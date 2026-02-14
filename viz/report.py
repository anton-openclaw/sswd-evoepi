"""Automated HTML summary report for SSWD-EvoEpi simulation runs.

Generates a self-contained HTML report with:
  a. Run metadata (config, seed, date, git hash)
  b. Executive summary (key numbers table)
  c. Population trajectories (all nodes)
  d. Disease dynamics (all nodes)
  e. Genetic evolution (allele freq plots, r̄ trajectory)
  f. Network dynamics (connectivity, dispersal stats)
  g. Anomalies/warnings (extinction, NaN, etc.)

Usage:
    python -m viz.report results/5node_epidemic_20yr/

Output:
    results/5node_epidemic_20yr/report.html
"""

from __future__ import annotations

import argparse
import base64
import datetime
import io
import json
import subprocess
import sys
from pathlib import Path
from typing import List, Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from viz.plot_utils import (
    SimulationData,
    NODE_COLORS,
    plot_population_all_nodes,
    plot_disease_all_nodes,
    plot_genetics_summary,
    plot_genetic_heatmap_all_nodes,
    plot_network_state,
    plot_allele_frequencies,
    plot_resistance_distribution,
)


# ═══════════════════════════════════════════════════════════════════════
# UTILITIES
# ═══════════════════════════════════════════════════════════════════════

def _fig_to_base64(fig: plt.Figure, dpi: int = 130) -> str:
    """Convert a matplotlib figure to base64-encoded PNG data URI."""
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight",
                facecolor="white", edgecolor="none")
    plt.close(fig)
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode("ascii")
    return f"data:image/png;base64,{b64}"


def _get_git_hash() -> str:
    """Get current git commit hash, or 'unknown'."""
    try:
        result = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, timeout=5,
        )
        return result.stdout.strip() if result.returncode == 0 else "unknown"
    except Exception:
        return "unknown"


def _detect_anomalies(data: SimulationData) -> List[dict]:
    """Scan simulation data for anomalies and warnings.

    Returns list of {severity, node, year, message} dicts.
    """
    anomalies = []

    for i in range(data.n_nodes):
        name = data.node_names[i]

        # Node extinction
        final_pop = data.yearly_pop[i, -1]
        if final_pop == 0:
            anomalies.append({
                "severity": "🔴 CRITICAL",
                "node": name,
                "year": int(np.argmin(data.yearly_pop[i] > 0)),
                "message": f"Population went EXTINCT (N=0)",
            })
        elif final_pop <= 5:
            yr_below = int(np.argmax(data.yearly_pop[i] <= 5))
            anomalies.append({
                "severity": "🟡 WARNING",
                "node": name,
                "year": yr_below,
                "message": f"Near-extinction: N={final_pop} at end of simulation",
            })

        # Population crash > 99%
        if data.disease_year < data.n_years:
            pre_pop = data.yearly_pop[i, data.disease_year]
            min_pop = data.yearly_pop[i, data.disease_year:].min()
            if pre_pop > 0:
                crash = 1 - min_pop / pre_pop
                if crash > 0.99:
                    anomalies.append({
                        "severity": "🟡 WARNING",
                        "node": name,
                        "year": int(data.disease_year + np.argmin(
                            data.yearly_pop[i, data.disease_year:])),
                        "message": f"Population crash > 99% ({crash:.1%})",
                    })

        # Check for NaN in key metrics
        for metric_name, metric_arr in [
            ("mean_resistance", data.yearly_mean_resistance[i]),
            ("ef1a_freq", data.yearly_ef1a_freq[i]),
            ("va", data.yearly_va[i]),
        ]:
            nan_count = np.isnan(metric_arr).sum()
            if nan_count > 0:
                first_nan = int(np.argmax(np.isnan(metric_arr)))
                anomalies.append({
                    "severity": "🔴 CRITICAL",
                    "node": name,
                    "year": first_nan,
                    "message": f"NaN detected in {metric_name} ({nan_count} years)",
                })

        # Negative resistance (should never happen)
        if np.any(data.yearly_mean_resistance[i] < 0):
            anomalies.append({
                "severity": "🔴 CRITICAL",
                "node": name,
                "year": int(np.argmax(data.yearly_mean_resistance[i] < 0)),
                "message": f"Negative mean resistance detected",
            })

        # EF1A fixation or loss
        ef1a = data.yearly_ef1a_freq[i]
        if np.any(ef1a >= 0.99):
            anomalies.append({
                "severity": "🟡 WARNING",
                "node": name,
                "year": int(np.argmax(ef1a >= 0.99)),
                "message": f"EF1A allele near fixation (q ≥ 0.99)",
            })
        if np.any(ef1a <= 0.01) and ef1a[0] > 0.01:
            anomalies.append({
                "severity": "🟡 WARNING",
                "node": name,
                "year": int(np.argmax(ef1a <= 0.01)),
                "message": f"EF1A allele nearly lost (q ≤ 0.01)",
            })

        # Zero recruits for 3+ consecutive years
        recruits = data.yearly_recruits[i]
        zero_run = 0
        max_zero_run = 0
        for yr in range(data.n_years):
            if recruits[yr] == 0:
                zero_run += 1
                max_zero_run = max(max_zero_run, zero_run)
            else:
                zero_run = 0
        if max_zero_run >= 3:
            anomalies.append({
                "severity": "🟡 WARNING",
                "node": name,
                "year": -1,
                "message": f"Zero recruits for {max_zero_run} consecutive years",
            })

    # Check total population
    total_pop = data.yearly_total_pop
    if total_pop[-1] == 0:
        anomalies.append({
            "severity": "🔴 CRITICAL",
            "node": "ALL",
            "year": int(np.argmin(total_pop > 0)),
            "message": "TOTAL population extinction across all nodes",
        })

    # No recovery check (monotonic decline post-epidemic)
    if data.disease_year + 3 < data.n_years:
        post_pop = total_pop[data.disease_year:]
        diffs = np.diff(post_pop)
        if np.all(diffs <= 0) and len(diffs) >= 5:
            anomalies.append({
                "severity": "🔵 INFO",
                "node": "ALL",
                "year": -1,
                "message": "Monotonic population decline post-epidemic (no recovery phase observed)",
            })

    return anomalies


# ═══════════════════════════════════════════════════════════════════════
# REPORT SECTIONS
# ═══════════════════════════════════════════════════════════════════════

def _section_metadata(data: SimulationData, result_dir: Path) -> str:
    """Section A: Run metadata."""
    git_hash = _get_git_hash()
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Try to load config if present
    config_path = result_dir / "config.yaml"
    config_summary = ""
    if config_path.exists():
        config_summary = f'<pre style="max-height:300px;overflow:auto;font-size:11px;background:#f8f8f8;padding:10px;border-radius:4px;">{config_path.read_text()[:3000]}</pre>'
    else:
        config_summary = '<p style="color:#999;">No config.yaml saved in results directory.</p>'

    return f"""
    <div class="section">
        <h2>A. Run Metadata</h2>
        <table class="meta-table">
            <tr><td>Simulation</td><td>{data.n_years} years, {data.n_nodes} nodes</td></tr>
            <tr><td>Disease introduced</td><td>Year {data.disease_year}</td></tr>
            <tr><td>Random seed</td><td>{data.seed}</td></tr>
            <tr><td>Results directory</td><td><code>{result_dir}</code></td></tr>
            <tr><td>Git commit</td><td><code>{git_hash}</code></td></tr>
            <tr><td>Report generated</td><td>{now}</td></tr>
            <tr><td>Node names</td><td>{", ".join(data.node_names)}</td></tr>
            <tr><td>Carrying capacities</td><td>{", ".join(str(int(k)) for k in data.node_K)}</td></tr>
        </table>
        <details>
            <summary>Config file</summary>
            {config_summary}
        </details>
    </div>
    """


def _section_executive_summary(data: SimulationData) -> str:
    """Section B: Executive summary table."""
    rows = []
    for i in range(data.n_nodes):
        K = int(data.node_K[i])
        pre_pop = int(data.yearly_pop[i, data.disease_year])
        final_pop = int(data.yearly_pop[i, -1])
        min_pop = int(data.yearly_pop[i, data.disease_year:].min())
        total_dis = int(data.yearly_disease_deaths[i].sum())
        peak_prev = float(data.peak_disease_prevalence[i])
        r_pre = float(data.yearly_mean_resistance[i, data.disease_year])
        r_post = float(data.yearly_mean_resistance[i, -1])
        delta_r = r_post - r_pre
        crash = (1 - min_pop / max(pre_pop, 1)) * 100

        crash_class = "critical" if crash >= 95 else "warning" if crash >= 80 else "ok"
        delta_class = "ok" if delta_r > 0 else "warning"

        rows.append(f"""
        <tr>
            <td style="color:{NODE_COLORS[i]};font-weight:bold">{data.node_names[i]}</td>
            <td>{K}</td><td>{pre_pop}</td><td>{final_pop}</td><td>{min_pop}</td>
            <td class="{crash_class}">{crash:.1f}%</td>
            <td>{total_dis}</td>
            <td>{peak_prev:.1%}</td>
            <td>{r_pre:.4f}</td><td>{r_post:.4f}</td>
            <td class="{delta_class}">{delta_r:+.4f}</td>
        </tr>""")

    # Totals
    total_pre = int(data.yearly_total_pop[data.disease_year])
    total_final = int(data.yearly_total_pop[-1])
    total_dis = int(data.yearly_disease_deaths.sum())
    total_crash = (1 - total_final / max(total_pre, 1)) * 100

    rows.append(f"""
    <tr class="total-row">
        <td>TOTAL</td>
        <td>{int(data.node_K.sum())}</td><td>{total_pre}</td>
        <td>{total_final}</td><td>—</td>
        <td>{total_crash:.1f}%</td>
        <td>{total_dis}</td>
        <td>—</td><td colspan="3">—</td>
    </tr>""")

    return f"""
    <div class="section">
        <h2>B. Executive Summary</h2>
        <table class="data-table">
            <thead>
                <tr>
                    <th>Node</th><th>K</th><th>Pre-N</th><th>Final N</th>
                    <th>Min N</th><th>Crash %</th><th>Dis. Deaths</th>
                    <th>Peak Prev.</th><th>r̄ pre</th><th>r̄ final</th><th>Δr̄</th>
                </tr>
            </thead>
            <tbody>{"".join(rows)}</tbody>
        </table>
    </div>
    """


def _section_population(data: SimulationData) -> str:
    """Section C: Population trajectory plots."""
    fig = plot_population_all_nodes(data)
    img = _fig_to_base64(fig, dpi=130)

    return f"""
    <div class="section">
        <h2>C. Population Trajectories</h2>
        <p>Stacked areas show juveniles/subadults (orange) and adults (blue).
        Annual recruits shown as green bars. Dashed red line = disease introduction.
        Dotted black line = carrying capacity K.</p>
        <img src="{img}" style="width:100%;max-width:1200px;" alt="Population trajectories">
    </div>
    """


def _section_disease(data: SimulationData) -> str:
    """Section D: Disease dynamics plots."""
    fig = plot_disease_all_nodes(data)
    img = _fig_to_base64(fig, dpi=130)

    return f"""
    <div class="section">
        <h2>D. Disease Dynamics</h2>
        <p>Stacked bars: gray = natural deaths, red = disease deaths.
        Disease mortality dominates after introduction year.</p>
        <img src="{img}" style="width:100%;max-width:1200px;" alt="Disease dynamics">
    </div>
    """


def _section_genetics(data: SimulationData) -> str:
    """Section E: Genetic evolution."""
    # Multi-panel summary
    fig_summary = plot_genetics_summary(data)
    img_summary = _fig_to_base64(fig_summary, dpi=130)

    # Heatmap
    fig_heatmap = plot_genetic_heatmap_all_nodes(data)
    img_heatmap = _fig_to_base64(fig_heatmap, dpi=130)

    # Per-node allele frequency plots
    allele_imgs = []
    for i in range(data.n_nodes):
        fig = plot_allele_frequencies(data, i)
        allele_imgs.append(_fig_to_base64(fig, dpi=110))

    allele_html = ""
    for i, img in enumerate(allele_imgs):
        allele_html += f"""
        <div style="display:inline-block;margin:5px;">
            <img src="{img}" style="width:100%;max-width:580px;" alt="Allele freq {data.node_names[i]}">
        </div>"""

    return f"""
    <div class="section">
        <h2>E. Genetic Evolution</h2>
        <h3>E.1 Summary (r̄, EF1A, Ne/N, V_A)</h3>
        <img src="{img_summary}" style="width:100%;max-width:1200px;" alt="Genetics summary">

        <h3>E.2 Allele Frequency Shifts (Δq = post − pre epidemic)</h3>
        <p>Heatmap: blue = decrease, red = increase. Gold line marks EF1A (locus 51).</p>
        <img src="{img_heatmap}" style="width:100%;max-width:1200px;" alt="Genetic heatmap">

        <h3>E.3 Top-3 Locus Trajectories (Per Node)</h3>
        {allele_html}
    </div>
    """


def _section_network(data: SimulationData) -> str:
    """Section F: Network dynamics."""
    # Network state at 3 key years
    years = [0, data.disease_year, min(data.disease_year + 5, data.n_years - 1), data.n_years - 1]
    years = sorted(set(years))  # deduplicate

    imgs = []
    for yr in years:
        fig = plot_network_state(data, yr)
        imgs.append((yr, _fig_to_base64(fig, dpi=110)))

    network_html = ""
    for yr, img in imgs:
        network_html += f"""
        <div style="display:inline-block;margin:5px;vertical-align:top;">
            <img src="{img}" style="width:100%;max-width:400px;" alt="Network year {yr}">
        </div>"""

    # Dispersal stats
    total_larvae = data.yearly_total_larvae
    dispersal_html = "<table class='meta-table'>"
    dispersal_html += "<tr><th>Year</th><th>Total Larvae Dispersed</th></tr>"
    for yr in range(data.n_years):
        dispersal_html += f"<tr><td>{yr}</td><td>{int(total_larvae[yr]):,}</td></tr>"
    dispersal_html += "</table>"

    return f"""
    <div class="section">
        <h2>F. Network Dynamics</h2>
        <h3>F.1 Network State Snapshots</h3>
        <p>Node size ∝ population, color = disease mortality fraction (green→red).
        Gray lines show connectivity.</p>
        {network_html}

        <h3>F.2 Larval Dispersal</h3>
        <details>
            <summary>Yearly total larvae dispersed (click to expand)</summary>
            {dispersal_html}
        </details>
    </div>
    """


def _section_anomalies(data: SimulationData) -> str:
    """Section G: Anomalies and warnings."""
    anomalies = _detect_anomalies(data)

    if not anomalies:
        return """
        <div class="section">
            <h2>G. Anomalies & Warnings</h2>
            <p style="color:#2ca02c;font-weight:bold;">✅ No anomalies detected.</p>
        </div>
        """

    rows = []
    for a in anomalies:
        yr_str = str(a["year"]) if a["year"] >= 0 else "—"
        rows.append(f"""
        <tr>
            <td>{a['severity']}</td>
            <td>{a['node']}</td>
            <td>{yr_str}</td>
            <td>{a['message']}</td>
        </tr>""")

    critical_count = sum(1 for a in anomalies if "CRITICAL" in a["severity"])
    warning_count = sum(1 for a in anomalies if "WARNING" in a["severity"])
    info_count = sum(1 for a in anomalies if "INFO" in a["severity"])

    return f"""
    <div class="section">
        <h2>G. Anomalies & Warnings</h2>
        <p>Found <strong>{len(anomalies)}</strong> items:
        <span style="color:red">{critical_count} critical</span>,
        <span style="color:orange">{warning_count} warnings</span>,
        <span style="color:blue">{info_count} info</span></p>
        <table class="data-table">
            <thead>
                <tr><th>Severity</th><th>Node</th><th>Year</th><th>Message</th></tr>
            </thead>
            <tbody>{"".join(rows)}</tbody>
        </table>
    </div>
    """


# ═══════════════════════════════════════════════════════════════════════
# REPORT ASSEMBLY
# ═══════════════════════════════════════════════════════════════════════

CSS = """
* { box-sizing: border-box; }
body {
    font-family: 'Segoe UI', 'Helvetica Neue', Arial, sans-serif;
    background: #f5f7fa;
    margin: 0;
    padding: 0;
    color: #333;
    line-height: 1.5;
}
.header {
    background: linear-gradient(135deg, #1a1a2e 0%, #16213e 50%, #0f3460 100%);
    color: white;
    padding: 30px 40px;
    text-align: center;
}
.header h1 { margin: 0 0 5px; font-size: 26px; }
.header .subtitle { font-size: 13px; opacity: 0.7; }
.nav {
    background: white;
    padding: 10px 20px;
    position: sticky;
    top: 0;
    z-index: 100;
    box-shadow: 0 2px 4px rgba(0,0,0,0.1);
    display: flex;
    gap: 8px;
    flex-wrap: wrap;
    justify-content: center;
}
.nav a {
    text-decoration: none;
    color: #0f3460;
    padding: 5px 14px;
    border-radius: 16px;
    font-size: 12px;
    font-weight: 600;
    border: 1px solid #ddd;
    transition: all 0.2s;
}
.nav a:hover { background: #0f3460; color: white; }
.section {
    background: white;
    margin: 16px 20px;
    padding: 20px 28px;
    border-radius: 6px;
    box-shadow: 0 1px 5px rgba(0,0,0,0.08);
}
.section h2 {
    color: #0f3460;
    border-bottom: 2px solid #eee;
    padding-bottom: 8px;
    margin-top: 0;
}
.section h3 { color: #555; margin-top: 20px; }
.section img { border-radius: 4px; display: block; margin: 10px auto; }
.meta-table {
    border-collapse: collapse;
    font-size: 13px;
    margin: 10px 0;
}
.meta-table td, .meta-table th {
    padding: 6px 14px;
    border-bottom: 1px solid #eee;
    text-align: left;
}
.meta-table td:first-child { font-weight: 600; color: #555; white-space: nowrap; }
.data-table {
    border-collapse: collapse;
    width: 100%;
    font-size: 12px;
    margin: 10px 0;
}
.data-table th {
    background: #333;
    color: white;
    padding: 8px;
    text-align: center;
}
.data-table td {
    padding: 6px 8px;
    text-align: center;
    border-bottom: 1px solid #eee;
}
.data-table .critical { color: #e41a1c; font-weight: bold; }
.data-table .warning { color: #ff7f0e; font-weight: bold; }
.data-table .ok { color: #2ca02c; }
.data-table .total-row {
    border-top: 2px solid #333;
    font-weight: bold;
    background: #f0f0f0;
}
details { margin: 10px 0; }
summary {
    cursor: pointer;
    color: #0f3460;
    font-weight: 600;
    font-size: 13px;
}
code {
    background: #f0f0f0;
    padding: 2px 6px;
    border-radius: 3px;
    font-size: 12px;
}
.footer {
    text-align: center;
    padding: 20px;
    color: #999;
    font-size: 11px;
}
"""


def generate_report(
    result_dir: str | Path,
    output_path: Optional[str | Path] = None,
) -> Path:
    """Generate a self-contained HTML report from simulation results.

    Args:
        result_dir: Directory with simulation_data.npz + metadata.json.
        output_path: Output HTML path (default: <result_dir>/report.html).

    Returns:
        Path to generated report.
    """
    result_dir = Path(result_dir)
    if output_path is None:
        output_path = result_dir / "report.html"
    output_path = Path(output_path)

    print(f"Loading simulation data from {result_dir}...")
    data = SimulationData.load(result_dir)

    print("Generating report sections...")

    print("  A. Metadata...")
    sec_meta = _section_metadata(data, result_dir)

    print("  B. Executive summary...")
    sec_exec = _section_executive_summary(data)

    print("  C. Population trajectories...")
    sec_pop = _section_population(data)

    print("  D. Disease dynamics...")
    sec_disease = _section_disease(data)

    print("  E. Genetic evolution...")
    sec_genetics = _section_genetics(data)

    print("  F. Network dynamics...")
    sec_network = _section_network(data)

    print("  G. Anomalies & warnings...")
    sec_anomalies = _section_anomalies(data)

    # Assemble
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    git_hash = _get_git_hash()

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>SSWD-EvoEpi Report — {data.n_years}yr {data.n_nodes}-node</title>
    <style>{CSS}</style>
</head>
<body>
    <div class="header">
        <h1>🌊 SSWD-EvoEpi Simulation Report</h1>
        <div class="subtitle">
            {data.n_years}-year, {data.n_nodes}-node epidemic simulation |
            Seed {data.seed} | {now}
        </div>
    </div>

    <div class="nav">
        <a href="#sec-a">📋 Metadata</a>
        <a href="#sec-b">📊 Summary</a>
        <a href="#sec-c">📈 Population</a>
        <a href="#sec-d">🦠 Disease</a>
        <a href="#sec-e">🧬 Genetics</a>
        <a href="#sec-f">🗺️ Network</a>
        <a href="#sec-g">⚠️ Anomalies</a>
    </div>

    <div id="sec-a">{sec_meta}</div>
    <div id="sec-b">{sec_exec}</div>
    <div id="sec-c">{sec_pop}</div>
    <div id="sec-d">{sec_disease}</div>
    <div id="sec-e">{sec_genetics}</div>
    <div id="sec-f">{sec_network}</div>
    <div id="sec-g">{sec_anomalies}</div>

    <div class="footer">
        SSWD-EvoEpi Automated Report | git {git_hash} |
        Generated {now} | Weertman Lab, UW FHL
    </div>
</body>
</html>"""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(html)

    size_kb = output_path.stat().st_size / 1024
    print(f"\n✅ Report written to {output_path}")
    print(f"   File size: {size_kb:.0f} KB")
    return output_path


# ═══════════════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Generate automated HTML report from SSWD-EvoEpi simulation results.",
        epilog="Example: python -m viz.report results/5node_epidemic_20yr/",
    )
    parser.add_argument("result_dir", help="Directory with simulation_data.npz + metadata.json")
    parser.add_argument("--output", "-o", default=None,
                        help="Output HTML path (default: <result_dir>/report.html)")

    args = parser.parse_args()
    result_dir = Path(args.result_dir)
    if not result_dir.exists():
        print(f"Error: {result_dir} does not exist", file=sys.stderr)
        sys.exit(1)

    generate_report(result_dir, args.output)


if __name__ == "__main__":
    main()
