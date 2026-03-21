#!/usr/bin/env python3
"""
Phase 4: Analyze Sobol results and generate publication-quality plots.

Collects batch results, runs Sobol analysis (if not already done),
generates figures and a markdown report.

Usage: python3 scripts/sensitivity/analyze_and_plot.py [--outdir results/sensitivity]
"""

import argparse
import json
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from scripts.sensitivity.param_spec import get_salib_problem, PARAM_SPEC
from scripts.sensitivity.spatial_runner import METRIC_NAMES

# Plotting
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


# ─── Dark theme ───────────────────────────────────────────────────────
DARK_BG = '#1a1a2e'
DARK_PANEL = '#16213e'
DARK_TEXT = '#e0e0e0'
ACCENT_1 = '#e94560'  # coral red
ACCENT_2 = '#0f3460'  # deep blue
ACCENT_3 = '#00d2ff'  # cyan
ACCENT_4 = '#ffd700'  # gold

plt.rcParams.update({
    'figure.facecolor': DARK_BG,
    'axes.facecolor': DARK_PANEL,
    'text.color': DARK_TEXT,
    'axes.labelcolor': DARK_TEXT,
    'xtick.color': DARK_TEXT,
    'ytick.color': DARK_TEXT,
    'axes.edgecolor': '#333366',
    'grid.color': '#333366',
    'grid.alpha': 0.3,
    'font.size': 10,
})


# ─── Friendly labels ─────────────────────────────────────────────────
PARAM_SHORT = {
    "disease.a_exposure": "Exposure rate",
    "disease.K_half": "Half-dose",
    "disease.sigma_1_eff": "Shedding I1",
    "disease.sigma_2_eff": "Shedding I2",
    "disease.sigma_D": "Shedding dead",
    "disease.rho_rec": "Recovery rate",
    "disease.mu_EI1_ref": "E→I1 rate",
    "disease.mu_I2D_ref": "I2→Death rate",
    "disease.P_env_max": "Background Vibrio",
    "disease.T_ref": "Vibrio T_opt",
    "population.F0": "Fecundity",
    "population.gamma_fert": "Fertilization γ",
    "population.settler_survival": "Settler survival",
    "population.alpha_srs": "SRS α",
    "population.senescence_age": "Senescence age",
    "population.k_growth": "Growth rate",
    "population.L_min_repro": "Min repro size",
    "genetics.n_resistance": "N resist",
    "genetics.n_tolerance": "N toler",
    "genetics.target_mean_t": "Mean t₀",
    "genetics.target_mean_c": "Mean c₀",
    "genetics.tau_max": "τ_max",
    "spawning.p_spontaneous_female": "Spawn rate ♀",
    "spawning.induction_female_to_male": "Cascade κ_fm",
    "disease.susceptibility_multiplier": "Immunosupp ×",
    "disease.T_vbnc": "VBNC T",
    "disease.s_min": "Salinity min",
}

METRIC_SHORT = {
    "pop_crash_pct": "Pop crash %",
    "final_pop_frac": "Final pop frac",
    "recovery": "Recovery (y/n)",
    "extinction": "Extinction",
    "peak_mortality": "Peak mortality",
    "time_to_nadir": "Time to nadir",
    "total_disease_deaths": "Total deaths",
    "resistance_shift_mean": "Resistance Δ (mean)",
    "resistance_shift_max": "Resistance Δ (max)",
    "va_retention_mean": "Va retention",
    "tolerance_shift_mean": "Toler Δ",
    "recovery_shift_mean": "Recov Δ",
    "n_extinct_nodes": "Extinct nodes",
    "north_south_mortality_gradient": "N→S mortality gradient",
    "fjord_protection_effect": "Fjord protection",
}


def short_name(name):
    return PARAM_SHORT.get(name, name.split(".")[-1])


def collect_batch_results(outdir):
    """Collect Y matrices from all batch result files."""
    # Load sample matrix
    sample_path = os.path.join(outdir, "sobol_samples.npz")
    if not os.path.exists(sample_path):
        raise FileNotFoundError(f"No sample matrix: {sample_path}")
    
    sdata = np.load(sample_path, allow_pickle=True)
    X = sdata["X"]
    param_names = list(sdata["param_names"])
    n_total = len(X)
    
    # Collect from batch files
    Y = np.full((n_total, len(METRIC_NAMES)), np.nan)
    
    # Try single-file first
    single = os.path.join(outdir, "sobol_batch0_results.npz")
    if os.path.exists(single):
        batch_files = sorted([
            f for f in os.listdir(outdir)
            if f.startswith("sobol_batch") and f.endswith("_results.npz")
        ])
        for bf in batch_files:
            data = np.load(os.path.join(outdir, bf))
            start = int(data["start_idx"])
            end = int(data["end_idx"])
            Y[start:end] = data["Y"]
            print(f"  Loaded {bf}: rows {start}-{end}")
    
    valid = np.sum(~np.isnan(Y[:, 0]))
    print(f"Total valid results: {valid}/{n_total} ({valid/n_total*100:.1f}%)")
    
    return X, Y, param_names


def plot_morris(outdir):
    """Plot Morris screening results."""
    path = os.path.join(outdir, "morris_screening.json")
    if not os.path.exists(path):
        print("No Morris results found — skipping")
        return
    
    with open(path) as f:
        data = json.load(f)
    
    param_names = data["param_names"]
    screened = set(data["screened_params"])
    
    # Pick a representative metric for the main plot
    # Use pop_crash_pct as primary
    n_metrics = len(data["morris_results"])
    fig, axes = plt.subplots(2, min(5, (n_metrics+1)//2), 
                              figsize=(20, 10), squeeze=False)
    fig.suptitle("Morris Screening: Elementary Effects (μ* vs σ)", 
                 fontsize=16, fontweight='bold', y=0.98)
    
    for idx, (metric, mdata) in enumerate(data["morris_results"].items()):
        ax = axes.flat[idx] if idx < len(axes.flat) else None
        if ax is None:
            break
        
        mu_star = np.array(mdata["mu_star"])
        sigma = np.array(mdata["sigma"])
        
        colors = [ACCENT_3 if p in screened else '#555555' for p in param_names]
        
        ax.scatter(mu_star, sigma, c=colors, s=50, alpha=0.8, edgecolors='white', linewidths=0.5)
        
        # Label top 5
        top5 = np.argsort(mu_star)[::-1][:5]
        for i in top5:
            ax.annotate(short_name(param_names[i]),
                       (mu_star[i], sigma[i]),
                       fontsize=7, color=DARK_TEXT,
                       xytext=(5, 5), textcoords='offset points')
        
        # Diagonal line (σ = μ* means non-linear / interactive)
        max_val = max(np.max(mu_star), np.max(sigma)) * 1.1
        ax.plot([0, max_val], [0, max_val], '--', color='#666', alpha=0.5, linewidth=1)
        
        ax.set_xlabel("μ* (mean absolute effect)")
        ax.set_ylabel("σ (interaction/non-linearity)")
        ax.set_title(METRIC_SHORT.get(metric, metric), fontsize=11)
    
    # Hide unused axes
    for idx in range(len(data["morris_results"]), len(axes.flat)):
        axes.flat[idx].set_visible(False)
    
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    out = os.path.join(outdir, "morris_screening.png")
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {out}")


def plot_sobol_bars(sobol_results, param_names, outdir):
    """Bar chart of S1 vs ST for each metric."""
    n_metrics = len(sobol_results)
    if n_metrics == 0:
        return
    
    n_cols = min(5, n_metrics)
    n_rows = (n_metrics + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4*n_cols, 4*n_rows), squeeze=False)
    fig.suptitle("Sobol Sensitivity Indices: S1 (first-order) vs ST (total)", 
                 fontsize=16, fontweight='bold', y=1.01)
    
    for idx, (metric, mdata) in enumerate(sobol_results.items()):
        ax = axes.flat[idx]
        
        S1 = np.array(mdata["S1"])
        ST = np.array(mdata["ST"])
        
        # Sort by ST
        order = np.argsort(ST)[::-1][:15]  # top 15
        
        y_pos = np.arange(len(order))
        labels = [short_name(param_names[i]) for i in order]
        
        ax.barh(y_pos + 0.15, ST[order], height=0.3, color=ACCENT_1, alpha=0.8, label='ST')
        ax.barh(y_pos - 0.15, S1[order], height=0.3, color=ACCENT_3, alpha=0.8, label='S1')
        
        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels, fontsize=8)
        ax.invert_yaxis()
        ax.set_xlabel("Index value")
        ax.set_title(METRIC_SHORT.get(metric, metric), fontsize=11)
        ax.legend(fontsize=8, loc='lower right')
        ax.set_xlim(-0.1, None)
        ax.axvline(x=0, color='#666', linewidth=0.5)
    
    for idx in range(len(sobol_results), len(axes.flat)):
        axes.flat[idx].set_visible(False)
    
    plt.tight_layout()
    out = os.path.join(outdir, "sobol_bars.png")
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {out}")


def plot_heatmap(sobol_results, param_names, outdir):
    """Heatmap: parameter × metric importance (ST values)."""
    metrics_with_data = [m for m in METRIC_NAMES if m in sobol_results]
    if len(metrics_with_data) < 2:
        return
    
    # Build matrix
    n_params = len(param_names)
    n_metrics = len(metrics_with_data)
    
    ST_matrix = np.zeros((n_params, n_metrics))
    for j, metric in enumerate(metrics_with_data):
        ST_matrix[:, j] = sobol_results[metric]["ST"]
    
    # Sort parameters by max ST across metrics
    max_ST = np.max(ST_matrix, axis=1)
    order = np.argsort(max_ST)[::-1]
    
    # Take top 20 (or all if fewer)
    n_show = min(20, n_params)
    order = order[:n_show]
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Custom colormap
    cmap = LinearSegmentedColormap.from_list("", [DARK_PANEL, "#0f3460", ACCENT_1, ACCENT_4])
    
    im = ax.imshow(ST_matrix[order], aspect='auto', cmap=cmap, vmin=0, vmax=max(0.5, np.max(ST_matrix)))
    
    ax.set_xticks(range(n_metrics))
    ax.set_xticklabels([METRIC_SHORT.get(m, m) for m in metrics_with_data],
                       rotation=45, ha='right', fontsize=9)
    ax.set_yticks(range(n_show))
    ax.set_yticklabels([short_name(param_names[i]) for i in order], fontsize=9)
    
    # Annotate cells
    for i in range(n_show):
        for j in range(n_metrics):
            val = ST_matrix[order[i], j]
            if val > 0.05:
                ax.text(j, i, f"{val:.2f}", ha='center', va='center',
                       fontsize=8, color='white' if val > 0.3 else DARK_TEXT)
    
    plt.colorbar(im, ax=ax, label="Total-order Sobol index (ST)")
    ax.set_title("Parameter × Metric Importance (ST)", fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    out = os.path.join(outdir, "sobol_heatmap.png")
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {out}")


def plot_interactions(sobol_results, param_names, outdir):
    """Plot interaction strength (ST - S1) for each metric."""
    metrics_with_data = [m for m in METRIC_NAMES if m in sobol_results]
    if len(metrics_with_data) < 2:
        return
    
    n_params = len(param_names)
    n_metrics = len(metrics_with_data)
    
    inter_matrix = np.zeros((n_params, n_metrics))
    for j, metric in enumerate(metrics_with_data):
        S1 = np.array(sobol_results[metric]["S1"])
        ST = np.array(sobol_results[metric]["ST"])
        inter_matrix[:, j] = np.maximum(0, ST - S1)  # clip negative to 0
    
    # Sort by max interaction
    max_inter = np.max(inter_matrix, axis=1)
    order = np.argsort(max_inter)[::-1][:15]
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    cmap = LinearSegmentedColormap.from_list("", [DARK_PANEL, "#1a4060", "#e94560"])
    
    im = ax.imshow(inter_matrix[order], aspect='auto', cmap=cmap, vmin=0)
    
    ax.set_xticks(range(n_metrics))
    ax.set_xticklabels([METRIC_SHORT.get(m, m) for m in metrics_with_data],
                       rotation=45, ha='right', fontsize=9)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([short_name(param_names[i]) for i in order], fontsize=9)
    
    for i in range(len(order)):
        for j in range(n_metrics):
            val = inter_matrix[order[i], j]
            if val > 0.02:
                ax.text(j, i, f"{val:.2f}", ha='center', va='center',
                       fontsize=8, color='white' if val > 0.15 else DARK_TEXT)
    
    plt.colorbar(im, ax=ax, label="Interaction strength (ST − S1)")
    ax.set_title("Parameter Interactions by Metric", fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    out = os.path.join(outdir, "sobol_interactions.png")
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {out}")


def plot_scatter_top_params(X, Y, param_names, sobol_results, outdir, top_n=5):
    """Scatter plots for top N parameters vs key metrics."""
    if not sobol_results:
        return
    
    # Find top parameters by max ST across all metrics
    n_params = len(param_names)
    max_ST = np.zeros(n_params)
    for metric, mdata in sobol_results.items():
        ST = np.array(mdata["ST"])
        max_ST = np.maximum(max_ST, ST)
    
    top_params = np.argsort(max_ST)[::-1][:top_n]
    
    key_metrics = ["pop_crash_pct", "resistance_shift", "recovery", "va_retention"]
    key_metrics = [m for m in key_metrics if m in sobol_results]
    
    n_cols = len(key_metrics)
    n_rows = top_n
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4*n_cols, 3*n_rows), squeeze=False)
    fig.suptitle(f"Top {top_n} Parameters vs Key Metrics", fontsize=14, fontweight='bold', y=1.01)
    
    for row, pi in enumerate(top_params):
        x_vals = X[:, pi]
        
        # Back-transform log-uniform params for display
        pname = param_names[pi]
        spec = PARAM_SPEC.get(pname, {})
        if spec.get("dist") == "loguniform":
            x_vals = 10.0 ** x_vals
            x_label = f"{short_name(pname)} (log)"
        elif spec.get("dist") == "discrete":
            vals = spec.get("values", [])
            x_vals_disc = np.array([vals[min(int(v), len(vals)-1)] for v in x_vals])
            x_vals = x_vals_disc
            x_label = short_name(pname)
        else:
            x_label = short_name(pname)
        
        for col, metric in enumerate(key_metrics):
            ax = axes[row, col]
            mi = METRIC_NAMES.index(metric)
            y_vals = Y[:, mi]
            
            valid = ~np.isnan(y_vals)
            ax.scatter(x_vals[valid], y_vals[valid], s=3, alpha=0.15, 
                      color=ACCENT_3, rasterized=True)
            
            # Add moving average trend line
            if np.sum(valid) > 20:
                sorted_idx = np.argsort(x_vals[valid])
                x_sorted = x_vals[valid][sorted_idx]
                y_sorted = y_vals[valid][sorted_idx]
                window = max(len(y_sorted) // 20, 5)
                y_smooth = np.convolve(y_sorted, np.ones(window)/window, mode='valid')
                x_smooth = x_sorted[window//2:window//2+len(y_smooth)]
                ax.plot(x_smooth, y_smooth, color=ACCENT_1, linewidth=2)
            
            if row == 0:
                ax.set_title(METRIC_SHORT.get(metric, metric), fontsize=10)
            if col == 0:
                ax.set_ylabel(x_label, fontsize=9)
            if row == n_rows - 1:
                ax.set_xlabel(METRIC_SHORT.get(metric, metric), fontsize=8)
    
    plt.tight_layout()
    out = os.path.join(outdir, "sobol_scatter.png")
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {out}")


def generate_report(outdir, param_names, sobol_results, morris_data=None):
    """Generate markdown report."""
    
    lines = [
        "# SSWD-EvoEpi Sensitivity Analysis Report",
        "",
        f"**Generated:** {__import__('datetime').datetime.now().strftime('%Y-%m-%d %H:%M PST')}",
        f"**Machine:** 16-core Ryzen 7 5700, 31GB RAM",
        "",
        "---",
        "",
    ]
    
    # Morris section
    if morris_data:
        lines.extend([
            "## Phase 1: Morris Screening",
            "",
            f"- **Parameters tested:** {morris_data['n_params']}",
            f"- **Trajectories:** {morris_data['n_trajectories']}",
            f"- **Total runs:** {morris_data['n_runs']}",
            f"- **Errors:** {morris_data['n_errors']}",
            f"- **Runtime:** {morris_data['elapsed_s']:.0f}s",
            "",
            f"### Screened parameters ({len(morris_data['screened_params'])}/{morris_data['n_params']}):",
            "",
        ])
        for p in morris_data['screened_params']:
            lines.append(f"- ✓ **{short_name(p)}** (`{p}`)")
        
        if morris_data['eliminated_params']:
            lines.extend(["", "### Eliminated (minimal effect):"])
            for p in morris_data['eliminated_params']:
                lines.append(f"- ✗ {short_name(p)} (`{p}`)")
        lines.append("")
    
    # Sobol section
    if sobol_results:
        lines.extend([
            "## Phase 2: Sobol Analysis",
            "",
            "### Top Parameters by Total-Order Index (ST)",
            "",
        ])
        
        for metric, mdata in sobol_results.items():
            S1 = np.array(mdata["S1"])
            ST = np.array(mdata["ST"])
            S1_conf = np.array(mdata["S1_conf"])
            ST_conf = np.array(mdata["ST_conf"])
            
            ranked = np.argsort(ST)[::-1][:10]
            
            lines.extend([
                f"#### {METRIC_SHORT.get(metric, metric)}",
                "",
                "| Rank | Parameter | S1 (±CI) | ST (±CI) | Interaction |",
                "|------|-----------|----------|----------|-------------|",
            ])
            
            for rank, idx in enumerate(ranked, 1):
                s1_str = f"{S1[idx]:.3f} ± {S1_conf[idx]:.3f}"
                st_str = f"{ST[idx]:.3f} ± {ST_conf[idx]:.3f}"
                inter = max(0, ST[idx] - S1[idx])
                lines.append(
                    f"| {rank} | {short_name(param_names[idx])} | {s1_str} | {st_str} | {inter:.3f} |"
                )
            lines.append("")
        
        # Summary: aggregate importance
        lines.extend([
            "### Aggregate Parameter Importance",
            "",
            "Ranked by maximum ST across all metrics:",
            "",
            "| Rank | Parameter | Max ST | Most Important For |",
            "|------|-----------|--------|-------------------|",
        ])
        
        max_ST = np.zeros(len(param_names))
        best_metric = [""] * len(param_names)
        
        for metric, mdata in sobol_results.items():
            ST = np.array(mdata["ST"])
            for i in range(len(param_names)):
                if ST[i] > max_ST[i]:
                    max_ST[i] = ST[i]
                    best_metric[i] = METRIC_SHORT.get(metric, metric)
        
        ranked = np.argsort(max_ST)[::-1]
        for rank, idx in enumerate(ranked, 1):
            lines.append(
                f"| {rank} | {short_name(param_names[idx])} | {max_ST[idx]:.3f} | {best_metric[idx]} |"
            )
        
        lines.extend([
            "",
            "### Key Findings",
            "",
            "*(Auto-generated — review for accuracy)*",
            "",
        ])
        
        # Auto-detect key findings
        top3_params = ranked[:3]
        lines.append(f"1. **Most influential parameters:** "
                    f"{', '.join(short_name(param_names[i]) for i in top3_params)}")
        
        # Check for high-interaction parameters
        max_inter = np.zeros(len(param_names))
        for metric, mdata in sobol_results.items():
            S1 = np.array(mdata["S1"])
            ST = np.array(mdata["ST"])
            inter = np.maximum(0, ST - S1)
            max_inter = np.maximum(max_inter, inter)
        
        high_inter = np.where(max_inter > 0.1)[0]
        if len(high_inter) > 0:
            lines.append(f"2. **High interaction effects:** "
                        f"{', '.join(short_name(param_names[i]) for i in high_inter)}")
        
        lines.extend([
            "",
            "### Figures",
            "",
            "- `morris_screening.png` — Morris μ* vs σ plots",
            "- `sobol_bars.png` — S1 vs ST bar charts per metric",
            "- `sobol_heatmap.png` — Parameter × metric importance matrix",
            "- `sobol_interactions.png` — Interaction strength heatmap",
            "- `sobol_scatter.png` — Top parameter scatter plots",
            "",
        ])
    
    report = "\n".join(lines)
    outpath = os.path.join(outdir, "SENSITIVITY_REPORT.md")
    with open(outpath, "w") as f:
        f.write(report)
    print(f"Report saved: {outpath}")
    
    return report


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", default="results/sensitivity")
    args = parser.parse_args()
    outdir = args.outdir
    
    # Load Morris data if available
    morris_data = None
    morris_path = os.path.join(outdir, "morris_screening.json")
    if os.path.exists(morris_path):
        with open(morris_path) as f:
            morris_data = json.load(f)
        print("Loaded Morris screening data")
        plot_morris(outdir)
    
    # Load Sobol results
    sobol_path = os.path.join(outdir, "sobol_results.json")
    sobol_results = None
    param_names = None
    X = None
    Y = None
    
    if os.path.exists(sobol_path):
        with open(sobol_path) as f:
            sdata = json.load(f)
        sobol_results = sdata["sobol_results"]
        param_names = sdata["param_names"]
        print(f"Loaded Sobol results: {len(sobol_results)} metrics analyzed")
    else:
        # Need to collect batch results and run analysis
        print("No sobol_results.json — collecting batch results...")
        try:
            X, Y, param_names = collect_batch_results(outdir)
            
            # Run Sobol analysis
            from SALib.analyze import sobol as sobol_analyze
            problem = get_salib_problem(param_names)
            
            sobol_results = {}
            for j, metric in enumerate(METRIC_NAMES):
                y = Y[:, j]
                valid = ~np.isnan(y)
                if np.sum(valid) < len(y) * 0.5:
                    continue
                y_clean = y.copy()
                if np.any(~valid):
                    y_clean[~valid] = np.nanmedian(y)
                if np.std(y_clean) < 1e-10:
                    continue
                
                try:
                    Si = sobol_analyze.analyze(problem, y_clean,
                        calc_second_order=False, num_resamples=1000,
                        conf_level=0.95, seed=54321)
                    sobol_results[metric] = {
                        "S1": Si["S1"].tolist(),
                        "S1_conf": Si["S1_conf"].tolist(),
                        "ST": Si["ST"].tolist(),
                        "ST_conf": Si["ST_conf"].tolist(),
                    }
                except Exception as e:
                    print(f"  {metric}: {e}")
            
            # Save
            with open(sobol_path, "w") as f:
                json.dump({"sobol_results": sobol_results, "param_names": param_names,
                          "metric_names": METRIC_NAMES}, f, indent=2)
        
        except Exception as e:
            print(f"Failed to collect/analyze batch results: {e}")
            import traceback
            traceback.print_exc()
    
    # Load X, Y if not already loaded (for scatter plots)
    if X is None:
        sample_path = os.path.join(outdir, "sobol_samples.npz")
        if os.path.exists(sample_path):
            sdata = np.load(sample_path, allow_pickle=True)
            X = sdata["X"]
        batch_path = os.path.join(outdir, "sobol_batch0_results.npz")
        if os.path.exists(batch_path):
            try:
                _, Y, _ = collect_batch_results(outdir)
            except:
                Y = None
    
    # Generate plots
    if sobol_results and param_names:
        plot_sobol_bars(sobol_results, param_names, outdir)
        plot_heatmap(sobol_results, param_names, outdir)
        plot_interactions(sobol_results, param_names, outdir)
        
        if X is not None and Y is not None:
            plot_scatter_top_params(X, Y, param_names, sobol_results, outdir)
    
    # Generate report
    generate_report(outdir, param_names or [], sobol_results or {}, morris_data)
    
    print("\nAnalysis complete!")


if __name__ == "__main__":
    main()
