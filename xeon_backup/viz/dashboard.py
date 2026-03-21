"""Interactive HTML dashboard for SSWD-EvoEpi simulation results.

Generates a self-contained HTML file with 6 panels using Plotly.js:
  1. Network Overview (map with time slider)
  2. Population Dynamics (per node, stacked areas)
  3. Disease Dynamics (deaths + Vibrio)
  4. Genetics (resistance, EF1A, Ne/N, Va, allele freq top loci)
  5. Genetic Heatmap (Δq across 52 loci, all nodes)
  6. Summary Statistics Table

Usage:
    python -m viz.dashboard results/5node_epidemic_20yr/ --output dashboard.html

"If I cannot see it, I cannot trust it." — Willem
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import List

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px

# ═══════════════════════════════════════════════════════════════════════
# DATA LOADING (standalone — no dependency on plot_utils at import time)
# ═══════════════════════════════════════════════════════════════════════

def load_data(result_dir: str | Path) -> dict:
    """Load simulation data into a plain dict."""
    result_dir = Path(result_dir)

    with open(result_dir / "metadata.json") as f:
        meta = json.load(f)

    d = np.load(result_dir / "simulation_data.npz", allow_pickle=True)

    return {
        "meta": meta,
        "yearly_pop": d["yearly_pop"],
        "yearly_adults": d["yearly_adults"],
        "yearly_recruits": d["yearly_recruits"],
        "yearly_natural_deaths": d["yearly_natural_deaths"],
        "yearly_disease_deaths": d["yearly_disease_deaths"],
        "yearly_mean_resistance": d["yearly_mean_resistance"],
        "yearly_vibrio_max": d["yearly_vibrio_max"],
        "yearly_ef1a_freq": d["yearly_ef1a_freq"],
        "yearly_va": d["yearly_va"],
        "yearly_ne_ratio": d["yearly_ne_ratio"],
        "yearly_allele_freq_top3": d["yearly_allele_freq_top3"],
        "yearly_total_pop": d["yearly_total_pop"],
        "yearly_total_larvae": d["yearly_total_larvae_dispersed"],
        "peak_disease_prevalence": d["peak_disease_prevalence"],
        "pre_epidemic_allele_freq": d["pre_epidemic_allele_freq"],
        "post_epidemic_allele_freq": d["post_epidemic_allele_freq"],
        "node_K": d["node_K"],
    }


# ═══════════════════════════════════════════════════════════════════════
# NODE POSITIONS (hardcoded for the 5-node test network)
# ═══════════════════════════════════════════════════════════════════════

NODE_POSITIONS = {
    "Sitka, AK": {"lat": 57.06, "lon": -135.34},
    "Howe Sound, BC": {"lat": 49.52, "lon": -123.25},
    "San Juan Islands, WA": {"lat": 48.53, "lon": -123.02},
    "Newport, OR": {"lat": 44.63, "lon": -124.05},
    "Monterey, CA": {"lat": 36.62, "lon": -121.90},
}

NODE_COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"]
NODE_SYMBOLS = ["circle", "diamond", "square", "triangle-up", "star"]


# ═══════════════════════════════════════════════════════════════════════
# PANEL BUILDERS
# ═══════════════════════════════════════════════════════════════════════

def build_network_panel(data: dict) -> go.Figure:
    """Panel 1: Interactive network map with time slider.

    Nodes on map (lat/lon), size = population, color = disease prevalence.
    Connectivity edges. Time slider scrubs through years.
    """
    meta = data["meta"]
    n_nodes = meta["n_nodes"]
    n_years = meta["n_years"]
    names = meta["node_names"]
    disease_year = meta["disease_year"]

    fig = go.Figure()

    # Pre-compute all frames
    frames = []
    for yr in range(n_years):
        pops = data["yearly_pop"][:, yr]
        dis_d = data["yearly_disease_deaths"][:, yr]
        nat_d = data["yearly_natural_deaths"][:, yr]
        total_d = dis_d + nat_d

        # Disease fraction for color
        dis_frac = np.where(total_d > 0, dis_d / np.maximum(total_d, 1), 0.0)

        lats = [NODE_POSITIONS.get(n, {}).get("lat", 0) for n in names]
        lons = [NODE_POSITIONS.get(n, {}).get("lon", 0) for n in names]

        # Size: proportional to population, min 10
        sizes = [max(10, p / max(data["node_K"].max(), 1) * 60) for p in pops]

        # Color: green to red based on disease fraction
        colors = [f"rgb({int(255*df)}, {int(255*(1-df))}, 50)" for df in dis_frac]

        hover_text = [
            f"<b>{names[i]}</b><br>"
            f"Population: {pops[i]}<br>"
            f"K: {data['node_K'][i]}<br>"
            f"Disease Deaths: {dis_d[i]}<br>"
            f"Natural Deaths: {nat_d[i]}<br>"
            f"Disease Fraction: {dis_frac[i]:.1%}<br>"
            f"Mean r̄: {data['yearly_mean_resistance'][i, yr]:.4f}"
            for i in range(n_nodes)
        ]

        node_trace = go.Scattergeo(
            lat=lats, lon=lons,
            mode="markers+text",
            marker=dict(size=sizes, color=colors, line=dict(width=1.5, color="black"),
                        opacity=0.9),
            text=[f"{n}\nN={p}" for n, p in zip(names, pops)],
            textposition="top center",
            textfont=dict(size=9),
            hovertext=hover_text,
            hoverinfo="text",
            name="Nodes",
        )

        # Edge traces
        edge_lats = []
        edge_lons = []
        for i in range(n_nodes):
            for j in range(i + 1, n_nodes):
                edge_lats += [lats[i], lats[j], None]
                edge_lons += [lons[i], lons[j], None]

        edge_trace = go.Scattergeo(
            lat=edge_lats, lon=edge_lons,
            mode="lines",
            line=dict(width=0.8, color="rgba(150,150,150,0.4)"),
            hoverinfo="skip",
            name="Connectivity",
        )

        frame_label = f"Year {yr}" + (" ← disease" if yr == disease_year else "")
        frames.append(go.Frame(
            data=[edge_trace, node_trace],
            name=str(yr),
            layout=go.Layout(title_text=f"Network State — Year {yr}")
        ))

    # Initial state (year 0)
    if frames:
        for trace in frames[0].data:
            fig.add_trace(trace)

    fig.frames = frames

    # Slider
    sliders = [dict(
        active=0,
        yanchor="top", xanchor="left",
        currentvalue=dict(prefix="Year: ", font=dict(size=14)),
        transition=dict(duration=100),
        pad=dict(b=10, t=50),
        len=0.9, x=0.05, y=0,
        steps=[dict(
            args=[[str(yr)], dict(frame=dict(duration=200, redraw=True),
                                  mode="immediate")],
            label=str(yr),
            method="animate",
        ) for yr in range(n_years)]
    )]

    # Play/Pause buttons
    updatemenus = [dict(
        type="buttons",
        showactive=False,
        y=0, x=0.05, xanchor="right", yanchor="top",
        pad=dict(t=87, r=10),
        buttons=[
            dict(label="▶ Play", method="animate",
                 args=[None, dict(frame=dict(duration=500, redraw=True),
                                  fromcurrent=True,
                                  transition=dict(duration=200))]),
            dict(label="⏸ Pause", method="animate",
                 args=[[None], dict(frame=dict(duration=0, redraw=True),
                                    mode="immediate",
                                    transition=dict(duration=0))]),
        ]
    )]

    fig.update_layout(
        title=dict(text="Panel 1: Network Overview", font=dict(size=18)),
        geo=dict(
            scope="north america",
            showland=True, landcolor="rgb(243, 243, 243)",
            showocean=True, oceancolor="rgb(204, 229, 255)",
            showcoastlines=True, coastlinecolor="rgb(100, 100, 100)",
            showlakes=True, lakecolor="rgb(204, 229, 255)",
            lonaxis=dict(range=[-140, -118]),
            lataxis=dict(range=[34, 60]),
            projection_type="mercator",
        ),
        height=650,
        sliders=sliders,
        updatemenus=updatemenus,
        showlegend=False,
        margin=dict(l=20, r=20, t=60, b=20),
    )

    return fig


def build_population_panel(data: dict) -> go.Figure:
    """Panel 2: Population dynamics per node (stacked areas + K)."""
    meta = data["meta"]
    n_nodes = meta["n_nodes"]
    n_years = meta["n_years"]
    names = meta["node_names"]
    disease_year = meta["disease_year"]
    years = list(range(n_years))

    fig = make_subplots(
        rows=n_nodes, cols=1,
        shared_xaxes=True,
        subplot_titles=[f"{names[i]} (K={data['node_K'][i]})" for i in range(n_nodes)],
        vertical_spacing=0.04,
    )

    for i in range(n_nodes):
        pop = data["yearly_pop"][i]
        adults = data["yearly_adults"][i]
        non_adult = pop - adults
        recruits = data["yearly_recruits"][i]
        K = int(data["node_K"][i])
        row = i + 1
        show_legend = (i == 0)

        # Non-adults (juveniles + subadults)
        fig.add_trace(go.Scatter(
            x=years, y=non_adult.tolist(),
            mode="lines", fill="tozeroy",
            fillcolor="rgba(252, 141, 98, 0.4)",
            line=dict(width=0, color="rgba(252, 141, 98, 0.4)"),
            name="Juv + Subadult",
            legendgroup="nonadult", showlegend=show_legend,
            hovertemplate="Juv+Sub: %{y}<extra></extra>",
        ), row=row, col=1)

        # Adults (stacked on top)
        fig.add_trace(go.Scatter(
            x=years, y=pop.tolist(),
            mode="lines", fill="tonexty",
            fillcolor="rgba(141, 160, 203, 0.6)",
            line=dict(width=0.5, color="rgba(141, 160, 203, 0.8)"),
            name="Adults",
            legendgroup="adult", showlegend=show_legend,
            hovertemplate="Total: %{y}<extra></extra>",
        ), row=row, col=1)

        # Recruits as bar-like scatter
        fig.add_trace(go.Bar(
            x=years, y=recruits.tolist(),
            marker_color="rgba(102, 194, 165, 0.7)",
            name="Recruits",
            legendgroup="recruits", showlegend=show_legend,
            hovertemplate="Recruits: %{y}<extra></extra>",
            width=0.5,
        ), row=row, col=1)

        # K line
        fig.add_trace(go.Scatter(
            x=[0, n_years - 1], y=[K, K],
            mode="lines",
            line=dict(dash="dot", color="black", width=1),
            name="K",
            legendgroup="K", showlegend=show_legend,
            hoverinfo="skip",
        ), row=row, col=1)

        # Disease year marker
        fig.add_vline(x=disease_year, line_dash="dash", line_color="red",
                      line_width=1.5, opacity=0.7, row=row, col=1)

    fig.update_layout(
        title=dict(text="Panel 2: Population Dynamics by Node", font=dict(size=18)),
        height=250 * n_nodes + 80,
        legend=dict(orientation="h", yanchor="bottom", y=1.01, xanchor="center", x=0.5),
        margin=dict(l=60, r=20, t=100, b=40),
        barmode="overlay",
    )
    fig.update_xaxes(title_text="Year", row=n_nodes, col=1)

    return fig


def build_disease_panel(data: dict) -> go.Figure:
    """Panel 3: Disease dynamics per node (mortality + Vibrio)."""
    meta = data["meta"]
    n_nodes = meta["n_nodes"]
    n_years = meta["n_years"]
    names = meta["node_names"]
    disease_year = meta["disease_year"]
    years = list(range(n_years))

    fig = make_subplots(
        rows=n_nodes, cols=1,
        shared_xaxes=True,
        subplot_titles=names,
        vertical_spacing=0.04,
        specs=[[{"secondary_y": True}] for _ in range(n_nodes)],
    )

    for i in range(n_nodes):
        nat_d = data["yearly_natural_deaths"][i]
        dis_d = data["yearly_disease_deaths"][i]
        vibrio = data["yearly_vibrio_max"][i]
        row = i + 1
        show_legend = (i == 0)

        # Natural deaths
        fig.add_trace(go.Bar(
            x=years, y=nat_d.tolist(),
            marker_color="rgba(136, 136, 136, 0.7)",
            name="Natural Deaths",
            legendgroup="natd", showlegend=show_legend,
            hovertemplate="Natural: %{y}<extra></extra>",
            width=0.6,
        ), row=row, col=1, secondary_y=False)

        # Disease deaths (stacked)
        fig.add_trace(go.Bar(
            x=years, y=dis_d.tolist(),
            marker_color="rgba(228, 26, 28, 0.8)",
            name="Disease Deaths",
            legendgroup="disd", showlegend=show_legend,
            hovertemplate="Disease: %{y}<extra></extra>",
            width=0.6,
        ), row=row, col=1, secondary_y=False)

        # Vibrio on secondary y
        fig.add_trace(go.Scatter(
            x=years, y=vibrio.tolist(),
            mode="lines+markers",
            line=dict(color="#ff7f0e", width=2),
            marker=dict(size=4),
            name="Peak Vibrio",
            legendgroup="vibrio", showlegend=show_legend,
            hovertemplate="Vibrio: %{y:.0f} bact/mL<extra></extra>",
        ), row=row, col=1, secondary_y=True)

        # Disease marker
        fig.add_vline(x=disease_year, line_dash="dash", line_color="red",
                      line_width=1.5, opacity=0.5, row=row, col=1)

    fig.update_layout(
        title=dict(text="Panel 3: Disease Dynamics by Node", font=dict(size=18)),
        height=250 * n_nodes + 80,
        barmode="stack",
        legend=dict(orientation="h", yanchor="bottom", y=1.01, xanchor="center", x=0.5),
        margin=dict(l=60, r=60, t=100, b=40),
    )
    fig.update_xaxes(title_text="Year", row=n_nodes, col=1)

    return fig


def build_genetics_panel(data: dict) -> go.Figure:
    """Panel 4: Genetics summary (4 subplots + allele freq)."""
    meta = data["meta"]
    n_nodes = meta["n_nodes"]
    n_years = meta["n_years"]
    names = meta["node_names"]
    disease_year = meta["disease_year"]
    years = list(range(n_years))

    fig = make_subplots(
        rows=3, cols=2,
        shared_xaxes=True,
        subplot_titles=[
            "Mean Resistance (r̄)", "EF1A Allele Frequency",
            "Ne / N Ratio", "Additive Genetic Variance (V_A)",
            "Top-3 Loci Allele Frequencies", "Heterozygosity Proxy (V_A × 2)",
        ],
        vertical_spacing=0.08,
        horizontal_spacing=0.08,
    )

    for i in range(n_nodes):
        show_legend = (i == 0)
        color = NODE_COLORS[i]

        # (1,1) Mean resistance
        fig.add_trace(go.Scatter(
            x=years, y=data["yearly_mean_resistance"][i].tolist(),
            mode="lines", line=dict(color=color, width=2),
            name=names[i], legendgroup=names[i], showlegend=show_legend,
            hovertemplate=f"{names[i]}: " + "%{y:.4f}<extra></extra>",
        ), row=1, col=1)

        # (1,2) EF1A
        fig.add_trace(go.Scatter(
            x=years, y=data["yearly_ef1a_freq"][i].tolist(),
            mode="lines", line=dict(color=color, width=2),
            name=names[i], legendgroup=names[i], showlegend=False,
        ), row=1, col=2)

        # (2,1) Ne/N
        fig.add_trace(go.Scatter(
            x=years, y=data["yearly_ne_ratio"][i].tolist(),
            mode="lines", line=dict(color=color, width=2),
            name=names[i], legendgroup=names[i], showlegend=False,
        ), row=2, col=1)

        # (2,2) Va
        fig.add_trace(go.Scatter(
            x=years, y=data["yearly_va"][i].tolist(),
            mode="lines", line=dict(color=color, width=2),
            name=names[i], legendgroup=names[i], showlegend=False,
        ), row=2, col=2)

        # (3,1) Top-3 loci (dashed lines per locus)
        freq_top3 = data["yearly_allele_freq_top3"][i]  # (n_years, 3)
        dashes = ["solid", "dash", "dot"]
        for j in range(3):
            fig.add_trace(go.Scatter(
                x=years, y=freq_top3[:, j].tolist(),
                mode="lines",
                line=dict(color=color, width=1.5, dash=dashes[j]),
                name=f"{names[i]} L{j+1}" if i == 0 else None,
                legendgroup=f"locus{j}", showlegend=(i == 0),
                hovertemplate=f"L{j+1}: " + "%{y:.4f}<extra></extra>",
            ), row=3, col=1)

        # (3,2) Heterozygosity proxy: 2 × Va
        het_proxy = (data["yearly_va"][i] * 2).tolist()
        fig.add_trace(go.Scatter(
            x=years, y=het_proxy,
            mode="lines", line=dict(color=color, width=2),
            name=names[i], legendgroup=names[i], showlegend=False,
        ), row=3, col=2)

    # Disease year markers
    for row in range(1, 4):
        for col in range(1, 3):
            fig.add_vline(x=disease_year, line_dash="dash", line_color="red",
                          line_width=1, opacity=0.5, row=row, col=col)

    fig.update_layout(
        title=dict(text="Panel 4: Genetics", font=dict(size=18)),
        height=900,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="center", x=0.5),
        margin=dict(l=60, r=20, t=120, b=40),
    )
    fig.update_xaxes(title_text="Year", row=3, col=1)
    fig.update_xaxes(title_text="Year", row=3, col=2)

    return fig


def build_heatmap_panel(data: dict) -> go.Figure:
    """Panel 5: Genetic heatmap — Δq across 52 loci for all nodes."""
    meta = data["meta"]
    n_nodes = meta["n_nodes"]
    names = meta["node_names"]

    # Build matrix: rows = nodes, columns = loci
    delta_q = data["post_epidemic_allele_freq"] - data["pre_epidemic_allele_freq"]
    # delta_q shape: (n_nodes, 52)

    locus_labels = [f"L{i}" for i in range(51)] + ["EF1A"]

    fig = go.Figure(data=go.Heatmap(
        z=delta_q.tolist(),
        x=locus_labels,
        y=names,
        colorscale="RdBu_r",
        zmid=0,
        zmin=-0.15, zmax=0.15,
        colorbar=dict(title=dict(text="Δq", side="right")),
        hovertemplate="<b>%{y}</b><br>%{x}<br>Δq: %{z:.4f}<extra></extra>",
    ))

    # Highlight EF1A column
    fig.add_vline(x=51, line_width=3, line_color="gold", opacity=0.8)
    fig.add_annotation(x=51, y=-0.3, text="EF1A ★", showarrow=False,
                       font=dict(size=10, color="goldenrod"), yref="paper")

    fig.update_layout(
        title=dict(text="Panel 5: Allele Frequency Shifts (Δq = post − pre epidemic)",
                   font=dict(size=18)),
        height=350,
        xaxis=dict(title="Locus", tickangle=0, dtick=5),
        yaxis=dict(title=""),
        margin=dict(l=120, r=20, t=60, b=60),
    )

    return fig


def build_summary_table(data: dict) -> str:
    """Panel 6: Summary statistics as an HTML table."""
    meta = data["meta"]
    n_nodes = meta["n_nodes"]
    names = meta["node_names"]
    disease_year = meta["disease_year"]

    rows_html = []
    for i in range(n_nodes):
        K = int(data["node_K"][i])
        pre_pop = int(data["yearly_pop"][i, disease_year])
        final_pop = int(data["yearly_pop"][i, -1])
        min_pop = int(data["yearly_pop"][i, disease_year:].min())
        total_dis_deaths = int(data["yearly_disease_deaths"][i].sum())
        total_nat_deaths = int(data["yearly_natural_deaths"][i].sum())
        peak_prev = float(data["peak_disease_prevalence"][i])
        r_pre = float(data["yearly_mean_resistance"][i, disease_year])
        r_post = float(data["yearly_mean_resistance"][i, -1])
        delta_r = r_post - r_pre
        ef1a_pre = float(data["yearly_ef1a_freq"][i, disease_year])
        ef1a_post = float(data["yearly_ef1a_freq"][i, -1])
        crash_pct = (1 - min_pop / max(pre_pop, 1)) * 100

        # Mortality color
        if crash_pct >= 95:
            crash_color = "#e41a1c"
        elif crash_pct >= 80:
            crash_color = "#ff7f0e"
        else:
            crash_color = "#2ca02c"

        # Delta r color
        dr_color = "#2ca02c" if delta_r > 0 else "#e41a1c"

        rows_html.append(f"""
        <tr>
          <td style="font-weight:bold; color:{NODE_COLORS[i]}">{names[i]}</td>
          <td>{K}</td>
          <td>{pre_pop}</td>
          <td>{final_pop}</td>
          <td>{min_pop}</td>
          <td style="color:{crash_color}; font-weight:bold">{crash_pct:.1f}%</td>
          <td>{total_dis_deaths}</td>
          <td>{total_nat_deaths}</td>
          <td>{peak_prev:.1%}</td>
          <td>{r_pre:.4f}</td>
          <td>{r_post:.4f}</td>
          <td style="color:{dr_color}">{delta_r:+.4f}</td>
          <td>{ef1a_pre:.3f}</td>
          <td>{ef1a_post:.3f}</td>
        </tr>""")

    # Totals row
    total_pre = int(data["yearly_total_pop"][disease_year])
    total_final = int(data["yearly_total_pop"][-1])
    total_dis = int(data["yearly_disease_deaths"].sum())
    total_nat = int(data["yearly_natural_deaths"].sum())
    total_crash = (1 - total_final / max(total_pre, 1)) * 100

    rows_html.append(f"""
    <tr style="border-top: 2px solid #333; font-weight: bold; background: #f0f0f0">
      <td>TOTAL</td>
      <td>{int(data['node_K'].sum())}</td>
      <td>{total_pre}</td>
      <td>{total_final}</td>
      <td>—</td>
      <td>{total_crash:.1f}%</td>
      <td>{total_dis}</td>
      <td>{total_nat}</td>
      <td>—</td>
      <td colspan="4">—</td>
      <td>—</td>
    </tr>""")

    return f"""
    <div style="padding: 20px;">
      <h2 style="font-family: Arial, sans-serif; color: #333;">
        Panel 6: Summary Statistics
      </h2>
      <p style="font-family: Arial, sans-serif; color: #666; font-size: 14px;">
        Simulation: {meta['n_years']} years, {n_nodes} nodes |
        Disease introduced: year {disease_year} |
        Seed: {meta['seed']}
      </p>
      <table style="border-collapse: collapse; width: 100%; font-family: Arial, sans-serif; font-size: 13px;">
        <thead>
          <tr style="background: #333; color: white;">
            <th style="padding: 8px; text-align: left;">Node</th>
            <th style="padding: 8px;">K</th>
            <th style="padding: 8px;">Pre-N</th>
            <th style="padding: 8px;">Final N</th>
            <th style="padding: 8px;">Min N</th>
            <th style="padding: 8px;">Crash %</th>
            <th style="padding: 8px;">Dis. Deaths</th>
            <th style="padding: 8px;">Nat. Deaths</th>
            <th style="padding: 8px;">Peak Prev.</th>
            <th style="padding: 8px;">r̄ pre</th>
            <th style="padding: 8px;">r̄ final</th>
            <th style="padding: 8px;">Δr̄</th>
            <th style="padding: 8px;">EF1A pre</th>
            <th style="padding: 8px;">EF1A final</th>
          </tr>
        </thead>
        <tbody>
          {"".join(rows_html)}
        </tbody>
      </table>
    </div>
    """


# ═══════════════════════════════════════════════════════════════════════
# DASHBOARD ASSEMBLY
# ═══════════════════════════════════════════════════════════════════════

def generate_dashboard(result_dir: str | Path, output_path: str | Path) -> Path:
    """Generate self-contained HTML dashboard.

    Args:
        result_dir: Directory containing simulation_data.npz + metadata.json.
        output_path: Output HTML file path.

    Returns:
        Path to the generated HTML file.
    """
    result_dir = Path(result_dir)
    output_path = Path(output_path)

    print(f"Loading data from {result_dir}...")
    data = load_data(result_dir)
    meta = data["meta"]

    print("Building Panel 1: Network Overview...")
    fig_network = build_network_panel(data)

    print("Building Panel 2: Population Dynamics...")
    fig_pop = build_population_panel(data)

    print("Building Panel 3: Disease Dynamics...")
    fig_disease = build_disease_panel(data)

    print("Building Panel 4: Genetics...")
    fig_genetics = build_genetics_panel(data)

    print("Building Panel 5: Genetic Heatmap...")
    fig_heatmap = build_heatmap_panel(data)

    print("Building Panel 6: Summary Table...")
    summary_html = build_summary_table(data)

    # Convert figures to HTML divs
    print("Assembling dashboard...")
    network_html = fig_network.to_html(full_html=False, include_plotlyjs=False)
    pop_html = fig_pop.to_html(full_html=False, include_plotlyjs=False)
    disease_html = fig_disease.to_html(full_html=False, include_plotlyjs=False)
    genetics_html = fig_genetics.to_html(full_html=False, include_plotlyjs=False)
    heatmap_html = fig_heatmap.to_html(full_html=False, include_plotlyjs=False)

    # Year-by-year total population sparkline for header
    total_pop = data["yearly_total_pop"]

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>SSWD-EvoEpi Dashboard — {meta['n_years']}-Year Epidemic Simulation</title>
    <script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
    <style>
        * {{ box-sizing: border-box; }}
        body {{
            font-family: 'Segoe UI', Arial, sans-serif;
            background: #f5f7fa;
            margin: 0;
            padding: 0;
            color: #333;
        }}
        .header {{
            background: linear-gradient(135deg, #1a1a2e 0%, #16213e 50%, #0f3460 100%);
            color: white;
            padding: 30px 40px;
            text-align: center;
        }}
        .header h1 {{
            margin: 0 0 8px;
            font-size: 28px;
            letter-spacing: 1px;
        }}
        .header .subtitle {{
            font-size: 14px;
            opacity: 0.8;
            margin-bottom: 15px;
        }}
        .header .stats {{
            display: flex;
            justify-content: center;
            gap: 30px;
            flex-wrap: wrap;
        }}
        .header .stat {{
            text-align: center;
        }}
        .header .stat .value {{
            font-size: 28px;
            font-weight: bold;
            display: block;
        }}
        .header .stat .label {{
            font-size: 11px;
            opacity: 0.7;
            text-transform: uppercase;
        }}
        .panel {{
            background: white;
            margin: 20px;
            border-radius: 8px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
            overflow: hidden;
        }}
        .panel-content {{
            padding: 10px 20px 20px;
        }}
        .nav {{
            background: white;
            padding: 10px 20px;
            position: sticky;
            top: 0;
            z-index: 100;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            display: flex;
            gap: 10px;
            flex-wrap: wrap;
            justify-content: center;
        }}
        .nav a {{
            text-decoration: none;
            color: #0f3460;
            padding: 6px 16px;
            border-radius: 20px;
            font-size: 13px;
            font-weight: 600;
            transition: all 0.2s;
            border: 1px solid #ddd;
        }}
        .nav a:hover {{
            background: #0f3460;
            color: white;
        }}
        .footer {{
            text-align: center;
            padding: 20px;
            color: #999;
            font-size: 12px;
        }}
        .quote {{
            font-style: italic;
            color: rgba(255,255,255,0.6);
            margin-top: 10px;
            font-size: 13px;
        }}
    </style>
</head>
<body>
    <div class="header">
        <h1>🌊 SSWD-EvoEpi Simulation Dashboard</h1>
        <div class="subtitle">
            Coupled Eco-Evolutionary Epidemiological Agent-Based Model
        </div>
        <div class="stats">
            <div class="stat">
                <span class="value">{meta['n_years']}</span>
                <span class="label">Years</span>
            </div>
            <div class="stat">
                <span class="value">{meta['n_nodes']}</span>
                <span class="label">Nodes</span>
            </div>
            <div class="stat">
                <span class="value">{meta['initial_total_pop']:,}</span>
                <span class="label">Initial Pop</span>
            </div>
            <div class="stat">
                <span class="value" style="color: #ff6b6b">{meta['final_total_pop']:,}</span>
                <span class="label">Final Pop</span>
            </div>
            <div class="stat">
                <span class="value" style="color: #ffd93d">Year {meta['disease_year']}</span>
                <span class="label">Disease Onset</span>
            </div>
            <div class="stat">
                <span class="value">{int(data['yearly_disease_deaths'].sum()):,}</span>
                <span class="label">Disease Deaths</span>
            </div>
        </div>
        <div class="quote">"If I cannot see it, I cannot trust it."</div>
    </div>

    <div class="nav">
        <a href="#panel1">🗺️ Network</a>
        <a href="#panel2">📊 Population</a>
        <a href="#panel3">🦠 Disease</a>
        <a href="#panel4">🧬 Genetics</a>
        <a href="#panel5">🔬 Heatmap</a>
        <a href="#panel6">📋 Summary</a>
    </div>

    <div id="panel1" class="panel">
        <div class="panel-content">
            {network_html}
        </div>
    </div>

    <div id="panel2" class="panel">
        <div class="panel-content">
            {pop_html}
        </div>
    </div>

    <div id="panel3" class="panel">
        <div class="panel-content">
            {disease_html}
        </div>
    </div>

    <div id="panel4" class="panel">
        <div class="panel-content">
            {genetics_html}
        </div>
    </div>

    <div id="panel5" class="panel">
        <div class="panel-content">
            {heatmap_html}
        </div>
    </div>

    <div id="panel6" class="panel">
        {summary_html}
    </div>

    <div class="footer">
        SSWD-EvoEpi v0.1 | Seed: {meta['seed']} |
        Generated by viz.dashboard | Weertman Lab, UW FHL
    </div>
</body>
</html>"""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(html)
    print(f"\n✅ Dashboard written to {output_path}")
    print(f"   File size: {output_path.stat().st_size / 1024:.0f} KB")

    return output_path


# ═══════════════════════════════════════════════════════════════════════
# STATIC PNG EXPORT
# ═══════════════════════════════════════════════════════════════════════

def export_static_pngs(result_dir: str | Path, output_dir: str | Path) -> List[Path]:
    """Export static PNG plots for each panel.

    Uses matplotlib (plot_utils) for publication-quality static output.

    Args:
        result_dir: Directory with simulation data.
        output_dir: Directory for PNG output.

    Returns:
        List of generated PNG file paths.
    """
    from viz.plot_utils import (
        SimulationData,
        plot_population_all_nodes,
        plot_disease_all_nodes,
        plot_genetics_summary,
        plot_genetic_heatmap_all_nodes,
        plot_network_state,
    )
    import matplotlib.pyplot as plt

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    data = SimulationData.load(result_dir)

    outputs = []

    # Network state at 3 key years
    for yr in [0, data.disease_year, data.n_years - 1]:
        path = output_dir / f"network_year{yr}.png"
        fig = plot_network_state(data, yr)
        fig.savefig(path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        outputs.append(path)
        print(f"  ✓ {path.name}")

    # Population
    path = output_dir / "population_all_nodes.png"
    fig = plot_population_all_nodes(data)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    outputs.append(path)
    print(f"  ✓ {path.name}")

    # Disease
    path = output_dir / "disease_all_nodes.png"
    fig = plot_disease_all_nodes(data)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    outputs.append(path)
    print(f"  ✓ {path.name}")

    # Genetics
    path = output_dir / "genetics_summary.png"
    fig = plot_genetics_summary(data)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    outputs.append(path)
    print(f"  ✓ {path.name}")

    # Heatmap
    path = output_dir / "genetic_heatmap.png"
    fig = plot_genetic_heatmap_all_nodes(data)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    outputs.append(path)
    print(f"  ✓ {path.name}")

    print(f"\n✅ {len(outputs)} static PNGs exported to {output_dir}/")
    return outputs


# ═══════════════════════════════════════════════════════════════════════
# CLI ENTRY POINT
# ═══════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Generate SSWD-EvoEpi interactive dashboard from simulation results.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python -m viz.dashboard results/5node_epidemic_20yr/
  python -m viz.dashboard results/5node_epidemic_20yr/ --output my_dashboard.html
  python -m viz.dashboard results/5node_epidemic_20yr/ --png
        """,
    )
    parser.add_argument("result_dir", help="Directory with simulation_data.npz + metadata.json")
    parser.add_argument("--output", "-o", default=None,
                        help="Output HTML path (default: <result_dir>/dashboard.html)")
    parser.add_argument("--png", action="store_true",
                        help="Also export static PNG plots")
    parser.add_argument("--png-dir", default=None,
                        help="Directory for PNG output (default: <result_dir>/plots/)")

    args = parser.parse_args()

    result_dir = Path(args.result_dir)
    if not result_dir.exists():
        print(f"Error: {result_dir} does not exist", file=sys.stderr)
        sys.exit(1)

    output = Path(args.output) if args.output else result_dir / "dashboard.html"

    generate_dashboard(result_dir, output)

    if args.png:
        png_dir = Path(args.png_dir) if args.png_dir else result_dir / "plots"
        print(f"\nExporting static PNGs to {png_dir}/...")
        export_static_pngs(result_dir, png_dir)


if __name__ == "__main__":
    main()
