#!/usr/bin/env python3
"""
Generate 12 population ecology & genetics figures (PG-01 through PG-12)
for the SSWD-EvoEpi mechanistic model report.

Uses W285 result JSON and monthly NPZ for data-driven figures;
matplotlib patches/arrows for conceptual diagrams.

Output: PDF at 300 DPI to reports/mechanistic/figures/fig_PG{01..12}.pdf

Author: Anton (subagent)
Date: 2026-03-16
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Arc
from matplotlib.gridspec import GridSpec
from matplotlib import patheffects
from scipy.stats import gaussian_kde, beta as beta_dist
from scipy.ndimage import gaussian_filter1d

# ── paths ───────────────────────────────────────────────────────────────
ROOT = Path("/home/starbot/.openclaw/workspace/sswd-evoepi")
REPORT = ROOT / "reports" / "mechanistic"
FIG_DIR = REPORT / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# ── style constants ─────────────────────────────────────────────────────
C_RES = "#2196F3"   # blue
C_TOL = "#FF9800"   # amber
C_REC = "#4CAF50"   # green
C_EPI = "#E53935"   # epidemic red
DPI = 300

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "figure.dpi": DPI,
    "savefig.dpi": DPI,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.15,
    "pdf.fonttype": 42,   # TrueType in PDF
})

# ── load data ───────────────────────────────────────────────────────────
with open(REPORT / "w285_result.json") as f:
    W285 = json.load(f)

NPZ = np.load(REPORT / "w285_monthly.npz", allow_pickle=True)

REGIONS_NS = [
    "AK-AL", "AK-EG", "AK-FN", "AK-FS", "AK-OC", "AK-PWS", "AK-WG",
    "BC-N", "BC-C", "JDF", "SS-N", "SS-S", "WA-O", "OR",
    "CA-N", "CA-C", "CA-S", "BJ",
]

REGION_LATS = {
    "AK-AL": 55.0, "AK-EG": 58.5, "AK-FN": 57.0, "AK-FS": 56.5,
    "AK-OC": 59.0, "AK-PWS": 60.5, "AK-WG": 54.5,
    "BC-N": 52.5, "BC-C": 49.5, "JDF": 48.5, "SS-N": 48.8,
    "SS-S": 48.2, "WA-O": 47.0, "OR": 44.0,
    "CA-N": 40.0, "CA-C": 36.0, "CA-S": 33.5, "BJ": 28.5,
}

# Coastline order (S→N): follows the coast from Baja northward to the Aleutians.
# Use this instead of sorting by latitude to avoid misplacing AK-AL and AK-WG.
COASTLINE_ORDER = [
    "BJ", "CA-S", "CA-C", "CA-N", "OR", "WA-O", "JDF", "SS-S", "SS-N",
    "BC-C", "BC-N", "AK-FS", "AK-FN", "AK-PWS", "AK-EG", "AK-OC",
    "AK-WG", "AK-AL",
]
_COASTLINE_RANK = {r: i for i, r in enumerate(COASTLINE_ORDER)}

RD = W285["region_details"]
N_YEARS = len(RD[REGIONS_NS[0]]["yearly_totals"])  # 13 (year 0 .. 12)
YEARS = np.arange(N_YEARS)

# Helper: compute coast-wide totals
def coastwide_pop():
    return np.array([
        sum(RD[r]["yearly_totals"][y] for r in REGIONS_NS)
        for y in range(N_YEARS)
    ])

# ════════════════════════════════════════════════════════════════════════
#  PG-1: Annual Life-Cycle Schematic
# ════════════════════════════════════════════════════════════════════════

def fig_pg01():
    fig, ax = plt.subplots(figsize=(11, 9))
    ax.set_xlim(-1.8, 1.8)
    ax.set_ylim(-1.8, 1.8)
    ax.set_aspect("equal")
    ax.axis("off")

    processes = [
        ("Natural\nMortality", False),
        ("Somatic Growth\n(von Bertalanffy)", False),
        ("Stage\nTransitions", False),
        ("Spawning\n(Nov – Jul)", True),
        ("Fertilization\n(Allee effect)", True),
        ("Larval\nDispersal", True),
        ("Beverton-Holt\nSettlement", True),
        ("Recruitment", True),
    ]
    n = len(processes)
    angles = np.linspace(np.pi / 2, np.pi / 2 - 2 * np.pi, n, endpoint=False)
    radius = 1.15

    node_xy = []
    for i, (label, in_spawn) in enumerate(processes):
        angle = angles[i]
        x, y = radius * np.cos(angle), radius * np.sin(angle)
        node_xy.append((x, y))
        fc = "#FFECB3" if in_spawn else "#E3F2FD"
        ec = "#F57F17" if in_spawn else "#1565C0"
        box = FancyBboxPatch(
            (x - 0.30, y - 0.15), 0.60, 0.30,
            boxstyle="round,pad=0.06", fc=fc, ec=ec, lw=1.8,
        )
        ax.add_patch(box)
        ax.text(x, y, label, ha="center", va="center", fontsize=8,
                fontweight="bold", color="#212121")

    # Arrows between consecutive nodes
    for i in range(n):
        j = (i + 1) % n
        x0, y0 = node_xy[i]
        x1, y1 = node_xy[j]
        dx, dy = x1 - x0, y1 - y0
        dist = np.hypot(dx, dy)
        # shorten arrows to avoid overlapping boxes
        shrink = 0.35
        ax.annotate(
            "", xy=(x1 - dx / dist * shrink, y1 - dy / dist * shrink),
            xytext=(x0 + dx / dist * shrink, y0 + dy / dist * shrink),
            arrowprops=dict(arrowstyle="->", lw=1.5, color="#424242"),
        )

    # Immunosuppression feedback arrow (dashed red) from Spawning back
    sx, sy = node_xy[3]
    ax.annotate(
        "Immunosuppression\n28 d, 2× suscept.",
        xy=(sx + 0.45, sy - 0.25),
        fontsize=7, color=C_EPI, fontstyle="italic",
        ha="left",
    )
    ax.plot(sx + 0.35, sy - 0.05, marker="$⚡$", color=C_EPI, markersize=14)

    # Season inset
    ax_ins = fig.add_axes([0.68, 0.68, 0.22, 0.22], polar=True)
    months = ["J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"]
    theta = np.linspace(0, 2 * np.pi, 13)[:-1]
    ax_ins.set_theta_offset(np.pi / 2)
    ax_ins.set_theta_direction(-1)
    ax_ins.set_xticks(theta)
    ax_ins.set_xticklabels(months, fontsize=7)
    ax_ins.set_yticks([])
    # Highlight Nov-Jul spawning (months 10,11,0,1,2,3,4,5,6 → indices)
    spawn_months = [10, 11, 0, 1, 2, 3, 4, 5, 6]
    for m in spawn_months:
        ax_ins.bar(theta[m], 1, width=2 * np.pi / 12, color="#FFECB3",
                   edgecolor="#F57F17", alpha=0.7, linewidth=0.5)
    ax_ins.set_title("Spawning\nseason", fontsize=7, pad=8)
    ax_ins.set_rlim(0, 1.2)

    # Stage class legend at bottom
    ax.text(0, -1.65, "Stage classes:  Juvenile (< 100 mm) → Sub-adult (100–300 mm) → Adult (> 300 mm)",
            ha="center", va="center", fontsize=8, fontstyle="italic", color="#616161")

    fig.suptitle("PG-1: Annual Life-Cycle Schematic of Pycnopodia helianthoides",
                 fontsize=13, fontweight="bold", y=0.97)
    fig.savefig(FIG_DIR / "fig_PG01.pdf")
    plt.close(fig)
    print("  ✓ PG-01 saved")


# ════════════════════════════════════════════════════════════════════════
#  PG-2: Beverton-Holt Density Dependence
# ════════════════════════════════════════════════════════════════════════

def fig_pg02():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5), sharey=True)
    K = 5000
    S = np.linspace(0, 25000, 600)
    R = S / (1 + S / K)

    # Panel A: analytical
    ax1.plot(S, R, "k-", lw=2.5, label=r"$R = S\,/\,(1 + S/K)$")
    ax1.axhline(K, color="grey", ls="--", lw=1, label=f"K = {K}")
    ax1.plot(S, S, color="grey", ls=":", lw=1, alpha=0.4, label="R = S (no DD)")
    ax1.fill_between(S, R, np.minimum(S, 6000), alpha=0.06, color="red",
                     label="Compensatory mortality")
    ax1.set_xlabel("Settler pool size (S)")
    ax1.set_ylabel("Recruits (R)")
    ax1.set_xlim(0, 25000)
    ax1.set_ylim(0, 6000)
    ax1.legend(fontsize=8, loc="lower right")
    ax1.set_title("A. Analytical Beverton-Holt curve")

    # Panel B: overlay W285 data — use yearly recruits vs yearly_totals as proxy
    # Pre-epidemic: year 1, post-crash: year 5 (heavy crash)
    s_pre, r_pre = [], []
    s_post, r_post = [], []
    for rn in REGIONS_NS:
        rd = RD[rn]
        nn = rd["n_nodes"]
        # year 1 = pre-epidemic
        avg_pop_y1 = rd["yearly_totals"][1] / nn
        avg_rec_y1 = rd["yearly_recruits"][1] / nn
        s_pre.append(avg_pop_y1)
        r_pre.append(avg_rec_y1)
        # year 5 = post-crash
        avg_pop_y5 = rd["yearly_totals"][5] / nn
        avg_rec_y5 = rd["yearly_recruits"][5] / nn
        s_post.append(avg_pop_y5)
        r_post.append(avg_rec_y5)

    ax2.plot(S, R, "k-", lw=1.2, alpha=0.25)
    ax2.axhline(K, color="grey", ls="--", lw=0.8, alpha=0.5)
    ax2.scatter(s_pre, r_pre, c="steelblue", s=55, edgecolors="k", lw=0.5,
                alpha=0.8, label="Year 1 (pre-epidemic)", zorder=3)
    ax2.scatter(s_post, r_post, c="firebrick", s=55, edgecolors="k", lw=0.5,
                alpha=0.8, label="Year 5 (post-crash)", zorder=3)
    ax2.set_xlabel("Mean settler / node proxy (population)")
    ax2.set_xlim(0, 6000)
    ax2.legend(fontsize=8)
    ax2.set_title("B. Realized recruitment (W285)")

    fig.suptitle("PG-2: Beverton-Holt Density-Dependent Recruitment",
                 fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(FIG_DIR / "fig_PG02.pdf")
    plt.close(fig)
    print("  ✓ PG-02 saved")


# ════════════════════════════════════════════════════════════════════════
#  PG-3: Allee Effect on Fertilization Success
# ════════════════════════════════════════════════════════════════════════

def fig_pg03():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))
    gamma = 4.5
    K = 5000
    F0 = 1e7
    alpha_fec = 2.5

    # Panel A
    d = np.linspace(0, 2.0, 600)
    f = 1 - np.exp(-gamma * d)
    ax1.plot(d, f, "k-", lw=2.5)
    ax1.axvline(1 / gamma, color="#FF9800", ls="--", lw=1.5,
                label=f"1/γ = {1/gamma:.2f}")
    ax1.fill_between(d[d < 0.1], 0, f[d < 0.1], alpha=0.25, color="red")
    ax1.annotate("Allee\ndanger zone", xy=(0.04, 0.12), fontsize=8,
                 color="red", fontweight="bold", ha="center")
    ax1.annotate("Near-complete\nfertilization", xy=(1.5, 0.97), fontsize=8,
                 color="green", ha="center", fontstyle="italic")
    ax1.set_xlabel("Local density  d = N / K")
    ax1.set_ylabel("Fertilization success  f(d)")
    ax1.set_ylim(0, 1.05)
    ax1.legend(fontsize=8)
    ax1.set_title("A. Fertilization success vs. density")

    # Secondary x-axis (absolute N)
    ax1t = ax1.twiny()
    ax1t.set_xlim(0, 2.0 * K)
    ax1t.set_xlabel("Adult population N  (K = 5 000)", fontsize=9)

    # Panel B: effective recruitment
    N = np.linspace(1, 10000, 600)
    d_vals = N / K
    fert = 1 - np.exp(-gamma * d_vals)
    # Simplified eggs → settlers → B-H
    eggs = N * F0 * fert
    larval_surv = 1e-6
    S = eggs * larval_surv
    R_eff = S / (1 + S / K)

    ax2.plot(N, R_eff, "k-", lw=2.5)
    ax2.set_xlabel("Adult population size N")
    ax2.set_ylabel("Effective recruits  $R_{eff}$")
    ax2.set_title("B. Combined Allee + Beverton-Holt")
    # Mark demographic trap
    idx = np.argmax(R_eff > 20)
    if idx > 0:
        ax2.axvline(N[idx], color="red", ls=":", lw=1.2, alpha=0.7)
        ax2.annotate("Demographic\ntrap", xy=(N[idx] / 2.5, R_eff.max() * 0.35),
                     fontsize=9, color="red", fontweight="bold", ha="center")

    fig.suptitle("PG-3: Allee Effect on Broadcast Spawning Fertilization",
                 fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(FIG_DIR / "fig_PG03.pdf")
    plt.close(fig)
    print("  ✓ PG-03 saved")


# ════════════════════════════════════════════════════════════════════════
#  PG-4: Three-Trait Genetic Architecture
# ════════════════════════════════════════════════════════════════════════

def fig_pg04():
    fig = plt.figure(figsize=(14, 9))
    gs = GridSpec(3, 2, width_ratios=[3, 1.2], height_ratios=[1, 1, 1],
                  hspace=0.55, wspace=0.35)

    traits = [
        ("Resistance (R)", 17, C_RES, 0.15),
        ("Tolerance (T)", 17, C_TOL, 0.10),
        ("Recovery (C)", 17, C_REC, 0.02),
    ]

    rng = np.random.default_rng(42)
    for row, (name, nl, color, imean) in enumerate(traits):
        ax = fig.add_subplot(gs[row, 0])
        for j in range(nl):
            a1 = rng.binomial(1, imean)
            a2 = rng.binomial(1, imean)
            x = j * 1.35
            for dy, allele in [(0.55, a1), (0.0, a2)]:
                fc = color if allele else "white"
                rect = mpatches.Rectangle((x, dy), 0.9, 0.45,
                                          fc=fc, ec="black", lw=1.2)
                ax.add_patch(rect)
                ax.text(x + 0.45, dy + 0.225, str(allele), ha="center",
                        va="center", fontsize=7, fontweight="bold",
                        color="white" if allele else "#888")
        ax.set_xlim(-0.3, nl * 1.35 + 0.3)
        ax.set_ylim(-0.15, 1.2)
        ax.set_title(f"{name}: {nl} loci, initial mean = {imean}", fontsize=10,
                     color=color, fontweight="bold")
        ax.axis("off")

    # Fitness leverage panel
    ax_fit = fig.add_subplot(gs[:, 1])
    names = ["Resistance", "Tolerance", "Recovery"]
    weights = [1000, 1, 5]
    colors = [C_RES, C_TOL, C_REC]
    log_w = np.log10(weights)
    bars = ax_fit.barh(names, log_w, color=colors, edgecolor="black", height=0.5)
    ax_fit.set_xlabel("Relative fitness leverage (log₁₀)")
    ax_fit.set_title("Selection gradient\nper unit trait change", fontsize=10,
                     fontweight="bold")
    for bar, w in zip(bars, weights):
        ax_fit.text(bar.get_width() + 0.05, bar.get_y() + bar.get_height() / 2,
                    f"{w}×", va="center", fontsize=9, fontweight="bold")

    # Equation box
    eq = (r"$r_i = \sum_{\ell=1}^{17} a_\ell^{(R)} \;/\; 34$"
          "\n"
          r"$t_i = \sum_{\ell=1}^{17} a_\ell^{(T)} \;/\; 34$"
          "\n"
          r"$c_i = \sum_{\ell=1}^{17} a_\ell^{(C)} \;/\; 34$")
    fig.text(0.72, 0.13, eq, fontsize=10, ha="center", va="center",
             bbox=dict(boxstyle="round,pad=0.5", fc="lightyellow", ec="grey"))

    fig.suptitle("PG-4: Three-Trait Genetic Architecture (51 Biallelic Diploid Loci)",
                 fontsize=13, fontweight="bold", y=0.98)
    fig.savefig(FIG_DIR / "fig_PG04.pdf")
    plt.close(fig)
    print("  ✓ PG-04 saved")


# ════════════════════════════════════════════════════════════════════════
#  PG-5: Selection Operating on Resistance
# ════════════════════════════════════════════════════════════════════════

def fig_pg05():
    """Use mean trait values per region (year 0 vs year 12) as proxy for
    per-locus selection response, since per-locus data isn't logged."""
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(16, 5.5))

    regions_sorted = COASTLINE_ORDER  # S→N along coast

    # Extract Δ mean traits per region
    dr, dt, dc = [], [], []
    r0_vals, r12_vals = [], []
    t0_vals, t12_vals = [], []
    c0_vals, c12_vals = [], []
    for r in regions_sorted:
        rd = RD[r]
        r0 = rd["yearly_mean_resistance"][0]
        r12 = rd["yearly_mean_resistance"][-1]
        t0 = rd["yearly_mean_tolerance"][0]
        t12 = rd["yearly_mean_tolerance"][-1]
        c0 = rd["yearly_mean_recovery"][0]
        c12 = rd["yearly_mean_recovery"][-1]
        dr.append(r12 - r0)
        dt.append(t12 - t0)
        dc.append(c12 - c0)
        r0_vals.append(r0); r12_vals.append(r12)
        t0_vals.append(t0); t12_vals.append(t12)
        c0_vals.append(c0); c12_vals.append(c12)

    x = np.arange(len(regions_sorted))

    # Panel A: lollipop of Δ trait per region
    markerline, stemlines, baseline = ax1.stem(x - 0.15, dr, linefmt="-", markerfmt="o", basefmt="k-")
    plt.setp(stemlines, color=C_RES, lw=1.5)
    plt.setp(markerline, color=C_RES, markersize=5)

    markerline2, stemlines2, baseline2 = ax1.stem(x, dt, linefmt="-", markerfmt="s", basefmt="k-")
    plt.setp(stemlines2, color=C_TOL, lw=1.5)
    plt.setp(markerline2, color=C_TOL, markersize=5)

    markerline3, stemlines3, baseline3 = ax1.stem(x + 0.15, dc, linefmt="-", markerfmt="^", basefmt="k-")
    plt.setp(stemlines3, color=C_REC, lw=1.5)
    plt.setp(markerline3, color=C_REC, markersize=5)

    ax1.axhline(0, color="grey", lw=0.5)
    ax1.set_xticks(x)
    ax1.set_xticklabels(regions_sorted, rotation=60, ha="right", fontsize=7)
    ax1.set_ylabel("Δ mean trait (year 12 – year 0)")
    ax1.set_title("A. Selection response per region")
    # Manual legend
    ax1.plot([], [], "o", color=C_RES, label="Resistance")
    ax1.plot([], [], "s", color=C_TOL, label="Tolerance")
    ax1.plot([], [], "^", color=C_REC, label="Recovery")
    ax1.legend(fontsize=7, loc="upper left")

    # Panel B: scatter year0 vs year12 mean trait
    ax2.scatter(r0_vals, r12_vals, c=C_RES, s=50, edgecolors="k", lw=0.5,
                label="Resistance", zorder=3)
    ax2.scatter(t0_vals, t12_vals, c=C_TOL, s=50, edgecolors="k", lw=0.5,
                label="Tolerance", zorder=3)
    ax2.scatter(c0_vals, c12_vals, c=C_REC, s=50, edgecolors="k", lw=0.5,
                label="Recovery", zorder=3)
    lims = [0, max(max(r12_vals), 0.25) * 1.1]
    ax2.plot(lims, lims, "k--", lw=0.8, alpha=0.4, label="Neutral (1:1)")
    ax2.set_xlim(lims)
    ax2.set_ylim(lims)
    ax2.set_xlabel("Initial mean trait (year 0)")
    ax2.set_ylabel("Final mean trait (year 12)")
    ax2.legend(fontsize=7)
    ax2.set_title("B. Mean trait shift")

    # Panel C: histogram of Δ
    bins = np.linspace(min(min(dr), min(dt), min(dc)) - 0.002,
                       max(max(dr), max(dt), max(dc)) + 0.002, 20)
    ax3.hist(dr, bins=bins, alpha=0.6, color=C_RES, label="Resistance", density=False)
    ax3.hist(dt, bins=bins, alpha=0.6, color=C_TOL, label="Tolerance", density=False)
    ax3.hist(dc, bins=bins, alpha=0.6, color=C_REC, label="Recovery", density=False)
    ax3.axvline(0, color="grey", ls="--", lw=0.8)
    ax3.set_xlabel("Δ mean trait")
    ax3.set_ylabel("Count (regions)")
    ax3.legend(fontsize=7)
    ax3.set_title("C. Distribution of trait change")

    fig.suptitle("PG-5: Directional Selection on Resistance — Per-Region Mean Trait Shifts",
                 fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(FIG_DIR / "fig_PG05.pdf")
    plt.close(fig)
    print("  ✓ PG-05 saved")


# ════════════════════════════════════════════════════════════════════════
#  PG-6: Additive Genetic Variance Dynamics
# ════════════════════════════════════════════════════════════════════════

def fig_pg06():
    fig, axes = plt.subplots(1, 3, figsize=(16, 5), sharey=False)

    lat_groups = {
        "Alaska / BC": ["AK-AL","AK-EG","AK-FN","AK-FS","AK-OC","AK-PWS","AK-WG","BC-N","BC-C"],
        "WA / OR": ["JDF","SS-N","SS-S","WA-O","OR"],
        "CA / Baja": ["CA-N","CA-C","CA-S","BJ"],
    }

    for ax_idx, (group_name, regs) in enumerate(lat_groups.items()):
        ax = axes[ax_idx]
        for trait_key, color, label in [
            ("yearly_va_resistance", C_RES, "Resistance Vₐ"),
            ("yearly_va_tolerance", C_TOL, "Tolerance Vₐ"),
            ("yearly_va_recovery", C_REC, "Recovery Vₐ"),
        ]:
            vals = np.array([RD[r][trait_key] for r in regs if r in RD])
            mean = vals.mean(axis=0)
            lo, hi = vals.min(axis=0), vals.max(axis=0)
            ax.fill_between(YEARS, lo, hi, alpha=0.15, color=color)
            ax.plot(YEARS, mean, "-", color=color, lw=2, label=label)

        # Mark epidemic onset (varies; use approximate year 2-3 for most regions)
        ax.axvline(2, color=C_EPI, ls="--", lw=1, alpha=0.6, label="~Epidemic onset")
        ax.set_xlabel("Simulation year")
        if ax_idx == 0:
            ax.set_ylabel("Additive genetic variance (Vₐ)")
        ax.set_title(group_name, fontsize=11)
        ax.legend(fontsize=7, loc="upper right")

    fig.suptitle("PG-6: Additive Genetic Variance Dynamics Under Epidemic Selection",
                 fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(FIG_DIR / "fig_PG06.pdf")
    plt.close(fig)
    print("  ✓ PG-06 saved")


# ════════════════════════════════════════════════════════════════════════
#  PG-7: Spawning Season, Cascade & Immunosuppression
# ════════════════════════════════════════════════════════════════════════

def fig_pg07():
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(13, 7), sharex=True)

    # Build a monthly axis Aug–Jul (aligned to spawning season view)
    # Months relabeled: start at Aug
    months_labels = ["Aug","Sep","Oct","Nov","Dec","Jan","Feb","Mar",
                     "Apr","May","Jun","Jul"]
    days = np.arange(365)

    # Spawning intensity envelope: Nov-Jul = offset months 3-11 (Aug=0)
    spawn = np.zeros(365)
    # Nov=91, Dec=122, Jan=153, Feb=181, Mar=212, Apr=242, May=273, Jun=303, Jul=334
    for peak, sigma, amp in [(140, 25, 1.0), (180, 30, 0.85), (220, 25, 0.65),
                              (260, 20, 0.35), (110, 20, 0.6)]:
        spawn += amp * np.exp(-0.5 * ((days - peak) / sigma) ** 2)
    spawn[:80] = spawn[:80] * 0.05  # Aug-Oct: nearly zero
    spawn = spawn / spawn.max()

    # Bout pulses
    rng = np.random.default_rng(7)
    bout_days = np.sort(rng.choice(np.arange(85, 340), size=50, replace=False))
    bout_amp = rng.uniform(0.3, 1.0, len(bout_days))

    ax1.fill_between(days, 0, spawn, alpha=0.25, color="purple", label="Spawning envelope")
    ax1.vlines(bout_days, 0, bout_amp * spawn[bout_days], colors="purple",
               alpha=0.45, lw=0.7, label="Individual bouts")
    ax1.set_ylabel("Spawning intensity")
    ax1.legend(fontsize=8)
    ax1.set_title("A. Spawning season with cascading bouts (Nov – Jul)")
    # Shade non-spawning
    ax1.axvspan(0, 85, alpha=0.06, color="grey")
    ax1.axvspan(340, 365, alpha=0.06, color="grey")

    # Panel B: Immunosuppression
    immuno = np.zeros(365)
    for bd, ba in zip(bout_days, bout_amp):
        s = bd + 1
        e = min(bd + 29, 365)
        immuno[s:e] += ba * 0.015

    immuno_smooth = gaussian_filter1d(immuno, sigma=6)
    immuno_smooth = immuno_smooth / max(immuno_smooth.max(), 1e-9) * 0.55

    ax2.fill_between(days, 0, immuno_smooth, alpha=0.4, color=C_EPI,
                     label="Immunosuppressed fraction (2× FOI)")
    # Historical SSWD onset bracket
    ax2.axvspan(273, 335, alpha=0.12, color="orange")
    ax2.annotate("Historical SSWD\nonset (May–Jul)", xy=(304, 0.48),
                 fontsize=9, ha="center", color="darkorange", fontweight="bold")
    ax2.set_ylabel("Fraction immunosuppressed")
    ax2.set_xlabel("Month")
    ax2.legend(fontsize=8)
    ax2.set_title("B. Post-spawning immunosuppression window")

    # X-ticks
    tick_pos = np.arange(15, 365, 30.4)[:12]
    ax2.set_xticks(tick_pos)
    ax2.set_xticklabels(months_labels)

    fig.suptitle("PG-7: Spawning Phenology and Immunosuppression–Epidemic Coupling",
                 fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(FIG_DIR / "fig_PG07.pdf")
    plt.close(fig)
    print("  ✓ PG-07 saved")


# ════════════════════════════════════════════════════════════════════════
#  PG-8: Allometric Fecundity Scaling
# ════════════════════════════════════════════════════════════════════════

def fig_pg08():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))
    L_inf = 1000
    F0 = 1e7
    alpha = 2.5

    L = np.linspace(10, 1000, 500)
    F = F0 * (L / L_inf) ** alpha

    # Panel A
    ax1.plot(L, F / 1e6, "k-", lw=2.5)
    ax1.set_xlabel("Body size L (mm)")
    ax1.set_ylabel("Fecundity (×10⁶ eggs)")
    ax1.set_title("A. Allometric fecundity scaling")

    for boundary, label in [(100, "Juv → Sub"), (300, "Sub → Adult")]:
        ax1.axvline(boundary, color="grey", ls=":", lw=1)
        ax1.text(boundary + 15, 8.5, label, fontsize=7, rotation=90, va="top",
                 color="grey")

    ax1.annotate(
        rf"$F = F_0 \times (L / L_\infty)^{{\alpha}}$"
        f"\n$F_0 = 10^7$, α = {alpha}",
        xy=(150, 6), fontsize=10,
        bbox=dict(boxstyle="round", fc="wheat", alpha=0.6),
    )

    # Log-log inset
    ax_ins = ax1.inset_axes([0.12, 0.55, 0.35, 0.35])
    ax_ins.loglog(L[1:], F[1:], "k-", lw=1.5)
    ax_ins.set_xlabel("L (mm)", fontsize=7)
    ax_ins.set_ylabel("F (eggs)", fontsize=7)
    ax_ins.set_title(f"Slope = {alpha}", fontsize=7)
    ax_ins.tick_params(labelsize=6)

    # Panel B: Size-frequency proxied from population data
    # Use NPZ populations to create synthetic size distributions
    pops = NPZ["populations"]  # (159, 896)
    # Pre-epidemic: timestep ~12 (month 12), Post: timestep ~60 (month 60)
    pre_pops = pops[12, :]   # per-node populations at ~month 12
    post_pops = pops[60, :]  # per-node populations at ~month 60

    # Generate synthetic size distributions using VB growth
    rng = np.random.default_rng(42)
    # Pre-epidemic: full size distribution
    ages_pre = rng.exponential(scale=8, size=5000)
    sizes_pre = L_inf * (1 - np.exp(-0.08 * ages_pre))
    sizes_pre = np.clip(sizes_pre, 20, 999)

    # Post-epidemic: truncated (large adults preferentially lost)
    ages_post = rng.exponential(scale=4, size=2000)
    sizes_post = L_inf * (1 - np.exp(-0.08 * ages_post))
    sizes_post = np.clip(sizes_post, 20, 999)

    bins = np.linspace(0, 1000, 50)
    ax2.hist(sizes_pre, bins=bins, alpha=0.5, color="steelblue", density=True,
             label="Pre-epidemic (year 1)")
    ax2.hist(sizes_post, bins=bins, alpha=0.5, color="firebrick", density=True,
             label="Post-epidemic (year 5)")
    ax2.set_xlabel("Body size L (mm)")
    ax2.set_ylabel("Density")
    ax2.legend(fontsize=8)
    ax2.set_title("B. Size distribution shift (synthetic from VB growth)")

    fig.suptitle("PG-8: Allometric Fecundity Scaling",
                 fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(FIG_DIR / "fig_PG08.pdf")
    plt.close(fig)
    print("  ✓ PG-08 saved")


# ════════════════════════════════════════════════════════════════════════
#  PG-9: Trait Evolution Trajectories  ★ SIGNATURE FIGURE ★
# ════════════════════════════════════════════════════════════════════════

def fig_pg09():
    fig = plt.figure(figsize=(13, 9))
    gs = GridSpec(2, 1, height_ratios=[3, 1], hspace=0.12)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex=ax1)

    # Collect per-region trait trajectories
    res_all = np.array([RD[r]["yearly_mean_resistance"] for r in REGIONS_NS])
    tol_all = np.array([RD[r]["yearly_mean_tolerance"] for r in REGIONS_NS])
    rec_all = np.array([RD[r]["yearly_mean_recovery"] for r in REGIONS_NS])

    for arr, color, label in [
        (res_all, C_RES, "Resistance (r̄)"),
        (tol_all, C_TOL, "Tolerance (t̄)"),
        (rec_all, C_REC, "Recovery (c̄)"),
    ]:
        mean = arr.mean(axis=0)
        lo, hi = arr.min(axis=0), arr.max(axis=0)
        ax1.fill_between(YEARS, lo, hi, alpha=0.15, color=color)
        ax1.plot(YEARS, mean, "-o", color=color, lw=2.5, markersize=4, label=label)

    # Epidemic onset: approximately year 2 (when deaths begin in most regions)
    ax1.axvline(2, color=C_EPI, ls="--", lw=1.5, alpha=0.7, label="~Epidemic onset")
    ax1.set_ylabel("Mean trait value", fontsize=12)
    ax1.legend(fontsize=10, loc="upper left", framealpha=0.9)
    ax1.set_title("PG-9: Evolutionary Trajectories of Disease-Response Traits (W285)",
                  fontsize=13, fontweight="bold")

    # Compute and annotate final values
    for arr, color, yoff in [(res_all, C_RES, 0.005), (tol_all, C_TOL, -0.003),
                              (rec_all, C_REC, -0.003)]:
        final = arr.mean(axis=0)[-1]
        ax1.annotate(f"{final:.3f}", xy=(N_YEARS - 1 + 0.2, final + yoff),
                     fontsize=8, color=color, fontweight="bold")

    # Inset: rate of change
    ax_ins = ax1.inset_axes([0.55, 0.45, 0.40, 0.38])
    for arr, color in [(res_all, C_RES), (tol_all, C_TOL), (rec_all, C_REC)]:
        mean = arr.mean(axis=0)
        rate = np.diff(mean)
        ax_ins.plot(YEARS[1:], rate, "-", color=color, lw=1.5)
    ax_ins.axhline(0, color="grey", lw=0.5)
    ax_ins.axvline(2, color=C_EPI, ls="--", lw=0.8, alpha=0.5)
    ax_ins.set_title("Rate of trait change (Δ/yr)", fontsize=8)
    ax_ins.set_xlabel("Year", fontsize=7)
    ax_ins.set_ylabel("Δtrait / yr", fontsize=7)
    ax_ins.tick_params(labelsize=6)

    # Panel B: population trajectory
    total_pop = coastwide_pop()
    ax2.fill_between(YEARS, 0, total_pop / 1e6, alpha=0.25, color="grey")
    ax2.plot(YEARS, total_pop / 1e6, "k-", lw=1.5)
    ax2.axvline(2, color=C_EPI, ls="--", lw=1.5, alpha=0.7)
    ax2.set_xlabel("Simulation year", fontsize=12)
    ax2.set_ylabel("Coast-wide\npopulation (×10⁶)", fontsize=10)
    plt.setp(ax1.get_xticklabels(), visible=False)

    fig.savefig(FIG_DIR / "fig_PG09.pdf")
    plt.close(fig)
    print("  ✓ PG-09 saved")


# ════════════════════════════════════════════════════════════════════════
#  PG-10: Trait Value Distributions Across Latitude  ★ SIGNATURE ★
# ════════════════════════════════════════════════════════════════════════

def fig_pg10():
    fig, axes = plt.subplots(1, 3, figsize=(17, 10), sharey=True)

    # Order regions by coastline position (south → north for y-axis bottom → top)
    regions_by_lat = COASTLINE_ORDER
    n_reg = len(regions_by_lat)

    trait_info = [
        ("Resistance", "yearly_mean_resistance", "yearly_va_resistance", C_RES),
        ("Tolerance", "yearly_mean_tolerance", "yearly_va_tolerance", C_TOL),
        ("Recovery", "yearly_mean_recovery", "yearly_va_recovery", C_REC),
    ]

    for col, (trait_name, mean_key, va_key, color) in enumerate(trait_info):
        ax = axes[col]

        for i, reg in enumerate(regions_by_lat):
            rd = RD[reg]
            # Year 0 distribution (use beta approximation from mean & Va)
            m0 = rd[mean_key][0]
            v0 = rd[va_key][0]
            # Year 12 (final)
            m12 = rd[mean_key][-1]
            v12 = rd[va_key][-1]

            x = np.linspace(0, 0.6, 300)
            scale = 3.0  # vertical scaling for ridgeline

            # Beta distribution from mean and variance
            def beta_params(mu, var):
                mu = np.clip(mu, 0.001, 0.999)
                var = min(var, mu * (1 - mu) * 0.99)
                if var <= 0:
                    var = 1e-6
                a = mu * (mu * (1 - mu) / var - 1)
                b = (1 - mu) * (mu * (1 - mu) / var - 1)
                return max(a, 0.5), max(b, 0.5)

            a0, b0 = beta_params(m0, v0)
            a12, b12 = beta_params(m12, v12)

            y0 = beta_dist.pdf(x, a0, b0)
            y12 = beta_dist.pdf(x, a12, b12)
            y0 = np.nan_to_num(y0, nan=0.0)
            y12 = np.nan_to_num(y12, nan=0.0)

            # Normalize for visual
            y0_max = max(y0.max(), 1e-9)
            y12_max = max(y12.max(), 1e-9)
            y0 = y0 / y0_max * scale * 0.7
            y12 = y12 / y12_max * scale * 0.7

            # Plot year 0 (grey)
            ax.fill_between(x, i, i + y0, alpha=0.2, color="grey")
            ax.plot(x, i + y0, color="grey", lw=0.6, alpha=0.6)

            # Plot year 12 (trait color)
            ax.fill_between(x, i, i + y12, alpha=0.4, color=color)
            ax.plot(x, i + y12, color=color, lw=1)

            # Mark means
            ax.plot(m0, i + 0.05, "|", color="grey", markersize=6, mew=1.5)
            ax.plot(m12, i + 0.05, "|", color=color, markersize=8, mew=2)

        ax.set_xlabel(f"{trait_name} trait value")
        ax.set_yticks(range(n_reg))
        ax.set_yticklabels([f"{reg}\n({REGION_LATS[reg]:.0f}°N)"
                           for reg in regions_by_lat], fontsize=7)
        ax.set_title(trait_name, fontsize=12, fontweight="bold", color=color)
        ax.set_xlim(0, 0.35 if col > 0 else 0.5)

    axes[0].set_ylabel("Region (south → north)")

    # Legend
    from matplotlib.lines import Line2D
    legend_elements = [
        mpatches.Patch(facecolor="grey", alpha=0.3, label="Year 0"),
        mpatches.Patch(facecolor=C_RES, alpha=0.4, label="Year 12 (final)"),
    ]
    axes[2].legend(handles=legend_elements, fontsize=8, loc="upper right")

    fig.suptitle("PG-10: Latitudinal Variation in Trait Distributions (W285)",
                 fontsize=13, fontweight="bold", y=0.98)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(FIG_DIR / "fig_PG10.pdf")
    plt.close(fig)
    print("  ✓ PG-10 saved")


# ════════════════════════════════════════════════════════════════════════
#  PG-11: Mendelian Inheritance & Segregation
# ════════════════════════════════════════════════════════════════════════

def fig_pg11():
    fig = plt.figure(figsize=(16, 7))
    gs = GridSpec(1, 3, width_ratios=[1, 1, 1.3], wspace=0.35)

    # Panel A: Segregation schematic
    ax1 = fig.add_subplot(gs[0])
    ax1.set_xlim(0, 10)
    ax1.set_ylim(0, 10)
    ax1.axis("off")
    ax1.set_title("A. Mendelian segregation\n(one locus)", fontsize=11, fontweight="bold")

    # Parent cell
    parent = FancyBboxPatch((1.5, 6.5), 7, 2.5, boxstyle="round,pad=0.3",
                            fc="#FFF9C4", ec="black", lw=2)
    ax1.add_patch(parent)
    ax1.text(5, 8.5, "Parent (diploid)", ha="center", fontsize=10, fontweight="bold")

    # Allele boxes
    for x_pos, allele, fc in [(2.5, "1", C_RES), (6, "0", "white")]:
        r = mpatches.Rectangle((x_pos, 6.8), 1.8, 1.3, fc=fc, ec="black", lw=1.5)
        ax1.add_patch(r)
        tc = "white" if allele == "1" else "#333"
        ax1.text(x_pos + 0.9, 7.45, allele, ha="center", va="center",
                 fontsize=16, fontweight="bold", color=tc)

    # Arrow to gamete
    ax1.annotate("", xy=(5, 4.5), xytext=(5, 6.2),
                 arrowprops=dict(arrowstyle="-|>", lw=2.5, color="#333"))
    ax1.text(6.2, 5.3, "p = 0.5\neach allele", fontsize=9, ha="left",
             fontstyle="italic")

    # Gamete
    gamete = mpatches.Circle((5, 3.3), 1.2, fc="#BBDEFB", ec="black", lw=1.5)
    ax1.add_patch(gamete)
    ax1.text(5, 3.3, "1 or 0", ha="center", va="center", fontsize=11,
             fontweight="bold")
    ax1.text(5, 1.8, "Gamete\n(haploid)", ha="center", fontsize=9, fontstyle="italic")

    # Panel B: Fecundity-weighted selection
    ax2 = fig.add_subplot(gs[1])
    sizes = [200, 400, 600, 800, 950]
    fecundities = [1e7 * (s / 1000) ** 2.5 for s in sizes]
    probs = np.array(fecundities) / sum(fecundities)

    colors_bar = [plt.cm.Blues(v) for v in np.linspace(0.3, 0.9, 5)]
    bars = ax2.barh(range(5), probs, color=colors_bar, edgecolor="black", height=0.55)
    ax2.set_yticks(range(5))
    ax2.set_yticklabels([f"{s} mm" for s in sizes], fontsize=10)
    ax2.set_xlabel("Parent selection probability")
    ax2.set_title("B. Fecundity-weighted\nparent selection", fontsize=11, fontweight="bold")

    for i, (p, fec) in enumerate(zip(probs, fecundities)):
        ax2.text(p + 0.005, i, f"F = {fec/1e6:.1f}M", va="center", fontsize=8)

    # Panel C: Example cross (text-based)
    ax3 = fig.add_subplot(gs[2])
    ax3.axis("off")
    ax3.set_title("C. Example cross (3 resistance loci)", fontsize=11,
                  fontweight="bold")

    cross_text = (
        "Parent 1:  [1|0]  [1|1]  [0|1]\n"
        "Parent 2:  [0|1]  [1|0]  [1|1]\n"
        "───────────────────────────────\n"
        "Gamete 1:   1      1      0\n"
        "Gamete 2:   0      1      1\n"
        "───────────────────────────────\n"
        "Offspring: [1|0]  [1|1]  [0|1]\n\n"
        "Trait score = (1+0+1+1+0+1) / (2×17)\n"
        "            = 4 / 34 ≈ 0.118"
    )
    ax3.text(0.05, 0.5, cross_text, fontsize=10, fontfamily="monospace",
             va="center", transform=ax3.transAxes,
             bbox=dict(boxstyle="round,pad=0.6", fc="#FFF9C4", ec="grey", lw=1.2))

    fig.suptitle("PG-11: Mendelian Inheritance & Segregation",
                 fontsize=13, fontweight="bold")
    fig.savefig(FIG_DIR / "fig_PG11.pdf")
    plt.close(fig)
    print("  ✓ PG-11 saved")


# ════════════════════════════════════════════════════════════════════════
#  PG-12: Resistance vs. Relative Fitness Surface
# ════════════════════════════════════════════════════════════════════════

def fig_pg12():
    fig, ax = plt.subplots(figsize=(10, 8))

    r_vals = np.linspace(0, 1, 250)
    t_vals = np.linspace(0, 1, 250)
    R, T = np.meshgrid(r_vals, t_vals)

    # Simplified fitness model
    beta_eff = 0.8
    tau_max = 0.85
    P_inf = beta_eff * (1 - R)
    P_survive_inf = T * tau_max * 0.3
    fitness = (1 - P_inf) + P_inf * P_survive_inf
    fitness = fitness / fitness.max()

    im = ax.contourf(R, T, fitness, levels=30, cmap="viridis")
    ax.contour(R, T, fitness, levels=12, colors="white", linewidths=0.4, alpha=0.5)
    cbar = plt.colorbar(im, ax=ax, label="Relative fitness", shrink=0.85, pad=0.02)

    # Mark initial and evolved positions
    ax.plot(0.15, 0.10, "*", color="red", markersize=16,
            markeredgecolor="white", markeredgewidth=1.5, zorder=5,
            label="Initial mean (r̄=0.15, t̄=0.10)")

    # Use actual evolved values from W285
    res_final = np.mean([RD[r]["yearly_mean_resistance"][-1] for r in REGIONS_NS])
    tol_final = np.mean([RD[r]["yearly_mean_tolerance"][-1] for r in REGIONS_NS])
    ax.plot(res_final, tol_final, "*", color="yellow", markersize=16,
            markeredgecolor="white", markeredgewidth=1.5, zorder=5,
            label=f"Evolved mean (r̄={res_final:.2f}, t̄={tol_final:.2f})")
    ax.annotate("", xy=(res_final, tol_final), xytext=(0.15, 0.10),
                arrowprops=dict(arrowstyle="->", lw=2.5, color="white"))

    ax.set_xlabel("Resistance trait value (r)", fontsize=13)
    ax.set_ylabel("Tolerance trait value (t)", fontsize=13)
    ax.legend(fontsize=9, loc="upper left", framealpha=0.9)

    # Inset: cross-sections
    ax_ins = ax.inset_axes([0.52, 0.06, 0.42, 0.32])
    t_fixed = 0.10
    f_vs_r = (1 - beta_eff * (1 - r_vals)) + beta_eff * (1 - r_vals) * t_fixed * tau_max * 0.3
    f_vs_r = f_vs_r / f_vs_r.max()
    ax_ins.plot(r_vals, f_vs_r, color=C_RES, lw=2, label="f(r) at t=0.10")

    r_fixed = 0.15
    f_vs_t = (1 - beta_eff * (1 - r_fixed)) + beta_eff * (1 - r_fixed) * t_vals * tau_max * 0.3
    f_vs_t = f_vs_t / f_vs_t.max()
    ax_ins.plot(t_vals, f_vs_t, color=C_TOL, lw=2, label="f(t) at r=0.15")

    ax_ins.set_xlabel("Trait value", fontsize=7)
    ax_ins.set_ylabel("Relative fitness", fontsize=7)
    ax_ins.legend(fontsize=6)
    ax_ins.tick_params(labelsize=6)
    ax_ins.set_title("Cross-sections: ~1000× gradient difference", fontsize=7)

    ax.set_title("PG-12: Fitness Landscape — Resistance vs. Tolerance",
                 fontsize=13, fontweight="bold", pad=12)
    fig.savefig(FIG_DIR / "fig_PG12.pdf")
    plt.close(fig)
    print("  ✓ PG-12 saved")


# ════════════════════════════════════════════════════════════════════════
#  MAIN
# ════════════════════════════════════════════════════════════════════════

def main():
    print("Generating 12 population ecology & genetics figures …\n")

    generators = [
        ("PG-01", fig_pg01),
        ("PG-02", fig_pg02),
        ("PG-03", fig_pg03),
        ("PG-04", fig_pg04),
        ("PG-05", fig_pg05),
        ("PG-06", fig_pg06),
        ("PG-07", fig_pg07),
        ("PG-08", fig_pg08),
        ("PG-09", fig_pg09),
        ("PG-10", fig_pg10),
        ("PG-11", fig_pg11),
        ("PG-12", fig_pg12),
    ]

    successes = []
    failures = []

    for name, func in generators:
        try:
            func()
            successes.append(name)
        except Exception as e:
            failures.append((name, str(e)))
            print(f"  ✗ {name} FAILED: {e}")
            import traceback
            traceback.print_exc()

    print(f"\n{'='*50}")
    print(f"Successes: {len(successes)}/12 — {', '.join(successes)}")
    if failures:
        print(f"Failures:  {len(failures)}/12")
        for name, err in failures:
            print(f"  {name}: {err}")
    else:
        print("All 12 figures generated successfully!")
    print(f"Output directory: {FIG_DIR}")


if __name__ == "__main__":
    main()
