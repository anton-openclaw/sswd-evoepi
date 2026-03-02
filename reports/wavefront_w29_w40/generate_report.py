#!/usr/bin/env python3
"""Generate wavefront calibration report for W29-W40."""

import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
from pathlib import Path

BASE = Path("/home/starbot/.openclaw/workspace/sswd-evoepi")
RESULTS = BASE / "results" / "calibration"
OUTDIR = BASE / "reports" / "wavefront_w29_w40"
FIGDIR = OUTDIR / "figures"
FIGDIR.mkdir(parents=True, exist_ok=True)

CONFIGS = {
    "W29": {"act": 50, "K_half": 200000, "T_vbnc": 12},
    "W30": {"act": 50, "K_half": 200000, "T_vbnc": 14},
    "W31": {"act": 50, "K_half": 400000, "T_vbnc": 12},
    "W32": {"act": 50, "K_half": 400000, "T_vbnc": 14},
    "W33": {"act": 50, "K_half": 600000, "T_vbnc": 12},
    "W34": {"act": 50, "K_half": 600000, "T_vbnc": 14},
    "W35": {"act": 100, "K_half": 200000, "T_vbnc": 12},
    "W36": {"act": 100, "K_half": 200000, "T_vbnc": 14},
    "W37": {"act": 100, "K_half": 400000, "T_vbnc": 12},
    "W38": {"act": 100, "K_half": 400000, "T_vbnc": 14},
    "W39": {"act": 100, "K_half": 600000, "T_vbnc": 12},
    "W40": {"act": 100, "K_half": 600000, "T_vbnc": 14},
}

TARGET_REGIONS = ["AK-PWS","AK-FN","AK-FS","BC-N","SS-S","JDF","OR","CA-N"]
TARGETS = {"AK-PWS":50,"AK-FN":50,"AK-FS":20,"BC-N":20,"SS-S":5,"JDF":2,"OR":0.25,"CA-N":0.1}

ALL_REGIONS = [
    "AK-WG","AK-AL","AK-EG","AK-PWS","AK-FN","AK-FS","AK-OC",
    "BC-N","BC-C","SS-N","SS-S","JDF","WA-O",
    "OR","CA-N","CA-C","CA-S","BJ"
]

cmap = plt.cm.coolwarm_r
REGION_COLORS = {r: cmap(i / (len(ALL_REGIONS)-1)) for i, r in enumerate(ALL_REGIONS)}

results = {}
for rnd in CONFIGS:
    rfile = RESULTS / rnd / "result_seed42.json"
    if rfile.exists():
        with open(rfile) as f:
            results[rnd] = json.load(f)
print(f"Loaded {len(results)} results")


def get_recovery_trajectory(data, region):
    rd = data.get("region_details", {}).get(region, {})
    yt = rd.get("yearly_totals", [])
    if not yt or yt[0] == 0:
        return None
    return [y / yt[0] for y in yt]


# ===== FIG 1: Recovery trajectories — T_vbnc=12 vs T_vbnc=14 comparison =====
print("Fig 1: T_vbnc comparison...")
fig, axes = plt.subplots(3, 2, figsize=(16, 18), sharex=True, sharey=True)

for row, khalf_label in enumerate(["200K", "400K", "600K"]):
    khalf = [200000, 400000, 600000][row]
    for col, tvbnc in enumerate([12, 14]):
        ax = axes[row, col]
        # Find matching round (act=50 for cleaner comparison)
        rnd = [r for r, c in CONFIGS.items() if c["K_half"]==khalf and c["T_vbnc"]==tvbnc and c["act"]==50][0]
        data = results.get(rnd)
        if not data:
            continue
        
        cal_years = [2012 + y for y in range(13)]
        for region in ALL_REGIONS:
            traj = get_recovery_trajectory(data, region)
            if traj is None:
                continue
            is_target = region in TARGETS
            lw = 2.5 if is_target else 0.7
            alpha = 1.0 if is_target else 0.25
            ls = '-' if is_target else '--'
            ax.plot(cal_years[:len(traj)], [t*100 for t in traj],
                    color=REGION_COLORS[region], lw=lw, alpha=alpha, ls=ls)
            if is_target and traj:
                final = traj[-1]*100
                if final > 1:
                    ax.annotate(f"{region} {final:.0f}%",
                               xy=(cal_years[len(traj)-1], final),
                               fontsize=6, color=REGION_COLORS[region], fontweight='bold',
                               xytext=(3, 0), textcoords='offset points')
        
        ax.set_title(f"{rnd}: K½={khalf_label}, T_vbnc={tvbnc}", fontsize=11)
        ax.set_ylim(-2, 105)
        ax.grid(True, alpha=0.2)
        if row == 2: ax.set_xlabel("Year")
        if col == 0: ax.set_ylabel("Recovery (%)")

handles = [Line2D([0],[0], color=REGION_COLORS[r], lw=2 if r in TARGETS else 0.7,
                  ls='-' if r in TARGETS else '--',
                  label=f"{r} (T:{TARGETS[r]}%)" if r in TARGETS else r)
           for r in ALL_REGIONS]
fig.legend(handles=handles, loc='center right', bbox_to_anchor=(1.11, 0.5), fontsize=7)
fig.suptitle("T_vbnc=12 (left) vs T_vbnc=14 (right) — Effect of VBNC Midpoint Shift\n"
             "act_thresh=50, D_P_max_range=500km", fontsize=13, fontweight='bold')
plt.tight_layout(rect=[0, 0, 0.89, 0.96])
fig.savefig(FIGDIR / "fig1_tvbnc_comparison.png", dpi=150, bbox_inches='tight')
plt.close()

# ===== FIG 2: K_half sensitivity — 200K vs 400K vs 600K =====
print("Fig 2: K_half sensitivity...")
fig, axes = plt.subplots(2, 3, figsize=(18, 10), sharex=True, sharey=True)

for row, tvbnc in enumerate([12, 14]):
    for col, (khalf, klab) in enumerate([(200000,"200K"),(400000,"400K"),(600000,"600K")]):
        ax = axes[row, col]
        rnd = [r for r, c in CONFIGS.items() if c["K_half"]==khalf and c["T_vbnc"]==tvbnc and c["act"]==50][0]
        data = results.get(rnd)
        if not data:
            continue
        
        cal_years = [2012 + y for y in range(13)]
        for region in ALL_REGIONS:
            traj = get_recovery_trajectory(data, region)
            if traj is None:
                continue
            is_target = region in TARGETS
            ax.plot(cal_years[:len(traj)], [t*100 for t in traj],
                    color=REGION_COLORS[region],
                    lw=2.5 if is_target else 0.7,
                    alpha=1.0 if is_target else 0.25,
                    ls='-' if is_target else '--')
        
        ax.set_title(f"{rnd}: K½={klab}, T_vbnc={tvbnc}", fontsize=10)
        ax.set_ylim(-2, 105)
        ax.grid(True, alpha=0.2)
        if row == 1: ax.set_xlabel("Year")
        if col == 0: ax.set_ylabel("Recovery (%)")

fig.suptitle("K_half Sensitivity: 200K → 400K → 600K\nact_thresh=50, D_P_max_range=500km",
             fontsize=13, fontweight='bold')
plt.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(FIGDIR / "fig2_khalf_sensitivity.png", dpi=150, bbox_inches='tight')
plt.close()

# ===== FIG 3: Heatmap — all 12 rounds =====
print("Fig 3: Summary heatmap...")
rounds = list(CONFIGS.keys())
regs = TARGET_REGIONS

matrix = np.zeros((len(rounds), len(regs)))
for i, rnd in enumerate(rounds):
    data = results.get(rnd, {})
    sc = data.get("scoring", {}).get("per_region", {})
    for j, reg in enumerate(regs):
        matrix[i, j] = sc.get(reg, {}).get("actual_pct", 0)

fig, ax = plt.subplots(figsize=(12, 8))
im = ax.imshow(matrix, aspect='auto', cmap='RdYlGn', vmin=0, vmax=100)

ax.set_xticks(range(len(regs)))
ax.set_xticklabels(regs, rotation=45, ha='right')
ax.set_yticks(range(len(rounds)))
ylabels = [f"{r} (act={CONFIGS[r]['act']}, K½={CONFIGS[r]['K_half']//1000}K, Tv={CONFIGS[r]['T_vbnc']})"
           for r in rounds]
ax.set_yticklabels(ylabels, fontsize=8)

for i in range(len(rounds)):
    for j in range(len(regs)):
        val = matrix[i, j]
        tgt = TARGETS[regs[j]]
        color = 'white' if val > 50 or val < 5 else 'black'
        ax.text(j, i, f"{val:.1f}\n({tgt}%)", ha='center', va='center', fontsize=6.5, color=color)

for j, reg in enumerate(regs):
    ax.text(j, -0.7, f"T:{TARGETS[reg]}%", ha='center', fontsize=8, fontweight='bold', color='blue')

plt.colorbar(im, ax=ax, label="Recovery %", shrink=0.8)
ax.set_title("W29-W40: Final Recovery % by Region\nact_thresh [50,100] × K_half [200K,400K,600K] × T_vbnc [12,14]",
             fontsize=12, fontweight='bold')
plt.tight_layout()
fig.savefig(FIGDIR / "fig3_heatmap.png", dpi=150, bbox_inches='tight')
plt.close()

# ===== FIG 4: Target matching — scatter of actual vs target =====
print("Fig 4: Target matching...")
fig, ax = plt.subplots(figsize=(10, 8))

for rnd in rounds:
    data = results.get(rnd, {})
    sc = data.get("scoring", {}).get("per_region", {})
    cfg = CONFIGS[rnd]
    marker = 'o' if cfg['T_vbnc'] == 12 else 's'
    size = {200000: 40, 400000: 80, 600000: 120}[cfg['K_half']]
    color = 'C0' if cfg['act'] == 50 else 'C1'
    
    for reg in TARGET_REGIONS:
        actual = sc.get(reg, {}).get("actual_pct", 0)
        target = TARGETS[reg]
        ax.scatter(target, actual, marker=marker, s=size, color=color, alpha=0.5,
                  edgecolors='black', linewidths=0.3)

ax.plot([0.01, 100], [0.01, 100], 'k--', alpha=0.3, label='Perfect match')
ax.plot([0.01, 100], [0.02, 200], 'k:', alpha=0.15, label='2× band')
ax.plot([0.01, 100], [0.005, 50], 'k:', alpha=0.15)
ax.set_xscale('symlog', linthresh=0.1)
ax.set_yscale('symlog', linthresh=0.1)
ax.set_xlabel("Target Recovery %", fontsize=11)
ax.set_ylabel("Actual Recovery %", fontsize=11)
ax.set_title("All W29-W40: Actual vs Target Recovery\n(○ T_vbnc=12, □ T_vbnc=14; size=K_half; blue=act50, orange=act100)",
             fontsize=11)
ax.legend()
ax.grid(True, alpha=0.2)
plt.tight_layout()
fig.savefig(FIGDIR / "fig4_target_matching.png", dpi=150, bbox_inches='tight')
plt.close()

# ===== RESULTS.md =====
print("Writing RESULTS.md...")
md = ["# Wavefront Calibration W29-W40 Results", "",
      "**Date:** 2026-03-02",
      "**Design:** activation_threshold [50, 100] × K_half [200K, 400K, 600K] × T_vbnc [12, 14]",
      "**Fixed:** k_vbnc=2.0, s0=0.002, P_env_max=2000, D_P=50, D_P_max_range=500km, 5 Channel Islands origins",
      "",
      "## Critical Finding: Wavefront Cannot Reach Alaska",
      "",
      "Across all 12 runs (and all 12 from W17-W28), disease **never reaches Alaska**.",
      "The wavefront propagates from Channel Islands to BC-N (~24-34 months) but stalls there.",
      "AK-PWS remains at 92-93% in every scenario — completely untouched by disease.",
      "",
      "### Why the Wavefront Stalls",
      "The BC coast north of BC-N has a long stretch (~1000-1500km) of remote coastline",
      "with widely spaced sites. Even at D_P_max_range=500km, pathogen concentrations",
      "dilute below any activation threshold before reaching Alaska. The dispersal",
      "kernel is exponential decay with distance — beyond ~200km, concentrations",
      "approach zero regardless of activation_threshold.",
      "",
      "## What DID Work: Mid-Coast Recovery Gradient",
      "",
      "### T_vbnc=14 is a Powerful Knob",
      "Shifting the VBNC midpoint from 12°C to 14°C dramatically increases recovery:",
      "- BC-N: 3.4% → 11.5% (K½=200K), 21.1% → 40.6% (K½=400K), 42.3% → 62.2% (K½=600K)",
      "- The VBNC sigmoid now gives meaningful winter disease suppression at mid-latitudes",
      "",
      "### K_half Controls Recovery Magnitude",
      "- K½=200K: BC-N 3-14%, SS-S 5-28%",
      "- K½=400K: BC-N 21-48%, SS-S 22-66%", 
      "- K½=600K: BC-N 42-65%, SS-S 39-73%",
      "",
      "### Best Target Matches (Ignoring Alaska)",
      "- **W31** (act=50, K½=400K, T_vbnc=12): BC-N=21.1% ✅, SS-S=22.1% (high), JDF=11.8% (high)",
      "- **W37** (act=100, K½=400K, T_vbnc=12): BC-N=29.3% (close), SS-S=22.7% (high)",
      "- **W30** (act=50, K½=200K, T_vbnc=14): BC-N=11.5% (low), SS-S=27.4% (high), JDF=21.3% (high)",
      "",
      "## Recommendation: Abandon Pure Wavefront for Alaska",
      "",
      "The waterborne dispersal mechanism cannot bridge the BC→AK gap with any",
      "reasonable parameters. Options:",
      "1. **Region-timed seeding**: Seed disease at observed arrival times per region",
      "   (Gravem 2021: CA-S Jun 2013, OR Jan 2014, WA/BC Aug 2015, AK Jan 2017)",
      "2. **Hybrid**: Keep wavefront for CA→BC, add timed seeding for AK",  
      "3. **Extended dispersal kernel**: Heavy-tailed (Lévy) instead of exponential",
      "",
      "Option 1 is most pragmatic and biologically defensible (multiple introduction",
      "pathways likely contributed to the actual northward spread).",
]

# Add summary table
md.extend(["", "## Summary Table", "",
           "| Round | act | K_half | T_vbnc | AK-PWS | BC-N | SS-S | JDF | OR | CA-N |",
           "|-------|-----|--------|--------|--------|------|------|-----|-----|------|"])
for rnd in rounds:
    cfg = CONFIGS[rnd]
    data = results.get(rnd, {})
    sc = data.get("scoring", {}).get("per_region", {})
    vals = [f"{sc.get(r, {}).get('actual_pct', 0):.1f}" for r in TARGET_REGIONS]
    md.append(f"| {rnd} | {cfg['act']} | {cfg['K_half']//1000}K | {cfg['T_vbnc']} | {' | '.join(vals)} |")
md.append(f"\n**Targets:** AK-PWS=50%, AK-FN=50%, AK-FS=20%, BC-N=20%, SS-S=5%, JDF=2%, OR=0.25%, CA-N=0.1%")

with open(OUTDIR / "RESULTS.md", 'w') as f:
    f.write('\n'.join(md))

print("Done!")
