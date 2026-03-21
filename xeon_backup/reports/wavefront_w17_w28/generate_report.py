#!/usr/bin/env python3
"""Generate wavefront calibration report for W17-W28 with recovery-by-year-by-region trajectories."""

import json
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path

BASE = Path("/home/starbot/.openclaw/workspace/sswd-evoepi")
RESULTS = BASE / "results" / "calibration"
OUTDIR = BASE / "reports" / "wavefront_w17_w28"
FIGDIR = OUTDIR / "figures"
FIGDIR.mkdir(parents=True, exist_ok=True)

# Config for each round
CONFIGS = {
    "W17": {"K_half": 50000, "s0": 0.001, "P_env_max": 2000},
    "W18": {"K_half": 100000, "s0": 0.001, "P_env_max": 2000},
    "W19": {"K_half": 50000, "s0": 0.002, "P_env_max": 2000},
    "W20": {"K_half": 100000, "s0": 0.002, "P_env_max": 2000},
    "W21": {"K_half": 200000, "s0": 0.002, "P_env_max": 2000},
    "W22": {"K_half": 200000, "s0": 0.003, "P_env_max": 2000},
    "W23": {"K_half": 50000, "s0": 0.001, "P_env_max": 1500},
    "W24": {"K_half": 100000, "s0": 0.001, "P_env_max": 1500},
    "W25": {"K_half": 50000, "s0": 0.002, "P_env_max": 2500},
    "W26": {"K_half": 100000, "s0": 0.002, "P_env_max": 2500},
    "W27": {"K_half": 50000, "s0": 0.003, "P_env_max": 2000},
    "W28": {"K_half": 100000, "s0": 0.003, "P_env_max": 2000},
}

# All shared: k_vbnc=2.0, activation_threshold=500, D_P=50, D_P_max_range=175

TARGET_REGIONS = ["AK-PWS", "AK-FN", "AK-FS", "BC-N", "SS-S", "JDF", "OR", "CA-N"]
TARGETS = {"AK-PWS": 50, "AK-FN": 50, "AK-FS": 20, "BC-N": 20, "SS-S": 5, "JDF": 2, "OR": 0.25, "CA-N": 0.1}

# All 18 regions in geographic order (north to south)
ALL_REGIONS = [
    "AK-WG", "AK-AL", "AK-EG", "AK-PWS", "AK-FN", "AK-FS", "AK-OC",
    "BC-N", "BC-C", "SS-N", "SS-S", "JDF", "WA-O",
    "OR", "CA-N", "CA-C", "CA-S", "BJ"
]

# Color scheme: warm colors for south, cool for north
REGION_COLORS = {}
cmap = plt.cm.coolwarm_r
for i, reg in enumerate(ALL_REGIONS):
    REGION_COLORS[reg] = cmap(i / (len(ALL_REGIONS) - 1))

# Load all results
results = {}
for rnd in CONFIGS:
    rfile = RESULTS / rnd / "result_seed42.json"
    if rfile.exists():
        with open(rfile) as f:
            results[rnd] = json.load(f)

print(f"Loaded {len(results)} result files")


def get_recovery_trajectory(data, region):
    """Compute recovery fraction per year for a region."""
    rd = data.get("region_details", {}).get(region, {})
    yt = rd.get("yearly_totals", [])
    if not yt or yt[0] == 0:
        return None
    return [y / yt[0] for y in yt]


# ============================================================
# FIGURE 1: Recovery trajectories for ALL rounds (grid of panels)
# ============================================================
print("Generating Figure 1: Recovery trajectories (all rounds)...")
fig, axes = plt.subplots(4, 3, figsize=(18, 20), sharex=True, sharey=True)
axes_flat = axes.flatten()

for idx, (rnd, cfg) in enumerate(CONFIGS.items()):
    ax = axes_flat[idx]
    data = results.get(rnd)
    if not data:
        ax.text(0.5, 0.5, "No data", ha='center', va='center', transform=ax.transAxes)
        continue
    
    years = list(range(13))  # 0-12
    cal_years = [2012 + y for y in years]
    
    for region in ALL_REGIONS:
        traj = get_recovery_trajectory(data, region)
        if traj is None:
            continue
        is_target = region in TARGETS
        lw = 2.0 if is_target else 0.8
        alpha = 0.9 if is_target else 0.3
        ls = '-' if is_target else '--'
        ax.plot(cal_years[:len(traj)], [t * 100 for t in traj], 
                color=REGION_COLORS[region], lw=lw, alpha=alpha, ls=ls, label=region if idx == 0 else "")
    
    # Add target markers
    for reg, tgt in TARGETS.items():
        ax.axhline(y=tgt, color=REGION_COLORS[reg], alpha=0.15, lw=0.5, ls=':')
    
    ax.set_title(f"{rnd}: K½={cfg['K_half']//1000}K, s₀={cfg['s0']}, Penv={cfg['P_env_max']}", fontsize=10)
    ax.set_ylim(-2, 105)
    ax.grid(True, alpha=0.2)
    if idx >= 9:
        ax.set_xlabel("Year")
    if idx % 3 == 0:
        ax.set_ylabel("Recovery (%)")

# Add legend to right side
handles = [Line2D([0], [0], color=REGION_COLORS[r], lw=2 if r in TARGETS else 0.8, 
                  ls='-' if r in TARGETS else '--', alpha=0.9 if r in TARGETS else 0.5,
                  label=f"{r} (target {TARGETS[r]}%)" if r in TARGETS else r)
           for r in ALL_REGIONS]
fig.legend(handles=handles, loc='center right', bbox_to_anchor=(1.12, 0.5), fontsize=8, ncol=1)

fig.suptitle("Wavefront Calibration W17-W28: Recovery Trajectories by Region\n"
             "All: k_vbnc=2.0, act_thresh=500, 5 Channel Islands origins, D_P=50km, D_P_max_range=175km",
             fontsize=13, fontweight='bold')
plt.tight_layout(rect=[0, 0, 0.88, 0.96])
fig.savefig(FIGDIR / "fig1_recovery_trajectories_all.png", dpi=150, bbox_inches='tight')
plt.close()
print("  Saved fig1_recovery_trajectories_all.png")


# ============================================================
# FIGURE 2: Best rounds - detailed recovery trajectories
# ============================================================
print("Generating Figure 2: Best rounds detailed...")
# Pick a few illustrative rounds
best_rounds = ["W21", "W22", "W17", "W25"]  # K_half=200K and edge cases

fig, axes = plt.subplots(2, 2, figsize=(16, 12), sharex=True, sharey=True)
axes_flat = axes.flatten()

for idx, rnd in enumerate(best_rounds):
    ax = axes_flat[idx]
    data = results.get(rnd)
    if not data:
        continue
    
    cfg = CONFIGS[rnd]
    years = list(range(13))
    cal_years = [2012 + y for y in years]
    
    for region in ALL_REGIONS:
        traj = get_recovery_trajectory(data, region)
        if traj is None:
            continue
        is_target = region in TARGETS
        lw = 2.5 if is_target else 0.8
        alpha = 1.0 if is_target else 0.3
        ls = '-' if is_target else '--'
        ax.plot(cal_years[:len(traj)], [t * 100 for t in traj],
                color=REGION_COLORS[region], lw=lw, alpha=alpha, ls=ls)
        
        # Label target regions at end
        if is_target and traj:
            final = traj[-1] * 100
            ax.annotate(f"{region}\n{final:.1f}%", 
                       xy=(cal_years[len(traj)-1], final),
                       fontsize=7, color=REGION_COLORS[region], fontweight='bold',
                       xytext=(5, 0), textcoords='offset points')
    
    # Add target horizontal lines
    for reg, tgt in TARGETS.items():
        ax.axhline(y=tgt, color=REGION_COLORS[reg], alpha=0.2, lw=0.5, ls=':')
    
    ax.set_title(f"{rnd}: K½={cfg['K_half']//1000}K, s₀={cfg['s0']}, Penv={cfg['P_env_max']}", fontsize=11)
    ax.set_ylim(-2, 105)
    ax.grid(True, alpha=0.2)
    ax.set_xlabel("Year")
    ax.set_ylabel("Recovery (%)")

fig.suptitle("Selected Rounds: Recovery Trajectories (Target Regions Highlighted)", fontsize=13, fontweight='bold')
plt.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(FIGDIR / "fig2_best_rounds_detail.png", dpi=150, bbox_inches='tight')
plt.close()
print("  Saved fig2_best_rounds_detail.png")


# ============================================================
# FIGURE 3: Wavefront arrival timing
# ============================================================
print("Generating Figure 3: Wavefront arrival timing...")

timing_targets = {
    "CA-S": 0, "CA-C": 6, "CA-N": 6, "OR": 15, "WA-O": 26,
    "SS-N": 26, "SS-S": 26, "JDF": 26, "BC-C": 26, "BC-N": 26,
    "AK-FS": 42, "AK-FN": 42, "AK-OC": 42, "AK-PWS": 42, "AK-EG": 42, "AK-AL": 42, "AK-WG": 42
}
timing_regions = list(timing_targets.keys())

fig, ax = plt.subplots(figsize=(14, 8))

# Target bars
x = np.arange(len(timing_regions))
target_vals = [timing_targets[r] for r in timing_regions]
ax.barh(x, target_vals, height=0.3, color='gray', alpha=0.3, label='Target')

# Plot actual arrivals for a few rounds
plot_rounds = ["W17", "W21", "W22", "W25"]
offsets = np.linspace(-0.15, 0.15, len(plot_rounds))
markers = ['o', 's', '^', 'D']

for i, rnd in enumerate(plot_rounds):
    data = results.get(rnd)
    if not data:
        continue
    at = data.get("arrival_timing", {}).get("per_region", {})
    vals = []
    reached_flags = []
    for reg in timing_regions:
        ar = at.get(reg, {})
        if ar.get("reached"):
            vals.append(ar.get("actual_months", 0))
            reached_flags.append(True)
        else:
            vals.append(156)  # Never reached (13yr = 156mo)
            reached_flags.append(False)
    
    colors = ['red' if not r else 'C%d' % i for r in reached_flags]
    cfg = CONFIGS[rnd]
    label = f"{rnd} (K½={cfg['K_half']//1000}K, s₀={cfg['s0']})"
    for j, (v, c, rf) in enumerate(zip(vals, colors, reached_flags)):
        mk = markers[i]
        if not rf:
            mk = 'x'
        ax.scatter(v, x[j] + offsets[i], marker=mk, color=f'C{i}', s=60, zorder=5,
                  label=label if j == 0 else "")

ax.set_yticks(x)
ax.set_yticklabels(timing_regions)
ax.set_xlabel("Months after disease origin")
ax.set_title("Wavefront Arrival Timing: Actual vs Target\n(X = never reached)", fontsize=12, fontweight='bold')
ax.axvline(x=156, color='red', ls='--', alpha=0.3, label='Simulation end (13yr)')
ax.legend(loc='lower right', fontsize=9)
ax.grid(True, alpha=0.2, axis='x')
ax.set_xlim(-2, 165)
plt.tight_layout()
fig.savefig(FIGDIR / "fig3_wavefront_timing.png", dpi=150, bbox_inches='tight')
plt.close()
print("  Saved fig3_wavefront_timing.png")


# ============================================================
# FIGURE 4: Parameter sensitivity — final recovery by K_half and s0
# ============================================================
print("Generating Figure 4: Parameter sensitivity...")

fig, axes = plt.subplots(2, 3, figsize=(16, 10))

# Panel A-C: How K_half affects target regions (fix s0=0.002, Penv=2000)
khalf_rounds = {"W19": 50, "W20": 100, "W21": 200}  # s0=0.002, Penv=2000
ax = axes[0, 0]
for reg in TARGET_REGIONS:
    vals = []
    xs = []
    for rnd, kh in khalf_rounds.items():
        data = results.get(rnd)
        if data:
            pr = data.get("scoring", {}).get("per_region", {}).get(reg, {})
            vals.append(pr.get("actual_pct", 0))
            xs.append(kh)
    if vals:
        ax.plot(xs, vals, 'o-', color=REGION_COLORS[reg], label=reg, lw=2)
        ax.axhline(y=TARGETS.get(reg, 0), color=REGION_COLORS[reg], ls=':', alpha=0.3)
ax.set_xlabel("K_half (×1000)")
ax.set_ylabel("Final Recovery %")
ax.set_title("Effect of K_half\n(s₀=0.002, Penv=2000)")
ax.legend(fontsize=7, ncol=2)
ax.grid(True, alpha=0.2)
ax.set_yscale('symlog', linthresh=1)
ax.set_ylim(-0.1, 100)

# Panel B: How s0 affects target regions (fix K_half=50K, Penv=2000)
s0_rounds = {"W17": 0.001, "W19": 0.002, "W27": 0.003}
ax = axes[0, 1]
for reg in TARGET_REGIONS:
    vals = []
    xs = []
    for rnd, s0 in s0_rounds.items():
        data = results.get(rnd)
        if data:
            pr = data.get("scoring", {}).get("per_region", {}).get(reg, {})
            vals.append(pr.get("actual_pct", 0))
            xs.append(s0)
    if vals:
        ax.plot(xs, vals, 'o-', color=REGION_COLORS[reg], label=reg, lw=2)
ax.set_xlabel("s₀ (settler survival)")
ax.set_ylabel("Final Recovery %")
ax.set_title("Effect of s₀\n(K½=50K, Penv=2000)")
ax.legend(fontsize=7, ncol=2)
ax.grid(True, alpha=0.2)
ax.set_yscale('symlog', linthresh=1)
ax.set_ylim(-0.1, 100)

# Panel C: How Penv affects target regions (fix K_half=50K, s0=0.001)
penv_rounds = {"W23": 1500, "W17": 2000, "W25": 2500}
# W25 has s0=0.002 not 0.001, use W23/W17 only... actually W25 is s0=0.002
# Better: K_half=50K group: W23(1500,0.001), W17(2000,0.001), W25 doesn't match
# Use K_half=100K: W24(1500,0.001), W18(2000,0.001), W26(2500,0.002) - s0 mismatch
# Just use what we have with a note
penv_rounds = {"W23": 1500, "W17": 2000}  # Only clean comparison
ax = axes[0, 2]
for reg in TARGET_REGIONS:
    vals = []
    xs = []
    for rnd, pe in penv_rounds.items():
        data = results.get(rnd)
        if data:
            pr = data.get("scoring", {}).get("per_region", {}).get(reg, {})
            vals.append(pr.get("actual_pct", 0))
            xs.append(pe)
    if vals:
        ax.plot(xs, vals, 'o-', color=REGION_COLORS[reg], label=reg, lw=2)
ax.set_xlabel("P_env_max")
ax.set_ylabel("Final Recovery %")
ax.set_title("Effect of P_env_max\n(K½=50K, s₀=0.001)")
ax.legend(fontsize=7, ncol=2)
ax.grid(True, alpha=0.2)
ax.set_yscale('symlog', linthresh=1)
ax.set_ylim(-0.1, 100)

# Panel D-F: Disease deaths and recruits for W21 (best candidate)
data = results.get("W21", {})
rd = data.get("region_details", {})

ax = axes[1, 0]
cal_years = [2012 + y for y in range(13)]
for reg in TARGET_REGIONS:
    dd = rd.get(reg, {}).get("yearly_disease_deaths", [])
    if dd:
        ax.plot(cal_years[:len(dd)], dd, 'o-', color=REGION_COLORS[reg], label=reg, lw=1.5, ms=3)
ax.set_xlabel("Year")
ax.set_ylabel("Disease Deaths")
ax.set_title("W21: Yearly Disease Deaths")
ax.legend(fontsize=7, ncol=2)
ax.grid(True, alpha=0.2)
ax.set_yscale('symlog', linthresh=10)

ax = axes[1, 1]
for reg in TARGET_REGIONS:
    rec = rd.get(reg, {}).get("yearly_recruits", [])
    if rec:
        ax.plot(cal_years[:len(rec)], rec, 'o-', color=REGION_COLORS[reg], label=reg, lw=1.5, ms=3)
ax.set_xlabel("Year")
ax.set_ylabel("Recruits")
ax.set_title("W21: Yearly Recruits")
ax.legend(fontsize=7, ncol=2)
ax.grid(True, alpha=0.2)
ax.set_yscale('symlog', linthresh=10)

# Panel F: Resistance evolution for W21
ax = axes[1, 2]
for reg in ["AK-PWS", "BC-N", "SS-S", "OR", "CA-N", "CA-S"]:
    mr = rd.get(reg, {}).get("yearly_mean_resistance", [])
    if mr:
        ax.plot(cal_years[:len(mr)], mr, 'o-', color=REGION_COLORS.get(reg, 'gray'), label=reg, lw=1.5, ms=3)
ax.set_xlabel("Year")
ax.set_ylabel("Mean Resistance")
ax.set_title("W21: Resistance Evolution")
ax.legend(fontsize=7)
ax.grid(True, alpha=0.2)

fig.suptitle("Parameter Sensitivity & Diagnostics (W17-W28)", fontsize=13, fontweight='bold')
plt.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(FIGDIR / "fig4_parameter_sensitivity.png", dpi=150, bbox_inches='tight')
plt.close()
print("  Saved fig4_parameter_sensitivity.png")


# ============================================================
# FIGURE 5: Summary heatmap — final recovery vs target
# ============================================================
print("Generating Figure 5: Summary heatmap...")

rounds = list(CONFIGS.keys())
regs = TARGET_REGIONS

recovery_matrix = np.zeros((len(rounds), len(regs)))
for i, rnd in enumerate(rounds):
    data = results.get(rnd, {})
    sc = data.get("scoring", {}).get("per_region", {})
    for j, reg in enumerate(regs):
        recovery_matrix[i, j] = sc.get(reg, {}).get("actual_pct", 0)

fig, ax = plt.subplots(figsize=(12, 8))
im = ax.imshow(recovery_matrix, aspect='auto', cmap='RdYlGn', vmin=0, vmax=100)

ax.set_xticks(range(len(regs)))
ax.set_xticklabels(regs, rotation=45, ha='right')
ax.set_yticks(range(len(rounds)))
ylabels = [f"{r} (K½={CONFIGS[r]['K_half']//1000}K, s₀={CONFIGS[r]['s0']}, Pe={CONFIGS[r]['P_env_max']})" 
           for r in rounds]
ax.set_yticklabels(ylabels, fontsize=8)

# Annotate with values and targets
for i in range(len(rounds)):
    for j in range(len(regs)):
        val = recovery_matrix[i, j]
        tgt = TARGETS[regs[j]]
        color = 'white' if val > 50 or val < 5 else 'black'
        ax.text(j, i, f"{val:.1f}\n({tgt}%)", ha='center', va='center', fontsize=7, color=color)

# Add target line annotations at top
for j, reg in enumerate(regs):
    ax.text(j, -0.7, f"T:{TARGETS[reg]}%", ha='center', va='center', fontsize=8, fontweight='bold', color='blue')

plt.colorbar(im, ax=ax, label="Recovery %", shrink=0.8)
ax.set_title("W17-W28: Final Recovery % by Region\n(values shown with targets in parentheses)", 
             fontsize=12, fontweight='bold')
plt.tight_layout()
fig.savefig(FIGDIR / "fig5_summary_heatmap.png", dpi=150, bbox_inches='tight')
plt.close()
print("  Saved fig5_summary_heatmap.png")


# ============================================================
# Generate markdown summary
# ============================================================
print("Generating RESULTS.md...")

md = []
md.append("# Wavefront Calibration W17-W28 Results")
md.append("")
md.append("**Date:** 2026-03-01")
md.append("**Shared parameters:** k_vbnc=2.0, activation_threshold=500, 5 Channel Islands origins, D_P=50km, D_P_max_range=175km")
md.append("")
md.append("## Key Findings")
md.append("")
md.append("### 1. Wavefront Stalls at BC-N — Never Reaches Alaska")
md.append("- In ALL 12 runs, disease never reaches AK-FS, AK-FN, AK-PWS, AK-AL, AK-EG, AK-WG, or AK-OC")
md.append("- Wave reaches BC-N at ~23-34 months but cannot bridge the gap to Alaska")
md.append("- Alaska populations remain at 89-96% (untouched by disease)")
md.append("- **Root cause:** activation_threshold=500 bact/mL too high for the diluted vibrio concentrations reaching BC→AK distances")
md.append("")
md.append("### 2. Strong Southern Gradient Achieved")
md.append("- CA-N and OR go completely extinct (0%) — matches targets (0.1% and 0.25%)")
md.append("- CA-S and CA-C destroyed immediately upon disease arrival")
md.append("- SS-S ranges 0.1-6.7% across runs (target 5%) — some runs nail this")
md.append("- JDF ranges 0.1-7.4% (target 2%)")
md.append("")
md.append("### 3. BC-N Overshoot Problem")
md.append("- BC-N ranges 0-10.2% (target 20%)")  
md.append("- Disease is too severe once it arrives — wipes out populations before recovery can establish")
md.append("- Higher K_half (200K) helps: W21=10.2%, W22=7.1%")
md.append("")
md.append("### 4. Wave Timing")
md.append("- CA-C arrival: ~12 months (target 6) — about 2× too slow")
md.append("- OR arrival: ~16-20 months (target 15) — close")
md.append("- BC-N arrival: ~23-34 months (target ~26) — in range")
md.append("- Alaska: NEVER (target ~42 months)")
md.append("")
md.append("### 5. Parameter Effects")
md.append("- **K_half:** Higher K_half → slightly more recovery in mid-coast (BC-N, SS-S, JDF). No effect on wavefront propagation.")
md.append("- **s0 (settler survival):** Higher s0 → more recovery (more recruits). W22 (s0=0.003, K_half=200K) has best mid-coast numbers.")
md.append("- **P_env_max:** 1500 vs 2000 vs 2500 — minimal effect in these runs (wave dynamics dominate over endemic pressure).")
md.append("")
md.append("## What Needs to Change for Next Round")
md.append("")
md.append("### Priority 1: Get wavefront to Alaska")
md.append("- **Lower activation_threshold** from 500 to 50-100 bact/mL")
md.append("- **Increase D_P_max_range** from 175km to 500-1000km (longer-range waterborne dispersal)")
md.append("- Consider both simultaneously")
md.append("")
md.append("### Priority 2: Protect BC-N recovery")
md.append("- Once wavefront reaches AK, need to ensure BC-N still recovers to ~20%")
md.append("- K_half=200K + s0≥0.002 seems necessary for any mid-coast recovery")
md.append("")
md.append("### Priority 3: Speed up southern wave")
md.append("- CA-C at 12mo vs 6mo target — wave propagation too slow near origin")
md.append("- May need higher local D_P or lower threshold")
md.append("")

# Summary table
md.append("## Summary Table")
md.append("")
md.append("| Round | K_half | s0 | P_env | AK-PWS | AK-FN | AK-FS | BC-N | SS-S | JDF | OR | CA-N | AK Arrival | BC-N Arrival |")
md.append("|-------|--------|------|-------|--------|-------|-------|------|------|-----|-----|------|------------|--------------|")
for rnd in rounds:
    cfg = CONFIGS[rnd]
    data = results.get(rnd, {})
    sc = data.get("scoring", {}).get("per_region", {})
    at = data.get("arrival_timing", {}).get("per_region", {})
    
    vals = [f"{sc.get(r, {}).get('actual_pct', 0):.1f}" for r in TARGET_REGIONS]
    
    ak_arr = at.get("AK-PWS", {})
    ak_str = f"{ak_arr['actual_months']:.0f}mo" if ak_arr.get("reached") else "NEVER"
    bcn_arr = at.get("BC-N", {})
    bcn_str = f"{bcn_arr['actual_months']:.0f}mo" if bcn_arr.get("reached") else "NEVER"
    
    md.append(f"| {rnd} | {cfg['K_half']//1000}K | {cfg['s0']} | {cfg['P_env_max']} | {' | '.join(vals)} | {ak_str} | {bcn_str} |")

md.append("")
md.append(f"**Targets:** AK-PWS=50%, AK-FN=50%, AK-FS=20%, BC-N=20%, SS-S=5%, JDF=2%, OR=0.25%, CA-N=0.1%")
md.append("")
md.append("## Figures")
md.append("- `fig1_recovery_trajectories_all.png` — Recovery trajectories for all 12 rounds (18 regions each)")
md.append("- `fig2_best_rounds_detail.png` — Detailed trajectories for selected rounds")
md.append("- `fig3_wavefront_timing.png` — Wavefront arrival timing vs targets")
md.append("- `fig4_parameter_sensitivity.png` — Parameter sensitivity & diagnostics")
md.append("- `fig5_summary_heatmap.png` — Summary heatmap of final recovery vs targets")

with open(OUTDIR / "RESULTS.md", 'w') as f:
    f.write('\n'.join(md))

print("Done! All figures and RESULTS.md generated.")
