#!/usr/bin/env python3
"""K Sensitivity Analysis — Multi-scale assessment of carrying capacity effects.

Analyzes how reducing K (for local diagnostic runs) affects:
1. Aggregate metrics (RMSLE, within-2x/5x)
2. Per-region recovery fractions and their variance
3. Population dynamics (trajectories, crash timing, nadir)
4. Evolutionary response (resistance gain, Va dynamics)
5. Spatial structure (wavefront timing, gradient)
6. Stochastic noise floor (inter-seed CV per region per K)

Key question: What is the minimum K where the model's behavior is
qualitatively faithful to K=5000 (Xeon runs)?
"""

import json
import os
import sys
import numpy as np
from pathlib import Path
from collections import defaultdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# ── Setup ───────────────────────────────────────────────────────────────
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

BASE = Path(__file__).resolve().parent.parent / 'results' / 'k_cv_sweep'
OUTDIR = BASE / 'analysis'
OUTDIR.mkdir(exist_ok=True)

K_VALUES = [500, 1000, 2000]
SEEDS = [42, 137, 256]

# Coastline ordering (S→N)
COASTLINE_ORDER = [
    "BJ", "CA-S", "CA-C", "CA-N", "OR", "WA-O", "JDF", "SS-S", "SS-N",
    "BC-C", "BC-N", "AK-FS", "AK-FN", "AK-PWS", "AK-EG", "AK-OC",
    "AK-WG", "AK-AL",
]

# 8 scored regions + targets
from sswd_evoepi.metrics import RECOVERY_TARGETS
SCORED = list(RECOVERY_TARGETS.keys())

# Colors per K
K_COLORS = {500: '#e74c3c', 1000: '#f39c12', 2000: '#2ecc71', 5000: '#3498db'}

# ── Load data ───────────────────────────────────────────────────────────
def load_all():
    """Load all result JSONs into nested dict: data[K][seed] = result_dict."""
    data = {}
    for K in K_VALUES:
        data[K] = {}
        for seed in SEEDS:
            path = BASE / f'K{K}_seed{seed}' / f'result_seed{seed}.json'
            if path.exists():
                with open(path) as f:
                    data[K][seed] = json.load(f)
    return data


def load_npz(K, seed):
    """Load monthly NPZ for a given K and seed."""
    path = BASE / f'K{K}_seed{seed}' / f'monthly_seed{seed}.npz'
    if path.exists():
        return np.load(path)
    return None


# ── 1. Aggregate metrics ───────────────────────────────────────────────
def analyze_aggregate(data):
    """RMSLE mean, std, CV across seeds for each K."""
    print("\n" + "="*70)
    print("1. AGGREGATE METRICS")
    print("="*70)

    rows = []
    for K in K_VALUES:
        rmsles = []
        w2x = []
        w5x = []
        for seed in SEEDS:
            if seed in data[K]:
                s = data[K][seed]['scoring']
                rmsles.append(s['rmsle'])
                w2x.append(s['within_2x'])
                w5x.append(s['within_5x'])
        if len(rmsles) >= 2:
            mean_r = np.mean(rmsles)
            std_r = np.std(rmsles, ddof=1)
            cv_r = std_r / mean_r * 100
            rows.append((K, mean_r, std_r, cv_r, np.mean(w2x), np.mean(w5x), rmsles))
            print(f"  K={K:5d}: RMSLE = {mean_r:.4f} ± {std_r:.4f} (CV={cv_r:.1f}%)"
                  f"  within-2x={np.mean(w2x):.0f}/8  within-5x={np.mean(w5x):.0f}/8"
                  f"  [{', '.join(f'{r:.3f}' for r in rmsles)}]")
        else:
            print(f"  K={K:5d}: only {len(rmsles)} seeds (need ≥2)")

    return rows


# ── 2. Per-region recovery ─────────────────────────────────────────────
def analyze_regional(data):
    """Per-region recovery: mean, std, CV across seeds for each K."""
    print("\n" + "="*70)
    print("2. PER-REGION RECOVERY (scored regions)")
    print("="*70)

    # Header
    print(f"  {'Region':>8s}  {'Target':>7s}", end="")
    for K in K_VALUES:
        print(f"  │ K={K:4d} mean   std    CV", end="")
    print()
    print("  " + "-"*90)

    region_data = {}  # region -> {K -> [recovery values]}
    for region in SCORED:
        target = RECOVERY_TARGETS[region]
        region_data[region] = {}
        print(f"  {region:>8s}  {target*100:6.1f}%", end="")
        for K in K_VALUES:
            vals = []
            for seed in SEEDS:
                if seed in data[K]:
                    rr = data[K][seed]['region_recovery']
                    vals.append(rr.get(region, 0))
            region_data[region][K] = vals
            if len(vals) >= 2:
                m = np.mean(vals) * 100
                s = np.std(vals, ddof=1) * 100
                cv = s / m * 100 if m > 0 else float('inf')
                print(f"  │ {m:6.2f}% {s:5.2f}% {cv:5.1f}%", end="")
            else:
                print(f"  │   {'N/A':>6s}  {'N/A':>5s}  {'N/A':>5s}", end="")
        print()

    return region_data


# ── 3. Population dynamics ─────────────────────────────────────────────
def analyze_dynamics(data):
    """Compare population trajectories: crash depth, nadir year, recovery shape."""
    print("\n" + "="*70)
    print("3. POPULATION DYNAMICS (coast-wide)")
    print("="*70)

    dyn = {}  # K -> {seed -> yearly coast-wide pop}
    for K in K_VALUES:
        dyn[K] = {}
        for seed in SEEDS:
            if seed not in data[K]:
                continue
            rd = data[K][seed]['region_details']
            n_years = len(rd[COASTLINE_ORDER[0]]['yearly_totals'])
            coast_pop = np.zeros(n_years)
            for r in COASTLINE_ORDER:
                if r in rd:
                    coast_pop += np.array(rd[r]['yearly_totals'])
            # Normalize by K to compare shapes
            total_K = sum(rd[r]['n_nodes'] for r in COASTLINE_ORDER if r in rd) * K
            dyn[K][seed] = coast_pop / total_K  # fraction of total capacity

    print(f"  {'K':>5s}  {'seed':>5s}  {'peak_frac':>9s}  {'nadir_frac':>10s}  {'nadir_yr':>8s}  {'final_frac':>10s}  {'crash%':>6s}")
    for K in K_VALUES:
        for seed in sorted(dyn[K].keys()):
            ts = dyn[K][seed]
            peak = ts[:3].max()  # pre-disease peak (first 3 years)
            nadir = ts.min()
            nadir_yr = np.argmin(ts)
            final = ts[-1]
            crash = (peak - nadir) / peak * 100
            print(f"  {K:5d}  {seed:5d}  {peak:9.4f}  {nadir:10.4f}  {nadir_yr:8d}  {final:10.4f}  {crash:5.1f}%")

    return dyn


# ── 4. Evolutionary response ──────────────────────────────────────────
def analyze_evolution(data):
    """Compare resistance evolution: Δr̄, rate of change, Va trajectories."""
    print("\n" + "="*70)
    print("4. EVOLUTIONARY RESPONSE")
    print("="*70)

    evo = {}  # K -> {seed -> {region -> (r_trajectory, va_trajectory)}}
    for K in K_VALUES:
        evo[K] = {}
        for seed in SEEDS:
            if seed not in data[K]:
                continue
            rd = data[K][seed]['region_details']
            evo[K][seed] = {}
            for r in COASTLINE_ORDER:
                if r in rd:
                    evo[K][seed][r] = {
                        'r': rd[r]['yearly_mean_resistance'],
                        't': rd[r]['yearly_mean_tolerance'],
                        'c': rd[r]['yearly_mean_recovery'],
                        'va_r': rd[r]['yearly_va_resistance'],
                    }

    # Coast-wide mean resistance shift
    print(f"\n  Coast-wide mean resistance: initial → final (Δr̄)")
    print(f"  {'K':>5s}  {'seed':>5s}  {'r̄₀':>7s}  {'r̄_final':>8s}  {'Δr̄':>7s}  {'max Va(r)':>9s}")
    for K in K_VALUES:
        for seed in sorted(evo[K].keys()):
            # Population-weighted mean across regions
            rd = data[K][seed]['region_details']
            r_trajectories = []
            va_trajectories = []
            weights = []
            for r in COASTLINE_ORDER:
                if r in rd:
                    r_trajectories.append(rd[r]['yearly_mean_resistance'])
                    va_trajectories.append(rd[r]['yearly_va_resistance'])
                    weights.append(np.mean(rd[r]['yearly_totals']))
            weights = np.array(weights)
            weights /= weights.sum()
            r_traj = np.average(r_trajectories, axis=0, weights=weights)
            va_traj = np.average(va_trajectories, axis=0, weights=weights)
            r0 = r_traj[0]
            r_final = r_traj[-1]
            delta_r = r_final - r0
            max_va = np.max(va_traj)
            print(f"  {K:5d}  {seed:5d}  {r0:7.4f}  {r_final:8.4f}  {delta_r:7.4f}  {max_va:9.6f}")

    return evo


# ── 5. Spatial gradient ────────────────────────────────────────────────
def analyze_gradient(data):
    """Compare the N-S recovery gradient across K values."""
    print("\n" + "="*70)
    print("5. SPATIAL GRADIENT (recovery by latitude)")
    print("="*70)

    # For each K, compute mean recovery per region (averaged over seeds)
    gradient = {}  # K -> {region -> mean_recovery}
    for K in K_VALUES:
        gradient[K] = {}
        for region in COASTLINE_ORDER:
            vals = []
            for seed in SEEDS:
                if seed in data[K]:
                    vals.append(data[K][seed]['region_recovery'].get(region, 0))
            if vals:
                gradient[K][region] = np.mean(vals)

    # Print gradient table
    print(f"  {'Region':>8s}", end="")
    for K in K_VALUES:
        print(f"  K={K:4d}", end="")
    print(f"  │ max_diff")
    for region in COASTLINE_ORDER:
        print(f"  {region:>8s}", end="")
        vals_across_K = []
        for K in K_VALUES:
            v = gradient[K].get(region, 0) * 100
            vals_across_K.append(v)
            print(f"  {v:5.1f}%", end="")
        diff = max(vals_across_K) - min(vals_across_K)
        print(f"  │ {diff:5.1f}pp")

    return gradient


# ── 6. Stochastic noise profile ───────────────────────────────────────
def analyze_noise(data):
    """Where does stochastic noise hit hardest? CV per region per K."""
    print("\n" + "="*70)
    print("6. STOCHASTIC NOISE PROFILE")
    print("="*70)

    print(f"  Regions with highest inter-seed CV at each K:")
    for K in K_VALUES:
        region_cvs = []
        for region in COASTLINE_ORDER:
            vals = []
            for seed in SEEDS:
                if seed in data[K]:
                    vals.append(data[K][seed]['region_recovery'].get(region, 0))
            if len(vals) >= 2 and np.mean(vals) > 0.001:
                cv = np.std(vals, ddof=1) / np.mean(vals) * 100
                region_cvs.append((region, cv, np.mean(vals)*100))
        region_cvs.sort(key=lambda x: -x[1])
        print(f"\n  K={K}:")
        for r, cv, mean in region_cvs[:5]:
            print(f"    {r:8s}: CV={cv:5.1f}% (mean={mean:5.2f}%)")


# ── 7. Wavefront timing ───────────────────────────────────────────────
def analyze_wavefront(data):
    """Compare disease arrival timing across K values."""
    print("\n" + "="*70)
    print("7. WAVEFRONT ARRIVAL TIMING")
    print("="*70)

    for K in K_VALUES:
        for seed in SEEDS:
            if seed not in data[K]:
                continue
            at = data[K][seed].get('arrival_timing', {})
            if at:
                mae = at.get('mae_months', 'N/A')
                pr = at.get('per_region', {})
                print(f"  K={K}, seed={seed}: MAE={mae}")


# ── Figures ─────────────────────────────────────────────────────────────
def plot_summary(data, dyn, evo, gradient):
    """Generate comprehensive summary figure."""
    fig = plt.figure(figsize=(20, 16))
    gs = GridSpec(3, 3, figure=fig, hspace=0.35, wspace=0.3)

    # Panel A: RMSLE by K (box/scatter)
    ax_a = fig.add_subplot(gs[0, 0])
    for i, K in enumerate(K_VALUES):
        rmsles = [data[K][s]['scoring']['rmsle'] for s in SEEDS if s in data[K]]
        ax_a.scatter([K]*len(rmsles), rmsles, color=K_COLORS[K], s=80, zorder=5, alpha=0.8)
        if len(rmsles) >= 2:
            ax_a.errorbar(K, np.mean(rmsles), yerr=np.std(rmsles, ddof=1),
                         color=K_COLORS[K], capsize=6, lw=2, zorder=4)
    ax_a.set_xlabel('K (carrying capacity per node)')
    ax_a.set_ylabel('RMSLE')
    ax_a.set_title('A. RMSLE vs K', fontweight='bold')
    ax_a.set_xticks(K_VALUES)

    # Panel B: Per-region recovery (scored regions) — grouped bars
    ax_b = fig.add_subplot(gs[0, 1:])
    x = np.arange(len(SCORED))
    width = 0.25
    for i, K in enumerate(K_VALUES):
        means = []
        errs = []
        for region in SCORED:
            vals = [data[K][s]['region_recovery'].get(region, 0) for s in SEEDS if s in data[K]]
            means.append(np.mean(vals) * 100)
            errs.append(np.std(vals, ddof=1) * 100 if len(vals) >= 2 else 0)
        ax_b.bar(x + i*width, means, width, yerr=errs, label=f'K={K}',
                color=K_COLORS[K], alpha=0.8, capsize=3)
    # Add target line markers
    for j, region in enumerate(SCORED):
        target = RECOVERY_TARGETS[region] * 100
        ax_b.plot([j - 0.15, j + len(K_VALUES)*width - 0.15], [target, target],
                 'k--', lw=1, alpha=0.5)
    ax_b.set_xticks(x + width)
    ax_b.set_xticklabels(SCORED, rotation=45, ha='right')
    ax_b.set_ylabel('Recovery (%)')
    ax_b.set_title('B. Regional Recovery by K (dashed = target)', fontweight='bold')
    ax_b.legend(fontsize=9)

    # Panel C: Normalized population trajectories (coast-wide)
    ax_c = fig.add_subplot(gs[1, 0])
    for K in K_VALUES:
        for seed in sorted(dyn[K].keys()):
            ts = dyn[K][seed]
            alpha = 0.3 if len(dyn[K]) > 1 else 0.8
            ax_c.plot(ts, color=K_COLORS[K], alpha=alpha, lw=1)
        # Mean trajectory
        if len(dyn[K]) >= 2:
            mean_ts = np.mean([dyn[K][s] for s in dyn[K]], axis=0)
            ax_c.plot(mean_ts, color=K_COLORS[K], lw=2.5, label=f'K={K}')
    ax_c.set_xlabel('Simulation Year')
    ax_c.set_ylabel('Population / Total Capacity')
    ax_c.set_title('C. Normalized Pop. Trajectories', fontweight='bold')
    ax_c.legend(fontsize=9)

    # Panel D: Resistance evolution trajectories
    ax_d = fig.add_subplot(gs[1, 1])
    for K in K_VALUES:
        for seed in SEEDS:
            if seed not in data[K]:
                continue
            rd = data[K][seed]['region_details']
            # Coast-wide weighted mean resistance
            r_traj = []
            weights = []
            for r in COASTLINE_ORDER:
                if r in rd:
                    r_traj.append(rd[r]['yearly_mean_resistance'])
                    weights.append(np.mean(rd[r]['yearly_totals']))
            weights = np.array(weights) / sum(weights)
            mean_r = np.average(r_traj, axis=0, weights=weights)
            alpha = 0.3 if len(data[K]) > 1 else 0.8
            ax_d.plot(mean_r, color=K_COLORS[K], alpha=alpha, lw=1)
        # Mean across seeds
        all_r = []
        for seed in SEEDS:
            if seed not in data[K]:
                continue
            rd = data[K][seed]['region_details']
            r_traj = []
            weights = []
            for r in COASTLINE_ORDER:
                if r in rd:
                    r_traj.append(rd[r]['yearly_mean_resistance'])
                    weights.append(np.mean(rd[r]['yearly_totals']))
            weights = np.array(weights) / sum(weights)
            all_r.append(np.average(r_traj, axis=0, weights=weights))
        if all_r:
            mean_all = np.mean(all_r, axis=0)
            ax_d.plot(mean_all, color=K_COLORS[K], lw=2.5, label=f'K={K}')
    ax_d.set_xlabel('Simulation Year')
    ax_d.set_ylabel('Mean Resistance (r̄)')
    ax_d.set_title('D. Resistance Evolution', fontweight='bold')
    ax_d.legend(fontsize=9)

    # Panel E: Va(resistance) dynamics
    ax_e = fig.add_subplot(gs[1, 2])
    for K in K_VALUES:
        all_va = []
        for seed in SEEDS:
            if seed not in data[K]:
                continue
            rd = data[K][seed]['region_details']
            va_traj = []
            weights = []
            for r in COASTLINE_ORDER:
                if r in rd:
                    va_traj.append(rd[r]['yearly_va_resistance'])
                    weights.append(np.mean(rd[r]['yearly_totals']))
            weights = np.array(weights) / sum(weights)
            all_va.append(np.average(va_traj, axis=0, weights=weights))
        if all_va:
            mean_va = np.mean(all_va, axis=0)
            ax_e.plot(mean_va, color=K_COLORS[K], lw=2.5, label=f'K={K}')
    ax_e.set_xlabel('Simulation Year')
    ax_e.set_ylabel('Additive Genetic Variance Va(r)')
    ax_e.set_title('E. Genetic Variance Dynamics', fontweight='bold')
    ax_e.legend(fontsize=9)

    # Panel F: Recovery gradient (mean across seeds)
    ax_f = fig.add_subplot(gs[2, 0:2])
    x = np.arange(len(COASTLINE_ORDER))
    for K in K_VALUES:
        vals = [gradient[K].get(r, 0) * 100 for r in COASTLINE_ORDER]
        ax_f.plot(x, vals, 'o-', color=K_COLORS[K], lw=2, markersize=5,
                 label=f'K={K}', alpha=0.8)
    ax_f.set_xticks(x)
    ax_f.set_xticklabels(COASTLINE_ORDER, rotation=60, ha='right', fontsize=7)
    ax_f.set_ylabel('Recovery (%)')
    ax_f.set_title('F. Recovery Gradient (S→N along coast)', fontweight='bold')
    ax_f.legend(fontsize=9)
    ax_f.axhline(0, color='gray', lw=0.5)

    # Panel G: Inter-seed CV by region (noise profile)
    ax_g = fig.add_subplot(gs[2, 2])
    for K in K_VALUES:
        cvs = []
        for region in COASTLINE_ORDER:
            vals = [data[K][s]['region_recovery'].get(region, 0)
                    for s in SEEDS if s in data[K]]
            if len(vals) >= 2 and np.mean(vals) > 0.001:
                cvs.append(np.std(vals, ddof=1) / np.mean(vals) * 100)
            else:
                cvs.append(0)
        ax_g.plot(range(len(COASTLINE_ORDER)), cvs, 'o-', color=K_COLORS[K],
                 lw=1.5, markersize=4, label=f'K={K}', alpha=0.8)
    ax_g.set_xticks(range(len(COASTLINE_ORDER)))
    ax_g.set_xticklabels(COASTLINE_ORDER, rotation=60, ha='right', fontsize=6)
    ax_g.set_ylabel('Inter-seed CV (%)')
    ax_g.set_title('G. Stochastic Noise by Region', fontweight='bold')
    ax_g.legend(fontsize=9)

    fig.suptitle('K Sensitivity Analysis — Carrying Capacity Effects on Model Behavior',
                 fontsize=16, fontweight='bold', y=0.98)
    plt.savefig(OUTDIR / 'k_sensitivity_summary.png', dpi=200, bbox_inches='tight')
    plt.close()
    print(f"\n  Saved: {OUTDIR / 'k_sensitivity_summary.png'}")


def plot_region_trajectories(data):
    """Per-region population trajectories for each K, all seeds overlaid."""
    fig, axes = plt.subplots(3, 3, figsize=(18, 14), sharex=True)
    axes = axes.flatten()

    for idx, region in enumerate(SCORED[:8]):
        ax = axes[idx]
        for K in K_VALUES:
            for seed in SEEDS:
                if seed not in data[K]:
                    continue
                rd = data[K][seed]['region_details']
                if region in rd:
                    pop = np.array(rd[region]['yearly_totals'])
                    # Normalize by K * n_nodes for comparability
                    n_nodes = rd[region]['n_nodes']
                    pop_norm = pop / (K * n_nodes)
                    ax.plot(pop_norm, color=K_COLORS[K], alpha=0.4, lw=1)
            # Mean
            all_pop = []
            for seed in SEEDS:
                if seed in data[K]:
                    rd = data[K][seed]['region_details']
                    if region in rd:
                        n_nodes = rd[region]['n_nodes']
                        all_pop.append(np.array(rd[region]['yearly_totals']) / (K * n_nodes))
            if all_pop:
                ax.plot(np.mean(all_pop, axis=0), color=K_COLORS[K], lw=2.5,
                       label=f'K={K}')

        target = RECOVERY_TARGETS.get(region, 0)
        ax.axhline(target, color='gray', ls='--', lw=1, alpha=0.5)
        ax.set_title(f'{region} (target={target*100:.1f}%)', fontweight='bold', fontsize=10)
        ax.set_ylabel('Pop / Capacity')
        if idx == 0:
            ax.legend(fontsize=7)

    # Use last subplot for legend
    axes[-1].axis('off')
    for K in K_VALUES:
        axes[-1].plot([], [], color=K_COLORS[K], lw=2.5, label=f'K={K}')
    axes[-1].legend(fontsize=12, loc='center')

    fig.suptitle('Per-Region Population Trajectories by K (normalized)',
                 fontsize=14, fontweight='bold')
    plt.savefig(OUTDIR / 'k_sensitivity_regions.png', dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {OUTDIR / 'k_sensitivity_regions.png'}")


# ── Summary verdict ────────────────────────────────────────────────────
def verdict(data, dyn):
    """Print a summary verdict on minimum acceptable K."""
    print("\n" + "="*70)
    print("VERDICT")
    print("="*70)

    for K in K_VALUES:
        rmsles = [data[K][s]['scoring']['rmsle'] for s in SEEDS if s in data[K]]
        if len(rmsles) < 2:
            continue
        cv = np.std(rmsles, ddof=1) / np.mean(rmsles) * 100

        # Check gradient preservation: ratio of AK-PWS to CA-N recovery
        ak_vals = [data[K][s]['region_recovery'].get('AK-PWS', 0) for s in SEEDS if s in data[K]]
        ca_vals = [data[K][s]['region_recovery'].get('CA-N', 0) for s in SEEDS if s in data[K]]
        ak_mean = np.mean(ak_vals) if ak_vals else 0
        ca_mean = np.mean(ca_vals) if ca_vals else 0
        ratio = ak_mean / ca_mean if ca_mean > 0 else float('inf')

        # Check dynamics shape preservation: crash depth
        crash_depths = []
        for seed in dyn[K]:
            ts = dyn[K][seed]
            peak = ts[:3].max()
            nadir = ts.min()
            crash_depths.append((peak - nadir) / peak * 100)

        print(f"\n  K={K}:")
        print(f"    RMSLE CV:         {cv:.1f}% {'✓' if cv < 5 else '⚠' if cv < 10 else '✗'}")
        print(f"    Gradient (AK/CA): {ratio:.1f}x {'✓' if 3 < ratio < 300 else '⚠'}")
        print(f"    Crash depth:      {np.mean(crash_depths):.1f}% ± {np.std(crash_depths):.1f}%")
        print(f"    Acceptable:       {'YES' if cv < 10 else 'MARGINAL' if cv < 20 else 'NO'}")


# ── Main ────────────────────────────────────────────────────────────────
def main():
    print("K Sensitivity Analysis")
    print(f"Base: {BASE}")
    print(f"Output: {OUTDIR}")

    data = load_all()

    # Check completeness
    for K in K_VALUES:
        n = len(data[K])
        print(f"  K={K}: {n}/{len(SEEDS)} seeds loaded")
        if n == 0:
            print(f"    WARNING: No data for K={K}, skipping")

    rows = analyze_aggregate(data)
    region_data = analyze_regional(data)
    dyn = analyze_dynamics(data)
    evo = analyze_evolution(data)
    gradient = analyze_gradient(data)
    analyze_noise(data)
    analyze_wavefront(data)

    # Figures
    print("\n" + "="*70)
    print("GENERATING FIGURES")
    print("="*70)
    plot_summary(data, dyn, evo, gradient)
    plot_region_trajectories(data)

    verdict(data, dyn)

    # Write machine-readable summary
    summary = {
        'k_values': K_VALUES,
        'seeds': SEEDS,
    }
    for K in K_VALUES:
        rmsles = [data[K][s]['scoring']['rmsle'] for s in SEEDS if s in data[K]]
        if rmsles:
            summary[f'K{K}'] = {
                'rmsle_mean': float(np.mean(rmsles)),
                'rmsle_std': float(np.std(rmsles, ddof=1)) if len(rmsles) >= 2 else 0,
                'rmsle_cv': float(np.std(rmsles, ddof=1) / np.mean(rmsles) * 100) if len(rmsles) >= 2 else 0,
                'rmsle_values': [float(r) for r in rmsles],
                'n_seeds': len(rmsles),
            }
    with open(OUTDIR / 'k_sensitivity_summary.json', 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"\n  Saved: {OUTDIR / 'k_sensitivity_summary.json'}")


if __name__ == '__main__':
    main()
