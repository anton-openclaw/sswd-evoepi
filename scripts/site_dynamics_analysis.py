#!/usr/bin/env python3
"""
Site-level population dynamics analysis for SSWD forecast simulations.

Key question: During regional recovery booms, do individual "oasis" sites
drive the recovery, or do many sites recover together?
"""

import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.colors import LogNorm
from collections import defaultdict
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Try seaborn style
try:
    plt.style.use('seaborn-v0_8-paper')
except:
    try:
        plt.style.use('seaborn-paper')
    except:
        pass

# Paths
BASE = Path('/home/starbot/.openclaw/workspace/sswd-evoepi')
DATA_DIR = BASE / 'results/calibration/F01'
OUT_DATA = BASE / 'reports/publishable/data'
OUT_FIG = BASE / 'reports/publishable/figures'
OUT_DATA.mkdir(parents=True, exist_ok=True)
OUT_FIG.mkdir(parents=True, exist_ok=True)

SEEDS = [42, 123, 999]
KEY_REGIONS = ['AK-PWS', 'AK-FN', 'BC-N', 'OR', 'CA-C']
ALL_REGIONS = ['AK-PWS', 'AK-FN', 'BC-N', 'OR', 'CA-C', 'CA-N']

# Color scheme
REGION_COLORS = {
    'AK-PWS': '#1b9e77',
    'AK-FN': '#d95f02',
    'BC-N': '#7570b3',
    'OR': '#e7298a',
    'CA-C': '#66a61e',
    'CA-N': '#e6ab02',
}

###############################################################################
# Data loading
###############################################################################

def load_seed(seed):
    """Load data for a given seed."""
    path = DATA_DIR / f'monthly_seed{seed}.npz'
    d = np.load(path, allow_pickle=True)
    return {
        'sim_days': d['sim_days'],
        'populations': d['populations'].astype(np.float64),
        'infected': d['infected'].astype(np.float64),
        'site_lats': d['site_lats'],
        'site_lons': d['site_lons'],
        'site_names': np.array([str(n) for n in d['site_names']]),
        'K': int(d['K']),
        'sst_start_year': int(d['sst_start_year']),
    }

def get_region_prefix(name):
    """Extract region prefix from site name."""
    parts = name.split('-')
    if len(parts) >= 3 and parts[0] == 'AK':
        return parts[0] + '-' + parts[1]
    elif len(parts) >= 3 and parts[0] == 'BC':
        return parts[0] + '-' + parts[1]
    elif len(parts) >= 3 and parts[0] == 'CA':
        return parts[0] + '-' + parts[1]
    elif len(parts) >= 3 and parts[0] == 'SS':
        return parts[0] + '-' + parts[1]
    elif len(parts) >= 2 and parts[0] == 'WA':
        return parts[0] + '-' + parts[1]
    else:
        return parts[0]

def build_region_indices(site_names):
    """Map region prefix -> list of site indices."""
    region_map = defaultdict(list)
    for i, name in enumerate(site_names):
        prefix = get_region_prefix(name)
        region_map[prefix].append(i)
    return dict(region_map)

def sim_days_to_dates(sim_days, start_year=2012):
    """Convert simulation days to approximate dates."""
    from datetime import datetime, timedelta
    base = datetime(start_year, 1, 1)
    return [base + timedelta(days=int(d)) for d in sim_days]

def sim_days_to_years(sim_days, start_year=2012):
    """Convert simulation days to fractional years."""
    return start_year + sim_days / 365.25


###############################################################################
# Analysis functions
###############################################################################

def gini_coefficient(values):
    """Compute Gini coefficient for an array of non-negative values."""
    v = np.array(values, dtype=np.float64)
    v = v[v >= 0]
    if len(v) == 0 or v.sum() == 0:
        return 0.0
    v = np.sort(v)
    n = len(v)
    index = np.arange(1, n + 1)
    return (2.0 * np.sum(index * v) / (n * np.sum(v))) - (n + 1.0) / n

def top_pct_concentration(values, pct=0.05):
    """Fraction of total held by top pct% of sites."""
    v = np.array(values, dtype=np.float64)
    if v.sum() == 0:
        return 0.0
    n_top = max(1, int(np.ceil(len(v) * pct)))
    sorted_v = np.sort(v)[::-1]
    return sorted_v[:n_top].sum() / v.sum()

def identify_booms(regional_pop, window=6, threshold=0.2):
    """
    Identify recovery booms: periods where regional population increases by >threshold
    within a `window`-timestep window (~6 months).
    
    Returns list of dicts with start_idx, end_idx, magnitude, rel_increase.
    """
    booms = []
    n = len(regional_pop)
    i = 0
    while i < n - window:
        start_pop = regional_pop[i]
        if start_pop <= 0:
            i += 1
            continue
        # Look for max increase in windows starting at i
        best_end = -1
        best_increase = 0
        for j in range(i + 1, min(i + window + 1, n)):
            end_pop = regional_pop[j]
            rel_increase = (end_pop - start_pop) / start_pop
            if rel_increase > best_increase:
                best_increase = rel_increase
                best_end = j
        if best_increase > threshold and best_end > 0:
            booms.append({
                'start_idx': int(i),
                'end_idx': int(best_end),
                'start_pop': float(start_pop),
                'end_pop': float(regional_pop[best_end]),
                'abs_increase': float(regional_pop[best_end] - start_pop),
                'rel_increase': float(best_increase),
            })
            i = best_end + 1  # Skip past this boom
        else:
            i += 1
    return booms

def decompose_boom(populations, site_indices, start_idx, end_idx):
    """
    Decompose a boom into site-level contributions.
    Returns sorted list of (site_idx, contribution, site_pop_start, site_pop_end).
    """
    contributions = []
    for si in site_indices:
        delta = populations[end_idx, si] - populations[start_idx, si]
        contributions.append({
            'site_idx': int(si),
            'delta': float(delta),
            'pop_start': float(populations[start_idx, si]),
            'pop_end': float(populations[end_idx, si]),
        })
    # Sort by contribution descending
    contributions.sort(key=lambda x: x['delta'], reverse=True)
    return contributions


###############################################################################
# Main analysis
###############################################################################

print("Loading data for all seeds...")
all_data = {}
for seed in SEEDS:
    all_data[seed] = load_seed(seed)
    print(f"  Seed {seed}: populations shape {all_data[seed]['populations'].shape}")

# Use seed 42 as primary
d = all_data[42]
site_names = d['site_names']
region_map = build_region_indices(site_names)
sim_days = d['sim_days']
years = sim_days_to_years(sim_days, d['sst_start_year'])
dates = sim_days_to_dates(sim_days, d['sst_start_year'])
populations = d['populations']

print(f"\nRegions found: {sorted(region_map.keys())}")
for r in KEY_REGIONS:
    print(f"  {r}: {len(region_map.get(r, []))} sites")

# ============================================================================
# Task 1: Regional population time series and boom identification
# ============================================================================
print("\n=== Task 1: Identifying boom-bust cycles ===")

regional_booms = {}
regional_pops = {}

for region in KEY_REGIONS:
    indices = region_map[region]
    reg_pop = populations[:, indices].sum(axis=1)
    regional_pops[region] = reg_pop
    booms = identify_booms(reg_pop, window=6, threshold=0.2)
    regional_booms[region] = booms
    print(f"\n  {region}: {len(booms)} booms detected")
    for b in booms[:5]:
        y_start = years[b['start_idx']]
        y_end = years[b['end_idx']]
        print(f"    {y_start:.1f}-{y_end:.1f}: +{b['rel_increase']*100:.0f}% "
              f"({b['start_pop']:.0f} -> {b['end_pop']:.0f})")

# ============================================================================
# Task 2: Decompose booms into site contributions
# ============================================================================
print("\n=== Task 2: Boom decomposition ===")

boom_decompositions = {}

for region in KEY_REGIONS:
    indices = region_map[region]
    booms = regional_booms[region]
    decomps = []
    for b in booms:
        contribs = decompose_boom(populations, indices, b['start_idx'], b['end_idx'])
        total_increase = b['abs_increase']
        
        # Sites with positive contributions
        positive_sites = [c for c in contribs if c['delta'] > 0]
        n_positive = len(positive_sites)
        n_total = len(indices)
        
        # Top-N concentration
        if total_increase > 0:
            top1_frac = sum(c['delta'] for c in contribs[:1]) / total_increase if contribs else 0
            top3_frac = sum(c['delta'] for c in contribs[:3]) / total_increase if contribs else 0
            top5_frac = sum(c['delta'] for c in contribs[:5]) / total_increase if contribs else 0
        else:
            top1_frac = top3_frac = top5_frac = 0
        
        decomps.append({
            'boom': b,
            'n_positive_sites': n_positive,
            'n_total_sites': n_total,
            'frac_sites_growing': n_positive / n_total if n_total > 0 else 0,
            'top1_frac': float(top1_frac),
            'top3_frac': float(top3_frac),
            'top5_frac': float(top5_frac),
            'top_sites': [
                {
                    'name': site_names[c['site_idx']],
                    'idx': c['site_idx'],
                    'delta': c['delta'],
                    'frac': c['delta'] / total_increase if total_increase > 0 else 0,
                }
                for c in contribs[:10]
            ],
        })
    boom_decompositions[region] = decomps
    
    if decomps:
        avg_top1 = np.mean([d['top1_frac'] for d in decomps])
        avg_top3 = np.mean([d['top3_frac'] for d in decomps])
        avg_top5 = np.mean([d['top5_frac'] for d in decomps])
        avg_frac_growing = np.mean([d['frac_sites_growing'] for d in decomps])
        print(f"\n  {region} ({len(decomps)} booms):")
        print(f"    Avg sites growing: {avg_frac_growing*100:.0f}%")
        print(f"    Avg top-1 site share: {avg_top1*100:.0f}%")
        print(f"    Avg top-3 sites share: {avg_top3*100:.0f}%")
        print(f"    Avg top-5 sites share: {avg_top5*100:.0f}%")

# ============================================================================
# Task 3: Identify oasis sites (repeatedly in top contributors)
# ============================================================================
print("\n=== Task 3: Oasis site identification ===")

oasis_sites = {}

for region in KEY_REGIONS:
    decomps = boom_decompositions[region]
    site_counts = defaultdict(int)
    site_total_contrib = defaultdict(float)
    
    for dec in decomps:
        for ts in dec['top_sites'][:5]:  # top 5 contributors per boom
            site_counts[ts['name']] += 1
            site_total_contrib[ts['name']] += ts['frac']
    
    n_booms = len(decomps)
    if n_booms == 0:
        continue
    
    # Sites appearing in top-5 in at least 30% of booms
    oasis_threshold = max(2, int(0.3 * n_booms))
    oases = [(name, count, site_total_contrib[name] / count)
             for name, count in site_counts.items()
             if count >= oasis_threshold]
    oases.sort(key=lambda x: x[1], reverse=True)
    
    oasis_sites[region] = oases
    print(f"\n  {region} oasis sites (appear in top-5 of ≥{oasis_threshold}/{n_booms} booms):")
    for name, count, avg_frac in oases[:10]:
        idx = np.where(site_names == name)[0][0]
        lat, lon = d['site_lats'][idx], d['site_lons'][idx]
        mean_pop = populations[:, idx].mean()
        print(f"    {name}: {count}/{n_booms} booms, avg share {avg_frac*100:.1f}%, "
              f"mean pop {mean_pop:.0f}, lat={lat:.2f}, lon={lon:.2f}")

# ============================================================================
# Task 4: Gini coefficient and concentration metrics
# ============================================================================
print("\n=== Task 4: Concentration metrics ===")

gini_trajectories = {}
top5pct_trajectories = {}

for region in ALL_REGIONS:
    indices = region_map[region]
    ginis = []
    top5pcts = []
    for t in range(len(sim_days)):
        site_pops = populations[t, indices]
        ginis.append(gini_coefficient(site_pops))
        top5pcts.append(top_pct_concentration(site_pops, pct=0.05))
    gini_trajectories[region] = ginis
    top5pct_trajectories[region] = top5pcts
    
    # Report early vs late
    early = np.mean(ginis[:24])  # first 2 years
    late = np.mean(ginis[-24:])  # last 2 years
    print(f"  {region}: Gini early={early:.3f}, late={late:.3f}, "
          f"change={late-early:+.3f}")

# ============================================================================
# Task 5: Cross-seed validation
# ============================================================================
print("\n=== Task 5: Cross-seed validation ===")

# For each seed, identify top oasis sites per region
seed_oases = {}

for seed in SEEDS:
    sd = all_data[seed]
    pop = sd['populations'].astype(np.float64)
    sn = np.array([str(n) for n in sd['site_names']])
    rm = build_region_indices(sn)
    
    seed_oases[seed] = {}
    for region in KEY_REGIONS:
        indices = rm[region]
        reg_pop = pop[:, indices].sum(axis=1)
        booms = identify_booms(reg_pop, window=6, threshold=0.2)
        
        site_counts = defaultdict(int)
        for b in booms:
            contribs = decompose_boom(pop, indices, b['start_idx'], b['end_idx'])
            for c in contribs[:5]:
                site_counts[sn[c['site_idx']]] += 1
        
        n_booms = len(booms)
        oasis_threshold = max(2, int(0.3 * n_booms))
        top_sites = sorted(site_counts.items(), key=lambda x: x[1], reverse=True)
        oasis_names = [name for name, count in top_sites if count >= oasis_threshold]
        seed_oases[seed][region] = set(oasis_names)

# Cross-seed consistency
print("\nCross-seed oasis consistency:")
oasis_consistency = {}
for region in KEY_REGIONS:
    sets = [seed_oases[s].get(region, set()) for s in SEEDS]
    union = sets[0] | sets[1] | sets[2]
    intersection = sets[0] & sets[1] & sets[2]
    if len(union) > 0:
        jaccard = len(intersection) / len(union)
    else:
        jaccard = 0
    
    # Count how many seeds each site appears in
    site_seed_count = {}
    for s_set, seed in zip(sets, SEEDS):
        for name in s_set:
            site_seed_count[name] = site_seed_count.get(name, 0) + 1
    
    consistent = [name for name, cnt in site_seed_count.items() if cnt >= 2]
    all_three = [name for name, cnt in site_seed_count.items() if cnt == 3]
    
    oasis_consistency[region] = {
        'union_size': len(union),
        'intersection_size': len(intersection),
        'jaccard': jaccard,
        'consistent_2of3': consistent,
        'consistent_3of3': all_three,
        'per_seed': {s: sorted(seed_oases[s].get(region, set())) for s in SEEDS},
    }
    
    print(f"\n  {region}:")
    print(f"    Union: {len(union)} sites, Intersection: {len(intersection)} sites")
    print(f"    Jaccard: {jaccard:.2f}")
    print(f"    In ≥2 seeds: {len(consistent)} sites: {sorted(consistent)[:10]}")
    print(f"    In all 3 seeds: {len(all_three)} sites: {sorted(all_three)[:10]}")


###############################################################################
# Save JSON data
###############################################################################
print("\n=== Saving JSON data ===")

json_output = {
    'metadata': {
        'seeds': SEEDS,
        'key_regions': KEY_REGIONS,
        'n_timesteps': int(len(sim_days)),
        'n_sites': int(len(site_names)),
        'start_year': int(d['sst_start_year']),
        'K': int(d['K']),
    },
    'booms': {},
    'boom_decompositions': {},
    'oasis_sites': {},
    'gini_trajectories': {},
    'top5pct_trajectories': {},
    'oasis_consistency': {},
    'years': [float(y) for y in years],
}

for region in KEY_REGIONS:
    json_output['booms'][region] = [
        {
            'start_timestep': b['start_idx'],
            'end_timestep': b['end_idx'],
            'start_year': float(years[b['start_idx']]),
            'end_year': float(years[b['end_idx']]),
            'magnitude': b['rel_increase'],
            'abs_increase': b['abs_increase'],
        }
        for b in regional_booms[region]
    ]
    
    json_output['boom_decompositions'][region] = [
        {
            'boom_start_year': float(years[dec['boom']['start_idx']]),
            'boom_end_year': float(years[dec['boom']['end_idx']]),
            'n_positive_sites': dec['n_positive_sites'],
            'frac_sites_growing': dec['frac_sites_growing'],
            'top1_frac': dec['top1_frac'],
            'top3_frac': dec['top3_frac'],
            'top5_frac': dec['top5_frac'],
            'top_contributing_sites': dec['top_sites'][:10],
        }
        for dec in boom_decompositions.get(region, [])
    ]
    
    json_output['oasis_sites'][region] = [
        {
            'name': name,
            'n_booms_in_top5': count,
            'avg_share': avg_frac,
            'lat': float(d['site_lats'][np.where(site_names == name)[0][0]]),
            'lon': float(d['site_lons'][np.where(site_names == name)[0][0]]),
        }
        for name, count, avg_frac in oasis_sites.get(region, [])
    ]

for region in ALL_REGIONS:
    # Subsample trajectories to avoid huge JSON (every 6 months)
    step = 6
    json_output['gini_trajectories'][region] = [
        float(g) for g in gini_trajectories[region][::step]
    ]
    json_output['top5pct_trajectories'][region] = [
        float(t) for t in top5pct_trajectories[region][::step]
    ]

# Convert sets to lists for JSON serialization
for region in KEY_REGIONS:
    oc = oasis_consistency.get(region, {})
    json_output['oasis_consistency'][region] = {
        'union_size': oc.get('union_size', 0),
        'intersection_size': oc.get('intersection_size', 0),
        'jaccard': oc.get('jaccard', 0),
        'consistent_2of3': sorted(oc.get('consistent_2of3', [])),
        'consistent_3of3': sorted(oc.get('consistent_3of3', [])),
        'per_seed': {str(k): v for k, v in oc.get('per_seed', {}).items()},
    }

with open(OUT_DATA / 'site_dynamics.json', 'w') as f:
    json.dump(json_output, f, indent=2)
print(f"Saved JSON to {OUT_DATA / 'site_dynamics.json'}")


###############################################################################
# FIGURE 1: Site decomposition heatmaps
###############################################################################
print("\n=== Generating Figure 1: Site Decomposition ===")

fig1_regions = ['AK-PWS', 'OR', 'CA-C', 'BC-N']
fig, axes = plt.subplots(len(fig1_regions), 1, figsize=(14, 16), sharex=True)

for ax, region in zip(axes, fig1_regions):
    indices = region_map[region]
    site_pops = populations[:, indices]  # (T, n_sites)
    n_sites = len(indices)
    
    # Sort sites by their mean population (descending) for visual clarity
    mean_pops = site_pops.mean(axis=0)
    sort_order = np.argsort(mean_pops)[::-1]
    site_pops_sorted = site_pops[:, sort_order]
    
    # Create heatmap
    im = ax.imshow(
        site_pops_sorted.T,
        aspect='auto',
        cmap='viridis',
        interpolation='nearest',
        extent=[years[0], years[-1], n_sites - 0.5, -0.5],
        vmin=0,
        vmax=d['K'],
    )
    
    ax.set_ylabel(f'{region}\n({n_sites} sites)', fontsize=11, fontweight='bold')
    ax.set_yticks([0, n_sites // 2, n_sites - 1])
    ax.set_yticklabels(['Top', 'Mid', 'Bottom'], fontsize=9)
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, label='Population', shrink=0.8, pad=0.02)
    cbar.ax.tick_params(labelsize=8)

axes[-1].set_xlabel('Year', fontsize=12)
axes[0].set_title('Site-Level Population Dynamics: Are Recoveries Concentrated or Distributed?',
                   fontsize=13, fontweight='bold', pad=15)

plt.tight_layout()
fig.savefig(OUT_FIG / 'fig_site_decomposition.png', dpi=300, bbox_inches='tight')
plt.close(fig)
print("  Saved fig_site_decomposition.png")


###############################################################################
# FIGURE 2: Oasis concentration (Gini + top-5%)
###############################################################################
print("\n=== Generating Figure 2: Oasis Concentration ===")

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

for region in ALL_REGIONS:
    color = REGION_COLORS.get(region, '#333333')
    ax1.plot(years, gini_trajectories[region], label=region, color=color, linewidth=1.5, alpha=0.85)
    ax2.plot(years, top5pct_trajectories[region], label=region, color=color, linewidth=1.5, alpha=0.85)

ax1.set_ylabel('Gini Coefficient', fontsize=12)
ax1.set_title('(a) Population Concentration: Gini Coefficient Over Time', fontsize=12, fontweight='bold')
ax1.legend(loc='upper left', fontsize=9, ncol=2, framealpha=0.9)
ax1.set_ylim(0, 1)
ax1.grid(True, alpha=0.3)
ax1.axhline(y=0.5, color='gray', linestyle='--', alpha=0.4, linewidth=0.8)

ax2.set_ylabel('Top 5% Site Share', fontsize=12)
ax2.set_xlabel('Year', fontsize=12)
ax2.set_title('(b) Fraction of Regional Population in Top 5% of Sites', fontsize=12, fontweight='bold')
ax2.legend(loc='upper left', fontsize=9, ncol=2, framealpha=0.9)
ax2.set_ylim(0, 1)
ax2.grid(True, alpha=0.3)

plt.tight_layout()
fig.savefig(OUT_FIG / 'fig_oasis_concentration.png', dpi=300, bbox_inches='tight')
plt.close(fig)
print("  Saved fig_oasis_concentration.png")


###############################################################################
# FIGURE 3: Boom decomposition bar charts
###############################################################################
print("\n=== Generating Figure 3: Boom Decomposition ===")

# Pick 2-3 interesting booms from different regions
selected_booms = []

for region in ['AK-PWS', 'OR', 'CA-C']:
    decomps = boom_decompositions.get(region, [])
    if decomps:
        # Pick the largest boom
        biggest = max(decomps, key=lambda d: d['boom']['abs_increase'])
        selected_booms.append((region, biggest))

n_panels = len(selected_booms)
fig, axes = plt.subplots(1, n_panels, figsize=(5 * n_panels, 6))
if n_panels == 1:
    axes = [axes]

for ax, (region, dec) in zip(axes, selected_booms):
    boom = dec['boom']
    y_start = years[boom['start_idx']]
    y_end = years[boom['end_idx']]
    
    # Get all site contributions
    indices = region_map[region]
    contribs = decompose_boom(populations, indices, boom['start_idx'], boom['end_idx'])
    total_increase = boom['abs_increase']
    
    # Show top 15 sites + "others"
    n_show = min(15, len(contribs))
    top_contribs = contribs[:n_show]
    other_delta = sum(c['delta'] for c in contribs[n_show:])
    
    names = [site_names[c['site_idx']].split('-')[-1] for c in top_contribs]
    fracs = [c['delta'] / total_increase if total_increase > 0 else 0 for c in top_contribs]
    
    if other_delta != 0 and total_increase > 0:
        names.append('Others')
        fracs.append(other_delta / total_increase)
    
    colors_bar = []
    for f in fracs:
        if f > 0:
            colors_bar.append(REGION_COLORS.get(region, '#1f77b4'))
        else:
            colors_bar.append('#cc4444')
    
    bars = ax.barh(range(len(names)), fracs, color=colors_bar, alpha=0.8, edgecolor='white', linewidth=0.5)
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel('Fraction of Regional Increase', fontsize=10)
    ax.set_title(f'{region}\n{y_start:.1f}–{y_end:.1f} (+{boom["rel_increase"]*100:.0f}%)',
                 fontsize=11, fontweight='bold')
    ax.axvline(x=0, color='black', linewidth=0.5)
    ax.grid(True, axis='x', alpha=0.3)
    
    # Annotate top1/3/5 concentration
    ax.text(0.98, 0.98, 
            f'Top 1: {dec["top1_frac"]*100:.0f}%\n'
            f'Top 3: {dec["top3_frac"]*100:.0f}%\n'
            f'Top 5: {dec["top5_frac"]*100:.0f}%\n'
            f'Growing: {dec["frac_sites_growing"]*100:.0f}%',
            transform=ax.transAxes, fontsize=8,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='wheat', alpha=0.8))

fig.suptitle('Boom Decomposition: Which Sites Drive Recovery?', fontsize=13, fontweight='bold', y=1.02)
plt.tight_layout()
fig.savefig(OUT_FIG / 'fig_boom_decomposition.png', dpi=300, bbox_inches='tight')
plt.close(fig)
print("  Saved fig_boom_decomposition.png")


###############################################################################
# FIGURE 4: Oasis consistency across seeds
###############################################################################
print("\n=== Generating Figure 4: Oasis Consistency ===")

# Find sites that are oases in at least 2 of 3 seeds
consistent_oases = {}
for region in KEY_REGIONS:
    oc = oasis_consistency[region]
    consistent_oases[region] = oc.get('consistent_2of3', [])

# Pick top oasis sites from 3 regions
plot_regions = ['AK-PWS', 'OR', 'CA-C']
sites_per_region = 4

fig, axes = plt.subplots(len(plot_regions), 1, figsize=(14, 4 * len(plot_regions)), sharex=True)

for ax, region in zip(axes, plot_regions):
    oases = consistent_oases.get(region, [])
    
    # If no consistent oases, use the top oasis sites from seed 42
    if not oases:
        oases = [o[0] for o in oasis_sites.get(region, [])[:sites_per_region]]
    else:
        oases = oases[:sites_per_region]
    
    if not oases:
        ax.text(0.5, 0.5, f'No oasis sites identified for {region}',
                transform=ax.transAxes, ha='center', va='center')
        continue
    
    linestyles = ['-', '--', ':']
    seed_colors = {42: '#1f77b4', 123: '#ff7f0e', 999: '#2ca02c'}
    
    for i, site_name in enumerate(oases):
        for seed in SEEDS:
            sd = all_data[seed]
            sn = np.array([str(n) for n in sd['site_names']])
            idx = np.where(sn == site_name)[0]
            if len(idx) == 0:
                continue
            idx = idx[0]
            pop = sd['populations'][:, idx].astype(float)
            
            label = f'{site_name} (seed {seed})' if i == 0 else (f'seed {seed}' if i == 1 and seed != 42 else None)
            alpha = 0.9 if seed == 42 else 0.5
            lw = 1.8 if seed == 42 else 1.0
            
            # Use consistent color per site, different linestyle per seed
            site_color = plt.cm.Set2(i / max(1, len(oases) - 1)) if len(oases) > 1 else REGION_COLORS.get(region, '#1f77b4')
            ax.plot(years, pop, color=site_color, linestyle=linestyles[SEEDS.index(seed)],
                    alpha=alpha, linewidth=lw)
    
    ax.set_ylabel('Population', fontsize=11)
    ax.set_title(f'{region} — Top Oasis Sites Across Seeds (solid=seed42, dashed=seed123, dotted=seed999)',
                 fontsize=11, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.set_ylim(bottom=0)
    
    # Add site name annotations
    for i, site_name in enumerate(oases):
        site_color = plt.cm.Set2(i / max(1, len(oases) - 1)) if len(oases) > 1 else REGION_COLORS.get(region, '#1f77b4')
        ax.annotate(site_name, xy=(0.01, 0.95 - 0.08 * i), xycoords='axes fraction',
                    fontsize=8, color=site_color, fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7))

axes[-1].set_xlabel('Year', fontsize=12)
fig.suptitle('Oasis Site Consistency: Same Sites Across Stochastic Replicates?',
             fontsize=13, fontweight='bold', y=1.01)
plt.tight_layout()
fig.savefig(OUT_FIG / 'fig_oasis_consistency.png', dpi=300, bbox_inches='tight')
plt.close(fig)
print("  Saved fig_oasis_consistency.png")


###############################################################################
# Summary
###############################################################################
print("\n" + "=" * 70)
print("ANALYSIS COMPLETE")
print("=" * 70)
print(f"\nJSON data: {OUT_DATA / 'site_dynamics.json'}")
print(f"Figures:   {OUT_FIG}/")
print(f"  - fig_site_decomposition.png")
print(f"  - fig_oasis_concentration.png")
print(f"  - fig_boom_decomposition.png")
print(f"  - fig_oasis_consistency.png")

# Quick summary of key finding
print("\n=== KEY FINDINGS ===")
for region in KEY_REGIONS:
    decomps = boom_decompositions.get(region, [])
    if decomps:
        avg_top3 = np.mean([d['top3_frac'] for d in decomps])
        avg_growing = np.mean([d['frac_sites_growing'] for d in decomps])
        oc = oasis_consistency.get(region, {})
        n_consistent = len(oc.get('consistent_3of3', []))
        print(f"\n  {region}:")
        print(f"    Recovery pattern: Top 3 sites drive {avg_top3*100:.0f}% of boom magnitude")
        print(f"    {avg_growing*100:.0f}% of sites show growth during booms")
        print(f"    {n_consistent} sites are oases across all 3 seeds")
