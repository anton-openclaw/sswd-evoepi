#!/usr/bin/env python3
"""Generate 10 spatial-temporal dynamics figures for W249-W268 calibration report.

All figures use W257 (best config) NPZ data and focus on model dynamics.
"""
import sys
sys.path.insert(0, '/home/starbot/.openclaw/workspace/sswd-evoepi')

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from matplotlib.ticker import MaxNLocator
from pathlib import Path

from viz.calibration import apply_pub_style
apply_pub_style()

# Bump fonts for publication quality
matplotlib.rcParams.update({
    'font.size': 12,
    'axes.labelsize': 13,
    'axes.titlesize': 14,
    'xtick.labelsize': 11,
    'ytick.labelsize': 11,
    'legend.fontsize': 10,
})

# ── Load data ──────────────────────────────────────────────────────────
NPZ = '/home/starbot/.openclaw/workspace/sswd-evoepi/reports/w249-w268-comprehensive/data/sweep/W257/monthly_seed42.npz'
OUTDIR = Path('/home/starbot/.openclaw/workspace/sswd-evoepi/reports/w249-w268-comprehensive/figures')
OUTDIR.mkdir(parents=True, exist_ok=True)

d = np.load(NPZ, allow_pickle=True)
sim_days   = d['sim_days']          # (159,)
populations = d['populations']      # (159, 896)
infected   = d['infected']          # (159, 896)
site_lats  = d['site_lats']         # (896,)
site_lons  = d['site_lons']         # (896,)
site_names = d['site_names']        # (896,) str
K          = int(d['K'])            # 5000
sst_start_year = int(d['sst_start_year'])  # 2012

years = sst_start_year + sim_days / 365.25
n_times, n_sites = populations.shape

# ── Region mapping ─────────────────────────────────────────────────────
regions = np.array([str(n).rsplit('-', 1)[0] for n in site_names])
region_order = ['AK-AL', 'AK-WG', 'AK-OC', 'AK-EG', 'AK-PWS', 'AK-FN', 'AK-FS',
                'BC-N', 'BC-C', 'SS-N', 'SS-S', 'JDF', 'WA-O', 'OR', 'CA-N', 'CA-C', 'CA-S', 'BJ']
scored_regions = ['AK-PWS', 'AK-FN', 'AK-FS', 'BC-N', 'SS-S', 'JDF', 'OR', 'CA-N']

# Latitude-based colors: cool (north) to warm (south)
region_colors = {r: cm.coolwarm(i / (len(region_order) - 1)) for i, r in enumerate(region_order)}

# Site indices per region
region_sites = {r: np.where(regions == r)[0] for r in region_order}

# Helper: find nearest time index for a given year
def tidx_for_year(target_year):
    return np.argmin(np.abs(years - target_year))

DPI = 300

# ═══════════════════════════════════════════════════════════════════════
# Fig 1: Population heatmap (time × region)
# ═══════════════════════════════════════════════════════════════════════
print("Fig 1: Population heatmap...")
fig, ax = plt.subplots(figsize=(12, 8))

# Build matrix: (n_regions, n_times) = mean relative population per region
heatmap_data = np.zeros((len(region_order), n_times))
for i, r in enumerate(region_order):
    idx = region_sites[r]
    if len(idx) > 0:
        # Sum pop across sites, divide by n_sites * K
        heatmap_data[i, :] = populations[:, idx].sum(axis=1) / (len(idx) * K)

im = ax.imshow(heatmap_data, aspect='auto', cmap='RdYlGn', vmin=0, vmax=1,
               origin='upper', interpolation='nearest')

# X-axis: years
year_ticks = np.arange(2012, 2026, 1)
tick_positions = [tidx_for_year(y) for y in year_ticks]
ax.set_xticks(tick_positions)
ax.set_xticklabels([str(int(y)) for y in year_ticks], rotation=45, ha='right')
ax.set_xlabel('Year')

# Y-axis: regions (north=top)
ax.set_yticks(range(len(region_order)))
ax.set_yticklabels(region_order)
ax.set_ylabel('Region (north → south)')

cbar = plt.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
cbar.set_label('Relative population (pop / K)', fontsize=12)
ax.set_title('Population dynamics: south-to-north crash wave', fontsize=14, fontweight='bold')

fig.tight_layout()
fig.savefig(OUTDIR / 'fig_st_population_heatmap.png', dpi=DPI, bbox_inches='tight', facecolor='white')
plt.close(fig)
print("  ✓ Saved fig_st_population_heatmap.png")


# ═══════════════════════════════════════════════════════════════════════
# Fig 2: Disease wavefront geographic snapshots
# ═══════════════════════════════════════════════════════════════════════
print("Fig 2: Wavefront snapshots...")
snapshot_years = [2013, 2015, 2019, 2025]
snapshot_labels = ['(a) Year 1 – Early spread', '(b) Year 3 – Peak southern crash',
                   '(c) Year 7 – Mid-recovery', '(d) Year 13 – Final state']

fig, axes = plt.subplots(2, 2, figsize=(14, 16))
axes = axes.flatten()

for panel_idx, (yr, label) in enumerate(zip(snapshot_years, snapshot_labels)):
    ax = axes[panel_idx]
    ti = tidx_for_year(yr)
    
    pop_t = populations[ti, :].astype(float)
    inf_t = infected[ti, :].astype(float)
    
    # Prevalence
    with np.errstate(divide='ignore', invalid='ignore'):
        prevalence = np.where(pop_t > 0, inf_t / pop_t, np.nan)
    
    # Sites with pop=0: small gray dots
    zero_mask = pop_t == 0
    live_mask = ~zero_mask
    
    # Size: proportional to population (scale for visibility)
    sizes_live = np.clip(pop_t[live_mask] / K * 80, 5, 80)
    
    # Plot zero-pop sites first (background)
    if zero_mask.any():
        ax.scatter(site_lons[zero_mask], site_lats[zero_mask], s=3, c='gray', alpha=0.4, zorder=1)
    
    # Plot live sites
    sc = ax.scatter(site_lons[live_mask], site_lats[live_mask], s=sizes_live,
                    c=prevalence[live_mask], cmap='Reds', vmin=0, vmax=0.3,
                    edgecolors='k', linewidths=0.2, alpha=0.8, zorder=2)
    
    ax.set_title(f'{label} ({yr})', fontsize=13, fontweight='bold')
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_xlim(site_lons.min() - 2, site_lons.max() + 2)
    ax.set_ylim(site_lats.min() - 1, site_lats.max() + 1)
    
    cbar = plt.colorbar(sc, ax=ax, shrink=0.7, pad=0.02)
    cbar.set_label('Infection prevalence')

fig.suptitle('Disease wavefront propagation', fontsize=16, fontweight='bold', y=0.98)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(OUTDIR / 'fig_st_wavefront_snapshots.png', dpi=DPI, bbox_inches='tight', facecolor='white')
plt.close(fig)
print("  ✓ Saved fig_st_wavefront_snapshots.png")


# ═══════════════════════════════════════════════════════════════════════
# Fig 3: Disease arrival timing map
# ═══════════════════════════════════════════════════════════════════════
print("Fig 3: Arrival timing map...")
fig, ax = plt.subplots(figsize=(10, 10))

# First month where infected > 0 for each site
arrival_month = np.full(n_sites, np.nan)
for s in range(n_sites):
    inf_col = infected[:, s]
    where_inf = np.where(inf_col > 0)[0]
    if len(where_inf) > 0:
        arrival_month[s] = where_inf[0]  # index = month number

never_infected = np.isnan(arrival_month)
ever_infected = ~never_infected

# Gray for never-infected
if never_infected.any():
    ax.scatter(site_lons[never_infected], site_lats[never_infected], 
               s=15, c='lightgray', edgecolors='gray', linewidths=0.3, alpha=0.6, label='Never infected', zorder=1)

sc = ax.scatter(site_lons[ever_infected], site_lats[ever_infected], s=25,
                c=arrival_month[ever_infected], cmap='plasma', 
                edgecolors='k', linewidths=0.2, alpha=0.85, zorder=2)

cbar = plt.colorbar(sc, ax=ax, shrink=0.7, pad=0.02)
cbar.set_label('Months after introduction', fontsize=12)

ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_title('Disease arrival timing', fontsize=14, fontweight='bold')
ax.legend(loc='lower left', fontsize=10)
ax.set_xlim(site_lons.min() - 2, site_lons.max() + 2)
ax.set_ylim(site_lats.min() - 1, site_lats.max() + 1)

fig.tight_layout()
fig.savefig(OUTDIR / 'fig_st_arrival_timing.png', dpi=DPI, bbox_inches='tight', facecolor='white')
plt.close(fig)
print("  ✓ Saved fig_st_arrival_timing.png")


# ═══════════════════════════════════════════════════════════════════════
# Fig 4: Population crash magnitude map
# ═══════════════════════════════════════════════════════════════════════
print("Fig 4: Crash magnitude map...")
fig, ax = plt.subplots(figsize=(10, 10))

initial_pop = populations[0, :].astype(float)
min_pop = populations.min(axis=0).astype(float)

# Crash magnitude: 1 - min_pop/initial_pop
crash_mag = np.where(initial_pop > 0, 1.0 - min_pop / initial_pop, np.nan)

# Size proportional to initial population
sizes = np.clip(initial_pop / K * 60, 5, 60)

valid = ~np.isnan(crash_mag)
sc = ax.scatter(site_lons[valid], site_lats[valid], s=sizes[valid],
                c=crash_mag[valid], cmap='Reds', vmin=0, vmax=1,
                edgecolors='k', linewidths=0.2, alpha=0.85)

cbar = plt.colorbar(sc, ax=ax, shrink=0.7, pad=0.02)
cbar.set_label('Max fractional population loss', fontsize=12)

ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_title('Population crash magnitude', fontsize=14, fontweight='bold')
ax.set_xlim(site_lons.min() - 2, site_lons.max() + 2)
ax.set_ylim(site_lats.min() - 1, site_lats.max() + 1)

fig.tight_layout()
fig.savefig(OUTDIR / 'fig_st_crash_magnitude.png', dpi=DPI, bbox_inches='tight', facecolor='white')
plt.close(fig)
print("  ✓ Saved fig_st_crash_magnitude.png")


# ═══════════════════════════════════════════════════════════════════════
# Fig 5: Final recovery map
# ═══════════════════════════════════════════════════════════════════════
print("Fig 5: Recovery map...")
fig, ax = plt.subplots(figsize=(10, 10))

final_pop = populations[-1, :].astype(float)
recovery_frac = np.where(initial_pop > 0, final_pop / initial_pop, np.nan)

# Clip for display
recovery_clipped = np.clip(recovery_frac, 0, 1.5)

sizes_final = np.clip(final_pop / K * 60, 3, 60)

valid = ~np.isnan(recovery_frac)
sc = ax.scatter(site_lons[valid], site_lats[valid], s=sizes_final[valid],
                c=recovery_clipped[valid], cmap='RdYlGn', vmin=0, vmax=1.2,
                edgecolors='k', linewidths=0.2, alpha=0.85)

cbar = plt.colorbar(sc, ax=ax, shrink=0.7, pad=0.02)
cbar.set_label('Recovery fraction (final / initial pop)', fontsize=12)

ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_title('Final recovery state (2025)', fontsize=14, fontweight='bold')
ax.set_xlim(site_lons.min() - 2, site_lons.max() + 2)
ax.set_ylim(site_lats.min() - 1, site_lats.max() + 1)

fig.tight_layout()
fig.savefig(OUTDIR / 'fig_st_recovery_map.png', dpi=DPI, bbox_inches='tight', facecolor='white')
plt.close(fig)
print("  ✓ Saved fig_st_recovery_map.png")


# ═══════════════════════════════════════════════════════════════════════
# Fig 6: Per-region infection prevalence time series
# ═══════════════════════════════════════════════════════════════════════
print("Fig 6: Prevalence curves...")
fig, ax = plt.subplots(figsize=(12, 7))

for r in scored_regions:
    idx = region_sites[r]
    if len(idx) == 0:
        continue
    total_pop = populations[:, idx].sum(axis=1).astype(float)
    total_inf = infected[:, idx].sum(axis=1).astype(float)
    prevalence = np.where(total_pop > 0, total_inf / total_pop, np.nan)
    ax.plot(years, prevalence, label=r, color=region_colors[r], linewidth=2)

ax.set_xlabel('Year')
ax.set_ylabel('Infection prevalence (infected / population)')
ax.set_title('Epidemic wave timing by region (scored regions)', fontsize=14, fontweight='bold')
ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), fontsize=10, frameon=True)
ax.set_xlim(years[0], years[-1])
ax.set_ylim(0, None)
ax.axhline(0, color='gray', linewidth=0.5, linestyle='--')

fig.tight_layout()
fig.savefig(OUTDIR / 'fig_st_prevalence_curves.png', dpi=DPI, bbox_inches='tight', facecolor='white')
plt.close(fig)
print("  ✓ Saved fig_st_prevalence_curves.png")


# ═══════════════════════════════════════════════════════════════════════
# Fig 7: North-south population cascade
# ═══════════════════════════════════════════════════════════════════════
print("Fig 7: Cascade panels...")
cascade_groups = [
    ('Alaska (PWS + FN + FS)', ['AK-PWS', 'AK-FN', 'AK-FS']),
    ('British Columbia (N + C)', ['BC-N', 'BC-C']),
    ('Salish Sea / WA', ['SS-N', 'SS-S', 'JDF', 'WA-O']),
    ('Oregon', ['OR']),
    ('California (N + C + S)', ['CA-N', 'CA-C', 'CA-S']),
]

fig, axes = plt.subplots(len(cascade_groups), 1, figsize=(12, 12), sharex=True)

for panel_idx, (title, group_regions) in enumerate(cascade_groups):
    ax = axes[panel_idx]
    
    # Combine all sites in the group
    all_idx = np.concatenate([region_sites[r] for r in group_regions if r in region_sites])
    if len(all_idx) == 0:
        ax.set_visible(False)
        continue
    
    total_pop = populations[:, all_idx].sum(axis=1).astype(float)
    initial = total_pop[0]
    relative = total_pop / initial if initial > 0 else total_pop * 0
    
    # Use average color of the group
    group_color_idx = np.mean([region_order.index(r) for r in group_regions])
    color = cm.coolwarm(group_color_idx / (len(region_order) - 1))
    
    ax.fill_between(years, 0, relative, alpha=0.3, color=color)
    ax.plot(years, relative, color=color, linewidth=2)
    ax.set_ylabel('Rel. pop.')
    ax.set_ylim(0, 1.15)
    ax.set_title(title, fontsize=12, fontweight='bold', loc='left')
    ax.axhline(1.0, color='gray', linewidth=0.5, linestyle='--', alpha=0.5)
    ax.axhline(0.0, color='gray', linewidth=0.5, linestyle='--', alpha=0.5)

axes[-1].set_xlabel('Year')
axes[-1].set_xlim(years[0], years[-1])

fig.suptitle('North-to-south population cascade', fontsize=15, fontweight='bold', y=0.99)
fig.tight_layout(rect=[0, 0, 1, 0.97])
fig.savefig(OUTDIR / 'fig_st_cascade_panels.png', dpi=DPI, bbox_inches='tight', facecolor='white')
plt.close(fig)
print("  ✓ Saved fig_st_cascade_panels.png")


# ═══════════════════════════════════════════════════════════════════════
# Fig 8: Population × infection overlay
# ═══════════════════════════════════════════════════════════════════════
print("Fig 8: Pop-infection overlay...")
fig, ax1 = plt.subplots(figsize=(12, 7))

total_pop = populations.sum(axis=1).astype(float)
total_inf = infected.sum(axis=1).astype(float)

color_pop = '#2166ac'
color_inf = '#b2182b'

ax1.plot(years, total_pop / 1e6, color=color_pop, linewidth=2.5, label='Total population')
ax1.set_xlabel('Year')
ax1.set_ylabel('Total population (millions)', color=color_pop, fontsize=13)
ax1.tick_params(axis='y', labelcolor=color_pop)
ax1.set_xlim(years[0], years[-1])
ax1.set_ylim(0, None)

ax2 = ax1.twinx()
ax2.fill_between(years, 0, total_inf / 1e3, alpha=0.25, color=color_inf)
ax2.plot(years, total_inf / 1e3, color=color_inf, linewidth=2, label='Total infected')
ax2.set_ylabel('Total infected (thousands)', color=color_inf, fontsize=13)
ax2.tick_params(axis='y', labelcolor=color_inf)
ax2.set_ylim(0, None)

# Combined legend
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right', fontsize=11, frameon=True)

ax1.set_title('Population and infection dynamics (W257, all sites)', fontsize=14, fontweight='bold')

fig.tight_layout()
fig.savefig(OUTDIR / 'fig_st_pop_infection_overlay.png', dpi=DPI, bbox_inches='tight', facecolor='white')
plt.close(fig)
print("  ✓ Saved fig_st_pop_infection_overlay.png")


# ═══════════════════════════════════════════════════════════════════════
# Fig 9: Seasonal dynamics (AK-PWS)
# ═══════════════════════════════════════════════════════════════════════
print("Fig 9: Seasonal dynamics...")
fig, ax1 = plt.subplots(figsize=(14, 7))

pws_idx = region_sites.get('AK-PWS', np.array([], dtype=int))
print(f"  AK-PWS sites: {len(pws_idx)}")

pws_pop = populations[:, pws_idx].sum(axis=1).astype(float)
pws_inf = infected[:, pws_idx].sum(axis=1).astype(float)

color_pop = '#2166ac'
color_inf = '#b2182b'

ax1.plot(years, pws_pop, color=color_pop, linewidth=1.8, label='Population')
ax1.set_xlabel('Year')
ax1.set_ylabel('Total population (AK-PWS)', color=color_pop, fontsize=13)
ax1.tick_params(axis='y', labelcolor=color_pop)
ax1.set_xlim(years[0], years[-1])

ax2 = ax1.twinx()
ax2.plot(years, pws_inf, color=color_inf, linewidth=1.5, label='Infected')
ax2.set_ylabel('Number infected (AK-PWS)', color=color_inf, fontsize=13)
ax2.tick_params(axis='y', labelcolor=color_inf)

# Spawning season shading (Nov-Mar of each year)
for yr in range(sst_start_year, sst_start_year + 14):
    # Nov of year yr to Mar of year yr+1
    nov_start = yr + 10/12  # November
    mar_end = yr + 1 + 3/12  # March
    ax1.axvspan(nov_start, mar_end, alpha=0.08, color='blue', zorder=0)

# Annotate one spawning band
ax1.annotate('Spawning\nseason', xy=(2015.9, pws_pop.max() * 0.9),
             fontsize=9, color='blue', alpha=0.6, ha='center')

lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right', fontsize=11, frameon=True)

ax1.set_title('Seasonal dynamics at monthly resolution (AK-PWS)', fontsize=14, fontweight='bold')

fig.tight_layout()
fig.savefig(OUTDIR / 'fig_st_seasonal_dynamics.png', dpi=DPI, bbox_inches='tight', facecolor='white')
plt.close(fig)
print("  ✓ Saved fig_st_seasonal_dynamics.png")


# ═══════════════════════════════════════════════════════════════════════
# Fig 10: Recovery vs latitude scatter
# ═══════════════════════════════════════════════════════════════════════
print("Fig 10: Recovery vs latitude...")
fig, ax = plt.subplots(figsize=(12, 8))

recovery = np.where(initial_pop > 0, final_pop / initial_pop, np.nan)

# Plot each region separately for legend
for r in region_order:
    idx = region_sites[r]
    if len(idx) == 0:
        continue
    valid_mask = ~np.isnan(recovery[idx])
    if not valid_mask.any():
        continue
    ax.scatter(site_lats[idx][valid_mask], recovery[idx][valid_mask],
               c=[region_colors[r]], s=20, alpha=0.6, label=r, edgecolors='none')

# Reference lines
ax.axhline(1.0, color='green', linewidth=1, linestyle='--', alpha=0.5, label='Full recovery')
ax.axhline(0.0, color='red', linewidth=1, linestyle='--', alpha=0.5)

# Moving average trend line (binned by latitude)
valid_all = ~np.isnan(recovery)
lat_sorted_idx = np.argsort(site_lats[valid_all])
lats_sorted = site_lats[valid_all][lat_sorted_idx]
rec_sorted = recovery[valid_all][lat_sorted_idx]

# Bin average
n_bins = 30
bin_edges = np.linspace(lats_sorted.min(), lats_sorted.max(), n_bins + 1)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
bin_means = np.zeros(n_bins)
for b in range(n_bins):
    mask = (lats_sorted >= bin_edges[b]) & (lats_sorted < bin_edges[b + 1])
    if mask.any():
        bin_means[b] = np.nanmean(rec_sorted[mask])
    else:
        bin_means[b] = np.nan

valid_bins = ~np.isnan(bin_means)
ax.plot(bin_centers[valid_bins], bin_means[valid_bins], 'k-', linewidth=3, alpha=0.7, label='Trend (binned mean)')

ax.set_xlabel('Latitude (°N)', fontsize=13)
ax.set_ylabel('Recovery fraction (final / initial pop)', fontsize=13)
ax.set_title('Site-level recovery vs latitude', fontsize=14, fontweight='bold')

# Use log scale if range is large
if np.nanmin(recovery[valid_all]) > 0:
    ax.set_yscale('log')
    ax.set_ylabel('Recovery fraction (log scale)', fontsize=13)
else:
    # Add small offset to see zeros - use symlog
    ax.set_yscale('symlog', linthresh=0.01)
    ax.set_ylabel('Recovery fraction (symlog scale)', fontsize=13)

ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), fontsize=9, frameon=True, ncol=1)

fig.tight_layout()
fig.savefig(OUTDIR / 'fig_st_recovery_vs_latitude.png', dpi=DPI, bbox_inches='tight', facecolor='white')
plt.close(fig)
print("  ✓ Saved fig_st_recovery_vs_latitude.png")


# ═══════════════════════════════════════════════════════════════════════
# Verify outputs
# ═══════════════════════════════════════════════════════════════════════
print("\n" + "="*60)
print("Output verification:")
expected = [
    'fig_st_population_heatmap.png',
    'fig_st_wavefront_snapshots.png',
    'fig_st_arrival_timing.png',
    'fig_st_crash_magnitude.png',
    'fig_st_recovery_map.png',
    'fig_st_prevalence_curves.png',
    'fig_st_cascade_panels.png',
    'fig_st_pop_infection_overlay.png',
    'fig_st_seasonal_dynamics.png',
    'fig_st_recovery_vs_latitude.png',
]
all_ok = True
for fname in expected:
    p = OUTDIR / fname
    if p.exists() and p.stat().st_size > 0:
        size_kb = p.stat().st_size / 1024
        print(f"  ✓ {fname:45s} {size_kb:8.1f} KB")
    else:
        print(f"  ✗ {fname:45s} MISSING OR EMPTY")
        all_ok = False

if all_ok:
    print("\nAll 10 figures generated successfully! ✅")
else:
    print("\nSome figures missing! ❌")
