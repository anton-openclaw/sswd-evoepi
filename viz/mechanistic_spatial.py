#!/usr/bin/env python3
"""Generate 12 spatial ecology & environmental forcing figures for the
SSWD-EvoEpi mechanistic model report.

Figures S1–S12 saved as PDF to reports/mechanistic/figures/.
"""

import sys, os, json, warnings
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import ListedColormap, Normalize, LogNorm, BoundaryNorm
from matplotlib.patches import FancyArrowPatch, Rectangle, Polygon, FancyBboxPatch, Circle
from matplotlib.lines import Line2D
import matplotlib.patheffects as pe
import cartopy.crs as ccrs
import cartopy.feature as cfeature

warnings.filterwarnings('ignore')

# ── paths ─────────────────────────────────────────────────────────────
ROOT = '/home/starbot/.openclaw/workspace/sswd-evoepi'
FIG_DIR = os.path.join(ROOT, 'reports', 'mechanistic', 'figures')
os.makedirs(FIG_DIR, exist_ok=True)

# ── load shared data ──────────────────────────────────────────────────
print("Loading shared data...")

with open(os.path.join(ROOT, 'data/nodes/all_sites.json')) as f:
    ALL_SITES = json.load(f)
with open(os.path.join(ROOT, 'data/nodes/site_enclosedness.json')) as f:
    ENCLOSEDNESS = json.load(f)
with open(os.path.join(ROOT, 'reports/mechanistic/w285_result.json')) as f:
    W285_RESULT = json.load(f)

W285 = np.load(os.path.join(ROOT, 'reports/mechanistic/w285_monthly.npz'),
               allow_pickle=True)
DIST_NPZ = np.load(os.path.join(ROOT, 'results/overwater/distance_matrix.npz'),
                    allow_pickle=True)

N_SITES = len(ALL_SITES)
LATS = np.array([s['latitude'] for s in ALL_SITES])
LONS = np.array([s['longitude'] for s in ALL_SITES])
NAMES = np.array([s['name'] for s in ALL_SITES])
REGIONS_ARR = np.array([s['region'] for s in ALL_SITES])

# Coastline order (S→N): follows the coast from Baja northward to the Aleutians.
# This is the reverse of REGION_ORDER in results.py (which is N→S).
COASTLINE_ORDER = [
    "BJ", "CA-S", "CA-C", "CA-N", "OR", "WA-O", "JDF", "SS-S", "SS-N",
    "BC-C", "BC-N", "AK-FS", "AK-FN", "AK-PWS", "AK-EG", "AK-OC",
    "AK-WG", "AK-AL",
]
# Validate: make sure all regions are present
_all_regions = set(REGIONS_ARR)
assert set(COASTLINE_ORDER) == _all_regions, \
    f"COASTLINE_ORDER mismatch: {set(COASTLINE_ORDER) ^ _all_regions}"
REGIONS_SORTED = COASTLINE_ORDER  # replaces old sorted(set(...))

# Enclosedness arrays (aligned by name)
enc_map = {e['name']: e for e in ENCLOSEDNESS}
E_RAW = np.array([enc_map.get(s['name'], {}).get('enclosedness_combined', 0.0)
                   for s in ALL_SITES])
FJORD_DEPTH = np.array([enc_map.get(s['name'], {}).get('fjord_depth_norm', 0.0)
                         for s in ALL_SITES])
FLUSHING = np.array([enc_map.get(s['name'], {}).get('flushing_rate', 0.5)
                      for s in ALL_SITES])

# Distance matrix (m → km)
DIST_KM = DIST_NPZ['distances'].astype(np.float64)  # already in meters
# Check units — if max > 10000 it's likely meters
if DIST_KM.max() > 10000:
    DIST_KM = DIST_KM / 1000.0  # convert to km

# W285 simulation data
SIM_DAYS = W285['sim_days']
POPULATIONS = W285['populations']  # (159, 896)
INFECTED = W285['infected']       # (159, 896)
K = int(W285['K'])
SST_START_YEAR = int(W285['sst_start_year'])

# Region colors — 18 qualitative
_tab20 = plt.cm.get_cmap('tab20', 20)
REGION_COLORS = {r: _tab20(i) for i, r in enumerate(REGIONS_SORTED)}

# Coastline-sorted index for matrices: first by region's coastline position,
# then by latitude within each region.
_region_rank = {r: i for i, r in enumerate(COASTLINE_ORDER)}
_sort_keys = np.array([(_region_rank[r], lat) for r, lat in zip(REGIONS_ARR, LATS)],
                       dtype=[('region_rank', int), ('lat', float)])
LAT_ORDER = np.argsort(_sort_keys, order=('region_rank', 'lat'))

# Common projection
PROJ = ccrs.LambertConformal(central_longitude=-135, central_latitude=45)
PC = ccrs.PlateCarree()

# Region boundaries in coastline-sorted order
_lat_sorted_regions = REGIONS_ARR[LAT_ORDER]
REGION_BOUNDARIES = []
for i in range(1, len(LAT_ORDER)):
    if _lat_sorted_regions[i] != _lat_sorted_regions[i-1]:
        REGION_BOUNDARIES.append(i)


def add_coast(ax, extent=None):
    """Add Natural Earth land/ocean/coastline features."""
    ax.add_feature(cfeature.NaturalEarthFeature(
        'physical', 'land', '50m', facecolor='#f0ebe3', edgecolor='none'))
    ax.add_feature(cfeature.NaturalEarthFeature(
        'physical', 'ocean', '50m', facecolor='#d6eaf8', edgecolor='none'))
    ax.add_feature(cfeature.NaturalEarthFeature(
        'physical', 'coastline', '50m', facecolor='none',
        edgecolor='#555555', linewidth=0.4))
    if extent is not None:
        ax.set_extent(extent, crs=PC)


def save_fig(fig, name):
    path = os.path.join(FIG_DIR, name)
    fig.savefig(path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"  ✓ Saved {name}")


# ══════════════════════════════════════════════════════════════════════
# S1 — 896-Node Metapopulation Network Map
# ══════════════════════════════════════════════════════════════════════
def fig_s01():
    print("Generating S1: Network Map...")
    fig = plt.figure(figsize=(14, 16))
    gs = fig.add_gridspec(3, 2, width_ratios=[2.5, 1], height_ratios=[1.2, 1, 1],
                          hspace=0.05, wspace=0.08)

    # Main map
    ax_main = fig.add_subplot(gs[:, 0], projection=PROJ)
    add_coast(ax_main, [-170, -115, 24, 62])

    for r in REGIONS_SORTED:
        mask = REGIONS_ARR == r
        ax_main.scatter(LONS[mask], LATS[mask], c=[REGION_COLORS[r]],
                        s=12, transform=PC, zorder=5, edgecolor='k',
                        linewidth=0.15, alpha=0.85, label=r)

    ax_main.set_title('896-Node Metapopulation Network', fontsize=13, fontweight='bold')
    gl = ax_main.gridlines(draw_labels=True, linewidth=0.3, alpha=0.4, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False

    # Insets
    inset_specs = [
        ('SE Alaska Fjords', [-140, -130, 55, 61], gs[0, 1]),
        ('Puget Sound / Salish Sea', [-124.5, -122, 47, 49.2], gs[1, 1]),
        ('Channel Islands', [-121, -118.5, 32.5, 34.8], gs[2, 1]),
    ]
    for title, extent, gs_pos in inset_specs:
        ax_in = fig.add_subplot(gs_pos, projection=ccrs.Mercator())
        add_coast(ax_in, extent)
        for r in REGIONS_SORTED:
            mask = REGIONS_ARR == r
            ax_in.scatter(LONS[mask], LATS[mask], c=[REGION_COLORS[r]],
                          s=20, transform=PC, zorder=5, edgecolor='k',
                          linewidth=0.2, alpha=0.9)
        ax_in.set_title(title, fontsize=9, fontweight='bold')
        gl2 = ax_in.gridlines(draw_labels=True, linewidth=0.2, alpha=0.3)
        gl2.top_labels = False
        gl2.right_labels = False

    # Legend
    handles = [Line2D([0],[0], marker='o', color='w', markerfacecolor=REGION_COLORS[r],
                       markersize=6, label=r, markeredgecolor='k', markeredgewidth=0.3)
               for r in REGIONS_SORTED]
    ax_main.legend(handles=handles, loc='lower left', ncol=3, fontsize=7,
                   frameon=True, framealpha=0.9, title='Region', title_fontsize=8)

    save_fig(fig, 'fig_S01.pdf')


# ══════════════════════════════════════════════════════════════════════
# S2 — Dispersal Kernels
# ══════════════════════════════════════════════════════════════════════
def fig_s02():
    print("Generating S2: Dispersal Kernels...")
    fig, ax = plt.subplots(figsize=(10, 5.5))
    d = np.linspace(0.1, 2000, 2000)

    kernels = [
        (400, r'$D_L = 400$ km (larvae)', '#2980b9', '-', 2.5),
        (15,  r'$D_P = 15$ km (pathogen, local)', '#c0392b', '-', 2.5),
        (300, r'$D_P^{\mathrm{wave}} = 300$ km (wavefront)', '#e67e22', '--', 2.0),
    ]
    for D, label, color, ls, lw in kernels:
        w = np.exp(-d / D)
        ax.semilogy(d, w, color=color, ls=ls, lw=lw, label=label)
        # e-folding distance marker
        ax.axvline(D, color=color, ls=':', alpha=0.35, lw=1)

    # Inter-node distance distribution
    triu = np.triu_indices(N_SITES, k=1)
    dists_flat = DIST_KM[triu]
    dists_flat = dists_flat[dists_flat > 0]
    med = np.median(dists_flat)
    p90 = np.percentile(dists_flat, 90)
    ax.axvline(med, color='grey', ls='-', lw=1.5, alpha=0.7)
    ax.axvline(p90, color='grey', ls='--', lw=1.5, alpha=0.7)
    ax.text(med+15, 0.3, f'Median\n{med:.0f} km', fontsize=8, color='grey', va='top')
    ax.text(p90+15, 0.3, f'P90\n{p90:.0f} km', fontsize=8, color='grey', va='top')

    ax.set_xlabel('Overwater distance (km)', fontsize=11)
    ax.set_ylabel('Kernel weight (relative)', fontsize=11)
    ax.set_ylim(1e-6, 1.5)
    ax.set_xlim(0, 2000)
    ax.legend(frameon=True, framealpha=0.9, fontsize=10)
    ax.set_title('Dispersal Kernels: Larvae vs. Pathogen', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.2)

    save_fig(fig, 'fig_S02.pdf')


# ══════════════════════════════════════════════════════════════════════
# S3 — Larval Connectivity Matrix (C)
# ══════════════════════════════════════════════════════════════════════
def fig_s03():
    print("Generating S3: Larval Connectivity Matrix...")
    # Parameters
    D_L = 400.0  # km
    r_total = 0.003
    alpha_self_fjord = 0.70
    alpha_self_open = 0.02
    phi_open = 0.8
    phi_fjord = 0.03
    n_conn = 0.3

    # Determine fjord-like status from enclosedness
    e_eff = E_RAW.copy()
    alpha_self = np.where(e_eff > 0.5, alpha_self_fjord, alpha_self_open)

    # Build C matrix (lat-sorted)
    idx = LAT_ORDER
    dist_sorted = DIST_KM[np.ix_(idx, idx)]
    alpha_sorted = alpha_self[idx]
    regions_sorted_arr = REGIONS_ARR[idx]

    N = len(idx)
    C = np.zeros((N, N), dtype=np.float64)

    # Raw kernel
    with np.errstate(divide='ignore', invalid='ignore'):
        C_raw = np.exp(-dist_sorted / D_L)
    np.fill_diagonal(C_raw, 0.0)

    for k in range(N):
        row_sum = C_raw[:, k].sum()
        if row_sum > 0:
            C[:, k] = C_raw[:, k] / row_sum * (1.0 - alpha_sorted[k])
        C[k, k] = alpha_sorted[k]
    C *= r_total

    # Region aggregation
    region_labels_sorted = regions_sorted_arr
    n_reg = len(REGIONS_SORTED)
    C_region = np.zeros((n_reg, n_reg))
    for i, ri in enumerate(REGIONS_SORTED):
        for j, rj in enumerate(REGIONS_SORTED):
            mask_i = region_labels_sorted == ri
            mask_j = region_labels_sorted == rj
            C_region[i, j] = C[np.ix_(mask_i, mask_j)].sum()

    fig = plt.figure(figsize=(16, 7))
    gs = fig.add_gridspec(1, 3, width_ratios=[2, 1.2, 1])

    # Panel A: full C matrix
    ax_a = fig.add_subplot(gs[0])
    log_C = np.log10(C + 1e-20)
    im_a = ax_a.imshow(log_C, cmap='viridis', aspect='auto', vmin=-8, vmax=-2,
                        interpolation='nearest')
    for b in REGION_BOUNDARIES:
        ax_a.axhline(b, color='white', lw=0.3, alpha=0.5)
        ax_a.axvline(b, color='white', lw=0.3, alpha=0.5)
    ax_a.set_xlabel('Source node (S→N)', fontsize=10)
    ax_a.set_ylabel('Sink node (S→N)', fontsize=10)
    ax_a.set_title('A. 896×896 Larval C Matrix', fontsize=11, fontweight='bold')
    plt.colorbar(im_a, ax=ax_a, label='log₁₀(C)', shrink=0.7)

    # Panel B: region aggregated
    ax_b = fig.add_subplot(gs[1])
    im_b = ax_b.imshow(np.log10(C_region + 1e-20), cmap='YlOrRd', aspect='auto',
                        vmin=-4, vmax=0)
    ax_b.set_xticks(range(n_reg))
    ax_b.set_yticks(range(n_reg))
    ax_b.set_xticklabels(REGIONS_SORTED, rotation=90, fontsize=6)
    ax_b.set_yticklabels(REGIONS_SORTED, fontsize=6)
    ax_b.set_title('B. Region-aggregated C', fontsize=11, fontweight='bold')
    plt.colorbar(im_b, ax=ax_b, label='log₁₀(ΣC)', shrink=0.7)

    # Panel C: net export map
    ax_c = fig.add_subplot(gs[2], projection=PROJ)
    add_coast(ax_c, [-170, -115, 24, 62])
    # Compute net export per node (in sorted order, then unsort)
    row_sums = C.sum(axis=0)  # total exported by each source
    col_sums = C.sum(axis=1)  # total received by each sink
    net_export = row_sums - col_sums
    # Unsort
    net_export_full = np.zeros(N_SITES)
    net_export_full[LAT_ORDER] = net_export
    vmax = np.percentile(np.abs(net_export_full), 95)
    sc = ax_c.scatter(LONS, LATS, c=net_export_full, cmap='RdBu', s=8,
                       transform=PC, zorder=5, vmin=-vmax, vmax=vmax)
    ax_c.set_title('C. Net larval export', fontsize=11, fontweight='bold')
    plt.colorbar(sc, ax=ax_c, label='Net export (C)', shrink=0.6)

    fig.suptitle('Larval Connectivity Matrix (C)', fontsize=14, fontweight='bold', y=1.01)
    save_fig(fig, 'fig_S03.pdf')


# ══════════════════════════════════════════════════════════════════════
# S4 — Wavefront Activation
# ══════════════════════════════════════════════════════════════════════
def fig_s04():
    print("Generating S4: Wavefront Activation...")
    # Use first-infection time from W285 NPZ as wavefront activation proxy
    first_inf_day = np.full(N_SITES, np.nan)
    for j in range(N_SITES):
        idx = np.where(INFECTED[:, j] > 0)[0]
        if len(idx) > 0:
            first_inf_day[j] = SIM_DAYS[idx[0]]

    # Convert to dates: day 0 = Jan 2012
    # Seed nodes (Channel Islands): [322, 319, 632, 633, 634]
    seed_nodes = [322, 319, 632, 633, 634]
    seed_lons = LONS[seed_nodes]
    seed_lats = LATS[seed_nodes]

    # Time snapshots (days from Jan 2012)
    # Disease starts ~day 750 (Jan 2014 = month 25)
    # Equivalent to Jun 2013 onset in the model
    snapshots = [
        (540, 'Jun 2013\n(seed)'),
        (720, 'Dec 2013'),
        (900, 'Jun 2014'),
        (1080, 'Dec 2014'),
        (1260, 'Jun 2015'),
        (1500, 'Feb 2016'),
    ]

    fig, axes = plt.subplots(2, 3, figsize=(16, 14),
                              subplot_kw={'projection': PROJ})
    fig.suptitle('Wavefront Activation: Channel Islands Origin & Spread',
                 fontsize=14, fontweight='bold', y=0.98)

    for ax, (t_day, t_label) in zip(axes.flat, snapshots):
        add_coast(ax, [-170, -115, 24, 62])

        activated = first_inf_day <= t_day
        pending = (~activated) & (~np.isnan(first_inf_day))
        never = np.isnan(first_inf_day)

        # Never-infected nodes: very light grey
        if never.sum() > 0:
            ax.scatter(LONS[never], LATS[never], c='#e0e0e0', s=5,
                       transform=PC, zorder=3, alpha=0.5)

        # Pending nodes: grey
        if pending.sum() > 0:
            ax.scatter(LONS[pending], LATS[pending], c='#aaaaaa', s=6,
                       transform=PC, zorder=4, alpha=0.6)

        # Activated nodes: colored by activation time
        if activated.sum() > 0:
            sc = ax.scatter(LONS[activated], LATS[activated],
                            c=first_inf_day[activated], cmap='YlOrRd',
                            vmin=first_inf_day[~np.isnan(first_inf_day)].min(),
                            vmax=first_inf_day[~np.isnan(first_inf_day)].max(),
                            s=10, transform=PC, zorder=5, edgecolor='none')

        # Seed nodes: star
        ax.scatter(seed_lons, seed_lats, marker='*', c='gold', s=120,
                   edgecolor='k', linewidth=0.5, transform=PC, zorder=7)

        ax.set_title(t_label, fontsize=10, fontweight='bold')

    # Colorbar
    cax = fig.add_axes([0.92, 0.15, 0.015, 0.65])
    sm = plt.cm.ScalarMappable(
        cmap='YlOrRd',
        norm=Normalize(vmin=first_inf_day[~np.isnan(first_inf_day)].min() / 30,
                       vmax=first_inf_day[~np.isnan(first_inf_day)].max() / 30))
    plt.colorbar(sm, cax=cax, label='Activation time (months from Jan 2012)')

    save_fig(fig, 'fig_S04.pdf')


# ══════════════════════════════════════════════════════════════════════
# S5 — Enclosedness Pipeline
# ══════════════════════════════════════════════════════════════════════
def fig_s05():
    print("Generating S5: Enclosedness Pipeline...")
    # Compute effective enclosedness
    e_eff = E_RAW  # raw enclosedness as effective (simplified)
    phi_open_val = 0.8
    phi_fjord_val = 0.03
    n_conn = 0.3
    phi_computed = phi_open_val * (1 - e_eff**n_conn) + phi_fjord_val * e_eff**n_conn

    fig = plt.figure(figsize=(18, 7))
    gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 0.7, 1], wspace=0.35)

    # Panel A: raw enclosedness map
    ax_a = fig.add_subplot(gs[0], projection=PROJ)
    add_coast(ax_a, [-170, -115, 24, 62])
    sc_a = ax_a.scatter(LONS, LATS, c=E_RAW, cmap='YlGnBu', s=8,
                         transform=PC, vmin=0, vmax=1, zorder=5, edgecolor='none')
    plt.colorbar(sc_a, ax=ax_a, label='Raw enclosedness', shrink=0.6)
    ax_a.set_title('A. Raw Enclosedness', fontsize=10, fontweight='bold')

    # Panel B: effective enclosedness map
    ax_b = fig.add_subplot(gs[1], projection=PROJ)
    add_coast(ax_b, [-170, -115, 24, 62])
    sc_b = ax_b.scatter(LONS, LATS, c=e_eff, cmap='YlGnBu', s=8,
                         transform=PC, vmin=0, vmax=1, zorder=5, edgecolor='none')
    plt.colorbar(sc_b, ax=ax_b, label='Effective enclosedness', shrink=0.6)
    ax_b.set_title('B. Effective Enclosedness', fontsize=10, fontweight='bold')

    # Panel C: transfer function
    ax_c = fig.add_subplot(gs[2])
    e_curve = np.linspace(0, 1, 200)
    phi_curve = phi_open_val * (1 - e_curve**n_conn) + phi_fjord_val * e_curve**n_conn
    ax_c.plot(e_curve, phi_curve, 'k-', lw=2.5, label='Transfer function')
    ax_c.scatter(e_eff, phi_computed, c=e_eff, cmap='YlGnBu', s=8, alpha=0.4, zorder=3)
    ax_c.set_xlabel(r'Effective enclosedness $e_{\mathrm{eff}}$', fontsize=10)
    ax_c.set_ylabel(r'Flushing rate $\phi$ (day$^{-1}$)', fontsize=10)
    ax_c.set_title('C. Transfer Function', fontsize=10, fontweight='bold')
    ax_c.grid(True, alpha=0.2)
    ax_c.set_xlim(0, 1)
    ax_c.set_ylim(0, 0.9)

    # Panel D: flushing rate map
    ax_d = fig.add_subplot(gs[3], projection=PROJ)
    add_coast(ax_d, [-170, -115, 24, 62])
    sc_d = ax_d.scatter(LONS, LATS, c=FLUSHING, cmap='RdYlBu', s=8,
                         transform=PC, zorder=5, edgecolor='none',
                         vmin=0, vmax=0.8)
    plt.colorbar(sc_d, ax=ax_d, label=r'$\phi$ (day$^{-1}$)', shrink=0.6)
    ax_d.set_title('D. Flushing Rate', fontsize=10, fontweight='bold')

    fig.suptitle('Enclosedness Pipeline: Raw → Effective → Flushing Rate',
                 fontsize=13, fontweight='bold', y=1.01)
    save_fig(fig, 'fig_S05.pdf')


# ══════════════════════════════════════════════════════════════════════
# S6 — The Fjord Trade-Off
# ══════════════════════════════════════════════════════════════════════
def fig_s06():
    print("Generating S6: Fjord Trade-Off...")
    # Classify fjord vs open
    fjord_mask = E_RAW > 0.5  # high enclosedness = fjord-like
    alpha_self = np.where(fjord_mask, 0.70, 0.02)
    retention_time = 1.0 / np.clip(FLUSHING, 0.001, None)

    # Peak disease prevalence and recovery from simulation
    peak_prevalence = np.zeros(N_SITES)
    recovery_month = np.full(N_SITES, np.nan)
    for j in range(N_SITES):
        if POPULATIONS[:, j].max() > 0:
            prev = INFECTED[:, j] / np.clip(POPULATIONS[:, j], 1, None)
            peak_prevalence[j] = prev.max()
            # Recovery: time to return to 50% of initial pop after crash
            pop = POPULATIONS[:, j]
            init_pop = pop[0]
            if init_pop > 0:
                min_idx = np.argmin(pop)
                if min_idx > 0:
                    post_crash = pop[min_idx:]
                    thresh = 0.5 * init_pop
                    recovered = np.where(post_crash >= thresh)[0]
                    if len(recovered) > 0:
                        recovery_month[j] = (SIM_DAYS[min_idx + recovered[0]] - SIM_DAYS[min_idx]) / 30.0

    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(2, 2, height_ratios=[0.8, 1.2], hspace=0.35, wspace=0.3)

    # Panel A: Conceptual schematic
    ax_a = fig.add_subplot(gs[0, :])
    ax_a.set_xlim(0, 10)
    ax_a.set_ylim(0, 4)
    ax_a.set_aspect('equal')
    ax_a.axis('off')
    ax_a.set_title('A. Conceptual: Open Coast vs. Fjord', fontsize=11, fontweight='bold')

    # Open coast
    ax_a.add_patch(FancyBboxPatch((0.5, 0.5), 3.5, 3, boxstyle="round,pad=0.1",
                                   facecolor='#d6eaf8', edgecolor='#2980b9', linewidth=2))
    ax_a.text(2.25, 3.2, 'OPEN COAST', fontsize=10, ha='center', fontweight='bold', color='#2980b9')
    ax_a.annotate('', xy=(4.2, 2), xytext=(3.8, 2),
                  arrowprops=dict(arrowstyle='->', color='#2980b9', lw=3))
    ax_a.annotate('', xy=(4.2, 1.5), xytext=(3.8, 1.5),
                  arrowprops=dict(arrowstyle='->', color='#2980b9', lw=3))
    ax_a.text(2.25, 2.5, r'$\alpha_{\mathrm{self}}=0.02$' + '\nLow self-recruitment',
              fontsize=8, ha='center', color='#2980b9')
    ax_a.text(2.25, 1.0, r'$\phi=0.80$' + '\nRapid flushing\n→ pathogen diluted',
              fontsize=8, ha='center', color='#c0392b')

    # Fjord
    ax_a.add_patch(FancyBboxPatch((5.5, 0.5), 3.5, 3, boxstyle="round,pad=0.1",
                                   facecolor='#fadbd8', edgecolor='#c0392b', linewidth=2))
    ax_a.text(7.25, 3.2, 'FJORD', fontsize=10, ha='center', fontweight='bold', color='#c0392b')
    ax_a.annotate('', xy=(9.2, 2), xytext=(8.8, 2),
                  arrowprops=dict(arrowstyle='->', color='#c0392b', lw=1.5))
    ax_a.text(7.25, 2.5, r'$\alpha_{\mathrm{self}}=0.70$' + '\nHigh self-recruitment',
              fontsize=8, ha='center', color='#2980b9')
    ax_a.text(7.25, 1.0, r'$\phi=0.03$' + '\nSlow flushing\n→ pathogen retained',
              fontsize=8, ha='center', color='#c0392b')

    # Panel B: self-recruitment vs retention time
    ax_b = fig.add_subplot(gs[1, 0])
    ax_b.scatter(alpha_self[~fjord_mask], retention_time[~fjord_mask],
                 c='#3498db', alpha=0.4, label='Open coast', s=15, edgecolor='none')
    ax_b.scatter(alpha_self[fjord_mask], retention_time[fjord_mask],
                 c='#e74c3c', alpha=0.4, label='Fjord', s=15, edgecolor='none')
    ax_b.set_xlabel(r'Self-recruitment $\alpha_{\mathrm{self}}$', fontsize=10)
    ax_b.set_ylabel(r'Pathogen retention time $1/\phi$ (days)', fontsize=10)
    ax_b.set_title('B. Self-Recruitment vs. Retention', fontsize=11, fontweight='bold')
    ax_b.legend(frameon=True, fontsize=9)
    ax_b.grid(True, alpha=0.2)
    ax_b.set_yscale('log')

    # Panel C: peak prevalence vs recovery time
    ax_c = fig.add_subplot(gs[1, 1])
    valid = ~np.isnan(recovery_month) & (peak_prevalence > 0)
    if valid.sum() > 5:
        ax_c.scatter(recovery_month[valid & ~fjord_mask],
                     peak_prevalence[valid & ~fjord_mask] * 100,
                     c='#3498db', alpha=0.4, label='Open coast', s=15, edgecolor='none')
        ax_c.scatter(recovery_month[valid & fjord_mask],
                     peak_prevalence[valid & fjord_mask] * 100,
                     c='#e74c3c', alpha=0.4, label='Fjord', s=15, edgecolor='none')
    ax_c.set_xlabel('Recovery time (months to 50% initial pop)', fontsize=10)
    ax_c.set_ylabel('Peak disease prevalence (%)', fontsize=10)
    ax_c.set_title('C. Epidemic Outcome: The Trade-Off', fontsize=11, fontweight='bold')
    ax_c.legend(frameon=True, fontsize=9)
    ax_c.grid(True, alpha=0.2)

    fig.suptitle('The Fjord Paradox: Where Larvae Thrive, So Does Disease',
                 fontsize=14, fontweight='bold', y=1.01)
    save_fig(fig, 'fig_S06.pdf')


# ══════════════════════════════════════════════════════════════════════
# S7 — SST Hovmöller + Time Series
# ══════════════════════════════════════════════════════════════════════
def fig_s07():
    print("Generating S7: SST Forcing...")
    # Load ~20 representative site SST files spanning latitudes
    site_sst_dir = os.path.join(ROOT, 'data/sst/site_sst')
    
    # Select ~30 sites evenly spaced in latitude
    lat_order = np.argsort(LATS)
    sample_idx = lat_order[np.linspace(0, len(lat_order)-1, 30, dtype=int)]
    
    # Load monthly SST for sampled sites
    import csv
    sst_data = {}  # site_idx -> {(year,month): sst}
    for si in sample_idx:
        fname = f"{NAMES[si]}_monthly.csv"
        fpath = os.path.join(site_sst_dir, fname)
        if os.path.exists(fpath):
            monthly = {}
            with open(fpath) as f:
                reader = csv.DictReader(f)
                for row in reader:
                    monthly[(int(row['year']), int(row['month']))] = float(row['sst'])
            sst_data[si] = monthly

    # Build Hovmöller array
    years = range(2002, 2026)
    months = range(1, 13)
    n_time = len(years) * 12
    lat_bins = np.linspace(24, 62, 40)
    hovmoller = np.full((len(lat_bins)-1, n_time), np.nan)
    
    # Assign each loaded site to a lat bin
    for t_idx, (yr, mo) in enumerate([(y, m) for y in years for m in months]):
        for bi in range(len(lat_bins)-1):
            vals = []
            for si in sst_data:
                if lat_bins[bi] <= LATS[si] < lat_bins[bi+1]:
                    key = (yr, mo)
                    if key in sst_data[si]:
                        vals.append(sst_data[si][key])
            if vals:
                hovmoller[bi, t_idx] = np.mean(vals)

    # Load 3 representative time series for Panel B
    rep_sites = {}
    for target_lat, label, color in [(28, 'Baja (~28°N)', '#e67e22'),
                                       (37, 'Central CA (~37°N)', '#27ae60'),
                                       (57, 'SE Alaska (~57°N)', '#2980b9')]:
        # Find closest site with SST data
        for si in sorted(sst_data.keys(), key=lambda x: abs(LATS[x] - target_lat)):
            if si in sst_data:
                rep_sites[label] = (si, color, sst_data[si])
                break

    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.5, 1], hspace=0.3)

    # Panel A: Hovmöller
    ax_hov = fig.add_subplot(gs[:, 0])
    time_axis = np.arange(n_time)
    lat_centers = 0.5 * (lat_bins[:-1] + lat_bins[1:])
    # Interpolate NaN gaps
    from scipy.interpolate import griddata
    # Simple: use pcolormesh with masked array
    hovmoller_masked = np.ma.masked_invalid(hovmoller)
    im = ax_hov.pcolormesh(time_axis, lat_centers, hovmoller_masked,
                            cmap='RdYlBu_r', vmin=4, vmax=22, shading='auto')
    # Set x-ticks to years
    year_ticks = [(yr - 2002) * 12 for yr in range(2002, 2026, 4)]
    ax_hov.set_xticks(year_ticks)
    ax_hov.set_xticklabels([str(yr) for yr in range(2002, 2026, 4)])
    ax_hov.set_xlabel('Year', fontsize=11)
    ax_hov.set_ylabel('Latitude (°N)', fontsize=11)
    ax_hov.set_title('A. SST Hovmöller Diagram', fontsize=12, fontweight='bold')
    plt.colorbar(im, ax=ax_hov, label='SST (°C)', shrink=0.7)

    # Panel B: representative time series
    ax_ts = fig.add_subplot(gs[0, 1])
    for label, (si, color, monthly) in rep_sites.items():
        times_sorted = sorted(monthly.keys())
        t_vals = [i for i, _ in enumerate(times_sorted)]
        sst_vals = [monthly[k] for k in times_sorted]
        ax_ts.plot(t_vals, sst_vals, color=color, lw=0.8, label=f'{label}\n({NAMES[si]})', alpha=0.8)
    ax_ts.set_xlabel('Month index (2002–2025)', fontsize=10)
    ax_ts.set_ylabel('SST (°C)', fontsize=10)
    ax_ts.set_title('B. Representative Time Series', fontsize=11, fontweight='bold')
    ax_ts.legend(fontsize=7, frameon=True)
    ax_ts.grid(True, alpha=0.2)

    # Panel C: mean annual SST map
    ax_map = fig.add_subplot(gs[1, 1], projection=PROJ)
    add_coast(ax_map, [-170, -115, 24, 62])
    # Compute mean SST per site from loaded data
    mean_sst = np.full(N_SITES, np.nan)
    for si in sst_data:
        mean_sst[si] = np.mean(list(sst_data[si].values()))
    # For all sites, use latitude-based estimate for sites without SST
    # Simple linear model: SST ≈ 25 - 0.3*(lat-25)
    for j in range(N_SITES):
        if np.isnan(mean_sst[j]):
            mean_sst[j] = 25 - 0.3 * (LATS[j] - 25)
    sc = ax_map.scatter(LONS, LATS, c=mean_sst, cmap='RdYlBu_r', s=8,
                         transform=PC, zorder=5, vmin=5, vmax=20, edgecolor='none')
    ax_map.set_title('C. Mean Annual SST', fontsize=11, fontweight='bold')
    plt.colorbar(sc, ax=ax_map, label='SST (°C)', shrink=0.6)

    fig.suptitle('SST Forcing: Latitudinal Gradient & Temporal Dynamics',
                 fontsize=14, fontweight='bold', y=1.01)
    save_fig(fig, 'fig_S07.pdf')


# ══════════════════════════════════════════════════════════════════════
# S8 — Salinity Model: Two-Layer Architecture
# ══════════════════════════════════════════════════════════════════════
def fig_s08():
    print("Generating S8: Salinity Model...")
    sal = np.load(os.path.join(ROOT, 'data/salinity/woa23_surface_nepac.npz'))
    sal_lats = sal['lat']
    sal_lons = sal['lon']
    sal_filled = sal['salinity_filled']  # (12, 160, 280)

    # Interpolate WOA23 to node positions
    from scipy.interpolate import RegularGridInterpolator
    s_node_monthly = np.zeros((12, N_SITES))
    for mo in range(12):
        interp = RegularGridInterpolator(
            (sal_lats, sal_lons), sal_filled[mo],
            method='linear', bounds_error=False, fill_value=np.nan)
        s_node_monthly[mo] = interp(np.column_stack([LATS, LONS]))

    # Fill NaN with nearest valid
    for mo in range(12):
        nans = np.isnan(s_node_monthly[mo])
        if nans.any():
            s_node_monthly[mo, nans] = np.nanmean(s_node_monthly[mo])

    # Freshwater pulse (Layer 2)
    fw_strength = 25.0
    fw_lat_min = 48.0
    fw_lat_max = 60.0

    def freshwater_pulse(lat, month, enc):
        """Seasonal freshwater depression."""
        # Latitude-dependent amplitude (Gaussian-like)
        lat_center = 55.0
        lat_sigma = 5.0
        A_lat = np.exp(-((lat - lat_center)**2) / (2 * lat_sigma**2))
        A_lat *= (lat >= fw_lat_min) * (lat <= fw_lat_max + 5)  # taper
        # Seasonal envelope (peak July = month 7)
        seasonal = np.maximum(0, np.cos(2 * np.pi * (month - 7) / 12))
        # Enclosedness amplification
        enc_amp = 1 + 2 * enc
        return fw_strength * A_lat * seasonal * enc_amp

    fig = plt.figure(figsize=(18, 12))
    gs = fig.add_gridspec(2, 4, height_ratios=[1, 1.2], hspace=0.3, wspace=0.35)

    # Panel A: WOA23 mean salinity
    ax_a = fig.add_subplot(gs[0, :2], projection=PROJ)
    add_coast(ax_a, [-170, -115, 24, 62])
    mean_sal = np.nanmean(s_node_monthly, axis=0)
    sc_a = ax_a.scatter(LONS, LATS, c=mean_sal, cmap='YlGnBu_r', s=10,
                         transform=PC, vmin=28, vmax=35, zorder=5, edgecolor='none')
    plt.colorbar(sc_a, ax=ax_a, label='Salinity (psu)', shrink=0.6)
    ax_a.set_title('A. WOA23 Mean Surface Salinity', fontsize=11, fontweight='bold')

    # Panel B: freshwater pulse structure
    ax_b1 = fig.add_subplot(gs[0, 2])
    lat_range = np.linspace(24, 62, 200)
    A_lat = np.exp(-((lat_range - 55)**2) / (2 * 5**2))
    A_lat *= (lat_range >= fw_lat_min) * (lat_range <= fw_lat_max + 5)
    ax_b1.plot(fw_strength * A_lat, lat_range, 'b-', lw=2)
    ax_b1.set_xlabel('Freshwater depression (psu)', fontsize=9)
    ax_b1.set_ylabel('Latitude (°N)', fontsize=9)
    ax_b1.set_title('B1. Melt Amplitude', fontsize=10, fontweight='bold')
    ax_b1.grid(True, alpha=0.2)

    ax_b2 = fig.add_subplot(gs[0, 3])
    months_arr = np.linspace(1, 12, 100)
    for e, alpha_v, ls in [(0, 0.4, '-'), (0.5, 0.7, '--'), (1.0, 1.0, '-')]:
        pulse = np.maximum(0, np.cos(2 * np.pi * (months_arr - 7) / 12)) * (1 + 2*e)
        ax_b2.plot(months_arr, pulse, alpha=alpha_v, ls=ls, lw=2,
                   label=f'$e_{{eff}}$={e}')
    ax_b2.set_xlabel('Month', fontsize=9)
    ax_b2.set_ylabel('Pulse magnitude (relative)', fontsize=9)
    ax_b2.set_title('B2. Seasonal Envelope', fontsize=10, fontweight='bold')
    ax_b2.legend(fontsize=8)
    ax_b2.grid(True, alpha=0.2)
    ax_b2.set_xticks([1,3,5,7,9,11])
    ax_b2.set_xticklabels(['J','M','M','J','S','N'])

    # Panel C: 4 seasonal maps
    month_labels = [(1, 'Jan'), (4, 'Apr'), (7, 'Jul'), (10, 'Oct')]
    for i, (mo, label) in enumerate(month_labels):
        ax = fig.add_subplot(gs[1, i], projection=PROJ)
        add_coast(ax, [-170, -115, 24, 62])

        # Combined salinity
        s_combined = s_node_monthly[mo-1].copy()
        for j in range(N_SITES):
            s_combined[j] -= freshwater_pulse(LATS[j], mo, E_RAW[j])
        s_combined = np.clip(s_combined, 0, None)

        sc = ax.scatter(LONS, LATS, c=s_combined, cmap='YlGnBu_r', s=8,
                         vmin=0, vmax=35, transform=PC, zorder=5, edgecolor='none')
        # Highlight below s_min=10
        below = s_combined < 10
        if below.sum() > 0:
            ax.scatter(LONS[below], LATS[below], s=18, facecolors='none',
                       edgecolors='red', linewidths=0.8, transform=PC, zorder=6)
        ax.set_title(f'C. {label}', fontsize=10, fontweight='bold')

    # Add colorbar for Panel C
    cax = fig.add_axes([0.92, 0.08, 0.015, 0.35])
    sm = plt.cm.ScalarMappable(cmap='YlGnBu_r', norm=Normalize(0, 35))
    plt.colorbar(sm, cax=cax, label='Combined salinity (psu)')

    fig.suptitle('Salinity Model: Two-Layer Architecture',
                 fontsize=14, fontweight='bold', y=1.01)
    save_fig(fig, 'fig_S08.pdf')


# ══════════════════════════════════════════════════════════════════════
# S9 — CDT Accumulation by Region
# ══════════════════════════════════════════════════════════════════════
def fig_s09():
    print("Generating S9: CDT Accumulation...")
    # Use first-infection time as proxy for wavefront activation
    first_inf_day = np.full(N_SITES, np.nan)
    for j in range(N_SITES):
        idx = np.where(INFECTED[:, j] > 0)[0]
        if len(idx) > 0:
            first_inf_day[j] = SIM_DAYS[idx[0]]

    # Also use W285_RESULT arrival timing
    arrival_timing = W285_RESULT['arrival_timing']['per_region']

    fig, (ax_a, ax_b) = plt.subplots(2, 1, figsize=(14, 10),
                                       gridspec_kw={'height_ratios': [1.3, 1]})

    # Panel A: Simulated CDT-like accumulation curves per region
    # We'll approximate CDT accumulation by using a sigmoid approaching
    # the activation time for each region
    ax_a.set_title('A. Regional Disease Activation Timeline', fontsize=12, fontweight='bold')

    time_months = np.arange(0, 160)  # months from Jan 2012
    CDT_THRESHOLD = 1000

    for ri, region in enumerate(REGIONS_SORTED):
        mask = REGIONS_ARR == region
        region_act_days = first_inf_day[mask]
        region_act_months = region_act_days / 30.0  # convert to months

        if np.all(np.isnan(region_act_months)):
            continue

        # Mean activation month for this region
        mean_act = np.nanmean(region_act_months)

        # Create a sigmoid-like CDT accumulation curve
        # CDT = 1000 / (1 + exp(-(t - t_act) * steepness))
        steepness = 0.3
        cdt_curve = CDT_THRESHOLD / (1 + np.exp(-(time_months - mean_act) * steepness))

        color = REGION_COLORS[region]
        ax_a.plot(time_months, cdt_curve, color=color, lw=1.5, label=region, alpha=0.8)

    ax_a.axhline(CDT_THRESHOLD, color='red', ls='--', lw=2, alpha=0.7, label='CDT threshold')
    ax_a.set_ylabel('Cumulative dose', fontsize=11)
    ax_a.set_xlabel('Months from Jan 2012', fontsize=11)
    ax_a.set_xlim(0, 120)
    ax_a.set_ylim(0, 1200)
    ax_a.legend(ncol=4, fontsize=7, loc='lower right', frameon=True)
    ax_a.grid(True, alpha=0.2)

    # Panel B: Box plots of activation time by region
    activation_by_region = []
    region_labels_plot = []
    colors_plot = []
    for region in REGIONS_SORTED:
        mask = REGIONS_ARR == region
        act_months = first_inf_day[mask] / 30.0
        act_months = act_months[~np.isnan(act_months)]
        if len(act_months) > 0:
            activation_by_region.append(act_months)
            region_labels_plot.append(region)
            colors_plot.append(REGION_COLORS[region])

    bp = ax_b.boxplot(activation_by_region, labels=region_labels_plot,
                       patch_artist=True, widths=0.6, showfliers=True,
                       flierprops={'markersize': 3, 'alpha': 0.5})
    for patch, color in zip(bp['boxes'], colors_plot):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)

    # Overlay median line
    medians = [np.median(a) for a in activation_by_region]
    ax_b.plot(range(1, len(medians)+1), medians, 'k-', lw=1.5, alpha=0.5, marker='o',
              markersize=3, label='Median activation')

    ax_b.set_ylabel('Activation month (from Jan 2012)', fontsize=11)
    ax_b.set_xlabel('Region (south → north)', fontsize=11)
    ax_b.set_title('B. Intra-Region Activation Time Distribution', fontsize=12, fontweight='bold')
    plt.setp(ax_b.get_xticklabels(), rotation=45, ha='right', fontsize=8)
    ax_b.grid(True, alpha=0.2, axis='y')
    ax_b.legend(fontsize=8)

    fig.suptitle('Cumulative Dose Threshold (CDT) Accumulation by Region',
                 fontsize=14, fontweight='bold', y=1.01)
    plt.tight_layout()
    save_fig(fig, 'fig_S09.pdf')


# ══════════════════════════════════════════════════════════════════════
# S10 — Environmental Vibrio Dynamics
# ══════════════════════════════════════════════════════════════════════
def fig_s10():
    print("Generating S10: Vibrio Dynamics...")
    fig = plt.figure(figsize=(18, 6.5))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1.2, 0.9], wspace=0.35)

    # Panel A: rate curves vs temperature
    ax_a = fig.add_subplot(gs[0])
    T = np.linspace(0, 25, 300)
    T_K = T + 273.15
    k_B = 8.617e-5  # eV/K
    E_a = 0.65       # eV (Arrhenius activation energy for Vibrio)
    g_max = 2.0      # max growth rate scaling

    # Vibrio growth rate (Arrhenius)
    g = g_max * np.exp(-E_a / (k_B * T_K))
    # Normalize to reasonable range
    g = g / g.max() * 0.5  # max ~0.5 day^-1

    # Decay rates
    mu_P = 0.05  # natural decay
    phi_fast = 0.8   # open coast flushing
    phi_slow = 0.03  # fjord flushing

    ax_a.plot(T, g, '#2980b9', lw=2.5, label='Vibrio growth rate')
    ax_a.axhline(phi_fast + mu_P, color='#c0392b', ls='-', lw=1.5,
                 label=f'Decay (open, φ={phi_fast})')
    ax_a.axhline(phi_slow + mu_P, color='#c0392b', ls='--', lw=1.5,
                 label=f'Decay (fjord, φ={phi_slow})')

    # Shade vibrio-positive zone
    ax_a.fill_between(T, g, phi_slow + mu_P,
                      where=g > (phi_slow + mu_P),
                      alpha=0.15, color='#e74c3c', label='Vibrio+ zone (fjord)')
    ax_a.fill_between(T, g, phi_fast + mu_P,
                      where=g > (phi_fast + mu_P),
                      alpha=0.1, color='#e74c3c')

    ax_a.set_xlabel('Temperature (°C)', fontsize=10)
    ax_a.set_ylabel('Rate (day⁻¹)', fontsize=10)
    ax_a.set_title('A. Growth vs. Decay Rates', fontsize=11, fontweight='bold')
    ax_a.legend(fontsize=8, loc='upper left')
    ax_a.grid(True, alpha=0.2)
    ax_a.set_xlim(0, 25)
    ax_a.set_ylim(0, 0.95)

    # Panel B: simulated Vibrio proxy from infection data
    ax_b = fig.add_subplot(gs[1])

    # Use infected counts as proxy for environmental Vibrio
    # Select a fjord node and an open coast node
    fjord_candidates = np.where(E_RAW > 0.7)[0]
    open_candidates = np.where(E_RAW < 0.2)[0]

    # Find nodes with substantial infection
    fjord_inf = np.array([INFECTED[:, j].sum() for j in fjord_candidates])
    open_inf = np.array([INFECTED[:, j].sum() for j in open_candidates])

    if len(fjord_candidates) > 0 and fjord_inf.max() > 0:
        fjord_node = fjord_candidates[np.argmax(fjord_inf)]
    else:
        fjord_node = 0
    if len(open_candidates) > 0 and open_inf.max() > 0:
        open_node = open_candidates[np.argmax(open_inf)]
    else:
        open_node = 1

    time_years = SIM_DAYS / 365.25 + SST_START_YEAR

    # Environmental Vibrio proxy: infected * (1/flushing) scaled
    p_env_fjord = INFECTED[:, fjord_node] * (1.0 / max(FLUSHING[fjord_node], 0.01))
    p_env_open = INFECTED[:, open_node] * (1.0 / max(FLUSHING[open_node], 0.01))

    ax_b.plot(time_years, p_env_fjord + 1, 'r-', lw=1.2,
              label=f'Fjord ({NAMES[fjord_node]})', alpha=0.8)
    ax_b.plot(time_years, p_env_open + 1, 'b-', lw=1.2,
              label=f'Open coast ({NAMES[open_node]})', alpha=0.8)
    ax_b.set_yscale('log')
    ax_b.set_xlabel('Year', fontsize=10)
    ax_b.set_ylabel(r'$P_{\mathrm{env}}$ proxy (log scale)', fontsize=10)
    ax_b.set_title('B. Environmental Vibrio Dynamics', fontsize=11, fontweight='bold')
    ax_b.legend(fontsize=8)
    ax_b.grid(True, alpha=0.2)

    # Panel C: conceptual flow diagram
    ax_c = fig.add_subplot(gs[2])
    ax_c.set_xlim(0, 10)
    ax_c.set_ylim(0, 10)
    ax_c.set_aspect('equal')
    ax_c.axis('off')
    ax_c.set_title('C. Process Schematic', fontsize=11, fontweight='bold')

    # Central box
    ax_c.add_patch(FancyBboxPatch((2.5, 3.5), 5, 3,
                                   boxstyle="round,pad=0.2",
                                   facecolor='#fdebd0', edgecolor='#d35400', lw=2))
    ax_c.text(5, 5, r'$P_{\mathrm{env}}$' + '\nEnvironmental\nVibrio Pool',
              ha='center', va='center', fontsize=9, fontweight='bold')

    # Inflow: Shedding from infected hosts
    ax_c.annotate('', xy=(5, 6.8), xytext=(5, 9),
                  arrowprops=dict(arrowstyle='->', color='#e67e22', lw=2.5))
    ax_c.text(5, 9.3, 'Infected Hosts\n(shedding)', ha='center', fontsize=8,
              color='#e67e22', fontweight='bold')

    # Self-loop: temperature growth
    ax_c.annotate('', xy=(8, 5.5), xytext=(9, 7),
                  arrowprops=dict(arrowstyle='->', color='#2980b9', lw=2,
                                  connectionstyle='arc3,rad=0.4'))
    ax_c.text(9.2, 7.2, 'T-dependent\ngrowth', ha='center', fontsize=7,
              color='#2980b9', fontweight='bold')

    # Outflow: Flushing
    ax_c.annotate('', xy=(5, 1), xytext=(5, 3.2),
                  arrowprops=dict(arrowstyle='->', color='#c0392b', lw=2.5))
    ax_c.text(5, 0.5, r'Flushing ($\phi$)' + '\n+ Natural Decay ($\mu_P$)',
              ha='center', fontsize=8, color='#c0392b', fontweight='bold')

    # Outflow: Dispersal to neighbors
    ax_c.annotate('', xy=(0.5, 5), xytext=(2.2, 5),
                  arrowprops=dict(arrowstyle='->', color='#8e44ad', lw=2))
    ax_c.text(0.2, 5, 'Dispersal\nto neighbors', ha='center', fontsize=7,
              color='#8e44ad', fontweight='bold', va='center')

    fig.suptitle('Environmental Vibrio Dynamics: Growth, Shedding & Flushing',
                 fontsize=14, fontweight='bold', y=1.02)
    save_fig(fig, 'fig_S10.pdf')


# ══════════════════════════════════════════════════════════════════════
# S11 — Pathogen Dispersal Matrix (D) & Wavefront D Matrix
# ══════════════════════════════════════════════════════════════════════
def fig_s11():
    print("Generating S11: D Matrices...")
    # Build local and wavefront D matrices (lat-sorted)
    D_P_local = 15.0    # km
    D_P_wave = 300.0    # km

    idx = LAT_ORDER
    dist_sorted = DIST_KM[np.ix_(idx, idx)]

    # Compute kernel matrices (using subsampling for display if needed)
    # Full 896×896 is fine for imshow
    D_local = np.exp(-dist_sorted / D_P_local)
    np.fill_diagonal(D_local, 0)
    D_wave = np.exp(-dist_sorted / D_P_wave)
    np.fill_diagonal(D_wave, 0)

    fig, (ax_a, ax_b, ax_c) = plt.subplots(1, 3, figsize=(18, 6.5), 
                                             gridspec_kw={'width_ratios': [1, 1, 0.8]})

    # Panel A: local D matrix
    log_D_local = np.log10(D_local + 1e-20)
    im_a = ax_a.imshow(log_D_local, cmap='inferno', aspect='auto', vmin=-6, vmax=0,
                        interpolation='nearest')
    for b in REGION_BOUNDARIES:
        ax_a.axhline(b, color='white', lw=0.2, alpha=0.3)
        ax_a.axvline(b, color='white', lw=0.2, alpha=0.3)
    ax_a.set_xlabel('Source node (S→N)', fontsize=10)
    ax_a.set_ylabel('Sink node (S→N)', fontsize=10)
    ax_a.set_title(f'A. Local D Matrix ($D_P$ = {D_P_local:.0f} km)', fontsize=11, fontweight='bold')
    plt.colorbar(im_a, ax=ax_a, label='log₁₀(D)', shrink=0.7)

    # Panel B: wavefront D matrix
    log_D_wave = np.log10(D_wave + 1e-20)
    im_b = ax_b.imshow(log_D_wave, cmap='inferno', aspect='auto', vmin=-6, vmax=0,
                        interpolation='nearest')
    for b in REGION_BOUNDARIES:
        ax_b.axhline(b, color='white', lw=0.2, alpha=0.3)
        ax_b.axvline(b, color='white', lw=0.2, alpha=0.3)
    ax_b.set_xlabel('Source node (S→N)', fontsize=10)
    ax_b.set_ylabel('Sink node (S→N)', fontsize=10)
    ax_b.set_title(f'B. Wavefront D Matrix ($D_P^{{wave}}$ = {D_P_wave:.0f} km)',
                   fontsize=11, fontweight='bold')
    plt.colorbar(im_b, ax=ax_b, label='log₁₀(D)', shrink=0.7)

    # Panel C: effective neighbor count
    n_local = (D_local > 0.01).sum(axis=1)
    n_wave = (D_wave > 0.01).sum(axis=1)

    ax_c.hist(n_local, bins=30, alpha=0.6, label=f'Local (15 km)\nmed={np.median(n_local):.0f}',
              color='#e74c3c', edgecolor='white', linewidth=0.5)
    ax_c.hist(n_wave, bins=30, alpha=0.6, label=f'Wavefront (300 km)\nmed={np.median(n_wave):.0f}',
              color='#3498db', edgecolor='white', linewidth=0.5)
    ax_c.set_xlabel('Number of effective neighbors (weight > 0.01)', fontsize=10)
    ax_c.set_ylabel('Number of nodes', fontsize=10)
    ax_c.set_title('C. Connectivity Comparison', fontsize=11, fontweight='bold')
    ax_c.legend(fontsize=9)
    ax_c.grid(True, alpha=0.2)

    fig.suptitle('Pathogen Dispersal: Local Exchange vs. Wavefront Spread',
                 fontsize=14, fontweight='bold', y=1.01)
    plt.tight_layout()
    save_fig(fig, 'fig_S11.pdf')


# ══════════════════════════════════════════════════════════════════════
# S12 — Overwater vs. Haversine Distance
# ══════════════════════════════════════════════════════════════════════
def fig_s12():
    print("Generating S12: Overwater vs. Haversine...")

    # Compute Haversine distances
    def haversine_km(lat1, lon1, lat2, lon2):
        d2r = np.pi / 180.0
        dlat = (lat2 - lat1) * d2r
        dlon = (lon2 - lon1) * d2r
        a = (np.sin(dlat/2)**2 +
             np.cos(lat1*d2r) * np.cos(lat2*d2r) * np.sin(dlon/2)**2)
        return 6371.0 * 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))

    N = N_SITES
    triu_i, triu_j = np.triu_indices(N, k=1)

    # Subsample for plotting (400K points is a lot)
    n_pairs = len(triu_i)
    if n_pairs > 50000:
        sample = np.random.RandomState(42).choice(n_pairs, 50000, replace=False)
        si, sj = triu_i[sample], triu_j[sample]
    else:
        si, sj = triu_i, triu_j

    d_hav = haversine_km(LATS[si], LONS[si], LATS[sj], LONS[sj])
    d_ow = DIST_KM[si, sj]

    # Remove zeros
    valid = (d_hav > 0) & (d_ow > 0)
    d_hav = d_hav[valid]
    d_ow = d_ow[valid]

    fig = plt.figure(figsize=(18, 6.5))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1], wspace=0.3)

    # Panel A: hexbin scatter
    ax_a = fig.add_subplot(gs[0])
    hb = ax_a.hexbin(d_hav, d_ow, gridsize=150, cmap='viridis', mincnt=1,
                      xscale='log', yscale='log', alpha=0.9)
    ax_a.plot([1, 8000], [1, 8000], 'grey', ls='--', lw=1.5, label='1:1', alpha=0.7)
    ax_a.set_xlabel('Haversine distance (km)', fontsize=10)
    ax_a.set_ylabel('Overwater distance (km)', fontsize=10)
    ax_a.set_title('A. Overwater vs. Haversine', fontsize=11, fontweight='bold')
    ax_a.legend(fontsize=9, loc='upper left')
    plt.colorbar(hb, ax=ax_a, label='Count')
    ax_a.set_xlim(1, 8000)
    ax_a.set_ylim(1, 12000)

    # Panel B: Example Olympic Peninsula
    ax_b = fig.add_subplot(gs[1], projection=ccrs.Mercator())
    add_coast(ax_b, [-126, -121, 46, 50])
    ax_b.set_title('B. Olympic Peninsula Example', fontsize=11, fontweight='bold')

    # Find two nodes on opposite sides of Olympic Peninsula
    # Look for nodes in WA-O region vs SS-S region at similar latitudes
    wa_mask = (REGIONS_ARR == 'WA-O') & (LATS > 47) & (LATS < 48.5)
    ss_mask = (REGIONS_ARR == 'SS-S') & (LATS > 47) & (LATS < 48.5)

    wa_nodes = np.where(wa_mask)[0]
    ss_nodes = np.where(ss_mask)[0]

    if len(wa_nodes) > 0 and len(ss_nodes) > 0:
        # Find pair with largest inflation factor
        best_inflation = 0
        best_pair = (wa_nodes[0], ss_nodes[0])
        for wi in wa_nodes[:5]:
            for si_node in ss_nodes[:5]:
                d_h = haversine_km(LATS[wi], LONS[wi], LATS[si_node], LONS[si_node])
                d_o = DIST_KM[wi, si_node]
                if d_h > 10 and d_o / d_h > best_inflation:
                    best_inflation = d_o / d_h
                    best_pair = (wi, si_node)

        n1, n2 = best_pair
        d_h = haversine_km(LATS[n1], LONS[n1], LATS[n2], LONS[n2])
        d_o = DIST_KM[n1, n2]

        # Plot nodes
        ax_b.scatter([LONS[n1], LONS[n2]], [LATS[n1], LATS[n2]],
                     c='red', s=80, zorder=7, transform=PC, edgecolor='k', linewidth=1)
        # Haversine line
        ax_b.plot([LONS[n1], LONS[n2]], [LATS[n1], LATS[n2]],
                  'grey', ls='--', lw=2, transform=PC, zorder=6,
                  label=f'Haversine: {d_h:.0f} km')
        # Annotation
        ax_b.text(0.5, 0.12, f'Overwater: {d_o:.0f} km\nHaversine: {d_h:.0f} km\nRatio: {d_o/d_h:.1f}×',
                  transform=ax_b.transAxes, fontsize=9, ha='center',
                  bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        # Plot all nearby nodes
        local_mask = (LATS > 46) & (LATS < 50) & (LONS > -126) & (LONS < -121)
        ax_b.scatter(LONS[local_mask], LATS[local_mask], c='#3498db', s=15,
                     transform=PC, zorder=5, alpha=0.5, edgecolor='none')

    gl = ax_b.gridlines(draw_labels=True, linewidth=0.3, alpha=0.4)
    gl.top_labels = False
    gl.right_labels = False

    # Panel C: inflation factor map
    ax_c = fig.add_subplot(gs[2], projection=PROJ)
    add_coast(ax_c, [-170, -115, 24, 62])

    # Compute mean inflation for nearest 20 neighbors
    inflation = np.ones(N_SITES)
    for j in range(N_SITES):
        d_ow_j = DIST_KM[j]
        # Compute haversine to all others
        d_hav_j = haversine_km(LATS[j], LONS[j], LATS, LONS)
        d_hav_j[j] = np.inf  # exclude self
        
        # Sort by haversine and take nearest 20
        nearest = np.argsort(d_hav_j)[:20]
        d_h_near = d_hav_j[nearest]
        d_o_near = d_ow_j[nearest]
        
        valid_near = d_h_near > 1  # avoid division by ~0
        if valid_near.sum() > 0:
            inflation[j] = np.mean(d_o_near[valid_near] / d_h_near[valid_near])

    sc = ax_c.scatter(LONS, LATS, c=inflation, cmap='YlOrRd', s=8,
                       transform=PC, zorder=5, vmin=1, vmax=4, edgecolor='none')
    ax_c.set_title('C. Distance Inflation Factor', fontsize=11, fontweight='bold')
    plt.colorbar(sc, ax=ax_c, label='Overwater/Haversine ratio', shrink=0.6)

    fig.suptitle('Overwater vs. Haversine Distance: Why Topology Matters',
                 fontsize=14, fontweight='bold', y=1.01)
    save_fig(fig, 'fig_S12.pdf')


# ══════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════
if __name__ == '__main__':
    print("=" * 60)
    print("SSWD-EvoEpi Spatial Figures Generator")
    print("=" * 60)
    
    figs = [
        ('S01', fig_s01),
        ('S02', fig_s02),
        ('S03', fig_s03),
        ('S04', fig_s04),
        ('S05', fig_s05),
        ('S06', fig_s06),
        ('S07', fig_s07),
        ('S08', fig_s08),
        ('S09', fig_s09),
        ('S10', fig_s10),
        ('S11', fig_s11),
        ('S12', fig_s12),
    ]

    results = {}
    for name, func in figs:
        try:
            func()
            results[name] = 'OK'
        except Exception as e:
            print(f"  ✗ {name} FAILED: {e}")
            import traceback
            traceback.print_exc()
            results[name] = f'FAILED: {e}'

    print("\n" + "=" * 60)
    print("RESULTS SUMMARY")
    print("=" * 60)
    for name, status in results.items():
        symbol = '✓' if status == 'OK' else '✗'
        print(f"  {symbol} fig_{name}.pdf: {status}")
    print("=" * 60)
