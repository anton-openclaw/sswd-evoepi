"""Shared utilities for W325 figure generation."""
import json
import numpy as np

REPO = "/home/starbot/.openclaw/workspace/sswd-evoepi"
RESULT = f"{REPO}/results/k1000_scaled_sweep/W325/result_seed42.json"
NPZ = f"{REPO}/results/k1000_scaled_sweep/W325/monthly_seed42.npz"
FIGDIR = f"{REPO}/reports/w325/figures"

COASTLINE_ORDER = [
    "BJ", "CA-S", "CA-C", "CA-N", "OR", "WA-O", "JDF", "SS-S", "SS-N",
    "BC-C", "BC-N", "AK-FS", "AK-FN", "AK-OC", "AK-PWS", "AK-EG", "AK-WG", "AK-AL"
]

SCORED_REGIONS = ["AK-PWS", "AK-FN", "AK-FS", "BC-N", "SS-S", "JDF", "OR", "CA-N"]

# Scored regions in coastline order (S→N) for bar charts
SCORED_COASTLINE = [r for r in COASTLINE_ORDER if r in SCORED_REGIONS]

TARGETS = {
    "AK-PWS": 0.20, "AK-FN": 0.20, "AK-FS": 0.20, "BC-N": 0.20,
    "SS-S": 0.05, "JDF": 0.02, "OR": 0.0025, "CA-N": 0.001,
}

K_PER_SITE = 1000
START_YEAR = 2012
YEARS = list(range(2012, 2025))  # 13 years of yearly data

PALETTE = {
    'disease': 'Reds',
    'population': 'YlGn',
    'latitude': 'coolwarm',
    'timing': 'plasma',
    'crash': 'RdYlGn_r',
    'heatmap': 'inferno',
    'evolution': 'viridis',
    'pathogen': 'magma',
}

def load_result():
    with open(RESULT) as f:
        return json.load(f)

def load_npz():
    data = np.load(NPZ, allow_pickle=True)
    return data

def get_regions_per_site(site_names):
    """Map site names to region names.
    
    Site names are either 'REGION-NNN' (e.g., 'BJ-001', 'JDF-023', 'OR-006')
    or 'REGION-SUB-NNN' (e.g., 'AK-PWS-001', 'CA-N-042').
    The region is everything before the final '-NNN' numeric suffix.
    """
    regions = []
    for n in site_names:
        parts = str(n).split("-")
        # Last part is always the numeric site ID — strip it
        if len(parts) >= 2 and parts[-1].isdigit():
            regions.append("-".join(parts[:-1]))
        else:
            # Fallback: take first two parts (shouldn't happen with valid data)
            regions.append("-".join(parts[:2]))
    return regions

def aggregate_to_regions(data_2d, regions_per_site, coastline_order=None):
    """Aggregate (time, sites) → (time, regions) by summing."""
    if coastline_order is None:
        coastline_order = COASTLINE_ORDER
    n_time = data_2d.shape[0]
    result = np.zeros((n_time, len(coastline_order)))
    for j, reg in enumerate(coastline_order):
        mask = np.array([r == reg for r in regions_per_site])
        if mask.any():
            result[:, j] = data_2d[:, mask].sum(axis=1)
    return result

def get_region_lats(site_names, site_lats):
    regions_per_site = get_regions_per_site(site_names)
    region_lats = {}
    for r, lat in zip(regions_per_site, site_lats):
        region_lats.setdefault(r, []).append(lat)
    return {r: np.mean(v) for r, v in region_lats.items()}

def get_region_site_counts(site_names):
    regions_per_site = get_regions_per_site(site_names)
    counts = {}
    for r in regions_per_site:
        counts[r] = counts.get(r, 0) + 1
    return counts

def setup_style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.style.use("default")
    plt.rcParams.update({
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "font.size": 10,
        "axes.titlesize": 12,
        "axes.labelsize": 11,
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.1,
    })
    return plt
