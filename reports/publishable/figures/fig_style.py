"""Shared style and utilities for all publishable figures."""
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import json, pathlib

# ── Paths ──────────────────────────────────────────────────────────────
ROOT = pathlib.Path(__file__).resolve().parent.parent.parent.parent  # sswd-evoepi/
DATA_JSON = ROOT / "reports" / "publishable" / "data" / "summary_data.json"
NPZ_DIR = ROOT / "results" / "calibration"
FIG_DIR = pathlib.Path(__file__).resolve().parent

# ── Style ──────────────────────────────────────────────────────────────
def apply_style():
    try:
        plt.style.use('seaborn-v0_8-paper')
    except OSError:
        plt.style.use('seaborn-paper')
    mpl.rcParams.update({
        'font.size': 10,
        'axes.labelsize': 10,
        'axes.titlesize': 11,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'legend.fontsize': 8,
        'figure.dpi': 300,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.1,
        'axes.spines.top': False,
        'axes.spines.right': False,
    })

# ── Scenario colours & labels ──────────────────────────────────────────
SCENARIO_COLORS = {
    'F01': 'tab:blue',
    'F02': 'tab:orange',
    'F03': 'tab:green',
    'F04': 'tab:red',
}
SCENARIO_LABELS = {
    'F01': 'SSP2-4.5, all evo ON',
    'F02': 'SSP2-4.5, all evo OFF',
    'F03': 'SSP2-4.5, thermal adpt ON',
    'F04': 'SSP5-8.5, all evo ON',
}

# ── Region ordering (N→S) ─────────────────────────────────────────────
REGIONS_NS = [
    'AK-AL', 'AK-WG', 'AK-OC', 'AK-EG', 'AK-PWS', 'AK-FN', 'AK-FS',
    'BC-N', 'BC-C', 'SS-N', 'SS-S', 'JDF', 'WA-O', 'OR',
    'CA-N', 'CA-C', 'CA-S', 'BJ',
]
KEY_REGIONS = ['AK-PWS', 'AK-FN', 'BC-N', 'OR', 'CA-C', 'CA-S']

SEEDS = [42, 123, 999]

# ── Helper: load summary JSON ─────────────────────────────────────────
_cache = {}
def load_summary():
    if 'summary' not in _cache:
        with open(DATA_JSON) as f:
            _cache['summary'] = json.load(f)
    return _cache['summary']

# ── Helper: years array from yearly data ───────────────────────────────
def years_from_yearly(n):
    """Return year array starting at 2012 with length n (one per simulation year)."""
    return np.arange(2012, 2012 + n)

# ── Helper: load monthly npz and compute region aggregates ─────────────
def load_monthly_regional(scenario, seed):
    """Return dict with region -> population time-series (monthly)."""
    npz_path = NPZ_DIR / scenario / f"monthly_seed{seed}.npz"
    d = np.load(npz_path, allow_pickle=True)
    names = d['site_names']
    pops = d['populations']
    sim_days = d['sim_days']
    # Build region -> indices
    region_idx = {}
    for i, name in enumerate(names):
        parts = str(name).split('-')
        # Region prefix is everything except the last numeric segment
        # e.g. AK-AL-001 -> AK-AL, OR-048 -> OR, JDF-005 -> JDF, BJ-031 -> BJ
        prefix = '-'.join(p for p in parts if not p.isdigit())
        region_idx.setdefault(prefix, []).append(i)
    result = {}
    for reg in REGIONS_NS:
        idx = region_idx.get(reg, [])
        if idx:
            result[reg] = pops[:, idx].sum(axis=1)
        else:
            result[reg] = np.zeros(len(sim_days))
    result['_sim_days'] = sim_days
    result['_years'] = 2012 + sim_days / 365.25
    return result

def savefig(fig, name):
    out = FIG_DIR / name
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"  ✓ {out.name}")
