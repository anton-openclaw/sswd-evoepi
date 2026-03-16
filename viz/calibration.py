#!/usr/bin/env python3
"""Calibration report visualization module.

Reusable plots for calibration sweep analysis. Uses publication style
from fig_style.py. Designed for multi-config × multi-seed sweeps.

Usage:
    from viz.calibration import CalibrationData, plot_rmse_ranking, ...

    data = CalibrationData.from_results_dir('results/calibration/', configs=['W197', 'W198'])
    fig = plot_rmse_ranking(data, baseline_rmse=0.608, baseline_label='W173')
    fig.savefig('figures/rmse_ranking.pdf')
"""
from __future__ import annotations

import json
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Patch
from pathlib import Path
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Any

# Ensure project root is on path for sswd_evoepi imports
_PROJECT_ROOT = str(Path(__file__).resolve().parent.parent)
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)

from sswd_evoepi.results import REGION_ORDER, SCORED_REGIONS   # noqa: E402
from sswd_evoepi.metrics import RECOVERY_TARGETS               # noqa: E402


# ── Style ──────────────────────────────────────────────────────────────

def apply_pub_style():
    """Apply publication-quality light style for PDF reports."""
    try:
        plt.style.use('seaborn-v0_8-paper')
    except OSError:
        plt.style.use('seaborn-paper')
    matplotlib.rcParams.update({
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


# ── Region definitions (from consolidated modules) ─────────────────────
# REGION_ORDER, SCORED_REGIONS imported from sswd_evoepi.results
# RECOVERY_TARGETS imported from sswd_evoepi.metrics

# Local aliases for backward compatibility
REGIONS_NS = REGION_ORDER
DEFAULT_TARGETS = RECOVERY_TARGETS


# ── Group coloring ─────────────────────────────────────────────────────

DEFAULT_GROUP_COLORS = {
    'default':    '#607D8B',
    'group_a':    '#2196F3',
    'group_b':    '#FF9800',
    'group_c':    '#4CAF50',
    'group_d':    '#9C27B0',
    'group_e':    '#F44336',
    'group_f':    '#FF5722',
    'group_g':    '#00BCD4',
    'group_h':    '#795548',
}


# ── Data container ─────────────────────────────────────────────────────

@dataclass
class CalibrationData:
    """Container for calibration sweep results.

    Attributes:
        configs: ordered list of config names (e.g. ['W197', 'W198'])
        seeds: list of seed values
        results: config_name -> seed -> result dict (from JSON)
        descriptions: config_name -> short description string
        groups: config_name -> group label (for coloring)
        group_colors: group_label -> hex color
    """
    configs: List[str]
    seeds: List[int]
    results: Dict[str, Dict[int, dict]]
    descriptions: Dict[str, str] = field(default_factory=dict)
    groups: Dict[str, str] = field(default_factory=dict)
    group_colors: Dict[str, str] = field(default_factory=lambda: dict(DEFAULT_GROUP_COLORS))

    @classmethod
    def from_results_dir(
        cls,
        results_dir: str | Path,
        configs: Optional[List[str]] = None,
        seeds: List[int] = (42, 123, 999),
    ) -> 'CalibrationData':
        """Load calibration results from a directory structure.

        Expected layout:
            results_dir/W199/result_seed42.json
            results_dir/W199/result_seed123.json
            ...
        """
        results_dir = Path(results_dir)
        if configs is None:
            configs = sorted([
                d.name for d in results_dir.iterdir()
                if d.is_dir() and d.name.startswith('W')
            ], key=lambda x: int(x[1:]))

        results = {}
        for config in configs:
            results[config] = {}
            for seed in seeds:
                fpath = results_dir / config / f'result_seed{seed}.json'
                if fpath.exists():
                    with open(fpath) as f:
                        results[config][seed] = json.load(f)

        return cls(configs=configs, seeds=list(seeds), results=results)

    def rmse_per_seed(self, config: str) -> List[float]:
        """Get RMSE values for each seed of a config."""
        vals = []
        for s in self.seeds:
            if s in self.results.get(config, {}):
                v = self.results[config][s]['scoring'].get('rmsle', self.results[config][s]['scoring'].get('rmse_log'))
                vals.append(v if v is not None and v != float('inf') else np.inf)
            else:
                vals.append(np.nan)
        return vals

    def mean_rmse(self, config: str) -> float:
        """Mean RMSE across seeds (inf-aware)."""
        vals = self.rmse_per_seed(config)
        finite = [v for v in vals if np.isfinite(v)]
        return np.mean(finite) if finite else np.inf

    def recovery(self, config: str, region: str, seed: Optional[int] = None) -> float:
        """Get recovery fraction for a config/region. Mean across seeds if seed=None."""
        if seed is not None:
            return self.results.get(config, {}).get(seed, {}).get('region_recovery', {}).get(region, 0.0)
        vals = []
        for s in self.seeds:
            v = self.results.get(config, {}).get(s, {}).get('region_recovery', {}).get(region, 0.0)
            vals.append(v)
        return np.mean(vals)

    def sorted_configs(self) -> List[str]:
        """Return configs sorted by mean RMSE (best first)."""
        return sorted(self.configs, key=lambda c: self.mean_rmse(c))

    def color(self, config: str) -> str:
        """Get color for a config based on its group."""
        group = self.groups.get(config, 'default')
        return self.group_colors.get(group, '#607D8B')

    def to_summary_json(self, path: str | Path):
        """Export data_summary.json for report skill."""
        summary = {
            'configs': {},
            'rankings': [],
            'targets': dict(DEFAULT_TARGETS),
            'seeds': self.seeds,
        }
        for config in self.sorted_configs():
            rmse_vals = self.rmse_per_seed(config)
            recovery = {}
            for reg in REGIONS_NS:
                recovery[reg] = self.recovery(config, reg)

            summary['configs'][config] = {
                'description': self.descriptions.get(config, ''),
                'group': self.groups.get(config, 'default'),
                'rmse_mean': float(self.mean_rmse(config)),
                'rmse_std': float(np.std([v for v in rmse_vals if np.isfinite(v)]) if any(np.isfinite(v) for v in rmse_vals) else 0),
                'rmse_per_seed': {str(s): float(v) for s, v in zip(self.seeds, rmse_vals)},
                'recovery': {reg: float(v) for reg, v in recovery.items()},
            }
            summary['rankings'].append(config)

        with open(path, 'w') as f:
            json.dump(summary, f, indent=2)
        print(f"Wrote {path}")


# ── NPZ loading ────────────────────────────────────────────────────────

def load_npz_regional(npz_path: str | Path) -> Dict[str, Any]:
    """Load monthly NPZ and return regional population time series.

    Returns dict with region -> population array, plus '_years' and '_sim_days'.
    """
    d = np.load(npz_path, allow_pickle=True)
    names = d['site_names']
    pops = d['populations']
    sim_days = d['sim_days']
    start_year = int(d.get('sst_start_year', 2012))

    region_idx = {}
    for i, name in enumerate(names):
        parts = str(name).split('-')
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
    result['_years'] = start_year + sim_days / 365.25
    return result


# ══════════════════════════════════════════════════════════════════════
# PLOT FUNCTIONS
# ══════════════════════════════════════════════════════════════════════

def plot_rmse_ranking(
    data: CalibrationData,
    baseline_rmse: Optional[float] = None,
    baseline_label: str = 'Baseline',
    max_rmse: float = 1.5,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Horizontal bar chart of configs ranked by mean RMSE.

    Args:
        data: CalibrationData instance
        baseline_rmse: optional reference line (e.g. previous best)
        baseline_label: label for reference line
        max_rmse: truncate x-axis at this value
        save_path: optional path to save figure
    """
    apply_pub_style()
    configs = data.sorted_configs()
    n = len(configs)

    fig, ax = plt.subplots(figsize=(10, max(5, 0.5 * n)))
    y_pos = np.arange(n)
    means = [data.mean_rmse(c) for c in configs]
    stds = [np.std([v for v in data.rmse_per_seed(c) if np.isfinite(v)]) for c in configs]
    colors = [data.color(c) for c in configs]

    ax.barh(y_pos, [min(m, max_rmse) for m in means], xerr=stds,
            color=colors, edgecolor='white', height=0.7, capsize=3)

    if baseline_rmse is not None:
        ax.axvline(x=baseline_rmse, color='red', linestyle='--', linewidth=1.5,
                   label=f'{baseline_label} ({baseline_rmse:.3f})')

    labels = []
    for c in configs:
        desc = data.descriptions.get(c, '')
        labels.append(f'{c}\n{desc}' if desc else c)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel('RMSE (log scale)')
    ax.set_title('Calibration Sweep: RMSE Rankings (mean ± σ)')
    ax.invert_yaxis()
    ax.set_xlim(0, max_rmse)

    # Value labels
    for i, (m, s) in enumerate(zip(means, stds)):
        if m < max_rmse:
            ax.text(min(m + s + 0.02, max_rmse - 0.05), i, f'{m:.3f}', va='center', fontsize=8)
        else:
            ax.text(max_rmse - 0.15, i, f'{m:.2f}', va='center', fontsize=8, color='red')

    # Group legend
    seen = {}
    for c in configs:
        g = data.groups.get(c, 'default')
        if g not in seen:
            seen[g] = data.color(c)
    legend_elements = [Patch(facecolor=col, label=g) for g, col in seen.items()]
    if baseline_rmse is not None:
        legend_elements.append(plt.Line2D([0], [0], color='red', linestyle='--', label=baseline_label))
    ax.legend(handles=legend_elements, loc='lower right', fontsize=7)

    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, bbox_inches='tight')
    return fig


def plot_region_heatmap(
    data: CalibrationData,
    regions: Optional[List[str]] = None,
    targets: Optional[Dict[str, float]] = None,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Heatmap of regional recovery % (configs × regions).

    Args:
        data: CalibrationData instance
        regions: list of regions to show (default: SCORED_REGIONS)
        targets: recovery targets dict (default: DEFAULT_TARGETS)
        save_path: optional path to save figure
    """
    apply_pub_style()
    if regions is None:
        regions = SCORED_REGIONS
    if targets is None:
        targets = DEFAULT_TARGETS

    configs = data.sorted_configs()
    n_configs = len(configs)
    n_regions = len(regions)

    matrix = np.zeros((n_configs, n_regions))
    for i, c in enumerate(configs):
        for j, reg in enumerate(regions):
            matrix[i, j] = data.recovery(c, reg) * 100

    fig, ax = plt.subplots(figsize=(max(8, n_regions * 1.2), max(5, n_configs * 0.5)))
    im = ax.imshow(matrix, cmap='YlOrRd', aspect='auto', vmin=0, vmax=50)

    ax.set_xticks(range(n_regions))
    ax.set_xticklabels(regions, rotation=45, ha='right')
    ax.set_yticks(range(n_configs))
    ax.set_yticklabels([f'{c} ({data.mean_rmse(c):.3f})' for c in configs], fontsize=9)

    # Annotate cells
    for i in range(n_configs):
        for j in range(n_regions):
            val = matrix[i, j]
            target_pct = targets.get(regions[j], 0) * 100
            color = 'white' if val > 25 else 'black'
            text = f'{val:.1f}' if val >= 1 else f'{val:.2f}'
            fontweight = 'bold' if val >= target_pct * 0.8 else 'normal'
            ax.text(j, i, text, ha='center', va='center', fontsize=7,
                    color=color, fontweight=fontweight)

    ax.set_title('Regional Recovery (%) — bold = within 80% of target')
    plt.colorbar(im, ax=ax, label='Recovery %', shrink=0.8)

    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, bbox_inches='tight')
    return fig


def plot_tradeoff(
    data: CalibrationData,
    x_region: str = 'AK-PWS',
    y_region: str = 'CA-N',
    x_target: Optional[float] = None,
    y_target: Optional[float] = None,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Scatter: recovery in one region vs another (the latitude tradeoff).

    Args:
        data: CalibrationData instance
        x_region: region for x-axis
        y_region: region for y-axis
        x_target: target for x-region (draws vertical line)
        y_target: target for y-region (draws horizontal line)
        save_path: optional path to save figure
    """
    apply_pub_style()
    if x_target is None:
        x_target = DEFAULT_TARGETS.get(x_region)
    if y_target is None:
        y_target = DEFAULT_TARGETS.get(y_region)

    fig, ax = plt.subplots(figsize=(9, 7))

    for c in data.configs:
        x_val = data.recovery(c, x_region) * 100
        y_val = data.recovery(c, y_region) * 100
        color = data.color(c)
        ax.scatter(x_val, y_val, color=color, s=100, zorder=5,
                   edgecolors='black', linewidth=0.5)
        ax.annotate(c, (x_val, y_val), textcoords='offset points',
                    xytext=(5, 5), fontsize=8)

    if x_target is not None:
        ax.axvline(x=x_target * 100, color='green', linestyle=':', alpha=0.5,
                   label=f'{x_region} target ({x_target*100:.0f}%)')
    if y_target is not None:
        ax.axhline(y=y_target * 100, color='red', linestyle=':', alpha=0.5,
                   label=f'{y_region} target ({y_target*100:.1f}%)')

    ax.set_xlabel(f'{x_region} Recovery (%)')
    ax.set_ylabel(f'{y_region} Recovery (%)')
    ax.set_title(f'Latitude Gradient Tradeoff\n{x_region} vs {y_region}')
    ax.legend()

    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, bbox_inches='tight')
    return fig


def plot_parameter_effects(
    data: CalibrationData,
    param_configs: Dict[str, List[Tuple[str, float]]],
    baseline_rmse: Optional[float] = None,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Multi-panel parameter sensitivity plot.

    Args:
        data: CalibrationData instance
        param_configs: dict of param_name -> [(config_name, param_value), ...]
            Example: {'s₀': [('W185', 0.2), ('W187', 0.5), ('W189', 1.0)]}
        baseline_rmse: optional reference line
        save_path: optional path to save figure
    """
    apply_pub_style()
    n_params = len(param_configs)
    cols = min(2, n_params)
    rows = (n_params + cols - 1) // cols

    fig, axes = plt.subplots(rows, cols, figsize=(6 * cols, 5 * rows))
    if n_params == 1:
        axes = np.array([axes])
    axes = axes.flatten()

    for idx, (param_name, config_vals) in enumerate(param_configs.items()):
        ax = axes[idx]
        x_vals = [v for _, v in config_vals]
        y_vals = [data.mean_rmse(c) for c, _ in config_vals]
        ax.plot(x_vals, y_vals, 'o-', markersize=10, color='#2196F3')
        for c, v in config_vals:
            ax.annotate(c, (v, data.mean_rmse(c)), textcoords='offset points',
                        xytext=(8, 5), fontsize=9)
        ax.set_xlabel(param_name)
        ax.set_ylabel('Mean RMSE')
        if baseline_rmse:
            ax.axhline(y=baseline_rmse, color='red', linestyle='--', alpha=0.5, label='Baseline')
            ax.legend(fontsize=8)

    # Hide unused axes
    for idx in range(n_params, len(axes)):
        axes[idx].set_visible(False)

    plt.suptitle('Parameter Sensitivity', fontsize=13)
    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, bbox_inches='tight')
    return fig


def plot_trajectories(
    data: CalibrationData,
    configs: Optional[List[str]] = None,
    npz_dir: Optional[str | Path] = None,
    regions: List[str] = ('AK-PWS', 'AK-FN', 'BC-N', 'OR', 'CA-N'),
    seed: int = 42,
    region_colors: Optional[Dict[str, str]] = None,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Population trajectories for selected configs across key regions.

    Args:
        data: CalibrationData instance
        configs: configs to plot (default: top 5 by RMSE)
        npz_dir: directory containing config/monthly_seed*.npz files
        regions: regions to plot
        seed: which seed to use for trajectories
        region_colors: optional color mapping for regions
        save_path: optional path to save figure
    """
    apply_pub_style()
    if configs is None:
        configs = data.sorted_configs()[:5]
    if region_colors is None:
        cmap = plt.cm.tab10
        region_colors = {reg: cmap(i) for i, reg in enumerate(regions)}

    n = len(configs)
    fig, axes = plt.subplots(n, 1, figsize=(12, 3.5 * n), sharex=True)
    if n == 1:
        axes = [axes]

    for idx, config in enumerate(configs):
        ax = axes[idx]
        npz_path = Path(npz_dir) / config / f'monthly_seed{seed}.npz' if npz_dir else None

        if npz_path is None or not npz_path.exists():
            ax.text(0.5, 0.5, f'{config}: NPZ not available',
                    transform=ax.transAxes, ha='center')
            continue

        regional = load_npz_regional(npz_path)
        years = regional['_years']

        for reg in regions:
            pop = regional.get(reg, np.zeros_like(years))
            if pop[0] > 0:
                pop_norm = pop / pop[0]
            else:
                pop_norm = pop
            ax.plot(years, pop_norm, color=region_colors[reg], label=reg, linewidth=1.5)

        rmse = data.mean_rmse(config)
        desc = data.descriptions.get(config, '')
        title = f'{config} — {desc} (RMSE {rmse:.3f})' if desc else f'{config} (RMSE {rmse:.3f})'
        ax.set_title(title, fontsize=11)
        ax.set_ylabel('Relative\nPopulation')
        ax.legend(loc='upper right', fontsize=8, ncol=len(regions))
        ax.set_ylim(-0.05, 1.1)
        ax.axhline(y=0, color='gray', linewidth=0.5)

    axes[-1].set_xlabel('Year')
    plt.suptitle(f'Population Trajectories (seed {seed})', fontsize=13)
    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, bbox_inches='tight')
    return fig


def plot_seed_variability(
    data: CalibrationData,
    baseline_rmse: Optional[float] = None,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Box plot of RMSE spread across seeds per config.

    Args:
        data: CalibrationData instance
        baseline_rmse: optional reference line
        save_path: optional path to save figure
    """
    apply_pub_style()
    configs = data.sorted_configs()

    fig, ax = plt.subplots(figsize=(max(8, len(configs) * 0.8), 5))

    box_data = []
    box_colors = []
    for c in configs:
        vals = [v for v in data.rmse_per_seed(c) if np.isfinite(v)]
        box_data.append(vals if vals else [0])
        box_colors.append(data.color(c))

    bp = ax.boxplot(box_data, tick_labels=[c for c in configs],
                    patch_artist=True, widths=0.6)
    for patch, color in zip(bp['boxes'], box_colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    if baseline_rmse:
        ax.axhline(y=baseline_rmse, color='red', linestyle='--', linewidth=1.5,
                   label=f'Baseline ({baseline_rmse:.3f})')
        ax.legend()

    ax.set_ylabel('RMSE')
    ax.set_title('Seed Variability — RMSE Spread')
    plt.xticks(rotation=45, ha='right')

    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, bbox_inches='tight')
    return fig


def plot_config_vs_targets(
    data: CalibrationData,
    config: str,
    targets: Optional[Dict[str, float]] = None,
    regions: Optional[List[str]] = None,
    config_color: str = '#2196F3',
    target_color: str = '#4CAF50',
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Side-by-side bar chart comparing a config's regional recovery to targets.

    Args:
        data: CalibrationData instance
        config: config name (e.g. 'W257')
        targets: recovery targets dict (default: DEFAULT_TARGETS)
        regions: list of regions to show (default: SCORED_REGIONS)
        config_color: color for actual recovery bars
        target_color: color for target bars
        save_path: optional path to save figure
    """
    apply_pub_style()
    if targets is None:
        targets = DEFAULT_TARGETS
    if regions is None:
        regions = SCORED_REGIONS

    x = np.arange(len(regions))
    width = 0.35

    target_vals = [targets.get(reg, 0) * 100 for reg in regions]
    actual_vals = [data.recovery(config, reg) * 100 for reg in regions]

    fig, ax = plt.subplots(figsize=(10, 5))
    bars_target = ax.bar(x - width / 2, target_vals, width, label='Target',
                         color=target_color, edgecolor='white', alpha=0.85)
    bars_actual = ax.bar(x + width / 2, actual_vals, width,
                         label=f'{config} actual', color=config_color,
                         edgecolor='white', alpha=0.85)

    # Value labels on bars
    for bar, val in zip(bars_target, target_vals):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.3,
                f'{val:.1f}%', ha='center', va='bottom', fontsize=8, color=target_color)
    for bar, val in zip(bars_actual, actual_vals):
        text = f'{val:.1f}%' if val >= 1 else f'{val:.2f}%'
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.3,
                text, ha='center', va='bottom', fontsize=8, color=config_color)

    ax.set_xticks(x)
    ax.set_xticklabels(regions, rotation=45, ha='right')
    ax.set_ylabel('Recovery (%)')
    desc = data.descriptions.get(config, '')
    rmse = data.mean_rmse(config)
    title = f'{config} vs Calibration Targets'
    if desc:
        title += f' — {desc}'
    title += f' (RMSE={rmse:.3f})'
    ax.set_title(title)
    ax.legend()
    ax.set_ylim(0, max(max(target_vals), max(actual_vals)) * 1.2)

    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, bbox_inches='tight')
    return fig


# ══════════════════════════════════════════════════════════════════════
# CONVENIENCE: Generate all standard calibration figures
# ══════════════════════════════════════════════════════════════════════

def generate_standard_figures(
    data: CalibrationData,
    output_dir: str | Path,
    npz_dir: Optional[str | Path] = None,
    baseline_rmse: Optional[float] = None,
    baseline_label: str = 'Baseline',
    param_configs: Optional[Dict] = None,
) -> List[Path]:
    """Generate all standard calibration report figures.

    Returns list of generated file paths.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    generated = []

    print("Generating standard calibration figures...")

    # 1. RMSE ranking
    fig = plot_rmse_ranking(data, baseline_rmse=baseline_rmse,
                            baseline_label=baseline_label,
                            save_path=output_dir / 'fig_rmse_ranking.pdf')
    plt.close(fig)
    generated.append(output_dir / 'fig_rmse_ranking.pdf')
    print(f"  ✓ fig_rmse_ranking.pdf")

    # 2. Region heatmap
    fig = plot_region_heatmap(data, save_path=output_dir / 'fig_region_heatmap.pdf')
    plt.close(fig)
    generated.append(output_dir / 'fig_region_heatmap.pdf')
    print(f"  ✓ fig_region_heatmap.pdf")

    # 3. Tradeoff scatter
    fig = plot_tradeoff(data, save_path=output_dir / 'fig_tradeoff.pdf')
    plt.close(fig)
    generated.append(output_dir / 'fig_tradeoff.pdf')
    print(f"  ✓ fig_tradeoff.pdf")

    # 4. Parameter effects (if provided)
    if param_configs:
        fig = plot_parameter_effects(data, param_configs, baseline_rmse=baseline_rmse,
                                     save_path=output_dir / 'fig_parameter_effects.pdf')
        plt.close(fig)
        generated.append(output_dir / 'fig_parameter_effects.pdf')
        print(f"  ✓ fig_parameter_effects.pdf")

    # 5. Trajectories (if NPZ available)
    if npz_dir:
        fig = plot_trajectories(data, npz_dir=npz_dir,
                                save_path=output_dir / 'fig_trajectories.pdf')
        plt.close(fig)
        generated.append(output_dir / 'fig_trajectories.pdf')
        print(f"  ✓ fig_trajectories.pdf")

    # 6. Seed variability
    fig = plot_seed_variability(data, baseline_rmse=baseline_rmse,
                                save_path=output_dir / 'fig_seed_variability.pdf')
    plt.close(fig)
    generated.append(output_dir / 'fig_seed_variability.pdf')
    print(f"  ✓ fig_seed_variability.pdf")

    # 7. Export data summary
    data.to_summary_json(output_dir / 'data_summary.json')
    generated.append(output_dir / 'data_summary.json')

    print(f"\nGenerated {len(generated)} files in {output_dir}")
    return generated
