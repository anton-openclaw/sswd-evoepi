#!/usr/bin/env python3
"""Plot scaling axes for SSWD-EvoEpi."""
import json, sys
sys.path.insert(0, '.')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# Load data
with open('results/performance/scaling_axes_data.json') as f:
    data = json.load(f)

# Dark theme
BG = '#1a1a2e'
FG = '#e6e6e6'
PANEL = '#16213e'
ACCENT1 = '#e94560'
ACCENT2 = '#0f3460'
ACCENT3 = '#533483'
ACCENT4 = '#00d2ff'
GRID = '#2a2a4e'

fig = plt.figure(figsize=(16, 12), facecolor=BG)
gs = GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.35)

def style_ax(ax, title, xlabel, ylabel):
    ax.set_facecolor(PANEL)
    ax.set_title(title, color=FG, fontsize=13, fontweight='bold', pad=10)
    ax.set_xlabel(xlabel, color=FG, fontsize=11)
    ax.set_ylabel(ylabel, color=FG, fontsize=11)
    ax.tick_params(colors=FG, labelsize=9)
    ax.grid(True, alpha=0.2, color=GRID)
    for spine in ax.spines.values():
        spine.set_color('#444')

# ── Panel 1: N scaling (linear) ──────────────────────────────────
ax1 = fig.add_subplot(gs[0, 0])
n_data = [(r['n'], r['wall_s']) for r in data['axis_N'] if r.get('ok')]
ns, ts = zip(*n_data)
ax1.plot(ns, ts, 'o-', color=ACCENT1, linewidth=2.5, markersize=8, zorder=5)

# Fit power law
log_ns, log_ts = np.log(ns), np.log(ts)
slope, intercept = np.polyfit(log_ns, log_ts, 1)
n_fit = np.linspace(min(ns), max(ns), 100)
t_fit = np.exp(intercept) * n_fit ** slope
ax1.plot(n_fit, t_fit, '--', color=ACCENT4, alpha=0.6, linewidth=1.5, label=f'O(N^{slope:.2f})')
ax1.legend(fontsize=10, facecolor=PANEL, edgecolor='#555', labelcolor=FG)
style_ax(ax1, 'Population Size Scaling', 'Agents per Node (N)', 'Wall Time (s)')

# ── Panel 2: N scaling (log-log) ─────────────────────────────────
ax2 = fig.add_subplot(gs[0, 1])
ax2.loglog(ns, ts, 'o-', color=ACCENT1, linewidth=2.5, markersize=8, zorder=5)
ax2.loglog(n_fit, t_fit, '--', color=ACCENT4, alpha=0.6, linewidth=1.5, label=f'O(N^{slope:.2f})')
ax2.legend(fontsize=10, facecolor=PANEL, edgecolor='#555', labelcolor=FG)
style_ax(ax2, 'Population Size (log-log)', 'Agents per Node (N)', 'Wall Time (s)')

# ── Panel 3: T scaling ───────────────────────────────────────────
ax3 = fig.add_subplot(gs[0, 2])
t_data = [(r['years'], r['wall_s']) for r in data['axis_T'] if r.get('ok')]
years, times = zip(*t_data)
ax3.plot(years, times, 's-', color=ACCENT3, linewidth=2.5, markersize=8, zorder=5)

# Fit
log_y, log_t = np.log(years), np.log(times)
slope_t, int_t = np.polyfit(log_y, log_t, 1)
y_fit = np.linspace(min(years), max(years), 100)
t_fit2 = np.exp(int_t) * y_fit ** slope_t
ax3.plot(y_fit, t_fit2, '--', color=ACCENT4, alpha=0.6, linewidth=1.5, label=f'O(T^{slope_t:.2f})')
ax3.legend(fontsize=10, facecolor=PANEL, edgecolor='#555', labelcolor=FG)
style_ax(ax3, 'Duration Scaling (200 agents)', 'Simulation Years (T)', 'Wall Time (s)')

# ── Panel 4: Spawning overhead ────────────────────────────────────
ax4 = fig.add_subplot(gs[1, 0])
sp_data = data['axis_spawning']
sp_ns = [r['n'] for r in sp_data]
sp_off = [r['no_spawn'] for r in sp_data]
sp_on = [r['spawn'] for r in sp_data]
x = np.arange(len(sp_ns))
w = 0.35
ax4.bar(x - w/2, sp_off, w, label='No Spawning', color=ACCENT2, alpha=0.85)
ax4.bar(x + w/2, sp_on, w, label='Spawning', color=ACCENT1, alpha=0.85)
ax4.set_xticks(x)
ax4.set_xticklabels(sp_ns)
ax4.legend(fontsize=9, facecolor=PANEL, edgecolor='#555', labelcolor=FG)
style_ax(ax4, 'Spawning Overhead (5yr)', 'Agents per Node', 'Wall Time (s)')

# ── Panel 5: Disease overhead ─────────────────────────────────────
ax5 = fig.add_subplot(gs[1, 1])
ds_data = data['axis_disease']
ds_ns = [r['n'] for r in ds_data]
ds_off = [r['no_disease'] for r in ds_data]
ds_on = [r['disease'] for r in ds_data]
x = np.arange(len(ds_ns))
ax5.bar(x - w/2, ds_off, w, label='No Disease', color=ACCENT2, alpha=0.85)
ax5.bar(x + w/2, ds_on, w, label='Disease (yr 3)', color=ACCENT3, alpha=0.85)
ax5.set_xticks(x)
ax5.set_xticklabels(ds_ns)
ax5.legend(fontsize=9, facecolor=PANEL, edgecolor='#555', labelcolor=FG)
style_ax(ax5, 'Disease Overhead (10yr)', 'Agents per Node', 'Wall Time (s)')

# ── Panel 6: Projections table ────────────────────────────────────
ax6 = fig.add_subplot(gs[1, 2])
ax6.set_facecolor(PANEL)
ax6.axis('off')

# Calculate projections
ref_n, ref_t = 2000, 3.185  # Largest measured point, 5yr
projections = []
for n, t in [(500, 20), (1000, 20), (2000, 20), (500, 50), (1000, 50), (2000, 50),
             (500, 100), (1000, 100)]:
    est = ref_t * (n / 2000) ** slope * (t / 5) ** slope_t
    projections.append((n, t, est))

# Table
col_labels = ['N', 'Years', 'Est. Time']
cell_text = []
for n, t, est in projections:
    if est < 60:
        tstr = f'{est:.1f}s'
    elif est < 3600:
        tstr = f'{est/60:.1f}min'
    else:
        tstr = f'{est/3600:.1f}hr'
    cell_text.append([str(n), str(t), tstr])

table = ax6.table(cellText=cell_text, colLabels=col_labels,
                   loc='center', cellLoc='center')
table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1, 1.5)

# Style table
for (row, col), cell in table.get_celld().items():
    cell.set_edgecolor('#444')
    if row == 0:
        cell.set_facecolor(ACCENT2)
        cell.set_text_props(color=FG, fontweight='bold')
    else:
        cell.set_facecolor(PANEL)
        cell.set_text_props(color=FG)

ax6.set_title('Time Projections', color=FG, fontsize=13, fontweight='bold', pad=10)

# Title
fig.suptitle('SSWD-EvoEpi: Scaling Analysis Along Major Axes',
             color=FG, fontsize=16, fontweight='bold', y=0.98)
fig.text(0.5, 0.94, f'N scaling: O(N^{slope:.2f})  |  T scaling: O(T^{slope_t:.2f})  |  Spawning: 1-5× overhead  |  Disease: net speedup (kills reduce N)',
         ha='center', color='#888', fontsize=10)

plt.savefig('results/performance/scaling_axes.png', dpi=150, facecolor=BG, bbox_inches='tight')
print('Saved results/performance/scaling_axes.png')
