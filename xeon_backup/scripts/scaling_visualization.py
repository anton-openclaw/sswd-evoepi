import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import json

# Load scaling data
with open('results/performance/scaling_data.json') as f:
    data = json.load(f)

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
fig.patch.set_facecolor('#1a1a2e')
for ax in axes:
    ax.set_facecolor('#16213e')
    ax.tick_params(colors='white')
    ax.spines['bottom'].set_color('#555')
    ax.spines['left'].set_color('#555')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.xaxis.label.set_color('white')
    ax.yaxis.label.set_color('white')
    ax.title.set_color('white')

# Panel 1: Wall time by population size (5yr)
data_5yr = [d for d in data if d['years'] == 5]
ns = [d['n'] for d in data_5yr]
times = [d['wall_s'] for d in data_5yr]
axes[0].plot(ns, times, 'o-', color='#e94560', linewidth=2, markersize=8)
axes[0].set_xlabel('Population Size')
axes[0].set_ylabel('Wall Time (s)')
axes[0].set_title('Scaling: 5-Year Simulation')

# Panel 2: Component breakdown (stacked bar)
labels = [d['label'] for d in data]
components_all = {}
for d in data:
    for k, v in d['components'].items():
        if k != '_total_s':
            if k not in components_all:
                components_all[k] = []
            components_all[k].append(v.get('total_s', 0))

# Pad missing entries
max_len = len(data)
for k in components_all:
    while len(components_all[k]) < max_len:
        components_all[k].append(0)

# Sort by total time
sorted_comps = sorted(components_all.items(), key=lambda x: -sum(x[1]))
colors = ['#e94560', '#0f3460', '#533483', '#16213e', '#1a1a2e', '#e94560']
bottom = [0] * max_len
x_pos = range(max_len)
for i, (name, vals) in enumerate(sorted_comps[:5]):
    color = colors[i % len(colors)]
    axes[1].bar(x_pos, vals, bottom=bottom, label=name, color=color, alpha=0.8)
    bottom = [b + v for b, v in zip(bottom, vals)]
axes[1].set_xticks(list(x_pos))
axes[1].set_xticklabels([d['label'].replace(' agents, ', '\n') for d in data], rotation=45, ha='right', fontsize=8, color='white')
axes[1].set_ylabel('Time (s)')
axes[1].set_title('Component Breakdown')
axes[1].legend(fontsize=8, facecolor='#16213e', edgecolor='#555', labelcolor='white')

plt.tight_layout()
plt.savefig('results/performance/scaling_analysis.png', dpi=150, facecolor='#1a1a2e')
print('Saved results/performance/scaling_analysis.png')