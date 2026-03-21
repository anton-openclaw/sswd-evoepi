import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import matplotlib.colors as mcolors

plt.style.use("seaborn-v0_8-whitegrid")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), gridspec_kw={"width_ratios": [1.2, 1]})

# ── Left panel: Resistance distributions ──────────────────────────────────
r = np.linspace(-0.1, 0.7, 500)

# Pre-disease
mu_pre, sd_pre = 0.15, 0.1612
y_pre = norm.pdf(r, mu_pre, sd_pre)

# Post-selection
mu_post, sd_post = 0.2119, 0.1454
y_post = norm.pdf(r, mu_post, sd_post)

ax1.fill_between(r, y_pre, alpha=0.35, color="#b0b0b0", label="Pre-disease (2012)")
ax1.plot(r, y_pre, color="#707070", lw=1.5)

ax1.fill_between(r, y_post, alpha=0.40, color="#6baed6", label="Post-selection (2024)")
ax1.plot(r, y_post, color="#2171b5", lw=2)

# Set limits first
ax1.set_xlim(0, 0.65)
ymax = max(y_pre.max(), y_post.max()) * 1.12
ax1.set_ylim(0, ymax)

# Threshold lines + labels (once)
for thresh, lbl in [(0.30, "r = 0.30"), (0.40, "r = 0.40"), (0.50, "r = 0.50")]:
    ax1.axvline(thresh, ls="--", color="#d62728", lw=1.2, alpha=0.7)
    ax1.text(thresh, ymax * 0.98, lbl, fontsize=8, color="#d62728",
             ha="left", va="top", fontweight="bold",
             bbox=dict(facecolor="white", edgecolor="none", alpha=0.7, pad=1))

ax1.set_title("Resistance Distribution in Washington State", fontsize=13, fontweight="bold", pad=10)
ax1.set_xlabel("Resistance (r) — probability of avoiding infection", fontsize=11)
ax1.set_ylabel("Probability density", fontsize=11)
ax1.legend(loc="upper left", fontsize=10, framealpha=0.9)

# ── Right panel: Sampling effort bar chart ────────────────────────────────
thresholds = ["r > 0.20", "r > 0.25", "r > 0.30", "r > 0.35",
              "r > 0.40", "r > 0.45", "r > 0.50"]
n_sample = [2, 3, 4, 6, 11, 20, 43]

# Color gradient: green → orange → red
cmap = mcolors.LinearSegmentedColormap.from_list("effort",
       ["#2ca02c", "#f0c929", "#e67e22", "#d62728"])
norm_c = plt.Normalize(vmin=min(n_sample), vmax=max(n_sample))
colors = [cmap(norm_c(n)) for n in n_sample]

y_pos = np.arange(len(thresholds))
bars = ax2.barh(y_pos, n_sample, color=colors, edgecolor="white", height=0.6)

ax2.set_yticks(y_pos)
ax2.set_yticklabels(thresholds, fontsize=10)
ax2.invert_yaxis()  # highest threshold at top

# Annotations
for i, (bar, n) in enumerate(zip(bars, n_sample)):
    ax2.text(bar.get_width() + 0.8, bar.get_y() + bar.get_height() / 2,
             str(n), va="center", ha="left", fontsize=10.5, fontweight="bold",
             color="#333333")

ax2.set_xlim(0, max(n_sample) * 1.18)
ax2.set_title("Sampling Effort for Breeding Program", fontsize=13, fontweight="bold", pad=10)
ax2.set_xlabel("Individuals to sample", fontsize=11)
ax2.set_ylabel("Resistance threshold", fontsize=11)

plt.tight_layout(w_pad=3)
out = "/home/starbot/.openclaw/workspace/sswd-evoepi/reports/progress_report/figures/fig_breeding_sampling.png"
plt.savefig(out, dpi=150, bbox_inches="tight", facecolor="white")
plt.close()

import os
size_kb = os.path.getsize(out) / 1024
print(f"Saved: {out}")
print(f"Size: {size_kb:.1f} KB")
