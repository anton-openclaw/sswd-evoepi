#!/usr/bin/env python3
"""Generate figures + LaTeX report for W95-W104 (pathogen thermal adaptation)."""
import json, os, math, subprocess
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = "/home/starbot/.openclaw/workspace/sswd-evoepi/reports/pathogen_adapt_w95_w104"
FIGS = os.path.join(OUT, "figures")
os.makedirs(FIGS, exist_ok=True)

# ── Load data from Xeon via pre-extracted file ──
regions = ["AK-PWS", "AK-FN", "AK-FS", "BC-N", "SS-S", "JDF", "OR", "CA-N"]
targets = {"AK-PWS": 50, "AK-FN": 50, "AK-FS": 20, "BC-N": 20, "SS-S": 5, "JDF": 2, "OR": 0.25, "CA-N": 0.1}
years = np.arange(2012, 2025)

region_colors = {
    "AK-PWS": "#1f77b4", "AK-FN": "#4a90d9", "AK-FS": "#7eb3e8",
    "BC-N": "#2ca02c", "SS-S": "#ff7f0e", "JDF": "#d62728",
    "OR": "#9467bd", "CA-N": "#e377c2"
}

# Parse the raw data
all_data = {}
with open("/tmp/w95_w104_raw.txt") as f:
    current = None
    for line in f:
        line = line.strip()
        if line.startswith("===W") and line.endswith("==="):
            current = line.replace("=", "")
            all_data[current] = {"regions": {}, "T_vbnc": {}}
        elif current and "|" in line and not line.startswith("wall"):
            parts = line.split("|")
            reg = parts[0]
            pct = float(parts[1])
            tvbnc = float(parts[2])
            yt = eval(parts[3]) if parts[3].strip() else []
            all_data[current]["regions"][reg] = {
                "pct": pct, "yearly_totals": yt
            }
            all_data[current]["T_vbnc"][reg] = tvbnc

# Run configs
configs = {
    "W90": {"floor": 300, "alpha": 0.20, "adapt": 0, "tmin": "-"},
    "W91": {"floor": 500, "alpha": 0.20, "adapt": 0, "tmin": "-"},
    "W95": {"floor": 300, "alpha": 0.20, "adapt": 0.0005, "tmin": 6},
    "W96": {"floor": 300, "alpha": 0.20, "adapt": 0.001, "tmin": 6},
    "W97": {"floor": 300, "alpha": 0.20, "adapt": 0.002, "tmin": 6},
    "W98": {"floor": 300, "alpha": 0.20, "adapt": 0.005, "tmin": 6},
    "W99": {"floor": 500, "alpha": 0.20, "adapt": 0.001, "tmin": 6},
    "W100": {"floor": 500, "alpha": 0.20, "adapt": 0.002, "tmin": 6},
    "W101": {"floor": 300, "alpha": 0.20, "adapt": 0.001, "tmin": 8},
    "W102": {"floor": 300, "alpha": 0.20, "adapt": 0.001, "tmin": 4},
    "W103": {"floor": 300, "alpha": 0.15, "adapt": 0.001, "tmin": 6},
    "W104": {"floor": 200, "alpha": 0.15, "adapt": 0.001, "tmin": 6},
}

def recovery_frac(yt):
    if not yt or len(yt) < 3:
        return []
    peak = max(yt[:3])
    if peak == 0:
        return []
    return [p / peak for p in yt[:13]]

# ══ FIG 1: Key runs 2×2 ══
fig, axes = plt.subplots(2, 2, figsize=(14, 10), sharex=True, sharey=True)
fig.suptitle("Recovery Time Series: Pathogen Adaptation Key Runs", fontsize=14, fontweight='bold')
key = [("W90", "W90 (no adapt, reference)"), ("W96", "W96 (adapt=0.001)"),
       ("W99", "W99 (floor=500, adapt=0.001) BEST"), ("W98", "W98 (adapt=0.005, highest)")]
for ax, (wn, title) in zip(axes.flat, key):
    d = all_data[wn]
    for reg in regions:
        if reg not in d["regions"]:
            continue
        frac = recovery_frac(d["regions"][reg]["yearly_totals"])
        if frac:
            ax.plot(years[:len(frac)], frac, '-o', ms=3, lw=1.5, color=region_colors[reg], label=reg)
    ax.set_title(title, fontsize=10)
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, alpha=0.3)
    ax.set_ylabel("Pop / Peak")
    ax.set_xlabel("Year")
handles, labels = axes[0,0].get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center', ncol=4, fontsize=9, bbox_to_anchor=(0.5, -0.02))
plt.tight_layout(rect=[0, 0.04, 1, 0.96])
plt.savefig(f"{FIGS}/fig1_key_runs.png", dpi=150, bbox_inches='tight')
plt.close()
print("Fig 1")

# ══ FIG 2: All runs grid ══
all_runs = ["W90", "W91"] + [f"W{i}" for i in range(95, 105)]
fig, axes = plt.subplots(3, 4, figsize=(18, 12), sharex=True, sharey=True)
fig.suptitle("All Runs: Recovery Trajectories", fontsize=14, fontweight='bold')
for idx, wn in enumerate(all_runs):
    row, col = divmod(idx, 4)
    ax = axes[row, col]
    d = all_data.get(wn, {})
    for reg in regions:
        if reg not in d.get("regions", {}):
            continue
        frac = recovery_frac(d["regions"][reg]["yearly_totals"])
        if frac:
            ax.plot(years[:len(frac)], frac, '-o', ms=2, lw=1, color=region_colors[reg],
                    label=reg if idx == 0 else None)
    c = configs.get(wn, {})
    ax.set_title(f"{wn} f={c.get('floor','')} α={c.get('alpha','')} r={c.get('adapt','')}", fontsize=8)
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, alpha=0.3)
handles, labels = axes[0,0].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper right', fontsize=8)
plt.tight_layout()
plt.savefig(f"{FIGS}/fig2_all_runs.png", dpi=150, bbox_inches='tight')
plt.close()
print("Fig 2")

# ══ FIG 3: Adapt rate effect (W95-W98) ══
fig, axes = plt.subplots(1, 4, figsize=(18, 5), sharey=True)
fig.suptitle("Adaptation Rate Effect (floor=300, α=0.20)", fontsize=13, fontweight='bold')
for ax, wn in zip(axes, ["W95", "W96", "W97", "W98"]):
    d = all_data[wn]
    for reg in ["AK-PWS", "BC-N", "SS-S", "JDF", "CA-N"]:
        frac = recovery_frac(d["regions"][reg]["yearly_totals"])
        if frac:
            ax.plot(years[:len(frac)], frac, '-o', ms=3, lw=1.5, color=region_colors[reg], label=reg)
            ax.axhline(y=targets[reg]/100, color=region_colors[reg], ls='--', alpha=0.3)
    c = configs[wn]
    t_ak = d["T_vbnc"].get("AK-PWS", 12)
    ax.set_title(f"{wn}: rate={c['adapt']}\nAK T_vbnc={t_ak:.2f}°C", fontsize=9)
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlabel("Year")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7, loc='upper right')
axes[0].set_ylabel("Pop / Peak")
plt.tight_layout()
plt.savefig(f"{FIGS}/fig3_adapt_rate.png", dpi=150, bbox_inches='tight')
plt.close()
print("Fig 3")

# ══ FIG 4: Bar chart W90 vs W99 vs targets ══
fig, ax = plt.subplots(figsize=(12, 6))
x = np.arange(len(regions))
width = 0.25
tgt = [targets[r] for r in regions]
w90 = [all_data["W90"]["regions"][r]["pct"] for r in regions]
w99 = [all_data["W99"]["regions"][r]["pct"] for r in regions]
ax.bar(x - width, tgt, width, label='Target', color='#2ecc71', alpha=0.8, edgecolor='black', lw=0.5)
ax.bar(x, w90, width, label='W90 (no adapt)', color='#3498db', alpha=0.8, edgecolor='black', lw=0.5)
ax.bar(x + width, w99, width, label='W99 (adapt, BEST)', color='#e74c3c', alpha=0.8, edgecolor='black', lw=0.5)
ax.set_xticks(x)
ax.set_xticklabels(regions)
ax.set_ylabel("Recovery %")
ax.set_title("Recovery: W90 (no adapt) vs W99 (adapt, RMSE=0.557) vs Targets", fontsize=13, fontweight='bold')
ax.legend()
ax.set_yscale('symlog', linthresh=1)
ax.set_ylim(0, 60)
ax.grid(True, alpha=0.3, axis='y')
plt.tight_layout()
plt.savefig(f"{FIGS}/fig4_bar_comparison.png", dpi=150, bbox_inches='tight')
plt.close()
print("Fig 4")

# ══ FIG 5: T_vbnc by region — THE key figure ══
fig, ax = plt.subplots(figsize=(14, 6))
adapt_runs = [f"W{i}" for i in range(95, 105)]
x = np.arange(len(regions))
width = 0.08
for i, wn in enumerate(adapt_runs):
    d = all_data[wn]
    vals = [d["T_vbnc"].get(r, 12.0) for r in regions]
    offset = (i - len(adapt_runs)/2 + 0.5) * width
    ax.bar(x + offset, vals, width, label=wn, alpha=0.8)
ax.axhline(y=12.0, color='red', ls='--', lw=2, label='Initial T_vbnc (12°C)')
ax.set_xticks(x)
ax.set_xticklabels(regions)
ax.set_ylabel("Final T_vbnc (°C)")
ax.set_title("Pathogen Thermal Adaptation: Final T_vbnc by Region\n(barely moves from 12.0°C)", fontsize=13, fontweight='bold')
ax.set_ylim(11.8, 12.05)
ax.legend(fontsize=7, ncol=5, loc='lower left')
ax.grid(True, alpha=0.3, axis='y')
plt.tight_layout()
plt.savefig(f"{FIGS}/fig5_tvbnc_adaptation.png", dpi=150, bbox_inches='tight')
plt.close()
print("Fig 5")

# ══ FIG 6: RMSE comparison ══
fig, ax = plt.subplots(figsize=(12, 5))
all_names = ["W90", "W91"] + [f"W{i}" for i in range(95, 105)]
rmses = []
for wn in all_names:
    d = all_data.get(wn, {})
    sq = []
    for r in regions:
        if r in d.get("regions", {}):
            actual = d["regions"][r]["pct"] / 100
            target = targets[r] / 100
            if actual > 0 and target > 0:
                sq.append((math.log10(actual) - math.log10(target))**2)
    rmses.append(math.sqrt(sum(sq)/len(sq)) if sq else 99)
colors = ['#3498db' if wn in ['W90','W91'] else '#e74c3c' for wn in all_names]
ax.bar(range(len(all_names)), rmses, color=colors, alpha=0.8, edgecolor='black', lw=0.5)
ax.set_xticks(range(len(all_names)))
ax.set_xticklabels(all_names, rotation=45, fontsize=9)
ax.set_ylabel("RMSE (log-space)")
ax.set_title("RMSE Comparison: Reference (blue) vs Adaptation Runs (red)", fontsize=13, fontweight='bold')
ax.axhline(y=min(rmses), color='green', ls='--', alpha=0.5, label=f'Best: {min(rmses):.3f}')
ax.legend()
ax.grid(True, alpha=0.3, axis='y')
plt.tight_layout()
plt.savefig(f"{FIGS}/fig6_rmse.png", dpi=150, bbox_inches='tight')
plt.close()
print("Fig 6")

# ══ FIG 7: T_vbnc_min irrelevance (W96=W101=W102) ══
fig, axes = plt.subplots(1, 3, figsize=(16, 5), sharey=True)
fig.suptitle("T_vbnc_min Irrelevance: W96 (Tmin=6) = W101 (Tmin=8) = W102 (Tmin=4)", fontsize=12, fontweight='bold')
for ax, wn in zip(axes, ["W96", "W101", "W102"]):
    d = all_data[wn]
    for reg in ["AK-PWS", "BC-N", "SS-S", "JDF", "CA-N"]:
        frac = recovery_frac(d["regions"][reg]["yearly_totals"])
        if frac:
            ax.plot(years[:len(frac)], frac, '-o', ms=3, lw=1.5, color=region_colors[reg], label=reg)
    c = configs[wn]
    ax.set_title(f"{wn}: T_vbnc_min={c['tmin']}°C", fontsize=10)
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlabel("Year")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7)
axes[0].set_ylabel("Pop / Peak")
plt.tight_layout()
plt.savefig(f"{FIGS}/fig7_tmin_irrelevance.png", dpi=150, bbox_inches='tight')
plt.close()
print("Fig 7")

# ══ FIG 8: Alpha effect (W96 vs W103 vs W104) ══
fig, axes = plt.subplots(1, 3, figsize=(16, 5), sharey=True)
fig.suptitle("α_env Effect with Adaptation: 0.20 vs 0.15", fontsize=13, fontweight='bold')
for ax, wn in zip(axes, ["W96", "W103", "W104"]):
    d = all_data[wn]
    for reg in regions:
        frac = recovery_frac(d["regions"][reg]["yearly_totals"])
        if frac:
            ax.plot(years[:len(frac)], frac, '-o', ms=3, lw=1.5, color=region_colors[reg], label=reg)
    c = configs[wn]
    ax.set_title(f"{wn}: f={c['floor']} α={c['alpha']}", fontsize=10)
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlabel("Year")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7, loc='upper right')
axes[0].set_ylabel("Pop / Peak")
plt.tight_layout()
plt.savefig(f"{FIGS}/fig8_alpha_effect.png", dpi=150, bbox_inches='tight')
plt.close()
print("Fig 8")

print("\nAll figures generated!")

# ══ Generate LaTeX ══
tex = r"""\documentclass[11pt]{article}
\usepackage[margin=1in]{geometry}
\usepackage{graphicx}
\usepackage{booktabs}
\usepackage{caption}
\usepackage{hyperref}
\usepackage{fancyhdr}
\pagestyle{fancy}
\lhead{W95--W104 Report}
\rhead{SSWD Calibration}

\title{Calibration Report: Pathogen Thermal Adaptation\\(W95--W104)}
\author{Anton Star --- SSWD-EvoEpi Project}
\date{March 5, 2026}

\begin{document}
\maketitle

\section{Executive Summary}
W95--W104 are the first calibration runs with \textbf{pathogen thermal adaptation} enabled --- a mechanism allowing each node's VBNC threshold ($T_\mathrm{vbnc}$) to evolve downward as cold-adapted pathogen strains are selected.

\textbf{Key finding:} Thermal adaptation is negligible at tested rates. Even at the highest adaptation rate (W98, $r=0.005$), $T_\mathrm{vbnc}$ shifted by only 0.14\textdegree C at AK-PWS (from 12.0 to 11.86). JDF and SS-S showed $<$0.02\textdegree C change. The mechanism does not operate on the 13-year simulation timescale.

\textbf{Best run: W99} (floor=500, $\alpha$=0.20, adapt=0.001) achieves RMSE=0.557 --- our best ever. However, this is driven by the floor/$\alpha$ combination, not adaptation. AK-PWS remains at 10.8\% (target 50\%).

\section{Methods}
All runs share: $K_{1/2}$=800K, CDT=1000, wfDP=300, wfRange=3000, dynamic $P_\mathrm{env}$, $\delta$=0.02, revert\_rate=0.0005.

\begin{table}[h]
\centering
\small
\begin{tabular}{lcccc}
\toprule
Run & Floor & $\alpha$ & Adapt Rate & $T_\mathrm{vbnc,min}$ \\
\midrule
W90 (ref) & 300 & 0.20 & --- & --- \\
W91 (ref) & 500 & 0.20 & --- & --- \\
\midrule
W95 & 300 & 0.20 & 0.0005 & 6 \\
W96 & 300 & 0.20 & 0.001 & 6 \\
W97 & 300 & 0.20 & 0.002 & 6 \\
W98 & 300 & 0.20 & 0.005 & 6 \\
W99 & 500 & 0.20 & 0.001 & 6 \\
W100 & 500 & 0.20 & 0.002 & 6 \\
W101 & 300 & 0.20 & 0.001 & 8 \\
W102 & 300 & 0.20 & 0.001 & 4 \\
W103 & 300 & 0.15 & 0.001 & 6 \\
W104 & 200 & 0.15 & 0.001 & 6 \\
\bottomrule
\end{tabular}
\end{table}

\section{Results}

\begin{table}[h]
\centering
\small
\begin{tabular}{lccccccccr}
\toprule
Run & AK-PWS & AK-FN & AK-FS & BC-N & SS-S & JDF & OR & CA-N & RMSE \\
 & (50) & (50) & (20) & (20) & (5) & (2) & (0.25) & (0.1) & \\
\midrule
W90 & 10.5 & 4.8 & 6.0 & 6.9 & 5.8 & 1.9 & 0.9 & 1.0 & 0.643 \\
W91 & 11.3 & 5.7 & 2.6 & 9.3 & 7.4 & 5.9 & 0.5 & 0.3 & 0.591 \\
\midrule
W95 & 10.6 & 4.6 & 4.2 & 8.4 & 9.3 & 19.7 & 1.0 & 1.6 & 0.790 \\
W96 & 11.3 & 5.3 & 6.4 & 8.5 & 4.0 & 9.2 & 1.4 & 0.6 & 0.649 \\
W97 & 9.9 & 7.0 & 12.0 & 7.7 & 5.4 & 2.0 & 3.0 & 0.8 & 0.654 \\
W98 & 12.1 & 6.0 & 3.6 & 8.6 & 4.8 & 15.4 & 1.4 & 0.8 & 0.709 \\
\textbf{W99} & \textbf{10.8} & \textbf{8.8} & \textbf{6.2} & \textbf{8.4} & \textbf{2.9} & \textbf{4.3} & \textbf{1.0} & \textbf{0.6} & \textbf{0.557} \\
W100 & 10.1 & 4.1 & 3.7 & 8.1 & 4.4 & 3.3 & 2.3 & 0.5 & 0.687 \\
W101 & 11.3 & 5.3 & 6.4 & 8.5 & 4.0 & 9.2 & 1.4 & 0.6 & 0.649 \\
W102 & 11.3 & 5.3 & 6.4 & 8.5 & 4.0 & 9.2 & 1.4 & 0.6 & 0.649 \\
W103 & 20.8 & 9.4 & 8.0 & 13.6 & 6.8 & 19.6 & 6.7 & 4.9 & 0.919 \\
W104 & 24.9 & 6.1 & 10.2 & 14.0 & 15.2 & 28.5 & 6.0 & 5.6 & 0.972 \\
\bottomrule
\end{tabular}
\caption{Recovery fractions (\%) by region. Targets in parentheses. W99 bolded as best RMSE.}
\end{table}

\subsection{Key Runs: Recovery Time Series}
\begin{figure}[h]
\centering
\includegraphics[width=\textwidth]{figures/fig1_key_runs.png}
\caption{Recovery trajectories for W90 (reference), W96 (adapt=0.001), W99 (best RMSE), W98 (highest adapt rate).}
\end{figure}

\subsection{All Runs Grid}
\begin{figure}[h]
\centering
\includegraphics[width=\textwidth]{figures/fig2_all_runs.png}
\caption{All 12 runs (2 references + 10 adaptation runs).}
\end{figure}

\subsection{Adaptation Rate Effect}
\begin{figure}[h]
\centering
\includegraphics[width=\textwidth]{figures/fig3_adapt_rate.png}
\caption{Increasing adaptation rate (W95--W98). Note final $T_\mathrm{vbnc}$ values barely change.}
\end{figure}

\subsection{Recovery Comparison}
\begin{figure}[h]
\centering
\includegraphics[width=\textwidth]{figures/fig4_bar_comparison.png}
\caption{W90 vs W99 vs targets. AK-PWS remains $\sim$5$\times$ below target.}
\end{figure}

\section{T\_vbnc Adaptation Analysis}

\begin{figure}[h]
\centering
\includegraphics[width=\textwidth]{figures/fig5_tvbnc_adaptation.png}
\caption{\textbf{The key negative result}: Final $T_\mathrm{vbnc}$ per region across all adaptation runs. All values cluster within 0.15\textdegree C of the initial 12.0\textdegree C. The mechanism is inert at these rates.}
\end{figure}

\textbf{Why adaptation fails:} The daily adaptation is $\Delta T = r \times \mathrm{prevalence} \times (T_\mathrm{vbnc} - T)$. With prevalence $\sim$1--5\%, temp gap $\sim$2\textdegree C, and $r=0.001$: $\Delta T \approx 0.00002$\textdegree C/day $\approx$ 0.09\textdegree C over 13 years. Reversion (0.0005/day) further counteracts.

\subsection{T\_vbnc\_min Irrelevance}
\begin{figure}[h]
\centering
\includegraphics[width=\textwidth]{figures/fig7_tmin_irrelevance.png}
\caption{W96/W101/W102 produce \textit{identical} results despite $T_\mathrm{vbnc,min}$ = 6/8/4\textdegree C. The floor is never approached.}
\end{figure}

\section{RMSE Comparison}
\begin{figure}[h]
\centering
\includegraphics[width=\textwidth]{figures/fig6_rmse.png}
\caption{RMSE across all runs. W99 (0.557) is best, but improvement over W91 (0.591) is modest and not from adaptation.}
\end{figure}

\section{$\alpha_\mathrm{env}$ Effect}
\begin{figure}[h]
\centering
\includegraphics[width=\textwidth]{figures/fig8_alpha_effect.png}
\caption{$\alpha=0.15$ (W103/W104) vs $\alpha=0.20$ (W96). Lower $\alpha$ reduces gradient, inflates southern recovery.}
\end{figure}

\section{Key Findings}
\begin{enumerate}
\item \textbf{Pathogen adaptation is negligible} at biologically plausible rates over 13-year timescale
\item \textbf{$T_\mathrm{vbnc,min}$ is irrelevant} --- adaptation never approaches the biophysical floor
\item \textbf{W99 achieves RMSE=0.557} (best ever) due to floor=500 + $\alpha$=0.20, not adaptation
\item \textbf{Alaska recovery remains the core problem}: 10--12\% across all runs vs 50\% target
\item \textbf{$\alpha_\mathrm{env}$ dominates}: $\alpha$=0.20 produces strong gradient; $\alpha$=0.15 destroys it
\end{enumerate}

\section{Implications}
The mechanism as designed (prevalence-scaled daily increments) is too slow. Two paths forward:
\begin{itemize}
\item \textbf{Redesign}: Model strain \textit{replacement} rather than gradual adaptation --- binary or threshold-based rather than continuous
\item \textbf{Refocus}: Accept that the north-south gradient comes from dynamic $P_\mathrm{env}$ (floor + $\alpha$) and focus on lifting Alaska recovery through other mechanisms (recruitment, $K_{1/2}$, settler survival)
\end{itemize}

\end{document}
"""

with open(f"{OUT}/main.tex", "w") as f:
    f.write(tex)
print("LaTeX written")

# Compile
os.chdir(OUT)
for _ in range(2):
    subprocess.run([os.path.expanduser("~/.TinyTeX/bin/x86_64-linux/pdflatex"),
                    "-interaction=nonstopmode", "main.tex"],
                   capture_output=True, timeout=60)
print("PDF compiled" if os.path.exists(f"{OUT}/main.pdf") else "PDF FAILED")
