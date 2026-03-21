#!/usr/bin/env python3
"""Generate 12 disease ecology figures for the mechanistic model report.

Figures D01–D12 covering:
  - Compartmental model schematic
  - Force of infection decomposition
  - Environmental Vibrio dynamics
  - VBNC dormancy sigmoid
  - Arrhenius temperature scaling
  - Vibrio decay rate
  - Pathogen thermal adaptation
  - Community virulence evolution
  - Wavefront disease spread
  - Salinity suppression
  - R0 sensitivity analysis
  - Single-node epidemic trajectory

Output: PDF figures at 300 DPI in reports/mechanistic/figures/
"""

from __future__ import annotations

import sys
import json
import math
import traceback
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle
from matplotlib import patheffects
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from sswd_evoepi.config import DiseaseSection, PathogenEvolutionSection
from sswd_evoepi.disease import (
    arrhenius, thermal_performance, salinity_modifier, size_susceptibility,
    force_of_infection, environmental_vibrio, vibrio_decay_rate,
    shedding_rate_I1, shedding_rate_I2, compute_R0,
    run_single_node_epidemic, adapt_pathogen_thermal, NodeDiseaseState,
    T_REF_K, T_OPT, T_MAX, BETA_L, L_BAR, SIGMA_L,
    XI_10C, XI_20C, CARCASS_SHED_DAYS,
)

# ═══════════════════════════════════════════════════════════════════════
# STYLE SETUP
# ═══════════════════════════════════════════════════════════════════════

FIGDIR = PROJECT_ROOT / 'reports' / 'mechanistic' / 'figures'
FIGDIR.mkdir(parents=True, exist_ok=True)

# Disease state colors
C_S = '#3498db'    # blue — Susceptible
C_E = '#f39c12'    # gold — Exposed
C_I1 = '#e67e22'   # orange — I1
C_I2 = '#e74c3c'   # crimson — I2
C_D = '#2c3e50'    # dark — Dead
C_R = '#2ecc71'    # green — Recovered
C_P = '#9b59b6'    # purple — Pathogen/Tolerance

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['DejaVu Sans', 'Helvetica', 'Arial'],
    'font.size': 10,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.15,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'legend.framealpha': 0.9,
    'legend.edgecolor': '0.8',
})

cfg = DiseaseSection()
pe_cfg = PathogenEvolutionSection()

# Apply W285 calibrated parameters
cfg.K_half = 1_500_000
cfg.P_env_max = 2000
cfg.T_vbnc = 12.0  # initial
cfg.T_vbnc_min = 10.0
cfg.k_vbnc = 2.0
cfg.pathogen_adapt_rate = 0.001
cfg.v_max_warm = 0.7
cfg.v_adapt_rate = 0.001
cfg.alpha_env = 0.18
cfg.delta_env = 0.02
cfg.activation_threshold = 50.0
cfg.cumulative_dose_threshold = 1000.0
cfg.wavefront_D_P = 300.0
cfg.seed_vibrio = 2000.0

results = {}


def save_fig(fig, name: str):
    """Save figure and close."""
    out = FIGDIR / name
    fig.savefig(str(out), format='pdf')
    plt.close(fig)
    print(f"  ✓ Saved {name}")


# ═══════════════════════════════════════════════════════════════════════
# FIGURE D01: SEIPD+R Compartmental Model — Flow Diagram
# ═══════════════════════════════════════════════════════════════════════

def fig_D01():
    """SEIPD+R compartmental model flow diagram."""
    fig, ax = plt.subplots(1, 1, figsize=(14, 8))
    ax.set_xlim(-0.5, 14.5)
    ax.set_ylim(-1.5, 8.5)
    ax.axis('off')
    ax.set_aspect('equal')

    # Compartment positions and colors
    compartments = {
        'S': (1.5, 5.5, C_S, 'Susceptible\n(S)'),
        'E': (4.5, 5.5, C_E, 'Exposed\n(E)'),
        'I₁': (7.5, 5.5, C_I1, 'Pre-symptomatic\n(I₁)'),
        'I₂': (10.5, 5.5, C_I2, 'Symptomatic\n(I₂)'),
        'D': (13.0, 5.5, C_D, 'Dead\n(D)'),
    }

    box_w, box_h = 2.0, 1.6

    # Draw boxes
    for name, (cx, cy, color, label) in compartments.items():
        rect = FancyBboxPatch(
            (cx - box_w/2, cy - box_h/2), box_w, box_h,
            boxstyle="round,pad=0.15", facecolor=color, alpha=0.25,
            edgecolor=color, linewidth=2.5
        )
        ax.add_patch(rect)
        ax.text(cx, cy, label, ha='center', va='center',
                fontsize=11, fontweight='bold', color=color)

    # Forward arrows (S→E→I1→I2→D)
    arrow_style = "Simple,tail_width=4,head_width=14,head_length=8"
    arrow_kw = dict(arrowstyle=arrow_style, connectionstyle="arc3,rad=0",
                    lw=1.5)

    transitions = [
        ((2.5, 5.5), (3.5, 5.5), '#555', 'λᵢ\n(FOI)', 0.3),
        ((5.5, 5.5), (6.5, 5.5), '#555', 'μ_EI₁(T)\nk=3', 0.3),
        ((8.5, 5.5), (9.5, 5.5), '#555', 'μ_I₁I₂(T)\nk=2', 0.3),
        ((11.5, 5.5), (12.0, 5.5), '#555', 'μ_I₂D(T)\nk=2', 0.3),
    ]

    for (x1, y1), (x2, y2), color, label, offset in transitions:
        ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle='->', color=color, lw=2.5,
                                    connectionstyle='arc3,rad=0'))
        ax.text((x1+x2)/2, y1 + offset + 0.3, label, ha='center', va='bottom',
                fontsize=8, color='#333', style='italic')

    # Recovery arrows (I1→S and I2→S) - curved green dashed
    for src_x, src_label in [(7.5, 'I₁'), (10.5, 'I₂')]:
        ax.annotate("", xy=(1.5, 5.5 - box_h/2 - 0.1),
                    xytext=(src_x, 5.5 - box_h/2 - 0.1),
                    arrowprops=dict(arrowstyle='->', color=C_R, lw=2,
                                    linestyle='dashed',
                                    connectionstyle=f'arc3,rad=0.4'))

    ax.text(5.5, 3.3, 'Recovery (cᵢ × ρ_rec)', ha='center', va='center',
            fontsize=9, color=C_R, fontweight='bold')

    # Genetic trait shields
    # Resistance shield on S→E arrow
    ax.text(3.0, 6.7, '[R]', ha='center', va='center', fontsize=11,
            color=C_S, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.3', facecolor=C_S, alpha=0.15))
    ax.text(3.0, 7.4, 'Resistance\n(1−rᵢ)', ha='center', va='center',
            fontsize=7, color=C_S)

    # Tolerance shield on I2→D arrow
    ax.text(11.8, 6.7, '[T]', ha='center', va='center', fontsize=11,
            color=C_P, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.3', facecolor=C_P, alpha=0.15))
    ax.text(11.8, 7.4, 'Tolerance\n(τ_max × tᵢ)', ha='center', va='center',
            fontsize=7, color=C_P)

    # Recovery shield on I1/I2→S
    ax.text(9.0, 3.8, '[C]', ha='center', va='center', fontsize=11,
            color=C_R, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.3', facecolor=C_R, alpha=0.15))

    # Behavioral modifier table
    table_data = [
        ['', 'S', 'E', 'I₁', 'I₂', 'D'],
        ['Speed', '1.0', '1.0', '0.5', '0.1', '0.0'],
        ['Feeding', '1.0', '1.0', '0.5', '0.0', '0.0'],
        ['Can spawn', '✓', '✓', '✗', '✗', '✗'],
    ]

    colors_row = ['white', C_S, C_E, C_I1, C_I2, C_D]

    table = ax.table(
        cellText=table_data, loc='bottom',
        bbox=[0.1, -0.12, 0.8, 0.22],
        cellLoc='center',
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor('#ccc')
        if row == 0:
            cell.set_text_props(fontweight='bold')
            if col > 0:
                cell.set_facecolor(colors_row[col])
                cell.set_alpha(0.2)

    # Carcass shedding annotation
    ax.text(13.0, 4.2, 'σ_D shedding\n(≤3 days)', ha='center', va='center',
            fontsize=8, color=C_D, style='italic',
            bbox=dict(boxstyle='round,pad=0.2', facecolor='#eee', alpha=0.7))

    # Title and parameter box
    ax.set_title('SEIPD+R Compartmental Disease Model with Three-Trait Host Genetics',
                 fontsize=14, fontweight='bold', pad=20)

    param_text = (
        'Parameters at T_ref = 20°C:\n'
        'μ_EI₁ = 0.233 d⁻¹ (4.3d)\n'
        'μ_I₁I₂ = 0.434 d⁻¹ (2.3d)\n'
        'μ_I₂D = 0.563 d⁻¹ (1.8d)\n'
        'ρ_rec = 0.05 d⁻¹'
    )
    ax.text(0.0, 1.8, param_text, fontsize=7.5, family='monospace',
            va='top', ha='left',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='#f8f8f8',
                      edgecolor='#ccc', alpha=0.9))

    save_fig(fig, 'fig_D01.pdf')


# ═══════════════════════════════════════════════════════════════════════
# FIGURE D02: Force of Infection — Component Decomposition
# ═══════════════════════════════════════════════════════════════════════

def fig_D02():
    """Force of infection: 4-panel decomposition."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('Force of Infection Decomposition\n'
                 'λᵢ = a × P/(K_half+P) × (1−rᵢ) × S_sal × f_size(Lᵢ)',
                 fontsize=14, fontweight='bold')

    # Panel A: Dose-response
    ax = axes[0, 0]
    P = np.logspace(0, 7.5, 500)
    K_half = cfg.K_half
    dr = P / (K_half + P)
    ax.semilogx(P, dr, color=C_I2, lw=2.5)
    ax.axvline(K_half, ls='--', color='gray', alpha=0.7, lw=1.5)
    ax.text(K_half * 1.3, 0.15, f'K_half = {K_half/1e6:.1f}M', fontsize=9,
            color='gray')
    ax.set_xlabel('Vibrio concentration P (bact/mL)')
    ax.set_ylabel('Dose-response P/(K_half+P)')
    ax.set_title('(A) Dose-Response (Michaelis-Menten)', fontweight='bold')
    ax.set_ylim(-0.05, 1.05)
    ax.fill_between(P, dr, alpha=0.1, color=C_I2)

    # Panel B: Resistance modifier
    ax = axes[0, 1]
    r = np.linspace(0, 1, 200)
    ax.plot(r, 1 - r, color=C_S, lw=2.5)
    # Population distribution overlay
    from scipy.stats import beta as beta_dist
    x_beta = np.linspace(0, 0.6, 200)
    pdf = beta_dist.pdf(x_beta, 2, 8)
    pdf_scaled = pdf / pdf.max() * 0.3
    ax.fill_between(x_beta, 0, pdf_scaled, alpha=0.2, color=C_S,
                     label='Pop. distribution\nBeta(2,8)')
    ax.plot(r, 1 - r, color=C_S, lw=2.5, label='(1 − rᵢ)')
    ax.set_xlabel('Resistance score rᵢ')
    ax.set_ylabel('Susceptibility (1 − rᵢ)')
    ax.set_title('(B) Resistance Modifier', fontweight='bold')
    ax.legend(loc='upper right', fontsize=8)
    ax.set_ylim(-0.05, 1.15)

    # Panel C: Salinity modifier
    ax = axes[1, 0]
    S = np.linspace(0, 36, 300)
    S_sal = np.array([salinity_modifier(s, cfg.s_min, cfg.s_full) for s in S])
    ax.plot(S, S_sal, color='#1abc9c', lw=2.5)
    ax.axvline(cfg.s_min, ls=':', color='gray', alpha=0.7)
    ax.axvline(cfg.s_full, ls=':', color='gray', alpha=0.7)
    ax.text(cfg.s_min + 0.3, 0.85, f's_min={cfg.s_min}', fontsize=8, color='gray')
    ax.text(cfg.s_full + 0.3, 0.85, f's_full={cfg.s_full}', fontsize=8, color='gray')
    ax.fill_between(S, S_sal, alpha=0.1, color='#1abc9c')
    # Shade freshwater zone
    ax.axvspan(0, cfg.s_min, alpha=0.1, color='blue', label='Freshwater (Vibrio unviable)')
    ax.set_xlabel('Salinity (psu)')
    ax.set_ylabel('Salinity modifier S_sal')
    ax.set_title('(C) Salinity Modifier (Quadratic)', fontweight='bold')
    ax.legend(loc='center right', fontsize=8)
    ax.set_ylim(-0.05, 1.15)

    # Panel D: Size susceptibility
    ax = axes[1, 1]
    L = np.linspace(100, 900, 300)
    f_size = np.exp(BETA_L * (L - L_BAR) / SIGMA_L)
    ax.plot(L, f_size, color=C_E, lw=2.5)
    ax.axvline(L_BAR, ls='--', color='gray', alpha=0.7, lw=1.5)
    ax.text(L_BAR + 10, ax.get_ylim()[0] + 0.1, f'L̄={L_BAR}mm', fontsize=9,
            color='gray')
    # Mark Eisenlord OR
    ax.annotate('OR=1.23\nper 10mm', xy=(500, np.exp(BETA_L * (500 - L_BAR) / SIGMA_L)),
                xytext=(600, 1.0), fontsize=8, color=C_E,
                arrowprops=dict(arrowstyle='->', color=C_E, lw=1))
    ax.fill_between(L, f_size, 1, where=f_size > 1, alpha=0.1, color=C_I2)
    ax.fill_between(L, f_size, 1, where=f_size < 1, alpha=0.1, color=C_S)
    ax.set_xlabel('Body size L (mm)')
    ax.set_ylabel('Relative susceptibility f_size')
    ax.set_title('(D) Size-Dependent Susceptibility', fontweight='bold')

    fig.tight_layout(rect=[0, 0, 1, 0.93])
    save_fig(fig, 'fig_D02.pdf')


# ═══════════════════════════════════════════════════════════════════════
# FIGURE D03: Environmental Vibrio Dynamics
# ═══════════════════════════════════════════════════════════════════════

def fig_D03():
    """Environmental Vibrio dynamics: schematic + simulation."""
    fig = plt.figure(figsize=(15, 5.5))
    gs = gridspec.GridSpec(1, 3, width_ratios=[1, 1.2, 1.2], wspace=0.35)

    # Panel A: Schematic of dP/dt
    ax_a = fig.add_subplot(gs[0])
    ax_a.set_xlim(-1, 11)
    ax_a.set_ylim(-1, 9)
    ax_a.axis('off')
    ax_a.set_title('(A) Vibrio Budget: dP/dt', fontweight='bold', fontsize=11)

    # Central P_k box
    rect = FancyBboxPatch((3.5, 3.5), 4, 2, boxstyle="round,pad=0.2",
                          facecolor=C_P, alpha=0.2, edgecolor=C_P, lw=2)
    ax_a.add_patch(rect)
    ax_a.text(5.5, 4.5, 'P_k\n(bact/mL)', ha='center', va='center',
              fontsize=11, fontweight='bold', color=C_P)

    # Incoming arrows
    inputs = [
        (5.5, 8.5, 5.5, 5.8, C_I2, 'Shedding\nσ₁I₁ + σ₂I₂ + σ_D·D'),
        (0.0, 4.5, 3.2, 4.5, C_R, 'P_env\n(reservoir)'),
        (0.0, 3.0, 3.2, 3.8, '#888', 'Dispersal\nΣ d_jk·P_j'),
    ]
    for x1, y1, x2, y2, color, label in inputs:
        ax_a.annotate("", xy=(x2, y2), xytext=(x1, y1),
                       arrowprops=dict(arrowstyle='->', color=color, lw=2))
        ax_a.text(x1, y1 + 0.2, label, ha='center', va='bottom',
                  fontsize=7.5, color=color)

    # Outgoing arrows
    outputs = [
        (7.8, 4.5, 10.5, 4.5, '#c0392b', 'Decay\nξ(T)·P'),
        (5.5, 3.2, 5.5, 0.5, '#2980b9', 'Flushing\nφ·P'),
    ]
    for x1, y1, x2, y2, color, label in outputs:
        ax_a.annotate("", xy=(x2, y2), xytext=(x1, y1),
                       arrowprops=dict(arrowstyle='->', color=color, lw=2))
        ax_a.text(x2, y2 - 0.3, label, ha='center', va='top',
                  fontsize=7.5, color=color)

    # Panel B: Simulated epidemic Vibrio time series
    ax_b = fig.add_subplot(gs[1])
    result = run_single_node_epidemic(
        n_individuals=500, T_celsius=15.0, salinity=32.0,
        phi_k=0.02, cfg=cfg, n_days=365, initial_infected=5,
        record_daily=True, seed=42, initial_vibrio=100.0,
    )
    days = np.arange(result.days)
    ax_b.semilogy(days, np.maximum(result.daily_P, 0.1), color=C_P, lw=2,
                  label='P_k (Vibrio)')
    ax_b2 = ax_b.twinx()
    ax_b2.plot(days, result.daily_I, color=C_I2, lw=1.5, alpha=0.7,
               label='Infected (I₁+I₂)')
    ax_b2.set_ylabel('Infected count', color=C_I2)
    ax_b.set_xlabel('Days')
    ax_b.set_ylabel('Vibrio (bact/mL)', color=C_P)
    ax_b.set_title('(B) Single-Node Epidemic: Vibrio Dynamics', fontweight='bold',
                    fontsize=11)
    # Combined legend
    lines1, labels1 = ax_b.get_legend_handles_labels()
    lines2, labels2 = ax_b2.get_legend_handles_labels()
    ax_b.legend(lines1 + lines2, labels1 + labels2, loc='upper right', fontsize=8)

    # Panel C: P_env_pool dynamics
    ax_c = fig.add_subplot(gs[2])
    # Simulate P_env_pool buildup/decay analytically
    n_days_sim = 365
    P_pool = np.zeros(n_days_sim)
    shedding_profile = np.zeros(n_days_sim)
    # Approximate shedding from epidemic curve
    for d in range(n_days_sim):
        n_inf = result.daily_I[d] if d < len(result.daily_I) else 0
        shed = shedding_rate_I2(15.0, cfg) * n_inf
        shedding_profile[d] = shed
        pool_input = cfg.alpha_env * shed
        pool_decay = cfg.delta_env * P_pool[max(0, d-1)]
        P_pool[d] = max(0, (P_pool[max(0, d-1)] if d > 0 else 0) + pool_input - pool_decay)

    ax_c.plot(days, P_pool, color='#8e44ad', lw=2.5, label='P_env_pool')
    ax_c.fill_between(days, 0, P_pool, alpha=0.15, color='#8e44ad')
    ax_c.set_xlabel('Days')
    ax_c.set_ylabel('P_env_pool (bact/mL)')
    ax_c.set_title('(C) Host-Amplified Environmental Pool', fontweight='bold',
                    fontsize=11)
    ax_c.text(0.95, 0.85, f'α_env={cfg.alpha_env}\nδ_env={cfg.delta_env}',
              transform=ax_c.transAxes, fontsize=8, ha='right', va='top',
              bbox=dict(boxstyle='round,pad=0.3', facecolor='#f8f8f8',
                        edgecolor='#ccc'))
    ax_c.legend(fontsize=8)

    fig.suptitle('Environmental Vibrio Concentration Dynamics',
                 fontsize=14, fontweight='bold', y=1.02)
    save_fig(fig, 'fig_D03.pdf')


# ═══════════════════════════════════════════════════════════════════════
# FIGURE D04: VBNC Dormancy — Sigmoid Activation Function
# ═══════════════════════════════════════════════════════════════════════

def fig_D04():
    """VBNC dormancy: sigmoid activation, thermal performance, combined P_env."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    T = np.linspace(0, 30, 300)
    k = cfg.k_vbnc

    # Panel A: VBNC sigmoid for different T_vbnc values
    ax = axes[0]
    T_vbnc_vals = [9, 10, 11, 12]
    colors_vbnc = ['#e74c3c', '#e67e22', '#f39c12', '#3498db']
    for T_vbnc, c in zip(T_vbnc_vals, colors_vbnc):
        sigmoid = 1.0 / (1.0 + np.exp(-k * (T - T_vbnc)))
        ax.plot(T, sigmoid, color=c, lw=2, label=f'T_vbnc={T_vbnc}°C')
    # Shade Alaska winter
    ax.axvspan(3, 5, alpha=0.15, color='#3498db', label='Alaska winter SST')
    ax.axvspan(18, 22, alpha=0.1, color='#e74c3c', label='SoCal summer SST')
    ax.set_xlabel('Temperature (°C)')
    ax.set_ylabel('VBNC Activation [0,1]')
    ax.set_title('(A) VBNC Resuscitation Sigmoid', fontweight='bold')
    ax.legend(fontsize=7.5, loc='center right')
    ax.set_ylim(-0.05, 1.1)

    # Panel B: Thermal performance curve
    ax = axes[1]
    g_peak = np.array([thermal_performance(3000.0, t) for t in T])
    ax.plot(T, g_peak, color=C_I1, lw=2.5)
    ax.fill_between(T, g_peak, alpha=0.1, color=C_I1)
    ax.axvline(T_OPT, ls='--', color='gray', alpha=0.5)
    ax.text(T_OPT + 0.3, 0.9, f'T_opt={T_OPT}°C', fontsize=8, color='gray')
    ax.axvline(T_MAX, ls=':', color='gray', alpha=0.5)
    ax.text(T_MAX - 3, 0.1, f'T_max={T_MAX}°C', fontsize=8, color='gray')
    ax.set_xlabel('Temperature (°C)')
    ax.set_ylabel('Thermal Performance g(T)')
    ax.set_title('(B) Thermal Performance Curve', fontweight='bold')

    # Panel C: Combined P_env
    ax = axes[2]
    for T_vbnc, c in zip(T_vbnc_vals, colors_vbnc):
        P_env = []
        for t in T:
            sig = 1.0 / (1.0 + np.exp(-k * (t - T_vbnc)))
            gp = thermal_performance(3000.0, t)
            sal = salinity_modifier(32.0, cfg.s_min, cfg.s_full)
            P_env.append(cfg.P_env_max * sig * gp * sal)
        ax.plot(T, P_env, color=c, lw=2, label=f'T_vbnc={T_vbnc}°C')
    ax.axvspan(3, 5, alpha=0.1, color='#3498db')
    ax.axvspan(18, 22, alpha=0.08, color='#e74c3c')
    ax.set_xlabel('Temperature (°C)')
    ax.set_ylabel('P_env (bact/mL/d)')
    ax.set_title('(C) Combined Environmental Vibrio Input', fontweight='bold')
    ax.legend(fontsize=7.5)
    ax.text(0.02, 0.95, f'P_env_max = {cfg.P_env_max}', transform=ax.transAxes,
            fontsize=8, va='top',
            bbox=dict(boxstyle='round', facecolor='#f8f8f8', edgecolor='#ccc'))

    fig.suptitle('VBNC Dormancy: Temperature-Dependent Vibrio Resuscitation',
                 fontsize=14, fontweight='bold')
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    save_fig(fig, 'fig_D04.pdf')


# ═══════════════════════════════════════════════════════════════════════
# FIGURE D05: Arrhenius Temperature Scaling of Disease Rates
# ═══════════════════════════════════════════════════════════════════════

def fig_D05():
    """Arrhenius temperature scaling of disease progression."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    T = np.linspace(5, 25, 200)

    rates_info = {
        'E→I₁': (cfg.mu_EI1_ref, cfg.Ea_EI1, C_E),
        'I₁→I₂': (cfg.mu_I1I2_ref, cfg.Ea_I1I2, C_I1),
        'I₂→D': (cfg.mu_I2D_ref, cfg.Ea_I2D, C_I2),
        'Shedding': (5.0, cfg.Ea_sigma, C_P),
    }

    # Panel A: Normalized rates
    ax = axes[0]
    for name, (ref, Ea, color) in rates_info.items():
        vals = np.array([arrhenius(ref, Ea, t) for t in T])
        ax.semilogy(T, vals / ref, color=color, lw=2, label=f'{name} (E_a/R={int(Ea)}K)')
    ax.axhline(1.0, ls=':', color='gray', alpha=0.5)
    ax.axvline(20, ls='--', color='gray', alpha=0.3, label='T_ref=20°C')
    ax.set_xlabel('Temperature (°C)')
    ax.set_ylabel('Rate / Rate_ref (dimensionless)')
    ax.set_title('(A) Normalized Arrhenius Rates', fontweight='bold')
    ax.legend(fontsize=7.5)

    # Panel B: Mean stage durations (absolute)
    ax = axes[1]
    stage_info = {
        'E (incubation)': (cfg.mu_EI1_ref, cfg.Ea_EI1, C_E),
        'I₁ (pre-symptomatic)': (cfg.mu_I1I2_ref, cfg.Ea_I1I2, C_I1),
        'I₂ (symptomatic)': (cfg.mu_I2D_ref, cfg.Ea_I2D, C_I2),
    }
    dur_total = np.zeros_like(T)
    for name, (ref, Ea, color) in stage_info.items():
        rates = np.array([arrhenius(ref, Ea, t) for t in T])
        duration = 1.0 / rates
        ax.plot(T, duration, color=color, lw=2, label=name)
        dur_total += duration

    ax.set_xlabel('Temperature (°C)')
    ax.set_ylabel('Mean Stage Duration (days)')
    ax.set_title('(B) Stage Durations vs. Temperature', fontweight='bold')
    ax.legend(fontsize=8)

    # Panel C: Total exposure-to-death time
    ax = axes[2]
    ax.plot(T, dur_total, color=C_D, lw=2.5)
    ax.fill_between(T, dur_total, alpha=0.1, color=C_D)
    # Mark empirical reference: 11.6d at 13°C
    ax.plot(13, 11.6, 'o', color='red', ms=8, zorder=5)
    ax.annotate('Prentice et al. 2025\n11.6d at 13°C', xy=(13, 11.6),
                xytext=(16, 18), fontsize=8, color='red',
                arrowprops=dict(arrowstyle='->', color='red', lw=1.2))
    ax.set_xlabel('Temperature (°C)')
    ax.set_ylabel('Total E → D Duration (days)')
    ax.set_title('(C) Total Exposure-to-Death Time', fontweight='bold')

    fig.suptitle('Arrhenius Temperature Dependence of Disease Progression Rates',
                 fontsize=14, fontweight='bold')
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    save_fig(fig, 'fig_D05.pdf')


# ═══════════════════════════════════════════════════════════════════════
# FIGURE D06: Vibrio Decay Rate — Temperature Dependence
# ═══════════════════════════════════════════════════════════════════════

def fig_D06():
    """Vibrio decay rate and half-life vs temperature."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))

    T = np.linspace(5, 25, 200)
    xi = np.array([vibrio_decay_rate(t) for t in T])
    half_life = np.log(2) / xi

    # Region shadings
    regions = [
        ('Alaska\nsummer', 5, 10, '#3498db', 0.12),
        ('PNW\nsummer', 10, 14, '#1abc9c', 0.12),
        ('SoCal\nsummer', 16, 22, '#e74c3c', 0.10),
    ]

    # Panel A: Decay rate
    ax = axes[0]
    ax.plot(T, xi, color='#c0392b', lw=2.5)
    ax.fill_between(T, xi, alpha=0.1, color='#c0392b')
    for name, t_lo, t_hi, color, alpha in regions:
        ax.axvspan(t_lo, t_hi, alpha=alpha, color=color)
        ax.text((t_lo + t_hi) / 2, ax.get_ylim()[0] + 0.03, name,
                ha='center', va='bottom', fontsize=7, color=color)
    # Mark known points
    ax.plot(10, XI_10C, 'o', color='#c0392b', ms=8, zorder=5)
    ax.plot(20, XI_20C, 'o', color='#c0392b', ms=8, zorder=5)
    ax.text(10.5, XI_10C + 0.03, f'ξ(10°C)={XI_10C}', fontsize=8)
    ax.text(17, XI_20C + 0.03, f'ξ(20°C)={XI_20C}', fontsize=8)
    ax.set_xlabel('Temperature (°C)')
    ax.set_ylabel('Decay rate ξ (d⁻¹)')
    ax.set_title('(A) Vibrio Decay Rate', fontweight='bold')

    # Panel B: Half-life
    ax = axes[1]
    ax.plot(T, half_life, color='#2980b9', lw=2.5)
    ax.fill_between(T, half_life, alpha=0.1, color='#2980b9')
    for name, t_lo, t_hi, color, alpha in regions:
        ax.axvspan(t_lo, t_hi, alpha=alpha, color=color)
    ax.set_xlabel('Temperature (°C)')
    ax.set_ylabel('Half-life ln(2)/ξ (days)')
    ax.set_title('(B) Vibrio Environmental Half-Life', fontweight='bold')

    # Annotation: warm → persists longer
    ax.annotate('Vibrio persists ~3× longer\nin warm water (positive\nfeedback with disease)',
                xy=(20, half_life[np.argmin(np.abs(T - 20))]),
                xytext=(14, 2.5), fontsize=8, color='#2c3e50',
                arrowprops=dict(arrowstyle='->', color='#2c3e50', lw=1.2),
                bbox=dict(boxstyle='round,pad=0.3', facecolor='#fff3cd',
                          edgecolor='#ffc107', alpha=0.8))

    fig.suptitle('Temperature-Dependent Vibrio Decay in Seawater',
                 fontsize=14, fontweight='bold')
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    save_fig(fig, 'fig_D06.pdf')


# ═══════════════════════════════════════════════════════════════════════
# FIGURE D07: Pathogen Thermal Adaptation
# ═══════════════════════════════════════════════════════════════════════

def fig_D07():
    """Pathogen thermal adaptation: T_vbnc evolution at different latitudes."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Panel A: Schematic — sigmoid shifting left
    ax = axes[0]
    T_plot = np.linspace(0, 25, 200)
    k = cfg.k_vbnc
    for T_vbnc, ls, alpha in [(12, '-', 0.4), (10.5, '--', 0.6), (9, '-', 1.0)]:
        sig = 1.0 / (1.0 + np.exp(-k * (T_plot - T_vbnc)))
        ax.plot(T_plot, sig, ls=ls, alpha=alpha, color=C_I2, lw=2,
                label=f'T_vbnc={T_vbnc}°C')
    ax.annotate('', xy=(8, 0.5), xytext=(12, 0.5),
                arrowprops=dict(arrowstyle='->', color='#2c3e50', lw=2.5))
    ax.text(10, 0.58, 'Adaptation\n(cold selection)', ha='center', fontsize=9,
            fontweight='bold', color='#2c3e50')
    ax.axhline(0.5, ls=':', color='gray', alpha=0.3)
    ax.set_xlabel('Temperature (°C)')
    ax.set_ylabel('VBNC Activation')
    ax.set_title('(A) Sigmoid Shift via Adaptation', fontweight='bold')
    ax.legend(fontsize=8)

    # Panel B: T_vbnc trajectories at different latitudes
    ax = axes[1]

    # Representative sites: (name, mean_SST, SST_amplitude, color)
    sites = [
        ('Alaska 57°N', 7.0, 4.0, '#3498db'),
        ('Haida Gwaii 53°N', 9.0, 3.5, '#1abc9c'),
        ('PNW 48°N', 10.5, 3.0, '#f39c12'),
        ('Central CA 37°N', 13.0, 3.0, '#e67e22'),
        ('SoCal 34°N', 16.0, 3.0, '#e74c3c'),
    ]

    n_years = 10
    n_days = n_years * 365
    for site_name, mean_sst, amp, color in sites:
        node = NodeDiseaseState(T_vbnc_local=12.0)
        T_vbnc_history = []
        for day in range(n_days):
            # Sinusoidal SST
            T_sst = mean_sst + amp * np.sin(2 * np.pi * (day - 45) / 365)
            adapt_pathogen_thermal(node, T_sst, P_env_pool=500.0, cfg=cfg)
            if day % 30 == 0:
                T_vbnc_history.append(node.T_vbnc_local)
        years = np.arange(len(T_vbnc_history)) * 30 / 365
        ax.plot(years, T_vbnc_history, color=color, lw=1.8, label=site_name)

    ax.axhline(cfg.T_vbnc_min, ls='--', color='gray', alpha=0.5,
               label=f'T_vbnc_min={cfg.T_vbnc_min}°C')
    ax.set_xlabel('Years')
    ax.set_ylabel('T_vbnc (°C)')
    ax.set_title('(B) T_vbnc Trajectories by Latitude', fontweight='bold')
    ax.legend(fontsize=7, loc='upper right')

    # Panel C: Latitudinal gradient of equilibrium T_vbnc
    ax = axes[2]
    lats = np.linspace(28, 61, 50)
    eq_T_vbnc = []
    mean_ssts = []
    for lat in lats:
        # Rough SST model: decreasing with latitude
        mean_sst = max(4, 25 - 0.35 * (lat - 28))
        mean_ssts.append(mean_sst)
        amp = 3.0 + 0.03 * (lat - 28)
        node = NodeDiseaseState(T_vbnc_local=12.0)
        for day in range(3650):
            T_sst = mean_sst + amp * np.sin(2 * np.pi * (day - 45) / 365)
            adapt_pathogen_thermal(node, T_sst, P_env_pool=500.0, cfg=cfg)
        eq_T_vbnc.append(node.T_vbnc_local)

    ax.plot(lats, eq_T_vbnc, 'o-', color=C_I2, ms=4, lw=1.5,
            label='Equilibrium T_vbnc')
    ax2 = ax.twinx()
    ax2.plot(lats, mean_ssts, '--', color='#2980b9', lw=1.5, alpha=0.7,
             label='Mean annual SST')
    ax2.set_ylabel('Mean SST (°C)', color='#2980b9')
    ax.set_xlabel('Latitude (°N)')
    ax.set_ylabel('Equilibrium T_vbnc (°C)', color=C_I2)
    ax.set_title('(C) Latitudinal T_vbnc Gradient', fontweight='bold')
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, fontsize=8, loc='upper right')

    fig.suptitle('Pathogen Thermal Adaptation: Cold-Tolerant Vibrio Strains',
                 fontsize=14, fontweight='bold')
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    save_fig(fig, 'fig_D07.pdf')


# ═══════════════════════════════════════════════════════════════════════
# FIGURE D08: Community Virulence Evolution
# ═══════════════════════════════════════════════════════════════════════

def fig_D08():
    """Community virulence evolution: tradeoff and optimum surface."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Panel A: Virulence optimum heatmap
    ax = axes[0]
    T_grid = np.linspace(5, 25, 100)
    d_grid = np.linspace(0.01, 1.0, 100)
    T_v_mid, T_v_width, v_max_warm = cfg.T_v_mid, cfg.T_v_width, cfg.v_max_warm

    v_opt = np.zeros((100, 100))
    for i, d in enumerate(d_grid):
        for j, T in enumerate(T_grid):
            tf = 1.0 / (1.0 + math.exp(-(T - T_v_mid) / max(T_v_width, 0.1)))
            v_opt[i, j] = v_max_warm * d * tf

    im = ax.pcolormesh(T_grid, d_grid, v_opt, cmap='YlOrRd', shading='auto')
    cs = ax.contour(T_grid, d_grid, v_opt, levels=[0.1, 0.3, 0.5], colors='k',
                     linewidths=1)
    ax.clabel(cs, fmt='v*=%.1f', fontsize=8)
    cbar = fig.colorbar(im, ax=ax, label='Optimal virulence v*', shrink=0.8)
    ax.set_xlabel('Temperature (°C)')
    ax.set_ylabel('Density ratio N/K')
    ax.set_title('(A) Virulence Optimum Surface', fontweight='bold')

    # Panel B: Tradeoff curves
    ax = axes[1]
    v = np.linspace(0, 1, 200)
    v_anchor = pe_cfg.v_anchor
    dv = v - v_anchor
    shed_mult = np.exp(pe_cfg.alpha_shed * dv)
    prog_mult = np.exp(pe_cfg.alpha_prog * dv)
    kill_mult = np.exp(pe_cfg.alpha_kill * dv)

    ax.plot(v, shed_mult, color=C_P, lw=2, label=f'Shedding (α={pe_cfg.alpha_shed})')
    ax.plot(v, prog_mult, color=C_I1, lw=2, label=f'Progression (α={pe_cfg.alpha_prog})')
    ax.plot(v, kill_mult, color=C_I2, lw=2, label=f'Kill rate (α={pe_cfg.alpha_kill})')
    ax.axvline(v_anchor, ls='--', color='gray', alpha=0.5)
    ax.text(v_anchor + 0.02, ax.get_ylim()[1] * 0.9, f'v_anchor={v_anchor}',
            fontsize=8, color='gray')
    ax.set_xlabel('Virulence v')
    ax.set_ylabel('Rate multiplier (relative to v_anchor)')
    ax.set_title('(B) Virulence-Transmission Tradeoff', fontweight='bold')
    ax.legend(fontsize=8)
    ax.set_yscale('log')

    # Panel C: Virulence trajectory at contrasting sites
    ax = axes[2]
    # Simulate two sites: warm/dense vs cold/sparse
    scenarios = [
        ('Warm + Dense (SoCal)', 16.0, 0.8, 5000, '#e74c3c'),
        ('Cool + Moderate (PNW)', 11.0, 0.5, 5000, '#f39c12'),
        ('Cold + Sparse (Alaska)', 7.0, 0.2, 5000, '#3498db'),
    ]
    for label, mean_T, density_ratio, K, color in scenarios:
        node = NodeDiseaseState(v_local=0.5, P_env_pool=500.0)
        v_history = []
        n_hosts = int(density_ratio * K)
        for day in range(3650):
            T_sst = mean_T + 3 * np.sin(2 * np.pi * (day - 45) / 365)
            # Simulate population crash and recovery
            if 365 < day < 730:
                n_hosts_now = int(n_hosts * 0.2)  # crash
            elif 730 <= day < 1460:
                frac = (day - 730) / 730
                n_hosts_now = int(n_hosts * (0.2 + 0.6 * frac))  # recovery
            else:
                n_hosts_now = n_hosts
            from sswd_evoepi.disease import adapt_community_virulence
            adapt_community_virulence(node, T_sst, n_hosts_now, K, cfg)
            if day % 30 == 0:
                v_history.append(node.v_local)
        years = np.arange(len(v_history)) * 30 / 365
        ax.plot(years, v_history, color=color, lw=1.8, label=label)

    # Mark crash period
    ax.axvspan(1, 2, alpha=0.1, color='gray')
    ax.text(1.5, 0.05, 'Population\ncrash', ha='center', fontsize=7, color='gray')
    ax.set_xlabel('Years')
    ax.set_ylabel('Community virulence v_local')
    ax.set_title('(C) Virulence Trajectories', fontweight='bold')
    ax.legend(fontsize=7.5, loc='upper right')

    fig.suptitle('Community Virulence Evolution: Density × Temperature Trade-off',
                 fontsize=14, fontweight='bold')
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    save_fig(fig, 'fig_D08.pdf')


# ═══════════════════════════════════════════════════════════════════════
# FIGURE D09: Wavefront Disease Spread
# ═══════════════════════════════════════════════════════════════════════

def fig_D09():
    """Wavefront disease spread from Channel Islands: map + distance scatter."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 6))

    # Panel A: Wavefront mechanism schematic
    ax = axes[0]
    ax.set_xlim(-0.5, 10.5)
    ax.set_ylim(-0.5, 8.5)
    ax.axis('off')
    ax.set_title('(A) Wavefront Mechanism', fontweight='bold')

    # Active node
    circ1 = Circle((2.5, 5), 1.2, facecolor=C_I2, alpha=0.3, edgecolor=C_I2, lw=2)
    ax.add_patch(circ1)
    ax.text(2.5, 5, 'Active\nNode', ha='center', va='center', fontsize=9,
            fontweight='bold', color=C_I2)

    # Inactive node
    circ2 = Circle((7.5, 5), 1.2, facecolor='#bdc3c7', alpha=0.4,
                    edgecolor='#7f8c8d', lw=2)
    ax.add_patch(circ2)
    ax.text(7.5, 5, 'Inactive\nNode', ha='center', va='center', fontsize=9,
            fontweight='bold', color='#7f8c8d')

    # Dispersal arrow
    ax.annotate('', xy=(6.0, 5), xytext=(4.0, 5),
                arrowprops=dict(arrowstyle='->', color=C_P, lw=2.5,
                                connectionstyle='arc3,rad=0'))
    ax.text(5, 5.7, 'Pathogen\ndispersal', ha='center', fontsize=8,
            color=C_P, fontweight='bold')
    ax.text(5, 4.0, f'D_P = {int(cfg.wavefront_D_P)} km', ha='center',
            fontsize=8, color=C_P)

    # Dose accumulation bar
    bar_y = 2.0
    ax.barh(bar_y, 0.7, height=0.6, left=5.5, color='#bdc3c7', edgecolor='#7f8c8d')
    ax.barh(bar_y, 0.7, height=0.6, left=5.5, color=C_I2, alpha=0.5,
            edgecolor=C_I2)
    ax.text(6.6, bar_y, '→', fontsize=14, ha='left', va='center', color=C_I2)
    ax.text(7.5, bar_y, f'CDT = {int(cfg.cumulative_dose_threshold)}', fontsize=8,
            ha='center', va='center', color='#555')
    ax.plot([8.8, 8.8], [bar_y - 0.5, bar_y + 0.5], '--', color='red', lw=2)
    ax.text(9.0, bar_y + 0.6, 'Threshold', fontsize=7, color='red')
    ax.text(7.5, 1.0, 'Cumulative dose\naccumulation', ha='center',
            fontsize=8, color='#555')

    # Load W285 data
    npz = np.load(str(PROJECT_ROOT / 'reports' / 'mechanistic' / 'w285_monthly.npz'),
                  allow_pickle=True)
    site_lats = npz['site_lats']
    site_lons = npz['site_lons']
    infected = npz['infected']
    sim_days_arr = npz['sim_days']

    # Compute first infection time for each site
    first_inf_step = np.full(896, -1, dtype=int)
    for s in range(896):
        inf_at_site = infected[:, s]
        has_inf = np.where(inf_at_site > 0)[0]
        if len(has_inf) > 0:
            first_inf_step[s] = has_inf[0]

    arrival_days = np.where(first_inf_step >= 0,
                            sim_days_arr[first_inf_step],
                            sim_days_arr[-1])

    # Channel Islands centroid
    origin_lat, origin_lon = 34.0, -119.5

    # Panel B: Map of disease arrival
    ax = axes[1]
    # Sort by arrival day so later arrivals are on top
    order = np.argsort(arrival_days)
    sc = ax.scatter(site_lons[order], site_lats[order],
                    c=arrival_days[order] / 30,  # convert to months
                    cmap='plasma', s=12, alpha=0.85, edgecolors='none')
    # Mark origin
    origin_mask = np.isin(np.arange(896), [319, 322, 632, 633, 634])
    ax.scatter(site_lons[origin_mask], site_lats[origin_mask],
               marker='*', s=120, c='lime', edgecolors='k', zorder=10,
               label='Origin (Channel Is.)')
    cbar = fig.colorbar(sc, ax=ax, label='Disease arrival (months)', shrink=0.8)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title('(B) Disease Arrival Map', fontweight='bold')
    ax.legend(loc='upper left', fontsize=8)
    ax.set_aspect(1.3)

    # Panel C: Distance vs arrival day
    ax = axes[2]

    def haversine_km(lat1, lon1, lat2, lon2):
        R = 6371
        dlat = np.radians(lat2 - lat1)
        dlon = np.radians(lon2 - lon1)
        a = (np.sin(dlat/2)**2 +
             np.cos(np.radians(lat1)) * np.cos(np.radians(lat2)) * np.sin(dlon/2)**2)
        return R * 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))

    distances = np.array([haversine_km(origin_lat, origin_lon, lat, lon)
                          for lat, lon in zip(site_lats, site_lons)])

    valid = first_inf_step >= 0
    sc2 = ax.scatter(distances[valid], arrival_days[valid] / 30,
                     c=site_lats[valid], cmap='coolwarm', s=8, alpha=0.7)
    fig.colorbar(sc2, ax=ax, label='Latitude (°N)', shrink=0.8)

    # Fit line for effective speed
    from numpy.polynomial import polynomial as P
    coeffs = P.polyfit(distances[valid], arrival_days[valid], 1)
    x_fit = np.linspace(0, distances[valid].max(), 100)
    y_fit = P.polyval(x_fit, coeffs)
    speed_km_per_month = 1 / (coeffs[1] / 30) if coeffs[1] > 0 else 0
    ax.plot(x_fit, y_fit / 30, '--', color='#555', lw=1.5,
            label=f'~{speed_km_per_month:.0f} km/month')

    ax.set_xlabel('Distance from Channel Islands (km)')
    ax.set_ylabel('Disease arrival (months)')
    ax.set_title('(C) Wavefront Speed', fontweight='bold')
    ax.legend(fontsize=8)

    fig.suptitle('Wavefront Disease Spread from Channel Islands Origin',
                 fontsize=14, fontweight='bold')
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    save_fig(fig, 'fig_D09.pdf')


# ═══════════════════════════════════════════════════════════════════════
# FIGURE D10: Salinity and Freshwater Suppression
# ═══════════════════════════════════════════════════════════════════════

def fig_D10():
    """Salinity suppression: modifier function + seasonal dynamics."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Panel A: Salinity modifier function
    ax = axes[0]
    S = np.linspace(0, 36, 300)
    S_sal = np.array([salinity_modifier(s, cfg.s_min, cfg.s_full) for s in S])
    ax.plot(S, S_sal, color='#1abc9c', lw=2.5)
    ax.fill_between(S, S_sal, alpha=0.1, color='#1abc9c')
    ax.axvline(cfg.s_min, ls=':', color='gray', alpha=0.7)
    ax.axvline(cfg.s_full, ls=':', color='gray', alpha=0.7)
    ax.text(cfg.s_min - 1.5, 0.5, f's_min\n{cfg.s_min}', fontsize=8,
            ha='center', color='gray')
    ax.text(cfg.s_full + 1, 0.5, f's_full\n{cfg.s_full}', fontsize=8,
            ha='center', color='gray')
    ax.axvspan(0, cfg.s_min, alpha=0.1, color='#3498db')
    ax.text(5, 0.85, 'Vibrio\nunviable', ha='center', fontsize=8, color='#3498db')
    ax.set_xlabel('Salinity (psu)')
    ax.set_ylabel('S_sal [0, 1]')
    ax.set_title('(A) Salinity Modifier Function', fontweight='bold')

    # Panel B: Two-layer salinity schematic
    ax = axes[1]
    ax.set_xlim(-1, 11)
    ax.set_ylim(-1, 9)
    ax.axis('off')
    ax.set_title('(B) Two-Layer Salinity Model', fontweight='bold')

    # Layer 1: WOA23
    rect1 = FancyBboxPatch((0.5, 5), 9, 2.5, boxstyle="round,pad=0.2",
                            facecolor='#2980b9', alpha=0.2, edgecolor='#2980b9', lw=2)
    ax.add_patch(rect1)
    ax.text(5, 6.25, 'Layer 1: WOA23 Monthly Climatology\n(0.25° resolution, baseline)',
            ha='center', va='center', fontsize=9, color='#2980b9', fontweight='bold')

    # Layer 2: Freshwater
    rect2 = FancyBboxPatch((0.5, 1.5), 9, 2.5, boxstyle="round,pad=0.2",
                            facecolor='#1abc9c', alpha=0.2, edgecolor='#1abc9c', lw=2)
    ax.add_patch(rect2)
    ax.text(5, 2.75, 'Layer 2: Fjord Freshwater Depression\n'
            'fw_strength × fjord_depth^exp × melt_pulse(day)\n'
            '× lat_factor(lat, 48°–60°N)',
            ha='center', va='center', fontsize=8, color='#1abc9c')

    # Arrow: combination
    ax.annotate('', xy=(5, 4.5), xytext=(5, 4.8),
                arrowprops=dict(arrowstyle='<->', color='#555', lw=2))
    ax.text(5, 4.65, 'S_local = Layer1 − Layer2', ha='center', fontsize=9,
            color='#555', fontweight='bold')

    # Melt pulse annotation
    ax.text(5, 0.7, 'Melt pulse: cos²(π(day−166)/365), peak ~June 15',
            ha='center', fontsize=7.5, style='italic', color='#7f8c8d')

    # Panel C: Seasonal salinity + S_sal for representative sites
    ax = axes[2]
    days = np.arange(365)

    sites = [
        ('Open coast (Monterey)', 33.5, 0.0, '#e74c3c'),
        ('Puget Sound', 29.5, 3.0, '#f39c12'),
        ('Glacial fjord (Howe Sound)', 24.0, 10.0, '#3498db'),
    ]

    ax2 = ax.twinx()
    for name, base_sal, fw_amp, color in sites:
        # Melt pulse: cosine squared centered at day 166
        melt = np.maximum(0, np.cos(np.pi * (days - 166) / 182)) ** 2
        sal_seasonal = base_sal - fw_amp * melt
        sal_seasonal = np.maximum(sal_seasonal, 5)
        S_sal_seasonal = np.array([salinity_modifier(s, cfg.s_min, cfg.s_full)
                                    for s in sal_seasonal])
        ax.plot(days, sal_seasonal, color=color, lw=1.8, label=f'{name}')
        ax2.plot(days, S_sal_seasonal, color=color, lw=1.2, ls='--', alpha=0.6)

    # Shade refugia periods
    ax2.axhspan(0, 0.5, alpha=0.08, color='#3498db')
    ax2.text(300, 0.25, 'Refugia\n(S_sal<0.5)', fontsize=7, color='#3498db',
             ha='center')

    ax.set_xlabel('Day of Year')
    ax.set_ylabel('Salinity (psu)')
    ax2.set_ylabel('S_sal (dashed)', color='gray')
    ax.set_title('(C) Seasonal Salinity at Representative Sites', fontweight='bold')
    ax.legend(fontsize=7.5, loc='lower left')
    ax.axhline(cfg.s_min, ls=':', color='gray', alpha=0.3)
    ax.axhline(cfg.s_full, ls=':', color='gray', alpha=0.3)

    fig.suptitle('Freshwater Refugia: Salinity-Dependent Vibrio Suppression',
                 fontsize=14, fontweight='bold')
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    save_fig(fig, 'fig_D10.pdf')


# ═══════════════════════════════════════════════════════════════════════
# FIGURE D11: R₀ Sensitivity Analysis
# ═══════════════════════════════════════════════════════════════════════

def fig_D11():
    """R0 sensitivity: temperature, density, flushing heatmap."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Panel A: R0 vs temperature at different flushing rates
    ax = axes[0]
    T = np.linspace(5, 25, 200)
    phi_vals = [0.01, 0.02, 0.05]
    phi_colors = ['#e74c3c', '#f39c12', '#3498db']
    for phi, color in zip(phi_vals, phi_colors):
        R0 = np.array([compute_R0(t, 200, phi, cfg) for t in T])
        ax.semilogy(T, R0, color=color, lw=2, label=f'φ={phi}')
    ax.axhline(1, ls='--', color='gray', alpha=0.7, lw=1.5)
    ax.text(6, 1.3, 'R₀ = 1 (epidemic threshold)', fontsize=8, color='gray')
    ax.set_xlabel('Temperature (°C)')
    ax.set_ylabel('R₀')
    ax.set_title('(A) R₀ vs. Temperature', fontweight='bold')
    ax.legend(fontsize=8)
    ax.set_ylim(0.01, 100)

    # Panel B: R0 vs host density
    ax = axes[1]
    N_vals = np.linspace(10, 500, 200)
    R0_N = np.array([compute_R0(15, int(n), 0.02, cfg) for n in N_vals])
    ax.plot(N_vals, R0_N, color=C_I2, lw=2.5)
    ax.axhline(1, ls='--', color='gray', alpha=0.7, lw=1.5)
    # Find threshold N
    threshold_idx = np.where(R0_N >= 1)[0]
    if len(threshold_idx) > 0:
        N_thresh = N_vals[threshold_idx[0]]
        ax.axvline(N_thresh, ls=':', color=C_I2, alpha=0.5)
        ax.text(N_thresh + 5, ax.get_ylim()[1] * 0.5,
                f'N_thresh≈{int(N_thresh)}', fontsize=8, color=C_I2)
    ax.fill_between(N_vals, R0_N, 1, where=R0_N > 1, alpha=0.1, color=C_I2)
    ax.set_xlabel('Host count N')
    ax.set_ylabel('R₀')
    ax.set_title('(B) R₀ vs. Host Density (T=15°C)', fontweight='bold')

    # Panel C: R0 heatmap
    ax = axes[2]
    T_grid = np.linspace(5, 25, 80)
    phi_grid = np.logspace(-2.3, -1, 60)
    R0_map = np.zeros((len(phi_grid), len(T_grid)))
    for i, phi in enumerate(phi_grid):
        for j, t in enumerate(T_grid):
            R0_map[i, j] = compute_R0(t, 200, phi, cfg)

    im = ax.pcolormesh(T_grid, phi_grid, np.log10(np.maximum(R0_map, 0.01)),
                       cmap='RdYlBu_r', shading='auto', vmin=-1.5, vmax=2)
    cs = ax.contour(T_grid, phi_grid, R0_map, levels=[1.0], colors='black',
                     linewidths=2.5)
    ax.clabel(cs, fmt='R₀=1', fontsize=9)
    ax.set_yscale('log')
    fig.colorbar(im, ax=ax, label='log₁₀(R₀)', shrink=0.8)
    ax.set_xlabel('Temperature (°C)')
    ax.set_ylabel('Flushing rate φ (d⁻¹)')
    ax.set_title('(C) R₀ Parameter Space', fontweight='bold')

    # Annotate epidemic window
    ax.text(18, 0.012, 'Epidemic\nwindow', fontsize=9, color='white',
            fontweight='bold', ha='center')

    fig.suptitle('Basic Reproduction Number R₀: Sensitivity Analysis',
                 fontsize=14, fontweight='bold')
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    save_fig(fig, 'fig_D11.pdf')


# ═══════════════════════════════════════════════════════════════════════
# FIGURE D12: Single-Node Epidemic Trajectory
# ═══════════════════════════════════════════════════════════════════════

def fig_D12():
    """Single-node epidemic: compartment dynamics + sensitivity."""
    fig, axes = plt.subplots(2, 2, figsize=(13, 10))

    # Reference epidemic
    ref = run_single_node_epidemic(
        n_individuals=500, T_celsius=15.0, salinity=32.0,
        phi_k=0.02, cfg=cfg, n_days=300, initial_infected=5,
        record_daily=True, seed=42, initial_vibrio=100.0,
    )
    days = np.arange(ref.days)

    # Panel A: Compartment time series
    ax = axes[0, 0]
    ax.plot(days, ref.daily_S, color=C_S, lw=2, label='S')
    ax.plot(days, ref.daily_I, color=C_I2, lw=2, label='I₁+I₂')
    ax.plot(days, ref.daily_D_cumul, color=C_D, lw=2, ls='--', label='Cumul. D')
    ax.fill_between(days, 0, ref.daily_S, alpha=0.08, color=C_S)
    ax.fill_between(days, 0, ref.daily_I, alpha=0.12, color=C_I2)
    ax.set_xlabel('Days')
    ax.set_ylabel('Count')
    ax.set_title('(A) Compartment Dynamics (N₀=500, T=15°C)', fontweight='bold')
    ax.legend(fontsize=8)
    ax.text(0.95, 0.6, f'Total deaths: {ref.total_deaths}\n'
            f'Mortality: {ref.mortality_fraction:.1%}\n'
            f'Peak prevalence: {ref.peak_prevalence:.1%}',
            transform=ax.transAxes, fontsize=8, ha='right', va='top',
            bbox=dict(boxstyle='round', facecolor='#f8f8f8', edgecolor='#ccc'))

    # Panel B: Vibrio concentration
    ax = axes[0, 1]
    ax.semilogy(days, np.maximum(ref.daily_P, 0.1), color=C_P, lw=2)
    ax.fill_between(days, 0.1, np.maximum(ref.daily_P, 0.1), alpha=0.1, color=C_P)
    peak_day = np.argmax(ref.daily_P)
    peak_val = ref.daily_P[peak_day]
    ax.plot(peak_day, peak_val, 'v', color=C_I2, ms=10, zorder=5)
    ax.annotate(f'Peak: {peak_val:.0f}\n(day {peak_day})',
                xy=(peak_day, peak_val),
                xytext=(peak_day + 30, peak_val * 2),
                fontsize=8, arrowprops=dict(arrowstyle='->', color=C_I2))
    ax.axhline(cfg.K_half, ls=':', color='gray', alpha=0.3)
    ax.text(200, cfg.K_half * 1.3, f'K_half={cfg.K_half/1e6:.1f}M',
            fontsize=7, color='gray')
    ax.set_xlabel('Days')
    ax.set_ylabel('Vibrio concentration (bact/mL)')
    ax.set_title('(B) Environmental Vibrio', fontweight='bold')

    # Panel C: Temperature sensitivity
    ax = axes[1, 0]
    temp_colors = {10: '#3498db', 15: '#f39c12', 20: '#e74c3c'}
    for T, color in temp_colors.items():
        r = run_single_node_epidemic(
            n_individuals=500, T_celsius=float(T), salinity=32.0,
            phi_k=0.02, cfg=cfg, n_days=300, initial_infected=5,
            record_daily=True, seed=42, initial_vibrio=100.0,
        )
        mort_frac = r.daily_D_cumul / 500
        ax.plot(np.arange(r.days), mort_frac, color=color, lw=2,
                label=f'T={T}°C (mort={r.mortality_fraction:.0%})')
    ax.set_xlabel('Days')
    ax.set_ylabel('Cumulative mortality fraction')
    ax.set_title('(C) Temperature Sensitivity', fontweight='bold')
    ax.legend(fontsize=8)
    ax.set_ylim(0, 1.05)

    # Panel D: Resistance sensitivity
    ax = axes[1, 1]
    res_colors = {0.05: '#e74c3c', 0.15: '#f39c12', 0.30: '#2ecc71'}
    for r_mean, color in res_colors.items():
        r = run_single_node_epidemic(
            n_individuals=500, T_celsius=15.0, salinity=32.0,
            phi_k=0.02, cfg=cfg, n_days=300, initial_infected=5,
            record_daily=True, seed=42, initial_vibrio=100.0,
            mean_resistance=r_mean,
        )
        mort_frac = r.daily_D_cumul / 500
        ax.plot(np.arange(r.days), mort_frac, color=color, lw=2,
                label=f'r̄={r_mean} (mort={r.mortality_fraction:.0%})')
    ax.set_xlabel('Days')
    ax.set_ylabel('Cumulative mortality fraction')
    ax.set_title('(D) Resistance Sensitivity', fontweight='bold')
    ax.legend(fontsize=8)
    ax.set_ylim(0, 1.05)

    fig.suptitle('Anatomy of a Single-Node SSWD Epidemic',
                 fontsize=14, fontweight='bold')
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    save_fig(fig, 'fig_D12.pdf')


# ═══════════════════════════════════════════════════════════════════════
# MAIN — Generate all figures
# ═══════════════════════════════════════════════════════════════════════

FIGURES = [
    ('D01', 'SEIPD+R Flow Diagram', fig_D01),
    ('D02', 'Force of Infection Decomposition', fig_D02),
    ('D03', 'Environmental Vibrio Dynamics', fig_D03),
    ('D04', 'VBNC Dormancy Sigmoid', fig_D04),
    ('D05', 'Arrhenius Temperature Scaling', fig_D05),
    ('D06', 'Vibrio Decay Rate', fig_D06),
    ('D07', 'Pathogen Thermal Adaptation', fig_D07),
    ('D08', 'Community Virulence Evolution', fig_D08),
    ('D09', 'Wavefront Disease Spread', fig_D09),
    ('D10', 'Salinity Freshwater Suppression', fig_D10),
    ('D11', 'R₀ Sensitivity Analysis', fig_D11),
    ('D12', 'Single-Node Epidemic Trajectory', fig_D12),
]

if __name__ == '__main__':
    print(f"Generating {len(FIGURES)} disease ecology figures...")
    print(f"Output directory: {FIGDIR}\n")

    successes = []
    failures = []

    for fig_id, title, func in FIGURES:
        print(f"[{fig_id}] {title}")
        try:
            func()
            successes.append(fig_id)
        except Exception as e:
            print(f"  ✗ FAILED: {e}")
            traceback.print_exc()
            failures.append((fig_id, str(e)))

    print(f"\n{'='*60}")
    print(f"RESULTS: {len(successes)}/{len(FIGURES)} figures generated")
    print(f"  ✓ Succeeded: {', '.join(successes)}")
    if failures:
        print(f"  ✗ Failed: {', '.join(f'{fid}: {msg[:50]}' for fid, msg in failures)}")
    print(f"\nFigures saved to: {FIGDIR}")
