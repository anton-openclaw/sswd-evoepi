#!/usr/bin/env python3
"""
Faithful single-site SSWD-EvoEpi disease simulation.

Mechanistically mirrors the actual model:
- Correlated random walk movement (Kay & Emlet 2002)
- InfectedDensityGrid spatial transmission
- Dynamic environmental pathogen pool (P_env)
- Arrhenius-scaled disease progression at 13°C
- Stochastic daily transitions

Seed = 42, N = 5000, habitat = 25,000 m², T = 13°C
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from scipy.ndimage import uniform_filter

# ── Seed ──
np.random.seed(42)

# ── Constants ──
N = 5000
HABITAT_AREA = 25_000.0  # m²
HABITAT_SIDE = np.sqrt(HABITAT_AREA)  # ~158.1 m
CELL_SIZE = 20.0  # m
NX = int(np.ceil(HABITAT_SIDE / CELL_SIZE))  # 8
NY = NX

# Disease states
S, E, I1, I2, D, R = 0, 1, 2, 3, 4, 5
STATE_NAMES = {S: "S", E: "E", I1: "I1", I2: "I2", D: "D", R: "R"}

# Movement
BASE_SPEED = 0.5  # m/min
SIGMA_TURN = 0.6  # radians
SUBSTEPS_PER_DAY = 24
DT_MINUTES = 1440.0 / SUBSTEPS_PER_DAY  # 60 min
SPEED_MODIFIERS = {S: 1.0, E: 1.0, I1: 0.5, I2: 0.1, D: 0.0, R: 1.0}

# Temperature & Arrhenius
T_CELSIUS = 13.0
T_KELVIN = T_CELSIUS + 273.15
T_REF_K = 293.15  # 20°C reference
EA = 5000.0  # activation energy (simplified)
T_VBNC = 12.0
K_VBNC = 2.0

# Shedding rates (Arrhenius-scaled from 20°C reference)
arrhenius_factor = np.exp(EA * (1.0 / T_REF_K - 1.0 / T_KELVIN))
SIGMA_1 = 5.0 * arrhenius_factor   # I1 shedding (~3.3)
SIGMA_2 = 50.0 * arrhenius_factor  # I2 shedding (~33)
SIGMA_D = 15.0  # fresh dead shedding (not Arrhenius-scaled in model)

# Pathogen dynamics
K_HALF = 87_000.0  # within-site dose-response
A_EXPOSURE = 0.75
DECAY_RATE = 0.1
PHI_K = 0.5  # flushing
P_ENV_FLOOR = 50.0
ALPHA_ENV = 0.18
DELTA_ENV = 0.05

# Disease progression at 13°C (Arrhenius-adjusted)
MU_EI1 = 0.167   # ~6 day mean
MU_I1I2 = 0.285  # ~3.5 day mean
MU_I2D = 0.475   # ~2.1 day mean
RHO_REC = 0.05   # recovery rate coefficient

# Size effect
BETA_L = 0.021
L_BAR = 300.0
SIGMA_L = 100.0

# VBNC activation
VBNC_ACTIVATION = 1.0 / (1.0 + np.exp(-K_VBNC * (T_CELSIUS - T_VBNC)))

# Simulation
N_DAYS = 45
N_INITIAL_INFECTED = 3

print(f"Arrhenius factor (13°C vs 20°C ref): {arrhenius_factor:.4f}")
print(f"SIGMA_1 = {SIGMA_1:.2f}, SIGMA_2 = {SIGMA_2:.2f}")
print(f"VBNC activation = {VBNC_ACTIVATION:.4f}")
print(f"Habitat side = {HABITAT_SIDE:.1f} m, Grid = {NX}x{NY} cells")
print(f"N = {N}, density = {N/HABITAT_AREA:.3f} ind/m²")


# ══════════════════════════════════════════════════════════════
# AGENT INITIALIZATION
# ══════════════════════════════════════════════════════════════

x = np.random.uniform(0, HABITAT_SIDE, N)
y = np.random.uniform(0, HABITAT_SIDE, N)
heading = np.random.uniform(0, 2 * np.pi, N)
size_mm = np.clip(np.random.normal(300, 100, N), 50, 600)
age = np.random.uniform(1, 12, N)
resistance = np.clip(np.random.normal(0.15, 0.05, N), 0, 1)
tolerance = np.clip(np.random.normal(0.10, 0.05, N), 0, 1)
recovery_ability = np.clip(np.random.normal(0.02, 0.01, N), 0, 1)
disease_state = np.full(N, S, dtype=int)
alive = np.ones(N, dtype=bool)

# Place 3 initial I1 near center
center = HABITAT_SIDE / 2.0
initial_infected_idx = np.argsort((x - center)**2 + (y - center)**2)[:N_INITIAL_INFECTED]
disease_state[initial_infected_idx] = I1
print(f"Initial I1 agents: {initial_infected_idx}")

# Track 3 agents for movement trails
tracked_agents = [initial_infected_idx[0], initial_infected_idx[1],
                  np.random.choice([i for i in range(N) if i not in initial_infected_idx])]
trail_history = {tid: [] for tid in tracked_agents}

# ══════════════════════════════════════════════════════════════
# BOUNDARY REFLECTION (vectorized)
# ══════════════════════════════════════════════════════════════

def reflect_boundary(coord):
    """Elastic reflection via modular arithmetic."""
    period = 2.0 * HABITAT_SIDE
    coord = coord % period
    mask = coord > HABITAT_SIDE
    coord[mask] = period - coord[mask]
    return coord


# ══════════════════════════════════════════════════════════════
# INFECTED DENSITY GRID
# ══════════════════════════════════════════════════════════════

def compute_exposure_grid(x, y, disease_state, alive):
    """
    Bin I1 and I2 into grid, apply 2 passes of 3x3 averaging with
    reflect padding, normalize so mean = 1.0.
    """
    grid = np.zeros((NY, NX), dtype=float)
    infected_mask = alive & ((disease_state == I1) | (disease_state == I2))
    
    if not np.any(infected_mask):
        return np.ones((NY, NX))
    
    ix = np.clip((x[infected_mask] / CELL_SIZE).astype(int), 0, NX - 1)
    iy = np.clip((y[infected_mask] / CELL_SIZE).astype(int), 0, NY - 1)
    
    # Vectorized binning using np.add.at
    np.add.at(grid, (iy, ix), 1.0)
    
    # 2 passes of 3x3 averaging with reflect padding
    for _ in range(2):
        grid = uniform_filter(grid, size=3, mode='reflect')
    
    # Normalize so mean = 1.0
    grid_mean = grid.mean()
    if grid_mean > 0:
        grid = grid / grid_mean
    else:
        grid = np.ones((NY, NX))
    
    return grid


def get_exposure_factor(x_arr, y_arr, grid):
    """Look up grid cell for each agent."""
    ix = np.clip((x_arr / CELL_SIZE).astype(int), 0, NX - 1)
    iy = np.clip((y_arr / CELL_SIZE).astype(int), 0, NY - 1)
    return grid[iy, ix]


# ══════════════════════════════════════════════════════════════
# SIMULATION LOOP
# ══════════════════════════════════════════════════════════════

# State tracking
P_k = 100.0  # initial environmental Vibrio
P_env_pool = 0.0
daily_counts = []
daily_P_k = []

# Snapshot storage: (day, x_copy, y_copy, state_copy, alive_copy, trails_copy)
snapshots = []


def take_snapshot(day):
    trail_snap = {tid: list(trail_history[tid]) for tid in tracked_agents}
    snapshots.append((day, x.copy(), y.copy(), disease_state.copy(),
                       alive.copy(), trail_snap))


take_snapshot(0)

# Pre-compute speed modifier array for vectorized movement
speed_mod_arr = np.zeros(6)
for st_key, mod_val in SPEED_MODIFIERS.items():
    speed_mod_arr[st_key] = mod_val

for day in range(1, N_DAYS + 1):
    # ── Movement (24 substeps per day) — fully vectorized ──
    for substep in range(SUBSTEPS_PER_DAY):
        # Compute per-agent speed modifier from disease state
        agent_speed_mod = speed_mod_arr[disease_state]
        # Only move alive agents with non-zero speed
        move_mask = alive & (agent_speed_mod > 0)
        n_move = np.sum(move_mask)
        if n_move > 0:
            heading[move_mask] += np.random.normal(0, SIGMA_TURN, n_move)
            speeds = BASE_SPEED * agent_speed_mod[move_mask]
            x[move_mask] += speeds * np.cos(heading[move_mask]) * DT_MINUTES
            y[move_mask] += speeds * np.sin(heading[move_mask]) * DT_MINUTES
        
        # Boundary reflection
        x = reflect_boundary(x)
        y = reflect_boundary(y)
    
    # Record trails (once per day after movement)
    for tid in tracked_agents:
        trail_history[tid].append((x[tid], y[tid]))
        # Keep last 5 days only
        if len(trail_history[tid]) > 5:
            trail_history[tid] = trail_history[tid][-5:]
    
    # ── Spatial exposure grid ──
    exposure_grid = compute_exposure_grid(x, y, disease_state, alive)
    
    # ── Count infected for shedding ──
    n_I1 = np.sum(alive & (disease_state == I1))
    n_I2 = np.sum(alive & (disease_state == I2))
    n_D_fresh = 0  # updated after transitions
    
    # ── Environmental pathogen dynamics ──
    shed = SIGMA_1 * n_I1 + SIGMA_2 * n_I2 + SIGMA_D * n_D_fresh
    P_env_pool += ALPHA_ENV * shed - DELTA_ENV * P_env_pool
    env_input = P_ENV_FLOOR * VBNC_ACTIVATION + P_env_pool
    dP = shed - DECAY_RATE * P_k - PHI_K * P_k + env_input
    P_k = max(0.0, P_k + dP)
    
    # ── Force of infection (transmission) — vectorized ──
    susceptible_mask = alive & (disease_state == S)
    if np.any(susceptible_mask) and P_k > 0:
        s_idx = np.where(susceptible_mask)[0]
        local_exposure = get_exposure_factor(x[s_idx], y[s_idx], exposure_grid)
        P_local = P_k * local_exposure
        dose_response = P_local / (K_HALF + P_local)
        resistance_mod = 1.0 - resistance[s_idx]
        size_mod = np.exp(BETA_L * (size_mm[s_idx] - L_BAR) / SIGMA_L)
        lambda_i = A_EXPOSURE * dose_response * resistance_mod * size_mod
        p_infection = 1.0 - np.exp(-lambda_i)
        
        # Stochastic infection
        rolls = np.random.random(len(s_idx))
        new_exposed = s_idx[rolls < p_infection]
        disease_state[new_exposed] = E
    
    # ── Disease progression (stochastic daily transitions) ──
    n_D_fresh = 0
    
    # E → I1
    e_mask = alive & (disease_state == E)
    if np.any(e_mask):
        e_idx = np.where(e_mask)[0]
        p_ei1 = 1.0 - np.exp(-MU_EI1)
        rolls = np.random.random(len(e_idx))
        disease_state[e_idx[rolls < p_ei1]] = I1
    
    # I1 → I2
    i1_mask = alive & (disease_state == I1)
    if np.any(i1_mask):
        i1_idx = np.where(i1_mask)[0]
        p_i1i2 = 1.0 - np.exp(-MU_I1I2)
        rolls = np.random.random(len(i1_idx))
        disease_state[i1_idx[rolls < p_i1i2]] = I2
    
    # I2 → D or R
    i2_mask = alive & (disease_state == I2)
    if np.any(i2_mask):
        i2_idx = np.where(i2_mask)[0]
        tol = tolerance[i2_idx]
        p_i2d = 1.0 - np.exp(-MU_I2D * (1.0 - 0.85 * tol))
        rolls_d = np.random.random(len(i2_idx))
        
        # Recovery (only if recovery_ability > 0.5, very rare)
        rec = recovery_ability[i2_idx]
        p_i2r = RHO_REC * rec
        p_i2r[rec <= 0.5] = 0.0
        rolls_r = np.random.random(len(i2_idx))
        
        # Apply recovery first (takes priority)
        recover_mask = rolls_r < p_i2r
        disease_state[i2_idx[recover_mask]] = R
        
        # Then death (only for those not recovering)
        death_mask = (~recover_mask) & (rolls_d < p_i2d)
        dying_idx = i2_idx[death_mask]
        disease_state[dying_idx] = D
        alive[dying_idx] = False
        n_D_fresh = len(dying_idx)
    
    # Update P_env with fresh dead shedding
    if n_D_fresh > 0:
        extra_shed = SIGMA_D * n_D_fresh
        P_env_pool += ALPHA_ENV * extra_shed
        P_k = max(0.0, P_k + extra_shed)
    
    # ── Record daily counts ──
    counts = {
        'S': int(np.sum(disease_state == S)),
        'E': int(np.sum(disease_state == E)),
        'I1': int(np.sum(disease_state == I1)),
        'I2': int(np.sum(disease_state == I2)),
        'D': int(np.sum(disease_state == D)),
        'R': int(np.sum(disease_state == R)),
    }
    daily_counts.append(counts)
    daily_P_k.append(P_k)
    
    take_snapshot(day)
    
    if day % 5 == 0 or day <= 5:
        print(f"Day {day:3d}: S={counts['S']:4d} E={counts['E']:4d} "
              f"I1={counts['I1']:4d} I2={counts['I2']:4d} "
              f"D={counts['D']:4d} R={counts['R']:4d} | P_k={P_k:.0f}")


# ══════════════════════════════════════════════════════════════
# ADAPTIVE TIMEPOINT SELECTION
# ══════════════════════════════════════════════════════════════

active_infected = [c['I1'] + c['I2'] for c in daily_counts]
total_infected = [c['E'] + c['I1'] + c['I2'] for c in daily_counts]
total_dead = [c['D'] for c in daily_counts]

peak_active_day = int(np.argmax(active_infected)) + 1  # 1-indexed
peak_total_day = int(np.argmax(total_infected)) + 1

print("\n── Summary ──")
print(f"Final: S={daily_counts[-1]['S']} E={daily_counts[-1]['E']} "
      f"I1={daily_counts[-1]['I1']} I2={daily_counts[-1]['I2']} "
      f"D={daily_counts[-1]['D']} R={daily_counts[-1]['R']}")
print(f"Peak active infected (I1+I2): {max(active_infected)} on day {peak_active_day}")
print(f"Peak total infected (E+I1+I2): {max(total_infected)} on day {peak_total_day}")
print(f"Final P_k: {daily_P_k[-1]:.0f}")

# Adaptive timepoints:
# 0: initial state
# 1: first visible cluster (~peak_active/4, early spread)
# 2: wave front visible (~peak_active/2)
# 3: peak infection
# 4: post-peak (~halfway between peak and end)
# 5: final state
t0 = 0
t1 = max(2, peak_active_day // 4)        # early spread
t2 = max(t1 + 2, peak_active_day // 2)   # wave front
t3 = peak_active_day                      # peak
t4 = min(N_DAYS, peak_active_day + (N_DAYS - peak_active_day) // 2)  # post-peak
t5 = N_DAYS                              # final

# Ensure strictly increasing and within bounds
timepoints = sorted(set([t0, t1, t2, t3, t4, t5]))
# If we have duplicates or too few, space them out
if len(timepoints) < 6:
    # Fill in with evenly spaced
    timepoints = [0]
    for i in range(1, 5):
        timepoints.append(int(round(i * N_DAYS / 5)))
    timepoints.append(N_DAYS)

timepoints = timepoints[:6]
print(f"Selected timepoints: {timepoints}")


# ══════════════════════════════════════════════════════════════
# FIGURE
# ══════════════════════════════════════════════════════════════

STATE_COLORS = {
    S: 'steelblue',
    E: 'gold',
    I1: 'orange',
    I2: 'crimson',
    D: 'gray',
    R: 'limegreen',
}

TRAIL_COLORS = ['#e74c3c', '#2ecc71', '#3498db']

plt.style.use('seaborn-v0_8-whitegrid')
fig, axes = plt.subplots(2, 3, figsize=(16, 11), dpi=150)
axes = axes.flatten()

for panel_idx, day in enumerate(timepoints):
    ax = axes[panel_idx]
    snap = snapshots[day]
    snap_day, snap_x, snap_y, snap_state, snap_alive, snap_trails = snap
    
    # Grid overlay (very faint for dense plots)
    for gx in np.arange(0, HABITAT_SIDE + CELL_SIZE, CELL_SIZE):
        ax.axvline(gx, color='lightgray', linewidth=0.3, alpha=0.1)
    for gy in np.arange(0, HABITAT_SIDE + CELL_SIZE, CELL_SIZE):
        ax.axhline(gy, color='lightgray', linewidth=0.3, alpha=0.1)
    
    # Movement trails for tracked agents
    for t_idx, tid in enumerate(tracked_agents):
        trail = snap_trails[tid]
        if len(trail) >= 2:
            tx = [p[0] for p in trail]
            ty = [p[1] for p in trail]
            seg_x, seg_y = [tx[0]], [ty[0]]
            for j in range(1, len(tx)):
                if abs(tx[j] - tx[j-1]) > HABITAT_SIDE * 0.4 or \
                   abs(ty[j] - ty[j-1]) > HABITAT_SIDE * 0.4:
                    if len(seg_x) >= 2:
                        ax.plot(seg_x, seg_y, color=TRAIL_COLORS[t_idx],
                                linewidth=0.3, alpha=0.4, zorder=2)
                    seg_x, seg_y = [tx[j]], [ty[j]]
                else:
                    seg_x.append(tx[j])
                    seg_y.append(ty[j])
            if len(seg_x) >= 2:
                ax.plot(seg_x, seg_y, color=TRAIL_COLORS[t_idx],
                        linewidth=0.3, alpha=0.4, zorder=2)
    
    # Plot agents by state (dead first, then alive)
    # Dead agents — low alpha since there will be many
    dead_mask = snap_state == D
    if np.any(dead_mask):
        ax.scatter(snap_x[dead_mask], snap_y[dead_mask], s=1,
                   c=STATE_COLORS[D], alpha=0.15, zorder=3, edgecolors='none')
    
    # Alive agents by state
    for st in [S, E, I1, I2, R]:
        mask = (snap_state == st) & snap_alive
        if np.any(mask):
            ax.scatter(snap_x[mask], snap_y[mask], s=3,
                       c=STATE_COLORS[st], alpha=0.85, zorder=4,
                       edgecolors='none')
    
    ax.set_xlim(-2, HABITAT_SIDE + 2)
    ax.set_ylim(-2, HABITAT_SIDE + 2)
    ax.set_aspect('equal')
    ax.set_title(f'Day {snap_day}', fontsize=13, fontweight='bold')
    
    # Count annotation
    counts = {st: int(np.sum(snap_state == st)) for st in [S, E, I1, I2, D, R]}
    n_infected = counts[I1] + counts[I2]
    count_text = (f"S:{counts[S]}  E:{counts[E]}  "
                  f"I:{n_infected}  D:{counts[D]}")
    if counts[R] > 0:
        count_text += f"  R:{counts[R]}"
    ax.text(0.02, 0.97, count_text, transform=ax.transAxes,
            fontsize=8, verticalalignment='top',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      alpha=0.85, edgecolor='gray', linewidth=0.5),
            fontfamily='monospace')
    
    # Axis labels
    if panel_idx >= 3:
        ax.set_xlabel('x (meters)', fontsize=9)
    if panel_idx % 3 == 0:
        ax.set_ylabel('y (meters)', fontsize=9)
    ax.tick_params(labelsize=7)

# Legend
legend_elements = [
    mpatches.Patch(facecolor='steelblue', label='Susceptible (S)'),
    mpatches.Patch(facecolor='gold', label='Exposed (E)'),
    mpatches.Patch(facecolor='orange', label='Pre-symptomatic (I1)'),
    mpatches.Patch(facecolor='crimson', label='Symptomatic (I2)'),
    mpatches.Patch(facecolor='gray', alpha=0.4, label='Dead (D)'),
    mpatches.Patch(facecolor='limegreen', label='Recovered (R)'),
    Line2D([0], [0], color='gray', linewidth=0.8, linestyle='-',
           alpha=0.5, label='Movement trail'),
]
fig.legend(handles=legend_elements, loc='lower center', ncol=7,
           fontsize=9, frameon=True, fancybox=True,
           edgecolor='gray', framealpha=0.9)

# Title
fig.suptitle('Within-Site Disease Dynamics (K = 5,000, Faithful Parameters)',
             fontsize=16, fontweight='bold', y=0.98)
fig.text(0.5, 0.945,
         f'T = {T_CELSIUS}°C,  K_half = {K_HALF:,.0f},  '
         f'a = {A_EXPOSURE},  base speed = {BASE_SPEED} m/min,  '
         f'habitat = {HABITAT_AREA:,.0f} m²',
         ha='center', fontsize=10, color='0.4', style='italic')

plt.tight_layout(rect=[0, 0.05, 1, 0.93])

outpath = '/home/starbot/.openclaw/workspace/sswd-evoepi/reports/progress_report/figures/fig_within_site_disease.png'
fig.savefig(outpath, dpi=150, bbox_inches='tight', facecolor='white')
plt.close()
print(f"\nFigure saved to {outpath}")
