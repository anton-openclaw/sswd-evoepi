"""Spawning dynamics visualizations + animated GIFs for SSWD-EvoEpi.

Phase 13: Deep visualization of the spawning overhaul (Phase 12).
Shows mass spawning events, spatial aggregation, cascade induction,
readiness waves, and density-dependent Allee effects.

Static plots (Part A):
  1. spawning_event_profile — daily spawning by sex stacked bar
  2. spawning_participation — per-node per-year bar chart
  3. readiness_cascade — timeline of readiness → trigger → cascade
  4. spawning_density_scatter — peak fraction vs density
  5. spawning_before_after — old vs new params comparison

Animated GIFs (Part B):
  1. Mass spawning event at high density
  2. Same event at low density (post-crash)
  3. Parameter comparison (old vs new)
  4. Full season overview

Author: Anton 🔬
"""

from __future__ import annotations

import matplotlib
matplotlib.use('Agg')

from pathlib import Path
from typing import Optional, List, Dict, Any, Tuple
from dataclasses import dataclass

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.animation as animation
import numpy as np

from sswd_evoepi.viz.style import (
    ACCENT_COLORS,
    DARK_BG,
    DARK_PANEL,
    GRID_COLOR,
    NODE_COLORS,
    TEXT_COLOR,
    apply_dark_theme,
    dark_figure,
    save_figure,
)

DAYS_PER_YEAR = 365

# ═══════════════════════════════════════════════════════════════════════
# AGENT STATE SNAPSHOT FOR ANIMATIONS
# ═══════════════════════════════════════════════════════════════════════

@dataclass
class DailySnapshot:
    """Captured agent state for one day (used by animation frames)."""
    day: int
    doy: int          # day of year (1-365)
    x: np.ndarray     # positions (alive only)
    y: np.ndarray
    sex: np.ndarray
    size: np.ndarray
    spawning_ready: np.ndarray
    has_spawned: np.ndarray
    last_spawn_day: np.ndarray
    spawn_gravity_timer: np.ndarray
    alive_count: int
    adult_count: int
    n_spawning_today: int
    n_ready: int
    n_recently_spawned: int  # spawned within last 3 days
    n_done: int              # reached max bouts


# ═══════════════════════════════════════════════════════════════════════
# SIMULATION RUNNER FOR FRAME CAPTURE
# ═══════════════════════════════════════════════════════════════════════

def run_simulation_with_snapshots(
    n_individuals: int = 500,
    carrying_capacity: int = 500,
    habitat_area: float = 10000.0,
    T_celsius: float = 10.5,
    n_years: int = 3,
    seed: int = 42,
    snapshot_start_doy: Optional[int] = None,
    snapshot_end_doy: Optional[int] = None,
    snapshot_year: int = 1,
    spawning_overrides: Optional[Dict[str, Any]] = None,
    disease_year: Optional[int] = None,
) -> Tuple[List[DailySnapshot], 'CoupledSimResult']:
    """Run a coupled simulation capturing daily agent snapshots.

    This is a specialized runner that captures full agent state each day
    during a specified window for animation purposes.

    Args:
        n_individuals: Starting population.
        carrying_capacity: K.
        habitat_area: Habitat area (m²).
        T_celsius: SST.
        n_years: Total years.
        seed: RNG seed.
        snapshot_start_doy: Start capturing at this DOY (1-365). None = entire season.
        snapshot_end_doy: Stop capturing at this DOY. None = entire season.
        snapshot_year: Which year to capture (0-indexed).
        spawning_overrides: Dict of SpawningSection field overrides.
        disease_year: Year to introduce disease. None = no disease.

    Returns:
        (snapshots, result) — list of DailySnapshot and CoupledSimResult.
    """
    from sswd_evoepi.config import default_config
    from sswd_evoepi.model import (
        initialize_population,
        annual_natural_mortality, annual_growth_and_aging,
        make_effect_sizes, settle_daily_cohorts,
        CoupledSimResult,
    )
    from sswd_evoepi.perf import PerfMonitor
    from sswd_evoepi.disease import NodeDiseaseState, daily_disease_update
    from sswd_evoepi.types import DiseaseState, Stage
    from sswd_evoepi.spawning import (
        spawning_step, in_spawning_season, reset_spawning_season,
    )
    from sswd_evoepi.movement import daily_movement

    config = default_config()

    # Apply spawning overrides
    if spawning_overrides:
        for key, val in spawning_overrides.items():
            if hasattr(config.spawning, key):
                setattr(config.spawning, key, val)

    # Disable disease if not wanted
    if disease_year is None:
        disease_year_actual = n_years + 10  # never triggers
    else:
        disease_year_actual = disease_year

    pop_cfg = config.population
    dis_cfg = config.disease
    spawning_cfg = config.spawning

    rng = np.random.default_rng(seed)
    effect_sizes = make_effect_sizes(config.genetics.effect_size_seed)

    max_agents = max(int(carrying_capacity * 2.5), n_individuals * 3)
    agents, genotypes = initialize_population(
        n_individuals=n_individuals,
        max_agents=max_agents,
        habitat_area=habitat_area,
        effect_sizes=effect_sizes,
        pop_cfg=pop_cfg,
        rng=rng,
        genetics_cfg=config.genetics,
    )

    node_disease = NodeDiseaseState(node_id=0)
    disease_active = False
    cumulative_disease_deaths = 0
    accumulated_cohorts = []
    previous_in_season = False
    hab_side = np.sqrt(max(habitat_area, 1.0))

    # Default snapshot window: entire spawning season
    if snapshot_start_doy is None:
        snapshot_start_doy = spawning_cfg.season_start_doy
    if snapshot_end_doy is None:
        snapshot_end_doy = spawning_cfg.season_end_doy

    snapshots = []

    # Tracking arrays for result
    total_days = n_years * DAYS_PER_YEAR
    daily_spawning_counts = np.zeros(total_days, dtype=np.int32)
    yearly_pop = np.zeros(n_years, dtype=np.int32)

    perf = PerfMonitor(enabled=False)

    for year in range(n_years):
        # Seed disease (simplified: just infect a few random agents)
        if year == disease_year_actual:
            alive_idx = np.where(agents['alive'])[0]
            n_infect = min(5, len(alive_idx))
            if n_infect > 0:
                infect_idx = rng.choice(alive_idx, size=n_infect, replace=False)
                agents['disease_state'][infect_idx] = DiseaseState.E
                agents['disease_timer'][infect_idx] = 7
            disease_active = True

        for day in range(DAYS_PER_YEAR):
            sim_day = year * DAYS_PER_YEAR + day
            doy = day + 1  # 1-based

            today_T = T_celsius

            # Movement (1 substep per day for speed; enough for viz)
            daily_movement(
                agents, hab_side,
                sigma_turn=0.5,
                base_speed=0.5,
                substeps=4,  # reduced for speed
                rng=rng,
                spawning_gravity_enabled=spawning_cfg.gravity_enabled,
                gravity_strength=spawning_cfg.gravity_strength,
                gravity_range=spawning_cfg.gravity_range,
                pre_spawn_gravity_days=spawning_cfg.pre_spawn_gravity_days,
                post_spawn_gravity_days=spawning_cfg.post_spawn_gravity_days,
            )

            # Disease step
            if disease_active:
                node_disease = daily_disease_update(
                    agents=agents,
                    node_state=node_disease,
                    T_celsius=today_T,
                    salinity=30.0,
                    phi_k=0.02,
                    dispersal_input=0.0,
                    day=sim_day,
                    cfg=dis_cfg,
                    rng=rng,
                )

            # Settlement
            if accumulated_cohorts:
                ready = [c for c in accumulated_cohorts
                         if (sim_day - c.spawn_day) >= c.pld_days]
                accumulated_cohorts = [c for c in accumulated_cohorts
                                       if (sim_day - c.spawn_day) < c.pld_days]
                if ready:
                    settle_daily_cohorts(
                        ready, agents, genotypes, carrying_capacity,
                        pop_cfg, rng, effect_sizes,
                        habitat_area=habitat_area,
                        sim_day=sim_day,
                    )

            # Spawning
            currently_in_season = in_spawning_season(
                doy, spawning_cfg.season_start_doy, spawning_cfg.season_end_doy
            )
            if previous_in_season and not currently_in_season:
                reset_spawning_season(agents)

            if currently_in_season:
                cohorts_today = spawning_step(
                    agents=agents,
                    genotypes=genotypes,
                    day_of_year=doy,
                    node_latitude=48.0,
                    spawning_config=spawning_cfg,
                    disease_config=dis_cfg,
                    rng=rng,
                    current_sim_day=sim_day,
                    current_sst=today_T,
                )
                if cohorts_today:
                    accumulated_cohorts.extend(cohorts_today)

                alive_spawned = (
                    agents['alive'] &
                    (agents['last_spawn_day'] == doy)
                )
                daily_spawning_counts[sim_day] = int(np.sum(alive_spawned))

            previous_in_season = currently_in_season

            # Capture snapshot if in window
            if year == snapshot_year:
                in_window = _doy_in_window(doy, snapshot_start_doy, snapshot_end_doy)
                if in_window:
                    snap = _capture_snapshot(
                        agents, day=sim_day, doy=doy,
                        spawning_cfg=spawning_cfg,
                    )
                    snapshots.append(snap)

        # Annual updates
        annual_natural_mortality(agents, pop_cfg, rng)
        annual_growth_and_aging(agents, pop_cfg, rng)
        yearly_pop[year] = int(np.sum(agents['alive']))

    # Build minimal result
    result = CoupledSimResult(
        n_years=n_years,
        yearly_pop=yearly_pop,
        yearly_adults=np.zeros(n_years, dtype=np.int32),
        yearly_recruits=np.zeros(n_years, dtype=np.int32),
        yearly_natural_deaths=np.zeros(n_years, dtype=np.int32),
        yearly_disease_deaths=np.zeros(n_years, dtype=np.int32),
        yearly_senescence_deaths=np.zeros(n_years, dtype=np.int32),
        yearly_mean_resistance=np.zeros(n_years, dtype=np.float64),
        yearly_fert_success=np.zeros(n_years, dtype=np.float64),
        daily_spawning_counts=daily_spawning_counts,
    )

    return snapshots, result


def _doy_in_window(doy: int, start_doy: int, end_doy: int) -> bool:
    """Check if DOY is within snapshot window (handles year wrapping)."""
    if start_doy <= end_doy:
        return start_doy <= doy <= end_doy
    else:
        return doy >= start_doy or doy <= end_doy


def _capture_snapshot(
    agents: np.ndarray,
    day: int,
    doy: int,
    spawning_cfg,
) -> DailySnapshot:
    """Capture current agent state as a DailySnapshot."""
    from sswd_evoepi.types import Stage

    alive = agents['alive'].astype(bool)
    adult = alive & (agents['stage'] == Stage.ADULT)

    # Count spawning categories
    n_spawning_today = int(np.sum(alive & (agents['last_spawn_day'] == doy)))
    n_ready = int(np.sum(alive & (agents['spawning_ready'] == 1)))

    # Recently spawned (within last 3 days, handling wrapping)
    last_spawn = agents['last_spawn_day'][alive]
    days_since = np.where(
        last_spawn == 0, 999,
        np.where(last_spawn <= doy, doy - last_spawn, (365 - last_spawn) + doy)
    )
    n_recently_spawned = int(np.sum((days_since > 0) & (days_since <= 3)))

    # Done for season
    females_done = alive & (agents['sex'] == 0) & (agents['has_spawned'] >= spawning_cfg.female_max_bouts)
    males_done = alive & (agents['sex'] == 1) & (agents['has_spawned'] >= spawning_cfg.male_max_bouts)
    n_done = int(np.sum(females_done | males_done))

    return DailySnapshot(
        day=day,
        doy=doy,
        x=agents['x'][alive].copy(),
        y=agents['y'][alive].copy(),
        sex=agents['sex'][alive].copy(),
        size=agents['size'][alive].copy(),
        spawning_ready=agents['spawning_ready'][alive].copy(),
        has_spawned=agents['has_spawned'][alive].copy(),
        last_spawn_day=agents['last_spawn_day'][alive].copy(),
        spawn_gravity_timer=agents['spawn_gravity_timer'][alive].copy(),
        alive_count=int(np.sum(alive)),
        adult_count=int(np.sum(adult)),
        n_spawning_today=n_spawning_today,
        n_ready=n_ready,
        n_recently_spawned=n_recently_spawned,
        n_done=n_done,
    )


# ═══════════════════════════════════════════════════════════════════════
# AGENT STATE CLASSIFICATION FOR ANIMATION
# ═══════════════════════════════════════════════════════════════════════

def classify_agents(snap: DailySnapshot, max_female_bouts: int = 2, max_male_bouts: int = 3) -> np.ndarray:
    """Classify each agent into a visual category.

    Categories (priority order, highest first):
      4 = spawning TODAY (orange)
      3 = recently spawned within 3 days (red)
      2 = done for season (green)
      1 = ready to spawn (blue)
      0 = not ready (gray)

    Returns:
        Integer array of categories, same length as snap.x.
    """
    n = len(snap.x)
    cats = np.zeros(n, dtype=np.int8)

    # 1: ready
    cats[snap.spawning_ready == 1] = 1

    # 2: done
    female_done = (snap.sex == 0) & (snap.has_spawned >= max_female_bouts)
    male_done = (snap.sex == 1) & (snap.has_spawned >= max_male_bouts)
    cats[female_done | male_done] = 2

    # 3: recently spawned (within 3 days, not today)
    days_since = np.where(
        snap.last_spawn_day == 0, 999,
        np.where(
            snap.last_spawn_day <= snap.doy,
            snap.doy - snap.last_spawn_day,
            (365 - snap.last_spawn_day) + snap.doy
        )
    )
    recently_spawned = (days_since > 0) & (days_since <= 3)
    cats[recently_spawned] = 3

    # 4: spawning today (overrides everything)
    spawning_today = snap.last_spawn_day == snap.doy
    cats[spawning_today] = 4

    return cats


# Category colors (dark background friendly)
CAT_COLORS = {
    0: '#555555',   # gray — not ready
    1: '#3498db',   # blue — ready
    2: '#2ecc71',   # green — done
    3: '#e74c3c',   # red — recently spawned
    4: '#f39c12',   # orange — spawning TODAY
}

CAT_LABELS = {
    0: 'Not ready',
    1: 'Ready',
    2: 'Done (season)',
    3: 'Recently spawned',
    4: 'Spawning TODAY',
}


# ═══════════════════════════════════════════════════════════════════════
# PART A: STATIC VISUALIZATIONS
# ═══════════════════════════════════════════════════════════════════════

# A1. Mass Spawning Event Profile
def plot_spawning_event_profile(
    snapshots: List[DailySnapshot],
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Daily spawning counts stacked by sex, showing mass event bursts.

    Caption: Mass spawning events concentrate reproduction into 1-2 day
    bursts driven by spawning induction cascade. Females (κ_fm=0.80) strongly
    trigger males; readiness induction accelerates the build-up. The sharp
    peak demonstrates that the cascade mechanism creates biologically
    realistic synchronous broadcast spawning from individual-level rules.
    """
    fig, ax = dark_figure(figsize=(14, 6))

    doys = [s.doy for s in snapshots]
    female_counts = []
    male_counts = []

    for snap in snapshots:
        spawning_today = snap.last_spawn_day == snap.doy
        female_spawning = spawning_today & (snap.sex == 0)
        male_spawning = spawning_today & (snap.sex == 1)
        female_counts.append(int(np.sum(female_spawning)))
        male_counts.append(int(np.sum(male_spawning)))

    female_arr = np.array(female_counts)
    male_arr = np.array(male_counts)

    ax.bar(doys, female_arr, color='#e94560', alpha=0.85, label='Females', width=1.0)
    ax.bar(doys, male_arr, bottom=female_arr, color='#3498db', alpha=0.85, label='Males', width=1.0)

    ax.set_xlabel('Day of Year', fontsize=12, color=TEXT_COLOR)
    ax.set_ylabel('Number Spawning', fontsize=12, color=TEXT_COLOR)
    ax.set_title('Mass Spawning Event Profile (Year 2, Pre-Disease)',
                 fontsize=14, fontweight='bold', color=TEXT_COLOR)
    ax.legend(facecolor=DARK_PANEL, edgecolor=GRID_COLOR,
              labelcolor=TEXT_COLOR, fontsize=10)

    # Annotate peak
    total = female_arr + male_arr
    peak_idx = np.argmax(total)
    if total[peak_idx] > 0:
        ax.annotate(
            f'Peak: {total[peak_idx]} spawners\n(DOY {doys[peak_idx]})',
            xy=(doys[peak_idx], total[peak_idx]),
            xytext=(15, 15), textcoords='offset points',
            fontsize=10, color='white', fontweight='bold',
            arrowprops=dict(arrowstyle='->', color='white', lw=1.5),
        )

    caption = (
        "Stacked daily spawning by sex. Mass events driven by cascade induction:\n"
        "females (κ_fm=0.80) strongly trigger males; readiness induction accelerates build-up.\n"
        "Sharp peaks = biologically realistic broadcast spawning from individual-level rules."
    )
    fig.text(0.5, -0.06, caption, ha='center', va='top', fontsize=9,
             color=TEXT_COLOR, style='italic')

    if save_path:
        save_figure(fig, save_path)
    return fig


# A2. Spawning Participation Rate
def plot_spawning_participation(
    snapshots_by_year: Dict[int, List[DailySnapshot]],
    n_individuals: int = 500,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Bar chart of fraction that spawned per year.

    Caption: Spawning participation rate shows the fraction of the
    adult population that successfully spawned at least once per season.
    At high density (pre-epidemic), participation is high due to effective
    cascade propagation. Post-epidemic, reduced density weakens the
    cascade, lowering participation — an emergent Allee effect on
    reproduction that compounds population decline.
    """
    fig, ax = dark_figure(figsize=(10, 6))

    years = sorted(snapshots_by_year.keys())
    participation_rates = []

    for yr in years:
        snaps = snapshots_by_year[yr]
        if not snaps:
            participation_rates.append(0.0)
            continue
        # Use end-of-season snapshot
        final_snap = snaps[-1]
        n_adults = final_snap.adult_count
        n_spawned = int(np.sum(final_snap.has_spawned > 0))
        rate = n_spawned / max(n_adults, 1)
        participation_rates.append(rate)

    colors = [ACCENT_COLORS[3] if r > 0.5 else ACCENT_COLORS[0] for r in participation_rates]
    bars = ax.bar(years, participation_rates, color=colors, edgecolor='white',
                  linewidth=0.5, alpha=0.85)

    for bar, rate in zip(bars, participation_rates):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.02,
                f'{rate:.0%}', ha='center', va='bottom', fontsize=11, color=TEXT_COLOR)

    ax.set_xlabel('Simulation Year', fontsize=12, color=TEXT_COLOR)
    ax.set_ylabel('Spawning Participation Rate', fontsize=12, color=TEXT_COLOR)
    ax.set_title('Spawning Participation by Year',
                 fontsize=14, fontweight='bold', color=TEXT_COLOR)
    ax.set_ylim(0, 1.1)
    ax.axhline(0.8, color=GRID_COLOR, linestyle='--', alpha=0.5, label='80% target')
    ax.legend(facecolor=DARK_PANEL, edgecolor=GRID_COLOR, labelcolor=TEXT_COLOR, fontsize=9)

    caption = (
        "Fraction of adults spawning at least once per season.\n"
        "High density → cascade propagation → high participation.\n"
        "Low density → cascade failure → Allee effect on reproduction."
    )
    fig.text(0.5, -0.06, caption, ha='center', va='top', fontsize=9,
             color=TEXT_COLOR, style='italic')

    if save_path:
        save_figure(fig, save_path)
    return fig


# A3. Readiness Induction Cascade Diagram
def plot_readiness_cascade(
    snapshots: List[DailySnapshot],
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Timeline showing readiness wave → trigger → cascade → mass event.

    Caption: The readiness cascade operates in distinct phases:
    (1) Seasonal readiness wave — individuals become reproductively ready
    following a Normal envelope around the peak DOY. (2) Spontaneous triggers —
    a few ready individuals spawn spontaneously. (3) Cascade induction —
    spawners trigger nearby ready individuals via chemical cues
    (κ_fm=0.80, κ_mf=0.60). (4) Mass event — positive feedback creates
    an exponential burst. (5) Exhaustion — pool of ready individuals depleted.
    """
    fig, axes = dark_figure(nrows=3, ncols=1, figsize=(14, 12))
    ax_ready, ax_spawn, ax_done = axes

    doys = [s.doy for s in snapshots]

    # Top: readiness accumulation
    n_ready = [s.n_ready for s in snapshots]
    ax_ready.fill_between(doys, 0, n_ready, color='#3498db', alpha=0.4)
    ax_ready.plot(doys, n_ready, color='#3498db', linewidth=2, label='Ready to spawn')
    ax_ready.set_ylabel('Count', fontsize=11, color=TEXT_COLOR)
    ax_ready.set_title('Phase 1: Readiness Wave', fontsize=13, fontweight='bold', color=TEXT_COLOR)
    ax_ready.legend(facecolor=DARK_PANEL, edgecolor=GRID_COLOR, labelcolor=TEXT_COLOR, fontsize=9)

    # Middle: daily spawning (spontaneous + cascade-induced)
    n_spawning = [s.n_spawning_today for s in snapshots]
    n_recent = [s.n_recently_spawned for s in snapshots]
    ax_spawn.bar(doys, n_spawning, color='#f39c12', alpha=0.85, label='Spawning today', width=1.0)
    ax_spawn.plot(doys, n_recent, color='#e74c3c', linewidth=1.5, alpha=0.7,
                  label='Recently spawned (3d window)')
    ax_spawn.set_ylabel('Count', fontsize=11, color=TEXT_COLOR)
    ax_spawn.set_title('Phase 2-4: Spontaneous → Cascade → Mass Event',
                       fontsize=13, fontweight='bold', color=TEXT_COLOR)
    ax_spawn.legend(facecolor=DARK_PANEL, edgecolor=GRID_COLOR, labelcolor=TEXT_COLOR, fontsize=9)

    # Bottom: cumulative done
    n_done = [s.n_done for s in snapshots]
    ax_done.fill_between(doys, 0, n_done, color='#2ecc71', alpha=0.4)
    ax_done.plot(doys, n_done, color='#2ecc71', linewidth=2, label='Done for season')
    ax_done.set_xlabel('Day of Year', fontsize=12, color=TEXT_COLOR)
    ax_done.set_ylabel('Count', fontsize=11, color=TEXT_COLOR)
    ax_done.set_title('Phase 5: Pool Exhaustion', fontsize=13, fontweight='bold', color=TEXT_COLOR)
    ax_done.legend(facecolor=DARK_PANEL, edgecolor=GRID_COLOR, labelcolor=TEXT_COLOR, fontsize=9)

    fig.suptitle('Spawning Cascade Dynamics: Readiness → Trigger → Burst → Quiet',
                 fontsize=15, fontweight='bold', color=TEXT_COLOR, y=1.01)

    caption = (
        "Three-phase cascade: (1) Seasonal readiness wave builds spawner pool.\n"
        "(2) Spontaneous triggers + cascade induction (κ_fm=0.80, κ_mf=0.60) create exponential burst.\n"
        "(3) Pool exhaustion as females reach single-spawn limit and males hit refractory."
    )
    fig.text(0.5, -0.02, caption, ha='center', va='top', fontsize=9,
             color=TEXT_COLOR, style='italic')

    if save_path:
        save_figure(fig, save_path)
    return fig


# A4. Spawning Intensity vs Density Scatter
def plot_spawning_vs_density(
    snapshots_list: List[Tuple[str, List[DailySnapshot], int]],
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Scatter: X = pop/K, Y = peak daily spawning fraction.

    Args:
        snapshots_list: List of (label, snapshots, K) tuples.

    Caption: Peak daily spawning intensity as a function of population
    density reveals the emergent density-dependent Allee effect.
    At high density (>0.5K), cascade propagation concentrates spawning into
    mass events. At low density, cascade fails — diffuse spawning, lower
    fertilization. This is the key mechanism linking population crash to
    delayed reproductive recovery.
    """
    fig, ax = dark_figure(figsize=(10, 8))

    for i, (label, snaps, K) in enumerate(snapshots_list):
        if not snaps:
            continue
        color = NODE_COLORS[i % len(NODE_COLORS)]

        # Calculate peak spawning fraction per "chunk" of ~30 days
        chunk_size = 30
        densities = []
        peak_fracs = []

        for start in range(0, len(snaps), chunk_size):
            chunk = snaps[start:start+chunk_size]
            if not chunk:
                continue
            avg_pop = np.mean([s.alive_count for s in chunk])
            peak_spawning = max(s.n_spawning_today for s in chunk)
            density = avg_pop / max(K, 1)
            peak_frac = peak_spawning / max(avg_pop, 1)
            densities.append(density)
            peak_fracs.append(peak_frac)

        ax.scatter(densities, peak_fracs, color=color, s=80, alpha=0.7,
                   edgecolors='white', linewidths=0.5, label=label, zorder=3)

    ax.set_xlabel('Population Density (fraction of K)', fontsize=12, color=TEXT_COLOR)
    ax.set_ylabel('Peak Daily Spawning Fraction', fontsize=12, color=TEXT_COLOR)
    ax.set_title('Spawning Intensity vs Population Density\n(Emergent Reproductive Allee Effect)',
                 fontsize=14, fontweight='bold', color=TEXT_COLOR)
    ax.legend(facecolor=DARK_PANEL, edgecolor=GRID_COLOR, labelcolor=TEXT_COLOR, fontsize=9)

    caption = (
        "Each point = 30-day chunk. High density → effective cascade → concentrated spawning.\n"
        "Low density → cascade failure → diffuse spawning → reproductive Allee effect."
    )
    fig.text(0.5, -0.06, caption, ha='center', va='top', fontsize=9,
             color=TEXT_COLOR, style='italic')

    if save_path:
        save_figure(fig, save_path)
    return fig


# A5. Before/After Spawning Overhaul
def plot_spawning_before_after(
    snapshots_old: List[DailySnapshot],
    snapshots_new: List[DailySnapshot],
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Two-panel comparison: old params vs new (Phase 12) params.

    Caption: LEFT: Old parameters (female single-bout, 21-day refractory,
    κ_mf=0.30) produce diffuse, spread-out spawning with weak cascade.
    RIGHT: Phase 12 parameters (female 2-bout, zero refractory, κ_mf=0.60,
    readiness induction) create sharper, tighter mass spawning events.
    The new system produces biologically realistic broadcast spawning bursts
    while maintaining extended seasonal coverage.
    """
    fig, (ax_old, ax_new) = dark_figure(nrows=1, ncols=2, figsize=(16, 6))

    for ax, snaps, title, color in [
        (ax_old, snapshots_old, 'Old Parameters\n(1 bout, 21d refrac, κ_mf=0.30)', ACCENT_COLORS[0]),
        (ax_new, snapshots_new, 'Phase 12 Parameters\n(2 bouts, 0d refrac, κ_mf=0.60)', ACCENT_COLORS[3]),
    ]:
        doys = [s.doy for s in snaps]
        counts = [s.n_spawning_today for s in snaps]
        ax.bar(doys, counts, color=color, alpha=0.85, width=1.0)
        ax.set_xlabel('Day of Year', fontsize=11, color=TEXT_COLOR)
        ax.set_ylabel('Spawners/Day', fontsize=11, color=TEXT_COLOR)
        ax.set_title(title, fontsize=12, fontweight='bold', color=TEXT_COLOR)

        # Stats annotation
        total = sum(counts)
        peak = max(counts) if counts else 0
        n_active_days = sum(1 for c in counts if c > 0)
        ax.text(0.98, 0.95,
                f'Total events: {total}\nPeak: {peak}/day\nActive days: {n_active_days}',
                transform=ax.transAxes, ha='right', va='top', fontsize=9,
                color=TEXT_COLOR, fontfamily='monospace',
                bbox=dict(boxstyle='round', facecolor=DARK_PANEL, edgecolor=GRID_COLOR, alpha=0.8))

    fig.suptitle('Spawning Parameter Comparison: Before vs After Phase 12 Overhaul',
                 fontsize=14, fontweight='bold', color=TEXT_COLOR, y=1.02)

    caption = (
        "LEFT: Old parameters produce diffuse spawning; weak cascade.\n"
        "RIGHT: Phase 12 params create sharper bursts via stronger induction + multi-bout females.\n"
        "Both maintain extended seasonal coverage (~270 days)."
    )
    fig.text(0.5, -0.06, caption, ha='center', va='top', fontsize=9,
             color=TEXT_COLOR, style='italic')

    if save_path:
        save_figure(fig, save_path)
    return fig


# ═══════════════════════════════════════════════════════════════════════
# PART B: ANIMATED GIFS
# ═══════════════════════════════════════════════════════════════════════

def _make_spawning_animation_frame(
    ax: plt.Axes,
    snap: DailySnapshot,
    hab_side: float,
    max_female_bouts: int = 2,
    max_male_bouts: int = 3,
    title_prefix: str = '',
) -> None:
    """Render a single animation frame onto an Axes.

    Clears the axes and redraws all agents with category-based coloring
    on a dark background.
    """
    ax.clear()
    ax.set_facecolor('#0d0d1a')  # very dark blue-black
    ax.set_xlim(0, hab_side)
    ax.set_ylim(0, hab_side)

    if len(snap.x) == 0:
        ax.text(hab_side / 2, hab_side / 2, 'No agents',
                ha='center', va='center', color=TEXT_COLOR, fontsize=14)
        return

    # Classify agents
    cats = classify_agents(snap, max_female_bouts, max_male_bouts)

    # Size: scale by body size (10-800mm → 5-50 points)
    sizes = np.clip(snap.size, 10, 800)
    point_sizes = 5 + (sizes - 10) / (800 - 10) * 45

    # Draw in category order (lowest first, so highest draws on top)
    for cat_val in [0, 1, 2, 3, 4]:
        mask = cats == cat_val
        if not np.any(mask):
            continue
        ax.scatter(
            snap.x[mask], snap.y[mask],
            s=point_sizes[mask],
            c=CAT_COLORS[cat_val],
            alpha=0.8 if cat_val >= 3 else 0.6,
            edgecolors='white' if cat_val == 4 else 'none',
            linewidths=1.0 if cat_val == 4 else 0,
            zorder=cat_val + 2,
        )

    # Frame info panel
    info_text = (
        f"DOY {snap.doy}  |  Pop: {snap.alive_count}  |  "
        f"Spawning: {snap.n_spawning_today}  |  Ready: {snap.n_ready}  |  "
        f"Done: {snap.n_done}"
    )
    ax.set_title(f'{title_prefix}{info_text}',
                 fontsize=10, color=TEXT_COLOR, fontweight='bold', pad=8)

    # Remove axis ticks for cleaner look
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_color('#2a2a4a')


def create_spawning_animation(
    snapshots: List[DailySnapshot],
    hab_side: float = 100.0,
    fps: int = 4,
    title_prefix: str = '',
    max_female_bouts: int = 2,
    max_male_bouts: int = 3,
    save_path: Optional[str] = None,
    dpi: int = 100,
) -> Optional[str]:
    """Create an animated GIF from daily snapshots.

    Args:
        snapshots: List of DailySnapshot objects (one per frame).
        hab_side: Habitat side length (m) for axis limits.
        fps: Frames per second.
        title_prefix: String prefix for frame titles.
        max_female_bouts: For category classification.
        max_male_bouts: For category classification.
        save_path: Path to save GIF. If None, returns None.
        dpi: Resolution.

    Returns:
        Save path if saved, else None.
    """
    if not snapshots:
        print("No snapshots to animate.")
        return None

    fig, ax = plt.subplots(figsize=(10, 10))
    fig.patch.set_facecolor('#0d0d1a')

    # Add legend
    legend_patches = [
        mpatches.Patch(color=CAT_COLORS[0], label=CAT_LABELS[0]),
        mpatches.Patch(color=CAT_COLORS[1], label=CAT_LABELS[1]),
        mpatches.Patch(color=CAT_COLORS[2], label=CAT_LABELS[2]),
        mpatches.Patch(color=CAT_COLORS[3], label=CAT_LABELS[3]),
        mpatches.Patch(color=CAT_COLORS[4], label=CAT_LABELS[4]),
    ]

    def animate(frame_idx):
        snap = snapshots[frame_idx]
        _make_spawning_animation_frame(
            ax, snap, hab_side, max_female_bouts, max_male_bouts, title_prefix
        )
        # Re-add legend after clear
        ax.legend(
            handles=legend_patches,
            loc='lower left',
            fontsize=8,
            facecolor='#1a1a2e',
            edgecolor='#2a2a4a',
            labelcolor=TEXT_COLOR,
            framealpha=0.9,
        )

    anim = animation.FuncAnimation(
        fig, animate, frames=len(snapshots), interval=1000 // fps, blit=False
    )

    if save_path:
        Path(save_path).parent.mkdir(parents=True, exist_ok=True)
        writer = animation.PillowWriter(fps=fps)
        anim.save(save_path, writer=writer, dpi=dpi)
        plt.close(fig)
        print(f"Saved animation: {save_path} ({len(snapshots)} frames)")
        return save_path

    plt.close(fig)
    return None


def create_comparison_animation(
    snapshots_a: List[DailySnapshot],
    snapshots_b: List[DailySnapshot],
    hab_side: float = 100.0,
    fps: int = 4,
    label_a: str = 'Old Parameters',
    label_b: str = 'New Parameters',
    max_female_bouts_a: int = 1,
    max_male_bouts_a: int = 3,
    max_female_bouts_b: int = 2,
    max_male_bouts_b: int = 3,
    save_path: Optional[str] = None,
    dpi: int = 100,
) -> Optional[str]:
    """Create side-by-side comparison animation GIF.

    Shows two simulations with different parameters running in parallel.
    Frames are synced by index (same DOY if simulations align).

    Args:
        snapshots_a: Snapshots from first simulation.
        snapshots_b: Snapshots from second simulation.
        hab_side: Habitat side length.
        fps: Frames per second.
        label_a: Title for left panel.
        label_b: Title for right panel.
        save_path: Output path.
        dpi: Resolution.

    Returns:
        Save path if saved.
    """
    n_frames = min(len(snapshots_a), len(snapshots_b))
    if n_frames == 0:
        print("No frames for comparison animation.")
        return None

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(18, 9))
    fig.patch.set_facecolor('#0d0d1a')

    legend_patches = [
        mpatches.Patch(color=CAT_COLORS[0], label=CAT_LABELS[0]),
        mpatches.Patch(color=CAT_COLORS[1], label=CAT_LABELS[1]),
        mpatches.Patch(color=CAT_COLORS[2], label=CAT_LABELS[2]),
        mpatches.Patch(color=CAT_COLORS[3], label=CAT_LABELS[3]),
        mpatches.Patch(color=CAT_COLORS[4], label=CAT_LABELS[4]),
    ]

    def animate(frame_idx):
        _make_spawning_animation_frame(
            ax_a, snapshots_a[frame_idx], hab_side,
            max_female_bouts_a, max_male_bouts_a, f'{label_a}  |  '
        )
        _make_spawning_animation_frame(
            ax_b, snapshots_b[frame_idx], hab_side,
            max_female_bouts_b, max_male_bouts_b, f'{label_b}  |  '
        )
        ax_b.legend(
            handles=legend_patches, loc='lower left', fontsize=7,
            facecolor='#1a1a2e', edgecolor='#2a2a4a',
            labelcolor=TEXT_COLOR, framealpha=0.9,
        )

    anim = animation.FuncAnimation(
        fig, animate, frames=n_frames, interval=1000 // fps, blit=False
    )

    if save_path:
        Path(save_path).parent.mkdir(parents=True, exist_ok=True)
        writer = animation.PillowWriter(fps=fps)
        anim.save(save_path, writer=writer, dpi=dpi)
        plt.close(fig)
        print(f"Saved comparison animation: {save_path} ({n_frames} frames)")
        return save_path

    plt.close(fig)
    return None
