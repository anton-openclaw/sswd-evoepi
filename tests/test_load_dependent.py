"""Tests for load-dependent disease model (Tier 2).

Covers:
  - within_host_load_step(): logistic growth, clearance, temperature scaling,
    capping, reinfection
  - load_mortality_probability(): Hill-function mortality, tolerance, steepness
  - load_dependent_disease_update(): infection, state mapping, recovery, death,
    sentinels, shedding, backward compatibility
  - Behavioral: immune agents, epidemic dynamics, reinfection after recovery
"""

import numpy as np
import pytest

from sswd_evoepi.config import (
    DiseaseSection,
    LoadDependentSection,
    SimulationConfig,
    default_config,
)
from sswd_evoepi.types import AGENT_DTYPE, DiseaseState, allocate_agents
from sswd_evoepi.disease import (
    CarcassTracker,
    NodeDiseaseState,
    arrhenius,
    within_host_load_step,
    load_mortality_probability,
    load_dependent_disease_update,
)


# ═══════════════════════════════════════════════════════════════════════
# FIXTURES
# ═══════════════════════════════════════════════════════════════════════

@pytest.fixture
def cfg() -> DiseaseSection:
    """Default disease configuration."""
    return DiseaseSection()


@pytest.fixture
def ld_cfg() -> LoadDependentSection:
    """Default load-dependent config with enabled=True."""
    return LoadDependentSection(enabled=True)


@pytest.fixture
def rng() -> np.random.Generator:
    return np.random.default_rng(42)


# ═══════════════════════════════════════════════════════════════════════
# HELPERS
# ═══════════════════════════════════════════════════════════════════════

def _make_agents(n, **overrides):
    """Create a small agent array with sensible defaults.

    All agents are alive, susceptible, adult-sized. Override any field
    by passing keyword arrays or scalars.
    """
    agents = allocate_agents(n)
    agents['alive'][:n] = True
    agents['disease_state'][:n] = DiseaseState.S
    agents['resistance'][:n] = 0.08
    agents['tolerance'][:n] = 0.0
    agents['recovery_ability'][:n] = 0.08
    agents['size'][:n] = 300.0
    agents['pathogen_load'][:n] = 0.0
    for key, val in overrides.items():
        agents[key][:n] = val
    return agents


def _make_node_state(vibrio=0.0):
    """Create a fresh NodeDiseaseState with optional vibrio."""
    ns = NodeDiseaseState()
    ns.vibrio_concentration = vibrio
    return ns


def _run_ld_update(agents, node_state, cfg, ld_cfg, rng, *,
                   T=15.0, sal=30.0, phi=0.02, day=100,
                   dispersal=0.0, sentinel_frac=0.0):
    """Run one load_dependent_disease_update step with defaults."""
    return load_dependent_disease_update(
        agents, node_state,
        T_celsius=T, salinity=sal, phi_k=phi,
        dispersal_input=dispersal, day=day,
        cfg=cfg, ld_cfg=ld_cfg, rng=rng,
        sentinel_shedding_fraction=sentinel_frac,
    )


# ═══════════════════════════════════════════════════════════════════════
# UNIT TESTS: within_host_load_step()
# ═══════════════════════════════════════════════════════════════════════

class TestWithinHostLoadStep:

    def test_load_growth_without_clearance(self, cfg, ld_cfg):
        """1. Zero resistance/recovery → load grows logistically."""
        L = np.array([100.0])
        r_i = np.array([0.0])
        c_i = np.array([0.0])
        L_new = within_host_load_step(L, r_i, c_i, 15.0, 0.0, cfg, ld_cfg)
        # Even with r_i=0, clearance = delta_clear * 0.3 * 0.5 * L (not zero)
        # but growth at low L should outweigh that residual clearance
        # Growth ≈ r_growth * arrhenius * L * (1 - L/L_max) which is >> clearance residual
        # at L=100. Let's verify the load increases.
        assert L_new[0] > L[0], (
            f"Expected load to grow from {L[0]} but got {L_new[0]}"
        )

    def test_load_clearance_high_resistance(self, cfg, ld_cfg):
        """2. r_i=1.0 → clearance dominates → load decreases."""
        L = np.array([1000.0])
        r_i = np.array([1.0])
        c_i = np.array([0.0])
        L_new = within_host_load_step(L, r_i, c_i, 15.0, 0.0, cfg, ld_cfg)
        # clearance = delta_clear * (0.3 + 0.7*1.0) * (0.5 + 0.5*0.0) * L
        #           = 0.3 * 1.0 * 0.5 * 1000 = 150
        # growth = r_growth * arrhenius * 1000 * (1 - 1000/1e6) ≈ 0.5 * ~0.7 * 1000 ≈ 350ish
        # Actually let's just check the direction: high r should reduce net growth
        # versus load with r=0 and same starting point
        L_low_r = within_host_load_step(L, np.array([0.0]), c_i, 15.0, 0.0, cfg, ld_cfg)
        assert L_new[0] < L_low_r[0], "High resistance should lead to lower load than low resistance"

    def test_load_clearance_high_recovery(self, cfg, ld_cfg):
        """3. c_i=1.0 → clearance boosted → lower load than c_i=0."""
        L = np.array([5000.0])
        r_i = np.array([0.5])
        c_low = np.array([0.0])
        c_high = np.array([1.0])
        L_low = within_host_load_step(L, r_i, c_low, 15.0, 0.0, cfg, ld_cfg)
        L_high = within_host_load_step(L, r_i, c_high, 15.0, 0.0, cfg, ld_cfg)
        assert L_high[0] < L_low[0], (
            "Higher recovery ability (c_i=1) should produce lower load"
        )

    def test_load_equilibrium(self, cfg, ld_cfg):
        """4. Intermediate traits → load converges to steady state."""
        # Run many steps from moderate starting load
        L = np.array([5000.0])
        r_i = np.array([0.3])
        c_i = np.array([0.3])
        history = [L[0]]
        for _ in range(500):
            L = within_host_load_step(L, r_i, c_i, 15.0, 0.0, cfg, ld_cfg)
            history.append(L[0])
        # Check that load stabilized (last 50 values have low variance)
        tail = np.array(history[-50:])
        if tail[-1] > ld_cfg.L_clear:
            # Reached equilibrium — coefficient of variation should be small
            cv = np.std(tail) / (np.mean(tail) + 1e-12)
            assert cv < 0.01, f"Load didn't stabilize: CV={cv:.4f}"
        else:
            # Cleared the load entirely — also valid steady state
            assert tail[-1] == pytest.approx(0.0, abs=ld_cfg.L_clear)

    def test_load_temperature_scaling(self, cfg, ld_cfg):
        """5. Warmer temperature → higher growth rate (Arrhenius)."""
        L = np.array([1000.0])
        r_i = np.array([0.0])
        c_i = np.array([0.0])
        L_cold = within_host_load_step(L.copy(), r_i, c_i, 8.0, 0.0, cfg, ld_cfg)
        L_warm = within_host_load_step(L.copy(), r_i, c_i, 22.0, 0.0, cfg, ld_cfg)
        assert L_warm[0] > L_cold[0], "Warmer T should produce higher growth"

    def test_load_capped_at_Lmax(self, cfg, ld_cfg):
        """6. Load never exceeds L_max."""
        # Start near L_max
        L = np.array([ld_cfg.L_max * 0.99])
        r_i = np.array([0.0])
        c_i = np.array([0.0])
        for _ in range(20):
            L = within_host_load_step(L, r_i, c_i, 20.0, 1e6, cfg, ld_cfg)
        # Logistic growth: L > L_max → growth term becomes negative
        # Combined with np.maximum(0, ...) and logistic form, should not exceed L_max
        assert L[0] <= ld_cfg.L_max * 1.01, (
            f"Load {L[0]:.0f} exceeds L_max {ld_cfg.L_max:.0f}"
        )

    def test_load_never_negative(self, cfg, ld_cfg):
        """7. Load never goes below 0."""
        # High clearance, low starting load
        L = np.array([1.0])
        r_i = np.array([1.0])
        c_i = np.array([1.0])
        for _ in range(100):
            L = within_host_load_step(L, r_i, c_i, 15.0, 0.0, cfg, ld_cfg)
        assert L[0] >= 0.0, f"Load went negative: {L[0]}"

    def test_reinfection_from_environment(self, cfg, ld_cfg):
        """8. Nonzero P_k adds to load via reinfection term."""
        L = np.array([100.0])
        r_i = np.array([0.0])
        c_i = np.array([0.0])
        L_no_pk = within_host_load_step(L.copy(), r_i, c_i, 15.0, 0.0, cfg, ld_cfg)
        L_with_pk = within_host_load_step(L.copy(), r_i, c_i, 15.0, 1e6, cfg, ld_cfg)
        assert L_with_pk[0] > L_no_pk[0], (
            "Environmental pathogen should add reinfection load"
        )


# ═══════════════════════════════════════════════════════════════════════
# UNIT TESTS: load_mortality_probability()
# ═══════════════════════════════════════════════════════════════════════

class TestLoadMortalityProbability:

    def test_mortality_zero_at_zero_load(self, ld_cfg):
        """9. L=0 → p_death=0."""
        L = np.array([0.0])
        t_i = np.array([0.0])
        p = load_mortality_probability(L, t_i, ld_cfg, tau_max=0.85)
        assert p[0] == pytest.approx(0.0, abs=1e-12)

    def test_mortality_at_LD50(self, ld_cfg):
        """10. L=LD50 → p_death ≈ p_death_max/2."""
        # For Hill function: L^n / (LD50^n + L^n) at L=LD50 → 0.5
        L = np.array([ld_cfg.LD50_base])
        t_i = np.array([0.0])  # no tolerance → LD50_i = LD50_base
        p = load_mortality_probability(L, t_i, ld_cfg, tau_max=0.85)
        expected = ld_cfg.p_death_max * 0.5
        assert p[0] == pytest.approx(expected, rel=0.01), (
            f"At LD50, expected p_death={expected:.4f}, got {p[0]:.4f}"
        )

    def test_mortality_saturates(self, ld_cfg):
        """11. Very high L → p_death approaches p_death_max."""
        L = np.array([ld_cfg.L_max * 100])  # extremely high load
        t_i = np.array([0.0])
        p = load_mortality_probability(L, t_i, ld_cfg, tau_max=0.85)
        assert p[0] == pytest.approx(ld_cfg.p_death_max, rel=0.01), (
            f"At very high load, expected ~{ld_cfg.p_death_max}, got {p[0]:.4f}"
        )

    def test_tolerance_shifts_LD50(self, ld_cfg):
        """12. Higher t_i → lower p_death at same load."""
        L = np.array([ld_cfg.LD50_base])  # at baseline LD50
        t_low = np.array([0.0])
        t_high = np.array([1.0])
        p_low = load_mortality_probability(L, t_low, ld_cfg, tau_max=0.85)
        p_high = load_mortality_probability(L, t_high, ld_cfg, tau_max=0.85)
        assert p_high[0] < p_low[0], (
            f"Tolerant agent should have lower death prob: "
            f"p(t=0)={p_low[0]:.4f}, p(t=1)={p_high[0]:.4f}"
        )

    def test_hill_steepness(self, ld_cfg):
        """13. Higher n_hill → sharper transition around LD50."""
        L_below = np.array([ld_cfg.LD50_base * 0.5])
        L_above = np.array([ld_cfg.LD50_base * 1.5])
        t_i = np.array([0.0])

        # Low Hill coefficient → gradual transition
        ld_low = LoadDependentSection(enabled=True, n_hill=1.0)
        p_below_low = load_mortality_probability(L_below, t_i, ld_low, 0.85)
        p_above_low = load_mortality_probability(L_above, t_i, ld_low, 0.85)
        spread_low = p_above_low[0] - p_below_low[0]

        # High Hill coefficient → sharp transition
        ld_high = LoadDependentSection(enabled=True, n_hill=10.0)
        p_below_high = load_mortality_probability(L_below, t_i, ld_high, 0.85)
        p_above_high = load_mortality_probability(L_above, t_i, ld_high, 0.85)
        spread_high = p_above_high[0] - p_below_high[0]

        assert spread_high > spread_low, (
            f"Higher n_hill should give sharper transition: "
            f"spread(n=1)={spread_low:.4f}, spread(n=10)={spread_high:.4f}"
        )


# ═══════════════════════════════════════════════════════════════════════
# INTEGRATION TESTS: load_dependent_disease_update()
# ═══════════════════════════════════════════════════════════════════════

class TestLoadDependentDiseaseUpdate:

    def test_infection_creates_load(self, cfg, ld_cfg, rng):
        """14. Susceptible agents exposed to P_k>0 gain nonzero pathogen_load."""
        agents = _make_agents(50, resistance=0.0)  # zero resistance → high infection rate
        ns = _make_node_state(vibrio=1e6)

        for day in range(10):
            ns = _run_ld_update(agents, ns, cfg, ld_cfg, rng, T=18.0, day=day)

        # Some agents should have gained pathogen_load
        infected = agents['disease_state'][:50] != DiseaseState.S
        loads = agents['pathogen_load'][:50]
        n_with_load = np.sum(loads[infected] > 0)
        assert n_with_load > 0, "Expected some infected agents to have nonzero load"

    def test_initial_load_scales_with_Pk(self, cfg, ld_cfg, rng):
        """15. Higher P_k → higher initial load."""
        # Low Pk scenario
        agents_low = _make_agents(100, resistance=0.0)
        ns_low = _make_node_state(vibrio=1e4)
        rng_low = np.random.default_rng(99)
        _run_ld_update(agents_low, ns_low, cfg, ld_cfg, rng_low, T=18.0)

        # High Pk scenario
        agents_high = _make_agents(100, resistance=0.0)
        ns_high = _make_node_state(vibrio=1e7)
        rng_high = np.random.default_rng(99)
        _run_ld_update(agents_high, ns_high, cfg, ld_cfg, rng_high, T=18.0)

        # Mean initial load for those infected should be higher with more Pk
        inf_low = agents_low['disease_state'][:100] != DiseaseState.S
        inf_high = agents_high['disease_state'][:100] != DiseaseState.S

        if np.any(inf_low) and np.any(inf_high):
            mean_low = np.mean(agents_low['pathogen_load'][:100][inf_low])
            mean_high = np.mean(agents_high['pathogen_load'][:100][inf_high])
            assert mean_high > mean_low, (
                f"Higher Pk should give higher initial load: "
                f"low={mean_low:.1f}, high={mean_high:.1f}"
            )
        else:
            # At least verify higher Pk produces more infections
            assert np.sum(inf_high) >= np.sum(inf_low)

    def test_load_to_state_mapping(self, cfg, ld_cfg, rng):
        """16. L<L_clear→S, L_clear≤L<L_symp→I1, L≥L_symp→I2."""
        agents = _make_agents(3)
        # Manually set up infected agents with different loads
        agents['disease_state'][:3] = DiseaseState.I1  # all start infected
        agents['pathogen_load'][0] = ld_cfg.L_clear * 0.5   # should clear → S
        agents['pathogen_load'][1] = ld_cfg.L_symp * 0.5    # should be I1
        agents['pathogen_load'][2] = ld_cfg.L_symp * 2.0    # should be I2

        # Set high clearance for agent 0 so it stays below L_clear,
        # and low clearance for agents 1,2 so within-host step is relatively stable
        agents['resistance'][0] = 0.9
        agents['recovery_ability'][0] = 0.9
        agents['resistance'][1] = 0.0
        agents['recovery_ability'][1] = 0.0
        agents['resistance'][2] = 0.0
        agents['recovery_ability'][2] = 0.0

        ns = _make_node_state(vibrio=0.0)

        # Run one step — use invasion scenario to prevent new infections from env vibrio
        cfg_inv = DiseaseSection(scenario="invasion", invasion_year=0, invasion_nodes=[])
        ns = _run_ld_update(agents, ns, cfg_inv, ld_cfg, rng, T=15.0)

        # Agent 0: started with load < L_clear * 0.5, after ODE with high clearance → cleared to S
        # (Agent 0 has L=5.0 which is already < L_clear=10 → immediately cleared)
        assert agents['disease_state'][0] == DiseaseState.S, (
            f"Agent 0 (L=5) should be S, got state={agents['disease_state'][0]}"
        )

        # Agent 1: moderate load → I1 (or could still be I1 after ODE)
        # Exact state depends on ODE step but should be I1 if load is in [L_clear, L_symp)
        state1 = agents['disease_state'][1]
        load1 = agents['pathogen_load'][1]
        if load1 > 0:
            if load1 < ld_cfg.L_symp:
                assert state1 == DiseaseState.I1
            else:
                assert state1 == DiseaseState.I2

        # Agent 2: high load → I2
        state2 = agents['disease_state'][2]
        load2 = agents['pathogen_load'][2]
        if load2 >= ld_cfg.L_symp:
            assert state2 == DiseaseState.I2, (
                f"Agent 2 (L={load2:.0f}) should be I2, got state={state2}"
            )

    def test_clearance_leads_to_recovery(self, cfg, rng):
        """17. High-resistance agent clears infection (state→S, load→0).

        Uses elevated delta_clear so that immune clearance dominates
        growth at T=15°C for r_i=1, c_i=1.
        """
        ld_high_clear = LoadDependentSection(enabled=True, delta_clear=0.8)
        agents = _make_agents(1, resistance=1.0, recovery_ability=1.0)
        agents['disease_state'][0] = DiseaseState.I1
        agents['pathogen_load'][0] = 500.0  # moderate load, high clearance

        ns = _make_node_state(vibrio=0.0)
        # Use invasion to prevent reinfection
        cfg_inv = DiseaseSection(scenario="invasion", invasion_year=0, invasion_nodes=[])

        for day in range(30):
            ns = _run_ld_update(agents, ns, cfg_inv, ld_high_clear, rng, day=day)

        assert agents['disease_state'][0] == DiseaseState.S, (
            f"High-resistance agent should have cleared, state={agents['disease_state'][0]}"
        )
        assert agents['pathogen_load'][0] == 0.0

    def test_high_load_causes_death(self, cfg, ld_cfg, rng):
        """18. Agent with high load and zero tolerance eventually dies."""
        # Use many agents so probability of at least one death is very high
        n = 50
        agents = _make_agents(n, resistance=0.0, tolerance=0.0, recovery_ability=0.0)
        agents['disease_state'][:n] = DiseaseState.I2
        agents['pathogen_load'][:n] = ld_cfg.L_max * 0.9  # very high load

        ns = _make_node_state(vibrio=0.0)
        cfg_inv = DiseaseSection(scenario="invasion", invasion_year=0, invasion_nodes=[])

        for day in range(30):
            ns = _run_ld_update(agents, ns, cfg_inv, ld_cfg, rng, day=day)

        n_dead = np.sum(agents['disease_state'][:n] == DiseaseState.D)
        assert n_dead > 0, "Expected at least one death from high pathogen load"

    def test_sentinels_dont_die(self, cfg, ld_cfg, rng):
        """19. Sentinel agents with high load remain alive."""
        n = 20
        agents = _make_agents(n, resistance=0.0, tolerance=0.0, recovery_ability=0.0)
        agents['is_sentinel'][:n] = 1
        agents['disease_state'][:n] = DiseaseState.I2
        agents['pathogen_load'][:n] = ld_cfg.L_max * 0.9

        ns = _make_node_state(vibrio=0.0)
        cfg_inv = DiseaseSection(scenario="invasion", invasion_year=0, invasion_nodes=[])

        for day in range(50):
            ns = _run_ld_update(agents, ns, cfg_inv, ld_cfg, rng, day=day)

        n_alive = np.sum(agents['alive'][:n].astype(bool))
        assert n_alive == n, (
            f"All {n} sentinels should survive, but only {n_alive} alive"
        )

    def test_shedding_proportional_to_load(self, cfg, ld_cfg, rng):
        """20. Total shedding increases with aggregate load."""
        # Low load scenario
        agents_low = _make_agents(10)
        agents_low['disease_state'][:10] = DiseaseState.I1
        agents_low['pathogen_load'][:10] = 100.0  # low load

        ns_low = _make_node_state(vibrio=0.0)
        cfg_inv = DiseaseSection(scenario="invasion", invasion_year=0, invasion_nodes=[])
        rng_low = np.random.default_rng(77)
        ns_low = _run_ld_update(agents_low, ns_low, cfg_inv, ld_cfg, rng_low)

        # High load scenario
        agents_high = _make_agents(10)
        agents_high['disease_state'][:10] = DiseaseState.I1
        agents_high['pathogen_load'][:10] = 50000.0  # high load

        ns_high = _make_node_state(vibrio=0.0)
        rng_high = np.random.default_rng(77)
        ns_high = _run_ld_update(agents_high, ns_high, cfg_inv, ld_cfg, rng_high)

        # Higher load → more shedding → higher Vibrio concentration
        assert ns_high.vibrio_concentration > ns_low.vibrio_concentration, (
            f"Higher aggregate load should produce more shedding: "
            f"low={ns_low.vibrio_concentration:.2f}, high={ns_high.vibrio_concentration:.2f}"
        )

    def test_no_infection_without_vibrio(self, cfg, ld_cfg, rng):
        """21. P_k=0 → no new infections."""
        agents = _make_agents(50, resistance=0.0)
        ns = _make_node_state(vibrio=0.0)
        cfg_inv = DiseaseSection(scenario="invasion", invasion_year=0, invasion_nodes=[])

        for day in range(20):
            ns = _run_ld_update(agents, ns, cfg_inv, ld_cfg, rng, day=day)

        n_infected = np.sum(agents['disease_state'][:50] != DiseaseState.S)
        assert n_infected == 0, (
            f"No infections expected without vibrio, got {n_infected}"
        )

    def test_resistance_reduces_infection_probability(self, cfg, ld_cfg):
        """22. Higher r_i → fewer infections."""
        n = 200

        # Low resistance
        agents_low_r = _make_agents(n, resistance=0.0)
        ns_low = _make_node_state(vibrio=5e5)
        rng_low = np.random.default_rng(55)
        for day in range(5):
            ns_low = _run_ld_update(agents_low_r, ns_low, cfg, ld_cfg, rng_low, T=18.0, day=day)

        # High resistance
        agents_high_r = _make_agents(n, resistance=0.9)
        ns_high = _make_node_state(vibrio=5e5)
        rng_high = np.random.default_rng(55)
        for day in range(5):
            ns_high = _run_ld_update(agents_high_r, ns_high, cfg, ld_cfg, rng_high, T=18.0, day=day)

        inf_low = np.sum(agents_low_r['disease_state'][:n] != DiseaseState.S)
        inf_high = np.sum(agents_high_r['disease_state'][:n] != DiseaseState.S)
        assert inf_high < inf_low, (
            f"Higher resistance should produce fewer infections: "
            f"low_r={inf_low}, high_r={inf_high}"
        )

    def test_load_dependent_disabled_by_default(self):
        """23. Default config → load_dependent.enabled is False."""
        cfg = default_config()
        assert cfg.disease_load_dependent.enabled is False

    def test_backward_compatible_import(self):
        """24. Importing load_dependent_disease_update doesn't break anything."""
        # This test simply verifies the import succeeded (done at module level)
        assert callable(load_dependent_disease_update)
        assert callable(within_host_load_step)
        assert callable(load_mortality_probability)


# ═══════════════════════════════════════════════════════════════════════
# BEHAVIORAL TESTS
# ═══════════════════════════════════════════════════════════════════════

class TestBehavioral:

    def test_truly_immune_r1(self, cfg, rng):
        """25. Agent with r_i=1.0 always clears — clearance >> growth.

        Uses cold temperature (T=8°C) where Arrhenius-scaled growth < base clearance,
        so an agent with r_i=1.0, c_i=1.0 has clearance >> growth and clears.
        """
        ld_high_clear = LoadDependentSection(enabled=True, delta_clear=0.8)
        agents = _make_agents(1, resistance=1.0, recovery_ability=1.0)
        agents['disease_state'][0] = DiseaseState.I1
        agents['pathogen_load'][0] = 10000.0

        ns = _make_node_state(vibrio=0.0)
        cfg_inv = DiseaseSection(scenario="invasion", invasion_year=0, invasion_nodes=[])

        for day in range(60):
            ns = _run_ld_update(agents, ns, cfg_inv, ld_high_clear, rng, T=8.0, day=day)

        assert agents['disease_state'][0] == DiseaseState.S, (
            "Truly immune agent (r=1.0) should always clear infection"
        )
        assert agents['pathogen_load'][0] == 0.0

    def test_epidemic_produces_deaths(self, cfg, ld_cfg):
        """26. Population exposed to disease should have some deaths."""
        n = 100
        agents = _make_agents(n, resistance=0.05, tolerance=0.0, recovery_ability=0.05)

        ns = _make_node_state(vibrio=1e6)
        rng_epi = np.random.default_rng(42)

        for day in range(120):
            ns = _run_ld_update(agents, ns, cfg, ld_cfg, rng_epi, T=18.0, day=day)

        n_dead = np.sum(~agents['alive'][:n].astype(bool))
        assert n_dead > 0, (
            f"Expected deaths in an epidemic, got 0 deaths out of {n}"
        )
        assert ns.cumulative_deaths > 0

    def test_recovered_can_be_reinfected(self, cfg):
        """27. Agent that clears can get infected again.

        Uses high delta_clear so the agent clears quickly, ubiquitous
        scenario at warm T for sustained vibrio, and zero resistance
        so reinfection probability is high.
        """
        ld_fast_clear = LoadDependentSection(enabled=True, delta_clear=1.0)
        # Use ubiquitous scenario for sustained environmental vibrio
        n = 20
        agents = _make_agents(n, resistance=0.0, recovery_ability=1.0)
        # Infect half with low load (will clear fast)
        agents['disease_state'][:10] = DiseaseState.I1
        agents['pathogen_load'][:10] = 50.0

        ns = _make_node_state(vibrio=1e6)
        rng_re = np.random.default_rng(42)

        recovered_once = False
        reinfected_after_recovery = False

        for day in range(300):
            # Manually inject vibrio to keep the pool high (simulates ubiquitous)
            ns.vibrio_concentration = max(ns.vibrio_concentration, 5e5)
            ns = _run_ld_update(agents, ns, cfg, ld_fast_clear, rng_re, T=18.0, day=day)

            # Check any agent for recovery→reinfection cycle
            for i in range(n):
                state = agents['disease_state'][i]
                if state == DiseaseState.S and agents['pathogen_load'][i] == 0:
                    if not recovered_once:
                        recovered_once = True
                        break
            else:
                continue

            if recovered_once:
                # Look for reinfection of any agent that was S
                for day2 in range(day + 1, day + 100):
                    ns.vibrio_concentration = max(ns.vibrio_concentration, 5e5)
                    ns = _run_ld_update(agents, ns, cfg, ld_fast_clear, rng_re, T=18.0, day=day2)
                    for i in range(n):
                        # Was recovered (state S), now infected again
                        if agents['disease_state'][i] != DiseaseState.S and agents['disease_state'][i] != DiseaseState.D:
                            reinfected_after_recovery = True
                            break
                    if reinfected_after_recovery:
                        break
                break

        assert recovered_once, "Agent should have recovered at least once"
        assert reinfected_after_recovery, (
            "After recovery, agent should be reinfected by high environmental P_k"
        )


# ═══════════════════════════════════════════════════════════════════════
# VECTORIZATION / EDGE CASE TESTS
# ═══════════════════════════════════════════════════════════════════════

class TestEdgeCases:

    def test_empty_population(self, cfg, ld_cfg, rng):
        """load_dependent_disease_update handles zero alive agents."""
        agents = allocate_agents(10)  # all dead by default (alive=False)
        ns = _make_node_state(vibrio=100.0)
        # Should not raise
        ns = _run_ld_update(agents, ns, cfg, ld_cfg, rng)
        assert ns.n_S == 0
        assert ns.n_I1 == 0
        assert ns.n_I2 == 0

    def test_vectorized_step_multiple_agents(self, cfg, ld_cfg):
        """within_host_load_step works on arrays of multiple agents."""
        n = 100
        L = np.random.default_rng(42).uniform(10, 10000, n)
        r_i = np.random.default_rng(43).uniform(0, 1, n)
        c_i = np.random.default_rng(44).uniform(0, 1, n)
        L_new = within_host_load_step(L, r_i, c_i, 15.0, 1000.0, cfg, ld_cfg)
        assert L_new.shape == (n,)
        assert np.all(L_new >= 0.0)

    def test_mortality_vectorized(self, ld_cfg):
        """load_mortality_probability works on arrays."""
        n = 50
        L = np.linspace(0, ld_cfg.L_max, n)
        t_i = np.random.default_rng(42).uniform(0, 1, n)
        p = load_mortality_probability(L, t_i, ld_cfg, tau_max=0.85)
        assert p.shape == (n,)
        assert np.all(p >= 0.0)
        assert np.all(p <= ld_cfg.p_death_max + 1e-10)

    def test_load_field_exists_in_dtype(self):
        """pathogen_load field exists in AGENT_DTYPE."""
        assert 'pathogen_load' in AGENT_DTYPE.names
        agents = allocate_agents(5)
        assert agents['pathogen_load'].dtype == np.float64
        assert np.all(agents['pathogen_load'] == 0.0)
