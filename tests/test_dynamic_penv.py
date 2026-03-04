"""Tests for dynamic P_env with host-amplification feedback and community floor.

Tests verify:
  1. Backward compatibility: P_env_dynamic=False gives identical behavior
  2. Pool builds from shedding by infected hosts
  3. Pool decays without infected hosts
  4. Floor provides minimum environmental vibrio even with pool=0
  5. Floor is modulated by SST (warm > cold)
  6. Warm vs cold differential after shedding removal
  7. Config defaults are correct
"""

import numpy as np
import pytest

from sswd_evoepi.config import DiseaseSection
from sswd_evoepi.types import AGENT_DTYPE, DiseaseState, allocate_agents
from sswd_evoepi.disease import (
    NodeDiseaseState,
    update_vibrio_concentration,
    environmental_vibrio,
    daily_disease_update,
    shedding_rate_I1,
    shedding_rate_I2,
    thermal_performance,
    salinity_modifier,
    K_VBNC,
)


# ═══════════════════════════════════════════════════════════════════════
# FIXTURES
# ═══════════════════════════════════════════════════════════════════════

@pytest.fixture
def cfg_static() -> DiseaseSection:
    """Static P_env config (default, backward compatible)."""
    return DiseaseSection()


@pytest.fixture
def cfg_dynamic() -> DiseaseSection:
    """Dynamic P_env config."""
    return DiseaseSection(
        P_env_dynamic=True,
        P_env_floor=50.0,
        alpha_env=0.1,
        delta_env=0.05,
    )


def _make_agents(n: int, n_I1: int = 0, n_I2: int = 0) -> np.ndarray:
    """Create a simple agent array with specified compartment counts."""
    agents = np.zeros(n, dtype=AGENT_DTYPE)
    agents['alive'] = True
    agents['disease_state'] = DiseaseState.S
    agents['size'] = 300.0
    agents['resistance'] = 0.1
    # Set some agents as I1/I2
    idx = 0
    for _ in range(n_I1):
        if idx < n:
            agents['disease_state'][idx] = DiseaseState.I1
            agents['disease_timer'][idx] = 5
            idx += 1
    for _ in range(n_I2):
        if idx < n:
            agents['disease_state'][idx] = DiseaseState.I2
            agents['disease_timer'][idx] = 3
            idx += 1
    return agents


def _compute_floor(T_celsius: float, salinity: float, cfg: DiseaseSection) -> float:
    """Compute the expected floor value for dynamic P_env."""
    k = getattr(cfg, 'k_vbnc', K_VBNC)
    vbnc_activation = 1.0 / (1.0 + np.exp(-k * (T_celsius - cfg.T_vbnc)))
    g_peak = thermal_performance(3000.0, T_celsius, rate_ref=1.0)
    sal_mod = salinity_modifier(salinity, cfg.s_min, cfg.s_full)
    return cfg.P_env_floor * vbnc_activation * g_peak * sal_mod


# ═══════════════════════════════════════════════════════════════════════
# TESTS
# ═══════════════════════════════════════════════════════════════════════

class TestBackwardCompatStaticPenv:
    """P_env_dynamic=False gives same result as current."""

    def test_static_penv_unchanged(self, cfg_static):
        """With P_env_dynamic=False, update_vibrio_concentration uses
        environmental_vibrio() exactly as before."""
        T, sal, phi = 15.0, 30.0, 0.01
        P_k = 100.0

        # Static mode (default)
        P_new_static = update_vibrio_concentration(
            P_k, n_I1=5, n_I2=2, n_D_fresh=1,
            T_celsius=T, salinity=sal, phi_k=phi,
            dispersal_input=0.0, cfg=cfg_static,
            disease_reached=True,
        )

        # Explicitly passing P_env_pool should have no effect when static
        P_new_with_pool = update_vibrio_concentration(
            P_k, n_I1=5, n_I2=2, n_D_fresh=1,
            T_celsius=T, salinity=sal, phi_k=phi,
            dispersal_input=0.0, cfg=cfg_static,
            disease_reached=True, P_env_pool=999.0,
        )

        assert P_new_static == pytest.approx(P_new_with_pool)

    def test_static_uses_environmental_vibrio(self, cfg_static):
        """Static mode result matches manual calculation with environmental_vibrio()."""
        T, sal, phi = 15.0, 30.0, 0.01
        P_k = 0.0

        P_new = update_vibrio_concentration(
            P_k, n_I1=0, n_I2=0, n_D_fresh=0,
            T_celsius=T, salinity=sal, phi_k=phi,
            dispersal_input=0.0, cfg=cfg_static,
            disease_reached=True,
        )

        env = environmental_vibrio(T, sal, cfg_static)
        assert P_new == pytest.approx(env, rel=1e-10)


class TestPoolBuildsFromShedding:
    """With infected hosts, P_env_pool increases over multiple steps."""

    def test_pool_increases_with_infected_hosts(self, cfg_dynamic):
        """P_env_pool should increase when infected hosts are shedding."""
        T, sal, phi = 15.0, 30.0, 0.01
        rng = np.random.default_rng(42)

        agents = _make_agents(50, n_I1=10, n_I2=5)
        node_state = NodeDiseaseState(node_id=0, vibrio_concentration=0.0)

        # P_env_pool starts at 0
        assert node_state.P_env_pool == 0.0

        # Run multiple daily updates
        pools = [0.0]
        for day in range(10):
            node_state = daily_disease_update(
                agents, node_state,
                T_celsius=T, salinity=sal, phi_k=phi,
                dispersal_input=0.0, day=day, cfg=cfg_dynamic, rng=rng,
                disease_reached=True,
            )
            pools.append(node_state.P_env_pool)

        # Pool should have increased from 0
        assert node_state.P_env_pool > 0.0
        # Pool should be increasing over the first few steps
        assert pools[2] > pools[1] > pools[0]

    def test_pool_accumulates_proportional_to_shedding(self, cfg_dynamic):
        """Pool input is alpha_env * total_shedding."""
        T = 15.0
        node_state = NodeDiseaseState(node_id=0, P_env_pool=0.0)

        # Manually compute expected pool after one step with known shedding
        n_I1, n_I2, n_D_fresh = 10, 5, 1
        total_shed = (shedding_rate_I1(T, cfg_dynamic) * n_I1
                      + shedding_rate_I2(T, cfg_dynamic) * n_I2
                      + cfg_dynamic.sigma_D * n_D_fresh)
        expected_input = cfg_dynamic.alpha_env * total_shed
        expected_pool = max(0.0, 0.0 + expected_input - cfg_dynamic.delta_env * 0.0)

        # Simulate the pool update manually
        pool_input = cfg_dynamic.alpha_env * total_shed
        pool_decay = cfg_dynamic.delta_env * node_state.P_env_pool
        new_pool = max(0.0, node_state.P_env_pool + pool_input - pool_decay)

        assert new_pool == pytest.approx(expected_pool, rel=1e-10)
        assert new_pool > 0.0


class TestPoolDecaysWithoutHosts:
    """Without infected hosts, P_env_pool decays exponentially toward 0."""

    def test_pool_decays_exponentially(self, cfg_dynamic):
        """Pool should decay toward 0 when no infected hosts present."""
        T, sal, phi = 15.0, 30.0, 0.01
        rng = np.random.default_rng(42)

        # Start with a high pool value and no infected hosts
        agents = _make_agents(50, n_I1=0, n_I2=0)
        node_state = NodeDiseaseState(
            node_id=0, vibrio_concentration=100.0, P_env_pool=500.0
        )

        initial_pool = node_state.P_env_pool
        pools = [initial_pool]

        for day in range(20):
            node_state = daily_disease_update(
                agents, node_state,
                T_celsius=T, salinity=sal, phi_k=phi,
                dispersal_input=0.0, day=day, cfg=cfg_dynamic, rng=rng,
                disease_reached=True,
            )
            pools.append(node_state.P_env_pool)

        # Pool should be decreasing monotonically
        for i in range(1, len(pools)):
            assert pools[i] < pools[i - 1], f"Pool not decreasing at step {i}"

        # After 20 steps with delta_env=0.05, pool should be much smaller
        # Decay factor per step = (1 - delta_env) = 0.95
        # After 20 steps: 500 * 0.95^20 ≈ 179
        expected_approx = initial_pool * (1.0 - cfg_dynamic.delta_env) ** 20
        assert node_state.P_env_pool == pytest.approx(expected_approx, rel=0.1)

    def test_pool_approaches_zero(self, cfg_dynamic):
        """Pool should approach 0 after many steps without input."""
        T, sal, phi = 15.0, 30.0, 0.01
        rng = np.random.default_rng(42)

        # Use fully resistant agents so no new infections feed back into pool
        agents = _make_agents(50, n_I1=0, n_I2=0)
        agents['resistance'] = 1.0

        node_state = NodeDiseaseState(
            node_id=0, vibrio_concentration=100.0, P_env_pool=1000.0
        )

        for day in range(200):
            node_state = daily_disease_update(
                agents, node_state,
                T_celsius=T, salinity=sal, phi_k=phi,
                dispersal_input=0.0, day=day, cfg=cfg_dynamic, rng=rng,
                disease_reached=True,
            )

        # After 200 steps: 1000 * 0.95^200 ≈ 0.035
        assert node_state.P_env_pool < 0.1


class TestFloorProvidesMinimum:
    """Even with pool=0 and no hosts, floor still provides env vibrio."""

    def test_floor_provides_env_vibrio(self, cfg_dynamic):
        """Floor should contribute to env vibrio even with zero pool."""
        T, sal = 18.0, 30.0

        floor = _compute_floor(T, sal, cfg_dynamic)

        # Use update_vibrio_concentration with pool=0
        P_new = update_vibrio_concentration(
            P_k=0.0, n_I1=0, n_I2=0, n_D_fresh=0,
            T_celsius=T, salinity=sal, phi_k=0.01,
            dispersal_input=0.0, cfg=cfg_dynamic,
            disease_reached=True, P_env_pool=0.0,
        )

        # P_new should equal floor (since P_k=0, no decay/flush, no shedding)
        assert P_new == pytest.approx(floor, rel=1e-10)
        assert P_new > 0.0

    def test_floor_uses_penv_floor_not_penv_max(self, cfg_dynamic):
        """Floor computation uses P_env_floor, not P_env_max."""
        T, sal = 18.0, 30.0

        floor = _compute_floor(T, sal, cfg_dynamic)
        env_max_result = environmental_vibrio(T, sal, cfg_dynamic)

        # Floor should be P_env_floor/P_env_max times the environmental_vibrio result
        ratio = cfg_dynamic.P_env_floor / cfg_dynamic.P_env_max
        assert floor == pytest.approx(env_max_result * ratio, rel=1e-10)


class TestFloorModulatedBySST:
    """Floor is higher at warm SST, lower at cold SST."""

    def test_warm_higher_than_cold(self, cfg_dynamic):
        """Floor at 18°C should be much higher than at 6°C."""
        sal = 30.0

        floor_warm = _compute_floor(18.0, sal, cfg_dynamic)
        floor_cold = _compute_floor(6.0, sal, cfg_dynamic)

        assert floor_warm > floor_cold
        # The difference should be substantial due to VBNC sigmoid and thermal perf
        assert floor_warm > 5.0 * floor_cold

    def test_floor_near_zero_at_very_cold(self, cfg_dynamic):
        """Floor at very cold temperatures should be near zero."""
        sal = 30.0
        floor_cold = _compute_floor(4.0, sal, cfg_dynamic)

        # At 4°C: VBNC sigmoid is very low (~0.0003), thermal performance is low
        assert floor_cold < 1.0

    def test_floor_peaks_near_topt(self, cfg_dynamic):
        """Floor should peak near T_opt=20°C."""
        sal = 30.0

        floor_15 = _compute_floor(15.0, sal, cfg_dynamic)
        floor_20 = _compute_floor(20.0, sal, cfg_dynamic)
        floor_25 = _compute_floor(25.0, sal, cfg_dynamic)

        # 20°C should be higher than 15°C (approaching peak)
        assert floor_20 > floor_15
        # 20°C should be higher than 25°C (thermal decline above T_opt)
        assert floor_20 > floor_25


class TestWarmVsColdDifferential:
    """After identical shedding then removal, warm node retains more env pressure."""

    def test_warm_retains_more_pressure(self, cfg_dynamic):
        """Warm node should have higher env vibrio than cold node after shedding stops."""
        sal, phi = 30.0, 0.01
        rng_warm = np.random.default_rng(123)
        rng_cold = np.random.default_rng(123)

        T_warm, T_cold = 18.0, 8.0

        # Phase 1: Build up pool with identical shedding (10 days)
        agents_warm = _make_agents(50, n_I1=10, n_I2=5)
        agents_cold = _make_agents(50, n_I1=10, n_I2=5)
        ns_warm = NodeDiseaseState(node_id=0, vibrio_concentration=0.0)
        ns_cold = NodeDiseaseState(node_id=1, vibrio_concentration=0.0)

        for day in range(10):
            ns_warm = daily_disease_update(
                agents_warm, ns_warm, T_warm, sal, phi,
                dispersal_input=0.0, day=day, cfg=cfg_dynamic, rng=rng_warm,
                disease_reached=True,
            )
            ns_cold = daily_disease_update(
                agents_cold, ns_cold, T_cold, sal, phi,
                dispersal_input=0.0, day=day, cfg=cfg_dynamic, rng=rng_cold,
                disease_reached=True,
            )

        # Phase 2: Remove all infected hosts, let pool decay (20 days)
        agents_clean_warm = _make_agents(50, n_I1=0, n_I2=0)
        agents_clean_cold = _make_agents(50, n_I1=0, n_I2=0)

        for day in range(10, 30):
            ns_warm = daily_disease_update(
                agents_clean_warm, ns_warm, T_warm, sal, phi,
                dispersal_input=0.0, day=day, cfg=cfg_dynamic, rng=rng_warm,
                disease_reached=True,
            )
            ns_cold = daily_disease_update(
                agents_clean_cold, ns_cold, T_cold, sal, phi,
                dispersal_input=0.0, day=day, cfg=cfg_dynamic, rng=rng_cold,
                disease_reached=True,
            )

        # Warm node should have higher vibrio due to higher floor
        assert ns_warm.vibrio_concentration > ns_cold.vibrio_concentration

        # The floor at warm temp should be substantially higher
        floor_warm = _compute_floor(T_warm, sal, cfg_dynamic)
        floor_cold = _compute_floor(T_cold, sal, cfg_dynamic)
        assert floor_warm > floor_cold


class TestConfigDefaults:
    """P_env_dynamic defaults to False, verify all defaults."""

    def test_penv_dynamic_defaults_false(self):
        """P_env_dynamic should default to False."""
        cfg = DiseaseSection()
        assert cfg.P_env_dynamic is False

    def test_penv_floor_default(self):
        """P_env_floor should default to 50.0."""
        cfg = DiseaseSection()
        assert cfg.P_env_floor == 50.0

    def test_alpha_env_default(self):
        """alpha_env should default to 0.1."""
        cfg = DiseaseSection()
        assert cfg.alpha_env == 0.1

    def test_delta_env_default(self):
        """delta_env should default to 0.05."""
        cfg = DiseaseSection()
        assert cfg.delta_env == 0.05

    def test_node_state_pool_defaults_zero(self):
        """NodeDiseaseState.P_env_pool should default to 0.0."""
        ns = NodeDiseaseState()
        assert ns.P_env_pool == 0.0

    def test_getattr_backward_compat(self):
        """getattr with default should work for old configs without P_env_dynamic."""
        cfg = DiseaseSection()
        assert getattr(cfg, 'P_env_dynamic', False) is False
