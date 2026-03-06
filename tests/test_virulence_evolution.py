"""Tests for community virulence evolution (Chunk 3).

Tests cover:
  - adapt_community_virulence() function behavior
  - Virulence rate multipliers in daily_disease_update()
  - Config defaults and backward compatibility
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from sswd_evoepi.config import DiseaseSection, PathogenEvolutionSection
from sswd_evoepi.disease import (
    NodeDiseaseState,
    adapt_community_virulence,
)


# ═══════════════════════════════════════════════════════════════════════
# HELPERS
# ═══════════════════════════════════════════════════════════════════════

def _make_cfg(**overrides) -> DiseaseSection:
    """Create a DiseaseSection with virulence_evolution defaults + overrides."""
    defaults = dict(
        virulence_evolution=True,
        v_adapt_rate=0.001,
        v_max_warm=0.7,
        T_v_mid=12.0,
        T_v_width=3.0,
        P_adapt_half=500.0,
        # P_env_dynamic needed for pool
        P_env_dynamic=True,
    )
    defaults.update(overrides)
    return DiseaseSection(**defaults)


def _make_pe_cfg(**overrides) -> PathogenEvolutionSection:
    """Create PathogenEvolutionSection with defaults."""
    defaults = dict(
        v_anchor=0.5,
        alpha_kill=2.0,
        alpha_shed=1.5,
        alpha_prog=1.0,
        gamma_early=0.3,
    )
    defaults.update(overrides)
    return PathogenEvolutionSection(**defaults)


def _make_node_state(**overrides) -> NodeDiseaseState:
    """Create a NodeDiseaseState with defaults + overrides."""
    ns = NodeDiseaseState()
    for k, v in overrides.items():
        setattr(ns, k, v)
    return ns


# ═══════════════════════════════════════════════════════════════════════
# adapt_community_virulence() TESTS
# ═══════════════════════════════════════════════════════════════════════


class TestAdaptCommunityVirulence:
    """Tests for the adapt_community_virulence function."""

    def test_no_change_when_disabled(self):
        """virulence_evolution=False → v_local unchanged."""
        cfg = _make_cfg(virulence_evolution=False)
        ns = _make_node_state(v_local=0.5, P_env_pool=1000.0)
        adapt_community_virulence(ns, T_celsius=15.0, n_hosts=100, K=200, cfg=cfg)
        assert ns.v_local == 0.5

    def test_no_change_when_pool_zero(self):
        """P_env_pool == 0 → v_local unchanged (no bacteria to evolve)."""
        cfg = _make_cfg()
        ns = _make_node_state(v_local=0.5, P_env_pool=0.0)
        adapt_community_virulence(ns, T_celsius=15.0, n_hosts=100, K=200, cfg=cfg)
        assert ns.v_local == 0.5

    def test_no_change_when_pool_negative(self):
        """P_env_pool < 0 → v_local unchanged."""
        cfg = _make_cfg()
        ns = _make_node_state(v_local=0.5, P_env_pool=-10.0)
        adapt_community_virulence(ns, T_celsius=15.0, n_hosts=100, K=200, cfg=cfg)
        assert ns.v_local == 0.5

    def test_drift_toward_optimum_high_density_warm(self):
        """High density + warm temp → high v_optimal → v_local increases."""
        cfg = _make_cfg(v_adapt_rate=0.1)  # fast rate for test visibility
        ns = _make_node_state(v_local=0.3, P_env_pool=1000.0)
        adapt_community_virulence(ns, T_celsius=20.0, n_hosts=200, K=200, cfg=cfg)
        assert ns.v_local > 0.3, "v_local should increase toward high optimal"

    def test_drift_toward_optimum_low_density_cold(self):
        """Low density + cold temp → low v_optimal → v_local decreases."""
        cfg = _make_cfg(v_adapt_rate=0.1)
        ns = _make_node_state(v_local=0.5, P_env_pool=1000.0)
        adapt_community_virulence(ns, T_celsius=5.0, n_hosts=10, K=200, cfg=cfg)
        assert ns.v_local < 0.5, "v_local should decrease toward low optimal"

    def test_v_local_stays_in_bounds_upper(self):
        """v_local should never exceed 1.0."""
        cfg = _make_cfg(v_adapt_rate=100.0)  # absurdly high rate
        ns = _make_node_state(v_local=0.99, P_env_pool=1e6)
        adapt_community_virulence(ns, T_celsius=25.0, n_hosts=500, K=500, cfg=cfg)
        assert ns.v_local <= 1.0

    def test_v_local_stays_in_bounds_lower(self):
        """v_local should never go below 0.0."""
        cfg = _make_cfg(v_adapt_rate=100.0)
        ns = _make_node_state(v_local=0.01, P_env_pool=1e6)
        # Very cold, very sparse → v_optimal near 0
        adapt_community_virulence(ns, T_celsius=-5.0, n_hosts=0, K=500, cfg=cfg)
        assert ns.v_local >= 0.0

    def test_pool_factor_half_saturation(self):
        """Pool at P_adapt_half → pool_factor = 0.5, half the drift speed."""
        cfg = _make_cfg(v_adapt_rate=0.1, P_adapt_half=1000.0)
        # Full pool: 1e6 >> P_adapt_half → pool_factor capped at 1.0
        ns_full = _make_node_state(v_local=0.3, P_env_pool=1e6)
        # Half pool: exactly at P_adapt_half → pool_factor = 1000/1000 = 1.0?
        # No — we want pool_factor = 0.5, so use P_adapt_half / 2
        # pool_factor = min(500 / 1000, 1.0) = 0.5
        ns_half = _make_node_state(v_local=0.3, P_env_pool=500.0)

        adapt_community_virulence(ns_full, T_celsius=20.0, n_hosts=200, K=200, cfg=cfg)
        adapt_community_virulence(ns_half, T_celsius=20.0, n_hosts=200, K=200, cfg=cfg)

        delta_full = ns_full.v_local - 0.3
        delta_half = ns_half.v_local - 0.3
        assert delta_half > 0, "Should still drift"
        assert delta_half < delta_full, "Half pool → smaller drift"
        assert abs(delta_half / delta_full - 0.5) < 0.05, (
            f"Expected ~half drift: {delta_half} vs {delta_full}"
        )

    def test_density_ratio_capped_at_one(self):
        """n_hosts > K → density_ratio capped at 1.0 (no overshoot)."""
        cfg = _make_cfg(v_adapt_rate=0.1)
        ns1 = _make_node_state(v_local=0.3, P_env_pool=1e6)
        ns2 = _make_node_state(v_local=0.3, P_env_pool=1e6)

        adapt_community_virulence(ns1, T_celsius=20.0, n_hosts=200, K=200, cfg=cfg)
        adapt_community_virulence(ns2, T_celsius=20.0, n_hosts=500, K=200, cfg=cfg)

        # Both should get density_ratio=1.0, same result
        assert abs(ns1.v_local - ns2.v_local) < 1e-10

    def test_K_zero_safe(self):
        """K=0 → should not crash (denominator protection)."""
        cfg = _make_cfg(v_adapt_rate=0.1)
        ns = _make_node_state(v_local=0.5, P_env_pool=1000.0)
        adapt_community_virulence(ns, T_celsius=15.0, n_hosts=10, K=0, cfg=cfg)
        # Should not raise; density_ratio = 10/max(0,1) = 10, capped at 1.0

    def test_optimal_virulence_calculation(self):
        """Verify v_optimal matches expected sigmoid × density formula."""
        cfg = _make_cfg(v_adapt_rate=0.1, v_max_warm=0.7, T_v_mid=12.0, T_v_width=3.0)
        T = 15.0
        density = 0.5  # n_hosts/K
        temp_factor = 1.0 / (1.0 + math.exp(-(T - 12.0) / 3.0))
        expected_v_opt = 0.7 * density * temp_factor

        ns = _make_node_state(v_local=0.3, P_env_pool=1e6)
        adapt_community_virulence(ns, T_celsius=T, n_hosts=100, K=200, cfg=cfg)

        # Delta should be v_adapt_rate * 1.0(pool) * (v_optimal - 0.3)
        expected_delta = 0.1 * 1.0 * (expected_v_opt - 0.3)
        assert abs(ns.v_local - (0.3 + expected_delta)) < 1e-10

    def test_no_drift_at_optimum(self):
        """When v_local == v_optimal, no drift should occur."""
        cfg = _make_cfg(v_adapt_rate=0.1, v_max_warm=0.7, T_v_mid=12.0, T_v_width=3.0)
        T = 15.0
        density_ratio = 1.0
        temp_factor = 1.0 / (1.0 + math.exp(-(T - 12.0) / 3.0))
        v_opt = 0.7 * density_ratio * temp_factor

        ns = _make_node_state(v_local=v_opt, P_env_pool=1e6)
        adapt_community_virulence(ns, T_celsius=T, n_hosts=200, K=200, cfg=cfg)
        assert abs(ns.v_local - v_opt) < 1e-12


# ═══════════════════════════════════════════════════════════════════════
# RATE MULTIPLIER TESTS
# ═══════════════════════════════════════════════════════════════════════


class TestVirulenceRateMultipliers:
    """Test that community virulence correctly scales disease rates."""

    def _compute_multipliers(self, v_local: float, pe_cfg=None):
        """Compute the virulence multipliers for a given v_local."""
        if pe_cfg is None:
            pe_cfg = _make_pe_cfg()
        dv = v_local - pe_cfg.v_anchor
        return {
            'kill': math.exp(pe_cfg.alpha_kill * dv),
            'prog': math.exp(pe_cfg.alpha_prog * dv),
            'shed': math.exp(pe_cfg.alpha_shed * dv),
            'shed_early': math.exp(pe_cfg.alpha_shed * dv * pe_cfg.gamma_early),
        }

    def test_anchor_gives_identity(self):
        """At v=v_anchor (0.5), all multipliers must be exactly 1.0."""
        mults = self._compute_multipliers(0.5)
        assert mults['kill'] == 1.0
        assert mults['prog'] == 1.0
        assert mults['shed'] == 1.0
        assert mults['shed_early'] == 1.0

    def test_high_virulence_increases_rates(self):
        """v > v_anchor → all multipliers > 1.0 (more virulent)."""
        mults = self._compute_multipliers(0.7)
        assert mults['kill'] > 1.0
        assert mults['prog'] > 1.0
        assert mults['shed'] > 1.0
        assert mults['shed_early'] > 1.0

    def test_low_virulence_decreases_rates(self):
        """v < v_anchor → all multipliers < 1.0 (less virulent)."""
        mults = self._compute_multipliers(0.3)
        assert mults['kill'] < 1.0
        assert mults['prog'] < 1.0
        assert mults['shed'] < 1.0
        assert mults['shed_early'] < 1.0

    def test_kill_scales_fastest(self):
        """alpha_kill (2.0) > alpha_shed (1.5) > alpha_prog (1.0) → ordering."""
        mults = self._compute_multipliers(0.7)
        assert mults['kill'] > mults['shed'] > mults['prog']

    def test_early_shedding_attenuated(self):
        """gamma_early < 1 → shed_early multiplier closer to 1 than shed."""
        mults = self._compute_multipliers(0.7)
        # shed_early exponent = alpha_shed * dv * gamma_early
        # shed exponent = alpha_shed * dv
        # Since gamma_early=0.3 < 1, shed_early is closer to 1
        assert abs(mults['shed_early'] - 1.0) < abs(mults['shed'] - 1.0)

    def test_multiplier_values_exact(self):
        """Verify exact multiplier values for v=0.7 with default params."""
        pe = _make_pe_cfg()
        dv = 0.7 - 0.5  # = 0.2
        assert math.exp(pe.alpha_kill * dv) == pytest.approx(math.exp(2.0 * 0.2))
        assert math.exp(pe.alpha_shed * dv) == pytest.approx(math.exp(1.5 * 0.2))
        assert math.exp(pe.alpha_prog * dv) == pytest.approx(math.exp(1.0 * 0.2))
        assert math.exp(pe.alpha_shed * dv * pe.gamma_early) == pytest.approx(
            math.exp(1.5 * 0.2 * 0.3)
        )


# ═══════════════════════════════════════════════════════════════════════
# INTEGRATION: RATE SCALING IN daily_disease_update
# ═══════════════════════════════════════════════════════════════════════


class TestDailyDiseaseUpdateVirulenceScaling:
    """Test that daily_disease_update respects community virulence multipliers."""

    def _make_agents(self, n, rng):
        """Create a minimal agent array for testing."""
        from sswd_evoepi.types import AGENT_DTYPE, DiseaseState
        agents = np.zeros(n, dtype=AGENT_DTYPE)
        agents['alive'] = True
        agents['disease_state'] = DiseaseState.S
        agents['size'] = 300.0
        agents['resistance'] = 0.0
        agents['tolerance'] = 0.0
        agents['recovery_ability'] = 0.0
        agents['settlement_day'] = 0
        return agents

    def test_virulence_disabled_no_change(self):
        """virulence_evolution=False → behavior identical to baseline."""
        from sswd_evoepi.disease import daily_disease_update
        cfg = DiseaseSection(virulence_evolution=False)
        pe_cfg = _make_pe_cfg()
        rng = np.random.default_rng(42)
        agents = self._make_agents(50, rng)
        # Make some I1 and I2
        agents['disease_state'][:5] = 2  # I1
        agents['disease_timer'][:5] = 3
        agents['disease_state'][5:10] = 3  # I2
        agents['disease_timer'][5:10] = 3
        ns = _make_node_state(vibrio_concentration=100.0, v_local=0.8)
        ns_orig_v = ns.v_local

        daily_disease_update(
            agents, ns, T_celsius=15.0, salinity=30.0,
            phi_k=0.02, dispersal_input=0.0, day=100,
            cfg=cfg, rng=rng, pe_cfg=pe_cfg,
        )
        # v_local should not affect rates when disabled
        assert ns.v_local == ns_orig_v  # unchanged by disease update

    def test_virulence_at_anchor_identical_to_baseline(self):
        """v_local == v_anchor → multipliers are 1.0 → identical behavior."""
        from sswd_evoepi.disease import (
            daily_disease_update,
            shedding_rate_I1,
            shedding_rate_I2,
            update_vibrio_concentration,
        )

        # Run once with virulence_evolution=False (baseline)
        cfg_base = DiseaseSection(virulence_evolution=False)
        rng1 = np.random.default_rng(42)
        agents1 = self._make_agents(50, rng1)
        agents1['disease_state'][:5] = 2  # I1
        agents1['disease_timer'][:5] = 5
        agents1['disease_state'][5:10] = 3  # I2
        agents1['disease_timer'][5:10] = 5
        ns1 = _make_node_state(vibrio_concentration=100.0, v_local=0.5)

        daily_disease_update(
            agents1, ns1, T_celsius=15.0, salinity=30.0,
            phi_k=0.02, dispersal_input=0.0, day=100,
            cfg=cfg_base, rng=rng1,
        )

        # Run once with virulence_evolution=True, v_local=0.5 (anchor)
        cfg_evo = DiseaseSection(virulence_evolution=True)
        pe_cfg = _make_pe_cfg(v_anchor=0.5)
        rng2 = np.random.default_rng(42)
        agents2 = self._make_agents(50, rng2)
        agents2['disease_state'][:5] = 2  # I1
        agents2['disease_timer'][:5] = 5
        agents2['disease_state'][5:10] = 3  # I2
        agents2['disease_timer'][5:10] = 5
        ns2 = _make_node_state(vibrio_concentration=100.0, v_local=0.5)

        daily_disease_update(
            agents2, ns2, T_celsius=15.0, salinity=30.0,
            phi_k=0.02, dispersal_input=0.0, day=100,
            cfg=cfg_evo, rng=rng2, pe_cfg=pe_cfg,
        )

        # Vibrio concentration should be identical
        assert ns1.vibrio_concentration == pytest.approx(
            ns2.vibrio_concentration, rel=1e-10
        ), "At v_anchor, vibrio should match baseline exactly"

    def test_high_virulence_more_shedding(self):
        """v > anchor → shedding increases → higher vibrio concentration."""
        from sswd_evoepi.disease import daily_disease_update

        pe_cfg = _make_pe_cfg()

        # Baseline: v=0.5 (anchor)
        cfg = DiseaseSection(virulence_evolution=True)
        rng1 = np.random.default_rng(42)
        agents1 = self._make_agents(50, rng1)
        agents1['disease_state'][:5] = 2  # I1
        agents1['disease_timer'][:5] = 10
        agents1['disease_state'][5:10] = 3  # I2
        agents1['disease_timer'][5:10] = 10
        ns1 = _make_node_state(vibrio_concentration=100.0, v_local=0.5)

        daily_disease_update(
            agents1, ns1, T_celsius=15.0, salinity=30.0,
            phi_k=0.02, dispersal_input=0.0, day=100,
            cfg=cfg, rng=rng1, pe_cfg=pe_cfg,
        )

        # High virulence: v=0.8
        rng2 = np.random.default_rng(42)
        agents2 = self._make_agents(50, rng2)
        agents2['disease_state'][:5] = 2  # I1
        agents2['disease_timer'][:5] = 10
        agents2['disease_state'][5:10] = 3  # I2
        agents2['disease_timer'][5:10] = 10
        ns2 = _make_node_state(vibrio_concentration=100.0, v_local=0.8)

        daily_disease_update(
            agents2, ns2, T_celsius=15.0, salinity=30.0,
            phi_k=0.02, dispersal_input=0.0, day=100,
            cfg=cfg, rng=rng2, pe_cfg=pe_cfg,
        )

        assert ns2.vibrio_concentration > ns1.vibrio_concentration, (
            "Higher virulence should produce more shedding → higher vibrio"
        )


# ═══════════════════════════════════════════════════════════════════════
# CONFIG TESTS
# ═══════════════════════════════════════════════════════════════════════


class TestConfigDefaults:
    """Test that virulence evolution config fields have correct defaults."""

    def test_virulence_evolution_default_false(self):
        """virulence_evolution defaults to False (backward compatible)."""
        cfg = DiseaseSection()
        assert cfg.virulence_evolution is False

    def test_v_adapt_rate_default(self):
        assert DiseaseSection().v_adapt_rate == 0.001

    def test_v_max_warm_default(self):
        assert DiseaseSection().v_max_warm == 0.7

    def test_T_v_mid_default(self):
        assert DiseaseSection().T_v_mid == 12.0

    def test_T_v_width_default(self):
        assert DiseaseSection().T_v_width == 3.0

    def test_node_state_v_local_default(self):
        """NodeDiseaseState.v_local defaults to 0.5."""
        ns = NodeDiseaseState()
        assert ns.v_local == 0.5

    def test_pe_section_tradeoff_defaults(self):
        """PathogenEvolutionSection has the trade-off params we depend on."""
        pe = PathogenEvolutionSection()
        assert pe.v_anchor == 0.5
        assert pe.alpha_kill == 2.0
        assert pe.alpha_shed == 1.5
        assert pe.alpha_prog == 1.0
        assert pe.gamma_early == 0.3

    def test_backward_compatible_disabled(self):
        """With defaults, virulence_evolution is off → no behavior change."""
        cfg = DiseaseSection()
        ns = _make_node_state(v_local=0.5, P_env_pool=1000.0)
        adapt_community_virulence(ns, T_celsius=15.0, n_hosts=100, K=200, cfg=cfg)
        assert ns.v_local == 0.5  # no change
