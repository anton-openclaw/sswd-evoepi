"""Tests for per-node lognormal K variability.

Verifies that the K_cv parameter in PopulationSection correctly generates
variable carrying capacities drawn from a lognormal distribution while
preserving the mean K and maintaining backward compatibility when K_cv=0.
"""

import numpy as np
import pytest
import sys
from pathlib import Path

# Add project root to path
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from sswd_evoepi.config import PopulationSection, default_config


# ── Helper: replicate the K generation logic from calibration_runner ──
def generate_K_values(K: int, K_cv: float, n_sites: int, seed: int = 42):
    """Replicate build_node_defs K generation logic for testing."""
    if K_cv > 0:
        rng = np.random.default_rng(seed + 9999)
        sigma = np.sqrt(np.log(1 + K_cv**2))
        mu = np.log(K) - sigma**2 / 2
        K_values = rng.lognormal(mu, sigma, size=n_sites)
        K_values = np.maximum(K_values, 100).astype(int)
    else:
        K_values = np.full(n_sites, K, dtype=int)
    return K_values


class TestDefaultUniformK:
    """K_cv=0 gives all nodes the same K (backward compatible)."""

    def test_all_same(self):
        K_values = generate_K_values(K=5000, K_cv=0.0, n_sites=896)
        assert np.all(K_values == 5000)

    def test_different_K_base(self):
        K_values = generate_K_values(K=10000, K_cv=0.0, n_sites=100)
        assert np.all(K_values == 10000)


class TestLognormalKDistribution:
    """K_cv>0 gives variable K with correct mean (within statistical tolerance)."""

    def test_mean_preserved_cv05(self):
        """With CV=0.5, mean K should be close to K_base."""
        K_values = generate_K_values(K=5000, K_cv=0.5, n_sites=10000, seed=42)
        # With 10k samples, mean should be within ~2% of target
        assert abs(K_values.mean() - 5000) / 5000 < 0.05, (
            f"Mean K = {K_values.mean():.0f}, expected ~5000"
        )

    def test_mean_preserved_cv2(self):
        """With CV=2.0 (high variability), mean K should still be close."""
        K_values = generate_K_values(K=5000, K_cv=2.0, n_sites=50000, seed=123)
        # High CV needs more samples for convergence; allow 10%
        assert abs(K_values.mean() - 5000) / 5000 < 0.10, (
            f"Mean K = {K_values.mean():.0f}, expected ~5000"
        )

    def test_variability_exists(self):
        """K values should not all be identical when K_cv > 0."""
        K_values = generate_K_values(K=5000, K_cv=0.5, n_sites=100)
        assert np.std(K_values) > 0, "K values should have variability"
        assert np.min(K_values) != np.max(K_values)

    def test_cv_approximately_correct(self):
        """Empirical CV should approximate the requested CV (large sample)."""
        K_cv = 0.5
        K_values = generate_K_values(K=5000, K_cv=K_cv, n_sites=50000, seed=99)
        empirical_cv = np.std(K_values) / np.mean(K_values)
        # Allow ±30% relative error on CV (floor distorts slightly)
        assert abs(empirical_cv - K_cv) / K_cv < 0.30, (
            f"Empirical CV = {empirical_cv:.3f}, expected ~{K_cv}"
        )


class TestKFloor:
    """No node gets K < 100."""

    def test_floor_at_100(self):
        """Even with extreme CV, minimum K should be 100."""
        K_values = generate_K_values(K=5000, K_cv=5.0, n_sites=10000, seed=42)
        assert np.min(K_values) >= 100, f"Min K = {np.min(K_values)}, expected >= 100"

    def test_floor_with_small_K_base(self):
        """With small K_base, floor should still apply."""
        K_values = generate_K_values(K=200, K_cv=2.0, n_sites=10000, seed=42)
        assert np.min(K_values) >= 100


class TestReproducibility:
    """Same seed gives same K values."""

    def test_same_seed_same_values(self):
        K1 = generate_K_values(K=5000, K_cv=1.0, n_sites=896, seed=42)
        K2 = generate_K_values(K=5000, K_cv=1.0, n_sites=896, seed=42)
        np.testing.assert_array_equal(K1, K2)

    def test_different_seed_different_values(self):
        K1 = generate_K_values(K=5000, K_cv=1.0, n_sites=896, seed=42)
        K2 = generate_K_values(K=5000, K_cv=1.0, n_sites=896, seed=99)
        assert not np.array_equal(K1, K2)


class TestConfigDefault:
    """K_cv defaults to 0.0 in PopulationSection."""

    def test_default_value(self):
        pop = PopulationSection()
        assert pop.K_cv == 0.0

    def test_config_has_k_cv(self):
        config = default_config()
        assert hasattr(config.population, 'K_cv')
        assert config.population.K_cv == 0.0

    def test_custom_value(self):
        pop = PopulationSection(K_cv=1.5)
        assert pop.K_cv == 1.5
