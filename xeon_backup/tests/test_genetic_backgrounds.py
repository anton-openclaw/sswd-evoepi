"""Tests for genetic background definitions (scripts/genetic_backgrounds.py).

Covers:
  - Trait values for each of the 6 canonical backgrounds
  - Monotonic resistance ordering across backgrounds
  - Allele frequency shape, bounds, and slice correctness
  - Lookup helpers (get_all_backgrounds, get_background_by_name)
  - Determinism of bred_generation under fixed seed

Authors: Anton 🔬 & Willem Weertman
"""

import numpy as np
import pytest

from scripts.genetic_backgrounds import (
    EFFECTS_R,
    EFFECTS_T,
    EFFECTS_C,
    N_RESISTANCE,
    N_TOLERANCE,
    N_RECOVERY,
    RES_SLICE,
    TOL_SLICE,
    REC_SLICE,
    bred_generation,
    get_all_backgrounds,
    get_background_by_name,
    optimal,
    pre_sswd,
    survivors_2019,
)
from sswd_evoepi.types import N_LOCI


# ═══════════════════════════════════════════════════════════════════════
# INDIVIDUAL BACKGROUND TRAIT VALUES
# ═══════════════════════════════════════════════════════════════════════


class TestPreSSWD:
    """Pre-2013 wild-type background."""

    def test_pre_sswd_traits(self):
        bg = pre_sswd()
        t = bg['expected_traits']
        assert t['resistance'] == pytest.approx(0.15, abs=0.005)
        assert t['tolerance'] == pytest.approx(0.10, abs=0.005)
        assert t['recovery'] == pytest.approx(0.02, abs=0.005)
        assert bg['name'] == 'pre_sswd'


class TestSurvivors2019:
    """Post-epidemic survivors ~2019."""

    def test_survivors_2019_traits(self):
        bg = survivors_2019()
        t = bg['expected_traits']
        assert t['resistance'] == pytest.approx(0.18, abs=0.005)
        assert t['tolerance'] == pytest.approx(0.14, abs=0.005)
        assert t['recovery'] == pytest.approx(0.03, abs=0.005)
        assert bg['name'] == 'survivors_2019'


class TestBredGenerations:
    """Selective breeding backgrounds (1, 2, 5 generations)."""

    def test_bred_1gen_improves_on_survivors(self):
        surv = survivors_2019()
        b1 = bred_generation(n_gen=1)
        assert b1['expected_traits']['resistance'] > surv['expected_traits']['resistance']
        assert b1['name'] == 'bred_1gen'

    def test_bred_2gen_improves_on_1gen(self):
        b1 = bred_generation(n_gen=1)
        b2 = bred_generation(n_gen=2)
        assert b2['expected_traits']['resistance'] > b1['expected_traits']['resistance']
        assert b2['name'] == 'bred_2gen'

    def test_bred_5gen_improves_on_2gen(self):
        b2 = bred_generation(n_gen=2)
        b5 = bred_generation(n_gen=5)
        assert b5['expected_traits']['resistance'] > b2['expected_traits']['resistance']
        assert b5['name'] == 'bred_5gen'


class TestOptimal:
    """Theoretical maximum resistance."""

    def test_optimal_max_resistance(self):
        bg = optimal()
        t = bg['expected_traits']
        # All 17 resistance loci fixed → E[r] = sum(effects) = 1.0
        assert t['resistance'] == pytest.approx(1.0, abs=0.001)
        assert bg['name'] == 'optimal'

    def test_optimal_resistance_loci_fixed(self):
        bg = optimal()
        r_freqs = bg['allele_freqs'][RES_SLICE]
        np.testing.assert_array_equal(r_freqs, 1.0)


# ═══════════════════════════════════════════════════════════════════════
# STRUCTURAL PROPERTIES
# ═══════════════════════════════════════════════════════════════════════


class TestAlleleFreqProperties:
    """Shape, bounds, and slice correctness for all backgrounds."""

    @pytest.fixture
    def all_backgrounds(self):
        return get_all_backgrounds()

    def test_allele_freqs_shape(self, all_backgrounds):
        for bg in all_backgrounds:
            assert bg['allele_freqs'].shape == (N_LOCI,), (
                f"{bg['name']}: expected shape ({N_LOCI},), "
                f"got {bg['allele_freqs'].shape}"
            )
            assert bg['n_loci'] == N_LOCI

    def test_allele_freqs_bounded(self, all_backgrounds):
        for bg in all_backgrounds:
            af = bg['allele_freqs']
            # All frequencies in [0, 1] (optimal uses exactly 1.0)
            assert np.all(af >= 0.0), f"{bg['name']}: negative allele freq"
            assert np.all(af <= 1.0), f"{bg['name']}: allele freq > 1.0"
            # Non-optimal backgrounds use clipped range [0.001, 0.999]
            if bg['name'] != 'optimal':
                assert np.all(af >= 0.001), f"{bg['name']}: freq below 0.001"
                assert np.all(af <= 0.999), f"{bg['name']}: freq above 0.999"

    def test_trait_slices_correct(self, all_backgrounds):
        bg = all_backgrounds[0]
        slices = bg['trait_slices']
        assert slices['resistance'] == slice(0, N_RESISTANCE)
        assert slices['tolerance'] == slice(N_RESISTANCE, N_RESISTANCE + N_TOLERANCE)
        assert slices['recovery'] == slice(
            N_RESISTANCE + N_TOLERANCE,
            N_RESISTANCE + N_TOLERANCE + N_RECOVERY,
        )
        # Full coverage
        total = N_RESISTANCE + N_TOLERANCE + N_RECOVERY
        assert total == N_LOCI


# ═══════════════════════════════════════════════════════════════════════
# CONVENIENCE FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════


class TestGetAllBackgrounds:
    """Tests for get_all_backgrounds()."""

    def test_get_all_backgrounds_count(self):
        bgs = get_all_backgrounds()
        assert len(bgs) == 6

    def test_get_all_backgrounds_ordering(self):
        """Resistance should increase monotonically across backgrounds."""
        bgs = get_all_backgrounds()
        r_vals = [bg['expected_traits']['resistance'] for bg in bgs]
        for i in range(len(r_vals) - 1):
            assert r_vals[i] < r_vals[i + 1], (
                f"Resistance not monotonically increasing: "
                f"{bgs[i]['name']} ({r_vals[i]:.4f}) >= "
                f"{bgs[i+1]['name']} ({r_vals[i+1]:.4f})"
            )


class TestGetBackgroundByName:
    """Tests for get_background_by_name()."""

    def test_get_background_by_name(self):
        for name in ['pre_sswd', 'survivors_2019', 'bred_1gen',
                      'bred_2gen', 'bred_5gen', 'optimal']:
            bg = get_background_by_name(name)
            assert bg['name'] == name

    def test_get_background_by_name_invalid(self):
        with pytest.raises(ValueError, match="Unknown background"):
            get_background_by_name('nonexistent')


# ═══════════════════════════════════════════════════════════════════════
# DETERMINISM
# ═══════════════════════════════════════════════════════════════════════


class TestDeterminism:
    """Breeding is deterministic under fixed seed."""

    def test_bred_generation_deterministic(self):
        b1a = bred_generation(n_gen=3, seed=99)
        b1b = bred_generation(n_gen=3, seed=99)
        np.testing.assert_array_equal(
            b1a['allele_freqs'], b1b['allele_freqs'],
        )
        assert b1a['expected_traits'] == b1b['expected_traits']
