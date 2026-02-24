"""Genetic background definitions for reintroduction experiments.

Generates per-locus allele frequency arrays for different genetic
backgrounds representing captive-bred populations at various stages
of selection and breeding.

Six backgrounds span the plausible range:
  1. pre_sswd      — Pre-2013 wild-type (baseline standing variation)
  2. survivors_2019 — Post-epidemic survivors (~2019, enriched for resistance)
  3. bred_1gen     — 1 generation of selective breeding from survivors
  4. bred_2gen     — 2 generations of selective breeding
  5. bred_5gen     — 5 generations of selective breeding (optimistic program)
  6. optimal       — Theoretical maximum resistance (all R-loci fixed)

Each function returns a dict:
  {
    'name': str,
    'description': str,
    'allele_freqs': np.ndarray (51,),
    'expected_traits': {'resistance': float, 'tolerance': float, 'recovery': float},
    'n_loci': int,
    'trait_slices': {'resistance': slice, 'tolerance': slice, 'recovery': slice},
  }

Usage:
    from scripts.genetic_backgrounds import get_all_backgrounds
    backgrounds = get_all_backgrounds()
    for bg in backgrounds:
        print(f"{bg['name']}: r={bg['expected_traits']['resistance']:.3f}")

Authors: Anton 🔬 & Willem Weertman
"""

from __future__ import annotations

from typing import Dict, List

import numpy as np

from sswd_evoepi.types import N_LOCI, trait_slices
from sswd_evoepi.genetics import initialize_trait_effect_sizes
from sswd_evoepi.model import make_effect_sizes


# ═══════════════════════════════════════════════════════════════════════
# SHARED CONSTANTS
# ═══════════════════════════════════════════════════════════════════════

N_RESISTANCE = 17
N_TOLERANCE = 17
N_RECOVERY = 17

RES_SLICE, TOL_SLICE, REC_SLICE = trait_slices(N_RESISTANCE, N_TOLERANCE, N_RECOVERY)

# Canonical effect sizes (deterministic from seed)
EFFECTS_R = make_effect_sizes(seed=12345, n_loci=N_RESISTANCE)
EFFECTS_T = initialize_trait_effect_sizes(
    np.random.default_rng(12346), N_TOLERANCE, total_weight=1.0,
)
EFFECTS_C = initialize_trait_effect_sizes(
    np.random.default_rng(12347), N_RECOVERY, total_weight=1.0,
)


def _make_result(
    name: str,
    description: str,
    allele_freqs: np.ndarray,
) -> Dict:
    """Build a standard background result dict."""
    # Compute expected trait values from allele frequencies
    # E[trait] = sum(e_l * q_l) for each trait block
    r_expected = float(np.dot(EFFECTS_R, allele_freqs[RES_SLICE]))
    t_expected = float(np.dot(EFFECTS_T, allele_freqs[TOL_SLICE]))
    c_expected = float(np.dot(EFFECTS_C, allele_freqs[REC_SLICE]))

    return {
        'name': name,
        'description': description,
        'allele_freqs': allele_freqs.copy(),
        'expected_traits': {
            'resistance': r_expected,
            'tolerance': t_expected,
            'recovery': c_expected,
        },
        'n_loci': N_LOCI,
        'trait_slices': {
            'resistance': RES_SLICE,
            'tolerance': TOL_SLICE,
            'recovery': REC_SLICE,
        },
    }


def _target_to_uniform_q(effects: np.ndarray, target_mean: float) -> np.ndarray:
    """Compute uniform allele frequencies to achieve a target trait mean.

    For a trait with effect sizes e_l, the expected trait value at
    uniform allele frequency q is: E[trait] = q * sum(e_l).
    So q = target / sum(e_l).

    Returns per-locus q array clipped to [0.001, 0.999].
    """
    total_weight = effects.sum()
    if total_weight <= 0:
        return np.full(len(effects), 0.01)
    q = target_mean / total_weight
    return np.clip(np.full(len(effects), q), 0.001, 0.999)


# ═══════════════════════════════════════════════════════════════════════
# BACKGROUND DEFINITIONS
# ═══════════════════════════════════════════════════════════════════════


def pre_sswd() -> Dict:
    """Pre-2013 wild-type genetic background.

    Represents the standing genetic variation before SSWD hit.
    Uses the model's default initialization targets:
      - resistance: 0.15 (low baseline immune exclusion)
      - tolerance: 0.10 (minimal damage limitation)
      - recovery: 0.02 (very rare pathogen clearance ability)
    """
    allele_freqs = np.zeros(N_LOCI, dtype=np.float64)
    allele_freqs[RES_SLICE] = _target_to_uniform_q(EFFECTS_R, 0.15)
    allele_freqs[TOL_SLICE] = _target_to_uniform_q(EFFECTS_T, 0.10)
    allele_freqs[REC_SLICE] = _target_to_uniform_q(EFFECTS_C, 0.02)

    return _make_result(
        name='pre_sswd',
        description='Pre-2013 wild-type: r=0.15, t=0.10, c=0.02',
        allele_freqs=allele_freqs,
    )


def survivors_2019() -> Dict:
    """Post-epidemic survivors circa 2019.

    Approximates the genetic composition of wild survivors ~6 years
    after the SSWD epidemic began. Values estimated from our
    R→S validation run (5K, 20yr, seed=42):
      - resistance: ~0.18 (enriched by selective mortality)
      - tolerance: ~0.14 (moderate enrichment)
      - recovery: ~0.03 (slight increase, weak selection signal)

    These are hardcoded approximations — running a full 907-node sim
    just for genetic backgrounds would be prohibitively expensive.
    The key insight from validation: resistance evolves fastest under
    R→S biology, recovery barely moves.
    """
    allele_freqs = np.zeros(N_LOCI, dtype=np.float64)
    allele_freqs[RES_SLICE] = _target_to_uniform_q(EFFECTS_R, 0.18)
    allele_freqs[TOL_SLICE] = _target_to_uniform_q(EFFECTS_T, 0.14)
    allele_freqs[REC_SLICE] = _target_to_uniform_q(EFFECTS_C, 0.03)

    return _make_result(
        name='survivors_2019',
        description='Post-epidemic survivors ~2019: r=0.18, t=0.14, c=0.03',
        allele_freqs=allele_freqs,
    )


def bred_generation(
    n_gen: int,
    n_founders: int = 200,
    n_select_frac: float = 0.2,
    seed: int = 42,
) -> Dict:
    """Simulate selective breeding from survivor genetics.

    Starts from survivors_2019 allele frequencies, then runs n_gen
    generations of truncation selection on a weighted index:
      w = 0.7 * resistance + 0.2 * tolerance + 0.1 * recovery

    Each generation:
      1. Generate n_founders individuals from current allele frequencies
      2. Compute selection index for each
      3. Select top fraction (n_select_frac) as parents
      4. Update allele frequencies from selected parents
      5. Apply one round of random mating (allele freq averaging)

    This is a simplified Wright-Fisher-style model that captures
    the direction and approximate magnitude of allele frequency
    change under selection without full pedigree tracking.

    Args:
        n_gen: Number of generations of selective breeding.
        n_founders: Effective population size per generation.
        n_select_frac: Fraction selected as parents (truncation).
        seed: Random seed.

    Returns:
        Background dict with updated allele frequencies.
    """
    rng = np.random.default_rng(seed)

    # Start from survivor genetics
    base = survivors_2019()
    current_freqs = base['allele_freqs'].copy()

    n_select = max(2, int(n_founders * n_select_frac))

    # Selection weights
    w_r, w_t, w_c = 0.7, 0.2, 0.1

    for gen in range(n_gen):
        # Generate founder genotypes from current allele frequencies
        genotypes = np.zeros((n_founders, N_LOCI, 2), dtype=np.int8)
        for l in range(N_LOCI):
            q = current_freqs[l]
            genotypes[:, l, 0] = (rng.random(n_founders) < q).astype(np.int8)
            genotypes[:, l, 1] = (rng.random(n_founders) < q).astype(np.int8)

        # Compute trait scores for all founders
        allele_means = genotypes.sum(axis=2).astype(np.float64) * 0.5
        r_scores = allele_means[:, RES_SLICE] @ EFFECTS_R
        t_scores = allele_means[:, TOL_SLICE] @ EFFECTS_T
        c_scores = allele_means[:, REC_SLICE] @ EFFECTS_C

        # Selection index
        fitness = w_r * r_scores + w_t * t_scores + w_c * c_scores

        # Truncation selection: keep top fraction
        selected_idx = np.argsort(fitness)[-n_select:]

        # Update allele frequencies from selected parents
        selected_geno = genotypes[selected_idx]
        current_freqs = (
            selected_geno.sum(axis=2).sum(axis=0).astype(np.float64)
            / (2.0 * n_select)
        )
        current_freqs = np.clip(current_freqs, 0.001, 0.999)

    gen_label = f'{n_gen}gen'
    return _make_result(
        name=f'bred_{gen_label}',
        description=(
            f'Selective breeding from survivors: {n_gen} generation(s), '
            f'top {n_select_frac:.0%} selected, '
            f'index = 0.7R + 0.2T + 0.1C'
        ),
        allele_freqs=current_freqs,
    )


def optimal() -> Dict:
    """Theoretical maximum resistance.

    All 17 resistance loci fixed at q=1.0 (homozygous derived),
    while tolerance and recovery remain at population mean.
    This represents a theoretical upper bound — unachievable in
    practice but useful as a reference for maximum possible effect.
    """
    allele_freqs = np.zeros(N_LOCI, dtype=np.float64)
    allele_freqs[RES_SLICE] = 1.0  # All resistance loci fixed
    allele_freqs[TOL_SLICE] = _target_to_uniform_q(EFFECTS_T, 0.10)
    allele_freqs[REC_SLICE] = _target_to_uniform_q(EFFECTS_C, 0.02)

    return _make_result(
        name='optimal',
        description='Theoretical max resistance: all 17 R-loci fixed at q=1.0',
        allele_freqs=allele_freqs,
    )


# ═══════════════════════════════════════════════════════════════════════
# CONVENIENCE: ALL BACKGROUNDS
# ═══════════════════════════════════════════════════════════════════════


def get_all_backgrounds() -> List[Dict]:
    """Return all 6 genetic backgrounds in order of increasing resistance.

    Returns:
        List of background dicts:
          [pre_sswd, survivors_2019, bred_1gen, bred_2gen, bred_5gen, optimal]
    """
    return [
        pre_sswd(),
        survivors_2019(),
        bred_generation(n_gen=1),
        bred_generation(n_gen=2),
        bred_generation(n_gen=5),
        optimal(),
    ]


def get_background_by_name(name: str) -> Dict:
    """Look up a genetic background by name.

    Valid names: pre_sswd, survivors_2019, bred_1gen, bred_2gen,
                 bred_5gen, optimal

    Args:
        name: Background name string.

    Returns:
        Background dict.

    Raises:
        ValueError: If name not recognized.
    """
    dispatch = {
        'pre_sswd': pre_sswd,
        'survivors_2019': survivors_2019,
        'bred_1gen': lambda: bred_generation(n_gen=1),
        'bred_2gen': lambda: bred_generation(n_gen=2),
        'bred_5gen': lambda: bred_generation(n_gen=5),
        'optimal': optimal,
    }
    if name not in dispatch:
        raise ValueError(
            f"Unknown background '{name}'. "
            f"Valid: {sorted(dispatch.keys())}"
        )
    return dispatch[name]()


# ═══════════════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    print("Genetic Backgrounds for Reintroduction Experiment")
    print("=" * 60)
    print()

    backgrounds = get_all_backgrounds()
    for bg in backgrounds:
        traits = bg['expected_traits']
        print(f"  {bg['name']:20s}  "
              f"r={traits['resistance']:.4f}  "
              f"t={traits['tolerance']:.4f}  "
              f"c={traits['recovery']:.4f}")
        # Show allele freq summary
        af = bg['allele_freqs']
        print(f"  {'':20s}  "
              f"q_R=[{af[RES_SLICE].min():.3f}-{af[RES_SLICE].max():.3f}]  "
              f"q_T=[{af[TOL_SLICE].min():.3f}-{af[TOL_SLICE].max():.3f}]  "
              f"q_C=[{af[REC_SLICE].min():.3f}-{af[REC_SLICE].max():.3f}]")
        print()

    # Show resistance progression
    r_vals = [bg['expected_traits']['resistance'] for bg in backgrounds]
    print(f"Resistance range: {min(r_vals):.4f} → {max(r_vals):.4f} "
          f"({max(r_vals)/min(r_vals):.1f}× increase)")
