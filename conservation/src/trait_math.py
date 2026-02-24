"""Analytical trait distribution calculations.

Implements equations from Section 2-3 of the conservation report:
- Additive genetic variance (Eq. 2.12)
- Trait mean and variance from allele frequencies (Eq. 4.2-4.3)
- Normal approximation for trait distributions
- Selection response (breeder's equation)
- Multi-generation allele frequency prediction
"""

import numpy as np
from scipy import stats
from typing import Optional


def trait_mean(effects: np.ndarray, freqs: np.ndarray) -> float:
    """Expected trait value given effect sizes and allele frequencies.
    
    E[τ] = Σ α_ℓ × q_ℓ  (Eq. 4.2)
    
    Args:
        effects: (n_loci,) effect sizes α_ℓ, summing to 1.0.
        freqs: (n_loci,) protective allele frequencies q_ℓ.
    
    Returns:
        Expected trait value.
    """
    return float(np.dot(effects, freqs))


def trait_variance(effects: np.ndarray, freqs: np.ndarray) -> float:
    """Trait variance from allele frequencies.
    
    σ²_τ = Σ (α_ℓ²/2) × q_ℓ(1-q_ℓ)  (Eq. 4.3)
    
    The factor of 1/2 arises because x_ℓ = (a₁+a₂)/2 with Var[a]=q(1-q).
    
    Args:
        effects: (n_loci,) effect sizes.
        freqs: (n_loci,) allele frequencies.
    
    Returns:
        Trait variance (phenotypic = additive under our model).
    """
    return float(np.sum(effects**2 / 2 * freqs * (1 - freqs)))


def additive_variance(effects: np.ndarray, freqs: np.ndarray) -> float:
    """Additive genetic variance V_A.
    
    V_A = Σ 2 q_ℓ(1-q_ℓ) α_ℓ²  (Eq. 2.12)
    
    Note: V_A = 4 × trait_variance because the trait uses allele means
    (x = (a₁+a₂)/2) while V_A uses allelic values directly.
    In our model with h²=1: V_P = trait_variance, and the 
    selection-relevant variance is trait_variance.
    """
    return float(np.sum(2 * freqs * (1 - freqs) * effects**2))


def selection_response(effects: np.ndarray, freqs: np.ndarray,
                       fraction_selected: float) -> float:
    """Expected selection response per generation (breeder's equation).
    
    R = h² × i × σ_P  (Eq. 3.5)
    
    With h² = 1 in our model: R = i × σ_τ.
    
    Args:
        effects: (n_loci,) effect sizes.
        freqs: (n_loci,) current allele frequencies.
        fraction_selected: Proportion of population selected (0 < p ≤ 1).
    
    Returns:
        Expected change in trait mean per generation.
    """
    sigma = np.sqrt(trait_variance(effects, freqs))
    if sigma == 0:
        return 0.0
    # Selection intensity for truncation selection (normal approximation)
    z = stats.norm.ppf(1 - fraction_selected)
    i = stats.norm.pdf(z) / fraction_selected
    return float(i * sigma)


def delta_q_per_locus(effects: np.ndarray, freqs: np.ndarray,
                      fraction_selected: float) -> np.ndarray:
    """Per-locus allele frequency change under truncation selection.
    
    Δq_ℓ ≈ i × α_ℓ × q_ℓ(1-q_ℓ) / σ_P  (Eq. 3.8)
    
    Args:
        effects: (n_loci,) effect sizes.
        freqs: (n_loci,) current allele frequencies.
        fraction_selected: Proportion selected.
    
    Returns:
        (n_loci,) array of frequency changes.
    """
    sigma = np.sqrt(trait_variance(effects, freqs))
    if sigma == 0:
        return np.zeros_like(freqs)
    z = stats.norm.ppf(1 - fraction_selected)
    i = stats.norm.pdf(z) / fraction_selected
    return i * effects * freqs * (1 - freqs) / sigma


def predict_generations(effects: np.ndarray, freqs_init: np.ndarray,
                        fraction_selected: float,
                        n_generations: int) -> dict:
    """Predict trait evolution over multiple generations.
    
    Iterates the breeder's equation with updating V_A.
    
    Args:
        effects: (n_loci,) effect sizes.
        freqs_init: (n_loci,) initial allele frequencies.
        fraction_selected: Proportion selected each generation.
        n_generations: Number of generations to simulate.
    
    Returns:
        Dictionary with per-generation: means, variances, allele freqs.
    """
    freqs = freqs_init.copy()
    history = {
        'generation': [],
        'mean': [],
        'variance': [],
        'va': [],
        'freqs': [],
        'response': [],
        'n_fixed': [],
    }
    
    for g in range(n_generations + 1):
        mu = trait_mean(effects, freqs)
        var = trait_variance(effects, freqs)
        va = additive_variance(effects, freqs)
        n_fixed = int(np.sum(freqs >= 0.999))
        
        history['generation'].append(g)
        history['mean'].append(mu)
        history['variance'].append(var)
        history['va'].append(va)
        history['freqs'].append(freqs.copy())
        history['n_fixed'].append(n_fixed)
        
        if g < n_generations:
            dq = delta_q_per_locus(effects, freqs, fraction_selected)
            history['response'].append(float(np.dot(effects, dq)))
            freqs = np.clip(freqs + dq, 0.0, 1.0)
        else:
            history['response'].append(0.0)
    
    return history


def generations_to_target(effects: np.ndarray, freqs_init: np.ndarray,
                          fraction_selected: float,
                          target_mean: float,
                          max_generations: int = 100) -> int:
    """How many generations to reach a target trait mean?
    
    Returns:
        Number of generations, or max_generations if not reached.
    """
    freqs = freqs_init.copy()
    for g in range(max_generations):
        mu = trait_mean(effects, freqs)
        if mu >= target_mean:
            return g
        dq = delta_q_per_locus(effects, freqs, fraction_selected)
        freqs = np.clip(freqs + dq, 0.0, 1.0)
    return max_generations


def exceedance_probability(threshold: float, effects: np.ndarray,
                           freqs: np.ndarray) -> float:
    """P(τ ≥ threshold) using normal approximation.
    
    Args:
        threshold: Trait value threshold.
        effects: Effect sizes.
        freqs: Allele frequencies.
    
    Returns:
        Probability of exceeding threshold.
    """
    mu = trait_mean(effects, freqs)
    sigma = np.sqrt(trait_variance(effects, freqs))
    if sigma == 0:
        return 1.0 if mu >= threshold else 0.0
    return float(1 - stats.norm.cdf(threshold, loc=mu, scale=sigma))


def expected_maximum(n_sample: int, effects: np.ndarray,
                     freqs: np.ndarray) -> float:
    """Expected best trait value from a sample of n individuals.
    
    E[τ_(n)] ≈ μ + σ × Φ⁻¹(n/(n+1))  (Eq. 4.6)
    
    Args:
        n_sample: Sample size.
        effects: Effect sizes.
        freqs: Allele frequencies.
    
    Returns:
        Expected maximum trait value.
    """
    mu = trait_mean(effects, freqs)
    sigma = np.sqrt(trait_variance(effects, freqs))
    if sigma == 0 or n_sample <= 1:
        return mu
    return mu + sigma * stats.norm.ppf(n_sample / (n_sample + 1))
