# Flushing & Residence Time Literature Summary

## Purpose
Inform the choice of `n_connectivity` exponent for mapping enclosedness → flushing rate (φ).

## Key Sources

### Liu et al. 2019 (GRL) — Global Coastal Residence Time
- CRT primarily controlled by **geometric enclosure** (χ = V/S, volume over cross-section area)
- CRT estimator using χ explains **73% of system-to-system variability**
- Relationship is **nonlinear** — CRT increases faster than linearly with χ
- Suggests n > 1 in a direct enclosedness→residence-time mapping

### Khangaonkar et al. 2022 — Salish Sea Empirical Flushing Times
| System | Flushing Time (days) |
|--------|---------------------|
| Open shelf | 1–7 |
| Whidbey Basin / South Sound (more open) | ~30 |
| Puget Sound overall | ~115 |
| Hood Canal (deep fjord-like basin) | ~138 |
| Georgia Basin | ~240 |

### Prince William Sound
- Flushing time **exceeds 1 year** (monthly transports 0.05–0.08 Sv)
- Massive enclosed volume with narrow straits connecting to Gulf of Alaska

## Empirical Scaling Summary

| Category | Residence Time (days) | Approximate φ (relative) |
|----------|----------------------|--------------------------|
| Open coast | 1–7 | 0.8 (high flushing) |
| Semi-enclosed bays | 30–60 | ~0.3–0.5 |
| Puget Sound | 60–115 | ~0.15–0.25 |
| Deep fjords / Hood Canal | 138–240 | ~0.05–0.10 |
| PWS | >365 | ~0.03 (low flushing) |

**Range**: ~50–500× from open coast to deep fjord.

## Implications for n_connectivity

The progression is **nonlinear** but the *direction* of nonlinearity depends on framing:

1. **Residence time vs. enclosedness**: The jump from open coast (1–7d) to semi-enclosed (30–60d) is proportionally huge, but from semi-enclosed to deep fjord (138–365d) is proportionally smaller. This suggests **n < 1** for our φ mapping — the sharpest drop in flushing happens at moderate enclosedness, then plateaus.

2. **Log-linear interpretation**: Residence time scales roughly exponentially with enclosedness, meaning φ (the inverse) drops exponentially. Equivalent to power law with **n ≈ 0.5–0.7**.

3. **Our sweep**: n = 0.5, 1.0, 2.0, 3.0 covers the empirically relevant range. **n = 0.5 may be most physically justified** — sharper flushing drop for even moderately enclosed sites.

## Formula Reference
```
effective = enclosedness × (0.5 + 0.5 × fjord_depth_norm)
effective_nl = effective ^ n_connectivity
φ = φ_open × (1 - effective_nl) + φ_fjord × effective_nl
```
Where φ_open = 0.8, φ_fjord = 0.03.

---
*Compiled Mar 7, 2026*
