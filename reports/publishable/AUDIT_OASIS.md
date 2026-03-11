# Audit: Section 4b — Site-Level Dynamics and the Oasis Question

**File audited:** `sections/sec4b_site_dynamics.tex`  
**Data sources:** `results/calibration/F01/monthly_seed{42,123,999}.npz`  
**Audit date:** 2026-03-11  
**Auditor:** automated (Claude, subagent oasis-audit)

---

## Numerical Claims

### 1. Oregon site counts: "17–21 of 59 sites retain population"

**CORRECT ✅**

| Seed | Alive sites | Total sites |
|------|-------------|-------------|
| 42   | 17          | 59          |
| 123  | 21          | 59          |
| 999  | 17          | 59          |

Range: 17–21. Matches claim exactly.

---

### 2. Oregon top 3 concentration: "39–79%"

**CORRECT ✅** (with rounding)

| Seed | Top 3 fraction |
|------|----------------|
| 42   | 79.0%          |
| 123  | 43.1%          |
| 999  | 38.9%          |

Actual range: 38.9–79.0%. The lower bound (38.9%) rounds to 39% at integer precision. Claim is valid.

---

### 3. Oregon Gini: "0.85–0.92"

**CORRECT ✅** (with rounding)

| Seed | Gini   |
|------|--------|
| 42   | 0.9153 |
| 123  | 0.8629 |
| 999  | 0.8450 |

Actual range: 0.845–0.915. At 2 decimal places: 0.85–0.92 (0.845 → 0.85 by round-half-up; 0.9153 → 0.92). Claim is valid.

---

### 4. AK-PWS site counts: "37–38 of 43 sites"

**CORRECT ✅**

| Seed | Alive sites | Total sites |
|------|-------------|-------------|
| 42   | 38          | 43          |
| 123  | 37          | 43          |
| 999  | 38          | 43          |

---

### 5. AK-PWS top 3 concentration: "25–29%"

**CORRECT ✅**

| Seed | Top 3 fraction |
|------|----------------|
| 42   | 28.6%          |
| 123  | 25.4%          |
| 999  | 26.8%          |

Range: 25.4–28.6%, which falls within the claimed 25–29%.

---

### 6. AK-PWS Gini: "0.65–0.71"

**CORRECT ✅**

| Seed | Gini   |
|------|--------|
| 42   | 0.7069 |
| 123  | 0.6516 |
| 999  | 0.6728 |

Range: 0.65–0.71 at 2 dp. Matches claim.

---

### 7. AK-FN site counts: "37–40 of 40 sites"

**CORRECT ✅**

| Seed | Alive sites | Total sites |
|------|-------------|-------------|
| 42   | 38          | 40          |
| 123  | 37          | 40          |
| 999  | 40          | 40          |

Range: 37–40. Matches claim exactly.

Note: The text states "(93–100%)" for the fraction — 37/40 = 92.5% ≈ 93%, 40/40 = 100%. Correct.

---

### 8. Cross-seed top 5: "13 different sites, only 2 in more than one seed"

**CORRECT ✅**

| Seed | Top 5 sites                              |
|------|------------------------------------------|
| 42   | OR-048, OR-049, OR-055, OR-052, OR-058   |
| 123  | OR-046, OR-047, OR-053, OR-055, OR-049   |
| 999  | OR-020, OR-050, OR-059, OR-044, OR-045   |

- **13 unique sites** across all three top-5 lists: OR-020, OR-044, OR-045, OR-046, OR-047, OR-048, OR-049, OR-050, OR-052, OR-053, OR-055, OR-058, OR-059
- **2 sites appear in >1 seed:** OR-049 (seeds 42, 123) and OR-055 (seeds 42, 123)

Matches claim exactly.

---

### 9. Southern Oregon zone: "42.0–42.8°N"

**CORRECT ✅** (with caveat acknowledged in text)

Top sites by seed and their latitudes:

**Seed 42:** All top 10 sites in 42.04–42.67°N ✅  
**Seed 123:** Top 5 in 42.04–42.73°N; sites at 44.68°N and 44.83°N appear at ranks 6–7 ✅  
**Seed 999:** **OR-020 (44.60°N) is the #1 site** — but this is explicitly noted in the text as the counterexample ("in seed 999, OR-020 (44.60°N)—a site 220 km to the north—is the dominant oasis"). The remaining 4 of top 5 are in 42.04–42.75°N.

The text's claim that "dominant survivors concentrate along the southern Oregon coast between 42.0°N and 42.8°N" is broadly correct — this zone dominates in all seeds except for the OR-020 outlier which is explicitly discussed.

---

### 10. Oasis sites OR-048 (42.67°N) and OR-020 (44.60°N)

**CORRECT ✅**

| Site   | Verified lat | Claimed lat |
|--------|-------------|-------------|
| OR-048 | 42.6744°N   | 42.67°N     |
| OR-020 | 44.5968°N   | 44.60°N     |

Additional claims verified:
- "In seed 42, OR-048 holds 37% of regional population" → 2790 / 7542 = **37.0%** ✅
- "In seed 999, OR-020 is the dominant oasis" → OR-020 has 3035, largest in seed 999 ✅
- "A site 220 km to the north" → Δlat ≈ 1.92° × 111 km/° ≈ 214 km. Claim of "220 km" is a reasonable approximation ✅

---

### 11. Initial conditions: "Gini = 0.00 at start, K = 5,000"

**CORRECT ✅**

- Oregon t=0 population: min=4997, max=5000, mean=4999 (effectively uniform at K ≈ 5000)
- Gini at t=0: 0.0001 (effectively 0.00)

---

### 12. AK-FN Gini: "0.63–0.66" (mentioned in text)

**MINOR ERROR ⚠️**

| Seed | Gini   |
|------|--------|
| 42   | 0.6298 |
| 123  | 0.6453 |
| 999  | 0.6583 |

Actual range at 2 dp: 0.63–0.66. Matches claim. However, seed 42 is 0.6298 which rounds to **0.63** at 2 dp, so the lower bound is borderline. Technically correct.

**CORRECT ✅** on closer examination.

---

## Biological Accuracy Checks

### Sweepstakes reproductive success (Hedgecock 2011)

**CORRECT ✅**

The citation is appropriate. Hedgecock & Pudovkin (2011) formalized the sweepstakes reproductive success hypothesis for marine broadcast spawners: at any given time, only a small fraction of adults contribute to successful recruitment due to matching of reproductive timing with favorable oceanographic conditions. This mechanism correctly explains why, at small population sizes, stochastic variation in which individuals happen to reproduce successfully creates enormous site-level variance.

### Allee effects logic

**CORRECT ✅**

The text argues that Allee effects (reduced per-capita fitness at low density) combine with sweepstakes reproduction to amplify stochastic variance in recovery outcomes. This is logically sound: below an Allee threshold, populations have reduced reproductive output, making successful recruitment events rarer and more stochastic. At the post-crash population sizes described (~400–2,000 across the region), Allee effects are biologically plausible for broadcast-spawning marine invertebrates.

### "Oasis identity is stochastic but oasis geography is deterministic"

**SUPPORTED ✅**

The data strongly supports this claim:
- **Geography is deterministic:** In all three seeds, the majority of the top 5 sites (12 of 13 unique sites, excluding OR-020) fall in the 42.0–42.8°N zone. This geographic consistency across independent stochastic realizations indicates a deterministic geographic advantage.
- **Identity is stochastic:** 13 different sites appear across three top-5 lists, with only 2 sites recurring. Zero sites appear in all three seeds' top 5. The specific "winners" within the favorable zone are not predictable.
- **The OR-020 exception** (seed 999, lat 44.6°N) shows that even geography is not absolute — a site well outside the primary oasis zone can dominate through a lucky reproductive event. The text honestly acknowledges this.

---

## Summary

| # | Claim | Verdict |
|---|-------|---------|
| 1 | OR alive sites: 17–21 of 59 | ✅ CORRECT |
| 2 | OR top 3 concentration: 39–79% | ✅ CORRECT (38.9% rounds to 39%) |
| 3 | OR Gini: 0.85–0.92 | ✅ CORRECT (0.845 rounds to 0.85) |
| 4 | AK-PWS alive: 37–38 of 43 | ✅ CORRECT |
| 5 | AK-PWS top 3: 25–29% | ✅ CORRECT |
| 6 | AK-PWS Gini: 0.65–0.71 | ✅ CORRECT |
| 7 | AK-FN alive: 37–40 of 40 | ✅ CORRECT |
| 8 | Cross-seed: 13 unique, 2 in >1 seed | ✅ CORRECT |
| 9 | Southern OR zone: 42.0–42.8°N | ✅ CORRECT (exception noted in text) |
| 10 | OR-048 @ 42.67°N, OR-020 @ 44.60°N | ✅ CORRECT |
| 11 | Initial Gini = 0.00, K = 5000 | ✅ CORRECT |
| 12 | AK-FN Gini: 0.63–0.66 | ✅ CORRECT |
| Bio | Sweepstakes reproduction | ✅ CORRECT |
| Bio | Allee effects | ✅ CORRECT |
| Bio | Stochastic identity / deterministic geography | ✅ SUPPORTED |

**Overall: All 15 claims verified. Zero errors found. Two claims (OR top 3 lower bound, OR Gini lower bound) rely on rounding from raw values that sit just below the stated boundary — consider adding one decimal place of precision or noting "approximately" in the text for maximum rigor.**
