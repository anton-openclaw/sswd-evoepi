# Numerical Audit Report

**Auditor:** Automated numerical verification  
**Date:** 2026-03-11  
**Data source:** `summary_data.json` (3-seed means from W154 production model)  
**Sections audited:** sec1_introduction, sec2_key_findings, sec3_scenario_design, sec4_results, sec6_discussion

---

## Summary

- **Total claims checked:** 127
- **CORRECT:** 96
- **ERRORS:** 18 (listed below)
- **WARNINGS:** 13 (listed below)

### Critical Check Results
| Check | Status |
|-------|--------|
| T_VBNC geography (North→9°C floor, South stays high) | ✅ PASS |
| Resistance described as hazard modifier (1−r) | ✅ PASS |
| "Immune exclusion" does NOT appear | ✅ PASS |
| RMSE 0.504 described as calibration RMSE | ✅ PASS |

---

## Error Summary

| # | Section | Claim | Report Value | Data Value | Severity |
|---|---------|-------|-------------|------------|----------|
| E1 | 2, 4 | CA-C final resistance | 0.478 | 0.425 | **HIGH** |
| E2 | 4 | OR final resistance | 0.356 | 0.378 | **HIGH** |
| E3 | 4 | AK-PWS final resistance | 0.261 | 0.253 | MEDIUM |
| E4 | 4 | AK-FN final resistance | 0.206 | 0.230 | **HIGH** |
| E5 | 4 | BC-N final resistance | 0.252 | 0.240 | MEDIUM |
| E6 | 2, 4 | CA-S described as "extinct" (r=0.000) | extinct | 52 indiv, r=0.115 | **HIGH** |
| E7 | 2, 4 | CA-C v_local | 0.244 | 0.235 | MEDIUM |
| E8 | 2, 4 | OR T_vbnc | 10.6°C / 10.62°C | 10.70°C | MEDIUM |
| E9 | 2, 4 | AK-PWS v_local | 0.086 | 0.081 | MEDIUM |
| E10 | 4 | BC-N v_local | 0.139 | 0.144 | MEDIUM |
| E11 | 2, 4 | BJ v_local | 0.343 | 0.335 | MEDIUM |
| E12 | 2 | CA-N "loses 0.9 pp" when evo removed | loses | gains 0.9 pp | **HIGH** |
| E13 | 2 | AK-FS called "Frederick Sound region" | Frederick Sound | Fairweather South | MEDIUM |
| E14 | 2 | CA-C resistance "more than tripling" | >3× | 2.83× (0.425/0.15) | MEDIUM |
| E15 | 2 | "52% penetrates" (from r=0.478) | 52% | ~57.5% (1−0.425) | MEDIUM |
| E16 | 2 | "nearly double the next-best region" | ~2× | 1.48× (17.6/11.9) | MEDIUM |
| E17 | 2→6 | Internal inconsistency: OR resistance | 0.356 (§4) vs ~0.38 (§6) | 0.378 | MEDIUM |
| E18 | 4 | AK-FN called "Funter Bay–northern SE AK" | Funter Bay | Fairweather North | LOW |

---

## Warning Summary

| # | Section | Claim | Report | Data | Note |
|---|---------|-------|--------|------|------|
| W1 | 2, 4 | OR v_local | 0.136 | 0.138 | Off by 0.002 |
| W2 | 4 | BC-N T_vbnc | 9.38°C | 9.39°C | Off by 0.01°C |
| W3 | 4 | CA-C T_vbnc | 11.80°C | 11.81°C | Off by 0.01°C |
| W4 | 2 | AK-FN climate gain | 3.1 pp | 3.0 pp (raw) | From rounded table: 9.9−6.8=3.1 |
| W5 | 2, 4 | OR climate loss | 7.2 pp | 7.1 pp (raw) | From rounded table: 17.6−10.4=7.2 |
| W6 | 4 | F04−F01 coastwide | +0.4 pp | +0.45 pp (raw) | From rounded: 5.8−5.4=0.4 |
| W7 | 2 | AK-PWS resistance "0.26" | 0.26 | 0.253 | Rounding: 0.253→0.25, not 0.26 |
| W8 | 6 | "v_local ≈ 0.35" for southern CA | 0.35 | CA-S: 0.296, BJ: 0.335 | Overestimates |
| W9 | 6 | "2–7% for Alaskan regions" | 2–7% | F04: AK-FN 9.9%, AK-FS 10.3% | Exceeds 7% under F04 |
| W10 | 2 | "three to four times the coastwide avg" | 3–4× | OR: 3.3×, CA-C: 2.2× | Only OR reaches 3× |
| W11 | 6 | CA-C resistance "0.43–0.45" | 0.43–0.45 | 0.425 (≈0.43) | Lower end only |
| W12 | 4 | Multiple diff computations from rounded table | — | — | Systematic; see note |
| W13 | 4 | CA-S table shows 0.0% (F01, F04) | 0.0% | 0.013%, 0.02% | Rounds down to 0.0 at 1 dp |

**Note on W12:** Many percentage-point differences (e.g., +4.2, −3.7, +2.0, +2.3, −1.1, −4.2) are computed from *rounded* table values rather than raw data. All are within ±0.1 pp of the raw difference. This is a consistent methodological choice, not random error, but readers should be aware that re-computing from the raw data yields slightly different deltas.

---

## Detailed Findings

### Section 1: Introduction

- [CORRECT] "90.6% range-wide decline" — external reference (Gravem et al. 2021), not model output
- [CORRECT] "RMSE of 0.504 against normalized regional abundance indices" — correctly described as calibration RMSE. Value not present in `summary_data.json` (which stores `mean_rmse_log = 1.184`, a different metric). Cannot independently verify the 0.504 figure from this data source, but description is appropriate.
- [CORRECT] "896 discrete coastal sites" — consistent with model configuration
- [CORRECT] "peak historical population across all scenarios was 4,447,241 individuals" — data: `total_peak_pop = 4,447,240.67` → rounds to 4,447,241 ✓
- [CORRECT] "three independent random seeds" — consistent with data structure
- [CORRECT] "2012 through 2025" calibration period — consistent
- [CORRECT] No problematic terminology (no "immune exclusion," no "avoidance probability")

### Section 2: Key Findings

#### Coastwide Population Numbers
- [CORRECT] "238,687 surviving individuals" — data F01: 238,686.67 → rounds to 238,687
- [CORRECT] "5.4% of the peak population of 4,447,241" — data: 5.367% → 5.4%
- [CORRECT] "261,632 individuals (5.9%)" — data F02: 261,632.0 and 5.883% → 5.9%
- [CORRECT] "224,638 individuals (5.1%)" — data F03: 224,638.33 → 224,638 and 5.051% → 5.1%
- [CORRECT] "258,760 individuals (5.8%)" — data F04: 258,760.0 and 5.818% → 5.8%
- [CORRECT] "5.1% to 5.9%" range — confirmed: F03=5.1%, F01=5.4%, F04=5.8%, F02=5.9%
- [CORRECT] "approximately 37,000 individuals" difference — F02−F03 = 261,632−224,638 = 36,994 ≈ 37,000
- [CORRECT] "nearly 4.5 million" — 4,447,241 ✓

#### Regional Recovery Rankings (§2.2)
- [CORRECT] "Oregon retains the highest fraction… at 17.6%" — data F01 OR: 17.58% → 17.6%, confirmed #1 across all 18 regions
- [CORRECT] "Central California follows at 11.9%" — data F01 CA-C: 11.95% → 11.9%, confirmed #2
- [CORRECT] "outer Washington coast at 8.8%" — data F01 WA-O: 8.76% → 8.8%, confirmed #3
- [CORRECT] "AK-FN at 6.8%" — data F01 AK-FN: 6.83% → 6.8%, confirmed #4
- [CORRECT] "southern Salish Sea (SS-S) at 6.5%" — data F01 SS-S: 6.53% → 6.5%, confirmed #5
- [CORRECT] "Baja California (BJ)… retaining 0.0%" — data F01 BJ: 0.0%
- [CORRECT] "southern California (CA-S)… retaining 0.0%" — data F01 CA-S: 0.013% → rounds to 0.0% at 1 dp
- [CORRECT] "Aleutians (AK-AL) at 0.4%" — data: 0.37% → 0.4%
- [CORRECT] "AK-OC at 1.8%" — data: 1.79% → 1.8%
- [CORRECT] "AK-EG at 1.7%" — data: 1.72% → 1.7%
- [CORRECT] "AK-PWS… just 2.4%" — data: 2.36% → 2.4%
- [CORRECT] "more than three times the coastwide average" — 17.6/5.4 = 3.26× ✓
- **[ERROR E16]** "nearly double the next-best region (CA-C at 11.9%)" — 17.6/11.9 = 1.48×. This is roughly 50% higher, not "nearly double." Should read "roughly 50% above."

#### Climate Paradox (§2.3)
- [CORRECT] "AK-FS gains 7.6 percentage points" — F04−F01: 10.27−2.70 = 7.57 → 7.6 pp
- [CORRECT] "BC-N gains 7.7 points" — F04−F01: 11.28−3.56 = 7.72 →