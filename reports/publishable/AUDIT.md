# Numerical Audit Report

**Auditor:** Automated numerical verification  
**Date:** 2026-03-11  
**Data source:** `summary_data.json` (3-seed means from W154 production model)  
**Sections audited:** sec1–sec4, sec6

---

## Summary

- **Total claims checked:** 127
- **CORRECT:** 96
- **ERRORS:** 18
- **WARNINGS:** 13

### Critical Check Results

| Check | Status |
|-------|--------|
| T_VBNC geography (North→9°C floor, South stays high) | ✅ PASS |
| Resistance described as hazard modifier (1−r) | ✅ PASS |
| "Immune exclusion" does NOT appear | ✅ PASS |
| RMSE 0.504 described as calibration RMSE | ✅ PASS |

---

## Errors

### E1 — CA-C final resistance [HIGH] (§2, §4)

**Report:** "mean host resistance rises from 0.15 to 0.48" (§2); "CA-C: 0.150 → 0.478" (§4)  
**Data:** F01 CA-C `yearly_resistance_mean[-1]` = **0.425**  
**Impact:** Cascades to E14 ("more than tripling" → 0.425/0.15 = 2.83×, not >3×) and E15 ("52% penetrates" → should be ~57.5%).  
**Note:** Resistance did reach ~0.49 in some earlier years (indices 33–35) but the final 2050 value is 0.425.

### E2 — OR final resistance [HIGH] (§4)

**Report:** "OR: 0.150 → 0.356"  
**Data:** F01 OR `yearly_resistance_mean[-1]` = **0.378**  
**Note:** The Discussion (§6) correctly states "r ≈ 0.38" — only §4's explicit value is wrong, creating an **internal inconsistency** (see E17).

### E3 — AK-PWS final resistance [MEDIUM] (§4)

**Report:** "AK-PWS: 0.150 → 0.261"  
**Data:** F01 AK-PWS `yearly_resistance_mean[-1]` = **0.253**

### E4 — AK-FN final resistance [HIGH] (§4)

**Report:** "AK-FN: 0.150 → 0.206"  
**Data:** F01 AK-FN `yearly_resistance_mean[-1]` = **0.230**

### E5 — BC-N final resistance [MEDIUM] (§4)

**Report:** "BC-N: 0.150 → 0.252"  
**Data:** F01 BC-N `yearly_resistance_mean[-1]` = **0.240**

### E6 — CA-S described as extinct [HIGH] (§2, §4)

**Report:** "CA-S: 0.150 → 0.000 (population extinct)" (§4); "the population goes extinct before meaningful selection can accumulate, with resistance dropping to 0.00" (§2)  
**Data:** F01 CA-S `final_pop_mean` = **51.67** (~52 individuals surviving); `yearly_resistance_mean[-1]` = **0.115**  
**Impact:** The population is functionally near extinction but NOT literally extinct. The resistance is 0.115, not 0.000.

### E7 — CA-C v_local [MEDIUM] (§2, §4)

**Report:** "a local virulence of 0.244" (§2); "CA-C: 0.244" (§4 table)  
**Data:** F01 CA-C `final_v_local_mean` = **0.235**

### E8 — OR T_vbnc [MEDIUM] (§2, §4)

**Report:** "a lower thermal threshold (10.6°C)" (§2); "OR: T_vbnc = 10.62°C" (§4)  
**Data:** F01 OR `final_T_vbnc_mean` = **10.70°C**  
**Note:** The Discussion (§6) correctly states "T_vbnc ≈ 10.7°C" — only §2 and §4 are wrong.

### E9 — AK-PWS v_local [MEDIUM] (§2, §4)

**Report:** "virulence (0.086)" (§2); "AK-PWS: 0.086" (§4 table)  
**Data:** F01 AK-PWS `final_v_local_mean` = **0.081**

### E10 — BC-N v_local [MEDIUM] (§4)

**Report:** "BC-N: 0.139" (§4 table)  
**Data:** F01 BC-N `final_v_local_mean` = **0.144**

### E11 — BJ v_local [MEDIUM] (§2, §4)

**Report:** "virulence of 0.343" (§2); "BJ: 0.343" (§4 table)  
**Data:** F01 BJ `final_v_local_mean` = **0.335**

### E12 — CA-N direction of evolution effect [HIGH] (§2)

**Report:** "But in Oregon, removing pathogen evolution reduces recovery by 4.8 percentage points… Northern California similarly loses 0.9 points."  
**Data:** F01 CA-N = 4.97%, F02 CA-N = 5.92%. Removing evolution **increases** CA-N recovery by +0.95 pp.  
**Impact:** The claim is backwards. CA-N *gains* under F02, it does not "similarly lose."

### E13 — AK-FS region name [MEDIUM] (§2)

**Report:** "the Frederick Sound region (AK-FS)"  
**Data label:** "Alaska – Fairweather South"  
**Note:** Frederick Sound is a different geographic feature in SE Alaska. The data region is Fairweather South.

### E14 — "more than tripling" [MEDIUM, cascading from E1] (§2)

**Report:** "mean host resistance rises from 0.15 to 0.48 by 2050—more than tripling"  
**Actual:** Using correct data value 0.425: 0.425/0.15 = 2.83×. Not more than tripling.

### E15 — "52% penetrates" [MEDIUM, cascading from E1] (§4)

**Report:** "52% of the pathogen's virulence force still penetrates host defenses"  
**Actual:** Using correct resistance 0.425: 1−0.425 = 0.575 = **57.5%**, not 52%.  
**Note:** §6 correctly states "57%" using r=0.43 — consistent with the data.

### E16 — "nearly double the next-best" [MEDIUM] (§2)

**Report:** "far above any other region… nearly double the next-best region (CA-C at 11.9%)"  
**Actual:** 17.6/11.9 = 1.48×, i.e., ~48% higher. "Nearly double" implies ~2×.

### E17 — Internal inconsistency in OR resistance [MEDIUM] (§4 vs §6)

**Report §4:** "OR: 0.150 → 0.356"  
**Report §6:** "Oregon populations have evolved resistance to r ≈ 0.38 by 2050"  
**Data:** 0.378. The Discussion is correct; §4 is wrong.

### E18 — AK-FN region name [LOW] (§4)

**Report:** "AK-FN (Funter Bay–northern Southeast Alaska)"  
**Data label:** "Alaska – Fairweather North"  
**Note:** Funter Bay is a different location. Minor geographic mislabeling.

---

## Warnings

### W1 — OR v_local: 0.136 vs 0.138 (§2, §4)

Report says 0.136; data says 0.138. Off by 0.002. Acceptable rounding.

### W2 — BC-N T_vbnc: 9.38°C vs 9.39°C (§4)

Off by 0.01°C. Likely a different rounding convention.

### W3 — CA-C T_vbnc: 11.80°C vs 11.81°C (§4)

Off by 0.01°C. Acceptable.

### W4 — AK-FN climate gain: 3.1 pp vs 3.0 pp raw (§2, §4)

Raw: 9.85−6.83 = 3.02 pp. Rounded table: 9.9−6.8 = 3.1 pp. Report uses rounded-table value.

### W5 — OR climate loss: 7.2 pp vs 7.1 pp raw (§2, §4)

Raw: 17.58−10.44 = 7.14 pp. Rounded table: 17.6−10.4 = 7.2 pp. Report uses rounded-table value.

### W6 — F04−F01 coastwide: +0.4 pp vs +0.45 pp raw (§4)

Raw: 5.818−5.367 = 0.451 pp. Rounded: 5.8−5.4 = 0.4 pp.

### W7 — AK-PWS resistance "0.26" (§2)

Report (§2) says "resistance reaches only 0.26"; data final value = 0.253. Rounds to 0.25 at 2 dp, not 0.26. Borderline between warning and error.

### W8 — "v_local ≈ 0.35" for southern CA (§6)

Report says "high evolved virulence (v_local ≈ 0.35)." Data: CA-S v_local = 0.296, BJ v_local = 0.335. Neither is 0.35. The tilde makes this approximate, but it overstates both.

### W9 — "2–7% recovery for Alaskan regions" (§6)

Under F04, AK-FN = 9.9% and AK-FS = 10.3%, exceeding the stated 7% upper bound. The claim holds for F01 but not "both climate scenarios" as stated.

### W10 — "three to four times the coastwide average" (§2)

Report: "Oregon and central California show recovery rates three to four times the coastwide average."  
Data: OR = 17.6/5.4 = 3.3×. CA-C = 11.9/5.4 = 2.2×. Only Oregon reaches 3×; CA-C is ~2.2×.

### W11 — CA-C resistance "0.43–0.45" (§6)

Data: final value = 0.425. The lower bound ~0.43 is roughly right, but 0.45 overstates. In some prior years resistance reached 0.49, so the range may reflect temporal variation rather than the endpoint.

### W12 — Systematic use of rounded-table differences

Multiple inter-scenario differences (e.g., +4.2, −3.7, +2.0, +2.3, −1.1, −4.2 pp) are computed from rounded 1-dp table values rather than raw data. All are within ±0.15 pp of raw values. Consistent methodological choice but may mislead readers who try to reproduce from raw data.

### W13 — CA-S table shows 0.0% under F01 and F04

Data: F01 CA-S = 0.013%, F04 CA-S = 0.02%. Both round to 0.0% at 1 dp. Technically correct rounding but obscures that ~52 individuals persist.

---

## Detailed Findings by Section

### Section 1: Introduction

| # | Claim | Verdict | Notes |
|---|-------|---------|-------|
| 1 | "90.6% range-wide decline" | CORRECT | External ref (Gravem et al. 2021) |
| 2 | "RMSE of 0.504" | UNVERIFIABLE | Not in summary_data.json; correctly described as calibration RMSE |
| 3 | "896 discrete coastal sites" | CONSISTENT | Matches model config |
| 4 | "peak … 4,447,241 individuals" | CORRECT | Data: 4,447,240.67 → rounds to 4,447,241 |
| 5 | "three independent random seeds" | CORRECT | Consistent with data |
| 6 | No "immune exclusion" terminology | CORRECT | Not present |

### Section 2: Key Findings

| # | Claim | Verdict | Notes |
|---|-------|---------|-------|
| 7 | F01: 238,687 survivors (5.4%) | CORRECT | 238,686.67→238,687; 5.367→5.4% |
| 8 | F02: 261,632 (5.9%) | CORRECT | 261,632.0; 5.883→5.9% |
| 9 | F03: 224,638 (5.1%) | CORRECT | 224,638.33→224,638; 5.051→5.1% |
| 10 | F04: 258,760 (5.8%) | CORRECT | 258,760.0; 5.818→5.8% |
| 11 | Range 5.1%–5.9% | CORRECT | Min F03=5.05%, Max F02=5.88% |
| 12 | ~37,000 difference best-worst | CORRECT | 261,632−224,638 = 36,994 |
| 13 | OR 17.6% = highest region | CORRECT | Data: 17.58%, rank #1 confirmed |
| 14 | CA-C 11.9%, WA-O 8.8%, AK-FN 6.8%, SS-S 6.5% | CORRECT | Rankings verified |
| 15 | BJ 0.0%, CA-S 0.0% | CORRECT | Rounding of 0.0% and 0.013% |
| 16 | AK-AL 0.4%, AK-OC 1.8%, AK-EG 1.7%, AK-PWS 2.4% | CORRECT | All match within rounding |
| 17 | "more than three times the coastwide average" | CORRECT | 17.6/5.4 = 3.26× |
| 18 | "nearly double the next-best region" | **ERROR E16** | 17.6/11.9 = 1.48× |
| 19 | AK-FS gains 7.6 pp under F04 | CORRECT | 10.27−2.70 = 7.57→7.6 |
| 20 | BC-N gains 7.7 pp | CORRECT | 11.28−3.56 = 7.72→7.7 |
| 21 | AK-OC gains 4.2 pp | CORRECT | 5.96−1.79 = 4.17→4.2 |
| 22 | AK-FN gains 3.1 pp | **WARNING W4** | Raw: 3.02→3.0; table: 3.1 |
| 23 | OR loses 7.2 pp | **WARNING W5** | Raw: 7.14→7.1; table: 7.2 |
| 24 | CA-C loses 4.3 pp | CORRECT | 11.95−7.63 = 4.32→4.3 |
| 25 | CA-N loses 0.7 pp | CORRECT | 4.97−4.31 = 0.66→0.7 |
| 26 | "Frederick Sound region (AK-FS)" | **ERROR E13** | Data: "Fairweather South" |
| 27 | CA-C resistance 0.15→0.48 | **ERROR E1** | Data: 0.425 |
| 28 | "more than tripling" | **ERROR E14** | 0.425/0.15 = 2.83× |
| 29 | OR resistance 0.36 | **ERROR E2** (approx) | Data: 0.378 |
| 30 | AK-PWS resistance "0.26" | **WARNING W7** | Data: 0.253 |
| 31 | CA-S r→0.00, "extinct" | **ERROR E6** | Data: r=0.115, 52 indiv |
| 32 | CA-C T_vbnc 11.8°C | CORRECT | Data: 11.81°C |
| 33 | CA-C v_local 0.244 | **ERROR E7** | Data: 0.235 |
| 34 | OR T_vbnc 10.6°C | **ERROR E8** | Data: 10.70°C |
| 35 | OR v_local 0.136 | **WARNING W1** | Data: 0.138 |
| 36 | PWS T_vbnc 9.0°C floor | CORRECT | Data: 9.001°C |
| 37 | PWS v_local 0.086 | **ERROR E9** | Data: 0.081 |
| 38 | BJ T_vbnc 12.0°C | CORRECT | Data: 11.996°C |
| 39 | BJ v_local 0.343 | **ERROR E11** | Data: 0.335 |
| 40 | "52% penetrates" | **ERROR E15** | 1−0.425 = 57.5% |
| 41 | F01−F03 = 0.3 pp coastwide | CORRECT | 5.37−5.05 = 0.32→0.3 |
| 42 | AK-AL gains 3.7 pp (F02−F01) | CORRECT | 4.07−0.37 = 3.70 |
| 43 | AK-WG gains 2.0 pp | CORRECT | 4.51−2.45 = 2.06; table 2.0 |
| 44 | BC-C gains 4.2 pp | CORRECT | 7.47−3.34 = 4.13; table 4.2 |
| 45 | OR loses 4.8 pp (F01→F02) | CORRECT | 17.58−12.77 = 4.81→4.8 |
| 46 | CA-N "similarly loses 0.9 pp" | **ERROR E12** | Actually gains 0.95 pp |
| 47 | F02 coastwide 5.9% | CORRECT | 5.883→5.9% |
| 48 | "three to four times" for OR+CA-C | **WARNING W10** | Only OR≥3× |
| 49 | "virulence increases fourfold" AK→BJ | CORRECT | ~4× even with slightly off values |

### Section 3: Scenario Design

| # | Claim | Verdict | Notes |
|---|-------|---------|-------|
| 50 | RMSE = 0.504 | UNVERIFIABLE | Calibration metric, not in summary data |
| 51 | K = 5,000 per cell | CORRECT | Data: K=5000 |
| 52 | r₀ = 0.150 | CORRECT | All regions start at ~0.150 |
| 53 | 896 cells, 2012–2050 | CONSISTENT | |
| 54 | Three seeds | CORRECT | |
| 55 | Initial T_vbnc = 12.0°C | CORRECT | F02 data (no evolution): all T_vbnc = 12.0 |
| 56 | +0.8°C (SSP2-4.5), +1.5°C (SSP5-8.5) | UNVERIFIABLE | Not in summary data |

### Section 4: Results

| # | Claim | Verdict | Notes |
|---|-------|---------|-------|
| 57 | Coast-wide table (F01–F04 survivors/%) | CORRECT | All 4 rows verified |
| 58 | Regional table (18 regions × 4 scenarios = 72 cells) | CORRECT | All values match within 1-dp rounding |
| 59 | Bold values = max per region | CORRECT | Spot-checked 8 regions |
| 60 | "AK-FN (Funter Bay–northern SE AK)" | **ERROR E18** | Label: "Fairweather North" |
| 61 | AK-FS +7.6 pp relative to F01 | CORRECT | |
| 62 | AK-OC +4.2 pp (F04−F01) | CORRECT | |
| 63 | AK-PWS: 0.150→0.261 | **ERROR E3** | Data: 0.253 |
| 64 | AK-FN: 0.150→0.206 | **ERROR E4** | Data: 0.230 |
| 65 | BC-N: 0.150→0.252 | **ERROR E5** | Data: 0.240 |
| 66 | OR: 0.150→0.356 | **ERROR E2** | Data: 0.378 |
| 67 | CA-C: 0.150→0.478 | **ERROR E1** | Data: 0.425 |
| 68 | CA-S: 0.150→0.000 (extinct) | **ERROR E6** | 52 indiv, r=0.115 |
| 69 | "52% penetrates" at r=0.478 | **ERROR E15** | Should be ~57.5% |
| 70 | AK-PWS T_vbnc = 9.00°C | CORRECT | Data: 9.001 |
| 71 | AK-FN T_vbnc = 9.00°C | CORRECT | Data: 9.002 |
| 72 | BC-N T_vbnc = 9.38°C | **WARNING W2** | Data: 9.39°C |
| 73 | OR T_vbnc = 10.62°C | **ERROR E8** | Data: 10.70°C |
| 74 | CA-C T_vbnc = 11.80°C | **WARNING W3** | Data: 11.81°C |
| 75 | BJ T_vbnc = 12.00°C | CORRECT | Data: 11.996 ≈ 12.00 |
| 76 | AK-PWS v_local = 0.086 | **ERROR E9** | Data: 0.081 |
| 77 | AK-FN v_local = 0.091 | CORRECT | Data: 0.091 |
| 78 | BC-N v_local = 0.139 | **ERROR E10** | Data: 0.144 |
| 79 | OR v_local = 0.136 | **WARNING W1** | Data: 0.138 |
| 80 | CA-C v_local = 0.244 | **ERROR E7** | Data: 0.235 |
| 81 | BJ v_local = 0.343 | **ERROR E11** | Data: 0.335 |
| 82 | F04−F01 coastwide +0.4 pp | **WARNING W6** | Raw: +0.45 pp |
| 83 | Northern beneficiaries list (5 items) | CORRECT | All verified |
| 84 | Southern losers list (4 items) | CORRECT | All verified |
| 85 | Breakeven ~48°N | CORRECT | Approximately correct |
| 86 | F02−F01 coastwide +0.5 pp | CORRECT | 5.88−5.37 = 0.51→0.5 |
| 87 | OR loses 4.8 pp under F02 | CORRECT | |
| 88 | WA-O loses 3.7 pp under F02 | CORRECT | Table-rounded: 8.8−5.1=3.7 |
| 89 | BC-C gains 4.2 pp, AK-AL 3.7, AK-WG 2.0, SS-S 2.3 | CORRECT | All from table-rounded |
| 90 | F03 effects: CA-N −2.1, OR −2.3, WA-O −4.6, AK-FN −1.1 | CORRECT | All verified |
| 91 | F01−F03 coastwide −0.3 pp | CORRECT | |
| 92 | Synthesis: +4.8% OR, −4.2% BC-C | CORRECT | OR exact, BC-C from rounded |

### Section 6: Discussion

| # | Claim | Verdict | Notes |
|---|-------|---------|-------|
| 93 | "94–95% below pre-disease levels" | CORRECT | Range 94.1–95.0% |
| 94 | F01: 94.6% decline | CORRECT | 100−5.367 = 94.633→94.6 |
| 95 | F04: 94.2% decline | CORRECT | 100−5.818 = 94.182→94.2 |
| 96 | OR 17.6% under F01 | CORRECT | |
| 97 | OR T_vbnc ≈ 10.7°C | CORRECT | Data: 10.70°C ✓ |
| 98 | OR v_local ≈ 0.14 | CORRECT | Data: 0.138 ≈ 0.14 |
| 99 | OR r ≈ 0.38 by 2050 | CORRECT | Data: 0.378 ≈ 0.38 |
| 100 | OR drops to 10.4% under F04 | CORRECT | Data: 10.44→10.4% |
| 101 | "2–7% recovery for Alaskan regions" | **WARNING W9** | F04: AK-FN 9.9%, AK-FS 10.3% |
| 102 | AK-FN 6.8%→9.9% = 45% improvement | CORRECT | (9.9−6.8)/6.8 = 45.6%→45% |
| 103 | CA-S ~0.01% recovery, ~52 individuals | CORRECT | Data: 0.013%, 51.67 |
| 104 | CA-C ~12% | CORRECT | Data: 11.95% |
| 105 | CA-C→7.6% under F04 | CORRECT | Data: 7.63→7.6% |
| 106 | v_local ≈ 0.35 for southern CA | **WARNING W8** | CA-S: 0.296, BJ: 0.335 |
| 107 | CA-C resistance "0.43–0.45" | **WARNING W11** | Data: 0.425 ≈ 0.43 |
| 108 | "57% penetrates" for r=0.43 | CORRECT | 1−0.43 = 0.57 ✓ |
| 109 | "nearly tripling" (0.15→0.43) | CORRECT | 0.43/0.15 = 2.87× ≈ "nearly 3×" |
| 110 | VA drops 88% (0.026→0.003) | CORRECT | (0.0262−0.00316)/0.0262 = 87.9% |
| 111 | "effective daily infection … well above 20%" | UNVERIFIABLE | Model-internal calculation |

---

## Highest-Priority Fixes

1. **All five resistance endpoint values in §4 evolutionary dynamics** (E1–E5) must be corrected. These propagate into the §2 narrative. The data values are:
   - AK-PWS: 0.253 (not 0.261)
   - AK-FN: 0.230 (not 0.206)
   - BC-N: 0.240 (not 0.252)
   - OR: 0.378 (not 0.356)
   - CA-C: 0.425 (not 0.478)

2. **CA-S is not extinct** (E6). It has ~52 individuals with r = 0.115. Recommend changing to "functionally extinct" or "near-extinction" and revising the resistance claim.

3. **CA-N direction of evolution effect** (E12). Section 2 says CA-N "similarly loses" 0.9 pp when evolution is removed. Data shows CA-N *gains* 0.95 pp (F02 > F01). This reverses the stated conclusion. Remove or correct.

4. **OR T_vbnc** (E8). Sections 2 and 4 say 10.6°C / 10.62°C, but data says 10.70°C. The Discussion (§6) correctly uses ~10.7°C. Harmonize to 10.7°C.

5. **Virulence values** (E7, E9–E11). The v_local table in §4 and narrative in §2 contain 4 incorrect values. Correct CA-C to 0.235, AK-PWS to 0.081, BC-N to 0.144, BJ to 0.335.

6. **Cascading claims** from E1: "more than tripling" (E14) → should say "nearly tripling"; "52% penetrates" (E15) → should be ~57.5%; §2 should match §6's correct values.

7. **Region names** (E13, E18): "Frederick Sound" → "Fairweather South"; "Funter Bay" → "Fairweather North."

---

## Note on Methodology

Many inter-scenario percentage-point differences in the report are computed from rounded table values (1 decimal place) rather than raw data. This is a defensible practice when the table is the primary reference, but it systematically introduces ±0.1 pp discrepancies. All such cases are flagged as warnings, not errors.

The §6 Discussion generally uses more accurate (correctly rounded) values than §2 and §4, suggesting the Discussion may have been written or revised with direct access to the data while §2 and §4 may have been generated from an earlier or intermediate data pipeline.
