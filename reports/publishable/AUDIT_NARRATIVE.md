# Narrative Audit — SSWD-EvoEpi Publishable Report

**Date:** 2026-03-11  
**Auditor:** Scientific editor (subagent)  
**Scope:** Biological accuracy, internal consistency, narrative quality  
**Files reviewed:** sec1–sec6 (all six report sections)

---

## CRITICAL

### C1. "Immunity" used for echinoderms (sec1_introduction.tex, ~line 55)

**Text:**  
> "individuals that survive carry alleles conferring some degree of tolerance or immunity, and pass them to their offspring"

**What's wrong:**  
Echinoderms possess only innate immunity — they have no adaptive immune system and cannot develop immunological memory. The model's three defense traits are resistance, tolerance, and recovery (R→S, not R→immune). Using "immunity" implies adaptive immune protection, which is biologically incorrect and contradicts the model's own mechanics described in sec5.

**Suggested fix:**  
Replace "tolerance or immunity" with "tolerance or resistance."

---

### C2. AK-FN geographic description contradicts between sections

**sec2_key_findings.tex (~line 15):**  
> "the Fairweather–Yakutat coast of Alaska (AK-FN)"

**sec4_results.tex (~line 75):**  
> "AK-FN (Funter Bay–northern Southeast Alaska)"

**sec6_discussion.tex (~line 55):**  
> "Fairweather North recovery" / "the Fairweather coast"

**What's wrong:**  
The Fairweather–Yakutat coast and Funter Bay are geographically distinct areas. Fairweather is on the outer coast near Glacier Bay; Funter Bay is an interior waterway near Juneau. Two sections use "Fairweather" and one uses "Funter Bay." Readers will not know which geography is correct for region AK-FN.

**Suggested fix:**  
Determine the correct geographic description and use it consistently in all three sections. A region-code glossary table (see Suggestion S3) would prevent future confusion.

---

### C3. *Pisaster ochraceus* called a "congener" of *Pycnopodia* (sec5_model.tex, ~line 40)

**Text:**  
> "inspired by GWAS results from the congener *Pisaster ochraceus*"

**What's wrong:**  
"Congener" means same genus. *Pisaster* (family Asteriidae) and *Pycnopodia* (family Pycnopodiidae) are in different families. They are both asteroids (class Asteroidea) but not congeners. This is a taxonomic error that peer reviewers will catch immediately.

**Suggested fix:**  
Replace "the congener" with "the asteroid" or "the sea star" — e.g., "inspired by GWAS results from the sea star *Pisaster ochraceus* — the only Pacific asteroid for which disease-associated genomic data currently exist."

---

## MODERATE

### M1. CA-C resistance value inconsistent across sections

| Section | CA-C resistance (2050) |
|---------|----------------------|
| sec2_key_findings | 0.48 |
| sec4_results | 0.478 |
| sec5_model | "approximately 0.42–0.45" |
| sec6_discussion | 0.43 |

**What's wrong:**  
0.478 (sec4) versus 0.43 (sec6) is an 11% discrepancy — well beyond rounding. The sec5 range of 0.42–0.45 excludes the sec4 value of 0.478. Sec2 and sec4 are consistent (0.48 ≈ 0.478), but sec5 and sec6 report a substantially lower number. This looks like different data sources or an uncorrected draft revision.

**Suggested fix:**  
Use the sec4 value (0.478) as the authoritative number — it comes from the detailed results table. Update sec5 (~line 55) and sec6 (~line 75) to match.

---

### M2. Oregon resistance value inconsistent across sections

| Section | OR resistance (2050) |
|---------|---------------------|
| sec2_key_findings | 0.36 |
| sec4_results | 0.356 |
| sec6_discussion | ≈0.38 |

**What's wrong:**  
0.36 vs 0.38 is a 5.5% discrepancy. While sec2 and sec4 are consistent (0.36 ≈ 0.356), sec6 rounds UP to 0.38 instead of down to 0.36. A reader cross-referencing sections will notice.

**Suggested fix:**  
Standardize to 0.36 (rounded from 0.356) in all sections, including sec6 (~line 25).

---

### M3. "Virulence force" used instead of "force of infection" (sec2, sec4)

**sec2_key_findings.tex (~line 62):**  
> "52% of the pathogen's virulence force still penetrates host defenses"

**sec4_results.tex (~line 55):**  
> "52% of the pathogen's virulence force still penetrates host defenses"

**sec6_discussion.tex (~line 70):**  
> "57% of the pathogen force of infection penetrates the host's defenses" *(correct terminology)*

**What's wrong:**  
Resistance (*r*) modifies the **infection hazard rate** by (1−*r*), not virulence. "Virulence force" conflates two distinct epidemiological concepts. The correct term is "force of infection" (as sec6 uses) or "infection hazard." Additionally, sec2/sec4 compute 52% from *r* = 0.478, while sec6 computes 57% from *r* = 0.43 — these are different *r* values for the same population (CA-C at 2050), making the discrepancy worse because it's both terminologically and numerically inconsistent.

**Suggested fix:**  
Replace "virulence force" with "force of infection" in sec2 and sec4. Align the *r* value used (see M1) so the percentage is consistent.

---

### M4. "Peak population" vs "initial population" — inconsistent terminology for 4,447,241

**sec1_introduction.tex (~line 85):** "peak historical population across all scenarios was 4,447,241"  
**sec2_key_findings.tex (~line 5):** "5.4% of the peak population of 4,447,241"  
**sec4_results.tex (~line 15):** "5.4% of the initial population of 4,447,241"  
**sec4 Table 1 header:** "Initial"

**What's wrong:**  
"Peak" and "initial" are used interchangeably. They happen to be the same number (the 2012 pre-disease population is both the starting point and the all-time maximum), but a reader may wonder whether these refer to different quantities.

**Suggested fix:**  
Choose one term and use it consistently. "Pre-disease peak" is the most informative — it conveys both that this is the maximum and that it occurred before SSWD. Update sec4 Table 1 header to match.

---

### M5. SSP abbreviation never spelled out (first use: sec1_introduction.tex, ~line 70)

**Text:**  
> "Moderate warming (SSP2-4.5)"

**What's wrong:**  
"Shared Socioeconomic Pathway" is never expanded. A conservation biologist or manager unfamiliar with IPCC scenarios may not know this acronym.

**Suggested fix:**  
On first use in sec1, write "Shared Socioeconomic Pathway 2-4.5 (SSP2-4.5)" and then use the abbreviation thereafter.

---

### M6. CMIP6 abbreviation never spelled out (first use: sec3_scenario_design.tex, ~line 45)

**Text:**  
> "a delta method applied to an ensemble of CMIP6 general circulation models"

**What's wrong:**  
"Coupled Model Intercomparison Project Phase 6" is never expanded. Same accessibility concern as M5.

**Suggested fix:**  
Spell out on first use.

---

## MINOR

### m1. GWAS never spelled out (sec5_model.tex, ~line 38; sec6_discussion.tex, ~line 100)

**Text:**  
> "inspired by GWAS results from…"

**What's wrong:**  
"Genome-wide association study" is not expanded. Most marine biologists will know the term, but conservation managers may not.

**Suggested fix:**  
Spell out on first use in sec5.

---

### m2. "Arrhenius-scaled" unexplained (sec5_model.tex, ~line 55)

**Text:**  
> "Transition rates between compartments are Arrhenius-scaled to local SST"

**What's wrong:**  
Arrhenius kinetics is a chemistry/physics concept. A biologist unfamiliar with temperature-dependent reaction rates may not understand this term.

**Suggested fix:**  
Add a brief parenthetical: "Arrhenius-scaled (i.e., exponentially increasing with temperature) to local SST."

---

### m3. Oregon T_vbnc rounding inconsistent across sections

| Section | OR T_vbnc value |
|---------|----------------|
| sec2 | 10.6°C |
| sec4 | 10.62°C |
| sec6 | ≈10.7°C |

**What's wrong:**  
10.6 and 10.62 are consistent, but 10.7 rounds in the wrong direction (10.62 → 10.6, not 10.7).

**Suggested fix:**  
Use 10.6°C consistently in approximate contexts (sec2, sec6), and 10.62°C when precision matters (sec4).

---

### m4. OISST not spelled out (sec3_scenario_design.tex, ~line 50)

**Text:**  
> "2012–2025 OISST climatology"

**What's wrong:**  
The full name "Optimum Interpolation SST" appears earlier in the same section, but "OISST" is used without introduction as an abbreviation.

**Suggested fix:**  
Either define "OISST" on first full-name usage or spell it out again here.

---

### m5. sec4 references "initial population" in regional table, while sec2 uses "peak population"

**sec4 Table 2 caption:**  
> "Regional recovery at 2050 (% of initial population surviving)"

**What's wrong:**  
Minor instance of the M4 issue. Should match whatever term is standardized.

---

### m6. sec6 uses \\citealt instead of \\citet (sec6_discussion.tex, ~line 45)

**Text:**  
> "\\citealt{gehman2025}"

**What's wrong:**  
`\citealt` omits parentheses, producing "Gehman 2025" without any delimiters. In context ("estimated ~50%; Gehman 2025"), this may be intentional, but `\citep` would be more standard for a parenthetical citation.

**Suggested fix:**  
Verify this is intentional. If not, use `\citep{gehman2025}`.

---

## SUGGESTIONS

### S1. Redundancy between Key Findings and Results

The Key Findings (sec2) and Results (sec4) sections contain substantial overlap in both data and mechanistic explanations. Examples:
- The "climate paradox" mechanism (warming helps north, hurts south) is explained in full in both sec2 §2.3 and sec4 §4.3.
- The "evolutionary arms race" is narrated in both sec2 §2.4 and sec4 §4.4 with nearly identical logic.
- The "virulence irrelevance" argument is made twice (sec2 §2.5, sec4 §4.5.4).

**Suggestion:** Sec2 could summarize results and defer mechanistic explanations to sec4, or sec4 could explicitly reference sec2's explanations instead of restating them. The current structure works for a standalone Key Findings reader, but creates noticeable repetition for anyone reading cover-to-cover.

---

### S2. Region code glossary

The 18 region codes (AK-AL, AK-WG, AK-OC, etc.) are introduced piecemeal across sections. A single reference table mapping codes to full geographic names, latitude ranges, and cluster membership would improve navigability — and would have prevented the AK-FN discrepancy (C2).

---

### S3. Consider a forward reference from Results to Model

The report's structure (Introduction → Key Findings → Scenario Design → Results → Model → Discussion) is effective for policy audiences, but sec4 uses terminology and concepts (VBNC sigmoid, Michaelis-Menten infection kinetics, sweepstakes reproduction) that are only explained in sec5. A brief note at the top of sec4 (similar to sec5's "Readers interested in conservation may proceed to Discussion") would help readers know that the Model section follows and will explain mechanics.

---

### S4. Known past errors — status

The following known past errors were checked and are **not present** in the current draft:

| Error | Status |
|-------|--------|
| "La Niña narrative" for 2022 bump | ✅ Not found — no La Niña claims about population bumps |
| "Immune exclusion" terminology | ✅ Not found — term not used anywhere |
| T_VBNC geography inverted | ✅ Correct — north evolves to 9°C floor, south unchanged |

---

### S5. The disease progression diagram omits Recovery

sec5 presents: S → E → I₁ → I₂ → Dead

But the recovery trait (*c*) allows individuals to clear infection and return to S. The diagram should include the R→S pathway to be complete. The text mentions it ("Cleared individuals return to susceptible status (R→S)"), but the displayed diagram doesn't show this branch.

---

### S6. "Candidate agent" hedging in sec5 may conflict with sec1's certainty

**sec5_model.tex (~line 48):**  
> "recently confirmed as a candidate agent of SSWD through fulfillment of Koch's postulates"

**sec1_introduction.tex (~line 30):**  
> "fulfilled Koch's postulates for the marine bacterium *Vibrio pectenicida*" (no hedging)

The sec5 phrasing "candidate agent" introduces uncertainty that sec1 does not. If Koch's postulates were fulfilled, it is the established etiological agent, not a "candidate." If the authors intend to hedge, this should be consistent across both sections.

---

## Summary

| Severity | Count |
|----------|-------|
| CRITICAL | 3 |
| MODERATE | 6 |
| MINOR | 6 |
| SUGGESTION | 6 |

**Most urgent fixes:** C1 (immunity language), C2 (AK-FN geography), C3 (congener taxonomy), M1 (CA-C resistance values), M3 (virulence force terminology).
