# SSWD-EvoEpi Build Review — Phase 14

**Date:** 2026-02-14  
**Reviewer:** Anton 🔬  
**Scope:** Full codebase (15 modules, 8 test files), CODE_ERRATA (CE-1 through CE-15), specs ERRATA (E1–E14, I1–I2), 5-node prototype results

---

## 1. What Works Well

### 1.1 Disease Module (disease.py) — ★★★★☆

The SEIPD+R implementation is the strongest module in the build. It correctly captures:

- **Temperature-dependent dynamics:** Arrhenius scaling produces biologically correct R₀ gradient (sub-threshold at 8°C, epidemic at 16°C). The ERRATA E1 correction (E_a/R for I₂→D from 6000 to 2000 K) was critical and properly implemented.
- **Force of infection:** The multiplicative λ_i = a × P/(K_half+P) × (1−r_i) × S_sal × f_size(L_i) architecture handles four susceptibility modifiers cleanly. Vectorized S→E transitions avoid per-agent Python loops for the most common operation.
- **Erlang-distributed stage durations:** Using countdown timers with shape k=2–3 produces realistic variance in progression times. This is better than fixed duration or exponential (memoryless) alternatives.
- **Dual etiological scenarios:** The ubiquitous/invasion config switch (CE-2) provides experimental flexibility.
- **Environmental pathogen dynamics:** The VBNC sigmoid, temperature-dependent decay (counter-intuitive: faster decay at cold T), and flushing all interact correctly.
- **R₀ computation:** Now includes carcass saprophytic contribution (CE-6 fix was critical — σ_D from 150 to 15 eliminated cold-temperature epidemic artifact).

**Patterns matching literature:**
- Latitudinal mortality gradient (Monterey 99.8% > Sitka 86%)
- Epidemic onset timing (within 1 year of introduction, matching 2013–2015 timeline)
- Peak prevalences in realistic range (5.7–66.7%)
- Continued population decline 7+ years post-epidemic (Hamilton 2021)

### 1.2 Genetics Module (genetics.py) — ★★★★☆

- **52-locus architecture** cleanly separates 51 additive + 1 overdominant (EF1A).
- **Exponential effect sizes** (CE-3, Lotterhos & Whitlock 2016) produce realistic heavy-tailed distribution — most loci contribute little, few contribute much.
- **Vectorized r_i computation** (compute_resistance_batch) is efficient: allele sum → dot product with effects → add overdominant bonus.
- **EF1A lethal elimination** correctly prevents ins/ins homozygotes from surviving (both at initialization and during reproduction).
- **Allele frequency tracking, V_A, heterozygosity, F_ST** — comprehensive diagnostics are computed per node per year. This is excellent for validation.
- **Selection signal detected:** 3/5 nodes show positive Δr̄ post-epidemic; direction is consistent even when magnitude is small (CE-10).
- **Genotype bank (Tier 2):** Implementation with SRS-weighted expansion (E7) is architecturally sound for future 150-node scaling.

### 1.3 Reproduction Module (reproduction.py) — ★★★★☆

- **SRS lottery with Pareto weights** — the core innovation. Correctly implements Hedgecock's sweepstakes: raw Pareto(α=1.35) weights × fecundity quality, then multinomial parent sampling. Produces realistically low Ne/N.
- **Mendelian inheritance** is clean and correct: per-locus, per-copy random choice from diploid parents.
- **Allee effects** work on three levels: fertilization kinetics (Lundquist & Botsford 2004 quadratic), settlement cue (Michaelis-Menten adult biofilm), and Beverton-Holt density-dependent recruitment. All three correctly suppress reproduction at low density.
- **BH formula correction** (CE-4) was the second-most critical bug fix: the original spec formula gave zero recruits.

### 1.4 Spatial Module (spatial.py) — ★★★☆☆

- **Two-matrix architecture** (C larval, D pathogen) correctly captures different dispersal scales (E5).
- **Haversine with tortuosity** — simple but functional coastline approximation.
- **Sill attenuation** for fjord pathogen barriers is a nice detail.
- **Network builder** cleanly constructs the full MetapopulationNetwork from node definitions.
- **5-node test network** is well-chosen: spans the latitudinal range, includes fjord (Howe Sound), semi-enclosed (SJI), and open coast.

### 1.5 Infrastructure

- **Configuration system** (config.py) — hierarchical YAML merge, dataclass sections, validation. Clean and extensible.
- **RNG hierarchy** (rng.py) — SeedSequence spawning ensures per-node independence and bit-exact replay. Excellent for reproducibility.
- **Types module** (types.py) — single source of truth for AGENT_DTYPE, enums, constants. Good architecture.
- **Full integration** (model.py) — 1,281 lines, but the main simulation loop is readable and correctly sequences daily disease → annual demographics → genetics tracking.

### 1.6 Verification Results

- **36/36 checks passed** in the final 5-node 20-year run with only 4 warnings.
- **21.5-second runtime** for a 5-node 20-year simulation — performant enough for development iteration.

---

## 2. What's Concerning

### 2.1 Monotonic Population Decline — No Recovery Phase

**Severity: HIGH**  
**Affects: All production run interpretations**

The most important behavioral concern. After SSWD introduction, all 5 nodes show monotonic decline with zero recovery inflection over 15 years. This is arguably consistent with field observations (Hamilton 2021), but it means the model has **no mechanism for population recovery** without external intervention.

The problem: endemic Vibrio + warm summers + Allee effects create a stable extinction vortex. Even as resistance evolves upward (Δr̄ ≈ +0.02), the population is too small for reproduction to exceed mortality. This is the core question the model is designed to answer, so it's working as intended — but it means 100+ year runs may just show slow extinction at all nodes.

**What could break this:** If the temperature-dependent force of infection is even slightly too high (e.g., a_exposure = 0.75 is too aggressive), populations can never escape the vortex. The parameter a_exposure is ★☆☆ confidence.

### 2.2 Per-Locus Allele Frequency Shifts Too Small (CE-10)

**Severity: MEDIUM**  
**Affects: Calibration against Schiebelhut 2018**

Schiebelhut 2018 reported Δq ≈ 0.08–0.15 at outlier loci after ~81% mortality in *Pisaster*. Our model produces Δq ≈ 0.01–0.03 at the top locus. This 3–5× discrepancy likely means:

1. The real genetic architecture has **fewer loci with larger effects** (10–20 loci rather than 51)
2. Linkage effects (hitchhiking) amplify apparent per-locus shifts in real data
3. SRS drift at small N masks the directional signal in our model

The phenotypic response (mean r̄ shift) is consistent with strong selection, so the model captures the *aggregate* signal. But if reviewers compare per-locus values to Schiebelhut, there's a mismatch.

### 2.3 Lab-to-Field Shedding Rate Scaling (ERRATA E2)

**Severity: MEDIUM**  
**Affects: Epidemic dynamics, R₀ calibration**

The ~3,000× lab-to-field scaling factor for shedding rates is the most uncertain number in the model. We have:
- Lab values (Lupo 2020): σ₂ ≈ 167,500 bact/mL/d/host in 25L tanks
- Field-effective: σ₂_eff = 50 bact/mL/d/host (calibrated to R₀ ≈ 2 at 16°C)

The scaling depends on an assumed effective mixing volume of ~200,000 L per individual. If the real effective volume is 50,000 L (more clumped habitat), shedding rates should be 4× higher and R₀ would be 4× higher. This would dramatically change epidemic dynamics.

**No empirical data exists for field-effective Vibrio shedding in sea stars.** This is a first-principles estimate.

### 2.4 Python For-Loops in Hot Paths

**Severity: MEDIUM (performance only)**  
**Affects: Scalability to 150 nodes**

Several critical paths use Python for-loops over individuals:
- `daily_disease_update`: disease progression iterates `for idx in diseased_indices`
- `annual_natural_mortality`: iterates `for idx in alive_idx`
- `annual_growth_and_aging`: iterates `for idx in alive_idx`
- `initialize_population`: iterates `for i in range(n_individuals)`
- Mendelian inheritance: `for l in range(N_LOCI)` per offspring batch

At K=500, these are fast. At K=10,000 with 150 nodes, the daily disease loop would process ~1.5M individuals per day × 365 days × 100 years = ~55 billion iterations. **This will not run in pure Python.**

Numba was explicitly removed (CE-12) because it wasn't installed. For production runs, either:
- Install Numba and @njit the hot loops
- Port disease progression to vectorized NumPy
- Use Cython for the inner loop

### 2.5 Oversimplified Environmental Forcing

**Severity: LOW-MEDIUM**  
**Affects: Climate scenario realism**

SST is a deterministic sinusoid + linear trend. Real SST has:
- **Marine heatwaves** (blob events, 2013–2016, 2021): +2–4°C anomalies lasting months
- **ENSO modulation**: El Niño years warmer, La Niña cooler
- **PDO multidecadal oscillation**: shifts the baseline by ~1°C
- **Interannual variability**: ±1°C year-to-year noise

The absence of environmental stochasticity means the model predicts a smooth, deterministic outcome. Real populations experience stochastic "bad years" that could trigger extinction at nodes hovering near viability. Conversely, cool years could provide recovery windows.

### 2.6 Conservation Module — Stub Only

**Severity: LOW (for v1 prototype)**  
**Affects: Practical policy relevance**

conservation.py is a stub. The model can't yet simulate captive breeding releases, which is the primary conservation application. This is acknowledged as a future phase, but it's the module Willem and stakeholders care most about.

**Update from web search (Feb 2026):** Sunflower Star Laboratory conducted the first experimental outplanting in Monterey in December 2025 — 47/48 juveniles survived 4 weeks. Cal Academy plans first spawning in 2026. The conservation module should model these real protocols.

### 2.7 Pathogen Dispersal Matrix D Effectively Zero

**Severity: LOW (for 5-node prototype)**  
**Affects: Spatial epidemic dynamics**

All pairwise waterway distances in the 5-node network exceed 50 km, so D is all zeros (CE-13). This means disease introduction is simultaneous at all nodes (by design in the simulation), not wave-like. The 150-node network with ~10 km spacing will fix this, but it means we haven't actually tested pathogen dispersal mechanics in a realistic simulation.

---

## 3. Uncertainties That Matter

### 3.1 Parameters We're Least Confident About

| Parameter | Current Value | Confidence | Impact if Wrong |
|-----------|--------------|------------|-----------------|
| `a_exposure` (exposure rate) | 0.75 d⁻¹ | ★☆☆ | Changes R₀ linearly; wrong by 2× changes epidemic dynamics entirely |
| `sigma_1_eff`, `sigma_2_eff` (shedding) | 5, 50 bact/mL/d | ★☆☆ | Changes R₀ linearly via lab-field scaling |
| `sigma_D` (carcass shedding) | 15 bact/mL/d | ★☆☆ | CE-6 reduced from 150; still uncertain |
| `K_half` (half-infective dose) | 87,000 bact/mL | ★★☆ | From Lupo 2020 challenge experiments; species-specific? |
| `alpha_srs` (Pareto SRS shape) | 1.35 | ★★☆ | Controls Ne/N ratio; Árnason 2023 gives 1.35 but for different taxon |
| `rho_rec` (recovery rate) | 0.05 d⁻¹ | ★☆☆ | No empirical recovery data for SSWD |
| `gamma_fert` (fertilization kinetics) | 4.5 m² | ★★☆ | From Lundquist & Botsford; not measured for Pycnopodia |
| `P_env_max` (background Vibrio) | 500 bact/mL/d | ★☆☆ | No field measurements of environmental V. pectenicida |
| Self-recruitment fraction (α_self) | 0.10–0.30 | ★☆☆ | E10: no empirical data for any NE Pacific asteroid |

### 3.2 Structural Assumptions Under Question

1. **52 loci vs fewer**: Schiebelhut 2018 found only 3 outlier loci. If resistance is truly oligogenic (3–10 major loci), per-locus selection is much stronger, genetic drift matters less, and evolutionary rescue is faster. The 52-locus assumption distributes selection too thinly.

2. **SEIPD+R vs simpler SIR**: The full compartmental model has 6 disease states with multiple rate parameters. A simpler SIR or SI model with temperature-dependent mortality might capture 80% of the dynamics with 50% of the parameters and much less calibration uncertainty.

3. **Constant vs evolving pathogen**: Vibrio is assumed static. If pathogen virulence evolves (as it may over decades), the co-evolutionary dynamics change fundamentally. The Clement et al. 2024 Tasmanian devil model showed that pathogen attenuation is critical for host-disease coexistence.

4. **No spatial heterogeneity in habitat quality**: Each node has uniform habitat. Real coastlines have microhabitats with different exposure, depth, flow — creating within-node spatial refugia that could buffer populations.

5. **Annual reproduction pulse vs extended season**: The model uses a single spawning day per year. In reality, Pycnopodia may spawn over weeks, and multiple cohorts could establish. This affects the Ne bottleneck — multiple spawning events would increase effective population size.

6. **No urchin feedback**: Pycnopodia is a keystone predator of urchins. Urchin barrens (from Pycnopodia collapse) prevent kelp recovery, which reduces habitat quality for Pycnopodia settlers. This positive feedback loop is absent.

### 3.3 What Would Change Results Most if Wrong

Ranked by impact on the key output (does the population recover?):

1. **Exposure rate `a_exposure`**: Controls whether disease remains endemic post-epidemic. If 0.3 instead of 0.75, populations likely recover spontaneously in ~20 years.
2. **Number of resistance loci**: If 10 instead of 51, per-locus selection is 5× stronger, evolutionary rescue occurs faster, potentially within 10–20 generations.
3. **Recovery rate `rho_rec`**: If 0.10 instead of 0.05, twice as many infected individuals survive, maintaining population and accelerating resistance evolution.
4. **Lab-to-field shedding scaling**: If the effective volume is 100,000 L instead of 200,000 L, R₀ doubles and endemic disease becomes more severe.
5. **Self-recruitment α_self**: If 0.40 (high self-recruitment), depleted nodes can self-rescue; if 0.05 (low), recovery requires immigration from healthy nodes.

---

## 4. What Could Use More Information

### 4.1 Empirical Data That Would Constrain Uncertain Parameters

1. **Field-effective Vibrio shedding rates**: Environmental DNA (eDNA) measurements of V. pectenicida concentration around infected sea stars in situ would directly calibrate σ₁_eff, σ₂_eff. The Prentice 2025 paper identified the pathogen; quantitative eDNA studies should follow.

2. **Recovery rates from SSWD**: Any sea star that was observed with wasting symptoms and subsequently recovered would provide an empirical estimate of ρ_rec. MARINe monitoring data may contain this.

3. **Pycnopodia-specific genomics**: The reference genome exists (NCBI GCA_032158295.1, published 2024). GWAS for SSWD resistance in surviving Pycnopodia populations would tell us the actual genetic architecture (number of loci, effect sizes, EF1A analog presence/absence). This would replace our borrowed *Pisaster* parameters.

4. **Larval dispersal distances**: Particle tracking models (ROMS-derived) specific to Pycnopodia PLD (49–146 days) would provide empirical C matrix entries. Some biophysical models exist for other species in the California Current.

5. **In-situ Vibrio concentration dynamics**: Time-series of V. pectenicida concentration at field sites during an outbreak (if any occurs) would calibrate the entire environmental pathogen module — shedding, decay, flushing, VBNC dynamics.

6. **Density-dependent fertilization success in Pycnopodia**: Spawning trials at controlled densities would directly parameterize gamma_fert and the Allee threshold.

### 4.2 Methodological Papers That Address Our Implementation Questions

1. **Clement et al. 2024** — "Coevolution promotes the coexistence of Tasmanian devils and a fatal, transmissible cancer" (Evolution, qpae143). The closest methodological analog to SSWD-EvoEpi: individual-based eco-evolutionary model of host-disease dynamics with GWAS-parameterized polygenic resistance. **Paywalled. Priority scrape request.**

2. **Prowse et al. 2016** — "An efficient protocol for the global sensitivity analysis of stochastic ecological models" (Ecosphere). Directly addresses our sensitivity analysis challenge: LHS with stochastic IBMs. Available in Wiley Open Access.

3. **Cosandey-Godin et al. 2015** — Individual-based eco-evolutionary marine models review (Proc R Soc B). Framework for validation approaches in genetically-explicit marine IBMs.

4. **Alexander et al. 2014** — "Evolutionary rescue: linking theory for conservation and medicine" (Evol Appl). Key reference for whether evolutionary rescue is plausible in Pycnopodia given its life history (long generation time, high fecundity, SRS).

### 4.3 Validation Data We Don't Have

1. **Population trajectories at specific sites**: Pre-SSWD and post-SSWD census data at our 5 model nodes (or equivalent sites) would let us compare predicted vs observed crash magnitude and timing. MARINe has some of this data.

2. **Allele frequency shifts in Pycnopodia survivors**: No equivalent of Schiebelhut 2018 exists for sunflower stars yet. If pre-SSWD tissue samples exist (museum collections?), comparison with surviving individuals could directly validate our genetics module.

3. **Between-site connectivity**: Realized larval transport or adult migration rates between our model nodes are unknown. Population genetic structure (F_ST) across the range could be used as an indirect calibration target.

4. **Captive breeding survival data**: The December 2025 outplanting (47/48 survived 4 weeks in Monterey) provides the first empirical datapoint for the conservation module's release survival parameter.

---

## 5. Module-by-Module Assessment

| Module | File | Lines | Tests | Confidence | Notes |
|--------|------|-------|-------|-----------|-------|
| Types | types.py | 157 | test_types.py | ★★★★★ | Clean, complete, single source of truth |
| Config | config.py | 289 | test_config.py | ★★★★★ | Validation, YAML merge, all sections |
| RNG | rng.py | 105 | test_rng.py | ★★★★★ | SeedSequence hierarchy, checkpoint/restore |
| Disease | disease.py | 480 | test_disease.py | ★★★★☆ | Strong, but hot-loop performance concern |
| Genetics | genetics.py | 470 | test_genetics.py | ★★★★☆ | Clean architecture, good diagnostics |
| Reproduction | reproduction.py | 560 | test_reproduction.py | ★★★★☆ | SRS lottery works, BH corrected |
| Population | population.py | 8 | — | ★★☆☆☆ | Stub — lifecycle split across model.py |
| Environment | environment.py | 145 | — | ★★★☆☆ | Functional but oversimplified |
| Spatial | spatial.py | 520 | test_spatial.py | ★★★☆☆ | Good architecture, D matrix untested at scale |
| Model | model.py | 1281 | test_integration.py | ★★★★☆ | Large but readable; performance bottleneck |
| Conservation | conservation.py | 8 | — | ★☆☆☆☆ | Stub |
| Recorder | recorder.py | — | — | ★★★☆☆ | Snapshot recording works |
| Utils | utils.py | — | — | ★★★☆☆ | Helper functions |
| Snapshots | snapshots.py | — | — | ★★★☆☆ | Data capture |

---

## 6. Critical Errata: All Resolved, But Some Need Monitoring

The 15 CODE_ERRATA entries were all resolved during the build. The two most critical (CE-4 Beverton-Holt, CE-6 carcass shedding) caught fundamental numerical bugs that would have produced nonsensical results. This speaks well of the incremental testing approach.

However, several spec ERRATA remain ⚠️ Pending:
- E4: disease_timer countdown semantics (clarification only)
- E6: Bimodal r_i distribution (SD target should be 0.08–0.10, not 0.04)
- E7: Tier 2 → Tier 1 must use SRS-weighted sampling (implemented correctly, spec needs update)
- E8: Nonlinear coupling at Tier 2 (not yet relevant — all nodes are Tier 1)
- E9: Hood Canal not a salinity refugium (needs scenario testing)
- E10: Self-recruitment fraction uncertain (needs sensitivity analysis)
- E12: Captive-bred resistance gap (key scientific result)
- E13: AGENT_DTYPE origin field (implemented, spec needs update)
- I2: E-state spawning exclusion (implemented correctly, spec inconsistency)

---

## 7. Overall Assessment

**The SSWD-EvoEpi prototype is a scientifically sound, well-architected model that correctly captures the qualitative dynamics of SSWD in Pycnopodia.** The 15-module build with 36/36 passing verification checks is a strong foundation.

**What gives us confidence:**
- Temperature gradient matches field observations
- Fjord protection mechanisms work
- Selection signal is detectable despite SRS noise
- Critical bugs (BH formula, carcass feedback) were caught and fixed during build
- Code is readable, well-documented, and reproducible

**What keeps us up at night:**
- The monotonic decline with no recovery inflection — is this biological reality or parameter artifacts?
- Lab-to-field scaling of shedding rates is a 3,000× guess with no empirical anchor
- 52-locus architecture may distribute selection too thinly compared to real system
- Performance bottleneck at 150 nodes without compiled inner loops
- Conservation module is a stub — the most policy-relevant component is unbuilt

**Recommendation:** This is ready for Willem's review and discussion. The key conversation topics should be:
1. Is the no-recovery behavior biologically correct, or should we tune parameters to allow spontaneous recovery?
2. Should we test alternative genetic architectures (10, 20, 51 loci)?
3. What performance target do we need for production runs?
4. Priority of conservation module vs environmental stochasticity for next build phase?
