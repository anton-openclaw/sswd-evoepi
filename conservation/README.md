# Conservation Genetics Module

Quantitative genetics framework for *Pycnopodia helianthoides* (sunflower sea star) conservation, built on the SSWD-EvoEpi model's three-trait genetic architecture.

## Purpose

Translate model predictions into actionable conservation tools:

1. **Predict** post-epidemic genetic state at wild populations
2. **Optimize** screening effort for founder selection
3. **Design** breeding programs balancing genetic gain and diversity
4. **Evaluate** reintroduction strategies via spatial simulation
5. **Synthesize** recommendations for practitioners

## Directory Structure

```
conservation/
├── README.md               ← You are here
├── analyses/               ← Analysis templates (run after calibration)
│   ├── params.yaml         ← Central configuration for all analyses
│   ├── 01_current_genetic_state.py
│   ├── 02_screening_effort.py
│   ├── 03_breeding_optimization.py
│   ├── 04_reintroduction_scenarios.py
│   └── 05_recommendations.py
├── src/                    ← Core library code (tested, validated)
│   ├── __init__.py
│   ├── trait_math.py       ← Analytical trait distributions (Eq. 2-4)
│   ├── screening.py        ← Screening theory & allocation (Eq. 4.1-4.10)
│   ├── breeding.py         ← Breeding simulation (Eq. 5.1-5.14)
│   ├── inbreeding.py       ← F coefficients, Ne, kinship (Eq. 6.2-6.8)
│   └── viz.py              ← Plotting functions
├── tests/                  ← Unit tests + validation reports
│   ├── test_trait_math.py
│   ├── test_screening.py
│   ├── test_breeding.py
│   ├── validation_trait_math.md
│   ├── validation_screening.md
│   └── validation_breeding.md
├── report/                 ← Theory report (LaTeX)
│   ├── main.tex
│   ├── main.pdf
│   ├── references.bib
│   ├── preamble.tex
│   └── sections/           ← Report sections (intro, theory, analysis plan)
├── results/                ← Analysis outputs (gitignored, generated)
├── ideas/                  ← Future directions parking lot
│   └── notes.md
└── specs/                  ← Design specs (in parent repo)
```

## How to Use

### 1. Run the tests (anytime)

```bash
cd sswd-evoepi
python -m pytest conservation/tests/ -v
```

### 2. Run analyses with default parameters (for pipeline testing)

```bash
cd conservation/analyses

# Analysis 1: Current genetic state (placeholder data)
python 01_current_genetic_state.py --seeds 5

# Analysis 2: Screening effort (requires Analysis 1)
python 02_screening_effort.py

# Analysis 3: Breeding optimization (requires Analysis 1)
python 03_breeding_optimization.py --founders 50 --generations 5 --replicates 10

# Analysis 4: Reintroduction scenarios (skeleton — dry run only)
python 04_reintroduction_scenarios.py --dry-run

# Analysis 5: Integrated recommendations (loads whatever results exist)
python 05_recommendations.py
```

### 3. Run with calibrated parameters (after ABC-SMC)

1. Complete ABC-SMC calibration (see main model README)
2. Update `analyses/params.yaml`:
   ```yaml
   calibrated_params_file: "../../results/calibration/abc_smc_posterior.json"
   ```
3. Rerun analyses 1-5 in order

### 4. Full production run (on Xeon)

```bash
# Analysis 1: ~4-8 hours (50 seeds × 11 sites × 13 years)
python 01_current_genetic_state.py --seeds 50

# Analysis 2: minutes
python 02_screening_effort.py

# Analysis 3: ~1-2 hours (4 strategies × 100 reps × 10 gen)
python 03_breeding_optimization.py --founders 200 --generations 10 --replicates 100

# Analysis 4: ~2-5 days (full scenario grid)
python 04_reintroduction_scenarios.py --seeds 50

# Analysis 5: minutes (synthesis)
python 05_recommendations.py
```

## Dependencies

- Python 3.10+
- NumPy, SciPy, matplotlib
- PyYAML
- sswd_evoepi model package (parent repository)

## Theory Report

The mathematical framework is documented in `report/main.pdf`:

- §2: Genetic architecture (additive model, Schiebelhut 51-locus basis)
- §3: Quantitative genetics (breeder's equation, selection response)
- §4: Screening theory (exceedance probabilities, optimal allocation)
- §5: Breeding programs (Mendelian crossing, selection schemes, OCS)
- §6: Inbreeding & diversity (genomic F, Ne, heterozygosity loss)
- §7: Reintroduction theory (release thresholds, allele spread)
- §8: Pycnopodia biology (life history, reproductive capacity)
- §9: Analysis plan (the analyses implemented here)

## Component Status

| Component | Code | Tests | Validation | Calibration |
|-----------|------|-------|------------|-------------|
| trait_math | ✅ | ✅ | ✅ | — |
| screening | ✅ | ✅ | ✅ | Needs calibrated allele frequencies |
| breeding | ✅ | ✅ | ✅ | Needs calibrated effect sizes |
| inbreeding | ✅ | ✅ | ✅ | — |
| viz | ✅ | — | — | — |
| Analysis 1 | ✅ template | — | — | **Blocked on ABC-SMC** |
| Analysis 2 | ✅ template | — | — | Blocked on Analysis 1 |
| Analysis 3 | ✅ template | — | — | Blocked on Analysis 1 |
| Analysis 4 | ⚠️ skeleton | — | — | Blocked on model extension + calibration |
| Analysis 5 | ✅ template | — | — | Blocked on Analyses 1-4 |

## Key Design Decisions

- **Uses actual model genetics** — not a separate implementation. Trait computations, allele frequencies, and Mendelian inheritance use `sswd_evoepi.genetics` directly.
- **Runnable without calibration** — all templates produce placeholder results with warnings, allowing pipeline validation before calibration is complete.
- **Central config** — `params.yaml` is the single source of truth for all analysis parameters. Change once, propagate everywhere.
- **Conservation report as spec** — the LaTeX theory report defines the equations; the code implements them with matching equation references.
