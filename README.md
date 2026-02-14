# SSWD-EvoEpi

**Coupled eco-evolutionary model of sea star wasting disease**

A spatially explicit, individual-based model (IBM) for simulating the dynamics of
sea star wasting disease (SSWD) in *Pycnopodia helianthoides* (sunflower sea star),
coupling:

- **Population dynamics** with sweepstakes reproductive success (SRS) and Allee effects
- **Disease dynamics** (SEIPD+R compartments for *Vibrio pectenicida* / SSWD)
- **Polygenic resistance evolution** (51 additive + 1 overdominant locus)
- **Metapopulation connectivity** via larval and pathogen dispersal
- **Conservation interventions** (captive breeding, assisted gene flow, release scenarios)

## Authors

- **Anton** 🔬 (AI research assistant)
- **Willem Weertman** — PhD candidate, UW Psychology / Neural Systems & Behavior / Sea Star Lab

## Status

🚧 **Active development** — Phase 0 infrastructure complete.

## Quick Start

```bash
# Install dependencies
pip install numpy scipy pyyaml pytest

# Run tests
python -m pytest tests/ -v

# Load config
python -c "from sswd_evoepi.config import load_config; c = load_config('configs/default.yaml'); print(c)"
```

## Project Structure

```
sswd-evoepi/
├── sswd_evoepi/          # Source code
│   ├── types.py          # Core data types (AGENT_DTYPE, enums)
│   ├── config.py         # Configuration system (YAML loading)
│   ├── rng.py            # Seeded RNG hierarchy
│   ├── population.py     # Population dynamics (stub)
│   ├── disease.py        # Disease dynamics (stub)
│   ├── genetics.py       # Genetics & evolution (stub)
│   ├── spatial.py        # Spatial connectivity (stub)
│   ├── conservation.py   # Conservation interventions (stub)
│   ├── environment.py    # Environmental forcing (stub)
│   ├── recorder.py       # Data recording (stub)
│   ├── model.py          # Simulation orchestrator (stub)
│   └── utils.py          # Utility functions
├── configs/
│   └── default.yaml      # Default parameters (all modules)
├── tests/                # Test suite
├── CODE_ERRATA.md        # Implementation errata tracker
└── README.md
```

## Design Documents

See `../sswd-literature/specs/` for detailed technical specifications:

- `integration-architecture-spec.md` — Master architecture document
- `population-dynamics-spec.md` — Population module spec
- `disease-module-spec.md` — Disease module spec
- `genetics-evolution-spec.md` — Genetics module spec
- `spatial-connectivity-spec.md` — Spatial module spec
- `conservation-module-spec.md` — Conservation module spec
- `data-parameterization-plan.md` — Parameter inventory & data pipeline

## Key Design Decisions

1. **No cost of resistance** — Removed per Willem's decision (CODE_ERRATA CE-1)
2. **Both etiological scenarios** — "ubiquitous" and "invasion" via config
3. **Exponential effect sizes** — Per Lotterhos & Whitlock 2016
4. **Individual-based genetics** — Diploid genotypes at 52 loci per agent
5. **SRS is the reproduction mechanism** — Not a post-hoc drift modifier

## License

CC-BY-4.0
