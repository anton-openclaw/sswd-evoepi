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

🚧 **Active development** — Core model calibrated, 1020 tests passing.
Reintroduction and forecast experiments complete; MHW (marine heatwave)
experiments pending. All computation runs locally on a Ryzen 7 / 32 GB
workstation (see [INFRASTRUCTURE.md](INFRASTRUCTURE.md)).

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
│   ├── population.py     # Population dynamics
│   ├── disease.py        # Disease dynamics (SEIPD+R)
│   ├── genetics.py       # Polygenic resistance evolution
│   ├── spatial.py        # Spatial connectivity & dispersal
│   ├── conservation.py   # Conservation interventions
│   ├── environment.py    # Environmental forcing (SST, MHW)
│   ├── recorder.py       # Data recording & output
│   ├── model.py          # Simulation orchestrator
│   └── utils.py          # Utility functions
├── configs/
│   └── default.yaml      # Default parameters (all modules)
├── tests/                # Test suite (1020 tests)
├── CODE_ERRATA.md        # Implementation errata tracker
├── INFRASTRUCTURE.md     # Compute environment & benchmarks
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
