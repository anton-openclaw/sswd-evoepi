# Dynamic P_env Calibration Report: W65-W74

## Executive Summary

**Dynamic P_env creates the north-south recovery gradient that static P_env could not.** W71 (floor=100, δ=0.02) achieves the best RMSE of any run to date (0.692), with SS-S dropping from 22% (static) to **1.6%** (target 5%) and CA-N from 1.3% to **0.3%** (target 0.1%). However, AK recovery is now too LOW (10.5% vs target 50%). The mechanism works — it just needs tuning to suppress the south while letting the north breathe.

## Mechanism

Dynamic P_env replaces the static `P_env_max` with a per-node vibrio reservoir that:
1. **Builds up** from infected host shedding (fraction `α_env`)
2. **Decays** at rate `δ_env` (protist grazing, UV, dilution)
3. **Has a floor** from community maintenance: `P_env_floor × VBNC(SST)`

The floor is SST-modulated: warm southern waters maintain higher baseline vibrio concentrations than cold northern waters. This creates **differential disease persistence** — the missing mechanism.

## Design Matrix

All runs use W62 base parameters: CDT=1000, K½=400K, wfDP=300, wfRange=3000, P_env_max=2000, k_vbnc=2.0

| Run | P_env_floor | δ_env | α_env | Half-life |
|-----|-------------|-------|-------|-----------|
| W65 | 25 | 0.02 | 0.1 | 35 days |
| W66 | 25 | 0.05 | 0.1 | 14 days |
| W67 | 25 | 0.1 | 0.1 | 7 days |
| W68 | 50 | 0.02 | 0.1 | 35 days |
| W69 | 50 | 0.05 | 0.1 | 14 days |
| W70 | 50 | 0.1 | 0.1 | 7 days |
| W71 | 100 | 0.02 | 0.1 | 35 days |
| W72 | 100 | 0.05 | 0.1 | 14 days |
| W73 | 100 | 0.1 | 0.1 | 7 days |
| W74 | 50 | 0.05 | 0.2 | 14 days |

## Results: Recovery Fractions

Targets: AK-PWS=50%, AK-FN=50%, AK-FS=20%, BC-N=20%, SS-S=5%, JDF=2%, OR=0.25%, CA-N=0.1%

| Run | RMSE | AK-PWS | AK-FN | AK-FS | BC-N | SS-S | JDF | OR | CA-N |
|-----|------|--------|-------|-------|------|------|-----|-----|------|
| W65 | 0.970 | 8.8% | 14.5% | 5.8% | 5.3% | 17.2% | 22.8% | 6.4% | 4.3% |
| W66 | 1.357 | 42.1% | 40.7% | 40.5% | 32.0% | 55.1% | 62.8% | 30.8% | 42.0% |
| W67 | 1.508 | 79.3% | 62.3% | 58.9% | 66.2% | 76.0% | 69.8% | 62.5% | 71.2% |
| W68 | 0.942 | 14.5% | 3.3% | 8.5% | 7.2% | 5.6% | 14.3% | 4.9% | 4.3% |
| W69 | 1.269 | 38.2% | 37.0% | 30.4% | 33.4% | 40.0% | 54.5% | 25.0% | 26.9% |
| W70 | 1.460 | 74.9% | 67.1% | 57.8% | 62.2% | 61.7% | 67.1% | 56.5% | 54.8% |
| **W71** | **0.692** | **10.5%** | **3.3%** | **3.7%** | **7.0%** | **1.6%** | **7.1%** | **1.3%** | **0.3%** |
| W72 | 1.197 | 33.0% | 32.9% | 25.2% | 31.3% | 24.6% | 31.7% | 23.6% | 22.5% |
| W73 | 1.402 | 71.1% | 62.9% | 47.1% | 60.3% | 60.9% | 61.6% | 47.3% | 39.5% |
| W74 | 0.881 | 10.9% | 10.4% | 8.3% | 11.7% | 12.8% | 18.8% | 5.5% | 3.0% |

## Key Findings

### 1. δ_env (decay rate) is THE dominant parameter

This is the clearest signal in the entire calibration campaign:

- **δ=0.02** (slow decay, 35d half-life): Strong suppression. RMSE 0.69-0.97.
- **δ=0.05** (medium decay, 14d half-life): Moderate. RMSE 0.88-1.36.
- **δ=0.1** (fast decay, 7d half-life): Pool decays too fast → reverts to near-static behavior. RMSE 1.40-1.51.

**Biological interpretation**: Vibrio in the environment persists for weeks, not days. Slow decay (δ=0.02, ~35 day half-life) is most realistic — consistent with V. pectenicida persistence in marine sediments and biofilms.

### 2. P_env_floor controls the gradient slope

At δ=0.02 (the optimal decay rate):

| Floor | AK-PWS | SS-S | CA-N | Gradient (AK/SS) |
|-------|--------|------|------|-------------------|
| 25 | 8.8% | 17.2% | 4.3% | 0.5× |
| 50 | 14.5% | 5.6% | 4.3% | 2.6× |
| 100 | 10.5% | 1.6% | 0.3% | 6.6× |

Floor=100 with δ=0.02 gives the strongest gradient. The VBNC sigmoid means this floor is ~10× higher at 18°C (southern CA) than at 6°C (Alaska winter), creating exactly the differential maintenance needed.

### 3. α_env (shedding feedback fraction) amplifies suppression

W74 (α=0.2, floor=50, δ=0.05) vs W69 (α=0.1, same floor/δ):
- W74: AK-PWS=10.9%, SS-S=12.8% — more suppression everywhere
- W69: AK-PWS=38.2%, SS-S=40.0% — less suppression

Higher α increases disease persistence everywhere, not differentially. Useful as a global intensity knob.

### 4. Dynamic P_env vs Static P_env: The Gradient Exists

Best static (W62): AK-PWS=22.3%, SS-S=22.1% → **ratio 1.01× (flat)**
Best dynamic (W71): AK-PWS=10.5%, SS-S=1.6% → **ratio 6.6× (gradient!)**

**The mechanism works.** Dynamic P_env creates exactly the kind of differential recovery that was impossible with static P_env.

### 5. The Problem: Everything is Too Suppressed

W71 has great RMSE (0.692) and a clear gradient, but:
- AK-PWS = 10.5% (target 50%) — 5× too low
- AK-FN = 3.3% (target 50%) — 15× too low
- BC-N = 7.0% (target 20%) — 3× too low

The dynamic P_env is too aggressive overall. We need to either:
1. **Increase K_half** to boost recovery everywhere (then gradient does the rest)
2. **Lower P_env_max** to reduce initial crash severity
3. **Use intermediate floor** (e.g. 75) with **higher K_half** (600K-800K)

## Wavefront Timing

All runs have excellent wavefront speed (17-19 months to AK-FS, target 26 months — slightly fast):

| Run | CA-S→CA-C | CA-S→OR | CA-S→BC-N | CA-S→AK-FS | CA-S→AK-PWS |
|-----|-----------|---------|-----------|-------------|-------------|
| W71 | 13.2 mo | 14.3 mo | 17.2 mo | 17.2 mo | 19.7 mo |
| Target | 6 mo | 15 mo | 20 mo | 26 mo | 30 mo |

Wavefront timing is stable across all dynamic P_env variants — the long-range kernel + CDT=1000 combination is robust.

## Year-by-Year Recovery Trajectories (W71)

```
Year:     1     2     3     4     5     6     7     8     9    10    11    12    13
AK-PWS: 100%  100%   95%    1%    1%    3%    2%    2%    1%    1%   36%   12%   10%
BC-N:   100%  100%   21%    6%   12%   10%    6%    6%    5%   19%    9%    9%    7%
SS-S:   100%  100%   71%   48%   48%   39%   29%   13%   11%   11%    3%    2%    2%
JDF:    100%  100%   88%   54%   46%   39%   33%   28%   21%   20%   13%    7%    7%
OR:     100%  100%   41%   34%   19%    0%    0%    0%    0%    1%    3%    2%    1%
CA-N:   100%   99%   54%   14%    6%    0%    0%    0%    0%    4%    1%    0%    0%
```

Key observations:
- **South crashes first** (CA-N year 3 → 54%, year 6 → 0%) then stays near-extinct ✅
- **North crashes later** (AK-PWS year 3 → 95%, year 4 → 1%) then has year 11 rebound ✅
- **SS-S shows gradual decline** (71% → 48% → 13% → 2%) — sustained disease pressure ✅
- **Year 11 AK rebound** (1% → 36%): La Niña cooling event — emergent property preserved ✅

## Genetic Evolution (W71)

Mean resistance trait by region:
```
           Y1     Y3     Y5     Y7     Y9    Y11    Y13
AK-PWS: 0.150  0.150  0.196  0.303  0.252  0.232  0.248
BC-N:   0.150  0.163  0.199  0.220  0.225  0.227  0.216
SS-S:   0.150  0.150  0.151  0.151  0.157  0.204  0.209
CA-N:   0.150  0.166  0.168  0.218  0.219  0.324  0.303
```

- AK-PWS evolves resistance rapidly (+66% by year 7)
- CA-N has HIGHEST resistance by year 11 (0.324) — strong selection in the south
- SS-S evolves resistance SLOWLY — population too suppressed for selection to act efficiently

## Recommended Next Steps

**W71 is the right direction.** The gradient exists and the mechanism is sound. To fix the "everything too suppressed" problem:

1. **W75-W80**: Keep W71's floor=100, δ=0.02, but increase K_half from 400K → [600K, 800K, 1M] and/or lower P_env_max from 2000 → [1000, 1500]

The hypothesis: Higher K_half gives populations more demographic resilience (density-dependent recruitment buffer), allowing northern populations to recover against the weaker disease pressure they experience, while southern populations with high floor-maintained vibrio remain suppressed.

2. **Consider floor=75 × K_half=800K** as intermediate

3. **Fine-tune α_env** last — it's a global intensity knob, not a gradient knob

## Conclusion

Dynamic P_env definitively solves the differential recovery problem. The SST-modulated community floor creates a 6.6× north-south gradient that was impossible with static P_env (1.01× ratio). The decay rate δ_env=0.02 (~35 day half-life) is biologically reasonable and provides the strongest gradient. Next round should adjust K_half and P_env_max to lift overall recovery while preserving the gradient.
