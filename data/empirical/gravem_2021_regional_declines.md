# Gravem et al. 2021 — Regional Population Declines of Pycnopodia helianthoides

## Source
Gravem, S.A., W.N. Heady, V.R. Saccomanno, K.F. Alvstad, A.L.M. Gehman, T.N. Frierson and S.L. Hamilton. 2021. Pycnopodia helianthoides. IUCN Red List of Threatened Species. e.T178290276A197818455.

Supplement: `reef.org/sites/default/files/pycnopodia_helianthoides_published_supplement.pdf`
Assessment: `reefcheck.org/wp-content/uploads/2021/10/Final_IUCNAssessment.pdf`

## Data Sources
- ~61,000 scientific surveys from 31 datasets (Supplement Table 2)
- Shallow nearshore roving-diver surveys (< 25m)
- Deep offshore trawl surveys (55-1280m)
- Belt transect surveys (British Columbia, 3-18m)

## Global Summary
- **Overall decline: 90.6%** (estimated ~5.75 billion deaths)
- Species listed **Critically Endangered**
- Range halved: eradicated from ~2,700 km southern half

## Regional Declines (from IUCN Assessment text + NOAA Status Review + J. Heredity 2024)

### Southern Regions (Baja California to Washington) — 97.2% to 99.8%
| Region | Estimated Decline | Source/Notes |
|--------|------------------|--------------|
| Baja California | ~99-100% | Functionally extirpated |
| Southern California | ~99-100% | Functionally extirpated |
| Central California | ~99-100% | Functionally extirpated |
| Northern California | ~99-100% | Functionally extirpated |
| Oregon | 100% (deep trawl) | No individuals collected in 2016 trawls |
| Washington outer coast | 99.2% (deep trawl) | 3.11 → 0.02 kg/10ha |

### Northern Regions (British Columbia to Alaska) — >87.8%
| Region | Estimated Decline | Source/Notes |
|--------|------------------|--------------|
| Salish Sea | 91.9% | Gravem et al. 2021 main text |
| Southeast Alaska | 94.7% | Gravem et al. 2021 main text |
| Coastal British Columbia | 94.9% | Gravem et al. 2021 main text |
| Central British Columbia | ~96% (belt transects) | Harvell et al. 2019 (0.57-0.93 → 0.01-0.07 kg/10m²) |
| Northern British Columbia | ~87.8-94.9% | Two estimates reported (Gravem 2021) |
| West Gulf of Alaska | >87.8% | Part of 87.8% aggregate for northern range |
| East Gulf of Alaska | >87.8% | Part of 87.8% aggregate for northern range |
| Aleutian Islands | >87.8% | Crash timing: Jan 1, 2017 (latest of all regions) |

### Aggregate Summary
- **Northern range (Aleutians to Cape Flattery, WA)**: 87.8% mortality
- **Southern range (Cape Flattery south to Baja)**: 99-100% mortality
- Gravem (2024 interview): "60-80% decline in Alaska in the last decade"

## Crash Timing by Region (from Supplement Figure 3)
- **Jan 1, 2014**: California (northern, central, southern) and Baja California
- **Aug 1, 2015**: Southeast Alaska, coastal, northern and central BC, Salish Sea, Washington outer coast, Oregon
- **Jan 1, 2017**: Aleutians, west and east Gulf of Alaska

## Deep Water Data (Harvell et al. 2019, Science Advances)
- 8,968 bottom trawls, 55–1,280m depth, 2004–2016
- California: 100% biomass decline (2.78 → 0 kg/10ha)
- Oregon: 100% biomass decline (1.73 → 0 kg/10ha)
- Washington: 99.2% biomass decline (3.11 → 0.02 kg/10ha)
- 2016: zero individuals in 692 trawls covering 1,264 ha

## Notes for Model Calibration
These empirical declines represent the **initial crash** (2013-2017), which our model simulates
as the acute SSWD epidemic phase. Our calibration targets (set by Willem's expert estimates
for ~2024, 10 years post-crash) include partial recovery in Alaska:
- AK-PWS: 50% recovery → ~50% of pre-crash levels
- AK-FN: 50% → ~50% of pre-crash
- AK-FS: 20% → ~20% of pre-crash
- BC-N: 20% → ~20% of pre-crash
- SS-S: 5% → ~5% of pre-crash
- JDF: 2% → ~2% of pre-crash
- OR: 0.25% → near-extirpated
- CA-N: 0.1% → functionally extirpated

The Gravem data confirms the model should produce:
1. Near-total elimination (>99%) south of Washington
2. Severe but incomplete crash (~88-95%) from BC to Alaska
3. A clear latitudinal gradient in severity
4. Earlier crash in the south (2013-14) vs later in the north (2015-17)
5. No deep-water refuge

## Data Extraction Notes
- Supplement Table 2 (PDF) could not be extracted as text — binary PDF format
- Values compiled from: search result snippets, IUCN main assessment text, NOAA Status Review
  Appendix A, Harvell et al. 2019 (PMC6353623), J. Heredity reference genome paper (2024),
  Wikipedia summary, UAF news article (2024)
- **NOT extracted yet**: individual region values from Supplement Table 2 columns (pre-crash pop,
  post-crash pop, area, survey effort per region)
