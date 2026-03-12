# K_cv Literature Review: Pycnopodia Density Variation Across Sites

## Published Pre-SSWD Density Estimates

### Montecino-Latorre et al. 2016 (PLoS ONE)
- **Salish Sea, strip transects, 6-18m depth**
- Pre-SSWD: 5.5 ind/100m² (= **0.055 ind/m²**)
- Post-SSWD: 0.18 ind/100m² (= **0.0018 ind/m²**) — 97% decline
- This is an average across San Juan Islands sites

### Casendino et al. 2023 (PLoS ONE)  
- **Port Madison, Puget Sound, trawl surveys, 10-70m depth**
- Washington outer coast: ~0.00025 ind/m² (converted from kg/ha, Harvell 2019 data)
- 20-year dataset (1997-2019), subtidal benthic
- High variability by depth: shallower = more Pycnopodia

### Gravem et al. 2021 (IUCN Red List / Supplement)
- **Figure 3**: Density map of Pycnopodia in shallow water (<25m)
- Historical (1976 to outbreak): densities shown by region
- Visually: highest in BC/AK, moderate in WA, very low in CA
- Cannot extract exact numbers from PDF, but supplement Table 5 has GLM results showing significant regional density differences
- "The species is most abundant in the waters off eastern Alaska and British Columbia"

### NOAA Status Review 2022 (Appendix A)
- "Very few density observations above zero individual per m²"
- Density varied with latitude AND depth
- Deep surveys (>25m, ROV): measured as ind/100m of transect, not area-based
- Hamilton et al. 2021 data: counts per site / total area surveyed

### Gehman et al. 2025 (Proc Royal Soc B)
- **BC fjords, biomass density (kg per 10m²)**
- High biomass in outer islands at higher temperatures
- High biomass in fjords at lower temperatures + higher salinity
- Suggests bimodal density distribution: productive open sites AND protected fjord sites

## What This Tells Us About K_cv

### Key observations
1. **Pre-SSWD densities varied by AT LEAST 1-2 orders of magnitude** across the range
   - AK/BC prime habitat: ~0.05-0.1 ind/m²
   - WA/OR intermediate: ~0.01-0.05 ind/m²  
   - CA marginal: ~0.001-0.01 ind/m²
   
2. **Within-region variation is also substantial**
   - Montecino-Latorre reports means across sites, but individual transects vary widely
   - Depth matters: shallow kelp forests >> deep subtidal
   - Habitat type matters: rocky reef/kelp >> sand/mud

3. **Marine benthic invertebrate densities typically follow lognormal distributions**
   - This is well-established in marine ecology (Taylor's power law)
   - CV of 0.5-1.5 is typical for sessile/slow-moving invertebrates

### Estimated CV from data
- If we take AK mean=0.055/m² and CA mean=0.005/m² as the full range:
  - Overall mean ≈ 0.03/m², overall SD ≈ 0.025/m²
  - CV ≈ 0.8
  
- But this confounds latitudinal gradient with site-level variation
- WITHIN a region (e.g., just Salish Sea), site-level CV is probably 0.5-0.8
- ACROSS the full range, including the latitudinal gradient: CV probably 0.8-1.2

### Recommended K_cv values for model
- **K_cv = 0.5**: Conservative. Mild site-level variation within regions.
- **K_cv = 0.8**: Moderate. Includes within-region habitat heterogeneity.
- **K_cv = 1.0**: Strong. Full range of variation including marginal sites.
- **K_cv > 1.0**: Extreme. Hard to justify without site-specific density data.

### Limitation
Our model applies K_cv uniformly — every region gets the same lognormal distribution.
In reality, the MEAN density varies by region (higher in AK, lower in CA).
A more realistic approach would be region-specific K_mean, but that's a bigger change.
K_cv captures within-region heterogeneity, which is the right first step.

## References
- Montecino-Latorre D et al. (2016) PLoS ONE 11(10): e0163190
- Casendino HR et al. (2023) PLoS ONE 18(6): e0286384
- Gravem SA et al. (2021) IUCN Red List Assessment + Supplement
- NOAA (2022) Status Review Report for Pycnopodia helianthoides
- Gehman AM et al. (2025) Proc Royal Soc B 292(2044): 20242770
