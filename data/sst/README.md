# SST Data Directory

Sea surface temperature data for the SSWD-EvoEpi model's 11-node
stepping-stone network.

## Node Coordinates

| Node          | Latitude | Longitude  | Region                |
|---------------|----------|------------|-----------------------|
| Sitka         | 57.05°N  | 135.33°W   | Southeast Alaska      |
| Ketchikan     | 55.34°N  | 131.64°W   | Southeast Alaska      |
| Haida Gwaii   | 53.25°N  | 132.07°W   | British Columbia      |
| Bella Bella   | 52.16°N  | 128.14°W   | Central BC Coast      |
| Howe Sound    | 49.38°N  | 123.23°W   | Southern BC           |
| SJI           | 48.53°N  | 123.01°W   | San Juan Islands, WA  |
| Westport      | 46.89°N  | 124.10°W   | Washington Coast      |
| Newport       | 44.63°N  | 124.05°W   | Oregon Coast          |
| Crescent City | 41.74°N  | 124.18°W   | Northern California   |
| Fort Bragg    | 39.43°N  | 123.80°W   | Northern California   |
| Monterey      | 36.60°N  | 121.90°W   | Central California    |

## SST Sources

The model supports three SST forcing modes (`sst_source` in config):

### 1. Sinusoidal (default)
- `sst_source = 'sinusoidal'`
- Analytical sinusoidal annual cycle + linear warming trend
- No data files needed; parameters set in config
- Good for quick runs and sensitivity analysis

### 2. Satellite Climatology
- `sst_source = 'satellite'`
- Uses `*_climatology.csv` files (365-day mean annual cycle)
- Tiled across simulation years with optional linear trend
- Real spatial pattern, no interannual variability

### 3. Monthly (Year-Specific)
- `sst_source = 'monthly'`
- Real year-specific SST from NOAA OISST observations (2002–2025)
- CMIP6 projections for future years (2026+)
- Daily values interpolated from monthly means

## Observational Data

**Source:** NOAA Optimum Interpolation SST v2.1 (OISST)
- Resolution: 0.25° × 0.25°, daily → aggregated to monthly
- Period: January 2002 – December 2025
- Reference: Huang et al. (2021), doi:10.1175/JCLI-D-20-0166.1

### Files

| Pattern                        | Description                              |
|--------------------------------|------------------------------------------|
| `{Node}_monthly.csv`          | Monthly mean SST (year, month, sst)      |
| `{Node}_climatology.csv`      | 365-day climatology (day_of_year, sst_mean, sst_std, n_years) |
| `all_nodes_climatology.csv`   | Combined climatology for all 11 nodes    |
| `climatology_curves.png`      | Visual validation of climatology curves  |
| `validation_plots.png`        | OISST vs climatology comparison plots    |
| `fetch_summary.json`          | Metadata from the OISST download process |

## Projection Data

**Directory:** `projections/`

CMIP6 multi-model ensemble projections, bias-corrected against OISST.

### Available Scenarios

| Scenario | Files                          | Source           | Status      |
|----------|--------------------------------|------------------|-------------|
| SSP2-4.5 | `{Node}_ssp245_monthly.csv`   | CMIP6 ensemble   | ✅ Real data |
| SSP5-8.5 | `{Node}_ssp585_placeholder_monthly.csv` | Trend extrapolation | ⚠️ Placeholder |

### SSP2-4.5 (Real CMIP6 Data)

- **Models:** GFDL-ESM4, IPSL-CM6A-LR (multi-model ensemble mean)
- **Period:** 2026–2100 (monthly)
- **Bias correction:** Per-node, per-month offset computed from 2015–2025
  overlap between CMIP6 and OISST. See `projections/bias_correction_report.txt`.
- **Source:** Pangeo CMIP6 catalog on Google Cloud Storage
- **Processing script:** `scripts/process_cmip6_sst.py`

### SSP5-8.5 (Placeholder)

⚠️ **These are NOT real CMIP6 projections.** They are linear trend
extrapolations from OISST observations, scaled by 2.5× to approximate
the SSP5-8.5 warming rate relative to SSP2-4.5. Replace with real CMIP6
data when available (see instructions below).

### Raw CMIP6 Downloads

**Directory:** `cmip6/`

Raw monthly SST extracted at the 11 nearest grid points from CMIP6 models.
See `cmip6/README.md` for details on available models and download
instructions for additional data.

## Configuration

To use year-specific SST in a simulation:

```python
config.simulation.sst_source = 'monthly'
config.simulation.sst_start_year = 2013       # first calendar year of simulation
config.simulation.sst_scenario = 'ssp245'     # 'observed_only', 'ssp245', 'ssp585'
config.simulation.sst_data_dir = 'data/sst'   # where *_monthly.csv files are
config.simulation.sst_projection_dir = 'data/sst/projections'  # where projection CSVs are
config.simulation.sst_obs_end_year = 2025      # last year of OISST data
```

### Scenario Options

| Value            | Behavior                                              |
|------------------|-------------------------------------------------------|
| `observed_only`  | Use OISST for available years; climatology fallback beyond 2025 |
| `ssp245`         | OISST through 2025, then SSP2-4.5 projections         |
| `ssp585`         | OISST through 2025, then SSP5-8.5 (currently placeholder) |
| `ssp126`/`ssp370`| Config validates these but data not yet available      |

### Fallback Hierarchy

For any year in the simulation, the SST source is resolved in order:
1. **OISST observations** (if all 12 months available and year ≤ obs_end_year)
2. **CMIP6 projections** (if scenario ≠ 'observed_only' and projection data exists)
3. **Satellite climatology** (from `*_climatology.csv`)
4. **Repeat last complete year** (last resort)

## How to Add New Scenarios or Extend the Time Range

### Adding Real SSP5-8.5 Data

1. Download CMIP6 SSP5-8.5 `tos` (monthly) from [ESGF](https://esgf-node.llnl.gov/search/cmip6/)
   for GFDL-ESM4 and IPSL-CM6A-LR (or other models)
2. Extract to `data/sst/cmip6/` as CSV with columns: `time, Sitka, Ketchikan, ..., Monterey`
3. Run the processing pipeline:
   ```bash
   python scripts/process_cmip6_sst.py \
       --cmip6-dir data/sst/cmip6 \
       --oisst-dir data/sst \
       --output-dir data/sst/projections \
       --scenarios ssp585
   ```

### Adding More GCMs to the Ensemble

1. Download additional model data from ESGF or Pangeo
2. Save as `{Model}_{scenario}_tos_monthly.csv` in `data/sst/cmip6/`
3. Re-run `process_cmip6_sst.py` — it automatically includes all models
   found in the directory and produces an ensemble mean

### Extending OISST Observations

When new OISST data becomes available:
1. Download from [NOAA OISST](https://www.ncei.noaa.gov/products/optimum-interpolation-sst)
2. Extract monthly means for the 11 node coordinates
3. Append rows to each `{Node}_monthly.csv` (format: `year,month,sst`)
4. Update `sst_obs_end_year` in config to match the new end year
5. Consider re-running bias correction if the overlap period changes

## File Naming Conventions

| Pattern | Example | Description |
|---------|---------|-------------|
| `{Node}_monthly.csv` | `Sitka_monthly.csv` | OISST observations |
| `{Node}_climatology.csv` | `Sitka_climatology.csv` | 365-day climatology |
| `{Node}_{scenario}_monthly.csv` | `Sitka_ssp245_monthly.csv` | Real CMIP6 projection |
| `{Node}_{scenario}_placeholder_monthly.csv` | `Sitka_ssp585_placeholder_monthly.csv` | Trend-based placeholder |
| `{Model}_{scenario}_tos_monthly.csv` | `GFDL-ESM4_ssp245_tos_monthly.csv` | Raw CMIP6 download |

Node names use underscores for spaces (e.g., `Howe_Sound`, `Crescent_City`, `Fort_Bragg`).
`SJI` is the abbreviation for San Juan Islands.
