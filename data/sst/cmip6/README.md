# CMIP6 SST Projection Data

## Overview

CMIP6 (Coupled Model Intercomparison Project Phase 6) sea surface temperature
projections for the 11 stepping-stone network nodes used in SSWD-EvoEpi.

## Available Data

### Downloaded (Phase 3, Feb 24 2026)

| Model         | Scenario | Years     | Grid  | Source           |
|---------------|----------|-----------|-------|------------------|
| GFDL-ESM4     | SSP2-4.5 | 2015-2100 | gr    | Pangeo GCS       |
| IPSL-CM6A-LR  | SSP2-4.5 | 2015-2100 | gr    | Pangeo GCS       |

**Variable:** `tos` (sea surface temperature), monthly, already extracted to
the 11 nearest grid points and converted to °C.

### Target Models (not yet downloaded)

The following models were targeted but not available on Pangeo GCS with
regridded (gr) grids, or had issues with curvilinear (gn) extraction:

- MRI-ESM2-0, CanESM5, MIROC6, CESM2, UKESM1-0-LL, NorESM2-LM, EC-Earth3

**SSP5-8.5** was not available for the downloaded models on Pangeo GCS.

## Manual Download Instructions

For additional models or scenarios, use one of these sources:

### Option A: ESGF (recommended for completeness)
1. Go to https://esgf-node.llnl.gov/search/cmip6/
2. Search for: variable=tos, frequency=mon, experiment_id=ssp245/ssp585
3. Download NetCDF files
4. Use `scripts/process_cmip6_sst.py` to process them, or adapt
   `scripts/download_cmip6_sst.py` which has the extraction logic

### Option B: Pangeo Catalog (code-based)
The `scripts/download_cmip6_sst.py` script queries the Pangeo CMIP6 catalog
on Google Cloud Storage. Requires: `intake`, `intake-esm`, `xarray`, `gcsfs`.

### Option C: Climate Explorer
1. Go to https://climexp.knmi.nl/
2. Select CMIP6 → SST (tos) → desired model/scenario
3. Download time series for coordinates of interest

## Processing Pipeline

```bash
# Process raw CMIP6 CSVs into model-ready format with bias correction
python scripts/process_cmip6_sst.py \
    --cmip6-dir data/sst/cmip6 \
    --oisst-dir data/sst \
    --output-dir data/sst/projections \
    --scenarios ssp245,ssp585
```

This:
1. Computes per-node, per-month bias correction (CMIP6 − OISST, 2015-2025 overlap)
2. Multi-model ensemble mean (if multiple models available)
3. Applies bias correction to projection period (2026-2100)
4. Outputs per-node CSVs matching the model's `*_monthly.csv` schema
5. For missing scenarios, generates trend extrapolation placeholders

## File Format

Raw CMIP6 CSVs (this directory):
```
time,Sitka,Ketchikan,...,Monterey
2015-01-16 12:00:00,2.4924,1.7288,...,11.9395
```

Processed projections (`data/sst/projections/`):
```
year,month,sst
2026,1,8.2172
```

## Bias Correction

See `data/sst/projections/bias_correction_report.txt` for per-node,
per-model bias statistics. Mean biases range from -1.4°C to +3.8°C
depending on model and node. IPSL-CM6A-LR is warm-biased at southern
nodes; GFDL-ESM4 is relatively unbiased at northern nodes.

The per-month correction preserves the seasonal structure of bias
(summer vs winter offsets can differ substantially).

## Placeholder Data

SSP5-8.5 files labeled `*_ssp585_placeholder_monthly.csv` are **not**
real CMIP6 projections. They are linear trend extrapolations from OISST
observations scaled by 2.5× to approximate the SSP5-8.5 warming rate.
Replace with real CMIP6 data when available.
