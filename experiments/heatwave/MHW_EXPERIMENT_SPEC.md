# MHW_EXPERIMENT_SPEC.md
## Marine Heatwave Experiment Technical Specification

### 1. Anomaly Extraction Methodology

#### 1.1 Baseline Climatology
- **Period**: 2002-2012 (11-year pre-Blob baseline)
- **Rationale**: Provides stable baseline excluding the Blob event while maximizing data availability
- **Computation**: For each node and month, calculate mean SST across 2002-2012

#### 1.2 Template Extraction
- **Template Period**: 2014-2016 (3-year lifecycle)
- **Components**:
  - 2014: Onset year (warming initiation)
  - 2015: Peak year (maximum intensity)
  - 2016: Dissipation year (gradual cooling)
- **Extraction Formula**: For each node `n`, month `m`, year `y`:
  ```
  anomaly[n,m,y] = SST[n,m,y] - baseline_climatology[n,m]
  ```

#### 1.3 Template Structure
Create a 36-month anomaly template capturing the full Blob lifecycle:
- Months 1-12: 2014 anomalies (onset)
- Months 13-24: 2015 anomalies (peak)
- Months 25-36: 2016 anomalies (dissipation)

### 2. Synthetic SST Generation Algorithm

#### 2.1 Temporal Envelope Options

**Option A: Full Lifecycle (36 months)**
- Preserves natural onset→peak→dissipation progression
- Most realistic but requires 3-year commitment in simulation time

**Option B: Compressed Lifecycle (24 months)**
- Month 1-6: Scaled 2014 onset (50% intensity)
- Month 7-18: Full 2015 peak (100% intensity)
- Month 19-24: Scaled 2016 dissipation (50% intensity)

**Option C: Peak-Only (12 months)**
- Use only 2015 peak year template
- Add 3-month linear ramp-up and ramp-down

**Recommendation**: Use Option B (24-month compressed) as default, with Option A for sensitivity analysis

#### 2.2 Intensity Scaling
Apply scaling factor `α` to all anomalies:
```
synthetic_SST[n,m] = baseline_SST[n,m] + α × anomaly_template[n,m]
```

Where:
- α = 0.5: Moderate heatwave (50% of Blob intensity)
- α = 1.0: Blob-equivalent heatwave
- α = 1.5: Extreme heatwave (150% of Blob intensity)

#### 2.3 Monthly Boundary Handling
- Use existing model's monthly→daily interpolation
- No additional smoothing needed (model handles this)

#### 2.4 Temperature Capping
Apply biological realism constraints:
```python
max_sst = 30.0  # °C, biological threshold
if synthetic_SST[n,m] > max_sst:
    synthetic_SST[n,m] = max_sst
```

### 3. Scenario Design Matrix

#### 3.1 Core Experimental Factors

**Factor 1: Intensity (3 levels)**
- Low: α = 0.5
- Medium: α = 1.0
- High: α = 1.5

**Factor 2: Frequency (4 levels)**
- Single event (year 2030)
- Every 10 years (2030, 2040, 2050)
- Every 7 years (2030, 2037, 2044)
- Every 5 years (2030, 2035, 2040, 2045, 2050)

**Factor 3: Duration (2 levels)**
- Standard: 24-month compressed lifecycle
- Extended: 36-month full lifecycle

#### 3.2 Scenario Allocation (Total: 36 scenarios)

**Primary Matrix (24 scenarios)**
- 3 intensities × 4 frequencies × 2 durations = 24 scenarios

**Additional Scenarios (12 scenarios)**
- 4 baseline/control runs (no heatwave, different random seeds)
- 4 "peak-only" runs (12-month duration, α = 1.0, different timings)
- 2 "compound event" runs (back-to-back heatwaves)
- 2 "ramping intensity" runs (increasing α over successive events)

### 4. Implementation Architecture

#### 4.1 Directory Structure
```
sswd-model/
├── data/
│   ├── sst/
│   │   ├── site_sst/              # Original OISST data
│   │   ├── site_sst_scenarios/    # Generated scenario data
│   │   │   ├── baseline/          # No-heatwave controls
│   │   │   ├── single_event/      # One-time heatwaves
│   │   │   ├── recurring/         # Periodic heatwaves
│   │   │   └── special/           # Compound/ramping scenarios
│   │   └── templates/             # Extracted anomaly templates
│   │       └── blob_anomalies_2014-2016.csv
├── scripts/
│   ├── extract_blob_template.py
│   ├── generate_mhw_scenarios.py
│   └── launch_mhw_experiments.py
└── configs/
    └── mhw_scenarios/
        ├── scenario_001_baseline.json
        ├── scenario_002_single_low.json
        └── ... (36 total)
```

#### 4.2 Generator Script Parameters
```python
# generate_mhw_scenarios.py
parameters = {
    "scenario_id": "recurring_high_5yr",
    "template_file": "templates/blob_anomalies_2014-2016.csv",
    "baseline_years": [2002, 2012],
    "events": [
        {"start_year": 2030, "intensity": 1.5, "duration_months": 24},
        {"start_year": 2035, "intensity": 1.5, "duration_months": 24},
        {"start_year": 2040, "intensity": 1.5, "duration_months": 24}
    ],
    "output_dir": "data/sst/site_sst_scenarios/recurring/high_5yr/"
}
```

#### 4.3 Config JSON Structure
```json
{
  "scenario_id": "recurring_high_5yr",
  "description": "High-intensity MHW every 5 years",
  "sst_data_dir": "data/sst/site_sst_scenarios/recurring/high_5yr/",
  "mhw_params": {
    "intensity_factor": 1.5,
    "frequency_years": 5,
    "duration_months": 24,
    "event_years": [2030, 2035, 2040, 2045, 2050]
  },
  "simulation": {
    "years": 40,
    "start_year": 2025,
    "K": 1000,
    "random_seed": 42
  }
}
```

#### 4.4 Launch Script Structure
```bash
#!/bin/bash
# launch_mhw_experiments.py

# Run 8 scenarios concurrently
parallel -j 8 python run_simulation.py --config {} ::: configs/mhw_scenarios/*.json

# Monitor and log progress
# Estimate: 36 scenarios × 2 hours / 8 parallel = 9 hours total
```

### 5. Edge Cases and Solutions

#### 5.1 Biologically Unrealistic Temperatures
- **Issue**: High-intensity scenarios in already-warm regions might exceed 30°C
- **Solution**: Apply temperature cap at 30°C (physiological limit for most marine life)
- **Implementation**: Check and cap in generator script after applying anomalies

#### 5.2 Salish Sea Buffering
- **Issue**: JDF/SS-S nodes show natural buffering in real data
- **Solution**: No special handling needed — the template already contains this spatial heterogeneity
- **Verification**: Confirm extracted anomalies show reduced magnitude for these nodes

#### 5.3 Post-2025 Climatology Interaction
- **Issue**: Model uses climatology repeat after 2025, potentially creating double-warming
- **Solution**: 
  1. Extract baseline climatology from same 2002-2012 period
  2. Apply anomalies to this consistent baseline
  3. Ensure generator uses same climatology calculation as main model

#### 5.4 Startup Transients
- **Issue**: Sudden heatwave onset might create model artifacts
- **Solution**: Always include 3-month linear ramp-up in synthetic events

### 6. Validation Checks

1. **Spatial Pattern Preservation**: Verify synthetic heatwaves maintain Blob spatial structure
2. **Magnitude Validation**: Ensure scaled intensities match expected ranges
3. **Temporal Smoothness**: Check monthly transitions are smooth after interpolation
4. **Biological Realism**: Verify no temperatures exceed 30°C threshold

### 7. Execution Timeline

1. **Week 1**: Extract Blob template, implement generator script
2. **Week 2**: Generate all scenario SST files, create config JSONs
3. **Week 3**: Run experiments (9-10 hours compute time)
4. **Week 4**: Initial analysis and validation

### 8. Code Snippets

#### 8.1 Template Extraction (Pseudocode)
```python
def extract_blob_template(sst_dir, baseline_start=2002, baseline_end=2012):
    # Load all node monthly data
    for site in sites:
        df = pd.read_csv(f"{sst_dir}/{site}_monthly.csv")
        
        # Calculate baseline climatology
        baseline = df[df.year.between(baseline_start, baseline_end)]
        climatology = baseline.groupby('month')['sst'].mean()
        
        # Extract 2014-2016 anomalies
        for year in [2014, 2015, 2016]:
            year_data = df[df.year == year]
            anomalies = year_data['sst'] - year_data['month'].map(climatology)
            save_anomalies(site, year, anomalies)
```

#### 8.2 Scenario Generation (Pseudocode)
```python
def generate_mhw_scenario(template, params, output_dir):
    for site in sites:
        # Load baseline SST
        baseline = load_baseline_sst(site)
        synthetic = baseline.copy()
        
        # Apply heatwave events
        for event in params['events']:
            start_idx = (event['start_year'] - 2025) * 12
            duration = event['duration_months']
            intensity = event['intensity']
            
            # Apply scaled anomalies
            for m in range(duration):
                if start_idx + m < len(synthetic):
                    anomaly = template[site][m % len(template)]
                    synthetic[start_idx + m] += intensity * anomaly
                    
        # Apply temperature cap
        synthetic = np.minimum(synthetic, 30.0)
        
        # Save to scenario directory
        save_synthetic_sst(synthetic, f"{output_dir}/{site}_monthly.csv")
```

### 9. Summary

This specification provides a complete framework for generating and executing marine heatwave experiments. The approach leverages the real Blob signature to create realistic synthetic events while maintaining computational feasibility within the 10-hour budget. The 36-scenario design balances comprehensive coverage of key factors with practical constraints.