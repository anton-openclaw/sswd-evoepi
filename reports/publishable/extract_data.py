#!/usr/bin/env python3
"""
Extract all F01-F04 forecast data into a clean JSON summary for report generation.
Reads combined_results.json from each scenario and produces summary_data.json.
"""
import json
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import numpy as np
import os

BASE = '/home/starbot/.openclaw/workspace/sswd-evoepi/results/calibration'
SCENARIOS = {
    'F01': 'SSP2-4.5 Baseline (all evolution)',
    'F02': 'SSP2-4.5 No pathogen evolution',
    'F03': 'SSP2-4.5 Thermal adapt ON, virulence evo OFF',
    'F04': 'SSP5-8.5 Baseline (all evolution)',
}

# CENTRALIZED: moved to sswd_evoepi.results
from sswd_evoepi.results import REGION_ORDER
# REGION_ORDER = [
#     'AK-AL', 'AK-WG', 'AK-OC', 'AK-EG', 'AK-PWS', 'AK-FN', 'AK-FS',
#     'BC-N', 'BC-C', 'SS-N', 'SS-S', 'JDF', 'WA-O', 'OR', 'CA-N', 'CA-C', 'CA-S', 'BJ'
# ]

REGION_LABELS = {
    'AK-AL': 'Alaska - Aleutians',
    'AK-WG': 'Alaska - Western Gulf',
    'AK-OC': 'Alaska - Outer Coast',
    'AK-EG': 'Alaska - Eastern Gulf',
    'AK-PWS': 'Alaska - Prince William Sound',
    'AK-FN': 'Alaska - Fairweather North',
    'AK-FS': 'Alaska - Fairweather South',
    'BC-N': 'British Columbia - North',
    'BC-C': 'British Columbia - Central',
    'SS-N': 'Salish Sea - North',
    'SS-S': 'Salish Sea - South',
    'JDF': 'Juan de Fuca Strait',
    'WA-O': 'Washington - Outer Coast',
    'OR': 'Oregon',
    'CA-N': 'California - North',
    'CA-C': 'California - Central',
    'CA-S': 'California - South',
    'BJ': 'Baja California',
}

# Approximate representative latitudes for each region
REGION_LATS = {
    'AK-AL': 53.0, 'AK-WG': 57.5, 'AK-OC': 59.0, 'AK-EG': 59.5,
    'AK-PWS': 60.5, 'AK-FN': 58.5, 'AK-FS': 57.0,
    'BC-N': 54.0, 'BC-C': 51.0, 'SS-N': 49.0, 'SS-S': 48.0,
    'JDF': 48.3, 'WA-O': 47.0, 'OR': 44.5,
    'CA-N': 40.5, 'CA-C': 36.5, 'CA-S': 34.0, 'BJ': 30.0,
}

def load_scenario(scenario_id):
    """Load all seed results for a scenario."""
    path = os.path.join(BASE, scenario_id, 'combined_results.json')
    with open(path) as f:
        data = json.load(f)
    
    # Also load individual seed results for multi-seed averaging
    seed_results = []
    for seed in [42, 123, 999]:
        rpath = os.path.join(BASE, scenario_id, f'result_seed{seed}.json')
        if os.path.exists(rpath):
            with open(rpath) as f:
                seed_results.append(json.load(f))
    
    return data, seed_results

def extract_regional_summary(scenario_id, data, seed_results):
    """Extract per-region summary stats, averaged across seeds."""
    summary = {}
    
    for region in REGION_ORDER:
        region_data = {'region': region, 'label': REGION_LABELS[region], 'lat': REGION_LATS[region]}
        
        # Collect from all seeds
        recovery_fracs = []
        final_pops = []
        peak_pops = []
        crash_pcts = []
        
        # Yearly time series (we'll average across seeds)
        yearly_totals_all = []
        yearly_recruits_all = []
        yearly_deaths_all = []
        yearly_resistance_all = []
        yearly_tolerance_all = []
        yearly_recovery_trait_all = []
        yearly_va_resistance_all = []
        final_tvbnc_all = []
        final_vlocal_all = []
        
        for sr in seed_results:
            rd = sr.get('region_details', {}).get(region, {})
            if not rd:
                continue
            recovery_fracs.append(rd.get('recovery_frac', 0))
            final_pops.append(rd.get('final_pop', 0))
            peak_pops.append(rd.get('peak_pop', 0))
            crash_pcts.append(rd.get('crash_pct', 100))
            yearly_totals_all.append(rd.get('yearly_totals', []))
            yearly_recruits_all.append(rd.get('yearly_recruits', []))
            yearly_deaths_all.append(rd.get('yearly_disease_deaths', []))
            yearly_resistance_all.append(rd.get('yearly_mean_resistance', []))
            yearly_tolerance_all.append(rd.get('yearly_mean_tolerance', []))
            yearly_recovery_trait_all.append(rd.get('yearly_mean_recovery', []))
            yearly_va_resistance_all.append(rd.get('yearly_va_resistance', []))
            final_tvbnc_all.append(rd.get('final_mean_T_vbnc', 12.0))
            final_vlocal_all.append(rd.get('final_mean_v_local', 0.0))
        
        if not recovery_fracs:
            # Fallback to combined_results
            cr = data['results'][0]['region_details'].get(region, {})
            recovery_fracs = [cr.get('recovery_frac', 0)]
            final_pops = [cr.get('final_pop', 0)]
            peak_pops = [cr.get('peak_pop', 0)]
            crash_pcts = [cr.get('crash_pct', 100)]
            yearly_totals_all = [cr.get('yearly_totals', [])]
            yearly_resistance_all = [cr.get('yearly_mean_resistance', [])]
            yearly_va_resistance_all = [cr.get('yearly_va_resistance', [])]
            final_tvbnc_all = [cr.get('final_mean_T_vbnc', 12.0)]
            final_vlocal_all = [cr.get('final_mean_v_local', 0.0)]
        
        region_data['recovery_frac_mean'] = float(np.mean(recovery_fracs))
        region_data['recovery_frac_std'] = float(np.std(recovery_fracs))
        region_data['recovery_pct'] = float(np.mean(recovery_fracs) * 100)
        region_data['final_pop_mean'] = float(np.mean(final_pops))
        region_data['peak_pop_mean'] = float(np.mean(peak_pops))
        region_data['crash_pct_mean'] = float(np.mean(crash_pcts))
        
        # Average yearly time series
        if yearly_totals_all and yearly_totals_all[0]:
            n_years = len(yearly_totals_all[0])
            region_data['yearly_pop_mean'] = [
                float(np.mean([s[i] for s in yearly_totals_all if i < len(s)]))
                for i in range(n_years)
            ]
        
        if yearly_resistance_all and yearly_resistance_all[0]:
            n_years = len(yearly_resistance_all[0])
            region_data['yearly_resistance_mean'] = [
                float(np.mean([s[i] for s in yearly_resistance_all if i < len(s)]))
                for i in range(n_years)
            ]
        
        if yearly_va_resistance_all and yearly_va_resistance_all[0]:
            n_years = len(yearly_va_resistance_all[0])
            region_data['yearly_va_resistance_mean'] = [
                float(np.mean([s[i] for s in yearly_va_resistance_all if i < len(s)]))
                for i in range(n_years)
            ]
        
        region_data['final_T_vbnc_mean'] = float(np.mean(final_tvbnc_all))
        region_data['final_v_local_mean'] = float(np.mean(final_vlocal_all))
        
        summary[region] = region_data
    
    return summary

def compute_overall_stats(regional_summary):
    """Compute coast-wide aggregate stats."""
    total_peak = sum(r['peak_pop_mean'] for r in regional_summary.values())
    total_final = sum(r['final_pop_mean'] for r in regional_summary.values())
    overall_recovery = total_final / total_peak if total_peak > 0 else 0
    
    return {
        'total_peak_pop': total_peak,
        'total_final_pop': total_final,
        'overall_recovery_pct': overall_recovery * 100,
        'overall_crash_pct': (1 - overall_recovery) * 100,
    }

def main():
    all_data = {}
    
    for scenario_id, description in SCENARIOS.items():
        print(f"Processing {scenario_id}: {description}")
        data, seed_results = load_scenario(scenario_id)
        
        regional = extract_regional_summary(scenario_id, data, seed_results)
        overall = compute_overall_stats(regional)
        
        all_data[scenario_id] = {
            'description': description,
            'regions': regional,
            'overall': overall,
            'mean_rmse_log': data.get('mean_rmse_log', None),
            'years': data.get('years', 38),
            'K': data.get('K', 5000),
        }
        
        print(f"  Overall: {overall['overall_recovery_pct']:.1f}% surviving, "
              f"peak={overall['total_peak_pop']:.0f}, final={overall['total_final_pop']:.0f}")
        print(f"  Top 3 regions: ", end='')
        ranked = sorted(regional.items(), key=lambda x: x[1]['recovery_pct'], reverse=True)[:3]
        print(', '.join(f"{r}: {d['recovery_pct']:.1f}%" for r, d in ranked))
    
    # Cross-scenario comparisons
    comparisons = {}
    
    # F01 vs F04 (climate effect)
    climate_diff = {}
    for region in REGION_ORDER:
        f01_pct = all_data['F01']['regions'][region]['recovery_pct']
        f04_pct = all_data['F04']['regions'][region]['recovery_pct']
        climate_diff[region] = {
            'f01_pct': f01_pct,
            'f04_pct': f04_pct,
            'diff_pct': f04_pct - f01_pct,
            'direction': 'warming_helps' if f04_pct > f01_pct else 'warming_hurts',
        }
    comparisons['climate_effect'] = climate_diff
    
    # F01 vs F02 (pathogen evolution effect)
    evo_diff = {}
    for region in REGION_ORDER:
        f01_pct = all_data['F01']['regions'][region]['recovery_pct']
        f02_pct = all_data['F02']['regions'][region]['recovery_pct']
        evo_diff[region] = {
            'f01_pct': f01_pct,
            'f02_pct': f02_pct,
            'diff_pct': f02_pct - f01_pct,
            'direction': 'evo_off_helps' if f02_pct > f01_pct else 'evo_off_hurts',
        }
    comparisons['pathogen_evo_effect'] = evo_diff
    
    # F02 vs F03 (virulence evo effect)  
    vir_diff = {}
    for region in REGION_ORDER:
        f02_pct = all_data['F02']['regions'][region]['recovery_pct']
        f03_pct = all_data['F03']['regions'][region]['recovery_pct']
        vir_diff[region] = {
            'f02_pct': f02_pct,
            'f03_pct': f03_pct,
            'diff_pct': f03_pct - f02_pct,
        }
    comparisons['virulence_only_effect'] = vir_diff
    
    all_data['comparisons'] = comparisons
    
    # Save
    outpath = '/home/starbot/.openclaw/workspace/sswd-evoepi/reports/publishable/data/summary_data.json'
    with open(outpath, 'w') as f:
        json.dump(all_data, f, indent=2)
    print(f"\nSaved to {outpath}")

if __name__ == '__main__':
    main()
