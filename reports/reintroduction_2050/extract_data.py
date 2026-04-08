#!/usr/bin/env python3
"""Extract all reintroduction experiment data into data_summary.json."""

import json
import os
import sys
import numpy as np

REPO = '/home/starbot/.openclaw/workspace/sswd-evoepi'
REINTRO_DIR = f'{REPO}/results/reintro_2050'
BASELINE_DIR = f'{REPO}/results/forecast_2050/F01_baseline'
OUT_DIR = f'{REPO}/reports/reintroduction_2050/data'

# Experiment design
LOCATIONS = {'L1': 'Monterey, CA', 'L2': 'Friday Harbor, WA'}
RESISTANCE = {
    'R1': {'label': 'CA-2025 survivors', 'r_mean': 0.31},
    'R2': {'label': 'Pre-epidemic naive', 'r_mean': 0.15},
    'R3': {'label': 'Mid-epidemic 2019', 'r_mean': 0.20},
    'R4': {'label': 'Artificially selected', 'r_mean': 0.50},
    'R5': {'label': 'Immune', 'r_mean': 1.00},
}
DENSITY = {
    'D1': {'label': '50/site', 'n': 50},
    'D2': {'label': '200/site', 'n': 200},
    'D3': {'label': '500/site', 'n': 500},
    'D4': {'label': '1000/site', 'n': 1000},
}

REGIONS_ORDERED = [
    'BJ', 'CA-S', 'CA-C', 'CA-N', 'OR', 'WA-O', 'JDF', 'SS-S', 'SS-N',
    'BC-C', 'BC-N', 'AK-FS', 'AK-FN', 'AK-OC', 'AK-PWS', 'AK-EG', 'AK-WG', 'AK-AL'
]

def load_result(path):
    with open(path) as f:
        return json.load(f)

def extract_scenario(sc_name, sc_dir):
    """Extract key metrics from a scenario."""
    result = load_result(f'{sc_dir}/result_seed42.json')
    
    parts = sc_name.split('_')
    loc = parts[1]
    res = parts[2]
    dens = parts[3]
    
    rd = result['region_details']
    
    # Extract per-region data
    regions = {}
    for reg in REGIONS_ORDERED:
        if reg in rd:
            d = rd[reg]
            regions[reg] = {
                'recovery_frac': d.get('recovery_frac', 0),
                'recovery_pct': d.get('recovery_frac', 0) * 100,
                'final_pop': d.get('final_pop', 0),
                'peak_pop': d.get('peak_pop', 0),
                'crash_pct': d.get('crash_pct', 0),
                'n_nodes': d.get('n_nodes', 0),
                'yearly_totals': d.get('yearly_totals', []),
            }
            # Check for genetics data
            for gkey in ['yearly_mean_resistance', 'yearly_mean_tolerance', 
                        'yearly_mean_recovery', 'yearly_va_resistance']:
                if gkey in d:
                    regions[reg][gkey] = d[gkey]
    
    return {
        'scenario': sc_name,
        'location': loc,
        'location_name': LOCATIONS[loc],
        'resistance': res,
        'resistance_label': RESISTANCE[res]['label'],
        'r_mean': RESISTANCE[res]['r_mean'],
        'density': dens,
        'density_label': DENSITY[dens]['label'],
        'density_n': DENSITY[dens]['n'],
        'rmsle': result['scoring']['rmsle'],
        'final_pop_frac': result['overall']['final_pop_frac'],
        'final_pop_pct': result['overall']['final_pop_frac'] * 100,
        'pop_crash_pct': result['overall']['pop_crash_pct'],
        'wall_time_s': result.get('wall_time_seconds', 0),
        'regions': regions,
        'scoring': result['scoring'],
    }

def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    
    # Load baseline
    bl_result = load_result(f'{BASELINE_DIR}/result_seed42.json')
    bl_rd = bl_result['region_details']
    baseline = {
        'rmsle': bl_result['scoring']['rmsle'],
        'final_pop_frac': bl_result['overall']['final_pop_frac'],
        'final_pop_pct': bl_result['overall']['final_pop_frac'] * 100,
        'regions': {}
    }
    for reg in REGIONS_ORDERED:
        if reg in bl_rd:
            baseline['regions'][reg] = {
                'recovery_frac': bl_rd[reg].get('recovery_frac', 0),
                'recovery_pct': bl_rd[reg].get('recovery_frac', 0) * 100,
                'final_pop': bl_rd[reg].get('final_pop', 0),
                'peak_pop': bl_rd[reg].get('peak_pop', 0),
                'yearly_totals': bl_rd[reg].get('yearly_totals', []),
            }
    
    # Load all scenarios
    scenarios = {}
    sc_dirs = sorted([d for d in os.listdir(REINTRO_DIR) if d.startswith('REINTRO_')])
    
    for sc_name in sc_dirs:
        sc_dir = f'{REINTRO_DIR}/{sc_name}'
        try:
            sc = extract_scenario(sc_name, sc_dir)
            scenarios[sc_name] = sc
        except Exception as e:
            print(f'ERROR loading {sc_name}: {e}', file=sys.stderr)
    
    # Check for genetics data
    genetics_available = False
    sample = list(scenarios.values())[0] if scenarios else None
    if sample:
        sample_reg = sample['regions'].get('CA-C', {})
        genetics_available = 'yearly_mean_resistance' in sample_reg
    
    # Compute summary tables
    # L1 heatmap: resistance × density → CA recovery metrics
    l1_heatmap = {}
    l2_heatmap = {}
    for sc in scenarios.values():
        key = f"{sc['resistance']}_{sc['density']}"
        entry = {
            'ca_s_pct': sc['regions'].get('CA-S', {}).get('recovery_pct', 0),
            'ca_c_pct': sc['regions'].get('CA-C', {}).get('recovery_pct', 0),
            'ca_n_pct': sc['regions'].get('CA-N', {}).get('recovery_pct', 0),
            'ca_total_pct': (
                sc['regions'].get('CA-S', {}).get('recovery_pct', 0) +
                sc['regions'].get('CA-C', {}).get('recovery_pct', 0) +
                sc['regions'].get('CA-N', {}).get('recovery_pct', 0)
            ) / 3,
            'jdf_pct': sc['regions'].get('JDF', {}).get('recovery_pct', 0),
            'or_pct': sc['regions'].get('OR', {}).get('recovery_pct', 0),
            'ak_pws_pct': sc['regions'].get('AK-PWS', {}).get('recovery_pct', 0),
            'bc_n_pct': sc['regions'].get('BC-N', {}).get('recovery_pct', 0),
            'final_pop_pct': sc['final_pop_pct'],
            'rmsle': sc['rmsle'],
        }
        if sc['location'] == 'L1':
            l1_heatmap[key] = entry
        else:
            l2_heatmap[key] = entry
    
    # Competitive dilution analysis
    dilution = {}
    for res in ['R1', 'R2', 'R3']:
        dilution[res] = {}
        for dens in ['D1', 'D2', 'D3', 'D4']:
            key = f'REINTRO_L1_{res}_{dens}'
            if key in scenarios:
                sc = scenarios[key]
                dilution[res][dens] = {
                    'ca_c_pct': sc['regions'].get('CA-C', {}).get('recovery_pct', 0),
                    'ca_n_pct': sc['regions'].get('CA-N', {}).get('recovery_pct', 0),
                    'density_n': sc['density_n'],
                }
    
    # Best scenarios
    best_l1_ca = max(
        [sc for sc in scenarios.values() if sc['location'] == 'L1'],
        key=lambda s: s['regions'].get('CA-C', {}).get('recovery_pct', 0)
    )
    best_l2_jdf = max(
        [sc for sc in scenarios.values() if sc['location'] == 'L2'],
        key=lambda s: s['regions'].get('JDF', {}).get('recovery_pct', 0)
    )
    
    summary = {
        'title': 'Reintroduction Experiment 2050',
        'n_scenarios': len(scenarios),
        'seed': 42,
        'K': 1000,
        'K_ref': 5000,
        'years': 38,
        'start_year': 2002,
        'release_year': 2025,
        'end_year': 2050,
        'genetics_available': genetics_available,
        'baseline': baseline,
        'scenarios': scenarios,
        'l1_heatmap': l1_heatmap,
        'l2_heatmap': l2_heatmap,
        'competitive_dilution': dilution,
        'best_l1_ca': {
            'scenario': best_l1_ca['scenario'],
            'ca_c_pct': best_l1_ca['regions'].get('CA-C', {}).get('recovery_pct', 0),
            'ca_n_pct': best_l1_ca['regions'].get('CA-N', {}).get('recovery_pct', 0),
        },
        'best_l2_jdf': {
            'scenario': best_l2_jdf['scenario'],
            'jdf_pct': best_l2_jdf['regions'].get('JDF', {}).get('recovery_pct', 0),
        },
        'design': {
            'locations': LOCATIONS,
            'resistance': RESISTANCE,
            'density': DENSITY,
            'regions_ordered': REGIONS_ORDERED,
        },
    }
    
    # Remove yearly_totals from scenarios to keep JSON manageable
    # (they're still in the full data but not in the summary)
    # Actually, keep them — they're needed for figure generation
    
    out_path = f'{OUT_DIR}/data_summary.json'
    
    # Custom serializer for numpy types
    class NpEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, np.integer):
                return int(obj)
            if isinstance(obj, np.floating):
                return float(obj)
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            return super().default(obj)
    
    with open(out_path, 'w') as f:
        json.dump(summary, f, indent=2, cls=NpEncoder)
    
    print(f'Wrote {out_path} ({os.path.getsize(out_path) / 1024:.0f} KB)')
    print(f'Scenarios: {len(scenarios)}')
    print(f'Genetics data: {genetics_available}')
    print(f'Best L1 (CA): {best_l1_ca["scenario"]} → CA-C {best_l1_ca["regions"]["CA-C"]["recovery_pct"]:.1f}%')
    print(f'Best L2 (JDF): {best_l2_jdf["scenario"]} → JDF {best_l2_jdf["regions"]["JDF"]["recovery_pct"]:.1f}%')

if __name__ == '__main__':
    main()
