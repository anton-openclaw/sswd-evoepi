#!/usr/bin/env python3
"""Analyze wavefront timing from W257 results and W249 (baseline) results.

Key questions:
1. When does disease arrive at each region?
2. When do secondary infections first appear at origin nodes?
3. How does the wavefront propagation speed compare to observed SSWD spread?
"""

import sys
sys.path.insert(0, '/home/starbot/.openclaw/workspace/sswd-evoepi')

import json
import numpy as np

# Load W257 and W249 results
RESULT_DIR = '/home/starbot/.openclaw/workspace/sswd-evoepi/results/calibration/W249-W268'

for run_id in ['W249', 'W257']:
    print(f"\n{'='*80}")
    print(f"  {run_id} ANALYSIS")
    print(f"{'='*80}")
    
    # Load result JSON
    result_path = f'{RESULT_DIR}/{run_id}/result_seed42.json'
    try:
        with open(result_path) as f:
            result = json.load(f)
    except FileNotFoundError:
        # Try seed 137
        result_path = f'{RESULT_DIR}/{run_id}/result_seed137.json'
        with open(result_path) as f:
            result = json.load(f)
    
    print(f"  RMSE: {result.get('rmse', 'N/A')}")
    print(f"  Recoveries: {result.get('region_recoveries', {})}")
    
    # Load monthly NPZ
    npz_path = f'{RESULT_DIR}/{run_id}/monthly_seed42.npz'
    try:
        data = np.load(npz_path, allow_pickle=True)
    except FileNotFoundError:
        npz_path = f'{RESULT_DIR}/{run_id}/monthly_seed137.npz'
        data = np.load(npz_path, allow_pickle=True)
    
    print(f"\n  NPZ keys: {list(data.keys())}")
    
    populations = data['populations']  # (months, sites)
    infected = data['infected']  # (months, sites)
    sim_days = data['sim_days']
    site_names = data['site_names'] if 'site_names' in data else None
    
    print(f"  Shape: populations={populations.shape}, infected={infected.shape}")
    print(f"  Sim days range: {sim_days[0]} - {sim_days[-1]}")
    
    # Origin nodes
    origins = [322, 319, 632, 633, 634]
    
    # For each origin node, find first month with infections
    print(f"\n  Origin node disease onset:")
    for idx in origins:
        inf = infected[:, idx]
        first_month = np.argmax(inf > 0)
        if inf[first_month] > 0:
            print(f"    Node {idx}: first infections at month {first_month} (day ~{sim_days[first_month]}), count={inf[first_month]:.0f}")
        else:
            print(f"    Node {idx}: NO infections in simulation")
    
    # Regional analysis: find first month each region has infections
    # Need region mapping
    # Load site regions from the network
    if site_names is not None:
        # Group by region prefix
        from collections import defaultdict
        region_sites = defaultdict(list)
        for i, name in enumerate(site_names):
            name_str = str(name)
            # Extract region from site name pattern
            # Sites are like "CA-S_xxx", "OR_xxx", "WA_xxx", "BC-N_xxx", etc.
            parts = name_str.split('_', 1)
            region = parts[0] if len(parts) > 1 else name_str
            region_sites[region].append(i)
        
        print(f"\n  Regional wavefront arrival (first site with infections):")
        regions_ordered = []
        for region, sites in sorted(region_sites.items()):
            first_month = None
            first_site = None
            for s in sites:
                inf = infected[:, s]
                fm = np.argmax(inf > 0) if np.any(inf > 0) else None
                if fm is not None and inf[fm] > 0:
                    if first_month is None or fm < first_month:
                        first_month = fm
                        first_site = s
            
            if first_month is not None:
                regions_ordered.append((first_month, region, first_site, len(sites)))
                day = sim_days[first_month]
                month_count = first_month
                # Peak infection
                region_inf = infected[:, sites].sum(axis=1)
                peak_month = np.argmax(region_inf)
                peak_inf = region_inf[peak_month]
                print(f"    {region:10s} ({len(sites):3d} sites): arrival month {month_count:3d} (day ~{day:5d}), peak inf={peak_inf:.0f} at month {peak_month}")
            else:
                print(f"    {region:10s} ({len(sites):3d} sites): NEVER reached")
        
        # Wavefront speed
        if regions_ordered:
            regions_ordered.sort()
            print(f"\n  Wavefront progression (ordered by arrival):")
            for fm, reg, fs, ns in regions_ordered:
                print(f"    Month {fm:3d}: {reg}")
    
    # Check total infected over time
    total_inf = infected.sum(axis=1)
    total_pop = populations.sum(axis=1)
    print(f"\n  Coastwide dynamics:")
    for m in range(0, len(sim_days), 6):
        if m < len(total_inf):
            pct = total_inf[m] / total_pop[m] * 100 if total_pop[m] > 0 else 0
            print(f"    Month {m:3d} (day {sim_days[m]:5d}): pop={total_pop[m]:,.0f}, inf={total_inf[m]:,.0f} ({pct:.1f}%)")
