#!/usr/bin/env python3
"""
Create a site distribution map using matplotlib.
"""

import json
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from collections import defaultdict
import matplotlib.colors as mcolors

def create_site_map():
    # Load merged site data
    data_dir = Path("/home/starbot/.openclaw/workspace/sswd-evoepi/data/nodes")
    sites_file = data_dir / "all_sites.json"
    
    with open(sites_file, 'r') as f:
        all_sites = json.load(f)
    
    print(f"Loaded {len(all_sites)} sites for mapping")
    
    # Extract coordinates and regions
    lats = []
    lons = []
    regions = []
    pycno_status = []
    
    for site in all_sites:
        lat = site.get('latitude') or site.get('lat')
        lon = site.get('longitude') or site.get('lon')
        region = site.get('region', 'Unknown')
        has_pycno = any(key in site and site[key] for key in ['pycnopodia_documented', 'pycnopodia_evidence']) or \
                   any(field in site and site[field] and 
                       ('pycnopodia' in str(site[field]).lower() or 'sunflower' in str(site[field]).lower())
                       for field in ['notes', 'description', 'location_description', 'habitat_type'])
        
        if lat is not None and lon is not None:
            lats.append(float(lat))
            lons.append(float(lon))
            regions.append(region)
            pycno_status.append(has_pycno)
    
    # Create region color map
    unique_regions = sorted(list(set(regions)))
    colors = plt.cm.Set3(np.linspace(0, 1, len(unique_regions)))
    region_colors = {region: colors[i] for i, region in enumerate(unique_regions)}
    
    # Create figure with dark background
    plt.style.use('dark_background')
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
    fig.patch.set_facecolor('#0d1b2a')
    
    # Map 1: All sites by region
    ax1.set_facecolor('#0d1b2a')
    for region in unique_regions:
        region_lats = [lats[i] for i in range(len(lats)) if regions[i] == region]
        region_lons = [lons[i] for i in range(len(lons)) if regions[i] == region]
        ax1.scatter(region_lons, region_lats, c=[region_colors[region]], 
                   label=f'{region} ({len(region_lats)})', alpha=0.7, s=30)
    
    ax1.set_xlabel('Longitude', fontsize=12, color='white')
    ax1.set_ylabel('Latitude', fontsize=12, color='white')
    ax1.set_title('SSWD-EvoEpi Site Network by Region\\n489 Total Sites', fontsize=14, color='white', pad=20)
    ax1.grid(True, alpha=0.3)
    ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=9)
    
    # Map 2: Sites by Pycnopodia documentation
    ax2.set_facecolor('#0d1b2a')
    pycno_lats = [lats[i] for i in range(len(lats)) if pycno_status[i]]
    pycno_lons = [lons[i] for i in range(len(lons)) if pycno_status[i]]
    no_pycno_lats = [lats[i] for i in range(len(lats)) if not pycno_status[i]]
    no_pycno_lons = [lons[i] for i in range(len(lons)) if not pycno_status[i]]
    
    ax2.scatter(no_pycno_lons, no_pycno_lats, c='gray', alpha=0.5, s=20, 
               label=f'No Pycnopodia docs ({len(no_pycno_lats)})')
    ax2.scatter(pycno_lons, pycno_lats, c='gold', alpha=0.8, s=40, 
               label=f'Pycnopodia documented ({len(pycno_lats)})')
    
    ax2.set_xlabel('Longitude', fontsize=12, color='white')
    ax2.set_ylabel('Latitude', fontsize=12, color='white')
    ax2.set_title('Sites by Pycnopodia Documentation\\n287/489 Sites (58.7%)', fontsize=14, color='white', pad=20)
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=11)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save map
    output_file = data_dir / "site_map.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight', 
               facecolor='#0d1b2a', edgecolor='none')
    print(f"✅ Site map saved: {output_file}")
    
    plt.close()
    
    # Create summary statistics
    region_stats = defaultdict(lambda: {"total": 0, "with_pycno": 0})
    for i, region in enumerate(regions):
        region_stats[region]["total"] += 1
        if pycno_status[i]:
            region_stats[region]["with_pycno"] += 1
    
    print("\nSITE MAP SUMMARY:")
    print("-" * 50)
    print("Region     Total  Pycno  Coverage")
    print("-" * 50)
    for region in sorted(region_stats.keys()):
        stats = region_stats[region]
        coverage = stats["with_pycno"] / stats["total"] * 100 if stats["total"] > 0 else 0
        print(f"{region:>6}  {stats['total']:>6}  {stats['with_pycno']:>5}  {coverage:>6.1f}%")
    
    total_sites = len(lats)
    total_pycno = sum(pycno_status)
    total_coverage = total_pycno / total_sites * 100 if total_sites > 0 else 0
    print("-" * 50)
    print(f"TOTAL   {total_sites:>6}  {total_pycno:>5}  {total_coverage:>6.1f}%")
    
    return output_file, region_stats

if __name__ == "__main__":
    create_site_map()