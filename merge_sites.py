#!/usr/bin/env python3
"""
Merge all regional site JSON files into a single all_sites.json with comprehensive stats.
"""

import json
import glob
import os
from pathlib import Path
from collections import defaultdict
import numpy as np

def main():
    # Find all site JSON files
    data_dir = Path("/home/starbot/.openclaw/workspace/sswd-evoepi/data/nodes")
    site_files = list(data_dir.glob("sites_*.json"))
    
    print(f"Found {len(site_files)} regional site files:")
    for f in sorted(site_files):
        print(f"  {f.name}")
    
    # Merge all sites
    all_sites = []
    region_stats = defaultdict(lambda: {"count": 0, "with_pycno": 0, "file_size": 0})
    
    for site_file in site_files:
        region = site_file.stem.replace("sites_", "")
        file_size = site_file.stat().st_size
        region_stats[region]["file_size"] = file_size
        
        try:
            with open(site_file, 'r') as f:
                sites = json.load(f)
            
            print(f"\n{region}: {len(sites)} sites ({file_size:,} bytes)")
            
            region_stats[region]["count"] = len(sites)
            
            # Check for Pycnopodia documentation
            pycno_count = 0
            for site in sites:
                # Add region to each site
                site["region"] = region
                
                # Check for Pycnopodia mentions
                has_pycno = False
                for field in ["species", "description", "notes", "taxa"]:
                    if field in site and site[field]:
                        if isinstance(site[field], str):
                            if "pycnopodia" in site[field].lower() or "sunflower" in site[field].lower():
                                has_pycno = True
                                break
                        elif isinstance(site[field], list):
                            for item in site[field]:
                                if isinstance(item, str) and ("pycnopodia" in item.lower() or "sunflower" in item.lower()):
                                    has_pycno = True
                                    break
                
                if has_pycno:
                    pycno_count += 1
                
                all_sites.append(site)
            
            region_stats[region]["with_pycno"] = pycno_count
            print(f"  → {pycno_count} sites with Pycnopodia documentation")
            
        except Exception as e:
            print(f"ERROR reading {site_file}: {e}")
            continue
    
    # Write merged file
    output_file = data_dir / "all_sites.json"
    with open(output_file, 'w') as f:
        json.dump(all_sites, f, indent=2)
    
    print(f"\n✅ Merged {len(all_sites)} total sites → {output_file}")
    print(f"File size: {output_file.stat().st_size:,} bytes")
    
    # Generate summary stats
    print("\n" + "="*60)
    print("REGIONAL BREAKDOWN:")
    print("="*60)
    total_pycno = 0
    for region in sorted(region_stats.keys()):
        stats = region_stats[region]
        print(f"{region:>6}: {stats['count']:>4} sites, {stats['with_pycno']:>3} with Pycnopodia ({stats['file_size']:>6,} bytes)")
        total_pycno += stats['with_pycno']
    
    print("-" * 60)
    print(f"TOTAL:  {len(all_sites):>4} sites, {total_pycno:>3} with Pycnopodia ({sum(s['file_size'] for s in region_stats.values()):>6,} bytes)")
    
    # Coordinate statistics
    if all_sites:
        lats = [s['lat'] for s in all_sites if 'lat' in s]
        lons = [s['lon'] for s in all_sites if 'lon' in s]
        
        print("\n" + "="*60)
        print("COORDINATE COVERAGE:")
        print("="*60)
        print(f"Latitude:  {min(lats):.3f}° to {max(lats):.3f}° ({max(lats)-min(lats):.1f}° span)")
        print(f"Longitude: {min(lons):.3f}° to {max(lons):.3f}° ({max(lons)-min(lons):.1f}° span)")
    
    return output_file, len(all_sites), total_pycno, region_stats

if __name__ == "__main__":
    main()