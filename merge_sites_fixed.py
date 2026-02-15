#!/usr/bin/env python3
"""
Merge all regional site JSON files into a single all_sites.json with comprehensive stats.
Handles both JSON formats: array of sites vs object with sites array.
"""

import json
import glob
import os
from pathlib import Path
from collections import defaultdict
import numpy as np

def extract_sites_from_json(data, region):
    """Extract sites array from various JSON formats."""
    if isinstance(data, list):
        # Format 1: Direct array of sites (e.g., BC-N)
        return data
    elif isinstance(data, dict) and 'sites' in data:
        # Format 2: Object with sites array (e.g., AK-WG)
        return data['sites']
    else:
        print(f"WARNING: Unknown JSON format for {region}")
        return []

def standardize_site(site, region):
    """Standardize site dictionary with consistent field names."""
    # Create a copy to avoid modifying original
    standardized = dict(site)
    standardized["region"] = region
    
    # Standardize coordinate field names
    if 'lat' in site and 'latitude' not in site:
        standardized['latitude'] = site['lat']
    if 'lon' in site and 'longitude' not in site:
        standardized['longitude'] = site['lon']
    
    return standardized

def has_pycnopodia_mentions(site):
    """Check if site has Pycnopodia or sunflower star mentions."""
    search_fields = ["species", "description", "notes", "taxa", "location_description", 
                    "pycnopodia_evidence", "habitat_type", "sswd_status"]
    
    # Check boolean flag first
    if site.get("pycnopodia_documented", False):
        return True
    
    # Check text fields
    for field in search_fields:
        if field in site and site[field]:
            value = site[field]
            if isinstance(value, str):
                text = value.lower()
                if "pycnopodia" in text or "sunflower" in text:
                    return True
            elif isinstance(value, list):
                for item in value:
                    if isinstance(item, str):
                        text = item.lower()
                        if "pycnopodia" in text or "sunflower" in text:
                            return True
    
    return False

def main():
    # Find all site JSON files
    data_dir = Path("/home/starbot/.openclaw/workspace/sswd-evoepi/data/nodes")
    site_files = list(data_dir.glob("sites_*.json"))
    
    print(f"Found {len(site_files)} regional site files:")
    for f in sorted(site_files):
        print(f"  {f.name}")
    
    # Merge all sites
    all_sites = []
    region_stats = defaultdict(lambda: {"count": 0, "with_pycno": 0, "file_size": 0, "status": "OK"})
    
    for site_file in sorted(site_files):
        region = site_file.stem.replace("sites_", "")
        file_size = site_file.stat().st_size
        region_stats[region]["file_size"] = file_size
        
        try:
            with open(site_file, 'r') as f:
                data = json.load(f)
            
            sites = extract_sites_from_json(data, region)
            print(f"\n{region}: {len(sites)} sites ({file_size:,} bytes)")
            
            region_stats[region]["count"] = len(sites)
            
            # Process each site
            pycno_count = 0
            for site in sites:
                try:
                    standardized_site = standardize_site(site, region)
                    
                    if has_pycnopodia_mentions(standardized_site):
                        pycno_count += 1
                    
                    all_sites.append(standardized_site)
                    
                except Exception as e:
                    print(f"  WARNING: Error processing site in {region}: {e}")
                    continue
            
            region_stats[region]["with_pycno"] = pycno_count
            print(f"  → {pycno_count} sites with Pycnopodia documentation")
            
        except Exception as e:
            print(f"ERROR reading {site_file}: {e}")
            region_stats[region]["status"] = f"FAILED: {e}"
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
    successful_regions = 0
    for region in sorted(region_stats.keys()):
        stats = region_stats[region]
        status_marker = "✅" if stats["status"] == "OK" else "❌"
        print(f"{status_marker} {region:>6}: {stats['count']:>4} sites, {stats['with_pycno']:>3} with Pycnopodia ({stats['file_size']:>6,} bytes)")
        if stats["status"] == "OK":
            total_pycno += stats['with_pycno']
            successful_regions += 1
        else:
            print(f"         Status: {stats['status']}")
    
    failed_regions = len(region_stats) - successful_regions
    print("-" * 60)
    print(f"TOTAL:  {len(all_sites):>4} sites, {total_pycno:>3} with Pycnopodia ({sum(s['file_size'] for s in region_stats.values()):>6,} bytes)")
    print(f"Status: {successful_regions}/{len(region_stats)} regions successful, {failed_regions} failed")
    
    # Coordinate statistics
    if all_sites:
        lats = []
        lons = []
        for site in all_sites:
            # Check both possible coordinate field names
            lat = site.get('latitude') or site.get('lat')
            lon = site.get('longitude') or site.get('lon')
            if lat is not None and lon is not None:
                lats.append(float(lat))
                lons.append(float(lon))
        
        if lats and lons:
            print("\n" + "="*60)
            print("COORDINATE COVERAGE:")
            print("="*60)
            print(f"Sites with coordinates: {len(lats)}/{len(all_sites)} ({len(lats)/len(all_sites)*100:.1f}%)")
            print(f"Latitude:  {min(lats):.3f}° to {max(lats):.3f}° ({max(lats)-min(lats):.1f}° span)")
            print(f"Longitude: {min(lons):.3f}° to {max(lons):.3f}° ({max(lons)-min(lons):.1f}° span)")
        else:
            print("\n❌ No valid coordinates found in sites!")
    
    return output_file, len(all_sites), total_pycno, region_stats

if __name__ == "__main__":
    main()