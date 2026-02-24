#!/usr/bin/env python3
"""
Interactive map of all 489 Pycnopodia survey/monitoring sites.
Opens in your default browser. Pan, zoom, click markers for details.

Usage: python3 scripts/site_map.py [--output path/to/map.html] [--open]
"""

import json
import argparse
import webbrowser
from pathlib import Path

import folium
from folium.plugins import MarkerCluster

BASE = Path(__file__).parent.parent
SITES_FILE = BASE / "data" / "nodes" / "all_sites.json"
DEFAULT_OUTPUT = BASE / "results" / "site_map_interactive.html"

# Region colors
REGION_COLORS = {
    'AK-WG': '#1f77b4',   # Western Gulf AK - dark blue
    'AK-AL': '#2ca02c',   # Aleutians - green
    'AK-EG': '#17becf',   # Eastern Gulf AK - cyan
    'AK-SE': '#9467bd',   # SE Alaska - purple
    'BC-N': '#ff7f0e',    # Northern BC - orange
    'BC-C': '#d62728',    # Central BC - red
    'SS': '#8c564b',      # Salish Sea - brown
    'WA-O': '#e377c2',    # WA outer coast - pink
    'OR': '#7f7f7f',      # Oregon - gray
    'CA-N': '#bcbd22',    # Northern CA - olive
    'CA-C': '#ff9896',    # Central CA - salmon
    'CA-S': '#c5b0d5',    # Southern CA - lavender
    'BJ': '#c49c94',      # Baja - tan
}

REGION_NAMES = {
    'AK-WG': 'Alaska - Western Gulf',
    'AK-AL': 'Alaska - Aleutians',
    'AK-EG': 'Alaska - Eastern Gulf',
    'AK-SE': 'Alaska - Southeast',
    'BC-N': 'British Columbia - North',
    'BC-C': 'British Columbia - Central',
    'SS': 'Salish Sea',
    'WA-O': 'Washington - Outer Coast',
    'OR': 'Oregon',
    'CA-N': 'California - North',
    'CA-C': 'California - Central',
    'CA-S': 'California - South',
    'BJ': 'Baja California',
}


def build_popup(site, idx):
    """Build HTML popup for a site marker."""
    name = site.get('name', 'Unknown')
    region = site.get('region', '?')
    region_full = REGION_NAMES.get(region, region)
    lat = site.get('latitude', site.get('lat', '?'))
    lon = site.get('longitude', site.get('lon', '?'))
    habitat = site.get('habitat', 'Unknown')
    sswd = site.get('sswd_status', 'Unknown')
    depth = site.get('depth_range_m', '?')
    notes = site.get('notes', '')
    pyc_doc = site.get('pycnopodia_documentation', site.get('pycnopodia_documented', ''))
    site_id = site.get('id', '')
    
    # Monitoring programs
    programs = site.get('monitoring_programs', [])
    prog_html = ', '.join(programs) if programs else 'None listed'
    
    html = f"""
    <div style="font-family: -apple-system, sans-serif; max-width: 320px;">
        <h3 style="margin: 0 0 8px 0; color: #1a1a2e;">{name}</h3>
        <table style="font-size: 13px; border-collapse: collapse; width: 100%;">
            <tr><td style="padding: 3px 8px 3px 0; font-weight: bold; color: #555;">Node Index</td>
                <td style="padding: 3px 0;">{idx}</td></tr>
            <tr><td style="padding: 3px 8px 3px 0; font-weight: bold; color: #555;">ID</td>
                <td style="padding: 3px 0;">{site_id}</td></tr>
            <tr><td style="padding: 3px 8px 3px 0; font-weight: bold; color: #555;">Region</td>
                <td style="padding: 3px 0;">{region_full} ({region})</td></tr>
            <tr><td style="padding: 3px 8px 3px 0; font-weight: bold; color: #555;">Coordinates</td>
                <td style="padding: 3px 0;">{lat:.4f}°N, {abs(lon) if isinstance(lon, (int,float)) else lon}°W</td></tr>
            <tr><td style="padding: 3px 8px 3px 0; font-weight: bold; color: #555;">Habitat</td>
                <td style="padding: 3px 0;">{habitat}</td></tr>
            <tr><td style="padding: 3px 8px 3px 0; font-weight: bold; color: #555;">Depth</td>
                <td style="padding: 3px 0;">{depth} m</td></tr>
            <tr><td style="padding: 3px 8px 3px 0; font-weight: bold; color: #555;">SSWD Status</td>
                <td style="padding: 3px 0;">{sswd}</td></tr>
            <tr><td style="padding: 3px 8px 3px 0; font-weight: bold; color: #555;">Pycnopodia</td>
                <td style="padding: 3px 0;">{pyc_doc}</td></tr>
            <tr><td style="padding: 3px 8px 3px 0; font-weight: bold; color: #555;">Monitoring</td>
                <td style="padding: 3px 0;">{prog_html}</td></tr>
        </table>
        {"<p style='font-size: 12px; color: #666; margin: 6px 0 0;'>" + notes + "</p>" if notes else ""}
    </div>
    """
    return folium.Popup(html, max_width=350)


def build_map(sites):
    """Build the interactive Folium map."""
    # Center on the Pacific NW coast
    m = folium.Map(
        location=[48.0, -135.0],
        zoom_start=4,
        tiles='CartoDB positron',
        control_scale=True,
    )
    
    # Add tile layers
    folium.TileLayer('OpenStreetMap', name='OpenStreetMap').add_to(m)
    folium.TileLayer(
        tiles='https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}',
        attr='Esri',
        name='Satellite',
    ).add_to(m)
    
    # Create feature groups per region (toggleable in layer control)
    region_groups = {}
    for region_code in REGION_COLORS:
        name = REGION_NAMES.get(region_code, region_code)
        fg = folium.FeatureGroup(name=f"{name} ({region_code})", show=True)
        region_groups[region_code] = fg
    
    # Special markers
    monterey_outplant = None
    
    for idx, site in enumerate(sites):
        lat = site.get('latitude', site.get('lat'))
        lon = site.get('longitude', site.get('lon'))
        if lat is None or lon is None:
            continue
            
        region = site.get('region', 'Unknown')
        color = REGION_COLORS.get(region, '#333333')
        name = site.get('name', f'Site {idx}')
        
        popup = build_popup(site, idx)
        tooltip = f"{name} (#{idx}, {region})"
        
        # Special icon for the Monterey outplanting site
        if idx == 283:
            marker = folium.Marker(
                location=[lat, lon],
                popup=popup,
                tooltip=tooltip,
                icon=folium.Icon(color='red', icon='star', prefix='fa'),
            )
            monterey_outplant = (lat, lon)
        else:
            marker = folium.CircleMarker(
                location=[lat, lon],
                radius=6,
                popup=popup,
                tooltip=tooltip,
                color=color,
                fill=True,
                fill_color=color,
                fill_opacity=0.7,
                weight=1,
            )
        
        fg = region_groups.get(region)
        if fg:
            marker.add_to(fg)
        else:
            marker.add_to(m)
    
    # Add all region groups to map
    for fg in region_groups.values():
        fg.add_to(m)
    
    # Layer control
    folium.LayerControl(collapsed=False).add_to(m)
    
    # Title
    title_html = """
    <div style="position: fixed; top: 10px; left: 60px; z-index: 1000;
                background: white; padding: 10px 16px; border-radius: 8px;
                box-shadow: 0 2px 6px rgba(0,0,0,0.2); font-family: -apple-system, sans-serif;">
        <h3 style="margin: 0; font-size: 16px;">🌊 Pycnopodia Network — 489 Sites</h3>
        <p style="margin: 4px 0 0; font-size: 12px; color: #666;">
            ⭐ = Monterey Outplanting Site (node 283) · Click markers for details · Toggle regions in layer panel
        </p>
    </div>
    """
    m.get_root().html.add_child(folium.Element(title_html))
    
    # Stats box
    from collections import Counter
    region_counts = Counter(s.get('region', '?') for s in sites)
    stats_rows = ''.join(
        f'<div style="display:flex;justify-content:space-between;padding:1px 0;">'
        f'<span style="color:{REGION_COLORS.get(r, "#333")};font-weight:bold;">●</span>'
        f'<span style="flex:1;margin:0 6px;">{REGION_NAMES.get(r,r)}</span>'
        f'<span style="font-weight:bold;">{c}</span></div>'
        for r, c in sorted(region_counts.items(), key=lambda x: -x[1])
    )
    stats_html = f"""
    <div style="position: fixed; bottom: 30px; left: 10px; z-index: 1000;
                background: white; padding: 10px 14px; border-radius: 8px;
                box-shadow: 0 2px 6px rgba(0,0,0,0.2); font-family: -apple-system, sans-serif;
                font-size: 11px; max-height: 300px; overflow-y: auto; width: 240px;">
        <h4 style="margin: 0 0 6px;">Regions ({len(sites)} sites)</h4>
        {stats_rows}
    </div>
    """
    m.get_root().html.add_child(folium.Element(stats_html))
    
    return m


def main():
    parser = argparse.ArgumentParser(description='Interactive map of Pycnopodia sites')
    parser.add_argument('--output', '-o', type=str, default=str(DEFAULT_OUTPUT),
                        help='Output HTML file path')
    parser.add_argument('--open', action='store_true', default=True,
                        help='Open in browser after generating (default: True)')
    parser.add_argument('--no-open', action='store_true',
                        help='Do not open in browser')
    args = parser.parse_args()
    
    print(f"Loading sites from {SITES_FILE}...")
    with open(SITES_FILE) as f:
        sites = json.load(f)
    print(f"Loaded {len(sites)} sites")
    
    print("Building map...")
    m = build_map(sites)
    
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    m.save(str(output))
    print(f"✓ Saved to {output}")
    print(f"  File size: {output.stat().st_size / 1024:.0f} KB")
    
    if args.open and not args.no_open:
        print("Opening in browser...")
        webbrowser.open(f'file://{output.resolve()}')


if __name__ == '__main__':
    main()
