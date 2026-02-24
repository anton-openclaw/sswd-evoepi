#!/usr/bin/env python3
"""
Interactive site editor for the Pycnopodia network.
Click map to add, click markers to edit, DRAG markers to reposition.
All changes logged to data/nodes/edit_history.jsonl for undo.

Usage: python3 scripts/site_editor.py [--port 5555]
"""

import json
import argparse
import webbrowser
from datetime import datetime, timezone
from pathlib import Path
from flask import Flask, request, jsonify, Response

BASE = Path(__file__).parent.parent
SITES_FILE = BASE / "data" / "nodes" / "all_sites.json"
HISTORY_FILE = BASE / "data" / "nodes" / "edit_history.jsonl"

app = Flask(__name__)


def load_sites():
    with open(SITES_FILE) as f:
        return json.load(f)


def save_sites(sites):
    with open(SITES_FILE, 'w') as f:
        json.dump(sites, f, indent=2)


def log_edit(action, index, before=None, after=None, note=""):
    """Append an edit record to the history log."""
    entry = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "action": action,
        "index": index,
        "note": note,
    }
    if before is not None:
        entry["before"] = before
    if after is not None:
        entry["after"] = after
    with open(HISTORY_FILE, 'a') as f:
        f.write(json.dumps(entry) + '\n')


REGION_COLORS = {
    'AK-WG': '#4363d8',  # royal blue
    'AK-AL': '#3cb44b',  # green
    'AK-EG': '#42d4f4',  # cyan
    'AK-SE': '#911eb4',  # purple
    'BC-N': '#f58231',   # orange
    'BC-C': '#e6194b',   # crimson
    'SS':   '#469990',   # teal
    'WA-O': '#f032e6',   # magenta
    'OR':   '#9A6324',   # brown
    'CA-N': '#bfef45',   # lime
    'CA-C': '#ffe119',   # yellow
    'CA-S': '#dcbeff',   # lavender
    'BJ':   '#800000',   # maroon
}

REGION_NAMES = {
    'AK-WG': 'Alaska - Western Gulf', 'AK-AL': 'Alaska - Aleutians',
    'AK-EG': 'Alaska - Eastern Gulf', 'AK-SE': 'Alaska - Southeast',
    'BC-N': 'British Columbia - North', 'BC-C': 'British Columbia - Central',
    'SS': 'Salish Sea', 'WA-O': 'Washington - Outer Coast',
    'OR': 'Oregon', 'CA-N': 'California - North',
    'CA-C': 'California - Central', 'CA-S': 'California - South',
    'BJ': 'Baja California',
}

HTML = r"""<!DOCTYPE html>
<html>
<head>
<title>Pycnopodia Site Editor</title>
<meta charset="utf-8" />
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css" />
<script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"></script>
<style>
  * { margin: 0; padding: 0; box-sizing: border-box; }
  body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif; }
  #map { position: absolute; top: 0; left: 0; right: 0; bottom: 0; z-index: 1; }

  .panel {
    position: fixed; top: 10px; right: 10px; z-index: 1000;
    background: white; border-radius: 10px; padding: 16px;
    box-shadow: 0 4px 12px rgba(0,0,0,0.15); width: 340px;
    max-height: calc(100vh - 20px); overflow-y: auto; font-size: 14px;
  }
  .panel h2 { font-size: 18px; margin-bottom: 12px; }
  .form-group { margin-bottom: 10px; }
  .form-group label { display: block; font-weight: 600; margin-bottom: 3px; color: #555; font-size: 12px; }
  .form-group input, .form-group select, .form-group textarea {
    width: 100%; padding: 7px 10px; border: 1px solid #ddd; border-radius: 6px;
    font-size: 13px; font-family: inherit;
  }
  .form-group textarea { resize: vertical; min-height: 50px; }
  .form-group input:focus, .form-group select:focus, .form-group textarea:focus {
    outline: none; border-color: #4a90d9; box-shadow: 0 0 0 2px rgba(74,144,217,0.2);
  }
  .coord-row { display: flex; gap: 8px; }
  .coord-row .form-group { flex: 1; }

  .btn {
    padding: 8px 16px; border: none; border-radius: 6px; cursor: pointer;
    font-size: 13px; font-weight: 600; margin-right: 6px; margin-top: 6px;
  }
  .btn-primary { background: #4a90d9; color: white; }
  .btn-primary:hover { background: #357abd; }
  .btn-danger { background: #e74c3c; color: white; }
  .btn-danger:hover { background: #c0392b; }
  .btn-secondary { background: #95a5a6; color: white; }
  .btn-secondary:hover { background: #7f8c8d; }
  .btn-undo { background: #f39c12; color: white; }
  .btn-undo:hover { background: #d68910; }

  .header {
    position: fixed; top: 10px; left: 60px; z-index: 1000;
    background: white; padding: 10px 16px; border-radius: 8px;
    box-shadow: 0 2px 6px rgba(0,0,0,0.2); border-left: 4px solid #3498db;
  }
  .header h3 { margin: 0; font-size: 16px; }
  .header p { margin: 4px 0 0; font-size: 12px; color: #666; }

  .stats {
    position: fixed; bottom: 30px; left: 10px; z-index: 1000;
    background: white; padding: 10px 14px; border-radius: 8px;
    box-shadow: 0 2px 6px rgba(0,0,0,0.2); font-size: 11px;
    max-height: 280px; overflow-y: auto; width: 220px;
  }
  .stats h4 { margin: 0 0 6px; }
  .stat-row { display: flex; justify-content: space-between; padding: 1px 0; }

  .history-panel {
    position: fixed; bottom: 30px; right: 10px; z-index: 1000;
    background: white; padding: 10px 14px; border-radius: 8px;
    box-shadow: 0 2px 6px rgba(0,0,0,0.2); font-size: 11px;
    max-height: 200px; overflow-y: auto; width: 300px; display: none;
  }
  .history-panel h4 { margin: 0 0 6px; display: flex; justify-content: space-between; }
  .hist-entry { padding: 3px 0; border-bottom: 1px solid #eee; }
  .hist-action { font-weight: bold; }
  .hist-add { color: #27ae60; }
  .hist-update { color: #2980b9; }
  .hist-move { color: #8e44ad; }
  .hist-delete { color: #e74c3c; }
  .hist-undo { color: #f39c12; }

  .toast {
    position: fixed; top: 60px; left: 50%; transform: translateX(-50%); z-index: 2000;
    background: #27ae60; color: white; padding: 10px 20px; border-radius: 8px;
    box-shadow: 0 4px 12px rgba(0,0,0,0.2); font-weight: 600;
    transition: opacity 0.3s; opacity: 0; pointer-events: none;
  }
  .toast.show { opacity: 1; }
  .toast.error { background: #e74c3c; }

  .coord-display {
    position: fixed; bottom: 10px; left: 50%; transform: translateX(-50%);
    z-index: 1000; background: rgba(0,0,0,0.7); color: white;
    padding: 4px 12px; border-radius: 4px; font-size: 12px; font-family: monospace;
  }

  .drag-hint {
    background: rgba(0,0,0,0.75); color: #fff; padding: 4px 8px;
    border-radius: 4px; font-size: 11px; white-space: nowrap;
  }
</style>
</head>
<body>

<div id="map"></div>

<div class="header">
  <h3>🌊 Pycnopodia Site Editor</h3>
  <p>Click to add (auto-names from coastline) · Drag to move · Click to edit · <a href="#" id="toggleHistory">History</a> · <kbd>Ctrl+Z</kbd> undo</p>
</div>

<div class="panel" id="panel" style="display:none;">
  <h2 id="panelTitle">Add New Site</h2>
  <form id="siteForm">
    <div class="form-group">
      <label>Name <span style="font-weight:normal;color:#999;">(auto-filled from location)</span></label>
      <input type="text" id="f_name" placeholder="Looking up location...">
    </div>
    <div class="coord-row">
      <div class="form-group">
        <label>Latitude</label>
        <input type="number" id="f_lat" step="0.0001" required>
      </div>
      <div class="form-group">
        <label>Longitude</label>
        <input type="number" id="f_lon" step="0.0001" required>
      </div>
    </div>
    <p style="font-size:11px;color:#888;margin:-6px 0 8px;">💡 Drag the marker on the map to reposition</p>
    <div class="form-group">
      <label>Region *</label>
      <select id="f_region"></select>
    </div>
    <div class="form-group">
      <label>Habitat</label>
      <input type="text" id="f_habitat" placeholder="e.g. rocky reef / kelp forest">
    </div>
    <div class="form-group">
      <label>Depth Range (m)</label>
      <input type="text" id="f_depth" placeholder="e.g. 5-20">
    </div>
    <div class="form-group">
      <label>SSWD Status</label>
      <input type="text" id="f_sswd" value="population severely impacted by SSWD starting 2013">
    </div>
    <div class="form-group">
      <label>Pycnopodia Documentation</label>
      <input type="text" id="f_pyc" placeholder="e.g. historically documented">
    </div>
    <div class="form-group">
      <label>Notes</label>
      <textarea id="f_notes" placeholder="Additional details..."></textarea>
    </div>
    <input type="hidden" id="f_editIdx" value="-1">
    <div>
      <button type="submit" class="btn btn-primary" id="btnSave">💾 Save</button>
      <button type="button" class="btn btn-danger" id="btnDelete" style="display:none;">🗑 Delete</button>
      <button type="button" class="btn btn-secondary" id="btnCancel">Cancel</button>
    </div>
  </form>
</div>

<div class="stats" id="stats"></div>

<div class="history-panel" id="historyPanel">
  <h4>
    <span>📝 Edit History</span>
    <button class="btn btn-undo" id="btnUndo" style="padding:4px 10px;margin:0;font-size:11px;">↩ Undo Last</button>
  </h4>
  <div id="historyList"></div>
</div>

<div class="coord-display" id="coords">Hover over map</div>
<div class="toast" id="toast"></div>

<script>
const REGION_COLORS = REGION_COLORS_JSON;
const REGION_NAMES = REGION_NAMES_JSON;

let sites = [];
let markers = [];
let map;
let activeMarker = null; // the draggable marker for editing/adding
let editingIdx = -1;
let panelOpen = false;

map = L.map('map', { zoomControl: true }).setView([48.0, -135.0], 4);

L.tileLayer('https://{s}.basemaps.cartocdn.com/light_all/{z}/{x}/{y}{r}.png', {
  attribution: '&copy; OpenStreetMap, &copy; CartoDB', maxZoom: 19
}).addTo(map);

const satellite = L.tileLayer('https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}', {
  attribution: 'Esri', maxZoom: 19
});
L.control.layers({
  'Map': L.tileLayer('https://{s}.basemaps.cartocdn.com/light_all/{z}/{x}/{y}{r}.png', {maxZoom: 19}),
  'Satellite': satellite
}, {}, { position: 'topleft' }).addTo(map);

map.on('mousemove', function(e) {
  document.getElementById('coords').textContent =
    e.latlng.lat.toFixed(4) + '°N, ' + Math.abs(e.latlng.lng).toFixed(4) + '°W';
});

// Reverse geocoding via Nominatim
let geocodeTimer = null;
async function reverseGeocode(lat, lng) {
  try {
    const resp = await fetch(
      `https://nominatim.openstreetmap.org/reverse?lat=${lat}&lon=${lng}&format=json&zoom=14&addressdetails=1&extratags=1`,
      { headers: { 'User-Agent': 'SSWD-EvoEpi-SiteEditor/1.0' } }
    );
    const data = await resp.json();
    return data;
  } catch (e) {
    console.warn('Geocode failed:', e);
    return null;
  }
}

function buildSiteName(geoData, lat, lng) {
  if (!geoData || geoData.error) return `Site at ${lat.toFixed(3)}°N, ${Math.abs(lng).toFixed(3)}°W`;
  const addr = geoData.address || {};
  const parts = [];
  // Prefer: natural feature > locality > town > county
  const feat = addr.natural || addr.beach || addr.bay || addr.cape || addr.peninsula
    || addr.water || addr.island || addr.archipelago || addr.reef;
  const locality = addr.hamlet || addr.village || addr.town || addr.city || addr.locality;
  const county = addr.county || addr.state_district;
  const state = addr.state;
  if (feat) parts.push(feat);
  if (locality && locality !== feat) parts.push(locality);
  if (!feat && !locality && geoData.name && geoData.name !== geoData.display_name) parts.push(geoData.name);
  if (parts.length === 0 && county) parts.push(county);
  if (state && parts.length < 2) parts.push(state);
  if (parts.length === 0) return `Site at ${lat.toFixed(3)}°N, ${Math.abs(lng).toFixed(3)}°W`;
  return parts.join(' - ');
}

async function autoFillFromLocation(lat, lng) {
  const nameField = document.getElementById('f_name');
  nameField.placeholder = 'Looking up location...';
  const geoData = await reverseGeocode(lat, lng);
  const name = buildSiteName(geoData, lat, lng);
  // Only auto-fill if user hasn't typed something
  if (!nameField.value || nameField.value.startsWith('Site at ') || nameField.value === nameField.dataset.autoName) {
    nameField.value = name;
    nameField.dataset.autoName = name;
  }
  nameField.placeholder = 'e.g. Tofino - Chesterman Beach';
}

// Region select
const regionSelect = document.getElementById('f_region');
Object.entries(REGION_NAMES).sort((a,b) => a[1].localeCompare(b[1])).forEach(([code, name]) => {
  const opt = document.createElement('option');
  opt.value = code; opt.textContent = name + ' (' + code + ')';
  regionSelect.appendChild(opt);
});

async function loadSites() {
  const resp = await fetch('/api/sites');
  sites = await resp.json();
  renderMarkers();
  renderStats();
}

function renderMarkers() {
  markers.forEach(m => map.removeLayer(m));
  markers = [];

  sites.forEach((site, idx) => {
    const lat = site.latitude || site.lat;
    const lon = site.longitude || site.lon;
    if (!lat || !lon) return;

    const region = site.region || '?';
    const color = REGION_COLORS[region] || '#333';
    const name = site.name || 'Site ' + idx;

    // All markers are draggable
    const marker = L.marker([lat, lon], {
      draggable: true,
      icon: idx === 283
        ? L.divIcon({ html: '<div style="font-size:20px;text-shadow:0 1px 3px rgba(0,0,0,0.5);">⭐</div>', iconSize: [24,24], iconAnchor: [12,12], className: '' })
        : L.divIcon({
            html: '<div style="width:12px;height:12px;border-radius:50%;background:'+color+';border:1.5px solid white;box-shadow:0 1px 3px rgba(0,0,0,0.3);"></div>',
            iconSize: [12,12], iconAnchor: [6,6], className: ''
          })
    });

    marker.bindTooltip(name + ' (#' + idx + ')', { direction: 'top', offset: [0, -8] });

    // Click to edit
    marker.on('click', function(e) {
      L.DomEvent.stopPropagation(e);
      editSite(idx);
    });

    // Drag to move
    marker.on('dragstart', function() {
      marker.closeTooltip();
    });
    marker.on('dragend', function(e) {
      const newPos = e.target.getLatLng();
      moveSite(idx, newPos.lat, newPos.lng);
    });

    marker.addTo(map);
    markers.push(marker);
  });
}

function renderStats() {
  const counts = {};
  sites.forEach(s => { const r = s.region || '?'; counts[r] = (counts[r]||0)+1; });
  const sorted = Object.entries(counts).sort((a,b) => b[1]-a[1]);
  let html = '<h4>' + sites.length + ' sites</h4>';
  sorted.forEach(([r, c]) => {
    html += '<div class="stat-row"><span style="color:'+(REGION_COLORS[r]||'#333')+'">●</span>'
         + '<span style="flex:1;margin:0 6px;">'+(REGION_NAMES[r]||r)+'</span>'
         + '<span style="font-weight:bold;">'+c+'</span></div>';
  });
  document.getElementById('stats').innerHTML = html;
}

// Move site via drag
async function moveSite(idx, newLat, newLng) {
  const resp = await fetch('/api/sites/' + idx + '/move', {
    method: 'POST',
    headers: {'Content-Type': 'application/json'},
    body: JSON.stringify({ latitude: newLat, longitude: newLng })
  });
  const result = await resp.json();
  if (result.ok) {
    showToast('Moved #' + idx + ' → ' + newLat.toFixed(4) + ', ' + newLng.toFixed(4));
    // Update local data without full reload
    sites[idx].latitude = newLat;
    sites[idx].longitude = newLng;
    if (sites[idx].lat !== undefined) sites[idx].lat = newLat;
    if (sites[idx].lon !== undefined) sites[idx].lon = newLng;
    // If edit panel is open for this site, update coords
    if (panelOpen && parseInt(document.getElementById('f_editIdx').value) === idx) {
      document.getElementById('f_lat').value = newLat.toFixed(4);
      document.getElementById('f_lon').value = newLng.toFixed(4);
    }
  }
}

// Click map to add
map.on('click', function(e) {
  if (panelOpen) return;
  showAddForm(e.latlng.lat, e.latlng.lng);
});

function showAddForm(lat, lng) {
  // Remove previous active marker
  if (activeMarker) { map.removeLayer(activeMarker); activeMarker = null; }

  activeMarker = L.marker([lat, lng], {
    draggable: true,
    icon: L.divIcon({
      html: '<div style="font-size:24px;filter:drop-shadow(0 2px 4px rgba(0,0,0,0.3));">📍</div>',
      iconSize: [24,24], iconAnchor: [12,24], className: ''
    })
  }).addTo(map);

  activeMarker.on('drag', function(e) {
    document.getElementById('f_lat').value = e.latlng.lat.toFixed(4);
    document.getElementById('f_lon').value = e.latlng.lng.toFixed(4);
    document.getElementById('f_region').value = guessRegion(e.latlng.lat, e.latlng.lng);
  });
  activeMarker.on('dragend', function(e) {
    // Re-geocode after drag
    autoFillFromLocation(e.target.getLatLng().lat, e.target.getLatLng().lng);
  });

  const guessed = guessRegion(lat, lng);
  document.getElementById('panelTitle').textContent = 'Add New Site';
  document.getElementById('f_name').value = '';
  document.getElementById('f_name').dataset.autoName = '';
  document.getElementById('f_lat').value = lat.toFixed(4);
  document.getElementById('f_lon').value = lng.toFixed(4);
  document.getElementById('f_region').value = guessed;
  document.getElementById('f_habitat').value = 'rocky reef / kelp forest';
  document.getElementById('f_depth').value = '5-25';
  document.getElementById('f_sswd').value = 'population severely impacted by SSWD starting 2013';
  document.getElementById('f_pyc').value = 'historically documented';
  document.getElementById('f_notes').value = '';
  document.getElementById('f_editIdx').value = '-1';
  document.getElementById('btnDelete').style.display = 'none';
  document.getElementById('btnSave').textContent = '💾 Add Site';
  document.getElementById('panel').style.display = 'block';
  panelOpen = true;

  // Auto-lookup name from location
  autoFillFromLocation(lat, lng);
}

function guessRegion(lat, lng) {
  if (lat > 58) return 'AK-EG';
  if (lat > 54 && lng < -160) return 'AK-AL';
  if (lat > 54 && lng < -140) return 'AK-SE';
  if (lat > 54) return 'BC-N';
  if (lat > 51) return 'BC-C';
  if (lat > 48.5 && lng > -125) return 'SS';
  if (lat > 48.5) return 'BC-C';
  if (lat > 46) return 'WA-O';
  if (lat > 42) return 'OR';
  if (lat > 39) return 'CA-N';
  if (lat > 35) return 'CA-C';
  if (lat > 30) return 'CA-S';
  return 'BJ';
}

function editSite(idx) {
  const site = sites[idx];
  if (activeMarker) { map.removeLayer(activeMarker); activeMarker = null; }

  document.getElementById('panelTitle').textContent = 'Edit Site #' + idx;
  document.getElementById('f_name').value = site.name || '';
  document.getElementById('f_lat').value = (site.latitude || site.lat || '');
  document.getElementById('f_lon').value = (site.longitude || site.lon || '');
  document.getElementById('f_region').value = site.region || '';
  document.getElementById('f_habitat').value = site.habitat || '';
  document.getElementById('f_depth').value = site.depth_range_m || '';
  document.getElementById('f_sswd').value = site.sswd_status || '';
  document.getElementById('f_pyc').value = site.pycnopodia_documentation || site.pycnopodia_documented || '';
  document.getElementById('f_notes').value = site.notes || '';
  document.getElementById('f_editIdx').value = idx;
  document.getElementById('btnDelete').style.display = 'inline-block';
  document.getElementById('btnSave').textContent = '💾 Update';
  document.getElementById('panel').style.display = 'block';
  panelOpen = true;

  const lat = site.latitude || site.lat;
  const lon = site.longitude || site.lon;
  if (lat && lon) map.panTo([lat, lon]);
}

// Save form
document.getElementById('siteForm').addEventListener('submit', async function(e) {
  e.preventDefault();
  const idx = parseInt(document.getElementById('f_editIdx').value);
  const lat = parseFloat(document.getElementById('f_lat').value);
  const lon = parseFloat(document.getElementById('f_lon').value);
  const nameVal = document.getElementById('f_name').value || 'Site at ' + lat.toFixed(3) + '°N, ' + Math.abs(lon).toFixed(3) + '°W';
  const siteData = {
    name: nameVal,
    latitude: lat,
    longitude: lon,
    region: document.getElementById('f_region').value,
    habitat: document.getElementById('f_habitat').value,
    depth_range_m: document.getElementById('f_depth').value,
    sswd_status: document.getElementById('f_sswd').value,
    pycnopodia_documented: true,
    pycnopodia_documentation: document.getElementById('f_pyc').value,
    notes: document.getElementById('f_notes').value,
  };

  let resp;
  if (idx >= 0) {
    resp = await fetch('/api/sites/' + idx, {
      method: 'PUT', headers: {'Content-Type': 'application/json'},
      body: JSON.stringify(siteData)
    });
  } else {
    resp = await fetch('/api/sites', {
      method: 'POST', headers: {'Content-Type': 'application/json'},
      body: JSON.stringify(siteData)
    });
  }

  const result = await resp.json();
  if (result.ok) {
    showToast(idx >= 0 ? 'Updated #' + idx : 'Added #' + result.index + ' ✓');
    closePanel();
    loadSites();
    loadHistory();
  } else {
    showToast('Error: ' + result.error, true);
  }
});

// Delete
document.getElementById('btnDelete').addEventListener('click', async function() {
  const idx = parseInt(document.getElementById('f_editIdx').value);
  if (idx < 0) return;
  if (!confirm('Delete #' + idx + ' (' + sites[idx].name + ')?\n\nThis will be logged in history and can be undone.')) return;

  const resp = await fetch('/api/sites/' + idx, { method: 'DELETE' });
  const result = await resp.json();
  if (result.ok) {
    showToast('Deleted #' + idx);
    closePanel();
    loadSites();
    loadHistory();
  }
});

document.getElementById('btnCancel').addEventListener('click', closePanel);

function closePanel() {
  document.getElementById('panel').style.display = 'none';
  panelOpen = false;
  if (activeMarker) { map.removeLayer(activeMarker); activeMarker = null; }
}

// History
document.getElementById('toggleHistory').addEventListener('click', function(e) {
  e.preventDefault();
  const hp = document.getElementById('historyPanel');
  hp.style.display = hp.style.display === 'none' ? 'block' : 'none';
  if (hp.style.display === 'block') loadHistory();
});

async function loadHistory() {
  const resp = await fetch('/api/history');
  const entries = await resp.json();
  const list = document.getElementById('historyList');
  if (entries.length === 0) {
    list.innerHTML = '<p style="color:#999;">No edits yet</p>';
    return;
  }
  // Show newest first, limit 30
  const recent = entries.slice(-30).reverse();
  list.innerHTML = recent.map(e => {
    const cls = 'hist-' + (e.action || 'update');
    const time = e.timestamp ? new Date(e.timestamp).toLocaleTimeString() : '?';
    const name = (e.after && e.after.name) || (e.before && e.before.name) || '#' + e.index;
    let desc = '';
    if (e.action === 'add') desc = 'Added ' + name;
    else if (e.action === 'delete') desc = 'Deleted ' + name;
    else if (e.action === 'move') desc = 'Moved ' + name;
    else if (e.action === 'update') desc = 'Updated ' + name;
    else if (e.action === 'undo') desc = 'Undid: ' + (e.note || '?');
    else desc = e.action + ' ' + name;
    return '<div class="hist-entry"><span class="hist-action '+cls+'">'+e.action+'</span> ' + desc + ' <span style="color:#aaa;float:right;">'+time+'</span></div>';
  }).join('');
}

// Undo
document.getElementById('btnUndo').addEventListener('click', async function() {
  const resp = await fetch('/api/undo', { method: 'POST' });
  const result = await resp.json();
  if (result.ok) {
    showToast('↩ ' + result.message);
    loadSites();
    loadHistory();
  } else {
    showToast(result.error || 'Nothing to undo', true);
  }
});

function showToast(msg, isError) {
  const t = document.getElementById('toast');
  t.textContent = msg;
  t.className = 'toast show' + (isError ? ' error' : '');
  setTimeout(() => t.className = 'toast', 2500);
}

document.addEventListener('keydown', function(e) {
  if (e.key === 'Escape') closePanel();
  if ((e.ctrlKey || e.metaKey) && e.key === 'z') {
    e.preventDefault();
    document.getElementById('btnUndo').click();
  }
});

loadSites();
</script>
</body>
</html>"""


@app.route('/')
def index():
    html = HTML.replace('REGION_COLORS_JSON', json.dumps(REGION_COLORS))
    html = html.replace('REGION_NAMES_JSON', json.dumps(REGION_NAMES))
    return Response(html, mimetype='text/html')


@app.route('/api/sites', methods=['GET'])
def get_sites():
    return jsonify(load_sites())


@app.route('/api/sites', methods=['POST'])
def add_site():
    sites = load_sites()
    site = request.json
    sites.append(site)
    save_sites(sites)
    log_edit('add', len(sites) - 1, after=site)
    return jsonify({'ok': True, 'index': len(sites) - 1, 'total': len(sites)})


@app.route('/api/sites/<int:idx>', methods=['PUT'])
def update_site(idx):
    sites = load_sites()
    if idx < 0 or idx >= len(sites):
        return jsonify({'ok': False, 'error': 'Invalid index'}), 404
    before = dict(sites[idx])
    updated = request.json
    sites[idx].update(updated)
    save_sites(sites)
    log_edit('update', idx, before=before, after=dict(sites[idx]))
    return jsonify({'ok': True, 'index': idx})


@app.route('/api/sites/<int:idx>/move', methods=['POST'])
def move_site(idx):
    """Drag-to-move endpoint — just updates lat/lon."""
    sites = load_sites()
    if idx < 0 or idx >= len(sites):
        return jsonify({'ok': False, 'error': 'Invalid index'}), 404
    data = request.json
    before = dict(sites[idx])
    sites[idx]['latitude'] = data['latitude']
    sites[idx]['longitude'] = data['longitude']
    if 'lat' in sites[idx]:
        sites[idx]['lat'] = data['latitude']
    if 'lon' in sites[idx]:
        sites[idx]['lon'] = data['longitude']
    save_sites(sites)
    log_edit('move', idx, before=before, after=dict(sites[idx]),
             note=f"{before.get('latitude','?')},{before.get('longitude','?')} → {data['latitude']:.4f},{data['longitude']:.4f}")
    return jsonify({'ok': True})


@app.route('/api/sites/<int:idx>', methods=['DELETE'])
def delete_site(idx):
    sites = load_sites()
    if idx < 0 or idx >= len(sites):
        return jsonify({'ok': False, 'error': 'Invalid index'}), 404
    removed = sites.pop(idx)
    save_sites(sites)
    log_edit('delete', idx, before=removed)
    return jsonify({'ok': True, 'removed': removed.get('name', '?'), 'total': len(sites)})


@app.route('/api/history', methods=['GET'])
def get_history():
    if not HISTORY_FILE.exists():
        return jsonify([])
    entries = []
    with open(HISTORY_FILE) as f:
        for line in f:
            line = line.strip()
            if line:
                try:
                    entries.append(json.loads(line))
                except json.JSONDecodeError:
                    pass
    return jsonify(entries)


@app.route('/api/undo', methods=['POST'])
def undo_last():
    """Undo the most recent edit."""
    if not HISTORY_FILE.exists():
        return jsonify({'ok': False, 'error': 'No history'})

    entries = []
    with open(HISTORY_FILE) as f:
        for line in f:
            line = line.strip()
            if line:
                try:
                    entries.append(json.loads(line))
                except json.JSONDecodeError:
                    pass

    if not entries:
        return jsonify({'ok': False, 'error': 'No history'})

    last = entries[-1]
    sites = load_sites()
    action = last.get('action')
    idx = last.get('index', -1)
    msg = ''

    if action == 'add':
        # Undo add = delete the last-added site
        if idx < len(sites):
            removed = sites.pop(idx)
            msg = f"Removed added site '{removed.get('name', '?')}'"
        else:
            return jsonify({'ok': False, 'error': 'Site index out of range'})

    elif action == 'delete':
        # Undo delete = re-insert the site
        before = last.get('before')
        if before:
            sites.insert(idx, before)
            msg = f"Restored deleted site '{before.get('name', '?')}'"
        else:
            return jsonify({'ok': False, 'error': 'No before state in history'})

    elif action in ('update', 'move'):
        # Undo update/move = restore before state
        before = last.get('before')
        if before and idx < len(sites):
            sites[idx] = before
            msg = f"Restored site #{idx} to previous state"
        else:
            return jsonify({'ok': False, 'error': 'Cannot restore'})

    else:
        return jsonify({'ok': False, 'error': f'Cannot undo action: {action}'})

    save_sites(sites)
    log_edit('undo', idx, note=f"Undid {action} on #{idx}")

    return jsonify({'ok': True, 'message': msg})


def main():
    parser = argparse.ArgumentParser(description='Interactive Pycnopodia site editor')
    parser.add_argument('--port', '-p', type=int, default=5555)
    parser.add_argument('--no-open', action='store_true')
    args = parser.parse_args()

    print(f"🌊 Pycnopodia Site Editor")
    print(f"   {len(load_sites())} sites loaded")
    print(f"   http://localhost:{args.port}")
    print(f"   Click to add · Click marker to edit · Drag to reposition")
    print(f"   Ctrl+Z to undo · All changes logged to edit_history.jsonl")
    print()

    if not args.no_open:
        import threading
        threading.Timer(1.0, lambda: webbrowser.open(f'http://localhost:{args.port}')).start()

    app.run(host='127.0.0.1', port=args.port, debug=False)


if __name__ == '__main__':
    main()
