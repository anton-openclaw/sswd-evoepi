#!/usr/bin/env python3
"""Generate a standalone HTML5 Canvas player for wildfire visualization.

Produces a single self-contained HTML file with embedded snapshot data,
play/pause/step/scrub controls, and the same dark marine theme as the
matplotlib version.

Usage:
    from sswd_evoepi.snapshots import SnapshotRecorder
    from viz.wildfire_player import generate_player

    snap = SnapshotRecorder.load("snapshots.npz")
    generate_player(snap, "wildfire_player.html",
                    node_names=["Sitka", "Howe Sound", ...],
                    habitat_sides={0: 223.6, 1: 141.4, ...})
"""

from __future__ import annotations

import json
import base64
import struct
from pathlib import Path
from typing import List, Optional, Dict

import numpy as np

from sswd_evoepi.snapshots import SnapshotRecorder
from sswd_evoepi.types import DiseaseState


def generate_player(
    recorder: SnapshotRecorder,
    output_path: str = "wildfire_player.html",
    node_names: Optional[List[str]] = None,
    habitat_sides: Optional[Dict[int, float]] = None,
    max_frames: int = 600,
    title: str = "SSWD-EvoEpi Wildfire Viewer",
):
    """Generate standalone HTML5 wildfire player.

    Args:
        recorder: SnapshotRecorder with data.
        output_path: Output HTML file path.
        node_names: List of node names indexed by node_id.
        habitat_sides: Dict node_id -> habitat side length (m).
        max_frames: Max unique days (subsamples if exceeded).
        title: Page title.
    """
    days = recorder.get_days()
    nodes = recorder.get_nodes()

    if not days or not nodes:
        print("No data to export.")
        return

    # Subsample days if needed
    if len(days) > max_frames:
        step = max(1, len(days) // max_frames)
        days = days[::step]

    if habitat_sides is None:
        habitat_sides = {}
    if node_names is None:
        node_names = [f"Node {n}" for n in nodes]

    # Pack snapshot data as JSON-friendly structure
    # For each day: { day, year, doy, nodes: { nid: { x[], y[], ds[], sz[], r[], n } } }
    frames = []
    for sim_day in days:
        year = sim_day // 365
        doy = sim_day % 365
        frame = {"day": sim_day, "year": year, "doy": doy, "nodes": {}}
        for nid in nodes:
            snap = recorder.get_snapshot(sim_day, nid)
            if snap is None or snap.n_alive == 0:
                frame["nodes"][str(nid)] = {"n": 0}
            else:
                frame["nodes"][str(nid)] = {
                    "x": snap.x.tolist(),
                    "y": snap.y.tolist(),
                    "ds": snap.disease_state.astype(int).tolist(),
                    "sz": np.round(snap.size, 1).tolist(),
                    "r": np.round(snap.resistance, 3).tolist(),
                    "n": snap.n_alive,
                }
        frames.append(frame)

    metadata = {
        "nodes": [{"id": nid, "name": node_names[nid] if nid < len(node_names) else f"Node {nid}",
                    "side": habitat_sides.get(nid, 100.0)} for nid in nodes],
        "n_frames": len(frames),
        "n_days": len(days),
    }

    # Estimate size
    json_data = json.dumps({"meta": metadata, "frames": frames}, separators=(',', ':'))
    size_mb = len(json_data) / 1024 / 1024
    print(f"Data size: {size_mb:.1f} MB ({len(frames)} frames, {len(nodes)} nodes)")

    html = _build_html(json_data, metadata, title)

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        f.write(html)

    file_size = Path(output_path).stat().st_size
    print(f"Player saved: {output_path} ({file_size / 1024:.0f} KB)")


def _build_html(json_data: str, metadata: dict, title: str) -> str:
    """Build the complete HTML document."""
    n_nodes = len(metadata["nodes"])

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>{title}</title>
<style>
* {{ margin: 0; padding: 0; box-sizing: border-box; }}
body {{
    background: #0d1b2a;
    color: #c8d6e5;
    font-family: -apple-system, 'Segoe UI', Roboto, Helvetica, sans-serif;
    display: flex;
    flex-direction: column;
    align-items: center;
    min-height: 100vh;
    padding: 12px;
}}
h1 {{
    font-size: 18px;
    font-weight: 600;
    margin: 8px 0 4px;
    color: #e0e8f0;
}}
#time-label {{
    font-size: 15px;
    color: #8899aa;
    margin-bottom: 8px;
    font-variant-numeric: tabular-nums;
}}
#canvas-container {{
    display: flex;
    flex-wrap: wrap;
    justify-content: center;
    gap: 6px;
    margin-bottom: 10px;
}}
canvas.node-canvas {{
    border: 1.5px solid #3d5a80;
    border-radius: 4px;
    background: #1b2838;
}}
.node-label {{
    text-align: center;
    font-size: 11px;
    font-weight: 600;
    color: #c8d6e5;
    margin-bottom: 2px;
}}
.node-stats {{
    text-align: center;
    font-size: 9px;
    color: #667788;
    font-family: 'SF Mono', 'Consolas', monospace;
    margin-top: 1px;
    height: 14px;
}}
.node-wrap {{ display: flex; flex-direction: column; align-items: center; }}
#controls {{
    display: flex;
    align-items: center;
    gap: 10px;
    margin: 6px 0;
    flex-wrap: wrap;
    justify-content: center;
}}
button {{
    background: #1b2838;
    color: #c8d6e5;
    border: 1px solid #3d5a80;
    border-radius: 4px;
    padding: 6px 14px;
    font-size: 13px;
    cursor: pointer;
    transition: background 0.15s;
}}
button:hover {{ background: #2d3f52; }}
button:active {{ background: #3d5a80; }}
button.active {{ background: #3d5a80; border-color: #66CCEE; }}
#scrubber {{
    width: 400px;
    max-width: 80vw;
    accent-color: #66CCEE;
}}
#speed-label {{ font-size: 11px; color: #667788; }}
#legend {{
    display: flex;
    gap: 14px;
    flex-wrap: wrap;
    justify-content: center;
    margin: 6px 0;
    font-size: 11px;
}}
.legend-item {{
    display: flex;
    align-items: center;
    gap: 4px;
}}
.legend-dot {{
    width: 10px;
    height: 10px;
    border-radius: 50%;
    display: inline-block;
}}
.legend-x {{
    font-size: 14px;
    font-weight: bold;
    line-height: 10px;
}}
#frame-counter {{
    font-size: 11px;
    color: #556677;
    font-variant-numeric: tabular-nums;
}}
</style>
</head>
<body>
<h1>{title}</h1>
<div id="time-label">Loading...</div>

<div id="canvas-container"></div>

<div id="controls">
    <button id="btn-prev" title="Previous frame">⏮</button>
    <button id="btn-play" title="Play/Pause" class="active">▶</button>
    <button id="btn-next" title="Next frame">⏭</button>
    <input type="range" id="scrubber" min="0" max="0" value="0">
    <button id="btn-slower" title="Slower">−</button>
    <span id="speed-label">1×</span>
    <button id="btn-faster" title="Faster">+</button>
    <span id="frame-counter"></span>
</div>

<div id="legend">
    <span class="legend-item"><span class="legend-dot" style="background:#4477AA"></span> Susceptible</span>
    <span class="legend-item"><span class="legend-dot" style="background:#CCBB44"></span> Exposed</span>
    <span class="legend-item"><span class="legend-dot" style="background:#EE7733"></span> Pre-symptomatic</span>
    <span class="legend-item"><span class="legend-dot" style="background:#EE6677"></span> Wasting</span>
    <span class="legend-item"><span class="legend-x" style="color:#888888">×</span> Dead</span>
    <span class="legend-item"><span class="legend-dot" style="background:#66CCEE"></span> Recovered</span>
</div>

<script>
const DATA = {json_data};
const META = DATA.meta;
const FRAMES = DATA.frames;

const DS_COLORS = {{
    0: '#4477AA',  // S
    1: '#CCBB44',  // E
    2: '#EE7733',  // I1
    3: '#EE6677',  // I2
    4: '#888888',  // D
    5: '#66CCEE',  // R
}};
const DS_NAMES = {{0:'S', 1:'E', 2:'I₁', 3:'I₂', 4:'D', 5:'R'}};
const MONTHS = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'];

// Canvas size per node
const CW = Math.min(280, Math.floor((window.innerWidth - 60) / Math.min(META.nodes.length, 5)));
const CH = CW;

// Create canvases
const container = document.getElementById('canvas-container');
const canvases = [];
const ctxs = [];
const statsEls = [];

META.nodes.forEach((node, i) => {{
    const wrap = document.createElement('div');
    wrap.className = 'node-wrap';
    const label = document.createElement('div');
    label.className = 'node-label';
    label.textContent = node.name;
    const cvs = document.createElement('canvas');
    cvs.className = 'node-canvas';
    cvs.width = CW * 2;  // retina
    cvs.height = CH * 2;
    cvs.style.width = CW + 'px';
    cvs.style.height = CH + 'px';
    const stats = document.createElement('div');
    stats.className = 'node-stats';
    wrap.appendChild(label);
    wrap.appendChild(cvs);
    wrap.appendChild(stats);
    container.appendChild(wrap);
    canvases.push(cvs);
    ctxs.push(cvs.getContext('2d'));
    statsEls.push(stats);
}});

// State
let frameIdx = 0;
let playing = false;
let speed = 1;
let timer = null;
const speeds = [0.25, 0.5, 1, 2, 4, 8];
let speedIdx = 2;

const scrubber = document.getElementById('scrubber');
scrubber.max = FRAMES.length - 1;

function renderFrame(idx) {{
    const frame = FRAMES[idx];
    const timeLabel = document.getElementById('time-label');
    const mo = MONTHS[Math.min(Math.floor(frame.doy / 30), 11)];
    timeLabel.textContent = `Year ${{frame.year}}  ·  ${{mo}} (day ${{frame.doy}})`;

    META.nodes.forEach((node, i) => {{
        const ctx = ctxs[i];
        const side = node.side;
        const scale = (CW * 2) / (side + 4);  // retina
        const offset = 2 * scale;

        // Clear
        ctx.fillStyle = '#1b2838';
        ctx.fillRect(0, 0, CW * 2, CH * 2);

        const nd = frame.nodes[String(node.id)];
        if (!nd || nd.n === 0) {{
            ctx.fillStyle = '#555555';
            ctx.font = 'italic 28px sans-serif';
            ctx.textAlign = 'center';
            ctx.fillText('Extinct', CW, CH);
            statsEls[i].textContent = 'N=0';
            return;
        }}

        const {{ x, y, ds, sz, r, n }} = nd;

        // Draw scale bar
        const barLen = Math.min(50, side * 0.25);
        const bx = 6;
        const by = CH * 2 - 10;
        ctx.strokeStyle = '#556677';
        ctx.lineWidth = 2;
        ctx.beginPath();
        ctx.moveTo(bx, by);
        ctx.lineTo(bx + barLen * scale, by);
        ctx.stroke();
        ctx.fillStyle = '#556677';
        ctx.font = '10px sans-serif';
        ctx.textAlign = 'center';
        ctx.fillText(barLen + 'm', bx + barLen * scale / 2, by - 4);

        // Draw agents
        for (let j = 0; j < x.length; j++) {{
            const px = offset + x[j] * scale;
            const py = (CH * 2) - offset - y[j] * scale;
            const state = ds[j];
            const color = DS_COLORS[state] || '#888888';
            const radius = Math.max(2, 1.5 + 5 * (sz[j] / 1000));
            const alpha = state === 4 ? 0.3 : (0.5 + 0.5 * Math.min(r[j], 1));

            if (state === 4) {{
                // Dead: × marker
                ctx.strokeStyle = color;
                ctx.globalAlpha = alpha;
                ctx.lineWidth = 1.2;
                ctx.beginPath();
                ctx.moveTo(px - 3, py - 3);
                ctx.lineTo(px + 3, py + 3);
                ctx.moveTo(px + 3, py - 3);
                ctx.lineTo(px - 3, py + 3);
                ctx.stroke();
                ctx.globalAlpha = 1;
            }} else {{
                // Alive: filled circle
                ctx.globalAlpha = alpha;
                ctx.fillStyle = color;
                ctx.beginPath();
                ctx.arc(px, py, radius, 0, Math.PI * 2);
                ctx.fill();
                ctx.globalAlpha = 1;
            }}
        }}

        // Stats
        const counts = [0,0,0,0,0,0];
        for (let j = 0; j < ds.length; j++) counts[ds[j]]++;
        statsEls[i].textContent = `N=${{n}}  S=${{counts[0]}}  E=${{counts[1]}}  I₁=${{counts[2]}}  I₂=${{counts[3]}}  R=${{counts[5]}}`;
    }});

    scrubber.value = idx;
    document.getElementById('frame-counter').textContent = `${{idx + 1}} / ${{FRAMES.length}}`;
}}

function play() {{
    if (timer) clearInterval(timer);
    playing = true;
    document.getElementById('btn-play').textContent = '⏸';
    document.getElementById('btn-play').classList.add('active');
    timer = setInterval(() => {{
        frameIdx++;
        if (frameIdx >= FRAMES.length) {{ frameIdx = 0; }}
        renderFrame(frameIdx);
    }}, Math.round(80 / speeds[speedIdx]));
}}

function pause() {{
    playing = false;
    if (timer) clearInterval(timer);
    timer = null;
    document.getElementById('btn-play').textContent = '▶';
    document.getElementById('btn-play').classList.remove('active');
}}

document.getElementById('btn-play').addEventListener('click', () => {{
    playing ? pause() : play();
}});

document.getElementById('btn-prev').addEventListener('click', () => {{
    pause();
    frameIdx = Math.max(0, frameIdx - 1);
    renderFrame(frameIdx);
}});

document.getElementById('btn-next').addEventListener('click', () => {{
    pause();
    frameIdx = Math.min(FRAMES.length - 1, frameIdx + 1);
    renderFrame(frameIdx);
}});

scrubber.addEventListener('input', (e) => {{
    pause();
    frameIdx = parseInt(e.target.value);
    renderFrame(frameIdx);
}});

document.getElementById('btn-slower').addEventListener('click', () => {{
    speedIdx = Math.max(0, speedIdx - 1);
    document.getElementById('speed-label').textContent = speeds[speedIdx] + '×';
    if (playing) play();
}});

document.getElementById('btn-faster').addEventListener('click', () => {{
    speedIdx = Math.min(speeds.length - 1, speedIdx + 1);
    document.getElementById('speed-label').textContent = speeds[speedIdx] + '×';
    if (playing) play();
}});

// Keyboard shortcuts
document.addEventListener('keydown', (e) => {{
    if (e.code === 'Space') {{ e.preventDefault(); playing ? pause() : play(); }}
    else if (e.code === 'ArrowLeft') {{ pause(); frameIdx = Math.max(0, frameIdx - 1); renderFrame(frameIdx); }}
    else if (e.code === 'ArrowRight') {{ pause(); frameIdx = Math.min(FRAMES.length - 1, frameIdx + 1); renderFrame(frameIdx); }}
    else if (e.code === 'Home') {{ pause(); frameIdx = 0; renderFrame(frameIdx); }}
    else if (e.code === 'End') {{ pause(); frameIdx = FRAMES.length - 1; renderFrame(frameIdx); }}
}});

// Initial render
renderFrame(0);
</script>
</body>
</html>"""
