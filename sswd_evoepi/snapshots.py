"""Optional individual-level snapshot recording.

Records (x, y, disease_state, resistance, size, alive) for every agent
at configurable intervals. Designed for the "wildfire" visualization —
watch the epidemic spread through individuals in real time.

Memory-conscious: stores compressed snapshots and can limit to specific
nodes or time windows.

Usage:
    recorder = SnapshotRecorder(
        enabled=True,
        interval_days=1,       # every day during epidemic
        nodes=[0, 1, 2, 3, 4],
    )
    
    # In simulation loop:
    recorder.capture(day, year, nodes)
    
    # After simulation:
    recorder.save("snapshots.npz")
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np


@dataclass
class IndividualSnapshot:
    """One snapshot: positions + states of all agents at one node at one time."""
    day: int        # simulation day (year * 365 + day_of_year)
    year: int
    node_id: int
    # Parallel arrays (only alive individuals)
    x: np.ndarray           # float32
    y: np.ndarray           # float32
    disease_state: np.ndarray  # int8
    resistance: np.ndarray  # float32
    size: np.ndarray        # float32
    n_alive: int


class SnapshotRecorder:
    """Records individual-level snapshots for wildfire visualization.
    
    When enabled=False, all methods are no-ops (zero overhead).
    """
    
    def __init__(
        self,
        enabled: bool = False,
        interval_days: int = 1,
        nodes: Optional[List[int]] = None,
        start_day: int = 0,
        end_day: int = 999999,
    ):
        """
        Args:
            enabled: Master switch. False = no-ops everywhere.
            interval_days: Capture every N days (1 = daily, 7 = weekly).
            nodes: List of node indices to capture (None = all).
            start_day: First simulation day to record.
            end_day: Last simulation day to record.
        """
        self.enabled = enabled
        self.interval_days = interval_days
        self.node_filter = set(nodes) if nodes is not None else None
        self.start_day = start_day
        self.end_day = end_day
        
        # Storage: dict of (sim_day, node_id) -> IndividualSnapshot
        self.snapshots: Dict[Tuple[int, int], IndividualSnapshot] = {}
    
    def should_capture(self, sim_day: int) -> bool:
        """Check if we should capture this day."""
        if not self.enabled:
            return False
        if sim_day < self.start_day or sim_day > self.end_day:
            return False
        return (sim_day % self.interval_days) == 0
    
    def capture_node(
        self,
        sim_day: int,
        year: int,
        node_id: int,
        agents: np.ndarray,
    ) -> None:
        """Capture a snapshot of one node's agents.
        
        Args:
            sim_day: Global simulation day (year * 365 + day_of_year).
            year: Current year.
            node_id: Node index.
            agents: Agent structured array.
        """
        if not self.enabled:
            return
        if self.node_filter is not None and node_id not in self.node_filter:
            return
        
        alive_mask = agents['alive'].astype(bool)
        n_alive = int(np.sum(alive_mask))
        
        snap = IndividualSnapshot(
            day=sim_day,
            year=year,
            node_id=node_id,
            x=agents['x'][alive_mask].copy(),
            y=agents['y'][alive_mask].copy(),
            disease_state=agents['disease_state'][alive_mask].copy(),
            resistance=agents['resistance'][alive_mask].copy(),
            size=agents['size'][alive_mask].copy(),
            n_alive=n_alive,
        )
        self.snapshots[(sim_day, node_id)] = snap
    
    def capture_all_nodes(self, sim_day: int, year: int, nodes: list) -> None:
        """Capture snapshots for all (filtered) nodes.
        
        Args:
            sim_day: Global simulation day.
            year: Current year.
            nodes: List of SpatialNode objects (must have .agents and .definition.node_id).
        """
        if not self.should_capture(sim_day):
            return
        for node in nodes:
            self.capture_node(sim_day, year, node.definition.node_id, node.agents)
    
    def get_days(self) -> List[int]:
        """Get sorted list of all captured simulation days."""
        return sorted(set(d for d, _ in self.snapshots.keys()))
    
    def get_nodes(self) -> List[int]:
        """Get sorted list of all captured node IDs."""
        return sorted(set(n for _, n in self.snapshots.keys()))
    
    def get_snapshot(self, sim_day: int, node_id: int) -> Optional[IndividualSnapshot]:
        """Get a specific snapshot."""
        return self.snapshots.get((sim_day, node_id))
    
    def save(self, path: str) -> None:
        """Save all snapshots to compressed npz file.
        
        Format: each snapshot stored as separate arrays with naming convention:
            d{day}_n{node}_x, d{day}_n{node}_y, etc.
        Plus metadata arrays: days, nodes, n_alive.
        """
        if not self.snapshots:
            return
        
        arrays = {}
        meta_days = []
        meta_nodes = []
        meta_years = []
        meta_n_alive = []
        
        for (day, nid), snap in sorted(self.snapshots.items()):
            prefix = f"d{day}_n{nid}"
            arrays[f"{prefix}_x"] = snap.x
            arrays[f"{prefix}_y"] = snap.y
            arrays[f"{prefix}_ds"] = snap.disease_state
            arrays[f"{prefix}_r"] = snap.resistance
            arrays[f"{prefix}_sz"] = snap.size
            meta_days.append(day)
            meta_nodes.append(nid)
            meta_years.append(snap.year)
            meta_n_alive.append(snap.n_alive)
        
        arrays['meta_days'] = np.array(meta_days, dtype=np.int32)
        arrays['meta_nodes'] = np.array(meta_nodes, dtype=np.int32)
        arrays['meta_years'] = np.array(meta_years, dtype=np.int32)
        arrays['meta_n_alive'] = np.array(meta_n_alive, dtype=np.int32)
        
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(path, **arrays)
    
    @classmethod
    def load(cls, path: str) -> 'SnapshotRecorder':
        """Load snapshots from npz file."""
        data = np.load(path)
        recorder = cls(enabled=False)  # Don't capture, just hold data
        
        days = data['meta_days']
        nodes = data['meta_nodes']
        years = data['meta_years']
        n_alive = data['meta_n_alive']
        
        for i in range(len(days)):
            d, n, y, na = int(days[i]), int(nodes[i]), int(years[i]), int(n_alive[i])
            prefix = f"d{d}_n{n}"
            snap = IndividualSnapshot(
                day=d, year=y, node_id=n,
                x=data[f"{prefix}_x"],
                y=data[f"{prefix}_y"],
                disease_state=data[f"{prefix}_ds"],
                resistance=data[f"{prefix}_r"],
                size=data[f"{prefix}_sz"],
                n_alive=na,
            )
            recorder.snapshots[(d, n)] = snap
        
        return recorder
    
    def memory_estimate_mb(self) -> float:
        """Estimate memory usage of stored snapshots."""
        total_bytes = 0
        for snap in self.snapshots.values():
            # 4 float32 arrays + 1 int8 array per alive individual
            total_bytes += snap.n_alive * (4 * 4 + 1)
        return total_bytes / (1024 * 1024)
