"""Lightweight monthly snapshot recorder for wavefront visualization.

Records per-node aggregate stats (population, infection counts) at
monthly intervals. Designed to be extremely low-overhead:
- Only samples every 30 days
- Stores just 2 int32 values per node per month
- Keeps data in memory, writes once at end

Usage:
    recorder = MonthlyRecorder(n_nodes=896)

    # In simulation daily loop:
    recorder.capture(sim_day, nodes)

    # After simulation:
    recorder.save("monthly_snapshots.npz", site_lats, site_lons, site_names, K)
"""

from __future__ import annotations
from pathlib import Path
from typing import List, Optional
import numpy as np


class MonthlyRecorder:
    """Captures monthly per-node population and infection counts.

    Infection = disease_state in {1 (E), 2 (I1), 3 (I2)}, alive only.
    """

    CAPTURE_INTERVAL = 30  # days between captures

    def __init__(self, n_nodes: int, sst_start_year: int = 2012):
        self.n_nodes = n_nodes
        self.sst_start_year = sst_start_year
        # Storage: list of (sim_day, pop_array, infected_array)
        self._frames: List[tuple] = []

    def capture(self, sim_day: int, nodes: list) -> None:
        """Capture one snapshot if it's a capture day.

        Args:
            sim_day: Global simulation day (year * 365 + day_of_year).
            nodes: List of SpatialNode objects with .agents structured array.
        """
        if sim_day % self.CAPTURE_INTERVAL != 0:
            return

        pop = np.zeros(self.n_nodes, dtype=np.int32)
        infected = np.zeros(self.n_nodes, dtype=np.int32)

        for i, node in enumerate(nodes):
            if i >= self.n_nodes:
                break
            agents = node.agents
            alive = agents['alive'].astype(bool)
            # Exclude sentinels from population counts
            _pyc_alive = alive
            if 'is_sentinel' in agents.dtype.names:
                _pyc_alive = alive & ~agents['is_sentinel'].astype(bool)
            pop[i] = int(np.sum(_pyc_alive))
            ds = agents['disease_state'][_pyc_alive]
            # E=1, I1=2, I2=3 are "sick"
            infected[i] = int(np.sum((ds >= 1) & (ds <= 3)))

        self._frames.append((sim_day, pop.copy(), infected.copy()))

    @property
    def n_frames(self) -> int:
        return len(self._frames)

    def save(
        self,
        path: str,
        site_lats: np.ndarray,
        site_lons: np.ndarray,
        site_names: list,
        K: int,
    ) -> None:
        """Save all snapshots to a compact .npz file.

        Stored arrays:
            sim_days: (n_frames,) int32
            populations: (n_frames, n_nodes) int32
            infected: (n_frames, n_nodes) int32
            site_lats: (n_nodes,) float64
            site_lons: (n_nodes,) float64
            site_names: (n_nodes,) str (as object array)
            K: int scalar
            sst_start_year: int scalar
        """
        if not self._frames:
            return

        n = len(self._frames)
        sim_days = np.array([f[0] for f in self._frames], dtype=np.int32)
        populations = np.stack([f[1] for f in self._frames])
        infected_arr = np.stack([f[2] for f in self._frames])

        Path(path).parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            path,
            sim_days=sim_days,
            populations=populations,
            infected=infected_arr,
            site_lats=np.asarray(site_lats, dtype=np.float64),
            site_lons=np.asarray(site_lons, dtype=np.float64),
            site_names=np.array(site_names, dtype=object),
            K=np.int32(K),
            sst_start_year=np.int32(self.sst_start_year),
        )

    @staticmethod
    def load(path: str) -> dict:
        """Load snapshot data from .npz file. Returns dict of arrays."""
        data = np.load(path, allow_pickle=True)
        return {
            'sim_days': data['sim_days'],
            'populations': data['populations'],
            'infected': data['infected'],
            'site_lats': data['site_lats'],
            'site_lons': data['site_lons'],
            'site_names': data['site_names'],
            'K': int(data['K']),
            'sst_start_year': int(data['sst_start_year']),
        }
