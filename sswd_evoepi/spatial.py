"""Spatial and connectivity module.

Defines the metapopulation network: nodes (sites), connectivity matrices
(larval C and pathogen D), and the spatial simulation loop that couples
per-node dynamics with inter-node dispersal.

Core classes:
  - NodeDefinition: static attributes of a site (coordinates, habitat, etc.)
  - SpatialNode: runtime state wrapping definition + population + disease
  - MetapopulationNetwork: collection of nodes + C/D matrices

Core functions:
  - construct_larval_connectivity: build C matrix (exponential kernel, D_L=400km)
  - construct_pathogen_dispersal: build D matrix (exponential kernel, D_P=15km)
  - pathogen_dispersal_step: daily inter-node pathogen exchange (D.T @ P)
  - distribute_larvae: annual larval dispersal via C matrix
  - haversine_km: geodesic distance between nodes

References:
  - spatial-connectivity-spec.md §§1–9
  - CODE_ERRATA CE-1 through CE-11
  - ERRATA E3: flushing range [0.007, 0.05] for fjords
  - ERRATA E5: two connectivity matrices (C larval, D pathogen)

Build target: Phase 9 (Spatial Network).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from sswd_evoepi.types import (
    AGENT_DTYPE,
    LarvalCohort,
    SettlerPacket,
    N_LOCI,
    allocate_agents,
    allocate_genotypes,
)
from sswd_evoepi.environment import (
    make_sst_timeseries,
    sst_with_trend,
    seasonal_flushing,
    salinity_modifier,
    vbnc_fraction,
)


# ═══════════════════════════════════════════════════════════════════════
# GEODESIC DISTANCE
# ═══════════════════════════════════════════════════════════════════════

_EARTH_RADIUS_KM = 6371.0


def haversine_km(lat1: float, lon1: float,
                 lat2: float, lon2: float) -> float:
    """Great-circle distance in km between two (lat, lon) points.

    Args:
        lat1, lon1, lat2, lon2: Decimal degrees.

    Returns:
        Distance in kilometres.
    """
    d2r = np.pi / 180.0
    rlat1 = lat1 * d2r
    rlat2 = lat2 * d2r
    dlat = (lat2 - lat1) * d2r
    dlon = (lon2 - lon1) * d2r
    a = (np.sin(dlat / 2.0) ** 2
         + np.cos(rlat1) * np.cos(rlat2) * np.sin(dlon / 2.0) ** 2)
    c = 2.0 * np.arctan2(np.sqrt(a), np.sqrt(1.0 - a))
    return _EARTH_RADIUS_KM * c


def compute_distance_matrix(lats: np.ndarray,
                            lons: np.ndarray,
                            tortuosity: Optional[np.ndarray] = None
                            ) -> np.ndarray:
    """Compute pairwise waterway-distance matrix.

    Uses Haversine × tortuosity factor per pair.  For v1 we use a
    uniform tortuosity (default 1.5 — between open coast 1.2 and
    fjord 2.5) or per-pair values if provided.

    Args:
        lats: (N,) node latitudes.
        lons: (N,) node longitudes.
        tortuosity: Scalar or (N, N) array. Default 1.5.

    Returns:
        (N, N) waterway distance matrix (km).
    """
    N = len(lats)
    dist = np.zeros((N, N), dtype=np.float64)
    for i in range(N):
        for j in range(i + 1, N):
            d = haversine_km(lats[i], lons[i], lats[j], lons[j])
            dist[i, j] = d
            dist[j, i] = d

    if tortuosity is not None:
        if np.ndim(tortuosity) == 0:
            dist *= float(tortuosity)
        else:
            dist *= tortuosity
    else:
        dist *= 1.5  # default tortuosity
    return dist


def load_overwater_distances(
    node_defs: List[NodeDefinition],
    npz_path: str = 'results/overwater/distance_matrix_489.npz',
    max_match_km: float = 50.0,
) -> np.ndarray:
    """Load precomputed overwater distances for the given nodes.

    Matches node_defs to the 489-site matrix by nearest lat/lon.
    If a node is >max_match_km from any site, falls back to Haversine×1.5.
    Handles inf (disconnected) pairs by setting a large but finite distance.

    Args:
        node_defs: List of NodeDefinition to match.
        npz_path: Path to precomputed distance matrix.
        max_match_km: Maximum distance for matching nodes (km).

    Returns:
        (N, N) distance matrix in km.
    """
    from pathlib import Path
    
    if not Path(npz_path).exists():
        print(f"Warning: {npz_path} not found, falling back to Haversine×1.5")
        lats = np.array([nd.lat for nd in node_defs])
        lons = np.array([nd.lon for nd in node_defs])
        return compute_distance_matrix(lats, lons, tortuosity=1.5)
    
    # Load the 489×489 matrix
    data = np.load(npz_path)
    matrix_coords = data['coordinates']  # (489, 2) [lat, lon]
    matrix_distances = data['distances']  # (489, 489)
    matrix_names = data['names']  # (489,) for logging
    
    N = len(node_defs)
    result = np.zeros((N, N), dtype=np.float64)
    matched_indices = []
    
    # Match each node to nearest site in the 489-site matrix
    for i, node in enumerate(node_defs):
        node_coord = np.array([node.lat, node.lon])
        
        # Find nearest site by haversine distance
        distances_to_sites = np.array([
            haversine_km(node.lat, node.lon, coord[0], coord[1])
            for coord in matrix_coords
        ])
        
        nearest_idx = np.argmin(distances_to_sites)
        nearest_dist = distances_to_sites[nearest_idx]
        
        if nearest_dist <= max_match_km:
            matched_indices.append(nearest_idx)
            print(f"Matched {node.name} → {matrix_names[nearest_idx]} ({nearest_dist:.1f} km)")
        else:
            matched_indices.append(None)
            print(f"Warning: {node.name} >50 km from any site, using Haversine fallback")
    
    # Build the distance matrix
    fallback_pairs = []
    for i in range(N):
        for j in range(i + 1, N):
            idx_i = matched_indices[i]
            idx_j = matched_indices[j]
            
            if idx_i is not None and idx_j is not None:
                # Both nodes matched — use overwater distance
                overwater_dist = matrix_distances[idx_i, idx_j]
                
                if np.isinf(overwater_dist):
                    # Disconnected pair — set to large finite distance
                    # Use 10x the maximum finite distance in the matrix
                    max_finite = np.max(matrix_distances[np.isfinite(matrix_distances)])
                    overwater_dist = max_finite + 1000.0
                    print(f"Warning: {node_defs[i].name} ↔ {node_defs[j].name} disconnected, set to {overwater_dist:.0f} km")
                
                result[i, j] = overwater_dist
                result[j, i] = overwater_dist
            else:
                # At least one node unmatched — use Haversine fallback
                haversine_dist = haversine_km(
                    node_defs[i].lat, node_defs[i].lon,
                    node_defs[j].lat, node_defs[j].lon
                ) * 1.5  # apply standard tortuosity
                result[i, j] = haversine_dist
                result[j, i] = haversine_dist
                fallback_pairs.append((node_defs[i].name, node_defs[j].name))
    
    if fallback_pairs:
        print(f"Used Haversine fallback for {len(fallback_pairs)} pairs")
    
    return result


# ═══════════════════════════════════════════════════════════════════════
# NODE DEFINITION
# ═══════════════════════════════════════════════════════════════════════

@dataclass
class NodeDefinition:
    """Static definition for a single spatial node.

    Set once at initialisation.  Environmental time series and
    population state live on the wrapping SpatialNode, not here.
    """
    node_id: int
    name: str
    lat: float               # decimal degrees N
    lon: float               # decimal degrees E (west is negative)
    subregion: str           # e.g. "AK-SE", "SS", "CA-BJ"
    habitat_area: float      # m² of suitable benthic habitat
    carrying_capacity: int   # K = floor(habitat_area × ρ_max)
    is_fjord: bool = False
    sill_depth: float = np.inf   # m; np.inf for open coast
    flushing_rate: float = 0.5   # mean annual φ_k (d⁻¹)
    mean_sst: float = 10.0       # baseline annual mean SST (°C)
    sst_amplitude: float = 3.0   # annual cycle half-range (°C)
    sst_trend: float = 0.02      # °C/yr warming trend
    salinity: float = 32.0       # effective mean salinity (psu)
    depth_range: Tuple[float, float] = (5.0, 60.0)  # min, max depth (m)


# ═══════════════════════════════════════════════════════════════════════
# SPATIAL NODE (runtime wrapper)
# ═══════════════════════════════════════════════════════════════════════

@dataclass
class SpatialNode:
    """Runtime state for a single node, wrapping its definition.

    Population arrays (agents, genotypes) are owned here.
    Environmental quantities are updated each time step.
    """
    definition: NodeDefinition

    # ── Population arrays ────────────────────────────────────────────
    agents: Optional[np.ndarray] = None       # AGENT_DTYPE
    genotypes: Optional[np.ndarray] = None    # (max_agents, N_LOCI, 2)

    # ── Disease ──────────────────────────────────────────────────────
    vibrio_concentration: float = 0.0

    # ── Current environment (updated per timestep) ───────────────────
    current_sst: float = 10.0
    current_salinity: float = 32.0
    current_flushing: float = 0.5

    # ── SST time series (generated at init) ──────────────────────────
    sst_timeseries: Optional[np.ndarray] = None  # daily SST

    # ── Diagnostics ──────────────────────────────────────────────────
    disease_active: bool = False
    disease_free_years: int = 0

    # ── Convenience properties ───────────────────────────────────────

    @property
    def node_id(self) -> int:
        return self.definition.node_id

    @property
    def name(self) -> str:
        return self.definition.name

    @property
    def K(self) -> int:
        return self.definition.carrying_capacity

    @property
    def n_alive(self) -> int:
        if self.agents is None:
            return 0
        return int(np.sum(self.agents['alive']))


# ═══════════════════════════════════════════════════════════════════════
# CONNECTIVITY MATRICES
# ═══════════════════════════════════════════════════════════════════════

def construct_larval_connectivity(
    nodes: List[NodeDefinition],
    distances: np.ndarray,
    D_L: float = 400.0,
    alpha_self: Optional[np.ndarray] = None,
    barriers: Optional[Dict[Tuple[int, int], float]] = None,
    r_total: float = 0.02,
) -> np.ndarray:
    """Construct larval connectivity matrix C.

    C[j, k] = probability that a larva from node j settles at node k.

    Exponential distance kernel with self-recruitment on diagonal,
    barrier attenuation for non-adjacent coastal segments, and
    row-normalisation to r_total.

    Args:
        nodes: List of NodeDefinition (length N).
        distances: (N, N) waterway distance matrix (km).
        D_L: Larval dispersal scale (km).  Default 400.
        alpha_self: (N,) per-node self-recruitment fraction.
            Default: 0.30 for fjord nodes, 0.10 for open coast.
        barriers: Dict mapping (i, j) pairs to attenuation factors [0, 1].
        r_total: Target row-sum (total settlement success fraction).

    Returns:
        (N, N) dense connectivity matrix (float64).
    """
    N = len(nodes)
    if alpha_self is None:
        alpha_self = np.array(
            [0.30 if n.is_fjord else 0.10 for n in nodes],
            dtype=np.float64,
        )

    C = np.zeros((N, N), dtype=np.float64)

    for j in range(N):
        for k in range(N):
            if j == k:
                C[j, k] = alpha_self[j]
            else:
                kernel = np.exp(-distances[j, k] / D_L)
                # Barrier attenuation
                b = 1.0
                if barriers:
                    if (j, k) in barriers:
                        b = barriers[(j, k)]
                    elif (k, j) in barriers:
                        b = barriers[(k, j)]
                C[j, k] = (1.0 - alpha_self[j]) * kernel * b

    # Row-normalise so each row sums to r_total
    for j in range(N):
        row_sum = C[j, :].sum()
        if row_sum > 0:
            C[j, :] *= r_total / row_sum

    return C


def construct_pathogen_dispersal(
    nodes: List[NodeDefinition],
    distances: np.ndarray,
    D_P: float = 15.0,
    f_out: float = 0.2,
    max_range: float = 50.0,
) -> np.ndarray:
    """Construct pathogen dispersal matrix D.

    D[j, k] = fraction of pathogen at node j reaching node k per day.

    Short-range exponential kernel modulated by source flushing rate
    and sill attenuation for fjord barriers.

    Args:
        nodes: List of NodeDefinition (length N).
        distances: (N, N) waterway distance matrix (km).
        D_P: Pathogen dispersal scale (km). Default 15.
        f_out: Fraction of flushed water reaching neighbours. Default 0.2.
        max_range: Maximum dispersal range (km). Default 50.

    Returns:
        (N, N) dense dispersal matrix (float64).
    """
    N = len(nodes)
    D = np.zeros((N, N), dtype=np.float64)

    for j in range(N):
        phi_j = nodes[j].flushing_rate
        for k in range(N):
            if j == k:
                continue
            d_jk = distances[j, k]
            if d_jk > max_range or d_jk <= 0:
                continue

            kernel = np.exp(-d_jk / D_P)

            # Sill attenuation
            sill = min(nodes[j].sill_depth, nodes[k].sill_depth)
            max_depth = max(nodes[j].depth_range[1], nodes[k].depth_range[1])
            S_jk = _sill_attenuation(sill, max_depth)

            D[j, k] = phi_j * f_out * kernel * S_jk

        # Cap total export at flushing rate
        row_sum = D[j, :].sum()
        if row_sum > phi_j and phi_j > 0:
            D[j, :] *= phi_j / row_sum

    return D


def _sill_attenuation(sill_depth: float, max_depth: float) -> float:
    """Pathogen dispersal attenuation through a fjord sill.

    S = (sill / max_depth)²  ∈ [0, 1].
    Returns 1.0 if no sill (sill_depth >= max_depth or inf).
    """
    if sill_depth >= max_depth or np.isinf(sill_depth):
        return 1.0
    if max_depth <= 0:
        return 1.0
    ratio = sill_depth / max_depth
    return min(1.0, ratio * ratio)


# ═══════════════════════════════════════════════════════════════════════
# PATHOGEN DISPERSAL STEP (daily)
# ═══════════════════════════════════════════════════════════════════════

def pathogen_dispersal_step(P: np.ndarray, D: np.ndarray) -> np.ndarray:
    """Compute pathogen dispersal input to each node.

    dispersal_input[k] = Σ_j D[j, k] × P[j]

    This is D.T @ P.

    Args:
        P: (N,) Vibrio concentration at each node.
        D: (N, N) pathogen dispersal matrix.

    Returns:
        (N,) dispersal input to each node.
    """
    return D.T @ P


# ═══════════════════════════════════════════════════════════════════════
# LARVAL DISPERSAL (annual)
# ═══════════════════════════════════════════════════════════════════════

def distribute_larvae(
    source_node_ids: List[int],
    source_n_larvae: List[int],
    source_genotypes: List[np.ndarray],
    C: np.ndarray,
    rng: np.random.Generator,
) -> Dict[int, List[Tuple[np.ndarray, int]]]:
    """Distribute larvae from source nodes to receiving nodes via C matrix.

    For each source node j with n_competent larvae:
      - Total settlement probability = sum(C[j, :])
      - Binomial draw for how many settle at all
      - Multinomial assignment to destination nodes
      - Genotypes carried to destination

    Args:
        source_node_ids: Source node indices.
        source_n_larvae: Number of competent larvae per source.
        source_genotypes: Genotype arrays per source; each (n, N_LOCI, 2).
        C: (N, N) larval connectivity matrix.
        rng: Random generator.

    Returns:
        Dict {dest_node_id: [(genotypes_array, source_node_id), ...]}.
    """
    N = C.shape[0]
    settlers: Dict[int, List[Tuple[np.ndarray, int]]] = {k: [] for k in range(N)}

    for idx, j in enumerate(source_node_ids):
        n = source_n_larvae[idx]
        geno = source_genotypes[idx]
        if n == 0 or geno is None or len(geno) == 0:
            continue

        probs = C[j, :].copy()
        total_p = probs.sum()
        if total_p < 1e-12:
            continue

        # How many larvae settle at all?
        n_settling = rng.binomial(n, min(total_p, 1.0))
        if n_settling == 0:
            continue

        # Conditional probabilities for destination
        probs_norm = probs / total_p

        # Multinomial assignment
        dest_counts = rng.multinomial(n_settling, probs_norm)

        for k in range(N):
            nk = dest_counts[k]
            if nk <= 0:
                continue
            # Sample genotypes with replacement from source
            idx_sample = rng.choice(len(geno), size=nk, replace=True)
            settlers[k].append((geno[idx_sample].copy(), j))

    return settlers


# ═══════════════════════════════════════════════════════════════════════
# METAPOPULATION NETWORK
# ═══════════════════════════════════════════════════════════════════════

@dataclass
class MetapopulationNetwork:
    """Collection of spatial nodes + connectivity matrices.

    This is the top-level spatial object that orchestrates the
    multi-node simulation.
    """
    nodes: List[SpatialNode]
    C: np.ndarray                  # (N, N) larval connectivity
    D: np.ndarray                  # (N, N) pathogen dispersal
    distances: np.ndarray          # (N, N) waterway distances (km)
    rng: np.random.Generator = field(
        default_factory=lambda: np.random.default_rng(42)
    )

    @property
    def n_nodes(self) -> int:
        return len(self.nodes)

    def get_node(self, idx: int) -> SpatialNode:
        return self.nodes[idx]

    def get_vibrio_array(self) -> np.ndarray:
        """Current Vibrio concentration at all nodes. Shape: (N,)."""
        return np.array([n.vibrio_concentration for n in self.nodes])

    def get_sst_array(self) -> np.ndarray:
        """Current SST at all nodes. Shape: (N,)."""
        return np.array([n.current_sst for n in self.nodes])

    def get_population_array(self) -> np.ndarray:
        """Current alive count at all nodes. Shape: (N,)."""
        return np.array([n.n_alive for n in self.nodes])

    def total_population(self) -> int:
        return sum(n.n_alive for n in self.nodes)

    def summary(self) -> str:
        """Human-readable summary of the network."""
        lines = [f"MetapopulationNetwork: {self.n_nodes} nodes"]
        for n in self.nodes:
            d = n.definition
            lines.append(
                f"  [{d.node_id}] {d.name}: pop={n.n_alive}/{d.carrying_capacity} "
                f"SST={n.current_sst:.1f}°C sal={n.current_salinity:.0f}psu "
                f"φ={n.current_flushing:.3f} V={n.vibrio_concentration:.1f} "
                f"{'fjord' if d.is_fjord else 'coast'}"
            )
        return "\n".join(lines)


# ═══════════════════════════════════════════════════════════════════════
# NETWORK BUILDER
# ═══════════════════════════════════════════════════════════════════════

def build_network(
    node_defs: List[NodeDefinition],
    D_L: float = 400.0,
    D_P: float = 15.0,
    r_total: float = 0.02,
    f_out: float = 0.2,
    barriers: Optional[Dict[Tuple[int, int], float]] = None,
    tortuosity: float = 1.5,
    seed: int = 42,
    overwater_npz: Optional[str] = None,
) -> MetapopulationNetwork:
    """Build a MetapopulationNetwork from node definitions.

    Computes distance matrix, C matrix, D matrix, and wraps everything
    in SpatialNode objects.

    Args:
        node_defs: List of NodeDefinition objects.
        D_L: Larval dispersal scale (km).
        D_P: Pathogen dispersal scale (km).
        r_total: Total larval settlement success.
        f_out: Pathogen dispersal export fraction.
        barriers: Barrier attenuation dict for C matrix.
        tortuosity: Coastline tortuosity multiplier.
        seed: RNG seed.
        overwater_npz: Optional path to precomputed overwater distances.
            If provided, uses load_overwater_distances instead of Haversine×tortuosity.

    Returns:
        Assembled MetapopulationNetwork.
    """
    # Choose distance computation method
    if overwater_npz is not None:
        distances = load_overwater_distances(node_defs, overwater_npz)
    else:
        lats = np.array([nd.lat for nd in node_defs])
        lons = np.array([nd.lon for nd in node_defs])
        distances = compute_distance_matrix(lats, lons, tortuosity=tortuosity)

    C = construct_larval_connectivity(
        node_defs, distances, D_L=D_L, barriers=barriers, r_total=r_total,
    )
    D = construct_pathogen_dispersal(
        node_defs, distances, D_P=D_P, f_out=f_out,
    )

    nodes = []
    for nd in node_defs:
        sn = SpatialNode(
            definition=nd,
            current_sst=nd.mean_sst,
            current_salinity=nd.salinity,
            current_flushing=nd.flushing_rate,
        )
        nodes.append(sn)

    rng = np.random.default_rng(seed)

    return MetapopulationNetwork(
        nodes=nodes, C=C, D=D, distances=distances, rng=rng,
    )


# ═══════════════════════════════════════════════════════════════════════
# 5-NODE TEST NETWORK
# ═══════════════════════════════════════════════════════════════════════

def make_5node_network(seed: int = 42, use_overwater: bool = False) -> MetapopulationNetwork:
    """Create the canonical 5-node test network.

    Nodes:
      0: Sitka, AK — open coast, cold, large population
      1: Howe Sound, BC — fjord refugium, low flushing
      2: San Juan Islands, WA — Salish Sea, moderate
      3: Newport, OR — open coast, moderate
      4: Monterey, CA — warm, near functional extinction

    Args:
        seed: Random seed for network construction.
        use_overwater: If True, use precomputed overwater distances.
            If False, use Haversine×1.5 (default for backward compatibility).

    Returns:
        MetapopulationNetwork ready for simulation.
    """
    node_defs = get_5node_definitions()
    overwater_npz = 'results/overwater/distance_matrix_489.npz' if use_overwater else None
    return build_network(node_defs, seed=seed, overwater_npz=overwater_npz)


def get_5node_definitions() -> List[NodeDefinition]:
    """Return the 5 canonical test node definitions."""
    return [
        NodeDefinition(
            node_id=0,
            name="Sitka, AK",
            lat=57.06, lon=-135.34,
            subregion="AK-SE",
            habitat_area=50000.0,
            carrying_capacity=1000,
            is_fjord=False,
            sill_depth=np.inf,
            flushing_rate=0.8,
            mean_sst=8.0,
            sst_amplitude=3.5,
            sst_trend=0.015,
            salinity=32.0,
            depth_range=(5.0, 60.0),
        ),
        NodeDefinition(
            node_id=1,
            name="Howe Sound, BC",
            lat=49.52, lon=-123.25,
            subregion="SS",
            habitat_area=20000.0,
            carrying_capacity=400,
            is_fjord=True,
            sill_depth=30.0,         # moderate sill
            flushing_rate=0.03,      # low flushing (fjord)
            mean_sst=10.0,
            sst_amplitude=4.0,
            sst_trend=0.02,
            salinity=22.0,           # freshwater influence
            depth_range=(5.0, 100.0),
        ),
        NodeDefinition(
            node_id=2,
            name="San Juan Islands, WA",
            lat=48.53, lon=-123.02,
            subregion="SS",
            habitat_area=40000.0,
            carrying_capacity=800,
            is_fjord=False,
            sill_depth=np.inf,
            flushing_rate=0.3,       # semi-enclosed
            mean_sst=10.0,
            sst_amplitude=4.0,
            sst_trend=0.02,
            salinity=30.0,
            depth_range=(5.0, 40.0),
        ),
        NodeDefinition(
            node_id=3,
            name="Newport, OR",
            lat=44.63, lon=-124.05,
            subregion="WA-OR",
            habitat_area=30000.0,
            carrying_capacity=600,
            is_fjord=False,
            sill_depth=np.inf,
            flushing_rate=1.0,       # open coast
            mean_sst=12.0,
            sst_amplitude=3.0,
            sst_trend=0.02,
            salinity=33.0,
            depth_range=(5.0, 50.0),
        ),
        NodeDefinition(
            node_id=4,
            name="Monterey, CA",
            lat=36.62, lon=-121.90,
            subregion="CA-BJ",
            habitat_area=35000.0,
            carrying_capacity=700,
            is_fjord=False,
            sill_depth=np.inf,
            flushing_rate=0.8,       # open coast
            mean_sst=14.0,
            sst_amplitude=2.5,
            sst_trend=0.025,
            salinity=33.5,
            depth_range=(10.0, 60.0),
        ),
    ]


# ═══════════════════════════════════════════════════════════════════════
# YAML I/O
# ═══════════════════════════════════════════════════════════════════════

def save_node_definitions_yaml(node_defs: List[NodeDefinition],
                               path: str) -> None:
    """Save node definitions to YAML file.

    Args:
        node_defs: List of NodeDefinition objects.
        path: Output file path.
    """
    import yaml

    data = {"nodes": []}
    for nd in node_defs:
        entry = {
            "node_id": nd.node_id,
            "name": nd.name,
            "lat": nd.lat,
            "lon": nd.lon,
            "subregion": nd.subregion,
            "habitat_area": nd.habitat_area,
            "carrying_capacity": nd.carrying_capacity,
            "is_fjord": nd.is_fjord,
            "sill_depth": float(nd.sill_depth) if not np.isinf(nd.sill_depth) else None,
            "flushing_rate": nd.flushing_rate,
            "mean_sst": nd.mean_sst,
            "sst_amplitude": nd.sst_amplitude,
            "sst_trend": nd.sst_trend,
            "salinity": nd.salinity,
            "depth_range": list(nd.depth_range),
        }
        data["nodes"].append(entry)

    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, 'w') as f:
        yaml.dump(data, f, default_flow_style=False, sort_keys=False)


def load_node_definitions_yaml(path: str) -> List[NodeDefinition]:
    """Load node definitions from YAML file.

    Args:
        path: Input file path.

    Returns:
        List of NodeDefinition objects.
    """
    import yaml

    with open(path) as f:
        data = yaml.safe_load(f)

    nodes = []
    for entry in data["nodes"]:
        nd = NodeDefinition(
            node_id=entry["node_id"],
            name=entry["name"],
            lat=entry["lat"],
            lon=entry["lon"],
            subregion=entry["subregion"],
            habitat_area=entry["habitat_area"],
            carrying_capacity=entry["carrying_capacity"],
            is_fjord=entry.get("is_fjord", False),
            sill_depth=entry["sill_depth"] if entry.get("sill_depth") is not None else np.inf,
            flushing_rate=entry.get("flushing_rate", 0.5),
            mean_sst=entry.get("mean_sst", 10.0),
            sst_amplitude=entry.get("sst_amplitude", 3.0),
            sst_trend=entry.get("sst_trend", 0.02),
            salinity=entry.get("salinity", 32.0),
            depth_range=tuple(entry.get("depth_range", [5.0, 60.0])),
        )
        nodes.append(nd)

    return nodes
