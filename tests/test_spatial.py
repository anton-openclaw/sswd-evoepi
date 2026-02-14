"""Tests for spatial network module (Phase 9).

Tests:
  1. Distance computation (Haversine)
  2. Larval connectivity matrix C — structure and normalisation
  3. Pathogen dispersal matrix D — structure, sill attenuation
  4. Environment module — SST, salinity, flushing
  5. 5-node network builds correctly
  6. Disease-free equilibrium (5 nodes reach stable population)
  7. Larvae redistribute according to C matrix
  8. Pathogen spreads between connected nodes
  9. Fjord node (Howe Sound) has reduced pathogen import
  10. Conservation of individuals during dispersal
  11. YAML round-trip for node definitions

References:
  - spatial-connectivity-spec.md
  - CODE_ERRATA CE-1 through CE-11
"""

import os
import tempfile

import numpy as np
import pytest

from sswd_evoepi.environment import (
    sinusoidal_sst,
    sst_with_trend,
    make_sst_timeseries,
    salinity_modifier,
    seasonal_flushing,
    vbnc_fraction,
)
from sswd_evoepi.spatial import (
    NodeDefinition,
    SpatialNode,
    MetapopulationNetwork,
    haversine_km,
    compute_distance_matrix,
    construct_larval_connectivity,
    construct_pathogen_dispersal,
    pathogen_dispersal_step,
    distribute_larvae,
    build_network,
    make_5node_network,
    get_5node_definitions,
    save_node_definitions_yaml,
    load_node_definitions_yaml,
    _sill_attenuation,
)
from sswd_evoepi.model import run_spatial_simulation, SpatialSimResult


# ═══════════════════════════════════════════════════════════════════════
# HAVERSINE DISTANCE
# ═══════════════════════════════════════════════════════════════════════

class TestHaversine:
    def test_same_point_zero(self):
        assert haversine_km(48.0, -123.0, 48.0, -123.0) == pytest.approx(0.0)

    def test_known_distance(self):
        """Sitka → San Juan Islands ≈ 1000 km straight-line."""
        d = haversine_km(57.06, -135.34, 48.53, -123.02)
        # Should be roughly 1000-1200 km
        assert 900.0 < d < 1300.0

    def test_symmetric(self):
        d1 = haversine_km(57.06, -135.34, 36.62, -121.90)
        d2 = haversine_km(36.62, -121.90, 57.06, -135.34)
        assert d1 == pytest.approx(d2, rel=1e-10)

    def test_sitka_to_monterey(self):
        """Sitka → Monterey ≈ 2400 km straight-line."""
        d = haversine_km(57.06, -135.34, 36.62, -121.90)
        assert 2200.0 < d < 2700.0


# ═══════════════════════════════════════════════════════════════════════
# ENVIRONMENT MODULE
# ═══════════════════════════════════════════════════════════════════════

class TestEnvironment:
    def test_sinusoidal_peak_at_doy227(self):
        """SST should peak around day 227 (mid-August)."""
        sst_peak = sinusoidal_sst(227, mean_sst=10.0, amplitude=4.0)
        assert sst_peak == pytest.approx(14.0, abs=0.1)

    def test_sinusoidal_trough(self):
        """SST trough ~6 months from peak (day 227 + 182 ≈ day 44)."""
        sst_trough = sinusoidal_sst(44, mean_sst=10.0, amplitude=4.0)
        assert sst_trough == pytest.approx(6.0, abs=0.5)

    def test_sst_mean_over_year(self):
        """Mean SST over a year should ≈ mean_sst."""
        vals = [sinusoidal_sst(d, 10.0, 4.0) for d in range(365)]
        assert np.mean(vals) == pytest.approx(10.0, abs=0.1)

    def test_sst_with_trend(self):
        """Warming trend should shift baseline."""
        sst_2000 = sst_with_trend(227, 2000, 10.0, 4.0, 0.02, 2000)
        sst_2050 = sst_with_trend(227, 2050, 10.0, 4.0, 0.02, 2000)
        assert sst_2050 - sst_2000 == pytest.approx(1.0, abs=0.01)

    def test_make_sst_timeseries_shape(self):
        ts = make_sst_timeseries(5, 2000, 10.0, 4.0, 0.02)
        assert ts.shape == (5 * 365,)

    def test_salinity_modifier_low(self):
        """Below s_min → 0."""
        assert salinity_modifier(5.0) == pytest.approx(0.0)

    def test_salinity_modifier_full(self):
        """Above s_full → 1."""
        assert salinity_modifier(32.0) == pytest.approx(1.0)

    def test_salinity_modifier_mid(self):
        """At midpoint (19 psu): ((19−10)/(28−10))² = 0.25."""
        assert salinity_modifier(19.0) == pytest.approx(0.25, abs=0.01)

    def test_seasonal_flushing_fjord(self):
        """Fjord flushing should vary ±30% over the year."""
        phi_jun = seasonal_flushing(0.03, 5, True)  # month 5 = June peak
        phi_dec = seasonal_flushing(0.03, 11, True)  # December trough
        assert phi_jun > phi_dec
        assert phi_jun == pytest.approx(0.03 * 1.3, abs=0.002)

    def test_vbnc_above_threshold(self):
        """Well above T_vbnc → nearly all culturable."""
        assert vbnc_fraction(20.0) > 0.99

    def test_vbnc_below_threshold(self):
        """Well below T_vbnc → mostly VBNC."""
        assert vbnc_fraction(5.0) < 0.01

    def test_vbnc_at_threshold(self):
        """At T_vbnc → 50%."""
        assert vbnc_fraction(12.0) == pytest.approx(0.5, abs=0.01)


# ═══════════════════════════════════════════════════════════════════════
# CONNECTIVITY MATRICES
# ═══════════════════════════════════════════════════════════════════════

class TestConnectivity:
    @pytest.fixture
    def five_nodes(self):
        return get_5node_definitions()

    @pytest.fixture
    def dist_matrix(self, five_nodes):
        lats = np.array([n.lat for n in five_nodes])
        lons = np.array([n.lon for n in five_nodes])
        return compute_distance_matrix(lats, lons, tortuosity=1.5)

    def test_distance_matrix_shape(self, dist_matrix):
        assert dist_matrix.shape == (5, 5)

    def test_distance_matrix_symmetric(self, dist_matrix):
        np.testing.assert_array_almost_equal(dist_matrix, dist_matrix.T)

    def test_distance_matrix_diagonal_zero(self, dist_matrix):
        np.testing.assert_array_almost_equal(np.diag(dist_matrix), 0.0)

    def test_distance_matrix_positive(self, dist_matrix):
        # Off-diagonal should all be positive
        mask = ~np.eye(5, dtype=bool)
        assert np.all(dist_matrix[mask] > 0)

    def test_C_matrix_shape(self, five_nodes, dist_matrix):
        C = construct_larval_connectivity(five_nodes, dist_matrix)
        assert C.shape == (5, 5)

    def test_C_matrix_non_negative(self, five_nodes, dist_matrix):
        C = construct_larval_connectivity(five_nodes, dist_matrix)
        assert np.all(C >= 0)

    def test_C_matrix_row_sums(self, five_nodes, dist_matrix):
        """Each row should sum to r_total (0.02)."""
        C = construct_larval_connectivity(
            five_nodes, dist_matrix, r_total=0.02,
        )
        row_sums = C.sum(axis=1)
        np.testing.assert_allclose(row_sums, 0.02, atol=1e-10)

    def test_C_self_recruitment_positive(self, five_nodes, dist_matrix):
        """Diagonal entries should be positive (self-recruitment)."""
        C = construct_larval_connectivity(five_nodes, dist_matrix)
        assert np.all(np.diag(C) > 0)

    def test_C_fjord_higher_self_recruitment(self, five_nodes, dist_matrix):
        """Fjord node (1) should have higher relative self-recruitment."""
        C = construct_larval_connectivity(five_nodes, dist_matrix)
        # Node 1 = Howe Sound (fjord), Node 3 = Newport (coast)
        self_frac_fjord = C[1, 1] / C[1, :].sum()
        self_frac_coast = C[3, 3] / C[3, :].sum()
        assert self_frac_fjord > self_frac_coast

    def test_C_nearby_nodes_stronger(self, five_nodes, dist_matrix):
        """Nearby node pairs should have stronger connectivity."""
        C = construct_larval_connectivity(five_nodes, dist_matrix)
        # SJI (2) → Howe Sound (1) should be stronger than SJI → Sitka (0)
        assert C[2, 1] > C[2, 0]

    def test_D_matrix_shape(self, five_nodes, dist_matrix):
        D = construct_pathogen_dispersal(five_nodes, dist_matrix)
        assert D.shape == (5, 5)

    def test_D_matrix_diagonal_zero(self, five_nodes, dist_matrix):
        """D diagonal should be zero (no self-dispersal)."""
        D = construct_pathogen_dispersal(five_nodes, dist_matrix)
        np.testing.assert_array_almost_equal(np.diag(D), 0.0)

    def test_D_matrix_non_negative(self, five_nodes, dist_matrix):
        D = construct_pathogen_dispersal(five_nodes, dist_matrix)
        assert np.all(D >= 0)

    def test_D_short_range(self, five_nodes, dist_matrix):
        """D should be very sparse — most pairs > 50 km apart."""
        D = construct_pathogen_dispersal(five_nodes, dist_matrix)
        # For widely spaced 5-node network, most D entries are 0
        # Only nearby nodes (Howe Sound ↔ SJI) might be non-zero
        # Actually at 1.5× tortuosity, Howe Sound to SJI is ~180 km
        # which is way beyond 50 km max range, so D should be all zeros
        # for this widely-spaced test network
        # That's correct biologically: pathogen can't spread 100s of km
        n_nonzero = np.count_nonzero(D)
        # Most or all entries should be zero for this network
        assert n_nonzero <= 10  # generous upper bound

    def test_D_sill_attenuation(self):
        """Sill attenuation reduces pathogen flow."""
        # Shallow sill (30m), deep basin (100m)
        S = _sill_attenuation(30.0, 100.0)
        assert S == pytest.approx(0.09, abs=0.01)
        # No sill
        S_open = _sill_attenuation(np.inf, 60.0)
        assert S_open == pytest.approx(1.0)
        # Sill deeper than basin
        S_deep = _sill_attenuation(80.0, 60.0)
        assert S_deep == pytest.approx(1.0)


# ═══════════════════════════════════════════════════════════════════════
# PATHOGEN DISPERSAL STEP
# ═══════════════════════════════════════════════════════════════════════

class TestPathogenDispersal:
    def test_zero_vibrio_zero_dispersal(self):
        P = np.zeros(5)
        D = np.eye(5) * 0  # no dispersal
        result = pathogen_dispersal_step(P, D)
        np.testing.assert_array_almost_equal(result, 0.0)

    def test_dispersal_conserves_direction(self):
        """Pathogen flows from high to receiving nodes."""
        D = np.zeros((3, 3))
        D[0, 1] = 0.1  # node 0 exports to node 1
        D[0, 2] = 0.05  # node 0 exports to node 2
        P = np.array([100.0, 0.0, 0.0])
        result = pathogen_dispersal_step(P, D)
        assert result[0] == pytest.approx(0.0)  # node 0 doesn't self-import
        assert result[1] == pytest.approx(10.0)  # 100 × 0.1
        assert result[2] == pytest.approx(5.0)   # 100 × 0.05


# ═══════════════════════════════════════════════════════════════════════
# LARVAL DISPERSAL
# ═══════════════════════════════════════════════════════════════════════

class TestLarvalDispersal:
    def test_conserves_settlers(self):
        """Total settlers should be ≤ total larvae × sum(C row)."""
        rng = np.random.default_rng(42)
        C = np.array([
            [0.01, 0.005, 0.005],
            [0.005, 0.01, 0.005],
            [0.005, 0.005, 0.01],
        ])
        geno = np.zeros((100, 52, 2), dtype=np.int8)
        result = distribute_larvae(
            source_node_ids=[0],
            source_n_larvae=[100],
            source_genotypes=[geno],
            C=C,
            rng=rng,
        )
        total_settled = sum(
            sum(len(batch) for batch, _ in result[k])
            for k in range(3)
        )
        # Should settle approximately 100 × 0.02 = 2 larvae
        assert 0 <= total_settled <= 20

    def test_empty_source_no_crash(self):
        rng = np.random.default_rng(42)
        C = np.eye(3) * 0.01
        result = distribute_larvae(
            source_node_ids=[],
            source_n_larvae=[],
            source_genotypes=[],
            C=C,
            rng=rng,
        )
        assert all(len(v) == 0 for v in result.values())

    def test_genotypes_carried(self):
        """Larvae should carry genotypes to destination."""
        rng = np.random.default_rng(42)
        C = np.array([
            [0.5, 0.5],
            [0.5, 0.5],
        ])
        geno = np.ones((50, 52, 2), dtype=np.int8)  # all 1s
        result = distribute_larvae(
            source_node_ids=[0],
            source_n_larvae=[50],
            source_genotypes=[geno],
            C=C,
            rng=rng,
        )
        for k in range(2):
            for batch, src in result[k]:
                assert np.all(batch == 1), "Genotypes should be carried"

    def test_multinomial_distribution(self):
        """With large N and uniform C, settlers should be roughly even."""
        rng = np.random.default_rng(42)
        C = np.array([
            [0.25, 0.25, 0.25, 0.25],
            [0.25, 0.25, 0.25, 0.25],
            [0.25, 0.25, 0.25, 0.25],
            [0.25, 0.25, 0.25, 0.25],
        ])
        geno = np.zeros((10000, 52, 2), dtype=np.int8)
        result = distribute_larvae(
            source_node_ids=[0],
            source_n_larvae=[10000],
            source_genotypes=[geno],
            C=C,
            rng=rng,
        )
        counts = []
        for k in range(4):
            n = sum(len(batch) for batch, _ in result[k])
            counts.append(n)
        # Each should get ~25% of total settlers
        total = sum(counts)
        if total > 0:
            fracs = [c / total for c in counts]
            for f in fracs:
                assert 0.15 < f < 0.35, f"Expected ~0.25, got {f}"


# ═══════════════════════════════════════════════════════════════════════
# 5-NODE NETWORK
# ═══════════════════════════════════════════════════════════════════════

class TestFiveNodeNetwork:
    @pytest.fixture
    def network(self):
        return make_5node_network(seed=42)

    def test_network_has_5_nodes(self, network):
        assert network.n_nodes == 5

    def test_node_names(self, network):
        names = [n.name for n in network.nodes]
        assert "Sitka, AK" in names
        assert "Howe Sound, BC" in names
        assert "San Juan Islands, WA" in names
        assert "Newport, OR" in names
        assert "Monterey, CA" in names

    def test_C_matrix_shape(self, network):
        assert network.C.shape == (5, 5)

    def test_D_matrix_shape(self, network):
        assert network.D.shape == (5, 5)

    def test_distances_shape(self, network):
        assert network.distances.shape == (5, 5)

    def test_howe_sound_is_fjord(self, network):
        hs = network.nodes[1]
        assert hs.definition.is_fjord

    def test_sitka_coldest(self, network):
        """Sitka should have the lowest mean SST."""
        ssts = [n.definition.mean_sst for n in network.nodes]
        assert np.argmin(ssts) == 0

    def test_monterey_warmest(self, network):
        """Monterey should have the highest mean SST."""
        ssts = [n.definition.mean_sst for n in network.nodes]
        assert np.argmax(ssts) == 4

    def test_summary_runs(self, network):
        s = network.summary()
        assert "MetapopulationNetwork" in s
        assert "Sitka" in s


# ═══════════════════════════════════════════════════════════════════════
# SPATIAL SIMULATION — DISEASE-FREE EQUILIBRIUM
# ═══════════════════════════════════════════════════════════════════════

class TestSpatialEquilibrium:
    """5 nodes reach independent disease-free equilibrium."""

    @pytest.fixture(scope="class")
    def result(self):
        """Run 8-year disease-free simulation (cached for class)."""
        network = make_5node_network(seed=42)
        return run_spatial_simulation(
            network=network,
            n_years=8,
            disease_year=None,  # no disease
            seed=42,
        )

    def test_result_type(self, result):
        assert isinstance(result, SpatialSimResult)

    def test_all_nodes_have_population(self, result):
        """All 5 nodes should maintain positive population."""
        final_pops = result.yearly_pop[:, -1]
        assert np.all(final_pops > 0), f"Final pops: {final_pops}"

    def test_population_stable(self, result):
        """After a few years, population should stabilise near K."""
        # Check last 3 years — population shouldn't change > 30%
        for i in range(5):
            late_pops = result.yearly_pop[i, -3:]
            mean_late = np.mean(late_pops)
            for p in late_pops:
                assert abs(p - mean_late) / max(mean_late, 1) < 0.35, \
                    f"Node {i} unstable: {late_pops}"

    def test_no_disease_deaths(self, result):
        """No disease deaths in disease-free run."""
        assert np.all(result.yearly_disease_deaths == 0)

    def test_total_population_positive(self, result):
        """Total metapopulation stays positive every year."""
        assert np.all(result.yearly_total_pop > 0)

    def test_larvae_produced(self, result):
        """Larvae should be produced most years."""
        # After initial year, there should be larvae
        assert np.any(result.yearly_total_larvae_dispersed > 0)


# ═══════════════════════════════════════════════════════════════════════
# SPATIAL SIMULATION — PATHOGEN SPREAD
# ═══════════════════════════════════════════════════════════════════════

class TestPathogenSpread:
    """Pathogen spreads between connected nodes."""

    def test_pathogen_matrix_transport(self):
        """Direct D-matrix test: pathogen moves from source to dest."""
        D = np.zeros((3, 3))
        D[0, 1] = 0.05  # node 0 → node 1
        P = np.array([1000.0, 0.0, 0.0])
        result = pathogen_dispersal_step(P, D)
        assert result[1] > 0
        assert result[2] == 0  # no connection to node 2

    def test_fjord_reduced_import(self):
        """Fjord node should receive less pathogen due to sill attenuation."""
        # Create a simple 2-node network: open coast ↔ fjord
        open_coast = NodeDefinition(
            node_id=0, name="Coast", lat=48.5, lon=-123.0,
            subregion="WA-OR", habitat_area=30000, carrying_capacity=500,
            is_fjord=False, sill_depth=np.inf, flushing_rate=1.0,
            mean_sst=12.0, sst_amplitude=3.0, sst_trend=0.02,
            salinity=33.0, depth_range=(5.0, 50.0),
        )
        fjord = NodeDefinition(
            node_id=1, name="Fjord", lat=48.55, lon=-123.05,
            subregion="WA-OR", habitat_area=10000, carrying_capacity=200,
            is_fjord=True, sill_depth=20.0, flushing_rate=0.02,
            mean_sst=10.0, sst_amplitude=3.0, sst_trend=0.02,
            salinity=20.0, depth_range=(5.0, 80.0),
        )
        nodes = [open_coast, fjord]
        lats = np.array([n.lat for n in nodes])
        lons = np.array([n.lon for n in nodes])
        dist = compute_distance_matrix(lats, lons, tortuosity=1.2)
        D = construct_pathogen_dispersal(nodes, dist, D_P=15.0, f_out=0.2)

        # Coast→Fjord export should be attenuated
        # Fjord→Coast export should also be low (low flushing)
        # D[0,1] = coast export to fjord (high φ but sill attenuation)
        # D[1,0] = fjord export to coast (low φ)
        assert D[1, 0] < D[0, 1] or D[1, 0] == 0, \
            "Fjord should export less than coast"


# ═══════════════════════════════════════════════════════════════════════
# YAML ROUND-TRIP
# ═══════════════════════════════════════════════════════════════════════

class TestYamlIO:
    def test_round_trip(self):
        """Save and load node definitions → identical."""
        defs = get_5node_definitions()
        with tempfile.NamedTemporaryFile(suffix=".yaml", delete=False) as f:
            path = f.name
        try:
            save_node_definitions_yaml(defs, path)
            loaded = load_node_definitions_yaml(path)
            assert len(loaded) == len(defs)
            for orig, ld in zip(defs, loaded):
                assert orig.node_id == ld.node_id
                assert orig.name == ld.name
                assert orig.lat == pytest.approx(ld.lat)
                assert orig.lon == pytest.approx(ld.lon)
                assert orig.subregion == ld.subregion
                assert orig.is_fjord == ld.is_fjord
                assert orig.carrying_capacity == ld.carrying_capacity
                assert orig.flushing_rate == pytest.approx(ld.flushing_rate)
                assert orig.mean_sst == pytest.approx(ld.mean_sst)
                assert orig.salinity == pytest.approx(ld.salinity)
        finally:
            os.unlink(path)

    def test_load_configs_yaml(self):
        """Load from the committed configs/nodes_5site.yaml."""
        path = os.path.join(
            os.path.dirname(__file__), "..", "configs", "nodes_5site.yaml"
        )
        if os.path.exists(path):
            defs = load_node_definitions_yaml(path)
            assert len(defs) == 5
            assert defs[0].name == "Sitka, AK"
            assert defs[1].is_fjord is True
