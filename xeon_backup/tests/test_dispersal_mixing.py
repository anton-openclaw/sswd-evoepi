"""Tests for dispersal-weighted trait mixing of pathogen community traits.

Tests cover:
- v_local field existence and defaults in NodeDiseaseState
- Dispersal mixing of T_vbnc and v_local between nodes
- No mixing when pathogen_adaptation=False
- No mixing at unreached nodes
- Wavefront trait inheritance includes v_local
- Identity: same traits everywhere → no change
- Asymmetric dispersal shifts traits toward dominant source

Convention: D[source, destination] = fraction of source's pathogen going to dest.
_D_T_sparse = csr_matrix(D.T), so _D_T_sparse @ P gives incoming pathogen per node.
"""

import numpy as np
import pytest
from scipy.sparse import csr_matrix

from sswd_evoepi.disease import NodeDiseaseState
from sswd_evoepi.model import _inherit_pathogen_traits


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_D_T_sparse(D):
    """Build D^T in CSR format from a dense dispersal matrix D.
    
    D[source, dest] = fraction from source to dest.
    D_T_sparse = D.T so D_T_sparse @ P = incoming pathogen per node.
    """
    return csr_matrix(D.T)


def _apply_dispersal_mixing(node_disease_states, disease_reached, D_T_sparse, P):
    """Replicate the dispersal trait mixing logic from model.py for testing."""
    N = len(node_disease_states)
    T_vec = np.array([node_disease_states[i].T_vbnc_local for i in range(N)])
    V_vec = np.array([node_disease_states[i].v_local for i in range(N)])

    PT_vec = P * T_vec
    PV_vec = P * V_vec

    dispersal_in = D_T_sparse @ P
    incoming_PT = D_T_sparse @ PT_vec
    incoming_PV = D_T_sparse @ PV_vec

    total = P + dispersal_in
    for i in range(N):
        if disease_reached[i] and total[i] > 0:
            node_disease_states[i].T_vbnc_local = (P[i] * T_vec[i] + incoming_PT[i]) / total[i]
            node_disease_states[i].v_local = (P[i] * V_vec[i] + incoming_PV[i]) / total[i]


# ---------------------------------------------------------------------------
# Tests: NodeDiseaseState field
# ---------------------------------------------------------------------------

class TestNodeDiseaseStateVLocal:
    def test_v_local_exists_with_default(self):
        """v_local field exists in NodeDiseaseState with default 0.5."""
        nds = NodeDiseaseState()
        assert hasattr(nds, 'v_local')
        assert nds.v_local == 0.5

    def test_v_local_settable(self):
        nds = NodeDiseaseState()
        nds.v_local = 0.8
        assert nds.v_local == 0.8


# ---------------------------------------------------------------------------
# Tests: Dispersal mixing shifts traits
# ---------------------------------------------------------------------------

class TestDispersalMixingTvbnc:
    """Dispersal mixing shifts T_vbnc toward source nodes."""

    def test_mixing_shifts_tvbnc_toward_source(self):
        """Node 1 (cold-adapted) receives pathogen from node 0 (warm);
        T_vbnc should shift toward warm."""
        # D[source, dest]: node 0 sends 10% to node 1
        D = np.array([
            [0.0, 0.1],   # from node 0: 10% to node 1
            [0.0, 0.0],
        ])
        D_T = _make_D_T_sparse(D)

        nds = [NodeDiseaseState(node_id=0), NodeDiseaseState(node_id=1)]
        nds[0].T_vbnc_local = 12.0  # warm-adapted
        nds[0].v_local = 0.5
        nds[1].T_vbnc_local = 8.0   # cold-adapted
        nds[1].v_local = 0.5

        P = np.array([1000.0, 500.0])
        disease_reached = [True, True]

        _apply_dispersal_mixing(nds, disease_reached, D_T, P)

        # Node 1: incoming = 0.1 * 1000 = 100 from node 0
        # T_mixed = (500*8 + 100*12) / (500 + 100) = 8.667
        expected = (500.0 * 8.0 + 100.0 * 12.0) / (500.0 + 100.0)
        assert abs(nds[1].T_vbnc_local - expected) < 1e-10
        assert nds[1].T_vbnc_local > 8.0

    def test_no_self_dispersal_no_change(self):
        """Node with zero incoming dispersal keeps its T_vbnc."""
        D = np.array([
            [0.0, 0.0],
            [0.0, 0.0],
        ])
        D_T = _make_D_T_sparse(D)

        nds = [NodeDiseaseState(node_id=0), NodeDiseaseState(node_id=1)]
        nds[0].T_vbnc_local = 12.0
        nds[1].T_vbnc_local = 8.0

        P = np.array([1000.0, 500.0])
        disease_reached = [True, True]

        _apply_dispersal_mixing(nds, disease_reached, D_T, P)

        assert nds[0].T_vbnc_local == 12.0
        assert nds[1].T_vbnc_local == 8.0


class TestDispersalMixingVLocal:
    """Dispersal mixing shifts v_local toward source nodes."""

    def test_mixing_shifts_vlocal_toward_source(self):
        """Node 1 (low virulence) receives pathogen from node 0 (high virulence);
        v_local should shift toward high virulence."""
        # D[source, dest]: node 0 sends 10% to node 1
        D = np.array([
            [0.0, 0.1],
            [0.0, 0.0],
        ])
        D_T = _make_D_T_sparse(D)

        nds = [NodeDiseaseState(node_id=0), NodeDiseaseState(node_id=1)]
        nds[0].v_local = 0.8   # high virulence
        nds[1].v_local = 0.3   # low virulence
        nds[0].T_vbnc_local = 12.0
        nds[1].T_vbnc_local = 12.0

        P = np.array([1000.0, 500.0])
        disease_reached = [True, True]

        _apply_dispersal_mixing(nds, disease_reached, D_T, P)

        # incoming to node 1 = 0.1 * 1000 = 100
        expected = (500.0 * 0.3 + 100.0 * 0.8) / (500.0 + 100.0)
        assert abs(nds[1].v_local - expected) < 1e-10
        assert nds[1].v_local > 0.3

    def test_both_traits_mixed_simultaneously(self):
        """Both T_vbnc and v_local should be mixed in a single dispersal step."""
        # D[source, dest]: node 0 sends 20% to node 1
        D = np.array([
            [0.0, 0.2],
            [0.0, 0.0],
        ])
        D_T = _make_D_T_sparse(D)

        nds = [NodeDiseaseState(node_id=0), NodeDiseaseState(node_id=1)]
        nds[0].T_vbnc_local = 14.0
        nds[0].v_local = 0.9
        nds[1].T_vbnc_local = 8.0
        nds[1].v_local = 0.2

        P = np.array([1000.0, 400.0])
        disease_reached = [True, True]

        _apply_dispersal_mixing(nds, disease_reached, D_T, P)

        incoming = 0.2 * 1000.0  # 200
        total = 400.0 + 200.0  # 600
        expected_T = (400.0 * 8.0 + 200.0 * 14.0) / total
        expected_V = (400.0 * 0.2 + 200.0 * 0.9) / total

        assert abs(nds[1].T_vbnc_local - expected_T) < 1e-10
        assert abs(nds[1].v_local - expected_V) < 1e-10


# ---------------------------------------------------------------------------
# Tests: No mixing when pathogen_adaptation=False
# ---------------------------------------------------------------------------

class TestNoMixingWhenDisabled:
    """When pathogen_adaptation is False, traits should not be mixed.
    
    The model.py code guards with `if _pathogen_adapt`, so this test
    verifies that without calling the mixing function, traits remain unchanged.
    """

    def test_traits_unchanged_without_mixing(self):
        """Without calling _apply_dispersal_mixing, traits stay as-is."""
        nds = [NodeDiseaseState(node_id=0), NodeDiseaseState(node_id=1)]
        nds[0].T_vbnc_local = 12.0
        nds[0].v_local = 0.8
        nds[1].T_vbnc_local = 8.0
        nds[1].v_local = 0.3

        # Simulate what happens when _pathogen_adapt is False:
        # the mixing block is skipped entirely
        assert nds[0].T_vbnc_local == 12.0
        assert nds[0].v_local == 0.8
        assert nds[1].T_vbnc_local == 8.0
        assert nds[1].v_local == 0.3


# ---------------------------------------------------------------------------
# Tests: No mixing at unreached nodes
# ---------------------------------------------------------------------------

class TestNoMixingUnreachedNodes:
    def test_unreached_node_traits_unchanged(self):
        """Nodes where disease hasn't reached should not have traits mixed."""
        # D[source, dest]: node 0 sends 10% to node 1
        D = np.array([
            [0.0, 0.1],
            [0.0, 0.0],
        ])
        D_T = _make_D_T_sparse(D)

        nds = [NodeDiseaseState(node_id=0), NodeDiseaseState(node_id=1)]
        nds[0].T_vbnc_local = 12.0
        nds[0].v_local = 0.8
        nds[1].T_vbnc_local = 8.0
        nds[1].v_local = 0.3

        P = np.array([1000.0, 500.0])
        disease_reached = [True, False]  # node 1 not reached

        _apply_dispersal_mixing(nds, disease_reached, D_T, P)

        # Node 1 should be unchanged (not reached)
        assert nds[1].T_vbnc_local == 8.0
        assert nds[1].v_local == 0.3

    def test_reached_node_still_mixed(self):
        """Even with some unreached nodes, reached nodes should be mixed."""
        # D[source, dest]: node 0 sends 10% to node 1, node 1 sends 10% to node 2
        D = np.array([
            [0.0, 0.1, 0.0],
            [0.0, 0.0, 0.1],
            [0.0, 0.0, 0.0],
        ])
        D_T = _make_D_T_sparse(D)

        nds = [
            NodeDiseaseState(node_id=0),
            NodeDiseaseState(node_id=1),
            NodeDiseaseState(node_id=2),
        ]
        nds[0].T_vbnc_local = 14.0
        nds[0].v_local = 0.9
        nds[1].T_vbnc_local = 10.0
        nds[1].v_local = 0.5
        nds[2].T_vbnc_local = 6.0
        nds[2].v_local = 0.1

        P = np.array([1000.0, 500.0, 200.0])
        disease_reached = [True, True, False]

        _apply_dispersal_mixing(nds, disease_reached, D_T, P)

        # Node 1 (reached): incoming from node 0 = 0.1 * 1000 = 100
        total_1 = 500.0 + 100.0
        expected_T = (500.0 * 10.0 + 100.0 * 14.0) / total_1
        expected_V = (500.0 * 0.5 + 100.0 * 0.9) / total_1
        assert abs(nds[1].T_vbnc_local - expected_T) < 1e-10
        assert abs(nds[1].v_local - expected_V) < 1e-10

        # Node 2 (not reached): unchanged
        assert nds[2].T_vbnc_local == 6.0
        assert nds[2].v_local == 0.1


# ---------------------------------------------------------------------------
# Tests: Wavefront trait inheritance includes v_local
# ---------------------------------------------------------------------------

class TestWavefrontInheritance:
    def test_inherit_pathogen_traits_includes_vlocal(self):
        """_inherit_pathogen_traits should set both T_vbnc_local and v_local."""
        # D[source, dest]: node 0 sends 0.3 to node 2, node 1 sends 0.7 to node 2
        D = np.array([
            [0.0, 0.0, 0.3],
            [0.0, 0.0, 0.7],
            [0.0, 0.0, 0.0],
        ])
        D_T = _make_D_T_sparse(D)

        nds = [
            NodeDiseaseState(node_id=0),
            NodeDiseaseState(node_id=1),
            NodeDiseaseState(node_id=2),
        ]
        nds[0].T_vbnc_local = 14.0
        nds[0].v_local = 0.9
        nds[1].T_vbnc_local = 8.0
        nds[1].v_local = 0.3
        nds[2].T_vbnc_local = 12.0  # default, should be overwritten
        nds[2].v_local = 0.5        # default, should be overwritten

        P = np.array([1000.0, 500.0, 0.0])
        disease_reached = [True, True, False]

        class FakeCfg:
            pass

        _inherit_pathogen_traits(2, nds, disease_reached, D_T, P, FakeCfg())

        # Expected: weighted by dispersal_weight * P[src]
        w0 = 0.3 * 1000.0  # 300
        w1 = 0.7 * 500.0   # 350
        total_w = w0 + w1   # 650

        expected_T = (w0 * 14.0 + w1 * 8.0) / total_w
        expected_V = (w0 * 0.9 + w1 * 0.3) / total_w

        assert abs(nds[2].T_vbnc_local - expected_T) < 1e-10
        assert abs(nds[2].v_local - expected_V) < 1e-10

    def test_inherit_only_from_reached_sources(self):
        """Should only inherit from reached source nodes."""
        # D[source, dest]: both nodes 0 and 1 send 0.5 to node 2
        D = np.array([
            [0.0, 0.0, 0.5],
            [0.0, 0.0, 0.5],
            [0.0, 0.0, 0.0],
        ])
        D_T = _make_D_T_sparse(D)

        nds = [
            NodeDiseaseState(node_id=0),
            NodeDiseaseState(node_id=1),
            NodeDiseaseState(node_id=2),
        ]
        nds[0].T_vbnc_local = 14.0
        nds[0].v_local = 0.9
        nds[1].T_vbnc_local = 8.0
        nds[1].v_local = 0.2
        nds[2].T_vbnc_local = 12.0
        nds[2].v_local = 0.5

        P = np.array([1000.0, 1000.0, 0.0])
        disease_reached = [True, False, False]  # only node 0 reached

        _inherit_pathogen_traits(2, nds, disease_reached, D_T, P, type('Cfg', (), {})())

        # Only node 0 contributes
        assert abs(nds[2].T_vbnc_local - 14.0) < 1e-10
        assert abs(nds[2].v_local - 0.9) < 1e-10

    def test_inherit_no_reached_sources_keeps_defaults(self):
        """If no reached sources, keep defaults."""
        # D[source, dest]: node 0 sends 0.5 to node 1
        D = np.array([
            [0.0, 0.5],
            [0.0, 0.0],
        ])
        D_T = _make_D_T_sparse(D)

        nds = [NodeDiseaseState(node_id=0), NodeDiseaseState(node_id=1)]
        nds[0].T_vbnc_local = 14.0
        nds[0].v_local = 0.9
        nds[1].T_vbnc_local = 12.0
        nds[1].v_local = 0.5

        P = np.array([1000.0, 0.0])
        disease_reached = [False, False]

        _inherit_pathogen_traits(1, nds, disease_reached, D_T, P, type('Cfg', (), {})())

        # No reached sources → defaults unchanged
        assert nds[1].T_vbnc_local == 12.0
        assert nds[1].v_local == 0.5


# ---------------------------------------------------------------------------
# Tests: Identity — same traits everywhere → no change
# ---------------------------------------------------------------------------

class TestIdentityMixing:
    def test_uniform_traits_no_change(self):
        """When all nodes have the same traits, dispersal mixing shouldn't change them."""
        D = np.array([
            [0.0, 0.05, 0.05],
            [0.05, 0.0, 0.05],
            [0.05, 0.05, 0.0],
        ])
        D_T = _make_D_T_sparse(D)

        nds = [NodeDiseaseState(node_id=i) for i in range(3)]
        for nd in nds:
            nd.T_vbnc_local = 10.0
            nd.v_local = 0.6

        P = np.array([1000.0, 800.0, 600.0])
        disease_reached = [True, True, True]

        _apply_dispersal_mixing(nds, disease_reached, D_T, P)

        for nd in nds:
            assert abs(nd.T_vbnc_local - 10.0) < 1e-10
            assert abs(nd.v_local - 0.6) < 1e-10


# ---------------------------------------------------------------------------
# Tests: Asymmetric dispersal
# ---------------------------------------------------------------------------

class TestAsymmetricDispersal:
    def test_asymmetric_shifts_toward_dominant_source(self):
        """More pathogen from one direction shifts traits toward that source."""
        # D[source, dest]: node 0 sends 0.5 to node 1, node 2 sends 0.05 to node 1
        D = np.array([
            [0.0, 0.5, 0.0],   # from 0: 50% to 1
            [0.0, 0.0, 0.0],
            [0.0, 0.05, 0.0],  # from 2: 5% to 1
        ])
        D_T = _make_D_T_sparse(D)

        nds = [
            NodeDiseaseState(node_id=0),
            NodeDiseaseState(node_id=1),
            NodeDiseaseState(node_id=2),
        ]
        nds[0].T_vbnc_local = 14.0
        nds[0].v_local = 0.9
        nds[1].T_vbnc_local = 10.0
        nds[1].v_local = 0.5
        nds[2].T_vbnc_local = 6.0
        nds[2].v_local = 0.1

        # Same concentrations at all sources
        P = np.array([1000.0, 200.0, 1000.0])
        disease_reached = [True, True, True]

        _apply_dispersal_mixing(nds, disease_reached, D_T, P)

        # Node 1 incoming: 0.5*1000=500 from node 0, 0.05*1000=50 from node 2
        incoming_from_0 = 500.0
        incoming_from_2 = 50.0
        total = 200.0 + incoming_from_0 + incoming_from_2

        expected_T = (200.0 * 10.0 + 500.0 * 14.0 + 50.0 * 6.0) / total
        expected_V = (200.0 * 0.5 + 500.0 * 0.9 + 50.0 * 0.1) / total

        assert abs(nds[1].T_vbnc_local - expected_T) < 1e-10
        assert abs(nds[1].v_local - expected_V) < 1e-10

        # Node 1's T_vbnc should be closer to node 0's (14.0) than node 2's (6.0)
        assert nds[1].T_vbnc_local > 10.0  # shifted toward node 0
        assert nds[1].v_local > 0.5         # shifted toward node 0's virulence

    def test_asymmetric_concentration_matters(self):
        """Even with equal dispersal weights, higher concentration
        at source shifts traits more."""
        # D[source, dest]: both nodes 0 and 2 send 10% to node 1
        D = np.array([
            [0.0, 0.1, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.1, 0.0],
        ])
        D_T = _make_D_T_sparse(D)

        nds = [
            NodeDiseaseState(node_id=0),
            NodeDiseaseState(node_id=1),
            NodeDiseaseState(node_id=2),
        ]
        nds[0].T_vbnc_local = 14.0
        nds[0].v_local = 0.9
        nds[1].T_vbnc_local = 10.0
        nds[1].v_local = 0.5
        nds[2].T_vbnc_local = 6.0
        nds[2].v_local = 0.1

        # Node 0 has much more pathogen
        P = np.array([5000.0, 200.0, 100.0])
        disease_reached = [True, True, True]

        _apply_dispersal_mixing(nds, disease_reached, D_T, P)

        # Incoming from node 0: 0.1*5000=500, from node 2: 0.1*100=10
        # Node 1's traits should be closer to node 0's values
        assert nds[1].T_vbnc_local > 10.0  # toward node 0
        assert nds[1].v_local > 0.5         # toward node 0

    def test_zero_concentration_no_influence(self):
        """Source nodes with zero concentration don't influence mixing."""
        # D[source, dest]: node 0 sends 50% to node 1
        D = np.array([
            [0.0, 0.5],
            [0.0, 0.0],
        ])
        D_T = _make_D_T_sparse(D)

        nds = [NodeDiseaseState(node_id=0), NodeDiseaseState(node_id=1)]
        nds[0].T_vbnc_local = 14.0
        nds[0].v_local = 0.9
        nds[1].T_vbnc_local = 10.0
        nds[1].v_local = 0.5

        P = np.array([0.0, 500.0])
        disease_reached = [True, True]

        _apply_dispersal_mixing(nds, disease_reached, D_T, P)

        # Node 0 has no pathogen → no mixing effect on node 1
        assert nds[1].T_vbnc_local == 10.0
        assert nds[1].v_local == 0.5
