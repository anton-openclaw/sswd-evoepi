"""Tests for sswd_evoepi.rng — seeded RNG hierarchy and checkpointing."""

import numpy as np
import pytest

from sswd_evoepi.rng import (
    create_rng_hierarchy,
    get_node_rng,
    restore_rng_state,
    rng_state_snapshot,
)


class TestCreateRngHierarchy:
    def test_returns_correct_keys(self):
        rngs = create_rng_hierarchy(42, n_nodes=5)
        assert 'global' in rngs
        assert 'larval' in rngs
        assert 'conservation' in rngs
        for i in range(5):
            assert f'node_{i}' in rngs
        assert len(rngs) == 5 + 3  # 5 nodes + 3 global streams

    def test_generators_are_independent(self):
        """Different streams produce different sequences."""
        rngs = create_rng_hierarchy(42, n_nodes=3)
        vals = {name: rng.random() for name, rng in rngs.items()}
        # All values should be unique (independent streams)
        unique_vals = set(vals.values())
        assert len(unique_vals) == len(vals), "RNG streams produced duplicate values"

    def test_reproducibility(self):
        """Same seed produces identical sequences."""
        rngs1 = create_rng_hierarchy(42, n_nodes=5)
        rngs2 = create_rng_hierarchy(42, n_nodes=5)

        for name in rngs1:
            v1 = rngs1[name].random(100)
            v2 = rngs2[name].random(100)
            np.testing.assert_array_equal(v1, v2)

    def test_different_seeds_differ(self):
        """Different seeds produce different sequences."""
        rngs1 = create_rng_hierarchy(42, n_nodes=3)
        rngs2 = create_rng_hierarchy(43, n_nodes=3)

        # At least one stream should differ (overwhelmingly likely)
        all_same = True
        for name in rngs1:
            v1 = rngs1[name].random(10)
            v2 = rngs2[name].random(10)
            if not np.array_equal(v1, v2):
                all_same = False
                break
        assert not all_same

    def test_zero_nodes(self):
        """Edge case: zero nodes is valid (only global streams)."""
        rngs = create_rng_hierarchy(42, n_nodes=0)
        assert len(rngs) == 3  # global, larval, conservation
        assert 'global' in rngs

    def test_node_independence(self):
        """Adding nodes doesn't change existing nodes' streams.

        This is the key SeedSequence property: spawning is deterministic
        based on position in the sequence, not total count.
        """
        rngs_5 = create_rng_hierarchy(42, n_nodes=5)
        rngs_10 = create_rng_hierarchy(42, n_nodes=10)

        # First 5 nodes should produce identical sequences
        for i in range(5):
            v5 = rngs_5[f'node_{i}'].random(50)
            v10 = rngs_10[f'node_{i}'].random(50)
            np.testing.assert_array_equal(v5, v10)


class TestGetNodeRng:
    def test_valid_node(self):
        rngs = create_rng_hierarchy(42, n_nodes=5)
        rng = get_node_rng(rngs, 3)
        assert isinstance(rng, np.random.Generator)
        # Should be the same object
        assert rng is rngs['node_3']

    def test_invalid_node(self):
        rngs = create_rng_hierarchy(42, n_nodes=5)
        with pytest.raises(KeyError):
            get_node_rng(rngs, 99)


class TestCheckpointing:
    def test_snapshot_and_restore(self):
        """Full round-trip: generate → snapshot → advance → restore → compare."""
        rngs = create_rng_hierarchy(42, n_nodes=3)

        # Advance each stream a bit
        for rng in rngs.values():
            rng.random(10)

        # Snapshot the state
        snapshot = rng_state_snapshot(rngs)

        # Record values AFTER snapshot point
        expected = {}
        for name, rng in rngs.items():
            expected[name] = rng.random(20)

        # Create fresh hierarchy and advance to same point
        rngs2 = create_rng_hierarchy(42, n_nodes=3)
        for rng in rngs2.values():
            rng.random(10)

        # Restore from snapshot
        restore_rng_state(rngs2, snapshot)

        # Values after restore should match
        for name, rng in rngs2.items():
            actual = rng.random(20)
            np.testing.assert_array_equal(actual, expected[name])

    def test_snapshot_keys(self):
        rngs = create_rng_hierarchy(42, n_nodes=2)
        snapshot = rng_state_snapshot(rngs)
        assert set(snapshot.keys()) == set(rngs.keys())

    def test_restore_unknown_stream_raises(self):
        rngs = create_rng_hierarchy(42, n_nodes=2)
        bad_states = {'nonexistent_stream': {}}
        with pytest.raises(KeyError, match="nonexistent_stream"):
            restore_rng_state(rngs, bad_states)


class TestPCG64Properties:
    def test_generator_type(self):
        """All generators use PCG64 bit generator."""
        rngs = create_rng_hierarchy(42, n_nodes=3)
        for rng in rngs.values():
            assert isinstance(rng.bit_generator, np.random.PCG64)

    def test_large_sample_uniformity(self):
        """Basic statistical check: uniform samples should have mean ~0.5."""
        rngs = create_rng_hierarchy(42, n_nodes=1)
        samples = rngs['node_0'].random(100_000)
        assert abs(samples.mean() - 0.5) < 0.01
        assert abs(samples.std() - (1.0 / 12**0.5)) < 0.01
