"""Seeded RNG factory for reproducible simulations.

Uses NumPy's SeedSequence → PCG64 hierarchy to guarantee:
  - Statistical independence between per-node streams
  - Bit-exact replay with the same master seed
  - Adding/removing nodes doesn't affect other nodes' streams

References:
  - integration-architecture-spec.md §5.3
  - NumPy docs: numpy.random.SeedSequence
"""

from __future__ import annotations

from typing import Dict

import numpy as np


def create_rng_hierarchy(
    master_seed: int,
    n_nodes: int,
) -> Dict[str, np.random.Generator]:
    """Create independent RNG streams for each node + global operations.

    Uses SeedSequence spawning to guarantee statistical independence
    between streams (no overlap in 2^128 period PCG64).

    Streams created:
      - 'global':       Global operations (initialization, etc.)
      - 'larval':       Larval dispersal and settlement
      - 'conservation': Conservation module (captive breeding, releases)
      - 'node_0' .. 'node_{n-1}': Per-node streams for disease, reproduction, mortality

    Args:
        master_seed: Master RNG seed (non-negative integer).
        n_nodes: Number of spatial nodes.

    Returns:
        Dictionary mapping stream names to numpy Generator instances.

    Example:
        >>> rngs = create_rng_hierarchy(42, n_nodes=15)
        >>> rngs['global'].random()  # reproducible
        >>> rngs['node_0'].integers(0, 100)
    """
    ss = np.random.SeedSequence(master_seed)
    # Spawn n_nodes + 3 child seeds: global, larval, conservation, then per-node
    child_seeds = ss.spawn(n_nodes + 3)

    rngs: Dict[str, np.random.Generator] = {
        'global': np.random.Generator(np.random.PCG64(child_seeds[0])),
        'larval': np.random.Generator(np.random.PCG64(child_seeds[1])),
        'conservation': np.random.Generator(np.random.PCG64(child_seeds[2])),
    }
    for i in range(n_nodes):
        rngs[f'node_{i}'] = np.random.Generator(
            np.random.PCG64(child_seeds[3 + i])
        )

    return rngs


def get_node_rng(
    rngs: Dict[str, np.random.Generator],
    node_id: int,
) -> np.random.Generator:
    """Get the RNG stream for a specific node.

    Args:
        rngs: RNG hierarchy from create_rng_hierarchy().
        node_id: Node index (0-based).

    Returns:
        Generator for the specified node.

    Raises:
        KeyError: If node_id doesn't have a stream.
    """
    key = f'node_{node_id}'
    if key not in rngs:
        raise KeyError(
            f"No RNG stream for node {node_id}. "
            f"Available nodes: 0–{max(int(k.split('_')[1]) for k in rngs if k.startswith('node_'))}"
        )
    return rngs[key]


def rng_state_snapshot(
    rngs: Dict[str, np.random.Generator],
) -> Dict[str, dict]:
    """Capture full RNG state for checkpointing.

    Returns a dict of {name: state_dict} that can be serialized (e.g. via pickle)
    and restored to resume a simulation exactly.

    Args:
        rngs: RNG hierarchy.

    Returns:
        Dictionary mapping stream names to their internal state dicts.
    """
    return {name: rng.bit_generator.state for name, rng in rngs.items()}


def restore_rng_state(
    rngs: Dict[str, np.random.Generator],
    states: Dict[str, dict],
) -> None:
    """Restore RNG state from a checkpoint snapshot.

    Args:
        rngs: RNG hierarchy (must have same keys as states).
        states: State snapshot from rng_state_snapshot().

    Raises:
        KeyError: If a stream in states doesn't exist in rngs.
    """
    for name, state in states.items():
        if name not in rngs:
            raise KeyError(f"Cannot restore RNG state for unknown stream '{name}'")
        rngs[name].bit_generator.state = state
