#!/usr/bin/env python3
"""Compare 5-node distances: Haversine vs Overwater."""

import numpy as np
from sswd_evoepi.spatial import get_5node_definitions, compute_distance_matrix, load_overwater_distances

def main():
    nodes = get_5node_definitions()
    node_names = [n.name for n in nodes]

    # Haversine distances
    lats = np.array([n.lat for n in nodes])
    lons = np.array([n.lon for n in nodes])
    haversine_dist = compute_distance_matrix(lats, lons, tortuosity=1.5)

    # Overwater distances
    overwater_dist = load_overwater_distances(nodes)

    print('5-Node Distance Comparison (km):')
    print('=' * 70)
    print(f'{"Pair":<30} {"Haversine*1.5":<15} {"Overwater":<15} {"Diff %":>10}')
    print('-' * 70)

    for i in range(5):
        for j in range(i+1, 5):
            h_dist = haversine_dist[i, j]
            o_dist = overwater_dist[i, j]
            diff_pct = ((o_dist - h_dist) / h_dist) * 100
            pair_name = f'{node_names[i][:12]} <-> {node_names[j][:12]}'
            print(f'{pair_name:<30} {h_dist:>12.0f} {o_dist:>12.0f} {diff_pct:>9.1f}%')

    print('=' * 70)
    print('Key: Positive % = overwater longer, Negative % = overwater shorter')
    print('     Overwater distances are more accurate (real coastal routing)')

if __name__ == '__main__':
    main()