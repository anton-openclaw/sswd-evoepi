#!/usr/bin/env python3
"""Minimal distance method comparison - matrices only."""

import sys
from pathlib import Path
import numpy as np

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from sswd_evoepi.spatial import make_5node_network, get_5node_definitions

def compare_matrices():
    print("=" * 60)
    print("DISTANCE METHOD COMPARISON - MATRICES ONLY")
    print("=" * 60)
    
    # Create networks
    print("Creating networks...")
    net_haversine = make_5node_network(seed=42, use_overwater=False)
    net_overwater = make_5node_network(seed=42, use_overwater=True)
    
    print("Networks created successfully!")
    print()
    
    # Get node names
    node_names = [nd.name for nd in get_5node_definitions()]
    
    # Distance matrices
    print("─── Distance Matrices (km) ───")
    print("Haversine × 1.5:")
    print_matrix(net_haversine.distances, node_names)
    print("\nOverwater:")
    print_matrix(net_overwater.distances, node_names)
    print("\nDifference (Overwater - Haversine):")
    diff = net_overwater.distances - net_haversine.distances
    print_matrix(diff, node_names)
    print(f"Max absolute difference: {np.max(np.abs(diff)):.1f} km")
    print()
    
    # Larval connectivity matrices
    print("─── Larval Connectivity Matrices (C_ij) ───")
    print("Haversine × 1.5:")
    print_matrix(net_haversine.C, node_names, format=".6f")
    print("\nOverwater:")
    print_matrix(net_overwater.C, node_names, format=".6f")
    print("\nDifference (Overwater - Haversine):")
    diff_c = net_overwater.C - net_haversine.C
    print_matrix(diff_c, node_names, format="+.6f")
    print()
    
    # Pathogen dispersal matrices
    print("─── Pathogen Dispersal Matrices (D_ij) ───")
    print("Haversine × 1.5:")
    print_matrix(net_haversine.D, node_names, format=".6f")
    print("\nOverwater:")
    print_matrix(net_overwater.D, node_names, format=".6f")
    print("\nDifference (Overwater - Haversine):")
    diff_d = net_overwater.D - net_haversine.D
    print_matrix(diff_d, node_names, format="+.6f")
    print()
    
    # Analysis
    print("─── Key Findings ───")
    print()
    
    # Self-recruitment changes
    print("Self-recruitment changes (diagonal of C matrix):")
    for i, name in enumerate(node_names):
        h_self = net_haversine.C[i, i]
        o_self = net_overwater.C[i, i]
        change = (o_self - h_self) / h_self * 100 if h_self > 0 else 0
        print(f"  {name:<25s}: {h_self:.6f} → {o_self:.6f} ({change:+5.1f}%)")
    print()
    
    # Connectivity strength changes
    print("Strongest connectivity changes:")
    flat_diff = diff_c.flatten()
    flat_idx = np.argsort(np.abs(flat_diff))[::-1]  # Sort by absolute change
    
    for idx in flat_idx[:10]:  # Top 10 changes
        i, j = divmod(idx, len(node_names))
        if i != j:  # Skip diagonal
            h_val = net_haversine.C[i, j]
            o_val = net_overwater.C[i, j]
            change = diff_c[i, j]
            pct_change = change / h_val * 100 if h_val > 0 else float('inf')
            print(f"  {node_names[i]:<12s} → {node_names[j]:<12s}: "
                  f"{h_val:.6f} → {o_val:.6f} ({change:+.6f}, {pct_change:+5.1f}%)")
    print()
    
    # Distance reduction analysis
    print("Distance reductions (negative = shorter overwater distance):")
    for i in range(len(node_names)):
        for j in range(i+1, len(node_names)):
            h_dist = net_haversine.distances[i, j]
            o_dist = net_overwater.distances[i, j]
            diff_dist = o_dist - h_dist
            pct_diff = diff_dist / h_dist * 100
            print(f"  {node_names[i]:<12s} → {node_names[j]:<12s}: "
                  f"{diff_dist:+6.0f} km ({pct_diff:+5.1f}%)")
    
    return net_haversine, net_overwater

def print_matrix(matrix, labels, format=".1f"):
    """Print a matrix with row/column labels."""
    n = len(labels)
    
    # Header row
    print("    " + "".join([f"{label[:8]:>10s}" for label in labels]))
    print("    " + "─" * (10 * n))
    
    # Data rows
    for i, row_label in enumerate(labels):
        row_str = f"{row_label[:8]:<4s}"
        for j in range(n):
            val = matrix[i, j]
            if format == ".6f":
                row_str += f"{val:>10.6f}"
            elif format == "+.6f":
                row_str += f"{val:>+10.6f}"
            else:  # Default to .1f
                row_str += f"{val:>10.1f}"
        print(row_str)

if __name__ == "__main__":
    compare_matrices()