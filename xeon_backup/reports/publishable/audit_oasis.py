#!/usr/bin/env python3
"""Audit numerical claims in sec4b_site_dynamics.tex against raw NPZ data."""

import numpy as np
import json

seeds = [42, 123, 999]
base = "/home/starbot/.openclaw/workspace/sswd-evoepi/results/calibration/F01"

def gini(x):
    """Compute Gini coefficient for array x."""
    x = np.array(x, dtype=float)
    if x.sum() == 0:
        return 0.0
    x = np.sort(x)
    n = len(x)
    index = np.arange(1, n + 1)
    return (2 * np.sum(index * x) - (n + 1) * np.sum(x)) / (n * np.sum(x))

results = {}

for seed in seeds:
    path = f"{base}/monthly_seed{seed}.npz"
    data = np.load(path, allow_pickle=True)
    
    pop = data['populations']       # [463, 896]
    names = data['site_names']      # [896]
    lats = data['site_lats']        # [896]
    
    # Convert names to strings if needed
    names = np.array([str(n) for n in names])
    
    print(f"\n=== Seed {seed} ===")
    print(f"  populations shape: {pop.shape}")
    print(f"  site_names shape: {names.shape}")
    print(f"  site_lats shape: {lats.shape}")
    
    # Region masks
    or_mask = np.array([n.startswith("OR-") for n in names])
    pws_mask = np.array([n.startswith("AK-PWS-") for n in names])
    fn_mask = np.array([n.startswith("AK-FN-") for n in names])
    
    n_or = or_mask.sum()
    n_pws = pws_mask.sum()
    n_fn = fn_mask.sum()
    
    print(f"  Oregon sites: {n_or}")
    print(f"  AK-PWS sites: {n_pws}")
    print(f"  AK-FN sites: {n_fn}")
    
    # Final timestep populations
    final_pop = pop[-1, :]
    
    # === Oregon ===
    or_final = final_pop[or_mask]
    or_names = names[or_mask]
    or_lats_r = lats[or_mask]
    
    or_alive = (or_final > 0).sum()
    or_total = or_final.sum()
    
    # Top 3 concentration
    or_sorted_idx = np.argsort(or_final)[::-1]
    or_top3_pop = or_final[or_sorted_idx[:3]].sum()
    or_top3_frac = or_top3_pop / or_total if or_total > 0 else 0
    
    # Top 5 sites
    or_top5_idx = or_sorted_idx[:5]
    or_top5_names = or_names[or_top5_idx]
    or_top5_lats = or_lats_r[or_top5_idx]
    or_top5_pops = or_final[or_top5_idx]
    
    # Gini
    or_gini = gini(or_final)
    
    print(f"\n  Oregon final timestep:")
    print(f"    Alive sites: {or_alive} of {n_or}")
    print(f"    Total pop: {or_total:.0f}")
    print(f"    Top 3 concentration: {or_top3_frac*100:.1f}%")
    print(f"    Gini: {or_gini:.4f}")
    print(f"    Top 5 sites:")
    for i in range(min(5, len(or_top5_names))):
        print(f"      {or_top5_names[i]}: pop={or_top5_pops[i]:.0f}, lat={or_top5_lats[i]:.4f}")
    
    # === AK-PWS ===
    pws_final = final_pop[pws_mask]
    pws_names = names[pws_mask]
    
    pws_alive = (pws_final > 0).sum()
    pws_total = pws_final.sum()
    
    pws_sorted_idx = np.argsort(pws_final)[::-1]
    pws_top3_pop = pws_final[pws_sorted_idx[:3]].sum()
    pws_top3_frac = pws_top3_pop / pws_total if pws_total > 0 else 0
    
    pws_gini = gini(pws_final)
    
    print(f"\n  AK-PWS final timestep:")
    print(f"    Alive sites: {pws_alive} of {n_pws}")
    print(f"    Total pop: {pws_total:.0f}")
    print(f"    Top 3 concentration: {pws_top3_frac*100:.1f}%")
    print(f"    Gini: {pws_gini:.4f}")
    
    # === AK-FN ===
    fn_final = final_pop[fn_mask]
    fn_names = names[fn_mask]
    
    fn_alive = (fn_final > 0).sum()
    fn_total = fn_final.sum()
    
    fn_sorted_idx = np.argsort(fn_final)[::-1]
    fn_top3_pop = fn_final[fn_sorted_idx[:3]].sum()
    fn_top3_frac = fn_top3_pop / fn_total if fn_total > 0 else 0
    
    fn_gini = gini(fn_final)
    
    print(f"\n  AK-FN final timestep:")
    print(f"    Alive sites: {fn_alive} of {n_fn}")
    print(f"    Total pop: {fn_total:.0f}")
    print(f"    Top 3 concentration: {fn_top3_frac*100:.1f}%")
    print(f"    Gini: {fn_gini:.4f}")
    
    results[seed] = {
        'or_alive': int(or_alive),
        'or_total': int(n_or),
        'or_top3_frac': float(or_top3_frac),
        'or_gini': float(or_gini),
        'or_top5_names': list(or_top5_names),
        'or_top5_lats': [float(x) for x in or_top5_lats],
        'or_top5_pops': [float(x) for x in or_top5_pops],
        'pws_alive': int(pws_alive),
        'pws_total': int(n_pws),
        'pws_top3_frac': float(pws_top3_frac),
        'pws_gini': float(pws_gini),
        'fn_alive': int(fn_alive),
        'fn_total': int(n_fn),
        'fn_top3_frac': float(fn_top3_frac),
        'fn_gini': float(fn_gini),
    }

# Cross-seed analysis: top 5 Oregon sites
print("\n\n=== Cross-seed Oregon Top 5 Analysis ===")
all_top5 = set()
top5_per_seed = {}
for seed in seeds:
    top5_per_seed[seed] = set(results[seed]['or_top5_names'])
    all_top5 |= top5_per_seed[seed]
    print(f"  Seed {seed}: {results[seed]['or_top5_names']}")

print(f"\n  Total unique sites in any top-5: {len(all_top5)}")
print(f"  Sites: {sorted(all_top5)}")

# Count how many seeds each site appears in
from collections import Counter
site_counts = Counter()
for seed in seeds:
    for s in top5_per_seed[seed]:
        site_counts[s] += 1

multi_seed = {s: c for s, c in site_counts.items() if c > 1}
print(f"  Sites in >1 seed's top 5: {multi_seed}")
print(f"  Count of sites in >1 seed: {len(multi_seed)}")

# Check specific sites: OR-048 and OR-020
print("\n\n=== Specific Site Verification ===")
for seed in seeds:
    path = f"{base}/monthly_seed{seed}.npz"
    data = np.load(path, allow_pickle=True)
    names = np.array([str(n) for n in data['site_names']])
    lats = data['site_lats']
    final_pop = data['populations'][-1, :]
    or_mask = np.array([n.startswith("OR-") for n in names])
    or_total = final_pop[or_mask].sum()
    
    for target in ['OR-048', 'OR-020']:
        idx = np.where(names == target)[0]
        if len(idx) > 0:
            i = idx[0]
            frac = final_pop[i] / or_total * 100 if or_total > 0 else 0
            print(f"  Seed {seed}: {target} lat={lats[i]:.4f}, pop={final_pop[i]:.0f}, frac={frac:.1f}%")
        else:
            print(f"  Seed {seed}: {target} NOT FOUND")

# Check latitude range of top sites
print("\n\n=== Southern Oregon Lat Range Check ===")
for seed in seeds:
    path = f"{base}/monthly_seed{seed}.npz"
    data = np.load(path, allow_pickle=True)
    names = np.array([str(n) for n in data['site_names']])
    lats = data['site_lats']
    final_pop = data['populations'][-1, :]
    or_mask = np.array([n.startswith("OR-") for n in names])
    or_final = final_pop[or_mask]
    or_names = names[or_mask]
    or_lats_r = lats[or_mask]
    
    or_sorted_idx = np.argsort(or_final)[::-1]
    print(f"  Seed {seed} top 10:")
    for i in range(min(10, len(or_sorted_idx))):
        idx = or_sorted_idx[i]
        print(f"    {or_names[idx]}: pop={or_final[idx]:.0f}, lat={or_lats_r[idx]:.4f}")

# Initial timestep check (should be uniform at K=5000)
print("\n\n=== Initial Timestep Check ===")
data = np.load(f"{base}/monthly_seed42.npz", allow_pickle=True)
pop = data['populations']
names = np.array([str(n) for n in data['site_names']])
or_mask = np.array([n.startswith("OR-") for n in names])
or_init = pop[0, or_mask]
print(f"  Oregon t=0: min={or_init.min():.0f}, max={or_init.max():.0f}, mean={or_init.mean():.0f}")
print(f"  Gini at t=0: {gini(or_init):.4f}")

# Summary for audit
print("\n\n========== AUDIT SUMMARY ==========")
print("\n1. Oregon site counts: claimed 17-21 of 59")
for seed in seeds:
    r = results[seed]
    print(f"   Seed {seed}: {r['or_alive']} of {r['or_total']}")

print("\n2. Oregon top 3 concentration: claimed 39-79%")
for seed in seeds:
    r = results[seed]
    print(f"   Seed {seed}: {r['or_top3_frac']*100:.1f}%")

print("\n3. Oregon Gini: claimed 0.85-0.92")
for seed in seeds:
    r = results[seed]
    print(f"   Seed {seed}: {r['or_gini']:.4f}")

print("\n4. AK-PWS site counts: claimed 37-38 of 43")
for seed in seeds:
    r = results[seed]
    print(f"   Seed {seed}: {r['pws_alive']} of {r['pws_total']}")

print("\n5. AK-PWS top 3 concentration: claimed 25-29%")
for seed in seeds:
    r = results[seed]
    print(f"   Seed {seed}: {r['pws_top3_frac']*100:.1f}%")

print("\n6. AK-PWS Gini: claimed 0.65-0.71")
for seed in seeds:
    r = results[seed]
    print(f"   Seed {seed}: {r['pws_gini']:.4f}")

print("\n7. AK-FN site counts: claimed 37-40 of 40")
for seed in seeds:
    r = results[seed]
    print(f"   Seed {seed}: {r['fn_alive']} of {r['fn_total']}")

print("\n8. Cross-seed top-5 unique: claimed 13 unique, 2 in >1 seed")
print(f"   Unique: {len(all_top5)}")
print(f"   In >1 seed: {len(multi_seed)} -> {multi_seed}")
