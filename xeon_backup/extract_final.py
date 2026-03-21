import json, os, ast

base = os.path.expanduser("~/projects/sswd-evoepi/results/calibration/W289-W308")
results = []

for i in range(289, 309):
    path = os.path.join(base, f"W{i}", "result_seed137.json")
    with open(path) as f:
        d = json.load(f)
    
    scoring = d["scoring"]
    if isinstance(scoring, str):
        scoring = ast.literal_eval(scoring)
    
    rr = d["region_recovery"]
    if isinstance(rr, str):
        rr = ast.literal_eval(rr)
    
    rmsle = scoring["rmsle"]
    
    # Try to get params from config file
    cfg_path = os.path.join(base, f"W{i}", "config.json")
    params = {}
    if os.path.isfile(cfg_path):
        with open(cfg_path) as f:
            params = json.load(f)
    else:
        # Try config.yaml
        cfg_path = os.path.join(base, f"W{i}", "config.yaml")
        if os.path.isfile(cfg_path):
            params = {"file": cfg_path}
    
    wall = d.get("wall_time_seconds", 0)
    
    results.append({
        "run": f"W{i}",
        "rmsle": rmsle,
        "wall_hrs": wall / 3600,
        "AK-PWS": rr.get("AK-PWS", 0),
        "AK-FN": rr.get("AK-FN", 0),
        "AK-FS": rr.get("AK-FS", 0),
        "BC-N": rr.get("BC-N", 0),
        "SS-S": rr.get("SS-S", 0),
        "JDF": rr.get("JDF", 0),
        "OR": rr.get("OR", 0),
        "CA-N": rr.get("CA-N", 0),
        "params": params
    })

results.sort(key=lambda x: x["rmsle"])

print("=== TOP 5 BY RMSLE ===\n")
for r in results[:5]:
    print(f"{r['run']}: RMSLE={r['rmsle']:.4f} ({r['wall_hrs']:.1f}h)")
    print(f"  AK-PWS={r['AK-PWS']:.1%}  AK-FN={r['AK-FN']:.1%}  AK-FS={r['AK-FS']:.1%}")
    print(f"  BC-N={r['BC-N']:.1%}  SS-S={r['SS-S']:.1%}  JDF={r['JDF']:.1%}")
    print(f"  OR={r['OR']:.1%}  CA-N={r['CA-N']:.1%}")
    if r["params"]:
        print(f"  Params: {json.dumps(r['params'], default=str)}")
    print()

print("=== FULL RANKING ===")
for r in results:
    print(f"{r['run']}: RMSLE={r['rmsle']:.4f}  AK-PWS={r['AK-PWS']:.1%}  OR={r['OR']:.1%}  CA-N={r['CA-N']:.1%}")

# Compare to W285 baseline
print(f"\n=== W285 BASELINE: RMSLE=0.666  AK-PWS=19.5%  OR=2.7%  CA-N=2.1% ===")
better = [r for r in results if r["rmsle"] < 0.666]
print(f"Runs beating W285: {len(better)}")
