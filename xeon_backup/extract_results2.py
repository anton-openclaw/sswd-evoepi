import json, os

base = os.path.expanduser("~/projects/sswd-evoepi/results/calibration/W289-W308")
results = []
for i in range(289, 309):
    path = os.path.join(base, f"W{i}", "result_seed137.json")
    if os.path.isfile(path):
        with open(path) as f:
            d = json.load(f)
        # Try multiple keys for RMSLE
        rmsle = d.get("rmsle", d.get("rmse_log", d.get("RMSLE", None)))
        rec = d.get("regional_recovery", d.get("recovery", {}))
        params = d.get("params", d.get("config", {}))
        
        # Debug: print all top-level keys for first result
        if i == 289:
            print(f"DEBUG W289 keys: {list(d.keys())}")
            print(f"DEBUG W289 snippet: { {k: str(v)[:80] for k,v in d.items()} }")
        
        results.append({
            "run": f"W{i}",
            "rmsle": rmsle,
            "rec": rec,
            "params": params,
            "raw": d
        })

# Filter out None RMSLE
valid = [r for r in results if r["rmsle"] is not None]
invalid = [r for r in results if r["rmsle"] is None]
print(f"\nValid RMSLE: {len(valid)}, Missing RMSLE: {len(invalid)}")
if invalid:
    for r in invalid:
        print(f"  {r['run']} missing RMSLE. Keys: {list(r['raw'].keys())}")

valid.sort(key=lambda x: x["rmsle"])

print("\n=== TOP 5 BY RMSLE ===")
for r in valid[:5]:
    run = r["run"]
    rmsle = r["rmsle"]
    rec = r["rec"]
    print(f"\n--- {run} (RMSLE={rmsle:.4f}) ---")
    for k in ["AK-PWS", "AK-FN", "BC-N", "SS-S", "JDF", "OR", "CA-N"]:
        v = rec.get(k, None)
        if v is not None:
            print(f"  {k}: {v:.1%}", end="")
        else:
            print(f"  {k}: N/A", end="")
    print()
    print(f"  Params: {json.dumps(r['params'], default=str)}")

print("\n=== FULL RANKING ===")
for r in valid:
    run = r["run"]
    rmsle = r["rmsle"]
    rec = r["rec"]
    ak = rec.get("AK-PWS")
    orv = rec.get("OR")
    can = rec.get("CA-N")
    ak_s = f"{ak:.1%}" if ak is not None else "N/A"
    or_s = f"{orv:.1%}" if orv is not None else "N/A"
    ca_s = f"{can:.1%}" if can is not None else "N/A"
    print(f"{run}: RMSLE={rmsle:.4f}  AK-PWS={ak_s}  OR={or_s}  CA-N={ca_s}")
