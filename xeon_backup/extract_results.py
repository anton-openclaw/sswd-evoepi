import json, os

base = os.path.expanduser("~/projects/sswd-evoepi/results/calibration/W289-W308")
results = []
for i in range(289, 309):
    path = os.path.join(base, f"W{i}", "result_seed137.json")
    if os.path.isfile(path):
        with open(path) as f:
            d = json.load(f)
        rmsle = d.get("rmsle", d.get("rmse_log", None))
        rec = d.get("regional_recovery", d.get("recovery", {}))
        ak_pws = rec.get("AK-PWS", None)
        ak_fn = rec.get("AK-FN", None)
        or_val = rec.get("OR", None)
        ca_n = rec.get("CA-N", None)
        bc_n = rec.get("BC-N", None)
        ss_s = rec.get("SS-S", None)
        jdf = rec.get("JDF", None)
        params = d.get("params", d.get("config", {}))
        results.append({
            "run": f"W{i}",
            "rmsle": rmsle,
            "AK-PWS": ak_pws,
            "AK-FN": ak_fn,
            "OR": or_val,
            "CA-N": ca_n,
            "BC-N": bc_n,
            "SS-S": ss_s,
            "JDF": jdf,
            "params": params
        })

results.sort(key=lambda x: x["rmsle"] if x["rmsle"] is not None else 999)

print("=== TOP 5 BY RMSLE ===")
for r in results[:5]:
    run = r["run"]
    rmsle = r["rmsle"]
    print(f"--- {run} (RMSLE={rmsle:.4f}) ---")
    vals = []
    for k in ["AK-PWS", "AK-FN", "BC-N", "SS-S", "JDF", "OR", "CA-N"]:
        v = r[k]
        if v is not None:
            vals.append(f"{k}={v:.1%}")
        else:
            vals.append(f"{k}=N/A")
    print("  " + "  ".join(vals))
    print(f"  Params: {json.dumps(r['params'], default=str)}")
    print()

print("=== FULL RANKING ===")
for r in results:
    run = r["run"]
    rmsle = r["rmsle"]
    ak = r["AK-PWS"]
    orv = r["OR"]
    can = r["CA-N"]
    ak_s = f"{ak:.1%}" if ak is not None else "N/A"
    or_s = f"{orv:.1%}" if orv is not None else "N/A"
    ca_s = f"{can:.1%}" if can is not None else "N/A"
    print(f"{run}: RMSLE={rmsle:.4f}  AK-PWS={ak_s}  OR={or_s}  CA-N={ca_s}")
