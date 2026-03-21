import json, os, ast

base = os.path.expanduser("~/projects/sswd-evoepi/results/calibration/W289-W308")

# First, look at scoring structure
path = os.path.join(base, "W289", "result_seed137.json")
with open(path) as f:
    d = json.load(f)

scoring = d["scoring"]
if isinstance(scoring, str):
    scoring = ast.literal_eval(scoring)
print("scoring keys:", list(scoring.keys()) if isinstance(scoring, dict) else type(scoring))
print("scoring:", json.dumps(scoring, indent=2, default=str)[:1000])

print("\nregion_recovery:")
rr = d["region_recovery"]
if isinstance(rr, str):
    rr = ast.literal_eval(rr)
print(json.dumps(rr, indent=2, default=str)[:500])
