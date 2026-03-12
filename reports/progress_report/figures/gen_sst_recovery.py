#!/usr/bin/env python3
"""Generate dual-panel SST vs Population Recovery figure."""

import json
import os
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

# ── Paths ──────────────────────────────────────────────────────────────
SST_DIR = "/home/starbot/.openclaw/workspace/sswd-evoepi/data/sst/site_sst/"
SITES_JSON = "/home/starbot/.openclaw/workspace/sswd-evoepi/data/nodes/all_sites.json"
POP_JSON = "/home/starbot/.openclaw/workspace/sswd-evoepi/results/calibration/W142/combined_results.json"
OUT_PATH = "/home/starbot/.openclaw/workspace/sswd-evoepi/reports/progress_report/figures/fig_sst_recovery.png"

REGIONS = ["AK-PWS", "AK-FN", "SS-S", "CA-N"]
REGION_LABELS = {
    "AK-PWS": "Alaska – Prince William Sound",
    "AK-FN": "Alaska – Fairweather/N. Gulf",
    "SS-S": "Salish Sea – South",
    "CA-N": "California – North",
}
COLORS = {
    "AK-PWS": "#1b9e77",   # teal green
    "AK-FN": "#d95f02",    # orange
    "SS-S":  "#7570b3",    # purple
    "CA-N":  "#e7298a",    # magenta/pink
}

YEAR_MIN, YEAR_MAX = 2012, 2024
POP_YEARS = list(range(2012, 2025))  # 13 years

# ── Load sites ─────────────────────────────────────────────────────────
with open(SITES_JSON) as f:
    all_sites = json.load(f)

region_site_names = {}
for r in REGIONS:
    region_site_names[r] = [s["name"] for s in all_sites if s["region"] == r]
    print(f"  {r}: {len(region_site_names[r])} sites")

# ── Load & aggregate SST ──────────────────────────────────────────────
region_sst = {}
for r in REGIONS:
    frames = []
    for name in region_site_names[r]:
        fp = os.path.join(SST_DIR, f"{name}_monthly.csv")
        if not os.path.exists(fp):
            continue
        df = pd.read_csv(fp)
        df = df[(df["year"] >= YEAR_MIN) & (df["year"] <= YEAR_MAX)]
        frames.append(df)
    big = pd.concat(frames, ignore_index=True)
    # Monthly mean across all sites in the region
    monthly = big.groupby(["year", "month"])["sst"].mean().reset_index()
    monthly["date"] = pd.to_datetime(
        monthly["year"].astype(str) + "-" + monthly["month"].astype(str).str.zfill(2) + "-15"
    )
    monthly = monthly.sort_values("date")
    region_sst[r] = monthly
    print(f"  SST {r}: {len(monthly)} months, range {monthly['sst'].min():.1f}–{monthly['sst'].max():.1f} °C")

# ── Load population data ──────────────────────────────────────────────
with open(POP_JSON) as f:
    pop_data = json.load(f)
rd = pop_data["results"][0]["region_details"]

region_pop = {}
for r in REGIONS:
    totals = np.array(rd[r]["yearly_totals"], dtype=float)
    # Normalize by year-0 (2012)
    region_pop[r] = totals / totals[0]
    print(f"  Pop {r}: initial={totals[0]:.0f}, min_frac={region_pop[r].min():.3f}, final_frac={region_pop[r][-1]:.3f}")

# ── Plot ──────────────────────────────────────────────────────────────
plt.style.use("seaborn-v0_8-whitegrid")
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), dpi=150, sharex=True,
                                gridspec_kw={"hspace": 0.12})

# ── Top panel: SST ────────────────────────────────────────────────────
for r in REGIONS:
    d = region_sst[r]
    ax1.plot(d["date"], d["sst"], color=COLORS[r], linewidth=1.0,
             label=REGION_LABELS[r], alpha=0.85)

# Threshold lines
ax1.axhline(12, color="firebrick", linestyle="--", linewidth=1.2, alpha=0.7)
ax1.text(pd.Timestamp("2024-04-01"), 12.3,
         r"$T_{\mathrm{vbnc}}$ = 12 °C (pathogen activation)",
         fontsize=8.5, color="firebrick", va="bottom", ha="right")

ax1.axhline(9, color="steelblue", linestyle="--", linewidth=1.2, alpha=0.7)
ax1.text(pd.Timestamp("2024-04-01"), 8.3,
         r"$T_{\mathrm{vbnc,min}}$ = 9 °C (thermal floor)",
         fontsize=8.5, color="steelblue", va="top", ha="right")

ax1.set_ylabel("Sea Surface Temperature (°C)", fontsize=11)
ax1.set_xlim(pd.Timestamp("2012-01-01"), pd.Timestamp("2025-01-01"))
ax1.legend(loc="upper left", fontsize=8.5, framealpha=0.9, ncol=2)
ax1.set_title("Regional Mean SST and Population Recovery — SSWD EvoEpi", fontsize=13, fontweight="bold")
ax1.tick_params(axis="x", labelsize=9)
ax1.tick_params(axis="y", labelsize=9)

# Format x-axis as years
ax1.xaxis.set_major_locator(mdates.YearLocator())
ax1.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))

# ── Bottom panel: Population ──────────────────────────────────────────
# Convert years to timestamps (mid-year) so x-axis aligns with SST panel
pop_dates = [pd.Timestamp(f"{y}-07-01") for y in POP_YEARS]
for r in REGIONS:
    ax2.plot(pop_dates, region_pop[r], color=COLORS[r], linewidth=2.0,
             marker="o", markersize=4, label=REGION_LABELS[r], alpha=0.85)

# Disease introduction line
ax2.axvline(pd.Timestamp("2013-06-01"), color="gray", linestyle="--", linewidth=1.2, alpha=0.7)
ax2.text(pd.Timestamp("2013-08-01"), 0.95, "Disease\nintroduction", fontsize=8.5, color="gray",
         va="top", ha="left")

ax2.set_ylabel("Population (fraction of initial)", fontsize=11)
ax2.set_xlabel("Year", fontsize=11)
ax2.legend(loc="upper right", fontsize=8.5, framealpha=0.9, ncol=2)
ax2.tick_params(axis="x", labelsize=9)
ax2.tick_params(axis="y", labelsize=9)

# Reference line at 1.0
ax2.axhline(1.0, color="black", linestyle=":", linewidth=0.8, alpha=0.4)

plt.tight_layout()
fig.savefig(OUT_PATH, dpi=150, bbox_inches="tight", facecolor="white")
plt.close()

# Report
size_kb = os.path.getsize(OUT_PATH) / 1024
print(f"\n✅ Saved: {OUT_PATH}")
print(f"   Size: {size_kb:.0f} KB")
