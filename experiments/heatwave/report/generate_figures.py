#!/usr/bin/env python3
"""Generate all visualization figures for the MHW temperature report.

Creates 8 publication-quality figures in experiments/heatwave/report/figures/.
"""

import json
import os
import sys
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns

# Setup
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, REPO)

TEMPLATE_PATH = os.path.join(REPO, "experiments", "heatwave", "templates", "blob_anomalies.npz")
SCENARIOS_DIR = os.path.join(REPO, "data", "sst", "scenarios")
SITES_JSON = os.path.join(REPO, "data", "nodes", "all_sites.json")
FIG_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")

os.makedirs(FIG_DIR, exist_ok=True)

# Region order S→N
REGION_ORDER = [
    "BJ", "CA-S", "CA-C", "CA-N", "OR", "WA-O", "JDF",
    "SS-S", "SS-N", "BC-C", "BC-N",
    "AK-FS", "AK-FN", "AK-OC", "AK-PWS", "AK-EG", "AK-WG", "AK-AL"
]

MONTHS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
          "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

# Style
plt.style.use("seaborn-v0_8-paper")
plt.rcParams.update({
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "figure.dpi": 150,
    "savefig.dpi": 300,
})


def load_template():
    """Load blob anomaly template."""
    tmpl = np.load(TEMPLATE_PATH, allow_pickle=True)
    return {
        "site_names": list(tmpl["site_names"]),
        "anomaly_2014": tmpl["anomaly_2014"],
        "anomaly_2015": tmpl["anomaly_2015"],
        "anomaly_2016": tmpl["anomaly_2016"],
        "climatology": tmpl["climatology"],
    }


def load_sites():
    """Load site info with regions."""
    with open(SITES_JSON) as f:
        return json.load(f)


def site_to_region(site_name):
    """Extract region from site name."""
    # Region is everything except the trailing -NNN
    parts = site_name.rsplit("-", 1)
    return parts[0]


def get_region_mean_anomalies(site_names, anomalies):
    """Compute mean anomaly per region. Returns dict region -> (12,) array."""
    region_data = {r: [] for r in REGION_ORDER}
    for i, name in enumerate(site_names):
        region = site_to_region(name)
        if region in region_data:
            region_data[region].append(anomalies[i])
    result = {}
    for r in REGION_ORDER:
        if region_data[r]:
            result[r] = np.mean(region_data[r], axis=0)
        else:
            result[r] = np.zeros(12)
    return result


def load_scenario_sst(scenario_id, site_names):
    """Load SST data for a scenario. Returns dict: site_name -> {(year,month): sst}."""
    d = os.path.join(SCENARIOS_DIR, scenario_id)
    data = {}
    for name in site_names:
        path = os.path.join(d, f"{name}_monthly.csv")
        site_data = {}
        with open(path) as f:
            reader = csv.DictReader(f)
            for row in reader:
                site_data[(int(row["year"]), int(row["month"]))] = float(row["sst"])
        data[name] = site_data
    return data


def compute_scenario_regional_anomalies(scenario_data, climatology, site_names, years=range(2025, 2051)):
    """Compute regional mean SST anomalies for a scenario over specified years.

    Returns (n_regions, n_years) array.
    """
    n_years = len(list(years))
    region_anomalies = np.zeros((len(REGION_ORDER), n_years))

    for ri, region in enumerate(REGION_ORDER):
        # Find sites in this region
        site_indices = [i for i, n in enumerate(site_names) if site_to_region(n) == region]
        if not site_indices:
            continue

        for yi, year in enumerate(years):
            annual_anomalies = []
            for si in site_indices:
                name = site_names[si]
                clim = climatology[si]
                for month in range(1, 13):
                    sst = scenario_data[name].get((year, month), clim[month-1])
                    anom = sst - clim[month-1]
                    annual_anomalies.append(anom)
            region_anomalies[ri, yi] = np.mean(annual_anomalies)

    return region_anomalies


# ============================================================
# Figure 1: Blob Template
# ============================================================
def fig_blob_template(tmpl):
    print("  Generating fig_blob_template...")
    site_names = tmpl["site_names"]

    fig, axes = plt.subplots(1, 3, figsize=(12, 5), sharey=True)

    for ax, year, key in zip(axes, [2014, 2015, 2016],
                              ["anomaly_2014", "anomaly_2015", "anomaly_2016"]):
        region_means = get_region_mean_anomalies(site_names, tmpl[key])
        data = np.array([region_means[r] for r in REGION_ORDER])

        vmax = 3.0
        im = ax.imshow(data, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                       interpolation="nearest")
        ax.set_xticks(range(12))
        ax.set_xticklabels(MONTHS, rotation=45, ha="right")
        ax.set_title(f"{year} Anomaly", fontweight="bold")
        ax.set_xlabel("Month")

    axes[0].set_yticks(range(len(REGION_ORDER)))
    axes[0].set_yticklabels(REGION_ORDER)
    axes[0].set_ylabel("Region (S → N)")

    fig.subplots_adjust(right=0.88)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label("SST Anomaly (°C)")

    fig.suptitle("Pacific Blob Anomaly Template (vs. 2002–2012 Climatology)", fontweight="bold", y=1.02)
    plt.savefig(os.path.join(FIG_DIR, "fig_blob_template.pdf"))
    plt.close()


# ============================================================
# Figure 2: Blob Spatial
# ============================================================
def fig_blob_spatial(tmpl, sites):
    print("  Generating fig_blob_spatial...")
    site_names = tmpl["site_names"]
    anom_2015 = tmpl["anomaly_2015"]

    # Get lat and annual mean anomaly per site
    site_lats = {}
    for s in sites:
        site_lats[s["name"]] = s.get("latitude", s.get("lat"))

    lats = []
    mean_anoms = []
    for i, name in enumerate(site_names):
        if name in site_lats:
            lats.append(site_lats[name])
            mean_anoms.append(np.mean(anom_2015[i]))

    lats = np.array(lats)
    mean_anoms = np.array(mean_anoms)

    fig, ax = plt.subplots(figsize=(8, 5))
    sc = ax.scatter(lats, mean_anoms, c=mean_anoms, cmap="RdBu_r",
                    vmin=-2, vmax=3, s=8, alpha=0.6, edgecolors="none")
    ax.axhline(0, color="gray", linewidth=0.5, linestyle="--")
    ax.set_xlabel("Latitude (°N)")
    ax.set_ylabel("Mean SST Anomaly 2015 (°C)")
    ax.set_title("Pacific Blob 2015: Spatial Gradient of Peak Anomaly", fontweight="bold")

    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label("SST Anomaly (°C)")

    # Add region labels
    for region in REGION_ORDER:
        region_sites = [i for i, n in enumerate(site_names) if site_to_region(n) == region and n in site_lats]
        if region_sites:
            rlats = [site_lats[site_names[i]] for i in region_sites]
            mean_lat = np.mean(rlats)
            ranoms = [np.mean(anom_2015[i]) for i in region_sites]
            mean_a = np.mean(ranoms)
            ax.annotate(region, (mean_lat, mean_a), fontsize=6, alpha=0.7,
                       ha="center", va="bottom")

    plt.savefig(os.path.join(FIG_DIR, "fig_blob_spatial.pdf"))
    plt.close()


# ============================================================
# Figure 3: Scenario Comparison (HERO figure)
# ============================================================
def fig_scenario_comparison(tmpl):
    print("  Generating fig_scenario_comparison...")
    site_names = tmpl["site_names"]
    climatology = tmpl["climatology"]

    scenarios = [
        "MHW_CTRL", "MHW_S1_LOW", "MHW_S1_MED", "MHW_S1_HIGH",
        "MHW_R10_MED", "MHW_R7_MED", "MHW_R5_MED", "MHW_R5_HIGH",
        "MHW_R3_MED", "MHW_RAMP", "MHW_CHRONIC", "MHW_COMPOUND"
    ]

    labels = [
        "Control", "Single Low (0.5×)", "Single Med (1.0×)", "Single High (1.5×)",
        "Recurring 10yr", "Recurring 7yr", "Recurring 5yr", "Recurring 5yr High",
        "Recurring 3yr", "Ramping", "Chronic +1°C", "Compound"
    ]

    years = list(range(2025, 2051))
    vmax = 2.5

    fig, axes = plt.subplots(4, 3, figsize=(14, 12), sharex=True, sharey=True)

    for idx, (scenario, label) in enumerate(zip(scenarios, labels)):
        ax = axes[idx // 3, idx % 3]

        print(f"    Loading {scenario}...")
        sdata = load_scenario_sst(scenario, site_names)
        regional = compute_scenario_regional_anomalies(sdata, climatology, site_names, years)

        im = ax.imshow(regional, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                       interpolation="nearest", origin="lower")

        ax.set_title(label, fontsize=8, fontweight="bold")

        if idx % 3 == 0:
            ax.set_yticks(range(len(REGION_ORDER)))
            ax.set_yticklabels(REGION_ORDER, fontsize=6)
        if idx >= 9:
            # Every 5 years
            tick_positions = [i for i, y in enumerate(years) if y % 5 == 0]
            ax.set_xticks(tick_positions)
            ax.set_xticklabels([years[i] for i in tick_positions], fontsize=7)

    fig.subplots_adjust(right=0.88, hspace=0.3, wspace=0.1)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.015, 0.7])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label("SST Anomaly (°C vs climatology)")

    fig.suptitle("Marine Heatwave SST Scenarios: Regional Anomalies 2025–2050",
                 fontweight="bold", fontsize=12, y=0.98)
    plt.savefig(os.path.join(FIG_DIR, "fig_scenario_comparison.pdf"))
    plt.close()


# ============================================================
# Figure 4: Intensity Comparison
# ============================================================
def fig_intensity_comparison(tmpl):
    print("  Generating fig_intensity_comparison...")
    site_names = tmpl["site_names"]
    climatology = tmpl["climatology"]

    scenarios = ["MHW_S1_LOW", "MHW_S1_MED", "MHW_S1_HIGH"]
    labels = ["Low (0.5×)", "Medium (1.0×)", "High (1.5×)"]
    years = list(range(2025, 2051))
    vmax = 2.5

    fig, axes = plt.subplots(1, 3, figsize=(14, 5), sharey=True)

    for ax, scenario, label in zip(axes, scenarios, labels):
        sdata = load_scenario_sst(scenario, site_names)
        regional = compute_scenario_regional_anomalies(sdata, climatology, site_names, years)

        im = ax.imshow(regional, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                       interpolation="nearest", origin="lower")
        ax.set_title(label, fontweight="bold")
        tick_positions = [i for i, y in enumerate(years) if y % 5 == 0]
        ax.set_xticks(tick_positions)
        ax.set_xticklabels([years[i] for i in tick_positions])
        ax.set_xlabel("Year")

    axes[0].set_yticks(range(len(REGION_ORDER)))
    axes[0].set_yticklabels(REGION_ORDER)
    axes[0].set_ylabel("Region (S → N)")

    fig.subplots_adjust(right=0.88)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label("SST Anomaly (°C)")

    fig.suptitle("Single-Event Intensity Comparison (2030 Event)", fontweight="bold", y=1.02)
    plt.savefig(os.path.join(FIG_DIR, "fig_intensity_comparison.pdf"))
    plt.close()


# ============================================================
# Figure 5: Frequency Comparison
# ============================================================
def fig_frequency_comparison(tmpl):
    print("  Generating fig_frequency_comparison...")
    site_names = tmpl["site_names"]
    climatology = tmpl["climatology"]

    scenarios = ["MHW_R10_MED", "MHW_R7_MED", "MHW_R5_MED", "MHW_R3_MED"]
    labels = ["Every 10 years", "Every 7 years", "Every 5 years", "Every 3 years"]
    years = list(range(2025, 2051))
    vmax = 2.5

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True, sharey=True)

    for idx, (scenario, label) in enumerate(zip(scenarios, labels)):
        ax = axes[idx // 2, idx % 2]
        sdata = load_scenario_sst(scenario, site_names)
        regional = compute_scenario_regional_anomalies(sdata, climatology, site_names, years)

        im = ax.imshow(regional, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                       interpolation="nearest", origin="lower")
        ax.set_title(label, fontweight="bold")

        if idx % 2 == 0:
            ax.set_yticks(range(len(REGION_ORDER)))
            ax.set_yticklabels(REGION_ORDER, fontsize=7)
        if idx >= 2:
            tick_positions = [i for i, y in enumerate(years) if y % 5 == 0]
            ax.set_xticks(tick_positions)
            ax.set_xticklabels([years[i] for i in tick_positions])

    fig.subplots_adjust(right=0.88, hspace=0.25)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label("SST Anomaly (°C)")

    fig.suptitle("Heatwave Recurrence Frequency Comparison (1.0× Intensity)",
                 fontweight="bold", y=0.98)
    plt.savefig(os.path.join(FIG_DIR, "fig_frequency_comparison.pdf"))
    plt.close()


# ============================================================
# Figure 6: Regional SST Time Series
# ============================================================
def fig_regional_sst_timeseries(tmpl):
    print("  Generating fig_regional_sst_timeseries...")
    site_names = tmpl["site_names"]

    regions = ["AK-PWS", "JDF", "OR", "CA-S"]
    scenarios = ["MHW_CTRL", "MHW_S1_MED", "MHW_R5_MED", "MHW_R5_HIGH"]
    scenario_labels = ["Control", "Single 1.0×", "Recurring 5yr", "Recurring 5yr 1.5×"]
    colors = ["#333333", "#2196F3", "#FF9800", "#E53935"]

    years = list(range(2012, 2051))

    # Load all scenario data
    scenario_data = {}
    for scen in scenarios:
        print(f"    Loading {scen}...")
        scenario_data[scen] = load_scenario_sst(scen, site_names)

    fig, axes = plt.subplots(4, 4, figsize=(16, 12), sharex=True)

    for ri, region in enumerate(regions):
        # Sites in this region
        region_sites = [i for i, n in enumerate(site_names) if site_to_region(n) == region]

        for si, (scen, slabel) in enumerate(zip(scenarios, scenario_labels)):
            ax = axes[ri, si]

            # Compute regional mean SST per year
            annual_sst = []
            for year in years:
                vals = []
                for site_i in region_sites:
                    name = site_names[site_i]
                    for month in range(1, 13):
                        sst = scenario_data[scen][name].get((year, month))
                        if sst is not None:
                            vals.append(sst)
                annual_sst.append(np.mean(vals) if vals else np.nan)

            ax.plot(years, annual_sst, color=colors[si], linewidth=1.2)
            ax.axvline(2025, color="gray", linewidth=0.5, linestyle="--", alpha=0.5)

            if ri == 0:
                ax.set_title(slabel, fontsize=9, fontweight="bold")
            if si == 0:
                ax.set_ylabel(f"{region}\nSST (°C)", fontsize=8)
            if ri == 3:
                ax.set_xlabel("Year")

            ax.tick_params(labelsize=7)

    fig.suptitle("Regional Mean SST Trajectories (2012–2050)", fontweight="bold", fontsize=12, y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(os.path.join(FIG_DIR, "fig_regional_sst_timeseries.pdf"))
    plt.close()


# ============================================================
# Figure 7: Peak Anomaly Map (bar chart)
# ============================================================
def fig_peak_anomaly_map(tmpl):
    print("  Generating fig_peak_anomaly_map...")
    site_names = tmpl["site_names"]
    climatology = tmpl["climatology"]
    anom_2015 = tmpl["anomaly_2015"]

    # Real Blob peak (2015): regional mean annual anomaly
    blob_peak = get_region_mean_anomalies(site_names, anom_2015)
    blob_annual = [np.mean(blob_peak[r]) for r in REGION_ORDER]

    # Simulated: load S1_MED, find max anomaly during 2030-2031
    print("    Loading MHW_S1_MED...")
    sdata = load_scenario_sst("MHW_S1_MED", site_names)

    sim_peak = []
    for region in REGION_ORDER:
        region_sites = [i for i, n in enumerate(site_names) if site_to_region(n) == region]
        max_anom = -999
        for year in [2030, 2031]:
            anoms = []
            for si in region_sites:
                name = site_names[si]
                clim = climatology[si]
                for month in range(1, 13):
                    sst = sdata[name].get((year, month), clim[month-1])
                    anoms.append(sst - clim[month-1])
            yr_mean = np.mean(anoms)
            max_anom = max(max_anom, yr_mean)
        sim_peak.append(max_anom)

    x = np.arange(len(REGION_ORDER))
    width = 0.35

    fig, ax = plt.subplots(figsize=(10, 5))
    bars1 = ax.bar(x - width/2, blob_annual, width, label="Real Blob (2015)", color="#1976D2", alpha=0.8)
    bars2 = ax.bar(x + width/2, sim_peak, width, label="Simulated (2030–2031)", color="#E53935", alpha=0.8)

    ax.set_xticks(x)
    ax.set_xticklabels(REGION_ORDER, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Mean SST Anomaly (°C)")
    ax.set_xlabel("Region (S → N)")
    ax.set_title("Peak Anomaly: Real Blob vs. Simulated MHW_S1_MED", fontweight="bold")
    ax.legend()
    ax.axhline(0, color="gray", linewidth=0.5)

    plt.savefig(os.path.join(FIG_DIR, "fig_peak_anomaly_map.pdf"))
    plt.close()


# ============================================================
# Figure 8: Chronic vs Pulse
# ============================================================
def fig_chronic_vs_pulse(tmpl):
    print("  Generating fig_chronic_vs_pulse...")
    site_names = tmpl["site_names"]
    climatology = tmpl["climatology"]

    scenarios = ["MHW_CHRONIC", "MHW_R5_MED"]
    labels = ["Chronic +1°C", "Pulsed (5yr recurrence)"]
    colors = ["#E53935", "#1976D2"]

    years = list(range(2025, 2051))

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Panel 1: Cumulative heat exposure per region
    ax = axes[0]
    for scen, label, color in zip(scenarios, labels, colors):
        print(f"    Loading {scen}...")
        sdata = load_scenario_sst(scen, site_names)

        # Compute cumulative anomaly over all years, averaged across regions
        cumulative = np.zeros(len(years))
        for yi, year in enumerate(years):
            total_anom = 0
            count = 0
            for i, name in enumerate(site_names):
                clim = climatology[i]
                for month in range(1, 13):
                    sst = sdata[name].get((year, month), clim[month-1])
                    total_anom += sst - clim[month-1]
                    count += 1
            cumulative[yi] = total_anom / count if count > 0 else 0

        cumulative_sum = np.cumsum(cumulative)
        ax.plot(years, cumulative_sum, label=label, color=color, linewidth=2)

    ax.set_xlabel("Year")
    ax.set_ylabel("Cumulative Anomaly (°C·years)")
    ax.set_title("Cumulative Heat Exposure", fontweight="bold")
    ax.legend()

    # Panel 2: Annual anomaly comparison
    ax = axes[1]
    for scen, label, color in zip(scenarios, labels, colors):
        sdata = load_scenario_sst(scen, site_names)
        annual_anom = []
        for year in years:
            anoms = []
            for i, name in enumerate(site_names):
                clim = climatology[i]
                for month in range(1, 13):
                    sst = sdata[name].get((year, month), clim[month-1])
                    anoms.append(sst - clim[month-1])
            annual_anom.append(np.mean(anoms))

        ax.plot(years, annual_anom, label=label, color=color, linewidth=2)
        ax.fill_between(years, 0, annual_anom, alpha=0.15, color=color)

    ax.set_xlabel("Year")
    ax.set_ylabel("Mean SST Anomaly (°C)")
    ax.set_title("Annual Mean Anomaly", fontweight="bold")
    ax.axhline(0, color="gray", linewidth=0.5, linestyle="--")
    ax.legend()

    fig.suptitle("Chronic vs. Pulsed Warming Comparison", fontweight="bold", y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, "fig_chronic_vs_pulse.pdf"))
    plt.close()


# ============================================================
# Main
# ============================================================
def main():
    print("Loading data...")
    tmpl = load_template()
    sites = load_sites()

    print("Generating figures...")
    fig_blob_template(tmpl)
    fig_blob_spatial(tmpl, sites)
    fig_scenario_comparison(tmpl)
    fig_intensity_comparison(tmpl)
    fig_frequency_comparison(tmpl)
    fig_regional_sst_timeseries(tmpl)
    fig_peak_anomaly_map(tmpl)
    fig_chronic_vs_pulse(tmpl)

    print("\nAll figures generated in:", FIG_DIR)
    for f in sorted(os.listdir(FIG_DIR)):
        size = os.path.getsize(os.path.join(FIG_DIR, f))
        print(f"  {f}: {size/1024:.1f} KB")


if __name__ == "__main__":
    main()
