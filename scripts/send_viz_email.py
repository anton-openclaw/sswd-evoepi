#!/usr/bin/env python3
"""Send visualization library email with all figures as attachments via JMAP."""

import json
import os
import sys
import base64
import mimetypes
import requests
from requests.adapters import HTTPAdapter
from pathlib import Path

# Force IPv4 to avoid "Network unreachable" on IPv6-broken hosts
import urllib3.util.connection
urllib3.util.connection.HAS_IPV6 = False

API_TOKEN = os.environ.get("FASTMAIL_API_TOKEN", "")
JMAP_SESSION_URL = "https://api.fastmail.com/jmap/session"

VIZ_DIR = Path(__file__).parent.parent / "results" / "viz_examples"

# --- JMAP helpers ---

def get_session():
    r = requests.get(JMAP_SESSION_URL, headers={
        "Authorization": f"Bearer {API_TOKEN}",
        "Content-Type": "application/json"
    })
    r.raise_for_status()
    session = r.json()
    api_url = session["apiUrl"]
    account_id = session["primaryAccounts"]["urn:ietf:params:jmap:mail"]
    upload_url = session["uploadUrl"].replace("{accountId}", account_id)
    return api_url, account_id, upload_url


def jmap_call(api_url, method_calls):
    r = requests.post(api_url, headers={
        "Authorization": f"Bearer {API_TOKEN}",
        "Content-Type": "application/json"
    }, json={
        "using": [
            "urn:ietf:params:jmap:core",
            "urn:ietf:params:jmap:mail",
            "urn:ietf:params:jmap:submission"
        ],
        "methodCalls": method_calls
    })
    r.raise_for_status()
    return r.json()


def upload_blob(upload_url, file_path):
    """Upload a file and return the blob ID."""
    mime = mimetypes.guess_type(str(file_path))[0] or "application/octet-stream"
    with open(file_path, "rb") as f:
        data = f.read()
    r = requests.post(upload_url, headers={
        "Authorization": f"Bearer {API_TOKEN}",
        "Content-Type": mime,
    }, data=data)
    r.raise_for_status()
    result = r.json()
    return result["blobId"], result["size"], mime


# --- Email content ---

EMAIL_HTML = """\
<html>
<head>
<style>
  body { font-family: Georgia, serif; color: #222; max-width: 900px; margin: 0 auto; padding: 20px; }
  h1 { color: #1a5276; border-bottom: 2px solid #1a5276; padding-bottom: 8px; }
  h2 { color: #2e86c1; margin-top: 40px; border-bottom: 1px solid #ccc; padding-bottom: 4px; }
  .figure { margin: 20px 0; padding: 12px; background: #f9f9f9; border: 1px solid #ddd; border-radius: 4px; }
  .fig-num { font-weight: bold; color: #1a5276; }
  .fig-file { font-family: monospace; font-size: 0.85em; color: #666; }
  .caption { margin-top: 6px; line-height: 1.5; }
  .category-count { color: #666; font-size: 0.9em; }
  .summary { background: #eaf2f8; padding: 16px; border-radius: 6px; margin: 20px 0; }
  .toc a { text-decoration: none; color: #2e86c1; }
  .toc li { margin: 4px 0; }
</style>
</head>
<body>

<h1>SSWD-EvoEpi: Visualization Library &mdash; Complete Figure Set with Captions</h1>

<div class="summary">
<strong>63 figures</strong> across 6 categories, generated from a 5-node, 20-year prototype simulation
(Sitka, Howe Sound, San Juan Islands, Newport, Monterey; seed=42).<br><br>
All figures use a dark publication theme with colorblind-safe palettes. Figures are attached as PNGs &mdash;
see inline captions below for interpretation guide.
</div>

<h3>Table of Contents</h3>
<ul class="toc">
  <li><a href="#population">1. Population &amp; Demographics</a> (10 figures)</li>
  <li><a href="#disease">2. Disease &amp; Epidemiology</a> (12 figures)</li>
  <li><a href="#genetics">3. Host Genetics &amp; Evolution</a> (12 figures)</li>
  <li><a href="#coevolution">4. Pathogen Co-Evolution</a> (8 figures)</li>
  <li><a href="#spatial">5. Spatial &amp; Metapopulation</a> (8 figures)</li>
  <li><a href="#dashboards">6. Dashboards &amp; Composite Views</a> (13 figures)</li>
</ul>

<hr>

<!-- ==================== POPULATION ==================== -->
<h2 id="population">1. Population &amp; Demographics <span class="category-count">(10 figures)</span></h2>

<div class="figure">
  <span class="fig-num">Figure 1.</span> <span class="fig-file">population/01_population_trajectory.png</span>
  <div class="caption">
    <strong>Total population over time.</strong> Y-axis: total agent count across all nodes; x-axis: simulation years.
    A horizontal dashed line marks aggregate carrying capacity (K). The vertical marker indicates year of first disease introduction.
    <em>Key takeaway:</em> Look for the magnitude and speed of the post-epidemic crash, and whether any stabilization or recovery plateau emerges.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 2.</span> <span class="fig-file">population/02_stage_composition.png</span>
  <div class="caption">
    <strong>Life-stage composition over time (stacked area).</strong> Y-axis: number of individuals; x-axis: years. Colors represent life stages (juvenile, sub-adult, adult).
    <em>Key takeaway:</em> Reveals whether disease disproportionately removes adults (expected for SSWD) and whether juvenile recruitment compensates. A healthy population shows stable stage ratios; a crashing one shows stage truncation.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 3.</span> <span class="fig-file">population/03_cause_of_death.png</span>
  <div class="caption">
    <strong>Annual deaths by cause (stacked bar).</strong> Y-axis: death count per year; x-axis: years. Colors encode cause: disease (SSWD), natural mortality, and senescence.
    <em>Key takeaway:</em> During epidemic years, disease deaths should dwarf background mortality. In post-epidemic years, look for whether natural mortality resumes dominance &mdash; a sign of population stabilization.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 4.</span> <span class="fig-file">population/04_age_size_pyramid.png</span>
  <div class="caption">
    <strong>Population pyramid at a snapshot year.</strong> Left panel: age distribution; right panel: size (disc diameter) distribution. Males extend left, females right.
    <em>Key takeaway:</em> A bottom-heavy pyramid indicates recent recruitment success; a truncated top indicates adult mortality from SSWD. Compare pre- and post-epidemic shapes to see demographic damage.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 5.</span> <span class="fig-file">population/05_population_heatmap.png</span>
  <div class="caption">
    <strong>Population heatmap: nodes &times; years.</strong> Y-axis: node names (ordered by latitude); x-axis: years. Color intensity represents population as a fraction of local K.
    <em>Key takeaway:</em> The north&ndash;south mortality gradient should be visible &mdash; northern nodes retain more color (higher N/K) than southern ones. Fjord-protected nodes (Howe Sound) should stand out as warmer colors persisting longer.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 6.</span> <span class="fig-file">population/06_recruitment_timeseries.png</span>
  <div class="caption">
    <strong>Annual recruitment (settlers) with reproductive adult overlay.</strong> Y-axis (left): new settlers per year; y-axis (right): number of reproductive adults. X-axis: years.
    <em>Key takeaway:</em> Recruitment should crash after adult population collapses (Allee effects, sweepstake reproduction). A lag between adult decline and recruitment decline reflects larval development time.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 7.</span> <span class="fig-file">population/07_survival_curves.png</span>
  <div class="caption">
    <strong>Kaplan&ndash;Meier-style survival curves by cohort.</strong> Y-axis: survival probability; x-axis: age (years). Each line is a birth-year cohort.
    <em>Key takeaway:</em> Pre-epidemic cohorts should show higher survival. Epidemic-era cohorts show steep early drop-offs. Post-epidemic cohorts may show intermediate survival if resistance alleles have increased in frequency.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 8.</span> <span class="fig-file">population/08_sex_ratio.png</span>
  <div class="caption">
    <strong>Proportion female over time.</strong> Y-axis: fraction female (0&ndash;1); x-axis: years. Dashed line at 0.5 marks parity.
    <em>Key takeaway:</em> Should remain near 0.5 unless disease has sex-biased effects. Deviations indicate potential issues with spawning dynamics (sex-asymmetric immunosuppression, differential mortality).
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 9.</span> <span class="fig-file">population/09_density_dependence.png</span>
  <div class="caption">
    <strong>Density dependence: population size vs. per-capita growth rate.</strong> X-axis: N/K (relative density); y-axis: per-capita growth rate r. Each point is one year.
    <em>Key takeaway:</em> Should show negative density dependence (downward slope) &mdash; confirming Beverton&ndash;Holt regulation is working. Points far below zero at low density indicate disease-driven decline overwhelming compensatory recruitment.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 10.</span> <span class="fig-file">population/10_node_comparison.png</span>
  <div class="caption">
    <strong>Initial vs. final population by node (grouped bars).</strong> Y-axis: population count; x-axis: nodes. Paired bars show year-0 vs. year-20 population, colored by subregion.
    <em>Key takeaway:</em> Quantifies per-node crash severity. Monterey and Newport should show >95% decline; Howe Sound (fjord) should retain the most. This is the key figure for spatial heterogeneity in outcomes.
  </div>
</div>

<!-- ==================== DISEASE ==================== -->
<h2 id="disease">2. Disease &amp; Epidemiology <span class="category-count">(12 figures)</span></h2>

<div class="figure">
  <span class="fig-num">Figure 11.</span> <span class="fig-file">disease/01_epidemic_curve.png</span>
  <div class="caption">
    <strong>Classic epidemic curve (all SEIPD+R compartments).</strong> Y-axis: number of individuals; x-axis: time (days or years). Stacked/overlaid areas show Susceptible, Exposed, Infected-stage-1 (I&#8321;), Infected-stage-2 (I&#8322;), Post-diseased, and Recovered.
    <em>Key takeaway:</em> The shape reveals epidemic speed, peak prevalence, and whether the population reaches an endemic equilibrium or crashes to near-zero. The recovery compartment (R) should be very small given SSWD's high CFR.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 12.</span> <span class="fig-file">disease/02_vibrio_concentration.png</span>
  <div class="caption">
    <strong>Environmental <em>Vibrio</em> concentration over time (log scale).</strong> Y-axis: log&#8321;&#8320; pathogen concentration in water column; x-axis: time.
    <em>Key takeaway:</em> Tracks environmental reservoir dynamics. Concentration should spike during epidemic peaks (mass shedding from infected hosts) and decay during low-prevalence periods. Temperature-dependent VBNC (viable but non-culturable) modulation may create seasonal oscillations.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 13.</span> <span class="fig-file">disease/03_force_of_infection.png</span>
  <div class="caption">
    <strong>Histogram of per-individual force of infection (&lambda;<sub>i</sub>) at a snapshot.</strong> X-axis: &lambda;<sub>i</sub> (daily infection probability); y-axis: count of individuals.
    <em>Key takeaway:</em> Reveals heterogeneity in exposure risk across the population. Resistant individuals should cluster at lower &lambda;<sub>i</sub> values. A long right tail indicates superspreader-driven dynamics.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 14.</span> <span class="fig-file">disease/04_R0_over_time.png</span>
  <div class="caption">
    <strong>Effective reproductive number R&#8320; over time.</strong> Y-axis: R&#8320; estimate; x-axis: years. A horizontal dashed line at R&#8320;=1 marks the epidemic threshold.
    <em>Key takeaway:</em> R&#8320; &gt; 1 means the epidemic is growing; R&#8320; &lt; 1 means it is declining. Watch for R&#8320; dropping below 1 as the susceptible population is depleted &mdash; this is depletion-driven, not immunity-driven, decline.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 15.</span> <span class="fig-file">disease/05_disease_mortality_by_node.png</span>
  <div class="caption">
    <strong>Cumulative disease mortality fraction per node.</strong> Y-axis: fraction of initial population killed by disease; x-axis: nodes (ordered by latitude). Color encodes latitude.
    <em>Key takeaway:</em> Should show the characteristic north&ndash;south gradient &mdash; warmer southern nodes suffer higher cumulative mortality (temperature-dependent <em>V. pectenicida</em> growth, T<sub>opt</sub>=20&deg;C).
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 16.</span> <span class="fig-file">disease/06_epidemic_wave_timing.png</span>
  <div class="caption">
    <strong>Epidemic wave propagation diagram.</strong> Y-axis: nodes (ordered by latitude or connectivity); x-axis: time. Horizontal bars show onset, peak, and end of the epidemic at each node.
    <em>Key takeaway:</em> Reveals the speed and directionality of spatial spread. If the epidemic starts at a single node and propagates via larval/adult movement, the timing sequence should correlate with network distance.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 17.</span> <span class="fig-file">disease/07_compartment_flow_sankey.png</span>
  <div class="caption">
    <strong>Compartment flow (Sankey-style) diagram.</strong> Nodes: disease compartments (S, E, I&#8321;, I&#8322;, D, R). Flow width proportional to total individual-transitions over the simulation.
    <em>Key takeaway:</em> Immediately shows what fraction of the population progresses to death vs. recovery. The S&rarr;E&rarr;I&#8321;&rarr;I&#8322;&rarr;D pathway should dominate massively over the recovery branch, consistent with SSWD's &gt;95% CFR.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 18.</span> <span class="fig-file">disease/08_shedding_timeseries.png</span>
  <div class="caption">
    <strong>Total pathogen shedding over time, decomposed by source.</strong> Y-axis: total shedding rate; x-axis: time. Colors distinguish I&#8321; (early) vs. I&#8322; (late-stage) shedding contributions.
    <em>Key takeaway:</em> Late-stage (I&#8322;) shedding should dominate total pathogen output due to the ~3,000&times; lab-to-field scaling factor. This drives the environmental reservoir and is the primary transmission mechanism.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 19.</span> <span class="fig-file">disease/09_disease_state_heatmap.png</span>
  <div class="caption">
    <strong>Disease state heatmap: nodes &times; years.</strong> Y-axis: nodes; x-axis: years. Color intensity represents disease prevalence (or mortality fraction) at each node-year.
    <em>Key takeaway:</em> Complements the population heatmap (Fig. 5) &mdash; bright cells indicate active outbreaks. Look for whether the epidemic is synchronous across nodes or shows traveling-wave structure.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 20.</span> <span class="fig-file">disease/10_immunosuppression_overlap.png</span>
  <div class="caption">
    <strong>Spawning season, immunosuppression, and VBNC overlap timeline.</strong> Three overlaid bands showing: spawning window (Nov&ndash;Jul), post-spawning immunosuppression (28-day window, 2&times; disease susceptibility), and warm-season <em>Vibrio</em> resurgence.
    <em>Key takeaway:</em> The critical biological hypothesis &mdash; spawning-induced immunosuppression coincides temporally with summer <em>Vibrio</em> emergence, creating a synergistic vulnerability window. This figure tests whether the model captures that overlap.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 21.</span> <span class="fig-file">disease/11_recovery_vs_resistance.png</span>
  <div class="caption">
    <strong>Scatter: host resistance (x) vs. recovery probability (y).</strong> Each point is an individual. Color indicates disease outcome (recovered vs. died).
    <em>Key takeaway:</em> Recovery probability p<sub>rec</sub> = &rho;<sub>rec</sub> &times; r<sub>i</sub>&sup2; (quadratic). Only individuals in the upper tail of resistance have meaningful recovery probability. This visualizes the selection mechanism &mdash; resistance doesn't prevent infection, it enables recovery.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 22.</span> <span class="fig-file">disease/12_cfr_over_time.png</span>
  <div class="caption">
    <strong>Case fatality rate (CFR) per year.</strong> Y-axis: CFR (disease deaths / total infections); x-axis: years.
    <em>Key takeaway:</em> Should start very high (~95&ndash;99%) and potentially decline slightly over time as selection increases mean resistance. If pathogen evolution is enabled, coevolutionary dynamics may further modify CFR trajectory (cf. myxomatosis analog).
  </div>
</div>

<!-- ==================== GENETICS ==================== -->
<h2 id="genetics">3. Host Genetics &amp; Evolution <span class="category-count">(12 figures)</span></h2>

<div class="figure">
  <span class="fig-num">Figure 23.</span> <span class="fig-file">genetics/01_resistance_trajectory.png</span>
  <div class="caption">
    <strong>Mean population resistance over time (&plusmn;1 SD band).</strong> Y-axis: mean resistance r&#772;; x-axis: years. Shaded region shows &plusmn;1 standard deviation.
    <em>Key takeaway:</em> The central evolutionary signal. A rising trajectory indicates directional selection for resistance. Compare the rate of increase (&Delta;r&#772;/year) with Schiebelhut et al.'s empirical &Delta;q &asymp; 0.08&ndash;0.15 per generation to assess biological plausibility.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 24.</span> <span class="fig-file">genetics/02_resistance_distribution.png</span>
  <div class="caption">
    <strong>Resistance distribution at multiple timepoints (overlaid KDE curves).</strong> X-axis: resistance value r<sub>i</sub>; y-axis: density. Each curve is a different year.
    <em>Key takeaway:</em> Watch for the distribution shifting rightward (increasing mean) AND narrowing (decreasing variance) over time. Variance depletion is the signature of directional selection eroding additive genetic variation &mdash; the fuel for evolutionary rescue.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 25.</span> <span class="fig-file">genetics/03_allele_freq_spaghetti.png</span>
  <div class="caption">
    <strong>Per-locus allele frequency trajectories (&ldquo;spaghetti plot&rdquo;).</strong> Y-axis: resistance allele frequency q; x-axis: years. Each line is one of the 51 additive loci, colored by effect size.
    <em>Key takeaway:</em> Large-effect loci (warm colors) should show faster frequency increases. Small-effect loci may be dominated by drift. Convergence of many loci toward q=1 would indicate a hard selective sweep; stochastic trajectories indicate drift-selection balance.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 26.</span> <span class="fig-file">genetics/04_additive_variance.png</span>
  <div class="caption">
    <strong>Additive genetic variance (V<sub>A</sub>) over time.</strong> Y-axis: V<sub>A</sub>; x-axis: years.
    <em>Key takeaway:</em> V<sub>A</sub> is the raw material for evolutionary rescue. It should initially increase (as rare resistance alleles increase) then decrease as alleles approach fixation. The rate of V<sub>A</sub> depletion determines how long selection can act &mdash; if V<sub>A</sub> &rarr; 0 before the population reaches a viable resistance level, evolutionary rescue fails.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 27.</span> <span class="fig-file">genetics/05_ef1a_dynamics.png</span>
  <div class="caption">
    <strong>EF1A overdominant locus: allele frequency and heterozygosity.</strong> Y-axis (left): allele frequency; y-axis (right): heterozygosity (2pq). X-axis: years.
    <em>Key takeaway:</em> The EF1A analog provides large resistance via heterozygote advantage (s<sub>het</sub>) but cannot fix &mdash; balancing selection maintains polymorphism. Heterozygosity should remain high even as additive loci shift. This locus alone contributes ~0.058 to mean resistance.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 28.</span> <span class="fig-file">genetics/06_selection_differential.png</span>
  <div class="caption">
    <strong>Selection differential (&Delta;r) per year.</strong> Y-axis: difference in mean resistance between survivors and the pre-selection population; x-axis: years.
    <em>Key takeaway:</em> Positive values indicate viability selection is operating &mdash; individuals who survive to reproduce have higher resistance than the cohort average. &Delta;r should be largest during peak epidemic years and diminish as the susceptible population is depleted.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 29.</span> <span class="fig-file">genetics/07_heritability.png</span>
  <div class="caption">
    <strong>Narrow-sense heritability (h&sup2; = V<sub>A</sub>/V<sub>P</sub>) over time.</strong> Y-axis: h&sup2;; x-axis: years.
    <em>Key takeaway:</em> h&sup2; determines how efficiently selection translates into evolutionary response (breeder's equation: R = h&sup2;S). Environmental variance is minimal in this model, so h&sup2; should be high (~0.5&ndash;0.9). Declining h&sup2; signals V<sub>A</sub> depletion.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 30.</span> <span class="fig-file">genetics/08_genotype_phenotype_map.png</span>
  <div class="caption">
    <strong>Genotype&ndash;phenotype map scatter.</strong> X-axis: genotypic score (sum of allele counts &times; effect sizes); y-axis: realized resistance r<sub>i</sub>. Each point is an individual.
    <em>Key takeaway:</em> Should show a tight linear relationship (high h&sup2;) with the EF1A locus creating discrete bands (heterozygotes shifted up). Spread around the line reflects the overdominant locus contribution, not environmental noise.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 31.</span> <span class="fig-file">genetics/09_effect_size_distribution.png</span>
  <div class="caption">
    <strong>Locus effect sizes (ranked bar chart).</strong> X-axis: locus index (ranked by effect size); y-axis: per-allele effect on resistance.
    <em>Key takeaway:</em> Shows the genetic architecture of resistance. Effect sizes are drawn from an exponential distribution (many small-effect, few large-effect loci). This architecture is critical &mdash; it determines whether evolution is driven by a few major loci or many minor ones.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 32.</span> <span class="fig-file">genetics/10_resistance_by_node_violin.png</span>
  <div class="caption">
    <strong>Resistance distributions per node (violin plot).</strong> X-axis: node name; y-axis: resistance r<sub>i</sub>. Violin width shows the density of individuals at each resistance level.
    <em>Key takeaway:</em> Spatial differentiation in resistance evolution. Southern nodes (stronger selection) should show higher mean resistance. Fjord-protected nodes may lag because weaker selection pressure preserves more low-resistance individuals.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 33.</span> <span class="fig-file">genetics/11_genetic_drift_null.png</span>
  <div class="caption">
    <strong>Drift null model vs. disease simulation (multi-seed).</strong> Multiple seed overlays comparing allele frequency trajectories with and without disease.
    <em>Key takeaway:</em> The critical control. Drift-only trajectories (no disease) should show random walks around initial frequencies. Disease trajectories should show consistent directional shifts above the drift envelope, confirming that observed allele changes are selection-driven, not stochastic.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 34.</span> <span class="fig-file">genetics/12_beta_init_visualization.png</span>
  <div class="caption">
    <strong>Beta(a,b) initialization visualization (3 panels).</strong> Left: Beta(2,8) PDF showing per-locus allele frequency distribution. Center: resulting per-locus q values across the 51 loci. Right: population resistance distribution after initialization.
    <em>Key takeaway:</em> Documents the new genetics initialization scheme. Beta(2,8) produces most loci with low initial frequency (q &lt; 0.3) but allows some at higher frequencies, yielding a realistic distribution of standing genetic variation with target mean r&#772; &asymp; 0.15.
  </div>
</div>

<!-- ==================== COEVOLUTION ==================== -->
<h2 id="coevolution">4. Pathogen Co-Evolution <span class="category-count">(8 figures)</span></h2>

<div class="figure">
  <span class="fig-num">Figure 35.</span> <span class="fig-file">coevolution/virulence_trajectory.png</span>
  <div class="caption">
    <strong>Mean pathogen virulence over time (&plusmn;SD band).</strong> Y-axis: mean virulence v&#772;; x-axis: years.
    <em>Key takeaway:</em> The pathogen-side evolutionary signal. Virulence should evolve toward an intermediate optimum set by the transmission&ndash;virulence tradeoff. Too high = kills hosts before sufficient shedding; too low = outcompeted by faster-replicating strains.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 36.</span> <span class="fig-file">coevolution/coevolution_phase_portrait.png</span>
  <div class="caption">
    <strong>Co-evolutionary phase portrait: mean host resistance (x) vs. mean pathogen virulence (y).</strong> Each point is a year; trajectory shows co-evolutionary dynamics.
    <em>Key takeaway:</em> THE key figure for co-evolutionary dynamics. Look for: (a) arms-race dynamics (both increasing), (b) convergence to a coexistence equilibrium, or (c) cycling. Compare with Clement et al. 2024 DFTD predictions where coevolution enabled host persistence.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 37.</span> <span class="fig-file">coevolution/virulence_distribution_over_time.png</span>
  <div class="caption">
    <strong>Virulence distribution at multiple timepoints (violin plots).</strong> X-axis: year; y-axis: virulence phenotype. Violin width shows density.
    <em>Key takeaway:</em> Reveals whether pathogen diversity narrows (selective sweep toward optimal virulence) or broadens (diversification into multiple strategies). Bimodal distributions would suggest coexistence of high- and low-virulence strains.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 38.</span> <span class="fig-file">coevolution/tradeoff_curve.png</span>
  <div class="caption">
    <strong>Virulence&ndash;transmission tradeoff curve.</strong> X-axis: virulence v; y-axis: transmission &beta;(v). The mechanistic constraint: shedding rate &times; duration &asymp; constant (sicker hosts shed faster but die sooner).
    <em>Key takeaway:</em> This is the fundamental constraint preventing unconstrained virulence evolution. The curve shape (concave = intermediate optimum exists) determines whether the pathogen evolves toward intermediate or maximal virulence. R&#8320;-maximizing virulence is marked.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 39.</span> <span class="fig-file">coevolution/R0_by_virulence.png</span>
  <div class="caption">
    <strong>R&#8320; as a function of virulence for different host densities.</strong> X-axis: virulence; y-axis: R&#8320;. Multiple curves for different population sizes (N).
    <em>Key takeaway:</em> Shows how the fitness landscape shifts as the host population declines. The R&#8320;-maximizing virulence may shift as N decreases, potentially selecting for LOWER virulence when hosts are scarce &mdash; the "prudent parasite" effect.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 40.</span> <span class="fig-file">coevolution/virulence_vs_host_density.png</span>
  <div class="caption">
    <strong>Mean virulence vs. total host population over time.</strong> X-axis: N (total population); y-axis: mean virulence. Points colored by time.
    <em>Key takeaway:</em> Tests the density-dependent virulence evolution hypothesis. If virulence decreases as N drops, the pathogen is adapting to host scarcity. This could create a natural stabilizing mechanism preventing complete host extinction.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 41.</span> <span class="fig-file">coevolution/strain_competition.png</span>
  <div class="caption">
    <strong>Strain competition: stacked area of virulence bins over time.</strong> Y-axis: proportion of infections; x-axis: time. Colors: low, medium, and high virulence bins.
    <em>Key takeaway:</em> Shows which virulence strategies are winning over evolutionary time. A shift from high to intermediate virulence would mirror the classic myxomatosis pattern. Complete dominance by one bin suggests a selective sweep.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 42.</span> <span class="fig-file">coevolution/coevolution_multi_seed.png</span>
  <div class="caption">
    <strong>Multi-seed co-evolution trajectories.</strong> Multiple simulation seeds overlaid: host resistance (one color) and pathogen virulence (another) over time.
    <em>Key takeaway:</em> Assesses robustness of co-evolutionary patterns across stochastic realizations. Consistent trajectories across seeds = deterministic signal; divergent trajectories = strong stochastic effects (especially likely in small populations post-crash).
  </div>
</div>

<!-- ==================== SPATIAL ==================== -->
<h2 id="spatial">5. Spatial &amp; Metapopulation <span class="category-count">(8 figures)</span></h2>

<div class="figure">
  <span class="fig-num">Figure 43.</span> <span class="fig-file">spatial/network_map.png</span>
  <div class="caption">
    <strong>Geographic map of the metapopulation network.</strong> Nodes plotted at lat/lon coordinates, sized by carrying capacity, colored by a metric (e.g., final population). Lines show connectivity (larval dispersal or adult movement corridors).
    <em>Key takeaway:</em> Provides geographic context for all spatial results. The 5-node prototype spans Sitka, AK to Monterey, CA &mdash; over 3,000 km of Pacific coastline.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 44.</span> <span class="fig-file">spatial/connectivity_heatmap.png</span>
  <div class="caption">
    <strong>Larval connectivity matrix (heatmap).</strong> Y-axis: source node; x-axis: destination node. Color intensity represents connectivity strength (C<sub>ij</sub>).
    <em>Key takeaway:</em> Shows which nodes supply larvae to which. Diagonal elements (self-recruitment) should be highest, especially for fjord nodes (&alpha;<sub>self,fjord</sub>=0.30 vs. &alpha;<sub>self,open</sub>=0.10). Off-diagonal elements reveal the rescue potential of connected populations.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 45.</span> <span class="fig-file">spatial/north_south_gradient_mortality.png</span>
  <div class="caption">
    <strong>North&ndash;south gradient in disease mortality.</strong> Y-axis: latitude; x-axis: cumulative mortality fraction (or other metric). Each point is a node.
    <em>Key takeaway:</em> THE validation figure. Hamilton et al. 2021 documented a strong north&ndash;south gradient in SSWD mortality across the NE Pacific. This should be reproduced: higher mortality at lower latitudes where warmer temperatures promote <em>V. pectenicida</em> growth.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 46.</span> <span class="fig-file">spatial/fjord_vs_open.png</span>
  <div class="caption">
    <strong>Fjord vs. open-coast paired comparison.</strong> Side-by-side panels comparing Howe Sound (fjord) against open-coast nodes on key metrics: population trajectory, disease prevalence, and resistance evolution.
    <em>Key takeaway:</em> Fjords serve as natural refugia due to higher self-recruitment, potential temperature buffering, and reduced larval connectivity (limiting pathogen import). Howe Sound should consistently outperform open-coast nodes.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 47.</span> <span class="fig-file">spatial/metapopulation_timeseries.png</span>
  <div class="caption">
    <strong>All nodes' population trajectories overlaid.</strong> Y-axis: population; x-axis: years. Each line is one node, colored/styled distinctly.
    <em>Key takeaway:</em> Shows synchrony vs. asynchrony across the metapopulation. Asynchronous dynamics are favorable for persistence (portfolio effect) &mdash; if all nodes crash simultaneously, rescue via recolonization is impossible.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 48.</span> <span class="fig-file">spatial/larval_flow_diagram.png</span>
  <div class="caption">
    <strong>Net larval dispersal flow diagram.</strong> Arrow diagram connecting nodes, with arrow width proportional to net larval flow. Direction shows source&rarr;sink relationships.
    <em>Key takeaway:</em> Identifies source populations (net larval exporters) and sink populations (net importers). Post-epidemic, the remaining source populations become critical for metapopulation recovery and should be prioritized for conservation.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 49.</span> <span class="fig-file">spatial/spatial_epidemic_timeline.png</span>
  <div class="caption">
    <strong>Spatial epidemic timeline (Gantt-chart style).</strong> Y-axis: nodes; x-axis: time. Horizontal bars color-coded by disease state (no disease / outbreak / endemic / cleared).
    <em>Key takeaway:</em> Visualizes the spatiotemporal pattern of epidemic progression across the network. Look for sequential vs. simultaneous outbreaks, and whether any nodes clear the epidemic before others &mdash; potential recolonization sources.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 50.</span> <span class="fig-file">spatial/node_fate_matrix.png</span>
  <div class="caption">
    <strong>Node fate matrix: one small multi-panel per node.</strong> Each sub-panel shows population (line), disease prevalence (fill), and mean resistance (overlay) for a single node.
    <em>Key takeaway:</em> The complete per-node story at a glance. Enables rapid comparison of trajectory shapes &mdash; which nodes crash-and-stabilize vs. crash-to-extinction. Resistance trends should be steepest at nodes with strongest selection (highest mortality).
  </div>
</div>

<!-- ==================== DASHBOARDS ==================== -->
<h2 id="dashboards">6. Dashboards &amp; Composite Views <span class="category-count">(13 figures)</span></h2>

<div class="figure">
  <span class="fig-num">Figure 51.</span> <span class="fig-file">dashboards/simulation_dashboard_pe.png</span>
  <div class="caption">
    <strong>Master simulation dashboard (pathogen evolution ON).</strong> 2&times;3 grid: (top-left) population over time, (top-center) epidemic curve, (top-right) mean resistance trajectory, (bottom-left) cumulative deaths by cause, (bottom-center) mean virulence trajectory, (bottom-right) co-evolutionary phase portrait.
    <em>Key takeaway:</em> The single-figure summary of a full co-evolutionary simulation. With pathogen evolution enabled, the bottom row reveals pathogen-side dynamics alongside host evolution.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 52.</span> <span class="fig-file">dashboards/simulation_dashboard_nope.png</span>
  <div class="caption">
    <strong>Master simulation dashboard (pathogen evolution OFF).</strong> Same 2&times;3 layout but the virulence panel is replaced with environmental <em>Vibrio</em> concentration.
    <em>Key takeaway:</em> The baseline comparison &mdash; host evolution only, static pathogen. Compare with Fig. 51 to see how pathogen co-evolution modifies outcomes.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 53.</span> <span class="fig-file">dashboards/spatial_dashboard.png</span>
  <div class="caption">
    <strong>Multi-node spatial overview dashboard.</strong> 2&times;3 grid: (top-left) per-node population trajectories, (top-center) disease mortality bar chart, (top-right) network map, (bottom-left) per-node resistance trajectories, (bottom-center) north&ndash;south gradient, (bottom-right) fjord vs. open comparison.
    <em>Key takeaway:</em> Integrates all spatial dimensions into one figure. Ideal for presentations &mdash; tells the complete geographic story of the epidemic.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 54.</span> <span class="fig-file">dashboards/scenario_comparison.png</span>
  <div class="caption">
    <strong>Scenario comparison grid (4 rows &times; N columns).</strong> Rows: population, resistance, virulence, deaths. Columns: different pathogen evolution scenarios (e.g., PE off, v&#8320;=0.3, v&#8320;=0.5).
    <em>Key takeaway:</em> Side-by-side comparison of how different assumptions about pathogen evolution change predictions. Essential for assessing whether co-evolution is a critical model ingredient or a second-order effect.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 55.</span> <span class="fig-file">dashboards/sensitivity_tornado_popcrash.png</span>
  <div class="caption">
    <strong>Sobol sensitivity tornado for population crash metric.</strong> Horizontal bars: top 15 parameters ranked by total-order Sobol index (S<sub>T</sub>). Inner bar = first-order (S<sub>1</sub>), full bar = total-order (S<sub>T</sub>). Color by parameter category.
    <em>Key takeaway:</em> The gap between S<sub>1</sub> and S<sub>T</sub> reveals interaction strength. susceptibility_multiplier and sigma_2_eff dominate, both with massive interaction effects. This confirms that population crash is driven by parameter interactions, not individual parameters alone.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 56.</span> <span class="fig-file">dashboards/sensitivity_tornado_resistance.png</span>
  <div class="caption">
    <strong>Sobol sensitivity tornado for mean resistance shift metric.</strong> Same format as Fig. 55 but for evolutionary outcome (&Delta;r&#772;).
    <em>Key takeaway:</em> mu_I2D_ref (I&#8322;&rarr;death rate) dominates &mdash; it controls who dies and therefore which genotypes survive. n_additive (number of loci) is #2, confirming that genetic architecture is a key uncertainty for evolutionary predictions.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 57.</span> <span class="fig-file">dashboards/sensitivity_heatmap.png</span>
  <div class="caption">
    <strong>Sobol total-order index heatmap: parameters &times; metrics.</strong> Y-axis: 20+ parameters (grouped by category); x-axis: 14 output metrics. Color = S<sub>T</sub> magnitude. Top 5 cells highlighted.
    <em>Key takeaway:</em> The complete sensitivity landscape in one figure. Reveals which parameters matter for which outcomes &mdash; e.g., a_exposure dominates fjord protection but not genetic metrics. Parameters with high S<sub>T</sub> across many metrics are priority calibration targets.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 58.</span> <span class="fig-file">dashboards/interaction_web_popcrash.png</span>
  <div class="caption">
    <strong>Parameter interaction matrix for population crash.</strong> Matrix plot: S<sub>T</sub>&minus;S<sub>1</sub> (approximate pairwise interaction strength) between parameters. Color intensity indicates interaction magnitude.
    <em>Key takeaway:</em> Confirms the sensitivity analysis finding that interactions dominate &mdash; extinction risk is ENTIRELY interaction-driven (S<sub>1</sub> &asymp; 0 for most parameters). The strongest interactions involve sigma_2_eff &times; mu_I2D_ref and a_exposure &times; susceptibility_multiplier.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 59.</span> <span class="fig-file">dashboards/morris_vs_sobol.png</span>
  <div class="caption">
    <strong>Morris &mu;* vs. Sobol S<sub>T</sub> scatter.</strong> Each point is a parameter; x-axis: Morris &mu;* (screening); y-axis: Sobol S<sub>T</sub> (variance decomposition). Color by category.
    <em>Key takeaway:</em> Tests whether the cheap screening method (Morris, 480 runs) agrees with the expensive decomposition (Sobol, 12,288 runs). Strong correlation = Morris is a reliable pre-screen. Outliers reveal parameters where nonlinear interactions make screening unreliable.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 60.</span> <span class="fig-file">dashboards/morris_vs_sobol_resistance.png</span>
  <div class="caption">
    <strong>Morris vs. Sobol scatter for resistance shift.</strong> Same format as Fig. 59 but for the evolutionary metric.
    <em>Key takeaway:</em> Evolutionary outcomes may show weaker Morris&ndash;Sobol agreement than demographic outcomes, because selection differentials depend on complex nonlinear interactions between mortality, genetic architecture, and population size.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 61.</span> <span class="fig-file">dashboards/evolutionary_rescue.png</span>
  <div class="caption">
    <strong>Evolutionary rescue assessment (4 panels).</strong> (a) Population vs. minimum viable threshold; (b) rate of resistance increase &Delta;r&#772;/year; (c) V<sub>A</sub> depletion trajectory; (d) estimated generations to recovery.
    <em>Key takeaway:</em> Synthesizes all evidence on whether evolutionary rescue can save Pycnopodia. Key question: does resistance increase fast enough relative to population decline? Our prototype suggests rescue is too slow (50&ndash;100+ generations) for conservation relevance &mdash; reinforcing the case for captive breeding.
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 62.</span> <span class="fig-file">dashboards/conservation_matrix.png</span>
  <div class="caption">
    <strong>Conservation scenario comparison matrix (3&times;3 grid).</strong> Scenarios: no intervention, captive breeding, gene flow enhancement. Metrics: population trajectory, resistance evolution, extinction probability.
    <em>Key takeaway:</em> Placeholder data currently &mdash; this visualization is designed to compare conservation interventions once the conservation module is built. The framework is ready for the outplanting protocol modeling (cf. December 2025 Monterey outplanting: 47/48 survived 4 weeks).
  </div>
</div>

<div class="figure">
  <span class="fig-num">Figure 63.</span> <span class="fig-file">dashboards/model_validation.png</span>
  <div class="caption">
    <strong>Model validation panel (3 panels).</strong> Model predictions vs. empirical data: (a) mortality range vs. Hamilton 2021, (b) recovery timeline vs. observed absence durations, (c) decline pattern shape comparison.
    <em>Key takeaway:</em> The critical reality check. Model should reproduce the observed 86&ndash;99.8% crash range (Hamilton 2021), the monotonic decline without recovery (no bounce-back observed empirically), and the north&ndash;south gradient. Discrepancies identify where model assumptions need revision.
  </div>
</div>

<hr>

<p><em>All 63 figures attached as PNGs. Generated from SSWD-EvoEpi visualization library (sswd_evoepi/viz/), 2026-02-18.</em></p>
<p><em>Questions, corrections, or requests for additional views &mdash; just reply.</em></p>
<p>&mdash; Anton 🔬</p>

</body>
</html>
"""

# Plain text version (fallback)
EMAIL_TEXT = """SSWD-EvoEpi: Visualization Library — Complete Figure Set with Captions

63 figures across 6 categories attached as PNGs. See HTML version for full captions.

Categories:
1. Population & Demographics (10 figures)
2. Disease & Epidemiology (12 figures)
3. Host Genetics & Evolution (12 figures)
4. Pathogen Co-Evolution (8 figures)
5. Spatial & Metapopulation (8 figures)
6. Dashboards & Composite Views (13 figures)

Generated from SSWD-EvoEpi viz library, 2026-02-18.
— Anton
"""


def collect_figures():
    """Collect all figure paths in order."""
    categories = [
        ("population", [
            "01_population_trajectory.png", "02_stage_composition.png",
            "03_cause_of_death.png", "04_age_size_pyramid.png",
            "05_population_heatmap.png", "06_recruitment_timeseries.png",
            "07_survival_curves.png", "08_sex_ratio.png",
            "09_density_dependence.png", "10_node_comparison.png",
        ]),
        ("disease", [
            "01_epidemic_curve.png", "02_vibrio_concentration.png",
            "03_force_of_infection.png", "04_R0_over_time.png",
            "05_disease_mortality_by_node.png", "06_epidemic_wave_timing.png",
            "07_compartment_flow_sankey.png", "08_shedding_timeseries.png",
            "09_disease_state_heatmap.png", "10_immunosuppression_overlap.png",
            "11_recovery_vs_resistance.png", "12_cfr_over_time.png",
        ]),
        ("genetics", [
            "01_resistance_trajectory.png", "02_resistance_distribution.png",
            "03_allele_freq_spaghetti.png", "04_additive_variance.png",
            "05_ef1a_dynamics.png", "06_selection_differential.png",
            "07_heritability.png", "08_genotype_phenotype_map.png",
            "09_effect_size_distribution.png", "10_resistance_by_node_violin.png",
            "11_genetic_drift_null.png", "12_beta_init_visualization.png",
        ]),
        ("coevolution", [
            "virulence_trajectory.png", "coevolution_phase_portrait.png",
            "virulence_distribution_over_time.png", "tradeoff_curve.png",
            "R0_by_virulence.png", "virulence_vs_host_density.png",
            "strain_competition.png", "coevolution_multi_seed.png",
        ]),
        ("spatial", [
            "network_map.png", "connectivity_heatmap.png",
            "north_south_gradient_mortality.png", "fjord_vs_open.png",
            "metapopulation_timeseries.png", "larval_flow_diagram.png",
            "spatial_epidemic_timeline.png", "node_fate_matrix.png",
        ]),
        ("dashboards", [
            "simulation_dashboard_pe.png", "simulation_dashboard_nope.png",
            "spatial_dashboard.png", "scenario_comparison.png",
            "sensitivity_tornado_popcrash.png", "sensitivity_tornado_resistance.png",
            "sensitivity_heatmap.png", "interaction_web_popcrash.png",
            "morris_vs_sobol.png", "morris_vs_sobol_resistance.png",
            "evolutionary_rescue.png", "conservation_matrix.png",
            "model_validation.png",
        ]),
    ]

    figures = []
    for cat, filenames in categories:
        for fname in filenames:
            path = VIZ_DIR / cat / fname
            if path.exists():
                figures.append(path)
            else:
                print(f"⚠️  Missing: {cat}/{fname}")
    return figures


def main():
    if not API_TOKEN:
        print("ERROR: FASTMAIL_API_TOKEN not set")
        sys.exit(1)

    api_url, account_id, upload_url = get_session()
    print(f"Session OK. Account: {account_id}")

    # Collect figures
    figures = collect_figures()
    print(f"\nFound {len(figures)} figures")

    if len(figures) < 60:
        print(f"ERROR: Only {len(figures)} figures found, expected 63. Aborting.")
        sys.exit(1)

    # Upload all blobs
    attachments = []
    for i, path in enumerate(figures):
        blob_id, size, mime = upload_blob(upload_url, path)
        attachments.append({
            "blobId": blob_id,
            "type": mime,
            "name": f"{path.parent.name}/{path.name}",
            "size": size,
        })
        print(f"  [{i+1}/{len(figures)}] Uploaded {path.parent.name}/{path.name} ({size:,} bytes)")

    # Get drafts mailbox
    resp = jmap_call(api_url, [
        ["Mailbox/query", {"accountId": account_id, "filter": {"role": "drafts"}}, "0"]
    ])
    drafts_id = resp["methodResponses"][0][1]["ids"][0]

    # Get sending identity
    id_resp = jmap_call(api_url, [
        ["Identity/get", {"accountId": account_id}, "0"]
    ])
    identities = id_resp["methodResponses"][0][1].get("list", [])
    identity_id = identities[0]["id"] if identities else None

    # Create and send email
    email_body = {
        "accountId": account_id,
        "create": {
            "draft1": {
                "mailboxIds": {drafts_id: True},
                "from": [{"name": "Anton 🔬", "email": "antonstar@fastmail.com"}],
                "to": [{"name": "Willem Weertman", "email": "wlweert@gmail.com"}],
                "subject": "SSWD-EvoEpi: Visualization Library — All 63 Figures with Captions",
                "htmlBody": [{"partId": "html", "type": "text/html"}],
                "textBody": [{"partId": "text", "type": "text/plain"}],
                "bodyValues": {
                    "html": {"value": EMAIL_HTML},
                    "text": {"value": EMAIL_TEXT},
                },
                "attachments": attachments,
                "keywords": {"$draft": True},
            }
        }
    }

    submission = {
        "emailId": "#draft1",
        "envelope": {
            "mailFrom": {"email": "antonstar@fastmail.com"},
            "rcptTo": [{"email": "wlweert@gmail.com"}],
        }
    }
    if identity_id:
        submission["identityId"] = identity_id

    resp = jmap_call(api_url, [
        ["Email/set", email_body, "0"],
        ["EmailSubmission/set", {
            "accountId": account_id,
            "create": {"sub1": submission},
            "onSuccessUpdateEmail": {
                "#sub1": {
                    f"mailboxIds/{drafts_id}": None,
                    "keywords/$draft": None,
                }
            }
        }, "1"]
    ])

    # Check result
    responses = resp["methodResponses"]
    for r in responses:
        if r[1] and r[1].get("notCreated"):
            print(f"ERROR: {json.dumps(r[1]['notCreated'], indent=2)}")
            sys.exit(1)

    print(f"\n✅ Email sent to wlweert@gmail.com with {len(attachments)} figures attached!")


if __name__ == "__main__":
    main()
