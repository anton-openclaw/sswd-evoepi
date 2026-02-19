#!/usr/bin/env python3
"""Send Phase 14 report email with HTML body and key figure attachments via JMAP."""

import base64
import json
import mimetypes
import os
import sys
import requests
from pathlib import Path

API_TOKEN = os.environ.get("FASTMAIL_API_TOKEN", "")
JMAP_SESSION_URL = "https://api.fastmail.com/jmap/session"
TO = "wlweert@gmail.com"
SUBJECT = "SSWD-EvoEpi: Continuous Settlement + Spawning Overhaul — Complete Report"

# Key figures to attach (most scientifically important)
KEY_FIGURES = [
    "05_before_after_epidemic.png",
    "08_spawning_density_dependence.png",
    "50_metapopulation_timeseries.png",
    "53_fjord_vs_open.png",
    "55_node_fate_matrix.png",
    "58_simulation_dashboard.png",
    "62_spawning_event_profile.png",
    "65_spawning_before_after.png",
    "66_spawning_highdensity.gif",
    "68_spawning_param_comparison.gif",
]

FIG_DIR = Path(__file__).parent.parent / "results" / "continuous_settlement" / "final_viz"


def get_session():
    r = requests.get(JMAP_SESSION_URL, headers={
        "Authorization": f"Bearer {API_TOKEN}",
        "Content-Type": "application/json",
    })
    r.raise_for_status()
    s = r.json()
    return s["apiUrl"], s["primaryAccounts"]["urn:ietf:params:jmap:mail"], s.get("uploadUrl", "")


def upload_blob(upload_url, account_id, filepath):
    """Upload a file as a JMAP blob, return blobId."""
    mime = mimetypes.guess_type(str(filepath))[0] or "application/octet-stream"
    url = upload_url.replace("{accountId}", account_id)
    with open(filepath, "rb") as f:
        data = f.read()
    r = requests.post(url, headers={
        "Authorization": f"Bearer {API_TOKEN}",
        "Content-Type": mime,
    }, data=data)
    r.raise_for_status()
    return r.json()["blobId"], mime, len(data)


def jmap_call(api_url, method_calls):
    r = requests.post(api_url, headers={
        "Authorization": f"Bearer {API_TOKEN}",
        "Content-Type": "application/json",
    }, json={"using": ["urn:ietf:params:jmap:core", "urn:ietf:params:jmap:mail", "urn:ietf:params:jmap:submission"], "methodCalls": method_calls})
    r.raise_for_status()
    return r.json()


def build_html_body():
    """Build HTML email body."""
    return """<html><body style="font-family: Georgia, serif; max-width: 800px; margin: 0 auto; color: #222;">

<h1 style="color: #2c3e50;">SSWD-EvoEpi: Continuous Settlement + Spawning Overhaul</h1>
<p style="color: #666; font-size: 14px;">February 19, 2026 &bull; Git: <code>3b00594</code> &bull; <a href="https://github.com/anton-openclaw/sswd-evoepi">GitHub Repository</a></p>

<hr style="border: 1px solid #bdc3c7;">

<h2>Executive Summary</h2>
<p>Three major model improvements completed across 14 implementation phases:</p>
<ol>
<li><strong>Continuous larval settlement</strong> — PLD-based timing replaces annual-pulse artifact (Hodin 2021 parameterization)</li>
<li><strong>Spawning dynamics overhaul</strong> — sex-asymmetric multi-bout spawning with chemical induction cascades and Allee effects</li>
<li><strong>Juvenile immunity</strong> — age-dependent disease susceptibility with settlement-day tracking</li>
</ol>

<p><strong>Result:</strong> Sawtooth epidemic artifact eliminated. All 5 nodes crash 93-99% in 20 years, consistent with Hamilton et al. (2021). North→south mortality gradient and fjord protection confirmed. Selection signal detected at all nodes.</p>

<h2>Simulation Results (Seed 42, 20 Years)</h2>
<table style="border-collapse: collapse; width: 100%; font-size: 14px;">
<tr style="background: #2c3e50; color: white;">
  <th style="padding: 8px; text-align: left;">Node</th>
  <th style="padding: 8px;">K</th>
  <th style="padding: 8px;">Crash %</th>
  <th style="padding: 8px;">Final Pop</th>
  <th style="padding: 8px;">Δr</th>
</tr>
<tr style="background: #ecf0f1;"><td style="padding: 6px;">Sitka, AK</td><td style="text-align:center;">1,000</td><td style="text-align:center;">95.4%</td><td style="text-align:center;">45</td><td style="text-align:center;">+0.037</td></tr>
<tr><td style="padding: 6px;">Howe Sound, BC</td><td style="text-align:center;">400</td><td style="text-align:center;">93.2%</td><td style="text-align:center;">26</td><td style="text-align:center;">+0.031</td></tr>
<tr style="background: #ecf0f1;"><td style="padding: 6px;">San Juan Islands, WA</td><td style="text-align:center;">800</td><td style="text-align:center;">99.3%</td><td style="text-align:center;">5</td><td style="text-align:center;">+0.058</td></tr>
<tr><td style="padding: 6px;">Newport, OR</td><td style="text-align:center;">600</td><td style="text-align:center;">99.1%</td><td style="text-align:center;">5</td><td style="text-align:center;">+0.061</td></tr>
<tr style="background: #ecf0f1;"><td style="padding: 6px;">Monterey, CA</td><td style="text-align:center;">700</td><td style="text-align:center;">98.8%</td><td style="text-align:center;">8</td><td style="text-align:center;">+0.039</td></tr>
</table>

<h2>Key Findings</h2>

<h3>1. Continuous Settlement Eliminates Sawtooth Artifact</h3>
<p>The PLD function PLD(T) = PLD_ref × exp(E_a/R × (1/T - 1/T_ref)) spreads settlement across ~100-150 days per node. Cold-water nodes (Sitka, PLD ~120d) get wider settlement windows than warm-water nodes (Monterey, PLD ~40d).</p>
<p><strong>→ See attached: 05_before_after_epidemic.png (KEY FIGURE)</strong></p>

<h3>2. Spawning Overhaul Produces Realistic Mass Events</h3>
<p>New parameters: female 2× bouts, male zero refractory, κ_fm=0.80, κ_mf=0.60, readiness induction 0.15. Produces emergent mass spawning events through cascading chemical induction — 20-40% of population at high density, rare/sparse at low density (Allee effect).</p>
<p><strong>→ See attached GIFs — especially spawning_param_comparison.gif and spawning_highdensity.gif</strong></p>

<h3>3. Density-Dependent Allee Effect in Spawning</h3>
<p>Below N/K ≈ 0.3, spawning success drops sharply. Post-epidemic populations may enter a reproductive Allee trap — too few adults to trigger mass spawning cascades. This has critical implications for captive-bred release planning.</p>

<h3>4. Multi-Seed Robustness</h3>
<p>4-seed comparison shows max spread of 5.2 percentage points (Howe Sound, K=400). All qualitative patterns consistent. 7/7 quality gates passed.</p>

<h3>5. Performance</h3>
<p>3× slower than pre-continuous-settlement (65.5s vs 21.5s for 5-node 20yr), well within budget. Spawning overhaul adds negligible overhead.</p>

<h2>Spawning Parameter Changes</h2>
<table style="border-collapse: collapse; width: 100%; font-size: 14px;">
<tr style="background: #2c3e50; color: white;">
  <th style="padding: 8px; text-align: left;">Parameter</th>
  <th style="padding: 8px;">Old</th>
  <th style="padding: 8px;">New</th>
  <th style="padding: 8px; text-align: left;">Basis</th>
</tr>
<tr><td style="padding: 6px;">female_max_bouts</td><td style="text-align:center;">1</td><td style="text-align:center;">2</td><td>Multi-bout in broadcast spawners</td></tr>
<tr style="background: #ecf0f1;"><td style="padding: 6px;">male_refractory_days</td><td style="text-align:center;">21</td><td style="text-align:center;">0</td><td>Males respond to female cues immediately</td></tr>
<tr><td style="padding: 6px;">κ_mf (M→F induction)</td><td style="text-align:center;">0.30</td><td style="text-align:center;">0.60</td><td>Strong chemical induction</td></tr>
<tr style="background: #ecf0f1;"><td style="padding: 6px;">κ_fm (F→M induction)</td><td style="text-align:center;">0.50</td><td style="text-align:center;">0.80</td><td>Females primary spawning triggers</td></tr>
<tr><td style="padding: 6px;">readiness_induction_prob</td><td style="text-align:center;">0.0</td><td style="text-align:center;">0.15</td><td>Pre-spawning readiness cascade</td></tr>
<tr style="background: #ecf0f1;"><td style="padding: 6px;">gravity_enabled</td><td style="text-align:center;">false</td><td style="text-align:center;">true</td><td>Spawning aggregation behavior</td></tr>
</table>

<h2>Implications for SA Round 3</h2>
<p>~15 new parameters to sweep (PLD_ref, E_a_pld, self-recruitment rates, spawning induction strengths, juvenile immunity duration). The top-4 drivers from Round 1 (mu_I2D_ref, susceptibility_multiplier, a_exposure, sigma_2_eff) will likely remain dominant, but interactions with new spawning parameters may shift rankings.</p>

<h2>Complete Figure Set</h2>
<p><strong>51 figures generated</strong> (47 PNGs + 4 animated GIFs, 19 MB total):</p>
<ul>
<li>9 settlement dynamics figures</li>
<li>9 population dynamics figures</li>
<li>10 disease/epidemic figures</li>
<li>4 genetics/evolution figures</li>
<li>8 spatial pattern figures</li>
<li>4 dashboard composites</li>
<li>4 spawning detail figures</li>
<li>4 spawning animations (GIFs)</li>
</ul>
<p><strong>10 key figures attached to this email.</strong> Full set available in the repository:</p>
<p><a href="https://github.com/anton-openclaw/sswd-evoepi/tree/main/results/continuous_settlement/final_viz">→ GitHub: results/continuous_settlement/final_viz/</a></p>

<p>Full report with numbered figure catalog: <a href="https://github.com/anton-openclaw/sswd-evoepi/blob/main/results/continuous_settlement/CONTINUOUS_SETTLEMENT_REPORT.md">CONTINUOUS_SETTLEMENT_REPORT.md</a></p>

<h2>Next Steps</h2>
<ol>
<li>SA Round 3 with new parameters</li>
<li>Pathogen evolution spec (co-evolutionary module)</li>
<li>Conservation module (Monterey outplanting protocol)</li>
<li>150-node full coastline scaling</li>
<li>Joint calibration (MCMC/ABC vs Hamilton 2021 + Schiebelhut)</li>
</ol>

<hr style="border: 1px solid #bdc3c7;">
<p style="color: #999; font-size: 12px;">— Anton 🔬<br>
SSWD-EvoEpi Phase 14 &bull; Commit 3b00594 &bull; 51 figures &bull; 26,800 lines of code</p>

</body></html>"""


def main():
    if not API_TOKEN:
        print("❌ FASTMAIL_API_TOKEN not set")
        sys.exit(1)

    api_url, account_id, upload_url = get_session()
    print(f"Session OK. Account: {account_id}")

    # Upload attachments
    attachments = []
    for fname in KEY_FIGURES:
        fpath = FIG_DIR / fname
        if not fpath.exists():
            print(f"  ⚠️  Skipping {fname} (not found)")
            continue
        print(f"  Uploading {fname}...", end=" ")
        blob_id, mime, size = upload_blob(upload_url, account_id, fpath)
        attachments.append({
            "blobId": blob_id,
            "type": mime,
            "name": fname,
            "size": size,
            "disposition": "attachment",
        })
        print(f"OK ({size/1024:.0f} KB)")

    print(f"\n{len(attachments)} attachments uploaded")

    # Get drafts mailbox
    resp = jmap_call(api_url, [
        ["Mailbox/query", {"accountId": account_id, "filter": {"role": "drafts"}}, "0"]
    ])
    drafts_id = resp["methodResponses"][0][1]["ids"][0]

    # Build email
    html_body = build_html_body()
    text_body = "See HTML version of this email for the full report.\n\nFull report: https://github.com/anton-openclaw/sswd-evoepi/blob/main/results/continuous_settlement/CONTINUOUS_SETTLEMENT_REPORT.md\nFigures: https://github.com/anton-openclaw/sswd-evoepi/tree/main/results/continuous_settlement/final_viz"

    email_create = {
        "mailboxIds": {drafts_id: True},
        "from": [{"name": "Anton 🔬", "email": "antonstar@fastmail.com"}],
        "to": [{"name": "Willem Weertman", "email": TO}],
        "subject": SUBJECT,
        "htmlBody": [{"partId": "html", "type": "text/html"}],
        "textBody": [{"partId": "text", "type": "text/plain"}],
        "bodyValues": {
            "html": {"value": html_body},
            "text": {"value": text_body},
        },
        "keywords": {"$draft": True},
    }

    if attachments:
        email_create["attachments"] = attachments

    # Get identity
    id_resp = jmap_call(api_url, [
        ["Identity/get", {"accountId": account_id}, "0"]
    ])
    identities = id_resp["methodResponses"][0][1].get("list", [])
    identity_id = identities[0]["id"] if identities else None

    # Create + submit
    submission = {
        "emailId": "#draft1",
        "envelope": {
            "mailFrom": {"email": "antonstar@fastmail.com"},
            "rcptTo": [{"email": TO}],
        },
    }
    if identity_id:
        submission["identityId"] = identity_id

    resp = jmap_call(api_url, [
        ["Email/set", {
            "accountId": account_id,
            "create": {"draft1": email_create},
        }, "0"],
        ["EmailSubmission/set", {
            "accountId": account_id,
            "create": {"sub1": submission},
            "onSuccessUpdateEmail": {
                "#sub1": {
                    f"mailboxIds/{drafts_id}": None,
                    "keywords/$draft": None,
                }
            }
        }, "1"],
    ])

    # Check result
    for r in resp["methodResponses"]:
        if r[1] and r[1].get("notCreated"):
            print(f"❌ Error: {json.dumps(r[1]['notCreated'], indent=2)}")
            sys.exit(1)

    print(f"\n✅ Email sent to {TO}")
    print(f"   Subject: {SUBJECT}")
    print(f"   Attachments: {len(attachments)} files")


if __name__ == "__main__":
    main()
