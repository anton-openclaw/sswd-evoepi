#!/usr/bin/env python3
"""Send the wavefront calibration R2 report via Fastmail JMAP."""

import os, json, base64, requests, sys, argparse

parser = argparse.ArgumentParser()
parser.add_argument("--pdf", default="/home/starbot/.openclaw/workspace/sswd-evoepi/reports/wavefront_calibration_r2/main.pdf")
parser.add_argument("--subject", default="SSWD-EvoEpi: Wavefront Calibration R2 Report (W05-W16)")
parser.add_argument("--to", default="wlweert@gmail.com")
args = parser.parse_args()

token = os.environ.get("FASTMAIL_API_TOKEN")
if not token:
    print("ERROR: FASTMAIL_API_TOKEN not set")
    sys.exit(1)

headers_json = {"Authorization": f"Bearer {token}", "Content-Type": "application/json"}

session_resp = requests.get("https://api.fastmail.com/.well-known/jmap",
                            headers={"Authorization": f"Bearer {token}"})
session_resp.raise_for_status()
session = session_resp.json()
jmap_url = session["apiUrl"]
upload_url_template = session["uploadUrl"]
account_id = list(session["accounts"].keys())[0]
upload_url = upload_url_template.replace("{accountId}", account_id)
identity_id = "172611807"

from_email = "antonstar@fastmail.com"

body_text = """Hi Willem,

Attached is the Wavefront Calibration Round 2 report for the SSWD-EvoEpi model (rounds W05-W16).

These are the first 12 runs with three major code fixes:
1. Configurable VBNC sigmoid steepness (k_vbnc)
2. Fixed Beverton-Holt recruitment pipeline (true settler counts)
3. Reduced settler survival (s0 = 0.001)

Key findings:
- Recovery gradient now working: 0% in California, 3-7% in Salish Sea, 89% in Alaska
- Best round: W15 (k=2.0, thresh=500, P_env=2000) with RMSE = 0.404
- All 8 target regions within 5x of targets (up from 4/8 in W01-W04)
- RMSE improved from 1.18 to 0.40 (66% reduction)
- P_env_max=5000 is definitively too aggressive (RMSE ~1.4-1.6)
- AK-PWS overshoots: 89% vs 50% target (disease doesn't reach Alaska)
- Wavefront timing still needs work: doesn't propagate to Alaska

Next batch (W17-W28) focuses on:
- Lowering K_half from 200K to 50-100K to reduce Alaska overshoot
- Exploring s0 = 0.002-0.003 to help CA-N reach 0.1% target

The report includes 6 figures, full results tables, and detailed analysis.

Best,
Anton (via OpenClaw)
"""

# Upload PDF
print(f"Uploading {args.pdf}...")
with open(args.pdf, "rb") as f:
    pdf_data = f.read()

upload_resp = requests.post(
    upload_url,
    headers={"Authorization": f"Bearer {token}", "Content-Type": "application/pdf"},
    data=pdf_data
)
upload_resp.raise_for_status()
blob_id = upload_resp.json()["blobId"]
blob_size = upload_resp.json()["size"]
print(f"Uploaded blob: {blob_id} ({blob_size} bytes)")

# Get Drafts mailbox
resp = requests.post(jmap_url, headers=headers_json, json={
    "using": ["urn:ietf:params:jmap:core", "urn:ietf:params:jmap:mail"],
    "methodCalls": [["Mailbox/query", {"accountId": account_id, "filter": {"role": "drafts"}}, "0"]]
})
resp.raise_for_status()
drafts_id = resp.json()["methodResponses"][0][1]["ids"][0]

email_create = {
    "mailboxIds": {drafts_id: True},
    "from": [{"email": from_email, "name": "Anton"}],
    "to": [{"email": args.to, "name": "Willem"}],
    "subject": args.subject,
    "bodyStructure": {
        "type": "multipart/mixed",
        "subParts": [
            {"type": "text/plain", "partId": "body"},
            {
                "type": "application/pdf",
                "blobId": blob_id,
                "name": "SSWD_EvoEpi_Wavefront_Calibration_R2_W05-W16.pdf",
                "disposition": "attachment",
                "size": blob_size,
            }
        ]
    },
    "bodyValues": {
        "body": {"value": body_text, "isEncodingProblem": False, "isTruncated": False}
    },
}

resp = requests.post(jmap_url, headers=headers_json, json={
    "using": ["urn:ietf:params:jmap:core", "urn:ietf:params:jmap:mail", "urn:ietf:params:jmap:submission"],
    "methodCalls": [
        ["Email/set", {"accountId": account_id, "create": {"draft1": email_create}}, "1"],
        ["EmailSubmission/set", {
            "accountId": account_id,
            "create": {"sub1": {"emailId": "#draft1", "identityId": identity_id}}
        }, "2"]
    ]
})
resp.raise_for_status()
result = resp.json()

email_resp = result["methodResponses"][0][1]
submit_resp = result["methodResponses"][1][1]

if email_resp.get("created") and "draft1" in email_resp["created"]:
    print(f"Email created: {email_resp['created']['draft1']['id']}")
else:
    print(f"Email creation failed: {json.dumps(email_resp.get('notCreated', email_resp), indent=2)}")
    sys.exit(1)

if submit_resp.get("created") and "sub1" in submit_resp["created"]:
    print(f"Email sent to {args.to}!")
    print(f"Subject: {args.subject}")
else:
    print(f"Submission failed: {json.dumps(submit_resp.get('notCreated', submit_resp), indent=2)}")
    sys.exit(1)

print("Done!")
