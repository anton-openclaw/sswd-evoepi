#!/usr/bin/env python3
"""Send the wavefront calibration report via Fastmail JMAP."""

import os, json, base64, requests, sys

# Configuration
token = os.environ.get("FASTMAIL_API_TOKEN")
if not token:
    print("ERROR: FASTMAIL_API_TOKEN not set")
    sys.exit(1)

headers_json = {"Authorization": f"Bearer {token}", "Content-Type": "application/json"}

# Discover JMAP session
session_resp = requests.get("https://api.fastmail.com/.well-known/jmap", 
                            headers={"Authorization": f"Bearer {token}"})
session_resp.raise_for_status()
session = session_resp.json()
jmap_url = session["apiUrl"]
upload_url_template = session["uploadUrl"]
account_id = list(session["accounts"].keys())[0]
upload_url = upload_url_template.replace("{accountId}", account_id)
print(f"Account: {account_id}")
print(f"API URL: {jmap_url}")
print(f"Upload URL: {upload_url}")

identity_id = "172611807"

pdf_path = "/home/starbot/.openclaw/workspace/sswd-evoepi/reports/wavefront_calibration/main.pdf"
to_email = "wlweert@gmail.com"
from_email = "antonstar@fastmail.com"
subject = "SSWD-EvoEpi: Wavefront Calibration Report (W01-W04)"

body_text = """Hi Willem,

Attached is the wavefront calibration report for rounds W01-W04 of the SSWD-EvoEpi model.

Key findings:
- Wavefront propagation works but is ~2x too fast (all regions reached in 12-20 months vs. target 6-42 months)
- D_P (50-200km) has surprisingly little effect on timing or recovery
- Recovery gradient remains flat (~21-46%) — same issue as pre-wavefront runs
- AK-PWS recovery is close to target (~46% vs 50%)
- RMSE nearly identical across rounds (~1.19)

Recommended next steps: increase activation_threshold to slow wavefront, explore T_VBNC for recovery gradient.

The report includes 6 figures with detailed analysis and data tables.

Best,
Anton (via OpenClaw)
"""

# Step 1: Upload the PDF blob
print("Uploading PDF...")
with open(pdf_path, "rb") as f:
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

# Step 2: Create the email with attachment and submit it
print("Creating and sending email...")

email_create = {
    "mailboxIds": {},  # Will be set to Drafts below
    "from": [{"email": from_email, "name": "Anton"}],
    "to": [{"email": to_email, "name": "Willem"}],
    "subject": subject,
    "bodyStructure": {
        "type": "multipart/mixed",
        "subParts": [
            {
                "type": "text/plain",
                "partId": "body",
            },
            {
                "type": "application/pdf",
                "blobId": blob_id,
                "name": "SSWD_EvoEpi_Wavefront_Calibration_W01-W04.pdf",
                "disposition": "attachment",
                "size": blob_size,
            }
        ]
    },
    "bodyValues": {
        "body": {
            "value": body_text,
            "isEncodingProblem": False,
            "isTruncated": False,
        }
    },
}

# First, get the Drafts mailbox
resp = requests.post(jmap_url, headers=headers_json, json={
    "using": ["urn:ietf:params:jmap:core", "urn:ietf:params:jmap:mail"],
    "methodCalls": [["Mailbox/query", {"accountId": account_id, "filter": {"role": "drafts"}}, "0"]]
})
resp.raise_for_status()
drafts_id = resp.json()["methodResponses"][0][1]["ids"][0]
print(f"Drafts mailbox: {drafts_id}")

email_create["mailboxIds"] = {drafts_id: True}

# Create email and submit in one call
resp = requests.post(jmap_url, headers=headers_json, json={
    "using": ["urn:ietf:params:jmap:core", "urn:ietf:params:jmap:mail", "urn:ietf:params:jmap:submission"],
    "methodCalls": [
        ["Email/set", {
            "accountId": account_id,
            "create": {"draft1": email_create}
        }, "1"],
        ["EmailSubmission/set", {
            "accountId": account_id,
            "create": {
                "sub1": {
                    "emailId": "#draft1",
                    "identityId": identity_id,
                }
            }
        }, "2"]
    ]
})
resp.raise_for_status()
result = resp.json()

# Debug: print full response
print("Response:")
print(json.dumps(result["methodResponses"], indent=2)[:3000])

email_resp = result["methodResponses"][0][1]
submit_resp = result["methodResponses"][1][1]

if email_resp.get("created") and "draft1" in email_resp["created"]:
    email_id = email_resp["created"]["draft1"]["id"]
    print(f"\nEmail created: {email_id}")
elif email_resp.get("notCreated"):
    print(f"\nEmail creation failed: {json.dumps(email_resp['notCreated'], indent=2)}")
    sys.exit(1)
else:
    print(f"\nEmail response unexpected: {json.dumps(email_resp, indent=2)}")
    sys.exit(1)

if submit_resp.get("created") and "sub1" in submit_resp["created"]:
    print(f"Email submitted successfully!")
    print(f"Sent to: {to_email}")
    print(f"Subject: {subject}")
elif submit_resp.get("notCreated"):
    print(f"Email submission failed: {json.dumps(submit_resp['notCreated'], indent=2)}")
    sys.exit(1)
else:
    print(f"Submission response unexpected: {json.dumps(submit_resp, indent=2)}")
    sys.exit(1)

print("Done!")
