#!/usr/bin/env python3
"""Send W75-W84 calibration report via Fastmail JMAP with PDF attachment."""

import json
import os
import sys
import base64
import requests

API_TOKEN = os.environ.get("FASTMAIL_API_TOKEN", "")
JMAP_SESSION_URL = "https://api.fastmail.com/jmap/session"
IDENTITY_ID = "172611807"
PDF_PATH = "/home/starbot/.openclaw/workspace/sswd-evoepi/reports/khalf_sweep_w75_w84/main.pdf"

def get_session():
    r = requests.get(JMAP_SESSION_URL, headers={
        "Authorization": f"Bearer {API_TOKEN}",
        "Content-Type": "application/json"
    })
    r.raise_for_status()
    session = r.json()
    api_url = session["apiUrl"]
    upload_url = session["uploadUrl"]
    account_id = session["primaryAccounts"]["urn:ietf:params:jmap:mail"]
    return api_url, upload_url, account_id

def upload_blob(upload_url, account_id, file_path):
    """Upload file as blob to JMAP."""
    url = upload_url.replace("{accountId}", account_id)
    with open(file_path, "rb") as f:
        data = f.read()
    r = requests.post(url, headers={
        "Authorization": f"Bearer {API_TOKEN}",
        "Content-Type": "application/pdf"
    }, data=data)
    r.raise_for_status()
    return r.json()["blobId"]

def get_drafts_mailbox(api_url, account_id):
    """Get Drafts mailbox ID."""
    resp = requests.post(api_url, headers={
        "Authorization": f"Bearer {API_TOKEN}",
        "Content-Type": "application/json"
    }, json={
        "using": ["urn:ietf:params:jmap:core", "urn:ietf:params:jmap:mail"],
        "methodCalls": [
            ["Mailbox/query", {"accountId": account_id, "filter": {"role": "drafts"}}, "0"]
        ]
    })
    resp.raise_for_status()
    return resp.json()["methodResponses"][0][1]["ids"][0]

def send_email(api_url, account_id, blob_id):
    """Create email with attachment and send."""
    drafts_id = get_drafts_mailbox(api_url, account_id)
    email_body = """Hi Willem,

Please find attached the calibration report for runs W75-W84 (K½/P_env_max sweep with dynamic P_env).

Key findings:
• P_env_max has ZERO effect — dynamic P_env fully dominates (W75=W78=W81, W76=W79=W82, etc.)
• K½ controls overall recovery level uniformly — no gradient improvement
• W76 (K½=800K): AK-PWS=57.1% ✅, AK-FS=21.2% ✅, but SS-S=44.1%, OR=36.8%, CA-N=34.9%
• Gradient collapsed from 35× (W71) to 1.6× (W76)
• The 10-run sweep effectively has only 4 distinct outcomes

The fundamental issue: uniform K½ cannot produce both correct AK recovery AND southern suppression. Need latitude-dependent K½ or stronger floor mechanism.

Best,
Anton (via OpenClaw)
"""

    resp = requests.post(api_url, headers={
        "Authorization": f"Bearer {API_TOKEN}",
        "Content-Type": "application/json"
    }, json={
        "using": [
            "urn:ietf:params:jmap:core",
            "urn:ietf:params:jmap:mail",
            "urn:ietf:params:jmap:submission"
        ],
        "methodCalls": [
            ["Email/set", {
                "accountId": account_id,
                "create": {
                    "draft1": {
                        "from": [{"name": "Anton", "email": "antonstar@fastmail.com"}],
                        "to": [{"name": "Willem", "email": "wlweert@gmail.com"}],
                        "subject": "SSWD Calibration Report: W75-W84 K½/P_env_max Sweep",
                        "textBody": [{"partId": "body", "type": "text/plain"}],
                        "attachments": [{
                            "blobId": blob_id,
                            "type": "application/pdf",
                            "name": "W75-W84_Calibration_Report.pdf"
                        }],
                        "bodyValues": {
                            "body": {"value": email_body}
                        },
                        "mailboxIds": {drafts_id: True},
                        "keywords": {"$draft": True}
                    }
                }
            }, "0"],
            ["EmailSubmission/set", {
                "accountId": account_id,
                "onSuccessUpdateEmail": {
                    "#sendIt": {
                        "mailboxIds/drafts": None,
                        "keywords/$draft": None
                    }
                },
                "create": {
                    "sendIt": {
                        "emailId": "#draft1",
                        "identityId": IDENTITY_ID
                    }
                }
            }, "1"]
        ]
    })
    resp.raise_for_status()
    result = resp.json()
    
    # Check for errors
    for method_resp in result["methodResponses"]:
        if "error" in str(method_resp[1]).lower():
            print(f"Warning: {method_resp}")
    
    return result

def main():
    if not API_TOKEN:
        print("ERROR: FASTMAIL_API_TOKEN not set")
        sys.exit(1)
    
    if not os.path.exists(PDF_PATH):
        print(f"ERROR: PDF not found at {PDF_PATH}")
        sys.exit(1)
    
    print(f"PDF size: {os.path.getsize(PDF_PATH)} bytes")
    
    print("Getting JMAP session...")
    api_url, upload_url, account_id = get_session()
    print(f"Account ID: {account_id}")
    
    print("Uploading PDF blob...")
    blob_id = upload_blob(upload_url, account_id, PDF_PATH)
    print(f"Blob ID: {blob_id}")
    
    print("Sending email...")
    result = send_email(api_url, account_id, blob_id)
    
    # Check submission result
    submissions = result["methodResponses"]
    for s in submissions:
        if s[0] == "EmailSubmission/set":
            created = s[1].get("created", {})
            if created:
                print("✓ Email sent successfully!")
            else:
                not_created = s[1].get("notCreated", {})
                print(f"✗ Email failed: {json.dumps(not_created, indent=2)}")
        elif s[0] == "Email/set":
            created = s[1].get("created", {})
            if created:
                print(f"✓ Draft created: {list(created.keys())}")
            else:
                not_created = s[1].get("notCreated", {})
                print(f"✗ Draft failed: {json.dumps(not_created, indent=2)}")

if __name__ == "__main__":
    main()
