#!/usr/bin/env python3
"""Enhanced JMAP client with attachment support for sending the calibration report."""

import json
import os
import sys
import requests
import base64
import mimetypes
from datetime import datetime, timezone
from pathlib import Path

API_TOKEN = os.environ.get("FASTMAIL_API_TOKEN", "")
JMAP_SESSION_URL = "https://api.fastmail.com/jmap/session"

def get_session():
    """Get JMAP session and return (api_url, account_id)."""
    r = requests.get(JMAP_SESSION_URL, headers={
        "Authorization": f"Bearer {API_TOKEN}",
        "Content-Type": "application/json"
    })
    r.raise_for_status()
    session = r.json()
    api_url = session["apiUrl"]
    account_id = session["primaryAccounts"]["urn:ietf:params:jmap:mail"]
    return api_url, account_id

def jmap_call(api_url, method_calls):
    """Make a JMAP API call."""
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

def upload_attachment(api_url, account_id, file_path):
    """Upload a file and return the blob_id for attachment."""
    file_path = Path(file_path)
    
    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")
    
    # Get upload URL
    session_resp = requests.get(JMAP_SESSION_URL, headers={
        "Authorization": f"Bearer {API_TOKEN}",
        "Content-Type": "application/json"
    })
    session_resp.raise_for_status()
    upload_url = session_resp.json()["uploadUrl"]
    
    # Upload the file
    with open(file_path, 'rb') as f:
        upload_resp = requests.post(
            upload_url.replace("{accountId}", account_id),
            headers={"Authorization": f"Bearer {API_TOKEN}"},
            files={"file": (file_path.name, f, mimetypes.guess_type(str(file_path))[0] or 'application/octet-stream')}
        )
    
    upload_resp.raise_for_status()
    upload_result = upload_resp.json()
    
    return upload_result["blobId"], upload_result["size"], upload_result["type"]

def send_email_with_attachment(to_addr, subject, body, attachment_path):
    """Send an email with attachment."""
    api_url, account_id = get_session()
    
    # Upload attachment
    print(f"Uploading attachment: {attachment_path}")
    blob_id, file_size, content_type = upload_attachment(api_url, account_id, attachment_path)
    print(f"Attachment uploaded: blob_id={blob_id}, size={file_size}, type={content_type}")
    
    # Get drafts mailbox
    resp = jmap_call(api_url, [
        ["Mailbox/query", {"accountId": account_id, "filter": {"role": "drafts"}}, "0"]
    ])
    drafts_id = resp["methodResponses"][0][1]["ids"][0]
    
    # Prepare email with attachment
    attachment_filename = Path(attachment_path).name
    
    email_body = {
        "accountId": account_id,
        "create": {
            "draft1": {
                "mailboxIds": {drafts_id: True},
                "from": [{"email": "antonstar@fastmail.com"}],
                "to": [{"email": to_addr}],
                "subject": subject,
                "bodyStructure": {
                    "type": "multipart/mixed",
                    "subParts": [
                        {
                            "partId": "text",
                            "type": "text/plain"
                        },
                        {
                            "partId": "attachment1",
                            "type": content_type,
                            "disposition": "attachment",
                            "name": attachment_filename,
                            "blobId": blob_id
                        }
                    ]
                },
                "textBody": [{"partId": "text", "type": "text/plain"}],
                "bodyValues": {
                    "text": {"value": body}
                },
                "keywords": {"$draft": True}
            }
        }
    }
    
    # Get sending identity
    id_resp = jmap_call(api_url, [
        ["Identity/get", {"accountId": account_id}, "0"]
    ])
    identities = id_resp["methodResponses"][0][1].get("list", [])
    identity_id = identities[0]["id"] if identities else None
    
    # Create email and submit
    submission = {
        "emailId": "#draft1",
        "envelope": {
            "mailFrom": {"email": "antonstar@fastmail.com"},
            "rcptTo": [{"email": to_addr}]
        }
    }
    if identity_id:
        submission["identityId"] = identity_id
    
    print("Sending email...")
    resp = jmap_call(api_url, [
        ["Email/set", email_body, "0"],
        ["EmailSubmission/set", {
            "accountId": account_id,
            "create": {"sub1": submission},
            "onSuccessUpdateEmail": {
                "#sub1": {
                    f"mailboxIds/{drafts_id}": None,
                    "keywords/$draft": None
                }
            }
        }, "1"]
    ])
    
    # Check for success
    print("Response:", json.dumps(resp, indent=2))
    responses = resp["methodResponses"]
    
    email_created = False
    submission_created = False
    
    for i, r in enumerate(responses):
        print(f"Response {i}: method={r[0]}, data type={type(r[1])}")
        if len(r) > 1 and r[1] is not None:
            if r[0] == "Email/set":
                if r[1].get("created", {}).get("draft1"):
                    print("✅ Email draft created")
                    email_created = True
                elif r[1].get("notCreated"):
                    print(f"❌ Email creation error: {json.dumps(r[1]['notCreated'], indent=2)}")
                    return False
            elif r[0] == "EmailSubmission/set":
                if r[1].get("created", {}).get("sub1"):
                    print(f"✅ Email submission created")
                    submission_created = True
                elif r[1].get("notCreated"):
                    print(f"❌ Email submission error: {json.dumps(r[1]['notCreated'], indent=2)}")
                    return False
        else:
            print(f"Response {i}: No data or null response")
    
    if email_created and submission_created:
        print(f"✅ Email sent successfully to {to_addr}")
        return True
    elif email_created:
        print(f"✅ Email created but submission unclear - likely sent to {to_addr}")
        return True
    else:
        print(f"❌ Email sending status unclear")
        return False

def main():
    """Send the calibration report via email."""
    # Email details
    to_email = "wlweert@gmail.com"
    subject = "SSWD-EvoEpi: W105-W114 Calibration Report (Multi-Trait Pathogen Evolution)"
    
    body = """Willem,

I've completed the comprehensive calibration analysis for W105-W114, exploring multi-trait pathogen community evolution in the SSWD-EvoEpi model.

Key findings:
• Best overall performance: W113 (control, no vir_evo) with RMSE 0.645
• Best virulence evolution run: W106 with RMSE 0.644  
• Realistic spatial gradients emerged: v_local 0.07→0.28 (AK→CA), T_vbnc 6.6→11.5°C
• Fast adaptation (v_adapt=0.005) was counterproductive - flattened gradients
• Alaska challenge persists: ~8% vs 50% target across all configs
• All regions still declining at year 13 - may need different recovery mechanisms

The attached PDF contains the complete analysis with 6 figures, detailed methods, results, and recommendations for the next calibration round (including raising T_vbnc_min to 9°C per literature).

Best regards,
Claude
"""
    
    attachment_path = "/home/starbot/.openclaw/workspace/sswd-evoepi/reports/virulence_evo_w105_w114/main.pdf"
    
    # Verify attachment exists
    if not Path(attachment_path).exists():
        print(f"❌ Error: PDF not found at {attachment_path}")
        return False
    
    # Check API token
    if not API_TOKEN:
        print("❌ Error: FASTMAIL_API_TOKEN not set")
        return False
    
    # Send email
    try:
        success = send_email_with_attachment(to_email, subject, body, attachment_path)
        return success
    except Exception as e:
        print(f"❌ Error sending email: {e}")
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)