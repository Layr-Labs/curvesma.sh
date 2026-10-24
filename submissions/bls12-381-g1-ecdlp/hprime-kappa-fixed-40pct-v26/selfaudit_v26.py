#!/usr/bin/env python3
"""v26 self-audit: laundry list, no-deferred grep, selector==staged path."""
import json, os, re, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
FORBID = re.compile(r"\b(deferred|TBD|to be derived|defer to v24|deferred to v24)\b", re.I)
hits = []
for root, dirs, files in os.walk(HERE):
    dirs[:] = [d for d in dirs if d not in ("__pycache__", ".git")]
    for f in files:
        p = os.path.join(root, f)
        if f == "selfaudit_v26.py":
            continue
        if f.endswith((".md", ".json", ".py", ".txt")):
            try:
                for i, line in enumerate(open(p, errors="ignore"), 1):
                    if FORBID.search(line):
                        hits.append((os.path.relpath(p, HERE), i, line.strip()[:100]))
            except OSError:
                pass
sub = json.load(open(os.path.join(HERE, "submission.json")))
pkg_bytes = 0
nfiles = 0
for root, dirs, files in os.walk(HERE):
    dirs[:] = [d for d in dirs if d not in ("__pycache__", ".git")]
    for f in files:
        nfiles += 1
        pkg_bytes += os.path.getsize(os.path.join(root, f))
rel = os.path.basename(HERE)
out = dict(
    forbidden_strings=hits,
    no_deferred=not hits,
    n_files=nfiles,
    package_bytes=pkg_bytes,
    under_30_files=bool(nfiles <= 30),
    under_256KiB=bool(pkg_bytes <= 262144),
    selector=rel,
    selector_matches_staged=bool(rel == "hprime-kappa-fixed-40pct-v26"),
    claim_text_chars=len(sub["claims"][0]["statement"]),
    claim_text_cap_ok=bool(len(sub["claims"][0]["statement"]) <= 7965),
    submission_keys=sorted(sub),
    submission_minimal=bool(set(sub) == {"version", "title", "claims", "experiments"}),
    evidence_count=len(sub["claims"][0]["evidence"]),
    all_evidence_present=all(os.path.exists(os.path.join(HERE, e)) for e in sub["claims"][0]["evidence"]),
)
print(json.dumps(out, indent=1))
ok = out["no_deferred"] and out["under_30_files"] and out["under_256KiB"] and out["claim_text_cap_ok"] and out["submission_minimal"] and out["all_evidence_present"]
raise SystemExit(0 if ok else 1)
