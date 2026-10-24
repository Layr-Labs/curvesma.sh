#!/usr/bin/env python3
"""v27 self-audit: intake-law v3 laundry, no-deferred grep, selector==staged path,
title/claim caps, package bytes, 30-file cap, no bytecode."""
import json, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
FORBID = re.compile(r"\b(deferred|TBD|to be derived|defer to v24|deferred to v24)\b", re.I)
hits = []
pkg_bytes = 0
nfiles = 0
pycache = []
for root, dirs, files in os.walk(HERE):
    dirs[:] = [d for d in dirs if d not in (".git",)]
    for d in dirs:
        if d == "__pycache__":
            pycache.append(os.path.relpath(os.path.join(root, d), HERE))
    for f in files:
        p = os.path.join(root, f)
        nfiles += 1
        pkg_bytes += os.path.getsize(p)
        if f.endswith(".pyc"):
            pycache.append(os.path.relpath(p, HERE))
        if f == "selfaudit_v27.py":
            continue
        if f.endswith((".md", ".json", ".py", ".txt")):
            try:
                for i, line in enumerate(open(p, errors="ignore"), 1):
                    if FORBID.search(line):
                        hits.append((os.path.relpath(p, HERE), i, line.strip()[:100]))
            except OSError:
                pass
sub = json.load(open(os.path.join(HERE, "submission.json")))
ev = sub["claims"][0]["evidence"]
rel = os.path.basename(HERE)
out = dict(
    forbidden_strings=hits,
    no_deferred=not hits,
    n_files=nfiles,
    package_bytes=pkg_bytes,
    under_30_files=bool(nfiles <= 30),
    under_262144=bool(pkg_bytes <= 262144),
    under_240KiB=bool(pkg_bytes <= 245760),
    no_bytecode=not pycache,
    pycache=pycache,
    selector=rel,
    selector_matches_staged=bool(rel == "hprime-derived-40pct-v27"),
    title_chars=len(sub["title"]),
    title_cap_ok=bool(len(sub["title"]) <= 300),
    claim_text_chars=len(sub["claims"][0]["statement"]),
    claim_text_cap_ok=bool(len(sub["claims"][0]["statement"]) <= 7965),
    submission_keys=sorted(sub),
    submission_minimal=bool(set(sub) == {"version", "title", "claims", "experiments"}),
    evidence_count=len(ev),
    all_evidence_present=all(os.path.exists(os.path.join(HERE, e)) for e in ev),
)
print(json.dumps(out, indent=1))
ok = (out["no_deferred"] and out["under_30_files"] and out["under_262144"]
      and out["no_bytecode"] and out["claim_text_cap_ok"] and out["submission_minimal"]
      and out["all_evidence_present"] and out["title_cap_ok"])
raise SystemExit(0 if ok else 1)
