#!/usr/bin/env python3
"""preflight.py -- v12: expand the shipped archive and ASSERT it is self-consistent.

Any failure is fatal.  Checks:
  1. every path referenced by submission.json / reproduction.json exists IN the package
     and, after expansion, INSIDE the expanded copy (no dangling refs);
  2. the ledger RUNS inside the expanded copy and its log2 W / claim rows match
     submission.json's score and appear verbatim in the public note;
  3. the ARCHIVE SOURCE (all editablePaths of the track, uncompressed, as intake
     measures it) is <= 262144 bytes, with a 240000 fail-safe, and
     tracks/<track>/selected-submission.txt points at THIS package (the intake kill
     that has taken three packages tonight);
  4. hygiene: <= 30 files, no __pycache__, no verdict files;
  5. field limits: title <= 300, claims[0].statement <= 7965.
"""
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
TRACK = "p256-ecdlp"
PY = "/opt/homebrew/opt/python@3.14/bin/python3.14"
fails = []


def chk(cond, msg):
    if not cond:
        fails.append(msg)
    return cond


files = []
for r, _d, fs in os.walk(HERE):
    if "__pycache__" in r or "/.git" in r:
        continue
    for f in fs:
        files.append(os.path.relpath(os.path.join(r, f), HERE))

chk(len(files) <= 30, "file count %d > 30" % len(files))
chk(not any("__pycache__" in f for f in files), "__pycache__ present")
chk(not any("verdict" in f.lower() for f in files), "verdict file inside submission")

sub = json.load(open(os.path.join(HERE, "submission.json")))
rep = json.load(open(os.path.join(HERE, "reproduction.json")))
refs = set(sub["claims"][0]["evidence"]) | {"public-note.md", "submission.json",
                                            "reproduction.json", "ledger_v12.py",
                                            "selfaudit.py", "assemble.py"}
for c in rep["commands"]:
    for m in re.finditer(r"[\w/]+\.py", c):
        refs.add(m.group(0))
for r in sorted(refs):
    chk(os.path.exists(os.path.join(HERE, r)), "referenced path missing: %s" % r)

tmp = tempfile.mkdtemp(prefix="v12preflight")
arc = os.path.join(tmp, "pkg")
shutil.copytree(HERE, arc, ignore=shutil.ignore_patterns("__pycache__"))
for r in sorted(refs):
    chk(os.path.exists(os.path.join(arc, r)), "missing INSIDE expanded archive: %s" % r)
size = sum(os.path.getsize(os.path.join(r_, f)) for r_, _d_, fs_ in os.walk(arc) for f in fs_)
chk(size <= 262144, "expanded size %d > 262144" % size)

REPO = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))  # <repo>/submissions/<track>/<pkg>
bpath = os.path.join(REPO, "benchmark.json")
arc_bytes = None
if os.path.exists(bpath):
    _b = json.load(open(bpath))
    paths = None
    for _t in _b.get("tracks", []):
        if _t.get("name") == TRACK:
            paths = _t.get("editablePaths")
    if paths:
        arc_bytes = 0
        for p_ in paths:
            a_ = os.path.join(REPO, p_)
            if os.path.isfile(a_):
                arc_bytes += os.path.getsize(a_)
                continue
            for _r, _d, _fs in os.walk(a_):
                if "__pycache__" in _r:
                    continue
                for _f in _fs:
                    arc_bytes += os.path.getsize(os.path.join(_r, _f))
        print("archive source (editablePaths, uncompressed) %d bytes -> %s"
              % (arc_bytes, "OK" if arc_bytes <= 240000 else "OVER 240000"))
        chk(arc_bytes <= 240000, "archive source %d > 240000 (intake 262144)" % arc_bytes)
        sel = os.path.join(REPO, "tracks", TRACK, "selected-submission.txt")
        want = os.path.relpath(HERE, REPO)
        chk(open(sel).read().strip() == want,
            "selected-submission.txt does not point at this package (%s)" % want)
        chk(os.path.isdir(os.path.join(REPO, want)), "selected path absent from archive")

env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
p = subprocess.run([PY, "ledger_v12.py"], cwd=arc, capture_output=True, text=True, env=env)
chk(p.returncode == 0, "ledger failed inside expanded archive: %s" % p.stderr[-400:])
out = p.stdout
m = re.search(r"log2 W\s+([\d.]+)\s+claim\s+([\d.]+)", out)
chk(m is not None, "ledger claim row not found")
if m:
    logw, claim = m.group(1), m.group(2)
    sc = sub["claims"][0]["score"]
    chk(abs(float(claim) - int(sc["numerator"]) / int(sc["denominator"])) < 1e-9,
        "score %s != ledger claim %s" % (sc, claim))
    note = open(os.path.join(HERE, "public-note.md")).read()
    chk(logw in note, "note does not quote ledger log2 W %s" % logw)
    chk(claim in note, "note does not quote ledger claim %s" % claim)
    chk("floor met True" in out, "ledger floor decision is not met")
    chk("NOT SUPPORTED" in out, "ledger lost the transfer disclosure row")
    for num in re.findall(r"\d+\.\d{6,}", note):
        chk(num in out or num in json.dumps(json.load(
            open(os.path.join(HERE, "evidence", "ledger-out.json")))),
            "note number %s has no ledger row" % num)

chk(len(sub.get("title", "")) <= 300, "title longer than 300")
chk(len(sub["claims"][0]["statement"]) <= 7965, "claims[0] too long")
chk(arc_bytes is not None, "archive source unmeasured -> refusing to guess")
shutil.rmtree(tmp, ignore_errors=True)
print("files %d  expanded %d bytes  refs %d" % (len(files), size, len(refs)))
if fails:
    print("PREFLIGHT FAILED:")
    for f in fails:
        print("  -", f)
    sys.exit(1)
print("PREFLIGHT OK (claim %s, log2 W %s, floor met, transfer disclosed)" % (claim, logw))
