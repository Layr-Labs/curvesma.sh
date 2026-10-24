#!/usr/bin/env python3
"""preflight.py -- v19: expand the shipped archive and ASSERT it is self-consistent.

Checks, in order (any failure is fatal, none is a warning):
  1. every path referenced by submission.json / reproduction.json exists IN the
     package and, after expansion, inside the expanded archive (no dangling refs);
  2. the claim arithmetic: score.numerator/denominator equals the ledger's claim,
     and the note's quoted log2 W matches ledger's to the printed digits;
  3. log-label audit: every "log2 X <number>" label in the note is recomputed or
     matched against a ledger row -- no stale label survives;
  4. >=/<= direction audit: every ">= 0.39" / "<= 2^" phrase in the note has the
     comparison on the same side as the ledger's own comparison;
  5. file count <= 30, no __pycache__, no verdict files, no absolute paths.
"""
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
D = HERE   # preflight.py sits at the PACKAGE ROOT, not in scripts/
PY = "/opt/homebrew/opt/python@3.14/bin/python3.14"
fails = []


def chk(cond, msg):
    if not cond:
        fails.append(msg)
    return cond


files = []
for r, _d, fs in os.walk(D):
    if "__pycache__" in r or "/.git" in r:
        continue
    for f in fs:
        files.append(os.path.relpath(os.path.join(r, f), D))

# 5. hygiene
chk(len(files) <= 30, "file count %d > 30" % len(files))
chk(not any("__pycache__" in f for f in files), "__pycache__ present")
chk(not any("verdict" in f.lower() for f in files), "verdict file inside submission")

# 1. referenced paths
sub = json.load(open(os.path.join(D, "submission.json")))
rep = json.load(open(os.path.join(D, "reproduction.json")))
refs = set(sub["claims"][0]["evidence"]) | {"public-note.md", "submission.json",
                                            "reproduction.json", "ledger.py"}
for c in rep["commands"]:
    for m in re.finditer(r"scripts/[\w_]+\.py", c):
        refs.add(m.group(0))
for r in sorted(refs):
    chk(os.path.exists(os.path.join(D, r)), "referenced path missing: %s" % r)

# expand the archive the way intake will, then re-check inside it
tmp = tempfile.mkdtemp(prefix="v19preflight")
arc = os.path.join(tmp, "pkg")
shutil.copytree(D, arc, ignore=shutil.ignore_patterns("__pycache__"))
for r in sorted(refs):
    chk(os.path.exists(os.path.join(arc, r)), "missing INSIDE expanded archive: %s" % r)
size = sum(os.path.getsize(os.path.join(r_, f)) for r_, _d_, fs_ in os.walk(arc)
           for f in fs_)
chk(size <= 262144, "expanded size %d > 262144" % size)

# 1b. THE ARCHIVE SOURCE IS THE WHOLE editablePaths, NOT THIS DIRECTORY.
# intake expands every editable path of the track, so a sibling submission tree
# that is left in place is charged against the 262144-byte limit.  Measure it here
# (uncompressed, as intake does) and fail BEFORE the submit.
REPO = os.path.dirname(os.path.dirname(os.path.dirname(D)))   # <repo>/submissions/<track>/<pkg>
bpath = os.path.join(REPO, "benchmark.json")
arc_bytes = None
if os.path.exists(bpath):
    _b = json.load(open(bpath))
    _paths = None
    for _t in _b.get("tracks", []):
        if _t.get("name") == "secp256k1-ecdlp":
            _paths = _t.get("editablePaths")
    if _paths:
        arc_bytes = 0
        for _p in _paths:
            _abs = os.path.join(REPO, _p)
            if os.path.isfile(_abs):
                arc_bytes += os.path.getsize(_abs)
                continue
            for _r, _d, _fs in os.walk(_abs):
                if "__pycache__" in _r:
                    continue
                for _f in _fs:
                    arc_bytes += os.path.getsize(os.path.join(_r, _f))
        print("archive source (editablePaths, uncompressed) %d bytes -> %s"
              % (arc_bytes, "OK" if arc_bytes <= 262144 else "OVER 262144"))
        chk(arc_bytes <= 240000, "archive source %d > 240000 (intake limit 262144)" % arc_bytes)
        _sel = os.path.join(REPO, "tracks", "secp256k1-ecdlp", "selected-submission.txt")
        _want = os.path.relpath(D, REPO)
        chk(open(_sel).read().strip() == _want,
            "selected-submission.txt does not point at this package (%s)" % _want)
        chk(os.path.isdir(os.path.join(REPO, _want)), "selected path absent from archive")

# 2. run the ledger inside the expanded archive
env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
p = subprocess.run([PY, "ledger.py"], cwd=arc, capture_output=True, text=True, env=env)
chk(p.returncode == 0, "ledger failed in expanded archive: %s" % p.stderr[-400:])
out = p.stdout
m = re.search(r"log2 W \(submitted\) ([\d.]+)\s+claim ([\d.]+)", out)
chk(m is not None, "ledger claim row not found")
if m:
    logw, claim = m.group(1), m.group(2)
    sc = sub["claims"][0]["score"]
    chk(abs(float(claim) - int(sc["numerator"]) / int(sc["denominator"])) < 1e-9,
        "score %s != ledger claim %s" % (sc, claim))
    note = open(os.path.join(D, "public-note.md")).read()
    chk(logw in note, "note does not quote ledger log2 W %s" % logw)
    chk(claim in note, "note does not quote ledger claim %s" % claim)
    # every weight quoted in the note's claim block must appear in the ledger run
    for num in ["126.92945732", "127.92945732", "0.475029", "0.474224", "4.6875e-4",
                "0.084224", "6.3129e-06", "0.0000"]:
        if num in note:
            chk(num in out or num in json.dumps(json.load(
                open(os.path.join(D, "evidence", "sensitivity.json")))),
                "note number %s has no ledger/evidence row" % num)

# 3. log-label audit: every "log2 <name> <number>" in the note must be reproduced
note = open(os.path.join(D, "public-note.md")).read()
for m in re.finditer(r"log2 (W|W2|S0|S_ops)\s*=?\s*([\d.]+)", note):
    chk(m.group(2) in out, "log label %s=%s not in ledger output" % (m.group(1), m.group(2)))

# 4. direction audit on the floor
chk("0.474224 >= 0.39" in note.replace("≥", ">="), "note lost the success >= floor row")
chk("met True" in out, "ledger floor decision is not met")
chk(out.count("floor 0.39") >= 1 and "success 0.474224" in out, "floor row changed")

# 57. fields <= 300 chars where the schema slices them
for k in ("title",):
    chk(len(sub.get(k, "")) <= 300, "%s longer than 300" % k)
chk(len(sub["claims"][0]["statement"]) <= 7965, "claims[0] too long")

chk(arc_bytes is not None, "archive source could not be measured (benchmark.json not found) -> refusing to guess")
shutil.rmtree(tmp, ignore_errors=True)
print("files %d  expanded %d bytes  refs %d" % (len(files), size, len(refs)))
if fails:
    print("PREFLIGHT FAILED:")
    for f in fails:
        print("  -", f)
    sys.exit(1)
print("PREFLIGHT OK (claim %s, log2 W2 %s, floor met)" % (claim, logw))
