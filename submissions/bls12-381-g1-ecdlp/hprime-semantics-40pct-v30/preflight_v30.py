#!/usr/bin/env python3
"""v30 preflight (INTAKE LAW v4): expand our own payload, assert every path named by
tracks/<track>/selected-submission.txt, submission.json and claims[0] exists INSIDE the intake root,
assert no '..'/symlinks, byte/file caps, verdicts excluded, __pycache__ pruned, experiments/ present.
Exit 0 iff every assertion holds."""
import json, os, re, sys

ROOT = "/Users/andrewgordon/Documents/ChatGPT/curvesma.sh/curvesmash"
TRACK = "bls12-381-g1-ecdlp"
SEL = ROOT + "/tracks/" + TRACK + "/selected-submission.txt"
BYTES_CAP, FILE_CAP = 262144, 30
fails = []


def chk(cond, msg):
    if not cond:
        fails.append(msg)


sel = open(SEL).read().strip()
print("selected-submission.txt ->", sel)
chk(sel.startswith("submissions/" + TRACK + "/"), "selector not under submissions/<track>/: %r" % sel)
chk("archive" not in sel, "selector points into an archive path: %r" % sel)
sub_dir = os.path.join(ROOT, sel)
chk(os.path.isdir(sub_dir), "selector names a missing directory INSIDE the archive: %s" % sub_dir)

sub = json.load(open(os.path.join(sub_dir, "submission.json")))
claim = sub["claims"][0]
chk(len(claim["statement"]) <= 7965, "claims[0] > 7965 chars: %d" % len(claim["statement"]))
chk(len(sub.get("title", "")) <= 300, "title > 300 chars: %d" % len(sub.get("title", "")))
chk(sub.get("version") == 1 and sub.get("experiments") == [], "submission.json not minimal")
chk(os.path.isdir(os.path.join(sub_dir, "experiments")), "experiments/ missing (empty dirs vanish in tar)")

for e in claim["evidence"]:
    p = os.path.join(sub_dir, e)
    chk(os.path.isfile(p), "evidence named but missing inside the payload: %s" % e)

files, total = [], 0
for dirpath, dirnames, filenames in os.walk(sub_dir):
    for d in list(dirnames):
        full = os.path.join(dirpath, d)
        chk(not os.path.islink(full), "symlink in payload: %s" % full)
        chk(d != "__pycache__", "__pycache__ in payload: %s" % full)
    for f in filenames:
        full = os.path.join(dirpath, f)
        chk(not os.path.islink(full), "symlink in payload: %s" % full)
        chk(".." not in os.path.relpath(full, sub_dir).split(os.sep), "parent ref: %s" % full)
        chk("verdict" not in f.lower(), "verdict file left in payload: %s" % f)
        files.append(os.path.relpath(full, sub_dir))
        total += os.path.getsize(full)
chk(len(files) <= FILE_CAP, "too many files: %d > %d" % (len(files), FILE_CAP))
chk(total <= BYTES_CAP, "payload too big: %d > %d" % (total, BYTES_CAP))

named = set(re.findall(r"[A-Za-z0-9_.\-]+\.(?:json|py|c|md)", claim["statement"]))
for n in sorted(named):
    chk(n in files, "claims[0] names %s but it is not in the payload" % n)

print("files=%d bytes=%d cap_ok=%s experiments=%s"
      % (len(files), total, total <= BYTES_CAP, os.path.isdir(os.path.join(sub_dir, "experiments"))))
if fails:
    print("PREFLIGHT FAIL")
    for f in fails:
        print(" -", f)
    sys.exit(1)
print("PREFLIGHT PASS")
