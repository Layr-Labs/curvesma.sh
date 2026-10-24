#!/usr/bin/env python3
"""manifest32.py -- evidence/manifest32.json: sha256 + byte count for every package
file.  The manifest cannot contain its own hash, so it lists every file EXCEPT
itself; cert32.py checks that set == the tree minus the manifest."""
import hashlib
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
OUT = os.path.join(PKG, "evidence", "manifest32.json")
SELF = os.path.relpath(OUT, PKG)

files = []
for r, ds, fs in os.walk(PKG):
    ds[:] = [d for d in ds if d != "__pycache__"]
    for f in fs:
        if f.startswith(".") or f == "__pycache__":
            continue
        p = os.path.join(r, f)
        rel = os.path.relpath(p, PKG)
        if rel == SELF:
            continue
        b = open(p, "rb").read()
        files.append({"path": rel, "bytes": len(b), "sha256": hashlib.sha256(b).hexdigest()})
files.sort(key=lambda x: x["path"])
json.dump({"files": files, "count": len(files), "excludes": SELF},
          open(OUT, "w"), indent=1)
print("manifest files", len(files), "total bytes", sum(f["bytes"] for f in files))
print("max file", max(f["bytes"] for f in files))
