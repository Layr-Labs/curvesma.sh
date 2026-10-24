#!/usr/bin/env python3
"""count_symbols.py -- instruction count of every label of the shipped listing."""
import json, os, re, sys
HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "walk_v15o.s")
out, cur = {}, None
for ln in open(path):
    m = re.match(r'^(_[A-Za-z0-9_.]+):', ln)
    if m:
        cur = m.group(1); out[cur] = 0; continue
    if cur and re.match(r'^\s+[a-z][a-z0-9.]*\s', ln) and not ln.strip().startswith('.'):
        out[cur] += 1
print(json.dumps({"kind": "counts-v15o", "listing": os.path.relpath(path, PKG),
                  "method": "instruction lines per labelled symbol", "symbols": out}, indent=1))
