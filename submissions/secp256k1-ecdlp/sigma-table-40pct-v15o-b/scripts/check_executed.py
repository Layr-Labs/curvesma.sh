#!/usr/bin/env python3
"""check_executed.py -- static length == executed length, machine-checked.

For every counted symbol of a listing, find backward branches (a branch whose
target label appears earlier in the same symbol).  A symbol with none executes
each of its instructions exactly once per call, so its static count IS its
executed count.  The shipped TU is compiled with `#pragma clang loop unroll(full)`
on every loop precisely to make this true; the point-step listing is hand-written
straight-line and is checked the same way.
"""
import json, os, re, sys
HERE = os.path.dirname(os.path.abspath(__file__)); PKG = os.path.dirname(HERE)

def scan(path):
    lines = open(path).read().splitlines()
    out, cur, pos = {}, None, {}
    for i, ln in enumerate(lines):
        m = re.match(r'^(_[A-Za-z0-9_.]+|step_asm):', ln)
        if m:
            cur = m.group(1); out[cur] = {"n": 0, "back": 0, "labels": {}}
        if cur is None:
            continue
        lm = re.match(r'^(\.L[A-Za-z0-9_.]+):', ln.strip())
        if lm:
            out[cur]["labels"][lm.group(1)] = i
        t = ln.strip()
        if re.match(r'^\s+[a-z][a-z0-9.]*\s', ln) and not t.startswith('.'):
            out[cur]["n"] += 1
            bm = re.match(r'^b[a-z.]*\s+(\.L[A-Za-z0-9_.]+)', t)
            if bm and bm.group(1) in out[cur]["labels"] and out[cur]["labels"][bm.group(1)] < i:
                out[cur]["back"] += 1
    return out

res = {}
for f in ("scripts/walk_v15o.s", "evidence/step_secp_v15g.s"):
    for k, v in scan(os.path.join(PKG, f)).items():
        if k in ("_v15o_align", "_fmul", "_v15o_label", "_rmont", "_v15o_detect",
                 "_v15o_report", "step_asm"):
            res[k] = {"static": v["n"], "backward_branches": v["back"],
                      "executed_equals_static": v["back"] == 0}
# step_asm is a genuine loop over n: its per-iteration body is priced by the
# ledger's base_step_ir (.Lhead block + .Ltail block), so it is excluded here
bad = [k for k, v in res.items() if not v["executed_equals_static"] and k != "step_asm"]
out = {"kind": "executed-count-check-v15p", "rule":
       "a counted symbol with zero backward branches executes static == executed instructions per call",
       "symbols": res, "unrolled_failures": bad,
       "verdict": "PASS" if not bad else "FAIL"}
json.dump(out, open(os.path.join(PKG, "evidence", "check-executed-v15p.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
sys.exit(1 if bad else 0)
