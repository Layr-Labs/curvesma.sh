#!/usr/bin/env python3
"""check_layout.py -- one layout: C struct vs the assembler constants.

Compiles scripts/layout_assert.c (whose _Static_asserts already fail the build on
any drift) and compares the compiler's sizeof/offsetof against the `.equ` values
the listing evidence/step_secp_v15g.s uses for its stride and field offsets.
Exits non-zero on any disagreement.
"""
import json, os, re, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__)); PKG = os.path.dirname(HERE)
exe = "/Users/andrewgordon/.hermes/cache/scratch/v15p-build/layout_assert"
os.makedirs(os.path.dirname(exe), exist_ok=True)
subprocess.run(["clang", "-O2", "-arch", "arm64", "-I", HERE,
                os.path.join(HERE, "layout_assert.c"), "-o", exe], check=True)
c = json.loads(subprocess.run([exe], capture_output=True, text=True, check=True).stdout)
asm = {}
for ln in open(os.path.join(PKG, "evidence", "step_secp_v15g.s")):
    m = re.match(r"\.equ\s+(REC_SIZE|OFF_IDX|OFF_LEN|OFF_ESC)\s*,\s*(\d+)", ln.strip())
    if m:
        asm[m.group(1)] = int(m.group(2))
want = {"REC_SIZE": c["sizeof_rec_t"], "OFF_IDX": c["off_idx"],
        "OFF_LEN": c["off_len"], "OFF_ESC": c["off_esc"]}
bad = {k: (asm.get(k), v) for k, v in want.items() if asm.get(k) != v}
out = {"kind": "layout-check-v15s", "compiler": c, "listing_equ": asm,
       "rule": ("sizeof(rec_t) and every field offset equal the listing's .equ constants; "
         "sizeof(worker_t) and its offsets are machine-derived with _Static_assert on every "
         "field, and the batch prefix array is WORKER_BATCH x 64 B"),
       "mismatches": bad, "verdict": "PASS" if not bad else "FAIL"}
json.dump(out, open(os.path.join(PKG, "evidence", "layout-check-v15s.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
sys.exit(1 if bad else 0)
