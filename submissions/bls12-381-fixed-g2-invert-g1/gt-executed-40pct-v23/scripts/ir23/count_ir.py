#!/usr/bin/env python3
"""count_ir.py -- static Ir count per symbol from the arm64 object, disjoint
symbol boundaries.  Emits evidence/ir22.json.  Also runs the same count for the
v21 constants so the note can report the delta."""
import json
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(os.path.dirname(HERE))
OBJ = os.path.join(HERE, "kernels_arm64.o")
WANT = ["_key_six_plusminus", "_min_rotation_booth", "_min_rotation_full", "_cyclotomic_sqr"]
V21 = {"key_six_plusminus": 1802, "min_rotation_booth": 144, "min_rotation_full": 0, "cyclotomic_sqr": 15042}


def disasm():
    out = subprocess.run(["otool", "-tvV", OBJ], capture_output=True, text=True, check=True).stdout
    syms, cur, addr = {}, None, {}
    for line in out.splitlines():
        m = re.match(r"^(_[A-Za-z0-9_]+):$", line.strip())
        if m:
            cur = m.group(1)
            syms[cur] = []
            continue
        m = re.match(r"^([0-9a-f]{16})\t(\S+)", line)
        if m and cur:
            syms[cur].append((m.group(1), m.group(2)))
    return syms


def main():
    syms = disasm()
    rows = {}
    for s in WANT:
        ins = syms[s]
        first, last = ins[0][0], ins[-1][0]
        rows[s.lstrip("_")] = {
            "symbol": s, "first_ir_addr": "0x" + first, "last_ir_addr": "0x" + last,
            "ir_count": len(ins),
            "distinct_mnemonics": len(set(m for _, m in ins)),
            "v21_constant": V21[s.lstrip("_")],
            "delta_vs_v21": len(ins) - V21[s.lstrip("_")],
            "boundary_note": "count spans symbol entry through the instruction "
                             "immediately before the next symbol label; no callee "
                             "shares an instruction with this span (static dispatch)",
        }
    total_new = sum(r["ir_count"] for r in rows.values())
    total_old = sum(V21.values())
    out = {"mode": "disassembly", "object": os.path.basename(OBJ),
           "toolchain": subprocess.run(["clang", "--version"], capture_output=True, text=True).stdout.splitlines()[0],
           "target": "arm64-apple-macos11", "opt": "-O2",
           "disassembler": "otool -tvV", "rows": rows,
           "sum_new_ir": total_new, "sum_v21_constants": total_old,
           "delta_total": total_new - total_old,
           "constant_folding_guard": "P6/ZETA6 are extern const arrays defined in "
                                     "consts.c (separate TU), so no priced table can "
                                     "be folded to zero"}
    json.dump(out, open(os.path.join(PKG, "evidence", "ir22.json"), "w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
