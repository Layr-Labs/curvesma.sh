#!/usr/bin/env python3
"""count_ir.py -- static Ir count per priced symbol from the arm64 object, with
disjoint symbol boundaries.  This is the ONLY source of the v25 kernel prices.

v25: the object is built from the CORRECTED kernels.c (f6mul cross terms repaired,
f6mulv made alias-safe, f12mul/f6inv/f2inv/f12inv added).  cyclotomic_sqr's span is
the squaring price; f12mul's span x 1200 + compare is the verification price;
key_six_plusminus and min_rotation_booth are the per-step key/rotation rows.
"""
import json
import os
import re
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(os.path.dirname(HERE))
OBJ = os.path.join(HERE, "kernels_arm64.o")
WANT = ["_key_six_plusminus", "_min_rotation_booth", "_min_rotation_full",
        "_cyclotomic_sqr", "_f12mul"]


def disasm(obj):
    out = subprocess.run(["otool", "-tvV", obj], capture_output=True, text=True,
                         check=True).stdout
    syms, cur = {}, None
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
    subprocess.run(["clang", "-O2", "-target", "arm64-apple-macos11", "-c",
                    os.path.join(HERE, "kernels.c"), "-o", OBJ], check=True)
    syms = disasm(OBJ)
    rows = {}
    for s in WANT:
        ins = syms[s]
        rows[s.lstrip("_")] = {
            "symbol": s, "first_ir_addr": "0x" + ins[0][0],
            "last_ir_addr": "0x" + ins[-1][0], "ir_count": len(ins),
            "distinct_mnemonics": len(set(m for _, m in ins)),
        }
    # lookup symbol: separate object
    lk_obj = os.path.join(HERE, "lookup_arm64.o")
    subprocess.run(["clang", "-O2", "-target", "arm64-apple-macos11", "-c",
                    os.path.join(HERE, "lookup_count.c"), "-o", lk_obj], check=True)
    lk = disasm(lk_obj)
    lsym = [s for s in lk if s.startswith("_lookup")]
    assert lsym, lk.keys()
    rows["lookup_probe"] = {
        "symbol": lsym[0], "first_ir_addr": "0x" + lk[lsym[0]][0][0],
        "last_ir_addr": "0x" + lk[lsym[0]][-1][0],
        "ir_count": len(lk[lsym[0]]),
        "distinct_mnemonics": len(set(m for _, m in lk[lsym[0]])),
    }
    out = {"mode": "disassembly", "object": os.path.basename(OBJ),
           "toolchain": subprocess.run(["clang", "--version"], capture_output=True,
                                        text=True).stdout.splitlines()[0],
           "target": "arm64-apple-macos11", "opt": "-O2",
           "disassembler": "otool -tvV", "rows": rows,
           "boundary_note": "count spans the symbol entry through the instruction "
                            "immediately before the next symbol label; the kernels "
                            "are static-dispatch, no callee shares an instruction",
           "constant_folding_guard": "P6/ZETA6/ONE6 are extern const arrays defined "
                                     "in a separate TU (consts.c)"}
    json.dump(out, open(os.path.join(PKG, "evidence", "ir25.json"), "w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
