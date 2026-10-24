#!/usr/bin/env python3
"""row_extract.py -- count the shipped kernel from its DISASSEMBLY, and prove it.

Two numbers matter and v17a shipped only the first (verdict 40b720e2 point 3):

  span_rows    the STATIC instruction span _cios4 .. _cios4_end
  exec_rows    the EXECUTED instruction count of ONE call: every row outside the
               outer loop runs once, every row inside the loop's body runs 4
               times (i = 0..3).  The outer loop is found as the backward branch
               whose region contains all other rows; the body is [target, branch].

It also builds and runs the SELFTEST (a,b < 2^192 < n, n odd 256-bit) and checks
every (n,a,b,out) tuple against a Python Montgomery reference
out = a*b*R^{-1} mod n, R = 2^256 -- so the counted code is the multiplication it
claims to be.  Then it re-imports the ledger and asserts ledger.M_insn == exec_rows.

    python3 scripts/row_extract.py
"""
import json
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.dirname(HERE)
KER = os.path.join(HERE, "kernel")
SRC = os.path.join(KER, "cios4.c")
# the object file is a BINARY: it is built in a temp dir and never enters the
# package (intake rejects non-UTF8 files outright).
import tempfile
OBJ = os.path.join(tempfile.mkdtemp(prefix="cios4-v17-"), "cios4.o")
OUT = os.path.join(D, "evidence", "cios-disasm.json")
LINE = re.compile(r"^([0-9a-f]{16})\t(\S+)\t?(.*)$")
BRANCH = re.compile(r"^b(\.|$)")
START, END = "_cios4:", "_cios4_end:"
R = 1 << 256


def rows_of(txt):
    on, rows = False, []
    for ln in txt.splitlines():
        if ln.strip() == START:
            on = True
            continue
        if ln.strip() == END:
            on = False
            continue
        if on:
            m = LINE.match(ln)
            if m:
                rows.append({"addr": int(m.group(1), 16), "insn": m.group(2),
                             "operands": m.group(3).strip()})
    return rows


def executed(rows):
    """exec_rows = once-rows + 4 * body-rows, body = outermost backward branch."""
    loops = []
    for i, r in enumerate(rows):
        if BRANCH.match(r["insn"]):
            tgt = re.match(r"0x([0-9a-f]+)", r["operands"])
            if not tgt:
                continue
            t = int(tgt.group(1), 16)
            if t < r["addr"]:
                loops.append((t, r["addr"], i))
    if not loops:
        return len(rows), None, rows
    t, a, idx = min(loops, key=lambda x: x[0])          # widest backward branch
    body = [r for r in rows if t <= r["addr"] <= a]
    once = [r for r in rows if not (t <= r["addr"] <= a)]
    return len(once) + 4 * len(body), {"loop_start": hex(t), "loop_branch": hex(a),
                                       "body_rows": len(body),
                                       "once_rows": len(once),
                                       "iterations": 4}, body


def selftest():
    exe = "/tmp/cios4_test_v17"
    subprocess.run(["clang", "-O2", "-DSELFTEST", SRC, "-o", exe], check=True)
    out = subprocess.run([exe], check=True, capture_output=True, text=True).stdout
    checked = ok = 0
    for ln in out.splitlines():
        p = ln.split()
        if len(p) != 4:
            continue
        n, a, b, o = (int(x, 16) for x in p)
        checked += 1
        ok += ((a * b * pow(R, -1, n)) % n == o)
    return {"cases": checked, "matched": ok, "domain": "a,b < 2^192 < n, n odd 256-bit",
            "reference": "a*b*R^{-1} mod n with R = 2^256"}


def main():
    subprocess.run(["clang", "-O2", "-c", SRC, "-o", OBJ], check=True)
    txt = subprocess.run(["otool", "-tv", OBJ], check=True,
                         capture_output=True, text=True).stdout
    rows = rows_of(txt)
    assert rows, "empty span: label names changed"
    ex, loop, body = executed(rows)
    st = selftest()
    assert st["matched"] == st["cases"] and st["cases"] > 0, st
    out = {"test": "cios4_disasm_span",
           "kernel": "scripts/kernel/cios4.c",
           "compile_cmd": "clang -O2 -c scripts/kernel/cios4.c -o cios4.o",
           "disasm_cmd": "otool -tv cios4.o (span %s .. %s)" % (START, END),
           "host": subprocess.run(["uname", "-m"], capture_output=True,
                                  text=True).stdout.strip(),
           "span_rows": len(rows),
           "exec_rows": ex,
           "loop": loop,
           "selftest": st,
           "asserted_equal_ledger_row": False,
           "listing": ["%s  %s %s" % (hex(r["addr"])[2:].zfill(16), r["insn"],
                                      r["operands"]) for r in rows]}
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(out, f, indent=1)
    src = "import sys; sys.path.insert(0, %r); import ledger; print(ledger.M_insn)" % D
    led = subprocess.run([sys.executable, "-c", src], capture_output=True, text=True)
    if led.returncode != 0:
        print(led.stderr[-600:], file=sys.stderr)
        raise SystemExit("ledger import failed")
    led_insn = int(led.stdout.strip().splitlines()[-1])
    out["ledger_M_insn"] = led_insn
    out["asserted_equal_ledger_row"] = (led_insn == ex)
    with open(OUT, "w") as f:
        json.dump(out, f, indent=1)
    print("span %d, exec %d (%s), selftest %d/%d, ledger %d, equal %s"
          % (len(rows), ex, loop, st["matched"], st["cases"], led_insn,
             out["asserted_equal_ledger_row"]))
    assert out["asserted_equal_ledger_row"], "disassembly and ledger disagree"


if __name__ == "__main__":
    main()
