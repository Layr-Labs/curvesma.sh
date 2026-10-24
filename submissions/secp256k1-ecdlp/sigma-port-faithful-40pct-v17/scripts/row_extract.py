#!/usr/bin/env python3
"""row_extract.py -- count the shipped kernel's instructions from its DISASSEMBLY.

Intake objection (d717195a point 4): "the asserted kernel disassembly measurements
are not supplied".  This script supplies them.  It compiles scripts/kernel/cios4.c,
disassembles the object with otool, extracts the instruction span between the
`_cios4` and `_cios4_end` labels, counts the rows, writes the listing and the count
to evidence/cios-disasm.json, and then verifies the count against the value the
ledger actually prices (ledger.M_insn, which reads that same file) by importing the
ledger in a subprocess.  If the two disagree the script exits non-zero -- the claim
and the artifact cannot drift apart.

    python3 scripts/row_extract.py

Artifact: evidence/cios-disasm.json  {span_rows, listing[], compile_cmd, ...}
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
OBJ = os.path.join(KER, "cios4.o")
OUT = os.path.join(D, "evidence", "cios-disasm.json")
LINE = re.compile(r"^([0-9a-f]{16})\t(\S+)\t?(.*)$")
START, END = "_cios4:", "_cios4_end:"


def build_and_disassemble():
    subprocess.run(["clang", "-O2", "-c", SRC, "-o", OBJ], check=True)
    txt = subprocess.run(["otool", "-tv", OBJ], check=True,
                         capture_output=True, text=True).stdout
    return txt


def extract(txt):
    on, rows = False, []
    for ln in txt.splitlines():
        if ln.strip() == START:
            on = True
            continue
        if ln.strip() == END:
            on = False
            continue
        if not on:
            continue
        m = LINE.match(ln)
        if m:
            rows.append({"addr": m.group(1), "insn": m.group(2),
                         "operands": m.group(3).strip()})
    return rows


def main():
    txt = build_and_disassemble()
    rows = extract(txt)
    assert rows, "empty span: label names changed"
    out = {"test": "cios4_disasm_span",
           "kernel": "scripts/kernel/cios4.c",
           "compile_cmd": "clang -O2 -c scripts/kernel/cios4.c -o cios4.o",
           "disasm_cmd": "otool -tv cios4.o (span %s .. %s)" % (START, END),
           "host": subprocess.run(["uname", "-m"], capture_output=True,
                                  text=True).stdout.strip(),
           "span_rows": len(rows),
           "asserted_equal_ledger_row": False,
           "listing": ["%s  %s %s" % (r["addr"], r["insn"], r["operands"])
                       for r in rows]}
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(out, f, indent=1)
    # second pass: does the LEDGER price this count?  (it reads this file)
    src = ("import sys; sys.path.insert(0, %r); import ledger; print(ledger.M_insn)"
           % D)
    led = subprocess.run([sys.executable, "-c", src], capture_output=True,
                         text=True)
    if led.returncode != 0:
        print(led.stderr[-800:], file=sys.stderr)
        raise SystemExit("ledger import failed; cannot verify the row")
    led_insn = int(led.stdout.strip().splitlines()[-1])
    out["ledger_M_insn"] = led_insn
    out["asserted_equal_ledger_row"] = (led_insn == len(rows))
    with open(OUT, "w") as f:
        json.dump(out, f, indent=1)
    print("span_rows %d, ledger M_insn %d, equal %s, host %s"
          % (len(rows), led_insn, out["asserted_equal_ledger_row"], out["host"]))
    assert out["asserted_equal_ledger_row"], "disassembly and ledger disagree"


if __name__ == "__main__":
    main()
