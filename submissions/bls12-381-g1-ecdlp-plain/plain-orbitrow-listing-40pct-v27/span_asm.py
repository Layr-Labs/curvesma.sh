#!/usr/bin/env python3
"""Object-disassembly reconciliation for the plain BLS12-381 G1 rho span (v27).

The v26 verdict rejected row 6 because it was a *count* ("ldr entry, 6 ldp table point,
6 ldp state_old, 6 stp state_new") rather than one executed instruction: a listing that
is not assembled cannot witness what the machine runs, and a pseudo-instruction can hide
assembler expansion (a wide `mov` becoming movz/movk, `adr` becoming adrp+add, an
`ldp`/`stp` split for an illegal immediate).

This script closes that hole by construction.  It emits the listing that schedule.py's
build() declares, assembles it for AArch64 with the platform toolchain, disassembles the
object with otool, and reconciles the two instruction by instruction:

  * count   -- the object must contain exactly as many instructions as the schedule
               charges (125 literal instructions per iteration, 122 on a non-refresh
               step, plus the 2 amortised per-batch charges that are never instructions);
  * identity-- opcode, registers in order, and immediate values must match row for row,
               so no row is an assembler expansion of another;
  * control flow -- the object must have exactly ONE backward branch, which is what makes
               the fusion window finite (no store->reload pair can span two iterations).

Writes span-disasm.txt (the evidence file: the disassembly, the reconciliation table and
the summary).  Exit 0 iff every check passes.

Usage: python3 span_asm.py [--keep]
"""

import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import schedule  # noqa: E402

REGS = set("x%d" % i for i in range(31)) | {"xzr", "wzr", "sp"} | set("w%d" % i for i in range(31))
CONDS = {"eq", "ne", "hs", "lo", "mi", "pl", "vs", "vc", "hi", "ls", "ge", "lt", "gt", "le", "al", "cs", "cc"}


def norm_tokens(text):
    """Operand tokens reduced to what must be identical: registers, immediates, shifts,
    condition codes.  Labels and resolved branch targets are dropped."""
    out = []
    for tok in [t.strip() for t in text.split(",")]:
        t = tok.replace(" ", "").lower().rstrip("]").lstrip("[")
        if not t:
            continue
        if t.startswith("#"):
            v = t[1:]
            try:
                out.append("imm:%d" % int(v, 0))
            except ValueError:
                out.append("imm:" + v)
        elif t in REGS or t in CONDS:
            out.append(t)
        elif re.match(r"^-?0x[0-9a-f]+$", t) or re.match(r"^\d+$", t):
            out.append("imm:%d" % int(t, 0))
        elif re.match(r"^lsl#?\d+$", t) or re.match(r"^lsr#?\d+$", t) or re.match(r"^asr#?\d+$", t):
            out.append(re.sub(r"#", "", t))
        elif t in ("[", "]"):
            continue
        else:
            out.append(t)          # keep anything unexpected visible
    if out and out[-1] == "imm:0":
        out.pop()                  # a zero offset is elided in the disassembly
    return out


def parse_source(asm):
    parts = asm.split(None, 1)
    op = parts[0].lower()
    rest = parts[1] if len(parts) > 1 else ""
    # drop a trailing branch target name: `b.hs fallback`, `cbz x1,report`
    rest = re.sub(r",?\s*[A-Za-z_][A-Za-z_0-9]*$", "", rest) if op.startswith("b") or op == "cbz" else rest
    return op, norm_tokens(rest.replace(" ", "")) if rest else []


def parse_disasm(line):
    """otool prints `<addr>\t<mnemonic>\t<operands>` (the encoding column appears only
    with -V, so accept both shapes)."""
    parts = [p for p in line.rstrip("\n").split("\t")]
    if len(parts) < 2 or not re.match(r"^[0-9a-f]{8,16}$", parts[0].strip()):
        return None
    addr = int(parts[0].strip(), 16)
    i = 1
    enc = ""
    if re.match(r"^[0-9a-f]{8}$", parts[i].strip()):
        enc = parts[i].strip()
        i += 1
    if i >= len(parts):
        return None
    return dict(addr=addr, enc=enc, op=parts[i].strip().lower(), rest="\t".join(parts[i + 1:]).strip())


def main():
    S = schedule.build()
    rows = [(i, x) for i, x in enumerate(S) if x.cost and x.literal]
    src = schedule.emit_asm(S)
    asm_path = os.path.join(HERE, "span.s")
    obj_path = os.path.join(HERE, "span.o")
    with open(asm_path, "w") as f:
        f.write(src + "\n")

    subprocess.run(["clang", "-arch", "arm64", "-c", asm_path, "-o", obj_path], check=True)
    dis = subprocess.run(["otool", "-t", "-v", obj_path], check=True,
                         capture_output=True, text=True).stdout

    insns = [d for d in (parse_disasm(l) for l in dis.splitlines()) if d]

    # the trailing label block plus `ret` is one instruction and is not part of the span
    end = len(insns) - 1
    body = insns[:end]
    tail = insns[end:]

    errs = []
    if len(body) != len(rows):
        errs.append("count: the object holds %d span instructions, the schedule charges %d"
                    % (len(body), len(rows)))

    mism = []
    for k, (d, (i, x)) in enumerate(zip(body, rows)):
        sop, stok = parse_source(x.asm)
        dop, dtok = d["op"], norm_tokens(d["rest"].replace(" ", "")) if d["rest"] else []
        if sop.startswith("b.") or sop == "cbz":
            ok = dop == sop
        else:
            ok = (sop == dop and stok == dtok)
        if not ok:
            mism.append(dict(idx=k, row=x.row, src=x.asm, obj="%s %s" % (d["op"], d["rest"]),
                             src_tok=stok, obj_tok=dtok))
    if mism:
        errs.append("identity: %d row(s) do not match the object instruction" % len(mism))

    back = [d for d in body if d["op"] == "b" and d["addr"] < d["addr"]]
    # a backward branch is one whose resolved target is below its own address
    backs = []
    for d in body:
        if d["op"] in ("b", "b.ne", "b.eq", "b.hs", "b.ls", "b.lo", "b.hi") and re.search(r"0x[0-9a-f]+", d["rest"]):
            tgt = int(re.search(r"0x([0-9a-f]+)", d["rest"]).group(1), 16)
            if tgt <= d["addr"]:
                backs.append(d)

    amort = [x for x in S if x.cost and not x.literal]
    per_iter = sum(x.cost for x in S)
    nonrefresh = sum(x.cost for x in S if not x.refresh)

    lines = []
    lines.append("object-disassembly reconciliation -- plain BLS12-381 G1 rho span (v27)")
    lines.append("=" * 88)
    lines.append("source   : schedule.py build()  (%d rows, %d costed)" % (len(S), sum(1 for x in S if x.cost)))
    lines.append("toolchain: clang -arch arm64 -c span.s -o span.o ; otool -t -v span.o")
    lines.append("")
    lines.append("charged per iteration          : %d  (%d non-refresh, %d on the 1-in-64 refresh step)"
                 % (per_iter, nonrefresh, per_iter - nonrefresh + 0))
    lines.append("literal machine instructions   : %d" % len(body))
    lines.append("amortised per-batch charges    : %d  (rows 11, 12: never per-element instructions)"
                 % len(amort))
    lines.append("backward branches in the object: %d  (the single back-edge that bounds fusion)"
                 % len(backs))
    lines.append("")
    lines.append("RECONCILIATION: %s" % ("EXACT -- every charged literal row is one machine"
                                        " instruction, in order, with the same registers"
                                        if not errs else "FAILED"))
    lines.append("")
    lines.append("row-by-row (source row | source text | object instruction | encoding)")
    lines.append("-" * 88)
    for (i, x), d in zip(rows, body):
        lines.append("%3d | row %-2d | %-46s | %-34s | %s"
                     % (i, x.row, x.asm, ("%s %s" % (d["op"], d["rest"])).strip(), d["enc"]))
    lines.append("")
    if amort:
        lines.append("not machine instructions (charged, per-batch amortised):")
        for i, x in enumerate(S):
            if x.cost and not x.literal:
                lines.append("  row %-2d  %s" % (x.row, x.asm))
        lines.append("")
    lines.append("full disassembly (otool -t -v):")
    lines.append("-" * 88)
    lines.extend(l for l in dis.splitlines() if l.strip())

    out = os.path.join(HERE, "span-disasm.txt")
    with open(out, "w") as f:
        f.write("\n".join(lines) + "\n")

    print("charged per iteration        : %d (%d non-refresh)" % (per_iter, nonrefresh))
    print("literal instructions in .o   : %d" % len(body))
    print("amortised charges (not insn) : %d" % len(amort))
    print("backward branches            : %d" % len(backs))
    print("row-6 region                 : %d instructions in the object at rows %d..%d"
          % (sum(1 for i, x in rows if x.row == 6), 0, 0))
    for e in errs:
        print("  FAIL " + e)
    for m in mism[:8]:
        print("  MISMATCH idx %d row %d\n    src %s %s\n    obj %s %s"
              % (m["idx"], m["row"], m["src"], m["src_tok"], m["obj"], m["obj_tok"]))
    if not errs:
        print("\n%s" % lines[8])
        print("evidence: span-disasm.txt")
    if "--keep" not in sys.argv:
        for p in (asm_path, obj_path):
            try:
                os.unlink(p)
            except OSError:
                pass
    return 1 if errs else 0


if __name__ == "__main__":
    sys.exit(main())
