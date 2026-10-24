#!/usr/bin/env python3
"""
Row-enumeration completeness proof for the shipped span (v32, verdict-v31 flaw 2).

(v33: see proof.md section 9 for the change log this file implements.)
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import schedule  # noqa: E402
from span_asm import parse_disasm, norm_tokens  # noqa: E402

C_INSTR = 2000
N_BATCH_MIN = 1 << 18
FERMAT = 232864


def rows_in_emit_order(S):
    """The enumeration: every costed literal row, in the order emit_asm writes it."""
    prologue = [x for x in S if x.literal and x.cost and x.row == 15]
    prefix = [x for x in S if x.literal and x.cost and x.row == 13]
    body = [x for x in S if x.literal and x.cost and x.row not in (13, 15)]
    return prologue + prefix + body


def line_for(x):
    """The assembler line emit_asm writes for row x (LABELS substitution included)."""
    a = x.asm
    if a in schedule.LABELS:
        op, lab = schedule.LABELS[a]
        if op == "cbz":
            return "cbz %s, %s" % (a.split(" ", 1)[1].split(",")[0], lab)
        return "%s %s" % (op, lab)
    return a.replace(",", ", ")


def emitted_lines(src):
    out = []
    for l in src.splitlines():
        s = l.strip()
        if not s or s.startswith(".") or s.startswith("//") or s.endswith(":"):
            continue
        out.append(s)
    return out


def split_synthetic(lines, want):
    """Walk the emitted lines against the row-derived lines; the gaps are the synthetic
    lines.  Row lines must appear in order (a missing or reordered row is an error).
    Returns the synthetic lines and the aligned listing (line, row index or None)."""
    synth, aligned, errs, k = [], [], [], 0
    for ln in lines:
        if k < len(want) and ln == want[k]:
            aligned.append((ln, k))
            k += 1
        else:
            synth.append(ln)
            aligned.append((ln, None))
    if k != len(want):
        errs.append("enumeration: %d/%d row lines found in the emitted listing" % (k, len(want)))
    return synth, aligned, errs


def build_object():
    src = schedule.emit_asm(schedule.build())
    asm = os.path.join(HERE, "span.s")
    obj = os.path.join(HERE, "span.o")
    with open(asm, "w") as f:
        f.write(src + "\n")
    subprocess.run(["clang", "-arch", "arm64", "-c", asm, "-o", obj], check=True)
    dis = subprocess.run(["otool", "-t", "-v", obj], check=True,
                         capture_output=True, text=True).stdout
    ins = [d for d in (parse_disasm(l) for l in dis.splitlines()) if d]
    hdr = subprocess.run(["otool", "-l", obj], check=True,
                         capture_output=True, text=True).stdout
    size = None
    for blk in hdr.split("sectname __text"):
        for l in blk.splitlines():
            ls = l.strip()
            if ls.startswith("size 0x"):
                size = int(ls.split()[1], 16)
                break
        if size is not None:
            break
    return src, asm, obj, ins, size


def identity(ins, lines):
    """1:1 opcode identity, in order (operand-normalised for non-branch instructions)."""
    bad = []
    for d, ln in zip(ins, lines):
        parts = ln.split(None, 1)
        op = parts[0]
        tok = parts[1] if len(parts) > 1 else ""
        if op.startswith("b.") or op == "cbz" or op == "b":
            if d["op"] != op:
                bad.append((ln, d["op"], d["rest"]))
            continue
        if d["op"] != op or norm_tokens(tok.replace(" ", "")) != norm_tokens(d["rest"].replace(" ", "")):
            bad.append((ln, "%s %s" % (d["op"], d["rest"]),
                        norm_tokens(tok.replace(" ", ""))))
    return bad


def check(rows, synth, ins, size, charge, aligned, quiet=False):
    errs = []
    want = [line_for(x) for x in rows]
    total = len(rows) + len(synth)
    if len(ins) != total:
        errs.append("count: object holds %d instructions, the enumeration names %d + %d = %d"
                    % (len(ins), len(rows), len(synth), total))
    if size is None:
        errs.append("span-end: no __text size in the object header")
    elif size != 4 * len(ins):
        errs.append("span-end: __text is 0x%x bytes, %d instructions x 4 = 0x%x (trailing data)"
                    % (size, len(ins), 4 * len(ins)))
    if ins and size is not None and ins[-1]["addr"] + 4 != size:
        errs.append("span-end: last instruction at 0x%x does not end __text (0x%x)"
                    % (ins[-1]["addr"], size))
    unattributed = max(0, len(ins) - total)
    if unattributed:
        errs.append("attribution: %d object instruction(s) have no row and no synthetic line"
                    % unattributed)
    bad = identity(ins, [ln for ln, _ in aligned])
    if bad:
        errs.append("identity: %d line(s) differ from the object (first: %r vs %r)"
                    % (len(bad), bad[0][0], bad[0][1]))
    # charge reconciliation
    literal = len(rows)
    refresh = sum(1 for x in rows if x.refresh)
    nonrefresh_literal = literal - refresh
    nonrefresh_charge = charge["nonrefresh"]
    refresh_charge = charge["total"] - nonrefresh_charge
    if nonrefresh_literal + charge["amortised"] != nonrefresh_charge:
        errs.append("charge: %d non-refresh literal + %d amortised != %d charged"
                    % (nonrefresh_literal, charge["amortised"], nonrefresh_charge))
    if literal + charge["amortised"] != charge["total"]:
        errs.append("charge: %d literal + %d amortised != %d charged"
                    % (literal, charge["amortised"], charge["total"]))
    if refresh_charge != refresh:
        errs.append("charge: the refresh step adds %d, the listing has %d refresh rows"
                    % (refresh_charge, refresh))
    # scaffolding budget against the reserve inside the charged price
    per_element = len(synth) / float(N_BATCH_MIN)
    reserve = 1.0 - FERMAT / float(N_BATCH_MIN)
    if not per_element < reserve:
        errs.append("budget: scaffolding %.3g/element is not below the %.4f reserve" % (per_element, reserve))
    # the level price is the exact 64-cycle mean of the two listings, not a per-step max
    if 63 * nonrefresh_charge + charge["total"] != 64 * 131.046875:
        errs.append("price: (63 x %d + %d)/64 is not 131.046875" % (nonrefresh_charge, charge["total"]))
    if not quiet:
        print("  enumeration   : %d rows (%d non-refresh + %d refresh)" % (literal, nonrefresh_literal, refresh))
        print("  synthetic     : %d lines (2 returns + 5 stubs x 3)" % len(synth))
        print("  object        : %d instructions, __text 0x%x bytes, last ends at 0x%x"
              % (len(ins), size or 0, (ins[-1]["addr"] + 4) if ins else 0))
        print("  charge        : %d non-refresh = %d literal + %d amortised; %d refresh = %d literal + %d amortised"
              % (nonrefresh_charge, nonrefresh_literal, charge["amortised"], charge["total"],
                 literal, charge["amortised"]))
        print("  scaffolding   : %.3g instructions/element, reserve %.4f/element"
              % (per_element, reserve))
    return errs


def main():
    S = schedule.build()
    rows = rows_in_emit_order(S)
    src, asm, obj, ins, size = build_object()
    lines = emitted_lines(src)
    synth, aligned, errs0 = split_synthetic(lines, [line_for(x) for x in rows])
    charge = dict(total=sum(x.cost for x in S),
                  nonrefresh=sum(x.cost for x in S if not x.refresh),
                  amortised=sum(1 for x in S if x.cost and not x.literal))
    print("completeness -- row enumeration vs the shipped object (v32)\n")
    errs = errs0 + check(rows, synth, ins, size, charge, aligned)

    # ---- controls: each must fire
    ctl = []
    r1 = rows[:-1]
    ctl.append(("C1 delete one row from the enumeration",
                bool(check(r1, synth, ins, size, charge, aligned[1:], quiet=True))))
    r2 = [x for x in rows if not x.refresh]
    ctl.append(("C2 delete the three refresh rows",
                bool(check(r2, synth, ins, size, charge, aligned, quiet=True))))
    ctl.append(("C3 shave one line off the synthetic set",
                bool(check(rows, synth[:-1], ins, size, charge, aligned, quiet=True))))
    for name, hit in ctl:
        print("  %-38s -> %s" % (name, "DETECTED as required" if hit else "MISSED = CONTROL DEAD"))
    if not all(h for _, h in ctl):
        errs.append("control: a planted enumeration defect was not detected")

    # ---- table
    out = []
    A = out.append
    A("row-enumeration completeness -- plain BLS12-381 G1 span (v32)")
    A("=" * 96)
    A("source    : schedule.py build() rows in emit_asm() order, not the disassembly")
    A("object    : clang -arch arm64 -c span.s -o span.o ; otool -t -v / -l")
    A("")
    A("| # | row | refresh | asm (source) | object instruction | encoding | addr |")
    A("|---:|---:|---|---|---|---|---:|")
    want = [line_for(x) for x in rows]
    idx = 0
    for x, ln in zip(rows, want):
        d = ins[idx] if idx < len(ins) else dict(op="?", rest="", enc="", addr=0)
        A("| %d | %d | %s | %s | %s %s | %s | 0x%x |"
          % (idx + 1, x.row, "yes" if x.refresh else "no", ln, d["op"], d["rest"].strip(),
             d["enc"], d["addr"]))
        idx += 1
    A("")
    A("synthetic lines appended by emit_asm() -- belong to no ledger row:")
    A("")
    A("| # | asm | object instruction | encoding | addr |")
    A("|---:|---|---|---|---:|")
    for j, ln in enumerate(synth):
        d = ins[len(rows) + j]
        A("| %d | %s | %s %s | %s | 0x%x |" % (j + 1, ln, d["op"], d["rest"].strip(),
                                               d["enc"], d["addr"]))
    A("")
    A("accounting closure")
    A("-" * 96)
    A("  enumerated rows                 : %d  (%d literal non-refresh + %d refresh)"
      % (len(rows), len(rows) - sum(1 for x in rows if x.refresh), sum(1 for x in rows if x.refresh)))
    A("  synthetic scaffolding           : %d  (2 returns + 5 stubs x 3)" % len(synth))
    A("  object instruction count        : %d" % len(ins))
    A("  __text bytes                    : 0x%x = 4 x %d" % (size or 0, len(ins)))
    A("  unattributed instructions       : %d" % max(0, len(ins) - len(rows) - len(synth)))
    A("  charged per non-refresh step    : %d = %d literal + %d amortised (meter, Fermat share)"
      % (charge["nonrefresh"], len(rows) - 3, charge["amortised"]))
    A("  charged on the refresh step     : %d = %d literal + %d amortised"
      % (charge["total"], len(rows), charge["amortised"]))
    A("  scaffolding budget              : %d per batch / 2^18 elements = %.3e per element"
      % (len(synth), len(synth) / float(N_BATCH_MIN)))
    A("  reserve inside the charged price: 1 - 232,864/2^18 = %.4f per element"
      % (1.0 - FERMAT / float(N_BATCH_MIN)))
    A("  the level price 1.0655234375    : (63 x 131 + 134)/64 = 131.046875 exactly, no per-step max")
    A("")
    A("RESULT: %s" % ("COMPLETE -- every instruction in __text is enumerated or charged"
                      if not errs else "FAILED"))
    for e in errs:
        A("  FAIL " + e)
    with open(os.path.join(HERE, "completeness.txt"), "w") as f:
        f.write("\n".join(out) + "\n")

    for e in errs:
        print("  FAIL " + e)
    if not errs:
        print("\nCOMPLETE: %d object instructions = %d rows + %d synthetic; __text ends at 0x%x"
              % (len(ins), len(rows), len(synth), size))
        print("evidence: completeness.txt")
    if "--keep" not in sys.argv:
        for p in (asm, obj):
            try:
                os.unlink(p)
            except OSError:
                pass
    return 1 if errs else 0


if __name__ == "__main__":
    sys.exit(main())
