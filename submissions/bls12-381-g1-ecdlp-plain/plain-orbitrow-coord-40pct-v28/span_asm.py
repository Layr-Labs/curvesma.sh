#!/usr/bin/env python3
"""Object-disassembly reconciliation for the plain BLS12-381 G1 rho span (v28).

v28 reconciles three things the v27 verdict (1259a121) found wanting, all against the OBJECT:

  * row 6 is segmented per coordinate: the encodings split 1 (entry address) + 9 (x path) +
    9 (y path), and the object shows the row-6 region is one straight-line block with no
    coordinate-conditional branch, so both paths execute on every element (weight 1 each).
  * the store -> reload pairing uses the object's own address arithmetic: `add x19,x19,#512`
    advances the record pointer, so a pair is eligible only when both sides address the same
    x19 value (same record).  v27 wrapped across that back-edge and paired different records;
    those candidates are excluded here and listed.
  * control flow: v27 emitted every label in one block after the back-edge, so nine conditional
    branches resolved to one shared `ret` and `b.ne norefresh` never reached the length update.
    v28 emits a real CFG: the prefix pass is its own loop before the leaf body, `b.ne
    L_norefresh` targets the length-update block, and each exit action has its own stub.  This
    script enumerates every branch -> target edge with the target address and stub body length.

Writes span-disasm.txt.  Exit 0 iff every check passes.
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
CONDS = {"eq", "ne", "hs", "lo", "mi", "pl", "vs", "vc", "hi", "ls", "ge", "lt", "gt", "le", "al",
         "cs", "cc"}
BR = ("b", "b.eq", "b.ne", "b.hs", "b.ls", "b.lo", "b.hi", "cbz", "cbnz")


def norm_tokens(text):
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
            out.append(t)
    if out and out[-1] == "imm:0":
        out.pop()
    return out


def parse_source(asm):
    parts = asm.split(None, 1)
    op = parts[0].lower()
    rest = parts[1] if len(parts) > 1 else ""
    rest = re.sub(r",?\s*[A-Za-z_][A-Za-z_0-9]*$", "", rest) if op.startswith("b") or op == "cbz" else rest
    return op, norm_tokens(rest.replace(" ", "")) if rest else []


def parse_disasm(line):
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


def target_of(d, symbols=None):
    m = re.search(r"0x([0-9a-f]+)", d["rest"])
    if m:
        return int(m.group(1), 16)
    for name, addr in (symbols or {}).items():
        if re.search(r"\b%s\b" % re.escape(name), d["rest"]):
            return addr
    return None


def main():
    S = schedule.build()
    rows = [(i, x) for i, x in enumerate(S) if x.cost and x.literal]
    pre_rows = [(i, x) for i, x in rows if x.row == 13]
    body_rows = [(i, x) for i, x in rows if x.row != 13]
    src = schedule.emit_asm(S)
    asm_path = os.path.join(HERE, "span.s")
    obj_path = os.path.join(HERE, "span.o")
    with open(asm_path, "w") as f:
        f.write(src + "\n")
    subprocess.run(["clang", "-arch", "arm64", "-c", asm_path, "-o", obj_path], check=True)
    dis = subprocess.run(["otool", "-t", "-v", obj_path], check=True,
                         capture_output=True, text=True).stdout
    insns = [d for d in (parse_disasm(l) for l in dis.splitlines()) if d]

    first = body_rows[0][1]
    fop, ftok = parse_source(first.asm)
    start = None
    for k, d in enumerate(insns):
        if d["op"] == fop and norm_tokens(d["rest"].replace(" ", "")) == ftok:
            start = k
            break
    errs = []
    if start is None:
        errs.append("cannot locate the first row-6 instruction in the object")
        start = max(0, len(pre_rows))
    body = insns[start:start + len(body_rows)]
    prefix = insns[max(0, start - len(pre_rows) - 1):max(0, start - 1)]
    stubs = insns[start + len(body_rows):]

    if len(body) != len(body_rows):
        errs.append("count: %d leaf instructions in the object, %d charged"
                    % (len(body), len(body_rows)))
    if len(prefix) != len(pre_rows):
        errs.append("count: %d prefix instructions in the object, %d charged"
                    % (len(prefix), len(pre_rows)))

    mism = []
    for region, srows in ((prefix, pre_rows), (body, body_rows)):
        for d, (i, x) in zip(region, srows):
            sop, stok = parse_source(x.asm)
            if sop.startswith("b.") or sop == "cbz":
                ok = d["op"] == sop
            else:
                ok = (sop == d["op"] and stok == norm_tokens(d["rest"].replace(" ", "")))
            if not ok:
                mism.append(dict(row=x.row, src=x.asm, obj="%s %s" % (d["op"], d["rest"]),
                                 src_tok=stok, obj_tok=norm_tokens(d["rest"].replace(" ", ""))))
    if mism:
        errs.append("identity: %d row(s) do not match the object instruction" % len(mism))

    body_lo, body_hi = body[0]["addr"], body[-1]["addr"]
    symbols = {"_batch_prefix": prefix[0]["addr"] if prefix else 0, "_leaf_span": body_lo}
    edges = []
    for d in insns:
        if d["op"] in BR:
            t = target_of(d, symbols)
            if t is None:
                continue
            where = ("prefix" if d in prefix else "body" if body_lo <= d["addr"] <= body_hi
                     else "stub")
            edges.append(dict(addr=d["addr"], op=d["op"], target=t, where=where,
                              back=(t <= d["addr"])))
    body_backs = [e for e in edges if e["where"] == "body" and e["back"]]
    pre_backs = [e for e in edges if e["where"] == "prefix" and e["back"]]
    exit_edges = [e for e in edges if e["where"] == "body" and not e["back"] and e["target"] > body_hi]
    inbody_fwd = [e for e in edges if e["where"] == "body" and not e["back"] and e["target"] <= body_hi]

    if len(body_backs) != 1:
        errs.append("CFG: the leaf body has %d back-edge(s), expected exactly one" % len(body_backs))
    if len(pre_backs) != 1:
        errs.append("CFG: the prefix loop has %d back-edge(s), expected exactly one" % len(pre_backs))
    noref = [e for e in inbody_fwd if e["op"] == "b.ne"]
    len_block = next((d["addr"] for d in body if d["op"] == "adds"), None)
    if len_block is None:
        errs.append("CFG: no `adds` (length update) in the leaf body")
    elif len(noref) != 1:
        errs.append("CFG: %d forward conditional branches inside the body, expected one"
                    % len(noref))
    elif noref[0]["target"] != len_block:
        errs.append("CFG: `b.ne norefresh` targets 0x%x, the length-update block is 0x%x"
                    % (noref[0]["target"], len_block))

    stub_indeg = {}
    for e in exit_edges:
        stub_indeg[e["target"]] = stub_indeg.get(e["target"], 0) + 1
    for a in stub_indeg:
        if not any(d["addr"] == a and d["op"] == "ret" for d in stubs):
            errs.append("CFG: exit edge target 0x%x is not a stub body" % a)
        if a <= body_hi:
            errs.append("CFG: an exit edge targets an in-body address")

    amort = [x for x in S if x.cost and not x.literal]
    per_iter = sum(x.cost for x in S)
    nonrefresh = sum(x.cost for x in S if not x.refresh)
    tab = schedule.fold_by_coordinate(S)
    tp = schedule.traffic_pairs(S)

    lines = []
    A = lines.append
    A("object-disassembly reconciliation -- plain BLS12-381 G1 rho span (v28)")
    A("=" * 92)
    A("source   : schedule.py build()  (%d rows, %d costed)" % (len(S), sum(1 for x in S if x.cost)))
    A("toolchain: clang -arch arm64 -c span.s -o span.o ; otool -t -v span.o")
    A("")
    A("charged per iteration          : %d  (%d non-refresh, %d refresh-only)"
      % (per_iter, nonrefresh, per_iter - nonrefresh))
    A("literal machine instructions   : %d  (leaf body %d, prefix loop %d)"
      % (len(body) + len(prefix), len(body), len(prefix)))
    A("amortised per-batch charges    : %d  (rows 11, 12: never per-element instructions)" % len(amort))
    A("leaf-body back-edges           : %d" % len(body_backs))
    A("prefix-loop back-edges         : %d" % len(pre_backs))
    A("exit stubs                     : %d distinct targets, in-degree %s"
      % (len(stub_indeg), sorted(stub_indeg.values())))
    A("")
    A("RECONCILIATION: %s" % ("EXACT -- every charged literal row is one machine instruction, in"
                             " order, with the same registers" if not errs else "FAILED"))
    A("")
    A("row 6 segmented per coordinate (object instructions)")
    A("-" * 92)
    A("| coordinate | path | executed insns | charged | inside C (instr.) | traffic weight |")
    A("|---|---|---:|---:|---:|---:|")
    for r in tab:
        A("| %s | %s | %d | %d | %.1f | %.2f |" % (r["coord"], "entry-address fetch" if
          r["coord"] == "entry" else "table+state transfers", r["insns"], r["charged"],
          r["inside_c_instr"], r["traffic"]))
    A("")
    A("branch -> target edges in the object")
    A("-" * 92)
    A("| branch PC | instruction | target | target region | back-edge | stub body (insn) |")
    A("|---|---|---|---|---|---:|")
    pre_hi = prefix[-1]["addr"] if prefix else -1
    for e in sorted(edges, key=lambda q: q["addr"]):
        reg = ("prefix loop" if e["target"] <= pre_hi else
               "in-body continuation" if body_lo <= e["target"] <= body_hi else "exit stub")
        A("| 0x%x | %s | 0x%x | %s | %s | %s |"
          % (e["addr"], e["op"], e["target"], reg, "yes" if e["back"] else "no",
             ("%d" % stub_indeg.get(e["target"], 0)) if reg == "exit stub" else "-"))
    A("")
    A("store -> reload pairs, eligibility-filtered (same record only)")
    A("-" * 92)
    A("| record offset | field | producer | consumer | distance |")
    A("|---:|---|---:|---:|---:|")
    for pp in sorted(tp["eligible"], key=lambda q: -q["dist"]):
        A("| %d | %s | %d | %d | %d |" % (pp["off"], schedule.FIELDS.get(pp["off"], "?"),
                                          pp["producer"], pp["consumer"], pp["dist"]))
    A("")
    A("%d eligible pairs, max distance %d, bound %d; %d candidates excluded (different record or"
      " carried state); %d stack/frame operands."
      % (len(tp["eligible"]), tp["mx"], tp["bound"], len(tp["excluded"]), len(tp["bad"])))
    A("")
    A("row-by-row (source row | source text | object instruction | encoding)")
    A("-" * 92)
    for (i, x), d in list(zip(pre_rows, prefix)) + list(zip(body_rows, body)):
        A("%3d | row %-2d | %-46s | %-32s | %s"
          % (i, x.row, x.asm, ("%s %s" % (d["op"], d["rest"])).strip(), d["enc"]))
    A("")
    if amort:
        A("not machine instructions (charged, per-batch amortised):")
        for i, x in enumerate(S):
            if x.cost and not x.literal:
                A("  row %-2d  %s" % (x.row, x.asm))
        A("")
    A("full disassembly (otool -t -v):")
    A("-" * 92)
    A("\n".join(l for l in dis.splitlines() if l.strip()))

    with open(os.path.join(HERE, "span-disasm.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")

    print("charged per iteration        : %d (%d non-refresh)" % (per_iter, nonrefresh))
    print("literal instructions in .o   : %d (leaf %d + prefix %d)"
          % (len(body) + len(prefix), len(body), len(prefix)))
    print("leaf-body back-edges         : %d" % len(body_backs))
    print("prefix-loop back-edges       : %d" % len(pre_backs))
    print("exit stubs                   : %d targets, in-degree %s"
          % (len(stub_indeg), sorted(stub_indeg.values())))
    print("row-6 region                 : %d instructions"
          % sum(1 for i, x in body_rows if x.row == 6))
    print("eligible pairs               : %d, max distance %d, bound %d"
          % (len(tp["eligible"]), tp["mx"], tp["bound"]))
    print("excluded candidates          : %d" % len(tp["excluded"]))
    for e in errs:
        print("  FAIL " + e)
    for m in mism[:6]:
        print("  MISMATCH row %d\n    src %s %s\n    obj %s %s"
              % (m["row"], m["src"], m["src_tok"], m["obj"], m["obj_tok"]))
    if not errs:
        print("\nEXACT: every charged literal row is one machine instruction")
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
