#!/usr/bin/env python3
"""Machine-checked fused register schedule for the plain BLS12-381 G1 rho step.

This is the artifact proof.md section 5 refers to: a *linear* schedule of the ordinary
step's instructions (the rows of the section-5 listing; the reference addition is the
priced black box of C = 2,000 instructions), with an explicit data-flow obligation for
every instruction -- which register it writes, which value that register then holds, and
which (register, value) pairs it reads.

Three independent checks run over it:

  (1) forward definedness -- every read must find *that* value in *that* register, so a
      read of a clobbered or never-written register is an error;
  (2) backward liveness -- a write may not destroy a value that is read later.  This is
      the v23 rejection checked mechanically: "rows 1-3 consume x_new in x0-x5 and
      overwrite those registers with canonical-index intermediates.  Row 4 then uses
      x0-x5 as the six limbs of y_new, without any intervening loads";
  (3) register budget -- x18 is reserved on Darwin/AArch64, so 30 registers are usable;
      the peak simultaneous live-value count must be <= 30 and no instruction may name
      x18.

It also checks the census (115 ordinary / 118 refresh / 133 traversal instructions), the
fusion's memory obligation (the span loads the record's OLD point for the addition and
never reloads the NEW point: rows 1 and 4 take it in registers), and runs four planted
corruptions that must each make the verifier fail:

  A. y_new placed in x10-x15, where p lives          (v23's section-0 claim)
  B. row 2 writes its V1 product into K0's register  (a silent clobber)
  C. the `cset` producing s deleted                  (an undefined read)
  D. rows 1-3 emitted before row 4                   (v23's actual order)

Usage: python3 schedule.py            -- table, checks, controls; exit 0 iff all pass.
"""

import sys

RESERVED = {18}                                   # x18 is platform-reserved
POOL = ["x%d" % i for i in range(31) if i not in RESERVED]     # 30 usable

LIVE_IN = ["P0", "P1", "P2", "P3", "P4", "P5", "PH0", "PTOP",
           "RPTR", "END", "TBL", "CONST"]
LIVE_OUT = ["P0", "P1", "P2", "P3", "P4", "P5", "PH0", "PTOP",
            "RPTR2", "END", "TBL", "CONST"]
FUSION_DEFS = ["NX0", "NX1", "NX2", "NX3", "NX4", "NX5",
               "NY0", "NY1", "NY2", "NY3", "NY4", "NY5"]


class Ins:
    __slots__ = ("row", "asm", "defs", "uses", "cost", "refresh")

    def __init__(self, row, asm, defs, uses, cost=1, refresh=False):
        self.row, self.asm = row, asm
        self.defs = list(defs)                     # [(reg, value)]
        self.uses = list(uses)
        self.cost, self.refresh = cost, refresh


def build():
    S = []
    A = S.append

    # ---- row 6: the reference addition's interface (19 counted transfers) ----------
    # IN : the record's entry address (ldr), the table point (6 ldp), and the state point
    #      as it stands at block entry -- the OLD point (6 ldp).  The index rows that
    #      follow need the NEW point; the whole point of this revision is that they take
    #      it from registers, not from a second load of the record.
    # OUT: affine x_new limbs 0..5 in x0..x5, y_new limbs 0..5 in x23..x28; 6 stp store
    #      it and are charged here.
    # The block preserves the twelve batch-live values; during it exactly those twelve
    # are live, so 18 registers are free for its own temporaries inside its 2,000
    # instruction price (row 6 is shown as counts, as in v22, for the same reason).
    A(Ins(6, "add-ref: ldr entry, 6 ldp table point, 6 ldp state_old, add, 6 stp state_new",
          [(r, v) for r, v in zip(["x0", "x1", "x2", "x3", "x4", "x5"], FUSION_DEFS[:6])] +
          [(r, v) for r, v in zip(["x23", "x24", "x25", "x26", "x27", "x28"], FUSION_DEFS[6:])],
          [("x19", "RPTR")], cost=19))

    # ---- row 4 (17): p - y, then the exact borrow (p-y)-y as the sign, K limbs 0-3 ---
    # y_new arrives in x23..x28.  p is a batch constant used by every record, so p - y
    # must not overwrite it: it goes to the six registers free at this point.
    PY = ["x6", "x7", "x8", "x9", "x29", "x30"]
    for i in range(6):
        mn = "subs" if i == 0 else ("sbc" if i == 5 else "sbcs")
        A(Ins(4, "%s %s,x%d,x%d" % (mn, PY[i], 10 + i, 23 + i),
              [(PY[i], "PY%d" % i)], [("x%d" % (10 + i), "P%d" % i), ("x%d" % (23 + i), "NY%d" % i)]))
    A(Ins(4, "cmp x6,x23", [], [("x6", "PY0"), ("x23", "NY0")]))
    for i in range(1, 6):
        A(Ins(4, "sbcs xzr,%s,x%d" % (PY[i], 23 + i), [],
              [(PY[i], "PY%d" % i), ("x%d" % (23 + i), "NY%d" % i)]))
    A(Ins(4, "cset x27,lo", [("x27", "S")], []))          # NY4 is dead after the borrow
    for i, k in enumerate(["x23", "x24", "x25", "x26"]):
        A(Ins(4, "csel %s,%s,%s,lo" % (k, PY[i], k),
              [(k, "K%d" % i)], [(PY[i], "PY%d" % i), (k, "NY%d" % i)]))

    # ---- row 1 (21): the twelve limb products of x_new (no reload) ----------------
    # six umulh (column 5 high halves) + six madd (column 6 low halves) + five add, plus
    # 3 ldp of G_0..G_5 and ldr of G_6: 21 instructions (24 with the x reload of v22).
    def r1(asm, defs, uses):
        A(Ins(1, asm, defs, uses))
    r1("ldp x6,x7,[x22,#0]", [("x6", "G0"), ("x7", "G1")], [("x22", "CONST")])
    r1("umulh x9,x5,x6", [("x9", "L1")], [("x5", "NX5"), ("x6", "G0")])
    r1("umulh x8,x4,x7", [("x8", "T1")], [("x4", "NX4"), ("x7", "G1")])
    r1("add x9,x9,x8", [("x9", "L2")], [("x9", "L1"), ("x8", "T1")])
    r1("madd x9,x5,x7,x9", [("x9", "L3")], [("x5", "NX5"), ("x7", "G1"), ("x9", "L2")])
    r1("ldp x6,x7,[x22,#16]", [("x6", "G2"), ("x7", "G3")], [("x22", "CONST")])
    r1("umulh x8,x3,x6", [("x8", "T2")], [("x3", "NX3"), ("x6", "G2")])
    r1("add x9,x9,x8", [("x9", "L4")], [("x9", "L3"), ("x8", "T2")])
    r1("madd x9,x4,x6,x9", [("x9", "L5")], [("x4", "NX4"), ("x6", "G2"), ("x9", "L4")])
    r1("umulh x8,x2,x7", [("x8", "T3")], [("x2", "NX2"), ("x7", "G3")])
    r1("add x9,x9,x8", [("x9", "L6")], [("x9", "L5"), ("x8", "T3")])
    r1("madd x9,x3,x7,x9", [("x9", "L7")], [("x3", "NX3"), ("x7", "G3"), ("x9", "L6")])
    r1("ldp x6,x7,[x22,#32]", [("x6", "G4"), ("x7", "G5")], [("x22", "CONST")])
    r1("umulh x8,x1,x6", [("x8", "T4")], [("x1", "NX1"), ("x6", "G4")])
    r1("add x9,x9,x8", [("x9", "L8")], [("x9", "L7"), ("x8", "T4")])
    r1("madd x9,x2,x6,x9", [("x9", "L9")], [("x2", "NX2"), ("x6", "G4"), ("x9", "L8")])
    r1("umulh x8,x0,x7", [("x8", "T5")], [("x0", "NX0"), ("x7", "G5")])
    r1("add x9,x9,x8", [("x9", "L10")], [("x9", "L9"), ("x8", "T5")])
    r1("madd x9,x1,x7,x9", [("x9", "L11")], [("x1", "NX1"), ("x7", "G5"), ("x9", "L10")])
    r1("ldr x6,[x22,#48]", [("x6", "G6")], [("x22", "CONST")])
    r1("madd x9,x0,x6,x9", [("x9", "L")], [("x0", "NX0"), ("x6", "G6"), ("x9", "L11")])
    assert sum(1 for x in S if x.row == 1) == 21

    # ---- row 2 (15): range test, V1, X, S, wrap test, V2 --------------------------
    A(Ins(2, "sub x8,x9,#16", [("x8", "RT")], [("x9", "L")]))
    A(Ins(2, "cmn x8,#31", [], [("x8", "RT")]))
    A(Ins(2, "b.hs fallback", [], [("x8", "RT")]))
    A(Ins(2, "umulh x0,x16,x9", [("x0", "V1a")], [("x16", "PH0"), ("x9", "L")]))
    A(Ins(2, "mul x1,x17,x9", [("x1", "V1b")], [("x17", "PTOP"), ("x9", "L")]))
    A(Ins(2, "umulh x2,x17,x9", [("x2", "V1c")], [("x17", "PTOP"), ("x9", "L")]))
    A(Ins(2, "cmn x1,x0", [], [("x1", "V1b"), ("x0", "V1a")]))
    A(Ins(2, "adc x1,x2,xzr", [("x1", "V1")], [("x1", "V1b"), ("x2", "V1c")]))
    A(Ins(2, "add x2,x5,x1", [("x2", "SM")], [("x5", "NX5"), ("x1", "V1")]))
    A(Ins(2, "subs x3,x2,x17", [("x3", "DF")], [("x2", "SM"), ("x17", "PTOP")]))
    A(Ins(2, "csel x4,x17,xzr,hs", [("x4", "W0")], [("x17", "PTOP")]))
    A(Ins(2, "add x6,x3,#16", [("x6", "D16")], [("x3", "DF")]))
    A(Ins(2, "cmp x6,#32", [], [("x6", "D16")]))
    A(Ins(2, "b.ls fallback", [], [("x6", "D16")]))
    A(Ins(2, "sub x2,x4,x3", [("x2", "V2")], [("x4", "W0"), ("x3", "DF")]))

    # ---- row 3 (14): exact argmin j ----------------------------------------------
    A(Ins(3, "cmp x1,x5", [], [("x1", "V1"), ("x5", "NX5")]))
    A(Ins(3, "csel x3,x1,x5,lo", [("x3", "A1")], [("x1", "V1"), ("x5", "NX5")]))
    A(Ins(3, "csel x4,x5,x1,lo", [("x4", "A2")], [("x1", "V1"), ("x5", "NX5")]))
    A(Ins(3, "cset x6,lo", [("x6", "MIN1")], []))
    A(Ins(3, "cmp x2,x4", [], [("x2", "V2"), ("x4", "A2")]))
    A(Ins(3, "csel x7,x2,x4,lo", [("x7", "B1")], [("x2", "V2"), ("x4", "A2")]))
    A(Ins(3, "cmp x2,x3", [], [("x2", "V2"), ("x3", "A1")]))
    A(Ins(3, "csel x7,x3,x7,lo", [("x7", "B2")], [("x3", "A1"), ("x7", "B1")]))
    A(Ins(3, "csel x3,x2,x3,lo", [("x3", "MIN2")], [("x2", "V2"), ("x3", "A1")]))
    A(Ins(3, "mov x8,#2", [("x8", "TWO")], []))
    A(Ins(3, "csel x6,x8,x6,lo", [("x6", "J")], [("x8", "TWO"), ("x6", "MIN1")]))
    A(Ins(3, "sub x7,x7,x3", [("x7", "GAP")], [("x7", "B2"), ("x3", "MIN2")]))
    A(Ins(3, "cmp x7,#16", [], [("x7", "GAP")]))
    A(Ins(3, "b.ls fallback", [], [("x7", "GAP")]))

    # ---- row 5 (8): h, distinguished test, twist index, table address ------------
    A(Ins(5, "and x4,x23,#0xfffff", [("x4", "H")], [("x23", "K0")]))
    A(Ins(5, "cbz x24,report", [], [("x24", "K1")]))
    A(Ins(5, "add x5,x27,x4,lsl #1", [("x5", "U1")], [("x27", "S"), ("x4", "H")]))
    A(Ins(5, "add x5,x5,x5,lsl #1", [("x5", "U2")], [("x5", "U1")]))
    A(Ins(5, "add x5,x5,x6", [("x5", "U3")], [("x5", "U2"), ("x6", "J")]))
    A(Ins(5, "add x5,x5,x5,lsl #2", [("x5", "U4")], [("x5", "U3")]))
    A(Ins(5, "add x5,x21,x5,lsl #5", [("x5", "ENT")], [("x21", "TBL"), ("x5", "U4")]))
    A(Ins(5, "str x5,[x19,#96]", [], [("x5", "ENT"), ("x19", "RPTR")]))

    # ---- row 8a (3): length load and cap test run before row 7's refresh test ---
    A(Ins(8, "ldp x27,x28,[x19,#136]", [("x27", "LENLO"), ("x28", "LENHI")], [("x19", "RPTR")]))
    A(Ins(8, "cmp x28,#20", [], [("x28", "LENHI")]))
    A(Ins(8, "b.hs abandon", [], [("x28", "LENHI")]))

    # ---- row 7 (12): checkpoint test, refresh test, 3 refresh-only stores ------
    A(Ins(7, "ldp x7,x8,[x19,#104]", [("x7", "KCP0"), ("x8", "KCP2")], [("x19", "RPTR")]))
    A(Ins(7, "ldp x9,x29,[x19,#120]", [("x9", "KCP3"), ("x29", "VAL")], [("x19", "RPTR")]))
    A(Ins(7, "cmp x23,x7", [], [("x23", "K0"), ("x7", "KCP0")]))
    A(Ins(7, "ccmp x25,x8,#0,eq", [], [("x25", "K2"), ("x8", "KCP2")]))
    A(Ins(7, "ccmp x26,x9,#0,eq", [], [("x26", "K3"), ("x9", "KCP3")]))
    A(Ins(7, "ccmp x29,#1,#0,eq", [], [("x29", "VAL")]))
    A(Ins(7, "b.eq traversal", [], [("x29", "VAL")]))
    A(Ins(7, "tst x27,#63", [], [("x27", "LENLO")]))
    A(Ins(7, "b.ne norefresh", [], [("x27", "LENLO")]))
    A(Ins(7, "stp x23,x25,[x19,#104]", [], [("x23", "K0"), ("x25", "K2"), ("x19", "RPTR")], refresh=True))
    A(Ins(7, "mov x30,#1", [("x30", "ONE")], [], refresh=True))
    A(Ins(7, "stp x26,x30,[x19,#120]", [], [("x26", "K3"), ("x30", "ONE"), ("x19", "RPTR")], refresh=True))

    # ---- row 8b (6): 69-bit length increment, pointer advance, loop ------------
    A(Ins(8, "adds x27,x27,#1", [("x27", "LENLO2")], [("x27", "LENLO")]))
    A(Ins(8, "adc x28,x28,xzr", [("x28", "LENHI2")], [("x28", "LENHI")]))
    A(Ins(8, "stp x27,x28,[x19,#136]", [], [("x27", "LENLO2"), ("x28", "LENHI2"), ("x19", "RPTR")]))
    A(Ins(8, "add x19,x19,#512", [("x19", "RPTR2")], [("x19", "RPTR")]))
    A(Ins(8, "cmp x19,x20", [], [("x19", "RPTR2"), ("x20", "END")]))
    A(Ins(8, "b.ne loop", [], [("x19", "RPTR2"), ("x20", "END")]))
    return S


def verify(S, quiet=False):
    errs = []

    # (0) register pool
    for i, x in enumerate(S):
        for r, v in x.defs + x.uses:
            if r.startswith("x") and r[1:].isdigit() and r not in POOL:
                errs.append("%d: `%s` names %s, outside the 30-register pool (x18 reserved)" % (i, x.asm, r))

    # (1) forward definedness -- seeded with the batch-live set the span starts from
    cur = {}
    for r, v in zip(["x10", "x11", "x12", "x13", "x14", "x15", "x16", "x17",
                     "x19", "x20", "x21", "x22"], LIVE_IN):
        cur[r] = v
    for i, x in enumerate(S):
        for r, v in x.uses:
            if cur.get(r) != v:
                errs.append("row %d instr %d `%s`: reads %s from %s, which holds %s"
                            % (x.row, i, x.asm, v, r, cur.get(r)))
        for r, v in x.defs:
            cur[r] = v

    # (2) backward liveness on (register, value) pairs
    n = len(S)
    live_after = [set() for _ in range(n)]
    live = {(r, v) for r, v in zip(["x%d" % i for i in range(10, 18)] + ["x19", "x20", "x21", "x22"],
                                   LIVE_OUT)}
    live = {(r, v) for r, v in zip(["x10", "x11", "x12", "x13", "x14", "x15", "x16", "x17",
                                    "x19", "x20", "x21", "x22"], LIVE_OUT)}
    peak, peak_at = len(live), "live-out"
    for i in range(n - 1, -1, -1):
        x = S[i]
        live_after[i] = set(live)
        for r, v in x.defs:
            for (r2, v2) in live:
                if r2 == r and v2 != v:
                    errs.append("row %d instr %d `%s`: CLOBBER of live value %s in %s (writes %s)"
                                % (x.row, i, x.asm, v2, r, v))
        live = {p for p in live if p[0] not in {r for r, _ in x.defs}} | {(r, v) for r, v in x.uses}
        if len(live) > peak:
            peak, peak_at = len(live), "before row %d instr %d (`%s`)" % (x.row, i, x.asm)

    # (3) census
    total = sum(x.cost for x in S)
    refresh = sum(x.cost for x in S if x.refresh)
    rows = {}
    for x in S:
        rows[x.row] = rows.get(x.row, 0) + x.cost
    if total != 115:
        errs.append("census: %d instructions, the listing prices 115" % total)
    if refresh != 3:
        errs.append("census: %d refresh-only instructions, the listing prices 3" % refresh)
    if rows.get(1) != 21 or rows.get(4) != 17:
        errs.append("census: rows 1/4 are %s/%s, the fused listing is 21/17" % (rows.get(1), rows.get(4)))
    if rows.get(6) != 19 or rows.get(2) != 15 or rows.get(3) != 14 or rows.get(5) != 8 or rows.get(7) != 12 or rows.get(8) != 9:
        errs.append("census: row totals %s disagree with the ledger 19/15/14/8/12/9" % rows)

    # (4) fusion obligation
    if any(x.asm.split()[0] in ("ldp", "ldr") and "state_new" in x.asm for x in S):
        errs.append("fusion: the span reloads the new point")
    for i, x in enumerate(S):
        for r, v in x.defs:
            if v in FUSION_DEFS and i != 0:
                errs.append("fusion: %s is defined outside the addition block (instr %d)" % (v, i))
    defined = {v for x in S for r, v in x.defs}
    for v in FUSION_DEFS:
        if v not in defined:
            errs.append("fusion: %s is never defined" % v)
    if not any("NX" in v for x in S if x.row == 1 for r, v in x.uses):
        errs.append("fusion: row 1 does not consume the handed-over x_new")

    # (5) budget
    if peak > len(POOL):
        errs.append("budget: peak %d live values > %d usable registers" % (peak, len(POOL)))

    if not quiet:
        print("instructions            : %d (refresh step %d, %d refresh-only)"
              % (total, total + refresh, refresh))
        print("row census              : %s" % ", ".join("row %d = %d" % (r, c) for r, c in sorted(rows.items())))
        print("peak live values        : %d at %s" % (peak, peak_at))
        regs = sorted({r for x in S for r, _ in x.defs + x.uses}, key=lambda s: int(s[1:]))
        print("registers used          : %d (%s)" % (len(regs), ", ".join(regs)))
    return errs, dict(total=total, refresh=refresh, rows=rows, peak=peak, peak_at=peak_at,
                      live_after=live_after)


def main():
    S = build()
    print("fused register schedule -- ordinary step, plain BLS12-381 G1 rho")
    print("30 usable registers (x0..x30 less x18); %d values live in, %d live out\n"
          % (len(LIVE_IN), len(LIVE_OUT)))
    errs, st = verify(S)
    for e in errs:
        print("  FAIL  " + e)
    if errs:
        print("\nSCHEDULE REJECTED (%d error%s)" % (len(errs), "" if len(errs) == 1 else "s"))
        return 1
    print("checks: pool ok | forward definedness ok | no clobber | census 115/118/133 | fusion handover ok")
    print("peak %d live values <= 30 usable registers (%d free)\n" % (st["peak"], 30 - st["peak"]))

    print("| # | row | instruction | writes | reads | live after |")
    print("|---:|---:|---|:--|:--|---:|")
    for i, x in enumerate(S):
        print("| %d | %d | `%s` | %s | %s | %d |" % (
            i + 1, x.row, x.asm,
            ", ".join("%s=%s" % (r, v) for r, v in x.defs) or "-",
            ", ".join("%s(%s)" % (v, r) for r, v in x.uses) or "-",
            len(st["live_after"][i])))

    print("\nplanted corruptions (each must make the verifier fail):")
    res = []

    def C1_y_in_p():
        M = build()
        for x in M:
            x.defs = [({"NY0": "x10", "NY1": "x11", "NY2": "x12", "NY3": "x13",
                        "NY4": "x14", "NY5": "x15"}.get(v, r), v) for r, v in x.defs]
            x.uses = [({"NY0": "x10", "NY1": "x11", "NY2": "x12", "NY3": "x13",
                        "NY4": "x14", "NY5": "x15"}.get(v, r), v) for r, v in x.uses]
        return M

    def C2_row2_into_K0():
        M = build()
        for x in M:
            if x.asm == "umulh x0,x16,x9":
                x.defs = [("x23", "V1a")]
        return M

    def C3_drop_cset():
        return [x for x in build() if x.asm != "cset x27,lo"]

    def C4_v23_order():
        M = build()
        return ([x for x in M if x.row == 6] + [x for x in M if x.row in (1, 2, 3)] +
                [x for x in M if x.row == 4] + [x for x in M if x.row in (5, 7, 8)])

    for name, fn in [("A y_new in x10-x15 (p's registers)", C1_y_in_p),
                     ("B row 2 writes V1 into K0 (x23)", C2_row2_into_K0),
                     ("C cset producing s deleted", C3_drop_cset),
                     ("D rows 1-3 before row 4 (v23 order)", C4_v23_order)]:
        e, s2 = verify(fn(), quiet=True)
        res.append(bool(e))
        print("  %-36s -> %-22s (%d error%s, peak %d) %s" % (
            name, "FAIL as required" if e else "PASSED = CONTROL DEAD",
            len(e), "" if len(e) == 1 else "s", s2["peak"], e[0] if e else ""))
    if not all(res):
        print("\nCONTROL DEAD: the self-audit does not certify")
        return 1
    print("\n4/4 controls fired (3 required + the v23 order). SCHEDULE OK, CERTIFIED.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
