#!/usr/bin/env python3
"""Machine-checked liveness/census contract for the plain BLS12-381 G1 rho step (v26 listing).

This is the artifact proof.md section 5.1 refers to: a *linear* schedule of the ordinary
step's instructions (rows 1-8 of the section-5 listing, plus the four one-per-element
batched-addition rows 10-13), with an explicit data-flow obligation for every instruction --
which register it writes, which value that register then holds, and which (register, value)
pairs it reads. The reference addition itself is the priced black box of C = 2,000
instructions (row 6 is its record interface: 19 counted transfers, shown as counts, as in
v22 and v25, because its 36 scalar transfers stream through the 18 free registers inside
the addition's own price).

Four independent checks run over it:

  (1) forward definedness -- every read must find *that* value in *that* register, so a
      read of a never-written register, of a clobbered register, or of the wrong register
      for a value is an error;
  (2) backward liveness -- a write may not destroy a value that is read later, which is the
      contract that killed v23 ("rows 1-3 consume x_new in x0-x5 and overwrite those
      registers ... row 4 then uses x0-x5 as the six limbs of y_new, without any
      intervening loads") and the contract the twelve batch-live values must survive;
  (3) register budget -- x18 is reserved on Darwin/AArch64, so 30 registers are usable; the
      peak simultaneous live-value count must be <= 30 and no instruction may name x18;
  (4) census -- 124 charged instructions on a non-refresh step, 3 of them refresh-only
      (127 on the 1-in-64 refresh step), with the row totals of the section-5.1 traffic
      table: 24/15/14/20/8/19/9+3/9/1/1/1/3.

It then runs the **three planted liveness-contract violations the v25 brief required and
that package did not ship** (0/3 fired), plus a planted census flip, each of which must
make the verifier fail:

  A. clobber-while-live     -- row 2 writes its V1 product into x23, where the new point's
                               y_0 limb (NY0) lives until row 4 reads it;
  B. missing-preservation   -- row 5 folds the table address into x21, destroying TBL, one
                               of the twelve batch-live values the span must hand back;
  C. wrong-output-register  -- row 1 takes x_2 from the truncated-product product instead
                               of from the record load, i.e. NX2 is written to x8 and the
                               product reads a register that never held it;
  D. census flip            -- one ledger row deleted, so the census no longer equals 124.

Usage: python3 schedule.py     # table, checks, controls; exit 0 iff all pass.
"""

import sys

RESERVED = {18}                                              # x18 is platform-reserved
POOL = ["x%d" % i for i in range(31) if i not in RESERVED]    # 30 usable

LIVE_IN = ["P0", "P1", "P2", "P3", "P4", "P5", "PH0", "PTOP", "RPTR", "END", "TBL", "CONST"]
LIVE_OUT = ["P0", "P1", "P2", "P3", "P4", "P5", "PH0", "PTOP", "RPTR2", "END", "TBL", "CONST"]
PINNED = ["x10", "x11", "x12", "x13", "x14", "x15", "x16", "x17", "x19", "x20", "x21", "x22"]


class Ins:
    __slots__ = ("row", "asm", "defs", "uses", "cost", "refresh")

    def __init__(self, row, asm, defs=(), uses=(), cost=1, refresh=False):
        self.row, self.asm = row, asm
        self.defs, self.uses = list(defs), list(uses)
        self.cost, self.refresh = cost, refresh


def build():
    S = []
    A = S.append
    NY = ["x23", "x24", "x25", "x26", "x27", "x28"]
    NX = ["x0", "x1", "x2", "x3", "x4", "x5"]
    PY = ["x6", "x7", "x9", "x29", "x30", "x5"]

    # ---- row 6 (19): the reference addition's record interface ---------------------
    # ldr entry address + 6 ldp table point + 6 ldp state (old) + 6 stp new state.
    # Shown as counts (its 36 scalar transfers stream through the addition's 18 free
    # registers); it reads the record pointer and must preserve the twelve pinned values.
    A(Ins(6, "add-ref: ldr entry, 6 ldp table point, 6 ldp state_old, 6 stp state_new",
          [("x0", "TX0")], [("x19", "RPTR")], cost=19))

    # ---- row 1 (24): the twelve limb products of x_new, reloaded from the record ----
    A(Ins(1, "ldp x0,x1,[x19,#0]", [("x0", "NX0"), ("x1", "NX1")], [("x19", "RPTR")]))
    A(Ins(1, "ldp x2,x3,[x19,#16]", [("x2", "NX2"), ("x3", "NX3")], [("x19", "RPTR")]))
    A(Ins(1, "ldp x4,x5,[x19,#32]", [("x4", "NX4"), ("x5", "NX5")], [("x19", "RPTR")]))
    A(Ins(1, "ldp x6,x7,[x22,#0]", [("x6", "G0"), ("x7", "G1")], [("x22", "CONST")]))
    A(Ins(1, "umulh x8,x5,x6", [("x8", "L1")], [("x5", "NX5"), ("x6", "G0")]))
    A(Ins(1, "umulh x9,x4,x7", [("x9", "T1")], [("x4", "NX4"), ("x7", "G1")]))
    A(Ins(1, "add x8,x8,x9", [("x8", "L2")], [("x8", "L1"), ("x9", "T1")]))
    A(Ins(1, "madd x8,x5,x7,x8", [("x8", "L3")], [("x5", "NX5"), ("x7", "G1"), ("x8", "L2")]))
    A(Ins(1, "ldp x6,x7,[x22,#16]", [("x6", "G2"), ("x7", "G3")], [("x22", "CONST")]))
    A(Ins(1, "umulh x9,x3,x6", [("x9", "T2")], [("x3", "NX3"), ("x6", "G2")]))
    A(Ins(1, "add x8,x8,x9", [("x8", "L4")], [("x8", "L3"), ("x9", "T2")]))
    A(Ins(1, "madd x8,x4,x6,x8", [("x8", "L5")], [("x4", "NX4"), ("x6", "G2"), ("x8", "L4")]))
    A(Ins(1, "umulh x9,x2,x7", [("x9", "T3")], [("x2", "NX2"), ("x7", "G3")]))
    A(Ins(1, "add x8,x8,x9", [("x8", "L6")], [("x8", "L5"), ("x9", "T3")]))
    A(Ins(1, "madd x8,x3,x7,x8", [("x8", "L7")], [("x3", "NX3"), ("x7", "G3"), ("x8", "L6")]))
    A(Ins(1, "ldp x6,x7,[x22,#32]", [("x6", "G4"), ("x7", "G5")], [("x22", "CONST")]))
    A(Ins(1, "umulh x9,x1,x6", [("x9", "T4")], [("x1", "NX1"), ("x6", "G4")]))
    A(Ins(1, "add x8,x8,x9", [("x8", "L8")], [("x8", "L7"), ("x9", "T4")]))
    A(Ins(1, "madd x8,x2,x6,x8", [("x8", "L9")], [("x2", "NX2"), ("x6", "G4"), ("x8", "L8")]))
    A(Ins(1, "umulh x9,x0,x7", [("x9", "T5")], [("x0", "NX0"), ("x7", "G5")]))
    A(Ins(1, "add x8,x8,x9", [("x8", "L10")], [("x8", "L9"), ("x9", "T5")]))
    A(Ins(1, "madd x8,x1,x7,x8", [("x8", "L11")], [("x1", "NX1"), ("x7", "G5"), ("x8", "L10")]))
    A(Ins(1, "ldr x6,[x22,#48]", [("x6", "G6")], [("x22", "CONST")]))
    A(Ins(1, "madd x9,x0,x6,x8", [("x9", "L")], [("x0", "NX0"), ("x6", "G6"), ("x8", "L11")]))

    # ---- row 2 (15): range test, V1, X, S, wrap test, V2 ---------------------------
    A(Ins(2, "sub x7,x9,#16", [("x7", "RT")], [("x9", "L")]))
    A(Ins(2, "cmn x7,#31", [], [("x7", "RT")]))
    A(Ins(2, "b.hs fallback", [], [("x7", "RT")]))
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

    # ---- row 3 (14): exact argmin j -----------------------------------------------
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
    A(Ins(3, "csel x8,x8,x6,lo", [("x8", "J")], [("x8", "TWO"), ("x6", "MIN1")]))
    A(Ins(3, "sub x7,x7,x3", [("x7", "GAP")], [("x7", "B2"), ("x3", "MIN2")]))
    A(Ins(3, "cmp x7,#16", [], [("x7", "GAP")]))
    A(Ins(3, "b.ls fallback", [], [("x7", "GAP")]))

    # ---- row 4 (20): y_new reload, p - y, borrow (p-y)-y as the sign, K limbs 0-3 ---
    A(Ins(4, "ldp x23,x24,[x19,#48]", [("x23", "NY0"), ("x24", "NY1")], [("x19", "RPTR")]))
    A(Ins(4, "ldp x25,x26,[x19,#64]", [("x25", "NY2"), ("x26", "NY3")], [("x19", "RPTR")]))
    A(Ins(4, "ldp x27,x28,[x19,#80]", [("x27", "NY4"), ("x28", "NY5")], [("x19", "RPTR")]))
    for i, (mn, dst) in enumerate(zip(["subs", "sbcs", "sbcs", "sbcs", "sbcs", "sbc"], PY)):
        A(Ins(4, "%s %s,x%d,x%d" % (mn, dst, 10 + i, 23 + i),
              [(dst, "PY%d" % i)], [("x%d" % (10 + i), "P%d" % i), ("x%d" % (23 + i), "NY%d" % i)]))
    A(Ins(4, "cmp x6,x23", [], [("x6", "PY0"), ("x23", "NY0")]))
    for i in range(1, 6):
        A(Ins(4, "sbcs xzr,%s,x%d" % (PY[i], 23 + i), [],
              [(PY[i], "PY%d" % i), ("x%d" % (23 + i), "NY%d" % i)]))
    A(Ins(4, "cset x4,lo", [("x4", "S")], []))
    for i, k in enumerate(["x0", "x1", "x2", "x3"]):
        A(Ins(4, "csel %s,%s,x%d,lo" % (k, PY[i], 23 + i),
              [(k, "K%d" % i)], [(PY[i], "PY%d" % i), ("x%d" % (23 + i), "NY%d" % i)]))

    # ---- row 5 (8): h, distinguished test, twist index, table address --------------
    A(Ins(5, "and x5,x0,#0xfffff", [("x5", "H")], [("x0", "K0")]))
    A(Ins(5, "cbz x1,report", [], [("x1", "K1")]))
    A(Ins(5, "add x6,x4,x5,lsl #1", [("x6", "U1")], [("x4", "S"), ("x5", "H")]))
    A(Ins(5, "add x6,x6,x6,lsl #1", [("x6", "U2")], [("x6", "U1")]))
    A(Ins(5, "add x6,x6,x8", [("x6", "U3")], [("x6", "U2"), ("x8", "J")]))
    A(Ins(5, "add x6,x6,x6,lsl #2", [("x6", "U4")], [("x6", "U3")]))
    A(Ins(5, "add x6,x21,x6,lsl #5", [("x6", "ENT")], [("x21", "TBL"), ("x6", "U4")]))
    A(Ins(5, "str x6,[x19,#96]", [], [("x6", "ENT"), ("x19", "RPTR")]))

    # ---- row 8a (3): length load and cap test run before row 7's refresh test ------
    A(Ins(8, "ldp x27,x28,[x19,#136]", [("x27", "LENLO"), ("x28", "LENHI")], [("x19", "RPTR")]))
    A(Ins(8, "cmp x28,#20", [], [("x28", "LENHI")]))
    A(Ins(8, "b.hs abandon", [], [("x28", "LENHI")]))

    # ---- row 7 (9 + 3 refresh): checkpoint test, refresh test, refresh stores ------
    A(Ins(7, "ldp x7,x8,[x19,#104]", [("x7", "KCP0"), ("x8", "KCP2")], [("x19", "RPTR")]))
    A(Ins(7, "ldp x9,x29,[x19,#120]", [("x9", "KCP3"), ("x29", "VAL")], [("x19", "RPTR")]))
    A(Ins(7, "cmp x0,x7", [], [("x0", "K0"), ("x7", "KCP0")]))
    A(Ins(7, "ccmp x2,x8,#0,eq", [], [("x2", "K2"), ("x8", "KCP2")]))
    A(Ins(7, "ccmp x3,x9,#0,eq", [], [("x3", "K3"), ("x9", "KCP3")]))
    A(Ins(7, "ccmp x29,#1,#0,eq", [], [("x29", "VAL")]))
    A(Ins(7, "b.eq traversal", [], [("x29", "VAL")]))
    A(Ins(7, "tst x27,#63", [], [("x27", "LENLO")]))
    A(Ins(7, "b.ne norefresh", [], [("x27", "LENLO")]))
    A(Ins(7, "stp x0,x2,[x19,#104]", [], [("x0", "K0"), ("x2", "K2"), ("x19", "RPTR")], refresh=True))
    A(Ins(7, "mov x26,#1", [("x26", "ONE")], [], refresh=True))
    A(Ins(7, "stp x3,x26,[x19,#120]", [], [("x3", "K3"), ("x26", "ONE"), ("x19", "RPTR")], refresh=True))

    # ---- row 8b (6): 69-bit length increment, pointer advance, loop ---------------
    A(Ins(8, "adds x27,x27,#1", [("x27", "LENLO2")], [("x27", "LENLO")]))
    A(Ins(8, "adc x28,x28,xzr", [("x28", "LENHI2")], [("x28", "LENHI")]))
    A(Ins(8, "stp x27,x28,[x19,#136]", [], [("x27", "LENLO2"), ("x28", "LENHI2"), ("x19", "RPTR")]))
    A(Ins(8, "add x19,x19,#512", [("x19", "RPTR2")], [("x19", "RPTR")]))
    A(Ins(8, "cmp x19,x20", [], [("x19", "RPTR2"), ("x20", "END")]))
    A(Ins(8, "b.ne loop", [], [("x19", "RPTR2"), ("x20", "END")]))

    # ---- rows 10-13 (6): the batched-addition overhead the section-5.1 table charges
    A(Ins(10, "cbz batch-product,zero-recover", [], [], cost=1))
    A(Ins(11, "add meter.lo,meter.lo,batch-size", [], [], cost=1))
    A(Ins(12, "share of the one Fermat inversion (232,864/N, N >= 2^18)", [], [], cost=1))
    A(Ins(13, "fwd prefix pass: ldr entry, subs counter, b.ne", [], [("x19", "RPTR2")], cost=3))
    return S


def verify(S, quiet=False):
    errs = []

    # (0) register pool
    for i, x in enumerate(S):
        for r, v in x.defs + x.uses:
            if r.startswith("x") and r[1:].isdigit() and r not in POOL:
                errs.append("%d: `%s` names %s, outside the 30-register pool (x18 reserved)" % (i, x.asm, r))

    # (1) forward definedness
    cur = dict(zip(PINNED, LIVE_IN))
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
    live = set(zip(PINNED, LIVE_OUT))
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

    # (3) register budget
    if peak > len(POOL):
        errs.append("budget: peak %d live values > %d usable registers" % (peak, len(POOL)))

    # (4) census, against the section-5.1 traffic table
    total = sum(x.cost for x in S)
    refresh = sum(x.cost for x in S if x.refresh)
    rows = {}
    for x in S:
        rows[x.row] = rows.get(x.row, 0) + x.cost
    want = {1: 24, 2: 15, 3: 14, 4: 20, 5: 8, 6: 19, 7: 12, 8: 9, 10: 1, 11: 1, 12: 1, 13: 3}
    if total - refresh != 124:
        errs.append("census: %d charged instructions, the section-5.1 table prices 124" % (total - refresh))
    if refresh != 3:
        errs.append("census: %d refresh-only instructions, the listing prices 3" % refresh)
    if rows != want:
        errs.append("census: row totals %s disagree with the traffic table %s" % (rows, want))
    if total != 127:
        errs.append("census: refresh step %d instructions, the listing prices 127" % total)

    if not quiet:
        print("charged instructions    : %d non-refresh (%d on the 1-in-64 refresh step, %d refresh-only)"
              % (total - refresh, total, refresh))
        print("row census              : %s" % ", ".join("row %d = %d" % (r, c) for r, c in sorted(rows.items())))
        print("peak live values        : %d at %s" % (peak, peak_at))
        regs = sorted({r for x in S for r, _ in x.defs + x.uses}, key=lambda s: int(s[1:]))
        print("registers touched       : %d (%s)" % (len(regs), ", ".join(regs)))
    return errs, dict(total=total, refresh=refresh, rows=rows, peak=peak, peak_at=peak_at,
                      live_after=live_after)


def main():
    S = build()
    print("register/liveness contract -- ordinary step, plain BLS12-381 G1 rho (v26 listing)")
    print("30 usable registers (x0..x30 less x18); %d values live in, %d live out\n"
          % (len(LIVE_IN), len(LIVE_OUT)))
    errs, st = verify(S)
    for e in errs:
        print("  FAIL  " + e)
    if errs:
        print("\nCONTRACT REJECTED (%d error%s)" % (len(errs), "" if len(errs) == 1 else "s"))
        return 1
    print("checks: pool ok | forward definedness ok | no clobber | census 124/127 | budget ok")
    print("peak %d live values <= 30 usable registers (%d free)\n" % (st["peak"], 30 - st["peak"]))

    print("| # | row | instruction | writes | reads | live after |")
    print("|---:|---:|---|:--|:--|---:|")
    for i, x in enumerate(S):
        print("| %d | %d | `%s` | %s | %s | %d |" % (
            i + 1, x.row, x.asm,
            ", ".join("%s=%s" % (r, v) for r, v in x.defs) or "-",
            ", ".join("%s(%s)" % (v, r) for r, v in x.uses) or "-",
            len(st["live_after"][i])))

    print("\nplanted contract violations (each must make the verifier fail):")
    res = []

    def A_clobber_while_live():            # row 2 overwrites NY0 (live until row 4)
        M = build()
        for x in M:
            if x.asm == "umulh x0,x16,x9":
                x.defs = [("x23", "V1a")]
        return M

    def B_missing_preservation():          # TBL, a batch-live value, is folded away
        M = build()
        for x in M:
            if x.asm == "add x6,x21,x6,lsl #5":
                x.defs = [("x21", "TBL2")]
        return M

    def C_wrong_output_register():         # NX2 is written to x8, product reads x2
        M = build()
        for x in M:
            if x.asm == "ldp x2,x3,[x19,#16]":
                x.defs = [("x8", "NX2"), ("x3", "NX3")]
        return M

    def D_census_flip():                   # one ledger row deleted
        return [x for x in build() if x.asm != "adds x27,x27,#1"]

    for name, fn in [("A clobber-while-live (row 2 -> NY0 in x23)", A_clobber_while_live),
                     ("B missing-preservation (row 5 -> TBL in x21)", B_missing_preservation),
                     ("C wrong-output-register (NX2 to x8, read x2)", C_wrong_output_register),
                     ("D census flip (length increment deleted)", D_census_flip)]:
        e, s2 = verify(fn(), quiet=True)
        res.append(bool(e))
        print("  %-46s -> %-22s (%d error%s, peak %d) %s" % (
            name, "FAIL as required" if e else "PASSED = CONTROL DEAD",
            len(e), "" if len(e) == 1 else "s", s2["peak"], e[0] if e else ""))
    if not all(res):
        print("\nCONTROL DEAD: the self-audit does not certify")
        return 1
    print("\n4/4 controls fired (3 required contract violations + the census flip). CONTRACT OK.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
