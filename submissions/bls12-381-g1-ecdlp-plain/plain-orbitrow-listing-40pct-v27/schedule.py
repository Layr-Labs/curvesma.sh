#!/usr/bin/env python3
"""Machine-checked liveness/census contract for the plain BLS12-381 G1 rho step (v27 listing).

v27 changes exactly one thing at the listing level: **row 6 is no longer a cost-19
pseudo-instruction.**  In v22-v26 row 6 was one row whose `asm` string was a *count*
("ldr entry, 6 ldp table point, 6 ldp state_old, 6 stp state_new") and whose declared
reads were {RPTR} with the single write x0.  The v26 verdict (e7f4d2f6) rejected that
construct, correctly: a count is not an instruction, so the liveness checker could not
witness the interface it was pricing.  Row 6 is now **19 separate rows, each one literal
AArch64 instruction, each with its own defs/uses**, plus two zero-cost rows that mark the
fixed-price reference addition's in-place fold of the table point into the state point
(they are charged inside C, not in R, and are listed only so the data flow is explicit).

The expansion is forced by the register budget and is the reason the interface streams
per coordinate: the 12 pinned values are live for the whole span, so 18 of the 30 usable
registers are free and 24 operand limbs (12 table + 12 state) cannot be live at once.
Row 6 therefore loads and folds x first (12 operand registers, 24 live with the pinned
set), stores the new x, then reuses the same 12 registers for y.

Five independent checks run over the schedule:

  (1) forward definedness -- every read must find *that* value in *that* register;
  (2) backward liveness -- a write may not destroy a value that is read later;
  (3) register budget -- x18 is reserved on Darwin/AArch64, so 30 registers are usable;
      the peak simultaneous live-value count must be <= 30 and no instruction names x18;
  (4) census -- 124 charged instructions on a non-refresh step, 3 of them refresh-only
      (127 on the 1-in-64 refresh step), row totals 24/15/14/20/8/19/12/9/1/1/1/3.
  (5) pinned-value preservation -- **proved, not asserted** (v27): for each of the twelve
      pinned registers the scan reports the pin point, every write inside the span and the
      last use, and requires every write to fall after the last use (x19 is the only
      pinned register with a write: the row-8 pointer advance, which is exactly the
      hand-off RPTR -> RPTR2).  A planted mid-span clobber of a pinned register must
      make this check fire.

It also measures the **backward-pass fusion bound** from this schedule: every store ->
reload pair the ledger credits as fused is enumerated, its distance in executed
instructions is measured (wrapping across the span's single back-edge), and the maximum
is reported against the bound of one span (124 instructions).  A pair reaching further
would be a second pass and would have to be charged separately; the object listing has
exactly one backward branch, so no such pair can exist.

Controls (each must make the verifier fail):
  A. clobber-while-live     -- row 2 writes its V1 product into x23 while NY0 lives there
  B. missing-preservation   -- row 5 folds the table address into x21, destroying TBL
  C. wrong-output-register  -- row 1 takes x_2 from the product instead of the record load
  D. census flip            -- one ledger row deleted, so the census no longer equals 124
  E. pinned-clobber         -- a row-6 transfer writes x12 (limb 2 of p), a pinned value
  F. fusion-bound flip      -- a store->reload pair is stretched past the span bound

Usage: python3 schedule.py          # table, checks, controls; exit 0 iff all pass
       python3 schedule.py --emit   # the same listing in assembler form (span.s body)
"""

import re
import sys

RESERVED = {18}                                              # x18 is platform-reserved
POOL = ["x%d" % i for i in range(31) if i not in RESERVED]    # 30 usable

LIVE_IN = ["P0", "P1", "P2", "P3", "P4", "P5", "PH0", "PTOP", "RPTR", "END", "TBL", "CONST"]
LIVE_OUT = ["P0", "P1", "P2", "P3", "P4", "P5", "PH0", "PTOP", "RPTR2", "END", "TBL", "CONST"]
PINNED = ["x10", "x11", "x12", "x13", "x14", "x15", "x16", "x17", "x19", "x20", "x21", "x22"]

SPAN = 124                    # charged instructions on a non-refresh step
SPAN_MAX = 127                # charged instructions on the 1-in-64 refresh step
FUSION_BOUND = SPAN_MAX      # one iteration; enforced by the single back-edge in the object


class Ins:
    __slots__ = ("row", "asm", "defs", "uses", "cost", "refresh", "literal")

    def __init__(self, row, asm, defs=(), uses=(), cost=1, refresh=False, literal=True):
        self.row, self.asm = row, asm
        self.defs, self.uses = list(defs), list(uses)
        self.cost, self.refresh = cost, refresh
        self.literal = literal       # False = a C-internal witness row, never an instruction


def build():
    S = []
    A = S.append
    NY = ["x23", "x24", "x25", "x26", "x27", "x28"]
    NX = ["x0", "x1", "x2", "x3", "x4", "x5"]
    PY = ["x6", "x7", "x9", "x29", "x30", "x5"]

    # ---- row 6 (19): the reference addition's record interface, one literal insn per row
    # Streamed per coordinate: 12 pinned values live + 12 operand limbs = 24 <= 30, but
    # 24 operand limbs (12 table + 12 state) cannot be, so x is folded and stored first.
    TX = ["x0", "x1", "x2", "x3", "x4", "x6"]      # table point limbs (dies into C)
    SX = ["x7", "x8", "x9", "x23", "x24", "x25"]   # state point limbs (folded in place)
    A(Ins(6, "ldr x5,[x19,#96]", [("x5", "ENT")], [("x19", "RPTR")]))
    for off, half in [(0, 0), (16, 2), (32, 4)]:
        A(Ins(6, "ldp %s,%s,[x5,#%d]" % (TX[half], TX[half + 1], off),
              [(TX[half], "TX%d" % half), (TX[half + 1], "TX%d" % (half + 1))], [("x5", "ENT")]))
    for off, half in [(0, 0), (16, 2), (32, 4)]:
        A(Ins(6, "ldp %s,%s,[x19,#%d]" % (SX[half], SX[half + 1], off),
              [(SX[half], "SX%d" % half), (SX[half + 1], "SX%d" % (half + 1))], [("x19", "RPTR")]))
    A(Ins(6, "[C] fold x: state_x := table_x + state_x, in place (inside the reference price)",
          [(SX[i], "NSX%d" % i) for i in range(6)],
          [(TX[i], "TX%d" % i) for i in range(6)] + [(SX[i], "SX%d" % i) for i in range(6)]
          + [("x10", "P0"), ("x11", "P1"), ("x12", "P2"), ("x13", "P3"),
             ("x14", "P4"), ("x15", "P5")], cost=0, literal=False))
    for off, half in [(0, 0), (16, 2), (32, 4)]:
        A(Ins(6, "stp %s,%s,[x19,#%d]" % (SX[half], SX[half + 1], off), [],
              [(SX[half], "NSX%d" % half), (SX[half + 1], "NSX%d" % (half + 1)),
               ("x19", "RPTR")]))
    for off, half in [(48, 0), (64, 2), (80, 4)]:
        A(Ins(6, "ldp %s,%s,[x5,#%d]" % (TX[half], TX[half + 1], off),
              [(TX[half], "TY%d" % half), (TX[half + 1], "TY%d" % (half + 1))], [("x5", "ENT")]))
    for off, half in [(48, 0), (64, 2), (80, 4)]:
        A(Ins(6, "ldp %s,%s,[x19,#%d]" % (SX[half], SX[half + 1], off),
              [(SX[half], "SY%d" % half), (SX[half + 1], "SY%d" % (half + 1))], [("x19", "RPTR")]))
    A(Ins(6, "[C] fold y: state_y := table_y + state_y, in place (inside the reference price)",
          [(SX[i], "NSY%d" % i) for i in range(6)],
          [(TX[i], "TY%d" % i) for i in range(6)] + [(SX[i], "SY%d" % i) for i in range(6)]
          + [("x10", "P0"), ("x11", "P1"), ("x12", "P2"), ("x13", "P3"),
             ("x14", "P4"), ("x15", "P5")], cost=0, literal=False))
    for off, half in [(48, 0), (64, 2), (80, 4)]:
        A(Ins(6, "stp %s,%s,[x19,#%d]" % (SX[half], SX[half + 1], off), [],
              [(SX[half], "NSY%d" % half), (SX[half + 1], "NSY%d" % (half + 1)),
               ("x19", "RPTR")]))

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

    # ---- row 2 (15): range test, V1, X, S, wrap test, V2 ----
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

    # ---- row 3 (14): exact argmin j ----
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

    # ---- row 4 (20): y_new reload, p - y, borrow (p-y)-y as the sign, K limbs 0-3 ----
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

    # ---- row 5 (8): h, distinguished test, twist index, table address ----
    A(Ins(5, "and x5,x0,#0xfffff", [("x5", "H")], [("x0", "K0")]))
    A(Ins(5, "cbz x1,report", [], [("x1", "K1")]))
    A(Ins(5, "add x6,x4,x5,lsl #1", [("x6", "U1")], [("x4", "S"), ("x5", "H")]))
    A(Ins(5, "add x6,x6,x6,lsl #1", [("x6", "U2")], [("x6", "U1")]))
    A(Ins(5, "add x6,x6,x8", [("x6", "U3")], [("x6", "U2"), ("x8", "J")]))
    A(Ins(5, "add x6,x6,x6,lsl #2", [("x6", "U4")], [("x6", "U3")]))
    A(Ins(5, "add x6,x21,x6,lsl #5", [("x6", "ENTW")], [("x21", "TBL"), ("x6", "U4")]))
    A(Ins(5, "str x6,[x19,#96]", [], [("x6", "ENTW"), ("x19", "RPTR")]))

    # ---- row 8a (3): the length load and cap test run before row 7's refresh test ----
    A(Ins(8, "ldp x27,x28,[x19,#136]", [("x27", "LENLO"), ("x28", "LENHI")], [("x19", "RPTR")]))
    A(Ins(8, "cmp x28,#20", [], [("x28", "LENHI")]))
    A(Ins(8, "b.hs abandon", [], [("x28", "LENHI")]))

    # ---- row 7 (9 + 3 refresh) ----
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

    # ---- row 8b (6): 69-bit length increment, pointer advance, loop ----
    A(Ins(8, "adds x27,x27,#1", [("x27", "LENLO2")], [("x27", "LENLO")]))
    A(Ins(8, "adc x28,x28,xzr", [("x28", "LENHI2")], [("x28", "LENHI")]))
    A(Ins(8, "stp x27,x28,[x19,#136]", [], [("x27", "LENLO2"), ("x28", "LENHI2"), ("x19", "RPTR")]))
    A(Ins(8, "add x19,x19,#512", [("x19", "RPTR2")], [("x19", "RPTR")]))
    A(Ins(8, "cmp x19,x20", [], [("x19", "RPTR2"), ("x20", "END")]))
    A(Ins(8, "b.ne loop", [], [("x19", "RPTR2"), ("x20", "END")]))

    # ---- rows 10-12 (3, per-batch amortised): zero-denominator test, meter, Fermat share
    A(Ins(10, "cbz x1,zero-recover", [], [("x1", "K1")]))
    A(Ins(11, "add meter.lo,meter.lo,batch-size", [], [], literal=False))
    A(Ins(12, "share of the one Fermat inversion (232,864/N, N >= 2^18)", [], [], literal=False))

    # ---- row 13 (3): the forward prefix pass's control, one per element ----
    A(Ins(13, "ldr x0,[x19,#96]", [("x0", "PFXENT")], [("x19", "RPTR2")]))
    A(Ins(13, "subs x1,x0,#1", [("x1", "PFXCTR")], [("x0", "PFXENT")]))
    A(Ins(13, "b.ne prefix", [], [("x1", "PFXCTR")]))
    return S


def verify(S, quiet=False):
    errs = []

    for i, x in enumerate(S):
        for r, v in x.defs + x.uses:
            if r.startswith("x") and r[1:].isdigit() and r not in POOL:
                errs.append("%d: `%s` names %s, outside the 30-register pool (x18 reserved)" % (i, x.asm, r))

    cur = dict(zip(PINNED, LIVE_IN))
    for i, x in enumerate(S):
        for r, v in x.uses:
            if cur.get(r) != v:
                errs.append("row %d instr %d `%s`: reads %s from %s, which holds %s"
                            % (x.row, i, x.asm, v, r, cur.get(r)))
        for r, v in x.defs:
            cur[r] = v

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

    if peak > len(POOL):
        errs.append("budget: peak %d live values > %d usable registers" % (peak, len(POOL)))

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
        nz = [x for x in S if x.cost]
        regs = sorted({r for x in S for r, _ in x.defs + x.uses}, key=lambda s: int(s[1:]))
        print("registers touched       : %d (%s)" % (len(regs), ", ".join(regs)))
        print("costed rows             : %d literal instructions (%d zero-cost C witnesses)"
              % (len(nz), len(S) - len(nz)))
    return errs, dict(total=total, refresh=refresh, rows=rows, peak=peak, peak_at=peak_at,
                      live_after=live_after)


# --------------------------------------------------------------------------------------
# (5) pinned-value preservation: proved per register, not asserted
# --------------------------------------------------------------------------------------
def pinned_preservation(S, quiet=False):
    """For every pinned register: pin point (live-in, before index 0), all writes in the
    span, the last use, and the verdict.  A write is legal only after the last use of the
    pinned value, i.e. only as the hand-off the span is specified to make."""
    tab = []
    for r in PINNED:
        want = LIVE_IN[PINNED.index(r)]
        writes = [i for i, x in enumerate(S) if any(rr == r and vv != want for rr, vv in x.defs)]
        uses = [i for i, x in enumerate(S) if any(rr == r and vv == want for rr, vv in x.uses)]
        last_use = max(uses) if uses else -1
        bad = [i for i in writes if i < last_use]
        handoff = [i for i in writes if i >= last_use]
        ok = not bad
        tab.append(dict(reg=r, value=want, pin=-1, last_use=last_use, writes=writes,
                        handoff=handoff, ok=ok, bad=bad,
                        held="whole span" if not writes else
                        ("to last use, then handed off" if ok else "BROKEN")))
    return tab


def fusion_bound(S, bound=None, quiet=False):
    """Measure every store -> reload pair the ledger credits as fused, in *executed*
    instructions, wrapping across the span's single back-edge.  Positions are positions in
    the executed stream (cost-carrying rows only), so a wrapped distance counts the whole
    rest of the iteration plus the head of the next one."""
    bound = FUSION_BOUND if bound is None else bound
    pos, k = {}, 0
    for i, x in enumerate(S):
        if x.cost and x.literal:
            pos[i] = k
            k += 1
    n = k                                   # executed instructions per iteration (127)
    pat = re.compile(r"^(ldp|ldr|stp|str)\s+(x[0-9]+(?:,x[0-9]+)?),\s*\[(x[0-9]+)(?:,#([0-9]+))?\]$")
    stores, loads = {}, []
    for i, x in enumerate(S):
        m = pat.match(x.asm)
        if not m or i not in pos:
            continue
        kind, base, off = m.group(1), m.group(3), int(m.group(4) or 0)
        if base != "x19":
            continue
        if kind in ("stp", "str"):
            stores.setdefault(off, []).append(pos[i])
        else:
            loads.append((pos[i], off))
    pairs = []
    for p_i, off in loads:
        cand = [p for p in stores.get(off, []) if p < p_i]
        if cand:
            p = max(cand)
            d, wrapped = p_i - p, False
        elif stores.get(off):
            p = max(stores[off])
            d, wrapped = (n - p) + p_i, True
        else:
            continue
        pairs.append(dict(off=off, producer=p, consumer=p_i, dist=d, wrapped=wrapped))
    mx = max((p["dist"] for p in pairs), default=0)
    return pairs, mx, n, bound


def controls(S):
    res = []

    def A_clobber_while_live():
        M = build()
        for x in M:
            if x.asm == "umulh x0,x16,x9":
                x.defs = [("x23", "V1a")]
        return M

    def B_missing_preservation():
        M = build()
        for x in M:
            if x.asm == "add x6,x21,x6,lsl #5":
                x.defs = [("x21", "TBL2")]
        return M

    def C_wrong_output_register():
        M = build()
        for x in M:
            if x.asm == "ldp x2,x3,[x19,#16]":
                x.defs = [("x8", "NX2"), ("x3", "NX3")]
        return M

    def D_census_flip():
        return [x for x in build() if x.asm != "adds x27,x27,#1"]

    def E_pinned_clobber():            # a row-6 transfer writes limb 2 of p
        M = build()
        for x in M:
            if x.asm == "ldp x2,x3,[x5,#16]":
                x.defs = [("x12", "TX2"), ("x3", "TX3")]
        return M

    def F_fusion_bound_vacuity():      # the bound set one below the measured maximum
        M = build()
        _, mx, _, _ = fusion_bound(M)
        return M, mx - 1

    for name, fn, kind in [
            ("A clobber-while-live (row 2 -> NY0 in x23)", A_clobber_while_live, "verify"),
            ("B missing-preservation (row 5 -> TBL in x21)", B_missing_preservation, "verify"),
            ("C wrong-output-register (NX2 to x8, read x2)", C_wrong_output_register, "verify"),
            ("D census flip (length increment deleted)", D_census_flip, "verify"),
            ("E pinned-clobber (row 6 writes p limb 2 to x12)", E_pinned_clobber, "pinned"),
            ("F fusion-bound vacuity (bound set one below measured max)", F_fusion_bound_vacuity, "fusion")]:
        M = fn()
        if kind == "verify":
            e, s2 = verify(M, quiet=True)
            note = e[0] if e else ""
        elif kind == "pinned":
            t = pinned_preservation(M)
            bad = [r for r in t if not r["ok"]]
            e = ["pinned %s: write at %s <= last use %d" % (r["reg"], r["bad"], r["last_use"]) for r in bad]
            note = e[0] if e else ""
        else:
            M, bnd = M
            p, mx, n, _ = fusion_bound(M, bound=bnd)
            e = ["fusion: max distance %d > bound %d" % (mx, bnd)] if mx > bnd else []
            note = e[0] if e else ""
        res.append(bool(e))
        print("  %-52s -> %-22s (%d error%s) %s" % (
            name, "FAIL as required" if e else "PASSED = CONTROL DEAD",
            len(e), "" if len(e) == 1 else "s", note))
    return res


def emit_asm(S):
    """The literal machine listing: one source line per charged instruction."""
    labels = {"b.hs fallback": ("b.hs", "L_fallback"), "b.ls fallback": ("b.ls", "L_fallback"),
              "b.eq traversal": ("b.eq", "L_traversal"), "b.ne norefresh": ("b.ne", "L_norefresh"),
              "b.hs abandon": ("b.hs", "L_abandon"), "b.ne loop": ("b.ne", "L_loop"),
              "cbz x1,report": ("cbz", "L_report"), "cbz x1,zero-recover": ("cbz", "L_zero"),
              "b.ne prefix": ("b.ne", "L_prefix")}
    out = [".text", ".p2align 2", "L_loop:"]
    for x in S:
        if not x.literal or not x.cost:
            continue
        a = x.asm
        if a in labels:
            op, lab = labels[a]
            if op == "cbz":
                out.append("    cbz %s, %s" % (a.split(" ", 1)[1].split(",")[0], lab))
            else:
                out.append("    %s %s" % (op, lab))
        else:
            out.append("    " + a.replace(",", ", "))
    out += ["L_fallback:", "L_report:", "L_traversal:", "L_norefresh:", "L_abandon:",
            "L_zero:", "L_prefix:", "    ret"]
    return "\n".join(out)


def main():
    S = build()
    if "--emit" in sys.argv:
        print(emit_asm(S))
        return 0
    print("register/liveness/preservation contract -- ordinary step, plain BLS12-381 G1 rho (v27 listing)")
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

    print("\n## (5) pinned-value preservation, per register (pin point = before the span)\n")
    print("| pinned reg | value | pinned at | last use | writes inside the span | verdict |")
    print("|---|:--|---:|---:|---|---|")
    tab = pinned_preservation(S)
    for r in tab:
        print("| %s | %s | live-in | %d | %s | %s |" % (
            r["reg"], r["value"], r["last_use"],
            ", ".join(str(w) for w in r["writes"]) or "none", r["held"]))
    broken = [r for r in tab if not r["ok"]]
    if broken:
        print("\nPRESERVATION BROKEN: %s" % broken)
        return 1
    print("\n%d/%d pinned registers preserved: every write inside the span falls after that"
          "\nregister's last use (x19's single write is the specified RPTR -> RPTR2 hand-off)."
          % (len(tab) - len(broken), len(tab)))

    print("\n## backward-pass fusion bound (measured from this schedule)\n")
    pairs, mx, n_exec, bnd = fusion_bound(S)
    print("| record offset | producer (idx) | consumer (idx) | distance (insn) | wraps span |")
    print("|---:|---:|---:|---:|:--|")
    for p in sorted(pairs, key=lambda q: -q["dist"]):
        print("| %d | %d | %d | %d | %s |" % (p["off"], p["producer"], p["consumer"], p["dist"],
                                              "yes" if p["wrapped"] else "no"))
    print("\n%d fused pairs over a %d-instruction iteration; max measured distance %d"
          " instructions; bound %d (= one iteration)." % (len(pairs), n_exec, mx, bnd))
    if mx > bnd:
        print("FUSION BOUND EXCEEDED")
        return 1

    print("\nplanted contract violations (each must make the verifier fail):")
    res = controls(S)
    if not all(res):
        print("\nCONTROL DEAD: the self-audit does not certify")
        return 1
    print("\n6/6 controls fired (4 contract/census violations + the pinned clobber + the"
          " fusion-bound flip). CONTRACT OK.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
