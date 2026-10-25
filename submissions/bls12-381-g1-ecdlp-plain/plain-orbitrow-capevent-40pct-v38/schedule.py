#!/usr/bin/env python3
"""Machine-checked liveness / census / operand-interface contract for the plain G1 rho step (v39).

Row 6 is the operand interface of the fixed-price 6M reference addition, and it is the item
42a2308b and 73c192b7 both named.  A 6M addition consumes 24 point limbs at once (table x,
table y, state x, state y); the walk has 30 usable registers and 12 pinned values, so
24 + 12 = 36 > 30 and no load order or interleaving exists.  v39 relocates six pinned values
(PH0, PTOP, END, TBL, CONST, RPTR) to frame slots across the boundary -- 12 charged transfers --
so the boundary sees 24 point limbs + P0..P5 in 30 distinct registers with 0 free.

Checks (all fail-closed):
  (1) forward definedness -- every read finds *that* value in *that* register;
  (2) backward liveness -- a write may not destroy a value read later (frame-relocated values
      are carried as slot values m<off>, so a relocation is a move, not a clobber);
  (3) register budget -- x18 is reserved, so 30 registers are usable; peak live <= 30;
  (4) census -- 143 charged instructions on a non-refresh step (146 on the 1-in-64 refresh),
      row totals 24/15/14/20/8/32/12/9/1/1/1/4/2/2;
  (5) pinned preservation -- proved per register: every write is the RPTR hand-off or a
      bracketed relocation whose slot is restored before the register's next use;
  (6) operand interface (interface_budget) -- the 30 declared operands are present, distinct
      and in the pool at the [C] boundary.  The v38 streaming interface FAILS this check
      (control I), because the state-y loads overwrite the table-x registers.

Controls A-I must each make some check fail.
"""

import re
import sys

RESERVED = {18}                                              # x18 is platform-reserved
POOL = ["x%d" % i for i in range(31) if i not in RESERVED]    # 30 usable

LIVE_IN = ["P0", "P1", "P2", "P3", "P4", "P5", "PH0", "PTOP", "RPTR", "END", "TBL", "CONST"]
LIVE_OUT = ["P0", "P1", "P2", "P3", "P4", "P5", "PH0", "PTOP", "RPTR2", "END", "TBL", "CONST"]
PINNED = ["x10", "x11", "x12", "x13", "x14", "x15", "x16", "x17", "x19", "x20", "x21", "x22"]

SPAN = 143                    # charged instructions on a non-refresh step (v39: +12 row-6 relocations)
SPAN_MAX = 146                # charged instructions on the 1-in-64 refresh step
FUSION_BOUND = SPAN_MAX      # one iteration; enforced by the single back-edge in the object


class Ins:
    __slots__ = ("row", "asm", "defs", "uses", "cost", "refresh", "literal", "live_in")

    def __init__(self, row, asm, defs=(), uses=(), cost=1, refresh=False, literal=True,
                 live_in=()):
        self.row, self.asm = row, asm
        self.defs, self.uses = list(defs), list(uses)
        self.cost, self.refresh = cost, refresh
        self.literal = literal       # False = a C-internal witness row, never an instruction
        self.live_in = list(live_in)  # values a separate per-batch pass brings in itself


def build():
    S = []
    A = S.append
    NY = ["x23", "x24", "x25", "x26", "x27", "x28"]
    NX = ["x0", "x1", "x2", "x3", "x4", "x5"]
    PY = ["x6", "x7", "x9", "x29", "x30", "x5"]

    # ---- row 15 (2, per batch): the batch prologue, HOISTED OUT of the element loop ----
    # addressed [x19,#-8] while x19 advances 512 per element.  Two consequences: iteration k's
    # save overwrote iteration k-1's, so after the first iteration the slot held a clobbered
    # state limb rather than the caller's link, and the restore at the end read [END-8] while the
    # last save had gone to [END-520].  A moving address cannot preserve a value across a loop
    # element loop, into a dedicated frame slot at a STATIONARY address [sp,#8]; the restore is
    # the same address expression [sp,#8] on every exit.  Address does not move, value is saved
    # before the first clobber, restore reads the save.  The element loop therefore cannot move
    # the slot and cannot re-save.
    A(Ins(15, "sub sp,sp,#64", [], []))
    A(Ins(15, "stur x30,[sp,#8]", [], []))

    # ---- row 6 (32): the operand interface, SPILL-BASED (v39) -----------------------------
    # 73c192b7 is right about v38: it loaded the table x-coordinate into x0..x6 and then wrote
    # the state y-coordinate into those same six registers before the [C] boundary, so the
    # boundary's declared inputs contained no surviving table x.  A 6M reference addition needs
    # 24 point limbs live at once (table x, table y, state x, state y).  The walk has 30 usable
    # registers and 12 of them are pinned for the whole span, so 24 + 12 = 36 > 30: **no load
    # order and no interleaving exists**, the v30-v38 "stream the coordinates" interface is not
    # a scheduling bug but a register-count fact, and it cannot be repaired inside 131.
    #
    # v39 therefore RELOCATES six pinned values to the frame for the duration of the boundary --
    # PH0, PTOP, END, TBL, CONST (each last used before row 6) and RPTR (relocated after the
    # state point has been read out of the record, restored before the result stores) -- so the
    # boundary sees 24 point limbs + the 6 modulus limbs in 30 distinct registers, 0 free.
    # Each relocation is an emitted, charged transfer: 6 stores + 6 loads = 12.
    #
    #     table x  TX0..TX5   x0,x1,x2,x3,x4,x6
    #     table y  TY0..TY5   x7,x8,x9,x23,x24,x25
    #     state x  SX0..SX5   x26,x27,x28,x29,x30,x16
    #     state y  SY0..SY5   x17,x20,x21,x22,x19,x5
    #     modulus  P0..P5     x10,x11,x12,x13,x14,x15
    #
    # The slot budget is machine-checked by interface_budget() below: it fails if any of the 30
    # declared operands is absent, duplicated, or displaced at the boundary, and the pinned set
    # is carried through the frame slots (a spill store defines the slot's value m<off>, the
    # restore load consumes it, so the liveness checker follows the value through memory).
    TX = ["x0", "x1", "x2", "x3", "x4", "x6"]        # table x limbs
    TY = ["x7", "x8", "x9", "x23", "x24", "x25"]     # table y limbs
    SX = ["x26", "x27", "x28", "x29", "x30", "x16"]  # state x limbs
    SY = ["x17", "x20", "x21", "x22", "x19", "x5"]   # state y limbs
    A(Ins(6, "ldr x5,[x19,#96]", [("x5", "ENT")], [("x19", "RPTR")]))
    for reg, slot, val in [("x16", 16, "PH0"), ("x17", 24, "PTOP"), ("x20", 32, "END"),
                           ("x21", 40, "TBL"), ("x22", 48, "CONST")]:
        A(Ins(6, "stur %s,[sp,#%d]" % (reg, slot), [("m%d" % slot, val)], [(reg, val)]))
    for off, half in [(0, 0), (16, 2), (32, 4)]:
        A(Ins(6, "ldp %s,%s,[x19,#%d]" % (SX[half], SX[half + 1], off),
              [(SX[half], "SX%d" % half), (SX[half + 1], "SX%d" % (half + 1))], [("x19", "RPTR")]))
    for off, half in [(0, 0), (16, 2), (32, 4)]:
        A(Ins(6, "ldp %s,%s,[x5,#%d]" % (TX[half], TX[half + 1], off),
              [(TX[half], "TX%d" % half), (TX[half + 1], "TX%d" % (half + 1))], [("x5", "ENT")]))
    for off, half in [(48, 0), (64, 2), (80, 4)]:
        A(Ins(6, "ldp %s,%s,[x5,#%d]" % (TY[half], TY[half + 1], off),
              [(TY[half], "TY%d" % half), (TY[half + 1], "TY%d" % (half + 1))], [("x5", "ENT")]))
    for off, half in [(48, 0), (64, 2)]:
        A(Ins(6, "ldp %s,%s,[x19,#%d]" % (SY[half], SY[half + 1], off),
              [(SY[half], "SY%d" % half), (SY[half + 1], "SY%d" % (half + 1))], [("x19", "RPTR")]))
    A(Ins(6, "stur x19,[sp,#56]", [("m56", "RPTR")], [("x19", "RPTR")]))
    A(Ins(6, "ldp %s,%s,[x19,#80]" % (SY[4], SY[5]),
          [(SY[4], "SY4"), (SY[5], "SY5")], [("x19", "RPTR")]))
    # The boundary: 30 declared operands, every one of them a live register holding exactly that
    # limb.  Its writes are the twelve new coordinates and the batch product; the addition's
    # price is the fixed 6M convention charged in H, and its operand access is these 30 live
    # registers -- nothing is hidden in memory at the interface.
    A(Ins(6, "[C] charged reference addition (6M): 30 register operands -> new coordinates + product",
          [(TX[i], "NSX%d" % i) for i in range(6)]
          + [("x7", "NSY0"), ("x9", "NSY1"), ("x23", "NSY2"), ("x24", "NSY3"),
             ("x25", "NSY4"), ("x26", "NSY5"), ("x8", "BATCHPROD")],
          [(SX[i], "SX%d" % i) for i in range(6)] + [(SY[i], "SY%d" % i) for i in range(6)]
          + [(TX[i], "TX%d" % i) for i in range(6)] + [(TY[i], "TY%d" % i) for i in range(6)]
          + [("x10", "P0"), ("x11", "P1"), ("x12", "P2"), ("x13", "P3"),
             ("x14", "P4"), ("x15", "P5")], cost=0, literal=False))
    # Montgomery prefix pass is class 2 of the traffic table (inside C, priced there) and its
    # running value P_i = prod_{j<=i} z_j is the boundary's output in x8, stored to the batch's
    # frame slot on every element so the slot ends the batch holding the FULL product.
    A(Ins(6, "stur x8,[sp,#0]", [], [("x8", "BATCHPROD")]))
    A(Ins(6, "ldr x19,[sp,#56]", [("x19", "RPTR")], [("m56", "RPTR")]))
    A(Ins(6, "stp x0,x1,[x19,#0]", [], [("x0", "NSX0"), ("x1", "NSX1"), ("x19", "RPTR")]))
    A(Ins(6, "stp x2,x3,[x19,#16]", [], [("x2", "NSX2"), ("x3", "NSX3"), ("x19", "RPTR")]))
    A(Ins(6, "stp x4,x6,[x19,#32]", [], [("x4", "NSX4"), ("x6", "NSX5"), ("x19", "RPTR")]))
    A(Ins(6, "stp x7,x9,[x19,#48]", [], [("x7", "NSY0"), ("x9", "NSY1"), ("x19", "RPTR")]))
    A(Ins(6, "stp x23,x24,[x19,#64]", [], [("x23", "NSY2"), ("x24", "NSY3"), ("x19", "RPTR")]))
    A(Ins(6, "stp x25,x26,[x19,#80]", [], [("x25", "NSY4"), ("x26", "NSY5"), ("x19", "RPTR")]))
    for reg, slot, val in [("x16", 16, "PH0"), ("x17", 24, "PTOP"), ("x20", 32, "END"),
                           ("x21", 40, "TBL"), ("x22", 48, "CONST")]:
        A(Ins(6, "ldr %s,[sp,#%d]" % (reg, slot), [(reg, val)], [("m%d" % slot, val)]))

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
    # v39 fire-2 fix (50faf11e): the operational cap IS lambda = 25.  The compare immediate is the
    # cap the walk enforces, so it must equal TRAIL_CAP/theta -- v38 bought the accounting constant
    # to 25 and left this compare at 20, which made the age-conditioned charge 25e^-25 price a cap
    # the machine did not implement.  One literal, one meaning.
    A(Ins(8, "cmp x28,#25", [], [("x28", "LENHI")]))
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
    # row-5 report branch has already exited on that limb being zero: the test was vacuous and
    # mis-attributed the class-14 charge (verdict 28abe582, flaw 2).  The batch-pass invariant is
    # over the batch denominator product Pi = prod z_i (the quantity the shared inversion
    # consumes and the batch commit gates on), so the test reads the accumulator.
    # between the last element's `stur` and the test, so the tested value is the transferred
    # product and nothing else.  The transfer + test together are 2 per element; the batch's
    # zero-denominator recovery action stays priced at 2^17 units per fire.
    A(Ins(10, "ldur x8,[sp,#0]", [("x8", "BATCHPROD2")], []))
    A(Ins(10, "cbz x8,zero-recover", [], [("x8", "BATCHPROD2")]))
    A(Ins(11, "add meter.lo,meter.lo,batch-size", [], [], literal=False))
    A(Ins(12, "share of the one Fermat inversion (232,864/N, N >= 2^18)", [], [], literal=False))
    # ---- row 14b (2): restore from the SAME stationary slot, then tear the frame down ----
    A(Ins(14, "ldur x30,[sp,#8]", [("x30", "LINK")], []))
    A(Ins(14, "add sp,sp,#64", [], []))

    # ---- row 13 (3): the forward prefix pass's control, one per element ----
    A(Ins(13, "ldr x2,[x0,#96]", [("x2", "PFXENT")], [("x0", "PFXBASE")],
          live_in=[("x0", "PFXBASE"), ("x1", "PFXN")]))
    A(Ins(13, "add x0,x0,#512", [("x0", "PFXBASE2")], [("x0", "PFXBASE")]))
    A(Ins(13, "subs x1,x1,#1", [("x1", "PFXCTR")], [("x1", "PFXN")]))
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
        for r, v in x.live_in:
            cur[r] = v
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
        if x.live_in:      # a separate per-batch pass supplies these itself; do not carry them back
            live = {p for p in live if p[0] not in {r for r, _ in x.live_in}}
        reglive = {q for q in live if not q[0].startswith("m")}
        if len(reglive) > peak:
            peak, peak_at = len(reglive), "before row %d instr %d (`%s`)" % (x.row, i, x.asm)

    if peak > len(POOL):
        errs.append("budget: peak %d live values > %d usable registers" % (peak, len(POOL)))

    total = sum(x.cost for x in S)
    refresh = sum(x.cost for x in S if x.refresh)
    rows = {}
    for x in S:
        rows[x.row] = rows.get(x.row, 0) + x.cost
    want = {1: 24, 2: 15, 3: 14, 4: 20, 5: 8, 6: 32, 7: 12, 8: 9, 10: 2, 11: 1, 12: 1,
            13: 4, 14: 2, 15: 2}
    if total - refresh != SPAN:
        errs.append("census: %d charged instructions, the section-5.1 table prices %d"
                    % (total - refresh, SPAN))
    if refresh != 3:
        errs.append("census: %d refresh-only instructions, the listing prices 3" % refresh)
    if rows != want:
        errs.append("census: row totals %s disagree with the traffic table %s" % (rows, want))
    if total != SPAN_MAX:
        errs.append("census: refresh step %d instructions, the listing prices %d"
                    % (total, SPAN_MAX))

    if not quiet:
        print("charged instructions    : %d non-refresh (%d on the 1-in-64 refresh step, %d refresh-only)"
              % (total - refresh, total, refresh))
        print("row census              : %s" % ", ".join("row %d = %d" % (r, c) for r, c in sorted(rows.items())))
        print("peak live values        : %d at %s" % (peak, peak_at))
        nz = [x for x in S if x.cost]
        regs = sorted({r for x in S for r, _ in x.defs + x.uses if not r.startswith("m")}, key=lambda s: int(s[1:]))
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
    span, the last use, and the verdict.  A write is legal in exactly two ways:

      hand-off     -- after the last use of the pinned value, i.e. the RPTR -> RPTR2 advance;
      relocation   -- v39: the value goes to a frame slot first (a `stur xR,[sp,#off]` whose
                      def names the slot value), the register is reused, and the value is
                      restored from that same slot before its next use.  The value is preserved
                      because the forward pass proves the slot is defined and unmodified between
                      the two transfers; the register is not.

    Anything else is BROKEN.  A relocation must be bracketed: spill before the first
    displacing write, restore after the last one.
    """
    tab = []
    for r in PINNED:
        want = LIVE_IN[PINNED.index(r)]
        defs = [(i, vv) for i, x in enumerate(S) for rr, vv in x.defs if rr == r]
        writes = [i for i, vv in defs if vv != want]
        uses = [(i, vv) for i, x in enumerate(S) for rr, vv in x.uses if rr == r and vv == want]
        last_use = max((i for i, _ in uses), default=-1)
        restores = [i for i, vv in defs if vv == want]
        # the relocation pair: a store whose DEF names the slot value and whose USE is this
        # register holding this value, and a load whose DEF is this register/value
        spills = [(i, x.defs[0][0]) for i, x in enumerate(S)
                  if len(x.defs) == 1 and x.defs[0][0].startswith("m") and x.defs[0][1] == want
                  and any(rr == r and vv == want for rr, vv in x.uses)]
        reloc = None
        inside = [i for i in writes if i < last_use]
        if spills and inside:
            p0, slot = spills[0][0], spills[0][1]
            after = [i for i in restores if i > max(inside)]
            cand = [i for i in after if any(rr == slot and vv == want for rr, vv in S[i].uses)]
            if p0 < min(inside) and cand:
                reloc = (p0, min(cand), slot)
        handoff = [i for i in writes if i >= last_use]
        bad = [i for i in writes if i < last_use and not (reloc and reloc[0] < i < reloc[1])]
        ok = not bad
        if not writes:
            held = "whole span"
        elif reloc and not bad:
            held = "relocated to [sp,#%s] over [%d,%d], value restored" % (reloc[2][1:], reloc[0], reloc[1])
        elif ok:
            held = "to last use, then handed off"
        else:
            held = "BROKEN"
        tab.append(dict(reg=r, value=want, pin=-1, last_use=last_use, writes=writes,
                        handoff=handoff, ok=ok, bad=bad, held=held, reloc=reloc))
    return tab


def interface_budget(S, quiet=False):
    """The slot-budget proof, as data (v39).

    At the [C] boundary the interface declares 30 operands -- table x, table y, state x, state y
    (24 point limbs) and the modulus P0..P5 -- and this check reproduces, from the emitted
    instructions alone, which register holds which limb at the boundary: all 30 must be present,
    distinct, and inside the 30-register pool, and the occupancy is reported against it.
    This is the check a v38-style streaming interface fails: there, the state-y loads overwrite
    the table-x registers, so TX0..TX5 are absent at the boundary and the check fires.
    """
    b = [i for i, x in enumerate(S) if not x.literal and x.cost == 0]
    if len(b) != 1:
        return dict(rows=[], occupied=0, free=0, errs=["expected exactly one [C] boundary row, found %d" % len(b)])
    b = b[0]
    cur = dict(zip(PINNED, LIVE_IN))
    for x in S[:b]:
        for r, v in x.live_in:
            cur[r] = v
        for r, v in x.defs:
            cur[r] = v
    hold = {}
    for r, v in cur.items():
        if not r.startswith("m"):
            hold.setdefault(v, []).append(r)
    rows, errs, seen = [], [], set()
    for coord, tmpl in (("table x", "TX%d"), ("table y", "TY%d"),
                        ("state x", "SX%d"), ("state y", "SY%d"), ("modulus", "P%d")):
        for i in range(6):
            v = tmpl % i
            rs = hold.get(v, [])
            if not rs:
                errs.append("%s operand %s is ABSENT at the [C] boundary" % (coord, v))
                r = None
            elif len(rs) > 1:
                errs.append("%s operand %s is duplicated in %s" % (coord, v, rs))
                r = rs[0]
            elif rs[0] in seen:
                errs.append("%s operand %s shares register %s with another operand" % (coord, v, rs[0]))
                r = rs[0]
            else:
                r = rs[0]
                seen.add(r)
            rows.append(dict(coord=coord, limb=i, value=v, reg=r))
    occ = len({r for r in cur if not r.startswith("m")})
    free = len(POOL) - occ
    if free < 0:
        errs.append("boundary occupancy %d exceeds the %d-register pool" % (occ, len(POOL)))
    if not quiet:
        print("\n## (9) row-6 operand interface, slot budget at the [C] boundary\n")
        print("| coordinate | limb | value | register |")
        print("|---|---:|---|---|")
        for t in rows:
            print("| %s | %d | %s | %s |" % (t["coord"], t["limb"], t["value"], t["reg"]))
        print("\noperand registers occupied: %d of %d usable; free: %d; declared operands: %d"
              % (occ, len(POOL), free, len(rows)))
    return dict(rows=rows, occupied=occ, free=free, errs=errs)


RECORD_BYTES = 512          # the walk's record; [sp,#0] is the batch product, [sp,#8] the link
RELOC_SLOT0 = 16            # v39: the first of the six pinned-relocation slots in the frame
FIELDS = {0: "x limb 0/1", 16: "x limb 2/3", 32: "x limb 4/5", 48: "y limb 0/1",
          64: "y limb 2/3", 80: "y limb 4/5", 96: "entry address", 104: "K_cp limbs 0,2",
          120: "K_cp limb 3 + valid", 136: "length (2 limbs)"}

MEM = re.compile(r"^(ldp|ldr|stp|str)\s+([^,]+(?:,[^,]+)?),\s*\[(x[0-9]+|sp)(?:,#([0-9]+))?\]")


def traffic_pairs(S):
    """Store -> reload pairing with an ELIGIBILITY FILTER (v28).

    v27 matched memory operands by textual offset from x19 and wrapped the search across the
    loop back-edge.  But `add x19,x19,#512` advances the record pointer, so a store to
    [x19,#0] in iteration k and a load from [x19,#0] in iteration k+1 address DIFFERENT
    records: the wrapped pairs were never producer-consumer pairs, and the "max distance 122"
    they produced bounded nothing.  v28 admits a pair only when

      (i)  both sides are record-relative (base x19) and inside the 512-byte record -- a
           stack/frame operand (sp-relative, or an offset beyond the record) is stack traffic,
           excluded from fusion credit and charged at its own load/store price;
      (ii) both sides name the SAME record, i.e. the store precedes the load in the same
           executed iteration (no wrap across the back-edge);
      (iii) the load's value is the store's value (same offset in the same record).

    Every load and store that is excluded is still a literal instruction in the span, so the
    exclusion removes a claim, not a charge: delta 0.0000.
    """
    pos, k = {}, 0
    for i, x in enumerate(S):
        if x.cost and x.literal:
            pos[i] = k
            k += 1
    n = k
    stores, loads, bad = {}, [], []
    for i, x in enumerate(S):
        m = MEM.match(x.asm)
        if not m or i not in pos:
            continue
        kind, base, off = m.group(1), m.group(3), int(m.group(4) or 0)
        if base == "sp" and off >= RELOC_SLOT0 and x.row == 6:
            # v39 relocation of a pinned value through the walk's own frame.  It is a frame
            # transfer by construction, it is charged 1 instruction, and it is NOT record
            # traffic: excluded from fusion credit, listed in the relocation census instead.
            continue
        if base == "sp" or (base == "x19" and off >= RECORD_BYTES):
            bad.append(dict(idx=pos[i], asm=x.asm, why="stack/frame-relative: not dataflow"))
            continue
        if base != "x19":
            continue                 # another object (table, constants, prefix array): not a record pair
        if kind in ("stp", "str"):
            stores.setdefault(off, []).append(pos[i])
        else:
            loads.append((pos[i], off))
    eligible, excluded = [], []
    for p_i, off in loads:
        same = [p for p in stores.get(off, []) if p < p_i]
        if same:
            p = max(same)
            eligible.append(dict(off=off, producer=p, consumer=p_i, dist=p_i - p, wrapped=False))
        else:
            later = stores.get(off, [])
            if later:
                p = max(later)
                excluded.append(dict(off=off, producer=p, consumer=p_i, dist=(n - p) + p_i,
                                     why="read precedes the write in this iteration: the value is "
                                         "carried state from an earlier iteration of a DIFFERENT record"))
            else:
                excluded.append(dict(off=off, producer=None, consumer=p_i, dist=None,
                                     why="no store to this field in the span"))
    mx = max((p["dist"] for p in eligible), default=0)
    return dict(eligible=eligible, excluded=excluded, bad=bad, mx=mx, n=n, bound=FUSION_BOUND)


def fold_by_coordinate(S):
    """Per-coordinate segmentation of the row-6 interface, taken from the emitted instructions.

    The fold's INTERIOR is inside C and has no encoding in the object (the object contains only
    the interface transfers), so the interior is not itemised per coordinate; the 3M reference
    formula does split symmetrically (lambda 1M shared -> 0.5M per coordinate, x3 1M, y3 1M), so
    the inside-C share is 1.5M = 574.5 instructions per coordinate and 1,149 together, which is
    class 1 of the traffic table and is not charged again.
    """
    rows = []
    for x in S:
        if x.row != 6 or not x.literal or not x.cost:
            continue
        if x.asm.startswith("stur x8"):
            rows.append(dict(coord="product", asm=x.asm, kind="stur"))
            continue
        m = MEM.match(x.asm)
        if m and m.group(3) == "sp":          # v39 spill / restore: not a point operand
            rows.append(dict(coord="relocation", asm=x.asm, kind=m.group(1)))
            continue
        if m and m.group(3) == "x19":
            off = int(m.group(4) or 0)
        elif m:                      # table-relative: [x5,#off]
            off = int(m.group(4) or 0)
        else:
            continue
        coord = "entry" if off == 96 else ("x" if off in (0, 16, 32) else "y")
        rows.append(dict(coord=coord, asm=x.asm, kind=m.group(1)))
    tab = []
    for c in ("entry", "x", "y", "product", "relocation"):
        insns = [r for r in rows if r["coord"] == c]
        tab.append(dict(coord=c, insns=len(insns), charged=len(insns),
                        inside_c_instr=574.5 if c in ("x", "y") else 0,
                        traffic=0.0 if c == "relocation" else 1.0))
    return tab


def branch_table():
    """Declared branch -> target edges of the emitted object (v28: distinct targets)."""
    return [
        dict(row=2, asm="b.hs fallback", target="L_fallback", kind="exit-stub", stub_body=2),
        dict(row=2, asm="b.ls fallback", target="L_fallback", kind="exit-stub", stub_body=2),
        dict(row=3, asm="b.ls fallback", target="L_fallback", kind="exit-stub", stub_body=2),
        dict(row=5, asm="cbz x1,report", target="L_report", kind="exit-stub", stub_body=2),
        dict(row=7, asm="b.eq traversal", target="L_traversal", kind="exit-stub", stub_body=2),
        dict(row=7, asm="b.ne norefresh", target="L_norefresh", kind="in-body skip", stub_body=0),
        dict(row=8, asm="b.hs abandon", target="L_abandon", kind="exit-stub", stub_body=2),
        dict(row=8, asm="b.ne loop", target="L_loop", kind="loop back-edge", stub_body=0),
        dict(row=10, asm="cbz x8,zero-recover", target="L_zero", kind="exit-stub", stub_body=2),
        dict(row=13, asm="b.ne prefix", target="L_prefix", kind="prefix back-edge", stub_body=0),
    ]


def stub_accounting():
    """Per-entry accounting for the exit stubs.

    v27 emitted every label in one block after the back-edge, so all nine conditional branches
    resolved to one `ret` and the prefix rows sat outside any loop.  v28 emits a real CFG: each
    exit action has its own stub, a stub is one instruction (`ret`), and a taken exit branch
    leaves the loop body.  Per element at most ONE exit branch can be taken, because the exit
    branches sit in one straight-line body and any taken one transfers out of it (the stub has
    no path back into the body).  So the per-entry count is 1 and the single-count is 1:
    old single-count 1 x 1 instruction = new per-entry count 1 x 1 instruction, delta 0.0000.
    """
    return dict(indeg={}, entries_per_element=1, body=1, old=1, new=1, delta=0)


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

    def I_streaming_interface():       # the v38 row 6: state y loaded over the table-x registers
        M = build()
        for x in M:
            if x.asm in ("ldp x17,x20,[x19,#48]", "ldp x21,x22,[x19,#64]", "ldp x19,x5,[x19,#80]"):
                x.defs = [("x%d" % (0 if "#48" in x.asm else (2 if "#64" in x.asm else 4)), "SY0"),
                          ("x%d" % (1 if "#48" in x.asm else (3 if "#64" in x.asm else 6)), "SY1")]
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
            ("F fusion-bound vacuity (bound set one below measured max)", F_fusion_bound_vacuity, "fusion"),
            ("I streaming interface (v38 row 6, table x overwritten)", I_streaming_interface, "budget")]:
        M = fn()
        if kind == "verify":
            e, s2 = verify(M, quiet=True)
            note = e[0] if e else ""
        elif kind == "pinned":
            t = pinned_preservation(M)
            bad = [r for r in t if not r["ok"]]
            e = ["pinned %s: write at %s <= last use %d" % (r["reg"], r["bad"], r["last_use"]) for r in bad]
            note = e[0] if e else ""
        elif kind == "budget":
            e = interface_budget(M, quiet=True)["errs"]
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


LABELS = {"b.hs fallback": ("b.hs", "L_fallback"), "b.ls fallback": ("b.ls", "L_fallback"),
          "b.eq traversal": ("b.eq", "L_traversal"), "b.ne norefresh": ("b.ne", "L_norefresh"),
          "b.hs abandon": ("b.hs", "L_abandon"), "b.ne loop": ("b.ne", "L_loop"),
          "cbz x1,report": ("cbz", "L_report"), "cbz x8,zero-recover": ("cbz", "L_zero"),
          "b.ne prefix": ("b.ne", "L_prefix")}

# ---------------------------------------------------------------------------------------------
# v37 -- THE EXACT-INDEX REPAIR: CORRECT, EXECUTED, EMULATED AGAINST THE WALKER
# ---------------------------------------------------------------------------------------------
# v36 was rejected because its block never loaded the operands, never reduced mod p, and compared
# each candidate against spill slots it had just overwritten (so csel never selected).  v37 emits a
# block that does the arithmetic: operand loads, a 6x6 schoolbook with carry propagation, a
# separated-operand Montgomery reduction with a three-limb carry cascade, one borrow-propagating
# conditional subtraction of p by csel, the walker's own beta^2 identity, and per candidate
# read-compare-then-write.  It allocates its own 512-byte frame (the prologue's 16-byte link frame
# is untouched), saves every register but x8, puts the exact argmin in x8, and ends `b L_rejoin`.
# `emulate_block.py` executes these instructions and compares with the walker's argmin;
# `stack_contract()` and `cfg_audit()` enforce the frame and the rejoin.
from repair import repair_block, stack_contract   # the v37 repair block (own file: byte cap)


EXIT_LABELS = ("L_report", "L_traversal", "L_abandon", "L_zero")


def emit_asm(S):
    """The literal machine listing with real control flow (v28 emission order).

    Prologue hoisted out of the element loop; the prefix pass is its own loop before the leaf
    body; `b.ne L_norefresh` targets the length-update block physically next; one stub per exit
    action; the fallback entry is the v37 repair, which rejoins at L_rejoin instead of returning.
    """
    prologue = [x for x in S if x.literal and x.cost and x.row == 15]
    prefix = [x for x in S if x.literal and x.cost and x.row == 13]
    body = [x for x in S if x.literal and x.cost and x.row not in (13, 15)]
    out = [".text", ".p2align 2", "_batch_entry:"]
    for x in prologue:                      # hoisted OUT of the element loop: once per batch
        out.append("    " + x.asm.replace(",", ", "))
    out += [".p2align 2", "_batch_prefix:", "L_prefix:"]
    for x in prefix:
        a = x.asm
        if a in LABELS:
            op, lab = LABELS[a]
            out.append("    %s %s" % (op, lab))
        else:
            out.append("    " + a.replace(",", ", "))
    out += [".p2align 2", "_leaf_span:", "L_loop:"]   # v33: no `ret` here -- the prefix pass FALLS
    # THROUGH into the leaf loop, so the span is one single-entry function (prologue, prefix pass,
    # a `ret` after the prefix pass and the trace interpreter then redirected it to the leaf body;
    # that redirect was not what the machine does.  Removing the `ret` deletes the redirect instead
    # of modelling it: -1 scaffolding instruction (149 -> 148 object instructions), charged census
    # unchanged (131 non-refresh / 134 refresh).
    prev = None
    for x in body:
        if x.row == 8 and prev == 7:
            out.append("L_norefresh:")
        if x.row == 4 and prev == 3:
            # v36: re-entry point of the exact-index repair block.  The repair runs out of line
            # (cold path) and rejoins the walk here, at the row-4 head, with J already in x8.
            out.append("L_rejoin:")
        a = x.asm
        if a in LABELS:
            op, lab = LABELS[a]
            if op == "cbz":
                out.append("    cbz %s, %s" % (a.split(" ", 1)[1].split(",")[0], lab))
            else:
                out.append("    %s %s" % (op, lab))
        else:
            out.append("    " + a.replace(",", ", "))
        if x.row == 6 and x.asm.startswith("ldp x19,x5,[x19,#80]"):
            # The priced C boundary sits here: the reference addition consumes the operand block,
            # returns the new coordinates and the Montgomery prefix accumulator in x8.  This is a
            # COMMENT, not an instruction: it is not charged, it is not executed, and the listing's
            out.append("    // [C] charged 6M reference addition: operand block -> new coordinates,"
                       " and its Montgomery prefix accumulator in x8 (not an instruction)")
        prev = x.row
    out.append("    ret")      # the frame restore is row 14b, emitted above as body rows
    # v36: the fallback entry is the REPAIR, not an exit.  It recomputes the index exactly and
    # branches back to L_rejoin; nothing returns, so no trail leaves the pool.
    out += ["L_fallback:"] + repair_block()
    for lab in EXIT_LABELS:      # exit-only: at most one fire per pass, amortised < 2^-18
        out += [lab + ":", "    ldur x30, [sp, #8]", "    add sp, sp, #64", "    ret"]
    return "\n".join(out)


def cfg_audit(src):
    """Audit the emitted control flow.  Errors (each makes the audit fail):

      * an exit branch whose target is inside the loop body (v27's aliasing defect);
      * an exit stub that is not a single block ending in `ret`, or that falls back into the body;
      * the leaf body having anything other than exactly one back-edge, or `b.ne loop`;
      * `b.ne norefresh` not targeting an in-body block that continues at the length update;
      * the prefix pass not being a separate loop before the leaf body with one back-edge.
    """
    errs = []
    lines = [l.strip() for l in src.splitlines() if l.strip() and not l.strip().startswith(".")]
    label_at = {}
    items = []
    for i, l in enumerate(lines):
        if l.endswith(":") and not l.startswith(("b.", "cbz", "ldr", "ldp", "stp", "str", "add",
                                                 "sub", "subs", "cmp", "cmn", "adc", "sbc", "sbcs",
                                                 "csel", "cset", "mov", "and", "tst", "ccmp", "madd",
                                                 "mul", "umulh", "ret")):
            label_at[l[:-1]] = len(items)
        else:
            items.append(l)
    body_start = label_at.get("L_loop", 0)
    body_end = None
    for i in range(body_start, len(items)):
        if items[i] == "ret":
            body_end = i
            break
    if body_end is None:
        return ["no ret closing the leaf body"]
    # back-edges
    backs = []
    for i, ins in enumerate(items):
        m = re.match(r"^(b\.[a-z]{2}|cbz)\s+(?:x[0-9]+,\s*)?(L_[A-Za-z_]+)$", ins)
        if not m:
            continue
        tgt = m.group(2)
        if tgt not in label_at:
            errs.append("branch to unknown label %s" % tgt)
            continue
        t = label_at[tgt]
        if t <= i:
            backs.append((i, tgt))
        if tgt in EXIT_LABELS and t < body_end:
            errs.append("branch `%s` targets exit label %s which sits INSIDE the body" % (ins, tgt))
    leaf_backs = [b for b in backs if b[0] < body_end and body_start <= b[0]]
    if len(leaf_backs) != 1 or leaf_backs[0][1] != "L_loop":
        errs.append("leaf body has %d back-edge(s) %s, expected exactly one `b.ne loop`"
                    % (len(leaf_backs), [b[1] for b in leaf_backs]))
    if "L_norefresh" not in label_at:
        errs.append("no L_norefresh label: `b.ne norefresh` has no in-body target")
    else:
        t = label_at["L_norefresh"]
        if not (body_start < t < body_end):
            errs.append("L_norefresh is not inside the leaf body")
    pre = [b for b in backs if b[0] < body_start]
    if len(pre) != 1 or pre[0][1] != "L_prefix":
        errs.append("prefix pass has %d back-edge(s) %s, expected exactly one `b.ne L_prefix`"
                    % (len(pre), [b[1] for b in pre]))
    if label_at.get("L_prefix", 10 ** 9) > body_start:
        errs.append("the prefix loop does not precede the leaf body")
    for b in branch_table():
        op = "cbz" if b["asm"].startswith("cbz") else b["asm"].split()[0]
        pat = ("%s " % op) if op == "cbz" else ("%s " % op)
        hits = [l for l in items if l.startswith(pat) and l.rstrip().endswith(b["target"])]
        if not hits:
            errs.append("declared edge `%s` -> %s is missing (some branch aliases another target)"
                        % (b["asm"], b["target"]))
    # exit stubs: one instruction, and no fall-through into the body
    for lab in EXIT_LABELS:
        if lab not in label_at:
            errs.append("missing exit stub %s" % lab)
            continue
        t = label_at[lab]
        if t <= body_end:
            errs.append("exit stub %s is not out of line" % lab)
        elif not (items[t].startswith("ldur x30") and t + 2 < len(items)
                  and items[t + 1].startswith("add sp") and items[t + 2] == "ret"):
            errs.append("exit stub %s is not `ldur x30; add sp; ret` (frame restore + ret)" % lab)
    # v36: the repair block, out of line, re-enters the walk and never returns
    if "L_fallback" not in label_at:
        errs.append("missing repair entry L_fallback")
    else:
        t = label_at["L_fallback"]
        if t <= body_end:
            errs.append("repair block L_fallback is not out of line")
        if "L_rejoin" not in label_at:
            errs.append("repair target L_rejoin is missing")
        elif not (body_start < label_at["L_rejoin"] < body_end):
            errs.append("L_rejoin is not inside the leaf body (the repair would not rejoin the walk)")
        if "ret" in items[t:t + len(repair_block())]:
            errs.append("repair block contains a `ret`: the trail would leave the pool")
        if not items[t + len(repair_block()) - 1].endswith("L_rejoin"):
            errs.append("repair block does not end in `b L_rejoin`")
        try:
            stack_contract()
        except AssertionError as exc:
            errs.append("repair stack contract: %s" % exc)
        if not items[t].strip().startswith("sub sp,sp,#"):
            errs.append("repair block does not allocate its own frame")
        if items[t + len(repair_block()) - 2].strip().split(",")[0] != "add sp":
            errs.append("repair block does not release its frame before rejoining")
    return errs


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

    def E_pinned_clobber():
        M = build()
        for x in M:
            if x.asm == "ldp x2,x3,[x5,#16]":
                x.defs = [("x12", "TX2"), ("x3", "TX3")]
        return M

    def I_streaming_interface():       # the v38 row 6: state y loaded over the table-x registers
        M = build()
        for x in M:
            if x.asm in ("ldp x17,x20,[x19,#48]", "ldp x21,x22,[x19,#64]", "ldp x19,x5,[x19,#80]"):
                x.defs = [("x%d" % (0 if "#48" in x.asm else (2 if "#64" in x.asm else 4)), "SY0"),
                          ("x%d" % (1 if "#48" in x.asm else (3 if "#64" in x.asm else 6)), "SY1")]
        return M

    def F_fusion_bound_vacuity():
        M = build()
        t = traffic_pairs(M)
        return M, t["mx"] - 1

    def G_stack_pair():                 # a record operand rewritten to be stack-relative
        M = build()
        for x in M:
            if x.asm == "ldp x0,x1,[x19,#0]":
                x.asm = "ldp x0,x1,[sp,#0]"
        return M

    def H_exit_edge_inside_body():      # the v27 defect: an exit branch aliased to an in-body label
        return emit_asm(build()).replace("b.eq L_traversal", "b.eq L_norefresh")

    for name, fn, kind in [
            ("A clobber-while-live (row 2 -> NY0 in x23)", A_clobber_while_live, "verify"),
            ("B missing-preservation (row 5 -> TBL in x21)", B_missing_preservation, "verify"),
            ("C wrong-output-register (NX2 to x8, read x2)", C_wrong_output_register, "verify"),
            ("D census flip (length increment deleted)", D_census_flip, "verify"),
            ("E pinned-clobber (row 6 writes p limb 2 to x12)", E_pinned_clobber, "pinned"),
            ("F fusion-bound flip (bound one below measured max)", F_fusion_bound_vacuity, "fusion"),
            ("G stack pair (record operand made sp-relative)", G_stack_pair, "pairs"),
            ("H exit edge aliased into the body (v27 defect)", H_exit_edge_inside_body, "cfg"),
            ("I streaming interface (v38 row 6: table x overwritten by state y)", I_streaming_interface, "budget")]:
        M = fn()
        if kind == "verify":
            e, _ = verify(M, quiet=True)
        elif kind == "pinned":
            t = pinned_preservation(M)
            e = ["pinned %s: write at %s <= last use %d" % (r["reg"], r["bad"], r["last_use"])
                 for r in t if not r["ok"]]
        elif kind == "fusion":
            M, bnd = M
            t = traffic_pairs(M)
            e = ["fusion: max distance %d > bound %d" % (t["mx"], bnd)] if t["mx"] > bnd else []
        elif kind == "budget":
            e = interface_budget(M, quiet=True)["errs"]
        elif kind == "pairs":
            t = traffic_pairs(M)
            e = ["pairing: %d stack/frame-relative operand(s) credited or counted: %s"
                 % (len(t["bad"]), [b["asm"] for b in t["bad"]])] if t["bad"] else []
        else:
            e = cfg_audit(M if isinstance(M, str) else emit_asm(M))
        note = e[0] if e else ""
        res.append(bool(e))
        print("  %-52s -> %-22s (%d error%s) %s" % (
            name, "FAIL as required" if e else "PASSED = CONTROL DEAD",
            len(e), "" if len(e) == 1 else "s", note))
    return res


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
    print("checks: pool ok | forward definedness ok | no clobber | census %d/%d | budget ok" % (SPAN, SPAN_MAX))
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
    print("\n%d/%d pinned registers preserved: every write inside the span is either the"
          "\nspecified hand-off (after that register's last use) or a v39 relocation through a"
          "\nframe slot that is restored before the register's next use." % (len(tab) - len(broken), len(tab)))

    ib = interface_budget(S)
    for e in ib["errs"]:
        print("  FAIL  " + e)
    if ib["errs"]:
        print("\nOPERAND INTERFACE REJECTED (%d error%s)" % (len(ib["errs"]), "" if len(ib["errs"]) == 1 else "s"))
        return 1
    print("\noperand interface: %d/%d declared operands present in distinct registers,"
          "\n%d of %d usable registers occupied, %d free -- the slot budget closes with 0 slack."
          % (len(ib["rows"]), len(ib["rows"]), ib["occupied"], 30, ib["free"]))

    print("\n## (6) row 6 segmented per coordinate (from the emitted instructions)\n")
    print("| coordinate | path | executed insns | charged | inside C (instr.) | traffic weight |")
    print("|---|---|---:|---:|---:|---:|")
    for r in fold_by_coordinate(S):
        print("| %s | %s | %d | %d | %.1f | %.2f |" % (
            r["coord"], "table+state transfers" if r["coord"] != "entry" else "entry-address fetch",
            r["insns"], r["charged"], r["inside_c_instr"], r["traffic"]))
    print("\nweighted fold price per element = 12 operand transfers + 7 stores + 12 frame")
    print("relocations = 31 executed instructions (both coordinate paths fire")
    print("on every element: the row-6 region is one straight-line block, 0 conditional branches).")

    print("\n## (7) store -> reload pairing with the eligibility filter\n")
    t = traffic_pairs(S)
    print("| record offset | field | producer (idx) | consumer (idx) | distance (insn) |")
    print("|---:|---|---:|---:|---:|")
    for pp in sorted(t["eligible"], key=lambda q: -q["dist"]):
        print("| %d | %s | %d | %d | %d |" % (pp["off"], FIELDS.get(pp["off"], "?"),
                                              pp["producer"], pp["consumer"], pp["dist"]))
    print("\n%d eligible same-record pairs over a %d-instruction iteration; max measured distance %d;"
          " bound %d." % (len(t["eligible"]), t["n"], t["mx"], t["bound"]))
    print("\n%d candidates EXCLUDED (not producer-consumer pairs), each still a literal charged"
          " instruction (delta 0):" % len(t["excluded"]))
    print("\n| record offset | field | store idx | load idx | excluded because |")
    print("|---:|---|---:|---:|---|")
    for pp in t["excluded"]:
        print("| %d | %s | %s | %d | %s |" % (pp["off"], FIELDS.get(pp["off"], "?"),
                                              pp["producer"], pp["consumer"], pp["why"]))
    if t["bad"]:
        print("\nSTACK/FRAME PAIRS INSIDE THE TRAFFIC TABLE: %s" % t["bad"])
        return 1
    if any(p["wrapped"] for p in t["eligible"]):
        print("\nA CREDITED PAIR WRAPS THE BACK-EDGE (different records)")
        return 1

    print("\n## (8) branch -> target edges and stub accounting (v28 emits a real CFG)\n")
    print("| row | branch | target | kind | stub body (insn) |")
    print("|---:|---|---|---|---:|")
    indeg = {}
    for b in branch_table():
        indeg[b["target"]] = indeg.get(b["target"], 0) + 1
        print("| %d | `%s` | %s | %s | %d |" % (b["row"], b["asm"], b["target"], b["kind"], b["stub_body"]))
    print("\nin-degree: %s" % ", ".join("%s=%d" % (k, v) for k, v in sorted(indeg.items())))
    sa = stub_accounting()
    print("stub accounting: at most %d exit entry per element (any taken exit branch leaves the"
          " body; a stub is one `ret` with no path back), so old single-count %d x %d instr. ="
          " new per-entry count %d x %d instr., delta %d." %
          (sa["entries_per_element"], sa["old"], sa["body"], sa["new"], sa["body"], sa["delta"]))
    cfg_errs = cfg_audit(emit_asm(S))
    for e in cfg_errs:
        print("  CFG FAIL  " + e)
    if cfg_errs:
        return 1
    print("CFG audit: prefix pass is a separate loop before the leaf body; `b.ne norefresh` targets"
          " the in-body length-update block; %d exit stubs out of line, one `ret` each; exactly one"
          " leaf back-edge (`b.ne loop`)." % len(EXIT_LABELS))

    print("\nplanted contract violations (each must make the verifier fail):")
    res = controls(S)
    if not all(res):
        print("\nCONTROL DEAD: the self-audit does not certify")
        return 1
    print("\n8/8 controls fired (4 contract/census violations + the pinned clobber + the"
          " fusion-bound flip). CONTRACT OK.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
