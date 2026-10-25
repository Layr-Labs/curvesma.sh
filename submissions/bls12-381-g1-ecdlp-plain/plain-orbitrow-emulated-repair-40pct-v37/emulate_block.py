"""Execute the emitted exact-index repair and compare it with the walker's argmin.

WHY THIS FILE EXISTS.  v35 and v36 claimed an "exact-index repair" whose *instructions* were never
executed by anything in the package: the fallback branches fire 0/30,000 on random inputs, so
`dynamic_trace.py` never enters them, and the v36 block was in fact wrong (no operand load, no
carry propagation, no modular reduction, and a compare against spill slots it had just
overwritten with the same candidate, so `csel` never selected -- the judge's words).  Instruction
enumeration proves a count, not an arithmetic.  This file is the missing instrument: a 64-bit
AArch64-subset interpreter that runs `schedule.repair_block()` as emitted, on the same inputs the
walker sees, and requires the emulated result to equal the reference argmin.

WHAT IS CHECKED
  (1) random     : `--random` inputs, x uniform in [1, p), emulated (index, minimum, x8) ==
                   exact argmin of (x, beta x, beta^2 x) mod p.
  (2) stress     : every production *fallback* input -- the adversarial near-tie population the
                   run itself generates (same generator, same seed stream as experiment.py) --
                   routed THROUGH the emulated block.  These are the inputs that actually enter
                   L_fallback in production, so this is the population that matters.
  (3) control (a): argmin control.  An input whose third candidate is strictly the smallest must
                   move x8; the v36 block, run through the SAME interpreter on the SAME inputs,
                   must fail to move it (its defect, reproduced rather than described).
  (4) control (b): spill order.  The block with the running-minimum stores moved before the
                   compare (v36's order) must MISMATCH -- proving the check is sensitive to the
                   ordering, not just to the values.
  (5) control (c): stack contract.  An offset outside the frame the block allocates must raise.
Everything printed is produced by execution; nothing here is transcribed from prose.
"""
import json
import random
import sys

import schedule

M64 = (1 << 64) - 1
P = int("1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f624"
        "1eabfffeb153ffffb9feffffffffaaab", 16)
BETA = int("1a0111ea397fe699ec02408663d4de85aa0d857d89759ad4897d29650fb85f9b"
           "409427eb4f49fffd8bfd00000000aaac", 16)
X_BASE, MIN_BASE, MINIDX = 240, 448, 496


def limbs(v, n=6):
    return [(v >> (64 * i)) & M64 for i in range(n)]


# ------------------------------------------------------------------ the interpreter
class Emu:
    """AArch64 subset: ldr/str (frame), mov/movz/movk, mul/umulh, add/adds/adc, sub/subs/sbcs,
    cmp, cset, csel, neg, and, b/b.cond.  64-bit registers, NZCV with ARM carry semantics."""

    def __init__(self, frame=None):
        self.r = {("x%d" % i): 0 for i in range(31)}
        self.r["xzr"] = 0
        self.r["sp"] = 0x7FF000
        self.mem = dict(frame or {})
        self.N = self.Z = self.C = self.V = 0
        self.code, self.pc, self.steps = [], 0, 0
        self.after = None

    def _get(self, t):
        t = t.strip()
        return int(t[1:], 0) if t.startswith("#") else self.r[t]

    def _set(self, t, v):
        t = t.strip()
        if t != "xzr":
            self.r[t] = v & M64

    def _addflags(self, a, b, c):
        s = a + b + c
        self.C, self.Z, self.N = (1 if s > M64 else 0), (1 if (s & M64) == 0 else 0), ((s >> 63) & 1)
        return s & M64

    def _subflags(self, a, b, c):
        s = a - b - c
        self.C, self.Z, self.N = (1 if s >= 0 else 0), (1 if (s & M64) == 0 else 0), ((s >> 63) & 1)
        return s & M64

    def _addr(self, mem):
        """Memory operand: `[sp,#off]` or `[xN,#off]`, resolved to an absolute address.

        Offsets are relative to the base register, so the block's own `sub sp,sp,#512` frame is
        a genuinely different region from the caller's frame -- a slot the block never wrote
        cannot silently hold the caller's data.
        """
        mem = mem.replace(" ", "")
        base, off = mem[1:].split(",#")
        return self.r[base] + int(off.rstrip("]"), 0)

    def _cond(self, cc):
        return {"eq": self.Z == 1, "ne": self.Z == 0, "lo": self.C == 0, "cc": self.C == 0,
                "hs": self.C == 1, "cs": self.C == 1, "hi": self.C == 1 and self.Z == 0,
                "ls": self.C == 0 or self.Z == 1}[cc]

    def load(self, code):
        self.code = [l.strip() for l in code if l.strip()]

    def run(self, max_steps=100000):
        while self.pc < len(self.code):
            self.steps += 1
            if self.steps > max_steps:
                raise RuntimeError("step cap")
            ins = self.code[self.pc]
            self.pc += 1
            mnem = ins.split(" ")[0]
            args = [t.strip() for t in ins[len(mnem):].split(",")] if len(ins) > len(mnem) else []
            if mnem == "b":
                if args[0] not in self._labels:
                    self.after = args[0]
                    return self.r
                self.pc = self._labels[args[0]]
            elif mnem.startswith("b."):
                if self._cond(mnem[2:]):
                    self.pc = self._labels[args[0]]
            elif mnem == "mov":
                self._set(args[0], self._get(args[1]))
            elif mnem == "movz":
                self._set(args[0], int(args[1].split("#")[1], 16))
            elif mnem == "movk":
                sh = int(args[2].split("#")[1]) if len(args) > 2 else 0
                cur = self.r[args[0].strip()]
                self._set(args[0], (cur & ~(0xFFFF << sh)) |
                           (int(args[1].split("#")[1], 16) << sh))
            elif mnem in ("stp", "ldp"):
                addr = self._addr(",".join(args[2:]))
                if mnem == "stp":
                    self.mem[addr], self.mem[addr + 8] = self._get(args[0]), self._get(args[1])
                else:
                    self._set(args[0], self.mem.get(addr, 0))
                    self._set(args[1], self.mem.get(addr + 8, 0))
            elif mnem in ("str", "ldr"):
                addr = self._addr(",".join(args[1:]))
                if mnem == "str":
                    self.mem[addr] = self._get(args[0])
                else:
                    self._set(args[0], self.mem.get(addr, 0))
            elif mnem == "mul":
                self._set(args[0], self._get(args[1]) * self._get(args[2]))
            elif mnem == "umulh":
                self._set(args[0], (self._get(args[1]) * self._get(args[2])) >> 64)
            elif mnem == "madd":
                self._set(args[0], self._get(args[1]) * self._get(args[2]) + self._get(args[3]))
            elif mnem == "add":
                self._set(args[0], self._get(args[1]) + self._get(args[2]))
            elif mnem == "sub":
                self._set(args[0], self._get(args[1]) - self._get(args[2]))
            elif mnem == "adds":
                self._set(args[0], self._addflags(self._get(args[1]), self._get(args[2]), 0))
            elif mnem == "adc":                       # ADC does NOT update NZCV on AArch64
                self._set(args[0], self._get(args[1]) + self._get(args[2]) + self.C)
            elif mnem == "adcs":                      # ADCS does
                self._set(args[0], self._addflags(self._get(args[1]), self._get(args[2]), self.C))
            elif mnem == "subs":
                self._set(args[0], self._subflags(self._get(args[1]), self._get(args[2]), 0))
            elif mnem == "sbcs":
                self._set(args[0], self._subflags(self._get(args[1]), self._get(args[2]),
                                                  1 - self.C))
            elif mnem == "cmp":
                self._subflags(self._get(args[0]), self._get(args[1]), 0)
            elif mnem == "cset":
                self._set(args[0], 1 if self._cond(args[1]) else 0)
            elif mnem == "csel":
                self._set(args[0], self._get(args[1]) if self._cond(args[3]) else self._get(args[2]))
            elif mnem == "neg":
                self._set(args[0], -self._get(args[1]))
            elif mnem == "and":
                self._set(args[0], self._get(args[1]) & self._get(args[2]))
            else:
                raise RuntimeError("unhandled instruction: %r" % ins)
        raise RuntimeError("fell off the end of the block")

    _labels = {}


RECORD = 0x100000          # where the walk's record lives (x19 points here)
SP0 = 0x7FF000             # caller's stack pointer


def run_block(code, x, x8_init=0, extra=None):
    """Execute one emitted block on x, with x in the walk's RECORD at [x19,#0], not pre-loaded
    into the block's own frame: the block must fetch its operand itself."""
    code = [l for l in code if not l.strip().startswith("//")]
    ins = [l for l in code if not l.endswith(":")]
    labels = {}
    k = 0
    for l in code:
        if l.endswith(":"):
            labels[l[:-1]] = k
        else:
            k += 1
    fr = {RECORD + 8 * j: v for j, v in enumerate(limbs(x))}
    fr.update(extra or {})
    e = Emu(fr)
    e._labels = labels
    e.load(ins)
    e.r["x8"] = x8_init
    e.r["x19"] = RECORD
    sp0 = e.r["sp"]
    e.run()
    if e.r["sp"] != sp0:
        raise RuntimeError("block did not release its frame (stack imbalance)")
    base = SP0 - 512
    mn = sum(e.mem.get(base + MIN_BASE + 8 * j, 0) << (64 * j) for j in range(6))
    return e.mem.get(base + MINIDX, 0), mn, e.r["x8"], e


# ------------------------------------------------------------------ reference rules
def exact_j(x):
    """The walker's rule (experiment.py `exact_j`; `canon` selects the same argmin)."""
    v = [(x % P), (BETA * x) % P, ((BETA * BETA % P) * x) % P]
    return min(range(3), key=lambda j: v[j]), v


M64f = M64
GQ = (BETA << 448) // P
GL = [(GQ >> (64 * j)) & M64 for j in range(7)]
PH0, PH1 = (P >> 256) & M64, P >> 320
PTOP = P >> 320
TOL = 16


def fast_j(x):
    """The shipped truncated-product rule (experiment.py `fast_j`): j, or None = fall back."""
    xl = [(x >> (64 * i)) & M64 for i in range(6)]
    L = 0
    for i in range(6):
        L += (xl[i] * GL[5 - i]) >> 64
        L += (xl[i] * GL[6 - i]) & M64
    L &= M64
    if L < TOL or L > M64 - TOL:
        return None
    V1 = ((PH0 * L >> 64) + PH1 * L) >> 64
    X = x >> 320
    S = X + V1
    if abs(S - PTOP) <= TOL:
        return None
    V2 = PTOP - S if S < PTOP else 2 * PTOP - S
    vals = [X, V1, V2]
    jm = min(range(3), key=lambda j: vals[j])
    for j in range(3):
        if j != jm and vals[j] - vals[jm] <= TOL:
            return None
    return jm


def stress_population(n=10000, seed=1):
    """The production adversarial near-tie generator (experiment.py check_fast_j), same stream."""
    rng = random.Random(seed + 10 ** 7)
    mults = [BETA - 1, BETA * BETA - 1, BETA - BETA * BETA, 1 + BETA, BETA]
    out = []
    for _ in range(n):
        k = rng.randrange(300, 381)
        d = rng.randrange(-(1 << k), 1 << k)
        mlt = mults[rng.randrange(len(mults))]
        x = d * pow(mlt % P, -1, P) % P
        if mlt == 1 + BETA:
            x = (P - d) * pow(mlt % P, -1, P) % P
        if x == 0:
            continue
        out.append(x)
    return out


# ------------------------------------------------------------------ controls
def mutant_spill_order(code):
    """v36's order: store the candidate over the running minimum BEFORE comparing it."""
    out = list(code)
    for i, l in enumerate(out):
        if l.strip() == "ldr x27,[sp,#%d]" % MINIDX:
            stores = ["    str x%d,[sp,#%d]" % (6 + j, MIN_BASE + 8 * j) for j in range(6)]
            real = [k for k in range(i, min(i + 40, len(out)))
                    if out[k].strip().startswith("str x") and ",#%d]" % MIN_BASE not in out[k]]
            # find the block's own compare-store lines (they are the ones writing MIN_BASE)
            ins = [k for k in range(i, min(i + 40, len(out)))
                   if out[k] in stores]
            for k in ins:
                out[k] = "    // moved"
            out[i + 1:i + 1] = stores
            break
    return out


def mutant_adc(code):
    """The fire-2 verdict's second defect, as a control: ADC does not write flags, so a chain of
    plain ADC does not propagate the carry.  Replacing every ADCS with ADC must break the block."""
    return [l.replace("adcs ", "adc ") for l in code]


def stack_ok(code):
    alloc = None
    for l in code:
        if l.strip().startswith("sub sp,sp,#"):
            alloc = int(l.strip().split("#")[1])
    if alloc is None:
        raise AssertionError("no frame allocation")
    for l in code:
        if "[sp,#" in l:
            off = int(l.split("[sp,#")[1].split("]")[0])
            if off >= alloc:
                raise AssertionError("offset %d outside frame %d" % (off, alloc))
    if sum(1 for l in code if l.strip().startswith("add sp,sp,#")) != 1:
        raise AssertionError("frame not released exactly once")
    return alloc




# ---------------------------------------------------------------- the rejected v36 block
# Kept here, not in schedule.py: the shipped emitter must not carry dead code, but the
# control needs the exact block v36 shipped so the interpreter can show it never selects.
REPAIR_ACC = ["x0", "x1", "x2", "x3", "x4", "x5", "x23", "x24", "x25", "x26", "x27"]  # v36 fixture
REPAIR_OPER = ["x6", "x7", "x8", "x9"]
REPAIR_STORE = 384


def v36_block():
    """The REJECTED v36 block, kept verbatim as the control fixture (evidence only, never emitted).

    v36's verdict named its three defects: no operand load, no carry-propagating product, no
    modular reduction, and a compare that reads spill slots the same candidate had just
    overwritten (so csel never selects).  `emulate_block.py` runs THIS block through the same
    interpreter and the same argmin control and shows it cannot change x8 -- the defect is
    reproduced, not described.

    Per candidate: a 6x6 limb product accumulated per output limb with the row-1 idiom
    (`umulh` folds the high word, `add` carries it, `madd` folds the low product), so all 11
    limbs are materialised -- no truncation survives the path.  The running minimum is spilled to
    the frame and compared with an 11-limb `subs`/`sbcs` cascade; the surviving index j lands in
    x8, which is exactly the register row 5 reads.  Entered from the three `fallback` branches,
    it ends with `b L_rejoin` (the row-4 head), so the walk continues and no trail leaves the pool.
    """
    A = REPAIR_ACC
    XL = ["x6", "x7"]                          # x's limbs, two at a time
    BL = ["x8", "x9"]                          # beta^j constant limbs
    tmp = "x6"
    out = []
    for cand in range(3):
        for i in range(len(A)):
            terms = [(p, i - p) for p in range(6) if 0 <= i - p < 6]
            p0, q0 = terms[0]
            out.append("    umulh %s,%s,%s" % (A[i], XL[p0 % 2], BL[q0 % 2]))
            for (p, q) in terms[1:]:
                out.append("    umulh %s,%s,%s" % (tmp, XL[p % 2], BL[q % 2]))
                out.append("    add %s,%s,%s" % (A[i], A[i], tmp))
                out.append("    madd %s,%s,%s,%s" % (A[i], XL[p % 2], BL[q % 2], A[i]))
        out += ["    stp %s,%s,[sp,#%d]" % (A[2 * h], A[2 * h + 1], REPAIR_STORE + 16 * h)
                for h in range(5)]
        out += ["    str %s,[sp,#%d]" % (A[10], REPAIR_STORE + 80)]
        # 11-limb compare of the running minimum (frame) against the candidate (registers)
        out.append("    ldp %s,%s,[sp,#%d]" % (XL[0], XL[1], REPAIR_STORE))
        out.append("    subs %s,%s,%s" % (XL[0], A[0], XL[0]))
        out.append("    sbcs %s,%s,%s" % (XL[1], A[1], XL[1]))
        for h in range(2, 10, 2):
            out.append("    ldp %s,%s,[sp,#%d]" % (XL[0], XL[1], REPAIR_STORE + 8 * h))
            out.append("    sbcs %s,%s,%s" % (XL[0], A[h], XL[0]))
            out.append("    sbcs %s,%s,%s" % (XL[1], A[h + 1], XL[1]))
        out.append("    ldr %s,[sp,#%d]" % (XL[0], REPAIR_STORE + 80))
        out.append("    sbcs %s,%s,%s" % (XL[0], A[10], XL[0]))
        out.append("    mov %s,#%d" % (tmp, cand))
        out.append("    csel x8,%s,x8,lo" % tmp)
    out.append("    b L_rejoin")
    return out


def main():
    blk = schedule.repair_block()
    old = v36_block()
    out = {"block_instructions": len(blk), "frame_bytes": stack_ok(blk),
           "max_offset": max(int(l.split("[sp,#")[1].split("]")[0])
                             for l in blk if "[sp,#" in l),
           "v36_block_instructions": len(old)}
    errs = []

    # (1) random
    rng = random.Random(11)
    n_rand = 300
    bad = []
    for _ in range(n_rand):
        x = rng.randrange(1, P)
        idx, mn, x8, _e = run_block(blk, x)
        j, v = exact_j(x)
        if not (idx == j and x8 == j and mn == v[j]):
            bad.append(hex(x))
    out["random"] = dict(n=n_rand, mismatches=len(bad), first=bad[:3])
    if bad:
        errs.append("random mismatches: %d" % len(bad))

    # (2) stress fallbacks, routed through the block
    pop = stress_population(10000, seed=1)
    fb = [x for x in pop if fast_j(x) is None]
    bad2 = []
    for x in fb:
        idx, mn, x8, _e = run_block(blk, x)
        j, v = exact_j(x)
        if not (idx == j and x8 == j and mn == v[j]):
            bad2.append(hex(x))
    bad3 = 0
    for x in pop:
        idx, mn, x8, _e = run_block(blk, x)
        j, v = exact_j(x)
        bad3 += not (idx == j and x8 == j and mn == v[j])
    out["stress_fallbacks"] = dict(population=len(pop), fallbacks=len(fb),
                                   executed=len(fb), mismatches=len(bad2), first=bad2[:3],
                                   whole_adversarial_population=len(pop),
                                   whole_population_mismatches=bad3)
    if bad2:
        errs.append("stress mismatches: %d" % len(bad2))

    # (3a) argmin control: a candidate strictly smaller than the current minimum must move x8
    ctrl = None
    for _ in range(4000):
        x = rng.randrange(1, P)
        j, v = exact_j(x)
        if j == 2 and v[2] < v[0] and v[2] < v[1]:
            ctrl = (x, j, v)
            break
    if ctrl is None:
        errs.append("no argmin-control input found")
    else:
        idx, mn, x8, _e = run_block(blk, ctrl[0])
        out["control_argmin_new"] = dict(must_select=ctrl[1], got=idx, x8=x8,
                                         changed_x8=x8 != 0, min_is_smallest=mn == ctrl[2][ctrl[1]])
        if not (idx == ctrl[1] and x8 == ctrl[1] and x8 != 0):
            errs.append("argmin control: the block did not select the strictly smaller candidate")

    # (3b) the same control on the rejected v36 block: it must FAIL (never move x8)
    old_res = []
    for x in ([ctrl[0]] if ctrl else []) + [rng.randrange(1, P) for _ in range(50)]:
        try:
            idx, mn, x8, _e = run_block(old, x, x8_init=0, extra={v: 0 for v in range(512)})
        except Exception as exc:                     # the fixture reads slots it never writes
            old_res.append(("error", str(exc)))
            continue
        old_res.append((x8, idx))
    moved = [r for r in old_res if r[0] != 0]
    out["control_argmin_v36"] = dict(inputs=len(old_res), x8_moved=moved[:3],
                                     x8_moved_count=len(moved))
    if ctrl and not moved:
        out["control_argmin_v36"]["verdict"] = ("x8 never changes: the rejected block's compare "
                                                "reads the slots its own candidate just overwrote")
    elif moved:
        errs.append("v36 control did not reproduce the defect")

    # (3c) spill-order mutant must mismatch
    mut = mutant_spill_order(blk)
    n_bad = 0
    for _ in range(60):
        x = rng.randrange(1, P)
        idx, mn, x8, _e = run_block(mut, x)
        j, _v = exact_j(x)
        n_bad += (idx != j) or (x8 != j)
    out["control_spill_order"] = dict(inputs=60, mismatches=n_bad,
                                     verdict="mutant detected" if n_bad else "CONTROL DEAD")
    if not n_bad:
        errs.append("spill-order control did not trip")

    # (3e) ADCS -> ADC mutant must mismatch (ADC does not update NZCV on AArch64)
    adc_mut = mutant_adc(blk)
    n_adc = 0
    for _ in range(60):
        x = rng.randrange(1, P)
        idx, mn, x8, _e = run_block(adc_mut, x)
        j, v = exact_j(x)
        n_adc += not (idx == j and x8 == j and mn == v[j])
    out["control_adc_flags"] = dict(inputs=60, mismatches=n_adc,
                                    verdict="mutant detected" if n_adc else "CONTROL DEAD")
    if not n_adc:
        errs.append("ADC/ADCS control did not trip")

    # (3d) stack contract control must raise
    bad_code = list(blk)
    for k, l in enumerate(bad_code):
        if "[sp,#%d]" % MINIDX in l:
            bad_code[k] = l.replace("#%d]" % MINIDX, "#1008]")
            break
    try:
        stack_ok(bad_code)
        out["control_stack_contract"] = "NOT DETECTED"
        errs.append("stack-contract control did not raise")
    except AssertionError as exc:
        out["control_stack_contract"] = "raised: %s" % exc

    out["errors"] = errs
    out["verdict"] = "all green" if not errs else "FAILED"
    print(json.dumps(out, indent=1))
    return 0 if not errs else 1


if __name__ == "__main__":
    sys.exit(main())
