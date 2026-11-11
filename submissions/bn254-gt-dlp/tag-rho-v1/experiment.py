"""Tag-traced Frobenius-class rho on BN254 GT (proof.md, machine.md). Standard library only.
K: AArch64 listings (leaves, fused Fp2, add4sums, c1dot) executed on a simulator vs Python integers.
A: BN254 identities (tower, Frobenius, keys, canonical element, least-key invariance, tag identity).
C: production walk in the real GT, tag-traced and direct in lock step, every step compared.
B: the complete algorithm in toy fields of the same tower.  F: useless-cycle law, escape charge.
D: exact ledger and score. Toy runs are evidence at their own scale only.
"""
import json
import math
import random
import sys
import time
from decimal import Decimal, ROUND_CEILING, getcontext
from fractions import Fraction
from pathlib import Path

T_START = time.time()
cfg = json.load(sys.stdin)
SEED = cfg.get("seed", 1)
KERNEL_CASES = cfg.get("kernel_cases", 300)
PROD_STEPS = cfg.get("production_steps", 3000)
CHUNK_L = 7
PROD_SECONDS = min(600, cfg.get("prod_seconds", 6))
Q_TOY = cfg.get("q", 151)
TRIALS = cfg.get("trials", 60)
R_TOY = cfg.get("R", 32)
SQ8 = cfg.get("sq8", 0)
DP_BITS = cfg.get("dp_bits", 3)
WALKS = cfg.get("walks", 4)
BUDGET = cfg.get("seconds", 6)
SHOW_LISTINGS = cfg.get("listings", False)
rng = random.Random(SEED)
rng_f = random.Random(SEED ^ 0x5EED)   # separate stream for the fused-leaf tests: the other parts
                                         # then replay the random stream of the base package exactly

base = Path("/snapshot/files")
if not (base / "target-profile.json").is_file():
    base = Path(__file__).resolve().parent
profile = json.loads((base / "target-profile.json").read_text())
p = int(profile["field_modulus"], 16)
r = int(profile["subgroup_order"], 16)
G_PROFILE = [int(c, 16) for c in profile["generator_g"]["coefficients"]]
M64 = (1 << 64) - 1
RMONT = 1 << 256
R_FACT = {2: 28, 3: 2, 13: 1, 29: 1, 983: 1, 11003: 1, 237073: 1, 405928799: 1,
          1670836401704629: 1, 13818364434197438864469338081: 1}

# =============================================================== Part K: AArch64 listings
# Register use: x0 = result pointer, x1 = a pointer, x2 = b pointer, x3..x17 scratch.
# Only caller-saved registers x0..x17 are written, so no spills are needed. xzr reads 0.
# Constants block at label K: p0, p1, p2, p3, n0 = -m^-1 mod 2^64 (one block per modulus).


def listing_fp_mul():
    """Montgomery product for a 4-limb modulus (machine.md section 2)."""
    L = []
    b = ["x3", "x4", "x5", "x6"]
    t = ["x7", "x8", "x9", "x10", "x11", "x12", "x13", "x14"]
    T = ["x15", "x16", "x17", "x2"]
    L += ["ldp x3, x4, [x2, #0]", "ldp x5, x6, [x2, #16]"]
    # row 0, a0 held in x14 (t7 is not live before row 3)
    L += ["ldr x14, [x1, #0]"]
    L += ["mul %s, x14, %s" % (t[j], b[j]) for j in range(4)]
    L += ["umulh %s, x14, %s" % (T[j], b[j]) for j in range(4)]
    L += ["adds x8, x8, x15", "adcs x9, x9, x16", "adcs x10, x10, x17", "adc x11, x2, xzr"]
    for i in (1, 2, 3):
        ai = "x14" if i < 3 else "x1"
        L += ["ldr %s, [x1, #%d]" % (ai, 8 * i)]
        L += ["mul %s, %s, %s" % (T[j], ai, b[j]) for j in range(4)]
        L += ["adds %s, %s, %s" % (t[i], t[i], T[0])]
        L += ["adcs %s, %s, %s" % (t[i + j], t[i + j], T[j]) for j in (1, 2, 3)]
        L += ["adc %s, xzr, xzr" % t[i + 4]]
        L += ["umulh %s, %s, %s" % (T[j], ai, b[j]) for j in range(4)]
        L += ["adds %s, %s, %s" % (t[i + 1], t[i + 1], T[0])]
        L += ["adcs %s, %s, %s" % (t[i + 1 + j], t[i + 1 + j], T[j]) for j in (1, 2)]
        L += ["adc %s, %s, %s" % (t[i + 4], t[i + 4], T[3])]
    # constants: p0..p3 -> x3..x6, n0 -> x1
    L += ["adr x2, K", "ldp x3, x4, [x2, #0]", "ldp x5, x6, [x2, #16]", "ldr x1, [x2, #32]"]
    pr = ["x3", "x4", "x5", "x6"]
    for i in range(4):
        m = t[i]
        L += ["subs xzr, %s, #1" % t[i]]              # C = (t_i != 0) = carry of t_i + lo(m*p0)
        L += ["mul %s, %s, x1" % (m, t[i])]           # m = t_i * n0 (flags untouched)
        L += ["mul %s, %s, %s" % (T[j], m, pr[j]) for j in (1, 2, 3)]
        L += ["adcs %s, %s, %s" % (t[i + j], t[i + j], T[j]) for j in (1, 2, 3)]
        L += ["umulh %s, %s, %s" % (T[j], m, pr[j]) for j in range(4)]
        L += ["adc x2, x2, xzr"]                      # hi3 + carry <= 2^64-1
        L += ["adds %s, %s, %s" % (t[i + 1], t[i + 1], T[0])]
        L += ["adcs %s, %s, %s" % (t[i + 1 + j], t[i + 1 + j], T[j]) for j in (1, 2, 3)]
        for k in range(i + 5, 8):
            L += ["adcs %s, %s, xzr" % (t[k], t[k])]
    U = t[4:]
    L += ["subs x15, %s, x3" % U[0], "sbcs x16, %s, x4" % U[1], "sbcs x17, %s, x5" % U[2],
          "sbcs x2, %s, x6" % U[3]]
    L += ["csel %s, %s, %s, cs" % (U[j], T[j], U[j]) for j in range(4)]
    L += ["stp %s, %s, [x0, #0]" % (U[0], U[1]), "stp %s, %s, [x0, #16]" % (U[2], U[3]), "ret"]
    return L


def listing_fp_add():
    """(a+b) mod m for a, b < m < 2^254 (a+b < 2^255, no carry out of 256 bits)."""
    return ["ldp x3, x4, [x1, #0]", "ldp x5, x6, [x1, #16]",
            "ldp x7, x8, [x2, #0]", "ldp x9, x10, [x2, #16]",
            "adds x3, x3, x7", "adcs x4, x4, x8", "adcs x5, x5, x9", "adc x6, x6, x10",
            "adr x11, K", "ldp x7, x8, [x11, #0]", "ldp x9, x10, [x11, #16]",
            "subs x7, x3, x7", "sbcs x8, x4, x8", "sbcs x9, x5, x9", "sbcs x10, x6, x10",
            "csel x3, x7, x3, cs", "csel x4, x8, x4, cs", "csel x5, x9, x5, cs", "csel x6, x10, x6, cs",
            "stp x3, x4, [x0, #0]", "stp x5, x6, [x0, #16]", "ret"]


def listing_fp_sub():
    """(a-b) mod m for a, b < m: subtract, then add m masked by the borrow."""
    return ["ldp x3, x4, [x1, #0]", "ldp x5, x6, [x1, #16]",
            "ldp x7, x8, [x2, #0]", "ldp x9, x10, [x2, #16]",
            "subs x3, x3, x7", "sbcs x4, x4, x8", "sbcs x5, x5, x9", "sbcs x6, x6, x10",
            "csetm x11, cc",
            "adr x12, K", "ldp x7, x8, [x12, #0]", "ldp x9, x10, [x12, #16]",
            "and x7, x7, x11", "and x8, x8, x11", "and x9, x9, x11", "and x10, x10, x11",
            "adds x3, x3, x7", "adcs x4, x4, x8", "adcs x5, x5, x9", "adc x6, x6, x10",
            "stp x3, x4, [x0, #0]", "stp x5, x6, [x0, #16]", "ret"]


CALLEE = ["x%d" % i for i in range(19, 29)]
ALLOWED = {"x%d" % i for i in range(18)} | set(CALLEE) | {"xzr", "sp"}


def simulate(listing, mem, regs, kaddr):
    """Execute a straight-line listing. Returns the executed instruction count."""
    C = 0
    count = 0

    def rd(x):
        if x not in ALLOWED:
            raise RuntimeError("register not allowed: " + x)
        return 0 if x == "xzr" else regs[x]

    def wr(x, v):
        if x not in ALLOWED:
            raise RuntimeError("register not allowed: " + x)
        if x != "xzr":
            regs[x] = v & M64

    def addr(s):
        s = s.strip("[]")
        bn, off = [u.strip() for u in s.split(",")]
        return rd(bn) + int(off.lstrip("#"))

    for line in listing:
        count += 1
        op, _, rest = line.partition(" ")
        args = [a.strip() for a in rest.replace(", [", "|[").split("|")]
        ops = [a.strip() for a in args[0].split(",")] if args[0] else []
        if op == "ret":
            break
        if op == "ldr":
            wr(ops[0], mem[addr(args[1])])
        elif op == "str":
            mem[addr(args[1])] = rd(ops[0])
        elif op == "ldp":
            a = addr(args[1])
            wr(ops[0], mem[a]); wr(ops[1], mem[a + 8])
        elif op == "stp":
            a = addr(args[1])
            mem[a] = rd(ops[0]); mem[a + 8] = rd(ops[1])
        elif op == "adr":
            wr(ops[0], kaddr)
        elif op == "mul":
            wr(ops[0], rd(ops[1]) * rd(ops[2]))
        elif op == "umulh":
            wr(ops[0], (rd(ops[1]) * rd(ops[2])) >> 64)
        elif op in ("adds", "adcs", "adc", "add"):
            cin = C if op in ("adcs", "adc") else 0
            bv = int(ops[2].lstrip("#")) if ops[2].startswith("#") else rd(ops[2])
            s = rd(ops[1]) + bv + cin
            wr(ops[0], s)
            if op in ("adds", "adcs"):
                C = s >> 64
        elif op in ("subs", "sbcs", "sub"):
            bv = int(ops[2].lstrip("#")) if ops[2].startswith("#") else rd(ops[2])
            bin_ = (1 - C) if op == "sbcs" else 0
            d = rd(ops[1]) - bv - bin_
            wr(ops[0], d)
            if op != "sub":
                C = 1 if d >= 0 else 0
        elif op == "csel":
            wr(ops[0], rd(ops[1]) if (C == 1) == (ops[3] == "cs") else rd(ops[2]))
        elif op == "csetm":
            wr(ops[0], M64 if (C == 1) == (ops[1] == "cs") else 0)
        elif op == "and":
            wr(ops[0], rd(ops[1]) & rd(ops[2]))
        else:
            raise RuntimeError("unsupported " + line)
    return count


def limbs(v):
    return [(v >> (64 * i)) & M64 for i in range(4)]


def run_kernel(listing, m, a, b):
    n0 = (-pow(m, -1, 1 << 64)) % (1 << 64)
    mem = {}
    A, B, O, K = 0x1000, 0x2000, 0x3000, 0x4000
    for i, w in enumerate(limbs(a)):
        mem[A + 8 * i] = w
    for i, w in enumerate(limbs(b)):
        mem[B + 8 * i] = w
    for i, w in enumerate(limbs(m) + [n0]):
        mem[K + 8 * i] = w
    regs = {"x%d" % i: rng.getrandbits(64) for i in range(18)}
    regs.update({"x0": O, "x1": A, "x2": B})
    n = simulate(listing, mem, regs, K)
    return sum(mem[O + 8 * i] << (64 * i) for i in range(4)), n


# ---- fused Fp2 leaves (lazy reduction inside one straight-line listing)
def mulu_core(A, B, t, T):
    """t[0..7] <- A*B (4x4 limbs, operand scanning by rows); T = 4 temporaries. 63 instructions."""
    L = ["mul %s, %s, %s" % (t[j], A[0], B[j]) for j in range(4)]
    L += ["umulh %s, %s, %s" % (T[j], A[0], B[j]) for j in range(4)]
    L += ["adds %s, %s, %s" % (t[1], t[1], T[0]), "adcs %s, %s, %s" % (t[2], t[2], T[1]),
          "adcs %s, %s, %s" % (t[3], t[3], T[2]), "adc %s, %s, xzr" % (t[4], T[3])]
    for i in (1, 2, 3):
        L += ["mul %s, %s, %s" % (T[j], A[i], B[j]) for j in range(4)]
        L += ["adds %s, %s, %s" % (t[i], t[i], T[0])]
        L += ["adcs %s, %s, %s" % (t[i + j], t[i + j], T[j]) for j in (1, 2, 3)]
        L += ["adc %s, xzr, xzr" % t[i + 4]]
        L += ["umulh %s, %s, %s" % (T[j], A[i], B[j]) for j in range(4)]
        L += ["adds %s, %s, %s" % (t[i + 1], t[i + 1], T[0])]
        L += ["adcs %s, %s, %s" % (t[i + 1 + j], t[i + 1 + j], T[j]) for j in (1, 2)]
        L += ["adc %s, %s, %s" % (t[i + 4], t[i + 4], T[3])]
    return L


def red_core(t, P, n0, T):
    """Montgomery reduction of an 8-limb t < m*2^256 to the canonical residue (82 instructions)."""
    L = []
    for i in range(4):
        m = t[i]
        L += ["subs xzr, %s, #1" % t[i]]              # C = (t_i != 0) = carry of t_i + lo(q*p0)
        L += ["mul %s, %s, %s" % (m, t[i], n0)]       # q = t_i * n0 (flags untouched)
        L += ["mul %s, %s, %s" % (T[j], m, P[j]) for j in (1, 2, 3)]
        L += ["adcs %s, %s, %s" % (t[i + j], t[i + j], T[j]) for j in (1, 2, 3)]
        L += ["umulh %s, %s, %s" % (T[j], m, P[j]) for j in range(4)]
        L += ["adc %s, %s, xzr" % (T[3], T[3])]       # hi3 + carry <= 2^64-1
        L += ["adds %s, %s, %s" % (t[i + 1], t[i + 1], T[0])]
        L += ["adcs %s, %s, %s" % (t[i + 1 + j], t[i + 1 + j], T[j]) for j in (1, 2, 3)]
        for k in range(i + 5, 8):
            L += ["adcs %s, %s, xzr" % (t[k], t[k])]
    U = t[4:]
    L += ["subs %s, %s, %s" % (T[0], U[0], P[0])]
    L += ["sbcs %s, %s, %s" % (T[j], U[j], P[j]) for j in (1, 2, 3)]
    L += ["csel %s, %s, %s, cs" % (U[j], T[j], U[j]) for j in range(4)]
    return L


def ld4(regs, base, off):
    return ["ldp %s, %s, [%s, #%d]" % (regs[0], regs[1], base, off),
            "ldp %s, %s, [%s, #%d]" % (regs[2], regs[3], base, off + 16)]


def st4(regs, base, off):
    return ["stp %s, %s, [%s, #%d]" % (regs[0], regs[1], base, off),
            "stp %s, %s, [%s, #%d]" % (regs[2], regs[3], base, off + 16)]


def listing_f2mul():
    """Fused Fp2 product with lazy Karatsuba (machine.md 2a): 3 unreduced products, 2 Montgomery reductions."""
    t = ["x11", "x12", "x13", "x14", "x15", "x16", "x17", "x19"]
    T = ["x20", "x21", "x22", "x23"]
    X = ["x3", "x4", "x5", "x6"]
    Y = ["x7", "x8", "x9", "x10"]
    Z = ["x24", "x25", "x26", "x27"]
    L = ["sub sp, sp, #208"]
    L += ["stp x%d, x%d, [sp, #%d]" % (19 + 2 * i, 20 + 2 * i, 16 * i) for i in range(5)]
    L += ld4(X, "x1", 0) + ld4(Y, "x2", 0)                      # a0, b0
    L += mulu_core(X, Y, t, T)
    L += st4(t[:4], "sp", 80) + st4(t[4:], "sp", 112)           # T0
    L += ld4(X, "x1", 32) + ld4(Y, "x2", 32)                    # a1, b1
    L += mulu_core(X, Y, t, T)
    L += st4(t[:4], "sp", 144) + st4(t[4:], "sp", 176)          # T1
    L += ld4(Z, "x1", 0)                                        # a0 again
    L += ["adds %s, %s, %s" % (X[0], X[0], Z[0])] + ["adcs %s, %s, %s" % (X[j], X[j], Z[j]) for j in (1, 2)]
    L += ["adc %s, %s, %s" % (X[3], X[3], Z[3])]               # s = a0 + a1 < 2p
    L += ld4(Z, "x2", 0)                                        # b0 again
    L += ["adds %s, %s, %s" % (Y[0], Y[0], Z[0])] + ["adcs %s, %s, %s" % (Y[j], Y[j], Z[j]) for j in (1, 2)]
    L += ["adc %s, %s, %s" % (Y[3], Y[3], Z[3])]               # t = b0 + b1 < 2p
    L += mulu_core(X, Y, t, T)                                  # T2
    U = ["x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10"]       # T0
    V = ["x24", "x25", "x26", "x27", "x28", "x1", "x2", "x20"]  # T1 (pointers a, b are dead)
    L += ld4(U[:4], "sp", 80) + ld4(U[4:], "sp", 112)
    L += ld4(V[:4], "sp", 144) + ld4(V[4:], "sp", 176)
    L += ["subs %s, %s, %s" % (t[0], t[0], U[0])] + ["sbcs %s, %s, %s" % (t[j], t[j], U[j]) for j in range(1, 8)]
    L += ["subs %s, %s, %s" % (t[0], t[0], V[0])] + ["sbcs %s, %s, %s" % (t[j], t[j], V[j]) for j in range(1, 8)]
    L += ["subs %s, %s, %s" % (U[0], U[0], V[0])] + ["sbcs %s, %s, %s" % (U[j], U[j], V[j]) for j in range(1, 8)]
    L += ["csetm x21, cc"]                                      # mask = all ones iff T0 < T1
    Pm = ["x24", "x25", "x26", "x27"]
    L += ["adr x28, K", "ldp x24, x25, [x28, #0]", "ldp x26, x27, [x28, #16]", "ldr x28, [x28, #32]"]
    Tm = ["x1", "x2", "x20", "x22"]
    L += ["and %s, %s, x21" % (Tm[j], Pm[j]) for j in range(4)]
    L += ["adds %s, %s, %s" % (U[4], U[4], Tm[0])] + ["adcs %s, %s, %s" % (U[4 + j], U[4 + j], Tm[j]) for j in (1, 2)]
    L += ["adc %s, %s, %s" % (U[7], U[7], Tm[3])]
    L += red_core(t, Pm, "x28", Tm)
    L += st4(t[4:], "x0", 32)                                   # c1
    L += red_core(U, Pm, "x28", Tm)
    L += st4(U[4:], "x0", 0)                                    # c0
    L += ["ldp x%d, x%d, [sp, #%d]" % (19 + 2 * i, 20 + 2 * i, 16 * i) for i in range(5)]
    L += ["add sp, sp, #208", "ret"]
    return L


def listing_f2sqr():
    """Fused Fp2 squaring (machine.md 2a): c1 = REDC(2 a0 a1), c0 = REDC((a0+a1)(a0-a1+p))."""
    t = ["x11", "x12", "x13", "x14", "x15", "x16", "x17", "x19"]
    T = ["x20", "x21", "x22", "x2"]
    X = ["x3", "x4", "x5", "x6"]
    Y = ["x7", "x8", "x9", "x10"]
    Pm = ["x23", "x24", "x25", "x26"]
    D = ["x27", "x28", "x1", "x20"]
    L = ["sub sp, sp, #80"]
    L += ["stp x%d, x%d, [sp, #%d]" % (19 + 2 * i, 20 + 2 * i, 16 * i) for i in range(5)]
    L += ld4(X, "x1", 0) + ld4(Y, "x1", 32)                     # a0, a1
    L += ["adr x2, K", "ldp x23, x24, [x2, #0]", "ldp x25, x26, [x2, #16]", "ldr x27, [x2, #32]"]
    L += mulu_core(X, Y, t, T)
    L += ["adds %s, %s, %s" % (t[0], t[0], t[0])] + ["adcs %s, %s, %s" % (t[j], t[j], t[j]) for j in range(1, 7)]
    L += ["adc %s, %s, %s" % (t[7], t[7], t[7])]
    L += red_core(t, Pm, "x27", T)
    L += st4(t[4:], "x0", 32)                                   # c1
    L += ["add x21, x27, xzr"]                                  # keep n0
    L += ["subs %s, %s, %s" % (D[0], X[0], Y[0])] + ["sbcs %s, %s, %s" % (D[j], X[j], Y[j]) for j in (1, 2, 3)]
    L += ["adds %s, %s, %s" % (D[0], D[0], Pm[0])] + ["adcs %s, %s, %s" % (D[j], D[j], Pm[j]) for j in (1, 2)]
    L += ["adc %s, %s, %s" % (D[3], D[3], Pm[3])]              # d = a0 - a1 + p (mod 2^256; true value in (0,2p))
    L += ["adds %s, %s, %s" % (X[0], X[0], Y[0])] + ["adcs %s, %s, %s" % (X[j], X[j], Y[j]) for j in (1, 2)]
    L += ["adc %s, %s, %s" % (X[3], X[3], Y[3])]               # s = a0 + a1
    T2 = ["x7", "x8", "x9", "x10"]
    L += mulu_core(X, D, t, T2)
    L += red_core(t, Pm, "x21", T2)
    L += st4(t[4:], "x0", 0)                                    # c0
    L += ["ldp x%d, x%d, [sp, #%d]" % (19 + 2 * i, 20 + 2 * i, 16 * i) for i in range(5)]
    L += ["add sp, sp, #80", "ret"]
    return L


def listing_f2mulxi():
    """(9 + u)(a0 + a1 u), canonical in and out (machine.md 2a)."""
    W = ["x3", "x4", "x5", "x6"]
    T = ["x7", "x8", "x9", "x10"]
    Pm = ["x11", "x12", "x13", "x14"]
    A = ["x15", "x16", "x17", "x2"]

    def addred(src):                # W <- (W + src) mod p, W, src < p
        return (["adds %s, %s, %s" % (W[0], W[0], src[0])]
                + ["adcs %s, %s, %s" % (W[j], W[j], src[j]) for j in (1, 2)]
                + ["adc %s, %s, %s" % (W[3], W[3], src[3])]
                + ["subs %s, %s, %s" % (T[0], W[0], Pm[0])]
                + ["sbcs %s, %s, %s" % (T[j], W[j], Pm[j]) for j in (1, 2, 3)]
                + ["csel %s, %s, %s, cs" % (W[j], T[j], W[j]) for j in range(4)])

    def nine(off):                  # A <- a_off, W <- 9 a_off mod p
        L = ld4(A, "x1", off)
        L += ["add %s, %s, xzr" % (W[j], A[j]) for j in range(4)]
        for _ in range(3):
            L += addred(W)
        return L + addred(A)

    L = ["sub sp, sp, #32", "adr x14, K", "ldp x11, x12, [x14, #0]", "ldp x13, x14, [x14, #16]"]
    L += nine(0)                                                # W = 9 a0
    L += ld4(T, "x1", 32)                                       # a1
    L += ["subs %s, %s, %s" % (W[0], W[0], T[0])] + ["sbcs %s, %s, %s" % (W[j], W[j], T[j]) for j in (1, 2, 3)]
    L += ["csetm x15, cc"]
    L += ["and %s, %s, x15" % (T[j], Pm[j]) for j in range(4)]
    L += ["adds %s, %s, %s" % (W[0], W[0], T[0])] + ["adcs %s, %s, %s" % (W[j], W[j], T[j]) for j in (1, 2)]
    L += ["adc %s, %s, %s" % (W[3], W[3], T[3])]
    L += ["stp x3, x4, [sp, #0]", "stp x5, x6, [sp, #16]"]     # c0 parked in the frame (out may alias a)
    L += nine(32)                                               # W = 9 a1
    L += ld4(A, "x1", 0)
    L += addred(A)                                              # + a0
    L += st4(W, "x0", 32)
    L += ["ldp x3, x4, [sp, #0]", "ldp x5, x6, [sp, #16]"]
    L += st4(W, "x0", 0)
    L += ["add sp, sp, #32", "ret"]
    return L


# ---- fused tag kernel c1dot (new in this package)
def mul_rows(Abase, Aoff, a, B, t, T):
    """t <- A*B, multiplier limb loaded per row into register a (mulu rows + 4 loads)."""
    L = ["ldr %s, [%s, #%d]" % (a, Abase, Aoff)]
    L += ["mul %s, %s, %s" % (t[j], a, B[j]) for j in range(4)]
    L += ["umulh %s, %s, %s" % (T[j], a, B[j]) for j in range(4)]
    L += ["adds %s, %s, %s" % (t[1], t[1], T[0]), "adcs %s, %s, %s" % (t[2], t[2], T[1]),
          "adcs %s, %s, %s" % (t[3], t[3], T[2]), "adc %s, %s, xzr" % (t[4], T[3])]
    for i in (1, 2, 3):
        L += ["ldr %s, [%s, #%d]" % (a, Abase, Aoff + 8 * i)]
        L += ["mul %s, %s, %s" % (T[j], a, B[j]) for j in range(4)]
        L += ["adds %s, %s, %s" % (t[i], t[i], T[0])]
        L += ["adcs %s, %s, %s" % (t[i + j], t[i + j], T[j]) for j in (1, 2, 3)]
        L += ["adc %s, xzr, xzr" % t[i + 4]]
        L += ["umulh %s, %s, %s" % (T[j], a, B[j]) for j in range(4)]
        L += ["adds %s, %s, %s" % (t[i + 1], t[i + 1], T[0])]
        L += ["adcs %s, %s, %s" % (t[i + 1 + j], t[i + 1 + j], T[j]) for j in (1, 2)]
        L += ["adc %s, %s, %s" % (t[i + 4], t[i + 4], T[3])]
    return L

def redc_rounds(t, P, n0, T):
    """four Montgomery word rounds (74 instructions) of red_core: t[4..7] <- (t + q p)/2^256"""
    L = []
    for i in range(4):
        m = t[i]
        L += ["subs xzr, %s, #1" % t[i]]
        L += ["mul %s, %s, %s" % (m, t[i], n0)]
        L += ["mul %s, %s, %s" % (T[j], m, P[j]) for j in (1, 2, 3)]
        L += ["adcs %s, %s, %s" % (t[i + j], t[i + j], T[j]) for j in (1, 2, 3)]
        L += ["umulh %s, %s, %s" % (T[j], m, P[j]) for j in range(4)]
        L += ["adc %s, %s, xzr" % (T[3], T[3])]
        L += ["adds %s, %s, %s" % (t[i + 1], t[i + 1], T[0])]
        L += ["adcs %s, %s, %s" % (t[i + 1 + j], t[i + 1 + j], T[j]) for j in (1, 2, 3)]
        for k in range(i + 5, 8):
            L += ["adcs %s, %s, xzr" % (t[k], t[k])]
    return L

def condsub(U, M, T):
    """U <- U - M if U >= M (8 instructions)"""
    return (["subs %s, %s, %s" % (T[0], U[0], M[0])] + ["sbcs %s, %s, %s" % (T[j], U[j], M[j]) for j in (1, 2, 3)]
            + ["csel %s, %s, %s, cs" % (U[j], T[j], U[j]) for j in range(4)])

# constant block K: p (0..24), n0 (32), 2p (40..64), 6p^2 (72..128, 8 limbs)
def kblock(m):
    n0 = (-pow(m, -1, 1 << 64)) % (1 << 64)
    w = [(m >> (64 * i)) & M64 for i in range(4)] + [n0] + [((2 * m) >> (64 * i)) & M64 for i in range(4)]
    w += [((6 * m * m) >> (64 * i)) & M64 for i in range(8)]
    return w

def listing_c1dot():
    """(Re, Im) of sum_{i<6} x_i Qt_i * 2^-256 mod p, canonical (machine.md 2b). x1 base block, x2 tag block (stride 96: e0, e1, e0+e1), x0 output."""
    acc = ["x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10", "x11"]
    t = ["x12", "x13", "x14", "x15", "x16", "x17", "x19", "x20"]
    T = ["x21", "x22", "x23", "x24"]
    B = ["x25", "x26", "x27", "x28"]
    a = "x0"
    L = ["sub sp, sp, #256"]
    L += ["stp x%d, x%d, [sp, #%d]" % (19 + 2 * i, 20 + 2 * i, 16 * i) for i in range(5)]
    L += ["str x0, [sp, #80]"]
    for part, off in enumerate((0, 32, 64)):
        for i in range(6):
            L += ["ldp %s, %s, [x2, #%d]" % (B[0], B[1], 96 * i + off), "ldp %s, %s, [x2, #%d]" % (B[2], B[3], 96 * i + off + 16)]
            if i == 0:
                L += mul_rows("x1", 96 * i + off, a, B, acc[:8], T)
                L += ["add %s, xzr, xzr" % acc[8]]
            else:
                L += mul_rows("x1", 96 * i + off, a, B, t, T)
                L += ["adds %s, %s, %s" % (acc[0], acc[0], t[0])] + ["adcs %s, %s, %s" % (acc[j], acc[j], t[j]) for j in range(1, 8)]
                L += ["adc %s, %s, xzr" % (acc[8], acc[8])]
        if part < 2:     # park S0 at sp+96, S1 at sp+176 (9 limbs each)
            base = 96 + 80 * part
            L += ["stp %s, %s, [sp, #%d]" % (acc[2 * j], acc[2 * j + 1], base + 16 * j) for j in range(4)]
            L += ["str %s, [sp, #%d]" % (acc[8], base + 64)]
    t9 = t + [T[0]]
    def load9(dst, base):
        return ["ldp %s, %s, [sp, #%d]" % (dst[2 * j], dst[2 * j + 1], base + 16 * j) for j in range(4)] + ["ldr %s, [sp, #%d]" % (dst[8], base + 64)]
    def sub9(dst, src):
        return ["subs %s, %s, %s" % (dst[0], dst[0], src[0])] + ["sbcs %s, %s, %s" % (dst[j], dst[j], src[j]) for j in range(1, 9)]
    L += load9(t9, 96) + sub9(acc, t9) + load9(t9, 176) + sub9(acc, t9)       # acc = S2 - S0 - S1
    P_, M2 = B, ["x13", "x14", "x15", "x16"]
    n0 = "x12"
    def reduce_out(src, outoff):
        R = ["adr x17, K", "ldp %s, %s, [x17, #0]" % (P_[0], P_[1]), "ldp %s, %s, [x17, #16]" % (P_[2], P_[3]),
             "ldr %s, [x17, #32]" % n0, "ldp %s, %s, [x17, #40]" % (M2[0], M2[1]), "ldp %s, %s, [x17, #56]" % (M2[2], M2[3])]
        R += redc_rounds(src, P_, n0, T)
        U = src[4:8]
        R += condsub(U, M2, T) + condsub(U, P_, T)
        R += ["ldr x0, [sp, #80]", "stp %s, %s, [x0, #%d]" % (U[0], U[1], outoff), "stp %s, %s, [x0, #%d]" % (U[2], U[3], outoff + 16)]
        return R
    L += reduce_out(acc[:8], 32)                                               # b
    L += load9(acc, 96) + load9(t9, 176) + sub9(acc, t9)                       # acc = S0 - S1
    C8 = t                                                                     # 6p^2
    L += ["adr x22, K"] + ["ldp %s, %s, [x22, #%d]" % (C8[2 * j], C8[2 * j + 1], 72 + 16 * j) for j in range(4)]
    L += ["adds %s, %s, %s" % (acc[0], acc[0], C8[0])] + ["adcs %s, %s, %s" % (acc[j], acc[j], C8[j]) for j in range(1, 8)]
    L += reduce_out(acc[:8], 0)                                                # a
    L += ["ldp x%d, x%d, [sp, #%d]" % (19 + 2 * i, 20 + 2 * i, 16 * i) for i in range(5)]
    L += ["add sp, sp, #256", "ret"]
    return L

def listing_add4sums():
    """Unreduced sums e0 + e1 (< 2p) at offset 96 i + 64 of the block at x0, i < 6 (machine.md 2b)."""
    L = []
    for i in range(6):
        o = 96 * i
        L += ["ldp x1, x2, [x0, #%d]" % o, "ldp x3, x4, [x0, #%d]" % (o + 16),
              "ldp x5, x6, [x0, #%d]" % (o + 32), "ldp x7, x8, [x0, #%d]" % (o + 48),
              "adds x1, x1, x5", "adcs x2, x2, x6", "adcs x3, x3, x7", "adc x4, x4, x8",
              "stp x1, x2, [x0, #%d]" % (o + 64), "stp x3, x4, [x0, #%d]" % (o + 80)]
    return L + ["ret"]


def run_add4sums(L, xs):
    mem = {}
    A = 0x10000
    for i, e in enumerate(xs):
        for j, v in enumerate(e):
            for k in range(4):
                mem[A + 96 * i + 32 * j + 8 * k] = (v >> (64 * k)) & M64
    regs = {"x%d" % i: rng_f.getrandbits(64) for i in range(29)}
    regs.update({"x0": A, "sp": 0x90000})
    n = simulate(L, mem, regs, 0)
    return [sum(mem[A + 96 * i + 64 + 8 * k] << (64 * k) for k in range(4)) for i in range(6)], n


def run_c1dot(L, xs, qs, xsums=None):
    """xs, qs: six Fp2 pairs each (stored Montgomery residues, canonical)."""
    mem = {}
    A, Bb, O, K, SP = 0x10000, 0x20000, 0x30000, 0x40000, 0x90000
    def put(base, off, v):
        for i in range(4): mem[base + off + 8 * i] = (v >> (64 * i)) & M64
    for i in range(6):
        for base, e in ((A, xs[i]), (Bb, qs[i])):
            put(base, 96 * i, e[0]); put(base, 96 * i + 32, e[1]); put(base, 96 * i + 64, e[0] + e[1])
        if xsums is not None:
            put(A, 96 * i + 64, xsums[i])
    for i, w in enumerate(kblock(p)): mem[K + 8 * i] = w
    regs = {"x%d" % i: rng_f.getrandbits(64) for i in range(29)}
    regs.update({"x0": O, "x1": A, "x2": Bb, "sp": SP})
    saved = {k: regs[k] for k in CALLEE + ["sp"]}
    before = dict(mem)
    n = simulate(L, mem, regs, K)
    assert all(regs[k] == saved[k] for k in saved)
    for k_ in mem:
        if k_ not in before or mem[k_] != before[k_]:
            assert O <= k_ < O + 64 or SP - 256 <= k_ < SP, hex(k_)
    get = lambda off: sum(mem[O + off + 8 * i] << (64 * i) for i in range(4))
    return (get(0), get(32)), n


def run_fused(listing, m, ins, nout, alias=False):
    """ins: list of field elements for operand a (x1) and optionally b (x2) as (alist, blist)."""
    n0 = (-pow(m, -1, 1 << 64)) % (1 << 64)
    mem = {}
    A, B, O, K, SP = 0x1000, 0x2000, 0x3000, 0x4000, 0x9000
    if alias:
        O = A
    for base_, vals in ((A, ins[0]), (B, ins[1])):
        for k, v in enumerate(vals):
            for i, w in enumerate(limbs(v)):
                mem[base_ + 32 * k + 8 * i] = w
    for i, w in enumerate(limbs(m) + [n0]):
        mem[K + 8 * i] = w
    regs = {"x%d" % i: rng_f.getrandbits(64) for i in range(29)}
    regs.update({"x0": O, "x1": A, "x2": B, "sp": SP})
    saved = {k: regs[k] for k in CALLEE + ["sp"]}
    n = simulate(listing, mem, regs, K)
    assert all(regs[k] == saved[k] for k in saved), "callee-saved register or sp not restored"
    assert all(SP - 256 <= a < SP or a < 0x5000 for a in mem), "write outside operands, constants and the frame"
    return [sum(mem[O + 32 * k + 8 * i] << (64 * i) for i in range(4)) for k in range(nout)], n



def part_K(cases):
    out = {}
    Ls = {"fp_mul": listing_fp_mul(), "fp_add": listing_fp_add(), "fp_sub": listing_fp_sub()}
    counts = set()
    for name, mod in (("p", p), ("r", r)):
        Rinv = pow(RMONT, -1, mod)
        edge = [0, 1, 2, mod - 2, mod - 1, (1 << 64) - 1, 1 << 64, mod >> 1, (mod >> 1) + 1,
                (1 << 192) - 1, 1 << 192, mod - (1 << 64)]
        pairs = [(x, y) for x in edge for y in edge]
        pairs += [(rng.randrange(mod), rng.randrange(mod)) for _ in range(cases)]
        # inputs whose limbs are all large exercise every carry chain
        pairs += [(mod - 1 - rng.getrandbits(32), mod - 1 - rng.getrandbits(32)) for _ in range(cases // 4)]
        fails = 0
        for x, y in pairs:
            want = {"fp_mul": x * y * Rinv % mod, "fp_add": (x + y) % mod, "fp_sub": (x - y) % mod}
            for k, L in Ls.items():
                got, n = run_kernel(L, mod, x, y)
                counts.add((k, n))
                fails += got != want[k]
        out["modulus_" + name] = {"cases": len(pairs), "failures": fails}
    # fused Fp2 leaves, modulus p only (they are used for field arithmetic, not for exponents)
    Fs = {"f2mul": listing_f2mul(), "f2sqr": listing_f2sqr(), "f2mulxi": listing_f2mulxi()}
    Rinv = pow(RMONT, -1, p)
    small = [0, 1, p - 1, (1 << 192) - 1]
    quads = [(a0, a1, b0, b1) for a0 in small for a1 in small for b0 in (0, p - 1) for b1 in (0, p - 1, p >> 1)]
    quads += [tuple(rng_f.randrange(p) for _ in range(4)) for _ in range(max(1, min(cases, 300) // 6))]
    quads += [tuple(p - 1 - rng_f.getrandbits(32) for _ in range(4)) for _ in range(max(1, min(cases, 300) // 12))]
    ffails = 0
    for a0, a1, b0, b1 in quads:
        want = {"f2mul": [(a0 * b0 - a1 * b1) * Rinv % p, (a0 * b1 + a1 * b0) * Rinv % p],
                "f2sqr": [(a0 * a0 - a1 * a1) * Rinv % p, 2 * a0 * a1 * Rinv % p],
                "f2mulxi": [(9 * a0 - a1) % p, (9 * a1 + a0) % p]}
        for k, L in Fs.items():
            for alias in (False, True):
                got, n = run_fused(L, p, ([a0, a1], [b0, b1]), 2, alias)
                counts.add((k, n))
                ffails += got != want[k]
    out["fused_leaves_modulus_p"] = {"cases": 2 * len(quads), "failures": ffails,
                                     "callee_saved_and_sp_restored": True}
    Ls.update(Fs)
    # fused tag kernel: (Re, Im) of sum_{i<6} x_i Qt_i * 2^-256 mod p
    L1 = listing_c1dot()
    Ls["c1dot"] = L1
    edge = [0, 1, p - 1, p - 2, (1 << 192) - 1, p >> 1, 1 << 64]
    tests = []
    for e in (0, p - 1):
        tests += [([(e, e)] * 6, [(p - 1, p - 1)] * 6), ([(e, 0)] * 6, [(0, p - 1)] * 6), ([(0, e)] * 6, [(p - 1, 0)] * 6)]
    tests += [([(rng_f.randrange(p), rng_f.randrange(p)) for _ in range(6)],
               [(rng_f.randrange(p), rng_f.randrange(p)) for _ in range(6)]) for _ in range(max(1, min(cases, 300) // 2))]
    tests += [([(rng_f.choice(edge), rng_f.choice(edge)) for _ in range(6)],
               [(rng_f.choice(edge), rng_f.choice(edge)) for _ in range(6)]) for _ in range(max(1, min(cases, 300) // 4))]
    cfails = 0
    for xs, qs in tests:
        want = (sum(x[0] * q[0] - x[1] * q[1] for x, q in zip(xs, qs)) * Rinv % p,
                sum(x[0] * q[1] + x[1] * q[0] for x, q in zip(xs, qs)) * Rinv % p)
        got, n = run_c1dot(L1, xs, qs)
        counts.add(("c1dot", n))
        cfails += got != want
    L4 = listing_add4sums()
    Ls["add4sums"] = L4
    afails = 0
    for xs, _ in tests:
        got, n = run_add4sums(L4, xs)
        counts.add(("add4sums", n))
        afails += got != [e[0] + e[1] for e in xs]
    out["unreduced_sums_add4sums"] = {"cases": len(tests), "failures": afails}
    out["tag_kernel_c1dot_modulus_p"] = {"cases": len(tests), "failures": cfails,
                                         "callee_saved_and_sp_restored": True,
                                         "writes_only_output_and_frame": True}
    out["executed_instructions"] = {k: n for k, n in sorted(counts)}
    out["listing_lengths"] = {k: len(L) for k, L in Ls.items()}
    import hashlib
    txt = "\n".join(L1) + "\n"
    out["c1dot_listing_sha256"] = hashlib.sha256(txt.encode()).hexdigest()
    lf = base / "c1dot.listing.txt"
    out["c1dot_listing_file_matches"] = lf.read_text() == txt if lf.is_file() else "file not present"
    if SHOW_LISTINGS:
        out["listings"] = {k: v for k, v in Ls.items() if k != "c1dot"}
    return out, out["listing_lengths"]


# =============================================================== counted tower
class Tower:
    """F_q2 = F_q[u]/(u^2+1), F_q12 = F_q2[w]/(w^6 - xi), xi = 9+u; optional operation counting."""

    def __init__(self, q):
        self.q = q
        self.count = None
        self.xi = (9 % q, 1)
        self.d = [self.f2pow(self.xi, m * (q - 1) // 6) for m in range(6)]
        z = self.f2pow(self.xi, (q * q - 1) // 6)
        assert z[1] == 0 and (z[0] * z[0] - z[0] + 1) % q == 0
        self.zeta = z[0]
        self.one = ((1, 0),) + ((0, 0),) * 5
        self._fk = {}

    # ---- base field, counted
    def c(self, k, n=1):
        if self.count is not None:
            self.count[k] = self.count.get(k, 0) + n

    def fadd(self, a, b):
        self.c("fp_add")
        return (a + b) % self.q

    def fsub(self, a, b):
        self.c("fp_sub")
        return (a - b) % self.q

    def fmul(self, a, b):
        self.c("fp_mul")
        return a * b % self.q

    # ---- Fp2, counted (add/sub/xi are inline pairs of leaf calls; mul/sqr are routines)
    def f2add(self, a, b):
        return (self.fadd(a[0], b[0]), self.fadd(a[1], b[1]))

    def f2sub(self, a, b):
        return (self.fsub(a[0], b[0]), self.fsub(a[1], b[1]))

    def f2mulxi(self, a):
        # fused leaf: (a0 + a1 u)(9 + u) = (9a0 - a1) + (9a1 + a0) u
        self.c("f2mulxi")
        q = self.q
        return ((9 * a[0] - a[1]) % q, (9 * a[1] + a[0]) % q)

    def f2mul(self, a, b):
        # fused leaf: lazy Karatsuba, two Montgomery reductions
        self.c("f2mul")
        q = self.q
        return ((a[0] * b[0] - a[1] * b[1]) % q, (a[0] * b[1] + a[1] * b[0]) % q)

    def f2sqr(self, a):
        # fused leaf: complex squaring, two Montgomery reductions
        self.c("f2sqr")
        q = self.q
        return ((a[0] + a[1]) * (a[0] - a[1]) % q, 2 * a[0] * a[1] % q)

    def f2conj(self, a):
        return (a[0], self.fsub(0, a[1]))

    # ---- uncounted helpers
    def f2mul_raw(self, a, b):
        q = self.q
        return ((a[0] * b[0] - a[1] * b[1]) % q, (a[0] * b[1] + a[1] * b[0]) % q)

    def f2pow(self, a, e):
        res = (1, 0)
        while e:
            if e & 1:
                res = self.f2mul_raw(res, a)
            a = self.f2mul_raw(a, a)
            e >>= 1
        return res

    # ---- Fp6 / Fp12 multiplication (Karatsuba), counted
    def f6mul(self, a, b, bs=None):
        """bs: the three sums (b1+b2, b0+b1, b0+b2) when b is a table operand (stored, not recomputed)."""
        self.c("routine")
        a0, a1, a2 = a
        b0, b1, b2 = b
        if bs is None:
            bs = (self.f2add(b1, b2), self.f2add(b0, b1), self.f2add(b0, b2))
        v0, v1, v2 = self.f2mul(a0, b0), self.f2mul(a1, b1), self.f2mul(a2, b2)
        t = self.f2mul(self.f2add(a1, a2), bs[0])
        c0 = self.f2add(v0, self.f2mulxi(self.f2sub(self.f2sub(t, v1), v2)))
        t = self.f2mul(self.f2add(a0, a1), bs[1])
        c1 = self.f2add(self.f2sub(self.f2sub(t, v0), v1), self.f2mulxi(v2))
        t = self.f2mul(self.f2add(a0, a2), bs[2])
        c2 = self.f2sub(self.f2add(self.f2sub(t, v0), v1), v2)
        return (c0, c1, c2)

    def presums(self, y):
        """Table-side Karatsuba sums of a multiplier (768 bytes per entry), computed when the table is built."""
        q = self.q
        add2 = lambda a, b: ((a[0] + b[0]) % q, (a[1] + b[1]) % q)
        s3 = lambda b: (add2(b[1], b[2]), add2(b[0], b[1]), add2(b[0], b[2]))
        Y0, Y1 = (y[0], y[2], y[4]), (y[1], y[3], y[5])
        Ys = tuple(add2(u, v) for u, v in zip(Y0, Y1))
        return (Ys, s3(Y0), s3(Y1), s3(Ys))

    def f6add(self, a, b):
        return tuple(self.f2add(x, y) for x, y in zip(a, b))

    def f6sub(self, a, b):
        return tuple(self.f2sub(x, y) for x, y in zip(a, b))

    def mul_counted(self, x, y, pre=None):
        """pre = presums(y) for a table operand; None recomputes the sums (counted)."""
        self.c("routine")
        X0, X1 = (x[0], x[2], x[4]), (x[1], x[3], x[5])
        Y0, Y1 = (y[0], y[2], y[4]), (y[1], y[3], y[5])
        if pre is None:
            Ys, s0, s1, ss = self.f6add(Y0, Y1), None, None, None
        else:
            Ys, s0, s1, ss = pre
        t0, t1 = self.f6mul(X0, Y0, s0), self.f6mul(X1, Y1, s1)
        vt1 = (self.f2mulxi(t1[2]), t1[0], t1[1])
        C0 = self.f6add(t0, vt1)
        t2 = self.f6mul(self.f6add(X0, X1), Ys, ss)
        C1 = self.f6sub(self.f6sub(t2, t0), t1)
        return (C0[0], C1[0], C0[1], C1[1], C0[2], C1[2])

    def mul(self, x, y):
        """uncounted schoolbook product, used for speed in setup and toy runs"""
        q = self.q
        acc = [[0, 0] for _ in range(11)]
        for i in range(6):
            a0, a1 = x[i]
            for j in range(6):
                b0, b1 = y[j]
                t = acc[i + j]
                t[0] += a0 * b0 - a1 * b1
                t[1] += a0 * b1 + a1 * b0
        out = []
        for m in range(6):
            c0, c1 = acc[m]
            if m + 6 <= 10:
                h0, h1 = acc[m + 6]            # w^6 = 9 + u
                c0 += 9 * h0 - h1
                c1 += 9 * h1 + h0
            out.append((c0 % q, c1 % q))
        return tuple(out)

    def pow(self, x, e):
        res = self.one
        while e:
            if e & 1:
                res = self.mul(res, x)
            x = self.mul(x, x)
            e >>= 1
        return res

    def frobk(self, y, k):
        k %= 12
        if k == 0:
            return y
        D = self._fk.get(k)
        if D is None:
            D = [self.f2pow(self.xi, m * (self.q ** k - 1) // 6) for m in range(6)]
            self._fk[k] = D
        q = self.q
        if k & 1:
            return tuple(self.f2mul_raw((c[0], (-c[1]) % q), D[m]) for m, c in enumerate(y))
        return tuple(self.f2mul_raw(c, D[m]) for m, c in enumerate(y))

    # ---- keys and cyclotomic squaring, counted
    def keys(self, y):
        """Re(w^1 coefficient of y^(q^k)), k = 0..11: 4 multiplications, 1 addition, 8 sub/neg."""
        a, b1 = y[1]
        B = self.fadd(self.fmul(a, self.d[1][0]), self.fmul(b1, self.d[1][1]))
        za, zb = self.fmul(a, self.zeta), self.fmul(B, self.zeta)
        A6 = [a, za, self.fsub(za, a), self.fsub(0, a), self.fsub(0, za), self.fsub(a, za)]
        B6 = [B, zb, self.fsub(zb, B), self.fsub(0, B), self.fsub(0, zb), self.fsub(B, zb)]
        out = [0] * 12
        for j in range(6):
            out[2 * j], out[2 * j + 1] = A6[j], B6[j]
        return out

    def keys_direct(self, y):
        return [self.frobk(y, k)[1][0] for k in range(12)]

    def fp4sqr(self, a, b):
        self.c("routine")
        a2, b2 = self.f2sqr(a), self.f2sqr(b)
        ab2 = self.f2sub(self.f2sub(self.f2sqr(self.f2add(a, b)), a2), b2)
        return self.f2add(a2, self.f2mulxi(b2)), ab2

    def cyc_sqr(self, y):
        """Granger-Scott squaring on the cyclotomic subgroup (proof.md F3)."""
        self.c("routine")
        c0, c1, c2, c3, c4, c5 = y
        A0, A1 = self.fp4sqr(c0, c3)
        B0, B1 = self.fp4sqr(c1, c4)
        C0, C1 = self.fp4sqr(c2, c5)
        xC1 = self.f2mulxi(C1)

        def minus(t, u):        # 3t - 2u = t + 2(t - u)
            dd = self.f2sub(t, u)
            return self.f2add(t, self.f2add(dd, dd))

        def plus(t, u):         # 3t + 2u = t + 2(t + u)
            ss = self.f2add(t, u)
            return self.f2add(t, self.f2add(ss, ss))

        return (minus(A0, c0), plus(xC1, c1), minus(B0, c2), plus(A1, c3), minus(C0, c4), plus(B1, c5))

    # ---- tag tracing (new): c_1 of a product from a stored tag block
    def tagvec(self, P_):
        """Qt with c_1(x * P_) = sum_i x_i Qt_i: (P1, P0, xi P5, xi P4, xi P3, xi P2)."""
        X = lambda c_: self.f2mul_raw(self.xi, c_)
        return (P_[1], P_[0], X(P_[5]), X(P_[4]), X(P_[3]), X(P_[2]))

    def c1_from_tag(self, x, Qt):
        q = self.q
        re = sum(a[0] * b[0] - a[1] * b[1] for a, b in zip(x, Qt)) % q
        im = sum(a[0] * b[1] + a[1] * b[0] for a, b in zip(x, Qt)) % q
        return (re, im)

    def keys_c1(self, c1):
        """the twelve keys from the w^1 coefficient alone (same counted operations as keys)"""
        return self.keys((None, c1))

    def rand_gt(self, rng_, rp):
        hc = (self.q ** 12 - 1) // rp
        while True:
            x = tuple((rng_.randrange(self.q), rng_.randrange(self.q)) for _ in range(6))
            y = self.pow(x, hc)
            if y != self.one:
                return y


def profile_to_tower(c, q):
    """Profile basis to tower coefficients (proof.md section 2)."""
    return tuple(((c[m] + 9 * c[m + 6]) % q, c[m + 6] % q) for m in range(6))


def poly_mul_profile(a, b, q):
    acc = [0] * 23
    for i in range(12):
        for j in range(12):
            acc[i + j] += a[i] * b[j]
    for k in range(22, 11, -1):
        h = acc[k]
        acc[k - 6] += 18 * h
        acc[k - 12] -= 82 * h
        acc[k] = 0
    return [v % q for v in acc[:12]]


def mult_order(a, n, fact):
    o = n - 1
    for pf, e in fact.items():
        for _ in range(e):
            if pow(a, o // pf, n) == 1:
                o //= pf
            else:
                break
    return o


# =============================================================== Part A
def part_A(checks):
    P = Tower(p)
    out = {}
    prod = 1
    for k, e in R_FACT.items():
        prod *= k ** e
    o2 = mult_order(2, r, R_FACT)
    op_ = mult_order(p % r, r, R_FACT)
    out.update({
        "p_mod_12": p % 12,
        "r_minus_1_factorisation_verified": prod == r - 1,
        "order_of_p_mod_r": op_,
        "p6_is_minus1_mod_r": pow(p, 6, r) == r - 1,
        "phi12_p_is_0_mod_r": (pow(p, 4, r) - pow(p, 2, r) + 1) % r == 0,
        "log2_order_of_2_mod_r": round(math.log2(o2), 3),
        "smallest_t_with_2^t_in_<p>_log2": round(math.log2(o2 // math.gcd(o2, 12)), 3),
        "zeta_in_Fp_primitive_6th_root": (P.zeta * P.zeta - P.zeta + 1) % p == 0,
        "xi_not_square_not_cube_in_Fp2": P.f2pow(P.xi, (p * p - 1) // 2) != (1, 0) and P.f2pow(P.xi, (p * p - 1) // 3) != (1, 0),
    })
    # profile basis -> tower: ring homomorphism check on random elements, and g of order r
    a = [rng.randrange(p) for _ in range(12)]
    b = [rng.randrange(p) for _ in range(12)]
    out["profile_to_tower_is_ring_map"] = profile_to_tower(poly_mul_profile(a, b, p), p) == P.mul(profile_to_tower(a, p), profile_to_tower(b, p))
    g = profile_to_tower(G_PROFILE, p)
    out["g_order_r"] = P.pow(g, r) == P.one and g != P.one
    ok = dict.fromkeys(["frob", "keys", "distinct", "canon", "inv", "sq", "mul", "trace_inv", "stored_canon", "minkey_inv", "tag_c1"], True)
    counts_mul, counts_sq, counts_keys = None, None, None
    Rp = RMONT % p
    for _ in range(checks):
        y = P.pow(g, rng.randrange(1, r))
        x2 = P.pow(g, rng.randrange(1, r))
        ok["frob"] &= P.frobk(y, 1) == P.pow(y, p)
        kd, kf = P.keys_direct(y), P.keys(y)
        ok["keys"] &= kd == kf
        ok["distinct"] &= len(set(kf)) == 12
        conj = [P.frobk(y, k) for k in range(12)]
        canon, canon_st = set(), set()
        for cj in conj:
            kk = P.keys(cj)
            canon.add(P.frobk(cj, kk.index(min(kk))))
            ks = [v * Rp % p for v in kk]          # stored Montgomery limbs
            canon_st.add(P.frobk(cj, ks.index(min(ks))))
        ok["canon"] &= len(canon) == 1
        ok["stored_canon"] &= len(canon_st) == 1
        ok["inv"] &= P.mul(conj[6], y) == P.one
        ok["trace_inv"] &= len({cj[0][0] for cj in conj}) == 1
        ok["minkey_inv"] &= len({min(v * Rp % p for v in P.keys(cj)) for cj in conj}) == 1
        ok["tag_c1"] &= P.c1_from_tag(y, P.tagvec(x2)) == P.mul(y, x2)[1]
        P.count = {}
        s = P.cyc_sqr(y)
        counts_sq = dict(P.count)
        P.count = {}
        m = P.mul_counted(y, x2, P.presums(x2))
        counts_mul = dict(P.count)
        P.count = None
        ok["mul"] &= P.mul_counted(y, x2) == m
        P.count = {}
        P.keys(y)
        counts_keys = dict(P.count)
        P.count = None
        ok["sq"] &= s == P.mul(y, y)
        ok["mul"] &= m == P.mul(y, x2)
    out.update({"random_GT_elements": checks, **{k + "_ok": v for k, v in ok.items()}})
    out["op_counts_fp12_mul"] = counts_mul
    out["op_counts_cyclotomic_square"] = counts_sq
    out["op_counts_keys"] = counts_keys
    return out, counts_mul, counts_sq, counts_keys


# =============================================================== Part C: production walk
def part_C(steps, rprod=32, chunk=CHUNK_L, kernel_checks=24):
    """Production walk in the real BN254 GT, tag-traced and direct in lock step (proof.md 3.2, 3.4, 7): compares every c_1, K, decision, materialisation and escape; runs c1dot on real stored values."""
    P = Tower(p)
    Rp = RMONT % p
    st = lambda v: v * Rp % p
    lam = [pow(p, k, r) for k in range(12)]
    g = profile_to_tower(G_PROFILE, p)
    x0 = rng.randrange(1, r)
    z0 = P.pow(g, x0)
    ah = [rng.randrange(r) for _ in range(rprod)]
    bh = [rng.randrange(r) for _ in range(rprod)]
    Mh = [P.mul(P.pow(g, ah[i]), P.pow(z0, bh[i])) for i in range(rprod)]
    item = {}
    for h in range(rprod):
        for t in range(12):
            item[12 * h + t] = (P.frobk(Mh[h], t), ah[h] * lam[t] % r, bh[h] * lam[t] % r)
    cache = {(): (P.one, 0, 0)}

    def entry(ms):                      # multiset product table entry (built on demand here)
        e = cache.get(ms)
        if e is None:
            par = entry(ms[:-1])
            it = item[ms[-1]]
            e = (P.mul(par[0], it[0]), (par[1] + it[1]) % r, (par[2] + it[2]) % r)
            cache[ms] = e
        return e

    L1 = listing_c1dot()
    L4 = listing_add4sums()
    stats = dict(steps=0, c1_mismatch=0, K_mismatch=0, decision_mismatch=0, materialised=0,
                 materialise_mismatch=0, two_cycle_escapes=0, checkpoint_escapes=0, ties=0,
                 kernel_runs=0, kernel_failures=0)

    def info(c1):
        ks = [st(v) for v in P.keys_c1(c1)]
        K = min(ks)
        stats["ties"] += ks.count(K) > 1
        return K, ks.index(K)

    class TagWalk:
        def __init__(s, y, a, b):
            s.y0, s.a0, s.b0, s.ms = y, a, b, ()
        def c1(s):
            if not s.ms:
                return s.y0[1]
            Qt = P.tagvec(entry(s.ms)[0])
            c1 = P.c1_from_tag(s.y0, Qt)
            if stats["kernel_runs"] < kernel_checks:
                xs = [(st(a), st(b)) for a, b in s.y0]
                sums, _ = run_add4sums(L4, xs)          # the charged routine builds the base sums
                got, _ = run_c1dot(L1, xs, [(st(a), st(b)) for a, b in Qt], sums)
                stats["kernel_runs"] += 1
                stats["kernel_failures"] += got != (st(c1[0]), st(c1[1]))
            return c1
        def handle(s):                   # what is needed to materialise the current element later
            return (s.y0, s.a0, s.b0, s.ms)
        @staticmethod
        def element(hd):
            y0, a0, b0, ms = hd
            E_ = entry(ms)
            return P.mul(y0, E_[0]), (a0 + E_[1]) % r, (b0 + E_[2]) % r
        def move(s, h, t):
            s.ms = tuple(sorted(s.ms + (12 * h + t,)))
            if len(s.ms) == chunk:
                E_ = entry(s.ms)
                s.y0 = P.mul(s.y0, E_[0])
                s.a0, s.b0 = (s.a0 + E_[1]) % r, (s.b0 + E_[2]) % r
                s.ms = ()
                return True
            return False
        def restart(s, y, a, b):
            s.__init__(y, a, b)

    class DirectWalk:
        def __init__(s, y, a, b):
            s.y, s.a, s.b = y, a, b
        def c1(s):
            return s.y[1]
        def handle(s):
            return (s.y, s.a, s.b)
        @staticmethod
        def element(hd):
            return hd
        def move(s, h, t):
            it = item[12 * h + t]
            s.y, s.a, s.b = P.mul(s.y, it[0]), (s.a + it[1]) % r, (s.b + it[2]) % r
            return False
        def restart(s, y, a, b):
            s.__init__(y, a, b)

    a_, b_ = rng.randrange(r), rng.randrange(r)
    y_ = P.mul(P.pow(g, a_), P.pow(z0, b_))
    walks = [TagWalk(y_, a_, b_), DirectWalk(y_, a_, b_)]
    ctl = [dict(hist=[], ck=None, n=0) for _ in walks]
    bins3, bins4 = [0] * 8, [0] * 32
    t0 = time.time()

    def escape(w, c, hd):            # square the element behind handle hd; new chunk base
        y, a, b = w.element(hd)
        w.restart(P.cyc_sqr(y), 2 * a % r, 2 * b % r)
        c.update(hist=[], ck=None, n=0)

    while stats["steps"] < steps and time.time() - t0 < PROD_SECONDS:
        stats["steps"] += 1
        c1s = [w.c1() for w in walks]
        stats["c1_mismatch"] += c1s[0] != c1s[1]
        infos = [info(c1) for c1 in c1s]
        stats["K_mismatch"] += infos[0][0] != infos[1][0]
        K = infos[0][0]
        bins3[K & 7] += 1
        bins4[(K >> 3) & 31] += 1
        events = []
        for w, c, (Kw, kap) in zip(walks, ctl, infos):
            hd = w.handle()
            if len(c["hist"]) >= 2 and c["hist"][-2][0] == Kw:
                _, hb = min(c["hist"][-1], (Kw, hd), key=lambda e_: e_[0])
                escape(w, c, hb)
                events.append("two")
                continue
            if c["ck"] is not None and c["ck"] == Kw:
                # cycle of length <= 64: traverse it once more, keep the least-K handle, square it
                best = (Kw, hd)
                for _ in range(64):
                    kk = info(w.c1())
                    w.move((kk[0] >> 3) % rprod, (12 - kk[1]) % 12)
                    k2 = info(w.c1())[0]
                    if k2 == Kw:
                        break
                    if k2 < best[0]:
                        best = (k2, w.handle())
                escape(w, c, best[1])
                events.append("ck")
                continue
            if c["n"] % 64 == 0:
                c["ck"] = Kw
            c["hist"] = (c["hist"] + [(Kw, hd)])[-2:]
            c["n"] += 1
            w.move((Kw >> 3) % rprod, (12 - kap) % 12)
            events.append("step")
        stats["decision_mismatch"] += len(set(events)) != 1
        if events[0] == "two":
            stats["two_cycle_escapes"] += 1
        if events[0] == "ck":
            stats["checkpoint_escapes"] += 1
        if events[0] != "step" or not walks[0].ms:
            stats["materialised"] += 1
            stats["materialise_mismatch"] += walks[0].element(walks[0].handle()) != walks[1].handle()
    yt = TagWalk.element(walks[0].handle())
    exp_ok = yt == walks[1].handle() and yt[0] == P.mul(P.pow(g, yt[1]), P.pow(z0, yt[2]))
    n = sum(bins3)
    chi3 = sum((c_ - n / 8) ** 2 / (n / 8) for c_ in bins3) if n else None
    chi4 = sum((c_ - n / 32) ** 2 / (n / 32) for c_ in bins4) if n else None
    return {"R_table": rprod, "chunk_length": chunk, **stats,
            "multiset_entries_built": len(cache) - 1,
            "exponent_invariant_and_final_element_agree": exp_ok,
            "chi2_low3_bits_of_stored_K_7dof": round(chi3, 2) if chi3 is not None else None,
            "chi2_index_bits_3_to_7_of_stored_K_31dof": round(chi4, 2) if chi4 is not None else None}


# =============================================================== Part B: toy algorithm
def toy_prime(q):
    n = q ** 4 - q * q + 1
    d, best = 2, 1
    while d * d <= n:
        while n % d == 0:
            best = max(best, d)
            n //= d
        d += 1
    return max(best, n)


def mix64(x):
    x &= M64
    x = ((x ^ (x >> 30)) * 0xBF58476D1CE4E5B9) & M64
    x = ((x ^ (x >> 27)) * 0x94D049BB133111EB) & M64
    return x ^ (x >> 31)


def classinfo(T, y):
    kf = T.keys(y)
    m = min(kf)
    ks = [k for k in range(12) if kf[k] == m]
    if len(ks) == 1:
        return ks[0], T.frobk(y, ks[0]), m
    can, kap = min((T.frobk(y, k), k) for k in ks)
    return kap, can, m


def toy_hash(can):
    h = 0x9E3779B97F4A7C15
    for c_ in can:
        h = mix64(h ^ (c_[0] * 0x100000001B3 + c_[1]))
    return h


def run_toy(T, q, rp, walks, Rm, sq8, dp_bits, cap_factor=40):
    qk = [pow(q, k, rp) for k in range(12)]
    g = T.rand_gt(rng, rp)
    xs = rng.randrange(1, rp)
    z0 = T.pow(g, xs)
    ah = [rng.randrange(rp) for _ in range(Rm)]
    bh = [rng.randrange(rp) for _ in range(Rm)]
    Mh = [T.mul(T.pow(g, ah[i]), T.pow(z0, bh[i])) for i in range(Rm)]
    Tab = [[T.frobk(Mh[i], k) for k in range(12)] for i in range(Rm)]
    At = [[ah[i] * qk[k] % rp for k in range(12)] for i in range(Rm)]
    Bt = [[bh[i] * qk[k] % rp for k in range(12)] for i in range(Rm)]
    mask = (1 << dp_bits) - 1

    def new_walk():
        a, b = rng.randrange(rp), rng.randrange(rp)
        return {"y": T.mul(T.pow(g, a), T.pow(z0, b)), "a": a, "b": b, "hist": [], "ck": None, "n": 0, "since": 0}

    def advance(y, a, b):
        kap, can, K = classinfo(T, y)
        u = toy_hash(can)
        if (u & 7) < sq8:
            return T.cyc_sqr(y), 2 * a % rp, 2 * b % rp
        h = (u >> 3) % Rm
        idx = (12 - kap) % 12
        return T.mul(y, Tab[h][idx]), (a + At[h][idx]) % rp, (b + Bt[h][idx]) % rp

    W = [new_walk() for _ in range(walks)]
    visited, dps = {}, {}
    steps, Tc, found = 0, None, None
    st = {"two_cycles": 0, "long_cycles": 0, "restarts": 0, "useless": 0, "sq": 0, "ord": 0}
    while found is None:
        for w in W:
            y, a, b = w["y"], w["a"], w["b"]
            kap, can, K = classinfo(T, y)
            u = toy_hash(can)
            ca, cb = a * qk[kap] % rp, b * qk[kap] % rp
            if Tc is None:
                old = visited.get(can)
                if old is None:
                    visited[can] = (ca, cb)
                elif old != (ca, cb):
                    Tc = steps
            if (u >> 13) & mask == 0:
                old = dps.get(can)
                if old is not None and old[1] != cb:
                    xr = (old[0] - ca) * pow(cb - old[1], -1, rp) % rp
                    if T.pow(g, xr) != z0:
                        raise RuntimeError("recovered logarithm failed verification")
                    found = xr
                    break
                if old is not None:
                    st["useless"] += 1
                dps[can] = (ca, cb)
                w["since"] = 0
            hist = w["hist"]
            key = (K, can)
            if len(hist) >= 2 and hist[-2][0][1] == can:
                st["two_cycles"] += 1
                _, yy, aa, bb = min(hist[-1], (key, y, a, b), key=lambda s_: s_[0])
                w["y"], w["a"], w["b"] = T.cyc_sqr(yy), 2 * aa % rp, 2 * bb % rp
                w["hist"], w["ck"] = [], None
                steps += 1
                continue
            if w["ck"] is not None and w["ck"][0] == can and w["ck"][1] < w["n"]:
                st["long_cycles"] += 1
                best = (key, y, a, b)
                yy, aa, bb = y, a, b
                for _ in range(w["n"] - w["ck"][1]):
                    yy, aa, bb = advance(yy, aa, bb)
                    steps += 1
                    k2, c2, K2 = classinfo(T, yy)
                    if (K2, c2) < best[0]:
                        best = ((K2, c2), yy, aa, bb)
                _, yy, aa, bb = best
                w["y"], w["a"], w["b"] = T.cyc_sqr(yy), 2 * aa % rp, 2 * bb % rp
                w["hist"], w["ck"] = [], None
                steps += 1
                continue
            if w["n"] % 64 == 0:
                w["ck"] = (can, w["n"])
            hist.append((key, y, a, b))
            if len(hist) > 2:
                hist.pop(0)
            st["ord"] += 1
            st["sq"] += (u & 7) < sq8
            w["y"], w["a"], w["b"] = advance(y, a, b)
            w["n"] += 1
            w["since"] += 1
            steps += 1
            if w["since"] > cap_factor << dp_bits:
                st["restarts"] += 1
                w.update(new_walk())
    return Tc, steps, st


def part_B(q, trials, Rm, sq8, dp_bits, walks, budget):
    T = Tower(q)
    rp = toy_prime(q)
    if not (pow(q, 12, rp) == 1 and all(pow(q, k, rp) != 1 for k in range(1, 12))):
        raise RuntimeError("toy field does not have embedding degree 12")
    N = (rp - 1) // 12
    s = sq8 / 8
    sumq2 = s * s + (1 - s) ** 2 / Rm
    Np = N / (1 - sumq2)
    Tcs, ev = [], {"two_cycles": 0, "long_cycles": 0, "restarts": 0, "useless": 0}
    t0 = time.time()
    allsteps = 0
    while len(Tcs) < trials and time.time() - t0 < budget:
        Tc, steps, st = run_toy(T, q, rp, walks, Rm, sq8, dp_bits)
        Tcs.append(Tc)
        allsteps += steps
        for k in ev:
            ev[k] += st[k]
    n = len(Tcs)
    sN = math.sqrt(N)
    S = {x: math.sqrt(2 * Np * math.log(1 / (1 - x))) for x in (0.40, 0.42, 0.43)}
    return {"q": q, "r_prime": rp, "classes_N": N, "R": Rm, "squaring_prob": s,
            "effective_N_over_N": round(1 / (1 - sumq2), 4), "trials_completed": n,
            "logs_recovered_and_verified": n,
            "mean_Tc_over_sqrtN": round(sum(Tcs) / n / sN, 4) if n else None,
            "model_mean_over_sqrtN": round(math.sqrt(math.pi * Np / 2) / sN, 4),
            **{"fraction_Tc_le_S_model_%.2f" % x: round(sum(t <= S[x] for t in Tcs) / max(1, n), 3) for x in S},
            "events": ev, "two_cycles_per_step": round(ev["two_cycles"] / max(1, allsteps), 6),
            "predicted_two_cycles_per_step": round((1 - s) ** 2 / (12 * Rm), 6)}


# =============================================================== Part F: useless-cycle law
def part_F(mmax=16, lmax=300, R=None):
    """Exact useless-cycle law G(l) of the R-adding class walk under H2 (proof.md 3.4): q(m) exact for m <= mmax, bound 1/(3 m^2) above; G(l) = l!/R^l [x^l] (sum q(m) x^m/m!)^R."""
    R = R or R_TABLE
    pw = [(1, 0, 0, 0)]
    for _ in range(11):
        a, b, c_, d = pw[-1]
        pw.append((-d, a, b + d, c_))
    dist = {(0, 0, 0, 0): 1}
    q = [Fraction(1)]
    for m in range(1, mmax + 1):
        nd = {}
        for v, c_ in dist.items():
            for s_ in pw:
                w = (v[0] + s_[0], v[1] + s_[1], v[2] + s_[2], v[3] + s_[3])
                nd[w] = nd.get(w, 0) + c_
        dist = nd
        q.append(Fraction(dist.get((0, 0, 0, 0), 0), 12 ** m))
    m2q_ok = all(m * m * q[m] <= Fraction(1, 3) for m in range(1, mmax + 1))
    qb = [float(q[m]) if m <= mmax else 1 / (3 * m * m) for m in range(lmax + 1)]
    NEG = -math.inf                                # log-domain coefficients (no under/overflow)
    f = [math.log(qb[m]) - math.lgamma(m + 1) if qb[m] > 0 else NEG for m in range(lmax + 1)]

    def sq(a):
        out = []
        for l in range(lmax + 1):
            t = [a[i] + a[l - i] for i in range(l + 1) if a[i] > NEG and a[l - i] > NEG]
            M = max(t) if t else NEG
            out.append(M + math.log(sum(math.exp(x - M) for x in t)) if t else NEG)
        return out
    Pw = f
    for _ in range(R.bit_length() - 1):          # R = 32 = 2^5: five squarings
        Pw = sq(Pw)
    G = [math.exp(Pw[l] + math.lgamma(l + 1) - l * math.log(R)) if Pw[l] > NEG else 0.0 for l in range(lmax + 1)]
    lG = [(Pw[l] + math.lgamma(l + 1) - l * math.log(R)) / math.log(2) if Pw[l] > NEG else NEG for l in range(lmax + 1)]
    esc = None
    two = G[2] * (2 + 4)                                     # 2-cycle: two steps plus the escape (< 4 steps)
    mid = sum(G[l] * (67 + 2 * l) for l in range(3, 65))     # checkpoint delay <= 63 + l, traverse l, escape
    Mt = max(lG[65:])
    tail_log2 = Mt + math.log2(sum(2 ** (x - Mt) for x in lG[65:]))
    n12 = 12 * R
    charged = 8 / n12 + 1000 / n12 ** 2
    return {"R": R, "q_m_exact": {str(m): str(q[m]) for m in range(1, mmax + 1)},
            "m2_q_m_at_most_one_third": m2q_ok,
            "G_2": G[2], "G_2_equals_1_over_12R": abs(G[2] - 1 / n12) < 1e-15,
            "log2_G": {str(l): round(lG[l], 3) for l in (3, 4, 5, 6, 8, 16, 32, 64, 100, 200, 300)},
            "sum_G_3_to_64": sum(G[3:65]), "log2_sum_G_65_to_%d" % lmax: round(tail_log2, 2),
            "G_nonincreasing_from_l6": all(lG[l + 1] <= lG[l] for l in range(6, lmax)),
            "escape_overhead_two_cycles": two, "escape_overhead_cycles_3_to_64": mid,
            "escape_overhead_total": two + mid, "eps_charged": charged, "eps_covers": two + mid <= charged}


# =============================================================== Part D: ledger
CALL_SITE = 4          # BL plus up to three pointer set-ups (ADD/MOV) per leaf call
ROUTINE = 32           # per non-leaf routine call: call site, frame, 4 callee-saved, epilogue, slack
COMMON = 300           # every step (itemised in machine.md section 4: 125)
KEYMIN = 160           # minimum of 12 four-limb keys, kappa, twist index (itemised 138, tie path included)
THETA_BITS = 70         # distinguished iff (K >> 13) mod 2^70 = 0
CHUNK_BOOK = 200       # per chunk: base swap, sums of the new base, multiset reset (itemised in machine.md)
R_TABLE = 32
CHUNK = CHUNK_L
P_MODEL = Decimal("0.43")
CLAIM = Decimal("123.02")


def ledger(lens, cmul, csq, ckeys):
    getcontext().prec = 60
    ln2 = Decimal(2).ln()
    lg = lambda v: Decimal(v).ln() / ln2
    leaf = {k: lens[k] + CALL_SITE for k in ("fp_mul", "fp_add", "fp_sub", "f2mul", "f2sqr", "f2mulxi", "c1dot")}
    leaf["routine"] = ROUTINE

    def price(cnt):
        return sum(leaf[k] * v for k, v in cnt.items())

    MULI, SQI = price(cmul), price(csq)
    KEYI = price(ckeys) + KEYMIN
    TAGI = leaf["c1dot"]
    # per chunk of CHUNK steps: CHUNK key computations, CHUNK - 1 tags (step 0 reads c_1 of the base),
    # one dense product, 6 fp_add for the new base's sums, 2 fp_add (mod r) for the exponents, bookkeeping
    SUMS = lens["add4sums"] + CALL_SITE          # unreduced base sums for c1dot (listing add4sums)
    CHUNKI = MULI + SUMS + 2 * leaf["fp_add"] + CHUNK_BOOK
    per_chunk = CHUNK * (KEYI + COMMON) + (CHUNK - 1) * TAGI + CHUNKI
    C = 20000
    cbar = Fraction(per_chunk, CHUNK * C)
    cb = Decimal(cbar.numerator) / Decimal(cbar.denominator)
    N = Fraction(r - 1, 12)
    sumq2 = Fraction(1, R_TABLE)
    Np = N / (1 - sumq2)
    Npd = Decimal(Np.numerator) / Decimal(Np.denominator)
    S = (2 * Npd * (1 / (1 - P_MODEL)).ln()).sqrt()
    n12 = 12 * R_TABLE
    eps = Decimal(8) / Decimal(n12) + Decimal(1000) / Decimal(n12 * n12)
    E = Decimal(2 ** 20) * 20 * Decimal(2 ** THETA_BITS)
    NONWALK = Decimal(2) ** Decimal("57")
    W = cb * S * (1 + eps) + cb * E + NONWALK
    lw = lg(W)
    ceil2 = (lw * 100).to_integral_value(rounding=ROUND_CEILING) / 100
    # escape price in instructions: materialise (dense product), cyclotomic square, exponents, new base sums
    ESC = MULI + SQI + 2 * leaf["fp_add"] + SUMS + COMMON
    from math import comb
    n_items = 12 * R_TABLE
    entries = [comb(n_items + k - 1, k) for k in range(1, CHUNK + 1)]
    table_bytes = sum(entries[:-1]) * (384 + 768 + 576 + 64) + entries[-1] * (384 + 768 + 64)
    build_units = Fraction(sum(entries) * (MULI + 24 * leaf["fp_add"] + 4 * leaf["f2mulxi"] + 6 * leaf["fp_add"] + 400), C)
    Wcap = Decimal(2) ** CLAIM
    return {"leaf_prices_incl_call_site": leaf,
            "fp12_mul_instr": MULI, "cyclotomic_square_instr": SQI,
            "keys_instr": KEYI, "tag_c1dot_instr": TAGI, "common_instr": COMMON,
            "chunk_instr": CHUNKI, "chunk_length": CHUNK, "per_chunk_instr": per_chunk,
            "R_table": R_TABLE, "items_12R": n_items,
            "cbar_units_per_step": float(cbar), "instr_per_step": float(Fraction(per_chunk, CHUNK)),
            "log2_N": float(lg(Decimal(N.numerator) / N.denominator)),
            "sum_q2": float(sumq2), "log2_N_eff": float(lg(Npd)), "log2_S": float(lg(S)), "eps": float(eps),
            "escape_instr": ESC, "escape_in_steps": float(Fraction(ESC * CHUNK, per_chunk)),
            "log2_E": float(lg(E)), "log2_nonwalk_bucket": 57, "theta_bits": THETA_BITS,
            "log2_reports_cap": float(lg(Decimal("1.1") * S / Decimal(2 ** THETA_BITS) + 2 ** 25)),
            "multiset_entries_by_size": entries, "log2_table_entries": float(lg(sum(entries))),
            "log2_table_bytes": float(lg(table_bytes)), "log2_table_build_units": float(lg(Decimal(build_units.numerator) / build_units.denominator)),
            "log2_W_nominal": float(lw), "least_hundredth_at_or_above": str(ceil2),
            "claim": str(CLAIM), "claim_minus_nominal_bits": float(CLAIM - lw),
            "claim_ok": CLAIM >= lw,
            "cap_slack_units_per_step": float((Wcap - cb * E - NONWALK) / (S * (1 + eps)) - cb)}


# =============================================================== main
out = {"seed": SEED}
kres, lens = part_K(KERNEL_CASES)
out["partK_kernels"] = kres
ares, cmul, csq, ckeys = part_A(max(1, cfg.get("production_checks", 2)))
out["partA_identities"] = ares
out["partD_ledger"] = ledger(lens, cmul, csq, ckeys)
out["partC_production_tag_walk"] = part_C(max(1, PROD_STEPS))
out["partF_useless_cycles"] = part_F()
# declared runs (seconds <= 14) stay inside the 20 s runner limit; longer budgets are offline only
budget_B = BUDGET if BUDGET > 14 else max(0.5, min(BUDGET, 15 - (time.time() - T_START)))
out["partB_toy"] = part_B(Q_TOY, TRIALS, R_TOY, SQ8, DP_BITS, WALKS, budget_B)
out["elapsed_seconds"] = round(time.time() - T_START, 2)
out["scope"] = ("Kernel listings are executed instruction by instruction on a simulator of the AArch64 "
                "subset they use; tower identities and the step rule run on the protected BN254 GT; "
                "collision statistics are toy-field evidence only. No production DLP is run.")
print(json.dumps(out, indent=1, default=str))
