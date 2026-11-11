"""AArch64 listings of the window re-base, an instruction-level simulator and an
encodability check.  Standard library only; loaded by experiment.py.

The four leaf listings fp_add, fp_sub, f2mul and f2mulxi (and their generators
mulu_core, red_core, ld4, st4) are copied unchanged from e-man07's submission
7e58d165-ebf6-41aa-8003-2f0a42370b6a (experiment.py, part K); fp_add and fp_sub
are, as stated there, the base-field leaves of 2ece9f04.  fr_add is fp_add
reading the constants of r (label KR) instead of p (label K).  The driver
`rebase` is new: Y0 <- Y0*T_e by the Karatsuba tower product with the table
operand's sums stored in the entry (Tower.mul_counted of 7e58d165), plus the two
label additions mod r, as straight-line calls of the leaves.
The simulator (class Machine) follows the one of 7e58d165, extended with bl,
ret, the link register x30 and per-function register checks.

Calling convention: x0 = output, x1, x2 = operands; leaves write only x0..x17,
the fused leaves save and restore the x19..x28 they use.  `rebase` keeps its
three pointers in x19..x21 and saves x19..x21 and the link register x30.
"""

M64 = (1 << 64) - 1


# ------------------------------------------------------------------ leaves (7e58d165)
def listing_fp_add(label="K"):
    """(a+b) mod m for a, b < m < 2^254 (a+b < 2^255, no carry out of 256 bits)."""
    return ["ldp x3, x4, [x1, #0]", "ldp x5, x6, [x1, #16]",
            "ldp x7, x8, [x2, #0]", "ldp x9, x10, [x2, #16]",
            "adds x3, x3, x7", "adcs x4, x4, x8", "adcs x5, x5, x9", "adc x6, x6, x10",
            "adr x11, %s" % label, "ldp x7, x8, [x11, #0]", "ldp x9, x10, [x11, #16]",
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
    """Montgomery reduction of the 8-limb value t < m*2^256: afterwards t[4..7] = t*2^-256 mod m,
    fully reduced. Four word rounds (74 instructions) and one conditional subtraction (8)."""
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
    """(c0 + c1 u) = (a0 + a1 u)(b0 + b1 u), u^2 = -1, canonical Montgomery inputs and outputs.
    Lazy Karatsuba: T0 = a0 b0, T1 = a1 b1, T2 = (a0+a1)(b0+b1) as unreduced 8-limb integers
    (T0, T1 < p^2, T2 < 4p^2 < p*2^256); c1 = REDC(T2 - T0 - T1), exact and in [0, 2p^2);
    c0 = REDC(T0 - T1 + [borrow] p*2^256), in [0, p*2^256). Frame: 80 bytes of callee-saved
    registers, 128 bytes for T0 and T1."""
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


def listing_f2mulxi():
    """(c0 + c1 u) = (9 + u)(a0 + a1 u) = (9 a0 - a1) + (9 a1 + a0) u, all values canonical.
    Each doubling/addition is followed by a conditional subtraction of p in registers."""
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


# ------------------------------------------------------------------ the re-base driver (new)
# Layouts (bytes).  Fp2 = 64 (c0, c1); Fp12 block = 384: the Fp2 coefficient a_i of w^i sits at
# 64*(3*(i%2) + i//2), i.e. Fp6 X0 = (a0, a2, a4) then X1 = (a1, a3, a5) -- the order of
# reference.py's BN254GT.coords and of the 48 limbs read by tt_run.
# Walk state: y at 0 (384), a at 384, b at 416 (mod r).  Table entry: F at 0 (768), next at 768
# (64), T at 832 (384), alpha at 1216, beta at 1248, then the stored Karatsuba sums at 1280:
# Ys = Y0+Y1 (3 Fp2), then for b in (Y0, Y1, Ys) the sums (b1+b2, b0+b1, b0+b2) (3 x 3 Fp2).
ENTRY_T, ENTRY_ALPHA, ENTRY_BETA, ENTRY_PRE, ENTRY_BYTES = 832, 1216, 1248, 1280, 2048
STATE_A, STATE_B = 384, 416
FRAME = 992


def coef_off(i):
    return 64 * (3 * (i % 2) + i // 2)


def listing_rebase():
    """rebase(x0 = out state, x1 = in state, x2 = entry): out.y = in.y * T_e (dense Fp12 product,
    Karatsuba over Fp6 and Fp2 as in 7e58d165, table-side sums loaded from the entry),
    out.a = in.a + alpha_e, out.b = in.b + beta_e (mod r).  out may equal in: every read of
    in.y precedes the first write to out."""
    L = ["sub sp, sp, #%d" % FRAME, "stp x19, x20, [sp, #0]", "stp x21, x30, [sp, #16]",
         "add x19, x0, #0", "add x20, x1, #0", "add x21, x2, #0"]
    OUT, IN, ENT, SP = "x19", "x20", "x21", "sp"
    fr = iter(range(32, FRAME, 64))
    tmp = {}

    def slot(name, n=1):            # n consecutive Fp2 slots in the frame
        if name not in tmp:
            tmp[name] = [(SP, next(fr)) for _ in range(n)]
        return tmp[name]

    def at(loc, d):
        return (loc[0], loc[1] + d)

    def call(fn, o, a, b=None):
        L.append("add x0, %s, #%d" % o)
        L.append("add x1, %s, #%d" % a)
        if b is not None:
            L.append("add x2, %s, #%d" % b)
        L.append("bl " + fn)

    def f2add(o, a, b):
        call("fp_add", o, a, b)
        call("fp_add", at(o, 32), at(a, 32), at(b, 32))

    def f2sub(o, a, b):
        call("fp_sub", o, a, b)
        call("fp_sub", at(o, 32), at(a, 32), at(b, 32))

    def f2mul(o, a, b):
        call("f2mul", o, a, b)

    def f2xi(o, a):
        call("f2mulxi", o, a)

    def f6mul(c, a, b, bs):
        """c = a*b in Fp6 (Karatsuba, 6 f2mul); bs = stored (b1+b2, b0+b1, b0+b2).
        c must not overlap a or the temporaries."""
        v = slot("v", 3)
        s, t, u = slot("s")[0], slot("t")[0], slot("u")[0]
        for k in range(3):
            f2mul(v[k], a[k], b[k])
        f2add(s, a[1], a[2]); f2mul(t, s, bs[0])
        f2sub(t, t, v[1]); f2sub(t, t, v[2]); f2xi(t, t); f2add(c[0], v[0], t)
        f2add(s, a[0], a[1]); f2mul(t, s, bs[1])
        f2sub(t, t, v[0]); f2sub(t, t, v[1]); f2xi(u, v[2]); f2add(c[1], t, u)
        f2add(s, a[0], a[2]); f2mul(t, s, bs[2])
        f2sub(t, t, v[0]); f2add(t, t, v[1]); f2sub(c[2], t, v[2])

    X0 = [(IN, coef_off(i)) for i in (0, 2, 4)]
    X1 = [(IN, coef_off(i)) for i in (1, 3, 5)]
    Y0 = [(ENT, ENTRY_T + coef_off(i)) for i in (0, 2, 4)]
    Y1 = [(ENT, ENTRY_T + coef_off(i)) for i in (1, 3, 5)]
    pre = [(ENT, ENTRY_PRE + 64 * k) for k in range(12)]
    Ys, s0, s1, ss = pre[0:3], pre[3:6], pre[6:9], pre[9:12]
    C0 = [(OUT, coef_off(i)) for i in (0, 2, 4)]
    C1 = [(OUT, coef_off(i)) for i in (1, 3, 5)]
    t0, t1, S = slot("t0", 3), slot("t1", 3), slot("S", 3)
    f6mul(t0, X0, Y0, s0)
    f6mul(t1, X1, Y1, s1)
    for k in range(3):
        f2add(S[k], X0[k], X1[k])                               # last read of in.y
    f6mul(C1, S, Ys, ss)
    for k in range(3):
        f2sub(C1[k], C1[k], t0[k]); f2sub(C1[k], C1[k], t1[k])
    x = slot("s")[0]
    f2xi(x, t1[2]); f2add(C0[0], t0[0], x)                      # C0 = t0 + v*t1
    f2add(C0[1], t0[1], t1[0]); f2add(C0[2], t0[2], t1[1])
    call("fr_add", (OUT, STATE_A), (IN, STATE_A), (ENT, ENTRY_ALPHA))
    call("fr_add", (OUT, STATE_B), (IN, STATE_B), (ENT, ENTRY_BETA))
    assert max(o for v in tmp.values() for _, o in v) + 64 <= FRAME
    L += ["ldp x19, x20, [sp, #0]", "ldp x21, x30, [sp, #16]", "add sp, sp, #%d" % FRAME, "ret"]
    return L


def program():
    return {"rebase": listing_rebase(), "f2mul": listing_f2mul(), "f2mulxi": listing_f2mulxi(),
            "fp_add": listing_fp_add("K"), "fp_sub": listing_fp_sub(), "fr_add": listing_fp_add("KR")}


def render(prog):
    out = []
    for name in ("rebase", "f2mul", "f2mulxi", "fp_add", "fp_sub", "fr_add"):
        out.append("%s:  // %d instructions" % (name, len(prog[name])))
        out += ["    " + s for s in prog[name]]
        out.append("")
    return "\n".join(out)


# ------------------------------------------------------------------ encodability
GPR = {"x%d" % i for i in range(31)}
CONDS = {"eq", "ne", "cs", "cc", "mi", "pl", "vs", "vc", "hi", "ls", "ge", "lt", "gt", "le"}


def check_encodable(line, labels):
    """Raise unless `line` is an encodable A64 instruction of the subset used.
    Immediate forms: ADD/SUB(S) imm12 0..4095 (register 31 = sp, except Rd of ADDS/SUBS = xzr);
    LDP/STP signed imm7*8 in [-512, 504]; LDR unsigned imm12*8 in [0, 32760].
    Register forms never name sp."""
    op, _, rest = line.partition(" ")
    if op == "ret":
        assert rest == "", line
        return
    if op in ("bl", "adr"):
        parts = [s.strip() for s in rest.split(",")]
        tgt = parts[-1]
        if op == "adr":
            assert parts[0] in GPR, line
        assert tgt in labels, line
        return
    if op in ("ldp", "stp", "ldr"):
        regs, _, mem = rest.partition("[")
        regs = [s.strip() for s in regs.strip().rstrip(",").split(",")]
        base, imm = [s.strip() for s in mem.rstrip("]").split(",")]
        assert base in GPR or base == "sp", line
        assert imm.startswith("#"), line
        v = int(imm[1:])
        assert all(r_ in GPR for r_ in regs), line
        if op == "ldr":
            assert len(regs) == 1 and v % 8 == 0 and 0 <= v <= 32760, line
        else:
            assert len(regs) == 2 and v % 8 == 0 and -512 <= v <= 504, line
            assert not (op == "ldp" and regs[0] == regs[1]), line
        return
    ops = [s.strip() for s in rest.split(",")]
    if op in ("add", "sub", "adds", "subs") and ops[-1].startswith("#"):
        v = int(ops[-1][1:])
        assert len(ops) == 3 and 0 <= v <= 4095, line
        rd_ok = (GPR | {"xzr"}) if op in ("adds", "subs") else (GPR | {"sp"})
        assert ops[0] in rd_ok and ops[1] in (GPR | {"sp"}), line
        return
    if op in ("add", "sub", "adds", "subs", "adc", "adcs", "sbcs", "and", "mul", "umulh"):
        assert len(ops) == 3 and all(o in GPR or o == "xzr" for o in ops), line
        return
    if op == "csel":
        assert len(ops) == 4 and all(o in GPR or o == "xzr" for o in ops[:3]) and ops[3] in CONDS, line
        return
    if op == "csetm":
        assert len(ops) == 2 and ops[0] in GPR and ops[1] in CONDS, line
        return
    raise AssertionError("not in the checked subset: " + line)


# ------------------------------------------------------------------ simulator
CALLER = ["x%d" % i for i in range(18)]
CALLEE = ["x%d" % i for i in range(19, 29)]


def _decode(line):
    op, _, rest = line.partition(" ")
    if op in ("ldp", "stp", "ldr"):
        regs, _, mem = rest.partition("[")
        regs = tuple(s.strip() for s in regs.strip().rstrip(",").split(","))
        base, imm = [s.strip() for s in mem.rstrip("]").split(",")]
        return (op, regs, base, int(imm[1:]))
    if op in ("ret",):
        return (op,)
    ops = tuple(s.strip() for s in rest.split(","))
    if op in ("add", "sub", "adds", "subs") and ops[-1].startswith("#"):
        return (op + "#", ops[0], ops[1], int(ops[2][1:]))
    return (op,) + ops


class Machine:
    """Executes the listings of `prog`.  Only the carry flag is modelled (the listings test no
    other flag).  bl writes a return address to x30 and ret jumps to the address held in x30,
    as on the hardware; the caller must preserve its own x30.  Every executed
    instruction, including bl and ret, is counted once.  Reading an unwritten address raises."""
    SENTINEL = 0xFFFF000000000000

    def __init__(self, prog, labels):
        self.code = {k: [_decode(s) for s in v] for k, v in prog.items()}
        self.labels = labels
        self.allowed = {k: set(CALLER) | set(CALLEE) | {"sp", "xzr"} for k in prog}
        self.allowed["rebase"] |= {"x30"}

    def run(self, entry, regs, mem):
        C, count = 0, 0
        fn, pc = entry, 0
        code, labels = self.code, self.labels
        rets = []
        regs["x30"] = self.SENTINEL
        allowed = self.allowed[fn]

        def rd(x):
            return 0 if x == "xzr" else regs[x]

        def wr(x, v):
            if x not in allowed:
                raise RuntimeError("%s writes %s" % (fn, x))
            if x != "xzr":
                regs[x] = v & M64

        while True:
            ins = code[fn][pc]
            pc += 1
            count += 1
            op = ins[0]
            if op == "ldp":
                a = rd(ins[2]) + ins[3]
                wr(ins[1][0], mem[a]); wr(ins[1][1], mem[a + 8])
            elif op == "stp":
                a = rd(ins[2]) + ins[3]
                mem[a] = rd(ins[1][0]); mem[a + 8] = rd(ins[1][1])
            elif op == "ldr":
                wr(ins[1][0], mem[rd(ins[2]) + ins[3]])
            elif op == "mul":
                wr(ins[1], rd(ins[2]) * rd(ins[3]))
            elif op == "umulh":
                wr(ins[1], (rd(ins[2]) * rd(ins[3])) >> 64)
            elif op in ("adds", "adcs", "adc", "add"):
                s = rd(ins[2]) + rd(ins[3]) + (C if op in ("adcs", "adc") else 0)
                wr(ins[1], s)
                if op != "adc" and op != "add":
                    C = s >> 64
            elif op in ("subs", "sbcs", "sub"):
                d = rd(ins[2]) - rd(ins[3]) - ((1 - C) if op == "sbcs" else 0)
                wr(ins[1], d)
                if op != "sub":
                    C = 1 if d >= 0 else 0
            elif op == "add#":
                wr(ins[1], rd(ins[2]) + ins[3])
            elif op == "sub#":
                wr(ins[1], rd(ins[2]) - ins[3])
            elif op == "subs#":
                d = rd(ins[2]) - ins[3]
                wr(ins[1], d)
                C = 1 if d >= 0 else 0
            elif op == "adds#":
                s = rd(ins[2]) + ins[3]
                wr(ins[1], s)
                C = s >> 64
            elif op == "csel":
                wr(ins[1], rd(ins[2]) if (C == 1) == (ins[4] == "cs") else rd(ins[3]))
            elif op == "csetm":
                wr(ins[1], M64 if (C == 1) == (ins[2] == "cs") else 0)
            elif op == "and":
                wr(ins[1], rd(ins[2]) & rd(ins[3]))
            elif op == "adr":
                wr(ins[1], labels[ins[2]])
            elif op == "bl":
                token = 0x7000000000 + 4 * len(rets)
                rets.append((fn, pc))
                regs["x30"] = token
                fn, pc = ins[1], 0
                allowed = self.allowed[fn]
            elif op == "ret":
                token = regs["x30"]
                if token == self.SENTINEL:
                    return count
                fn, pc = rets[(token - 0x7000000000) // 4]
                allowed = self.allowed[fn]
            else:
                raise RuntimeError("unsupported " + str(ins))


def count_by_function(prog):
    """Static call structure of rebase: number of bl to each leaf."""
    calls = {}
    for s in prog["rebase"]:
        if s.startswith("bl "):
            calls[s[3:]] = calls.get(s[3:], 0) + 1
    return calls
