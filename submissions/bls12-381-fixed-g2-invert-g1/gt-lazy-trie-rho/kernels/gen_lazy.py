"""Emit the two key kernels of the lazily evaluated walk (straight-line scalar AArch64, same subset,
register map, Montgomery form and reduction routine as gen.py):

  lazy.s  the class key of y' = x * Pi for a trie node Pi, WITHOUT forming y': only the w^1 coefficient
          c_1(x * Pi) = sum_k x_k gamma_k is computed, by Winograd's pairing
          c_1 = sum_{i<3} (x_2i + g_2i+1)(x_2i+1 + g_2i) - xi_x - xi_g, xi_x = sum_i x_2i x_2i+1 (stored per base
          by root.s), xi_g = sum_i g_2i g_2i+1 (stored per node): 3 Fp2 products, Karatsuba, 9 products of
          384 x 384 bits; then a = Re c_1, b = Im c_1, the twelve keys and their minimum K with its index kappa
          (proof.md sections F2, F5).
  root.s  the same key from an explicit x (a, b are coordinates 6 and 7 of x), plus what lazy.s reads for the
          rest of the block: the six Karatsuba pre-sums x_k.re + x_k.im and xi_x (9 products, 2 reductions).

Calling convention: x0 = out, x1 = walk state (x: 12 x 48 bytes, then the 6 pre-sums: 6 x 48 bytes, then
xi_x.re, xi_x.im: 2 x 48 bytes), x2 = node entry (lazy.s: 18 vectors gamma_k.re, gamma_k.im, gamma_k.re + gamma_k.im,
k = 0..5, then C_re = 4p - xi_g.re and C_im = 2p - xi_g.im), x3 = constant pool (gen.py's pool, plus Re d_1, Im d_1,
zeta in Montgomery form). All stored residues (x, gamma, xi_x, xi_g) are in [0, p).
out: [0] a, [48] b, [96] B, [144] a zeta, [192] B zeta, [240] K (6 limbs), [288] kappa (one word).
"""
import json, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gen as G
import model as M

P, R = M.P, M.R
POOL = dict(G.POOL)
POOL.update({"Rd1": 2048, "Id1": 2096, "zeta": 2144})
RM = R % P


def f2pow(a, e):
    r = (1, 0)
    while e:
        if e & 1:
            r = ((r[0] * a[0] - r[1] * a[1]) % P, (r[0] * a[1] + r[1] * a[0]) % P)
        a = ((a[0] * a[0] - a[1] * a[1]) % P, (2 * a[0] * a[1]) % P)
        e >>= 1
    return r


D1 = f2pow((1, 1), (P - 1) // 6)
ZETA = f2pow((1, 1), (P * P - 1) // 6)
assert ZETA[1] == 0 and (ZETA[0] ** 2 - ZETA[0] + 1) % P == 0
ZETA = ZETA[0]
# Ranges (proof.md F5), every operand and sum an exact non-negative integer. lazy.s: S0 and S1 operands
# x.re + g.re, x.im + g.im < 2p, S2 operands (pre-sum) + (node sum) < 4p < 2^384, so S0, S1 < 12p^2 and
# S2 < 48p^2 < 2^768. Im = S2 - S0 - S1 (exactly sum_i u0 v1 + u1 v0, in [0, 24p^2)) + (C_im - xi_x.im) 2^384 with
# C_im - xi_x.im in (0, 2p]; Re = S0 - S1 (in (-12p^2, 12p^2)) + (C_re - xi_x.re) 2^384 with C_re - xi_x.re in
# (2p, 4p] and 2p 2^384 >= 12p^2, so Re is in [0, 12p^2 + 4p 2^384). root.s: operands < 2p, S0, S1 < 3p^2,
# S2 < 12p^2; Im = S2 - S0 - S1 in [0, 6p^2), Re = S0 - S1 + p 2^384 in [0, 3p^2 + p 2^384). Montgomery reduction
# of T returns r < T/2^384 + 1 + p, below 8p, the range of redc_store's q-table subtraction. The in-register
# multiply-accumulate folds two carries into hi(A_5 b_j) <= A_5 - 1, safe while A_5 <= 2^64 - 2 (A < 4p).
assert 4 * P < (1 << 384) and (4 * P - 1) >> 320 <= (1 << 64) - 2
assert 48 * P * P < (1 << 768) and 2 * R >= 12 * P
assert 24 * P * P + 2 * P * R < (1 << 768) and 12 * P * P + 4 * P * R < (1 << 768)
assert (24 * P * P + 2 * P * R) // R + 1 + P <= 8 * P and (12 * P * P + 4 * P * R) // R + 1 + P <= 8 * P
assert 3 * P * P + P * R < (1 << 768) and (3 * P * P + P * R) // R + 1 + P <= 8 * P and 6 * P * P // R + 1 + P <= 8 * P
assert (2 * P * P) // R + 1 + P <= 8 * P


def fidx(m, comp):
    """blst flat index of component comp of the w^m coefficient (m = 2j + i -> 6i + 2j + comp)."""
    return 6 * (m % 2) + 2 * (m // 2) + comp


XS = 576          # offset of the six pre-sums in the walk state
OUT = {"a": 0, "b": 48, "B": 96, "az": 144, "Bz": 192, "K": 240, "kap": 288}


def load_consts(a):
    a.ins("mem", "ldp x0, x3, [sp, #96]")
    a.load6(G.PR, "x3", POOL["p"])
    a.ins("mem", f"ldr {G.N0R}, [x3, #{POOL['n0']}]")


def mont_mul(a, src1, src2, dst, acc_with=None, scratch=None):
    """dst (out slot) = REDC(src1 * src2 [+ acc_with product]); src = (base, off). Clobbers Z, A, B, T."""
    a.set_pool(["x2"], shift=400)
    G.product(a, [(1,) + src1], *src2)
    if acc_with is not None:
        a.store6(G.Z[:6], "sp", scratch)
        a.store6(G.Z[6:], "sp", scratch + 48)
        G.product(a, [(1,) + acc_with[0]], *acc_with[1])
        a.acc12(1, "sp", scratch)
    load_consts(a)
    G.redc_store(a, dst // 48)


def keys_tail(a, scratch):
    """B, a zeta, B zeta, the twelve keys and their minimum K with its index kappa; a, b in out[0], out[48]."""
    mont_mul(a, ("x0", OUT["a"]), ("x3", POOL["Rd1"]), OUT["B"], acc_with=(("x0", OUT["b"]), ("x3", POOL["Id1"])),
             scratch=scratch)
    mont_mul(a, ("x0", OUT["a"]), ("x3", POOL["zeta"]), OUT["az"])
    mont_mul(a, ("x0", OUT["B"]), ("x3", POOL["zeta"]), OUT["Bz"])
    # keys: family f (0: a, 1: B): v0 = value, v1 = value*zeta;
    # key_{2j+f} for j = 0..5: v0, v1, v1 - v0, -v0, -v1, v0 - v1 (mod p; negation as p - x).
    MN = G.Z[:6]           # running minimum
    IDX = G.Z[6]           # its index
    IMM = G.Z[7]
    C = G.T                # candidate
    PR = G.PR              # p (loaded)
    V0, V1 = G.Z[8:12] + [G.B[0], G.B[1]], None
    a.comment("keys and their minimum (strict <, so the first key attaining it in this order is kept)")
    a.load6(MN, "x0", OUT["a"])
    a.ins("keys", f"mov {IDX}, #0")

    def consider(idx):
        for i in range(6):
            a.ins("keys", f"{'subs' if i == 0 else 'sbcs'} xzr, {C[i]}, {MN[i]}")
        for i in range(6):
            a.ins("keys", f"csel {MN[i]}, {C[i]}, {MN[i]}, lo")
        a.ins("keys", f"mov {IMM}, #{idx}")
        a.ins("keys", f"csel {IDX}, {IMM}, {IDX}, lo")

    for f, (s0, s1) in enumerate(((OUT["a"], OUT["az"]), (OUT["B"], OUT["Bz"]))):
        # V0 <- v0 (kept in registers for the differences)
        a.load6(V0, "x0", s0)
        if f == 1:
            a.load6(C, "x0", s0)
            consider(1)
        # v1
        a.load6(C, "x0", s1)
        consider(2 + f)
        # -v1 = p - v1
        for i in range(6):
            a.ins("keys", f"{'subs' if i == 0 else 'sbcs'} {C[i]}, {PR[i]}, {C[i]}")
        consider(8 + f)
        # -v0 = p - v0
        for i in range(6):
            a.ins("keys", f"{'subs' if i == 0 else 'sbcs'} {C[i]}, {PR[i]}, {V0[i]}")
        consider(6 + f)
        # v1 - v0 mod p: C <- v1 - v0; if borrow add p
        a.load6(C, "x0", s1)
        for i in range(6):
            a.ins("keys", f"{'subs' if i == 0 else 'sbcs'} {C[i]}, {C[i]}, {V0[i]}")
        # V0 is dead now: V0 <- C + p, select on the borrow of the subtraction
        for i in range(6):
            a.ins("keys", f"{'adds' if i == 0 else 'adcs'} {V0[i]}, {C[i]}, {PR[i]}")
        # the flags now belong to the add; redo the borrow test on the stored values: v1 < v0 iff
        # the add above carried out (C - v0 wrapped exactly when the sum C + p is >= 2^384)
        for i in range(6):
            a.ins("keys", f"csel {C[i]}, {V0[i]}, {C[i]}, hs")
        consider(4 + f)
        # v0 - v1 = p - (v1 - v0) mod p
        for i in range(6):
            a.ins("keys", f"{'subs' if i == 0 else 'sbcs'} {C[i]}, {PR[i]}, {C[i]}")
        consider(10 + f)
    a.store6(MN, "x0", OUT["K"])
    a.ins("mem", f"str {IDX}, [x0, #{OUT['kap']}]")


BV = ["x23", "x24", "x25", "x26", "x27", "x28"]   # b-side operand held in registers
TT, TH, THN = "x29", "x30", "x18"                  # product temporaries


def load_sum(a, dst, terms):
    """dst = sum of the 6-limb vectors at terms [(base, off)], exact; the carry chain runs across the ldp pairs."""
    a.load6(dst, *terms[0])
    for base, off in terms[1:]:
        for i in range(0, 6, 2):
            a.ins("mem", f"ldp {TT}, {TH}, {a.ref(base, off + 8 * i)}")
            a.ins("addsub", f"{'adds' if i == 0 else 'adcs'} {dst[i]}, {dst[i]}, {TT}")
            a.ins("addsub", f"adcs {dst[i + 1]}, {dst[i + 1]}, {TH}")


def product_rr(a, mac):
    """Z = A * BV (mac False, 143 instructions) or Z += A * BV (mac True, 155), both operands in registers, by rows
    of BV. In the MAC form the carry out of row j's lo chain and of row j - 1 (both at limb j + 6) are folded into
    hi(A_5 b_j) <= A_5 - 1, which cannot overflow because A_5 < 2^63 (A < 4p, asserted); the carry out of the last
    row (limb 12) is 0 because the running sum stays below the asserted bound of the whole sum (< 2^768)."""
    A, Z, t = G.A, G.Z, TT
    if not mac:
        b = BV[0]
        for i in range(6):
            a.ins("prod", f"mul {Z[i]}, {A[i]}, {b}")
        for i in range(5):
            a.ins("prod", f"umulh {t}, {A[i]}, {b}")
            a.ins("prod", f"{'adds' if i == 0 else 'adcs'} {Z[i + 1]}, {Z[i + 1]}, {t}")
        a.ins("prod", f"umulh {t}, {A[5]}, {b}")
        a.ins("prod", f"adc {Z[6]}, {t}, xzr")
        for j in range(1, 6):
            b = BV[j]
            for i in range(6):
                a.ins("prod", f"mul {t}, {A[i]}, {b}")
                a.ins("prod", f"{'adds' if i == 0 else 'adcs'} {Z[j + i]}, {Z[j + i]}, {t}")
            a.ins("prod", f"adc {Z[j + 6]}, xzr, xzr")
            for i in range(5):
                a.ins("prod", f"umulh {t}, {A[i]}, {b}")
                a.ins("prod", f"{'adds' if i == 0 else 'adcs'} {Z[j + 1 + i]}, {Z[j + 1 + i]}, {t}")
            a.ins("prod", f"umulh {t}, {A[5]}, {b}")
            a.ins("prod", f"adc {Z[j + 6]}, {Z[j + 6]}, {t}")
        return
    h, hn = TH, THN
    a.ins("prod", f"umulh {h}, {A[5]}, {BV[0]}")
    for j in range(6):
        b = BV[j]
        for i in range(6):
            a.ins("prod", f"mul {t}, {A[i]}, {b}")
            a.ins("prod", f"{'adds' if i == 0 else 'adcs'} {Z[j + i]}, {Z[j + i]}, {t}")
        a.ins("prod", f"adc {h}, {h}, xzr")
        if j < 5:
            a.ins("prod", f"umulh {hn}, {A[5]}, {BV[j + 1]}")
        for i in range(5):
            a.ins("prod", f"umulh {t}, {A[i]}, {b}")
            a.ins("prod", f"{'adds' if i == 0 else 'adcs'} {Z[j + 1 + i]}, {Z[j + 1 + i]}, {t}")
        a.ins("prod", f"adcs {Z[j + 6]}, {Z[j + 6]}, {h}")
        if j < 5:
            a.ins("prod", f"adc {hn}, {hn}, xzr")
            h, hn = hn, h


def sum_of_products(a, pairs):
    """Z = sum over pairs of (sum of A terms) * (sum of B terms), the first product fresh, the rest accumulated."""
    for n, (aterms, bterms) in enumerate(pairs):
        load_sum(a, G.A, aterms)
        load_sum(a, BV, bterms)
        product_rr(a, n > 0)


def high_add(a, sign, base, off):
    """Z[6:12] +-= the 6-limb vector at base+off (adds c 2^384 to the 768-bit sum; carry out of limb 11 dropped)."""
    for i in range(0, 6, 2):
        a.ins("mem", f"ldp {TT}, {TH}, {a.ref(base, off + 8 * i)}")
        op = ("adds" if i == 0 else "adcs") if sign > 0 else ("subs" if i == 0 else "sbcs")
        op2 = "adcs" if sign > 0 else "sbcs"
        a.ins("addsub", f"{op} {G.Z[6 + i]}, {G.Z[6 + i]}, {TT}")
        a.ins("addsub", f"{op2} {G.Z[7 + i]}, {G.Z[7 + i]}, {TH}")


def wino_pairs(xa, xb, ga, gb):
    """Operand lists of the three Karatsuba sums S0, S1, S2 over the pairs (x_{2i} + g_{2i+1}, x_{2i+1} + g_{2i}),
    i < 3; xa(k, c), ga(k, c) give (base, off) of component c (0 re, 1 im, 2 re + im) of x_k, g_k (ga None: none)."""
    out = []
    for c in range(3):
        pr = []
        for i in range(3):
            at = [xa(2 * i, c)] + ([ga(2 * i + 1, c)] if ga else [])
            bt = [xb(2 * i + 1, c)] + ([gb(2 * i, c)] if gb else [])
            pr.append((at, bt))
        out.append(pr)
    return out


XI = XS + 288      # offset of the stored xi_x = sum_i x_{2i} x_{2i+1} (re, im) in the walk state
NODE_C = 864       # offset of the node constants C_re = 3p - xi_g.re, C_im = 2p - xi_g.im


def lazy_kernel():
    a = G.Asm()
    acc = {0: G.SAVE, 1: G.SAVE + 96}
    scratch = G.SAVE + 192
    frame = (scratch + 96 + 15) // 16 * 16
    a.comment("BLS12-381 lazy key: c_1(x * Pi) = sum_i (x_2i + g_2i+1)(x_2i+1 + g_2i) - xi_x - xi_g (Winograd), Karatsuba\n"
              "over Fp2 (9 products), then the twelve keys, their minimum K and index kappa; generated by gen_lazy.py\n"
              "args x0 out, x1 walk state (x, pre-sums, xi_x), x2 node entry (18 vectors, C_re, C_im), x3 constant pool")
    G.prologue(a, frame)
    a.set_pool(["x0", "x3"], shift=400)
    xo = lambda k, c: ("x1", XS + 48 * k if c == 2 else 48 * fidx(k, c))
    go = lambda k, c: ("x2", 48 * (3 * k + c))
    S = wino_pairs(xo, xo, go, go)
    for s in (0, 1, 2):
        sum_of_products(a, S[s])
        if s < 2:
            a.store6(G.Z[:6], "sp", acc[s])
            a.store6(G.Z[6:], "sp", acc[s] + 48)
    # Im = S2 - S0 - S1 + (C_im - xi_x.im) 2^384, in [0, 24p^2 + 2p 2^384)
    a.acc12(-1, "sp", acc[0])
    a.acc12(-1, "sp", acc[1])
    high_add(a, 1, "x2", NODE_C + 48)
    high_add(a, -1, "x1", XI + 48)
    # Re = S0 - S1 + (C_re - xi_x.re) 2^384, in [0, 12p^2 + 3p 2^384); x1, x2 still needed: reduce Im afterwards
    a.store6(G.Z[:6], "sp", acc[1] + 96)
    a.store6(G.Z[6:], "sp", acc[1] + 144)
    a.load6(G.Z[:6], "sp", acc[0])
    a.load6(G.Z[6:], "sp", acc[0] + 48)
    a.acc12(-1, "sp", acc[1])
    high_add(a, 1, "x2", NODE_C)
    high_add(a, -1, "x1", XI)
    load_consts(a)
    a.set_pool(["x2"], shift=400)
    G.redc_store(a, OUT["a"] // 48)
    a.load6(G.Z[:6], "sp", acc[1] + 96)
    a.load6(G.Z[6:], "sp", acc[1] + 144)
    G.redc_store(a, OUT["b"] // 48)
    keys_tail(a, scratch)
    G.epilogue(a, frame)
    return a, frame


def root_kernel():
    a = G.Asm()
    acc = {0: G.SAVE, 1: G.SAVE + 96}
    scratch = G.SAVE + 192
    frame = (scratch + 96 + 15) // 16 * 16
    a.comment("BLS12-381 root key from an explicit x: a, b = coordinates 6, 7; the six pre-sums x_k.re + x_k.im and\n"
              "xi_x = sum_i x_2i x_2i+1 (Karatsuba over Fp2, 9 products) for lazy.s; then the keys as lazy.s; generated by gen_lazy.py")
    G.prologue(a, frame)
    a.set_pool(["x2"], shift=400)
    for k in range(6):
        a.load6(G.A, "x1", 48 * fidx(k, 0))
        a.load6(G.T, "x1", 48 * fidx(k, 1))
        a.chain("addsub", 1, G.A, G.A, G.T)
        a.store6(G.A, "x1", XS + 48 * k)
    xo = lambda k, c: ("x1", XS + 48 * k if c == 2 else 48 * fidx(k, c))
    S = wino_pairs(xo, xo, None, None)
    for s in (0, 1, 2):
        sum_of_products(a, S[s])
        if s < 2:
            a.store6(G.Z[:6], "sp", acc[s])
            a.store6(G.Z[6:], "sp", acc[s] + 48)
    # Im = S2 - S0 - S1 in [0, 6p^2); Re = S0 - S1 + p 2^384 in [0, 3p^2 + p 2^384)
    a.acc12(-1, "sp", acc[0])
    a.acc12(-1, "sp", acc[1])
    load_consts(a)
    a.ins("addr", f"add x0, x1, #{XI}")
    G.redc_store(a, 1)
    a.load6(G.Z[:6], "sp", acc[0])
    a.load6(G.Z[6:], "sp", acc[0] + 48)
    a.acc12(-1, "sp", acc[1])
    a.chain("addsub", 1, G.Z[6:], G.Z[6:], G.PR, first=True)
    G.redc_store(a, 0)
    a.ins("mem", "ldr x0, [sp, #96]")
    for comp, slot in ((0, OUT["a"]), (1, OUT["b"])):
        a.load6(G.A, "x1", 48 * fidx(1, comp))
        a.store6(G.A, "x0", slot)
    keys_tail(a, scratch)
    G.epilogue(a, frame)
    return a, frame


def pool_words():
    """Constant pool image (64-bit words by byte offset) for the emulator and for the setup charge."""
    w = {}

    def put(off, v, n=6):
        for j in range(n):
            w[off + 8 * j] = (v >> (64 * j)) & ((1 << 64) - 1)
    for q in range(8):
        put(POOL["qt"] + 64 * q, q * P)
    w[POOL["n0"]] = G.N0
    w[POOL["c"]] = G.QC
    put(POOL["p"], P)
    put(POOL["Rd1"], D1[0] * RM % P)
    put(POOL["Id1"], D1[1] * RM % P)
    put(POOL["zeta"], ZETA * RM % P)
    return w


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    rep = {}
    for name, fn in (("lazy", lazy_kernel), ("root", root_kernel)):
        a, frame = fn()
        with open(os.path.join(here, f"{name}.s"), "w") as f:
            f.write(a.text())
        cnt = a.count()
        rep[name] = {"executed": sum(cnt.values()), "count": cnt, "frame": frame}
        print(name, sum(cnt.values()), cnt)
    with open(os.path.join(here, "lazy.layout.json"), "w") as f:
        json.dump({"pool": POOL, "out": OUT, "xs": XS, "kernels": rep}, f, indent=1)
    return rep


if __name__ == "__main__":
    main()
