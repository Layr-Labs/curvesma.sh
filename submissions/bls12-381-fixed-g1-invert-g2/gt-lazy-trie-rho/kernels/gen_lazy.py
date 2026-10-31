"""Emit the two key kernels of the lazily evaluated walk (straight-line scalar AArch64, same subset,
register map, Montgomery form and reduction routine as gen.py):

  lazy.s  the class key of y' = x * Pi for a trie node Pi, WITHOUT forming y': only the w^1 coefficient
          c_1(x * Pi) = sum_k x_k gamma_k (6 Fp2 products, Karatsuba, 18 products of 384 x 384 bits)
          is computed, from the node's stored vector gamma(Pi); then a = Re c_1, b = Im c_1, the twelve
          keys and their minimum K with its index kappa (proof.md section F2).
  root.s  the same key from an explicit x (a, b are coordinates 6 and 7 of x), plus the six Karatsuba
          pre-sums x_k.re + x_k.im that lazy.s reads for the rest of the block.

Calling convention: x0 = out, x1 = walk state (x: 12 x 48 bytes, then the 6 pre-sums: 6 x 48 bytes),
x2 = node entry (lazy.s: 18 vectors gamma_k.re, gamma_k.im, gamma_k.re + gamma_k.im, k = 0..5), x3 = constant
pool (gen.py's pool, plus Re d_1, Im d_1, zeta in Montgomery form and the 768-bit constant 6p^2).
out: [0] a, [48] b, [96] B, [144] a zeta, [192] B zeta, [240] K (6 limbs), [288] kappa (one word).
"""
import json, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gen as G
import model as M

P, R = M.P, M.R
POOL = dict(G.POOL)
POOL.update({"Rd1": 2048, "Id1": 2096, "zeta": 2144, "K6": 2192})
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
K6 = 6 * P * P
assert K6 < (1 << 768)
# Ranges (proof.md F5): stored inputs are in [0, p); a pre-sum x.re + x.im and a node sum g.re + g.im are below
# 2p < 2^384 (one 6-limb operand); S0, S1 < 6 p^2, S2 < 24 p^2; the reduced inputs Im = S2 - S0 - S1 (exactly the
# non-negative sum of x.re g.im + x.im g.re) and Re + 6p^2 = S0 - S1 + 6p^2 lie in [0, 12 p^2); B's input a Rd1 +
# b Id1 < 2 p^2 and a single product < p^2. Montgomery reduction of T < 12 p^2 returns r < T/2^384 + 1 + p, below
# 8p, the range of redc_store's q-table subtraction.
assert 2 * P < (1 << 384) and 24 * P * P < (1 << 768)
assert (12 * P * P) // R + 1 + P <= 8 * P and (2 * P * P) // R + 1 + P <= 8 * P

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

def lazy_kernel():
    a = G.Asm()
    acc = {0: G.SAVE, 1: G.SAVE + 96, 2: G.SAVE + 192}
    scratch = G.SAVE + 288
    frame = (scratch + 96 + 15) // 16 * 16
    a.comment("BLS12-381 lazy key: c_1(x * Pi) from the trie node vector gamma(Pi), Karatsuba over Fp2,\n"
              "then the twelve keys, their minimum K and index kappa; generated by gen_lazy.py\n"
              "args x0 out, x1 walk state (x, pre-sums), x2 node entry (18 vectors), x3 constant pool")
    G.prologue(a, frame)
    a.set_pool(["x0", "x3"], shift=400)
    # S0 = sum x_k.re gamma_k.re, S1 = sum x_k.im gamma_k.im, S2 = sum (x_k.re + x_k.im)(gamma_k.re + gamma_k.im)
    for s in (0, 1, 2):
        for k in range(6):
            if s < 2:
                xo = 48 * fidx(k, s)
            else:
                xo = XS + 48 * k
            go = 48 * (3 * k + s)
            G.product(a, [(1, "x1", xo)], "x2", go)
            if k == 0:
                a.store6(G.Z[:6], "sp", acc[s])
                a.store6(G.Z[6:], "sp", acc[s] + 48)
            else:
                a.acc12(1, "sp", acc[s])
                if not (s == 2 and k == 5):
                    a.store6(G.Z[:6], "sp", acc[s])
                    a.store6(G.Z[6:], "sp", acc[s] + 48)
    # Im = S2 - S0 - S1 (exactly the non-negative sum of x_re g_im + x_im g_re, below 12 p^2)
    a.acc12(-1, "sp", acc[0])
    a.acc12(-1, "sp", acc[1])
    load_consts(a)
    a.set_pool(["x2"], shift=400)
    G.redc_store(a, OUT["b"] // 48)
    # Re = S0 - S1 + 6p^2 (in [0, 12 p^2))
    a.load6(G.Z[:6], "sp", acc[0])
    a.load6(G.Z[6:], "sp", acc[0] + 48)
    a.acc12(-1, "sp", acc[1])
    a.acc12(1, "x3", POOL["K6"])
    G.redc_store(a, OUT["a"] // 48)
    keys_tail(a, scratch)
    G.epilogue(a, frame)
    return a, frame

def root_kernel():
    a = G.Asm()
    scratch = G.SAVE
    frame = (scratch + 96 + 15) // 16 * 16
    a.comment("BLS12-381 root key from an explicit x: a, b = coordinates 6, 7; the six pre-sums\n"
              "x_k.re + x_k.im for lazy.s; then the keys as lazy.s; generated by gen_lazy.py")
    G.prologue(a, frame)
    a.set_pool(["x2"], shift=400)
    for k in range(6):
        a.load6(G.A, "x1", 48 * fidx(k, 0))
        a.load6(G.T, "x1", 48 * fidx(k, 1))
        a.chain("addsub", 1, G.A, G.A, G.T)
        a.store6(G.A, "x1", XS + 48 * k)
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
    put(POOL["K6"], K6, n=12)
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
