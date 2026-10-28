"""Emit the two key kernels of the lazily evaluated walk (straight-line scalar AArch64, same subset,
register map, Montgomery form and reduction routine as gen.py):

  lazy.s  the class key of y' = x * Pi for a trie node Pi, without forming y': only the w^1 coefficient
          c_1(x * Pi) = sum_k x_k gamma_k, by Winograd's inner-product identity over Fp2 (pairs (u, v) in PAIRS)
            c_1 = sum (x_u + gamma_v)(x_v + gamma_u) - xi(x) - eta(Pi), xi = sum x_u x_v, eta = sum gamma_u gamma_v,
          xi per base (root.s), eta per node (stored); 3 Fp2 products by Karatsuba = 9 products of 384 x 384 bits,
          each operand the sum of two stored values in [0, p); then a = Re c_1, d = a - b, the twelve keys and
          their minimum K with its index kappa (proof.md section F2).
  root.s  the same key from an explicit x (a, b = coordinates 6, 7), plus the per-base data lazy.s reads for the
          rest of the block: the six pre-sums x_k.re + x_k.im mod p, and Xa, Xd (768 bits).

Im d_1 = -Re d_1 (asserted), so the odd-Frobenius base B = a Re d_1 + b Im d_1 = Re d_1 (a - b): one constant
product of d, and B zeta = d (Re d_1 zeta) a second.

Calling convention: x0 = out, x1 = walk state (x: 12 x 48 bytes, the pre-sums 6 x 48, Xa, Xd 2 x 96), x2 = node
entry (18 vectors gamma_k.re, gamma_k.im, gamma_k.re + gamma_k.im mod p, then Ea, Ed), x3 = constant pool (gen.py's,
plus Re d_1, Re d_1 zeta and zeta in Montgomery form). out: [0] a, [48] b, [96] B (d until B overwrites it),
[144] a zeta, [192] B zeta, [240] K (6 limbs), [288] kappa (one word).
"""
import json, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gen as G
import model as M

P, R = M.P, M.R
POOL = dict(G.POOL)
POOL.update({"Rd1": 2048, "Rd1z": 2096, "zeta": 2144})
RM = R % P
RINV = pow(R, -1, P)

def f2pow(a, e):
	r = (1, 0)
	while e:
		if e & 1:
			r = ((r[0] * a[0] - r[1] * a[1]) % P, (r[0] * a[1] + r[1] * a[0]) % P)
		a = ((a[0] * a[0] - a[1] * a[1]) % P, (2 * a[0] * a[1]) % P)
		e >>= 1
	return r

D1 = f2pow((1, 1), (P - 1) // 6)
assert (D1[0] + D1[1]) % P == 0
ZETA = f2pow((1, 1), (P * P - 1) // 6)
assert ZETA[1] == 0 and (ZETA[0] ** 2 - ZETA[0] + 1) % P == 0
ZETA = ZETA[0]
# Ranges (proof.md F5), asserted on integers. Stored values are in [0, p); pre-sums are reduced mod p (an
# unreduced one would let a pre-sum product reach 9 p^2 and d's reduction pass 8p). A lazy operand is at most
# 2p - 2 < 2^384, so S0, S1, S2 (three products each) are in [0, 3 (2p - 2)^2]. Per base (stored x): Xa = pR -
# (s0 - s1), Xd = pR - (2 s0 - s2), s_j in [0, 3 (p - 1)^2]. Per node: Ea, Ed in [p, 2p). lazy.s reduces a_T = S0 - S1
# + Xa + Ea 2^384 = a 2^384 mod p and d_T = 2 S0 - S2 + Xd + Ed 2^384 = d 2^384 mod p; Montgomery reduction returns
# r <= T/2^384 + p < 8p, inside the range of redc_store's q-table subtraction; a constant product is at most (p-1)^2.
OP = 2 * (P - 1)
SL = 3 * OP ** 2
SR = 3 * (P - 1) ** 2
XA = (P * R - SR, P * R + SR)
XD = (P * R - 2 * SR, P * R + SR)
EC = (P, 2 * P - 1)
TA = (-SL + XA[0] + EC[0] * R, SL + XA[1] + EC[1] * R)
TD = (-SL + XD[0] + EC[0] * R, 2 * SL + XD[1] + EC[1] * R)
assert OP < (1 << 384) and 2 * P < (1 << 384) and min(XA[0], XD[0], TA[0], TD[0]) >= 0
assert max(XA[1], XD[1], TA[1], TD[1]) < (1 << 768) and EC[1] < (1 << 384)
assert all(t // R + 1 + P <= 8 * P for t in (TA[1], TD[1], (P - 1) ** 2))

PAIRS = ((0, 1), (2, 3), (4, 5))   # Winograd pairs (u, v): operands x_u + gamma_v and x_v + gamma_u

def node_vector(pairs):
	"""Node entry from the six stored (re, im) of gamma_k, each in [0, p): 18 vectors (re, im, re + im mod p) and
    Ea, Ed with Ea 2^384 = -Re(eta) 2^768, Ed 2^384 = -(Re - Im)(eta) 2^768 mod p (stored products), in [p, 2p)."""
	vec = []
	for re, im in pairs:
		assert 0 <= re < P and 0 <= im < P
		vec += [re, im, (re + im) % P]
	e = [sum(vec[3 * u + s] * vec[3 * v + s] for u, v in PAIRS) for s in range(3)]
	return vec + [-(e[0] - e[1]) * RINV % P + P, -(2 * e[0] - e[2]) * RINV % P + P]

def block_vectors(xflat):
	"""Per-base data root.s writes after x (stored flat x in [0, p)): the six pre-sums (x_k.re + x_k.im mod p),
    then Xa and Xd as 768-bit integers, each as two 6-limb halves (low, high)."""
	w = [(xflat[fidx(k, 0)], xflat[fidx(k, 1)]) for k in range(6)]
	assert all(0 <= v < P for v in xflat)
	xs = [(re + im) % P for re, im in w]
	s = [sum(w[u][j] * w[v][j] for u, v in PAIRS) for j in range(2)] + [sum(xs[u] * xs[v] for u, v in PAIRS)]
	xa, xd = P * R - (s[0] - s[1]), P * R - (2 * s[0] - s[2])
	assert XA[0] <= xa <= XA[1] and XD[0] <= xd <= XD[1]
	return xs + [xa % R, xa // R, xd % R, xd // R]

def fidx(m, comp):
	"""blst flat index of component comp of the w^m coefficient (m = 2j + i -> 6i + 2j + comp)."""
	return 6 * (m % 2) + 2 * (m // 2) + comp

XS = 576          # walk state: the six pre-sums, then Xa (96 bytes) and Xd (96 bytes)
XA_OFF, XD_OFF = XS + 288, XS + 384
EA_OFF, ED_OFF = 864, 912          # node entry: Ea, Ed after the 18 vectors
OUT = {"a": 0, "b": 48, "B": 96, "az": 144, "Bz": 192, "K": 240, "kap": 288}

def load_consts(a):
	a.ins("mem", "ldp x0, x3, [sp, #96]")
	a.load6(G.PR, "x3", POOL["p"])
	a.ins("mem", f"ldr {G.N0R}, [x3, #{POOL['n0']}]")

def add_hi(a, off):
	"""Z += (the 6-limb multiple of p at pool offset off) * 2^384, on the upper six limbs only."""
	a.load6(G.T, "x3", off)
	a.chain("addsub", 1, G.Z[6:], G.Z[6:], G.T)

def mont_mul(a, src1, src2, dst):
	"""dst (out slot) = REDC(src1 * src2); src = (base, off). Clobbers Z, A, B, T."""
	a.set_pool(["x2"], shift=400)
	G.product(a, [(1,) + src1], *src2)
	load_consts(a)
	G.redc_store(a, dst // 48)

def sub_mod(a, dst, x, y):
	"""dst = x - y mod p for x, y in [0, p) (6 registers each, p in G.PR); clobbers G.T."""
	a.chain("keys", -1, dst, x, y)
	# x - y + p carries out of 2^384 exactly when x - y borrowed (2p < 2^384)
	a.chain("keys", 1, G.T, dst, G.PR)
	for i in range(6):
		a.ins("keys", f"csel {dst[i]}, {G.T[i]}, {dst[i]}, hs")

def keys_tail(a):
	"""a in out[a], d = a - b in out[B]: a zeta, B zeta, B, the twelve keys and their minimum K with its index
    kappa."""
	mont_mul(a, ("x0", OUT["a"]), ("x3", POOL["zeta"]), OUT["az"])
	mont_mul(a, ("x0", OUT["B"]), ("x3", POOL["Rd1z"]), OUT["Bz"])
	mont_mul(a, ("x0", OUT["B"]), ("x3", POOL["Rd1"]), OUT["B"])
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

def wprod(a, xa, ga, xb, gb):
	"""Z = (X_a + G_a)(X_b + G_b), operands (base, off) of stored values in [0, p): the first sum in A, the second
    in Z[7..11] and B[0] (b_j in Z[j + 7] for j < 5), then gen.product's rows."""
	BR = G.Z[7:12] + [G.B[0]]
	for dst, x, g in ((G.A, xa, ga), (BR, xb, gb)):
		a.load6(dst, *x)
		a.load6(G.T, *g)
		a.chain("addsub", 1, dst, dst, G.T)
	G.product(a, None, None, None, bregs=BR)

def xoff(k, s):
	"""Walk-state offset of x_k's component s (0 re, 1 im, 2 the stored pre-sum)."""
	return 48 * fidx(k, s) if s < 2 else XS + 48 * k

def store12(a, off, base="sp", load=False):
	for h in (0, 6):
		(a.load6 if load else a.store6)(G.Z[h:h + 6], base, off + 8 * h)

def sums(a, acc, prod, consts):
	"""S_s (s = 1, 2, then 0) over the three pairs i, prod(a, s, i) leaving each product in Z. S1 and S2 are
    stored at acc[s] minus the 768-bit constant and the upper-half constant consts[s] (mod 2^768); S0 stays in Z."""
	for s in (1, 2, 0):
		for i in range(3):
			prod(a, s, i)
			if i:
				a.acc12(1, "sp", acc[s])
			elif s in consts:
				(b12, o12), (b6, o6) = consts[s]
				a.acc12(-1, b12, o12)
				a.load6(G.T, b6, o6)
				a.chain("addsub", -1, G.Z[6:], G.Z[6:], G.T)
			if s or i < 2:
				store12(a, acc[s])

def lazy_kernel():
	a = G.Asm()
	acc = {0: G.SAVE, 1: G.SAVE + 96, 2: G.SAVE + 192}
	frame = (G.SAVE + 288 + 15) // 16 * 16
	a.comment("BLS12-381 lazy key: c_1(x * Pi) by Winograd's identity over Fp2 (3 Fp2 products, Karatsuba, 9 Fp\n"
			 "products) from x, its per-base data and the trie node entry, then the twelve keys, their minimum K and\n"
			 "index kappa; generated by gen_lazy.py; args x0 out, x1 walk state, x2 node entry, x3 constant pool")
	G.prologue(a, frame)
	a.set_pool(["x0", "x3"], shift=400)

	def prod(a, s, i):
		u, v = PAIRS[i]
		wprod(a, ("x1", xoff(u, s)), ("x2", 48 * (3 * v + s)), ("x1", xoff(v, s)), ("x2", 48 * (3 * u + s)))
	# S1 - Xa - Ea 2^384 and S2 - Xd - Ed 2^384 at acc[1], acc[2]; S0 in Z
	sums(a, acc, prod, {1: (("x1", XA_OFF), ("x2", EA_OFF)), 2: (("x1", XD_OFF), ("x2", ED_OFF))})
	# d_T = 2 S0 - (S2 - Xd - Ed 2^384), exact; S0 kept for a_T
	store12(a, acc[0])
	a.chain("addsub", 1, G.Z, G.Z, G.Z)
	a.acc12(-1, "sp", acc[2])
	load_consts(a)
	a.set_pool(["x2"], shift=400)
	G.redc_store(a, OUT["B"] // 48)
	# a_T = S0 - (S1 - Xa - Ea 2^384), exact
	store12(a, acc[0], load=True)
	a.acc12(-1, "sp", acc[1])
	G.redc_store(a, OUT["a"] // 48)
	# b = a - d mod p (a still in Z[:6], p in PR)
	a.load6(G.T, "x0", OUT["B"])
	sub_mod(a, G.Z[6:], G.Z[:6], G.T)
	a.store6(G.Z[6:], "x0", OUT["b"])
	keys_tail(a)
	G.epilogue(a, frame)
	return a, frame

def root_kernel():
	a = G.Asm()
	acc = {0: G.SAVE, 1: G.SAVE + 96, 2: G.SAVE + 192}
	frame = (G.SAVE + 288 + 15) // 16 * 16
	a.comment("BLS12-381 root key from an explicit x: a, b = coordinates 6, 7; the per-base data for lazy.s (six\n"
			 "pre-sums x_k.re + x_k.im mod p; Xa = pR - Re xi, Xd = pR - (Re - Im) xi, xi = sum_i x_2i x_2i+1, 9\n"
			 "products); d = a - b; then the keys as lazy.s; generated by gen_lazy.py")
	G.prologue(a, frame)
	a.set_pool(["x2"], shift=400)
	PZ = G.Z[6:]
	a.load6(PZ, "x3", POOL["p"])
	for k in range(6):
		# s = re + im < 2p < 2^384; s - p borrows exactly when s < p
		a.load6(G.A, "x1", xoff(k, 0))
		a.load6(G.T, "x1", xoff(k, 1))
		a.chain("addsub", 1, G.A, G.A, G.T)
		a.chain("addsub", -1, G.T, G.A, PZ)
		for i in range(6):
			a.ins("addsub", f"csel {G.A[i]}, {G.T[i]}, {G.A[i]}, hs")
		a.store6(G.A, "x1", xoff(k, 2))

	def prod(a, s, i):
		G.product(a, [(1, "x1", xoff(PAIRS[i][0], s))], "x1", xoff(PAIRS[i][1], s))
	sums(a, acc, prod, {})
	# Xd = pR - (2 s0 - s2), Xa = pR - (s0 - s1), exact
	store12(a, acc[0])
	a.chain("addsub", 1, G.Z, G.Z, G.Z)
	for s, off in ((2, XD_OFF), (1, XA_OFF)):
		if s == 1:
			store12(a, acc[0], load=True)
		for h in (0, 6):
			a.load6(G.T, "sp", acc[s] + 8 * h)
			a.chain("addsub", -1, G.Z[h:h + 6], G.T, G.Z[h:h + 6], first=not h)
		add_hi(a, POOL["p"])
		store12(a, off, "x1")
	a.load6(G.PR, "x3", POOL["p"])
	a.load6(G.Z[:6], "x1", 48 * fidx(1, 0))
	a.store6(G.Z[:6], "x0", OUT["a"])
	a.load6(G.Z[6:], "x1", 48 * fidx(1, 1))
	a.store6(G.Z[6:], "x0", OUT["b"])
	sub_mod(a, G.Z[6:], G.Z[:6], G.Z[6:])
	a.store6(G.Z[6:], "x0", OUT["B"])
	keys_tail(a)
	G.epilogue(a, frame)
	return a, frame

def pool_words():
	"""Constant pool image (64-bit words by byte offset) for the emulator and for the setup charge."""
	w = {}

	def put(off, v, n=6):
		for j in range(n):
			w[off + 8 * j] = (v >> (64 * j)) & ((1 << 64) - 1)
	for q in range(G.QMAX + 1):
		put(POOL["qt"] + 64 * q, q * P)
	w[POOL["n0"]] = G.N0
	w[POOL["c"]] = G.QC
	put(POOL["p"], P)
	put(POOL["Rd1"], D1[0] * RM % P)
	put(POOL["Rd1z"], D1[0] * ZETA * RM % P)
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
