"""Emit the key kernels of the lazily evaluated walk (straight-line scalar AArch64, same subset, register map,
Montgomery form and reduction routine as gen.py; proof.md F5):

  lazy.s  the class key of x * Pi for a trie node Pi without forming it: c_1(x Pi) = sum_k x_k gamma_k equals, on GT
          (W(x) = 2 x0 x3 - x1 x2 - xi x4 x5 = -x3), (U0 + A)(x3 + B) + (x1 + C)(x2 + D) + (U4 + E)(x5 + F) - eta,
          per-base U0 = -2 x0, U4 = xi x4, node constants A..F and eta: 3 Fp2 products by Karatsuba (9 products of
          384 x 384 bits), no per-base product; then a = Re c_1, d = a - b, the twelve keys, K and kappa. Off GT it
          returns the same integer function of its inputs.
  root.s  the key from an explicit x (a, b = coordinates 6, 7) and the per-base data lazy.s reads for the rest of
          the block: the six operand pre-sums (re + im mod p), then U0 and U4, each mod p.
  node.s  the node build's pre-sums and Ea, Ed (setup).

Im d_1 = -Re d_1 (asserted), so B = a Re d_1 + b Im d_1 = Re d_1 (a - b) and B zeta = d (Re d_1 zeta).
Calling convention: x0 = out, x1 = walk state (x: 12 x 48 bytes, then the per-base data: 10 x 48), x2 = node
entry (A..F as re, im, re + im mod p, then Ea, Ed), x3 = constant pool (gen.py's, plus Re d_1, Re d_1 zeta and zeta
in Montgomery form). out: [0] a, [48] b, [96] B (d until B overwrites it), [144] a zeta, [192] B zeta, [240] K
(6 limbs), [288] kappa (one word).
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
# Ranges (proof.md F5), asserted on integers. Every stored value (x, the per-base forms U0, U4 and pre-sums, the
# node's constants and pre-sums) is in [0, p), so a lazy operand is at most 2p - 2 < 2^384 and S0, S1, S2 (three
# products each) are in [0, 3 (2p - 2)^2]. Per node: Ea, Ed in [2p, 3p) (2p 2^384 > 3 (2p - 2)^2 keeps the sums
# non-negative). lazy.s reduces a_T = S0 - S1 + Ea 2^384 and d_T = 2 S0 - S2 + Ed 2^384; Montgomery reduction
# returns r <= T/2^384 + p < 8p, the range of redc_store's q-table subtraction; a constant product is at most (p-1)^2.
OP = 2 * (P - 1)
SL = 3 * OP ** 2
EOFF = 2 * P
EC = (EOFF, 3 * P - 1)
TA = (-SL + EC[0] * R, SL + EC[1] * R)
TD = (-SL + EC[0] * R, 2 * SL + EC[1] * R)
assert OP < (1 << 384) and EC[1] < (1 << 384) and min(TA[0], TD[0]) >= 0 and max(TA[1], TD[1]) < (1 << 768)
# node.s: Ta = pR - (e0 - e1), Td = pR - (2 e0 - e2), e_s = eta's three products on component s in [0, 3 (p - 1)^2].
SR = 3 * (P - 1) ** 2
NA, ND = (P * R - SR, P * R + SR), (P * R - 2 * SR, P * R + SR)
assert min(NA[0], ND[0]) >= 0 and max(NA[1], ND[1]) < (1 << 768)
assert all(t // R + 1 + P <= 8 * P for t in (TA[1], TD[1], NA[1], ND[1], (P - 1) ** 2))

# Operands n = 2i + j of pair i (j = 0 left, 1 right): c_1 = (U0 + A)(x3 + B) + (x1 + C)(x2 + D) + (U4 + E)(x5 + F)
# - eta on GT, U0 = -2 x0, U4 = xi x4 (per-base forms), A..F node constants, eta = AB + CD + EF (proof.md F2).
OPND = (("U", 0), ("x", 3), ("x", 1), ("x", 2), ("U", 1), ("x", 5))

def node_vector(g):
	"""Node entry from the six stored gamma_k = (re, im), each in [0, p): constants A = g3 - 1, B = -g0/2, C = g2,
    D = g1, E = g5, F = g4/xi as 18 vectors (re, im, re + im mod p), then entry()'s Ea, Ed."""
	assert all(0 <= v < P for c in g for v in c)
	h = lambda v: v * (P + 1) // 2 % P
	c = [((g[3][0] - RM) % P, g[3][1]), (h(-g[0][0]), h(-g[0][1])), g[2], g[1], g[5],
		(h(g[4][0] + g[4][1]), h(g[4][1] - g[4][0]))]
	return entry([v for re, im in c for v in (re, im, (re + im) % P)])

def entry(vec):
	"""The 18 vectors and Ea, Ed in [2p, 3p): Ea 2^384 = -Re(eta) 2^768, Ed 2^384 = -(Re - Im)(eta) 2^768 mod p,
    eta = sum_i (constant 2i)(constant 2i + 1) on the stored values."""
	e = [sum(vec[6 * i + s] * vec[6 * i + 3 + s] for i in range(3)) for s in range(3)]
	return vec + [-(e[0] - e[1]) * RINV % P + EOFF, -(2 * e[0] - e[2]) * RINV % P + EOFF]

def block_vectors(xflat):
	"""Per-base data root.s writes after x (stored flat x in [0, p)): the six operand pre-sums, then U0, U4."""
	assert all(0 <= v < P for v in xflat)
	x = [(xflat[fidx(k, 0)], xflat[fidx(k, 1)]) for k in range(6)]
	u = [((-2 * x[0][0]) % P, (-2 * x[0][1]) % P), ((x[4][0] - x[4][1]) % P, (x[4][0] + x[4][1]) % P)]
	ops = [u[v] if k == "U" else x[v] for k, v in OPND]
	return [(re + im) % P for re, im in ops] + [v for c in u for v in c]

def fidx(m, comp):
	"""blst flat index of component comp of the w^m coefficient (m = 2j + i -> 6i + 2j + comp)."""
	return 6 * (m % 2) + 2 * (m // 2) + comp

XS = 576          # walk state: the six operand pre-sums, then U0, U4 (re, im)
UF = XS + 288
EA_OFF, ED_OFF = 864, 912          # node entry: Ea, Ed after the 18 vectors
OUT = {"a": 0, "b": 48, "B": 96, "az": 144, "Bz": 192, "K": 240, "kap": 288}

def load_consts(a):
	a.ins("mem", "ldp x0, x3, [sp, #96]")
	a.load6(G.PR, "x3", POOL["p"])
	a.ins("mem", f"ldr {G.N0R}, [x3, #{POOL['n0']}]")

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

def woff(n, s):
	"""Walk-state offset of operand n's component s (0 re, 1 im, 2 its stored pre-sum)."""
	k, v = OPND[n]
	return XS + 48 * n if s == 2 else (48 * fidx(v, s) if k == "x" else UF + 96 * v + 48 * s)

def store12(a, off, base="sp", load=False):
	for h in (0, 6):
		(a.load6 if load else a.store6)(G.Z[h:h + 6], base, off + 8 * h)

def sums(a, acc, prod, sub=None):
	"""S_s (s = 1, 2, then 0) over the three pairs, prod(s, i) leaving each product in Z; S1 and S2 stored at acc[s]
    (minus the node offset at x2 + sub[s - 1], times 2^384, when sub is given); S0 stays in Z."""
	for s in (1, 2, 0):
		for i in range(3):
			prod(s, i)
			if i:
				a.acc12(1, "sp", acc[s])
			elif s and sub:
				a.load6(G.T, "x2", sub[s - 1])
				a.chain("addsub", -1, G.Z[6:], G.Z[6:], G.T)
			if s or i < 2:
				store12(a, acc[s])

ACC = {0: G.SAVE, 1: G.SAVE + 96, 2: G.SAVE + 192}
FRAME = (G.SAVE + 288 + 15) // 16 * 16

def lazy_kernel():
	a = G.Asm()
	a.comment("BLS12-381 lazy key: c_1(x * Pi) on GT as three Fp2 products of (per-base operand + node constant)\n"
			 "(Karatsuba, 9 Fp products) minus the node's eta, then the twelve keys, their minimum K and index kappa;\n"
			 "generated by gen_lazy.py; args x0 out, x1 walk state, x2 node entry, x3 constant pool")
	G.prologue(a, FRAME)
	a.set_pool(["x0", "x3"], shift=0)
	sums(a, ACC, lambda s, i: wprod(a, ("x1", woff(2 * i, s)), ("x2", 288 * i + 48 * s), ("x1", woff(2 * i + 1, s)),
									("x2", 288 * i + 144 + 48 * s)), (EA_OFF, ED_OFF))
	# d_T = 2 S0 - (S2 - Ed 2^384), exact; S0 kept for a_T
	store12(a, ACC[0])
	a.chain("addsub", 1, G.Z, G.Z, G.Z)
	a.acc12(-1, "sp", ACC[2])
	load_consts(a)
	a.set_pool(["x2"], shift=400)
	G.redc_store(a, OUT["B"] // 48)
	# a_T = S0 - (S1 - Ea 2^384), exact
	store12(a, ACC[0], load=True)
	a.acc12(-1, "sp", ACC[1])
	G.redc_store(a, OUT["a"] // 48)
	# b = a - d mod p (a still in Z[:6], p in PR)
	a.load6(G.T, "x0", OUT["B"])
	sub_mod(a, G.Z[6:], G.Z[:6], G.T)
	a.store6(G.Z[6:], "x0", OUT["b"])
	keys_tail(a)
	G.epilogue(a, FRAME)
	return a, FRAME

def node_kernel():
	"""The trie-node build's key-entry sequence (proof.md section 5), x0 = x1 = node entry holding the constants'
    (re, im): their six pre-sums mod p, eta's sums e_s, Ed = REDC(pR - (2 e0 - e2)) + 2p, Ea = REDC(pR - (e0 - e1)) + 2p."""
	a = G.Asm()
	G.prologue(a, FRAME)
	a.set_pool(["x2"], shift=400)
	PZ = G.Z[6:]
	a.load6(PZ, "x3", POOL["p"])
	for n in range(6):
		# s = re + im < 2p < 2^384; s - p borrows exactly when s < p
		a.load6(G.A, "x1", 144 * n)
		a.load6(G.T, "x1", 144 * n + 48)
		a.chain("addsub", 1, G.A, G.A, G.T)
		a.chain("addsub", -1, G.T, G.A, PZ)
		for i in range(6):
			a.ins("addsub", f"csel {G.A[i]}, {G.T[i]}, {G.A[i]}, hs")
		a.store6(G.A, "x1", 144 * n + 96)
	sums(a, ACC, lambda s, i: G.product(a, [(1, "x1", 288 * i + 48 * s)], "x1", 288 * i + 144 + 48 * s))
	store12(a, ACC[0])
	a.chain("addsub", 1, G.Z, G.Z, G.Z)
	for s, off in ((2, ED_OFF), (1, EA_OFF)):
		if s == 1:
			store12(a, ACC[0], load=True)
		for h in (0, 6):
			a.load6(G.T, "sp", ACC[s] + 8 * h)
			a.chain("addsub", -1, G.Z[h:h + 6], G.T, G.Z[h:h + 6], first=not h)
		a.load6(G.T, "x3", POOL["p"])
		a.chain("addsub", 1, G.Z[6:], G.Z[6:], G.T)
		load_consts(a)
		G.redc_store(a, off // 48)
		a.load6(G.T, "x3", POOL["2p"])
		a.chain("addsub", 1, G.Z[:6], G.Z[:6], G.T)
		a.store6(G.Z[:6], "x0", off)
	G.epilogue(a, FRAME)
	return a, FRAME

def root_kernel():
	a = G.Asm()
	frame = (G.SAVE + 15) // 16 * 16
	a.comment("BLS12-381 root key from an explicit x: a, b = coordinates 6, 7; the per-base data for lazy.s (U0 = -2 x0,\n"
			 "U4 = xi x4 and the six operand pre-sums, each mod p); d = a - b; then the keys as lazy.s; generated by\n"
			 "gen_lazy.py")
	G.prologue(a, frame)
	a.set_pool(["x2"], shift=400)
	X, Y, PZ = G.A, G.Z[:6], G.Z[6:]
	a.load6(PZ, "x3", POOL["p"])

	def modop(sign, x, y, off):
		"""Y = x +- y mod p (x, y in [0, p); y may be Y), stored at walk-state offset off: T = x +- y exactly, then
        Y = T -+ p, kept when it did not borrow (sign +) or carried out of 2^384 (sign -, T borrowed)."""
		a.chain("addsub", sign, G.T, x, y)
		a.chain("addsub", -sign, Y, G.T, PZ)
		for i in range(6):
			a.ins("addsub", f"csel {Y[i]}, {Y[i]}, {G.T[i]}, hs")
		a.store6(Y, "x1", off)

	wx = lambda k, s: 48 * fidx(k, s)
	for n, (k, v) in enumerate(OPND):
		if k == "x":
			a.load6(X, "x1", wx(v, 0))
			a.load6(Y, "x1", wx(v, 1))
			modop(1, X, Y, woff(n, 2))
	# U4 = (x4.re - x4.im, x4.re + x4.im), its pre-sum 2 x4.re
	a.load6(X, "x1", wx(4, 0))
	for sign, comp in ((-1, 0), (1, 1)):
		a.load6(Y, "x1", wx(4, 1))
		modop(sign, X, Y, woff(4, comp))
	modop(1, X, X, woff(4, 2))
	# U0 = -2 x0 = 0 - (2 x0 mod p) per component, its pre-sum U0.re + U0.im
	for comp in (0, 1):
		a.load6(X, "x1", wx(0, comp))
		modop(1, X, X, woff(0, comp))
		modop(-1, ["xzr"] * 6, Y, woff(0, comp))
	a.load6(X, "x1", woff(0, 0))
	modop(1, X, Y, woff(0, 2))
	a.load6(G.PR, "x3", POOL["p"])
	a.load6(G.Z[:6], "x1", wx(1, 0))
	a.store6(G.Z[:6], "x0", OUT["a"])
	a.load6(G.Z[6:], "x1", wx(1, 1))
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
	for q in range(8):
		put(POOL["qt"] + 64 * q, q * P)
	w[POOL["n0"]] = G.N0
	w[POOL["c"]] = G.QC
	put(POOL["p"], P)
	put(POOL["2p"], 2 * P)
	put(POOL["Rd1"], D1[0] * RM % P)
	put(POOL["Rd1z"], D1[0] * ZETA * RM % P)
	put(POOL["zeta"], ZETA * RM % P)
	return w

def main():
	here = os.path.dirname(os.path.abspath(__file__))
	rep = {}
	for name, fn in (("lazy", lazy_kernel), ("root", root_kernel), ("node", node_kernel)):
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
