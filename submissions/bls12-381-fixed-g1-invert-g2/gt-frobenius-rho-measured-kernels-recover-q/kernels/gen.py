"""Emit straight-line scalar AArch64 for two BLS12-381 Fp12 kernels (blst's tower, Montgomery form,
6 x 64-bit limbs per Fp element, outputs fully reduced):

  mul.s  out = x * m, m a fixed table entry stored as 45 Fp-linear forms of m reduced mod p (Toom-3 over Fp6)
  sqr.s  out = x^2 for x in the cyclotomic subgroup (Granger-Scott, as blst's cyclotomic_sqr_fp12)

Calling convention of both: x0 = out (12 x 48 bytes), x1 = x (12 x 48 bytes), x2 = table entry (mul),
x3 = constant pool (layout.json). General registers only, no branches, no calls; ends in ret.

Lazy reduction: each output coordinate is one 768-bit integer T = K + sum(+-products) (+ 2 a R terms
for the squaring); every operand is an exact sum of stored values, each with a proven bound (inputs and entry
values in [0, p), mul's reduced instance-0 forms in [0, p) or [0, 2p)), so T is an integer polynomial whose
range model.py bounds monomial by monomial (K: least multiple of p 2^384 making T >= 0 for mul, of p for sqr);
one Montgomery reduction, then subtract q p and p.
"""
import json, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import model as M

P, R = M.P, M.R
N0 = (-pow(P, -1, 1 << 64)) % (1 << 64)
Z = [f"x{i}" for i in range(4, 16)]
A = ["x16", "x17", "x19", "x20", "x21", "x22"]
B = ["x23", "x24"]
T = [f"x{i}" for i in range(25, 31)]
PR, N0R, MR = A, "x23", "x24"          # registers of the reduction phase
SAVE = 112                               # callee-saved x19..x30 (96 bytes) and x0, x3 (16 bytes)

# constant pool: the q*p table at 64q (q = 0..9), n0, c, p, 2p, 4p, then the per-kernel offsets K_k
QMAX = 9
POOL = {"qt": 0, "n0": 48, "c": 56, "p": 64, "2p": 128, "4p": 256, "K": 64 * (QMAX + 1)}
QC = (1 << 64) // ((P >> 320) + 1)      # q = r5 * QC >> 64 is monotone in r: check both ends of each [kp, (k+1)p)
assert QC == 9 and 9 * P < (1 << 384) < 10 * P
assert all(k - 1 <= (k * P >> 320) * QC >> 64 and ((k + 1) * P - 1 >> 320) * QC >> 64 <= k for k in range(10))

class Asm:
	def __init__(self):
		self.lines, self.cats, self.pool, self.addr, self.clock = [], [], [], {}, 0

	def ins(self, cat, text):
		op = text.split()[0]
		if op in ("ldp", "stp", "ldr", "str"):
			cat = "mem"
		self.lines.append("    " + text)
		self.cats.append(cat)

	def comment(self, text):
		for t in text.split("\n"):
			self.lines.append("// " + t)
			self.cats.append(None)

	def set_pool(self, regs, shift=None):
		self.pool, self.addr = list(regs), {}
		self.shift = shift if shift is not None else getattr(self, "shift", 48)

	def ref(self, base, off, pair=True):
		lo, hi = (-512, 504) if pair else (0, 32760)
		if lo <= off <= hi:
			return f"[{base}, #{off}]"
		self.clock += 1
		for r, (b, o, _) in list(self.addr.items()):
			if b == base and lo <= off - o <= hi:
				self.addr[r] = (b, o, self.clock)
				return f"[{r}, #{off - o}]"
		free = [r for r in self.pool if r not in self.addr]
		r = free[0] if free else min(self.addr, key=lambda k: self.addr[k][2])
		new = off + self.shift
		if new > 4095:
			self.ins("addr", f"add {r}, {base}, #{new >> 12}, lsl #12")
			if new & 4095:
				self.ins("addr", f"add {r}, {r}, #{new & 4095}")
		else:
			self.ins("addr", f"add {r}, {base}, #{new}")
		self.addr[r] = (base, new, self.clock)
		return f"[{r}, #{off - new}]"

	def load6(self, regs, base, off):
		for i in range(0, 6, 2):
			self.ins("mem", f"ldp {regs[i]}, {regs[i + 1]}, {self.ref(base, off + 8 * i)}")

	def store6(self, regs, base, off):
		for i in range(0, 6, 2):
			self.ins("mem", f"stp {regs[i]}, {regs[i + 1]}, {self.ref(base, off + 8 * i)}")

	def chain(self, cat, sign, dst, x, y, first=True):
		"""dst = x +- y over len(dst) limbs; first=False continues the carry already in C."""
		for i in range(len(dst)):
			if sign > 0:
				op = "adds" if (i == 0 and first) else "adcs"
			else:
				op = "subs" if (i == 0 and first) else "sbcs"
			self.ins(cat, f"{op} {dst[i]}, {x[i]}, {y[i]}")

	def acc12(self, sign, base, off, cat="addsub"):
		"""Z (12 limbs) +-= the 768-bit value at base+off, through T in two halves (flags survive ldp)."""
		self.load6(T, base, off)
		self.chain(cat, sign, Z[:6], Z[:6], T)
		self.load6(T, base, off + 48)
		self.chain(cat, sign, Z[6:], Z[6:], T, first=False)

	def text(self):
		return "\n".join(self.lines) + "\n"

	def count(self):
		out = {}
		for c in self.cats:
			if c is not None:
				out[c] = out.get(c, 0) + 1
		return out

def prologue(a, frame):
	left = frame
	while left:
		step = min(left, 4080)
		a.ins("glue", f"sub sp, sp, #{step}")
		left -= step
	regs = ["x19", "x20", "x21", "x22", "x23", "x24", "x25", "x26", "x27", "x28", "x29", "x30"]
	for i in range(0, 12, 2):
		a.ins("mem", f"stp {regs[i]}, {regs[i + 1]}, [sp, #{8 * i}]")
	a.ins("mem", "stp x0, x3, [sp, #96]")

def epilogue(a, frame):
	regs = ["x19", "x20", "x21", "x22", "x23", "x24", "x25", "x26", "x27", "x28", "x29", "x30"]
	for i in range(0, 12, 2):
		a.ins("mem", f"ldp {regs[i]}, {regs[i + 1]}, [sp, #{8 * i}]")
	left = frame
	while left:
		step = min(left, 4080)
		a.ins("glue", f"add sp, sp, #{step}")
		left -= step
	a.ins("glue", "ret")

def load_lin(a, dst, terms):
	"""dst = sum of sign * vector over `terms` [(sign, base, off)], exact (range checked by the caller)."""
	s0, b0, o0 = terms[0]
	assert s0 > 0
	a.load6(dst, b0, o0)
	for s, b, o in terms[1:]:
		a.load6(T, b, o)
		a.chain("addsub", s, dst, dst, T)

def product(a, aterms, bbase, boff, double_a=False, bregs=None):
	"""Z = (sum of aterms) * (vector at bbase+boff), schoolbook by rows of b (143 arithmetic instructions);
    with bregs, A and b are already in A and bregs (b_j must survive rows 0..j: row j writes Z[j..j+6], T)."""
	if not bregs:
		load_lin(a, A, aterms)
	if double_a:
		a.chain("addsub", 1, A, A, A)
	for jj in range(0, 6, 2):
		if not bregs:
			a.ins("mem", f"ldp {B[0]}, {B[1]}, {a.ref(bbase, boff + 8 * jj)}")
		for j, b in (zip((jj, jj + 1), bregs[jj:jj + 2]) if bregs else ((jj, B[0]), (jj + 1, B[1]))):
			if j == 0:
				for i in range(6):
					a.ins("prod", f"mul {Z[i]}, {A[i]}, {b}")
				for i in range(6):
					a.ins("prod", f"umulh {T[i]}, {A[i]}, {b}")
				a.chain("prod", 1, Z[1:6], Z[1:6], T[:5])
				a.ins("prod", f"adc {Z[6]}, {T[5]}, xzr")
			else:
				for i in range(6):
					a.ins("prod", f"mul {T[i]}, {A[i]}, {b}")
				a.chain("prod", 1, Z[j:j + 6], Z[j:j + 6], T)
				a.ins("prod", f"adc {Z[j + 6]}, xzr, xzr")
				for i in range(6):
					a.ins("prod", f"umulh {T[i]}, {A[i]}, {b}")
				a.chain("prod", 1, Z[j + 1:j + 6], Z[j + 1:j + 6], T[:5])
				a.ins("prod", f"adc {Z[j + 6]}, {Z[j + 6]}, {T[5]}")

def redc_store(a, k):
	"""Montgomery-reduce Z (T < 2^768) to r = (T + M p)/2^384 < T/2^384 + p, then subtract the
    estimated multiple q p (r < 2^384 < 10p) and p branch-free, and store coordinate k. Window registers
    rotate; the low limb of each round becomes the next top (its sum with lo(m p0) is 0 mod 2^64,
    carry = [limb != 0])."""
	W = Z[:6]
	for _ in range(6):
		a.ins("redc", f"mul {MR}, {W[0]}, {N0R}")
		for j in range(1, 6):
			a.ins("redc", f"mul {T[j]}, {MR}, {PR[j]}")
		a.ins("redc", f"subs xzr, {W[0]}, #1")
		for j in range(1, 6):
			a.ins("redc", f"adcs {W[j]}, {W[j]}, {T[j]}")
		a.ins("redc", f"adc {W[0]}, xzr, xzr")
		for j in range(6):
			a.ins("redc", f"umulh {T[j]}, {MR}, {PR[j]}")
		a.chain("redc", 1, W[1:6], W[1:6], T[:5])
		a.ins("redc", f"adc {W[0]}, {W[0]}, {T[5]}")
		W = W[1:] + W[:1]
	a.chain("redc", 1, W, W, Z[6:])
	assert W == Z[:6]                             # r < 10p; q = floor(r5 * QC / 2^64) in [floor(r/p) - 1, floor(r/p)]
	q, ad = Z[6], Z[7]
	a.ins("csub", f"ldr {q}, [x3, #{POOL['c']}]")
	a.ins("csub", f"umulh {q}, {W[5]}, {q}")
	a.ins("csub", f"add {ad}, x3, {q}, lsl #6")
	a.load6(T, ad, POOL["qt"])
	a.chain("csub", -1, W, W, T)
	a.chain("csub", -1, T, W, PR)
	for j in range(6):
		a.ins("csub", f"csel {W[j]}, {W[j]}, {T[j]}, lo")
	a.store6(W, "x0", 48 * k)

def out_bounds(polys, vmax):
	"""K = least multiple of p with K + lo >= 0."""
	lohi = [M.poly_range(poly, vmax) for poly in polys]
	Ks = [-(-(-lo) // P) * P if lo < 0 else 0 for lo, hi in lohi]
	rows = []
	for k, ((lo, hi), K) in enumerate(zip(lohi, Ks)):
		top = K + hi
		assert K + lo >= 0 and top < (1 << 768)
		bres = top // R + 1 + P            # r < T/R + p
		assert bres <= 8 * P
		rows.append((k, lo, hi, K, bres))
	return rows

def bound_table(rows, cap="8p"):
	t = ["output  lo/p^2    hi/p^2    K/p^2   T_max/p^2  r bound/p",
		"(T = K + polynomial in [K + lo, K + hi], lo/hi summed monomial by monomial, inputs in [0, p))"]
	for k, lo, hi, K, bres in rows:
		t.append(f"c{k:<6} {lo / P / P:+8.3f} {hi / P / P:+9.3f} {K / P / P:8.3f} {(K + hi) / P / P:10.3f} "
				f"{bres / P:10.4f}")
	t.append(f"(each r < {cap}: subtract q p, q = umulh(r5, 9) read from the q*p table, then p unless that borrows)")
	return "\n".join(t)

# -- multiplication by a table entry: Toom-3 over Fp6 at the Gaussian points t = 0, 1, -1, -u, inf
# Fp6 y n (degree 2 each in v) is interpolated from the five Fp2 products y(t) n(t), the point inf reading the
# top coefficients. The entry stores B_t = n(t)/lam_t (lam = -1, -2-2u, -2-2u, -2, -1), so every weight COEF[o][t]
# of W_t = y(t) B_t in output o of y n mod (v^3 - xi) is a Gaussian integer of norm at most 2 (each part -1, 0
# or 1). GPOW[t]: t^0, t^1, t^2. PAIRS: S = W_x - W_y, built once per Fp12 output half and read by several outputs.
GPOW = (((1, 0), (0, 0), (0, 0)), ((1, 0), (1, 0), (1, 0)), ((1, 0), (-1, 0), (1, 0)), ((1, 0), (0, -1), (-1, 0)), None)
COEF = (((0, -1), (-1, -1), (-1, 1), (-1, 1), (1, -1)), ((0, 1), (0, -1), (1, 0), (0, -1), (-1, 0)),
		((1, 0), (-1, -1), (-1, -1), (0, 0), (1, 0)))
PAIRS = (((4, 0), (2, 0)), ((4, 1), (2, 1)), ((1, 0), (0, 1)), ((1, 0), (0, 0)), ((3, 0), (1, 1)))
PARTIAL = {(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (3, 0)}  # instance-0 forms (t, k) left in [0, 2p), not [0, p)

def at(y, t):
	"""y(t), exact integer Lins, for an Fp6 y given as three Fp2 Lin pairs."""
	if GPOW[t] is None:
		return y[2]
	re, im = {}, {}
	sc = lambda L, k: {n: k * c for n, c in L.items() if k * c}
	for (a, b), (f0, f1) in zip(GPOW[t], y):
		re, im = M.ladd(re, M.ladd(sc(f0, a), sc(f1, b), -1)), M.ladd(im, M.ladd(sc(f1, a), sc(f0, b)))
	return re, im

def entry_program():
	"""The 45 entry vectors of a table value from its 12 raw coordinates m, by the straight-line program the table
  build runs: ops[0] counts modular additions, subtractions and negations, ops[1] modular halvings. Instance
  m-sides n = m0, v m1 - m0, m1 - m0; per point B_t = n(t)/lam_t (1/lam = -1, (u-1)/4, (u-1)/4, -1/2, -1),
  stored as B1, B0 + B1, B0 - B1. Returns [instance][point] triples of Lins with coefficients mod p."""
	ops = [0, 0]
	def ad(a, b, s=1):
		ops[0] += 1
		return {n: c % P for n, c in M.ladd(a, b, s).items() if c % P}
	def hv(a):
		ops[1] += 1
		return {n: c * (P + 1) // 2 % P for n, c in a.items()}
	ng = lambda a: ad({}, a, -1)
	sub = lambda v, w: (ad(v[0], w[0], -1), ad(v[1], w[1], -1))
	qt = lambda x0, x1: (hv(hv(ad(ng(x0), x1, -1))), hv(hv(ad(x0, x1, -1))))
	m0, m1 = M.fp6_vars("m", 0), M.fp6_vars("m", 6)
	vm1 = [(ad(*m1[2], -1), ad(*m1[2])), m1[0], m1[1]]
	out = []
	for (a0, a1), (b0, b1), (c0, c1) in (m0, [sub(v, w) for v, w in zip(vm1, m0)], [sub(v, w) for v, w in zip(m1, m0)]):
		e0, e1 = ad(a0, c0), ad(a1, c1)
		Bt = [(ng(a0), ng(a1)), qt(ad(e0, b0), ad(e1, b1)), qt(ad(e0, b0, -1), ad(e1, b1, -1)),
			 (hv(ad(ad(c0, a0, -1), b1, -1)), hv(ad(ad(b0, c1), a1, -1))), (ng(c0), ng(c1))]
		out.append([(B1, ad(B0, B1), ad(B0, B1, -1)) for B0, B1 in Bt])
	return out, ops

def mul_plan():
	"""Fp12 x m: k1 = (X0+X1) M0, k2 = X1 (v M1 - M0), k3 = X0 (M1 - M0), c0 = k1 + k2, c1 = k1 + k3. Product
  q = 9t + 3i + j is instance i's y(t) side (A0 + A1, A0, A1) against entry vector q (B1, B0 + B1, B0 - B1), so
  P_t = (P1 - P0, P2 + P0) = y(t) B_t. c0, c1 are summed at the points: W_t = P_t(k1) + P_t(k2 or k3), output
  (c, o, k) = part k of sum_t COEF[o][t] W_t."""
	x0, x1 = M.fp6_vars("x", 0), M.fp6_vars("x", 6)
	ys = [M.fp6_lsum(x0, x1), x1, x0]
	ent, ops = entry_program()
	prods = []
	for t in range(5):
		for i in range(3):
			f0, f1 = at(ys[i], t)
			prods += zip((M.ladd(f0, f1), f0, f1), ent[i][t])
	wterms, outs = [], []
	for c in range(2):
		for o in range(3):
			for k in range(2):
				w = {}
				for t, (a, b) in enumerate(COEF[o]):
					for kk, s in ((0, a), (1, -b)) if k == 0 else ((1, a), (0, b)):
						if s:
							w[t, kk] = s
				d = {}
				for (t, kk), s in w.items():
					for q in (9 * t, 9 * t + 3 + 3 * c):
						for qq, v in ((q + 1, 1), (q, -1)) if kk == 0 else ((q + 2, 1), (q, 1)):
							d[qq] = d.get(qq, 0) + s * v
				wterms.append(w)
				outs.append({q: v for q, v in d.items() if v})
	return prods, outs, wterms, ops, ys

def mul_kernel():
	a = Asm()
	prods, outs, wterms, ops, ys = mul_plan()
	assert len(prods) == 45 and ops == [104, 30]  # entry build: 104 modular add/sub/neg and 30 halvings per table value
	mod = lambda q: {k: v % P for k, v in q.items() if v % P}
	fprods, fouts = M.mul_model()             # the reference Karatsuba model, coefficient by coefficient mod p
	assert all(mod(M.expand(prods, c)) == mod(M.expand(fprods, f)) for c, f in zip(outs, fouts))
	vx = {f"x{i}": P - 1 for i in range(12)}
	lone = lambda L: list(L.values()) == [1] and "1" not in L
	xf = {}                                   # instances 1, 2: y(t) + k p, the least k making it non-negative
	for i in (1, 2):
		for t in range(5):
			for k, L in enumerate(at(ys[i], t)):
				lo = M.lin_range(L, vx)[0]
				xf[i, t, k] = M.ladd(L, {"1": -(lo // P) * P}) if lo < 0 else L
	wt = [dict(w) for w in wterms]
	for w, raw in zip(wt, wterms):
		for s, (x, y) in enumerate(PAIRS):
			if x in w and y in w and w[x] * w[y] == -1:
				w["S", s] = w.pop(x)
				w.pop(y)
		back = {}
		for key, s in w.items():
			for kk, v in ((PAIRS[key[1]][0], 1), (PAIRS[key[1]][1], -1)) if key[0] == "S" else ((key, 1),):
				back[kk] = back.get(kk, 0) + s * v
		assert {k: v for k, v in back.items() if v} == raw and all(abs(s) == 1 for s in w.values())
	# frame: tmp and the instance-0 point product, then per point its stored x-side forms, W slots and pair sums
	tmp, atom, off = SAVE, SAVE + 96, SAVE + 288
	st, wslot, sslot = {}, {}, {}
	for t in range(5):
		for i in (1, 2, 0):
			for k in range(2):
				if i == 0 or not lone(xf[i, t, k]):
					st[i, t, k], off = off, off + 48
		for c in range(2):
			for k in range(2):
				wslot[c, t, k], off = off, off + 96
				for s, (x, y) in enumerate(PAIRS):
					if x == (t, k):
						sslot[c, s], off = off, off + 96
	frame = (off + 15) // 16 * 16
	prologue(a, frame)

	def operand(i, t, k):
		if (i, t, k) in st:
			return "sp", st[i, t, k]
		(n, _), = xf[i, t, k].items()
		return "x1", 48 * int(n[1:])

	def build_form(L):
		"""A = L exactly: inputs with coefficients +-1, then + k p from Z[6:]."""
		terms = sorted(((n, c) for n, c in L.items() if n != "1"), key=lambda nc: -nc[1])
		kp = L.get("1", 0) // P
		assert L.get("1", 0) == kp * P and terms[0][1] == 1 and all(abs(c) == 1 for _, c in terms)
		for j, (n, c) in enumerate(terms):
			a.load6(T if j else A, "x1", 48 * int(n[1:]))
			if j:
				a.chain("addsub", c, A, A, T)
		for _ in range(kp):
			a.chain("addsub", 1, A, A, Z[6:])

	def reduce_to_p(hi, final):
		"""A < hi < 8p: subtract q p (q = umulh(A5, QC), the q*p table) to [0, 2p); then p unless that borrows.
    Returns the bound the emitted code leaves, which the range check reads."""
		if hi >= 2 * P:
			a.ins("xred", f"umulh {B[1]}, {A[5]}, {B[0]}")
			a.ins("xred", f"add {B[1]}, x3, {B[1]}, lsl #6")
			a.load6(T, B[1], POOL["qt"])
			a.chain("xred", -1, A, A, T)
		if not final:
			return 2 * P - 1
		a.chain("xred", -1, T, A, Z[6:])
		for j in range(6):
			a.ins("xred", f"csel {A[j]}, {A[j]}, {T[j]}, lo")
		return P - 1

	a.set_pool(["x0"])                        # x-side forms; p in Z[6:], QC in B[0]
	a.load6(Z[6:], "x3", POOL["p"])
	a.ins("mem", f"ldr {B[0]}, [x3, #{POOL['c']}]")
	for i in (1, 2):
		for k in range(2):                      # e = y0 + y2 once: y(1) = e + y1, y(-1) = e - y1 + p
			(n0, _), (n1, _), (n2, _) = (list(ys[i][j][k].items())[0] for j in range(3))
			e = M.lin(n0, n2)
			assert xf[i, 1, k] == M.ladd(e, M.lin(n1)) and xf[i, 2, k] == M.ladd(M.ladd(e, M.lin(n1), -1), {"1": P})
			a.load6(Z[:6], "x1", 48 * int(n0[1:]))
			a.load6(T, "x1", 48 * int(n2[1:]))
			a.chain("addsub", 1, Z[:6], Z[:6], T)
			a.load6(T, "x1", 48 * int(n1[1:]))
			a.chain("addsub", 1, A, Z[:6], T)
			a.store6(A, "sp", st[i, 1, k])
			a.chain("addsub", -1, A, Z[:6], T)
			a.chain("addsub", 1, A, A, Z[6:])
			a.store6(A, "sp", st[i, 2, k])
			for t in (0, 3, 4):
				if (i, t, k) in st:
					build_form(xf[i, t, k])
					a.store6(A, "sp", st[i, t, k])
	vb = {}                                   # instance 0: y1(t) + y2(t) reduced, the bound the code leaves
	for t in range(5):
		for k in range(2):
			a.load6(A, *operand(1, t, k))
			a.load6(T, *operand(2, t, k))
			a.chain("addsub", 1, A, A, T)
			hi = M.lin_range(M.ladd(xf[1, t, k], xf[2, t, k]), vx)[1]
			assert hi < 8 * P
			vb[f"r{2 * t + k}"] = reduce_to_p(hi, (t, k) not in PARTIAL)
			a.store6(A, "sp", st[0, t, k])
	# ranges: every stored x-side form is one variable bounded as emitted, entry values reduced, inputs in [0, p)
	vm = dict(vx, **vb, **{f"e{q}": P - 1 for q in range(45)})
	def var(i, t, k):
		if i == 0:
			return {f"r{2 * t + k}": 1}
		L = xf[i, t, k]
		if lone(L):
			return dict(L)
		lo, vm[f"f{i}_{t}_{k}"] = M.lin_range(L, vx)
		assert lo >= 0
		return {f"f{i}_{t}_{k}": 1}
	rx = []
	for q in range(45):
		t, i, j = q // 9, q // 3 % 3, q % 3
		f0, f1 = var(i, t, 0), var(i, t, 1)
		rx.append((M.ladd(f0, f1), f0, f1)[j])
	rows = []
	for k, oc in enumerate(outs):
		lo, hi = M.poly_range(M.expand([(xs, {f"e{q}": 1}) for q, xs in enumerate(rx)], oc), vm)
		K = -(lo // (P * R)) * P * R if lo < 0 else 0
		bres = (K + hi) // R + 1 + P            # r < T/R + p
		assert K + lo >= 0 and K + hi < (1 << 768) and bres <= (1 << 384) and K <= QMAX * P * R
		rows.append((k, lo, hi, K, bres))
	opmax = max(M.lin_range(xs, vm)[1] for xs in rx)
	assert opmax <= 6 * P and max(r[4] for r in rows) * 10000 < 96411 * P < (1 << 384) * 10000
	a.lines[:0] = ["// " + s for s in (
		"BLS12-381 Fp12 multiplication by a fixed table entry, Montgomery form, outputs in [0, p)\n"
		"generated by gen.py; args x0 out, x1 x, x2 entry, x3 constant pool\n"
		"Toom-3 over Fp6 at t = 0, 1, -1, -u, inf: 45 Fp products, c0 and c1 summed at the points\n"
		f"x-side operands at most {opmax / P:.2f}p < 2^384; entry values reduced mod p\n"
		f"entry: {len(prods)} vectors x 48 bytes = {48 * len(prods)} bytes per table entry\n"
		+ bound_table(rows, "2^384 < 10p")).split("\n")]
	a.cats[:0] = [None] * (len(a.lines) - len(a.cats))

	a.set_pool(["x0", "x3"])                  # products, summed at the points

	def prod(i, t, j):
		terms = [(1,) + operand(i, t, 0), (1,) + operand(i, t, 1)] if j == 0 else [(1,) + operand(i, t, j - 1)]
		product(a, terms, "x2", 48 * (9 * t + 3 * i + j))

	def store12(off):
		a.store6(Z[:6], "sp", off)
		a.store6(Z[6:], "sp", off + 48)

	for t in range(5):
		prod(0, t, 0)
		store12(tmp)
		for k in range(2):
			prod(0, t, 1 + k)
			a.acc12(1 if k else -1, "sp", tmp)
			store12(atom + 96 * k)
		for c in range(2):
			prod(1 + c, t, 0)
			store12(tmp)
			for k in range(2):
				prod(1 + c, t, 1 + k)
				a.acc12(1 if k else -1, "sp", tmp)
				a.acc12(1, "sp", atom + 96 * k)
				store12(wslot[c, t, k])
				fresh = True
				for s, (x, y) in enumerate(PAIRS):
					if x == (t, k):
						if not fresh:
							a.load6(Z[:6], "sp", wslot[c, t, k])
							a.load6(Z[6:], "sp", wslot[c, t, k] + 48)
						a.acc12(-1, "sp", wslot[(c,) + y])
						store12(sslot[c, s])
						fresh = False
	a.ins("mem", "ldp x0, x3, [sp, #96]")    # outputs: one lazy sum per coordinate, K = j p 2^384, reduction
	a.load6(PR, "x3", POOL["p"])
	a.ins("mem", f"ldr {N0R}, [x3, #{POOL['n0']}]")
	a.set_pool(["x1", "x2"])
	for k, w in enumerate(wt):
		c = k // 6
		terms = sorted(w.items(), key=lambda kv: -kv[1])
		offs = [sslot[c, key[1]] if key[0] == "S" else wslot[(c,) + key] for key, _ in terms]
		a.load6(Z[:6], "sp", offs[0])
		a.load6(Z[6:], "sp", offs[0] + 48)
		if terms[0][1] < 0:                     # every term negative: negate the first
			a.ins("addsub", f"subs {Z[0]}, xzr, {Z[0]}")
			for j in range(1, 12):
				a.ins("addsub", f"sbcs {Z[j]}, xzr, {Z[j]}")
		for (_, s), o in zip(terms[1:], offs[1:]):
			a.acc12(s, "sp", o)
		if rows[k][3]:
			a.load6(T, "x3", 64 * (rows[k][3] // (P * R)))
			a.chain("addsub", 1, Z[6:], Z[6:], T)
		redc_store(a, k)
	epilogue(a, frame)
	return a, {"entry": [e for _, e in prods], "xs": [xs for xs, _ in prods], "outs": outs, "frame": frame,
				"K": {}}, rows

# -- cyclotomic squaring
def sqr_kernel():
	a = Asm()
	vmax = {f"a{i}": P - 1 for i in range(12)}
	prods, outs = M.sqr_model()
	polys = [M.expand(prods, c) for c in outs]
	rows = out_bounds(polys, vmax)

	# stored vectors: the Fp4 sums s = a + b (two inputs, coefficient 1, no offset) first, then every other
	# B operand sum (the Fp4 sums are inside B operands of the t2 squarings)
	keys = [tuple(sorted(L.items())) for Aop, Bop in prods for L in (Aop, Bop)
			if len(L) == 2 and "1" not in L and all(c == 1 for c in L.values())]
	keys += [tuple(sorted(Bop.items())) for _, Bop in prods if len(Bop) > 1]
	allv, scr = {}, SAVE
	for key in keys:
		if key not in allv:
			allv[key], scr = scr, scr + 48

	def decompose(L):
		L = dict(L)
		terms = []
		if "1" in L:
			c = L.pop("1")
			mult = c // P
			assert c == mult * P and mult in (1, 2, 4)
			terms.append((1, "x3", POOL["p" if mult == 1 else f"{mult}p"]))
		neg = []
		for key in sorted(allv, key=lambda k: (len(k), k), reverse=True):
			if len(key) == 2 and all(L.get(n, 0) >= c > 0 for n, c in key):
				for n, c in key:
					L[n] -= c
				terms.append((1, "sp", allv[key]))
			elif len(key) == 2 and all(L.get(n, 0) <= -c < 0 for n, c in key):
				for n, c in key:
					L[n] += c
				neg.append((-1, "sp", allv[key]))
		L = {n: c for n, c in L.items() if c}
		for n, c in sorted(L.items(), key=lambda kv: -kv[1]):
			for _ in range(abs(c)):
				(terms if c > 0 else neg).append((1 if c > 0 else -1, "x1", 48 * int(n[1:])))
		return terms + neg

	prod_off = scr
	scr += 96 * len(prods)
	frame = (scr + 15) // 16 * 16
	opmax = 0
	for Aop, Bop in prods:
		for L in (Aop, Bop):
			lo, hi = M.lin_range(L, vmax)
			assert lo >= 0 and hi < (1 << 384)
			opmax = max(opmax, hi)
	a.comment("BLS12-381 cyclotomic squaring (Granger-Scott, blst's cyclotomic_sqr_fp12 formula), Montgomery\n"
			 "form, outputs in [0, p); generated by gen.py; args x0 out, x1 x, x3 constant pool\n"
			 f"every product operand is an exact non-negative sum below {opmax / P:.2f}p < 2^384\n"
			 + bound_table(rows))
	prologue(a, frame)
	a.set_pool(["x0", "x2"])
	for key, off in list(allv.items()):          # stored sums (x3 still the constant pool here)
		L = dict(key)
		keep = allv.pop(key)
		terms = decompose(L)
		allv[key] = keep
		load_lin(a, A, terms)
		a.store6(A, "sp", off)
	a.set_pool(["x0", "x2", "x3"])
	for i, (Aop, Bop) in enumerate(prods):
		bkey = tuple(sorted(Bop.items()))
		if bkey in allv:
			bb, bo = "sp", allv[bkey]
		else:
			(n, c), = Bop.items()
			assert c == 1
			bb, bo = "x1", 48 * int(n[1:])
		double = all(v == 2 for v in Aop.values())
		aterms = decompose({n: 1 for n in Aop}) if double else decompose(Aop)
		product(a, aterms, bb, bo, double_a=double)
		a.store6(Z[:6], "sp", prod_off + 96 * i)
		a.store6(Z[6:], "sp", prod_off + 96 * i + 48)
	a.ins("mem", "ldp x0, x3, [sp, #96]")
	a.load6(PR, "x3", POOL["p"])
	a.ins("mem", f"ldr {N0R}, [x3, #{POOL['n0']}]")
	a.set_pool(["x2"])
	consts = {}
	for k, c in enumerate(outs):
		prodterms = sorted((pi, v // 3) for pi, v in c.items() if not isinstance(pi, tuple))
		assert all(v * 3 == c[pi] for pi, v in prodterms)
		hi_terms = [(key[1], v) for key, v in c.items() if isinstance(key, tuple)]
		first = next(i for i, (pi, s) in enumerate(prodterms) if s > 0)
		pi, s = prodterms.pop(first)
		a.load6(Z[:6], "sp", prod_off + 96 * pi)
		a.load6(Z[6:], "sp", prod_off + 96 * pi + 48)
		for pi, s in prodterms:
			for _ in range(abs(s)):
				a.acc12(1 if s > 0 else -1, "sp", prod_off + 96 * pi)
		a.ins("addsub", f"add {T[0]}, {Z[0]}, {Z[0]}")   # 3Z = Z + 2Z, 2Z by extr (no flags)
		for i in range(12):
			if i < 11:
				a.ins("addsub", f"extr {T[1 - i % 2]}, {Z[i + 1]}, {Z[i]}, #63")
			a.ins("addsub", f"{'adcs' if i else 'adds'} {Z[i]}, {Z[i]}, {T[i % 2]}")
		for name, v in hi_terms:                   # +-2 a in the upper half
			a.load6(T, "x1", 48 * int(name[1:]))
			for _ in range(abs(v)):
				a.chain("addsub", 1 if v > 0 else -1, Z[6:], Z[6:], T)
		K = rows[k][3]
		if K:
			consts[k] = K
			a.acc12(1, "x3", POOL["K"] + 96 * k)
		redc_store(a, k)
	epilogue(a, frame)
	return a, {"K": {str(k): v for k, v in consts.items()}, "frame": frame}, rows

def main():
	here = os.path.dirname(os.path.abspath(__file__))
	rep = {}
	for name, fn in (("mul", mul_kernel), ("sqr", sqr_kernel)):
		a, layout, rows = fn()
		with open(os.path.join(here, f"{name}.s"), "w") as f:
			f.write(a.text())
		layout["count"] = a.count()
		layout["executed"] = sum(layout["count"].values())
		layout["pool"] = {"n0": N0, "qc": QC, "offsets": POOL}
		with open(os.path.join(here, f"{name}.layout.json"), "w") as f:
			json.dump(layout, f, indent=1, default=str)
		rep[name] = (layout["executed"], layout["count"])
		print(name, layout["executed"], layout["count"])
	return rep

if __name__ == "__main__":
	main()
