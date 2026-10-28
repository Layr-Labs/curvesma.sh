"""Run lazy.s and root.s (gen_lazy.py) in emu.py against a big-integer reference on real BLS12-381 GT
elements: lazy.s must return exactly the stored a, b, B, a zeta, B zeta of y' = x * Pi, the minimum K of the twelve
stored keys and kappa, which must select the minimum-key conjugate among the twelve explicit Frobenius conjugates;
root.s the same for y' = x and its per-base data (six pre-sums, Xa, Xd), equal to gen_lazy.block_vectors (walk state
pre-filled with random words). Every call must execute the same count. Besides the sampled cases, a fixed set runs every stored input at p - 1 or 0 (the operand and accumulator
extremes: S0 = 0 with S1 and S2 at their maximum, S0 at its maximum with S1 and S2 near 0) and keys placed by
solving for node coordinates: a = 0, b = 0, d = a - b = 0, exact ties, and key pairs that agree on every limb but
one. With the argument "mutants", each deliberately broken listing (a line edit, or a generator constant changed
for one kernel only) must fail.
Usage: python3 -I check_lazy.py <n_cases> <seed> [mutants]"""
import json, os, random, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import emu
import check as C
import gen_lazy as L

P = C.P
RM = C.RM
RINV = pow(RM, -1, P)
OUT, XB, EB, PB, SP = 0x10000, 0x20000, 0x30000, 0x40000, 0x80000

def wflat(a):
	"""blst flat vector -> w^m coefficients (Fp2 pairs), m = 0..5."""
	v = C.flat(a)
	return [(v[L.fidx(m, 0)], v[L.fidx(m, 1)]) for m in range(6)]

def frob_keys_direct(a):
	"""Re(w^1 coefficient) of a^(p^k), k = 0..11, by explicit Frobenius maps (the definition)."""
	out = []
	y = a
	for k in range(12):
		out.append(wflat(y)[1][0])
		y = C.f12pow(y, P)
	return out

def gamma(m):
	"""Node entry of Pi = m: gamma_0 = m_1, gamma_1 = m_0, gamma_k = xi m_{7-k} (k = 2..5), via gen_lazy.node_vector."""
	w = wflat(m)
	xi = lambda c: ((c[0] - c[1]) % P, (c[0] + c[1]) % P)
	g = [w[1], w[0]] + [xi(w[7 - k]) for k in range(2, 6)]
	return L.node_vector([(re * RM % P, im * RM % P) for re, im in g])

def from_stored(v):
	"""The element whose twelve stored (Montgomery) flat coordinates are v."""
	return C.unflat([t * RINV % P for t in v])

def node_from_stored(pairs):
	"""The node element m whose gamma_k has stored (re, im) = pairs[k]."""
	g = [(re * RINV % P, im * RINV % P) for re, im in pairs]
	half = pow(2, -1, P)
	w = [g[1], g[0], None, None, None, None]
	for k in range(2, 6):          # xi (c0, c1) = (c0 - c1, c0 + c1) = g_k
		t0, t1 = g[k]
		w[7 - k] = ((t0 + t1) * half % P, (t1 - t0) * half % P)
	v = [0] * 12
	for m in range(6):
		v[L.fidx(m, 0)], v[L.fidx(m, 1)] = w[m]
	return C.unflat(v)

FUNC = {"a": lambda a, b: a, "b": lambda a, b: b, "d": lambda a, b: (a - b) % P}

def place(x, m, conds):
	"""m adjusted on the first len(conds) flat coordinates (that give an invertible system) so that y = x * m
    has the plain w^1 functionals of conds (name, plain value) exactly."""
	f = lambda e: [FUNC[n](*wflat(C.f12mul(x, e))[1]) for n, _ in conds]
	base, mv = f(m), C.flat(m)
	unit = lambda j: C.unflat([int(i == j) for i in range(12)])
	cols = [f(unit(j)) for j in range(12)]
	for js in ([0], [1], [2], [3]) if len(conds) == 1 else ([0, 1], [0, 2], [1, 2], [2, 3]):
		A = [[cols[j][i] for j in js] for i in range(len(conds))]
		r = [(t - base[i]) % P for i, (_, t) in enumerate(conds)]
		det = A[0][0] if len(js) == 1 else (A[0][0] * A[1][1] - A[0][1] * A[1][0]) % P
		if det:
			if len(js) == 1:
				lam = [r[0] * pow(det, -1, P) % P]
			else:
				inv = pow(det, -1, P)
				lam = [(r[0] * A[1][1] - r[1] * A[0][1]) * inv % P, (A[0][0] * r[1] - A[1][0] * r[0]) * inv % P]
			for j, l in zip(js, lam):
				mv[j] = (mv[j] + l) % P
			return C.unflat(mv)
	raise AssertionError("no invertible coordinate pair")

def run(ops, pool, xv, entry, rng):
	m = emu.Machine()
	for i in range(31):
		m.r[i] = rng.getrandbits(64)
	saved = list(m.r[19:31])
	m.mem.update(pool)
	xm = [(v * RM) % P for v in C.flat(xv)]
	C.put(m.mem, XB, xm)
	if entry is not None:
		C.put(m.mem, XB + L.XS, L.block_vectors(xm))
		C.put(m.mem, EB, entry)
	else:
		C.put(m.mem, XB + L.XS, [rng.getrandbits(384) for _ in range(10)])
	m.r[0], m.r[1], m.r[2], m.r[3], m.r[31] = OUT, XB, EB, PB, SP
	n = emu.run(ops, m)
	assert m.r[31] == SP and m.r[19:31] == saved, "callee-saved state or sp not restored"
	rd = lambda off: sum(m.mem[OUT + off + 8 * j] << (64 * j) for j in range(6))
	res = {k: rd(L.OUT[k]) for k in ("a", "b", "B", "az", "Bz", "K")}
	res["kap"] = m.mem[OUT + L.OUT["kap"]]
	xs = [sum(m.mem[XB + L.XS + 48 * k + 8 * j] << (64 * j) for j in range(6)) for k in range(10)]
	return n, res, xs, xm

def reference(y):
	w = wflat(y)
	a, b = w[1]
	d1 = L.D1
	B = (a * d1[0] + b * d1[1]) % P
	st = lambda v: v * RM % P
	keys_int = frob_keys_direct(y)
	keys = [st(v) for v in keys_int]
	K = min(keys)
	return {"a": st(a), "b": st(b), "B": st(B), "az": st(a * L.ZETA % P), "Bz": st(B * L.ZETA % P),
			"K": K, "kap_set": [k for k in range(12) if keys[k] == K]}

def fixed_cases(rng, gens):
	"""(kind, x, m) with the operand extremes and placed keys; root.s runs on x, lazy.s on x * m."""
	top, zero = [P - 1] * 12, [0] * 12
	out = [("max", from_stored(top), node_from_stored([(P - 1, P - 1)] * 6)),
		 ("S0_zero_S1_S2_max", from_stored([0, P - 1] * 6), node_from_stored([(0, P - 1)] * 6)),
		 ("S0_max_S1_S2_low", from_stored([P - 1, 1] * 6), node_from_stored([(P - 1, 1)] * 6)),
		 ("S0_S2_max_S1_zero", from_stored([P - 1, 0] * 6), node_from_stored([(P - 1, 0)] * 6)),
		 ("zero_x", from_stored(zero), gens[1]), ("zero_node", gens[0], node_from_stored([(0, 0)] * 6))]
	# S0 = 0 < S1 with (S1 / p mod 2^384) p < S1: the offsets in Xa and Ea keep a_T >= 0
	R, n = 1 << 384, 0
	while True:
		r2 = random.Random(n)
		v, h = [r2.randrange(P) for _ in range(6)], [r2.randrange(P) for _ in range(6)]
		s1 = sum(i * j for i, j in zip(v, h))
		if s1 * pow(P, -1, R) % R * P < s1:
			break
		n += 1
	xv = [0] * 12
	for k in range(6):
		xv[L.fidx(k, 1)] = v[k]
	out.append(("a_T_neg", from_stored(xv), node_from_stored([(0, j) for j in h])))
	pl = lambda v: v * RINV % P
	dpl = lambda Bst: pl(Bst) * pow(L.D1[0], -1, P) % P
	X = 3 << 64 | 7
	placed = {"a0": [("a", 0)], "b0": [("b", 0)], "d0": [("d", 0)], "a_p-1": [("a", pl(P - 1))],
			 "tie": [("a", pl(X)), ("d", dpl(X))], "B_low_limb": [("a", pl(X)), ("d", dpl(X - 5))],
			 "B_up_limb1": [("a", pl(X)), ("d", dpl(X + (1 << 64)))],
			 "a_low_limb": [("a", pl(X - 5)), ("d", dpl(X))], "B_hi": [("a", pl(X)), ("d", dpl(X + (1 << 320)))]}
	for name, conds in placed.items():
		x, m = gens[0], C.rand12(rng)
		m = place(x, m, conds)
		out.append(("lazy:" + name, x, m))
		xv = C.flat(C.rand12(rng))
		y = C.f12mul(x, m)
		xv[6], xv[7] = C.flat(y)[6], C.flat(y)[7]
		out.append(("root:" + name, C.unflat(xv), m))
	return out

def check_all(progs, pool, n_cases, rng, gens, stop_at_first=False):
	counts = {"lazy": set(), "root": set()}
	fails = {"lazy": 0, "root": 0}
	kinds = {}
	x, m = gens[0], gens[1]
	cases = []
	for i in range(n_cases):
		kind = ("gt", "gt", "rand", "edge")[i % 4]
		if kind == "gt":
			x = C.f12mul(x, gens[2])
			m = C.f12mul(m, gens[1])
		elif kind == "rand":
			x, m = C.rand12(rng), C.rand12(rng)
		else:
			x = C.unflat([rng.choice((0, 1, P - 1, P - 2, rng.randrange(P))) for _ in range(12)])
			m = C.unflat([rng.choice((0, 1, P - 1, P - 2, rng.randrange(P))) for _ in range(12)])
		cases.append((kind, x, m))
	cases = fixed_cases(rng, gens) + cases
	for kind, x, m in cases:
		kinds[kind] = kinds.get(kind, 0) + 1
		for name in ("lazy", "root"):
			if kind.startswith("lazy:") and name == "root" or kind.startswith("root:") and name == "lazy":
				continue
			try:
				if name == "lazy":
					n, got, _, _ = run(progs[name], pool, x, gamma(m), rng)
					y = C.f12mul(x, m)
				else:
					n, got, xs, xm = run(progs[name], pool, x, None, rng)
					y = x
					if xs != L.block_vectors(xm):
						fails[name] += 1
			except (KeyError, ValueError, AssertionError):
				fails[name] += 1
				if stop_at_first:
					return fails, counts, kinds
				continue
			ref = reference(y)
			counts[name].add(n)
			ok = all(got[k] == ref[k] for k in ("a", "b", "B", "az", "Bz", "K")) and got["kap"] in ref["kap_set"]
			if not ok:
				fails[name] += 1
				if stop_at_first:
					return fails, counts, kinds
				if fails[name] <= 3:
					print("MISMATCH", kind, name, {k: (got[k], ref.get(k)) for k in got if got[k] != ref.get(k)})
	return fails, counts, kinds

# Deliberately broken listings: (name, kernels, line transform). Each must fail check_all.
CMP = re.compile(r"sbcs xzr, (x2[5-9]|x30), (x[4-9]|x1[0-5])$")
MUTANTS = [
	("min_top_limb_only", ("lazy", "root"), lambda s: CMP.sub(r"subs xzr, \1, \2", s)),
	("min_z_flag", ("lazy", "root"), lambda s: re.sub(r"(csel x(?:[4-9]|1[0-5]), x(?:2[5-9]|30), x(?:[4-9]|1[0-5])), lo$",
													 r"\1, ls", s)),
	("min_dropped_csel", ("lazy", "root"), "drop:csel x4, x25, x4, lo"),
	("b_select_flipped", ("lazy", "root"), lambda s: re.sub(r"(csel (x1[0-5]), x(?:2[5-9]|30), \2), hs$", r"\1, lo", s)),
	("presum_unreduced", ("root",), lambda s: re.sub(r"csel (x1[67]|x19|x2[0-2]), x(?:2[5-9]|30), \1, hs$",
													r"csel \1, \1, \1, hs", s)),
	("operand_carry_dropped", ("lazy",), lambda s: re.sub(r"adcs (x2[01]), \1, (x2[89])$", r"add \1, \1, \2", s)),
	("d_not_doubled", ("lazy",), "drop:adds x4, x4, x4"),
	("xa_read_as_xd", ("lazy",), {"XA_OFF": L.XD_OFF}),
	("ea_read_as_ed", ("lazy",), {"EA_OFF": L.ED_OFF}),
	("pairs_regrouped", ("lazy",), {"PAIRS": ((0, 2), (1, 3), (4, 5))}),
]

def mutate(text, how, kern):
	if isinstance(how, dict):
		keep = {k: getattr(L, k) for k in how}
		try:
			for k, v in how.items():
				setattr(L, k, v)
			return getattr(L, kern + "_kernel")()[0].text()
		finally:
			for k, v in keep.items():
				setattr(L, k, v)
	lines = text.split("\n")
	if callable(how):
		return "\n".join(how(l) for l in lines)
	pat = how.split(":", 1)[1]
	del lines[[l.strip() for l in lines].index(pat)]
	return "\n".join(lines)

def main():
	n_cases, seed = int(sys.argv[1]), int(sys.argv[2])
	rng = random.Random(seed)
	pool = L.pool_words()
	pool = {PB + k: v for k, v in pool.items()}
	text = {name: open(os.path.join(HERE, f"{name}.s")).read() for name in ("lazy", "root")}
	progs = {name: emu.compile_prog(emu.parse(t)) for name, t in text.items()}
	e = (P ** 12 - 1) // C.Rr
	gens = [C.f12pow(C.rand12(rng), e) for _ in range(3)]
	fails, counts, kinds = check_all(progs, pool, n_cases, rng, gens)
	res = {"cases": kinds, "executed_counts": {k: sorted(v) for k, v in counts.items()}, "fails": fails}
	if sys.argv[3:] == ["mutants"]:
		caught = {}
		for name, kerns, how in MUTANTS:
			bad = dict(progs)
			for k in kerns:
				t = mutate(text[k], how, k)
				assert t != text[k], name
				bad[k] = emu.compile_prog(emu.parse(t))
			f, _, _ = check_all(bad, pool, 4, random.Random(seed + 1), gens, stop_at_first=True)
			caught[name] = sum(f.values()) > 0
		res["mutants_caught"] = caught
		assert all(caught.values()), caught
	print(json.dumps(res))
	assert all(len(v) == 1 for v in counts.values()) and not any(fails.values())
	print("PASS")

if __name__ == "__main__":
	main()
