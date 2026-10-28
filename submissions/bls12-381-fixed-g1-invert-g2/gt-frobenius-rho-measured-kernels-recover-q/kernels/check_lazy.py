"""Run lazy.s, root.s and node.s (gen_lazy.py) in emu.py against big-integer references. root.s on any x, and lazy.s
for x in GT and any node Pi, must return exactly the stored a, b, B, a zeta, B zeta of y' = x * Pi, the minimum K of
the twelve stored keys and kappa, which must select the minimum-key conjugate among the twelve explicit Frobenius
conjugates; root.s must write gen_lazy.block_vectors (walk state pre-filled with random words). lazy.s uses W(x) = -x3, true on GT only, so on every case it is checked against the integer
function it computes, REDC of S0 - S1 + Ea 2^384 and 2 S0 - S2 + Ed 2^384 (ident), and on GT ident against the
definition; node.s must produce gen_lazy.entry. Fixed cases: raw stored inputs at 0 and p - 1 with Ea, Ed at both
ends of [2p, 3p) and the offset's own extreme (a_T_neg), root.s at x0 = 0, x4 extremes and all coordinates 0 or
p - 1, keys placed by solving for node coordinates (a = 0, b = 0, d = 0, exact ties, keys equal on all limbs but one).
Every call must execute the same count. With "mutants", each broken listing or generator constant must fail.
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

def run(ops, pool, xm, bd, entry, rng):
	m = emu.Machine()
	for i in range(31):
		m.r[i] = rng.getrandbits(64)
	saved = list(m.r[19:31])
	m.mem.update(pool)
	C.put(m.mem, XB, xm)
	C.put(m.mem, XB + L.XS, bd if entry is not None else [rng.getrandbits(384) for _ in range(10)])
	if entry is not None:
		C.put(m.mem, EB, entry)
	m.r[0], m.r[1], m.r[2], m.r[3], m.r[31] = OUT, XB, EB, PB, SP
	n = emu.run(ops, m)
	assert m.r[31] == SP and m.r[19:31] == saved, "callee-saved state or sp not restored"
	rd = lambda off: sum(m.mem[OUT + off + 8 * j] << (64 * j) for j in range(6))
	res = {k: rd(L.OUT[k]) for k in ("a", "b", "B", "az", "Bz", "K")}
	res["kap"] = m.mem[OUT + L.OUT["kap"]]
	xs = [sum(m.mem[XB + L.XS + 48 * k + 8 * j] << (64 * j) for j in range(6)) for k in range(10)]
	return n, res, xs

def node(ops, pool, vec, rng):
	"""node.s on a node entry holding the constants' (re, im) of vec (pre-sums and Ea, Ed pre-filled with random words):
  the entry must become gen_lazy.entry of those constants with their pre-sums mod p."""
	vec = [v if i % 3 < 2 else (vec[i - 2] + vec[i - 1]) % P for i, v in enumerate(vec)]
	m = emu.Machine()
	m.mem.update(pool)
	C.put(m.mem, EB, [v if i % 3 < 2 else rng.getrandbits(384) for i, v in enumerate(vec)] + [rng.getrandbits(384)] * 2)
	m.r[0], m.r[1], m.r[3], m.r[31] = EB, EB, PB, SP
	n = emu.run(ops, m)
	got = [sum(m.mem[EB + 48 * i + 8 * j] << (64 * j) for j in range(6)) for i in range(20)]
	assert m.r[31] == SP and got == L.entry(vec), "node.s"
	return n

def outputs(a, b):
	"""Stored a, b -> the stored outputs and the twelve keys in Tower.keys order (all Fp-linear in a, b)."""
	B = (a - b) * L.D1[0] % P
	az, Bz = a * L.ZETA % P, B * L.ZETA % P
	keys = []
	for j, f in enumerate(((lambda v, z: v), (lambda v, z: z), (lambda v, z: z - v), (lambda v, z: -v),
												 (lambda v, z: -z), (lambda v, z: v - z))):
		keys += [f(a, az) % P, f(B, Bz) % P]
	K = min(keys)
	return {"a": a, "b": b, "B": B, "az": az, "Bz": Bz, "K": K, "kap_set": [k for k in range(12) if keys[k] == K]}

def ident(xm, bd, ent):
	"""lazy.s as an integer function of its stored inputs (any x): stored a, d = REDC(S0 - S1 + Ea 2^384),
  REDC(2 S0 - S2 + Ed 2^384), S_s = sum_i (operand 2i + node 2i)(operand 2i + 1 + node 2i + 1) on component s."""
	def opv(n, s):
		k, v = L.OPND[n]
		return bd[n] if s == 2 else (xm[L.fidx(v, s)] if k == "x" else bd[6 + 2 * v + s])
	vals = [opv(n, s) for n in range(6) for s in range(3)] + ent
	assert all(0 <= v < P for v in vals[:-2] + ent[:-2]) and all(L.EC[0] <= v <= L.EC[1] for v in ent[-2:])
	S = [sum((opv(2 * i, s) + ent[6 * i + s]) * (opv(2 * i + 1, s) + ent[6 * i + 3 + s]) for i in range(3)) for s in range(3)]
	aT, dT = S[0] - S[1] + ent[18] * L.R, 2 * S[0] - S[2] + ent[19] * L.R
	assert L.TA[0] <= aT <= L.TA[1] and L.TD[0] <= dT <= L.TD[1]
	a, d = aT * RINV % P, dT * RINV % P
	return outputs(a, (a - d) % P)

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

def raw(xv, bv, nv, off=None):
	"""Raw stored lazy inputs: x (12), per-base data (10) and node vectors (18) from per-component values
  (re, im, pre-sum) for every operand and constant; Ea, Ed = off, or entry()'s rule."""
	xm, bd = [0] * 12, [0] * 10
	for n, (k, v) in enumerate(L.OPND):
		bd[n] = xv[2]
		for s in (0, 1):
			if k == "x":
				xm[L.fidx(v, s)] = xv[s]
			else:
				bd[6 + 2 * v + s] = xv[s]
	vec = list(nv) * 6
	return xm, bd, (vec + [off, off] if off else L.entry(vec))

def fixed_cases(rng, gens):
	"""(kind, x, m) for root.s and lazy.s on x * m, and ("raw:...", xm, bd, ent) for lazy.s alone."""
	q, z = P - 1, 0
	out = [("raw:max", *raw((q, q, q), None, (q, q, q), L.EC[1])), ("raw:zero", *raw((z, z, z), None, (z, z, z), L.EC[0])),
		 ("raw:S0_zero_S1_S2_max", *raw((z, q, q), None, (z, q, q), L.EC[0])),
		 ("raw:S0_max_S1_S2_zero", *raw((q, z, z), None, (q, z, z), L.EC[1]))]
	# S0 = 0 < S1 with Ea (entry()'s rule) below S1 / 2^384 - p: a_T >= 0 only through the offset 2p
	n = 0
	while True:
		v = random.Random(n).randrange(P)
		case = raw((z, q, q), None, (z, v, v))
		if case[2][18] - L.EOFF + P < 3 * (q + v) ** 2 // L.R:
			break
		n += 1
	out.append(("raw:a_T_neg",) + case)
	X0 = [rng.randrange(P) for _ in range(12)]
	X0[0] = X0[1] = 0
	X0[L.fidx(4, 0)] = X0[L.fidx(4, 1)] = q
	X4 = list(X0)
	X4[L.fidx(4, 0)], X4[L.fidx(4, 1)] = 0, q
	out += [("root:" + k, from_stored(v), None) for k, v in (("max", [q] * 12), ("zero", [z] * 12), ("x0_zero", X0),
			("x4_edges", X4))]
	out.append(("zero_node", gens[0], C.unflat([0] * 12)))
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
	"""lazy.s against ident on every case and, for x in GT (gt, zero_node, lazy:*), ident against the definition;
  root.s against the definition and block_vectors."""
	counts = {"lazy": set(), "root": set(), "node": set()}
	fails = {"lazy": 0, "root": 0, "ident_vs_definition": 0}
	kinds = {}
	gx, gm = gens[0], gens[1]
	cases = []
	for i in range(n_cases):
		kind = ("gt", "gt", "rand", "edge")[i % 4]
		if kind == "gt":
			gx, gm = C.f12mul(gx, gens[2]), C.f12mul(gm, gens[1])
			x, m = gx, gm
		elif kind == "rand":
			x, m = C.rand12(rng), C.rand12(rng)
		else:
			x = C.unflat([rng.choice((0, 1, P - 1, P - 2, rng.randrange(P))) for _ in range(12)])
			m = C.unflat([rng.choice((0, 1, P - 1, P - 2, rng.randrange(P))) for _ in range(12)])
		cases.append((kind, x, m))
	for case in fixed_cases(rng, gens) + cases:
		kind = case[0]
		kinds[kind] = kinds.get(kind, 0) + 1
		for name in ("lazy", "root"):
			if kind.startswith(("lazy:", "raw:")) and name == "root" or kind.startswith("root:") and name == "lazy":
				continue
			try:
				xm = case[1] if kind.startswith("raw:") else [v * RM % P for v in C.flat(case[1])]
				if name == "lazy":
					bd, ent = (case[2], case[3]) if kind.startswith("raw:") else (L.block_vectors(xm), gamma(case[2]))
					ref = ident(xm, bd, ent)
					if kind in ("gt", "zero_node") or kind.startswith("lazy:"):
						d = reference(C.f12mul(case[1], case[2]))
						if any(ref[k] != d[k] for k in ("a", "b", "B", "az", "Bz", "K", "kap_set")):
							fails["ident_vs_definition"] += 1
					n, got, _ = run(progs[name], pool, xm, bd, ent, rng)
					counts["node"].add(node(progs["node"], pool, ent[:18], rng))
				elif name == "root":
					n, got, xs = run(progs[name], pool, xm, None, None, rng)
					ref = reference(case[1])
					if xs != L.block_vectors(xm):
						fails[name] += 1
			except (KeyError, ValueError, AssertionError):
				fails[name] += 1
				if stop_at_first:
					return fails, counts, kinds
				continue
			counts[name].add(n)
			ok = all(got[k] == ref[k] for k in ("a", "b", "B", "az", "Bz", "K")) and got["kap"] in ref["kap_set"]
			if not ok:
				fails[name] += 1
				if stop_at_first:
					return fails, counts, kinds
				if fails[name] <= 3:
					print("MISMATCH", kind, name, {k: (got[k], ref.get(k)) for k in got if got[k] != ref.get(k)})
	return fails, counts, kinds

def u0_unreduced():
	"""root.s's negation in U0 keeps T + p (p, not 0, for x0 = 0): its six csel take the uncorrected sum."""
	state = {"left": 0}
	def f(s):
		if s.strip() == "subs x25, xzr, x4":
			state["left"] = 18
		if state["left"]:
			state["left"] -= 1
			return re.sub(r"csel (x[4-9]), \1, (x2[5-9]|x30), hs$", r"csel \1, \1, \1, hs", s)
		return s
	return f

# Deliberately broken listings or generator constants: (name, kernels, line transform or constants). Each must fail.
CMP = re.compile(r"sbcs xzr, (x2[5-9]|x30), (x[4-9]|x1[0-5])$")
MUTANTS = [
	("min_top_limb_only", ("lazy", "root"), lambda s: CMP.sub(r"subs xzr, \1, \2", s)),
	("min_z_flag", ("lazy", "root"), lambda s: re.sub(r"(csel x(?:[4-9]|1[0-5]), x(?:2[5-9]|30), x(?:[4-9]|1[0-5])), lo$",
													 r"\1, ls", s)),
	("min_dropped_csel", ("lazy", "root"), "drop:csel x4, x25, x4, lo"),
	("b_select_flipped", ("lazy", "root"), lambda s: re.sub(r"(csel (x1[0-5]), x(?:2[5-9]|30), \2), hs$", r"\1, lo", s)),
	("presum_unreduced", ("root",), lambda s: re.sub(r"csel (x[4-9]), \1, (x2[5-9]|x30), hs$", r"csel \1, \2, \2, hs", s)),
	("u0_negation_unreduced", ("root",), "u0"),
	("operand_carry_dropped", ("lazy",), lambda s: re.sub(r"adcs (x2[01]), \1, (x2[89])$", r"add \1, \1, \2", s)),
	("d_not_doubled", ("lazy",), "drop:adds x4, x4, x4"),
	("ea_read_as_ed", ("lazy",), {"EA_OFF": L.ED_OFF}),
	("operands_regrouped", ("lazy", "root"), {"OPND": (("U", 0), ("x", 2), ("x", 1), ("x", 3), ("U", 1), ("x", 5))}),
	("node_offset_p", (), {"EOFF": P}),
]

def mutate(text, how, kern):
	if isinstance(how, dict):
		return getattr(L, kern + "_kernel")()[0].text()
	lines = text.split("\n")
	if how == "u0":
		how = u0_unreduced()
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
	text = {name: open(os.path.join(HERE, f"{name}.s")).read() for name in ("lazy", "root", "node")}
	progs = {name: emu.compile_prog(emu.parse(t)) for name, t in text.items()}
	e = (P ** 12 - 1) // C.Rr
	gens = [C.f12pow(C.rand12(rng), e) for _ in range(3)]
	fails, counts, kinds = check_all(progs, pool, n_cases, rng, gens)
	res = {"cases": kinds, "executed_counts": {k: sorted(v) for k, v in counts.items()}, "fails": fails}
	if sys.argv[3:] == ["mutants"]:
		caught = {}
		for name, kerns, how in MUTANTS:
			keep = {k: getattr(L, k) for k in how} if isinstance(how, dict) else {}
			try:
				for k, v in keep.items():
					setattr(L, k, how[k])
				bad = dict(progs)
				for k in kerns:
					t = mutate(text[k], how, k)
					assert t != text[k], name
					bad[k] = emu.compile_prog(emu.parse(t))
				f, _, _ = check_all(bad, pool, 4, random.Random(seed + 1), gens, stop_at_first=True)
			finally:
				for k, v in keep.items():
					setattr(L, k, v)
			caught[name] = sum(f.values()) > 0
		res["mutants_caught"] = caught
		assert all(caught.values()), caught
	print(json.dumps(res))
	assert all(len(v) == 1 for v in counts.values()) and not any(fails.values())
	print("PASS")

if __name__ == "__main__":
	main()
