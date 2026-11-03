"""evidence16: the experiment of proof.md (six parts, see proof.md Claim); stdout is one JSON object, progress on stderr."""
import base64, json, lzma, math, os, sys, time
from decimal import Decimal, localcontext
from fractions import Fraction as F
from random import Random
_HERE = os.path.dirname(os.path.abspath(__file__))
_DIR = os.path.join(_HERE, "files") if os.path.exists(os.path.join(_HERE, "files", "model16.py")) else _HERE
sys.path.insert(0, _DIR)
import common12 as C
import model16 as M
import owner10 as O
O.install(M)
import descent16 as DS
import store16
from threshold_model import Retention
import ideal_canonical as IC
import lattice_count_bounds as LCB
import embedding_bounds as EB
import filtsim
P = 4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787
U = C.U
H = C.H
LN2 = math.log(2)
CLAIM = 120.750399548
# ---- Murphy alpha over the published prime range (see proof.md section 2) ----
import retain16
ALPHA_PMAX = 2000
ALPHA_F, ALPHA_G = 4.868096, 0.843643   # alpha_contrib summed over p <= 2000 (42e80f80 --alpha; primitive-correct p=2 double-root term)
ALPHA_G_PMAX = 50000
ALPHA_G_50K = 0.731544   # alpha_contrib summed over p < 50000 (g side)
ALPHA_F_50K = 5.085036   # alpha_contrib summed over p < 50000 (f side; 42e80f80's --alpha50k routine, 280 s); the charged 5.15 is above it
def alpha_contrib(p):
	"""Expected p-adic valuation of f and g over the prime ideals above p, from
	the Kummer-Dedekind structure: alpha = sum_p log2(p)(1/(p-1) - sum_P d1*v_P),
	conditional on base-primitivity (counted relations are base-primitive).
	Simple roots use q/(q^2-1) = (q/(q+1))*(1/(q-1)): b is a unit with
	probability q/(q+1) given primitivity (DESCENT.md), else the homogeneous
	value is a unit; given b a unit the Hensel lift gives 1/(q-1).
	Repeated roots at p > 2 use 1/(q+1) (conservative: below the two certain
	tails 2/(q+1), so alpha is overestimated).
	At p = 2, g(x) = x^2 + U*y has double roots (mu = 2) above q = 4 (d1 = 2)
	and q = 16 (d1 = 4); v2(U) = 16 and y is a unit mod any prime above 2
	(H(0) = -1). Given b a unit (prob q/(q+1)), v = v_P(a^2 + U*y*b^2) with
	v_P(U) = 16: for j <= 16, v >= j iff a = 0 mod P^ceil(j/2), exact tails
	q/(q+1)*q^-ceil(j/2); j > 16 bounded by the multiplicity bound
	2*q/(q+1)*q^-ceil(j/2) (retain16.local_pdf tail). Total per prime at most
	2q/(q^2-1)*(1+q^-8). This is the primitive-correct valuation (the v6
	2/(q-1) omitted the q/(q+1) factor)."""
	vf = vg = 0.0
	for d1, e, mf, mg in retain16.structure(p):
		q = float(p ** d1)
		for side, mults in enumerate((mf, mg)):
			v = 0.0
			for mu in mults:
				if mu == 1: v += q / (q * q - 1.0)
				elif p == 2 and side == 1 and mu == 2: v += (2.0 * q / (q * q - 1.0)) * (1.0 + q ** -8)
				else: v += 1.0 / (q + 1.0)
			if side == 0: vf += d1 * v
			else: vg += d1 * v
	t = math.log(p) / math.log(2.0)
	return t * (1.0 / (p - 1.0) - vf), t * (1.0 / (p - 1.0) - vg)
DECLARED = {'ledger': ('<=', 120.750399548), 'relcol': ('<=', 120.348119967), 'la': ('<=', 118.440377117), 'log2_N': ('<=', 59.4521),
	'log2_rels': ('>=', 63.8902), 'rels_2se': ('>=', 63.68), 'margin': ('>=', 0.3014), 'work_se': ('<=', 0.01961),
	'two_log': ('<=', 115.3704), 'markov': ('<=', 0.577703), 'log2_EN': ('<=', 18.12), 'over': ('<=', -7.38),
	'rep': ('<=', -14.30), 'mrep': ('<=', -8.58), 'cheb': ('<=', -30), 'boot_share': ('<=', -26.07), 'boot_phi': ('<=', 0.001763641),
	'success': ('>=', 0.410103), 'gen': ('<=', 107.644448746)}
PARTS = ("core", "ledger", "descent", "lattices", "retention", "losses")
RELS_EST, YSE = 63.8903, 0.1014   # sub-region retained-relation estimate (log2) and its pool standard error ([4]); the descent input is RELS_EST - 2 YSE floored to 63.68
_REPORT = []
T0 = time.time()
def load_stats():
	"""stats16.b85: LZMA-compressed base85 JSON of the measured statistics."""
	import lzma, base64
	return json.loads(lzma.decompress(base64.b85decode(open(os.path.join(_DIR, "stats16.b85")).read())))
def say(*args):
	line = " ".join(str(a) for a in args)
	_REPORT.append(line)
	print(line, file=sys.stderr, flush=True)
def f6_mulp(A, B):
	return [c % P for c in C.mul([c % P for c in A], [c % P for c in B])]
def f6_pow(A, n):
	r, base = [1], [c % P for c in A]
	while n:
		if n & 1:
			r = f6_mulp(r, base)
		base = f6_mulp(base, base)
		n >>= 1
	return (r + [0] * 6)[:6]
def check_resultant():
	F0 = [0] * 6
	for d, c in C.TERMS:
		term = (c + [0] * 6)[:6]
		for _ in range(d // 2):
			term = C.mul(term, [0, -U])
		F0 = [x + y for x, y in zip(F0, term)]
	r = C.norm(F0)
	ok = r != 0 and r % P == 0
	say(f"[0] Res_y(h, f(x^2 -> uy)): {r.bit_length()} bits, 0 mod p {'ok' if ok else 'FAIL'}")
	return ok
def check_tower():
	yp = f6_pow([0, 1], P)
	cur, ok_irr = [0, 1, 0, 0, 0, 0], True
	for k in range(1, 7):
		acc = [0] * 6
		for coef in reversed(cur):
			acc = f6_mulp(acc, yp)
			acc[0] = (acc[0] + coef) % P
		cur = acc
		is_y = cur == [0, 1, 0, 0, 0, 0]
		if (k in (1, 2, 3) and is_y) or (k == 6 and not is_y):
			ok_irr = False
	uy = f6_mulp([(-U) % P], [0, 1])
	chi = f6_pow(uy, (P ** 6 - 1) // 2)
	ok_sq = chi == [P - 1, 0, 0, 0, 0, 0]
	say(f"[0b] h irreducible mod p {'ok' if ok_irr else 'FAIL'}; (u y)^((p^6-1)/2) = -1 {'ok' if ok_sq else 'FAIL'}")
	cert = IC.field_certificate(P)
	ok_c = cert['squarefree'] and cert['irreducible_mod_p'] and cert['maximal_power_basis'] and cert['discriminant'] == 58673
	say(f"[0c] Disc(h) = {cert['discriminant']} ({cert['discriminant_factors']}), O_h = Z[y] {'ok' if ok_c else 'FAIL'}")
	return ok_irr and ok_sq and ok_c
def check_rho_and_semismooth():
	ref = {2: 0.30685281944, 3: 0.048608388, 4: 0.0049109256, 6: 1.9649696e-5, 8: 3.2320693e-8, 10: 2.7701718e-11}
	err = max(abs(DS.rho_series(u) - v) / v for u, v in ref.items())
	ok_br = all(DS.rho_lo(u) <= DS.rho_series(u) <= DS.rho_hi(u) + 1e-300 for u in (1.5, 2.3, 7.7, 15.2, 33.3, 41.0))
	rows, h = DS.window_rows(80.0, 69.0, 48, 48)
	ok_w = abs(sum(rows[(1, 0)]) - math.log(80 / 69)) < 1e-12
	say(f"[1] rho: series vs reference {err:.1e} {'ok' if err < 1e-6 else 'FAIL'}; brackets {'ok' if ok_br else 'FAIL'}; weight sum = ln(n/base) "
		f"{'ok' if ok_w else 'FAIL'}")
	return err < 1e-6 and ok_br and ok_w
def uniform_half_ellipsoid(q, inverse_q, radius2, rng):
	cap = []
	for j in range(6):
		c = F(radius2) * inverse_q[j][j]
		cap.append(math.isqrt(c.numerator // c.denominator))
	while True:
		a = [rng.randrange(-c, c + 1) for c in cap[:-1]] + [rng.randrange(1, cap[-1] + 1)]
		if C.quadratic(q, a) <= radius2:
			return a, C.quadratic(q, a)
def regenerate(cfg, seed, index):
	qa, da, ia = C.integer_metric(cfg['matrix_a']); qb, db, ib = C.integer_metric(cfg['matrix_b'])
	A2, theta = cfg['production_A_squared'], F(cfg['production_theta'])
	rng = Random(f"v11-production-ideal-key:{seed}:{index}")
	a, aq = uniform_half_ellipsoid(qa, ia, int(A2 * theta * da), rng)
	b, bq = uniform_half_ellipsoid(qb, ib, (A2 - F(aq, da)) * db, rng)
	return a, b, aq, bq
def coprime_pair(a, b):
	da, db = C.norm(a), C.norm(b)
	d = math.gcd(da, db)
	if d == 1:
		return True
	rows = [C.mul(v, [0] * i + [1]) for v in (a, b) for i in range(6)]
	key = IC.row_hnf(rows, d)
	return all(key[i][i] == 1 for i in range(6))
def verify_pool(cfg, pool, samples, disks):
	qa, da, ia = C.integer_metric(cfg['matrix_a']); qb, db, ib = C.integer_metric(cfg['matrix_b'])
	A2 = cfg['production_A_squared']; theta_reg = F(cfg['production_theta_region']); n_out = 0
	index = {(r['seed'], r['index']): r for r in pool}
	n_can = n_prim = 0; bad = []
	t0 = time.time()
	idx = sorted(set([r['index'] for r in pool if r['seed'] == 33][:samples] + list(range(samples))))
	for i in idx:
		a, b, aq, bq = regenerate(cfg, 33, i)
		can, work, witness = IC.canonical_a(a, cfg['matrix_a'], 100000)
		row = index.get((33, i))
		if can and row is None and F(aq, da) / A2 > theta_reg:
			n_out += 1; continue   # canonical but outside the collected sub-region: not stored (pool12 header 'subregion')
		if can != (row is not None):
			bad.append(('canonical flag', i, can)); continue
		if not can:
			continue
		n_can += 1
		qa_frac = F(aq, da) / A2
		if round(qa_frac * 65536) != round(row['qa'] * 65536):
			bad.append(('Q_a', i)); continue
		nf, ng = C.norms(a, b)
		if abs(math.log2(nf) - row['nf']) > 1 / 16 or abs(math.log2(ng) - row['ng']) > 1 / 16:
			bad.append(('norms', i)); continue
		if coprime_pair(a, b) != row['primitive']:
			bad.append(('primitive', i)); continue
		n_prim += row['primitive']
		radius = (A2 - F(aq, da)) * db
		t = sum(qb[0][j] * b[j] for j in range(1, 6))
		c = sum(qb[j][k] * b[j] * b[k] for j in range(1, 6) for k in range(1, 6)) - radius
		iv = C.integer_interval(qb[0][0], t, c)
		blo = (b[0] // 64) * 64; lo, hi = max(iv.start, blo), min(iv.stop - 1, blo + 63)
		if 16 * max(1, (iv.stop - iv.start) // 16) != row['rlen'] or hi - lo + 1 != row['blen']:
			bad.append(('row/block length', i)); continue
		est = EB.NormBounds(a, disks)
		leaves = est.adaptive_row(a, b[1:], lo, hi, singleton_policy='reject', width=128)
		at = [l for l in leaves if l[0] <= b[0] <= l[1]][0]
		rej = isinstance(at[2], dict)
		if est.counts['block_evaluations'] != row['evals'] or len(leaves) != row['leaves'] or rej != row['rejected']:
			bad.append(('nbound statistics', i)); continue
		if not rej:
			(fl, fu), (gl, gu) = at[2]
			if not (fl <= 16 * math.log2(nf) <= fu and gl <= 16 * math.log2(ng) <= gu):
				bad.append(('bound does not contain the exact norm', i)); continue
	say(f"[2] {len(idx)} seed-33 indices regenerated ({time.time()-t0:.1f}s): {n_can} canonical in the sub-region ({n_out} canonical outside it, not stored), {n_prim} primitive; every stored "
		f"field {'ok' if not bad else 'FAIL ' + str(bad[:3])}")
	return not bad, n_can
def poly_roots_count(poly, p):
	poly = [c % p for c in poly]
	while poly and poly[-1] == 0:
		poly.pop()
	if len(poly) <= 1:
		return 0
	def pmod(a, m):
		a = a[:]
		inv = pow(m[-1], -1, p)
		while len(a) >= len(m):
			if a[-1]:
				z = a[-1] * inv % p; s = len(a) - len(m)
				for j, v in enumerate(m):
					a[s + j] = (a[s + j] - z * v) % p
			a.pop()
			while a and a[-1] == 0:
				a.pop()
		return a
	def pmulmod(a, b, m):
		r = [0] * (len(a) + len(b) - 1)
		for i, x in enumerate(a):
			if x:
				for j, y in enumerate(b):
					r[i + j] = (r[i + j] + x * y) % p
		return pmod(r, m)
	r, base, e = [1], [0, 1], p
	while e:
		if e & 1:
			r = pmulmod(r, base, poly)
		base = pmulmod(base, base, poly); e >>= 1
	r = (r + [0, 0])[:max(len(r), 2)]; r[1] = (r[1] - 1) % p
	a, b = poly, pmod(r, poly) if len(r) >= len(poly) else [c % p for c in r]
	while b and any(b):
		a, b = b, pmod(a, b)
	return len(a) - 1
def nbar(lo_bits, hi_bits):
	from threshold_model import primes_upto
	ps = [q for q in primes_upto(2 ** hi_bits) if q >= 2 ** lo_bits]
	tot = 0
	rng = Random(1)
	for p in ps:
		for s in sorted(set(roots_modp(H, p, rng))):   # the roots of h mod p
			tot += poly_roots_count([1 + s - s * s - s ** 3, 0, s ** 5, 0, 0, 0, 2 * s ** 3, 0, 0, 0, -2 * s, 0, 1], p) + \
				poly_roots_count([U * s, 0, 1], p)
	return tot / len(ps), len(ps)
def _asm(name):
	"""compiled kernel listings: kernels_asm.b85 = base85(lzma(text)), sections '@@ la_mac.s' and '@@ bw_dot.s'"""
	t = lzma.decompress(base64.b85decode("".join(open(os.path.join(_DIR, "kernels_asm.b85")).read().split()))).decode()
	return t.split("@@ " + name + "\n")[1].split("\n@@ ")[0].splitlines()
def check_la_kernel(name="la_mac.s", start="LBB9_3:"):
	on, body = False, []
	for ln in _asm(name):
		s = ln.strip()
		if s.startswith(start):
			on = True; continue
		if on:
			if s.startswith("LBB") or (s.startswith(";") and "bb." in s):
				break
			if s and not s.startswith(";") and not s.startswith("."):
				body.append(s.split()[0])
	return len(body), [m for m in body if m.startswith("b.") or m == "b"]
def check_la_func(fn):
	"""static count of a loop-free function (all branches forward); calls are listed"""
	L = [x.strip() for x in _asm("la_mac.s")]
	i = L.index(fn + ":") + 1; n, seen, fwd, calls = 0, set(), True, []
	while i < len(L) and not L[i].startswith(("_", ";")):
		x = L[i]; i += 1
		if x.startswith("LBB"): seen.add(x[:-1])
		elif x:
			n += 1; fwd = fwd and x.split()[-1] not in seen
			if x.split()[0] == "bl": calls.append(x.split()[-1])
	return n, fwd, calls
def la_row_tail():
	"""row tail per row and product and one 9-limb reduction, as itemised in STORAGE.md"""
	f = {k: check_la_func(k) for k in ("_send_pair", "_combine", "_finish", "_reduce9", "_fold9", "_reduce6")}
	ok = all(v[1] for v in f.values()) and f["_finish"][2] == ["_fold9", "_reduce6"] and f["_reduce9"][2] == ["_fold9", "_reduce6"] \
		and not any(f[k][2] for k in ("_send_pair", "_combine", "_fold9", "_reduce6"))
	sub = f["_fold9"][0] + f["_reduce6"][0]
	k = O.ACTIVE_OWNERS
	row = (k - 1) * (f["_send_pair"][0] + 18) + (k - 1) * (f["_combine"][0] + 2) + f["_finish"][0] + sub + 2 + 10 * k + 40 + 168   # k owners (owner10.py, STORAGE.md)
	red9 = f["_reduce9"][0] + sub + 2
	return row, red9, {k: v[0] for k, v in f.items()}, ok
def ball_counts_check():
	"""exact N_12(R^2) by the integer-point recurrence for R <= 8 against CELLS;"""
	import functools
	@functools.lru_cache(None)
	def Nd(d, s):
		if s < 0:
			return 0
		if d == 0:
			return 1
		r = math.isqrt(s)
		return sum(Nd(d - 1, s - x * x) for x in range(-r, r + 1))
	ok = True
	for R in (4, 6, 8):
		n12 = Nd(12, R * R)
		ok &= (n12 - 1) // 2 == DS.CELLS[R]
	V12 = math.pi ** 6 / 720
	for R, half in DS.CELLS.items():
		n = 2 * half + 1
		ok &= V12 * (R - math.sqrt(3)) ** 12 <= n <= V12 * (R + math.sqrt(3)) ** 12
	return ok
def roots_modp(poly, p, rng):
	poly = [c % p for c in poly]
	while poly and poly[-1] == 0:
		poly.pop()
	if len(poly) <= 1:
		return []
	def pmod(a, m):
		a = a[:]
		inv = pow(m[-1], -1, p)
		while len(a) >= len(m):
		  if a[-1]:
		    z = a[-1] * inv % p; s = len(a) - len(m)
		    for j, v in enumerate(m):
		      a[s + j] = (a[s + j] - z * v) % p
		  a.pop()
		  while a and a[-1] == 0:
		    a.pop()
		return a
	def pmul(a, b, m):
		r = [0] * (len(a) + len(b) - 1)
		for i, x in enumerate(a):
		  if x:
		    for j, y in enumerate(b):
		      r[i + j] = (r[i + j] + x * y) % p
		return pmod(r, m)
	def powmod(a, e, m):
		r = [1]
		while e:
		  if e & 1:
		    r = pmul(r, a, m)
		  a = pmul(a, a, m); e >>= 1
		return r
	xp = powmod([0, 1], p, poly)
	xp = xp + [0] * max(0, 2 - len(xp)); xp[1] = (xp[1] - 1) % p
	a, b = poly, pmod(xp, poly)
	while b and any(b):
		a, b = b, pmod(a, b)
	if not a or not any(a):
		return []
	inv = pow(a[-1], -1, p); g = [x * inv % p for x in a]
	def split(g):
		if len(g) <= 1:
		  return []
		if len(g) == 2:
		  return [(-g[0] * pow(g[1], -1, p)) % p]
		while True:
		  c = [rng.randrange(p) for _ in range(len(g) - 2)] + [rng.randrange(1, p), 1]
		  t = powmod(c, (p - 1) // 2, g)
		  t = t + [0] * max(0, 1 - len(t)); t = [(t[0] - 1) % p] + t[1:]
		  h = pmod(t, g) if len(t) >= len(g) else t
		  aa, bb = g, h
		  while bb and any(bb):
		    aa, bb = bb, pmod(aa, bb)
		  if aa and 1 < len(aa) < len(g):
		    iv = pow(aa[-1], -1, p); hh = [x * iv % p for x in aa]
		    qq = g[:]; v = [0] * (len(g) - len(hh) + 1)
		    for k in range(len(v) - 1, -1, -1):
		      v[k] = qq[k + len(hh) - 1]
		      for j, x in enumerate(hh):
		        qq[k + j] = (qq[k + j] - v[k] * x) % p
		    assert not any(qq)
		    return split(hh) + split(v)
	return sorted(split(g))
def repeat_bound(dres):
	"""[7] node counter N_MAX = 2^NODE_CAP_LOG2 over both logs: (a) Markov P(abort) <= E[N]/N_MAX; (b) special-q repeats <= C(N_MAX, 2) 2^6/|smallest bin|;
	(c) missing-ideal repeats under the shared missing set <= C(N_MAX, 2) sum_I m_I wh_I^2; (d) Chebyshev exponent of the zero-occurrence mass (descent16.miss_shared). Rounded UP to 0.01 bit."""
	import counts16
	EN, _ = counts16.expected_nodes(DS.P(), dres)
	pmin = min(math.log2(2 * (2 ** DS.log2li(b[1]) - 2 ** DS.log2li(b[0]))) for b in dres['bins'])
	over = math.ceil((math.log2(EN) - M.NODE_CAP_LOG2) * 100) / 100
	coll = math.ceil((2 * M.NODE_CAP_LOG2 - 1 + 6 - pmin) * 100) / 100
	mw, cheb = DS.miss_shared(DS.P())
	mrep = math.ceil((2 * M.NODE_CAP_LOG2 - 1 + mw) * 100) / 100
	return math.log2(EN), over, coll, pmin, mw, mrep, math.ceil(cheb * 100) / 100
def audit(rng_seed=31001):
	import descent_witness as W
	assert (W.P ** 6 - 1) % W.R != 0 and ((W.P ** 12 - 1) // W.R) % W.R != 0
	rng = Random(rng_seed)
	while True:
		x = ([rng.randrange(W.P) for _ in range(6)], [rng.randrange(W.P) for _ in range(6)])
		if any(x[0] + x[1]):
			break
	w0 = W.pow12(x, (W.P ** 12 - 1) // W.R)
	tt = rng.randrange(1, W.R)
	z = W.pow12(w0, tt)
	B, swaps = W.boot_lattice(z)
	bst = W.boot_ball_stats(z, B, rng_seed + 1, points=6, radii=(96,))
	boot_bins_max = max(x[1] for x in DS.P.BOOT_BINS)
	sq = []
	Pm = DS.P()
	cfmax = max(b[1] for b in Pm.QBINS); cgmax = max(b[2] for b in Pm.QBINS)
	env_ok = True
	for bits, side, seed in ((260, 'f', 13300),):
		r = W.specialq_ball_stats(bits, side, seed, points=6, radii=(96,))
		for st_ in r['stats']:
			env_ok &= st_['cf_max'] <= cfmax and st_['cg_max'] <= cgmax
		sq.append(dict(bits=bits, side=side, seed=seed, swaps=r['swaps'], stats=[{k: (round(v, 2) if isinstance(v, float) else v) for k, v in s.items()} for s in r['stats']]))
		if side == 'f':
			units = W.small_units(2)
			um = W.unit_mates(r['basis'], seed + 7, 96, 8, units)
	return dict(boot=dict(swaps=swaps, stats=[{k: (round(v, 2) if isinstance(v, float) else v) for k, v in s.items()} for s in bst],
						 within_model_bins=all(s['cg_max'] <= boot_bins_max for s in bst)),
				specialq=sq, specialq_within_bins=env_ok,
				unit_mates=dict(points=um['points'], units=um['units'], counts=um['counts'], distinct_over_raw=round(um['distinct_over_raw'], 3)),
				storage=store16.check())
def check_fields():
	"""[13] r is unramified in K_f and K_g (premise of H3, proof.md section 5)."""
	import retain16 as RT, descent_witness as W
	R = W.R; Hr = [x % R for x in H]
	def res(a, b):
		"""Res(a, b) over F_r by the Euclidean recurrence."""
		a = RT._trim([x % R for x in a]); b = RT._trim([x % R for x in b]); s = 1
		while b and len(b) > 1:
			da, db = len(a) - 1, len(b) - 1; rem = RT.pmod(a, b, R)
			if not rem:
				return 0
			s = s * (-1) ** (da * db) * pow(b[-1], da - (len(rem) - 1), R) % R; a, b = b, rem
		return s * pow(b[0], len(a) - 1, R) % R if b else 0
	def norm_poly(terms, deg):
		"""N_{K_h/Q}(s(X)) mod r for s = sum_d c_d(y) X^d monic of degree deg: values at 6 deg + 1 points, Newton interpolation."""
		n = 6 * deg + 1; xs = list(range(n))
		ys = [res(Hr, [sum(c[i] * pow(k, d, R) for d, c in terms if i < len(c)) % R for i in range(6)]) for k in xs]
		for j in range(1, n):
			for i in range(n - 1, j - 1, -1):
				ys[i] = (ys[i] - ys[i - 1]) * pow(xs[i] - xs[i - j], -1, R) % R
		poly = [ys[-1]]
		for i in range(n - 2, -1, -1):
			poly = [(-xs[i] * poly[0] + ys[i]) % R] + [(poly[k - 1] - xs[i] * poly[k]) % R for k in range(1, len(poly))] + [poly[-1]]
		return poly
	out = {}
	for name, terms, deg in (('f', C.TERMS, 12), ('g', RT.G_TERMS, 2)):
		poly = norm_poly(terms, deg)
		out['deg_' + name] = len(poly) - 1
		out['sqf_' + name] = poly[-1] == 1 and len(poly) == 6 * deg + 1 and len(RT.pgcd(poly, RT.pderiv(poly, R), R)) == 1
	out['unramified'] = out['sqf_f'] and out['sqf_g']
	out['disc_h'] = C.norm([-1, 2, 3, 0, 0, 6])                     # |disc h| = |N(h'(y))|
	out['sig'] = signatures()
	return out
def signatures():
	"""[13] Sturm counts of real embeddings of K_f, K_f^+ = K_h(x^2), K_g; sigma-odd unit rank."""
	def tr(p):
		while p and p[-1] == 0: p = p[:-1]
		return p
	def sturm(p):
		s = [tr([F(c) for c in p])]; s.append([i * c for i, c in enumerate(s[0])][1:])
		while len(s[-1]) > 1:
			a, b = s[-2][:], s[-1]
			while len(a) >= len(b):
				q = a[-1] / b[-1]; d = len(a) - len(b); a = tr([x - q * b[i - d] if i >= d else x for i, x in enumerate(a)][:-1])
			if not a: break
			s.append([-x for x in a])
		return s
	def ch(s, x):
		v = [t for t in (sum(c * x ** i for i, c in enumerate(q)) for q in s) if t]
		return sum((a > 0) != (b > 0) for a, b in zip(v, v[1:]))
	nre = lambda p, B=F(10) ** 60: (lambda s: ch(s, -B) - ch(s, B))(sturm(p))
	sh = sturm(H); ys = []
	def iso(a, b):
		n = ch(sh, a) - ch(sh, b)
		if n == 1 and b - a < F(1, 10 ** 30): ys.append((a + b) / 2)
		elif n: iso(a, (a + b) / 2); iso((a + b) / 2, b)
	iso(F(-10), F(10))
	rf = rfp = rg = 0
	for y in ys:
		fx = [0] * 13
		for d, c in C.TERMS: fx[d] = sum(F(ci) * y ** i for i, ci in enumerate(c))
		rf += nre(fx); rfp += nre(fx[0::2]); rg += nre([U * y, 0, 1])
	rk = lambda n, r1: r1 + (n - r1) // 2 - 1
	return dict(f=(rf, 36 - rf // 2), fp=(rfp, 18 - rfp // 2), g=(rg, 6 - rg // 2), odd=rk(72, rf) - rk(36, rfp) + rk(12, rg) - rk(6, len(ys)))
def success_of(two_log, cap, bt, over, coll, mrep, share, log2_N):
	"""success lower bound (proof.md section 7): one minus the descent terms, the Block-Wiedemann/main-phase 2^-191 and 2/r."""
	markov = 2 ** (two_log - cap)
	p_desc = 1 - 2 * bt['phi'] - markov - 2.0 ** over - 2.0 ** coll - 2.0 ** mrep - 2.0 ** share - 2.0 ** -30 - 2.0 ** -48 - 2.0 ** -191 - 2.0 / R_ORDER   # 2^-30: the zero-occurrence concentration bound (DESCENT.md); 2^-191: Wiedemann <= 2^-192 plus main phase >= 1 - 2^-192
	return p_desc, markov
R_ORDER = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
def fields_line():
	t0 = time.time(); fl = check_fields()
	say(f"[13] number fields ({time.time()-t0:.1f}s): N(f(X)) (degree {fl['deg_f']}) and N(g(X)) (degree {fl['deg_g']}) squarefree mod r {'ok' if fl['unramified'] else 'FAIL'}; "
		f"disc h = {fl['disc_h']}; signatures K_f {fl['sig']['f']}, K_f+ {fl['sig']['fp']}, K_g {fl['sig']['g']}: odd unit rank {fl['sig']['odd']} = SM columns {M.N_SM} {'ok' if fl['sig']['odd'] == M.N_SM else 'FAIL'}")
	return fl['unramified'] and fl['sig']['odd'] == M.N_SM
def main():
	cfg_in = json.load(sys.stdin) if not sys.stdin.isatty() else {}
	samples = int(cfg_in.get("samples", 8)); seed = int(cfg_in.get("seed", 1)); part = str(cfg_in.get("part", "core"))
	if part not in PARTS:
		raise SystemExit(f"part must be one of {PARTS}")
	say(f"evidence: stnfs (bls12-381-fixed-g2-invert-g1); part={part} samples={samples} seed={seed}; unit {M.UNIT:.0f} instructions")
	ok0 = check_resultant(); ok0b = check_tower(); ok1 = check_rho_and_semismooth()
	cfg = json.load(open(os.path.join(_DIR, "config.json")))
	hdr, pool = M.load_pool(open(os.path.join(_DIR, "pool12-1.b85")).read() + open(os.path.join(_DIR, "pool12-2.b85")).read())
	n_total = sum(hdr['samples'].values())
	theta_pool = F(cfg['production_theta']); theta_reg = F(cfg['production_theta_region'])
	assert theta_reg <= theta_pool == M.OP['theta']
	pool_reg = [r for r in pool if r['qa'] <= float(theta_reg)]
	n_eff = n_total * float((theta_reg / theta_pool) ** 3)
	rmin = math.isqrt(int(theta_reg * cfg['production_A_squared']))
	corr = LCB.correction(cfg, rmin)
	geo = M.Geometry(cfg, theta_reg, F(corr['safe_yield_factor']), F(corr['safe_work_population_factor']))
	say(f"[3] collected region theta={theta_reg} (pool at theta={theta_pool}; the sub-region's rows are a uniform sample of it), a-radius {rmin}; "
		f"eligible a 2^{math.log2(geo.va):.3f}, largest half-slice 2^{math.log2(geo.vb):.3f} cells; certified count products "
		f"{float(corr['safe_yield_factor']):.9f} (yield), {float(corr['safe_work_population_factor']):.9f} (work) at every radius >= {rmin}; "
		f"{len(pool_reg)} canonical rows of {len(pool)} in the region ({sum(r['primitive'] for r in pool_reg)} primitive), "
		f"effective sample {n_eff:.3f} of {n_total} draws")
	op = M.OP
	log2_U = math.log2((M.NBAR / 2) * 2 ** M.log2_li(op['lb']) + 2 ** 31)
	la = M.linear_algebra(log2_U, op['excess'])
	cap = M.DESCENT_CAP
	checks = dict(resultant=ok0, tower=ok0b, rho_semismooth=ok1)
	vals = {}
	extra = {}
	full = part in ("ledger", "descent")
	if part in ("core", "ledger"):
		if part == "ledger":
			disks = EB.root_disks()
			ok2, n_can = verify_pool(cfg, pool, samples, disks)
			checks['pool'] = ok2
		retain = Retention(int(op['lb']), int(op['lx2']), step=.25, max_bits=850, tilt=.05)
		cost = Retention(850, int(op['lx2']), step=.25, max_bits=851, tilt=.005)
		t0 = time.time()
		e = M.evaluate(pool_reg, n_eff, geo, op, retain, cost, return_work=True)
		wc = e['work_contribs']; ncan = len(wc); mw = sum(wc) / ncan
		sw = math.sqrt(sum((x - mw) ** 2 for x in wc) / (ncan - 1))
		cond_se = sw / math.sqrt(ncan) / mw
		pcan = ncan / n_total
		work_se = math.sqrt(cond_se ** 2 + (1 - pcan) / ncan)
		allow_need = 1 + M.WORK_SE_K * work_se
		ok_allow = M.CANON_WORK_ALLOWANCE >= allow_need
		say(f"[4a] work estimator: {ncan} contributing draws of {n_total} ({n_eff:.1f} effective); row SE {100*cond_se:.3f}%, "
			f"incidence {100*math.sqrt((1-pcan)/ncan):.3f}%, total {100*work_se:.3f}%; allowance "
			f"{M.CANON_WORK_ALLOWANCE:.4f} >= 1 + {M.WORK_SE_K:.2f} SE = {allow_need:.5f} {'ok' if ok_allow else 'FAIL'}")
		ok_relmeans = e['nf_rel'] >= 69 * DS.P.MAIN_UF and e['ng_rel'] >= 69 * DS.P.MAIN_UG
		rels_2se = e['log2_rels'] - 2 * e['yield_se_log2']
		ok_rels_input = rels_2se >= DS.P.LOG2_RELS
		rf = retain.at(e['nf_mean'], op['slack_f']); rg = retain.at(e['ng_mean'], op['slack_g'])
		say(f"[4] B=2^{op['lb']:.0f}, X2=2^{op['lx2']:.0f}, slacks ({op['slack_f']}, {op['slack_g']}), excess {op['excess']}, "
			f"width {op['width_bits']}, blocks of {op['L0']}, alpha ({op['alpha_f']}, {op['alpha_g']}); mean norms "
			f"{e['nf_mean']:.2f}/{e['ng_mean']:.2f} (relation-weighted {e['nf_rel']:.2f}/{e['ng_rel']:.2f} >= descent inputs {69*DS.P.MAIN_UF:.1f}/{69*DS.P.MAIN_UG:.1f} {'ok' if ok_relmeans else 'FAIL'}); retention at the means {rf['retention']:.4f}/{rg['retention']:.4f}; "
			f"rels 2^{e['log2_rels']:.4f} vs need "
			f"2^{e['log2_need']:.4f}, margin {math.floor(e['margin']*1000)/1000:.3f} >= {M.MARGIN} ({'feasible' if e['feasible'] else 'INFEASIBLE'}); pool SE {e['yield_se_log2']:.4f} bits; "
			f"estimate - 2 SE = {rels_2se:.4f} >= descent LOG2_RELS input {DS.P.LOG2_RELS:.3f} {'ok' if ok_rels_input else 'FAIL'} ({time.time()-t0:.1f}s)")
		say(f"[5] population 2^{e['log2_population']:.3f} cells at {e['per_cell']:.1f} instructions per cell (survivor cofactorisation 2^{e['c_surv_log2']:.4f}): relation collection 2^{e['log2_relcol']:.6f} units; "
			f"N/U {la['ratio']:.5f} x pad {M.PAD_N} x weight-200 ratio {M.W200_RATIO} x trap bloat (1 + {M.TRAP_BLOAT}) -> N = 2^{la['log2_N']:.4f}; linear algebra 2^{la['log2_la']:.6f} units ({M.C_ENTRY:.1f} instr/entry), generator 2^{la['log2_gen']:.6f}")
		if part == "ledger":
			t0 = time.time()
			Us, res = filtsim.run(op['lb'], e['nf_rel'] / op['lb'], e['ng_rel'] / op['lb'], op['excess'], seed, 12, kmax=64, targets=(100,), verbose=False)
			r100 = res[100]
			say(f"[6] filtering (scale 2^12, seed {seed}, {time.time()-t0:.1f}s): N/U = {1/r100[0]:.4f} at weight {r100[3]:.0f} (max |coefficient| {r100[2]}); "
				"2^19 basis at excess 1.1 (seed, U/N): " + ", ".join(f"({s_}, {un:.3f})" for (s_, un, mc, w_) in M.FILT_TABLE[1.1]) + f"; max coef {max(mc for (s_, un, mc, w_) in M.FILT_TABLE[1.1])}, weight <= {max(w_ for (s_, un, mc, w_) in M.FILT_TABLE[1.1])}")
			t0 = time.time()
			rec_ab = filtsim.missing_recon(seed, uf=e['nf_rel'] / op['lb'], ug=e['ng_rel'] / op['lb'], excess=op['excess'], addback=True)
			rec = rec_ab['bands']
			ok_rec = all(mf <= md for (lo, hi, n, mf, md, h3) in rec)
			say(f"[6b] missing logs after purge and back-substitution (2^12, seed {seed}, {time.time()-t0:.1f}s) per lambda band (occurring ideals, measured unknown fraction among them, uncalibrated P(Poisson <= 2), lost with >= 3): "
				+ "; ".join(f"[{lo:g},{hi:g}): {n}, {mf:.4f}, {md:.4f}, {h3}" for (lo, hi, n, mf, md, h3) in rec) + f" {'ok' if ok_rec else 'FAIL'}")
			ok_trap = rec_ab['lost3_after'] == 0 and rec_ab['added'] / rec_ab['U'] <= M.TRAP_SU_COVER
			say(f"[6t] trap add-back (same 2^12 run): {rec_ab['trapped']} trapped occ>=3 ideals, {rec_ab['added']} rows appended (s/U = {rec_ab['added']/rec_ab['U']:.5f} <= cover {M.TRAP_SU_COVER}), occ>=3 still unknown afterwards {rec_ab['lost3_after']} {'ok' if ok_trap else 'FAIL'}")
			cal = []
			for (lo, hi, n, mf, md, h3) in rec:
				c = next((c_ for lo_, hi_, c_ in DS.CALIB if lo_ <= lo < hi_), None)
				if c is None:
					cal.append((lo, hi, n, mf, None, None, True)); continue
				cal.append((lo, hi, n, mf, c, c, mf <= c + 1e-12))
			ok_tab = all(c >= 1.5 * DS.CALIB_MAX[(lo, hi)] for lo, hi, c in DS.CALIB if (lo, hi) in DS.CALIB_MAX)
			ok_cal = all(x[6] for x in cal) and ok_tab
			say(f"[7c] missing-log calibration rerun (2^12, seed {seed}): per band (lambda lo, measured fraction of OCCURRING ideals still unknown <= cover c_r, no allowance): "
				+ "; ".join((f"[{lo:g},{hi:g}) {mf:.4f} <= {c:.3f} {'ok' if ok_ else 'FAIL'}" if c is not None else f"[{lo:g},{hi:g}) {mf:.4f} (uncalibrated)") for lo, hi, n, mf, c, al, ok_ in cal)
				+ "; covers " + ", ".join(f"[{lo:g},{hi:g}) {c}" for lo, hi, c in DS.CALIB) + f" >= 1.5 x the 16-run maxima {DS.CALIB_MAX} {'ok' if ok_tab else 'FAIL'}; the never-occurring mass is charged in full (descent16.p_missing); overall {'ok' if ok_cal else 'FAIL'}")
			nb, np_ = nbar(8, 12)
			say(f"[9] degree-1 ideals per prime, both sides, {np_} primes in [2^8, 2^12): {nb:.3f} (model NBAR = {M.NBAR})")
			ok9 = nb <= M.NBAR * 1.01
			la_n, la_br = check_la_kernel()
			tail, red9, tail_n, tail_ok = la_row_tail()
			say(f"[10] row tail (la_mac.s): {tail_n} instructions (calls ok: {tail_ok}); {tail} per row <= C_ROW {M.C_ROW:.0f}; 9-limb reduction {red9} <= C_RED9 {M.C_RED9:.0f}")
			charged = M.ratio_N_over_U(op['excess']) * M.PAD_N
			worst14 = max(max(v) for k, v in M.SCALE_NU11.items() if k >= 14)
			w200max = max(max(v) for v in M.W200_11.values())
			nu_ok = charged >= worst14 and w200max <= M.W200_RATIO and all(1 / r_[1] <= charged / M.PAD_N * (1 + 1e-12) for r_ in M.FILT_TABLE[1.1])
			say(f"[10] weight-100 N/U at excess 1.1, max over seeds: " + ", ".join(f"2^{k} {max(v):.4f} (n={len(v)})" for k, v in sorted(M.SCALE_NU11.items()))
				+ f"; charged {charged:.5f} = {M.PAD_N} x 2^19 worst {charged / M.PAD_N:.5f} >= every value from 2^14 up ({worst14:.5f}); weight-200/100 ratio at excess 1.1 max {w200max:.4f} <= charged {M.W200_RATIO} (42e80f80's excess-1.2 data: max {max(max(v) for v in M.W200.values()):.4f}): {'ok' if nu_ok else 'FAIL'}")
			say(f"[10] la_mac.s group loop: {la_n} instructions {'ok' if la_n == 33 and la_br == ['b.ne'] else 'MISMATCH'}")
			bw_n, bw_br = check_la_kernel("bw_dot.s"); bf_n, bf_br = check_la_kernel("bw_dot.s", "LBB1_2:")
			say(f"[10] bw_dot.s bw_full inner loop: {bf_n} instructions (charged {M.C_FULL:.0f}) {'ok' if bf_n == 97 and bf_n <= M.C_FULL and bf_br == ['b.ne'] else 'MISMATCH'}")
			ok13 = fields_line()
			_t0 = time.time(); al_s = {p: alpha_contrib(p) for p in (2, 19, 23)}
			ok_alpha = op['alpha_f'] >= ALPHA_F_50K and op['alpha_g'] >= ALPHA_G_50K and DS.P.ALPHA_G >= ALPHA_G_50K and DS.P.ALPHA_F >= ALPHA_F_50K and all(round(v[0], 3) >= 0 for v in al_s.values())
			say(f"[14] Murphy alpha: computed over p <= {ALPHA_PMAX} ({ALPHA_F:.6f}, {ALPHA_G:.6f}) and over p < {ALPHA_G_PMAX} ({ALPHA_F_50K:.6f}, {ALPHA_G_50K:.6f}) <= charged ({op['alpha_f']}, {op['alpha_g']}) in collection and ({DS.P.ALPHA_F}, {DS.P.ALPHA_G}) in the descent {'ok' if ok_alpha else 'FAIL'}; "
				f"in-run exact checks at p=2/19/23 (f,g): " + ", ".join(f"{p}: ({v[0]:.6f}, {v[1]:.6f})" for p, v in sorted(al_s.items())) + f" ({time.time()-_t0:.1f}s)")
			checks.update(nbar=ok9, la_loop=la_n == 33 and bf_n == 97 and bf_br == ['b.ne'] and check_la_kernel('la_mac.s', 'LBB9_1:') == (11, ['b.ne']) and tail_ok and tail <= M.C_ROW and red9 <= M.C_RED9 and nu_ok and M.C_ENTRY >= 4 + 3 + M.UNIT_FRAC_CHARGED * 11 + (1 - M.UNIT_FRAC_CHARGED) * 33,
						  alpha_ok=ok_alpha, missing_logs=ok_rec, calibration=ok_cal, fields=ok13, trap_addback=ok_trap)
		tot = M.total_log2(e['log2_relcol'], la['log2_la'], cap, la['log2_gen'])
		say(f"ledger (log2 units): relation collection {e['log2_relcol']:.9f}, linear algebra {la['log2_la']:.9f}, generator {la['log2_gen']:.9f}, "
			f"descent cap {cap:g}, canonical list {M.LOG2_CANON:.0f}, Schirokauer maps {M.LOG2_SM:.0f}, "
			f"factor base {M.LOG2_FBC:.0f}, misc {M.LOG2_MISC:.0f}, verify {M.LOG2_VERIFY:.0f}, certificates 50")
		tot_decl = M.total_log2(DECLARED['relcol'][1], DECLARED['la'][1], cap, DECLARED['gen'][1])
		say(f"ledger total {tot:.9f} (unrounded rows); declared rows {DECLARED['relcol'][1]:.9f} + {DECLARED['la'][1]:.9f} + {cap:g} + {DECLARED['gen'][1]:.9f} + side budgets = {tot_decl:.9f}; claimed {CLAIM:.9f}: {'ok' if tot <= tot_decl <= CLAIM else 'EXCEEDS CLAIM'}; reserve {CLAIM - tot:.2e} bits")
		vals.update(ledger=tot, relcol=e['log2_relcol'], la=la['log2_la'], gen=la['log2_gen'], log2_N=la['log2_N'], log2_rels=e['log2_rels'], rels_2se=rels_2se, margin=e['margin'], work_se=work_se)
		checks.update(owner_schedule=O.check(M, la), rels_input=ok_rels_input and ok_relmeans, work_allowance=ok_allow, feasible=e['feasible'], within_claim=tot <= tot_decl <= CLAIM)
		extra.update(log2_relcol=e['log2_relcol'], log2_rels=e['log2_rels'], rels_minus_2se=rels_2se, margin_bits=e['margin'], work_se=work_se, work_allowance=M.CANON_WORK_ALLOWANCE, yield_se_log2=e['yield_se_log2'],
					 log2_N=la['log2_N'], log2_la=la['log2_la'], log2_gen=la['log2_gen'], total_log2W_full_ledger=tot, total_declared_rows=tot_decl, claimed=CLAIM)
	if part in ("core", "descent"):
		t0 = time.time()
		ok_counts = ball_counts_check()
		Pm = DS.P()
		dres = DS.evaluate(Pm, DS.PLAN16, DS.BOOT16)
		two_log = dres['two_log_units']
		bt = dres['boot']
		log2_EN, over, coll, pmin, mw, mrep, cheb = repeat_bound(dres)
		share = math.ceil((math.log2(math.comb(int(bt['nrels']), 2)) - math.log2(bt['balls'])) * 100) / 100
		success, markov = success_of(two_log, cap, bt, over, coll, mrep, share, la['log2_N'])
		rows = DS.report(dres)
		sv = [dres['DET'][k] for k in dres['DET'] if dres['DET'][k]['how'] == 'sieve']
		hmax = max(DS.CELLS[d['R']] for d in sv)
		ok_width = hmax < 1 << 103 and (hmax - 1) >> 40 < 1 << 63 and all(16 * max(d['nf'], d['ng']) < 1 << 16 for d in sv)
		mu_min = min(r['mu'] for r in rows); psi_max = max(max(r['psi'] for r in rows), bt['psi']); rmin_plan = min(r['R'] for r in rows)
		miss_f = [d['miss'] for k, d in dres['DET'].items() if k[1] == 'f']; miss_g = [d['miss'] for k, d in dres['DET'].items() if k[1] == 'g']
		say(f"[7] descent ({time.time()-t0:.1f}s): ball counts {'ok' if ok_counts else 'FAIL'}; {len(rows)} nodes (R = {rmin_plan}..{max(r['R'] for r in rows)}), min mu {math.floor(mu_min*10)/10:.1f}, max phi {max(r['phi'] for r in rows):.4f}, max psi {math.ceil(psi_max*1000)/1000:.3f}, retention floor {Pm.RET_FLOOR}; "
			f"missing ideals per relation (calibrated) f {min(miss_f):.2f}-{max(miss_f):.2f}, g {min(miss_g):.2f}-{max(miss_g):.2f}, boot {bt['miss']:.2f}; "
			f"boot n1={bt['n1']} nA={bt['nA']:.2f} kA={bt['kA']} R={bt['R']} relations={bt['nrels']} (m = {bt['m']:.0f} bits): 2^{bt['trials_log2']:.2f} trials at 2^{bt['trial_log2']:.2f} over 2^{math.log2(bt['balls']):.2f} lattices, "
			f"own 2^{bt['own_log2']:.2f} + subtree 2^{bt['sub_log2']:.2f} x retry {bt['retry']:.3f}; one log 2^{bt['one_log_log2']:.3f}; two-log work 2^{two_log:.4f} units; "
			f"cap 2^{cap:g}: Markov <= {markov:.6f}; P(boot fails) <= psi^{bt['nrels']} = {bt['psi']:.4f}^{bt['nrels']} < {math.ceil(bt['phi']*1e8)/1e8:.3e}; shared boot lattice <= C({bt['nrels']},2)/2^{math.log2(bt['balls']):.2f} = 2^{share:.2f}; "
			f"node counter: E[N] <= 2^{log2_EN:.2f} (counts16), N_MAX 2^{M.NODE_CAP_LOG2:g}, P(abort) <= 2^{over:.2f}; special-q repeats (smallest bin 2^{pmin:.2f} ideals) <= 2^{coll:.2f}; "
			f"missing-ideal repeats under the shared missing set (sum_I m_I wh_I^2 = 2^{mw:.2f}) <= 2^{mrep:.2f}; zero-occurrence mass: Chebyshev P(S_0 > (1 + 2^-8) E S_0) <= 2^{cheb:.2f}, charged 2^-30; all charged; success = {math.floor(success*1e6)/1e6:.6f} (floor); record widths {'ok' if ok_width else 'FAIL'}")
		cons = DS.miss_consistency(Pm)
		say(f"[7u] one ideal universe: the shared-set population sum_I w_I(u) m_I agrees with the recurrence's per-relation missing count to {cons:.1e} over the plan's norm sizes {'ok' if cons < 1e-9 else 'FAIL'}")
		ok_plan = mu_min >= Pm.MIN_RELS and psi_max <= Pm.PSI_MAX and max(r['phi'] for r in rows) <= 0.05 and coll <= -4 and rmin_plan >= 192 and all(r['how'] == 'sieve' for r in rows) and cons < 1e-9
		if part == "descent":
			import descent_witness as W0
			t0 = time.time()
			rc = W0.certify_forms(); ok_roots = rc['disjoint'] and rc['max_rel_radius'] < 2.0 ** -60
			say(f"[7r] root enclosures of the linear forms (Weierstrass disks, 60 digits): max relative radius {rc['max_rel_radius']:.1e}, min modulus {rc['min_modulus']:.3f}, disjoint {rc['disjoint']} {'ok' if ok_roots else 'FAIL'} ({time.time()-t0:.1f}s)")
			t0 = time.time()
			sens = []
			def variant(tag, relax=False, **kw):
				class PV(DS.P):
					pass
				for k_, v_ in kw.items():
					setattr(PV, k_, v_)
				if relax:
					PV.MIN_RELS = 1e-9; PV.PSI_MAX = 1.0 - 1e-12
				try:
					d_ = DS.evaluate(PV, DS.PLAN16, DS.BOOT16)
					b_ = d_['boot']; EN_, ov_, cl_, _, _, mr_, _ = repeat_bound(d_)
					s_, mk_ = success_of(d_['two_log_units'], cap, b_, ov_, cl_, mr_, share, la['log2_N'])
					viol = sum(1 for dd in d_['DET'].values() if dd['mu'] < Pm.MIN_RELS or dd['psi'] > Pm.PSI_MAX or dd['phi'] > 0.05)
					sens.append(f"{tag}: two-log 2^{d_['two_log_units']:.4f}, boot phi {b_['phi']:.2e}, boot psi {b_['psi']:.3f}, success at the shipped cap {s_:.4f}" + (f" (plan conditions violated at {viol} nodes)" if relax else ""))
				except ValueError as ex:
					if not relax:
						return variant(tag, relax=True, **kw)
					sens.append(f"{tag}: infeasible ({ex})")
			variant(f"estimate - 1 SE (LOG2_RELS {RELS_EST - YSE:.3f})", LOG2_RELS=RELS_EST - YSE)
			variant(f"point estimate (LOG2_RELS {RELS_EST:.3f})", LOG2_RELS=RELS_EST)
			variant("fitted NB shapes (r 0.6/1.86)", NB_R1=0.6, NB_R2=1.86)
			variant("UDUP 0.283 (smallest measured relation-weighted value)", UDUP=0.283)
			saved = DS.CALIB; DS.CALIB = tuple((lo, hi, 1.0) for lo, hi, c in saved); DS.miss_side.cache_clear()
			variant("c_r = 1 (every occurring ideal with lambda < 8 unknown)")
			DS.CALIB = saved; DS.miss_side.cache_clear()
			markov_max = 1 - 0.39 - (2 * bt['phi'] + 2.0 ** over + 2.0 ** coll + 2.0 ** mrep + 2.0 ** share + 2.0 ** -30 + 2.0 ** -48 + 2.0 ** -191 + 2.0 / R_ORDER)
			sens.append(f"headroom: E[W] may rise {math.log2(markov_max / markov):.3f} bits ({100 * (markov_max / markov - 1):.1f}%) before success reaches 0.39")
			say(f"[7s] sensitivity at the shipped plan and meter ({time.time()-t0:.1f}s): shipped: two-log 2^{two_log:.4f}, success {success:.4f}; " + "; ".join(sens))
			checks['roots'] = ok_roots
			extra['sensitivity'] = sens
		vals.update(two_log=two_log, markov=markov, log2_EN=log2_EN, over=over, rep=coll, mrep=mrep, cheb=cheb, boot_share=share, boot_phi=bt['phi'], success=success)
		checks.update(ball_counts=ok_counts, record_widths=ok_width, descent_plan_feasible=ok_plan, success_ok=success >= 0.39)
		extra.update(descent_two_logs_log2_units=two_log, descent_cap=cap, markov_abort=markov, log2_expected_nodes=log2_EN, node_cap_log2=M.NODE_CAP_LOG2,
					 counter_abort_log2=over, repeat_log2=coll, missing_repeat_log2=mrep, boot_share_log2=share, boot_phi=bt['phi'], success_exact=success, modelled_success_floor6=math.floor(success * 1e6) / 1e6,
					 descent_nodes=len(rows), descent_boot_log2=round(bt['one_log_log2'], 3), min_radius=rmin_plan)
	if part == "lattices":
		t0 = time.time()
		au = audit(31001 + seed)
		ok11 = au['boot']['within_model_bins'] and au['specialq_within_bins']
		say(f"[11] lattices ({time.time()-t0:.1f}s): certified sqrt(U)-weighted boot lattice ({au['boot']['swaps']} swaps), R = 96 cg max {max(s['cg_max'] for s in au['boot']['stats'])} "
			f"<= {max(x[1] for x in DS.P.BOOT_BINS):.1f} {'ok' if au['boot']['within_model_bins'] else 'FAIL'}; special-q lattices "
			+ "; ".join(f"{s['side']}@2^{s['bits']}: cf max {max(x['cf_max'] for x in s['stats'])}, cg max {max(x['cg_max'] for x in s['stats'])}" for s in au['specialq'])
			+ f" within bins (cf <= {max(b[1] for b in DS.P.QBINS):.1f}, cg <= {max(b[2] for b in DS.P.QBINS):.1f}) {'ok' if au['specialq_within_bins'] else 'FAIL'}; "
			f"unit mates f@2^260 ({au['unit_mates']['points']} points x {au['unit_mates']['units']} units, unconditioned): counts {au['unit_mates']['counts']}, "
			f"distinct/raw {au['unit_mates']['distinct_over_raw']} (charge {DS.P.UDUP}); largest line <= 2^{au['storage']['line_cells_log2_max']} cells; "
			f"chunk protocol {au['storage']['routing_stable']}")
		import ecm16 as EC
		t0 = time.time()
		toyecm = EC.toy(60, seed, sizes=((36, 1024, 50), (40, 2048, 50)))
		nf69 = math.ceil(72 * ((75 + math.log2(DS.CELLS[384])) / 12) + max(b[1] for b in DS.P.QBINS) - 69)
		sch = [EC.schedule(69.0, nf69), EC.schedule(DS.P.LX + DS.BOOT16[1] * (DS.BOOT16[0] - DS.P.LX) / DS.P.MESH, int(math.ceil(12 * ((DS.P.LOG2P6 + math.log2(DS.CELLS[DS.BOOT16[3]])) / 12) + DS.P.LOG2U6 + max(x[1] for x in DS.P.BOOT_BINS))))]
		ok12 = all(t['empirical'] >= t['model_corrected'] - 2 * math.sqrt(t['model_corrected'] * (1 - t['model_corrected']) / t['trials']) for t in toyecm)
		say(f"[12] ECM schedules ({time.time()-t0:.1f}s): " + "; ".join(
			f"{x['fbits']:.2f}-bit primes in {x['nbits']}-bit integers: B1 2^{math.log2(x['B1']):.0f}, B2 {x['B2'] // x['B1']} B1, {x['curves']} curves at p = {x['p_curve']:.2e}, 2^{x['log2_instructions']:.2f} instructions" for x in sch)
			+ "; toy detection (bits, trials, stage 1, stage 2, empirical, P1 raw/corrected): "
			+ "; ".join(f"({t['bits']}, {t['trials']}, {t['stage1']}, {t['stage2']}, {t['empirical']:.3f}, {t['model']:.3f}, {t['model_corrected']:.3f})" for t in toyecm)
			+ f" {'ok' if ok12 else 'FAIL'}")
		ok13 = fields_line()
		checks.update(audit=ok11, ecm=ok12, fields=ok13)
	if part == "retention":
		import descent_witness as W
		import retain16 as RT
		t0 = time.time(); st = load_stats(); RT.load_table(st)
		nb = st['nb_empty_window']
		rr = [(b, s_, R_) + W.retention(b, s_, R_, seed + b, 40) for b, s_, R_ in ((69, 'g', DS.PLAN16[(69, 'g')][4]), (150, 'g', DS.PLAN16[(150, 'g')][4]))]
		nb_ok = all(z <= zp for b_, s_, R_, lb_, m_, r_, z, zp in nb['rows'] if m_ >= 7.4)
		tab = st['prime_structure']['table']
		special = [p for p in RT.primes(2, 512) if any(e != 1 or any(m != 1 for m in mf + mg) for d1, e, mf, mg in tab[str(p)])]
		spot = sorted(set([p for i, p in enumerate(RT.primes(2, 512)) if i % 7 == seed % 7] + special))
		ok_tab = len(tab) == len(RT.primes(2, 512)) and all(RT.structure(p) == tab[str(p)] for p in spot)
		sample = {}
		for p in (2, 23):
			rs = W.residue_valuations(p, 300, seed + p)
			for side in 'fg':
				bound = 1 - RT.local_pdf(p, tab[str(p)], side, False, True)[0].get(0, 0.0)
				frac = 1 - rs[side].get(0, 0) / rs['n']
				sample[f"{p}{side}"] = (round(frac, 4), round(bound, 4), frac <= bound + 3 * math.sqrt(bound * (1 - bound) / rs['n']) + 1e-12)
		ok_sample = all(v[2] for v in sample.values())
		ok11b = nb_ok and all(x[3] for x in rr) and ok_tab and ok_sample
		say(f"[11b] retention inputs ({time.time()-t0:.1f}s): prime-ideal table recomputed at {len(spot)} primes {'ok' if ok_tab else 'FAIL'}; residue samples (300 base-primitive classes mod p^K, P(p | N) sample/bound): " + ", ".join(f"{k} {v[0]}/{v[1]}" for k, v in sample.items()) + f" {'ok' if ok_sample else 'FAIL'}"
			+ f"; base-primitive {RT.base_primitive():.4f}; certified cell bounds >= |N| at 40 cells: " + "; ".join(f"{b}{s_} R={R_}: {ok_}, max gap {g:.1e} bits" for b, s_, R_, ok_, g, k in rr)
			+ "; empty-window toys (bits, side, R, lb, mean, r, zeros, at r=0.5): " + "; ".join(str(tuple(x[:2]) + tuple(x[4:])) for x in nb['rows'])
			+ f" {'ok' if nb_ok else 'FAIL'}")
		checks.update(retention=ok11b)
	if part == "losses":
		import retain16 as RT
		t0 = time.time(); RT.load_table(load_stats())
		pf, ordf = RT.omitted_pdf('f', True); pg, ordg = RT.omitted_pdf('g', True)
		ex = lambda pdf: sum(x * v for x, v in pdf.items()) / RT.GRID
		counted = RT.counted_fraction()
		Pm = DS.P()
		t1 = time.time()
		rets = []
		cfmax = max(b[1] for b in Pm.QBINS); cgmax = max(b[2] for b in Pm.QBINS)
		for (lo, hi) in DS.bins_of(Pm):
			for side in ('f', 'g'):
				n, iA, kA, kB, R, thin = DS.PLAN16[(lo, side)]
				c = (hi + math.log2(DS.CELLS[R])) / 12   # the node's largest norms, as in descent16.node_cost
				nf = 72 * c + cfmax - (lo if side == 'f' else 0); ng = 12 * c + Pm.LOG2U6 + cgmax - (lo if side == 'g' else 0)
				ret, lf, lg = RT.retention(nf, ng, float(n), iA, kA, kB, Pm.SLACK)
				rets.append((lo, side, ret, lf, lg))
		ret_min = min(r[2] for r in rets); loss_max = max(max(r[3], r[4]) for r in rets)
		ok_floor = ret_min >= Pm.RET_FLOOR and counted >= Pm.COUNTED
		say(f"[11c] retention ({time.time()-t0:.1f}s): base-primitive {RT.base_primitive():.6f}, ordinary {ordf:.4f}/{ordg:.4f}, counted fraction {counted:.6f} >= pinned COUNTED {Pm.COUNTED} {'ok' if counted >= Pm.COUNTED else 'FAIL'}; "
			f"omitted mass means {ex(pf):.2f}/{ex(pg):.2f} bits, P(X > 96) {RT.cell_loss('f', 96):.2e}/{RT.cell_loss('g', 96):.2e}; exact retention at all {len(rets)} nodes at their largest norms ({time.time()-t1:.1f}s): min {ret_min:.6f} >= pinned RET_FLOOR {Pm.RET_FLOOR} {'ok' if ret_min >= Pm.RET_FLOOR else 'FAIL'}, largest one-side loss {loss_max:.2e}; per node (bin, side, retention, loss_f, loss_g): "
			+ "; ".join(f"({b}, {s_}, {r:.6f}, {lf:.1e}, {lg:.1e})" for b, s_, r, lf, lg in rets))
		checks.update(retention_floor=ok_floor)
	dec_ok = {k: (vals[k] <= v if op_ == '<=' else vals[k] >= v) for k, (op_, v) in DECLARED.items() if k in vals}
	if part in ("core", "ledger"):
		dec_ok['work_factor'] = corr['safe_work_population_factor'] <= F('1.029409690') and corr['safe_yield_factor'] >= F('0.969673713')  # certified products vs the stated 9-decimal bounds
	if dec_ok:
		say("declared bounds (unrounded value vs stated bound): " + "; ".join(f"{k} {vals[k]!r} {DECLARED[k][0]} {DECLARED[k][1]} {'ok' if dec_ok[k] else 'FAIL'}" for k in dec_ok if k in DECLARED) + (f"; count products {'ok' if dec_ok.get('work_factor', True) else 'FAIL'}" if 'work_factor' in dec_ok else ""))
		checks['declared_bounds'] = all(dec_ok.values())
	say(f"done in {time.time()-T0:.1f}s; part {part}; checks: {checks}")
	out = {"checks_ok": all(checks.values()), "part": part, "checks": checks, "declared": {k: v for k, v in DECLARED.items() if k in vals}}
	out.update(extra)
	print(json.dumps(out), flush=True)
if __name__ == "__main__":
	main()
