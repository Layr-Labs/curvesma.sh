"""stnfs declared experiment (bls12-381-fixed-g2-invert-g1). Stdlib only; stdout
is exactly one JSON object, progress on stderr. Sections [0]-[13]: field and
tower certificates, Dickman rho, pool regeneration, region metrics, yield/cost"""
import json, math, os, sys, time
from decimal import Decimal, localcontext
from fractions import Fraction as F
from random import Random
_HERE = os.path.dirname(os.path.abspath(__file__))
_DIR = os.path.join(_HERE, "files") if os.path.exists(os.path.join(_HERE, "files", "model16.py")) else _HERE
sys.path.insert(0, _DIR)
import common12 as C
import model16 as M
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
CLAIM = 122.248
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
	A2 = cfg['production_A_squared']
	index = {(r['seed'], r['index']): r for r in pool}
	n_can = n_prim = 0; bad = []
	t0 = time.time()
	for i in range(samples):
		a, b, aq, bq = regenerate(cfg, 1, i)
		can, work, witness = IC.canonical_a(a, cfg['matrix_a'], 100000)
		row = index.get((1, i))
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
	say(f"[2] seed-1 indices 0..{samples-1} regenerated ({time.time()-t0:.1f}s): {n_can} canonical, {n_prim} primitive; every stored "
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
def check_la_kernel(name="la_mac.s", start="LBB0_2:"):
	lines = open(os.path.join(_DIR, name)).read().splitlines()
	on, body = False, []
	for ln in lines:
		s = ln.strip()
		if s.startswith(start):
			on = True; continue
		if on:
			if s.startswith("LBB") or (s.startswith(";") and "bb." in s):
				break
			if s and not s.startswith(";") and not s.startswith("."):
				body.append(s.split()[0])
	return len(body), [m for m in body if m.startswith("b.") or m == "b"]
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
	"""[7] worst case under the meter: at most 2^(cap)/min own cost nodes in the two-log run, all placed in the
	smallest bin; two nodes draw the same ideal with probability <= 2^6/|bin| (1/p weighting within a 6-bit bin)."""
	nmax = M.DESCENT_CAP + math.log2(M.UNIT) - math.log2(min(d['own'] for d in dres['DET'].values()))
	pmin = min(math.log2(2 * (2 ** DS.log2li(b[1]) - 2 ** DS.log2li(b[0]))) for b in dres['bins'])
	return nmax, 2 * nmax - 1 + 6 - pmin
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
	B, swaps = W.lll(W.boot_basis(z))
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
			um = W.unit_mates(r['basis'], seed + 7, 96, 2, units)
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
def main():
	cfg_in = json.load(sys.stdin) if not sys.stdin.isatty() else {}
	samples = int(cfg_in.get("samples", 16)); seed = int(cfg_in.get("seed", 1))
	say(f"evidence: stnfs (bls12-381-fixed-g2-invert-g1); samples={samples} seed={seed}; unit {M.UNIT:.0f} instructions")
	ok0 = check_resultant(); ok0b = check_tower(); ok1 = check_rho_and_semismooth()
	cfg = json.load(open(os.path.join(_DIR, "config.json")))
	hdr, pool = M.load_pool(open(os.path.join(_DIR, "pool12-1.b85")).read() + open(os.path.join(_DIR, "pool12-2.b85")).read())
	n_total = sum(hdr['samples'].values())
	disks = EB.root_disks()
	ok2, n_can = verify_pool(cfg, pool, samples, disks)
	corr = LCB.correction(cfg)
	geo = M.Geometry(cfg, M.OP['theta'], F(corr['safe_yield_factor']), F(corr['safe_work_population_factor']))
	say(f"[3] eligible a 2^{math.log2(geo.va):.3f}, largest half-slice 2^{math.log2(geo.vb):.3f} cells; count factors "
		f"{corr['safe_yield_factor']} (yield), {corr['safe_work_population_factor']} (work); {len(pool)} canonical rows of {n_total} "
		f"({sum(r['primitive'] for r in pool)} primitive)")
	op = M.OP
	retain = Retention(int(op['lb']), int(op['lx2']), step=.25, max_bits=850, tilt=.05)
	cost = Retention(850, int(op['lx2']), step=.25, max_bits=851, tilt=.005)
	t0 = time.time()
	e = M.evaluate(pool, n_total, geo, op, retain, cost)
	ok_relmeans = e['nf_rel'] >= 69 * DS.P.MAIN_UF and e['ng_rel'] >= 69 * DS.P.MAIN_UG   # descent missing-log inputs
	rf = retain.at(e['nf_mean'], op['slack_f']); rg = retain.at(e['ng_mean'], op['slack_g'])
	say(f"[4] B=2^{op['lb']:.0f}, X2=2^{op['lx2']:.0f}, slacks ({op['slack_f']}, {op['slack_g']}), excess {op['excess']}, "
		f"width {op['width_bits']}, blocks of {op['L0']}, alpha ({op['alpha_f']}, {op['alpha_g']}); mean norms "
		f"{e['nf_mean']:.1f}/{e['ng_mean']:.1f} (relation-weighted {e['nf_rel']:.1f}/{e['ng_rel']:.1f}); retention at the means {rf['retention']:.4f}/{rg['retention']:.4f}; "
		f"no block selection (every block enumerated and charged): rels 2^{e['log2_rels']:.3f} vs need "
		f"2^{e['log2_need']:.3f}, margin {e['margin']:.3f} >= {M.MARGIN} ({'feasible' if e['feasible'] else 'INFEASIBLE'}); pool SE {e['yield_se_log2']:.3f} bits ({time.time()-t0:.1f}s)")
	la = M.linear_algebra(e['log2_U'], op['excess'])
	say(f"[5] population 2^{e['log2_population']:.3f} cells at {e['per_cell']:.1f} instructions per cell (routing {M.REC_ROUTE:.0f}/record, "
		f"{M.CELL_ROUTE:.0f}/cell; apply {M.C_APPLY:.2f}, recovery {M.C_RECOVER:.3f}): relation collection 2^{e['log2_relcol']:.4f} units; "
		f"N/U {la['ratio']:.4f} x pad {M.PAD_N} -> N = 2^{la['log2_N']:.4f}; linear algebra 2^{la['log2_la']:.4f} units ({M.C_ENTRY:.0f} instr/entry)")
	t0 = time.time()
	Us, res = filtsim.run(op['lb'], e['nf_mean'] / op['lb'], e['ng_mean'] / op['lb'], op['excess'],
						 seed, 12, kmax=64, targets=(100,), verbose=False)
	r100 = res[100]
	say(f"[6] filtering (scale 2^12, seed {seed}, {time.time()-t0:.1f}s): N/U = {1/r100[0]:.4f} at weight {r100[3]:.0f} (max |coefficient| {r100[2]}); "
		"embedded 2^17 runs (excess: seed, U/N, N, max coef): " + "; ".join(f"{ex}: " + ", ".join(f"({s}, {un:.2f}, {n}, {mc})" for (s, un, n, mc, w) in v)
								 for ex, v in sorted(M.FILT_TABLE.items())))
	t0 = time.time()
	rec = filtsim.missing_recon(seed)
	ok_rec = all(mf <= md for (lo, hi, n, mf, md, h3) in rec)
	say(f"[6b] missing logs after purge and back-substitution (scale 2^12, seed {seed}, {time.time()-t0:.1f}s) per lambda band "
		"(ideals, measured, model P(Poisson <= 2), lost with >= 3 occurrences): "
		+ "; ".join(f"[{lo:g},{hi:g}): {n}, {mf:.4f}, {md:.4f}, {h3}" for (lo, hi, n, mf, md, h3) in rec) + f" {'ok' if ok_rec else 'FAIL'}")
	t0 = time.time()
	ok_counts = ball_counts_check()
	import retain16 as RT0
	RT0.load_table(load_stats())
	Pm = DS.P()
	dres = DS.evaluate(Pm, DS.PLAN16, DS.BOOT16)
	two_log = dres['two_log_units']
	bt = dres['boot']
	cap = M.DESCENT_CAP
	markov = 2 ** (two_log - cap)
	nmax, coll = repeat_bound(dres); REP = 2.0 ** (math.ceil(coll * 100) / 100)
	p_desc = 1 - 2 * bt['phi'] - markov - REP
	p_main = 1.0 - 2.0 ** (la['log2_N'] + 2 - 254)   # under H3 the main phase fails only by the Wiedemann projections (proof.md section 5)
	success = p_main * p_desc
	rows = DS.report(dres)
	sv = [dres['DET'][k] for k in dres['DET'] if dres['DET'][k]['how'] == 'sieve']
	hmax = max(DS.CELLS[d['R']] for d in sv)
	ok_width = hmax < 1 << 103 and (hmax - 1) >> 40 < 1 << 63 and all(16 * max(d['nf'], d['ng']) < 1 << 16 for d in sv)
	mu_min = min(r['mu'] for r in rows); psi_max = max(max(r['psi'] for r in rows), bt['psi'])
	say(f"[7] descent ({time.time()-t0:.1f}s): ball counts {'ok' if ok_counts else 'FAIL'}; {len(rows)} nodes, min mu {math.floor(mu_min*10)/10:.1f}, max psi {math.ceil(psi_max*1000)/1000:.3f}; "
		f"boot n1={bt['n1']} nA={bt['nA']:.2f} kA={bt['kA']} R={bt['R']}: 2^{bt['trials_log2']:.2f} trials at 2^{bt['trial_log2']:.2f} (pr 2^{bt['pr_log2']:.1f}), "
		f"own 2^{bt['own_log2']:.2f} + subtree 2^{bt['sub_log2']:.2f} x retry {bt['retry']:.3f}; one log 2^{bt['one_log_log2']:.3f}; two-log work 2^{two_log:.3f} units; "
		f"cap 2^{cap:g}: Markov <= 2^{math.ceil((two_log-cap)*1000)/1000:.3f}; P(boot fails) < {math.ceil(bt['phi']*1e7)/1e7:.2e}; special-q repeats: <= 2^{nmax:.2f} nodes under the cap, P <= 2^{math.ceil(coll*100)/100:.2f}, charged; success = {success:.4f}; record widths: half count <= "
		f"2^{math.log2(hmax):.2f}, owner < 2^{math.log2((hmax - 1) >> 40):.2f}, max score {16 * max(max(d['nf'], d['ng']) for d in sv):.0f} {'ok' if ok_width else 'FAIL'}")
	nb, np_ = nbar(8, 12)
	say(f"[9] degree-1 ideals per prime, both sides, {np_} primes in [2^8, 2^12): {nb:.3f} (model NBAR = {M.NBAR})")
	ok9 = nb <= M.NBAR * 1.01
	la_n, la_br = check_la_kernel()
	say(f"[10] la_mac.s inner loop: {la_n} instructions, branches {la_br} {'ok' if la_n == 33 and la_br == ['b.ne'] else 'MISMATCH'}")
	bw_n, bw_br = check_la_kernel("bw_dot.s"); bf_n, bf_br = check_la_kernel("bw_dot.s", "LBB1_2:")
	say(f"[10] bw_dot.s bw_full (full-field multiply-accumulate) inner loop: {bf_n} instructions (charged {M.C_FULL:.0f}) {'ok' if bf_n == 97 and bf_n <= M.C_FULL and bf_br == ['b.ne'] else 'MISMATCH'}")
	t0 = time.time()
	au = audit()
	ok11 = au['boot']['within_model_bins'] and au['specialq_within_bins']
	say(f"[11] audit ({time.time()-t0:.1f}s): certified boot lattice ({au['boot']['swaps']} swaps), boot cg max {max(s['cg_max'] for s in au['boot']['stats'])} "
		f"<= {max(x[1] for x in DS.P.BOOT_BINS):.1f} {'ok' if au['boot']['within_model_bins'] else 'FAIL'}; special-q lattices "
		+ "; ".join(f"{s['side']}@2^{s['bits']}: cf max {max(x['cf_max'] for x in s['stats'])}, cg max {max(x['cg_max'] for x in s['stats'])}" for s in au['specialq'])
		+ f" within bins (cf <= {max(b[1] for b in DS.P.QBINS):.1f}, cg <= {max(b[2] for b in DS.P.QBINS):.1f}) {'ok' if au['specialq_within_bins'] else 'FAIL'}; "
		f"unit mates f@2^260 ({au['unit_mates']['points']} points x {au['unit_mates']['units']} units): counts {au['unit_mates']['counts']}, "
		f"distinct/raw {au['unit_mates']['distinct_over_raw']} (model {DS.P.UDUP}); largest line <= 2^{au['storage']['line_cells_log2_max']} cells; "
		f"chunk protocol {au['storage']['routing_stable']}")
	import descent_witness as W
	import retain16 as RT
	import ecm16 as EC
	t0 = time.time(); st = load_stats()
	nb = st['nb_empty_window']
	rr = [(b, s_, R_) + W.retention(b, s_, R_, seed + b, 40) for b, s_, R_ in [(b, s_, dres['DET'][(bb, s_)]['R']) for b, s_ in ((69, 'g'), (150, 'g')) for bb in dres['bins'] if bb[0] == b]]
	nb_ok = all(z <= zp for b_, s_, R_, lb_, m_, r_, z, zp in nb['rows'] if m_ >= 7.4)
	tab = st['prime_structure']['table']
	special = [p for p in RT.primes(2, 512) if any(e != 1 or any(m != 1 for m in mf + mg) for d1, e, mf, mg in tab[str(p)])]
	spot = sorted(set([p for i, p in enumerate(RT.primes(2, 512)) if i % 7 == seed % 7] + special))
	ok_tab = len(tab) == len(RT.primes(2, 512)) and all(RT.structure(p) == tab[str(p)] for p in spot)
	sample = {}
	for p in (2, 23):
		rs = W.residue_valuations(p, 400, seed + p)
		for side in 'fg':
			bound = 1 - RT.local_pdf(p, tab[str(p)], side, False, True).get(0, 0.0)
			frac = 1 - rs[side].get(0, 0) / rs['n']
			sample[f"{p}{side}"] = (round(frac, 4), round(bound, 4), frac <= bound + 3 * math.sqrt(bound * (1 - bound) / rs['n']) + 1e-12)
	ok_sample = all(v[2] for v in sample.values())
	pf, ordf = RT.omitted_pdf('f', True); pg, ordg = RT.omitted_pdf('g', True)
	ex = lambda pdf: sum(x * v for x, v in pdf.items()) / RT.GRID
	sieve_rows = [r for r in rows if r['how'] == 'sieve']
	ret_min = min(r['ret'] for r in sieve_rows); loss_max = max(max(v['loss_f'], v['loss_g']) for v in dres['DET'].values())
	ok11b = nb_ok and all(x[3] for x in rr) and ok_tab and ok_sample
	say(f"[11b] retention ({time.time()-t0:.1f}s): prime-ideal table (degree, exponent, root multiplicities of f, g) recomputed at {len(spot)} primes {'ok' if ok_tab else 'FAIL'}; "
		"repeated roots or ramification at " + "; ".join(f"{p}: {tab[str(p)]}" for p in special)
		+ "; residue-class samples (400 base-primitive classes mod p^K, P(p | N) sample/bound): " + ", ".join(f"{k} {v[0]}/{v[1]}" for k, v in sample.items()) + f" {'ok' if ok_sample else 'FAIL'}"
		+ f"; base-primitive {RT.base_primitive():.4f}, ordinary {ordf:.4f}/{ordg:.4f}; omitted mass means {ex(pf):.2f}/{ex(pg):.2f} bits, "
		f"P(X > 96) {RT.cell_loss('f', 96):.2e}/{RT.cell_loss('g', 96):.2e}; loss <= {loss_max:.1e}, retention >= {ret_min:.4f}; bounds >= |N|: " + "; ".join(f"{b}{s_} R={R_}: {ok_}, max gap {g:.1e} bits" for b, s_, R_, ok_, g, k in rr)
		+ "; empty-window toys (bits, side, R, lb, mean, r, zeros, zeros at r=0.5): " + "; ".join(str(tuple(x)) for x in nb['rows'])
		+ f" {'ok' if nb_ok else 'FAIL'}")
	t0 = time.time()
	toyecm = EC.toy(60, seed, sizes=((36, 1024, 50), (40, 2048, 50)))
	sch = [EC.schedule(69.0, int(math.ceil(dres['DET'][((69, 75), 'f')]['nf'])) ), EC.schedule(bt['nA'], int(math.ceil(bt['m'])))]
	ok12 = all(t['empirical'] >= t['model_corrected'] - 2 * math.sqrt(t['model_corrected'] * (1 - t['model_corrected']) / t['trials']) for t in toyecm)
	say(f"[12] ECM schedules ({time.time()-t0:.1f}s): " + "; ".join(
		f"{x['fbits']:.2f}-bit primes in {x['nbits']}-bit integers: B1 2^{math.log2(x['B1']):.0f}, B2 {x['B2'] // x['B1']} B1, wheel {x['wheel']}, "
		f"{x['curves']} curves at p = {x['p_curve']:.2e}, {x['stage1_M']:.3g} + {x['stage2_M']:.3g} M per curve, 2^{x['log2_instructions']:.2f} instructions" for x in sch)
		+ "; toy detection (bits, trials, stage 1, stage 2, empirical, P1 without/with the power correction): "
		+ "; ".join(f"({t['bits']}, {t['trials']}, {t['stage1']}, {t['stage2']}, {t['empirical']:.3f}, {t['model']:.3f}, {t['model_corrected']:.3f})" for t in toyecm)
		+ f" {'ok' if ok12 else 'FAIL'}")
	t0 = time.time(); fl = check_fields()
	say(f"[13] number fields ({time.time()-t0:.1f}s): N(f(X)) of degree {fl['deg_f']} and N(g(X)) of degree {fl['deg_g']} squarefree mod r {'ok' if fl['unramified'] else 'FAIL'} (r unramified in K_f, K_g); "
		f"disc h = {fl['disc_h']}; signatures K_f {fl['sig']['f']}, K_f+ {fl['sig']['fp']}, K_g {fl['sig']['g']}: odd unit rank {fl['sig']['odd']} = SM columns {M.N_SM} {'ok' if fl['sig']['odd'] == M.N_SM else 'FAIL'}")
	ok13 = fl['unramified'] and fl['sig']['odd'] == M.N_SM
	tot = M.total_log2(e['log2_relcol'], la['log2_la'], cap, la['log2_gen'])
	say(f"ledger (log2 units): relation collection {e['log2_relcol']:.4f}, linear algebra {la['log2_la']:.4f}, generator {la['log2_gen']:.2f}, "
		f"descent cap {cap:g} (expected {two_log:.3f}), canonical list {M.LOG2_CANON:.0f}, Schirokauer maps {M.LOG2_SM:.0f}, "
		f"factor base {M.LOG2_FBC:.0f}, misc {M.LOG2_MISC:.0f}, verify {M.LOG2_VERIFY:.0f}")
	say(f"total_log2W_full_ledger = {tot:.6f} (claimed {CLAIM:.6f}: {'ok' if tot <= CLAIM else 'EXCEEDS CLAIM'}); "
		f"reserve {CLAIM - tot:.3f} bits")
	checks = dict(resultant=ok0, tower=ok0b, rho_semismooth=ok1, pool=ok2, nbar=ok9, la_loop=la_n == 33 and bf_n == 97 and bf_br == ['b.ne'],
				 ball_counts=ok_counts and ok_relmeans, missing_logs=ok_rec, record_widths=ok_width, descent_plan_feasible=mu_min >= Pm.MIN_RELS and psi_max <= Pm.PSI_MAX and max(r['phi'] for r in rows) <= 0.05 and coll <= -4,
				 audit=ok11, retention=ok11b, ecm=ok12, fields=ok13, feasible=e['feasible'], within_claim=tot <= CLAIM, success_ok=success >= 0.39)
	say(f"done in {time.time()-T0:.1f}s; checks: {checks}")
	print(json.dumps({"checks_ok": all(checks.values()), "checks": checks,
					 "total_log2W_full_ledger": round(tot, 6), "claimed": CLAIM,
					 "log2_relcol": round(e['log2_relcol'], 4), "log2_rels": round(e['log2_rels'], 4),
					 "margin_bits": round(e['margin'], 4), "log2_N": round(la['log2_N'], 4),
					 "log2_la": round(la['log2_la'], 4), "descent_two_logs_log2_units": round(two_log, 3),
					 "descent_cap": cap, "markov_abort": markov, "modelled_success": round(success, 4),
					 "descent_nodes": len(rows), "descent_boot_log2": round(bt['one_log_log2'], 3)}), flush=True)
if __name__ == "__main__":
	main()
