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
import filtsim
P = 4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787
U = C.U
H = C.H
LN2 = math.log(2)
CLAIM = 120.254755005279
import retain16
ALPHA_PMAX = 2000
ALPHA_F, ALPHA_G = 4.868096, 0.843643 # via --alpha (p=2 double-root term)
ALPHA_G_PMAX = 50000
ALPHA_G_50K = 0.731543 # 50k sweep (judge bound 0.73154)
ALPHA_F_PMAX = 50000
ALPHA_F_50K = 5.085035 # 50k sweep (exact Kummer-Dedekind)
def alpha_contrib(p):
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
def alpha_sweep(pmax):
	af = ag = 0.0
	for p in retain16.primes(2, pmax):
		a, b = alpha_contrib(p); af += a; ag += b
	return af, ag
DECLARED = {'ledger': ('<=', 120.254755005279), 'relcol': ('<=', 119.9317), 'la': ('<=', 117.7346), 'log2_rels': ('>=', 63.6890),
	'margin': ('>=', 0.1002), 'work_se': ('<=', 0.00825), 'two_log': ('<=', 113.158), 'markov': ('<=', 0.279),
	'log2_EN': ('<=', 17.674), 'over': ('<=', -9.62), 'rep': ('<=', -10.70), 'boot_phi': ('<=', 1.0e-04),
	'success': ('>=', 0.479)} # v22s: S* per-node exceptional mass in the recurrence; floor 0.03, BOOT_RELS 24
_REPORT = []
T0 = time.time()
def load_stats():
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
def verify_pool(cfg, pool, samples, seed, sample_counts, anchors):
	index = {(r['seed'], r['index']): r for r in pool}
	vrng = Random(f"pool-verify:{seed}")
	picks = [(1, i) for i in vrng.sample(range(int(sample_counts['1'])), samples // 2)]
	picks += [(2, i) for i in vrng.sample(range(int(sample_counts['2'])), samples - samples // 2)]
	n_can = n_prim = 0; bad = []
	t0 = time.time()
	for sd, i in picks:
		a, b, aq, bq = C.regenerate_pair(cfg, sd, i)
		try:
			can, work, witness = IC.canonical_a(a, cfg['matrix_a'], 100000)
		except RuntimeError:
			can, work, witness = IC.canonical_a(a, cfg['matrix_a'], 1000000) # sampler fallback (pool header records the retry)
		row = index.get((sd, i))
		if can != (row is not None):
			bad.append(('canonical flag', sd, i, can)); continue
		if not can:
			continue
		n_can += 1
		n_prim += row['primitive']
	n_anch = 0
	for a in anchors:
		row = pool[a['position']]
		if not (round(row['qa'] * 65536) == a['qa16'] and round(row['nf'] * 16) == a['nf16']
				and round(row['ng'] * 16) == a['ng16'] and row['primitive'] == bool(a['flags'] & 1)
				and row['rejected'] == bool(a['flags'] & 2) and not row['nonfinite']
				and row['evals'] == a['evals'] and row['leaves'] == a['leaves']
				and row['blen'] == a['blen'] and row['rlen'] == 16 * a['rlen16']):
			bad.append(('anchor fields', a['position'])); continue
		n_anch += 1
	say(f"[2] {samples} membership picks across seeds 1-2 plus {len(anchors)} field anchors ({time.time()-t0:.1f}s): {n_can} canonical, {n_prim} primitive, "
		f"{n_anch} anchors match; list integrity and regeneration {'ok' if not bad else 'FAIL ' + str(bad[:3])}")
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
		for s in sorted(set(roots_modp(H, p, rng))):
			tot += poly_roots_count([1 + s - s * s - s ** 3, 0, s ** 5, 0, 0, 0, 2 * s ** 3, 0, 0, 0, -2 * s, 0, 1], p) + \
				poly_roots_count([U * s, 0, 1], p)
	return tot / len(ps), len(ps)
def _asm(name):
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
	f = {k: check_la_func(k) for k in ("_send_pair", "_combine", "_finish", "_reduce9", "_fold9", "_reduce6")}
	ok = all(v[1] for v in f.values()) and f["_finish"][2] == ["_fold9", "_reduce6"] and f["_reduce9"][2] == ["_fold9", "_reduce6"] \
		and not any(f[k][2] for k in ("_send_pair", "_combine", "_fold9", "_reduce6"))
	sub = f["_fold9"][0] + f["_reduce6"][0]
	row = 9 * (f["_send_pair"][0] + 18) + 9 * (f["_combine"][0] + 2) + f["_finish"][0] + sub + 2 + 100 + 40 + 168
	red9 = f["_reduce9"][0] + sub + 2
	return row, red9, {k: v[0] for k, v in f.items()}, ok
def ball_counts_check():
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
	import counts16
	EN, _ = counts16.expected_nodes(DS.P(), dres)
	pmin = min(math.log2(2 * (2 ** DS.log2li(b[1]) - 2 ** DS.log2li(b[0]))) for b in dres['bins'])
	over = math.ceil((math.log2(EN) - M.NODE_CAP_LOG2) * 100) / 100
	coll = math.ceil((2 * M.NODE_CAP_LOG2 - 1 + 6 - pmin) * 100) / 100
	return math.log2(EN), over, coll, pmin
def audit(rng_seed=31001, npr=None, bpr=None):
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
	B, swaps = W.boot_lll(W.boot_basis(z))
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
	units = W.gen_units(8)
	assert len(units) == 2442
	n1, iAb, kAb, _ = DS.BOOT16
	sb = sum(DS.ss_lower2(lng, Pm.LX, float(n1), iAb, kAb, 1, Pm.MESH, Pm.ALPHA_G, Pm.TRUNC) / (1 + m)
		for _, lng, m in W.unit_mates(B, rng_seed + 33, 4, 96, units, boot=True)['triples'])
	srng = Random(rng_seed + 77)
	dp = [1] + [0] * 16
	for _ in range(12):
		ndp = [0] * 17
		for s in range(17):
			for v in range(-4, 5):
				if s + v * v <= 16:
					ndp[s + v * v] += dp[s]
		dp = ndp
	tot = sum(dp[1:]); sh = [0] * 17; okshape = True
	for _ in range(600):
		c = W.int_ball_point(srng, 4)
		n2 = sum(x * x for x in c)
		okshape &= 0 < n2 <= 16
		for x in reversed(c):
			if x:
				okshape &= x > 0
				break
		sh[n2] += 1
	chi2 = sum((sh[s] - 600 * dp[s] / tot) ** 2 / (600 * dp[s] / tot) for s in range(1, 17))
	return dict(boot=dict(swaps=swaps, stats=[{k: (round(v, 2) if isinstance(v, float) else v) for k, v in s.items()} for s in bst],
						 within_model_bins=all(s['cg_max'] <= boot_bins_max for s in bst)),
				specialq=sq, specialq_within_bins=env_ok,
				relw=dict(boot_yield_ratio=(sb / 96) / 2.0 ** bpr),
				sampler=dict(chi2=chi2, ok=bool(okshape) and chi2 < 60),
				storage=store16.check())
def check_fields():
	import retain16 as RT, descent_witness as W
	R = W.R; Hr = [x % R for x in H]
	def res(a, b):
		a = RT._trim([x % R for x in a]); b = RT._trim([x % R for x in b]); s = 1
		while b and len(b) > 1:
			da, db = len(a) - 1, len(b) - 1; rem = RT.pmod(a, b, R)
			if not rem:
				return 0
			s = s * (-1) ** (da * db) * pow(b[-1], da - (len(rem) - 1), R) % R; a, b = b, rem
		return s * pow(b[0], len(a) - 1, R) % R if b else 0
	def norm_poly(terms, deg):
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
	out['disc_h'] = C.norm([-1, 2, 3, 0, 0, 6]) # |disc h| = |N(h'(y))|
	out['sig'] = signatures()
	return out
def signatures():
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
def sstar(dres, st):
	import itertools as it, descent_witness as W
	S = st['sstar']; Pm = DS.P(); bins, D, PH, bt = dres['bins'], dres['DET'], dres['PHI'], dres['boot']; X = 1e-9; LX = Pm.LX
	pmx = lambda sd, n: max([PH[(b, sd)] for b in bins if b[0] < n - X] + [0.0])
	ch, todo, seen = {}, [(bt['nA'], LX, 'g'), (bt['n1'], bt['nA'], 'g')], set()
	while todo: # reachability from the boot; charge = min over parent slots
		top, lo, sides = todo.pop()
		for k in [(b, sd) for b in bins for sd in sides if b[0] < top - X and b[1] > lo + X]:
			ch[k] = min(ch.get(k, 1.0), pmx(k[1], top))
			if k not in seen:
				seen.add(k); d = D[k]
				if d['kA'] and d['nA'] > LX + X: todo.append((d['nA'], LX, 'fg'))
				if d['kB'] and d['n'] > d['nA'] + X: todo.append((d['n'], d['nA'], 'fg'))
	b85 = next(b for b in bins if b[0] == 85); pm = max(PH[k] for k in PH if k[0][0] <= Pm.THIN_MAX)
	cl = {'%d%s' % (k[0][0], k[1]): (D[k], c) for k, c in ch.items()}; cl.update({'comp' + sd: (D[(b85, sd)], pm) for sd in 'fg'})
	rows = {r[0]: r[1:] for r in S['rows']}; a = S['alpha'] / len(rows); res = []
	for c, (d, g) in cl.items():
		n, K, b, m, k = rows[c]; lo, hi = 0.0, 1.0
		for _ in range(80): # Clopper-Pearson upper, level 1 - a
			p = (lo + hi) / 2
			if sum(math.comb(n, i) * p ** i * (1 - p) ** (n - i) for i in range(k + 1)) > a: lo = p
			else: hi = p
		s = d['mu'] * (1 - d['psi']); r = Pm.NB_R1 if (d['kA'] if d['nA'] > LX + X else 0) + (d['kB'] if d['n'] > d['nA'] + X else 0) <= 1 else Pm.NB_R2
		res.append(max((g - (1 + x * s / r) ** -r - hi / (1 - math.exp(-3 * K * (m - x) ** 2 / (2 * b * (2 * x + m)))), hi, x, c) for x in [m * j / 400 for j in range(4, 400)]))
	U = ([0, -1, 0, 0, 0, 0], [1, 0, -1, -1, 1, -1], [-2, 0, 1, 1, -1, 1]), ([1, -1, -1, 0, 0, -1], [0, -1, -1, -1, 0, -1], [0, -1, -1, 1, -1, 0])
	z = W._croots(list(C.H)); z = sorted([q for q in z if abs(q.imag) < X], key=lambda q: q.real) + sorted([q for q in z if q.imag > X], key=lambda q: q.real)
	ev = lambda u, q: sum(u[j] * q ** j for j in range(6)); L = [[math.log(abs(ev(u, q))) for u in U[0]] for q in z]; err = 0.0; M = 2 ** 64
	d3 = lambda Q: Q[0][0] * (Q[1][1] * Q[2][2] - Q[1][2] * Q[2][1]) - Q[0][1] * (Q[1][0] * Q[2][2] - Q[1][2] * Q[2][0]) + Q[0][2] * (Q[1][0] * Q[2][1] - Q[1][1] * Q[2][0])
	for c, N, s0, t0, x, ys in S['replay']: # campaign replay: lattice, splitmix64 cells, complete orbits
		d = cl[c][0]; R = d['R']; w = [pow(s0, j, N) for j in range(6)] + [-t0 * pow(s0, j, N) % N for j in range(6)]
		B, _ = W.lll([[N] + [0] * 11] + [[-w[j] % N] + [int(i == j) for i in range(1, 12)] for j in range(1, 12)])
		Bi = C.inverse(B); lq = math.log2(N); Mh = [[R * math.sqrt(sum(abs(ev(B[i][6 * h:6 * h + 6], q)) ** 2 for i in range(12))) for q in z] for h in (0, 1)]
		for y0 in ys:
			cc = []
			while not (any(cc) and sum(u * u for u in cc) <= R * R):
				cc = []
				while len(cc) < 12:
					x = (x + 0x9E3779B97F4A7C15) % M; y = (x ^ x >> 30) * 0xBF58476D1CE4E5B9 % M; y = (y ^ y >> 27) * 0x94D049BB133111EB % M; y = (y ^ y >> 31) & (1 << (2 * R).bit_length()) - 1
					if y <= 2 * R: cc.append(y - R)
			if [u for u in cc if u][-1] < 0: cc = [-u for u in cc]
			v = [sum(cc[i] * B[i][j] for i in range(12)) for j in range(12)]
			sv = math.prod(DS.ss_lower2(math.ceil(math.log2(nn) - lq * (c[-1] == h) + X), LX, d['n'], d['iA'], d['kA'], d['kB'], Pm.MESH, al, Pm.TRUNC) for nn, h, al in zip(C.norms(v[:6], v[6:]), 'fg', (Pm.ALPHA_F, Pm.ALPHA_G)))
			t = [math.log(min(Mh[0][k] / abs(ev(v[:6], q)), Mh[1][k] / abs(ev(v[6:], q))) * 1.02) + X for k, q in enumerate(z)]
			V = [[d3([[t[i] if j == e else L[i][j] for j in range(3)] for i in range(4) if i != q]) / d3([L[i] for i in range(4) if i != q]) for e in range(3)] for q in range(4)]
			m2 = 0
			for e in it.product(*[range(math.floor(min(p[j] for p in V)) - 1, math.ceil(max(p[j] for p in V)) + 2) for j in range(3)]):
				if any(e) and all(sum(L[k][j] * e[j] for j in range(3)) <= t[k] + 1e-7 for k in range(4)):
					u = [1, 0, 0, 0, 0, 0]
					for j in range(3):
						for _ in range(abs(e[j])): u = C.mul(u, U[e[j] < 0][j])
					y = [sum(q * Bi[i][j] for i, q in enumerate(C.mul(u, v[:6]) + C.mul(u, v[6:]))) for j in range(12)]
					m2 += all(q.denominator == 1 for q in y) and sum(q * q for q in y) <= R * R
			err = max(err, abs(sv / (1 + m2) / (2.0 ** d['pr_log2'] * Pm.UDUP) / y0 - 1))
	return dict(same=set(cl) == set(rows), n=len(rows), lat=sum(r[0] for r in rows.values()), k=sum(r[4] for r in rows.values()), worst=min(res), err=err, nrep=sum(len(r[5]) for r in S['replay']))
def main():
	cfg_in = json.load(sys.stdin) if not sys.stdin.isatty() else {}
	samples = int(cfg_in.get("samples", 16)); seed = int(cfg_in.get("seed", 1))
	say(f"evidence: stnfs (bls12-381-fixed-g1-invert-g2); samples={samples} seed={seed}; unit {M.UNIT:.0f} instructions")
	ok0 = check_resultant(); ok0b = check_tower(); ok1 = check_rho_and_semismooth()
	cfg = json.load(open(os.path.join(_DIR, "config.json")))
	hdr, pool = M.load_pool(open(os.path.join(_DIR, "pool12.b85")).read(), cfg)
	n_total = sum(hdr['samples'].values())
	ok2, n_can = verify_pool(cfg, pool, samples, seed, hdr['samples'], hdr['anchors'])
	theta_pool = F(cfg['production_theta']); theta_reg = F(cfg['production_theta_region'])
	assert theta_reg <= theta_pool == M.OP['theta']
	pool_reg = [r for r in pool if r['qa'] <= float(theta_reg)]
	n_eff = n_total * float((theta_reg / theta_pool) ** 3)
	rmin = math.isqrt(int(theta_reg * cfg['production_A_squared']))
	corr = LCB.correction(cfg, rmin)
	geo = M.Geometry(cfg, theta_reg, F(corr['safe_yield_factor']), F(corr['safe_work_population_factor']))
	say(f"[3] region theta={theta_reg} (pool theta={theta_pool}), a-radius {rmin}; "
		f"eligible a 2^{math.log2(geo.va):.3f}, half-slice 2^{math.log2(geo.vb):.3f} cells; certified count products "
		f"{float(corr['safe_yield_factor']):.9f} (yield), {float(corr['safe_work_population_factor']):.9f} (work) at every radius >= {rmin}; "
		f"{len(pool_reg)}/{len(pool)} canonical rows in region ({sum(r['primitive'] for r in pool_reg)} primitive), "
		f"effective sample {n_eff:.3f} of {n_total} draws")
	op = M.OP
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
	say(f"[4a] work estimator: {ncan}/{n_total} contributing draws; row SE {100*cond_se:.3f}%, "
		f"incidence {100*math.sqrt((1-pcan)/ncan):.3f}%, total SE {100*work_se:.3f}%; allowance "
		f"{M.CANON_WORK_ALLOWANCE:.4f} >= 1 + {M.WORK_SE_K:.2f} SE = {allow_need:.5f} {'ok' if ok_allow else 'FAIL'}")
	ok_relmeans = e['nf_rel'] >= 69 * DS.P.MAIN_UF and e['ng_rel'] >= 69 * DS.P.MAIN_UG # descent missing-log inputs
	ok_rels_input = e['log2_rels'] >= DS.P.LOG2_RELS
	rf = retain.at(e['nf_mean'], op['slack_f']); rg = retain.at(e['ng_mean'], op['slack_g'])
	say(f"[4] B=2^{op['lb']:.0f}, X2=2^{op['lx2']:.0f}, slacks ({op['slack_f']}, {op['slack_g']}), excess {op['excess']}, "
		f"width {op['width_bits']}, blocks {op['L0']}, alpha ({op['alpha_f']}, {op['alpha_g']}); mean norms "
		f"{e['nf_mean']:.1f}/{e['ng_mean']:.1f} (rel-weighted {e['nf_rel']:.1f}/{e['ng_rel']:.1f}); retention {rf['retention']:.4f}/{rg['retention']:.4f}; "
		f"no block selection: rels 2^{e['log2_rels']:.3f} vs need "
		f"2^{e['log2_need']:.3f}, margin {math.floor(e['margin']*1000)/1000:.3f} >= max({M.MARGIN_FLOOR}, {M.MARGIN_SE:g} x SE {e['yield_se_log2']:.4f}) = {e['margin_need']:.4f} ({'feasible' if e['feasible'] else 'INFEASIBLE'}); "
		f"central-core cut qa < {M.Q_CORE} drops {e['n_core_cut']} rows from the yield (cells still charged); "
		f"descent LOG2_RELS input {DS.P.LOG2_RELS:.3f} <= supplied {e['log2_rels']:.4f} {'ok' if ok_rels_input else 'FAIL'} ({time.time()-t0:.1f}s)")
	la = M.linear_algebra(e['log2_U'], op['excess'])
	say(f"[5] 2^{e['log2_population']:.3f} cells at {e['per_cell']:.1f} instr/cell (routing {M.REC_ROUTE:.0f}/record, "
		f"{M.CELL_ROUTE:.0f}/cell; apply {M.C_APPLY:.2f}, recovery {M.C_RECOVER:.3f}): collection 2^{e['log2_relcol']:.4f} units; "
		f"N/U {la['ratio']:.4f} x pad {M.PAD_N} x w200 {M.W200_RATIO} -> N = 2^{la['log2_N']:.4f}; LA 2^{la['log2_la']:.4f} units ({M.C_ENTRY:.0f}/entry)")
	t0 = time.time()
	Us, res = filtsim.run(op['lb'], e['nf_mean'] / op['lb'], e['ng_mean'] / op['lb'], op['excess'],
						 seed, 12, kmax=64, targets=(100,), verbose=False)
	r100 = res[100]
	say(f"[6] filtering (2^12, seed {seed}, {time.time()-t0:.1f}s): N/U = {1/r100[0]:.4f} at weight {r100[3]:.0f} (max |coef| {r100[2]}); "
		"embedded 2^17 runs (excess: seed, U/N, N, max coef): " + "; ".join(f"{ex}: " + ", ".join(f"({s}, {un:.2f}, {n}, {mc})" for (s, un, n, mc, w) in v)
								 for ex, v in sorted(M.FILT_TABLE.items())))
	t0 = time.time()
	rec = filtsim.missing_recon(seed, excess=op['excess'])
	ok_rec = all(mf <= md for (lo, hi, n, mf, md, r1, ot, mg) in rec)
	say(f"[6b] missing logs after purge+clique+merge+backsub (2^12, seed {seed}, {time.time()-t0:.1f}s) per band "
		"(n, measured, model, r1/other/merge): "
		+ "; ".join(f"[{lo:g},{hi:g}): {n}, {mf:.4f}, {md:.4f}, {r1}/{ot}/{mg}" for (lo, hi, n, mf, md, r1, ot, mg) in rec) + f" {'ok' if ok_rec else 'FAIL'}")
	ok_zero = all(math.exp(-l) * (1 + l + l * l / 2) >= math.exp(-l) for l in [x / 10 for x in range(1, 400)])
	say(f"[6b-zero] proxy coef-1 >= e^-lambda on grid 0.1..39.9 {'ok' if ok_zero else 'FAIL'}")
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
	log2_EN, over, coll, pmin = repeat_bound(dres); OVER = 2.0 ** over; REP = 2.0 ** coll
	t0 = time.time()
	Rrels = 2.0 ** Pm.LOG2_RELS
	maxbits = max(r['nf'] + r['ng'] for r in pool)
	WMAX = 1500
	ok_wmax = maxbits <= 1400
	wmax = 400.0 / Rrels
	sum_c2 = Rrels * (2 * WMAX * wmax) ** 2
	T_ABS = 0.01
	logdelta_mcd = -2 * T_ABS * T_ABS / sum_c2
	ok_mcd = logdelta_mcd + math.log(63) <= math.log(1e-9) # per-node McDiarmid + union over 63
	union_hi = 0.0
	for b10 in range(80, 690): # v22r1: full range (b<40: e^-lambda underflows to 0.0)
		b = b10 / 10.0
		for um in (Pm.MAIN_UF, Pm.MAIN_UG):
			lam = 2.0 ** (Pm.LOG2_RELS - b) * 2 * DS.kappa(um, b / 69.0)
			if lam > 100:
				union_hi += 2.0 ** (b + 6) * math.exp(-lam) * (1 + lam + lam * lam / 2)
	ok_union = union_hi <= 1e-9
	u4b = math.ceil(4 * bt['m'] / 69.0)
	dom_ok = True
	krat_ok = True
	for k, d in dres['DET'].items():
		for m, um in ((d['nf'], Pm.MAIN_UF), (d['ng'], Pm.MAIN_UG)):
			u4 = math.ceil(4 * m / 69.0)
			for b10 in range(80, 690):
				if DS.kappa(u4 / 4.0, b10 / 690.0) > DS.kappa(u4b / 4.0, b10 / 690.0):
					dom_ok = False
				if DS.kappa(u4 / 4.0, b10 / 690.0) > 4 * DS.kappa(um, b10 / 690.0):
					krat_ok = False
	def cm_rate(lam): # v22r1: top band covers lam >= 8 (old code fell through past 1e9)
		for lo, hi, r in M.CM_RATES:
			if lo <= lam < hi:
				return r
		return M.CM_RATES[0][2] if lam < 1 else M.CM_RATES[-1][2]
	E_boot = 0.0
	for b10 in range(80, 690):
		b = b10 / 10.0
		lnw = math.log((b + 0.1) / b)
		for um in (Pm.MAIN_UF, Pm.MAIN_UG):
			lam = 2.0 ** (Pm.LOG2_RELS - b) * 2 * DS.kappa(um, b / 69.0)
			p_ge3 = 1 - math.exp(-lam) * (1 + lam + lam * lam / 2)
			E_boot += DS.kappa(u4b / 4.0, b / 69.0) * lnw * p_ge3 * cm_rate(lam)
	ok_cover = all(lost / occ <= M.CM_RATES[i][2] for scale in M.CM_BIG for i, (occ, lost) in enumerate(scale))
	tail0 = 0.0
	u4s = {u4b} | {math.ceil(4 * m / 69.0) for d in dres['DET'].values() for m in (d['nf'], d['ng'])}
	for LOG2R in (Pm.LOG2_RELS, 63.5): # shipped input + floor (tail decreases in lambda)
		for u4 in u4s:
			for b10 in range(80, 400):
				b = b10 / 10.0
				for um in (Pm.MAIN_UF, Pm.MAIN_UG):
					lam = 2.0 ** (LOG2R - b) * 2 * DS.kappa(um, b / 69.0)
					tail0 += DS.kappa(u4 / 4.0, b / 69.0) * math.log((b + 0.1) / b) * math.exp(-lam) * (1 + lam + lam * lam / 2)
	ok_tail0 = tail0 == 0.0
	T_M = Pm.T_RES - T_ABS
	delta_M = E_boot / T_M
	ok6c = ok_wmax and ok_mcd and ok_union and dom_ok and krat_ok and ok_cover and ok_tail0 and T_M > 0 and delta_M + 1e-9 + union_hi <= M.DELTA_TABLE
	say(f"[6c] table-good event ({time.time()-t0:.1f}s): pool norm bits <= {maxbits:.0f} -> WMAX {WMAX}; McDiarmid+union log-delta {logdelta_mcd:.1f}+log63 (<= ln 1e-9 {'ok' if ok_mcd else 'FAIL'}); high-lambda union {union_hi:.2e} {'ok' if ok_union else 'FAIL'}; boot-weight dominance {'ok' if dom_ok else 'FAIL'}; kappa ratio <= 4 {'ok' if krat_ok else 'FAIL'}; cover {'ok' if ok_cover else 'FAIL'}; tail0 {tail0!r} {'ok' if ok_tail0 else 'FAIL'}; E[boot residual] {E_boot:.4f}, T_M {T_M:.2f}, Markov delta {delta_M:.4f} + 1e-9 + union <= DELTA {M.DELTA_TABLE} {'ok' if ok6c else 'FAIL'}")
	p_desc = 1 - 2 * bt['phi'] - markov - OVER - REP - M.DELTA_TABLE
	p_main = 1.0 - 2.0 ** (la['log2_N'] + 2 - 254) # H3: main phase fails only by Wiedemann projections
	success = p_main * p_desc
	rows = DS.report(dres)
	sv = [dres['DET'][k] for k in dres['DET'] if dres['DET'][k]['how'] == 'sieve']
	hmax = max(DS.CELLS[d['R']] for d in sv)
	ok_width = hmax < 1 << 103 and (hmax - 1) >> 40 < 1 << 63 and all(16 * max(d['nf'], d['ng']) < 1 << 16 for d in sv)
	mu_min = min(r['mu'] for r in rows); psi_max = max(max(r['psi'] for r in rows), bt['psi'])
	say(f"[7] descent ({time.time()-t0:.1f}s): ball counts {'ok' if ok_counts else 'FAIL'}; {len(rows)} nodes, min mu {math.floor(mu_min*10)/10:.1f}, max psi {math.ceil(psi_max*1000)/1000:.3f}; "
		f"boot n1={bt['n1']} nA={bt['nA']:.2f} kA={bt['kA']} R={bt['R']}: 2^{bt['trials_log2']:.2f} trials at 2^{bt['trial_log2']:.2f} (pr 2^{bt['pr_log2']:.1f}), "
		f"own 2^{bt['own_log2']:.2f} + subtree 2^{bt['sub_log2']:.2f} x retry {bt['retry']:.3f}; one log 2^{bt['one_log_log2']:.3f}; two-log work 2^{two_log:.3f} units; "
		f"cap 2^{cap:g}: Markov <= 2^{math.ceil((two_log-cap)*1000)/1000:.3f}; P(boot fails) < {math.ceil(bt['phi']*1e7)/1e7:.2e}; E[N] <= 2^{math.ceil(log2_EN * 100) / 100:.2f} (counts16), N_MAX 2^{M.NODE_CAP_LOG2:g}, P(abort) <= 2^{over:.2f}; repeats in 2^{pmin:.2f}-ideal bin: P <= 2^{coll:.2f}, charged; success = {math.floor(success*10000)/10000:.4f} (floor); widths: half count <= "
		f"2^{math.log2(hmax):.2f}, owner < 2^{math.log2((hmax - 1) >> 40):.2f}, max score {16 * max(max(d['nf'], d['ng']) for d in sv):.0f} {'ok' if ok_width else 'FAIL'}")
	t0 = time.time(); sens = []
	for k_se in (0, 1, 2, 3):
		Ps = DS.P(); Ps.LOG2_RELS = math.floor((e['log2_rels'] - k_se * e['yield_se_log2']) * 1000) / 1000; DS.miss_side.cache_clear()
		sens.append((k_se, Ps.LOG2_RELS, DS.evaluate(Ps, DS.PLAN16, DS.BOOT16)['two_log_units']))
	DS.miss_side.cache_clear()
	ok_rels_2se = DS.P.LOG2_RELS <= e['log2_rels'] - 2 * e['yield_se_log2'] and abs(sens[2][2] - two_log) < 1e-9
	say(f"[7s] input sensitivity ({time.time()-t0:.1f}s), two-log at LOG2_RELS = est - k SE (SE {e['yield_se_log2']:.4f}): " + "; ".join(f"k={k}: {v:.3f} -> 2^{w:.4f}" for k, v, w in sens) + f"; shipped {DS.P.LOG2_RELS} <= {e['log2_rels'] - 2 * e['yield_se_log2']:.4f} {'ok' if ok_rels_2se else 'FAIL'}")
	nb, np_ = nbar(8, 12)
	say(f"[9] deg-1 ideals/prime, both sides, {np_} primes in [2^8, 2^12): {nb:.3f} (model {M.NBAR})")
	ok9 = nb <= M.NBAR * 1.01
	la_n, la_br = check_la_kernel()
	tail, red9, tail_n, tail_ok = la_row_tail()
	say(f"[10] row tail (la_mac.s): {tail_n} instr (loop-free, calls {tail_ok}); {tail}/row <= C_ROW {M.C_ROW:.0f}; 9-limb reduction {red9} <= C_RED9 {M.C_RED9:.0f}")
	nu_ok = max(M.W200.values()) <= M.W200_RATIO and all(max(M.SCALE_NU[a]) >= max(M.SCALE_NU[b]) for a, b in zip(sorted(M.SCALE_NU), sorted(M.SCALE_NU)[1:])) and M.ratio_N_over_U(M.OP['excess']) * M.PAD_N >= max(max(v) for k, v in M.SCALE_NU.items() if k >= 15)
	say(f"[10] N/U by scale (max over seeds): " + ", ".join(f"2^{k} {max(v):.4f}" for k, v in sorted(M.SCALE_NU.items())) + f"; non-increasing, charged {M.ratio_N_over_U(M.OP['excess']) * M.PAD_N:.5f} >= every 2^15+ value; w200 max {max(M.W200.values()):.3f} <= {M.W200_RATIO}: {nu_ok}")
	say(f"[10] la_mac.s group loop: {la_n} instr, branches {la_br} {'ok' if la_n == 33 and la_br == ['b.ne'] else 'MISMATCH'}")
	bw_n, bw_br = check_la_kernel("bw_dot.s"); bf_n, bf_br = check_la_kernel("bw_dot.s", "LBB1_2:")
	say(f"[10] bw_dot.s bw_full inner loop: {bf_n} instr (charged {M.C_FULL:.0f}) {'ok' if bf_n == 97 and bf_n <= M.C_FULL and bf_br == ['b.ne'] else 'MISMATCH'}")
	t0 = time.time()
	au = audit(bpr=bt['pr_log2'])
	ok11 = au['boot']['within_model_bins'] and au['specialq_within_bins']
	say(f"[11] audit ({time.time()-t0:.1f}s): boot lattice ({au['boot']['swaps']} swaps), boot cg max {max(s['cg_max'] for s in au['boot']['stats'])} "
		f"<= {max(x[1] for x in DS.P.BOOT_BINS):.1f} {'ok' if au['boot']['within_model_bins'] else 'FAIL'}; special-q lattices "
		+ "; ".join(f"{s['side']}@2^{s['bits']}: cf max {max(x['cf_max'] for x in s['stats'])}, cg max {max(x['cg_max'] for x in s['stats'])}" for s in au['specialq'])
		+ f" within bins (cf <= {max(b[1] for b in DS.P.QBINS):.1f}, cg <= {max(b[2] for b in DS.P.QBINS):.1f}) {'ok' if au['specialq_within_bins'] else 'FAIL'}; "
		f"relw E[shat/m]/pr (live): boot {au['relw']['boot_yield_ratio']:.3f} (> 0.70, 96 cells, {'ok' if au['relw']['boot_yield_ratio'] > 0.70 else 'FAIL'}) (UDUP {DS.P.UDUP}); "
		f"sampler shell-chi2 {au['sampler']['chi2']:.1f} (< 60 {'ok' if au['sampler']['ok'] else 'FAIL'}); "
		f"largest line <= 2^{au['storage']['line_cells_log2_max']} cells; chunk protocol {au['storage']['routing_stable']}")
	ok11r = au['relw']['boot_yield_ratio'] > 0.70
	ok11s = au['sampler']['ok']
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
	say(f"[11b] retention ({time.time()-t0:.1f}s): prime-ideal table recomputed at {len(spot)} primes {'ok' if ok_tab else 'FAIL'}; "
		"repeated roots or ramification at " + "; ".join(f"{p}: {tab[str(p)]}" for p in special)
		+ "; residue samples (400 classes, sample/bound): " + ", ".join(f"{k} {v[0]}/{v[1]}" for k, v in sample.items()) + f" {'ok' if ok_sample else 'FAIL'}"
		+ f"; base-primitive {RT.base_primitive():.4f}, ordinary {ordf:.4f}/{ordg:.4f}; omitted mass means {ex(pf):.2f}/{ex(pg):.2f} bits, "
		f"P(X > 96) {RT.cell_loss('f', 96):.2e}/{RT.cell_loss('g', 96):.2e}; loss <= {loss_max:.1e}, retention >= {ret_min:.4f}; bounds >= |N|: " + "; ".join(f"{b}{s_} R={R_}: {ok_}, max gap {g:.1e} bits" for b, s_, R_, ok_, g, k in rr)
		+ "; empty-window toys: " + "; ".join(str(tuple(x)) for x in nb['rows'])
		+ f" {'ok' if nb_ok else 'FAIL'}")
	uo = st['unit_orbit']
	ok_orb = math.floor(uo['boot_exact_ci_lower_charged'] * 100) / 100 >= uo['boot_claim'] >= DS.P.UDUP
	say(f"[11u] unit orbit ({uo['method']}): boot {uo['boot_points']} cells yr {uo['boot_exact_yr_charged']} [{uo['boot_exact_ci_lower_charged']}, {uo['boot_exact_ci_upper_charged']}] vs charged pr (claim {uo['boot_claim']}, UDUP {DS.P.UDUP}; nodes: [11s]) {'ok' if ok_orb else 'FAIL'}")
	t0 = time.time(); ss = sstar(dres, st); ok_ss = ss['same'] and ss['worst'][0] >= 0 and ss['err'] < 1e-6
	say(f"[11s] S* ({time.time()-t0:.1f}s): {ss['n']} classes, live == campaign set {ss['same']}; {ss['lat']} lattices x 200 cells, {ss['k']} flagged; tightest {ss['worst'][3]}: P(G < {ss['worst'][2]:.2f}) <= CP {ss['worst'][1]:.4f}/(1-beta), slack {ss['worst'][0]:+.4f}; replay {ss['nrep']} cells, rel err {ss['err']:.1e} {'ok' if ok_ss else 'FAIL'}")
	t0 = time.time()
	toyecm = EC.toy(60, seed, sizes=((36, 1024, 50), (40, 2048, 50)))
	sch = [EC.schedule(69.0, int(math.ceil(dres['DET'][((69, 75), 'f')]['nf'])) ), EC.schedule(bt['nA'], int(math.ceil(bt['m'])))]
	ok12 = all(t['empirical'] >= t['model_corrected'] - 2 * math.sqrt(t['model_corrected'] * (1 - t['model_corrected']) / t['trials']) for t in toyecm)
	say(f"[12] ECM schedules ({time.time()-t0:.1f}s): " + "; ".join(
		f"{x['fbits']:.2f}-bit primes in {x['nbits']}-bit integers: B1 2^{math.log2(x['B1']):.0f}, B2 {x['B2'] // x['B1']} B1, wheel {x['wheel']}, "
		f"{x['curves']} curves at p = {x['p_curve']:.2e}, {x['stage1_M']:.3g} + {x['stage2_M']:.3g} M per curve, 2^{x['log2_instructions']:.2f} instructions" for x in sch)
		+ "; toy detection: "
		+ "; ".join(f"({t['bits']}, {t['trials']}, {t['stage1']}, {t['stage2']}, {t['empirical']:.3f}, {t['model']:.3f}, {t['model_corrected']:.3f})" for t in toyecm)
		+ f" {'ok' if ok12 else 'FAIL'}")
	t0 = time.time(); fl = check_fields()
	say(f"[13] number fields ({time.time()-t0:.1f}s): N(f),N(g) of degree {fl['deg_f']},{fl['deg_g']} squarefree mod r {'ok' if fl['unramified'] else 'FAIL'} (r unramified); "
		f"disc h = {fl['disc_h']}; odd unit rank {fl['sig']['odd']} = SM {M.N_SM} {'ok' if fl['sig']['odd'] == M.N_SM else 'FAIL'}")
	ok13 = fl['unramified'] and fl['sig']['odd'] == M.N_SM
	import time as _t
	_t0 = _t.time(); al_s = {p: alpha_contrib(p) for p in (2, 19, 23)}
	ok_alpha = op['alpha_f'] >= ALPHA_F and op['alpha_f'] >= ALPHA_F_50K and op['alpha_g'] >= ALPHA_G and op['alpha_g'] >= ALPHA_G_50K and DS.P.ALPHA_G >= ALPHA_G and DS.P.ALPHA_G >= ALPHA_G_50K and DS.P.ALPHA_F >= ALPHA_F and DS.P.ALPHA_F >= ALPHA_F_50K and all(round(v[0], 3) >= 0 for v in al_s.values())
	say(f"[14] alpha (p <= {ALPHA_PMAX}, p < 50000, exact roots): computed ({ALPHA_F:.6f}, {ALPHA_G:.6f}; 50k ({ALPHA_F_50K:.6f}, {ALPHA_G_50K:.6f})) <= charged ({op['alpha_f']}, {op['alpha_g']}) {'ok' if ok_alpha else 'FAIL'}; "
		f"exact checks at p=2/19/23: " + ", ".join(f"{p}: ({v[0]:.6f}, {v[1]:.6f})" for p, v in sorted(al_s.items())) + f" ({_t.time()-_t0:.1f}s)")
	tot = M.total_log2(e['log2_relcol'], la['log2_la'], cap, la['log2_gen'])
	say(f"ledger: collection {e['log2_relcol']:.4f}, LA {la['log2_la']:.4f}, generator {la['log2_gen']:.2f}, "
		f"descent cap {cap:g} (exp {two_log:.3f}), canon list {M.LOG2_CANON:.0f}, Schirokauer {M.LOG2_SM:.0f}, "
		f"factor base {M.LOG2_FBC:.0f}, misc {M.LOG2_MISC:.0f}, verify {M.LOG2_VERIFY:.0f}")
	say(f"total = {tot:.6f} (claimed {CLAIM:.6f}: {'ok' if tot <= CLAIM else 'EXCEEDS CLAIM'}); "
		f"reserve {CLAIM - tot:.3f} bits")
	vals = dict(ledger=tot, relcol=e['log2_relcol'], la=la['log2_la'], log2_rels=e['log2_rels'], margin=e['margin'],
				work_se=work_se, two_log=two_log, markov=markov, log2_EN=log2_EN, over=math.log2(OVER), rep=math.log2(REP),
				boot_phi=bt['phi'], success=success)
	dec_ok = {k: (vals[k] <= v if op == '<=' else vals[k] >= v) for k, (op, v) in DECLARED.items()}
	dec_ok['work_factor'] = 786 ** 12 * 10 ** 9 <= 784 ** 12 * 1031045427 and float(corr['safe_work_population_factor']) <= 1.031045427 # (786/784)^12 ceiling at a-radius 784
	ok_declared = all(dec_ok.values())
	say("declared bounds: " + "; ".join(f"{k} {vals[k]!r} {op} {v} {'ok' if dec_ok[k] else 'FAIL'}" for k, (op, v) in DECLARED.items()))
	checks = dict(owner_schedule=O.check(M, la), resultant=ok0, tower=ok0b, rho_semismooth=ok1, pool=ok2, nbar=ok9, la_loop=la_n == 33 and bf_n == 97 and bf_br == ['b.ne'] and check_la_kernel('la_mac.s', 'LBB9_1:') == (11, ['b.ne']) and tail_ok and tail <= M.C_ROW and red9 <= M.C_RED9 and nu_ok and M.C_ENTRY >= 4 + 3 + M.UNIT_FRAC_CHARGED * 11 + (1 - M.UNIT_FRAC_CHARGED) * 33,
				 ball_counts=ok_counts and ok_relmeans, rels_input=ok_rels_input and ok_rels_2se, work_allowance=ok_allow, declared_bounds=ok_declared, alpha_ok=ok_alpha, missing_logs=ok_rec and ok_zero, table_concentration=ok6c, record_widths=ok_width, descent_plan_feasible=mu_min >= Pm.MIN_RELS and psi_max <= Pm.PSI_MAX and max(r['phi'] for r in rows) <= 0.05 and coll <= -4,
				 audit=ok11, retention=ok11b, unit_orbit=ok_orb, sstar=ok_ss, ecm=ok12, fields=ok13, feasible=e['feasible'], within_claim=tot <= CLAIM, success_ok=success >= 0.39, udup_relw=ok11r, sampler_uniform=ok11s)
	say(f"done in {time.time()-T0:.1f}s; checks: {checks}")
	print(json.dumps({"checks_ok": all(checks.values()), "checks": checks,
					 "total_log2W_full_ledger": math.ceil(tot * 1000000) / 1000000, "claimed": CLAIM,
					 "log2_relcol_ceil4": math.ceil(e['log2_relcol'] * 10000) / 10000, "log2_rels_floor4": math.floor(e['log2_rels'] * 10000) / 10000,
					 "margin_bits_floor4": math.floor(e['margin'] * 10000) / 10000, "work_se_ceil5": math.ceil(work_se * 100000) / 100000, "work_allowance": M.CANON_WORK_ALLOWANCE, "yield_se_log2": round(e['yield_se_log2'], 4), "log2_N": round(la['log2_N'], 4),
					 "log2_la_ceil4": math.ceil(la['log2_la'] * 10000) / 10000, "descent_two_logs_log2_units_ceil4": math.ceil(two_log * 10000) / 10000,
					 "descent_cap": cap, "markov_abort": markov, "log2_expected_nodes_ceil4": math.ceil(log2_EN * 10000) / 10000, "node_cap_log2": M.NODE_CAP_LOG2, "counter_abort": OVER, "repeat": REP, "success_exact": success, "modelled_success_floor4": math.floor(success * 10000) / 10000,
					 "descent_nodes": len(rows), "descent_boot_log2": round(bt['one_log_log2'], 3)}), flush=True)
if __name__ == "__main__":
	if "--alpha" in sys.argv:
		print(json.dumps(dict(zip(("alpha_f", "alpha_g"), alpha_sweep(ALPHA_PMAX))), sort_keys=True)); sys.exit(0)
	if "--alpha50k" in sys.argv:
		print(json.dumps(dict(zip(("alpha_f", "alpha_g"), alpha_sweep(ALPHA_G_PMAX))), sort_keys=True)); sys.exit(0)
	main()
