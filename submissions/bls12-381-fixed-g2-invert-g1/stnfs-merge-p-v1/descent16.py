"""stnfs descent model (DESCENT.md): nodes, lower sums, survivors, the work recurrence."""
import math, functools
from decimal import Decimal, localcontext
LN2 = math.log(2); UNIT = 20000.0
KMAX = 64
RSTEP = 1.0 / 1024
def _rho_coefficients(kmax=KMAX, nt=600):
	with localcontext() as ctx:
		ctx.prec = 200
		zero, one = Decimal(0), Decimal(1)
		coefs = [None, [one] + [zero] * (nt - 1)]
		for k in range(2, kmax + 1):
			previous, c = coefs[k - 1], [zero] * nt
			for m in range(nt - 1):
				c[m + 1] = (previous[m] + m * c[m]) / (k * (m + 1))
			c[0] = previous[0] - sum(c[1:])
			coefs.append(c)
		out = [None]
		for c in coefs[1:]:
			big = max(abs(x) for x in c)
			while len(c) > 1 and abs(c[-1]) < big * Decimal('1e-25'):
				c = c[:-1]
			out.append([float(x) for x in c])
		return out
RHO_COEFS = _rho_coefficients()
def rho_series(u):
	if u <= 1:
		return 1.0
	k = math.ceil(u)
	if k > KMAX:
		return 0.0
	x, value = k - u, 0.0
	for c in reversed(RHO_COEFS[k]):
		value = value * x + c
	return max(value, 0.0)
RHO_TAB = [rho_series(i * RSTEP) for i in range(int(KMAX / RSTEP) + 2)]
def rho_lo(u):
	"""lower bracket of rho(u): the table value at the next grid point."""
	if u <= 1:
		return 1.0
	i = math.ceil(u / RSTEP)
	return RHO_TAB[i] if i < len(RHO_TAB) else 0.0
def rho_hi(u):
	"""upper bracket: the previous grid point;"""
	if u <= 1:
		return 1.0
	i = math.floor(u / RSTEP)
	if i + 1 >= len(RHO_TAB):
		return 1.0 / math.factorial(math.floor(u))
	return RHO_TAB[i]
def kappa(u, t):
	return rho_series(u - t) / rho_series(u) if rho_series(u) > 0 else 0.0
CALIB = ((1.0, 2.0, 0.172), (2.0, 3.0, 0.117), (3.0, 5.0, 0.057), (5.0, 8.0, 0.02))   # (lambda lo, hi, c_r): charged unknown fraction of OCCURRING ideals after purge and back-substitution ([7c])
CALIB_MAX = {(1.0, 2.0): 0.1141, (2.0, 3.0): 0.0775, (3.0, 5.0): 0.0362, (5.0, 8.0): 0.0071}   # H5c-mean premise means: the largest of 16 filtsim.missing_recon runs at excess 1.1 (DESCENT.md)
ZERO_EPS = 2.0 ** -8   # relative allowance on the never-occurring mass S_0 (two-sided Chebyshev, table_tails)
TAIL_A, TAIL_K = 0.5, 0.45   # lambda >= 8: occurring-ideal loss charged at max(P(N <= 2 | N >= 1), TAIL_A e^(-TAIL_K lambda)) (DESCENT.md)
COVER_K = 1.5          # H5c-mean: charged occurring fraction >= COVER_K x the premise mean
XI_PREM = 2.0 ** 14    # H5c-ES: N_s n E[Phi_s^2] / Lambda_s^2 <= XI_PREM (DESCENT.md)
DELTA_CONC = 2.0 ** -30   # charged against success: >= the union of the S_0 and S_1 table bounds (table_tails)
U4_RANGE = range(24, 161)   # relation sizes u = u4/4 of the table event ([7] checks the plan's)
def c_occ(lam):
	"""charged fraction of occurring ideals left without a logarithm (non-increasing)."""
	if lam < 1.0:
		return 1.0
	for lo, hi, c in CALIB:
		if lo <= lam < hi:
			return c
	p0 = math.exp(-lam)
	return max(p0 * (lam + lam * lam / 2) / (1 - p0), TAIL_A * math.exp(-TAIL_K * lam))
def p_band(lam_lo, lam_hi):
	"""charged P(no logarithm), lambda in [lam_lo, lam_hi]: never occurring plus occurring but unknown, each at its worst lambda."""
	return min(1.0, (1.0 + ZERO_EPS) * math.exp(-lam_lo) + (1.0 - math.exp(-lam_hi)) * c_occ(lam_lo))
@functools.lru_cache(None)
def kap_edges(u):
	"""kappa(u, t_j), t_j = (40 + j/10)/69."""
	r = rho_series(u)
	return tuple(rho_series(u - (400 + j) / 690.0) / r for j in range(291))
SUB = 10   # sub-bands per 1/10-bit band: main-phase lambda and 1/q are bounded on 1/100-bit sub-bands
@functools.lru_cache(None)
def band_table(log2_rels, u_main):
	"""one side's ideal universe (DESCENT.md): per 1/10-bit band, sums over 1/100-bit sub-bands k (n_k sigma-orbits, lambda in
	[lo_k, hi_k]): (b, n, lam_lo, lam_hi, A, A0, B0, V0, C1, Q2), A = sum n_k 2^-b_k P_miss_k, A0/B0 = upper/lower never-occurring
	sums, V0 = sum n_k 4^-b_k e^-lo_k, C1 = sum n_k 2^-b_k (1 - e^-hi_k) c_occ(lo_k), Q2 = sum n_k P_miss_k 4^-b_k."""
	r = rho_series(u_main); kk = [rho_series(u_main - (4000 + i) / 6900.0) / r for i in range(2901)]
	out = []
	for j in range(290):
		acc = [0.0] * 7; n_j = 0.0; lmin = 1e300; lmax = 0.0
		for i in range(SUB * j, SUB * j + SUB):
			b = (4000 + i) / 100.0; b1 = (4001 + i) / 100.0
			n = (2 ** log2li(b1) - 2 ** log2li(b)) / 2; q = 2.0 ** -b
			lo = 2.0 ** (log2_rels - b1) * 2 * kk[i]; hi = 2.0 ** (log2_rels - b) * 2 * kk[i + 1]
			pm = p_band(lo, hi); elo = math.exp(-lo); ehi = math.exp(-hi)
			for t, v in enumerate((n * q * pm, n * q * elo, n * q * 2.0 ** -0.01 * ehi, n * q * q * elo, n * q * (1 - ehi) * c_occ(lo), n * pm * q * q)):
				acc[t] += v
			n_j += n; lmin = min(lmin, lo); lmax = max(lmax, hi)
		out.append(((400 + j) / 10.0, n_j, lmin, lmax) + tuple(acc[:6]))
	return tuple(out)
@functools.lru_cache(None)
def miss_side(u4, log2_rels, u_main):
	"""charged missing-log children of one side of size u4/4."""
	k = kap_edges(u4 / 4.0)
	return sum(2 * k[j + 1] * row[4] for j, row in enumerate(band_table(log2_rels, u_main)))
def log2_n_side():
	"""log2 of one side's orbits with t >= 0.8 (the toy's range)."""
	return math.log2((2 ** log2li(69.0) - 2 ** log2li(55.2)) / 2)
@functools.lru_cache(None)
def table_tails(log2_rels, u_main_f, u_main_g, xi=XI_PREM):
	"""the realized-table event (DESCENT.md): log2 of (a) the S_0 two-sided Chebyshev union, (b) the S_1 Efron-Stein union,
	4 (G/g)^2 xi/N_s per (side, u4), (c) sum_I m_I wh_I^2 over both sides, wh_I = max_u w_I^+(u)/(B_0 - ZERO_EPS A_0)."""
	d0 = d1 = pair = 0.0; ns = 2.0 ** log2_n_side()
	for um in (u_main_f, u_main_g):
		tab = band_table(log2_rels, um); km = kap_edges(um); wh = [0.0] * len(tab)
		for u4 in U4_RANGE:
			k = kap_edges(u4 / 4.0); A0 = B0 = V0 = 0.0
			for j, row in enumerate(tab):
				kp = 2 * k[j + 1]; A0 += kp * row[5]; B0 += 2 * k[j] * row[6]; V0 += kp * kp * row[7]
			g = [k[j] / km[j] for j in range(291)]
			d0 += V0 / (ZERO_EPS * A0) ** 2
			d1 += (max(g) / min(g)) ** 2 * xi / ns / (COVER_K - 1) ** 2
			zlo = B0 - ZERO_EPS * A0
			for j in range(len(tab)):
				w = 2 * k[j + 1] / zlo
				if w > wh[j]: wh[j] = w
		pair += sum(w * w * row[9] for row, w in zip(tab, wh))
	return math.log2(d0), math.log2(d1), math.log2(pair)
def clear_caches():
	for f in (band_table, miss_side, table_tails):
		f.cache_clear()
def universe_check(Pm, u4=60):
	"""[7u] miss_side against a direct 1/100-bit sum; log2 sum n_k against log2 (Li(2^69) - Li(2^40))/2."""
	out = []
	for um in (Pm.MAIN_UF, Pm.MAIN_UG):
		direct = 0.0; n_tot = 0.0
		for i in range(2900):
			b = (4000 + i) / 100.0; b1 = b + 0.01; n = (2 ** log2li(b1) - 2 ** log2li(b)) / 2; n_tot += n
			lo = 2.0 ** (Pm.LOG2_RELS - b1) * 2 * kappa(um, b / 69.0); hi = 2.0 ** (Pm.LOG2_RELS - b) * 2 * kappa(um, b1 / 69.0)
			direct += n * 2 * kappa(u4 / 4.0, (400 + i // SUB + 1) / 690.0) * 2.0 ** -b * p_band(lo, hi)
		out.append((miss_side(u4, Pm.LOG2_RELS, um), direct, math.log2(n_tot), math.log2((2 ** log2li(69.0) - 2 ** log2li(40.0)) / 2)))
	return out
def missing(Pm, nf, ng):
	"""missing-log ideals per relation: f side at nf bits, g side at ng bits."""
	f = lambda m, um: miss_side(math.ceil(4 * m / 69.0), Pm.LOG2_RELS, um)
	return (f(nf, Pm.MAIN_UF) if nf else 0.0) + (f(ng, Pm.MAIN_UG) if ng else 0.0)
KA_MAX = 16
def _conv(a, b):
	"""out[k] = sum_i a[i] b[k - i] over the nonzero supports, one dot product per k."""
	from operator import mul as _mul
	la, lb = len(a), len(b)
	out = [0.0] * (la + lb - 1)
	A = la
	while A > 0 and a[A - 1] == 0.0:
		A -= 1
	B0 = 0
	while B0 < lb and b[B0] == 0.0:
		B0 += 1
	B1 = lb
	while B1 > B0 and b[B1 - 1] == 0.0:
		B1 -= 1
	if A == 0 or B0 >= B1:
		return out
	bb = b[B0:B1][::-1]; L = B1 - B0
	for k in range(B0, A + B1 - 1):
		kk = k - B0
		i0 = max(0, kk - L + 1); i1 = min(kk, A - 1)
		out[k] = sum(map(_mul, a[i0:i1 + 1], bb[L - 1 - kk + i0:L - kk + i1]))
	return out
@functools.lru_cache(None)
def window_rows(n, base, iA, mesh, kmax=KA_MAX, kB=1):
	"""rows[(jA, jB)]: measure/(jA! jB!) of jA A-type and jB B-type interval products by index sum; h = (n - base)/mesh; jB rows only when the B window is used."""
	h = (n - base) / mesh
	wA = [math.log1p(h / (base + i * h)) if i < iA else 0.0 for i in range(mesh)]
	wB = [math.log1p(h / (base + i * h)) if i >= iA else 0.0 for i in range(mesh)]
	PA = [[1.0]]
	for j in range(1, kmax + 1):
		PA.append([x / j for x in _conv(PA[-1], wA)])
	PB = [[1.0], wB]
	rows = {}
	for jA in range(kmax + 1):
		for jB in range(2 if kB else 1):
			rows[(jA, jB)] = PA[jA] if jB == 0 else _conv(PA[jA], PB[1])
	return rows, h
@functools.lru_cache(None)
def ss_cum2(m, base, n, iA, mesh, alpha, trunc, kmax=KA_MAX, kB=1):
	"""lower sums L[jA][jB] (cumulative in jA) at integer m."""
	rows, h = window_rows(n, base, iA, mesh, kmax, kB)
	out = [[0.0, 0.0] for _ in range(kmax + 1)]
	for (jA, jB), row in rows.items():
		j = jA + jB
		imax = len(row) - 1
		ifeas = min(imax, int(math.floor((m - j * h - j * base) / h))) if h > 0 else imax
		if ifeas < 0:
			continue
		imin = 0 if not trunc else max(0, ifeas - int(trunc * base / h) - 1)
		tab, ntab, ceil = RHO_TAB, len(RHO_TAB), math.ceil
		k0 = (m - j * base + alpha) / (base * RSTEP); hk = h / (base * RSTEP)
		s = 0.0
		for i in range(imin, ifeas + 1):
			weight = row[i]
			if weight:
				k = ceil(k0 - i * hk)
				s += weight * (1.0 if k <= 0 else (tab[k] if k < ntab else 0.0))
		out[jA][jB] = s
	cum = [[0.0, 0.0] for _ in range(kmax + 1)]
	for jB in range(2):
		acc = 0.0
		for jA in range(kmax + 1):
			acc += out[jA][jB]
			cum[jA][jB] = acc
	return tuple(tuple(x) for x in cum)
def ss_lower2(m, base, n, iA, kA, kB, mesh, alpha, trunc=3.0):
	"""P(all primes <= 2^n, <= kA in (2^base, 2^nA], <= kB in (2^nA, 2^n]) lower sum."""
	m = math.ceil(m)
	if n <= base:
		return rho_lo((m + alpha) / base)
	cum = ss_cum2(m, base, n, iA, mesh, alpha, trunc, kA, 1 if kB else 0)
	return cum[kA][0] + (cum[kA][1] if kB else 0.0)
def surv_upper(m, base, T):
	"""P(product of the prime factors above 2^base of an m-bit integer <= 2^T) <= (T/base) rho((m - T)/base) for base <= T < m."""
	if T <= 0:
		return 0.0
	if T >= m:
		return 1.0
	return min(1.0, max(1.0, T / base) * rho_hi((math.floor(m) - T) / base))
MONT5 = 367.76               # measured five-limb Montgomery multiplication (ecm12.c, KERNELS.md)
def ecm_instr(fbits, nbits, success):
	"""instructions of the ecm16 schedule finding every prime of <= fbits bits in an nbits-bit integer with probability >= success per prime."""
	import ecm16
	return ecm16.ecm_instr(fbits, nbits, success)
def log2li(l):
	x = l * LN2
	return l - math.log2(x) + math.log2(1 + 1 / x + 2 / x / x)
CELLS = {4:13315520,6:1573359464,8:47982017060,10:687193400292,12:6073556677752,14:38414493101628,
		 16:190058478539564,20:2754504316057176,24:24505531366261200,28:155610653675916932,
		 32:771909602591773000,40:11221331971811657968,48:99995402971660463340,
		 56:635616373364584622852,64:3155035383949694198996,80:45900097950907999286464,
		 96:409193506988359583551012,112:2601677688965738307448548,128:12916156661164168796080328,160:187942845067274337246511324,192:1675659512501707896389090868,
		 224:10654623198171176305341831664,256:52897551248440914796084637156,288:217401149762189752879411886484,
		 320:769748406668707279529221545272,384:6863095027560050329435265332068}
class P:
	"""the parameter set (proof.md / DESCENT.md)."""
	LB = 69.0; LX = 69.0; PMIN = 8.0; SLACK = 96.0
	ALPHA_F = 5.15; ALPHA_G = 1.0; LOG2U6 = 6 * 63.714   # Murphy alpha charged above the p < 50,000 sums (5.085, 0.732; [14])
	QBINS = [(1 / 6, -26.4 + 23, 14.5 + 6), (1 / 6, -11.2 + 23, 14.5 + 6), (1 / 6, -0.1 + 23, 14.5 + 6),
			 (1 / 6, 10.5 + 23, 14.5 + 6), (1 / 6, 23.2 + 23, 14.5 + 6), (1 / 6, 86.6 + 23, 14.5 + 6)]
	BOOT_BINS = [(0.5, -191.5 + 4), (0.4, -185.9 + 4), (0.1, -175.0)]   # weighted boot basis, uniform points: p50/p90 + 4; top decile at the max of 21,504 norms (DESCENT.md, [11])
	QBINS_LO = [(-140.05 - 23, -54.56 - 6), (-26.4 - 23, -54.56 - 6), (-11.2 - 23, -54.56 - 6),
				(-0.1 - 23, -54.56 - 6), (10.5 - 23, -54.56 - 6), (23.2 - 23, -54.56 - 6)]
	LOG2_RELS = 63.68; MAIN_UF = 741.2 / 69; MAIN_UG = 498.7 / 69   # pool estimate - 2 SE (floor), means at or below the relation-weighted ones; [4] checks both
	THIN_MAX = 90               # nodes with q <= 2^95 discard relations touching a missing log
	UDUP = 0.25                 # distinct relations per raw relation (unit duplicates; stratified relation-conditioned minimum 0.283, DESCENT.md)
	MIN_RELS = 16.0; PSI_MAX = 0.8; PHI_FLOOR = 0.02
	TABLE_SLACK = 2.0 ** -16    # relative margin kept on the boot failure (DESCENT.md)
	RET_FLOOR = 0.8592          # retention floor of every sieve node (part=losses)
	COUNTED = 0.8592             # the boot's counted fraction (part=losses)
	NB_R1 = 0.5; NB_R2 = 1.0    # negative-binomial shapes: windows of <= 1 / >= 2 primes (toy fits 0.60 / 1.86-12.2)
	BOUND_CELL = 2000.0         # certified per-cell norm bounds, both sides (84 forms; ~1500 measured, KERNELS.md)
	HITS_PER_CELL = 2 * math.log(69.0 / 8.0)
	SIEVE_HIT = 850.0; SCAN = 10.0; ZERO = 2.0   # f and g score arrays (757 enumeration + 84 transfer + 9 apply; KERNELS.md)
	PAIR_SETUP = 2.0 ** 25.1; LLL_NODE = 2.0 ** 36; NORMS = 2.0 ** 26; TDIV = 2.0 ** 24.5; IDEAL = 2.0 ** 34; PRIME = 2.0 ** 24
	ECM_SUCCESS = 0.999
	RADII = (4, 6, 8, 10, 12, 14, 16, 20, 24, 28, 32, 40, 48, 56, 64, 80, 96, 112, 128, 160, 192, 224, 256, 288, 320, 384)
	MESH = 384; TRUNC = 0.0
	LOG2P6 = 6 * 380.99700
	LMAX = 225; BIN = 5
def pairs_count(Pm):
	return 2.0 * 2 ** log2li(Pm.LX)
def node_cost(Pm, l_lo, l_hi, side, n, iA, kA, kB, R, thin=None):
	"""own cost (instructions) and statistics of one node plan;"""
	pairs = pairs_count(Pm)
	cfmax = max(b[1] for b in Pm.QBINS); cgmax = max(b[2] for b in Pm.QBINS)
	S = math.log2(CELLS[R]); c = (l_hi + S) / 12
	h = (n - Pm.LX) / Pm.MESH if n > Pm.LX else 0.0
	nA = Pm.LX + iA * h if n > Pm.LX else Pm.LX
	T = kA * nA + kB * n + Pm.SLACK
	c_lo = (l_lo + S) / 12
	pr = 0.0; sv = 0.0
	for (w, cf, cg), (cf_lo, cg_lo) in zip(Pm.QBINS, Pm.QBINS_LO):
		nf = 72 * c + cf - (l_lo if side == 'f' else 0)
		ng = 12 * c + Pm.LOG2U6 + cg - (l_lo if side == 'g' else 0)
		pf = ss_lower2(nf, Pm.LX, n, iA, kA, kB, Pm.MESH, Pm.ALPHA_F, Pm.TRUNC)
		pg = ss_lower2(ng, Pm.LX, n, iA, kA, kB, Pm.MESH, Pm.ALPHA_G, Pm.TRUNC)
		pr += w * pf * pg
		nf_lo = 72 * c_lo + cf_lo - (l_hi if side == 'f' else 0)
		ng_lo = 12 * c_lo + Pm.LOG2U6 + cg_lo - (l_hi if side == 'g' else 0)
		sv += w * surv_upper(nf_lo, Pm.LX, T) * surv_upper(ng_lo, Pm.LX, T)
	if pr <= 0:
		return None
	nf_max = 72 * c + cfmax - (l_lo if side == 'f' else 0)
	ng_max = 12 * c + Pm.LOG2U6 + cgmax - (l_lo if side == 'g' else 0)
	det = 1.0 - (1 - Pm.ECM_SUCCESS) * ((nf_max + ng_max) / 20.0 + 2 * (kA + kB))
	miss = missing(Pm, nf_max, ng_max)
	if thin is None or l_lo <= Pm.THIN_MAX:
		thin = thin or l_lo <= Pm.THIN_MAX
	mu = CELLS[R] * pr * Pm.UDUP * det * (math.exp(-miss) if thin else 1.0)
	if mu < Pm.MIN_RELS:
		return None
	surv_cost = (Pm.NORMS + Pm.TDIV + ecm_instr(Pm.LX, nf_max, Pm.ECM_SUCCESS) + ecm_instr(Pm.LX, ng_max, Pm.ECM_SUCCESS)
				 + 2 * (ecm_instr(nA, T, Pm.ECM_SUCCESS) + Pm.PRIME) + Pm.IDEAL)
	survivors = CELLS[R] * sv + CELLS[R] * pr
	sieve = CELLS[R] * (Pm.HITS_PER_CELL * Pm.SIEVE_HIT + Pm.SCAN + Pm.ZERO + Pm.BOUND_CELL) + pairs * Pm.PAIR_SETUP
	trials = CELLS[R] * surv_cost
	if trials < sieve + survivors * surv_cost:
		own = Pm.LLL_NODE + trials; how = 'ecm'; survivors = CELLS[R]
		ret, loss_f, loss_g = 1.0, 0.0, 0.0
	else:
		own = Pm.LLL_NODE + sieve + survivors * surv_cost; how = 'sieve'
		ret, loss_f, loss_g = Pm.RET_FLOOR, None, None
		mu *= ret                 # retention floor (part=losses recomputes every node's exact value)
		if mu < Pm.MIN_RELS:
			return None
	return dict(own=own, n=n, iA=iA, nA=nA, kA=kA, kB=kB, R=R, S=S, mu=mu, how=how, miss=miss, ret=ret, loss_f=loss_f, loss_g=loss_g,
				miss_children=0.0 if thin else miss, thin=bool(thin), pr_log2=math.log2(pr),
				surv_log2=math.log2(survivors), own_log2=math.log2(own), nf=nf_max, ng=ng_max, det=det, T=T,
				sieve_log2=math.log2(sieve), survcost_log2=math.log2(survivors * surv_cost))
def bins_of(Pm):
	bins = [(int(Pm.LB), int(Pm.LB) + Pm.BIN + 1)]
	l = int(Pm.LB) + Pm.BIN + 1
	while l < Pm.LMAX:
		bins.append((l, l + Pm.BIN)); l += Pm.BIN
	return bins
def phi_nb(mean, Pm, window):
	"""P(no distinct relation) for a negative-binomial count of the given mean."""
	r = Pm.NB_R1 if window <= 1 else Pm.NB_R2
	return max((1.0 + mean / r) ** (-r), Pm.PHI_FLOOR)
def evaluate(Pm, PLAN, BOOT):
	"""evaluate the shipped plan: per-bin W (instructions), phi, details, and the two-log expected work;"""
	shared = ('LX', 'MESH', 'TRUNC', 'ALPHA_F', 'ALPHA_G')
	if any(getattr(Pm, name) != getattr(P, name) for name in shared):
		raise ValueError('descent and retention parameters must agree')
	bins = bins_of(Pm)
	W = {}; PHI = {}; DET = {}
	def wmax(side, n):
		return max([W[(b, side)] for b in bins if b[0] < n - 1e-9] + [0.0])
	def phimax(side, n):
		return max([PHI[(b, side)] for b in bins if b[0] < n - 1e-9] + [0.0])
	def miss_terms():
		"""failure probability and work of a missing-log child Q (prime special-q node in its own lattice L_Q): the larger base-bin values over both sides (DESCENT.md)."""
		base = [k for k in W if k[0][0] == bins[0][0]]
		return (max([PHI[k] for k in base] + [0.0]), max([W[k] for k in base] + [0.0]))
	for b in bins:
		for side in ('f', 'g'):
			n, iA, kA, kB, R, thin = PLAN[(b[0], side)]
			d = node_cost(Pm, b[0], b[1], side, float(n), iA, kA, kB, R, bool(thin))
			if d is None:
				raise ValueError('infeasible plan entry %s %s' % (b, side))
			nA = d['nA']
			pm, wm = miss_terms()
			psi = kA * (phimax('f', nA) + phimax('g', nA)) + kB * (phimax('f', n) + phimax('g', n)) + d['miss_children'] * pm
			if psi > Pm.PSI_MAX:
				raise ValueError('psi too large at %s %s' % (b, side))
			retry = 1.0 / (1.0 - psi)
			sub = kA * (wmax('f', nA) + wmax('g', nA)) + kB * (wmax('f', n) + wmax('g', n)) + d['miss_children'] * wm
			W[(b, side)] = d['own'] + retry * sub
			window = (kA if nA > Pm.LX + 1e-9 else 0) + (kB if n > nA + 1e-9 else 0)
			PHI[(b, side)] = phi_nb(d['mu'] * (1 - psi), Pm, window)
			DET[(b, side)] = dict(d, psi=psi, retry=retry, sub_log2=math.log2(sub) if sub > 0 else None,
								 W_log2=math.log2(W[(b, side)]), phi=PHI[(b, side)])
	n1, iA, kA, Rb, nrels = BOOT   # nrels distinct boot relations are kept and tried in turn
	if n1 > Pm.LMAX:
		raise ValueError('boot window above the last bin')
	Sb = math.log2(CELLS[Rb]); cb = (Pm.LOG2P6 + Sb) / 12
	h = (n1 - Pm.LX) / Pm.MESH; nA = Pm.LX + iA * h
	m_max = 12 * cb + Pm.LOG2U6 + max(x[1] for x in Pm.BOOT_BINS)
	pr = sum(w * ss_lower2(12 * cb + Pm.LOG2U6 + cg, Pm.LX, float(n1), iA, kA, 1, Pm.MESH, Pm.ALPHA_G, Pm.TRUNC)
			 for w, cg in Pm.BOOT_BINS)
	det = 1.0 - (1 - Pm.ECM_SUCCESS) * (m_max / 20.0 + kA + 1)
	trials = nrels / (pr * det * Pm.UDUP * Pm.COUNTED)
	balls = 1.0 + trials / CELLS[Rb]
	trial = Pm.NORMS + ecm_instr(nA, m_max, Pm.ECM_SUCCESS) + Pm.PRIME + Pm.IDEAL
	pm, wm = miss_terms(); mb = missing(Pm, 0, m_max)
	psi = kA * phimax('g', nA) + phimax('g', float(n1)) + mb * pm
	if psi > Pm.PSI_MAX:
		raise ValueError('boot psi too large')
	retry = 1 / (1 - psi)
	sub = kA * wmax('g', nA) + wmax('g', float(n1)) + mb * wm
	one_log = Pm.LLL_NODE * balls + trials * trial + retry * sub
	phi_boot = psi ** nrels * (1 + Pm.TABLE_SLACK)   # all nrels kept relations fail, each with probability <= psi given any table in the event E (table_tails; DESCENT.md)
	boot = dict(n1=n1, nA=nA, kA=kA, R=Rb, nrels=nrels, S=Sb, m=m_max, pr_log2=math.log2(pr), trials_log2=math.log2(trials),
				trial_log2=math.log2(trial), balls=balls, miss=mb, own_log2=math.log2(trials * trial), sub_log2=math.log2(sub),
				psi=psi, retry=retry, phi=phi_boot, one_log_log2=math.log2(one_log))
	return dict(W=W, PHI=PHI, DET=DET, bins=bins, boot=boot, two_log_units=math.log2(2 * one_log / UNIT),
				one_log_instr=one_log)
PLAN16 = {(int(k[:-1]), k[-1]): tuple(float(x) if "." in x else int(x) for x in v.split(","))   # (bin, side): (n, iA x 8, kA, kB, R, discard); the plan of 6895f2e2 re-searched at these inputs (DESCENT.md)
					for k, v in (e.split(":") for e in (
 "69f:69,368,0,0,288,1 69g:69,384,0,0,320,1 75f:74,268,1,0,256,1 75g:75,336,1,0,256,1 80f:75,328,1,0,320,0 80g:76,228,2,0,192,0 85f:81,188,2,0,192,1 "
 "85g:75,368,2,0,224,1 90f:77,280,2,0,256,1 90g:76,324,2,0,320,0 95f:81,192,3,0,288,1 95g:80,200,2,1,256,1 100f:85,92,1,1,320,1 100g:85,136,2,1,288,1 "
 "105f:84,140,2,1,320,1 105g:85,136,3,1,320,0 110f:90,148,2,1,288,1 110g:94,156,2,1,320,1 115f:92,180,2,1,320,1 115g:100,196,2,1,288,1 "
 "120f:94,228,2,1,320,0 120g:94,372,3,0,384,0 125f:100,192,2,1,384,1 125g:100,192,3,1,384,1 130f:100,192,3,1,320,1 130g:101,192,4,1,384,1 "
 "135f:103,216,3,1,320,1 135g:110,236,3,1,320,0 140f:105,260,3,1,320,1 140g:110,276,3,1,384,0 145f:115,208,3,1,320,1 145g:115,248,4,1,288,1 "
 "150f:110,272,3,1,384,0 150g:110,336,4,1,384,0 155f:115,216,4,1,288,0 155g:124,236,4,1,384,1 160f:115,216,4,1,384,0 160g:128,236,4,1,384,0 "
 "165f:124,196,4,1,384,0 165g:130,268,4,1,384,0 170f:124,224,4,1,384,0 170g:135,288,4,1,384,1 175f:136,240,4,1,288,1 175g:150,240,4,1,384,0 "
 "180f:130,288,4,1,320,0 180g:130,272,6,1,320,0 185f:140,240,4,1,384,0 185g:137,244,6,1,384,0 190f:145,240,4,1,384,0 190g:157,188,6,1,384,0 "
 "195f:150,240,4,1,384,0 195g:161,316,4,1,384,0 200f:150,288,4,1,320,0 200g:150,280,6,1,288,0 205f:165,240,4,1,384,1 205g:160,240,6,1,384,0 "
 "210f:160,288,4,1,320,0 210g:173,220,6,1,384,0 215f:165,288,4,1,320,0 215g:170,288,6,1,224,0 220f:170,288,4,1,384,1 220g:170,288,6,1,288,0 ").split())}
BOOT16 = (210, 356, 13, 4, 16)
def report(res):
	"""compact per-bin table for the experiment report."""
	rows = []
	for b in res['bins']:
		for side in ('f', 'g'):
			d = res['DET'][(b, side)]
			rows.append(dict(bin=b[0], side=side, n=int(d['n']), nA=round(d['nA'], 1), kA=d['kA'], kB=d['kB'], R=d['R'], ret=round(d['ret'], 4),
							 loss_f=d['loss_f'], loss_g=d['loss_g'], how=d['how'], mu=round(d['mu'], 2), own_log2=round(d['own_log2'], 3),
							 W_log2=round(d['W_log2'], 3), psi=round(d['psi'], 4), phi=round(d['phi'], 6), nf=round(d['nf'], 1), ng=round(d['ng'], 1), T=round(d['T'], 1)))
	return rows
