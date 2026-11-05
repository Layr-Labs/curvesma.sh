"""stnfs ledger model: relation collection (chunk-routed records, every block
charged; STORAGE.md), linear algebra per entry, per row and per reduction (STORAGE.md)
and the metered descent cap DESCENT_CAP (log2 units; descent16.evaluate gives"""
import math
from fractions import Fraction as F
from common12 import rho, integer_metric, quadratic, regenerate_pair, norms, mul, norm, integer_interval
import ideal_canonical as IC
import embedding_bounds as EB
import filtsim
from threshold_model import Retention
LN2 = math.log(2)
UNIT = 20000.0
OP = dict(lb=69.0, lx2=64.0, excess=1.05, slack_f=295, slack_g=226, width_bits=8, L0=64,
   theta=F(5, 16), alpha_f=5.15, alpha_g=0.85)
MARGIN = 0.06
MARGIN_SE = 2.0
Q_CORE = 0.02
DESCENT_CAP = 113.156048
NODE_CAP_LOG2 = 27.3
CANON_WORK_ALLOWANCE = 1.0047
WORK_SE_K = 2.0
Q_F_SAFETY = 100.0
ECM_FIND = 0.999
REC_ROUTE = 10.01
CELL_ROUTE = 8.0
LINE_TAB = [(18.81, 17.2), (13.82, 24.1), (9.82, 38.99), (5.82, 88.32), (1.83, 352.1)]
LINE_SE_K = 2.0
REC_TAB = [(18.91, 44.67), (10.91, 77.42), (3.91, 250.54), (0.91, 810.82), (-1.09, 1711.20)]
C_COSET_EMPTY = 1600.0
C_APPLY = 7.00
C_RECOVER = 7.003
WIDE_FROM, WIDE_REC, WIDE_SETUP = 63.0, 1.01, 1.0
REC_TAB_W = [(18.91, 23449665 / 522256), (10.91, 3216127 / 41121), (3.91, 1745730 / 6539), (0.91, 3721383 / 3967), (-1.09, 6474056 / 2901)]
LOG2_WIDE_LAT = 70.6
REC7_TAB = [(21.99, 63.72), (19.99, 71.42), (17.99, 82.58), (15.99, 92.07), (14.99, 94.65), (13.99, 104.03), (12.99, 121.05), (11.99, 136.87), (10.99, 143.12), (9.99, 177.37), (8.99, 192.16), (7.99, 255.03), (6.99, 334.57), (5.99, 417.68), (4.99, 417.68)]
P7, EXTRA7, RUN_ALLOW, STEP7 = 52.0, 2.0, 0.01, True
LOG2_PREP7 = 68.3
def step7(l):
	return next(c for h, c in REC7_TAB if l >= h)
C_SCAN_CELL, C_SCAN_LEAF = 8.0, 21.0
C_BUILD_E, C_BUILD_SLICE, C_BUILD_SEG = 28.17, 2.0 ** 21, 2.0 ** 17
C_ZERO = 0.625
C_SETUP6 = 406.87 + 72.0 + 21.13
C_NB_EVAL = 1508.2
C_NB_ROW = 4600.0
C_NB_SLICE = 17500.0
MONT5 = 367.76
C_NORMS = 2.0 ** 22
C_TDIV = 2.0 ** 24.5
NBAR = 2.0
LX0, P_MID = 8.0, 48.0
BAND_FBITS = 48.0
BAND_NBITS = (832, 576)
BAND_FMAX = (41, 28)
BAND_Q = 0.999
BAND_ROUNDS = 2
P_MIN = 256
ECM_TAB = {20: (11e3, 86), 25: (5e4, 214), 30: (25e4, 430), 35: (1e6, 904), 40: (3e6, 2350), 45: (11e6, 4480),
   50: (43e6, 7553), 55: (11e7, 17769), 60: (26e7, 42017), 65: (85e7, 69408), 70: (29e8, 102212),
   75: (76e8, 188056), 80: (25e9, 265557)}
def ecm_instr(fbits, nbits, success=ECM_FIND):
	import ecm16
	return ecm16.ecm_instr(fbits, nbits, success)
def lam_top(u, t0, m=40):
	r = filtsim.rho(u)
	if r <= 0 or t0 >= 1:
		return 0.0
	return sum(filtsim.rho(u - t) / r / t * ((1 - t0) / m) for t in (t0 + (1 - t0) * (i + 0.5) / m for i in range(m)))
def log2_li(log2x):
	x = log2x * LN2
	return log2x - math.log2(x) + math.log2(1 + 1 / x + 2 / x / x)
def tab(t, log2h):
	if log2h >= t[0][0]:
		return t[0][1]
	for (h1, c1), (h2, c2) in zip(t, t[1:]):
		if log2h >= h2:
			return c2 + (c1 - c2) * (log2h - h2) / (h1 - h2)
	return t[-1][1]
def coprime_pair(a, b):
	d = math.gcd(norm(a), norm(b))
	if d == 1:
		return True
	rows = [mul(v, [0] * k + [1]) for v in (a, b) for k in range(6)]
	key = IC.row_hnf(rows, d)
	return all(key[k][k] == 1 for k in range(6))
_CTX = {}
def pool_row(cfg, sd, i):
	if not _CTX:
		_CTX.update(qa=integer_metric(cfg['matrix_a']), qb=integer_metric(cfg['matrix_b']), disks=EB.root_disks())
	(qa, da, ia), (qb, db, ib), disks = _CTX['qa'], _CTX['qb'], _CTX['disks']
	A2 = cfg['production_A_squared']
	a, b, aq, bq = regenerate_pair(cfg, sd, i)
	qa16 = round(F(aq, da) / A2 * 65536)
	nf, ng = norms(a, b)
	nf16 = round(16 * math.log2(nf)); ng16 = round(16 * math.log2(ng))
	prim = coprime_pair(a, b)
	radius = (A2 - F(aq, da)) * db
	t = sum(qb[0][j] * b[j] for j in range(1, 6))
	c = sum(qb[j][k] * b[j] * b[k] for j in range(1, 6) for k in range(1, 6)) - radius
	iv = integer_interval(qb[0][0], t, c)
	blo = (b[0] // 64) * 64; lo, hi = max(iv.start, blo), min(iv.stop - 1, blo + 63)
	rlen = 16 * max(1, (iv.stop - iv.start) // 16); blen = hi - lo + 1
	est = EB.NormBounds(a, disks)
	box = box_rows(qb, radius, b, blo)
	leaves = EB.box_leaves(est, [(est.row(list(s)), l, h) for s, l, h in box], min(x[1] for x in box), max(x[2] for x in box))
	leaves = [l for l in leaves if l[0] <= hi and l[1] >= lo and l[3] in (None, 0)]
	at = [l for l in leaves if l[0] <= b[0] <= l[1]][0]
	return dict(seed=sd, index=i, qa=qa16 / 65536.0, nf=nf16 / 16.0, ng=ng16 / 16.0, primitive=prim, rejected=isinstance(at[2], dict),
	   nonfinite=False, evals=est.counts['block_evaluations'], bcells=sum(h - l + 1 for _, l, h in box), leaves=len(leaves), blen=blen, rlen=rlen)
SUM_KEYS = ('rows', 'ysum', 'y2', 'nfy', 'ngy', 'work', 'work2', 'w', 'nfw', 'ngw', 'cells')
def add_sums(parts):
	return {k: sum(p[k] for p in parts) for k in SUM_KEYS}
def pool_totals(s, n_eff, n_draws, geo):
	n, ys, y2, w, w2 = s['rows'], s['ysum'], s['y2'], s['work'], s['work2']
	mw = w / n
	cond = math.sqrt(max(0.0, (w2 - n * mw * mw) / (n - 1)) / n) / mw
	return dict(log2_rels=math.log2(geo.va / n_eff * geo.vb * geo.lower * ys),
	   yield_se_log2=math.sqrt(max(0.0, y2 / n_eff - (ys / n_eff) ** 2) / n_eff) / (ys / n_eff) / LN2,
	   nf_rel=s['nfy'] / ys, ng_rel=s['ngy'] / ys, nf_mean=s['nfw'] / s['w'], ng_mean=s['ngw'] / s['w'],
	   log2_population=math.log2(geo.va / n_eff * s['cells']), per_cell=w / s['cells'], log2_est=math.log2(geo.va / n_eff * w / UNIT),
	   cond_se=cond, work_se=math.sqrt(cond ** 2 + (1 - n / n_draws) / n))
BOX = (2, 2, 4)
def box_rows(qb, radius, b, blo):
	import itertools
	base = [(b[j] // s) * s for j, s in zip((2, 3, 4), BOX)]; out = []
	for d in itertools.product(*(range(s) for s in BOX)):
		x = (b[1], base[0] + d[0], base[1] + d[1], base[2] + d[2], b[5])
		t = sum(qb[0][j] * x[j - 1] for j in range(1, 6))
		c = sum(qb[j][k] * x[j - 1] * x[k - 1] for j in range(1, 6) for k in range(1, 6)) - radius
		iv = integer_interval(qb[0][0], t, c); lo, hi = max(iv.start, blo), min(iv.stop - 1, blo + 63)
		if lo <= hi:
			out.append((x, lo, hi))
	out.sort(key=lambda r: r[0] != tuple(b[1:]))
	return out
class Geometry:
	def __init__(self, cfg, theta, lower=F(2401, 2500), upper=F(2601, 2500)):
		self.A2 = cfg['production_A_squared']
		self.qa, self.da, self.ia = integer_metric(cfg['matrix_a'])
		self.qb, self.db, self.ib = integer_metric(cfg['matrix_b'])
		from ideal_canonical import determinant
		det_a = F(determinant(self.qa), self.da ** 6); det_b = F(determinant(self.qb), self.db ** 6)
		v6 = math.pi ** 3 / 6
		self.va = v6 * (float(theta) * self.A2) ** 3 / (2 * math.sqrt(float(det_a)))
		self.vb = v6 * self.A2 ** 3 / (2 * math.sqrt(float(det_b)))
		self.lower, self.upper = float(lower), float(upper)
		self.rb = [math.sqrt(float(self.ib[j][j]) * self.A2 * self.db) for j in range(6)]
		self.inv55 = float(self.ib[5][5]) * self.db
	def slice_cells(self, qa_frac):
		return self.vb * (1 - qa_frac) ** 3 * self.upper
	def lines(self, qa_frac):
		return max(1.0, math.floor(math.sqrt((1 - qa_frac) * self.A2 * self.inv55)))
	def entries_upper(self, qa_frac):
		s = math.sqrt(1 - qa_frac); L4 = 1.0
		for j in range(1, 5):
			L4 *= 2 * math.floor(self.rb[j] * s) + 3
		return (math.floor(math.sqrt((1 - qa_frac) * self.A2 * self.inv55)) + 2) * L4
	def suffix_entries(self, qa_frac):
		s = math.sqrt(1 - qa_frac); L4 = 1.0
		for j in range(1, 5):
			L4 *= 2 * math.floor(self.rb[j] * s) + 1
		return L4
def evaluate(rows, n_total, geo, op=OP, retain=None, cost=None, scale_consts=1.0, verbose=False, return_work=False):
	lb, lx2, excess = op['lb'], op['lx2'], op['excess']
	Sf, Sg, W = op['slack_f'], op['slack_g'], op['width_bits']
	retain = retain or Retention(int(lb), int(lx2), step=.25, max_bits=850, tilt=.05)
	cost = cost or Retention(850, int(lx2), step=.25, max_bits=851, tilt=.005)
	wsum = 0.0; ysum = 0.0; y2 = 0.0; rej_w = 0.0; nf_m = ng_m = 0.0; nf_r = ng_r = 0.0
	work_sum = work2 = 0.0
	work_list = [] if return_work else None
	c_ecm_f, c_ecm_g = ecm_instr(lb, Sf + 8), ecm_instr(lb, Sg + 8)
	c_band_f = BAND_ROUNDS * ecm_instr(BAND_FBITS, BAND_NBITS[0], BAND_Q)
	c_band_g = BAND_ROUNDS * ecm_instr(BAND_FBITS, BAND_NBITS[1], BAND_Q)
	band_keep = (1 - BAND_FMAX[0] * (1 - BAND_Q) ** BAND_ROUNDS) * (1 - BAND_FMAX[1] * (1 - BAND_Q) ** BAND_ROUNDS)
	pairs_all = (NBAR / 2) * 2 ** log2_li(lx2) + ((WIDE_SETUP - 1) * (NBAR / 2) * (2 ** log2_li(lx2) - 2 ** log2_li(WIDE_FROM)) if lx2 > WIDE_FROM else 0.0)
	for r in rows:
		w = (1 - r['qa']) ** 3
		wsum += w
		if r['primitive'] and not r['rejected'] and r['qa'] >= Q_CORE:
			rf = retain.at(r['nf'], Sf)['retention']; rg = retain.at(r['ng'], Sg)['retention']
			uf_r, ug_r = (r['nf'] + op['alpha_f']) / lb, (r['ng'] + op['alpha_g']) / lb
			miss = ECM_FIND ** (lam_top(uf_r, lx2 / lb) + lam_top(ug_r, lx2 / lb))
			y = rho(uf_r) * rho(ug_r) * rf * rg * miss * band_keep
			ysum += w * y; y2 += w * w * y * y
			nf_r += w * y * r['nf']; ng_r += w * y * r['ng']
		if r['rejected']:
			rej_w += w
		nf_m += w * r['nf']; ng_m += w * r['ng']
		S = geo.slice_cells(r['qa']); K = geo.lines(r['qa'])
		p_all = dict(build=(C_BUILD_E * geo.entries_upper(r['qa']) + C_BUILD_SLICE) / S + C_BUILD_SEG * 4 / 2 ** 32, zero=0.0, setup6=pairs_all * C_SETUP6 / S,
		    thr=C_NB_EVAL * r['evals'] / r['bcells'] + C_NB_ROW / r['rlen'] + C_NB_SLICE / S)
		lp = LX0; line = 0.0
		while lp < P_MID:
			lp2 = min(lp + 0.25, P_MID); dsum = math.log(lp2 / lp)
			line += NBAR * dsum * tab(LINE_TAB, math.log2(S / K) - 0.5 * (lp + lp2)); lp = lp2
		lp = P_MID; rec = 0.0; empty = 0.0; save7 = 0.0; lg = math.log2(r.get('gain', 1.0))
		while lp < lx2:
			lp2 = min(lp + 0.25, lx2); dsum = math.log(lp2 / lp)
			log2h = math.log2(S) - 0.5 * (lp + lp2)
			if log2h >= REC_TAB[-1][0]:
				c6 = NBAR * dsum * (WIDE_REC * tab(REC_TAB_W, log2h) if lp >= WIDE_FROM else tab(REC_TAB, log2h)); rec += c6
				if P7 <= lp < WIDE_FROM and log2h + lg >= REC7_TAB[-1][0]:
					dp = (NBAR / 2) * (2 ** log2_li(lp2) - 2 ** log2_li(lp))
					c7 = NBAR * dsum * ((step7(log2h + lg) if STEP7 else tab(REC7_TAB, log2h + lg)) + EXTRA7) + dp * C_SETUP6 * (2 ** -lg - 1) / S
					save7 += max(0.0, c6 - c7)
			else:
				dpairs = (NBAR / 2) * (2 ** log2_li(lp2) - 2 ** log2_li(lp))
				empty += dpairs * C_COSET_EMPTY / S * (WIDE_REC if lp >= WIDE_FROM else 1.0)
			lp = lp2
		rec -= (1 - RUN_ALLOW) * save7
		hits = NBAR * math.log(lx2 / P_MID); hits_w = NBAR * math.log(lx2 / max(P_MID, min(lx2, WIDE_FROM)))
		sf = cost.survivor_upper(r['nf'], Sf, W); sg = cost.survivor_upper(r['ng'], Sg, W)
		q_f = min(1.0, Q_F_SAFETY * rho((r['nf'] + op['alpha_f']) / lb) / cost.survivor_upper(r['nf'], Sf, 0))
		cbf = BAND_ROUNDS * ecm_instr(BAND_FBITS, min(BAND_NBITS[0], math.ceil(r['nf']) + 1), BAND_Q)
		cbg = BAND_ROUNDS * ecm_instr(BAND_FBITS, min(BAND_NBITS[1], math.ceil(r['ng']) + 1), BAND_Q)
		c_surv = C_NORMS + C_TDIV + cbf + c_ecm_f + q_f * (cbg + c_ecm_g)
		p_keep = dict(line=line, rec=rec, empty=empty,
		    apply=(hits + (WIDE_REC - 1) * hits_w) * C_APPLY, recover=hits * C_RECOVER,
		    route=CELL_ROUTE + REC_ROUTE * hits,
		    scan=C_SCAN_CELL + C_SCAN_LEAF * r['leaves'] / r['blen'],
		    zero=C_ZERO, cof=sf * sg * c_surv)
		per_cell = sum(p_all.values()) + sum(p_keep.values())
		per_cell = per_cell * scale_consts if scale_consts != 1.0 else per_cell
		contrib = S * per_cell
		work_sum += contrib; work2 += contrib * contrib
		if return_work: work_list.append(contrib)
	n_a_upper = geo.va / n_total
	population = n_a_upper * sum(geo.slice_cells(r['qa']) for r in rows)
	relcol_instr = n_a_upper * work_sum * CANON_WORK_ALLOWANCE
	per_cell_mean = relcol_instr / population
	rels = n_a_upper * geo.vb * geo.lower * ysum
	yield_se = math.sqrt(max(0.0, y2 / n_total - (ysum / n_total) ** 2) / n_total) / (ysum / n_total) / LN2 if ysum > 0 else float('inf')
	nf_m /= wsum; ng_m /= wsum
	log2_U = math.log2((NBAR / 2) * 2 ** log2_li(lb) + 2 ** 31)
	log2_need = log2_U + math.log2(excess)
	log2_rels = math.log2(rels) if rels > 0 else float('-inf')
	feasible = log2_rels >= log2_need + max(MARGIN, MARGIN_SE * yield_se)
	log2_relcol = math.log2(relcol_instr / UNIT)
	out = dict(op={k: (str(v) if isinstance(v, F) else v) for k, v in op.items()},
	   log2_population=math.log2(population), per_cell=per_cell_mean,
	   log2_relcol=log2_relcol,
	   log2_rels=log2_rels, log2_U=log2_U, log2_need=log2_need,
	   margin=log2_rels - log2_need, yield_se_log2=yield_se, feasible=feasible,
	   rejected_weight=rej_w / wsum, nf_mean=nf_m, ng_mean=ng_m, nf_rel=nf_r / ysum, ng_rel=ng_r / ysum,
	   records_per_cell=NBAR * math.log(lx2 / P_MID),
	   c_surv_log2=math.log2(C_NORMS + C_TDIV + c_band_f + c_band_g + c_ecm_f + c_ecm_g), mean_slice_log2=math.log2(population / (n_a_upper * len(rows))))
	out['sums'] = dict(rows=len(rows), ysum=ysum, y2=y2, nfy=nf_r, ngy=ng_r, work=work_sum, work2=work2, w=wsum, nfw=nf_m * wsum, ngw=ng_m * wsum,
	       cells=sum(geo.slice_cells(r['qa']) for r in rows))
	if return_work: out['work_contribs'] = work_list
	return out
FILT_RUNS = {
 12: ((1,375,263,99),(2,374,262,103),(3,361,250,142),(4,387,269,163),(5,405,282,97),(6,386,265,113),(7,469,322,97),(8,405,275,141),(9,489,344,143),(10,369,258,185)),
 13: ((1,516,358,161),(2,557,394,152),(3,424,297,151),(4,447,308,129),(5,511,354,142),(6,652,450,164)),
 14: ((1,1074,733,140),(2,1096,755,123),(3,991,678,133),(4,938,637,151),(5,1019,693,131),(6,1170,810,171)),
 15: ((1,1762,1211,188),(2,1921,1312,186),(3,1831,1262,250),(4,1837,1259,165),(5,1856,1293,126),(6,1970,1380,191)),
 16: ((1,3843,2638,199),(2,3554,2423,141),(3,4127,2834,185),(4,3672,2508,196),(5,3719,2512,162),(6,3825,2622,155),(7,3908,2666,142),(8,4073,2798,166),(9,3903,2676,134),(10,3891,2685,150),(11,4169,2834,154),(12,4146,2856,175)),
 17: ((1,7799,5288,200),(2,7309,4976,186),(3,7380,5046,161),(4,7026,4822,143),(5,7212,4925,179),(6,7226,4949,192),(7,7916,5426,212),(8,7762,5300,166),(9,7880,5385,181),(10,7025,4826,186),(11,7650,5191,187),(12,7277,5029,186),(13,7596,5206,188),(14,7251,4965,176),(15,7056,4823,158),(16,7120,4862,169)),
 18: ((1,14271,9743,202),(2,13785,9398,224),(3,13521,9299,197),(4,14147,9651,180),(5,13754,9387,197),(6,13558,9275,179),(7,13957,9527,173),(8,14279,9703,206)),
}
NU300 = 0.03943
W_ROW = 301.0
C_ENTRY = 20.9
C_ENTRY_K = 13.21
C_ROW = 2710.0
C_RED9 = 420.0
UNIT_FRAC_CHARGED = 0.87
BW_M = 28
BW_N = 5
C_PROJ = 60.0
N_SM = 21
C_FULL = 97.0
LOG2_SM = 88.0
LOG2_FBC = 72.0
LOG2_CANON = 90.0
LOG2_MISC = 80.0
LOG2_VERIFY = 11.0
def linear_algebra(log2_U, excess=None):
	ratio = NU300
	log2_N = log2_U + math.log2(NU300)
	c_prod = W_ROW * C_ENTRY + C_ROW + N_SM * C_FULL
	c_kry = W_ROW * C_ENTRY_K + C_ROW + N_SM * C_FULL
	per_n2 = c_kry * (BW_N / BW_M + 1) + (1 / BW_M + 1 / BW_N) * (c_prod + BW_N * C_FULL + C_RED9) + C_PROJ * (BW_M + BW_N)
	per_n1 = (BW_M + BW_N) * C_RED9
	log2_la = math.log2((2 ** (2 * log2_N) * per_n2 + 2 ** log2_N * per_n1) / UNIT)
	return dict(ratio=ratio, log2_N=log2_N, log2_la=log2_la, log2_gen=generator(log2_N))
def generator(log2_N):
	k = BW_M + BW_N
	log2_L = math.log2((2 ** log2_N + N_SM) * (1 / BW_M + 1 / BW_N) + 2)
	return math.log2(4 * k ** 3 * (C_FULL + C_RED9) / UNIT) + math.log2(3) * log2_L + math.log2(log2_L)
def total_log2(log2_relcol, log2_la, log2_desc=DESCENT_CAP, log2_gen=0.0):
	return math.log2(2 ** log2_relcol + 2 ** log2_la + 2 ** log2_desc + 2 ** log2_gen + 2 ** LOG2_SM + 2 ** LOG2_FBC
	     + 2 ** LOG2_CANON + 2 ** LOG2_MISC + 2 ** LOG2_VERIFY + 2 ** LOG2_WIDE_LAT + 2 ** LOG2_PREP7 + 2 ** 50)
