"""stnfs ledger model: collection, linear algebra, the metered descent cap."""
import math, json, struct, zlib, lzma, base64
from fractions import Fraction as F
from common12 import rho, integer_metric, quadratic
import filtsim
from threshold_model import Retention
LN2 = math.log(2)
UNIT = 20000.0
OP = dict(lb=69.0, lx2=64.0, excess=1.1, slack_f=295, slack_g=224, width_bits=8, L0=64,
		 theta=F(5, 16), alpha_f=5.15, alpha_g=0.85)   # pool theta (9ff4fc85's pool); collected sub-region: config production_theta_region
MARGIN = 0.06      # feasibility: margin >= max(MARGIN, MARGIN_SE x pool SE) (proof.md 3)
MARGIN_SE = 2.9
Q_CORE = 0.02      # rows with Q_a(a)/A^2 < Q_CORE credit no relations (cells still charged)
POOL_CHUNKS = 16   # part=pool regenerates one chunk of the canonical list and checks its sums
Q_F_SAFETY = 100.0   # early abort: g-side charged at 100x the H1 P(f factors) (proof.md 3)
PAD_N = 1.45   # x the 2^19 worst: >= 1.03 x every weight-100 N/U from 2^14 up, both means (proof.md 4)
DESCENT_CAP = 116.495   # meter: E[two-log] 114.7514 + Markov gap for modelled success 0.70 (proof.md 7)
NODE_CAP_LOG2 = 28.5    # node counter for both logs (DESCENT.md)
CANON_WORK_ALLOWANCE = 1.0167   # >= 1 + WORK_SE_K x work-estimator rel SE ([4a])
WORK_SE_K = 2.0
ECM_FIND = 0.999
REC_ROUTE = 10.01        # record route (STORAGE.md)
CELL_ROUTE = 8.0        # instr/cell (bookkeeping, sentinels)
LINE_TAB = [(13.82, 24.10), (9.82, 31.07), (5.82, 87.57), (1.82, 298.78)]
REC_TAB = [(18.91, 44.67), (10.91, 77.42), (3.91, 250.54), (0.91, 810.82), (-1.09, 1711.20)]
C_COSET_EMPTY = 1600.0
C_APPLY = 7.00
C_RECOVER = 7.003
WIDE_FROM, WIDE_FACTOR = 63.0, 2.0  # STORAGE.md
C_SCAN_CELL, C_SCAN_LEAF = 8.0, 21.0
C_BUILD = 12.325
C_ZERO = 0.625
C_SETUP6 = 406.87 + 72.0 + 21.13
C_NB_EVAL = 1508.2
C_NB_ROW = 4300.0
C_NB_SLICE = 17500.0
MONT5 = 367.76
C_NORMS = 2.0 ** 22
C_TDIV = 2.0 ** 24.5
NBAR = 2.0
LX0, P_MID = 8.0, 48.0   # line sieve to 2^48 (hits/coset stays inside the measured line table); record walk above
BAND_FBITS = 48.0   # survivor band strip recovers line-credited (2^20, 2^48] primes
BAND_NBITS = (832, 576)   # unstripped cofactor caps (bits)
BAND_FMAX = (41, 28)   # max band factors per side: floor(nbits/20)
BAND_Q = 0.999   # per-round per-prime find probability
BAND_ROUNDS = 2
def ecm_instr(fbits, nbits, success=ECM_FIND):
	"""the ecm16 schedule ([12]): cheapest (B1,B2) with per-prime success >= success."""
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
def pool_lists(text):
	"""pool12.b85 (9ff4fc85): per-seed canonical-index delta-varint lists plus 64 field anchors;
	row fields regenerate from (seed, index) (evidence16.regen_row)."""
	raw = lzma.decompress(base64.b85decode(''.join(text.split())), format=lzma.FORMAT_RAW,
						 filters=[{"id": lzma.FILTER_LZMA2}])
	hdr, blob = raw.split(b'\n', 1)
	hdr = json.loads(hdr)
	def varints(buf):
		out = []; n = 0; s = 0
		for byte in buf:
			n |= (byte & 0x7F) << s
			if byte & 0x80: s += 7
			else: out.append(n); n = 0; s = 0
		assert s == 0 and n == 0
		return out
	o = 0; items = []
	for sd in hdr['seeds']:
		prev = -1
		for d in varints(blob[o:o + hdr[f'seed{sd}_bytes']]):
			prev += d; items.append((sd, prev))
		o += hdr[f'seed{sd}_bytes']
	fmt = '<HHHHBBBBB'; sz = struct.calcsize(fmt); ab = blob[o:o + hdr['anchor_bytes']]
	assert len(ab) == hdr['anchor_rows'] * sz and o + hdr['anchor_bytes'] == len(blob) and len(items) == hdr['canonical_rows']
	hdr['anchors'] = [dict(zip(('position', 'qa16', 'nf16', 'ng16', 'flags', 'evals', 'leaves', 'blen', 'rlen16'),
							   struct.unpack(fmt, ab[k * sz:(k + 1) * sz]))) for k in range(hdr['anchor_rows'])]
	return hdr, items
def chunk_items(items, c):
	L = -(-len(items) // POOL_CHUNKS)
	return items[c * L:(c + 1) * L]
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
	def suffix_entries(self, qa_frac):
		s = math.sqrt(1 - qa_frac); L4 = 1.0
		for j in range(1, 5):
			L4 *= 2 * math.floor(self.rb[j] * s) + 1
		return L4
SUM_KEYS = ('n', 'wsum', 'ysum', 'y2', 'rej_w', 'nf_m', 'ng_m', 'nf_r', 'ng_r', 'work', 'work2', 'q_sum', 'slices')
def accumulate(rows, geo, op=OP, retain=None, cost=None):
	"""per-row sums of the collection estimator (SUM_KEYS); rows with Q_a/A^2 < Q_CORE credit no relations."""
	lb, lx2 = op['lb'], op['lx2']
	Sf, Sg, W = op['slack_f'], op['slack_g'], op['width_bits']
	retain = retain or Retention(int(lb), int(lx2), step=.25, max_bits=850, tilt=.05)
	cost = cost or Retention(850, int(lx2), step=.25, max_bits=851, tilt=.005)
	t = dict.fromkeys(SUM_KEYS, 0.0)
	c_ecm_f, c_ecm_g = ecm_instr(lb, Sf + 8), ecm_instr(lb, Sg + 8)
	c_band_f = BAND_ROUNDS * ecm_instr(BAND_FBITS, BAND_NBITS[0], BAND_Q)
	c_band_g = BAND_ROUNDS * ecm_instr(BAND_FBITS, BAND_NBITS[1], BAND_Q)
	q_cache = {}
	band_keep = (1 - BAND_FMAX[0] * (1 - BAND_Q) ** BAND_ROUNDS) * (1 - BAND_FMAX[1] * (1 - BAND_Q) ** BAND_ROUNDS)
	pairs_all = (NBAR / 2) * 2 ** log2_li(lx2) + ((WIDE_FACTOR - 1) * (NBAR / 2) * (2 ** log2_li(lx2) - 2 ** log2_li(WIDE_FROM)) if lx2 > WIDE_FROM else 0.0)
	for r in rows:
		w = (1 - r['qa']) ** 3
		t['n'] += 1; t['wsum'] += w
		if r['primitive'] and not r['rejected'] and r['qa'] >= Q_CORE:
			rf = retain.at(r['nf'], Sf)['retention']; rg = retain.at(r['ng'], Sg)['retention']
			uf_r, ug_r = (r['nf'] + op['alpha_f']) / lb, (r['ng'] + op['alpha_g']) / lb
			miss = ECM_FIND ** (lam_top(uf_r, lx2 / lb) + lam_top(ug_r, lx2 / lb))
			y = rho(uf_r) * rho(ug_r) * rf * rg * miss * band_keep
			t['ysum'] += w * y; t['y2'] += w * w * y * y
			t['nf_r'] += w * y * r['nf']; t['ng_r'] += w * y * r['ng']
		if r['rejected']:
			t['rej_w'] += w
		t['nf_m'] += w * r['nf']; t['ng_m'] += w * r['ng']
		S = geo.slice_cells(r['qa']); K = geo.lines(r['qa'])
		p_all = dict(build=C_BUILD, zero=0.0, setup6=pairs_all * C_SETUP6 / S,
					 thr=C_NB_EVAL * r['evals'] / r['blen'] + C_NB_ROW / r['rlen'] + C_NB_SLICE / S)
		lp = LX0; line = 0.0
		while lp < P_MID:
			lp2 = min(lp + 0.25, P_MID); dsum = math.log(lp2 / lp)
			line += NBAR * dsum * tab(LINE_TAB, math.log2(S / K) - 0.5 * (lp + lp2)); lp = lp2
		lp = P_MID; rec = 0.0; empty = 0.0
		while lp < lx2:
			lp2 = min(lp + 0.25, lx2); dsum = math.log(lp2 / lp)
			log2h = math.log2(S) - 0.5 * (lp + lp2)
			if log2h >= REC_TAB[-1][0]:
				rec += NBAR * dsum * tab(REC_TAB, log2h) * (WIDE_FACTOR if lp >= WIDE_FROM else 1.0)
			else:
				dpairs = (NBAR / 2) * (2 ** log2_li(lp2) - 2 ** log2_li(lp))
				empty += dpairs * C_COSET_EMPTY / S
			lp = lp2
		hits = NBAR * math.log(lx2 / P_MID)
		sf = cost.survivor_upper(r['nf'], Sf, W); sg = cost.survivor_upper(r['ng'], Sg, W)
		if r['nf'] not in q_cache: q_cache[r['nf']] = min(1.0, Q_F_SAFETY * rho((r['nf'] + op['alpha_f']) / lb) / cost.survivor_upper(r['nf'], Sf, 0))
		q_f = q_cache[r['nf']]
		c_surv = C_NORMS + C_TDIV + c_band_f + c_ecm_f + q_f * (c_band_g + c_ecm_g)
		t['q_sum'] += q_f
		p_keep = dict(line=line, rec=rec, empty=empty,
					 apply=hits * C_APPLY, recover=hits * C_RECOVER,
					 route=CELL_ROUTE + REC_ROUTE * hits,
					 scan=C_SCAN_CELL + C_SCAN_LEAF * r['leaves'] / r['blen'],
					 zero=C_ZERO, cof=sf * sg * c_surv)
		contrib = S * (sum(p_all.values()) + sum(p_keep.values()))
		t['work'] += contrib; t['work2'] += contrib * contrib; t['slices'] += S
	return t
def finalize(t, n_eff, n_draws, geo, op=OP):
	"""relation collection from the pool sums: relations, yield SE, work, work SE, feasibility."""
	lb, lx2, excess = op['lb'], op['lx2'], op['excess']
	c_ecm_f, c_ecm_g = ecm_instr(lb, op['slack_f'] + 8), ecm_instr(lb, op['slack_g'] + 8)
	c_band_f = BAND_ROUNDS * ecm_instr(BAND_FBITS, BAND_NBITS[0], BAND_Q)
	c_band_g = BAND_ROUNDS * ecm_instr(BAND_FBITS, BAND_NBITS[1], BAND_Q)
	n = t['n']; ysum = t['ysum']
	n_a_upper = geo.va / n_eff
	population = n_a_upper * t['slices']
	relcol_instr = n_a_upper * t['work'] * CANON_WORK_ALLOWANCE
	rels = n_a_upper * geo.vb * geo.lower * ysum
	yield_se = math.sqrt(max(0.0, t['y2'] / n_eff - (ysum / n_eff) ** 2) / n_eff) / (ysum / n_eff) / LN2 if ysum > 0 else float('inf')
	mw = t['work'] / n; sw2 = max(0.0, (t['work2'] - n * mw * mw) / (n - 1))
	cond_se = math.sqrt(sw2 / n) / mw; pcan = n / n_draws
	work_se = math.sqrt(cond_se ** 2 + (1 - pcan) / n)
	log2_U = math.log2((NBAR / 2) * 2 ** log2_li(lb) + 2 ** 31)
	log2_need = log2_U + math.log2(excess)
	log2_rels = math.log2(rels) if rels > 0 else float('-inf')
	margin = log2_rels - log2_need
	return dict(op={k: (str(v) if isinstance(v, F) else v) for k, v in op.items()},
				log2_population=math.log2(population), per_cell=relcol_instr / population,
				log2_relcol=math.log2(relcol_instr / UNIT), log2_rels=log2_rels, log2_U=log2_U, log2_need=log2_need,
				margin=margin, yield_se_log2=yield_se, feasible=margin >= max(MARGIN, MARGIN_SE * yield_se),
				rejected_weight=t['rej_w'] / t['wsum'], nf_mean=t['nf_m'] / t['wsum'], ng_mean=t['ng_m'] / t['wsum'],
				nf_rel=t['nf_r'] / ysum, ng_rel=t['ng_r'] / ysum, records_per_cell=NBAR * math.log(lx2 / P_MID),
				c_surv_log2=math.log2(C_NORMS + C_TDIV + c_band_f + c_ecm_f + (t['q_sum'] / n) * (c_band_g + c_ecm_g)),
				q_f_mean=t['q_sum'] / n, mean_slice_log2=math.log2(t['slices'] / n), rows=int(n),
				work_cond_se=cond_se, work_se=work_se)
FILT_TABLE = {   # 2^19 basis, excess 1.1 (seed, U/N, max coef, weight)
	1.1: [(1, 18.923, 40, 100.003), (2, 19.223, 40, 100.005), (3, 19.205, 41, 100.007), (4, 18.955, 40, 100.007), (5, 18.900, 42, 100.008), (6, 19.364, 40, 100.004), (7, 19.259, 46, 100.007), (8, 18.483, 34, 100.005), (9, 19.458, 42, 100.009), (10, 18.836, 40, 100.002), (11, 18.480, 37, 100.003), (12, 18.856, 54, 100.004), (13, 19.115, 38, 100.001), (14, 19.324, 37, 100.002), (15, 18.694, 40, 100.001), (16, 19.073, 50, 100.009)],
}
SCALE_NU11 = {   # N/U at weight 100, excess 1.1, by log2 U
	15: (0.05865, 0.06180, 0.06326, 0.06168, 0.05591, 0.06189, 0.06354, 0.07056, 0.06589, 0.06113, 0.07022, 0.06522, 0.05469, 0.05893, 0.05853, 0.06046, 0.06360, 0.06503, 0.06357, 0.06339, 0.06757, 0.06445, 0.06299, 0.06870, 0.05713, 0.06363, 0.05652, 0.06546, 0.06348, 0.06363,),
	16: (0.06375, 0.05461, 0.06076, 0.05553,),
	17: (0.05865, 0.05482,),
	18: (0.05404, 0.05314,),
	19: (0.05284, 0.05202, 0.05207, 0.05276, 0.05291, 0.05164, 0.05192, 0.05410, 0.05139, 0.05309, 0.05411, 0.05303, 0.05231, 0.05175, 0.05349, 0.05243,),
	20: (0.05212, 0.04971, 0.04972,),
}
SCALE_NU11_SHIP = {14: (0.06445, 0.06342, 0.06000, 0.05896, 0.06128, 0.06519,), 15: (0.05713, 0.06000, 0.06110, 0.05719, 0.05475, 0.05859, 0.06726, 0.06717, 0.06467, 0.06204, 0.07401, 0.06223, 0.05884, 0.06100, 0.05804, 0.06030, 0.06116, 0.06375, 0.05792, 0.05750, 0.06711, 0.06500, 0.06354, 0.06628, 0.05386, 0.06027, 0.05325, 0.06113, 0.06445, 0.06543,), 16: (0.05913, 0.05484, 0.05907, 0.05489,), 17: (0.05515, 0.05424,), 18: (0.05332, 0.05239,), 19: (0.05180, 0.04967, 0.04980, 0.05207,)}   # N/U at weight 100, excess 1.1, at the shipped relation-weighted means 740.34/498.98 (seeds as SCALE_NU11)
W200_11_SHIP = {14: (0.7756, 0.7786, 0.7772, 0.7836, 0.7819, 0.7800,), 15: (0.7788, 0.7752, 0.7857, 0.7855, 0.7731, 0.7776, 0.7718, 0.7746, 0.7739, 0.7757, 0.7720, 0.7724, 0.7733, 0.7779, 0.7787, 0.7910, 0.7879, 0.7731, 0.7713, 0.7781, 0.7740, 0.7737, 0.7733, 0.7744, 0.7807, 0.7777, 0.7765, 0.7698, 0.7770, 0.7682,), 16: (0.7706, 0.7760, 0.7737, 0.7815,), 17: (0.7727, 0.7744,), 18: (0.7759, 0.7748,), 19: (0.7719, 0.7716, 0.7715, 0.7730,)}   # N(200)/N(100) at the shipped means
W200_WEIGHT_MAX = 201.30   # largest mean row weight after the weight-200 merge in any run (shipped means)
W200_11 = {15: (0.7815, 0.7827, 0.7795, 0.7798, 0.7817, 0.7781, 0.7757, 0.7790, 0.7767, 0.7703,), 16: (0.7736, 0.7776, 0.7727,), 17: (0.7733, 0.7727,)}   # N(200)/N(100), excess 1.1
W_ROW = 202.0   # row weight bound (target 200; >= W200_WEIGHT_MAX): entries = N W_ROW
W200_RATIO = 0.80   # charged; max measured 0.7827
W200 = {12: (0.7831, 0.7994, 0.7698, 0.7993, 0.7752), 13: (0.7749, 0.7759, 0.8061, 0.7927, 0.7844), 14: (0.7891, 0.7927, 0.7873, 0.7798, 0.7825), 15: (0.7816, 0.7866, 0.7805, 0.7828, 0.7812), 16: (0.7780, 0.7860, 0.7843, 0.7864, 0.7866), 17: (0.7813, 0.7817, 0.7830, 0.7860, 0.7805), 18: (0.7795, 0.7783, 0.7786)}   # N(200)/N(100), excess 1.2
C_ENTRY = 20.7   # 4+3+0.88*11+0.12*33 = 20.64; single-vector (Horner)
C_ENTRY_K = 17.9   # Krylov, five vectors/pass: 4+3+0.88*39/5+0.12*33 = 17.82
C_PROJ = 60.0   # projections, 128-bit X: bw_dot2, counted 59
C_ROW = 1750.0   # eleven-owner row tail, counted 1,745
C_RED9 = 420.0   # 9-limb reduction + add, counted 385
UNIT_FRAC_CHARGED = 0.88   # H4u, below every measured value (0.917..0.948)
BW_M = 28
BW_N = 5
N_SM = 21                   # Schirokauer columns: odd unit ranks 18 + 3 ([13])
C_FULL = 97.0
LOG2_SM = 88.0
LOG2_FBC = 72.0
LOG2_CANON = 90.0
LOG2_MISC = 80.0
LOG2_VERIFY = 66.0   # checks, map, pairings, [a]P0, output
def ratio_N_over_U(excess):
	xs = sorted(FILT_TABLE)
	excess = min(max(excess, xs[0]), xs[-1])
	for a, b in zip(xs, xs[1:]):
		if excess <= b:
			ra = max(1 / r[1] for r in FILT_TABLE[a]); rb = max(1 / r[1] for r in FILT_TABLE[b])
			return ra + (rb - ra) * (excess - a) / (b - a)
	return max(1 / r[1] for r in FILT_TABLE[xs[-1]])
def linear_algebra(log2_U, excess, pad_n=PAD_N):
	ratio = ratio_N_over_U(excess)
	log2_N = log2_U + math.log2(ratio * pad_n * W200_RATIO)
	# block Wiedemann on A = [M | S] (STORAGE.md)
	c_prod = W_ROW * C_ENTRY + C_ROW + N_SM * C_FULL
	c_kry = W_ROW * C_ENTRY_K + C_ROW + N_SM * C_FULL   # five Krylov vectors share each pass
	per_n2 = c_kry * (BW_N / BW_M + 1) + (1 / BW_M + 1 / BW_N) * (c_prod + BW_N * C_FULL + C_RED9) + C_PROJ * (BW_M + BW_N)
	per_n1 = (BW_M + BW_N) * C_RED9
	log2_la = math.log2((2 ** (2 * log2_N) * per_n2 + 2 ** log2_N * per_n1) / UNIT)
	return dict(ratio=ratio, log2_N=log2_N, log2_la=log2_la, log2_gen=generator(log2_N))
def generator(log2_N):
	# (m+n) x (m+n) matrix generator (STORAGE.md)
	k = BW_M + BW_N
	log2_L = math.log2((2 ** log2_N + N_SM) * (1 / BW_M + 1 / BW_N) + 2)
	return math.log2(4 * k ** 3 * (C_FULL + C_RED9) / UNIT) + math.log2(3) * log2_L + math.log2(log2_L)
def total_log2(log2_relcol, log2_la, log2_desc=DESCENT_CAP, log2_gen=0.0):
	return math.log2(2 ** log2_relcol + 2 ** log2_la + 2 ** log2_desc + 2 ** log2_gen + 2 ** LOG2_SM + 2 ** LOG2_FBC
					 + 2 ** LOG2_CANON + 2 ** LOG2_MISC + 2 ** LOG2_VERIFY + 2 ** 50)
