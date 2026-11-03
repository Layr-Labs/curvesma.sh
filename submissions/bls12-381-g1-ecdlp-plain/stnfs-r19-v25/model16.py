"""stnfs ledger model: collection (band strip)"""
import math, json, struct, zlib, lzma, base64
from fractions import Fraction as F
from common12 import rho, integer_metric, quadratic
import filtsim
from threshold_model import Retention
LN2 = math.log(2)
UNIT = 20000.0
OP = dict(lb=69.0, lx2=64.0, excess=1.3, slack_f=297, slack_g=233, width_bits=8, L0=64,
  theta=F(19, 50), alpha_f=5.0, alpha_g=0.85)
MARGIN = 0.30
DESCENT_RELS_MIN = 63.97   # = descent16.P.LOG2_RELS
N_POOL = 40000   # pool draws on the 19/50 region (42e80f80 seeds 33/34)
YIELD_SE_K = 2.0   # descent relation input <= retained count x (1 - 2 relative...)
PAD_N = 1.10   # covers every SCALE_NU value from 2^15 up (STORAGE.md)
DESCENT_CAP = 115.3115   # meter: success >= 0.39 at RHO_CAP 0.05 (audit_g1.py)
NODE_CAP_LOG2 = 25.7    # explicit special-q node counter for both logs (proof.md...)
CANON_WORK_ALLOWANCE = 1.0312   # >= 1 + WORK_SE_K x the incidence-inclusive relative SE of the...
WORK_SE_K = 2.0
ECM_FIND = 0.999   # main-phase per-prime ECM success (this revision)
REC_ROUTE = 10.3        # per record: emit 7 + owner copy 13/4 + page hand-off <= 0.01...
CELL_ROUTE = 8.0        # instructions per cell (row bookkeeping, sentinels)
LINE_TAB = [(13.82, F(909617, 37746)), (9.82, F(1914335, 61528)), (5.82, F(846656, 9587)), (1.82, F(7448729, 23974))]
REC_TAB = [(18.91, F(23449665, 522256)), (10.91, F(3216127, 41121)), (3.91, F(1745730, 6539)), (0.91, F(3721383, 3967)), (-1.09, F(6474043, 2901))]
# exact logged quotients (audit re-parses measurements12.txt)
BOUND_SHIFT = 0.01
C_COSET_EMPTY = 1600.0
C_APPLY = 7.00
C_RECOVER = 7.003
CHUNK_CELLS = 2.0 ** 30   # chunk size (cells); recovery is not gated
C_SCAN_CELL, C_SCAN_LEAF = 8.0, 21.0
C_BUILD = 12.325
C_ZERO = 0.625
C_SETUP6 = 406.87 + 72.0 + 21.13
WIDE_LP = 63.0; WIDE_REC = 2.0; WIDE_SETUP = 2.0   # [2^63, 2^64) charged 2x (KERNELS.md)
C_NB_EVAL = 1508.2
C_NB_ROW = 4300.0
C_NB_SLICE = 17500.0
MONT5 = 367.76
C_NORMS = 2.0 ** 22
NBAR = 2.0
LX0, P_MID = 11.0, 48.0   # line sieve (scores only) 2^11..2^48
P_TDIV = 20.0          # survivors: trial division below 2^20
BAND_Q, BAND_ROUNDS = 0.999, 2   # band strip: two independent rounds of the 48-bit ecm16...
TDIV_PER_PRIME = 512.0  # per-prime trial-division bound per side (STORAGE.md)
def _tdiv_prime_count(log2p):   # Rosser-Schoenfeld pi(x) < 1.25506 x/ln x (x >= 17)
 return math.ceil(1.25506 * 2.0 ** log2p / (log2p * math.log(2))) + 1
N_TDIV_PRIMES = _tdiv_prime_count(P_TDIV)
C_TDIV = 2.0 * TDIV_PER_PRIME * N_TDIV_PRIMES
P_MIN = 2048
def ecm_instr(fbits, nbits, success=ECM_FIND):
 """Explicit ecm16 schedule ([12])"""
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
_TF = {}
def tab(t, log2h):   # Fraction op float is float(Fraction) op float
 log2h = log2h - BOUND_SHIFT
 if log2h >= t[0][0]:
  return t[0][1]
 seg = _TF.get(id(t)) or _TF.setdefault(id(t), [(h1, h2, float(c2), float(c1 - c2)) for (h1, c1), (h2, c2) in zip(t, t[1:])])
 for h1, h2, c2, d in seg:
  if log2h >= h2:
   return c2 + d * (log2h - h2) / (h1 - h2)
 return t[-1][1]
def load_pool(text):
 raw = lzma.decompress(base64.b85decode(''.join(text.split())), format=lzma.FORMAT_RAW,
      filters=[{"id": lzma.FILTER_LZMA2}])
 hdr, blob = raw.split(b'\n', 1)
 hdr = json.loads(hdr)
 fmt = '<BHHHHBBBBB'; m = hdr['rows']; cols = []; o = 0
 for ch in fmt[1:]:                     # column-major, index as differences mod 2^16 (hdr['layout'])
  k = struct.calcsize('<' + ch) * m; cols.append(struct.unpack('<' + ch * m, blob[o:o + k])); o += k
 idx = [cols[1][0]]
 for d in cols[1][1:]: idx.append((idx[-1] + d) & 0xFFFF)
 cols[1] = idx
 rows = []
 for seed, index, qa16, nf16, ng16, flags, evals, leaves, blen, rlen16 in zip(*cols):
  rows.append(dict(seed=seed, index=index, qa=qa16 / 65536.0, nf=nf16 / 16.0, ng=ng16 / 16.0,
      primitive=bool(flags & 1), rejected=bool(flags & 2), nonfinite=bool(flags & 4),
      evals=evals, leaves=leaves, blen=blen, rlen=16 * rlen16))
 assert len(rows) == hdr['rows']
 return hdr, rows
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
def evaluate(rows, n_total, geo, op=OP, retain=None, cost=None, scale_consts=1.0, verbose=False, return_work=False):
 """relation collection: chunk-routed records, every block charged."""
 lb, lx2, excess = op['lb'], op['lx2'], op['excess']
 Sf, Sg, W = op['slack_f'], op['slack_g'], op['width_bits']
 retain = retain or Retention(int(lb), int(lx2), prime_min=P_MIN, step=.25, max_bits=850, tilt=.05)
 cost = cost or Retention(850, int(lx2), prime_min=P_MIN, step=.25, max_bits=851, tilt=.005)
 wsum = 0.0; ysum = 0.0; y2 = 0.0; rej_w = 0.0; nf_m = ng_m = 0.0; nf_r = ng_r = 0.0
 work_sum = 0.0
 work_list = [] if return_work else None; wy_list = []
 c_ecm = ecm_instr(lb, Sf + 8) + ecm_instr(lb, Sg + 8)
 c_surv = C_NORMS + C_TDIV + c_ecm; c_surv_max = c_surv
 pairs_all = (NBAR / 2) * 2 ** log2_li(lx2)
 for r in rows:
  w = (1 - r['qa']) ** 3
  wsum += w
  nb = [64 * math.ceil((r[k] + 1) / 64) for k in ('nf', 'ng')]   # band strip sized by this cell's exact norms (whole 64-bit...)
  c_band = BAND_ROUNDS * sum(ecm_instr(P_MID, n, BAND_Q) for n in nb)
  band_keep = math.prod(1 - (n // int(P_TDIV)) * (1 - BAND_Q) ** BAND_ROUNDS for n in nb)   # a missed band prime loses the relation
  c_surv_max = max(c_surv_max, c_surv + c_band)
  if r['primitive'] and not r['rejected']:
   rf = retain.at(r['nf'], Sf)['retention']; rg = retain.at(r['ng'], Sg)['retention']
   uf_r, ug_r = (r['nf'] + op['alpha_f']) / lb, (r['ng'] + op['alpha_g']) / lb
   miss = ECM_FIND ** (lam_top(uf_r, lx2 / lb) + lam_top(ug_r, lx2 / lb))
   y = rho(uf_r) * rho(ug_r) * rf * rg * miss * band_keep
   ysum += w * y; y2 += w * w * y * y
   if return_work: wy_list.append((w * y, r['nf'], r['ng']))
   nf_r += w * y * r['nf']; ng_r += w * y * r['ng']
  if r['rejected']:
   rej_w += w
  nf_m += w * r['nf']; ng_m += w * r['ng']
  S = geo.slice_cells(r['qa']); K = geo.lines(r['qa'])
  p_all = dict(build=C_BUILD, zero=0.0, setup6=(pairs_all + (WIDE_SETUP - 1) * max(0.0, pairs_all - (NBAR / 2) * 2 ** log2_li(min(lx2, WIDE_LP)))) * C_SETUP6 / S,
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
    rec += NBAR * dsum * tab(REC_TAB, log2h) * (WIDE_REC if lp >= WIDE_LP else 1.0)
   else:
    dpairs = (NBAR / 2) * (2 ** log2_li(lp2) - 2 ** log2_li(lp))
    empty += dpairs * C_COSET_EMPTY / S * (WIDE_REC if lp >= WIDE_LP else 1.0)
   lp = lp2
  hits = NBAR * math.log(lx2 / P_MID); hits_wide = NBAR * math.log(lx2 / max(P_MID, min(lx2, WIDE_LP)))
  rec_hits = hits   # line hits are scores only
  sf = cost.survivor_upper(r['nf'], Sf, W); sg = cost.survivor_upper(r['ng'], Sg, W)
  p_keep = dict(line=line, rec=rec, empty=empty,
     apply=(rec_hits + (WIDE_REC - 1) * hits_wide) * C_APPLY,
     recover=rec_hits * C_RECOVER,
     route=CELL_ROUTE + REC_ROUTE * rec_hits,
     scan=C_SCAN_CELL + C_SCAN_LEAF * r['leaves'] / r['blen'],
     zero=C_ZERO, cof=sf * sg * (c_surv + c_band))
  per_cell = sum(p_all.values()) + sum(p_keep.values())
  per_cell = per_cell * scale_consts if scale_consts != 1.0 else per_cell
  contrib = S * per_cell
  work_sum += contrib
  if return_work: work_list.append(contrib)
 n_a_upper = geo.va / n_total
 population = n_a_upper * sum(geo.slice_cells(r['qa']) for r in rows)
 relcol_instr = n_a_upper * work_sum * CANON_WORK_ALLOWANCE
 per_cell_mean = relcol_instr / population
 rels = n_a_upper * geo.vb * geo.lower * ysum
 # fixed-pool estimator over all N_POOL original draws (zero...)
 # relative variance (Q/S^2 - 1/N) N/(N - 1), in bits
 yield_se = math.sqrt(max(0.0, y2 / ysum ** 2 - 1.0 / N_POOL) * N_POOL / (N_POOL - 1)) / LN2 if ysum > 0 else float('inf')
 nf_m /= wsum; ng_m /= wsum
 log2_U = math.log2((NBAR / 2) * 2 ** log2_li(lb) + 2 ** 31)
 log2_need = log2_U + math.log2(excess)
 log2_rels = math.log2(rels) if rels > 0 else float('-inf')
 feasible = log2_rels >= log2_need + MARGIN and log2_rels + math.log2(1 - YIELD_SE_K * yield_se * LN2) >= DESCENT_RELS_MIN
 log2_relcol = math.log2(relcol_instr / UNIT)
 out = dict(op={k: (str(v) if isinstance(v, F) else v) for k, v in op.items()},
    log2_population=math.log2(population), per_cell=per_cell_mean,
    log2_relcol=log2_relcol,
    log2_rels=log2_rels, log2_U=log2_U, log2_need=log2_need,
    margin=log2_rels - log2_need, yield_se_log2=yield_se, feasible=feasible,
    rejected_weight=rej_w / wsum, nf_mean=nf_m, ng_mean=ng_m, nf_rel=nf_r / ysum, ng_rel=ng_r / ysum,
    records_per_cell=NBAR * math.log(lx2 / P_MID),
    c_surv_log2=math.log2(c_surv_max), mean_slice_log2=math.log2(population / (n_a_upper * len(rows))))
 if return_work:
  out['work_contribs'] = work_list
  for k, j in (('nf_rel_se', 1), ('ng_rel_se', 2)):   # ratio-estimator SE of the relation-weighted mean norm over...
   R = out[k[:6]]; out[k] = math.sqrt(N_POOL / (N_POOL - 1) * sum((a * (t[j] - R)) ** 2 for t in wy_list for a in (t[0],))) / ysum
 return out
FILT_TABLE = {   # signed sigma-quotient system at scale 2^17, excess 1.3
 1.3: [(1, 31.713, 4133, 71, 300.06), (2, 32.850, 3990, 56, 300.06), (3, 31.015, 4226, 53, 300.02)],
}
W_ROW = 304.0   # mean row weight bound (max measured 302.4)
SCALE_NU = {   # N/U at mean weight 300
 12: (0.04175, 0.05567, 0.04346, 0.04493, 0.04542, 0.05567, 0.0647, 0.05103, 0.06519, 0.042),
 13: (0.03895, 0.0431, 0.03614, 0.03431, 0.03443, 0.04029),
 14: (0.03608, 0.03455, 0.03675, 0.03388, 0.03577, 0.03602),
 15: (0.03299, 0.0348, 0.03431), 16: (0.03274, 0.03162, 0.03335), 17: (0.03154, 0.03045, 0.03225), 18: (0.03038, 0.02855),
}
C_ENTRY = 20.7   # per entry 4 + 3 + 0.88*11 + 0.12*33 = 20.64 (STORAGE.md)
LA_OWNERS = 7    # active vector-block owners: N + 21 < 7 x 2^56 (STORAGE.md)
C_ROW = 1320.0   # per row and product: 7-owner row tail incl
C_RED9 = 420.0   # one 9-limb accumulator mod r plus a modular add
UNIT_FRAC_CHARGED = 0.88   # H4u, below every measured value (audit_g1.py unit_entries)
UNIT_FRAC = UNIT_FRAC_CHARGED
BW_M = 22
BW_N = 4
N_SM = 21                   # dense Schirokauer columns of A = [M | S]
C_FULL = 97.0
LOG2_SM = 88.0
LOG2_FBC = 72.0
LOG2_CANON = 90.0
LOG2_MISC = 80.0
LOG2_VERIFY = 11.0
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
 log2_N = log2_U + math.log2(ratio * pad_n)
 # Block Wiedemann on A = [M | S] (STORAGE.md)
 c_prod = W_ROW * C_ENTRY + C_ROW + N_SM * C_FULL
 per_n2 = c_prod * (BW_N / BW_M + 1) + (1 / BW_M + 1 / BW_N) * (c_prod + BW_N * C_FULL + C_RED9) + C_FULL * (BW_M + BW_N)
 per_n1 = (BW_M + BW_N) * C_RED9
 log2_la = math.log2((2 ** (2 * log2_N) * per_n2 + 2 ** log2_N * per_n1) / UNIT)
 return dict(ratio=ratio, log2_N=log2_N, log2_la=log2_la, log2_gen=generator(log2_N))
def generator(log2_N):
 # (m+n)x(m+n) generator of the length-L sequence by approximant...
 # 4(m+n)^3 L^log2(3) log2(L) multiplications at C_FULL + C_RED9
 k = BW_M + BW_N
 log2_L = math.log2((2 ** log2_N + N_SM) * (1 / BW_M + 1 / BW_N) + 2)
 return math.log2(4 * k ** 3 * (C_FULL + C_RED9) / UNIT) + math.log2(3) * log2_L + math.log2(log2_L)
def total_log2(log2_relcol, log2_la, log2_desc=DESCENT_CAP, log2_gen=0.0):
 return math.log2(2 ** log2_relcol + 2 ** log2_la + 2 ** log2_desc + 2 ** log2_gen + 2 ** LOG2_SM + 2 ** LOG2_FBC
     + 2 ** LOG2_CANON + 2 ** LOG2_MISC + 2 ** LOG2_VERIFY + 2 ** 50)
