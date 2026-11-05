"""stnfs ledger model: collection, linear algebra (STORAGE.md) and the metered descent cap (descent16)."""
import math, json, struct, zlib, lzma, base64
from fractions import Fraction as F
from common12 import rho, integer_metric, quadratic
import filtsim
from threshold_model import Retention
LN2 = math.log(2)
UNIT = 20000.0
OP = dict(lb=69.0, lx2=64.0, excess=1.3, slack_f=297, slack_g=233, width_bits=8, L0=64,
     theta=F(1, 2), alpha_f=5.0, alpha_g=0.844)
MARGIN = 0.30
DESCENT_RELS_MIN = 63.97   # = descent16.P.LOG2_RELS
N_POOL = 20000   # original pool draws (pool header samples; audit checks)
POOL_QA_CUT = F(191, 500)   # pool files keep canonical rows with stored qa <= this (SB1)
YIELD_SE_K = 1.0   # descent relation input <= retained count x (1 - 1 relative pool SE): linear lower confidence bound
PAD_N = 1.08   # covers every SCALE_NU value from 2^15 up (max 0.03480 <= 0.03225 x 1.08; STORAGE.md, linear algebra)
DESCENT_CAP = 118.3128
NODE_CAP_LOG2 = 27.3    # explicit special-q node counter for both logs (proof.md section 7)
CANON_WORK_ALLOWANCE = 1.04432   # >= 1 + WORK_SE_K x the fresh-incidence relative SE of the ratio work estimator ([4a]); verified from canon12.b85
FRESH_SEED = 4   # fresh canonicity sample (canon12.b85), header pinned by the audit
FRESH_N = 100000
FRESH_THETA = "24700/65536"
WORK_SE_K = 2.0
ECM_FIND = 0.999   # main-phase per-prime ECM success (R1)
Q_F_SAFETY = 100.0   # early-abort g-side probability charged at 100x the H1 model estimate (LEAP1; proof.md section 3)
REC_ROUTE = 10.3        # per record: scatter 7 + owner copy 13/4 + hand-off <= 0.01 (STORAGE.md, R4)
CELL_ROUTE = 8.0        # instructions per cell (row bookkeeping, sentinels)
LINE_TAB = [(13.82, F(909617, 37746)), (9.82, F(1914335, 61528)), (5.82, F(846656, 9587)), (1.82, F(7448729, 23974))]
REC_TAB = [(18.91, F(23449665, 522256)), (10.91, F(3216127, 41121)), (3.91, F(1745730, 6539)), (0.91, F(3721383, 3967)), (-1.09, F(6474043, 2901))]
# R11: exact quotients of measurements12.txt (audit re-parses); BOUND_SHIFT = 2x the two-decimal boundary rounding (H4, convex)
BOUND_SHIFT = 0.01
C_COSET_EMPTY = 1600.0
C_APPLY = 7.00
C_RECOVER = 7.003
CHUNK_CELLS = 2.0 ** 30   # chunk size (cells); recovery is not gated (R5 withdrawn)
PAGE_RECS = 2.0 ** 16; C_PAGE_FREE = 8.0; C_CHUNK_GATE = 16.0
C_SCAN_CELL, C_SCAN_LEAF = 8.0, 21.0
C_BUILD = 12.325
C_ZERO = 0.625
C_SETUP6 = 406.87 + 72.0 + 21.13
WIDE_LP = 63.0; WIDE_REC = 1.01; WIDE_SETUP = 1.0   # WIDE1 (KERNELS.md, wide12.txt): primes in [2^63, 2^64) at measured prices; 1.01 on record walk, empty cosets, apply
LOG2_WIDE_LAT = 70.6   # WIDE1 factor-base lattices: <= 6 pi(2^64) (p, s) at 2 x 5,615,924 Ir, separate ledger row
C_NB_EVAL = 1508.2
C_NB_ROW = 4300.0
C_NB_SLICE = 17500.0
MONT5 = 367.76
C_NORMS = 2.0 ** 22
NBAR = 2.0
LX0, P_MID = 11.0, 47.75   # line sieve (scores only) 2^11..2^47.75; record walk above
P_TDIV = 20.0          # survivors: trial division below 2^20, band strip for (2^20, 2^47.75], then ECM to 2^69
BAND_Q, BAND_ROUNDS = 0.999, 2   # two independent 47-bit ecm16 rounds at 99.9% per prime
TDIV_PER_PRIME = 512.0  # per-prime trial-division bound per side (STORAGE.md)
def _tdiv_prime_count(log2p):   # Rosser-Schoenfeld pi(x) < 1.25506 x/ln x (x >= 17); +1 covers float rounding
  return math.ceil(1.25506 * 2.0 ** log2p / (log2p * math.log(2))) + 1
N_TDIV_PRIMES = _tdiv_prime_count(P_TDIV)
C_TDIV = 2.0 * TDIV_PER_PRIME * N_TDIV_PRIMES
P_MIN = 2048
def ecm_instr(fbits, nbits, success=ECM_FIND):
  """Explicit ecm16 schedule ([12]): cheapest grid point finding every prime
  <= 2^fbits with probability >= success, at the measured MONT5."""
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
  log2h = log2h - BOUND_SHIFT   # R11: charge at a slightly smaller h (cost decreases in h)
  if log2h >= t[0][0]:
    return t[0][1]
  for (h1, c1), (h2, c2) in zip(t, t[1:]):
    if log2h >= h2:
      return c2 + (c1 - c2) * (log2h - h2) / (h1 - h2)
  return t[-1][1]
def load_pool(text):
  raw = lzma.decompress(base64.b85decode(''.join(text.split())), format=lzma.FORMAT_RAW,
             filters=[{"id": lzma.FILTER_LZMA2}])
  hdr, blob = raw.split(b'\n', 1)
  hdr = json.loads(hdr)
  fmt = '<BHHHHBBBBBB'; m = hdr['rows']; cols = []; o = 0
  for ch in fmt[1:]:                     # column-major, index as differences mod 2^16 (hdr['layout'])
    k = struct.calcsize('<' + ch) * m; cols.append(struct.unpack('<' + ch * m, blob[o:o + k])); o += k
  idx = [cols[1][0]]
  for d in cols[1][1:]: idx.append((idx[-1] + d) & 0xFFFF)
  cols[1] = idx
  rows = []
  for seed, index, qa16, nf16, ng16, flags, evals, leaves, blen, rlen16, wpack in zip(*cols):
    # WIDTH1 packed widths: wpack = 9*ceil(f/16)+ceil(g/16)
    rows.append(dict(seed=seed, index=index, qa=qa16 / 65536.0, nf=nf16 / 16.0, ng=ng16 / 16.0,
             primitive=bool(flags & 1), rejected=bool(flags & 2), nonfinite=bool(flags & 4),
             evals=evals, leaves=leaves, blen=blen, rlen=16 * rlen16, wf=16 * (wpack // 9), wg=16 * (wpack % 9)))
  assert len(rows) == hdr['rows']
  return hdr, rows
def load_fresh(text):
  """Fresh canonicity sample (canon12.b85): header JSON + first canonical draw
  index <I then successive differences <H, lzma-RAW + base85 like the pool."""
  raw = lzma.decompress(base64.b85decode(''.join(text.split())), format=lzma.FORMAT_RAW,
           filters=[{"id": lzma.FILTER_LZMA2}])
  hdr, blob = raw.split(b'\n', 1)
  hdr = json.loads(hdr)
  n = (len(blob) - 4) // 2
  first = struct.unpack('<I', blob[:4])[0]
  deltas = struct.unpack('<' + 'H' * n, blob[4:4 + 2 * n])
  assert 4 + 2 * n == len(blob)
  idx = [first]
  for d in deltas: idx.append(idx[-1] + d)
  assert len(idx) == hdr['C'] and all(b > a for a, b in zip(idx, idx[1:])) and idx[-1] < hdr['N']
  return hdr, idx
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
def row_count_factors(cfg, a_radius, rows):
  """Per-row b-side certs at guarded brad=isqrt((65535-qa16)*32)
  (+ a-side at a_radius); same half_count_ratios theorem. -1 ulp guards
  every rounding mode (qa16 within 1/2 ulp); isqrt floors."""
  import lattice_count_bounds as L
  A2 = cfg['production_A_squared']; t = L.correction(cfg, a_radius, 1024)['cube_radius']
  a = L.half_count_ratios(cfg['matrix_a'], a_radius, t[0]); B0 = []; B1 = []
  for r in rows:
    b = L.half_count_ratios(cfg['matrix_b'], math.isqrt((65535 - round(r['qa'] * 65536)) * 32), t[1])
    B0.append(float(b[0])); B1.append(float(b[1]))
  assert a[0] >= F(98, 100) and min(B0) >= 0.98 and a[1] <= F(102, 100) and max(B1) <= 1.02
  return (float(a[0]), float(a[1]), B0, B1)
def evaluate(rows, n_total, geo, op=OP, retain=None, cost=None, scale_consts=1.0, verbose=False, return_work=False, p_fresh=None, row_ab=None):
  """relation collection: chunk-routed records, every block charged.
  p_fresh (fresh-sample canonical rate) selects the ratio work estimator
  va x p_fresh x mean-row-work = expansion point x p_fresh/(rows/n_total);
  the yield estimator stays the frozen pool expansion in both cases.
  row_ab from row_count_factors; None uses geo.lower/geo.upper."""
  lb, lx2, excess = op['lb'], op['lx2'], op['excess']
  Sf, Sg = op['slack_f'], op['slack_g']   # WIDTH1: per-row width cover (r['wf']/wg)
  retain = retain or Retention(int(lb), int(lx2), prime_min=P_MIN, step=.25, max_bits=850, tilt=.05)
  cost = cost or Retention(850, int(lx2), prime_min=P_MIN, step=.25, max_bits=851, tilt=.005)
  if row_ab is None:
    a0, a1 = geo.lower, geo.upper
    Brow0 = Brow1 = [1.0] * len(rows)
  else:
    a0, a1, Brow0, Brow1 = row_ab
  wsum = 0.0; ysum = 0.0; y2 = 0.0; rej_w = 0.0; nf_m = ng_m = 0.0; nf_r = ng_r = 0.0
  work_sum = 0.0
  work_list = [] if return_work else None
  c_ecm_f, c_ecm_g = ecm_instr(lb, Sf + 8), ecm_instr(lb, Sg + 8)
  c_surv_max = 0.0; qf_max = 0.0; qf_sum = 0.0
  pairs_all = (NBAR / 2) * 2 ** log2_li(lx2)
  for idx, r in enumerate(rows):
    w = (1 - r['qa']) ** 3
    wsum += w
    nb = [64 * math.ceil((r[k] + 1) / 64) for k in ('nf', 'ng')]   # band sized by exact norms (whole limbs)
    c_band_f = BAND_ROUNDS * ecm_instr(P_MID, nb[0], BAND_Q)
    c_band_g = BAND_ROUNDS * ecm_instr(P_MID, nb[1], BAND_Q)
    band_keep = math.prod(1 - (n // int(P_TDIV)) * (1 - BAND_Q) ** BAND_ROUNDS for n in nb)   # a missed band prime loses the relation
    # LEAP1 early abort: g-side band+ECM runs only if f fully factors; charge q_f = min(1, 100 rho_f / P(omitted_f <= slack_f))
    q_f = min(1.0, Q_F_SAFETY * rho((r['nf'] + op['alpha_f']) / lb) / cost.survivor_upper(r['nf'], Sf, 0))
    qf_max = max(qf_max, q_f); qf_sum += q_f
    c_surv_row = C_NORMS + C_TDIV + c_band_f + c_ecm_f + q_f * (c_band_g + c_ecm_g)
    c_surv_max = max(c_surv_max, c_surv_row)
    if r['primitive'] and not r['rejected']:
      rf = retain.at(r['nf'], Sf)['retention']; rg = retain.at(r['ng'], Sg)['retention']
      uf_r, ug_r = (r['nf'] + op['alpha_f']) / lb, (r['ng'] + op['alpha_g']) / lb
      miss = ECM_FIND ** (lam_top(uf_r, lx2 / lb) + lam_top(ug_r, lx2 / lb))
      y = rho(uf_r) * rho(ug_r) * rf * rg * miss * band_keep
      z = w * y * Brow0[idx]
      ysum += z; y2 += z * z
      nf_r += z * r['nf']; ng_r += z * r['ng']
    if r['rejected']:
      rej_w += w
    nf_m += w * r['nf']; ng_m += w * r['ng']
    S = geo.vb * (1 - r['qa']) ** 3 * a1 * Brow1[idx]; K = geo.lines(r['qa'])
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
    sf = cost.survivor_upper(r['nf'], Sf, r['wf'] / 16.0); sg = cost.survivor_upper(r['ng'], Sg, r['wg'] / 16.0)
    p_keep = dict(line=line, rec=rec, empty=empty,
           apply=(rec_hits + (WIDE_REC - 1) * hits_wide) * C_APPLY,
           recover=rec_hits * C_RECOVER,
           route=CELL_ROUTE + REC_ROUTE * rec_hits,
           scan=C_SCAN_CELL + C_SCAN_LEAF * r['leaves'] / r['blen'],
           zero=C_ZERO, cof=sf * sg * c_surv_row)
    per_cell = sum(p_all.values()) + sum(p_keep.values())
    per_cell = per_cell * scale_consts if scale_consts != 1.0 else per_cell
    contrib = S * per_cell
    work_sum += contrib
    if return_work: work_list.append(contrib)
  n_a_upper = geo.va / n_total
  population = n_a_upper * sum(geo.vb * (1 - r['qa']) ** 3 * a1 * Brow1[i] for i, r in enumerate(rows))
  ratio_factor = p_fresh * n_total / len(rows) if p_fresh is not None else 1.0
  relcol_instr = n_a_upper * work_sum * CANON_WORK_ALLOWANCE * ratio_factor
  per_cell_mean = relcol_instr / population
  rels = n_a_upper * geo.vb * a0 * ysum
  # fixed-pool estimator over all N_POOL original draws (zero outside the canonical primitive sub-region):
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
        log2_relcol=log2_relcol, work_ratio_factor=ratio_factor,
        log2_rels=log2_rels, log2_U=log2_U, log2_need=log2_need,
        margin=log2_rels - log2_need, yield_se_log2=yield_se, feasible=feasible,
        rejected_weight=rej_w / wsum, nf_mean=nf_m, ng_mean=ng_m, nf_rel=nf_r / ysum, ng_rel=ng_r / ysum,
        qf_max=qf_max, qf_mean=qf_sum / len(rows),
        records_per_cell=NBAR * math.log(lx2 / P_MID),
        c_surv_log2=math.log2(c_surv_max), mean_slice_log2=math.log2(population / (n_a_upper * len(rows))))
  if return_work: out['work_contribs'] = work_list
  return out
FILT_TABLE = {   # signed sigma-quotient system at scale 2^17, excess 1.3: (seed, U/N, N, max coef, weight); R17 W=300
  1.3: [(1, 31.713, 4133, 71, 300.06), (2, 32.850, 3990, 56, 300.06), (3, 31.015, 4226, 53, 300.02)],
}
W_ROW = 304.0   # mean row weight bound (R17: max measured 302.4)
SCALE_NU = {   # N/U at mean weight 300, excess 1.3, by log2 U (filtsim.run, kmax 64); R17
  12: (0.04175, 0.05567, 0.04346, 0.04493, 0.04542, 0.05567, 0.0647, 0.05103, 0.06519, 0.042),
  13: (0.03895, 0.0431, 0.03614, 0.03431, 0.03443, 0.04029),
  14: (0.03608, 0.03455, 0.03675, 0.03388, 0.03577, 0.03602),
  15: (0.03299, 0.0348, 0.03431), 16: (0.03274, 0.03162, 0.03335), 17: (0.03154, 0.03045, 0.03225), 18: (0.03038, 0.02855),
}
C_ENTRY = 20.7   # per entry 4 + 3 + 0.88*11 + 0.12*33 = 20.64 (STORAGE.md)
LA_OWNERS = 7    # active vector-block owners: N + 21 < 7 x 2^56 (STORAGE.md)
C_ROW = 1320.0   # per row and product: 7-owner row tail incl. reduction (STORAGE.md), counted 1317 by [10]
C_RED9 = 420.0   # one 9-limb accumulator mod r plus a modular add, counted 385 by [10]
UNIT_FRAC_CHARGED = 0.88   # H4u, below every measured value (audit_g1.py unit_entries)
UNIT_FRAC = UNIT_FRAC_CHARGED
BW_M = 22
BW_N = 4
N_SM = 21                   # dense Schirokauer columns of A = [M | S]: sigma-odd unit ranks 18 + 3 ([13])
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
  # Block Wiedemann on A = [M | S] (STORAGE.md): Krylov + Horner + projections.
  c_prod = W_ROW * C_ENTRY + C_ROW + N_SM * C_FULL
  per_n2 = c_prod * (BW_N / BW_M + 1) + (1 / BW_M + 1 / BW_N) * (c_prod + BW_N * C_FULL + C_RED9) + C_FULL * (BW_M + BW_N)
  per_n1 = (BW_M + BW_N) * C_RED9
  log2_la = math.log2((2 ** (2 * log2_N) * per_n2 + 2 ** log2_N * per_n1) / UNIT)
  return dict(ratio=ratio, log2_N=log2_N, log2_la=log2_la, log2_gen=generator(log2_N))
def generator(log2_N):
  # (m+n)x(m+n) generator of the length-L sequence by approximant bases: <=
  # 4(m+n)^3 L^log2(3) log2(L) multiplications at C_FULL + C_RED9; log2 units.
  k = BW_M + BW_N
  log2_L = math.log2((2 ** log2_N + N_SM) * (1 / BW_M + 1 / BW_N) + 2)
  return math.log2(4 * k ** 3 * (C_FULL + C_RED9) / UNIT) + math.log2(3) * log2_L + math.log2(log2_L)
def total_log2(log2_relcol, log2_la, log2_desc=DESCENT_CAP, log2_gen=0.0):
  return math.log2(2 ** log2_relcol + 2 ** log2_la + 2 ** log2_desc + 2 ** log2_gen + 2 ** LOG2_SM + 2 ** LOG2_FBC + 2 ** LOG2_WIDE_LAT
           + 2 ** LOG2_CANON + 2 ** LOG2_MISC + 2 ** LOG2_VERIFY + 2 ** 50)
