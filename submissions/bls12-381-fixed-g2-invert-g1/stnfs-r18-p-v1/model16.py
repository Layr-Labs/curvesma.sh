"""stnfs ledger model: relation collection, linear algebra, generator, metered descent cap (log2 units)."""
import math, json, struct, zlib, lzma, base64
from fractions import Fraction as F
from common12 import rho, integer_metric, quadratic
import filtsim
from threshold_model import Retention
LN2 = math.log(2)
UNIT = 20000.0
OP = dict(lb=69.0, lx2=64.0, excess=1.2, slack_f=293, slack_g=224, width_bits=8, L0=64,
     theta=F(1, 2), alpha_f=5.0, alpha_g=0.844)
MARGIN = 0.30
PAD_N = 1.25   # 1.25x the 2^19 worst-seed N/U; covers all sims at scales >= 2^15 (R11)
DESCENT_CAP = 118.76810   # meter: E 118.03418 + gap 0.73393 (R18)
NODE_CAP_LOG2 = 27.3    # explicit special-q node counter for both logs (proof.md section 7)
CANON_WORK_ALLOWANCE = 1.0772   # >= 1 + WORK_SE_K x the incidence-inclusive relative SE of the sampled work estimator ([4a]; R18)
WORK_SE_K = 2.0
ECM_FIND = 0.999   # main-phase per-prime ECM success (this revision; was 0.99): the same ecm16 grid, more curves per survivor
BAND_FBITS = 48.0   # R16 band-strip: recover (2^20, 2^48] primes credited by the R15 line path
BAND_NBITS = (832, 576)   # unstripped cofactor caps: pool maxima 796.375/513.0 bits, up to 13/9 limbs
BAND_FMAX = (41, 28)   # max band factors per side: floor(nbits/20)
BAND_Q = 0.999   # per-round per-prime band-schedule find probability
BAND_ROUNDS = 2   # fixed rounds per side; failure <= FMAX*(1-Q)^ROUNDS
REC_ROUTE = 10.005      # measured 10.0015 (STORAGE.md), 0.035% margin (R7)
CELL_ROUTE = 8.0        # instructions per cell (row bookkeeping, sentinels)
LINE_TAB = [(13.82, 24.09836803899751), (9.82, 31.113232999609934), (5.82, 88.31292375091269), (1.82, 310.7003003253525)]  # exact raw-count quotients (R5)
REC_TAB = [(18.91, 44.90070961367605), (10.91, 78.21130322706159), (3.91, 266.9720140694296), (0.91, 938.0849508444669), (-1.09, 2231.659083074802)]  # exact raw-count quotients (R5)
C_COSET_EMPTY = 1600.0
C_APPLY = 7.00
C_RECOVER = 7.003
C_SCAN_CELL, C_SCAN_LEAF = 8.0, 21.0
C_BUILD = 12.325
C_ZERO = 0.625
C_SETUP6 = 406.87 + 72.0 + 21.13
C_NB_EVAL = 1508.1812874481486  # exact raw-count quotient (R5)
C_NB_ROW = 4300.0
C_NB_SLICE = 17500.0
MONT5 = 367.76
C_NORMS = 2.0 ** 22
C_TDIV = 2.0 ** 24.5
NBAR = 2.0
LX0, P_MID = 8.0, 48.0   # R15: line sieve to 2^48 (production hits/coset >= 2^6.79 stays inside the measured line-tab range); record walk above
P_MIN = 256
def ecm_instr(fbits, nbits, success=ECM_FIND):
  """Cheapest ecm16 (B1,B2) grid point finding every prime <= 2^fbits in an nbits-bit integer with probability >= success ([12])."""
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
  retain = retain or Retention(int(lb), int(lx2), step=.25, max_bits=850, tilt=.05)
  cost = cost or Retention(850, int(lx2), step=.25, max_bits=851, tilt=.005)
  wsum = 0.0; ysum = 0.0; y2 = 0.0; rej_w = 0.0; nf_m = ng_m = 0.0; nf_r = ng_r = 0.0
  work_sum = 0.0
  comp_sum = {}
  work_list = [] if return_work else None
  c_ecm = ecm_instr(lb, Sf + 8) + ecm_instr(lb, Sg + 8)
  c_band = BAND_ROUNDS * (ecm_instr(BAND_FBITS, BAND_NBITS[0], BAND_Q) + ecm_instr(BAND_FBITS, BAND_NBITS[1], BAND_Q))
  c_surv = C_NORMS + C_TDIV + c_ecm + c_band
  band_keep = (1 - BAND_FMAX[0] * (1 - BAND_Q) ** BAND_ROUNDS) * (1 - BAND_FMAX[1] * (1 - BAND_Q) ** BAND_ROUNDS)
  pairs_all = (NBAR / 2) * 2 ** log2_li(lx2)
  for r in rows:
    w = (1 - r['qa']) ** 3
    wsum += w
    if r['primitive'] and not r['rejected']:
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
        rec += NBAR * dsum * tab(REC_TAB, log2h)
      else:
        dpairs = (NBAR / 2) * (2 ** log2_li(lp2) - 2 ** log2_li(lp))
        empty += dpairs * C_COSET_EMPTY / S
      lp = lp2
    hits = NBAR * math.log(lx2 / P_MID)
    sf = cost.survivor_upper(r['nf'], Sf, W); sg = cost.survivor_upper(r['ng'], Sg, W)
    p_keep = dict(line=line, rec=rec, empty=empty,
           apply=hits * C_APPLY, recover=hits * C_RECOVER,
           route=CELL_ROUTE + REC_ROUTE * hits,
           scan=C_SCAN_CELL + C_SCAN_LEAF * r['leaves'] / r['blen'],
           zero=C_ZERO, cof=sf * sg * c_surv)
    per_cell = sum(p_all.values()) + sum(p_keep.values())
    per_cell = per_cell * scale_consts if scale_consts != 1.0 else per_cell
    contrib = S * per_cell
    work_sum += contrib
    for _k, _v in list(p_all.items()) + list(p_keep.items()):
      comp_sum[_k] = comp_sum.get(_k, 0.0) + S * _v
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
  feasible = log2_rels >= log2_need + MARGIN
  log2_relcol = math.log2(relcol_instr / UNIT)
  out = dict(op={k: (str(v) if isinstance(v, F) else v) for k, v in op.items()},
        comp_log2={k: (math.log2(n_a_upper * v / UNIT) if v > 0 else float('-inf')) for k, v in comp_sum.items()},
        log2_population=math.log2(population), per_cell=per_cell_mean,
        log2_relcol=log2_relcol,
        log2_rels=log2_rels, log2_U=log2_U, log2_need=log2_need,
        margin=log2_rels - log2_need, yield_se_log2=yield_se, feasible=feasible,
        rejected_weight=rej_w / wsum, nf_mean=nf_m, ng_mean=ng_m, nf_rel=nf_r / ysum, ng_rel=ng_r / ysum,
        records_per_cell=NBAR * math.log(lx2 / P_MID),
        c_surv_log2=math.log2(c_surv), c_band_log2=math.log2(c_band), band_keep=band_keep,
        mean_slice_log2=math.log2(population / (n_a_upper * len(rows))))
  if return_work: out['work_contribs'] = work_list
  return out
FILT_TABLE = {   # 2^19 filtering basis at excess 1.2: (seed, U/N, N, max coef, weight); U/N rounded down
  1.2: [(1, 18.206, 28796, 43, 100.001), (2, 18.638, 28129, 35, 100.008), (3, 18.831, 27841, 37, 100.004), (4, 18.208, 28794, 37, 100.004), (5, 18.150, 28886, 36, 100.001), (6, 18.709, 28022, 45, 100.004), (7, 18.387, 28513, 47, 100.004), (8, 17.986, 29149, 41, 100.016)],
}
W_ROW = 100.7   # covers all filtering sims at scales >= 2^13 (max 100.61; R6)
C_ENTRY = 73.0   # 33+4+3+20+13 (STORAGE.md, R7)
BW_M = 22
BW_N = 4   # min over integer block sizes (REVISION.md R7)
N_SM = 21                   # dense Schirokauer columns of A = [M | S]: sigma-odd unit ranks 18 + 3 ([13])
C_FULL = 97.0
C_RED = 2 * MONT5 + 40      # mod-r reduction of a row accumulator: two five-limb Montgomery products + fixups
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
  # Block Wiedemann on square A = [M | S] with random m x n blocks: Krylov products, full-field projections, one Horner reconstruction
  c_prod = W_ROW * C_ENTRY + N_SM * C_FULL
  per_n2 = c_prod * (BW_N / BW_M + 1) + (1 / BW_M + 1 / BW_N) * (c_prod + BW_N * C_FULL) + C_FULL * (BW_M + BW_N)
  per_n2 += (BW_M + BW_N) * C_RED / 2 ** log2_N   # one reduction per projection inner product of length N
  log2_la = 2 * log2_N + math.log2(per_n2 / UNIT)
  return dict(ratio=ratio, log2_N=log2_N, log2_la=log2_la, log2_gen=generator(log2_N))
def generator(log2_N):
  # (m+n) x (m+n) matrix generator via divide-and-conquer approximant bases (Karatsuba); log2 units
  k = BW_M + BW_N
  log2_L = math.log2((2 ** log2_N + N_SM) * (1 / BW_M + 1 / BW_N) + 2)
  return math.log2(4 * k ** 3 * (C_FULL + C_RED) / UNIT) + math.log2(3) * log2_L + math.log2(log2_L)
def total_log2(log2_relcol, log2_la, log2_desc=DESCENT_CAP, log2_gen=0.0):
  return math.log2(2 ** log2_relcol + 2 ** log2_la + 2 ** log2_desc + 2 ** log2_gen + 2 ** LOG2_SM + 2 ** LOG2_FBC
           + 2 ** LOG2_CANON + 2 ** LOG2_MISC + 2 ** LOG2_VERIFY + 2 ** 50)
