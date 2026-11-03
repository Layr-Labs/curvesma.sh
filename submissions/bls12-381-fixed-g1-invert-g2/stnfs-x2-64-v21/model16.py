"""stnfs ledger model: relation collection (chunk-routed records, every block"""
import math, json, struct, zlib, lzma, base64
from fractions import Fraction as F
from common12 import rho, integer_metric, quadratic, regenerate_pair, norms, mul, norm, integer_interval
import ideal_canonical as IC
import embedding_bounds as EB
import filtsim
from threshold_model import Retention
LN2 = math.log(2)
UNIT = 20000.0
OP = dict(lb=69.0, lx2=64.0, excess=1.1, slack_f=295, slack_g=226, width_bits=8, L0=64,
     theta=F(5, 16), alpha_f=5.15, alpha_g=0.85)
MARGIN = 0.11   # legacy; feasibility uses MARGIN_FLOOR/MARGIN_SE below
MARGIN_FLOOR = 0.06   # absolute margin floor (bits over required excess)
MARGIN_SE = 2.0       # margin >= max(floor, MARGIN_SE x SE)
Q_CORE = 0.02         # core cut: qa < Q_CORE rows give no yield (cells charged);
                      # dropping a non-negative stratum only lowers the estimate (proof 3)
PAD_N = 1.1192   # covers SCALE_NU 2^15+ at excess 1.1 (proof 4)
DESCENT_CAP = 113.214
NODE_CAP_LOG2 = 27.3    # node counter cap (proof 7)
CANON_WORK_ALLOWANCE = 1.0165   # >= 1 + WORK_SE_K x work-estimator relative SE ([4a])
WORK_SE_K = 2.0
Q_F_SAFETY = 100.0   # early-abort prob charged at 100x the model estimate (proof 3)
ECM_FIND = 0.999
REC_ROUTE = 10.01        # measured record route (STORAGE.md)
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
LX0, P_MID = 8.0, 48.0   # line sieve to 2^48; record walk above
BAND_FBITS = 48.0   # band strip recovers line-credited (2^20, 2^48] primes
BAND_NBITS = (832, 576)   # unstripped cofactor caps
BAND_FMAX = (41, 28)   # max band factors/side: floor(nbits/20)
BAND_Q = 0.999   # per-round per-prime find prob
BAND_ROUNDS = 2
P_MIN = 256
ECM_TAB = {20: (11e3, 86), 25: (5e4, 214), 30: (25e4, 430), 35: (1e6, 904), 40: (3e6, 2350), 45: (11e6, 4480),
      50: (43e6, 7553), 55: (11e7, 17769), 60: (26e7, 42017), 65: (85e7, 69408), 70: (29e8, 102212),
      75: (76e8, 188056), 80: (25e9, 265557)}
def ecm_instr(fbits, nbits, success=ECM_FIND):
  """The explicit ecm16 schedule (the algorithm actually run, [12]): the cheapest"""
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
def load_pool(text, cfg):
  """Compact pool: per-seed canonical-index delta-varint lists plus a field-anchor"""
  raw = lzma.decompress(base64.b85decode(''.join(text.split())), format=lzma.FORMAT_RAW,
             filters=[{"id": lzma.FILTER_LZMA2}])
  hdr, blob = raw.split(b'\n', 1)
  hdr = json.loads(hdr)
  def varints(buf):
    out = []; n = 0; s = 0
    for byte in buf:
      n |= (byte & 0x7F) << s
      if byte & 0x80:
        s += 7
      else:
        out.append(n); n = 0; s = 0
    assert s == 0 and n == 0
    return out
  o = 0
  lists = {}
  for sd, key in ((1, 'seed1_bytes'), (2, 'seed2_bytes')):
    idx = []; prev = -1
    for d in varints(blob[o:o + hdr[key]]):
      prev += d; idx.append(prev)
    lists[sd] = idx; o += hdr[key]
  fmt = '<HHHHBBBBB'; sz = struct.calcsize(fmt)
  ab = blob[o:o + hdr['anchor_bytes']]
  assert len(ab) == hdr['anchor_rows'] * sz and o + hdr['anchor_bytes'] == len(blob)
  anchors = []
  for k in range(hdr['anchor_rows']):
    a = struct.unpack(fmt, ab[k * sz:(k + 1) * sz])
    anchors.append(dict(position=a[0], qa16=a[1], nf16=a[2], ng16=a[3], flags=a[4],
                 evals=a[5], leaves=a[6], blen=a[7], rlen16=a[8]))
  hdr['anchors'] = anchors
  qa, da, ia = integer_metric(cfg['matrix_a']); qb, db, ib = integer_metric(cfg['matrix_b'])
  A2 = cfg['production_A_squared']
  disks = EB.root_disks()
  rows = []
  for sd in (1, 2):
    for i in lists[sd]:
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
      leaves = est.adaptive_row(a, b[1:], lo, hi, singleton_policy='reject', width=128)
      at = [l for l in leaves if l[0] <= b[0] <= l[1]][0]
      rej = isinstance(at[2], dict)
      rows.append(dict(seed=sd, index=i, qa=qa16 / 65536.0, nf=nf16 / 16.0, ng=ng16 / 16.0,
               primitive=prim, rejected=rej, nonfinite=False,
               evals=est.counts['block_evaluations'], leaves=len(leaves), blen=blen, rlen=rlen))
  assert len(rows) == hdr['canonical_rows']
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
  work_list = [] if return_work else None
  # early abort: f side first; g-side ECM only if f factored, prob q_f = 100 rho_f / P(omitted_f <= slack_f)
  c_ecm_f, c_ecm_g = ecm_instr(lb, Sf + 8), ecm_instr(lb, Sg + 8)
  # band strip: two 48-bit rounds per side (f before its 69-bit ECM, g after f factored); a miss loses the relation
  c_band_f = BAND_ROUNDS * ecm_instr(BAND_FBITS, BAND_NBITS[0], BAND_Q)
  c_band_g = BAND_ROUNDS * ecm_instr(BAND_FBITS, BAND_NBITS[1], BAND_Q)
  band_keep = (1 - BAND_FMAX[0] * (1 - BAND_Q) ** BAND_ROUNDS) * (1 - BAND_FMAX[1] * (1 - BAND_Q) ** BAND_ROUNDS)
  pairs_all = (NBAR / 2) * 2 ** log2_li(lx2) + ((WIDE_FACTOR - 1) * (NBAR / 2) * (2 ** log2_li(lx2) - 2 ** log2_li(WIDE_FROM)) if lx2 > WIDE_FROM else 0.0)
  n_cut = 0
  for r in rows:
    w = (1 - r['qa']) ** 3
    wsum += w
    if r['qa'] < Q_CORE:
      n_cut += 1   # core cut: no yield credit, work still charged
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
    q_f = min(1.0, Q_F_SAFETY * rho((r['nf'] + op['alpha_f']) / lb) / cost.survivor_upper(r['nf'], Sf, 0))
    c_surv = C_NORMS + C_TDIV + c_band_f + c_ecm_f + q_f * (c_band_g + c_ecm_g)
    p_keep = dict(line=line, rec=rec, empty=empty,
           apply=hits * C_APPLY, recover=hits * C_RECOVER,
           route=CELL_ROUTE + REC_ROUTE * hits,
           scan=C_SCAN_CELL + C_SCAN_LEAF * r['leaves'] / r['blen'],
           zero=C_ZERO, cof=sf * sg * c_surv)
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
  yield_se = math.sqrt(max(0.0, y2 / n_total - (ysum / n_total) ** 2) / n_total) / (ysum / n_total) / LN2 if ysum > 0 else float('inf')
  nf_m /= wsum; ng_m /= wsum
  log2_U = math.log2((NBAR / 2) * 2 ** log2_li(lb) + 2 ** 31)
  log2_need = log2_U + math.log2(excess)
  log2_rels = math.log2(rels) if rels > 0 else float('-inf')
  margin_need = max(MARGIN_FLOOR, MARGIN_SE * yield_se)
  feasible = log2_rels >= log2_need + margin_need
  log2_relcol = math.log2(relcol_instr / UNIT)
  out = dict(op={k: (str(v) if isinstance(v, F) else v) for k, v in op.items()},
        log2_population=math.log2(population), per_cell=per_cell_mean,
        log2_relcol=log2_relcol,
        log2_rels=log2_rels, log2_U=log2_U, log2_need=log2_need,
        margin=log2_rels - log2_need, yield_se_log2=yield_se, feasible=feasible, margin_need=margin_need,
        n_core_cut=n_cut,
        rejected_weight=rej_w / wsum, nf_mean=nf_m, ng_mean=ng_m, nf_rel=nf_r / ysum, ng_rel=ng_r / ysum,
        records_per_cell=NBAR * math.log(lx2 / P_MID),
        c_surv_log2=math.log2(C_NORMS + C_TDIV + c_band_f + c_band_g + c_ecm_f + c_ecm_g), mean_slice_log2=math.log2(population / (n_a_upper * len(rows))))
  if return_work: out['work_contribs'] = work_list
  return out
FILT_TABLE = {   # signed sigma-quotient system at scale 2^17, excess 1.1: (seed, U/N, N, max coef, weight)
  1.1: [(1, 18.181, 7209, 44, 100.05), (2, 19.078, 6870, 36, 100.02), (3, 18.401, 7123, 41, 100.03), (4, 19.571, 6697, 31, 100.00)],
}
W_ROW = 201.0   # mean row weight bound (merge target 200): total entries = N W_ROW
W200_RATIO = 0.80   # N(200)/N(100) charged (filtsim.run maxima)
W200 = {   # max N(200)/N(100) per log2 U at excess 1.1
  12: 0.796, 13: 0.793, 14: 0.793, 15: 0.788, 16: 0.785, 17: 0.778, 18: 0.777, }
SCALE_NU = {   # N/U at weight 100, excess 1.1, by log2 U (filtsim.run, kmax 64, rel means)
  12: (0.08472, 0.08472, 0.07397, 0.08643, 0.07544, 0.09033, 0.09668, 0.08374, 0.10913, 0.08569),
  13: (0.06213, 0.0675, 0.06665, 0.0553, 0.06396, 0.07227),
  14: (0.06207, 0.06421, 0.06201, 0.05872, 0.06189, 0.07104),
  15: (0.06, 0.05756, 0.0585, 0.05872, 0.05164, 0.06155),
  16: (0.058, 0.05525, 0.05991, 0.05351, 0.05653, 0.05766),
  17: (0.055, 0.05241, 0.05434, 0.05109),
  18: (0.05187, 0.05138),
}
C_ENTRY = 20.7   # per entry 4 + 3 + 0.88*11 + 0.12*33 = 20.64 (STORAGE.md)
C_ROW = 2710.0   # row tail incl. reduction, counted 2708 by [10]
C_RED9 = 420.0   # 9-limb accumulator + add, counted 385 by [10]
UNIT_FRAC_CHARGED = 0.88   # H4u, below every measured value at weight 100 and 200 (unit_entries.py)
BW_M = 19
BW_N = 4
N_SM = 21                   # dense Schirokauer columns (unit ranks 18+3; [13])
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
  log2_N = log2_U + math.log2(ratio * pad_n * W200_RATIO)
  # Block Wiedemann on square A = [M | S] (kernel = virtual-log line under H3), random m x n blocks:
  # Krylov N (n/m + 1) products (W_ROW sparse at C_ENTRY, N_SM dense at C_FULL per row);
  # projections (m + n) N^2 terms at C_FULL; one Horner reconstruction;
  # every product reduces each row in C_ROW; Horner/projection entries reduce at C_RED9
  c_prod = W_ROW * C_ENTRY + C_ROW + N_SM * C_FULL
  per_n2 = c_prod * (BW_N / BW_M + 1) + (1 / BW_M + 1 / BW_N) * (c_prod + BW_N * C_FULL + C_RED9) + C_FULL * (BW_M + BW_N)
  per_n1 = (BW_M + BW_N) * C_RED9
  log2_la = math.log2((2 ** (2 * log2_N) * per_n2 + 2 ** log2_N * per_n1) / UNIT)
  return dict(ratio=ratio, log2_N=log2_N, log2_la=log2_la, log2_gen=generator(log2_N))
def generator(log2_N):
  # (m+n) x (m+n) matrix generator (L = (N+N_SM)(1/m+1/n) + 2) by approximant bases:
  # <= 4 (m+n)^3 L^log2(3) log2(L) mults at C_FULL + C_RED9; log2 units
  k = BW_M + BW_N
  log2_L = math.log2((2 ** log2_N + N_SM) * (1 / BW_M + 1 / BW_N) + 2)
  return math.log2(4 * k ** 3 * (C_FULL + C_RED9) / UNIT) + math.log2(3) * log2_L + math.log2(log2_L)
def total_log2(log2_relcol, log2_la, log2_desc=DESCENT_CAP, log2_gen=0.0):
  return math.log2(2 ** log2_relcol + 2 ** log2_la + 2 ** log2_desc + 2 ** log2_gen + 2 ** LOG2_SM + 2 ** LOG2_FBC
           + 2 ** LOG2_CANON + 2 ** LOG2_MISC + 2 ** LOG2_VERIFY + 2 ** 50)
