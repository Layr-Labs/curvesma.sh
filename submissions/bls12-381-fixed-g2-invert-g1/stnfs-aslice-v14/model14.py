"""stnfs-aslice-v14 ledger model: relation collection with lean record routing
and certified block selection, LA at 88 instructions/entry, descent cap 2^111.

Changes vs the v13 model12.py (all constants and structure documented in
proof.md/STORAGE.md/DESCENT.md of this package):
 - record routing: per-line bounded buffers + two-level page directory (STORAGE.md),
   20 instructions per retained record + 8 per cell (apply/recovery separate), replacing 256+8;
 - block selection: blocks whose certified upper bound (embedding_bounds, 64-cell
   blocks) exceeds the selection cut are dropped before sieving; the all-cell
   costs (norm bounds, layout, pair setups) are still paid for them. The retained
   cells/yield fractions are the measured curve of stats14.json (regenerated
   per-row from the pool by the declared experiment at reduced scope).
 - LA: 88 instructions/entry (16 blocks, four-hop partial sums), BW=3, W=100.
 - descent: cap 2^119; the recurrence is descent14.py (window-sieved descent).
"""
import math, json, struct, zlib, base64
from fractions import Fraction as F
from common12 import rho, integer_metric, quadratic
import filtsim
from threshold_model import Retention

LN2 = math.log(2)
UNIT = 20000.0

OP = dict(lb=69.0, lx2=63.0, excess=1.3, slack_f=292, slack_g=220, width_bits=8, L0=64,
          theta=F(1, 2), alpha_f=5.0, alpha_g=1.0)
MARGIN = 0.30
PAD_N = 2.0
DESCENT_CAP = 119.0
ECM_FIND = 0.99

# --- lean record routing (STORAGE.md) and block selection (proof.md 4.7) ---
REC_ROUTE = 20.0        # instructions per retained record (append 6 + page 4 +
                        # directory 2 + sentinel/dispatch 8; apply/recovery are charged
                        # separately at 7.00/7.003 per record)
CELL_ROUTE = 8.0        # instructions per cell (row bookkeeping, sentinels)
SEL_CUT = 1285.0        # selection cut on the certified (log2 N_f + log2 N_g) upper bound
# measured selection curve at 64-cell blocks, width 128 (stats14.json):
SEL_CELLS = 0.8078      # fraction of cells kept at cut 1285 (stats14.json)
SEL_YIELD = 0.9617      # fraction of retained yield kept at cut 1285 (stats14.json)

LINE_TAB = [(13.82, 24.10), (9.82, 31.11), (5.82, 88.31), (1.82, 310.70)]
REC_TAB = [(18.91, 44.90), (10.91, 78.21), (3.91, 266.97), (0.91, 938.08), (-1.09, 2231.66)]
C_COSET_EMPTY = 1600.0
C_APPLY = 7.00
C_RECOVER = 7.003
C_SCAN_CELL, C_SCAN_LEAF = 8.0, 21.0
C_BUILD = 12.325
C_ZERO = 0.625
C_SETUP6 = 406.87 + 72.0 + 21.13
C_NB_EVAL = 5900.0
C_NB_ROW = 4200.0
C_NB_SLICE = 15300.0
MONT5 = 367.76
C_NORMS = 2.0 ** 22
C_TDIV = 2.0 ** 24.5
NBAR = 2.0
LX0, P_MID = 8.0, 20.0
P_MIN = 256

ECM_TAB = {20: (11e3, 86), 25: (5e4, 214), 30: (25e4, 430), 35: (1e6, 904), 40: (3e6, 2350), 45: (11e6, 4480),
           50: (43e6, 7553), 55: (11e7, 17769), 60: (26e7, 42017), 65: (85e7, 69408), 70: (29e8, 102212),
           75: (76e8, 188056), 80: (25e9, 265557)}


def ecm_instr(fbits, nbits, success=ECM_FIND):
    d = fbits / 3.3219
    ks = sorted(ECM_TAB)
    if d <= ks[-1]:
        k = min(x for x in ks if x >= d); b1, cur = ECM_TAB[k]
    else:
        b1, cur = ECM_TAB[ks[-1]]; e = (d - ks[-1]) / 5; b1 *= 3.3 ** e; cur *= 1.5 ** e
    curves = cur * math.log(1 / (1 - success))
    L = max(5, math.ceil(nbits / 64)); mulmod = MONT5 * (L / 5.0) ** 2
    return curves * 2 * 1.44 * b1 * 11 * mulmod


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
    raw = zlib.decompress(base64.b85decode(''.join(text.split())))
    hdr, blob = raw.split(b'\n', 1)
    hdr = json.loads(hdr)
    fmt = '<BHHHHBBBBB'; n = struct.calcsize(fmt)
    rows = []
    for i in range(0, len(blob), n):
        seed, index, qa16, nf16, ng16, flags, evals, leaves, blen, rlen16 = struct.unpack(fmt, blob[i:i + n])
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


def evaluate(rows, n_total, geo, op=OP, retain=None, cost=None, scale_consts=1.0, verbose=False):
    """v14 relation collection: lean routing + block selection.

    Cost split: the all-cell part (norm-bound blocks, layout, pair setups, zero)
    is paid for every cell of the region; the sieve part (line path, record walk,
    apply/recover, routing, scan, cofactorisation) only for cells of retained
    blocks. The retention fractions SEL_CELLS/SEL_YIELD are the measured curve.
    """
    lb, lx2, excess = op['lb'], op['lx2'], op['excess']
    Sf, Sg, W = op['slack_f'], op['slack_g'], op['width_bits']
    retain = retain or Retention(int(lb), int(lx2), step=.25, max_bits=850, tilt=.05)
    cost = cost or Retention(850, int(lx2), step=.25, max_bits=851, tilt=.005)

    wsum = 0.0; ysum = 0.0; y2 = 0.0; rej_w = 0.0; nf_m = ng_m = 0.0
    work_sum = 0.0
    c_ecm = ecm_instr(lb, Sf + 8) + ecm_instr(lb, Sg + 8)
    c_surv = C_NORMS + C_TDIV + c_ecm
    pairs_all = (NBAR / 2) * 2 ** log2_li(lx2)
    for r in rows:
        w = (1 - r['qa']) ** 3
        wsum += w
        if r['primitive'] and not r['rejected']:
            rf = retain.at(r['nf'], Sf)['retention']; rg = retain.at(r['ng'], Sg)['retention']
            uf_r, ug_r = (r['nf'] + op['alpha_f']) / lb, (r['ng'] + op['alpha_g']) / lb
            miss = ECM_FIND ** (lam_top(uf_r, lx2 / lb) + lam_top(ug_r, lx2 / lb))
            y = rho(uf_r) * rho(ug_r) * rf * rg * miss
            ysum += w * y; y2 += w * w * y * y
        if r['rejected']:
            rej_w += w
        nf_m += w * r['nf']; ng_m += w * r['ng']

        S = geo.slice_cells(r['qa']); K = geo.lines(r['qa'])
        # all-cell part
        p_all = dict(build=C_BUILD, zero=0.0, setup6=pairs_all * C_SETUP6 / S,
                     thr=C_NB_EVAL * r['evals'] / r['blen'] + C_NB_ROW / r['rlen'] + C_NB_SLICE / S)
        # kept-cell part (line, records, routing, scan, cofactorisation, zero)
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
        per_cell = sum(p_all.values()) + SEL_CELLS * sum(p_keep.values())
        per_cell = per_cell * scale_consts if scale_consts != 1.0 else per_cell
        work_sum += S * per_cell

    n_a_upper = geo.va / n_total
    population = n_a_upper * sum(geo.slice_cells(r['qa']) for r in rows)
    relcol_instr = n_a_upper * work_sum
    per_cell_mean = relcol_instr / population
    # retained yield at the selection cut
    rels = n_a_upper * geo.vb * geo.lower * ysum * SEL_YIELD
    yield_se = math.sqrt(max(0.0, y2 / n_total - (ysum / n_total) ** 2) / n_total) / (ysum / n_total) / LN2 if ysum > 0 else float('inf')
    nf_m /= wsum; ng_m /= wsum

    log2_U = math.log2((NBAR / 2) * 2 ** log2_li(lb) + 2 ** 31)
    log2_need = log2_U + math.log2(excess)
    log2_rels = math.log2(rels) if rels > 0 else float('-inf')
    feasible = log2_rels >= log2_need + MARGIN
    log2_relcol = math.log2(relcol_instr / UNIT)
    out = dict(op={k: (str(v) if isinstance(v, F) else v) for k, v in op.items()},
               log2_population=math.log2(population), per_cell=per_cell_mean,
               log2_relcol=log2_relcol, sel_cells=SEL_CELLS, sel_yield=SEL_YIELD,
               log2_rels=log2_rels, log2_U=log2_U, log2_need=log2_need,
               margin=log2_rels - log2_need, yield_se_log2=yield_se, feasible=feasible,
               rejected_weight=rej_w / wsum, nf_mean=nf_m, ng_mean=ng_m,
               records_per_cell=NBAR * math.log(lx2 / P_MID),
               c_surv_log2=math.log2(c_surv), mean_slice_log2=math.log2(population / (n_a_upper * len(rows))))
    return out


FILT_TABLE = {
    1.2: [(1, 16.258, 8062, 26, 100.02), (2, 17.251, 7598, 25, 100.02)],
    1.3: [(1, 18.600, 7047, 26, 100.03), (2, 19.253, 6808, 21, 100.00)],
    1.4: [(1, 20.789, 6305, 29, 100.04), (2, 20.687, 6336, 29, 100.02)],
    1.6: [(1, 23.072, 5681, 34, 100.03), (2, 24.068, 5446, 23, 100.00)],
    1.8: [(1, 24.019, 5457, 28, 100.05), (2, 24.504, 5349, 24, 100.01)],
}
W_ROW = 100.0
C_ENTRY = 88.0
BW = 3.0
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
    log2_la = math.log2(BW) + 2 * log2_N + math.log2(W_ROW * C_ENTRY / UNIT)
    return dict(ratio=ratio, log2_N=log2_N, log2_la=log2_la)


def total_log2(log2_relcol, log2_la, log2_desc=DESCENT_CAP):
    return math.log2(2 ** log2_relcol + 2 ** log2_la + 2 ** log2_desc + 2 ** LOG2_SM + 2 ** LOG2_FBC
                     + 2 ** LOG2_CANON + 2 ** LOG2_MISC + 2 ** LOG2_VERIFY)


def memory_report(e):
    S = 2 ** e['largest_slice_log2'] if 'largest_slice_log2' in e else 2 ** 67.5
    seg = 2 ** 30
    cells = seg; meta = 8 * e.get('L4_largest', 2 ** 49); recs = 16 * e['records_per_cell'] * seg
    return dict(cells_log2=math.log2(cells), meta_log2=math.log2(meta), records_log2=math.log2(recs),
                job_log2=math.log2(cells + meta + recs), segments='per-line bounded, 2^30-cell sum arrays')
