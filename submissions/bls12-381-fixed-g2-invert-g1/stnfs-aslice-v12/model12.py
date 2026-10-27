"""stnfs-aslice-v12 cost and yield model (stdlib only; imported by evidence12.py): population, retained relation orbits, per-cell
cost per pool row from the measured kernel constants (KERNELS.md), relation collection, filtering -> N, linear algebra, ledger."""
import math, json, struct, zlib, base64
from fractions import Fraction as F
from common12 import rho, integer_metric, quadratic
import filtsim                 # tabulated Dickman rho for the many evaluations of lam_top
from threshold_model import Retention
LN2 = math.log(2)
UNIT = 20000.0                 # one dense Fp12 multiplication (tracks/classical-group-work-v2.json)

# ------------------------------------------------------------------ operating point ------------------------------
OP = dict(lb=69.0, lx2=63.0, excess=1.3, slack_f=292, slack_g=220, width_bits=8, L0=64, theta=F(1, 2), alpha_f=5.0, alpha_g=1.0)
MARGIN = 0.30                  # required log2 margin of the modelled retained relations over excess x U
PAD_N = 2.0                    # pad on the simulated filtered dimension (proof.md section 5)
DESCENT_CAP = 106.0            # metered cap on both descents (log2 units)
ECM_FIND = 0.99                # per-prime ECM success at the full curve budget; the yield is charged 0.99^(expected number of primes in (X2, B] on both sides)

# ------------------------------------------------------------------ measured constants (KERNELS.md) -------------
# line path (p < 2^20): Ir per bucket update INCLUDING the per-(pair, line, coset) setup, by log2 h, h = points per line coset
LINE_TAB = [(13.82, 24.10), (9.82, 31.11), (5.82, 88.31), (1.82, 310.70)]
# record path (2^20 <= p < 2^63): Ir per hit (6-dim walk, fold, record store, per-coset recentring/centre), h = points per slice coset
REC_TAB = [(18.91, 44.90), (10.91, 78.21), (3.91, 266.97), (0.91, 938.08), (-1.09, 2231.66)]
C_COSET_EMPTY = 1600.0         # per (slice, pair) when h < 2^-1.09: the exact enumeration of an (almost always) empty coset (measured 943-1444 at h = 2^-1..2^-3; 2231.66 x 2^-1.09 = 1050)
C_APPLY = 7.00                 # per record: apply pass (load cell, add batch log)
C_RECOVER = 7.003              # per record: recovery pass (load cell, sentinel test)
C_SCAN_CELL, C_SCAN_LEAF = 8.0, 21.0   # scan: 9.313/cell at 16-cell blocks and 13.250 at 4-cell blocks = 8.0 per cell + 21 per threshold block (leaf)
C_BUILD = 12.325               # per cell, whole-program build of the packed layout (allocation, metadata, rows, first zero pass) at the benchmark geometry
C_ZERO = 0.625                 # per cell, one zero pass over both sum arrays (charged again on top of the build's)
C_SETUP6 = 406.87 + 72.0 + 21.13   # per (slice, pair): a(s) mod p (Horner/Barrett), c = a(s) r^-1 (max of three measured primes) + 36 streamed words + rounding to 500
C_NB_EVAL = 5900.0             # per norm-bound block evaluation (nbound12: both sides, four places, outward intervals, adaptive bookkeeping; measured 5862.7 on 16-blocks, 5833.0 on 64/128 blocks)
C_NB_ROW = 4200.0              # per row: nbound12 row_setup (suffix evaluation at the four places; measured 4120)
C_NB_SLICE = 15300.0           # per slice: nbound12 setup_a (a(v) and its powers at the four places; measured 19342.8 - 4120 per case)
MONT5 = 367.76                 # ecm12 mont_mul, 5 x 64-bit limbs (cofactors < 2^300)
C_NORMS = 2.0 ** 22            # exact norms of one survivor (two resultants over O_h, <= 1302 and 520 bits)
C_TDIV = 2.0 ** 24.5           # trial division of both norms by the 82,025 primes below 2^20 (<= 21 limbs x ~5 Ir per limb per prime x 2), valuations of the recorded ideals, degree-1 checks
NBAR = 2.0                     # degree-1 prime ideals per rational prime, both sides: Chebotarev mean 1 per side; evidence12 part [8] measures 1.996 over [2^8, 2^12) (v1-v10 charged 1.947)
LX0, P_MID = 8.0, 20.0         # primes below 2^8 are not sieved (omitted mass); the line path ends and the record path starts at 2^20
P_MIN = 256

ECM_TAB = {20: (11e3, 86), 25: (5e4, 214), 30: (25e4, 430), 35: (1e6, 904), 40: (3e6, 2350), 45: (11e6, 4480),
           50: (43e6, 7553), 55: (11e7, 17769), 60: (26e7, 42017), 65: (85e7, 69408), 70: (29e8, 102212), 75: (76e8, 188056), 80: (25e9, 265557)}


def ecm_instr(fbits, nbits, success=ECM_FIND):
    """ECM to 2^fbits on an nbits-bit cofactor, `success` per prime: GMP-ECM B1/curve table x ln(1/(1-success)) curves, stage 1 =
    1.44 B1 bits x 11 modular multiplications, stage 2 as a second stage 1, one multiplication = MONT5 (L/5)^2 Ir (L >= 5 limbs)."""
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
    """expected number of prime factors with log_B size in (t0, 1] of a B-smooth integer with u = log N / log B (Dickman-conditioned density kappa(u, t)/t)"""
    r = filtsim.rho(u)
    if r <= 0 or t0 >= 1:
        return 0.0
    return sum(filtsim.rho(u - t) / r / t * ((1 - t0) / m) for t in (t0 + (1 - t0) * (i + 0.5) / m for i in range(m)))


def log2_li(log2x):
    x = log2x * LN2
    return log2x - math.log2(x) + math.log2(1 + 1 / x + 2 / x / x)


def tab(t, log2h):
    """piecewise-linear in log2 h; clamped to the dense end (the per-update cost decreases with h), extended at the sparse end"""
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
        # rlen = 16 floor(row length / 16) <= the true row length: the per-row setup charge 4200 / rlen is an upper bound
        rows.append(dict(seed=seed, index=index, qa=qa16 / 65536.0, nf=nf16 / 16.0, ng=ng16 / 16.0, primitive=bool(flags & 1),
                         rejected=bool(flags & 2), nonfinite=bool(flags & 4), evals=evals, leaves=leaves, blen=blen, rlen=16 * rlen16))
    assert len(rows) == hdr['rows']
    return hdr, rows


class Geometry:
    """Exact metrics of config.json; continuous counts of the eligible a and of the b half-slices, with the certified count factors."""

    def __init__(self, cfg, theta, lower=F(2401, 2500), upper=F(2601, 2500)):
        self.A2 = cfg['production_A_squared']
        self.qa, self.da, self.ia = integer_metric(cfg['matrix_a'])
        self.qb, self.db, self.ib = integer_metric(cfg['matrix_b'])
        from ideal_canonical import determinant
        det_a = F(determinant(self.qa), self.da ** 6); det_b = F(determinant(self.qb), self.db ** 6)
        v6 = math.pi ** 3 / 6
        self.va = v6 * (float(theta) * self.A2) ** 3 / (2 * math.sqrt(float(det_a)))     # eligible a (a5 > 0, Q_a <= theta A^2), continuous
        self.vb = v6 * self.A2 ** 3 / (2 * math.sqrt(float(det_b)))                        # largest half-slice (Q_a = 0), continuous
        self.lower, self.upper = float(lower), float(upper)
        self.rb = [math.sqrt(float(self.ib[j][j]) * self.A2 * self.db) for j in range(6)]    # box radii of the largest slice (b5: lines)
        self.inv55 = float(self.ib[5][5]) * self.db

    def slice_cells(self, qa_frac):            # upper count of cells of the half-slice of a (Q_a = qa_frac A^2)
        return self.vb * (1 - qa_frac) ** 3 * self.upper

    def lines(self, qa_frac):                  # K = floor(sqrt((A^2 - Q_a) (Q_b^-1)_55))
        return max(1.0, math.floor(math.sqrt((1 - qa_frac) * self.A2 * self.inv55)))

    def suffix_entries(self, qa_frac):         # L4 = prod_{i=1..4} (2 r_i + 1) of the slice (dense metadata per line)
        s = math.sqrt(1 - qa_frac); L4 = 1.0
        for j in range(1, 5):
            L4 *= 2 * math.floor(self.rb[j] * s) + 1
        return L4


def evaluate(rows, n_total, geo, op=OP, retain=None, cost=None, scale_consts=1.0, pad_n=PAD_N, verbose=False):
    """The model at an operating point over the pool rows (n_total = pool denominator 20,000)."""
    lb, lx2, excess = op['lb'], op['lx2'], op['excess']
    Sf, Sg, W = op['slack_f'], op['slack_g'], op['width_bits']
    retain = retain or Retention(int(lb), int(lx2), step=.25, max_bits=850, tilt=.05)
    cost = cost or Retention(850, int(lx2), step=.25, max_bits=851, tilt=.005)
    # ---- population and yield ----
    wsum = 0.0; ysum = 0.0; y2 = 0.0; rej_w = 0.0; nf_m = ng_m = 0.0
    work_sum = 0.0
    c_ecm = ecm_instr(lb, Sf + 8) + ecm_instr(lb, Sg + 8)          # cofactors are <= slack + W bits (the preliminary test)
    c_surv = C_NORMS + C_TDIV + c_ecm
    pairs_all = (NBAR / 2) * 2 ** log2_li(lx2)                      # sigma-pairs with p < X2 (record and line paths)
    parts = dict(build=0.0, zero=0.0, line=0.0, rec=0.0, apply=0.0, recover=0.0, empty=0.0, setup6=0.0, scan=0.0, thr=0.0, cof=0.0, meta=0.0)
    for r in rows:
        w = (1 - r['qa']) ** 3
        wsum += w
        # yield: retained relation orbits per cell of this slice (zero for non-primitive pairs and rejected sample points)
        if r['primitive'] and not r['rejected']:
            rf = retain.at(r['nf'], Sf)['retention']; rg = retain.at(r['ng'], Sg)['retention']
            uf_r, ug_r = (r['nf'] + op['alpha_f']) / lb, (r['ng'] + op['alpha_g']) / lb
            miss = ECM_FIND ** (lam_top(uf_r, lx2 / lb) + lam_top(ug_r, lx2 / lb))   # E[0.99^K] >= 0.99^E[K] (convexity)
            y = rho(uf_r) * rho(ug_r) * rf * rg * miss
            ysum += w * y; y2 += w * w * y * y
        if r['rejected']:
            rej_w += w
        nf_m += w * r['nf']; ng_m += w * r['ng']
        # cost per cell of this slice
        S = geo.slice_cells(r['qa']); K = geo.lines(r['qa']); L4 = geo.suffix_entries(r['qa'])
        p = dict(build=C_BUILD, zero=C_ZERO, meta=0.0)
        # line path: h = line points / p; NBAR / p updates per cell
        lp = LX0; line = 0.0
        while lp < P_MID:
            lp2 = min(lp + 0.25, P_MID); dsum = math.log(lp2 / lp)         # sum of 1/p over (2^lp, 2^lp2] = ln(lp2/lp) (Mertens)
            line += NBAR * dsum * tab(LINE_TAB, math.log2(S / K) - 0.5 * (lp + lp2)); lp = lp2
        p['line'] = line
        # record path: h = slice points / p
        lp = P_MID; rec = 0.0; empty = 0.0
        while lp < lx2:
            lp2 = min(lp + 0.25, lx2); dsum = math.log(lp2 / lp)
            log2h = math.log2(S) - 0.5 * (lp + lp2)          # h = cells / p: one bit below the benchmark's convention (2 cells / p), i.e. charged at a sparser regime
            if log2h >= REC_TAB[-1][0]:
                rec += NBAR * dsum * tab(REC_TAB, log2h)
            else:                                                       # per (slice, pair) empty-coset enumeration
                dpairs = (NBAR / 2) * (2 ** log2_li(lp2) - 2 ** log2_li(lp))
                empty += dpairs * C_COSET_EMPTY / S
            lp = lp2
        hits = NBAR * math.log(lx2 / P_MID)                             # records per cell (sum of NBAR / p over the record-path primes)
        p['rec'] = rec; p['empty'] = empty; p['apply'] = hits * C_APPLY; p['recover'] = hits * C_RECOVER
        p['setup6'] = pairs_all * C_SETUP6 / S
        p['scan'] = C_SCAN_CELL + C_SCAN_LEAF * r['leaves'] / r['blen']
        p['thr'] = C_NB_EVAL * r['evals'] / r['blen'] + C_NB_ROW / r['rlen'] + C_NB_SLICE / S
        sf = cost.survivor_upper(r['nf'], Sf, W); sg = cost.survivor_upper(r['ng'], Sg, W)
        p['cof'] = sf * sg * c_surv
        per_cell = sum(v for k, v in p.items() if k != 'cof') * scale_consts + p['cof']
        work_sum += S * per_cell
        for k, v in p.items():
            parts[k] += S * v
    # expected counts: eligible a = va (x count factors); canonical fraction over the pool = (#rows)/n_total with weights
    # the certified count factors (lattice_count_bounds: a-count x b-count in [2401/2500, 2601/2500]) are applied once, on the
    # cells of a slice (upper, work) and on the relation count (lower, yield); every pool sample represents va / n_total eligible a
    n_a_upper = geo.va / n_total
    population = n_a_upper * sum(geo.slice_cells(r['qa']) for r in rows)
    relcol_instr = n_a_upper * work_sum
    per_cell_mean = relcol_instr / population
    rels = n_a_upper * geo.vb * geo.lower * ysum
    yield_se = math.sqrt(max(0.0, y2 / n_total - (ysum / n_total) ** 2) / n_total) / (ysum / n_total) / LN2 if ysum > 0 else float('inf')
    nf_m /= wsum; ng_m /= wsum
    # unknowns: sigma-orbits of the prime ideals of relative degree 1 over K_h with norm <= B on both sides: (NBAR/2) Li(B) over the
    # degree-1 primes of K_h plus those over K_h-primes of degree >= 2 (norm p^f <= B: < 2 Li(B^(1/2)) < 2^31 of them, charged as 2^31)
    log2_U = math.log2((NBAR / 2) * 2 ** log2_li(lb) + 2 ** 31)
    log2_need = log2_U + math.log2(excess)
    log2_rels = math.log2(rels)
    feasible = log2_rels >= log2_need + MARGIN
    log2_relcol = math.log2(relcol_instr / UNIT)
    out = dict(op={k: (str(v) if isinstance(v, F) else v) for k, v in op.items()}, log2_population=math.log2(population), per_cell=per_cell_mean,
               parts={k: v * n_a_upper / population for k, v in parts.items()}, log2_relcol=log2_relcol,
               log2_rels=log2_rels, log2_U=log2_U, log2_need=log2_need, margin=log2_rels - log2_need, yield_se_log2=yield_se,
               feasible=feasible, rejected_weight=rej_w / wsum, nf_mean=nf_m, ng_mean=ng_m, uf=nf_m / lb, ug=ng_m / lb,
               c_surv_log2=math.log2(c_surv), mean_slice_log2=math.log2(population / (n_a_upper * len(rows))),
               largest_slice_log2=math.log2(geo.slice_cells(0.0)), lines_largest=geo.lines(0.0), L4_largest=geo.suffix_entries(0.0),
               pairs_log2=math.log2(pairs_all), records_per_cell=NBAR * math.log(lx2 / P_MID))
    return out


# filtering simulation results (filtsim.py, this profile u_f = 11.12, u_g = 7.29, B = 2^69, scale 2^17 sampled columns,
# kmax 64, dense columns propagated): excess -> [(seed, U/N at W=100, N, max |coef|, full row weight), ...]   @@FILT_TABLE@@
FILT_TABLE = {
    1.2: [(1, 16.258, 8062, 26, 100.02), (2, 17.251, 7598, 25, 100.02)],
    1.3: [(1, 18.600, 7047, 26, 100.03), (2, 19.253, 6808, 21, 100.00)],
    1.4: [(1, 20.789, 6305, 29, 100.04), (2, 20.687, 6336, 29, 100.02)],
    1.6: [(1, 23.072, 5681, 34, 100.03), (2, 24.068, 5446, 23, 100.00)],
    1.8: [(1, 24.019, 5457, 28, 100.05), (2, 24.504, 5349, 24, 100.01)],
}
W_ROW = 100.0                  # charged row weight of the filtered matrix
C_ENTRY = 36.0                 # la_mac.s: 33 per entry (64-bit column indices) + row finish and 4-block combination amortised over W_ROW
BW = 3.0                       # block Wiedemann: ~3 N single-vector products (sequence, generator application)
LOG2_SM = 88.0                 # Schirokauer maps: <= 2^37 Ir per relation x 2^64.5 relations / 2^14.29
LOG2_FBC = 72.0                # factor base: 2^64.4 ideals x (roots mod p 2^20 + 6-dim reduction 2^13 + tables) Ir / 2^14.29
LOG2_CANON = 90.0              # canonical a-list: bulk HNF-key/segmented-sort route, analytical bound (CANONICAL.md / proof.md 4.2)
LOG2_MISC = 64.0               # relation store, dedupe, purge/merge passes, matrix storage, challenge pairing, randomness
LOG2_VERIFY = 11.0


def ratio_N_over_U(excess):
    """max over seeds of N/U at W = 100 from FILT_TABLE at this excess (interpolated linearly in excess)"""
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
    return math.log2(2 ** log2_relcol + 2 ** log2_la + 2 ** log2_desc + 2 ** LOG2_SM + 2 ** LOG2_FBC + 2 ** LOG2_CANON + 2 ** LOG2_MISC + 2 ** LOG2_VERIFY)


def memory_report(e):
    """peak storage of one a-slice job (bytes, log2): cells (2 x uint16), suffix metadata (int64 per entry per line), records (16 bytes)"""
    S = 2 ** e['largest_slice_log2']
    cells = 4 * S; meta = 8 * e['L4_largest'] * e['lines_largest']; recs = 16 * e['records_per_cell'] * S
    return dict(cells_log2=math.log2(cells), meta_log2=math.log2(meta), records_log2=math.log2(recs), job_log2=math.log2(cells + meta + recs),
                segments=math.ceil(cells / 2 ** 62))
