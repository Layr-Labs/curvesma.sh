"""bn1 cost model: per-cell collection work and yield on exact pool rows; batch tests; small rows.  UNIT = 20,000 instructions."""
import math
from math import ceil, exp, floor, gcd, log, log2, prod, sqrt
from fractions import Fraction as F
from common12 import rho, integer_metric, regenerate_pair, norms, mul, norm, integer_interval, sieve_coords, sieve_metric
import ideal_canonical as IC
import embedding_bounds as EB
from threshold_model import Retention
LN2 = log(2); UNIT = 20000.
OP = dict(lb=56., lx2=56., excess=1.095, slack_f=60., width_bits=32, theta=F(3, 8), alpha_f=-.40, alpha_g=2.76)
MARGIN = .06; MARGIN_SE = 2.; Q_CORE = .02
KAP2 = (.07, .07)                              # >= kappa2(u) on every yielding row ([14k])
def a2(r, s, op=OP):
    """Bits added to a norm whose exact 2-part 2^v is known: the odd part n - v is priced as a random integer of size n - v + 1 + kappa
    (a random integer has v_2 = k with probability 2^-(k+1)), plus alpha over the odd primes."""
    return op['alpha_' + 'fg'[s]] + 1. + KAP2[s] - r['v2' + 'fg'[s]]
def kappa2(u, lb=56.):
    lo, hi, r0 = .0, 1., rho(u)
    for _ in range(40):
        m = (lo + hi) / 2
        if sum(2. ** -(k + 1) * rho(u - (k - 1 - m) / lb) for k in range(60)) > r0: lo = m
        else: hi = m
    return hi
# --- kernel prices (KERNELS.md): carried tables are read only inside their measured density ranges
CANON_WORK_ALLOWANCE = 1.02; WORK_SE_K = 2.; Q_F_SAFETY = 100.
REC_ROUTE, REC_ROUTE7 = 10.01, 3.01; CELL_ROUTE = .02
LINE_TAB = [(18.78, 6.06), (13.82, 9.04), (9.82, 14.26), (5.82, 36.), (1.83, 174.93)]
REC_TAB = [(18.91, 23331674 / 522256), (10.91, 3183663 / 41121), (3.91, 1638258 / 6539), (.91, 3216517 / 3967), (-1.09, 4964184 / 2901)]
C_COSET_EMPTY = 1600.; C_APPLY = 5.76
REC7_TAB = [(25., 33.17), (24.01, 37.65), (19., 37.65), (18.01, 38.33), (16.02, 38.33), (16., 43.76), (14.49, 49.66), (14.01, 54.79), (13.01, 54.95), (12.98, 60.33), (12.03, 60.33), (12., 62.88), (11.02, 62.88), (11.01, 73.74), (10.55, 73.74), (10.31, 82.85), (10., 86.39), (9.7, 87.26), (9.55, 87.26), (8.65, 127.16), (7.47, 127.16), (7., 167.75), (5.5, 243.83), (4.99, 243.83)]
REC7_K = [(4.99, 1.)]                         # set from the BN254 campaign (KERNELS.md): (stratum lower edge, mean + 2 SE of Ir / (H x ref))
P7, X7, EXTRA7, RUN_ALLOW, RUN_T = 40., 56., 2., .006, (16, 1024)
COUNT = .9905, 1.002                              # lattice points per volume of the half-ellipsoids: (lower a x b, upper a x b), [3]
C_SCAN_CELL, C_SCAN_LEAF = 1.125, 42.88; C_BUILD_E, C_BUILD_SLICE, C_BUILD_SEG = 28.57, 2. ** 21, 2. ** 17; C_SETUP6 = 406.87 + 72. + 21.13
NB_ALLOW = 1.25                                # allowance on the norm-bound kernel prices (the kernel runs with BN254 coefficients, KERNELS.md)
C_NB_EVAL = 1178. * NB_ALLOW; C_NB_ROW = 378. * NB_ALLOW; C_NB_SLICE = 18500. * NB_ALLOW
MONT5 = 367.76; C_NORMS = 27000.; NBAR = 2.; SIEVE_SIDES = .5
LX0, P_MID = 14., 40.; PRIME_MIN, EXPLICIT_BITS = 2 ** 14, 14
NORM_BITS = 576                                # both norms are below 2^576 on the region ([5t]): 9 limbs
TD_BITS, BAND_TOP = 20, 40                     # trial division below 2^20; strip 1 on (2^20, 2^40]; strip 2 on (2^40, B]
MAXBITS = 600
def hc_sup(t, xm):
    if xm < t[-1][0]: return .0
    cand = [xm] + [h for h, c in t if h <= xm]
    for (h1, c1), (h2, c2) in zip(t, t[1:]):
        s = (c1 - c2) / (h1 - h2); x = h2 + (-s / log(2) - c2) / s
        if h2 <= x <= min(h1, xm): cand.append(x)
    return max(2 ** x * tab(t, x) for x in cand)
def contam_sup(geo):
    """Upper bound of the slice-walk work of one slice on [2^P7, 2^X7) (band maxima of hits x price, empty cosets, setup6)."""
    S = geo.slice_cells(.0); lp = P7; d = .0
    while lp < X7:
        lp2 = min(lp + .25, X7); m = .5 * (lp + lp2); dp = NBAR / 2 * (2 ** log2_li(lp2) - 2 ** log2_li(lp))
        d += NBAR * log(lp2 / lp) * 2 ** m * hc_sup(REC_TAB, log2(S) - m) + dp * (C_COSET_EMPTY + C_SETUP6); lp = lp2
    return d
def step7(l, m=1): return tab(REC7_TAB, l) * (next(k for e, k in REC7_K if l >= e) if m else 1)
def log2_fb(lb): return log2(NBAR / 2 * 2 ** log2_li(lb) + 2 ** 31) - 1
def log2_li(log2x): x = log2x * LN2; return log2x - log2(x) + log2(1 + 1 / x + 2 / x / x)
def tab(t, log2h):
    if log2h >= t[0][0]: return t[0][1]
    for (h1, c1), (h2, c2) in zip(t, t[1:]):
        if log2h >= h2: return c2 + (c1 - c2) * (log2h - h2) / (h1 - h2)
    return t[-1][1]
def coprime_pair(a, b):
    d = gcd(norm(a), norm(b))
    if d == 1: return True
    rows = [mul(v, [0] * k + [1] + [0] * (5 - k)) for v in (a, b) for k in range(6)]; key = IC.row_hnf(rows, d); return all(key[k][k] == 1 for k in range(6))
_CTX = {}
def _ctx(cfg):
    if not _CTX:
        qb, db, ib = integer_metric(cfg['matrix_b']); _CTX.update(qa=integer_metric(cfg['matrix_a']), qb=(sieve_metric(qb), db), disks=EB.root_disks())
    return _CTX
def v2(n): return (n & -n).bit_length() - 1
def pool_row(cfg, sd, i, full=True):
    """One pool draw: exact norms, primitivity, 2-parts; with full: the two-level threshold counts of its row and 256-cell block."""
    c = _ctx(cfg); (qa, da, ia), (qb, db), disks = c['qa'], c['qb'], c['disks']; A2 = cfg['production_A_squared']
    a, b, aq, bq = regenerate_pair(cfg, sd, i); nf, ng = norms(a, b)
    row = dict(seed=sd, index=i, qa=round(F(aq, da) / A2 * 65536) / 65536., nf=round(16 * log2(nf)) / 16., ng=round(16 * log2(ng)) / 16., v2f=v2(nf), v2g=v2(ng), primitive=coprime_pair(a, b), rejected=False)
    if full: row.update(threshold_counts(cfg, a, b, aq))
    return row
def threshold_counts(cfg, a, b, aq):
    c = _ctx(cfg); (qa, da, ia), (qb, db), disks = c['qa'], c['qb'], c['disks']; A2 = cfg['production_A_squared']
    radius = (F(cfg['production_region_lam']) * A2 - F(aq, da)) * db; b = sieve_coords(b)
    t = sum(qb[0][j] * b[j] for j in range(1, 6)); cc = sum(qb[j][k] * b[j] * b[k] for j in range(1, 6) for k in range(1, 6)) - radius; iv = integer_interval(qb[0][0], t, cc)
    blo = b[0] // 256 * 256; bhi = blo + 255; lo, hi = max(iv.start, blo), min(iv.stop - 1, bhi); rlen = 16 * max(1, (iv.stop - iv.start) // 16); blen = hi - lo + 1
    est = EB.NormBounds(a, disks); box = box_rows(qb, radius, b)
    leaves = [l for l in two_level(est, [(est.row(list(s)), l, h) for s, l, h, _ in box], [q for *_, q in box], blo) if l[0] <= hi and l[1] >= lo and (l[3] is None or 0 in l[3])]
    at = [l for l in leaves if l[0] <= b[0] <= l[1]][0]
    return dict(rejected=isinstance(at[2], dict), evals=est.counts['block_evaluations'], bcells=sum(min(h, bhi) - max(l, blo) + 1 for _, l, h, _ in box if l <= bhi and h >= blo),
                leaves=sum(1 if y - x > 254 else y // 64 - x // 64 + 1 for x, y in ((max(l[0], lo), min(l[1], hi)) for l in leaves)), blen=blen, rlen=rlen)
S0, L1 = 64, 16
def two_level(est, rows, keys, blo):
    """Thresholds of one aligned 256-cell block for the rows of an 8 x 8 x 8 box: level 1 on the hull of the 2 x 2 x 4 sub-box hulls (64-cell nodes,
    halved down to 16 cells), level 2 per sub-box (EB.box_leaves); one merged leaf when the whole block is within the width."""
    W = 16 * OP['width_bits']; ok = lambda bd: bd[0] is not None and bd[0][1] - bd[0][0] <= W; A = min(a for _, a, _ in rows); B = max(b for *_, b in rows); sub = {}
    P = [(x, x + S0 - 1) for x in range(blo + 256 - S0, blo - 1, -S0)]; L = []; l2 = 0
    for k, q in enumerate(keys): sub.setdefault(q, []).append(k)
    H = EB.hull([EB.hull([rows[k][0] for k in ks]) for ks in sub.values()])
    while P:
        lo, hi = P.pop()
        if hi < A or lo > B: continue
        bd = est.block(H, lo, hi)
        if ok(bd): L.append((lo, hi, bd, None))
        elif hi - lo >= L1: m = (lo + hi) // 2; P += [(m + 1, hi), (lo, m)]
        else:
            l2 = 1
            for ks in sub.values(): L += [(x, y, b2, tuple(ks) if k2 is None else (ks[k2],)) for x, y, b2, k2 in EB.box_leaves(est, [rows[k] for k in ks], lo, hi, W)]
    if not l2:
        lo, hi = min(l[2][0][0] for l in L), max(l[2][0][1] for l in L)
        if hi - lo <= W: L = [(blo, blo + 255, [[lo, hi], None], None)]
    return L
BOX = 8, 8, 8
def box_rows(qb, radius, b):
    import itertools
    base = [b[j] // s * s for j, s in zip((2, 3, 4), BOX)]; out = []
    for d in itertools.product(*(range(s) for s in BOX)):
        x = b[1], base[0] + d[0], base[1] + d[1], base[2] + d[2], b[5]; t = sum(qb[0][j] * x[j - 1] for j in range(1, 6))
        c = sum(qb[j][k] * x[j - 1] * x[k - 1] for j in range(1, 6) for k in range(1, 6)) - radius; iv = integer_interval(qb[0][0], t, c)
        if iv: out.append((x, iv.start, iv.stop - 1, (d[0] // 2, d[1] // 2, d[2] // 4)))
    out.sort(key=lambda r: r[0] != tuple(b[1:])); return out
class Geometry:
    """Volumes of the pool: a in {a5 >= 1, Q_a <= theta A^2} (va), b in {c5 >= 1, Q_b <= (lam - qa) A^2} (vb (lam - qa)^3), sieve basis for b."""
    def __init__(self, cfg, theta, lower=COUNT[0], upper=COUNT[1]):
        self.A2 = cfg['production_A_squared']; self.lam, self.kap = float(F(cfg['production_region_lam'])), float(F(cfg['production_region_kap']))
        self.qa, self.da, self.ia = integer_metric(cfg['matrix_a']); qb, self.db, _ = integer_metric(cfg['matrix_b']); self.qb = sieve_metric(qb)
        from common12 import inverse
        self.ib = inverse(self.qb); det_a = F(IC.determinant(self.qa), self.da ** 6); det_b = F(IC.determinant(self.qb), self.db ** 6); v6 = math.pi ** 3 / 6
        self.va = v6 * (float(theta) * self.A2) ** 3 / (2 * sqrt(float(det_a))); self.vb = v6 * self.A2 ** 3 / (2 * sqrt(float(det_b))); self.lower, self.upper = float(lower), float(upper)
        self.rb = [sqrt(float(self.ib[j][j]) * self.A2 * self.db) for j in range(6)]; self.inv55 = float(self.ib[5][5]) * self.db
    def rad(self, q): return self.lam - self.kap * q
    def slice_cells(self, q): return self.vb * self.rad(q) ** 3 * self.upper
    weight = slice_cells
    def lines(self, q): return max(1., floor(sqrt(self.rad(q) * self.A2 * self.inv55)))
    def entries_upper(self, q):
        s = sqrt(self.rad(q)); L4 = 1.
        for j in range(1, 5): L4 *= 2 * floor(self.rb[j] * s) + 3
        return (floor(sqrt(self.rad(q) * self.A2 * self.inv55)) + 2) * L4
def retentions(op=OP):
    """(yield law given B-smooth, survivor law over all integers): Euler-product models without the prime 2, read at n - v_2."""
    kw = dict(prime_min=PRIME_MIN, explicit_bits=EXPLICIT_BITS, step=.25, odd=1)
    return Retention(op['lb'], op['lx2'], max_bits=MAXBITS, tilt=.05, **kw), Retention(MAXBITS, op['lx2'], max_bits=MAXBITS + 1, tilt=.005, **kw)
def row_cost(r, geo, op, cost, D7, pairs_all, c_chain):
    """Instructions per cell of the slice of pool row r, by component."""
    S = geo.slice_cells(r['qa']) * geo.lower / geo.upper; K = geo.lines(r['qa']); lx2 = op['lx2']
    p = dict(build=(C_BUILD_E * geo.entries_upper(r['qa']) + C_BUILD_SLICE) / S + C_BUILD_SEG * 4 / 2 ** 32, setup6=pairs_all * C_SETUP6 / S,
             thr=C_NB_EVAL * r['evals'] / r['bcells'] + C_NB_ROW / r['rlen'] + C_NB_SLICE / S)
    lp = LX0; line = .0
    while lp < P_MID: lp2 = min(lp + .25, P_MID); line += NBAR * log(lp2 / lp) * tab(LINE_TAB, log2(S / K) - .5 * (lp + lp2)); lp = lp2
    lp = P_MID; rec = .0; empty = .0; u7 = h7 = .0
    lg = log2(geo.vb * geo.upper) + r['lrs'] if RUN_T[0] <= r.get('T', 0) <= RUN_T[1] else None; l7 = lg and lg + 1 + log2(geo.lower / geo.upper)
    while lp < lx2:
        lp2 = min(lp + .25, lx2); dsum = log(lp2 / lp); m = .5 * (lp + lp2); log2h = log2(S) - m; dp = NBAR / 2 * (2 ** log2_li(lp2) - 2 ** log2_li(lp))
        if lg is not None and P7 <= lp < X7 and l7 - m >= REC7_TAB[-1][0]: rec += NBAR * dsum * (step7(l7 - m) + EXTRA7); h7 += NBAR * dsum; u7 += dp * (C_SETUP6 * 2 ** -l7 - C_SETUP6 / S)
        elif log2h >= REC_TAB[-1][0]: rec += NBAR * dsum * tab(REC_TAB, log2h)
        else: empty += dp * C_COSET_EMPTY / S
        lp = lp2
    rec += RUN_ALLOW * D7 / S; p['setup6'] = max(.0, p['setup6'] + u7); hits = NBAR * log(lx2 / P_MID) * SIEVE_SIDES; W = op['width_bits']
    no = r['nf'] - r['v2f']; sl = op['slack_f'] - r['v2f']; sf = cost.survivor_upper(no, sl, W); s0 = cost.survivor_upper(no, sl, .25)
    q_g = min(1., Q_F_SAFETY * rho((r['ng'] + a2(r, 1, op)) / op['lb']))
    p.update(line=line * SIEVE_SIDES, rec=rec * SIEVE_SIDES, empty=empty * SIEVE_SIDES, apply=hits * C_APPLY, route=CELL_ROUTE + REC_ROUTE * hits - (REC_ROUTE - REC_ROUTE7) * h7 * SIEVE_SIDES,
             scan=C_SCAN_CELL + C_SCAN_LEAF * r['leaves'] / r['blen'], cof=sf * C_NORMS + s0 * (c_chain[0] + q_g * c_chain[1]))
    return p
def chain(op=OP):
    """(per rechecked survivor, per g-smooth survivor): g batch; then the f side: trial division, two strips."""
    lb = op['lb']; return C_GBATCH, C_TDIV + C_BATCH + C_BATCH2
def evaluate(rows, geo, op=OP, retain=None, cost=None, parts=None):
    """Sums over pool rows (canonical a, inside the region): yield (primitive rows, qa >= Q_CORE) and work (every row)."""
    lb = op['lb']; Sf = op['slack_f']
    if retain is None: retain, cost = retentions(op)
    wsum = ysum = y2 = nf_m = ng_m = nf_r = ng_r = work_sum = work2 = cells = .0; wrows = 0; pairs_all = NBAR / 2 * 2 ** log2_li(op['lx2']); D7 = contam_sup(geo); cc = chain(op)
    for r in rows:
        w = geo.rad(r['qa']) ** 3; wsum += w
        if r['primitive'] and not r['rejected'] and r['qa'] >= Q_CORE:
            rf = retain.at(r['nf'] - r['v2f'], Sf - r['v2f'])['retention']; y = rho((r['nf'] + a2(r, 0, op)) / lb) * rho((r['ng'] + a2(r, 1, op)) / lb) * rf
            ysum += w * y; y2 += w * w * y * y; nf_r += w * y * r['nf']; ng_r += w * y * r['ng']
        nf_m += w * r['nf']; ng_m += w * r['ng']
        if 'lrs' in r:
            pc = row_cost(r, geo, op, cost, D7, pairs_all, cc); contrib = geo.weight(r['qa']) * sum(pc.values()); work_sum += contrib; work2 += contrib * contrib; cells += geo.weight(r['qa']); wrows += 1
            if parts is not None:
                for k, v in pc.items(): parts[k] = parts.get(k, .0) + geo.weight(r['qa']) * v
    return dict(rows=len(rows), ysum=ysum, y2=y2, nfy=nf_r, ngy=ng_r, work=work_sum, work2=work2, w=wsum, nfw=nf_m, ngw=ng_m, cells=cells, wrows=wrows)
SUM_KEYS = 'rows', 'ysum', 'y2', 'nfy', 'ngy', 'work', 'work2', 'w', 'nfw', 'ngw', 'cells', 'wrows'
def add_sums(parts): return {k: sum(p[k] for p in parts) for k in SUM_KEYS}
def pool_totals(s, n_draws, geo):
    """Estimates from the sums of n_draws uniform draws (rows = canonical draws inside the region)."""
    n, ys, y2 = s['rows'], s['ysum'], s['y2']; out = dict(log2_rels=log2(geo.va / n_draws * geo.vb * geo.lower * ys), yield_se_log2=sqrt(max(.0, y2 / n_draws - (ys / n_draws) ** 2) / n_draws) / (ys / n_draws) / LN2,
                                                   nf_rel=s['nfy'] / ys, ng_rel=s['ngy'] / ys, nf_mean=s['nfw'] / s['w'], ng_mean=s['ngw'] / s['w'])
    if s['cells']:
        # work rows: the region rows carrying run and threshold counts (a fixed-rank subsample of the region rows, fraction wrows / rows):
        # work = (region rows / draws) x mean row work; relative variance = (1 - n/N)/n (count of region rows) + cond^2 (mean over the work rows)
        w, w2, nw = s['work'], s['work2'], s['wrows']; nd = n_draws * nw / n; mw = w / nw; cond = sqrt(max(.0, (w2 - nw * mw * mw) / (nw - 1)) / nw) / mw
        out.update(log2_population=log2(geo.va / nd * s['cells']), per_cell=w / s['cells'], log2_est=log2(geo.va / nd * w / UNIT), cond_se=cond, work_se=sqrt(cond ** 2 + (1 - n / n_draws) / n))
    return out
# --- batch tests (exact scaled remainder tree, proof.md section 3; H-FFT)
MR_ROUNDS = 40
FFT_A, FFT_B, FFT_ROOT, FFT_COPY, FFT_LEAF, FFT_XFER, FFT_LOCAL = 8.5, 1200., 6., 3., 2. ** 15, 504., 2 ** 26
def fft_F(N): return FFT_A * N * log2(max(N, 2)) ** 2 + FFT_B + (FFT_XFER * N if N > FFT_LOCAL else .0)
def gbatch_cost(pbits, leaf=9, L=None, c_leaf=FFT_LEAF):
    """Instructions per leaf of one wave of 2^L leaves of `leaf` limbs against a product of pbits bits."""
    SP = ceil(pbits / 64); L = L if L is not None else ceil(log2((SP + 1) / leaf)); G = L + 2; n = lambda j: leaf * 2 ** j
    PT = sum(2 ** (L - j) * (fft_F(2 * n(j - 1)) + FFT_COPY * 2 * n(j - 1)) for j in range(1, L + 1)); RT = sum(2 ** (L - j + 1) * (fft_F(n(j) + G + n(j - 1)) + FFT_COPY * (n(j - 1) + G)) for j in range(1, L + 1))
    return (PT + RT + FFT_ROOT * fft_F(2 * (max(SP, n(L)) + G))) / 2 ** L + c_leaf
def pbits(lo, hi): return 1.000028 * (2. ** hi - (2. ** lo if lo else 0)) / LN2 if lo == 0 else 1.000028 * 2. ** hi / LN2
C_GBATCH = 2. ** 24.07; C_BATCH = 2. ** 22.7; C_BATCH2 = 2. ** 24.11; C_TDIV = 2. ** 23.45      # >= gbatch_cost / tdiv_cost at these sizes ([5t])
def k2(z, r, e=10): v = (z & -z).bit_length() - 1; z >>= v; return gcd(z, pow(r, 2 ** e, z)) << v
def part(z, P):
    m = 1; g = gcd(z, P)
    while g > 1: z //= g; m *= g; g = gcd(z, g)
    return m
def ptree(zs):
    T = [list(zs)]
    while len(T[-1]) > 1: a = T[-1]; T.append([a[i] * a[i + 1] if i + 1 < len(a) else a[i] for i in range(0, len(a), 2)])
    return T
def srt(P, zs, G=None, dx=0):
    """Exact scaled remainder tree: P mod z for every leaf z (Newton quotient at the root offset by dx to exercise the correction)."""
    T = ptree(zs); G = G if G is not None else len(T) + 1; w = lambda Z: 64 * (-(-Z.bit_length() // 64) + G); Z = T[-1][0]; b = w(Z); s = b + P.bit_length() + 8; x = (1 << s) // Z + dx; q = P * x >> s - b; r = (P << b) - q * Z; fix = 0
    while r < 0: q -= 1; r += Z; fix += 1
    while r >= Z: q += 1; r -= Z; fix += 1
    Y = [q % (1 << b)]
    for l in range(len(T) - 2, -1, -1):
        C = T[l]; Yn = []
        for i, z in enumerate(C): v = i // 2; Zv = T[l + 1][v]; sib = C[i ^ 1] if i ^ 1 < len(C) else 1; bv, bc = w(Zv), w(z); Yn.append(Y[v] * sib % (1 << bv) >> bv - bc)
        Y = Yn
    out = []
    for z, y in zip(zs, Y): b = w(z); r = z * y + (1 << b) - 1 >> b; out.append(0 if r == z else r)
    return out, fix
def srt_toy(k, seed, dxs=(-2, 2), G=None):
    import random
    from threshold_model import primes_upto
    ps = [p for p in primes_upto(2 ** 21) if p > 2 ** 20]; P = ptree(ps)[-1][0]; g = random.Random(seed); zs = []
    for _ in range(2 ** k):
        z = g.getrandbits(NORM_BITS - 64) | 1
        for _ in range(g.randrange(3)): z *= g.choice(ps) ** g.randrange(1, 3)
        zs.append(z)
    R = [srt(P, zs, G, dx) for dx in dxs]
    return sum(a != P % z for r in R for a, z in zip(r[0], zs)) + sum(k2(z, a) != part(z, P) for a, z in zip(R[0][0], zs)), max(r[1] for r in R)
def c_mr(nbits, mulc): return MR_ROUNDS * 2 * nbits * mulc
def ecm_instr(fbits, nbits, success=.999):
    import ecm16
    return ecm16.ecm_instr(fbits, nbits, success)
def tdiv_cost():
    """Trial division of both norms by every prime below 2^TD_BITS (5 instructions per limb and prime, 22 per side, 4 loop; divisions out)."""
    from threshold_model import primes_upto
    n = len(primes_upto(2 ** TD_BITS)); L = NORM_BITS // 64; return n * (5 * 2 * L + 2 * 22 + 4) + 12 * 2 * L * NORM_BITS, n
SPLIT_LIMB, SPLIT_SIDE, SPLIT_LOOP = 6, 22, 4; SPLIT_DIV, SPLIT_IDEAL = 2. ** 10, 2. ** 26
def split_cost(lb=OP['lb']):
    """Instructions to factor one kept relation: trial division by every prime of (2^20, 2^BAND_TOP] on both sides, divisions to full power,
    ECM (ecm16 schedule, 99.9% per prime, repeated until found) for the primes of (2^BAND_TOP, B], Miller-Rabin, ideal lookups."""
    L = NORM_BITS // 64; np_ = 1.25506 * 2. ** BAND_TOP / (BAND_TOP * LN2); kb = floor(NORM_BITS / TD_BITS); kt = floor(NORM_BITS / BAND_TOP)
    td = np_ * (SPLIT_LIMB * 2 * L + 2 * SPLIT_SIDE + SPLIT_LOOP) + 2 * kb * SPLIT_DIV
    ecm = 2 * kt * (1.002 * ecm_instr(lb, NORM_BITS) + SPLIT_IDEAL) + 2 * (kt + 1) * c_mr(NORM_BITS, MONT5 * (L / 5.) ** 2)
    return td + ecm + 2 * kb * SPLIT_IDEAL
