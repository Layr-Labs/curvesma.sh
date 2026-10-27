"""stnfs v16 declared experiment (bls12-381-fixed-g2-invert-g1). Stdlib only;
stdout is exactly one JSON object; progress on stderr.  Sections:
 [0]  exact field/tower certificates;  [1] Dickman rho and the semismooth lower sums
      (constant-integrand identity, independent midpoint rule, exact small-integer
      oracle);  [2] canonical-pool regeneration and verification;  [3] exact region
      metrics;  [4] relation-collection yield/retention/selection ledger with the
      chunk-routed record price;  [5] filtering dimension and linear algebra;
      [6] filtering simulation;  [7] the descent: exact ball counts, the shipped plan
      evaluated (mu, survivors, costs, phi/psi, expected work), the cap and the
      success arithmetic;  [8] toy descents (fixed balls vs the model, a complete
      toy tree);  [9] degree-one ideal count;  [10] la_mac.s recount;  [11] production
      audit: certified boot lattice with ball norms, actual special-q lattices with
      ball norms against the model bins, unit mates, chunk protocol checks."""
import json, math, os, sys, time
from decimal import Decimal, localcontext
from fractions import Fraction as F
from random import Random

_HERE = os.path.dirname(os.path.abspath(__file__))
_DIR = os.path.join(_HERE, "files") if os.path.exists(os.path.join(_HERE, "files", "model16.py")) else _HERE
sys.path.insert(0, _DIR)

import common12 as C
import model16 as M
import descent16 as DS
import toy16 as TY
import store16
from threshold_model import Retention
import ideal_canonical as IC
import lattice_count_bounds as LCB
import embedding_bounds as EB
import filtsim

P = 4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787
U = C.U
H = C.H
LN2 = math.log(2)
CLAIM = 123.5
_REPORT = []
T0 = time.time()


def say(*args):
    line = " ".join(str(a) for a in args)
    _REPORT.append(line)
    print(line, file=sys.stderr, flush=True)


def f6_mulp(A, B):
    return [c % P for c in C.mul([c % P for c in A], [c % P for c in B])]


def f6_pow(A, n):
    r, base = [1], [c % P for c in A]
    while n:
        if n & 1:
            r = f6_mulp(r, base)
        base = f6_mulp(base, base)
        n >>= 1
    return (r + [0] * 6)[:6]


def check_resultant():
    F0 = [0] * 6
    for d, c in C.TERMS:
        term = (c + [0] * 6)[:6]
        for _ in range(d // 2):
            term = C.mul(term, [0, -U])
        F0 = [x + y for x, y in zip(F0, term)]
    r = C.norm(F0)
    ok = r != 0 and r % P == 0
    say(f"[0] Res_y(h, f(x^2 -> uy)) has {r.bit_length()} bits and is {'0' if ok else 'NOT 0'} mod p: {'ok' if ok else 'FAIL'}")
    return ok


def check_tower():
    yp = f6_pow([0, 1], P)
    cur, ok_irr = [0, 1, 0, 0, 0, 0], True
    for k in range(1, 7):
        acc = [0] * 6
        for coef in reversed(cur):
            acc = f6_mulp(acc, yp)
            acc[0] = (acc[0] + coef) % P
        cur = acc
        is_y = cur == [0, 1, 0, 0, 0, 0]
        if (k in (1, 2, 3) and is_y) or (k == 6 and not is_y):
            ok_irr = False
    uy = f6_mulp([(-U) % P], [0, 1])
    chi = f6_pow(uy, (P ** 6 - 1) // 2)
    ok_sq = chi == [P - 1, 0, 0, 0, 0, 0]
    say(f"[0b] h irreducible mod p: {'ok' if ok_irr else 'FAIL'}; (u y)^((p^6-1)/2) = -1: {'ok' if ok_sq else 'FAIL'} -> F_p[y,x]/(h,g) = F_p12")
    cert = IC.field_certificate(P)
    ok_c = cert['squarefree'] and cert['irreducible_mod_p'] and cert['maximal_power_basis'] and cert['discriminant'] == 58673
    say(f"[0c] field certificate: Disc(h) = {cert['discriminant']} ({cert['discriminant_factors']}), O_h = Z[y]: {'ok' if ok_c else 'FAIL'}")
    return ok_irr and ok_sq and ok_c


def check_rho_and_semismooth():
    ref = {2: 0.30685281944, 3: 0.048608388, 4: 0.0049109256, 6: 1.9649696e-5, 8: 3.2320693e-8, 10: 2.7701718e-11}
    err = max(abs(DS.rho_series(u) - v) / v for u, v in ref.items())
    # independent integrator: rho(u) = (1/u) int_{u-1}^{u} rho(t) dt by the implicit
    # trapezoid rule in 300-digit decimals (no cancellation), steps 1/200 and 1/400,
    # Richardson-extrapolated; compared at every integer u <= 64
    def trapezoid(n, umax=64):
        with localcontext() as ctx:
            ctx.prec = 300
            h = Decimal(1) / n; r = [Decimal(1)] * (n + 1); win = Decimal(1); vals = {}
            for i in range(n + 1, umax * n + 1):
                win -= (r[-n - 1] + r[-n]) * h / 2
                ri = (win + h / 2 * r[-1]) / (i * h - h / 2)
                win += h / 2 * (r[-1] + ri); r.append(ri)
                if i % n == 0:
                    vals[i // n] = ri
                if len(r) > 3 * n:
                    del r[:n]
            return vals
    lo, hi = trapezoid(200), trapezoid(400)
    deq = max(abs(DS.rho_series(u) / float(hi[u] + (hi[u] - lo[u]) / 3) - 1) for u in range(2, 65))
    ok_rho = err < 1e-6 and deq < 1e-5
    # brackets: rho_lo <= rho <= rho_hi on the grid, monotone
    ok_br = all(DS.rho_lo(u) <= DS.rho_series(u) <= DS.rho_hi(u) + 1e-300 for u in (1.5, 2.3, 7.7, 15.2, 33.3, 41.0))
    # constant-integrand identity: with rho == 1 the one-prime lower sum must be ln(n/base)
    # (checked by evaluating the window rows: sum of weights = ln(n/base))
    rows, h = DS.window_rows(80.0, 69.0, 48, 48)
    wsum = sum(rows[(1, 0)])
    ok_w = abs(wsum - math.log(80 / 69)) < 1e-12
    # independent two-dimensional midpoint rule (no convolution) vs the lower sum, base 5, n 8, m 25, k 2
    base, n, m, N = 5.0, 8.0, 25, 128
    hh = (n - base) / N
    ts = [base + (i + .5) * hh for i in range(N)]
    mid = DS.rho_series(m / base) + sum(DS.rho_series((m - t) / base) * hh / t for t in ts) \
        + sum(DS.rho_series((m - t - u) / base) * hh * hh / (2 * t * u) for t in ts for u in ts)
    lo = DS.ss_lower2(m, base, n, 48, 2, 0, 48, 0.0, 0.0)
    ok_mid = lo <= mid
    # exact small-integer oracle: every integer <= 2^16, B = 32, Y = 128: all primes <= Y,
    # at most j primes in (B, Y]; the asymptotic lower sum must not exceed the finite frequency
    X = 1 << 16
    spf = list(range(X + 1))
    for p in range(2, math.isqrt(X) + 1):
        if spf[p] == p:
            for v in range(p * p, X + 1, p):
                if spf[v] == v:
                    spf[v] = p
    counts = [0] * 4
    for value in range(1, X + 1):
        v = value; large = 0; ok = True
        while v > 1:
            p = spf[v]; v //= p
            if p > 128:
                ok = False; break
            if p > 32:
                large += 1
        if ok:
            for k in range(large, 4):
                counts[k] += 1
    exact = [c / X for c in counts]
    lows = [DS.ss_lower2(16, 5.0, 7.0, 48, k, 0, 48, 0.0, 0.0) for k in range(4)]
    ok_or = all(l <= e for l, e in zip(lows, exact))
    say(f"[1] Dickman rho (series, 200-digit coefficients): max relative error {err:.1e} vs reference values to u=10, "
        f"independent 300-digit delay-equation integrator at every integer u <= 64: max relative error {deq:.1e} {'ok' if ok_rho else 'FAIL'}; table brackets {'ok' if ok_br else 'FAIL'}; one-prime weight sum = ln(n/base) "
        f"{'ok' if ok_w else 'FAIL'}; two-prime lower sum {lo:.5f} <= independent midpoint rule {mid:.5f} "
        f"{'ok' if ok_mid else 'FAIL'}; exact oracle (all integers <= 2^16, B=32, Y=128): finite frequencies "
        f"{[round(x, 4) for x in exact]} vs asymptotic lower sums {[round(x, 4) for x in lows]} "
        f"({'lower sums below the finite frequencies: ok' if ok_or else 'FAIL'})")
    return ok_rho and ok_br and ok_w and ok_mid and ok_or


def uniform_half_ellipsoid(q, inverse_q, radius2, rng):
    cap = []
    for j in range(6):
        c = F(radius2) * inverse_q[j][j]
        cap.append(math.isqrt(c.numerator // c.denominator))
    while True:
        a = [rng.randrange(-c, c + 1) for c in cap[:-1]] + [rng.randrange(1, cap[-1] + 1)]
        if C.quadratic(q, a) <= radius2:
            return a, C.quadratic(q, a)


def regenerate(cfg, seed, index):
    qa, da, ia = C.integer_metric(cfg['matrix_a']); qb, db, ib = C.integer_metric(cfg['matrix_b'])
    A2, theta = cfg['production_A_squared'], F(cfg['production_theta'])
    rng = Random(f"v11-production-ideal-key:{seed}:{index}")
    a, aq = uniform_half_ellipsoid(qa, ia, int(A2 * theta * da), rng)
    b, bq = uniform_half_ellipsoid(qb, ib, (A2 - F(aq, da)) * db, rng)
    return a, b, aq, bq


def coprime_pair(a, b):
    da, db = C.norm(a), C.norm(b)
    d = math.gcd(da, db)
    if d == 1:
        return True
    rows = [C.mul(v, [0] * i + [1]) for v in (a, b) for i in range(6)]
    key = IC.row_hnf(rows, d)
    return all(key[i][i] == 1 for i in range(6))


def verify_pool(cfg, pool, samples, disks):
    qa, da, ia = C.integer_metric(cfg['matrix_a']); qb, db, ib = C.integer_metric(cfg['matrix_b'])
    A2 = cfg['production_A_squared']
    index = {(r['seed'], r['index']): r for r in pool}
    n_can = n_prim = 0; bad = []
    t0 = time.time()
    for i in range(samples):
        a, b, aq, bq = regenerate(cfg, 1, i)
        can, work, witness = IC.canonical_a(a, cfg['matrix_a'], 100000)
        row = index.get((1, i))
        if can != (row is not None):
            bad.append(('canonical flag', i, can)); continue
        if not can:
            continue
        n_can += 1
        qa_frac = F(aq, da) / A2
        if round(qa_frac * 65536) != round(row['qa'] * 65536):
            bad.append(('Q_a', i)); continue
        nf, ng = C.norms(a, b)
        if abs(math.log2(nf) - row['nf']) > 1 / 16 or abs(math.log2(ng) - row['ng']) > 1 / 16:
            bad.append(('norms', i)); continue
        if coprime_pair(a, b) != row['primitive']:
            bad.append(('primitive', i)); continue
        n_prim += row['primitive']
        radius = (A2 - F(aq, da)) * db
        t = sum(qb[0][j] * b[j] for j in range(1, 6))
        c = sum(qb[j][k] * b[j] * b[k] for j in range(1, 6) for k in range(1, 6)) - radius
        iv = C.integer_interval(qb[0][0], t, c)
        blo = (b[0] // 64) * 64; lo, hi = max(iv.start, blo), min(iv.stop - 1, blo + 63)
        if 16 * max(1, (iv.stop - iv.start) // 16) != row['rlen'] or hi - lo + 1 != row['blen']:
            bad.append(('row/block length', i)); continue
        est = EB.NormBounds(a, disks)
        leaves = est.adaptive_row(a, b[1:], lo, hi, singleton_policy='reject', width=128)
        at = [l for l in leaves if l[0] <= b[0] <= l[1]][0]
        rej = isinstance(at[2], dict)
        if est.counts['block_evaluations'] != row['evals'] or len(leaves) != row['leaves'] or rej != row['rejected']:
            bad.append(('nbound statistics', i)); continue
        if not rej:
            (fl, fu), (gl, gu) = at[2]
            if not (fl <= 16 * math.log2(nf) <= fu and gl <= 16 * math.log2(ng) <= gu):
                bad.append(('bound does not contain the exact norm', i)); continue
    say(f"[2] seed-1 indices 0..{samples-1} regenerated in {time.time()-t0:.1f}s: {n_can} canonical, {n_prim} primitive; "
        f"flags, Q_a, exact norms, primitivity, row lengths, norm-bound statistics and containment all "
        f"{'match: ok' if not bad else 'FAIL ' + str(bad[:3])}")
    return not bad, n_can


def poly_roots_count(poly, p):
    poly = [c % p for c in poly]
    while poly and poly[-1] == 0:
        poly.pop()
    if len(poly) <= 1:
        return 0

    def pmod(a, m):
        a = a[:]
        inv = pow(m[-1], -1, p)
        while len(a) >= len(m):
            if a[-1]:
                z = a[-1] * inv % p; s = len(a) - len(m)
                for j, v in enumerate(m):
                    a[s + j] = (a[s + j] - z * v) % p
            a.pop()
            while a and a[-1] == 0:
                a.pop()
        return a

    def pmulmod(a, b, m):
        r = [0] * (len(a) + len(b) - 1)
        for i, x in enumerate(a):
            if x:
                for j, y in enumerate(b):
                    r[i + j] = (r[i + j] + x * y) % p
        return pmod(r, m)

    r, base, e = [1], [0, 1], p
    while e:
        if e & 1:
            r = pmulmod(r, base, poly)
        base = pmulmod(base, base, poly); e >>= 1
    r = (r + [0, 0])[:max(len(r), 2)]; r[1] = (r[1] - 1) % p
    a, b = poly, pmod(r, poly) if len(r) >= len(poly) else [c % p for c in r]
    while b and any(b):
        a, b = b, pmod(a, b)
    return len(a) - 1


def nbar(lo_bits, hi_bits):
    from threshold_model import primes_upto
    ps = [q for q in primes_upto(2 ** hi_bits) if q >= 2 ** lo_bits]
    tot = 0
    rng = Random(1)
    for p in ps:
        for s in sorted(set(TY.roots_modp(H, p, rng))):   # the roots of h mod p
            tot += poly_roots_count([1 + s - s * s - s ** 3, 0, s ** 5, 0, 0, 0, 2 * s ** 3, 0, 0, 0, -2 * s, 0, 1], p) + \
                poly_roots_count([U * s, 0, 1], p)
    return tot / len(ps), len(ps)


def check_la_kernel():
    lines = open(os.path.join(_DIR, "la_mac.s")).read().splitlines()
    on, body = False, []
    for ln in lines:
        s = ln.strip()
        if s.startswith("LBB0_2:"):
            on = True; continue
        if on:
            if s.startswith("LBB") or (s.startswith(";") and "bb." in s):
                break
            if s and not s.startswith(";") and not s.startswith("."):
                body.append(s.split()[0])
    return len(body), [m for m in body if m.startswith("b.") or m == "b"]


def ball_counts_check():
    """exact N_12(R^2) by the integer-point recurrence for R <= 8 against CELLS; the
    continuous brackets V_12 (R -+ sqrt 3)^12 for every shipped radius."""
    import functools

    @functools.lru_cache(None)
    def Nd(d, s):
        if s < 0:
            return 0
        if d == 0:
            return 1
        r = math.isqrt(s)
        return sum(Nd(d - 1, s - x * x) for x in range(-r, r + 1))
    ok = True
    for R in (4, 6, 8):
        n12 = Nd(12, R * R)
        ok &= (n12 - 1) // 2 == DS.CELLS[R]
    V12 = math.pi ** 6 / 720
    for R, half in DS.CELLS.items():
        n = 2 * half + 1
        ok &= V12 * (R - math.sqrt(3)) ** 12 <= n <= V12 * (R + math.sqrt(3)) ** 12
    return ok


def toy_section():
    lb, n, k = 10, 13, 2
    model = lambda nf, ng: (DS.ss_lower2(math.log2(nf), float(lb), float(n), 48, k, 0, 48, 0.0, 0.0) *
                            DS.ss_lower2(math.log2(ng), float(lb), float(n), 48, k, 0, 48, 0.0, 0.0))
    t0 = time.time()
    rows = TY.fixed_balls(1000, 14, 'g', 2, lb, n, k, 2, model)
    fb = [dict(q=r['q'], cells=r['cells'], raw=r['raw'], distinct=r['distinct'], duplicates=r['duplicates'],
               model_mean=round(r['model_mean'], 2)) for r in rows]
    say(f"[8] toy fixed balls (toy tower |u|=4, deg-4 f; 14-bit g-side ideals, R=2 ball of {rows[0]['cells']} cells, "
        f"B=2^{lb}, window <= {k} primes <= 2^{n}; {time.time()-t0:.1f}s): per ideal (raw, distinct, model mean) "
        + "; ".join(f"({r['raw']}, {r['distinct']}, {r['model_mean']})" for r in fb)
        + f"; total distinct {sum(r['distinct'] for r in fb)} vs model x UDUP {sum(r['model_mean'] for r in fb) * DS.P.UDUP:.1f}")
    t0 = time.time()
    plan = lambda bits: (2, n if bits > n else max(lb, bits - 3), 2 if bits > lb + 3 else 1)
    tr = TY.Tree(6001, lb, plan, max_nodes=8, order='children_first')
    rng = TY.Random(6001)
    q, s, t = TY.random_ideal(rng, 14, 'g')
    ok = tr.descend(q, 'g')
    nodes = [dict(bits=nd['bits'], side=nd['side'], depth=nd['depth'], R=nd['R'], cells=nd['cells'], raw=nd['raw'],
                  distinct=nd['distinct'], tried=nd['tried'], success=nd['success']) for nd in tr.nodes]
    say(f"[8b] toy complete descent from a 14-bit g-side ideal to B=2^{lb} ({time.time()-t0:.1f}s): "
        f"{'success' if ok else 'FAILED'} after {len(tr.nodes)} nodes, {tr.cells} cells; nodes (bits, side, depth, "
        f"distinct relations, relations tried): " + "; ".join(f"({nd['bits']}, {nd['side']}, {nd['depth']}, {nd['distinct']}, {nd['tried']})" for nd in nodes))
    return dict(fixed_balls=fb, tree=dict(success=ok, nodes=nodes, cells=tr.cells)), ok and max(nd['depth'] for nd in nodes) >= 1 and all(r['distinct'] > 0 for r in rows)


def audit(rng_seed=31001):
    import descent_witness as W
    assert (W.P ** 6 - 1) % W.R != 0 and ((W.P ** 12 - 1) // W.R) % W.R != 0
    rng = Random(rng_seed)
    while True:
        x = ([rng.randrange(W.P) for _ in range(6)], [rng.randrange(W.P) for _ in range(6)])
        if any(x[0] + x[1]):
            break
    w0 = W.pow12(x, (W.P ** 12 - 1) // W.R)
    tt = rng.randrange(1, W.R)
    z = W.pow12(w0, tt)
    B, swaps = W.lll(W.boot_basis(z))
    bst = W.boot_ball_stats(z, B, rng_seed + 1, points=12, radii=(14, 96))
    boot_bins_max = max(x[1] for x in DS.P.BOOT_BINS)
    sq = []
    Pm = DS.P()
    cfmax = max(b[1] for b in Pm.QBINS); cgmax = max(b[2] for b in Pm.QBINS)
    env_ok = True
    for bits, side, seed in ((260, 'f', 13300), (150, 'g', 13400)):
        r = W.specialq_ball_stats(bits, side, seed, points=12, radii=(16, 96))
        for st_ in r['stats']:
            env_ok &= st_['cf_max'] <= cfmax and st_['cg_max'] <= cgmax
        sq.append(dict(bits=bits, side=side, seed=seed, swaps=r['swaps'], stats=[{k: (round(v, 2) if isinstance(v, float) else v) for k, v in s.items()} for s in r['stats']]))
        if side == 'f':
            units = W.small_units(2)
            um = W.unit_mates(r['basis'], seed + 7, 96, 2, units)
    return dict(boot=dict(swaps=swaps, stats=[{k: (round(v, 2) if isinstance(v, float) else v) for k, v in s.items()} for s in bst],
                          within_model_bins=all(s['cg_max'] <= boot_bins_max for s in bst)),
                specialq=sq, specialq_within_bins=env_ok,
                unit_mates=dict(points=um['points'], units=um['units'], counts=um['counts'], distinct_over_raw=round(um['distinct_over_raw'], 3)),
                storage=store16.check())


def main():
    cfg_in = json.load(sys.stdin) if not sys.stdin.isatty() else {}
    samples = int(cfg_in.get("samples", 16)); seed = int(cfg_in.get("seed", 1))
    say(f"evidence: stnfs v16 (bls12-381-fixed-g2-invert-g1); samples={samples} seed={seed}; 1 unit = {M.UNIT:.0f} instructions")
    ok0 = check_resultant(); ok0b = check_tower(); ok1 = check_rho_and_semismooth()
    cfg = json.load(open(os.path.join(_DIR, "config.json")))
    hdr, pool = M.load_pool(open(os.path.join(_DIR, "pool12-1.b85")).read() + open(os.path.join(_DIR, "pool12-2.b85")).read())
    n_total = sum(hdr['samples'].values())
    disks = EB.root_disks()
    ok2, n_can = verify_pool(cfg, pool, samples, disks)

    corr = LCB.correction(cfg)
    geo = M.Geometry(cfg, M.OP['theta'], F(corr['safe_yield_factor']), F(corr['safe_work_population_factor']))
    say(f"[3] exact metrics: eligible a 2^{math.log2(geo.va):.3f}, largest half-slice 2^{math.log2(geo.vb):.3f} cells; "
        f"count factors yield x {corr['safe_yield_factor']}, work x {corr['safe_work_population_factor']}; "
        f"{len(pool)} canonical rows of {n_total} ({sum(r['primitive'] for r in pool)} primitive)")

    op = M.OP
    retain = Retention(int(op['lb']), int(op['lx2']), step=.25, max_bits=850, tilt=.05)
    cost = Retention(850, int(op['lx2']), step=.25, max_bits=851, tilt=.005)
    t0 = time.time()
    e = M.evaluate(pool, n_total, geo, op, retain, cost)
    rf = retain.at(e['nf_mean'], op['slack_f']); rg = retain.at(e['ng_mean'], op['slack_g'])
    say(f"[4] B=2^{op['lb']:.0f}, X2=2^{op['lx2']:.0f}, slacks ({op['slack_f']}, {op['slack_g']}), excess {op['excess']}, "
        f"width {op['width_bits']}, blocks of {op['L0']}, alpha ({op['alpha_f']}, {op['alpha_g']}); mean norms "
        f"{e['nf_mean']:.1f}/{e['ng_mean']:.1f}; retention at the means {rf['retention']:.4f}/{rg['retention']:.4f}; "
        f"block selection cut {M.SEL_CUT:.0f} keeps {M.SEL_CELLS:.4f} of cells and {M.SEL_YIELD:.4f} of yield "
        f"(stats16.json): rels 2^{e['log2_rels']:.3f} vs need 2^{e['log2_need']:.3f}: margin {e['margin']:.3f} "
        f">= {M.MARGIN} ({'feasible' if e['feasible'] else 'INFEASIBLE'}); pool SE {e['yield_se_log2']:.3f} bits; "
        f"[{time.time()-t0:.1f}s]")
    la = M.linear_algebra(e['log2_U'], op['excess'])
    say(f"[5] population 2^{e['log2_population']:.3f} cells; instructions per cell {e['per_cell']:.1f} "
        f"(all-cell bounds/layout/setups + {M.SEL_CELLS:.4f} x kept-cell sieve machinery at "
        f"{M.REC_ROUTE:.0f}/record + {M.CELL_ROUTE:.0f}/cell routing, {M.C_APPLY:.2f}/{M.C_RECOVER:.3f} apply/recovery); "
        f"relation collection 2^{e['log2_relcol']:.4f} units; filtering N/U {la['ratio']:.4f} x pad {M.PAD_N} -> N = 2^{la['log2_N']:.4f}; "
        f"linear algebra 2^{la['log2_la']:.4f} units (20 blocks, five hops, {M.C_ENTRY:.0f} instr/entry)")

    t0 = time.time()
    Us, res = filtsim.run(op['lb'], e['nf_mean'] / op['lb'], e['ng_mean'] / op['lb'], op['excess'],
                         seed, 12, kmax=64, targets=(100,), verbose=False)
    r100 = res[100]
    say(f"[6] filtering simulation (scale 2^12, seed {seed}, {time.time()-t0:.1f}s): N/U = {1/r100[0]:.4f} at "
        f"full weight {r100[3]:.0f} (max |coefficient| {r100[2]}); embedded 2^17 runs (excess: seed, U/N, N, "
        "max coef): " + "; ".join(f"{ex}: " + ", ".join(f"({s}, {un:.2f}, {n}, {mc})" for (s, un, n, mc, w) in v)
                                  for ex, v in sorted(M.FILT_TABLE.items())))

    t0 = time.time()
    ok_counts = ball_counts_check()
    Pm = DS.P()
    dres = DS.evaluate(Pm, DS.PLAN16, DS.BOOT16)
    two_log = dres['two_log_units']
    bt = dres['boot']
    cap = M.DESCENT_CAP
    markov = 2 ** (two_log - cap)
    p_desc = 1 - 2 * bt['phi'] - markov
    p_main = 0.9
    success = p_main * p_desc
    rows = DS.report(dres)
    mu_min = min(r['mu'] for r in rows); psi_max = max(max(r['psi'] for r in rows), bt['psi'])
    say(f"[7] descent ({time.time()-t0:.1f}s): exact ball counts (R<=8 recomputed; brackets for all radii) "
        f"{'ok' if ok_counts else 'FAIL'}; shipped plan: {len(rows)} (bin, side) nodes, all mu >= {Pm.MIN_RELS:.0f} "
        f"(min {mu_min:.1f}), psi <= {Pm.PSI_MAX} (max {psi_max:.3f}), phi floor {Pm.PHI_FLOOR}; boot: n1={bt['n1']} "
        f"nA={bt['nA']:.1f} kA={bt['kA']} R={bt['R']} (ball 2^{bt['S']:.1f}, {bt['trials_log2']:.1f} bits of trials at "
        f"2^{bt['trial_log2']:.1f} each, pr 2^{bt['pr_log2']:.1f}), boot own 2^{bt['own_log2']:.2f} + subtree 2^{bt['sub_log2']:.2f} "
        f"instructions x retry {bt['retry']:.3f}; one log 2^{bt['one_log_log2']:.3f} instructions; expected two-log work "
        f"2^{two_log:.3f} units; cap 2^{cap:.0f}: Markov P(cap exceeded) <= 2^{two_log-cap:.2f}; P(boot has no relation) "
        f"<= {bt['phi']:.1e} each; modelled success >= {p_main} x (1 - 2 x {bt['phi']:.1e} - 2^{two_log-cap:.2f}) = {success:.4f}")
    for r in rows:
        if r['bin'] in (70, 100, 150, 200, 250, 305):
            say(f"    node ({r['bin']},{r['side']}): n={r['n']} nA={r['nA']} kA={r['kA']} kB={r['kB']} R={r['R']} S={r['S']} {r['how']} "
                f"mu={r['mu']} pr=2^{r['pr_log2']} surv=2^{r['surv_log2']} own 2^{r['own_log2']} W 2^{r['W_log2']} psi={r['psi']} nf={r['nf']} ng={r['ng']} T={r['T']}")

    toy, ok8 = toy_section()

    nb, np_ = nbar(8, 12)
    say(f"[9] degree-1 prime ideals per rational prime, both sides, over the {np_} primes in [2^8, 2^12): {nb:.3f} "
        f"(the model charges NBAR = {M.NBAR})")
    ok9 = nb <= M.NBAR * 1.01

    la_n, la_br = check_la_kernel()
    say(f"[10] la_mac.s inner loop: {la_n} instructions per entry, branches {la_br} "
        f"({'ok' if la_n == 33 and la_br == ['b.ne'] else 'MISMATCH'})")

    t0 = time.time()
    au = audit()
    ok11 = au['boot']['within_model_bins'] and au['specialq_within_bins']
    say(f"[11] production audit ({time.time()-t0:.1f}s): certified boot lattice (Gram det p^12, LLL conditions; "
        f"{au['boot']['swaps']} swaps), ball norms at R=14 and 96: cg max {max(s['cg_max'] for s in au['boot']['stats'])} "
        f"<= model max {max(x[1] for x in DS.P.BOOT_BINS):.1f} {'ok' if au['boot']['within_model_bins'] else 'FAIL'}; "
        f"actual special-q lattices " + "; ".join(f"{s['side']}@2^{s['bits']}: cf max {max(x['cf_max'] for x in s['stats'])}, cg max {max(x['cg_max'] for x in s['stats'])}" for s in au['specialq'])
        + f" within the model bins (cf <= {max(b[1] for b in DS.P.QBINS):.1f}, cg <= {max(b[2] for b in DS.P.QBINS):.1f}) {'ok' if au['specialq_within_bins'] else 'FAIL'}; "
        f"unit mates on the f@2^260 lattice ({au['unit_mates']['points']} ball points x {au['unit_mates']['units']} height-2 units): "
        f"counts {au['unit_mates']['counts']}, distinct/raw {au['unit_mates']['distinct_over_raw']} (model charges {DS.P.UDUP}); "
        f"largest line <= 2^{au['storage']['line_cells_log2_max']} cells (exact-determinant bound); chunk protocol checks passed: {au['storage']['routing_stable']}")

    tot = M.total_log2(e['log2_relcol'], la['log2_la'], cap)
    say(f"ledger (log2 units): relation collection {e['log2_relcol']:.4f}, linear algebra {la['log2_la']:.4f}, "
        f"descent cap {cap:.0f} (expected {two_log:.3f}), canonical list {M.LOG2_CANON:.0f}, Schirokauer maps {M.LOG2_SM:.0f}, "
        f"factor base {M.LOG2_FBC:.0f}, misc {M.LOG2_MISC:.0f}, verify {M.LOG2_VERIFY:.0f}")
    say(f"total_log2W_full_ledger = {tot:.6f} (claimed {CLAIM:.6f}: {'ok' if tot <= CLAIM else 'EXCEEDS CLAIM'}); "
        f"reserve {CLAIM - tot:.3f} bits")
    # every relation-collection constant x1.5 scales the per-cell cost linearly
    rc15 = e['log2_relcol'] + math.log2(1.5)
    say(f"    sensitivity [sieve constants x1.5]: relcol {rc15:.4f}, total {M.total_log2(rc15, la['log2_la'], cap):.4f}")
    checks = dict(resultant=ok0, tower=ok0b, rho_semismooth=ok1, pool=ok2, nbar=ok9, la_loop=la_n == 33,
                  ball_counts=ok_counts, descent_plan_feasible=mu_min >= Pm.MIN_RELS and psi_max <= Pm.PSI_MAX,
                  toy=ok8, audit=ok11, feasible=e['feasible'], within_claim=tot <= CLAIM, success_ok=success >= 0.39)
    say(f"done in {time.time()-T0:.1f}s; checks: {checks}")
    print(json.dumps({"checks_ok": all(checks.values()), "checks": checks,
                      "total_log2W_full_ledger": round(tot, 6), "claimed": CLAIM,
                      "log2_relcol": round(e['log2_relcol'], 4), "log2_rels": round(e['log2_rels'], 4),
                      "margin_bits": round(e['margin'], 4), "log2_N": round(la['log2_N'], 4),
                      "log2_la": round(la['log2_la'], 4), "descent_two_logs_log2_units": round(two_log, 3),
                      "descent_cap": cap, "markov_abort": markov, "modelled_success": round(success, 4),
                      "descent_boot": {k: (round(v, 4) if isinstance(v, float) else v) for k, v in bt.items()},
                      "descent_plan": rows, "toy": toy, "production_audit": {k: v for k, v in au.items() if k != 'storage'},
                      "storage_check": au['storage'], "report": _REPORT}, indent=1), flush=True)


if __name__ == "__main__":
    main()
