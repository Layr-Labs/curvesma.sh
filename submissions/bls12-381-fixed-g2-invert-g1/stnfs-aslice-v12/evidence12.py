"""stnfs-aslice-v12 evidence (bls12-381-fixed-g2-invert-g1): one stdlib-only sandbox run; stdin {"samples", "seed"}; progress on stderr,
one JSON object on stdout.  Parts [0] resultant, [0b] tower = F_{p^12}, [0c] field certificate (O_h = Z[y]), [1] Dickman rho, [2] pool
regeneration from per-index seeds (canonical oracle, Q_a, exact norms, primitivity, certified norm-bound block statistics and bounds),
[3] geometry and lattice-count factors, [4] retention/yield vs excess x U, [5] cost model and ledger, [6] filtering simulation,
[7] descent DP, [8] degree-1 ideals per prime, [9] la_mac.s loop count.  proof.md references the parts as [k]."""
import json, math, os, sys, time
from fractions import Fraction as F
from random import Random
_HERE = os.path.dirname(os.path.abspath(__file__))
_DIR = os.path.join(_HERE, "files") if os.path.exists(os.path.join(_HERE, "files", "model12.py")) else _HERE
sys.path.insert(0, _DIR)
import common12 as C                       # noqa: E402
import model12 as M                        # noqa: E402
from threshold_model import Retention      # noqa: E402
import ideal_canonical as IC               # noqa: E402
import lattice_count_bounds as LCB         # noqa: E402
import embedding_bounds as EB              # noqa: E402
import filtsim, descent12               # noqa: E402

P = 4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787
R = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
U = C.U
H = C.H
LN2 = math.log(2)
CLAIM = 122.343600   # @@CLAIM@@
_REPORT = []


def say(*args):
    line = " ".join(str(a) for a in args)
    _REPORT.append(line)
    print(line, file=sys.stderr, flush=True)


# ------------------------------------------------------------------ [0] polynomials, tower --------------------------
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
    """f is even in x and g = x^2 + |u| y is monic in x, so Res_x(f, g) = F0(y)^2 with F0 = f(x^2 -> -|u| y);
    Res_y(h, F0) = N_{K_h}(F0(y)) must vanish mod p for f and g to share a root in F_{p^12}."""
    F0 = [0] * 6
    for d, c in C.TERMS:
        term = (c + [0] * 6)[:6]
        for _ in range(d // 2):
            term = C.mul(term, [0, -U])
        F0 = [x + y for x, y in zip(F0, term)]
    r = C.norm(F0)
    ok = r != 0 and r % P == 0
    say(f"[0] Res_y(h, f(x^2 -> -|u| y)) has {r.bit_length()} bits and is {'0' if ok else 'NOT 0'} mod p: {'ok' if ok else 'FAIL'}")
    return ok


def check_tower():
    yp = f6_pow([0, 1], P)
    cur, ok_irr = [0, 1, 0, 0, 0, 0], True
    for k in range(1, 7):                       # cur = y^(p^k) by composition cur(y) -> cur(y^p)
        acc = [0] * 6
        for coef in reversed(cur):
            acc = f6_mulp(acc, yp)
            acc[0] = (acc[0] + coef) % P
        cur = acc
        is_y = cur == [0, 1, 0, 0, 0, 0]
        if (k in (1, 2, 3) and is_y) or (k == 6 and not is_y):
            ok_irr = False
    uy = f6_mulp([(-U) % P], [0, 1])             # u y with u = -|u|
    chi = f6_pow(uy, (P ** 6 - 1) // 2)
    ok_sq = chi == [P - 1, 0, 0, 0, 0, 0]
    say(f"[0b] h irreducible mod p (y^(p^6) = y, y^(p^k) != y for k <= 3): {'ok' if ok_irr else 'FAIL'}; p mod 4 = {P % 4}, h(0) = {H[0]}; "
        f"(u y)^((p^6-1)/2) = -1 in F_p[y]/h: {'ok' if ok_sq else 'FAIL'} -> F_p[y, x]/(h, g) = F_{{p^12}}")
    cert = IC.field_certificate(P)
    ok_c = cert['squarefree'] and cert['irreducible_mod_p'] and cert['maximal_power_basis'] and cert['discriminant'] == 58673
    say(f"[0c] field certificate: Disc(h) = {cert['discriminant']} = {cert['discriminant_factors']}, squarefree {cert['squarefree']}, "
        f"irreducible mod p {cert['irreducible_mod_p']}, O_h = Z[y]: {cert['maximal_power_basis']}: {'ok' if ok_c else 'FAIL'}")
    return ok_irr and ok_sq and ok_c


def check_rho():
    ref = {2: 0.30685281944, 3: 0.048608388, 4: 0.0049109256, 6: 1.9649696e-5, 8: 3.2320693e-8, 10: 2.7701718e-11}
    err = max(abs(C.rho(u) - v) / v for u, v in ref.items())
    say(f"[1] Dickman rho: max relative error vs reference values {err:.2e} {'ok' if err < 1e-6 else 'FAIL'}")
    return err < 1e-6


# ------------------------------------------------------------------ [2] pool regeneration ---------------------------
def uniform_half_ellipsoid(q, inverse_q, radius2, rng, label_counts=None):
    cap = []
    for j in range(6):
        c = F(radius2) * inverse_q[j][j]
        cap.append(math.isqrt(c.numerator // c.denominator))
    while True:
        a = [rng.randrange(-c, c + 1) for c in cap[:-1]] + [rng.randrange(1, cap[-1] + 1)]
        value = C.quadratic(q, a)
        if value <= radius2:
            return a, value


def regenerate(cfg, seed, index):
    """the pool sampler (v11 yield_experiment.sample): per-index seeded RNG; uniform eligible a (a5 > 0, Q_a <= theta A^2) then
    uniform b in its half-slice (b5 >= 1).  Returns a, b, Q_a numerator (den da), Q_b numerator (den db)."""
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
    n_can = n_prim = n_nb = 0; bad = []
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
        # certified norm-bound blocks (nbound12's Python reference): the aligned 64-block of the row containing b0, width 128
        radius = (A2 - F(aq, da)) * db
        t = sum(qb[0][j] * b[j] for j in range(1, 6)); c = sum(qb[j][k] * b[j] * b[k] for j in range(1, 6) for k in range(1, 6)) - radius
        iv = C.integer_interval(qb[0][0], t, c)
        blo = (b[0] // 64) * 64; lo, hi = max(iv.start, blo), min(iv.stop - 1, blo + 63)
        if 16 * max(1, (iv.stop - iv.start) // 16) != row['rlen'] or hi - lo + 1 != row['blen']:
            bad.append(('row/block length', i)); continue
        est = EB.NormBounds(a, disks)
        leaves = est.adaptive_row(a, b[1:], lo, hi, singleton_policy='reject', width=128)
        at = [l for l in leaves if l[0] <= b[0] <= l[1]][0]
        rej = isinstance(at[2], dict)
        if est.counts['block_evaluations'] != row['evals'] or len(leaves) != row['leaves'] or rej != row['rejected']:
            bad.append(('nbound statistics', i, est.counts['block_evaluations'], row['evals'], len(leaves), row['leaves'])); continue
        if not rej:
            (fl, fu), (gl, gu) = at[2]
            if not (fl <= 16 * math.log2(nf) <= fu and gl <= 16 * math.log2(ng) <= gu):
                bad.append(('bound does not contain the exact norm', i)); continue
        n_nb += 1
    say(f"[2] seed-1 indices 0..{samples - 1} regenerated in {time.time() - t0:.1f}s: {n_can} canonical (embedded pool: "
        f"{sum(1 for r in pool if r['seed'] == 1 and r['index'] < samples)} in this range), {n_prim} primitive; canonical flags, Q_a, exact norms, "
        f"primitivity, norm-bound statistics and bounds all match: {'ok' if not bad else 'FAIL ' + str(bad[:3])}; {n_nb} rows with bounds containing the exact norms")
    return not bad, n_can


# ------------------------------------------------------------------ [8] NBAR ------------------------------------------
def poly_roots_count(poly, p):
    """number of distinct roots in F_p of poly (coefficient list, low degree first): deg gcd(x^p - x, poly)"""
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
    # x^p mod poly
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
    """mean over the primes in [2^lo, 2^hi) of the number of degree-1 prime ideals of K_f and K_g: roots s of h mod p, then roots of
    f(x, s) = x^12 - 2 s x^10 + 2 s^3 x^6 + s^5 x^2 + (1 + s - s^2 - s^3) and of g(x, s) = x^2 + |u| s"""
    from threshold_model import primes_upto
    ps = [q for q in primes_upto(2 ** hi_bits) if q >= 2 ** lo_bits]
    tot = 0
    for p in ps:
        for s in (s for s in range(p) if sum(c * pow(s, j, p) for j, c in enumerate(H)) % p == 0):
            tot += poly_roots_count([1 + s - s * s - s ** 3, 0, s ** 5, 0, 0, 0, 2 * s ** 3, 0, 0, 0, -2 * s, 0, 1], p) + poly_roots_count([U * s, 0, 1], p)
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


# ------------------------------------------------------------------ main -------------------------------------------------
def main():
    cfg_in = json.load(sys.stdin) if not sys.stdin.isatty() else {}
    samples = int(cfg_in.get("samples", 100)); seed = int(cfg_in.get("seed", 1))
    t_start = time.time()
    say(f"evidence: stnfs-aslice-v12 (bls12-381-fixed-g2-invert-g1); samples={samples} seed={seed}; 1 unit = {M.UNIT:.0f} instructions")
    ok0 = check_resultant(); ok0b = check_tower(); ok1 = check_rho()
    cfg = json.load(open(os.path.join(_DIR, "config.json")))
    hdr, pool = M.load_pool(open(os.path.join(_DIR, "pool12-1.b85")).read() + open(os.path.join(_DIR, "pool12-2.b85")).read())
    n_total = sum(hdr['samples'].values())
    disks = EB.root_disks()
    ok2, n_can = verify_pool(cfg, pool, samples, disks)
    # [3] geometry
    corr = LCB.correction(cfg)
    geo = M.Geometry(cfg, M.OP['theta'], F(corr['safe_yield_factor']), F(corr['safe_work_population_factor']))
    say(f"[3] exact metrics: eligible a 2^{math.log2(geo.va):.3f} (continuous), largest half-slice 2^{math.log2(geo.vb):.3f} cells; count factors for every radius >= 1024: "
        f"yield x {corr['safe_yield_factor']}, work x {corr['safe_work_population_factor']} (a in [{corr['a_count_ratio'][0]:.4f}, {corr['a_count_ratio'][1]:.4f}], "
        f"b in [{corr['b_count_ratio_uniform'][0]:.4f}, {corr['b_count_ratio_uniform'][1]:.4f}]); pool {len(pool)} canonical rows of {n_total} ({sum(r['primitive'] for r in pool)} primitive, "
        f"{sum(r['rejected'] for r in pool)} rejected points); lines of the largest slice {geo.lines(0):.0f}, suffix-box entries per line 2^{math.log2(geo.suffix_entries(0)):.2f}")
    # [4]+[5] the model at the operating point
    op = M.OP
    retain = Retention(int(op['lb']), int(op['lx2']), step=.25, max_bits=850, tilt=.05)
    cost = Retention(850, int(op['lx2']), step=.25, max_bits=851, tilt=.005)
    e = M.evaluate(pool, n_total, geo, op, retain, cost)
    rf = retain.at(e['nf_mean'], op['slack_f']); rg = retain.at(e['ng_mean'], op['slack_g'])
    say(f"[4] B=2^{op['lb']:.0f}, X2=2^{op['lx2']:.0f}, slacks ({op['slack_f']}, {op['slack_g']}), excess {op['excess']}, width {op['width_bits']} bits, blocks of {op['L0']}, alpha ({op['alpha_f']}, {op['alpha_g']}); "
        f"mean log2 norms {e['nf_mean']:.1f} / {e['ng_mean']:.1f} (u_f {e['uf']:.2f}, u_g {e['ug']:.2f})")
    say(f"    retention at the mean norms: f {rf['retention']:.4f} (rounding budget {rf['rounding_budget_bits']:.2f} bits, omitted mass <= {rf['omitted_mass_limit_bits']:.1f}), g {rg['retention']:.4f}; "
        f"retained relation orbits 2^{e['log2_rels']:.3f} (x {M.ECM_FIND}^(expected primes in (X2, B]), lower count factor, rejected weight {e['rejected_weight']:.4f}) vs excess x U = 2^{e['log2_need']:.3f} "
        f"(U = 2^{e['log2_U']:.3f}): margin {e['margin']:.3f} bits >= {M.MARGIN} ({'feasible' if e['feasible'] else 'INFEASIBLE'}); pool standard error {e['yield_se_log2']:.3f} bits")
    la = M.linear_algebra(e['log2_U'], op['excess'])
    parts = e['parts']
    say(f"[5] population 2^{e['log2_population']:.3f} cells (mean slice 2^{e['mean_slice_log2']:.2f}, largest 2^{e['largest_slice_log2']:.2f}; pairs below X2 2^{e['pairs_log2']:.2f}; "
        f"{e['records_per_cell']:.2f} records per cell); instructions per cell {e['per_cell']:.1f} = layout {parts['build']:.1f} + zero {parts['zero']:.2f} + line {parts['line']:.1f} + "
        f"record {parts['rec']:.1f} + apply {parts['apply']:.1f} + recover {parts['recover']:.1f} + empty cosets {parts['empty']:.2f} + pair setups {parts['setup6']:.2f} + "
        f"scan {parts['scan']:.1f} + norm bounds {parts['thr']:.1f} + cofactorisation {parts['cof']:.1f} (2^{e['c_surv_log2']:.1f} per survivor)")
    say(f"    relation collection 2^{e['log2_relcol']:.3f} units; filtering N/U = {la['ratio']:.4f} (max over seeds, 2^17 runs, excess {op['excess']}, W = {M.W_ROW:.0f}) x pad {M.PAD_N} -> "
        f"N = 2^{la['log2_N']:.2f}; linear algebra 2^{la['log2_la']:.3f} units")
    # [6] filtering simulation at reduced scale
    t = time.time()
    Us, res = filtsim.run(op['lb'], e['uf'], e['ug'], op['excess'], seed, 13, kmax=64, targets=(100, 150), verbose=False)
    r100 = res[100]
    say(f"[6] filtering simulation (production occurrence profile, dense columns propagated, scale 2^13 sampled columns, seed {seed}, {time.time() - t:.1f}s): "
        f"N/U = {1 / r100[0]:.4f} at full row weight {r100[3]:.0f} (max |coefficient| {r100[2]}); embedded 2^17 runs (excess: seed, U/N, N, max coef): "
        + "; ".join(f"{ex}: " + ", ".join(f"({s}, {un:.2f}, {n}, {mc})" for (s, un, n, mc, w) in v) for ex, v in sorted(M.FILT_TABLE.items())))
    # [7] descent
    t = time.time()
    Pm = descent12.P(); Pm.lB = op['lb']; Pm.alpha_f, Pm.alpha_g = op['alpha_f'], op['alpha_g']; Pm.LMAX = 265; Pm.BOOT_NORM = 2725.0
    ld, bestb, detail = descent12.solve_fast(Pm, s_step=2)
    tot_d, lB1, own, sub, cand, nk = bestb
    log2_desc = ld - math.log2(M.UNIT)
    say(f"[7] descent DP ({time.time() - t:.1f}s; ECM multiplication = measured ecm12 mont_mul {M.MONT5} Ir x (L/5)^2): boot B1 = 2^{lB1}, 2^{math.log2(cand):.1f} candidates, "
        f"boot 2^{math.log2(own):.1f} instr, tree 2^{math.log2(sub):.1f} instr; per logarithm 2^{log2_desc:.2f} units, both 2^{log2_desc + 1:.2f}; "
        f"cap 2^{M.DESCENT_CAP:.0f} -> Markov P(abort) <= {2 ** (log2_desc + 1 - M.DESCENT_CAP):.4f}")
    for b in (69, 81, 151, 261):
        for side in ("f", "g"):
            key = (b - ((b - int(math.floor(op['lb']))) % 2), side)
            if key in detail:
                tt, n, c, S, nf, ng, keep, how, nkids = detail[key]
                say(f"    node 2^{key[0]}-2^{key[0] + 2} {side}: target 2^{n:g} via {how}, 2^{S:.0f} candidates, N_f {nf:.0f} N_g {ng:.0f} bits, keep {keep:.2f}, own 2^{math.log2(c):.1f}, "
                    f"children {nkids:.1f}, subtree 2^{math.log2(tt):.1f} instr")
    # [8] NBAR
    nb, np_ = nbar(8, 12)
    say(f"[8] degree-1 prime ideals per rational prime, both sides, over the {np_} primes in [2^8, 2^12): {nb:.3f} (the sieved ideals; model charges NBAR = {M.NBAR} for hits and for the degree-1 part of U)")
    ok8 = nb <= M.NBAR * 1.01
    # [9] LA kernel
    la_n, la_br = check_la_kernel()
    say(f"[9] la_mac.s inner loop: {la_n} instructions per entry, branches {la_br} ({'ok' if la_n == 33 and la_br == ['b.ne'] else 'MISMATCH'})")
    # ledger
    tot = M.total_log2(e['log2_relcol'], la['log2_la'])
    mem = M.memory_report(e)
    say(f"ledger (log2 units): relation collection {e['log2_relcol']:.3f}, linear algebra {la['log2_la']:.3f}, descent cap {M.DESCENT_CAP:.0f}, canonical list {M.LOG2_CANON:.0f}, "
        f"Schirokauer maps {M.LOG2_SM:.0f}, factor base {M.LOG2_FBC:.0f}, misc {M.LOG2_MISC:.0f}, verify {M.LOG2_VERIFY:.0f}")
    say(f"total_log2W_full_ledger = {tot:.6f}  (claimed {CLAIM:.6f}: {'ok' if tot <= CLAIM else 'EXCEEDS CLAIM'})")
    say(f"memory of one a-slice job (bytes, log2): cells {mem['cells_log2']:.2f} ({mem['segments']} segments of <= 2^62), suffix metadata {mem['meta_log2']:.2f}, records {mem['records_log2']:.2f}, "
        f"total {mem['job_log2']:.2f}")
    for label, kw in (("sieve constants x1.5", dict(scale_consts=1.5)), ("matrix pad 4 instead of 2", dict(pad_n=4.0))):
        e2 = M.evaluate(pool, n_total, geo, op, retain, cost, scale_consts=kw.get('scale_consts', 1.0))
        la2 = M.linear_algebra(e2['log2_U'], op['excess'], pad_n=kw.get('pad_n', M.PAD_N))
        say(f"    sensitivity [{label}]: relcol {e2['log2_relcol']:.3f}, LA {la2['log2_la']:.3f}, total {M.total_log2(e2['log2_relcol'], la2['log2_la']):.3f}")
    for sf, sg in ((op['slack_f'] - 4, op['slack_g'] - 4), (op['slack_f'] + 4, op['slack_g'] + 4)):
        e2 = M.evaluate(pool, n_total, geo, dict(op, slack_f=sf, slack_g=sg), retain, cost)
        say(f"    neighbour slacks ({sf}, {sg}): relcol {e2['log2_relcol']:.3f}, margin {e2['margin']:.3f}")
    checks = dict(resultant=ok0, tower=ok0b, rho=ok1, pool=ok2, nbar=ok8, la_loop=la_n == 33, feasible=e['feasible'], within_claim=tot <= CLAIM)
    say(f"done in {time.time() - t_start:.1f}s; checks: {checks}")
    print(json.dumps({"checks_ok": all(checks.values()), "checks": checks, "total_log2W_full_ledger": round(tot, 6), "claimed": CLAIM,
                      "operating_point": e['op'], "log2_population": round(e['log2_population'], 4), "instructions_per_cell": round(e['per_cell'], 2),
                      "log2_relcol": round(e['log2_relcol'], 4), "log2_rels": round(e['log2_rels'], 4), "log2_U": round(e['log2_U'], 4), "margin_bits": round(e['margin'], 4),
                      "log2_N": round(la['log2_N'], 4), "log2_la": round(la['log2_la'], 4), "descent_two_logs_log2_units": round(log2_desc + 1, 3),
                      "regenerated_canonical_rows": n_can, "memory_log2_bytes_per_job": round(mem['job_log2'], 3), "report": _REPORT}, indent=1), flush=True)


if __name__ == "__main__":
    main()
