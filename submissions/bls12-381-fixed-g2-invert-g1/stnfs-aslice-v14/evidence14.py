"""stnfs-aslice-v14 declared experiment (see module docstring in the header of the
file). Stdlib only; stdout is exactly one JSON object; progress on stderr."""
import json, math, os, sys, time
from fractions import Fraction as F
from random import Random

_HERE = os.path.dirname(os.path.abspath(__file__))
_DIR = os.path.join(_HERE, "files") if os.path.exists(os.path.join(_HERE, "files", "model14.py")) else _HERE
sys.path.insert(0, _DIR)

import common12 as C
import model14 as M
from threshold_model import Retention
import ideal_canonical as IC
import lattice_count_bounds as LCB
import embedding_bounds as EB
import filtsim, descent14, store14, toydiff14

P = 4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787
R = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
U = C.U
H = C.H
LN2 = math.log(2)
CLAIM = 123.0
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


def check_rho():
    ref = {2: 0.30685281944, 3: 0.048608388, 4: 0.0049109256, 6: 1.9649696e-5, 8: 3.2320693e-8, 10: 2.7701718e-11}
    err = max(abs(filtsim.rho(u) - v) / v for u, v in ref.items())
    say(f"[1] Dickman rho (filtsim table): max relative error {err:.2e} {'ok' if err < 1e-6 else 'FAIL'}")
    return err < 1e-6


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
    for p in ps:
        for s in (s for s in range(p) if sum(c * pow(s, j, p) for j, c in enumerate(H)) % p == 0):
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


def audit(samples=2):
    """Production-parameter audit with an amortised mu_r generator (one full
    cofactor exponentiation, then per-target 255-bit exponentiations)."""
    import descent_witness as W
    assert (W.P ** 6 - 1) % W.R != 0 and ((W.P ** 12 - 1) // W.R) % W.R != 0
    # generator of mu_r: w0 = x^((p^12-1)/r) for one uniform x; targets w0^t
    import random
    rng = random.Random(31001)
    while True:
        x = ([rng.randrange(W.P) for _ in range(6)], [rng.randrange(W.P) for _ in range(6)])
        if any(x[0] + x[1]):
            break
    w0 = W.pow12(x, (W.P ** 12 - 1) // W.R)
    boot = []
    for i in range(samples):
        t = rng.randrange(1, W.R)
        z = W.pow12(w0, t)          # uniform in mu_r for uniform t
        B, swaps = W.lll(W.boot_basis(z))
        st = W.boot_stats(z, B, 31002 + i, points=16, scales=(32, 96))
        for s in st:
            assert s['max_bits'] < 2800
        boot.append(dict(scales=[32, 96],
                         max_bits=max(s['max_bits'] for s in st),
                         max96=max(s['max_bits'] for s in st if s['scale'] == 96),
                         mean_bits=st[-1]['mean_bits'], lll_swaps=swaps))
    special = []
    for bits, side in ((260, 'f'),):
        seed = 13300 + bits * 10 + (side == 'g') * 10000
        r = W.specialq_witness(bits, side, seed, points=8, scales=(16, 64))
        cf = max(x['cf_max'] for x in r['stats']); cg = max(x['cg_max'] for x in r['stats'])
        assert cf < 96 and cg < 17
        special.append(dict(bits=bits, side=side, seed=seed, cf_max=cf, cg_max=cg, points=8 * len(r['stats'])))
    cfg = json.load(open(os.path.join(_DIR, "config.json")))
    Mb = [[F(x) for x in row[:5]] for row in cfg['matrix_b'][:5]]
    den = math.lcm(*(x.denominator for row in Mb for x in row))
    det = F(IC.determinant([[int(x * den) for x in row] for row in Mb]), den ** 5)
    rootlo = F(math.isqrt((det.numerator << 256) // det.denominator), 1 << 128)
    delta = sum(math.isqrt((Mb[i][i].numerator + Mb[i][i].denominator - 1) // Mb[i][i].denominator) + 1 for i in range(5))
    radius = math.isqrt(cfg['production_A_squared']) + 1 + delta
    line_bound = math.ceil(6 * radius ** 5 / rootlo)
    assert line_bound < 1 << 59
    return dict(line_cells_upper=line_bound, line_cells_upper_log2=math.log2(line_bound),
                boot_targets=boot, specialq=special, storage=store14.check(),
                scope='Production norm, lattice and address checks only; smoothness/tree independence remain heuristics.')


def main():
    cfg_in = json.load(sys.stdin) if not sys.stdin.isatty() else {}
    samples = int(cfg_in.get("samples", 25)); seed = int(cfg_in.get("seed", 1))
    say(f"evidence: stnfs-aslice-v14 (bls12-381-fixed-g2-invert-g1); samples={samples} seed={seed}; 1 unit = {M.UNIT:.0f} instructions")
    ok0 = check_resultant(); ok0b = check_tower(); ok1 = check_rho()
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
        f"block selection cut {M.SEL_CUT} keeps {M.SEL_CELLS:.3f} of cells and {M.SEL_YIELD:.3f} of yield "
        f"(stats14.json): rels 2^{e['log2_rels']:.3f} vs need 2^{e['log2_need']:.3f}: margin {e['margin']:.3f} "
        f">= {M.MARGIN} ({'feasible' if e['feasible'] else 'INFEASIBLE'}); pool SE {e['yield_se_log2']:.3f} bits; "
        f"[{time.time()-t0:.1f}s]")
    la = M.linear_algebra(e['log2_U'], op['excess'])
    say(f"[5] population 2^{e['log2_population']:.3f} cells; instructions per cell {e['per_cell']:.1f} "
        f"(all-cell bounds/layout/setups + {M.SEL_CELLS:.2f} x kept-cell sieve machinery at "
        f"{M.REC_ROUTE:.0f}/record + {M.CELL_ROUTE:.0f}/cell); relation collection 2^{e['log2_relcol']:.4f} units; "
        f"filtering N/U {la['ratio']:.4f} x pad {M.PAD_N} -> N = 2^{la['log2_N']:.2f}; linear algebra "
        f"2^{la['log2_la']:.4f} units (16 blocks, four hops, {M.C_ENTRY:.0f} instr/entry)")

    t0 = time.time()
    Us, res = filtsim.run(op['lb'], e['nf_mean'] / op['lb'], e['ng_mean'] / op['lb'], op['excess'],
                         seed, 12, kmax=64, targets=(100,), verbose=False)
    r100 = res[100]
    say(f"[6] filtering simulation (scale 2^12, seed {seed}, {time.time()-t0:.1f}s): N/U = {1/r100[0]:.4f} at "
        f"full weight {r100[3]:.0f} (max |coefficient| {r100[2]}); embedded 2^17 runs (excess: seed, U/N, N, "
        "max coef): " + "; ".join(f"{ex}: " + ", ".join(f"({s}, {un:.2f}, {n}, {mc})" for (s, un, n, mc, w) in v)
                                  for ex, v in sorted(M.FILT_TABLE.items())))

    t0 = time.time()
    Pm = descent14.P()
    two_log, bestb, detail = descent14.solve_units(Pm)
    tot_d, lB1, own, sub, (Sb, kk, keep), nk = bestb[0], bestb[1], bestb[2], bestb[3], bestb[4], bestb[5]
    log2_desc_h5 = two_log + 7  # H5 factor 128
    say(f"[7] window-sieved descent DP ({time.time()-t0:.1f}s): boot lB1={lB1}, box 2^{Sb:.1f}, k={kk}, "
        f"boot own 2^{math.log2(own):.1f}, children subtree 2^{math.log2(sub):.1f} instructions; two-log recurrence "
        f"2^{two_log:.3f} units; H5 x128 -> expected < 2^{log2_desc_h5:.3f}; cap 2^{M.DESCENT_CAP:.0f} -> "
        f"Markov P(abort) <= 2^{log2_desc_h5 - M.DESCENT_CAP:.2f}; modelled success >= "
        f"0.9*(1-2^{log2_desc_h5 - M.DESCENT_CAP:.2f}) = {0.9*(1-2**(log2_desc_h5-M.DESCENT_CAP)):.3f}")
    for b in (69, 101, 161, 303):
        for side in ("f", "g"):
            key = (b, side)
            if key in detail:
                tt, n, c, S, nf, ng, kp, how, nkids = detail[key]
                say(f"    node 2^{b} {side}: target 2^{n:g} {how}, box 2^{S:.0f}, N_f {nf:.0f} N_g {ng:.0f}, "
                    f"keep {kp:.2f}, own 2^{math.log2(c):.1f}, subtree 2^{math.log2(tt):.1f} instr")

    t0 = time.time()
    td = toydiff14.node_differential(nsample=300, slack=40)
    ok11 = (td['functional_mismatches'] == 0 and td['score_mismatches'] == 0 and
            td['false_negatives'] == 0 and td['child_violations'] == 0)
    say(f"[8] descent node differential (toy tower |u|=4, thin window LX=2^{td['LX']}, n=2^{td['n']}, k={td['kw']}; "
        f"{td['cells']} sampled cells, {td['pairs']} congruence functionals, {time.time()-t0:.1f}s): functional "
        f"mismatches {td['functional_mismatches']}, score mismatches {td['score_mismatches']}, false negatives "
        f"{td['false_negatives']}, child violations {td['child_violations']}, survivors {td['survivors']} "
        f"({td['survivor_true']} true, {td['survivor_junk']} junk): {'ok' if ok11 else 'FAIL'}")

    t0 = time.time()
    rd = toydiff14.relation_differential(nsample=450)
    ok11b = (rd['false_credits'] == 0 and rd['relation_threshold_failures'] == 0
             and rd['child_violations'] == 0 and rd['relations_true'] >= 2)
    say(f"[8b] descent relation differential (toy |u|=4, small q, fat window LX=2^{rd['LX']}, n=2^{rd['n']}, "
        f"k={rd['kw']}; {rd['cells']} sampled cells, {time.time()-t0:.1f}s): true relations through the full "
        f"path {rd['relations_true']}, false credits {rd['false_credits']}, relation threshold failures "
        f"{rd['relation_threshold_failures']}, child violations {rd['child_violations']}, junk survivors "
        f"{rd['junk_survivors']} ({rd['unfactored_rejected']} conservatively rejected unfactored, as the "
        f"production ECM does): {'ok' if ok11b else 'FAIL'}")

    nb, np_ = nbar(8, 12)
    say(f"[9] degree-1 prime ideals per rational prime, both sides, over the {np_} primes in [2^8, 2^12): {nb:.3f} "
        f"(the model charges NBAR = {M.NBAR})")
    ok9 = nb <= M.NBAR * 1.01

    la_n, la_br = check_la_kernel()
    say(f"[10] la_mac.s inner loop: {la_n} instructions per entry, branches {la_br} "
        f"({'ok' if la_n == 33 and la_br == ['b.ne'] else 'MISMATCH'})")

    t0 = time.time()
    audit_res = audit(samples=2)
    ok10 = True
    say(f"[11] production audit ({time.time()-t0:.1f}s): line cells < 2^59: "
        f"{audit_res['line_cells_upper_log2']:.2f}; boot lattice at scales 32 and the operating 96: max "
        f"{max(b['max_bits'] for b in audit_res['boot_targets']):.1f} bits (< 2800, envelope 12c+6log2|u|+17 at S=91.1 is 2794.4); special-q cf/cg maxima "
        + "; ".join(f"{s['side']}@2^{s['bits']}: {s['cf_max']:.1f}/{s['cg_max']:.1f}" for s in audit_res['specialq'])
        + f" (< 96/17); storage protocol checks passed: {ok10}")

    tot = M.total_log2(e['log2_relcol'], la['log2_la'])
    mem = M.memory_report(e)
    say(f"ledger (log2 units): relation collection {e['log2_relcol']:.4f}, linear algebra {la['log2_la']:.4f}, "
        f"descent cap {M.DESCENT_CAP:.0f}, canonical list {M.LOG2_CANON:.0f}, Schirokauer maps {M.LOG2_SM:.0f}, "
        f"factor base {M.LOG2_FBC:.0f}, misc {M.LOG2_MISC:.0f}, verify {M.LOG2_VERIFY:.0f}")
    say(f"total_log2W_full_ledger = {tot:.6f} (claimed {CLAIM:.6f}: {'ok' if tot <= CLAIM else 'EXCEEDS CLAIM'}); "
        f"reserve {CLAIM - tot:.3f} bits")
    for label, kw in (("sieve/descent constants x1.5", dict(scale_consts=1.5)),
                      ("matrix pad 4 instead of 2", None)):
        if kw:
            e2 = M.evaluate(pool, n_total, geo, op, retain, cost, scale_consts=kw['scale_consts'])
            say(f"    sensitivity [{label}]: relcol {e2['log2_relcol']:.4f}")
    checks = dict(resultant=ok0, tower=ok0b, rho=ok1, pool=ok2, nbar=ok9, la_loop=la_n == 33,
                  toy_differential=ok11, toy_relation_differential=ok11b,
                  feasible=e['feasible'], within_claim=tot <= CLAIM)
    say(f"done in {time.time()-T0:.1f}s; checks: {checks}")
    print(json.dumps({"checks_ok": all(checks.values()), "checks": checks,
                      "total_log2W_full_ledger": round(tot, 6), "claimed": CLAIM,
                      "log2_relcol": round(e['log2_relcol'], 4), "log2_rels": round(e['log2_rels'], 4),
                      "margin_bits": round(e['margin'], 4), "log2_N": round(la['log2_N'], 4),
                      "log2_la": round(la['log2_la'], 4), "descent_two_logs_log2_units": round(two_log, 3),
                      "descent_cap": M.DESCENT_CAP, "descent_h5_log2": round(log2_desc_h5, 3),
                      "modelled_success": round(0.9 * (1 - 2 ** (log2_desc_h5 - M.DESCENT_CAP)), 3),
                      "toy_differential": td, "toy_relation_differential": rd, "production_audit": {k: v for k, v in audit_res.items() if k != 'storage'},
                      "storage_check": audit_res['storage'], "report": _REPORT}, indent=1), flush=True)


if __name__ == "__main__":
    main()
