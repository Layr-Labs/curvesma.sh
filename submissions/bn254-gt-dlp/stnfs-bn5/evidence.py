import hashlib, json, math, sys, time
from fractions import Fraction as F
from math import log2, sqrt, floor, ceil
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE / "files" if (HERE / "files").is_dir() else HERE; sys.path.insert(0, str(ROOT))
import common12 as C, model16 as M, alpha12 as A, runs16 as RU, ideal_canonical as IC, embedding_bounds as EB, countmc as CM, lattice_count_bounds as LC, ledger16 as L
T0 = time.time(); LINES = []; CHECKS = {}
def say(*a):
    line = " ".join(str(x) for x in a); LINES.append(line); print(line, file=sys.stderr, flush=True)
def ok(name, cond):
    CHECKS[name] = bool(cond) and CHECKS.get(name, True); return "ok" if cond else "FAIL"
def tower():
    U = C.U; P, R = L.P, L.R
    def pm(a, b):
        c = [0] * 11
        for i, x in enumerate(a):
            if x:
                for j, y in enumerate(b): c[i + j] += x * y
        for i in range(10, 5, -1):
            if c[i]:
                for j in range(6): c[i - 6 + j] -= c[i] * C.H[j]
        return [x % P for x in c[:6]]
    def pw(a, e):
        r = [1, 0, 0, 0, 0, 0]
        while e:
            if e & 1: r = pm(r, a)
            a = pm(a, a); e >>= 1
        return r
    c0 = P == 36 * U ** 4 + 36 * U ** 3 + 24 * U ** 2 + 6 * U + 1 == int(L.PROFILE_P, 16) and R == 36 * U ** 4 + 36 * U ** 3 + 18 * U ** 2 + 6 * U + 1 == int(L.PROFILE_R, 16)
    c1 = pow(P, 12, R) == 1 and all(pow(P, k, R) != 1 for k in (1, 2, 3, 4, 6)) and (P ** 4 - P ** 2 + 1) % R == 0
    say(f"[0] p, r from u = {U} equal the protected profile; ord_r(p) = 12, r | Phi_12(p) {ok('params', c0 and c1)}")
    y = [0, 1, 0, 0, 0, 0]; fr = [y]
    for _ in range(6): fr.append(pw(fr[-1], P))
    irr = fr[6] == y and all(fr[k] != y for k in (1, 2, 3)); uw = [(U * x) % P for x in C.W]; gi = pw(uw, (P ** 6 - 1) // 2) == [P - 1, 0, 0, 0, 0, 0]
    acc = [0] * 6
    for k, ck in enumerate(C.FC): acc = [(a + ck * b) % P for a, b in zip(acc, pm(pw([x % P for x in C.W], k), pw(uw, 4 - k)))]
    say(f"[0a] h irreducible mod p (y^(p^6) = y, y^(p^k) != y for k = 1, 2, 3); (u w)^((p^6-1)/2) = -1 so g is irreducible over F_p6; f(x^2 -> u w) = w^4 F(u) = 0 mod p {ok('tower', irr and gi and acc == [0] * 6)}")
    tot = [0] * 6; ac = C.ONE
    for ck in C.H: tot = [t + ck * a for t, a in zip(tot, ac)]; ac = C.mul(ac, C.SIGP)
    aut = tot == [0] * 6 and C.sigma(C.SIGP) == [0, 1, 0, 0, 0, 0]; cu = C.mul(C.CU, C.sigma(C.CU)) == [-1, 0, 0, 0, 0, 0] and C.sigma(C.W) == C.mul(C.mul(C.CU, C.CU), C.W)
    fb = fr[3] == [x % P for x in C.SIGP] and pw(uw, (P ** 3 - 1) // 2) == [x % P for x in C.CU]; chi = pow(P, 3, R)
    say(f"[0b] sigma is an automorphism of order 2, c sigma(c) = -1, sigma(w) = c^2 w (so tau maps f to c^8 f, g to c^2 g); sigma(y) = y^(p^3), (u w)^((p^3-1)/2) = c: tau = Frobenius^3; chi^2 = -1 mod r {ok('tau', aut and cu and fb and (chi * chi + 1) % R == 0)}")
    d = L.disc_h(); ded = L.maximal_order()
    say(f"[0c] disc h = {d} = 2^4 5^3 7^2; Dedekind at 2, 5, 7: Z[y] is the maximal order; N(w) = {C.norm(C.W)}; units of norm 1 {ok('order', d == 98000 and ded and C.norm(C.W) == 5 and all(C.norm(u) == 1 for u in C.UNITS))}")
def pool(st, seed):
    cfg, pd = st['config'], st['pool_d']; qa, da, _ = C.integer_metric(cfg['matrix_a']); QA = [[F(x) for x in r] for r in cfg['matrix_a']]; A2 = cfg['production_A_squared']; th = F(cfg['production_theta_region']); lim = int(th * A2 * da)
    good = True; n = 0
    for sd, (cnt, lst) in sorted(pd['replay'].items())[seed % 2::2]:
        got = [i for i in range(cnt) if IC.canonical_a(C.draw_a(cfg, int(sd), i)[0], QA)[0]]; good &= got == lst; n += cnt
    say(f"[2] canonical flags of the first {n} draws of seed {sorted(pd['replay'])[seed % 2]} equal the stored list {ok('pool', good)}")
    ul = RU.unit_list(); geo = M.Geometry(cfg, F(cfg['production_theta'])); ret, cost = M.retentions(); D7 = M.contam_sup(geo); pa = M.NBAR / 2 * 2 ** M.log2_li(M.OP['lx2']); cc = M.chain(); res = []
    pick = [pd['anchors'][seed % len(pd['anchors'])]]
    for an in pick:
        r = M.pool_row(cfg, an[0], an[1]); a, aq = C.draw_a(cfg, an[0], an[1]); lo, hi = RU.run(a, ul, qa, lim); r.update(T=hi - lo + 1, lrs=round(RU.lrs(a, lo, hi, qa, da, A2, F(cfg['production_region_lam']), 1), 6))
        w = sum(M.row_cost(r, geo, M.OP, cost, D7, pa, cc).values())
        got = [an[0], an[1], round(r['qa'] * 65536), round(r['nf'] * 16), round(r['ng'] * 16), r['v2f'], r['v2g'], int(r['primitive']), r['T'], lo, hi, r['lrs'], r['evals'], r['bcells'], r['leaves'], r['blen'], r['rlen']]
        can = [x[5] >= 1 and C.quadratic(qa, x) <= lim and IC.canonical_a(x, QA)[0] for x in ([a[0] + t] + a[1:] for t in range(lo - 1, hi + 2))]
        res.append((got == an[:17] and math.isclose(w, an[17], rel_tol=1e-6), all(can[1:-1]) and not can[0] and not can[-1], r['T']))
        est = EB.NormBounds(a, M._ctx(cfg)['disks']); b, _ = C.draw_b(cfg, an[0], an[1], aq); c = C.sieve_coords(b); bd = est.block(est.row(c[1:]), c[0], c[0]); nf, ng = C.norms(a, b)
        res[-1] += (bd[0] is not None and bd[0][0] <= C.scaled_log(nf) and C.scaled_log(nf, upper=True) <= bd[0][1] + 1,)
    say(f"[2a] {len(pick)} anchor rows (norms, 2-parts, primitivity, run, threshold counts, per-cell work) recomputed {ok('pool', all(x[0] for x in res))}; [2r] their runs (T = {[x[2] for x in res]}): every a-value canonical, both neighbours outside {ok('runs', all(x[1] for x in res))}; [2b] exact N_f inside the cell bound {ok('bounds', all(x[3] for x in res))}")
    return geo
def rebuild(st):
    """Arithmetic only: emit draft declarations, never claim qualification."""
    o = L.compute(st); d = dict(L.DECL)
    up = lambda x, scale: math.ceil(math.nextafter(x * scale, math.inf)) / scale
    rows = {k: format(up(v, 10**6), '.6f') for k, v in o['rows'].items()}
    cap = math.ceil(float(rows['two_logs'])) + 20
    t7 = o['t7']
    tot = L.totals(rows, cap, dict(log2_N=o['la']['log2_N'], boot_fail=t7['boot_phi'], abort_log2=t7['abort_log2'], repeat_log2=t7['repeat_log2']))
    success = math.floor(math.nextafter(tot['success'] * 10**6, -math.inf)) / 10**6
    if success <= .39: raise ValueError('insufficient conditional success for gate')
    ex = up(tot['exp_run'], 10**9)
    sk = format(success, '.6f'); ek = format(ex, '.9f')
    K, claim = L.gate(sk, ek)
    score = up(claim, 10**9)
    d.update(rows=rows, cap=cap, success=sk, exp_run=ek,
             split_cap=112,
             gate_K=K, claim=format(score, '.9f'))
    d.pop('run_bound', None)
    return dict(status='draft_arithmetic_only', declared=d,
                summary=dict(tot, cap=cap, gate_K=K, claim=score),
                descent=o['t7'])

def main():
    req = json.load(sys.stdin) if not sys.stdin.isatty() else {}
    part = req.get("part", "ledger"); seed = req.get("seed", 1)
    if part == 'bw':
        import bw7probe as BWR
        cfg=int(req.get('config',0));cases=[dict(kind='diagonal',seed=1),dict(kind='diagonal',seed=2),dict(kind='filtered',seed=1),dict(kind='diagonal',seed=1,bad_pin=True)]
        print(json.dumps(BWR.main(cases[cfg])));return
    if part == 'adaptive_stats':
        import adaptive_replay as ADS
        print(json.dumps(ADS.run({},ROOT)));return
    if part == 'adaptive':
        import adaptive_toy_order as ADT
        r=ADT.run(dict(side=req.get('side','f'),seed=req.get('root_seed',0),order=req.get('order','deterministic')));print(json.dumps(dict(experiment_completed=True,tree=r)));return
    if part == 'pgf':
        if set(req)-{'part','seed'}: raise ValueError('pgf input')
        import pgf_replay as PG
        result=PG.run({'grid':True},ROOT)
        print(json.dumps(dict(all_pass=result['all_pass'],pgf=result)))
        return
    if part == 'porttoy':
        if set(req)-{'part','seed','kind','config','parent','child'}: raise ValueError('porttoy input')
        import childrecompute as CR
        result=CR.run(req,ROOT)
        print(json.dumps(result))
        return
    if set(req) - {"part", "seed"} or part not in ("ledger", "pool", "alpha", "filter", "descent", "laws", "rebuild", "all") or type(seed) is not int or not 0 <= seed <= 1000000: raise ValueError("input")
    if part == 'rebuild':
        print(json.dumps(rebuild(C.stats())))
        return
    st = C.stats(); say(f"part={part} seed={seed}; unit {M.UNIT:.0f} instructions"); out = {}
    if part in ("ledger", "all"): tower()
    if part in ("pool", "all"): pool(st, seed)
    out = L.check_all(st, part, seed, say, ok)
    if part == 'ledger':
        import pgf_replay as PG
        pg=PG.run({'grid':False},ROOT);ok('pgf_calibration',pg['all_pass']);out['pgf_calibration']=pg
        import adaptive_replay as AR
        ad=AR.run({},ROOT);ok('adaptive_calibration',ad['all_pass']);out['adaptive_calibration']=ad
    good = all(CHECKS.values()); say(f"checks {sum(CHECKS.values())}/{len(CHECKS)} {'ALL PASS' if good else 'FAILED: ' + ', '.join(k for k, v in CHECKS.items() if not v)}; {time.time() - T0:.1f} s")
    print(json.dumps(dict(checks=CHECKS, all_pass=good, **out, output_sha256=hashlib.sha256("\n".join(LINES[:-1]).encode()).hexdigest())))
if __name__ == "__main__": main()
