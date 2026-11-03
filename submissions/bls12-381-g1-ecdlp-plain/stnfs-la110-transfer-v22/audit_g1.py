"""Initial audit for the G1 transfer (stdlib only; one JSON object): collection, SE/allce, LA, prices, recovery,
success, ledger, G1 conversion. {"descent": true} reruns descent; {"unit_entries": [k, seed]}: H4u."""
import json, math, os, sys, time
from decimal import Decimal, localcontext
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "files") if os.path.exists(os.path.join(HERE, "files", "model16.py")) else HERE
sys.path.insert(0, ROOT)

import evidence16 as EV
import model16 as M
import descent16 as DS
import lattice_count_bounds as LCB
import retain16
from threshold_model import Retention

CORE_CLAIM = 121.2637      # upper bound on the complete GT-core ledger, C = 20,000 units
G1_CLAIM = "124.587070069061269"       # declared score, C = 2,000 units (least 1e-12 ceiling of the upper model)
CORE_DECLARED = {'ledger': ('<=', 121.2637), 'relcol': ('<=', 120.9128), 'la': ('<=', 117.7327), 'log2_rels': ('>=', 64.1310),
                 'margin': ('>=', 0.3012), 'work_se': ('<=', 0.03444), 'two_log': ('<=', 117.5873), 'markov': ('<=', 0.60478),
                 'log2_EN': ('<=', 16.8719), 'over': ('<=', -10.42), 'rep': ('<=', -10.70), 'boot_phi': ('<=', 1.780e-03),
                 'success': ('>=', 0.390320)}


def unit_entries(k, seed):
    """H4u (STORAGE.md): unit share of row weight at mean weight 100; filtsim.merge is copied unchanged
    except that each row also carries its integer combination of original relations."""
    import inspect, filtsim as FS
    src = inspect.getsource(FS.merge).replace("  rows = [dict(r) for r in rows]\n", "  rows = [dict(r) for r in rows]; ORIG = [{i: 1} for i in range(len(rows))]\n")
    src = src.replace("totw += pr.get(-1, 0.0)\n", "totw += pr.get(-1, 0.0)\n      O = ORIG[i]\n      for q_, v_ in ORIG[p].items():\n        O[q_] = O.get(q_, 0) - m * v_\n        if not O[q_]: del O[q_]\n")
    src = src.replace("      results[tlist[ti]] = (nrows, ncols, avg, maxc, wgt)", "      results[tlist[ti]] = (nrows, ncols, avg, maxc, wgt); SNAP[tlist[ti]] = [(rows[j], ORIG[j]) for j in range(len(rows)) if alive[j]]")
    ns = dict(FS.__dict__, SNAP={}); exec(src, ns); FS.merge = ns['merge']
    L, us = 69.0, (747.5 / 69, 498.9 / 69)
    FS.run(L, us[0], us[1], 1.3, seed, k, targets=(100,), verbose=False)
    lsmall = max(8.0, math.log2(round(1.3 * 2 ** k) / 25.0)); S1 = S2 = 0.0
    for u in us:
        for i in range(4000):
            lb = lsmall + (0.8 * L - lsmall) * (i + 0.5) / 4000; x = 2.0 ** lb
            p = min(1.0, 2 * FS.kappa(u, lb / L) / x); d = x * (0.8 * L - lsmall) / 4000 / lb
            S2 += p * p * d; S1 += 2 * p / x * d
    wmid = FS.build_system(L, us, 1.3, 1, k)[3]; W = U1 = 0.0
    for row, orig in ns['SNAP'][100]:
        W += len(row) - 1 + row.get(-1, 0.0); n = len(orig); n1 = sum(1 for v in orig.values() if abs(v) == 1)
        U1 += sum(1 for c, v in row.items() if c >= 0 and abs(v) == 1) + max(0.0, n1 * wmid - n * (n - 1) * S2 - n1 * S1)
    return dict(scale=k, seed=seed, rows=len(ns['SNAP'][100]), unit_fraction=U1 / W, charged=M.UNIT_FRAC, ok=U1 / W >= M.UNIT_FRAC)


def main():
    text = sys.stdin.read()
    req = json.loads(text) if text.strip() else {}
    if set(req) == {"source"}:
        print(json.dumps({"section": req["source"], "text": "\n".join(EV._asm(req["source"]))})); return
    if set(req) == {"unit_entries"}:
        print(json.dumps(unit_entries(*map(int, req["unit_entries"])))); return
    if req not in ({}, {"descent": True}, {"descent": False}):
        raise ValueError('input: {} (fast), {"descent": true} (also reruns the descent) or {"unit_entries": [k, seed]}')
    full = bool(req.get("descent"))
    t0 = time.time()
    checks = {}
    op = M.OP
    cfg = json.load(open(os.path.join(ROOT, "config.json")))
    checks["parameters"] = (op['slack_f'] == 297 and op['slack_g'] == 233 and M.P_MIN == 2048 and M.LX0 == 11.0 and M.P_MID == 45.5 and M.P_TDIV == 30.0 and op['lx2'] == 64.0 and M.WIDE_REC == 2.0 and M.WIDE_SETUP == 2.0 and op['lb'] == 69.0
                            and M.ECM_FIND == 0.999 and cfg['production_theta_region'] == "190/500"
                            and M.CANON_WORK_ALLOWANCE == 1.06905
                            and M.CHUNK_CELLS == 2 ** 30 and M.DESCENT_CAP == 118.3128 and DS.BOOT16 == (214, 320, 16, 4)
                            and not hasattr(M, "C_LR_CHECK") and not hasattr(M, "C_SETUP_LINE") and not hasattr(M, "N_LINE_PRIMES")
                            and M.TDIV_PER_PRIME == 512.0 and M.N_TDIV_PRIMES == M._tdiv_prime_count(M.P_TDIV)
                            and M.C_TDIV == 2.0 * M.TDIV_PER_PRIME * M.N_TDIV_PRIMES
                            and DS.P.LOG2_RELS == 63.97 and abs(DS.P.MAIN_UF * 69 - 745.9) < 1e-9 and abs(DS.P.MAIN_UG * 69 - 498.9) < 1e-9
                            and DS.P.ECM_SUCCESS == 0.999)
    # ---- R14: RS prime-count formula shape-check at sieved pi(2^20) ----
    _sv = bytearray(b'\x01') * 2 ** 20; _sv[0:2] = b'\x00\x00'
    for _p in range(2, 1025):
        if _sv[_p]:
            _sv[_p * _p:: _p] = b'\x00' * ((2 ** 20 - 1 - _p * _p) // _p + 1)
    checks["tdiv_prime_bound"] = M._tdiv_prime_count(20.0) >= sum(_sv) == 82025
    # ---- record routing: direct chunk scatter at emission (STORAGE.md listing, REVISION.md R2) ----
    src = open(os.path.join(ROOT, "STORAGE.md")).read().split("emit_hot:")[1].split("emit_end:")[0]
    hot = [l.split("//")[0].strip() for l in src.splitlines()]
    hot = [l for l in hot if l]
    rv = open(os.path.join(ROOT, "STORAGE.md")).read().split("recv_loop:")[1].split("recv_end:")[0]
    rv = [l.strip() for l in rv.splitlines() if l.strip()]
    checks["route_hot_path"] = len(hot) == 7 and len(rv) == 13 and rv[-1].startswith("b.") and M.REC_ROUTE >= 7 + len(rv) / 4 + 0.01
    tail, red9, tail_n, tail_ok = EV.la_row_tail()   # row tail, reductions and +-1 loop in kernels_asm.b85 (STORAGE.md, linear algebra)
    nu_ok = all(max(M.SCALE_NU[a]) >= max(M.SCALE_NU[b]) for a, b in zip(sorted(M.SCALE_NU), sorted(M.SCALE_NU)[1:])) and M.ratio_N_over_U(1.3) * M.PAD_N >= max(max(v) for k, v in M.SCALE_NU.items() if k >= 15)
    checks["la_prices"] = (tail_ok and tail <= M.C_ROW and red9 <= M.C_RED9 and EV.check_la_kernel('la_mac.s', 'LBB9_1:') == (11, ['b.ne'])
                           and EV.check_la_kernel() == (33, ['b.ne']) and M.C_ENTRY >= 4 + 3 + M.UNIT_FRAC_CHARGED * 11 + (1 - M.UNIT_FRAC_CHARGED) * 33
                           and (M.BW_M, M.BW_N) == (19, 4) and M.C_ROW == 2750 and M.C_RED9 == 420 and M.C_ENTRY == 21)
    checks["matrix_pad"] = nu_ok and M.PAD_N == 1.10
    # ---- R11: exact prices re-parsed from measurements12.txt; table shape ----
    import re
    mlog = "\n".join(EV._asm("measurements12.txt"))
    pairs = {(k, int(p)): F(int(i2) - int(i1), int(h2) - int(h1)) for k, p, i1, i2, h1, h2
             in re.findall(r"(line|record) p=(\d+): Ir (\d+) -> (\d+), hits (\d+) -> (\d+)", mlog)}
    lp_ = (257, 4099, 65537, 1048583); rp_ = (257, 65537, 8388617, 67108879, 268435459)
    lq_ = (24.10, 31.11, 88.31, 310.70); rq_ = (44.90, 78.21, 266.97, 938.08, 2231.66)
    su = re.search(r"setup6 p=1048583: Ir (\d+) -> (\d+), units (\d+) -> (\d+)", mlog)
    nb = re.search(r"bench 1/3 reps: (\d+) (\d+) ; benchsetup 1/3: (\d+) (\d+) ; evaluations/rep (\d+)", mlog)
    su = F(int(su.group(2)) - int(su.group(1)), int(su.group(4)) - int(su.group(3))) if su else None
    nb = F((int(nb.group(2)) - int(nb.group(1))) - (int(nb.group(4)) - int(nb.group(3))), 2 * int(nb.group(5))) if nb else None
    checks["exact_prices"] = bool(
        len(pairs) == 9 and su is not None and nb is not None and su <= 406.87 <= M.C_SETUP6 and nb <= M.C_NB_EVAL
        and all(M.LINE_TAB[i][1] == pairs[("line", p)] for i, p in enumerate(lp_))
        and all(M.REC_TAB[i][1] == pairs[("record", p)] for i, p in enumerate(rp_))
        and all(abs(float(pairs[("line", p)]) - q) < 0.01 for p, q in zip(lp_, lq_))
        and all(abs(float(pairs[("record", p)]) - q) < 0.01 for p, q in zip(rp_, rq_)))

    def _sl(t):
        return [(t[i + 1][1] - t[i][1]) / (F(str(t[i + 1][0])) - F(str(t[i][0]))) for i in range(len(t) - 1)]
    sl, sr = _sl(M.LINE_TAB), _sl(M.REC_TAB)
    checks["price_shape"] = (all(x[1] < x[0] for x in zip(sl, sl[1:])) and all(x[1] < x[0] for x in zip(sr, sr[1:]))
                             and all(M.LINE_TAB[i + 1][1] > M.LINE_TAB[i][1] for i in range(3))
                             and all(M.REC_TAB[i + 1][1] > M.REC_TAB[i][1] for i in range(4)))
    import store16
    checks["store_protocol"] = bool(store16.check()["routing_stable"])
    # ---- R15: no rerun; STORAGE.md must not contain the retired lr_loop snippet ----
    checks["no_rerun"] = "lr_loop:" not in open(os.path.join(ROOT, "STORAGE.md")).read()
    # ---- relation collection on the theta = 190/500 uniform sub-region (as evidence16 [3]-[5]) ----
    hdr, pool = M.load_pool(open(os.path.join(ROOT, "pool12-1.b85")).read() + open(os.path.join(ROOT, "pool12-2.b85")).read())
    n_total = sum(hdr['samples'].values())
    th_pool = F(cfg['production_theta']); th = F(cfg['production_theta_region'])
    pool_reg = [r for r in pool if r['qa'] <= float(th)]
    n_eff = n_total * float((th / th_pool) ** 3)
    rmin = math.isqrt(int(th * cfg['production_A_squared']))
    corr = LCB.correction(cfg, rmin)
    geo = M.Geometry(cfg, th, F(corr['safe_yield_factor']), F(corr['safe_work_population_factor']))
    retain = Retention(int(op['lb']), int(op['lx2']), prime_min=M.P_MIN, step=.25, max_bits=850, tilt=.05)
    cost = Retention(850, int(op['lx2']), prime_min=M.P_MIN, step=.25, max_bits=851, tilt=.005)
    e = M.evaluate(pool_reg, n_eff, geo, op, retain, cost, return_work=True)
    wc = e.pop('work_contribs'); ncan = len(wc); mw = sum(wc) / ncan
    sw = math.sqrt(sum((x - mw) ** 2 for x in wc) / (ncan - 1))
    cond_se = sw / math.sqrt(ncan) / mw
    pcan = ncan / n_total
    work_se = math.sqrt(cond_se ** 2 + (1 - pcan) / ncan)
    checks["work_allowance"] = M.CANON_WORK_ALLOWANCE >= 1 + M.WORK_SE_K * work_se
    checks["feasible_margin"] = bool(e['feasible'])
    checks["rels_input"] = e['log2_rels'] >= DS.P.LOG2_RELS
    checks["pool_draws"] = n_total == M.N_POOL   # yield SE uses all original draws (fixed-pool estimator)
    checks["yield_lower_bound_covers_descent"] = e['log2_rels'] + math.log2(1 - M.YIELD_SE_K * e['yield_se_log2'] * math.log(2)) >= DS.P.LOG2_RELS == M.DESCENT_RELS_MIN
    checks["mean_norm_inputs"] = e['nf_rel'] >= 69 * DS.P.MAIN_UF and e['ng_rel'] >= 69 * DS.P.MAIN_UG
    checks["count_products_in_range"] = 0.97 <= float(corr['safe_yield_factor']) <= 1 <= float(corr['safe_work_population_factor']) <= 1.03
    import ecm16
    sch = [ecm16.schedule(op['lb'], s + 8, M.ECM_FIND) for s in (op['slack_f'], op['slack_g'])]
    # ---- R11: shift coverage of production log2h, and sub-half-ulp residual bound ----
    def _wins(a, b):
        w = []; lp = a
        while lp < b:
            lp2 = min(lp + 0.25, b); w.append((lp, lp2)); lp = lp2
        return w
    lh_min = min(math.log2(geo.slice_cells(r['qa']) / geo.lines(r['qa'])) - 0.5 * (lp + lp2)
                 for r in pool_reg for lp, lp2 in _wins(M.LX0, M.P_MID))
    rh_min = min(math.log2(geo.slice_cells(r['qa'])) - 0.5 * (lp + lp2)
                 for r in pool_reg for lp, lp2 in _wins(M.P_MID, 64.0))
    checks["tab_shift_ranges"] = (M.BOUND_SHIFT == 0.01
                                  and M.tab(M.LINE_TAB, 45.0) == M.LINE_TAB[0][1]
                                  and M.tab(M.REC_TAB, 18.91) > float(M.REC_TAB[0][1])
                                  and (M.P_MID - M.LX0) / 0.25 == int((M.P_MID - M.LX0) / 0.25)
                                  and (64.0 - M.P_MID) / 0.25 == int((64.0 - M.P_MID) / 0.25)
                                  and lh_min - M.BOUND_SHIFT >= M.LINE_TAB[-1][0]
                                  and rh_min - M.BOUND_SHIFT >= M.REC_TAB[-2][0])
    _saved = (M.C_APPLY, M.C_RECOVER, M.C_BUILD, M.C_ZERO, ecm16.MONT5)
    M.C_APPLY, M.C_RECOVER, M.C_BUILD, M.C_ZERO = 7.005, 7.0035, 12.3255, 0.6255
    ecm16.MONT5 = 367.765; ecm16.schedule.cache_clear()
    try:
        e2 = M.evaluate(pool_reg, n_eff, geo, op, retain, cost)
    finally:
        M.C_APPLY, M.C_RECOVER, M.C_BUILD, M.C_ZERO = _saved[:4]
        ecm16.MONT5 = _saved[4]; ecm16.schedule.cache_clear()
    rr = 2.0 ** (e2['log2_relcol'] - e['log2_relcol'])
    checks["residual_bound"] = rr <= 1.000029   # R16: the manifest states 1.000029x; the audit now checks the stated multiplier, not just 1.0001
    la = M.linear_algebra(e['log2_U'], op['excess'])
    # ---- descent at the new inputs (as evidence16 [7]); the fast run takes the declared descent bounds ----
    cap = M.DESCENT_CAP
    if full:
        retain16.load_table(EV.load_stats())
        Pm = DS.P()
        dres = DS.evaluate(Pm, DS.PLAN16, DS.BOOT16)
        two_log = dres['two_log_units']; bt = dres['boot']
        markov = 2 ** (two_log - cap)
        log2_EN, over, coll, pmin = EV.repeat_bound(dres)
        rows = DS.report(dres)
        mu_min = min(r['mu'] for r in rows); psi_max = max(max(r['psi'] for r in rows), bt['psi'])
        checks["descent_plan_feasible"] = mu_min >= Pm.MIN_RELS and psi_max <= Pm.PSI_MAX and max(r['phi'] for r in rows) <= 0.05 and coll <= -4
    else:
        D = {k: v for k, (o, v) in CORE_DECLARED.items()}
        two_log, log2_EN, over, coll, markov = D['two_log'], D['log2_EN'], D['over'], D['rep'], D['markov']; bt = {'phi': D['boot_phi']}
    p_desc = 1 - 2 * bt['phi'] - markov - 2.0 ** over - 2.0 ** coll
    p_main = 1.0 - 2.0 ** (la['log2_N'] + 2 - 254)
    success = p_main * p_desc
    tot = M.total_log2(e['log2_relcol'], la['log2_la'], cap, la['log2_gen'])
    checks["residual_ledger"] = 1 + 2.0 ** (e['log2_relcol'] - tot) * (rr - 1) <= 1.000023   # stated complete-ledger multiplier, checked (R16)
    vals = dict(ledger=tot, relcol=e['log2_relcol'], la=la['log2_la'], log2_rels=e['log2_rels'], margin=e['margin'],
                work_se=work_se, two_log=two_log, markov=markov, log2_EN=log2_EN, over=over, rep=coll,
                boot_phi=bt['phi'], success=success)
    dec = {k: (vals[k] <= v if o == '<=' else vals[k] >= v) for k, (o, v) in CORE_DECLARED.items()}
    checks["core_declared_bounds"] = all(dec.values())
    checks["core_within_claim"] = tot <= CORE_CLAIM
    up4 = lambda v: math.ceil(v * 10000) / 10000
    rows_up = M.total_log2(up4(e['log2_relcol']), up4(la['log2_la']), cap, up4(la['log2_gen']))
    checks["rounded_rows_within_claim"] = rows_up <= CORE_CLAIM   # sum of the upper-rounded ledger rows <= stated core bound
    # ---- ordinary-G1 transfer arithmetic (proof.md) ----
    with localcontext() as ctx:
        ctx.prec = 80
        two = Decimal(2)
        p2 = lambda x: (Decimal(x) * two.ln()).exp()
        core = Decimal(1001) / Decimal(1000) * p2(str(CORE_CLAIM))
        g1 = Decimal(10) * core + p2(80) / Decimal(2000)
        g1_log2 = g1.ln() / two.ln()
        s_g1 = (1 - two ** -128) * Decimal(str(CORE_DECLARED['success'][1])) - two ** -48
        checks["g1_budget"] = g1 < p2(G1_CLAIM)
        checks["g1_success"] = s_g1 > Decimal("0.39")
    toy = {}
    for r in (3, 5, 7, 11, 17, 31):
        image = {(b, (b * a) % r) for b in range(1, r) for a in range(r)}
        toy[r] = len(image) == r * (r - 1)
    checks["toy_distribution_bijections"] = all(toy.values())
    ok = all(checks.values())
    print(json.dumps({
        "checks_ok": ok, "checks": checks,
        "core": {"theta_region": str(th), "a_radius": rmin, "canonical_rows": ncan,
                 "count_products": [float(corr["safe_yield_factor"]), float(corr["safe_work_population_factor"])],
                 "main_ecm": [{k: x[k] for k in ('fbits', 'nbits', 'B1', 'B2', 'wheel', 'curves', 'log2_instructions')} for x in sch],
                 "c_surv_log2": e['c_surv_log2'], "per_cell": e['per_cell'], "log2_population": e['log2_population'],
                 "rels_log2": e['log2_rels'], "need_log2": e['log2_need'], "margin": e['margin'], "yield_se_log2": e['yield_se_log2'],
                 "tab_shift_min_log2h": [lh_min, rh_min], "residual_ratio": rr,
                 "nf_rel": e['nf_rel'], "ng_rel": e['ng_rel'], "work_se": work_se, "allowance": M.CANON_WORK_ALLOWANCE,
                 "relcol_log2": e['log2_relcol'], "la_log2": la['log2_la'], "gen_log2": la['log2_gen'], "log2_N": la['log2_N'],
                 "two_log_log2": two_log, "descent_cap": cap, "markov": markov, "boot_phi": bt['phi'],
                 "log2_EN": log2_EN, "over": over, "rep": coll, "success": success, "ledger_log2": tot,
                 "declared": {k: [vals[k], o, v, dec[k]] for k, (o, v) in CORE_DECLARED.items()}},
        "g1": {"log2_upper_model": str(g1_log2), "declared": G1_CLAIM, "success_lower_model": str(s_g1)},
        "descent": "recomputed" if full else "declared bounds used; rerun with {\"descent\": true}",
        "seconds": round(time.time() - t0, 1),
        "scope": "Recomputes the claim; inherited certificates rerun only by evidence16.py. Conditional on H1-H5."}))
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
