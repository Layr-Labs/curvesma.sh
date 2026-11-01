"""Bounded initial audit for the ordinary-G1 transfer (stdlib only; stdout is one JSON object).

Recomputes, from the shipped model code and data, every quantity changed by this revision:
the theta = 207/500 relation-collection model with the 99.9%-per-prime main-phase ECM
schedule and slacks (300, 228), its incidence-inclusive work standard error and charged
allowance, the linear algebra at the unit/row-tail entry price, the success bound and the full ledger; with
{"descent": true} also the complete descent recurrence at the new missing-log inputs and the
node-counter and repeat terms (the fast run uses their declared bounds).  {"unit_entries": [k, seed]}: H4u.  It then checks the
G1 unit conversion (10 x 101/100 x core + 2^80/2000) and the finite (b, x) -> (h^b, h^x)
bijections.  It does NOT rerun the unchanged inherited certificates (pool regeneration, tower,
fields, alpha, retention toys, ECM toys, filtering at 2^12); `evidence16.py` does all of that
(about 19 s, near the 20 s runner limit, so it is not declared; run it locally).  Nothing here executes a production attack.
"""
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

CORE_CLAIM = 121.6287      # upper bound on the complete GT-core ledger, C = 20,000 units
G1_CLAIM = "124.9650"       # declared score, C = 2,000 units
CORE_DECLARED = {'ledger': ('<=', 121.6287), 'relcol': ('<=', 121.1926), 'la': ('<=', 118.7428),
                 'log2_rels': ('>=', 64.1310), 'margin': ('>=', 0.3012), 'work_se': ('<=', 0.03000),
                 'two_log': ('<=', 117.9122), 'markov': ('<=', 0.60589), 'log2_EN': ('<=', 17.60),
                 'over': ('<=', -9.70), 'rep': ('<=', -10.70), 'boot_phi': ('<=', 1.7864e-04),
                 'success': ('>=', 0.391949)}


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
    if set(req) == {"unit_entries"}:
        print(json.dumps(unit_entries(*map(int, req["unit_entries"])))); return
    if req not in ({}, {"descent": True}, {"descent": False}):
        raise ValueError('input: {} (fast), {"descent": true} (also reruns the descent) or {"unit_entries": [k, seed]}')
    full = bool(req.get("descent"))
    t0 = time.time()
    checks = {}
    op = M.OP
    cfg = json.load(open(os.path.join(ROOT, "config.json")))
    checks["parameters"] = (op['slack_f'] == 300 and op['slack_g'] == 228 and op['lx2'] == 63.0 and op['lb'] == 69.0
                            and M.ECM_FIND == 0.999 and cfg['production_theta_region'] == "207/500"
                            and M.CANON_WORK_ALLOWANCE == 1.0600 and M.CHUNK_CELLS == 2 ** 30 and M.DESCENT_CAP == 118.635
                            and DS.P.LOG2_RELS == 64.131 and abs(DS.P.MAIN_UF * 69 - 747.5) < 1e-9
                            and abs(DS.P.MAIN_UG * 69 - 498.9) < 1e-9 and DS.P.ECM_SUCCESS == 0.999)
    # ---- record routing: direct chunk scatter at emission (STORAGE.md listing, REVISION.md R2) ----
    src = open(os.path.join(ROOT, "STORAGE.md")).read().split("emit_hot:")[1].split("emit_end:")[0]
    hot = [l.split("//")[0].strip() for l in src.splitlines()]
    hot = [l for l in hot if l]
    rv = open(os.path.join(ROOT, "STORAGE.md")).read().split("recv_loop:")[1].split("recv_end:")[0]
    rv = [l.strip() for l in rv.splitlines() if l.strip()]
    checks["route_hot_path"] = len(hot) == 7 and len(rv) == 13 and rv[-1].startswith("b.") and M.REC_ROUTE >= 7 + len(rv) / 4 + 0.01
    lp = EV.la_entry_price()   # unit/general group loops and combined row tail in la_mac.s (STORAGE.md, REVISION.md R3)
    checks["la_entry_price"] = bool(lp['ok']) and lp['tail'] == 2090 and M.C_ENTRY == 44 and (M.BW_M, M.BW_N) == (18, 4)
    import store16
    checks["store_protocol"] = bool(store16.check()["routing_stable"])
    # ---- relation collection on the theta = 207/500 uniform sub-region (as evidence16 [3]-[5]) ----
    hdr, pool = M.load_pool(open(os.path.join(ROOT, "pool12-1.b85")).read() + open(os.path.join(ROOT, "pool12-2.b85")).read())
    n_total = sum(hdr['samples'].values())
    th_pool = F(cfg['production_theta']); th = F(cfg['production_theta_region'])
    pool_reg = [r for r in pool if r['qa'] <= float(th)]
    n_eff = n_total * float((th / th_pool) ** 3)
    rmin = math.isqrt(int(th * cfg['production_A_squared']))
    corr = LCB.correction(cfg, rmin)
    geo = M.Geometry(cfg, th, F(corr['safe_yield_factor']), F(corr['safe_work_population_factor']))
    retain = Retention(int(op['lb']), int(op['lx2']), step=.25, max_bits=850, tilt=.05)
    cost = Retention(850, int(op['lx2']), step=.25, max_bits=851, tilt=.005)
    e = M.evaluate(pool_reg, n_eff, geo, op, retain, cost, return_work=True)
    wc = e.pop('work_contribs'); ncan = len(wc); mw = sum(wc) / ncan
    sw = math.sqrt(sum((x - mw) ** 2 for x in wc) / (ncan - 1))
    cond_se = sw / math.sqrt(ncan) / mw
    pcan = ncan / n_total
    work_se = math.sqrt(cond_se ** 2 + (1 - pcan) / ncan)
    checks["work_allowance"] = M.CANON_WORK_ALLOWANCE >= 1 + M.WORK_SE_K * work_se
    checks["feasible_margin"] = bool(e['feasible'])
    checks["rels_input"] = e['log2_rels'] >= DS.P.LOG2_RELS
    checks["mean_norm_inputs"] = e['nf_rel'] >= 69 * DS.P.MAIN_UF and e['ng_rel'] >= 69 * DS.P.MAIN_UG
    checks["count_products_in_range"] = 0.97 <= float(corr['safe_yield_factor']) <= 1 <= float(corr['safe_work_population_factor']) <= 1.03
    import ecm16
    sch = [ecm16.schedule(op['lb'], s + 8, M.ECM_FIND) for s in (op['slack_f'], op['slack_g'])]
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
    vals = dict(ledger=tot, relcol=e['log2_relcol'], la=la['log2_la'], log2_rels=e['log2_rels'], margin=e['margin'],
                work_se=work_se, two_log=two_log, markov=markov, log2_EN=log2_EN, over=over, rep=coll,
                boot_phi=bt['phi'], success=success)
    dec = {k: (vals[k] <= v if o == '<=' else vals[k] >= v) for k, (o, v) in CORE_DECLARED.items()}
    checks["core_declared_bounds"] = all(dec.values())
    checks["core_within_claim"] = tot <= CORE_CLAIM
    # ---- ordinary-G1 transfer arithmetic (proof.md) ----
    with localcontext() as ctx:
        ctx.prec = 80
        two = Decimal(2)
        p2 = lambda x: (Decimal(x) * two.ln()).exp()
        core = Decimal(101) / Decimal(100) * p2(str(CORE_CLAIM))
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
                 "nf_rel": e['nf_rel'], "ng_rel": e['ng_rel'], "work_se": work_se, "allowance": M.CANON_WORK_ALLOWANCE,
                 "relcol_log2": e['log2_relcol'], "la_log2": la['log2_la'], "gen_log2": la['log2_gen'], "log2_N": la['log2_N'],
                 "two_log_log2": two_log, "descent_cap": cap, "markov": markov, "boot_phi": bt['phi'],
                 "log2_EN": log2_EN, "over": over, "rep": coll, "success": success, "ledger_log2": tot,
                 "declared": {k: [vals[k], o, v, dec[k]] for k, (o, v) in CORE_DECLARED.items()}},
        "g1": {"log2_upper_model": str(g1_log2), "declared": G1_CLAIM, "success_lower_model": str(s_g1)},
        "descent": "recomputed" if full else "declared bounds used; rerun with {\"descent\": true}",
        "seconds": round(time.time() - t0, 1),
        "scope": "Recomputes the changed collection, allowance, linear-algebra, descent, success and ledger quantities and the G1 conversion; inherited certificates are rerun only by evidence16.py. Conditional on H1-H5; not a production run."}))
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
