"""Bounded initial audit for the recover-P package. Stdlib only; stdout one JSON object."""
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

CORE_CLAIM = 121.24261     # GT-core ledger bound (R22)
P_CLAIM = "121.243331167243655"        # recover-P score (R22)
CORE_DECLARED = {'ledger': ('<=', 121.24261), 'relcol': ('<=', 120.61756), 'la': ('<=', 118.714572),
                 'log2_rels': ('>=', 64.0155), 'margin': ('>=', 0.3011), 'work_se': ('<=', 0.03856),
                 'two_log': ('<=', 118.0192), 'markov': ('<=', 0.60138), 'log2_EN': ('<=', 17.77),
                 'over': ('<=', -9.53), 'rep': ('<=', -10.7), 'boot_phi': ('<=', 2.49e-03),
                 'success': ('>=', 0.3917)}


def main():
    text = sys.stdin.read()
    req = json.loads(text) if text.strip() else {}
    if req not in ({},):
        raise ValueError("the audit takes an empty object")
    t0 = time.time()
    checks = {}
    op = M.OP
    cfg = json.load(open(os.path.join(ROOT, "config.json")))
    checks["parameters"] = (op['slack_f'] == 293 and op['slack_g'] == 224 and op['lx2'] == 64.0 and op['lb'] == 69.0
                            and op['excess'] == 1.2 and M.PAD_N == 1.2
                            and M.ECM_FIND == 0.999 and cfg['production_theta_region'] == "719/2000"
                            and M.P_MID == 48.0 and M.CANON_WORK_ALLOWANCE == 1.0772 and M.DESCENT_CAP == 118.75286
                            and M.BAND_FBITS == 48.0 and M.BAND_NBITS == (832, 576) and M.BAND_FMAX == (41, 28)
                            and M.BAND_Q == 0.999 and M.BAND_ROUNDS == 2 and M.BW_M == 21 and M.BW_N == 4 and M.C_ENTRY == 67.0
                            and DS.P.LOG2_RELS == 63.925 and abs(DS.P.MAIN_UF * 69 - 744.8) < 1e-9
                            and abs(DS.P.MAIN_UG * 69 - 498.9) < 1e-9 and DS.P.ECM_SUCCESS == 0.999)
    # record route (R3 lines)
    src = open(os.path.join(ROOT, "STORAGE.md")).read().split("emit_hot:")[1].split("emit_end:")[0]
    hot = [l for l in (x.strip() for x in src.splitlines()) if l]
    per = {}
    for l in open(os.path.join(ROOT, "measurements12.txt")):
        if l.startswith("R3 mode"):
            w = l.split(); n = int(w[4]); per[(w[2], int(w[3]), n)] = (int(w[6]) - int(w[5])) / n
    route = per[("p16", 1, 262144)] - 7 - (per[("p16", 3, 262144)] - per[("p16", 2, 262144)])
    page = (per[("p10", 1, 65536)] - per[("p16", 1, 65536)]) / (2 ** -10 - 2 ** -16)
    checks["route_measured"] = len(hot) == 7 and abs(per[("p16", 0, 262144)] - 9) < 1e-3 and route + page * 2 ** -16 <= M.REC_ROUTE
    # LA finishing budget per row
    checks["la_finish_budget"] = M.C_RED >= 2 * M.MONT5 + 40 and 100 + 24 * 12 + M.C_RED <= 12 * M.W_ROW
    checks["la_transfer_budget"] = 24 * 56 + 20 <= 15 * M.W_ROW
    import store16
    checks["store_protocol"] = bool(store16.check()["routing_stable"])
    # relation collection (as [3]-[5])
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
    bsch = [ecm16.schedule(M.BAND_FBITS, n, M.BAND_Q) for n in M.BAND_NBITS]
    checks["band_strip"] = all(x['B1'] == 8192 and x['curves'] == 92 for x in bsch) and e['band_keep'] > 1 - 1e-4
    la = M.linear_algebra(e['log2_U'], op['excess'])
    checks["la_occupancy"] = 12 * 2**56 >= 2 ** la['log2_N']
    # descent (as [7])
    retain16.load_table(EV.load_stats())
    Pm = DS.P()
    dres = DS.evaluate(Pm, DS.PLAN16, DS.BOOT16)
    two_log = dres['two_log_units']; bt = dres['boot']; cap = M.DESCENT_CAP
    markov = 2 ** (two_log - cap)
    log2_EN, over, coll, pmin = EV.repeat_bound(dres)
    p_desc = 1 - 2 * bt['phi'] - markov - 2.0 ** over - 2.0 ** coll
    p_main = 1.0 - 2.0 ** (la['log2_N'] + 2 - 254)
    success = p_main * p_desc
    rows = DS.report(dres)
    mu_min = min(r['mu'] for r in rows); psi_max = max(max(r['psi'] for r in rows), bt['psi'])
    checks["descent_plan_feasible"] = mu_min >= Pm.MIN_RELS and psi_max <= Pm.PSI_MAX and max(r['phi'] for r in rows) <= 0.05 and coll <= -4
    tot = M.total_log2(e['log2_relcol'], la['log2_la'], cap, la['log2_gen'])
    vals = dict(ledger=tot, relcol=e['log2_relcol'], la=la['log2_la'], log2_rels=e['log2_rels'], margin=e['margin'],
                work_se=work_se, two_log=two_log, markov=markov, log2_EN=log2_EN, over=over, rep=coll,
                boot_phi=bt['phi'], success=success)
    dec = {k: (vals[k] <= v if o == '<=' else vals[k] >= v) for k, (o, v) in CORE_DECLARED.items()}
    checks["core_declared_bounds"] = all(dec.values())
    checks["core_within_claim"] = tot <= CORE_CLAIM
    # recover-P total
    with localcontext() as ctx:
        ctx.prec = 80
        two = Decimal(2)
        p2 = lambda x: (Decimal(x) * two.ln()).exp()
        tot_p = Decimal(2001) / Decimal(2000) * p2(str(CORE_CLAIM)) + p2(80) / Decimal(20000)
        p_log2 = tot_p.ln() / two.ln()
        s_p = Decimal(str(CORE_DECLARED['success'][1])) - two ** -48
        checks["p_budget"] = tot_p < p2(P_CLAIM)
        checks["p_success"] = s_p > Decimal("0.39")

    # exact prices; rounding residual < 5e-5
    _mlog = open(os.path.join(ROOT, "measurements12.txt")).read()
    def _q(tag, p):
        for _l in _mlog.splitlines():
            if _l.startswith(tag + " p=%d:" % p):
                _w = _l.replace("->", " ").replace(",", " ").split()
                return (int(_w[4]) - int(_w[3])) / (int(_w[7]) - int(_w[6]))
    _ok = [abs(M.LINE_TAB[i][1] - _q("line", p)) < 1e-9 for i, p in enumerate((257, 4099, 65537, 1048583))]
    _ok += [abs(M.REC_TAB[i][1] - _q("record", p)) < 1e-9 for i, p in enumerate((257, 65537, 8388617, 67108879, 268435459))]
    for _l in _mlog.splitlines():
        if _l.startswith("bench 1/3 reps:"):
            _w = _l.replace(";", " ").replace(":", " ").split()
            _ok.append(abs(M.C_NB_EVAL - ((int(_w[4]) - int(_w[3])) - (int(_w[8]) - int(_w[7]))) / (2 * int(_w[10]))) < 1e-9)
    checks["prices_exact"] = len(_ok) == 10 and all(_ok)
    _smin = min(math.log2(geo.slice_cells(r['qa'])) for r in pool_reg)
    checks["empty_never"] = _smin - op['lx2'] >= M.REC_TAB[-1][0] and e['comp_log2']['empty'] == float('-inf')
    _lmin = min(math.log2(geo.slice_cells(r['qa']) / geo.lines(r['qa'])) for r in pool_reg)
    checks["line_measured"] = _lmin - M.P_MID >= 5.82   # production line-h stays in the measured tab range
    def _sh(_k):
        _v = e['comp_log2'][_k]
        return (2 ** (_v - tot)) * M.CANON_WORK_ALLOWANCE if _v > -1e300 else 0.0
    residual = ((0.005 / 7.0) * _sh('apply') + (0.0005 / 7.003) * _sh('recover') + (0.0005 / 12.325) * _sh('build')
                + (0.0005 / 0.625) * _sh('zero') + (0.005 / 367.76) * (_sh('cof') + 2 ** (la['log2_la'] - tot) + 2 ** (la['log2_gen'] - tot))
                + 0.0003 * _sh('setup6') + 1e-12)
    checks["reserve_residual"] = residual < 0.00005
    ok = all(checks.values())
    print(json.dumps({
        "checks_ok": ok, "checks": checks,
        "core": {"theta_region": str(th), "a_radius": rmin, "canonical_rows": ncan,
                 "count_products": [float(corr["safe_yield_factor"]), float(corr["safe_work_population_factor"])],
                 "main_ecm": [{k: x[k] for k in ('fbits', 'nbits', 'B1', 'B2', 'wheel', 'curves', 'log2_instructions')} for x in sch],
                 "band_ecm": [{k: x[k] for k in ('fbits', 'nbits', 'B1', 'B2', 'wheel', 'curves', 'log2_instructions')} for x in bsch],
                 "c_band_log2": e['c_band_log2'], "band_keep": e['band_keep'],
                 "route_per_record": route, "route_page": page, "c_red": M.C_RED, "c_surv_log2": e['c_surv_log2'], "per_cell": e['per_cell'], "log2_population": e['log2_population'],
                 "rels_log2": e['log2_rels'], "need_log2": e['log2_need'], "margin": e['margin'], "yield_se_log2": e['yield_se_log2'],
                 "nf_rel": e['nf_rel'], "ng_rel": e['ng_rel'], "work_se": work_se, "allowance": M.CANON_WORK_ALLOWANCE,
                 "relcol_log2": e['log2_relcol'], "la_log2": la['log2_la'], "gen_log2": la['log2_gen'], "log2_N": la['log2_N'],
                 "two_log_log2": two_log, "descent_cap": cap, "markov": markov, "boot_phi": bt['phi'], "rounding_residual": residual,
                 "log2_EN": log2_EN, "over": over, "rep": coll, "success": success, "ledger_log2": tot,
                 "declared": {k: [vals[k], o, v, dec[k]] for k, (o, v) in CORE_DECLARED.items()}},
        "recover_p": {"log2_upper_model": str(p_log2), "declared": P_CLAIM, "success_lower_model": str(s_p)},
        "seconds": round(time.time() - t0, 1),
        "scope": "Changed quantities + recover-P total (H1-H5)."}))
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
