#!/usr/bin/env python3
"""Re-derive every number in the claim and check every recorded run.

Checks, all recomputed here rather than quoted:
  * the itemisation sums (fold = sum of its rows, step totals)
  * the step budget S from the single-source constants, as an exact integer
  * W_raw and W_stress as exact Fractions, and the claim-window rules
  * bit-disjointness of the routing fields and the recorded collapse test
  * every recorded check from part_a / sim / toy / routing
  * single-constant discipline: no module restates a walk constant or a price
"""
import json, os, re, sys
from fractions import Fraction
from math import isqrt, log, log2
import constants as K
import ledger as L
import routing

HERE = os.path.dirname(os.path.abspath(__file__))
CERT = []


def chk(name, ok, detail=""):
    CERT.append({"check": name, "result": "PASS" if ok else "FAIL", "detail": detail})
    return ok


def load(fn):
    p = os.path.join(HERE, fn)
    if not os.path.exists(p):
        return None
    return json.load(open(p))


def main():
    r = L.rows()

    # ---- itemisation sums
    parts = ('selector_ab', 'selector_zeta', 'fold_adds', 'fold_cmp', 'hashing',
             'dp_test', 'detector', 'hist_slot', 'loop')
    chk("itemisation_fold_sum", r['fold'] == sum(r[k] for k in parts),
        "fold = %d" % r['fold'])
    chk("itemisation_mul_step", r['mul_step'] == K.C_REF + r['fold'])
    chk("itemisation_sqr_step", r['sqr_step'] == r['cyclo_sqr'] + r['fold'])
    chk("itemisation_cyclo_row", r['cyclo_sqr'] == 24*K.FP_MUL + 122*K.FP_ADD,
        "cyclo = %d" % r['cyclo_sqr'])

    # ---- exact budget and ledgers
    S = K.step_budget()
    s2 = (2*float(K.MU_WALK)*K.N_CLASSES*log(1/(1-float(K.TARGET)))
          / (1 - float(K.FRAC_SQR)**2))
    chk("step_budget_is_least", (S*S >= s2) and ((S-1)*(S-1) < s2), "S = %d" % S)
    Lg = L.ledgers()
    chk("ledger_S_matches_constants", Lg['S'] == S)
    chk("ledger_c_mul", Lg['c_mul'] == Fraction(r['mul_step'], K.C_REF))
    chk("ledger_c_sqr", Lg['c_sqr'] == Fraction(r['sqr_step'], K.C_REF))
    chk("ledger_w_raw", Lg['w_raw'] == S*((1-K.FRAC_SQR)*Lg['c_mul'] + K.FRAC_SQR*Lg['c_sqr']),
        "W_raw = 2^%.5f" % log2(float(Lg['w_raw'])))
    chk("ledger_w_stress", Lg['w_stress'] > Lg['w_raw'],
        "W_stress = 2^%.5f" % log2(float(Lg['w_stress'])))

    # ---- claim window
    claim = float(K.CLAIM)
    raw = log2(float(Lg['w_raw'])); st = log2(float(Lg['w_stress']))
    chk("claim_strictly_below_head", claim < 125.90, "claim %.2f < head 125.90" % claim)
    chk("claim_at_least_raw_plus_0.10", claim >= raw + 0.10, "pad %.5f" % (claim-raw))
    chk("claim_above_stress_plus_0.005", claim >= st + 0.005, "pad %.5f" % (claim-st))
    chk("claim_below_stress_plus_0.06", claim <= st + 0.06, "pad %.5f" % (claim-st))

    # ---- success floor
    success = 1 - pow(2.718281828459045,
                      -S*S*(1-float(K.FRAC_SQR)**2)/(2*float(K.MU_WALK)*K.N_CLASSES))
    chk("budgeted_success_at_mu", abs(success - float(K.TARGET)) < 1e-4,
        "model %.4f" % success)
    chk("success_floor_above_requirement",
        K.TARGET - K.LOSSES > K.SUCCESS,
        "%.4f - %.4f = %.4f > %.2f" % (K.TARGET, K.LOSSES, float(K.TARGET-K.LOSSES), float(K.SUCCESS)))

    # ---- routing
    chk("routing_fields_disjoint", routing.disjoint())
    chk("routing_layout_production_width",
        routing.TABLE_BITS == 14 and routing.DP_BITS == 80 and routing.STEP_BITS == 3)
    ro = load("routing-out.json")
    chk("routing_histogram_present", ro is not None and ro.get('keys', 0) >= 100000)
    if ro:
        chk("routing_no_collapse", ro['collapse'] == 0,
            "table_distinct %d, dp_distinct %d" % (ro['table_distinct'], ro['dp_distinct']))
        chk("routing_table_spread", ro['table_distinct'] > 0.95*(1 << routing.TABLE_BITS))

    # ---- recorded runs
    pa = load("part_a-out.json")
    chk("part_a_present", pa is not None)
    if pa:
        chk("part_a_all_checks", pa['part_a'] and pa['checks_n'] >= 20,
            "%d checks" % pa['checks_n'])
        for n, v in pa['checks'].items():
            chk("part_a_" + n, bool(v))
    sm = load("sim-out.json")
    chk("sim_present", sm is not None)
    if sm:
        n_groups = sum(x['groups'] for x in sm['rows'])
        chk("sim_group_count", n_groups >= 4000, "%d groups" % n_groups)
        chk("sim_pooled_above_requirement", sm['pooled_p'] >= float(K.SUCCESS),
            "pooled %.4f +- %.4f" % (sm['pooled_p'], sm['pooled_se_binomial']))
        chk("sim_pooled_within_noise_of_model",
            abs(sm['pooled_p'] - sm['model_pooled']) <= 3*sm['pooled_se_binomial'],
            "measured %.4f model %.4f se %.4f" % (sm['pooled_p'], sm['model_pooled'],
                                                  sm['pooled_se_binomial']))
        chk("sim_uses_production_constants",
            abs(sm['f'] - float(K.FRAC_SQR)) < 1e-12 and abs(sm['mu'] - float(K.MU_WALK)) < 1e-12)
        chk("sim_injectivity_penalty_is_conservative",
            sm['injectivity_ratio'] <= sm['predicted_ratio'],
            "measured %.4f charged %.4f" % (sm['injectivity_ratio'], sm['predicted_ratio']))
    ty = load("toy-out.json")
    chk("toy_present", ty is not None)
    if ty:
        chk("toy_end_to_end_verified", ty['all_verified'],
            "%d/%d recovered and verified" % (ty['verified'], ty['recovered']))
        chk("toy_secret_match", ty['secret_match'] == ty['recovered'] > 0)
        chk("toy_cap_refusal", ty['cap_probe_refused'])
        chk("toy_no_alignment_failure", ty['aligned_fail'] == 0)
        chk("toy_production_shape", ty['N'] >= 500000 and ty['p'] % 12 == 7)

    # ---- single-constant discipline: numeric literals only, so prose and
    # comments are ignored; every price, the target, the mixture and the claim
    # must come from constants.py.
    import tokenize, io
    forbidden_numbers = {20000, 383, Fraction(392, 1000), Fraction(3, 8),
                         Fraction(12589, 100), 0.392, 0.375, 125.89}
    viol = []
    for fn in sorted(os.listdir(HERE)):
        if not fn.endswith(".py") or fn in ("constants.py", "verify_certificate.py"):
            continue
        src = open(os.path.join(HERE, fn), "rb").read()
        for tok in tokenize.tokenize(io.BytesIO(src).readline):
            if tok.type != tokenize.NUMBER:
                continue
            try:
                val = Fraction(tok.string).limit_denominator(10**6) if "/" in tok.string or "." in tok.string else Fraction(int(tok.string, 0))
            except Exception:
                continue
            if val in forbidden_numbers:
                viol.append("%s:%d %s" % (fn, tok.start[0], tok.string))
    chk("single_constant_import", not viol, "; ".join(viol[:3]) or "no restated constant")

    ok = all(c['result'] == "PASS" for c in CERT)
    cert = {"package": "gt-rho-tight-40pct-v4",
            "certificate": "PASS" if ok else "FAIL",
            "checks": CERT,
            "numbers": {
                "claim": str(K.CLAIM), "claim_float": claim,
                "W_raw_log2": raw, "W_stress_log2": st,
                "S": S, "S_log2": log2(S),
                "c_mul": float(Lg['c_mul']), "c_sqr": float(Lg['c_sqr']),
                "c_raw": float(Lg['c_raw']),
                "fold_instructions": r['fold'], "cyclo_sqr_instructions": r['cyclo_sqr'],
                "raw_pad": claim - raw, "stress_pad": claim - st,
                "budgeted_success": success,
                "success_floor": float(K.TARGET - K.LOSSES),
                "N_classes_log2": log2(K.N_CLASSES),
                "Fp12_multiplies": "the reference operation, C = %d" % K.C_REF,
            }}
    json.dump(cert, open(os.path.join(HERE, "certificate.json"), "w"), indent=1)
    npass = sum(c['result'] == "PASS" for c in CERT)
    print("CERTIFICATE", cert['certificate'], "%d/%d checks" % (npass, len(CERT)))
    for c in CERT:
        if c['result'] != "PASS":
            print("  FAIL", c['check'], c['detail'])
    print("numbers: claim %.2f  W_raw 2^%.5f  W_stress 2^%.5f  S 2^%.4f  c_mul %.4f  c_sqr %.4f"
          % (claim, raw, st, log2(S), Lg['c_mul'], Lg['c_sqr']))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
