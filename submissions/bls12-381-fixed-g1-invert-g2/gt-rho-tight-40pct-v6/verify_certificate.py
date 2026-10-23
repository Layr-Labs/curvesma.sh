#!/usr/bin/env python3
"""Re-derive every number in the claim and check every recorded run.

Checks, all recomputed here rather than quoted:
  * the itemisation sums (fold = sum of its ten rows, step totals)
  * the step budget S from the single-source constants, as an exact integer
  * W_raw and W_stress as exact Fractions, and the claim-window rules, including
    W_raw + 2e-5 <= claim < live gate as an exact rational inequality
  * the corrected-recovery row: the measured pooled recovery as an exact Fraction must
    be at or above the 0.39 requirement and inside three binomial standard errors of
    the model, in both the f = 3/8 walk and the f = 0 control
  * the equivariance row: the exponent-domain probe must agree on every sample under the
    v6 rule, and the v5 rule must demonstrably fail it on the same samples (the defect
    is reproduced, not asserted)
  * bit-disjointness of the routing fields, the recorded collapse test, the matched
    density row (realised steps and reports against the budget)
  * every recorded check from part_a / sim / toy_core / routing
  * single-constant discipline: no module restates a walk constant or a price
"""
import json, os, sys
from fractions import Fraction
from math import log, log2
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
             'dp_test', 'detector', 'hist_slot', 'loop', 'witness_add')
    chk("itemisation_fold_sum", r['fold'] == sum(r[k] for k in parts),
        "fold = %d from %d rows" % (r['fold'], len(parts)))
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

    # ---- claim window (exact rationals; the gate is the live 125.90 head)
    claim = K.CLAIM
    gate = Fraction(12590, 100)
    w_raw = Lg['w_raw']; w_st = Lg['w_stress']
    raw = log2(float(w_raw)); st = log2(float(w_st))
    chk("claim_strictly_below_gate", claim < gate,
        "claim %.2f < live gate 125.90" % float(claim))
    raw_f = Fraction.from_float(raw).limit_denominator(10**9)
    st_f = Fraction.from_float(st).limit_denominator(10**9)
    chk("claim_above_raw_plus_2e_5", claim - raw_f >= Fraction(2, 100000),
        "raw pad %.6f" % (float(claim) - raw))
    chk("claim_above_stress_plus_0.005", claim - st_f >= Fraction(5, 1000),
        "stress pad %.6f" % (float(claim) - st))
    chk("claim_at_most_stress_plus_0.06", claim - st_f <= Fraction(6, 100),
        "stress pad %.6f" % (float(claim) - st))

    # ---- success floor
    success = 1 - pow(2.718281828459045,
                      -S*S*(1-float(K.FRAC_SQR)**2)/(2*float(K.MU_WALK)*K.N_CLASSES))
    chk("budgeted_success_at_mu", abs(success - float(K.TARGET)) < 1e-4,
        "model %.4f" % success)
    chk("success_floor_above_requirement", K.TARGET - K.LOSSES > K.SUCCESS,
        "%.4f - %.4f = %.4f > %.2f" % (K.TARGET, K.LOSSES, float(K.TARGET-K.LOSSES),
                                       float(K.SUCCESS)))

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
        n_ctrl = sum(x['groups'] for x in sm['control'])
        chk("sim_group_count", n_groups >= 4000 and n_ctrl >= 4000,
            "%d mixed groups, %d control groups" % (n_groups, n_ctrl))

        # corrected-recovery row, measured pooled recovery as an exact Fraction
        rec = Fraction(int(round(sm['pooled_p']*10**6)), 10**6)
        chk("recovery_row_above_requirement", rec >= K.SUCCESS,
            "measured %s = %.6f >= 0.39" % (rec, sm['pooled_p']))
        chk("recovery_row_within_3se_of_model",
            abs(sm['pooled_p'] - sm['model']) <= 3*sm['pooled_se'],
            "measured %.4f model %.4f se %.4f" % (sm['pooled_p'], sm['model'],
                                                  sm['pooled_se']))
        ctl = Fraction(int(round(sm['control_p']*10**6)), 10**6)
        chk("recovery_control_above_requirement", ctl >= K.SUCCESS,
            "f = 0 control %.4f >= 0.39" % sm['control_p'])
        chk("recovery_control_within_3se_of_model",
            abs(sm['control_p'] - sm['model_control']) <= 3*sm['control_se'],
            "measured %.4f model %.4f se %.4f" % (sm['control_p'], sm['model_control'],
                                                  sm['control_se']))

        # the v5 rule, run as a control on the same fields, tables and seeds
        chk("legacy_rule_reproduces_the_defect",
            sm['legacy_p'] < 0.5*sm['pooled_p']
            and sm['legacy_p'] < sm['model_legacy'] - 3*sm['legacy_se'],
            "v5 rule %.4f vs v6 %.4f, model %.4f" % (sm['legacy_p'], sm['pooled_p'],
                                                     sm['model_legacy']))

        # equivariance row
        eq = sm['equivariance']
        chk("equivariance_v6_all_samples",
            eq['v6_agree'] == eq['samples'] and eq['samples'] >= 1000,
            "v6 %d/%d agree" % (eq['v6_agree'], eq['samples']))
        chk("equivariance_v5_rule_fails", eq['v5_agree'] < 0.6*eq['samples'],
            "v5-style fixed index %d/%d" % (eq['v5_agree'], eq['samples']))

        # matched-density disclosure row
        chk("matched_density_steps", 0.95 <= sm['steps_over_S'] <= 1.05,
            "realised steps / S = %.4f" % sm['steps_over_S'])
        chk("matched_density_reports",
            abs(sm['reports_per_group'] - sm['reports_expected']) <= 0.05*sm['reports_expected'],
            "reports/group %.2f vs expected %.2f" % (sm['reports_per_group'],
                                                     sm['reports_expected']))
        chk("sim_uses_production_constants",
            abs(sm['f'] - float(K.FRAC_SQR)) < 1e-12 and abs(sm['mu'] - float(K.MU_WALK)) < 1e-12)
        chk("sim_squaring_fraction_is_f", abs(sm['sqr_fraction'] - float(K.FRAC_SQR)) <= 0.002,
            "measured %.4f" % sm['sqr_fraction'])
        chk("sim_control_squaring_fraction_zero", sm['control_sqr_fraction'] == 0.0,
            "measured %.4f" % sm['control_sqr_fraction'])
        chk("sim_injectivity_penalty_is_conservative",
            abs(sm['injectivity_ratio'] - sm['injectivity_model']) <= 0.005*sm['injectivity_model'],
            "measured %.4f charged %.4f (difference %.4f bits)"
            % (sm['injectivity_ratio'], sm['injectivity_model'],
               abs(sm['injectivity_ratio'] - sm['injectivity_model'])/sm['injectivity_model']/0.693))

    ty = load("toy-out.json")
    chk("toy_present", ty is not None)
    if ty:
        chk("toy_end_to_end_verified", ty['all_verified'],
            "%d/%d recovered and verified" % (ty['verified'], ty['recovered']))
        chk("toy_secret_match", ty['secret_match'] == ty['recovered'] > 0)
        chk("toy_cap_refusal", ty['cap_probe_refused'])
        chk("toy_no_alignment_failure", ty['aligned_fail'] == 0)
        chk("toy_production_shape", ty['N'] >= 500000 and ty['p'] % 12 == 7)
        cr = ty['detector_cycle_test']
        chk("toy_cycle_regression_448", cr['fired'] == 448 and cr['all_phases_fire'],
            "%d/%d (lengths 2..8 x 64 phases)" % (cr['fired'], 7*64))
        ca = ty['detector_cycle_adversarial']
        chk("toy_cycle_adversarial_all_fire", ca['all_phases_fire'],
            str(ca['phases_by_length']))
        chk("toy_equivariance_all_twelve",
            ty['equivariance_probe']['all_twelve_agree'] == ty['equivariance_probe']['samples']
            and ty['equivariance_probe']['samples'] > 0,
            "%d/%d samples" % (ty['equivariance_probe']['all_twelve_agree'],
                               ty['equivariance_probe']['samples']))

    # ---- single-constant discipline: numeric literals only, so prose and comments are
    # ignored; every price, the target, the mixture and the claim come from constants.py.
    import tokenize, io
    forbidden_numbers = {20000, 383, Fraction(392, 1000), Fraction(3, 8),
                         Fraction(12587, 100), 0.392, 0.375, 125.87}
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
    raw = log2(float(w_raw)); st = log2(float(w_st))
    cert = {"package": "gt-rho-tight-40pct-v6",
            "certificate": "PASS" if ok else "FAIL",
            "checks": CERT,
            "numbers": {
                "claim": str(claim), "claim_float": float(claim),
                "W_raw_log2": raw, "W_stress_log2": st,
                "S": S, "S_log2": log2(S),
                "c_mul": float(Lg['c_mul']), "c_sqr": float(Lg['c_sqr']),
                "c_raw": float(Lg['c_raw']),
                "fold_instructions": r['fold'], "cyclo_sqr_instructions": r['cyclo_sqr'],
                "raw_pad": float(claim) - raw, "stress_pad": float(claim) - st,
                "budgeted_success": success,
                "success_floor": float(K.TARGET - K.LOSSES),
                "measured_recovery": sm['pooled_p'] if sm else None,
                "measured_recovery_se": sm['pooled_se'] if sm else None,
                "measured_recovery_control": sm['control_p'] if sm else None,
                "measured_recovery_model": sm['model'] if sm else None,
                "measured_recovery_legacy_rule": sm['legacy_p'] if sm else None,
                "equivariance_v6": sm['equivariance']['v6_agree'] if sm else None,
                "equivariance_samples": sm['equivariance']['samples'] if sm else None,
                "equivariance_v5_rule": sm['equivariance']['v5_agree'] if sm else None,
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
          % (float(claim), raw, st, log2(S), Lg['c_mul'], Lg['c_sqr']))
    if sm:
        print("recovery: v6 %.4f +- %.4f (model %.4f)  control %.4f  v5 rule %.4f"
              % (sm['pooled_p'], sm['pooled_se'], sm['model'], sm['control_p'], sm['legacy_p']))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
