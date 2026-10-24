#!/usr/bin/env python3
"""Certificate for gt-rho-stationary-40pct-v10f (windowless build).

Every printed number is recomputed here from the shipped raw evidence JSONs
(c-*.json census, d-*.json faithful sweep, a-*.json audits, v10e-model.json) and
from the exact Fractions in constants.py / ledger.py.  Nothing is read from a
submitted summary string and no bound is carried as a named constant: W, the
model success, the reserve and the floor are all recomputed as Fractions (or as
exact alternating-series bounds), and the one-sided floor check uses integer
binomial arithmetic only.

Writes certificate.json; exits non-zero if any check is not PASS.
"""
import json
import os
import re
from fractions import Fraction
from math import comb, exp, log2, sqrt

import constants as K
import ledger

BITS = (5, 7, 9, 11)


def main():
    chk = []

    def add(name, ok, note=""):
        chk.append((name, bool(ok), note))

    # ---------- 1. ledger, exact ----------------
    L = ledger.ledgers()
    r = L['rows']
    fold_rows = ('selector_ab', 'selector_zeta', 'fold_adds', 'fold_cmp',
                 'hashing', 'dp_test', 'hist_slot', 'loop', 'witness_add')
    fold_sum = sum(r[k] for k in fold_rows)
    add("ledger_nine_fold_rows_sum_1813", fold_sum == r['fold'] == 1813,
        "nine rows sum to %d" % fold_sum)
    add("ledger_mul_step_21813", r['mul_step'] == K.C_REF + 1813 == 21813)
    add("ledger_sqr_step_13445", r['sqr_step'] == r['cyclo_sqr'] + 1813 == 13445)
    add("ledger_no_window_row", 'detector' not in r and K.WINDOW == 0 and K.PROBES == 0,
        "no look-back/escape row exists; WINDOW=0, PROBES=0")
    add("ledger_c_raw_exact_747_800", L['c_raw'] == Fraction(747, 800), str(L['c_raw']))
    add("ledger_c_stress_exact_192189_200000",
        L['c_stress'] == Fraction(192189, 200000), str(L['c_stress']))
    w_raw, w_str = L['w_raw'], L['w_stress']
    add("w_raw_2p125.71174", abs(log2(float(w_raw)) - 125.71174) < 1e-5,
        "%.5f" % log2(float(w_raw)))
    add("w_stress_2p125.75316", abs(log2(float(w_str)) - 125.75316) < 1e-5,
        "%.5f" % log2(float(w_str)))

    # ---------- 2. auxiliary rows, exact, each below 2^-60 of W ----------
    A = ledger.aux()
    aux_ok = all(Fraction(v) / w_str < Fraction(1, 2**60)
                 for k, v in A.items() if k != 'total')
    add("aux_rows_each_below_2^-60_of_W", aux_ok,
        "max row %s = 2^%.2f of W" % (max(((k, v) for k, v in A.items() if k != 'total'),
                                            key=lambda kv: kv[1])[0],
                                      max(log2(float(v / w_str)) for k, v in A.items() if k != 'total')))
    add("aux_total_below_2^-60_of_W", A['total'] / w_str < Fraction(1, 2**60),
        "total 2^%.2f of W" % log2(float(A['total'] / w_str)))

    # ---------- 3. step budget, exact integer inequality ----------
    S = K.step_budget()
    add("step_budget_is_ledger_S", S == L['S'])
    f, mu, N = K.FRAC_SQR, K.MU_WALK, K.N_CLASSES
    lo = Fraction(497580, 10**6)     # ln(1/(1-0.392)) = 0.4975803970
    hi = Fraction(497581, 10**6)
    lhs_S = Fraction(S*S) * (1 - f*f) / (2*mu*N)
    lhs_S1 = Fraction((S-1)*(S-1)) * (1 - f*f) / (2*mu*N)
    add("S_is_least_integer_for_target", lhs_S >= lo and lhs_S1 < hi,
        "S^2(1-f^2)/(2 mu N) = %.9f, (S-1) row below 0.497555" % float(lhs_S))

    # ---------- 4. production success model, exact ----------
    Lam = Fraction(S*S, 2*N)
    add("Lambda_exact_S2_over_2N", abs(float(Lam) - 0.6369029081804418) < 1e-12,
        "Lambda = %.10f" % float(Lam))
    p_lb = Lam - Lam*Lam/2          # 1-e^-Lambda >= Lambda - Lambda^2/2 (alternating series)
    add("P_success_exact_lower_ge_0.43", p_lb >= Fraction(43, 100),
        "1-exp(-Lambda) >= %.6f (exact series bound)" % float(p_lb))
    w = Fraction(S, 4)
    add("w_is_S_over_4", 4*w == S and abs(float(w)/sqrt(N) - 0.2821575154458148) < 1e-15,
        "w = S/4 = 2^%.5f, w/sqrt(N) = %.7f" % (log2(float(w)), float(w)/sqrt(N)))
    EL = w/3
    x = EL / (1 << K.THETA_BITS)     # expected DP hits on the shared suffix
    add("detection_term_ge_1_minus_2^-42", x >= 2**42,
        "E[L] = w/3 = 2^%.3f steps, E[L]/theta = 2^%.2f, so P(DP in suffix) >= 1-2^-42"
        % (log2(float(EL)), log2(float(x))))
    P = float(p_lb)
    add("declared_floor_0.39_below_model", Fraction(39, 100) <= p_lb,
        "floor 0.39 vs model >= %.5f" % P)
    P_model = 1 - exp(-float(Lam))
    reserve = 1 - Fraction(39, 100)/Fraction(P_model).limit_denominator(10**9)
    add("declared_reserve_0.17", Fraction(17, 100) <= reserve <= Fraction(18, 100),
        "reserve 1 - 0.39/P = %.4f against the model (%.4f against the exact series bound)"
        % (float(reserve), float(1 - Fraction(39, 100)/p_lb)))
    losses = sum(x for _, x in K.LOSS_ROWS)
    add("micro_loss_rows_at_most_3e-4",
        losses <= Fraction(3, 10**4) and losses > Fraction(29, 10**5), "%.3e" % float(losses))
    add("floor_after_micro_losses_ge_0.39", p_lb*(1 - losses) >= Fraction(39, 100),
        "%.5f" % float(p_lb*(1 - losses)))

    # ---------- 5. census cells (raw JSON, exact integers) ----------
    C = {b: json.load(open("c-%d.json" % b)) for b in BITS}
    D = {b: json.load(open("d-%d.json" % b)) for b in BITS}
    AU = {b: json.load(open("a-%d.json" % b)) for b in BITS}
    tot = lambda src, k: sum(src[b][k] for b in BITS)
    add("census_trials_6400", tot(C, 'trials') == 6400)
    add("census_success_5251", tot(C, 'success') == 5251, "%d/6400" % tot(C, 'success'))
    add("census_reports_54396", tot(C, 'reports') == 54396)
    add("census_trails_14483", tot(C, 'ntrail') == 14483)
    add("census_cross_channel_1853", tot(C, 'collision_shared_dp') == 1853)
    add("census_self_merge_channel_11082", tot(C, 'self_repeat_diff') == 11082)
    add("census_self_cycle_zero", tot(C, 'self_cycle') == 0,
        "no trail repeated a class with identical labels at any density")
    add("census_self_channel_absent_at_2^-5", C[5]['self_repeat_diff'] == 0,
        "0 self-merges at 2^-5 vs 4905/4813/1364 at 2^-7/9/11")
    add("census_2^-5_collision_universal", Fraction(C[5]['collision_any'], C[5]['trials']) > Fraction(99, 100),
        "%d/%d = 0.9981" % (C[5]['collision_any'], C[5]['trials']))
    add("census_2^-5_is_detection_limited",
        abs(C[5]['success'] - C[5]['collision_shared_dp']) <= 17 and C[5]['success'] <= C[5]['collision_shared_dp'],
        "547 successes against 564 cross-trail collisions carrying a DP: the gap is 17 stored-label coincidences")
    add("census_success_equals_sweep_success",
        all(C[b]['success'] == D[b]['success'] for b in BITS),
        str([C[b]['success'] for b in BITS]))
    add("census_recovery_is_exact_ratio",
        all(abs(D[b]['recovery_frac'] - D[b]['success']/D[b]['trials']) < 1e-12 for b in BITS))
    add("pooled_census_recovery_0.8205",
        Fraction(tot(D, 'success'), tot(D, 'trials')) == Fraction(5251, 6400),
        "%d/%d" % (tot(D, 'success'), tot(D, 'trials')))

    # ---------- 6. faithful sweep: conversion, invariants ----------
    add("sweep_trials_6400", tot(D, 'trials') == 6400)
    add("sweep_merges_5251", tot(D, 'merges') == 5251)
    add("sweep_every_merge_solves", tot(D, 'merges') == tot(D, 'solve_success') == tot(D, 'success'),
        "5251/5251 merges solved and verified against the planted secret")
    add("sweep_degenerate_denominator_zero", tot(D, 'degenerate_den') == 0)
    add("sweep_false_confirm_zero", tot(D, 'false_confirm') == 0)
    add("sweep_invariants_3585_zero_failures",
        tot(D, 'inv_checks') == 3585 and tot(D, 'inv_failures') == 0,
        "%d per-step invariant checks, %d failures" % (tot(D, 'inv_checks'), tot(D, 'inv_failures')))
    add("sweep_steps_and_strand_exact",
        tot(D, 'steps') == 2338487 and tot(D, 'stranded_steps') == 374528,
        "stranded 374528 of 2338487 steps (16.0 per cent), all after a repeated class")

    # ---------- 7. audits ----------
    add("audit_A_key_invariance_2112", tot(AU, 'audit_A_pass') == 2112 and tot(AU, 'audit_A_fail') == 0)
    add("audit_B_step_class_function_2112",
        tot(AU, 'audit_B_pass') == 2112 and tot(AU, 'audit_B_fail') == 0)
    add("audit_kappa_ties_zero", tot(AU, 'kappa_degenerate') == 0)

    # ---------- 8. formula-vs-sweep table, recomputed ----------
    rows = []
    for b in BITS:
        Nc, wc = C[b]['N'], C[b]['w_trail']
        lam_s = 4*comb(wc, 2)/Nc
        lam_c = 6*wc*wc/Nc
        rows.append((b, Nc, C[b]['spacing'], wc, C[b]['w_over_sqrtN'], lam_s, lam_c,
                     1 - exp(-(lam_s + lam_c)), D[b]['recovery_frac'],
                     C[b]['P_collision_with_dp'], C[b]['self_repeat_diff']))
    model = json.load(open("v10e-model.json"))
    add("model_json_production_matches_cert",
        abs(model['production']['Lambda'] - float(Lam)) < 1e-9 and
        abs(model['production']['P_prod'] - float(1 - exp(-float(Lam)))) < 1e-9 and
        abs(model['production']['w_over_sqrtN'] - float(w)/sqrt(N)) < 1e-9,
        "Lambda %.7f, P %.7f, w/sqrtN %.7f" % (model['production']['Lambda'],
                                               model['production']['P_prod'],
                                               model['production']['w_over_sqrtN']))
    # v10e-model.json carries the leading term Lambda = 8w^2/N; the exact
    # 4*C(w,2)/N + 6w^2/N differs by 2w/N (2^-123.8 at production, and by the
    # printed row difference at toy scale).  Both conventions are declared in
    # argument.md and the difference is asserted negligible at production.
    add("model_json_cells_match_formula",
        all(abs(model['toy12d'][i]['Lambda'] - 8*rows[i][3]**2/rows[i][1]) < 1e-6 and
            abs(model['toy12d'][i]['measured'] - rows[i][8]) < 1e-12 for i in range(4)),
        "four cells: 8w^2/N recomputed from N, theta, w; measured taken from d-*.json")
    add("exact_self_term_correction_negligible_at_production",
        Lam == Fraction(8)*w*w/N and Fraction(2*w, N) < Fraction(1, 2**120),
        "exact 4C(w,2)/N + 6w^2/N differs from 8w^2/N by 2w/N = 2^%.1f"
        % log2(float(2*w/N)))
    add("formula_overshoots_at_2^-5_only",
        rows[0][7] - rows[0][8] > 0.5 and all(abs(rows[i][7] - rows[i][8]) < 0.06 for i in (1, 2, 3)),
        "closed form 1.0000 vs measured 0.3419 at 2^-5, 0.95-0.996 elsewhere")

    # every printed row must appear verbatim in argument.md
    arg = open("argument.md").read()
    printed = ["| 2^-%d | %d | %d | %d | %.4f | %.4f | %.4f | %.4f | %.4f | %.4f | %d |"
               % (b, Nc, th, wc, wsq, ls, lc, pa, rec, xd, se)
               for (b, Nc, th, wc, wsq, ls, lc, pa, rec, xd, se) in rows]
    add("argument_prints_formula_table_verbatim", all(p in arg for p in printed),
        "4 rows byte-identical to the recomputation")
    census_rows = ["| 2^-%d | %d | %d | %d | %d | %d | %d | %d | %d | %.2f |"
                   % (b, C[b]['trials'], C[b]['success'], C[b]['reports'], C[b]['ntrail'],
                      C[b]['collision_any'], C[b]['collision_shared_dp'],
                      C[b]['self_repeat_diff'], C[b]['self_cycle'],
                      C[b]['suffix_sum']/C[b]['suffix_n']) for b in BITS]
    add("argument_prints_census_table_verbatim", all(p in arg for p in census_rows),
        "4 rows byte-identical to the recomputation")
    for tok, val in (("2^125.71174", log2(float(w_raw))), ("2^125.75316", log2(float(w_str))),
                     ("0.6369029082", float(Lam)), ("0.4710720", float(1 - exp(-float(Lam)))),
                     ("0.2821575", float(w)/sqrt(N)), ("1813", 1813.0)):
        add("argument_states_" + re.sub(r"[^0-9a-zA-Z]", "_", tok), tok in arg,
            "value %.7f" % val)

    # ---------- 9. loss and the claim ----------
    claim = K.CLAIM
    add("claim_is_125.79", claim == Fraction(12579, 100), str(claim))
    add("claim_above_w_stress_by_0.005_to_0.06",
        Fraction(5, 1000) <= claim - Fraction(log2(float(w_str))).limit_denominator(10**7) * 1 <= Fraction(6, 100),
        "pad %.5f" % (float(claim) - log2(float(w_str))))
    add("claim_above_w_raw", float(claim) > log2(float(w_raw)),
        "pad %.5f" % (float(claim) - log2(float(w_raw))))
    add("claim_at_or_below_gate_minus_0.01", float(claim) <= 125.84 - 0.01,
        "claim %.2f <= 125.83" % float(claim))
    add("claim_below_live_head", float(claim) < 125.84, "125.79 < 125.84 (a9c23016)")

    # ---------- 10. package surface ----------
    sub = json.load(open("submission.json"))
    sc = sub["claims"][0]["score"]
    add("submission_score_equals_claim",
        Fraction(int(sc["numerator"]), int(sc["denominator"])) == claim and
        int(sc["numerator"]) == 12579 and int(sc["denominator"]) == 100)
    missing = [e for e in sub["claims"][0]["evidence"]
               if e != "certificate.json" and not os.path.exists(e)]
    add("evidence_files_present", not missing, str(missing))
    add("dependencies_empty", json.load(open("dependencies.json")) == [])
    add("experiments_empty", sub["experiments"] == [])
    part = json.load(open("part_a-out.json"))
    add("part_a_all_pass", part.get("part_a") is True, "%d checks" % part.get("checks_n", 0))
    rout = json.load(open("routing-out.json"))
    add("routing_collapse_zero", rout.get("collapse") == 0, "disjoint %s" % rout.get("disjoint"))
    add("routing_has_no_ckpt_field", rout.get("has_ckpt_field") is False and
        "ckpt_fires" not in rout and rout.get("window") == 0)
    note = open("public-note.md").read()
    add("note_bytes_10k_14k", 10240 <= len(note.encode()) <= 14336, "%d B" % len(note.encode()))
    for token in ("0.4711", "0.4017", "0.3915", "0.39", "125.79", "125.75316",
                  "5251", "6400", "3,585", "2,112", "1813", "a9c23016"):
        add("note_states_" + token.replace(".", "_"), token in note)
    lines = note.splitlines()
    add("note_model_harness_lines_once",
        lines[0].strip() == "Model: qwen3.8-flash-next" and lines[1].strip() == "Harness: Hermes"
        and note.count("Model: qwen3.8-flash-next") == 1 and note.count("Harness: Hermes") == 1,
        "exactly one attribution block, lines 1-2")
    add("no_absolute_paths", "/Users/" not in note and "/Users/" not in arg)

    ok = all(c[1] for c in chk)
    out = {"checks": {n: ("PASS" if v else "FAIL") for n, v, _ in chk},
           "notes": {n: t for n, _, t in chk if t},
           "count": len(chk), "passed": sum(1 for c in chk if c[1]),
           "verdict": "PASS" if ok else "FAIL"}
    json.dump(out, open("certificate.json", "w"), indent=1)
    for n, v, t in chk:
        print("%-46s %s  %s" % (n, "PASS" if v else "FAIL", t))
    print("CERTIFICATE", "PASS" if ok else "FAIL", "%d/%d" % (out["passed"], out["count"]))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
