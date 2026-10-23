"""Recompute the STNFS ledger for fixed-G1 pairing inversion into G2.

This is not a run of the attack. It takes the published simulation parameters of
Guillevic-Singh (Mathematical Cryptology 1(1), 2021, Table 10, curve BLS-381),
recomputes the relation-collection and block-Wiedemann costs, derives the
individual-logarithm (descent) cost from p^12 with the published L[1/3,1.44]
constant, adds the ECM cofactorization, applies this submission's instruction
conversion (1024 scalar instructions per GS21 model operation), its capped
per-attempt budget, its charged retry factor and its Markov truncation factor,
and finally the organizer's fixed reference price and the explicit verification
charge. The result is the submitted hard work cap.

Every number that appears in proof.md or submission.json is recomputed here from
the inputs below; nothing is hard-coded except the published constants that the
inputs name.
"""
import json
import math
import sys

P_HEX = ("0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf67"
         "30d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab")


def recompute(cfg):
    # --- published GS21 Table 10 parameters -------------------------------
    A = float(cfg["sieve_bound"])                    # special-q bound
    eta = int(cfg["eta"])                            # exponent of the volume
    log2_B = float(cfg["log2_smoothness_bound"])     # B = 2^65.316
    log2_fb = float(cfg["log2_factor_base"])         # #FB = 2^60.848
    wt = int(cfg["matrix_weight"])
    flt = float(cfg["filtering_factor"])
    words = int(cfg["la_words"])

    # --- this submission's conversions (each justified in proof.md) -------
    instr_per_op = float(cfg["instructions_per_model_op"])   # 1024
    success_floor = float(cfg["success_floor"])              # organizer floor 0.39
    retry = 1.0 / success_floor                              # E[#attempts] <= 1/p
    cap_factor = float(cfg["cap_factor"])                    # 2
    C = float(cfg["reference_price"])                        # 20000
    descent_coeff = float(cfg["descent_coefficient"])        # 1.44
    ecm_log2 = float(cfg["cofactorization_log2"])            # 87.6
    verification_log2 = float(cfg["verification_reference_ops_log2"])  # 14
    claimed_score = float(cfg["claimed_score_log2"])
    bytes_per_relation = int(cfg["relation_bytes_per_entry"])
    log2_relations = float(cfg["log2_relations"])
    merge_passes = int(cfg["merge_passes"])
    compare_word_ops_log2 = float(cfg["compare_word_ops_log2"])
    published_total = float(cfg["published_total_log2"])     # 126
    p_int = int(cfg["field_prime_hex"], 16)

    # --- relation collection (GS21 Sec. 6.1) ------------------------------
    space_log2 = 2 * eta * math.log2(2 * A + 1) - 1.0
    lnB = log2_B * math.log(2.0)
    sieve_log2 = space_log2 + math.log2(math.log(lnB))

    # --- block-Wiedemann (GS21 Eq. 6.5) -----------------------------------
    la_log2 = math.log2(words * wt) + 2.0 * (log2_fb - math.log2(flt))
    core = 2.0 ** sieve_log2 + 2.0 ** la_log2
    core_log2 = math.log2(core)

    # --- descent from p^12, L[1/3, c] ------------------------------------
    ln_N = 12.0 * math.log(p_int)
    theta = ln_N ** (1.0 / 3.0) * math.log(ln_N) ** (2.0 / 3.0)
    descent_log2 = descent_coeff * theta / math.log(2.0)
    descent = 2.0 ** descent_log2
    ecm = 2.0 ** ecm_log2

    total = core + descent + ecm
    total_log2 = math.log2(total)

    # --- cap ledger -------------------------------------------------------
    single_log2 = total_log2 + math.log2(instr_per_op)      # B, hard cap
    expected_log2 = single_log2 + math.log2(retry)          # E[T]
    cap_log2 = expected_log2 + math.log2(cap_factor)        # R_cap = 4*B
    score = math.log2(2.0 ** (cap_log2 - math.log2(C)) + 2.0 ** verification_log2)
    score_rounded = math.ceil(score * 1e6) / 1e6

    # --- storage schedule (Sec. 7) ----------------------------------------
    relation_bytes_log2 = log2_relations + math.log2(bytes_per_relation)
    matrix_bytes_log2 = float(cfg["matrix_bytes_log2"])
    sort_ops_log2 = math.log2(merge_passes) + compare_word_ops_log2

    # --- first-principles sieving instruction count (Sec. 4.2) -----------
    # A sieved pair receives one update per (prime, root) hit.  Summed over the
    # factor base up to B the expected hit count per pair is 2*ln ln B, and a
    # hit is a load/add/store on an 8-bit log array, i.e. 5 scalar
    # instructions; an untouched pair still costs one initialisation byte.
    hits_per_pair = 2.0 * math.log(lnB)
    instr_per_pair_first_principles = 8.0 + 5.0 * hits_per_pair
    charged_instr_per_pair = 2.0 ** (sieve_log2 - space_log2) * instr_per_op
    per_pair_margin = charged_instr_per_pair / instr_per_pair_first_principles

    checks = {}
    checks["space_matches_published_2_124_08"] = abs(space_log2 - 124.0774) <= 0.05
    checks["core_matches_published_2_126"] = abs(core_log2 - published_total) <= 0.25
    checks["descent_base_is_2_122_6"] = (abs(descent_log2 - 122.6418) <= 0.01
                                        and descent < core and descent > 2.0 ** 120)
    checks["descent_is_subdominant_but_charged"] = (
        descent / core < 0.5 and descent / core > 2.0 ** -8)
    checks["conversion_dominates_published_calibration"] = (
        instr_per_op >= 8.0 * (2.0 ** 3.5))
    checks["conversion_is_1024"] = instr_per_op == 1024.0
    checks["retry_is_inverse_success_floor"] = (
        abs(retry - 1.0 / success_floor) < 1e-12)
    checks["success_floor_is_organizer_floor"] = abs(success_floor - 0.39) < 1e-12
    checks["cap_is_two_expected"] = (
        abs(cap_log2 - (expected_log2 + 1.0)) < 1e-12)
    checks["expected_is_inverse_floor_times_B"] = (
        abs(expected_log2 - (single_log2 + math.log2(1.0 / 0.39))) < 1e-12)
    checks["charged_per_pair_dominates_first_principles"] = per_pair_margin >= 50.0
    checks["score_rounded_up_is_claimed"] = score_rounded == claimed_score
    checks["claimed_score_matches_rational"] = (
        abs(claimed_score - 124.339147) < 1e-12)
    checks["relation_file_is_2_72_39"] = abs(relation_bytes_log2 - 72.389800) <= 1e-6
    checks["sort_is_negligible"] = sort_ops_log2 < total_log2 - 40.0
    checks["memory_peak_is_relation_file"] = relation_bytes_log2 > matrix_bytes_log2
    checks["all_terms_positive_and_finite"] = all(
        x > 0.0 and math.isfinite(x)
        for x in (space_log2, sieve_log2, la_log2, descent_log2, total_log2,
                  single_log2, expected_log2, cap_log2, score))
    # this revision charges 8x the previous one and is deliberately weaker
    checks["score_is_weaker_than_the_128_instruction_revision"] = score > 120.98
    return {
        "curve": "BLS12-381, Fp12 STNFS, GS21 Table 10 BLS-381 parameters",
        "log2_field_prime": math.log2(p_int),
        "log2_p_pow_12": ln_N / math.log(2.0),
        "theta_lnN_third_lnlnN_two_thirds": theta,
        "sieve_bound_A": A,
        "log2_smoothness_bound_B": log2_B,
        "log2_sieving_space_recomputed": space_log2,
        "log2_sieve_cost": sieve_log2,
        "log2_linear_algebra_cost": la_log2,
        "log2_core_model_cost": core_log2,
        "descent_coefficient": descent_coeff,
        "log2_descent_cost": descent_log2,
        "log2_cofactorization_cost": ecm_log2,
        "log2_total_model_cost": total_log2,
        "instructions_per_model_op": instr_per_op,
        "log2_single_attempt_instructions_B": single_log2,
        "log2_expected_total_instructions_E_T": expected_log2,
        "log2_hard_cap_instructions_R_cap": cap_log2,
        "success_floor": success_floor,
        "retry_factor": retry,
        "hits_per_sieved_pair": hits_per_pair,
        "first_principles_instr_per_pair": instr_per_pair_first_principles,
        "charged_instr_per_pair": charged_instr_per_pair,
        "per_pair_margin": per_pair_margin,
        "reference_price_C": C,
        "verification_log2_reference_ops_H": verification_log2,
        "score_log2": score,
        "score_rounded_up_6dp": score_rounded,
        "preprocessing_log2": score,
        "log2_relation_file_bytes": relation_bytes_log2,
        "log2_matrix_bytes": matrix_bytes_log2,
        "log2_sort_word_ops": sort_ops_log2,
        "checks": checks,
        "scope": ("arithmetic reproduction of a published cost model under this "
                  "submission's hard-capped conversion; not a production attack run"),
    }


def main():
    cfg = json.load(sys.stdin)
    out = recompute(cfg)
    failed = [k for k, v in out["checks"].items() if not v]
    if failed:
        raise SystemExit("check(s) failed: " + json.dumps(failed))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
