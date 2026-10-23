"""Recompute the published STNFS cost for BLS12-381 and the hard-capped ledger.

This is not a run of the attack. It takes the published simulation parameters of
Guillevic-Singh (Mathematical Cryptology 1(1), 2021, Table 10, curve BLS-381),
recomputes the relation-collection and block-Wiedemann costs of that model, then
applies this submission's instruction conversion (16 cycle-equivalents per model
operation at 8 scalar instructions per cycle), its charged retry factor and its
Markov truncation factor, and finally the fixed reference price and the explicit
verification charge. The result is the submitted hard work cap.
"""
import json
import math
import sys


def main():
    cfg = json.load(sys.stdin)
    A = float(cfg["sieve_bound"])
    eta = int(cfg["eta"])
    log2_B = float(cfg["log2_smoothness_bound"])
    log2_fb = float(cfg["log2_factor_base"])
    log2_rel = float(cfg["log2_relations"])
    wt = int(cfg["matrix_weight"])
    flt = float(cfg["filtering_factor"])
    words = int(cfg["la_words"])
    ipc = float(cfg["ipc_bound"])
    cycles_per_model_op = float(cfg["cycles_per_model_op"])
    retry = float(cfg["retry_factor"])
    cap_factor = float(cfg["cap_factor"])
    C = float(cfg["reference_price"])
    published_total = float(cfg["published_total_log2"])
    verification_log2 = float(cfg["verification_reference_ops_log2"])
    claimed_score = float(cfg["claimed_score_log2"])
    entries_per_relation = int(cfg["relation_entries_bound"])

    # Relation-collection model (GS21 Sec. 6.1): volume/2 * log(log(B)).
    space_log2 = 2 * eta * math.log2(2 * A + 1) - 1.0
    lnB = log2_B * math.log(2.0)
    sieve_log2 = space_log2 + math.log2(math.log(lnB))

    # Linear algebra (GS21 Eq. 6.5): words * wt * (#B/flt)^2.
    la_log2 = math.log2(words * wt) + 2.0 * (log2_fb - math.log2(flt))

    total = 2.0 ** sieve_log2 + 2.0 ** la_log2
    total_log2 = math.log2(total)

    # Instruction conversion: cycles_per_model_op must dominate the worst
    # measured model/hardware ratio published for this model family (BD19
    # Table 6 reports log2(c_sieve) up to 3.5), and ipc must dominate the widest
    # scalar retirement width of a current AArch64 core (8).
    instr_per_op = cycles_per_model_op * ipc
    single_log2 = total_log2 + math.log2(instr_per_op)
    expected_log2 = single_log2 + math.log2(retry)      # E[T], Markov input only
    cap_log2 = expected_log2 + math.log2(cap_factor)   # enforced hard cap
    score = math.log2(2.0 ** (cap_log2 - math.log2(C)) + 2.0 ** verification_log2)
    # The field-only core dominates; the full charged total is reported as the
    # per-instance preprocessing bound (no reuse credit is taken).
    preproc = score

    # Peak storage: the filtered matrix and, if held, the relation file.
    # 40 bytes = 8-byte index + four 64-bit value words per nonzero.
    rows = 2.0 ** log2_fb / flt
    matrix_bytes = rows * wt * 40.0
    relation_bytes = (2.0 ** log2_rel) * entries_per_relation * 8.0
    memory_bytes = max(matrix_bytes, relation_bytes) + 2.0 ** 20

    checks = {
        "space_matches_published_2_124_08": abs(space_log2 - 124.08) < 0.05,
        "total_within_0_25_bits_of_published": abs(total_log2 - published_total) < 0.25,
        "calibration_dominates_worst_measured_ratio": cycles_per_model_op >= 2.0 ** 3.5,
        "ipc_dominates_widest_aarch64_retirement": ipc >= 8.0,
        "cap_is_two_times_expected": abs(cap_log2 - expected_log2 - 1.0) < 1e-9,
        "score_matches_claim_within_1e_6": abs(score - claimed_score) < 1e-6,
    }
    if not all(checks.values()):
        raise SystemExit("ledger check failed: " + json.dumps(checks))

    print(json.dumps({
        "curve": "BLS12-381, Fp12 STNFS, GS21 Table 10 BLS-381 parameters",
        "sieve_bound_A": A,
        "log2_smoothness_bound_B": log2_B,
        "log2_sieving_space_recomputed": space_log2,
        "log2_sieve_cost": sieve_log2,
        "log2_linear_algebra_cost": la_log2,
        "log2_total_model_cost": total_log2,
        "cycles_per_model_op": cycles_per_model_op,
        "instructions_per_cycle": ipc,
        "instructions_per_model_op": instr_per_op,
        "log2_single_attempt_instructions": single_log2,
        "log2_expected_total_instructions_E_T": expected_log2,
        "log2_hard_cap_instructions_R_cap": cap_log2,
        "retry_factor": retry,
        "cap_factor": cap_factor,
        "verification_log2_reference_ops": verification_log2,
        "score_log2": score,
        "preprocessing_log2": preproc,
        "memory_log2_bytes": math.log2(memory_bytes),
        "checks": checks,
        "scope": "arithmetic reproduction of a published cost model with this "
                 "submission's hard-capped conversion; not a production attack run",
    }, indent=1))


if __name__ == "__main__":
    main()
