"""Recompute the published STNFS cost for BLS12-381 and the submitted ledger.

This is not a run of the attack. It takes the published simulation parameters
of Guillevic-Singh (Mathematical Cryptology 1(1), 2021, Table 10, curve
BLS-381), recomputes the relation-collection and block-Wiedemann costs of the
same model, and applies the submission's conservative instruction-equivalent
conversion. It is evidence for the arithmetic and for the published-input
reproduction only; the model assumptions are stated in proof.md.
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
    overhead = float(cfg["overhead_factor"])
    C = float(cfg["reference_price"])
    published_total = float(cfg["published_total_log2"])
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

    # Conservative conversion: model operations are calibrated cycles; at most
    # 8 scalar instructions retire per cycle, and a factor 16 covers the worst
    # measured sieving constant (27 cycles/tuple = 7.08 cycles per model
    # operation) plus a 1.25x allowance for descent, ECM, retries and
    # verification. Then divide by the fixed price C.
    instr_log2 = total_log2 + math.log2(ipc * overhead)
    score = instr_log2 - math.log2(C)
    # The field-only core dominates; report the full charged total as the
    # per-instance preprocessing upper bound (no reuse credit is taken).
    preproc = score

    # Reportable peak storage: the filtered matrix and, if held, the relation
    # file. 40 bytes = 8-byte index + four 64-bit value words per nonzero.
    rows = 2.0 ** log2_fb / flt
    matrix_bytes = rows * wt * 40.0
    relation_bytes = (2.0 ** log2_rel) * entries_per_relation * 8.0
    memory_bytes = max(matrix_bytes, relation_bytes) + 2.0 ** 20

    checks = {
        "space_matches_published_2_124_08": abs(space_log2 - 124.08) < 0.05,
        "total_within_0_25_bits_of_published": abs(total_log2 - published_total) < 0.25,
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
        "instruction_budget_log2": instr_log2,
        "score_log2": score,
        "preprocessing_log2": preproc,
        "memory_log2_bytes": math.log2(memory_bytes),
        "checks": checks,
        "scope": "arithmetic reproduction of a published cost model; not a production attack run",
    }))


if __name__ == "__main__":
    main()
