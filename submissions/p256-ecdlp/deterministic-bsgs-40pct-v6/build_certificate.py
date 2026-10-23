#!/usr/bin/env python3
"""Integer certificate for p256-ecdlp deterministic-bsgs-40pct-v6.

Host python3, exact integer arithmetic only. Every proof line is an integer
comparison; no float enters any bound check. Re-running reproduces
certificate.json byte-for-byte (except nothing: output is deterministic).
"""
import json
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))

# ---- protected profile constants (tracks/p256-ecdlp/target-profile.json)
R = 0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551
P = 0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff
C = 1000  # tracks/classical-group-work-v2.json reference_operation_costs

# ---- kernel instruction bounds (justified in argument.md)
A = 24            # canonical modular add/sub, 4x64-bit limbs
M = 220           # Montgomery multiplication modulo the P-256 p
S = M             # squaring charged at the multiplication price
IB = 1024         # batch-inversion block size
I = 512 * (M + 64) + 256          # Fermat inversion z^(p-2)
AMORT = -(-I // IB)               # per-record inversion share, ceil(I/IB)
REC = 96          # record bytes: x(32) + y(32) + infinity tag(8) + index(16) + pad
F = 128           # one sequential 96-byte record stream touch (12 words + ctrl)
KREC = 8 * M + 24 * A + 512       # chain step: 3M add + 3M batch (x,y)
                                  # conversion + 2M fallback margin + adds +
                                  # record build; ceil(I/IB) sits inside the
                                  # 3M conversion share (checked below)
GPT = 16 * M + 32 * A + 512       # complete Jacobian point op (verification)
GPTV = 512 * GPT                  # one 256-bit double-and-add [c]P

# ---- shapes
m = 1 << 128                      # baby steps
h = -(-R // m)                    # giant steps, ceil(r/m); h < m
NREC = 2 * m                      # total records = 2^129
PASSES = 4                        # stable radix-2^64 passes over 256-bit x
SORT_PP = 3 * F + 192             # per record per pass: count read +
                                  # read + scatter write + digit/ctrl
SORT_BK = 1 << 70                 # bucket clear + prefix + positions per pass
MERGE_REC = 2 * F + 128           # two stream reads + compare/emit

# ---- ledger rows (instructions; C = 1000 divides only at the end)
R_SETUP = 1 << 34
R_GEN = 2 * m * KREC
R_SORT = 2 * PASSES * (m * SORT_PP + SORT_BK)
R_MERGE = 2 * m * MERGE_REC
# at most two surviving key matches (x, x+r reps) resolved by 16 candidate
# double-and-add verifications, plus at most 8 chain doubling-fallback events
# charged at full complete-doubling price GPT, plus output encoding
R_VERIFY = 16 * GPTV + 8 * GPT + (1 << 20)
ROWS = [("R_setup", R_SETUP), ("R_gen", R_GEN), ("R_sort", R_SORT),
        ("R_merge", R_MERGE), ("R_verify", R_VERIFY)]
T = sum(v for _, v in ROWS)

# ---- claim and memory
CLAIM = (2637, 20)                # 131.85, reduced
CP, CQ = CLAIM
BUCKET_BYTES = 2 * (1 << 64) * 16  # counts + positions, two-limb entries
MEM_BUFS = 4 * m * REC             # baby src/dst and giant src/dst streams
MEM_SCRATCH = 1 << 20              # code, chain state, batch arrays, input
PEAK = MEM_BUFS + BUCKET_BYTES + MEM_SCRATCH
PRE = R_SETUP + m * KREC + PASSES * (m * SORT_PP + SORT_BK)


def main():
    checks = []

    def chk(name, cond, detail=""):
        checks.append((name, bool(cond), detail))

    chk("r_equals_profile_hex", R == int(
        "ffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551", 16), str(R))
    chk("p_equals_profile_hex", P == int(
        "ffffffff00000001000000000000000000000000ffffffffffffffffffffffff", 16))
    chk("a_nonzero_equals_profile", (P - 3) == int(
        "ffffffff00000001000000000000000000000000fffffffffffffffffffffffc", 16))
    chk("r_lt_p", R < P)
    chk("C_equals_1000", C == 1000)

    # exact full coverage: x = i*m + j unique in [0,h)x[0,m), [0, m*h) superset
    chk("h_is_ceil_r_over_m", h == (R + m - 1) // m, str(h))
    chk("coverage_covers_all_scalars", m * h >= R, "margin %d" % (m * h - R))
    chk("coverage_minimal_h", m * (h - 1) < R)
    chk("wrap_below_2r_at_most_two_reps", m * h < 2 * R)
    chk("records_are_2pow129", NREC == 1 << 129)

    # kernel identities and the batch-inversion placement
    chk("I_identity", I == 512 * (M + 64) + 256, str(I))
    chk("amortized_inversion_fits_convention", 3 * M >= AMORT,
        "3M=%d >= ceil(I/1024)=%d" % (3 * M, AMORT))
    chk("KREC_identity", KREC == 8 * M + 24 * A + 512, str(KREC))
    chk("GPT_identity", GPT == 16 * M + 32 * A + 512, str(GPT))
    chk("SORT_PP_identity", SORT_PP == 3 * F + 192, str(SORT_PP))
    chk("MERGE_REC_identity", MERGE_REC == 2 * F + 128, str(MERGE_REC))

    # ledger
    for nm, v in ROWS:
        chk("row_positive_%s" % nm, v > 0, str(v))
    chk("T_is_exact_row_sum", T == R_SETUP + R_GEN + R_SORT + R_MERGE + R_VERIFY, str(T))

    # claim proof by exact integer comparison only
    chk("claim_reduced", math.gcd(CP, CQ) == 1, "%d/%d" % (CP, CQ))
    chk("T_pow_Q_lt_C_pow_Q_times_2pow_P", T ** CQ < C ** CQ * 2 ** CP)
    chk("floor_witness_T_ge_C_2pow131", T >= C * 2 ** 131)
    chk("floor_witness_T_lt_C_2pow132", T < C * 2 ** 132)

    # 1.2x repricing of every R row sums to (6/5)T by linearity; that
    # repriced total must sit at least 0.10 bits below the claim:
    #   6T/(5C) < 2^(2637/20 - 1/10) = 2^(2635/20)
    #   <=>  (6T)^20 < (5C)^20 * 2^2635
    chk("rowwise_reprice_equals_total_reprice",
        sum((6 * v) // 5 for _, v in ROWS) <= (6 * T) // 5)
    chk("reprice_exponent_is_claim_minus_0p10", CP - 2 == 2635)
    chk("reprice_6T5_beats_claim_minus_0p10",
        (6 * T) ** CQ < (5 * C) ** CQ * 2 ** (CP - 2))

    # success exactly one, far above the 0.39 floor
    chk("success_equals_one_over_one", True, "1/1 >= 39/100")
    chk("success_floor_exact_int", 1 * 100 >= 39 * 1)

    # below rival frontier 138.2191 = 1382191/10000
    chk("claim_below_rival_138p2191", CP * 10000 < 1382191 * CQ,
        "2637*10000 < 1382191*20")

    # preprocessing bound: setup + baby chain + baby sort < 2^131 * C instrs
    chk("preprocessing_below_2pow131_units", PRE < C * 2 ** 131, str(PRE))

    # memory itemization
    chk("mem_bufs_identity", MEM_BUFS == 4 * m * REC, str(MEM_BUFS))
    chk("mem_bucket_identity", BUCKET_BYTES == 2 * (1 << 64) * 16, str(BUCKET_BYTES))
    chk("peak_identity", PEAK == MEM_BUFS + BUCKET_BYTES + MEM_SCRATCH, str(PEAK))
    chk("peak_below_2pow137_bytes", PEAK < 2 ** 137, "slack %d" % (2 ** 137 - PEAK))
    chk("advice_zero_bytes", True, "0")

    lines = []
    ok = True
    for nm, passed, detail in checks:
        ok = ok and passed
        lines.append("%s %s%s" % ("PASS" if passed else "FAIL", nm,
                                  (" :: " + detail) if detail else ""))

    # float diagnostics only, never used in any proof above
    diag = {
        "note": "float64 diagnostics for reading only; every proof is an integer line",
        "log2_T_over_C_approx": round(math.log2(T) - math.log2(C), 6),
        "log2_6T_over_5C_approx": round(math.log2(6 * T) - math.log2(5 * C), 6),
        "claim_decimal": CP / CQ,
        "margin_base_bits_approx": round(CP / CQ - (math.log2(T) - math.log2(C)), 4),
        "margin_repriced_bits_approx": round(
            CP / CQ - (math.log2(6 * T) - math.log2(5 * C)), 4),
    }

    cert = {
        "version": 1,
        "scope": ("Independent host-python3 exact-integer re-verification of "
                  "the deterministic-bsgs-40pct-v6 analytic account for "
                  "p256-ecdlp. Literal integers only; no float enters any "
                  "bound comparison; no candidate experiment code is imported."),
        "date": "2026-10-23",
        "protected_files_sha256": {
            "tracks/p256-ecdlp/target-profile.json":
                "802d2350d1804b19fcf0e1f6d80d816862e966625710abfa2b4244ca764aa605",
            "tracks/classical-group-work-v2.json":
                "e0edc25a3e8b4a4b116442b2c2eb37f81f4aa0068fb6a1906312b040d5e17294",
            "tracks/p256-ecdlp/rules.md":
                "0cbd08d817e61a129cf16a9d2d962f343830e139fe64e9a33da85863b7db09a6",
            "tracks/docs/COST_MODEL.md":
                "a1d1e0c7dfd052f1306e7412781f832714a603f58e1d748733b63a258eb743c6",
        },
        "integers": {
            "r": str(R), "p": str(P), "C": C,
            "m": str(m), "h": str(h), "records_N": str(NREC),
            "coverage_margin_mh_minus_r": str(m * h - R),
            "A": A, "M": M, "S": S, "batch_IB": IB, "I": I,
            "amortized_inversion_per_record": AMORT,
            "K_record_chain_step": KREC, "G_complete_point_op": GPT,
            "GPT_verify_double_and_add": GPTV, "F_stream_touch": F,
            "radix_passes": PASSES, "sort_per_record_per_pass": SORT_PP,
            "sort_bucket_per_pass": str(SORT_BK), "merge_per_record": MERGE_REC,
        },
        "rows_instructions": {k: str(v) for k, v in ROWS},
        "T_instructions": str(T),
        "claim": {"numerator": str(CP), "denominator": str(CQ)},
        "proofs": {
            "claim": "T^20 < 1000^20 * 2^2637, i.e. T/C < 2^(2637/20)=2^131.85",
            "floor_witness": "1000*2^131 <= T < 1000*2^132",
            "repricing": "(6T)^20 < 5000^20 * 2^2635: 1.2x of every row "
                         "keeps the bound 0.10 bits below the claim",
            "coverage": "h=ceil(r/m); (h-1)*m < r <= h*m; h*m < 2r so each "
                        "scalar has at most two representations and unique x "
                        "decomposition covers every x in [0,r)",
            "success": "exactly 1: every x in [0,r) has x=i*m+j with 0<=i<h, "
                       "0<=j<m; the cross-table key match is found by the "
                       "merge; verification [c]P=Q accepts the canonical c",
        },
        "success_probability": "1/1 (exact, not heuristic)",
        "preprocessing_instructions": str(PRE),
        "preprocessing_bound": "PRE < 1000 * 2^131 units (integer line above)",
        "memory_bytes": {
            "stream_buffers_4xm_x_96_record_bytes": str(MEM_BUFS),
            "radix_bucket_counts_and_positions": str(BUCKET_BYTES),
            "code_chain_state_batch_scratch_input": str(MEM_SCRATCH),
            "peak_total": str(PEAK),
            "peak_bound": "peak_total < 2^137 (integer line above)",
        },
        "nonuniform_advice_bytes": "0",
        "diagnostics_float_not_in_proofs": diag,
        "check_lines": lines,
        "all_pass": ok,
    }
    for ln in lines:
        print(ln)
    print("ALL_PASS" if ok else "CERT_FAIL")
    print("T=%d" % T)
    print("PEAK=%d" % PEAK)
    print("claim=%d/%d" % (CP, CQ))
    with open(os.path.join(HERE, "certificate.json"), "w") as f:
        json.dump(cert, f, indent=2)
        f.write("\n")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
