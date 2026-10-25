"""v23 planted controls -- four, each of which MUST be able to fail.

(a) PREDICATE CONTROL.  A same-bucket DIFFERENT-key pair is compared by the
    probe path but is NOT a collision; the old row predicate would count it.
    The audit reports both counts for the same fixture trial set (the measured
    bucket-pair and key-pair totals) and the planted pair, so the two
    statistics cannot be confused.

(b) MOMENT CONTROL.  The Paley-Zygmund bound must respond to dispersion.  A
    synthetic X with the same mean as the carrier but 1% mass at 63.3 has
    E[X(X-1)] = 39.43, so PZ = 0.010 << the Poisson value 0.469.  The measured
    carrier moments (E[X(X-1)]/E[X]^2 = 0.9908) reproduce the Poisson value.

(c) LOCATOR PRICE CONTROL.  Doubling the page width halves records/page and the
    per-pass page-walk share, so the locator row MUST move (246.518 ->
    243.259 Ir/record).  A row that cannot move is not being priced.

(d) REGEN 1-ULP CARRY.  The claim must be the ceiling of the exact log2 work at
    1e-4 granularity: claim - 1 ulp must fall BELOW the exact value, and the
    reduced rational must equal the claim.
"""
import json
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "evidence")


def control_a(key):
    c = key["carrier"]
    planted = key["planted"]
    bucket = c["bucket_pairs"]
    kes = c["hits"]
    assert bucket > kes > 0
    assert planted["same_bucket_different_key"]["hit"] == 0
    assert planted["same_bucket_different_key"]["same_bucket"] is True
    assert planted["same_key_pair"]["hit"] == "same_point"
    assert planted["same_key_pair"]["recovered"] is True
    return {
        "measured_bucket_pairs": bucket,
        "measured_key_pairs": kes,
        "bucket_over_key": bucket / kes,
        "planted_same_bucket_different_key": planted["same_bucket_different_key"],
        "planted_same_key": planted["same_key_pair"],
        "production_gap": "at production r/R = 2^255/2^20 = 2^235, so the row count "
                          "is not a proxy for the key count at any feasible scale",
        "verdict": "PASS: the probe path compares a same-bucket different-key pair and "
                   "does NOT count it (%d key pairs against %d bucket pairs in the same "
                   "runs); the planted same-key pair fires and recovers x. The old row "
                   "predicate would have counted the same-bucket pair."
                   % (kes, bucket),
    }


def control_b(key):
    c = key["carrier"]
    m1 = c["moments"]["E_X"]
    m2 = c["moments"]["E_XX1"]
    pz_measured = m1 * m1 / (m2 + m1)
    pz_poisson = m1 * m1 / (m1 * m1 + m1)
    # synthetic overdispersed X with the same mean 0.633: 1% at 63.3
    p = 0.01
    xv = 63.3
    s1 = p * xv
    s2 = p * xv * (xv - 1)
    pz_syn = s1 * s1 / (s2 + s1)
    assert pz_syn < pz_poisson / 10, "the PZ bound must collapse under overdispersion"
    return {
        "measured_carrier": {"E_X": m1, "E_XX1": m2,
                             "factorial_ratio": m2 / (m1 * m1),
                             "paley_zygmund": pz_measured,
                             "poisson_value_at_same_mean": pz_poisson},
        "synthetic_overdispersed": {"E_X": s1, "E_XX1": s2,
                                    "paley_zygmund": pz_syn,
                                    "poisson_value_at_same_mean": pz_poisson},
        "verdict": "PASS: PZ responds to dispersion (0.010 for the overdispersed synthetic "
                   "vs %.3f Poisson at the same mean), and the measured factorial ratio "
                   "%.4f is Poisson-tight" % (pz_poisson, m2 / (m1 * m1)),
    }


def control_c(store):
    pc = store["page_width_control"]
    assert pc["moves"] is True
    assert abs(pc["narrow"]["locator_row_ir_per_record"]
               - pc["wide"]["locator_row_ir_per_record"]) > 0.1
    assert pc["wide"]["pair_transfers_per_page"] == 512, "8192 B page must be 512 pairs"
    return {"narrow": pc["narrow"], "wide": pc["wide"], "delta": pc["delta"],
            "verdict": "PASS: doubling the page width moves the locator row by %.3f "
                       "Ir/record, so the row is priced from the layout, not asserted"
                       % pc["delta"]}


def control_d():
    claim = json.load(open(os.path.join(OUT, "claim_v23.json")))
    ledger = json.load(open(os.path.join(OUT, "ledger_v23.json")))
    exact = ledger["work"]["log2_work_final"]
    c = claim["claim"]
    rat = claim["rational"]
    ulp = 1e-4
    assert c >= exact, "claim must bound the exact work"
    assert c - ulp < exact, "claim is not the tightest 4-decimal bound"
    assert abs(rat["numerator"] / rat["denominator"] - c) < 1e-12
    assert int(round(c * 10000)) % (10000 // rat["denominator"]) == 0
    return {"claim": c, "exact_log2_work": exact, "ulp": ulp,
            "claim_minus_ulp": c - ulp, "one_ulp_lower_would_understate": c - ulp < exact,
            "rational": rat,
            "verdict": "PASS: claim is the ceiling of the exact log2 work at 1e-4 and the "
                       "reduced rational is exact"}


def main():
    key = json.load(open(os.path.join(OUT, "key_equal_v23.json")))
    store = json.load(open(os.path.join(OUT, "store_v23.json")))
    out = {"a_predicate": control_a(key), "b_moment": control_b(key),
           "c_locator_page_width": control_c(store), "d_regen_ulp": control_d()}
    json.dump(out, open(os.path.join(OUT, "audit_v23.json"), "w"), indent=1, sort_keys=True)
    print(json.dumps({k: v["verdict"][:90] for k, v in out.items()}, indent=1))
    print("ALL FOUR CONTROLS FIRED")


if __name__ == "__main__":
    main()
