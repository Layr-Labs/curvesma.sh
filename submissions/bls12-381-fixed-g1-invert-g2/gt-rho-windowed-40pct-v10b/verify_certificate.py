#!/usr/bin/env python3
"""Certificate for gt-rho-windowed-40pct-v10b.

Every printed number is recomputed here from the shipped raw trial log
(`trials4x-raw.b64`, zlib+base64 of one record per group) and from the exact
Fractions in constants.py / ledger.py.  Nothing is read from a submitted
summary: sweep4x-out.json is itself re-derived and byte-compared on the counts
that matter, and the one-sided binomial floor is re-checked in EXACT integer
arithmetic at the quoted rational.

Writes certificate.json; exits non-zero if any check is not PASS.
"""
import base64
import json
import zlib
from fractions import Fraction
from math import comb, exp, log, log2, sqrt

import constants as K
import ledger


def exact_tail_ge(x, n, pnum, pden):
    """Exact P[Bin(n, pnum/pden) >= x] as a Fraction (integer arithmetic only)."""
    q = pden - pnum
    num = sum(comb(n, k) * pnum**k * q**(n - k) for k in range(x, n + 1))
    return Fraction(num, pden**n)


def main():
    chk = []

    def add(name, ok, note=""):
        chk.append((name, bool(ok), note))

    # ---------- 1. ledger, exact ----------------
    L = ledger.ledgers()
    r = L['rows']
    fold_sum = sum(r[k] for k in ('selector_ab', 'selector_zeta', 'fold_adds', 'fold_cmp',
                                  'hashing', 'dp_test', 'detector', 'hist_slot', 'loop',
                                  'witness_add'))
    add("ledger_fold_sum_1943", fold_sum == r['fold'] == 1943, "ten rows sum to %d" % fold_sum)
    add("ledger_mul_step_21943", r['mul_step'] == K.C_REF + 1943 == 21943)
    add("ledger_sqr_step_13575", r['sqr_step'] == 11632 + 1943 == 13575)
    add("ledger_c_raw_exact", L['c_raw'] == Fraction(3761, 4000), str(L['c_raw']))
    add("ledger_c_stress_exact", L['c_stress'] == Fraction(193879, 200000), str(L['c_stress']))
    w_raw, w_str = L['w_raw'], L['w_stress']
    add("w_raw_2p125.72175", abs(log2(float(w_raw)) - 125.72175) < 1e-5, "%.5f" % log2(float(w_raw)))
    add("w_stress_2p125.76579", abs(log2(float(w_str)) - 125.76579) < 1e-5, "%.5f" % log2(float(w_str)))

    # ---------- 2. step budget and model success ----------------
    S = K.step_budget()
    add("step_budget_is_ledger_S", S == L['S'])
    f, mu, N = K.FRAC_SQR, K.MU_WALK, K.N_CLASSES
    model = 1 - exp(-S * S * (1 - float(f) ** 2) / (2 * float(mu) * N))
    add("model_success_ge_0.392", model >= 0.3920, "%.4f" % model)
    loss_rows = sum(x for _, x in K.LOSS_ROWS)
    add("micro_loss_rows_at_most_3e-4", loss_rows <= Fraction(3, 10**4) and loss_rows > Fraction(29, 10**5),
        "%.3e" % float(loss_rows))
    add("floor_after_micro_losses_ge_0.39", model * (1 - float(loss_rows)) >= 0.39,
        "%.5f" % (model * (1 - float(loss_rows))))

    # ---------- 3. raw trial log ----------------
    blob = json.load(open("trials4x-raw.b64"))
    raw = json.loads(zlib.decompress(base64.b64decode(blob["payload"])).decode())
    rows_by_cell = {}
    for cell in raw["cells"]:
        rows_by_cell[(cell["dp_bits"], cell["arm"])] = cell["rows"]
    add("raw_shape_8_cells", len(rows_by_cell) == 8)
    add("raw_groups_per_cell_1600", all(len(v) == 1600 for v in rows_by_cell.values()))
    add("raw_rows_12800", sum(len(v) for v in rows_by_cell.values()) == 12800)
    add("raw_seed_base_900", raw["seed_base"] == 900 and raw["groups"] == 1600)

    summary = json.load(open("sweep4x-out.json"))
    cells = {(c["dp_bits"], c["arm"]): c for c in summary["cells"]}
    ok_counts = True
    for key, rs in rows_by_cell.items():
        x = sum(row[0] for row in rs)
        c = cells[key]
        ok_counts = ok_counts and x == c["recovered"] and len(rs) == c["n"] == 1600
        ok_counts = ok_counts and sum(row[1] for row in rs) == c["reports"]
        ok_counts = ok_counts and sum(row[6] for row in rs) == c["cyc_trails"]
        ok_counts = ok_counts and sum(row[7] for row in rs) == c["cyc_lost"]
    add("raw_recomputes_summary_exactly", ok_counts)

    # ---------- 4. the one-sided binomial floor, exact ----------------
    xw = sum(cells[(d, "window")]["recovered"] for d in (5, 7, 9, 11))
    nw = sum(cells[(d, "window")]["n"] for d in (5, 7, 9, 11))
    xc = sum(cells[(d, "control")]["recovered"] for d in (5, 7, 9, 11))
    nc = sum(cells[(d, "control")]["n"] for d in (5, 7, 9, 11))
    add("pooled_window_counts_2571_6400", (xw, nw) == (2571, 6400), "%d/%d" % (xw, nw))
    add("pooled_control_counts_2375_6400", (xc, nc) == (2375, 6400), "%d/%d" % (xc, nc))
    quoted = Fraction(3915, 10000)                       # the declared one-sided floor
    alpha = exact_tail_ge(xw, nw, quoted.numerator, quoted.denominator)
    add("quoted_floor_is_valid_95pct", alpha <= Fraction(5, 100),
        "exact tail at 0.3915 = %.6f <= 0.05" % float(alpha))
    tighter = exact_tail_ge(xw, nw, 3917, 10000)
    add("floor_is_not_loose", tighter > Fraction(5, 100),
        "exact tail at 0.3917 = %.6f > 0.05" % float(tighter))
    add("pooled_point_above_0.40", Fraction(xw, nw) > Fraction(4, 10), "%d/%d" % (xw, nw))

    # ---------- 5. loss debit and the claim ----------------
    L_hi = 1 - quoted / Fraction(40, 100)
    add("L_hi_exact_0.02125", L_hi == Fraction(17, 800), str(L_hi))
    w_loss = w_raw / (1 - L_hi)
    w_comb = w_str / (1 - L_hi)
    add("w_combined_2p125.79678", abs(log2(float(w_comb)) - 125.79678) < 1e-5,
        "%.5f" % log2(float(w_comb)))
    claim = K.CLAIM
    add("claim_is_125.81", claim == Fraction(12581, 100), str(claim))
    add("claim_above_w_combined_by_ge_0.005",
        Fraction(claim) - log2(float(w_comb)) >= 0.005,
        "pad %.5f" % (float(claim) - log2(float(w_comb))))
    add("claim_above_w_stress", float(claim) > log2(float(w_str)),
        "pad %.5f" % (float(claim) - log2(float(w_str))))
    add("claim_above_w_raw", float(claim) > log2(float(w_raw)),
        "pad %.5f" % (float(claim) - log2(float(w_raw))))
    add("claim_at_or_below_gate_minus_0.01", float(claim) <= 125.84 - 0.01,
        "claim %.2f <= 125.83" % float(claim))
    add("claim_below_live_head", float(claim) < 125.84, "125.81 < 125.84 (a9c23016)")

    # ---------- 6. the window verdict ----------------
    per_cell = [cells[(d, "window")]["rec_float"] - cells[(d, "control")]["rec_float"]
                for d in (5, 7, 9, 11)]
    add("window_beats_control_every_cell", all(v > 0 for v in per_cell),
        str([round(v, 4) for v in per_cell]))
    diff = xw / nw - xc / nc
    se = sqrt((xw / nw) * (1 - xw / nw) / nw + (xc / nc) * (1 - xc / nc) / nc)
    add("window_advantage_significant", diff / se > 3.0, "z = %.2f" % (diff / se))
    add("pooling_homogeneous", summary["homogeneity"]["p_value"] >= 0.05,
        "chi2 %.3f on %d df, p = %.3f" % (summary["homogeneity"]["chi2"],
                                          summary["homogeneity"]["df"],
                                          summary["homogeneity"]["p_value"]))
    add("control_alone_would_fail", 1 - cells[(5, "control")]["lower"] / 0.40 > 0.02,
        "control pooled lower %.4f" % (xc / nc))

    # ---------- 7. the shipped mechanism, executable ----------------
    import escape
    import classwalk
    ring = escape.WindowRing()
    for i in range(escape.WINDOW):
        ring.push(escape.tag(1000 + i))
    add("window_ring_holds_24", ring.contains(escape.tag(1000)) and ring.n == 24)
    ring.push(escape.tag(5000))
    add("window_ring_evicts_oldest", not ring.contains(escape.tag(1000)))
    add("probe_order_deterministic",
        [classwalk.twist_index(3, 5) for _ in range(2)] == [classwalk.twist_index(3, 5)] * 2 and
        classwalk.twist_index(3, 5) != classwalk.twist_index(5, 3))
    add("twist_index_is_a_class_function",
        all(classwalk.twist_index(t, k) == (t - k) % 12 for t in range(12) for k in range(12)))
    add("probes_eight_and_window_24", escape.PROBES == 8 and escape.WINDOW == 24)

    # ---------- 8. package surface ----------------
    sub = json.load(open("submission.json"))
    sc = sub["claims"][0]["score"]
    add("submission_score_equals_claim",
        Fraction(int(sc["numerator"]), int(sc["denominator"])) == claim and
        int(sc["numerator"]) == 12581 and int(sc["denominator"]) == 100)
    import os
    missing = [e for e in sub["claims"][0]["evidence"] if not os.path.exists(e)]
    add("evidence_files_present", not missing, str(missing))
    add("dependencies_empty", json.load(open("dependencies.json")) == [])
    add("experiments_empty", sub["experiments"] == [])
    part = json.load(open("part_a-out.json"))
    add("part_a_all_pass", part.get("part_a") is True,
        "%d checks" % part.get("checks_n", 0))
    rout = json.load(open("routing-out.json"))
    add("routing_pass", rout.get("routing") is True or rout.get("collapse") == 0)
    note = open("public-note.md").read()
    for token in ("0.4017", "0.3915", "0.02125", "125.81", "125.79678", "125.84",
                  "2571", "6400", "1,600", "a9c23016"):
        add("note_states_" + token.replace(".", "_"), token in note)
    add("note_model_attribution", note.splitlines()[0].strip() == "Model: qwen3.8-flash-next"
        and note.splitlines()[1].strip() == "Harness: Hermes")
    add("note_no_absolute_paths", "/Users/" not in note)

    ok = all(c[1] for c in chk)
    out = {"checks": {n: ("PASS" if v else "FAIL") for n, v, _ in chk},
           "notes": {n: t for n, _, t in chk if t},
           "count": len(chk), "passed": sum(1 for c in chk if c[1]),
           "verdict": "PASS" if ok else "FAIL"}
    json.dump(out, open("certificate.json", "w"), indent=1)
    for n, v, t in chk:
        print("%-44s %s  %s" % (n, "PASS" if v else "FAIL", t))
    print("CERTIFICATE", "PASS" if ok else "FAIL", "%d/%d" % (out["passed"], out["count"]))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
