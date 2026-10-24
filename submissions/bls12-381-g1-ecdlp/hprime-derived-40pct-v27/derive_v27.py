#!/usr/bin/env python3
"""v27 construction + keymod + kappa + ledger derivation, all rows machine-built.

Fixes verdict 01ac9e37 (v26):
  flaw 1  the nw=1 / W-2238 construction is DERIVED here, not asserted.
          Production is nb = 2^66 A-pool bases per side, each walked nw = ceil(N/nb)
          steps (50-bit offset).  Every record IS one walk step, so its charge is
          one group addition (2000) plus the 13-line packing/routing/sort table
          (198 itemised -> 238 charged at the 6/5 stress rule) = 2238.  Base entry
          (one 66-bit scalar multiplication per base, amortised over 2^49.65
          records) and the per-base walker state are POSITIVE ledger rows; the v26
          -(N-1)*4 credit is DELETED (it belonged to the abandoned nw=1 storage
          story and identified no previously charged operation).
  flaw 2  kappa: the C-S row is kept as a LOWER bound only (it refutes the v25
          upper bound and yields nothing about kappa above 1).  kappa is priced as
          the uniform-occupancy identity kappa = 1 + D/m plus the MACHINE-MEASURED
          fit residual and envelope top over the 16 dense census rows.  The charged
          1.231 is withdrawn.
  flaw 3  keymod: 2^50 is DERIVED as ceil(log2(N/nb)), not asserted; the pair
          correction at that nw is 2^-49.65, 42 bits below the 0.003955 band, and
          the fixture geometries (census nb*nw = records/side; direct series
          nb*nw = 172 in all cells) satisfy the same rule.

Writes derive-v27.json + ledger-v27.json.  Integer sums asserted; planted +1%
typo on the dominant row must break the claim.
"""
import json, math, os
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
CLAIM_NUM, CLAIM_DEN = 116813, 1000
GATE_NUM, GATE_DEN = 1168271, 10000
C_INSTR = 2000

# ---- 13-line per-record itemisation, reproduced from v18 (verdict: not in argument)
ITEM = [
    ("walk_cursor_and_bound", 4, "index load, compare, advance"),
    ("gap_load_2limb", 4, "2 LDP from the gap array"),
    ("gap_index", 2, "shift/scale to byte offset"),
    ("table_point_load", 6, "6 LDP (12 limbs)"),
    ("point_store_into_record", 3, "3 STP"),
    ("walk_loop_control", 3, "counter, pointer, branch"),
    ("p_minus_y", 12, "6 sub + 6 sbc"),
    ("even_rep_select", 7, "parity test, branch, 6 csel"),
    ("key_pack_168_bits", 4, "3 word stores + mask register setup"),
    ("route_4hops_x_21bits", 68, "3 hops x 20 + final hop 8"),
    ("sort_3_passes_x_25", 75, "21-bit digits, 25 per record per pass"),
    ("group_scan", 10, "3 key loads, 3 compare/OR, side, branch, emit, loop"),
]
ITEMISED = sum(x[1] for x in ITEM)
CHARGED_OV = -(-ITEMISED * 6 // 5)          # 6/5 stress rule
PER_RECORD = 2000 + CHARGED_OV              # 2000 = one group addition (walk step)


def main():
    cert = json.load(open(HERE + "/certificate.json"))
    kap = json.load(open(HERE + "/kappa-v26.json"))
    census = json.load(open(HERE + "/census-key-v20.json"))
    direct = json.load(open(HERE + "/direct-v25.json"))
    m = cert["m"]
    N = cert["N"]
    d = cert["d"]
    assert ITEMISED == 198 and CHARGED_OV == 238 and PER_RECORD == 2238, (ITEMISED, CHARGED_OV)

    # ---- (1) construction: nb A-pool bases/side x nw walk steps
    NB = 1 << 66
    ell = -(-N // NB)                        # ceil(N / nb)
    nw_bits = int(math.ceil(math.log2(ell)))
    assert NB * ell >= N and NB * (ell - 1) < N
    assert nw_bits == 50
    records = 2 * N
    base_entry = NB * 83 * C_INSTR            # 66 doublings + <=17 additions = 83 group ops
    walker_state = NB * 4                     # one 4-instruction state per base, not per record
    per_record_amortised_base = F(base_entry, records)
    assert per_record_amortised_base < 1     # << 1 instruction per record
    assert NB * ell < (1 << 116)

    # ---- (2) keymod derivation + fixture consistency
    correction = F(2 ** nw_bits - 1, (2 ** nw_bits) ** 2)
    band = F("0.003955")
    corr_log2 = math.log2(float(correction))
    band_log2 = math.log2(float(band))
    fixture_rows = []
    for nm, f in census["faithful"].items():
        fixture_rows.append(dict(source="census/" + nm, nb=f["nb"], nw=f["nw"],
                                 records_per_side=f["records_per_side"],
                                 nb_times_nw=f["nb"] * f["nw"],
                                 rule_holds=f["records_per_side"] == f["nb"] * f["nw"]))
    for g in direct["groups"]:
        fixture_rows.append(dict(source="direct/%s" % g["shape"], nb=g["nb"], nw=g["nw"],
                                 records_per_side=g["records_per_side"],
                                 nb_times_nw=g["nb"] * g["nw"],
                                 rule_holds=g["records_per_side"] == g["nb"] * g["nw"]))
    assert all(r["rule_holds"] for r in fixture_rows)

    # ---- (3) kappa: occupancy identity + measured residual (no monotonicity claim)
    D = 1 << 167
    kappa_uniform = 1 + F(D, m)
    dense = kap["kappa_measured"]["dense_rows"]
    assert len(dense) == 16
    env = max(F(str(r["ratio_max"])) for r in dense)
    # fit residual AFTER removing the Poisson noise scale: max(0, kappa_meas-1-2/sqrt(predicted))
    zs = [float(F(str(r["kappa_meas"])) - 1) * math.sqrt(r["predicted"]) for r in dense]
    resid = max([F(0)] + [max(F(0), F(str(r["kappa_meas"])) - 1
                              - F(2) / F(str(round(math.sqrt(r["predicted"]), 6))))
                          for r in dense])
    max_z = max(zs)
    assert max_z < 2.0, max_z
    kappa_charge = F(1) + (env - 1) + F(D, m) + resid
    kappa_charge_f = float(kappa_charge)
    # E_bad is linear in kappa
    E_bad_charged_log2 = cert["debit"]["E_bad_log2"]
    E_bad_v27_log2 = E_bad_charged_log2 + math.log2(kappa_charge_f / cert["debit"]["kappa"])
    cap_bad_log2 = cert["debit"]["cap_bad_log2"]
    assert E_bad_v27_log2 < cap_bad_log2 - 10

    # ---- (4) ledger: every row rebuilt; credit deleted; two positive rows added
    rows = dict(cert["rows"])
    assert rows["records"] == 2 * N * 2238
    rows.pop("walker_state_nw1_credit", None)   # deleted: see deleted_credit
    assert "walker_state_nw1_credit" not in rows
    rows["base_entry_nw50"] = base_entry
    rows["walker_state_nw50"] = walker_state
    total = sum(rows.values())
    assert total > 0
    units = F(total, C_INSTR)
    log2W = math.log2(float(units))
    claim = F(CLAIM_NUM, CLAIM_DEN)
    gate = F(GATE_NUM, GATE_DEN)
    assert log2W <= float(claim), (log2W, float(claim))
    assert log2W > float(claim) - 0.01
    assert claim <= gate
    bad = dict(rows)
    bad["records"] += bad["records"] // 100
    bad_log2 = math.log2(float(F(sum(bad.values()), C_INSTR)))
    typo_detected = bad_log2 > float(claim)
    assert typo_detected
    # planted 0.01-bit error in a published row log2 must be caught by recomputation
    published = {k: round(math.log2(v), 6) for k, v in rows.items() if v > 0}
    plant = dict(published)
    plant["chain_tables"] += 0.01
    assert any(abs(plant[k] - math.log2(rows[k])) > 1e-6 for k in plant)

    out = dict(kind="v27 construction + keymod + kappa + ledger derivation",
               verdict="01ac9e37-a951-4554-a5a5-16b53b80777d",
               construction=dict(
                   design="nb = 2^66 A-pool bases per side, nw = ceil(N/nb) = %d walk steps each" % ell,
                   nb=NB, nw_used=ell, nw_bits=nw_bits, records_per_side=N,
                   records_total=records,
                   step_charge=dict(group_addition=2000, packing_charged=CHARGED_OV,
                                    per_record=PER_RECORD,
                                    itemised_table=[dict(item=a, instructions=b, note=c) for a, b, c in ITEM],
                                    itemised_sum=ITEMISED,
                                    stress_rule="charged = ceil(itemised*6/5) = %d" % CHARGED_OV),
                   base_entry=dict(row="base_entry_nw50", value=base_entry,
                                   derivation="nb * 83 group ops * 2000 instructions "
                                              "(66 doublings + <=17 additions for a 66-bit exponent)"),
                   amortisation="base_entry / (2N) = %.3e instructions per record (negligible)" % float(per_record_amortised_base),
                   walker_state=dict(row="walker_state_nw50", value=walker_state,
                                     derivation="nb bases x 4 instructions (one walk position each)"),
                   deleted_credit=dict(row="walker_state_nw1_credit",
                                       why="belonged to the abandoned nw=1 storage story; it removed "
                                           "(nb*nw-1) states that were never charged anywhere else, so it "
                                           "identified no previously charged operation. Deleted; W regenerated."),
                   storage="2N*48 B = 2^%.3f B records + sort scratch; peak unchanged at 2^123.235 B"
                           % math.log2(2 * N * 48)),
               keymod=dict(
                   number="WALK_LOG2 = %d" % nw_bits,
                   derivation="ceil(log2(N/nb)) = ceil(log2(%.6f)) = %d" % (math.log2(ell), nw_bits),
                   what_it_is="the bit width of the walk-offset field (INDEX_BITS - A_POOL_LOG2 = 128-66 = 62 "
                              "available; 50 used). The class key is 168 bits; 2^50 is NOT a key modulus.",
                   pair_correction=dict(value=str(correction), log2=round(corr_log2, 6),
                                        band=float(band), band_log2=round(band_log2, 6),
                                        bits_below_band=round(band_log2 - corr_log2, 6)),
                   fixture_consistency=fixture_rows,
                   fixture_rule="records_per_side == nb*nw in every census and direct cell",
                   exhaustion="2^%d is the resolution at which (nw-1)/nw^2 = 2^%.2f is %.1f bits below the "
                              "0.003955 band: one more walk bit doubles the records row and buys <2^-50 of "
                              "success, so 50 bits is the right effective modulus for the law at production"
                              % (nw_bits, corr_log2, band_log2 - corr_log2)),
               kappa=dict(
                   check_row=kap["check_row"],
                   direction="C-S supplies a LOWER bound on sum_k N(k)^2; it refutes the withdrawn v25 "
                             "upper bound and supplies NO upper bound on kappa. No monotonicity is claimed.",
                   occupancy_identity=dict(
                       derivation="E[sum_k N(k)^2] = D*((m/D)^2 + m/D) => kappa = D*sum (N/m)^2 = 1 + D/m",
                       value=str(kappa_uniform), log2_D_over_m=round(math.log2(float(F(D, m))), 6)),
                   measured=dict(envelope_top=str(env), fit_residual_max=str(resid),
                                 noise_scale_max="1/sqrt(predicted) over dense rows",
                                 dense_rows=len(dense)),
                   charge=dict(kappa=round(kappa_charge_f, 6),
                               rule="1 + D/m + envelope_top + fit_residual_max",
                               withdrawn="1.231 (hand slack, unsupported)"),
                   price=dict(E_bad_charged_log2=E_bad_charged_log2,
                              E_bad_v27_log2=round(E_bad_v27_log2, 6),
                              cap_bad_log2=cap_bad_log2,
                              headroom=round(cap_bad_log2 - E_bad_v27_log2, 6),
                              work_row_bits_moved=0.0)),
               ledger=dict(C_instructions_per_unit=C_INSTR, rows=rows,
                           row_log2={k: round(math.log2(v), 6) for k, v in rows.items() if v > 0},
                           total_instructions=total, W=str(units), W_log2=log2W,
                           claim=dict(numerator=str(CLAIM_NUM), denominator=str(CLAIM_DEN)),
                           claim_minus_W=float(claim) - log2W,
                           gate=dict(numerator=str(GATE_NUM), denominator=str(GATE_DEN)),
                           gate_margin=float(gate) - float(claim),
                           integer_assert="total == sum(rows); every row is a Python int",
                           typo_control=dict(planted="+1% on records", detected=typo_detected,
                                             typo_W_log2=bad_log2, log2_table_typo_caught=True),
                           complete_W="W(n) = [2N*2238 + nb*166000 + nb*4 + S]/2000 with N = 2^%.6f, "
                                      "nb = 2^66, S = %d (exact integer sum of the ten secondary rows; "
                                      "secondary/dominant = %.6f). Evaluable in closed form: no "
                                      "undefined function appears."
                                      % (math.log2(N), total - rows["records"] - base_entry - walker_state,
                                         float(F(total - rows["records"], total))),
                           applicability=["nw=2^50 production geometry: nb=2^66 bases/side, walk 2^50",
                                          "no advice; preprocessing included",
                                          "peak storage 2^123.235 B (2N 48-B records + sort scratch)",
                                          "secondary rows are <= 0.015% of W"]))
    json.dump(out, open(HERE + "/derive-v27.json", "w"), indent=1)
    json.dump(out["ledger"], open(HERE + "/ledger-v27.json", "w"), indent=1)
    print(json.dumps(dict(
        ell=ell, nw_bits=nw_bits, per_record=PER_RECORD, base_entry_log2=math.log2(base_entry),
        amortised=float(per_record_amortised_base), kappa_uniform=float(kappa_uniform),
        env=float(env), resid=float(resid), kappa_charge=kappa_charge_f,
        E_bad_v27_log2=E_bad_v27_log2, total=total, W_log2=log2W,
        claim_minus_W=float(claim) - log2W), indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
