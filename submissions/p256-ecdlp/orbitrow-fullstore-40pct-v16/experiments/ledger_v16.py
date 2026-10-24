"""v16 ledger: regenerated from the 48-byte config, the executed radix sort, the
structural/additive split of the cycle channel, and dependent-replay pricing.

Changes from v15:
  * radix_sort_detect 128 (unexplained constant) -> the EXECUTED LSD radix sort in
    experiments/radix_sort_v16.c: 16 passes x (count 7.33 + scatter 14.0) Ir per
    48-byte record = 341.33 Ir, from the shipped arm64 -O2 build and objdump.
  * replay priced as two DEPENDENT trails: each replayed step pays the full
    extended-Euclid inversion (388 fe_mul = 57,586 Ir) instead of the 2^20-batch
    share (444 Ir), i.e. 40.2x, and only 2 trails are re-run, not 2^20 walkers.
  * cycle channel: structural [A,A,B,B] at R^-2 and additive 4-sum at R^-1, both
    measured at R=32/64/128; only the residual is charged at 1/r.
  * storage: PEAK bound = records 48V + two sort arrays 2x48V = 144V bytes.
"""
import json, math, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v16kit as K

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "evidence")
LOG2_V = 127.67
WALKERS = 1 << 20
STEPS_PER_WALKER = 2 ** (LOG2_V - 20)
LOG2_R = 256.0


def main():
    price = json.load(open(os.path.join(OUT, "insn_price.json")))
    fe = price["prices_ir"]["fe_mul"]
    fa = price["prices_ir"]["fe_add"]
    fs = price["prices_ir"]["fe_sub"]
    cyc = json.load(open(os.path.join(OUT, "cycle_channel_v16.json")))
    cen = json.load(open(os.path.join(OUT, "census_v16.json")))
    rad = json.load(open(os.path.join(OUT, "radix_v16.json")))
    rel_cycle = cyc["pricing"]["relative_work_loss"]

    batch = 1 << 20
    sort_ir = rad["ir_per_record_sort"]
    rows = {
        "affine_add": 3 * fe + 4 * fa + 2 * fs,          # 2M+1S + chord/coord add-subs
        "batch_invert_share": 3 * fe,                     # 3 muls per record, one chain per batch
        "batch_chain_share": (385 * fe) / batch,          # 385-mul chain amortised over the batch
        "table_read_index": 10,
        "canon_sign_rule": 20,
        "record_append_48B": 12,
        "radix_sort_detect": sort_ir,
    }
    per_visit = sum(rows.values())
    units = per_visit / K.NORM_C
    log2_rows = LOG2_V + K.log2(units)

    # dependent replay: two trails, full per-step inversion, no batch amortisation
    full_step_ir = 388 * fe + 4 * fa + 2 * fs + 30     # 388 fe_mul incl. own inverse chain
    bulk_step_ir = rows["affine_add"] + rows["batch_invert_share"] + 30
    infl = full_step_ir / bulk_step_ir
    replay_rel = 2.0 * 2 ** 107.67 * infl / (WALKERS * 2 ** (LOG2_V - 20))
    pre = WALKERS * 2 * 3 * fe
    log2_work = log2_rows + K.log2(1.0 + rel_cycle) + K.log2(1.0 + replay_rel)
    log2_work += K.log2(1.0 + pre / (2 ** log2_work))
    claim = math.ceil(log2_work * 1e4) / 1e4

    # PEAK storage: live records + two full sort arrays (48B each), plus batch queue
    peak_bytes = K.RECORD_BYTES + 2 * K.RECORD_BYTES
    queue_bytes_per_walker = 32
    store_bytes_log2 = LOG2_V + K.log2(peak_bytes)
    queue_log2 = K.log2(WALKERS * queue_bytes_per_walker)
    peak_log2 = K.log2(2 ** store_bytes_log2 + 2 ** queue_log2)
    storage_bound = math.ceil(peak_log2 * 1e4) / 1e4
    claim_bound = min(claim, storage_bound)
    assert claim >= log2_work
    assert claim_bound >= log2_work

    # success: production-shape census (prime order, three walk shapes)
    shapes = cen["shapes"]
    lo_min = min(s["cp95_lower"] for s in shapes)
    p_min = min(s["measured_success"] for s in shapes)
    lam = 2 ** (2 * LOG2_V - LOG2_R)
    p_success = 1.0 - math.exp(-lam)
    ev = {
        "rows_ir_per_visit": {k: round(v, 6) for k, v in rows.items()},
        "instructions_per_visit": per_visit,
        "units_per_visit": units,
        "price_sources": {
            "fe_mul_ir": {"value": fe, "owner": "evidence/insn_price.json"},
            "fe_add_ir": {"value": fa, "owner": "evidence/insn_price.json"},
            "fe_sub_ir": {"value": fs, "owner": "evidence/insn_price.json"},
            "radix_sort": {"ir_per_record": sort_ir, "owner": "evidence/radix_v16.json",
                           "basis": "shipped C, arm64 -O2, objdump inner loops x executed "
                                    "16 passes; 2^20 records sorted, 805306368 bytes moved"},
            "inversion": {"chain_muls": 385, "formula": "256 squarings + popcount(p-2)=128 "
                          "muls + 1 fixup, 388 executed fe_mul"},
        },
        "work": {"log2_visits": LOG2_V, "walkers": WALKERS,
                 "steps_per_walker": STEPS_PER_WALKER,
                 "log2_rows": log2_rows, "restart_relative_loss": rel_cycle,
                 "replay_relative_loss": replay_rel,
                 "replay_inflation_factor": infl,
                 "replay_full_step_ir": round(full_step_ir, 6),
                 "replay_bulk_step_ir": round(bulk_step_ir, 6),
                 "log2_work_final": log2_work, "claim_time_log2": claim},
        "bounds": {"work_log2": claim, "storage_log2_peak": storage_bound,
                   "storage_peak_bytes_per_visit": peak_bytes,
                   "storage_records_log2": store_bytes_log2,
                   "storage_queue_log2": queue_log2,
                   "claim": claim_bound,
                   "note": "peak = 48V live records + 2 x 48V sort arrays + 2^20 x 32B "
                           "batch queue; MAX of live sets, not a sum of averages"},
        "success": {
            "lambda": lam, "birthday_lower": p_success, "floor": 0.39,
            "census_shapes": [{"shape": f"{s['walkers']}x{s['steps_per_walker']}",
                               "trials": s["trials"], "measured": s["measured_success"],
                               "cp95_lower": s["cp95_lower"],
                               "collision_ratio": s["collision_ratio"]} for s in shapes],
            "census_min_cp95_lower": lo_min, "census_min_measured": p_min,
            "meets_floor": p_success >= 0.39 and lo_min >= 0.39,
            "basis": "prime-order fixture (999979) at the production walk shape; the "
                     "collision count matches the birthday prediction, so the measured "
                     "success is the birthday success of this table+selector, not a "
                     "recurrence-limited value",
        },
        "record": {"bytes": K.RECORD_BYTES,
                   "layout": "32B class key || 16B identity (20-bit trail || 108-bit step)",
                   "peak_bytes_per_visit": peak_bytes,
                   "sort": "LSD radix, 16 passes of 4 bits over 4 x 64-bit limbs"},
        "claim_time_log2": claim_bound,
        "claim_reduced_rational": None,
    }
    num = int(round(claim_bound * 10000))
    from math import gcd
    g = gcd(num, 10000)
    ev["claim_reduced_rational"] = {"numerator": num // g, "denominator": 10000 // g}
    assert ev["bounds"]["claim"] == claim_bound
    json.dump(ev, open(os.path.join(OUT, "ledger_v16.json"), "w"), indent=1, sort_keys=True)
    print(json.dumps({"per_visit": per_visit, "units": units, "log2_work": log2_work,
                      "claim": claim, "storage_peak": storage_bound,
                      "replay_rel": replay_rel, "cycle_rel": rel_cycle,
                      "p_success": p_success, "census_min_cp95": lo_min}, indent=1))


if __name__ == "__main__":
    main()
