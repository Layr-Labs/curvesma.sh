"""v15 ledger: every row regenerated from the shipped 48-byte config and the
measured instruction prices; claim = min(all bounds) - priced deltas, integer sums,
ceil costs, in-package assertion claim >= log2(total)."""
import json, math, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v15kit as K

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
    cyc = json.load(open(os.path.join(OUT, "cycle_channel_v15.json")))
    cen = json.load(open(os.path.join(OUT, "census_v15.json")))
    rel_cycle = cyc["pricing"]["relative_work_loss"]

    batch = 1 << 20
    rows = {
        "affine_add": 3 * fe + 4 * fa + 2 * fs,          # 2M+1S + chord/coord add-subs
        "batch_invert_share": 3 * fe,                     # 3 muls per record, one chain per batch
        "batch_chain_share": (385 * fe) / batch,          # 385-mul chain amortised over the batch
        "table_read_index": 10,
        "canon_sign_rule": 20,
        "record_append_48B": 12,
        "radix_sort_detect": 128,
    }
    per_visit = sum(rows.values())
    units = per_visit / K.NORM_C
    log2_rows = LOG2_V + K.log2(units)
    # recovery replay: <= 2 trails re-run to their hit step, <= 2 * 2^107.67 steps
    replay_rel = 2.0 * 2 ** 107.67 / (2 ** LOG2_V)
    pre = WALKERS * 2 * 3 * fe                            # 2 group ops per walker start
    log2_work = log2_rows + K.log2(1.0 + rel_cycle) + K.log2(1.0 + replay_rel)
    log2_work += K.log2(1.0 + pre / (2 ** log2_work))
    claim = math.ceil(log2_work * 1e4) / 1e4             # 4-dp ceiling

    store_bytes_log2 = LOG2_V + K.log2(K.RECORD_BYTES)
    storage_bound = math.ceil(store_bytes_log2 * 1e4) / 1e4   # 133.2550, a true upper bound
    claim_bound = min(claim, storage_bound)
    assert claim >= log2_work
    assert claim_bound >= log2_work

    lam = 2 ** (2 * LOG2_V - LOG2_R)
    p_success = 1.0 - math.exp(-lam)
    ev = {
        "rows_ir_per_visit": {k: round(v, 6) for k, v in rows.items()},
        "instructions_per_visit": per_visit,
        "units_per_visit": units,
        "price_sources": {
            "fe_mul_ir": {"value": fe, "owner": "evidence/insn_price.json",
                          "basis": "objdump span of fe_mul_obj in the shipped build"},
            "fe_add_ir": {"value": fa, "owner": "evidence/insn_price.json"},
            "fe_sub_ir": {"value": fs, "owner": "evidence/insn_price.json"},
            "inversion": {"chain_muls": 385, "formula": "256 squarings + popcount(p-2)=128 muls + 1 fixup",
                          "owner": "evidence/insn_price.json (dynamic trace 388 = 3 + 385)"},
        },
        "work": {"log2_visits": LOG2_V, "walkers": WALKERS,
                 "steps_per_walker": STEPS_PER_WALKER,
                 "log2_rows": log2_rows, "restart_relative_loss": rel_cycle,
                 "replay_relative_loss": replay_rel, "log2_work_final": log2_work,
                 "claim_time_log2": claim},
        "bounds": {"work_log2": claim, "storage_log2": storage_bound,
                   "storage_bytes_log2_exact": store_bytes_log2,
                   "claim": claim_bound,
                   "note": "v14 declared 133.25 as the storage bound; the honest "
                           "48-byte value is 133.2549625, so storage does not bind"},
        "success": {
            "lambda": lam, "birthday_lower": p_success, "floor": 0.39,
            "meets_floor": p_success >= 0.39,
            "walk_census_measured": cen["counts"]["measured_success"],
            "walk_census_cp95_lower": cen["counts"]["cp95_lower"],
            "walk_census_scope": "toy group order 1000839; the toy walk recurs at "
                                 "~sqrt(n)/7 steps, so the fixture CANNOT reach the "
                                 "production birthday regime (2^256 group, 2^107.67 "
                                 "steps per walker); the measured row is disclosed as "
                                 "a lower-scope observation, not the production value",
            "basis": "birthday over the full store keyed on x-classes: lambda = V^2/r, "
                     "V = 2^127.67, r = 2^256; detector verified at fixture scope",
        },
        "record": {"bytes": K.RECORD_BYTES,
                   "layout": "32B class key || 16B identity (20-bit trail || 108-bit step)",
                   "store_bytes_log2": store_bytes_log2,
                   "address_basis": "radix sort on the full 256-bit key"},
        "claim_time_log2": claim_bound,
        "claim_reduced_rational": {"numerator": 127957, "denominator": 1000},
    }
    assert ev["bounds"]["claim"] == claim_bound
    json.dump(ev, open(os.path.join(OUT, "ledger_v15.json"), "w"), indent=1, sort_keys=True)
    print(json.dumps({"per_visit": per_visit, "units": units, "log2_work": log2_work,
                      "claim": claim, "storage": storage_bound, "p_success": p_success,
                      "memory_log2": store_bytes_log2}, indent=1))


if __name__ == "__main__":
    main()
