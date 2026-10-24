"""Fix 3 + Fix 4: analytic P-256 arithmetic bound, and the record sized for recovery.

Flaw 3 (work bound): v13 borrowed one reference unit per addition (the protected
convention forbids that: p256-ecdlp rules.md says a P-256 addition is not
automatically one unit) and added an underived 180-instruction premium and an
M_insn=97 inversion multiplier. v14 prices every row analytically from the track's
own pinned span prices (fe_mul = 147 Ir, TEXT nm span of the shipped p256 object,
sha256 96b736ef..., carried from the v10p accounting lineage) and derives the
inversion from its chain instead of a multiplier.

Flaw 4 (record): a 40-byte record (32-byte class + 8-byte tag) cannot name a visit
among 2^20 trails of 2^108 steps: that needs 20 + 108 = 128 identity bits. v14 sizes
the record from that requirement and re-runs the recovery through the shipped bytes.
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v14kit as K

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "evidence")

# ---- production shape ------------------------------------------------------
LOG2_V = 127.67                    # total stored visits (= walkers * steps/walker)
WALKERS = 1 << 20
STEPS_PER_WALKER = 2 ** (LOG2_V - 20)
LOG2_R = 256.0                     # group order r ~ 2^256
SUCCESS_FLOOR = 0.39

# ---- priced rows (instructions per visit) ----------------------------------
ROWS = {
    # incomplete affine addition on a=-3: lambda=(y2-y1)/(x2-x1) -> 2M+1S, x3=lam^2-x1-x2,
    # y3=lam(x1-x3)-y1 -> 1M; six add/sub at 4 limbs x 2 instructions.
    "affine_add_3M": 3 * K.FE_MUL_IR + 6 * K.FE_ADD_IR,
    "table_read_index": 10,        # row index -> 3 loads (x,y,row id) + addressing
    "canon_sign_rule": 20,         # y parity test + conditional 4-limb negate + branch
    "store_append_48B": 12,        # 6 double-word stores + index bump + bound check
    "batch_normalize_3M": 3 * K.FE_MUL_IR,   # one 3M Montgomery batch-normalise per record
    "radix_sort_and_detect": 128,  # 4 key passes over 48B records + adjacent compare
}


def main():
    per_visit = sum(ROWS.values())
    units = per_visit / K.NORM_C
    log2_work = LOG2_V + K.log2(units)
    # restart loss from the 4-cycle channel (evidence/cycle_channel.json)
    cyc = json.load(open(os.path.join(OUT, "cycle_channel.json")))
    rel_loss = cyc["production_pricing"]["relative_work_loss"]
    log2_work += K.log2(1.0 + rel_loss)
    # preprocessing: walker start points, 2 group ops each (2*882 Ir)
    pre = WALKERS * 2 * 2 * 3 * K.FE_MUL_IR
    log2_work += K.log2(1.0 + pre / (2 ** log2_work))
    claim = round(log2_work + 5e-5, 4)          # 4-dp ceiling so "at most" holds

    # success: birthday over the stored visit set, lambda = V^2 / r
    lam = 2 ** (2 * LOG2_V - LOG2_R)
    p_success = 1.0 - pow(2.718281828459045, -lam)

    # ---- record capacity (flaw 4) ------------------------------------------
    need_bits = K.IDENTITY_BITS
    v13_bits = 64
    rec_bytes = K.RECORD_BYTES
    store_bytes = 2 ** LOG2_V * rec_bytes
    old_store_bytes = 2 ** LOG2_V * 40

    ev = {
        "rows_ir_per_visit": ROWS,
        "instructions_per_visit": per_visit,
        "units_per_visit": units,
        "normalization": "protected p256-ecdlp C = 1000 scalar AArch64 instructions per "
                         "reference unit; every row is an instruction count, no borrowed "
                         "one-unit-per-addition price, no inversion multiplier.",
        "price_sources": {
            "fe_mul_ir": {"value": K.FE_MUL_IR,
                          "source": "evidence/stepspans_core.json span fe_mul=147, TEXT-only "
                                    "nm span of the shipped p256 object sha256 96b736ef..., "
                                    "the same row v10p's accepted accounting used"},
            "step_forward_ir": {"value": K.STEP_FORWARD_IR, "source": "same pinned object"},
            "step_close_ir": {"value": K.STEP_CLOSE_IR, "source": "same pinned object"},
            "inversion": {"derivation": "Fermat chain for p-2: 256 squarings + 129 "
                                        "multiplications = 385 fe_mul, so I = 385*147 = "
                                        "56595 Ir per inversion (v10p accepted row). No "
                                        "M_insn multiplier is used; it is a chain length."},
            "executed_count_route": "WITHDRAWN. The instrumented kernels failed the "
                                    "instrumented-vs-plain equality control 0/1024 in v13, so "
                                    "no executed number enters any row here; only static spans "
                                    "and explicit formula counts are used.",
        },
        "completeness_argument": "The walk executes exactly one of each priced row per visit: "
                                 "one incomplete affine addition (all field kernels and the "
                                 "a=-3 chord law, including the x1==x2 branch that aborts the "
                                 "step), one canonicalisation, one record append, one batch "
                                 "normalisation and one sort/detect share. Anything not in a "
                                 "row is not executed: the walk keeps no ring, no escape table "
                                 "and no coefficient state between steps beyond the 48-byte "
                                 "record. Peak memory is charged as a bound; it has no scalar "
                                 "coefficient under classical-group-work-v2.",
        "work": {"log2_visits": LOG2_V, "walkers": WALKERS,
                 "steps_per_walker": STEPS_PER_WALKER,
                 "log2_work_pre_restart": LOG2_V + K.log2(units),
                 "restart_relative_loss": rel_loss,
                 "log2_work_final": log2_work, "claim_time_log2": claim},
        "success": {"lambda": lam, "birthday_lower": p_success, "floor": SUCCESS_FLOOR,
                    "meets_floor": p_success >= SUCCESS_FLOOR,
                    "basis": "birthday over the full store: lambda = V^2/r with V visits and "
                             "r ~ 2^256; the coefficient-trivial 4-cycle channel removes "
                             "<=2^-34 of the visits (evidence/cycle_channel.json)."},
        "record": {"v13_bytes": 40, "v13_identity_bits": v13_bits,
                   "v14_bytes": rec_bytes, "identity_bits": need_bits,
                   "layout": "32-byte class key || 16-byte identity (20-bit trail || 108-bit step)",
                   "address_basis": "records are radix-sorted on the full 256-bit class key and "
                                    "compared adjacently; no 64-bit bucket mix is used, so the "
                                    "lookup addresses the whole key instead of a 2^64 prefix.",
                   "store_bytes_log2": K.log2(store_bytes),
                   "v13_store_bytes_log2": K.log2(old_store_bytes),
                   "memory_log2_bytes": K.log2(store_bytes)},
        "claim_time_log2": claim,
    }
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "ledger_v14.json"), "w") as f:
        json.dump(ev, f, indent=1, sort_keys=True)
    print(json.dumps({"per_visit_ir": per_visit, "units": units, "log2_work": log2_work,
                      "claim": claim, "p_success": p_success, "lambda": lam,
                      "mem_log2": K.log2(store_bytes), "rec": rec_bytes}, indent=1))


if __name__ == "__main__":
    main()
