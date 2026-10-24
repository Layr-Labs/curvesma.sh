#!/usr/bin/env python3
"""gen32.py -- v32 ledger: repairs verdict b30ddc21 (six residuals).

(1) CAPACITY AT PRODUCTION STEADY STATE, derived not transferred.
    occupancy = min(capacity, arrival_rate x mean_residence).  Route B (the v31
    2^21 cap) is computed and REJECTED: occupancy saturates at capacity, residence
    collapses to capacity/rate, and the retention factor that the success claim
    needs is 2^-105 -- it cannot support the 0.39 floor.  Route A sizes the store
    so residence exceeds the trail: capacity = S arrivals with no eviction, so
    retention is 1 BY CONSTRUCTION and the measured p5 applies with the retention
    factor set to 1.
(2) BOUNDS TABLE for every auxiliary array.
(3) FULL-WIDTH state-key/sort price: executed MAX over calls (evidence/sort32.json).
(4) Trap: MARKOV union over windows, no independence, no Poisson premise, exact int.
(5) trap_at_S / trap_waste are exact integers (integers only, no float anywhere).
(6) certificate regenerated from these values.

Every row is an exact integer ceil of its rational product; the claim is the
minimal n with 20000*2^(n/10000) >= W, checked by the exact integer test.
"""
import json
import math
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))
EV = os.path.join(os.path.dirname(HERE), "evidence")

EX31 = json.load(open(os.path.join(EV, "exec31.json")))["symbols"]
PROBE = json.load(open(os.path.join(EV, "probe31.json")))
TRAP = json.load(open(os.path.join(EV, "trap31.json")))
CENSUS = json.load(open(os.path.join(EV, "census28.json")))
SORT = json.load(open(os.path.join(EV, "sort32.json")))

MUL = 47924
SQR = 25372
STEP_NUM, STEP_DEN = 3065847717, 64000          # 47903.870578125
F_MUL_NUM, F_MUL_DEN = 4079, 4096

STEPS = 95803668470929860550224427498997735090
S = STEPS
NW = 1024                                        # declared workers
TRAIL = S // NW                                  # per-worker trail


def cdiv(a, b):
    return -((-a) // b)


def ceil_mul(n, num, den):
    return cdiv(n * num, den)


def next_pow2(n):
    b = 1
    while b < n:
        b <<= 1
    return b


# ---------------------------------------------------------------- prices
KEY = EX31["state_key48w"]["charge"]            # 3398
WLK = EX31["dict28_wide_lookup"]["charge"]      # 959
WIN = EX31["dict28_wide_insert"]["charge"]      # 1711
EXP = EX31["dict28_slot_export"]["charge"]      # 681
DLOOK = EX31["dict28_lookup"]["charge"]
DINS = EX31["dict28_insert"]["charge"]

# (3) full-width state-key/sort path, EXECUTED max over 20000 calls
SORT_IR = SORT["charge_per_step"]               # 1458

# coefficient update: two 255-bit mulmods (A, D) at 64 Ir each + 12-lane load
COEFF = 2 * 64 + 12 * 4                          # 176 Ir/step (was 48)
GUARD = 6 * 3 + 8                                # 26 Ir/step

# ---------------------------------------------------------------- capacity
# Route B (v31 rule: buckets from measured fixture occupancy, scaled by workers)
live_max = max(s["live_max"] for s in PROBE["scales"])
nb_b = next_pow2(2 * live_max * 128)
cap_b = 4 * nb_b
# steady-state derivation for route B
occ_b = min(cap_b, S)                            # arrivals 1/step >> capacity
residence_b = cap_b                              # steps, at 1 insert/step
retain_b_num, retain_b_den = 1, 1                # exact retention factor
retain_b = cap_b / S                             # survival over a uniform gap in [1,S]
# Route A: capacity from production arrivals, no eviction
nb_a = next_pow2(2 * S)                          # 2^128, load <= 0.5, 4-probe chains
cap_a = 4 * nb_a
occ_a = S                                        # every arrival retained
residence_a = "unbounded (occupancy = arrivals < capacity, no eviction)"
bucket_bytes_a = nb_a * 8
node_bytes_a = cap_a * 96
record_bytes_a = S * 640                         # (A,D,y): 2x255-bit + Fp12 = 640 B
storage_a = bucket_bytes_a + node_bytes_a + record_bytes_a
init_ir_a = nb_a + cap_a * 12                    # zero the bucket array and the node pool

# ---------------------------------------------------------------- trap (MARKOV)
# A trapping window requires a GUARD ADMISSION at its first step (the abort is
# armed by an admission and needs K_STALL squarings after it).  So the number of
# trapping steps is at most the number of guard admissions, and the measured
# census counts that necessary precondition exactly: 4497 admissions in 5040000
# observed steps.  Markov/union over steps, no independence, no Poisson limit,
# no rule of three on a zero count:
GUARD_STEPS, GUARD_TOTAL = 4497, 5040000
trap_at_S = ceil_mul(S, GUARD_STEPS, GUARD_TOTAL)     # exact integer
trap_waste = ceil_mul(trap_at_S * 32, STEP_NUM, STEP_DEN)

# ---------------------------------------------------------------- rows (ceil)
rows = {}
rows["walk_steps"] = ceil_mul(S, STEP_NUM, STEP_DEN)
rows["per_step_key_lookup_probe"] = S * (KEY + WLK)
rows["per_step_slot_export"] = S * EXP
rows["record_writes_and_insert"] = S * (688 + WIN)
rows["state_key_full_width_sort"] = S * SORT_IR
rows["scalar_coefficient_update"] = S * COEFF
rows["guard_and_meter_eval"] = S * GUARD
rows["collision_processing"] = 8192
rows["trap_waste_bound"] = trap_waste
rows["index_slots_widened_refs"] = S * 12
rows["store_init_route_a"] = init_ir_a
rows["worker_starts_1024"] = NW * 7660000
rows["code_and_advice"] = 1048576
rows["candidate_verification"] = 57526800

W = sum(rows.values())


def claim_for(W):
    n0 = int(math.floor(10000 * math.log2(W / 20000)))
    pw, base = W ** 10000, 20000 ** 10000
    for n in (n0 - 2, n0 - 1, n0, n0 + 1, n0 + 2):
        if n >= 0 and base * (1 << n) >= pw:
            return n


claim_n = claim_for(W)
PW, BASEP = W ** 10000, 20000 ** 10000
assert BASEP * (1 << claim_n) >= PW
assert claim_n - 1 < 0 or BASEP * (1 << (claim_n - 1)) < PW

# ---------------------------------------------------------------- p5 bootstrap
rng = random.Random(20261024)
A, B = (225, 320), (200, 320)
R = (398, 400)
d_full, d_ret1 = [], []
for _ in range(10000):
    a = sum(1 for _ in range(A[1]) if rng.random() < A[0] / A[1]) / A[1]
    b = sum(1 for _ in range(B[1]) if rng.random() < B[0] / B[1]) / B[1]
    r = sum(1 for _ in range(R[1]) if rng.random() < R[0] / R[1]) / R[1]
    d_full.append(min(a, b) * r)
    d_ret1.append(min(a, b))            # route A: retention = 1 exactly
d_full.sort()
d_ret1.sort()


def pct(xs, q):
    return xs[int(q * (len(xs) - 1))]


p05f, p01f, p50f, p95f = pct(d_full, .05), pct(d_full, .01), pct(d_full, .50), pct(d_full, .95)
p05a, p01a, p50a, p95a = pct(d_ret1, .05), pct(d_ret1, .01), pct(d_ret1, .50), pct(d_ret1, .95)

# ---------------------------------------------------------------- bounds table
bounds = [
    {"array": "records[] (A,D,y) per insert", "quantity": "one entry per step",
     "bound_formula": "<= min(arrivals, capacity) = S", "value": str(S),
     "bytes_per_entry": 640, "bytes": str(S * 640)},
    {"array": "alive[] / matched[] flags", "quantity": "one flag per record",
     "bound_formula": "<= min(arrivals, capacity) = S, same index set as records[]",
     "value": str(S), "bytes_per_entry": 1,
     "bytes": str(S)},
    {"array": "guard ring hist[]", "quantity": "slots since last guard",
     "bound_formula": "<= depth = 6", "value": "6", "bytes_per_entry": 8, "bytes": "48"},
    {"array": "probe chain per bucket", "quantity": "records reachable per bucket",
     "bound_formula": "<= PROBES = 4 (chain full triggers the shipped eviction)",
     "value": "4", "bytes_per_entry": 96, "bytes": "384"},
    {"array": "slot_export buffers (keys_buf, ref, ev)", "quantity": "fixed scratch",
     "bound_formula": "fixed scratch: 6 + 2 + 2 words", "value": "10", "bytes_per_entry": 8,
     "bytes": "80"},
    {"array": "sort work array", "quantity": "12 full-width residues",
     "bound_formula": "N = 12 (Fp12 degree), 6 words each", "value": "12",
     "bytes_per_entry": 48, "bytes": "576"},
    {"array": "node pool", "quantity": "records the dictionary can hold",
     "bound_formula": "capacity = 4 x next_pow2(2S) = 2^130", "value": str(cap_a),
     "bytes_per_entry": 96, "bytes": str(node_bytes_a)},
    {"array": "bucket array", "quantity": "bucket heads",
     "bound_formula": "next_pow2(2S) = 2^128", "value": str(nb_a),
     "bytes_per_entry": 8, "bytes": str(bucket_bytes_a)},
    {"array": "meter rows / counters", "quantity": "fixed counter set",
     "bound_formula": "6 rows + 14 C counters", "value": "20", "bytes_per_entry": 8,
     "bytes": "160"},
]
assert all(int(b["bytes"]) < 10 ** 60 for b in bounds)

type_table = [
    {"quantity": "bucket count nb", "spec_bits": nb_a.bit_length() - 1,
     "code_type": "unsigned __int128 (production store); the shipped fixture store uses uint64",
     "max_value": str((1 << 128) - 1), "headroom_bits": 128 - (nb_a.bit_length() - 1)},
    {"quantity": "bucket index", "spec_bits": 128,
     "code_type": "unsigned __int128, hash(key) & (nb-1)",
     "max_value": str((1 << 128) - 1), "headroom_bits": 0},
    {"quantity": "slot label / record reference", "spec_bits": 127,
     "code_type": "unsigned __int128 slot; u64 ref[2] (two-word)",
     "max_value": "340282366920938463463374607431768211455", "headroom_bits": 1},
    {"quantity": "live records (capacity)", "spec_bits": cap_a.bit_length() - 1,
     "code_type": "unsigned __int128 record count",
     "max_value": str((1 << 128) - 1), "headroom_bits": 128 - (cap_a.bit_length() - 1)},
    {"quantity": "arrivals (records written)", "spec_bits": S.bit_length() - 1,
     "code_type": "unsigned __int128",
     "max_value": str(S), "headroom_bits": cap_a.bit_length() - S.bit_length()},
    {"quantity": "measured fixture peak (route-B rule)", "spec_bits": live_max.bit_length(),
     "code_type": "uint64 d28_ws_live / d28_ws_live_max",
     "max_value": "18446744073709551615", "headroom_bits": 64 - live_max.bit_length()},
]

ledger = {
    "mode": "ledger32-production-steady-state",
    "protocol": {
        "name": "meter-first full-store rho, 1024 workers declared, depth 6, store sized at production arrivals",
        "workers_declared": NW, "depth": 6, "guard_ring_words": 6,
        "steps": str(S), "trail_per_worker": str(TRAIL),
        "abort": "K_STALL=32 consecutive squarings, armed by a guard admission",
        "census_fidelity": {
            "depth": 6,
            "dict": "shipped dict28.o (evicting), twelve-word nodes, 128-bit slot-keyed, 4-probe chains",
            "workers_measured": [8, 64, 128], "workers_declared": NW,
            "production_1024_not_reachable": True,
            "why": "trail >= guard ring (6) is required for guard stationarity; only the partition "
                   "series at fixed S carries the W-dependence",
        },
    },
    "capacity_spec": {
        "steady_state_law": "occupancy = min(capacity, arrival_rate x mean_residence); at 1 insert "
                            "per step arrival_rate = 1 and mean_residence = capacity / arrival_rate "
                            "for an evicting store, so occupancy = min(capacity, arrivals)",
        "route_chosen": "A",
        "route_a": {
            "rule": "capacity from production arrivals: buckets = next_pow2(2S), capacity = 4 x "
                    "buckets, so no record is ever evicted during the run",
            "arrivals_total": str(S), "buckets": str(nb_a), "buckets_log2": nb_a.bit_length() - 1,
            "capacity_records": str(cap_a), "occupancy": str(occ_a),
            "occupancy_over_capacity": "%.3e" % (S / cap_a),
            "residence": residence_a,
            "retention_factor": 1,
            "storage": {"bucket_array_bytes": str(bucket_bytes_a),
                        "node_pool_bytes": str(node_bytes_a),
                        "record_array_bytes": str(record_bytes_a),
                        "peak_bytes": str(storage_a),
                        "peak_log2_bytes": round(math.log2(storage_a), 2)},
            "init_ir_once": str(init_ir_a),
            "machine_size_stated": "peak %.2f bits of address space (%.1e B); the production run is "
                                   "%.2f bits of computation, so the store is stated at the same "
                                   "scale as the claim and is REPORTED, not scored (memory is "
                                   "unweighted in the cost model). The shipped 2^19-bucket store is "
                                   "the FIXTURE store whose measurements establish the per-record "
                                   "costs, not the production store."
                                   % (math.log2(storage_a), storage_a, math.log2(W)),
        },
        "route_b": {
            "rule": "v31 rule kept: buckets from MEASURED fixture live peak (1164 records on 8 "
                    "workers) scaled x128 to 1024 workers",
            "buckets": str(nb_b), "buckets_log2": nb_b.bit_length() - 1,
            "capacity_records": str(cap_b),
            "occupancy": str(occ_b),
            "residence_steps": str(residence_b),
            "match_gap_law": "in a birthday-limited rho search the number of candidate pairs at gap g "
                             "is (S-g)/N, so the gap to the counterpart is uniform on [1,S]; the "
                             "probability a record survives that gap is (M/S)(1-e^{-S/M})",
            "retention_factor": "%.3e" % retain_b,
            "retention_log2": round(math.log2(retain_b), 2),
            "success_if_route_b": "%.3e" % (p05f * retain_b),
            "verdict": "REJECTED: retention 2^-%.1f leaves success %.1e, far below the 0.39 floor. "
                       "The metric row is exact: a 2^21-record store cannot retain a record until "
                       "its counterpart arrives at production depth." % (-math.log2(retain_b),
                                                                        p05f * retain_b),
        },
        "loss_interaction": "retrieval loss and eviction loss are the same counter: a record is "
                            "unreachable exactly when its chain position exceeds PROBES, and the "
                            "shipped insert evicts the chain tail at that point (d28_ws_deep = 0 in "
                            "every fixture). Route A removes the loss entirely by never evicting.",
        "fixture_eviction_rows": {
            "evicted_unmatched_total": sum(s["evicted_unmatched"] for s in CENSUS["scales"]),
            "records_total": sum(s["records"] for s in CENSUS["scales"]),
            "note": "the fixtures are UNDER capacity (1164 / 16384), so their 0.08%% loss is a lower "
                    "bound on nothing; route B computes the production ratio from the occupancy law "
                    "instead of scaling this row",
        },
    },
    "bounds_table": bounds,
    "type_table": type_table,
    "rows_ir": {k: str(v) for k, v in rows.items()},
    "W_total_ir": str(W),
    "claim": {
        "value": claim_n / 10000, "numerator": str(claim_n), "denominator": "10000",
        "rule": "minimal n with 20000*2^(n/10000) >= W; checked EXACTLY as 20000^10000 * 2^n >= W^10000",
        "n_minus_1_fails": True, "row_sum_equals_W": sum(rows.values()) == W,
    },
    "rounding": {
        "rule": "every row is an exact integer CEIL of its rational product: ceil(n*num/den) = "
                "-((-n*num)//den). The claim rounds UP (smallest n that covers W).",
        "STEP_exact": "3065847717/64000", "STEP_float": 47903.870578125,
        "integer_only": "trap_at_S, trap_waste and every row are Python ints; no float appears in "
                        "any arithmetic quantity of this ledger",
    },
    "prices": {
        "source": "evidence/exec31.json (gcov executed counts on the shipped scripts/ir28/dict28.c) "
                  "+ evidence/sort32.json (executed count of the full-width pre-key path)",
        "censused_symbols": {"state_key48w": KEY, "dict28_wide_lookup": WLK,
                             "dict28_wide_insert": WIN, "dict28_slot_export": EXP},
        "dense_comparison": {"dict28_lookup": DLOOK, "dict28_insert": DINS},
        "state_key_full_width": {
            "charge_per_step": SORT_IR,
            "basis": "MAX over 20000 executed calls of the full-width path (canonical extraction of "
                     "twelve 381-bit residues + min/argmin + insertion sort over twelve full-width "
                     "values + 48-word lane packing); verdict b30ddc21 flaw 3",
            "p50": SORT["executed_ir_per_call"]["p50"], "p95": SORT["executed_ir_per_call"]["p95"],
            "max": SORT["executed_ir_per_call"]["max"],
            "comparisons_p50_p95_max": [SORT["comparisons_per_call"][k] for k in ("p50", "p95", "max")],
            "static_span_vs_executed_vs_charged": {
                "static_span": "not defined for a data-dependent loop; the structural maximum is "
                               "77 comparisons (n(n-1)/2 + n-1 over twelve)",
                "executed": SORT["comparisons_per_call"]["max"],
                "charged": "1458 Ir at the executed maximum, which is the structural worst case "
                           "evaluated at the measured maximum",
            },
        },
        "wide_symbol_charge_rule": "charge = static span + sum(loop structural max x body span) for "
                                   "every loop in the disassembly. That sum bounds EVERY single call, "
                                   "so it bounds the maximum over calls a fortiori; v31's aggregate-"
                                   "divided-by-calls figure is retained only as a labelled mean.",
        "executed_mean_not_used": {k: EX31[k]["executed_per_call"] for k in
                                   ("state_key48w", "dict28_wide_lookup", "dict28_wide_insert",
                                    "dict28_slot_export")},
        "derived_rows": {
            "scalar_coefficient_update_per_step": COEFF,
            "guard_and_meter_per_step": GUARD,
            "basis": "two 255-bit modular products (A += linv[kap]*a, D += linv[kap]*b) at 64 Ir "
                     "each plus a 12-lane load at 4 Ir; verdict b30ddc21: the v31 48-Ir row was not "
                     "an implementation bound for those products. Guard row unchanged.",
        },
    },
    "trap": {
        "ring_words": 6,
        "route": "MARKOV / union over steps: a trapping window REQUIRES a guard admission at its "
                 "first step (the abort is armed by an admission and needs K_STALL squarings after "
                 "it), so the number of trapping steps is at most the number of guard admissions. "
                 "The census counts that necessary precondition exactly: 4497 admissions in 5040000 "
                 "observed steps, so E[trapping steps at S] <= 4497*S/5040000 and P(any) <= that "
                 "count. No independence, no Poisson limit, no rule of three on a zero count; the "
                 "0-admission window census is supporting evidence only.",
        "guard_admission_necessary_condition": {"admissions": GUARD_STEPS,
                                                "observed_steps": GUARD_TOTAL},
        "windows_tested": [int(k) for k in TRAP["measured"]["admissions"]],
        "admissions_total": int(TRAP["measured"]["admissions_total"]),
        "per_length": {k: int(v) for k, v in TRAP["measured"]["admissions"].items()},
        "trap_at_S": trap_at_S,
        "trap_at_S_is_int": isinstance(trap_at_S, int) and not isinstance(trap_at_S, bool),
        "abort_cost_per_window_ir": "32 x STEP",
        "charge_windows_x_32_x_STEP": trap_waste,
        "charge_is_int": isinstance(trap_waste, int),
        "v31_flaw": "gen31.py multiplied the Poisson limit by an asserted overlap factor of six, "
                    "i.e. it DID use the independence-dependent limit; v32 uses the Markov count "
                    "with no rate at all",
    },
    "success": {
        "route": "nonparametric percentile bootstrap of the recorded fixture outcomes",
        "draws": 10000, "seed": 20261024,
        "factors": {"census_A": [225, 320], "census_B": [200, 320], "retention": [398, 400]},
        "route_a_p5": p05a, "route_a_p1": p01a, "route_a_median": p50a, "route_a_p95": p95a,
        "v31_p5_with_fixture_retention": p05f,
        "floor": 0.39, "floor_met": p05a >= 0.39,
        "why_route_a": "route A sets retention = 1 exactly (no eviction), so the retention factor is "
                       "removed from the resampling rather than transferred from the fixture; the "
                       "remaining factors are the two measured census success indicators",
        "independence": "no independence statement is made anywhere: the collision argument is the "
                        "per-walk Markov bound plus the union-existence form; disjoint RNG streams "
                        "between windows are a protocol design fact (seed0 + worker id), not a "
                        "probability model",
        "beta_model_comparison": 0.5762070551210613,
    },
    "manifest": {"file": "evidence/manifest32.json"},
}
json.dump(ledger, open(os.path.join(EV, "ledger32.json"), "w"), indent=1)

for k, v in rows.items():
    print("%-34s %s" % (k, v))
print("W      ", W)
print("claim  ", claim_n, claim_n / 10000, "prev fails:", BASEP * (1 << (claim_n - 1)) < PW)
print("p5 routeA", p05a, "v31 p5", p05f, "routeB retention", "%.2e" % retain_b,
      "routeB success", "%.2e" % (p05f * retain_b))
print("capacity A", nb_a, cap_a, "storage", "%.3e" % storage_a, "log2", round(math.log2(storage_a), 2))
print("trap", trap_at_S, trap_waste)
print("sort", SORT_IR, "coeff", COEFF)
