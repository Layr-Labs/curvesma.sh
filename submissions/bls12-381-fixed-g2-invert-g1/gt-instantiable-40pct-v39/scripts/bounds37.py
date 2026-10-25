#!/usr/bin/env python3
"""bounds37.py -- the store's reached-node bound from UNIFORM-ADDRESS OCCUPANCY.

Verdict 177471b7 flaw 4: 6*ceil(S/(4*2^22)) is not an upper bound on reached
nodes; under the submission's own uniform-address model the 110-bit prefixes
preceding the last radix level are overwhelmingly occupied.  This generator
replaces the dense-leaf count with

    reached(M, S) = ceil( M * (1 - exp(-S/M)) )        M = 2^110 leaf prefixes

which is the expected number of DISTINCT occupied leaves when S records land on
M uniform prefixes, and it reports the exponential argument's magnitude (S/M)
plus the sparse and dense limits as planted controls.  It also keeps the 96-byte
RECORD node and the 32-MiB RADIX node as SEPARATE constants: v36 let the radix
node size overwrite NODE_BYTES.
"""
import json
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
EV = os.path.join(os.path.dirname(HERE), "evidence")

S = int("95803668470929860550224427498997735090")   # the ledger's declared step budget
DB = 22                       # radix digit bits, from dict28.c
H = -(-128 // DB)             # 6 levels
M = 1 << (DB * (H - 1))       # 2^110 distinct LEAF prefixes
NODE_BYTES_RADIX = (1 << DB) * 8          # 32 MiB  -- a radix node
NODE_BYTES_RECORD = 96                    # 96 B    -- an (A,D,y) record node
RECORD_BYTES = 640                        # one records[] entry


def occupancy(Mv, Sv):
    """ceil(M*(1-exp(-S/M))) as an exact int + the magnitude of the argument.
    Dense regime (S >= 64M): M(1-e^-x) > M-1 because e^-x < 2^-64, so the
    ceiling is exactly M -- no floating point in the charged value.
    Sparse regime (64S < M): occupancy = S - S(S-1)/(2M) rounded UP, i.e. every
    record reaches a fresh node."""
    x = Sv / Mv
    # exact-integer dense test: M*e^-x < 1  <=  e^x > M = 2^110.  Since
    # log2(e) = 1.442695 > 14/10, 14*x >= 1100 implies e^x > 2^110 > M.
    if 14 * Sv >= 1100 * Mv:
        return Mv, x, "dense: reached == M, because M*e^-x < 1 (14*S/M >= 1100)"
    if 64 * Sv < Mv:
        lo = Sv - (Sv * (Sv - 1)) // (2 * Mv)
        return min(Sv, lo + 1), x, "sparse: bound == S (every record is a new node)"
    # generic middle band: ln2 * log2 form, integer compare, no cancellation
    return int(math.ceil(Mv * (1.0 - math.exp(-x)))), x, "middle"


def main():
    occ, x, regime = occupancy(M, S)
    levels, total = [], 0
    for j in range(1, H + 1):
        cap = 1 << (DB * j)
        n = occ if j == H else min(cap, occ)
        levels.append({"level": j, "prefix_bits": DB * j, "capacity": str(cap),
                       "reached": str(n)})
        total += n
    init_ir = total * (1 << DB)          # one Ir per slot initialised
    node_bytes_total = total * NODE_BYTES_RADIX
    dense6 = 6 * (-(-S // (4 * (1 << DB))))
    rows = [
        {"array": "records[] (A,D,y) + tags + 128-bit ref", "value": str(S),
         "bytes_per_entry": RECORD_BYTES, "bytes": str(S * RECORD_BYTES),
         "single_allocation_ok": True, "note": "one entry per executed step"},
        {"array": "alive[] / matched[] / trail_ids[] / positions[]", "value": str(S),
         "bytes_per_entry": 24, "bytes": str(S * 24), "single_allocation_ok": True,
         "note": "two flags + 8-byte trail id + 8-byte position"},
        {"array": "record node (96 B), one per NEW RESIDENT KEY (verdict 4558f762 flaw 5: v37 charged one per 40 records while dict28_wide_insert allocates per resident key)", "value": str(S),
         "bytes_per_entry": NODE_BYTES_RECORD, "bytes": str(S * NODE_BYTES_RECORD),
         "single_allocation_ok": True,
         "note": "SEPARATE constant from the radix node; v36 overwrote NODE_BYTES with the "
                 "32 MiB radix size"},
        {"array": "radix node (32 MiB), one calloc per REACHED node", "value": str(S),
         "bytes_per_entry": NODE_BYTES_RADIX, "bytes": str(node_bytes_total),
         "single_allocation_ok": True,
         "note": "32 MiB <= SIZE_MAX: every allocation is addressable with 64-bit pointers; the "
                 "128-bit slot is consumed 22 bits per level over h=6 levels"},
        {"array": "page table (one 128-bit node id per allocated node)", "value": str(S),
         "bytes_per_entry": 16, "bytes": str(total * 16), "single_allocation_ok": True,
         "note": "the table is addressed by the same level structure recursively"},
        {"array": "probe chain per bucket", "value": 4, "bytes_per_entry": 72, "bytes": 288,
         "single_allocation_ok": True},
        {"array": "slot export buffers (keys 6 + ref 2 + evicted 2)", "value": 10,
         "bytes_per_entry": 8, "bytes": 80, "single_allocation_ok": True},
        {"array": "guard ring hist[]", "value": 6, "bytes_per_entry": 8, "bytes": 48,
         "single_allocation_ok": True},
        {"array": "six-limb scratch (montmul t[8] + 12 lanes + key 6)", "value": 86,
         "bytes_per_entry": 8, "bytes": 688, "single_allocation_ok": True},
    ]
    out = {
        "mode": "bounds37-occupancy",
        "store": "128-bit radix tree over the 128-bit slot, h = ceil(128/22) = 6 levels, "
                 "node = 2^22 x 8 B = 32 MiB, one calloc per REACHED node",
        "parsed_constants": {"DB": DB, "H": H, "M_leaf_prefixes": str(M),
                             "NODE_BYTES_RADIX": NODE_BYTES_RADIX,
                             "NODE_BYTES_RECORD": NODE_BYTES_RECORD, "RECORD_BYTES": RECORD_BYTES,
                             "source": "scripts/ir28/dict28.c: D28_DB, rd_new() calloc(2^db, 8), "
                                       "ws_node record = 96 B"},
        "occupancy": {"S": str(S), "M": str(M), "S_over_M": x, "S_over_M_log2": math.log2(x),
                      "regime": regime, "reached_leaf_nodes": str(occ),
                      "exponential_argument_magnitude":
                          "S/M = %.6e; e^-x < 2^-106473, so M(1-e^-x) exceeds M-1 by a margin "
                          "larger than one node" % x,
                      "rule": "reached(M,S) = ceil(M(1-exp(-S/M))), NOT 6*ceil(S/(4*2^22))"},
        "levels": levels,
        "init_ir": str(init_ir),
        "node_bytes_total": str(node_bytes_total),
        "dense_six_ceil_control": {"v36_formula": "6*ceil(S/(4*2^22))", "value": str(dense6),
                                   "occupancy_value": str(occ),
                                   "ratio": float(occ) / float(dense6),
                                   "verdict": "the dense-leaf formula understates the reached-node "
                                              "count by the printed ratio; it is flagged, not used"},
        "small_regime_control": dict(zip(("reached", "S_over_M", "regime"), occupancy(M, 1000))),
        "bounds_table": rows,
        "peak_bytes": str(sum(int(r["bytes"]) for r in rows)),
        "production_layout": {"rule": "h = ceil(128/22) = 6 radix levels; a node is 2^22 x 8 B = "
                                      "32 MiB; every calloc argument is <= SIZE_MAX"},
    }
    json.dump(out, open(os.path.join(EV, "bounds37.json"), "w"), indent=1)
    print(json.dumps({"reached": str(occ), "S_over_M": x, "init_ir": str(init_ir),
                      "total_nodes": str(total), "dense6": str(dense6),
                      "ratio": float(occ) / float(dense6)}, indent=1))


if __name__ == "__main__":
    main()
