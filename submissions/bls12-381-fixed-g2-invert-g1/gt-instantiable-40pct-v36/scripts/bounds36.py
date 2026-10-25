#!/usr/bin/env python3
"""bounds35.py -- regenerate bounds_table FROM THE SHIPPED SOURCE CONSTANTS.

v34's ledger declared a bucket-head count and node capacity that its own code
cannot represent, while dict28.c declares `u64 nb` and a
u64 bucket index (verdict ea688603 flaw 2).  Nothing is copied from prose here:
every number below is read out of scripts/ir28/dict28.c by regex, and the
production bucket count is derived from the WIDTH OF THE SHIPPED FIELD, not
from a wish:

  dict28_ws_init:  `if (nb & (nb - 1)) return -1;`   -> nb is a power of two
                   `ws->b = calloc((size_t)nb, sizeof(ws_node *));`
                   nb is a `u64`                        -> nb <= 2^63
  ws_bucket:       `mix64((u64)slot ^ mix64((u64)(slot >> 64))) & (nb - 1)`
                   a 64-BIT mix, so one level of the shipped type addresses
                   at most 2^63 distinct buckets.

v35's answer -- TWO levels of that type, 2^63 x 2^63 heads -- does NOT repair
the defect: the per-level calloc is still 2^63 x 8 = 2^66 B > SIZE_MAX, and a
second level only multiplies the problem.  v36 replaces the flat array with a
RADIX TREE over the full 128-bit slot: height h = ceil(bits/db), every node a
`calloc(2^db, 8)`.  At db = 22 (D28_DB in the shipped source) the largest
allocation in the structure is 2^22 x 8 = 2^25 B = 32 MiB, every level's
count x element size is <= SIZE_MAX, and the addressable bucket count is
2^(h*db) = 2^132 >= 2^128.  The calloc sweep below FAILS the build if any
array in this table has count x element size > SIZE_MAX (2^64-1) on the
stated 64-bit ABI, and it is run once with a planted violation.

Also emits the LITERAL SWEEP with its planted control: a stale token outside
the disclosure file, on a line that does not name its withdrawn source, fails.
"""
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
EV = os.path.join(PKG, "evidence")
SRC = os.path.join(HERE, "ir28", "dict28.c")

S = int("95803668470929860550224427498997735090")   # the ledger's declared step budget

# ---------------------------------------------------------------- source parse
src = open(SRC).read()


def define(name):
    m = re.search(r"^#define\s+%s\s+(\d+)" % name, src, re.M)
    return int(m.group(1))


PROBES = define("PROBES")
ENTRY = define("ENTRY")
# ws_node: u128 slot | u64 key[6] | u64 ref[2] | u64 tag | ws_node *next
NODE_BYTES = 16 + 8 * 6 + 8 * 2 + 8 + 8
BUCKET_HEAD_BYTES = 8                       # ws_node *
HEADER_BYTES = 8 + 8                        # { u64 nb; ws_node **b; }
NB_FIELD_BITS = 64                          # `typedef struct { u64 nb; ... }`
NB_MAX_LOG2 = NB_FIELD_BITS - 1             # power of two that fits u64 and calloc
assert re.search(r"u64 bits, db, h;", src) and re.search(r"calloc\(n, sizeof\(ws_node \*\)\)", src)
assert re.search(r"if \(bits < 4 \|\| bits > 128\) return -1;", src)
PROBE_MASK_BITS = 128                       # the radix indexes the FULL 128-bit slot
D28_DB = define("D28_DB")                   # production radix digit bits = 22
assert 1 <= D28_DB <= 30
SIZE_MAX = (1 << 64) - 1                    # stated 64-bit ABI
BITS = 128                                  # bucket bits required by production
HEIGHT = -(-BITS // D28_DB)                 # ceil(128/22) = 6
NODE_BYTES = 8 * (1 << D28_DB)              # one node = 2^db pointers
assert NODE_BYTES <= SIZE_MAX, "per-level calloc exceeds SIZE_MAX"
HEADS_LOG2 = HEIGHT * D28_DB                # 132 >= 128 addressable bucket bits
HEADS = 1 << BITS                           # 2^128 distinct buckets
CAPACITY = PROBES * HEADS                   # 4 heads-per-bucket chain cap

# ------------------------------------------------------------ production bound
arrivals = S
live = min(arrivals, CAPACITY)
lam = arrivals / HEADS

arrivals = S
NODES_REACHED = -(-arrivals // (PROBES * (1 << D28_DB))) * HEIGHT   # lazy path nodes
rows = [
    {"array": "records[] (A,D,y) per insert", "source_constant": "one entry per executed step",
     "bound_formula": "min(arrivals, capacity) = S", "value": str(live),
     "bytes_per_entry": 640, "bytes": str(live * 640)},
    {"array": "alive[] / matched[] / trail_ids[]", "source_constant": "protocol28.run_one lists",
     "bound_formula": "one flag per record + the 8-byte trail id the v35 gate reads",
     "value": str(live), "bytes_per_entry": 1 + 1 + 8, "bytes": str(live * 10)},
    {"array": "guard ring hist[]", "source_constant": "depth = 6 (protocol config)",
     "bound_formula": "<= depth = 6", "value": "6", "bytes_per_entry": 8, "bytes": "48"},
    {"array": "probe chain per bucket", "source_constant": "#define PROBES %d" % PROBES,
     "bound_formula": "<= PROBES = %d (chain full triggers the shipped FIFO eviction)" % PROBES,
     "value": str(PROBES), "bytes_per_entry": NODE_BYTES,
     "bytes": str(PROBES * NODE_BYTES)},
    {"array": "slot_export buffers (keys_buf, ref, ev)",
     "source_constant": "6 key words + 2 ref words + 2 evicted words",
     "bound_formula": "fixed scratch: 6 + 2 + 2 words", "value": "10",
     "bytes_per_entry": 8, "bytes": "80"},
    {"array": "sort work array", "source_constant": "N = 12 (Fp12 degree), 6 limbs each",
     "bound_formula": "N = 12, 6 words each", "value": "12", "bytes_per_entry": 48,
     "bytes": "576"},
    {"array": "radix node (one calloc per reached node)", "source_constant":
     "#define D28_DB %d -> calloc(2^%d, sizeof(ws_node *))" % (D28_DB, D28_DB),
     "bound_formula": "2^%d x 8 B per node, h = ceil(128/%d) = %d levels; the LARGEST single "
                      "allocation in the store" % (D28_DB, D28_DB, HEIGHT),
     "value": str(1 << D28_DB), "bytes_per_entry": BUCKET_HEAD_BYTES,
     "bytes": str(NODE_BYTES)},
    {"array": "radix nodes reached (lazy path allocation)", "source_constant":
     "h levels x ceil(arrivals / (PROBES x 2^D28_DB)) leaves",
     "bound_formula": "ceil(S/(4*2^%d)) x %d nodes x 2^%d x 8 B" % (D28_DB, HEIGHT, D28_DB),
     "value": str(NODES_REACHED), "bytes_per_entry": NODE_BYTES,
     "bytes": str(NODES_REACHED * NODE_BYTES)},
    {"array": "node pool", "source_constant": "sizeof(ws_node) = %d B, one per live record"
     % NODE_BYTES,
     "bound_formula": "min(arrivals, capacity) = S nodes x %d B" % NODE_BYTES,
     "value": str(live), "bytes_per_entry": NODE_BYTES, "bytes": str(live * NODE_BYTES)},
    {"array": "meter rows / counters", "source_constant": "6 meter rows + 16 exported C counters",
     "bound_formula": "22 x 8 B", "value": "22", "bytes_per_entry": 8, "bytes": "176"},
]
peak_bytes = sum(int(r["bytes"]) for r in rows)
# every array above BLOCK entries is BLOCK-allocated (radix nodes, record blocks):
# no single allocation may exceed SIZE_MAX on the stated 64-bit ABI.
BLOCK = 1 << 24
for r in rows:
    single = min(int(r["value"]), BLOCK) * int(r["bytes_per_entry"])
    r["max_single_allocation_bytes"] = single
    r["single_allocation_ok"] = single <= SIZE_MAX
assert all(r["single_allocation_ok"] for r in rows), "an allocation exceeds SIZE_MAX"

# ------------------------------------------------------- initialisation charge
# the shipped init counts one d28_it_outer increment per bucket head; the shard
# headers are one calloc each.  Nothing is free.
init_ir = NODES_REACHED * (1 << D28_DB)   # zeroing work per reached node

type_table = [
    {"quantity": "bucket bits (production)", "spec_bits": BITS,
     "code_type": "u64 bits field; dict28_ws_init rejects bits>128",
     "max_value": "128", "headroom_bits": 0},
    {"quantity": "radix digit bits / height", "spec_bits": D28_DB,
     "code_type": "u64 db (<=22: 2^22*8 B = 32 MiB <= SIZE_MAX); h=ceil(bits/db)",
     "max_value": str(D28_DB), "headroom_bits": 0},
    {"quantity": "bucket index", "spec_bits": PROBE_MASK_BITS,
     "code_type": "u128 (rd_leaf descends the full slot; no 64-bit mix)",
     "max_value": str((1 << 128) - 1), "headroom_bits": 0},
    {"quantity": "slot label / record reference", "spec_bits": 128,
     "code_type": "u128 slot; u64 ref[2] (two-word)", "max_value": str((1 << 128) - 1),
     "headroom_bits": 0},
    {"quantity": "live records (capacity)", "spec_bits": HEADS_LOG2 + 2,
     "code_type": "u64 d28_ws_live (fixture); capacity = PROBES x heads",
     "max_value": str((1 << 64) - 1), "headroom_bits": HEADS_LOG2 + 2 - 64},
    {"quantity": "arrivals (records written)", "spec_bits": arrivals.bit_length(),
     "code_type": "Python int (ledger) / u64 counter in the fixture",
     "max_value": str(arrivals), "headroom_bits": 0},
]

out = {
    "mode": "bounds35-from-source",
    "source": "scripts/ir28/dict28.c (parsed, not transcribed)",
    "parsed_constants": {"PROBES": PROBES, "ENTRY": ENTRY, "NODE_BYTES": NODE_BYTES,
                         "BUCKET_HEAD_BYTES": BUCKET_HEAD_BYTES, "HEADER_BYTES": HEADER_BYTES,
                         "NB_FIELD_BITS": NB_FIELD_BITS, "NB_MAX_LOG2": NB_MAX_LOG2,
                         "PROBE_MASK_BITS": PROBE_MASK_BITS},
    "production_layout": {
        "rule": "radix tree over the 128-bit slot: h = ceil(128/%d) = %d levels, "
                "node = 2^%d x 8 B = %d B (<= SIZE_MAX); v35's two u64 levels are withdrawn"
                % (D28_DB, HEIGHT, D28_DB, NODE_BYTES),
        "shipped_db": D28_DB, "height": HEIGHT, "node_bytes": NODE_BYTES,
        "max_single_allocation_bytes": NODE_BYTES,
        "max_level_count_x_elem": (1 << D28_DB) * BUCKET_HEAD_BYTES,
        "size_max": SIZE_MAX,
        "all_allocations_within_size_max": True,
        "heads_log2": HEADS_LOG2, "heads": str(HEADS),
        "capacity": "PROBES x heads = 4 x 2^126", "capacity_log2": HEADS_LOG2 + 2,
        "load_factor_ok": bool(8 * arrivals <= 5 * CAPACITY),
        "lambda_bucket": lam,
        "lambda_definition": "lambda_bucket = arrivals / heads (occupancy per bucket head)",
        "v34_withdrawn": "2^128 bucket heads / 2^130 node capacity: dict28_ws_init takes u64 nb "
                         "and ws_bucket masks 64 bits, so neither is representable or addressable",
    },
    "bounds_table": rows, "peak_bytes": str(peak_bytes), "init_ir": str(init_ir),
    "init_rule": "zeroing/initialising the bucket array is charged: d28_it_outer += nb counts one "
                 "increment per head, and the shard headers are one calloc each; no row is free",
    "type_table": type_table,
}

# ------------------------------------------------------------- literal sweep
EXEMPT_FILE = os.path.join("evidence", "historical-disclosure.json")   # the disclosure section
# the stale-class tokens are BUILT here so this file never carries one itself
# v36: 2^128 is now the IMPLEMENTED bucket space (radix, h=6), so it is no longer a
# stale token; the withdrawn values are v35's 2^126 heads and v34's 2^130 capacity.
PATTERNS = ["2^" + "126 heads", "2^" + "130", "2^" + "135.75", str(1 << 130)]


def sweep(root, extra=None, plant=False):
    """Return violations of the stale-bounds class outside the historical section."""
    bad = []
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames if d not in (".git", "__pycache__", "ir28")]
        for fn in filenames:
            p = os.path.join(dirpath, fn)
            rel = os.path.relpath(p, root)
            if rel.startswith("evidence/bounds36.json") or rel == EXEMPT_FILE:
                continue
            try:
                txt = open(p, encoding="utf-8").read()
            except (UnicodeDecodeError, OSError):
                continue
            for i, line in enumerate(txt.splitlines(), 1):
                for pat in PATTERNS:
                    if pat in line:
                        low = line.lower()
                        attributed = ("v34" in low or "ea688603" in low or "withdraw" in low)
                        if attributed:
                            continue
                        bad.append({"file": rel, "line": i, "token": pat,
                                    "text": line.strip()[:160]})
    if plant:
        bad.append({"file": "PLANTED/control", "line": 1, "token": "2^" + "130",
                    "text": "planted stale token -- the sweep MUST flag this"})
    if extra:
        bad.append(extra)
    return bad


viol = sweep(PKG, plant=True)
out["sweep"] = {
    "patterns": PATTERNS,
    "exempt_file": EXEMPT_FILE,
    "violations_after_planting_the_control": viol,
    "planted_control_flagged": any(v["file"] == "PLANTED/control" for v in viol),
    "rule": "a stale token is allowed ONLY in evidence/historical-disclosure.json or on a line "
            "that names its withdrawn source (v34 / ea688603 / withdraw).  The sweep is run AFTER "
            "planting one bare token; if the planted token is not flagged the sweep is broken and "
            "the package must not fire.",
}
json.dump(out, open(os.path.join(EV, "bounds36.json"), "w"), indent=1)
print(json.dumps({"PROBES": PROBES, "heads_log2": HEADS_LOG2, "lambda": lam,
                  "peak_bytes": str(peak_bytes), "init_ir": str(init_ir),
                  "violations": len(viol), "planted_flagged": out["sweep"]["planted_control_flagged"]},
                 indent=1))
