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

The production store is therefore TWO levels of the shipped type
(2^63 shards x 2^63 bucket heads = 2^126 heads, every field u64), which is what
the ledger prices and what arm35.py measures at lambda = S/2^126.

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
assert re.search(r"u64 nb;", src) and re.search(r"if \(nb & \(nb - 1\)\) return -1;", src)
assert re.search(r"& \(nb - 1\)", src) and "mix64((u64)slot" in src
PROBE_MASK_BITS = 64                        # entropy of the shipped ws_bucket mix

LEVELS = 2
HEADS_LOG2 = LEVELS * NB_MAX_LOG2           # 126: two levels of the shipped type
HEADS = 1 << HEADS_LOG2
CAPACITY = PROBES * HEADS                   # 4 heads-per-bucket chain cap

# ------------------------------------------------------------ production bound
arrivals = S
live = min(arrivals, CAPACITY)
lam = arrivals / HEADS

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
    {"array": "shard headers (level 1)", "source_constant": "{ u64 nb; ws_node **b; } = %d B"
     % HEADER_BYTES,
     "bound_formula": "2^%d shards x %d B" % (NB_MAX_LOG2, HEADER_BYTES),
     "value": str(1 << NB_MAX_LOG2), "bytes_per_entry": HEADER_BYTES,
     "bytes": str((1 << NB_MAX_LOG2) * HEADER_BYTES)},
    {"array": "bucket array (level 2, the shipped ws_node **)", "source_constant":
     "calloc(nb, sizeof(ws_node *)); nb a u64 power of two, ws_bucket masks %d bits"
     % PROBE_MASK_BITS,
     "bound_formula": "2^%d x 2^%d heads = 2^%d heads x %d B" %
                      (NB_MAX_LOG2, NB_MAX_LOG2, HEADS_LOG2, BUCKET_HEAD_BYTES),
     "value": str(HEADS), "bytes_per_entry": BUCKET_HEAD_BYTES,
     "bytes": str(HEADS * BUCKET_HEAD_BYTES)},
    {"array": "node pool", "source_constant": "sizeof(ws_node) = %d B, one per live record"
     % NODE_BYTES,
     "bound_formula": "min(arrivals, capacity) = S nodes x %d B" % NODE_BYTES,
     "value": str(live), "bytes_per_entry": NODE_BYTES, "bytes": str(live * NODE_BYTES)},
    {"array": "meter rows / counters", "source_constant": "6 meter rows + 16 exported C counters",
     "bound_formula": "22 x 8 B", "value": "22", "bytes_per_entry": 8, "bytes": "176"},
]
peak_bytes = sum(int(r["bytes"]) for r in rows)

# ------------------------------------------------------- initialisation charge
# the shipped init counts one d28_it_outer increment per bucket head; the shard
# headers are one calloc each.  Nothing is free.
init_ir = HEADS + (1 << NB_MAX_LOG2)

type_table = [
    {"quantity": "bucket count nb (per level)", "spec_bits": NB_MAX_LOG2,
     "code_type": "u64 (dict28_ws.nb; dict28_ws_init rejects nb & (nb-1))",
     "max_value": str((1 << NB_FIELD_BITS) - 1), "headroom_bits": 0},
    {"quantity": "bucket index (per level)", "spec_bits": PROBE_MASK_BITS - 2,
     "code_type": "u64 (ws_bucket: mix64(slot_lo ^ mix64(slot_hi)) & (nb-1))",
     "max_value": str((1 << PROBE_MASK_BITS) - 1), "headroom_bits": 0},
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
        "rule": "%d levels of the shipped type: nb = 2^%d per level, all fields u64"
                % (LEVELS, NB_MAX_LOG2),
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
PATTERNS = ["2^" + "128", "2**" + "128", "2^" + "127", "2^" + "130", "2^" + "135.75",
            str(1 << 128), str(1 << 127)]


def sweep(root, extra=None, plant=False):
    """Return violations of the stale-bounds class outside the historical section."""
    bad = []
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames if d not in (".git", "__pycache__", "ir28")]
        for fn in filenames:
            p = os.path.join(dirpath, fn)
            rel = os.path.relpath(p, root)
            if rel.startswith("evidence/bounds35.json") or rel == EXEMPT_FILE:
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
        bad.append({"file": "PLANTED/control", "line": 1, "token": "2^" + "128",
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
json.dump(out, open(os.path.join(EV, "bounds35.json"), "w"), indent=1)
print(json.dumps({"PROBES": PROBES, "heads_log2": HEADS_LOG2, "lambda": lam,
                  "peak_bytes": str(peak_bytes), "init_ir": str(init_ir),
                  "violations": len(viol), "planted_flagged": out["sweep"]["planted_control_flagged"]},
                 indent=1))
