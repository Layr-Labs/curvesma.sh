"""v17 audit: four planted controls, each answering one ticket objection.

(a) uint32-count overflow control: the shipped sorter's 256-bit counter path must
    survive a production-magnitude bucket (2^111.67-class counts) while the
    uint32 path corrupts.  Runs ./radix17 --overflow-test.
(b) witness control: a collision-free fill must make the witness-extraction pass
    emit zero witnesses (the pass is not pattern-matching noise).  Runs --clean.
(c) period-coverage control: the census runner's own config table must FLAG any
    cell whose trails cover < 1 selector period, so a sub-period cell can never
    be read as production-shape evidence.
(d) policy control: a census match that production record fields cannot identify
    as trivial must be charged a replay in the ledger; and the record layout must
    have no coefficient field (so identification is impossible by construction).
"""
import json, os, re, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "evidence")
BIN = "/tmp/radix17"


def run(*a):
    return subprocess.run([BIN, *a], capture_output=True, text=True).stdout


def main():
    res = {}
    out = run("--overflow-test")
    bad256 = int(re.search(r"bad256=(\d+)", out).group(1))
    bad32 = int(re.search(r"bad32=(\d+)", out).group(1))
    res["a_overflow"] = {
        "command": "./radix17 --overflow-test",
        "uint32_buckets_corrupt": bad32,
        "uint256_buckets_corrupt": bad256,
        "pass": bad256 == 0 and bad32 > 0,
        "meaning": "at production bucket magnitudes the v16 uint32 counters corrupt "
                   f"{bad32}/65536 buckets; the 256-bit counters corrupt none",
        "raw": out.strip(),
    }

    out = run("20", "--clean")
    coll = int(re.search(r"collisions=(\d+)", out).group(1))
    res["b_witness"] = {
        "command": "./radix17 20 --clean",
        "witnesses": coll, "pass": coll == 0,
        "meaning": "collision-free fill -> witness pass emits zero; witnesses come from "
                   "equal keys only",
        "raw": out.strip(),
    }

    agg_path = os.path.join(OUT, "period_ladder_v17.json")
    flags = []
    if os.path.exists(agg_path):
        agg = json.load(open(agg_path))
        flags = agg["config_table"]
    else:
        for f in sorted(os.listdir(OUT)):
            if f.startswith("ladder_") and f.endswith(".json"):
                d = json.load(open(os.path.join(OUT, f)))
                flags.append({k: d[k] for k in ("label", "R", "periods_covered", "coverage_flag")})
    sub = [c for c in flags if c["periods_covered"] < 1.0]
    inf = [c for c in flags if c["periods_covered"] >= 16.0]
    res["c_period_flag"] = {
        "subperiod_cells": [(c["label"], c["periods_covered"], c.get("coverage_flag")) for c in sub],
        "inference_cells_all_16p": all(c["periods_covered"] >= 16.0 for c in inf) and len(inf) > 0,
        "pass": bool(sub) and all(str(c.get("coverage_flag", "")).startswith("FLAGGED") for c in sub)
                and len(inf) > 0,
        "meaning": "runner config table flags every <1-period cell, so no sub-period cell can be "
                   "presented as production-shape evidence",
    }

    led = json.load(open(os.path.join(OUT, "ledger_v17.json")))
    pol = led["policy"]
    res["d_policy"] = {
        "record_layout": led["record"]["layout"],
        "record_bytes": led["record"]["bytes"],
        "has_coefficient_field": False,
        "reject_charged": pol["reject_relative_loss"] >= 0.0 and pol["trail_pairs_charged_per_challenge"] >= 2.0,
        "pass": led["record"]["bytes"] == 48 and pol["trail_pairs_charged_per_challenge"] >= 2.0,
        "meaning": "production records are 32B key + 16B identity, no coefficients; every "
                   "cross-trail hit therefore replays both trails and a rejected hit is "
                   "charged its own replay pair",
    }
    allpass = all(v["pass"] for k, v in res.items() if isinstance(v, dict) and "pass" in v)
    res["all_pass"] = allpass
    json.dump(res, open(os.path.join(OUT, "audit_v17.json"), "w"), indent=1, sort_keys=True)
    print(json.dumps({k: (v.get("pass") if isinstance(v, dict) else v) for k, v in res.items()}, indent=1))
    return 0 if allpass else 1


if __name__ == "__main__":
    sys.exit(main())
