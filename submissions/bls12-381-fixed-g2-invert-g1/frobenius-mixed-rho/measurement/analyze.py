"""Per-step instruction split for the blst GT walk, from a measure.sh output directory.

  python3 -I analyze.py <out-dir> --s <square probability> --label <name>

Costs (Callgrind Ir, differenced as in the organizer's method):
  total  = (Ir(step, 128 rounds) - Ir(step, 64 rounds)) / (64 rounds * 64 walks)
  mul    = (Ir(mode mul, 8192) - Ir(mode mul, 4096)) / 4096   one blst_fp12_mul (the reference op)
  sqr    = (Ir(mode sqr, 8192) - Ir(mode sqr, 4096)) / 4096   one blst_fp12_cyclotomic_sqr
Event fractions come from the kernel's final stderr line of the dump run (native-128.stderr).
  other  = total - adds_fraction * mul - squares_fraction * sqr
           (key, table lookup, exponent updates, DP test, cycle filter, state copies)
  c2 = other / C,  kappa = sqr / C,  ref_units_per_add = mul / C,  C = 20000.
The REF_FUNCS pass of measure.sh is not used: toggling on ref_mul also counts the callees of
blst_fp12_mul, so its totals equal the full run and cannot isolate the reference operation.
"""
import argparse
import json
import os
import re

C = 20000
WALKS = 64


def summary(path):
    for line in open(path):
        if line.startswith("summary:"):
            return int(line.split()[1])
    raise SystemExit("no summary line in " + path)


def events(path):
    lines = [l for l in open(path).read().splitlines() if re.fullmatch(r"\w+=\S+( \w+=\S+)*", l.strip())]
    out = {}
    for tok in lines[-1].split():
        k, _, v = tok.partition("=")
        out[k] = int(v)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out_dir")
    ap.add_argument("--s", type=float, required=True)
    ap.add_argument("--ra", type=int, default=1024)
    ap.add_argument("--label", required=True)
    args = ap.parse_args()
    d = args.out_dir
    per_step = 64 * WALKS
    total = (summary(f"{d}/step-128.callgrind") - summary(f"{d}/step-64.callgrind")) / per_step
    mul = (summary(f"{d}/mode-mul-8192.callgrind") - summary(f"{d}/mode-mul-4096.callgrind")) / 4096
    sqr = (summary(f"{d}/mode-sqr-8192.callgrind") - summary(f"{d}/mode-sqr-4096.callgrind")) / 4096
    ev = events(f"{d}/native-128.stderr")
    steps = ev["steps"]
    add_frac = ev["adds"] / steps
    sqr_frac = ev["sqrs"] / steps
    ref = add_frac * mul
    sq = sqr_frac * sqr
    other = total - ref - sq
    res = {
        "label": args.label,
        "s_requested": args.s,
        "S_square_slots": round(args.s * args.ra),
        "ra": args.ra,
        "C": C,
        "per_step_instructions": {
            "total_full_step": round(total, 3),
            "reference_blst_fp12_mul_per_call": round(mul, 3),
            "cyclotomic_sqr_per_call": round(sqr, 3),
            "ref_part_add_fraction_times_mul": round(ref, 3),
            "sqr_part_sqr_fraction_times_sqr": round(sq, 3),
            "other_key_table_exponents_dp_cycles": round(other, 3),
        },
        "fractions": {
            "adds": round(add_frac, 6),
            "squares": round(sqr_frac, 6),
            "distinguished_points": round(ev["dps"] / steps, 6),
            "cycle2_events": ev["cyc2"],
            "cycle_long_events": ev["cycl"],
            "steps": steps,
        },
        "c2_other_over_C": round(other / C, 6),
        "kappa_sqr_over_C": round(sqr / C, 6),
        "ref_mul_over_C": round(mul / C, 6),
        "raw": {
            "step_callgrind": ["step-64.callgrind", "step-128.callgrind"],
            "mul_callgrind": ["mode-mul-4096.callgrind", "mode-mul-8192.callgrind"],
            "sqr_callgrind": ["mode-sqr-4096.callgrind", "mode-sqr-8192.callgrind"],
            "dump": "dump-128.txt",
            "native_stderr": "native-128.stderr",
        },
    }
    with open(f"{d}/analysis.json", "w") as f:
        json.dump(res, f, indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
