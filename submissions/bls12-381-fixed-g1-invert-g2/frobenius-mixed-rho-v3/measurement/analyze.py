"""Per-step instruction split for the blst GT walk, from a measure.sh output directory.

  python3 -I analyze.py <out-dir> --n1 256 --n2 512 --s 0.447265625 --label <name>

Callgrind Ir, differenced as in the organizer's method, over the same interval for the walk
and for the event counts:
  total = (Ir(step, n2) - Ir(step, n1)) / ((n2 - n1) * WALKS)        instructions per main step
  mul   = (Ir(mode mul, MODE_N2) - Ir(mode mul, MODE_N1)) / (MODE_N2 - MODE_N1)
  sqr   = (Ir(mode sqr, MODE_N2) - Ir(mode sqr, MODE_N1)) / (MODE_N2 - MODE_N1)
  counts = kernel event counts of the step-n2 run minus those of the step-n1 run (same inputs)
  a = adds / steps, q = sqrs / steps   group operations per main step, including the extra
      operations of cycle checks and escapes
  other = total - a * mul - q * sqr    (key, table reads, exponents, DP ring, cycle filter)
  c2 = other / C,  kappa = sqr / C,  ref_mul = mul / C,  C = 20000.
The production cost per main step is a * 1 + q * kappa + c2; accounting.py takes a and q from
the uninstrumented runs (native_run.sh), which have far more steps.
"""
import argparse
import json
import os
import re

C = 20000
WALKS = 64
COUNT_KEYS = ["steps", "adds", "sqrs", "dps", "cyc2", "cycl", "restarts", "cap_restarts"]


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
    ap.add_argument("--n1", type=int, required=True)
    ap.add_argument("--n2", type=int, required=True)
    ap.add_argument("--mode-n1", type=int, default=4096)
    ap.add_argument("--mode-n2", type=int, default=8192)
    ap.add_argument("--s", type=float, required=True)
    ap.add_argument("--ra", type=int, default=1024)
    ap.add_argument("--label", required=True)
    args = ap.parse_args()
    d = args.out_dir
    steps_per_size = WALKS
    total = (summary(f"{d}/step-{args.n2}.callgrind") - summary(f"{d}/step-{args.n1}.callgrind")) / (
        (args.n2 - args.n1) * steps_per_size)
    mdiff = args.mode_n2 - args.mode_n1
    mul = (summary(f"{d}/mode-mul-{args.mode_n2}.callgrind") - summary(f"{d}/mode-mul-{args.mode_n1}.callgrind")) / mdiff
    sqr = (summary(f"{d}/mode-sqr-{args.mode_n2}.callgrind") - summary(f"{d}/mode-sqr-{args.mode_n1}.callgrind")) / mdiff
    e1 = events(f"{d}/step-{args.n1}.stderr")
    e2 = events(f"{d}/step-{args.n2}.stderr")
    ev = {k: e2[k] - e1[k] for k in COUNT_KEYS}
    steps = ev["steps"]
    a = ev["adds"] / steps
    q = ev["sqrs"] / steps
    other = total - a * mul - q * sqr
    # restart alone (mode restartbench, optional): one restart_walk call per iteration
    restart_path = f"{d}/mode-restartbench-{args.mode_n2}.callgrind"
    restart = None
    if os.path.exists(restart_path):
        restart = (summary(restart_path) - summary(f"{d}/mode-restartbench-{args.mode_n1}.callgrind")) / mdiff
    res = {
        "label": args.label,
        "s": args.s,
        "S_square_slots": round(args.s * args.ra),
        "ra": args.ra,
        "C": C,
        "interval": {"n1": args.n1, "n2": args.n2, "walks": WALKS, "steps": steps},
        "per_main_step_instructions": {
            "total": round(total, 3),
            "reference_blst_fp12_mul_per_call": round(mul, 3),
            "cyclotomic_sqr_per_call": round(sqr, 3),
            "adds_times_mul": round(a * mul, 3),
            "sqrs_times_sqr": round(q * sqr, 3),
            "other": round(other, 3),
            **({"restart_per_call": round(restart, 3)} if restart is not None else {}),
        },
        "interval_counts": ev,
        "per_step_fractions": {"adds_a": round(a, 9), "sqrs_q": round(q, 9),
                               "cycle_events_per_step": round((ev["cyc2"] + ev["cycl"]) / steps, 9)},
        "c2_other_over_C": round(other / C, 6),
        "kappa_sqr_over_C": round(sqr / C, 6),
        "ref_mul_over_C": round(mul / C, 6),
        "raw": {
            "step_callgrind": [f"step-{args.n1}.callgrind", f"step-{args.n2}.callgrind"],
            "step_stderr": [f"step-{args.n1}.stderr", f"step-{args.n2}.stderr"],
            "mul_callgrind": [f"mode-mul-{args.mode_n1}.callgrind", f"mode-mul-{args.mode_n2}.callgrind"],
            "sqr_callgrind": [f"mode-sqr-{args.mode_n1}.callgrind", f"mode-sqr-{args.mode_n2}.callgrind"],
            **({"restart_callgrind": [f"mode-restartbench-{args.mode_n1}.callgrind",
                                      f"mode-restartbench-{args.mode_n2}.callgrind"]} if restart is not None else {}),
            "dump": f"dump-{args.n2}.txt",
            "native_stderr": f"native-{args.n2}.stderr",
        },
    }
    with open(f"{d}/analysis.json", "w") as f:
        json.dump(res, f, indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
