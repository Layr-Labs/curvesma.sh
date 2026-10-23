"""Per-step instruction split of the blst GT walk from a measure.sh output directory.

  python3 -I analyze.py <out-dir> --n1 256 --n2 512 --s S --label NAME

total = (Ir(step, n2) - Ir(step, n1)) / ((n2 - n1) * 64); mul, sqr, restart, store and prf are per call from the
mode runs (MODE_N1, MODE_N2); counts are the event-count differences of the two step runs. The restart comb count is
the comb_mults field of the restartbench stderr. C = 20000 instructions per charge unit.
"""
import argparse
import json
import os
import re

C = 20000
WALKS = 64
COUNT_KEYS = ["steps", "adds", "sqrs", "dps", "cyc2", "cycl", "restarts", "cap_restarts", "chained", "escapes",
              "escape_dps", "comb_mults"]


def summary(path):
    for line in open(path):
        if line.startswith("summary:"):
            return int(line.split()[1])
    raise SystemExit("no summary line in " + path)


def events(path):
    lines = [l for l in open(path).read().splitlines() if re.fullmatch(r"\w+=\S+( \w+=\S+)*", l.strip())]
    return {k: int(v) for k, v in (t.split("=") for t in lines[-1].split())}


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
    d, md = args.out_dir, args.mode_n2 - args.mode_n1

    def per_call(mode):
        p2, p1 = f"{d}/mode-{mode}-{args.mode_n2}.callgrind", f"{d}/mode-{mode}-{args.mode_n1}.callgrind"
        return (summary(p2) - summary(p1)) / md if os.path.exists(p2) else None

    total = (summary(f"{d}/step-{args.n2}.callgrind") - summary(f"{d}/step-{args.n1}.callgrind")) / (
        (args.n2 - args.n1) * WALKS)
    mul, sqr = per_call("mul"), per_call("sqr")
    e1, e2 = events(f"{d}/step-{args.n1}.stderr"), events(f"{d}/step-{args.n2}.stderr")
    ev = {k: e2[k] - e1[k] for k in COUNT_KEYS}
    steps = ev["steps"]
    a, q = ev["adds"] / steps, ev["sqrs"] / steps
    other = total - a * mul - q * sqr
    per = {"total": round(total, 3), "reference_blst_fp12_mul_per_call": round(mul, 3),
           "cyclotomic_sqr_per_call": round(sqr, 3), "other": round(other, 3)}
    restart, store, prf = per_call("restartbench"), per_call("storebench"), per_call("prfbench")
    if restart is not None:
        per["restart_per_call"] = round(restart, 3)
        per["restart_comb_mults_per_call"] = round(
            events(f"{d}/mode-restartbench-{args.mode_n2}.stderr")["comb_mults"] / args.mode_n2, 4)
    if store is not None:
        per["store_per_call"] = round(store, 3)
    if prf is not None:
        per["prf_per_call"] = round(prf, 3)
    res = {"label": args.label, "s": args.s, "S_square_slots": round(args.s * args.ra), "ra": args.ra, "C": C,
           "interval": {"n1": args.n1, "n2": args.n2, "walks": WALKS, "steps": steps},
           "per_main_step_instructions": per, "interval_counts": ev,
           "c2_other_over_C": round(other / C, 6), "kappa_sqr_over_C": round(sqr / C, 6),
           "ref_mul_over_C": round(mul / C, 6)}
    with open(f"{d}/analysis.json", "w") as f:
        json.dump(res, f, indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
