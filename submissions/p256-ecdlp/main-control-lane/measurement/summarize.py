"""Per-step instruction counts from a run.sh output directory:
  python3 -I summarize.py <out-dir> --steps-per-n 65536 --C 1000 [--n1 4 --n2 8] [--mode single:1]
Per-step counts are (Ir(n2) - Ir(n1)) / ((n2 - n1) * steps_per_n), the organizer's reference-cost
method; the kernel's final stderr line of key=value event counts is reported for n2.
"""
import argparse
import json
import re


def summary(path):
    with open(path) as f:
        for line in f:
            if line.startswith("summary:"):
                return int(line.split()[1])
    raise SystemExit("no summary line in " + path)


def events(path):
    try:
        lines = [l for l in open(path).read().splitlines()
                 if re.fullmatch(r"\w+=\S+( \w+=\S+)*", l.strip())]
    except FileNotFoundError:
        return {}
    if not lines:
        return {}
    out = {}
    for tok in lines[-1].split():
        k, _, v = tok.partition("=")
        try:
            out[k] = int(v)
        except ValueError:
            out[k] = v
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("work")
    ap.add_argument("--steps-per-n", type=int, required=True)
    ap.add_argument("--C", type=int, required=True)
    ap.add_argument("--n1", type=int, default=64)
    ap.add_argument("--n2", type=int, default=128)
    ap.add_argument("--mode", action="append", default=[],
                    help="name:steps_per_n for MODES runs (sizes MODE_N1/MODE_N2 via --mode-n1/--mode-n2)")
    ap.add_argument("--mode-n1", type=int, default=4096)
    ap.add_argument("--mode-n2", type=int, default=8192)
    a = ap.parse_args()
    w, n1, n2 = a.work, a.n1, a.n2
    total = {n: summary("%s/step-%d.callgrind" % (w, n)) for n in (n1, n2)}
    extra = (n2 - n1) * a.steps_per_n
    per_total = (total[n2] - total[n1]) / extra
    out = {"Ir_total": total, "events_n2": events("%s/step-%d.stderr" % (w, n2)),
           "per_step_total_Ir": round(per_total, 3), "C": a.C}
    out["per_step_units_all_instructions"] = round(per_total / a.C, 6)
    for spec in a.mode:
        name, _, spn = spec.partition(":")
        m = {n: summary("%s/mode-%s-%d.callgrind" % (w, name, n)) for n in (a.mode_n1, a.mode_n2)}
        out["mode_" + name] = {"Ir": m, "per_unit_Ir": round((m[a.mode_n2] - m[a.mode_n1]) / ((a.mode_n2 - a.mode_n1) * int(spn or 1)), 3)}
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
