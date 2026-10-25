"""v20 flatness / drift analysis over the stability ladder.

The question the review left open (d32b4f73, flaw 1) is whether the measured
success rate can be carried from 16 selector periods, a 49-bit group and one
walker to 2^87.67 periods, a 256-bit group and 2^20 walkers.  This script tests
the two axes that can be moved independently in the fixture:

  * COVERAGE axis (period count at fixed R, fixed W, fixed lambda),
  * WALKER axis (W at fixed R, fixed coverage, fixed lambda),
  * R axis (R at fixed W, fixed coverage, fixed lambda).

Verdicts are mechanical: two cells' 95% Clopper-Pearson intervals are computed
exactly (bisection on the binomial tail, no normal approximation), and an axis is
called DRIFT only when the interval of the extreme cells is DISJOINT -- the
conservative test, which reports FLAT whenever the data cannot tell the
difference.  A weighted slope in log2(axis value) is reported alongside.

Two planted controls make sure the detector is not a rubber stamp: counts that
drift must be flagged, counts that are flat must not be.
"""
import json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "evidence")


def cp_interval(k, n, alpha=0.025):
    if k == 0:
        lo = 0.0
    else:
        def tail(p):
            lp, lq = math.log(p), math.log1p(-p)
            s = 0.0
            lg = math.lgamma
            for i in range(k, n + 1):
                s += math.exp(lg(n + 1) - lg(i + 1) - lg(n - i + 1) + i * lp + (n - i) * lq)
            return s
        a, b = 0.0, 1.0
        for _ in range(60):
            m = 0.5 * (a + b)
            if tail(m) < alpha:
                a = m
            else:
                b = m
        lo = a
    if k == n:
        hi = 1.0
    else:
        def tail2(p):
            lp, lq = math.log(p), math.log1p(-p)
            s = 0.0
            lg = math.lgamma
            for i in range(0, k + 1):
                s += math.exp(lg(n + 1) - lg(i + 1) - lg(n - i + 1) + i * lp + (n - i) * lq)
            return s
        a, b = 0.0, 1.0
        for _ in range(60):
            m = 0.5 * (a + b)
            if tail2(m) > alpha:
                a = m
            else:
                b = m
        hi = a
    return lo, hi


def points(k, n):
    lo, hi = cp_interval(k, n)
    return {"wins": k, "trials": n, "p_hat": k / n, "cp95_lower": lo, "cp95_upper": hi}


def slope(cells):
    xs = [math.log2(c["axis"]) for c in cells]
    ys = [c["p_hat"] for c in cells]
    ws = [c["trials"] for c in cells]
    sw = sum(ws)
    mx = sum(w * x for w, x in zip(ws, xs)) / sw
    my = sum(w * y for w, y in zip(ws, ys)) / sw
    num = sum(w * (x - mx) * (y - my) for w, x, y in zip(ws, xs, ys))
    den = sum(w * (x - mx) ** 2 for w, x in zip(ws, xs))
    return num / den if den else 0.0


def axis(name, rows, field="cross"):
    """rows: list of (axis_value, cell_label, wins, trials) already in the chosen mode."""
    rows = [r for r in rows if r is not None]
    cells = []
    for v, label, k, n in rows:
        p = points(k, n)
        p.update({"axis": v, "label": label})
        cells.append(p)
    first, last = cells[0], cells[-1]
    up = last["cp95_lower"] > first["cp95_upper"]
    down = last["cp95_upper"] < first["cp95_lower"]
    verdict = "DRIFT_UP" if up else ("DRIFT_DOWN" if down else "FLAT")
    return {
        "axis": name, "mode": field, "cells": cells,
        "slope_per_log2_axis": slope(cells),
        "extreme_intervals_disjoint": bool(up or down),
        "verdict": verdict,
        "test": "exact CP-95 intervals of the extreme cells; FLAT unless disjoint",
    }


def main():
    def load(label):
        f = os.path.join(OUT, "ladder20_%s.json" % label)
        if not os.path.exists(f):
            return None
        return json.load(open(f))

    def row(label, axisval, mode):
        d = load(label)
        if d is None:
            return None
        if mode == "cross":
            return (axisval, label, d["wins_cross"], d["trials"])
        return (axisval, label, d["wins"], d["trials"])

    axes = []
    # R axis at fixed W=64, 16 periods (cross-only: the production collision mix)
    axes.append(axis("R@W64_p16", [row(l, r, "cross") for r, l in
                                   ((256, "X_R256_W64_p16"), (1024, "X_R1024_W64_p16"),
                                    (4096, "X_R4096_W64_p16"))]))
    # coverage axis at fixed R=256, W=64
    axes.append(axis("periods@R256_W64", [row(l, p, "cross") for p, l in
                                          ((16, "X_R256_W64_p16"), (64, "X_R256_W64_p64"),
                                           (256, "X_R256_W64_p256"))]))
    # coverage axis at fixed R=1024, W=16
    axes.append(axis("periods@R1024_W16", [row(l, p, "cross") for p, l in
                                           ((16, "P_R1024_p16"), (64, "P_R1024_p64"),
                                            (256, "P_R1024_p256"), (1024, "P_R1024_p1024"))]))
    # coverage axis at fixed R=4096, W=4
    axes.append(axis("periods@R4096_W4", [row(l, p, "cross") for p, l in
                                          ((16, "P_R4096_p16"), (64, "P_R4096_p64"),
                                           (256, "P_R4096_p256"), (1024, "P_R4096_p1024"))]))
    # coverage axis at fixed R=65536, W=1 (W=1 has no cross term, so pooled)
    axes.append(axis("periods@R65536_W1", [row(l, p, "pooled") for p, l in
                                           ((16, "P_R65536_p16"), (64, "P_R65536_p64"),
                                            (256, "P_R65536_p256"), (1024, "P_R65536_p1024"))]))
    # walker axis at fixed R=4096, 256 periods, pooled (W=1 has no cross term)
    axes.append(axis("walkers@R4096_p256", [row(l, w, "pooled") for w, l in
                                            ((1, "W_R4096_W1"), (16, "W_R4096_W16"),
                                             (64, "W_R4096_W64"))]))

    # planted controls for the detector itself
    def synth(name, counts, n):
        rows = [(i + 1, "synthetic_%d" % i, k, n) for i, k in enumerate(counts)]
        r = axis(name, rows, "pooled")
        return r

    ctl_drift = synth("control_drifting", [96, 84, 62, 40], 200)
    ctl_flat = synth("control_flat", [95, 92, 99, 93], 200)
    ctl_mild = synth("control_mild_decline", [95, 92, 89, 86], 200)

    out = {
        "why": "verdict d32b4f73 flaw 1: does success per trial stay put as coverage, "
               "walker count and R move, at fixed lambda = 0.634 and the production selector",
        "axes": axes,
        "detector_controls": {
            "drifting_counts_must_flag": ctl_drift,
            "flat_counts_must_not_flag": ctl_flat,
            "mild_decline_sensitivity": ctl_mild,
            "detected_drift": ctl_drift["verdict"] != "FLAT",
            "did_not_flag_flat": ctl_flat["verdict"] == "FLAT",
            "sensitivity_note": "the test is deliberately blunt: a monotone decline smaller "
                                "than the CP width at these trial counts is reported FLAT, "
                                "which is why the claim uses the carrier cell's own CP lower "
                                "bound and not a fitted trend",
        },
        "all_axes_flat": all(a["verdict"] == "FLAT" for a in axes),
    }
    json.dump(out, open(os.path.join(OUT, "drift_v20.json"), "w"), indent=1, sort_keys=True)
    for a in axes:
        print(a["axis"], a["verdict"], "slope %.4f" % a["slope_per_log2_axis"],
              " ".join("%s %.3f[%.3f,%.3f]" % (c["label"].split("_", 1)[-1], c["p_hat"],
                                               c["cp95_lower"], c["cp95_upper"])
                       for c in a["cells"]))
    print("controls: drift=%s flat=%s mild=%s" % (ctl_drift["verdict"], ctl_flat["verdict"],
                                                 ctl_mild["verdict"]))


if __name__ == "__main__":
    main()
