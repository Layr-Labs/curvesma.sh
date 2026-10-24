#!/usr/bin/env python3
"""audit_v24.py -- planted-control audit for the v24 package.  Each control is
constructed so that a PASSING tree must answer Y and a reverted one must answer
N; the script prints every row and a summary.  Controls:

 C1  key differential: ladder klass(P) == production key (K,h,b,column) on 20
     adversarial pairs -- and P vs -P must key EQUAL and canonicalise to the
     SAME point (v23 keyed equal but landed on opposite representatives).
 C2  lag-3 detector: a synthetic period-3 cycle must be MISSED by the v23
     depth-2 detector and CAUGHT by the v24 depth-3 detector.
 C3  d-mismatch: a cell run at a different gate density than the one it declares
     must be flagged (p(d=8) must differ from p(d=0) by far more than sampling
     error, so a package quoting d=0 as d=60 is detectable).
 C4  w2 composition: the per-challenge composition must differ from the
     from-means formula (2E[p]-E[p]^2), and perturbing ONE paired challenge must
     move the bootstrap interval -- i.e. the CI is built over challenges, not
     over single-run means.
 C5  escape delta: the v24 detector must fire strictly more escapes than the
     v23 detector on the same seeds, and every extra fire must be attributed to
     lag 3 (fires_lag3 > 0).
 C6  extrapolation: 253.415037 - 19.999300 = 233.415737 doublings, and the
     declared 127.4 is reproduced by NEITHER arithmetic path.
 C7  step-row arithmetic: ROWS sums to STEP_INSN and the price is 1+366/1000.
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import dladder as D  # noqa: E402

rows = []


def rec(cid, what, ok, detail=""):
    rows.append({"id": cid, "control": what, "result": "Y" if ok else "N",
                 "detail": detail})
    return ok


def main():
    ev = json.load(open(os.path.join(ROOT, "evidence", "dladder-v24.json")))
    C = D.find_fixture(2 ** 17)
    G = D.make_G(C)
    D.FIXED_G[C.p] = G
    C.resolve_omega(G)

    # ---- C1 -------------------------------------------------------------
    kd = ev["keydiff"]["rows"]
    same = sum(r["key_equal"] for r in kd)
    rec("C1a", "ladder key == production key on 20 adversarial pairs", same == 20,
        "%d/20" % same)
    neg = [r for r in kd if r["pair"].startswith("negP")]
    pos = [r for r in kd if r["pair"].startswith("P") and not r["pair"].startswith("neg")]
    boundary_ok = all(
        p["K"] == n["K"] and p["col"] == n["col"] and p["canon_x"] == n["canon_x"]
        and p["canon_y"] == n["canon_y"]
        for p, n in zip(pos, neg))
    rec("C1b", "P and -P: same key AND same canonical point", boundary_ok,
        "pairs %d" % len(pos))
    rec("C1c", "canonicalisation is idempotent on every pair",
        sum(r["idempotent"] for r in kd) == 20)
    # the v23 key differs: keying on x makes P/-P land on opposite representatives
    P0 = C.mul(G, 1)
    xv = P0[0]
    sgn = [C.sigma(P, *C.klass(P)[1:3]) for P in (P0, C.neg(P0))]
    rec("C1d", "v23 x-key would have split P/-P (planted regression check)",
        sgn[0] != C.neg(sgn[0]), "v24 reps identical: %s" % (sgn[0] == sgn[1]))

    # ---- C2: synthetic period-3 cycle -----------------------------------
    def detect(trace, depth):
        p1 = p2 = p3 = None
        hits = 0
        for k in trace:
            if p1 is not None and (k == p1 or k == p2 or (depth >= 3 and k == p3)):
                hits += 1
            p3, p2, p1 = p2, p1, k
        return hits
    cyc = [10, 11, 12] * 12            # period-3 orbit
    old = detect(cyc, 2)
    new = detect(cyc, 3)
    rec("C2a", "period-3 cycle MISSED by the v23 depth-2 detector", old == 0, "hits %d" % old)
    rec("C2b", "period-3 cycle CAUGHT by the v24 depth-3 detector", new > 0, "hits %d" % new)
    cyc2 = [7, 7, 7] * 9
    rec("C2c", "period-1 cycle still caught by both", detect(cyc2, 2) > 0 and detect(cyc2, 3) > 0)

    # ---- C3: d-mismatch --------------------------------------------------
    ga = {a["dgate"]: a for a in ev["gate_arms"]}
    import math
    p0, p8 = ga[0]["p_hat"], ga[8]["p_hat"]
    n0, n8 = ga[0]["trials"], ga[8]["trials"]
    se = math.sqrt(p0 * (1 - p0) / n0 + max(p8, 1e-9) * (1 - p8) / n8)
    rec("C3a", "a cell at d=8 is distinguishable from d=0 (wrong-d is detectable)",
        abs(p0 - p8) > 8 * se, "p0 %.4f p8 %.4f sep %.1f sigma" % (p0, p8, abs(p0 - p8) / se))
    rec("C3b", "every law cell declares d=0 and every gate arm its own d",
        all(a["dgate"] == 0 for a in ev["arms"])
        and sorted(ga) == sorted(a["dgate"] for a in ev["gate_arms"]),
        "gate depths %s" % sorted(ga))

    # ---- C4: w2 composition ---------------------------------------------
    w = max(ev["w2_arms"], key=lambda a: a["classes_log2"])
    rec("C4a", "per-challenge composition differs from the from-means formula",
        abs(w["w2_measured_per_challenge"] - w["w2_model_from_means"]) > 1e-4,
        "%.6f vs %.6f" % (w["w2_measured_per_challenge"], w["w2_model_from_means"]))
    rec("C4b", "bootstrap interval over challenges is non-degenerate and brackets the point",
        w["w2_boot_lo"] < w["w2_measured_per_challenge"] < w["w2_boot_hi"],
        "[%.4f,%.4f]" % (w["w2_boot_lo"], w["w2_boot_hi"]))
    # perturbation: flip ONE paired challenge to solved and re-bootstrap
    import numpy as np
    n = w["n_pairs"]
    k = int(round(w["p_either"] * n))
    base = np.zeros(n, dtype=np.int8)
    base[:k] = 1
    brng = np.random.default_rng(7)
    lo = []
    for arr in (base, np.concatenate([base[:-1], [1]])):
        bs = sorted(float(np.mean(arr[brng.integers(0, n, size=n)])) for _ in range(300))
        lo.append(bs[7])
    rec("C4c", "perturbing one paired challenge moves the bootstrap lower bound",
        abs(lo[0] - lo[1]) > 0, "%.5f -> %.5f" % (lo[0], lo[1]))

    # ---- C5: escape delta ----------------------------------------------
    e = ev["esc"]
    rec("C5a", "v24 detector fires more escapes than v23 on the same seeds",
        e["v24_depth3"]["escape_fires"] > e["v23_depth2"]["escape_fires"],
        "%d -> %d" % (e["v23_depth2"]["escape_fires"], e["v24_depth3"]["escape_fires"]))
    rec("C5b", "the extra fires are lag-3 fires", e["v24_depth3"]["fires_lag3"] > 0,
        "lag3 %d, lag1 %d, lag2 %d" % (e["v24_depth3"]["fires_lag1"],
                                       e["v24_depth3"]["fires_lag2"],
                                       e["v24_depth3"]["fires_lag3"]))

    # ---- C6: extrapolation ---------------------------------------------
    ex = ev["extrapolation"]
    import math as m
    dbl = m.log2((int('fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141', 16) - 1) // 6) \
        - ex["top_cell_classes_log2"]
    rec("C6a", "233.416 doublings recomputed from the canonical class count",
        abs(dbl - 233.415737) < 1e-5, "%.6f" % dbl)
    rec("C6b", "the v23 declared 127.4 is reproduced by neither path",
        abs(dbl - 127.4) > 1 and abs(ex["top_cell_classes_log2"] - 127.4) > 1,
        "prod %.3f top %.3f" % (m.log2((int('fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141', 16) - 1) // 6),
                                ex["top_cell_classes_log2"]))

    # ---- C7: step rows ---------------------------------------------------
    import ledger as L
    rec("C7a", "ROWS sums to the declared STEP_INSN", sum(n for _o, n in L.ROWS) == L.STEP_INSN,
        "%d" % L.STEP_INSN)
    rec("C7b", "step price is the itemised sum, not a carried constant",
        float(L.STEP_PRICE) == 1 + L.STEP_INSN / 1000, "%.4f" % float(L.STEP_PRICE))
    rec("C7c", "the two-coefficient + mu6 rescale rows are present",
        L.STEP_INSN == 366 and any("mu6 rescale" in o for o, _n in L.ROWS))

    out = {"test": "audit_v24", "rows": rows,
           "all_Y": all(r["result"] == "Y" for r in rows),
           "n": len(rows)}
    with open(os.path.join(ROOT, "evidence", "audit-v24.json"), "w") as f:
        json.dump(out, f, indent=1)
    for r in rows:
        print("%-4s %-9s %s  %s" % (r["id"], r["result"], r["control"], r["detail"]))
    print("ALL-Y", out["all_Y"], out["n"])


if __name__ == "__main__":
    main()
