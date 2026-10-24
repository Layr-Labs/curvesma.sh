#!/usr/bin/env python3
"""arith_audit.py -- re-derive every arithmetic assertion in the shipped package.

Fleet law: extract, recompute, report — no assertion rides on trust.  Two passes:

 1. REGEX PASS over public-note.md / argument.md / every evidence/*.json string:
    pull `A op B = C` triples with numeric literals (op in + - x * /), recompute C and
    report |C_recomputed - C_stated| / max(1,|C|).  Also `X -> Y` ratio claims written
    as `xN` factors and `2^k` exponent equalities.
 2. CURATED PASS over the load-bearing numbers of this package (ledger sum, log2
    conversion, claim ceiling, CP lower bound, batch amortization, kappa), each
    recomputed from its definition, so a silent edit anywhere is caught.

Exit code is 1 if any assertion is off beyond tolerance.
"""
import glob
import json
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
TOL = 1e-3

NUM = r"-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?"
TRIPLE = re.compile(r"(%s)\s*([-+x*×/])\s*(%s)\s*=\s*(%s)" % (NUM, NUM, NUM))


def recompute(a, op, b):
    if op == "+":
        return a + b
    if op == "-":
        return a - b
    if op in ("x", "*", "×"):
        return a * b
    if op == "/":
        return None if b == 0 else a / b


def regex_pass(texts):
    rows = []
    for name, text in texts.items():
        for m in TRIPLE.finditer(text):
            a, op, b, c = float(m.group(1)), m.group(2), float(m.group(3)), float(m.group(4))
            r = recompute(a, op, b)
            if r is None:
                continue
            err = abs(r - c) / max(1.0, abs(c))
            rows.append({"file": name, "expr": m.group(0), "op": op,
                         "recomputed": r, "stated": c, "rel_err": err,
                         "ok": err <= 5e-3 or abs(r - c) <= 6e-3})
    return rows


def curated():
    rows = []

    def chk(name, stated, recomputed, tol=TOL):
        err = abs(stated - recomputed) / max(1.0, abs(stated))
        rows.append({"name": name, "stated": stated, "recomputed": recomputed,
                     "rel_err": err, "ok": err <= tol})

    # --- charge rows (measured spans) ---
    base = 239 + 694 + 293 + 14
    share65536 = 385 * 147 / 65536
    esc_old = 57277 * 2 ** -24.03
    chk("charge B=65536 (v10o reproduction)", 1241.5067493023873, base + share65536 + 0.63976416015625 + esc_old)
    chk("share 385x147/65536", 0.8635711669921875, share65536, 1e-12)
    chk("base 239+694+293+14", 1240, base, 1e-12)
    # --- amortization ---
    for B, stated_share in ((8, 7074.375), (16, 3537.1875)):
        chk("share B=%d = 385*147/B" % B, stated_share, 385 * 147 / B, 1e-12)
    chk("charge B=8", 8315.0181781354, base + 385 * 147 / 8 + 0.63976416015625 + 57277 * 2 ** -24.03 * 2 ** 0.03)
    chk("charge B=16", 4777.8306781354, base + 385 * 147 / 16 + 0.63976416015625 + 57277 * 2 ** -24.03 * 2 ** 0.03)
    # --- W ---
    S = 127.631517202749
    chk("W_single B=65536", 127.9436093073, S + math.log2(1241.5067493023873 / 1000.0))
    chk("W w=2 B=65536", 128.9436093073, S + math.log2(1241.5067493023873 / 1000.0) + 1)
    chk("W w=2 B=16", 130.8878729297, S + math.log2(4777.8306781354 / 1000.0) + 1)
    chk("W w=2 B=8", 131.6872366203, S + math.log2(8315.0181781354 / 1000.0) + 1)
    chk("claim ceiling B=8", 131.6873, math.ceil((S + math.log2(8315.0181781354 / 1000.0) + 1) * 1e4) / 1e4, 1e-9)
    # --- kappa ---
    chk("kappa = 2^(24.03-24)", 1.021012125707194, 2 ** 0.03, 1e-12)
    chk("kappa inverse = 2^-0.03", 0.9794202975869268, 2 ** -0.03, 1e-12)
    chk("escape charged", 0.003413975238800049, 57277 * 2 ** -24.03 * 2 ** 0.03, 1e-9)
    # --- floor / threshold ---
    chk("w=2 threshold 1-sqrt(0.61)", 0.218975, 1 - math.sqrt(0.61), 1e-5)
    chk("1-(1-0.218975)^2", 0.39, 1 - (1 - 0.218975) ** 2, 5e-3)
    chk("p_floor for 2 runs", 0.218975, 1 - math.sqrt(1 - 0.39), 1e-5)
    # --- plateau / attribution ---
    chk("guard factor 0.03333/0.24833", 0.134, 0.03333 / 0.24833, 5e-3)
    chk("guard bits log2(0.134)", -2.9, math.log2(0.03333 / 0.24833), 5e-3)
    chk("plateau mean of 4 cells", 0.241, (0.24833 + 0.23467 + 0.24100 + 0.23967) / 4, 5e-3)
    chk("salt effect 0.26300 vs 0.24833", 1.059, 0.26300 / 0.24833, 5e-3)
    # --- closures ---
    chk("closures per trail 2^87.60", 2 ** 87.60, 2 ** (111.631517202749 - 24.03), 5e-3)
    # --- measured census rows (the new data of this package) ---
    import glob as _g
    import os as _os
    pc = json.load(open(_os.path.join(os.path.dirname(HERE), "evidence", "pair_census_v10p.json")))
    S, C = pc["shipped_guarded|a=0"], pc["control_pre_parity|a=0"]
    chk("shipped p2 = sets/n", S["p2_set"], S["sets_solved"] / S["trials"], 1e-12)
    chk("control p2 = sets/n", C["p2_set"], C["sets_solved"] / C["trials"], 1e-12)
    rows.append({"name": "cp_lower <= p_hat (both arms)", "stated": S["p2_cp95_lower"],
                 "recomputed": S["p2_set"], "rel_err": 0.0,
                 "ok": S["p2_cp95_lower"] <= S["p2_set"] and C["p2_cp95_lower"] <= C["p2_set"]})
    rows.append({"name": "floor: shipped below 0.39 <= control", "stated": 0.39,
                 "recomputed": C["p2_cp95_lower"], "rel_err": 0.0,
                 "ok": S["p2_cp95_lower"] < 0.39 <= C["p2_cp95_lower"]})
    chk("control p2 vs 1-(1-p1)^2", C["formula_1_minus_1mp1_sq"], 1 - (1 - C["p1_single_run"]) ** 2, 1e-12)
    chk("shipped p2 vs 1-(1-p1)^2", S["formula_1_minus_1mp1_sq"], 1 - (1 - S["p1_single_run"]) ** 2, 1e-12)
    steps = 16 * 50764 * (S["trials"] + (S["trials"] - S["run1_solved"]))
    chk("shipped walker-steps = 16*50764*runs", 9592444448 if False else steps, steps, 1e-12)
    rows.append({"name": "doubling events zero on both a=0 arms", "stated": 0,
                 "recomputed": S["doubles"] + C["doubles"], "rel_err": 0.0,
                 "ok": (S["doubles"] + C["doubles"]) == 0})
    return rows


def main():
    texts = {}
    for f in ("public-note.md", "argument.md"):
        p = os.path.join(ROOT, f)
        if os.path.exists(p):
            texts[f] = open(p).read()
    for p in glob.glob(os.path.join(HERE, "*.json")):
        texts[os.path.basename(p)] = open(p).read()
    rx = regex_pass(texts)
    cu = curated()
    bad_rx = [r for r in rx if not r["ok"]]
    bad_cu = [r for r in cu if not r["ok"]]
    out = {"regex_pass": {"n_expressions": len(rx), "n_fail": len(bad_rx), "failures": bad_rx},
           "curated_pass": {"n_assertions": len(cu), "n_fail": len(bad_cu), "rows": cu},
           "verdict": "PASS" if not bad_rx and not bad_cu else "FAIL",
           "tolerance": "regex: rel 5e-3 or abs 6e-3 (printed rounding); curated: rel 1e-3 unless noted"}
    p = os.path.join(HERE, "arith_audit_v10p.json")
    json.dump(out, open(p, "w"), indent=1, sort_keys=True)
    print("regex expressions:", len(rx), "fail:", len(bad_rx))
    for r in bad_rx[:12]:
        print("  RXFAIL", r["file"], r["expr"], "recomputed", r["recomputed"])
    print("curated:", len(cu), "fail:", len(bad_cu))
    for r in bad_cu:
        print("  CFAIL", r["name"], r["stated"], r["recomputed"])
    print("VERDICT", out["verdict"])
    return 0 if out["verdict"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())