#!/usr/bin/env python3
"""run_ladder_v24.py -- emits evidence/dladder-v24.json (all v24 ladder arms).

Arms:
  law   : 2^14 / 2^17 / 2^20 class space, n=1200 per cell, production key,
          lag-3 detector, d=0 (all visits recorded).
  w2    : n=1200 SAME-CHALLENGE PAIRS at 2^20 and 2^17 (two runs, fresh coins),
          per-challenge composition + bootstrap over challenges.
  gate  : 2^17 at d = 0,4,8,12,16 -- the DP gate law in RECORDS.
  esc   : 2^17, v23 detector (depth 2) vs v24 (depth 3), matched seeds.
  keydiff: 20 adversarial pairs, ladder key vs production key.
"""
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dladder as D

T = int(sys.argv[1]) if len(sys.argv) > 1 else 1200
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "evidence", "dladder-v24.json")

# ---- C1: the key differential (fixture 2^17, the middle cell) -----------------
C = D.find_fixture(2 ** 17)
G = D.make_G(C)
D.FIXED_G[C.p] = G
C.resolve_omega(G)
KEYDIFF = D.key_diff(C, G)

# ---- LAW ---------------------------------------------------------------------
law = [D.cell(2 ** 14, T, 11, L=16, R=32),
       D.cell(2 ** 17, T, 12, L=16, R=32),
       D.cell(2 ** 20, T, 13, L=16, R=32)]

# ---- W2 (same-challenge pairs) ----------------------------------------------
w2 = [D.cell(2 ** 20, T, 13, L=16, R=32, pairs=True),
      D.cell(2 ** 17, T, 12, L=16, R=32, pairs=True)]

# ---- GATE --------------------------------------------------------------------
gate = D.gate_scan(2 ** 17, 12, max(300, T // 2))

# ---- ESCAPE DELTA ------------------------------------------------------------
o, n = D.escape_delta(2 ** 17, 12, max(300, T // 2), pairs=True)

# ---- derived rows -----------------------------------------------------------
xs = [a["classes_log2"] for a in law]
ps = [a["p_hat"] for a in law]
mx = sum(xs) / len(xs)
my = sum(ps) / len(ps)
sxx = sum((x - mx) ** 2 for x in xs)
sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ps))
slope = sxy / sxx if sxx else 0.0
resid = [y - (my + slope * (x - mx)) for x, y in zip(xs, ps)]
s2 = sum(r * r for r in resid) / max(1, len(xs) - 2)
se = math.sqrt(s2 / sxx) if sxx else 0.0
# pooled CP-95 over all single-run trials
n_tot = sum(a["trials"] for a in law)
k_tot = sum(a["successes"] for a in law)
pooled = k_tot / n_tot
pooled_cp = D.cp_lower(k_tot, n_tot)

# extrapolation arithmetic (verdict point 2: 253.415 - 19.999 = 233.416)
PROD_LOG2 = math.log2((int('fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141', 16) - 1) // 6)
TOP = law[-1]["classes_log2"]
DOUBLINGS = PROD_LOG2 - TOP
DOUBLINGS_WRONG = PROD_LOG2 - 19.999

doc = {
    "test": "transfer_ladder_w2_v24",
    "mechanism": {
        "arithmetic": "affine elliptic-curve addition over F_p on real curves, batched "
                      "inversion (Montgomery's trick, one field inversion per step for "
                      "the whole lane fleet)",
        "table": "orbit rows J0(m)=[alpha_m]G+[d_m]Q, six mu6 slots per column, "
                 "m = (K*R)>>(bitlen(p)-1) with K = min(y,p-y) (PRODUCTION KEY)",
        "escape": "lag-1 and lag-3 class compare (prev1/prev2/prev3), alternate row "
                  "(m+1) mod R, K_ESC_MAX=8 then lane restart from fresh coins",
        "gate": "SplitMix fingerprint of the class key K, threshold 2^-d",
        "record": {"success": "two records in one class, mu6 solve over all six channels, "
                              "verified [x]G==Q",
                   "channels_solved": "1, omega, omega^2, -1, -omega, -omega^2"},
        "key_fix": "klass(P) = (K=min(y,p-y), h=argmin_h beta^h x, b=[y>p-y]); the "
                   "canonical representative is sigma(P,h,b) = (beta^h x, K).  P and -P "
                   "map to the SAME point (control C1).  v23 keyed on x only, so P and -P "
                   "shared a key but landed on opposite representatives.",
        "only_shrunk": "the class space: (N-1)/6 = 16415 / 131253 / 1048057 against "
                       "production (r-1)/6",
    },
    "keydiff": {"rows": KEYDIFF, "equal": sum(r["key_equal"] for r in KEYDIFF),
                "n": len(KEYDIFF),
                "idempotent": sum(r["idempotent"] for r in KEYDIFF),
                "note": "key_equal compares (K, h, b, column m) against the production "
                        "key function on the same point"},
    "arms": law,
    "slope_per_doubling": slope,
    "slope_se": se,
    "slope_ci95": [slope - 1.96 * se, slope + 1.96 * se],
    "span_doublings_measured": law[-1]["classes_log2"] - law[0]["classes_log2"],
    "pooled": {"n": n_tot, "k": k_tot, "p_hat": pooled, "cp95": pooled_cp},
    "w2_arms": w2,
    "gate_arms": gate,
    "esc": {"v23_depth2": o, "v24_depth3": n,
            "delta_fires": n["escape_fires"] - o["escape_fires"],
            "delta_lag3_fires": n["fires_lag3"],
            "delta_p_hat": n["p_hat"] - o["p_hat"]},
    "extrapolation": {
        "log2_prod_classes": PROD_LOG2,
        "top_cell_classes_log2": TOP,
        "doublings": DOUBLINGS,
        "doublings_as_declared_v23": 127.4,
        "arith_line": "log2((r-1)/6) = %.6f ; top cell = %.4f ; %.6f - %.4f = %.3f doublings"
                      % (PROD_LOG2, TOP, PROD_LOG2, TOP, DOUBLINGS),
        "wrong_233_check": "%.6f - 19.999 = %.3f (the verdict's recomputation against the "
                           "declared 127.4: the declared figure was the PRODUCTION-side "
                           "log2 minus a mis-typed cell exponent)" % (PROD_LOG2, DOUBLINGS_WRONG),
    },
    "protocol": {
        "trials_per_cell": T,
        "pairs_per_cell_w2": T,
        "same_challenge_pairing": "both runs of a pair use the SAME Q; fresh walker coins",
        "bootstrap": "400 resamples over the paired challenge rows",
        "key": "byte-identical to the production key function (selector.cell + twist_index)",
    },
}
with open(OUT, "w") as f:
    json.dump(doc, f, indent=1)
print(json.dumps({k: doc[k] for k in ("slope_per_doubling", "slope_se", "pooled",
                                      "extrapolation")}, indent=1)[:1200])
for a in law:
    print(a["p_hat"], a["successes"], a["trials"], a["classes_log2"],
          a["escape_fires"], a["fires_lag3"], a["elapsed_s"])
for a in w2:
    print("w2", a["classes_log2"], a["p_either"], a["w2_boot_lo"], a["w2_boot_hi"],
          a["w2_model_from_means"], a["p_run_pooled"])
for a in gate:
    print("gate", a["dgate"], a["p_hat"], a["visits"], a["successes"])
print("esc", o["p_hat"], n["p_hat"], o["escape_fires"], n["escape_fires"], n["fires_lag3"])
print("keydiff", doc["keydiff"]["equal"], "/", doc["keydiff"]["n"])
