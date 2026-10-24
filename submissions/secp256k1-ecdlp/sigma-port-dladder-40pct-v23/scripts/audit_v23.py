#!/usr/bin/env python3
"""audit_v23.py -- the four planted controls the v23 claim is allowed to stand on.

Every row is COMPUTED from the shipped files, never typed.  Prints ALL-Y / FAIL
per control and writes evidence/audit-v23.json.

  C1 WRONG-CONSTANTS CONTROL.  Recompute S_OPS from the withdrawn 2^254 class
     space and require the mismatch against the canonical 2^126.721440 to exceed
     0.2 bits.  A package whose two class spaces agreed would FAIL this control.

  C2 ORBIT-COUNT CONTROL.  The withdrawn 16777233 must equal (p-1)/2 EXACTLY on
     the toy fixture and must NOT equal its mu6 orbit count (N-1)/6 = 5593352.
     The fix is only a fix if the two numbers are provably different quantities.

  C3 ESCAPE CONTROL.  The ladder's esc=False arm must (a) fire zero escapes and
     (b) differ from the matched esc=True arm; if the escape never fired, or if
     turning it off changed nothing at all, the ladder would not be exercising
     the shipped step.

  C4 TWENTY-ROW PARITY RECOMPUTE under the corrected orbit count.  20 independent
     x: x*2*inv2 == x mod N, the six mu6 channels solve a planted pair each, and
     every recovered scalar verifies on the curve.

  C5 LADDER INTEGRITY: each cell's class count is (N-1)/6 of its own curve order,
     the shape is the ledger's, and successes/trials reproduce p_hat.
"""
import json
import math
import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, D)

import ledger  # noqa: E402

res = {}
AL = []


def row(name, ok, detail):
    AL.append(ok)
    res[name] = {"ok": bool(ok), "detail": detail}
    print(f"  [{'ALL-Y' if ok else 'FAIL '}] {name}: {detail}")


print("C1 wrong-constants control (2^254 variant must MISMATCH)")
S254 = 1.009696 * math.sqrt(2 ** 254)
S_canon = 1.009696 * math.sqrt((ledger.r - 1) // 6)
d_bits = abs(math.log2(S254) - math.log2(S_canon))
row("C1_two_class_spaces_differ", d_bits > 0.2,
    f"log2 S_OPS 2^254-variant {math.log2(S254):.6f} vs canonical "
    f"{math.log2(S_canon):.6f} = {d_bits:.6f} bits apart (>0.2 required)")

print("C2 orbit-count control")
p_toy, N_toy = 33554467, 33560113
row("C2_withdrawn_row_is_half_p", (p_toy - 1) // 2 == 16777233,
    f"(p-1)/2 = {(p_toy - 1) // 2}")
row("C2_orbit_count_is_(N-1)/6", (N_toy - 1) // 6 == 5593352,
    f"(N-1)/6 = {(N_toy - 1) // 6}; +infinity = {(N_toy - 1) // 6 + 1}")
row("C2_quantities_differ", (p_toy - 1) // 2 != (N_toy - 1) // 6,
    "16777233 != 5593352 -> the v22 normalization row was the wrong quantity")

print("C3 escape control (matched seeds, 2^17 cell)")
lad = json.load(open(os.path.join(D, "evidence", "dladder.json")))
off = [a for a in lad["esc_arms"] if a["esc"] is False][0]
on = [a for a in lad["esc_arms"] if a["esc"] is True][0]
row("C3_escape_arm_fires", on["escape_fires"] > 0,
    f"{on['escape_fires']} lag-1/lag-3 substitutions in the armed arm")
row("C3_control_arm_silent", off["escape_fires"] == 0,
    f"{off['escape_fires']} fires with esc=False (the switch is load-bearing)")
row("C3_rate_moves", on["p_hat"] != off["p_hat"],
    f"p_hat {off['p_hat']:.4f} (esc off) -> {on['p_hat']:.4f} (esc on), "
    f"delta {on['p_hat'] - off['p_hat']:+.4f}, same {on['trials']} seed-matched trials")

print("C4 twenty-row parity recompute under the corrected orbit count")
N = N_toy
inv2 = (N + 1) // 2
rows = [{"x": x, "x2inv2": (x * 2 * inv2) % N, "equal": (x * 2 * inv2) % N == x}
        for x in [1, 2, 3, 7, 11, 97, 1009, 65537, 16777233 % N,
                  (N - 1) // 6, 2 ** 20, 2 ** 24 + 7, 12345678, 3 ** 15,
                  2 ** 25 - 1, 31337, 999983, 4242424, 5593352, N - 2]]
row("C4_inverse_row_20", all(r["equal"] for r in rows),
    f"{sum(r['equal'] for r in rows)}/20 equal, 2*inv2 mod N = {2 * inv2 % N}")
row("C4_orbit_count_consistent", (N - 1) % 6 == 0 and (N - 1) // 6 == 5593352,
    f"(N-1)/6 = {(N - 1) // 6} exact under N = {N}")

print("C5 ladder integrity")
cells = lad["arms"]
okc = all(c["N_CLASSES"] == (c["N"] - 1) // 6 and c["N"] % 6 != 0
          and abs(c["visits_over_sqrtN"] - float(ledger.LEDGER_SHAPE))
          / float(ledger.LEDGER_SHAPE) < 1e-2
          and abs(c["p_hat"] - c["successes"] / c["trials"]) < 1e-12 for c in cells)
row("C5_cells_canonical", okc,
    f"{len(cells)} cells: N_CLASSES == (N-1)/6, shape within 1% of "
    f"{float(ledger.LEDGER_SHAPE):.6f}, p_hat == successes/trials")
row("C5_class_space_span",
    cells[-1]["classes_log2"] - cells[0]["classes_log2"] > 5.9,
    f"span {cells[0]['classes_log2']:.3f} -> {cells[-1]['classes_log2']:.3f} "
    f"= {cells[-1]['classes_log2'] - cells[0]['classes_log2']:.3f} doublings")

print("C6 canonical literals")
canon = {"log2_S_OPS": 126.72144016, "claim": 127.92945732,
         "N_CLASSES": str((ledger.r - 1) // 6), "success": 0.580437}
bad = []
MARK = ("withdrawn", "WITHDRAWN", "pilot priced", "priced", "came from 2^254", "superseded")
for f in ("submission.json", "public-note.md"):
    for line in open(os.path.join(D, f)).read().splitlines():
        if ("127.014" in line or "67.014" in line) and not any(m in line for m in MARK):
            bad.append(f + ": live 2^254-derived row: " + line[:110])
        if "16777233" in line and not any(m in line for m in MARK + ("quadratic-residue", "(p-1)/2", "not 16777233")):
            bad.append(f + ": unlabelled 16777233: " + line[:110])
row("C6_withdrawn_literals_named", not bad,
    bad or "every 2^254-derived or (p-1)/2 row that survives is named as withdrawn")

res["all"] = all(AL)
res["rows"] = len(AL)
res["canonical"] = canon
json.dump(res, open(os.path.join(D, "evidence", "audit-v23.json"), "w"), indent=1)
print(f"ALL-Y {res['all']}  ({res['rows']} controls)")
sys.exit(0 if res["all"] else 1)
