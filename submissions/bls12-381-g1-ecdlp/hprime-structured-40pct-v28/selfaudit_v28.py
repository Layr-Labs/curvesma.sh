#!/usr/bin/env python3
"""v28 self-audit: re-derive the ledger from certificate.json, assert every printed
integer identity, re-check kappa, re-check the arithmetic audit and the ladder rows.
Exit 1 on any failure.  Run from the package directory."""
import json, math, os, sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
FAILS = []


def chk(cond, msg):
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


def main():
    cert = json.load(open(HERE + "/certificate.json"))
    led = json.load(open(HERE + "/ledger-v28.json"))
    der = json.load(open(HERE + "/derive-v28.json"))
    kap = json.load(open(HERE + "/kappa-v28.json"))
    lad = json.load(open(HERE + "/structured-ladder.json"))
    aud = json.load(open(HERE + "/arithmetic-audit.json"))
    sub = json.load(open(HERE + "/submission.json"))

    N, nb = cert["N"], 1 << 66
    ell = -(-N // nb)
    chk(N % 2 == 1, "N is odd (the v27 equality was the defect)")
    chk(nb * (ell - 1) < N <= nb * ell, "nb*(ell-1) < N <= nb*ell")
    chk(der["construction"]["ell"] == ell, "derive ell matches")
    chk(der["construction"]["records_per_side"] == nb * ell, "records/side == nb*ell (integer identity)")
    chk(led["rows"]["records"] == 2 * nb * ell * 2238, "records row == 2*nb*ell*2238")
    chk(sum(led["rows"].values()) == led["total_instructions"], "integer sum == total")
    chk(led["S_written"] == sum(led["rows"][k] for k in led["rows"]
                               if k not in ("records", "base_entry_nw50",
                                            "walker_state_nw50", "gap_point_tables")),
        "printed S equals the nine secondary rows")
    chk(led["S_written"] != led["S_v27_argument_written"], "stale v27 S replaced")
    W = math.log2(float(F(led["total_instructions"], 2000)))
    chk(abs(W - led["W_log2"]) < 1e-12, "W_log2 recomputes exactly")
    claim, gate = F(116813, 1000), F(1168271, 10000)
    chk(W <= float(claim), "W <= claim 116.8130")
    chk(float(claim) <= float(gate), "claim <= gate 116.8271")
    chk(led["claim_minus_W"] > 0, "claim covers W (headroom > 0)")
    chk(led["rows"]["gap_point_tables"] > 0 and led["rows"]["base_entry_nw50"] > 0,
        "gap-point and pool-entry rows are positive charges")
    chk(kap["kappa_charge"] >= kap["kappa_measured_max"], "kappa charged >= measured max")
    chk(kap["price"]["E_bad_v28_log2"] < kap["price"]["cap_bad_log2"] - 10,
        "E_bad at least 10 bits under the candidate cap")
    chk(aud["failures"] == 0, "arithmetic audit: %d patterns, 0 failures" % aud["patterns_checked"])
    chk(all(r["K"] == lad["K"] for r in lad["ladder"]), "all ladder rungs at the same density K")
    chk(lad["primary"]["floor_met"], "primary rung CP-95 lower >= 0.39")
    chk(0.5093 < lad["K"] < 0.5095, "ladder density matches production 0.5094")
    txt = sub["claims"][0]["statement"]
    chk(len(txt) <= 7965, "claims[0] length %d <= 7965" % len(txt))
    chk(sub["claims"][0]["score"] == {"numerator": "116813", "denominator": "1000"},
        "submission score == 116813/1000")
    total = sum(os.path.getsize(os.path.join(HERE, f)) for f in os.listdir(HERE)
                if os.path.isfile(os.path.join(HERE, f)))
    chk(total <= 262144, "package bytes %d <= 262144" % total)
    print("\n%d failures" % len(FAILS))
    return 1 if FAILS else 0


if __name__ == "__main__":
    raise SystemExit(main())
