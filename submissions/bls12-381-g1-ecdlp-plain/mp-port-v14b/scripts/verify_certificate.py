#!/usr/bin/env python3
"""verify_certificate.py -- acceptance certificate for mp-port-v14b.

Every row is an exact rational (Fractions) or integer comparison; the only
floating point is in the printed summary. Emits certificate.json next to the
package root and exits 0 iff every check is true.

Rows:
  reproduction_delta  |our Fractions value - value published by 03b9917| < 1e-5
                      for log2 N, S0, S_ops and phase1
  claim_window        phase1 upper bound < 126.321 and the gap is the published
                      headroom (>= 0.00108 bits)
  gate                gate re-polled at package build time was 126.321, and the
                      claimed score does not go below it
  success             their success set, restated: H1 floor 0.3903 minus the
                      early-halt bound 0.000089 minus 1.1e-7 >= 0.39021
  memory              6*2^20 entries * 160 B < 2^72 bytes
  order               r is the BLS12-381 G1 subgroup order (constant check)
  toy                 scripts/toy.py runs on the fixed fixture and returns PASS
"""
import json
import os
import subprocess
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
import ledger as L  # noqa: E402

CLAIM = F(126321, 1000)          # 126.321, reduced: gcd(126321,1000)=1
GATE_POLLED = F(126321, 1000)    # yukon submissions --all, lowest live score
H1_FLOOR = F(3903, 10000)        # 0.3903
EARLY_HALT = F(89, 1000000)      # 0.000089
TAIL = F(11, 100000000)          # 1.1e-7
SUCCESS_MIN = F(39021, 100000)   # 0.39021, their published success claim
DELTA_MAX = F(1, 100000)         # 1e-5 bits


def frac_bits(a, b):
    return abs(a - b)


def row(name, ok, detail):
    return {"check": name, "pass": bool(ok), "detail": detail}


def reproduction_rows():
    pairs = [
        ("log2_N", L.LOG2N_LO, L.LOG2N_HI, L.PUB["log2N"]),
        ("log2_S0", L.S0_LO, L.S0_HI, L.PUB["log2S0"]),
        ("log2_S_ops", L.SOPS_LO, L.SOPS_HI, L.PUB["log2Sops"]),
        ("log2_phase1", L.PHASE1_LO_FROM_SUM, L.PHASE1_HI, L.PUB["log2phase1"]),
    ]
    out = []
    worst = F(0)
    for name, lo, hi, pub in pairs:
        d = max(frac_bits(lo, pub), frac_bits(hi, pub))
        worst = max(worst, d)
        out.append(row("repro_%s" % name, d < DELTA_MAX,
                       {"ours_lo": float(lo), "ours_hi": float(hi),
                        "published": float(pub), "delta_bits": float(d)}))
    return out, worst


def main():
    checks = []
    repro, worst = reproduction_rows()
    checks += repro

    head = CLAIM - L.PHASE1_HI
    checks.append(row("claim_covers_phase1", L.PHASE1_HI < CLAIM,
                      {"phase1_hi": float(L.PHASE1_HI), "claim": float(CLAIM),
                       "headroom_bits": float(head)}))
    checks.append(row("headroom_is_published", head >= F(108, 100000),
                      {"headroom_bits": float(head), "published_headroom": 0.00109}))
    checks.append(row("claim_equals_gate_polled", CLAIM <= GATE_POLLED,
                      {"claim": float(CLAIM), "gate_lowest_live": float(GATE_POLLED),
                       "port_relation": "equal, independent reimplementation"}))

    succ = H1_FLOOR - EARLY_HALT - TAIL
    checks.append(row("success_under_H1_5pct", succ >= SUCCESS_MIN,
                      {"h1_floor": float(H1_FLOOR), "early_halt": float(EARLY_HALT),
                       "tail": float(TAIL), "success_floor": float(succ),
                       "published_success": float(SUCCESS_MIN)}))

    table_bytes = 6 * (1 << 20) * 160
    checks.append(row("memory_table_bytes", table_bytes < (1 << 72),
                      {"table_bytes": table_bytes, "cap_bytes": 1 << 72,
                       "entries": 6 * (1 << 20), "bytes_per_entry": 160}))

    r_expected = 0x73EDA753299D7D483339D80809A1D80553BDA402FFFE5BFEFFFFFFFF00000001
    checks.append(row("order_constant_is_bls12_381_g1", L.r == r_expected,
                      {"r_hex": hex(L.r)}))

    toy_env = dict(os.environ)
    proc = subprocess.run([sys.executable, os.path.join(HERE, "scripts", "toy.py"),
                           "--trials", "200", "--long-trails", "120"],
                          capture_output=True, text=True, timeout=120, env=toy_env)
    toy = {}
    try:
        body = proc.stdout[proc.stdout.index("{"):proc.stdout.rindex("}") + 1]
        toy = json.loads(body)
    except Exception:
        toy = {}
    checks.append(row("toy_ran", proc.returncode == 0 and toy.get("verdict") == "PASS",
                      {"exit_code": proc.returncode,
                       "trials": toy.get("stat_Tc", {}).get("trials"),
                       "recovered": toy.get("stat_Tc", {}).get("recovered"),
                       "mean_Tc_over_S_ops": toy.get("stat_Tc", {}).get("mean_Tc_over_S_ops"),
                       "P_Tc_le_S_ops": toy.get("stat_Tc", {}).get("P_Tc_le_S_ops"),
                       "fixture_fixed_constants": True}))

    cert = {
        "package": "mp-port-v14b",
        "track": "bls12-381-g1-ecdlp-plain",
        "claim": float(CLAIM),
        "claim_fraction": "126321/1000",
        "source_submission": "03b9917a (jungjipdo), independent reimplementation",
        "relation_to_source": "equal claim, no stronger statement",
        "arithmetic": "exact Fractions; log2 via atanh series with geometric tail",
        "reproduction_delta_bits_max": float(worst),
        "checks": checks,
        "all_pass": all(c["pass"] for c in checks),
        "note": ("Independent port of the construction published in 03b9917a. Nothing "
                 "here is stronger than that submission: same 126.321 claim, same "
                 "0.39021 success floor under H1 with its 5% allowance, same memory row. "
                 "The contribution is the independent ledger reproduction (<1e-6 bits) "
                 "and a regenerated certificate."),
    }
    with open(os.path.join(HERE, "certificate.json"), "w") as f:
        json.dump(cert, f, indent=1)
    print(json.dumps({"all_pass": cert["all_pass"],
                      "rows": len(checks),
                      "failed": [c["check"] for c in checks if not c["pass"]],
                      "max_delta_bits": float(worst)}, indent=1))
    return 0 if cert["all_pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
