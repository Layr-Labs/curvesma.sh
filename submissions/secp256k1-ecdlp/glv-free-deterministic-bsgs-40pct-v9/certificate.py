#!/usr/bin/env python3
"""certificate.py - exact-integer machine verification of the v9 ledger claim.
Host integer arithmetic only (no floats, no secp group operations). Exit 0 iff
every check passes. Re-derived independently of certificate.json. Usage:
python3 certificate.py"""
import sys
r = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
C = 1000
m = 2**128
h = -(-r // m)
G = 78037248
H = m + (h - 1)
R = {
 "baby_write": m*64, "baby_reduce": m*32, "baby_scatter": 6*m*64,
 "giant_write": h*64, "giant_reduce": h*32, "giant_scatter": 6*h*64,
 "merge_scan": (m+h)*64, "heads": 6*((2**48)*4 + 2**40),
 "constants": 512*G, "exceptions": (2**20)*G, "cand_verify": 4*1024*G,
 "fixed_setup": 2**40, "output": 2**40,
}
Rtot = sum(R.values())
T = 1000*H + Rtot                 # total instruction equivalents; W = T/C
Tr = 10000*H + 12*Rtot            # 1.2x repricing of every R row; W_r = Tr/10000
q = 100000
cb = 12982492                          # claim*q for the submitted claim 3245623/25000
fr = 13670178                     # frontier 136.70178 scaled by q
checks = [
 ("h_eq_ceil_r_over_m_eq_2pow128_minus_1", h == 2**128 - 1),
 ("full_coverage_mh_ge_r", m*h >= r),
 ("giants_minimal_m_hm1_lt_r", m*(h-1) < r),
 ("at_most_two_matches_mh_lt_2r", m*h < 2*r),
 ("claim_proof_T_q_le_1000_q_2_cb", pow(T, q) <= (C**q) << cb),
 ("repriced_1p2x_margin_ge_0p10", pow(Tr, q) <= (10000**q) << (cb - 10000)),
 ("claim_below_frontier_by_ge_0p10", pow(T, q) << 100000 <= (C**q) << (fr - 100000)),
 ("peak_bytes_lt_2pow137", 2*(m+h)*64 + 6*(2**48)*8 + 2**24 < 2**137),
 ("claim_is_reduced", __import__("math").gcd(3245623, 25000) == 1),
 ("claim_times_q_eq_cb", 3245623 * q == cb * 25000),
]
bad = [n for n, ok in checks if not ok]
for n, ok in checks:
    print(("PASS" if ok else "FAIL"), n)
sys.exit(1 if bad else 0)
