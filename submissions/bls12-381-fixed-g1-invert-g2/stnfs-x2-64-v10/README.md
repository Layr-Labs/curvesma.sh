# STNFS with 100k pool, both-bars alpha (5.15, 0.85), theta 2901/8000, slacks (295, 226), UDUP 0.6: 121.15796

Track bls12-381-fixed-g1-invert-g2: return Q in G2 with e(P,Q)=z for nonidentity P in G1.

Both-bars alpha (5.15, 0.85) (f 50k 5.085035 + 0.065 margin, g 50k
0.731543 + 0.118 margin with q/(q+1) p = 2 fix), 100k pool, theta 2901/8000
(4518 rows, margin 0.1201), slacks (295, 226) (collection 2^120.9506),
inputs 63.949, 746.9/499.8 (two-log 2^115.8579), UDUP 0.6, meter 2^116.59.
Ledger <= 2^121.15796, success >= 0.395327. Analytic, not a full-size run.

Entry points: audit_refinement.py, unit_entries.py, evidence16.py.
