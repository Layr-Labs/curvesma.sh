# STNFS with 150k pool, both-bars alpha (5.15, 0.85), theta 2881/8000, slacks (295, 226), UDUP 0.65: 121.083537139968

Track bls12-381-fixed-g1-invert-g2: return Q in G2 with e(P,Q)=z for nonidentity P in G1.

Both-bars alpha (5.15, 0.85) (f 50k 5.085035 + 0.065 margin, g 50k
0.731543 + 0.118 margin with q/(q+1) p = 2 fix), 150k pool, theta 2881/8000
(6628 rows, margin 0.1104), slacks (295, 226) (collection 2^120.9254),
inputs 63.834 (= 63.940 - 2 SE), 746.1/499.9 (two-log 2^113.242), UDUP 0.65, 12-relation
boot, meter 2^113.975. Ledger <= 2^121.083537139968, success >= 0.3952. Analytic, not a full-size run.

Entry points: audit_refinement.py, unit_entries.py, evidence16.py.
