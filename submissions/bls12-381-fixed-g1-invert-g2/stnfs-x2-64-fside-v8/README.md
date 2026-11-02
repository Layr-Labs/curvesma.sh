# STNFS with both-bars alpha (5.15, 0.85), theta 19/50, slacks (295, 224): 121.336508

Track bls12-381-fixed-g1-invert-g2: return Q in G2 with e(P,Q)=z for nonidentity P in G1.

Both-bars alpha (5.15, 0.85) (f 50k 5.085035 + 0.065 margin, g 50k
0.731543 + 0.118 margin with q/(q+1) p = 2 fix), theta 19/50
(1040 rows, margin 0.3014), slacks (295, 224) (collection 2^121.1553),
inputs 64.131, 746.0/498.9 (two-log 2^115.8581), meter 2^116.58.
Ledger <= 2^121.336508, success >= 0.391337. Analytic, not a full-size run.

Entry points: audit_refinement.py, unit_entries.py, evidence16.py.
