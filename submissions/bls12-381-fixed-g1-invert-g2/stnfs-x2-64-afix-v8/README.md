# STNFS with alpha (5.25, 0.844), theta 762/2000, slacks (294, 224): 121.343021689718

Track bls12-381-fixed-g1-invert-g2: return Q in G2 with e(P,Q)=z for nonidentity P in G1.

Exact 400k alpha sweeps (f 5.114285, g 0.760075 + 0.13/0.08 tails, 12 spot
checks; charges (5.25, 0.844)), theta 759/2000->762/2000
(1052 rows, margin 0.3022), slacks (294, 224) (collection 2^121.1626),
inputs 64.132, 746.1/498.9 (two-log 2^115.8593; 4 plan moves), meter 2^116.582.
Ledger <= 2^121.343021689718, success >= 0.391658. Analytic, not a full-size run.

Entry points: audit_refinement.py, unit_entries.py, evidence16.py.
