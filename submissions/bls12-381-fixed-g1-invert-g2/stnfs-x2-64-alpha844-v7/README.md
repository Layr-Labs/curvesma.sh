# STNFS with corrected alpha, theta 759/2000, slacks (294, 224): 121.327036929110

Track bls12-381-fixed-g1-invert-g2: return Q in G2 with e(P,Q)=z for nonidentity P in G1.

Review-driven alpha fix: p = 2 double-root tails corrected to exact primitive
2q(1-q^-7)/(q^2-1) (primitivity included), charge alpha_g 0.844 (50k partial
0.731625, 15% margin); theta 759/2000 (1034 rows, margin 0.3007, allowance
1.0674), slacks (294, 224) (collection 2^121.1452), inputs 64.130, 746.0/499.0
(two-log 2^115.8429, meter 2^116.5647). Ledger <= 2^121.327036929110, success
>= 0.391369. Analytic, not a full-size run.

Entry points: audit_refinement.py, unit_entries.py, evidence16.py.
