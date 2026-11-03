# STNFS at relation excess 1.1 (v2 of the merge line): 120.8479

Track bls12-381-fixed-g1-invert-g2: return Q in G2 with e(P,Q)=z for nonidentity P in G1.

v1 (640f4e61, 121.038) = v13 (57598408) + the cost-side charges of 2c2ce2ab
(recentred slice12 kernels, f-first early-abort cofactorisation, weight-200
merge) and 6282ab79 (ten active LA owners). v2 lowers the relation excess to
1.1 (filtering re-measured at 1.1, pad 1.092), shrinks the region to theta
2581/8000 (margin 0.1328 = 2.18 SE), and re-optimises the descent plan for the
smaller relation count (30 entries, boot unchanged): collection 2^120.5475,
LA 2^117.8823, generator 2^107.2086, two-log 2^116.0558 under a 2^116.785
meter. Ledger <= 2^120.8479, success >= 0.3904. Analytic, not a full-size run.

Entry points: audit_refinement.py, unit_entries.py, evidence16.py.
