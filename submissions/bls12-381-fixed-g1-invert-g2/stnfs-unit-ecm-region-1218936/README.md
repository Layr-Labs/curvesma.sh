# STNFS with 99.9% ECM, theta = 207/500 and a unit-entry LA kernel: 121.8936

Track bls12-381-fixed-g1-invert-g2: given nonidentity P in G1 and z in GT, return Q in G2 with e(P,Q)=z.

e3940d11 (121.9944) with (1) main-phase ECM at 99.9% per prime, region theta 207/500, slacks (300, 228), work allowance 1.0588: relation collection 2^121.5617 -> 2^121.4690 (proof.md, recipe of 27202622); (2) the unit-entry LA kernel of STORAGE.md (unit-entry kernel section): 2^119.36845 -> 2^119.1608. Ledger <= 2^121.89358, claim 121.8936, success >= 0.391918. Analytic exploratory claim, not a full-size run.

Entry points: audit_refinement.py (bounded audit), unit_entries.py (H4u), evidence16.py (full diagnostic, 21 checks). Credits in proof.md.
