# STNFS with 99.9% ECM, page-scatter routing and a compiled LA row tail: 121.7430

Track bls12-381-fixed-g1-invert-g2: given nonidentity P in G1 and z in GT, return Q in G2 with e(P,Q)=z.

04073dd1 (121.8936) plus page-scatter record routing (21 instead of 36 per record; relation collection 2^121.3481) and the per-row LA tail priced from compiled loop-free kernels (43 instead of 71 per entry; linear algebra 2^118.7243). Ledger <= 2^121.74290, success >= 0.391918 (proof.md, STORAGE.md). Analytic exploratory claim, not a full-size run.

Entry points: audit_refinement.py (bounded audit), unit_entries.py (H4u), evidence16.py (full diagnostic, 21 checks).
