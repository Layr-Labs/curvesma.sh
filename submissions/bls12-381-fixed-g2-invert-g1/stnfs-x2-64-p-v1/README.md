# Recover P via the GT STNFS core with 99.9% ECM, X2 = 2^64, a measured record route and fully reduced LA kernels: 121.4953522

Track bls12-381-fixed-g2-invert-g1: return P in G1 with e(P,Q)=z for the supplied nonidentity Q in G2.

Read proof.md section 1 (reduction g = e(P0,Q), P = [a]P0, field-basis conversion, wrapper), then the core (sections 2-10, STORAGE.md, DESCENT.md, REFINEMENT.md). Core: 711418d3 from the recover-Q track (relation collection 2^121.155782 with every kernel price rounded up, linear algebra 2^117.732694, two-log meter 2^118.617). Upper-rounded ledger 2^121.4953522, success >= 0.392133. Analytic, not a full-size run.

Entry points: audit_refinement.py, unit_entries.py, evidence16.py.
