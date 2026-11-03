# Recover P by STNFS on a directly sampled region: 121.400000

Track bls12-381-fixed-g2-invert-g1: return P in G1 with e(P,Q)=z for nonidentity Q in G2.

Read proof.md (claim, reduction, ledger, success), then DESCENT.md, STORAGE.md,
KERNELS.md and REFINEMENT.md. Initial experiment: `audit_refinement.py` (empty
input, a few seconds). Optional: `evidence16.py` {"samples": 16, "seed": 1}
recomputes every figure (about 10 s natively); `unit_entries.py` reruns the
+-1 entry share. Analytic exploratory claim under H1-H5, not a full-size run.
