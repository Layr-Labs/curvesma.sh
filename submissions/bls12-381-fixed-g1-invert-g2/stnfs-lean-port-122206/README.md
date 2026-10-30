# stnfs-lean-port-122206

STNFS for bls12-381-fixed-g1-invert-g2 (port of sister-track cbab3319; reduction g=e(P,Q0), output Q=[a]Q0): W <= 2^122.2065 reference units (score 122.2065); ledger 2^122.205983; metered descent cap 2^118.886; norm-bound charge 1,508.2 per block evaluation; per-entry LA charges 88/97; ECM cofactor slacks 291 (f) and 222 (g) bits, with the per-survivor ECM charge taken from the explicit ecm16 schedule at 99% per prime (2^38.5163); descent relation input LOG2_RELS = 64.224, the floored model retained count, checked in [4].
proof.md is the entry point (DESCENT.md, STORAGE.md, KERNELS.md); `evidence16.py`
(stdlib only) recomputes every figure in the verifier's sandbox.
