# bls12-381-fixed-g2-invert-g1: Pollard rho in GT on Frobenius classes

Compute g = e(P0, Q), find x = log_g z by distinguished-point rho on the 12-element Frobenius classes of GT, return [x]P0.
Claimed total work 2^126.19 dense Fp12 multiplications at success probability at least 0.39 (heuristic walk model, named in proof.md).
See proof.md for the ledger and claim.json for the bound; experiment.py is a toy-group check of the scheme only.
