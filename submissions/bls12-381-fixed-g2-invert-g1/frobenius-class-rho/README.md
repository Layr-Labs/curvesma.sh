# bls12-381-fixed-g2-invert-g1: Pollard rho in GT on Frobenius classes

Compute g = e(P0, Q), find x = log_g z by distinguished-point rho on the 12-element Frobenius classes of GT, return [x]P0.
Claimed total work 2^126.19 dense Fp12 multiplications at success probability at least 0.39 and peak memory 2^53.7 bytes (heuristic walk model, named in proof.md).
See proof.md for the ledger and claim.json for the bound; experiment.py is a toy-group check of the scheme only.

Revision 2 (2026-10-22): resubmitted unchanged in substance after the challenge's evaluation workflow could not be dispatched during a GitHub Actions outage; the claim, proof and score are identical to revision 1.
Revision note (2026-10-22 16:40Z): resubmitted after another evaluation-dispatch failure during the GitHub Actions outage; substance unchanged.
Revision 3 (2026-10-23): the cycle escape is now a squaring of the state at the cycle's least-key member (a table-step escape can be undone by the next step); the cycle-waste factor, the source of the index and distinguished-test bits (now a charged mix of the key) and the memory bound (now 2^53.7 bytes) are rederived in proof.md. The score 126.19 is unchanged.
