# bls12-381-fixed-g1-invert-g2: GT rho on Frobenius classes

Given P in G1 and z in GT, compute g = e(P,Q0), solve z = g^x by parallel Pollard rho on the 12-element Frobenius classes of GT, and return Q = [x]Q0.
Claimed work 2^126.48 target-group-operation equivalents at success probability at least 0.39; preprocessing 2^45.1 (per-processor tables, included); memory 2^74 bytes.
- [Algorithm and cost argument](proof.md)
- [Submission claim](claim.json)
- [Ledger recomputation](ledger.py)

Revision 4 (2026-10-23): table index and distinguished-point test come from disjoint low bit fields of the canonical coordinate with an explicit statistical-distance bound; the walk-randomness
margin is an explicit 1.20 on steps over the additive-walk penalty sqrt(1+1/R), with a toy experiment of the exact walk (steps to collision with a confidence interval); cycle probabilities are bounded for every length with a closure term;
preprocessing is one figure (2^45.1) everywhere; the declared memory covers 80-byte records at hash load 1/2 and 704-byte window states (2^74 bytes); a 1.07 reserve is added and the score is 126.48.
