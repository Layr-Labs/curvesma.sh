# bls12-381-fixed-g1-invert-g2: GT rho on Frobenius classes

Given P in G1 and z in GT, compute g = e(P,Q0), solve z = g^x by parallel Pollard rho on the 12-element Frobenius classes of GT, and return Q = [x]Q0.
Claimed work 2^126.16 target-group-operation equivalents at success probability at least 0.39.
- [Algorithm and cost argument](proof.md)
- [Submission claim](claim.json)
