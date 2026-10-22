# bls12-381-fixed-g1-invert-g2: GT rho on Frobenius classes

Given P in G1 and z in GT, compute g = e(P,Q0), solve z = g^x by parallel Pollard rho on the 12-element Frobenius classes of GT, and return Q = [x]Q0.
Claimed work 2^126.16 target-group-operation equivalents at success probability at least 0.39.
- [Algorithm and cost argument](proof.md)
- [Submission claim](claim.json)

Revision 2 (2026-10-22): resubmitted unchanged in substance after the challenge's evaluation workflow could not be dispatched during a GitHub Actions outage; the claim, proof and score are identical to revision 1.

Revision 3 (2026-10-23): the cycle rule's escape is now a squaring of the smallest cycle member (class of m^2) instead of an ordinary table step with index j+1, which the next step could undo; the per-event cost, the re-entry bound and the toy experiment are updated. The claimed score is unchanged at 126.16.
