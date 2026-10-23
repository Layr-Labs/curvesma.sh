# bls12-381-fixed-g2-invert-g1: truncated BSGS with a fingerprint join

Selected package for the fixed-G2 pairing-inversion track. `submission.json`
carries one scored claim, `127878612/1000000 = 2^127.878612` target-group
operations; `claim.json` carries the same quantities as machine-readable
metadata; `proof.md` holds the reduction, the fixed-width representations, the
algorithm, the correctness and success arguments, the complete work ledger and
the peak-storage bound.

The attack reduces `e(P,Q)=z` for the supplied nonidentity `Q` to the discrete
logarithm of `z` to base `g=e(P0,Q)` in GT, searches the truncated grid
`a=im+j` with `m^2 >= 2r/5`, and joins the two `m`-element lists with a
randomized 256-bit fingerprint radix sort. Record references, radix counters
and prefix sums are explicit 128-bit two-limb values, and every fingerprint
match is confirmed by an exact comparison of all 72 words of the two stored GT
elements. Success is the exact covered-exponent count `m^2/r > 0.4`, less a
`2^-38` Markov term. No heuristic collision constant is used.
