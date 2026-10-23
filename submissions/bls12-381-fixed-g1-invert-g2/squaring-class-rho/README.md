# Squaring-class mixed add/square rho

Candidate for `bls12-381-fixed-g1-invert-g2`. The supplied `P` and `z` reduce
exactly to one discrete logarithm in `GT`: solve `z = e(P,Q0)^a` and return
`Q = [a]Q0`.

The walk replaces one dense Fp12 multiplication by a block of one dense
multiplication plus three cyclotomic squarings (`1.75` reference units) that
emits eight group elements: `X, X^2, X^4, X^8` and their free conjugates
`X^-1, X^-2, X^-4, X^-8`. Every emitted element is stored and compared, so the
fourteen disjoint ratio classes `{+/-2^d : d in [-3,3]}` between two blocks
are fourteen equalities the walk actually visits. The budget is set so that
the birthday bound gives success probability `0.45`, above the `0.39` floor.

Declared score: `time_log2 = 126.4911`, memory `2^134.402` bytes,
preprocessing `2^20`, advice `0`.

See `proof.md` for the reduction, the disjointness and solvability proof, the
complete instruction ledger, and the priced-out full Frobenius orbit.
`experiment.py` checks the collision counting and the birthday model at small
scope in a prime-order subgroup of `F_p^*`.
