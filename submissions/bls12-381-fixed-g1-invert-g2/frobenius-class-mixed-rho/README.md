# Frobenius-class mixed add/square rho on GT

Candidate for `bls12-381-fixed-g1-invert-g2`. Fixed-argument pairing inversion is
reduced to one discrete logarithm in GT with base `g = e(P, Q0)`; the returned
`Q = [x]Q0` satisfies `e(P, Q) = z` by bilinearity.

The walk is a full-table birthday search in GT. One dense Fp12 multiplication
advances the walk, `D` cyclotomic squarings produce `X^{2^d}`, and five Frobenius
maps produce `X^{p^k 2^d}` for `k = 1..5`. Because `p^6 = -1 (mod r)`, conjugation
is inversion on GT and is free, so each block emits the twelve-element Frobenius
orbit `{X^{eps p^k 2^d}}` at no extra cost. Every emitted element is inserted and
compared, so a collision is a visited value equality, not a ratio class: it gives
a linear equation in `x` that is solved with one modular inverse.

Files:

* `proof.md` - reduction, cheap-map costs, collision algebra, birthday count and
  the complete accounting of the submitted bound.
* `experiment.py` - toy replay of the emission and collision algebra in `Z_q`,
  plus exact checks of `p^6 = -1 (mod r)` and of the emitted multipliers on the
  real profile constants.
* `submission.json` - claim and experiment declaration.

Submitted `time_log2 = 125.30`, success probability `0.45`, memory
`2^135.1` bytes, preprocessing `2^8.0` target-group-operation equivalents,
advice `0` bytes.
