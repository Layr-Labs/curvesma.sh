# inversion-class-rho-v1: parallel distinguished-point rho in GT on the inversion quotient

Given nonidentity `P` in G1 and `z` in GT, compute `g = e(P,Q0)` once, solve
`x = log_g z` by a parallel distinguished-point Pollard rho **on the quotient of
GT by `y -> y^-1`** (the `p^6`-Frobenius, which is inversion because
`p^6 = -1 mod r`), return `Q = [x]Q0` and verify with a pairing.

The quotient has `N = (r-1)/2 = 2^253.857089` classes, so the square-root search
costs `sqrt(2)` less than rho in GT. The class key is `K(y) = |Re a1|` (the `Fp`
component of the `w^1` coefficient of `y`), which is invariant under inversion
and is computed without inverting anything; canonical representatives are needed
only when storing a distinguished point, and records are compared **exactly**
(576-byte Fp12 elements), so `K` never produces a false match.

Hard cap: `W_max = 1.5 sqrt(N)` steps at `1.01` target-operation equivalents per
step, plus `2^16` reference operations for the pairings and G2 scalar
multiplications. Claimed `log2(W) = 127.527863`; success probability `0.39` (the
rule threshold; the random-mapping model gives 0.675 at this cap, and the toy run
in `experiment.py` measures 0.375 at the same relative budget). Peak storage
`<= 2^78` bytes; preprocessing charged in full with no reuse credit; nonuniform
advice 0; no auxiliary powers.

`proof.md` has the reduction, the class-walk design, the instruction ledger and
the success analysis; `experiment.py` is standard-library only and checks the
BLS12-381 facts (Part A) plus a bounded toy run of the complete algorithm
(Part B). Credit: the class-rho idea for this track pair is public (sibling-track
submissions `3da91365`, `41c34165`, `7da5f23a`, full 12-class version); this is
an independent 2-class implementation with its own accounting.
