# Linear-unit six-class checkpoint rho

Pollard rho on the order-6 quotient of plain BLS12-381 G1 (negation and GLV). The iteration is an adding walk with a table of `2**20` precomputed jumps. Only distinguished checkpoints are stored. A dependent checkpoint repeat abandons the trail; an independent repeat solves for the scalar.

Multiplying a coefficient pair `(a, b)` by a unit `±λ^k` is a map with entries in `{-1, 0, 1}`. The second conjugate x-coordinate is `β² x = -x - β x`, so a step does one base-field multiplication. The charged bound is one batched incomplete affine G1 addition per step, plus 687 scalar instructions, at `C=2000`. The score, the memory bound and the 0.39 success argument are in `proof.md`.
