# Six-class checkpoint rho

Pollard rho on the order-6 quotient of plain BLS12-381 G1 (negation and GLV). The iteration is an adding walk with a table of `2**20` precomputed jumps. Only distinguished checkpoints are stored. A dependent checkpoint repeat abandons the trail; an independent repeat solves for the scalar.

The charged bound is one batched incomplete affine G1 addition per step, plus the scalar instruction count in `reference.py`, at `C=2000`. The score, the memory bound and the 0.39 success argument are in `proof.md`.
