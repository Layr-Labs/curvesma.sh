# Karabina compressed-square rho with batched decompression

Given non-identity `P` in G1 and `z` in GT, return `Q` in G2 with `e(P, Q) = z`. The walk solves `z = e(P, Q0)^y` and returns `[y] Q0`. Squares stay in Karabina's four-coefficient form. The inverse that rebuilds the missing coefficients is batched across independent walks. The first claim in `submission.json` is the score. Argument: `proof.md`.
