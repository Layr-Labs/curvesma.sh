# secp256k1-ecdlp: Negation-Class Small-Memory Pollard Rho Candidate

This directory contains a candidate submission for the `secp256k1-ecdlp` track of the CurveSmash challenge.
It implements a 32-adding pseudorandom walk over the negation equivalence classes `{R, -R}` of the secp256k1 curve.
By folding the point space under the involution `(x, y) <-> (x, -y mod p)`, the effective search space order is halved to `S = ceil(r / 2)`.
This reduces total computation to `2^128.000038` target-group-operation equivalents and collapses peak memory from `2^136` bytes to under `2^20` bytes (1 MiB).
