# Inversion-batched parallel double-or-add walk

Deterministic-budget Pollard walk on secp256k1. 4096 walks advance together so each round's denominators are inverted as one Montgomery batch and each step is priced as one reference batched incomplete affine addition (C=1000) plus a short instruction residual. The walk stops at the birthday budget S(S-1) >= 2r, which is the acceptance floor 0.39 under the random-mapping heuristic, not a full expected-time discrete log. No negation map and no GLV endomorphism.
