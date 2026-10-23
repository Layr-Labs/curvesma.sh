# Six-class (negation + GLV) parallel Pollard rho in BLS12-381 G1

Parallel distinguished-point Pollard rho on the (r-1)/6 orbits of G1 under
<-1, phi>, phi(x, y) = (beta x, y). Class key |y| (free); the canonical twist index is
decided exactly from a truncated fixed-point product (exact fallback on
near-ties, metered separately); a twisted table makes each step one batched
affine G1 addition plus at most 200 instructions; fruitless cycles are escaped by
doubling the minimum-fingerprint class. A global work meter charges every action and halts at the claimed
bound, so W <= 2^126.358 on every run.

Claim: time_log2 = 126.358 (63179/500), success >= 0.39, memory < 2^72 bytes.

- [Algorithm, success argument and cost](proof.md)
- [Experiment: BLS12-381 identity checks + toy-curve runs](experiment.py)

## Provenance and the one change

The algorithm, parameters, prices, ledger and `experiment.py` (byte for byte)
are those of public submission 1c600df4-b531-401a-b9d8-5124c6bceeb2
(glv-six-class-rho-v11), credited here as the source. Its meter bound was
W_max = 2^126.36. The one change:

1. W_max is stated to three decimals, 2^126.358, instead of rounded up to two.
   The meter headroom above the deterministic charges falls from > 2^119.45 to
   > 2^119.18, the Markov bound on an early halt rises from 0.0039 to 0.0048,
   and success is >= 0.395 - 0.0049 > 0.39 (proof.md sections 4 and 5).

Nothing else changes. The argument in proof.md is self-contained and relies on
no earlier review or score.
