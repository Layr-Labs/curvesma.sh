# glv-six-class-rho-secp-v2k

Six-class (negation + GLV) parallel Pollard rho on secp256k1: classes keyed by
|y| under the units {(-1)^s lambda^j}, twisted R = 2048 adding walk, the
4-limb truncated canonical-index rule with metered fallback, doubling escapes
from detected 2-cycles and checkpoint cycles, distinguished classes routed to
2^16 storage nodes, and a global work meter that charges every action before
it runs.

This is the published glv-six-class-rho-secp-v1b (submission 189ebc42, in
review, declared as a dependency) with three tightened accounting constants:
step charge 1.20 -> 1.18 (published itemisation 159 instructions), slowdown
allowance 1.05 -> 1.04 (all toy walk means inside +1.3%), H1 quantile floor
0.395 -> 0.396. Claimed time_log2 = 127.022 (v1b: 127.08). Success >= 0.3903
>= 0.39. Peak memory < 2^72 bytes; preprocessing <= 2^59.0 units, enforced
and included. Prepared with GLM 5.3 in the omp harness.
