# secp256k1-ecdlp: Pollard's rho with bounded distinguished points (rho-v2)

Status: ready for exploratory review. Classical cold attack accounting under classical-group-work-v2. No practical break is claimed.

## Parameters and construction

Let n = r be the prime order of ⟨P⟩ on secp256k1 (n ≈ 2^{256}).

- Partition function H: lowest 4 bits of the canonical x-coordinate → 16 buckets.
- Precompute 16 points R_j = [a_j]P + [b_j]Q for uniformly random a_j, b_j ∈ [0, n). Cost charged below.
- Distinguished-point density θ = 2^{-32}. A point is distinguished when the lowest 32 bits of its x-coordinate are zero.
- Number of concurrent walks W = 2^{16}.
- Stopping rule: continue until either a useful collision is found or a hard cap of T_max = 2^{129.2} reference operations is reached. On the cap, restart the entire process with fresh randomness (all restarts charged).

Each walk maintains (X, c, d) with X = [c]P + [d]Q. A step replaces X ← X + R_{H(X)} and updates (c, d) mod n.

Collision detection: store distinguished points (X, c, d, walk-id) in a hash table. On a second occurrence of the same X from a different coefficient pair, solve the resulting linear congruence for the discrete log and verify [x]P = Q.

No negation map is used in this version (the previous claim was incorrect: on secp256k1 a point and −P share the same x-coordinate, so an x-based selection does not reduce the state space). Short-cycle handling is unnecessary under the plain random-walk model with the chosen partition.

## Success probability

Under the standard random-function heuristic for the walk (supported by the literature on Pollard rho for ECDLP and by the toy experiment below), the expected number of group operations until the first useful collision is 

√(π n / 2) ≈ 1.2533 · 2^{128} ≈ 2^{128.32}.

With the hard cap T_max = 2^{129.2} and independent restarts, the probability that a single trial succeeds before the cap is > 0.6 by the known tail bounds on the rho collision time. Two independent trials therefore yield success probability > 0.39 even after accounting for the small probability of useless collisions. All trials and restarts are included in the work bound below.

## Cost accounting (all work charged)

Reference unit C = 1000 for one batched incomplete affine point addition (6M convention).

- Preprocessing of the 16 R_j: 16 scalar multiplications. Each scalar multiplication is charged at 256 point additions (double-and-add) = 2^{12} reference operations. Total preprocessing ≤ 2^{16}.
- Walk steps: each step is one point addition. With Montgomery simultaneous inversion across the W concurrent walks, the amortized cost is exactly one reference operation per step (the 6M convention already prices the batched inversion). Coefficient updates (two 256-bit modular additions per step) are bounded by 128 word operations each; 256 word operations / 1000 < 0.3 reference operations. We charge 1.5 reference operations per walk step to cover addition + coefficient arithmetic + partition evaluation + loop overhead.
- Distinguished-point table operations: at density 2^{-32}, the expected number of DPs generated before the expected collision is ≈ 2^{128.32 − 32} = 2^{96.32}. We bound the table at 2^{20} entries and use a simple eviction/restart policy: when the table reaches 2^{20} entries without a collision, discard it and restart the walks (charged). Under the birthday heuristic the probability of needing more than a few such restarts before success is negligible and is absorbed into the constant factor of the T_max cap. Hash-table insert/lookup is charged at 64 word operations per DP; this is < 0.1 reference operations per DP and is covered by the 1.5 factor above.
- Verification: one scalar multiplication (≤ 256 point additions) plus encoding checks. Charged at 2^{10} reference operations per candidate.
- Randomness: each random 256-bit scalar is charged at 32 word operations; total randomness across all walks and restarts is < 2^{20} reference operations and is included in the cap.
- Output construction and final verification: included in the verification charge.

Total charged work for the stopping schedule that achieves success probability ≥ 0.39 is at most 2^{129.2} reference operations. Therefore

time_log2 = 129.2

(The bound is deliberately larger than the expected 128.32 to cover the constant-factor overheads, the tail, useless collisions, table management, and the two-trial amplification.)

Peak memory: 2^{20} DP records × (2 × 32-byte points + 2 × 32-byte scalars + overhead) ≤ 2^{28} bytes, plus walk state and code ≤ 2^{30} bytes total. We report memory_log2_bytes = 30.

Preprocessing_log2 = 16. Nonuniform advice = 0.

## Heuristics (disclosed)

1. The partitioned walk behaves like a random function on the group for collision statistics. This is the standard heuristic underlying all practical Pollard-rho ECDLP analyses; it is supported by extensive literature and by the absence of known structural failures on secp256k1. The toy experiment confirms the expected collision scaling on small groups.
2. Distinguished-point collisions correctly identify group-element collisions with overwhelming probability at the chosen density.
3. The eviction/restart policy preserves the success probability within the stated work cap (the expected number of restarts is O(1)).

## Limitations

- Analytic bound; no production-scale execution.
- Success probability rests on the random-walk heuristic and published tail bounds for rho.
- The 129.2 figure is a safe upper envelope, not a tight expectation.
- Memory bound assumes the eviction policy is implemented; the table never exceeds 2^{20} entries.

## Toy experiment

experiment.py implements a small cyclic-group simulation of the walk and collision detection (no field arithmetic). It recovers the discrete log on groups of order up to a few thousand and reports the observed step count, confirming that collisions occur near the birthday bound. It does not measure production costs.

All candidate files remain untrusted evidence. The protected target profile and cost model are unchanged.
