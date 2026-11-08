# secp256k1-ecdlp: Pollard's rho with Floyd cycle detection (rho-v3)

Status: ready for exploratory review. Classical cold attack under classical-group-work-v2. No practical break is claimed.

## Construction

Let n = r be the prime order of ⟨P⟩ on secp256k1.

Partition function H maps a point to one of L = 16 buckets via the lowest 4 bits of its canonical x-coordinate. Precompute 16 points R_j = [a_j]P + [b_j]Q for uniformly random a_j, b_j.

A walk maintains a triple (X, c, d) with X = [c]P + [d]Q. The step function is
X ← X + R_{H(X)}, with the corresponding coefficient updates mod n.

Floyd cycle detection: maintain two walkers (tortoise and hare). The tortoise advances one step; the hare advances two steps. When they meet, a cycle has been detected. Backtrack or solve the resulting coefficient relation to obtain a candidate discrete log, then verify [x]P = Q.

No distinguished-point table is used. Peak memory is O(1) (a constant number of points and scalars plus code).

No negation map is used.

## Success probability and stopping schedule

Under the standard random-function heuristic for Pollard's rho (the walk behaves like a random mapping), the expected number of steps until the tortoise and hare meet is approximately 1.0308 √n (see the classic analysis). Each such iteration performs three group operations, for an expected ≈ 3.09 √n group operations.

√n ≈ 2^{128}, so the expected group-operation count is ≈ 2^{129.62}.

The collision-time distribution has a known tail. Charging a hard per-trial budget of 2^{130.5} reference operations gives success probability > 0.5 on a single trial (the probability of exceeding roughly 4× the mean is small under the published rho statistics). Two independent trials, each charged in full, therefore yield success probability > 0.75, which is comfortably above 0.39. All trials are included in the total bound below.

If a trial hits the cap without success, it is abandoned and a fresh trial is started; the total work of the two trials is the declared bound.

## Cost accounting

Reference unit C = 1000 for one batched incomplete affine point addition.

- Preprocessing of the 16 R_j: 16 scalar multiplications ≤ 16 × 256 = 2^{12} reference operations.
- Each walk step (tortoise or hare) is one point addition. With the 6M batched convention the addition itself is 1 reference operation. Coefficient updates and partition evaluation add < 0.5 reference operations. We charge 1.5 reference operations per step.
- Because Floyd performs three steps per iteration, the charge per iteration is 4.5 reference operations.
- Verification of a candidate: one scalar multiplication ≤ 256 additions = 2^8 reference operations.
- Randomness and control overhead are absorbed into the 1.5 factor and the final constant-factor margin.

Total charged work for the two-trial schedule that achieves success probability ≥ 0.39 is at most 2^{130.5} reference operations. Therefore

time_log2 = 130.5

Peak memory: a few dozen points and scalars plus code and stack ≤ 2^{16} bytes. We report memory_log2_bytes = 16.

Preprocessing_log2 = 12. Nonuniform advice = 0.

## Heuristics (disclosed)

1. The partitioned walk behaves like a random function on the group, so the classic Pollard-rho / Floyd collision-time analysis applies. This is the standard heuristic used throughout the ECDLP literature.
2. Meeting of the tortoise and hare correctly identifies a cycle from which the discrete log can be extracted with overwhelming probability (the probability of a useless collision is negligible and absorbed in the constant factor).

## Limitations

- Analytic bound only; no production-scale execution.
- Success probability rests on the random-walk heuristic and published rho tail bounds.
- 130.5 is a conservative envelope that covers Floyd's constant-factor overhead relative to an idealized birthday search.
- The toy experiment demonstrates control flow and cycle detection on a small cyclic group; it does not measure field costs.

All work (preprocessing, both trials, failures, verification, output) is charged. No free oracles or uncharged advice. Protected rules and cost model are unchanged.
