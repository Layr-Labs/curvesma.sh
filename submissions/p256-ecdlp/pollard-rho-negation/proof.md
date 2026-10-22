# p256-ecdlp: parallel Pollard rho with the negation map

Status: ready for independent exploratory review; no qualification or emitted
score is claimed by this file. This is a generic best-known-style attack
accounted under `classical-group-work-v2`, not a cryptanalytic advance.
`claim.json` binds `tracks/p256-ecdlp/target-profile.json` and
`tracks/classical-group-work-v2.json` by SHA-256. No auxiliary input, advice,
oracle or quantum step is used.

## Why a generic method

P-256 has prime order r (cofactor 1), r != p (not anomalous), embedding degree
far above any feasible pairing transfer, and a = -3 gives no known efficiently
computable endomorphism other than negation. The only automorphisms available
to a generic walk are {1, -1}. The method below therefore uses the negation map
(factor sqrt(2) in steps) and nothing more.

## Algorithm

Notation: p, r, P, a, b from the protected profile; input Q. Constants of the
run: partition count R = 2048, walk count w = 2^40 (2^30 processor cores each
running a batch of B = 1024 walks), distinguished-point probability
theta = 2^-32, useful-step budget k = 27*2^123, per-walk length
L = k/w + 2^37 = 27*2^83 + 2^37, total steps K = w*L.

1. Input handling. Decode Q canonically, check it lies on the curve and that
   [r]Q = infinity (cofactor 1 makes this a validity check). If Q = infinity
   output 0. If Q = P output 1 (any trivial short-cut is optional; it does not
   affect the bound).
2. Partition table. For i in 0..R-1 draw c_i, d_i uniformly in [0,r) by
   rejection sampling 256-bit words (rejection probability < 2^-32 per draw)
   and compute M_i = [c_i]P + [d_i]Q in affine form. If any M_i is infinity,
   redraw it.
3. Classes. Work on classes {X, -X}. The canonical representative of X=(x,y)
   is (x, min(y, p-y)). The partition index is h(X) = (x >> 8) mod R and
   depends on x only, hence is a class function.
4. Walk. From a class representative X with i = h(X), set Y = canon(X + M_i).
   Look-ahead (Wiener-Zuccherato): if h(Y) = i, replace by
   Y = canon(X + M_{(i+1) mod R}), repeating with the next index while the
   index repeats (at most R tries, then the walk is restarted). This removes
   2-cycles except with probability O(R^-3) per step.
5. Longer fruitless cycles. Every 4096 steps a walk performs a cycle check
   over its next 128 steps: it records the representative with minimal x seen
   in the first 64 of these steps and tests whether that representative recurs
   in the following 64. On recurrence the walk escapes by replacing the
   minimal representative X_min with canon(2*X_min) and continuing. Cycles of
   length <= 64 are therefore escaped; longer fruitless cycles occur with
   negligible probability for R = 2048 (Bernstein-Lange-Schwabe).
6. Start points. Walk j draws seed coefficients (u_j, v_j) uniformly in
   [0,r) and starts from canon([u_j]P + [v_j]Q). Coefficients are not tracked
   during the main walk; only the seed (walk id) is retained.
7. Batching. Each core advances its 1024 walks in lockstep, using Montgomery's
   simultaneous-inversion trick so that one field inversion serves 1024
   affine additions.
8. Distinguished points. A representative is distinguished when the low 32
   bits of x are zero. The DP (x, walk id, step index) is inserted into a
   shared open-addressing hash table of capacity 2*ceil_cap where ceil_cap =
   2*floor(K/2^32) is the DP cap. If another walk id already stored the same x,
   go to step 9 and restart the newer walk from a fresh seed.
9. Collision solving. Re-walk both colliding walks from their seeds to the
   shared DP, now tracking coefficients: each walk state carries (alpha, beta)
   with X = [alpha]P + [beta]Q, updated by adding (c_i, d_i) and negated when
   canonicalisation negates. At the DP, [alpha1]P + [beta1]Q = s*([alpha2]P +
   [beta2]Q) with s = +1 or -1 determined by comparing the actual y. If
   beta1 - s*beta2 is nonzero mod r, output
   x = (s*alpha2 - alpha1) / (beta1 - s*beta2) mod r, verify [x]P = Q and
   return it. At most 4 collisions are processed; after the 4th, or when all
   walks have run L steps, or the DP cap is reached, stop with failure.

## Success argument

Heuristic H1 (random-walk heuristic): on the N = (r+1)/2 classes the R-adding
walk with look-ahead behaves, for collision purposes, like a random mapping up
to a factor (1 - 1/R)^(-1/2) in expected collision time (Teske; Brent-Pollard;
Bernstein-Lange-Schwabe for the negation-map variant). Heuristic H2: fruitless
cycles longer than 64 steps are negligible for R = 2048, as analysed and
measured by Bernstein-Lange-Schwabe.

Consider the first k/w steps of every walk, k = 27*2^123 points in total.
Under H1 the probability that these contain no repeated class is at most
exp(-k(k-1)/(2N)) with k^2/(2N) = k^2/(r+1) = 0.71191... Hence a repeated
class occurs with probability >= 1 - e^(-0.7119) = 0.5093. The R-adding walk
correction (1 - 1/2048) changes this by less than 0.0002. Even if the true
walk behaved 1.4 times worse than random in the collision exponent, the
bound would still be 1 - e^(-0.508) = 0.398 >= 0.39. That leaves a stated
margin for H1.

A repeat between two different walks merges them, since the step function is
defined on classes. Each merged walk still has at least 2^37 steps to reach a
DP; it fails to do so with probability (1 - 2^-32)^(2^37) = e^(-32), so both
reach the same DP with probability >= 1 - 2e^(-32). A repeat inside one walk
needs a self-intersection of a walk of length < 2^89 in a set of 2^255
classes, probability < 2^-60 per walk and < 2^-20 in total by a union bound.
The collision is useless only if beta1 = +-beta2 mod r; seeds and table
coefficients are uniform, so this has probability about 2/r. The DP cap
2*floor(K/2^32) is exceeded with probability below e^(-2^90) by a Chernoff
bound. Hash probing longer than 64 at load factor 1/2 can only drop a DP; this
delays detection with probability < 2^-50 per DP.
Removing these losses (all below 2^-19 in total) gives success probability
>= 0.509 - 2^-19 > 0.39 under H1/H2, over uniform x and fresh coins.

## Instruction ledger (scalar AArch64, 64-bit GPRs)

Field elements are 4 limbs in Montgomery form. For P-256, -p^{-1} mod 2^64 = 1,
and p has words (2^64-1, 2^32-1, 0, 2^64-2^32+1). That makes each reduction
round shift-and-add work plus one MUL/UMULH pair.

| Kernel | Bound (instructions) | Justification |
| --- | ---: | --- |
| Field mult/square M | 256 | 4x4 schoolbook: 16 MUL + 16 UMULH + <= 64 ADDS/ADCS; 4 Montgomery rounds <= 16 each (64); final conditional subtract 12; 8 loads + 4 stores; call/return 8: total 192 |
| Field add/sub A | 32 | 4 ADDS/ADCS, 4 SUBS/SBCS against p, 4 CSEL, 8 loads, 4 stores, call: <= 28 |
| Field inversion I | 76,800 = 300M | Fermat x^(p-2) by a fixed chain of 255 squarings and <= 13 multiplications, plus loop control |
| Unbatched affine add/double | 78,080 = I + 4M + 8A | slope (with inversion), x3, y3, case tests for infinity/opposite/doubling |
| Scalar multiplication | 512 * 78,080 | double-and-add over 256 bits |

Per main-walk step (batched), S = 2304:

| Item | Instructions |
| --- | ---: |
| Montgomery batch trick, 3M per element | 768 |
| Slope multiply, lambda^2, lambda*(x1-x3): 2M + 1S | 768 |
| 6 field subtractions | 192 |
| Amortised inversion I/1024 = 75, plus its bookkeeping | 128 |
| Canonical negation (p-y, compare, select) | 64 |
| Partition index and 8-word table load | 32 |
| Walk state and batch scratch loads/stores | 64 |
| Distinguished test | 8 |
| Look-ahead fallback: probability R^-1/(1-R^-1) of an unbatched add | 64 (needs 38.2) |
| Subtotal 2088; reserve | 216 |

Cycle checks and escapes (step 5) are charged a further 1/16 of every step,
which is 256 instructions per 4096 steps covering 128 compare-and-record
steps and an amortised doubling. The charged per-step cost is
S' = 2304 * 17/16 = 2448 instructions, applied to all K steps, including
the 2^37-step detection tail of every walk.

Other components, all charged in full:

- Setup: 2^24 for decoding, curve/order checks of Q and constants, plus
  2048 * (2 scalar multiplications + 1 addition) for the table and 2048 * 1024
  for its randomness and rejection sampling.
- Walk starts: per walk, 2 scalar multiplications, 1 addition and 1024 for
  seed randomness.
- DP handling: per DP, 2^14 covering hashing, at most 64 probes of a 64-byte
  record compare (16 instructions/word + 64), and insertion; charged for the
  full cap of 2*floor(K/2^32) DPs.
- Re-walks: at most 4 collisions, two walks each of at most L steps,
  unbatched with coefficient tracking at <= 2^17 per step (78,080 plus two
  mod-r additions, negations and look-ahead).
- Solving, verification and output: 2^28 covering mod-r inversion, the [x]P=Q
  check, and serialization.

## Submitted computation bound

| Component | Instructions |
| --- | ---: |
| main walk K * 2448 | 702853228875198766214901053726844700852224 |
| setup | 163924410368 |
| walk starts | 87997240494075674624 |
| DP handling | 2190500237194381558707807745736704 |
| re-walks | 273812529649297694838475968217088 |
| solving, verification, output | 268435456 |

K = 287113247089541979662949776849201266688 steps; L =
261127977036760039175487488.

Total T = **702853231339511533058668304513786683326464** instructions.
Dividing by the fixed C = 1000 gives log2(T/1000) = 129.0464910...,
which rounds up to submitted time_log2 **129.046492**. Preprocessing (setup,
independent of later walks) is 163924410368 instructions, giving
preprocessing_log2 **27.288456**, included in the total. All point arithmetic is
charged through its instruction count; no addition is priced as one unit.

## Memory and advice

The DP hash table has capacity 2 * 133697524242821140057849593856 records of
64 bytes (x, walk id, step index). Walk states, batch scratch and seeds take
2^40 * 128 bytes. Code, the 2048-entry table and constants fit in a 2^30-byte
reserve, covered by the 2^48 term. Peak memory <= 2*ndp*64 + 2^48 =
17113283103081106208879724724224 bytes, memory_log2_bytes **103.754888**.
Nonuniform advice: none (2^0 bytes).

## Evidence and limitations

The cost and success arguments are analytic. `experiment.py` runs the same
negation-map walk, with look-ahead and distinguished points, on a toy
prime-order curve. It is evidence only at that scope. H1 and H2 are standard,
well-supported heuristics but are not theorems. The instruction bounds are
analytical upper bounds, not Callgrind measurements of a P-256 kernel.

References: van Oorschot-Wiener, "Parallel collision search with cryptanalytic
applications" (J. Cryptology 1999); Wiener-Zuccherato, "Faster attacks on
elliptic curve cryptosystems" (SAC 1998); Teske, "Speeding up Pollard's rho
method for computing discrete logarithms" (ANTS 1998); Bernstein-Lange-Schwabe,
"On the correct use of the negation map in the Pollard rho method" (PKC 2011).
