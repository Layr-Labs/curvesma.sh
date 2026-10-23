# GT Frobenius-class walk: shared specification (revision v9)

Used by `curvesmash-core/cscore/gt.py` (toy model, Python), the blst step kernel
(`curvesmash-core/examples/blst-gt-rho/step.c`) and its Python replica. Written
2026-10-22; revision v7 2026-10-23 (sections 6, 8, 9; the distinguished-point table and the
classification are now specified in full; revision v8 (2026-10-23) makes each trail end at its first DP and
restart at the walk's start element, section 6 and 9). Any change here must be applied to all three.

## 1. Group and Frobenius

- G = cyclic group of prime order m (production: m = r, the BLS12-381 subgroup order;
  toy: m = ℓ, a prime dividing Φ12(q) = q^4 − q^2 + 1 with ℓ ∤ 12).
- σ(x) = x^q, with q = p in production and a small prime in the toy model. λ = q mod m.
  Then σ(g) = g^λ, λ has order 12 and λ^6 = −1, so σ^6(x) = x^{-1}.
- Frobenius conjugates of x: σ^n(x) for n = 0..11. The orbit has 12 elements when x ≠ 1.
- Classes: the orbits. N = (m − 1)/12 classes (the identity excluded).

## 2. Scalar value and key

- For x in GT, V_n(x) is a fixed scalar of σ^n(x) (production: the real part of the
  w^1 coefficient in the blst basis Fp2[w]/(w^6 − ξ), ξ = 1+u; toy: the base-q integer
  encoding of the coefficient vector of σ^n(x) in F_q[x]/f).
- Key: K(x) = min_n V_n(x) and k(x) = argmin_n V_n(x) (smallest n on ties).
- K and the orbit are class invariants: K(σ^j x) = K(x). Ties are exceptional (probability
  about 2^-250 per step in production; they are resolved by smallest n and verified by the
  replica).
- Production implementation of the 12 values: with a = Re c1 and B = Re(conj(c1)·d1),
  where d1 = γ_{1,1} = ξ^{(p−1)/6}, the 12 values are ±a·ζ^j and ±B·ζ^j for ζ a primitive
  sixth root of unity in Fp, j = 0, 1, 2. Four Fp multiplications (B: two; aζ and Bζ: one
  each; aζ² = aζ − a and Bζ² = Bζ − B need none). The replica checks this formula against
  direct computation of V_n(σ^n(x)) on real GT elements.

## 3. Exponent invariant

- A walk state is (x, a, b) with x = g^{a + b·y} (mod m), where z = g^y is the target and
  g generates G (production: g = e(P, Q0), z = e(P, Q0)^y).
- The state x need NOT be the canonical element. The canonical element is σ^k(x) =
  g^{A + B·y} with (A, B) = (λ^k a, λ^k b), where k = k(x).
- Two states whose canonical elements are equal give y = (A2 − A1)/(B1 − B2) mod m, when
  B1 ≠ B2. Every candidate y is checked by g^y = z before it is output.

## 4. Twisted adder table (no Frobenius map per step)

- r_a = 2^adder_bits adders. Adder h has exponents (c_h, d_h), chosen at random, with
  R_h = g^{c_h} z^{d_h}. Adder elements are built by the simulation shortcut in the toy
  model; the walk never reads y.
- For k in 0..11: T[h][k] = σ^{-k}(R_h) = g^{λ^{-k} c_h} z^{λ^{-k} d_h}. Store the pairs
  E_c[h][k] = λ^{-k} c_h mod m and E_d[h][k] = λ^{-k} d_h mod m, with the element.
  Table size 12 · r_a elements (production: 12 · 1024 · 576 bytes ≈ 7 MiB).

## 5. One step (mixed add/square walk)

Let S = round(s · r_a), with s ∈ [0, 1) the square probability. For the add-only walk
s = 0.

1. Compute K = K(x) and k = k(x) (section 2).
2. h = K mod r_a.
3. If h < S (square slot): x ← x²; (a, b) ← (2a, 2b). Squaring commutes with σ, so no
   twist is needed. Production uses blst_fp12_cyclotomic_sqr (charged at its instructions).
4. Else (add slot): x ← x · T[h][k]; a ← a + E_c[h][k]; b ← b + E_d[h][k] (mod m).
   Production uses one blst_fp12_mul (one reference unit).

Correctness of step 4: σ^k(x') = σ^k(x)·R_h, so the class of x' is the class of the
canonical element times the adder, as required. Step 3 is exact.

## 6. Cycles, trails, distinguished points and solving

- Trails (v9). A walk's trail begins at its start element and ends at its first distinguished point (DP). At the DP
  the walk stores its record, then restarts with an adder chosen by a per-walk counter: let t be the walk's restart
  counter, incremented before use (so t = 1 at its first restart), and w the walk's global index (0, 1, ...). The
  restart index is h = mix((w + 1) * GAMMA + t) mod r_a, with GAMMA = 0x9E3779B97F4A7C15 and
  mix(x) = (z3) where z = x; z ^= z >> 30; z *= 0xBF58476D1CE4E5B9; z ^= z >> 27; z *= 0x94D049BB133111EB; z ^= z >> 31
  (all arithmetic mod 2^64). The restart is the additive step with index h, i.e. x <- x * T[h][k(x)],
  (a, b) <- (a + E_c[h][k], b + E_d[h][k]) mod r, so the invariant x = g^(a + b y) holds. It costs one multiplication
  (charged as an adder step) and clears the walk's detector state (previous key, anchor, trail age). The restarted state
  is itself a DP candidate: if it is a DP, its record is stored and the walk restarts again with the next counter value
  (chained restarts; probability 2^-60 per restart in production). A restart is not a walk step: it does not advance the
  anchor schedule, which refreshes every 64 walk steps (rounds, one step per walk per round). A walk stores no
  start element after loading. The counter makes the restart depend on the walk's history, so a walk that returns to
  an earlier DP does not repeat its earlier trail (v8 used a fixed start element, which could).
- 2-cycle (checked every step, per walk): the new key equals the key two steps back and the following step returns
  to the older key. Escape: square the least-key state of the cycle (commutes with sigma). The escape is a group
  operation and is counted.
- Longer cycles up to length 64 (per walk): an anchor key is recorded every 64 of the walk's own steps. When the new
  key equals the anchor key, walk the cycle (at most 64 steps) and escape from its least-key state if the walk closes.
  A cycle of length L <= 64 entered after tail length t is detected within 127 steps of t (the walk may traverse it
  several times before the anchor is on the cycle; those traversals are revisit steps, at most 127 per closure).
- Trail cap (per walk): if a trail reaches CAP = 64 * 2^dp_bits = 2^66 steps without a DP, the walk restarts as at a
  DP but stores no record (charged the same).
- Distinguished point: (K >> adder_bits) & (2^dp_bits − 1) = 0 (production dp_bits = 60). Store (fingerprint, A, B),
  80 bytes: the fingerprint is key bits 70..197 of the Montgomery integer of K (128 bits; bits 10..69 are zero for
  every distinguished point), and A, B are the canonical exponents A = lambda^kn a, B = lambda^kn b mod r, where kn is
  the argmin index. The record is stored before the restart.
- Duplicates. If a walk's trail reaches a state of another walk's trail or of its own earlier trail (a join), it follows that trail deterministically
  to that trail's DP, stores one duplicate of that record, and restarts with its next counter value (so the two trails
  diverge). A walk joins at most once per trail. A join is a collision: the two records have
  equal fingerprints and, generically, different B, so the solver finds y from the pair.
- Solving: records with equal fingerprints are compared pairwise. A pair from the same class with B1 != B2 gives
  y = (A2 - A1)/(B1 - B2) mod r; a degenerate pair (B1 = B2) is discarded and costs one comparison. Every candidate is
  verified by g^y = z before output. No trail is replayed.
- Walks: WALKS = 2^20 walks run in parallel. Each walk is a separate state with its own start element, detector and
  trail age. M is the number of main steps (of all walks together) that the charge bounds cover (section 9).

## 7. Parameters

- Production: r_a = 1024, adder_bits = 10, S = 458 square slots (s = S / r_a = 0.447265625).
  S is chosen by minimising the per-step charge over sqrt(V) (adder and square charges of section 8); the optimum
  is within a few parts in 10^6 of S = 458, which is the measured count.
- Toy: adder_bits and s are inputs of the experiment. The success statistic uses the partition
  V = 1 − s^2 − (r_a − S)/r_a^2 with S = s r_a (the r_a − S additive adders are each selected with
  probability 1/r_a). The v2 value, 1 − s^2 − (1 − s)^2/r_a, is withdrawn.

## 8. Measured quantities required by the accounting

- Per-type step costs: the adder step (reference operation blst_fp12_mul, charged 1, plus the other instructions
  divided by C = 20,000) and the square step (blst_fp12_cyclotomic_sqr at kappa, plus the other instructions divided by C).
  Every main step is charged the adder cost (the larger of the two; section 9.3). The step code is straight-line
  (comparisons over all limbs, branch-free modular addition and minimum), so each type has one instruction count,
  measured with every step of one type (forced inputs). The forced counts of revision v8 are in
  measurement/RESULT.json of the package (they supersede the v7 values after the restart change).
- Restart: one multiplication by the start element and the key derivation, charged as an adder step (A_UNITS).
- Distinguished-point costs: the store of one record (the dp path of the kernel, measured as the difference between the
  every-step-DP and the no-DP runs minus the restart; RESULT.json) and the bitonic sorting network of section 9.2:
  c_cmp instructions per compare-exchange and c_scan per slot (measured in the standalone kernel sortnet.c;
  measurement/sortnet/RESULT-sortnet.json).
- kappa: instructions of one blst_fp12_cyclotomic_sqr, divided by C.
- The operation counts per main step and the cycle-event counts come from uninstrumented runs of the same kernel
  (measurement/native_run.sh) over 2^20 rounds per seed (evidence/native-counts.json).
- Measurement: measurement/reproduce.sh (measure/measure.sh with LIB=blst, N1 = 256, N2 = 512, two seeds), and
  replica.py for the dump, distinguished-point, restart and count checks.

## 9. Procedure, work cutoff, capacity and the distinguished-point table

### 9.1 Work counter and cutoff

The procedure keeps a work counter W in units (charges as in section 8). Each operation adds its charge: main steps
(adder or square, each at the adder charge), restarts (one multiplication each, at DP and at the trail cap), detected-cycle
work (at most 66 extra operations per detection: walk-through, check and escape), distinguished points (section 9.2 and
9.4), walk starts (2 exponentiations each), verified fingerprint comparisons (one exponentiation each), the adder table,
and the final solve (2^29 instructions plus one pairing).

The walks run until W reaches W_MAX − reserve, where reserve is the solve and output charge. Then the procedure stops and
reports failure if no verified y was found. W_MAX is the claimed computation.

### 9.2 Distinguished-point table (buckets, v9)

- Buckets. The table is B = 2^37 independent buckets, each of n = 2^31 slots (64-bit indices), so each bucket fits a 64-bit
  address space. A record goes to bucket b = fp[1] >> 27 (the top 37 bits of the 128-bit fingerprint). Equal fingerprints
  therefore share a bucket, and matching needs no global sort. Total slots 2^68 (padding included).
- Store. A distinguished point appends its 80-byte record to its bucket (a cursor per bucket). If a bucket is full the
  procedure stops and reports failure (overflow rule); if the table holds D records before a store the same rule applies.
- Pad. Before sorting, the unused slots of each bucket are filled with sentinels (fingerprint all ones).
- Sort. Each bucket is sorted independently by a bitonic network on its n = 2^31 slots, ascending by the 128-bit fingerprint
  (fp[1] the high word, unsigned): for k = 2, 4, ..., n; for j = k/2, ..., 1; for i = 0..n−1 with l = i XOR j and l > i:
  compare-exchange(i, l), ascending if (i AND k) = 0. The compare-exchange is branch-free (mask compare, masked XOR swap of the
  80-byte records through one temporary). Each bucket has n k (k+1)/4 compare-exchanges, k = 31.
- Match. Equal neighbouring fingerprints in a bucket are enumerated over the equal run, each pair charged as a fingerprint
  match (section 9.3). A pair with B1 = B2 is discarded at one comparison.
- Scratch and memory. A bucket sorts in place; its scratch is one 80-byte temporary. Total storage 2^68 slots x 80 bytes
  (2^74.3 bytes), spread over the buckets; each bucket is 2^31 x 80 bytes, below the 64-bit address space.
- Cost. Per record: store (measured), one restart (charged as an adder step), bucket scatter (measured, scatter64),
  transfer and index widening (an explicit assumption, 100,000 instructions per record, far above a local copy), and a
  10-instruction overflow check. Per slot: the bucket sort (c_cmp64 per compare-exchange times 248 per slot, k = 31) and the
  scan (c_scan64 per slot). The table is charged for its full slots (padding included) and its records at capacity D.

### 9.3 Charges and reserves for the first M main steps (accounting.py, build())

Premise H7 (v9). Every new state's key is uniform given the past (the adder index and the square decision of the next step
are uniform); a restart index is a deterministic function of the walk's counter (section 6), treated as uniform. Escapes
are deterministic functions of the history and are counted as part of their closure event.

Unconditional means. Let p = event_rate_bound (section 9.5 and accounting.py) be the unconditional probability per main
step that a closure event starts (its first revisit of a window state of its own trail, lengths 2 to 64). Its expectation
over the first M steps is at most p M (a union over the windows, each closing with the probability of its formal or non-formal
pattern given fresh keys, by the tower property; the real process and the i.i.d.-key process agree until the first revisit).

- Closure events (Markov). The number of closure events E satisfies P(E > K p M) <= 1/K. The reserve is K = 50 (failure
  1/50 = 0.02, in delta). Each event is charged 66 adder-type operations (walk-through, check, escape).
- Revisit steps. Each closure event is followed by at most 127 revisit steps before its detection (the anchor fires within
  127 steps of entry, the walk being periodic on classes after the first revisit). These steps are main steps. So the
  revisit steps are at most 127 E <= 127 K p M, and the reserve fraction is w = 127 K p (+ the join suffix term below).
- Joins (Markov). A join is a first landing of a walk on a state of another trail (its own earlier trail included). Its
  expectation is at most M^2 / N (H1: a step's next state is a visited state with probability at most (visited)/N). With
  the allowance J_MAX = 1000, P(J > J_MAX) <= E[J] / J_MAX (failure about 2.8e-3 at the production M). A join is followed by at
  most CAP revisit steps (the joined trail's remainder), so the join suffixes are at most J_MAX CAP steps, a fraction of M of
  order 2^-50.
- Detection delay. A collision is detected when both walks reach a DP. Each walk's residual to its DP after the collision is
  geometric with mean 2^60 (under H7, DP flags of fresh keys have probability 2^-60). The allowance is T_DELAY = 40 mean DP
  waits per walk, so lag = WALKS * 40 * 2^60 steps, and the residual exceeds the allowance with probability at most
  e^-40 (in delta, once for the colliding pair). The collision rate is multiplied by (1 - e^-40) in the success model.
- Square fraction (Azuma). Square decisions of new states are conditionally Bernoulli(s) (H7), so the square count over the
  first M steps is within 10^-6 M of s M with probability at least 1 - exp(-2 10^-12 M). The collision rate V is evaluated at s
  (the fresh decisions); the revisit steps do not change it, they only reduce the fresh steps (reserve above).
- Fresh steps. The fresh steps are M minus the revisit steps; M is chosen so that M (1 - w) − J_MAX CAP >= T + lag, where w is the
  reserve fraction (127 K p) and J_MAX CAP covers the join suffixes. The success model needs T fresh steps.
- Distinguished points: at capacity (9.4). Restarts: one multiplication each (charged A_UNITS); cap resets are at most M/CAP.
- Fingerprint matches: metered (one exponentiation each); the metered total is bounded by 2^60 times the expected count
  (Markov, failure 2^-60).
- Non-formal closures of length at most 64: expectation 12 * 64 M / r (Markov).
- Trap events (undetected cycles longer than 64): the trap term of section 9.5.

### 9.4 Distinguished-point capacity

Records are fresh DP flags X (each new fresh key flags with probability 2^-60; X has mean mu = M 2^-60, and P(X > 2 mu) <= exp(-mu/3)
by Chernoff under H7) plus duplicates, one per join, so records <= X + J pathwise. The capacity is D = ceil(2 mu) + J_MAX.
Overflow needs X > 2 mu (failure exp(-mu/3)), or J > J_MAX (failure E[J]/J_MAX, 9.3), or a bucket above 2^31 records. Bucket
occupancy is binomial with mean D/2^37 (about 2^30.4 at the production M) against capacity 2^31; the Chernoff bound per bucket
and the union over 2^37 buckets give a negligible failure (accounting.py, bucket_failure_log2). The procedure stops when a
bucket or the total capacity is exceeded (9.2), so only these events bind.

### 9.5 Formal cancellations and window lengths

A closed window is a window of steps whose composition returns the walk to an earlier class or state. Its exponent coefficient is
the sum of the step contributions. If the window contains a square, the coefficient carries a factor 2^q with q >= 1 that is not
a 12th root of unity (section 9.6), so it is nonzero: such closures are non-formal and bounded by the 1/r term (H5). If the
window is all additive, the coefficient of adder h is sum_k n_{h,k} lambda^(-k) in F_r, where n_{h,k} counts the steps of adder h
with rotation k. Lemma (norm). For windows of length L < 2^63, this is zero in F_r only if sum_k n_{h,k} zeta^(-k) = 0 in
Z[zeta12]. Proof: the image of alpha = sum n_{h,k} zeta^(-k) in F_r is zero, so alpha lies in a prime 𝔭 of Z[zeta12] above r
(r = 1 mod 12, so r splits), and a nonzero alpha in 𝔭 has norm divisible by r with |N(alpha)| <= L^4 < 2^252 < r, so alpha = 0.

Proposition G (vanishing sums of 12th roots of unity). Let n in Z_{>=0}^12 with sum_k n_k zeta^k = 0 for a primitive 12th root of
unity zeta. Then n is a nonnegative integer combination of the indicator vectors of the 2-gons {k, k+6} and the 3-gons {k, k+4, k+8}
(indices mod 12). The proof is in the package's proof.md, section 2 (minimal polynomial x^4 − x^2 + 1; four coordinate
conditions C1-C4; signed antipodal differences; induction on the weight), and is stated here without repetition.

Per adder, the closed window's rotation multiset is a nonnegative solution, so Proposition G gives its decomposition. The union bound
of accounting.py (union_log2, lengths 65 to 400) counts the gon decompositions of each length. For windows longer than 400 the
probability that every step is additive, (1 − s)^401/s, is used.

### 9.6 Order of 2 and the smallest q

r − 1 = 2^32 · 3 · 11 · 19 · 10177 · 125527 · 859267 · 906349^2 · 2508409 · 2529403 · 52437899 · 254760293^2 (cscore/number_theory.py:
the product, primality of each factor by Miller–Rabin with the first twelve primes as bases, deterministic below 3.18·10^23, and r
as a probable prime with the same bases). The multiplicative order of 2 modulo r is o = 3277242198445386904965483781761622864855659531282977363912728668746161324032,
251 bits, certified by 2^o = 1 and 2^(o/p) != 1 for each prime divisor p of o. Since 2^q is a 12th root of unity iff o | 12q iff q0 | q,
with q0 = o/gcd(o,12) (gcd(q0, 12/gcd(o,12)) = 1), the smallest q >= 1 is q0, a 248-bit number larger than 2^247. Hence 2^q is not a
12th root of unity for 1 <= q < 2^247.

### 9.7 Failure terms

Each term is computed in accounting.py (build()) and summed into the failure total subtracted from the success probability:
the trap term (9.5), the fingerprint Markov term 2^-60, the square-fraction Azuma term (9.3), the closure-event Markov term 1/K
(9.3), the non-formal closure expectation, the join Markov term (9.3), the fresh-flag overflow exp(-mu/3) and bucket overflow
(9.4), and the detection delay e^-40 (9.3). Success is then the model with the collision factor (1 - e^-40) minus the total.
