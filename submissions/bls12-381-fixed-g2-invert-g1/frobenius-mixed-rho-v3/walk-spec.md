# GT Frobenius-class walk: shared specification (revision v8)

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

- Trails (v8). A walk's trail begins at its start element (x0, a0, b0) and ends at its first distinguished
  point (DP). At the DP the walk stores its record and restarts: x <- x * X0, (a, b) <- (a + A0, b + B0) mod r, where
  (X0, A0, B0) is the walk's own start element (fixed per walk; X0 = g^(A0 + B0 y)), so the invariant x = g^(a + b y)
  holds. The restart is one multiplication (charged as an adder step) and clears the walk's detector state (previous
  key, anchor, trail age).
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
- Duplicates. If a walk's trail reaches a state of another walk's trail (a join), it follows that trail deterministically
  to the other trail's DP, stores one duplicate of that record, and restarts at its own X0 (different from the other
  walk's, so the two trails diverge). A walk joins at most once per trail. A join is a collision: the two records have
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

### 9.2 Distinguished-point table (concrete procedure)

- Store. A distinguished point appends its 80-byte record (fingerprint, A, B) to the table array of n slots, where n is the
  smallest power of two not less than the capacity D (section 9.4). If the table holds D records before a store, the
  procedure stops and reports failure (overflow rule).
- Pad. After the walks stop, the slots c..n−1 are filled with sentinel records whose fingerprint is all ones, so they sort last.
- Sort. A bitonic sorting network on the n slots, ascending by the 128-bit fingerprint (fp[1] the high word, fp[0] the low
  word, unsigned). Standard iterative form: for k = 2, 4, ..., n; for j = k/2, ..., 1; for i = 0..n−1 with l = i XOR j and
  l > i: compare-exchange(i, l), ascending if (i AND k) = 0, descending otherwise. Each compare-exchange computes the comparison
  as an all-ones or all-zeros mask and swaps the two 80-byte records with masked XOR through one 80-byte temporary. The network
  has n k (k + 1)/4 compare-exchanges, k = log2 n; its control and index arithmetic depend only on n.
- Match. For i = 0..c−2, equal neighbouring fingerprints are enumerated over the equal run and charged as fingerprint matches
  (section 9.3). A pair with B1 = B2 is discarded at one comparison.
- Scratch memory. The network sorts in place: the scratch is the 80-byte temporary and the padding (n − c records at most),
  both inside the n-slot table.
- Cost. Instructions: store (measured) and restart (charged as an adder step) per record, c_cmp per compare-exchange, c_scan
  per slot, and 10 per record for the overflow test. The table is charged for its capacity D.

### 9.3 Charges for the first M main steps (accounting.py, build())

Premise H7 (conditional uniformity). Given the history F, the key of a fresh state (one that no walk has visited) is uniform,
so its adder index h = K mod r_a and rotation are uniform, and its distinguished-point flag has probability 2^-60. Escapes and
restarts are deterministic functions of the history (a restart uses the walk's start element). H7 is a stated premise.

- Main steps. Every main step (adder or square) is charged A_UNITS, the larger of the two step charges. No operation-mix
  concentration enters the charge.
- Closure events. An event E_i is a closure at step i: the walk revisits a state of its own trail inside its window of 64
  steps (2-cycle or anchor match, length 2 to 64). Claim: P(E_i | F_{i−64}) <= p_ev, where p_ev = (2-cycle term (1 − s)/(12 r_a))
  + (formal union over lengths 3 to 64, section 9.5) + (non-formal 12 * 64/r), times 1 + 10^-3. Reason: the window's states
  from its first state onward are fresh until the closure, so their keys are uniform given the history before the window
  (H7), and the closure probability over those uniform keys is the union bound p_ev. Splitting the steps into 64 residue
  classes mod 64, the events of a class are dominated by independent Bernoulli(p_ev) variables (conditional bound given the
  history up to the start of the class's previous window), so Chernoff on each class and a union bound give failure at most
  64 exp(−(10^-3)^2 (p_ev M / 64)/3). Each event is charged 66 adder-type operations (walk-through, check, escape).
- Revisit steps and the square fraction. A closure is followed by at most 127 revisit steps before its detection (the anchor
  fires within 127 steps of entry). The revisit steps repeat earlier square decisions, so the realised square fraction over
  the first M steps differs from s by at most 10^-6 (Hoeffding over the fresh decisions, which are Bernoulli(s) given the
  history, failure exp(−2 · 10^-12 M)) plus 127 p_ev (1 + 10^-3) plus the join suffixes (9.4 below, at most 2^-50 of M).
  The success model uses V_low = min V(s ± eps_s) with eps_s = 10^-6 + 127 p_ev (1 + 10^-3) + 2^-50.
- Joins. Each join adds a suffix of at most CAP revisit steps (the joined trail's remainder, bounded by the trail cap), so
  the join suffixes are at most JOIN * CAP steps where JOIN is the number of joins, bounded as in 9.4.
- Non-formal detections of length at most 64: expected count 12 * 64 M / r (Markov).
- Distinguished points: 9.4.
- Fingerprint matches are metered when they occur (one exponentiation each); the metered total is bounded by 2^60 times the
  expected count (Markov, failure 2^-60).
- Trap events (cycles longer than 64 that are not detected, and their revisit steps up to the trail cap): the trap term of
  section 9.5 bounds their expected number by M times the per-step trap probability, so they occur with probability at most
  the trap term.

### 9.4 Distinguished-point capacity (pathwise)

Each record is produced at a DP step. A DP step is either a fresh flag (its state is visited for the first time by any walk)
or a duplicate (the walk has joined another trail and reaches that trail's DP, which that trail has already recorded). A walk
joins at most once per trail (after the join it follows the other trail to its DP and restarts), so duplicates are at most the
joins J. Hence the number of records is at most X + J, pathwise, where X is the number of fresh DP flags.

- X is dominated by Binomial(M, 2^-60) (each fresh state flags with probability 2^-60 given the past, H7).
- J is dominated by Binomial(M, M/N): a step's next state is a visited state of another trail with probability at most
  (visited states)/N <= M/N (H1). Adder-level form: a given admissible additive index is chosen with probability 1/r_a (the index
  is uniform over the r_a values; the square decision is the same key bits, so the index is not multiplied by 1 - s), and with the
  adder elements independent of the path (H5) a given state is an adder neighbour of the current state with probability about
  r_a/N, which gives (visited)/N.

The capacity is D = ⌈2 μ⌉ + 64 with μ = M 2^-60, n = the smallest power of two not less than D. Overflow requires X > 2μ (Chernoff,
failure exp(−μ/3)) or J > 64 (failure at most (e μ_J / 64)^64 with μ_J = M^2/N). The procedure stops when the table is full
(section 9.2), so only these events bind. The restart count among the first M steps is at most D + M/CAP (charged).

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

Each failure term is computed in accounting.py (build()): the trap term (formal union over 65..400, additive tail, and the non-formal
closures under H5), the fingerprint Markov term 2^-60, the square-fraction Hoeffding term (section 9.3), the closure-count
concentration (section 9.3), the non-formal closure expectation, the DP fresh-flag overflow and the join allowance (section 9.4).
Their sum is subtracted from the success probability.
