# GT Frobenius-class walk: shared specification (revision v7)

Used by `curvesmash-core/cscore/gt.py` (toy model, Python), the blst step kernel
(`curvesmash-core/examples/blst-gt-rho/step.c`) and its Python replica. Written
2026-10-22; revision v7 2026-10-23 (sections 6, 8, 9; the distinguished-point table and the
classification are now specified in full). Any change here must be applied to all three.

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

- 2-cycle (checked every step, per walk): the new key equals the key two steps back and the
  following step returns to the older key. Escape: square the least-key state of the cycle
  (commutes with sigma). The escape is a group operation and is counted.
- Longer cycles up to length 64 (per walk): an anchor key is recorded every 64 of the walk's own
  steps. When the new key equals the anchor key, walk the cycle (at most 64 steps) and escape from
  its least-key state if the walk closes. A cycle of length L <= 64 entered after tail length t is
  found within 127 steps of t.
- Undetected cycles and the trail cap (per walk): each walk keeps a trail age, the main steps since
  its last distinguished point. When the age reaches CAP = 64 * 2^dp_bits = 2^66 steps, the walk is
  restarted from a fresh random start (2 exponentiations, charged), and its age is reset. A detected
  cycle wastes at most 127 main steps, charged in the budget.
- Distinguished point: (K >> adder_bits) & (2^dp_bits − 1) = 0 (production dp_bits = 60). Store
  (fingerprint, A, B), 80 bytes: the fingerprint is key bits 70..197 of the Montgomery integer of K
  (128 bits; bits 10..69 are zero for every distinguished point), and A, B are the canonical exponents
  A = λ^kn a, B = λ^kn b mod r, where kn is the argmin index.
- Solving: records with equal fingerprints are compared pairwise. A pair from the same class with
  B1 ≠ B2 gives y = (A2 − A1)/(B1 − B2) mod r; a degenerate pair (B1 = B2) is discarded and costs one
  comparison. Every candidate is verified by g^y = z before output. No trail is replayed.
- Walks: WALKS = 2^20 walks run in parallel. Each walk is a separate state with its own anchor and age.
  M is the number of main steps (of all walks together) that the charge bounds cover (section 9).

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
  measured with every step of one type (forced inputs): 26,730.2 (adder, with the anchor refresh charged) and
  13,511.4 (square) instructions.
- Distinguished-point costs: the store of one record, 1,479 instructions (measured, the dp path of the kernel); the
  bitonic sorting network of section 9.2, c_cmp instructions per compare-exchange and c_scan instructions per slot
  (measured with modes `sortnet` and `sortscan` of step.c, measurement/sortnet/RESULT-sortnet.json).
- kappa: instructions of one blst_fp12_cyclotomic_sqr, divided by C.
- The operation counts per main step (adds, squarings, both including escapes and cycle-check
  steps) and the cycle-event counts come from uninstrumented runs of the same kernel
  (measurement/native_run.sh) over 2^26 main steps per seed.
- Measurement: measurement/reproduce.sh (measure/measure.sh with LIB=blst, N1 = 256, N2 = 512, two
  seeds), and replica.py for the dump, distinguished-point and count checks.

## 9. Procedure, work cutoff, capacity and the distinguished-point table

### 9.1 Work counter and cutoff

The procedure keeps a work counter W in units (charges as in section 8). Each operation adds its charge:
main steps (adder or square, each at the adder charge), detected-cycle work (at most 66 adder-type operations per
detection: walk-through, check and escape), distinguished points (section 9.2), walk starts and restarts
(2 exponentiations each), verified fingerprint comparisons (one exponentiation each), the adder table, and the
final solve (2^29 instructions plus one pairing).

The walks run until W reaches W_MAX − reserve, where reserve is the solve and output charge. Then the
procedure stops and reports failure if no verified y was found. W_MAX is the claimed computation.

### 9.2 Distinguished-point table (concrete procedure)

- Store. A distinguished point appends its 80-byte record (fingerprint, A, B) to the table array of n slots, where
  n is the smallest power of two not less than the capacity D (section 9.4). Slot c is written and c is incremented.
  If the table holds D records before a store, the procedure stops and reports failure (overflow rule).
- Pad. After the walks stop, the slots c..n−1 are filled with sentinel records whose fingerprint is all ones, so they
  sort last. Padding is part of the table's memory (section 9.4).
- Sort. A bitonic sorting network on the n slots, ascending by the 128-bit fingerprint (fp[1] the high word, fp[0] the
  low word, unsigned). Standard iterative form: for k = 2, 4, ..., n; for j = k/2, ..., 1; for i = 0..n−1 with
  l = i XOR j and l > i: compare-exchange(i, l), ascending if (i AND k) = 0, descending otherwise. Each compare-exchange
  computes the comparison as an all-ones or all-zeros mask from the two fingerprints and swaps the two 80-byte records
  with masked XOR through one 80-byte temporary (no branch on data). The network has n k (k + 1) / 4 compare-exchanges,
  k = log2 n. Its control and index arithmetic depend only on n and are fixed.
- Match. For i = 0..c−2, a fingerprint equal to the next slot's (computed without a data-dependent branch, a constant-cost
  pass) identifies a pair of equal fingerprints, which is enumerated over the equal run and charged as a fingerprint
  match (section 9.3). Records with equal fingerprints but B1 = B2 are discarded at one comparison.
- Scratch memory. The network sorts in place. Scratch is the 80-byte temporary and the padding slots (n − c records),
  both inside the n-slot table counted in section 9.4. No other scratch is used.
- Cost. Instructions: store 1,479 per record (measured), c_cmp per compare-exchange (measured), c_scan per slot
  (measured, sortscan minus sortnet), and 10 per record for the overflow test. The table is charged for its capacity D:
  DP_inst(D) = 1,479 D + c_cmp n k(k+1)/4 + c_scan n + 10 D, with n the power of two of section 9.2. Charging the
  capacity, not the number of stored records, keeps the charge independent of the random count.

### 9.3 Charges for the first M main steps (accounting.py, build())

Premise H7 (conditional uniformity). Given the history F up to step i − 1, the adder index h = K(x) mod r_a of step i is
uniform on 0..r_a − 1, the distinguished-point and fingerprint bits of a fresh state are uniform, and the state of an
escape or an anchor is a deterministic function of F. A revisit step repeats the state of an earlier step and thus
repeats its decision bits.

- Operation mix (no concentration needed for the charge). Every main step is charged A_UNITS, the larger of the two step
  charges, so the charge does not depend on the square count. The square count enters only the success model through
  the realised fraction of square decisions. Over the fresh steps the decisions are conditionally Bernoulli(s) (H7), so
  Hoeffding gives deviation 10^-6 with failure exp(−2 · 10^-12 M). Revisit steps repeat an earlier decision, and there
  are at most as many revisit steps as detected cycles (section below). The realised square fraction is therefore within
  ε_s = 10^-6 + p_ev · (1 + 10^-3) of s, and the success model uses the collision rate V_low = min V(s ± ε_s).
- Detected cycles. A detection at step i needs the window ending at i (of length ℓ ≤ 64, from an anchor or a 2-cycle) to
  close. Its probability given the history F_{i−64} is at most p_ev = (2-cycle term) + (formal union over lengths 3..64,
  section 9.5) + (non-formal 12 · 64 / r), times 1 + 10^-3. Splitting the steps into 64 residue classes mod 64, the
  indicators along each class are dominated by independent Bernoulli(p_ev) variables (each has conditional probability at
  most p_ev given the history 64 steps earlier). Chernoff on each class and a union over the 64 classes give failure at most
  64 exp(−(10^-3)^2 (p_ev M / 64) / 3). Each detection is charged 66 adder-type operations: 66 A_UNITS per detection, so the
  charge is 66 A_UNITS p_ev per step on average, bounded as a count.
- Non-formal detections of length at most 64 have expected count 12 · 64 M / r (Markov).
- Distinguished points: section 9.4.
- Fingerprint matches are metered when they occur (one exponentiation each). The metered total is bounded by 2^60 times the
  expected count (Markov, failure 2^-60); above that bound the procedure stops with failure.
- Restarts: a walk restarts at age 2^66 (charged 2 exponentiations); the restarts among the first M steps number at most
  M/2^66 + WALKS.

### 9.4 Distinguished-point capacity

Let μ = M 2^-60 and c = (WINDOW + 1)(1 − s)/r_a = 0.0691 (WINDOW = 127). The capacity is

    D = ⌈2 (μ + 64) / (1 − c)⌉ + 2 WALKS,

and n = the smallest power of two ≥ D (padding included in memory).

Justification. Let X count the distinguished-point flags of fresh states. Under H7 each is an independent-looking test with
probability 2^-60 given the history, so X is dominated by Binomial(M, 2^-60) with mean μ. Let Y count the flags of revisits.
A revisit lands on a state of the last WINDOW + 1 steps (the trail window and the anchor), and an adder step reaches a given
earlier state with probability at most (1 − s)/r_a (the adder index is uniform and at most one index reaches the state).
So a flagged state is revisited with expected count at most c per step during which it is visible, and E[Y | past] ≤ c D_tot.
Cross-walk collisions (a step landing on a state of another walk) flag at most 64 records with failure at most
(e μ_c / 64)^64, μ_c = M²/(2N) ≈ 0.8. Hence D_tot ≤ (X + 64)/(1 − c). The overflow failure is at most exp(−μ/3) (Chernoff on X)
plus the collision term. The 2 WALKS term is an allowance for the detection lag of each walk (the walk's last record before
its collision is seen).

The overflow rule (section 9.2) stops the procedure when the table is full, so the only event that binds is D_tot > D.

### 9.5 Formal cancellations and window lengths

A closed window is a window of steps whose composition returns the walk to an earlier class or state. Its exponent coefficient
is the sum of the step contributions. If the window contains a square, the coefficient carries a factor 2^q with q ≥ 1 that is
not a 12th root of unity (section 9.7), so it is nonzero: such closures are non-formal and bounded by the 1/r term (H5). If the
window is all additive, the coefficient of adder h is Σ_k n_{h,k} λ^{-k} in F_r, where n_{h,k} counts the steps of adder h with
rotation k. Lemma (norm). For windows of length L < 2^63, this is zero in F_r only if Σ_k n_{h,k} ζ^{-k} = 0 in Z[ζ12]. Proof:
the image of α = Σ n_{h,k} ζ^{-k} in F_r is zero, so α lies in a prime 𝔭 of Z[ζ12] above r (r ≡ 1 mod 12 splits). A nonzero α in
𝔭 has norm divisible by r, and |N(α)| ≤ L^4 < 2^252 < r, so α = 0.

Proposition G (vanishing sums of 12th roots of unity). Let n ∈ Z_{≥0}^{12} with Σ_k n_k ζ^k = 0 for a primitive 12th root of
unity ζ. Then n is a nonnegative integer combination of the indicator vectors of the 2-gons {k, k+6} and the 3-gons {k, k+4, k+8}
(indices mod 12).

Proof. (i) Minimal polynomial. Φ(x) = x^4 − x^2 + 1 is the minimal polynomial of ζ over Q. It vanishes at ζ, since ζ^6 = −1 gives
0 = ζ^6 + 1 = (ζ^2 + 1)(ζ^4 − ζ^2 + 1) with ζ^2 ≠ −1. It is irreducible: a factorisation over Z (Gauss's lemma, Φ monic) has no
linear factor since Φ(±1) = 1; a quadratic factorisation (x^2 + ax + b)(x^2 − ax + c) gives b + c − a^2 = −1, a(c − b) = 0, bc = 1,
so b = c = 1 (then a^2 = 3) or b = c = −1 (then a^2 = −1), both impossible.
Hence Σ n_k ζ^k = 0 iff the remainder of f(x) = Σ n_k x^k modulo Φ is zero.

(ii) Coordinates. Modulo Φ: x^4 ≡ x^2 − 1, x^5 ≡ x^3 − x, x^6 ≡ −1, and x^{k+6} ≡ −x^k. The coordinates of x^k (k = 0..11), in the basis
1, x, x^2, x^3, are v_0 = (1,0,0,0), v_1 = (0,1,0,0), v_2 = (0,0,1,0), v_3 = (0,0,0,1), v_4 = (−1,0,1,0), v_5 = (0,−1,0,1),
v_6 = (−1,0,0,0), v_7 = (0,−1,0,0), v_8 = (0,0,−1,0), v_9 = (0,0,0,−1), v_10 = (1,0,−1,0), v_11 = (0,1,0,−1). The remainder is zero iff

  (C1) n_0 − n_4 − n_6 + n_10 = 0,  (C2) n_1 − n_5 − n_7 + n_11 = 0,
  (C3) n_2 + n_4 − n_8 − n_10 = 0,  (C4) n_3 + n_5 − n_9 − n_11 = 0.

(iii) Signed antipodal differences. Let p_j = n_j − n_{j+6} for j = 0..5. Then (C1)–(C4) read p_0 = p_4, p_1 = p_5, p_2 = −p_4 and p_3 = −p_5.
With a = p_4 and b = p_5 this gives p = (a, b, −a, −b, a, b).

(iv) Induction on the total weight. The zero vector is the empty combination. Let n ≠ 0 be a nonnegative solution. If some pair has
n_j > 0 and n_{j+6} > 0, the 2-gon {j, j+6} lies in the support; subtracting its indicator keeps n nonnegative and a solution (the gon
vector is a solution, since v_j + v_{j+6} = 0). Otherwise n_j n_{j+6} = 0 for all j, so p_j = n_j or p_j = −n_{j+6}, and n ≠ 0 forces
(a, b) ≠ (0, 0). If a > 0, then n_0 = a, n_4 = a and n_8 = a (from p_0 = a, p_4 = a, p_2 = −a with the other entry of each pair zero),
so the 3-gon {0, 4, 8} lies in the support. If a < 0, the 3-gon {2, 6, 10} lies in the support (n_2 = n_6 = n_10 = −a). If a = 0 and b > 0,
the 3-gon {1, 5, 9} lies in the support; if a = 0 and b < 0, the 3-gon {3, 7, 11} does. The 3-gon indicator is a solution, since
ζ^r(1 + ζ^4 + ζ^8) = 0 (ζ^4 is a primitive cube root of unity). Subtracting it keeps n nonnegative and a solution. Induction on the weight
completes the proof. ∎

Per adder, the closed window's rotation multiset is a nonnegative solution, so Proposition G gives its gon decomposition. The union bound of
accounting.py (union_log2, lengths 65..400) counts the gon decompositions of each length; the classification needs no enumeration. For
windows longer than 400 the all-additive probability (1 − s)^401/s is used, and no classification is needed.

### 9.6 Order of 2 and the smallest q

r − 1 = 2^32 · 3 · 11 · 19 · 10177 · 125527 · 859267 · 906349^2 · 2508409 · 2529403 · 52437899 · 254760293^2 (number_theory.py: the
product, primality of each factor by Miller–Rabin with the first twelve primes as bases, deterministic below 3.18·10^23, and r itself as a
probable prime with the same bases). The multiplicative order of 2 modulo r is o = 3277242198445386904965483781761622864855659531282977363912728668746161324032,
251 bits, certified by 2^o = 1 and 2^(o/p) ≠ 1 for each of the twelve primes p dividing o. The smallest q ≥ 1 with 2^q a 12th root of unity is
q0 = o/gcd(o, 12) (2^q is a 12th root of unity iff o | 12q iff q0 | q, using gcd(q0, 12/gcd(o,12)) = 1), a 248-bit number larger than 2^247.
Consequently 2^q is not a 12th root of unity for 1 ≤ q < 2^247. The test 2^((r−1)/p) for primes p | r − 1 is reported in number_theory.py
for information only (it shows that 2 is not a primitive root, since 2^((r−1)/2) = 1); the order certificate uses the exact order.

### 9.7 Failure terms

Each failure term is computed in accounting.py (build()): the trap term (formal union over 65..400, additive tail, and the non-formal
closures under H5), the fingerprint Markov term 2^-60, the square-fraction Hoeffding term (section 9.3), the detection concentration
(section 9.3), the non-formal detection expectation, the distinguished-point overflow (section 9.4) and the cross-walk collision allowance
(section 9.4). Their sum is the failure total subtracted from the success probability.
