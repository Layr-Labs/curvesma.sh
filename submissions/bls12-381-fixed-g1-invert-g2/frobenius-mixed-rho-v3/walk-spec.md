# GT Frobenius-class walk: shared specification

Used by `curvesmash-core/cscore/gt.py` (toy model, Python), the blst step kernel
(`curvesmash-core/examples/blst-gt-rho/step.c`) and its Python replica. Written
2026-10-22 for the rho package (RESEARCH.md section 7, item 1). Any change here must be
applied to all three.

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
  cycle wastes at most 127 main steps, charged in the budget. A trapped walk therefore loses at most
  CAP steps before the restart.
- Distinguished point: (K >> adder_bits) & (2^dp_bits - 1) = 0 (production dp_bits = 60). Store
  (fingerprint, A, B): the fingerprint is key bits 70..197 of the Montgomery integer of K (128 bits;
  bits 10..69 are zero for every distinguished point, so they are not used), and A, B are the
  canonical exponents A = lambda^kn a, B = lambda^kn b mod r, where kn is the argmin index (sigma^kn(x)
  is the canonical element, so its exponent is lambda^kn (a + b y)). Each record is 80 bytes.
- Solving: records with equal fingerprints are all compared pairwise. A pair from the same class with
  B1 != B2 gives y = (A2 - A1)/(B1 - B2) mod r; a degenerate pair (B1 = B2) is discarded and costs one
  comparison. Every candidate is verified by g^y = z before output. No trail is replayed.
- The walk is parallel over M = 2^20 walks. Each walk is a separate state with its own anchor and age.

## 7. Parameters

- Production: r_a = 1024, adder_bits = 10, S = 458 square slots (s = S / r_a = 0.447265625).
  S is chosen by minimising the per-step charge over sqrt(V) (adder and square charges of section 8); the optimum
  is within a few parts in 10^6 of S = 458, which is the measured count.
- Toy: adder_bits and s are inputs of the experiment. The success statistic uses the partition
  V = 1 - s^2 - (r_a - S)/r_a^2 with S = s r_a (the r_a - S additive adders are each selected with
  probability 1/r_a). The v2 value, 1 - s^2 - (1 - s)^2/r_a, is withdrawn.

## 8. Measured quantities required by the accounting

- Per-type step costs: the adder step (reference operation blst_fp12_mul, charged 1, plus the other instructions
  divided by C = 20,000) and the square step (blst_fp12_cyclotomic_sqr at kappa, plus the other instructions divided by C).
  The step code is straight-line (comparisons over all limbs, branch-free modular addition and minimum), so each type has one
  instruction count, measured with every step of one type (forced inputs), at 26,730.2 (adder, with the anchor refresh charged)
  and 13,511.4 (square) instructions. The distinguished-point store is measured (1,479 instructions) and the sort per
  record is 67 levels of a measured 48-instruction exchange.
- kappa: instructions of one blst_fp12_cyclotomic_sqr, divided by C.
- The operation counts per main step (adds, squarings, both including escapes and cycle-check
  steps) and the cycle-event counts come from uninstrumented runs of the same kernel
  (measurement/native_run.sh) over 2^26 main steps per seed.
- Measurement: measurement/reproduce.sh (measure/measure.sh with LIB=blst, N1 = 256, N2 = 512, two
  seeds), and replica.py for the dump, distinguished-point and count checks.

## 9. Procedure, hard work cutoff and capacity

- The procedure keeps a work counter W in units (charges as in section 8). Each operation adds its charge:
  main steps (2 kinds: adder and square), detected-cycle work (at most 66 extra operations per detection:
  walk-through, check and escape), distinguished points (2^20 instructions each), walk starts and restarts
  (2 exponentiations each), verified fingerprint comparisons (one exponentiation each), the adder table, and the
  final solve (2^29 instructions plus one pairing).
- The walks run until W reaches W_MAX - reserve, where reserve is the solve and output charge. Then the
  procedure stops and reports failure if no verified y was found. W_MAX is the claimed computation.
- Charge bound for the first M main steps (accounting.py, build()): under H7 (the adder and square choices are i.i.d.
  uniform), the square count is at least (s - 1e-6) M with probability at least 1 - exp(-2e-12 M) (Hoeffding). Each
  detected cycle of length at most 64 is a function of at most 64 consecutive choices, so the indicators split into 64
  residue classes mod 64 whose members are independent; Chernoff on each class and a union bound over the 64 classes
  give failure at most 64 exp(-(1e-3)^2 (p M / 64) / 3), with p the per-step detection bound. Non-formal detections
  have expected count 12 * 64 M / r. Each detection is charged at most 66 adder-type steps.
  The distinguished-point count is bounded by the table capacity (Chernoff, failure exp(-M 2^-60 / 3)).
- Fingerprint matches are metered: each match is charged when it occurs (one exponentiation). The metered total is
  bounded by 2^60 times the expected count (Markov, failure 2^-60); above that bound the procedure stops with failure.
- Distinguished-point capacity: the table holds D = 2 M 2^-60 + 2 WALKS records of 80 bytes. When the table is
  full, the procedure stops and reports failure (overflow rule). The overflow probability is at most
  exp(-M 2^-60 / 3), which is below 2^-(10^19) at the production M.
- Trail cap: each walk's age (main steps since its last DP) is kept as a 128-bit counter. At 2^66 the walk
  restarts (charged as one restart per DP, which bounds the restarts by M / 2^66 + WALKS).
- Formal cancellations (classification used in the trap bound): a closed window with exact coefficient zero
  is a window of additive steps whose multiset of (adder, rotation) pairs decomposes into 2-gons (k, k+6)
  and 3-gons (k, k+4, k+8) per adder (Lam-Leung; checked by enumeration up to total count 7). The order of 2
  modulo r is computed exactly from the factorisation of r - 1, so 2^q is not a 12th root of unity for
  q < 2^247 and every window containing a square has a nonzero coefficient.
