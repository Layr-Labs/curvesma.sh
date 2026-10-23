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

## 6. Cycles, distinguished points and solving

- 2-cycle: the new state has the same key as the previous state, and the next step returns
  to the old state. Escape: square the least-key state of the cycle (commutes with σ).
- Short cycles (≤ 64): an anchor key is recorded every 64 global steps. On a match, walk the
  cycle and escape from its least-key state, as in cscore/rho.py.
- Distinguished point: (K >> adder_bits) & (2^dp_bits − 1) = 0. Store (K, A, B) from the
  canonical exponents. A repeated DP key gives a candidate y; verify g^y = z.
- The walk is parallel over M walks (production M = 2^20 with theta = 2^-60) in the
  accounting. Each walk is a separate state.

## 7. Parameters

- Production: r_a = 1024, adder_bits = 10, s chosen from the measured c2 and κ (RESEARCH.md
  section 3.4; s* ≈ 0.46 at c2 = 0.12, κ = 0.49).
- Toy: adder_bits and s are inputs of the experiment. The success statistic uses
  V = 1 − s² − (1 − s)²/r_a (RESEARCH.md section 3.4).

## 8. Measured quantities required by the accounting

- c2: instructions per step of work other than the reference multiplication (key,
  argmin, table lookup, exponent updates, DP test, cycle filter), divided by C = 20,000.
- κ: instructions of one blst_fp12_cyclotomic_sqr (and its wrapper), divided by C.
- Both are per step averages over the walk mixture, measured by Callgrind on two input
  sets (measure/measure.sh, LIB=blst, REF_FUNCS=ref_mul).
