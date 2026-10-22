# Fixed-G1 pairing inversion via Frobenius-orbit parallel rho in GT

Track: `bls12-381-fixed-g1-invert-g2`.

* Claimed score: **time_log2 = 126.32** (3158/25) target-group-operation
  equivalents. The unit is the dense Fp12 multiplication, at C = 20,000 scalar
  AArch64 instruction equivalents each.
* Success probability: at least 0.39, under a stated random-mapping heuristic
  supported by experiments.

This is generic accounting of the standard reduction combined with the standard
automorphism speedup for rho (Duursma–Gaudry–Morain 1999; Gallant–Lambert–Vanstone
2000; Wiener–Zuccherato 1998). It is not a new attack.

## 1. Reduction

Input: a nonidentity P in G1 (order r) and z in GT. Q0 is the public G2
generator from the profile.

1. Validate P and z as specified in the profile.
2. If z = 1, output Q = infinity (the only valid answer).
3. Compute g = e(P, Q0) with the pinned map (`pairing(Q0, P)` in py_ecc
   argument order). The pinned map is bilinear and nondegenerate on G1 x G2,
   and P and Q0 have prime order r, so g != 1 generates GT (prime order r).
4. Find a in [0, r) with z = g^a (section 2).
5. Output Q = [a]Q0. Then e(P, Q) = e(P, Q0)^a = g^a = z, and Q is in G2.
   Optionally recompute e(P, Q) and compare all 12 coefficients.

P is never replaced. This is the admissible route described in `rules.md`.

## 2. Discrete logarithm in GT on Frobenius orbits

### 2.1 Frobenius structure

The tower is the blst / py_ecc one:

* Fp2 = Fp[u]/(u^2 + 1)
* Fp6 = Fp2[v]/(v^3 - xi), with xi = 1 + u
* Fp12 = Fp6[w]/(w^2 - v), so w^6 = xi

Write x = (a0 + a1 v + a2 v^2) + (b0 + b1 v + b2 v^2) w, with ai, bi in Fp2.
The profile's flat w-basis is converted to this tower at input and output only.

Let phi(x) = x^p. phi conjugates Fp2 coefficients and maps w to
w · xi^((p-1)/6). Hence the w-coefficient of phi^k(x) is

    conj^k(b0) · delta_k,   delta_k = xi^((p^k - 1)/6) in Fp2   (public constants).

We checked these constants numerically for this p:

* delta_k lies in Fp for every even k.
* delta_6 = -1, and delta_(k+6) = -delta_k.
* lambda = p mod r has multiplicative order exactly 12, and lambda^6 = -1 mod r.

On GT, phi(x) = x^lambda, so phi is a group automorphism of order 12. Every
x != 1 has an orbit of exactly 12 elements, and 1 is fixed. The orbit count is

    N = (r - 1)/12 + 1.

### 2.2 Class key and representative

Define the key K(y) = Re(w-coefficient of y), an Fp integer in [0, p). For a
state u, the 12 keys K(phi^k(u)) for k = 0..11 are:

| k | Key | Fp multiplications |
| --- | --- | --- |
| 0 | Re(b0) | free |
| 6 | p - Re(b0) | free |
| 2, 4 | Re(b0) · delta_k (delta_k in Fp) | 1 each |
| 1, 3, 5 | Re(b0)·Re(delta_k) + Im(b0)·Im(delta_k) | 2 each |
| k + 6 | negation of the k value, since delta_(k+6) = -delta_k | free |

That is 8 Fp multiplications in total. The class representative is phi^k*(u),
where k* = argmin_k K(phi^k(u)).

* Ties between two conjugates occur with probability about 2^-374 per step.
  They are broken deterministically by comparing the full canonical encodings,
  which costs at most 12 Frobenius maps (about 4 units), i.e. negligible.
* b0 = 0 for a nonidentity element of GT has probability about 2^-762. In that
  case the same rule is applied to b1.

The minimal key, and hence the partition index and the distinguished-point
property derived from it, is a class invariant.

### 2.3 Twisted-state walk (no canonicalization of the state)

Precompute a table of R = 2^16 entries M_j = g^(alpha_j) z^(beta_j), each stored
with all 12 Frobenius twists T[j][k] = phi^(-k)(M_j). The exponents of T[j][k]
are alpha_j · lambda^(-k) and beta_j · lambda^(-k) mod r.

The state is an arbitrary element u of its class, together with exponents
(A, B) such that u = g^A z^B. The class representative is phi^k*(u). One step:

1. Read j = H(min key), 16 bits of the minimal key.
2. Set u' = u · T[j][k*] (one dense Fp12 multiplication), A' = A + alpha_j ·
   lambda^(-k*), B' = B + beta_j · lambda^(-k*).
3. Compute the 12 keys of u' and the new k*.

Correctness: phi^k*(u) · M_j = phi^k*(u · phi^(-k*)(M_j)) = phi^k*(u'). So u'
lies in the class of (representative · M_j), which is exactly the class walk
"next = rep · M_j".

Well-definedness: replace u by u2 = phi^s(u) from the same class. Then
k2* = k* - s, and u2 · T[j][k2*] = phi^s(u'), which is in the same class. The
walk is therefore a function on classes, and the state is never canonicalized.

When a distinguished point is stored, the representative's exponents are
(A · lambda^k*, B · lambda^k*) mod r. This is computed once per distinguished
point.

### 2.4 Fruitless cycles

A 2-cycle needs y1 = phi^s(y0 · M_j) followed by H(y1) = j and a return to the
class of y0. In additive notation, the class of y2 is that of
lambda^(t+s) y0 + (lambda^(t+s) + lambda^t) m_j. For generic m_j this returns to
the class of y0 only when lambda^s = -1, i.e. s = 6. So 2-cycles have
probability about 1/(12R) before mitigation.

Mitigations:

* Wiener–Zuccherato look-ahead: if H(new class) = j, retry with j+1, j+2, and
  so on. This removes the (j, j) 2-cycles, leaving rates around 1/R^2 or below
  for longer cycles. A retry costs one extra multiplication with probability
  about 1/R = 2^-16.
* With R = 2^16, we conservatively assume a fruitless-cycle entry rate of at most
  2^-30 per step. This is about 12 times the negation-map 4-cycle rate
  1/(4R^2) = 2^-34.
* Walks run for 20/θ steps with θ = 2^-20. A walk that reaches 20 · 2^20 steps
  without a distinguished point is abandoned (its steps are charged) and
  restarted. The probability of getting stuck before a distinguished point is
  about 2^20 · 2^-30 = 2^-10, and the waste is at most 20 · 2^20 steps per stuck
  walk. Expected waste is about 2%; we charge 3%.

### 2.5 Collision search

Run w ≤ 2^40 walks from random starts g^A0 z^B0 (van Oorschot–Wiener). A point is
distinguished when 20 designated bits of the minimal key are zero. Store
(256-bit fingerprint of the representative, A, B), 96 bytes, in a hash table,
and restart the walk after each distinguished point.

A fingerprint match is confirmed by recomputing the representatives; false
matches total about 2^-43. Equal representatives give g^A z^B = g^A' z^B', so
a = (A - A')/(B' - B) mod r when B != B'. B = B' has probability about 1/r.

## 3. Steps and success probability

Model the class walk as a random mapping on N classes. Let

    k0 = sqrt(2 N ln(1/0.61)),   log2 k0 = 125.627787

At k0 steps the random-mapping collision probability is 0.39. We budget
**1.25 · k0 useful steps**, where the random-mapping success is
1 - 0.61^1.5625 = 0.538. The margin also covers a walk up to about 25% less
efficient than random, and detection lag.

Support for the heuristic:

* Teske (2001): r-adding walks with r ≥ 16 behave close to random.
* Duursma–Gaudry–Morain and Bos et al. report that the full automorphism-group
  speedup (√m) is realized in practice.
* Our toy experiment `experiment.py` uses exactly this twisted-state algorithm,
  in the prime-order-ℓ subgroup of F_q* with an order-12 automorphism x -> x^lam
  (12 | ℓ - 1). It uses 64 partitions, look-ahead and stuck-walk restarts, and
  counts every step. All runs were local, at seeds as shown.

| ℓ bits | walks / θ | trials | mean steps / √(πN/2) | success within 1.25 k0 |
| ---: | ---: | ---: | ---: | ---: |
| 24 | 8 / 2^-4 | 100 (seed 41) | 0.986 | 0.49 |
| 30 | 8 / 2^-4 | 8 × 100 (seeds 2–9) | 0.926–1.060 (mean about 0.99) | 0.47–0.59 (mean 0.54) |
| 30 | 8 / 2^-4 | 100 (seed 1) | 0.822 | 0.69 |

These agree with the random-mapping prediction 0.538, so the √12 class speedup
is realized at toy scale. In the toy, the fruitless-cycle rate (R = 64) is far
higher than at production size and distinguished-point lag is a larger share of
the work, so the toy is pessimistic on both effects. It remains toy-scale
evidence only.

## 4. Cost accounting (dense Fp12 multiplication = 1 unit)

### Per step

* One dense Fp12 multiplication u · T[j][k*]: 1 unit (the reference operation).
* Other work, itemized in scalar instructions, using the pinned blst_fp_mul cost
  of 383 instructions:

| Item | Instructions |
| --- | ---: |
| 8 Fp multiplications for the 12 keys | 3,064 |
| 3 modular additions and 6 negations (6 limbs each) | about 200 |
| 12 six-limb comparisons and argmin | about 250 |
| Partition index and distinguished-point test | about 10 |
| Table-entry load (72 words) | about 80 |
| Two 4-limb modular exponent additions (pre-twisted exponents) | about 80 |
| State stores and loop control | about 100 |
| Look-ahead extra multiplication, probability 2^-16 | negligible |
| **Total** | **about 3,800** |

We charge 5,000 instructions = 0.25 units, so the per-step cost is
**s = 1.25 units**.

### Steps

1.25 · k0 useful steps, times 1.03 for abandoned walks.

### Additive work (all charged)

| Item | Bound |
| --- | --- |
| Pairing g = e(P, Q0) | ≤ 2^14 units: under 126 Miller steps at under 64 units each, plus Fermat inversion and the final exponentiation by binary powering (≤ 2 × 4572 products each) |
| Table: 2^16 entries, each 2 exponentiations of 255 bits (≤ 1020 units) and 11 Frobenius maps (≤ 1 unit each) | ≤ 2^26.1 |
| Walk starts: ≤ 2^40 walks × (1020 + 1.25) units | ≤ 2^50 |
| Distinguished-point lag (each walk's tail, and detection delay) | ≤ w/θ · 1.25 ≤ 2^60.4 |
| Hash-table inserts and representative exponents: about 2^106.3 distinguished points × 2^13 instructions / 20,000 | ≈ 2^105, below 2^-21 of the main term |
| Output: a mod r, [a]Q0 in G2 (< 2^22 instructions), verification pairing (2^14 units) | negligible |

### Total

    W = 1.25 · 1.03 · 1.25 · k0 + (additive terms below 2^-20 of the main term)
    log2 W = 125.627787 + log2(1.609375) + (< 2^-19)
           = 125.627787 + 0.686500 + ...
           = 126.3143

The submitted bound is **time_log2 = 126.32**, rounded up.

## 5. Memory, preprocessing, advice

* Peak memory is about 2^106.32 distinguished points × 96 bytes, plus the
  twisted table (2^16 · 12 · 640 bytes ≈ 503 MB) and 2^40 walk states of about
  768 bytes. That totals **2^112.9 bytes**; claimed memory_log2_bytes = 113.
  Memory is reported only.
* Preprocessing (pairing and table) is at most 2^26.2 units; preprocessing_log2
  ≤ 27. It is included in W.
* Nonuniform advice: none. The Frobenius constants delta_k and lambda are
  derived from the public p and r; each costs at most one 4572-bit
  exponentiation, which is negligible and included in the setup allowance.

## 6. Heuristics and limitations

* **Heuristic:** the 2^16-adding orbit walk with look-ahead behaves like a random
  mapping on the N classes, to within the 25% step margin. Support: the
  literature above, and the toy runs at a ratio of about 1.0. It is not proven
  at production size.
* **Heuristic:** fruitless-cycle entry rates are at most 2^-30 per step at
  R = 2^16. A misestimate costs only charged waste, because walks are abandoned
  after 20/θ steps. Even a rate of 2^-26 would add at most about 30% to the
  step count, i.e. under 0.4 bits.
* The key uses one coordinate, which is enough to choose a canonical conjugate.
  The full representative is never materialized during the walk.
* Not claimed: TNFS/SexTNFS for DLP in Fp12 (p^12 is about 4572 bits). Published
  estimates for BLS12-381 are around 2^120–2^126 in those papers' own cost units.
  A ledger converted to this cost model is future work.
* Analytic claim only; no production-size execution.

References:

* J. Pollard, Math. Comp. 1978.
* P. van Oorschot and M. Wiener, J. Cryptology 1999.
* M. Wiener and R. Zuccherato, SAC 1998.
* I. Duursma, P. Gaudry and F. Morain, ASIACRYPT 1999.
* R. Gallant, R. Lambert and S. Vanstone, Math. Comp. 2000.
* E. Teske, Math. Comp. 2001.
* J. Bos, T. Kleinjung and A. Lenstra, ANTS 2010.
* D. Bernstein, T. Lange and P. Schwabe, PKC 2011.
