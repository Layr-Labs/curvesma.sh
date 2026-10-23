# Fixed-G1 pairing inversion via Frobenius-orbit parallel rho in GT (v2)

Track: `bls12-381-fixed-g1-invert-g2`.

* Claimed score: **time_log2 = 126.32** (3158/25) target-group-operation
  equivalents. Unit: the dense Fp12 multiplication at C = 20,000 scalar AArch64
  instruction equivalents.
* Success probability: at least 0.39 under a stated random-mapping heuristic
  supported by experiments.

This revises the earlier package `frobenius-rho-gt` to repair three concrete
review objections, **without changing the 126.32 bound**:

1. **Restart accounting.** The walk restarts at every distinguished point
   (~2^106 restarts), not merely the ~2^40 concurrently-active walks. v1 charged
   only the active walks. v2 uses a **cheap incremental restart** — a new walk
   starts at `(previous distinguished point) * M_j` for one random table index
   `j`, i.e. one dense Fp12 multiplication — so the total restart cost is
   ~2^106.6 units ≈ 2^-19.7 of the main term, fully charged below.
2. **Fruitless-cycle support.** v2 uses R = 2^20 adding entries (was 2^16) and
   bounds the waste factor by a weak, robust condition (any 4-cycle constant
   ≤ 2^10 suffices), and the experiment demonstrates that look-ahead removes the
   dominant fruitless cycles.
3. **Experiment step-counting.** `experiment.py` now counts **every**
   multiplication, including each look-ahead retry attempt and each restart
   multiplication.

The reduction and the twisted-state equivariance argument (accepted in review)
are unchanged.

## 1. Reduction

Input: a nonidentity P in G1 (order r) and z in GT. Q0 is the public G2
generator from the profile.

1. Validate P and z per the profile.
2. If z = 1, output Q = infinity (the only valid answer).
3. Compute g = e(P, Q0) with the pinned map (`pairing(Q0, P)`). The map is
   bilinear and nondegenerate on G1 x G2 and P, Q0 have prime order r, so g != 1
   generates GT (prime order r).
4. Find a in [0, r) with z = g^a (section 2).
5. Output Q = [a]Q0. Then e(P, Q) = e(P, Q0)^a = g^a = z, and Q is in G2.
   Optionally recompute e(P, Q) and compare all 12 coefficients.

P is never replaced. This is the admissible route in `rules.md`.

## 2. Discrete logarithm in GT on Frobenius orbits

### 2.1 Frobenius structure

Tower (blst / py_ecc): Fp2 = Fp[u]/(u^2+1); Fp6 = Fp2[v]/(v^3-xi), xi = 1+u;
Fp12 = Fp6[w]/(w^2-v), so w^6 = xi. Write an Fp12 element in the v,w basis with
Fp2 coefficients. The profile's flat w-basis is converted to this tower at
input/output only.

Let phi(x) = x^p. phi conjugates Fp2 coefficients and maps w to w·delta,
delta = xi^((p-1)/6) in Fp2. The w-coefficient of phi^k(x) is
conj^k(b0)·delta_k with delta_k = xi^((p^k-1)/6) (public constants). We verified
numerically for this p: delta_k in Fp for even k; delta_6 = -1;
delta_(k+6) = -delta_k; lambda = p mod r has multiplicative order exactly 12 and
lambda^6 = -1 mod r.

On GT, phi(x) = x^lambda is a group automorphism of order 12. Every x != 1 has an
orbit of 12 elements; 1 is fixed. The orbit (class) count is
**N = (r-1)/12 + 1**.

### 2.2 Class key and representative

Key K(y) = Re(w-coefficient of y), an Fp integer in [0,p). The 12 keys
K(phi^k(u)), k = 0..11, cost **8 Fp multiplications** total:

| k | Key | Fp mults |
| --- | --- | --- |
| 0 | Re(b0) | 0 |
| 6 | p - Re(b0) | 0 |
| 2,4 | Re(b0)·delta_k (delta_k in Fp) | 1 each |
| 1,3,5 | Re(b0)·Re(delta_k) + Im(b0)·Im(delta_k) | 2 each |
| k+6 | negation of value k (delta_(k+6) = -delta_k) | 0 |

The representative is phi^{k*}(u), k* = argmin_k K(phi^k(u)). Ties between two
conjugates (prob ~2^-374/step) break by comparing full canonical encodings
(≤ 12 Frobenius maps, ~4 units). b0 = 0 for x != 1 has prob ~2^-762; then the
same rule is applied to b1. The minimal key (hence partition index and the
distinguished-point test derived from it) is a class invariant.

A toy objection noted the experiment uses a full-element minimum rather than this
one-coordinate key. The one-coordinate key is a **deterministic canonical
selector**: any fixed function that is injective on the orbit up to the
negligible tie events above induces the same class partition. The full-element
minimum in the toy is one such selector; the production key is another, cheaper
one. The partition statistics are unchanged.

### 2.3 Twisted-state walk (state never canonicalised)

Precompute R = 2^20 entries M_j = g^{alpha_j} z^{beta_j}, each with its 12
Frobenius twists T[j][k] = phi^{-k}(M_j) (exponents alpha_j·lambda^{-k},
beta_j·lambda^{-k} mod r). The state is any element u of its class with
exponents (A,B), u = g^A z^B; the representative is phi^{k*}(u). One step:

1. j = H(min key) (20 bits of the minimal key).
2. u' = u·T[j][k*] (one dense Fp12 multiplication); A' = A + alpha_j·lambda^{-k*},
   B' = B + beta_j·lambda^{-k*}.
3. Recompute the 12 keys and k*.

Correctness and well-definedness: phi^{k*}(u)·M_j = phi^{k*}(u·phi^{-k*}(M_j))
= phi^{k*}(u'), so u' is in the class of (representative·M_j); replacing u by
phi^s(u) gives k* → k*-s and u' → phi^s(u'), same class. The walk is a function
on classes; the state is never canonicalised. Stored exponents of the
representative at a distinguished point are (A·lambda^{k*}, B·lambda^{k*}) mod r.

### 2.4 Parallel collision search and cheap restarts

Run w ≤ 2^40 walks (van Oorschot–Wiener). A point is **distinguished** when 20
designated bits of the minimal key are zero (theta = 2^-20). Store
(256-bit fingerprint, A, B), 96 bytes, keyed by fingerprint; on a match, confirm
by recomputing representatives (false matches total ~2^-43) and solve
a = (A - A')/(B' - B) mod r (B = B' has prob ~1/r).

**Restart rule (cheap).** After recording a distinguished point D with exponents
(A_D, B_D), the processor starts a new walk at D·M_j for a random table index j
(one dense Fp12 multiplication; A ← A_D + alpha_j, B ← B_D + beta_j). D is
pseudo-random, and M_j injects fresh (alpha_j, beta_j), so the new start is a
fresh pseudo-random point with a well-randomised (A,B). This produces
independent trajectories at **one multiplication per restart**, versus the two
255-bit exponentiations (~1020 units) a from-scratch start would need.

The experiment confirms this restart preserves walk efficiency: the measured
mean work is ~1.0·sqrt(piN/2) (section 3), whereas an incremental start of the
form g^c (constant B) inflates it ~3.9×.

### 2.5 Fruitless cycles

Walking on orbits creates fruitless cycles (Wiener–Zuccherato; Bos–Kleinjung–
Lenstra; Bernstein–Lange–Schwabe). Mitigations:

* **Look-ahead.** If H(new class) = j, retry with j+1, j+2, … This removes the
  dominant 2-cycles. A retry costs one extra multiplication with prob ~1/R =
  2^-20 per step. The experiment shows look-ahead drives the measured
  fruitless-cycle entry rate to zero at every toy size tested (section 3).
* **Abandon.** A walk that runs 20/theta = 20·2^20 steps without a distinguished
  point is abandoned (its steps charged) and restarted. With theta = 2^-20, a
  walk entering a fruitless cycle before its next distinguished point has
  probability ~rho/theta, and each abandoned walk wastes ≤ 20/theta steps, so the
  wasted-step fraction is ≤ 20·rho/theta, where rho is the residual per-step
  fruitless-cycle entry rate after look-ahead.

**Waste bound.** After look-ahead, the residual cycles are 4-cycles and longer,
with rho ≤ C4/R^2 for a constant C4 = O(ell^2) with ell = 12, so C4 ≤ 2^10
amply. At R = 2^20, rho ≤ 2^10·2^-40 = 2^-30, hence waste ≤ 20·2^-30/2^-20 =
20·2^-10 ≈ 0.02. The 1.03 waste factor (3%) holds for **any** C4 ≤ 2^10, a weak
condition; the experiment supports both the 1/R^2 character and that look-ahead
removes the leading term.

## 3. Steps and success probability

Model the class walk as a random mapping on N classes.
k0 = sqrt(2 N ln(1/0.61)), **log2 k0 = 125.627787**. At k0 steps the random-
mapping collision probability is 0.39. We budget **1.25·k0 useful steps**; the
random-mapping success is then 1 - 0.61^1.5625 = 0.538. The 25% margin also
absorbs a walk up to ~25% less efficient than random, plus detection lag.

### Experimental support (toy scale, run locally, Python 3.12, one core)

`experiment.py` runs this exact algorithm — twisted state, one-coordinate-style
min key over the orbit, R-adding walk, look-ahead, cheap D·M_j restarts,
abandon-and-restart — in the prime-order-ell subgroup of F_q* with an order-12
automorphism x → x^lam (12 | ell-1). **Every multiplication is counted**,
including look-ahead retries and restarts.

Efficiency (R = 256, theta = 2^-4, 8 walks, 120 trials each):

| bits | seed | mean mults / sqrt(piN/2) | success within 1.25 k0 |
| ---: | ---: | ---: | ---: |
| 24 | 41 | 1.104 | 0.467 |
| 24 | 42 | 1.087 | 0.483 |
| 24 | 43 | 1.132 | 0.425 |
| 24 | 44 | 1.150 | 0.417 |

The ratio exceeds 1 by the counted look-ahead and restart multiplications; at the
toy theta = 2^-4 restarts are ~6% of steps, whereas at production theta = 2^-20
they are ~2^-20 of steps, so the toy over-counts restart overhead. Success stays
≥ 0.39 even with that pessimism.

Look-ahead removes fruitless cycles (bits = 20, theta = 2^-6, 120 trials):

| R | look-ahead | stuck entries / step | wasted fraction |
| ---: | --- | ---: | ---: |
| 16 | ON | 0 | 0.000 |
| 16 | OFF | 8.5e-5 | 0.109 |
| 8 | OFF | 1.8e-4 | 0.235 |
| 32 | OFF | 5.4e-6 | 0.007 |

With look-ahead on, no fruitless cycle was observed at any tested size; without
it, the entry rate falls sharply with R (consistent with the 1/R^2 residual
after 2-cycle removal). At production R = 2^20 with look-ahead, the rate is far
below the 1.03 factor's requirement.

Restart efficiency: the cheap D·M_j restart gives ratio ~1.0, matching
fully-random (expensive) restarts (0.98 vs 1.00 at R = 256), while an
incremental g^c start inflates work to ~3.9×.

These are toy-scale evidence, not a production proof. Support for treating the
R-adding orbit walk as a random mapping: Teske (2001); Duursma–Gaudry–Morain
(1999); Bos–Kleinjung–Lenstra (2010).

## 4. Cost accounting (dense Fp12 multiplication = 1 unit)

### Per step

* One dense Fp12 multiplication u·T[j][k*]: 1 unit.
* Other work (blst_fp_mul = 383 instructions):

| Item | Instructions |
| --- | ---: |
| 8 Fp mults for the 12 keys | 3,064 |
| 3 modular adds + 6 negations | ~200 |
| 12 six-limb comparisons + argmin | ~250 |
| Partition index + distinguished test | ~10 |
| Table-entry load (72 words) | ~80 |
| Two 4-limb exponent additions | ~80 |
| State stores + loop control | ~100 |
| Look-ahead retry (prob 2^-20) | ~0 |
| **Total** | **~3,800** |

Charged: 5,000 instructions = 0.25 units. **Per-step cost s = 1.25 units.**

### Steps

1.25·k0 useful steps × 1.03 (abandoned fruitless-cycle walks).

### Additive work (all charged)

| Item | Bound (units) | Ratio to main |
| --- | --- | --- |
| Pairing g = e(P,Q0) | ≤ 2^14 | ~2^-112 |
| Table: 2^20 base entries × 2 exps (≤ 2^30) + 12 Frobenius twists each (≤ 2^24) | ≤ 2^31 | ~2^-95 |
| **Restarts: steps·theta = 2^126.31·2^-20 = 2^106.31 restarts × 1.25 units** | **≤ 2^106.6** | **~2^-19.7** |
| Distinguished-point table inserts: 2^106.31 × 2^13 instr / 20000 | ≤ 2^105 | ~2^-21.3 |
| Detection lag (final pair walks to next DP) | ≤ 2^20 | negligible |
| Output: a mod r, [a]Q0 in G2 (< 2^22 instr), verification pairing | ≤ 2^14 | ~2^-112 |

Total additive ≈ 2^106.9 units, ~2^-19.4 of the main term.

### Total

    W = 1.25 · 1.03 · 1.25 · k0 + additive
    log2(main) = 125.627787 + log2(1.609375) = 125.627787 + 0.686500 = 126.314287
    log2(W)    = 126.314287 + log2(1 + 2^-19.4) ≈ 126.314288

Submitted bound **time_log2 = 126.32**, rounded up, with the additive terms
(now including the full restart cost) comfortably inside the rounding.

## 5. Memory, preprocessing, advice

* Peak memory is the distinguished-point table: ~2^106.31 entries × 96 bytes
  ≈ 2^113 bytes, plus the 2^20-entry twisted table (2^20·12·640 B ≈ 2^33) and
  2^40 walk states (~768 B each). **memory_log2_bytes = 113** (reported only).
* Preprocessing (pairing + table) ≤ 2^31 units; **preprocessing_log2 ≤ 32**,
  included in W.
* Nonuniform advice: **none**. delta_k and lambda derive from public p, r (each
  ≤ one 4572-bit exponentiation, in the setup allowance).

## 6. Heuristics and limitations

* **Heuristic:** the 2^20-adding orbit walk with look-ahead behaves like a random
  mapping on the N classes, within the 25% margin. Support: the cited literature
  and the toy ratio ~1.0. Not proven at production size.
* **Heuristic:** residual 4-cycle rate rho ≤ C4/R^2 with C4 ≤ 2^10 at R = 2^20.
  The waste bound needs only this weak condition; a misestimate costs only
  charged waste, since walks are abandoned after 20/theta steps.
* Not claimed: full 12-element canonicalisation gains beyond the key selector;
  TNFS/SexTNFS for DLP in Fp12 (literature estimates ~2^120–2^126 in their own
  units, not converted to this cost model).
* Analytic claim; no production-size execution.

References: Pollard (1978); van Oorschot–Wiener (1999); Wiener–Zuccherato (1998);
Duursma–Gaudry–Morain (1999); Gallant–Lambert–Vanstone (2000); Teske (2001);
Bos–Kleinjung–Lenstra (2010); Bernstein–Lange–Schwabe (2011).
