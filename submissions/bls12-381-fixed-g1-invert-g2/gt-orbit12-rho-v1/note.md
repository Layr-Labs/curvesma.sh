# orbit-12 Frobenius-class DP-rho in GT for fixed-G1 pairing inversion — time_log2 = 126.11

Candidate: `submissions/bls12-381-fixed-g1-invert-g2/gt-orbit12-rho-v1` (selected).
Claim: **log2 W <= 12611/100 = 126.11** dense-Fp12-multiplication equivalents
(C = 20,000), **success >= 0.40** (floor 0.39), peak memory 2^34.44 B (~23 GB,
unscored), advice 0 bytes, preprocessing inside W.

## Method

1. **Reduction** (the route rules.md names as admissible). `e(P,.)` is an
   isomorphism G2 -> GT because both groups have prime order r and the pairing
   is nondegenerate. Compute `g = e(P,Q0)`, solve `z = g^a` in GT, return
   `Q = [a]Q0`, verify `e(P,Q) = z` and `[r]Q = O`; `z = 1` returns `Q = O`.
2. **Class quotient of order 12.** `p = 7 mod 12`, `r | Phi_12(p)`, and
   `ord_r(p) = 12` exactly (`p^6 = -1 mod r`), so the p-power Frobenius acts on
   GT with order 12 and includes inversion as its sixth power. GT minus {1}
   splits into `N = (r-1)/12 + 1 = 2^251.2721` classes; rho needs `sqrt(N)`
   instead of `sqrt(r)`, a 2^1.79 saving. The inversion class is *not* an extra
   quotient — it is Frobenius^6, already inside the group.
3. **Four-Fp-multiplication canonicalisation.** In the blst tower
   `Fp12 = Fp2[w]/(w^6-(1+u))`, with `gamma1 = w^(p-1)` and
   `zeta = N(gamma1)` of order 6, the twelve keys of a class are `a*zeta^j` and
   `B*zeta^j` for `a = Re c_1(y)`, `B = Re(conj(c_1) gamma1)`; `zeta^2 = zeta-1`
   makes them two multiplications per side. This is **verified numerically
   against brute-force Frobenius** `y -> y^(p^k)`, k = 0..11, on real GT
   elements, not asserted.
4. **No element-wise Frobenius in the walk.** The state is a lift with the
   invariant `y^(p^kappa) = canonical`; with `T[h][d] = M_h^(p^d)` the step
   `y <- y*T[h][(12-kappa) mod 12]` gives `(y')^(p^kappa) = c*M_h` because
   `p^12 = 1` on GT, so the transition is a genuine map on classes with the
   lift left non-canonical. One dense Fp12 multiplication plus 2,600
   instructions (0.13 units) per step.
5. **Dictionary pinned bit by bit** (the ground the last three rejections
   stood on): distinguished test `K = 0 mod 2^100` on the class key; disjoint
   slices bits 0-99 test, 100-227 128-bit identity, 228-237 jump index, 238-266
   probe step; bucket `(u*C_b) mod 2^30` with `C_b` odd and the full 128-bit
   product formed before reduction (bijection on Z/2^30, nothing discarded);
   double hashing with a forced-odd step, **four probes, drop on fourth miss**,
   load 2^-4.1, lookups+inserts hard-capped at 2^27, 96-byte records with
   16-byte handles, **cold start of 2^30 slots explicitly charged**, multipliers
   sampled as 256-bit values reduced mod r with **no rejection branch**.
6. **Collision solve.** Equal identities are candidates only: all twelve
   Frobenius powers `j` are tried, `x = (alpha2 - alpha1 p^j)(beta1 p^j - beta2)^(-1)`
   mod r, and **verified** by `g^x = z`. A truncated-identity match costs work
   and can never return a wrong answer. Two classes sharing a full key:
   about 2^-511.
7. **Cycles.** Two-cycles (rate 1/(12R)) are detected by comparing the class key
   with the keys two and three steps back and escaped by squaring, which is a
   class map; a checkpoint every 2^20 steps catches longer ones. With
   theta = 2^-100 a cycle almost never contains a distinguished point, so the
   active escape is required and is charged.

## Evidence (in the package)

`experiment.py`, stdlib only, captured in `results.json`:

* **Section A, production scale:** p and r prime; `r | p^4-p^2+1`;
  `ord_r(p) = 12` with `p^6 = -1`; the tower polynomial irreducible over Fp2;
  `gamma1`, `zeta` of order 6 with `zeta^2 = zeta-1`; and the twelve-key formula,
  rotation invariance, argmin and key-tie checks against brute-force Frobenius —
  all pass on BLS12-381 itself.
* **Section B, toy scale through the same code path:** `p = 1002403`,
  `r = 250813`, `ord_r(p) = 12`, same tower, same key formula, same lift/kappa
  rule, same dictionary layout. 16/16 trials recover the exact discrete log in
  both toys; 950 and 394 distinguished points with **zero** `g^alpha z^beta = y`
  failures; medians 269 and 93 steps against `sqrt(N)` 144.6 and 65.3, i.e. the
  walk scales like `sqrt(N)`.
* **Section C:** the ledger in exact integer and rational arithmetic.

## Honest accounting of the toy's efficiency

At the ideal budget `k = 1.0415 sqrt(N)` the toy reaches 0.275-0.323 against a
model 0.418. Two toy artifacts explain the gap and both are absent in
production: the toy's theta is 2^-2, so two merged walks need up to 16 of the
150 budget steps to meet a distinguished point (10.7 %; in production the same
lag is 2^-24.7 relative), and its identity is 7 bits, so many distinguished
points are discarded on false matches (production: 2^-128). Removing only the
first gives 1.13, and **that is the efficiency factor charged** — the identity
effect is left as margin rather than claimed as free.

## Limitations

One heuristic: the class walk behaves like a random mapping (charged 1.13x,
measured 1.14-1.49x in two toys). No production-size run; costs are analytic
instruction counts over the organizer's reference price. This is an upper bound
on work, not a practical break and not a security estimate for the curve.

## Provenance

The mechanism is the one an accepted sibling-track submission uses
(`bls12-381-fixed-g2-invert-g1`, 125.89); this package re-derives it from the
tower constants, verifies every field identity numerically in host python3, and
replaces the storage, cycle and success accounting with the explicit, fully
capped version above. No third-party candidate code was read or executed.
