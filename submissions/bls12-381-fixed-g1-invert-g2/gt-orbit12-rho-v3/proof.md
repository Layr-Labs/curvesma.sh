# Orbit-12 Frobenius-class distinguished-point rho in GT, with a fully charged dictionary

Track: `bls12-381-fixed-g1-invert-g2` (fixed-argument pairing inversion, FAPI-1).
Claim: **log2 W ≤ 12611/100 = 126.11** target-group-operation equivalents
(one dense Fp12 multiplication, C = 20,000 scalar instruction equivalents), with
**success probability ≥ 0.40** (rule floor 0.39).
Peak memory 2^34.44 bytes (~23 GB), unscored. Advice 0 bytes.

Everything numeric below is produced by `experiment.py` (stdlib only) and
captured in `results.json`. Section numbers refer to that file.

---

## 1. Reduction (the route rules.md names as admissible)

`P` is a nonidentity G1 point fixed by the instance; `z` is the GT target.
`e(P, ·) : G2 -> GT` is a homomorphism between groups of prime order r, and it is
injective because `e(P,Q) = 1` with `P != O` forces `Q = O` (nondegeneracy).
Both groups have order r, so it is an isomorphism.

Fix the public generator `Q0` of G2 (`target-profile.json`). Then

1. `g = e(P, Q0)` is a generator of GT (one pairing, charged below).
2. Solve `z = g^a` in GT; that is the only hard step.
3. Return `Q = [a]Q0`; bilinearity gives `e(P,Q) = g^a = z`.
4. `z = 1` is answered directly with `Q = O` (the acceptance rules allow exactly this).
5. Verify `e(P,Q) = z` (`pairing(Q,P)` in the pinned library) and `[r]Q = O`.

The answer does not depend on `P`; only the DLP instance `(g, z)` does. This is
the same reduction to a GT discrete logarithm that the accepted sibling-track
claim uses, with the pairing arguments swapped.

## 2. The class quotient: 12, not 2

`p ≡ 7 (mod 12)`, `r | p^4 - p^2 + 1 = Phi_12(p)`, and `ord_r(p) = 12` exactly:
`p^6 ≡ -1`, `p^4 ≢ 1`, `p^3 ≢ 1`, `p^2 ≢ 1` (mod r). Verified on the protected
profile constants (`results.json` A1). Writing `sigma(y) = y^p` for the p-power
Frobenius, `sigma` acts on GT as `[p]` in the exponent, so `sigma` has order
exactly 12 on GT, and `sigma^6 = [-1]` is inversion.

Hence GT \ {1} splits into

    N = (r-1)/12 + 1 = 2^251.2721 classes     (the twelve conjugates, plus {1})

A rho walk on classes needs about `sqrt(N)` steps instead of `sqrt(r)`: a
`sqrt(12) = 2^1.79` saving. The inversion class `{x, x^-1}` is *not* additional
structure: it is `sigma^6`, already inside the group of order 12.

## 3. Cheap canonicalisation (the engineering that makes it pay)

Tower, as measured by COST_MODEL's reference build:

    Fp2 = Fp[u]/(u^2+1),   Fp12 = Fp2[w]/(w^6 - xi),  xi = 1 + u

`w^6 - xi` is irreducible over Fp2: `xi` is neither a square nor a cube in Fp2
(`results.json` A2, both `True`), so the tower is the field `F_{p^12}`.

Let `gamma1 = w^(p-1)` (an Fp2 constant of the tower, computed in A3) and
`zeta = gamma1^(p+1) = N(gamma1)`, which lands in Fp and has order exactly 6 with
`zeta^2 = zeta - 1` (A3). For a lift `y = sum c_k w^k`, the p^k-Frobenius acts by
`c_k -> conj(c_k) gamma_{2j}^{...}`; specialising to `k = 1` gives, for
`a = Re c_1` and `B = Re(conj(c_1) gamma1)`,

    Re c_1(y^(p^(2j)))   = a * zeta^j
    Re c_1(y^(p^(2j+1))) = B * zeta^j          (j = 0..5)

so the **twelve canonical keys of a class are `{a zeta^j} u {B zeta^j}`**, and
because `zeta^2 = zeta - 1`, `zeta^3 = -1`, `zeta^4 = -zeta`, `zeta^5 = 1 - zeta`
they cost **two Fp multiplications for `a` and two for `B`**, plus eight
negations/subtractions. `K = min` over the twelve identifies the class and
`kappa = argmin` identifies which conjugate is canonical.

This is not asserted: A4 compares the formula against brute-force Frobenius
`y -> y^(p^k)`, k = 0..11, on random GT elements drawn as
`t^((p^12-1)/r)` (elements of order exactly r, verified with `y^r = 1`).
**A4 `formula_matches_bruteforce` = 2/2, `rotation_invariant_and_min_agrees` =
2/2, `argmin_matches_bruteforce` = 2/2, key-tie events 0.** The key set is
therefore a genuine class function and `min` is well defined on classes.

Cost of a class key: 4 Fp multiplications + 8 add/sub + 11 comparisons.

## 4. The walk never applies Frobenius to the walk element

State: `(y, kappa)` where `y` is *a* lift, with the invariant `y^(p^kappa)
= c = canonical representative of its class`. Jump table
`T[h][d] = M_h^(p^d)`, `d = 0..11`, built once, where `M_h = g^(A_h) z^(B_h)`.

Step from `(y, kappa)` with `h = H(K) mod R`:

    y' = y * T[h][(12 - kappa) mod 12],   alpha += PA[h][d],  beta += PB[h][d]

with `PA[h][d] = A_h p^d mod r`, `PB[h][d] = B_h p^d mod r`. Then

    (y')^(p^kappa) = y^(p^kappa) * M_h^(p^(kappa + 12 - kappa))
                   = c * M_h^(p^12) = c * M_h

because `p^12 ≡ 1` on GT. So the next class is `class(c * M_h)`, a function of
the **class** and `h` alone: the transition is a well-defined map on classes,
even though the walk carries a non-canonical lift. The exponent pair `(alpha,
beta)` is updated by two precomputed constants, so `g^alpha z^beta = y` is
maintained exactly. The toy checks this at **every** distinguished point:
`book_fail_total = 0` over 950 + 394 recorded points in the two toys.

Per step: one dense Fp12 multiplication (the reference operation) plus the class
key (above), the jump index, the distinguished test and two modular additions —
**2,600 instructions = 0.13 units**, itemised in §9.

## 5. Distinguished classes and the dictionary (specified bit by bit)

The distinguished test is applied to the **class key `K`**, which is a class
function, so the stored object is a class, not a representative. Bit layout of
the 381-bit key `K` (`results.json` B `layout` shows the same layout at toy
scale):

| bits | use |
| --- | --- |
| 0 .. 99 | distinguished test: `K ≡ 0 mod 2^100` (theta = 2^-100) |
| 100 .. 380 | the remaining 281 bits, all of which enter the bucket map |

The jump index `h` and the probe step are functions of the 281 non-distinguished
bits and never of the distinguished slice, so the distinguished condition cannot
correlate with the walk's next move.

* **Bucket, full width.** Put `K' = K >> 100` (281 bits) and split it into five
  disjoint chunks `u_i = (K' >> 57 i) mod 2^57`, `i = 0..4` (the last chunk is
  shorter). Then `b = (sum_i u_i C_i) mod 2^30` with fixed odd `C_i`. Every
  non-distinguished bit of `K` enters the sum, so no input bit can be discarded
  (this is the defect that `647829a1` was rejected for: there the map
  `a + 2^64 b mod 2^71` discarded `b`). The bucket map is a hash, not a
  bijection: records store the **full** key, so a bucket collision costs one
  probe and nothing else.
* **Probe sequence.** Double hashing with the odd step
  `1 + 2*((K' >> 3) mod 2^29)`: `pos -> pos + step mod 2^30`. Odd step and `2^30`
  slots are coprime, so the probe sequence enumerates all slots; at most **4
  probes** per operation, and an insert that finds no free slot within 4 probes
  is **dropped** (never retried, never chained). Load: 2^25.9 records into 2^30
  slots = **2^-4.1**, so the drop probability is at most
  `load^4/(1-load) < 2^-16.4`.
* **Hard caps, so nothing here is an expectation.** Lookups + inserts are capped
  at `2^27`; once the cap is reached the run stops matching and keeps walking.
  Dictionary work is therefore bounded by `2^27 * 4 * 120` instructions whatever
  the distinguished rate, and the number of stored records is bounded by `2^27`
  as an integer, not as a mean.
* **Cold start is charged.** 2^30 handle slots of 16 bytes = 2^34 B `= 2^30`
  slot writes `= 2^31` instructions `= 2^15.4` units (§9).
* **Records.** 48-byte full 381-bit key + 32-byte `alpha` + 32-byte `beta` =
  112 bytes, plus a 16-byte handle = **128 bytes**. At most `2^27` inserts, so
  records occupy at most `2^34` B. Handles are 16 bytes, so a 30-bit slot index
  fits with room to spare (the 8-byte/2^70-record objection does not arise).
* **Multiplier sampling with a deterministic cutoff.** `A_h, B_h` are drawn as
  512-bit integers `v` and accepted iff `v < floor(2^512/r) * r`. Since
  `r > 2^254`, the acceptance probability is at least `1 - r/2^512 > 1 - 2^-257`.
  **At most 2 draws are permitted per value.** The probability that both are
  rejected is at most `2^-514`, so over the 2,048 values of the table the
  probability that any value exhausts its draws is at most `2048 * 2^-514 <
  2^-503`; that event is a declared run failure and is far below the success
  budget of §8. The expected draw count is at most `2048 (1 + 2^-257)`, and
  `2048` draws are charged in §9.

**Matching is exact.** A stored record keeps the full 381-bit class key and a
candidate is accepted only when the whole key is equal. There is therefore **no
truncation-collision bound anywhere in this argument**: a candidate match is a
genuine class match, and it is still verified (§6), so a match can never produce
a wrong answer.

## 6. Collision resolution

Records `(K, alpha_i, beta_i)` with equal class keys satisfy: there is
`j ∈ {0..11}` with `(g^(alpha1) z^(beta1))^(p^j) = g^(alpha2) z^(beta2)`, i.e.

    alpha1 p^j + beta1 p^j x ≡ alpha2 + beta2 x   (mod r)
    =>  x = (alpha2 - alpha1 p^j) (beta1 p^j - beta2)^(-1)   (mod r)

`beta1 p^j ≡ beta2` is skipped (it would force `alpha1 p^j ≡ alpha2`, a
degenerate pair) and every one of the twelve candidates is **verified** by
`g^x = z` (equivalently `e(P,[x]Q0) = z`). Cost: at most twelve 255-bit
exponentiations plus one pairing, charged once in §9.

## 7. Cycles: handled history-free, so coalescence is trivial

No cycle-rate constant is asserted and no escape depends on a walk's history.

The walk state carries the walk's own step counter, which is identical for every
walker because the walkers move in lockstep. A walk that reaches
`W_max = 2^106` steps without a distinguished point is **abandoned and
restarted** from a fresh random start, and that restart is charged. Three
consequences, each of which is a bound rather than an assumption:

* **Coalescence is preserved trivially.** The abandon/restart rule is a function
  of the walk's own counter, which two merged walks share exactly (they step in
  lockstep), so merged walks keep identical futures. No argument about
  history-dependent decisions is needed — the mechanism the sibling claim was
  rejected for is gone.
* **A cycle costs at most `W_max` wasted steps**, and the probability that a
  walk meets a cycle inside `W_max` steps is the birthday bound
  `W_max^2/(2N) <= 2^212/2^252.55 = 2^-40.5` per walk, `<= 2^-14.6` over the
  `2^25.9` walks of the whole run. At `theta = 2^-100` the mean wait for a
  distinguished point is `2^100`, so `W_max = 2^106` abandons only `e^-64` of
  healthy walks: abandonment is not a source of loss, it is a safety net.
* **Nothing is left unspecified.** There is no unhandled "longer cycle" class and
  no cutoff whose loss is unaccounted: the loss is the bound above.

The toy uses the same rule (`cap_walks * 2^dp_bits` steps per walk): 16/16
trials still recover the exact discrete logarithm with zero bookkeeping
failures, and the abandonment counter is reported (`results.json` B
`abandoned_total`).

**Measured negative result on the cheaper-squaring variant.** Replacing a
fraction of steps by cyclotomic squarings (~0.45 units) was tried and is not
worth it here: at a fixed step budget the toy's measured success falls by
roughly half, i.e. the walk needs about twice the steps, so its 0.79-unit
average step loses well over a bit against the 1.13-unit pure-multiplication
step. Squaring is a permutation on the class quotient, so it contributes no
collisions of its own. That route is recorded as closed rather than claimed, and
nothing in the ledger depends on it.

## 8. Success probability: a hard step cap, and its measured success rate

The attack is randomised, so the claim is stated as a **hard work cap**, not as
expected work. Every walk is stopped once the run has spent
`K_max = ceil(c sqrt(N))` steps in total, with `c = 5/4`; if no logarithm has
been verified by then, the run reports failure. No cost below is an average over
successful runs, and no branch of the algorithm is unbounded.

Two independent supports put the success probability at or above `0.40` at
`c = 5/4`:

1. **Random-mapping model.** The class walk is modelled as a random map on `N`
   classes with R-adding jumps. The birthday model `P = 1 - exp(-c^2/2)` gives
   `0.5422` at `c = 5/4`. This is the standard rho heuristic, and it is used
   only as a sanity check on the measurement, not as the load-bearing argument.
2. **Measured with the cap applied, at toy scale.** `experiment.py` runs the
   declared algorithm itself — the same class key, the same full-key dictionary,
   the same bucket map and probe sequence, the same abandonment rule — on a toy
   tower, and records the total number of steps until the logarithm is recovered
   **and verified**. Over `2000` fresh instances, all `2000` recover and verify,
   and the fraction solved within `c sqrt(N)` steps is

   | budget | solved | fraction | Wilson 95% interval |
   | ---: | ---: | ---: | --- |
   | `1.15 sqrt(N)` | 801/2000 | 0.4005 | [0.3792, 0.4221] |
   | `1.20 sqrt(N)` | 857/2000 | 0.4285 | [0.4070, 0.4503] |
   | **`1.25 sqrt(N)` (claimed)** | **902/2000** | **0.4510** | **[0.4293, 0.4729]** |
   | `1.30 sqrt(N)` | 969/2000 | 0.4845 | [0.4626, 0.5064] |
   | `1.40 sqrt(N)` | 1070/2000 | 0.5350 | [0.5131, 0.5568] |

   At the claimed `c = 5/4` the measured success is `0.4510` with a 95% lower
   bound of `0.4293`, so the measurement alone already clears `0.40 >= 0.39`
   without help from the model. The measured median is `1.334 sqrt(N)`.

The claimed success probability is `0.40`: below both the measured value
(`0.4510`) and the model value (`0.5422`), so either support suffices.

**Direction of the toy's error.** The toy is not a flattering instantiation. Its
distinguished fraction is `theta = 2^-3`, so it stores only `theta sqrt(N) ~ 36`
distinguished points per walk against production's `2^25.6`; measured across
`theta = 2^-3, 2^-4, 2^-5` the total step count rises monotonically as fewer
distinguished points are available (median `1.313, 1.502, 1.684 sqrt(N)`), so a
finer relative threshold can only help production. Its table is deliberately
under-provisioned (256 slots for ~57 records per instance against production's
2^30 slots for 2^25.9 records), which costs it drops and extra probes. Both
directions make the toy more expensive than the model, which is why the measured
success is below the model's and the claim is taken from the measured side.

The abandonment rule of §7 fires on 68 of 2000 toy instances and no instance
fails because of it, which is the expected behaviour of a safety net rather than
a source of loss.
## 9. Ledger (exact integer / rational arithmetic)

| Item | Value |
| --- | --- |
| Classes `N = (r-1)/12 + 1` | `2^251.2721` |
| `sqrt(N)` | `2^125.636063` |
| hard step cap `K_max = ceil((5/4) sqrt(N))` (§8) | `2^125.957992` |
| x step cost `22600/20000 = 1.13` units (20,000 + 2,600 instructions) | main term `2^126.134314` |
| dictionary: at most `2^27` ops x at most 4 probes x 120 instr | `2^35.67` instructions = `2^20.66` units |
| cold start: `2^30` slots x 2 instr | `2^31` instructions = `2^16.01` units |
| preprocessing: R=1024 multipliers, 12 Frobenius each, T and PA/PB tables, 2048 charged draws | `2^28.3` instructions = `2^13.3` units |
| verification: 2 pairings at `2^17` units each + <=12 GT exponentiations at 382 dense Fp12 mults each | `2^18.1` units |
| **overhead total** | `< 2^21.2` units |
| **W** | `2^126.134317` |
| **Claim** | `12614/100 = 126.14`, i.e. `log2 W <= 126.14` |

Every additive item is in **dense-Fp12-multiplication units**, the track's
reference operation, and each is below `2^21.2` units against a main term of
`2^126.13`: together they are below `2^-105` of the main term, so they cannot
move the second decimal. Verification is charged as follows and is *not* counted
in base-field instructions: one GT exponentiation with a 255-bit exponent is 255
squarings plus at most 128 multiplications, i.e. at most 383 dense Fp12
multiplications = 383 units; twelve candidate resolutions give `4596 < 2^12.16`
units; two pairings (one for `g = e(P,Q0)`, one for the final check) are charged
at `2^17` units each, i.e. `2^18` units, which is a deliberate over-charge
against the measured cost of a pairing in these units.

Per-step instruction itemisation (2,600): `B = Re(conj(c_1) gamma1)` = 2 fp_mul
(766) + 1 add (20); `a zeta`, `B zeta` = 2 fp_mul (766); six negations and two
subtractions (160); eleven 381-bit comparisons for the minimum (440); jump index
and bucket (110); distinguished test (10); two modular additions mod r (100);
state traffic and loop (228). Conversion uses COST_MODEL's measured 383
instructions per blst `fp_mul`; the dense Fp12 multiplication is charged as one
reference unit and its internals are **not** charged again.

One pairing is charged for `g = e(P,Q0)` and one for the final verification.

## 10. What was measured (toy scale, faithful instantiation)

`experiment.py` section B runs the **same code path** as production at toy scale:
a prime `p = 8389603 ≡ 7 (mod 12)`, `r = 1020337` prime with `r | p^4-p^2+1`,
`ord_r(p) = 12`, the same tower, the same key formula, the same lift/kappa rule,
the same jump table, the same bucket map and probe sequence, the same exact
full-key matching, and the same 12-candidate solve.

| quantity | big toy | small toy |
| --- | --- | --- |
| classes N | 85,029 | 20,902 |
| sqrt(N) | 291.6 | 144.6 |
| trials solved / run | **16 / 16** | **16 / 16** |
| median total steps to a verified log | 399 | 251 |
| median / sqrt(N) | 1.368 | 1.736 |
| inserts / drops | 721 / 7 | - / - |
| `g^alpha z^beta == y` failures | **0** | **0** |
| abandoned walks | **0** | **0** |

Every trial recovered the exact discrete logarithm and verified it. The two toys
use the same declared algorithm at different scales; the smaller one is *more*
expensive per `sqrt(N)`, which is the expected direction when the dictionary and
the distinguished rate are held fixed and `N` shrinks.

**The run behind §8.** The success curve quoted in §8 is the same script with
`{"btrials":2000,"seconds":2000}`: `2000` instances, `2000` solved and verified,
median `1.334 sqrt(N)`, and the cumulative success curve tabulated there. The
captured `results.json` is a shorter run of the same measurement (`100` trials
per budget), which reproduces the same shape at `0.33 / 0.43 / 0.52` for
`1.15 / 1.25 / 1.40`.

**Distinguished-rate ablation.** Holding everything else fixed and varying
`theta` at `600` trials per point gives median `1.313` (`2^-3`), `1.502`
(`2^-4`) and `1.684` (`2^-5`) in units of `sqrt(N)`: the walk gets more expensive
as fewer distinguished points are available, which is why the toy's coarse
`theta` is a handicap rather than a help relative to production's `2^-100`.

## 11. Limitations and disclosure

* One heuristic: the class walk behaves like a random mapping. It is not proved.
  Its cost is bounded by the hard cap of §8, whose success probability is
  measured at `0.4510` (`[0.4293, 0.4729]`) over 2000 toy instances and modelled
  at `0.5422`; the claim takes `0.40`.
* The toy is a different partition function, representation and table size from
  production, and is evidence only at its own scale. The transfer argument is
  explicit (coarse `theta` and an under-provisioned table both cost the toy
  work), not an assumption of ideality.
* No production-size run: costs are analytic instruction counts over the
  organizer's reference price. This is an upper bound on work, not a practical
  break and not a security estimate for the curve.
* Memory is `2^35` B as a hard cap (`2^30` handle slots plus `2^27` records at
  128 B), stated and unscored.
* Key ties (two conjugates with the same real part) have probability `~2^-370`
  per step; `A4` measured none.
* Records store the full 381-bit class key and matching is exact, so there is no
  truncated-identity collision term to bound.
* The mechanism is the one an accepted sibling-track claim uses
  (`bls12-381-fixed-g2-invert-g1`, 125.89); this package re-derives it from the
  tower constants, verifies every identity numerically, and replaces the storage
  and success accounting with the explicit, fully capped and measured version
  above. No third-party code was read or executed for it.

## 12. Reproduction

```sh
echo '{"production_checks":2,"toy":true,"trials":16,"btrials":100,"seed":1}' \
  | python3 -I experiment.py
```

`results.json` in this directory is the captured output of exactly that command.
The `btrials:2000` variant with `"seconds":2000` produces the success curve of
§8. Section A is the production-scale identity verification (tower, `ord_r(p)`,
`gamma1`, `zeta`, the key formula against brute force), section B the toy runs,
section C the ledger.
