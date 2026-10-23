# Proof: O(1)-memory Brent rho in GT for fixed-G1 pairing inversion

Candidate `submissions/bls12-381-fixed-g1-invert-g2/gt-brent-min`.
Claim: total work `log2 <= 12815/100 = 128.15` in reference dense-Fp12-multiplier
equivalents (`C = 20,000`), with success `>= 39/100` over uniform full-production
challenges and fresh algorithm coins.

The structure is the one accepted on `e444842b` (plain Floyd rho in GT, 129.20).
This candidate changes exactly one thing: the cycle detection is Brent's, with a
single saved point, instead of Floyd's tortoise-and-hare. The detection time is
then an identity about the loop, verified by direct simulation, rather than an
inequality or a schedule. No fingerprint, checkpoint grid, comparison window,
distinguished-point set, dictionary, cycle-length bound or positivity lemma is
present.

---

## 1. The problem and the reduction

The instance supplies a nonidentity `P` in `G1` and `z` in `GT`; the task is `Q` in
`G2` with `e(P,Q) = z` for the pinned public pairing.

1. Compute `g = e(P,Q0) = pairing(Q0,P)` once. The pairing is nondegenerate and
   `r` is prime, `GT = mu_r`, so `g` generates `GT` whenever `P` is nonidentity.
   (Section A1/A2/A3 check on BLS12-381 itself: `p` and `r` prime,
   `r | Phi_12(p)`, `ord_r(p) = 12` with `p^6 = -1`, the blst tower
   `Fp2 = Fp[u]/(u^2+1)`, `Fp12 = Fp2[w]/(w^6-(1+u))` is a field, and a random
   `t^((p^12-1)/r)` is a nonidentity element of order `r`.)
2. Find `a = log_g z` by the rho of section 2.
3. Return `Q = [a]Q0`. Then `e(P,Q) = e(P,Q0)^a = g^a = z`. For `z = 1` this is
   `a = 0`, `Q = O` (the identity of `G2`), which is a valid answer.

The collision-to-logarithm identity, checked on real 381-bit GT elements
(section A4, 2/2): write `X_i = g^(a_i) z^(b_i)`. If `X_1 = X_2` with
`b_1 != b_2 (mod r)`, then `g^(a_1-a_2) = z^(b_2-b_1) = g^(x(b_2-b_1))`, so
`x = (a_1-a_2)(b_2-b_1)^(-1) mod r`. If `b_1 = b_2` the relation carries nothing
and is counted, never divided (section 6).

## 2. The walk

A public table `Th[h] = g^(A_h) z^(B_h)` for `h = 0..2^20-1` is built once per
instance. Each entry costs two 255-bit fixed-exponent GT exponentiations.

State `(X, a, b)` with `X = g^a z^b`. One step:

    h = H(X)
    X <- X * Th[h]          (one dense Fp12 multiplication = 1 reference unit)
    a <- a + A_h mod r
    b <- b + B_h mod r

`H(X)` reads ONE 64-bit limb of the reduced tower representation and applies a
fixed odd-multiplier xor-shift, keeping the top 20 bits. It is a fixed positional
read: no canonical-form conversion, no encoding export, no rejection, no
derivation from the position. The next state is therefore a function of the group
element alone -- the walk is memoryless, which is exactly the hypothesis the
random-map heuristic needs.

`A_h, B_h` are the run's uniform coins reduced modulo `r`, 512 bits per draw, no
rejection sampling; the statistical distance from uniform on `Z_r` is below
`2^-257`. The step offsets range over all of `Z_r`: there is no small-coefficient
window and no lemma about one.

## 3. Brent's cycle detection and its exact law

### 3.1 The loop

    power = lam = 1
    saved = X_0                      # position 0
    hare  = X_1                      # position 1
    while saved != hare:
        if power == lam:
            saved = hare             # save the hare's position
            power *= 2
            lam = 0
        hare = next(hare)
        lam += 1

One group multiplication and ONE full element comparison per step. The comparison
is a 576-byte equality (72 eight-byte words), so a comparison cannot be a false
positive and no work is ever aborted.

**Lemma 1 (saved positions).** The saved positions are `p = 2^j - 1` for
`j = 0, 1, 2, ...`, and after a save at `p` the loop saves next at `2p + 1`.

*Proof.* At the start the saved position is `0 = 2^0 - 1` and `power = lam = 1`.
Suppose a save happens at position `p` with `power = lam = P`. Immediately after
the save `power = 2P`, `lam = 0`, and `lam` increments once per subsequent hare
advance, so at hare position `p + k` we have `lam = k`. The next save occurs at the
first position where `lam = power = 2P`, i.e. at `p + 2P`. By induction
`P = p + 1` (base: `p = 0`, `P = 1`), so the next save is at
`p + 2(p+1) = 2p+1`; and `2(2^j-1)+1 = 2^(j+1)-1`. `[]`

Note the convention: the comparison is evaluated *before* the save test, so the
position at which the hare is saved is still compared against the previous saved
point.

### 3.2 The comparison window

**Lemma 2 (window).** A point saved at `p` is compared with the hare exactly at
positions `p+1, p+2, ..., 2p+1`.

*Proof.* Immediate from Lemma 1 and the loop order: after the save at `p` the hare
advances from `p` and is compared at each following position; the next save is at
`2p+1`, and the comparison at `2p+1` still uses the point saved at `p`. `[]`

### 3.3 The exact detection law

Let the trajectory have tail `mu` and cycle `lam` (rho shape: positions
`c < mu` are tail, positions `c >= mu` lie on a cycle of length `lam`). Then
positions `c < m` carry the same element iff `m = c`, or `c >= mu` and
`m = c mod lam`.

**Theorem (first hit).** The first position `m > 0` at which the hare equals the
saved point is

    m = p_min + lam,     p_min = min{ 2^j - 1 : 2^j - 1 >= max(mu, lam-1) }.

*Proof.* Two directions.

*(a) The hit occurs no later.* Take `p = p_min`. Then `p >= mu`, so the saved
element is on the cycle, and `p >= lam - 1`, i.e. `lam <= p + 1`, so the position
`p + lam` lies in the window `[p+1, 2p+1]` of Lemma 2. Since `p` and `p+lam` are
congruent modulo `lam` and both on the cycle, they carry the same element, so the
comparison at `p + lam` succeeds.

*(b) No earlier comparison succeeds.* Consider any earlier position `m` and the
point saved at that time, at some `p' < p_min`. If `p' < mu` the saved element is
in the tail, and a tail position is never revisited (the trajectory leaves it at
`p'` and never returns), so no later comparison can match it. If `p' >= mu` then
`p' < p_min` forces `p' < lam - 1`, so the window ends at `2p' + 1 < p' + lam`:
the only comparison at a position congruent to `p'` modulo `lam` is at `p' + lam`
itself (the next one, `p' + 2*lam`, is larger still), and it lies beyond the
window. Hence that saved point is dropped before it can be hit.

Combining (a) and (b), the first hit is at `p_min + lam`. `[]`

The `lam - 1` (not `lam`) condition is the whole content of the `-1` in the grid:
it is what makes small cycles detectable at all. For `mu = 0, lam = 1`
(fixed point) the law gives `p_min = 0`, `m = 1`; with `p >= lam` it would give
`m = 2`, which is false. The theorem also reproduces the textbook worst case: for
`mu = 0, lam = 9` it gives `p_min = 15` and `m = 24`, so the naive
`2*lam + 3 = 21` bound is false (as the review of `624fef4` pointed out). This
candidate claims **no** deterministic bound; `m <= 3*rho` is mentioned only to say
that the cap below is not derived from it.

### 3.4 Verification against the loop (section B1)

`brent_law_simulate` builds the rho-shaped functional graph for `(mu, lam)`, runs
the textbook loop of 3.1 with no formula, and reports the position of the hare at
the first successful comparison. `check_law` compares it with the theorem on every
`(mu, lam)` with `mu >= 0`, `lam >= 1`, `mu + lam <= 61`:

| pairs scanned | mismatches |
| ---: | ---: |
| 1,891 | **0** |

The scan is exhaustive on that range, so the `-1` in `p_min`, the window end
`2p+1`, and the tail/cycle distinction are all exercised, not assumed.

## 4. The hard cap, derived from the detection law

The only modelling assumption is the random-map heuristic for `H`
(Flajolet-Odlyzko): the trajectory has
`Pr[rho >= t] = exp(-t^2/(2N))` with `rho = mu + lam` and, given `rho`, the entry
point `mu` is uniform on `[0, rho)`. This is the heuristic the review has accepted
on this track, and it is the only one used.

Pushing that model through the exact law of 3.3 by deterministic quadrature
(section D1, 600 x 200 cells on `u = rho/sqrt(N)`, carried mass 0.999981) gives,
at `N = r`:

| quantity of `m/sqrt(r)` | value |
| --- | ---: |
| `c_q39` | 1.4366 |
| `c_q40` | 1.4602 |
| `c_q45` | 1.5906 |
| `c_q50` | 1.7441 |
| mean | 1.9830 |

Two independent checks:

* **D2, Monte-Carlo over the simulated loop.** At `N = 2^18` with 2000 draws of
  `(rho, mu)` and `m` taken from `brent_law_simulate` (not from the theorem):
  `c_q39 = 1.5410`, `c_q40 = 1.5527`, mean 1.9707. The quadrature at the same
  scale gives `c_q39 = 1.5338`, `c_q40 = 1.5508`, mean 1.9792. Agreement to about
  `0.01 sqrt(N)`.
* **D3, dyadic alignment.** `sqrt(N_t) = sqrt(r) * 2^(t/16)` for `t = 0..15`
  sweeps every fractional part of `log2 sqrt(N)` at spacing `1/16`. `c_q39` ranges
  over `[1.4362, 1.5496]` and `c_q40` peaks at 1.5834, so a cap fitted only to the
  alignment `r` happens to occupy would not be robust.

The declared cap is the smallest hundredth that keeps the model success at or
above 0.40 at **every** alignment:

    c = 160/100 = 1.60,    B = ceil(c * sqrt(r)).

`B` is computed exactly as the least integer with `(100 B)^2 >= 160^2 r`:
`B = 366382096244239324903765723136630613606 = 2^128.106617`, and `B - 1` fails
the inequality (`cap_tight` in the ledger). A run takes at most `B` steps; a run
that has not hit by then is censored and counted as a failure of the declared
bound.

At the cap, the model success is **0.4534** at `N = r` and **0.4097** at the worst
alignment.

## 5. The ledger

Units are one dense Fp12 multiplication, `C = 20,000` instructions
(`tracks/docs/COST_MODEL.md`, blst `fp12_mul` at 22,996 measured instructions).
Per-step instruction counts are 64-bit scalar equivalents, as the cost model
requires (arithmetic, loads/stores, addressing and control together).

| per step | instructions |
| --- | ---: |
| index mix, load, multiply, xor-shift | 8 |
| table entry load, element plus two 4-limb coefficients | 12 |
| accumulator `a` add-reduce, 4 limbs | 24 |
| accumulator `b` add-reduce, 4 limbs | 24 |
| Brent bookkeeping: `power == lam` test, branch, `lam += 1` | 8 |
| one full element comparison, 72 words x load/load/compare/branch | 288 |

The five small items sum to 76 and are charged at four times their itemised value
(304); the comparison count is exact, and there is exactly one comparison per step
because Brent keeps one saved point. So

    step = (20000 + 304 + 288) / 20000 = 1287/1250 = 2^0.042084 reference units.

Main term: `B * 1287/1250 = 2^128.148701`.

Additive items, all charged:

| item | units |
| --- | ---: |
| table build: `2^20` entries x 2 exponentiations x 512 dense multiplications | 2^30.0004 |
| one pairing for `g = e(P,Q0)` | 2^17 |
| final `[a]Q0` and `e(P,[a]Q0) = z` | 2^17 + 2^12 |

Total `2^128.148701`; the claim is its ceiling at the next hundredth,
`12815/100 = 2563/20 = 128.15`.

Nothing else is charged because nothing else runs: no fingerprint, no abort
allowance, no distinguished-point set, no dictionary, no sorting, no rejection
sampling, no checkpoint structure, no second pass. The final verification is
charged even though it is not needed for the algorithm to be correct.

## 6. Success accounting

Success is over uniform `P`, uniform `z` and fresh algorithm coins, as rules.md
requires.

* The walk hits within the cap unless the trajectory's first hit exceeds `B`, i.e.
  with probability `1 - P(m > B)`, which the model above prices at >= 0.4097 for
  every alignment.
* **Zero-denominator charge.** A hit with `b_hare = b_saved (mod r)` carries no
  logarithm: it is counted (`degenerate`), never divided, and the walk continues,
  since the saved point moves on at the next power of two. The run makes exactly
  `B` comparisons and each gap coefficient sum is uniform in `Z_r`, so
  `Pr[a usable relation is lost this way] <= B/r = 2^-126.75`.
* Every comparison is exact (a full 576-byte equality), so there is no
  false-positive rate, no abort probability and no premise about the joint
  distribution of checkpoint/hare pairs. Nothing else touches the success claim.

Hence

    success >= 0.4097 - 2^-126.75 > 0.39     at every dyadic alignment.

The code path for the unusable branch is exercised directly (section B3): a
constructed pair with an identical exponent pair is recognised as unusable, and an
injected identity table entry makes the walk take the branch 64 times in a
64-step run, count each one, and never divide.

## 7. Toy evidence and code-path identity

The toy is a prime-order subgroup of `F_p^*` with `p = 2q+1` and `q ~ 2^24`; the
table has `R = 1024` entries fixed for the whole walk; coefficients are uniform
mod `q`; the cap is `ceil(1.60 sqrt(q)) = 6554`. The toy and the production
instantiation share `brent_walk`: the same function, with the group operations
injected (`toy_ops` vs `gt_ops`), so the loop, the single saved point, the
one-comparison-per-step rule and the full-element equality are literally the same
code. Section A6 runs that function on real 381-bit GT elements and reports
`full_compares == steps`, `capped` at the cap, and the saved-point identity
recomputed from the saved exponent pair.

600 trials at `q ~ 2^24`:

| aggregate | value |
| --- | --- |
| solved | 253 / 600 = 0.4217, Wilson 95% [0.3828, 0.4616] |
| model at that scale | 0.4266 |
| full element comparisons | 3,277,315 |
| bookkeeping failures | 0 |
| degenerate relations | 0 |
| solved + capped = trials | true |

Detection CDF at cap factors 0.8 / 1.0 / 1.2 / 1.4 / 1.6, empirical vs model:
0.135/0.137, 0.212/0.217, 0.305/0.287, 0.350/0.344, 0.422/0.427. The agreement is
over the whole distribution, not just at the cap. A second size, `q ~ 2^22`, gives
0.4150 over 200 trials against the model's 0.4269.

Two honest caveats. (i) The toy's table is small relative to `sqrt(q)`
(`R/sqrt(q) = 0.25` here against `R/sqrt(r) = 2^-107.5` at production), so its
increments repeat far more often than production's and the toy is a harsher test
of the walk than the production instance, not a softer one. (ii) The toy's coins
are a BLAKE2b counter surrogate and make no statistical-distance claim; the toy
numbers are evidence at toy scope only and are not used as a production complexity
proof. The production claim rests on the theorem of 3.3 and the model quantiles of
section 4.

## 8. Memory

| object | bytes |
| --- | ---: |
| jump table `2^20` x (576-byte GT element + 32-byte `A_h` + 32-byte `B_h`) | 671,088,640 |
| one saved point: element 576 + `a` 32 + `b` 32 + position 8 | 648 |
| hare state: element 576 + `a` 32 + `b` 32 + step 8 | 648 |
| base `g` | 576 |
| control: `power`, `lam`, counters, RNG state, pointers | 256 |
| code, stack and scratch | 65,536 |
| **total** | **671,156,304** = 2^29.32 |

One native address space, `O(1)` walk state, one stored point at peak, no
segmented or distributed memory and no external memory traffic beyond the table
entry each step (charged in the per-step ledger).

## 9. Disclosures and limitations

* The random-map model is the only heuristic, and it is the one the review has
  accepted on this track. The detection time itself is not heuristic: it is the
  theorem of 3.3, verified exhaustively against the loop on the scanned range.
* No deterministic cycle-length or work bound is claimed. The cap is a quantile of
  the model; a run that exceeds it is a declared failure, not a retry.
* Toy evidence is toy-scope evidence, with the caveats of section 7.
* The claim is a bound on total work, not a demonstration of a production-scale
  run; rules.md explicitly does not require one.
* The `2^20`-entry table is built per instance and charged in full; no advice is
  reused across instances and nothing is amortised.

## 10. Reproduction

    python3 experiment.py < inputs.json

`inputs.json`:

    {"production_checks": 2, "table_checks": 2, "gt_walk_checks": 1,
     "trials": 600, "toy_bits": 25, "R": 1024, "seconds": 60,
     "seed": 1, "model_draws": 2000}

Standard library only, about 8 s on one core. `results.json` is exactly that run.
Section A: production facts on BLS12-381. A6: the walk function on real GT
elements. B1: the law against the loop. B3: the toy, with the zero-denominator
branch. B4: a second toy size. C: the ledger. D: the model success probability,
its Monte-Carlo cross-check and the alignment sweep.
