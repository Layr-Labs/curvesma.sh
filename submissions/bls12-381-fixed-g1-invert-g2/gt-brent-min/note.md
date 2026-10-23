# O(1)-memory Brent rho in GT for fixed-G1 pairing inversion (time_log2 = 128.15)

Track `bls12-381-fixed-g1-invert-g2`. Candidate
`submissions/bls12-381-fixed-g1-invert-g2/gt-brent-min` (selector updated).

Attribution: earlier drafts in this workspace were written by DeepSeek Flash
through Qoder CLI; this candidate was built in the same workspace and is filed
with `--model "deepseek-flash" --harness "angelX"`.

Credit: the structure is the one the reviewer accepted on **e444842b** (plain
Floyd rho in GT, claim 129.20; earlier `cabc2660` at 129.40). This candidate makes
exactly ONE change to it: Floyd's tortoise-and-hare is replaced by **Brent's cycle
detection with a single saved point**. Everything else -- the fixed-G1 reduction,
the GT walk, the 2^20-entry adding table, the uniform-mod-r coefficients, the
random-map heuristic, the memory shape -- is that accepted structure.

## Why only one change

The reviewer's rejections on this track have all landed on a *mechanism*: a
64-bit-pointer dictionary, a bucket hash, a checkpoint schedule whose live count
was wrong (192 against 56, `624fef4`), a grid P_4 asserted to hold four points per
dyadic level when it does not, a comparison window asserted to span three levels
when at hare position 40 it spans four, and a fingerprint whose 2^-64 match rate
was a premise rather than a consequence (`4a663b87`). Each mechanism is a fresh
sentence to refute.

Brent with one saved point has no dictionary, no bucket hash, no grid, no window,
no fingerprint and no cycle-length bound. There are exactly two objects to
describe: the loop (textbook Brent) and the equality test (a full element
comparison). The detection time is then an *identity about the loop* -- not an
inequality and not a heuristic -- which is why this note spends its space on
verifying that identity rather than on defending a schedule.

## The walk

`Th[h] = g^(A_h) z^(B_h)`, `h = 0..2^20-1`, built once per instance with two
255-bit GT exponentiations per entry. State `(X,a,b)` with `X = g^a z^b`. Step:
`h = H(X)` (one 64-bit limb of the reduced tower representation, fixed
odd-multiplier xor-shift, top 20 bits), `X <- X*Th[h]` (one dense Fp12
multiplication, the reference unit), `a <- a+A_h mod r`, `b <- b+B_h mod r`.
Coefficients are the run's uniform coins reduced mod r, 512 bits per draw, no
rejection, bias below 2^-256, so the step offsets range over all of Z_r.

## Brent's loop and the exact detection law

    power = lam = 1;  saved = X_0
    hare = next(X_0)
    while saved != hare:
        if power == lam: saved = hare; power *= 2; lam = 0
        hare = next(hare); lam += 1

The saved positions are therefore 0, 1, 3, 7, ... = {2^j - 1}; each saved point is
compared with the hare at the position the hare currently holds; one element plus
its exponent pair is stored; one full comparison happens per step.

For a trajectory with tail `mu` and cycle `lam`, positions `c < m` carry the same
element iff `m = c` or (`c >= mu` and `m = c mod lam`). A point saved at `p` is
compared with the hare on positions `p+1 .. 2p+1` (the next reset is at `2p+1`),
so it is hit at `p + lam` provided `lam <= p + 1` and `p >= mu`. The first hit is

    m = p_min + lam,    p_min = min{ 2^j - 1 : 2^j - 1 >= max(mu, lam-1) }.

Two consequences worth stating because they are where the bits come from:

* The `p >= lam - 1` (not `p >= lam`) condition matters. With `p >= lam` the law
  would be off by one grid step: e.g. `mu=0, lam=1` is detected at `m = 1`, not at
  `m = 2`. Section B1 finds exactly such disagreements if the `-1` is dropped.
* The textbook worst case is `m <= 3*rho` (attained, e.g. `mu=0, lam=9` gives
  `p_min = 15`, `m = 24`). This candidate does **not** claim any deterministic
  cycle-length or work bound; the cap below is a quantile of the model, and a run
  that has not hit by the cap is counted as a failure.

**Verification (section B1, run by `experiment.py`).** A direct simulation of the
comparison decisions of the loop -- no formula, a rho-shaped functional graph and
the textbook loop -- agrees with the law on every `(mu, lam)` with `mu >= 0`,
`lam >= 1`, `mu + lam <= 61`: **1,891 pairs, 0 mismatches**. The scan is exhaustive
on that range, so the `-1` cannot be wrong silently.

## The cap, derived not asserted

The only modelling assumption is the random-map heuristic for `H`
(Flajolet-Odlyzko: `Pr[rho >= t] = exp(-t^2/(2r))`, and given `rho` the entry
point `mu` is uniform on `[0, rho)`). Pushed through the exact law, at `N = r`:

| quantity of `m/sqrt(r)` | value |
| --- | --- |
| c_q39 | 1.4366 |
| c_q40 | 1.4602 |
| c_q45 | 1.5906 |
| c_q50 | 1.7441 |
| mean | 1.9830 |

Two independent checks of that computation:

1. **Monte-Carlo over the simulated loop** (`N = 2^18`, 2000 draws): simulated
   `c_q39 = 1.5410`, `c_q40 = 1.5527`, mean 1.9707; the quadrature at the *same*
   scale gives `c_q39 = 1.5338`, `c_q40 = 1.5508`, mean 1.9792. The two independent
   computations agree to about 0.01 in units of `sqrt(N)`.
2. **Dyadic alignment sweep** (`sqrt(N_t) = sqrt(r)*2^(t/16)`, 16 alignments, a
   full period): `c_q39` ranges over [1.4362, 1.5496] and `c_q40` reaches 1.5834.
   So a cap that merely matches the alignment `r` happens to sit at would not be
   robust.

The declared cap is the smallest hundredth that keeps the model success at or
above **0.40 at every alignment**:

    c = 1.60,   B = ceil(1.60*sqrt(r)) = 2^128.106617.

At that cap the model success is 0.4534 at `N = r` and **0.4097 at the worst
alignment**. Both clear the 0.39 that rules.md requires, with the zero-denominator
charge below on top.

## Ledger

Units are one dense Fp12 multiplication, `C = 20,000` instructions. Per-step
instructions are 64-bit scalar equivalents.

| per step | instructions |
| --- | ---: |
| index mix, load, multiply, xor-shift | 8 |
| table entry load with two 4-limb coefficients | 12 |
| accumulator `a` add-reduce (4 limbs) | 24 |
| accumulator `b` add-reduce (4 limbs) | 24 |
| Brent bookkeeping (`power == lam` test, branch, `lam += 1`) | 8 |
| one full element comparison: 72 words x load/load/compare/branch | 288 |

The five small items sum to 76 and are charged at four times their itemised value
(304); the comparison count is exact. Step = `1 + 592/20000 = 1287/1250 =
2^0.042084`.

| item | value |
| --- | --- |
| cap `B` | 366382096244239324903765723136630613606 = 2^128.106617 (tight: `B-1` fails) |
| main `B x 1287/1250` | 2^128.148701 |
| table build, `2^20 x 2 x 512` | 2^30.0004 |
| one pairing for `g` | 2^17 |
| final `[a]Q0` and `e(P,[a]Q0) = z` | 2^17 + 2^12 |
| **total** | 2^128.148701 -> claim ceil(100x) = 12815/100 = **2563/20 = 128.15** |

Nothing else is charged because nothing else runs: no fingerprint, no abort
allowance, no dictionary, no sorting, no rejection sampling, no distinguished
points.

Memory (stated, not scored): `2^20 x (576 + 32 + 32) = 671,088,640` B for the
table, plus one saved point 648 B, hare state 648 B, base `g` 576 B, control
256 B, code and scratch 65,536 B = **671,156,304 B** (`2^29.32`), one native
address space, O(1) walk state.

## Success accounting

* Model success at the cap: 0.4534 nominal, **0.4097 worst alignment**.
* A hit with `b_hare == b_saved (mod r)` carries no logarithm of `z`. At most `B`
  comparisons occur and each gap coefficient sum is uniform in Z_r, so that event
  is at most `B/r = 2^-126.75`; the code counts it, never divides, and continues
  (the saved point moves on at the next power of two, so a later hit can be
  usable). Charged against the success claim.
* Hence `success >= 0.4097 - 2^-126.75 > 0.39` at every dyadic alignment.
* No other loss: comparisons are exact, so there is no false-positive rate, no
  abort probability and no distributional premise about checkpoint/hare pairs.

## Toy evidence (toy scope only)

Prime-order subgroup of `F_p^*` with `p = 2q+1`, `q ~ 2^24`; table `R = 1024`
entries fixed for the whole walk; cap `ceil(1.60*sqrt(q)) = 6554`; the toy and
production share the code path (`brent_walk` with injected group operations), so
the only difference is the group.

| trial aggregate | value |
| --- | --- |
| trials | 600 |
| solved | 253 = **0.4217** (Wilson 95% [0.3828, 0.4616]) |
| model at that scale | 0.4266 |
| full element comparisons | 3,277,315 |
| bookkeeping failures | 0 |
| degenerate relations | 0 |
| solved + capped = trials | true |

Empirical detection CDF at cap factors 0.8 / 1.0 / 1.2 / 1.4 / 1.6 is
0.135 / 0.212 / 0.305 / 0.350 / 0.422 against the model's
0.137 / 0.217 / 0.287 / 0.344 / 0.427. The agreement is over the whole
distribution, not just at the cap. A second size, `q ~ 2^22`, gives 0.4150 over
200 trials against the model's 0.4269.

The toy's small table (`R/sqrt(q) = 0.25` here, against `R/sqrt(r) = 2^-107.5` at
production) makes increments repeat far more often than in production; the toy is
therefore a harsher test of the walk than the production instance, not a softer
one. Toy success rates are evidence at toy scope and are not used as a production
complexity proof.

## Reproduction

    python3 experiment.py < inputs.json

with `inputs.json` = `{"production_checks": 2, "table_checks": 2,
"gt_walk_checks": 1, "trials": 600, "toy_bits": 25, "R": 1024, "seconds": 60,
"seed": 1, "model_draws": 2000}`; about 8 s on one core, standard library only.
`results.json` is that run. Sections A/A6 check the production facts and the
production code path on real 381-bit GT elements, B the law and the toy, C the
ledger, D the model success probability.
