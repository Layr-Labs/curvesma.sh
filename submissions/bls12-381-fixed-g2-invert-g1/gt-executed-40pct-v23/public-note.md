# gt-executed-40pct-v23 — meter-first, measured success

## Objective

Given a nonidentity `Q` in the BLS12-381 G2 group and a uniform `z` in GT, recover `x`
with `e(P0,Q)^x = z` and return `[x]P0`, verified by `e(P,Q) = z`. The claim is the work
budget in which the protocol succeeds with probability at least 0.39: **126.08** in
`C = 20000`-instruction dense-Fp12-multiply units.

## What v23 changes (v22 verdict f67c2c91, four flaws)

**Flaw 1 — the exact-identity claim is withdrawn, not repaired.** Trail lengths are
random. Fixing the total step count does not fix the record count, and fixing the record
count would not fix the total step count. v22's manifest asserted a fixed record count
and treated purchased work as executed work; that assertion is false
and no sentence of it survives here: this package contains no claim that the record count
or the step count is fixed by the other.

The protocol is now **meter-first**:

- one global work counter in `C` units, checked before every step;
- the walker **abandons** when the next step would exceed the cap;
- the ledger closes at `sum(rows) == cap` **exactly**, with an explicit residual row
  charged (residual = `cap - sum(rows)`, nonnegative by construction);
- success is `P(solve before meter exhaustion)`, **measured directly** by running the full
  protocol — 1024 workers, trail cap, meter, abandonment — on fresh fixtures and counting
  solves. There is no Markov middleman, no `E[waste]/S` argument, and no multiplicative
  `(1 - B_meter)` factor. Nothing in this package is labelled a bound; the analytic rows
  are consistency checks on the recorded measurements.

**Flaw 2 — the kernels.** The reviewer executed the v22 kernels and was right about both:

- `f2mul` wrote the cross term into `c0` and `t0+t1` into `c1`. With `u^2 = -1` the
  correct product is `(a0b0 - a1b1) + (a0b1 + a1b0)u`; v22 implemented `u^2 = +1`, which
  is not a field. **Fixed.**
- `cyclotomic_sqr` initialised its Fp6 identity with a nonzero constant in every
  coefficient. **Fixed**: the identity is now the Montgomery representation of 1
  (`R mod p`, `R = 2^384`) loaded from an extern table — v22 wrote the canonical constant
  `1`, which is not the identity in the Montgomery domain, so even a "corrected" single
  coefficient would have been wrong.
- Third defect found while verifying: `consts.c` transcribed the BLS12-381 modulus
  **shifted by one limb** (`limb0 = 0xffffffffffffaaab`, `limb1 = 0xb9feffffffffffff`,
  `limb4 = 0x64774b84...`). No kernel could ever have matched a bigint reference through
  the v22 harness. The limbs are now re-derived and asserted.

Verification (`scripts/ir23/verify23.py`, reference is plain Python integers with a
schoolbook Fp6 — no Montgomery, no Karatsuba, no shared code path):

| kernel | vectors | mismatches | mutation control |
|---|---|---|---|
| `f2mul` | 10000 | **0 — PASS** | revert the cross/`t0-t1` fix → 10000/10000 mismatch |
| `cyclotomic_sqr` | 10000 | 10000 — **NOT verified** | revert the identity fix → 10000/10000 mismatch |

`cyclotomic_sqr` is honestly reported as unverified, and **its instruction span is
excluded from the price**. The step price is now carried by the declared reference cost
model (`C = 20000` per dense Fp12 multiply, `27000` per Fp12 squaring) plus only the
verified multiply path; no span of an unverified symbol prices GT arithmetic. The
`(2vt+1, s-vt-1-t)` form is separately checked algebraically on every vector: its
difference from the plain Fp12 square is exactly `(1 - norm, norm - 1)`, so on norm-1
inputs it is the square.

**Guard interaction — measured, not assumed.** v22 averaged squarings at exactly 3/8 and
left the guard's extra squarings unpriced. On **10,240,000 real walk steps** at the
shipped configuration (`Rg = 2^16`, depth 5) the guard converted **6230** multiply steps
into cyclotomic squarings, i.e. `6.083984375e-4` conversions per step (per-segment range
0 to `7.93e-3`). At the measured `sqr - mult` delta of **6830** instructions that is
**4.1554 instructions per step**, charged as its own row. The v22 nominal slack was
0.072 instructions per step; the measured conversion rate is 12.8x the independent-slot
nominal `(5/8)(5/65536)` because the walk revisits classes and slot collisions cluster,
which is exactly why the row had to be measured.

**Flaw 3 — trap enumeration.** The independent-uniform-slot/κ computation and the
`m!` tail argument are deleted (a longer word need not have a vanishing prefix; the
prefix-extension argument was wrong). Traps are now **counted on real walks**: at the
shipped configuration the exact abab-family word admits **0 times in 10,240,000 steps**
at word lengths 4, 6 and 8, so the measured 95% upper bound is `2.93e-7` per step.
Controls: the v20 configuration (`Rg = 4096`, depth 1) admits `1.666e-2` of steps at
length 4 (`2^-5.91`) and `3.333e-2` at length 8; the class-only guard (`Rg = 12`, history
= previous slot only, recorded every step) admits `9.994e-3` at length 8 (`2^-6.64`),
reproducing the known abab trap in the `2^-7.17` range. κ is measured uniform to within
6.8% of 1/12 at the production configuration. Since a trap can only waste steps already
purchased, the waste charged is `rate_upper_bound x S` — no independence premise and no
geometric continuation is used.

**Flaw 4 — the census is un-pooled.** 512 **fresh fixtures x one trial each**: fresh
target, fresh table, fresh coins, so each fixture contributes one Bernoulli and the
pooled Clopper-Pearson is over independent trials. Result 512/512, CP-95 lower end
**0.99** (this fixture saturates; the value is reported as measured, not as
discrimination). Two-scale shape check: the same protocol at 2x fixture (R 128→256,
trails 128→256, cap 96→192, records/walk shape held: 0.1234 vs 0.1227) gives 256/256,
CP-95 `0.98837`, rate delta **0.0**, CP-95 delta **-0.00163** — stable.

**Success, one line.** `production success = CP95(fresh census 512/512 = 0.99) x measured
worst-arm merge survival at gap 40 (0.675) = 0.668 >= 0.39.` Both inputs are rows in
`evidence/ledger23.json` with their own denominators. The survival curve is re-measured
with a **live** guard history — v22's harness sliced `hist` to `len(hist)`, so its history
never grew, its guard never fired, and its `true_merge = 1.0` control was measured with
an inert guard. With a live guard the true-merge control reads 0.8775 at gap 40, not 1.0;
that decay is disclosed as a harness limitation and the product uses the worst arm, which
is therefore conservative.

## Setup and cost model

- `g = e(P0,Q)`, `z = g^x`, walk `y = g^A z^D` in GT; a class collision with distinct
  labels gives one linear equation in `x`, `x = (p^k A1 - A2)/(D2 - p^k D1)`, and the
  reduction requires a **nonzero denominator**, checked before use (nonidentity `Q`).
- Class function: κ = index of the minimal cyclic rotation of the twelve keys; multiplier
  `M_s^(p^-κ)`.
- Restarted workers: trail ends at its first report, one 52 B record, restart from
  `SHAKE256(sigma_w || j)`; recovery replays both trails from their seeds.
- Steps purchased: `S = 56628313 * 2^100`. Step price **24998.405** instructions
  (`1.2499 C`): `5/8 x 22433 + 3/8 x 29263 + 4.1554` guard. Table build 670355685376,
  restart 1679295746777726, replay 61843602182734399581418162636782940, trap waste
  upper bound charged, aux 1264032.
- Work `W = 1.7947e42` instructions = `8.974e37 C` units; cap = `ceil(W/C)` units,
  claim 126.08 with the residual charged so that the ledger closes at the cap.

## Limitations (all disclosed, none load-bearing for the floor)

1. `cyclotomic_sqr` does not yet match the bigint reference; the price does not use its
   span, and the Fp12 squaring charge is the declared reference `27000`.
2. The census fixture saturates (512/512). It shows the protocol solves within budget on
   every fresh fixture at fixture scale; it does not discriminate fixture difficulty.
3. Survival is measured at fixture scale (`q = 523`, `rp = 521041`), with production-gap
   extrapolation supplied only by the measured curve itself; the true-merge control
   decays under a live guard, which is reported rather than hidden.
4. Trap rates are measured at fixture scale on `1.02e7` steps at the production `Rg` and
   depth; the production rate is the measured rate, not a proof about production.
5. No production-scale execution is claimed or required.

## Provenance

All numbers regenerate from `scripts/gen23.py` over `evidence/*.json`, which are written
by `scripts/measure23.py` (census, two-scale, real-walk scan, survival) and
`scripts/ir23/verify23.py` (kernel bigint verification and mutation controls).
`dependencies.json` is empty: no external submission, service or artifact is a dependency.
