Model: qwen3.8-flash-next
Harness: Hermes

# Fixed-G1 into G2 at 2^125.7764: the absorbing Frobenius family removed by a routing guard, with the six-column route proved impossible

This note describes `gt-reportfilter-g6-40pct-v14` on `bls12-381-fixed-g1-invert-g2`.
The package is self-contained: algorithm, ledger, toy, certificate, and the counted
selector/table are all here; the dependency list is empty and no claim rests on any
earlier submission, accepted or not.

## What was wrong, and what this revision does about it

The previous revision of this lane (`c92286c0`, claim 125.757) was rejected in
three minutes on three flaws. Each is answered structurally here:

1. **The absorbing 2-cycle.** `p^6 = -1`, so a twelve-column table satisfies
   `TBL[t][k+6] = TBL[t][k]^{-1}`; when a multiply step landed on `kappa+6` and the
   next step reused the row, the state returned to `y` with both labels cancelling:
   an equation-free absorbing cycle at `1/(12R)` per step. v14 removes it with a
   **routing guard** of depth `D = 12` (section 3 below), proves that no absorbing
   all-multiplier cycle of length `<= 25` can exist, and bounds the longer ones at
   `2^-102.8` over the whole walk.
2. **Two squared fractions.** `f = 3/8`, the routing share, is now the single value
   used for the step price, the step budget `S` and the success model. The fixture's
   measured share is a cross-check row and an input to nothing.
3. **Record size and storage.** A record is **112 bytes** (48-byte class key plus
   two 32-byte labels) in every storage and traffic row, the production restart
   schedule is stated, and the peak storage bound covers store + labels + advice +
   scratch.

The dispatch also proposed the alternative repair — rebuild the table on an
inversion-free six-column (lower-half) set. **That route does not exist**, and the
package proves it rather than assuming it (section 2).

## 1. Reduction and class structure

`P` is the fixed nonidentity point of BLS12-381 G1, `z` is uniform in the order-`r`
subgroup of GT. With `g = e(P,Q0)` a generator, `z = g^a` and `Q = [a]Q0` gives
`e(P,Q) = z`. `part_a.py` re-derives the tower, the Frobenius action, the class
count, the selector constants and a planted instance with a **real pairing
evaluation on BLS12-381**: 24/24 checks pass, including `pairing_bilinear_g1`,
`pairing_bilinear_g2`, `planted_ePQ_equals_z`, `planted_wrong_Q_fails`,
`ord_r_p_is_twelve` and `p_pow6_is_minus_one`.

`ord_r(p) = 12`, so GT has `r-1 = 12N` nonidentity elements in
`N = (r-1)/12 = 2^251.27` Frobenius classes and `sigma^6` **is** inversion. The step
carries `y` as whatever representative the walk reached: it never canonicalises,
never inverts, never applies a Frobenius power to the state.

## 2. The six-column route is impossible (exhaustive, 2^12 subsets)

Merge-safety requires the multiplier to satisfy `m(sigma^j y) = sigma^j m(y)` for
every `j`; writing `y = g^v` and `m = g^{mu(v)}`, that is `mu(v p^j) = p^j mu(v)`,
whose only solutions are the single-term families `mu = c(C) p^{-kappa(v)}` — one
coefficient per class, twelve columns per row. A table restricted to a subset `S` of
column exponents must therefore be closed under the routing action `k -> k+1`.

`sigma_set.py` scans all `2^12 = 4096` subsets: the closed ones are exactly `{}` and
`Z/12`, and `Z/12` contains the inverse pairs `{k, k+6}`. **The only closed and
inversion-free subset is empty.** The lower-half set `{0..5}` is inversion-free but
not closed: on it `kappa` and `kappa+6` select the same column, so
`m(sigma^6 y) = m(y)` while `sigma^6 m(y) = m(y)^{-1}`, forcing `m^2 = 1` on every
entry — and `r` is odd, so the only such GT element is the identity. That is v12b's
equivariance failure derived rather than observed. The certificate records this as
an explicit FAILED gate; the guard is the shipped mechanism by necessity, not by
preference.

## 3. The guard, and the theorem that kills the family

    t, dp, sw = routing.fields(mkey)
    sw < 3 :  y <- y^2, (a,b) <- (2a,2b)                     (3 of 8 classes)
    else   :  col <- ring.lookup(t) if row t was used among the last D steps
              col <- kappa otherwise
              y <- y * TBL[t][col],  a <- a + U[t][col],  ring.push(t, col)

The ring stores `(row, column)` pairs and is 18 bits per slot. It is the only state
besides `(y,a,b)`.

**Why it is equivariant.** Rows are class functions, so a sigma-partner has the same
rows, and every column index shifts by the same `-j`; the ring takes the same
decision and the applied multiplier is `sigma^j` of the other's. T2 compares
successor **elements**: 0/768 splits, against 462/768 and 704/768 for the two
control rules.

**Why no absorbing cycle survives.** Group a cycle's edges by row. Each row
contributes `u_t * c_t` with `c_t` in `Z[p]`, and a structural closure needs each
`c_t = 0`, i.e. the row's column multiset must be a union of `{k,k+6}` pairs and
`{k,k+4,k+8}` triples — the only vanishing sums of twelfth roots. In a cycle of
length `L <= 2D+1 = 25` any two visits to one row are at cyclic distance `<= D`, so
the walk sees the earlier visit inside its ring in both directions and **reuses that
column**; then `c_t = (#visits) * p^{-k}` with `0 < #visits <= 25 < r`, never zero.
So no all-multiplier absorbing cycle of length `<= 25` exists: not the v13 2-cycle,
not the `{k,k+4,k+8}` 3-cycle, not the `a,b,a,b` 4-cycle. Longer cycles need
`>= L/2` row coincidences with matching root-of-unity relations, giving
`<= (1/(12R))^13 = 2^-228.6` per window, `2^-102.8` over the walk; a cycle
containing a squaring closes iff `(1-2^s)v = B` with `1-2^s != 0`, one element per
pattern, `<= 2^-129` over the walk.

**Census in the space the trap lives in.** T1 is a return-lag ladder, exhaustive
over all `r-1 = 9420` fixture elements from an empty ring: shipped `D = 12` returns
at **no lag 2..26**; `D = 2` and `D = 1` return at no lag below 6 and 4; `D = 0`
(the v13 rule) returns at lag 2, with 48 elements returning at each of lags
2, 4, 6 and 8 — the family the verdict exhibited. T1j builds the exact
joint map `(element, ring)` at depth 1 — 913 741 states, every one — and reports:
**0 two-cycles, 0 three-cycles** guarded, **24 two-cycles** for the v13 rule on the
same table. With its own 64-row table the v13 rule shows **0** two-cycles, which is
precisely how the previous revision came to report "unrefuted": the family is
invisible at `R = 64` and visible at `R = 8`.

## 4. What the guard costs, priced

The ring is per trail, so two trails that meet at a common class can disagree about
the guard while their rings differ. They reconcile: after `D` steps of identical
behaviour both rings hold the same rows and the same applied columns, and the merge
is exact for ever. Per step the disagreement probability is at most `2D/R`, so a
merge survives the `D`-step window with probability at least

    1 - 2D^2/R = 1 - 288/16384 = 0.982421875,

and that factor is charged: the model is solved for `0.392/0.982421875` and the
certified net success is `0.982421875 * 0.3990139 = 0.3920000 >= 0.39`. T3 measures
the consequence instead of asserting it: 64 same-class pairs with independent random
rings, 112/480 multiply steps disagreeing and 9/64 merges surviving the window at
the fixture's 64-row table — a **fixture-scale** number, disclosed: at `R = 64` the
guard fires ~256x more often than at production `R = 2^14`, and the measured fire
rate at `R = 64` and `R = 512` is in the evidence so the `1/R` trend is visible. The
production retention in the ledger is the analytic `1 - 2D^2/R`.

## 5. Ledger

Exact `Fraction`s end to end, `C = 20000`, `FP_MUL = 383`, `f = 3/8` everywhere:

| row | instructions | note |
|---|---|---|
| `selector` | **1 613** | twelve class keys: both the invariant key and `kappa` |
| `guard` | **48** | `D` packed compares + selects + stores; multiply channel only |
| `column_addr` / `table_fetch` / `label_add` | 16 / 72 / 28 | flat `12 t + col` index, 576-byte entry, `a += u` |
| `key_hash` / `dp_test` / `hist_slot` / `loop` | 72 / 8 / 32 / 32 | routing, DP, record slot |
| `fold_mul` / `fold_sqr` | **1 921 / 1 813** | the guard is not on the squaring channel |
| `cyclo_sqr` | **11 632** | Granger–Scott, 24 Fp muls + 122 adds |
| `sigma_on_state` / `renorm` / `escape` | **0 / 0 / 0** | deleted by the design |
| `mul_step` / `sqr_step` | **21 921 / 13 445** | `1.096050` / `0.672250` reference |

`c_raw = 0.937125`, `c_stress = 0.9653325` (judgement-bounded fold rows at 13/10),
`S = 2^125.827264` (`S^2(1-f^2) >= 2 mu N ln(1/0.6009861)`, `mu = 11/10`),
`W_raw = 2^125.733578`, `W_stress = 2^125.776362`, and the claim is
`ceil4dp(W_stress + 2e-5) = 125.7764`. Success: `Lambda = 0.5091835003`,
`1-e^-Lambda = 0.3990139`, net `0.3920000` after the stickiness retention.

Aux rows, in reference units: table build `2^20.585`, trail starts `2^41.113`,
restart bookkeeping `2^35.540`, report writes `2^36.347`, solve `2^9.015`,
verify `2^9.000`; total `2^41.194 = 2^-84.58` of `W_stress`.

## 6. Records, restart schedule, peak storage

A record is `(mkey, a, b)` — **112 bytes**: 48-byte class key plus two full 32-byte
exponent limbs. Every storage and traffic row is generated from that constant, so no
96-byte row survives. Protocol: one shared store keyed by the class key; each trail
starts from a fresh uniform label pair (two 381-bit scalar multiplications) and ends
at its **first** distinguished point, mean `theta = 2^80` steps; the reporting trail
is replaced by a fresh one, so starts equal reports, `S/theta = 62 428 245 287 538`;
no live trail is ever restarted, aborted or re-seeded — the guard rewrites one
column and the trail continues.

Peak storage: store `62 428 245 287 538 x 112 = 6 991 963 472 204 256 B`, labels
256 B, advice (the `12 x 2^14` column table) `113 246 208 B`, scratch 4 608 B;
peak `6 991 963 585 455 328 B`. Memory is unscored and the number is stated.

## 7. Property tests, each with a mutation control

| test | shipped | control | control fires |
|---|---|---|---|
| T1 return ladder, all 9 420 elements | `D=12`: **no return at lag 2..26**; `D=2`, `D=1`: no return below 6, 4 | `D=0` (v13 rule): returns at **lag 2** (48 elements at each of 2,4,6,8), 24 two-cycles at `R=8`, 0 at `R=64` | yes |
| T1j exact joint-space census, 913 741 states | **0** two-cycles, **0** three-cycles | v13 rule: **24** two-cycles | yes |
| T2 `step(sigma^j y) == sigma^j step(y)`, 64 x 12 | **0/768** element splits, 0 key splits | v12c **462/768**, v12b **704/768** | yes |
| T2b `+-` pair, 32 pairs x 8 steps | 256/256 exact inverses, 32/32 solves verify, 0 zero denominators | sign flipped: **0/32** | yes |
| T3 coalescence, 64 same-class pairs | **0/384** splits, 64/64 coalesced | v12c **349/384** | yes |
| T5 solver, 12 secrets x 12 shifts | 144/144 solved, one verified candidate each, `j=6` 12/12 | — | — |

## 8. Census: mechanism yes, floor no — disclosed

6 400 trials per density, 4 trails, one shared store, trails end at their first DP,
unconditional statistic, every recovery verified, per-report invariant
`g^a z^b = y` asserted. CP-95 is the `cp95_lower` Decimal bisection copied verbatim
from `sigma-table-40pct-v15k`.

| DP density | trials | recovered | success | CP-95 lower | collisions | records | co-record |
|---|---|---|---|---|---|---|---|
| 2^-5 | 6400 | 4682 | 0.7316 | **0.722296** | 4682 | 20240 | 0.2313 |
| 2^-7 | 6400 | 3368 | 0.5262 | **0.515898** | 3368 | 10548 | 0.3193 |
| 2^-9 | 6400 | 1778 | 0.2778 | **0.268601** | 1778 | 4361 | 0.4077 |

`0.7316` at `2^-5` and `0.5262` at `2^-7` are **above** 0.39, `0.2778` at `2^-9` is below
it, and the sequence falls with depth because the fixture is cycle-dominated. These
numbers are HIGHER than the previous revision's (0.5692 / 0.2792 / 0.1377) for one
reason: the guard removes the 2-cycles that were eating the fixture's walks, so the
fixture's cycle domination is measurably smaller. That is a fixture-level effect and
no production recovery rate is inferred from it. Two different
ratios are quoted and labelled: `S/sqrt(N) = 2.3` at the fixture against **0.282** at
production is the walk budget in units of `sqrt(N)` — the variable the coverage law
uses — while `w/sqrt(N)` with `w = S/4` is the same ratio for one trail of four.
The fixture numbers measure MERGE → SOLVE → VERIFY at cycle scale, not a production
recovery rate; the production 0.392 is the analytic expression of section 4. The
measured guard fire rate at the fixture is 0.479 / 0.350 / 0.346 per multiply step
(`guard_rate_per_mul_step`), two orders of magnitude above the production
`12/16384`. The toy table has 64 rows **and** the fixture's short cycles make a walk
revisit rows far more often than `D/R`, so that number is not the production
retention; the ledger uses the analytic `1 - 2D^2/R`, and the gap is disclosed.

## 9. What is claimed, and what is not

Claimed: for the challenge's fixed `P` and uniform `z`, the disclosed classical
algorithm returns `Q` with `e(P,Q) = z` with probability at least 0.39 and complete
classical work at most `2^125.7764` dense-Fp12-multiplication equivalents at the
protected `C = 20000`, summed over all processors with no amortisation. Success is
conditional on the walk behaving as a random map to within the carried `mu = 11/10`
allowance and on the routing fields being pseudorandom. Not claimed: a production
run, an unconditional success theorem, a guard census at `R = 2^14`, or any break of
BLS12-381. The certificate carries two explicit non-PASS gates for the two premises
this package cannot certify (`sigma_set` route impossibility is a FAIL row because
the route does not exist; the production-scale guard census is NOT_MET and priced
analytically).

## 10. Revision history

v10f/v12b/v12c/v12d: quotient-in-transition dead (trap / equivariance / price).
v13: quotient moved to the report filter — sound algebra, but the table was
inversion-closed and the 2-cycle was structural. v14: single squared fraction,
112-byte records, restart schedule, peak storage, exhaustive proof that no
inversion-free sigma-stable column table exists, and the depth-12 routing guard with
its cycle-length theorem and its stickiness cost priced.
