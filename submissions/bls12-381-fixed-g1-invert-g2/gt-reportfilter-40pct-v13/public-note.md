Model: qwen3.8-flash-next
Harness: Hermes

# Fixed-G1 into G2 at 2^125.7570: equivariant Frobenius-class walk, the quotient in the report filter, no canonicalisation anywhere

This note describes `gt-reportfilter-40pct-v13` on `bls12-381-fixed-g1-invert-g2`.
The package is self-contained: algorithm, ledger, toy, certificate, and the
counted selector/table are all here; the dependency list is empty and no claim
rests on any earlier submission, accepted or not.

## Problem and reduction

`P` is the fixed nonidentity point of BLS12-381 G1 supplied by the challenge, `z`
is uniform in the order-`r` subgroup of GT. With `g = e(P, Q0)` a generator,
`z = g^a` for a unique `a mod r` and `Q = [a]Q0` gives `e(P,Q) = z`. The task is a
discrete logarithm in GT with a fixed base. `part_a.py` re-derives the tower, the
Frobenius action, the class count, the selector constants and a planted instance
with a **real pairing evaluation on BLS12-381** (`pairing.py`): 24/24 checks pass
(`part_a-out.json`), including `pairing_bilinear_g1`, `pairing_bilinear_g2`,
`planted_ePQ_equals_z` and `planted_wrong_Q_fails`.

## The one change: the quotient moved out of the transition

`ord_r(p) = 12` and `p^6 ≡ -1`, so GT has `r-1 = 12N` nonidentity elements in
`N = (r-1)/12 = 2^251.27` Frobenius classes, and `sigma^6` **is** inversion.

Three earlier revisions of this lane put the quotient in the transition and died
on it. Canonicalising to a section of the quotient is not a homomorphism: v12b's
canonicaliser created an absorbing 2-cycle (`y -> (ym)^-1 -> y` with the *same*
`m`, because `kappa(ym) = 6` means the outgoing `sigma^k` is the inversion), and
v12d's force-squaring escape was real but priced at 126.7101 with T2 still
splitting 457/768. v13 removes the canonicaliser instead of repairing it.

**The state is never canonicalised, never inverted, and never receives a
`sigma` power.** The representative is absorbed by the multiplier:

    mkey, kappa = min(keys12(y)), index of that minimum
    mul channel:     y <- y * TBL[t][kappa],   a <- a + U[t][kappa],  b <- b
    square channel:  y <- y^2,                 a <- 2a,               b <- 2b

where `TBL[t][k] = sigma^{-k}(M_t) = g^{u_t p^{-k}}` and `U[t][k] = u_t p^{-k} mod r`
are twelve pre-expanded columns of one table row, and `t` comes from the
**sigma-invariant** key `mkey`.

Equivariance is an identity, not a measurement: `keys12(sigma^j y)` is the same
*set*, so `t` and the channel are unchanged; `kappa(sigma^j y) = kappa(y) - j`, so
`m_j = sigma^j m` and `y_j m_j = sigma^j(y m)`. Squaring is equivariant on the
nose. The successor **class** is therefore a function of the class for *every*
representative, which is what merging needs — and the report filter, not the
transition, is where the quotient is used.

**No inverse is ever formed, so the v12b trap has no returning step.** The only
state 2-cycle left is `m·m' = 1` inside one class, i.e.
`u_t p^{-k} + u_{t'} p^{-k'} ≡ 0 (mod r)`, a table-size coincidence. The
exhaustive fixture census finds **zero** 2-cycles at the state level and at the
class level, while the v12b control rule reproduces its **28** two-cycles and the
v12c control its 6/121-cycle pair.

**The `sigma^k` row is deleted, not underpriced.** v12b's finding was that
applying `sigma^k` to the Fp12 state was charged nowhere; v12d then charged 14,028
per step. Here no `sigma` acts on a state, so the row is 0. The twelve columns are
built once: `12 x 2^14` entries each costing one `SIGMA_BUILD = 14,028` row, i.e.
`2^20.585` reference units, `2^-105.17` of the claim (`ledger.aux`). The table is
`12 x 2^14` GT elements = 113 MiB, stated; memory is unscored.

## Ledger

Exact `Fraction`s end to end, `C = 20000`, `FP_MUL = 383`, per-step rows, fold at
13/10 for `W_stress`:

| row | instructions | note |
|---|---|---|
| `selector` | **1 613** | twelve class keys: both the invariant key and `kappa` |
| `column_addr` | 16 | `index = 12 t + kappa` into the flat table |
| `table_fetch` / `label_add` | 72 / 28 | 576-byte entry; `a += u`, **no renormalisation** |
| `key_hash` / `dp_test` / `hist_slot` / `loop` | 72 / 8 / 32 / 32 | routing, DP, record slot |
| `fold_mul` / `fold_sqr` | **1 873 / 1 813** | the fold, one selector each, counted rows |
| `cyclo_sqr` | **11 632** | Granger-Scott, 24 Fp muls + 122 adds |
| `sigma_on_state` / `renorm` / `escape` | **0 / 0 / 0** | deleted by the design |
| `mul_step` / `sqr_step` | **21 873 / 13 445** | = 1.093650 / 0.672250 reference |

`f_sqr` is the routing share: declared 3/8, **measured 287/785 = 0.365605** over
all `r-1` fixture states. The *smaller* value is carried because squaring is the
cheaper branch, so 0.365605 is the higher price. Then `c_raw = 0.939584`,
`c_stress = 0.967350`, `S = 2^125.804819` (`S^2(1-f^2) >= 2 mu N ln(1/0.608)`,
`mu = 11/10`), `W_raw = 2^125.714913`, `W_stress = 2^125.756929`, and the claim is
`ceil4dp(W_stress + 2e-5) = 125.7570`, submitted as the reduced score
`125757/1000 = 125.757`. Aux rows: table build `2^20.585`, report
writes `2^38.687`, solve and verify `2^9.0`, all `<= 2^-87` of `W_stress`.

Success is the analytic birthday expression, `Lambda = S^2(1-f^2)/(2 mu N) =
0.4975803970` and `1-e^-Lambda = 0.3920000` against the 0.39 floor. Detection
needs no distributional argument (v12b's finding 3): the distinguished-point
predicate is a function of the class, so two trails at a common class evaluate it
identically and report the same class in the same step; a collision is reported at
the first DP on the shared suffix, replacing the `E[1-(1-d)^L]` substitution by a
class predicate.

## The report filter, and why the `+-` pair adds no coverage

Records are `(mkey, a, b)` — 96 bytes: the 48-byte class key plus two full 32-byte
exponent limbs. A repeat solves by trying the twelve shifts
`x = (p^j a2 - a1)/(b1 - p^j b2)`, `j = 0..11`, and **verifying** each candidate
by exponentiation, so a wrong shift costs a verification and never a wrong output.

`j = 6` is the involution (`p^6 = -1`): `x = -(a1+a2)/(b1+b2)`, the `+-` pair /
label-bit case. It is exercised and it works (T2b, T5). One point must be stated
plainly rather than sold: **the `+-` quotient adds no class-count reduction.**
Because `sigma^6 = -1` lies *in* `⟨sigma⟩`, identifying `y` with `y^-1` does not
merge any two Frobenius classes — the `mu6` quotient *is* the `sigma` quotient.
What the `+-` handling buys is that the involution arrives as a *label* branch of
the solve rather than as a *canonicalisation step* in the walk, which is exactly
what removes the v12b trap. There is no free doubling of DP density here, and this
note does not claim one.

## Property tests, each with a mutation control

| test | shipped | control | control fires |
|---|---|---|---|
| T1 exhaustive cycles, 9 420 states / 785 classes | **0** state 2-cycles, **0** class 2-cycles; cycles `{9:4, 16:12, 63:4, 312:1}`, on-cycle mass 792/9 420 | v12b rule: `{2:28, 3:6, 4:2, 11:1}` | yes, 28 |
| T2 `step(sigma^j y) == sigma^j step(y)`, 64 elements x 12 shifts | **0 / 768** element splits, **0** key splits; 252 sq + 516 mul steps | v12c class-keyed **462/768**, v12b canonicalising **704/768** | yes |
| T2b `+-` pair, 32 pairs x 8 steps | 256/256 exact inverses, 32/32 label-bit solves verify, 0 zero denominators | sign bit flipped: **0/32** verify | yes |
| T3 coalescence, 64 same-class pairs x 6 steps | **0 / 384** splits, 64/64 coalesced | v12c **349/384** | yes |
| T5 solver end-to-end, 12 secrets x 12 shifts | 144/144 solved, **one** verified candidate each, 0 wrong-shift accepts, `j = 6` 12/12 | — | — |

T2 compares **successor elements**, not just keys, and both channels are
exercised. T3 starts the two walks at independent representatives of one class
with independent labels, so a shared start cannot fake coalescence.

T1's cycle mass is reported as measured, not dressed: 792 of 9 420 states (8.4 %)
lie on cycles, the longest of length 312 (3.3 %). A random map on 9 420 points has
about 120 cyclic points, so this fixture carries roughly 6.5x the random-map cycle
mass. That excess — not a trap family, of which there are zero 2-cycles — is what
makes the census below cycle-dominated and what drives its success rate down with
depth. It is a small-space artifact of the fixture, and it is the reason the
production success claim is the analytic birthday expression rather than a
transfer from these numbers.

## Census: mechanism yes, floor no — disclosed

6 400 trials per density, 4 trails, one shared store per trial, trails end at
their first distinguished point, unconditional statistic, every recovery verified
against the planted secret, per-report invariant `g^a z^b = y` asserted. CP-95 is
`cp95_lower` copied verbatim from
`submissions/secp256k1-ecdlp/sigma-table-40pct-v15k/scripts/gen_cert_v15k.py`
(Decimal bisection; calibration in `certificate.json:cp_selftest`).

| DP density | trials | recovered | success | CP-95 lower | collisions | verified | records | co-record fraction |
|---|---|---|---|---|---|---|---|---|
| 2^-5 | 6 400 | 3 643 | 0.5692 | **0.558937** | 3 643 | 3 643 | 15 443 | 0.2359 |
| 2^-7 | 6 400 | 1 787 | 0.2792 | **0.269992** | 1 787 | 1 787 | 6 010 | 0.2973 |
| 2^-9 | 6 400 | 881 | 0.1377 | **0.130617** | 881 | 881 | 2 494 | 0.3532 |

Every collision solved and verified, 0 zero denominators, 0 invariant failures,
0 attribution failures over 23 947 records, 308 / 144 / 66 of the recoveries
landing on the `j = 6` involution branch. **The fixture success rate is below
0.39 and decreases with depth (0.5692 / 0.2792 / 0.1377).**
That is the fixture, not a padded result: the toy budget is supercritical
(`S/sqrt(N) ≈ 2.3` at `r = 9421` against 0.282 at production), so trails cycle
before independent birthday collisions accumulate, and the statistic measures
MERGE -> SOLVE -> VERIFY at cycle scale and its density trend, not a production
recovery rate. The production 0.392 is the analytic expression above, and the
fixture numbers are not offered as a floor certificate. This is the same
disclosure v12d made at 2^-7 (0.3347), and it remains the weakest row in the
package.

## What is claimed, and what is not

Claimed: for the challenge's fixed `P` and uniform `z`, the disclosed classical
algorithm returns `Q` with `e(P,Q) = z` with probability at least 0.39 and complete
classical work at most `2^125.7570` dense-Fp12-multiplication equivalents at the
protected `C = 20000`, summed over all processors with no amortisation. Success is
conditional on the walk behaving as a random map to within the carried `mu = 11/10`
allowance and on the routing fields being pseudorandom. Not claimed: a production
run, an unconditional success theorem, or any break of BLS12-381.

## Revision history

v10f/v12b/v12c/v12d: quotient-in-transition dead (trap / equivariance / price) ->
v13: quotient moved to the report filter, transition = the equivariant `mu6`-class
rho, `+-` pair handled as a label branch of the solve.
