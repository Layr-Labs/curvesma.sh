# argument.md — derivation, itemisation, and scope for `gt-reportfilter-40pct-v13`

Every number here is emitted by `ledger.py`, `gen_cert.py` and `toy_census.py`; the
certificate recomputes each of them from the raw evidence JSONs.

## 1. Reduction (checked on the real curve, not asserted)

`P` fixed nonidentity in G1, `z` uniform in the order-`r` subgroup of GT,
`g = e(P,Q0)` a generator by non-degeneracy, `a = log_g z`, `Q = [a]Q0`, so
`e(P,Q) = g^a = z`. `part_a.py` instantiates the ate pairing over the Fp12 tower
in `tower.py`/`pairing.py` and checks 24 named booleans on BLS12-381 itself:
`pairing_nondegenerate`, `pairing_order_r`, `pairing_bilinear_g1`,
`pairing_bilinear_g2`, `planted_ePQ_equals_z`, `planted_wrong_Q_fails`,
`ord_r_p_is_twelve`, `p_pow6_is_minus_one`, `n_classes_exact`, `key_rule`,
`keys_distinct`, `kappa_ok`, `min_ident`, `cyclo_sqr_identity` and the selector
identities. All 24 pass; `part_a-out.json:part_a` is `true`.

## 2. Class structure and the key

`ord_r(p) = 12`, `p^6 ≡ -1 (mod r)`, `N = (r-1)/12 = 2^251.27`, and the twelve
Frobenius images of `y` form one class. `keys12(y)` returns the twelve values
`{A zeta^j} ∪ {B zeta^j}`, `j = 0..5`, computed from the `c1[0]` coefficient of the
Fp12 element; `part_a` checks (`key_rule`) that this set equals the twelve
conjugate keys, (`keys_distinct`) that the twelve are distinct, and (`min_ident`)
that the minimum equals `-max` of the six absolute values — the identity the
selector row is counted against.

Under `sigma`, the *set* is preserved and the *index* of its minimum shifts:
`keys12(sigma^j y)` is the same set, `kappa(sigma^j y) = kappa(y) - j (mod 12)`.
Consequently `mkey = min keys12(y)` is `sigma`-invariant (a class function) and
`kappa` is a section coordinate. Both come out of the same twelve-way scan.

## 3. Equivariance of the transition — an identity

Step: `t, dp, sw = routing.fields(mkey)`;
if `sw < 3`, `y' = y^2`, `(a,b) -> (2a, 2b)`;
else `m = TBL[t][kappa] = sigma^{-kappa}(M_t)`, `y' = y·m`,
`(a,b) -> (a + U[t][kappa], b)` with `U[t][k] = u_t p^{-k} mod r`.

Let `y_j = sigma^j y`. `mkey` is unchanged, so `t` and `sw` are unchanged.
`kappa_j = kappa - j`, so

    m_j = sigma^{-(kappa - j)}(M_t) = sigma^j m,
    y_j · m_j = sigma^j y · sigma^j m = sigma^j (y · m) = sigma^j y'.

The squaring channel is `(sigma^j y)^2 = sigma^j(y^2)`. Hence the successor class
is a function of the class, for every representative, with no canonicalisation and
no section. `toy_census.t2_equivariance` compares the successor **elements**
(not keys) for 64 elements × 12 shifts, 768 checks, 0 splits. The two controls
both split: the v12c class-keyed rule 462/768 and the v12b canonicalising rule
704/768, so the test is shown to detect the defect it is meant to detect. The
shipped channel mix in that run is 252 squarings and 516 multiplies, so both
branches are exercised.

`t3_coalescence` starts two walks at independent representatives of one class with
independent labels and compares class keys at every checkpoint: 0 splits in 384
comparisons, 64/64 pairs coalesced at the end, against 349/384 for the v12c
control. A synchronized start cannot produce that number, because the two starts
are different elements.

## 4. Why no absorbing cycle exists (the v12b trap cannot be built)

v12b's cycle needs the step to **adopt** `(y m)^{-1}` and then apply the same `m`
again. In v13 no branch applies `sigma` to a state, so `sigma^6` (the inversion) is
never applied to a state, and the returning step does not exist. The only state
2-cycle available is `m·m' = 1` within one class, i.e.
`u_t p^{-kappa} + u_{t'} p^{-kappa'} ≡ 0 (mod r)`; with independent table exponents
this is a `1/r`-scale coincidence, not a structural family.

`t1_cycles` censuses the whole fixture (9 420 states, 785 classes): shipped rule
**0** state 2-cycles and **0** class 2-cycles; cycle-length histogram
`{9:4, 16:12, 63:4, 312:1}`, on-cycle mass 792/9 420 = 8.4 % (a random-map-sized
tail, not a trap basin). Mutation controls on the same census: the v12b rule
reproduces `{2:28, 3:6, 4:2, 11:1}` with **28** two-cycles, and the v12c rule its
`{6:1, 121:1}`. The trap is therefore measured, not argued away.

## 5. Sigma billing — the row is deleted, and the deletion is paid for

v12b's second finding was that applying `sigma^k` to the full Fp12 state was not
charged anywhere; v12d then charged `SIGMA_APPLY = 14 028` on every step (nine Fp2
constant products = 36 Fp muls in the odd-`k` worst case, six conjugations, six
accumulations). v13 never transforms the state, so the per-step row is
**0 by construction**: the `sigma` powers are pre-applied to the twelve table
columns.

The build cost is charged: `12 × 2^14` entries, each one `SIGMA_BUILD = 14 028`
plus one exponent multiply (`381 × FP_MUL`), i.e.
`2^20.585` reference units (`ledger.aux:table_build`), which is `2^-105.17` of the
claim. The table is `12 × 2^14` GT elements = 113 MiB; memory is unscored and
stated. `gen_cert` asserts the aux row is at least the `12 × 2^14 × SIGMA_BUILD`
product, so the row cannot silently shrink.

## 6. Per-step rows, itemised

| row | derivation | instructions |
|---|---|---|
| `selector` | `SELECTOR_AB = 1 MUL + 1 ADD`; `SELECTOR_ZETA = 2 MUL` (zeta^2 forms free); `FOLD_ADDS = 9 ADD`; `FOLD_CMP = 11 CMP6` | 1 613 |
| `column_addr` | `index = 12·t + kappa`, one multiply-add | 16 |
| `table_fetch` | 576-byte entry = 72 eight-byte loads | 72 |
| `label_add` | `a += u mod r`, one 255-bit add with conditional subtract | 28 |
| `key_hash` / `dp_test` / `hist_slot` / `loop` | one 256-bit mix into three fields; DP word test; 48-byte record slot; addressing | 72 / 8 / 32 / 32 |
| `fold_mul` | the eight rows above | **1 873** |
| `fold_sqr` | the same minus the fetch and the add, plus `2 × LABEL_ADD` for the doubling | **1 813** |
| `cyclo_sqr` | Granger–Scott: 24 Fp muls + 122 Fp adds | **11 632** |
| `mul_step` / `sqr_step` | `20000 + fold` / `11632 + fold` | **21 873 / 13 445** |

`c_mul = 1.093650`, `c_sqr = 0.672250`. With `f = 287/785 = 0.365605` measured
(declared 3/8; the smaller value is carried because squaring is cheaper, so it is
the higher price): `c_raw = 0.939584`.

`W_stress` reprices every judgement-bounded fold row at 13/10 and leaves the
counted selector and the counted field rows alone: `c_stress = 0.967350`,
`W_raw = 2^125.714913`, `W_stress = 2^125.756929`, `delta = 0.029551` (+2.07 %).

## 7. Step budget, success, detection, solve

`S` is the least integer with `S^2 (1-f^2) >= 2 mu N ln(1/(1-0.392))` at
`mu = 11/10`: `S = 2^125.804819`. `Lambda = 0.4975803970`,
`1 - e^-Lambda = 0.3920000 >= 0.39`. Claim `= ceil4dp(W_stress + 2e-5)`
(`ledger.claim_floor`, exact integer arithmetic, both directions asserted).

Detection: the DP predicate is a function of `mkey`, i.e. of the class. Two trails
at a common class report the same class in the same step, so no
`E[1-(1-d)^L]`/`E[L]` substitution is used (v12b's third finding): the collision is
observed at the first DP on the shared suffix.

Solve: two records `(mkey, a1, b1)`, `(mkey, a2, b2)` with `y1 = sigma^j y2` give
`a1 + b1 x = p^j (a2 + b2 x)`, hence `x = (p^j a2 - a1)/(b1 - p^j b2)`, tried for
`j = 0..11` and **verified by exponentiation**. `j = 6` is the involution and reads
`x = -(a1+a2)/(b1+b2)`. `t5_solver_end_to_end` builds, for each of 12 planted
secrets, a genuine `sigma^j`-partner record with an **independent** label pair (so
`b2` is not `p^j b1` and the denominator does not degenerate): 144/144 solved, one
verified candidate each, 0 wrong-shift accepts, 0 zero denominators, `j = 6`
covered 12/12.

`t2b_pair_algebra` is the `+-`-pair test proper: two walks start at `y` and
`y^{-1}` with independent labels and stay exact inverses for 8 steps (256/256),
32/32 label-bit solves verify, and the sign-bit-flipped control verifies 0/32.

## 8. What the `+-` quotient does and does not buy

`sigma^6 = -1`, so `±` is a subgroup of the `sigma` action: identifying `y` with
`y^{-1}` does not merge two distinct Frobenius classes. The `mu6` quotient is the
`sigma` quotient, and the class count is `N` either way. There is **no** free
doubling of effective DP density from the `+-` identification, and this package
does not claim one. What the `+-`/label-bit treatment buys is that the involution
arrives as a *branch of the solve* (`j = 6`) instead of as a *canonicalisation step
in the walk*, which is precisely what deletes the v12b trap. That is the whole of
the design change, and it is a deletion plus a re-derived fold, not a new
mechanism.

## 9. Census scope

`toy_census.py` runs the shipped step, the shipped routing and the shipped solver.
6 400 trials per density, 4 trails, one shared store per trial, trails end at their
first DP, unconditional statistic, per-report invariant `g^a z^b = y` asserted,
every recovery verified by exponentiation, record keys re-derived at
`(trial, trail, step)` as a planted attribution check. Results and the mandatory
disclosure that the fixture success rate (0.5692 at 2^-5, 0.2792 at 2^-7) is below
the 0.39 floor and falls with depth — because the fixture is cycle-dominated,
`S/sqrt(N) = 2.3` against 0.282 at production — are in `public-note.md` and
`certificate.json:census`. The fixture measures the mechanism (every collision
solves and verifies; 0 invariant failures; 0 zero denominators over 21 453
records) and the density trend, not a production recovery rate.

## 10. Legality note

Claim 125.7570 against a live lowest head of 125.84 polled at dispatch. Nothing in
this package depends on the identity or the validity of any other submission; the
comparison is context only.
