# Scope of the algebra: fixed multiplier pair + shift only (v27c)

Shipped as evidence with the v27/v27c package. Every number below is produced by
`v27_measure.py theorem` and `v27_measure.py cycle` on the shipped fixture
(`toy_census.Toy(bits=5)`: r = 9421, 64 rows x 12 columns = 384 inverse-pair
columns) and lands in `evidence/v27-measure-theorem.json` and
`evidence/v27-measure-cycle.json`. Nothing here is sampled.

**Scope, stated first.** What follows proves a statement about a FIXED multiplier
pair and shift: exactly one return value exists per (multiplier pair, shift) in the
order-r subgroup. It is NOT a production mixing bound, NOT a statement about the
full walk on 2^251 classes, and NOT a minimality theorem for step design. The only
minimality-adjacent statement made anywhere in this package is empirical: among the
step rules actually tested on the shipped fixture (plain multiply, inversion twist,
column permutation, row partition, squaring), squaring is the one that admits no
inverse-pair family; untested rules are not covered.

## 1. The identity

On the order-r subgroup of GT, sigma^6 = inversion: p^6 = -1 mod r, so the sixth
Frobenius power sends y to y^(p^6) = y^(-1).

Measured on the shipped fixture: `sigma6_is_inversion = True`,
`inverse_pair_columns = 384 / 384` (in EVERY row t and every k,
`TBL[t][k+6] * TBL[t][k] = 1`), and `orbit_columns = 768 / 768` (every column is
`sigma^(-k)(TBL[t][0])`).

## 2. Why the inverse-pair family cannot be partitioned away

Class-level merging needs `successor class = f(class)`. For sigma-partners
`y1 = sigma^j y2`, that forces `m(kappa - j) = sigma^j m(kappa)` for every j,
i.e. the multiplier map must be `m(k) = sigma^(-k)(A_t)` for a single per-row
`A_t`. Any such table therefore has

    m(k+6) = sigma^(-k-6)(A_t) = sigma^(-k)(sigma^(-6) A_t) = m(k)^(-1).

The inverse pair is an identity of the sigma-orbit structure, not a table defect.
A per-row column permutation that separated antipodal pairs would be non-affine,
would break `m(kappa - j) = sigma^j m(kappa)`, and would destroy class
determinism (measured cost in v22: retention 0.808333, success 0.3523, sub-floor).
Row partitioning cannot help either: both offending multipliers are columns of the
SAME row, by definition of the family. The family is representative-independent
(`kappa_delta_partner_invariant`), with rate
`P(t' = t) * P(kappa' - kappa = 6) = 1/(12 R) = 2^-17.585` per multiply step at
R = 2^14; a trail of the shipped cap 2^87 accumulates about 2^69 expected
absorptions, so the guard-free history-free walk cannot reach even one
distinguished point at production width. v23's trap row (theta^2/2N = 2^-92.3)
prices the ordinary rho coincidence, not this family, and is off by about 2^75.

## 3. The family-free step and its return equation

Make the step `y <- (y * m(kappa))^2`. It is still equivariant:

    (sigma^j y * sigma^j m)^2 = sigma^j((y m)^2),

so the successor class is a function of the class alone -- no prev_row, no ring,
no history -- and a class repeat IS a state repeat, so the cycle abort is sound
(this is exactly what v24's `step_v24` broke and its verdict names).

With `sigma^j(y) = y^(p^j)` the two-step return is

    y^(4 - p^j) = (m1^4 m2^2)^(-1),        not   y^3 = sigma^j / (m1^4 m2^2).

The exponent `e_j = (4 - p^j) mod r` is invertible on the subgroup order r
(`gcd(e_j, r) = 1` for j = 0..3), so the return has EXACTLY ONE solution per
(fixed multiplier pair, shift). Dynamically on the shipped object: 4 j-arms x
10,000 samples = 40,000 exact solves, 0 misses, 0 subgroup failures; the j=0 arm
reproduces the old cubing exactly (regression kept, 10,000/10,000) and the old
cubing returns 0 of 30,000 for j != 0. The v26 cubing form quoted above is
WITHDRAWN for j != 0, and with it the sentence "no structural family exists" --
what is established is uniqueness for a fixed pair and shift, not a mixing bound.

Price, measured on executed counts rather than estimated: the multiply step gains
exactly one Granger-Scott cyclotomic squaring,
`cyclo_sqr = 24*FP_MUL + 122*FP_ADD = 11632` instructions, so `mul_step` goes
21889 -> 33521 and `c_stress` 0.964032 -> 1.346817, i.e. **+0.4824 bits** (v25's
estimate was +0.473; the measured delta is larger because the squaring fraction is
priced at the executed pooled 0.3473967898 = 14566/41929 rather than the declared
3/8, both reported in the ledger). Claim: **126.2672** = 78917/625, which is 0.4472
bits BEHIND the live board minimum 125.82. That rank loss is disclosed, not
hidden: the squared step is a validity play, not a rank play.

## 4. The column scan has teeth (planted control)

A scan that always reports "clean" is worthless. `theorem_scan(plant=(7,3))`
mutates ONE table entry, `TBL[7][3]`, and re-runs the same column scan. Result:
`clean = False`, `inverse_pair_columns 383/384` with failure `[7, 3, 9]` and
`orbit_columns 767/768` with failure `[7, 3]` -- the mutated column and its inverse
partner -- while the other 383 pairs stay clean. The clean scan on the unmutated
table reports 384/384 and 768/768. Both rows are in
`evidence/v27-measure-theorem.json` (`planted_flagged: true`).

## 5. Cycle census: exact map vs class-level map

`full_map_census` builds the step map on ALL r-1 = 9420 elements of the fixture and
returns the full cycle-length histogram. Measured (`evidence/v27-measure-cycle.json`):

* v27 squared step, EXACT-element map: 4 cycles, lengths {216: 2, 24: 1, 60: 1};
  exact 2-cycles = 0.
* v13 rule (control), exact map: 21 cycles, lengths {16: 12, 312: 1, 63: 4, 9: 4};
  exact 2-cycles = 0.

Three honest consequences, none argued away:

1. The exact-element map has no 2-cycle, but the CLASS-level (collapsed) map still
   reports ONE 2-cycle, disclosed in `evidence/v27-audit.json`. A class 2-cycle
   lifts to an exact cycle of period `2*ord(sigma^J)`, so the absence of exact
   two-cycles is NOT offered as absence of class-level absorption. The class-level
   abandon rate is what the trap row prices; it is not claimed zero.
2. v26 collapsed 9420 exact states onto 785 class keys by `dlog[class_key]`, a
   factor 12.00; v27 keys by the full twelve-value signature, which is why the
   cycle lengths above differ from the v26 record.
3. The v13 control did not reproduce the family as an exact 2-cycle at this
   fixture (0 measured, 12.27 predicted by 1/(12R) at R = 64): the exact-return
   condition also requires `t(y') = t(y)`, which the routing decorrelates at this
   table size. The 2-cycle census is therefore the wrong instrument for this
   family and is reported as a FAILED control, not used as load-bearing evidence.

The load-bearing evidence for family absence is the algebraic uniqueness of
section 3 (`y^(4-p^j)` has one solution per fixed pair and shift) plus the 384/384
identity and its planted mutation in section 4. The v23 family itself was measured
over long trails (rate 2^-17.585 per multiply step), not as a 2-cycle, and that
measurement is carried in the v25 record.

## 6. Consequence

At R = 2^14 the class-merge rho family on this track is either unsound
(history-free), rejected (prev_row history), or 0.4472 bits over the 125.82 board
minimum (squared step). The walk family, not the packaging, is the blocker. What
the N=785 fixture cannot do is measure production retention: it is a toy at
`S/theta <= 16` against production's `2^45.8`, and no production rate is claimed
from it anywhere in this package.
