# Impossibility of the table-geometry route, and why the squared step is the minimal sound fix

Shipped as evidence with the v26 package. Every number below is produced by
`v26_measure.py theorem` and `v26_measure.py cycle` on the shipped fixture
(`toy_census.Toy(bits=5)`: r = 9421, 64 rows x 12 columns = 384 inverse-pair
columns). Nothing here is sampled.

## 1. The identity

On the order-r subgroup of GT, sigma^6 = inversion: p^6 = -1 mod r, so the sixth
Frobenius power sends y to y^(p^6) = y^(-1).

Measured on the shipped fixture: `sigma6_is_inversion = True`,
`inverse_pair_columns = 384 / 384` (in EVERY row t and every k,
`TBL[t][k+6] * TBL[t][k] = 1`), and `orbit_columns = 768 / 768` (every column is
`sigma^(-k)(TBL[t][0])`).

## 2. The impossibility

Class-level merging needs `successor class = f(class)`. For sigma-partners
`y1 = sigma^j y2`, that forces `m(kappa - j) = sigma^j m(kappa)` for every j,
i.e. the multiplier map must be `m(k) = sigma^(-k)(A_t)` for a single per-row
`A_t`. Any such table therefore has

    m(k+6) = sigma^(-k-6)(A_t) = sigma^(-k)(sigma^(-6) A_t) = m(k)^(-1).

The inverse pair is not a table defect that can be partitioned away: it is an
identity of the sigma-orbit structure. A per-row column permutation pi_t that
separated antipodal pairs would be non-affine, would break
`m(kappa - j) = sigma^j m(kappa)`, and would destroy class determinism (measured
cost in v22: retention 0.808333, success 0.3523, sub-floor). Row partitioning
cannot help either: both offending multipliers are columns of the SAME row, by
definition of the family.

The family is representative-independent (`kappa_delta_partner_invariant`), so it
is a property of the class map with rate

    P(t' = t) * P(kappa' - kappa = 6) = 1/(12 R) = 2^-17.585 per multiply step at R = 2^14.

A trail of the shipped cap 2^87 accumulates about 2^69 expected absorptions, so
the guard-free history-free walk cannot reach even one distinguished point at
production width. v23's trap row (theta^2/2N = 2^-92.3) prices the ordinary rho
coincidence, not this family, and is off by about 2^75.

## 3. The only family-free class-merge step, and its price

Make the step `y <- (y * m(kappa))^2`. It is still equivariant:

    (sigma^j y * sigma^j m)^2 = sigma^j((y m)^2),

so the successor class is a function of the class alone -- no prev_row, no ring,
no history -- and a class repeat IS a state repeat, so the cycle abort is sound
(this is exactly what v24's `step_v24` broke and its verdict names).

The two-step return equation becomes `y^3 = sigma^j / (m1^4 m2^2)`, which has
EXACTLY ONE solution per multiplier pair (cubing is a bijection of the order-r
group), a 2^-219 fraction of the space: no structural family exists.

Price, measured on executed counts rather than estimated: the multiply step gains
exactly one Granger-Scott cyclotomic squaring, `cyclo_sqr = 24*FP_MUL + 122*FP_ADD
= 11632` instructions, so `mul_step` goes 21889 -> 33521 and
`c_stress` 0.964032 -> 1.346817, i.e. **+0.4824 bits** (v25's estimate was +0.473;
the measured delta is +0.4824 because the squaring fraction is priced at the
executed 0.353722 rather than the declared 3/8, both reported in the ledger).
Claim: **126.2508**, which is 0.4308 bits BEHIND the live board minimum 125.82.
That rank loss is disclosed, not hidden: the squared step is a validity play, not
a rank play.

Cheaper equivariant twists do not work: `y <- y^(-1) m` just moves the family to
column Delta = 0 at the same 2^-17.6 rate; `sigma^a` twists are affine column
relabelings and change nothing.

## 4. The theorem has teeth (planted control)

A scan that always reports "clean" is worthless. `theorem_scan(plant=(7,3))`
mutates ONE table entry, `TBL[7][3]`, and re-runs the same column scan. Result:
the scan reports `clean = False` and flags exactly column pair `[7, 3, 9]` -- the
mutated column and its inverse partner -- while the other 383 pairs stay clean.
The clean scan on the unmutated table reports 384/384 and 768/768. Both rows are
in `evidence/v26-measure-theorem.json`.

## 5. Cycle census: what it does and does not show (reported as measured)

`full_map_census` builds the step map on ALL r-1 = 9420 elements of the fixture and
returns the full cycle-length histogram. Measured:

* v26 squared step: cycles 3, lengths {2: 1, 5: 1, 36: 1}; 2-cycles = 1.
* v13 rule (control): cycles 4, lengths {3: 1, 16: 1, 21: 1, 26: 1}; 2-cycles = 0.
* Random functional graph on 9420 nodes: expected 2-cycles = 1/2, expected cycle
  nodes ~ sqrt(pi n / 2) = 122. Both maps sit at the random expectation.

Two honest consequences, and neither is argued away:

1. The v26 two-cycle count is ONE, not zero. It equals the random-map expectation
   and is 2^-13.2 of the space; it is NOT evidence of a structural family, but it
   is also not a zero, and this document does not claim a zero.
2. The v13 CONTROL DID NOT REPRODUCE THE FAMILY as an exact 2-cycle at this
   fixture (0 measured, 12 predicted by 1/(12R) at R = 64). The exact-return
   condition also requires t(y') = t(y), which the routing decorrelates at this
   table size, so the 2-cycle census is the wrong instrument for this family and
   it is reported as a FAILED control rather than used as load-bearing evidence.

The load-bearing evidence for family absence is therefore not a cycle count. It is
the algebraic uniqueness in section 3 (`y^3` has one solution per multiplier
pair, so no multiplier pair admits a return for more than a 2^-219 fraction of
states) plus the 384/384 identity and its planted mutation in section 4. The
v23 family itself was measured over long trails (rate 2^-17.585 per multiply
step), not as a 2-cycle, and that measurement is carried in the v25 record.

## 6. Consequence

At R = 2^14 the class-merge rho family on this track is either unsound
(history-free), rejected (prev_row history), or 0.48 bits over the 125.82 board
minimum (squared step). The walk family, not the packaging, is the blocker. What
v26 ships is the sound one.
