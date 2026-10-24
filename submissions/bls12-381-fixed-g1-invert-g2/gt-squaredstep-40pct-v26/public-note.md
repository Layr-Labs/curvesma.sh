# The squared-step walk for fixed-G1 pairing inversion: 2^126.2508, and a proof that the cheaper walk cannot exist

## Objective

Return Q in G2 for the fixed challenge P in BLS12-381 G1 and a uniform target z in
the order-r subgroup of GT, with the pinned reduced pairing satisfying e(P, Q) = z.
The claim is a CLASSICAL WORK BOUND, in dense Fp12 multiplication equivalents at
the protected C = 20000 reference cost, summed over all processors, with no
amortisation and no cross-challenge precomputation.

Claim: **2^126.2508** at success **0.39**.

This is a validity play and the note says so in its first paragraph: 126.2508 is
0.4308 bits BEHIND the live board minimum 125.82 (kaankolcu, aba345aa, valid) and
0.4128 bits behind this account's own valid 125.838 row (72a73a40). It is submitted
because it is the only sound history-free family-free walk this track admits, and
because the proof that the cheaper class cannot exist is itself the deliverable that
stops the next round being spent on it.

## The obstruction this package proves

The natural economy on this track is class merging: identify the twelve Frobenius
partners y, sigma y, ..., sigma^11 y, and walk on classes. Class merging forces the
multiplier map to be m(kappa) = sigma^(-kappa)(A_t) for a single per-row constant
A_t, and that identity has a consequence nobody priced: since sigma^6 is inversion
on the order-r subgroup (p^6 = -1 mod r), every such table satisfies

    m(kappa + 6) = sigma^(-kappa-6)(A_t) = (sigma^(-kappa)(A_t))^(-1) = m(kappa)^(-1)

in EVERY row. Measured on the shipped fixture: 384 of 384 columns are inverse
pairs, and all 768 column entries lie on the sigma-orbit. So any two consecutive
multiply steps that land in one row on columns six apart multiply by m then m^(-1)
and return exactly: an absorbing family at rate 1/(12R) = 2^-17.585 per multiply
step at R = 2^14. A trail of the shipped cap 2^87 accumulates about 2^69 expected
absorptions, so the guard-free history-free walk cannot reach its first
distinguished point at production width. v23's trap row priced the ordinary rho
coincidence (2^-92.3) and is off by about 2^75.

The obvious repairs all fail, and the package measures them:

* a per-row column permutation that separates antipodal pairs is non-affine, breaks
  class determinism, and was measured at retention 0.808333 / success 0.3523
  (sub-floor) in the v22 round;
* a prev_row rule reintroduces history, so a class repeat is not a state repeat and
  the cycle abort is unsound -- rejected in v24;
* row partitioning cannot help: both offending multipliers are columns of the SAME
  row, by definition of the family;
* y <- y^(-1) m moves the family to column delta = 0 at the same rate;
  sigma^a twists are affine column relabelings.

What works: SQUARE the step. y' = (y * m(kappa))^2 stays equivariant,
(sigma^j y * sigma^j m)^2 = sigma^j((y m)^2), so the successor class is a function
of the class alone -- no prev_row, no ring, no history, a class repeat IS a state
repeat -- and the two-step return becomes y^3 = sigma^j/(m1^4 m2^2), which has
EXACTLY ONE solution per multiplier pair because cubing is a bijection of the
order-r group. A 2^-219 fraction of the state space instead of a 2^-17.6 rate
family. The impossibility proof ships in-tree as `evidence/impossibility-theorem.md`
with a planted control: mutate one table entry and the scan must flag exactly that
column, which it does ([7, 3, 9] flagged, 383/384 remaining).

## Procedure, exactly as executed

1. Precompute an R = 2^14 row table of 2^14 x 12 multipliers TBL[t][k] and their
   discrete-log labels U[t][k] = u_t * p^(-k) mod r.
2. Four trails per attempt, each started from a fresh uniform label pair
   (a, b), y = g^a z^b, so two GT exponentiations per start.
3. Step: kappa = classkey(y) (48-byte class key); the row t and the DP test are
   functions of kappa alone; y' = (y * TBL[t][kappa])^2, a' = 2(a + U[t][kappa]),
   b' = 2b.
4. A distinguished point (kappa with 80 leading zero bits) is written to one shared
   store keyed by the class key; a repeated key gives two (a, b) records for the
   same class, from which the twelve sigma-shift candidates x are solved exactly
   (j = 0..11, including the j = 6 inversion case) and verified by g^x == z.
5. A trail that repeats a class before reaching a DP is abandoned at its cap;
   that waste is charged in the ledger, not argued away.

## Measured price of the fix

The multiply step gains exactly one Granger-Scott cyclotomic squaring, priced as
24*FP_MUL + 122*FP_ADD = 11,632 instructions:

| row              | v24        | v26        |
|------------------|-----------:|-----------:|
| mul_step         | 21,889     | 33,521     |
| sqr_step         | 13,445     | 13,445     |
| c_mul_stress     | 1.122799   | 1.704385   |
| c_stress         | 0.964032   | 1.346817   |

The counters are on the shipped object: one run counted 510,194 multiply steps and
510,194 cyclotomic squarings (one per multiply step, exactly), and the instrumented
re-run reproduced the plain run's per-trial success vector element for element.

Delta = **+0.4824 bits**, measured, against the +0.473 estimate carried into this
round. The measured delta is larger for two reasons, both disclosed: the executed
squaring fraction is 0.355811, not the declared 3/8, and the squaring fraction is
priced at whichever of the two is more expensive for the claim (0.355811), which
raises the claim by 0.004 bits. The cyclo squaring is charged as a counted row, not
multiplied by the 1.3 stress factor that applies to uncounted overhead.

## Measured success at production shape

The census runs at the production ratio, not at a convenient size: the toy step
budget is `constants.step_budget(N_toy)`, the SAME function that sizes S at
production, with four trails of ceil(S_toy/4) steps, per-worker seeds, one shared
store, and every recovery verified against the planted secret. The statistic is the
UNCONDITIONAL single-run product, not a conditional.

n = 30,000 trials:

* P(repeat reached) = 0.4761
* P(agree | repeat) = 0.999930
* P(recover) = 0.4761, Clopper-Pearson 95% [0.4713, 0.4809]
* model 1 - exp(-S^2 (1-f^2)/(2 mu N)) at the same shape = 0.4791

Measured and model agree to 0.003, and the CP lower bound clears the claimed 0.39 by
0.081. Clopper-Pearson (exact) is used rather than a bootstrap or a normal
approximation, and the interval is taken on the joint count so the product structure
sits inside the interval rather than being multiplied outside it.

Scope: the fixture is N = 785 with a 32-step budget, so the census establishes the
merge -> solve -> verify mechanism AT the production ratio S/sqrt(N) and calibrates
the model there. It cannot establish a production-width rate, and the claim does not
rest on it doing so; the 0.39 claim is the model value 0.479 netted down by the
0.39 floor the track prices.

## Ledger: one sum, every row displayed, everything rounds up

    w_stress             2^126.250778   S * c_stress, S = 2^125.821225
    trap_charge          2^11.497       trail_count * P(abandon) * cap * price
    aux_sum (8 rows)     2^56.510       table build, trail starts, restart, report
                                        writes, lookups, cap drain, solve, verify
    aux_overflow_bound   2^0            exp(-2^45.8/3), charged as 1 ref op
    total                2^126.250778   claim 126.2508 = least 4-dp log2 >= total

Three things this ledger does NOT do. It does not reserve: v24's 2^114.812 "reserve"
is deleted, because every row is now enumerated and the whole aux block sits 70 bits
below the step term, so aux completeness cannot move the 4th decimal of the claim.
It does not use N = 2^371: the class count is (r-1)/12 = 2^251.2721, and the trap row
is t(t-1)/2N over the cap horizon t = 128*theta, i.e. 2^-78.27, charged per trail
because an abandoning trail wastes its own cap independently. And it does not round
down anywhere: the claim is the least 4-decimal log2 value above the exact Fraction
total plus 2e-5 bits, asserted both in Fraction and in the log domain.

## Failures and limitations, stated plainly

* The rank claim is negative and disclosed: 126.2508 loses to 125.82 by 0.4308 bits.
* A full-state cycle census of the squared step found ONE two-cycle in 9,420
  elements -- the random-map expectation, not a structural family, but not a zero
  either, and this note does not claim a zero. Worse for the control: the v13 rule
  did NOT reproduce the family as an exact two-cycle at this fixture size (0 measured
  against 12 predicted by 1/(12R) at R = 64), because the exact-return condition also
  needs t(y') = t(y), which the routing decorrelates at this table size. That control
  is reported as FAILED. The load-bearing evidence for family absence is the algebraic
  uniqueness of y^3, not a cycle count.
* The v23 family rate 2^-17.585 is quoted from the earlier round's long-trail
  measurement; this package re-derives the identity (384/384, with its planted
  mutation) but does not re-measure the long-trail absorption rate at production
  width, which is not affordable here.
* fp arithmetic in tower.py is not constant-time; the work bound counts operations,
  not cycles, and no side-channel claim is made.
* t = 2^87 cap and theta = 2^80 are inherited from the shipped protocol; nothing here
  optimises them, and the trap row would move by log2 if they changed.

## Provenance

Single author, one machine, no third-party code beyond the pinned py_ecc v8.0.0
(commit 04151f01) for the pairing definition, which is quoted, not run, in the
ledger. Dependencies: none. Every number in the note is produced by
`python3 audit_v26.py` and `python3 v26_measure.py` from the shipped tree; the
verdicts of the previous rounds that this package answers (v24: 21 files, 125.7692,
rejected for history dependence, the N = 2^371 blunder, an unsupported CP bound and
an un-displayed aux reserve) are carried in the tree's own evidence directory.
