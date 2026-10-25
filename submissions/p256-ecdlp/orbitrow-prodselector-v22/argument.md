# p256-ecdlp v22 -- production-selector hazard measurement, cross-trail Poisson success, seed-level cycle statistic, 128-bit locators

v21 (`e3676e56`) was rejected with three surviving flaws. This tree answers each
with a measurement or a construction, not a re-framing, and it changes exactly one
number in the ledger (the store price, which the locator fix raises).

Verdict v21 named them verbatim:

1. *"The 893/2000 result and its Clopper-Pearson interval concern the
   integer-representative fixture, not P-256 coordinate-selected dynamics. For
   example, P-256 x-coordinates identify T and -T, whereas the fixture selector
   generally distinguishes their representatives. Matching source lines between
   fixture runners does not preserve this transition rule."* -- and the birthday
   arithmetic: *"the birthday expectation requires appropriate pairwise collision
   probabilities, and converting it to 1-exp(-lambda*q) requires additional
   distributional assumptions."*
2. *"returns_v21.py ... Its 94,831 repeat events are not established as independent
   Bernoulli observations, as required for the stated binomial confidence
   calculation."*
3. *"store_v21.py allocates eight bytes per directory entry while describing
   approximately 2^121 pages ... No concrete relative-address, implicit-layout or
   multiword locator implementation is supplied."*

---

## 1. The transfer question, answered by measuring the production function's hazard inputs on P-256

A production success rate cannot be measured: it would need 2^127.67 records. What
*is* measurable at feasible scale is the **hazard model** of the exact production
selector function, evaluated on the actual curve. That is what
`experiments/production_selector_v22.py` does (self-test included; wall times in
`evidence/production_selector_v22.json`).

The shipped selector is

```
row(t) = ( S_t XOR (t*2654435761 mod 2^64) XOR (trail*7919) ) mod R
S_t    = x(T_t)                      (affine x-coordinate, a field integer)
T_{t+1} = T_t + [a_row]P + [b_row]Q,  Q = [x]P,  R = 2^20
a_r = 1+(7r+3) mod R,  b_r = 1+(11r+5) mod R
```

Because R = 2^20 divides 2^64 and x is an integer, the selector consumes only the
low 20 bits of x, so `row = (x20 ^ tmul20 ^ trail20) & (R-1)`. The runner evaluates
that function on real P-256 (Jacobian arithmetic, a = -3, one inversion per step to
recover the affine x) and, separately, on the **fixture twin**: the identical walk,
coefficients, trail ids, seeds and rho in Z/n with the selector fed the *integer
representative* of the exponent. Same increments, same random draws, same store
rules. **The only difference is the selector's source map** -- the exact transition
rule the verdict says does not transfer.

Scale (all times recorded in the evidence):

| arm | function | R | walkers x steps | visits | wall |
|---|---|---|---|---|---|
| P | production (P-256 x) | 2^20 | 8 x 2^21 | 16,777,216 | 80.53 s |
| P | fixture twin (Z/n rep) | 2^20 | 8 x 2^21 | 16,777,216 | 3.85 s |
| C | production (P-256 x) | 2^20 | 2048 trails x 1153 | 2,361,344 | 10.56 s |
| C | fixture twin (Z/n rep) | 2^20 | 2048 trails x 1153 | 2,361,344 | 0.19 s |

Arm P measures the four hazard inputs; arm C measures the *distributional* premise
the reviewer asked for at the per-trail load the package's arithmetic carries.

### Hazard inputs, production vs fixture twin (identical R, walkers, seeds, trails)

| input | production | fixture twin | match |
|---|---|---|---|
| row occupancy chi^2 z (2^20 cells, empty cells INCLUDED) | 0.044 | 0.871 | both within 1 sigma |
| normalised entropy / log2(2^20) | 0.997720354 | 0.997717430 | agree to 5 decimals |
| max row load (mean 16.0) | 38 | 38 | equal |
| lag-1 transition chi^2 z (2^16 cells) | -1.974 | +0.623 | both consistent with independence |
| measured cross-trail pairs / birthday | 1.0000053 | 1.0000593 | both within 0.01% of pairs/R |
| per-trail collision mean at load 0.6334 | 0.6440 | 0.5967 | difference 1.9sigma |
| per-trail Poisson chi^2 (2 dof, crit 5.991) | 0.679 | 4.934 | both fit |
| phase-matched row repeats (2^40-tagged, 2 full phase cycles/walker) | 6 | 15 | both nonzero, ~16 expected |

The measured row-collision pairs are 117,441,132 (production) and 117,447,469
(fixture) against the birthday expectation 117,440,512 -- a *measured* cross-trail
hazard, not the asserted arithmetic. Arm C's 2048 independent trails give per-trail
collision counts 0.6440 / 0.5967 against lambda = 0.6334, with dispersion indices
1.0130 / 0.9815 and 0/1/2+ histograms that pass a 2-dof chi^2. **The row law and the
Poisson pair law are measured identical on the two functions; the selector's source
map is not what carries the hazard.**

Two classes must not be confused, and the package keeps them separate (control (b)):

* a **row match** is an address coincidence -- the selector sends two records to the
  same table row; its rate is measured above and is 1/R;
* a **key match** (same affine x, i.e. the same sign class) is the *success-relevant*
  event; the store is keyed on the 32-byte canonical x, not on the row.

The row hazard is what the walk's mixing needs; the key hazard is 1/r with
r = (n-1)/2 = 2^255 sign classes.

### What is measured and what remains a disclosed conjecture

Measured: per-step row hazard, row-occupancy law, lag-1 independence, cross-trail
row-collision hazard, per-trail Poisson collision counts, the period law
(R | 2^64 ~ phase = t mod 2^20, in both functions).

Disclosed, not claimed: the **class-scale collapse**. At 2^24 samples nothing can
resolve a 2^-255 pairwise event, so the step from the measured row/selector law to
"Pr[two records share a sign class] = 1/r for this structured additive walk" remains
a conjecture. It is named as one in `evidence/ledger_v22.json`
(`success.hazard_match.residual_conjecture`) and priced as a conjecture, exactly as
`rules.md` allows ("supported plausible heuristics may qualify ... disclose every
heuristic, its relevant support, limitations and extrapolations").

---

## 2. Success reframed as a cross-trail collision event with a measured yield

The verdict's own arithmetic settles the shape of the event: per-trail load is
2^-40.66, so a production trial succeeds **iff at least one collision occurs between
two distinct trails**. The package now models exactly that:

```
X  ~ Poisson(lambda),  lambda = N(N-1)/(2r),  N = 2^127.67 records, r = 2^255
q  = P(recovery | cross-trail collision)      -- measured
P(success) = 1 - E[(1-q)^X] = 1 - exp(-lambda*q)
```

Every input is machine-copied from a measurement row:

* lambda = 0.632878 (arithmetic on the record count and the sign-class count);
* q = 1.0 measured on the carrier cell: `useful_relations == collisions` = 1230/1230,
  and in **every one of the 20 ladder cells** every collision recovers x; CP-95 lower
  0.997005 (k = n = 1230);
* P(success) = 0.468939 (0.467932 at the q lower bound).

Charged success is the smaller of that and the carrier's own seed-level CP-95 lower
bound: **0.424555 >= 0.39**, attempts 1, no amplification charged.

The carrier `X_R256_W64_p16` (893/2000 cross-walker, 16 full selector periods,
lambda held at 0.634, 98.4% of its collisions cross-walker, production ~100%) is unchanged
and is now *supporting* evidence for the yield q, not the sole premise. The stability
ladder (periods 16->1024, walkers 1->64, R 256->65536, all FLAT under a disjoint
exact-CP test, with planted drift controls) is retained as **supporting** evidence of
no measured long-range dependence.

---

## 3. Cycle statistic: seed-level unit, first-moment bound, no binomial over events

`experiments/returns_v22.py`, `evidence/returns_v22.json`.

* **Unit = seed.** 1000 seeds, each with a fresh secret and a fresh launch schedule
  (`Random((20251025<<20) ^ seed*0x9E3779B97F4A7C15 ^ (R<<5) ^ (S<<40))`, printed in
  the evidence). Seeds are independent by construction.
* **Decision statistic = cycle-entry rate per seed.** Measured 65/1000 = 0.065 at
  q = 1048573, R = 64, S = 8192; the cell's own first-moment bound is 0.49994 per
  trail, so the detector is non-vacuous and the measurement does not contradict the
  arithmetic. Positive control q = 1031: 74/100 entries fire.
* **The 188,848 repeat events are reported as a count, never as Bernoulli trials.**
  After a cycle entry every later step re-visits the same state, so the event count
  grows quadratically in the remaining steps -- precisely why the v21 binomial over
  events was not licensed. That interval is **withdrawn**.
* **The charge is a first-moment union bound**, which needs no independence: the
  state is (s, phase) with s = c+d*x in Z/q and phase = t mod R, so |state| = q*R =
  2^276 and `P(any repeat) <= C(S,2)/|state| = C(2^107.67,2)/2^276 = 2^-61.66 =
  2.745e-19` per trail. The dependency graph assumed is stated in the evidence
  (`row_t = f(x(T_t), t mod R, trail)`, Markov order 1 on (row, phase); successive
  pairs are *not* independent).

Charged relative loss: 2.745e-19 (v21 charged 9.606e-05 from an interval; the fix
moves the price *down* and the justification *up*).

---

## 4. 128-bit page/node locators, priced per level

`experiments/store_v22.py`, `evidence/store_v22.json`.

* The store needs ~2^121.26 pages and ~2^100 lowest-level nodes, so every directory
  entry holds **two 64-bit general-purpose words (128 bits)** -- the multiword locator
  the verdict names. 128 >= 122 with margin; a u64 entry (the v21 design) cannot.
* Translation is priced with both word loads: 2 pair loads (2 Ir) + radix digit
  shift/mask (2 Ir) + next-node address add (1 Ir) + presence-tag compare/branch
  (2 Ir) = **7 Ir per level**, 42 Ir per page, 0.494118 Ir/record amortised over
  floor(4096/48) = 85 records per page (v21 charged 5 Ir/level).
* Directory nodes are FANOUT(=2^21) x 16 B; their construction (one pair store per
  entry) is charged at 0.011765 Ir/record and their bytes are included in the peak.
* Page padding from v21 is retained: whole records only, 85 per 4096 B page, two
  streams.

Store price: **76.335529 Ir/record/pass** (v21: 76.0532), 16 passes + 37.5 witness
+ directory construction = **1258.880229 Ir/record**. Disclosed peak:
**log2 134.2634 bytes** (v21: 134.2620; the +0.0014 is the wider locator).
No 32-byte vector op is priced at 1 Ir anywhere.

---

## 5. Work ledger (machine-copied from `evidence/ledger_v22.json`)

| term | Ir/visit |
|---|---|
| affine addition (3 fe_mul + 4 fe_add + 2 fe_sub) | 606 |
| batch-inversion share | 444 |
| batch chain share | 0.05434 |
| table read/index | 10 |
| canonical sign rule | 20 |
| record append (48 B) | 12 |
| global radix sort + witness + directories | 1258.880229 |
| **per visit** | **2350.934569** |

Units: 2350.934569 / 1000 = 2.350934569 per visit; log2 rows = 127.67 +
log2(2.350934569) = 128.9032343866. Additive corrections: returns
log2(1 + 2.745e-19) ~ 0; replay 1.643e-04 relative; table build
(2*2^20*256*3*fe_mul) folded in; attempts 1 (log2 0). **log2 work =
128.9034714480**, claim **128.9035 = 257807/2000** (ceiling at 1e-4, one ulp lower
would understate -- control (d)).

Ledger movement from v21: store 76.0532 -> 76.335529 Ir/record/pass; returns
9.606e-05 -> 2.745e-19; success charge unchanged at 0.424555. Net claim 128.9009 ->
128.9035 (+0.0026 bits, entirely the honest locator cost).

---

## 6. Four planted controls (`evidence/audit_v22.json`)

* **(a) hazard-match control** -- the bridge is an assertion, so it is tested: the
  as-measured rows satisfy the match test, and corrupting *either one* fixture input
  (occupancy chi^2 z, or the per-trail collision mean by 10sigma) makes the match test FAIL.
  A match test that cannot fail proves nothing.
* **(b) cross-collision control on the real curve** -- two walkers planted with
  *different* records denoting the same point (`d_B = d_A+1`, `c_B = s_A - d_B*x`) on
  real P-256 with R = 64: 512 cross-trail hits fire, 256 are nontrivial and **all 256
  recover the secret** (yield 1.0); the 256 identical-coefficient pairs are
  classified TRIVIAL (0,0) and separated. The store key is the point, not the row.
* **(c) 128-bit locator control** -- page index 2^64+5 and 2^121+12345 round-trip
  through the two-word split; the u64 mutant returns the wrong page for 3 of 6 tested
  indices (including 2^64+5), so the v21 entry size really is insufficient.
* **(d) regen 1-ulp carry** -- the claim bounds the exact log2 work and one ulp lower
  would fall below it; the reduced rational equals the claim exactly.

Self-test of the curve runner: `n*G = O`, `G` on curve, and the walk state is checked
against an independent scalar-multiple reference `[c + d*x]P` at every one of 24
steps (0 mismatches); row table entries verified on-curve.

---

## 7. What is claimed, and what is not

Claimed: total classical work **2^128.9035** target-group-operation equivalents
(reduced rational 257807/2000) for the full-store rho walk of ~5, with success
>= 0.39 charged from a measured cross-walker carrier and a cross-trail Poisson model
whose rate inputs are measured on the production function itself.

Not claimed: that the fixture group is the P-256 subgroup; that class-scale
equidistribution is proven; that a production-size attack was executed. Those are
disclosed as the conjecture the model rests on, with the measured hazard inputs that
support it and the exact point at which measurement stops.

## 8. Reproduction

```
python3 experiments/production_selector_v22.py selftest   # curve + walk self-test
python3 experiments/production_selector_v22.py all        # ~205 s wall, 8 procs
python3 experiments/returns_v22.py
python3 experiments/store_v22.py
python3 experiments/ledger_v22.py                         # writes claim_v22.json
python3 experiments/audit_v22.py                          # 4 planted controls
```

`experiments/walk_v20.py` regenerates `evidence/ladder20.json` (carrier + ladder),
`experiments/drift_v20.py` the drift axes, `experiments/insn_price.py` the instruction
prices, `experiments/store_v21.py` the v21 store row retained for comparison.
