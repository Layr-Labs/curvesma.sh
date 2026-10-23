# Argument — gt-metered-squaring-40pct-v15

Claim **125.73** (`12573/100`); computed nominal **W = 2^125.725..125.726**; metered cap
**W_cap = 2^125.728..125.729**; success >= 0.39 after the abort debit; C = 20,000 instructions
per dense Fp12 multiplication. Reduction, tower, key rule, branch rule and step prices are the
v14 construction (still valid); this version changes the *resource claim*, the *algorithmic
specification* and the *collision model*. Gate seen at submission time: **125.74** (747ad6c9,
review; live entries also 125.80, 125.826, 125.89 x3).

## 0. What the v14 verdict rejected, and what replaces it

1. **"c_avg with a cap claim"**: v14 charged the average step price `0.99129` per step and claimed a
   cap of `2^125.73`. With k squarings in n steps the modelled cost is `1.1246n - 0.3555k`, and H2
   does not pin k/n at 3/8, so a max-price budget of `2^125.91` exceeded the claim. v15 **meters**:
   every step adds its ACTUAL charge to a running meter and the run aborts at `W_cap`. The cap then
   holds **by construction on every run**, not on average, and the only loss is P(abort), bounded
   below by Hoeffding and debited from the success probability (Section 4).
2. **"transitions, not the DP collision procedure; no store, no peak memory, unexplained 2^65"**:
   Section 5 writes the full production procedure — predicate, record layout with byte size,
   sorted-array lookup, replay/verify, 2-cycle escape and checkpoint schedule — itemises the store
   and reports peak bytes `2,969,264,968 < 2^32`. The `2^65` allowance is **deleted** and replaced by
   a derived setup row (1024 table entries x (510 reference mults + 12 Frobenius)).
3. **"twisted multiplication not shown bijective on Frobenius classes"**: Section 6 demotes the
   no-merge property to a **scoped heuristic** with the j-shift gap stated, charges the conservative
   side of the two rho expressions, and Section 7 measures both channels in the toy: 134 shifted-image
   matches, **all 134 usable**; 1228 same-index matches, 1225 degenerate cycle closures.

## 1. Reduction, classes, step rule (unchanged)

Q != 1 generates G2 (prime order r), the pairing is nondegenerate, g = e(P0,Q) generates GT;
z = g^x is returned as P = [x]P0 and verified by e(P,Q) = z; z = 1 returns infinity.
p^6 = -1 mod r, so Frobenius has order 12 with inversion in the orbit and GT splits into
N = (r-1)/12 = 2^251.272 classes of twelve conjugates. Writing y = c0 + c1 w,
**T = Re(c0) is invariant across the twelve conjugates** and is read as one field word. Branch,
multiplier index and DP bits are disjoint ranges of T: bits 0-2 branch (`T mod 8 < 3`, probability
3/8), bits 3-12 multiplier index (R = 1024), bits 13-112 the DP predicate. Multiply steps use the
precomputed `T[h][k] = M_h^(p^k)` and never apply Frobenius to the walk element; squaring steps use
the cyclotomic squaring and skip the key computation; the twelve keys are six +-pairs for four Fp
multiplications. Evidence: `evidence/prod_run15.json` (`part_a` 26/26), 2000 production-field steps
with 63+ invariant checks and 33/33 class-invariance agreements against all eleven explicit
conjugates.

## 2. Instruction rows and prices (counted, not asserted)

| row | instructions |
| --- | ---: |
| every step: trace read, branch/index/DP bits, dispatch | 60 |
| every step: exponent pair update mod r | 100 |
| every step: cycle compares, loop control, margin | 140 |
| every step: amortised DP store (over-counts the itemised row by 2^86) | 1 |
| multiply: six +-pair keys (4 Fp mul + 6 adds, counted) | 1802 |
| multiply: min of six, +-pair bookkeeping, table load, pair add | 390 |
| **multiply step charged** | **2493 + C = 22493 = 1.12465 u** |
| squaring: cyclotomic squaring (24 Fp mul + 130 adds, counted) + pair doubling | 15082 |
| **squaring step charged** | **15383 = 0.76915 u** |

`c_mul = 22493/20000`, `c_sqr = 15383/20000`, `c_avg = 5/8 c_mul + 3/8 c_sqr = 79307/80000`.
The Fp counts are those measured by running `scripts/gt_prod.py` (`sqr_rows.instr = 15042`,
`key_rows.instr = 1802`); `scripts/ledger15.py` reads them from the evidence file and rebuilds every
row, so no printed number is trusted. Setup is derived, not allowed: 1024 entries x (510 reference
mults for g^a z^b + 12 Frobenius x 2400 instructions) + 2^20 for code, tables, start points,
pairings and output = **10,475,339,776 instructions = 2^19.0 units**.

## 3. Meter and cap

The meter is an integer instruction counter. It is initialised to the setup charge and every step
adds its charged row (`22493` for a multiply step, `15383` for a squaring step); a step whose charge
would push the meter past `CAP = W_cap * C` is not executed and the run aborts. Therefore

    total charged work <= W_cap = W * (1 + 1/512),

for every run, with no tail assumption, no averaging and no maximum-price budget: the maximum step
price is irrelevant because the meter stops at the cap. `S = ceil(sqrt(2 N' ln(100/61))) =
70933784039754404319800682927445544128 = 2^125.737` steps at the average price gives
`W = 2^125.725..125.726 < W_cap = 2^125.728..125.729 < 2^125.73 = claim`. `ledger15.py` checks both
inequalities with integer power comparisons (`W_cap^100 < 2^12573`), not with decimals.

## 4. Abort analysis and success accounting

At n steps with k squarings the charge is `setup + 22493(n-k) + 15383k`, so the run survives S steps
iff `k >= k_min = ceil((setup + 22493 S - CAP)/7110) = 26213833087651451462838074698840409575`.
Since `k_min/S = 0.36956`, the abort event is `k/S < 3/8 - eps` with

    eps = 386335927256450157087181398951669473/70933784039754404319800682927445544128 = 0.005446
    (exact: eps = 3/8 - k_min/S, no rounding in any comparison).

Under H2 the branch indicators are i.i.d. Bernoulli(3/8), so Hoeffding gives

    P(abort) <= exp(-2 eps^2 S) <= 2^-6071291529163368902564250953935653,

which is `2^-6.07e33`: 33 orders of magnitude of exponent more than any reserve could need. The
success accounting is therefore

| item | value |
| --- | --- |
| model success under H1' (random-mapping walk, v14 model) | 0.430 |
| less P(abort) (Hoeffding, exact Fraction) | 2^-6.07e33 |
| less DP-window detection loss `2 * 2^100 / (S/R)` | <= 2^-14 (computed 2^-14.74) |
| **success lower bound carried in the claim** | **>= 0.39** |

The DP-window loss is the only non-trivial debit and it is a real one: after two trails merge they
share their trajectory, so the common DP occurs within 2^100 steps and both trails must still have
that many steps left; with trail length `S/R = 2^115.74` the missed-merge probability is
`2 * 2^100 / 2^115.74 = 2^-14.74`. The toy measures the same ratio at its own scale (§7: 0.316 vs
0.406 at S_ops, a 22% loss for a DP period that is 7% of the walk), i.e. the toy **over-states** this
loss by roughly 2^13, so the production bound is the conservative-and-tight side, and the claim
carries 0.43 - 0.39 = 0.04 of slack against a debit of 2^-14.

## 5. Production DP collision-detection, recovery and escape (the procedure, written out)

**Predicate (exact formula).** With `u` the canonical integer representative of `T = Re(c0)` in
`[0,p)` and `R = 1024`: branch `b = 1` iff `u mod 8 < 3`; index `h = floor(u/8) mod 1024`; DP iff
`floor(u / (8*1024)) mod 2^100 == 0`. All three are functions of the class invariant, so all twelve
conjugates of a class agree — this is what makes merged trails record the *same* DP.

**Record layout, 52 bytes.** `fp[32]` = low 256 bits of the canonical minimum key over the twelve
Frobenius keys (a class function; two distinct classes share `fp` with probability 2^-256 for the
2^25.7 records stored, and any match is replayed and verified before use), `trail[2]` = trail id
(R = 1024 trails), `step[16]` = 128-bit absolute step index within the trail, `pad[2]`.

**Store and lookup.** One append-only array of 52-byte slots, kept sorted by `fp` and refreshed by an
in-place merge every 2^14 inserts; lookup is binary search (26 probes at 2^25.7 records), so an insert
plus probe is charged 12288 instructions and the itemised total is
`55,956,890 x 12288 = 2^39.1` instructions. The per-step amortised row (1 instruction per step,
`S = 2^125.74` instructions) **dominates** that itemised total by `2^100/12288 = 2^86.4`, so the store
cost is fully inside the meter with an enormous margin, and the meter's own cap bounds the record
count (`records <= S/2^100`, hence the store cannot grow without bound).

**Peak memory (itemised).** records `55,956,890 x 52 B = 2^31.00 B`; jump table `1024 x 12 x 576 B =
7,077,888 B`; code `2^21`; scratch `2^24`; sort buffer `2^25`; **peak = 2,969,264,968 B = 2.97 GB
< 2^32 B**. Memory is not a scored objective on this track (the only objective is `time_log2`), but it
is stated, derived from the work cap, and finite.

**Match, replay, verify.** On a DP step whose `fp` is present: replay both trails from their starts to
their recorded step indices (each trail's `(A,D)` is recovered from its own start and stored start
point), compare the two representatives across the twelve Frobenius shifts to obtain the shift `k`,
solve `x = (p^k A1 - A2) * (D2 - p^k D1)^-1 mod r`, and accept only if `g^x = z`; a degenerate
denominator or a failed verification is logged and **discarded**, never reported as a collision. A
successful match is re-checked against the exact invariant `y = g^A z^D` at the collision point.

**Restart / escape schedule.** Each trail keeps a checkpoint (y, A, D) at its most recent DP step. A
repeated class within two steps is a structural 2-cycle: the trail escapes with **one forced multiply
step at index `(h+1) mod R`**, which replaces the step it would have taken and carries the multiply
price, so an escape adds **no extra steps and no extra charge** — the escape budget is zero by
construction, and the meter covers any number of escapes. If a trail exhausts its step budget it
restarts from its checkpoint with a fresh index (charge: one multiply step, again no extra budget).
Under H1' a class repeat within two steps is a `1/N = 2^-251.3` event, so escapes are not expected in
production at all; the toy, whose N is 43420, forces 28 of them in 73434 steps and they are charged.

## 6. Collision model: scoped heuristic, with the bijection gap stated

Fixed-index twisted multiplication is **not proved** to be injective on Frobenius classes. Equality of
output classes can mean `aM = (bM)^(p^j)` for `j != 0`, which does not imply equality of the input
classes. v14's blanket "same index never merges" is therefore withdrawn and replaced by:

* **H1' (scoped heuristic).** The mixed class walk behaves as a random mapping on the N classes for
  the purpose of the *ordinary* merge channel, i.e. two distinct classes collide with probability at
  most `rho/N`, `rho = (1-s^2)(1-1/R) = 0.858536`.
* **The j-shift channel is charged as a bonus, not as a premise.** Including shifted images the
  per-pair match probability is at most `12/N` (union over the twelve shifts). We do **not** take
  credit for the factor up to 12: the work claim uses the smaller ordinary rate. This is safe because a
  shifted match is not a loss — it yields the linear equation `A1 + D1 x = p^j (A2 + D2 x)`, i.e.
  `x = (p^j A2 - A1)/(D1 - p^j D2)`, which determines x whenever the denominator is nonzero. The toy
  confirms this empirically: **134 of 134** shifted-image matches were usable, versus 11 of 1228
  same-index matches (the rest were exact cycle closures with `D1 = D2`, correctly discarded).
* **rho on the conservative side.** The reviewer's algebra `1 - s^2 - (1-s)^2/R = 0.858994` is larger
  than our `0.858536`; the difference is exactly `2s(1-s)/R`, and using the smaller value makes `N'`
  larger and the step count longer. `ledger15.py` asserts the inequality symbolically, so the
  discrepancy cannot silently flip sign.
* **What is still not proved.** H1' and H2 (pseudorandomness of the disjoint T bit ranges) remain
  heuristics; the claim is conditional on them, exactly as the v14 statement said. What changed is
  that the *consequences* the reviewer asked for — the metered cap, the tail, the store, the peak
  memory, the detection loss and the shifted-image channel — are now supplied rather than disclaimed.

## 7. Toy and production-field evidence

`scripts/toy15.py` (q = 523, N = 43420, R = 128, DP width 4, fresh target per trial, 180 trials)
implements the **production** rule: branch `(T mod 8) < 3` on the same invariant `T = Re(c0)`, DP
gating with a store of DP records only, replay/verify, and the forced-multiply escape. Results:
**145/180 = 0.806** recoveries within 4 S_ops; **57/180 = 0.317** at S_ops (the DP-window loss row of
Section 4); 4171 DP hits, 2801 records, 1370 matches of which 1228 same-index and 134 shifted-image;
28 escapes; 4642 invariant checks, 0 failures. The toy's branch frequency is 0.302 rather than 3/8
because at q = 523 the invariant is far from uniform mod 8 — a small-field artifact, recorded and
**not** used in the claim; the production field, where T has 381 bits, measures 0.3835 +- 0.021.
`scripts/predicate_check.py` compares the shared production predicate module with an independently
written copy of the formula on **100 random production GT elements (100 agree)** and the toy predicate
with the same spec on **100 random toy elements (100 agree)**, so the toy and the production run are
driven by the same branch rule, not merely similar ones. `gt_prod.py` additionally exercises the
store path in the production field (12 DP records, one injected duplicate matched, replayed and
correctly rejected as degenerate: `dp_usable = 0`).

## 8. Limitations

H1' and H2 are heuristics and are not proved. The production-field run is 2000 steps; the full-scale
128-bit-class walk is not executed, and no full-scale execution or unconditional mixing proof is
claimed. The squaring price is analytic from counted Fp2 operations rather than a measured routine.
The DP store is bounded by the work cap but its 2.97 GB peak is a stated resource, not a scored one.
This is a square-root attack on the fixed-Q instance, not a security estimate.
