# Six-class parallel Pollard rho on BLS12-381 G1 (negation + GLV)

**Claim.** time_log2 = 126.2037: the algorithm recovers x from T = [x]P0 with
probability >= 0.39, and on every run its total work is <= W_max = 2^126.2037
units (one unit = one batched affine G1 addition = C = 2,000 scalar AArch64
instructions). The work bound is unconditional (global work meter, section 3).
Memory < 2^72 bytes; preprocessing (everything before the first walk step:
jump table, report-table initialisation, the first 2^22 trail starts) <= 2^57.2
units, enforced by the meter and included in W; no advice. Success
rests on the random-mapping heuristic H1 (section 4), supported by toy runs.

## 1. Algebra (all checked on BLS12-381 in `experiment.py`, part A)

p = field prime, r = subgroup order (log2 r = 254.85709), E: y^2 = x^3 + 4.

* beta = 0x1a0111ea...00000000aaac (primitive cube root of 1 in Fp),
  phi(x, y) = (beta x, y) = [lambda] on G1 with
  lambda = 0xac45a4010001a40200000000ffffffff (phi(P0) = [lambda]P0 checked).
* Units u_t = (-1)^s lambda^j (t = 3s + j) act by [u_t](x, y) = (beta^j x, (-1)^s y).
  No nonidentity G1 point is fixed by a unit != 1, so G1 \ {O} has exactly
  N = (r-1)/6 classes of size 6 (log2 N = 252.27213).
* **Parity key** K = y if y is even, else p - y. p is odd, so exactly one of
  y, p - y is even: "take the even one" is invariant under y -> p - y and
  unchanged by phi (which keeps y), so K is constant on a class. Equal keys give
  y' = +-y and x'^3 = x^3, so x' = beta^j x: K is a complete class invariant,
  and bit 0 of K is identically 0.
* **Canonical member**: j minimises beta^j x mod p, s = [y odd] (bit 0 of y);
  c = (beta^j x, K) = [u_t]P. The exact rule costs one Fp multiplication
  (beta x; beta^2 x = p - ((x + beta x) mod p)); the walk needs only the index
  t, which section 1.1 decides exactly without that multiplication.

### 1.1 Exact canonical index from three truncated fixed-point fractions

Only j = argmin(v0, v1, v2) is needed, where v0 = x, v1 = beta x mod p,
v2 = beta^2 x mod p (distinct integers in [1, p)). Let
f_j = floor(2^64 v_j / p) in [0, 2^64), so argmin f_j = argmin v_j (ties between
f_j are the fallback case). The walk forms three 64-bit integer under-estimates
F_j = f_j - e_j and decides j from them, never computing v1, v2 or an inverse.

* **F1 (for v1).** With G = floor(2^448 beta/p) (7 limbs), xG = 2^448 floor(beta x/p)
  + 2^448 v1/p - eps, 0 <= eps < x < 2^381. So the window bits 384..447 of xG equal
  floor(2^64 v1/p - eps/2^384) in {f1 - 1, f1} (eps/2^384 < 1/8). F1 sums only the
  high halves of the six column-5 limb products x_i G_(5-i) and the low halves of the
  six column-6 products x_i G_(6-i), mod 2^64 (6 umulh + 6 madd + 5 add). The dropped
  columns 0..4 and the column-5 low halves contribute a carry into bit 384 of at most
  12, all non-negative, so **0 <= f1 - F1 <= 13** (true bits in {f1-1,f1}, minus a
  carry deficit <= 12).
* **F0 (for v0 = x).** With A = floor(2^448/p) = 9*2^64 + CF (A has 68 bits, so its
  top limb is 9 and CF = A - 9*2^64 < 2^64), and x5 = floor(x/2^320),
  F0 = 9*x5 + floor(x5*CF/2^64) = floor(x5*A/2^64) - {x5 CF/2^64}-rounding
  (1 umulh + 1 madd). Then 2^64 x/p = x5*(2^384/p) + 2^64 xlo/p with xlo < 2^320:
  x5*2^384/p = x5*A/2^64 + x5*eta', 0 <= x5*eta' < 1 (eta' = 2^384/p - A/2^64 < 2^-64),
  the floor of x5*A/2^64 drops a fraction < 1, and 2^64 xlo/p < 2^384/p < 9.6, all
  non-negative, so **0 <= f0 - F0 <= 11**.
* **F2 (for v2), from the sum identity.** 1 + beta + beta^2 = 0 mod p and the v_j
  are nonzero (x != 0; the identity T = O is handled separately and returns 0), so
  v0 + v1 + v2 = k p with k in {1, 2} and f0 + f1 + f2 = 2^64 k - Phi, Phi in {1, 2}.
  Hence f2 = 2^64 k - Phi - f0 - f1, and with F2 = (2^64 - 25) - F0 - F1 mod 2^64
  (1 add + 1 sub) we get, mod 2^64, **f2 - F2 = 25 - Phi - (f0-F0) - (f1-F1)**, i.e.
  e2 = 25 - Phi - e0 - e1 in **[25 - 2 - 11 - 13, 25 - 1 - 0 - 0] = [-1, 24]**.
* **Decision (THR = 26).** Sort F0, F1, F2; let mn, sd be the two smallest and mx the
  largest, jm = argmin. **Fall back** (a separately metered action, 0.25 units) if
  mx >= 2^64 - 26, or mn < 26, or sd - mn < 26; otherwise output jm. Correctness:
  every e_j in [-1, 24], so |e_a - e_b| <= 25 < 26; if sd - mn >= 26 then
  f_sd - f_mn = (sd - mn) + (e_sd - e_mn) >= 26 - 25 = 1 > 0 and no reordering of the
  two smallest is possible, so argmin F = argmin f = argmin v. The two guards remove
  the only wrap cases: F0, F1 are pure under-estimates, so when v0 or v1 is near 0 the
  corresponding F under-wraps to >= 2^64 - 11 (caught by mx >= 2^64 - 26); F2 can
  over-estimate f2 by 1, so when v2 is near p (f2 = 2^64 - 1) F2 wraps to 0 (caught by
  mn < 26), and a genuinely small min (v near 0) is also sent to the exact fallback.
  The gap, min and max tests are one `cmp`/`ccmp`/`ccmn`/`b.hs` chain in the section-5
  listing. The fallback windows cover < 2^7 of the > 2^63 range per fraction, so the
  fallback rate is < 2^-52 per step if beta x mod p is equidistributed (this affects
  only the meter's early-halt probability, never correctness).

Part A compares this rule (full 64-bit limb/flag emulation of the rows-1-3 listing)
with the exact rule j = argmin(v0, v1, v2) on 200,000 random and 150,000 adversarial
near-tie/near-wrap production inputs: 0 mismatches, 0 listing-vs-rule mismatches, and
the measured error ranges e0 in [0, 10.84], e1 in [1.5, 8.2], e2 in [6.9, 22.8] all
lie inside the proved [-1, 24] (none of the 200,000 random inputs fell to the
fallback; the adversarial inputs are built to be near-ties and do).

## 2. Walk

Parameters: R = 2^21; distinguished probability theta = 2^-64 (limb 1 of K is
zero); M = 2^22 walks; L_max = 20 * 2^64; checkpoint period 64.

* **Table.** c_h, d_h uniform in [0, r) (fresh coins, rejection sampling),
  M_h = [c_h]P0 + [d_h]T (redraw if O), Tab[h][t] = [u_t^{-1}]M_h
  = (-1)^s phi^{-j}(M_h) with coefficients u_t^{-1}(c_h, d_h) mod r.
* **Step.** F(C) = class(c(C) + M_h), h = bits 1..21 of K (bit 0 is always 0). The walk keeps any member P
  and sets P <- P + Tab[h][t]; since c + M_h = [u_t](P + Tab[h][t]) this is
  F. Each ordinary step is exactly one batched affine addition (reference
  operation) plus at most 93.0561 instructions (93 listed per step and a share
  below 0.0561 of the per-batch work, section 5).
* **Trails.** Start S = [a]P0 + [b]T, a, b uniform from fresh coins, stored
  in full (redraw if S = O). At each class arrival, in order:

```
if K distinguished: report (K, a, b); end trail
if length >= L_max: abandon trail
if fp(K) == K_cp: traverse the cycle once more (<= 64 steps, else abandon),
     resume at the member of its minimum-fp class remembered during the
     traversal (no extra operation), ESC; clear history; continue
if length mod 64 == 0: K_cp <- fp(K)
P <- P + Tab[h][t]
```

  Cycles are detected only by the checkpoint (the separate 2-cycle test of
  earlier versions is dropped: at R = 2^21 cycles are rare, see section 4).
  "Clear history" invalidates K_cp by writing 0 into its limb-1 word (also 0
  at every trail start). A class that reaches the checkpoint test is not
  distinguished (the distinguished test, limb 1 of K = 0, ends the trail
  first), so its limb 1 is nonzero and a K_cp whose limb-1 word is 0 never
  matches; a refresh stores limb 1 of a non-distinguished K, so a refreshed
  K_cp is valid. No separate valid word is kept.
  ESC is the **doubling escape** P <- [2]P (batched affine doubling). Doubling
  commutes with every unit, so it maps class C to class(2c(C)): a class map.
  fp(K) = limbs 0, 1, 2 of K, i.e. K mod 2^192 (192 bits, 191 free: bit 0 is
  always 0; limb 1 is nonzero for every compared class), stored as K_cp and
  compared as a 192-bit integer (limb 2 most significant, then limb 1, then
  limb 0): the
  checkpoint test and the traversal minimum use this one fixed total order
  on classes (distinct classes have distinct fingerprints except for the
  spurious-match failure event of section 4); full keys are never compared.
  `length` counts every operation (step,
  doubling, traversal step), so a trail has <= L_max + 65 operations. If an addition has equal
  x-coordinates (P = +-Tab), the trail is abandoned (< 2^-114 in total).
* **Merging lemma.** Let trails A and B arrive at the same class C, with any
  histories. Until one of them escapes, both apply F and visit the same
  classes. A trail escapes only after K = K_cp, i.e. a repeat within 64
  ordinary steps (K_cp is set at an arrival and the test precedes the next
  refresh 64 operations later), so only when it is on a cycle Z of F of
  length <= 64 (a spurious fingerprint match is a failure event, section 4);
  Z then contains no distinguished class (the trail would have ended there
  first). **Detection delay:** a trail on Z reaches a length that is a
  multiple of 64 within 63 arrivals, which sets K_cp to a class of Z; it is
  back at that class |Z| arrivals later and the test fires (a K_cp set before
  the trail entered Z is never matched, barring a non-fruitless return). The
  traversal then visits all of Z once and the trail doubles the remembered
  member of min Z (the minimum-fp class of Z), going to
  E(Z) = class(2 c(min Z)), which depends on Z only. The detour takes at most
  63 + 2|Z| + 1 <= 192 operations; fixed points of F (|Z| = 1) take the same
  path. If A escapes from Z, then B, following the same F path, also enters Z,
  detects it within 128 arrivals and escapes to the same E(Z). So A and B are
  again at a common class after at most 192 extra operations of one of them,
  and by induction they reach the same distinguished class, unless one of
  them first hits L_max or an exceptional addition (failure events, section
  4). K_cp is invalidated at every doubling and is invalid at starts, so
  K = K_cp always witnesses a real cycle of F.
* **No structural cycles through escapes.** A trail can return to a class c
  after m >= 1 doublings only if u 2^m c + (sum of unit and power-of-2 multiples
  of table points M_h met on the way) = c for a unit u. The order of 2 mod r is
  2^250.86 (checked in part A from the full factorisation of r - 1), so
  1 - u 2^m != 0 mod r for every m < 2^248: the relation cannot hold by
  cancellation of the table terms (unlike an additive escape multiplier, where
  negation pairs cancel identically); it pins c to one specific point. This
  excludes identities, not coincidences: that such returns occur no more
  often than ordinary self-collisions of a trail is part of heuristic H1
  (section 4). The
  experiment's long-trail test checks this with many escapes per trail
  (section 6).
* **Solving.** On a key match with another trail's report, re-walk both from
  their stored (a, b), tracking P = [A]P0 + [B]T with every value reduced
  mod r after every operation (each step adds the table entry's coefficients;
  a doubling doubles A and B). Canonical coefficients
  (e, f) = u_t (A, B). If the canonical points agree and f_1 != f_2:
  x = (e_2 - e_1)/(f_1 - f_2) mod r (sign checked on production data); verify
  [x]P0 = T and output x. T = O outputs 0.

## 3. Phases and the global work meter

* **Meter.** Every action (random words, sampling attempts, starts and redraws,
  steps, doublings, traversals, reports, hash probes, table build and hash
  initialisation, re-walks, solving, verification) is charged its price
  (section 5) before it runs; each price upper-bounds its instructions / C.
  If a charge would exceed W_max = 2^126.2037 the attack halts with failure.
  So W <= 2^126.2037 on every run; no loop needs its own limit. Preprocessing
  (all charges before the first walk step) has its own meter cap 2^57.2: if a
  preprocessing charge would exceed it, the attack halts with failure. The meter is a
  192-bit fixed-point counter; each of at most 2^22 processors (one per walk)
  reserves blocks of 2^38 units by an atomic update (unused remainders count
  as spent: < 2 * 2^22 * 2^38 = 2^61 units, included in the hard extras).
* **Phase 1** runs until a global operation counter (steps, doublings,
  traversal steps of phase-1 trails; re-walk operations are separate and not
  counted; at most 2^22 processors reserve blocks of 2^38, and every reserved
  block is completed before phase 1 stops) reaches S_ops. Each processor can be
  inside one reserved block of 2^38 operations when the counter crosses S_ops,
  so the overshoot is at most 2^22 * 2^38 = 2^60: at least S_ops
  and fewer than S_ops + 2^60 operations are performed. Ordinary steps are
  priced 1.04652805, traversal steps 1.05902805 (they may copy and compare the
  remembered member and maintain its minimum and counter), doublings 1.49902805; a refresh of K_cp is a separate
  0.0045-unit action (section 5), and a trail makes at most one refresh per 64
  of its operations plus one: refreshes happen only at arrivals whose length is
  a multiple of 64, and a trail's length grows by one per operation, so all
  trails together make at most (operations)/64 + (trails) refreshes.
  **Phase 2** finishes open trails, no new
  trails. Report table capacity D_cap = 2^63.25 and at most 8 processed
  matches (overflow = failure); these bound memory and re-walk work.
* **Batches.** Ready walks (next action a step or doubling) run as one batch of
  2^22 per Montgomery inversion; a smaller batch pays its inversion as a
  separate 118-unit action. In phase 1 every walk is ready at every batch: each
  action a walk branches to from the per-step loop either completes before the
  loop moves on (fallback, report and restart, refresh, abandonment and
  restart, block reservation) or leaves the walk's next action a step or a
  doubling (a traversal, the escape doubling), so a phase-1 batch is the whole
  array of 2^22 walk records. The meter is charged once per batch
  (batch size x price, one 192-bit addition, before the batch runs), and the
  zero-denominator test is done once per batch on the batch product of the
  denominators (zero iff some denominator is zero; then a 2^17-unit recovery
  action scans the batch, at most 2^22 denominators at <= 12 instructions each,
  and abandons the offending walks, and the other walks run as a new batch,
  charged as usual; probability < 2^-114 in total).

## 4. Success probability

N' = N * R/(R-1). T_c = the number of class arrivals of all trails (trail
starts, steps, doublings, traversal steps, and repeated arrivals inside
fruitless cycles: nothing is excluded) up to the first arrival at a class
already visited by another trail. This is exactly what the toy measures.

**H1 (heuristic, at the budget).** P(T_c <= S_ops) >= 1 - exp(-(S_ops/kappa)^2
/(2N')) with kappa = 1.006 (a 0.6% slowdown allowance relative to the
random-mapping law). Toy slowdowns are a harness effect that scales with
(mean trail length)/sqrt(N): a trail's own classes are not useful collision
targets, which delays T_c by about one trail. A trail-length sweep on three
toy curves (section 6) fits slowdown = a + b (trail/sqrt(N)) with
a = +0.24% +- 0.18% (pooled sweep, 168,000 trials); in production
trail/sqrt(N) = 2^64/2^126.1 = 2^-62, so the 0.6% allowance is exactly 2.0
standard errors above the extrapolated value. (Earlier
fixed-trail runs gave -0.2% to +2.5% at trail ratios 0.9%-1.8%; a
no-automorphism control is within 0.6% of the law. The trail-length
explanation and a first sweep are due to submission 11084040, section 6.) Also part of H1: partition
and distinguished bits of K (bits 1..21 and 64..127; the constant bit 0 is used by
neither) are pseudo-random functions of the class; and,
**excluding the fruitless cycles of F of length <= 64** (2-cycles at rate
1/(6R) per step and the rarer short cycles below, which the walk detects and
escapes, and whose repeated arrivals are counted inside T_c), a trail makes a
*non-fruitless return* (an arrival at a class it visited earlier, either
before one of its escapes or along an F-cycle longer than 64) with probability
at most L_max^2/(2N) per trail, the random-mapping self-collision rate
(supported by the long-trail test of section 6; F-cycles of length 65..200
themselves need every partition index repeated, rate
sum_L sum_b S2(L, b) R^-(L-b) < 2^-206 per step, < 2^-137 per trail, and
longer ones are rarer still). Finally, **the canonical twist index of a new
class is uniform on the 6 units and independent of its partition index**
(used for the 2-cycle rate below).

S0 = sqrt(2 N' ln(1/0.6097)) = 2^126.12851; S_ops = 1.006 S0 = 2^126.13714
(H1 gives >= 0.3903; pure model 0.3940); 1.0465983625 (S_ops + 2^60) < 2^126.20285
(1.0465983625 = 1.04652805 per operation plus 0.0045/64 for the refreshes, section 5).

**Phase 1 covers the budget.** Phase 1 performs at least S_ops operations
(hard counter) and each completed operation is one class arrival (an
operation aborted by an exceptional addition is a failure event), so phase 1
makes at least S_ops arrivals: if T_c <= S_ops, the first inter-trail collision happens
in phase 1. No statement about wasted visits is needed: fruitless-cycle
arrivals are counted inside T_c.

**Doubling and traversal cost.** Only the surcharges of doublings
(1.49902805 - 1.04652805 = 0.4525 each) and traversal steps (0.0125 each) are not fixed
by the operation count. A doubling follows a detected cycle:
2-cycles need the next class to have the same partition index and canonical
twist -1, probability 1/(6R) per step; other fruitless cycles of length <= 64
have rate <= sum_L sum_b S2(L, b) R^-(L-b) = 2^-37.999 (every partition index
must occur >= 2 times, which includes identities such as one index thrice with
twists 1, lambda, lambda^2). So the expected number of doublings is
< 2^-23.584 (S_ops + 2^60) (1/(6R) = 2^-23.585; linearity of expectation over the pseudo-random
partition and twist of each new class; no independence across steps is
used), and their expected surcharge is < 0.4525 * 2^-23.584 * 2^126.14 <
2^101.41. The expected number of traversal
steps is < (2/(6R) + 64 * 2^-37.99)(S_ops + 2^60) < 2^103.56; at 0.0125 units
each the expected traversal surcharge is < 2^97.24. The checkpoint-only detection delay (at most
63 + |Z| ordinary steps per cycle) is ordinary work counted in S_ops and its
arrivals are counted in T_c; at R = 2^21 it is < 2^-17.4 of all operations.
(The large table, R = 2^21, is chosen to make cycles rare.)

**Failure events** (union bound with unconditional probabilities; if
T_c <= S_ops, by the merging lemma the two trails reach the same report and
the attack succeeds unless one of these occurs): trail cap hit first (Gamma(2) tail, < 1e-7);
a non-fruitless return inside one of the colliding trails (H1: L_max^2/(2N)
per trail, < 2^-53 in total; detected fruitless cycles are not failures: they
are escaped and the merging lemma covers them);
report overflow (under H1 the distinguished bits of distinct classes are
independent 2^-64 coins, so new-class reports are Binomial(n, 2^-64) with
n <= S_ops + 2^60 + 2^90.4 + starts, mean < 2^62.16; D_cap = 2^63.25 is 2.16
times the mean; Chernoff: < exp(-2^60.9); repeated reports of a class are
matches, at most 8 processed);
f_1 = f_2 (b uniform given the start point: 1/r); exceptional additions
(< 2^-114); 8 processed matches all failing (the attack stops at the first
verified x, so this needs 8 spurious or degenerate matches: < 2^-60); spurious fingerprint match
(fp = limbs 0-2 of K, 192 bits of which 191 are free; both compared classes are
non-distinguished, so under H1 limb 1 of each is uniform on the 2^64 - 1 nonzero
words and two distinct compared classes share a fingerprint with probability
2^-127/(2^64 - 1) < 2^-190.99; union bound over < 2^127.5 fingerprint
comparisons: < 2^-63.49 < 2^-63); meter halting
early (expected surcharges 2^101.41 + 2^97.24 plus other random charges
< 2^76.5, together < 2^101.49, against headroom > 2^115.48: < 2^-13.99 < 0.000062
by Markov, section 5); preprocessing cap (< 2^-32).
**Success >= 0.3903 - 0.000062 - 1.1e-7 > 0.39023 > 0.39.** A doubling and
traversal rate up to 3 times the 1/(6R) value gives a meter bound < 0.00019
and success >= 0.3901 (section 5).

## 5. Cost

**Per ordinary step** (instructions besides the reference addition; blst
Fp multiplication = 383):

| Work | Instr. |
| --- | ---: |
| rows 1-3 (section 1.1), F1: 3 `ldp` (limbs of x; the first is pre-indexed, `ldp x0,x1,[x19,#504]!`, and also advances the walk-record pointer), 6 `umulh` (column 5: x5G0, x4G1, x3G2, x2G3, x1G4, x0G5), 5 `add` (the later column-5 high halves), 6 `madd` (column 6: x5G1..x0G6, each adding a low half into F1) | 20 |
| rows 1-3, F0 and F2: `umulh`, `madd` (F0 = 9 x5 + hi(x5 CF)); `add`, `sub` (F2 = (2^64 - 25) - F0 - F1) | 4 |
| rows 1-3 decision: min/mid/max network (4 `cmp`, 9 `csel`, `cset`, `sub` for a, b, jm, mn, t, mx, sd and the gap; the constant 2 for j is held in x29), then the three-way fallback (`cmp`, `ccmp`, `ccmn`, `b.hs`) | 17 |
| row 4: `ldp` + `ldr` (limbs 0-2 of p from the constants block into x3,x4,x5), `ldp` + `ldr` (limbs 0-2 of y), p - y limbs 0-2 (`subs`, 2 `sbcs`), `ands` (s = bit 0 of y, sets the flags), 3 `csel` (limbs 0-2 of K) | 11 |
| row 5: h, distinguished test, table address (`ubfx` (bits 1..21 of K), `cbz` (limb 1 of K), 3 `add` with shifted operand (8h + 4s + j, then table + 256 times it), `str` of the entry address) | 6 |
| row 6: loads of table point and state, store of state (`ldr` of the entry address, 6 `ldp` table point, 6 `ldp` state, 6 `stp` state) | 19 |
| row 7: checkpoint test (`ldp`, `ldr`, `cmp`, 2 `ccmp`, `b.eq`: 6), refresh test (`tst` length, `b.eq` to the separately priced refresh action: 2) | 8 |
| row 8: length increment (69-bit), cap test, loop (`ldp`, `cmp`, `b.hs`, `adds`, `adc`, `stp`; `cmp`, `b.ne`; the pointer advance is row 1's pre-indexed load) | 8 |
| per-batch work at its share per step: zero-denominator test on the batch product (<= 9), batch meter charge (<= 10), one Fermat inversion (380 squarings + 228 multiplications = 608 x 383 = 232,864) and its exponent scan (< 2,000), per-batch constants and pointers (<= 16): <= 234,898 instructions per batch of 2^22 elements (every phase-1 batch is the whole walk array, section 3), and 234,898 < 0.0561 x 2^22 = 235,323.8 | < 0.0561 |
| **total, a per-step maximum, charged without allowance** | **93.0561** |

The rare fallback of section 1.1 is not in this price (separate action), and
neither is the refresh of K_cp (a separate action, priced below).

Register use (no two resident constants share a register; verified in part A by
emulating rows 1-5 on one register file and counting any write to a resident
register). The resident per-batch constants, set once per batch and never written
in the loop body, are the seven limbs G_0..G_6 (x14-x17, x24-x26), CF (x11), the
constant 9 (x12), the constant 2^64 - 25 (x13), the constant 2 (x29) and the table
base (x21). The pointers x19 (walk-record, advanced by row 1's pre-indexed load),
x20 (batch's last record), x22 (constants block), x23 (meter block) are also held.
**Limbs 0-2 of p are not resident**: row 4 loads them from the constants block into
the temporaries x3, x4, x5 each step, so no constant register is reused for p (an
earlier draft placed CF and 9 in x11, x12 and then reused x11, x12 as p limbs in
row 4, which was a defect; p now lives in temporaries). The G limbs, CF, 9, 2^64 - 25,
the constant 2 and the table base are reloaded only per batch (the "per-batch
constants and pointers" item, <= 16 instructions); the limbs of x (row 1), of p and
of y (row 4) are reloaded every step. x0-x10, x27, x28 are temporaries; x18 and x30
are not used. The twist index jm ends row 3 in x6, which row 4's loads of p (x3,x4,x5),
y (x0,x1,x2), the p - y chain (x7,x8,x9) and s (x28) do not touch, so it survives into
row 5; row 5 reuses x4 and x5 (free after the p - y subtraction) for h and the
address. The decision quantities mn (x3), gap (x5), mx (x4) stay in registers to the
fallback branch. A taken branch (fallback, report, traversal, refresh, abandonment,
block reservation, end of an operation block) leaves the loop for a separately priced
action, which saves and restores what it clobbers inside its own price and returns
into the loop.

Walk record (504 bytes, the largest pre-index offset of `ldp`): x at offset 0, y at
48, the next table-entry address at 96, K_cp at 104 (limbs 0, 1) and 120 (limb 2),
length (two limbs) at 128, the start coefficients a and b at 144 and 176 (32 bytes
each), the remembered member (x, y) of a traversal at 208 and its fingerprint at
304, the traversal counter at 328, and the batched addition's denominator and
prefix product at 336 and 384: 432 bytes used, every offset a multiple of 8 and
every `ldp`/`stp` offset <= 504. x19 starts one record before the batch's first
record and x20 holds the address of its last record, so row 1's pre-indexed load
moves x19 to the current record and row 8's `cmp x19,x20; b.ne loop` ends the loop
after the last one. Table entry (h, t), t = 3s + j, is at table + 256 (8h + 4s + j)
(each entry uses 160 bytes: x, y and the two coefficients; the slots with
4s + j = 3 or 7 are unused). The length load and the cap test of row 8 run before
row 7, whose refresh test reads the length's low limb.

Rows 1-8 as instruction sequences (every operand load and store is listed; `ldp`
and `stp` move two adjacent limbs of a field element in one instruction):

```
row 1 (20)  ldp x0,x1,[x19,#504]!; ldp x2,x3,[x19,#16]; ldp x4,x5,[x19,#32]   (x19 += 504 first; x0..x5 = limbs of x)
            umulh x9,x5,x14; umulh x8,x4,x15; add x9,x9,x8; madd x9,x5,x15,x9
            umulh x8,x3,x16; add x9,x9,x8; madd x9,x4,x16,x9
            umulh x8,x2,x17; add x9,x9,x8; madd x9,x3,x17,x9
            umulh x8,x1,x24; add x9,x9,x8; madd x9,x2,x24,x9
            umulh x8,x0,x25; add x9,x9,x8; madd x9,x1,x25,x9
            madd x9,x0,x26,x9                                    (x9 = F1)
row 2 (4)   umulh x6,x5,x11; madd x6,x5,x12,x6                   (x6 = F0 = 9*x5 + hi(x5*CF))
            add x7,x6,x9; sub x7,x13,x7                          (x7 = F2 = (2^64-25) - F0 - F1)
row 3 (17)  cmp x6,x9; csel x0,x6,x9,lo; csel x1,x9,x6,lo; cset x6,hs   (a=min(F0,F1)->x0, b=max->x1, pair idx->x6)
            cmp x7,x0; csel x3,x7,x0,lo; csel x6,x29,x6,lo       (mn=min(F2,a)->x3; jm=2 if F2<a else pair idx; x29=2)
            cmp x7,x1; csel x8,x7,x1,lo; csel x4,x7,x1,hs        (t=min(F2,b)->x8; mx=max(F2,b)->x4)
            cmp x8,x0; csel x5,x8,x0,hs; sub x5,x5,x3            (sd=max(t,a)->x5; gap=sd-mn->x5)
            cmp x3,#26; ccmp x5,#26,#0,hs; ccmn x4,#26,#2,hs; b.hs fallback
row 4 (11)  ldp x3,x4,[x22,#56]; ldr x5,[x22,#72]                (limbs 0-2 of p into temporaries)
            ldp x0,x1,[x19,#48]; ldr x2,[x19,#64]                (limbs 0-2 of y)
            subs x7,x3,x0; sbcs x8,x4,x1; sbcs x9,x5,x2          (limbs 0-2 of p - y)
            ands x28,x0,#1                                       (x28 = s = [y odd]; ne iff y odd)
            csel x0,x7,x0,ne; csel x1,x8,x1,ne; csel x2,x9,x2,ne (limbs 0-2 of K)
row 5 (6)   ubfx x4,x0,#1,#21; cbz x1,report                     (h = bits 1..21 of K; limb 1 = 0: distinguished)
            add x5,x6,x28,lsl #2; add x5,x5,x4,lsl #3            (8h + 4s + j)
            add x5,x21,x5,lsl #8; str x5,[x19,#96]               (table + 256 (8h + 4s + j))
row 6 (19)  ldr x5,[x19,#96]                                  (entry address)
            6 ldp: table point x, y (12 limbs)         6 ldp: state point x, y (12 limbs)
            6 stp: new state point x, y (12 limbs)
row 7 (8)   ldp x7,x8,[x19,#104]; ldr x9,[x19,#120]              (K_cp limbs 0, 1, 2)
            cmp x0,x7; ccmp x1,x8,#0,eq; ccmp x2,x9,#0,eq; b.eq traversal
            tst x27,#63; b.eq refresh                            (length = 0 mod 64: refresh action)
row 8 (8)   ldp x27,x28,[x19,#128]; cmp x28,#20; b.hs abandon    (length >= 20 * 2^64 = L_max)
 norefresh: adds x27,x27,#1; adc x28,x28,xzr; stp x27,x28,[x19,#128]
            cmp x19,x20; b.ne loop                               (x20 = last record of the batch)
refresh action (9 listed, priced 9):
            ldp x7,x8,[x23]; ldr x9,[x22,#104]; subs x7,x7,x9; sbcs x8,x8,xzr; b.lo reserve
            stp x7,x8,[x23]; stp x0,x1,[x19,#104]; str x2,[x19,#120]; b norefresh
```

The constants block at x22 holds G_0..G_6 at offsets 0-48, limbs 0-2 of p at 56, 64, 72
(loaded fresh each step by row 4, never left resident), CF at 80, the constant 9 at 88,
the constant 2^64 - 25 at 96, and the refresh price at 104; the per-batch init loads CF,
9 and 2^64 - 25 into x11, x12, x13 and G_0..G_6 into x14-x17, x24-x26 once per batch.

Row 1's 12 products are all present (six `umulh` for the column-5 high halves,
six `madd` for the column-6 low halves), and its sum of the 12 terms uses 11
additions (five `add` plus six accumulated by `madd`), one fewer than a sum of
12 terms needs to start from zero. **The three-way fallback test.** After row 3's
network x3 = mn, x5 = gap = sd - mn and x4 = mx. `cmp x3,#26` leaves `lo` iff
mn < 26. `ccmp x5,#26,#0,hs` runs `cmp x5,#26` only when `hs` (mn >= 26), leaving
`lo` iff gap < 26, and otherwise forces NZCV = 0000 (C = 0, `lo`); so after it `lo`
holds iff (mn < 26 or gap < 26). `ccmn x4,#26,#2,hs` runs `cmn x4,#26` only when
`hs` (i.e. the previous two tests both failed), setting C iff mx + 26 >= 2^64 i.e.
mx >= 2^64 - 26, and otherwise forces NZCV = 0010 (C = 1, `hs`); so `b.hs fallback`
is taken iff (mn < 26 or gap < 26 or mx >= 2^64 - 26), exactly the three section-1.1
fallback conditions, and the immediate 26 fits the 5-bit field of `ccmp`/`ccmn`.
**Part A runs rows 1-5 as one composed leaf loop body on a single 64-bit register
file** with NZCV flags, so the index rows (1-3) and the key/address rows (4-5) are
validated together, not in isolation. It compares the branch, the emitted index jm,
the key limbs 0-2 of K, the sign s, the partition h, the distinguished test and the
table address 8h + 4s + j with the exact rules, on 200,000 random and 150,000
adversarial inputs (and a mix of random, boundary and distinguished-key y): 0
mismatches. It also counts every write to a resident-constant register (G_0..G_6, CF,
9, 2^64 - 25, the constant 2, the table base) and asserts those registers are
distinct: **0 clobbers**, confirming the register allocation above is conflict-free
(in particular, p is loaded into x3,x4,x5 in row 4 and never reuses CF's x11 or 9's
x12). The key row is exact: p is odd, so y and p - y have opposite parities, K is
p - y exactly when y is odd, and s = [y odd] is bit 0 of y (`ands`, which also sets
the flags the three `csel` read). A borrow propagates only upward, so the three
`subs`/`sbcs` give limbs 0-2 of p - y exactly and limbs 3-5 of y are not loaded; the
row has no branch and no fallback. Only limbs 0-2 of K are selected, because a step
uses only those (h from bits 1..21 of limb 0, the distinguished test from limb 1, the
fingerprint from limbs 0-2); the report action recomputes all six limbs inside its
0.5-unit price. The checkpoint test needs no valid word: a K_cp whose limb-1 word
is 0 cannot equal limbs 0-2 of K, whose limb 1 is nonzero after row 5's `cbz`.
A separate boundary-stressed check (`check_parity_rule`) emulates rows 4 and 5 limb
by limb against the exact K mod 2^192, s, h, the distinguished test and the table
address on random, boundary and distinguished-key y, and verifies the address map is
injective into the 2^24 table slots. Row 6 keeps the scope of the earlier 36-transfer
row (one load of each limb of the table point and of the state, one store of each
limb of the new state) and pairs adjacent limbs; it is the reference addition's
own operand traffic, so it is shown as counts of `ldp`/`stp` rather than as a
register schedule that would collide with the addition's six multiplications.
Its 36 scalar transfers become 19 instructions. **The refresh action** charges its
own price first (the processor's meter-block remainder at [x23], 128-bit fixed
point, less the refresh price rounded up; a borrow branches to the separately
priced block reservation), then stores limbs 0-2 of K as K_cp; since limb 1 of K is
nonzero here, the stored K_cp is valid. It is taken only at arrivals whose length
is a multiple of 64 (`tst x27,#63` reads the length before the increment), as in
the walk of section 2.

**Prices** (units): ordinary step 1 + 93.0561/2000 = 1.04652805 (the per-step
maximum, charged without an allowance); refresh of K_cp 0.0045 (9 instructions = the 9 listed: the fall-through path is exactly the 9 listed instructions and the borrow path leaves to the separately priced block reservation); traversal step 1.05902805 (93.0561 + 25 instructions for the
traversal action: copying the remembered member when it changes (6 `ldp` + 6
`stp` = 12), loading the running minimum fingerprint (`ldp` + `ldr` = 2),
comparing it (`cmp`, 2 `ccmp`, `b.hs` = 4), storing the new minimum on the
taken path (`stp` + `str` = 2), and maintaining and testing the traversal
counter (`ldr`, `add`, `str`, `cmp`, `b.eq` = 5): 118.0561 instructions =
0.05902805 x 2000; this is the maximum path, charged without claiming savings
from omitted ordinary-loop instructions); doubling 1.49902805 (7 x 383 = 2,681 + 7 x
32 = 224 + 93.0561 = 2,998.0561 instructions = 1.49902805 x 2000; a doubling adds no
table point, so the six table-point `ldp` of row 6 are not executed and cover its
`str xzr,[x19,#112]`, which writes 0 into the limb-1 word of K_cp); section-1.1
fallback 0.25 (383 + 32 + 15 + 24 = 454 <= 500 instr.);
random 64-bit word (RNDR) 0.005; sampling attempt 0.025; trail start
[a]P0 + [b]T (two Jacobian scalar multiplications + inversion, about 4.8e6
instr., including the `str xzr` that leaves K_cp invalid) 2^12; report 0.5; hash
probe 0.015; hash initialisation 2^57.1 (once); table 2^33.2 (once: 2^21
entries, plus its sampling); phase-2 or re-walk operation (unbatched, incl.
inversion, coefficient updates and a possible refresh) 119 (232,864 + 2,000 +
1,149 + 192 + 93.0561 + 9 + 100 = 236,414.0561 <= 238,000); small-batch inversion
118 (232,864 + 2,000 = 234,864 <= 236,000); solving 2^9; verification 2^13.

**Budget.** Phase-1 operations and their refreshes: every operation is charged
1.04652805 plus, for doublings and traversal steps, its surcharge, and the refreshes
number at most (operations)/64 + (trails) (section 3). The deterministic part is
(1.04652805 + 0.0045/64) (S_ops + 2^60) = 1.0465983625 (S_ops + 2^60) < 2^126.20285; the
one-per-trail refreshes are counted with the trail starts below, and the doubling
and traversal surcharges are random (expectation < 2^101.41 + 2^97.24).
Hard extras: meter-block remainders < 2^61, block reservations (an atomic update
of a global counter, < 64 instructions, at most once per 2^38 units or 2^38
operations of a processor: < 2^89.26 of them, < 2^84.3 units), phase 2
(2^22 (L_max + 65) * 119 = 2^97.217), re-walks (8 * 2 * ((L_max + 65) * 119 +
2^12) = 2^79.217), reports (D_cap * 0.5 = 2^62.25), hash init 2^57.1, table and
verification < 2^33.3: together < 2^97.3 (2^97.2169). Random extras (expected):
starts and their first refreshes < 2^64.4 * (2^12 + 0.0045) < 2^76.49,
words/attempts/probes < 2^63, small batches < 2^8 (a phase-1 batch holds all 2^22
walks), fallbacks < 2^72.3, zero-batch recoveries and their wasted batch charges
< 2^-114 (2^17 + 2^23) < 2^-90: together < 2^76.5 (2^76.4819). The random charges
are together < 2^101.41 + 2^97.24 + 2^76.5 < 2^101.49. Headroom
2^126.2037 - 2^126.20285 - 2^97.3 > 2^115.48, so the meter halts early with probability
< 2^101.49/2^115.48 = 2^-13.99 < 0.000062
(Markov), and success is at least 0.3903 - 0.000062 - 1.1e-7 > 0.39023. With a doubling and traversal rate up to 3 times the 1/(6R) value
the bound is (3 * 2^101.41 + 3 * 2^97.24 + 2^76.5)/2^115.48 < 0.00019 and the
success is still >= 0.3903 - 0.00019 - 1.1e-7 > 0.3901.

**Preprocessing:** jump table (priced 2^33.2, plus sampling),
report-table initialisation 2^57.1, first 2^22 starts (2^34, plus sampling):
expected < 2^57.11, cap 2^57.2 (the expected random sampling part is < 2^21,
far below the cap gap, so the cap halts a run with probability < 2^-32,
Markov); all of it is inside W.

**W <= 2^126.2037 on every run.** Memory: D_cap records x 112 bytes at load 1/2
= 2^71.06, plus table (2^21 x 8 slots x 256 bytes = 2^32 bytes) and walk states
(2^22 x 504 bytes = 2^31 bytes): < 2^72.

## 6. Experiments and limitations

`experiment.py` (stdlib). **Part A** (production BLS12-381; all true): curve,
order, beta/lambda pairing, orbit and unit action, parity-key invariance (one
key on all six members, bit 0 zero, a bijection with |y| = min(y, p - y), so the
same classes), the boundary-stressed row-4/row-5 limb rule against the exact K
(limbs 0-2), s, h and distinguished test (40,000
random, 4,013 boundary and 1,000 distinguished-key y; with limbs 0-2, the
`cbz` and the table address table + 256 (8h + 4s + j) for every j; 0 mismatches),
**the composed rows 1-5 run instruction by instruction on one register file**, with
the branch, the emitted index jm, the key K, s, h, the distinguished test and the
table address all checked against the exact rules, and every write to a resident
constant register counted (200,000 random and 150,000 adversarial near-tie/near-wrap
inputs, 0 mismatches of any output, 0 resident-register clobbers, residents distinct;
the measured error ranges e0 in [0, 10.84], e1 in
[1.5, 8.2], e2 in [6.9, 22.8] all inside the proved [-1, 24]), unique canonical member, twisted-step class identity, doubling class map, order of 2
mod r, table
entries, coefficient tracking on production steps, and recovery of a planted
x from two production trails forced to merge. **Part B** (toy y^2 = x^3 + 4
over F_q, q = 58493047 or 861101287): the per-trail algorithm of section 2
(parity key with h and the distinguished test read above bit 0, exact canonical
rule, which the section-1.1 three-fraction rule always reproduces on production inputs, escapes, caps, restarts, re-walk solving). Default mode
reports T_c; `"pipeline": true` also applies the scaled production rules
(phase-1 budget 1.006 S0 class arrivals, phase 2, report capacity, 8-match
cap; the toy capacity is about 2.7 times the mean instead of 2.06, i.e.
looser) and reports success. Trails run sequentially (same visit distribution
under H1); [a]G + [b]T is evaluated as [a + bx]G for speed (same point); toy
keys are below 2^64, so the fingerprint (limbs 0, 1, 2) of a toy key is the key
itself and the toy's equality tests and minima follow the production rule
exactly. A toy key has limb 1 = 0 (its distinguished test reads the bits above h
instead), so the toy marks an invalid K_cp with an explicit flag; the rule is the
same as in production: an invalid K_cp never matches.

The sweep and scaled-pipeline figures below were measured with the toy walk keyed
by |y| = min(y, q - y) (h = |y| mod R, distinguished test on the next bits), which
`"class_key": "abs"` reruns exactly. The two keys induce the same classes and
determine each other (K = |y| if |y| is even, else q - |y|), so they differ only
in which pseudo-random label each class carries; H1 assumes only that these labels
are pseudo-random. The declared run uses the parity key.

**Trail-length sweep** (checkpoint-only walk, R = 2048, T_c mode, seeds
12001-12028 and 14001-14021; every recovered x verified; toy curve q34:
q = 10146855109, r' = 3382230043, N = 2^29.07; S = kappa sqrt(2 N' ln(1/0.6097))):

| toy | dp_bits | mean trail / sqrt(N) | trials | mean T_c vs model | P[T_c <= S], kappa = 1.006 |
| --- | ---: | ---: | ---: | --- | ---: |
| q34 | 2 | 0.017% | 12000 | +0.67% +- 0.47% | 0.392 +- 0.005 |
| q34 | 4 | 0.067% | 39000 | +0.40% +- 0.27% | 0.397 +- 0.003 |
| q34 | 5 | 0.135% | 21000 | +0.38% +- 0.36% | 0.396 +- 0.003 |
| q34 | 6 | 0.27% | 51000 | +0.11% +- 0.23% | 0.400 +- 0.002 |
| q30 | 4 | 0.23% | 24000 | +0.60% +- 0.34% | 0.397 +- 0.003 |
| q30 | 6 | 0.93% | 24000 | +1.18% +- 0.34% | 0.391 +- 0.003 |

Weighted fit slowdown = a + b (trail/sqrt(N)) over these 171,000 trials:
a = +0.24% +- 0.18%, b = 0.86 +- 0.46 (consistent with 11084040's q26/q30
sweep, a = +0.12% +- 0.16%, b = 1.10 +- 0.10). Every row is above the H1
floor 0.3903 at kappa = 1.006 (the lowest, 0.391, is the largest trail ratio,
as the trail law predicts). **Scaled full pipeline at kappa = 1.006** (budget
1.006 S0 class arrivals, phase 2, report capacity, 8-match cap; seeds
14022-14029 and 15001-15014): q34 dp 6 0.406 +- 0.004 (16,000 trials), q26
dp 3 0.392 +- 0.004 (18,000), q30 dp 4 0.392 +- 0.003 (24,000); pooled
0.396 +- 0.002 over 58,000 trials against the claimed 0.39. Earlier fixed-trail
runs of the same walk (dp_bits 5-6, trail ratios 0.9%-1.8%) gave +0.2% to
+2.5%, as the trail law predicts; tables R <= 256 are not used (their
checkpoint delay shows in the toy; at R = 2^21 it is < 2^-17.4). The schema
keeps toy trails short (dp_bits <= 6); a 14-second run has standard error
~0.02-0.05.

**Long-trail trap test** (`long_trails`; q30, N = 2^25.5, walks without
distinguished points, 150 trails of 2,048 operations each; repeats inside a
fruitless cycle of length <= 64 (same escape epoch, <= 64 arrivals after the
previous visit) are counted separately; every other repeat is a non-fruitless return;
ordinary self-collisions predict non-fruitless returns in about
1 - exp(-L'^2/(2N)) of trails, L' = new-class arrivals per trail):

| R | escape | escapes per trail | trails with non-fruitless returns |
| ---: | --- | ---: | ---: |
| 4 | doubling (the attack) | 29 | 1 / 150 |
| 8 | doubling (the attack) | 23 | 1 / 150 |
| 4 | additive M_esc (control) | 30 | 46 / 150 |
| 8 | additive M_esc (control) | 25 | 14 / 150 |

(With checkpoint-only detection most operations at R = 4 and 8 are cycle
repeats, so a trail makes only about 670 and 1,090 new-class arrivals; the
self-collision expectation from those is about 0.7 and 1.9 trails of 150.)

Doubling escapes behave like ordinary self-collisions even with many escapes
per trail; the additive escape traps trails structurally, as expected.

Limitations: H1 is a heuristic (random mapping with Brent-Pollard correction;
van Oorschot-Wiener; automorphism classes as in Wiener-Zuccherato, GLV,
Duursma-Gaudry-Morain; fruitless-cycle escapes as in Bos-Kleinjung-Lenstra and
Bernstein-Lange-Schwabe), supported only at toy scale; instruction counts are
analytic; generic attack (no TNFS).
