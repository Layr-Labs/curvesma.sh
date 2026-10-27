# Lockstep six-class Pollard rho on BLS12-381 G1 (negation + GLV)

**Claim.** time_log2 = 126.1868: the algorithm recovers x from T = [x]P0 with
probability >= 0.39, and on every run its total work is <= W_max = 2^126.1868
units (one unit = one batched affine G1 addition = C = 2,000 scalar AArch64
instructions). The work bound is unconditional (global work meter, section 3).
Memory < 2^62.3 bytes (every array is below 2^63 bytes, so all addressing is
ordinary 64-bit); preprocessing (everything before the first walk step:
jump table, report-table initialisation, the first 2^22 trail starts)
<= 2^56 units, enforced by the meter and included in W; no advice. Success
rests on the random-mapping heuristic H1 (section 4), supported by toy runs.

The attack is one sequential program (one processor); parallel execution would
change elapsed time only, and the score sums work. It is the six-class walk of
the earlier v21/v31 candidates (submission 9b9712fa; same classes, parity key,
twisted table, doubling escape, distinguished points, re-walk solving) with
five changes, each removing listed per-step work rather than re-pricing it:

1. the canonical twist index is decided in 64-bit fraction arithmetic
   (section 1.1): rows 2-3 are 19 instructions instead of 26;
2. the 7 limbs of G stay in registers for the whole batch (row 1: 20 instead of 24);
3. all walks run in lockstep batches with one global batch index t, so the
   checkpoint refresh and the trail cap are scheduled by t (refresh batches)
   and no per-walk length counter or per-step refresh test exists;
4. R = 2^40 makes fruitless cycles rare (2^-42.585 per step), so the checkpoint
   is compared only in test batches (every 256th) and a detected cycle is
   traversed and escaped synchronously by a separately metered action;
5. ordinary batches need only limbs 0-1 of the class key.

An ordinary step costs one batched affine addition plus at most 72.0561 listed
instructions (section 5); test and refresh batches add 3 and 12 instructions
per walk.

## 1. Algebra (all checked on BLS12-381 in `experiment.py`, part A)

p = field prime, r = subgroup order (log2 r = 254.85709), E: y^2 = x^3 + 4.

* beta = 0x1a0111ea...00000000aaac (primitive cube root of 1 in Fp),
  phi(x, y) = (beta x, y) = [lambda] on G1 with
  lambda = 0xac45a4010001a40200000000ffffffff (phi(P0) = [lambda]P0 checked).
* Units u_t = (-1)^s lambda^j (t = 3s + j) act by [u_t](x, y) = (beta^j x, (-1)^s y).
  No nonidentity G1 point is fixed by a unit != 1, so G1 \ {O} has exactly
  N = (r-1)/6 classes of size 6 (log2 N = 252.27213).
* **Parity key** K = y if y is even, else p - y. p is odd, so exactly one of
  y, p - y is even: K is invariant under y -> p - y and unchanged by phi (which
  keeps y), so K is constant on a class. Equal keys give y' = +-y and
  x'^3 = x^3, so x' = beta^j x: K is a complete class invariant, and bit 0 of K
  is identically 0.
* **Canonical member**: j minimises v_j = beta^j x mod p, s = [y odd] (bit 0 of
  y); c = (beta^j x, K) = [u_t]P. The exact rule costs one Fp multiplication
  (beta x; beta^2 x = p - ((x + beta x) mod p)); the walk needs only the index
  t, which section 1.1 decides exactly without that multiplication.

### 1.1 Exact canonical index in 64-bit fraction arithmetic

Only j = argmin(v0, v1, v2) is needed (v0 = x, v1 = beta x mod p,
v2 = beta^2 x mod p, distinct integers in [1, p)). Write F_j = 2^64 v_j / p
(reals in (0, 2^64)). Because x + beta x + beta^2 x = (1 + beta + beta^2) x = 0
mod p, v0 + v1 + v2 is p or 2p, so **F0 + F1 + F2 is exactly 2^64 or 2^65**:
F2 is determined by F0 and F1 modulo 2^64. The rule computes 64-bit estimates:

* **F1 = L'.** With G = floor(2^448 beta/p) (7 limbs), frac(xG/2^448) = v1/p - e,
  0 <= e < 2^-67. L' = (sum of the high halves of the six column-5 limb
  products x_i G_(5-i) and the low halves of the six column-6 products
  x_i G_(6-i)) mod 2^64: 12 mul/umulh. Everything dropped from xG mod 2^448
  (columns 0..4, column-5 low halves) is < 12 * 2^384, so
  2^64 frac(xG/2^448) lies in [L', L' + 12) modulo 2^64, and
  **F1 - L' lies in [0, 12.125)** (modulo 2^64).
* **F0 = 9 X + umulh(X, Dfrac)**, X = limb 5 of x, where
  floor(2^448/p) = 9 * 2^64 + Dfrac (Dfrac = 0xd835d2f3cc9e45ce; 2^384/p =
  2^3.2993 = 9.85). With x = X 2^320 + rho, 0 <= rho < 2^320: the estimate is
  <= X 2^384/p <= F0 (floors only drop) and > X 2^384/p - 2 > F0 - 9.85 - 2,
  so **F0 - est0 lies in [0, 11.85)**; it never wraps (0 <= est0 < 2^64).
* **F2 = -(F0 + L') mod 2^64.** Since F0 + F1 + F2 = 0 mod 2^64,
  est2 - F2 = (F0 - est0) + (F1 - L') = e2 lies in [0, 24) (mod 2^64).
* **Tests.** (A) If L' >= 2^64 - 14, fall back (only then can L' have wrapped
  below 0, i.e. F1 < 12.125). (B) If est2 < 31, fall back (only then can est2
  have wrapped above 2^64, i.e. F2 > 2^64 - 24; 31 is the largest `ccmp`
  immediate). Otherwise the three estimates
  are the true real values minus errors d0 in [0, 11.85), d1 in [0, 12.125),
  and plus e2 = d0 + d1. Any two errors differ by less than
  2 * 12.125 + 11.85 < 37 (the largest spread is est2 against est1:
  d0 + 2 d1). (C) If the smallest estimate is below both others by at least
  41, it is the exact argmin (a gap >= 41 > 37 cannot be reversed); otherwise
  fall back.
* **Fallback** recomputes j exactly with one Fp multiplication: a separately
  metered action (0.25 units). The fallback windows (A: 14, B: 31, C: a
  pair of true values within 41 + 37 = 78 of each other, < 3 * 2 * 78 values) cover < 2^9 values out
  of 2^64, so the fallback rate is < 2^-55 per step if beta x mod p is
  equidistributed (this affects only the meter's early-halt probability, never
  correctness).

Part A compares this rule and the instruction-level listing of rows 1-3 with
the exact rule on 24,000 random and 8,000 adversarial near-tie production
inputs (0 mismatches; 0 mismatches of the merged fallback branch against the
OR of tests A-C; observed error ranges d0 <= 10.9, d1 <= 8.2, e2 <= 18.3,
inside the bounds above; offline larger runs in section 6).

## 2. Walk

Parameters: R = 2^40 partition indices; distinguished probability
theta = 2^-73 (limb 1 of K is zero and bits 55..63 of limb 0 are zero:
the loop tests limb 1 with one `cbz`, and the rare taken branch (2^-64 per
step) tests the 9 further bits in a separately priced `dpcheck` action);
M = 2^22 walks; trail cap L_max = 20 * 2^73 operations; test period k = 256 batches; refresh period
P = 2^14 = 64 k batches (every refresh batch is also a test batch); L_max is a
multiple of P.

* **Table.** c_h, d_h uniform in [0, r) (fresh coins, rejection sampling),
  M_h = [c_h]P0 + [d_h]T (redraw if O), Tab[h][t] = [u_t^{-1}]M_h
  = (-1)^s phi^{-j}(M_h) with coefficients u_t^{-1}(c_h, d_h) mod r, for
  h < 2^40 and t < 6. Entry (h, t), t = 3s + j, is at table + 256 (8h + 4s + j)
  (160 bytes used per 256-byte slot: x, y and the two coefficients; slots with
  4s + j = 3 or 7 unused): 2^43 slots, 2^51 bytes.
* **Step.** F(C) = class(c(C) + M_h), h = bits 1..40 of K (bit 0 is always 0).
  The walk keeps any member P and sets P <- P + Tab[h][t]; since
  c + M_h = [u_t](P + Tab[h][t]), this is F.
* **Lockstep batches.** The program keeps the M walk records in one array and
  processes them in batches indexed by a global 128-bit batch counter t. In
  every batch each walk performs exactly one batched affine addition (its
  step). Before the addition, the canonicalisation loop of section 5 processes
  the walk's current point (its *arrival*) in this order:

```
if K distinguished:            report (K, a, b, t_start); start a new trail; redo this arrival
if t = 0 mod k  (test batch):  if limb 0 of K == limb 0 of K_cp: detect action (below)
if t = 0 mod P  (refresh):     if t >= expiry: abandon; start a new trail; redo this arrival
                               K_cp <- limbs 0-2 of K
store the table-entry address for the batched addition (P <- P + Tab[h][t])
```

  "Redo this arrival" means the action ends by moving the record pointer back
  one record, so the loop body runs again on the same record (for the new
  point) in the same batch; that re-run is part of the action's price. A trail
  started in batch t_start stores a, b and expiry = t_start + L_max (2 limbs);
  a trail therefore makes at most L_max + P operations. K_cp is *invalid* at
  trail starts and after escapes: its limb 0 is set to 1, which never equals
  limb 0 of a key (bit 0 of K is 0).
* **Detect action.** It first compares limbs 1-2 of K with limbs 1-2 of K_cp
  (limb 2 of K recomputed); if they differ, it was a spurious 64-bit match and
  the loop continues with the refresh test. Otherwise (a 192-bit match) the
  walk's class has repeated since the refresh batch that stored K_cp, with no
  escape in between, so the walk is on a cycle Z of F of length dividing a
  number <= P. The action **traverses Z synchronously**: it applies F with
  unbatched operations (each with its own inversion, priced 119 units, section
  5) from the current class until it returns to that class (full-key
  comparison), at most P operations (more is impossible except after a
  spurious 192-bit match; the trail is then abandoned), remembering the
  member of the class of Z with the smallest full key K. It then **escapes** by
  doubling that member (unbatched, 119 units), sets K_cp invalid, and redoes
  the arrival for the escape point. Doubling commutes with every unit, so the
  escape goes to E(Z) = class(2 c(min Z)), which depends on Z only.
* **Detection delay.** If a walk is on a cycle Z with |Z| <= 64, it is detected
  within 2P batches: the next refresh batch t0 stores K_cp = a class of Z;
  at the test batches t0 + 256 m the walk is 256 m steps further along Z, so
  the test matches for m = |Z|/gcd(|Z|, 256) <= 64, i.e. no later than batch
  t0 + P, whose test precedes its refresh. A trail on a cycle longer than 64
  may stay undetected; such a trail is abandoned at its cap. A trail that
  re-enters a class it visited earlier is necessarily on a cycle of F (F is a
  function on classes and no escape intervened), so a 192-bit match never
  arises from a transit class.
* **Merging lemma.** Let trails A and B arrive at the same class C, at any
  batches. Each applies F once per batch, so both visit the same classes until
  one of them escapes. A trail escapes only after a detected repeat, i.e. on a
  cycle Z of F; then the other, following the same F path, also enters Z. If
  |Z| <= 64 both detect Z within 2P batches of entering it and both escape to
  the same class E(Z). So A and B are again at a common class, and by
  induction they reach the same distinguished class, unless one of them first
  hits its cap, an exceptional addition, or a fruitless cycle longer than 64
  (failure events, section 4). Z contains no distinguished class (the trail
  would have ended there).
* **No structural cycles through escapes.** A trail can return to a class c
  after m >= 1 doublings only if u 2^m c + (sum of unit and power-of-2
  multiples of table points M_h met on the way) = c for a unit u. The order of
  2 mod r is 2^250.86 (checked in part A from the full factorisation of r - 1),
  so 1 - u 2^m != 0 mod r for every m < 2^248: the relation cannot hold by
  cancellation of the table terms; it pins c to one specific point. This
  excludes identities, not coincidences: that such returns occur no more often
  than ordinary self-collisions of a trail is part of heuristic H1 (section 4).
* **Solving by replay.** On a key match with an earlier report, re-walk both
  trails from their stored (a, b, t_start), executing the same per-walk rules
  with batch indices t_start, t_start + 1, ... (the rules of a walk depend only
  on its own state and on t, and the canonical index is exact, so the replay
  repeats the original trail operation by operation, including its detections
  and escapes), tracking P = [A]P0 + [B]T with every value reduced mod r after
  every operation (each step adds the table entry's coefficients, a doubling
  doubles A and B). Canonical coefficients (e, f) = u_t (A, B). If the
  canonical points agree and f_1 != f_2: x = (e_2 - e_1)/(f_1 - f_2) mod r
  (sign checked on production data); verify [x]P0 = T and output x. T = O
  outputs 0.
* **Reports: log and index.** A report appends the 128-byte record
  (K: 48 bytes, a, b: 32 each, t_start: 16) to an append-only **report log**
  of capacity D_cap = 2^54.35 records (2^61.35 bytes; a full log is a failure
  event) and inserts (fingerprint, log index) into the **report index**, an
  open-addressing table of m = 2^57 slots of 16 bytes (2^61 bytes) with
  linear probing: home slot i = bits 0..56 of limb 2 of K, fingerprint =
  limb 3 of K (empty slot: index field all ones). A lookup probes slots i,
  i + 1, ... (mod m) until an empty slot (no match: insert there) or a slot
  whose fingerprint equals limb 3 of K, in which case the log record is read
  and its full K compared (equal: a match; different: continue probing).
  At most 256 probes are made; reaching that cap is a failure event. The load
  is at most 2^54.35/2^57 = 0.159. Under H1 the index bits of distinct classes
  are uniform and independent, so a lookup needs more than k probes only if
  some run of L >= k consecutive slots containing its home slot receives
  >= L keys; by the Chernoff bound P(Bin(n, L/m) >= L) <= (e^(1-a) a)^L with
  a = 0.159, e^(1-a) a = 0.369, so P(more than 256 probes) <= sum_{L>=256}
  L 0.369^L < 2^-359 per lookup, and E[probes] <= 1 + sum_L L^2 0.369^L <
  3.02. Repeated reports of a class are matches (the class is stored once);
  at most 8 matches are processed. Every array (table 2^51 bytes, index 2^61,
  log 2^61.35, walk records 2^30.3) is below 2^63 bytes, so all addresses are
  ordinary 64-bit addresses.

## 3. Phases and the global work meter

* **Meter.** Every action (random words, sampling attempts, table build, index
  initialisation, trail starts and redraws, batches, dpcheck branches, reports, index probes,
  fallbacks, detect actions and each of their unbatched operations, phase-2
  and re-walk operations, solving, verification) is charged its price (section
  5) before it runs; each price upper-bounds its instructions / C. A batch is
  charged once, before it runs, as (number of walks) x (the step price of its
  variant). If a charge would exceed W_max = 2^126.1868 the attack halts with
  failure. So W <= 2^126.1868 on every run. Preprocessing (all charges before
  the first walk step) has its own cap 2^56: if a preprocessing charge would
  exceed it, the attack halts with failure. The meter is a 192-bit fixed-point
  counter in memory, updated by the one sequential program; there are no
  per-processor reservations.
* **Phase 1** runs B1 = ceil(S_ops / M) full batches, so it performs
  B1 M operations, at least S_ops and fewer than S_ops + 2^22. Every phase-1
  batch is the whole array of M walk records: every action that a walk
  branches to from the loop (fallback, report and restart, detect with its
  synchronous traversal and escape, abandonment and restart) completes before
  the loop moves on and leaves the walk with a current point whose arrival has
  been processed, so each walk performs exactly one batched addition per batch.
* **Batch variants.** The loop body exists in three straight-line variants:
  ordinary (section 5 listing), test (t = 0 mod 256: the ordinary body plus
  `ldr x2,[x19,#104]; cmp x0,x2; b.eq detect`, 3 instructions per walk) and
  refresh (t = 0 mod 2^14: row 4 replaced by row 4', which also produces limb 2
  of K (+3), the test (+3), the cap test `ldp x2,x3,[x19,#128]; cmp x13,x2;
  sbcs xzr,x16,x3; b.hs abandon` against the 128-bit batch index held in
  x13:x16 (+4), and `stp x0,x1,[x19,#104]; str x7,[x19,#120]` storing K_cp
  (+2): 12 instructions per walk). The variant is chosen once per batch from t
  (the per-batch dispatch item of section 5).
* **Phase 2** finishes open trails (no new trails), continuing the same rules
  and batch indices; every phase-2 operation, batched or not, is charged the
  unbatched price 119. **Report log** capacity D_cap = 2^54.35, at most 256
  index probes per lookup and at most 8 processed matches (overflow = failure);
  these bound memory, lookup work and re-walk work.
* **Batches.** One Montgomery inversion per batch of 2^22 walks; the
  zero-denominator test is done once per batch on the batch product of the
  denominators (zero iff some denominator is zero; then a 2^17-unit recovery
  action scans the batch, at most 2^22 denominators at <= 12 instructions each,
  and abandons and restarts the offending walks (restart charged), and the
  batch runs again, charged again; probability < 2^-114 in total).

## 4. Success probability

N' = N * R/(R-1). T_c = the number of class arrivals of all trails (trail
starts, batched steps, the detect action's traversal steps and escape
doublings, every repeated arrival: nothing is excluded) up to the first
arrival at a class already visited by another trail. This is exactly what the
toy measures.

**H1 (heuristic, at the budget).** P(T_c <= S_ops) >= 1 - exp(-(S_ops/kappa)^2
/(2N')) with kappa = 1.005 (a 0.5% slowdown allowance relative to the
random-mapping law). Toy slowdowns are dominated by a harness effect that
scales with (mean trail length)/sqrt(N): a trail's own classes are not useful
collision targets. In production trail/sqrt(N) = 2^73/2^126.1 = 2^-53; the
remaining production effects are the detection-delay repeats below
(a fraction < 2^-27.5 of all operations) and the r-adding-walk correction,
which N' already contains. The toy evidence (section 6) is a sweep over trail
ratios on three toy curves with the lockstep production schedule, fitted as
slowdown = a + b (trail/sqrt(N)), and the earlier v31 sweep of the same class
walk (a = +0.24% +- 0.18%). Also part of H1: the partition and distinguished
bits of K (bits 1..40 and 55..127), limb 0 of K and the index bits (limb 2)
and fingerprint (limb 3) are pseudo-random functions of the class; excluding the fruitless cycles of F of length <= 64 (2-cycles
at rate 1/(6R) per step, and the rarer short cycles below, which are detected
and escaped and whose repeated arrivals are counted inside T_c), a trail makes
a *non-fruitless return* (an arrival at a class it visited earlier, either
before one of its escapes or along an F-cycle longer than 64) with probability
at most (L_max + P)^2/(2N) per trail, the random-mapping self-collision rate;
the canonical twist index of a new class is uniform on the 6 units and
independent of its partition index; and beta x mod p is equidistributed along
the walk (these last two affect only the cycle and fallback rates, hence only
the meter's early-halt probability).

S0 = sqrt(2 N' ln(1/0.6097)) = 2^126.12851 (R = 2^40); S_ops = 1.005 S0 =
2^126.13570 (H1 gives >= 0.3903; the pure model 0.3939);
pi (S_ops + 2^22) < 2^126.186773 (pi = 1.0360342755859375, section 5).

**Phase 1 covers the budget.** Phase 1 performs at least S_ops batched
operations and each completed operation is one class arrival (an operation
aborted by an exceptional addition is a failure event), so phase 1 makes at
least S_ops arrivals: if T_c <= S_ops, the first inter-trail collision
happens in phase 1. Arrivals in actions (starts, traversals, escapes) only add
to the count. No statement about wasted visits is needed: repeated arrivals
are counted inside T_c.

**Cycle rates and detection delay.** A 2-cycle needs the next class to have
the same partition index and canonical twist -1: probability 1/(6R) =
2^-42.585 per step. Other fruitless cycles of length <= 64 need every
partition index to occur >= 2 times: rate <= sum_L sum_b S2(L, b) R^-(L-b)
= 2^-77.0 (computed in the experiment), and lengths 65..2^14 are smaller
still (each term at most S2(L, L/2) R^-(L/2) < 2^-900). So the expected
number of cycle entries is < 2^-42.58 * 2^126.19 = 2^83.6 (linearity of
expectation over the pseudo-random partition and twist of each new class).
Each costs at most 2P = 2^15 batched operations of detection delay (counted in
T_c and in S_ops, not extra work): together < 2^98.6 operations, a fraction
< 2^-27.5 of phase 1; at R = 2^40 this is negligible, which is why the
checkpoint is tested only every 256th batch.

**Failure events** (union bound with unconditional probabilities; if
T_c <= S_ops, by the merging lemma the two trails reach the same report and
the attack succeeds unless one of these occurs): a colliding trail hits its
cap first (each colliding trail's remaining length to a distinguished class
is geometric with mean 2^73; P(> 20 * 2^73) < e^-20; together with the
Gamma(2) length-biased tail < 1e-7); a non-fruitless return inside one of
the colliding trails (H1: (L_max + P)^2/(2N) = 2^-98.6 per trail, < 2^-97 in total);
a colliding trail entering a fruitless cycle longer than 64 (rate < 2^-900
per step, times 2^77.4 steps); report-log overflow (under H1 the distinguished
bits of distinct classes are independent 2^-73 coins, so new-class reports
are Binomial(n, 2^-73) with n <= 2^126.19 + 2^22 (L_max + P), mean
< 2^53.19; D_cap = 2^54.35 is 2.16 times the mean; Chernoff:
< exp(-2^51.9)); a lookup exceeding 256 probes (< 2^-359 per lookup,
< 2^-300 in total);
f_1 = f_2 (b uniform given the start point: 1/r); exceptional additions
(< 2^-114); 8 processed matches all failing (< 2^-60); a spurious 192-bit
K_cp match between distinct classes (< 2^-190 per test under H1, < 2^118
tests: < 2^-72; spurious 64-bit matches are rejected by the detect action and
cost only its price); meter halting early (< 0.0000034, section 5);
preprocessing cap (< 2^-18). The escape minimum uses full keys, so it is
well defined. **Success >= 0.3903 - 0.0000034 - 2^-18 - 1.1e-7 > 0.39029 >
0.39.**

## 5. Cost

**Scope of the reference operation.** The protected convention prices one
batched incomplete affine G1 addition at one unit and states that the 6M
convention "amortises the batch work and has no standalone per-addition
inversion term". This candidate charges every step one unit for that
operation and, in addition and conservatively, (i) its operand traffic
(row 6: loading the selected table point and the walk's point, storing the
new point) and (ii) the per-batch Montgomery inversion and batch bookkeeping
(last row), as v21/v31 did. Everything the attack does besides the batched
addition is listed instruction by instruction below.

**Per ordinary step** (instructions besides the reference addition; blst Fp
multiplication = 383):

| Row | Work | Instr. |
| --- | --- | ---: |
| 1 | 1.1 est1 = L': 3 `ldp` (limbs of x; the first is pre-indexed, `ldp x0,x1,[x19,#304]!`, and also advances the record pointer), 6 `umulh` (column-5 high halves), 5 `add`, 6 `madd` (column-6 low halves); G0..G6 are resident in x22..x28 | 20 |
| 2 | 1.1 est0 = 9 X + umulh(X, Dfrac) (`umulh`, `madd`; X = limb 5 of x already in x5), est2 = -(est0 + L') (`add`, `neg`) | 4 |
| 3 | 1.1 argmin and gap (3 `cmp`, 6 `csel`, `cset`, `sub`), then the one merged fallback test of tests A-C (`cmn`, `ccmp`, `ccmp`, `b.lo`) | 15 |
| 4 | limbs 0-1 of y (`ldp`), limbs 0-1 of p - y (`subs`, `sbcs`), s = bit 0 of y (`ands`, sets the flags), limbs 0-1 of K (2 `csel`) | 6 |
| 5 | h = bits 1..40 of K (`ubfx`), first-stage distinguished test (`cbz` limb 1 of K), table address table + 256 (8h + 4s + j) (3 `add` with shifted operand), `str` of the entry address | 6 |
| 6 | operand traffic of the batched addition: `ldr` of the entry address, 6 `ldp` table point, 6 `ldp` state point, 6 `stp` new state point | 19 |
| 7 | loop control of the canonicalisation loop (`cmp x19,x20; b.ne loop`) | 2 |
| 8 | per-batch work at its share per step: zero-denominator test on the batch product (<= 9), batch meter charge (<= 10), one Fermat inversion (380 squarings + 228 multiplications = 608 x 383 = 232,864) and its exponent scan (< 2,000), loading the 19 resident registers and pointers (<= 32), batch-index increment and variant dispatch (<= 8): <= 234,922 instructions per batch of 2^22 walks, and 234,922 < 0.0561 x 2^22 = 235,323.8 | < 0.0561 |
| | **total, a per-step maximum, charged without allowance** | **72.0561** |

Rows 1-5 and 7 as instruction sequences (every operand load and store is listed):

```
row 1 (20)  ldp x0,x1,[x19,#304]!; ldp x2,x3,[x19,#16]; ldp x4,x5,[x19,#32]   (x19 += 304 first)
            umulh x9,x5,x22; umulh x8,x4,x23; add x9,x9,x8; madd x9,x5,x23,x9
            umulh x8,x3,x24; add x9,x9,x8; madd x9,x4,x24,x9
            umulh x8,x2,x25; add x9,x9,x8; madd x9,x3,x25,x9
            umulh x8,x1,x26; add x9,x9,x8; madd x9,x2,x26,x9
            umulh x8,x0,x27; add x9,x9,x8; madd x9,x1,x27,x9
            madd x9,x0,x28,x9                                        (x9 = L' = est1)
row 2 (4)   umulh x0,x5,x29; madd x0,x5,x30,x0                        (x0 = est0; x29 = Dfrac, x30 = 9)
            add x2,x0,x9; neg x2,x2                                   (x2 = est2)
row 3 (15)  cmp x9,x0; csel x3,x9,x0,lo; csel x4,x0,x9,lo; cset x6,lo
            cmp x2,x4; csel x7,x2,x4,lo; cmp x2,x3; csel x7,x3,x7,lo; csel x3,x2,x3,lo
            csel x6,x15,x6,lo; sub x7,x7,x3                           (x6 = j, x7 = gap; x15 = 2)
            cmn x9,#14; ccmp x2,#31,#0,cc; ccmp x7,x14,#0,hs; b.lo fallback   (x14 = 41)
row 4 (6)   ldp x0,x1,[x19,#48]                                       (limbs 0-1 of y)
            subs x2,x10,x0; sbcs x3,x11,x1                            (limbs 0-1 of p - y; x10, x11 = p0, p1)
            ands x8,x0,#1                                             (x8 = s = [y odd]; ne iff y odd)
            csel x0,x2,x0,ne; csel x1,x3,x1,ne                        (limbs 0-1 of K)
row 5 (6)   ubfx x4,x0,#1,#40; cbz x1,dpcheck                         (h; limb 1 = 0: stage 2)
            add x5,x6,x8,lsl #2; add x5,x5,x4,lsl #3                  (8h + 4s + j)
            add x5,x21,x5,lsl #8; str x5,[x19,#96]                    (table + 256 (8h + 4s + j))
row 7 (2)   cmp x19,x20; b.ne loop                                    (x20 = last record of the batch)
```

Row 1's 12 products are all present (six `umulh` for the column-5 high halves
x5G0, x4G1, x3G2, x2G3, x1G4, x0G5 and six `madd` for the column-6 low halves
x5G1, x4G2, x3G3, x2G4, x1G5, x0G6); the 12 terms are summed with 11
additions (five `add`, six accumulated by `madd`). Row 2 uses x5 = X, which
row 1 does not overwrite. **The merged fallback test.** `cmn x9,#14` sets C iff
L' + 14 >= 2^64 (test A fails). Only if C is clear (`cc`) does
`ccmp x2,#31,#0,cc` compare est2 with 31 (C set iff est2 >= 31, test B
passes); otherwise it sets NZCV = 0000 (C clear). Only if C is set (`hs`) does
`ccmp x7,x14,#0,hs` compare the gap with 41 (C set iff gap >= 41, test C
passes); otherwise NZCV = 0000. So `b.lo fallback` (C clear) is taken iff test
A, B or C fails, exactly the OR of the three section-1.1 conditions. Part A
executes rows 1-3 instruction by instruction with NZCV flags and compares the
branch with the OR of tests A-C evaluated on exact integers, and j with the
exact rule. **Row 4**: p is odd, so y and p - y have opposite parities, K is
p - y exactly when y is odd, and s = [y odd] is bit 0 of y (`ands`, which also
sets the flags the two `csel` read). A borrow propagates only upward, so
`subs`/`sbcs` give limbs 0-1 of p - y exactly; higher limbs are not needed in
ordinary batches (h uses bits 1..40 of limb 0, the distinguished test limb 1,
the test-batch compare limb 0). Part A emulates rows 4, 4' and 5 limb by limb
against the exact K on random, boundary and distinguished-key y.
**Row 6** keeps the scope of the v31 row (one load of each limb of the table
point and of the walk's point, one store of each limb of the new point, plus
the entry-address load); it is the addition pass's operand traffic, so it is
given as counts of `ldp`/`stp` rather than as a register schedule. The entry
address is handed from the canonicalisation loop to the next batch's addition
pass through the walk record (`str` in row 5, `ldr` in row 6).

**Registers.** For the whole canonicalisation loop of a batch: x10, x11 = limbs
0-1 of p, x14 = 41, x15 = 2, x19 = record pointer (starts one record before
the batch's first record), x20 = address of the last record, x21 = table base,
x22..x28 = G0..G6, x29 = Dfrac, x30 = 9; in refresh batches also x12 = limb 2
of p and x13:x16 = the 128-bit batch index t. These 19 values are loaded once
per batch (row 8, "loading the 19 resident registers and pointers"; the
loop calls nothing, so x29/x30 are free as data registers and are saved and
restored around the loop inside that item). Temporaries are x0-x9 only; x17
and x18 are unused. A taken branch (fallback, dpcheck and report, detect, abandon) leaves
the loop for a separately priced action, which saves and restores what it
clobbers inside its own price and returns into the loop (after a restart or an
escape it re-runs the loop body for the same record, `sub x19,x19,#304`
before re-entry, inside its price).

**Test and refresh variants** (whole batches, chosen once per batch from t):

```
test batch, after row 5 (+3):     ldr x2,[x19,#104]; cmp x0,x2; b.eq detect
refresh batch, row 4' replaces row 4 (+3):
            ldp x0,x1,[x19,#48]; ldr x7,[x19,#64]
            subs x2,x10,x0; sbcs x3,x11,x1; sbcs x9,x12,x7; ands x8,x0,#1
            csel x0,x2,x0,ne; csel x1,x3,x1,ne; csel x7,x9,x7,ne     (limbs 0-2 of K)
            then row 5, the test (+3), and (+6):
            ldp x2,x3,[x19,#128]; cmp x13,x2; sbcs xzr,x16,x3; b.hs abandon   (t >= expiry)
            stp x0,x1,[x19,#104]; str x7,[x19,#120]                  (K_cp <- limbs 0-2 of K)
```

The `b.eq detect` compares limb 0 only; K_cp's limb 0 is 1 when invalid (never
equal to limb 0 of a key, whose bit 0 is 0). The detect action checks limbs
1-2 before doing anything else (section 2).

**Walk record** (304 bytes = the pre-index stride): x at 0, y at 48, next
table-entry address at 96, K_cp limbs 0-1 at 104 and limb 2 at 120, trail
expiry (two limbs) at 128, start coefficients a and b at 144 and 176 (32 bytes
each), batch denominator and prefix product at 208 and 256; every offset a
multiple of 8 and every `ldp`/`stp` offset <= 504. The start batch index is
expiry - L_max.

**Prices** (units): ordinary step 1 + 72.0561/2000 = 1.03602805 (the per-step
maximum, charged without an allowance); a test-batch step 1.03752805 (+3
instructions); a refresh-batch step 1.04202805 (+12). Over phase 1 the
test and refresh batches number at most B1/256 + 1 and B1/2^14 + 1, so the
phase-1 charge is at most pi (S_ops + 2^22) + 2^22 (15/2000) with
**pi = 1.03602805 + (3/256 + 12/16384)/2000 = 1.0360342755859375** (this
counts the refresh batches' 12 instructions on top of their test-batch 3,
i.e. 3 instructions per refresh-batch walk too many). Fallback of section 1.1
0.25 (one Fp multiplication and an exact three-way compare, <= 500
instructions); random 64-bit word (RNDR) 0.005; sampling attempt 0.025;
trail start or restart [a]P0 + [b]T (two Jacobian scalar multiplications and
an inversion, about 4.8e6 instructions; writes a, b, expiry, the invalid K_cp
and re-runs the loop body) 2^12; dpcheck (`tst x0,#0xff80000000000000; b.ne` back into row 5, else `b report`;
entered with probability 2^-64 per step) 0.01; report (completes limbs 2-5
of K, appends the 128-byte log record, writes the index slot) 0.5; index
probe (load a 16-byte slot, compare, advance mod m; a fingerprint hit also
reads and compares the log record's K) 0.015, at most 256 per report; detect action 0.1 fixed (limb-2
recomputation and the limbs 1-2 comparison, the loop re-run) plus 119 per
unbatched operation it performs; **unbatched operation** (one step or
doubling with its own inversion, used by the detect action, phase 2 and
re-walks: 232,864 + 2,000 (inversion and scan) + 5 x 383 (a doubling needs 4M,
an addition 3M) + 192 (field additions) + 98 (canonical index, key, table
address) + 100 (full-key compare and minimum, or coefficient updates mod r)
+ 84 (record and loop) = 237,253 <= 238,000) 119; zero-batch recovery 2^17;
index initialisation (2^57 slots set to empty with one `stp` each plus loop
control, < 1.25 instructions per slot: 2^46.4) (once); table 2^52 (2^40 entries at <= 2^12 units
each, including the six twisted copies, plus their sampling, once); solving
2^9; verification 2^13.

**Budget.** S_ops = 2^126.13570 (section 4). Deterministic part:
pi (S_ops + 2^22) < 2^126.186773.
Hard extras: phase 2 (2^22 (L_max + P) * 119 = 2^106.2167), re-walks
(8 * 2 * ((L_max + P) * 119 + 2^12) = 2^88.22), reports with their maximal
256 probes (D_cap * (0.5 + 256 * 0.015) < 2^56.5), index init 2^46.4, table
2^52, first 2^22 starts 2^34, boundary test/refresh batches 2^14.9, solving
and verification: together < 2^106.22.
Random extras (expected): detect actions (expected detections: cycle-entry
rate (1/(6R) + 2^-77.0) times < 2^126.19 operations of phase 1, phase 2
and re-walks together, plus spurious 64-bit
matches 2^-63 per test, < 2^83.6; each 2-cycle costs 0.1 + 3 * 119 units, and
longer cycles are rarer by factors of R) < 2^92.09; restarts (expected
number of trails < 2^53.2) < 2^65.2; expected index probes (< 3.02 per
report, section 2) < 2^48; dpcheck branches < 2^55.6; fallbacks < 2^69.2;
zero-batch recoveries < 2^-90: together **< 2^92.2**.
Headroom 2^126.1868 - 2^126.186773 - 2^106.22 > 2^110.40, so the meter
halts early with probability < 2^92.2/2^110.40 < 0.0000034 (Markov), and
success is at least 0.3903 - 0.0000034 - 2^-18 - 1.1e-7 > 0.39029. With all
random charges 3 times larger, the bound is < 0.0000100 and success is still
> 0.39028.

**Preprocessing:** jump table (2^52, plus sampling), index initialisation
2^46.4, first 2^22 starts (2^34, plus sampling): the priced part is
< 2^52.03; the cap is 2^56, so the gap exceeds 2^55.9 while the expected
random sampling part (two rejection-sampled scalars per table entry and per
start, < 0.1 units each) is < 2^37: the cap halts a run with probability
< 2^-18 (Markov). All of it is inside W.

**W <= 2^126.1868 on every run.** Memory: report index 2^57 x 16 bytes =
2^61, report log D_cap x 128 bytes = 2^61.35, table 2^43 slots x 256 bytes =
2^51 bytes, walk records 2^22 x 304 bytes < 2^31 bytes, code and constants
< 2^20: < 2^62.3 bytes, every array below 2^63 bytes.

**Encodability.** Every instruction of rows 1-5 and 7, of the test and
refresh variants and of the dpcheck branch is listed as text in
`experiment.py` (`LISTING`) and checked against the AArch64 operand ranges
it uses (`ldp`/`stp` offsets multiples of 8 in [-512, 504], `ldr`/`str`
offsets multiples of 8 up to 32760, add/sub/cmp/cmn 12-bit immediates,
`ccmp`/`ccmn` 5-bit immediates and 4-bit flags, logical bitmask immediates,
`ubfx` fields, shift amounts); negative controls (`ccmp x2,#32,...`,
`ands x8,x0,#0`, an offset of 508, an out-of-range `ubfx`, `cmp #4096`) are
rejected by the checker. The listed ordinary-step count it reports is 72.

## 6. Experiments and limitations

`experiment.py` (stdlib). **Part A** (production BLS12-381; all true in the
declared run): curve, order, beta/lambda pairing, x = 0 points are 3-torsion
(not in G1, so v0 = x >= 1), orbit and unit action, parity-key invariance,
unique canonical member (and the fast rule agrees with it on every orbit
member), the section-1.1 rule and the rows 1-3 listing executed instruction by
instruction with NZCV flags (24,000 random and 8,000 adversarial near-tie
inputs: 0 mismatches of j, 0 mismatches of the merged branch against the OR of
tests A-C, 0 violations of the error bounds d0 in [0, 11.85), d1 in
[0, 12.125), e2 = d0 + d1; of the adversarial inputs, 463, 435 and 683 fail
tests A, B and C), rows 4, 4' and 5 limb by limb (30,000 random, 3,013
boundary and 750 distinguished-key y; key limbs 0-1 and 0-2, parity, h, the
distinguished test, the table address for every j, and the invalid marker never
matching: 0 mismatches), twisted-step class identity, doubling class map,
order of 2 mod r from the full factorisation of r - 1, table entries,
coefficient tracking on production steps, recovery of a planted x from two
production trails forced to merge, and the encodability check of every listed
instruction (section 5). The cost block recomputes every number of
section 5 (ledger, per-batch work, prices, cycle rates by the exact Stirling
recurrence, hard and random extras, headroom, early-halt bound, success) and
checks that 126.1868 is the smallest 4-decimal claim meeting the success
margin.

**Part B** (toy y^2 = x^3 + 4 over F_q, q = 58493047, 861101287 or
10146855109): the complete lockstep pipeline of sections 2-3 at toy scale: M
walks in lockstep batches, test batches every k and refresh batches every 64k
batches (toy k = 1: at toy R the cycle rate is 2^13 times the production rate,
so the production spacing would only add toy-specific repeats), synchronous
traversal and doubling escape, cap in refresh batches, restarts, phase 1 of
ceil(1.005 S0 / M) batches, phase 2, the report log (capacity 2.16 times the
mean) with its linear-probing index at load <= 0.149 and the 256-probe cap
(toy index bits from a multiplicative hash of the short toy key), 8-match cap, replay re-walks from (a, b, start batch) and verification. It
reports P[T_c <= S1] (S1 = 1.005 S0), the H1 quantity, and the pipeline
success rate. The toy pipeline success rate is not comparable with 0.39: at
toy scale trail starts add 2^-dp_bits arrivals per operation and phase 2
completes M open trails of mean length 2^dp_bits, both a large fraction of
S0 (in production, 2^-73 and 2^-31.1 of S_ops). Declared run (q26, R 2048,
dp_bits 2, 16 walks, 14 s): a few hundred trials, so P[T_c <= S1] has a
binomial standard error of about 0.02-0.03 (two local runs: 266 trials, 0.391;
583 trials, 0.377; model 0.393, sweep value at this trail ratio 0.392); every
reported x was correct (a wrong x is rejected by the verification step).

**Lockstep sweep** (offline; `tcsim.c` in this directory is a C
re-implementation of the Part B walk that measures T_c only: same parity key,
exact canonical index, twisted table, partition and distinguished bits, lockstep
schedule with k = 1, synchronous traversal and doubling escape, cap and
restarts; fresh table and target per trial; 16 walks; 1,840,000 trials):

| toy | R | dp_bits | mean trail / sqrt(N) | trials | mean T_c vs model | P[T_c <= S1] |
| --- | ---: | ---: | ---: | ---: | --- | ---: |
| q26 | 2048 | 0 | 0.056% | 200000 | +0.25% +- 0.12% | 0.3914 |
| q26 | 2048 | 1 | 0.111% | 200000 | +0.09% +- 0.12% | 0.3928 |
| q26 | 2048 | 2 | 0.222% | 200000 | +0.26% +- 0.12% | 0.3917 |
| q26 | 2048 | 3 | 0.444% | 200000 | +0.39% +- 0.12% | 0.3911 |
| q26 | 2048 | 4 | 0.888% | 200000 | +0.67% +- 0.12% | 0.3881 |
| q26 | 2048 | 5 | 1.775% | 200000 | +1.18% +- 0.12% | 0.3851 |
| q26 | 2048 | 6 | 3.550% | 200000 | +1.75% +- 0.12% | 0.3808 |
| q26 | 16384 | 3 | 0.444% | 100000 | +0.40% +- 0.17% | 0.3894 |
| q30 | 2048 | 0 | 0.015% | 60000 | +0.01% +- 0.21% | 0.3937 |
| q30 | 2048 | 1 | 0.029% | 60000 | +0.30% +- 0.21% | 0.3911 |
| q30 | 2048 | 3 | 0.116% | 60000 | +0.09% +- 0.21% | 0.3937 |
| q30 | 2048 | 5 | 0.463% | 60000 | +0.62% +- 0.21% | 0.3891 |
| q30 | 2048 | 6 | 0.925% | 60000 | +0.67% +- 0.21% | 0.3871 |
| q34 | 2048 | 2 | 0.017% | 20000 | +0.15% +- 0.37% | 0.3921 |
| q34 | 2048 | 5 | 0.135% | 20000 | -0.36% +- 0.37% | 0.3992 |

Weighted fit slowdown = a + b (trail/sqrt(N)) over all 15 rows:
**a = +0.17% +- 0.05%**, b = 0.47 +- 0.04 (chi^2 = 10.0 on 13 degrees of
freedom); over the 13 rows with trail ratio <= 1%: a = +0.11% +- 0.07%. This
agrees with the v31 sweep of the per-trail walk (a = +0.24% +- 0.18%,
b = 0.86 +- 0.46). The
production ratio 2^-62 sits at the intercept; the 0.5% allowance of H1 is
6.6 standard errors above the fitted a. The toy intercept also contains the
toy's detection-delay repeats (cycle-detection rate about 6e-5 per operation
at R = 2048, each followed by up to about 64 repeated arrivals), which at
R = 2^40 are below 2^-27.5 of all operations; the R = 16384 row
(detections 8e-6 per operation) matches the R = 2048 row at the same trail
ratio. Toy rows with large trail ratios fall below 0.3903, as the trail law
predicts; H1 is a statement at the production ratio.

Limitations: H1 is a heuristic (random mapping with Brent-Pollard correction;
van Oorschot-Wiener parallel collision search; automorphism classes as in
Wiener-Zuccherato, GLV, Duursma-Gaudry-Morain; fruitless-cycle handling as in
Bos-Kleinjung-Lenstra and Bernstein-Lange-Schwabe), supported only at toy
scale; instruction counts are analytic listings, not measurements; the
reference-operation scope follows the protected convention as stated at the
top of section 5; generic attack (no TNFS).
