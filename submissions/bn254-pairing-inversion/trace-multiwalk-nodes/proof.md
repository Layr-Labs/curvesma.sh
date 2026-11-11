# bn254-pairing-inversion: sixteen random-stride trace walks matched across all pairs, on a node machine

**time_log2 = 119.979552**

Exact rational: 3749361/31250 (in lowest terms).

This is a hard cap on every run: every row of section 6 is enforced by a loop bound or a counter. No run
gate is used, and no row is an average.

- *Output:* points P in G1 and Q in G2 with e(P, Q) = z: P = [x]P0 and Q = Q0 for the exponent x found
  in GT (section 2, step 8). The identity z = 1 gives P = infinity (null) and Q = Q0.
- *Success:* at least 0.39 over uniform x and the algorithm's ordinary fresh independent coins: an exact
  coverage inequality minus three loss terms, two of them bounded under the disclosed heuristics H1 and H2
  (section 5).
- *Memory:* at most 2^132.418039 bytes, over nodes that each hold below 2^63 bytes.
- Every transfer between nodes is executed by the receiving node (section 4).
- *Setup:* at most 2^58.0393 units, included in the total.
- Setup means coins, table building, constants and zeroing, done in every run.
- *Not used:* advice, amortisation, a production-size run. A generic square-root method with explicit
  scalar AArch64 accounting; not a sieve attack and not a security estimate.
- *Notation:* p, r, P0, Q0 and g = e(P0, Q0) are the protected modulus, group order and generators;
  C = 20,000; lambda = p mod r; W = 16 walks, numbered w = 0 to 15; walk w uses the exponent w + 1 of z.

## 1. Facts used (mode `algebra`)

- **Orders.** p and r are prime and p has order exactly 12 modulo r, so lambda has order 12 and
  lambda^6 = -1 (mod r). On GT the Frobenius sigma(y) = y^p acts as y -> y^lambda; the experiment checks
  sigma^k(y) = y^(lambda^k) for k = 1, 5, 6, 11 on a random y. An element of GT fixed by sigma^j for some
  0 < j < 12 lies in a proper subfield of Fp12, which meets GT only in 1 (r divides no p^d - 1, d < 12).
- **Trace.** Tr(y) = sum_{k<12} sigma^k(y) lies in Fp and equals 12 Re(e_0) in the pinned tower
  (Fp2 = Fp[u]/(u^2+1), xi = 9+u, Fp6 = Fp2[v]/(v^3-xi), Fp12 = Fp6[w]/(w^2-v)); the experiment checks this by
  summing the 12 conjugates. Tr is constant on each Frobenius class {y^(lambda^k)}.
- **Recurrence.** For y in GT, y != 1, the characteristic polynomial f_y = prod_k (X - sigma^k(y)) lies in
  Fp[X], is monic with constant term 1 and is palindromic (the roots are closed under inversion). So for any
  y0, t_n = Tr(y0 y^n) satisfies t_{n+12} = -(t_n + sum over i = 1 to 5 of c_i (t_{n+i} + t_{n+12-i}) + c_6 t_{n+6}),
  with c_i from the power sums Tr(y^k), k <= 12 (Newton). Checked on the production walks.
- **Denominators.** For walks w != v and every k < 12, and for w = v and every 0 < k < 12, the residue
  (w + 1) - lambda^k (v + 1) is nonzero mod r. For w < v and k != k', the inverses of
  (w + 1) - lambda^k (v + 1) and (w + 1) - lambda^k' (v + 1) differ (they are the coefficients of a coin in
  section 3). Mode `algebra` checks both in exact arithmetic:

Denominators checked: 3,056, all nonzero; coefficient pairs checked: 7,920, all distinct.

## 2. Algorithm

Parameters: N is the least N with mu - mu^2/2 >= 0.3901 for mu = 6 W (W - 1) N^2/r = 12 C(16, 2) N^2/r.
Mode `algebra` checks the inequality in exact rationals for N and its failure for N - 1.

N = 2841460719550568469467457918666427110

Nodes: 16 walkers, 2^36 relays, 2^72 homes; blocks of K = 256 records.

0. **Input and identity.** Every walker holds z (input loading is in the fixed per-node row). Walker 0
   compares z with 1 on all 12 coefficients (inside the output row); if they are equal it outputs
   P = infinity (null) and Q = Q0 and the run ends. Every node zeroes its cursor tables (listing `zero`).
1. **Coins and setup messages.** For each walk w, walker 0 draws delta_w uniform in [1, r) and beta_w
   uniform in [0, r), each by rejection: 4 coin words give a 254-bit value, accepted when it lies in the
   range, at most 128 attempts per draw (failure loss in section 5). It draws salts s1, s2 (2 words). It
   sends every delta_w, beta_w and the salts to each of walkers 1 to 15, one message per walker
   (`setup_send`, handler `h_setup`).
2. **Tables.** Walker w (w >= 1 after waiting for its setup message, `wait`) computes D_w = g^(-delta_w),
   its start y_w = z^(w + 1) g^(beta_w) and c_1 to c_6 of f_(D_w), and fills 30 tables of 2^51 entries: table
   t = 5(i - 1) + d holds -c_i v 2^(51d) mod p at row v (listing `tabbuild`), plus a 32-entry fold table
   v 2^254 mod p.
3. **Walks.** Walker w produces S_(w,n) = Tr(y_w D_w^n), 0 <= n < N; the exponent of y_w D_w^n is
   e_(w,n) = (w + 1) x + beta_w - delta_w n mod r. Each trace step (`tstep`, normal form, machine.md section 2)
   is followed by `emit`: two salted hash words w1, w2; home (d1, d2) with d1 = top 36 bits of w1 and
   d2 = next 28 bits of w1 then top 8 bits of w2; the 40-byte record (t mod 2^128, d2 | w 2^60, n) goes into
   the walker's block for relay d1. A full block is sent (`callsend_w`, `send`).
4. **Relays.** Handler `h_relay` copies each record into its block for home d2 and sends full blocks.
5. **End of input.** After its walk a walker runs `flush_w`: one END message to every relay carrying the
   remaining records and the parameters (z, every delta_w and beta_w). A relay that has received END from
   all 16 walkers sends one END message to each of its homes with its remaining records and the parameters.
6. **Homes.** Handler `h_home` appends each record to the home's one list. A record past the cap sets the
   overflow flag and is not stored.

List cap: 19254486782305944 records per home.

   On END the handler stores the parameters and wakes the home's main processor. If the overflow flag is
   set, the home skips sorting and matching and reports not found (the run then fails unless another home
   finds x). Otherwise it radix-sorts the list (cursor at most the cap) on the 128-bit key in three passes,
   scans it for adjacent equal keys (`scan`), expands each run of equal keys into all its pairs and rechecks
   the pairs in order, at most 2^10 per home, stopping at the first kept x. For records (w, n) and (v, n') it
   computes u = z^(w + 1) g^(beta_w - delta_w n) and y = z^(v + 1) g^(beta_v - delta_v n'), finds the first
   k < 12 with u = sigma^k(y) (12 Frobenius maps), sets
   x = (lambda^k (beta_v - delta_v n') - beta_w + delta_w n) / ((w + 1) - lambda^k (v + 1)) mod r
   and keeps x only if g^x = z. Pairs of records from one walk are rechecked the same way.
7. **Return route.** Each home sends one report (found flag, x) to its relay (`report`); the relay's report
   handler keeps any reported x and counts reports; after its last home it wakes the relay's main processor,
   which reports to walker 0 (`relay_up`). Walker 0 counts the relay reports and keeps x when one was
   found (`output`); otherwise the run fails. There is no restart.
8. **Output points.** Walker 0 computes P = [x]P0 by left-to-right double-and-add from the top bit of x in
   Jacobian coordinates (`outstep.py`: a doubling uses 7 products, a mixed addition of P0 uses 11),
   converts P to affine with one inversion by exponentiation to p - 2, and outputs P and Q = Q0. Before
   each step the accumulator is [k]P0 with 1 <= k < x < r, so a doubling never meets infinity and the added
   P0 is never equal or opposite to it. As attack-side verification it evaluates the reduced pairing
   e(P, Q0) (optimal ate Miller loop over the profile's 6x + 2, then the final exponentiation by binary
   powering of (p^12 - 1)/r) and compares it with z on all 12 coefficients. They agree because the home
   checked g^x = z; a mismatch would end the run as a failure.

**Correctness.** A kept x satisfies g^x = z, checked at the home on all 12 coefficients. By bilinearity
e([x]P0, Q0) = e(P0, Q0)^x = g^x = z; P lies in G1 (cofactor 1) and Q0 is the profile generator of G2, so
no wrong output is possible. z = 1 is answered by step 0 with the infinity output that the rules permit
exactly for z = 1.

**Why a covered x is found.** Suppose e_(w,n) = lambda^k e_(v,n') for walks w < v (a cross triple, section
3). Then y_w D_w^n = sigma^k(y_v D_v^n'), so the two traces are equal; equal traces give equal hash words
and keys, so both records reach the same home and lie in one run of equal keys after the sort, and the pair
is rechecked. If u != 1, k is the only match, since u = sigma^k(y) = sigma^k'(y) would make y fixed by
sigma^(k - k'). If u = y = 1, both exponents are 0 and the equation of step 6 holds for x at every k. The
denominator is nonzero (section 1), so the equation has one solution, which is x. The same holds for a
conjugate pair inside one walk (then k != 0). So every equal-key pair whose elements are conjugates yields x.

## 3. Coverage (exact)

A cross triple is t = (w, n, v, n', k) with w < v, n and n' below N and k < 12; there are 12 C(16, 2) N^2 of
them. Let x_t be the one solution of (w + 1) x + beta_w - delta_w n = lambda^k ((v + 1) x + beta_v - delta_v n'):

  x_t = (lambda^k beta_v - beta_w + delta_w n - lambda^k delta_v n') / ((w + 1) - lambda^k (v + 1)).

x is covered when x = x_t for some cross triple. x_t depends only on the coins and x is uniform and
independent of them, so P(x = x_t) = 1/r. For distinct t, t', x_t - x_t' is an affine function of the coins
with a nonzero coefficient:

- *different walk pairs:* some walk belongs to one pair only, and its beta appears in only one of x_t, x_t',
  with a nonzero coefficient;
- *the same pair, k != k':* the coefficients of beta_w are the negated inverses of the two denominators,
  which differ (section 1);
- *the same pair and k, (n, n') != (m, m'):* the beta terms cancel and the difference is
  (delta_w (n - m) - lambda^k delta_v (n' - m')) / ((w + 1) - lambda^k (v + 1)), where n - m or n' - m' is
  nonzero and below N < r in absolute value.

A nonzero coefficient on a beta (uniform on [0, r)) gives P(x_t = x_t') = 1/r; on a delta (uniform on
[1, r), independent of the other coins) at most 1/(r - 1). So P(x = x_t = x_t') <= 1/(r (r - 1)), and
Bonferroni gives, with T the number of cross triples and mu = T/r < 1,

  P(covered) >= T/r - C(T, 2)/(r (r - 1)) >= mu - mu^2/2 >= 0.3901

(since T <= r, (T - 1)/(r - 1) <= T/r). This is an average over the challenge and the coins, as the claim
states. No bound on overlaps of Frobenius images is needed: each walk has its own delta. Conjugate pairs
inside one walk are rechecked and can only add solved targets; they are not counted.

## 4. Machine model, transport and phase order

- **Nodes.** Each node has a main processor and one delivery processor per inbound link class, all scalar
  AArch64, sharing the node's flat address space; every node holds below 2^63 bytes (section 7).
- **Fixed wiring.** Each link (walker to relay d1, relay to home d2, home to its relay, relay to walker 0,
  walker 0 to walker w) is fixed for the whole run: the block address the sender passes to `send` selects
  the receiver, and no executed lookup or routing table is involved beyond the address arithmetic counted
  in the listings.
- **Receive window.** A delivery processor has one read-only receive window at a fixed address, which the
  record handlers load with MOVZ on entry. While it delivers a message,
  the window maps exactly that message's block in the sender's memory (64-byte header, K records of 40 bytes,
  a 1,408-byte parameter area; a setup message is 64 + 8 x 130 bytes); reads stay inside it and stores to
  it are forbidden (the simulation asserts both). No link holds a copy: the receiver's executed loads are
  the transfer.
- **Completion.** `send` writes the count and END words and then the sequence word with STLR. That store
  retires only after the receiver's handler returns, so a block is never rewritten in flight. The handler's
  first window load is LDAR of the sequence word. The sequence word and the sent and handled counters
  (STADD) are only incremented and stored; no listing compares their values or uses them as addresses, so
  their wrap modulo 2^64 at production size has no effect.
- **Arbitration.** Senders to the same delivery processor are served one at a time, lowest link index first;
  a waiting sender stalls at its completing store. Handlers commute: list contents are multisets, END and
  report counters are sums, and a kept x is always verified. The simulation runs two schedules and requires
  identical home contents and the same output.
- **No deadlock.** Walker main processors wait only on relay delivery processor 1, which waits only on home
  delivery processors. Home handlers, report handlers and `h_setup` never send. Main processors of homes
  and relays wait on their own wake word, then send upstream to delivery processors whose handlers never
  wait; a walker's setup wait depends only on walker 0, which sends its setup messages before it walks.
  The waits-for graph has no cycle.
- **Waits.** `wait` is LDAXR, CBNZ, WFE, B on the node's wake word at state + 4096, the only word written
  in its aligned 2^11-byte block (2^11 bytes is the architectural maximum of the exclusives reservation
  granule; CTR_EL0.ERG gives log2 of the granule size in words, at most 9: 512 words), with the event
  stream off and interrupts masked: only stores to that word wake it, plus at most one stale event. Each
  wake word is written once (the simulation records every store into each
  node's wake block), so a wait runs at most 3 passes: 11 instructions, charged on every wait.
- **Calls.** BL writes x30 and every RET must return to it; the interpreter asserts this and that each
  routine returns with x30 and sp unchanged. A delivery processor's registers persist between messages:
  the simulation sets its state register once, checks it before every message, and overwrites the window
  register with a junk value before each record handler. Each delivery processor's x30 holds its handler's
  entry address, set once; every handler leaves x30 unchanged (the record handlers save and restore it),
  so its final RET branches back to its own entry and waits for the next message. No instruction runs between messages besides the handler's own.
  The calls of the main routines between listings (walker: tables, `walk_all`,
  flush, wait, output; relay: wait, `relay_up`) are in the fixed per-node row; a home's are in its driver row.
- **Phase order:** identity test; coins and setup messages; tables and walks; END from all walkers; relay
  flushes; home sort, scan and rechecks; home reports; relay reports; output points and pairing check.
  Every step after the coins is triggered by a message or a wake word, so no barrier is needed.

Every send, handler, flush, report and wait is a listing executed by the experiment (counts in
machine.md); the simulation (mode `transport`) checks that the executed totals of every processor equal the
cost functions used below.

## 5. Success: losses and heuristics

| | Cause | Bound |
| --- | --- | ---: |
| L1 | one of the 32 coin draws fails 128 times: 32 (1 - (r - 1)/2^254)^128 | < 2^-255.5 |
| L2 | a list exceeds its cap (overflow) | < 2^-48.7 |
| L3 | the home of a covering pair uses its 2^10 rechecks on foreign pairs first | < 2^-30.8 |

- **L1.** (r - 1)/2^254 is the smaller of the two acceptance probabilities (delta: r - 1 values, beta: r).
- **L2.** The cap is twice the ceiling of the mean list size M/2^72, M = 16 N records. Repeated trace
  values come from equal-trace pairs of records: an H1 coincidence (at most 2^-250 per pair) or a pair whose
  elements are conjugates, which for uniform x has probability at most 12/r per pair (one event of
  probability 1/r per k). Mode `ledger` checks in exact rationals that the expected number,
  C(M, 2)(2^-250 + 12/r), is within this line:

Expected repeats in all lists: at most 13M^2/r < 2.

  So at least 2^50 repeats has probability below 2^-49 (Markov). Under H2 the distinct values in a list are
  binomial with mean M/2^72, and the Chernoff bound for reaching the cap minus 2^50 (still almost twice
  that mean), times 2^72 lists, is far smaller.
- **L3.** A foreign pair has equal 128-bit keys but elements that are not conjugates; every other equal-key
  pair yields x (section 2). Under H1 and H2 the expected number of foreign pairs at a home is at most
  C(cap, 2)(2^-128 + 2^-250):

Expected foreign pairs at a home: below 2^-20.8.

  A covering pair is reached unless at least 2^10 foreign pairs come first (Markov).
- **Result.** Mode `ledger` computes it with exact rationals and high-precision decimals before rounding
  each exponent down:

Success >= 0.3901 - 2^-255.5 - 2^-48.7 - 2^-30.8 >= 0.3900999994 > 0.39.

- **H1 (trace coincidences).** For elements in different Frobenius classes, the 128 key bits of the trace
  behave like uniform values; equal full traces of different classes are no more likely than 2^-250 per
  pair. Used only through Markov bounds, so an error of a large factor leaves the result intact.
- **H2 (salted hash).** For distinct traces, the home bits of the salted hash under fresh salts behave like
  independent uniform values. Used for L2 and L3.

**Toy evidence** (mode `toy`, toy scale only): the whole algorithm (16 walks with fresh coins, routing over
relays and homes, caps, sort, run expansion, 12-conjugate recheck) in the BN-family fields q = p(-15)
(small) and q = p(-41) (large), with the toy's own least N:

Walk length in the toy: N = 26 (small), N = 192 (large).

For every target without overflow, solved equals the exact ground truth (some two records whose exponents
lie in one Frobenius class, computed from x and the coins), and no output is wrong. As a negative control
the same run also keeps only x found from pairs inside one walk; that count must equal the targets whose
class pair lies inside one walk, so the cross-walk pairs carry the coverage. The schema's `toy_cap` plants a
small cap: those runs overflow and fail. The toy counts a target with any overflowing list as failed, while
production skips only the overflowing home, so the toy is the stricter test.

- *Run fraction.* The cross-walk covered targets of one run, over toy_targets = n, are binomial; the spread
  is at most 0.5/sqrt(n) for each allowed n:

| Targets | Spread at most |
| ---: | ---: |
| 20 | 0.112 |
| 40 | 0.080 |
| 60 | 0.065 |

  So one run can print a fraction below 0.39 although its mean is above. The claim does not rest on this
  fraction; it is printed with its spread.
- *Exact coverage.* Every toy run lists all cross triples of the small field for fresh coin draws and
  counts the covered x exactly (no sampling); each draw's fraction must reach the Bonferroni bound of
  section 3 for the small field's N (printed rounded down to 4 decimals):

Exact coverage, small field: r' = 1,704,961, N = 26, 973,440 cross triples, 2 coin draws, Bonferroni bound 0.4079.

## 6. Work ledger (instructions; every row a hard cap)

Records: 16 N = 45463371512809095511479326698662833760.

Full blocks at hop 1: at most 177591294971910529341716119916651680.

Full blocks at hop 2: at most 177591294971910529341716119916651694.

These are 16 floor(N/K) and the floor of the records over K; each full block holds K records.

Per record: **576** executed instructions.

The rows are the cost functions of the experiment evaluated at production size; mode `transport` checks the
same functions against executed totals (its walkers call `walk` directly; `walk_all` is executed in mode
`kernels`).

| Row | Every-run bound | log2 units, rounded up |
| --- | ---: | ---: |
| walkers | 20,920,668,927,379,328,147,768,156,202,587,273,808,846 | 119.6544 |
| relay handlers | 873,038,806,081,913,130,329,005,678,016,548,856,518 | 115.0717 |
| home handlers | 687,811,085,426,210,519,061,092,763,759,139,025,422 | 114.7276 |
| radix sort | 3,365,192,951,601,258,184,525,232,276,498,076,559,296 | 117.0182 |
| scan | 363,706,972,102,511,525,275,926,007,637,216,686,848 | 113.8084 |
| rechecks | 104,813,394,661,667,060,512,388,974,379,008 | 92.0819 |
| reports | 273,897,256,011,112,346,814,548 | 63.5703 |
| home drivers | 19,342,813,113,834,066,795,298,816 | 69.7123 |
| setup | 5,923,494,519,332,151,664,240 | 58.0393 |
| node fixed | 4,951,760,157,213,578,693,651,202,048 | 77.7123 |
| output | 65,536 | 1.7123 |
| point output | 4,194,304 | 7.7123 |
| pairing check | 268,435,456 | 13.7123 |
| **total** | **26,210,418,847,409,587,551,417,551,384,092,446,990,886** | **119.979552** |

- *walkers:* 16 `walk_init`; per walker `walk_all` (7 + 4 q) and its q + 1 `walk` calls (3 + 460 n +
  floor((n - 1)/12) for n steps; `walk` counts steps in one 64-bit register, so q = floor((N - 13)/c)
  calls take c = 3 x 2^60 steps and the last takes the other N - 12 - q c), 21 per hop-1 full block
  (`callsend_w` with `send`, which replaces the RET of `emit`), 16 flushes (6 plus 207 per relay over
  2^36 relays), `setup_send` (15 setup messages), 16 capped waits (15 setup waits, walker 0's report
  wait) and `output`.
- *relay handlers:* (hop-1 full blocks + 16 END messages per relay) x fixed part + 19 per record + the
  full-block extra per hop-2 full block + the END extra per END message + per relay one flush, 2 plus 205
  per home over 2^36 homes.
- *home handlers:* (hop-2 full blocks + one END per home) x fixed part + 15 per record + END extra per home.
- *radix sort:* per record the three passes, (9 + 15) + (10 + 16) + (9 + 15); per home and pass, `zero` over
  2^43 counters (5 per 4 counters, 1 return), `prefix` over 2^43 counters (6 per counter, 1 return) and the
  returns of `count` and `scatter` (2).
- *scan:* 8 per record plus 8 per call, at most 2^10 + 2 calls per home.
- *rechecks:* per home 2^10 pairs x (1,554 Fp12 products + 2^10 products mod r at the executed `fr_mul`
  count + 2^12 + 64): u and y (each a power of g with at most 508 products, a power z^(w + 1) with at
  most 8 and one product), 12 Frobenius maps (each at most one product), g^x (508); the products mod r
  cover the two exponents of the records, the inverse of the denominator by exponentiation to r - 2 (at
  most 2 x 254 products), x (lambda^k is one of 12 public constants lambda^k mod r stored with the code),
  and the conversions into and out of Montgomery form; 2^12 for sums, compares and the loop over k, 64 for
  the pair routine.
- *reports:* every `report`, report handler (25 each, 27 for the last at each receiver), `relay_up` and
  capped wait on homes and relays, walker 0's report handlers, and the 15 `h_setup` messages.
- *home drivers:* 2^12 per home: the overflow test, the calls of `wait`, the radix passes, `scan` and
  `report`, register set-up and buffer swaps. The call, test and branch back around a `scan` call that
  returns a pair are inside that pair's 64; those of the other calls are inside the 2^12.
- *setup:* coins 524,320 (32 draws x 128 attempts x (4 words x 16 + 64) + 2 salt words x 16), per walker
  4,096 Fp12 products + 2^24 (D_w, y_w, power sums, Newton, ring start), each walker's 30 tables at 19 per
  entry for 2^51 - 1 entries plus 1 return and 2^10 per table, 2^10 per fold table, and cursor zeroing of
  16 walkers and 2^36 relays (5 per 4 words over 2^36 words, 1 return).
- *node fixed:* 2^20 per node (code loading, stack, state block, input, the main routine's calls, and each
  delivery processor's re-entry up to its blocking LDAR after its last message).
- *output:* 2^16, including the identity test of step 0.
- *point output:* 2^22 for step 8's scalar multiplication: at most 2 x 253 Jacobian steps for x < 2^254,
  each at most 11 products and 14 sums, every field operation charged at the executed `fp_mul` count, plus
  128 per step for the bit test, branch and call; the inversion (at most 2 x 254 products) and 4 products
  of the affine conversion; 4,096 for encoding and control:

Scalar multiplication: 2 x 253 x ((11 + 14) x 155 + 128) + (2 x 254 + 4) x 155 + 4,096 = 2,108,974 <= 2^22.

- *pairing check:* 2^28 for the attack-side check of step 8. Fp12 products, each charged at the executed
  Fp12 product: the Miller loop over the 65 bits of 6x + 2 (per bit a squaring and a line product, per
  set bit one more line product, two final line products for the Frobenius images of Q0) and the final
  exponentiation by binary powering of (p^12 - 1)/r, 2,790 bits, with 128 per exponent bit for the bit
  test, branch and call. G2 steps (a doubling per bit, an addition per set bit, two final additions), each
  in affine coordinates with its line coefficients: at most one Fp2 inversion (two squarings, one sum, an
  Fp inversion by exponentiation, two products, one negation), 12 Fp2 products (4 products and 2 sums
  each, enough for the Frobenius images of Q0 too) and 40 Fp2 sums (2 sums each), every field operation
  charged at `fp_mul`, plus 4,096 per step for embedding the line and control; 4,096 for the comparison
  with z:

One G2 step: 2 x 254 + 6 + 12 x 6 + 40 x 2 = 666 field operations.

Pairing check: (3 x 65 + 2 + 2 x 2,790) x 13,843 + 2,790 x 128 + (2 x 65 + 2) x (666 x 155 + 4,096) + 4,096 = 94,499,259 <= 2^28.

W = total/20,000. The claim **119.979552** is log2 W rounded up at the sixth decimal; mode `ledger` prints
the unrounded value and the exact integers.

## 7. Memory per node (every mode checks each below 2^63 and the address ranges disjoint)

| Role | Count | Bytes |
| --- | ---: | ---: |
| walker | 16 | 2,162,533,213,406,266,368 |
| relay | 2^36 | 805,392,268,394,496 |
| home | 2^72 | 1,540,429,311,329,701,760 |

Peak total: 2^132.418039 bytes (rounded up).

Layout in machine.md section 4.

## 8. Checks and limitations

- `experiment.py` (with `progs.py`, `outstep.py` and `toy.py`) runs every check in every mode, except the
  toy, which runs in modes `toy` and `all`; it exits with status 1 when any check fails and names the failed
  checks.
- Every mode compares the build-filled lines of this file and machine.md (the figure lines, the tables, the
  cross-triple counts, the output-step operation counts and bounds, and the walkers and rechecks items of
  section 6) with the computed values (mode `ledger` prints the comparison). Any other numeral in the two files must equal a computed value, a code
  constant or a listed structural constant; that check is by value, so it does not tie such a numeral to
  its line. The claim may print no other figures.
- It runs a mutation test: a table step one instruction longer or of 1000 instructions, a relay record, home
  record, scan step or second-pass scatter one instruction longer, a swapped stated bound, a swapped memory
  figure, a swapped program count and a swapped claim figure must each make the comparison fail.
- It loads the shipped copy of the protected target profile, checks its sha256 against the protected value

649e2c45757c3fd789e168fe1b7ea58a40b924640476551586a177ebea3c1ed0

  and asserts p, r, x, g, P0, Q0, the ate loop count and the Fp12 basis. It probes every immediate, label and condition form the
  generator emits with out-of-range values and x18:

43 probes rejected, 20 in-range forms accepted.

- The success bound uses H1 and H2. Counts are executed counts of generated programs on the interpreter, not
  timings. The node machine (windows, stalling completion store, wake-only-on-store waits) is a modelling
  premise. The simulation sorts with Python (the radix listings run in mode `kernels`), and expands
  equal-key runs and computes the recheck with Python big integers (both charged in the rechecks row).
  Nothing runs at production size; sieve methods for the degree-12 field are expected to be cheaper.
