# bn254-gt-dlp: trace baby-step giant-step with a random stride on a node machine

**time_log2 = 120.534383**

Exact rational: 120534383/1000000.

This is a hard cap on every run: every row of section 6 is enforced by a loop bound or a counter. No run
gate is used, and no row is an average.

- *Success:* at least 0.39 over uniform x and the algorithm's ordinary fresh independent coins: an exact
  coverage inequality minus three loss terms, two of them bounded under the disclosed heuristics H1 and H2
  (section 5).
- *Memory:* at most 2^132.871418 bytes, over nodes that each hold below 2^63 bytes.
- Every transfer between nodes is executed by the receiving node (section 4).
- *Setup:* at most 2^58.0349 units, included in the total.
- Setup means coins, table building, constants and zeroing, done in every run.
- *Not used:* advice, amortisation, a production-size run. A generic square-root method with explicit
  scalar AArch64 accounting; not a sieve attack and not a security estimate.
- *Notation:* p, r and g are the protected modulus, group order and generator; C = 20,000; lambda = p mod r.

## 1. Facts used (mode `algebra`)

- **Orders.** p and r are prime and p has order exactly 12 modulo r, so lambda^6 = -1 (mod r). On GT the
  Frobenius sigma(y) = y^p acts as y -> y^lambda; the experiment checks sigma^k(y) = y^(lambda^k) for
  k = 1, 5, 6, 11 on a random y.
- **Trace.** Tr(y) = sum_{k<12} sigma^k(y) lies in Fp and equals 12 Re(e_0) in the pinned tower
  (Fp2 = Fp[u]/(u^2+1), xi = 9+u, Fp6 = Fp2[v]/(v^3-xi), Fp12 = Fp6[w]/(w^2-v)); the experiment checks this by
  summing the 12 conjugates. Tr is constant on each Frobenius class {y^(lambda^k)}.
- **Recurrence.** For y in GT, y != 1, the characteristic polynomial f_y = prod_k (X - sigma^k(y)) lies in
  Fp[X], is monic with constant term 1 and is palindromic (the roots are closed under inversion). So for any
  y0, t_n = Tr(y0 y^n) satisfies t_{n+12} = -(t_n + sum over i = 1 to 5 of c_i (t_{n+i} + t_{n+12-i}) + c_6 t_{n+6}),
  with c_i from the power sums Tr(y^k), k <= 12 (Newton). Checked on the production walks.
- **Disjoint images.** For k = 1 to 11, k != 6, the lattice {(a, b) : b = lambda^k a mod r} has no nonzero
  vector of sup-norm at most m. Two-dimensional reduction bounds the sup-norm from below by half the
  shortest Euclidean length:

Sup-norm lower bound: 2^125.798; log2 m <= 124.55.

  Since lambda^6 = -1 and 2m < r, the 12 sets lambda^k [1, m] are pairwise disjoint and miss 0. So the
  baby set B = {j lambda^k mod r : 0 <= j <= m, k < 12} has exactly 12m + 1 elements.

## 2. Algorithm

Parameters: m = h, the least m with mu - mu^2/2 >= 0.3901 for mu = h(12m + 1)/r. Mode `algebra` checks the
inequality in exact rationals for m and its failure for m - 1.

m = h = 31126642647254138422725958960262852045

Nodes: 2 walkers, 2^36 relays, 2^72 homes; blocks of K = 256 records.

0. **Input and identity.** Both walkers hold z (input loading is in the fixed per-node row). Walker 0
   compares z with 1 on all 12 coefficients (inside the output row); if they are equal it outputs 0 and the
   run ends. Every node zeroes its cursor tables (listing `zero`).
1. **Coins and setup message.** Walker 0 draws delta uniform in [1, r) by rejection: 4 coin words give a
   254-bit value, accepted when 1 <= value < r, at most 128 attempts (failure loss in section 5). It draws
   salts s1, s2 (2 words). It sends (delta, s1, s2) to walker 1 (`setup_send`, handler `h_setup`).
2. **Tables.** Walker 0 uses y = g; walker 1 waits for the setup message (`wait`) and uses y = D = g^-delta.
   Each computes c_1 to c_6 of f_y and fills 30 tables of 2^51 entries: table t = 5(i - 1) + d holds
   -c_i v 2^(51d) mod p at row v (listing `tabbuild`), plus a 32-entry fold table v 2^254 mod p.
3. **Walks.** Walker 0 produces T_j = Tr(g^j), 0 <= j <= m; walker 1 produces S_i = Tr(z D^i) = Tr(z g^(-delta i)),
   0 <= i < h. Each trace step (`tstep`, normal form, machine.md section 2) is followed by `emit`: two
   salted hash words w1, w2; home (d1, d2) with d1 = top 36 bits of w1 and d2 = next 28 bits of w1 then top
   8 bits of w2; the 40-byte record (t mod 2^128, d2 | list bit, index) goes into the walker's block for
   relay d1. A full block is sent (`callsend_w`, `send`).
4. **Relays.** Handler `h_relay` copies each record into its block for home d2 and sends full blocks.
5. **End of input.** After its walk a walker runs `flush_w`: one END message to every relay carrying the
   remaining records and the parameters (z, delta). A relay that has received END from both walkers sends
   one END message to each of its homes with its remaining records and the parameters.
6. **Homes.** Handler `h_home` appends each record to the baby or giant list. A record past the cap sets
   the overflow flag and is not stored.

List cap: 13182645929817536 records per list.

   On END the handler stores the parameters and wakes the home's main processor. If the overflow flag is
   set, the home skips sorting and matching and reports not found (the run then fails unless another home
   finds x). Otherwise it radix-sorts both lists (cursors at most the cap) on the 128-bit key, merges them
   and rechecks every equal-key pair (j, i), at most 2^10 pairs per home: it computes u = g^j and y = z D^i,
   finds k < 12 with y = sigma^k(u) (12 Frobenius maps), sets x = (delta i + j lambda^k) mod r and keeps x
   only if g^x = z.
7. **Return route.** Each home sends one report (found flag, x) to its relay (`report`); the relay's report
   handler keeps any reported x and counts reports; after its last home it wakes the relay's main processor,
   which reports to walker 0 (`relay_up`). Walker 0 counts the relay reports and outputs x when one was
   found (`output`); otherwise the run fails. There is no restart.

**Correctness.** An output x != 0 satisfies g^x = z, checked at the home on all 12 coefficients, so no
wrong output is possible; z = 1 is answered by step 0.

**Why a covered x is found.** Suppose x - delta i is in B for some i < h, say x - delta i = j lambda^k. Then
z D^i = g^(j lambda^k) = sigma^k(g^j), whose trace equals T_j. Equal traces give equal hash words and keys,
so both records reach the same home; the merge stops at their run and the recheck finds k.

## 3. Coverage (exact)

Let A_i = delta i + B (mod r). x is covered when it lies in the union of A_0 to A_{h-1}. Average over uniform
x and delta uniform in [1, r): P(x in A_i) = (12m + 1)/r. For i != i' and fixed b != b' in B,
P(delta (i' - i) = b - b') = 1/(r - 1), so P(x in A_i and A_i') <= (12m + 1)^2/r^2. Bonferroni gives

  P(covered) >= mu - mu^2/2 >= 0.3901, mu = h(12m + 1)/r.

This is an average over the challenge and the coin delta, as the claim states; it needs |B| = 12m + 1
(section 1) and delta uniform on the nonzero residues (rejection sampling gives exactly that).

## 4. Machine model, transport and phase order

- **Nodes.** Each node has a main processor and one delivery processor per inbound link class, all scalar
  AArch64, sharing the node's flat address space; every node holds below 2^63 bytes (section 7).
- **Fixed wiring.** Each link (walker to relay d1, relay to home d2, home to its relay, relay to walker 0,
  walker 0 to walker 1) is fixed for the whole run: the block address the sender passes to `send` selects
  the receiver, and no executed lookup or routing table is involved beyond the address arithmetic counted
  in the listings.
- **Receive window.** A delivery processor has one read-only receive window at a fixed address, which the
  record handlers load with MOVZ on entry. While it delivers a message,
  the window maps exactly that message's block in the sender's memory (64-byte header, K records of 40 bytes,
  a 416-byte parameter area); reads stay inside it and stores to it are forbidden (the simulation asserts
  both). No link holds a copy: the receiver's executed loads are the transfer.
- **Completion.** `send` writes the count and END words and then the sequence word with STLR. That store
  retires only after the receiver's handler returns, so a block is never rewritten in flight. The handler's
  first window load is LDAR of the sequence word.
- **Arbitration.** Senders to the same delivery processor are served one at a time, lowest link index first;
  a waiting sender stalls at its completing store. Handlers commute: list contents are multisets, END and
  report counters are sums, and a kept x is always verified. The simulation runs two schedules and requires
  identical home contents and the same output.
- **No deadlock.** Walker main processors wait only on relay delivery processor 1, which waits only on home
  delivery processors. Home handlers, report handlers and `h_setup` never send. Main processors of homes
  and relays wait on their own wake word, then send upstream to delivery processors whose handlers never
  wait. The waits-for graph has no cycle.
- **Waits.** `wait` is LDAXR, CBNZ, WFE, B on the node's wake word at state + 4096, the only word written
  in its aligned 2^11-byte block (2^11 bytes is the architectural maximum of the exclusives reservation
  granule), with the event stream off and interrupts masked: only stores to that word wake it, plus at
  most one stale event. Each wake word is written once (the simulation records every store into each
  node's wake block), so a wait runs at most 3 passes: 11 instructions, charged on every wait.
- **Calls.** BL writes x30 and every RET must return to it; the interpreter asserts this and that each
  routine returns with x30 and sp unchanged. A delivery processor's registers persist between messages:
  the simulation sets its state register once, checks it before every message, and overwrites the window
  register with a junk value before each record handler. Each delivery processor's x30 holds its handler's
  entry address, set once; the handler saves and restores x30, so its final RET branches back to its own
  entry and waits for the next message. No instruction runs between messages besides the handler's own.
  The calls of the main routines between listings (walker: tables, walks,
  flush, wait, output; relay: wait, `relay_up`) are in the fixed per-node row; a home's are in its driver row.
- **Phase order:** identity test; coins and setup message; tables and walks; END from both walkers; relay
  flushes; home sort, merge and rechecks; home reports; relay reports; output. Every step after the coins
  is triggered by a message or a wake word, so no barrier is needed.

Every send, handler, flush, report and wait is a listing executed by the experiment (counts in
machine.md); the simulation (mode `transport`) checks that the executed totals of every processor equal the
cost functions used below.

## 5. Success: losses and heuristics

| | Cause | Bound |
| --- | --- | ---: |
| L1 | delta sampling fails 128 times: (1 - (r - 1)/2^254)^128 | < 2^-260.5 |
| L2 | a list exceeds its cap (overflow) | < 2^-49.3 |
| L3 | the home of a covering pair uses its 2^10 rechecks on foreign pairs first | < 2^-30.9 |

- **L2.** The cap is twice the ceiling of the larger mean list size. Repeated trace values come from
  equal-trace pairs in one list: a baby pair only by an H1 coincidence (at most 2^-250 per pair); a giant
  pair by an H1 coincidence or when z D^i' is a conjugate sigma^k(z D^i), k != 0, which for uniform x has
  probability 1/r per (i, i', k), so 11/r per pair. Mode `ledger` checks in exact rationals that the
  expected number of such pairs, C(m + 1, 2) 2^-250 + C(h, 2)(2^-250 + 11/r), is within this line:

Expected repeats in all lists: at most 18(m + 1)^2/r < 1.

  So at least 2^50 repeats has probability below 2^-50 (Markov). Under H2 the distinct values in a list are
  binomial with mean at most (m + 1)/2^72, and the Chernoff bound for reaching the cap minus 2^50 (still
  almost twice that mean), times 2^73 lists, is far smaller.
- **L3.** A foreign pair has equal 128-bit keys but traces that are not Frobenius conjugates. Under H1 and
  H2 its expected number at a home is at most cap^2 (2^-128 + 2^-250):

Expected foreign pairs at a home: below 2^-20.9.

  A covering pair is reached unless at least 2^10 foreign pairs come first (Markov).
- **Result.** Mode `ledger` computes it in exact arithmetic before rounding each exponent down:

Success >= 0.3901 - 2^-260.5 - 2^-49.3 - 2^-30.9 >= 0.3900999995 > 0.39.

- **H1 (trace coincidences).** For elements in different Frobenius classes, the 128 key bits of the trace
  behave like uniform values; equal full traces of different classes are no more likely than 2^-250 per
  pair. Used only through Markov bounds, so an error of a large factor leaves the result intact.
- **H2 (salted hash).** For distinct traces, the home bits of the salted hash under fresh salts behave like
  independent uniform values. Used for L2 and L3.

**Toy evidence** (mode `toy`, toy scale only): the whole algorithm (random delta, routing over relays and
homes, caps, sort, merge, 12-conjugate recheck) in the BN-family fields q = p(-15) (small) and q = p(-41)
(large). For every target without overflow, solved equals exact membership of x in the union of the A_i,
and no output is wrong. The schema's `toy_cap` plants a small cap as a negative control: those runs
overflow and fail. The toy counts a target with any overflowing list as failed, while production skips
only the overflowing home, so the toy is the stricter test.

- *Run fraction.* The covered targets of one run, over toy_targets = n, are binomial; the spread is at most
  0.5/sqrt(n) for each allowed n:

| Targets | Spread at most |
| ---: | ---: |
| 20 | 0.112 |
| 40 | 0.080 |
| 60 | 0.065 |

  So one run can print a fraction below 0.39 although its mean is above. The claim does not rest on this
  fraction; it is printed with its spread.
- *Coverage estimate.* Every toy run also tests fresh (x, delta) pairs by exact membership in the union of
  the A_i, with no group operations, and requires the fraction to be at least the Bonferroni bound of
  section 3 for the toy's own m, h and |B|. It prints the margin in spreads:

Coverage estimate, small field: 100,000 pairs, spread at most 0.0016, Bonferroni bound 0.3906.

Coverage estimate, large field: 30,000 pairs, spread at most 0.0029, Bonferroni bound 0.3902.

## 6. Work ledger (instructions; every row a hard cap)

Records: (m + 1) + h = 62253285294508276845451917920525704091.

Full blocks at hop 1: at most 243176895681672956427546554377053530.

Full blocks at hop 2: at most 243176895681672956427546554377053531.

These are floor((m + 1)/K) + floor(h/K) and the floor of the records over K; each full block holds K
records.

Per record: **618** executed instructions.

The rows are the cost functions of the experiment evaluated at production size; mode `transport` checks the
same functions against executed totals.

| Row | Every-run bound | log2 units, rounded up |
| --- | ---: | ---: |
| walkers | 28,646,805,724,057,664,837,396,648,391,493,251,889,629 | 120.1078 |
| relay handlers | 1,195,457,619,171,104,607,975,305,089,322,808,861,406 | 115.5251 |
| home handlers | 1,252,604,189,656,297,814,126,542,794,124,981,546,529 | 115.5925 |
| radix sort | 6,225,329,705,830,584,847,457,803,832,821,791,836,300 | 117.9057 |
| merge | 1,182,812,420,595,710,556,691,712,107,305,870,150,785 | 115.5098 |
| rechecks | 102,907,880,949,087,203,548,700,400,943,104 | 92.0554 |
| reports | 273,897,256,011,112,346,812,430 | 63.5703 |
| home drivers | 19,342,813,113,834,066,795,298,816 | 69.7123 |
| setup | 5,905,525,155,615,322,903,530 | 58.0349 |
| node fixed | 4,951,760,157,213,578,693,636,521,984 | 77.7123 |
| output | 65,536 | 1.7123 |
| **total** | **38,503,009,762,224,214,995,508,324,343,257,206,830,049** | **120.534383** |

- *walkers:* two `walk_init`, two `walk` calls (3 + 460 n + floor((n - 1)/12) for n steps; n is m + 1 - 12
  and h - 12), 21 per hop-1 full block (`callsend_w` with `send`, which replaces the RET of `emit`), two
  flushes (6 plus 77 per relay over 2^36 relays), `setup_send`, two capped waits and `output`.
- *relay handlers:* (hop-1 full blocks + 2 END messages per relay) x fixed part + 19 per record + the
  full-block extra per hop-2 full block + the END extra per END message + per relay one flush, 2 plus 75 per
  home over 2^36 homes.
- *home handlers:* (hop-2 full blocks + one END per home) x fixed part + 20 per record + END extra per home.
- *radix sort:* 4 passes x (9 + 16) per record; per home, list and pass, `zero` over 2^32 counters (5 per 4
  counters, 1 return), `prefix` over 2^32 counters (6 per counter, 1 return) and the returns of `count` and
  `scatter` (2).
- *merge:* 19 per record plus 11 per call, at most 2^10 + 2 calls per home.
- *rechecks:* per home 2^10 pairs x (1,537 Fp12 products + 2^12 + 64): g^j, D^i (at most 508 products
  each), the product with z, 12 Frobenius maps (each at most one product), g^x (508); 2^12 for scalar
  arithmetic and compares, 64 for the pair routine.
- *reports:* every `report`, report handler (25 each, 27 for the last at each receiver), `relay_up` and
  capped wait on homes and relays, walker 0's report handler, and `h_setup`.
- *home drivers:* 2^12 per home: the overflow test, the calls of `wait`, the 8 radix passes, `merge` and
  `report`, register set-up and buffer swaps. The call, test and branch back around a `merge` call that
  returns a pair are inside that pair's 64; those of the other calls are inside the 2^12.
- *setup:* coins 16,416 (128 attempts x (4 words x 16 + 64) + 2 salt words x 16), per walker 4,096 Fp12
  products + 2^24 (D, power sums, Newton, ring start), the 60 tables at 19 per entry for 2^51 - 1 entries
  plus 1 return and 2^10 per table, 2^10 per fold table, and cursor zeroing of 2 walkers and 2^36 relays
  (5 per 4 words over 2^36 words, 1 return).
- *node fixed:* 2^20 per node (code loading, stack, state block, input, the main routine's calls).
- *output:* 2^16, including the identity test of step 0.

W = total/20,000. The claim **120.534383** is log2 W rounded up at the sixth decimal; mode `ledger` prints
the unrounded value and the exact integers.

## 7. Memory per node (every mode checks each below 2^63 and the address ranges disjoint)

| Role | Count | Bytes |
| --- | ---: | ---: |
| walker | 2 | 2,162,465,043,685,311,488 |
| relay | 2^36 | 737,222,547,472,384 |
| home | 2^72 | 2,109,223,383,131,592,704 |

Peak total: 2^132.871418 bytes (rounded up).

Layout in machine.md section 4.

## 8. Checks and limitations

- `experiment.py` (with `progs.py`) runs every check in every mode, except the toy, which runs in modes
  `toy` and `all`; it exits with status 1 when any check fails and names the failed checks.
- Every mode compares every line above that carries a computed figure, and the tables of machine.md,
  with the computed values (mode `ledger` prints the comparison); every other numeral in this file and machine.md must be a listed structural
  constant, and the claim may print no other figures.
- It runs a mutation test: a table step one instruction longer or of 1000 instructions, a relay record one
  instruction longer, a swapped stated bound, a swapped memory figure, a swapped program count and a swapped
  claim figure must each make the comparison fail.
- It loads the shipped copy of the protected target profile, checks its sha256 against the protected value

b9c9b8f4a214f4c3327f27913caa78efa7bf260e7244e7dc6ba165e2172ce708

  and asserts p, r, x, g and the Fp12 basis. It probes every immediate, label and condition form the
  generator emits with out-of-range values and x18:

43 probes rejected, 20 in-range forms accepted.

- The success bound uses H1 and H2. Counts are executed counts of generated programs on the interpreter, not
  timings. The node machine (windows, stalling completion store, wake-only-on-store waits) is a modelling
  premise. The simulation sorts with Python (the radix listings run in mode `kernels`), and expands
  equal-key runs and computes the recheck with Python big integers (both charged in the rechecks row).
  Nothing runs at production size; sieve methods for the degree-12 field are expected to be cheaper.
