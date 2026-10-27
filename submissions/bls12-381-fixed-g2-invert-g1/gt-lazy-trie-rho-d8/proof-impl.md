# proof.md sections 0-i to 0-iv (implementation details)

This file holds sections 0-i to 0-iv of `proof.md`, split out to keep each file below the verifier's 64 KiB
per-file limit. Section numbers are unchanged; "section N" anywhere in the package refers to the same numbered
section, in `proof.md` or here.

### 0-i Counters and the work meter: widths and charged updates

Every counter of section 3.4 has a stated width and a priced update. The stuck-walk counter must
reach 20/theta = 20 * 2^72, which needs more than one word, so it is a
**2-word (128-bit) counter** in a register pair, incremented with `adds x, x, #1; adc y, y, xzr`
so that it cannot wrap. The report count must reach the cap D = 2^55 only: it is **one 64-bit
word**, incremented with `add x, x, #1` and compared with the cap by `cmp x10, x11; b.hs` (3 instructions, charged as 5 inside the report item of section 5), with no carry. The work meter must hold a whole run in instructions,
2^124.087 * 20,000 < 2^138.38 instructions, so it is a **3-word (192-bit) counter** in three
registers, with a headroom factor of at least 2^53.6 (192 - 138.38 bits).

**Pre-debit rule.** The meter counts the instructions charged since the walk phase began. Three registers
hold a fixed 192-bit threshold L = B_walk - A_max, where B_walk is the walk-phase budget of section 3.4 and
A_max = 2^25 instructions is above the charge of any single walk-phase action (a step at most 19,563 (a lazy step
with its flush), an escape 24,000, a checkpoint escape with its at most 64-step traversal charged as one action of
at most 24,000 + 64 x 19,563 + 65 x 64 = 1,280,192 instructions (asserted below A_max), a report 4,000, a switch 600, a walk start or restart (including the fresh start that follows every report) with its meter update and restart counter at
most 2^10.1 units plus 43 instructions, the meter adding that fixed allowance; `experiment.py` asserts each, and also the sum of a step, a report and an escape, and of a report and a fresh start, which bound what one iteration adds between two compares). **Before** every charged action the
meter is compared with L; the action runs only if the meter is below L, and **after** it the meter adds
the action's charge. So the meter stays below L + A_max = B_walk. When the compare finds the meter at or above L, the walk phase
ends and nothing more is charged to the meter: that failed compare, the hand-over to the matcher and all
later work are paid from the stopping allowance of section 3.4.

Per step, every kind, the meter is charged **8 instructions** and the counter 5 in the table of section 5.
At the top of the step, before the key kernel runs, the meter is compared with L by the full three-word compare
`cmp x0, x5; sbcs xzr, x1, x6; sbcs xzr, x2, x7; b.hs end_walk_phase` (4 instructions), which branches exactly
when meter >= L as 192-bit unsigned integers: after a borrow chain the carry flag C is the full-width "no borrow"
result, so a carry-only condition (`b.hs`, `b.lo`) is exact, while Z reflects only the last word, so no condition that
reads Z (`b.hi`, `b.ls`, `b.eq`, `b.ne`) follows any chained compare here. `flags_check.py`, run by `experiment.py`,
executes this and every other multiword listing of 0-i, 0-ii and 0-iv as written under the AArch64 NZCV rules on edge and
random cases. **This exact compare runs before every step of every walk**, with no cheaper pre-test in front of it. The step's kind is known from the node depth before the kernel runs, and
the step's charge (1,734 instructions for a root step, 5,137 for a lazy step, 19,563 for the lazy step whose child
is a leaf and which therefore flushes: exactly the prices the ledger of section 5 charges) is materialised in a register (`mov x4, #imm`, 1 instruction) and added
with `adds x0, x0, x4; adcs x1, x1, xzr; adc x2, x2, xzr` (3 instructions). A report, an escape and a walk
start or restart (including the fresh start after a report) are separate charged actions under the same rule; each pays its own compare and add
inside its own charge (section 5). The 2-word stuck-walk
counter (low word x8, high word x9) costs 5 instructions per step:
`adds x8, x8, #1; adc x9, x9, xzr; mov x16, #5120; cmp x9, x16; b.hs restart`. Only the high word is compared
because the threshold is an exact multiple of 2^64: 20/theta = 20 * 2^72 = 5120 * 2^64, so the counter
has reached it exactly when x9 >= 5120. 5120 is not an immediate of `cmp` (12 bits, optionally shifted left by
12), so it is moved into the scratch register x16 on every step, a charged instruction.

The stuck-walk counter is **reset at every walk start, restart and distinguished-point report**; since every
report ends its trail and is followed by a fresh start, it counts the steps of the current trail, and the
stuck-walk rule fires at most once per trail (a trail is a report gap), not once per walk. The
reset at a report is two moves into the register pair (`mov x8, xzr; mov x9, xzr`, 2 instructions); it is charged 4 instructions inside the
4,000-instruction report item of section 5 (a listed line of the item; the total stays 4,000),
at all 2^55 reports (the fresh start that follows resets it again inside its own 40, a redundant charged reset). At a walk start or restart the reset is part of the 40 instructions of meter bookkeeping
below, and the fresh start itself is charged in section 5 (at most 2^10.1 units each: 2^46.1 for the 2^36
stuck-walk restarts of the cap, at most 2^65.1 for the at most 2^55 fresh starts after reports, below 2^65.11 with their bookkeeping).

**Schedule.** All M = 2^20 walks run on one processor, round-robin: it advances one walk for a chunk of
K = 2^63 consecutive step iterations, counted down in one register (`subs x12, x12, #1; b.eq switch`, 2 of the
itemised 17 of loop control of section 5), then a **switch** stores that walk's state (at most 2 KB, 256 stores), loads the
next walk's (256 loads), sets the index mod M and its address (4) and the chunk counter (1), and pays its own
compare and add (8) and 75 for branch, call and margin: **600 instructions**, a charged action. A restart
replaces the walk's state in place and does not end the chunk (the fresh start after a report likewise); the phase opens with a switch. There are at most n_max/K + M + 2^36 + 1 < 2^62.72
switches (n_max of 0-iii; an over-count), **below
2^57.65 units**, a term of the ledger (section 5). The meter is that processor's one 192-bit register triple for the whole
walk phase; no switch saves or loads it and nothing runs concurrently, so the pre-debit compare sees all work
charged before every action. Each start or restart pays 40 instructions (compare and add 8, counter reset 2, loads and moves 30):
2^20 * 40 < 2^11.04 units, 2^36 * 40 < 2^27.04 units for stuck-walk restarts and 2^55 * 40 < 2^46.04 units for the fresh starts after reports, all in section 5's non-walk bucket.

### 0-ii The report log and the matcher: one 64-bit address space, no initialisation, worst-case matching

**Report log.** Reports go into an append-only log: an array of exactly **REPORT_CAP = 2^55 record
slots of 80 bytes** with no pointer field, **2^55 * 80 = 2^61.322 bytes**, which is the whole memory of the
store. A record holds 16 bytes of fingerprint F (the 112 bits of F, zero-extended) and the reporting trail's **start pair**, 32 bytes for alpha0 and 32 for beta0,
each 32-byte field as four 64-bit words, most significant first; both are residues in [0, r) as the sampler
returned them (section 0-iv) and are copied from the walk state, where they have been kept since the trail's fresh
start. No exponent of the reporting element is computed or stored. The slot of a report is the value of the one-word report counter
c before its increment, at byte offset 80 c = 64 c + 16 c (`add x2, x0, x1, lsl #6; add x2, x2, x1, lsl #4`).
Every byte offset is below 2^61.322 < 2^63, one native 64-bit register with the sign bit clear; there is no
sharding and no multiword address (`experiment.py` part D asserts each bound).
A report compares c with the cap; at c = 2^55 the walk phase ends and **the run fails**, a cap failure
bounded in section 4 (below exp(-2^52) under H2). There is no per-walk report cap: the global cap and its failure rule are
the only capacity rule.

**Report charge (4,000 instructions, at most 2^55 reports, 2^52.68 units).** Itemised as in the parent: the
fingerprint F = bits 136..247 of the stored K (K_2 is in a register from the DP test; `ldr` of K_3, `extr` for bits
136..199, `ubfx` for bits 200..247) 3; loading the trail's start pair (four ldp) 4; the pre-debit compare and the charge add 8; the report counter
compare, branch and increment 5; the slot address base + 80 c 2; appending the 80-byte record 12; resetting the 2-word
stuck-walk counter 4; the operating system's zero-fill of the log's fresh 2 MiB pages 13 (below); margin for loads,
moves, call and return 3,949. No keys and no products mod r are computed at a report.

**Nothing reads a slot that was not written.** The sort, the dedup and the scan below cover the slots
[0, c) only, every one of which a report wrote. So the log needs no clearing pass, no used bit and no free-slot
test. The walk states, counters and histories
of section 0-iv are still initialised, and charged, before they are read.

**Fresh pages.** The log's 2^61.322 bytes are reserved in 2 MiB pages and never zeroed by the algorithm. If the
operating system zero-fills a page when a report first writes into it, the page costs at most 2 MiB / 16 =
131,072 sixteen-byte stores, at most one loop-control instruction per store and one page fault, charged at
65,536 instructions for its entry, page-table update and return. A page holds 26,214 whole 80-byte records, so
this is (2 x 131,072 + 65,536) / 26,214 = 12.5002 instructions per record (5.0001 for the stores), charged as
**13 inside the report's 4,000** (section 5; `experiment.py` asserts it). The first n reports touch at most
80n/2^21 + 1 pages, so what 13n leaves is at most one page, 327,680 instructions (16.4 units), inside the
setup bound's slack (setup items sum below 2^30.11). Every other array is written before it is read; its at
most 1,045 fresh 2 MiB pages are zero-filled the same way, below 2^14.1 units.

**Matching runs once, after the walk phase**, on the n <= 2^55 written records, in place:

* *Heapsort of the whole 80-byte record* (F first, then alpha0, then beta0): Floyd's heap construction, then
  n - 1 extractions, with one 80-byte temporary and no other memory. The construction performs at most n
  sift-down levels and the extractions at most sum_{k=2..n} floor(log2 k), so at most
  **54 x 2^55 + 57 levels for n = 2^55** (`experiment.py` sums the formula exactly in power-of-two blocks). A
  record is ten words w0..w9 (offset 8k, w0 most significant, each field high word first); both compares
  are carry-only chains from w9. Child compare (20 `ldr`, `cmp` and 9 `sbcs`, 11 `csel`; x12 := the larger child, a tie
  picks left): 41. Parent compare (10 `ldr`, `cmp`, 9 `sbcs` from w9, `b.hs`): 21. A swap moves one record
  (10 loads, 10 stores, 20); index arithmetic and loop control 50: **132 instructions** (the exact listings
  are run in `flags_check.py`). Each of the at most n/2 + n - 1 outer iterations is charged 40 more. That is 7,188 instructions per record, **below 2^53.53 units** for 2^55 records.
  Heap indices are below 2^56 and their byte offsets below 2^63, so each fits one word.
* *Duplicates*: one adjacent 80-byte comparison per record (100 instructions) drops a record identical to
  its predecessor and marks where a run of equal F begins. (A trail reports once, at its end, so two identical
  records need two trails with the same start pair.)
* *Replay.* The exponents of a record are recovered by **replaying its trail**: start at y = g^alpha0 z^beta0 with
  (alpha, beta) = (alpha0, beta0) and execute the walk rule of section 3 iteration by iteration, now also updating
  (alpha, beta) mod r on every step (doubled at an escape's squaring, plus the scalar-table
  entry (a_h p^t, b_h p^t) on every step; the replay also keeps the pair of the remembered previous element,
  so a 2-cycle escape that squares that element doubles its pair), until the first distinguished class; there compute the 12 keys for
  kappa and the canonical pair (alpha^, beta^) = p^kappa (alpha, beta). The walk rule reads only the trail's own
  state (its base buffers, trie node, its two-step K history, its checkpoint and its step count mod 64, all set at the fresh start,
  section 0-iv), never the schedule, the meter or another walk, so the replay repeats the trail exactly and ends at
  the reported class (`experiment.py` part B replays every matched pair and asserts that the replayed F equals the
  record's). A reported trail has at most 20/theta steps (otherwise the stuck-walk rule would have restarted it
  before it reported); a 2-cycle escape needs two steps since the previous escape or the start, and a checkpoint
  escape at least one (an escape resets the two-step history and the checkpoint, as in part B's toy), so a trail has at
at most 2 * 20/theta + 3 iterations. A replayed iteration repeats at most one step (at most 19,563 instructions plus a
  64-instruction exponent update) and one escape; the costliest escape, a checkpoint escape with its at most 64-step
  traversal and all its exponent updates, 24,000 + 64 x 19,563 + 65 x 64 = 1,280,192 instructions plus 64 x 64 of exponent updates, is below A_max = 2^25
  (`experiment.py` asserts it). So an iteration is below 2 * 2^25 instructions, and one replay, with its start
  exponentiations (at most 1,021 units, far inside the per-iteration slack), is below (2 * 20/theta + 3) * 2^26
  instructions, **below 2^89.04 units**, whatever the walk did.
* *Pairs*: inside each maximal run of records with equal F, every pair is taken in turn. Both records of the pair
  are replayed. A pair with beta^1 = beta^2 gives no candidate; otherwise the candidate
  x = (alpha^1 - alpha^2)/(beta^2 - beta^1) mod r is checked by g^x = z. A pair costs at most two replays, two key
  computations with their two products mod r (at most 4,000 instructions each) and 600 units (counter, test, one
  inversion, the check): **below 2^90.04 units**. The matcher **halts at the first pair whose candidate verifies**.
  Two distinct records of one class are two trails with different start pairs ending at the same class, so
  g^alpha^1 z^beta^1 = g^alpha^2 z^beta^2; they verify unless beta^1 = beta^2 (probability about 1/r for independent
  fresh starts, section 4). So every pair charged before the halting one is a false-key pair of two distinct
  classes or such a degenerate pair, and the pairs charged are at most those plus one.
* *Hard pair budget*: a pair counter stops the matcher, and **the run fails**, when it would charge pair
  2^16 + 1. By linearity of expectation, with no independence assumption, the expected number of false-key
  pairs is at most (2^55)^2/2 = 2^109 record pairs times 2^-112 (the chance under H2 that two distinct classes share
  the 112 bits of F): 2^-3. Markov's inequality bounds the overflow probability by 2^-19, a loss term of section 4.

Only the stuck-walk counter (two words) and the work meter (three words) exceed 64 bits; they are
arithmetic counters, not addresses (section 0-i). The trie's transient leaf id (below 2^65.36, section 3) also
lives in a two-register pair, only between its computation at a depth-7 push and the flush's banked read; it is
an index into the banked leaf array, never a byte offset, and is never stored.

### 0-iii One restart cap

The restart cap is **2^36 stuck-walk restarts** (the fresh start that follows every report is not a stuck-walk
restart and is bounded by the report cap 2^55). Derivation. A walk is restarted
when its stuck-walk counter reaches 20/theta = 20*2^72 steps, and the counter is reset at every report and fresh start
(section 0-i). The rule therefore fires at most once per **trail** (a gap: the steps from a fresh start to the trail's
first distinguished point), and each trail is a separate chance. A trail ends at its
first distinguished point; each step is one at density theta = 2^-72, so a trail exceeds
20/theta steps with probability (1 - theta)^(20/theta) <= e^-20 < 2.07e-9. The number of trails is the number of
reports plus restarts, expected about theta times the steps. The walk phase runs until its budget is spent (section 3.4); every block of k iterations, with the escape or the
report and fresh start that ends it, is charged at least k c, so the walk phase takes at most
n_max = B_walk/(c C) + 6M < **2^125.71** iterations (section 3.4), and the expected number of stuck-walk restarts over
the whole run is at most n_max theta e^-20 < **2^24.85** (each walk has about 2^34 trails, so this is not a
per-walk figure). A walk caught in a useless closed walk of more than 64 steps
with no distinguished class (rate below 2^-105.6 per step, charged as 2^-104, section 3.3 and part F) is the other source: at most
n_max 2^-104 < **2^21.71** in expectation.
A cycle that contains a distinguished class cannot trap a trail: the trail ends at that class the first time it
reaches it, is recorded, and the walk restarts from a fresh start. So, as in the parent, there is no failed-run
event for long cycles through distinguished classes and no loss term for them.
The expected number of stuck-walk restarts of both kinds is below **2^25.01**, a factor 2^10.99 below the cap 2^36.
Under H2 the restart events of distinct trails are independent, so the Chernoff bound puts the probability of more
than 2^36 restarts below exp(-2^38) (section 4); reaching the cap fails the run. The
loss sum of section 4 charges the cap's own worst case as steps spent on restarts, not as their expected count: the
2^36 restarts' one-full-allowance worst case is 2^36 * 20/theta < 2^112.33 steps; part D computes its fraction of
the walk's iteration maximum in 60-digit decimal arithmetic, 9.37e-5, and subtracts that value; no other restart term is subtracted. The meter supports the cap: 2^36
restarts cost 2^36 * 40 < 2^41.33 instructions of meter bookkeeping (below 2^27.04 units, section 0-i), 2^36 fresh starts of at most 2^10.1
units each (2^46.1 units, section 5) and, at worst, one full stuck-walk allowance each, 2^36 * 20/theta < 2^112.33 steps,
which section 5 charges in the step bucket and whose nominal cost is below 2^110.81 units. The cap
therefore neither breaks the meter nor adds to the success allowance.

### 0-iv Initial values

Nothing the design reads is assumed to start at zero or at a known value; each is set and charged:

* the report log: not initialised, because nothing reads a slot before a report writes it (section 0-ii); an
  operating system's zero-fill of its fresh 2 MiB pages is charged, 13 instructions per record inside the 4,000;
* the 3-word work meter (set to zero), the 3-word threshold triple (loaded with L = B_walk - A_max), the
  one-word report, restart and pair counters and walk index (each set to zero; each switch sets the chunk counter): at
  most 32 instructions once, before the walks, inside the at most 2^31 units of the setup items (section 5);
* the restart counter: incremented and compared with the cap 2^36 at each restart (`add`, `cmp`, `b.hs`,
  3 instructions), 2^36 * 3 / 20,000 below 2^23.30 units at the cap, a named term of the non-walk sum
  of section 5; the pair counter: the same 3 instructions per pair, inside its at most 600 units;
* the kernels' constant pools: set in setup (section 5);
* the stuck-walk counter: set to zero at each walk start (not once before the walks), and reset at every
  restart, report and fresh start (section 0-i);
* the multiplier table T[h][t] and the scalar tables a_h p^t, b_h p^t mod r (read only by the replays of section 0-ii): computed in setup, the table item
  of section 5 (24 * (1,021 + 11) units, 48 sampler draws and 2 * 288 products mod r at 800
  instructions each, and the 288 x 54 Karatsuba pre-sums of the table entries, each a sum of at most 8
  reduced coordinates at 200 instructions; below 2^14.7 units);
* the prefix-product trie of section 3: built in setup level by level, each node from its parent by one `mul.s` call
  with the table entry T[h][t], then its 42 Karatsuba pre-sums and (internal nodes) its vector gamma, charged one unit
  per node (16,068 instructions plus the banked addressing, itemised in part D), 2^57.19 + 2^65.36 nodes held in
  banks of 2^32 nodes (section 3); its 2^76.71 bytes are written before they are read, and an operating system's
  zero-fill of its fresh 2 MiB pages is charged at 0.15625 instructions per byte (2^58.85 units); below 2^65.393
  units with the rest of setup; the bank tables (2^36.36 bytes) are written in setup before the walks;
* the sampler of section 3: the analysis is carried out in the ideal model in which the bit stream is an
  ideal random source; the implementation realizes that ideal mapping with a ChaCha20 generator (RFC 8439's
  round function), and this package assumes, as a disclosed heuristic of the same kind as H1 and H2, that
  the realization is indistinguishable from the ideal mapping. No proved distributional or independence
  bound on the generated scalars is claimed: the key and nonce define a finite family of streams, and the
  success probability is a heuristic over the ideal mapping, not a total-variation statement about that
  family. One attempt takes 255 bits of the stream (one 64-byte block) as X, accepts iff X < r by the
  four-word carry chain `cmp x0,x4; sbcs xzr,x1,x5; sbcs xzr,x2,x6; sbcs xzr,x3,x7; b.lo accept`, and
  otherwise retries (20 instructions with the retry counter and loop branch). A block is charged 2,300
  instructions, an itemised upper bound: 80 quarter-rounds of 12 instructions (4 additions, 4 XORs and 4
  rotates) = 960, the 16-word feed-forward add 16, loading and storing the 16-word state 32, the 64-bit
  block counter's `adds`/`adc` 2, address set-up and loop control 10, and 250 for register spills, the
  prologue and return, 1,270 in all; the charge keeps 2,300, a 1.8x margin. The counter is 64 bits (the
  original ChaCha20 64-bit-counter/64-bit-nonce state layout; RFC 8439's 96-bit-nonce instantiation has
  only a 32-bit counter, which would wrap). A run draws at most 48 table scalars and two per walk start,
  stuck-walk restart and fresh start after a report, 48 + 2 (2^20 + 2^36 + 2^55) < 2^56.0001 samples, with
  at most 2^7 attempts (blocks) each: below 2^63.0001 < 2^64 blocks, so the 64-bit counter cannot wrap and its
  update is inside the block's 2,300 (this is why R_MAX is 128 here, not the parent's 256: with the fresh
  starts after reports, 256 attempts per sample could exceed 2^64 blocks). It charges a fixed **2,320 instructions per attempt and R_MAX = 128
  attempts, 296,960 instructions = 14.848 units per sample**, independent of the draw. X is 255 bits and
  2^254 < r < 2^255, so acceptance is above 1/2; after 128 rejects the sampler returns the fixed residue
  1 (the exhaustion rule). Under the ideal-source model the exhaustion event has probability at most
  2^-128 per sample and at most 2^-71.99 over a run's below 2^56.0001 samples (the union bound); it is
  carried as an exhaustion allowance of 2^-71.99 in section 4's budget and changes no digit. The
  256-bit key, the 64-bit nonce (the operating system's random source), the 64-bit block counter and the
  retry counter are initialised in setup, 2^12 instructions inside the 2^31-unit bound. The table draws
  two scalars per multiplier and each walk start, restart and fresh start two; nothing else samples. `experiment.py`'s
  `sample_below` runs the same compare, cap and rule over its seeded bit source (part S);
* each walk's state y and its trail's start pair (alpha0, beta0): a fresh start y = g^alpha0 z^beta0 from fresh uniform
  (alpha0, beta0) at every walk start, stuck-walk restart and after every report (at most 2^10.1 units each, section 5; two
  exponentiations of at most 1,021 units, two sampler draws and 43 instructions of bookkeeping take 1,050.7, inside the bound);
  the walk keeps (alpha0, beta0) unchanged for the record and carries no running exponents;
* each walk's histories, the two remembered K slots, the checkpoint K, the node (root), the depth and the step count mod 64 of section 3.3:
  set at every walk start, restart and fresh start, the K fields to the six-limb all-ones value 2^384 - 1, which no stored K
  (a residue below p < 2^381) equals, so the first steps cannot read a false repeat, and the count to zero.
  These few dozen instructions sit inside the fresh start's 2^10.1 units, whose two exponentiations take at
  most 1,021 units and, with the two draws and 43 instructions, leave more than 46 units (asserted: more than 16).

