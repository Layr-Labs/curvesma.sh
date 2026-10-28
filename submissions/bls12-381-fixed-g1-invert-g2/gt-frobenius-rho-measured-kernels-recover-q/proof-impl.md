# proof.md sections 0-i to 0-iv (implementation details)

This file holds sections 0-i to 0-iv of `proof.md`, split out to keep each file below the verifier's 64 KiB
per-file limit. Section numbers are unchanged; "section N" anywhere in the package refers to the same numbered
section, in `proof.md` or here.

### 0-i Counters and the work meter: widths and charged updates

Every counter of section 3.4 has a stated width and a priced update. The stuck-walk counter must reach
28/theta = 28 * 2^86: a **2-word (128-bit) counter** in a register pair, which cannot wrap (28 * 2^86 < 2^128). The
report count reaches the cap D = 2^55 only: **one 64-bit word**, `add x, x, #1` and `cmp x10, x11; b.hs`
(3 instructions, charged as 5 inside the report item of section 5), no carry. The work meter must hold a whole run,
2^123.826 * 20,000 < 2^138.114 instructions: a **3-word (192-bit) counter** in three registers, headroom at least
2^53.88 (192 - 138.114 bits).

**Pre-debit rule.** The meter counts the instructions charged since the walk phase began. Three registers hold a
fixed threshold L = B_walk - A_max, B_walk the walk-phase budget of section 3.4 and A_max = 2^25 instructions, above
the charge of any single walk-phase action (a step at most 16,342, an escape 24,000, a checkpoint escape with its at
most 64-step traversal at most 24,000 + 64 x 16,342 + 65 x 64 = 1,074,048, a report 4,000, a switch 600, a walk start
or restart, the fresh start after a report included, with its bookkeeping at most 2^10.1 units plus 43 instructions;
`experiment.py` asserts each, and the sums step + report + escape and report + fresh start, which bound what one
iteration adds between two compares). **Before** every charged action the meter is compared with L and the action
runs only if the meter is below L; **after** it the meter adds the action's charge. So the meter stays below
L + A_max = B_walk. When the compare finds the meter at or above L the walk phase ends and nothing more is charged to
the meter: that compare, the hand-over to the matcher and all later work are paid from the stopping allowance of
section 3.4.

Per step, every kind, the meter costs **8 instructions** and the counter 5 (section 5). At the top of every step,
before the key kernel and with no cheaper pre-test in front of it, `cmp x0, x5; sbcs xzr, x1, x6; sbcs xzr, x2, x7;
b.hs end_walk_phase` branches exactly when meter >= L as 192-bit unsigned integers: after a borrow chain the carry
flag C is the full-width "no borrow" result, so a carry-only condition (`b.hs`, `b.lo`) is exact, while Z reflects
only the last word, so no condition that reads Z (`b.hi`, `b.ls`, `b.eq`, `b.ne`) follows any chained compare here.
`flags_check.py`, run by `experiment.py`, executes every multiword listing of 0-i, 0-ii and 0-iv as written under the
AArch64 NZCV rules on edge and random cases. The step's kind is known from the node depth before the kernel runs; its
charge (1,789 for a root step, 3,589 for a lazy step, 16,342 for the lazy step whose child is a leaf and which
therefore flushes: the ledger's prices) is set by `mov x4, #imm` and added by `adds x0, x0, x4; adcs x1, x1, xzr;
adc x2, x2, xzr`. A report, an escape and a walk start or restart are separate charged actions under the same rule,
each paying its own compare and add inside its charge (section 5). The stuck-walk counter (low word x8, high word
x9) adds 1 at every iteration: `adds x8, x8, #1; adc x9, x9, xzr; mov x16, #0x7000000; cmp x9, x16; b.hs restart`.
Only the high word is compared because 28/theta = 0x7000000 * 2^64 exactly, so the counter has reached it exactly
when x9 >= 0x7000000; that value is not a `cmp` immediate (12 bits, optionally shifted by 12) but is one MOVZ
(0x700 << 16, 16 bits at shift 16), moved into x16 on every step, a charged instruction.

**Credits at escapes.** An escape ends an *episode* (the iterations since the trail's start or previous escape) and
lowers the counter by at most the episode's revisits (0-iii). A 2-word **mark** holds the counter at the episode's
start ([x3, #16]; the step count mod 64 is at [x3, #8]). A 2-cycle escape credits 2 (`subs x8, x8, #2;
sbc x9, x9, xzr; stp x8, x9, [x3, #16]`; its episode has at least 3 iterations). A checkpoint escape of a cycle of
length l = ((s - 1) mod 64) + 1, s the step count mod 64, sets the counter to max(counter - 64, mark + l), a credit of
min(64, counter - mark - l), by the 18-instruction listing `CREDIT_CHECKPOINT` of `flags_check.py` (carry-only
conditions). Both store the new counter as the mark, so an escape leaves the counter at least 1 above the mark and
credits at most 64. The checkpoint escape's traversal to the cycle's minimum is part of the escape action and does not
touch the counter. `flags_check.py` runs both listings and the reset on edge cases (a low word below 64, where
counter - 64 wraps and is not selected; a credit of exactly 64; borrows and carries across the word boundary; l = 1,
l = 64; mark 0) and random ones; the 18 are priced in the escape's 24,000. The counter and mark are **reset at every
walk start, restart and distinguished-point report** (`mov x8, xzr; mov x9, xzr; stp xzr, xzr, [x3, #16]`), so the
rule fires at most once per trail (a report gap). At a report the reset is charged 4 inside the 4,000-instruction
report item of section 5 (the fresh start that follows resets again inside its own 40); at a walk start or restart it
is part of the 40 below. The fresh start itself is charged in section 5 (at most 2^10.1 units each: 2^30.1 for the
2^20 restarts, 2^65.1 for the at most 2^55 fresh starts after reports, below 2^65.11 with their bookkeeping).

**Schedule.** All M = 2^20 walks run on one processor, round-robin: it advances one walk for a chunk of K = 2^63
consecutive step iterations, counted down in one register (`subs x12, x12, #1; b.eq switch`, 2 of the itemised 17
of loop control of section 5), then a **switch** stores that walk's state (at most 2 KB, 256 stores), loads the
next walk's (256 loads), sets the index mod M and its address (4) and the chunk counter (1), and pays its own
compare and add (8) and 75 for branch, call and margin: **600 instructions**, a charged action. A restart or a
fresh start replaces the walk's state in place and does not end the chunk; the phase opens with a switch. There are
at most n_max/K + M + 2^20 + 1 < 2^62.712 switches (an over-count), **below 2^57.653 units**, a term of the ledger
(section 5). The meter is that processor's one 192-bit register triple for the whole walk phase; no switch saves or
loads it and nothing runs concurrently, so the pre-debit compare sees all work charged before every action. Each
start or restart pays 40 instructions (compare and add 8, counter and mark reset 3, loads and moves 29):
2^20 * 40 < 2^11.04 units for walk starts and for stuck-walk restarts, 2^55 * 40 < 2^46.04 units for the fresh
starts after reports, all in section 5's non-walk bucket.

### 0-ii The report log and the matcher: one 64-bit address space, no initialisation, worst-case matching

**Report log.** Reports go into an append-only log of exactly **REPORT_CAP = 2^55 record slots of 80 bytes**, no
pointer field, **2^55 * 80 = 2^61.322 bytes**, the whole memory of the store. A record holds the 128 bits of F
(16 bytes) and the reporting trail's **start pair**, alpha0 and beta0 as 32-byte fields of four 64-bit words, most
significant first: residues in [0, r) as the sampler returned them (0-iv), copied from the walk state, which keeps
them from the trail's fresh start. No exponent of the reporting element is computed or stored. A report's slot is
the one-word report counter c before its increment, at byte offset 80 c = 64 c + 16 c
(`add x2, x0, x1, lsl #6; add x2, x2, x1, lsl #4`), below 2^61.322 < 2^63: one native 64-bit register with the sign
bit clear, no sharding and no multiword address (`experiment.py` part D asserts each bound). A report compares c
with the cap; at c = 2^55 the walk phase ends and **the run fails**, a cap failure charged as its own loss, at most
2^-15.28 by Markov's inequality (proof.md sections 3.4 and 4); the global cap is the only capacity rule.

**Report charge (4,000 instructions, at most 2^55 reports, 2^52.68 units).** Itemised as in the sibling package: the
fingerprint F = bits 150..277 of the stored K (`ldp K2,K3`, `ldr K4`, `extr F_lo, K3, K2, #22` for F bits 0..63 and `extr F_hi, K4, K3, #22` for bits 64..127) 5; loading the trail's start pair (four ldp) 4; the pre-debit compare and the charge add 8; the report counter
compare, branch and increment 5; the slot address base + 80 c 2; appending the 80-byte record 12; resetting the
stuck-walk counter and its mark 4; the operating system's zero-fill of the log's fresh 2 MiB pages 13 (below); margin for loads,
moves, call and return 3,947. No keys and no products mod r are computed at a report.

**Nothing reads a slot that was not written.** The sort, the dedup and the scan below cover the slots
[0, c) only, every one of which a report wrote. So the log needs no clearing pass, no used bit and no free-slot
test. The walk states, counters and histories
of section 0-iv are still initialised, and charged, before they are read.

**Fresh pages.** The log's 2^61.322 bytes are reserved in 2 MiB pages, never zeroed by the algorithm. If the
operating system zero-fills a page when a report first writes into it, the page costs at most 2 MiB / 16 = 131,072
sixteen-byte stores, at most one loop-control instruction per store, and one page fault charged 65,536 instructions.
A page holds 26,214 whole 80-byte records, so this is (2 x 131,072 + 65,536) / 26,214 = 12.5002 instructions per
record (5.0001 for the stores), charged as **13 inside the report's 4,000** (asserted). The first n reports touch at
most 80n/2^21 + 1 pages, so what 13n leaves is at most one page, 327,680 instructions (16.4 units), inside the setup
bound's slack (2^52.21 units against setup items below 2^52.2098, a slack above 2^39.5 units). Every other array is
written before it is read; its at most 4,100 fresh 2 MiB pages (below 2^33 bytes, section 5) are zero-filled the same
way, below 2^16.1 units, inside the same slack (`experiment.py` asserts both pages' sum inside it).

**Matching runs once, after the walk phase**, on the n <= 2^55 written records, in place:

* *Heapsort of the whole 80-byte record* (F first, then alpha0, then beta0): Floyd's heap construction, then n - 1
  extractions, with one 80-byte temporary and no other memory: at most n sift-down levels for the construction and
  sum_{k=2..n} floor(log2 k) for the extractions, so at most **54 x 2^55 + 57 levels for n = 2^55** (`experiment.py`
  sums the formula exactly in power-of-two blocks). A record is ten words w0..w9 (offset 8k, w0 most significant,
  each field high word first); both compares are carry-only chains from w9. Child compare (20 `ldr`, `cmp` and 9
  `sbcs`, 11 `csel`; x12 := the larger child, a tie picks left): 41. Parent compare (10 `ldr`, `cmp`, 9 `sbcs` from w9,
  `b.hs`): 21. A swap moves one record (10 loads, 10 stores, 20); index arithmetic and loop control 50: **132
  instructions** (the exact listings run in `flags_check.py`). Each of the at most n/2 + n - 1 outer iterations is
  charged 40 more: 7,188 instructions per record, **below 2^53.53 units** for 2^55 records. Heap indices are below
  2^56 and their byte offsets below 2^63, so each fits one word.
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
  record's). A trail has at most I_max = 65 x 28/theta iterations and at most 28/theta escapes (0-iii). A replayed
  iteration repeats one step (at most 16,342 instructions plus a 64-instruction exponent update); an escape costs at
  most 24,000 + 64 x 16,342 + 65 x 64 = 1,074,048 (a checkpoint escape with its traversal; below A_max = 2^25,
  asserted) plus 64 x 64 of exponent updates. So one replay, with its fresh start (at most 2^10.1 units), is below
  28/theta (65 x 16,406 + 1,078,144) instructions plus that start, **below 2^97.552 units** (asserted), whatever the walk
  did; the stopping allowance keeps its charge per replay, (2 x 28/theta + 3) x 2^26 instructions (below 2^103.52 units).
* *Pairs*: inside each maximal run of records with equal F, every pair is taken in turn and both records are
  replayed. A pair with beta^1 = beta^2 gives no candidate; otherwise x = (alpha^1 - alpha^2)/(beta^2 - beta^1) mod r
  is checked by g^x = z. A pair costs at most two replays, two key computations with their products mod r (at most
  4,000 instructions each) and 600 units (counter, test, one inversion, the check): **below 2^104.52 units**
  (asserted). The matcher **halts at the first pair whose candidate verifies**. Two distinct records of one class are
  two trails with different start pairs ending at that class, so g^alpha^1 z^beta^1 = g^alpha^2 z^beta^2: they verify
  unless beta^1 = beta^2 (about 1/r for independent fresh starts, section 4). So the pairs charged before the halting
  one are false-key pairs of two distinct classes or such degenerate pairs.
* *Hard pair budget*: a pair counter stops the matcher, and **the run fails**, when it would charge pair 2^4 + 1. By
  linearity of expectation, with no independence assumption, the pairs charged before the halting one number at most
  (2^55 choose 2) (2^-128 + 2^-254) < 2^-19 in expectation: (2^55 choose 2) = 2^109 - 2^54 record pairs, each with
  chance 2^-128 under H2 that two distinct classes share the 128 bits of F, plus 2^-254 for beta^1 = beta^2 (a fresh
  start pair is uniform); the product is 2^-19 - 2^-74 + 2^-145 - 2^-200 < 2^-19 (asserted, with r > 2^254). Markov's
  inequality, valid for any dependence, gives P(the budget is reached) <= 2^-19 / 2^4 = 2^-23 (part D asserts and
  prints it).

Only the stuck-walk counter and its mark (two words each) and the work meter (three words) exceed 64 bits; they are
arithmetic counters, not addresses (section 0-i).

### 0-iii One restart cap: the trail tail in first visits

The cap is **2^20 stuck-walk restarts** (a fresh start after a report is not one). A walk restarts when its counter
reaches 28/theta = 28*2^86, at most once per **trail** (from a fresh start to its first distinguished class; 0-i).

**What the counter counts.** A *revisit* is a visit to a class some trail visited before; a first visit tests a DP
bit no visit tested. The step maps classes to classes (section 3.1), so inside an episode a revisit means a cycle.
If its length l is at most 64 the episode ends at its escape: a 2-cycle (or a 1-cycle, seen one iteration later)
after at most 2 revisits; a cycle of 3 <= l <= 64 after at most 64, since the first checkpoint stored on it, at most
63 steps after the walk entered it, is met again l steps later, and the episode also holds the cycle's l first
visits. So the credit, 2 or min(64, counter - mark - l), covers every revisit of the episode. Any other revisit (a
merge into another trail's path, a class of an earlier episode, a cycle longer than 64) is *uncredited*. Until its
first uncredited revisit a trail's counter is at most its first visits after each escape, and at most 64 more in
between; a trail restarted with none has made at least 28/theta - 64 first visits, none distinguished. Under H2 the
bits of K are a random function of the class, so a class's DP bit at its first visit is a new uniform draw whatever
came before; by the chain rule each trail, whatever the history before it, restarts that way with probability at
most x = (1 - theta)^(28/theta - 64) <= e^(-28 + 2^-80) < 6.92e-13. Part B runs this counter in the toy and checks it
against per-trail sets of visited classes (evidence.md).

**Expected restarts** (linearity, no independence): with no uncredited revisit, at most x times the expected trails
(2^20 walk starts, at most 2^20 restarts, the reports; below 2^39.712), below 2^-0.68; useless closed walks longer
than 64 with no distinguished class (below 2^-122.6 per step, charged 2^-121; section 3.3, part F), at most n_max
2^-121 < 2^4.712 over the at most n_max < **2^125.712** iterations (section 3.4); any other uncredited revisit (one
restart at most) is a collision of the class sequence, at most n_max^2/N < 2^0.151 under H1, read as "a class the
walk reaches by a step, an escape or a fresh start equals a given earlier class with probability at most 1/N". A
cycle through a distinguished class traps no trail. In all below **2^4.804**, a factor above 2^15.196 (37,550) below
the cap; Markov's inequality, which assumes no independence, gives P(2^20 restarts) <= 2^-15.196, a loss line of
proof.md section 4.

**Trail length.** An episode's iterations are at most its credit (at most 64) plus the counter's rise over it (at
least 1): a trail has at most I_max = 65 x 28/theta iterations and 28/theta escapes, whatever the walk does. Credited
iterations are revisits in fruitless cycles, at most the steps X charges for them (2 per 2-cycle, 64 + l per
checkpoint escape; section 3.3), so a restarted trail has at most 28/theta iterations outside X; the cap's worst
case, 2^20 * 28/theta < 2^110.81 such steps (asserted), is charged as a fraction below 3.27e-5 of n_max (section 4).
The ledger pays 2^20 restarts' bookkeeping 2^20 * 40 < 2^25.33 instructions (below 2^11.04 units), their fresh starts
(2^30.1 units) and 28/theta steps each (below 2^109 units; all asserted).

### 0-iv Initial values

Nothing the design reads is assumed to start at zero or at a known value; each is set and charged:

* the report log: not initialised (nothing reads a slot before a report writes it, 0-ii); a zero-fill of its fresh
  2 MiB pages is charged, 13 instructions per record inside the 4,000;
* the 3-word meter (zero), the threshold triple (L = B_walk - A_max), the one-word report, restart and pair counters
  and walk index (zero; each switch sets the chunk counter): at most 32 instructions once, inside setup (section 5);
* the restart counter: incremented and compared with the cap 2^20 at each restart (`add`, `cmp`, `b.hs`,
  3 instructions), 2^20 * 3 / 20,000 below 2^7.30 units at the cap, a named term of the non-walk sum
  of section 5; the pair counter: the same 3 instructions per pair, inside its at most 600 units;
* the kernels' constant pools: set in setup (section 5);
* the stuck-walk counter and its mark: set to zero at each walk start (not once before the walks), and reset at
  every restart, report and fresh start (section 0-i);
* the multiplier table T[h][t] and the scalar tables a_h p^t, b_h p^t mod r (read only by the replays of section 0-ii): computed in setup, the table item
  of section 5 (32 * (1,021 + 11) units, 64 sampler draws and 2 * 384 products mod r at 800
  instructions each, and the 54 table-side forms of each of the 384 entries, 270 modular additions per entry at 60
  instructions; below 2^15.07 units);
* the prefix-product trie of section 3: built in setup level by level, each node from its parent by one `mul.s` call
  (the table entry's 12 coordinates against the parent's stored forms), then the 54 forms of its own entry and
  (internal nodes) its key entry (A..F, pre-sums, Ea, Ed), charged 1.6 units per node (31,985 instructions itemised in part D),
  2^42.929 + 2^51.510 nodes; its 2^62.855 bytes are written before
  they are read, and an operating system's zero-fill of its fresh 2 MiB pages is charged at 0.15625 instructions per
  byte (below 2^45.9 units); below 2^52.2098 units with the rest of setup;
* the sampler of section 3: the analysis is in the ideal model in which the bit stream is an ideal random source;
  the implementation realizes it with a ChaCha20 generator (RFC 8439's round function), assumed, as a disclosed
  heuristic of the same kind as H1 and H2, indistinguishable from the ideal mapping. No proved distributional or
  independence bound on the generated scalars is claimed: the key and nonce define a finite family of streams, and
  success is a heuristic over the ideal mapping, not a total-variation statement about that family. One attempt
  takes 255 bits of the stream (one 64-byte block) as X, accepts iff X < r by `cmp x0,x4; sbcs xzr,x1,x5;
  sbcs xzr,x2,x6; sbcs xzr,x3,x7; b.lo accept`, and otherwise retries (20 instructions with the retry counter and
  loop branch). A block is charged 2,300 instructions against an itemised 1,270 (80 quarter-rounds of 12
  instructions, 4 additions, 4 XORs and 4 rotates, 960; feed-forward 16; state load and store 32; the 64-bit block
  counter's `adds`/`adc` 2; addressing and loop control 10; spills, prologue and return 250), a 1.8x margin. The
  counter is 64 bits (ChaCha20's original 64-bit-counter/64-bit-nonce layout; RFC 8439's 96-bit nonce leaves a
  32-bit counter, which would wrap). A run draws at most 64 table scalars and two per walk start, stuck-walk restart
  and fresh start after a report, 64 + 2 (2^20 + 2^20 + 2^55) < 2^56.0001 samples, at most 2^7 attempts (blocks)
  each: below 2^63.0001 < 2^64 blocks, so the counter cannot wrap (hence R_MAX = 128, not the previous package's
  256). The charge is a fixed **2,320 instructions per attempt and R_MAX = 128 attempts, 296,960 instructions =
  14.848 units per sample**, whatever the draw. X is 255 bits and 2^254 < r < 2^255, so acceptance is above 1/2;
  after 128 rejects the sampler returns the fixed residue 1. Under the ideal-source model that has probability at
  most 2^-128 per sample and 2^-71.99 over a run (union bound), carried as an exhaustion allowance of 2^-71.99 in
  section 4. The 256-bit key, the 64-bit nonce (from the operating system's random source), the block counter and
  the retry counter are initialised in setup (2^12 instructions inside the 2^52.21 units). Nothing but the table
  and the starts samples. `experiment.py`'s `sample_below` runs the same compare, cap and rule (part S);
* each walk's state y and its trail's start pair (alpha0, beta0): a fresh start y = g^alpha0 z^beta0 from fresh uniform
  (alpha0, beta0) at every walk start, stuck-walk restart and after every report (at most 2^10.1 units each, section 5; two
  exponentiations of at most 1,021 units, two sampler draws and 43 instructions of bookkeeping take 1,050.7, inside the bound);
  the walk keeps (alpha0, beta0) unchanged for the record and carries no running exponents;
* each walk's histories, the two remembered K slots, the checkpoint K, the node (root), the depth and the step count
  mod 64 of section 3.3: set at every walk start, restart and fresh start, the K fields to the all-ones 2^384 - 1,
  which no stored K (below p < 2^381) equals, so the first steps cannot read a false repeat, and the count to zero;
  these few dozen instructions sit inside the fresh start's 2^10.1 units, of which the two exponentiations (at most
  1,021 units), the two draws and 43 instructions leave more than 46 units (asserted).

