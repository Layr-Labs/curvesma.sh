# Appendix: machine details (proof.md appendix.md sections A-i..A-iv)

This file carries the four machine-detail subsections of proof.md section 0 verbatim (renumbered A-i..A-iv), so proof.md stays under the verifier per-file byte cap. All cross-references to `section 0-x` in the package mean `appendix.md section A-x`.

## A-i Counters and the work meter: widths and charged updates

Every counter of section 3.4 has a stated width and a priced update. The stuck-walk counter must
reach 20/theta = 20 * 2^72, which needs more than one word, so it is a
**2-word (128-bit) counter** in a register pair, incremented with `adds x, x, #1; adc y, y, xzr`
so that it cannot wrap. The report count must reach the cap D = 2^55 only: it is **one 64-bit
word**, incremented with `add x, x, #1` and compared with the cap by `cmp x10, x11; b.hs` (3 instructions, charged as 5 inside the report item of section 5), with no carry. The work meter must hold a whole run in instructions,
2^125.238 * 20,000 < 2^139.58 instructions, so it is a **3-word (192-bit) counter** in three
registers, with a headroom factor of at least 2^52.42 (192 - 139.58 bits).

**Pre-debit rule.** The meter counts the instructions charged since the walk phase began. Three registers
hold a fixed 192-bit threshold L = B_walk - A_max, where B_walk is the walk-phase budget of section 3.4 and
A_max = 2^25 instructions is above the charge of any single walk-phase action (a step at most 17,837, an
escape at most 3 units, a report 4,000, a switch 600, a walk start or restart with its meter update and restart counter at
most 2^10.1 units plus 43 instructions; `experiment.py` asserts each, and also their sum, which bounds what a step
that reports or escapes adds between two compares). **Before** every charged action the
meter is compared with L; the action runs only if the meter is below L, and **after** it the meter adds
the action's charge. So the meter stays below L + A_max = B_walk. When the compare finds the meter at or above L, the walk phase
ends and nothing more is charged to the meter: that failed compare, the hand-over to the matcher and all
later work are paid from the stopping allowance of section 3.4.

Per step, both branches, the meter is charged **8 instructions** and the counter 5 in the table of section 5.
At the top of the step, before T is read, the meter is compared with L by the full three-word compare
`cmp x0, x5; sbcs xzr, x1, x6; sbcs xzr, x2, x7; b.hs end_walk_phase` (4 instructions), which branches exactly
when meter >= L as 192-bit unsigned integers: after a borrow chain the carry flag C is the full-width "no borrow"
result, so a carry-only condition (`b.hs`, `b.lo`) is exact, while Z reflects only the last word, so no condition that
reads Z (`b.hi`, `b.ls`, `b.eq`, `b.ne`) follows any chained compare here. `flags_check.py`, run by `experiment.py`,
executes this and every other multiword listing of 0-i, 0-ii and 0-iv as written under the AArch64 NZCV rules on edge and
random cases. **This exact compare runs before every step of every walk**, with no cheaper pre-test in front of it. Once the branch is known,
the step's charge (17,837 instructions for a multiply step, 7,452 for a squaring step: exactly the prices
the ledger of section 5 charges) is materialised in a register (`mov x4, #imm`, 1 instruction) and added
with `adds x0, x0, x4; adcs x1, x1, xzr; adc x2, x2, xzr` (3 instructions). A report, an escape and a walk
start or restart are separate charged actions under the same rule; each pays its own compare and add
inside its own charge (section 5). The 2-word stuck-walk
counter (low word x8, high word x9) costs 5 instructions per step:
`adds x8, x8, #1; adc x9, x9, xzr; mov x16, #5120; cmp x9, x16; b.hs restart`. Only the high word is compared
because the threshold is an exact multiple of 2^64: 20/theta = 20 * 2^72 = 5120 * 2^64, so the counter
has reached it exactly when x9 >= 5120. 5120 is not an immediate of `cmp` (12 bits, optionally shifted left by
12), so it is moved into the scratch register x16 on every step, a charged instruction.

The stuck-walk counter is **reset at every walk start, restart and distinguished-point report**, so it
counts the steps since the last report (or start) and the stuck-walk rule fires once per report gap, not
once per walk. The
reset at a report is two moves into the register pair (`mov x8, xzr; mov x9, xzr`, 2 instructions); it is charged 4 instructions inside the
4,000-instruction report item of section 5 (a listed line of the item; the total stays 4,000),
at all 2^55 reports. At a walk start or restart the reset is part of the 40 instructions of meter bookkeeping
below, and the fresh start itself is charged in section 5 (at most 2^10.1 units each, 2^46.1 for the 2^36
restarts of the cap).

**Schedule.** All M = 2^20 walks run on one processor, round-robin: it advances one walk for a chunk of
K = 2^63 consecutive step iterations, counted down in one register (`subs x12, x12, #1; b.eq switch`, 2 of the
itemised 19 of loop control of section 5), then a **switch** stores that walk's state (100 words: y 72, alpha/beta 8, T current/previous/checkpoint 18,
step count and history flags 2 — 50 stp of the 1 KB state bound), loads the next walk's (50 ldp), sets the
index mod M and its address (4) and the chunk counter (1), and pays its own compare and add (8) and 75 for
branch, call and margin: **188 instructions**, a charged action. A restart
replaces the walk's state in place and does not end the chunk; the phase opens with a switch. There are at most N_max/K + M + 2^36 + 1 < 2^63.71
switches (N_max of appendix.md section A-iii; an over-count), **below
2^56.93 units**, a term of the ledger (section 5). The meter is that processor's one 192-bit register triple for the whole
walk phase; no switch saves or loads it and nothing runs concurrently, so the pre-debit compare sees all work
charged before every action. Each start or restart pays 40 instructions (compare and add 8, counter reset 2, loads and moves 30):
2^20 * 40 < 2^11.04 units and 2^36 * 40 < 2^27.04 units for restarts, both in section 5's non-walk bucket.

## A-ii The report log and the matcher: one 64-bit address space, no initialisation, worst-case matching

**Report log.** Reports go into an append-only log: an array of exactly **REPORT_CAP = 2^55 record
slots of 80 bytes** with no pointer field, **2^55 * 80 = 2^61.322 bytes**, which is the whole memory of the
store. A record holds 16 bytes of fingerprint F, 32 bytes for the p^kappa-scaled alpha and 32 for beta,
each 32-byte field as four 64-bit words, most significant first; stored exponents are
residues in [0, r) (each of the report's two multiplications by p^kappa mod r ends with its conditional
subtraction, inside its 800 instructions). The slot of a report is the value of the one-word report counter
c before its increment, at byte offset 80 c = 64 c + 16 c (`add x2, x0, x1, lsl #6; add x2, x2, x1, lsl #4`).
Every byte offset is below 2^61.322 < 2^63, one native 64-bit register with the sign bit clear; there is no
sharding and no multiword address (`experiment.py` part D asserts each bound).
A report compares c with the cap; at c = 2^55 the walk phase ends and **the run fails**, a cap failure
inside the 0.01 of H1+ (section 3.4). There is no per-walk report cap: the global cap and its failure rule are
the only capacity rule, and the text needs no other.

**Nothing reads a slot that was not written.** The sort, the dedup and the scan below cover the slots
[0, c) only, every one of which a report wrote. So the log needs no clearing pass, no used bit and no free-slot
test. The walk states, counters and histories
of appendix.md section A-iv are still initialised, and charged, before they are read.

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

* *Heapsort of the whole 80-byte record* (F first, then alpha, then beta): Floyd's heap construction, then
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
  its predecessor and marks where a run of equal F begins.
* *Pairs*: inside each maximal run of records with equal F, every pair is taken in turn. A pair with
  beta1 = beta2 is two different classes (equal canonical element and equal beta force equal alpha, an identical
  record the dedup removed); otherwise the candidate x = (alpha1 - alpha2)/(beta2 - beta1) mod r is checked by
  g^x = z. A pair costs at most 600 units (counter, test, one inversion, the check). The matcher **halts at the
  first pair whose candidate verifies**. Two distinct records of one class always verify, so every pair charged
  before that is a false-key pair of two distinct classes, and the pairs charged are at most the false-key
  pairs plus one.
* *Hard pair budget*: a pair counter stops the matcher, and **the run fails**, when it would charge pair
  2^16 + 1. By linearity of expectation, with no independence assumption, the expected number of false-key
  pairs is at most (2^55)^2/2 = 2^109 record pairs times 2^-128 (the chance under H1 that two distinct classes share
  the 128 bits of F): 2^-19. Markov's inequality bounds the overflow probability by 2^-35; it is a cap failure
  inside the 0.01 of H1+.

The distinguished-point rate theta = 2^-72 (bits 16 to 87 of T, section 3) is what makes the log fit: the report cap scales with theta, so at theta = 2^-60 (the record's rate) it would be 2^55 * 2^12 =
2^67 reports occupying 80*2^67 bytes, which no 64-bit offset can address.

Memory (section 5): below 2^61.33 bytes, stated as peak memory below 2^63 bytes.
Memory is reported, not scored (section 6).

Only the stuck-walk counter (two words) and the work meter (three words) exceed 64 bits; they are
arithmetic counters, not addresses (appendix.md section A-i).

## A-iii One restart cap

The restart cap is **2^36 restarts**, quoted with these digits everywhere (sections 3.3, 3.4, 4 and 5,
README, note). Derivation. A walk is restarted
when its stuck-walk counter reaches 20/theta = 20*2^72 steps, and the counter is reset at every report
(appendix.md section A-i). The rule therefore fires per **gap**, a run of steps between two consecutive reports (or
between a start or restart and the first report), and each gap is a separate chance. A gap ends at the
first distinguished point; each step is one at density theta = 2^-72, so a gap exceeds
20/theta steps with probability (1 - theta)^(20/theta) <= e^-20 < 2.07e-9. The number of gaps is the number of
reports, expected theta times the steps. The walk phase runs until its budget is spent (section 3.4); the budget is
below 2^125.238 units and every step costs at least 0.3726 units, so the walk phase takes at most
N_max < **2^126.67** steps, and the expected number of stuck-walk restarts over the whole run is at most
N_max theta e^-20 < **2^25.81** (the experiment's figure is below 2^25.8085; each walk has about 2^35 gaps, so this is not a
per-walk figure). A walk caught in a fruitless cycle of more than 64 steps
with no distinguished class (rate below 2^-100 per step, section 3.3) is the other source: at most
N_max 2^-100 < **2^26.67** in expectation.
A long cycle that contains a distinguished class is never restarted, since each pass reports and resets the counter: its walk
appends one identical record per pass (at most 2^55 records, the cap) until the walk phase ends or the report cap fails the run,
and the meter bounds its work. Each such event is charged as a failed run. A cycle of l distinct classes contains a distinguished
one with chance at most l theta, so by the section 3.3 rates their expected number is
N_max theta (2^-288 sum_{l=65..147} l + sum_{l>=148} l (27/64)^l) < **2^-121.60** (experiment: 2^-121.607), which bounds
their probability; section 4 subtracts 2e-11.
The expected number of restarts of both kinds is below **2^27.30** (the experiment's 2^27.298), which is below 2^-8.65 of the cap 2^36.
No tail bound is taken: more than 2^36 restarts is a cap failure inside the 0.01 of H1+ (section 3.4), like every other cap. The
loss sum of section 4 charges the cap's own worst case as steps spent on restarts, not as their expected count: the
2^36 restarts' one-full-allowance worst case is 2^36 * 20/theta < 2^112.33 steps; part D computes the exact fraction < 4.83e-5 of
the walk's hard maximum (part D), subtracted as 4.83e-5; no other restart term is subtracted. The meter supports the cap: 2^36
restarts cost 2^36 * 40 < 2^41.33 instructions of meter bookkeeping (below 2^27.04 units, appendix.md section A-i), 2^36 fresh starts of at most 2^10.1
units each (2^46.1 units, section 5) and, at worst, one full stuck-walk allowance each, 2^36 * 20/theta < 2^112.33 steps,
which section 5 charges in the step bucket and whose nominal cost is below 2^111.58 units. The cap
therefore neither breaks the meter nor adds to the success allowance, and it is the only restart cap in this document.

## A-iv Initial values

Nothing the design reads is assumed to start at zero or at a known value; each is set and charged:

* the report log: not initialised, because nothing reads a slot before a report writes it (appendix.md section A-ii); an
  operating system's zero-fill of its fresh 2 MiB pages is charged, 13 instructions per record inside the 4,000;
* the 3-word work meter (set to zero), the 3-word threshold triple (loaded with L = B_walk - A_max), the
  one-word report, restart and pair counters and walk index (each set to zero; each switch sets the chunk counter): at
  most 32 instructions once, before the walks, inside the at most 2^31 units of the setup items (section 5);
* the restart counter: incremented and compared with the cap 2^36 at each restart (`add`, `cmp`, `b.hs`,
  3 instructions), 2^36 * 3 / 20,000 **below 2^23.30 units** at the cap, a named term of the non-walk sum
  of section 5; the pair counter: the same 3 instructions per pair, inside its at most 600 units;
* the stuck-walk counter: set to zero at each walk start (not once before the walks), and reset at every
  restart and report (appendix.md section A-i);
* the multiplier table T[h][t] and the scalar tables a_h p^t, b_h p^t mod r: computed in setup, the table item
  of section 5 (1,024 * (1,021 + 11) units, 2,048 sampler draws and 2 * 12,288 products mod r at 800
  instructions each, and the 12,288 x 54 Karatsuba pre-sums of the table entries, each a sum of at most 8
  reduced coordinates at 200 instructions, 6,635.52 units; below 2^20.2 units; `experiment.py` asserts it);
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
  only a 32-bit counter, which a run's below 2^45.01 blocks -- 2,048 table draws and two per start and
  restart, below 2^37.01 samples at at most 2^8 attempts each -- would wrap), so it cannot wrap and its
  update is inside the block's 2,300. It charges a fixed **2,320 instructions per attempt and R_MAX = 256
  attempts, 593,920 instructions = 29.696 units per sample**, independent of the draw. X is 255 bits and
  2^254 < r < 2^255, so acceptance is above 1/2; after 256 rejects the sampler returns the fixed residue
  **1** (the exhaustion rule). Under the ideal-source model the exhaustion event has probability at most
  **2^-256** per sample and at most 2^-218.99 over a run's below 2^37.01 samples (the union bound); it is
  carried as an exhaustion allowance of **2^-218.99** in section 4's budget and changes no digit. The
  256-bit key, the 64-bit nonce (the operating system's random source), the 64-bit block counter and the
  retry counter are initialised in setup, 2^12 instructions inside the 2^31-unit bound. The table draws
  two scalars per multiplier and each start and restart two; nothing else samples. `experiment.py`'s
  `sample_below` runs the same compare, cap and rule over its seeded bit source (part S);
* each walk's state y, alpha, beta: a fresh start g^alpha z^beta from fresh uniform (alpha, beta) at every
  walk start and restart (at most 2^10.1 units each, section 5; two sampler draws and 43 instructions of
  bookkeeping take 1,080.4, inside the bound);
* each walk's histories, the remembered previous T, the checkpoint T and the step count mod 64 of section 3.3:
  set at every walk start and restart, both T fields to the six-limb all-ones value 2^384 - 1, which no stored T
  (a residue below p < 2^381) equals, so the first steps cannot read a false repeat, and the count to zero.
  These few dozen instructions sit inside the fresh start's 2^10.1 units, whose two exponentiations take at
  most 1,021 units and leave more than 16 units (320,000 instructions; asserted).

