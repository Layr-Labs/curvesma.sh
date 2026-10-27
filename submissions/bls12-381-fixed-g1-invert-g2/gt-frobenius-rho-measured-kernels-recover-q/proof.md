# Frobenius-class rho (2^20 round-robin walks) with a cyclotomic-squaring branch, priced by its measured instruction count

**Claim.** For uniform nonidentity P in G1 and uniform z in GT the algorithm below
returns Q in G2 with e(P,Q) = z for that exact P with probability at least 0.39, using at most
**2^125.014** dense-Fp12-multiplication equivalents (C = 20,000), including setup, the two
pairings, table construction, walk starts, distinguished-point (DP) handling, cycle
escapes, the DP detection delay, matching (the sort, duplicate removal and pair checks), output and verification. Peak memory is below 2^63
bytes, preprocessing (included) at most 2^31 units (the setup items of section 5), nonuniform advice 0 bytes. The success
bound rests on the usual random-mapping heuristic for rho walks, adjusted for the branch
mix (section 4), and is supported by toy-field experiments and by running the production
branch rule (with a 16-entry table and no distinguished-point test) in the real GT. The success probability is computed in the ideal model, in which the ChaCha20
stream that supplies the table and start scalars is indistinguishable from an ideal random source; that is
a disclosed heuristic, and no total-variation or independence bound on the realized scalars is claimed.
No production-size run is claimed.

The class construction (Frobenius classes, the 4-multiplication key, the twisted multiplier
table) is published in submission 7da5f23a-6633-46b1-9c03-11f6cfedeeec (mpjunior92). The
cyclotomic-squaring branch, the branch-mix collision law, the report and memory caps, and the
branch, multiplier and DP bits read from T = Re(c_0) (so a squaring step needs no key) are
later revisions of that design, the historical source of section 0, which lists every difference from it.

## 0. What this revision changes

This package is a changed copy of the historical source's package (`gt-tracecounted-exact-v25`, claim
2^125.72, submission 56fc259f), which solves the mirror problem on the fixed-G2 track. Here 'record' means
that historical source only, not a baseline of this track; changed quantities are this package's.

0. *Target.* This package solves this track: P in G1 is the supplied nonidentity input, Q0 the fixed G2
   generator, g = e(P, Q0) (pinned order pairing(Q0, P)) and the output is Q = [x]Q0 in G2 (section 1). The GT
   walk is the same map on GT whatever the base. The G2 output is charged at its counted hard maximum,
   12,809 units, not the record's 510 G1 operations; both pairings at 2^14 units each (counted at 10,331).

Prices (items 1 to 3 and 9 below):

| Work | record's charge | this package |
| --- | ---: | ---: |
| every step: 270 itemised (section 5) + 13 meter and counter (section 0-i) | 270 = 0.0135 units | **283 = 0.01415 units** |
| multiply step (keys 2,119 itemised, unchanged) | 1.11945 units | **15,137 = 0.75685 units** |
| cyclotomic squaring routine | 11,369 (analytic) | **6,786** (kernel, section 5) |
| squaring step in all | 0.58195 units | **7,069 = 0.35345 units** |
| squaring fraction s (branch rule) | 3/8 (T mod 8 < 3) | **37/64** (T mod 64 < 37, a 6-bit field) |

**Changes from the record** (record value -> this package, where the record states one):

1. *Squaring price.* 0.58195 -> **0.35345** units: the kernel `kernels/sqr.s` (6,786 executed; blst's
   measured 9,803 in `measure/` is the baseline) replaces the analytic 11,369.
2. *Multiply-step price.* 1.11945 -> **0.75685** units: the multiply kernel's 12,735 executed instructions
   replace the reference 20,000, plus 2,119 of keys and 283 of every-step work.
3. *Every-step work.* 270 -> 283: the 192-bit meter's exact 3-word pre-debit compare and add/carry update
   (8) and the 2-word stuck-walk counter with its compare (5), widths and updates in section 0-i, plus 40
   per walk start or restart.
4. *Average step cost, ledger and claim.* c 0.9178875 -> **0.523634375** units; W_nominal 2^125.7075 ->
   below **2^125.0125964**; claim **125.72 -> 125.014**, the least thousandth at or above the ledger with a
   0.05-point drift design margin (drift moves success only), checked on integers (section 3.4).
5. *Distinguished-point rate.* theta 2^-60 -> **2^-72** (DP bits 16 to 87 of T, fingerprint bits 88 to 215), so the report cap is D = 2^55, the stuck-walk threshold
   20/theta = 20*2^72 and the restart cap **2^36** (the record's 2^50; section 0-iii). The ledger's step
   bucket outside the walk is c (M + 2^36) 20/theta (the record's c M 20/theta) and the non-walk bucket
   2^54.2 units (the record's 2^64.3).
6. *Report store and matching.* The record's chained multimap (load factor 2/3, below 2^73 bytes, 6,000
   instructions per report, matched on arrival) is replaced by a log of exactly 2^55 slots of 80 bytes (peak
   memory below 2^63 bytes) at 4,000 instructions per report, matched once after the walk phase: an in-place
   heapsort, adjacent dedup and a hard budget of 2^16 pairs whose overflow fails the run (section 0-ii).
7. *Success chain.* The store loses nothing before its cap, whose failure is inside the 0.01 of H1+. The
   model value is 0.4001 (record 0.43) and the losses sum to at most 0.0100503, below 0.010051, so success
   after every loss is at least 0.390049, 0.000049 over 0.39 (section 4); eps = 3.2e-4 (record 6.5e-4, section 3.3).
7a. *Work-cap rule.* The record failed a run once its metered work passed the bound. Here a pre-debit meter
   runs the walk phase until its budget, net of setup and a reserved stopping allowance, is spent, so total work is at most W_cap on every run; the M walks run round-robin, K = 2^63 steps a
   chunk (sections 0-i and 3.4).
8. *Code.* Section 7 (evidence.md) runs these rules and prices; part D recomputes every ledger figure.
9. *Squaring fraction.* s = 3/8 -> **37/64** (T mod 64 < 37; the record's T mod 8 < 3). At these prices more squaring lowers the average step price
    (c = 0.523634375) by more than the collision law raises the step count: sum q^2 = 0.334402 (record
    0.141006), log2 N' at most 251.859405 (record 251.4914), S at most 2^125.9454 (record 2^125.8102). At theta = 2^-71
    the expected reports would be at most 0.9642 of the cap, at 2^-72 more than 2.07 below. Section 7's tables are
    at s = 37/64.

### 0-i Counters and the work meter: widths and charged updates

Every counter of section 3.4 has a stated width and a priced update. The stuck-walk counter must
reach 20/theta = 20 * 2^72, which needs more than one word, so it is a
**2-word (128-bit) counter** in a register pair, incremented with `adds x, x, #1; adc y, y, xzr`
so that it cannot wrap. The report count must reach the cap D = 2^55 only: it is **one 64-bit
word**, incremented with `add x, x, #1` and compared with the cap by `cmp x10, x11; b.hs` (3 instructions, charged as 5 inside the report item of section 5), with no carry. The work meter must hold a whole run in instructions,
2^125.014 * 20,000 < 2^139.31 instructions, so it is a **3-word (192-bit) counter** in three
registers, with a headroom factor of at least 2^52.69 (192 - 139.31 bits).

**Pre-debit rule.** The meter counts the instructions charged since the walk phase began. Three registers
hold a fixed 192-bit threshold L = B_walk - A_max, where B_walk is the walk-phase budget of section 3.4 and
A_max = 2^25 instructions is above the charge of any single walk-phase action (a step at most 15,137, an
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
the step's charge (15,137 instructions for a multiply step, 7,069 for a squaring step: exactly the prices
the ledger of section 5 charges) is materialised in a register (`mov x4, #imm`, 1 instruction) and added
with `adds x0, x0, x4; adcs x1, x1, xzr; adc x2, x2, xzr` (3 instructions). A report, an escape and a walk
start or restart are separate charged actions under the same rule; each pays its own compare and add
inside its own charge (section 5). The 2-word stuck-walk
counter (low word x8, high word x9) costs 5 instructions per step:
`adds x8, x8, #1; adc x9, x9, xzr; mov x16, #5120; cmp x9, x16; b.hs restart`. Only the high word is compared
because the threshold is an exact multiple of 2^64: 20/theta = 20 * 2^72 = 5120 * 2^64, so the counter
has reached it exactly when x9 >= 5120. 5120 is not an immediate of `cmp` (12 bits, optionally shifted left by
12), so it is moved into the scratch register x16 on every step, a charged instruction. The kernels restore x19-x30 and sp and overwrite x0-x17, so the step
routine saves the meter, threshold and counters around each kernel call, inside its 21 frame and 28
load-and-move instructions (section 5).

The stuck-walk counter is **reset at every walk start, restart and distinguished-point report**, so it
counts the steps since the last report (or start) and the stuck-walk rule fires once per report gap, not
once per walk. The
reset at a report is two moves into the register pair (`mov x8, xzr; mov x9, xzr`, 2 instructions); it is charged 4 instructions inside the
4,000-instruction report item of section 5 (a listed line of the item; the total stays 4,000),
at all 2^55 reports. At a walk start or restart the reset is part of the 40 instructions of meter bookkeeping
below; the fresh start itself is charged in section 5.

**Schedule.** All M = 2^20 walks run on one processor, round-robin: it advances one walk for a chunk of
K = 2^63 consecutive step iterations, counted down in one register (`subs x12, x12, #1; b.eq switch`, 2 of the
itemised 30 of loop control of section 5), then a **switch** stores that walk's state (at most 2 KB, 256 stores), loads the
next walk's (256 loads), sets the index mod M and its address (4) and the chunk counter (1), and pays its own
compare and add (8) and 75 for branch, call and margin: **600 instructions**, a charged action. A restart
replaces the walk's state in place and does not end the chunk; the phase opens with a switch. There are at most N_max/K + M + 2^36 + 1 < 2^63.64
switches (N_max of 0-iii; an over-count), **below
2^58.58 units**, a term of the ledger (section 5). The meter is that processor's one 192-bit register triple for the whole
walk phase; no switch saves or loads it and nothing runs concurrently, so the pre-debit compare sees all work
charged before every action. Each start or restart pays 40 instructions (compare and add 8, counter reset 2, loads and moves 30),
charged in section 5's non-walk bucket.

### 0-ii The report log and the matcher: one 64-bit address space, no initialisation, worst-case matching

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
the only capacity rule.

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
arithmetic counters, not addresses (section 0-i).

### 0-iii One restart cap

The restart cap is **2^36 restarts**, quoted with these digits everywhere (sections 3.3, 3.4, 4 and 5). Derivation. A walk is restarted
when its stuck-walk counter reaches 20/theta = 20*2^72 steps, and the counter is reset at every report
(section 0-i). The rule therefore fires per **gap**, a run of steps between two consecutive reports (or
between a start or restart and the first report), and each gap is a separate chance. A gap ends at the
first distinguished point; each step is one at density theta = 2^-72, so a gap exceeds
20/theta steps with probability (1 - theta)^(20/theta) <= e^-20 < 2.07e-9. The number of gaps is the number of
reports, expected theta times the steps. The walk phase runs until its budget is spent (section 3.4); the budget is
below 2^125.014 units and every step costs at least 0.35345 units, so the walk phase takes at most
N_max < **2^126.64** steps, and the expected number of stuck-walk restarts over the whole run is at most
N_max theta e^-20 < **2^25.78** (the experiment's 2^25.7796, rounded up; each walk has about 2^35 gaps, so this is not a
per-walk figure). A walk caught in a fruitless cycle of more than 64 steps
with no distinguished class (rate below 2^-100 per step, section 3.3) is the other source: at most
N_max 2^-100 < **2^26.64** in expectation.
A long cycle that contains a distinguished class is never restarted, since each pass reports and resets the counter: its walk
appends one identical record per pass (at most 2^55 records, the cap) until the walk phase ends or the report cap fails the run,
and the meter bounds its work. Each such event is charged as a failed run. A cycle of l distinct classes contains a distinguished
one with chance at most l theta, so by the section 3.3 rates their expected number is
N_max theta (2^-288 sum_{l=65..147} l + sum_{l>=148} l (27/64)^l) < **2^-121.63** (experiment: 2^-121.6361), which bounds
their probability; section 4 subtracts 2e-11.
The expected number of restarts of both kinds is below **2^27.27** (the experiment's 2^27.269), which is below 2^-8.73 of the cap 2^36.
No tail bound is taken: more than 2^36 restarts is a cap failure inside the 0.01 of H1+ (section 3.4), like every other cap. The
loss sum of section 4 charges the cap's own worst case (2^36 full stuck-walk allowances) as steps spent on restarts, not
their expected count. The cap neither breaks the meter nor adds to the success
allowance, and it is the only restart cap in this document.

### 0-iv Initial values

Nothing the design reads is assumed to start at zero or at a known value; each is set and charged:

* the report log: not initialised, because nothing reads a slot before a report writes it (section 0-ii); an
  operating system's zero-fill of its fresh 2 MiB pages is charged, 13 instructions per record inside the 4,000;
* the 3-word work meter (set to zero), the 3-word threshold triple (loaded with L = B_walk - A_max), the
  one-word report, restart and pair counters and walk index (each set to zero; each switch sets the chunk counter): at
  most 32 instructions once, before the walks, inside the at most 2^31 units of the setup items (section 5);
* the restart counter: incremented and compared with the cap 2^36 at each restart (`add`, `cmp`, `b.hs`,
  3 instructions), 2^36 * 3 / 20,000 **below 2^23.30 units** at the cap, a named term of the non-walk sum
  of section 5; the pair counter: the same 3 instructions per pair, inside its at most 600 units;
* the stuck-walk counter: set to zero at each walk start (not once before the walks), and reset at every
  restart and report (section 0-i);
* the kernels' constant pools: set in setup (section 5);
* the multiplier table T[h][t] and the scalar tables a_h p^t, b_h p^t mod r: computed in setup, the table item
  of section 5 (1,024 * (1,021 + 11) units, 2,048 sampler draws, 2 * 12,288 products mod r at 800
  instructions each and 12,288 x 270 m-side additions at 60 each, below 2^20.2 units;
  `experiment.py` asserts it);
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
  original 64-bit-counter/64-bit-nonce layout; RFC 8439's 32-bit counter would wrap over a run's below
  2^45.01 blocks: 2,048 table draws and two per start and restart, below 2^37.01 samples of at most 2^8
  attempts), so it cannot wrap and its
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

## 0a. Representation: which limbs the rule reads

A walk element is held as the six-limb Montgomery vectors that the pinned blst kernels store,
x~ = x*2^384 mod p. The decision rule is stated on those stored limbs, which is what a limb-level
implementation of the design reads at zero cost (reading an integer residue instead would need a
Montgomery reduction per T read and twelve per multiply step, work no charge here contains):

* T is the six limbs of Re(c_0) **as stored**. The branch bits are (stored T) mod 64, the table index is
  (stored T >> 6) mod R, and the distinguished-point predicate is read from the same word.
* kappa is the minimum of the twelve **stored** key products a*zeta^j and d*zeta^j, and the twist index
  is the smallest index attaining it.

This is the rule on ordinary integers composed with a fixed bijection of Z/p, so every structural
property carries over: the twelve stored keys are distinct, their multiset is the same for all twelve
conjugates (Frobenius permutes it), so kappa is class-determined and the twist index is covariant;
stored T is invariant across conjugates because Re(c_0) is. The selected multiplier index is a
pseudorandom function of the class under either rule, which is all H1 uses. No conversion is charged.

**The squaring branch returns the same stored form as the multiply branch.** The walk runs the two
straight-line kernels `kernels/mul.s` and `kernels/sqr.s`, not the blst routines. `kernels/check.py` compares
every output limb of both kernels with the canonical value x*2^384 mod p of a big-integer reference (0
mismatches in `kernels/check_kernels.py`), so both return the unique stored limb vector of their result and an
element reached by either branch is stored identically. The shipped Python walk computes in ordinary integers
and applies the Montgomery factor only to the values it reads (part C, stored rule), which gives the same bits.

Evidence, part C, both rules, 2,000 production-field steps each: the stored rule
gives 0 class-check failures in the first 150 steps and the exponent invariant y = g^alpha z^beta holds
at the end; its squaring fraction is 0.5824 over 1,997 branch-rule decisions (2,000 steps less 3
two-cycle escapes), 0.42 points above 37/64, z-score 0.38 (section 3.4), chi-square of the low 6 bits
59.56 (63 dof) and of the next 4 bits 16.57 (15 dof); it visits 1,997 distinct T values, the last new one
at step 1,999, and shows no repeating period at the end. The integer rule gives 0 failures and the
exponent invariant, with squaring fraction 0.5877 over 1,996 decisions, z-score 0.86, chi-squares 57.26
(63 dof) and 15.43 (15 dof), 1,996 distinct T values, the last new one at step 1,999, and no repeating
period. The selected key index kappa differs from the integer rule in 755 of the stored rule's 834
multiply steps, so the two rules are different walks and the claim is about the one the charges can
pay for.
`experiment.py` reports both blocks (`partC_production_walk`, `partC_production_walk_integer_rule`).

## 1. Reduction

The input is a nonidentity P in G1 (order r) and z in GT. Compute g = e(P, Q0) with the pinned convention
pairing(Q0, P). The reduced pairing is bilinear and nondegenerate and Q0 generates G2, so for every permitted
P, g != 1 and g generates GT (prime order r); as P = [b]P0 with b uniform in [1, r-1], g = e(P0, Q0)^b is
uniform over GT minus 1, independent of z. z = g^x for a unique x in [0, r). If z = 1 output Q = infinity
(x = 0). Otherwise find x by the walk below, output Q = [x]Q0 in G2, and check e(P, Q) = z; bilinearity gives
e(P, [x]Q0) = g^x = z, for the supplied P. `reduction.py` (part F of `experiment.py`; its checks are listed in evidence.md) runs this on
BLS12-381: the pinned pairing operation by operation (10,331 counted operations, each charged one unit;
charged 2^14) and Q = [x]Q0 by double-and-add with the pinned projective formulas over Fp2 and one Fp2
inversion. The ladder's hard
maximum takes the most expensive branch of the pinned addition (28 operations, the equal-points branch) and
doubling (20) at each of at most 254 steps plus the inversion (617): **12,809 units**, every Fp2 operation
charged a full dense Fp12 multiplication, an overcount.

## 2. Arithmetic facts (checked in `experiment.py`, parts A and C)

Tower: Fp2 = Fp[u]/(u^2+1), Fp12 = Fp2[w]/(w^6 - xi), xi = 1+u (the blst tower). With
u = w^6 - 1 this is the profile field, w^12 - 2w^6 + 2 = 0, and conversion of coefficients
is free: c_m = a_m - b_m, c_{m+6} = b_m.

**F1 Frobenius classes.** (sum c_m w^m)^p = sum conj(c_m) d_m w^m with d_m = xi^{m(p-1)/6}.
p has order exactly 12 modulo r, p^6 = -1 and 1 + p^4 + p^8 = 0 (mod r), so the maps
y -> y^{p^k} form a group of order 12 acting on GT, containing inversion. Every y != 1
lies in a class [y] of exactly 12 elements; N = (r-1)/12 classes, log2 N = 251.2721.

**F2 Keys.** Put zeta = xi^{(p^2-1)/6}; it lies in Fp and is a primitive 6th root of unity
(zeta^2 = zeta - 1). Let c_1 = a + b u be the w^1 coefficient of y and B = Re(conj(c_1) d_1)
= a Re(d_1) + b Im(d_1). The w^1 coefficient of y^{p^(k+1)} is conj(c^{(k)}) d_1, and
conj(d_1) d_1 = zeta, hence
`key_k := Re(w^1 coeff of y^{p^k}) = a zeta^j (k = 2j), B zeta^j (k = 2j+1)`.
With zeta^j in {1, zeta, zeta-1, -1, -zeta, 1-zeta} all 12 keys cost **4 Fp multiplications**
(a zeta, a Re d_1, b Im d_1, B zeta) plus 8 subtractions/negations and one addition. The
12 keys are distinct unless two of them coincide by a linear accident (probability <= 66/p
per element; ties are charged as failure, below 2^-247 over a run, section 4). K([y]) = min_k key_k is a class invariant,
kappa(y) = argmin_k key_k, and y^{p^kappa} is the same element for all 12 members of the
class (the canonical representative y^). K is the minimum of 12 values that behave like
independent uniform residues (they come in negation pairs, so K lies in [0, p/2]); its low 200
bits are uniform up to statistical distance about 2^200/p = 2^-181 (the density of the minimum is smooth).

**F2b Trace invariant.** The w^0 coefficient of y^{p^k} is conj^k(c_0) (no twist), so T([y]) := Re(c_0) is the same
for all 12 conjugates, and 12 Re(c_0) is the trace of y to Fp. T is read directly from the coefficients (no
multiplication) and, for distinct classes, is an essentially uniform residue (collisions of T between distinct classes
have probability about 1/p per pair). It plays the role of the class identity for the bits below; K is needed only to
find kappa. `experiment.py` part A checks that Re(c_0) is equal over all 12 conjugates of random GT elements.

**F3 Cyclotomic squaring (Granger-Scott).** GT lies in the cyclotomic subgroup, where
squaring is cheaper than a product. Write y = A + B w + C w^2 with A=(c0,c3), B=(c1,c4),
C=(c2,c5) in Fp4 = Fp2[v]/(v^2-xi), v = w^3, and (A0,A1) = A^2, (B0,B1) = B^2,
(C0,C1) = C^2 computed in Fp4 with 3 Fp2 squarings each. Then
`c0' = 3A0 - 2c0, c3' = 3A1 + 2c3, c1' = 3 xi C1 + 2c1, c4' = 3C0 - 2c4, c2' = 3B0 - 2c2, c5' = 3B1 + 2c5.`
`experiment.py` verifies y^2 = formula on real BLS12-381 GT elements and counts the
operations of the implementation: **9 Fp2 squarings and 34 Fp2 additions/subtractions**
(the multiplications by xi are counted as one Fp2 addition).

blst's `blst_fp12_cyclotomic_sqr` (`src/fp12_tower.c` lines 630-660) is this map: three `sqr_fp4` calls
and the 3t - 2u / 3t + 2u combinations, on blst's Fp12 = Fp6[w]/(w^2 - v) layout, in which our coefficient
c_m of w^m is the Fp2 element a[m mod 2][m div 2]. The walk computes it with the kernel `kernels/sqr.s`, which
section 5 prices; blst's measured count is the baseline.

**F4 Orders.** r - 1 = 2^32 * 3 * 11 * 19 * 10177 * 125527 * 859267 * 906349^2 * 2508409 *
2529403 * 52437899 * 254760293^2 (verified by multiplication in `experiment.py`). The order
of 2 modulo r is (r-1)/16 (about 2^250.9), so 2^t lies in the order-12 group <p> only if t is a
multiple of about 2^247.3.

## 3. Algorithm

Parameters: R = 1024 table multipliers, squaring probability s = 37/64, theta = 2^-72 (DP
density), M = 2^20 walks round-robin on one processor (section 0-i). Trails do not depend on the order walks are
advanced, so parallel collision search's analysis (M walks, distinguished points) and section 4 apply: the schedule
changes which walk's steps occur first, not the step budget or any trail, and a collision of two walks is found at
matching when both trails reached the distinguished point within the walk phase. Bits of the class invariant T = Re(c_0) (read as the six stored Montgomery limbs of section 0a):
T mod 64 selects
the branch (squaring iff T mod 64 < 37), h = (T >> 6) mod 1024 selects the multiplier, the class
is distinguished iff (T >> 16) mod 2^72 = 0 (bits 16 to 87 of T). The fingerprint of a class used for matching reports is
F = (T >> 88) mod 2^128 (bits 88 to 215 of T), the 128 bits above the branch, multiplier and DP bits, modeled as independent under H2;
T has 381 bits, so every bit used is a bit of the stored T and none is shared between the roles.

Precompute M_h = g^{a_h} z^{b_h} for h < R (fresh uniform a_h, b_h mod r, drawn by the sampler of
section 0-iv) and
T[h][t] = M_h^{p^t}, t = 0..11, is the T table (1,024 x 12 = 12,288 entries of 576 bytes, 7.1 MB);
the scalar tables a_h p^t, b_h p^t mod r hold 2 x 12,288 = 24,576 entries of 32 bytes (0.8 MB in all).

A walk keeps some member y of the current class together with (alpha, beta) mod r such that
y = g^alpha z^beta. Each walk starts at fresh uniform (alpha, beta), two sampler draws of section 0-iv. One step:

```
1  T <- Re(c_0) of y (no multiplication); branch, h and DP bits from T
2  if class distinguished: compute the 12 keys (4 Fp mult) -> kappa; append (F, p^kappa*alpha, p^kappa*beta)
   to the report log of section 0-ii
3  2-cycle / checkpoint tests on T (section 3.3)
4  if T mod 64 < 37 (squaring branch):  y <- y^2 (Granger-Scott);  alpha,beta <- 2*alpha, 2*beta
   [no keys, no twist]
   else (multiply branch): compute the 12 keys (4 Fp mult) -> kappa; h <- (T>>6) mod 1024;
        t <- (12 - kappa) mod 12;  y <- y * T[h][t]   (ONE dense Fp12 multiplication);
        alpha += a_h p^t;  beta += b_h p^t
```

**3.1 The step is a map on classes.** Write y^ = y^{p^kappa}. Multiply branch:
y * M_h^{p^{-kappa}} = (y^ * M_h)^{p^{-kappa}} lies in the class [y^ M_h]. Squaring branch:
y^2 = (y^^2)^{p^{-kappa}} lies in [y^^2], because squaring commutes with Frobenius. Both
depend only on the class (T is a class function; kappa and y^ enter only through the conjugation
absorbed by the twist), so walks that reach the same class stay together and reach the same DP.

**3.2 Collision.** The log keeps every report until the cap (records are appended, never replaced or
dropped), so a useful collision is lost to the store only by a cap failure. After the walk phase the
matcher of section 0-ii sorts the log and takes the pairs inside each run of equal F. Two records with equal F
are the same class, up to
a false match of probability 2^-128 per pair of distinct classes: with D <= 2^55 reports there are at most
D^2/2^129 = 2^-19 false-key pairs in expectation; every pair's candidate is verified (section 0-ii's 600-unit charge). Their canonical exponents (alpha^, beta^) = p^kappa (alpha, beta) satisfy
alpha^1 + x beta^1 = alpha^2 + x beta^2 mod r, so x = (alpha^1 - alpha^2)/(beta^2 - beta^1).
beta^1 = beta^2 happens with probability 1/r. Output Q = [x]Q0 in G2 after checking g^x = z.

**3.3 Fruitless cycles.** A closed walk of the class graph is *useless* if the canonical
exponent pair returns unchanged. After l steps the pair is Lambda (alpha^_0, beta^_0) + tau
where Lambda = 2^t p^m (t squaring steps, m a sum of twists) and tau = sum_h c_h (a_h, b_h)
with c_h in Z[p]. Because (alpha^_0, beta^_0) is uniform and independent of the table, the pair
returns unchanged with probability about 1/r unless Lambda = 1 and every c_h = 0.
* **No squaring step can occur in a useless cycle.** Lambda = 1 with t >= 1 needs
  2^t in <p> (F4), i.e. t >= 2^247.
* A useless cycle made of multiply steps needs every table index used to occur at least twice
  with coefficients p^{e_1} + p^{e_2} + ... = 0 in Z[p]/(Phi_12(p)): for two occurrences
  e_2 = e_1 + 6 (probability 1/12), for three {e, e+4, e+8}, and so on; for any number of occurrences
  the probability is at most 1/12 (the last twist must hit one value). For 2-cycles this gives the
  model rate (27/64)^2/(12 R) = 2^-16.07 per step at the production R = 1024; at the
  census's R = 64 the prediction is 2.32e-4 and the measured rates are 3.27e-4 (q = 67, 3,000 runs) and
  2.21e-4 (q = 163, 900 runs), with 2.65e-4 in the declared 40-trial run. For l = 2k indices drawn independently and
  uniformly (heuristic H2) the probability is at most (2k)!/(2^k k!) (12R)^-k (pairings dominate), which
  is below 2^-288 for 65 <= l <= 147. For l >= 148 all l consecutive steps must be multiply steps, which
  has probability (27/64)^l < 2^-184 (the branch bits are bits of T, H2; the charged allowance keeps the
  conservative 2^-100 rate). Such cycles are never met in
  practice; if one is and holds no distinguished class, the stuck-walk rule restarts the walk after 20/theta = 20*2^72 steps, which adds
  less than 2^-13 relative work even at the restart cap of section 0-iii (2^-13.62, section 5); one that holds a distinguished
  class is not restarted, a failed run (section 0-iii).
* **Detection.** The step before last is remembered; if its T equals the current T the two
  classes form a 2-cycle and the walk replaces the representative of the class with the smaller
  T by its square (y -> y^2, exponents doubled), cost <= 3 units. This depends on the cycle
  only, so merged walks stay merged. A checkpoint T is stored every 64 steps and compared at
  every step; a repeat means a cycle of length <= 64, which is traversed once more to its minimum-T
  class that is then squared (waste <= 200 steps). Every other walk whose stuck-walk
  counter reaches 20/theta = 20*2^72 steps is restarted from a fresh uniform pair (stuck-walk rule). The counter is reset to zero
  at every walk start, restart and report (section 0-i), so a walk is restarted only after 20/theta steps with no report, and the
  rule is checked once per report gap.
* Total overhead factor charged on the step count: eps = 1/(12R) + 250/R^2 < 3.2e-4 (the record charged 5/(12R) + 250/R^2). Derivation for this walk: a 2-cycle wastes at most two multiply steps and one escape, 2 * 0.75685 + 3 = 4.5137 units, 8.62 average steps (c = 0.523634375); a cycle of 3 to 64 steps needs every index used twice, rate below 10 (12R)^-2 (the l = 4 pairing term is 3 (12R)^-2), at most 250 steps each. At the model 2-cycle rate the overhead is below 1.42e-4; at the R = 64 census's largest measured rate 3.27e-4 scaled by 64/R it is below 1.93e-4, and eps exceeds that by a factor above 1.65 (`experiment.py` part D asserts both).

**3.4 Caps and the work cap.** The run has three phases: setup (at most 2^31 units, section 5), the walk
phase, and the stopping phase (matching, output, verification, halting). The stopping phase is bounded in
advance: the heapsort of at most 2^55 records, the dedup, at most 2^16 pairs at 600 units, [x]Q0 and the
verification pairing, and 100 instructions for the failed compare, the hand-over and halting, together the
**stopping allowance R_stop, below 2^53.55 units** (`experiment.py` part D). The walk phase has the budget
**B_walk = W_cap C - 2^31 C - R_stop C instructions**, with W_cap = 2^125.014 the submitted bound and
W_cap C rounded down to an integer. The walk phase meters its own charged work with the prices of section 5
(0.75685 units per multiply step, 0.35345 per squaring step, the stated costs of escapes, switches, starts,
restarts and reports) in the 3-word (192-bit) meter of section 0-i under its **pre-debit rule**: before each
charged action the meter is compared with L = B_walk - A_max, the action runs only if the meter is below L,
and the meter then adds the action's charge, at most A_max. So the walk phase charges below B_walk, and
when a compare finds the meter at or above L the walk phase ends and the stopping phase begins. Every run therefore does at most
2^31 C + B_walk + R_stop C <= W_cap C instructions, a deterministic inequality that `experiment.py` part D
checks on integers. The walk phase can always spend L (asserted); the probability that the threshold binds early is part
of the failure budget below.

*Heuristic H1+ (no drift of the walk statistics).* The cost bound does not rest on H1+: the pre-debit meter
enforces W_cap on every run, whatever squaring fraction the run realises, so drift moves success only. A
multiply step costs 0.75685 units and a squaring step 0.35345 (difference 0.4034), so a drop of D points of
the squaring fraction raises the average step price by a factor 1 + delta, delta = 0.004034 D /
0.523634375, at most 0.77039% per point. The cap carries a drift design margin of **0.05 points** (delta at
most 0.038520%), the rounding floor: the nominal ledger with every step priced that much higher is below
2^125.0131520, the least thousandth at or above it is 125.014, and the exact integer inequality
(numerator^1000 <= denominator^1000 * 2^125014, that ledger rounded up to 40 digits) is asserted in part D.
The cap exceeds W_nominal < 2^125.0125964 by at least 0.0973% (0.0014 bits); the walk budget's per-step slack
exceeds 0.000509 units, so the threshold does not end the walk phase before the ledger's walk is done while
the squaring fraction averaged over the whole run stays within **0.1263 points** of 37/64 (asserted). That
average is over n = S(1+eps) > 2^125.94 steps; under H2 it is a mean of n near-independent Bernoulli(37/64)
decisions, standard deviation sqrt(s(1-s)/n) < 2^-63.99, so 0.1263 points is more than 2^54.3 standard
deviations (both asserted). Short samples test H2, not this average: part C (stored rule) reads 0.5824 over
1,997 decisions (z-score 0.38); part B's toy walks at most 0.66 points below 37/64 (evidence.md). Drift beyond the tolerance is charged to success: the walk phase ends
early, the budget buys n/(1+delta) steps of the nominal n, 1/(1+eps) of them useful, under the section-4
collision law (rounded down; loss up):

| drop of the squaring fraction | success (lower bound) | loss from 0.4001 (upper bound) |
| ---: | ---: | ---: |
| up to 0.1263 points (the cap's tolerance; design margin 0.05) | 0.4001 | 0 |
| 0.5 points | 0.3983 | 0.0018 |
| 1 point | 0.3960 | 0.0041 |
| 1.5 points | 0.3936 | 0.0065 |
| 2 points | 0.3913 | 0.0088 |
| 3 points | 0.3868 | 0.0133 |
| 5 points | 0.3779 | 0.0222 |

H1+ states that over the whole run the branch, 2-cycle and DP frequencies stay inside a drift of this size
(the 0.01 budget is not exhausted until about 2.28 points, more than 2^58 whole-run standard deviations;
asserted) and the reports, false candidates and restarts
stay within the caps below, with total probability of any cap failure or early threshold at most 0.01. This is a much weaker
statement than independence of the steps; it excludes only a drift by percentage points over 2^125 steps.
Support: (i) H2 models the low six bits as uniform, giving nominal squaring probability 37/64; exact uniformity over production classes is not proved;
(ii) in the toy runs the fluctuation of the mean charge scales as the independent-steps prediction:
rms of (relative deviation x sqrt(steps)) is 0.4395 (q = 67 census) and 0.4199 (q = 163) against 0.4313 predicted at the earlier step prices (the
declared 40-trial run reads 0.3439 against 0.3805 at this package's), so there is no excess correlation of branch labels at that scale; (iii) the part-C and toy
fractions above agree. The toy runs cannot test
concentration at 2^125 steps; H1+ extrapolates it, exactly as H1 extrapolates the collision law. The
success analysis (section 4) charges 0.01 for it.

*Caps.* reports D <= 2^55, counted in a one-word counter (section 0-i); the expected number of reports is theta
times the steps the walk budget buys, below 2^53.95 at the nominal average price, more than a factor 2.07 below the cap, and
reaching the cap fails the run; at most 2^16 charged pairs in the matcher (section 0-ii: overflow
probability at most 2^-35, a failed run; the pair work is at most 2^16 * 600 units < 2^25.24); **at most 2^36 restarts** (section 0-iii; their starts and bookkeeping, 2^27.04 units, and their
stuck-walk allowances are charged in section 5); the total number of steps
that the ledger needs is S(1+eps) + (M + 2^36) 20/theta. The walk phase runs until its budget is spent,
so the detection window and restart allowance are steps the budget buys,
and the steps beyond the ledger's come from the cap's slack over the nominal ledger, at least
2^115.00 units (0.0973%). Fruitless cycles
longer than 64 steps occur at a rate below 2^-100 per step (section 3.3; below 2^26.64
events in expectation; one with no distinguished class wastes at most 20/theta = 20*2^72 steps before the stuck-walk restart, below 2^102.96 steps
in total in expectation; one through a distinguished class is a failed run, section 0-iii).
These caps bound the work and memory stated below.

## 4. Success probability

**H1 (random mapping with branch mix).** The class sequence behaves like iterates of a
random mapping on N classes in which two distinct classes have the same image with
probability (1 - sum_b q_b^2)/N, where q_b are the branch probabilities. Justification: for
c1 != c2 the equation F(c1) = F(c2) has, for each branch pair (b1,b2) with b1 != b2, about one
solution among the N classes (12 candidate representatives, each canonical with probability
1/12); for b1 = b2 = squaring it has none (squaring is a bijection of GT of odd order that
commutes with Frobenius), and for b1 = b2 = a table index it has about 11/12 of one. We
use the conservative value sum q^2 = s^2 + (1-s)^2/R = 0.334402 for all pairs; the exact
value (s^2 + (1-s)^2/(12R)) is smaller by 1.6e-4. Part E of `experiment.py` tests this law directly (section 7). This is the same correction as the Brent-Pollard
factor R/(R-1) of the pure adding walk (s = 0), and for Pollard's original three-branch walk
it is the mechanism that makes Pollard's original square-and-add walk slower than an r-adding walk
(Teske 1998). Hence the effective class count is
N' = N / (1 - sum q^2), log2 N' at most 251.859405.

Under H1 the probability that the first S steps contain a class collision is
1 - exp(-S^2/(2N')). Choose P_model = 0.4001:
`S = sqrt(2 N' ln(1/0.5999)) < 2^125.9454.`
Losses before the first useful collision is detected: selected-trail detection failure
20 e^-20/(1 - e^-20) < 4.13e-8 for each of the two selected trails, the colliding and the collided-into one,
8.25e-8 in all (length-biased; the fresh-trail e^-20 = 2.07e-9 does not apply), a key tie (66/p per key
computation over at most N_max of them, below 2^-247), a
collision inside a single walk (fraction about 1/M = 1e-6), beta^1 = beta^2 (1/r, below 2^-254.8), cap
failures and the threshold ending the walk phase early (at most 0.01 together, H1+, section 3.4; the report cap and
the pair budget are among these caps), steps spent on restarts (the restart cap's own worst case: 2^36
restarts at one full stuck-walk allowance of 20/theta = 20*2^72 steps each are 2^36 * 20/theta < 2^112.33
steps; part D computes the exact fraction < 4.92e-5 of the walk's hard maximum, subtracted as 4.92e-5;
restart counts and costs, section 0-iii), a long cycle through a
distinguished class, a failed run, probability below 2^-121.63, subtracted as 2e-11, and the sampler's
exhaustion allowance, at most 2^-218.99 under the ideal-source model of section 0-iv, subtracted as
2^-218.99. The report log loses no report before the cap (section 0-ii), so the store adds
no loss term. The sum is the 0.01 of cap failures plus the restart and small terms: at most 0.0100503, **below 0.010051**
(experiment.py part D). The success after every loss is therefore at least 0.4001 - 0.010051 = 0.390049,
leaving **0.000049 absolute slack** over the required 0.39 for the inaccuracy of H1. The value 0.4001 (the record's 0.43) is set against the pre-registered census400 of evidence.md section 7
(240,000 toy runs at q = 67, 163 and 211; q = 379 was not in it): at every tested q the fraction of runs with a useful
collision within S(0.40), minus 4 binomial standard errors, must be at least 0.391, a rule fixed before any count was
read. S(0.4001) > S(0.40), so that fraction is a lower bound on the fraction within S(0.4001). In the completed census p - 4 SE is at least 0.3932 at every
tested q (toy support, not a production proof). The earlier census405 (S(0.405), a longer budget) is provenance only.
Individual settings of a few hundred runs scatter by +-0.04.

## 5. Cost (units of one dense Fp12 multiplication = 20,000 instructions)

Counts are executed scalar AArch64 instructions as the cost model counts them
(`tracks/docs/COST_MODEL.md`): Valgrind Callgrind, Debian 13.7, GCC 14.2.0, blst f7a88067 built `-O3
-fno-builtin -fPIC -fno-tree-vectorize -mgeneral-regs-only`, collection around the call loop, marginal
(Ir(65)-Ir(1))/64; `measure/` holds the record and part D reads it.

| routine | per-call marginal | status |
| --- | ---: | --- |
| `blst_fp_mul` | 383 | measured; equals the cost model's reference count |
| `blst_fp12_mul` | 22,996 | measured; equals the cost model's reference count |
| `blst_fp12_cyclotomic_sqr` | 9,803 | measured, on real GT elements, identical at five slice offsets |

| walk kernel | executed instructions | status |
| --- | ---: | --- |
| `kernels/mul.s` (x * fixed entry) | 12,735 | executed count, branch-free, 0 mismatches in `kernels/check_kernels.py` |
| `kernels/sqr.s` (Granger-Scott) | 6,786 | executed count, branch-free, 0 mismatches in `kernels/check_kernels.py` |

**Table multiply:** the product y*T[h][t] uses `kernels/mul.s`, **12,735 executed instructions** (0.63675
units), below blst's dense `blst_fp12_mul` reference operation (22,996).

**Squaring branch:** the Granger-Scott map of F3 is what blst's `blst_fp12_cyclotomic_sqr` (`src/fp12_tower.c` lines 630-660) computes. The straight-line kernel `kernels/sqr.s` computes it in
**6,786 executed instructions** (0.3393 units), charged in full; blst's measured 9,803 for the same map is
the baseline it replaces. The squaring update needs no keys; DP report keys are charged separately.

**Kernels.** `kernels/gen.py` emits both: general-purpose registers only, no NEON, no branch, no call
(`ldp`/`stp`/`extr` allowed), Montgomery form, 6 x 64-bit limbs, outputs fully reduced. The multiply kernel uses
that m is fixed: every m-side linear form is stored reduced mod p, so the kernel adds only x-side operands. On the
layout of section 2, c0 = k1 + k2 and c1 = k1 + k3 with k1 = (x0 + x1) m0, k2 = x1 (v m1 - m0), k3 = x0 (m1 - m0);
each Fp6 product is six Fp2 products and each Fp2 product three Fp products (`gen.py` `mul_plan`), 54 in all, 12
per output coordinate. `gen.py` asserts each output equals the reference model mod p, every x-side operand is an
exact sum below 2^384, and each output's 768-bit sum is in [0, 2^768) with stored values as independent inputs in
[0, p) (a multiple K of p 2^384 from the 0..7p table is added to the real parts); the squaring kernel
triples a sum as Z + 2Z with an `extr` doubling chain. No branch and no call, so executed equals emitted on
every call. `kernels/emu.py`
interprets that subset. `gen.py` asserts r < 8p for every reduced r and floor(r/p) - 1 <= q <= floor(r/p) for
q = umulh(r5, 9), r5 the top limb, at both ends of each [kp, (k+1)p) (exhaustive: q is monotone in r); the kernel
subtracts q p from a table of 0..7p, then p unless that borrows, so every output is in [0, p).
`kernels/check_kernels.py` regenerates both listings and checks a sample, asserting 12,735 and 6,786 executed;
`experiment.py` runs a smaller sample. The 54 stored values per entry (2,592 bytes; 270 additions at 60 instructions, above
the 45 priced below) are built in the table construction, inside its 2^20.2-unit item. The listings are regenerated, not shipped.

**Per-step work.** Every step (both branches), analytic:

| Work | Instructions |
| --- | ---: |
| read Re(c_0) (6 limbs), branch bits, DP test, table address | 50 |
| (alpha, beta) update mod r (two 4-limb add + conditional subtract, or doubling) | 100 |
| 2-cycle and checkpoint comparisons on T, history | 60 |
| loop control: back-branch (1), chunk counter subs/b.eq (2), branch select on the 6-bit field, and #63/cmp #37/b.lo (3), DP test branch (1), 2-cycle and checkpoint test branches (4), call and return of the step routine with frame save and restore (21) | 32 |
| margin for loads and moves | 28 |
| **subtotal: the record's 270, with the loop-control row itemised** | **270** |
| 192-bit work meter: step charge in a register (1), adds/adcs/adc (3), exact 3-word compare with the cap, cmp/sbcs/sbcs/b.hs, on every step (4) (section 0-i) | 8 |
| 2-word stuck-walk counter: adds/adc increment (2), mov of 5120 into x16 (1), cmp of the high word with x16 and b.hs (2), exact as 20/theta = 5120 * 2^64 (section 0-i) | 5 |
| **total counted 283, charged** | **283 = 0.01415 units** |

Multiply steps in addition (the keys give kappa), analytic, built on the measured blst_fp_mul = 383 and
a modular add/sub/neg of 6 limbs called as a routine at 45 instructions (6 add/sub, 6 conditional-subtract,
6 select, 12 loads, 6 stores, 9 for call, return and address set-up):

| Work | Instructions |
| --- | ---: |
| 4 Fp multiplications (a zeta, a Re d_1, b Im d_1, B zeta) | 1,532 |
| 1 Fp addition (B) + 8 Fp subtractions/negations, 45 each | 405 |
| minimum of the 12 keys (12 top-limb loads, 11 x 3 compare/select, tie test, index) | 150 |
| load of the w^1 coefficient, kappa, twist index | 32 |
| **total counted 2,119, charged** | **2,119 = 0.10595 units** |

So a multiply step costs (12,735 + 283 + 2,119) / 20,000 = **0.75685 units** and a squaring step
(6,786 + 283) / 20,000 = **0.35345** units.
Average step cost c = (27/64) * 0.75685 + (37/64) * 0.35345 = **0.523634375 units**.

Setup (at most 2^31 units, the preprocessing bound): pairing g = e(P, Q0) (2^14, section 1), check z^r = 1 (<= 510), the
table (1024 multipliers x (2 exponentiations <= 1,021 units + 11 Frobenius maps) + 2,048 sampler draws
x 29.696 units, the scalar tables and the 12,288 x 54 stored m-side forms, <= 2^20.2), the sampler's state (<= 2^12 instructions), M = 2^20
starts g^alpha z^beta (<= 2^30.1, each with two sampler draws), the initial counter values (at most
32 instructions, section 0-iv) and the kernels' constant pools (2 x 194 words: 0..7p, n0, 9, at most twelve 12-word K,
1,940 instructions); their sum is below 2^30.11 units. DP handling: at most 2^55 reports (the cap of section 3.4) x **4,000 instructions**,
itemised as: 12 keys for kappa 2,119; two 255-bit modular multiplications by p^kappa mod r 1,600; the pre-debit
compare of the meter with the threshold and the charge's 3-word add 8; the one-word report counter compare
against the cap, branch and increment 5; the slot address base + 80 c 2; appending the 80-byte record (ten 8-byte
stores and address set-up) 12; reset of the 2-word stuck-walk counter 4; operating-system zero-fill of the
log's 2 MiB pages on first write (131,072 stores and one fault per 26,214 records) 13; margin for loads, moves,
call and return 237. That is 2^55 * 4,000 / 20,000 < 2^52.679 units. Meter bookkeeping (section 0-i) adds 2^20 * 40 / 20,000 < 2^11.04 units for the walk starts and 2^36 * 40 / 20,000 < 2^27.04 for the
restarts, the restart counter adds below 2^23.30 units (section 0-iv), and the restarted walks' fresh starts add
2^36 * 2^10.1 = 2^46.1 units (each start's two sampler draws are inside its 2^10.1 bound). Matching, output and halting, the stopping allowance of section 3.4, at their worst
case: the heapsort of 2^55 records (7,188 instructions per record, section 0-ii, below 2^53.53 units), the dedup (2^55 x 100 instructions), at most
2^16 pairs at 600 units (below 2^25.24 units), [x]Q0 in G2 (at most 12,809 units, section 1), the verification
pairing (2^14) and 100 instructions of hand-over and halting, together below 2^53.55 units. All non-walk work
other than the switches (section 0-i, charged separately below) is therefore below
2^53.55 + 2^52.679 + 2^46.1 + 2^31 + 2^27.04 + 2^23.30 + 2^11.04 < 2^54.2 units (`experiment.py` part D
sums the exact terms, below 2^54.181, and asserts the sum below 2^54.2).

Detection delay: after the first useful collision the ledger charges 20/theta further
steps for each of the M walks, E = M * 20 / theta < 2^96.33 steps. With one stuck-walk allowance per
restart under the cap of section 0-iii, the step bucket outside the walk is (2^20 + 2^36) * 20/theta < 2^112.33 steps
(2^-13.62 of S). As 20/theta = 10,240 K, the colliding walk gets its 20/theta further steps while each
other walk gets at most 10,240 chunks, 20/theta, so the window holds; switches add below 2^58.58 units (section 0-i).

`W = c * S * (1 + eps) + c * (2^20 + 2^36) * 20/theta + 2^54.2 + 2^58.58`, evaluated by `experiment.py` part D with the exact S = 2^125.945389..., 1 + eps = 1.000319798... and (2^20 + 2^36) * 20/theta = 2^112.321...: **W < 2^125.0125964**.

Submitted bound: **time_log2 = 125.014** (the least thousandth at or above the ledger recomputed with the
average step price raised by the 0.05-point drift factor of section 3.4, below 2^125.0131520), at least 0.0014 bits above the
nominal ledger and below the cap headroom of section 3.4. `experiment.py` (part D)
recomputes this ledger, and two earlier price vectors (11,400/300/2,200 and the record's 11,369/270/2,119:
2^125.7147 and 2^125.7075, at their own s = 3/8, model 0.43, theta = 2^-60 and non-walk input 2^63.8) through the same code path, with exact decimal arithmetic, and checks the 125.014 inequality on integers.

Memory: the report log of section 0-ii is 2^55 slots x 80 bytes (128-bit fingerprint, two 256-bit
exponents, no pointer field) < 2^61.322 bytes, sorted in place with one 80-byte temporary, plus the T table
(12,288 entries of 576 bytes, 7.1 MB), its 54 stored m-side forms (2,592 bytes per entry, 31.9 MB), the
scalar tables (24,576 entries of 32 bytes, 0.8 MB) and 2^20 walk states (at most 2 KB each, including the
remembered previous state, 2^31 bytes): below 2^61.33 bytes, and the claim states below 2^63 bytes. Memory
is reported, not scored.

## 6. Heuristics and limitations

* H1 above (random-mapping behaviour of the mixed walk on classes), H2 (pseudorandom branch and DP bits
  of T), independence of table indices (used for cycle probabilities). The success probability is computed
  in the ideal model: a uniformly random mapping on the classes (H1), realized by table multipliers and walk
  starts drawn from the ChaCha20 generator; the analysis assumes that realization is indistinguishable from
  the ideal mapping, discloses this as a heuristic of the same class as H1 and H2, and claims no proved
  distributional or independence bound on the sampled scalars. Under the ideal-source model the bounded
  rejection is exact uniform conditioned on acceptance, and the exhaustion rule is the event charged at
  2^-218.99 in section 4. The experiment draws the same bounded reduction from Python's seeded Mersenne
  Twister, so its part S exercises the cap, compare and exhaustion rule, not the ChaCha20 block. These are
  standard in rho analyses (Teske 1998/2001; van Oorschot-Wiener 1999; Duursma-Gaudry-Morain 1999;
  Bos-Kleinjung-Lenstra 2010), supported here by toy and production-field runs, not proved.
* Scope of the prices. The two walk kernels' counts are executed counts from the shipped emulator, equal to
  the emitted counts because the programs have no branch or call; `kernels/check_kernels.py` checks a
  sample. The other charges are analytic bounds: the 283 of every-step work, the 2,119 of key work, the
  4,000 per report and the 132 per heapsort level. The blst prices (383, 22,996, 9,803) are measured with
  the cost model's toolchain and are the baseline; the walk itself was not compiled and run as a blst
  program; the equal stored limbs of the two branches rest on the kernel check (section 0a).
* Collision statistics are validated in toy fields only, at R = 64 (section 4, evidence.md); the declared
  40-trial run is a single small sample and no agreement is claimed for it.
* Memory is unscored, and the claim is an instruction-count bound: no physical machine maps 2^61.322 bytes.
* This is a generic square-root attack. It does not use index-calculus (TNFS) methods in Fp12.
  It is a valid upper bound, not a security estimate.

## 7. Experiments (`experiment.py`, standard library only)

The experiments of parts A to E and S and their results, the offline toy aggregates and the pair-collision
tables, are in `evidence.md` (this section, moved whole to keep this file under the size limit).

References: Pollard 1978; van Oorschot & Wiener, J. Cryptology 12 (1999); Teske, Math. Comp.
70 (2001); Duursma, Gaudry & Morain, ASIACRYPT 1999; Granger & Scott, PKC 2010; Aranha et
al., EUROCRYPT 2011; Bos, Kleinjung & Lenstra, ANTS 2010; mpjunior92, submission 7da5f23a.
