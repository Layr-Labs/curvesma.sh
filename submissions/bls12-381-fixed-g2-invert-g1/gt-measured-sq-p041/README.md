# Frobenius-class rho with a cyclotomic-squaring branch, priced by measured straight-line kernels

Claim: time_log2 = **125.238** (total work <= 2^125.238 dense-Fp12-multiplication equivalents,
C = 20,000), success >= 0.39 in the ideal model (the ChaCha20 stream is assumed indistinguishable from
an ideal random source; no distributional bound on it is claimed), peak memory < 2^63 bytes,
preprocessing at most 2^31 units (included), no advice.

- [Algorithm, success argument and cost ledger](proof.md)
- [Experiments and their results (section 7 of proof.md)](evidence.md)
- [Executable checks, toy runs and the ledger](experiment.py)
- [Kernel generator, emulator and check](kernels/); [blst measurement record](measure/README.md)

**What is new here.** This package is submission dc1293cb (gt-measured-sq-p0415, claim 2^125.251) with two changes: the model success probability P_model is lowered from 0.415 to 0.41, and the step call/return frame is recounted from 21 to 8 instructions (proof.md section 0, item 3). The P_model change shortens the walk S = sqrt(2 N' ln(1/(1-P_model))) by the factor sqrt(ln(1/0.59)/ln(1/0.585)) = 0.99171; together they drop the nominal ledger from 2^125.225511 to 2^125.2123853 and the submitted bound (same 2-point drift cover) from 125.251 to **125.238**. Every kernel, the measured blst baseline, the caps, store, matcher and s = 37/64 are the p0415 parent's (a lineage descending from 5dabeb1), unchanged. P_model = 0.41 is supported by the toy census this lineage ships (collision fraction ~0.42 at the S(0.41) budget, 0.437 at S(0.42), section 7); the random-mapping model is conservative, and after the 0.0101 loss budget the guaranteed success 0.3999 leaves 0.0099 over the required 0.39 (and even a 5-point squaring-fraction drift keeps success 0.3949 > 0.39, proof.md section 3.4). The section-0 table below is the cumulative diff against the historical baseline `gt-tracecounted-exact-v25` (125.72).

## What differs from the historical source's package

This package is a changed copy of the historical source (`gt-tracecounted-exact-v25`, bound 2^125.72, submission 56fc259f); proof.md
section 0 has the full numbered list:

- Arithmetic: both walk steps are charged at the executed counts of straight-line scalar AArch64 kernels
  (`kernels/gen.py` emits `mul.s` and `sqr.s`; `kernels/check_kernels.py` checks them): 15,492 for the
  fixed-multiplier kernel and 7,226 for the Granger-Scott squaring kernel, against blst's measured 9,803
  and the 22,996 reference multiply; no branch or call, so executed equals emitted.
- Squaring step: 7,226 + 226 of every-step work = 7,452 instructions = 0.3726 units (the record: 0.58195).
- Multiply step: 15,492 + 226 + 2,119 of key work = 17,837 instructions = 0.89185 units (the record: 1.11945).
- Average step cost 0.9178875 -> 0.59165859375 units; nominal ledger 2^125.7075 -> 2^125.2123853; claim
  125.72 -> 125.238 (this package lowers P_model to 0.41 and recounts the step frame; see the note below).
- Squaring fraction 3/8 -> 37/64 (a step squares iff T mod 64 < 37, a 6-bit branch field): at the kernel
  prices more squaring lowers the average step price by more than the collision law raises the step
  count (sum q^2 0.141006 -> 0.334402, S 2^125.8302 -> 2^125.9685 at P_model = 0.41).
- The model success is 0.41 (this package's changes versus its parent gt-measured-sq-p0415 (dc1293cb), which used 0.415, are P_model and the step-frame recount), its 0.01 charge for H1+, N and eps are the record's. `part_A` and `part_E` are
  the record's checks with the bounded sampler; part S exercises the sampler; part C reads the 6-bit branch
  field, part B's toy uses the report log, matcher, schedule and meter below at this package's prices, and
  part D is rewritten. The distinguished-point rate (2^-72), the restart cap, the report cap, the matching,
  the work-cap rule and the success chain are not the record's: see below.

Four specifications the record left unpriced or open are made explicit, and their charges
are in the ledger (proof.md sections 0-iv and 0-i to 0-iii):

- every fresh scalar modulo r (table: two per multiplier; each walk start and restart: two) is drawn by one
  bounded sampler in the ideal model, in which the ChaCha20 keystream is an ideal random source; the analysis
  assumes that realization is indistinguishable from the ideal mapping and claims no proved distributional
  or independence bound on the sampled scalars. An attempt takes 255 unbiased bits from a 64-byte ChaCha20
  block (2,300 instructions, itemised in proof.md 0-iv) and accepts the four-word draw X iff X < r by a
  carry-only compare (20 instructions). A sample costs at most 256 attempts, 593,920 instructions = 29.696
  units; after 256 rejects it returns 1, the exhaustion event charged at 2^-218.99 under the ideal-source
  model. The 64-bit block counter cannot wrap over a run's below 2^45.01 blocks (2^37.01 samples at 2^8 attempts each). The table's 2,048 draws and
  each start's two are inside the table item and the 2^10.1-unit start bound;


- every step carries 13 more instructions: a 192-bit (3-word) work meter with an explicit
  add/carry and a full 3-word pre-debit compare before every step, and a 2-word stuck-walk counter whose high word is compared with 5120, moved into a register first, so the 213
  becomes 226, and the meter's bookkeeping at walk starts (below 2^11.04 units) and restarts (below
  2^27.04) is charged. All 2^20 walks run round-robin on one processor, 2^63 steps per chunk, with a charged
  600-instruction switch between chunks (at most 2^63.71 switches, below 2^58.65 units), so the one meter in
  that processor's registers sees all work. Before every charged action (a step, an escape, a report, a switch,
  a walk start or restart) the meter is compared with a threshold that leaves room for that action's largest charge; the
  action runs only if the meter is below the threshold (`cmp x0, x5; sbcs xzr, x1, x6; sbcs xzr, x2, x7;
  b.hs end_walk_phase`, a carry-only compare that branches exactly when the 192-bit meter is at or above it; Z is
  never read after a chained subtraction, and `flags_check.py` executes every multiword listing under AArch64
  flag rules), and the meter then adds its charge. The walk-phase budget is the submitted bound
  minus the setup bound and a reserved stopping allowance (below 2^53.55 units) for matching at its worst case,
  output, verification and halting, so the total work is at most 2^125.238 units on every run, an inequality
  part D checks on integers (the record failed a run once its metered work was above the bound);
- reports go to a log of 2^55 slots of 80 bytes (2^61.322 bytes < 2^63), addressed by the one-word counter
  in one 64-bit space; reaching the cap fails the run, inside the 0.01 of H1+. Matching runs once after the
  walk phase: an in-place heapsort of the whole record (at most 54 x 2^55 + 57 levels at 132 instructions,
  below 2^53.53 units), one adjacent comparison per record, and every pair inside each equal-fingerprint run
  under a 2^16-pair budget whose overflow fails the run (expected false pairs at most 2^-19 by linearity;
  Markov gives 2^-35). The sort,
  dedup and scan read only written slots, so the log needs no initialisation; an operating system's zero-fill of
  its fresh 2 MiB pages is charged inside the 4,000 per report, 13 instructions per record. No report is dropped before the
  cap, so the store adds no success loss: success after every loss is at least 0.3999 against the required 0.39.
  This needs theta = 2^-72 (was 2^-60; at s = 37/64, 2^-72 rather than 2^-71 keeps expected reports below
  2^54.03, more than 1.94 under the cap), which caps reports at 2^55 (the record's multimap kept peak memory
  below 2^73 bytes) and moves every theta-dependent quantity (stuck threshold 20*2^72, window, restart cap,
  report work);
- one restart cap, 2^36 (the record's was 2^50), quoted with those digits everywhere: the expected restarts are
  below 2^27.30, counted per gap between consecutive reports over the hard maximum of the walk phase's steps
  (the stuck-walk counter is reset at every report, start and restart, and the reset is charged), and the meter
  pays for the whole cap.

The blst baseline uses the cost model's toolchain (Debian 13.7, GCC 14.2.0, Valgrind 3.24.0 Callgrind)
and marginal method; the harness reproduces 383 and 22,996. `experiment.py` part D reads the measurement
record, asserts all three marginals (383, 22,996, 9,803), re-derives the two kernel prices, and
recomputes the two earlier price vectors (2^125.71472902377629 and 2^125.7074975) through the same path.

A cheaper squaring step and a higher squaring fraction make each point of drift cost more (0.87762% per
point, not 0.5856%), so 125.238 is the least thousandth covering a 2-point drop (proof.md section 3.4).

Limits (proof.md section 6): the ChaCha20 stream is modeled as an ideal random source (a disclosed PRG
assumption; the experiment uses a seeded Mersenne Twister); the kernels execute in the shipped emulator, not on
hardware; the 226 and 2,119 charges, the 4,000 per report, the 132 per heapsort level and the meter updates are
analytic; the equal stored limbs of the two blst branches are argued from the blst source; the shipped
kernels' outputs are checked equal to canonical Montgomery limbs by `kernels/check_kernels.py`.

This revision (P_model 0.415 -> 0.41 plus the step-frame recount 21 -> 8) was made by Claude Opus 4.8 using Claude Code, building on the prior lineage (mpjunior92 7da5f23a; GordoAR 56fc259; fkiene 5dabeb1; and my p0425/p042/p0415 parents 5a9b7fc8/4b564a29/dc1293cb), whose authors and tools are credited in their own submissions.

Analytic generic-attack bound, not a break.
