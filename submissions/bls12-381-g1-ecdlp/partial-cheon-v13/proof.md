# bls12-381-g1-ecdlp: partial-coverage Cheon attack v13, at n=2^20

Status: ready for independent exploratory review. Classical single-shot
attack using the supplied auxiliary power A_d, d = 2^20; no quantum
steps, no uncharged advice, no amortization over targets. Cheon's
discrete-logarithm-with-auxiliary-input attack is standard cryptanalysis.

v13 responds to the tenth-round review, which accepted the reduction,
the parity-key injectivity, the coverage argument, the fixed-parameter
factorization, and the toy result, and rejected (1) a scratch bound left
inconsistent across files (v12's claim text still said 2^32 bytes where
the README itemized 304 B/record) and (2) a gcd instruction price
asserted without an implementation. This document is regenerated whole
so that every figure is the current one, stated once, consistently.
Every instruction allowance below is either an executed construction
(the binary gcd, the merge sort, the binary search, the comb) or carries
full limb-level line items; every resource figure is computed
programmatically; the consistency of all bounds across the five files of
this submission is enforced by a grep gate before submission.

## Claimed bounds

| Field | Value |
| --- | --- |
| time_log2 (score) = log2(W(2^20)) | **123.836032** (= 1934938/15625, rounded up) |
| success_probability (supported lower bound) | 0.39 (exact per-offset minimum (m-1)m/N >= 0.4 - 2^-118; no probabilistic quantity anywhere) |
| memory_log2_bytes | 124.8 (global peak = the sort phase, 256m = 2^124.767581 bytes) |
| preprocessing_log2 (included in total) | 21.0 (preprocessing total <= 2^20.01 units) |
| nonuniform_advice_log2_bytes | 0 |

W(n) for n >= 2^20 is constant: the attack reads only A_1 = T and
A_d = [s^d]P0 with d = 2^20 <= n (at the scoring n = 2^20, A_d is the
last element, index d-1). The bound does not apply for n < 2^20.

## Conventions

Protected target: BLS12-381 G1, generator P0 of prime order r,
challenge s uniform in [1, r-1), T = A_1 = [s]P0, powers
A_i = [s^i mod r]P0 for i = 1..n, scored at n = 2^20. Cost model
`classical-group-work-v2`: cost = H + R/2000; H = reference operations
(batched incomplete affine G1 addition, 6M, S = M, A = 0; one reference
operation = 1 unit), R = other scalar 64-bit AArch64 instruction
equivalents. Pinned base-field multiplication: 383 instructions (blst
`blst_fp_mul`, mod p). Mod-r multiplication (r 255 bits): 5x52-limb
schoolbook (25 limb products of 52x52 -> 104 bits via 64x64->128
instructions, plus carry accumulation) followed by Montgomery
reduction of the same shape: charged **512 instructions** (0.256
units). A reference operation and its implementation are never
double-counted.

## Number-theoretic setup (verified by experiment.py Part A)

- The complete prime factorization of r - 1 is verified in the
  experiment (R_M1_FACTORS, hex literals; every factor Miller-Rabin
  prime, product equality asserted exactly):
  r - 1 = 2^32 * 3 * 11 * 19 * 10177 * 125527 * 859267 * 906349^2 *
  2508409 * 2529403 * 52437899 * 254760293^2; every prime factor is
  < 2^28.
- 2^20 | r-1, and d = 2^20 is the largest divisor of r-1 not exceeding
  the scoring n (any divisor with an odd part is at most 2^18 * 3 < 2^20).
- zeta = 7 generates Z_r*; eta = zeta^(2^20) has exact order
  N = (r-1)/2^20 (log2 N = 234.857); zeta_e = zeta^N has exact order
  2^20.

## Why the algebra works

s uniform in [1, r-1) lies in Z_r*. The d-th power map on Z_r* (order
r-1, d | r-1) is exactly d-to-1 onto the unique order-N subgroup, so
s^d is uniform there and s^d = eta^k1 for unique k1 uniform in [0, N);
and s = zeta^(k1 + t*N) for unique t in [0, 2^20). Hence
A_d = [eta^k1]P0, Q' = [zeta^(-k1)]A_1 = [zeta_e^t]P0.

## Algorithm v13 (single shot; the only coin is the window offset)

0. Read A_1 and A_d (charged).
1. Runtime factorization of r-1 (the bounded construction below);
   construct zeta, eta, zeta_e; build fixed-base comb window tables
   for P0, A_d, A_1 (8,160 points each, <= 2^13 additions per table,
   ~3 MiB each; <= 2^15.579 units total).
2. Phase 1, partial-coverage BSGS with a contiguous window:
   m1 = ceil(sqrt(N)), M = ceil(N/m1) (= m1 here), and
   m = isqrt(floor(2N/5)) + 1 = 141430898198098033869460368287841487
   (pure integer arithmetic; m < m1).
   - Baby scalars: eta^i computed sequentially (one 512-instruction
     mod-r multiplication per step), never stored as an m-sized list:
     the streaming comb consumes each block of scalars as it is
     produced.
   - Baby points: fixed-base comb multiplication (the construction
     below), one 64-unit charge per point, streamed block by block
     into the m-entry record array (the record array is the table by
     design and is the only m-sized object of this phase).
   - Baby table: sort the m records (key = (x, y mod 2), payload = the
     baby index) by bottom-up merge sort with an auxiliary buffer:
     117 passes, at most m * 117 key compares (15 instructions each:
     six limb loads + compares + branch), at most m * 118 record moves
     (50 instructions each: 16 loads + 16 stores + 4 index updates +
     control; includes one possible final copy-back), and <= 16
     instructions per record of one-time key preparation: <= 7671
     instructions per baby, an upper bound asserted in the experiment.
     The auxiliary merge buffer is freed after the sort.
   - Window: uniform offset o in [0, M); the m consecutive blocks
     o..o+m-1 (mod M) form at most two contiguous runs (run 2 starts
     at block 0 with scalar 1); each giant scalar is one sequential
     mod-r multiplication; at most two window-start exponentiations
     (charged).
   - Giants, STREAMING: giant scalars are generated and consumed in
     blocks of 2^24 records (never materialized as an m-sized list);
     each block's accumulator array is allocated for exactly the
     block size; each block's points are binary-searched immediately
     and then discarded. The sorted baby table stays resident
     throughout - the only m-sized live object of this phase.
   - Lookups: binary search of the sorted array for the key
     (x(Y), y(Y) mod 2): at most floor(log2 m) + 1 = 117 probes (the
     candidate range at most halves per probe; enforced per lookup in
     the experiment), each <= 32 instructions (2-word midpoint 6,
     pointer + record loads 4, six-limb key load/compare 12, parity
     bit 2, interval update + termination + branches 8): <= 3744
     instructions = 1.872 units per lookup. On a hit at baby index i:
     k1 = (i + m1*j) mod N - exact, because the key match is an exact
     point match (the parity key is unique per point; see below).
3. Phase 2 (deterministic, full coverage): Q' = [zeta^(-k1)]A_1; BSGS
   for t in [0, 2^20) with 2^11 comb multiplications and a 2^10-entry
   sorted baby table searched identically (all charged).
4. s = zeta^(k1 + t*N) mod r; verify [s]P0 = T; output. No phase-1 hit
   -> fail; no restart or amplification.

### The fixed-base comb construction (executed by the toy)

Per multiplication: one window-table hit assignment plus at most
ceil(bits/winbits) - 1 = 32 window additions for a 255-bit scalar at
8-bit windows, each window addition exactly one reference
batched-incomplete-affine operation by the cost model's definition
(inversion amortization included in the reference), with every window
layer's additions within a block sharing ONE Montgomery batched
inversion; window addressing and control <= 304 instructions per
multiplication. Construction total ~32.15 units; charged 64 units
(2x). Live working storage per block, itemized (all listed objects
coexist; no in-place aliasing is claimed or needed): block scalars
32 B + affine accumulators 96 B + Montgomery prefix products 48 B +
inverse-output array 48 B + pending-operation entries (index, two
point references, denominator) 72 B + giant indices 8 B = 304 B per
record, charged at <= 512 B per record x 2^24 = <= 2^33 bytes. The
toy EXECUTES this construction: per-block live accumulator and scalar
array lengths are recorded AND asserted (live_accs = live_scalars = 64
at toy scale, 5 blocks), and equality with plain double-and-add is
asserted in full for EVERY point - all m babies and every giant of
every trial.

### The sort key: (x, y mod 2)

On this curve y^2 = x^3 + 4 over odd Fp, the (at most two) points
sharing an x have parities y and p - y, which differ, so (x, y mod 2)
uniquely determines the affine point. Baby and giant records are both
raw affine points (no canonicalization anywhere in the Cheon attack);
a key match is therefore an exact point match and k1 = (i + m1*j) mod N
is exact with no sign ambiguity. The uniqueness property is proved in
one line, asserted on the baby set (m distinct keys), and sampled on
negation pairs in Part A.

### Runtime factorization: bounded, deterministic, executed

The input r-1 is the fixed public BLS12-381 parameter and the seed
sequence is fixed, so the computation is DETERMINISTIC: the hosted run
reproduces the measured 111,717 modular multiplications, 310 gcds and
9 rho cycles (Part A reports the counts). No probabilistic model is
used anywhere in this claim. The hard caps bound the worst case for
adversarial inputs:

- Each rho cycle checks its remaining iteration budget BEFORE each
  consumed batch (halving rounds, 128-step mini-batches, and the
  backtrack scan), so a COMPLETED cycle provably used at most the
  2^17-iteration cap; a hard 12-cycle guard bounds the total.
- Every modular multiplication costs 512 instructions (the explicit
  5x52 schoolbook + Montgomery construction); the inner loop performs
  two per iteration plus 32 instructions of loop, index and batching
  work: 1056 instructions per main-loop iteration.
- Every gcd is an EXECUTED round-counted binary gcd (Stein's
  algorithm, implemented in the experiment in place of math.gcd,
  verified against math.gcd on 201 sampled pairs including
  (r-1, 254760293^2)): each round is a 4-limb compare (<= 10
  instructions), a conditional subtraction (<= 10), multiword
  trailing-zero elimination (<= 12), a multiword right shift (<= 10)
  and loop control (<= 6): <= 48 instructions per round, and each round
  strips >= 1 bit of a < 2^256 total, so <= 512 rounds per call
  (asserted per call; measured maximum 153): <= 24,576 instructions
  per gcd.
- The backtrack scan is bounded by the 128-step mini-batch size per
  cycle; each backtrack step is one 512-instruction multiplication
  plus one 24,576-instruction binary gcd = 25,088 instructions.

Joint worst case, three bounded terms:

| Term | Bound | Instructions |
| --- | --- | ---: |
| main loop | 12 cycles x 2^17 iterations x 1056 | 1,660,944,384 = 2^30.629357 |
| batched gcds | 12 x 2^10 gcds x 24,576 | 301,989,888 = 2^28.169925 |
| backtrack | 12 x 128 steps x 25,088 | 38,535,168 = 2^25.199672 |
| **total** | | **2,001,469,440 = 2^30.898412** |

= 2^19.932628 units, plus <= 2^10.60 units of Miller-Rabin
certification (16 certifications x ~380 mulmods x 512 instructions)
and <= 2^15.579 units for the three comb tables: preprocessing total
<= 2^20.01 units, inside the 2^21 budget. A cap abort fails the run;
for the fixed production input it demonstrably does not fire.

## Correctness

A key match means [eta^(i + m1*j)]P0 and A_d are the same point; both
are powers of eta acting on P0, and eta has exact order N on P0, so
equality of points implies equality of exponents mod N: k1 =
(i + m1*j) mod N is exact including every wrap case. Phase 2 recovers
t exactly; verification closes the chain. The toy reproduces this
with zero wrong extractions, every extraction checked against the true
scalar.

## Success probability (exact integer arithmetic, no heuristic)

k1 is uniform over [0, N) (the d-th power map is exactly d-to-1 onto
the order-N subgroup). The run covers k1 iff some baby index i < m and
some block j in the window satisfy k1 = (m1*j + i) mod N - the
wrap-around hits near the end of the k1 range are included. Block j's
covered set is exactly the m values {(m1*j + i) mod N : i < m} (every
baby index is valid; the arc is injective), and consecutive blocks'
arcs are disjoint because their starts differ by m1 > m and no start
collision m1*Delta = 0 mod N can occur for window spans (m1*Delta < N).
Therefore EVERY window offset covers at least (m-1)*m values (dropping
the one block pair that may overlap through the mod-N wrap), i.e. at
least

  (m-1)*m / N >= 0.4 - 2^-118 >= 0.39

k1 values - an exact per-offset statement with no expectation
argument. The expectation over the uniform offset is >= the same
value; the toy brute-forces the exact wrap-inclusive expectation
(1334179/3316697 = 14675969/36483667 = 0.4022613461525126, matching
the third-round review's computation), and Part A crosschecks the
closed form against an exhaustive every-offset brute force on five
smaller parameter sets. Claimed lower bound: **0.39**.

## Work ledger (all charged, worst-case upper bounds)

| Component | Count | Units |
| --- | ---: | ---: |
| Baby scalars eta^i (one 512-instruction mod-r multiplication each, sequential) | m | m * 0.256 |
| Baby points [eta^i]P0 (comb: 32 window additions = 32 reference ops + 0.152 control; charged 2x = 64) | m | m * 64 |
| Baby sort: <= 117m compares + <= 118m moves (incl. copy-back) + key prep | m | m * 3.8355 |
| Giant scalars (one 512-instruction mod-r multiplication each, sequential) | m | m * 0.256 |
| Giant points [eta^(-m1*j)]A_d (comb: same construction, streamed) | m | m * 64 |
| Giant lookups: binary search, <= 117 probes x 32 instr | m | m * 1.872 |
| Window-start exponentiations (<= 2) | <= 2 | <= 196 |
| Phase-2 BSGS (comb mults + tiny sorted table) | 2^11 | 2^17 |
| Endgame: [zeta^(-k1)]A_1, verification, output | <= 4 mults + 1 exponentiation | <= 600 |
| Preprocessing: factor r-1 (three-term bound above), zeta/eta/zeta_e, three comb tables | one-off | <= 2^21 |

Per baby: 64 + 0.256 + 3.8355 = 68.0915 units.
Per giant: 64 + 0.256 + 1.872 = 66.128 units.

W(2^20) = m * 68.0915 + m * 66.128 + <= 2^21.9
        = m * 134.2195 + <= 2^21.9
        = 2^123.836032 (rounded up; exact rational 1934938/15625).

Preprocessing <= 2^20.01 units -> preprocessing_log2 = **21.0**,
included.

## Memory, by phase

- baby generation (streaming): window table ~3 MiB + block scratch
  <= 2^33 bytes + the m-entry record array (the table itself, by
  design);
- sort phase: main array 128m + merge buffer 128m = 256m =
  2^124.767581 bytes - the GLOBAL PEAK - plus code and scratch;
- giant generation + searches (streaming): sorted baby table 128m +
  block scratch <= 2^33 bytes + window table ~1 MiB: peak
  <= 128m x (1 + 2^-90), strictly below the sort phase; the toy asserts
  the block-local live arrays (live_accs = live_scalars = 64).

The global peak is the sort phase: 256*m bytes, strictly above
2^124.7675 including comb tables, code and scratch (all positive),
reported rounded up as memory_log2_bytes = **124.8**. Nonuniform
advice: **0 bytes** (everything constructed and charged at runtime).

## Alternatives and sensitivity

- r+1 variant: largest divisor of r+1 <= 2^20 is 484486 < 2^20.
  Rejected.
- d = 2^21 divides r-1 but A_{2^21} is not supplied at n = 2^20.
- Polynomial combinations of the supplied powers: p(s) uniform over
  Z_r for nonconstant p; birthday over them needs ~sqrt(r) samples.
  No gain.
- Sensitivity: charging full double-and-add (384 units) per point
  multiplication everywhere gives ~2^126.3; the score does not depend
  on the comb, but the comb construction (32 window additions + 0.152
  units control; charged 64) is stated so the charge is checkable.

## Evidence and reproduction

- Part A (production, exact): the deterministic hard-capped counted
  Brent-Pollard rho factorization (pre-batch-enforced 2^17-iteration
  cap, 12-cycle guard; measured 111,717 modular multiplications, 310
  gcds at a measured maximum of 153 binary-gcd rounds, 9 cycles)
  executes on the production r-1 and is asserted equal to the supplied
  factor list in both directions; the binary gcd is verified against
  math.gcd on 201 sampled pairs including (r-1, 254760293^2);
  2^20 | r-1; all factors < 2^28, Miller-Rabin prime, product equality;
  zeta = 7 generates Z_r*; eta has exact order N; zeta_e has exact
  order 2^20; d-th powers lie in the order-N subgroup; the per-offset
  minimum coverage is computed in exact integers and asserted
  >= 39/100 for every offset; the closed-form coverage bound is
  crosschecked against an exhaustive every-offset brute force on five
  smaller parameter sets; the parity-key uniqueness property is
  asserted on the baby set and sampled on negation pairs.
- Part B (toy, complete two-phase attack, d = 6 | n-1, N_toy = 174563,
  n-1 factored by the same capped construction): the declared run
  `{"trials": 200, "mode": "partial"}` reports 88/200 = 0.44 vs the
  exact wrap-inclusive expectation 0.4022613461525126, with 0 wrong
  extractions (every extraction checked against the true scalar); the
  toy runs the streaming comb (5 blocks of 64; live_accs =
  live_scalars = 64 asserted), the explicit merge sort (compare/move
  upper bounds asserted), the explicit binary search (probe budget
  floor(log2 m_toy)+1 = 9 asserted per lookup), and asserts equality
  with plain double-and-add for EVERY point (all babies and every
  giant of every trial). Full-coverage control: 200/200, 0 wrong
  extractions.
- Reproduce:
  `echo '{"trials": 200, "mode": "partial"}' | python3 submissions/bls12-381-g1-ecdlp/partial-cheon-v13/experiment.py`.

Small experiments are evidence only at their actual scope.

## Heuristics, limitations, objections

- No heuristic anywhere in the success argument: exact per-offset
  integer coverage including wrap, minimum 0.4 - 2^-118; the
  factorization is deterministic on the fixed public input; no
  probabilistic quantity appears anywhere in the claim.
- The table machinery is deterministic with worst-case upper bounds
  enforced by assertion in the experiment (merge-sort compares/moves
  including the copy-back; binary-search probe budget; block-local
  live array lengths; binary-gcd round bound per call; rho iteration
  and cycle caps before each consumed batch).
- Every instruction allowance is either an executed construction
  (binary gcd, merge sort, binary search, comb) or carries explicit
  limb-level line items (512-instruction mod-r multiplication, 15-,
  50-, 16-, 32-instruction table operations, 48-instruction gcd
  rounds, 1056-instruction rho main-loop iterations); every allowance
  carries a >= 2x margin or an enforced bound.
- Soft numbers are instruction allowances with stated constructions;
  the sensitivity row shows the claim survives 384-unit point
  multiplications.
- Single shot at the 0.39 threshold; deterministic full coverage costs
  m1/m times more in the leading term (about +0.86 log2).
- The bound applies only for n >= 2^20.
- Standard cryptanalysis (Cheon); the contribution is the sequential
  window structure, the wrap-exact coverage computation, and a
  hashless, move-charged, worst-case-exact, fully executed sorted-array
  design under this cost model.

References: `tracks/docs/COST_MODEL.md`; J. H. Cheon, "Discrete
logarithm problems with auxiliary inputs" (ASIACRYPT 2008 / J.
Cryptology 2010); the track's protected rules.md; HAC ch. 3.6 (BSGS).
Self-contained above.
