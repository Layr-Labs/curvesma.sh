# bls12-381-g1-ecdlp: partial-coverage Cheon attack v11, at n=2^20

Status: ready for independent exploratory review. Classical single-shot
attack using the supplied auxiliary power A_d, d = 2^20; no quantum steps,
no uncharged advice, no amortization over targets. Cheon's
discrete-logarithm-with-auxiliary-input attack is standard cryptanalysis.
v9 is v7 resubmitted: v7 and its twin v8 both failed server-side during
a platform incident (workflow dispatch; judgment null on both; other
solvers' submissions failed simultaneously; no review occurred).
v8 additionally fixed a malformed score-table header in the note. v8 responds to
the sixth-round review of v6 (which accepted the reduction, the parity-key
uniqueness, the wrap-aware coverage argument supporting the 0.39
threshold, and the production parameter checks, and rejected the
remaining resource accounting):

1. **Bounded, deterministic factorization construction.** v10 states
   the construction the implementation actually runs. The input r-1 is
   the fixed public BLS12-381 parameter and the seed sequence is fixed,
   so the computation is DETERMINISTIC: the hosted run reproduces the
   measured 111,717 modular multiplications and 310 gcds over 9 rho
   cycles (Part A reports the counts); no probabilistic model is used
   anywhere in this claim (the v7-v9 birthday-tail subtraction is
   withdrawn). The hard caps bound the worst case for adversarial
   inputs and are now enforced exactly as stated: each cycle checks its
   remaining iteration budget BEFORE consuming each batch (halving
   rounds, 128-step mini-batches, and the backtrack scan), so a
   COMPLETED cycle provably used at most the cap - the seventh-round
   review correctly noted the v9 assertion fired only after a whole
   doubling round; and a hard 12-cycle guard bounds the total. The
   per-iteration charge is corrected to the implementation's actual
   work: TWO modular multiplications per inner-loop iteration
   (2 x 512 instructions) plus 32 instructions of loop, index and
   batching work = 1056 instructions (the v9 allowance of 520
   undercounted the double-mulmod inner loop, as the review noted).
   Worst case: 12 cycles x 2^17 iterations x 1056 instructions, plus
   <= 128 backtrack steps per cycle at <= 1536 instructions and <= 16
   Miller-Rabin certifications: <= 2^30.63 instructions = 2^19.63 units,
   inside the 2^21 preprocessing budget. A cap abort fails the run; for
   the fixed production input it demonstrably does not fire.
2. **Instruction-level line items, asserted abstract counts.** The
   experiment asserts the abstract counts (compares, moves, probes,
   additions, inversions); every instruction allowance below now
   carries a complete line-item breakdown:
   - Binary-search probe (32 instructions): 2-word midpoint computation
     (6: add, shift, carry), record pointer load + record load (4),
     six-limb key load and compare with early exit (12), parity-bit
     compare (2), interval update + termination test + branches (8).
     Worst case exactly floor(log2 m) + 1 = 117 probes: <= 3744
     instructions = 1.872 units per lookup, enforced per lookup in the
     experiment.
   - Merge-sort compare (15): six 8-byte limb loads + compares + branch;
     record move (50): sixteen loads + sixteen stores + four index
     updates + control; key preparation (16): limb loads, parity extract,
     pack. All upper bounds asserted.
   - Comb control: window table addressing (<= 192 instructions), loop
     control (<= 96), accumulator initialization (<= 16): <= 304
     instructions = 0.152 units per multiplication against the 32 units
     of slack inside the 64-unit charge (2x the 32.152-unit
     construction).
   - The toy no longer multiplies points by plain double-and-add: the
     baby and giant point sets are computed by the charged construction
     itself - windowed additions with every window layer's affine
     additions across ALL points sharing ONE batched inversion
     (Montgomery's trick), equality with plain double-and-add asserted
     in full for every baby point (once, outside the trial loop), in
     full for every giant point on trial 0, and on 8 random giants every
     subsequent trial, and per-multiplication addition counts and
     per-layer inversion counts reported (baby side: 499 additions, 2
     inversions for all m babies; giant side: ~507 additions and exactly
     2 inversions per trial). An intermediate accumulator hitting
     infinity (window-prefix scalar ~ 0 mod r, probability <= 2^-112 per
     run) is handled exactly: the comb continues from the identity.
3. **Exact-count claims corrected to upper bounds.** The review correctly
   observed that a merge of two nonempty lists uses strictly fewer
   comparisons than their combined length, so "exactly 117*m key
   comparisons" was false; that ordinary binary search over m slots needs
   at most floor(log2 m) + 1 = 117 probes, not 118; and that the v5 toy
   performed an uncounted final copy-back. v6 states CONSERVATIVE UPPER
   BOUNDS throughout, charges the copy-back, and enforces every bound by
   assertion in the experiment:
   - The comb generation now STREAMS (the eighth-round review correctly
     noted that v10's toy still allocated the accumulator array for all
     m scalars and materialized every giant scalar): giant scalars are
     generated and consumed block by block (never materialized as an
     m-sized list), each block's accumulator array is allocated for
     EXACTLY the block size, and each block's giant points are
     binary-searched immediately and then discarded - so the only
     m-sized live object during the giant phase is the sorted baby
     table. Baby scalars stream the same way into the (by design
     m-sized) record array. Live working storage per block: block
     scalars 32 B, affine accumulators 96 B, Montgomery prefix
     products 48 B, pending-operation entries ~72 B: <= 256 B per
     record x 2^24 = <= 2^32 bytes. The merge buffer is freed after the
     sort. The global memory peak is therefore exactly the sort phase
     (see Memory). The toy DEMONSTRATES the O(block) property: the
     per-block live accumulator and scalar array lengths are recorded
     and asserted (5 blocks of 64; live_accs = live_scalars = 64;
     max_pending 63), with equality to double-and-add asserted for
     every baby point and trial 0's giants plus 8 random giants per
     block on later trials.
   - Merge sort (bottom-up, auxiliary buffer): ceil(log2 m) = 117 passes;
     per pass each record is compared at most once (15 instructions) and
     moves at most once (one 128-byte load + store, 50 instructions), and
     the sort may end with at most one full-array copy-back (m more
     moves): total <= m * 117 * 15 + m * 118 * 50 instructions, plus
     <= 16 instructions per record for the one-time (x, parity) key
     preparation (six limb loads, parity extract, pack). Charged
     7671 instructions = 3.8355 units per baby, an upper bound asserted
     in the experiment.
   - Binary search: worst case EXACTLY floor(log2 m) + 1 = 117 probes
     (the candidate range at most halves each probe and empties on the
     last; for m < 2^117 this gives 117); enforced per lookup in the
     experiment. Each probe: 2-word mid computation (6), record pointer
     load + record load (4), six-limb key load + compare with early
     exit (12), parity bit (2), interval update + termination test +
     branches (8) = 32 instructions: <= 3744 instructions = 1.872
     units per lookup, worst case, including multiword (117-bit) index
     arithmetic and addressing.
2. **Memory: phase analysis.** The seventh-round review correctly
   observed that during giant generation the sorted baby table stays
   resident while the unblocked layer-batched comb held m-sized giant
   scalars (32m), affine accumulators (96m) and Montgomery prefix
   products (48m) - 304m bytes ~ 2^125.02, exceeding the claimed
   2^124.8. v10 blocks the comb generation, so the m-sized live objects
   are exactly the two arrays during the sort and the single sorted
   table afterwards:
   - baby generation: window table ~3 MiB + block scratch <= 2^32 bytes;
   - sort phase: main array 128m + merge buffer 128m = 256m =
     2^124.767581 bytes (the global peak), plus code and scratch;
   - giant generation + searches (streaming): sorted baby table 128m +
     block scratch <= 2^32 bytes + window table ~1 MiB: peak
     <= 128m x (1 + 2^-92), strictly below the sort phase; the toy
     asserts the block-local live arrays (live_accs = live_scalars =
     64 at toy scale).
   The global peak is the sort phase: 256*m bytes, strictly above
   2^124.7675 including comb tables, code and scratch (all positive),
   reported rounded up as memory_log2_bytes = 124.8 - the same figure as
   v7-v9, now with the phase analysis that sustains it.
3. **Complete constructions for the charged allowances.**
   - Fixed-base comb: each window addition is EXACTLY ONE reference
     operation by construction of the cost model (the pinned reference is
     the batched incomplete affine G1 addition at 6M, S = M, A = 0; the
     inversion batching is part of the reference definition - the m baby
     multiplications proceed in lockstep, so the additions of each
     window layer are independent and amortize their inversions by
     Montgomery's trick exactly as the reference op assumes). Per
     multiplication: 32 window additions (32 reference units) + window
     table addressing and loop control (<= 8 units by count); charged 64
     units (2x). The point-key preparation is charged separately (16
     instructions per record, above).
   - Mod-r multiplication: 512 instructions (explicit 5x52-limb
     schoolbook + Montgomery reduction construction, unchanged).
   - Runtime factorization of r-1, a hard-capped, counted construction,
     performed at production scale by the experiment: Brent-Pollard rho
     with Miller-Rabin certification; each cycle capped at 2^17
     iterations (<= 520 instructions per iteration including batched
     gcd work), <= 12 cycles: worst case <= 2^29.61 instructions =
     2^18.64 units, plus <= 2^10.3 units of certification - inside the
     preprocessing 2^21 budget. The implementation counts its own
     modular multiplications, enforces the cap, and reports the
     measured count (111,717 for the production r-1). Part A EXECUTES
     this construction on the production r-1 and asserts the result
     equals the supplied factor list in both directions; the toy factors
     n-1 with the same routine. Failure model: birthday tail
     exp(-beta^2/2), beta = 6.57 at the worst factor; <= 2^-27.5 over
     the run, subtracted from the success probability.
4. Everything accepted in earlier rounds is unchanged: the reduction, the
   contiguous window, the sequential giants, the parity-keyed sorted
   array (same-x pairs have opposite y-parities, so a key match is an
   exact point match and k1 = (i + m1*j) mod N is exact), the
   wrap-inclusive exact coverage, and the every-offset crosschecks on
   five smaller parameter sets.

## Claimed bounds

| Field | Value |
| --- | --- |
| time_log2 (score) = log2(W(2^20)) | **123.836032** (= 1934938/15625, rounded up) |
| success_probability (supported lower bound) | 0.39 (exact per-offset minimum (m-1)m/N >= 0.4 - 2^-118) |
| memory_log2_bytes | 124.8 (main array + merge buffer) |
| preprocessing_log2 (included in total) | 21.0 |
| nonuniform_advice_log2_bytes | 0 |

W(n) for n >= 2^20 is constant: the attack reads only A_1 = T and
A_d = [s^d]P0 with d = 2^20 <= n (at the scoring n = 2^20, A_d is the last
element, index d-1). The bound does not apply for n < 2^20.

## Conventions

Protected target: BLS12-381 G1, generator P0 of prime order r, challenge
s uniform in [1, r-1), T = A_1 = [s]P0, powers A_i = [s^i mod r]P0 for
i = 1..n, scored at n = 2^20. Cost model `classical-group-work-v2`:
cost = H + R/2000; H = reference operations (batched incomplete affine G1
addition, 6M, one reference operation = 1 unit), R = other scalar 64-bit
AArch64 instruction equivalents. Pinned base-field multiplication: 383
instructions (blst `blst_fp_mul`, mod p). Mod-r multiplication (r 255
bits): 5x52-limb schoolbook (25 limb products of 52x52 -> 104 bits via
64x64->128 instructions, plus carry accumulation) followed by Montgomery
reduction of the same shape: charged **512 instructions** (0.256 units),
used once per baby/giant scalar update and in at most three
exponentiations per run (each <= 255 squarings + 128 multiplications plus
loop work, charged). A reference operation and its implementation are
never double-counted.

## Number-theoretic setup (verified by experiment.py Part A)

- The complete prime factorization of r - 1 is verified in the
  experiment (R_M1_FACTORS, hex literals; every factor Miller-Rabin
  prime, product equality asserted exactly):
  r - 1 = 2^32 * 3 * 11 * 19 * 10177 * 125527 * 859267 * 906349^2 *
  2508409 * 2529403 * 52437899 * 254760293^2; every prime factor is
  < 2^28, so the runtime re-derivation (trial division/Pollard rho,
  <= 2^21 units, charged) is cheap.
- 2^20 | r-1, and d = 2^20 is the largest divisor of r-1 not exceeding the
  scoring n (any divisor with an odd part is at most 2^18 * 3 < 2^20).
- zeta = 7 generates Z_r*; eta = zeta^(2^20) has exact order
  N = (r-1)/2^20 (log2 N = 234.857); zeta_e = zeta^N has exact order 2^20.

## Why the algebra works

s uniform in [1, r-1) lies in Z_r*. The d-th power map on Z_r* (order
r-1, d | r-1) is exactly d-to-1 onto the unique order-N subgroup, so s^d
is uniform there and s^d = eta^k1 for unique k1 uniform in [0, N); and
s = zeta^(k1 + t*N) for unique t in [0, 2^20). Hence
A_d = [eta^k1]P0, Q' = [zeta^(-k1)]A_1 = [zeta_e^t]P0.

## Algorithm v5 (single shot; the only coin is the window offset)

0. Read A_1 and A_d (charged).
1. Runtime factorization of r-1; construct zeta, eta, zeta_e; build
   fixed-base comb tables for P0, A_d, A_1 (each 8,160 points,
   <= 2^14 reference additions, ~3 MiB, charged).
2. Phase 1, partial-coverage BSGS with a contiguous window:
   m1 = ceil(sqrt(N)), M = ceil(N/m1) (= m1 here), and
   m = isqrt(floor(2N/5)) + 1 = 141430898198098033869460368287841487
   (pure integer arithmetic; m < m1).
   - Baby scalars: for i = 0..m-1 compute eta^i sequentially (one mod-r
     multiplication, 512 instructions, per step) and the point
     [eta^i]P0 (one comb multiplication, 64 units).
   - Baby table: sort the m records (key = (x, y mod 2); payload = the
     baby index i) by bottom-up merge sort with an auxiliary buffer:
     117 passes, at most m * 117 key compares, at most m * 118 record
     moves (including one possible final copy-back) and <= 16
     instructions per record of key preparation: <= 7671 instructions
     per baby, an upper bound asserted in the experiment.
   - Window: uniform offset o in [0, M); the m consecutive blocks
     o..o+m-1 (mod M) are processed as at most two contiguous runs
     (run 2 starts at block 0 with scalar 1); each giant scalar is one
     sequential mod-r multiplication; at most two window-start
     exponentiations (charged).
   - Giants: point [eta^(-m1*j)]A_d per block j in the window (one comb
     multiplication); binary search of the sorted baby array for the key
     (x(Y), y(Y) mod 2): at most floor(log2 m) + 1 = 117 slot probes
     (exact worst case, enforced per lookup in the experiment), each
     <= 30 instructions including multiword index arithmetic and
     addressing: <= 3510 instructions per lookup. On a hit at baby index
     i: k1 = (i + m1*j) mod N - exact, because the key match is an exact
     point match (see above).
3. Phase 2 (deterministic, full coverage): Q' = [zeta^(-k1)]A_1; BSGS
   for t in [0, 2^20) with 2^11 comb multiplications and a 2^10-entry
   sorted baby table searched identically (all charged).
4. s = zeta^(k1 + t*N) mod r; verify [s]P0 = T; output. No phase-1 hit ->
   fail; no restart or amplification.

Correctness: a key match means [eta^(i + m1*j)]P0 = s^d-mod-r component
as scalars in Z_r* acting on P0 (both sides are powers of eta times P0;
equality of points implies equality of exponents mod N since eta has
exact order N on P0), so k1 = (i + m1*j) mod N is exact including every
wrap case; phase 2 recovers t exactly; verification closes the chain.
The toy reproduces this with zero wrong extractions.

## Success probability (exact integer arithmetic, no heuristic)

k1 is uniform over [0, N). The run covers k1 iff some baby index i < m
and some block j in the window satisfy k1 = (m1*j + i) mod N - the
wrap-around hits near the end of the k1 range are included. Block j's
covered set is exactly the m values {(m1*j + i) mod N : i < m} (every
baby index is valid; the arc is injective), and consecutive blocks' arcs
are disjoint because their starts differ by m1 > m and no start collision
m1*Delta = 0 mod N can occur for window spans (m1*Delta < N). Therefore
EVERY window offset covers at least (m-1)*m values (dropping the one
block pair that may overlap through the mod-N wrap), i.e. at least

  (m-1)*m / N >= 0.4 - 2^-118

k1 values - an exact per-offset statement with no expectation argument.
The expectation over the uniform offset is >= the same value; the toy
brute-forces the exact wrap-inclusive expectation
(1334179/3316697 = 14675969/36483667 = 0.4022613461525126, matching the
third-round review's computation), and Part A crosschecks the closed form
against an exhaustive every-offset brute force on five smaller parameter
sets - a fully deterministic statement, with no probabilistic quantity
anywhere in the claim:

  0.4 - 2^-118 >= 0.39.

Claimed lower bound: **0.39**.

## Work ledger (all charged, worst-case exact)

| Component | Count | Units |
| --- | ---: | ---: |
| Baby scalars eta^i (one 512-instruction mod-r multiplication each) | m | m * 0.256 |
| Baby points [eta^i]P0 (comb: 32 window additions = 32 reference ops + 0.152 control; charged 2x = 64) | m | m * 64 |
| Baby sort: merge sort, <= 117m compares + <= 118m moves (incl. copy-back) + key prep | m | m * 3.8355 |
| Giant scalars (one 512-instruction mod-r multiplication each, sequential) | m | m * 0.256 |
| Giant points [eta^(-m1*j)]A_d (comb: same construction, exercised by the toy) | m | m * 64 |
| Giant lookups: binary search, <= 117 probes x 32 instr (mid 6, loads 4, key 12, parity 2, loop 8) | m | m * 1.872 |
| Window-start exponentiations (<= 2) | <= 2 | <= 196 |
| Phase-2 BSGS (comb mults + tiny sorted table) | 2^11 | 2^17 |
| Endgame: [zeta^(-k1)]A_1, verification, output | <= 4 mults + 1 exponentiation | <= 600 |
| Preprocessing: factor r-1 (hard-capped counted Brent-rho, executed in Part A), zeta/eta/zeta_e, three comb tables | one-off | <= 2^21 |

Per baby: 64 + 0.256 + 3.8355 = 68.0915 units.
Per giant: 64 + 0.256 + 1.872 = 66.128 units.

W(2^20) = m * 68.0915 + m * 66.128 + <= 2^21.9
        = m * 134.2195 + <= 2^21.9
        = 2^123.836032 (rounded up; exact rational 1934938/15625).

Preprocessing <= 2^21 units -> preprocessing_log2 = **21.0**, included.

## Memory

Main sorted array: m entries x 128 bytes (key: x 48 bytes + 1 parity bit
padded to 56, payload index 56, metadata 16) plus the merge buffer: m
entries x 128 bytes: together exactly 256 * m = 2^124.767581 bytes; comb
tables ~3 MiB, code and scratch positive on top. Peak storage is strictly
above 2^124.7675 bytes, reported rounded up as memory_log2_bytes =
**124.8**. Nonuniform advice: **0 bytes** (everything constructed and
charged at runtime).

## Alternatives and sensitivity

- r+1 variant: largest divisor of r+1 <= 2^20 is 484486 < 2^20. Rejected.
- d = 2^21 divides r-1 but A_{2^21} is not supplied at n = 2^20.
- Polynomial combinations of the supplied powers: p(s) uniform over Z_r
  for nonconstant p; birthday over them needs ~sqrt(r) samples. No gain.
- Sensitivity: charging full double-and-add (384 units) per point
  multiplication everywhere gives ~2^126.3; the score does not depend on
  the comb, but the comb construction (32 window additions, each exactly
  one reference operation by the cost model's definition, plus <= 8 units
  control; charged 64) is stated so the charge is checkable. The v6 score
  contains every sort/lookup/prepare instruction charge demanded by the
  fourth- and fifth-round reviews, as upper bounds.

## Evidence and reproduction

- Part A (production, exact): the runtime factorization CONSTRUCTION
  (deterministic, hard-capped, counted Brent-Pollard rho + Miller-Rabin;
  111,717 measured modular multiplications, 310 gcds, 9 cycles against
  the pre-batch-enforced 2^17-per-cycle cap and the 12-cycle guard) is
  executed on the production r-1 and asserted equal to the supplied
  factor list in both directions;
  2^20 | r-1; all factors < 2^28 and Miller-Rabin prime with product
  equality; zeta = 7 generates Z_r*; eta has exact order N; zeta_e has
  exact order 2^20; d-th powers of random scalars lie in the order-N
  subgroup; the per-offset minimum coverage is computed in exact
  integers and asserted >= 39/100 for every offset; the closed-form
  coverage bound is additionally crosschecked against an exhaustive
  every-offset brute force on five smaller (N, m1, m) parameter sets; the
  parity-key uniqueness property is asserted on the baby set and sampled
  on negation pairs.
- Part B (toy, complete two-phase attack, d = 6 | n-1, N_toy = 174563,
  production-shaped sequential window giants computed by the charged
  layer-batched comb construction itself - every window layer's affine
  additions share ONE batched inversion, equality with plain
  double-and-add asserted for every point, addition and inversion counts
  reported - wrap-inclusive exact coverage, parity-keyed sorted array
  with the explicit merge sort and binary search, the merge compare/move
  upper bounds and the floor(log2 m_toy)+1 = 9 probe budget enforced
  per lookup, and n-1 factored by the same capped Brent-rho
  construction): declared run `{"trials": 200, "mode": "partial"}`:
  88/200 = 0.44 vs the exact wrap-inclusive expectation
  0.4022613461525126 (equal to the third-round review's corrected
  value); 0 wrong extractions; every extraction verified by brute force.
  Full-coverage control: 200/200, 0 wrong extractions.
- Reproduce:
  `echo '{"trials": 200, "mode": "partial"}' | python3 submissions/bls12-381-g1-ecdlp/partial-cheon-v11/experiment.py`.

Small experiments are evidence only at their actual scope.

## Heuristics, limitations, objections

- No heuristic anywhere in the success argument: exact per-offset integer
  coverage including wrap, minimum 0.4 - 2^-118.
- The table machinery is deterministic with worst-case upper bounds
  asserted in the experiment (merge sort compares/moves incl. the final
  copy-back; binary-search probe budget floor(log2 m)+1; parity-key
  uniqueness theorem); no hash function, no randomness, no probe
  distribution, no uncharged moves, addressing or key preparation.
- Every instruction allowance carries a complete line-item breakdown
  and a >= 2x margin or an enforced upper bound: 64-unit comb charge
  (2x the stated 32-addition + 0.152-control construction, exercised by
  the toy), 512-instruction mod-r multiplication (5x52 schoolbook +
  Montgomery), 15-instruction key compares, 50-instruction record moves,
  32-instruction binary-search probes (mid 6 + loads 4 + key 12 +
  parity 2 + interval update/termination/branches 8), 16-instruction
  key preparation, and the hard-capped counted Brent-rho factorization
  (worst case 2^18.64 units; measured 111,717 mulmods at production
  scale; cap-failure probability <= 2^-27.5 subtracted from the success
  probability).
- Single shot at the 0.39 threshold; deterministic full coverage costs
  m1/m times more in the leading term (about +0.86 log2).
- The bound applies only for n >= 2^20.
- Standard cryptanalysis (Cheon); the contribution is the sequential
  window structure, the wrap-exact coverage computation, and a
  hashless, move-charged, worst-case-exact sorted-array design under this
  cost model.

References: `tracks/docs/COST_MODEL.md`; J. H. Cheon, "Discrete logarithm
problems with auxiliary inputs" (ASIACRYPT 2008 / J. Cryptology 2010);
the track's protected rules.md; HAC ch. 3.6 (BSGS). Self-contained above.