# bls12-381-g1-ecdlp: partial-coverage Cheon attack v5, at n=2^20

Status: ready for independent exploratory review. Classical single-shot
attack using the supplied auxiliary power A_d, d = 2^20; no quantum steps,
no uncharged advice, no amortization over targets. Cheon's
discrete-logarithm-with-auxiliary-input attack is standard cryptanalysis.
v5 responds to the fourth-round review of v4 (which accepted the reduction,
the contiguous-window construction, the sequential giants, the sorted-array
replacement for hashing, the wrap-exact coverage, and the 88/200 toy run,
and rejected the resource accounting):

1. **Sort/lookup ledger completeness.** The v4 ledger charged only key
   comparisons. Merge sort also moves every record once per pass and
   performs addressing and control; binary search performs index
   arithmetic, addressing and loads. v5 charges all of it explicitly:
   - Merge sort (bottom-up, with an auxiliary buffer): ceil(log2 m) = 117
     passes; each record moves exactly once per pass (one 128-byte load +
     one 128-byte store, charged 50 instructions) and is compared at most
     once per pass during the merges that touch it, charged 15
     instructions per key compare: exactly 117 * (50 + 15) = 7605
     instructions = 3.8025 units per baby, worst case.
   - Binary search: exactly ceil(log2(2m)) = 118 slot probes, each one
     381-bit key compare (15 instructions) plus index arithmetic,
     addressing and the load (10 instructions): 2950 instructions = 1.475
     units per lookup, worst case.
   - Fixed-base comb multiplication: 32 window additions of 381-bit affine
     points (each one batched incomplete affine reference addition, 6M,
     plus window addressing and control), i.e. 32-40 units by count;
     charged 64 units (2x). The v4 charge was the same number; v5 states
     the underlying operation count.
   - Mod-r multiplication: 512 instructions (explicit 5x52-limb
     schoolbook + Montgomery reduction construction, unchanged).
   The v5 score is HIGHER than v4's (123.831402 vs 123.792992) precisely
   because these costs are now in the ledger.
2. **Merge-sort memory.** Ordinary buffered merge sort needs linear
   auxiliary storage. v5 allocates the merge buffer explicitly: the main
   sorted array (m entries x 128 bytes) plus the merge buffer (m entries x
   128 bytes) = 2^124.77 bytes, plus comb tables (~3 MiB), code and
   scratch; the peak is strictly above 2^124.77 bytes, reported rounded up
   as memory_log2_bytes = 124.8.
3. **The equal-x collision.** Equal x-coordinates establish equality only
   up to point negation (an x determines y up to sign), so a v4 x-keyed
   lookup could hit the negated baby [eta^{i+N/2}]P0, breaking the exact
   k1 formula. v5 keys every record (baby and giant) by the pair
   (x, y mod 2): on this curve y^2 = x^3 + 4 over odd Fp, the (at most
   two) points sharing an x have y and p - y, opposite parities, so
   (x, y mod 2) uniquely determines the point. Baby and giant records are
   both raw affine points (no canonicalization machinery is involved
   anywhere in the Cheon attack); a key match is therefore an exact point
   match and k1 = (i + m1*j) mod N is exact. The uniqueness property is
   proved in one line in the experiment docstring, asserted at production
   parameters in Part A (the baby set has m distinct keys), and exercised
   by the toy run, whose table is now the parity-keyed sorted array built
   by the explicit merge sort and searched by the explicit binary search
   (the v4 toy used a full-(x,y) dictionary, which did not exercise the
   key).
4. **Coverage, independently validated.** The v4 Part A computed the
   coverage formula's value but did not independently validate per-offset
   behavior. v5 Part A additionally brute-forces EVERY window offset of
   five smaller (N, m1, m) parameter sets (same construction as
   production) and asserts the closed-form bound (m-1)*m/N holds at every
   offset of each; the production claim itself remains the exact
   per-offset statement (m-1)*m/N >= 0.4 - 2^-118 >= 0.39.

## Claimed bounds

| Field | Value |
| --- | --- |
| time_log2 (score) = log2(W(2^20)) | **123.831402** (= 61915701/500000, rounded up) |
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
     exactly m * ceil(log2 m) key compares and m * ceil(log2 m) record
     moves in 117 passes, 7605 instructions per baby, all charged.
   - Window: uniform offset o in [0, M); the m consecutive blocks
     o..o+m-1 (mod M) are processed as at most two contiguous runs
     (run 2 starts at block 0 with scalar 1); each giant scalar is one
     sequential mod-r multiplication; at most two window-start
     exponentiations (charged).
   - Giants: point [eta^(-m1*j)]A_d per block j in the window (one comb
     multiplication); binary search of the sorted baby array for the key
     (x(Y), y(Y) mod 2): exactly 118 slot probes, 2950 instructions,
     worst case. On a hit at baby index i: k1 = (i + m1*j) mod N -
     exact, because the key match is an exact point match (see 3 above).
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

  (m-1)*m / N >= 0.4 - 2^-118 >= 0.39

k1 values - an exact per-offset statement with no expectation argument.
The expectation over the uniform offset is >= the same value; the toy
brute-forces the exact wrap-inclusive expectation
(1334179/3316697 = 14675969/36483667 = 0.4022613461525126, matching the
third-round review's computation), and Part A crosschecks the closed form
against an exhaustive every-offset brute force on five smaller parameter
sets. Claimed lower bound: **0.39**.

## Work ledger (all charged, worst-case exact)

| Component | Count | Units |
| --- | ---: | ---: |
| Baby scalars eta^i (one 512-instruction mod-r multiplication each) | m | m * 0.256 |
| Baby points [eta^i]P0 (comb multiplications, 32-40 units by count) | m | m * 64 |
| Baby sort: merge sort, 117 passes x (50 move + 15 compare) instr | m | m * 3.8025 |
| Giant scalars (one 512-instruction mod-r multiplication each, sequential) | m | m * 0.256 |
| Giant points [eta^(-m1*j)]A_d (comb multiplications) | m | m * 64 |
| Giant lookups: binary search, 118 probes x 25 instr | m | m * 1.475 |
| Window-start exponentiations (<= 2) | <= 2 | <= 196 |
| Phase-2 BSGS (comb mults + tiny sorted table) | 2^11 | 2^17 |
| Endgame: [zeta^(-k1)]A_1, verification, output | <= 4 mults + 1 exponentiation | <= 600 |
| Preprocessing: factor r-1, zeta/eta/zeta_e, three comb tables | one-off | <= 2^21 |

Per baby: 64 + 0.256 + 3.8025 = 68.0585 units.
Per giant: 64 + 0.256 + 1.475 = 65.731 units.

W(2^20) = m * 68.0585 + m * 65.731 + <= 2^21.9
        = m * 133.7895 + <= 2^21.9
        = 2^123.831402 (rounded up; exact rational 61915701/500000).

Preprocessing <= 2^21 units -> preprocessing_log2 = **21.0**, included.

## Memory

Main sorted array: m entries x 128 bytes (key: x 48 bytes + 1 parity bit
padded to 56, payload index 56, metadata 16) = 2^123.77 bytes; merge
buffer: m entries x 128 bytes = 2^123.77 bytes; comb tables ~3 MiB; code
and scratch positive. Peak storage is strictly above 2^124.77 bytes,
reported rounded up as memory_log2_bytes = **124.8**. Nonuniform advice:
**0 bytes** (everything constructed and charged at runtime).

## Alternatives and sensitivity

- r+1 variant: largest divisor of r+1 <= 2^20 is 484486 < 2^20. Rejected.
- d = 2^21 divides r-1 but A_{2^21} is not supplied at n = 2^20.
- Polynomial combinations of the supplied powers: p(s) uniform over Z_r
  for nonconstant p; birthday over them needs ~sqrt(r) samples. No gain.
- Sensitivity: charging full double-and-add (384 units) per point
  multiplication everywhere gives ~2^126.3; the score does not depend on
  the comb, but the comb count (32-40 window additions) is stated so the
  charge is checkable. The v5 score already contains every sort/lookup
  instruction charge demanded by the fourth-round review.

## Evidence and reproduction

- Part A (production, exact): 2^20 | r-1; factorization verified
  (Miller-Rabin + product equality, all factors < 2^28); zeta = 7
  generates Z_r*; eta has exact order N; zeta_e has exact order 2^20;
  d-th powers of random scalars lie in the order-N subgroup; the
  per-offset minimum coverage is computed in exact integers and asserted
  >= 39/100 for every offset; the closed-form coverage bound is
  additionally crosschecked against an exhaustive every-offset brute
  force on five smaller (N, m1, m) parameter sets; the parity-key
  uniqueness property is asserted on the baby set and sampled on
  negation pairs.
- Part B (toy, complete two-phase attack, d = 6 | n-1, N_toy = 174563,
  production-shaped sequential window giants, wrap-inclusive exact
  coverage, parity-keyed sorted array with the explicit merge sort and
  binary search): declared run `{"trials": 200, "mode": "partial"}`:
  88/200 = 0.44 vs the exact wrap-inclusive expectation 0.4022613461525126
  (equal to the third-round review's corrected value); 0 wrong
  extractions; every extraction verified by brute force. The toy binary
  search enforces its probe budget (ceil(log2(2m_toy)) = 10) as a hard
  assert per lookup. Full-coverage control: 200/200,
  0 wrong extractions.
- Reproduce:
  `echo '{"trials": 200, "mode": "partial"}' | python3 submissions/bls12-381-g1-ecdlp/partial-cheon-v5/experiment.py`.

Small experiments are evidence only at their actual scope.

## Heuristics, limitations, objections

- No heuristic anywhere in the success argument: exact per-offset integer
  coverage including wrap, minimum 0.4 - 2^-118.
- The table machinery is deterministic and worst-case exact (merge sort
  pass count, binary-search probe count, parity-key uniqueness theorem);
  no hash function, no randomness, no probe distribution, no uncharged
  moves or addressing.
- Soft numbers are instruction sub-costs, each with an explicit
  construction: 64-unit comb charge (2x the 32-40 window-addition count),
  512-instruction mod-r multiplication (5x52 schoolbook + Montgomery),
  15-instruction 381-bit compares, 50-instruction record moves,
  25-instruction binary-search probes. The sensitivity row shows the
  claim survives 384-unit point multiplications.
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