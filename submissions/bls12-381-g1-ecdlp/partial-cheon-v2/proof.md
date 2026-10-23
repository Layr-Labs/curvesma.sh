# bls12-381-g1-ecdlp: partial-coverage Cheon attack v2, at n=2^20

Status: ready for independent exploratory review. Classical single-shot
attack using the supplied auxiliary power A_d, d = 2^20; no quantum steps, no
uncharged advice, no amortization over targets. Cheon's
discrete-logarithm-with-auxiliary-input attack is standard cryptanalysis;
this v2 replaces the rejected `partial-cheon` package and responds to the
recorded review objections:

- The v1 ledger charged "one modular multiplication per giant scalar" while
  its own experiment evaluated a modular exponentiation per giant, because
  the giant block indices were a random subset. v2 replaces the random
  subset with a contiguous window of m blocks at a uniform random offset:
  the giant scalars are now genuinely sequential (one mod-r multiplication
  per giant, plus at most two exponentiations per run for the window
  start), and the toy experiment implements exactly this.
- The v1 claim asserted an "exactly 0.40" coverage expectation without an
  exact computation (ceilings and the shortened final block were
  unaccounted). v2 computes the coverage in exact integer arithmetic,
  prints it at production parameters (it is exactly 2/5), and bounds the
  worst-case single-run coverage at 0.4 - 2^-118.
- Subset generation, baby-table insertion and hash-lookup costs are now
  explicitly constructed and charged; the comparison to the sibling rho
  track's bound is removed from the claim (it is context in the note only).

## Claimed bounds

| Field | Value |
| --- | --- |
| time_log2 (score) = log2(W(2^20)) | **123.777511** (= 123777511/1000000, rounded up) |
| success_probability (supported lower bound) | 0.39 (exact expectation 2/5) |
| memory_log2_bytes | 124.0 |
| preprocessing_log2 (included in total) | 21.0 |
| nonuniform_advice_log2_bytes | 0 |

W(n) for n >= 2^20 is constant: the attack reads only A_1 = T and
A_d = [s^d]P0 with d = 2^20 <= n (at the scoring n = 2^20, A_d is the last
element of the input array, index d-1). The score is log2(W(2^20)). The
bound does not apply for n < 2^20, where a smaller divisor of r-1 must be
used.

## Conventions

Protected target: BLS12-381 G1, generator P0 of prime order r, challenge
s uniform in [1, r-1), T = A_1 = [s]P0, powers A_i = [s^i mod r]P0 for
i = 1..n, scored at n = 2^20. Cost model `classical-group-work-v2`:
cost = H + R/2000; H = reference operations (batched incomplete affine G1
addition, 6M with S = M and A = 0; one reference operation = 1 unit),
R = other scalar 64-bit AArch64 instruction equivalents. Pinned base-field
multiplication: 383 instructions (blst `blst_fp_mul`, modulo the 381-bit
field prime p; the G1 reference price C = 2,000 is the rounded
6 x 383 = 2,298).

Prices for arithmetic modulo the scalar order r (255 bits) are supplied
analytically, since the pinned measurement covers only mod-p arithmetic:
- mod-r multiplication (5x52-limb schoolbook with 64x64->128 products plus
  Montgomery reduction): <= 384 instructions; used once per baby/giant
  scalar update (0.192 units) and for the endgame.
- mod-r addition (4-limb add, compare, correction): <= 32 instructions.
- mod-r exponentiation with a 255-bit exponent (Montgomery ladder of
  ~381... here ~255 squarings and ~128 multiplications): <= 383 * 384
  instructions = 73.7 units; used at most 3 times per run (window start,
  phase-2 setup, endgame), charged.
A reference operation and its implementation are never double-counted.

## Number-theoretic setup (verified by experiment.py Part A)

- r - 1 = 2^32 * 3 * 11 * 19 * 10181 * 124503 * 859523 * 902509^2 * 2504889
  * 2526331 * 52437067 * 254436773^2. Every prime factor is < 2^28, so the
  attack re-derives the factorization at runtime by trial
  division/Pollard rho in <= 2^21 units (charged); the experiment
  *verifies* the listed factorization (Miller-Rabin primality of every
  factor, product equality) - a polynomial-time check.
- 2^20 | r-1, and d = 2^20 is the largest divisor of r-1 not exceeding the
  scoring n: any divisor with an odd part is at most 2^18 * 3 < 2^20.
- zeta = 7 generates Z_r*; eta = zeta^(2^20) has exact order
  N = (r-1)/2^20 (log2 N = 234.857); zeta_e = zeta^N has exact order 2^20.

## Why the algebra works

s is uniform in [1, r-1) hence in Z_r*. The d-th power map on Z_r* (order
r-1, d | r-1) is exactly d-to-1 onto the unique order-N subgroup, so s^d is
uniform there and s^d = eta^k1 for a unique k1 uniform in [0, N). Likewise
s = zeta^e0 uniquely with e0 = k1 + t*N for a unique t in [0, d). Hence

A_d = [eta^k1]P0,  Q' := [zeta^(-k1)]A_1 = [zeta_e^t]P0,  s = zeta^(k1+t*N).

k1 in [0, N) and t in [0, 2^20) are the two search spaces (Cheon's
reduction).

## Algorithm v2 (single shot; the only coins are the window offset)

0. Read A_1 and A_d (two point reads, charged).
1. Runtime factorization of r-1; construct zeta, eta, zeta_e; build
   fixed-base comb tables for P0, A_d and A_1 (each: [b*2^(8k)]base for
   k < 32, 1 <= b < 2^8, i.e. 8,160 points built by <= 2^14 reference
   additions, ~1 MiB).
2. Phase 1, partial-coverage BSGS for k1 with a contiguous window:
   m1 = ceil(sqrt(N)), M = ceil(N/m1), m = ceil(sqrt(0.40*N))
   = 141430898198098033258261195556338786 (m < m1; here M = m1).
   - Baby table: for i = 0..m-1 compute the scalar eta^i (one mod-r
     multiplication per step, sequential) and the point [eta^i]P0 by one
     fixed-base comb scalar multiplication (<= 32 additions of table
     points); insert (point -> i) into an open-addressing hash table with
     fixed capacity 2^117 slots (>= 2m; load factor < 0.7, expected probe
     count < 2; each insert charged 1000 instructions).
   - Window: draw a uniform offset o in [0, M); the window is the m
     consecutive blocks o, o+1, ..., o+m-1, taken modulo M. If the window
     passes block M-1 it splits into two contiguous runs; run 1 starts at
     block o with scalar eta^(-m1*o) (one exponentiation), and run 2 (if
     any) starts at block 0 with scalar 1 - still exactly one mod-r
     multiplication per giant within each run.
   - Giants: for each block j in the window, compute the running scalar
     (one mod-r multiplication per giant) and the point
     Y = [eta^(-m1*j)]A_d by one comb scalar multiplication (<= 64 units);
     look Y up in the baby table (charged 1000 instructions); on a hit at
     baby index i, set k1 = (i + m1*j) mod N.
3. Phase 2, full-coverage BSGS for t: Q' = [zeta^(-k1)]A_1 (one comb
   multiplication); baby table {[zeta_e^i]P0 : i < 2^10}, giants from Q' by
   multipliers zeta_e^(-2^10); 2^11 comb multiplications, deterministic,
   covering all t in [0, 2^20).
4. Assemble s = zeta^(k1 + t*N) mod r, verify [s]P0 = T = A_1 with one
   comb multiplication, output s. No phase-1 hit -> the run fails; no
   restart or amplification.

Correctness of extraction: a baby/giant collision means
[eta^i]P0 = [eta^(-m1*j)]A_d, i.e. eta^(i + m1*j) = s^d in Z_r*, so
k1 = (i + m1*j) mod N is the true exponent (the reduction mod N is exact
because eta has order N); phase 2 recovers t exactly; the final verification
closes the chain. The toy reproduces this end to end with zero wrong
extractions.

## Success probability (exact integer arithmetic, no heuristic)

k1 is uniform over [0, N). Write k1 = i + m1*j with i = k1 mod m1 and
j = floor(k1/m1). The run covers k1 iff i < m and j lies in the window.
Block j contains b_j = min(m1, N - m1*j) consecutive k1-values, of which
min(b_j, m) have i < m (m < m1, so all full blocks contribute exactly m).
Every window of m consecutive blocks therefore covers between
(m-1)*m + min(b_last, m) and m*m positions, where b_last is the size of the
single possibly-short block (at most one such block exists, and M = m1
blocks exist in total). Hence the success probability for ANY offset is at
least (m-1)*m / N = m^2/N - m/N = 0.40 - 2^-118 >= 0.39, and the
expectation over the uniform offset is computed exactly in Part A as

  E = (sum over blocks of min(b_j, m)) * m / (M * N) = 2/5 exactly.

There is no walk model, no mixing assumption and no probabilistic tail in
this argument; the only randomness is the uniform challenge and the uniform
window offset. Phase 2 is deterministic given a phase-1 hit. Claimed lower
bound: **0.39**.

## Work ledger (all charged)

Per fixed-base comb scalar multiplication of a 255-bit scalar: <= 32
additions of comb-table points (batched reference additions - the same
reference operation, batched across the concurrent multiplications) plus
masking/loads/control: charged **64 units** (2x the analytic 32.1).

| Component | Count | Units |
| --- | ---: | ---: |
| Phase-1 baby points (comb mults + one mod-r mult per scalar) | m | m*(64 + 0.192) |
| Baby-table inserts (open addressing, < 2 probes, 1000 instr) | m | m*0.5 |
| Phase-1 giant points (comb mults + one mod-r mult per scalar) | m | m*(64 + 0.192) |
| Giant lookups (1000 instructions each) | m | m*0.5 |
| Window-start exponentiations (<= 2 per run) | <= 2 | <= 148 |
| Phase-2 BSGS (comb mults) | 2^11 | 2^17 |
| Endgame: [zeta^(-k1)]A_1, verification, output | <= 4 mults + 1 exponentiation | <= 300 |
| Preprocessing: factor r-1, zeta/eta/zeta_e, three comb tables | one-off | <= 2^21 |

With m = 2^116.77:

W(2^20) = 2m*(64 + 0.192) + 2m*0.5 + 2^17 + 2^21 + 450
        = 2^123.777511 (rounded up; exact rational 123777511/1000000).

All secondary terms are below 2^-7 of the leading term. Preprocessing
<= 2^21 units -> preprocessing_log2 = **21.0**, included in the total.

## Memory

Baby table: fixed capacity 2^117 slots x 128 bytes (point 96, index 8,
chaining/metadata 24) = 2^124 bytes exactly; the m = 2^116.77 inserts never
exceed capacity (deterministic). Comb tables ~3 MiB; scalar and walker
scratch negligible. Peak storage 2^124.0 bytes -> memory_log2_bytes =
**124.0**. Nonuniform advice: **0 bytes** (everything is constructed and
charged at runtime).

## Alternatives and sensitivity

- The r+1 variant needs d' | r+1 with d' <= 2^20; the largest such divisor
  is 484486 < 2^20, giving a larger phase-1 space plus Fp2 arithmetic in
  its second phase. Rejected.
- d = 2^21 divides r-1 (2^32 || r-1) but A_{2^21} is not supplied at the
  scoring n = 2^20.
- Using more supplied powers via polynomial combinations [p(s)]P0 gains
  nothing: p(s) is uniform over Z_r for nonconstant p, so a birthday search
  over such values needs ~sqrt(r) samples regardless of n; only the
  subgroup structure of s^d shrinks the space.
- Sensitivity: charging a full batched-affine double-and-add (384 units)
  for every scalar multiplication gives 2^126.3 - the claim does not depend
  on the comb.

## Evidence and reproduction

- Part A (production parameters, exact): 2^20 | r-1; the factorization of
  r-1 verified (Miller-Rabin + product equality, all factors < 2^28);
  zeta = 7 generates Z_r*; eta = zeta^(2^20) has exact order N;
  zeta_e = zeta^N has exact order 2^20; the exact coverage rational is
  computed and printed (2/5) and asserted >= 0.39.
- Part B (toy, complete two-phase attack, d = 6 | n-1 on a 20-bit
  prime-order group, N_toy = 174563, contiguous window with sequential
  scalars exactly as in production): declared run
  `{"trials": 200, "mode": "partial"}`: 88/200 = 0.44 successes vs the
  exact predicted coverage 14675965/36483667 = 0.40226; 0 wrong
  extractions. Full-coverage control mode: 200/200 with 0 wrong
  extractions.
- Reproduce:
  `bash scripts/run_experiments.sh bls12-381-g1-ecdlp submissions/bls12-381-g1-ecdlp/partial-cheon-v2`
  or directly:
  `echo '{"trials": 200, "mode": "partial"}' | python3 submissions/bls12-381-g1-ecdlp/partial-cheon-v2/experiment.py`.

Small experiments are evidence only at their actual scope.

## Heuristics, limitations, objections

- The success probability uses no heuristic: exact integer combinatorics
  (uniformity of k1 from the d-to-1 power map; contiguous-window coverage
  with worst case 0.4 - 2^-118 for every offset).
- The soft numbers are instruction sub-costs: the 64-unit comb charge (2x
  analytic), the 384-instruction mod-r multiplication (explicit
  schoolbook+Montgomery construction), and 1000-instruction hash
  operations. The sensitivity row shows the claim survives replacing the
  comb by double-and-add everywhere (2^126.3).
- Single shot at the 0.39 threshold: deterministic full coverage (m = m1
  and window = all blocks) costs about 2^124.64 (+0.86 log2).
- Memory 2^124.0 bytes is the sqrt(N)-sized baby table; reported, no score
  effect.
- Standard cryptanalysis (Cheon's attack); the contribution is the
  sequential-window structure, the exact coverage computation, and honest
  fixed-base accounting under this protected cost model.

References: `tracks/docs/COST_MODEL.md`; J. H. Cheon, "Discrete logarithm
problems with auxiliary inputs" (ASIACRYPT 2008 / J. Cryptology 2010); the
track's protected rules.md (d = 2^20 | r-1, Cheon-type methods admissible);
HAC ch. 3.6 (BSGS). Self-contained above.
