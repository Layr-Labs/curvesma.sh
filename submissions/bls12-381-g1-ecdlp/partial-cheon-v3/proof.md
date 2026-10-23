# bls12-381-g1-ecdlp: partial-coverage Cheon attack v3, at n=2^20

Status: ready for independent exploratory review. Classical single-shot
attack using the supplied auxiliary power A_d, d = 2^20; no quantum steps,
no uncharged advice, no amortization over targets. Cheon's
discrete-logarithm-with-auxiliary-input attack is standard cryptanalysis.
v3 responds to the second-round review of v2 (the review accepted the
reduction, the contiguous-window construction and the sequential giant
updates, and rejected the accounting):

1. **The ledger now matches the score.** v2's own ledger table summed to
   129.384m while its script computed 128.884m; the review's arithmetic
   (log2(129.384m) = 123.783096 before extras) was correct and exceeded
   the submitted 123.777499. v3 recomputes the total from the table exactly
   and submits that value: **123.783097**.
2. **No "exactly 2/5" claim.** The exact coverage rational is extremely
   close to, but not equal to, 2/5 (v2's prose wrongly asserted equality;
   only the float display was 0.4). v3 reports the exact rational as
   printed by the experiment and an exact integer verification that it is
   >= 39/100, plus the worst-single-offset bound (m-1)m/N = 0.4 - 2^-118.
   m is now computed by pure integer arithmetic
   (m = isqrt(floor(2N/5)) + 1), removing the floating-point rounding the
   review identified in v2's int(0.40*N) expression.
3. **Hash-table parameters stated honestly.** With m = 2^116.7676 inserts,
   a 2^117-slot table has load 0.851 (v2 wrongly asserted < 0.7 and
   "capacity >= 2m"). v3 allocates 2^118 slots: load 0.4256, expected
   unsuccessful-search probes 1/(1-0.4256) = 1.74 < 2, so the charged
   1000-instruction insert/lookup covers the probe sequence with margin.
4. **Memory and index width.** v2 claimed exactly 2^124 bytes while
   listing positive additional storage, and an 8-byte index that cannot
   represent 117-bit baby indices. v3: slots are 128 bytes with a 16-byte
   index; the table is exactly 2^125 bytes and the total peak is strictly
   above that, reported as **125.1** rounded up.

## Claimed bounds

| Field | Value |
| --- | --- |
| time_log2 (score) = log2(W(2^20)) | **123.783097** (= 123783097/1000000, rounded up) |
| success_probability (supported lower bound) | 0.39 (exact rational expectation printed; >= 39/100 by integer verification) |
| memory_log2_bytes | 125.1 |
| preprocessing_log2 (included in total) | 21.0 |
| nonuniform_advice_log2_bytes | 0 |

W(n) for n >= 2^20 is constant: the attack reads only A_1 = T and
A_d = [s^d]P0 with d = 2^20 <= n (at the scoring n = 2^20, A_d is the last
element of the input array, index d-1). The bound does not apply for
n < 2^20.

## Conventions

Protected target: BLS12-381 G1, generator P0 of prime order r, challenge
s uniform in [1, r-1), T = A_1 = [s]P0, powers A_i = [s^i mod r]P0 for
i = 1..n, scored at n = 2^20. Cost model `classical-group-work-v2`:
cost = H + R/2000; H = reference operations (batched incomplete affine G1
addition, 6M, one reference operation = 1 unit), R = other scalar 64-bit
AArch64 instruction equivalents. Pinned base-field multiplication: 383
instructions (blst `blst_fp_mul`, mod p). Mod-r prices (r 255 bits) are
analytic constructions: multiplication <= 384 instructions (5x52-limb
schoolbook + Montgomery reduction), used once per baby/giant scalar update
(0.192 units) and at most three times per run as exponentiations
(<= 383*384 instructions each, charged); addition <= 32 instructions. A
reference operation and its implementation are never double-counted.

## Number-theoretic setup (verified by experiment.py Part A)

- r - 1 = 2^32 * 3 * 11 * 19 * 10181 * 124503 * 859523 * 902509^2 * 2504889
  * 2526331 * 52437067 * 254436773^2; every prime factor is < 2^28, so the
  runtime re-derivation (trial division/Pollard rho, <= 2^21 units, charged)
  is cheap; the experiment verifies the listed factorization (Miller-Rabin
  primality of every factor plus product equality).
- 2^20 | r-1, and d = 2^20 is the largest divisor of r-1 not exceeding the
  scoring n (any divisor with an odd part is at most 2^18 * 3 < 2^20).
- zeta = 7 generates Z_r*; eta = zeta^(2^20) has exact order
  N = (r-1)/2^20 (log2 N = 234.857); zeta_e = zeta^N has exact order 2^20.

## Why the algebra works

s uniform in [1, r-1) lies in Z_r*. The d-th power map on Z_r* (order
r-1, d | r-1) is exactly d-to-1 onto the unique order-N subgroup, so s^d is
uniform there and s^d = eta^k1 for unique k1 uniform in [0, N); and
s = zeta^(k1 + t*N) for unique t in [0, 2^20). Hence
A_d = [eta^k1]P0, Q' = [zeta^(-k1)]A_1 = [zeta_e^t]P0.

## Algorithm v3 (single shot; the only coin is the window offset)

0. Read A_1 and A_d (charged).
1. Runtime factorization of r-1; construct zeta, eta, zeta_e; build
   fixed-base comb tables for P0, A_d, A_1 (each 8,160 points,
   <= 2^14 reference additions, ~1 MiB).
2. Phase 1, partial-coverage BSGS with a contiguous window:
   m1 = ceil(sqrt(N)), M = ceil(N/m1) (= m1 here), and
   m = isqrt(floor(2N/5)) + 1 = 141430898198098033869460368287841487
   (pure integer arithmetic; m < m1).
   - Baby table: for i = 0..m-1: scalar eta^i (one mod-r multiplication
     per step, sequential), point [eta^i]P0 (one comb multiplication),
     insert (point -> i) into an open-addressing table with fixed
     capacity 2^118 slots (load 0.4256, expected probes < 2, charged
     1000 instructions per insert). Slot: 128 bytes = point 96 + index 16
     + chaining/metadata 16.
   - Window: uniform offset o in [0, M); the m consecutive blocks
     o, ..., o+m-1 (mod M) are processed as at most two contiguous runs
     (run 2 starts at block 0 with scalar 1); each giant scalar is one
     sequential mod-r multiplication; at most two window-start
     exponentiations per run (charged).
   - Giants: point [eta^(-m1*j)]A_d (one comb multiplication) per block j
     in the window; lookup in the baby table (1000 instructions); on a hit
     at index i: k1 = (i + m1*j) mod N.
3. Phase 2 (deterministic, full coverage): Q' = [zeta^(-k1)]A_1; BSGS for
   t in [0, 2^20) with 2^11 comb multiplications.
4. s = zeta^(k1 + t*N) mod r; verify [s]P0 = T; output. No phase-1 hit ->
   fail; no restart or amplification.

Correctness: a baby/giant collision means eta^(i + m1*j) = s^d in Z_r*, so
k1 = (i + m1*j) mod N is exact (eta has order N); phase 2 recovers t
exactly; verification closes the chain. The toy reproduces this with zero
wrong extractions.

## Success probability (exact integer arithmetic, no heuristic)

k1 is uniform over [0, N). Write k1 = i + m1*j, i = k1 mod m1,
j = floor(k1/m1). The run covers k1 iff i < m and j is in the window.
Block j holds b_j = min(m1, N - m1*j) consecutive values, of which
min(b_j, m) have i < m (m < m1). Every window of m consecutive blocks
covers at least (m-1)*m + min(b_last, m) >= (m-1)m values and at most
m^2, where b_last is the single possibly-short block; so for EVERY offset
the success probability is at least (m-1)m/N >= 0.4 - 2^-118 >= 0.39.
The expectation over the uniform offset is the exact rational
E = (sum_j min(b_j, m)) * m / (M*N), printed by the experiment and
verified >= 39/100 by exact integer comparison. No walk model and no
probabilistic tail is used; the only randomness is the uniform challenge
and the uniform offset. Claimed lower bound: **0.39**.

## Work ledger (all charged)

Per comb scalar multiplication of a 255-bit scalar: <= 32 additions of
comb-table points (the same batched reference operation) plus
masking/loads/control: charged **64 units** (2x the analytic 32.1).

| Component | Count | Units |
| --- | ---: | ---: |
| Baby points: comb mult + sequential scalar mult | m | m*(64 + 0.192) |
| Baby-table inserts (open addressing, < 2 probes, 1000 instr) | m | m*0.5 |
| Giant points: comb mult + sequential scalar mult | m | m*(64 + 0.192) |
| Giant lookups (1000 instructions) | m | m*0.5 |
| Window-start exponentiations (<= 2 per run) | <= 2 | <= 148 |
| Phase-2 BSGS | 2^11 | 2^17 |
| Endgame: [zeta^(-k1)]A_1, verification, output | <= 4 mults + 1 exponentiation | <= 300 |
| Preprocessing: factor r-1, zeta/eta/zeta_e, three comb tables | one-off | <= 2^21 |

W(2^20) = 2m*(64 + 0.192) + 2m*0.5 + 2^17 + 2^21 + 450
        = 129.384m + <= 2^21.9
        = 2^123.783097 (rounded up; exact rational 123783097/1000000).

Preprocessing <= 2^21 units -> preprocessing_log2 = **21.0**, included.

## Memory

Baby table: fixed capacity 2^118 slots x 128 bytes (point 96, 16-byte
index, chaining/metadata 16) = 2^125 bytes exactly; the m = 2^116.77
inserts never exceed capacity (deterministic). Comb tables ~3 MiB; scalar
and scratch state negligible. The total peak is strictly above 2^125 and
far below 2^126: memory_log2_bytes = **125.1** (rounded up). Nonuniform
advice: **0 bytes** (everything constructed and charged at runtime).

## Alternatives and sensitivity

- r+1 variant: largest divisor of r+1 <= 2^20 is 484486 < 2^20 -> larger
  phase-1 space plus Fp2 arithmetic. Rejected.
- d = 2^21 divides r-1 but A_{2^21} is not supplied at n = 2^20.
- Polynomial combinations of the supplied powers: p(s) uniform over Z_r
  for nonconstant p; birthday over them needs ~sqrt(r) samples. No gain.
- Sensitivity: charging full double-and-add (384 units) per multiplication
  everywhere gives 2^126.3; the claim does not depend on the comb.

## Evidence and reproduction

- Part A (production, exact): 2^20 | r-1; factorization verified
  (Miller-Rabin + product equality, all factors < 2^28); zeta = 7
  generates Z_r*; eta has exact order N; zeta_e has exact order 2^20; the
  exact coverage rational is printed (it is 0.4 to 23+ decimal places but
  NOT asserted to be exactly 2/5) and verified >= 39/100 in integers.
- Part B (toy, complete two-phase attack, d = 6 | n-1, N_toy = 174563,
  production-shaped sequential window giants): declared run
  `{"trials": 200, "mode": "partial"}`: 88/200 = 0.44 vs the exact
  predicted coverage 14675965/36483667 = 0.40226; 0 wrong extractions.
  Full-coverage control: 200/200, 0 wrong extractions.
- Reproduce:
  `echo '{"trials": 200, "mode": "partial"}' | python3 submissions/bls12-381-g1-ecdlp/partial-cheon-v3/experiment.py`.

Small experiments are evidence only at their actual scope.

## Heuristics, limitations, objections

- No heuristic in the success argument: exact integer combinatorics with
  the worst-single-offset bound 0.4 - 2^-118.
- Soft numbers are instruction sub-costs (64-unit comb charge at 2x
  analytic; 384-instruction mod-r multiplication by explicit
  construction; 1000-instruction hash operations at load 0.4256); the
  sensitivity row shows the claim survives 384-unit multiplications.
- Single shot at the 0.39 threshold; deterministic full coverage costs
  about 2^124.64 (+0.86 log2).
- The bound applies only for n >= 2^20.
- Standard cryptanalysis (Cheon); the contribution is the sequential
  window structure, the exact integer coverage computation, and honest
  fixed-base accounting under this protected cost model.

References: `tracks/docs/COST_MODEL.md`; J. H. Cheon, "Discrete logarithm
problems with auxiliary inputs" (ASIACRYPT 2008 / J. Cryptology 2010); the
track's protected rules.md; HAC ch. 3.6 (BSGS). Self-contained above.