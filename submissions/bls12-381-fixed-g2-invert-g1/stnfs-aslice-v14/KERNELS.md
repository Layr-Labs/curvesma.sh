# Scope of the retained and new kernel measurements

v14 keeps the inherited bounded v12/v13 arithmetic/walk measurements
(measurements12.txt: callgrind Ir, Debian12/linux-arm64, GCC12.2, valgrind 3.19,
flags -O3 -fno-tree-vectorize -fno-tree-loop-distribute-patterns, differential
two-repetition method) for the main-phase slice walk, certified norm bounds, ECM
and the la_mac.s loop, and adds the descent-sieve kernel
(descent_sieve.c + measurements14.txt): the 12-dimensional hit-lattice walk,
record append, apply and scan of the descent nodes. The v13 production record
allocator and its fourteen-level routing are withdrawn; STORAGE.md specifies the
lean line-local protocol with an executable reference (store14.py) and an
analytical allowance. No claimed production bound is inferred solely from a
small correctness check.

## Descent sieve (new in v14)

descent_sieve.c implements, per sieving prime p and root-pair functional
w (12 coefficients mod p): the systematic basis of the hit lattice
L_p = {c : w.c = 0 mod p} (determinant p), a 12-dimensional integer LLL
reduction (float Gram-Schmidt for decisions, exact __int128 row updates), the
Schnorr-Euchner enumeration of L_p's points in the box c in [-k,k]^12, the
mixed-radix box index of each hit, the append of (segment = idx>>16, offset,
scaled log2 p) records into per-segment bounded buckets, and the per-segment
apply (sum array of 2^16 entries) and threshold scan. Its `check` subcommand
verifies against brute force that the walk enumerates exactly the congruence
solutions (no duplicates, correct count within 5% at k>=8, exact at k<=1), and
its `bench` subcommands are measured at p near 2^20, 2^40 and 2^62 and
k = 16..128 (three trials, max charged). The measured instruction counts are in measurements14.txt. The model charges
PAIR_SETUP = 2^24.8 per (prime, root pair) for the measured sound reduction plus
the functional (the double-decision build at 2^18.65 is unsound at 62-bit and is
not charged), and SIEVE_CELL = 2^21.64 per cell: the measured 2946 Ir per leaf
times the smallest measured leaves-per-hit (254 at p = 2^20; range 254..1663),
4.32 hits per cell, plus 10.9 zero/apply/scan and 1200 per cell for the c-space
certified screen (twice the measured a-space evaluation cost at twice its rate:
the interval-(a,b) variant is specified in DESCENT.md, not separately measured).
The production leaves-per-hit extrapolation and both-end sensitivities are in
DESCENT.md; beyond the measured p sizes the reduction grows about 10% per log2 p.

## Slice walk (retained)

slice12.c: bounded packed-row benchmark, line and record paths, exact recentring,
Barrett setups, apply/scan/recovery; the reallocating test driver is excluded from
production. (Its historical differential-check harness slice12_check.c and the
raw logs stay in the v12/v13 packages; the retained measurements and the digest
references are in measurements12.txt.) Line costs by hits per line coset:
(13.82,24.10),(9.82,31.11),(5.82,88.31),(1.82,310.70); record-path costs by hits
per whole-slice coset: (18.91,44.90),(10.91,78.21),(3.91,266.97),(0.91,938.08),
(-1.09,2231.66); empty-coset allowance 1600; apply 7.00 and recovery 7.003 per
record (before the STORAGE.md routing supplement of 36 per record); layout
12.325 per cell, zero 0.625, scan 8 per cell plus 21 per leaf. These are the
constants in model14.py.

## Certified norm bounds (retained; extended to c-space in v14)

The C port nbound12.c (exact rational Rouche disks, outward one-ulp binary64
intervals, dyadic log table, recursive splits to width 128) is retained in the
v12/v13 packages; its measured constants (5900 instructions per evaluation,
4200 per row, 15300 per slice, mean 0.043 evaluations per cell) are in
measurements12.txt, and the executable Python reference embedding_bounds.py is
shipped here. v14 uses the same machinery in the descent's c-space along the
shortest reduced-basis direction, where a and b are affine in the line
coordinate with integer coefficients; the sandbox audit [11] checks actual
special-q norms against the envelopes, and the block-selection study
re-measured 2.4-3.3 evaluations per 64-cell block over all 2401 pool rows.

## ECM and linear algebra (retained)

ecm12.c: five-limb Montgomery multiplication, checked on 20,000 products,
367.76 instructions; ECM curve budgets and limb scaling remain analytical
models. la_mac.s: 33-instruction inner loop with 64-bit column indices and
five-limb accumulator ([10] recounts it); the v14 price of 88 per entry covers
the loop plus block selection, signs, four-hop partial transfers over 16 blocks
and Barrett finishing (STORAGE.md); counting the loop alone would not establish
that price.

The shared verifier's sandbox reads these C/assembly sources but does not
compile them. All executable evidence runs under the credential-free
llm-verifier runner, never a candidate process on the host.
