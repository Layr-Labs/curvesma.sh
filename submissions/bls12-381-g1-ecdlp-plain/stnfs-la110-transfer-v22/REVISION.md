# Revisions to the GT STNFS core (source: e3940d11, core ledger 2^121.99438)

Polynomials, factor base, pool, filtering and H1-H5 are unchanged.

## R1-R9 (as in our AI-valid 3eccc6c9, 124.9584)

ECM_FIND 0.999 (403 curves); chunk-scatter route, measured 10.26 charged 10.3;
R3 replaced by R10; R5 withdrawn (re-read 7.003); [2^63, 2^64) record walk 2x;
line sieve from 2^11; R8 yield guard; R9 plan + LOG2_RELS 63.97, theta 190/500,
slacks (297, 233), allowance 1.0676; collection <= 121.1117, two-log 2^118.2732,
meter 2^118.996, success >= 0.390320.

## R10 (from our recover-Q 711418d3, AI-valid): fully priced linear algebra, pad 1.10

Loop-free kernels (row tail 2,750/row, red9 420, 21/entry, H4u); pad 1.61 ->
**1.10** on the measured N/U trend (H3). LA 2^118.7428 -> **2^117.7327**,
generator 2^105.5088 -> 2^107.2082; ledger <= 121.6221 -> **<= 121.5198**.

## R11 (submitted as aa2a2ba5, failed on packaging oversize only): exact price ceilings, 0.01 table shift, reserve 1001/1000

The nine line/record prices become the exact logged quotients (audit
`exact_prices` re-parses the Ir pairs; setup6/norm-bound verified ceilings).
Both tables are checked decreasing and convex (exact slopes; `price_shape`),
and `tab()` charges 0.01 below its argument, covering boundary rounding
twofold; production log2h is checked covered (`tab_shift_ranges`). Collection
at +half-ulp on every quoted decimal moves the ledger 1.000036x
(`residual_bound`), so the reserve drops 101/100 -> **1001/1000**. Collection
<= 121.1117 -> **<= 121.1127**; ledger <= 121.5198 -> **<= 121.5206**; G1
124.85609 -> **124.843970069062**.

## R12 (submitted as 494e4872, rejected: unpriced line-factor recovery)

The line/record split moved P_MID 20 -> 48.75, but line-sieved primes in
[2^20, 2^48.75) were neither trial-divided nor recorded, leaving their
cofactor removal unpriced. R12b attempted the revision; R12c below completes it.

## R12b (submitted as 43f5ea7e, rejected: step-p recovery unjustified; superseded by R12c)

Priced a second line pass with a step-p flag loop (6/multiple) plus setups,
P_MID 47.5. The review found hits are not step-p in the packed cell index, the
consecutive-integer toy did not test the packed layout, and row-wise or
rerun alternatives were unpriced. Claimed collection <= 120.8569, ledger <=
121.3321, G1 124.655470069062 (withdrawn).

## R13 (submitted as 16253187, validation cancelled: it inherited the R12b gap; content re-verified here)

SLACK 96 -> **60** in T (cutoff s follows; meter pinned). Tighter T cuts
node work faster than retention loss rises (ret min 0.8575; base mu 401/441;
node phi <= 0.0386, psi <= 0.8): two-log 2^118.2732 -> **2^118.2303**. Meter
re-solved (SLACK scan flat at 58-65): cap **2^118.9532**, Markov **0.60590**,
boot phi **7.858e-04**, E[N] 2^18.005, success claim **0.390320**.

## R12c (as in our AI-valid 0f2da92, 124.754970069062): recovery reruns the coset enumerator, P_MID 45.5

The second pass reruns the first-pass coset enumerator over every (line, p,
pair) and flag-checks each hit (STORAGE.md): identical stepping code, so
stepping costs at most the measured update `tab(h)` (the skipped score path
is >= 0); flag ops LDRB+CBZ = 2, recounted (`line_recov_loop`); emit is
<= 2^97 total (absorbed); setups recompute the same coset structure (488.24
per (line, p, pair); K counts lines because setups are per-line objects).
Completeness: slice12's own line-path-vs-brute-force checks (digests MATCH)
plus a packed-layout lattice-enumeration toy with root-pair congruences
(`line_recov_toy`: recorded = brute-force divisors, no duplicates). P_MID
re-optimised on the 0.25 grid: **45.5** (neighbours worse; work SE 3.406%,
allowance 1.06905 holds). Collection <= 121.1127 -> **<= 121.0039** (387.8
per cell); ledger <= 121.5206 -> **<= 121.4316** (unrounded 121.4315135, rows
121.4315234); G1 124.843970069062 -> **124.754970069062**. Yield, descent,
success and all other parameters are untouched.

## R14 (as in our submitted v20tdiv, 124.697570069062): trial division to 2^30, P_TDIV 20 -> 30

Rerun stepping cost 40.3 per cell while trial division is per survivor
(3.2e-11 of cells), so the recovery boundary moves up: trial division covers
every prime below 2^30 on both sides (512 per prime per side, at most
64,806,364 primes; C_TDIV = 2^35.95 per survivor, STORAGE.md) and the R12c
rerun covers [2^30, 2^45.5). Completeness and yield are unchanged.
Work SE 3.406% -> 3.439%; 1.06905 still covers 1 + 2.00 SE (needs 1.06877).
Collection <= 121.0039 -> **<= 120.9262**; ledger <= 121.4316 -> **<= 121.3742**
(unrounded 121.3741113, rows 121.3741570); G1 -> **124.697570069062**.

## R15 (as in our submitted v21, 124.687770069062): the first pass records line hits at/above 2^30; the rerun is deleted

Per-hit recovery by re-stepping (tab(h) + 2 flag ops, plus 1.0 flag stores
per cell and per-(line, p, pair) setups) costs more than the shared record
path (emit 7 inside REC_ROUTE 10.3, apply 7.00, recovery 7.003), so the
second pass is gone: the first-pass line enumerator emits a 16-byte record
per hit with p in [2^30, 2^45.5) through the identical chunk-scatter path,
and chunk owners apply, scan and recover line records exactly like walk
records (STORAGE.md). The per-cell record rate becomes NBAR ln(2^64/2^30);
no new prices, no new heuristic. Completeness: the same complete
enumerator slice12 checks against brute force, plus an unconditional store
per enumerated hit. Yield, descent, success and every other parameter are
untouched. Work SE 3.439% -> 3.444% (declared bound 0.03444); 1.06905
still covers 1 + 2.00 SE (needs 1.06887). Collection <= 120.9262 ->
**<= 120.9128** (364.1 per cell); ledger <= 121.3742 -> **<= 121.3644**
(unrounded 121.3643154, rows 121.3643460); G1 -> **124.687770069062**.

## Evidence

`audit_g1.py` recomputes every changed quantity (`la_prices`, `matrix_pad`,
`exact_prices`, `price_shape`, `tab_shift_ranges`, `residual_bound`,
`residual_ledger`, `no_rerun`, `tdiv_prime_bound`, the BOOT16 pin;
proof.md, Evidence); `evidence16.py` [10] recounts the kernels.

## R16 (this package): descent boot + plan re-searched at slack 60, meter re-solved; stated residual multipliers corrected

The shipped boot/plan was coordinate-searched before R13 moved SLACK 96 ->
60. Re-searching at 60 under the same model and constraints (mu >= 16, phi
<= 0.05, psi <= 0.8, iA on the x8 mesh lattice, R in RADII) moves BOOT16
(215,344,15,4) -> **(214,320,16,4)** (nA 189.83 hugs the 190 bin top, so the
16 A-children stay in cheaper bins) and twelve plan entries (85f, 85g, 95f,
110f, 110g, 115f, 120g, 125f, 125g, 165f, 165g, 185g). No model, heuristic,
price or collection parameter moves. Expected two-log work 2^118.2303 ->
**2^117.5873** (one log 2^130.875; boot psi 0.6732, phi 1.780e-03; min mu
32.3; E[N] 2^16.87, P(abort) 2^-10.42); cap re-solved to **2^118.3128**
(Markov 0.60478), re-binding success at the 0.390320 claim (model 0.3903323).
Ledger <= 121.3644 -> **<= 121.2637** (unrounded 121.2636532, rows
121.2636860); G1 -> **124.587070069062**.

Also fixed: v21 was rejected only for stating the v19 residual figure
(1.000014x) while the computed collection multiplier is 1.0000286x (ledger
1.0000224x); the verdict raised no objection to any model number or
mechanism. This package states the computed figures (collection **1.000029x**,
ledger **1.000023x**) and the audit now checks the stated multipliers, not
just 1.0001.
