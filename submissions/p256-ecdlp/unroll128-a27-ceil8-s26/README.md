# unroll128-table-sched-finvw-a27 (128-copy jump-table blocks)

P-256 sigma rho, declared **126.73132979** (12673132979/100000000 at C=1000; the same
exact ledger total W, log2(W) = 126.731329784942..., now rounded up at the eighth decimal
instead of the sixth, 5.06e-9 above log2(W), certified by 60-digit Decimal logarithms in
`accounting.py`), on the anchor-skip A=2^27 kernel (126.732475). Both batched blocks hold
128 labelled copies; the backward block enters by a 128-entry signed-offset jump table and
the forward block by a fixed branch into copy 2. Each end test runs once per 128 lanes.
Measured on all three arms at B=2^16 and B=2^20, this subtracts 491,517
instructions per 2^20-walk round (15B/32 - 3), with a local trampoline keeping the 2^60 cap
test in range. The kernel also keeps the windowed `finv` (11,814 per round), the
one-instruction cap, the rare-path repricing and the three-arm anchor schedule.

PLAIN_FIXED 605,090,025; C_COMMON 576.108991; CORR 14.210693; reserve 0.0075889;
PER_STEP 590.327272. Preprocessing 2^84.961190 units, peak memory below 2^77 bytes.
Paired counts and identical walk dumps are in `evidence/unroll128.json`; the schedule
rows are in `evidence/anchor-schedule.json`; the argument is in `proof.md` and
`unroll.md`. The base submissions carry `branchy.json`, `forced/summary.json`,
`sign-law.json` and `gen55.py`; this package omits them for the 256 KiB limit, and also
omits the redundant `claim.json` metadata duplicate (same content in `submission.json`
and `proof.md`). See `accounting.py` and `submission.json`.
