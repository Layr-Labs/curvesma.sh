# unroll256-table-sched-finvw-a27 (256-copy jump-table blocks)

P-256 sigma rho, declared **126.731292** (31682823/250000 at C=1000), on the
anchor-skip A=2^27 kernel. Both batched blocks hold 256 labelled copies, entered by a
256-entry signed-offset jump table; each end test runs once per 256 lanes. The measured
128-copy term at B=2^16 and B=2^20 subtracts 491,517 instructions per 2^20-walk round
(4B/8 - 4B/128 - 3). Halving the block count adds the analytic term
4*(B/128 - B/256) = 16,384, giving 507,901 per round (4B/8 - 4B/256 - 3); the six-instruction
dispatch is unchanged and the local trampoline keeps the 2^60 cap test in range. The kernel
also keeps the windowed `finv`, the one-instruction cap, the rare-path repricing and the
three-arm anchor schedule.

PLAIN_FIXED 605,073,641; C_COMMON 576.093366; CORR 14.210693; reserve 0.0075889;
PER_STEP 590.311647. Preprocessing 2^84.961190 units, peak memory below 2^77 bytes.
Measured 128-copy rows are in `evidence/unroll128.json`; the analytic 128 -> 256 extension
is in `evidence/unroll256.json`; schedule rows in `evidence/anchor-schedule.json`; argument
in `proof.md` and `unroll.md`. Base submissions carry `branchy.json`, `forced/summary.json`,
`sign-law.json` and `gen55.py`; this package omits them for the size limit.
