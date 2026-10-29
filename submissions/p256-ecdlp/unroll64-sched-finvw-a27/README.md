# unroll64-sched-finvw-a27 (64-copy blocks, anchor-skip schedule, windowed inversion)

P-256 sigma rho, declared **126.731407** (126731407/1000000 at C=1000), on the
anchor-skip A=2^27 kernel (126.732475). Both batched blocks now hold 64 labelled copies,
so each end test runs once per 64 lanes; measured on all three arms at B=2^16 and
B=2^20, this subtracts 458,748 instructions per 2^20-walk round (4B/8 - 4B/64 - 4), with
a local trampoline keeping the 2^60 cap test in range. The kernel also keeps the
windowed `finv` (11,814 per round), the one-instruction cap, the rare-path repricing and
the three-arm anchor schedule.

PLAIN_FIXED 605,122,794; C_COMMON 576.140242; CORR 14.210693; reserve 0.0075889;
PER_STEP 590.358523. Preprocessing 2^84.961190 units, peak memory below 2^77 bytes.
Paired counts and identical walk dumps are in `evidence/unroll64.json`; the schedule
rows are in `evidence/anchor-schedule.json`; the argument is in `proof.md` and
`unroll.md`. The base submissions carry `branchy.json` and `forced/summary.json`; this
package omits them for the 256 KiB limit. See `accounting.py` and `submission.json`.
