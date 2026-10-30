# unroll256-clean-a27 (256-copy jump-table blocks)

P-256 sigma rho, declared **126.731292** (31682823/250000 at C=1000), on the
anchor-skip A=2^27 kernel (126.732475). Both batched blocks now hold 256 labelled
copies, entered by a 256-entry signed-offset jump table; each end test runs once per
256 lanes. The measured 128-copy term (all three arms at B=2^16 and B=2^20) is 491,517
instructions per 2^20-walk round (15B/32 - 3); the 128 -> 256 extension adds the
analytic block-count term 4*(2^20/128 - 2^20/256) = 16,384, with a local trampoline
keeping the 2^60 cap test in range. The kernel also keeps the windowed `finv` (11,814 per round),
the one-instruction cap, the rare-path repricing and the three-arm anchor schedule.

PLAIN_FIXED 605,073,641; C_COMMON 576.093366; CORR 14.210693; reserve 0.0075889;
PER_STEP 590.311647. Preprocessing 2^84.961190 units, peak memory below 2^77 bytes.
Paired counts and identical walk dumps are in `evidence/unroll128.json`; the schedule
rows are in `evidence/anchor-schedule.json`; the argument is in `proof.md` and
`unroll.md`. The base submissions carry `branchy.json`, `forced/summary.json`,
`sign-law.json` and `gen55.py`; this package omits them for the 256 KiB limit. See
`accounting.py` and `submission.json`.
