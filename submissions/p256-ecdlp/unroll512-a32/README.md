# unroll512-a32 (A=2^32 retune of the 512-copy jump-table blocks)

P-256 sigma rho, declared **126.731244** (31682811/250000 at C=1000), on the
anchor-skip kernel retuned from A=2^27 to A=2^32 (126.732475 baseline). The
backward block holds 512 labelled copies entered through a 512-entry signed-offset
jump table (six dispatch instructions once per round); the forward block is entered
by a direct branch to copy 2 and runs B-2 steps. Each end test runs once per 512
lanes. The measured 128-copy term (all three arms at B=2^16 and B=2^20) is 491,517
instructions per 2^20-walk round (15B/32 - 3); the 128 -> 512 extension adds the
analytic block-count term 4*(2^20/128 - 2^20/512) = 24,576, with the table and the
2^60 cap test both inside the local branch range (see unroll.md). The kernel also
keeps the windowed `finv` (11,814 per round), the one-instruction cap, the
rare-path repricing and the three-arm anchor schedule.

PLAIN_FIXED 605,065,449; C_COMMON 576.085553; CORR 14.210693; reserve 0.000292113;
PER_STEP 590.296538. Preprocessing 2^89.961190 units, peak memory below 2^77 bytes.
The A=2^32 re-derivation is in proof.md section "The adder table at A=2^32":
two-cycle rate 2^-33, epsilon 4.563481e-8, slot shift 32 (M=2^164,
DELTA_HASH < 2.9463e-8), DELTA_EXCEPTIONAL 2^-95, DELTA_ADDER_CLASS 1.7748e-29;
accounting.py recomputes all of them (the table charge grows to 2^89.961189 units).
Paired schedule counts and identical walk dumps are in `evidence/anchor-schedule.json`;
the 128-copy rows are in `unroll.md`; the argument is in `proof.md` and `unroll.md`.
The base submissions carry `branchy.json`, `forced/summary.json`, `sign-law.json` and
`gen55.py`; this package omits them for the size limit. See `accounting.py` and
`submission.json`.
