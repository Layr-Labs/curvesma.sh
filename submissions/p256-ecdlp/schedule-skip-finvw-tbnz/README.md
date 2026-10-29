# schedule-skip-finvw-tbnz (three-arm anchor schedule, A=2^27, windowed inversion)

P-256 sigma rho, declared **126.732475** (5069299/40000 at C=1000), on the top-limb
distinguished-test package (126.735172). Five changes: `finv` is 4-bit windowed (253
squarings, 39 multiplications), saving 11,814 instructions per round; the 2^60 anchor cap
is one `tbnz x12,#60`; every unbatched rare-path addition shares `finv`, so F_MAX
41,190 -> 29,376 and wL 2,789,700 -> 2,021,790; the adder table grows to A=2^27, which
cuts the two-cycle reserve to 0.0075319 per step; and the walk compares against its anchor
only in rounds `t mod 64 = 0` and the plain rounds `33..63`, using a compare-free tail for
offsets 1..32 (2*(B-1)-2 fewer instructions per skip round).

PLAIN_FIXED 605,581,542; anchor extra 3*(2^20-1)+9; skip arm 605,581,542-2*(2^20-1)+2;
C_COMMON 576.577738; CORR 14.210693; reserve 0.0075889; total step 590.796019.
Preprocessing 2^84.961190 units, peak memory below 2^77 bytes.

Paired organizer-image counts: plain 39,222,158 -> 39,210,344 (11,814); anchor
39,485,165 -> 39,407,816 (77,349); anchor-skip deltas at B=2^16 (rounds 1, 33, 64) and
B=2^20 (round 1) reproduce the arm model exactly, with identical correction meters and
byte-identical one-round walk dumps. The referenced branchy.json checkpoints are in the
base submissions; this package omits that file for the 256 KiB limit. See proof.md,
submission.json, evidence/anchor-schedule.json and accounting.py.
