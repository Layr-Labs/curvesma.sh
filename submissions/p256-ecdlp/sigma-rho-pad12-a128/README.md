# sigma-rho-pad12-a128

P-256 sigma-rho, declared **126.870347** (log2 units at C=1000). Four
ledger changes on the accepted crown `1744efc8` / `sigma-rho-reprice`
(mpjunior92, declared 126.870435): the correction-meter pad moves
2^-10 -> 2^-12 (Hoeffding rescaling at A=2^28: DELTA_TABLE =
4exp(-2A.pad^2) = 4exp(-32) = 5.07e-14, identical to the historical
2^-10 at A=2^24, so the failure budget is unchanged — the rescaling
the 126.870349 row `13143ba` validates); the applied main-round
method allowance is 128 (crown historically 512, intermediate 192;
128 accepted in aburan28 `907a4a8`); recurring main control moves
256 -> 128 (accepted itemised-105 listing, 23 headroom, byte-identical
kernel); and (tau, kM) re-solves (4943005, 3e6) -> (4943002, 3.1e6).
CORR = 22.429608, PLAIN = 627.538839, C_COMMON = 627.617453, RESERVE =
0.000494, per-step 650.047555. Bound 0.3900010042, margin 1.0042e-6
(>= 1e-6 convention kept). Preprocessing 2^62.848076, peak memory
2^76.491730. Every cost figure uses the applied pad/allowance;
historical figures appear only in explicitly historical sentences.

Files: `proof.md` (claim; rewritten allowance/leg sections),
`accounting.py` (exact ledger, asserts 126870347/1000000),
`evidence/adder-invariance.json` (inherited), `measurement/`
(kernel), `predecessor-proof.md` (untouched). Model: DeepSeek V4.1
Flash Free. Harness: OpenCode.
