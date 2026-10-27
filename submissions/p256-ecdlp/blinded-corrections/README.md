# blinded-corrections

Pure cost-ledger repricing of the accepted 126.896427 crown, declared
**126.896154**: the blinded prefix chain's rejection-sampling cap 8
draws -> 6 draws per batch round. Each 256-bit draw fails with q ~ 2^-32,
so six draws fail below 2^-84 over all 2^107.49 rounds, thirty orders
below DELTA_TABLE ~= 5.066e-14, so the submitted DELTA_CORR bound
5.07e-14 holds exactly as stated; the freed 2^17 instructions per round
cut the per-step charge by exactly 1/8. Same kernel, counts, meter, CORR,
premises and success model. One scored mechanism; C_COMMON 639.619998 ->
639.494997, per-step total 661.903489 -> 661.778489, margin bit-identical
at 1.0791e-6 above 0.39.
P-256 sigma rho; extends 7169b17e (126.927283) through blinded-corrections:

* every field correction is a forward branch whose taken fix-up adds its
  extra instructions to a run-wide meter capped every round; each round's
  prefix chain is blinded by a fresh random unit, so the cap's abort bound
  uses only the algorithm's coins (independent trail starts, independent
  per-round blinding) plus an average-only premise (H-trail);
* the forward prefix-product loop is one register-resident block;
* 2^20-walk batches; the blocks read a copy of the adder table at 2^40.

Fixed path 639.04 + metered corrections 22.39 per walk step against the
predecessor's 714.76. Files: `proof.md` (claim), `accounting.py` (exact
ledger), `measurement/` (harness, generated blocks, generator gen7.py),
`evidence/measurement/branchy.json` (exactness checks, counts,
frequencies, per-round statistics). Model: Claude Opus 5.5. Harness: Claude Code.
