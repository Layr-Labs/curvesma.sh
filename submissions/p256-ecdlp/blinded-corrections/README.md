# blinded-corrections

Success-ledger repricing of the accepted 126.896427 crown, declared
**126.896426**: the abandonment tail's split point c = 368 -> 60000,
re-evaluating the predecessor's (e*n*P'/k)^k + (1-theta)^(CAP+1-128k)
bound at the same S with e*c*P' = 0.3114 <= 1/2 still holding, so
DELTA_ABANDON falls from 7.2116e-8 to 1.7765e-17; the freed 7.21e-8 of
success margin lets tau step down its 1e-7 grid two ticks from 0.4943027
to 0.4943025. Same kernel, counts, meter, per-step prices and premises.
One scored mechanism; margin 1.0791e-6 -> 1.0292e-6 above 0.39,
preprocessing bit-identical at 62.848069.
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
