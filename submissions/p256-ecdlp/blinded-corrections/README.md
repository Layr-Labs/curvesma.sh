# blinded-corrections

Pure cost-ledger repricing of the accepted 126.896669 crown, declared
**126.896427**: CORR per-site padding 2^-8 -> 2^-9, still twice the
2^-10 adder-mean deviation the crown supports by Hoeffding, with CORR
22.277997 staying 0.11 above the highest measured mean (22.17). Same
kernel, counts, meter, premises and success model (tau 0.4943027, delta
2^-24 kept). One scored mechanism; per-step total 662.014817 ->
661.903489 (exact 57/512 cut), p_claim and all abort terms
bit-identical. Accepted promotion 136b8848, restored onto the
birthday-tau engine.
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
