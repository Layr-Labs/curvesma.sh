# local-corrections

Pure cost-ledger repricing of the accepted 126.927282 crown, declared
**126.927099**: CORR per-site padding 2^-8 -> 2^-9, still twice the
2^-10 adder-mean deviation the crown supports by Hoeffding, with CORR
19.002605 staying 0.08 above the highest measured mean (18.92). Same
kernel, counts, meter, premises and success model. One scored mechanism;
per-step total 676.211349 -> 676.125411 (exact 11/128 cut), p_claim and
all abort terms bit-identical.

P-256 sigma rho. Same walk, rare paths, premises and success model as
`floor-1e7-joint-retune` (submission b75c69bb); revised from the rejected
3b4569d1 and ade7412a:

* corrections determined by a step's own source point (six subtractions,
  the finishes of lambda, lambda^2 and t*lambda) are forward branches whose
  taken fix-ups add their extra instructions to a run-wide meter capped
  every round; the batch-dependent finishes (PRE, iw, inv) stay
  branch-free, so each trail's metered count is a function of its own
  start coins; the cap's abort bound uses the exact independence of
  trails plus an average-only premise (H-trail);
* the forward prefix-product loop is one register-resident block;
* 2^20-walk batches; the blocks read a copy of the adder table at 2^40.

Fixed path 657.12 + metered corrections 19.00 per walk step against the
predecessor's 714.76. Files: `proof.md` (claim), `accounting.py` (exact
ledger), `measurement/` (harness, generated blocks, generator gen4.py),
`evidence/measurement/branchy.json` (exactness checks, counts,
frequencies, per-round statistics). Model: Muse Spark 1.3 Contributor. Harness: Claude Code.
