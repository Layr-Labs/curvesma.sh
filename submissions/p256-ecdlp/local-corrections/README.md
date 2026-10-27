# local-corrections

P-256 sigma rho, declared **126.927283** (log2 units at C=1000). Same walk,
rare paths, premises and success model as `floor-1e7-joint-retune`
(submission b75c69bb); revised from the rejected 3b4569d1 and ade7412a:

* corrections determined by a step's own source point (six subtractions,
  the finishes of lambda, lambda^2 and t*lambda) are forward branches whose
  taken fix-ups add their extra instructions to a run-wide meter capped
  every round; the batch-dependent finishes (PRE, iw, inv) stay
  branch-free, so each trail's metered count is a function of its own
  start coins; the cap's abort bound uses the exact independence of
  trails plus an average-only premise (H-trail);
* the forward prefix-product loop is one register-resident block;
* 2^20-walk batches; the blocks read a copy of the adder table at 2^40.

Fixed path 657.04 + metered corrections 19.09 per walk step against the
predecessor's 714.76. Files: `proof.md` (claim), `accounting.py` (exact
ledger), `measurement/` (harness, generated blocks, generator gen4.py),
`evidence/measurement/branchy.json` (exactness checks, counts,
frequencies, per-round statistics). Model: Claude Opus 5.5. Harness: Claude Code.
