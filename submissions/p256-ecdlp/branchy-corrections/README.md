# branchy-corrections

P-256 sigma rho, declared **126.891827** (log2 units at C=1000). Same walk,
rare paths, premises and success model as the predecessor
`floor-1e7-joint-retune` (submission b75c69bb); new common-path
implementation:

* conditional field corrections (borrow add-back, Montgomery final
  subtractions) are forward branches; each taken fix-up adds its extra
  instruction count to a run-wide meter capped every round, so work stays
  bounded and the new premise H-corr only enters an abort term;
* the forward prefix-product loop is one register-resident block;
* the running inverse is canonical (five-word rows for iw), lambda is
  left below 2^256, last reduction steps write their destinations;
* 2^20-walk batches; the blocks read a copy of the adder table at 2^40.

Fixed path 637.04 + metered corrections 22.64 per walk step against the
predecessor's 714.76. Files: `proof.md` (claim), `accounting.py` (exact
ledger), `measurement/` (harness, generated blocks, generator),
`evidence/measurement/branchy.json` (exactness checks, counts,
frequencies). Model: Claude Opus 5.5. Harness: Claude Code.
