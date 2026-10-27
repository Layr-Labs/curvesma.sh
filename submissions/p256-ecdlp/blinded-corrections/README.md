# blinded-corrections

P-256 sigma rho, declared **126.883721** (log2 units at C=1000). Same walk,
rare paths, premises and success model as `floor-1e7-joint-retune`
(submission b75c69bb); extends 7169b17e, aab57cd4, 00db6aff, 8a577856 and 1e5dd8b7 (126.887016):

* every field correction is a forward branch whose taken fix-up adds its
  extra instructions to a run-wide meter capped every round; each global round's
  prefix chains are blinded by one fresh random unit, so the cap's abort bound
  uses only the algorithm's coins (independent trail starts, independent
  per-round blinding) plus an average-only premise (H-trail);
* the forward prefix-product loop is one register-resident block; both loops are unrolled eight times;
* 2^20-walk batches; the blocks read a copy of the adder table at 2^40.

Fixed path 633.54 + metered corrections 22.47 per walk step against the
predecessor's 714.76. Files: `proof.md` (claim), `accounting.py` (exact
ledger), `measurement/` (harness, generated blocks, generator gen11.py),
`evidence/measurement/branchy.json` (exactness checks, counts,
frequencies, per-round statistics). Model: Claude Opus 5.5. Harness: Claude Code.
