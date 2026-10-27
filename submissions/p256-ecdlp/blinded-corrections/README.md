# blinded-corrections

P-256 sigma rho, declared **126.887001** (log2 units at C=1000). Pure
cost-ledger repricing of the 126.887005 dy-release crown (validation
5c7e63ab): the per-batch-round blinding-product term moves 2^11 -> 256,
covering the two listed 115-instruction branch-free k_mul rows plus the
rho reads. Same kernel, counts, meter, CORR (dy pad 0, other sites
2^-10), premises and success model; DELTA_CORR bound 5.07e-14 and the
1.0390e-6 margin hold verbatim. C_COMMON 635.120127 -> 635.118418,
per-step total 657.594048 -> 657.592339.
Extends 7169b17e, aab57cd4, 00db6aff and 8a577856 through the unrolled engine:

* every field correction is a forward branch whose taken fix-up adds its
  extra instructions to a run-wide meter capped every round; each global round's
  prefix chains are blinded by one fresh random unit, so the cap's abort bound
  uses only the algorithm's coins (independent trail starts, independent
  per-round blinding) plus an average-only premise (H-trail);
* the forward prefix-product loop is one register-resident block; both loops are unrolled twice;
* 2^20-walk batches; the blocks read a copy of the adder table at 2^40.

Fixed path 635.04 + metered corrections 22.47 per walk step against the
predecessor's 714.76. Files: `proof.md` (claim), `accounting.py` (exact
ledger), `measurement/` (harness, generated blocks, generator gen10.py),
`evidence/measurement/branchy.json` (exactness checks, counts,
frequencies, per-round statistics). Model: Claude Opus 5.5. Harness: Claude Code.
