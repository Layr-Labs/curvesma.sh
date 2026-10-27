# blinded-corrections

P-256 sigma rho, declared **126.887005** (log2 units at C=1000). Pure
cost-ledger repricing of the 126.887016 twice-unrolled crown (validation
ecafd1c3): the dy correction site's per-site padding is released
(2^-10 -> 0) because dy borrows with probability exactly 1/2 by L-ori
sign symmetry, so there is no adder-mean deviation for padding to cover;
all other sites keep 2^-10. CORR 22.473310 -> 22.468426 (exact 5/1024
cut), per-step total 657.598931 -> 657.594048. Same kernel, counts,
meter, premises and success model; DELTA_CORR bound 5.07e-14 and the
1.0390e-6 margin hold verbatim. Extends 7169b17e, aab57cd4, 00db6aff and 8a577856 through the unrolled engine:

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
