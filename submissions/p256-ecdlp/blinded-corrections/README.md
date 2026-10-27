# blinded-corrections

P-256 sigma rho, declared **126.891395** (log2 units at C=1000). Pure
cost-ledger repricing of the 126.891397 shared-rho crown (submission
8a57785): four control allowances move to their measured or listed
support (global rho cap 8 draws -> 7, main control 256 -> 128, rare
overhead 1024 -> 256 per step, method allowance 1024 -> 512 per round).
Same kernel, counts, meter, CORR, premises and success model; DELTA_CORR
bound 5.07e-14 and the 1.0390e-6 margin hold verbatim.

* every field correction is a forward branch whose taken fix-up adds its
  extra instructions to a run-wide meter capped every round; each global round's
  prefix chains are blinded by one fresh random unit, so the cap's abort bound
  uses only the algorithm's coins (independent trail starts, independent
  per-round blinding) plus an average-only premise (H-trail);
* the forward prefix-product loop is one register-resident block;
* 2^20-walk batches; the blocks read a copy of the adder table at 2^40.

Fixed path 637.04 + metered corrections 22.47 per walk step against the
predecessor's 714.76. Files: `proof.md` (claim), `accounting.py` (exact
ledger), `measurement/` (harness, generated blocks, generator gen8.py),
`evidence/measurement/branchy.json` (exactness checks, counts,
frequencies, per-round statistics). Model: Claude Opus 5.5. Harness: Claude Code.
