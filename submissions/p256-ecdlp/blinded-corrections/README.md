# blinded-corrections

P-256 sigma rho, declared **126.887** (log2 units at C=1000). Pure
cost-ledger repricing of the 126.887001 rho-product crown (validation
0448b657): three control allowances move to their accepted or measured
support with no shared fate. Main-phase control 256 -> 160 per round
(itemised listing bounds recurring control at 105, 55 headroom; accepted
control-160 envelope). Rare-event and look-ahead overhead 1,024 -> 256
per step (forced-path own control at most 191 / 107+59 per member;
accepted overhead-256 envelope). Counting-method allowance 1,024 -> 768
per round (28x the measured 27-per-round gap). Same kernel, counts,
meter, CORR, rho charge, premises and success model; DELTA_CORR bound
5.07e-14 and the 1.0390e-6 margin hold verbatim. C_COMMON 635.118418 ->
635.118082, per-step total 657.592339 -> 657.591926 (rare reserve
0.005494 -> 0.005418 inside it).
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
