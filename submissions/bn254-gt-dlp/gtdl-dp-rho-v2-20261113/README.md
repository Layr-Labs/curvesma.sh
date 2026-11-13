# gtdl-dp-rho-v2-20261113

Gated segmented distinguished-point Pollard rho on the BN254 target group, with an
analytic scalar AArch64 price for every instruction charged.

- `submission.json` — claim, exact score, experiment entry point
- `proof.md` — complete argument, success support, heuristics, limits
- `machine.md` — cell-level instruction accounting and the measured cross-check
- `evidence.py` — executable checks, production-GT recovery/statistics, toy
  campaign, ledger, rebuild
- `kernel.c` — portable C instance of the tower; modes `g` (g^r = 1), `s`
  (2,000,000 production walk steps), `t` (fixed product), plus the counting loops
- `kernel-counts.json`, `gt-stats.json` — measured records
- `declared.json` — rounded declarations reproduced by `part=rebuild`

The claim is a generic square-root attack priced end to end. It improves on this
benchmark's promoted baseline (2^135.216266) and is explicitly dominated by the
published STNFS submissions on this track (about 2^88.4), which this entry does not
reproduce or assess.
