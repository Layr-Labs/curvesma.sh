# gtdl-dp-rho-v3-20261113

Gated segmented distinguished-point Pollard rho on the BN254 target group, with an
analytic scalar AArch64 price for every instruction charged.

Revision notes: the two earlier candidates on this track were rejected. This one
removes their defects — no borrowed reference price, the 2,000,000-step C walk
starts at `g^a z^b` and reads canonical-form windows, tower additions counted (287),
the segment budget clamped, and fresh-segment exponentiations charged in a 2^40
one-time term. See `proof.md` section 0.

Files: `submission.json` (claim and experiment), `proof.md` (argument), `machine.md`
(accounting), `evidence.py` (checks, GT statistics, toys, ledger), `kernel.c` (C
instance; modes `g`, `s`, `t`, plus the counting loops), `kernel-counts.json`,
`gt-stats.json`, `declared.json`, `target-profile.json`.

The claim is a generic square-root attack priced end to end. It improves on this
benchmark's promoted baseline (2^135.216266) and is explicitly dominated by the
published STNFS submissions on this track (about 2^88.4), which it does not
reproduce or assess.
