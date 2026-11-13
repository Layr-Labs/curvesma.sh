# gtdl-dp-rho-20261113

Gated segmented distinguished-point Pollard rho on the BN254 target group.
`submission.json` declares the claim and the experiment entry point;
`proof.md` is the complete argument, `machine.md` the instruction accounting,
`kernel.c` + `kernel-counts.json` the executable arithmetic and its measurement.
`declared.json` holds the rounded declarations that `evidence.py part=rebuild`
reproduces.

The claim is a generic square-root attack priced end to end. It improves on this
benchmark's promoted baseline (2^135.216266) and is explicitly dominated by the
published STNFS submissions on this track (about 2^88.4), which this entry does
not reproduce or assess.
