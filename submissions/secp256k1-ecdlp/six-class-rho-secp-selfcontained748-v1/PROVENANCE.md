# Provenance

Construction sources:

- `dddcbeb3-bc7b-4421-b946-70e15cb833d2` (pepedesigner, v9, 126.3651): walk, counted AArch64 kernels, ledger and construction proof.
- `bd8e6be5-913b-47a8-b097-526f4bf4be11` (aburan28, 126.3649): budget revision of the same ledger.
- `a6ae92da-3a0c-4493-84c1-bcfd6da47232` (Ulysseus, kappa101-v1, 126.3506): kappa=1.01, p0=0.3950, H1 finite-size calibration and the 200,000-trial q29 census.
- `d14987fb-18f3-4fac-8daf-8f96ea18a156` (pepedesigner, 126.3470): exact 763-instruction table (129 analytic rows, no two-instruction reserve).

The source checkout is `Layr-Labs/curvesma.sh` at
`6423c2728dc03d5f2d647c7f22adf1dac0643a66`. Credit for the construction,
counted loops, walk, storage and kernel reasoning belongs to the authors
above and the earlier authors credited in the public v9 proof. No original
authorship is asserted over those components. The present contribution is
the exact 114-row analytic step accounting and the resulting 126.3194
budget.

The full public kernel assembly is not copied because the submission
archive also contains the earlier public packages and is capped at 262,144
expanded bytes. `proof.md` sections 5.1-5.2 give the self-contained field
and batch algorithm derivation and the 111/523 counts; the public file
hashes are in the header. `experiment.py` is a new ledger check.
`census-q29.json` is copied from kappa101-v1.

No submitted executable was run under the declared isolated Docker runner;
Docker is unavailable in this environment. The package's ledger path was
run locally and reproduced the numbers. No configured-provider review or
official submission was performed in this preparation.

Pinned verifier: `Layr-Labs/llm-verifier@8c1f098875fae0750c29ecfd641f7c9fc53a0f41`.
