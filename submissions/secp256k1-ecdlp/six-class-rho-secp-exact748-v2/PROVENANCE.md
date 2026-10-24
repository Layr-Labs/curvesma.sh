# Provenance and verification status

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

This package is deliberately small because the submission archive also
contains the earlier public packages and is capped at 262,144 expanded
bytes. The public kernel assembly is not copied; `kernel-counts.md` pins
its SHA-256, per-iteration counts and opcode histogram. `experiment.py` is
a new standard-library ledger check. `census.py` and `census-q29.json` are
copied from kappa101-v1 and are supplementary toy evidence, not the
declared experiment.

No submitted executable was run under the declared isolated Docker runner;
Docker is unavailable in this environment. The package's own ledger path
was run locally and reproduced the numbers in `proof.md`. No
configured-provider review or official submission has been performed in
this preparation.

Pinned verifier: `Layr-Labs/llm-verifier@8c1f098875fae0750c29ecfd641f7c9fc53a0f41`.

| Path | SHA-256 |
| --- | --- |
| `proof.md` | `986ab3da966dd69967107f0d6c0ad2cafdd40ecfeef8f77d6a01e9f5b7acbb4b` |
| `experiment.py` | `3ad0da34a21c88007fc04593635dc123beffc834fea58c5b3ad4217505d8a3fa` |
| `kernel-counts.md` | `2c20a0aa84f674dbb8eb7947b877e4786863f149ed5dc82590892c1dac4e4269` |
| `census.py` | `4e87a10efb5dd9dea4e61ea61bdb17a6ea5524fae3c7af0a719ca4fb133df16f` |
| `census-q29.json` | `546a7c0d1bde79df1f3b48669e0d2c3cc27405fb892b71042f8caa70ee055ce1` |
