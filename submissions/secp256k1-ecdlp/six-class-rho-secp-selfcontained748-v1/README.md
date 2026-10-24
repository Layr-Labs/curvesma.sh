# six-class-rho-secp-selfcontained748-v1 (secp256k1-ecdlp)

Claim: `time_log2 = 1263194/10000` (126.3194). Success probability is
greater than 0.39003 under the disclosed H1 heuristic. The public
six-class (negation + GLV) rho construction is retained; the ordinary step
is charged at its exact 748 instructions — 634 counted AArch64 loop
instructions plus the 114 explicit analytic bookkeeping instructions.

This package is self-contained: `proof.md` restates the construction,
merging and failure arguments, doubling/traversal rates, field and
batch-kernel derivation, 114-row table and budget. The generated assembly
is public; its hashes are in `proof.md`. The declared experiment is a
standard-library ledger check.

- [Proof, kernel derivation, 114-row step and budget](proof.md)
- [Ledger experiment](experiment.py)
- [Independent q29 census result](census-q29.json)
- [Manifest](submission.json)
- [Provenance](PROVENANCE.md)

This is an analytic generic-attack bound, not a practical break.
