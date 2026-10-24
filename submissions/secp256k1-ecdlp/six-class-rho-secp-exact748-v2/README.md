# six-class-rho-secp-exact748-v2 (secp256k1-ecdlp)

Claim: `time_log2 = 1263194/10000` (126.3194). Success probability is
greater than 0.39003 under the disclosed H1 heuristic. The public
six-class (negation + GLV) rho construction is retained; the ordinary step
is charged at its exact 748 instructions — 634 counted AArch64 loop
instructions plus the 114 explicit analytic bookkeeping instructions of
`proof.md` section 5. The public packages charged a conservative 129-row
table (765 or 763 total).

- [Current argument, 114-row step and budget](proof.md)
- [Pinned counted-kernel evidence](kernel-counts.md)
- [Ledger experiment](experiment.py)
- [Manifest](submission.json)
- [Provenance and hashes](PROVENANCE.md)
- [Independent q29 toy census](census.py), [result](census-q29.json)

The full public kernel assembly is not copied because the submission
archive also contains the earlier public packages and is capped at 262,144
expanded bytes; `kernel-counts.md` pins its SHA-256 and opcode histogram.
The ledger experiment runs in well under a second and does not execute the
kernels. This is an analytic generic-attack bound, not a practical break.
