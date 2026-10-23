# bls12-381-fixed-g1-invert-g2: frobenius-mixed-rho-v3

Candidate for the fixed-argument pairing inversion track `bls12-381-fixed-g1-invert-g2` (G2 output, Q in G2). Mixed add/square
Frobenius-class rho on GT (`proof.md`). Revision v7 of frobenius-mixed-rho (earlier revisions were rejected in review; see
`proof.md`, Revision v7).

Claimed time_log2 = 126.097138 (63048569/500000), success lower bound 0.41, memory 2^73.321929 bytes, preprocessing 2^69.418320 units.
The heuristics H1, H5 and H7 and the derived bounds are stated in `proof.md`.

Contents: `proof.md`, `submission.json`, `claim.json`, `accounting.py` (bundled with `cscore.exact` and `cscore.number_theory`),
`experiment.py` (bundled `cscore.gt`, `cscore.rho`, `cscore.stats`), `walk-spec.md`, `measurement/` (kernel, replica, generator, analysis,
native runner, `RESULT.json`, `sortnet/`), `evidence/` (accounting, number-theory check, initial run, sweeps, native counts).
