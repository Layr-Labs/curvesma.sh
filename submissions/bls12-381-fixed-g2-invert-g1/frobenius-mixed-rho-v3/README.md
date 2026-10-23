# bls12-381-fixed-g2-invert-g1: frobenius-mixed-rho-v3

Candidate for the fixed-argument pairing inversion track `bls12-381-fixed-g2-invert-g1` (G1 output, P in G1). Mixed add/square
Frobenius-class rho on GT (`proof.md`). Revision v6 of frobenius-mixed-rho (earlier revisions were rejected in review; see
`proof.md`, Revision v5).

Claimed time_log2 = 125.789684 (31447421/250000), success lower bound 0.41, memory 2^73.166564 bytes, preprocessing 2^76.418283 units.
The heuristics H1, H5 and H7 and the derived bounds are stated in `proof.md`.

Contents: `proof.md`, `submission.json`, `claim.json`, `accounting.py` (bundled with `cscore.exact`), `experiment.py`
(bundled `cscore.gt`, `cscore.rho`, `cscore.stats`), `walk-spec.md`, `measurement/` (kernel, replica, generator, analysis,
native runner, `RESULT.json`), `evidence/` (accounting, initial run, sweeps, native counts).
