# bls12-381-fixed-g1-invert-g2: frobenius-mixed-rho

Candidate for the fixed-argument pairing inversion track `bls12-381-fixed-g1-invert-g2`. Mixed add/square Frobenius-class
rho on GT (`proof.md`). Revision v2 (of v1 frobenius-mixed-rho, claim 125.700213). Claimed time_log2 = 125.747292 (31436823/250000), success lower bound
0.41, memory 2^73.200717 bytes, preprocessing 2^75.452436 units. Heuristics H1-H4 are in `proof.md`.

Status: not submitted. Nothing in this package has been posted or pushed. Attribution for any
later submission: model "Claude Haiku 5.5", harness "Claude Code".

Contents: `proof.md`, `submission.json`, `claim.json`, `accounting.py` (bundled with `cscore.exact`),
`experiment.py` (bundled `cscore.gt`, `cscore.rho`, `cscore.stats`), `walk-spec.md`, `measurement/`
(kernel, replica, generator, analysis, `RESULT.json`), `evidence/` (accounting, initial run, sweeps).
