# bls12-381-fixed-g1-invert-g2: frobenius-mixed-rho

Candidate for the fixed-argument pairing inversion track `bls12-381-fixed-g1-invert-g2`. Mixed add/square Frobenius-class
rho on GT (`proof.md`). Claimed time_log2 = 125.700213 (125700213/1000000), success lower bound
0.39, memory 2^73.153638 bytes, preprocessing 2^75.405357 units. Heuristics H1-H4 are in `proof.md`.

Status: not submitted. Nothing in this package has been posted or pushed. Attribution for any
later submission: model "Claude Haiku 5.5", harness "Claude Code".

Contents: `proof.md`, `submission.json`, `claim.json`, `accounting.py` (bundled with `cscore.exact`),
`experiment.py` (bundled `cscore.gt`, `cscore.rho`, `cscore.stats`), `walk-spec.md`, `measurement/`
(kernel, replica, generator, analysis, `RESULT.json`), `evidence/` (accounting, initial run, sweeps).
