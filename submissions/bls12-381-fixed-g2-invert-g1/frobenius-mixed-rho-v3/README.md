# bls12-381-fixed-g2-invert-g1: frobenius-mixed-rho-v3

Candidate for the fixed-argument pairing inversion track `bls12-381-fixed-g2-invert-g1` (G1 output, P in G1). Mixed add/square
Frobenius-class rho on GT (`proof.md`). Revision v3 of frobenius-mixed-rho (v2 was rejected by the reviewer on cycle
detection, the collision rate, the c2 extraction and the memory record; see `proof.md`, Revision v3).

Claimed time_log2 = 125.719085 (25143817/200000), success lower bound 0.41, memory 2^73.166751 bytes, preprocessing 2^75.418470 units.
Heuristics H1-H4 are in `proof.md`.

Contents: `proof.md`, `submission.json`, `claim.json`, `accounting.py` (bundled with `cscore.exact`), `experiment.py`
(bundled `cscore.gt`, `cscore.rho`, `cscore.stats`), `walk-spec.md`, `measurement/` (kernel, replica, generator, analysis,
native runner, `RESULT.json`), `evidence/` (accounting, initial run, sweeps, native counts).
