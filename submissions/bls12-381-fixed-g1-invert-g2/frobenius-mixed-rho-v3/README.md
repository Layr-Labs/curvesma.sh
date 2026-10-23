# bls12-381-fixed-g1-invert-g2: frobenius-mixed-rho-v3

Candidate for the fixed-argument pairing inversion track `bls12-381-fixed-g1-invert-g2` (G2 output, Q in G2). Mixed add/square
Frobenius-class rho on GT (`proof.md`). Revision v4 of frobenius-mixed-rho (v3 was rejected in round 3 on long-cycle trapping,
the fingerprint bit count and the missing work cutoff; see `proof.md`, Revision v4).

Claimed time_log2 = 125.719870 (12571987/100000), success lower bound 0.41, memory 2^73.166749 bytes, preprocessing 2^76.418468 units.
Heuristics are listed in `proof.md`.

Contents: `proof.md`, `submission.json`, `claim.json`, `accounting.py` (bundled with `cscore.exact`), `experiment.py`
(bundled `cscore.gt`, `cscore.rho`, `cscore.stats`), `walk-spec.md`, `measurement/` (kernel, replica, generator, analysis,
native runner, `RESULT.json`), `evidence/` (accounting, initial run, sweeps, native counts).
