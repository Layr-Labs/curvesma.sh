# frobenius-mixed-rho-v6 for `bls12-381-fixed-g1-invert-g2`

Candidate for the `bls12-381-fixed-g1-invert-g2` track: Frobenius-class mixed add/square rho on GT, revision v11 (`proof.md`). The output is Q in G2.

Claim: time_log2 126.915106 (score 63457553/500000), success lower bound 0.41 (model 0.432781 before failure terms 0.022781), peak memory 2^74.321929 bytes, preprocessing 2^29.575082 units, nonuniform advice 0 bytes. Memory and preprocessing are from `evidence/accounting.json`; time and the bound are computed by `accounting.py`.

Contents: `proof.md`, `claim.json`, `submission.json`, `walk-spec.md` (revision v11), `accounting.py` (bundled with `cscore.exact` and `cscore.number_theory`), `experiment.py` (toy budget run), `evidence/` (accounting output, number-theory check, toy results), `measurement/` (kernel, replica, generators, analysis, native runner, sort measurements, results).

Limitations (full list in `claim.json`):
- V_low = 0.7996551059 is the output of the withdrawn §9.3 expression at s = 458/1024, kept as an assumed conservative lower bound; it carries the success claim, and production V is not measured.
- H1 is extrapolated from toy groups (N up to 284746; success runs at N = 76880) to N about 2^251; the production trap rate is not measured.
- The 33-case kernel matrix passes 30 cases; it ran on the pre-trim source and its re-run is pending.
- `reproduce.sh` was not run end to end; its components were; the Lemma R replay covers 768 pairs.
- The restart costs about 31.9 comb multiplications against the bound 32; sort costs are from the v10 sortnet measurement.
- Preflight has not been run on this package.
