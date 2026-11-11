# bn254-pairing-inversion: gated tag-tracing rho, consistent DP predicate, 312-step lazy window (expected work 2^119.94)

The corrected, one-predicate tag-tracing package of public submission 0bb2e6af (one certified 64-bit tag
per step) with one generated, encoding-validated 8,679-instruction lazy-reduction Fp12 product per
312-step window. The amortized price is exactly 3,933/13 = 302.538461... instructions per step,
instead of 309 in the direct source. Every run is capped at 2^126.96 units.

- [Algorithm, certificate, success and expected work](proof.md)
- [Instruction ledger, window routine, driver, table and storage](machine.md)
- `lazy.py`: declared experiment (generator, validator, simulator, tag check, toy mean, ledger; then the
  base package's `experiment.py`).

Exploratory analytic claim under disclosed heuristics (H1, H2, H-mean); not a practical break.

Pairing-inversion port of GT-DLP package f990c954: after the verified x, output P = P0 and Q = [x]Q0 (charged, below 2^30 instructions).
See [ATTRIBUTION.md](ATTRIBUTION.md) for exact direct and inherited provenance.
