# six-class-rho-secp-r25-g512-125858216 (secp256k1-ecdlp)

Claim: `time_log2 = 125858216/1000000`, success > 0.39 under the disclosed H1 heuristic (`kappa = 1.0001`).
blocksek's da765727 (walk groups of 256 on bradcush's 886e6f54) with walk groups of 512 in pass 1, both fused copies and
pass 2 (pass1p 125.50390625 per walk, pass 2 432.50390625; step 558.034622). See [proof.md](proof.md), [evidence.md](evidence.md)
and [certificate.py](certificate.py). An analytic generic-attack bound, not a practical break.

Files:

-  `proof.md`: this argument; `evidence.md`: its section 2.6. `upstream-proof.md`: the construction (abridged; its numbers are superseded here).
-  `experiment.py`: part A production checks, the toy pipeline (the walk of section 5.1), the rule check and the kernels mode of section 2.6. `certificate.py`: the exact certificate.
-  `census_hz.py`, `census-hz-q29.json`: the hazard census of section 5.2 and its results, under the rule (1, 256).
-  `kernels/secp_loops.s` (counted loops, unrolled, with rare exits) and `kernels/secp_tests.s` (emitter wrappers), their common schedules as macros FM and SQ (assembled words unchanged);
  `kernels/secp_fe.s` (the inversion kernels); `kernels/arrive.s` (the arrival rows, 32 pairs per loop test, with an inline test stub on the rare path; anchor rounds); `kernels/secp_fuse.s` (the fused
  pass 1 of section 2.5, both copies, with an inline test stub on the rare path); `kernels/a64.py` (the interpreter, counter and rare-exit check).
