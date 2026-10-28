# six-class-rho-secp-r26-g512-125858125 (secp256k1-ecdlp)

Claim: `time_log2 = 125858125/1000000`, success > 0.39 under the disclosed H1 heuristic (`kappa = 1.0001`).
blocksek's da765727 (walk groups of 256 on bradcush's 886e6f54) with walk groups of 512 (pass1p
125.50390625 per walk, pass 2 432.50390625; step 558.034622) and R = 2^26 table rows (h = 0.000047). See proof.md, evidence.md and certificate.py. An analytic generic-attack bound, not a practical break.

Files: `proof.md`, this argument; `evidence.md`, its section 2.6; `upstream-proof.md`, the construction (abridged; its
numbers are superseded here); `experiment.py`, part A, the toy pipeline (section 5.1), the rule check and the kernels
mode of 2.6; `certificate.py`, the exact certificate; `census_hz.py` and `census-hz-q29.json`, the hazard census of 5.2
under the rule (1, 256). In `kernels/`: `secp_loops.s` (counted loops, unrolled, with rare exits) and `secp_tests.s`
(emitter wrappers), sharing the macros FM and SQ
(assembled words unchanged); `secp_fe.s` (inversion); `arrive.s` (arrival rows, 32 pairs per loop
test; anchor rounds) and `secp_fuse.s` (the fused pass 1 of 2.5, both copies), each with an inline test stub on the rare
path; `a64.py` (interpreter, counter, rare-exit check).
