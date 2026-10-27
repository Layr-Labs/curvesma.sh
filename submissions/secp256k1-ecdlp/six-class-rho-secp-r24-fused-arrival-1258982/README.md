# six-class-rho-secp-r24-fused-arrival-1258982 (secp256k1-ecdlp)

Claim: `time_log2 = 1258982/10000` (125.8982). Success probability is greater
than `0.39` under the disclosed H1 heuristic (`kappa = 1.0001`, `p0 = 0.3900503`,
Markov headroom `h = 0.00005`).

It is 4c8222cc's package with each walk's arrival fused into the next round's
pass 1 (`kernels/secp_fuse.s`, proof.md section 2.5): step 573.5713 instructions.

- [Argument](proof.md), [certificate](certificate.py), [census](census_hz.py) and [its results](census-hz-q29.json)

The hosted initial input is submission.json's `initial_inputs[0]`.

This is an analytic generic-attack bound, not a practical break.
