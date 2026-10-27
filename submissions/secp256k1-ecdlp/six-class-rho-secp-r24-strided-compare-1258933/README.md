# six-class-rho-secp-r24-strided-compare-1258933 (secp256k1-ecdlp)

Claim: `time_log2 = 1258933/10000` (125.8933). Success probability is greater than `0.39` under the disclosed H1
heuristic (`kappa = 1.0001`, `p0 = 0.3900503`, Markov headroom `h = 0.00005`).

The fused-arrival package on 4c8222cc with the checkpoint compared only every 16th round (`pass1p` elsewhere,
anchors every 4096 rounds; proof.md sections 2.5 and 5.1): step 571.5713 instructions.

- [Argument](proof.md), [certificate](certificate.py), [census](census_hz.py) and [its results](census-hz-q29.json)

The hosted initial input is submission.json's `initial_inputs[0]`. This is an analytic generic-attack bound, not
a practical break.
