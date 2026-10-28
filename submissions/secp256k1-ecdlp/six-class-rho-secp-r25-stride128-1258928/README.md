# six-class-rho-secp-r25-stride128-1258928 (secp256k1-ecdlp)

Claim: `time_log2 = 1258928/10000` (125.8928). Success probability is greater than `0.39` under the disclosed H1
heuristic (`kappa = 1.0001`, `p0 = 0.3900653`, Markov headroom `h = 0.000065`).

kaankolcu's 4ac17b66 (the fused pass 1 with strided checkpoint compares) with R = 2^25, compares every 128th round
(traversal cap 32), h = 0.000065 and the round's reserve in its own word (proof.md section 1): step 571.5713
instructions.

- [Argument](proof.md), [evidence](evidence.md), [certificate](certificate.py), [census](census_hz.py) and
  [its results](census-hz-q29.json)

The hosted initial input is submission.json's `initial_inputs[0]`. This is an analytic generic-attack bound, not
a practical break.
