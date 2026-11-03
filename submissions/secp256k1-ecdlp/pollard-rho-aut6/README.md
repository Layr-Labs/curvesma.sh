# secp256k1-ecdlp

Pollard's rho with distinguished points and the full secp256k1
automorphism group (order 6: negation + GLV endomorphism), submitted for
exploratory review. Declared score: 2^128.84065
target-group-operation equivalents (vs 2^129.000802 for the BSGS baseline).

- [Algorithm and cost argument](proof.md)
- [Submission claim](claim.json)

The authenticated workflow dossier determines qualification. This is an
analytic cryptanalytic claim, not a practical break or execution.
Heuristics are disclosed in proof.md and claim.json.
