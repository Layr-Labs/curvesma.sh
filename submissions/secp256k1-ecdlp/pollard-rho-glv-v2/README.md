# secp256k1-ecdlp — Pollard rho over GLV automorphism equivalence classes (v2)

Candidate for exploratory review. Pollard rho with distinguished points over
the order-6 automorphism classes of secp256k1 (negation + GLV endomorphism):
claimed total work 2^127.22 target-group-operation equivalents at success
probability 0.405, under the disclosed random-mapping heuristic H1.

Revision of the v1 candidate in the sibling folder: same algorithm and walk
analysis; the accounting changes are a single-multiplication
canonicalisation (`beta^2*x = -(x + beta*x) mod p`, exact) and a success
exponent recomputed for probability 0.405 (threshold 0.39 + 0.015 margin
after all deductions). See proof.md Sections 7 and 9.

- [Algorithm, probability and cost argument](proof.md)
- [Submission claim](claim.json)

The authenticated workflow dossier determines qualification. This is an
analytic attack estimate under the fixed cost model, not a practical break
or an executed production attack.
