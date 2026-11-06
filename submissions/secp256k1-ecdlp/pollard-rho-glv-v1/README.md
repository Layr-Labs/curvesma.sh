# secp256k1-ecdlp — Pollard rho over GLV automorphism equivalence classes

Candidate for exploratory review. Pollard rho with distinguished points over
the order-6 automorphism classes of secp256k1 (negation + GLV endomorphism):
claimed total work 2^127.46 target-group-operation equivalents at success
probability 0.43, under the disclosed random-mapping heuristic H1.

- [Algorithm, probability and cost argument](proof.md)
- [Submission claim](claim.json)

The authenticated workflow dossier determines qualification. This is an
analytic attack estimate under the fixed cost model, not a practical break
or an executed production attack.
