# secp256k1-ecdlp: batched-inversion BSGS

Candidate for exploratory review. Deterministic scalar recovery on the full
`production-v2` secp256k1 target with claimed `time_log2 = 128.35`, against the
organizer's 129.000802 accounting control.

- [Argument and cost ledger](proof.md)
- [Manifest](submission.json)
- [Toy-scale correctness evidence](experiment.py)

The improvement is implementation cost, not cryptanalysis: affine addition with
Montgomery batch inversion (one inversion per 2^16 additions instead of one per
addition) and an open-addressing hash table instead of a sorted table. The search
structure is unchanged classical BSGS with `m = 2^128`, `h = 2^128 - 1`, success
probability 1, 2^133.32 bytes of table memory and no advice.

Not a practical break, not a best-known-attack estimate, no full-size execution.
