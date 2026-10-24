# orbit-canonical-glv-bsgs

Deterministic baby-step/giant-step for the plain BLS12-381 G1 discrete logarithm.
Baby points are `a*P0` and giant points are `T-b*phi(P0)`, where `phi` is the
order-3 GLV endomorphism and `phi(P0)=[lambda]P0` with `lambda^2+lambda+1=r`.
Each point is stored once, under the minimum compressed encoding of its
order-6 orbit (the endomorphism together with negation).

- [Claim and score](submission.json)
- [Algorithm, exact coverage argument and work ledger](proof.md)
- [Production parameter check](reference.py)
- [Toy coverage and generator identity](experiment.py)

The recoverable set contains an explicit set `S` of cardinality
`(m-1)*(3*m-1)`, so the success probability is at least `|S|/r >= 39/100`.
Submitted `time_log2` is 127.863728, peak memory `2^133.126762` bytes,
preprocessing `2^126.845142` target-group-operation equivalents (included in
total work), advice 0 bytes.

This is an analytic generic accounting control. It is not a practical break,
a best-known attack estimate or a cryptanalytic advance, and no full-size run
is claimed.
