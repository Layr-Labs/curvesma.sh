# negation-expanded-exact-interval-bsgs

Deterministic baby-step/giant-step for the plain BLS12-381 G1 discrete logarithm.
Nonzero baby points are stored under both compressed encodings; the negative key
is a sign-flag copy. Giant steps use stride `2m-1`.

- [Claim and score](submission.json)
- [Algorithm, exactness proof and complete work ledger](proof.md)
- [Small-scale reference implementation](reference.py)
- [Toy coverage experiment](experiment.py)

The recoverable set is exactly `[0, h*s-m]` together with `[r-(m-1), r)`, so the
success probability is exactly `h*s/r >= 39/100`. At most one candidate is
emitted and verified. Submitted `time_log2` is 127.995855, peak memory
`2^133.982508` bytes, preprocessing `2^126.959413` target-group-operation
equivalents (included in total work), advice 0 bytes.

This is an analytic generic accounting control. It is not a practical break, a
best-known attack estimate or a cryptanalytic advance, and no full-size run is
claimed.
