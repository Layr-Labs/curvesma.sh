# exact-compressed-key-truncated-bsgs

Deterministic truncated baby-step/giant-step for the plain BLS12-381 G1 discrete
logarithm, keyed on the 48-byte compressed affine encoding.

- [Claim and score](submission.json)
- [Algorithm, exactness proof and complete work ledger](proof.md)
- [Small-scale reference implementation](reference.py)
- [Toy coverage experiment](experiment.py)

The recoverable set is exactly the integer interval `[0, m^2)`, so the success
probability is exactly `m^2/r >= 39/100`; at most one candidate is emitted and
verified. Submitted `time_log2` is 128.295287, peak memory `2^134.071246` bytes,
preprocessing `2^127.271374` target-group-operation equivalents (included in
total work), advice 0 bytes.

This is an analytic generic accounting control. It is not a practical break, a
best-known attack estimate or a cryptanalytic advance, and no full-size run is
claimed.
