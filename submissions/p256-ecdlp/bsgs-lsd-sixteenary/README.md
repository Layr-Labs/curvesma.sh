# p256-ecdlp

Deterministic odd-lift BSGS. Baby x-coordinates are ordered by a global LSD
byte distribution. Lookup is an implicit sixteen-ary descent on the sorted
array. No heap, no per-bucket directory and no port queue.

- [Algorithm and instruction bound](proof.md)
- [Claim](claim.json)

Success probability is 1 on the profile distribution. Submitted time_log2 is
1329/10. The toy experiment does not measure that bound.
