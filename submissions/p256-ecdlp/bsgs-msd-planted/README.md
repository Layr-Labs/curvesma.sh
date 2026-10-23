# p256-ecdlp

Deterministic odd-lift BSGS. Baby x-coordinates are ordered by a cross-page
planted-base MSD byte distribution. No heap and no port queue.

- [Algorithm and instruction bound](proof.md)
- [Claim](claim.json)

Success probability is 1 on the profile distribution. Submitted time_log2 is 133.
The toy experiment does not measure that bound.
