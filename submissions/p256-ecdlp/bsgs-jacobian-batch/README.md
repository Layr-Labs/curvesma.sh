# p256-ecdlp: BSGS with projective mixed additions and block batch inversion

Recover x from Q=[x]P on NIST P-256 (secp256r1, a=-3).

Deterministic baby-step giant-step over [0,r) with m=2^128 and h=ceil(r/m).
Both walks run in Jacobian coordinates with mixed additions, so no inversion is
needed inside a loop; every 1024 steps the buffered Z values are batch-inverted
with one inversion and the block is converted to affine. The baby table is
open-addressed over 2^132 slots keyed by the full affine pair. Coverage is
exact, including the identity target, and the only probabilistic element is the
hash multiplier (failure at most 2^-128, charged at the 64-probe cap).

Files:

- `submission.json` - claim and the toy experiment schema.
- `claim.json` - profile/cost-model bindings and declared metrics.
- `proof.md` - construction, coverage proof, batching schedule, instruction
  ledger and exact integer totals.
- `experiment.py` - the same update on the real P-256 curve at toy size.

Declared: time_log2 132.078610, preprocessing_log2 131.105679,
memory_log2_bytes 138.321929, advice 0. This is a conservative generic
accounting control, not a practical attack; the exponential table is mandatory
for the construction.
