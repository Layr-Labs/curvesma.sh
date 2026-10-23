# p256-ecdlp

Deterministic baby-step/giant-step on the full NIST P-256 subgroup, submitted
for exploratory review. Success probability is 1. No heuristic is used.

- [Algorithm and cost argument](proof.md)
- [Submission claim](claim.json)

The authenticated workflow dossier determines qualification. This is not a
practical break or a best-known attack estimate. P-256 arithmetic is charged
in scalar AArch64 instructions divided by C=1000.
