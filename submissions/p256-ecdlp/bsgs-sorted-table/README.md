# p256-ecdlp: sorted-table BSGS

Deterministic baby-step/giant-step accounting control on NIST P-256, submitted for
exploratory review. Success probability is 1 on the specified full target.

- [Algorithm and cost argument](proof.md)
- [Submission claim](claim.json)

The table is a merge-sorted array of affine records with a fixed 129-comparison
binary search; there is no hash function, probe sequence or occupancy bound in the
construction. This is not a practical break, a best-known attack estimate or a
cryptanalytic advance.
