# p256-ecdlp

Parallel Pollard rho on the r group points (no negation map), submitted for
exploratory review. Success probability is at least 0.64 under the disclosed
random-mapping heuristic M1 on the specified full target.

- [Algorithm and cost argument](proof.md)
- [Submission claim](claim.json)

The authenticated workflow dossier determines qualification. This is not a
practical break, best-known attack estimate or cryptanalytic advance.
P-256 arithmetic is charged in scalar AArch64 instructions divided by C=1000.
The negation map was deliberately declined; see proof.md for the reason.
