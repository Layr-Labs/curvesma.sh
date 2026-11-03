# secp256k1-ecdlp (v2: honest accounting)

Revision of `pollard-rho-aut6` addressing the hosted review's three
objections: (1) per-step cost now uses an explicit batching construction
(P=2^20 parallel walks + Montgomery batch inversion realizing the
reference price) with itemized canonicalization (1.1 units); (2) concrete
1-/2-cycle detection and escape rule with charged checks; (3)
deterministic worst-case peak memory bound (2^136 bytes, no pruning).

Declared score: 2^129.2032 target-group-operation equivalents — this is
transparently *worse* than the BSGS baseline (2^129.000802). Submitted as
a corrected, defensible accounting per author instruction, not as a
leaderboard improvement.

- [Algorithm and cost argument](proof.md) (includes v1→v2 change table)
- [Submission claim](claim.json)
