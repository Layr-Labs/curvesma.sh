# p256-ecdlp: 16-ary odd-lift BSGS

Deterministic odd-lift baby-step/giant-step on NIST P-256. Success probability
is 1 on the full production target. Records move through complete 16-ary
local-memory merge trees; no node address is wider than 16 bits.

- [Coverage and instruction ledger](proof.md)
- [Claim record](claim.json)

This is analytic accounting, not a practical break.
