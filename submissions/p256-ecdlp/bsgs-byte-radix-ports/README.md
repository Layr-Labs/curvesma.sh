# p256-ecdlp: byte-radix odd-lift BSGS

Deterministic odd-lift baby-step/giant-step on NIST P-256. Success probability
is 1 on the full production target. Each 16 KiB page holds 32 records and
orders them with a stable byte-radix shuffle through 256 ports of capacity 4.
Sorted pages then meet in a 16-ary local-memory merge. No page address is a
scaled 128-bit index.

- [Coverage and instruction ledger](proof.md)
- [Claim record](claim.json)

This is analytic accounting, not a practical break.
