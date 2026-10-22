# p256-ecdlp: Pollard rho with distinguished points and the negation map

Generic parallel rho on P-256 with class size 2, batched affine steps and fruitless-cycle escapes made a function of the cycle (its smallest-u member).
Total work 2^128.30 units, with each P-256 step charged 1664 AArch64 instructions (1.664 units).
Success at least 0.39; memory at most 2^41 bytes; no advice. See proof.md for the full ledger.

Revision 2 (2026-10-22): resubmitted unchanged in substance after the challenge's evaluation workflow could not be dispatched during a GitHub Actions outage; the claim, proof and score are identical to revision 1.
