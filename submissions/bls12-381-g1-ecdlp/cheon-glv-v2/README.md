BLS12-381 G1 auxiliary-powers track: a deterministic Cheon attack priced with
GLV scalar multiplication and explicit batching/table accounting.

The protected rules supply A_i=[s^i]P0 for i=1..n with scoring n=2^20 and
state that d=2^20 divides r-1 and that Cheon's attack using A_d is admissible.
This package uses A_1=T and A_d only. Phase 1 solves the order-M discrete
logarithm s^d=zeta_hat^(k0) by baby-step giant-step over [zeta_hat^k0]P0;
phase 2 solves the order-d logarithm s*zeta^(-k0)=zeta_check^k1. Each scalar
multiplication is a GLV multiplication on the production curve with the
verified endomorphism pair (beta,lambda); half-scalars are below 2^127.

Accounting: additions are the protected batched incomplete affine reference
(1 unit) in batches of 2^64 live chains with the 240-byte/chain workspace
counted; doublings are 7M/2000 = 1.3405 units; a 16-unit per-multiplication
overhead is itemised at <=3584 instructions; the baby table is a sorted array
of m 112-byte records (96-byte point + 16-byte index) with merge sort charged
5 units/entry over 118 passes and binary search 4 units/lookup. Total:
(m+m_hat)*320.265 + m*5 + m_hat*4 + 2048*320.265 + small, with
m=m_hat=223621884914696853578958571250384896, gives 2^126.7717970576 and the
exact score 63385899/500000 = 126.771798.

Success probability is 1 (deterministic). Peak memory is 2^125.2359 bytes
(two m-record sort buffers plus batch workspace), preprocessing is
2^125.774016 units (included), advice zero. Every stated exponent bound is at
least the ledger's exact value.

Read proof.md for the attack, the accounting and the limitations.
experiment.py recomputes the ledger and production constants and runs the
two-phase attack end-to-end on toy groups.
