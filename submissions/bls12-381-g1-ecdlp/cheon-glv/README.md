BLS12-381 G1 auxiliary-powers track: a deterministic Cheon attack priced with
GLV scalar multiplication.

The protected rules supply A_i=[s^i]P0 for i=1..n with scoring n=2^20 and
state that d=2^20 divides r-1 and that Cheon's attack using A_d is admissible.
This package uses A_1=T and A_d only.  Phase 1 solves the order-M discrete
logarithm s^d=zeta_hat^(k0) by baby-step giant-step over the geometric
progression [zeta_hat^k0]P0; phase 2 solves the order-d logarithm
s*zeta^(-k0)=zeta_check^k1.  Each scalar multiplication is a GLV
multiplication on the production curve, with the endomorphism pair
(beta,lambda) verified on P0 and half-scalars below 2^127.

Ledger: (m+m_hat)*321.265 + 2048*320.265 + 8 + 2^40/2000 units, with
m=m_hat=223621884914696853578958571250384896, gives
2^126.7561647122 and the exact score 25351233/200000 = 126.756165.
Success probability is 1 (deterministic).  Memory is 2^124.2359 bytes,
preprocessing 2^125.7562 units (included), advice zero.

Read proof.md for the attack, the cost conventions and the limitations.
experiment.py recomputes the ledger and production constants and runs the
two-phase attack end-to-end on toy groups.
