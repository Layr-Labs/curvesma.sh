# SexTNFS discrete logarithm in GT for fixed-G1 pairing inversion

Reduces e(P,Q)=z to the discrete logarithm z=e(P,Q0)^a in the order-r
subgroup of F_{p^12}^*. That logarithm is solved with the special extended
tower NFS, and the attack returns Q=[a]Q0. The claimed time_log2 is 120,
justified analytically in proof.md. This is a heuristic estimate and has no
production run.
