# rho-auto-v3

Automorphism-aware Pollard rho on BLS12-381 G1, revision 3. The walk is on the
six-element orbits of negation and the GLV endomorphism. Each step is one
batched incomplete affine addition (the track reference unit, C = 2000) plus
935 instructions: two base-field multiplications by the public cube roots of
unity, the finite-batch inversion excess at batch 4096, and the limb
comparisons that pick the canonical representative.

`proof.md` derives the step ceiling and the score
127016229/1000000 = log2(W) in target-group-operation equivalents.
`experiment.py` recomputes that ceiling and runs the same canonical walk on
prime orders 10009, 100003, 1000003 and 5000011, plus the order-4093 subgroup
of y^2 = x^3 + 4 over F_12253. `claim.json` binds the protected profile and
cost-model hashes.
