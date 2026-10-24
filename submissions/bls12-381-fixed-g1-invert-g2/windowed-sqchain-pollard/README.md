# Windowed squaring-backbone Pollard walk

Fixed-G1 pairing inversion. Compute γ = e(P, Q0), walk in ⟨γ⟩ with a shared
cyclotomic squaring backbone, and return Q = [a]Q0 once z = γ^a.

The score is log2 of total target-group-operation equivalents under
`classical-group-work-v2` with C = 20000. The accounting and the
random-image hypothesis are in `proof.md`. The toy walker in `experiment.py`
checks only collision algebra on a prime-order group.
