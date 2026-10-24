# Commuting-square Frobenius-class rho

Fixed-G1 pairing inversion. Compute γ = e(P, Q0), walk on Frobenius orbits in ⟨γ⟩, and return Q = [a]Q0 once z = γ^a.

Square steps are Granger–Scott cyclotomic squares. They commute with Frobenius, so they do not canonicalise. Multiply steps, half of the steps under the limb heuristic, select the orbit index from eleven real limbs of one Fp2 slot and multiply by the matching conjugate of γ or z.

The score is log2 of total target-group-operation equivalents under `classical-group-work-v2` with C = 20000. Accounting and the orbit heuristic are in `proof.md`.
