# Limb-arithmetic accounting control for BLS12-381 G1 DLOG

Deterministic baby-step/giant-step recovery of `s` from `T = [s]P0` on BLS12-381
G1, charged on a scalar 64-bit instruction basis with the fixed G1 reference
price `C = 2000` (batched incomplete affine 6M addition). The bound is an
analytic accounting control, not a practical execution; see `proof.md`.
