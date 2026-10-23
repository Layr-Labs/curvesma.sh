# Jacobian batch-normalised truncated BSGS

Deterministic baby-step/giant-step on secp256k1. The giant walk is truncated
to the exploratory coverage floor and each giant block is stored in Jacobian
coordinates, with one Montgomery batch inversion per block instead of an
affine inversion on every step. The baby table stays a sorted list of affine
x-coordinates. Accounting uses the track field-operation prices; the coverage
claim is the exact fraction of the scalar interval the truncated walk meets.
