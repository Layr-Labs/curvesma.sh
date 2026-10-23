# adding-walk

Floyd cycle-finding Teske r-adding walk (R = 32 jumps) on BLS12-381 G1 with
canonical affine state. `proof.md` specifies the attack and derives the
integer instruction bound W = 514,694,484,759,148,014,223,122,347,276,647,627,
415,912,096, scored as time_log2 = 137562769/1000000 target-group-operation
equivalents (C = 2000, classical-group-work-v2). `experiment.py` recomputes
the ledger, asserts the score is a ceiling of log2(W/2000), and runs the exact
walk function on toy groups (prime additive orders 10007/100003/1000003 and a
prime-order-3187 subgroup of a short Weierstrass curve over F_9613).
`claim.json` mirrors the scored fields. Field multiplication is charged at
the pinned 383-instruction `blst_fp_mul` measurement.
