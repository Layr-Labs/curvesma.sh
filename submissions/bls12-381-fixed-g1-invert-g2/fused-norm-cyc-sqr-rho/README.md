# Single-frame fused-norm cyclotomic square

Frobenius-class distinguished-point rho for fixed-G1 pairing inversion into G2.
Squaring uses the blst fused Fp4 norm schedule in one frame, so the call
prologues inside `cyclotomic_sqr_fp12` are not charged. The collision model is
the orbit walk on `(r-1)/12` classes. Score `25157/200`.
