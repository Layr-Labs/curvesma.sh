# unroll512-a32-anch16384-km330k

P-256 sigma rho, declared **126.728713** = 126728713/1000000 at C=1000. 512-copy
backward jump-table block, direct-branch forward block; anchor period 16384, tail
U=16334 (16384 mod l < 50 for all l<=64); 2P+64=32832 persistent states;
CAP_PERIODS=2^52. PLAIN_FIXED 605,065,449; C_COMMON 575.044965; reserve 0.002694847;
PER_STEP at most 589.258353.

Against the inherited 126.728714 package this copy changes only the free
residual-Markov factor kM (375000 -> 330000) and replaces the 2.5e-7 floating-point
success guard with an exact-rational certificate (margin 2.422e-8, floor 1e-9).
See proof.md, submission.json, accounting.py.
