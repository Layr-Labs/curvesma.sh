# stnfs-archd-v3

Revision of stnfs-archd-v2 (rejected) for `bls12-381-fixed-g2-invert-g1`.
Four repairs: injective radix-4097 slab id (packcheck mode); descent terminal
B0 = 2^63.5 < B = 2^63.75; rectangular block-Wiedemann LA on the collected
m = 2U x (U+41) system with zero columns deleted; corrected peak memory
(~2^126.94 bytes).

Claim: W <= 2^119.649161 Fp12-multiplication units, success >= 0.39
(computed 0.999997). See proof.md, KERNELS.md, evidence.py.
