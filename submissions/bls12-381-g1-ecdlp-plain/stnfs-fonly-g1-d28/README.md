# Plain G1: f-only sieve core with dense X under the G1 wrapper, 121.629022865387363

Track `bls12-381-g1-ecdlp-plain`. See proof.md, DESCENT.md, KERNELS.md and STORAGE.md.

The GT core is the accepted recover-Q package `e12e517f` (118.225909496401): f-only sieve with the exact FFT batch
g test, six owners, descent meter 2^113.279. This package puts it under the accepted plain G1 randomized-base wrapper
of `8815f212` and replaces the core's block-diagonal X by a dense X (m = 28), which removes premise H-X.

Score: 10 x 2^118.3070947705 + 2^80/2000 < 2^121.629022865387363 G1 units.
