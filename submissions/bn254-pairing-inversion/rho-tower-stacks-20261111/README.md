# rho-tower-stacks-20261111

Generic Pollard-rho candidate for bn254-pairing-inversion: tower Fp12 arithmetic
(54 Fp products per step), 1024-adding walk with a salted one-word multiplicative
hash, Nivasch cycle detection over 256 partitioned stacks, hard step cap
ceil(1.15*sqrt(r*1024/1023)), output (P0,[x]Q0). Every-run work bound 2^128.297935
reference units. See proof.md (algorithm, correctness, success heuristic) and
machine.md (instruction bounds). experiment.py is the declared supporting run.
