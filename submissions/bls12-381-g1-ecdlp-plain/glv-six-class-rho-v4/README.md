# Six-class (negation + GLV) parallel Pollard rho in BLS12-381 G1

Parallel distinguished-point Pollard rho on the (r-1)/6 orbits of G1 under
<-1, phi>, phi(x, y) = (beta x, y). Class key |y| (free); canonical member via
one Fp multiplication; a twisted table makes each step one batched affine G1
addition plus at most 622 instructions; fruitless cycles are escaped by
doubling. A global work meter charges every action and halts at the claimed
bound, so W <= 2^126.64 on every run.

Claim: time_log2 = 126.64, success >= 0.40, memory < 2^72 bytes.

- [Algorithm, success argument and cost](proof.md)
- [Experiment: BLS12-381 identity checks + toy-curve runs](experiment.py)
