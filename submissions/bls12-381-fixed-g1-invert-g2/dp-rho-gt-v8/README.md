# Fixed-G1 pairing inversion: single-trail distinguished-point rho in GT

Solve z = e(P,Q0)^a in GT with a single-trail distinguished-point Pollard rho
(ONE dense Fp12 multiplication per step -- no tortoise/hare factor), then return
Q = [a]Q0. The partition and DP predicate read a bit-slice of the first tower
coefficient in the kernel's NATIVE Montgomery representation (mont(c)=c*2^384 mod p,
no export). The ~2^-64 DP density and routing balance (joint, and conditional on
the DP event) follow from a representation-independent subgroup exponential-sum
(BGK / Bourgain-Chang) equidistribution bound -- GT has order r | p^4-p^2+1 and,
by embedding degree 12, lies in no proper subfield -- and are confirmed at a
feasible 12-bit scale by the experiment on real GT elements. Distinguished points
are routed to 2^16 separate storage nodes (each < 2^64 bytes, routing charged);
exponents are sampled from the rules' fresh uniform algorithm coins. Hard step
budget 1.3*k0; success >= 0.39. Claimed time_log2 = 127.85.

- [Algorithm, success and cost argument](proof.md)
- [Montgomery-form DP-bit uniformity experiment on real GT](experiment.py)

Generic square-root accounting; not a practical break.
