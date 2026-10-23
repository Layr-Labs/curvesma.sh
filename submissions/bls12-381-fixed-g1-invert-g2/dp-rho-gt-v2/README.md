# Fixed-G1 pairing inversion: parallel distinguished-point rho in GT

Solve z = e(P,Q0)^a in GT with a van Oorschot-Wiener parallel distinguished-point
Pollard rho (ONE dense Fp12 multiplication per step -- no tortoise/hare factor),
then return Q = [a]Q0. Distinguished points are routed to 2^16 separate storage
nodes (each < 2^64 bytes, routing charged); the partition and DP predicate read
the reference multiplication's own reduced tower representation (no export);
exponents are sampled from the rules' fresh uniform algorithm coins. Hard step
budget 1.3*k0; success >= 0.39. Claimed time_log2 = 127.85.

- [Algorithm, success and cost argument](proof.md)
- [Toy parallel-DP experiment](experiment.py)

Generic square-root accounting; not a practical break.
