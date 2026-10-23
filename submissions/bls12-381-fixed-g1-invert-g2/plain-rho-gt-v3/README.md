# Fixed-G1 pairing inversion: bounded distinguished-point rho in GT (v3)

Solve z = e(P,Q0)^a in GT with a fixed-budget distinguished-point Pollard rho
(no automorphism), then return Q = [a]Q0. The attack has a hard global work
cutoff W_max = ceil(1.25*k0) and per-walk truncation; exceptional
distinguished-point-free trajectories are allowed and bounded (~2^-92). Claimed
time_log2 = 127.80; success >= 0.39 within W_max under the random-mapping
heuristic (the single assumption).

Supersedes plain-rho-gt (v1/v2): drops the false "no fruitless cycles" claim,
adds the explicit charged-work cutoff, and analyses success allowing exceptional
unsuccessful trajectories.

- [Algorithm, success and cost argument](proof.md)
- [Toy bounded experiment](experiment.py)

Generic square-root accounting; not a practical break.
