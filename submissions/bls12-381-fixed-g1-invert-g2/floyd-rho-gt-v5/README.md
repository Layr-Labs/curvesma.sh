# Fixed-G1 pairing inversion: O(1)-memory Floyd rho in GT (v5, budget 1.1)

Solve z = e(P,Q0)^a in GT with a deterministic-budget Pollard rho using Floyd
tortoise-hare cycle detection (tortoise one step, hare two; meet -> solve), then
return Q = [a]Q0. Working set is a ~2^30-byte adding table plus O(1) walk state,
all native-addressable (no large/segmented/distributed memory). Claimed
time_log2 = 129.20; success >= 0.39 within the hard step budget under the
random-mapping heuristic.

- [Algorithm, success and cost argument](proof.md)
- [Toy O(1)-memory experiment](experiment.py)

Generic square-root accounting; not a practical break.
