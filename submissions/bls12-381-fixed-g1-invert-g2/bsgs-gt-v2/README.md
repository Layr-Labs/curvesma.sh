# Fixed-G1 pairing inversion: deterministic sort-based BSGS in GT (v2)

Solve z = e(P,Q0)^a in GT with deterministic baby-step giant-step whose collision
is found by RADIX SORT + adjacent scan (no hash table), then return Q = [a]Q0.
Success probability exactly 1; hard assumption-free work bound. Claimed
time_log2 = 128.60. Exponential memory (reported).

- [Algorithm and hard cost bound](proof.md)
- [Toy experiment: sort-based deterministic recovery](experiment.py)

Generic deterministic control; not a best-known attack or a practical break.
