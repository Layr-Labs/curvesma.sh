# Fixed-G1 pairing inversion: deterministic exact-element sort BSGS in GT (v3)

Solve z = e(P,Q0)^a in GT with deterministic baby-step giant-step whose collision
is found by RADIX SORT on the exact element + adjacent scan (no hash, no fingerprint), then return Q = [a]Q0.
Success probability exactly 1; hard assumption-free work bound. Claimed
time_log2 = 129.25. Exponential memory (reported).

- [Algorithm and hard cost bound](proof.md)
- [Toy experiment: sort-based deterministic recovery](experiment.py)

Generic deterministic control; not a best-known attack or a practical break.
