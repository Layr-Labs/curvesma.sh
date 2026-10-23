# Fixed-G1 pairing inversion: deterministic BSGS in GT

Solve z = e(P,Q0)^a in GT with deterministic baby-step giant-step, then return
Q = [a]Q0. Success probability exactly 1 (no heuristic, no cycles, no success
analysis). Claimed time_log2 = 128.47 dense-Fp12-multiply equivalents; hard
deterministic work bound 2*sqrt(r) operations. Exponential memory (reported).

A clean deterministic accounting control for the first frontier on this track;
rho/Frobenius (~126.3) and TNFS are the directions to beat it.

- [Algorithm and hard cost bound](proof.md)
- [Toy experiment: deterministic recovery, 2*sqrt(l) ops](experiment.py)

Generic deterministic control; not a best-known attack or a practical break.
