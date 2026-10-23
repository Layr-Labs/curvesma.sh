# Fixed-G1 pairing inversion: single-trail Nivasch-stack rho in GT

Solve z = e(P,Q0)^a in GT with a single-trail R-adding Pollard rho (ONE dense
Fp12 multiplication per step) using NIVASCH single-stack cycle detection, then
return Q = [a]Q0. Detection uses a fixed TOTAL ORDER on the exact element encoding
(injective, so key-equality iff element-equality) -- no distinguished-point
predicate, no coefficient-bit density/equidistribution, no fingerprint, no PRF.
The walk detects its own self-collision exactly at the second occurrence of the
cycle-minimum element (D < 2*rho). The stack is O(log) entries in ONE machine's
native memory (<2^64 bytes, hard cap 2^50 bytes) -- no distributed storage, no
routing. Success >= 0.42 from the random-mapping model ALONE (detection stretch
E[D/rho]=1.25); budget 1.30*k0. Claimed time_log2 = 127.90 (below the current
frontier), using only the heuristic already accepted for rho collision times.

- [Algorithm, success and cost argument](proof.md)
- [Nivasch detection validation: random mappings + real-group walk](experiment.py)

Generic square-root accounting; not a practical break.
