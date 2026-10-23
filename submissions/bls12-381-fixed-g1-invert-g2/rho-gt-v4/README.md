# Fixed-G1 pairing inversion: bounded single-processor rho in GT (v4)

Solve z = e(P,Q0)^a in GT with a strictly-bounded single-processor (serial van
Oorschot-Wiener) distinguished-point Pollard rho, then return Q = [a]Q0. Hard
global work cutoff W_max = ceil(1.25*k0); 384-bit fingerprints with resolution;
exceptional events (truncation, detection lag, false match) bounded and allowed.
Claimed time_log2 = 127.80; success >= 0.39 within W_max under the random-mapping
heuristic (the single assumption).

- [Algorithm, success and cost argument](proof.md)
- [Strictly-bounded toy experiment](experiment.py)

Generic square-root accounting; not a practical break.
