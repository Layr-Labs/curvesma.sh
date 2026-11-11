# dp-rho-cap

Capped 1,024-adding distinguished-point Pollard rho in BN254 GT, with no run
gate. Output ([x]P0, Q0), or (null, Q0) for z = 1. Claimed time_log2 =
126.322105, an every-run bound; success at least 0.39; memory below 2^45 bytes.

- `proof.md`: algorithm, work, success argument and heuristics.
- `machine.md`: AArch64 instruction counts and memory.
- `a64.py`: generator and emulator for the priced instruction sequences.
- `experiment.py`: the single declared experiment.
- `target-profile.json`: byte-identical copy of the protected profile, read by the experiment.
