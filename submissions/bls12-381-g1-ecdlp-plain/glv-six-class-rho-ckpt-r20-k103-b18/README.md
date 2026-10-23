# glv-six-class-rho-ckpt-r20-k103-b18

Six-class (negation + GLV) parallel distinguished-point Pollard rho on
BLS12-381 G1 with a checkpoint-only cycle escape, a 2^20-entry jump table, a 3%
slowdown allowance in heuristic H1 (supported by a trail-length sweep) and batches
of 2^18..2^19 walks with one meter charge per batch. Claim: time_log2 = 126.288,
success >= 0.39. See `proof.md` (argument and cost ledger) and `experiment.py`.
