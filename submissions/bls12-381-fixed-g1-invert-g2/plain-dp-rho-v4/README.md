# plain-dp-rho-v1: parallel distinguished-point rho in GT

Given nonidentity `P` in G1 and `z` in GT, compute `g = e(P,Q0)` once, solve
`x = log_g z` with a parallel distinguished-point Pollard rho in the full group
GT, return `Q = [x]Q0` and verify with a pairing.

Hard cap `W_max = 1.5 sqrt(r) = 2^128.013507` steps at `1.01`
target-operation equivalents per step (one dense Fp12 multiplication plus 200
charged scalar instructions), plus an overhead charged by exact summation for the setup pairing, the
`2^20`-multiplier table, the `2^20` walk starts, verification and the capped
insert cost (total `< 2^67.001`). Claimed `log2(W) = 128.027863`.

Records store the exact 576-byte Fp12 element (no fingerprints, so no false
matches; bucket collisions are resolved by exact comparison). Two independent
hash slices are used: `H_1` drives the walk index and `H_2` the distinguished
test, so distinguished points are not confined to one bucket. Distinguished
density `theta = 2^-60`; the dictionary is open-addressed (2^71 slots of 16-byte record handles,
linear probing over <=4 slots, drop-on-fourth-miss) and the toy run measures
0.750 at the same relative budget with dictionary work hard-capped at 2^70
lookups/inserts = 2^66.678072 units, `2^20` walks in lockstep,
hard storage capacity `2^70` records enforced by a stop-inserting policy, with
every memory component added up to `<= 2^80` bytes (records 2^79.321928 B +
multiplier table 2^29.2 B + index/state/scratch < 2^40 B). Success probability claimed: `0.39` (the rule threshold;
the random-mapping model gives 0.675 at this cap, and the toy run measures 0.775
at the same relative budget). Preprocessing charged in full with no reuse credit;
nonuniform advice 0; no auxiliary powers.

No group quotient is used. An earlier two-class inversion-quotient draft of mine
was withdrawn after review correctly showed that its transition is not a function
on inversion classes; the `sqrt(2)` saving is not claimed here.

`proof.md` has the reduction, the search design, the instruction ledger, the
success analysis and the storage policy; `experiment.py` is standard-library only
(Part A: BLS12-381 group-order facts; Part B: bounded toy run of the complete
algorithm with the uncensored collision statistic).
