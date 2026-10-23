# plain-dp-rho-v1: parallel distinguished-point rho in GT

Given nonidentity `P` in G1 and `z` in GT, compute `g = e(P,Q0)` once, solve
`x = log_g z` with a parallel distinguished-point Pollard rho in the full group
GT, return `Q = [x]Q0` and verify with a pairing.

Hard cap `W_max = 1.5 sqrt(r) = 2^128.013507` steps at `1.01`
target-operation equivalents per step (one dense Fp12 multiplication plus 200
charged scalar instructions), plus a flat `2^31` reference-operation overhead for
the setup pairing, the `2^20`-multiplier table, the `2^20` walk starts and
verification. Claimed `log2(W) = 128.027863`.

Records store the exact 576-byte Fp12 element (no fingerprints, so no false
matches), distinguished-point density `theta = 2^-60`, `2^20` walks in lockstep,
hard storage capacity `2^70` records (`<= 2^79` bytes) enforced by a
stop-inserting policy. Success probability claimed: `0.39` (the rule threshold;
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
