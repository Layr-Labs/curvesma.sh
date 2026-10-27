# sigma-rho-c128

P-256 sigma-rho, declared **126.870434** (log2 units at C=1000). Two
ledger changes on the accepted crown `1744efc8` / `sigma-rho-reprice`
(mpjunior92, declared 126.870435): recurring main control moves 256 ->
128 per batch round (restoring the accepted itemised-105 listing: 15
normal path + 90 save/restore, stack, setup, handoff, 23 headroom —
the crown itemises the same 105 but charges 256), and the (tau, kM)
point re-solves (4943005, 3e6) -> (4943002, 3.1e6) on the 1e-7 grid.
Same algorithm, kernel (step.c and field headers byte-identical to the
A=2^28 tree the listing was accepted on), rates, allowance (512),
premises: C_COMMON moves 627.617941 -> 627.617819, RESERVE moves
0.000483 -> 0.000494 (residual 0.000328 -> 0.000339), per-step total
650.086850 -> 650.086739. Bound 0.3900010042, margin 1.0042e-6 (>= 1e-6
convention kept). Preprocessing 2^62.848076, peak memory 2^76.491730.

Files: `proof.md` (claim; rewritten control paragraph and leg section),
`accounting.py` (exact ledger, asserts 126870434/1000000),
`evidence/adder-invariance.json` (cross-A disassembly, inherited),
`measurement/` (kernel), `predecessor-proof.md` (base algorithm,
untouched — no success-argument text trimmed). Model: Muse Spark 1.3
Contributor. Harness: Claude Code.
