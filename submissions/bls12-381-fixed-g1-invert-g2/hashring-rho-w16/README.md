# hashring-rho-w16

Fixed-argument pairing inversion into G2 for BLS12-381 by solving `z = e(P,Q0)^a` with a distinguished-point rho
walk on the 12-element Frobenius classes of GT (`x ~ x^p`, which acts as `x -> x^lambda`, `lambda` of order 12).
One dense Fp12 multiplication plus a 3-multiplication Frobenius-covariant class selector per step, a look-back
window that removes fruitless short cycles, trails that end at distinguished points (density 2^-80), and recovery by
replaying the two colliding trails.

* Claim: `time_log2 = 125.817` (125817/1000), success probability at least 0.39 under disclosed heuristics H1-H3.
* `proof.md`: the algorithm, the exact work ledger, the success argument, the trap lemma, memory, scope of the experiments.
* `experiment.py`: standard library only; modes `ledger`, `algebra`, `toy`, `sim` (see `submission.json`).
* Selected through `tracks/bls12-381-fixed-g1-invert-g2/selected-submission.txt`.

This package is the accepted `hashring-rho-sel1exact` (125.82) with one algorithmic change: the look-back window
shrinks from 24 to 16 class keys. The 24-slot scan costs 72 instructions per slot; 16 slots cost 48, saving 24
instructions per multiplication slot with no other row touched. The trap lemma is re-proved for the new threshold
(patterns of length 17 and above are undetected; their total is still negligible), the coalescence loss `L2` drops
from `2*24*nu` to `2*16*nu`, and the claim is stated at thousandths to capture the gain. The walk, the
3-multiplication selector, the `2^-80` distinguished-point rule, the budget `S_max`, the `mu = 11/10` allowance,
the 9-instruction explicit margin (kept), the `2^47`-record store, replay-based recovery and every additive item
are unchanged in structure; only the window size, the unit, and the claim precision move.
