# classfold-rho-sel3

Fixed-argument pairing inversion into G2 for BLS12-381 by solving `z = e(P,Q0)^a` with a distinguished-point rho
walk on the 12-element Frobenius classes of GT (`x ~ x^p`, which acts as `x -> x^lambda`, `lambda` of order 12).
One dense Fp12 multiplication plus a 3-multiplication Frobenius-covariant class selector per step, a look-back
window that removes fruitless short cycles, trails that end at distinguished points (density 2^-80), and recovery by
replaying the two colliding trails.

* Claim: `time_log2 = 125.86` (12586/100), success probability at least 0.39 under disclosed heuristics H1-H3.
* `proof.md`: the algorithm, the exact work ledger, the success argument, the trap lemma, memory, scope of the experiments.
* `experiment.py`: standard library only; modes `ledger`, `algebra`, `toy`, `sim` (see `submission.json`).
* Selected through `tracks/bls12-381-fixed-g1-invert-g2/selected-submission.txt`.

The design is the one in `classfold-rho-tight` unchanged: the same `2^14` twisted-multiplier table, the same
24-class look-back window on the same keys with the same escape rule, the same `2^-80` distinguished-point rule,
the same `2^47`-record store of full 48-byte class keys, the same replay-based recovery, and the same walk-quality
allowance `mu = 11/10` with its simulation evidence. What changed is one line of the selector's accounting, and it is
exact: in this tower the Frobenius constant is `gamma = g0*(1-u)` (so `g1 = -g0`, and equivalently `g0^6 = -1/8` and
`zeta = 2 g0^2`), hence the selector's second seed `d = Re(gamma*conj(c1)) = g0*a + g1*b` is `g0*(a-b)`, one
multiplication of the difference `a-b` instead of two. With `zeta*a` and `zeta*d` the twelve keys `{zeta^j a}` and
`{zeta^j d}` therefore cost 3 Fp multiplications, not 4, and the per-slot itemisation drops from 2,624 to 2,241
instructions (charged 2,253 with the same 0.5% margin), so each multiplication slot costs `22253/20000` instead of
`11319/10000`. The success target, the accepted-step budget `S_max`, the memory ledger and every other row of the
ledger are unchanged; only the total work moves, from `2^125.85144` against the earlier `2^125.87618`.

`dependencies.json` credits the AI-valid package whose design this is.
