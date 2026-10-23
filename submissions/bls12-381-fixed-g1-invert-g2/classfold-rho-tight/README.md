# classfold-rho-tight

Fixed-argument pairing inversion into G2 for BLS12-381 by solving `z = e(P,Q0)^a` with a distinguished-point rho
walk on the 12-element Frobenius classes of GT (`x ~ x^p`, which acts as `x -> x^lambda`, `lambda` of order 12).
One dense Fp12 multiplication plus a 4-multiplication Frobenius-covariant class selector per step, a look-back
window that removes fruitless short cycles, trails that end at distinguished points (density 2^-80), and recovery by
replaying the two colliding trails.

* Claim: `time_log2 = 125.90` (12590/100), success probability at least 0.39 under disclosed heuristics H1-H3.
* `proof.md`: the algorithm, the exact work ledger, the success argument, the trap lemma, memory, scope of the experiments.
* `experiment.py`: standard library only; modes `ledger`, `algebra`, `toy`, `sim` (see `submission.json`).
* Selected through `tracks/bls12-381-fixed-g1-invert-g2/selected-submission.txt`.

This package is `classfold-rho-tuned` (125.93) with two accountings tightened and the design unchanged: the success
target is the required 0.39 (main event 392/1000 at `S_max`) instead of 0.40, and the per-slot non-multiplication work
is itemised as hard per-item bounds totalling 2,865 instructions and charged at 2,900, instead of 3,000. The sampler,
the table, the store (full 48-byte class keys), the key representation, the walk constant `mu = 11/10` and the memory
arrays are exactly those of `classfold-rho-tuned`.
