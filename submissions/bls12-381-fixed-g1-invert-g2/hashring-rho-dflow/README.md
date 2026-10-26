# hashring-rho-dflow

Fixed-argument pairing inversion into G2 for BLS12-381 by solving `z = e(P,Q0)^a` with a distinguished-point rho
walk on the 12-element Frobenius classes of GT (`x ~ x^p`, which acts as `x -> x^lambda`, `lambda` of order 12).
One dense Fp12 multiplication plus a 3-multiplication Frobenius-covariant class selector per step, a look-back
window that removes fruitless short cycles, trails that end at distinguished points (density 2^-80), and recovery by
replaying the two colliding trails.

* Claim: `time_log2 = 125.819` (125819/1000), success probability at least 0.39 under disclosed heuristics H1-H3.
* `proof.md`: the algorithm, the exact work ledger, the success argument, the trap lemma, memory, scope of the experiments.
* `experiment.py`: standard library only; modes `ledger`, `algebra`, `toy`, `sim` (see `submission.json`).
* Selected through `tracks/bls12-381-fixed-g1-invert-g2/selected-submission.txt`.

This package is the accepted `hashring-rho-sel1exact` (125.82) with one correctness fix and two carriers. The fix
closes the composed-loop dataflow gap the review of the sibling no-margin line identified: no listed operation
transferred the fresh coefficient `a` (and `b`) from the dense-kernel output buffer into the selector's `V` array,
while `selector_tail_prog` and block (i) loaded their operands from `V`. The repair repoints those existing loads
to the product buffer through the `y'` base register `x29` — the same loads, different immediates, zero added
instructions — and the tests now populate a buffer model instead of prepopulating `V`. The carriers: the claim is
stated at thousandths (the fix costs nothing, so the total is the parent's and the least thousandth above it is
125.819), and `proof.md` is lightly trimmed to fit the 65,536-byte intake cap. The walk, the 3-multiplication
selector, the 24-class hash-ring window, the `2^-80` distinguished-point rule, the budget `S_max`, the `mu = 11/10`
allowance, the per-slot charge of 1,749 (9-instruction margin kept), the `2^47`-record store, replay-based recovery
and every additive item are otherwise unchanged.
