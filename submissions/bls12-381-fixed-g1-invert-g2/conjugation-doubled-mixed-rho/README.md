# Conjugation-doubled mixed-map rho on GT

Selected candidate for `bls12-381-fixed-g1-invert-g2`. Given nonidentity `P` in
`G1` and `z` in `GT`, the attack computes `g = e(P,Q0)` with one pairing, runs
`2^20` hash-driven walks in `GT` whose steps are cyclotomic squarings, Frobenius
maps, and rare table additions, emits each state together with its free
conjugate, and solves the first collision as a linear equation in the discrete
logarithm. It returns `Q = [a]Q0` after checking `e(P,Q) = z`.

* Claimed `time_log2`: **125.73** (target-group-operation equivalents).
* Peak memory: at most `2^138.3` bytes (birthday table of `2^127.9205` elements).
* Preprocessing: at most `2^20.1`, included in the total. Advice: none.

The success argument counts distinct linear conditions rather than element
pairs; the rare additions are what prevent the ratio collapse that a pure
power-map orbit suffers. Details, prices and limitations are in `proof.md`.
`experiment.py` checks only the collision algebra on a toy prime-order group.
