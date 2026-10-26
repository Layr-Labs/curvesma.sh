# Frobenius-class rho with an exactly-priced walk kernel: 125.69

Claim: time_log2 = **125.69** (total work <= 2^125.69 dense-Fp12-multiplication equivalents,
C = 20,000), success probability >= 0.39, peak memory < 2^73 bytes, preprocessing <= 2^31 units
(included), no advice.

- [Algorithm, success argument and cost ledger](proof.md)
- [Executable checks and toy runs](experiment.py)

## Lineage of this revision

The walk, branch mix (s = 3/8), model point (0.43), escapes, caps and heuristics are those of the
live gate package `gt-tracecounted-exact-v25` (GordoAR, submission `56fc259f`, 125.72), a revision
of `747ad6c9` (125.74), `6d4e57d1` (125.80) and of the Frobenius-class rho of `7da5f23a`
(mpjunior92). The walk dynamics are byte-identical, so every collision statistic carries over.

What changes is the pricing integrity of the kernel. The parent prices every 6-limb modular
add/sub/neg at 45 instructions "as a called routine" — including 9 instructions of call, return
and address set-up per add — even inside the straight-line walk kernel where the adds are
necessarily inlined. This revision prices inlined adds at their counted 36 instructions (6
add/sub, 6 conditional-subtract, 6 select, 12 loads, 6 stores), keeps the measured blst_fp_mul =
383 for the multiplications, and regenerates the ledger: Fp2 squaring 901 -> 874, Fp2 addition 90
-> 72, Granger-Scott 11,369 -> 10,514, multiply-extra 2,119 -> 2,158 (see below), squaring step
0.58195 -> 0.5392, multiply step 1.11945 -> 1.1214, average step cost c 0.9178875 -> **0.903075**
(= 36123/40000), W 2^125.7075 -> 2^125.6841, claim 125.72 -> **125.69**. Two correctness items
ride with the repricing so the meter is exact on every run: the twelve-key minimum is the full
stored-key minimum with a specified lower-limb fallback charged at its 270-instruction worst case
on every multiply step (the +120 over the parent's top-limb line is inside the 2,158), and the
chained multimap is replaced by append-plus-sort matching with a fully bounded work account in
the 2^64.3 non-walk bucket. `experiment.py` part D reproduces the parent's ledger and this one
through the same code path, and the new part F verifies sort matching finds exactly the
dict-matching collisions on toy reports.

This is an analytic generic-attack bound, not a practical break.
