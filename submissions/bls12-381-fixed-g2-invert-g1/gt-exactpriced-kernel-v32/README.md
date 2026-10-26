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

What changes is the pricing integrity of the kernel, stated as exactly one tariff. The parent
prices every 6-limb modular add/sub/neg at 45 instructions "as a called routine" — including 9
instructions of call, return and address set-up per add — even inside the straight-line walk
kernel where the adds are necessarily inlined. This revision prices inlined adds at their counted
36 instructions (6 add/sub, 6 conditional-subtract, 6 select, 12 loads, 6 stores), keeps the
measured blst_fp_mul = 383 for the multiplications, and regenerates the ledger: Fp2 squaring 874,
Fp2 addition 72, Granger-Scott **10,514** instructions (0.5257 units), multiply-extra **2,158**,
common 270. The resulting step prices are multiply 1.1214, squaring 0.5392, average
c = **0.903075** (= 36123/40000); S = 2^125.8302, W = 2^125.6840259, claim **125.69**. Two
correctness items ride with the repricing so the meter is exact on every run: the twelve-key
minimum is the full stored-key minimum with a specified lower-limb fallback charged at its
270-instruction worst case on every multiply step (inside the 2,158), and the chained multimap is
replaced by append-plus-sort matching — O(1) appends plus one final mergesort by fingerprint and
a neighbour-only scan — with the sort work charged in the 2^65.2 non-walk bucket. The parent's
superseded figures appear in this package only as labelled historical values that Part D
reproduces as a check (2^125.7074975); the claim asserts one tariff only. `experiment.py` part D
derives both ledgers from one pricing function, and the new part F runs the specified
neighbour-only scan on toy reports and checks it finds the same useful pairs as dict matching
within the verification cap.

This is an analytic generic-attack bound, not a practical break.
