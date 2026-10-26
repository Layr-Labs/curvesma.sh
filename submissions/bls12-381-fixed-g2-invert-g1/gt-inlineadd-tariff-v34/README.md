# Frobenius-class rho with an inline-add tariff: 125.68

Claim: time_log2 = **125.68** (total work <= 2^125.68 dense-Fp12-multiplication equivalents,
C = 20,000), success probability >= 0.39, peak memory < 2^73 bytes, preprocessing <= 2^31 units
(included), no advice.

- [Algorithm, success argument and cost ledger](proof.md)
- [Executable checks and toy runs](experiment.py)

## Lineage of this revision

The machine, algebra, walk, branch mix (s = 3/8), model point (0.43), escapes, caps, heuristics,
evidence, report multimap and memory layout are those of the live gate package
`gt-tracecounted-exact-v25` (GordoAR, submission `56fc259f`, 125.72), a revision of `747ad6c9`
(125.74), `6d4e57d1` (125.80) and of the Frobenius-class rho of `7da5f23a` (mpjunior92). The walk
dynamics are byte-identical, so every collision statistic carries over.

The parent prices every 6-limb modular add/sub/neg at 45 instructions "as a called routine" —
including 9 instructions of call, return and address set-up per add — even inside the
straight-line walk kernel where the adds are necessarily inlined. This revision prices inlined
adds at their counted 36 instructions (6 add/sub, 6 conditional-subtract, 6 select, 12 loads, 6
stores), keeps the measured blst_fp_mul = 383 and every other line of the parent's listing
unchanged, and regenerates the ledger: Fp2 squaring 874, Fp2 addition 72, Granger-Scott 10,514
instructions (0.5257 units), multiply-extra 2,038, common 270. Step prices: multiply 1.1154,
squaring 0.5392, average c = 0.899325 (= 35973/40000); S = 2^125.8302,
W = 2^125.6780227, claim **125.68**. `experiment.py` part D recomputes the ledger from the single
pricing function; section 3.4 re-derives the abort headroom.

This is an analytic generic-attack bound, not a practical break.
