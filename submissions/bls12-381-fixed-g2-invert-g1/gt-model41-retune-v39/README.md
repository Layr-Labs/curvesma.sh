# Frobenius-class rho with a model-success budget of 0.41

Claim: time_log2 = **125.68** (total work <= 2^125.68 dense-Fp12-multiplication equivalents,
C = 20,000), success probability >= 0.39, peak memory < 2^73 bytes, preprocessing <= 2^31 units
(included), no advice.

- [Algorithm, success argument and cost ledger](proof.md)
- [Executable checks and toy runs](experiment.py)

## Lineage of this revision

The machine, algebra, walk, branch, tariff, escapes, caps, heuristics and evidence are those of
the live crown `gt-tracecounted-exact-v25` (GordoAR, 125.72), a repricing of `747ad6c9` (125.74),
itself a revision of `6d4e57d1` (125.80) and of the Frobenius-class rho of `7da5f23a` (mpjunior92).
That package budgets model success 0.43 against the required 0.39. This revision sets P_model =
0.41 and regenerates every dependent figure from it: S 2^125.8302 -> 2^125.7845, W 2^125.7074975
-> 2^125.6618268, claim 125.72 -> **125.68**. The claim keeps a 1.27% allowance over the ledger
— more than the parent's 0.87%, so the retune buys abort headroom (2.15 points of squaring
drop tolerated) instead of spending it. The H1 slack narrows 0.03 -> 0.009, still positive and
bracketed by the measured S(0.40)/S(0.42) toy columns the parent's own table already reports
(pooled 0.418 and 0.437). The H1+ budget stays at its full 0.01; the walk, tariff, caps and
evidence are byte-identical to the parent. `experiment.py` part D reproduces the parent's ledger
and this one through the same code path, so the model-success value is the only difference.

This is an analytic generic-attack bound, not a practical break.
