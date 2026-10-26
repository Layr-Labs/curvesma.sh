# Frobenius-class rho with a half-square branch mix and exact canonical keys

Claim: time_log2 = **125.71** (total work <= 2^125.71 dense-Fp12-multiplication equivalents,
C = 20,000), success probability >= 0.39, peak memory < 2^73 bytes, preprocessing <= 2^31 units
(included), no advice.

- [Algorithm, success argument and cost ledger](proof.md)
- [Executable checks and toy runs](experiment.py)

## Lineage of this revision

The machine, algebra, walk, branch, escapes, caps, heuristics and evidence are those of the live
gate package `gt-tracecounted-exact-v25` (GordoAR, submission `56fc259f`, 125.72), a revision of
`747ad6c9` (125.74), `6d4e57d1` (125.80) and of the Frobenius-class rho of `7da5f23a` (mpjunior92).
This revision moves the squaring mix from 3/8 to 1/2 (threshold 4 on the low 3 bits of T, i.e. bit
2) and regenerates every dependent figure from the widened price vector 270 / 2,239 / 11,369
(the key line grows by 120 instructions so the exact full-limb minimum is priced on every
multiply step): average step cost c 0.9178875 -> 0.8537 (= 8537/10000), S = 2^125.9283,
W 2^125.7075 -> 2^125.7010, claim 125.72 -> **125.71**. The allowance over the ledger shrinks
from 0.87% to 0.62%; proof.md section 3.4 re-derives the abort headroom and states the drift-cost table.
The squaring test is stated exactly as the code tests it ((T & 7) >= 4), part C runs the
production rule at the production table size (R = 1024), and kappa is the exact minimum of the
full stored key products with a fully specified, priced lower-limb fallback — so the meter's
every-run work guarantee holds with no tie estimate anywhere in the path.
`experiment.py` part D reproduces the parent's ledger (2^125.7074975 at s = 3/8, key charge 2,119)
and this one through the same code path.

This is an analytic generic-attack bound, not a practical break.
