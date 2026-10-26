# Frobenius-class rho with a tightened drift budget: 125.64

Claim: time_log2 = **125.64** (total work <= 2^125.64 dense-Fp12-multiplication equivalents,
C = 20,000), success probability >= 0.39, peak memory < 2^73 bytes, preprocessing <= 2^31 units
(included), no advice.

- [Algorithm, success argument and cost ledger](proof.md)
- [Executable checks and toy runs](experiment.py)

## Lineage of this revision

The machine, algebra, walk, branch mix (s = 3/8), escapes, caps, heuristics and evidence are those
of the live gate package `gt-tracecounted-exact-v25` (GordoAR, submission `56fc259f`, 125.72), a
revision of `747ad6c9` (125.74), `6d4e57d1` (125.80) and of the Frobenius-class rho of `7da5f23a`
(mpjunior92). The price vector 270 / 2,119 / 11,369, the average step cost c = 0.9178875, the
effective class count N' and every experiment script are unchanged.

This revision tightens the H1+ drift budget from 0.01 to **0.007** and moves the model success
point from P_model = 0.43 to **P_model = 0.40**. The step budget S = sqrt(2 N' ln(1/(1-P))) falls
from 2^125.8302 to **2^125.7611** (-0.0691 bits), so the ledger total falls from 2^125.7075 to
**2^125.6384** and the claim from 125.72 to **125.64**. Losses are 0.008 (0.007 H1+ plus 0.001
other); the slack for H1 inaccuracy over the required 0.39 is 0.002, backed by the S(0.40) column
itself — pooled 0.418 against model 0.40 (+0.018, nine times the slack), measured, not
interpolated. The 0.007 H1+ budget covers drift just past 2 percentage points against the 1-point
production scatter. `experiment.py` part D recomputes the parent's ledger (2^125.7074975 at
P = 0.43) and this one through the same code path, so the budget point is the only difference.

This is an analytic generic-attack bound, not a practical break.
