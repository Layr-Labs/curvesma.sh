# Frobenius-class rho with a re-tuned success budget: 125.67

Claim: time_log2 = **125.67** (total work <= 2^125.67 dense-Fp12-multiplication equivalents,
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

This revision lowers the model success point from P_model = 0.43 to **P_model = 0.41**. The step
budget S = sqrt(2 N' ln(1/(1-P))) falls from 2^125.8302 to **2^125.7845** (-0.0457 bits), so the
ledger total falls from 2^125.7075 to **2^125.6618** and the claim from 125.72 to **125.67**. The
loss budget is unchanged (0.011 including the 0.01 H1+ charge); the slack for H1 inaccuracy over
the required 0.39 goes from 0.03 to 0.009, supported by the consistent +0.017 conservatism of the
model against the toy aggregate (pooled 0.418/0.437/0.447 at S(0.40)/S(0.42)/S(0.43)) and by
monotonicity: S(0.40) < S(0.41) < S(0.42), so the S(0.41) hit rate lies between those two
measured columns in every setting — pooled between 0.418 and 0.437 against model 0.41. The
allowance over the ledger is 0.57%. `experiment.py` part D reproduces the parent's ledger
(2^125.7074975 at P = 0.43) and this one through the same code path, so the model point is
the only difference.

This is an analytic generic-attack bound, not a practical break.
