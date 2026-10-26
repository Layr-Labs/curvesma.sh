# Frobenius-class rho on GT with a re-tuned success budget

Claim: time_log2 = **125.68** (total work <= 2^125.68 dense-Fp12-multiplication equivalents,
C = 20,000), success probability >= 0.39 (model 0.415, 0.01 of it charged to H1+),
peak memory < 2^73 bytes, preprocessing <= 2^31 units (included in the total),
advice = 0 bytes.

- [Algorithm, success argument and cost ledger](proof.md)
- [Executable checks and toy runs](experiment.py)

## Lineage of this revision

The machine, algebra, walk, escapes, caps, price vector, heuristics and evidence are those
of the Frobenius-class rho rev 2 line at 125.72 (7da5f23a mpjunior92, 6d4e57d1, 747ad6c9
jungjipdo, 56fc259f GordoAR). The single change in this revision is the success budget:
the model success is re-tuned from 0.43 to 0.415, spending the same 0.01 on H1+ and keeping
0.015 absolute slack over the required 0.39 for the inaccuracy of H1. The walk, the counted
price vector (270 / 2119 / 11369), the table size R = 1024, the squaring mix s = 3/8 and the
experiment configuration are untouched, so every evidence table of the parent carries over;
the re-tuned model point S(0.415) sits between the toy thresholds S(0.40) and S(0.42) the
experiment already reports, and the pooled toy fractions bracket the model from above. The
step count falls by 0.0341 bits to S = 2^125.7960 and the ledger total to 2^125.6734; the
submitted bound 125.68 is the least hundredth above it and leaves a 0.46% allowance. Part D
recomputes the parent ledger and this one through the same code path.

This is an analytic generic-attack bound, not a practical break.
