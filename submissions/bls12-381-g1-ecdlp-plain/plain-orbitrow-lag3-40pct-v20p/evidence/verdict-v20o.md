# Verdict v20o (b170149d-92a1-449c-bf76-a365ba18ef61)

- status: rejected
- judge.valid: False

## Reason (verbatim)

Reject the complete original claim under exploratory-v1. The required success probability is at least 0.39 over uniform full-production challenges and fresh algorithm coins. The reported 0.2825 toy success rate and 0.2694 lower confidence bound do not meet that requirement. The 0.3995 birthday-model prediction is not a supported production lower bound for this predecessor-dependent repaired walk. The experiments use different column counts, walker counts and distinguished-point densities; no supported transfer or budgeted amplification closes the gap. This is insufficient support, not proof that actual production success is below 0.39.

The corrected canon=cn assignments and shared starting states address the earlier stationary-loop and unmatched-control problems. The reported paired comparison is relevant local evidence for the escape mechanism, but several original measurement claims still have concrete defects.

First, the cap-loss statistic is not the claimed conditional probability. match_control.py divides losses by all 3,200 trials, rather than the 941 would-solve trials. For zero losses among 941, the one-sided 95% binomial upper bound would be approximately 0.00318, not 0.0009357. More fundamentally, meeting times are collected only from the 1,687-step trace: none can exceed its full-length cap by construction. Testing shorter prefixes does not establish success lost beyond that observation window or through the production expenditure meter. The claimed 0.26915 adjustment is therefore unsupported.

Second, census_v20n.py's decomposition computes maximum gaps between flattened trace indices, including arrivals from different walkers. These are not cycle periods of the full predecessor-dependent state. With 16 walkers, a two-step repeat within one lane gives an index gap of 32, which the code classifies as >=4. Thus the claimed cycle-length split is not established by this census.

Third, the lag-2/3 extension calls _apply with the pre-step canon but coefficients already updated by qwalk.step, generally breaking a+c*x=canon. It also retains the distinguished-point hit mask computed before that substitution. Its reported +0.000625 therefore does not validate the claimed invariant-preserving extension and reporting protocol.

The original claim's propagated CP value 0.2784 also conflicts with the supplied calculation, which reports 0.2658668 using an assumed detection factor of 0.9. Neither that factor nor the production extrapolation is justified as a lower bound. The smaller absorption experiments provide evidence at their measured sizes, not a production success-loss bound.

The reported log2 work 126.26594404 rounds upward to the submitted 126.266 in the stipulated units, but ledger arithmetic alone does not establish the complete resource-and-success claim. Its SUCCESS calculation still consumes walk-census.json's larger-budget arm at visits/sqrt(N)=1.134177, rather than the current census near 1.0096963, and its step accounting describes the lag-1 repair rather than the claimed shipped lag-2/3 extension. Consequently the work and success assertions do not jointly qualify. Citation audits do not resolve these mathematical and statistical defects. Scores remain null under the protected binding configuration.
