# Verdict v20n (54681c8b-ff0c-4150-b8c4-4e634b07d039)

- status: rejected
- judge.valid: False

## Reason (verbatim)

Reject the complete original claim under exploratory-v1. The declared success lower bound 0.3739 is below the required 0.39. Neither the raw toy success rate nor the separately labelled 0.3995 model prediction supplies a supported production guarantee meeting that floor. This is insufficient support, not proof that the true production success probability is below 0.39.

There is also a concrete mismatch within the declared claim: the 0.373911 figure comes from walk-census.json's 64-walker arm at visits/sqrt(N)=1.134177, not the claimed 16-walker experiment within 1e-3 of 1.0096963. The actual at-budget 16-walker CP lower bounds in walk-census-exact.json are approximately 0.334008, 0.359947 and 0.341671 before debits. Consequently the claimed statistical support is misidentified, and these small scalar-universe experiments do not establish the stated full-production bound.

The strongest offered response is the birthday-model decomposition with production walker and distinguished-point factors. Its transfer to this predecessor-dependent repaired walk remains unsupported, particularly for cycle effects and success lost to the expenditure cap. The missing matched control cannot resolve this: static inspection shows that both census_v20n.py's walk_arm and match_control.py's run never assign the returned cn back to canon inside their loops. They update coefficients while retaining the initial class, rather than advancing the stated walk. The reported zero-escape failure therefore does not justify the note's diagnosis that the starting-state block is responsible.

The supplied absorption measurements and restart repricing address local issues but do not close these gaps. In particular, the large absorption fixture has p_ge3=1; approximately 0.000732 is only the exact-three-cycle share, not a bound on all residual cycles or production success loss. Two toy sizes and a citation audit cannot establish that extrapolation.

The ledger reports log2(W)=126.26594404, rounded to the submitted 126.266 in the required units. Arithmetic aggregation alone does not establish recovery within that budget: the submission expressly supplies no cap-loss bound, and its success statistic is not for the complete capped production attack. Thus the resource and success claims cannot jointly qualify. Scores remain null as required by the protected binding configuration.
