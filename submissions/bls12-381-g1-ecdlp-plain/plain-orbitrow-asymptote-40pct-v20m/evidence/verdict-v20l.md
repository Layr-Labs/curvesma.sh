# Verdict v20l (2d98d2a8)

- id: 2d98d2a8-e5c8-4253-9ecf-200e1baf4c18
- status: rejected
- judge.valid: False

## Reason (verbatim)

Reject the entire declared claim under exploratory-v1. The submission does not support the required production success probability at the stated work, and several specific claims conflict with the supplied implementation.

1. The declared supported success bound, approximately 0.373911, is below the mandatory 0.39. The current ledger and public note acknowledge this. A fixture Clopper–Pearson lower bound below 0.39 does not prove that actual success is below 0.39, but neither does it establish the required guarantee. Removing the oversized drain window is a useful correction, not a substitute for that guarantee.

2. The census is not at the ledger's normalized work budget. From ledger.py, S_ops/sqrt(N) is approximately 1.015*sqrt(2*(-ln(0.6097))*N_CELLS/(N_CELLS-1)) = 1.0097, whereas the cited large fixture uses 1.1342 visits/sqrt(N). It therefore supplies about 12.3% more normalized walking work. Moreover, two increasing fixture estimates do not justify a production monotonicity assumption or make the larger fixture a conservative production bound. The small scalar-group experiments and their statistical intervals remain evidence only for their tested regimes.

3. The absorption experiment does not implement the claimed predecessor-dependent process. In absorb_direct.arm, the returned predecessor pv is ignored: the assignment canon, prev = cn, prev_in leaves prev fixed at the initial class. The escape consequently tests returns to that initial class, not to the preceding class. The call also sets kcap=None, disabling the claimed restart threshold, and the repetition key omits the escape-run counter. Thus the reported counts do not establish absorption probabilities for the shipped process as claimed.

4. Restart statistics and pricing remain unsupported at their asserted scope. walk_census.run initializes restarts=0 but never updates it, so its reported zero is not an independently accumulated restart count. Even genuine zero observations in these short fixtures would not transfer the rule-of-three trial bound to production-length attacks without supporting reasoning. Separately, restart_count.py and restart-count.json explicitly price each of the 255 doublings and 64 window additions at one reference operation, contrary to the claim that doublings were independently instruction-priced rather than treated as affine additions. Executing the known-scalar toy fresh_state does not validate production scalar-multiplication instruction costs. A deterministic work cap alone also supplies no bound on success lost when that cap halts execution.

These are concrete support and implementation failures, not a demand for an exponentially large production run. The submitted computation bound and complete success/absorption/resource claim cannot be accepted exactly as written.
