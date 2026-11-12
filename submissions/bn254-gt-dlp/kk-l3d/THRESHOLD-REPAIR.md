# Count tails and sequential availability

This repairs the existing traversal without count/descendant independence. H-PGF remains E[t^M]<=G(t)=(1+beta(1-t))^-r for0<=t<=.8, beta=.5mu/r. It applies to the specified pre-entry reached-ideal ensemble, not a completely fixed ideal.

Let F_i be failure of relation i's children and A_i={M>=i,F_1,...,F_(i-1)}. The existing sequential hypothesis is explicitly P(F_i|A_i,prehistory)<=psi whenever P(A_i)>0. Here psi is the sum of its child failure bounds. This conditional hazard and the separate conditional per-attempt child-work envelope remain disclosed transfer heuristics; neither is inferred from count PGFs alone.

A_(i+1) is contained in A_i intersect F_i, so P(A_(i+1))<=psi P(A_i). Induction gives P(M>=k,F_1,...,F_k)<=psi^k. A failed node either has M<k or fails its first k relations. Therefore, for every integer k>=1,

    P(node fails)<=P(M<k)+psi^k.

Since t^M>=t^(k-1) on M<k, Markov gives P(M<k)<=G(t)/t^(k-1). For k1 use t0 and P(M0)<=G0. We use max(.02,min(1,min_(k=1..64)[T_k+psi^k])), where T_1=G0 and, for k>1,

    t_k=min(.8,(k-1)(1+1/beta)/(k-1+r)),
    T_k=min(1,G(t_k)/t_k^(k-1)).

Differentiating log(T) gives r*beta/(1+beta(1-t))-(k-1)/t; its numerator is increasing linear, with zero t_k before clipping. Every finite k/t choice is valid; global optimization is unnecessary. mu0 gives bound1. Replay adds absolute1e-10 outward slack; it is a floating diagnostic of the exact formulas, not an interval proof.

E[attempts]=sum_i P(A_i)<=1/(1-psi). Together with the separate conditional child-work bound this preserves own_full_ball+child_work/(1-psi). All downstream losses and work are repriced using the fixed repaired radii.

The old inference P(fail)<=PGF(psi) was invalid: M1 with probability.1 and a failing sole relation, otherwise M100 with immediate success, passes the shape1 rawmu60 PGF but has failure.1>G(.1)=1/28. The new bound separates the correlated low-count event. The default run checks this counterexample (threshold_counterexample); no stronger joint-M law is introduced.
