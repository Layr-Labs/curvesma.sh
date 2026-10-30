This file replaces the predecessor abandonment and selected-merge probability bounds and supplies
the current report-overflow argument. Its M1 is the shared per-class random environment of predecessor
Section 2.0, not per-trail resampling. H1/M1/L-ori remain heuristic premises for real P-256.
The rare meter uses worst-case debits checked BEFORE each operation. Report processing terminates
at the first duplicate without draining a batch of duplicates. These are production procedure invariants.

# All-abandonment envelope and report overflow

**All-abandonment envelope, replacing Lemma M's selected-tail bound.** Put R=TRAIL_MAX, B=DP_MAX+2K, J=67R+1, k=floor(CAP/716800), P_e=(64/(2A))(1+2^-31), and d=CAP-32768k-1 (P=16384: 2P=32768 steps per short-cycle entry). We bound any abandonment on G, including after a merge. Fix target, table and M1's entire per-class random environment. Lemma R gives a common first-DP/short-cycle-escape path for every anchor phase. Analytically expose its first class discoveries, testing all members before escaping a short cycle; stop at the first per-trail G violation. Each short cycle has one fresh-source closing entry. Its repeated live traversals cost at most 2P-1 <= 32767 processed steps by E2; its least-class doubling has a new depth by D1. False look-ahead exposes this same path, not a new path. Therefore abandonment under any phase implies CAP<=1+Q+(2P-1)E_c, where Q counts distinct tested non-start classes and E_c short-cycle entries. A processed step generates at most one ordinary, 65 look-ahead and one escape point, so J bounds all discovery opportunities, including speculative ones.

E1/M1/L-ori give conditional entry probability at most P_e at a source's first exposure. The small pad covers key residue rounding. Zero-pad after termination or the G stop, without conditioning on G. Factorial moments give Pr(E_c>=k)<=binom(J,k)P_e^k<=(eJP_e/k)^k. If E_c<k, at least d distinct DP tests on the common discovery prefix failed. Their M1 first-exposure probability is at most (1-theta)^d<=exp(-theta*d). Previously exposed cycle/look-ahead classes supply no additional tests. A union bound gives q_*<=(eJP_e/k)^k+exp(-theta*d). This bounds the common envelope for all phases, not a separate tail conditioned on a selected phase. Here eJP_e/k<0.98, B*0.98^k<2^-100 and B*exp(-theta*d)<1.2413e-6.

For environment E let q_E be the fresh-uniform-start probability of this envelope before the local G stop and before any conditioned accepted-adder hit, which is paid separately. The single-trail M1 calculation gives E_E q_E<=q_*. Given E and any earlier history, each next start still has independent fresh uniform coins. Its adaptive phase is covered by the common envelope. Hence for the ith of at most B starts, Pr(i starts and abandons on G)<=E q_E<=q_*, and Pr(any abandonment,G,no accepted-adder hit)<=Bq_*<1.3e-6. We averaged over the shared environment after using independent start coins; we did not resample a known prefix or assert a conditional geometric law. Conditioned accepted-adder hits are stopped and paid separately. The M1 per-class random-environment premise remains heuristic. No adaptive-set H1 hit law is assumed.

On G, without abandonment or another charged abort, P/R ensure the first merge reports a common DP and resolves. Thus DELTA_ANY_ABANDON=13/10^7 replaces both DELTA_ABANDON and DELTA_MERGE. It also covers excess-abandonment failure. The original selected-history tail statement is deleted.

## Finite key-index coupling

M1 assigns each class a key uniform on [0,p), not an exactly uniform residue modulo A.
Write p=q*A+r0 with 0<=r0<A. Generate a block H: each of 0..q-1 has probability A/p,
and the last partial block H=q has probability r0/p. In a coupled ideal experiment the low
index I is independent uniform on 0..A-1 in EVERY block. In the real experiment it is the same I
on a full block and a fresh uniform residue 0..r0-1 in the last block. Couple orientations identically.
The experiments agree unless a partial block is exposed, probability r0/p<A/p per fresh class.
The ideal DP indicator is 1{H<2^196/A}, independent of I and of orientation; its mean is exactly
2^196/p, since the threshold is a full-block boundary. Conditioning on DP=false therefore leaves
I exactly uniform in the ideal model. This is a coupling, not exact real-key index uniformity.

The actual report analysis exposes at most n classes as below. Each of at most B started trails'
common-phase envelope exposes at most J additional candidates; include these hypothetical discovery
paths too. Stopped first exposures and a union bound give coupling failure at most
(n+B*J)*A/p < 2^-80, checked exactly in accounting.py. Repeated keys add no new coupling trial.
All exact-uniform index/window calculations in the retained E1/F and two-cycle analysis refer to
this ideal experiment, stopped at the coupling failure, not to the raw p-element key model.
The 2^-80 charge is included in the 2^-59 small bucket. No history-dependent bad-set hit law is added.
The per-class environment and fresh-start averaging use the same coupled environment throughout.

## Report overflow

**Report overflow.** Stop at the first G violation, duplicate resolution, resource/meter abort or attempted count overflow. Use the original space with zero padding. The rare meter uses worst-case primitive debits and checks before operations. Every auxiliary generated point, tested or not, consumes at least ten credits within those debits, so its count is at most floor(B_rare/10). This includes early-exit look-ahead and escape, not floor(B_rare/W_ENTRY). At most B=D+2K trails start before the stops, D=DP_MAX. Including all four point attempts per capped adder trial/start and K+1 extra handoff opportunities, put n=S+floor(B_rare/10)+4(A*2^40+B)+K+1. Short trails and identity retries are included. Replay begins after the duplicate stop and contributes no report trials.

At each first class-key exposure set Y=1{key<2^196}; M1 gives conditional mean theta. Repeated classes reuse Y. Before a duplicate, report keys are distinct, so reports<=sum Y+1, including a terminal duplicate even if overflow is tested first. A positive class exposed at a start or speculative point is charged to that earlier first exposure. Serial first-match termination prevents a batch drain of repeated positive reports. For lambda>0, conditional expectation gives E exp(lambda sum Y)<=exp(mu(exp(lambda)-1)), mu=n*theta. Here mu<1.01*S*theta and (D-1)/mu>1.9829. Thus Pr(overflow,G,before other aborts)<=Pr(sum Y>=D-1)<=exp(-mu/10)<=exp(-S*theta/10)<2^-60, using eta=1/2. No Chernoff assertion is made for repeated walk steps.

Allocate 2^-60 for this overflow, 1/r for the degenerate first match, and 2^-95 for exceptional additions. The latter safely bounds n*A/N under the same fixed-adder-set H1 accounting already used for DELTA_ADDER_CLASS; a zero denominator/unhandled ordinary addition requires an adder-class hit, while nonidentity doubling cannot give identity at odd prime order. The hit event is after setup; including construction in n only enlarges its numeric opportunity ceiling. Sampling/identity retry exhaustion stays in DELTA_CAP. Together with the 2^-80 finite-index coupling loss, their sum is below 2^-59, replacing the old unnamed 2^-20 bucket. DELTA_CUT, CAP, PROBE, HASH, CORR and the rare-meter losses retain their separate charges. Keeping DELTA_ADDER_CLASS separately is conservative double payment. No abandonment loss remains in this bucket.

## Current success allocation

The conditional success lower bound is

    1-exp(-tau_exact*(1-epsilon))
      -1/KM-2^-100-2^-59-DELTA_CORR-DELTA_ADDER_CLASS
      -DELTA_ANY_ABANDON-DELTA_CUT-DELTA_CAP-DELTA_PROBE-DELTA_HASH.

With DELTA_ANY_ABANDON=13/10^7 it is greater than 0.39000025; at A=2^32 and P=16384 the numerical evaluation is approximately 0.390000279549271. The ledger checks a 2.5e-7 margin. The work and memory formulas and DP_MAX are recomputed at A=2^32; the declared score is 63364359/500000 (126.728718).
