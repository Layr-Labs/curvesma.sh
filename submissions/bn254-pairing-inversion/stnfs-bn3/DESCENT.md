# Finite-family sieve descent and independently randomized boot

hybrid_model.py is the executable cost model. It keeps the original82 node windows and radii in descent16.PLAN16 and the sieve amortization. It supersedes the earlier negative-binomial counts, lazy cofactor discount, DISC, BALL_PAIR, subtree independence, node-counter and birthday-repeat calculations. The only run abort comes from the complete-work meters in proof.md.

## Fixed auxiliary families

For each side, enumerate the unramified degree-one root ideals (p,y-s,x-t) with128<=p<256, in lexicographic(p,s,t) order, and retain the first64. These small primes exclude the ramified primes2,5,7 and the w-prime5. The table contains their known factor-base logarithms; its construction and validation are charged in preprocessing. Let L be its actual nonzero size: L_f=12 and L_g=36, both at most64. Every node draws uniformly and independently with replacement from this SAME fixed table. The algorithm checks known-log availability, distinctness and the side-polynomial/base-polynomial congruences. It does not claim64 different rational primes.

For an input ideal q, intersect its lattice with the selected auxiliary q' lattice by CRT. The determinant is qq'. Reduce and certify under QN. Sieve the complete radius-R coefficient ball modulo sign, using the row's original R. All primes in[2^8,2^56] and relative-degree-one roots are sieved; the original2400-instruction norm-bound pass and score threshold apply. Compute surviving norms exactly, remove the known q and q' factors, and use finite cofactorization to verify every candidate relation. The ball's deterministic qualifying relation count includes primitivity, ordinary-ideal conditions, norm/score gates, unit duplicates and missing-log rejection. Factorization randomness is treated separately below.

q' has at most8 determinant bits but at least7 removed bits. Relative to the original single-q norm model, the conservative residual-norm shifts are(+25,+8) for f-special nodes and(+32,+1) for g-special nodes. hybrid_model shifts all norm-bin upper/lower offsets consistently, recomputes retention and surviving-cofactor work, and charges the full ball. No negative-binomial lazy fraction reduces cofactor work.

If a row rejects relations containing unlogged factor-base ideals, that rejection is included in its qualifying count. The other rows recursively solve unlogged small ideals using the original known-log approximately20-bit auxiliary construction; the effective ideal is at most77 bits and is charged as the most expensive discarding helper row at or below that limit. The helper rows themselves reject missing-log relations. This prevents a same-size helper recursion. All other children obey the row's strictly decreasing prime-size windows. The branch sum in hybrid_model includes these missing-log helper expectations.

## Finite-family probability argument

Fix the input q and the previously computed known-log data. For each auxiliary index j, let N_j be the exact number of qualifying SMOOTH relations in that auxiliary's complete ball, before randomized ECM detection. This is a nonnegative integer fixed by the public arithmetic, input ideal and deterministic sieve predicates. Let m=(sum_j N_j)/L. The deterministic identity

E[N_j^2]=(sum_j N_j^2)/L <= (sum_j N_j)^2/L = L m^2

holds for every finite family: no pairwise independence or empirical count-distribution shape is assumed. In particular, if m>0, at least one auxiliary ball has a relation and a uniform auxiliary hits a nonempty ball with probability at least1/L.

For that ball fix its first qualifying relation in the deterministic traversal order. The finite ECM/primality procedure detects and verifies this particular relation with probability at least delta, the separate entire-factorization H-ECM bound delta=1/2 specified below. Thus a fresh ball trial succeeds with probability at least delta/L, conditional on the input ideal and known-log data. Independent auxiliary and schedule coins give expected trials at most L/delta; the evaluator uses the actual side-specific L and fixed delta=1/2. This argument does not replace actual N_j by a smoothness prediction in the algebraic identity.

The remaining H1 premise is therefore explicit: for EVERY input ideal reached by this attack, the fixed finite family has positive actual qualifying relation count. The model's positive predicted family mean, recomputed norm laws and retention estimates are supporting heuristic evidence for that premise; a sum of predicted Dickman probabilities is not an actual relation-count certificate. The finite-family norm campaign measures exact CRT geometry and local count proxies for multiple auxiliaries at the SAME q; its scope and extrapolation to production semismoothness are disclosed in the accompanying evidence. If the family has no qualifying relation, the attack reaches its global meter and fails; H1 does not hide this case. No independent second-moment assumption is asserted.

## Cost and recursive solution

Each ball charges its complete sieve/bound pass and ALL predicted cofactor survivors, with no early-success or negative-binomial discount, plus construction overhead. The expected own work is at most L/delta times that per-ball work. After one relation is found, solve every child recursively; there is no Bernoulli child-success event. By induction over the decreasing windows, helpers as above, and the positive-family premise, all children finish almost surely before application of the global meter. Linearity of expectation gives

W_node <=(L/delta) W_ball + kA(max W_f(nA)+max W_g(nA)) + kB(max W_f(n)+max W_g(n)) + W_missing_helpers.

The model evaluates this finite recurrence in increasing bin order. Repeated ideals merely cause charged redundant work. They are not a probability-of-collision failure, and no caching is needed. The complete-work meter accounts for every retry, failed schedule and recursive child.

## Boot: one cell per fresh exponent

For every boot trial choose an independent uniform nonzero exponent e, form T'=T^e, build and certify its p^6-determinant boot lattice under QB, and draw ONE uniform cell modulo sign from the radius4 coefficient ball. Compute its g-norm and perform finite cofactorization. Accept only the original windows: primes at most260 bits, at most12 between56 and132.5 bits and at most one above132.5. Check primitivity/ordinary conditions and all required ideal data.

A rejected trial starts with a fresh exponent and lattice. This samples the original single-cell boot population independently; it does not traverse whole balls or postulate a within-ball relation-count law. Its probability is the original measured boot-bin Dickman/Murphy model, finite ECM detection, conservative0.65 duplicate factor and counted-cell fraction. The measured boot data comprise6144 exact norms on16 certified lattices at radii4,6,8, with the original shifted decile/sextile envelopes. These support the stated single-trial norm/smoothness hypothesis at their actual scope. They do not measure the production rare-event rate.

The trial price includes construction and exponentiation every time. Expected own work is trial_price/p_trial. Once a relation qualifies, solve all its children as above and undo e. The boot never applies the old eighteen-relation failure probability, BALL_PAIR or subtree-independence assumption. Two logarithms are charged independently in work only; success uses the global expected-work meter, not independence of their outcomes.

## Memory and scope

Unlike the independent-cell fallback, this attack preserves the original production-size sieve balls. The largest radius128 node uses four bytes per cell pair, below2^85.42 bytes; it dominates the existing main-phase storage. Memory is reported separately and is not scored. Norms, windows, auxiliary-family actual positivity, detection, expected work and full-work caps all refer to the full-production BN254 instance. This is a conditional exploratory claim, not a production factorization demonstration.

## Fixed-original-modulus factor extraction

No factorization forest is used. Each finite ECM schedule runs every curve modulo the original composite N throughout. On every curve charge both the ordinary sweep and an unconditional replay. Collect gcd outputs only at stage-1 prime-power prefixes and stage-2 giant-block products. Refine the divisor partition from these outputs; no gcd on every baby/giant difference is assumed. For the exact selected schedule, hybrid_model.fixed_modulus_ecm charges curves*(2*ops*mul + 2*gcd + 5*(pi_upper(B1)+G)*gcd). ops includes the original 10% arithmetic allowance. The five-gcd allowance includes output gcd, remainder formation and partition refinement; refinement reduces against the current component product tree, whose total input bit lengths are bounded by the original modulus. All scratch, control and setup also belong to the per-trial 2^40 allowance. A partial or composite partition is rejected; correctness is checked before accepting a relation.

H-ECM is a separate extraction hypothesis: on any qualifying relation in the specified windows, the complete bounded schedules, prefix/block replay and partition procedure fully extract the required factors with probability at least delta=1/2. This is NOT inferred from the 99.9% single-prime smooth-order hit model. It includes separation failures, setup failures and all component factors on both sides. Node retries are L/delta, with L=12 or36, and the boot single-trial probability uses delta=1/2. 

compound-summary.json records 400 complete, valid planted-product partitions over four predetermined seeds, 64 curves per product, 2-4 factors of10-18bits. These experiments test original-modulus stage1 and unconditional prefix replay only; stage2 extraction and production separation remain extrapolations under H-ECM. Replay is available as {"part":"ecm","seed":1}. count-summary.json records all80 finite-family proxy cases and20,480 production norm pairs, with two base ideals per stratum and four predetermined auxiliary choices. Cutoffs257/4093 gave9/22 paired one-large-prime probable-residual events. Those are permissive proxies, not the production two-window events or a uniform finite-family positivity theorem. The actual-positive-family premise H1 remains explicit. Replay is available as {"part":"family","seed":1}; auxiliary_pools.json gives the complete finite pools.
