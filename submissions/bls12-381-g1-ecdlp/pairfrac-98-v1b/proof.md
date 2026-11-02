# Exact kernel body and retuned coverage budget for the powers-of-tau target

## Complete claim and source scope

For the protected BLS12-381 powers-of-tau target, at n=2^20, the supplied paired class-walk Cheon construction has total work at most 2^116.308545645980436 G1-reference units (C=2000), conditional on the explicitly disclosed H1-H5 of source-proof.md. The exact score is 116308545645980436/1000000000000000. Success is at least 0.39 under those hypotheses. Memory remains below 2^121.5854 bytes, preprocessing below 2^115.312 units and included in total work, and advice is zero. Only A_1 and A_(2^20) are accessed, with all accesses charged; the same bound W(n) applies for 2^20<=n<r. No claim is made for smaller n, and no auxiliary powers are transferred to an ordinary track.

The complete algorithm, uniform-challenge reduction, class sieves, paired walk, planner, bounded transport, collision handling, root recovery, verification and H1-H5 support/limitations are supplied in source-proof.md. That file is the immutable source argument from submission 9c1b2d5d-98fb-42a4-99ab-593f6406a59f, commit d0a6e4cfea64902889b71847dbc4409709b073fc. prior-proof.md supplies the published parity-key amendment a38cb931-1c42-4bc2-b209-c429ce92dfb1. prior2-proof.md supplies the published lazy-fused side-output amendment 49e06a3d-e5c5-4fb4-8302-2106dc93104a. All three are included as evidence, not prior verdicts. This document replaces only their kernel per-output reserve, coverage-budget constant, score grid, and numerical ledger displays; all their target, algorithm, cap, probability and memory obligations remain. The predecessors' old kernel prices, old coverage constant, old totals and old declarations are historical, not additional claims by this package.

The source uses r as protected, d=2^20 dividing r-1, m=(r-1)/d, six-element GLV/negation classes, and minimum pair fraction 49/50. The coverage target is retuned to K=6N^2/m>=24757/50000 (0.49514). N=64239578809815195260251725009983698 is the least integer satisfying the retuned coverage budget. The scheme and all failure checks are unchanged. No previous hosted decision proves H1-H5.

## 1. Exact kernel body charge

prior2-proof.md proves a 1025-instruction worst-path side-output body and then charges 1026 per side output, allocating one instruction beyond the body. That extra instruction is pure slack: fixed call costs are covered separately and explicitly. The ledger reserves 512 instructions per possible invocation, and asserts the reserve dominates the bound:

- fixed-call bound: 512*(BASES+2)*LANE_OPS = 1066509928891664157045105880989696,
- fixed-call reserve: F*N+PADOPS = 1574286273568044487852289600231039001/25, about 6.30e34 (about 59 times the bound).

The assertion `fixed_call_reserve >= fixed_call_bound` is evaluated by evidence.py and holds. The 1025 body is itself already a worst-path upper bound (prior2-proof.md section 1; the native audit in native_check.py verifies the exact 26-per-lane saving on both final-subtraction branches, and no branch-frequency assumption is used). Charging the exact worst-path body therefore removes no coverage: every executed instruction remains charged, either in the body or in the separate fixed-call reserve.

A pair costs 2000+1025+53=3078 instructions for two points; a root costs 2050. The primary 2000 reference charge is retained untouched, including its batch-inversion amortization, which a side-output formula alone does not replace. Pair bookkeeping (53) and root bookkeeping (50), including their seven-instruction worst-path parity keys, and the average record transport 81/2 are unchanged. With the unchanged pair fraction F=49/50, the weighted per-point charge is 39743/25 = 1589.72 instructions.

## 2. Retuned coverage budget

The success floor is 1-exp(-K)-0.0005-2^-20-4*2^-40-2^-111, where K=6N^2/m is the realized coverage parameter and the subtracted terms are, respectively, the H1 production-coverage allowance, the modelling budget, the H3/H5 cap-tail union, and the diagonal event (source-proof.md section 7). The H1 allowance (0.0005) and every other subtracted term are retained exactly; only the target K changes, from 619/1250 (0.4952) to 24757/50000 (0.49514).

At K=0.49514 the exact floor is 1-exp(-0.49514)-0.0005-2^-20-4*2^-40-2^-111 = 0.39001347..., asserted by evidence.py to exceed 0.39001. The declared success remains 0.39, with arithmetic margin 0.000013 over the required threshold. This margin is pure exact arithmetic on 1-exp(-K); no heuristic allowance was narrowed to obtain it. N is the least integer with 6N^2>=KT*m (asserted both directions), so the realized K is the minimal feasible value and no coverage is wasted.

## 3. Exact complete ledger

The unchanged source dimensions, quotas and capped schedules are evaluated by evidence.py at the retuned (KT, F, KERNEL). Every row is an exact integer, including rounding at the original locations. The complete instruction total is 205774230904364042770467320491569263820, divided by 2000 for W. Preprocessing is 103130242993948898393733029380434086487 instructions, included. Its baby-phase transport is charged at 45 per record rather than the combined 40.5 average, exactly as in the source.

- base interval sieve, characters, assignment: 14778874133111650197864576 instructions.
- candidate delivery and replies: 5770237022568448 instructions.
- candidate re-walks: 1257895200023828017403475787776 instructions.
- candidate tests: 13858116485374300651520 instructions.
- character sieve, class lists, their transfer: 30225869176339564179694288896 instructions.
- comparisons at the per-owner cap: 175596271737596302496818302434672640 instructions.
- gap tables: 127161744382312499152959732592009770 instructions.
- global barriers (wakes, sweeps, releases): 6630756600764919022560032736 instructions.
- listing allowance: 8 instructions per point: 1027833260957043124164027600159739168 instructions.
- owner initialisation: 83076749736557315843464236105728000 instructions.
- padding ops at the pair price: 51290896618689753540182191330295808 instructions.
- per-wave link words, memory barriers: 15209023373991571020430912520192 instructions.
- plan distribution (copy tree): 121734115224496661965030475956224 instructions.
- planner (scans, matching, emission): 487728773955200596712816258318336 instructions.
- planner descendant traversals (adjacency lists): 10161011842454401296372011630592 instructions.
- records: 2N points: 204245886451078824418254744565742568770 instructions.
- root recovery and input distribution: 73872550223301050368 instructions.
- segment restarts (scalar multiplications): 62741744178240868946713197064620000 instructions.
- table strand starts and base points: 7657983763737776997925928960000 instructions.

Memory and all per-node layout assertions are retained: total below 2^121.5854 bytes and every node below 2^63 bytes. No free precomputation, future-target amortization or wall-clock parallel discount is introduced. The exponentially large node/memory assumptions are analytical accounting conventions, not feasible hardware resources.

The score grid is 10^-15: evidence.py reports the least fifteen-decimal ceiling 116308545645980436/1000000000000000 of the exact log2 W = 116.308545645980435809... (90-digit Decimal enclosure; the on-grid assertion holds). log_bound.py encloses ln(x) by 64 rational terms of 2*sum(t^(2j+1)/(2j+1)), t=(x-1)/(x+1), and bounds the remaining positive geometric tail. It similarly encloses ln2 and divides in the safe direction. Applied to the exact total/2000, it proves the declared fifteen-decimal ceiling. The source rounded memory/preprocessing bounds are retained/recomputed separately and are not new optimization objectives.

## 4. Success and honest limitations

The actual points, canonical stored keys, all routing and candidate decisions, and recovered scalars are identical to the preceding algorithm. Therefore neither the kernel-body exact charge nor the coverage retune changes any source success event or cap distribution, except that N (and hence the emitted sets) is recomputed at the retuned KT through the same construction. H1-H5 keep their precise definitions, evidence and failure budgets from source-proof.md section 7 and its planner sections. Those hypotheses remain extrapolations requiring review.

The recomputed lower floor is 1-exp(-0.49514)-0.0005-2^-20-4*2^-40-2^-111 > 0.39001347, and the declared success is only 0.39. This is conditional on the disclosed model, not an empirical production success rate or reviewer-confidence statement. All failures and capped work remain in the complete ledger. The exact uniform target and auxiliary-power limit remain unchanged.

The full combined test reruns the production-parameter integer ledger, planner controls, reduced-order curve recovery, synchronous transport including overflow controls, toy census, and 12142 parity-key cases, plus the 458-lane native audit. Tests exercise finite cases only; the algebra and accounting arguments above carry the production extrapolation. No full-size attack is run and no previous judgment is used as proof.
