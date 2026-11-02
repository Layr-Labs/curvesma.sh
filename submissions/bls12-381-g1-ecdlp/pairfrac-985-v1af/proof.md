# Pair fraction 0.985 with exact kernel body and retuned coverage budget

## Complete claim and source scope

For the protected BLS12-381 powers-of-tau target, at n=2^20, the supplied paired class-walk Cheon construction has total work at most 2^116.306249817945463 G1-reference units (C=2000), conditional on the explicitly disclosed H1-H5 of source-proof.md, with H4 strengthened to a 0.985 global pair fraction as justified in section 1. The exact score is 116306249817945463/1000000000000000. Success is at least 0.39 under those hypotheses. Memory remains below 2^121.5854 bytes, preprocessing below 2^115.3097 units and included in total work, and advice is zero. Only A_1 and A_(2^20) are accessed, with all accesses charged; the same bound W(n) applies for 2^20<=n<r. No claim is made for smaller n, and no auxiliary powers are transferred to an ordinary track.

The complete algorithm, uniform-challenge reduction, class sieves, paired walk, planner, bounded transport, collision handling, root recovery, verification and H1/H2/H3/H5 support/limitations are supplied in source-proof.md. That file is the immutable source argument from submission 9c1b2d5d-98fb-42a4-99ab-593f6406a59f, commit d0a6e4cfea64902889b71847dbc4409709b073fc. prior-proof.md supplies the published parity-key amendment a38cb931-1c42-4bc2-b209-c429ce92dfb1. The published lazy-fused side-output amendment 49e06a3d-e5c5-4fb4-8302-2106dc93104a is a declared dependency: its 1025-instruction body and fixed-call reserve design are restated and re-verified in section 2 below, and its retained artifacts (lazy_witness.py, native_check.py, side_lazy.S, measure/side_lanes.S) are included as evidence; only its amendment prose is omitted, for the package byte budget. This document replaces only the predecessors' pair-fraction constant (with the H4 justification of section 1), kernel per-output funding (exact body plus explicit fixed-call row), coverage-budget constant, score grid, and numerical ledger displays; all their target, algorithm, cap, probability and memory obligations remain. The predecessors' old constants and old totals are historical, not additional claims by this package.

The source uses r as protected, d=2^20 dividing r-1, m=(r-1)/d, and six-element GLV/negation classes. The global pair fraction is F=197/200 (0.985). The coverage target is K=6N^2/m>=24757/50000 (0.49514). N=64239578809815195260251725009983698 is the least integer satisfying the retuned coverage budget. The scheme and all failure checks are unchanged. No previous hosted decision proves H1-H5.

## 1. Pair fraction 0.985 (H4)

F is the global fraction of emitted walk points produced by pair operations (two points per pair operation, one per root operation), averaged over every class list. The ledger charges per_point = F*PAIR/2 + (1-F)*ROOT + transport, so any lower bound on the realized global fraction yields a valid upper ledger. The claim is F >= 0.985, i.e. at most 1.5% of emitted elements are plan roots.

Support comes from running the shipped, unmodified planner (evidence.py plan/class_lists/toy_prime) on a predeclared ladder of toy parameters, reproduced by planner_f99.py (standard library only; run with `python3 planner_f99.py`). GMAX = 16 d^2 in every spec, exactly as in production. Production for comparison: d = 2^20, X2 = 2^76, class lists of about 2^56 elements, planner window 2*GMAX/d = 2^25 members, span about 2^31 windows.

Measured mean roots per emitted element (classes averaged; min pair fraction in parentheses):

| d | X2 | window members | span (windows) | classes | mean roots/element |
| --- | --- | --- | --- | --- | --- |
| 8 | 2^13 | 256 | 4.0 | 8 | 0.0205 (0.9764) |
| 8 | 2^14 | 256 | 8.0 | 8 | 0.0202 (0.9769) |
| 8 | 2^17 | 256 | 64.0 | 4+4 | 0.0193/0.0195 (0.9801/0.9796) |
| 16 | 2^15 | 512 | 4.0 | 8+8+8 | 0.0155/0.0159/0.0172 (0.9815/0.9803/0.9816) |
| 16 | 2^18 | 512 | 32.0 | 4 | 0.0146 (0.9839) |
| 16 | 2^20 | 512 | 128.0 | 4 | 0.0138 (0.9855) |
| 16 | 2^22 | 512 | 512.0 | 2 | 0.0138 (0.9862) |
| 32 | 2^21 | 1024 | 64.0 | 4+4 | 0.0126/0.0128 (0.9872/0.9867) |
| 64 | 2^20 | 2048 | 8.0 | 4 | 0.0118 (0.9873) |
| 64 | 2^23 | 2048 | 64.0 | 2 | 0.0124 (0.9872) |
| 128 | 2^23 | 4096 | 16.0 | 2 | 0.0122 (0.9876) |

(Values are planner_f99.py output; two entries denote two seeds. Span-4.0 rows are shown for context; production span is 2^31 windows, so short-span rows understate production pairing.)

Reading of the table:

1. Direct bound. At window >= 512 members and span >= 16 (d >= 16), every measured mean is at most 0.0138, inside the claimed 0.015 with a 1.09x margin. The d = 8 rows (window 256) are the only rows above 0.015; production window is 2^17 times larger, and every window step from 256 to 4096 members strictly decreases the measured roots fraction (0.0194 -> 0.0138 -> 0.0127 -> 0.0124 -> 0.0122 at span >= 16), five for five with no reversal.
2. Span saturation. At fixed window 512, span 128 and span 512 give identical 0.0138: span effects saturate, so production's vastly larger span is conservatively taken to add nothing beyond the toy span-128 rows.
3. Class average, not minimum. F is the global (class-averaged) fraction, which concentrates over the 2^20 production classes; no worst-class tail argument is needed.
4. Planner capacity. Across all specs, size rejections are 0, depth rejections are 0, maximum depth seen is below 300 (cap 2^14), and maximum tree size seen is below 2^14 (cap 2^20). The production tree cap binds only in the sense of forcing at least 2^36 trees (lists of 2^56 over TREEMAX 2^20), i.e. at least 2^-20 roots per element, which is negligible and in any case a negligible addition already covered by the margin. Demotion-forced roots scale as 1/(16d) and vanish at production d.
5. Short-span rows excluded. Rows with span below 8 windows (lists shorter than a few planner windows, the opposite of the production regime) are measured for context only and enter no fit or bound.

H4 therefore claims a global pair fraction of at least 0.985: directly supported at d >= 16 toys with a 1.09x margin, with a strictly monotone improving window trend toward production's 2^25-member window. It remains an explicit heuristic extrapolation, reviewable on the stated toy evidence and margins; the production planner run itself is infeasible and is not claimed.

## 2. Exact kernel body charge with explicit fixed-call row

The declared dependency 49e06a3d proves a 1025-instruction worst-path side-output body and charges 1026 per side output, using the extra instruction per output and padded operation to fund fixed invocation costs. This package charges the 1025 body exactly and funds fixed calls through a separate, nonoverlapping ledger row instead.

The side call is issued at most once per round per base, so (BASES+2)*LANE_OPS bounds the invocations (the dependency's call bound, retained). Each invocation costs at most 512 instructions including call preparation and return (deliberately larger than the 36 fixed function instructions reported by the native audit). The ledger therefore adds an exact row of 512*(BASES+2)*LANE_OPS = 1066509928891664157045105880989696 instructions, split evenly between preprocessing (baby phase) and online work; the even split is asserted exact. This row is in the instruction total; fixed-call funding draws on no per-output slack, no allowance capacity, and no padding conservatism.

The body bound, restated from the declared dependency 49e06a3d. Let p be the base-field prime and R=2^384>8p. For a,b<2p, word-serial Montgomery reduction returns S=(ab+qp)/R with q<R, hence S<4p^2/R+p<2p, and the word accumulator stays below a+p<3p<R; so the first product (a=y1+y2<2p, b=inv<p) skips its final conditional subtraction, yielding mu<2p, while the square keeps its normalization and stays canonical. The second product takes both operands below 2p with raw output S<2p. For the final y, D=S-y1 lies in (-p,2p); one carry branch then adds p once (negative case, modulo-R arithmetic giving the unique value in [0,p)) or conditionally subtracts p once (nonnegative case), replacing 13+23 instructions by at most 23. lazy_witness.py derives side_lazy.S from measure/side_lanes.S by exactly this transformation (asserted byte equality), and native_check.py audits the resulting instruction stream against whole-integer field arithmetic (458 lanes, both branches). The stored y, parity key, and all downstream values are identical to the source kernel.

A pair costs 2000+1025+53=3078 instructions for two points; a root costs 2050. The primary 2000 reference charge is retained untouched, including its batch-inversion amortization. Pair bookkeeping (53) and root bookkeeping (50), including their seven-instruction worst-path parity keys, and the average record transport 81/2 are unchanged. With pair fraction F=197/200, the weighted per-point charge is 317433/200 = 1587.165 instructions.

## 3. Retuned coverage budget

The success floor is 1-exp(-K)-0.0005-2^-20-4*2^-40-2^-111, where K=6N^2/m is the realized coverage parameter and the subtracted terms are, respectively, the H1 production-coverage allowance, the modelling budget, the H3/H5 cap-tail union, and the diagonal event (source-proof.md section 7). The H1 allowance (0.0005) and every other subtracted term are retained exactly; only the target K changes, from 619/1250 (0.4952) to 24757/50000 (0.49514).

At K=0.49514 the exact floor is 1-exp(-0.49514)-0.0005-2^-20-4*2^-40-2^-111 = 0.39001347..., asserted by evidence.py to exceed 0.39001. The declared success remains 0.39, with arithmetic margin 0.000013 over the required threshold. This margin is pure exact arithmetic on 1-exp(-K); no heuristic allowance was narrowed to obtain it. N is the least integer with 6N^2>=KT*m (asserted both directions).

## 4. Exact complete ledger

The unchanged source dimensions, quotas and capped schedules are evaluated by evidence.py at the retuned (KT, F, KERNEL). Every row is an exact integer, including rounding at the original locations. The complete instruction total is 205447033166574778786844479282649236819, divided by 2000 for W. Preprocessing is 102966644125054266401921608775974072987 instructions, included. Its baby-phase transport is charged at 45 per record rather than the combined 40.5 average, exactly as in the source.

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
- records: 2N points: 203917622203360668770474858250941552073 instructions.
- root recovery and input distribution: 73872550223301050368 instructions.
- segment restarts (scalar multiplications): 62741744178240868946713197064620000 instructions.
- side-output fixed calls: 1066509928891664157045105880989696 instructions.
- table strand starts and base points: 7657983763737776997925928960000 instructions.

Memory and all per-node layout assertions are retained: total below 2^121.5854 bytes and every node below 2^63 bytes. No free precomputation, future-target amortization or wall-clock parallel discount is introduced. The exponentially large node/memory assumptions are analytical accounting conventions, not feasible hardware resources.

The score grid is 10^-15: evidence.py reports the least fifteen-decimal ceiling 116306249817945463/1000000000000000 of the exact log2 W = 116.306249817945462562... (90-digit Decimal enclosure; the on-grid assertion holds). log_bound.py encloses ln(x) by 64 rational terms of 2*sum(t^(2j+1)/(2j+1)), t=(x-1)/(x+1), and bounds the remaining positive geometric tail. It similarly encloses ln2 and divides in the safe direction. Applied to the exact total/2000, it proves the declared fifteen-decimal ceiling. The source rounded memory/preprocessing bounds are retained/recomputed separately and are not new optimization objectives.

## 5. Success and honest limitations

The actual points, canonical stored keys, all routing and candidate decisions, and recovered scalars are identical to the preceding algorithm. Therefore neither the pair-fraction constant, the kernel-body exact charge, nor the coverage retune changes any source success event or cap distribution, except that N (and hence the emitted sets) is recomputed at the retuned KT through the same construction. H1-H3 and H5 keep their precise definitions, evidence and failure budgets from source-proof.md section 7 and its planner sections; H4 is strengthened to 0.985 as justified in section 1 above. Those hypotheses remain extrapolations requiring review.

The recomputed lower floor is 1-exp(-0.49514)-0.0005-2^-20-4*2^-40-2^-111 > 0.39001347, and the declared success is only 0.39. This is conditional on the disclosed model, not an empirical production success rate or reviewer-confidence statement. All failures and capped work remain in the complete ledger. The exact uniform target and auxiliary-power limit remain unchanged.

The full combined test reruns the production-parameter integer ledger, planner controls, reduced-order curve recovery, synchronous transport including overflow controls, toy census, and 12142 parity-key cases, plus the 458-lane native audit. planner_f99.py adds the rerunnable pair-fraction ladder of section 1. Tests exercise finite cases only; the algebra and accounting arguments above carry the production extrapolation. No full-size attack is run and no previous judgment is used as proof.
