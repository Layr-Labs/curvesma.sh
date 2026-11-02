# Exact kernel body and retuned coverage budget for the powers-of-tau target

## Complete claim and source scope

For the protected BLS12-381 powers-of-tau target, at n=2^20, the supplied paired class-walk Cheon construction has total work at most 2^116.303717121705438 G1-reference units (C=2000), conditional on the explicitly disclosed H1-H6 of source-proof.md and section 6 below. The exact score is 116303717121705438/1000000000000000. Success is at least 0.39 under those hypotheses. Memory remains below 2^121.5854 bytes, preprocessing below 2^115.3125 units and included in total work, and advice is zero. Only A_1 and A_(2^20) are accessed, with all accesses charged; the same bound W(n) applies for 2^20<=n<r. No claim is made for smaller n, and no auxiliary powers are transferred to an ordinary track.

The complete algorithm, uniform-challenge reduction, class sieves, paired walk, planner, bounded transport, collision handling, root recovery, verification and H1-H5 support are in source-proof.md, the immutable source argument from submission 9c1b2d5d-98fb-42a4-99ab-593f6406a59f, commit d0a6e4cfea64902889b71847dbc4409709b073fc. prior-proof.md supplies the published parity-key amendment a38cb931-1c42-4bc2-b209-c429ce92dfb1. Declared dependency 49e06a3d-e5c5-4fb4-8302-2106dc93104a supplies the 1025-instruction body and fixed-call design restated in section 1, with retained artifacts (lazy_witness.py, native_check.py, side_lazy.S, measure/side_lanes.S); only its amendment prose is omitted, for the byte budget. The direct base c5eeef54-4766-47d0-a0cd-4429cd61b2a1 supplies the explicit rows, aggregate residual bounds, 33-bit repair, retuned budget and grid, retained verbatim. The rejected 0.985 attempt 8d7d6d7c is cited for its H4-mean ladder design only (remeasured independently by h6_ladder.py); its H6 reasoning is replaced by section 6. Ancestor-only Callgrind build/run shell recipes are omitted for the budget and retained in ancestor packages; all machine-audited kernel artifacts and the measurement log are included. This document changes only the allowance level (5 to 3, section 3), the priced pair fraction (49/50 to 197/200 pooled, section 6) and numerical displays. Predecessors' old constants and totals are historical, not additional claims.

The source uses r as protected, d=2^20 dividing r-1, m=(r-1)/d, six-element GLV/negation classes, and minimum pair fraction 49/50. The coverage target is K=6N^2/m>=24757/50000 (0.49514). N=64239578809815195260251725009983698 is the least integer satisfying the retuned budget. The scheme and all failure checks are unchanged. No previous hosted decision proves H1-H6.

## 1. Exact kernel body charge with explicit fixed-call row

Declared dependency 49e06a3d proves a 1025-instruction worst-path side-output body. This package charges the 1025 body exactly and funds fixed calls through a separate, nonoverlapping ledger row.

The side call is issued at most once per round per base, so (BASES+2)*LANE_OPS bounds the invocations (the dependency's call bound, retained). Each invocation costs at most 512 instructions including call preparation and return (deliberately larger than the 36 fixed function instructions of the native audit). The ledger adds an exact row of 512*(BASES+2)*LANE_OPS = 1066509928891664157045105880989696 instructions, split evenly between preprocessing and online work; the even split is asserted exact. This row is in the instruction total; fixed-call funding draws on no per-output slack, no allowance capacity, and no padding conservatism.

The body bound, restated from 49e06a3d. Let p be the base-field prime and R=2^384>8p. For a,b<2p, word-serial Montgomery reduction returns S=(ab+qp)/R with q<R, hence S<4p^2/R+p<2p, and the word accumulator stays below a+p<3p<R; so the first product (a=y1+y2<2p, b=inv<p) skips its final conditional subtraction, yielding mu<2p, while the square keeps its normalization and stays canonical. The second product takes both operands below 2p with raw output S<2p. For the final y, D=S-y1 lies in (-p,2p); one carry branch adds p once (negative case, modulo-R arithmetic giving the unique value in [0,p)) or conditionally subtracts p once (nonnegative case), replacing 13+23 instructions by at most 23. lazy_witness.py derives side_lazy.S from measure/side_lanes.S by exactly this transformation (asserted byte equality), and native_check.py audits the instruction stream against whole-integer field arithmetic (458 lanes, both branches). The stored y, parity key, and all downstream values are identical to the source kernel.

A pair costs 2000+1025+53=3078 instructions for two points; a root costs 2050. The primary 2000 reference charge is retained untouched, including its batch-inversion amortization. Pair bookkeeping (53) and root bookkeeping (50), including seven-instruction worst-path parity keys, and the average record transport 81/2 are unchanged. With the pooled pair fraction F=197/200 of section 6, the weighted per-point charge is 317433/200 = 1587.165 instructions.

## 2. Retuned coverage budget

The success floor is 1-exp(-K)-0.0005-2^-20-4*2^-40-2^-111, where K=6N^2/m is realized coverage and the subtracted terms are the H1 allowance, the modelling budget, the H3/H5 cap-tail union, and the diagonal event (source-proof.md section 7). The H1 allowance (0.0005) and every other term are retained exactly; only the target K changes, from 619/1250 (0.4952) to 24757/50000 (0.49514).

At K=0.49514 the exact floor is 1-exp(-0.49514)-0.0005-2^-20-4*2^-40-2^-111 = 0.39001347..., asserted by evidence.py to exceed 0.39001. The declared success remains 0.39, with arithmetic margin 0.000013. This margin is pure exact arithmetic on 1-exp(-K); no heuristic allowance was narrowed. N is the least integer with 6N^2>=KT*m (asserted both directions), so realized K is minimal feasible and no coverage is wasted.

## 3. Explicit pool-cap charge and allowance 3

The source charges an 8-instructions-per-point allowance row for unlisted bookkeeping, naming its largest known residual: the owner pool-cap compare, two instructions per baby record (source-proof.md section 5.2). The base charges it explicitly: every baby record is inserted once into its owner's pool, and the insert path performs one pool-counter bound compare plus its not-taken branch (subs/cmp plus b.hi, two instructions), for an exact row of N*2 = 128479157619630390520503450019967396 instructions. Giant records only read their head slot and walk the chain (REC_LOOK); they never insert and need no pool-cap check. The pool-cap row sits entirely in the baby phase and is included in preprocessing in full.

With the pool-cap residual and the head-index combine explicit, the allowance row drops from 5 to 3 instructions per point (2N*3 = 385437472858891171561510350059902188), removing 4N instructions (2 per point) while keeping the dominant residuals auditable.

Aggregate residual bounds, retained from the base (verdicts cited as procedural context only). (i) Pool-cap compares: charged exactly (N x 2 above), residual zero. (ii) Primary-key reloads: zero. The side kernel stores the key to the record slot and the generator reloads it with the listed ldp+ldr pair (side records; loop control for main records); the key then travels inside the 40-byte record through already-charged loads, and no listing reloads a key anywhere else. (iii) Head-index word combine: charged exactly. The specified 57-bit index (key bits 60..116) straddles two 64-bit words on the stated record layout, so one extraction cannot produce it; the least sequence is lsr plus bfi (two instructions where the listing charges one). Every baby insert and giant lookup performs exactly one extraction: 2N = 128479157619630390520503450019967396 instructions, split exactly N baby (preprocessing) and N giant (online). The index is consumed by the head load on the next instruction, so no re-extraction path exists and the residual beyond the new row is zero. (iv) Comparison-cap enforcement: charged exactly by one down-counter check per comparison (cmp_cap x OWNERS = 25085181676799471785259757490667520, online only, no preprocessing term); the listed per-comparison loop control is retained, so any overlap double-charges conservatively. (v) Plan, batch and output pass control: below 2^115.2 instructions, over 8x below the retained margin. Rows with listed loop control (walk, transport, comparison, restart rows) need nothing further. Owner-init loops dominate the remainder: per owner 2^56 charged stps iterate with subs+b.ne control, unrolled fourfold at identical stp count (setup inside the listed +64 per owner), giving below 2^115.01 including the outer owner loop; even without unrolling, below 2^117.1. All other bulk rows (plan distribution, planner scans and traversals, sieves, candidate delivery/tests, barriers, table strands, root recovery) charge below 2^109 instructions in sum; every loop iteration at every nesting level performs at least one charged operation with at most subs+b.ne control per level, and at most four levels, so their aggregate control is below 8x the charged sum (below 2^112). Total (v) is below 2^115.2 against the retained 3 x 2N = 385437472858891171561510350059902188 (2^118.21). With (i)-(v) bounded or charged, the retained allowance covers only purely-unknown misses.

## 4. Exact complete ledger

The unchanged source dimensions, quotas and capped schedules are evaluated by evidence.py at (KT, F, KERNEL). Every row is an exact integer, including rounding at the original locations. The complete instruction total is 205086680875392687087068228690080002151, divided by 2000 for W. Preprocessing is 103166429215152791659180991640755122287 instructions, included. Its baby-phase transport is charged at 45 per record rather than the combined 40.5 average, and its baby-phase walk part at the conservative floor F_PRE=39/40 of section 6, not at the pooled F.

- base interval sieve, characters, assignment: 14778874133111650197864576 instructions.
- candidate delivery and replies: 5770237022568448 instructions.
- candidate re-walks: 1257895200023828017403475787776 instructions.
- candidate tests: 13858116485374300651520 instructions.
- character sieve, class lists, their transfer: 30225869176339564179694288896 instructions.
- comparison cap enforcement: 25085181676799471785259757490667520 instructions.
- comparisons at the per-owner cap: 175596271737596302496818302434672640 instructions.
- gap tables: 127161744382312499152959732592009770 instructions.
- global barriers (wakes, sweeps, releases): 6630756600764919022560032736 instructions.
- head-index word combine: 128479157619630390520503450019967396 instructions.
- listing allowance: 3 instructions per point: 385437472858891171561510350059902188 instructions.
- owner initialisation: 83076749736557315843464236105728000 instructions.
- owner pool-cap compares: 128479157619630390520503450019967396 instructions.
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

The score grid is 10^-15: evidence.py reports the least fifteen-decimal ceiling 116303717121705438/1000000000000000 of the exact log2 W = 116.303717121705437899... (90-digit Decimal enclosure; the on-grid assertion holds). log_bound.py encloses ln(x) by 64 rational terms of 2*sum(t^(2j+1)/(2j+1)), t=(x-1)/(x+1), and bounds the remaining positive geometric tail. It similarly encloses ln2 and divides in the safe direction. Applied to the exact total/2000, it proves the declared fifteen-decimal ceiling.

## 5. Success and honest limitations

The actual points, canonical stored keys, all routing and candidate decisions, and recovered scalars are identical to the preceding algorithm. Neither the kernel-body exact charge, the coverage retune, nor the pooled-fraction repricing changes any source success event or cap distribution, except that N is recomputed at the retuned KT through the same construction. H1-H5 keep their definitions, evidence and failure budgets from source-proof.md section 7 and its planner sections; H6 is defined in section 6. The allowance change is exact accounting with a stated residual enumeration. Those hypotheses remain extrapolations requiring review.

The recomputed floor is 1-exp(-0.49514)-0.0005-2^-20-4*2^-40-2^-111 > 0.39001347, and the declared success is only 0.39. This is conditional on the disclosed model, not an empirical production success rate or reviewer-confidence statement. All failures and capped work remain in the complete ledger.

The full combined test reruns the production-parameter integer ledger, planner controls, reduced-order curve recovery, synchronous transport including overflow controls, toy census, and 12142 parity-key cases, plus the 458-lane native audit. Tests exercise finite cases only; the algebra and accounting above carry the production extrapolation. No full-size attack is run and no previous judgment is used as proof.

## 6. Pooled pair fraction 197/200: H4 class mean plus jointly measured H6 weighting

The scored records row prices all 2N emitted points at the pooled pair fraction F=197/200. Two premises. First, the production class-mean pair fraction is at least 0.9862 (H4 mean). Second, the production selected-base-weighted pooled fraction is at least 0.985 (new disclosed H6). The per-class floor 49/50 of source H4 is retained alongside both, and bounds the worst case.

H4 mean. The shipped planner was run at production degree parameter GMAX=16d^2 on real chi-class lists (h6_ladder.py mean-only specs, rerunnable; the declared hosted integrated run recomputes two small joint specs (h6_hosted output) and the table below is quoted from the same code offline, as with the source planner tables):

| d | span X2/2GMAX | classes | mean f | min f |
| --- | --- | --- | --- | --- |
| 16 | 16.0 | 4 | 0.98651 | 0.98601 |
| 32 | 8.0 | 4 | 0.98689 | 0.98589 |
| 64 | 4.0 | 2 | 0.98646 | 0.98502 |

Every row sits above 0.9862; the source tables show the same level (0.9861 at d=16 span 16, 0.9874 at d=64 span 16) with means rising in span and in d (uncapped controls 0.9794, 0.9834, 0.9858, 0.9859, 0.9879 at d=8..128). Production (span 2^31, d=2^20) lies far beyond every measured row along both rising trends, and every deficit source (boundary GMAX/span, cycle 1/d, demotions) vanishes at the production ratios. H4-mean is the disclosed extrapolation that the measured level and trends continue; the production run is infeasible and is not claimed.

H6 (joint weighting measurement). Production emitted-point weights depend jointly on class lengths, counts of selected base primes of the opposite character, and the final partial walks. h6_ladder.py measures exactly this joint process at toy scale: real chi-class lists, the shipped planner on every class, bases as consecutive primes above X2 in order, baby taking the first N_t elements and giant the next (last base partial), each base walking the plan of its opposite-character class. Twenty joint specs (d=8..64, 7..21 walked classes, 20..54 bases, uncapped plus production-ratio-capped max_tree=d controls) decompose pooled-minus-class-mean into selection, length-weighting, base-count-weighting and partial-walk steps. The full 20-row table is reproduced by `python3 h6_ladder.py < offline-specs.json` (about one minute, standard library); summary:

- pooled-minus-walked-mean over all 20 rows: range [-0.00185, +0.00282], mean +0.00014;
- pooled-minus-class-mean over all 20 rows: range [-0.00214, +0.00233], mean +0.00030;
- in the 13 rows with at least 12 walked classes, |pooled-minus-walked-mean| <= 0.00044;
- the 3 rows above 0.0005 all walk at most 8 classes with up to 31% of emitted weight on a single class;
- length-weighting step over all 20 rows: range [-0.00001, +0.00032], nonnegative in 18 of 20 rows (the other two are -0.00001, rounding-scale);
- base-count step: [-0.00133, +0.00191]; partial-walk step: [-0.00082, +0.00060];
- capped controls (cap binding, fractions 0.82-0.94): weighting deltas -0.00043, +0.00023, +0.00282, same scale.

Production weights are at least 10^4x more diffuse than the most concentrated measured row: one class contributes at most 2^42 bases x 2^57 elements = 2^99 of N ~= 2^115.63 emitted points per side, i.e. max single-class weight 2^-16.6, against 1/7 in the worst measured row; all 2^20 classes are walked (about 2^39 bases per character per side), so the selection step vanishes. The measured deltas center on zero with spread shrinking in diffuseness and no systematic negative component across 20 diverse fixed weightings. The H6 margin, 0.9862-0.985=0.0012, covers 2.7x the worst drag measured at production-like diffuseness (0.00044). The partial-walk contribution is additionally bounded analytically: at most one remainder per side, under 2^57 elements each, with pair-fraction deficit under 0.3 against N ~= 2^115.63 points, i.e. under 2^-60 pooled drag. Adversarial alignment of weights with low-fraction classes would break H6, but is measured-absent and mechanism-opposed (H2 equidistribution over 2^20 diffuse classes); this residual risk is disclosed, not hidden.

Phase-safe accounting. The combined records row is exactly linear in F (2N x (F x PAIR/2 + (1-F) x ROOT + 81/2)), so the pooled H6 bound covers any baby/giant split with no per-phase fraction asserted for the scored total. Preprocessing is funded separately and conservatively: its baby-phase walk part uses F_PRE=39/40=0.975, below the retained H4 per-class floor, with no weighting claim at all (103166429215152791659180991640755122287 instructions, below 2^115.3125 units).

## 7. Segment-index encoding repair (33 bits)

The source specifies a 32-bit segment index in the record chain id (source-proof.md section 4.2) while allowing classes with up to 2^33 segments (section 3 under the H2 length bounds), so the index overflows on long classes. This package amends the encoding explicitly (verdict cited as procedural context only): the segment-index field is 33 bits in both the pool record (62-bit base index + 33-bit segment index + 56-bit link = 226 of 256 bits, within the 31 spare bits) and the 128-bit transport chain id (62-bit base index + 33-bit segment index + spare = 95 of 128 bits). The segment directory already holds 2^33+1 eight-byte entries and needs no change, and the 62-bit base index needs no change (at most 2^61.6 bases). The repair costs zero instructions: stores write the same words, loads read the same words, and the candidate re-walk's index extract widens an immediate only (ubfx width 32 to 33, still one instruction); record sizes, memory, and all occupancy and cap analysis are unchanged.
