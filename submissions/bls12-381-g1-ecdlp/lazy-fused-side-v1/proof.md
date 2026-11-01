# Lazy and fused side-output kernel for the powers-of-tau target

## Complete claim and source scope

For the protected BLS12-381 powers-of-tau target, at n=2^20, the supplied paired class-walk Cheon construction has total work at most2^116.309074519993 G1-reference units (C=2000), conditional on the explicitly disclosed H1-H5 of source-proof.md. The exact score is116309074519993/1000000000000. Success is at least0.39 under those hypotheses. Memory remains below2^121.5854 bytes, preprocessing below2^115.3125 units and included in total work, and advice is zero. Only A_1 and A_(2^20) are accessed, with all accesses charged; the same bound W(n) applies for2^20<=n<r. No claim is made for smaller n, and no auxiliary powers are transferred to an ordinary track.

The complete algorithm, uniform-challenge reduction, class sieves, paired walk, planner, bounded transport, collision handling, root recovery, verification and H1-H5 support/limitations are supplied in source-proof.md. That file is the immutable source argument from submission9c1b2d5d-98fb-42a4-99ab-593f6406a59f, commitd0a6e4cfea64902889b71847dbc4409709b073fc. prior-proof.md supplies our published parity-key amendment a38cb931-1c42-4bc2-b209-c429ce92dfb1. Both are included as evidence, not prior verdicts. This document replaces only their side-kernel implementation/pricing, fixed-call accounting and numerical ledger displays; all their target, algorithm, cap, probability and memory obligations remain. The source's old1054/1051 kernel prices and old totals are historical, not additional claims by this package.

The source uses r as protected, d=2^20 dividing r-1, m=(r-1)/d, six-element GLV/negation classes, K=6N^2/m>=619/1250 and minimum pair fraction49/50. N=64243470898890460540588974107897450 is the least integer satisfying the source coverage budget. The scheme and all failure checks are unchanged. Neither the arithmetic gain nor a previous hosted decision proves H1-H5.

## 1. Universal arithmetic transformation

Let p be the protected base-field prime and R=2^384. R>8p. The source side-output formula first forms y1+y2, multiplies by a shared inverse, squares to derive the opposite x coordinate, then multiplies again and subtracts y1 to obtain its y coordinate. The pinned original field bodies and the exact old parity-modified function are supplied as measure/side_lanes.S; lazy_witness.py deterministically derives side_lazy.S.

For nonnegative a,b<2p, ordinary Montgomery reduction returns S=(ab+qp)/R with0<=q<R, hence S<p+4p^2/R<3p/2<2p. The word-serial accumulator also fits the original radix: if T<a+p before a step, then (T+a*b_j+q_j*p)/2^64<a+p, since b_j,q_j<2^64. Thus T<3p<R throughout. For the specialized square, the exact product of a<2p fits the double-width accumulator and the same final Montgomery bound applies. These inequalities justify the wider input representatives; a comparison on finite random data alone would not.

The first product receives a=y1+y2<2p and b<p. Its final13-instruction conditional-p subtraction can be omitted. Its raw result mu<2p is passed directly to the square. A bijective swap of its input/output register banks makes the raw output arrive in the consuming registers without any added MOVs. The square's original final normalization remains and yields the identical canonical square; the x-coordinate formula, storage and all downstream operands are unchanged modulo p and canonical where the original formula requires them. The second multiplication accepts its two representatives below2p, so its raw output S is again below2p.

For the final y coordinate, instead of first canonicalizing S then executing the23-instruction modular subtraction, form the six-word difference D=S-y1. As y1 is canonical, -p<D<2p. The carry branch distinguishes D<0 from D>=0. In the negative branch add p once; although the machine subtraction wrapped modulo R, adding p modulo R gives the unique value D+p in[0,p). In the nonnegative branch conditionally subtract p once, giving the unique value in[0,p). No second subtraction is needed. This replaces13+23 instructions by at most23, saving another13; the negative branch saves an additional5, which is not assumed for pricing.

The canonical y coordinate is therefore exactly identical to the source. The parity orbit key reads the same canonical y and remains identical. Inputs and all retained outputs match, while SP/LR and the original callee-save contract remain intact. The entire loop and caller-visible effects are unchanged except for fewer executed instructions. lazy_witness.py records the precise register/instruction transformation, and native_check.py interprets the actual resulting scalar instruction stream against whole-integer field/point calculations.

## 2. Complete kernel and call accounting

The original parity-modified worst-path body costs1051 instructions. The two proven transformations save at least26 per lane, giving1025. The function has36 fixed instructions; the full ledger reserves a deliberately larger512 per possible invocation including call preparation/return. The side call is issued at most once per round per base; (BASES+2)*LANE_OPS bounds such calls, including the partial terminal bases. It is not assumed that every call has2^18 lanes.

Charge1026 per side output, allocating one instruction beyond the1025 body. The existing padding term is charged at the same updated pair price. Its explicit reserve F*N+PADOPS is62975266199621165091686708678568685, exceeding512*(BASES+2)*LANE_OPS=1066574545623588727886416258793472. Thus the globally charged reserve covers all fixed calls even for short batches. If the actual paired fraction exceeds F=49/50, pairing is cheaper than two root operations, so using the minimum fraction in the full upper ledger remains conservative. No branch-frequency assumption is introduced.

The primary2000 reference charge is retained. In particular it includes batch-inversion amortization and must not be replaced by the cheaper side-point formula. Pair bookkeeping remains53 instructions and root bookkeeping50, including their seven-instruction worst-path parity keys. A pair costs2000+1026+53=3079; a root costs2050. Including the source average record transport81/2, the weighted per-point charge is159021/100. All additional costs below remain charged separately. This is an analytic instruction delta from a pinned source measurement, not a newly measured Callgrind or hardware throughput result.

## 3. Exact complete ledger

The unchanged source dimensions, quotas and capped schedules are evaluated by evidence.py. Every row is an exact integer, including rounding at the original locations. The complete instruction total is205849678999786419940317141650955669379, divided by2000 for W. Preprocessing is103167973916871702166031792963436304095 instructions, included. Its baby-phase transport is charged at45 per record rather than the combined40.5 average, exactly as in the source.

- base interval sieve, characters, assignment: 14778874146939128294669440 instructions.
- candidate delivery and replies: 5770237022568448 instructions.
- candidate re-walks: 1257895200023828017403475787776 instructions.
- candidate tests: 13858116485374300651520 instructions.
- character sieve, class lists, their transfer: 30225869176339564179694288896 instructions.
- comparisons at the per-owner cap: 175617550116043315285975015794671616 instructions.
- gap tables: 127169448742230040386870230497541850 instructions.
- global barriers (wakes, sweeps, releases): 6630756600764919022560032736 instructions.
- listing allowance: 8 instructions per point: 1027895534382247368649423585726359200 instructions.
- owner initialisation: 83076749736557315843464236105728000 instructions.
- padding ops at the pair price: 51310668903513872919393768661057536 instructions.
- per-wave link words, memory barriers: 15209023373991571020430912520192 instructions.
- plan distribution (copy tree): 121734115224496661965030475956224 instructions.
- planner (scans, matching, emission): 487728773955200596712816258318336 instructions.
- planner descendant traversals (adjacency lists): 10161011842454401296372011630592 instructions.
- records: 2N points: 204321219716249198512499985032239207929 instructions.
- root recovery and input distribution: 73872550223301050368 instructions.
- segment restarts (scalar multiplications): 62745545517876083201867020314828720 instructions.
- table strand starts and base points: 7658447738681224922543308800000 instructions.

Memory and all per-node layout assertions are retained: total below2^121.5854 bytes and every node below2^63 bytes. The replacement uses the same register-save footprint and the existing generous code allocation. No free precomputation, future-target amortization or wall-clock parallel discount is introduced. The exponentially large node/memory assumptions are analytical accounting conventions, not feasible hardware resources.

log_bound.py encloses ln(x) by64 rational terms of2*sum(t^(2j+1)/(2j+1)), t=(x-1)/(x+1), and bounds the remaining positive geometric tail. It similarly encloses ln2 and divides in the safe direction. Applied to the exact total/2000, it proves the declared ceiling116309074519993/1000000000000. Decimal116.3090745199924178... is explanatory only. The source rounded memory/preprocessing bounds are retained/recomputed separately and are not new optimization objectives.

## 4. Success and honest limitations

The actual points, canonical stored keys, all routing and candidate decisions, and recovered scalars are identical to the preceding algorithm. Therefore the arithmetic transformation does not change any source success event or cap distribution. H1 is the production coverage allowance, H2 the class/base distribution and runtime bounds, H3 the key/load model, H4 the planner's pair fraction, and H5 the candidate-multiplicity/tail assumptions. Their precise definitions, evidence and failure budgets remain in source-proof.md section7 and its planner sections. In particular the toy coverage deficit is not hidden or called a proof of production coverage. Those hypotheses remain extrapolations requiring review.

The inherited lower floor is1-exp(-0.4952)-0.0005-2^-20-4*2^-40-2^-111>0.390050041, and the declared success is only0.39. This is conditional on the disclosed model, not an empirical production success rate or reviewer-confidence statement. All failures and capped work remain in the complete ledger. The exact uniform target and auxiliary-power limit remain unchanged.

The full combined test reruns the production-parameter integer ledger, planner controls, reduced-order curve recovery, synchronous transport including overflow controls, toy census, and12142 parity-key cases. The independent new native audit adds458 field lanes,32 genuine production-curve point cases, random register contents, read-only input and callee-save checks, loose representatives including mu=p, and both final-subtraction branches. It verifies the exact26-per-lane saving plus5 on each negative final difference. Tests exercise finite cases only; the algebra and accounting arguments above carry the production extrapolation. No full-size attack is run and no previous judgment is used as proof.
