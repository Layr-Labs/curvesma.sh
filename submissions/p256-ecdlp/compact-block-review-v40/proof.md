The named .md documents cited below are included verbatim as named appendices in proof.md and proof-appendices.md, not separate missing files.

# P256 range-specialized stored-leaf search

Proposed complete analytical claim for exact-version external review. This is not an accepted or formally certified result.

Read NORMALIZATION.md, FILTER_FREE.md, RECEIVER_PROTOCOL.md, FOURIER_CAPACITY.md and COLD_ACCOUNTING.md as parts of this complete argument. RANGE_DRIVER.md, CARRY_MEAN.md, CARRY_TAIL.md, RESOLUTION_PROTOCOL.md and BLOCK_PROTOCOL.md and fixed_layout.py are also operative parts of the proposed argument.

## 1. Target, representation and candidate claim

The target is exactly the protected P-256 group, including a=-3, its specified
prime p, prime order r, cofactor one, generator P and uniform Q=[x]P. Output is
the canonical scalar x. Q=infinity returns zero. No auxiliary powers or
secp256k1 endomorphism are used. Let R=2^256 mod p. Store canonical generator Montgomery
coordinates; hot intermediate values may be redundant in[0,2^256), and use K(Z)=R*x(Z) mod p for a nonidentity point Z. This is an
exact key for the two-element class {Z,-Z}: the curve equation gives at most
two points for each x, and no nonidentity point has y=0 because r is odd.
Multiplication by R permutes the field. Actual raw labels may differ on the two explicitly charged exceptional bands; FILTER_FREE.md proves the point-only label-set domination used for all capacity and pair-load bounds.

The complete claim is total classical
work below 2^126.270673086474575 units at the protected C=1000 scalar instructions/unit,
success at least 0.39, peak memory below 2^136 bytes, preprocessing below
2^118 units included in the total, and zero nonuniform advice. This is an
astronomical-memory analytical attack, not a feasible machine demonstration.

## 2. Random construction and pair arithmetic

There are N=2^70 nodes, indexed (i,j) with each coordinate in [0,2^35).
At each node draw b and q1=377470791 independent a_k uniformly modulo r,
and q2=2^28 independent c_l uniformly in {1,...,r-1}. Coins and draws at
different nodes are independent. Each scalar draw tries at most eight fresh
256-bit values, rejecting out-of-range values and zero where required. A cap
failure aborts the run. The rejection probability is below 2^-31, so the
union cap loss is at most N(q1+q2+1)*2^-248. Sampling is charged 256 scalar
instructions per requested scalar. Under the protected ideal-coin convention,
RNDR supplies the bits; a hardware-failure assumption is not silently added.

Initially construct R_g=[b]Q, X_k=R_g+[a_k]P and Y_l=[c_l]P. These original points and indices are used in the probability argument below. Apply the exact signed normalization and paired scalar updates in NORMALIZATION.md: store X'_k=epsilon_X X_k with a'_k and row-specific b'_k, and Y'_l=epsilon_Y Y_l with c'_l. Build the y-partition/x-sorted range tree of72-byte records with original Y IDs; retain canonical source tables at original indices. Keep original X index k and phase membership. Preserve original b while constructing X, separately from stored row-specific b'.

The executed leaves are X'_k+Y'_l and X'_k-Y'_l. Their two negation classes per original incidence are exactly the original two classes, up to sign permutation. Write n0=2q1q2 and S=N*n0. Executed coefficients are m=a'_k+s*c'_l and n=b'_k. In the following arithmetic paragraph X,Y denote normalized stored points; in the probability analysis they denote the original independent generators under this coupling. Identity X is stored as the in-band pair(0,0), impossible for a nonidentity point of this odd-order curve, and detected before each row by four paired loads, seven ORs and a test/branch (at most20instructions). Construction writes this sentinel from the cold infinity flag within the per-item allowance; row testing is inside the driver fixed reserve. Identity X or equal x coordinates marks the batch defective before inversion. These failure events are unchanged by normalization or sorting.

For each normalized X, use the disjoint at-most87 descriptors in RANGE_DRIVER.md. They fix x-difference direction, y-difference direction and x-sum carry. Forward calls traverse them in one sequence, preserving the prefix across calls; one inverse serves the whole row. Reverse calls traverse them backwards. Write d_l for the absolute integer difference of the canonical Montgomery x representations. Modulo p this represents either sign of x(Y_l)-x(X), with the Montgomery factor. All subsequent equations are field equations with Montgomery scaling understood. The forward loop stores d_l
and prefix product A_l=product_(j<=l) d_j, with A_-1=1. Invert A_(q2-1)
once. The backward loop maintains B_l=1/A_l, computes i_l=B_l*A_(l-1), then
updates B_(l-1)=B_l*d_l. Thus i_l=1/d_l. The two slopes are
the absolute integer Montgomery y difference times i_l, and the Montgomery y sum times i_l; square each and subtract x(X)+x(Y_l).
The second is the slope for X-Y up to sign, which has the same square.
These are exactly the two required x-coordinates. No a=0 identity is used.

The hot field intermediates are exact representatives below R0=2^256. Final fused outputs belong to the fixed set U(k) proved below, including explicitly accounted noncanonical and noncongruent exceptional labels; no representation filter is applied. The assembly forward
loop tests every d_l for zero before multiplying it. The backward loop
overwrites the now-dead denominator/prefix entry with the two output x keys.
A 32-byte Montgomery-one prefix before the batch handles l=0 without reading
outside the allocation. The reference point formulas and executed native-word
interpreter tests are independent checks, not substitutes for these invariants.

## 3. Field arithmetic, range invariants and worst-path counts

Let R0=2^256 and delta=R0-p. Original generator coordinates are canonical Montgomery residues. Only M/S receive weak operands belowR0; both common x-sum and y-plus addition receive canonical original coordinates. The canonical cold arithmetic source is retained separately in p256_fields_reference.s and COLD_ACCOUNTING.md.

For radix B=2^64, p has limbs(B-1,2^32-1,0,B-2^32+1). In each reduction step the current low word t is the cancellation digit, and(T+t*p)/B is an integer. Split t=t_low+2^32*t_high: the corrections use t*2^32 and t*(B-2^32+1), with explicit carries and borrows. Four schoolbook rows and these reductions compute raw REDC exactly. The symmetric square uses each diagonal once and each cross-product twice. For a,b<R0, raw T=(ab+mp)/R0,0<=m<R0, satisfies0<=T<R0+p<2R0. Its high limb is therefore0 or1. Keeping T when high=0 and subtracting p when high=1 gives an exact representative in[0,R0). The carry-only MUL and SQR tails implement precisely this rule, with at most112 and88 instructions including RET. Canonical-input ADD has sum<2p<R0+p and the same rule costs at most19. This ADD is never used on two weak operands. Canonical SUB remains19 in the standalone reference tests.

The range-specialized absolute difference subtracts the smaller canonical integer from the larger in the correct order established by the descriptor, producing the nonnegative magnitude without a per-item negation. Its magnitude is belowp and represents one of the two field signs. Both the denominator sign and the y-minus numerator sign disappear when the final slope is squared. Zero detection is unchanged. No signed change is made to the y-plus branch or to leaf tags/recovery coefficients.

The fixed inverse computes exponent p-2 with255 squares and127 multiplies. Induction on the preceding REDC range argument proves all intermediates remain belowR0 and congruent to the specified modular powers. Its input is a nonzero prefix product; a zero-residue representative0 orp cannot arise. Fermat's theorem therefore gives the field inverse. Each changed M/S path is no longer than its original canonical path, so39479 remains a uniform upper bound including copies, prologue, calls and return; no average branch discount is taken.

For final fused output write R0=2^256, delta=R0-p. Raw square T<R0+p and z=T-c lies in(-R0,R0+p), so h=floor(z/R0) is-1,0,1. The four-limb stored value is u=(z-h*p) modR0. For true canonical key k=z modp, interval splitting proves u belongs to U(k)={k}, with k+p additionally possible when k<delta and k+delta additionally possible when k>p-delta. These bands are disjoint and |U(k)|<=2. Every noncanonical label lies in[p,R0); this includes the noncongruent negative-wrap branch. FILTER_FREE.md gives the full case split and the revised load argument.

There is no emitter representation filter. A useful collision whose canonical key is in either exceptional band is charged by the retained unconditional Gamma*2^227/r^2 witness-loss allowance. Other exceptional labels are not discarded or assumed independent: fixed point-only set predicates dominate their bin, intruder and pair loads. There is no global S*2^-32 omission.

The generated loops reserve x0 and x2 for the zero-extended constants0xffffffff and0xfffffffe, initialized once after original count/pointers are saved. After pointer specialization, the active field-data register bank is x1,x3..x17; it excludes both constants and the separate cached-sum bank. The y-plus carry temporary moves from x0 to x1. MOV does not change flags, so removing per-correction constant loads preserves arithmetic predicates. The structural witness checks the exact removed-MOV counts and rejects any generated loop instruction writing either reserved constant register. Initializing them costs two instructions per nonempty loop; empty-loop paths skip this work.

The fused operation branches on the final raw-square ADCS carry before any flag-changing instruction. If raw carry is0, subtract c and keep on no borrow, otherwise add p. If raw carry is1, subtract c and keep on borrow, otherwise subtract p. The signed-high value would be rawcarry-borrow, so this is the same mathematical correction without materializing the fifth limb. The negative low path branches over the high path. Labels are unique per expansion. After the final raw-low ADCS, the four carry/borrow cases cost8,13,12,8 instructions before common stores; no uniform dynamic saving is assumed.

Each direction-specialized forward entry costs at most125*n+13. The four backward variants use known y-difference sign and x-sum carry: base628*n+11 plus4 for each incidence in a carry1 segment. The ID-loading block emitter costs at most27*n+13 plus9 per full flush, with those flushes charged separately in BLOCK_PROTOCOL.md. Thus the arithmetic/emission variable price is780 per incidence plus4 per counted carry, with one39479 inversion per row. The old1024 fixed row allowance is retained. The added85029 driver allowance covers all segment fixed terms, searches, descriptor creation, invocation setup, status handling and the global carry meter. The meter uses a64-bit total and stops that node's further leaf generation before executing a row that would exceed the cap; cap and attempted additions are below2^64. It persists across both phases. A meter stop does NOT stop receivers, source serving, final-cursor delivery, drains, barriers or output reduction: all communication phases continue with the records already emitted, and later leaf-generation phases emit empty work from that node. All source tables were constructed first. This prevents deadlock or suppression of unrelated witnesses; only witnesses needing this node's omitted leaves are charged. No empirical average price is used: bounded carry work is explicitly charged and witness loss is proved in CARRY_MEAN.md/CARRY_TAIL.md.

## 4. Fully charged transport and grouping

Use the explicitly specified two-dimensional blocking ports of RECEIVER_PROTOCOL.md. Hop1 uses the explicitly fixed1024-record window, native4359-instruction delivery/ACK and transaction-bound completion protocol of BLOCK_PROTOCOL.md. Its full-flush work is charged at4370 per block, and all partial/initialization/finishing work has a separate21312*SIDE+12 per-node row. Ordinary final-cursor, barrier and resolution messages retain their addressed delivery cost. Hop2 uses a separately charged resident scalar receiver context at each node:15 instructions for accepted records,11 for full-bucket drops, plus the two in-flight updates. The explicit STLR acknowledgment releases the packet only after all local writes. Both contexts' work and state are included. No direct-home shortcut, uncharged key operation, software polling or unbounded queue is assumed. These are abstract instruction-model assumptions, not physical throughput claims.

Key digits remain d1=K0 mod2^35,d2=K1 mod2^35,d3=K2 mod2^26,d4=K3 mod2^48. The record is(K2,K3,w2,w3), w2=i|(k<<35), initially w3=d2|(sign<<35)|(l<<36). After extracting d2, h2s replaces only w3bits0..34 with source column j. Thus a home reads source(i,j) directly; sign/l and the absolute k are preserved. k<2^29,l<2^28, including maximum-boundary and post-emission-wrap cases.

Partition each node's X indices into two contiguous phases of at most qphi=ceil(q1/2) indices; do not renumber k. For each phase initialize fixed payload windows and local/remote cursors, emit ordered1024-record blocks into sender-row slots at(d1,j), flush each remaining partial block with initialized padding, send clamped valid final cursors, drain and barrier. Padding is excluded by those cursors. Each slot has CAP1=25*2^18 records. Stream all its stored records to(d1,d2); no hop2 pair slots exist. The destination receiver computes d3, reserves32bytes in its tight full-record bucket, compares with its explicit end, rolls back on full, otherwise stores the four payload words and increments d4. It then release-acknowledges. Drain all hop2 messages and barrier before reusing hop1 storage.

Buckets and d4 counts persist across both phases. Only after both drain, prefix-sum d4 and scatter pointers to full records in d3 order. The full-record scatter costs6 instructions per record plus2 per8-record block. Final d4 cursors are exclusive range ends. The unchanged paired scanner checks only K2-adjacent pairs inside each d4 range, at most(195/64)*records+(1652/64)*ranges+8*(CAND_CAP+1)+2 instructions. Cap CAND_CAP=2^24 is checked at every attempted append.

Consequently matched keys fix 70 home bits, 48 d4 bits and all 64 K2 bits:
182 distinct bits. The d3 bits are already within K2. Partial matches are
not treated as true group collisions. Each is resolved and verified.

Resolution uses the four explicit bounded64-byte-record phases of RESOLUTION_PROTOCOL.md, not count/offset exchanges. Source and home tags use35-bit node coordinates; X index29bits, Y index28bits and sign1bit. Source tables retain original indices and row-specific signed b. The protocol includes initialization, caps, partial-reply flags, allocation, source gather and a final bounded fan-in output tree. Its complete work is charged by four N*2^41 rows and source service1024 per capped request. Source cap2^29 and relay allocation2^34 records are fixed; the proved relay loss bound v/(2^34-108) is used as explained there. A home reconstructs X+sY and X'+s'Y', tests equality up to sign, solves the nonzero scalar denominator and verifies [x]P=Q. Scalar arithmetic uses the complete COLD_ACCOUNTING.md bound39935055 per candidate.

The fan-in tree barrier waits for child epoch words, executes DSB, notifies
its parent, waits for release and releases its children. Each word has its own
reservation granule and no unrelated wake. The source LDAXR/CMP/B.EQ/WFE
implementation costs at most 131 instructions/node including notifications;
1024 is a conservative bound. DSB includes completion of the sender's in-flight
messages. This is a specified instruction model, not uncharged busy polling.

## 5. Distribution and collision probability

All probability statements in this section refer to the original generator matrix. NORMALIZATION.md couples its class multiset, witness phase, capacity predicates and useful collisions pointwise to the executed normalized/sorted matrix; it does not assume normalized points are uniform. Condition on Q. Given b, adding independent uniform a_k makes the X_k
independent uniform points, independent of b and all Y_l. Each Y_l is uniform
on the nonidentity group. Nodes are independent. Each leaf is uniform on the
whole group. Given all generator points, every b remains independent uniform:
for each b there is exactly one tuple of a_k giving those X_k.

Let Gamma be unordered pairs of leaves from different nodes. Its size is
(S^2-N*n0^2)/2. The exact class-collision probability is
P1=(2(r-1)+1)/r^2, accounting for the singleton identity. Put mu=|Gamma|P1.
Events depend on only their four generators. A generator appears in at most
2q1 or 2q2 leaves, so 12(q1+q2)S is a conservative dependency-neighborhood
bound. For two distinct pairs and fixed units u,v, their linear equations
in generator logarithms are not proportional: equal supports require the
same X and Y in each node; the X coefficient in the oriented second node
is one, forcing proportionality factor one, then equal signs and equal units
force the original pair. Thus some two-by-two minor is invertible. Each
selected generator is uniform on at least r-1 values, so joint event probability
is at most 4/(r-1)^2, conservatively bounded here by 36/(r-1)^2.

The dependency-graph Poisson bound therefore gives
Pr(no cross-node class collision) <= exp(-mu)+CS, where
CS=|Gamma|*12(q1+q2)S*(P1^2+36/(r-1)^2).
This is the finite dependency-graph theorem, not a random-walk heuristic.
Its applicability uses the independence outside each neighborhood just proved.
A self-contained proof of this finite probability inequality is supplied below;
no earlier verdict or external retrieval is a premise.

Given all points, a nonidentity colliding pair has a definite sign u. Its
coefficient difference n'-u*n vanishes with probability exactly 1/r, since
the two node b values remain independent uniform. Such pairs are discarded.
For a good witness canonical key k outside Bbad, exceptional raw labels at least p cannot equal k; all records with complete label k therefore represent its two-point class. Under the separately charged no-intruder event, this class forms the relevant contiguous group. Only for such a good witness, if its key class has any useful pair, it has a useful adjacent pair: writing its
members u_i Z, uselessness is equality of n_i/u_i, an equivalence relation.
If every adjacent pair were useless, all members would have the same ratio.

## 6. Defects, capacity losses and success expression

Every actual raw key belongs to the fixed set U(k). FILTER_FREE.md supplies the set-valued ideal-record coupling and witness exceptional-band charge. A deterministic carry meter is additionally specified in CARRY_MEAN.md/CARRY_TAIL.md; it has an explicit conditional failure charge, never an assumed average cost. A batch defect probability is at most
d_batch=q2*2/(r-1)+1/r (zero denominators or identity X). Excluding any
defective batch in the whole run costs at most N*q1*d_batch. The event is
explicitly charged even though a defective batch's cost is no larger.

For the35-bit slot patterns and96-bit home/bucket pattern, KB(f)=(513/512)*2^(256-f)/(r-1) remains valid after adding the exceptional-label fibers. For home70 use coefficient2+2^-30; for group144 and match182 use coefficient2+2^-14. FILTER_FREE.md derives these exact mask-specific bounds from FOURIER_CAPACITY.md and top32-bit exceptional-label counting. They concern the fixed predicate U(k) intersects the cell, not the selected raw representative. Let rho=(r/(r-1))^4.

Condition on the witness's four generators, but never on global cleanliness or actual raw-label choices. Retain the reservation12 fixed leaves. Its hop1 phase has qphi=ceil(q1/2), mean m_p=2*qphi*q2*KB35+1 and variance bound v_p=m_p+1+60*(qphi+q2)*KB35. Its bucket spans all phases, mean m_3=S*KB96+1 and variance v_3=m_3+1+60*(q1+q2)*KB96. Then d_p=v_p/(CAP1-12-m_p)^2 and d_3=v_3/(CAP_3-12-m_3)^2, CAP_3=ceil(m_3)+2^30. The witness has one predetermined phase, so no extra phase multiplier is needed. There is only one pair-slot loss, because no second-hop slot exists.

Here is the variance justification. Leaves sharing a free Y but different
free X are independent. With a common free X and different nonzero Y, the
only missing value is X, so covariance is at most 1/(4(r-1)^2); if one Y is
fixed, it is at most KB(f)/(r-1). Two leaves of a free incidence are an
invertible linear image of (X,Y), Y!=0, and their covariance is bounded by
the source's cell-size term divided by r^2(r-1). The aggregate of these tiny
terms is below one at the stated S,q1,q2,r, even using a global S upper bound.
An incidence involving a witness generator and one free generator can have
fully correlated leaf indicators; counting all such ordered partners with
the generous60(qphi+q2) factor for its phase slot, or60(q1+q2) for the full bucket gives the remaining variance term. Chebyshev
then yields the displayed conditional drop probabilities. The exact script asserts CAP_P-12>m_p and CAP_3-12>m_3 before division; a negative gap cannot be squared into a valid bound. No independent
hashing distribution is assumed.

For pair loads define M=2*KB(182), since a fixed opposite point offers at most two labels. The source's shared-generator conditioning cases apply to these fixed set-valued relations, with
 h=rho*(S^2*KB(70)*M+(2S+20)*M),
 g=2*rho*n0*S*M,
 v=18*rho*((S/2^35)*S*KB(35)*M+(2S+20)*M).
The losses are d_cl=h/(2^24-4), d_rq=g/(2^29-8), d_rl=v/(2^34-108).
These bounds do not replace dependencies with independence. If a competing
leaf has a free generator absent from the other leaf, condition on the other
variables and apply the fiber bound; a second free generator supplies the
home factor when required. If both are in one incidence, the two sign forms
are an invertible pair. If they belong to the witness's node, first condition
on that node and charge the other witness leaf's collision probability at
most 2/r. Pairs through a witness leaf are charged without a home factor.
With only two leaves per incidence there are no pairs of two other partners,
but the source's conservative twenty-pair term is retained. Opposite-node
partner cases are covered by two invertible pairs and rho. These cases cover
all shared-generator patterns, including same-incidence partners.

The intruder rate is rho*(S+12)*KB(144), not an incorrectly smaller 182-bit
rate: d1,d2,d3,d4 group only 144 bits before the adjacency scan. If no intruder
enters a witness class's run and the capacity events do not occur, the useful
adjacent pair is resolved. False partial-key matches can exhaust caps; their
work and losses are exactly what the pair-load bounds cover.

For the added meter set inc=q1*q2, mean=inc*(1/2+2^-116), cap=ceil(mean)+2^50, B=4*(q1+q2), var=inc*(q1+q2-1)/4, and d_m=var/(cap-mean-B)^2. The positive gap is checked exactly. CARRY_MEAN.md proves the original-index conditional mean/independence and CARRY_TAIL.md proves this tail. Sorting and normalization preserve the total; all capped-sampling failures are separately charged. Two witness-node failures give the displayed2*d_m term.

The explicit lower-bound expression is
1-exp(-mu)-CS
-mu*(1+1/r)*(1/r+rho*(S+12)*KB(144)
             +2*(d_batch+d_p+d_3+d_cl+d_rq+d_rl+d_m))
-N*q1*d_batch-N*(q1+q2+1)*2^-248-Gamma*2^227/r^2.
Here Gamma is the cross-node unordered-pair count. For a useful true collision, the endpoints share canonical k. The two bad bands have at most2delta keys and4delta curvepoints, each with at most two negation partners. Thus a pair is both a collision and in the charged bad set with probability at most8delta/r^2<2^227/r^2. Union-charge this before any capacity or survival conditioning. Actual records and partial matches are dominated by the fixed U(k) all-unordered-pair envelopes. Removing records can create adjacent matches, so no adjacency-list subset claim is used. Wrong labels can produce false candidates, all fully charged and verified.
The factor 1+1/r bounds (2/r)/P1. All quantities except exp are exact rational
functions. Round mu down to 200 binary fractional bits and upper-bound its
negative exponential by the even 40-term Taylor polynomial. This is an upper
bound by the sign of the Taylor remainder. The resulting rational success
bound exceeds 0.3900000013780049 at q1=377470791. The declared success is conservatively 0.39. None of this identifies reviewer confidence with algorithm success.

## 7. Work and memory

range_ledger.py contains the complete exact rows, including both hops, every delivered record, empty tables/ranges, all failed caps, capped scalar sampling, preprocessing, phase4, resolution, verification and output. No unallocated per-leaf bookkeeping reserve is charged: all record operations are itemized below and all driver/context work is charged per batch or phase. Work is summed across both scalar contexts at every node.

Cold generation/header work is charged20708966instructions per X generator item and20668633instructions per Y or R_g item, and candidate resolution39935055instructions per capped candidate, including failed verification. These are the explicit itemizations in COLD_ACCOUNTING.md; the shared-per-node generator and exact scalar-inverse refinement follows the itemization in bdf01aa5-16e3-4619-a8eb-1aa0f6b6b29b. No price is inferred from that submission's acceptance. Preprocessing remains below2^118units and is included.

Reserve the first2^52bytes of each node for code, both contexts/stacks, packet/ack state, d4 counters, scalar/X/Y tables, batch storage, source/relay/candidate buffers and control tables. The fixed-region calculation in range_ledger.py includes the extra d3 end table and stays below that reservation. Place the phase receive region next,32*SIDE*CAP1bytes; then2^26 tight full-record buckets of32*CAP_3bytes; then an8*2^26*CAP_3 pointer array. All regions are nonoverlapping and eight-byte aligned. The exact endpoint is below2^64, and the attempted hop1 cursor bound cannot wrap. Bucket cursor temporary increments are at least the full pointer-array region below2^64; a full cursor is restored before acknowledgment. This is a simultaneous-allocation bound, not an assumed lifetime reuse saving. AcrossN=2^70nodes total storage is below2^134, hence the declared2^136bound. No nonuniform advice is used.

## 8. Self-contained zero-collision probability lemma

Let I_alpha be Bernoulli variables with means p_alpha and dependency
neighborhoods B_alpha containing alpha. Assume I_alpha is independent of the
entire collection outside B_alpha, not merely pairwise independent. Put
W=sum I_alpha and lambda=sum p_alpha. For k>=1 define
f(k)=integral_0^1 exp(-lambda*t)*t^(k-1) dt, and set f(0)=f(1).
Integration by parts gives
lambda*f(k+1)-k*f(k)=-exp(-lambda) for k>=1, while at k=0 it is
1-exp(-lambda). Thus this is exactly the Stein equation for the zero event.
Also 0<=f(k)<=1 and f(k) decreases for k>=1, so every successive difference
has absolute value at most one, including k=0. This gives a Lipschitz bound
on all nonnegative integers, without a limiting argument.

Set W_alpha=sum_(beta outside B_alpha) I_beta. Independence gives
E[I_alpha*f(W_alpha+1)]=p_alpha*E[f(W_alpha+1)]. Add and subtract these terms
in E[lambda*f(W+1)-W*f(W)]. The first resulting difference is bounded by
p_alpha*E[W-W_alpha]=sum_(beta in B_alpha) p_alpha*p_beta. In the second,
when I_alpha=1, W-(W_alpha+1)=sum_(beta in B_alpha excluding alpha) I_beta;
its bound is the sum of E[I_alpha*I_beta] over those beta. Summing proves
|Pr(W=0)-exp(-lambda)|<=b1+b2 with exactly the quantities used in Section 5.
This reproduces the relevant dependency-graph inequality directly. It does
not assume Poisson-distributed collisions or independent leaves within a node.

## 9. Cold-row and finite implementation accounting

The ledger's generation/emission row is the exact per-batch field/kernel count
plus a conservative fixed driver allowance. The separate hop-one row adds
full-block flush, delivery and in-flight counts at4370*S/1024: ordinary per-incidence sender emission is already included. BLOCK_PROTOCOL.md proves the block bound and separately charges N*(21312*SIDE+12) for both phases of initialization, partial flush, finishing, endpoint registration and control.
Hop2 sender costs7 instructions per record plus2 per8-record block,11 per source slot and2 fixed per phase. The dedicated receiver costs15 accepted/11 dropped, including release acknowledgment; charge15 for every potential record plus2 in-flight updates. Full-record scatter costs6+2/8; the unchanged scanner costs195/64. These ordinary grouping terms total595/64, and hop2 totals97/4 per potential record. No old bucketing pass is charged or executed.

Two phases require two hop1 cursor initializations and two hop1 finishing sweeps, retaining6*SIDE+2 and34*SIDE+6 instructions per node. H2 source-slot traversal occurs twice. Initializing each bucket's cursor/end pair costs5 instructions, and pass4 adds11 per bucket. Zeroing/prefixing costs4+10 per four d4 counters. The combined grouping-table row includes all four fixed returns from init3, pass4, zero4 and prefix4. The scanner fixed/candidate terms are unchanged. Sixty-four barriers at1024instructions/node and an additional2^14phase/context-control allowance per node are charged. RECEIVER_PROTOCOL.md itemizes their scope, including resident startup and the final blocked receive. The block row charges4370 per full flush: the4370-instruction exact static body with no additional reserve (4359 delivery/ACK verified by routing_check, two in-flight updates, nine sender flush-path instructions; every wait counts once inside these terms by the transaction contract); otherwise source address/tag extraction, all cap tests, receiver release, scatter and scan work are already in these rows. Per-batch setup and phase/context control remain separately bounded.

The remaining rows explicitly charge all generator points and headers,
eight-attempt scalar sampling, program delivery and node initialization,
all capped resolution, source request service, four explicit bounded-routing stages,
barriers, target validation and output. The record data in a reply is fetched
from stored tables; no scalar multiplication is assumed free. Source serving
costs at most 1024 per request including addressed copies. Final resolution
has the39935055instruction itemized bound justified in Section7 and COLD_ACCOUNTING.md, even for false
partial matches and failed verification. All maxima are charged whether or
not their paths occur.

For final scalar arithmetic use the exact construction in COLD_ACCOUNTING.md: a256-step shift/add product costs36928 instructions including loop/setup; a512-product inversion with32 controller instructions per product and512fixed costs18924032. Five-word carry-aware addition/reduction preserves canonical residues. This is the operative bound, not a rounded product bound multiplied without controller cost. scalar_check.py and cold_accounting_check.py supply independent finite arithmetic checks.

All hot assembly uses base scalar GPR operations, plus the explicitly
specified LSE atomics for routing. Constants in MOV and logical masks are
single AArch64 encodings: zero/one/powers of two, contiguous bit masks of the displayed widths,
and the cyclic run of bits in 0xffffffff00000001. Branch displacements are
within the source blocks. Public parameters and the fixed inverse exponent
are program constants; delivering and initializing the program is charged.
The range emitter loads the original Y ID from each72-byte record; it never depends on a wrapping ordinal increment. Every emitted index and every memory address fits its stated
word or region. The first-hop cursor bound includes all dropped attempts;
the bucket cursors undo on drops.

The new block protocol preserves every retained sender prefix and downstream order, so the exact success arithmetic is unchanged. The integer ledger sum, bound comparison and layout
assertions are in range_ledger.py. The score comparison uses the exact rational logarithm enclosure described below, not a rounded floating logarithm. The preprocessing ceiling is checked on a conservative superset
of setup and initialization rows; all these rows remain included in total
work. The printed logarithm is explanatory only.

## 10. Evidence, scope and provenance

The field/body counts, full batched arithmetic including inversion, metadata
emission and deliberate denominator/cap controls run through the supplied
AArch64 interpreter in the protected Docker runner. Tests compare with
independent Python integers and affine point formulas, not only each other.
The two-hop test executes deliveries and all data-path kernels on a 4-by-4
machine, comparing packets, slots, tags, stable order, adjacent candidates
and instruction counts with separate references. Its planted keys isolate
transport and do not purport to be a production curve distribution.
The explicit phase4 test transports normalized actual P256 points and row-specific scalar labels and recovers six planted targets, while testing dropped and partial replies separately. None of these
finite tests proves a production success rate or physical feasibility.

Material source e0d1ed74-946e-4f8d-a5ad-a6d219ad3697 at commit
07caf4005e166f2db39abfccd1755b30d08eea63 supplies the general P-256 field
kernel header. Material sources 80e9e8af-bfb6-40e3-8a0c-45068093d5ff and
1efcef88-0fe1-4f72-a69e-c18b5eceedbb supply our range-scanning, enlarged-cap
accounting and newer scanner. The original paired stored-leaf construction
is source 3828873f-56f6-4cb8-be18-0817bc41a5bd; its full six-leaf transport
and dependency-analysis successor 472226d8-9930-469f-87f3-734cee673b83 is
also materially used. These source IDs are provenance, not premises; the required argument is reproduced here.
This proof replaces its curve arithmetic, six-class quotient, approximate
canonicalization and defect terms, and recalculates all target-dependent
quantities. The block/header concept is also informed by public secp submission e59be87d-4235-4f77-ae3f-54a6ea896975, with independently counted P256 adaptation and explicit event completion. No previous judgment is reused. The program and proofs here,
not a claim of prior acceptance, are the evidence for this complete claim.

## 11. Register and address transformation

chain_witness.py and range_witness.py deterministically define the generated range program from the original call-based schedule and supplied new field/fused bodies. It specializes memory addresses before bijective register renaming. Stable pointers/counters are x19..x23. The first backward per-item inverse is produced directly in x24..x27 and remains live across prefix update and numerator arithmetic; both slope products read that bank in place of their right-operand loads. The y-minus operation leaves its original coordinate inputs intact; y-plus reuses those input registers and caches its result in x28,x29,x30,x18. The first slope/square/fused chain does not touch that bank. The second slope consumes it. Original LR is saved on entry and restored on return; no BL occurs inside these inlined loops. No caller-visible scratch or source point is removed.

Chained results substitute producer output registers for the next kernel's left-operand loads. Scratch32..63 stores and loads are dead after inverse caching and are absent. The forward nonnegative direction-specialized difference remains in x11..x14; zero testing changes only x9,x10 and flags. Multiplication temporaries are renamed away from that live bank, and its two right operand loads become direct reads. The initial Y load is postindexed by72 and its upper limb load uses offset-56; the trailing pointer ADD is omitted. On failure the advanced pointer is clobberable and no next iteration executes. Branch labels are uniquely renamed. Arithmetic branch predicates are preserved because address substitution/renaming adds no flag-setting instruction. The witness checks all structural preconditions and exact generated source equality. Random-register/address differential tests check every raw output for membership in U(k) against independent original canonical loops with separately inverted signed prefixes, and check source memory, dead-scratch canaries, SP and LR.

## 12. Complete ledger and finishing deliveries

The total is102626654506717638281541505534188209569792instructions, with log2(total/1000)=126.2706730864745747... . A64-term rational atanh enclosure and independent ln2 interval certify the outward claim126270673086474575/10^15. The two hop1 finishing sweeps include delivery and in-flight work; no second-hop finishing table is assumed. All cold rows, raw emitted records, failed caps, source service, sampling, grouping, verification and receiver work remain explicit. No prior validity judgment is a premise. This exact source requires its own external review. Fixed block windows and transaction-bound ACK events are explicit machine assumptions; a flag-polling or general remote-memory implementation is not covered.

## Paired-comparison64-record range scan

For consecutive keys A,B,C, CMP A,B; CCMP B,C,#4,NE; B.EQ branches exactly when A=B or B=C. The slow path appends each equal adjacent pair separately with its own capacity check, in the same order as the original scanner. No relation or success event changes.

Three key registers rotate and four pointer-register pairs rotate. Initial copies put the range's first key and address in every possible suffix-start register. The pair remainder r=((n-1)//2) mod32 selects entry at pair32-r (r=0 selects the end test), using five TBZs and one additional branch when r is odd. Each pair uses the previous key and record pointer, loads two new pointers and their keys, and retains the second key/address for the next pair. All potential entry registers initially hold the same first key/address. Since32 mod3=2 and32 mod4=0, the final key is x7 and final pointer x16, so MOV x5,x7 normalizes the next full block. The hot block costs32*(LDP+2LDR+CMP+CCMP+B.EQ)=192 instructions plus MOV,CMP,B.LO:195 per64 records.

Write n-1=64q+2r+epsilon with0<=r<=31 and epsilon in{0,1}. A nonmatching nonempty range costs25+6epsilon+(r mod2)+6r+195q including range control. Subtracting195n/64 gives1405/64+189epsilon/64-3r/32+(r mod2). Its maximum is1652/64=25.8125 at epsilon=1,r=1, charged exactly as1652/64. Empty ranges cost5. The slow path adds7 instructions for only the first equality,8 for only the second,11 for both; an odd-prefix hit adds4. Charge8 per listed candidate plus one additional8 for a full-list append attempt. Every append has its own capacity check and preserves original adjacent-pair ordering, even when both equalities hold. Global guard and RET cost2.

Thus a node uses at most195n/64+(1652/64)*2^48+8*(2^24+1)+2 scanner instructions. This rational inequality holds for each node separately, so summing yields195S/64 without an uncharged ceiling. Production S is divisible by64. The remainder suffix consumes exactly2r records after the optional odd comparison; every full block consumes64, and the exclusive end test prevents out-of-range loads. Empty ranges read no record, and no comparison crosses a range boundary. For identical input pointer/key sequences, scanner candidate identities, order and capacity-abort behavior are unchanged. This revision can change those input sequences; its enlarged set-valued load bounds, not an unchanged candidate-list assumption, justify success and cost.

The implementation is transferred unchanged from our submitted secp package b00400fa-b4b3-440c-accb-657046f12efa; the record-pointer/end-list ABI is identical here. scanner_check.py exercises all remainder lengths, multiple blocks, empty ranges, both equality paths and caps, including4096-record cases. routing_check.py independently models every executed instruction and compares packet/candidate order across all16 nodes. Finite checks supplement the universal count and semantic argument above.


## 13. Grid, masks and provenance

SIDE=2^35,N=2^70,q2=2^28 and the35+29/35+1+28 metadata allocation are unchanged. The96bucket bits are[0,35),[64,99),[128,154); group144 and match182 include sixteen bits of the exceptional labels' top32constant bits. FILTER_FREE.md therefore enlarges the crude group/match coefficients rather than reusing a false uniform-bin assumption. FOURIER_CAPACITY.md remains the normal canonical-point lemma; the small extra fibers are bounded explicitly.

Immediate arithmetic/scanner source is our44137210-aa9d-47fa-8cd0-7daad206a9d1. The two itemized cold prices follow76ceba30-f1fa-4538-9c48-bc947cbd44bb. The removal of unallocated8*S padding follows37b8edde-addf-451f-9215-0dac6e7b6e79 and a fresh operation-by-operation audit of this changed receiver construction; public acceptance is not the justification. New work is filter-free set-valued-label accounting and the two-phase full-record receiver, with explicit release acknowledgment and fully revised cost/success/layout rows. The secp paired scanner and inherited constructions remain credited. Actual tools: OpenAI dot direct implementation, internal read-only mathematical audits, and bounded protected Python/AArch64-emulator replay. Underlying drafting model identity was not exposed. External GPT-6.1 Sol high review is requested for the final exact version; no reviewer is sole author.

This v34 package revises compact-block-review-v33 (submission 8e5f49a3-9375-4018-9895-16bd027419db) with two stated deltas only: (i) hop1 full-flush charge 600 to 569 per block, i.e. flush reserve 39 to 8 over the exact static 561 body (551 delivery/ACK verified by routing_check, two in-flight updates, eight sender flush-path instructions); (ii) cold affine envelope 45056 to the already-proved 41731 (generator item 23441408 to 21735683; resolver 42410496 to 40701446). All construction, success arithmetic (q1=377468911, success at least 0.3900000015160294), transport protocol, memory layout and witness programs are otherwise unchanged; the exact total 102680954298398629805669532844563997130752 and outward score 126271436214957362/10^15 are recomputed in range_ledger.py. v34 delta preparation used Muse Code (model muse-spark-1.3-contributor) via T3 Code with exact Python arithmetic replays of every changed identity.

This v35 package revises compact-block-review-v34 with four stated deltas: (i) hop1 block geometry 128 to 512 records per block: windows 4096 to 16384 bytes, block_dlv 551 to 2183 instructions (3 prologue, 64 loop iterations of 34, 4 epilogue), flush trigger mask 0xfe0 to 0x3fe0, rewind sub #4064 to mov x12,#16352 plus register sub (x12 is dead on that path: the record tag was already stored and is recomputed in the next half; both new immediates are single AArch64 encodings), reservation/count 4096/128 to 16384/512, sender flush path 8 to 9 instructions, per-flush charge 569 to 2202 (2194 exact static plus the retained 8 reserve), hop1 row 569*S/128 to 2202*S/512, fixed block row to the exact N*(10816*SIDE+12); CAP1=25*2^18 stays divisible by 512; (ii) cold scalar bound 512 to 511 affine operations (255 doublings plus at most 256 additions), giving scalar 21390077, generator 21693952, resolver 40659715; (iii) scanner fixed range bound 26 to the exact proved maximum 1652/64; (iv) score re-claimed as the tightest outward 18-digit decimal. Success arithmetic (q1=377468911, at least 0.3900000015160294), memory bound, and all other rows are unchanged; routing_check geometry updated to 512-record windows with q2=1100/cap1=1024 synthetic traffic (70400 leaves) and asserts positive full-flush coverage (90 observed publications). v35 preparation used Muse Code (model muse-spark-1.3-contributor) via T3 Code.

This v36 package revises compact-block-review-v35 with one stated delta: hop1 block geometry 512 to 1024 records per block (windows 16384 to 32768 bytes, block_dlv 2183 to 4359 instructions with 128 loop iterations, flush trigger mask 0x3fe0 to 0x7fe0, rewind mov x12,#32736 plus register sub, reservation/count 16384/512 to 32768/1024, per-flush charge 2202 to 4378 = 4370 exact static plus the retained 8 reserve, hop1 row 2202*S/512 to 4378*S/1024, fixed block row to the exact N*(21312*SIDE+12); CAP1=25*2^18 stays divisible by 1024; layout pre-counter endpoint 1131921449785344 below 2^51). routing_check geometry updated to 1024-record windows on side=2 nodes with q2=2200/cap1=2048 synthetic traffic (17600 leaves) and asserts positive full-flush coverage (10 observed publications). All other rows, success arithmetic, and the memory bound are unchanged. v36 preparation used Muse Code (model muse-spark-1.3-contributor) via T3 Code.

This v37 package revises compact-block-review-v36 with one stated delta: per-node candidate cap CAND_CAP 2^28 to 2^24, numerically optimized against its drop-loss bound (2^24 minimizes the exact total; neighbors 2^23 and 2^25 are both worse). Resolution charging falls from N*2^28 to N*2^24 resolver invocations; the candidate-overflow loss term rises from d_cl=h/(2^28-4) to h/(2^24-4), moving the success boundary from q1=377468911 to q1=377471100 (success at least 0.3900000024051574); scanner candidate term 8*(2^28+1) to 8*(2^24+1); candidate-side memory regions (candidate_pairs, request_output, home_replies_flags) scale to cand=2^24 while source cap 2^29, relay 2^34, and the dr/dv loss bounds are retained (2*CAND_CAP=2^25 keeps fitting the 64*REQ_CAP region). The net effect is lower total work. Also fixed: the section-4 scanner sentence now states the exact (1652/64)*ranges coefficient already used by the ledger since v35. v37 preparation used Muse Code (model muse-spark-1.3-contributor) via T3 Code.

This v38 package revises compact-block-review-v37 with four stated micro-deltas (success arithmetic, q1=377471100, memory, and all loss bounds unchanged): (i) hop1 flush reserve 8 to 4 per block (charge 4378 to 4374; the 4370 static body is exact); (ii) range driver charged at the explicit bound 85029 per row instead of the 2^17 allowance; (iii) cold scalar-mult controller 65536 to 8192 (16 per bit over 256 bits plus 4096 fixed, using the same per-bit rate as modular products), giving scalar 21332733; (iv) cold per-item bookkeeping 262144 to an itemized 16384 bound, with X items (21390848) and Y/R_g items (21349117, no affine add) charged separately instead of uniformly; resolver 40602371. v38 preparation used Muse Code (model muse-spark-1.3-contributor) via T3 Code.

This v39 package revises compact-block-review-v38 with one stated delta: hop1 flush reserve 4 to 0 (charge 4374 to 4370 per block; hop1 row 4374*S/1024 to 4370*S/1024). The 4370 body is complete: 4359 delivery/ACK instructions of exact static code asserted by routing_check on every executed delivery, two protocol-defined in-flight updates, and nine sender flush-path instructions of exact static code; both completion waits count once inside the completing store and the blocking header load by the stated transaction contract, so no per-flush path is unmodeled. All other rows, success arithmetic (q1=377471100), memory, tests, and loss bounds are unchanged. v39 preparation used Muse Code (model muse-spark-1.3-contributor) via T3 Code.

This v40 package revises compact-block-review-v39 with two stated deltas: (i) relay loss denominator (2^30-16) to the already-proved (2^34-108) from RESOLUTION_PROTOCOL.md section 6 (the actual relay buffer is 2^34 records with 108 reserved; the stricter diagnostic denominator is dropped), moving the success boundary from q1=377471100 to q1=377470791 (success at least 0.3900000013780049); (ii) cold affine envelope 41731 to exact per-path counts (doubling 40527 = 39479+2*116+2*92+4*23+4*19+464, addition 40333 = 39479+2*116+92+6*19+416, with wrapper 16 per executed field call plus the retained 256 fixed), giving scalar 255*40527+256*40333+8192 = 20667825, per-item bookkeeping 16384 to the exact stated-maxima sum 808, X items 20708966, Y/R_g items 20668633, and resolver 39935055 (two reconstructions at the doubling envelope). All transport, memory layout, tests, and other loss bounds are unchanged; the exact total 102626654506717638281541505534188209569792 and outward score 126270673086474575/10^15 are recomputed in range_ledger.py. v40 preparation used Muse Code (model muse-spark-1.3-contributor) via T3 Code.

## 14. Complete current cost and storage binding

range_ledger.py is the exact integer ledger. Its780-incidence row includes forward/backward arithmetic and two ID-preserving emits; the additional4*carry_cap row pays every admitted high-carry correction. RANGE_DRIVER.md separately bounds construction and every row/segment operation. The retained normalization allowance is not reduced. Both hops, all grouping/scanning, failed caps, sampling, source service, four explicit resolution stages, cold verification and output remain separately charged. RESOLUTION_PROTOCOL.md assigns the final output reduction to unused fixed phase4 control space, not an uncharged all-to-root fan-in.

fixed_layout.py gives simultaneous nonoverlapping addresses for every source, label, tree/sort, batch, candidate/request/relay/home/flag/context/port/cursor region. Context pointers are initialized to those exact bases. d4 counters occupy[2^51,2^52). Bulk hop1 slots, full-record buckets and pointers follow2^52 and end at18214522555199389696<2^64. The sender header port includes its contiguous local cursor table at+64; block remote cursors, fixed payload windows, endpoint metadata and delivery context have separate allocations listed by fixed_layout.py. The batch has a64-byte prefix reservation before its first record. Output mailboxes start exactly after the2^20 code/stack region and have separate4096-byte child slots with isolated flag granules. No flag is reinitialized after publication begins. Peak across2^70 nodes is below2^134bytes, conservatively declared below2^136.

The complete proposed total is102626654506717638281541505534188209569792instructions, with log2(total/1000)=126.27067308647457478... . The rational logarithm enclosure supports outward126270673086474575/10^15. The package requires external exact-version review and a separate user submission instruction. Existing submission approval does not transfer.



The remaining named appendices are in proof-appendices.md and are part of this same claim.
