The named .md documents cited below are included verbatim as named appendices in proof.md and proof-appendices.md, not separate missing files.

# P256 range-specialized stored-leaf search

Proposed complete analytical claim for exact-version external review. This is not an accepted or formally certified result.

NORMALIZATION.md, FILTER_FREE.md, RECEIVER_PROTOCOL.md, FOURIER_CAPACITY.md, COLD_ACCOUNTING.md, RANGE_DRIVER.md, CARRY_MEAN.md, CARRY_TAIL.md, RESOLUTION_PROTOCOL.md, BLOCK_PROTOCOL.md and fixed_layout.py are operative parts of this argument.

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
work below 2^126.206384844364137 units at the protected C=1000 scalar instructions/unit,
success at least 0.39, peak memory below 2^136 bytes, preprocessing below
2^118 units included in the total, and zero nonuniform advice. This is an
astronomical-memory analytical attack, not a feasible machine demonstration.

## 2. Random construction and pair arithmetic

There are N=2^70 nodes, indexed (i,j) with each coordinate in [0,2^35).
At each node draw b and q1=377468911 independent a_k uniformly modulo r,
and q2=2^28 independent c_l uniformly in {1,...,r-1}. Coins and draws at
different nodes are independent. Each scalar draw tries at most eight fresh
256-bit values, rejecting out-of-range values and zero where required. A cap
failure aborts the run. The rejection probability is below 2^-31, so the
union cap loss is at most N(q1+q2+1)*2^-248. Sampling is charged 256 scalar
instructions per requested scalar. Under the protected ideal-coin convention,
RNDR supplies the bits; a hardware-failure assumption is not silently added.

Initially construct R_g=[b]Q, X_k=R_g+[a_k]P and Y_l=[c_l]P. These original points and indices are used in the probability argument below. Apply the exact signed normalization and paired scalar updates in NORMALIZATION.md: store X'_k=epsilon_X X_k with a'_k and row-specific b'_k, and Y'_l=epsilon_Y Y_l with c'_l. Build the y-partition/x-sorted range tree of72-byte records with original Y IDs; retain canonical source tables at original indices. Keep original X index k and phase membership. Preserve original b while constructing X, separately from stored row-specific b'.

The executed leaves are X'_k+Y'_l and X'_k-Y'_l. Their two negation classes per original incidence are exactly the original two classes, up to sign permutation. Write n0=2q1q2 and S=N*n0. Executed coefficients are m=a'_k+s*c'_l and n=b'_k. In the next arithmetic paragraph X,Y denote normalized stored points; in the probability analysis they denote the original independent generators under this coupling. Identity X is stored as the in-band pair(0,0), impossible for a nonidentity point of this odd-order curve, and detected before each row by four paired loads, seven ORs and a test/branch (at most20instructions). Construction writes this sentinel from the cold infinity flag within the per-item allowance; row testing is inside the driver fixed reserve. Identity X or equal x coordinates marks the batch defective before inversion. Normalization and sorting do not change these failure events.

For each normalized X, use the disjoint at-most87 descriptors in RANGE_DRIVER.md. They fix x-difference direction, y-difference direction and x-sum carry. Forward calls traverse them in one sequence, preserving the prefix across calls; one inverse serves the whole row. Reverse calls traverse them backwards. Write d_l for the absolute integer difference of the canonical Montgomery x representations. Modulo p this represents either sign of x(Y_l)-x(X), with the Montgomery factor. All later equations are field equations with Montgomery scaling understood. The forward loop stores d_l
and prefix product A_l=product_(j<=l) d_j, with A_-1=1. Invert A_(q2-1)
once. The backward loop maintains B_l=1/A_l, computes i_l=B_l*A_(l-1), then
updates B_(l-1)=B_l*d_l. Thus i_l=1/d_l. The two slopes are
the absolute integer Montgomery y difference times i_l, and the Montgomery y sum times i_l; square each and subtract x(X)+x(Y_l).
The second is the slope for X-Y up to sign, which has the same square.
These are exactly the two required x-coordinates. No a=0 identity is used.

Hot field intermediates are exact representatives below R0=2^256. Final fused outputs belong to the fixed set U(k) proved below, including explicitly accounted noncanonical and noncongruent exceptional labels; no representation filter is applied. The assembly forward
loop multiplies every d_l with no zero test; defective rows complete as bounded garbage. The backward loop
overwrites the now-dead denominator/prefix entry with the two output x keys.
A 32-byte Montgomery-one prefix before the batch handles l=0 without reading
outside the allocation. The reference point formulas and executed native-word
interpreter tests check these invariants independently; they do not replace them.

## 3. Field arithmetic, range invariants and worst-path counts

Let R0=2^256 and delta=R0-p. Original generator coordinates are canonical Montgomery residues. Only M/S receive weak operands belowR0; common x-sum and y-plus addition receive canonical original coordinates. Canonical cold arithmetic follows COLD_ACCOUNTING.md.

For radix B=2^64, p has limbs(B-1,2^32-1,0,B-2^32+1). In each reduction step the current low word t is the cancellation digit, and(T+t*p)/B is an integer. Split t=t_low+2^32*t_high: the corrections use t*2^32 and t*(B-2^32+1), with explicit carries and borrows. Four schoolbook rows and these reductions compute raw REDC exactly. The symmetric square uses each diagonal once and each cross-product twice. For a,b<R0, raw T=(ab+mp)/R0,0<=m<R0, satisfies0<=T<R0+p<2R0. Its high limb is0 or1. Keeping T when high=0 and subtracting p when high=1 gives an exact representative in[0,R0). The carry-only MUL and SQR tails implement precisely this rule, with at most110 and88 instructions including RET for the new standalone MUL and retained SQR. Canonical-input ADD has sum<2p<R0+p and the same rule costs at most19. This ADD is never used on two weak operands. Canonical SUB stays19 in the standalone reference tests.

The range-specialized absolute difference subtracts the smaller canonical integer from the larger in the descriptor's order, producing the nonnegative magnitude without a per-item negation. Its magnitude is belowp and represents one of the two field signs. Both the denominator sign and the y-minus numerator sign disappear when the final slope is squared. No zero test remains; defective rows complete as bounded garbage within the same envelopes, detailed below in this section. The y-plus branch and leaf tags/recovery coefficients are unchanged.

Direction specialization is machine-checked limb by limb, not by eyeballing one long line. The pair_fwd calling convention passes the Y stream in x1 (copied to x20) and the X point in x3 (copied to x22), so the inlined four-limb subtraction reads the Y minuend into x3-x6 and the X subtrahend into x7-x10. The unreversed pair_fwd_hi therefore computes Y-X and is selected exactly when Y.x>=X.x; pair_fwd_lo swaps all four limbs (chain_witness.py normalize_transform, one replace per limb) to compute X-Y and is selected exactly when Y.x<X.x. range_witness.assert_limb_specialization asserts on the generated program that the lo section holds exactly the four swapped limbs and none of the forward ones, the hi section exactly the four forward limbs and none of the swapped ones, and each of the four backward y-partition variants exactly its four expected limbs. The probe additionally executes the generated pair_fwd_lo on X.x=2^128 with Y.x=1 (plus its mirror and eight further high-half/boundary operand pairs) and asserts the stored denominator equals the exact absolute difference in every case (2^128-1 for the headline pair).

The fixed inverse computes exponent p-2 with255 squares and127 multiplies. Induction on the preceding REDC range argument proves all intermediates stay belowR0 and congruent to the specified modular powers. The fixed inverse instruction sequence is retained byte-for-byte after macro expansion. It calls external k_mul/k_sqr; the experiment links the changed hot routines. Its conservative39479 allowance is retained without crediting their smaller counts. Canonical cold point operations instead use the separately specified namespaced canonical routines. On non-defective rows its input is a nonzero prefix product and no zero-residue representative0 orp arises; defective rows complete as bounded garbage within the same envelopes, detailed below in this section. Fermat's theorem therefore gives the field inverse. Each changed M/S path is no longer than its original canonical path, so39479 stays a uniform upper bound including copies, prologue, calls and return; no average branch discount is taken.

For final fused output write R0=2^256, delta=R0-p. Raw square T<R0+p and z=T-c lies in(-R0,R0+p), so h=floor(z/R0) is-1,0,1. The four-limb stored value is u=(z-h*p) modR0. For true canonical key k=z modp, interval splitting proves u belongs to U(k)={k}, with k+p additionally possible when k<delta and k+delta additionally possible when k>p-delta. These bands are disjoint and |U(k)|<=2. Every noncanonical label lies in[p,R0), including the noncongruent negative-wrap branch. FILTER_FREE.md gives the full case split and the revised load argument.

There is no emitter representation filter. A useful collision whose canonical key is in either exceptional band is charged by the retained unconditional Gamma*2^227/r^2 witness-loss allowance. Other exceptional labels are neither discarded nor assumed independent: fixed point-only set predicates dominate their bin, intruder and pair loads. There is no global S*2^-32 omission.

The generated loops reserve x0 and x2 for the zero-extended constants0xffffffff andMC=0xffffffff00000001, initialized once after original count/pointers are saved. After pointer specialization, the active field-data register bank is x1,x3..x17; it excludes both constants and the separate cached-sum bank. The y-plus carry temporary moves from x0 to x1. MOV does not change flags, so removing per-correction constant loads preserves arithmetic predicates. For OLD=0x00000000fffffffe=~MC modulo2^64, ADC(y,OLD,C) and SBC(y,MC,C) have identical result and NZCV when corresponding flag-setting forms are used; conversely SBC(y,OLD,C)=ADC(y,MC,C). Every original CONST_B consumer is transformed this way before MUL MC uses are introduced. The MC constant is a single encodable logical-immediate MOV, hoisted once per nonempty loop. The structural witness checks the exact removed-MOV counts and rejects any generated loop instruction writing either reserved constant register. Initializing them costs two instructions per nonempty loop; empty-loop paths skip this work.

The fused operation branches on the final raw-square ADCS carry before any flag-changing instruction. If raw carry is0, subtract c and keep on no borrow, otherwise add p. If raw carry is1, subtract c and keep on borrow, otherwise subtract p. The signed-high value would be rawcarry-borrow, so this is the same mathematical correction without materializing the fifth limb. The negative low path branches over the high path. Labels are unique per expansion. After the final raw-low ADCS, the four carry/borrow cases cost8,13,12,8 before common stores; no uniform dynamic saving is assumed.

Each direction-specialized forward entry costs at most441/4*n+46; the per-item zero test is removed, so defective rows complete with bounded garbage instead of aborting. Main-loop items chain the prefix through alternating register banks (no A-loads); the at most 7 pinned tail items keep memory loads and are covered by the fixed premium. Pair emission is fused into the backward kernel: each incidence emits its two records inline to block windows (16/17 instructions plus one flush-test branch), so the fused backward bound is2533/4*n+33 and there is no separate emitter pass; the v40 base 616*n+11 plus4 per carry-1 incidence is retained inside it. Every new generated eightfold loop keeps the same control skeleton in the same order; forward mains alternate even/odd bank-chained copies. A nonempty backward/emission call costs exactly4-14*floor(t/8) more than v40, a nonempty forward call exactly6-30*floor(t/8) more, while empty calls are unchanged. Because14*floor(t/8)>=7*t/4-49/4, calls use2533/4*n+33 fused-backward; forward uses441/4*n+46 (main items 110 worst, tail singles 114 worst, prologue/epilogue/tail overhead 19, so 46-19-3.75*r>=0.75 for every remainder r<8). Sixty-three fixed per segment (46 forward incl. tail premium and bank preload, 17 backward) are charged for all87 segments, 63*87 total. The row still pays every garbage/empty/small-call case: defective rows complete within the same envelopes (all kernels preserve below-R0 on zero inputs, the Fermat inverse has no division, all addressing is sequential, garbage candidates fail final verification, and the defective event keeps its already-charged d_batch exclusion); empty calls retain the prior fixed bound. No empirical branch average is used.

Full-flush publishes go through the shared flush_publish_b subroutine (15 taken / 12 dropped plus 3 call-site, unit-tested on 32 cases), charged at 18 per firing in the hop1 row (600 to 610 per 128 leaves). The arithmetic/emission variable price is now2974/4 per incidence plus4 per counted carry. One39479 inversion and old1024 fixed row allowance remain. The per-pair 63*87=5481 segment-fixed allowance is removed: proof-appendices.md
states that those terms and dispatch are covered by the additional2^17 per-row driver
charge, which the driver row retains in full as(1<<17)*hi per pair. Charging both pays
twice for the same work. Carry-meter semantics, failure handling, cold construction and all original phase arguments below are unchanged. Hop2 completion now uses DONE-flag signaling (no per-message acknowledgment or in-flight updates); see RECEIVER_PROTOCOL.md. These constants give exactly98153882402615006327661479627458466545664 instructions; range_ledger.py verifies the outward rational score126206384844364137/10^15.

The fused transform first splices inline window emission plus a shared flush subroutine into the four rolled backward bodies (frame spills, machine-checked liveness peak 30 with publish scratch free on all variants), then the loop transformer unrolls both forward entries and all four fused backward variants: it checks that the body never reads or writes count x19, makes eight uniquely labelled copies and uses a biased unsigned count. Initially subtract8 and skip to the tail on borrow. After each group subtract8 and repeat on no borrow. At exit add8 to recover a remainder0..7, then run the unchanged single-item loop. A zero initial count uses the unchanged early return. Internal arithmetic branches are renamed per copy; each item's pointer increments and output order are unchanged. The full finite replay compares every tested segment's entire memory and exact count against frozen v40, with a routing differential on flush/no-flush cases. This finite test supports, but does not replace, the control-flow induction.

The generated arithmetic listing has104036 encoded bytes, checked by parsing its scalar instructions. All supplied .s helpers are independently macro-expanded and counted by the same parser. Their combined encoded bytes plus104036 generated, 1860 unrolled-block and65536 stack bytes are asserted within the existing1MiB code_stacks allocation. Raw source length is NOT used as an encoded-size bound. The emitted scalar program is a public algorithm constant; per-node initialization can materialize each32-bit instruction with two immediate moves, one store and one pointer increment, costing at most the combined encoded-byte count in scalar instructions. Thus the existing1MiB code initialization charge covers delivery/materialization, with no runtime Python generator or uncharged expansion assumed, and no new lookup or dispatch.

## 4. Fully charged transport and grouping

Use the explicitly specified two-dimensional blocking ports of RECEIVER_PROTOCOL.md. Hop1 uses the explicitly fixed128-record window, native551-instruction delivery/ACK and transaction-bound completion protocol of BLOCK_PROTOCOL.md. Its full-flush work is charged at600 per block, and all partial/initialization/finishing work has a separate4096*SIDE per-node reserve. Ordinary final-cursor, barrier and resolution messages retain their addressed delivery cost. Hop2 uses a separately charged resident scalar receiver context at each node:16 instructions for accepted records,12 for full-bucket drops, with no per-message acknowledgment or in-flight updates. Each sender closes its phase stream with one DONE marker per destination; the receiver counts markers and signals per-phase completion with a single STLR to its DONE cell, which the main context awaits with the same LDAXR/WFE completion primitive as barrier epoch words. Both contexts' work and state are included. The hop2 receive window is released only after the final flag load of the fixed five-word packet, when all four payload words and the flag are already in receiver registers; the two earlier payload loads never release it, so a publication cannot overwrite a packet being read. No direct-home shortcut, uncharged key operation, software polling or unbounded queue is assumed. These are abstract instruction-model assumptions, not physical throughput claims.

Key digits stay d1=K0 mod2^35,d2=K1 mod2^35,d3=K2 mod2^26,d4=K3 mod2^48. The record is(K2,K3,w2,w3), w2=i|(k<<35), initially w3=d2|(sign<<35)|(l<<36). After extracting d2, h2s replaces only w3bits0..34 with source column j. Thus a home reads source(i,j) directly; sign/l and the absolute k are preserved. k<2^29,l<2^28, including maximum-boundary and post-emission-wrap cases.

Partition each node's X indices into two contiguous phases of at most qphi=ceil(q1/2) indices; do not renumber k. For each phase initialize fixed payload windows and local/remote cursors, emit ordered128-record blocks into sender-row slots at(d1,j), flush each remaining partial block with initialized padding, send clamped valid final cursors, drain and barrier. Padding is excluded by those cursors. Each slot has CAP1=25*2^18 records. Stream all its stored records to(d1,d2); no hop2 pair slots exist. The destination receiver computes d3, reserves32bytes in its tight full-record bucket, compares with its explicit end, rolls back on full, otherwise stores the four payload words and increments d4. After streaming, each sender publishes one DONE marker per destination; each main context awaits its receiver's DONE signal, then all nodes barrier before reusing hop1 storage.

Buckets and d4 counts persist across both phases. Only after both drain, prefix-sum d4 and scatter pointers to full records in d3 order. The full-record scatter costs6 per record plus2 per8-record block. Final d4 cursors are exclusive range ends. The unchanged paired scanner checks only K2-adjacent pairs inside each d4 range, at most(195/64)*records+26*ranges+8*(CAND_CAP+1)+2 instructions. Cap CAND_CAP=2^28 is checked at every attempted append.

Consequently matched keys fix 70 home bits, 48 d4 bits and all 64 K2 bits:
182 distinct bits. The d3 bits are already within K2. Partial matches are
not treated as true group collisions. Each is resolved and verified.

Resolution uses the four explicit bounded64-byte-record phases of RESOLUTION_PROTOCOL.md, not count/offset exchanges. Source and home tags use35-bit node coordinates; X index29bits, Y index28bits and sign1bit. Source tables keep original indices and row-specific signed b. The protocol includes initialization, caps, partial-reply flags, allocation, source gather and a final bounded fan-in output tree. Its complete work is charged by four N*2^41 rows and source service1024 per capped request. Source cap2^29 and relay allocation2^34 records are fixed; the original stricter logical relay loss bound is retained as explained there. A home reconstructs X+sY and X'+s'Y', tests equality up to sign, solves the nonzero scalar denominator and verifies [x]P=Q. Scalar arithmetic uses the complete COLD_ACCOUNTING.md bound42410496 per candidate.

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
whole group. Given all generator points, every b stays independent uniform:
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
the two node b values stay independent uniform. Such pairs are discarded.
For a good witness canonical key k outside Bbad, exceptional raw labels at least p cannot equal k; all records with complete label k therefore represent its two-point class. Under the separately charged no-intruder event, this class forms the relevant contiguous group. Only for such a good witness, if its key class has any useful pair, it has a useful adjacent pair: writing its
members u_i Z, uselessness is equality of n_i/u_i, an equivalence relation.
If every adjacent pair were useless, all members would have the same ratio.

## 6. Defects, capacity losses and success expression

Every actual raw key belongs to the fixed set U(k). FILTER_FREE.md supplies the set-valued ideal-record coupling and witness exceptional-band charge. A deterministic carry meter is also specified in CARRY_MEAN.md/CARRY_TAIL.md; it has an explicit conditional failure charge, never an assumed average cost. A batch defect probability is at most
d_batch=q2*2/(r-1)+1/r (zero denominators or identity X). Excluding any
defective batch in the whole run costs at most N*q1*d_batch. The event is
explicitly charged even though a defective batch's cost is no larger.

For the35-bit slot patterns and96-bit home/bucket pattern, KB(f)=(513/512)*2^(256-f)/(r-1) stays valid after adding the exceptional-label fibers. For home70 use coefficient2+2^-30; for group144 and match182 use coefficient2+2^-14. FILTER_FREE.md derives these exact mask-specific bounds from FOURIER_CAPACITY.md and top32-bit exceptional-label counting. They concern the fixed predicate U(k) intersects the cell, not the selected raw representative. Let rho=(r/(r-1))^4.

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
The losses are d_cl=h/(2^28-4), d_rq=g/(2^29-8), d_rl=v/(2^30-16).
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
negative exponential by the even 40-term Taylor polynomial, an upper
bound by the sign of the Taylor remainder. The resulting rational success
bound exceeds 0.3900000015160294 at q1=377468911. The declared success is conservatively 0.39. None of this identifies reviewer confidence with algorithm success.

## 7. Work and memory

range_ledger.py contains the complete exact rows, including both hops, every delivered record, empty tables/ranges, all failed caps, capped scalar sampling, preprocessing, phase4, resolution, verification and output. No unallocated per-leaf bookkeeping reserve is charged: all record operations are itemized below and all driver/context work is charged per batch or phase. Work is summed across both scalar contexts at every node.

Cold generation/header work is charged23441408instructions per requested generator item, and candidate resolution42410496instructions per capped candidate, including failed verification. These are the explicit itemizations in COLD_ACCOUNTING.md; the shared-per-node generator and exact scalar-inverse refinement follows the itemization in bdf01aa5-16e3-4619-a8eb-1aa0f6b6b29b (its acceptance is not a premise). Preprocessing stays below2^118units and is included.

Reserve the first2^52bytes of each node for code, both contexts/stacks, packet/ack state, DONE cell/granule, d4 counters, scalar/X/Y tables, batch storage, source/relay/candidate buffers and control tables. The fixed-region calculation in range_ledger.py includes the extra d3 end table and DONE granule and stays below that reservation. Place the phase receive region next,32*SIDE*CAP1bytes; then2^26 tight full-record buckets of32*CAP_3bytes; then an8*2^26*CAP_3 pointer array. All regions are nonoverlapping and eight-byte aligned. The exact endpoint is below2^64, and the attempted hop1 cursor bound cannot wrap. Bucket cursor temporary increments are at least the full pointer-array region below2^64; a full cursor is restored before the next packet is consumed. This is a simultaneous-allocation bound, not an assumed lifetime reuse saving. AcrossN=2^70nodes total storage is below2^134, hence the declared2^136bound. No nonuniform advice is used.

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
This reproduces the relevant dependency-graph inequality directly. It
does not assume Poisson-distributed collisions or independent leaves in a node.

## 9. Cold-row and finite implementation accounting

The ledger's generation/emission row is the exact per-batch field/kernel count
plus a conservative fixed driver allowance. The separate hop-one row adds
full-block flush, delivery and in-flight counts at600*S/128: ordinary per-incidence sender emission is already included. BLOCK_PROTOCOL.md proves the block bound and separately charges4096*N*SIDE for both phases of initialization, partial flush, finishing, endpoint registration and control.
Hop2 sender costs7 instructions per record plus2 per8-record block,11 per source slot and2 fixed per phase. The dedicated receiver costs16 accepted/12 dropped, with no per-message acknowledgment or in-flight updates; charge16 for every potential record. DONE-marker publication costs3+6*SIDE per sender per phase, DONE handling costs8 per non-final and10 per final marker, and the DONE signal/wait/reset costs are charged per node per phase; all sit in a separate fixed row. Full-record scatter costs6+2/8; the unchanged scanner costs195/64. These ordinary grouping terms total595/64, and hop2 totals93/4 per potential record. No old bucketing pass is charged or executed.

Two phases need two hop1 cursor initializations and two hop1 finishing sweeps, retaining6*SIDE+2 and34*SIDE+6 instructions per node. H2 source-slot traversal occurs twice. Initializing each bucket's cursor/end pair costs5, and pass4 adds11 per bucket. Zeroing/prefixing costs4+10 per four d4 counters. The combined grouping-table row includes all four fixed returns from init3, pass4, zero4 and prefix4. The scanner fixed/candidate terms are unchanged. Sixty-four barriers at1024instructions/node and another2^14phase/context-control allowance per node are charged. RECEIVER_PROTOCOL.md itemizes their scope, including resident startup and the final blocked receive. The block row includes a39-instruction margin above its561-instruction full-block body; otherwise source address/tag extraction, all cap tests, receiver release, scatter and scan work are already in these rows. Per-batch setup and phase/context control stay separately bounded.

The remaining rows explicitly charge all generator points and headers,
eight-attempt scalar sampling, program delivery and node initialization,
all capped resolution, source request service, four explicit bounded-routing stages,
barriers, target validation and output. Reply record data is fetched
from stored tables; no scalar multiplication is assumed free. Source serving
costs at most 1024 per request including addressed copies. Final resolution
has the42410496instruction itemized bound justified in Section7 and COLD_ACCOUNTING.md, even for false
partial matches and failed verification. All maxima are charged whether or
not their paths occur.

For final scalar arithmetic use the exact construction in COLD_ACCOUNTING.md: a256-step shift/add product costs36928 instructions including loop/setup; a512-product inversion with32 controller instructions per product and512fixed costs18924032. Five-word carry-aware addition/reduction preserves canonical residues. This is the operative bound, not a rounded product bound multiplied without controller cost. scalar_check.py and cold_accounting_check.py supply independent finite arithmetic checks.

All hot assembly uses base scalar GPR operations, plus the explicitly
specified LSE atomics for routing. Constants in MOV and logical masks are
single AArch64 encodings: zero/one/powers of two, contiguous bit masks of the displayed widths,
and the cyclic run of bits in 0xffffffff00000001. Branch displacements are
within the source blocks. Public parameters and the fixed inverse exponent
are program constants; delivering and initializing the program is charged.
The range emitter loads the original Y ID from each72-byte record; it never depends on a wrapping ordinal increment. Every emitted index and memory address fits its stated
word or region. The first-hop cursor bound includes all dropped attempts;
the bucket cursors undo on drops.

The new block protocol preserves every retained sender prefix and downstream order, so the exact success arithmetic is unchanged. The integer ledger sum, bound comparison and layout
assertions are in range_ledger.py. The score comparison uses the exact rational logarithm enclosure described below, not a rounded floating logarithm. The preprocessing ceiling is checked on a conservative superset
of setup and initialization rows; all these rows stay included in total
work. The printed logarithm is explanatory only.

## 10. Evidence, scope and provenance

Field/body counts, full batched arithmetic with inversion, metadata
emission and denominator/cap controls run through the supplied
AArch64 interpreter in the protected runner. Tests compare against
independent Python integers and affine point formulas.
The two-hop test executes deliveries and data-path kernels on a 4-by-4
machine against separate references; its planted keys isolate
transport only.
The phase4 test transports normalized actual P256 points and row-specific labels, recovers six planted targets, and tests dropped and partial replies. None of these
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

chain_witness.py and range_witness.py deterministically define the generated range program from the original call-based schedule and supplied new field/fused bodies. It specializes memory addresses before bijective register renaming. Stable pointers/counters are x19..x23, except that the forward loop no longer keeps the X point or previous-prefix pointers live: X is loaded once in the prologue into x7-x10 (which the renamed multiply never writes) and the prefix is read relative to x21 (pinned tail) or the ping-pong banks (mains), so x22/x23 are unused there and their setup MOVs are deleted. The first backward per-item inverse is produced directly in x24..x27 and stays live across prefix update and numerator arithmetic; both slope products read that bank in place of their right-operand loads. The y-minus operation leaves its original coordinate inputs intact; y-plus reuses those input registers and caches its result in x28,x29,x30,x18. The first slope/square/fused chain does not touch that bank. The second slope consumes it. Original LR is saved on entry and restored on return; no BL occurs inside these inlined loops. No caller-visible scratch or source point is removed.

Chained results substitute producer output registers for the next kernel's left-operand loads. Scratch32..63 stores and loads are dead after inverse caching and are absent. The forward nonnegative direction-specialized difference stays in x11..x14; no zero test remains; defective rows complete as garbage. Multiplication temporaries are renamed away from that live bank, and its two right operand loads become direct reads. The initial Y load is postindexed by72 and its upper limb load uses offset-56; the trailing pointer ADD is omitted. The advanced pointer always completes the iteration. Branch labels are uniquely renamed. Arithmetic branch predicates are preserved because address substitution/renaming adds no flag-setting instruction. The witness checks all structural preconditions and exact generated source equality. Random-register/address differential tests check every raw output for membership in U(k) against independent original canonical loops with separately inverted signed prefixes, plus source memory, dead-scratch canaries, SP and LR.

## 12. Complete ledger and finishing deliveries

The total is98153882402615006327661479627458466545664instructions, with log2(total/1000)=126.2063848802652355... . A64-term rational atanh enclosure and independent ln2 interval certify the outward claim126206384844364137/10^15. The two hop1 finishing sweeps include delivery and in-flight work; no second-hop finishing table is assumed. All cold rows, raw emitted records, failed caps, source service, sampling, grouping, verification and receiver work remain explicit. No prior validity judgment is a premise. This exact source requires its own external review. Fixed block windows, hop1 transaction-bound ACK events and hop2 DONE-flag completion are explicit machine assumptions; a flag-polling or general remote-memory implementation is not covered.

## Paired-comparison64-record range scan

For consecutive keys A,B,C, CMP A,B; CCMP B,C,#4,NE; B.EQ branches exactly when A=B or B=C. The slow path appends each equal adjacent pair separately with its own capacity check, in the same order as the original scanner. No relation or success event changes.

Three key registers rotate and four pointer-register pairs rotate. Initial copies put the range's first key and address in every possible suffix-start register. The pair remainder r=((n-1)//2) mod32 selects entry at pair32-r (r=0 selects the end test), using five TBZs and one added branch when r is odd. Each pair uses the previous key and record pointer, loads two new pointers and their keys, and retains the second key/address for the next pair. All potential entry registers initially hold the same first key/address. Since32 mod3=2 and32 mod4=0, the final key is x7 and final pointer x16, so MOV x5,x7 normalizes the next full block. The hot block costs32*(LDP+2LDR+CMP+CCMP+B.EQ)=192 instructions plus MOV,CMP,B.LO:195 per64 records.

Write n-1=64q+2r+epsilon with0<=r<=31 and epsilon in{0,1}. A nonmatching nonempty range costs25+6epsilon+(r mod2)+6r+195q including range control. Subtracting195n/64 gives1405/64+189epsilon/64-3r/32+(r mod2). Its maximum is1652/64=25.8125 at epsilon=1,r=1, bounded by26. Empty ranges cost5. The slow path adds7 for only the first equality,8 for only the second,11 for both; an odd-prefix hit adds4. Charge8 per listed candidate plus one added8 for a full-list append attempt. Every append has its own capacity check and preserves original adjacent-pair ordering, even when both equalities hold. Global guard and RET cost2.

Thus a node uses at most195n/64+26*2^48+8*(2^28+1)+2 scanner instructions. This rational inequality holds per node separately, so summing yields195S/64 without an uncharged ceiling. Production S is divisible by64. The remainder suffix consumes exactly2r records after the optional odd comparison; every full block consumes64, and the exclusive end test prevents out-of-range loads. Empty ranges read no record, and no comparison crosses a range boundary. For identical input pointer/key sequences, scanner candidate identities, order and capacity-abort behavior are unchanged. This revision can change those input sequences; its enlarged set-valued load bounds, not an unchanged candidate-list assumption, justify success and cost.

The implementation is transferred unchanged from our submitted secp package b00400fa-b4b3-440c-accb-657046f12efa; the record-pointer/end-list ABI is identical here. scanner_check.py exercises all remainder lengths, multiple blocks, empty ranges, both equality paths and caps, including4096-record cases. routing_check.py independently models every executed instruction and compares packet/candidate order across all16 nodes. Finite checks supplement the universal count and semantic argument above.


## 13. Grid, masks and provenance

SIDE=2^35,N=2^70,q2=2^28 and the35+29/35+1+28 metadata allocation are unchanged. The96bucket bits are[0,35),[64,99),[128,154); group144 and match182 include sixteen bits of the exceptional labels' top32constant bits. FILTER_FREE.md therefore enlarges the crude group/match coefficients rather than reusing a false uniform-bin assumption. FOURIER_CAPACITY.md stays the normal canonical-point lemma; the small extra fibers are bounded explicitly.

Immediate arithmetic/scanner source is our44137210-aa9d-47fa-8cd0-7daad206a9d1. The two itemized cold prices follow76ceba30-f1fa-4538-9c48-bc947cbd44bb. The removal of unallocated8*S padding follows37b8edde-addf-451f-9215-0dac6e7b6e79 and a fresh operation-by-operation audit of this changed receiver construction; public acceptance is not the justification. New work is filter-free set-valued-label accounting and the two-phase full-record receiver, with DONE-flag completion and fully revised cost/success/layout rows. The secp paired scanner and inherited constructions stay credited. Tools: direct implementation, read-only audits, bounded protected replay; no reviewer is sole author.

## 14. Complete current cost and storage binding

range_ledger.py is the exact integer ledger. Its2974/4-incidence row includes forward arithmetic and fused backward arithmetic with inline pair emission; the added4*carry_cap row pays every admitted high-carry correction. RANGE_DRIVER.md separately bounds construction and every row/segment operation. The retained normalization allowance is not reduced. Both hops, all grouping/scanning, failed caps, sampling, source service, four explicit resolution stages, cold verification and output stay separately charged. RESOLUTION_PROTOCOL.md assigns the final output reduction to unused fixed phase4 control space, not an uncharged all-to-root fan-in.

fixed_layout.py gives simultaneous nonoverlapping addresses for every source, label, tree/sort, batch, candidate/request/relay/home/flag/context/port/cursor region. Context pointers are initialized to those exact bases. d4 counters occupy[2^51,2^52). Bulk hop1 slots, full-record buckets and pointers follow2^52 and end at18214482104660525056<2^64. The sender header port includes its contiguous local cursor table at+64; block remote cursors, fixed payload windows, endpoint metadata and delivery context have separate allocations listed by fixed_layout.py. The batch has a64-byte prefix reservation before its first record. Output mailboxes start exactly after the2^20 code/stack region and have separate4096-byte child slots with isolated flag granules. No flag is reinitialized after publication begins. Peak across2^70 nodes is below2^134bytes, conservatively declared below2^136.

The complete proposed total is98153882402615006327661479627458466545664instructions, with log2(total/1000)=126.206384844364137... . The rational logarithm enclosure supports outward126206384844364137/10^15. The package requires external exact-version review and a separate user submission instruction. Existing submission approval does not transfer.



The remaining named appendices are in proof-appendices.md and are part of this same claim.


## Short Montgomery multiplication and exact cost change

The new scalar multiplication core is adapted from Amazon s2n-bignum commit56152cdb51f47ef9f63b179cda0a8b7518379240, arm/p256/bignum_montmul_p256_alt.S. Copyright and Apache-2.0 OR ISC OR MIT-0 licensing are preserved in p256_fields.s and the repository LICENSE. Primary source and its corresponding arm/proofs/bignum_montmul_p256_alt.ml are available at https://github.com/awslabs/s2n-bignum/tree/56152cdb51f47ef9f63b179cda0a8b7518379240 . This attribution is not a claim that the upstream HOL proof automatically certifies our patched tails or register transformations.

The core first forms the full four-by-four schoolbook product with MUL/UMULH and explicit carry propagation. Its low limbs are reduced four times using the cancellation digit t and p's limbs from section3. Writing MC=2^64-2^32+1, the seven-instruction reduction adds t shifted32 to the next limb, carries t>>32 onward, and adds the low/high halves of t*MC to the remaining two limbs with carry. Each step is exact division of the current integer plus t*p by2^64. Four steps give raw T=(ab+m*p)/R0 with0<=m<R0; uniqueness of the cancellation m modulo R0 shows this T equals the previous multiply's raw value for all a,b<R0. The final high carry is consumed immediately by B.CC. The weak tail keeps T if T<R0 and returns T-p otherwise, always in[0,R0). Old and new correction predicates agree pathwise.

Including RET and standalone tail constants, old MUL had112 instructions when correcting and106 otherwise; the new one has110 and104. The generated loops remove the same two tail materializations and hoist the one MC materialization; every executed MUL is exactly3 instructions shorter than the previous generated version, on either branch. A complete incidence executes one forward and four backward MULs, yielding15 saved instructions. Before unrolling, this proves forward122 and backward616+4*carry from the previous125 and628+4*carry bounds, not an average branch extrapolation. For a zero denominator, the full forward pass completes with bounded garbage within the 441/4-instruction incidence allowance, keeping its fixed46 allowance; the defective event keeps its already-charged d_batch exclusion, so no separate abort path is counted. No backward saving is claimed for unexecuted work. All attempted rows stay charged to the uniform full-row upper bound, so failed paths are not omitted.

The new RHS limbs occupy x7,x8,x13,x14 and are read-only after loads. MC occupies a distinct scratch register standalone and the reserved x2 in generated loops. The forward specialization maps the live canonical difference into that RHS bank and conflicting scratch into x24..x27. Backward slope products substitute cached inverse x24..x27 for all four RHS registers, without assigning MC to a cached limb. The original fused square and its raw-label case split are retained. A proposed alternate square offers no hot gain after MC hoisting and is not included. Canonical cold routines, fixed inverse, memory map, hop1/block transport, success probabilities, caps and all other ledger rows are unchanged; hop2 completion is re-derived above.

The inherited v40 generated listing is pinned to SHA25607f98f57a417d080a5ccdf49fc8d1381b501b174251f89606e5ea69a48b95a29 (v40 base with the zero test deleted in chain_witness.generate, X hoisted to the forward prologue, and all per-item pointer updates folded into post-indexed memory ops). Its short-multiply improvement supplies the reference variable coefficient122+616+27=765. This version additionally applies the exact eightfold transformation in section3, giving2974/4 per incidence, plus the separately charged4 per carry and5481 fixed per row. The new generated listing SHA256 is32fa3f312d18679d2cdb3b018749588d0107f45681b411465dc894c4889f2cd3 and its encoded size is104036 bytes. b5-b11 built the eightfold-unrolled fused-emission kernels with garbage completion, X-hoist, post-indexed pointers and emission pair-loads (449/4*n+30 forward, 2533/4*n+33 fused backward at 34*87 fixed; rolled base untouched, pin holds). b12 (now 441/4*n+46 forward, 2533/4*n+33 fused backward, 63*87 fixed) chains the forward prefix through alternating banks (A=x28,x29,x17,x18; B=x15,x16,x24,x25) on main items, deleting both A-loads (-2); even copies read A/write B, odd copies swap; the tail stays pinned (validated b11 bytes) and the prologue preloads bank A (+2 fixed). Banks dodge X/d/consts/LR by construction (mapping injectivity asserted); A-stores stay for bwd. b14 replaces hop2 per-message acknowledgment plus sender in-flight updates with DONE-flag completion (address-word flag test +2, per-message STLR -1, one DONE marker per sender, single per-phase STLR awaited via the barrier LDAXR/WFE primitive). Hop2 falls97/4 to93/4 per potential record; DONE publication/handling/wait are a separate fixed row. All scalar sampling, cold work, failed work, remaining transport and output remain charged.


Control-flow supplement. The unchanged v40 forward and backward bodies do not inspect x19. Execute8 copies with unique local branch labels, decrementing a biased count8 at group boundaries; tail restores count and runs the original loop for0..7 remaining items. No fail label remains; defective rows complete as bounded garbage within the same envelopes. Fused bodies keep the same loop tail, so unrolling them obeys the same law with the fused per-item base. Every item keeps the same order and pointer updates. Every nonempty backward/emission call changes exact dynamic count by4-14*floor(n/8), every nonempty forward call by6-30*floor(n/8), independently of arithmetic branches. Empty calls unchanged. Asymptotic saving is15/4 per forward incidence and7/4 per backward incidence; small forward calls can lose6, small backward4. The proof funds63 fixed instructions for each of87 possible descriptors, keeps all prior driver/cold charges, and checks encoded program size. This supplement alone is not the full claim; see proof.md and range_ledger.py. V40 stays immutable and its review does not authorize this version. Full finite replay retains source/label/route/recovery checks and exact differential against v40 backward memory for every tested segment.

Interpreter compaction removes only interpreter comments and function docstrings, verified by inert AST comparison. Module provenance retained; full new replay required. V42/V43 exceeded the archive cap and were never submitted.
