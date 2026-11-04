The named .md documents cited below are included verbatim as named appendices in proof.md and proof-appendices.md, not separate missing files.

# P256 range-specialized stored-leaf search

Proposed complete analytical claim for exact-version external review; not certified.

NORMALIZATION.md, FILTER_FREE.md, RECEIVER_PROTOCOL.md, FOURIER_CAPACITY.md, COLD_ACCOUNTING.md, RANGE_DRIVER.md, CARRY_MEAN.md, CARRY_TAIL.md, RESOLUTION_PROTOCOL.md, BLOCK_PROTOCOL.md, section15 (grouped early abort, hash chains, expected-work ledger) and fixed_layout.py are operative parts of this argument.

## 1. Target, representation and candidate claim

The target is exactly the protected P-256 group, including a=-3, its specified
prime p, prime order r, cofactor one, generator P and uniform Q=[x]P. Output is
the canonical scalar x. Q=infinity returns zero. No auxiliary powers or
secp256k1 endomorphism are used. Let R=2^256 mod p. Store canonical generator Montgomery
coordinates; hot intermediate values may be redundant in[0,2^256), and use K(Z)=R*x(Z) mod p for a nonidentity point Z. This is an
exact key for the two-element class {Z,-Z}: the curve equation gives at most
two points for each x, and no nonidentity point has y=0 because r is odd.
Multiplication by R permutes the field. Actual raw labels may differ on the two explicitly charged exceptional bands; FILTER_FREE.md proves the point-only label-set domination used for all capacity and pair-load bounds.

The claim is expected classical work below 2^125.937710102075634 units at the
protected C=1000 scalar instructions/unit (expectation over the challenge and
algorithm coins of section15's grouped early-abort run; every run costs at most
the disclosed 2^178.3505 units, with all512 groups run and every uncapped hash
walk visiting its chain), success at least 0.39, peak memory below2^136 bytes,
preprocessing below2^119 units included, and zero nonuniform advice.

## 2. Random construction and pair arithmetic

There are N=2^70 nodes, indexed (i,j) with each coordinate in [0,2^35).
At each node draw b and q1=377471100 independent a_k uniformly modulo r,
and q2=2^28 independent c_l uniformly in {1,...,r-1}. Coins and draws at
different nodes are independent. Each scalar draw tries at most eight fresh
256-bit values, rejecting out-of-range values and zero where required. A cap
failure aborts the run. The rejection probability is below 2^-31, so the
union cap loss is at most N(q1+q2+1)*2^-248. Sampling is charged 256 scalar
instructions per requested scalar. Under the protected ideal-coin convention,
RNDR supplies the bits; a hardware-failure assumption is not silently added.

Initially construct R_g=[b]Q, X_k=R_g+[a_k]P and Y_l=[c_l]P. These original points and indices are used in the probability argument below. Apply the exact signed normalization and paired scalar updates in NORMALIZATION.md: store X'_k=epsilon_X X_k with a'_k and row-specific b'_k, and Y'_l=epsilon_Y Y_l with c'_l. Build the y-partition/x-sorted range tree of72-byte records with original Y IDs; retain canonical source tables at original indices. Keep original X index k and group membership (section15). Preserve original b while constructing X, separately from stored row-specific b'.

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

Direction specialization is machine-checked limb by limb: pair_fwd_hi computes Y-X and is selected exactly when Y.x>=X.x; pair_fwd_lo swaps all four limbs (chain_witness.py normalize_transform) to compute X-Y and is selected exactly when Y.x<X.x. range_witness.assert_limb_specialization asserts the exact limb sets of both forward sections and of the four backward variants on the generated program, and the probe executes pair_fwd_lo on ten boundary operand pairs (X.x=2^128 with Y.x=1 and its mirror among them) and asserts the stored denominator equals the exact absolute difference.

The fixed inverse computes exponent p-2 with255 squares and127 multiplies. Induction on the preceding REDC range argument proves all intermediates stay belowR0 and congruent to the specified modular powers. The fixed inverse instruction sequence is retained byte-for-byte after macro expansion; it calls external k_mul/k_sqr, and its conservative39479 allowance is retained without crediting their smaller counts. Canonical cold point operations instead use the separately specified namespaced canonical routines. On non-defective rows its input is a nonzero prefix product and no zero-residue representative0 orp arises; defective rows complete as bounded garbage within the same envelopes, detailed below in this section. Fermat's theorem therefore gives the field inverse. Each changed M/S path is no longer than its original canonical path, so39479 stays a uniform upper bound including copies, prologue, calls and return; no average branch discount is taken.

For final fused output write R0=2^256, delta=R0-p. Raw square T<R0+p and z=T-c lies in(-R0,R0+p), so h=floor(z/R0) is-1,0,1. The four-limb stored value is u=(z-h*p) modR0. For true canonical key k=z modp, interval splitting proves u belongs to U(k)={k}, with k+p additionally possible when k<delta and k+delta additionally possible when k>p-delta. These bands are disjoint and |U(k)|<=2. Every noncanonical label lies in[p,R0), including the noncongruent negative-wrap branch. FILTER_FREE.md gives the full case split and the revised load argument.

There is no emitter representation filter. A useful collision whose canonical key is in either exceptional band is charged by the retained unconditional Gamma*2^227/r^2 witness-loss allowance. Other exceptional labels are neither discarded nor assumed independent: fixed point-only set predicates dominate their bin, intruder and pair loads. There is no global S*2^-32 omission.

The generated loops reserve x0 and x2 for the zero-extended constants0xffffffff andMC=0xffffffff00000001, initialized once after original count/pointers are saved. After pointer specialization, the active field-data register bank is x1,x3..x17; it excludes both constants and the separate cached-sum bank. The y-plus carry temporary moves from x0 to x1. MOV does not change flags, so removing per-correction constant loads preserves arithmetic predicates. For OLD=0x00000000fffffffe=~MC modulo2^64, ADC(y,OLD,C) and SBC(y,MC,C) have identical result and NZCV when corresponding flag-setting forms are used; conversely SBC(y,OLD,C)=ADC(y,MC,C). Every original CONST_B consumer is transformed this way before MUL MC uses are introduced. The MC constant is a single encodable logical-immediate MOV, hoisted once per nonempty loop. The structural witness checks the exact removed-MOV counts and rejects any generated loop instruction writing either reserved constant register. Initializing them costs two instructions per nonempty loop; empty-loop paths skip this work.

The fused operation branches on the final raw-square ADCS carry before any flag-changing instruction. If raw carry is0, subtract c and keep on no borrow, otherwise add p. If raw carry is1, subtract c and keep on borrow, otherwise subtract p. The signed-high value would be rawcarry-borrow, so this is the same mathematical correction without materializing the fifth limb. The negative low path branches over the high path. Labels are unique per expansion. After the final raw-low ADCS, the four carry/borrow cases cost8,13,12,8 before common stores; no uniform dynamic saving is assumed.

Each direction-specialized forward entry costs at most441/4*n+46; the per-item zero test is removed, so defective rows complete with bounded garbage instead of aborting. Main-loop items chain the prefix through alternating register banks (no A-loads); the at most 7 pinned tail items keep memory loads and are covered by the fixed premium. Pair emission is fused into the backward kernel: each incidence emits its two records inline to block windows (9 instructions for the first record and 7 for the second, including the one-instruction TBZ full-window test, in place of the two output stores), so the fused backward bound is2489/4*n+33 and there is no separate emitter pass; the v40 base 616*n+11 plus4 per carry-1 incidence is retained inside it. Every new generated eightfold loop keeps the same control skeleton in the same order; forward mains alternate even/odd bank-chained copies. A nonempty backward/emission call costs exactly4-14*floor(t/8)-2t more than v40, a nonempty forward call exactly6-30*floor(t/8) more, while empty calls are unchanged. Because14*floor(t/8)>=7*t/4-49/4, calls use2489/4*n+33 fused-backward; forward uses441/4*n+46 (main items 110 worst, tail singles 114 worst, prologue/epilogue/tail overhead 19, so 46-19-3.75*r>=0.75 for every remainder r<8). Sixty-three fixed instructions per segment (46 forward incl. tail premium and bank preload, 17 backward) over at most87 segments are covered by the2^17 per-row driver charge below; no separate63*87 term is charged. Garbage/empty/small-call cases stay inside these envelopes (bounded garbage completion, sequential addressing, the d_batch exclusion of section2); empty calls retain the prior fixed bound. No empirical branch average is used.

Full-flush publishes go through the shared flush_publish_b subroutine (16 taken / 13 dropped; plus 3 call-site), charged at 19 per firing in the hop1 row (584 per 128 leaves: 561 delivery with ACK and in-flight bracketing, the 19-instruction fused flush, 4 reserve; the reference emitter's inline flush costs 8). The arithmetic/emission variable price is now2930/4 per incidence plus4 per counted carry. b16 renames the first backward multiply's accumulator limbs to x1,x17 so its low operand x3,x4 survives; the second multiply drops the `ldp x3,x4` reload, one instruction per backward incidence. b18 extends the reuse to the high pair: the first multiply's top accumulator limb goes to x28 (x6 is never written there and x28 is otherwise free, both asserted), so the second multiply also drops its `ldp x5,x6` reload, one more instruction per backward incidence. b17 shortens each inline record emission by three instructions: the full-window test `tst x14,#0xfe0; b.ne` becomes `tbz x15,#12` (payload windows sit at an8192-byte stride from an8192-byte-aligned base, so bit12 of a slot address is set exactly on the128th slot, where bits5..11 are zero); the d2/sign/ID word is formed by `bfi x1,x8,#0,#35` into the loaded ID word (whose low36 bits are zero) instead of `and`+`orr`; and the cursor-table address is formed by `bfi x13,x7,#3,#35` into the2^38-byte-aligned cursor-table base instead of `and`+`add`. b19 carries the ID word (x1) and the LDADD addend (x17) from the first emission to the second, and b21 carries the cursor-table base (x13) and the W2 word (x14) the same way (all asserted), so the second emission drops its `ldp`, `ldr` and `mov`, three instructions per backward incidence (2489/4 per incidence). Every stored word and address is unchanged. All three forms assemble under clang as single AArch64 instructions; no shifted carry form is used. The probe checks `ct = foct + 4 - 14*floor(n/8) - 2n` against the rolled fused base and identical final memory, and a dedicated full-window test prefills every emitted window to 127 records so all 280 fused emissions across its two trials fire the real subroutine (140 taken with exact publish metadata and cursor advance, 140 dropped with cursor reset and no publish), covering all four variants with per-segment counts inside2489/4*n+33+19F+16D. One39479 inversion and old1024 fixed row allowance remain. The per-pair 63*87=5481 segment-fixed allowance is not charged separately: the2^17 per-row driver charge, retained in full, covers those terms and dispatch. Carry-meter semantics, failure handling, cold construction and all original phase arguments below are unchanged. Hop2 completion now uses DONE-flag signaling (no per-message acknowledgment or in-flight updates); see RECEIVER_PROTOCOL.md. These constants enter the per-group rows of the expected-work ledger in section15; range_ledger.py verifies the outward rational score125937710102075634/10^15.

The fused transform first splices inline window emission plus a shared flush subroutine into the four rolled backward bodies (frame spills, machine-checked liveness peak30 with publish scratch free on all variants), then the loop transformer unrolls both forward entries and all four fused backward variants: it checks that the body never reads or writes count x19, makes eight uniquely labelled copies and uses a biased unsigned count (subtract8, skip to the tail on borrow, repeat on no borrow, add8 back for the remainder; zero uses the early return). Internal branches are renamed per copy; pointer increments and output order are unchanged. The finite replay compares every tested segment's memory and exact count against frozen v40; it supports, but does not replace, the control-flow induction.

The generated arithmetic listing has102456 encoded bytes, checked by parsing its scalar instructions; all supplied .s helpers are macro-expanded and counted by the same parser, and their encoded bytes plus102456 generated,1860 unrolled-block and65536 stack bytes are asserted within the1MiB code_stacks allocation (raw source length is not used). The program is a public constant; per-node initialization materializes each32-bit instruction with two immediate moves, one store and one pointer increment, within the1MiB code initialization charge, with no runtime generator, lookup or dispatch.

## 4. Fully charged transport and grouping

Use the explicitly specified two-dimensional blocking ports of RECEIVER_PROTOCOL.md. Hop1 uses the explicitly fixed128-record window, native551-instruction delivery/ACK and transaction-bound completion protocol of BLOCK_PROTOCOL.md. Its full-flush work is charged at584 per block (561 delivery with ACK and in-flight bracketing, the 19-instruction fused flush, 4 reserve); all partial/initialization/finishing work is charged per group in the worst-case group_fixed row. Ordinary final-cursor, barrier and resolution messages retain their addressed delivery cost. Hop2 uses a separately charged resident scalar receiver context at each node: it appends each accepted record to the node's single home record array and immediately joins it into the hash chains (15 instructions per accepted record plus6 per chain visit and4 per K2 match;7 for a full-array drop), with no per-message acknowledgment or in-flight updates. Each sender closes its group stream with one DONE marker per destination; the receiver counts markers and signals per-group completion with a single STLR to its DONE cell, which the main context awaits with the same LDAXR/WFE completion primitive as barrier epoch words. Both contexts' work and state are included. The hop2 receive window is released only after the final flag load of the fixed five-word packet, when all four payload words and the flag are already in receiver registers; the two earlier payload loads never release it, so a publication cannot overwrite a packet being read. No direct-home shortcut, uncharged key operation, software polling or unbounded queue is assumed.

Key digits stay d1=K0 mod2^35,d2=K1 mod2^35,d4=K3 mod2^48; no d3 bucket digit is used. The record is(K2,K3,w2,w3), w2=i|(k<<35), initially w3=d2|(sign<<35)|(l<<36). After extracting d2, h2s replaces only w3bits0..34 with source column j. Thus a home reads source(i,j) directly; sign/l and the absolute k are preserved. k<2^29,l<2^28, including maximum-boundary and post-emission-wrap cases.

Partition each node's X indices into G=512 contiguous groups of at most qg=ceil(q1/512) indices (sizes differ by at most one; range_ledger.py group_sizes); do not renumber k. Every node runs the same group at the same time, and the group loop tests the latched stop flag before each group (section15). For each group initialize fixed payload windows and local/remote cursors, emit ordered128-record blocks into sender-row slots at(d1,j), flush each remaining partial block with initialized padding, send clamped valid final cursors, drain and barrier. Padding is excluded by those cursors. Each per-group slot has CAP_PG=5*2^17 records; the slot region is reused by the next group only after the barrier that follows this group's hop2 drain. Stream all its stored records to(d1,d2); no hop2 pair slots exist. The destination receiver compares its register-resident record cursor with the last admissible position and drops the record if the array is full; otherwise it stores the four payload words at the next32-byte position, prepends the record to the hash chain of its K2 index and walks that chain (section12). After streaming, each sender publishes one DONE marker per destination; each main context awaits its receiver's DONE signal, then all nodes barrier. The group's new candidates (the list entries between the previous and the newly published list end) are then resolved by the four phase-4 stages of RESOLUTION_PROTOCOL.md; a resolver-confirmed x latches the stop flag.

The home record array and its cursor, the hash heads, next words and the candidate list persist across all groups. No d3 buckets, d4 counters, prefix sums, scatter pass, separate join pass or adjacency scan exist. Cap CAND_CAP=2^24 is checked at every attempted append.

A listed pair agrees on the 70 home bits, all 64 K2 bits (the chain test) and the 48 d4 bits
(the verified compare): 182 distinct bits. The hash-index bits are within K2. Partial matches are
not treated as true group collisions. Each is resolved and verified.

Resolution uses the four bounded64-byte-record phases of RESOLUTION_PROTOCOL.md, not count/offset exchanges. Tags use35-bit node coordinates; X index29bits, Y index28bits and sign1bit. Source tables keep original indices and signed b. The protocol includes initialization, caps, partial-reply flags, allocation, source gather and a final fan-in tree. Its work is charged by four N*2^41 rows per group and service1024 per capped request; the cumulative list has its2^24 cap, so at most2^24 candidates and2^25 requests are resolved. Source cap2^29 and relay allocation2^34 are fixed; the stricter logical relay loss bound is retained. A home reconstructs X+sY and X'+s'Y', tests equality up to sign, solves the nonzero denominator and verifies [x]P=Q. Scalar arithmetic uses the COLD_ACCOUNTING.md bound20759040 per candidate.

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

For the35-bit slot patterns, KB(f)=(513/512)*2^(256-f)/(r-1) stays valid after adding the exceptional-label fibers, and the same lemma gives KBH=(513/512)*2^186/(r-1) for the70-bit home pattern of the home record array. The pair-load term h keeps the crude home70 coefficient2+2^-30 (KB70); for group144 and match182 use coefficient2+2^-14. FILTER_FREE.md derives these exact mask-specific bounds from FOURIER_CAPACITY.md and top32-bit exceptional-label counting. They concern the fixed predicate U(k) intersects the cell, not the selected raw representative. Let rho=(r/(r-1))^4.

Condition on the witness's four generators, but never on global cleanliness or actual raw-label choices. Retain the reservation12 fixed leaves. Its hop1 group slot has qg=ceil(q1/G), mean m_p=2*qg*q2*KB35+1 and variance bound v_p=m_p+1+60*(qg+q2)*KB35. Its home record array spans all groups, mean m_h=S*KBH+1 and variance v_h=m_h+1+60*(q1+q2)*KBH. Then d_p=v_p/(CAP_PG-12-m_p)^2 and d_h=v_h/(CAP_H-12-m_h)^2, CAP_H=ceil(m_h)+2^50. The witness's group is predetermined by its X index k, so the uniform bound over the equal-sized groups applies and no group multiplier is needed. There is only one pair-slot loss, because no second-hop slot exists.

Here is the variance justification. Leaves sharing a free Y but different
free X are independent. With a common free X and different nonzero Y, the
only missing value is X, so covariance is at most 1/(4(r-1)^2); if one Y is
fixed, it is at most KB(f)/(r-1). Two leaves of a free incidence are an
invertible linear image of (X,Y), Y!=0, and their covariance is bounded by
the source's cell-size term divided by r^2(r-1). The aggregate of these tiny
terms is below one at the stated S,q1,q2,r, even using a global S upper bound.
An incidence involving a witness generator and one free generator can have
fully correlated leaf indicators; counting all such ordered partners with
the generous60(qg+q2) factor for its group slot, or60(q1+q2) for the home array gives the remaining variance term. Chebyshev
then yields the displayed conditional drop probabilities. The exact script asserts CAP_P-12>m_p and CAP_H-12>m_h before division; a negative gap cannot be squared into a valid bound. No independent
hashing distribution is assumed.

For pair loads define M=2*KB(182), since a fixed opposite point offers at most two labels. The source's shared-generator conditioning cases apply to these fixed set-valued relations, with
 h=rho*(S^2*KB(70)*M+(2S+20)*M),
 g=2*rho*n0*S*M,
 v=18*rho*((S/2^35)*S*KB(35)*M+(2S+20)*M).
The losses are d_cl=h/(2^24-4), d_rq=g/(2^29-8), d_rl=v/(2^30-16).
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

The intruder rate is rho*(S+12)*KB(182): a hash walk stops only at a prior record
agreeing with the new record on all 182 bits (home, K2 and d4), so only a 182-bit-agreeing
record of another class arriving at the home between the two witness records can split
their class's chain (section15); every other record is merely visited. If no such intruder
arrives and the capacity events do not occur, the useful chain edge is listed and resolved. False partial-key matches can exhaust caps; their
work and losses are exactly what the pair-load bounds cover.

For the added meter set inc=q1*q2, mean=inc*(1/2+2^-116), cap=ceil(mean)+2^50, B=4*(q1+q2), var=inc*(q1+q2-1)/4, and d_m=var/(cap-mean-B)^2. The positive gap is checked exactly. CARRY_MEAN.md proves the original-index conditional mean/independence and CARRY_TAIL.md proves this tail. Sorting and normalization preserve the total; all capped-sampling failures are separately charged. Two witness-node failures give the displayed2*d_m term.

The explicit lower-bound expression is
1-exp(-mu)-CS
-mu*(1+1/r)*(1/r+rho*(S+12)*KB(182)
             +2*(d_batch+d_p+d_h+d_cl+d_rq+d_rl+d_m))
-N*q1*d_batch-N*(q1+q2+1)*2^-248-Gamma*2^227/r^2.
Here Gamma is the cross-node unordered-pair count. For a useful true collision, the endpoints share canonical k. The two bad bands have at most2delta keys and4delta curvepoints, each with at most two negation partners. Thus a pair is both a collision and in the charged bad set with probability at most8delta/r^2<2^227/r^2. Union-charge this before any capacity or survival conditioning. Actual records and partial matches are dominated by the fixed U(k) all-unordered-pair envelopes. Removing records can create adjacent matches, so no adjacency-list subset claim is used. Wrong labels can produce false candidates, all fully charged and verified.
The factor 1+1/r bounds (2/r)/P1. All quantities except exp are exact rational
functions. Round mu down to 200 binary fractional bits and upper-bound its
negative exponential by the even 40-term Taylor polynomial, an upper
bound by the sign of the Taylor remainder. The resulting rational success
bound exceeds 0.3900000016983508 at q1=377464334 (the least q1 with the bound at least0.390000001). The declared success is conservatively 0.39. None of this identifies reviewer confidence with algorithm success.

## 7. Work and memory

range_ledger.py contains the complete exact rows: the E-averaged per-group rows of section15 (arithmetic, emission, per-row driver, both hops including the receiver's append and hash insert of every delivered record, and the expected main-pair walks), the global same-row walk, K2-verification and defective-walk rows, and the worst-case rows for all failed caps, generation, capped scalar sampling, normalization/sorting, the carry meter cap, preprocessing, per-group fixed work, phase4, resolution, verification and output. No unallocated per-leaf bookkeeping reserve is charged: all record operations are itemized below and all driver/context work is charged per batch or group. Work is summed across both scalar contexts at every node.

Cold generation/header work is charged1789952instructions per requested generator item, and candidate resolution20759040instructions per capped candidate, including failed verification. These are the explicit itemizations in COLD_ACCOUNTING.md; the shared-per-node generator and exact scalar-inverse refinement follows the itemization in bdf01aa5-16e3-4619-a8eb-1aa0f6b6b29b (its acceptance is not a premise). Preprocessing (node initialization, generation, scalar sampling, the one-time hash-head zeroing and the range tree) stays below2^119units and is included; range_ledger.py asserts that bound on a superset of the setup rows.

Reserve the first2^52bytes of each node for code, both contexts/stacks, packet/ack state, DONE cell/granule, the receiver's list-end cell, the unused former d4-counter range, scalar/X/Y tables, batch storage, source/relay/candidate buffers and control tables. The fixed-region calculation in range_ledger.py includes the DONE granule and stays below that reservation. Place the per-group hop1 slot region next,32*SIDE*CAP_PGbytes; then the2^60 hash head words (2^63bytes); then the home record array of32*CAP_H bytes; then the8*CAP_H-byte next-array (one8-byte next word per record position, section12). All regions are nonoverlapping and eight-byte aligned. The exact endpoint18115304426367778224 is below2^64, and the attempted hop1 cursor bound cannot wrap. The receiver tests its cursor against the last admissible position before any store, so no temporary overshoot exists. This is a simultaneous-allocation bound, not an assumed lifetime reuse saving. AcrossN=2^70nodes total storage is below2^134, hence the declared2^136bound. No nonuniform advice is used.

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
full-block flush, delivery and in-flight counts at584*S/128 in total: ordinary per-incidence sender emission is already included. BLOCK_PROTOCOL.md proves the block bound; per-group initialization, partial flush, finishing, endpoint registration and control (1472*SIDE+6) sit in the worst-case group_fixed row.
The dedicated receiver costs15 per accepted record, including the append, the hash insert and the chain-head test (the window-release contract of RECEIVER_PROTOCOL.md is unchanged), and7 per full-array drop, with no per-message acknowledgment or in-flight updates; charge15 for every potential record. Its chain walks cost6 per visit and4 more per K2 match, charged by pair counting in section15; a listing replaces the three-instruction advance by a three-instruction append, so it adds nothing. DONE-marker publication costs3+6*SIDE per sender per group, DONE handling costs8 per non-final and11 per final marker (the final one also publishes the list end), and the DONE signal/wait/reset and the main context's list-end load are charged per node per group inside the group_fixed row; hop2 totals89/4 per potential record. No bucketing pass, separate join pass, d4 counting, scatter or scan is charged or executed.

Every group needs its own hop1 cursor initialization (3*SIDE+1), finishing sweep (17*SIDE+3), block work (1472*SIDE+6, BLOCK_PROTOCOL.md), h2s slot traversal (11*SIDE+2), h2done (3+6*SIDE), DONE handling (8*SIDE+3), DONE wait/reset (5), list-end load (1) and stop test (5): at most1517*SIDE+29, charged as1536*SIDE+64 per group per node for all512 groups in the worst-case group_fixed row. The receiver's resident registers are set once inside the2^20 phase/context-control allowance; the2^60 hash head words (64-byte aligned) are zeroed once by zero4096 at1026 per4096 words plus its RET: 512 DC ZVA zeroes with explicit addresses plus the loop. Cache lines are at least 64 bytes (specified geometry, true of every ARMv8-A implementation), so each ZVA zeroes at least its aligned 64-byte block and counts as one instruction like the other single-instruction memory operations. 2048 barriers at1024instructions/node and a2^20phase/context-control allowance per node (33 per kernel call, at most20 calls per group, the stop test and group loop) are charged. RECEIVER_PROTOCOL.md itemizes their scope, including resident startup and the final blocked receive. The block row charges584 per full block: the561-instruction full-block body with ACK and in-flight bracketing, the 19-instruction fused flush, and a4-instruction margin; source address/tag extraction, all cap tests, receiver release and the hash join are in the rows above. Per-batch setup and group/context control stay separately bounded.

The remaining rows charge all generator points and headers, eight-attempt scalar sampling, program delivery and node initialization, all capped resolution (20759040 per candidate, COLD_ACCOUNTING.md, even for false partial matches and failed verification), source request service (at most1024 per request including addressed copies), the four bounded-routing stages per group, barriers, target validation and output. All maxima are charged whether or not their paths occur.

For final scalar arithmetic use the exact construction in COLD_ACCOUNTING.md: a256-step shift/add product costs36928 instructions including loop/setup; a512-product inversion with32 controller instructions per product and512fixed costs18924032. Five-word carry-aware addition/reduction preserves canonical residues. This is the operative bound, not a rounded product bound multiplied without controller cost. scalar_check.py and cold_accounting_check.py supply independent finite arithmetic checks.

All hot assembly uses base scalar GPR operations, plus the explicitly
specified LSE atomics for routing. Constants in MOV and logical masks are
single AArch64 encodings: zero/one/powers of two, contiguous bit masks of the displayed widths,
and the cyclic run of bits in 0xffffffff00000001. Branch displacements are
within the source blocks. Public parameters and the fixed inverse exponent
are program constants; delivering and initializing the program is charged.
The range emitter loads the original Y ID from each72-byte record; it never depends on a wrapping ordinal increment. Every emitted index and memory address fits its stated
word or region. The first-hop cursor bound includes all dropped attempts;
the home-array cursor is tested against its end before any store.

The new block protocol preserves every retained sender prefix and downstream order, so the exact success arithmetic is unchanged. The integer ledger sum, bound comparison and layout
assertions are in range_ledger.py. The score comparison uses the exact rational logarithm enclosure described below, not a rounded floating logarithm. The preprocessing ceiling is checked on a conservative superset
of setup and initialization rows; all these rows stay included in total
work. The printed logarithm is explanatory only.

## 10. Evidence, scope and provenance

Field/body counts, full batched arithmetic with inversion, metadata
emission and denominator/cap controls run through the supplied
AArch64 interpreter in the protected runner. Tests compare against
independent Python integers and affine point formulas.
The two-hop test executes deliveries, data-path kernels and the receiver's hash join over two executed groups of three planned (a toy stop rule latches after the second) on a 4-by-4 machine against separate references, with the exact receiver path model and full path coverage; its planted keys isolate transport and grouping only.
None of these finite tests proves a production success rate or physical feasibility.

Material source e0d1ed74-946e-4f8d-a5ad-a6d219ad3697 supplies the general P-256 field kernel header. Material sources 80e9e8af-bfb6-40e3-8a0c-45068093d5ff and
1efcef88-0fe1-4f72-a69e-c18b5eceedbb supply our enlarged-cap accounting. The original paired stored-leaf construction
is source 3828873f-56f6-4cb8-be18-0817bc41a5bd; its full six-leaf transport
and dependency-analysis successor 472226d8-9930-469f-87f3-734cee673b83 is
also materially used. These source IDs are provenance, not premises; the required argument is reproduced here.
This proof replaces its curve arithmetic, six-class quotient, approximate canonicalization and defect terms, and recalculates all target-dependent quantities. The block/header concept is informed by public secp submission e59be87d-4235-4f77-ae3f-54a6ea896975, with independently counted P256 adaptation. No previous judgment is reused; the program and proofs here are the evidence for this claim.

## 11. Register and address transformation

chain_witness.py and range_witness.py deterministically define the generated range program from the original call-based schedule and supplied new field/fused bodies. It specializes memory addresses before bijective register renaming. Stable pointers/counters are x19..x23, except that the forward loop no longer keeps the X point or previous-prefix pointers live: X is loaded once in the prologue into x7-x10 (which the renamed multiply never writes) and the prefix is read relative to x21 (pinned tail) or the ping-pong banks (mains), so x22/x23 are unused there and their setup MOVs are deleted. The first backward per-item inverse is produced directly in x24..x27 and stays live across prefix update and numerator arithmetic; both slope products read that bank in place of their right-operand loads. The y-minus operation leaves its original coordinate inputs intact; y-plus reuses those input registers and caches its result in x28,x29,x30,x18. The first slope/square/fused chain does not touch that bank. The second slope consumes it. Original LR is saved on entry and restored on return; no BL occurs inside these inlined loops. No caller-visible scratch or source point is removed.

Chained results substitute producer output registers for the next kernel's left-operand loads. Scratch32..63 stores and loads are dead after inverse caching and are absent. The forward nonnegative direction-specialized difference stays in x11..x14; no zero test remains; defective rows complete as garbage. Multiplication temporaries are renamed away from that live bank, and its two right operand loads become direct reads. The initial Y load is postindexed by72 and its upper limb load uses offset-56; the trailing pointer ADD is omitted. The advanced pointer always completes the iteration. Branch labels are uniquely renamed. Arithmetic branch predicates are preserved because address substitution/renaming adds no flag-setting instruction. The witness checks all structural preconditions and exact generated source equality. Random-register/address differential tests check every raw output for membership in U(k) against independent original canonical loops with separately inverted signed prefixes, plus source memory, dead-scratch canaries, SP and LR.

## 12. Complete ledger and finishing deliveries

The expected total is81475742315604668696621117609394939345234instructions (the integer ceiling of the exact rational E of section15), with log2(E/1000) is below 125.937710102075634 . A64-term rational atanh enclosure and independent ln2 interval certify the outward claim125937710102075634/10^15 as the least multiple of10^-15 above the enclosure. Fixed windows, hop1 ACK events, hop2 DONE-flag completion and the latched stop flag are explicit machine assumptions; flag polling or general remote memory is not covered.

## Receiver hash join and the next-array

The receiver (route.s receive) joins every accepted record into its home's hash chains at arrival; there is no separate join pass. Each home owns2^60 head words (index K2 mod2^60 via the resident mask register, head address x9+8*index), zeroed once by zero4096 before the first group, one home record array with a register-resident cursor and end, and a next-array with one8-byte word per record position: the next word of the record at address A is at (A>>2)+C, C=next_array_base-(record_base>>2) (fixed_layout.py next_constant; the map is injective and cannot wrap below2^64). Records are appended in arrival order, so every accepted record is stored and inserted exactly once; records dropped at the full array are neither stored nor inserted.

Accepted path (15 including RET): the two payload LDPs and the flag LDR, which releases the window, with its branch (4); CMP/B.HS of the cursor against the last admissible position (2); a pre-indexed STP of (K2,K3) advancing the cursor and an STP of (w2,w3) (2); AND with the hash mask, head LDR and STR of the new record into the head, both register-offset (3); LSR and a register-offset STR of the old head into the new record's next word (2); CBZ on an empty chain and RET (2). Per visit: prior K2 LDR, CMP, B.EQ, then LSR, register-offset LDR of the prior's next word and CBNZ (6). On a K2 match: prior K3 LDR, EOR, TST of the low48 bits (d4) and B.NE (4); on a full match the advance is replaced by CMP/B.HS against the list end and, unless full, the post-indexed STP of (prior,new) (3, or2 when full), then RET. A full-array drop costs7. The receive row is15 per potential record; walks are charged6 per visit and4 per K2 match (section15); a listing adds nothing.

Same-index records form one arrival-order chain, newest first, every member reachable from the head: the new record is prepended (its next word is the old head) before the walk, which starts at the old head. The walk is uncapped and runs to the first prior agreeing on K2 and d4 or to the chain end; chains are finite (one home's records, at most CAP_H), so every walk terminates on every run and no overflow event exists. Each new record is therefore linked to its most recent182-bit-agreeing predecessor on every run, and the listed pairs are a subset of the182-bit-matching pairs at their home, so the candidate, request and relay rows (d_cl, d_rq, d_rl) apply unchanged. The candidate list is cumulative with its2^24 cap checked at every attempted append; the final DONE marker of each group stores the list cursor to the receiver's list-end cell before the release STLR, and the main context loads it after acquiring the DONE cell, so the group's new candidates are exactly the entries between the previous and the published end. routing_check.py executes the receiver on the true layout (persistent resident registers, a toy head table at the production addressing forms, planted cross-group same-key and same-K2/other-d4 pairs, full home arrays and full lists): heads, next words, listed pairs, the published list end and the exact per-call count equal the path model, every head-following reaches exactly its index's members newest first, and every same-key class is connected by listed edges while the list is not full.

## 13. Grid, masks and provenance

SIDE=2^35,N=2^70,q2=2^28 and the35+29/35+1+28 metadata allocation are unchanged. The home bits are[0,35),[64,99) and the hash-index bits are K2 bits[0,60); group144 and match182 include sixteen bits of the exceptional labels' top32constant bits. FILTER_FREE.md therefore enlarges the crude group/match coefficients rather than reusing a false uniform-bin assumption. FOURIER_CAPACITY.md stays the normal canonical-point lemma; the small extra fibers are bounded explicitly.

Arithmetic source44137210-aa9d-47fa-8cd0-7daad206a9d1, cold prices76ceba30-f1fa-4538-9c48-bc947cbd44bb and unallocated-padding removal37b8edde-addf-451f-9215-0dac6e7b6e79 are re-audited here; no acceptance is a premise. Inherited from7b724c27-6821-41d2-8524-dc95dbd928d6: filter-free labels, the DONE receiver, grouped early abort, hash chains and the section15 E-ledger. Inherited from the rjoin line: the home-array receiver with join at arrival (section12), ported from secp46385d36-072f-4ec3-bddf-b58376f50805 (via7c8264bf,03438c81,b63af001) and re-derived for P-256. New in this version: G=512 groups, both-pair operand reuse with carried inline emission, and the tightened hop1 row (section3).

## 14. Complete current cost and storage binding

range_ledger.py is the exact rational ledger. Its per-group2930/4-incidence rows include forward arithmetic and fused backward arithmetic with inline pair emission, weighted by the non-stopping bounds of section15; the4*carry_cap row pays every admitted high-carry correction on every run. RANGE_DRIVER.md separately bounds construction and every row/segment operation. The retained normalization allowance is not reduced. Both hops, the receiver's hash join, failed caps, sampling, source service, four explicit resolution stages, cold verification and output stay separately charged. RESOLUTION_PROTOCOL.md assigns the final output reduction to unused fixed phase4 control space, not an uncharged all-to-root fan-in.

fixed_layout.py gives simultaneous nonoverlapping addresses for every source, label, tree/sort, batch, candidate/request/relay/home/flag/context/port/cursor region, and context pointers are initialized to those exact bases. The former d4-counter range[2^51,2^52) stays reserved and unused; per-group hop1 slots, the hash heads, the home record array and the next-array follow2^52 and end at18115304426367778224<2^64. The sender header port and its local cursor table are separate allocations: the cursor-table base is a multiple of2^38 so the emitter forms a cursor address by inserting d1 into bits3..37, and the payload windows occupy an8192-byte stride per destination from an8192-byte-aligned base. The batch has a64-byte prefix reservation; output mailboxes follow the2^20 code/stack region with separate4096-byte child slots and isolated flag granules, and no flag is reinitialized after publication begins. Peak across2^70 nodes is below2^134bytes, conservatively declared below2^136.

The expected total and the outward claim are those of section12. The package requires external exact-version review.



The remaining named appendices are in proof-appendices.md and are part of this same claim.


## 15. Grouped early abort, stopping time and the expected-work ledger

**Groups.** The q1 X-rows of every node are partitioned once, before the run, into G=512 contiguous groups of sizes ceil(q1/G) or floor(q1/G) (range_ledger.py group_sizes; q1=377464334 gives737236 and737235). All nodes run group g at the same time, separated by barriers: hop1 of the group's rows into fresh per-group slots (CAP_PG=5*2^17 records, BLOCK_PROTOCOL.md), hop2 into the cumulative home record arrays with the receiver's hash join at arrival, the four resolution stages for the group's new candidates, and the stop test. The stop flag is a latched word set only by a home whose resolver has verified [x]P=Q, broadcast with the barrier release and tested by every node before the next group (5 instructions, phase-control row); a set flag halts the run with that x as output. No group is skipped for any other reason. Rows, scalars, points, normalization, the range tree, coins and caps are those of sections2 to6; only the execution order and the detection pass change.

**Stopping time.** Let T be the first group whose resolution confirms an x (T=G+1 if none). Group g's arithmetic, emission, hop1, hop2 (with the receiver's append and insert) and per-group fixed work are deterministic given g and are executed iff T>g-1. The success event is {T<=G}, and the section6 expression is evaluated for the full run: the chain detection does not depend on group boundaries (every record walks the chain of all its priors at arrival, every accepted record is inserted exactly once, the candidate list is cumulative, and a group's candidates are resolved before its stop test). Thus Pr(T<=G)>=0.390000001 at q1=377464334, the least such q1 (range_ledger.py).

**Non-stopping bounds.** Write m_g for the number of rows in groups1..g-1 (m_1=0). Pr(T>g-1) is at most the probability that no useful collision among the first m_g rows of every node is resolved by the end of group g-1. The section6 expression re-evaluated at q1=m_g (same N and q2; CAP_PG, CAP_H, the carry cap and the candidate, request and relay caps fixed at their full-run values; calc(m,caps) in range_ledger.py) lower-bounds that resolution probability: the Chen-Stein bound and every defect, capacity and usefulness row condition only on coins and records of prefix rows, prefix means are at most full means, the fixed caps keep every Chebyshev gap at least its prefix value, and later groups never unlink a prefix record from its chain nor remove a listed pair. Hence PHI[g]=min(1,1-lo(m_g)), rounded up to a multiple of2^-120, bounds Pr(T>g-1); PHI[1]=1, nonincreasing, PHI[512]<0.61118; the group E-factor sum_g PHI[g]*q_g/q1 is below0.85743.

**Expected work.** By linearity over the indicators1[T>g-1], E[work]<=W_rest+sum_g PHI[g]*W_g+W_sameX+W_K2+W_def, where W_g is the group's per-row and per-record work (arithmetic, emission, the39479 inversion, the1024 row allowance and the2^17 driver per row; hop1, hop2 and the receiver's append and insert per record) plus its main-pair walk bound (below), W_rest holds the node-initialization, hash-head-zeroing and range-tree rows, the per-group fixed work of all512 groups, the four phase-4 stages per group, generation, sampling, normalization/sorting, the full carry meter cap (4*carry_cap) and the resolution, source-service, barrier, control and output rows, all at their worst case, and the global rows are charged without PHI. No independence between stopping and the chains is claimed.

**Walk means (history-conditional domination with a validity split).** Let H_g=(all rows of groups<g at every node, all Y generators, every node's b and R_g, Q). The event {T>g-1} is a function of H_g: records, chains, listed pairs, requests, replies, the resolver's sign test, the usefulness test and [x]P=Q are deterministic functions of those points, scalars and Q. Group g's walk cost at a home is sum over its inserts u of the visits V_u, and V_u is at most the number of priors v in u's hash chain (every visit examines one distinct chain member, section12). Visits are counted against the130-bit cell (home70 bits plus the60 hash-index bits of K2, bits[0,35),[64,99),[128,188)): every prior in u's chain lies in that cell, and one agreeing on the index but not on K2's top four bits or on d4 is a plain visit or a K2-only check. Split the pairs (u,v): main pairs, whose leaves come from different X-rows (in earlier groups or in group g itself), and same-X pairs (same row, hence the same group and node). For a main pair condition on H_g and on every group-g row except u's; v's row is then conditioned whether it lies in an earlier group or in group g, and v's actual stored label, a deterministic function of its own row's batch (X_v and all Y) and the earlier transport, is one fixed label. If u's row is not defective, u's record label lies in U(k) for the exact key k of the class of X_u+sY_u, where X_u=R_g+[a_k]P with a_k fresh, uniform and independent of everything conditioned on: given Y_u and the sign, X_u+sY_u is uniform on the group, so the probability that U(k) meets the one fixed130-bit cell of v's label is at most CELL130=(2*2^126+4*2^94)/(r-1) (the crude two-point fiber bound with the exceptional-label fibers of FILTER_FREE.md, section3, h=0 fixed positions in[224,256)); X_u=O (probability1/r) leaves a fixed leaf. So Pr(v in u's chain, u valid | conditioning)<=PAIRM=CELL130+1/r uniformly, whatever v's own validity. Which representative of U(k) u's batch arithmetic emits is irrelevant: the event uses only membership in U(k). If u's row is defective (unconditional d_batch per row), u's label is garbage and its walk is at most CAP_H visits (the home array holds all of its chain), each possibly K2-verified: charged globally as DEF=(6+4)*d_batch*CAP_H*S. Hence E[1[T>g-1]*walks_g^main]<=PHI[g]*6*PAIRM*(s_g*p_g+s_g(s_g-1)/2-sameX_g)+DEF, with s_g the group's leaves over all nodes, p_g the earlier groups' leaves and the pair count deterministic (arrivals never exceed sends): the conditional bound is established before multiplying by PHI[g]. Same-X pairs are charged globally without PHI and WITHOUT fixing either actual label: the two records share their row's batch-prefix arithmetic, so neither raw representative is fixed when the other's generator is exposed. Both actual labels lie in the fixed sets U(k_u) and U(k_v) of their ideal classes, so the records share a chain only if some130-bit cell meets both sets, and Pr<=sum over cells C of Pr(U(k_u) meets C and U(k_v) meets C). For two leaves with different Y, condition on X: the ideal points are independent uniform (Y_u, Y_v fresh), so the sum is at most max_C Pr(U(k_v) meets C) times the expected number of cells met by U(k_u), at most2*CELL130. For the two leaves of one incidence, (X+Y,X-Y) is an invertible linear image of(X,Y) with Y nonzero, uniform on pairs of distinct points, so given the first leaf the second is uniform off one point and the same bound2*CELL130 holds. With X=O (probability1/r) the two leaves of an incidence share a class: PAIRM_SAME=2*CELL130+1/r (SXV in the ledger; the same two-cell union FILTER_FREE.md uses as M=2*KB182). K2 verifications drop the stopping indicator (E[1*X]<=E[X]) and are charged globally over ALL unordered pairs, main and same-row alike, at the two-cell union2*KB(134), KB(134)=(2*2^122+4*2^90)/(r-1): EK2=4*PAIRS*2*KB(134); for main pairs the single cell would suffice by the conditioning above, and the union dominates both cases. A visit costs its probe and advance (6); a listing replaces the advance at no extra cost and the list-full case costs one less, so no listing row is needed.

**Equivalence with the ungrouped success analysis.** Same records (same kernels and row order; groups are contiguous row ranges; hop1 cursors reset per group as formerly per phase); same home-array contents in arrival order (cumulative cursor, cap and drops); every accepted record inserted exactly once at arrival; same-key records chain-connected unless a182-bit intruder arrives between them (the intr row); listed pairs resolved by the same protocol; a useful class edge exists whenever a useful pair exists (section5); stopping only skips groups after a confirmed x. So the section6 bound holds for the grouped run, and the E-ledger charges confirmed stopping only.

**Rows.** range_ledger.py uses the protected subgroup order r (asserted below the field modulus) in P1, rho, every KB bound, the prefix success bounds and PHI. It evaluates, exactly in rationals, W_g=NODES*q_g*(2930/4*Q2+39479+1024+2^17)+s_g*(584/128+29/4+15)+6*PAIRM*(main pairs of g); globally without PHI, SXV=6*PAIRM_SAME*(same-row pairs), EK2 and DEF as above; the worst-case group_fixed row NODES*512*(1536*SIDE+64), hash-head zeroing2^48*1026+1, and the unchanged worst-case rows of section9 (four phase-4 stages per group,2048 barriers, phase control2^20, the per-group stop broadcast in the output row). The expected total E has integer ceiling81475742315604668696621117609394939345234, log2(E/1000) is below 125.937710102075634, and the claim125937710102075634/10^15 is the least15-decimal value at or above the rational enclosure (score_bounds). With PHI=1 everywhere the rows give less than2^126.1586; with every walk at its CAP_H worst case and all groups run, the disclosed per-run worst case is2^178.3505 units (an outward4-decimal ceiling of a certified logarithm upper bound, asserted in range_ledger.py). Preprocessing (node init, generators, sampling, zeroing, range tree) is below2^119 units, included.

## Short Montgomery multiplication and exact cost change

The new scalar multiplication core is adapted from Amazon s2n-bignum commit56152cdb51f47ef9f63b179cda0a8b7518379240, arm/p256/bignum_montmul_p256_alt.S. Copyright and Apache-2.0 OR ISC OR MIT-0 licensing are preserved in p256_fields.s and the repository LICENSE. Primary source and its corresponding arm/proofs/bignum_montmul_p256_alt.ml are available at https://github.com/awslabs/s2n-bignum/tree/56152cdb51f47ef9f63b179cda0a8b7518379240 . This attribution is not a claim that the upstream HOL proof automatically certifies our patched tails or register transformations.

The core first forms the full four-by-four schoolbook product with MUL/UMULH and explicit carry propagation. Its low limbs are reduced four times using the cancellation digit t and p's limbs from section3. Writing MC=2^64-2^32+1, the seven-instruction reduction adds t shifted32 to the next limb, carries t>>32 onward, and adds the low/high halves of t*MC to the remaining two limbs with carry. Each step is exact division of the current integer plus t*p by2^64. Four steps give raw T=(ab+m*p)/R0 with0<=m<R0; uniqueness of the cancellation m modulo R0 shows this T equals the previous multiply's raw value for all a,b<R0. The final high carry is consumed immediately by B.CC. The weak tail keeps T if T<R0 and returns T-p otherwise, always in[0,R0). Old and new correction predicates agree pathwise.

Including RET and standalone tail constants, MUL costs110 instructions when correcting and104 otherwise (probe short_mul); the generated loops hoist the one MC materialization, so before unrolling the rolled bounds are forward122 and backward616+4*carry, not an average branch extrapolation. For a zero denominator the full forward pass completes with bounded garbage within the441/4-instruction incidence allowance and its fixed46; the defective event keeps its already-charged d_batch exclusion. All attempted rows stay charged to the uniform full-row upper bound.

The new RHS limbs occupy x7,x8,x13,x14 and are read-only after loads. MC occupies a distinct scratch register standalone and the reserved x2 in generated loops. The forward specialization maps the live canonical difference into that RHS bank and conflicting scratch into x24..x27. Backward slope products substitute cached inverse x24..x27 for all four RHS registers. The original fused square and its raw-label case split are retained.

The inherited rolled listing is pinned to SHA25607f98f57a417d080a5ccdf49fc8d1381b501b174251f89606e5ea69a48b95a29; this version applies the exact eightfold transformation of section3, giving2930/4 per incidence plus4 per carry, with the driver's2^17 row covering segment fixed terms. The generated listing SHA256 is899b087806f02cb67232fc350b27fe05760f6f33cd98b3ad7be1110e731cfcee, encoded size102456 bytes.



