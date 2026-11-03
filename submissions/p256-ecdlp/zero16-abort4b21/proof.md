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

The complete claim is expected total classical
work below 2^126.075215663167288 units at the protected C=1000 scalar instructions/unit (section 15: four X-groups, early abort on a verified x; every run stays below 2^126.232872),
success at least 0.39, peak memory below 2^136 bytes, preprocessing below
2^118 units included in the total, and zero nonuniform advice. This is an
astronomical-memory analytical attack, not a feasible machine demonstration.

## 2. Random construction and pair arithmetic

There are N=2^70 nodes, indexed (i,j) with each coordinate in [0,2^35).
At each node draw b and q1=377471015 independent a_k uniformly modulo r,
and q2=2^28 independent c_l uniformly in {1,...,r-1}. Coins and draws at
different nodes are independent. Each scalar draw tries at most eight fresh
256-bit values, rejecting out-of-range values and zero where required. A cap
failure aborts the run. The rejection probability is below 2^-31, so the
union cap loss is at most N(q1+q2+1)*2^-248. Sampling is charged 256 scalar
instructions per requested scalar. Under the protected ideal-coin convention,
RNDR supplies the bits; a hardware-failure assumption is not silently added.

Initially construct R_g=[b]Q, X_k=R_g+[a_k]P and Y_l=[c_l]P. These original points and indices are used in the probability argument below. Apply the exact signed normalization and paired scalar updates in NORMALIZATION.md: store X'_k=epsilon_X X_k with a'_k and row-specific b'_k, and Y'_l=epsilon_Y Y_l with c'_l. Build the y-partition/x-sorted range tree of72-byte records with original Y IDs; retain canonical source tables at original indices. Keep original X index k and X-group membership. Preserve original b while constructing X, separately from stored row-specific b'.

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
outside the allocation. The reference point formulas and the executed interpreter tests check these invariants independently.

## 3. Field arithmetic, range invariants and worst-path counts

Let R0=2^256 and delta=R0-p. Original generator coordinates are canonical Montgomery residues. Only M/S receive weak operands belowR0; common x-sum and y-plus addition receive canonical original coordinates. Canonical cold arithmetic follows COLD_ACCOUNTING.md.

For radix B=2^64, p has limbs(B-1,2^32-1,0,B-2^32+1). In each reduction step the current low word t is the cancellation digit, and(T+t*p)/B is an integer. Split t=t_low+2^32*t_high: the corrections use t*2^32 and t*(B-2^32+1), with explicit carries and borrows. Four schoolbook rows and these reductions compute raw REDC exactly. The symmetric square uses each diagonal once and each cross-product twice. For a,b<R0, raw T=(ab+mp)/R0,0<=m<R0, satisfies0<=T<R0+p<2R0. Its high limb is0 or1. Keeping T when high=0 and subtracting p when high=1 gives an exact representative in[0,R0). The carry-only MUL and SQR tails implement precisely this rule, with at most110 and88 instructions including RET for the new standalone MUL and retained SQR. Canonical-input ADD has sum<2p<R0+p and the same rule costs at most19. This ADD is never used on two weak operands. Canonical SUB stays19 in the standalone reference tests.

The range-specialized absolute difference subtracts the smaller canonical integer from the larger in the descriptor's order, producing the nonnegative magnitude without a per-item negation. Its magnitude is belowp and represents one of the two field signs. Both the denominator sign and the y-minus numerator sign disappear when the final slope is squared. No zero test remains; defective rows complete as bounded garbage within the same envelopes, detailed below in this section. The y-plus branch and leaf tags/recovery coefficients are unchanged.

Direction specialization is machine-checked limb by limb. The pair_fwd calling convention passes the Y stream in x1 (copied to x20) and the X point in x3 (copied to x22), so the inlined four-limb subtraction reads the Y minuend into x3-x6 and the X subtrahend into x7-x10. The unreversed pair_fwd_hi therefore computes Y-X and is selected exactly when Y.x>=X.x; pair_fwd_lo swaps all four limbs (chain_witness.py normalize_transform, one replace per limb) to compute X-Y and is selected exactly when Y.x<X.x. range_witness.assert_limb_specialization asserts on the generated program that the lo section holds exactly the four swapped limbs and none of the forward ones, the hi section exactly the four forward limbs and none of the swapped ones, and each of the four backward y-partition variants exactly its four expected limbs. The probe additionally executes the generated pair_fwd_lo on X.x=2^128 with Y.x=1 (plus its mirror and eight further high-half/boundary operand pairs) and asserts the stored denominator equals the exact absolute difference in every case (2^128-1 for the headline pair).

The fixed inverse computes exponent p-2 with255 squares and127 multiplies. Induction on the preceding REDC range argument proves all intermediates stay belowR0 and congruent to the specified modular powers. The fixed inverse instruction sequence is retained byte-for-byte after macro expansion. It calls external k_mul/k_sqr; the experiment links the changed hot routines. Its conservative39479 allowance is retained without crediting their smaller counts. Canonical cold point operations instead use the separately specified namespaced canonical routines. On non-defective rows its input is a nonzero prefix product and no zero-residue representative0 orp arises; defective rows complete as bounded garbage within the same envelopes, detailed below in this section. Fermat's theorem therefore gives the field inverse. Each changed M/S path is no longer than its original canonical path, so39479 stays a uniform upper bound including copies, prologue, calls and return; no average branch discount is taken.

For final fused output write R0=2^256, delta=R0-p. Raw square T<R0+p and z=T-c lies in(-R0,R0+p), so h=floor(z/R0) is-1,0,1. The four-limb stored value is u=(z-h*p) modR0. For true canonical key k=z modp, interval splitting proves u belongs to U(k)={k}, with k+p additionally possible when k<delta and k+delta additionally possible when k>p-delta. These bands are disjoint and |U(k)|<=2. Every noncanonical label lies in[p,R0), including the noncongruent negative-wrap branch. FILTER_FREE.md gives the full case split and the revised load argument.

There is no emitter representation filter. A useful collision whose canonical key is in either exceptional band is charged by the retained unconditional Gamma*2^227/r^2 witness-loss allowance. Other exceptional labels are neither discarded nor assumed independent: fixed point-only set predicates dominate their bin, intruder and pair loads. There is no global S*2^-32 omission.

The generated loops reserve x0 and x2 for the zero-extended constants0xffffffff andMC=0xffffffff00000001, initialized once after original count/pointers are saved. After pointer specialization, the active field-data register bank is x1,x3..x17; it excludes both constants and the separate cached-sum bank. The y-plus carry temporary moves from x0 to x1. MOV does not change flags, so removing per-correction constant loads preserves arithmetic predicates. For OLD=0x00000000fffffffe=~MC modulo2^64, ADC(y,OLD,C) and SBC(y,MC,C) have identical result and NZCV when corresponding flag-setting forms are used; conversely SBC(y,OLD,C)=ADC(y,MC,C). Every original CONST_B consumer is transformed this way before MUL MC uses are introduced. The MC constant is a single encodable logical-immediate MOV, hoisted once per nonempty loop. The structural witness checks the exact removed-MOV counts and rejects any generated loop instruction writing either reserved constant register. Initializing them costs two instructions per nonempty loop; empty-loop paths skip this work.

The fused operation branches on the final raw-square ADCS carry before any flag-changing instruction. If raw carry is0, subtract c and keep on no borrow, otherwise add p. If raw carry is1, subtract c and keep on borrow, otherwise subtract p. The signed-high value would be rawcarry-borrow, so this is the same mathematical correction without materializing the fifth limb. The negative low path branches over the high path. Labels are unique per expansion. After the final raw-low ADCS, the four carry/borrow cases cost8,13,12,8 before common stores; no uniform dynamic saving is assumed.

Each direction-specialized forward entry costs at most441/4*n+46; the per-item zero test is removed, so defective rows complete with bounded garbage instead of aborting. Main-loop items chain the prefix through alternating register banks (no A-loads); the at most 7 pinned tail items keep memory loads and are covered by the fixed premium. Pair emission is fused into the backward kernel: each incidence emits its two records inline to block windows (9 instructions for the first record and 10 for the second, including the one-instruction TBZ full-window test, in place of the two output stores), so the fused backward bound is2489/4*n+33 and there is no separate emitter pass; the v40 base 616*n+11 plus4 per carry-1 incidence is retained inside it. Every new generated eightfold loop keeps the same control skeleton in the same order; forward mains alternate even/odd bank-chained copies. A nonempty backward/emission call costs exactly4-14*floor(t/8) more than v40, a nonempty forward call exactly6-30*floor(t/8) more, while empty calls are unchanged. Because14*floor(t/8)>=7*t/4-49/4, calls use2489/4*n+33 fused-backward; forward uses441/4*n+46 (main items 110 worst, tail singles 114 worst, prologue/epilogue/tail overhead 19, so 46-19-3.75*r>=0.75 for every remainder r<8). Sixty-three fixed instructions per segment (46 forward incl. tail premium and bank preload, 17 backward) over at most87 segments are covered by the2^17 per-row driver charge below; no separate63*87 term is charged. The row still pays every garbage/empty/small-call case: defective rows complete within the same envelopes (all kernels preserve below-R0 on zero inputs, the Fermat inverse has no division, all addressing is sequential, garbage candidates fail final verification, and the defective event keeps its already-charged d_batch exclusion); empty calls retain the prior fixed bound.

Full-flush publishes go through the shared flush_publish_b subroutine (16 taken / 13 dropped, including its own d1 mask, plus 3 call-site, unit-tested on 32 cases), charged at 19 per firing in the hop1 row; that row charges2213 per 512-leaf flush (2183 block_dlv+10 in-flight+19 publish+1 reserve, each measured or itemized), i.e.2213/512 per leaf. The arithmetic/emission variable price is now2930/4 per incidence plus4 per counted carry. b16 renames the first backward multiply's accumulator limbs to x1,x17 so its low operand x3,x4 survives; the second multiply drops the `ldp x3,x4` reload, one instruction per backward incidence. b18 extends this to the high operand pair: in the first product x6 is never written and x5 is first written only as the top accumulator limb (`umulh x5,x6,x14`), after every operand read of x5, so renaming that accumulator role to x28 (dead there: the y-plus cache x28,x29,x30,x18 is produced later in the item and consumed before the next item) keeps x5,x6 intact and the second product drops its `ldp x5,x6,[x23,#16]` reload, another instruction per backward incidence; the generator asserts that x6 is never a destination, that x28 is unused in the first product and that no x5/x6/x28 access precedes the deleted load. b17 shortens each inline record emission by three instructions (2489/4 per incidence with b18, b19 and b21): the full-window test `tst x14,#0xfe0; b.ne` becomes `tbz x14,#14` (b20: payload windows hold512 records and sit at a32768-byte stride from a32768-byte-aligned base, so bit14 of a slot address is set exactly on the512th slot, where bits5..13 are zero); the d2/sign/ID word is formed by `bfi x15,x8,#0,#35` into the loaded ID word (whose low36 bits are zero) instead of `and`+`orr`; and the cursor-table address is formed by `bfi x13,x7,#3,#35` into the2^38-byte-aligned cursor-table base instead of `and`+`add`. Both inserted fields are the same35-bit low halves the old masks selected, so every stored word and address is unchanged. The unmasked d1 register is masked once inside the shared flush subroutine (+1 per firing, charged above). All three forms assemble under clang as single AArch64 instructions; no shifted carry form is used. b19 carries two values from the first record's emission to the second: the loaded ID word stays in x1 (its low35 bits, overwritten by BFI, carry d2 and are overwritten again), and the LDADD addend32 stays in x17, so the second emission drops its `ldr` and its `mov #32` (two instructions per incidence). Between the emissions the slope product (index8) never touches x1 or x17, and the fused square (index9) has its x1/x17 scratch renamed to x24/x25, the inverse-bank registers that are dead once that product has consumed them; the generator asserts both facts on every variant. b21 also carries the cursor-table base (x13) and the W2 word (x14) across: the first emission loads both with one `ldp x13,x14,[sp]` (LDADD result in x15, w3 in x1), the second re-inserts its d1 into x13 by the same BFI (field bits overwritten, base bits above bit37 untouched) and forms its w3 in x9, dropping its `ldp` (one per incidence); the flush subroutine takes the slot in x15 and its remote-cursor result in x3, writing neither x13 nor x14, and the index-9 square's x13/x14 scratch becomes x26/x27. The generator asserts index8 touches none of x1,x13,x14,x17 and x24..x27 are absent from index9 before the rename; fused_check asserts on the fused listing that between the `ldp` and the second BFI nothing writes x14 and only the two BFIs write x13. The probe checks `ct = foct + 4 - 14*floor(n/8) - 2n` against the rolled fused base and identical final memory. One39479 inversion and old1024 fixed row allowance remain. Carry-meter semantics, failure handling and cold construction are unchanged. Hop2 completion uses DONE-flag signaling (no per-message acknowledgment or in-flight updates); see RECEIVER_PROTOCOL.md. Section 15 combines these constants into the exact totals.

The fused transform first splices inline window emission plus a shared flush subroutine into the four rolled backward bodies (frame spills, machine-checked liveness peak 30 with publish scratch free on all variants), then the loop transformer unrolls both forward entries and all four fused backward variants: it checks that the body never reads or writes count x19, makes eight uniquely labelled copies and uses a biased unsigned count. Initially subtract8 and skip to the tail on borrow. After each group subtract8 and repeat on no borrow. At exit add8 to recover a remainder0..7, then run the unchanged single-item loop. A zero initial count uses the unchanged early return. Internal arithmetic branches are renamed per copy; each item's pointer increments and output order are unchanged. The full finite replay compares every tested segment's entire memory and exact count against frozen v40, with a routing differential on flush/no-flush cases.

The generated arithmetic listing has102456 encoded bytes, counted by parsing its scalar instructions; all .s helpers are macro-expanded and counted by the same parser, and their encoded bytes plus the102456 generated, 1860 unrolled-block and65536 stack bytes are asserted within the1MiB code_stacks allocation. The emitted program is a public algorithm constant; per-node initialization materializes each32-bit instruction with two immediate moves, one store and one pointer increment, inside the1MiB code initialization charge, with no runtime generator, uncharged expansion, lookup or dispatch.

## 4. Fully charged transport and grouping

Use the explicitly specified two-dimensional blocking ports of RECEIVER_PROTOCOL.md. Hop1 uses the explicitly fixed512-record window, native2183-instruction delivery/ACK and transaction-bound completion protocol of BLOCK_PROTOCOL.md. Its full-flush work is charged at2213 per block, and all partial/initialization/finishing work has a separate8192*SIDE per-node reserve per X-group. Ordinary final-cursor, barrier and resolution messages retain their addressed delivery cost. Hop2 uses a separately charged resident scalar receiver context at each node:14 instructions for accepted records,11 for full-bucket drops, with no per-message acknowledgment or in-flight updates; its resident cursor-table base register (a multiple of2^30) and counter base register (2^51) receive d3 and d4 by one BFI each, so no mask or address add is executed. Each sender closes its group stream with one DONE marker per destination; the receiver counts markers and signals per-group completion with a single STLR to its DONE cell, which the main context awaits with the same LDAXR/WFE completion primitive as barrier epoch words. Both contexts' work and state are included. The hop2 receive window is released only after the final flag load of the fixed five-word packet, when all four payload words and the flag are already in receiver registers; the two earlier payload loads never release it, so a publication cannot overwrite a packet being read. No direct-home shortcut, uncharged key operation, software polling or unbounded queue is assumed.

Key digits stay d1=K0 mod2^35,d2=K1 mod2^35,d3=K2 mod2^26,d4=K3 mod2^48. The record is(K2,K3,w2,w3), w2=i|(k<<35), initially w3=d2|(sign<<35)|(l<<36). After extracting d2, h2s replaces only w3bits0..34 with source column j. Thus a home reads source(i,j) directly; sign/l and the absolute k are preserved. k<2^29,l<2^28, including maximum-boundary and post-emission-wrap cases.

Partition each node's X indices into the four contiguous X-groups of section 15, at most qg=123857676 indices each; do not renumber k. For each group initialize fixed payload windows and local/remote cursors, emit ordered512-record blocks into sender-row slots at(d1,j), flush each remaining partial block with initialized padding, send clamped valid final cursors, drain and barrier. Padding is excluded by those cursors. Each slot has CAP1=25*2^18 records. Stream all its stored records to(d1,d2); no hop2 pair slots exist. The destination receiver computes d3, reserves32bytes in its tight full-record bucket, compares with its explicit end, rolls back on full, otherwise stores the four payload words and increments d4. After streaming, each sender publishes one DONE marker per destination; each main context awaits its receiver's DONE signal, then all nodes barrier before reusing hop1 storage.

Buckets and d4 counts persist across groups; only receivers write them. After each group drains, prefix4c writes d4 range starts into a separate cursor copy and pass4 scatters pointers to all records so far, in d3 order, through it. The full-record scatter costs6 per record plus2 per8-record block. Final d4 cursors are exclusive range ends. The unchanged paired scanner checks only K2-adjacent pairs inside each d4 range, at most(195/64)*records+26*ranges+8*(CAND_CAP+1)+2 instructions. Cap CAND_CAP=2^24 is checked at every attempted append.

Consequently matched keys fix 70 home bits, 48 d4 bits and all 64 K2 bits:
182 distinct bits. The d3 bits are already within K2. Partial matches are
not treated as true group collisions. Each is resolved and verified.

Each scan's list is resolved by the four explicit bounded64-byte-record phases of RESOLUTION_PROTOCOL.md, not count/offset exchanges. Source and home tags use35-bit node coordinates; X index29bits, Y index28bits and sign1bit. Source tables keep original indices and row-specific signed b. The protocol includes initialization, caps, partial-reply flags, allocation, source gather and a final bounded fan-in output tree. Its complete work is charged by four N*2^41 rows and source service1024 per capped request. Source cap2^29 and relay allocation2^34 records are fixed; the original stricter logical relay loss bound is retained as explained there. A home reconstructs X+sY and X'+s'Y', tests equality up to sign, solves the nonzero scalar denominator and verifies [x]P=Q. Scalar arithmetic uses the complete COLD_ACCOUNTING.md bound20592111 per candidate.

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

Condition on the witness's four generators, but never on global cleanliness or actual raw-label choices. Retain the reservation12 fixed leaves. Its hop1 group has at most qg indices, mean m_p=2*qg*q2*KB35+1 and variance bound v_p=m_p+1+60*(qg+q2)*KB35. Its bucket spans all groups, mean m_3=S*KB96+1 and variance v_3=m_3+1+60*(q1+q2)*KB96. Then d_p=v_p/(CAP1-12-m_p)^2 and d_3=v_3/(CAP_3-12-m_3)^2, CAP_3=ceil(m_3)+2^30. The witness has one predetermined group, so no extra group multiplier is needed. There is only one pair-slot loss, because no second-hop slot exists.

Here is the variance justification. Leaves sharing a free Y but different
free X are independent. With a common free X and different nonzero Y, the
only missing value is X, so covariance is at most 1/(4(r-1)^2); if one Y is
fixed, it is at most KB(f)/(r-1). Two leaves of a free incidence are an
invertible linear image of (X,Y), Y!=0, and their covariance is bounded by
the source's cell-size term divided by r^2(r-1). The aggregate of these tiny
terms is below one at the stated S,q1,q2,r, even using a global S upper bound.
An incidence involving a witness generator and one free generator can have
fully correlated leaf indicators; counting all such ordered partners with
the generous60(qg+q2) factor for its group slot, or60(q1+q2) for the full bucket gives the remaining variance term. Chebyshev
then yields the displayed conditional drop probabilities. The exact script asserts CAP_P-12>m_p and CAP_3-12>m_3 before division; a negative gap cannot be squared into a valid bound. No independent
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
bound exceeds 0.3900000021509504 at q1=377471015 with qg in d_p (the candidate cap2^24 raises d_cl; the smaller largest group lowers d_p). The declared success is conservatively 0.39. None of this identifies reviewer confidence with algorithm success.

## 7. Work and memory

range_ledger.py contains the complete exact rows, including both hops, every delivered record, empty tables/ranges, all failed caps, capped scalar sampling, preprocessing, phase4, resolution, verification and output. No unallocated per-leaf bookkeeping reserve is charged: all record operations are itemized below and all driver/context work is charged per batch or group (section 15). Work is summed across both scalar contexts at every node.

Cold generation/header work is charged1623023instructions per X generator item (1577967 per Y item) and candidate resolution20592111instructions per capped candidate, including failed verification; the scalar multiplication inside both is the itemized Jacobian double-and-add ladder of COLD_ACCOUNTING.md (1315823, one inversion per point) in place of the former 512-affine-operation envelope. These are the explicit itemizations in COLD_ACCOUNTING.md; the shared-per-node generator and exact scalar-inverse refinement follows the itemization in bdf01aa5-16e3-4619-a8eb-1aa0f6b6b29b (its acceptance is not a premise). Preprocessing stays below2^118units and is included.

Reserve the first2^52bytes of each node for code, both contexts/stacks, packet/ack state, DONE cell/granule, per-scan mailboxes and stop words, d4 counters, scalar/X/Y tables, batch storage, source/relay/candidate buffers and control tables. The fixed-region calculation in range_ledger.py includes the extra d3 end table and DONE granule and stays below that reservation. Place the phase receive region next,32*SIDE*CAP1bytes; then2^26 tight full-record buckets of32*CAP_3bytes; then an8*2^26*CAP_3 pointer array; then the2^51-byte scan cursor copy. All regions are nonoverlapping and eight-byte aligned. The exact endpoint is below2^64, and the attempted hop1 cursor bound cannot wrap. Bucket cursor temporary increments are at least the full pointer-array region below2^64; a full cursor is restored before the next packet is consumed. This is a simultaneous-allocation bound, not an assumed lifetime reuse saving. AcrossN=2^70nodes total storage is below2^134, hence the declared2^136bound. No nonuniform advice is used.

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
full-block flush, delivery, in-flight and publish counts at2213*S/512 in total: ordinary per-incidence sender emission is already included. BLOCK_PROTOCOL.md proves the block bound and separately charges8192*N*SIDE per group for initialization, partial flush, finishing, endpoint registration and control.
Hop2 sender costs7 instructions per record plus2 per8-record block,11 per source slot and2 fixed per phase. The dedicated receiver costs14 accepted/11 dropped, with no per-message acknowledgment or in-flight updates; charge14 for every potential record. DONE-marker publication costs3+6*SIDE per sender per phase, DONE handling costs8 per non-final and10 per final marker, and the DONE signal/wait/reset costs are charged per node per group; all sit in the per-group term GRP. Full-record scatter costs6+2/8; the unchanged scanner costs195/64. These ordinary grouping terms total595/64, and hop2 totals85/4 per potential record. No old bucketing pass is charged or executed.

Each group needs one hop1 cursor initialization and one finishing sweep, charged3*SIDE+1 and17*SIDE+3 per node. H2 source-slot traversal occurs once per group. Initializing each bucket's cursor/end pair costs5 once, and pass4 adds11 per bucket per scan. Zeroing costs4 per four d4 counters once; prefix4c costs10 per four per scan. The fixed returns of init3, pass4, zero4 and prefix4c are charged. The scanner fixed/candidate terms are unchanged. Sixty-four barriers at1024instructions/node and another2^14 phase/context-control allowance per node are charged per group. RECEIVER_PROTOCOL.md itemizes their scope, including resident startup and the final blocked receive. The block row includes a1-instruction margin above its2212-instruction full flush; otherwise source address/tag extraction, all cap tests, receiver release, scatter and scan work are already in these rows. Per-batch setup and phase/context control stay separately bounded.

The remaining rows explicitly charge all generator points and headers,
eight-attempt scalar sampling, program delivery and node initialization,
all capped resolution, source request service, four explicit bounded-routing stages,
barriers, target validation and output. Reply record data is fetched
from stored tables; no scalar multiplication is assumed free. Source serving
costs at most 1024 per request including addressed copies. Final resolution
has the20592111instruction itemized bound justified in Section7 and COLD_ACCOUNTING.md, even for false
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

The new block protocol preserves every retained sender prefix and downstream order, so the exact success arithmetic is unchanged. range_ledger.py holds the exact integer sums, the rational logarithm comparison (not a floating logarithm), a preprocessing check on a superset of setup rows, all included in total work, and the layout assertions.

## 10. Evidence, scope and provenance

Field/body counts, full batched arithmetic with inversion, metadata
emission and denominator/cap controls run through the supplied
AArch64 interpreter in the protected runner. Tests compare against
independent Python integers and affine point formulas.
The two-hop test executes deliveries and data-path kernels on a 2-by-2
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

chain_witness.py and range_witness.py deterministically define the generated range program from the original call-based schedule and supplied new field/fused bodies. It specializes memory addresses before bijective register renaming. Stable pointers/counters are x19..x23, except that the forward loop no longer keeps the X point or previous-prefix pointers live: X is loaded once in the prologue into x7-x10 (which the renamed multiply never writes) and the prefix is read relative to x21 (pinned tail) or the ping-pong banks (mains), so x22/x23 are unused there and their setup MOVs are deleted. The first backward per-item inverse is produced directly in x24..x27 and stays live across prefix update and numerator arithmetic; both slope products read that bank in place of their right-operand loads. The y-minus operation leaves its original coordinate inputs intact; y-plus reuses those input registers and caches its result in x28,x29,x30,x18 (x28 also serves as the first backward product's top accumulator limb earlier in the same item, before the cache is written). The first slope/square/fused chain does not touch that bank. The second slope consumes it. Original LR is saved on entry and restored on return; no BL occurs inside these inlined loops. No caller-visible scratch or source point is removed.

Chained results substitute producer output registers for the next kernel's left-operand loads. Scratch32..63 stores and loads are dead after inverse caching and are absent. The forward nonnegative direction-specialized difference stays in x11..x14; no zero test remains; defective rows complete as garbage. Multiplication temporaries are renamed away from that live bank, and its two right operand loads become direct reads. The initial Y load is postindexed by72 and its upper limb load uses offset-56; the trailing pointer ADD is omitted. The advanced pointer always completes the iteration. Branch labels are uniquely renamed. Arithmetic branch predicates are preserved because address substitution/renaming adds no flag-setting instruction. The witness checks all structural preconditions and exact generated source equality. Random-register/address differential tests check every raw output for membership in U(k) against independent original canonical loops with separately inverted signed prefixes, plus source memory, dead-scratch canaries, SP and LR.

## 12. Complete ledger and finishing deliveries

Section 15 states the exact expected and worst totals; a64-term rational atanh enclosure and independent ln2 interval certify the outward claim. The hop1 finishing sweeps include delivery and in-flight work; no second-hop finishing table is assumed. All cold rows, raw emitted records, failed caps, source service, sampling, grouping, verification and receiver work remain explicit. No prior validity judgment is a premise.

## Paired-comparison64-record range scan

For consecutive keys A,B,C, CMP A,B; CCMP B,C,#4,NE; B.EQ branches exactly when A=B or B=C. The slow path appends each equal adjacent pair separately with its own capacity check, in the same order as the original scanner. No relation or success event changes.

Three key registers rotate and four pointer-register pairs rotate. Initial copies put the range's first key and address in every possible suffix-start register. The pair remainder r=((n-1)//2) mod32 selects entry at pair32-r (r=0 selects the end test), using five TBZs and one added branch when r is odd. Each pair uses the previous key and record pointer, loads two new pointers and their keys, and retains the second key/address for the next pair. All potential entry registers initially hold the same first key/address. Since32 mod3=2 and32 mod4=0, the final key is x7 and final pointer x16, so MOV x5,x7 normalizes the next full block. The hot block costs32*(LDP+2LDR+CMP+CCMP+B.EQ)=192 instructions plus MOV,CMP,B.LO:195 per64 records.

Write n-1=64q+2r+epsilon with0<=r<=31 and epsilon in{0,1}. A nonmatching nonempty range costs24+6epsilon+(r mod2)+6r+195q including range control (the dispatch index is one `ubfx x17,x17,#1,#5`, exactly LSR#1 then AND#31). Subtracting195n/64 gives1341/64+189epsilon/64-3r/32+(r mod2). Its maximum is1588/64=24.8125 at epsilon=1,r=1, bounded by25. Empty ranges cost5. The slow path adds7 for only the first equality,8 for only the second,11 for both; an odd-prefix hit adds4. Charge8 per listed candidate plus one added8 for a full-list append attempt. Every append has its own capacity check and preserves original adjacent-pair ordering, even when both equalities hold. Global guard and RET cost2.

Thus a node uses at most195n/64+25*2^48+8*(2^24+1)+2 scanner instructions. This rational inequality holds per node separately, so summing yields195S/64 without an uncharged ceiling. Production S is divisible by64. The remainder suffix consumes exactly2r records after the optional odd comparison; every full block consumes64, and the exclusive end test prevents out-of-range loads. Empty ranges read no record, and no comparison crosses a range boundary. For identical input pointer/key sequences, scanner candidate identities, order and capacity-abort behavior are unchanged. This revision can change those input sequences; its enlarged set-valued load bounds, not an unchanged candidate-list assumption, justify success and cost.

The implementation is transferred from our secp package b00400fa-b4b3-440c-accb-657046f12efa with one change, the UBFX dispatch index above; the record-pointer/end-list ABI is identical. scanner_check.py exercises all remainder lengths, multiple blocks, empty ranges, both equality paths and caps, including4096-record cases. routing_check.py independently models every executed instruction and compares packet/candidate order across all4 nodes. Finite checks supplement the universal count and semantic argument above.

## 13. Grid, masks and provenance

SIDE=2^35,N=2^70,q2=2^28 and the35+29/35+1+28 metadata allocation are unchanged. The96bucket bits are[0,35),[64,99),[128,154); group144 and match182 include sixteen bits of the exceptional labels' top32constant bits. FILTER_FREE.md therefore enlarges the crude group/match coefficients rather than reusing a false uniform-bin assumption. FOURIER_CAPACITY.md stays the normal canonical-point lemma; the small extra fibers are bounded explicitly.

Immediate arithmetic/scanner source is our44137210-aa9d-47fa-8cd0-7daad206a9d1. The two itemized cold prices follow76ceba30-f1fa-4538-9c48-bc947cbd44bb. The removal of unallocated8*S padding follows37b8edde-addf-451f-9215-0dac6e7b6e79 and a fresh operation-by-operation audit of this changed receiver construction; public acceptance is not the justification. New work is filter-free set-valued-label accounting and the full-record receiver, with DONE-flag completion and fully revised cost/success/layout rows. The secp paired scanner and inherited constructions stay credited.

## 14. Complete current cost and storage binding

fixed_layout.py gives simultaneous nonoverlapping addresses for every source, label, tree/sort, batch, candidate/request/relay/home/flag/context/port/cursor region, and context pointers are initialized to those exact bases. d4 counters occupy[2^51,2^52); bulk hop1 slots, full-record buckets, pointers and the2^51-byte scan cursor copy follow2^52 and end at18216779176113864704<2^64. The sender header port and its local cursor table are separate allocations: the cursor-table base is a multiple of2^38 so the emitter forms a cursor address by inserting d1 into bits3..37, the hop2 d3 cursor/end table is a multiple of2^30 and the d4 counters start at2^51 so the resident receiver inserts d3 and d4 by BFI, and the512-record payload windows occupy a32768-byte stride per destination from a32768-byte-aligned base. The batch has a64-byte prefix reservation; four per-scan output mailbox regions and stop words follow the2^20 code/stack region with separate4096-byte child slots and isolated flag granules, and no flag is reinitialized after publication begins. Peak across2^70 nodes is below2^134bytes, conservatively declared below2^136.

The expected total is at most89623443111300696274737432003854996954054 instructions, log2(E/1000)=126.0752156631672870...; the rational logarithm enclosure supports outward126075215663167288/10^15.

The remaining named appendices are in proof-appendices.md and are part of this same claim.

## Short Montgomery multiplication and exact cost change

The new scalar multiplication core is adapted from Amazon s2n-bignum commit56152cdb51f47ef9f63b179cda0a8b7518379240, arm/p256/bignum_montmul_p256_alt.S. Copyright and Apache-2.0 OR ISC OR MIT-0 licensing are preserved in p256_fields.s and the repository LICENSE. Primary source and its corresponding arm/proofs/bignum_montmul_p256_alt.ml are available at https://github.com/awslabs/s2n-bignum/tree/56152cdb51f47ef9f63b179cda0a8b7518379240 . This attribution is not a claim that the upstream HOL proof automatically certifies our patched tails or register transformations.

The core first forms the full four-by-four schoolbook product with MUL/UMULH and explicit carry propagation. Its low limbs are reduced four times using the cancellation digit t and p's limbs from section3. Writing MC=2^64-2^32+1, the seven-instruction reduction adds t shifted32 to the next limb, carries t>>32 onward, and adds the low/high halves of t*MC to the remaining two limbs with carry. Each step is exact division of the current integer plus t*p by2^64. Four steps give raw T=(ab+m*p)/R0 with0<=m<R0; uniqueness of the cancellation m modulo R0 shows this T equals the previous multiply's raw value for all a,b<R0. The final high carry is consumed immediately by B.CC. The weak tail keeps T if T<R0 and returns T-p otherwise, always in[0,R0). Old and new correction predicates agree pathwise.

Including RET and standalone tail constants the MUL has110 instructions when correcting and104 otherwise; the generated loops remove the two tail materializations and hoist the MC materialization, so every executed MUL is3 shorter than the previous generated version on either branch (one forward and four backward MULs per incidence: 15 saved). Before unrolling this proves forward122 and backward616+4*carry from the previous125 and628+4*carry, not an average branch extrapolation. Defective rows complete as bounded garbage within these bounds (section3).

The RHS limbs occupy x7,x8,x13,x14 and are read-only after loads; MC occupies a distinct scratch register standalone and the reserved x2 in generated loops. The forward specialization maps the live canonical difference into that RHS bank and conflicting scratch into x24..x27; backward slope products substitute the cached inverse x24..x27 for all four RHS registers. The original fused square and its raw-label case split are retained; an alternate square offers no hot gain after MC hoisting.

## 15. Grouped early abort and expected work

Schedule. Each node's X indices form G=4 contiguous X-groups ending at floor(q1*e/64), e=21,36,50,64, run by all nodes in lockstep. Group g generates, samples and normalizes its X's, runs their batches and does hop1/hop2 as one transport phase of sections 2-4 and the appendices (there, a transport phase is an X-group and per-phase transport charges are per group; the four resolution phases run once per scan); the witness's group has at most qg=123857676 indices. Each node then scans all records received so far: prefix4c writes d4 range starts into the cursor copy, pass4 scatters, scan_ranges lists K2-adjacent pairs, and the resolution phases resolve the list with fresh caps (2^24 candidates, 2^29 requests, 2^34 relay records; lists, cursors and home flags reset per scan inside the2^37 per-stage allocation). The output tree uses the scan's own mailbox region. A root holding a verified x ([x]P=Q checked) latches and outputs it; after scans 1-3 it broadcasts stop or continue down the tree through per-scan stop words (own2048-byte granule, zeroed at setup, one parent STLR: 1 continue, 2 stop; awaited with the barrier LDAXR/CMP/B.EQ/WFE primitive; forwarded to at most four children; 1024 per node per scan). A node reading 2 halts; no node starts group g+1 before reading continue.

State. Only receivers write buckets, d3 cursor/end pairs and d4 counts; they persist, sized for the whole run. Source tables keep every X, Y and scalar at its original index all run, so any group's candidate is served identically. A scan only reads that state and writes its own scratch (cursor copy, pointers, list, phase-4 buffers, mailboxes), each rewritten before use. Scans start after the group's DONE signals and barrier; the next group's hop2 starts after its hop1 barrier, which each node reaches only after its scan and resolution barriers. Receivers are thus idle during scans, and each flag and stop word has one publisher and one use.

Success. Without a stop the run is the four-group form of sections 2-6 plus read-only scans: the final scan sees exactly the records, buckets and counts of the ungrouped four-phase run and fails only on the section 6 events, with qg in d_p. A stop needs a verified x, the answer since P has prime order r. Success is at least 0.3900000021509504 (q1=377471015).

Stopping time. Let T be the first group whose scan stops (5 if none). Group g and its scan run iff T>g-1, an event decided by earlier groups, and cost at most the deterministic bound W_g; linearity gives E[work]<=ONCE+sum_g P(T>g-1)W_g, no independence assumed (expectations over the uniform challenge and fresh coins, as for success). P(T>g-1) is at most the probability that the scan after group g-1 does not stop, at most PHI[g]=1-lo(b), b that group's end, lo the section 6 expression with q1=b (Gamma, mu, CS, intruder, pair loads, defects, meter at the prefix) and qg in d_p. Its bucket and meter caps are no larger than production and its list, request and relay caps are production, so every prefix drop term is conservative. PHI=(1,0.948168471,0.855210063,0.739555396), exact and nonincreasing in range_ledger.py.

Rows per node. X, per generated X: 2930/4*q2+39479+1024 arithmetic, 1623023 generation, 256 sampling, 2^17 driver, 1280 normalization, 2*q2*(2213/512+85/4) hops. GRP, per group: 8192*SIDE blocks, 11*SIDE+2 h2s, 14*SIDE+10 DONE, 3*SIDE+1 and 17*SIDE+3 cursor setup and finish, 64*1024 barriers, 2^14 control. SCAN(g): 595/64 per record so far (2*q2 per X), 11*2^26+1 pass4, (5/2)*2^48+1 prefix4c, 25*2^48+8*(2^24+1)+2 scanner, 4*2^41 routing, 2^24*20592111 resolution, 2^39 source service, 1024 stop word. ONCE: 2^48+1 zero4, 5*2^26+1 init3, 2^20 setup, 1577967*q2+1623023 Y and R_g, (q2+1)*256 sampling, 8192*q2 tree, 3968*q2+4352 sorting, 1 preset, 4*carry_cap (full meter cap, unweighted), 64 output; plus 2^28 global output. W_g=(group size)X+GRP+SCAN(g). With two groups, one scan and one barrier/control allowance these rows give the previous total 96224478515501901751686351617947303149568 exactly at q1=377471100.

Totals. E=N*(ONCE+sum_g PHI[g]W_g)+2^28 has ceiling 89623443111300696274737432003854996954054 instructions; the rational enclosure certifies log2(E/1000)<=126.075215663167288. Every run, all groups and scans included, does at most 2^126.232872 units (disclosed worst case). Preprocessing stays below 2^118 units. Memory adds the 2^51-byte cursor copy after the pointer array, and four mailbox regions and stop words in the fixed region; the node endpoint is 18216779176113864704<2^64.

Evidence. grouped_check.py runs the actual receive, zero4, init3, prefix4c, pass4 and scan_ranges kernels on four homes over four groups with full-bucket drops: a stop case halts after group 2 on a planted cross-group witness, a full case runs all four. Scans leave counts and bucket words unchanged, pointers and candidates equal the single-pass model of the records so far, and every count matches its formula. route.s with prefix4c assembles under clang (aarch64, armv8.1-a). The abort ledger adapts our secp design 266e82b9-a818-4077-a211-ef36749de8da to rescanning without hash chains; expected work to success is the usual collision-search cost (van Oorschot and Wiener, J. Cryptology 12, 1999).
