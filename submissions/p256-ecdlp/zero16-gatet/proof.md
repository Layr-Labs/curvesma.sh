The named .md documents cited below are included verbatim as named appendices in proof.md and proof-appendices.md, not separate missing files; sections16 (correction tails) and17 (run gate) are proof-tails.md.

# P256 range-specialized stored-leaf search

Proposed complete analytical claim for exact-version external review. This is not an accepted or formally certified result.

NORMALIZATION.md, FRACTION_LABEL.md, RECEIVER_PROTOCOL.md, COLD_ACCOUNTING.md, RANGE_DRIVER.md, RESOLUTION_PROTOCOL.md, BLOCK_PROTOCOL.md, section15 (grouped early abort, hash chains, expected-work ledger) and fixed_layout.py are operative parts of this argument.

## 1. Target, representation and candidate claim

The target is exactly the protected P-256 group, including a=-3, its specified
prime p, prime order r, cofactor one, generator P and uniform Q=[x]P. Output is
the canonical scalar x. Q=infinity returns zero. No auxiliary powers or
secp256k1 endomorphism are used. Let R=2^256 mod p. Store canonical generator Montgomery
coordinates; hot intermediate values may be redundant in[0,2^256). For a nonidentity point Z let T(Z)=floor(2^192*x(Z)/p); x is an
exact invariant of the two-element class {Z,-Z}: the curve equation gives at most
two points for each x, and no nonidentity point has y=0 because r is odd.
A record's key is bits16..191 of the computed fraction (176 bits); the kernel computes T within[-19,10] (section3), so equal classes give equal keys except on a charged label-error set, and FRACTION_LABEL.md proves the point-only cell bounds used for all capacity and pair-load bounds.

The complete claim is expected total classical
work below 2^125.025952461019009 units at the protected C=1000 scalar instructions/unit
(the expectation over the challenge and the algorithm coins of the section17 gate and, if open, the grouped early-abort run of
section15, under the disclosed correction-tail heuristic H_T of section16; on every run the work is at most the disclosed worst case 2^180.6801 units, all526 groups run
and every uncapped hash walk visiting its whole chain), success at least 0.39, peak memory below
2^136 bytes, preprocessing below 2^119 units included in the total, and zero nonuniform advice.

Gate (section17): a fresh 64-bit coin u opens the run iff u<K=7768471761051508417. Success>=(K/2^64)*s_run>=0.39 with s_run>=0.92608 (section6); E=64+(K/2^64)*E_run with E_run the section15 total.

## 2. Random construction and pair arithmetic

There are N=2^70 nodes, indexed (i,j) with each coordinate in [0,2^35).
At each node draw b and q1=893476688 independent a_k uniformly modulo r,
and q2=2^28 independent c_l uniformly in {1,...,r-1}. Coins and draws at
different nodes are independent. Each scalar draw tries at most eight fresh
256-bit values, rejecting out-of-range values and zero where required. A cap
failure aborts the run. The rejection probability is below 2^-31, so the
union cap loss is at most N(q1+q2+1)*2^-248. Sampling is charged 256 scalar
instructions per requested scalar. Under the protected ideal-coin convention,
RNDR supplies the bits; a hardware-failure assumption is not silently added.

Initially construct R_g=[b]Q, X_k=R_g+[a_k]P and Y_l=[c_l]P. These original points and indices are used in the probability argument below. Apply the exact signed normalization and paired scalar updates in NORMALIZATION.md: store X'_k=epsilon_X X_k with a'_k and row-specific b'_k, and Y'_l=epsilon_Y Y_l with c'_l. Build the y-partition/x-sorted range tree of72-byte records with original Y IDs; retain canonical source tables at original indices. Keep original X index k and group membership (section15). Preserve original b while constructing X, separately from stored row-specific b'.

The executed leaves are X'_k+Y'_l and X'_k-Y'_l. Their two negation classes per original incidence are exactly the original two classes, up to sign permutation. Write n0=2q1q2 and S=N*n0. Executed coefficients are m=a'_k+s*c'_l and n=b'_k. In the next arithmetic paragraph X,Y denote normalized stored points; in the probability analysis they denote the original independent generators under this coupling. Identity X is stored as the in-band pair(0,0), impossible for a nonidentity point of this odd-order curve, and detected before each row by four paired loads, seven ORs and a test/branch (at most20instructions). Construction writes this sentinel from the cold infinity flag within the per-item allowance; row testing is inside the driver fixed reserve. Identity X or equal x coordinates marks the batch defective before inversion. Normalization and sorting do not change these failure events.

For each normalized X, use the disjoint at-most87 descriptors in RANGE_DRIVER.md. They fix the x-difference direction, which selects the forward variant; one backward kernel serves all. Forward calls traverse them in one sequence, preserving the prefix across calls; one inverse serves the whole row. Reverse calls traverse them backwards. Write d_l for the absolute integer difference of the canonical Montgomery x representations. Modulo p this represents either sign of x(Y_l)-x(X), with the Montgomery factor. All later equations are field equations with Montgomery scaling understood. The forward loop stores d_l
and prefix product A_l=product_(j<=l) d_j, with A_-1=1. Invert A_(q2-1)
once. The backward loop maintains B_l=1/A_l, computes i_l=B_l*A_(l-1), then
updates B_(l-1)=B_l*d_l. Thus i_l=1/d_l. Then w_l=i_l^2 and M1_l=y(Y_l)*w_l, and both leaves'
labels follow from x(X-+Y)=(y_Y^2+y_X^2)w-+2y_X(y_Yw)-x(X)-x(Y) by two truncated products against tabled fraction rows (section3).
No a=0 identity is used.

Hot field intermediates are exact representatives below R0=2^256; labels are fractions of true residues (section3), so no representation filter is applied. The assembly forward
loop multiplies every d_l with no zero test; defective rows complete as bounded garbage. The backward loop
writes no batch entry; it emits two records per incidence.
A 32-byte Montgomery-one prefix before the batch handles l=0 without reading
outside the allocation. The reference point formulas and executed native-word
interpreter tests check these invariants independently; they do not replace them.

## 3. Field arithmetic, fraction labels and worst-path counts

Let R0=2^256 and delta=R0-p. Original generator coordinates are canonical Montgomery residues. Only M/S receive weak operands belowR0. Canonical cold arithmetic follows COLD_ACCOUNTING.md.

For radix B=2^64, p has limbs(B-1,2^32-1,0,B-2^32+1). In each reduction step the current low word t is the cancellation digit, and(T+t*p)/B is an integer. Split t=t_low+2^32*t_high: the corrections use t*2^32 and t*(B-2^32+1), with explicit carries and borrows. Four schoolbook rows and these reductions compute raw REDC exactly. The symmetric square uses each diagonal once and each cross-product twice. For a,b<R0, raw T=(ab+mp)/R0,0<=m<R0, satisfies0<=T<R0+p<2R0. Its high limb is0 or1. Keeping T when high=0 and subtracting p when high=1 gives an exact representative in[0,R0). The carry-only MUL and SQR tails implement precisely this rule, with at most110 and88 instructions including RET for the standalone MUL and SQR.

The forward absolute difference subtracts the smaller canonical integer from the larger in the descriptor's order, producing the nonnegative magnitude d_l<p without a per-item negation; it represents one of the two field signs of x(Y)-x(X), and the sign disappears in w=i^2. No zero test remains; defective rows complete as bounded garbage within the same envelopes. Direction specialization is machine-checked limb by limb: pair_fwd_hi computes Y-X and is selected exactly when Y.x>=X.x; pair_fwd_lo swaps all four limbs (chain_witness.py normalize_transform) to compute X-Y and is selected exactly when Y.x<X.x. range_witness.assert_limb_specialization asserts the exact limb sets on the rolled base program, and the probe executes pair_fwd_lo on ten boundary operand pairs and asserts the stored denominator equals the exact absolute difference.

The fixed inverse computes exponent p-2 with255 squares and127 multiplies. Induction on the preceding REDC range argument proves all intermediates stay belowR0 and congruent to the specified modular powers. The fixed inverse instruction sequence is retained byte-for-byte after macro expansion; its conservative39479 allowance is retained without crediting the smaller counts. On non-defective rows its input is a nonzero prefix product, so Fermat's theorem gives the field inverse.

**Fraction labels.** For a nonidentity point Z let T(Z)=floor(2^192*x(Z)/p), a 192-bit integer; the record key of Z is bits16..191 of T (176 bits). For constants c let row_k(c)=floor(2^256*(c*2^(64k) mod p)/p), k=0..3, four words each, and for an integer multiplier m=sum m_k 2^(64k)<2^256 let TP(m;rows)=sum_k floor(m_k*row_k/2^64) mod 2^192: the words of each 320-bit product m_k*row_k at2^64..2^255, with its lowest word dropped. Since sum_k m_k*(c*2^(64k) mod p) is congruent to c*m, frac(sum_k m_k*row_k/2^256) differs from (c*m mod p)/p only by the row floors: 2^192*(c*m mod p)/p-TP lies in[0,4+4e) modulo2^192 when every row is low by less than e units of2^-256 (four limb floors below one unit each, and m_k<2^64). In Montgomery form w=i^2*R and M1=y(Y)*w*R as stored; the per-Y entry (256 bytes) holds row_k(y(Y)^2*R^-1) at+32k, T(x(Y)) at+128, a zero word at+152 and the stored y(Y) at+160, and the per-X area (512 bytes) row_k(2*y(X)*R^-1) at+32k, row_k(y(X)^2*R^-1) at+128+32k and T(x(X)) at+256, so that the kernel reads both tables in whole64-byte lines by LD64B: two rows per line, and the Y entry's third line holds T(x(Y)), the zero word and y(Y) (section4). The entry's y(Y) is the normalized stored y of Y'_l, the same value (same epsilon_Y, NORMALIZATION.md) that the Y records and the source tables hold, so M1 and the sign of each emitted record are exactly those of 008389f5. With a fresh row sum sigma_k=row_k(y_Y^2R^-1)+row_k(y_X^2R^-1) mod2^256 (low by less than2 units), TP1=TP(w;sigma) is low by[0,12) and TP2=TP(M1;row(2y_XR^-1)) by[0,8) against the fractions of (y_Y^2+y_X^2)w and2y_Xy_Yw; T(x(X)), T(x(Y)) are low by[0,1). As x(X-+Y)=(y_Y^2+y_X^2)w-+2y_Xy_Yw-x(X)-x(Y) in the field (i=1/d, any X,Y with x(X)!=x(Y); no a=0 identity), L+-=TP1-+TP2-T(x(X))-T(x(Y)) mod2^192 satisfy L-T(X-+Y) in[-19,10] (the four deficits and one final floor). Hence the computed key of a leaf equals floor(T/2^16) of its point unless T mod2^16 is among the29 residues[0,19) and[2^16-10,2^16): the label-error row of section6. The label is a function of the true residue (no filter). FRACTION_LABEL.md proves the cell bounds used in section6. frac_witness.py is a bit-exact model of the kernel; probe.py and raw_differential.py compare every emitted record against it and every label against the exact T of the true sum.

**Backward kernel pair_bwd_f.** One backward kernel serves every descriptor. Per incidence it runs the retained batch-inversion pair i_l=M(B_l,A_(l-1)), B_(l-1)=M(B_l,d_l) (with B stored as before), then: w=M(i,i) (75 plus the B.CC tail test); the ID word by one pre-indexed LDR that also steps the Y-record pointer by -72; the Y entry address FB+(ID>>28) by one ADD (x30 holds the entry base FB through the loop, restored after every BL by the existing x16 save; entries are256 bytes at original index l); TP1 with the16-word row sum (72: the entry's first two64-byte lines by LD64B into x10..x17 with one ADD of x18 between them, eight LDP of the X area's rows into x3,x4 added in place (LDP leaves the flags), sixteen sum instructions, twenty-four multiplies and twenty-one accumulations); U=TP1-T(x(Y))-T(x(X)) (an ADD and one LD64B of the entry's third line, which brings T(x(Y)), the zero word and y(Y) into x10..x17, an LDP and an LDR of T(x(X)), six instructions); M1=M(y(Y),w) (95 plus the tail test: the base stream with y(Y) read from x14..x17, its temporaries x15,x16 renamed x11,x12 and its second output limb x27); TP2 (48: the X area's2y(X)R^-1 rows by LD64B x10,[x22], one ADD and a second LD64B, then the same forty-five multiply and accumulate instructions); L+=U-TP2 (3) and its record (9); L-=U+TP2 (3) and its record (8). A record is[K2,K3,W2,ID'] with K2=L1, K3=EXTR(L3,L2,#35) (L2 bits35..63, then L3 bits0..34), ID'=ID|sign<<35|L3>>29 (d2, one ORR of a shifted register) and d1=L2 mod2^35 (BFI into the cursor address); W2 and the stores, the LDADD slot cursor and the TBZ full-window test are unchanged. The key is therefore K2, the routing digits d1, d2 and K3 bits0..28 and45..63 (bits29..44 of K3 are fraction bits0..15, below the key, and the receiver masks them). A tail-free call of n>=1 incidences costs exactly529n+21-14*floor(n/8), plus4 for each taken tail among its four sites per incidence (i, B, w, M1), checked call by call against the model; with all tails at most2173/4*n+34; an empty call costs15. Forward entries cost at most441/4*n+46 (unchanged). Eighty fixed instructions per segment (46 forward,34 backward) over at most87 segments (6960) are covered by the2^17 per-row driver charge; no separate term is charged.

The generated loops reserve x0 and x2 for the zero-extended constants0xffffffff andMC=0xffffffff00000001, initialized once per nonempty loop. MOV does not change flags; for OLD=0x00000000fffffffe=~MC modulo2^64, ADC(y,OLD,C) and SBC(y,MC,C) have identical result and NZCV, and conversely SBC(y,OLD,C)=ADC(y,MC,C); every original CONST_B consumer is transformed this way. The structural witness checks the removed-MOV counts and rejects any generated loop instruction writing either reserved constant register.

Full flushes call the shared flush_publish_b subroutine (16 taken /13 dropped, plus3 call-site), charged19 per firing in the hop1 row (291 per128 leaves:268 LS64 delivery with ACK and in-flight bracketing, the19-instruction flush,4 reserve); L2 is placed in x7, the subroutine's d1 word, and every value live across the call (U, TP2, ID, x13, x14, x17, x30 via x16) is outside the subroutine's written set. The worst-case arithmetic/emission price is2614/4 per incidence (2534/4 with all five correction tails skipped, section16). Backward operand reuse in the batch-inversion pair is retained (the first multiply writes x1,x17,x28 in place of x3,x4,x5, so the second keeps x3..x6 and drops two reloads). The cursor base, W2 word, ID word and LDADD addend are carried from the first record to the second. All forms assemble under clang as single AArch64 instructions (ADD and ORR with an LSR-shifted register, EXTR, TST with the rotated-run mask; probe.py's ISA block executes each form against Python integers). A dedicated full-window test prefills every emitted window to127 records so all280 emissions across its two trials fire the real subroutine (140 taken with exact publish metadata and cursor advance,140 dropped with cursor reset), with exact per-call counts529n+21-14floor(n/8)+4tails+19F+16D. These constants enter section15; range_ledger.py verifies the outward rational score125025952461019009/10^15.

frac_witness.py derives pair_bwd_f from the fused y0_c0 variant of the hash-pinned rolled base (range_witness.fused_base with the operand-reuse pass): it keeps the prologue, the batch-inversion pair and its B store, replaces the tail after the B store by the code above (the square and multiply bodies are the base's own instruction streams with renamed registers), moves the frame set-up before the empty-call test, and range_witness.unroll8 makes eight uniquely labelled copies with a biased unsigned count (subtract8, skip to the tail on borrow, repeat on no borrow, add8 back for the remainder). The finite replay supports, but does not replace, the control-flow induction.

The generated arithmetic listing has28048 encoded bytes, checked by parsing its scalar instructions; all supplied .s helpers are macro-expanded and counted by the same parser, and their encoded bytes plus28048 generated,2688 unrolled-block and65536 stack bytes are asserted within the1MiB code_stacks allocation. The program is a public constant; per-node initialization materializes each32-bit instruction with two immediate moves, one store and one pointer increment, within the1MiB code initialization charge, with no runtime generator, lookup or dispatch.

## 4. Fully charged transport and grouping

Use the explicitly specified two-dimensional blocking ports of RECEIVER_PROTOCOL.md. Hop1 uses the explicitly fixed128-record window, native266-instruction LS64 delivery/ACK and transaction-bound completion protocol of BLOCK_PROTOCOL.md. Its full-flush work is charged at291 per block (268 delivery with ACK and in-flight bracketing, the 19-instruction fused flush, 4 reserve); all partial/initialization/finishing work is charged per group in the worst-case group_fixed row. Ordinary final-cursor, barrier and resolution messages retain their addressed delivery cost. Hop2 uses a separately charged resident scalar receiver context at each node: it appends each accepted record to the node's single home record array and immediately joins it into the hash chains (12 instructions per accepted record plus4 per chain visit and5 per K2 match;5 for a full-array drop), with no per-message acknowledgment or in-flight updates. Each sender closes its group stream with one DONE marker per destination; the receiver counts markers and signals per-group completion with a single STLR to its DONE cell, which the main context awaits with the same LDAXR/WFE completion primitive as barrier epoch words. No direct-home shortcut, uncharged key operation, software polling or unbounded queue is assumed.

**LS64 (machine premise).** Each node implements FEAT_LS64 (SCTLR_EL1.EnALS=1): LD64B/ST64B Xt,[Xn] move Xt..Xt+7 in one single-copy-atomic64-byte access; t even, at most22; no offset or writeback; Xn64-byte aligned (Arm DDI0487). Supported locations are IMPLEMENTATION DEFINED: here exactly the hop1 payload windows as mapped into the delivery context (Device; LD64B), the hop1 slot region (Normal Non-cacheable; LD64B, ST64B), the hop2 row-payload port (Device; ST64B), the receiver packet window (Device; LD64B) and the fraction tables, that is the per-Y fraction entries and the per-X fraction area (Normal Non-cacheable; LD64B; written by ordinary stores before they are read: the Y entries once during preprocessing, the X area at the start of each X row), nowhere else; a64.py asserts every rule on each executed LS64 instruction and probe.py checks that violations are rejected. Each counts one step like every instruction. The slot region starts at2^52+32, so block_dlv's63 slot pairs align on both sides; h2s and h2done send one ST64B packet per record or DONE marker, read by one LD64B (layout and the bit63 DONE test: RECEIVER_PROTOCOL.md). pair_bwd_f reads each fraction entry by three LD64B (x18 at the entry, then two ADDs of64) and the X area's2y(X) rows by two (x22, then x22+64); both tables are64-byte aligned (fixed_layout.py). Without FEAT_LS64 the fraction tables would be read by LDP as in008389f5 (12 more per incidence: the claim would be0.025140161631570 bits higher, exactly) and the transport would revert to the LDP/STP moves of ba0f040b-7f86-43e7-9fc2-c72543fb3a1b (about0.028 bits more, as measured on that line).

Key digits are d1=L2 mod2^35 and d2=L3>>29 (fraction bits[64,99) and[29,64)); no d3 bucket digit is used. The record is(K2,K3,w2,w3) with K2=L1 and K3=EXTR(L3,L2,#35) (section3), w2=i|((k mod2^29)<<35) (section17), initially w3=d2|(sign<<35)|(l<<36). After extracting d2, h2s replaces only w3bits0..34 with source column j. Thus a home reads source(i,j) directly; sign/l and k mod2^29 are preserved. l<2^28, including maximum-boundary and post-emission-wrap cases.

Partition each node's X indices into G=526 contiguous groups of at most qg indices (range_ledger.py group_sizes; section15); do not renumber k. Every node runs the same group at the same time, and the group loop tests the latched stop flag before each group (section15). For each group initialize fixed payload windows and local/remote cursors, emit ordered128-record blocks into sender-row slots at(d1,j), flush each remaining partial block with initialized padding, send clamped valid final cursors, drain and barrier. Padding is excluded by those cursors. Each per-group slot has CAP_PG=5*2^17 records; the slot region is reused by the next group only after the barrier that follows this group's hop2 drain. Stream all its stored records to(d1,d2); no hop2 pair slots exist. The destination receiver compares its register-resident record cursor with the last admissible position and drops the record if the array is full; otherwise it stores the record at the next32-byte position, prepends the record to the hash chain of its K2 index and walks that chain (section12). After streaming, each sender publishes one DONE marker per destination; each main context awaits its receiver's DONE signal, then all nodes barrier. The group's new candidates (the list entries between the previous and the newly published list end) are then resolved by the four phase-4 stages of RESOLUTION_PROTOCOL.md; a resolver-confirmed x latches the stop flag.

The home record array and its cursor, the hash heads, next words and the candidate list persist across all groups. Cap CAND_CAP=2^24 is checked at every attempted append.

A listed pair agrees on the 70 home bits, all 64 K2 bits and the 42 further key bits of the masked K3
(one K3' compare in a chain of equal K2 bits0..57): the 176 key bits. The hash-index bits are within K2. Partial matches are
not treated as true group collisions. Each is resolved and verified.

Resolution uses the four explicit bounded64-byte-record phases of RESOLUTION_PROTOCOL.md, not count/offset exchanges. Source and home tags use35-bit node coordinates; X index29bits, Y index28bits and sign1bit. Source tables keep original indices and row-specific signed b. The protocol includes initialization, caps, partial-reply flags, allocation, source gather and a final bounded fan-in output tree. Its complete work is charged by four N*2^41 rows per group and source service1024 per capped request; the candidate list is cumulative with its2^24 cap, so over all groups at most2^24 candidates are resolved and2^25 requests served, inside the single resolution and source-service rows. Source cap2^29 and relay allocation2^34 records are fixed; the original stricter logical relay loss bound is retained as explained there. A home reconstructs X+sY and X'+s'Y', tests equality up to sign, solves the nonzero scalar denominator and verifies [x]P=Q. Scalar arithmetic uses the complete COLD_ACCOUNTING.md bound42410496 per candidate.

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
For a good witness, whose class fraction T avoids the29 label-error residues modulo2^16 (section6), both witness records carry the key floor(T/2^16); records of other classes with that key are the charged intruders and false candidates. Under the separately charged no-intruder event, this class forms the relevant contiguous group. Only for such a good witness, if its key class has any useful pair, it has a useful adjacent pair: writing its
members u_i Z, uselessness is equality of n_i/u_i, an equivalence relation.
If every adjacent pair were useless, all members would have the same ratio.

## 6. Defects, capacity losses and success expression

Every computed key of a point Z lies in the set {floor((T(Z)+e)/2^16): e in[-19,10]}, at most two consecutive values (section3); FRACTION_LABEL.md gives the point-only cell bounds. There is no carry meter: the backward kernel forms no x-sum. A batch defect probability is at most
d_batch=q2*2/(r-1)+1/r (zero denominators or identity X). Excluding any
defective batch in the whole run costs at most N*q1*d_batch. The event is
explicitly charged even though a defective batch's cost is no larger.

For the35-bit slot pattern (d1, fraction bits[64,99)), KB(35)=(513/512)*2^221/(r-1), and KBH=(513/512)*2^186/(r-1) for the70-bit home pattern (bits[29,99)), by the interval lemma of FRACTION_LABEL.md. The pair-load term h uses the crude home coefficient2+2^-23 (KB70); the176-bit key uses KB(176)=(2+58*2^-16)*2^80/(r-1); fixed cells whose lowest fixed bit is at least29 use2+2^-23 (KBC). These concern the fixed set of keys a point can receive, not the selected value. Let rho=(r/(r-1))^4.

Condition on the witness's four generators, but never on global cleanliness or actual raw-label choices. Retain the reservation12 fixed leaves. Its hop1 group slot has qg=max group size, mean m_p=2*qg*q2*KB35+1 and variance bound v_p=m_p+1+60*(qg+q2)*KB35. Its home record array spans all groups, mean m_h=S*KBH+1 and variance v_h=m_h+1+60*(q1+q2)*KBH. Then d_p=v_p/(CAP_PG-12-m_p)^2 and d_h=v_h/(CAP_H-12-m_h)^2, CAP_H=ceil(m_h)+2^50. The witness's group is predetermined by its X index k, so the uniform bound over the equal-sized groups applies and no group multiplier is needed. There is only one pair-slot loss, because no second-hop slot exists.

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

For pair loads define M=2*KB(176), since a fixed opposite point offers at most two keys. The source's shared-generator conditioning cases apply to these fixed set-valued relations, with
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

The intruder rate is rho*(S+12)*KB(176): a hash walk stops only at a prior record
agreeing with the new record on all 176 key bits (home, K2 and the masked K3), so only a 176-bit-agreeing
record of another class arriving at the home between the two witness records can split
their class's chain (section15); every other record is merely visited. If no such intruder
arrives and the capacity events do not occur, the useful chain edge is listed and resolved. False partial-key matches can exhaust caps; their
work and losses are exactly what the pair-load bounds cover.

The explicit lower-bound expression is
1-exp(-mu)-CS
-mu*(1+1/r)*(1/r+rho*(S+12)*KB(176)
             +2*(d_batch+d_p+d_h+d_cl+d_rq+d_rl))
-N*q1*d_batch-N*(q1+q2+1)*2^-248-Gamma*4*29*2^240/(r-1)^2.
Here Gamma is the cross-node unordered-pair count. For a useful true collision the endpoints share x, hence T. At most29*2^176 values of T are bad (T mod2^16 among the29 label-error residues), each taken by at most2^64 values of x, each x by at most two points with at most two negation partners, so a pair is both a collision and bad with probability at most4*29*2^240/(r-1)^2. Union-charge this before any capacity or survival conditioning. Actual records and partial matches are dominated by the point-only key cells of FRACTION_LABEL.md. Removing records can create adjacent matches, so no adjacency-list subset claim is used. Wrong labels can produce false candidates, all fully charged and verified.
The factor 1+1/r bounds (2/r)/P1. All quantities except exp are exact rational
functions. Round mu down to 200 binary fractional bits and upper-bound its
negative exponential by the even 40-term Taylor polynomial, an upper
bound for every mu>0 by the sign of the Taylor remainder (the ledger asserts mu<4). The resulting rational success
bound exceeds 0.9260804969152573 at q1=893476688 (the largest q1 whose home records end below the heads at7*2^61; section17 identifies rows k>=2^29). This is the success of one run; the section17 gate multiplies it by K/2^64, and the declared overall success is0.39. None of this identifies reviewer confidence with algorithm success.

## 7. Work and memory

range_ledger.py contains the complete exact rows: the E-averaged per-group rows of section15 (arithmetic, emission, per-row driver, both hops including the receiver's append and hash insert of every delivered record, and the expected main-pair walks), the global same-row walk, K2-verification, defective-walk, correction_tails and first_prefix_tail (section16) rows, and the worst-case rows for all failed caps, generation, capped scalar sampling, normalization/sorting, the fraction tables, preprocessing, per-group fixed work, phase4, resolution, verification and output. No unallocated per-leaf bookkeeping reserve is charged: all record operations are itemized below and all driver/context work is charged per batch or group. Work is summed across both scalar contexts at every node.

Cold generation/header work is charged23441408instructions per requested generator item, and candidate resolution42410496instructions per capped candidate, including failed verification. Each X row's fraction area costs at most2^17 and each Y fraction entry at most2^16 (fraction_tables row, COLD_ACCOUNTING.md). These are the explicit itemizations in COLD_ACCOUNTING.md; the shared-per-node generator and exact scalar-inverse refinement follows the itemization in bdf01aa5-16e3-4619-a8eb-1aa0f6b6b29b (its acceptance is not a premise). Preprocessing (node initialization, generation, scalar sampling, the one-time hash-head zeroing and the range tree) stays below2^119units and is included; range_ledger.py asserts that bound on a superset of the setup rows.

Reserve the first2^52bytes of each node for code, both contexts/stacks, packet/ack state, DONE cell/granule, the receiver's list-end cell, the unused former d4-counter range, scalar/X/Y tables, the Y fraction entries (256 bytes per original Y index) and the512-byte X fraction area, batch storage, source/relay/candidate buffers and control tables. The fixed-region calculation in range_ledger.py includes the DONE granule and stays below that reservation. Place the per-group hop1 slot region next,32*SIDE*CAP_PGbytes from2^52+32 (slot bases32 mod64); then the home record array of32*CAP_H bytes (section12), ending below7*2^61; then the2^58 hash head words (2^61bytes) at[7*2^61,2^64). All regions are nonoverlapping and eight-byte aligned; the heads end exactly at2^64 (last byte2^64-1), and the attempted hop1 cursor bound cannot wrap. The receiver tests its cursor against the last admissible position before any store, so no temporary overshoot exists. This is a simultaneous-allocation bound, not an assumed lifetime reuse saving. AcrossN=2^70nodes total storage is at most2^134, hence the declared2^136bound. No nonuniform advice is used.

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
full-block flush, delivery and in-flight counts at291*S/128 in total: ordinary per-incidence sender emission is already included. BLOCK_PROTOCOL.md proves the block bound; per-group initialization, partial flush, finishing, endpoint registration and control (1187*SIDE+6) sit in the worst-case group_fixed row.
The dedicated receiver costs12 per accepted record, including the LD64B packet read, the append, the BFI/SWP hash insert and the chain-head test (window release as in RECEIVER_PROTOCOL.md), and5 per full-array drop, with no per-message acknowledgment or in-flight updates; charge12 for every potential record. DONE-marker publication costs2+4*SIDE per sender per group, DONE handling costs6 per non-final and9 per final marker (the final one also publishes the list end), and the DONE signal/wait/reset and the main context's list-end load are charged per node per group inside the group_fixed row; hop2 totals481/32 per potential record (h2s129/32 plus11).

Per-group work (itemized in BLOCK_PROTOCOL.md) is at most1236*SIDE+27, charged at that bound per group per node for all526 groups in the worst-case group_fixed row. The receiver's resident registers are set once inside the2^20 phase/context-control allowance; the2^58 hash head words (at7*2^61, 256-byte aligned) are zeroed once by zero16384 at1026 per16384 words plus its RET: 512 DC ZVA zeroes with explicit addresses (ADD #256) plus the loop. The machine's DC ZVA block is DCZID_EL0.BS=6, 256 bytes (IMPLEMENTATION DEFINED; as on Fujitsu A64FX), with DZP=0 and SCTLR_EL1.DZE=1, on cacheable Normal memory, never a port; each ZVA zeroes its aligned256-byte block and counts as one instruction like the other single-instruction memory operations (a64.py asserts the alignment). 2104 barriers (4 per group) at1024instructions/node and a2^20phase/context-control allowance per node (33 per kernel call, at most20 calls per group, the stop test and group loop) are charged. Per-batch setup and group/context control stay separately bounded.

The remaining rows charge all generator points and headers, eight-attempt scalar sampling, program delivery and node initialization, all capped resolution (42410496 per candidate, COLD_ACCOUNTING.md, even for false partial matches and failed verification), source request service (at most1024 per request including addressed copies), the four bounded-routing stages per group, barriers, target validation and output. All maxima are charged whether or not their paths occur.

For final scalar arithmetic use the exact construction in COLD_ACCOUNTING.md: a256-step shift/add product costs36928 instructions including loop/setup; a512-product inversion with32 controller instructions per product and512fixed costs18924032. Five-word carry-aware addition/reduction preserves canonical residues. This is the operative bound, not a rounded product bound multiplied without controller cost. scalar_check.py and cold_accounting_check.py supply independent finite arithmetic checks.

All hot assembly uses base scalar GPR operations, plus the explicitly
specified LSE atomics for routing and the FEAT_LS64 moves of section4. Constants in MOV and logical masks are
single AArch64 encodings: zero/one/powers of two, contiguous bit masks of the displayed widths,
and the cyclic runs of bits in 0xffffffff00000001 and in the receiver's K3 mask0xffffe0001fffffff (clang-checked). Branch displacements are
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
AArch64 interpreter in the protected runner. Tests compare against Python integers and affine formulas.
The two-hop test executes deliveries, data-path kernels and the receiver's hash join over two executed groups of three planned (a toy stop rule latches after the second) on a 4-by-4 machine against separate references, with the exact receiver path model and full path coverage; its planted keys isolate transport and grouping only.
None of these finite tests proves a production success rate or physical feasibility.

Provenance, not premises (the argument is reproduced here): e0d1ed74, 80e9e8af, 1efcef88, 3828873f, 472226d8; block/header concept after secp e59be87d. Target-dependent quantities are recalculated; no previous judgment is reused.

## 11. Register and address transformation

chain_witness.py and range_witness.py deterministically define the generated range program from the original call-based schedule and the supplied field/fused bodies, specializing memory addresses before bijective register renaming. Stable pointers/counters are x19..x23; the forward loop loads X once into x7-x10 (never written by the renamed multiply) and reads the prefix relative to x21 (pinned tail) or the ping-pong banks (mains), so x22/x23 and their setup MOVs are gone there. The first backward per-item inverse is produced in x24..x27 and stays live across the prefix update; the square reads that bank. x30 holds the Y-entry base through the backward loop; LR is saved on entry and restored on return, and around the flush call x30 is held in x16. The backward Y pointer is stepped by the pre-indexed ID load ([x20,#-72]!), with y(Y) read at offsets-32 and-16. No caller-visible scratch or source point is removed.

Chained results substitute producer output registers for the next kernel's left-operand loads; scratch32..63 stores and loads are dead after inverse caching and absent. The forward direction-specialized difference stays in x11..x14 with no zero test (defective rows complete as garbage); multiplication temporaries avoid that bank and its two right-operand loads become direct reads. The initial Y load is postindexed by72, its upper limb load uses offset-56 and the trailing pointer ADD is omitted. Labels are uniquely renamed; address substitution/renaming adds no flag-setting instruction, so arithmetic branch predicates are preserved. The witness checks all structural preconditions and exact generated source equality; random-register tests compare every emitted record with the bit-exact model of frac_witness.py and every label with the exact T of the true sum, plus source memory, SP and LR.

## 12. Complete ledger and finishing deliveries

The expected total is43307381650312544483650326063708894449647instructions (integer ceiling of the gated E of section17, under H_T of section16), with log2(E/1000)=125.02595246101900821... . A64-term rational atanh enclosure and independent ln2 interval certify the outward claim125025952461019009/10^15 as the least multiple of10^-15 above the enclosure. No prior validity judgment is a premise. Fixed block windows, hop1 transaction-bound ACK events, hop2 DONE-flag completion and the latched stop flag are explicit machine assumptions; a flag-polling or general remote-memory implementation is not covered.

## Receiver hash join and in-record next words

The receiver (route.s receive) joins every accepted record into its home's hash chains at arrival; there is no separate join pass. Each home owns2^58 head words at[7*2^61,2^64): the head of index K2 mod2^58 is at7*2^61+8*index, formed by one BFI of K2's low58 bits into bits3..60 of the resident register x12 (bits61..63 set, bits0..2 clear; no other instruction writes x12), zeroed once by zero16384 before the first group, and one home record array of32-byte records below7*2^61 with a register-resident cursor and end. A record pointer p addresses the record's w2 word: its next word at p-16, K3' at p-8, (w2,w3) at p; heads, next words and list entries hold such pointers (zero ends a chain), and make_requests reads (w2,w3) at p. SWP (FEAT_LSE, like LDADD) is one load-and-store; after the zeroing only the receiver writes its heads and records. Records are appended in arrival order, so every accepted record is stored and inserted exactly once; records dropped at the full array are neither stored nor inserted.

Accepted path (12 including RET): the LD64B of the whole packet, which releases the window, and the DONE-bit branch (2); CMP/B.HS of the cursor against the last admissible position (2); a pre-indexed STP of (w2,w3) advancing the cursor by32 (1); the BFI head address (1); LSR and BFI forming K3', the packet's K3 with bits29..44 (fraction bits below the key) replaced by K2's bits58..63 and zeros (2); an SWP that stores the new pointer into the head and returns the old head (1); an STP of (old head,K3') (1); CBZ on an empty chain and RET (2). Every record of a chain has the same K2 bits0..57 (its head index), so K3' equality is exactly agreement on K2 and the masked K3, i.e. on all176 key bits. Per visit: one LDP of the prior's (next,K3'), CMP, B.EQ and a CBZ/CBNZ on the next pointer (4); the walk alternates two pointer registers. On a K3' match: CMP of the list cursor against the list end and B.HS back into the walk (2); with room the post-indexed STP of (prior,new) replaces the visit's pointer test, then RET. A full-array drop costs5. The receive row is12 per potential record; walks are charged4 per visit and, conservatively,5 per K2 match (section15; a full-key match costs2).

Same-index records form one arrival-order chain, newest first, every member reachable from the head: the new record is prepended (its next word is the old head) before the walk, which starts at the old head. The walk is uncapped and runs to the first prior agreeing on K2 and the masked K3 while the list has room, else to the chain end; chains are finite (one home's records, at most CAP_H), so every walk terminates on every run and no overflow event exists. Each new record is therefore linked to its most recent176-bit-agreeing predecessor on every run, and the listed pairs are a subset of the176-bit-matching pairs at their home, so the candidate, request and relay rows (d_cl, d_rq, d_rl) apply unchanged. The candidate list is cumulative with its2^24 cap checked at every attempted append; the final DONE marker of each group stores the list cursor to the receiver's list-end cell before the release STLR, and the main context loads it after acquiring the DONE cell, so the group's new candidates are exactly the entries between the previous and the published end. routing_check.py executes the receiver on the true layout (persistent resident registers, x12 checked after every call, heads at7*2^61, toy K2 below2^6 except top bits so distinct K2 share chains, planted cross-group same-key and same-K2/other-masked-K3 pairs, full home arrays and full lists): heads, next words, listed pairs, the published list end and the exact per-call count equal the path model, every head-following reaches exactly its index's members newest first, and every same-key class is connected by listed edges while the list is not full.

## 13. Grid, masks and provenance

SIDE=2^35,N=2^70,q2=2^28 and the35+29/35+1+28 metadata allocation are unchanged. The home digits are fraction bits[29,99) (d2 then d1), K2 is bits[128,192) with hash index K2 bits[0,58), and the key is bits[16,192); FRACTION_LABEL.md bounds every cell.

Lineage (provenance only, re-derived here; no acceptance is a premise): 44137210, 76ceba30, 37b8edde, 7b724c27, secp46385d36, a0d1722d, ba0f040b, 7e5b1fb9, secp ae21db39, e8c1f92c, 26e42f93 and ad76eaf8, aburan28 9ecd4df5 and ecf155d9, 008389f5, faf905b6, our gated 1c475fc0 and 6dd6d114. This version: section17 gate, q1=893476688, row split at2^29, 58-bit hash index (heads at7*2^61), 32-byte records with K2 bits58..63 in K3'.

## 14. Complete current cost and storage binding

range_ledger.py is the exact rational ledger. Its per-group2534/4-incidence rows include forward arithmetic and the fraction-label backward arithmetic with inline pair emission on the tail-free path, weighted by the non-stopping bounds of section15, plus the PHI-weighted correction_tails row (section16); the fraction_tables row pays every per-X area and per-Y entry. RANGE_DRIVER.md separately bounds construction and every row/segment operation. The retained normalization allowance is not reduced. Both hops, the receiver's hash join, failed caps, sampling, source service, four explicit resolution stages, cold verification and output stay separately charged. RESOLUTION_PROTOCOL.md assigns the final output reduction to unused fixed phase4 control space, not an uncharged all-to-root fan-in.

fixed_layout.py gives simultaneous nonoverlapping addresses for every source, label, tree/sort, batch, candidate/request/relay/home/flag/context/port/cursor region, and context pointers are initialized to those exact bases. The former d4-counter range[2^51,2^52) stays reserved and unused; per-group hop1 slots, the home record array (below7*2^61) and the hash heads at[7*2^61,2^64) follow2^52 (the slots from2^52+32). The sender header port and its local cursor table are separate allocations: the cursor-table base is a multiple of2^38 so the emitter forms a cursor address by inserting d1 into bits3..37, and the payload windows occupy an8192-byte stride per destination from an8192-byte-aligned base. The batch has a64-byte prefix reservation; output mailboxes follow the2^20 code/stack region with separate4096-byte child slots and isolated flag granules, and no flag is reinitialized after publication begins. Peak across2^70 nodes is at most2^134bytes, conservatively declared below2^136.

The expected total and the outward claim are those of section12.

The remaining named appendices are in proof-appendices.md and are part of this same claim.

## 15. Grouped early abort, stopping time and the expected-work ledger

**Groups.** The q1 X-rows of every node are partitioned once, before the run, into G=526 contiguous groups: G1=316 groups of1698958 or1698959 rows cover exactly rows[0,2^29), the other210 groups rows[2^29,q1) (sizes1698122 or1698123; range_ledger.py group_sizes, qg=1698959). All nodes run group g at the same time, separated by barriers: hop1 of the group's rows into fresh per-group slots (CAP_PG=5*2^17 records, BLOCK_PROTOCOL.md), hop2 into the cumulative home record arrays with the receiver's hash join at arrival, the four resolution stages for the group's new candidates, and the stop test. The stop flag is a latched word set only by a home whose resolver has verified [x]P=Q, broadcast with the barrier release and tested by every node before the next group (5 instructions, phase-control row); a set flag halts the run with that x as output. No group is skipped for any other reason. Rows, scalars, points, normalization, the range tree, coins and caps are those of sections2 to6; only the execution order and the detection pass change.

**Stopping time.** Let T be the first group whose resolution confirms an x (T=G+1 if none). Group g's arithmetic, emission, hop1, hop2 (with the receiver's append and insert) and per-group fixed work are deterministic given g and are executed iff T>g-1. The success event is {T<=G}, and the section6 expression is evaluated for the full run: the chain detection does not depend on group boundaries (every record walks the chain of all its priors at arrival, every accepted record is inserted exactly once, the candidate list is cumulative, and a group's candidates are resolved before its stop test). Thus Pr(T<=G)>=0.92608 at q1=893476688 (range_ledger.py); section17 gates the run.

**Non-stopping bounds.** Write m_g for the number of rows in groups1..g-1 (m_1=0). Pr(T>g-1) is at most the probability that no useful collision among the first m_g rows of every node is resolved by the end of group g-1. The section6 expression re-evaluated at q1=m_g (same N and q2; CAP_PG, CAP_H and the candidate, request and relay caps fixed at their full-run values; calc(m,caps) in range_ledger.py) lower-bounds that resolution probability: the Chen-Stein bound and every defect, capacity and usefulness row condition only on coins and records of prefix rows, prefix means are at most full means, the fixed caps keep every Chebyshev gap at least its prefix value, and later groups never unlink a prefix record from its chain nor remove a listed pair. Hence PHI[g]=min(1,1-lo(m_g)), rounded up to a multiple of2^-120, bounds Pr(T>g-1); PHI[1]=1, nonincreasing, PHI[526]<0.07451; the group E-factor sum_g PHI[g]*q_g/q1 is below0.52607.

**Expected work.** By linearity over the indicators1[T>g-1], E[work]<=W_rest+sum_g PHI[g]*W_g+W_sameX+W_K2+W_def+W_tails, where W_tails is the two PHI-weighted tail rows of section16 and W_g is the group's per-row and per-record work (arithmetic, emission, the39479 inversion, the1024 row allowance and the2^17 driver per row; hop1, hop2 and the receiver's append and insert per record) plus its main-pair walk bound (below), W_rest holds the node-initialization, hash-head-zeroing and range-tree rows, the per-group fixed work of all526 groups, the four phase-4 stages per group, generation, sampling, normalization/sorting, the fraction tables and the resolution, source-service, barrier, control and output rows, all at their worst case, and the global rows are charged without PHI. No independence between stopping and the chains is claimed.

**Walk means (history-conditional domination with a validity split).** Let H_g=(all rows of groups<g at every node, all Y generators, every node's b and R_g, Q). The event {T>g-1} is a function of H_g: records, chains, listed pairs, requests, replies, the resolver's sign test, the usefulness test and [x]P=Q are deterministic functions of those points, scalars and Q. Group g's walk cost at a home is sum over its inserts u of the visits V_u, and V_u is at most the number of priors v in u's hash chain (every visit examines one distinct chain member, section12). Visits are counted against the128-bit cell (home70 bits plus the58 hash-index bits of K2, bits[0,35),[64,99),[128,186)): every prior in u's chain lies in that cell, and one agreeing on the index but not on K3' is a plain visit. Split the pairs (u,v): main pairs, whose leaves come from different X-rows (in earlier groups or in group g itself), and same-X pairs (same row, hence the same group and node). For a main pair condition on H_g and on every group-g row except u's; v's row is then conditioned whether it lies in an earlier group or in group g, and v's actual stored label, a deterministic function of its own row's batch (X_v and all Y) and the earlier transport, is one fixed label. If u's row is not defective, u's record key lies in the at-most-two-value key set of the class of X_u+sY_u, where X_u=R_g+[a_k]P with a_k fresh, uniform and independent of everything conditioned on: given Y_u and the sign, X_u+sY_u is uniform on the group, so the probability that this set meets the one fixed128-bit cell of v's label is at most CELL130=(2+2^-23)*2^128/(r-1) (FRACTION_LABEL.md; the cell's lowest fixed bit is29); X_u=O (probability1/r) leaves a fixed leaf. So Pr(v in u's chain, u valid | conditioning)<=PAIRM=CELL130+1/r uniformly, whatever v's own validity. Which of its key values u's arithmetic emits is irrelevant: the event uses only the set. If u's row is defective (unconditional d_batch per row), u's label is garbage and its walk is at most CAP_H visits (the home array holds all of its chain), each possibly K2-verified: charged globally as DEF=(4+5)*d_batch*CAP_H*S. Hence E[1[T>g-1]*walks_g^main]<=PHI[g]*4*PAIRM*(s_g*p_g+s_g(s_g-1)/2-sameX_g)+DEF, with s_g the group's leaves over all nodes, p_g the earlier groups' leaves and the pair count deterministic (arrivals never exceed sends): the conditional bound is established before multiplying by PHI[g]. Same-X pairs are charged globally without PHI and WITHOUT fixing either actual label: the two records share their row's batch-prefix arithmetic, so neither raw representative is fixed when the other's generator is exposed. Both actual keys lie in the fixed key sets K_u and K_v of their ideal classes, so the records share a chain only if some128-bit cell meets both sets, and Pr<=sum over cells C of Pr(K_u meets C and K_v meets C). For two leaves with different Y, condition on X: the ideal points are independent uniform (Y_u, Y_v fresh), so the sum is at most max_C Pr(K_v meets C) times the expected number of cells met by K_u, at most2*CELL130. For the two leaves of one incidence, (X+Y,X-Y) is an invertible linear image of(X,Y) with Y nonzero, uniform on pairs of distinct points, so given the first leaf the second is uniform off one point and the same bound2*CELL130 holds. With X=O (probability1/r) the two leaves of an incidence share a class: PAIRM_SAME=2*CELL130+1/r (SXV in the ledger; the same two-cell union as M=2*KB176). K2 verifications drop the stopping indicator (E[1*X]<=E[X]) and are charged globally over ALL unordered pairs, main and same-row alike, at the two-cell union2*KB(134), KB(134)=(2+2^-23)*2^122/(r-1): EK2=5*PAIRS*2*KB(134); for main pairs the single cell would suffice by the conditioning above, and the union dominates both cases. A visit costs its probe and pointer test (4), a full-key match its list-cap test (2, charged as5 per K2 match); a listing replaces the pointer test and with the list full the walk continues, so no listing row is needed.

**Equivalence with the ungrouped success analysis.** Same records (same kernels and row order; groups are contiguous row ranges; hop1 cursors reset per group as formerly per phase); same home-array contents in arrival order (cumulative cursor, cap and drops); every accepted record inserted exactly once at arrival; same-key records chain-connected unless a176-bit intruder arrives between them (the intr row); listed pairs resolved by the same protocol; a useful class edge exists whenever a useful pair exists (section5); stopping only skips groups after a confirmed x. So the section6 bound holds for the grouped run, and the E-ledger charges confirmed stopping only.

**Rows.** range_ledger.py uses the protected subgroup order r (asserted below the field modulus) in P1, rho, every KB bound, the prefix success bounds and PHI. It evaluates, exactly in rationals, W_g=NODES*q_g*(2534/4*Q2+39479+1024+2^17)+s_g*(291/128+129/32+11)+4*PAIRM*(main pairs of g); globally without PHI, SXV=4*PAIRM_SAME*(same-row pairs), EK2 and DEF; with PHI the tail rows of section16; the worst-case group_fixed row NODES*526*(1236*SIDE+27), hash-head zeroing2^44*1026+1, and the unchanged worst-case rows of section9 (four phase-4 stages per group,2104 barriers (4 per group), phase control2^20, the per-group stop broadcast in the output row). The gated total E=64+(K/2^64)*E_run (section17) has integer ceiling43307381650312544483650326063708894449647, log2(E/1000)=125.0259524610190082158..., and the claim125025952461019009/10^15 is the least15-decimal value at or above the rational enclosure (score_bounds). With PHI=1 everywhere the run rows give less than2^127.2291; with every walk at its CAP_H worst case, every correction tail taken and all groups run, the disclosed per-run worst case is2^180.6801 units (an outward4-decimal ceiling of a certified logarithm upper bound, asserted in range_ledger.py). Preprocessing (node init, generators, sampling, hash-head zeroing, range tree) is below2^119 units and inside the total.

## Short Montgomery multiplication and exact cost change

The scalar multiplication core is adapted from Amazon s2n-bignum commit56152cdb51f47ef9f63b179cda0a8b7518379240, arm/p256/bignum_montmul_p256_alt.S, licensed Apache-2.0 OR ISC OR MIT-0 (see the repository LICENSE). The upstream HOL proof does not certify our patched tails or register transformations.

The core forms the full four-by-four schoolbook product with MUL/UMULH and explicit carries, then reduces the low limbs four times with the cancellation digit t (MC=2^64-2^32+1: add t shifted32 to the next limb, carry t>>32 onward, add the low/high halves of t*MC to the remaining two limbs with carry), each step an exact division of the current integer plus t*p by2^64. This gives raw T=(ab+m*p)/R0 with0<=m<R0, equal to the previous multiply's raw value for all a,b<R0 by uniqueness of m modulo R0. B.CC consumes the final high carry: the weak tail keeps T if T<R0 and returns T-p otherwise, always in[0,R0). Old and new correction predicates agree pathwise.

Including RET and standalone tail constants, MUL costs110 instructions when correcting and104 otherwise (probe short_mul); with the MC materialization hoisted, the rolled worst-case bounds are forward122 and backward616+4*carry. For a zero denominator the forward pass completes with bounded garbage within the441/4-instruction incidence allowance and its fixed46, under the already-charged d_batch exclusion. RHS limbs occupy x7,x8,x13,x14 and are read-only after loads; MC is a distinct scratch register standalone and the reserved x2 in generated loops. The forward specialization maps the live canonical difference into that RHS bank and conflicting scratch into x24..x27; the backward square reads the cached inverse x24..x27.

The inherited rolled listing is pinned to SHA25607f98f57a417d080a5ccdf49fc8d1381b501b174251f89606e5ea69a48b95a29; this version derives pair_bwd_f from it and applies the exact eightfold transformation of section3, giving2614/4 per incidence, with the driver's2^17 row covering segment fixed terms. The generated listing SHA256 is1236f74b9c333702bd57a9c5981e4503c8f39a59272a5ac657f8929275f48fdd, encoded size28048 bytes.

## 16. Correction tails

proof-tails.md (part of this claim) gives the correction_tails row and the disclosed heuristic H_T.
