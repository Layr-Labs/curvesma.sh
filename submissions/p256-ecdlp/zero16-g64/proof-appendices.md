# Appendix: NORMALIZATION.md

# Signed normalization with original-ID range traversal

Normalization changes signs and traversal, not the underlying class multiset.

## 1. Exact point and scalar representation

Generate the original independent scalars exactly as in the baseline: X_k=[a_k]P+[b]Q, Y_l=[c_l]P. Preserve b during all X construction. Before hot work, for each nonidentity point A choose epsilon(A)=+1 if canonical Montgomery y_A lies in [1,(p-1)/2], else -1. Odd prime order excludes nonidentity y=0. Replace A by epsilon(A)A. For X store a'_k=epsilon(X_k)a_k mod r and b'_k=epsilon(X_k)b mod r; for Y store c'_l=epsilon(Y_l)c_l mod r. Identity X keeps its pre-row detection/defective-marking handling and never enters the arithmetic kernel. Normalization takes canonical p-y, not a redundant hot representation.

Keep source Y coordinates and scalar labels indexed by original l. Range-tree records carry that original l explicitly; sorting permutes only arithmetic traversal, never probabilistic anchors. X keeps original k and its group. Tree and descriptor construction are in RANGE_DRIVER.md.

For every original incidence and s in {+1,-1},

X'_k+sY'_l = epsilon(X_k) [ X_k + s epsilon(X_k)epsilon(Y_l) Y_l ].

Thus each incidence's two negation classes are exactly unchanged as a multiset, and so is the full ideal leaf matrix indexed by original generators. Analyze independence and collision counts on those original generators, not on falsely uniform normalized points. Sorting permutes leaves within each full row; normalization may exchange its two signs. Weak raw labels may change with prefix order, but the deterministic U(K) bound in FILTER_FREE.md is representation-history independent and still dominates each label. All good-witness, load and pair-envelope statements therefore need only this unchanged class matrix and U(K), not execution-order independence. Capacity predicates are fixed point-only predicates; sorted-prefix order is not independent randomness.

The event x_X=x_Y and identity X is invariant. Batch-defect marking stays before inversion; defective rows complete as bounded garbage. Intruder/false-pair bounds cover all pairs, so reordering equal-key runs cannot increase their envelope. For good canonical witness labels, the adjacent useful-pair lemma applies to any ordering: if every adjacent coefficient pair has zero resolver denominator, transitivity makes every pair useless.

## 2. Hot kernels

RANGE_DRIVER.md partitions each row into at most87 segments with known x-difference direction, y-difference direction and x-sum carry. range_witness.generated_assembly defines the fixed loops. Across segments use a single prefix chain and inverse, then reverse the descriptors. Forward costs441/4 per incidence with garbage completion (no zero test; X hoisted, post-incremented stores, ping-pong prefix banks on mains, pinned tail in fixed); fused backward2505/4 plus4 for a high-carry incidence with inline pair emission and shared-subroutine flushes charged per firing (post-decremented pointers, fused pair-loads, b16 first-operand reuse, b17 TBZ/BFI emission). Segment fixed terms and dispatch are covered by the additional2^17 per-row driver charge. Normalized y sum remains below p. All arithmetic/source/ID invariants are retained. Planted algebraic fixtures exercise rare raw labels; actual curve fixtures are separate.

## 3. Resolution and uselessness

An emitted normalized plus/minus leaf has coefficients (a'_k+s c'_l,b'_k). Fetch row-specific b'_k rather than the node's global b. The reply remains224 bytes: X, Y, a', b', c'. A collision Z2=epsilon Z1 yields
x=(epsilon A1-A2)/(B2-epsilon B1) mod r,
followed by [x]P=Q verification. All operations and fallback/infinity cases remain charged.

Condition on original generator points and Q. Since a is uniform modulo r, original b is uniform and independent of X. The normalization signs are fixed by these conditioned points, so signed b'_k is uniform too. Across distinct nodes the coefficient pair stays independent. More strongly, if original collision sign is u and the corresponding normalized leaves have signs epsilon_X1 and epsilon_X2, then u'=epsilon_X2*u*epsilon_X1. Both resolver numerator and denominator become epsilon_X2 times the original ones. Vanishing and the recovered scalar are thus preserved pointwise, keeping the cross-node1/r uselessness bound. Same-node loads use the original matrix coupling, not independence of row-specific b' within one node.

## 4. Cold work and memory

Retain the previous1280q1+3968q2+4352 normalization/sort/search allowance in full, even where the new driver makes parts redundant. Additionally charge8192q2+2^17 q1 per node for the range tree and driver. No saving is taken from the old reserve. RANGE_DRIVER.md itemizes the added construction. fixed_layout.py allocates source tables, labels, all tree records/nodes/sort scratch and batch storage simultaneously below the fixed-memory boundary; nothing relies on unproved reuse. Totals are in range_ledger.py and proof.md.


# Appendix: FILTER_FREE.md

# Filter-free labels and complete load domination

Write R0=2^256, delta=R0-p and k for the canonical Montgomery x-key of an ideal nonidentity leaf. Set U(O)=empty. Hot arithmetic and its denominator/identity exception charge are unchanged.

## 1. Exact output set

The raw Montgomery square satisfies 0<=T<R0+p. The common-sum representation c lies in[0,R0), hence z=T-c lies in(-R0,R0+p), and h=floor(z/R0) is -1,0 or1. The assembly stores u=(z-h*p) modulo R0. Its canonical field answer is k=z modulo p.

- For h=-1 and z>=-p, u=z+p=k. If z<-p, then u=z+p+R0=k+delta, with k=z+2p>p-delta.
- For h=0, u=z. If z<p this is k; otherwise u=k+p and k<delta.
- For h=1, u=z-p lies in[delta,R0). It is k if below p; otherwise it is k+p with k<delta.

Thus every output belongs to U(k)={k}, with k+p added when k<delta and k+delta added when k>p-delta. Since delta<p/2, the exceptional bands are disjoint and |U(k)|<=2. Every extra label is in[p,R0), not arbitrary machine garbage. This covers the noncongruent negative-wrap case and redundant congruent residues. The emitter outputs u without a representation test, subject only to its fully charged slot cap.

## 2. Witness exception loss

Let Bbad=[0,delta) union(p-delta,p). A true negation collision shares its canonical k. This set has at most2delta keys and4delta curvepoints, each with at most two negation-collision partners. For two independent cross-node uniform points, a collision with k in Bbad has probability at most8delta/r^2<2^227/r^2. Summing over Gamma gives the retained Gamma*2^227/r^2 allowance, applied unconditionally before any capacity/survival conditioning. No global S*2^-32 term is omitted: exceptional competitors are handled by the set-valued load bounds below, not excluded globally.

## 3. Single-cell bounds

For a fixed f-bit pattern let h be the number of its fixed positions in[224,256). Because p>R0-2^224, every extra label has all top32bits set. The containing interval[R0-2^224,R0) has at most2^(224-f+h) integers satisfying the pattern. Each extra raw label has at most two inverse translated canonical keys, and each canonical key at most two curvepoints. Extra point preimages are therefore at most4*2^(224-f+h).

For the35-bit slot masks and96-bit home/bucket mask, h=0. FOURIER_CAPACITY.md bounds canonical point counts. In the96-bit case its error is below4721808*2^128, while the new extra fibers add at most4*2^128. Their sum4721812*2^128 is below2^151, the spare1/512 fraction of2^160. For35bits, canonical error plus4*2^189 is also below2^221/512. Both refined coefficients remain513/512.

For home70, the crude canonical count2*2^186 plus4*2^154 gives coefficient2+2^-30. The group144 and matched182 masks include sixteen fixed positions in[224,256), so their coefficient is2+2^-14. Divide these counts by r-1 where a free generator excludes one group point. These are mask-specific bounds, not claims about arbitrary f-bit subsets.

Use the fixed indicator A_C(P)=1[U(k(P)) intersects C]. It is deterministic in the ideal point, even when the actual representation depends on a whole batch. Every actual record in C is dominated by this indicator. Define all ideal variables under the original independent generator distribution; do not condition on clean execution or chosen raw representatives.

## 4. Pair loads and dependence

For a fixed opposite ideal point P', U(k(P')) has at most two labels. A free point P whose possible labels can match P' on182bits lies in a union of at most two fixed cells. Hence the matched-degree bound is M=2*KB182. Home and relay-cell factors use the single-cell KB bounds. With rho=(r/(r-1))^4, retain the source's conservative shared-generator cases with

h=rho*(S^2*KB70*M+(2S+20)*M),
g=2*rho*n0*S*M,
v=18*rho*((S/SIDE)*S*KB35*M+(2S+20)*M).

The witness is charged separately if bad; otherwise its actual key is its unique canonical label and its home/group are fixed by the witness generators. Thus its intruder rate is rho*(S+12)*KB144, without another factor2. Slot and bucket indicators remain fixed point predicates, so the source covariance argument applies: independent incidences/free Y cases, the one-excluded-point bound for a common X, invertible two-sign forms for a shared incidence, and the conservative60-partner allowance for witness-generator incidences. The phase slot uses qphi; the persistent bucket uses all q1.

Actual partial/full matches can be false, including between two exceptional labels. They are still members of the ideal all-unordered-pair envelopes above. Caps and verification account for them. Removing records can create adjacent pairs, so no subset claim about an earlier adjacency list is used. Source labels, signs, generator dependencies and exact final verification are unchanged.

## 5. Work and evidence

The block emitter follows BLOCK_PROTOCOL.md; its full-record envelope already charges up to two records per incidence, so formerly filtered outputs are not free new work. Exact success and cost rows, including the enlarged h/g/v bounds, are recomputed in range_ledger.py. Boundary and differential tests assert membership in U(k), not equality of every raw output with k; they supplement this argument and do not prove production success.


# Appendix: FOURIER_CAPACITY.md

# Fixed-bit capacity lemma

Let p be the protected P-256 field prime, r its prime curve order, R0=2^256, and rho=(R0 mod p)^(-1) mod p. For canonical Montgomery keyk, the number of affine curvepoints is1+chi((rho*k)^3-3rho*k+b). Define g(k) as this quadratic character for0<=k<p and0 forp<=k<R0.

Primary reference for the mixed character bound: Todd Cochrane, Exponential sums modulo prime powers, Acta Arithmetica101.2(2002), introduction equation(1.2), author PDF https://www.math.k-state.edu/~cochrane/research/expprpowa.pdf . A nonzero additive frequency has polynomial degree1 and the multiplicative cubic degree3, characterorder2, giving3sqrt(p). At zero frequency Hasse gives2sqrt(p). The cubic is nonsingular and not a square.

An arbitrary realfrequency exponential over0..p-1 expands in the p additive characters. The geometric coefficients satisfy min(1,1/(2p||theta-t/p||)). Counting at mosttwo frequencygridpoints per distanceannulus bounds their normalizedL1 sum by2+H_floor(p/2)<3+lnp; the deliberately looser4+2lnp is also safe. Hence every R0-DFT coefficient of zero-extendedg has magnitude at most3sqrt(p)(4+2lnp).

For an interval of lengthL in Z/q, its normalized FourierL1 norm is at most2+lnL: use|D_L(t/q)|<=min(L,q/(2|t|)) and the decreasing-sum integral. Fixing a bitblock[a,b) is a periodic interval of length2^a and period2^b; periodic lifting preserves normalizedL1. The lowest-bit fixedblock has norm exactly1. The complete bitmask is a POINTWISE PRODUCT of these indicators, so the Fourier convolution inequality bounds its normalizedL1 by the product of these norms. Do not use the false normalizedL1 product claim for pointwise Fourier products of AP convolutions.

For SIDE35 and d3=26, bucket bits are[0,35),[64,99),[128,154). There are2^160 total integerkeys in the cell. Its indicator has norm at most(2+64ln2)(2+128ln2). Pairing againstg gives pointcount <=2^160+3sqrt(p)(4+2lnp)(2+64ln2)(2+128ln2). Zero extension avoids wrap or proper-GAP assumptions. The count of canonical integers in the cell is at most2^160.

Safe rational constants: ln2<7/10, lnp<180, sqrt(p)<2^128,2+64ln2<47,2+128ln2<92. Error<3*364*47*92*2^128=4721808*2^128<2^151. Thus use KB96=(513/512)*2^160/(r-1). For either35-bit slot, norm<=47, error<3*364*47*2^128<2^144, far below2^221/512; the same513/512 coefficient is safe. Retain crude two-point fiber bounds for all other patterns, particularly144and182bits. This explicit lemma is used by proof.md; its argument and cited character-sum theorem are material evidence.

Conditional leaves with one free independent generator are uniform or uniform excluding one group point; denominatorr-1 is retained. Recheck every dependency/variance step and every positive capgap when applying the new KBs. Do not merelysquare a negative gap.


# Appendix: RECEIVER_PROTOCOL.md

# Grouped direct-bucket receiver machine

This specifies the changed transport component. It is an abstract scalar AArch64 instruction-work construction, not a claim about measured network throughput. All processors' instructions and the stated delivery/in-flight work are charged. No direct-home shortcut or uncharged key oracle is used.

## 1. Contexts, channels and bounded state

Each algorithmic node has its main scalar context and one dedicated scalar receiver context sharing that node's coherent private memory. The receiver is resident, not invoked by a per-packet software dispatcher. Its code and register parameters are fixed before hop2 begins and persist across both phases. It executes independently while the main context stalls on a send or waits for completion. Both contexts' work is included.

Final-cursor, barrier and phase4 messages keep the original addressed delivery handler. Hop1 follows BLOCK_PROTOCOL.md with its fixed payload windows and completion contract. Hop2 payloads use a separate fixed row-payload endpoint/port bound to the receiver window; control messages cannot enter it (a fixed channel distinction, not content-dependent classification). Extra context, stack,48-byte packet window, DONE cell in its own4096-byte granule, done-count register and payload-port state need less than8KiB per node, inside the fixed-region reservation. Hop2 senders keep no per-message in-flight state; hop1/block ports and theirs are unchanged and charged separately.

Use the source blocking-port convention: a sender's completing STP waits if its bounded port/window cannot accept another message; the receiver's first LDP waits until a complete packet is atomically published. A waiting instruction counts once, not repeatedly as a poll. Arbitration is fair, backpressured and the same transport mechanism as the source delivery construction. There is one bounded receiver window, no hidden unbounded software queue, and no software context-switch loop. The packet is immutable until all five words are loaded. The two payload LDPs never release the window; only the final flag LDR, the last of the fixed five-word read, does, after every word is already in receiver registers. A later publication therefore cannot corrupt in-progress processing. Publication and the final blocking load have acquire ordering; the window is an ordered device endpoint, not an ordinary unsynchronized buffer.

## 2. Receiver body and DONE-flag completion

Resident inputs are the packet address, d3 mask, cursor/end-table base, DONE-cell address, done count, sender count, and constants1,32,-32. The body reads the four payload words in two LDPs, then loads the address word and branches on its flag: data always carries zero there (the streaming publisher stores xzr), DONE markers carry one. A data packet computes d3, addresses the interleaved cursor/end pair, reserves32bytes with LDADD, loads the end and compares. On a full bucket it undoes the increment. Otherwise it stores the two payload pairs; no d4 counter exists in this version.

Both data exits drop the old per-message STLR and branch to the blocking first load. The finite callable test uses RET for that final branch at the same one-instruction price. The accepted path costs13 (the former d4 mask, counter address and STADD are deleted; flag load and flag branch retained, no STLR); a full-bucket drop costs12. Charge13 for every potential record, including ones that never reach hop2 through earlier drops. A DONE marker instead increments the done count and compares against the sender count: non-final costs8, final costs10 (STLR of1 to the DONE cell plus count reset). STLR is a base scalar release store with an unindexed address; its encoding is restricted in the interpreter.

Hop2 has no per-message acknowledgment or in-flight update. Delivery completion still means processing completion, not inbox arrival: the receiver is strictly sequential, so consuming a sender's DONE marker implies every earlier packet from that sender is fully processed, and reaching the sender count implies every stream is fully processed. The main context awaits the DONE cell with the same LDAXR/CMP/B.EQ/WFE primitive as barrier epoch words (4 per phase, explicitly charged); the DONE cell sits in its own4096-byte granule with no unrelated wake, exactly as for barrier words. The main context then resets the cell with one store before the next phase. This DONE/completion contract is part of the machine definition.

There is no per-message call, argument setup, register save/restore or polling: parameters stay in the dedicated context. Startup and public parameter/address setup are bounded by the separate2^14 per-node phase/context-control allowance. A final issued blocking receive after the last packet costs at most one more instruction and is included there; with no further messages it executes no further work. Program, context and port initialization remain within the charged2^20 node initialization row. DONE publication/handling, the DONE signal, the main-context wait and both resets sit in the separate hop2 DONE fixed row, not hidden in those allowances. These bounds do not excuse an implementation adding a software scheduler; that would need its own extra ledger.

## 3. Phases and progress

Partition absolute X indices into G=64 contiguous groups of at mostceil(q1/64) indices (proof.md section15). Keep absolute k in w2. Per group initialize the block windows and local/remote cursor tables as in BLOCK_PROTOCOL.md; emit full and partial blocks and final valid cursors; drain and barrier. Then h2s streams those local slots through the row-payload port, and h2done publishes one DONE marker to every row destination, including empty ones. Extract destination d2 before BFI replaces w3bits0..34 with original source column j. Sign/l and w2 are untouched.

The receiver never sends another payload or waits on main-context work while processing a packet: it only performs bounded local operations, and its single STLR targets its own node's DONE cell, which never blocks. Fair execution therefore drains finite sender workloads even if main contexts block on their ports. After h2s and h2done, each main context awaits its receiver's DONE signal and all nodes barrier. Only then may hop1 slot storage be reused by the next group. Every sender publishes exactly one DONE marker per destination per group by construction, so the count reaches the sender count exactly when all streams are processed; empty senders still publish, so no receiver hangs. Buckets, their cursors and read cursors, the hash heads, next words and the candidate list persist through all groups. The hash join hchg of a group starts only after that group's hop2 fully drains and the barrier completes.

Every group's first-hop cursor initialization and finishing sweep are charged. h2s costs7n+2floor(n/8)+11*SIDE+2 per group; h2done costs3+6*SIDE per group. DONE handling costs8 per non-final and10 per final marker per receiver per group; the main-context wait costs4, the cell reset1, and the initial preset1 per node. Bucket cursor/end initialization costs5*2^26+1 and read-cursor initialization4*2^26+1, once. After each group's drain, hchg costs13 per bucket plus11 per new record plus2 per8-record pass plus its walks, verifications and listings (proof.md section15). All per-group fixed terms are charged in the worst-case group_fixed row of range_ledger.py for all64 groups,1536*SIDE+14*2^26+64 per group per node.2048 global barriers at1024instructions/node and2^20 phase/context-control instructions per node conservatively cover the finite public group transitions, kernel arguments, stop tests and completion setup. Receiver record work is not hidden inside those fixed allowances.

## 4. Addresses and resource obligations

Reserve[0,2^52) for fixed memory; add the d3 end table, the4096-byte DONE granule and8KiB extra receiver state explicitly. Place the per-group hop1 slots, tight32*CAP_3 full-record buckets, the8-byte-per-record next-array and the2^59 hash head words in disjoint increasing regions. range_ledger.py and fixed_layout.py check the exact endpoint below2^64, fixed reservation, positive cap gaps and attempted hop1 cursor bound. The bucket increment/rollback cannot wrap because its end is strictly below the subsequent next-array region and the final address limit. No second-hop R2 receive array or final-cursor table is needed. The per-group hop1 slot region is reused only after the barrier that follows the group's hop2 drain; nothing else relies on reusing dead storage.

Finite tests run actual emit, first-hop delivery, final cursors, sender, DONE markers, receiver DONE signaling, the hash-chain join and a toy stop rule over sixteen nodes and three executed groups of four planned, with filter-free high raw labels, slot/bucket/list caps and independent reference contents/counts. They establish the stated finite behavior and instruction counts, not a physical asynchronous implementation. The protocol and all ideal-machine assumptions remain explicit review obligations.


# Appendix: COLD_ACCOUNTING.md

# Conservative complete cold allowances

This analytic construction refines unused reserves, not measured production throughput. It uses the exact P-256 field primitives already supplied and checked: multiply116, square92, add23, subtract19, inverse39479, each including its return. Scalars and field values are four64-bit words. No arbitrary-width operation is a unit-cost primitive.

## Affine add/double, including exceptional points

Cold temporary frames represent a point by two canonical Montgomery field coordinates and one infinity flag. Persistent source_X stores infinity in-band as(0,0); nonidentity points cannot have y=0 in the protected odd-order group. Encoding/decoding this sentinel is charged in the existing construction/row control allowances. Check the flags, then coordinate equality. Copy the other point for an identity input; return infinity for opposite points or doubling with y=0. Otherwise use the ordinary slope formulas. Addition uses two differences, one inverse, two multiplies, one square and the remaining coordinate differences. Doubling uses slope=(3*x^2-3)/(2*y): two squares, two multiplies, one inverse and at most eight additions/subtractions. The constants1 and3 are stored in Montgomery representation. Thus either nontrivial path uses at most16 field calls, at most one an inverse. The other fifteen may all be charged at116, deliberately exceeding the actual mix.

A direct wrapper places at most16 temporary fields in a fixed frame, keeps live pointers/counters in frame slots across calls, and uses no dynamic allocation. Charge16 non-field instructions per call (pointer loads/address formation, BL and live-state management),128 for flag/coordinate comparisons, branches and copies, and128 for frame/constants/result handling. Four-word equality/zero tests use bounded loads, XOR/OR and branches; coordinates and frame offsets are below4096. These512 wrapper instructions are a conservative explicit construction allowance. Required field loads/stores and RET are already in each field primitive count and are not omitted.

One complete affine operation therefore costs at most39479+15*116+512=41731, below45056. The larger45056 allowance also covers exceptional branches, which skip arithmetic. Storage is bounded by the fixed frame and does not alter the published node-memory bound.

## Scalar multiplication and generation

Left-to-right binary multiplication over at most256 scalar bits performs at most256 doublings and256 additions, starting from infinity. Allow65536 further instructions for bit extraction, loop control, initialization, conversion/copies and output, separately from the point wrappers. Its work is at most512*45056+65536=23134208.

Compute R_g=[b]Q once per node for this target, retaining original b and R_g. Each X then needs one [a_k]P and one affine addition with R_g; each Y needs one [c_l]P. The row is itemized: X items (q1) at 23134208+45056+262144=23441408, Y items (q2) at
23134208+262144=23396352 since they need no affine addition, R_g at 23441408, i.e.
NODES*(q1*23441408+q2*23396352+23441408). No work is removed. It does not amortize across targets, and row-specific signed b'_k never replaces original b during construction. Sampling is separate. Header/canonical conversions and fixed bookkeeping stay within262144 per item.

## Scalar-order arithmetic

For modulus r, canonical four-word modular addition can add into five words, subtract r extended by a zero word, and select the canonical result. Even charging all loads/stores and control,64 scalar instructions suffice: twelve input/modulus word loads, five additions, five subtractions, four conditional selects, four stores and at most16 pointer/control instructions total46; the64 allowance is larger. The same construction doubles a canonical value. The fifth word handles overflow; inputs below r give a result below2r and one subtraction suffices.

A256-step shift-and-add modular product doubles the accumulator each step and conditionally adds the multiplicand. Charge16 more instructions per bit for extraction/branch/loop state and64 fixed instructions. Worst work is64+256*(2*64+16)=36928, strictly below2^16. Each addition is canonical, so induction gives the exact product modulo r. Binary exponentiation to r-2 uses at most512 such products. Including32 controller instructions per product and512 fixed instructions costs at most512*(36928+32)+512=18924032, below2^25. Fermat applies because the protected scalar group has prime order r. A zero solve denominator is rejected before inversion.

## Candidate resolution

The source-service and four explicit bounded-routing rows stay separately charged in full. The resolver receives the stored generator points and scalar labels, reconstructs its two leaf points with at most two affine operations, tests the negation relation and nonzero scalar denominator, solves for x modulo r and verifies [x]P=Q. Allocate:
- scalar inverse:18924032;
- verification scalar multiplication:23134208;
- two affine reconstructions:90112;
- all remaining scalar operations, comparisons, canonical output, local message preparation and bookkeeping:2^18=262144.

The last row can cover four full modular products at36928 each, sixteen modular additions at64 each, and more than100000 instructions for the remaining fixed-width control/copy/output operations. The solve itself needs only one product after the inverse and a few modular differences/sign changes. Transport/service counts are neither reduced nor charged implicitly to this row.

The sum is42410496, below2^26. Charge42410496 for every allowed candidate, including false partial matches and failed verifications. Early exits cannot increase this bound. No claim depends on a prior submission's acceptance, or on the sharper54413312 resolution allowance in public submission1d6d2b74.

The exact integer identities are reproduced in range_ledger.py; scalar_check.py independently checks word-level products, carry-aware halving and inversion. Canonical cold field counts (multiply116, square92, add23, subtract19, inverse39479) are the COLD_ACCOUNTING.md analytical allowances, distinct from the weak hot routines; cold operations and inversion are charged at those allowances. (The never-executed canonical exhibit file was dropped to meet the archive cap; no probe behavior changes.) Expanded program delivery stays within2^20 per node. Four coordinate decodings, sign/equality checks and exceptional dispatch are included in the262144 miscellaneous resolver allowance; no second charge is silently removed. Those finite tests supplement, not replace, the word-level and loop-invariant arguments above. A full protected package replay and external exact-version review are still required before publication.


# Appendix: RANGE_DRIVER.md

# Range-tree driver decomposition

2026-11-02. Analytical construction with tested native merge/search helpers. Python finite orchestration is not itself a counted production implementation.

Build an array of72-byte records (x,y,original ID), stable-sort by(y,ID), and use its fixed balanced binary tree. At each tree node store the same records sorted by(x,ID). Build node lists bottom-up by merging child lists, not by independently sorting every node. This avoids the extra logarithm of the naive construction. Keep IDs unchanged and original scalar tables separately. Each level stores exactly n records; there are29 levels for n=2^28. Total tree records72*n*29; a64-byte descriptor per node and two72*n sorting buffers bring the explicit allocation below2^40 bytes. fixed_layout.py inserts these regions without overlap into the operative address map.

A scalar-instruction budget is128 per merged record plus64 per merge. A fixed five-word lexicographic comparison, selected72-byte record copy and pointer/count branches fit this allowance, The supplied merge.s has at most64instructions per record plus8fixed, leaving added space for caller setup in the128/64 allowance. A y sort and an x merge tree each cost at most128*28*n+64*n. Add64 per tree node and64 per leaf copy: total7488*n-64, within8192*n. No free sorting oracle, per-target amortization or construction of advice is assumed.

At query rank r, a single root-to-leaf split path leaves at most29 disjoint fully contained tree nodes covering all leaves, with a fixed y-side bit on each. In each node list, binary searches at X.x and2^256-X.x split it into at most3 segments. If X.x=0, the latter threshold is a sentinel beyond all field coordinates, not a wrapped zero. Thus at most87 descriptors suffice. Each record occurs exactly once because node intervals are disjoint and each local split is a partition. Empty segments are omitted.

There are at most29*(1+2*29)=1711 coordinate comparisons: one global rank search and two searches per cover node. With a 32-instruction comparison/search-step budget (search.s uses at most16 per step plus10 fixed),256 per segment for descriptor writes/argument setup/call overhead,64 per visited path/cover node, and4096 fixed initialization, the per-row driver bound is84768, below2^17. Each64-byte descriptor stores node-list base, start/count, batch offset, x/y direction, carry and original-ID cursor. At most48 memory/address instructions construct it; caller spills/restores and arguments cost at most128, and all three routine fixed terms total37. These fit256 per descriptor. Merge callers use at most64 instructions for public loop/address calculations, six argument moves, BL and bookkeeping. Tree nodes use at most64 for eight field stores, address arithmetic and traversal. These reserves are deliberately separate from per-item arithmetic; call overhead already charged elsewhere must not be subtracted twice. The stated constants give the bound directly; BLOCK_PROTOCOL.md conservatively allows the additional261 per row.

Forward calls consume descriptors in order while carrying the prefix across segment boundaries. Reverse calls traverse the same descriptors backwards with the one row inverse. A zero denominator marks the entire row defective before any key emission; the row completes as bounded garbage within the same envelopes. Carry count is the sum of lengths of high-carry segments; the global cap must be checked before emission, and cannot reset per group. Slot admission and source-ID recovery are separate obligations.

The tree is generated by explicit bounded loops over levels and contiguous runs, with two reusable y-sort buffers and a separately allocated list per tree node; no library sort or uncharged allocator is used. Child/node offsets come from level/range indices; at most eight64-bit descriptor words are stored; stack depth is at most29. The x-key search uses the node base, the y-rank search the y-sorted array base+32 with the same72-byte stride; the second x threshold search runs only for X.x!=0. All addresses fit the explicit layout. The instruction budgets above include these operations, including empty searches/segments, rather than relying on finite timing. The full attack proof separately treats carry caps, phase4 and source recovery.

Control suballocation, relative to the4096-byte-aligned65536-byte contexts_ports_barriers base: query descriptors[0,5568); traversal/call stack[6144,14336); resident contexts and endpoints[14336,22528); source scratch[22528,23296); eight isolated2048-byte barrier-word granules[24576,40960); DONE-cell granule[40960,45056). The gaps are reserved padding. Barrier words are2048-byte aligned and no unrelated stores share their granules. Epoch reuse follows the drain/barrier protocol. All regions are disjoint and end below65536. The caller retains original IDs and treats sentinel(0,0) X as a defective row before search.


# Appendix: RESOLUTION_PROTOCOL.md

# Explicit bounded resolution routing

This protocol fully replaces the earlier count/offset exchange. It does not assume payload conservation establishes grouping or allocation work. There is no count exchange, sender-side sorting, offset reply, per-destination linked list, or undisclosed allocator.

## 1. Record formats and endpoints

A phase-4 packet occupies an 80-byte endpoint window: destination index and unused zero address word (16 bytes), then an eight-word, 64-byte record. Publication occurs only when all five STPs have completed. The first receiving LDP blocks until the complete record is published; the window stays unchanged until STLR releases it. A busy sender store blocks. In-flight state is incremented once on acceptance and decremented only after the receiver's release. Both updates are charged.

All node coordinates are 35-bit integers. A request record contains:

- word 0: source row | (X index k << 35);
- word 1: source column | (leaf sign c << 35) | (Y index l << 36);
- word 2: home row | (candidate index j << 35);
- word 3: home column | (candidate side << 35);
- words 4 through 7: zero.

The fields fit: k < q1 < 2^29, l < 2^28, c is 0 or 1, j < CAND_CAP = 2^24, and side is 0 or 1. These come from already bounded leaf and scanner records, not untrusted network input. Each request is uniquely identified by (home row, home column, j, side).

A reply record contains:

- word 0: home row | (j << 35);
- word 1: home column | (side << 35) | (chunk << 36);
- word 2: valid flag 1, separate from every routing field;
- word 3: zero;
- words 4 through 7: one 256-bit data chunk.

There are exactly eight chunks numbered 0 through 7. Chunks 0 through 6 carry X.x, X.y, Y.x, Y.y, a'_k, b'_k and c'_l; chunk 7 is zero padding. The leaf type stays in the home's candidate record, so it need not be retransmitted.

Each phase has a distinct input endpoint and resident receiving context. All four contexts are initialized once and never reconfigured while blocked. Phases 1 and 3 have distinct cursor words although they reuse the same relay storage after a drain. The phase-2 source cursor and phase-4 indexed-home handler are separate. The earlier leaf receiver has its own endpoint. Ordinary 48-byte fixed-address delivery stays a different message class; it does not silently acquire a new ACK requirement.

Every context performs the blocking payload load BEFORE reserving a record. Thus an idle receiver cannot leave a phantom cursor advance after its last packet. After completion, its next blocked load performs no reservation or further work; that one attempted load and all context setup are in the fixed control allowance. The output/failure decision is made after all required data drains, without waiting for another packet at an idle context.

## 2. Four phases and buffer lifetimes

At a home, `make_requests` reads each candidate's two leaf-record pointers and their source metadata, then writes two 64-byte requests. It costs exactly 17 instructions per candidate plus 3. The output fits the existing 64*REQ_CAP region because 2*CAND_CAP <= REQ_CAP. This region is not the relay region.

1. Request hop 1: home (h_r,h_c) sends each request to relay (g_r,h_c), using word 0's low 35 bits. The relay appends to its 64*RELAY_CAP region.
2. Request hop 2: after all hop-1 senders drain and a barrier completes, each relay scans its accepted records once and sends to source (g_r,g_c), using word 1's low 35 bits. The source appends to the 64*REQ_CAP region. The old home request output is dead by this phase and may be overwritten. Relay and source regions are disjoint.
3. Reply hop 1: after phase 2 drains, each source retains its accepted requests while processing them sequentially. It gathers the seven requested 256-bit values, builds eight tagged replies in a 512-byte temporary and sends them to relay (h_r,g_c). The phase-1 relay region is now dead and reused, with its separate phase-3 cursor initialized at the base.
4. Reply hop 2: after phase 3 drains, each relay scans each accepted reply exactly once and sends it to home (h_r,h_c). Home placement is direct indexed storage, not append allocation.

Each phase ends with all sender ports drained and a global barrier before any receiving buffer is read or reused. No source discards its request before all eight replies have been copied into outgoing packet windows. The 512-byte output temporary can then be reused because the endpoint holds its own packet bytes, not a pointer into that temporary. At most eight replies per request are ever produced. Forwarding does not duplicate records.

The append receiver accepts the first CAP records in arrival order. It reserves 64 bytes by LDADD, compares the returned address with the end, and either stores all eight words or rolls the cursor back. Both paths release the packet and continue. Excess records are dropped individually; no whole node is aborted. The cursor stays at or below the end between packets, and its temporary end+64 address is within the reserved fixed-space padding.

## 3. Reply association and completeness

At home the receiving kernel decodes j, side and chunk from the first two words and computes

slot = ((2*j + side)*8 + chunk),
address = reply_base + 64*slot.

This is injective on the valid tag tuple. Since each request is produced once and each accepted request produces each chunk once, no two legitimate writes collide. Two candidate leaves from the same source still differ by side. Different candidates requesting the same source leaf still differ by j or home.

Before any request is emitted, all 16*CAND_CAP flag words are zeroed. Only reply records write their flag to 1. After the final drain/barrier, the `complete` kernel checks all 16 flags for a candidate before any payload is used. A missing request or any missing reply chunk therefore drops that candidate. A final chunk cannot falsely certify the earlier chunks. All flags are from this run; there is no epoch reuse without initialization. Tags have bounded widths: with j < CAND_CAP = 2^24 the largest slot is 2^28-1 and the home reply region is exactly 16*CAND_CAP*64 bytes (the kernel decodes the full 28-bit j field; the finite test exercises j = 2^28-1).

The home reconstructs the two leaves from X,Y and the retained leaf types, rejects either reconstructed identity leaf, tests equality up to negation, selects a matching sign with nonzero denominator, performs at most one scalar inversion and one target verification for that candidate, and otherwise rejects. For valid nonidentity witness leaves the matching sign is unique; exceptional or invalid candidate data may be rejected without further inversions. It verifies [x]P=Q before output. The transport returns precisely the seven values stored at the requested source indices. Source X stride is64 bytes, with x and y at offsets 0 and 32. Source Y stride is64 bytes, with x at0 and y at32. Scalar arrays have 32-byte entries; b'_k comes from its row-specific32-byte scalar slot; original b is never substituted.

## 4. Executed instruction counts

The supplied phase4.s and resolution_check.py execute these scalar routines:

- `send64`: 14 per record plus 2. It loads the chosen routing word, masks the destination, writes the routing header, copies all eight payload words, advances the input pointer and loops.
- `append64`: 13 on acceptance, 7 on a full-buffer drop, including acquire-before-reserve, cursor/cap handling, payload copies or rollback, release and the loop return/branch.
- `home64`: 16 per record, including tag decoding, indexed address calculation, all payload/flag stores and release.
- `clear_flags`: 3 per flag plus 2, using a postindexed store, SUBS and B.NE.
- `complete`: at most 67 per candidate, including all 16 loads and checks and its result.
- `make_requests`: 17 per candidate plus 3, including both leaf metadata loads and every request store.
- `gather224`: 31 per served request, including all seven 256-bit reads, corresponding writes and padding initialization.
- `build_replies`: 84 per request for all eight complete tagged records.

The finite code returns to a sentinel for counting. A production resident receiver branches back to its blocking first load instead of RET; this is the same one instruction. The model's blocking/acquire/release behavior is explicit, not something the finite emulator proves about real hardware.

For source service, decoding and checking k,l,type and locating X,Y,a,b,c need at most 64 instructions (shifts/masks, comparisons, fixed-stride MUL/ADD and constant loads). Gathering takes 31 and constructing replies 84. Allow another 256 for argument setup, spilled loop state, request iteration and scratch management. Total 435 is below the retained 1024/request allowance. Output transmission is separately charged below. Source processing uses only 256 bytes of gather scratch and 512 bytes of reply scratch, not an unallocated 8*REQ_CAP reply array.

## 5. Complete phase-4 work bound

Per phase, every input record is visited once by `send64`. Intermediate phases read contiguous append arrays; no sort or directory lookup is needed. Request creation and final completeness checking are explicit routines. Per record, the executed sender (14), receiving path (at most 16) and in-flight updates (2) total at most 32. Charge 256 per record, which also dominates per-request send-call fixed costs, record iteration and any fixed-width address/tag moves not already charged. Empty calls and context initialization are in fixed work.

The fixed allocation per node per phase is 2^37 instructions, covering conservatively request creation (17*CAND_CAP+3), flag initialization (48*CAND_CAP+2), all candidate completeness checks (67*CAND_CAP), cursor resets, scalar argument/context initialization, buffer-length calculations, state spills and phase exits. Even 512*CAND_CAP = 2^33 is below 2^37; the explicitly listed per-candidate terms sum to 132*CAND_CAP, leaving substantial fixed-control space. The separate source-service allowance is not used twice to justify missing transport work. A final constant-degree tree reduction selects one verified scalar or reports failure. The parent of (i,j) with i>0 is ((i-1)//2,j); for i=0,j>0 it is (0,(j-1)//2). Each node has at most four children with public child-slot addresses in a separate16384-byte output-mailbox region. Each of four slots has a2048-byte data granule and a distinct2048-byte flag granule, so data writes cannot wake a flag wait. All four flags are initialized to zero at node setup BEFORE the first global barrier or any possible publication, never cleared lazily when waiting begins. They are disjoint from all earlier barrier words. Each child sends one 64-byte completion record with its selected scalar and a final flag (1=found,2=empty); send_completion writes the parent axis index and public child-slot address in its routing header, then all eight payload words, in10instructions. For a row child i>0, slot=(i-1)&1; for a column child i=0,j>0, slot=2+((j-1)&1). Choose OUTPUT_BASE=2^20, exactly after the2^20-byte code/stack reservation; it is4096-byte aligned. Reserve those16384bytes before all later fixed regions, so no added alignment padding is needed. The public destination address is OUTPUT_BASE+4096*slot. Parent and slot use shifts/parity from the public row/column, at most64setup instructions. publish64 stores the seven data words and release-publishes the last packet word at child_slot+2048 in12instructions. The parent waits on that isolated flag address. Each flag has exactly one post-initialization publisher. Parent waits acquire the flag before reading the scalar; every child writes once, so the existing bounded single-wakeup wait construction applies. Each node chooses at most one valid result and forwards once. 4096instructions/node covers four child reads/selections/waits, one send/delivery and fixed addressing, inside the unused fixed allowance (132*CAND_CAP+2^20<2^37), not the global OUTPUT constant, which root serialization alone uses. Many simultaneous successful candidates therefore cause no uncharged all-to-root fan-in.

Let M be the total attempted records in a phase across all nodes, including those later dropped. Requests have M <= N*2*CAND_CAP; replies have M <= N*8*REQ_CAP. Drops can only reduce subsequent phase inputs. Therefore each phase costs at most

N*2^37 + 256*M <= N*(2^37 + 256*2^32) = N*9*2^37 < N*2^41.

The existing ledger's four EXCH rows, each N*2^41, thus bound this complete replacement. This does not assert that conservation alone establishes work: the routines and fixed-cost calculation above establish the work per record and per node; conservation then bounds how many times they execute. No per-candidate arithmetic slack is needed to cover this protocol.

## 6. Memory and success-loss compatibility

Existing fixed-space regions already include 64*RELAY_CAP relay bytes, 64*REQ_CAP source/request-output bytes, and 16*CAND_CAP*64 home reply bytes. The candidate-pointer list is separate. The relay is reused only after phase 2 drains. Source/request-output reuse occurs only after phase 1 drains. Home flags occupy their own metadata words inside reply records.

Reserve 8192 bytes, replacing the earlier 2048, for the leaf and four phase-4 resident contexts, their separate 80-byte windows and cursor/ACK words, phase-specific sender windows, 768 bytes of source scratch, and bounded stack/control state. Even six 31-register contexts (including fixed completion delivery), twelve 80-byte endpoint windows,768bytes of scratch and1024bytes of stack/control fit. Together with the separate16384-byte output-mailbox region, this increases the explicit F sum by22528; fixed_layout.py explicitly allocates these and a larger64KiB context/control region below2^52. The reserved bulk endpoint E and total memory bound do not change.

Candidate cap2^24 (reduced from2^28; its Markov loss h/(2^24-4) is charged in the success ledger and q1 re-optimized) and source cap2^29 are retained. Allocate2^34 relay records per phase. For each witness endpoint put B=rho*((S/SIDE)*S*KB35*M+(2S+20)*M), where M=2*KB182 and SIDE=2^35. A request traverses relay(source row,home column); its two endpoint requests have envelope2B. A reply traverses relay(home row,source column); its sixteen candidate chunks have envelope16B. Each direction uses the corresponding35-bit marginal bound, identical by the symmetric bit-block bounds in FOURIER_CAPACITY.md/FILTER_FREE.md. The shared-generator cases in proof.md apply separately to these fixed marginals. Both dominate ALL unordered matching pairs, including pairs newly adjacent after deletions; no adjacency-subset assertion is used. Their union has conditional nonfixed-load expectation at most18B=v.

Four anchored original generators fix at most four sign-leaves (two incidences), hence at most six fully fixed unordered pairs. Reserve108=6*18 possible request/reply records, not merely the historical16. Markov therefore bounds either required-relay failure by v/(2^34-108), which is at most the retained numerical v/(2^30-16). The stricter diagnostic denominator stays safe; it is not the actual buffer size. Drops can only remove traffic. Missing home flags arise only from one of these capped request/reply paths. No added completion probability or whole-node abort is introduced. Source-request and candidate-load bounds stay as in proof.md with their existing explicit reservations. Full work is bounded by native per-record routines plus fixed initialization and only THEN by packet-count conservation.

## 7. Target-specific resolver and evidence

COLD_ACCOUNTING.md supplies the complete42410496 P256 bound, including18924032 scalar inverse,23134208 verification multiplication, two45056 affine reconstructions and262144 miscellaneous control/scalar/coordinate-decode work. Payload coordinates are canonical Montgomery values; decode them with four field products or execute the same canonical Montgomery affine routines. Signed coefficients use only negation, additions and one final product after inversion. Zero denominators and reconstructed identities are rejected before inversion; false candidates still get the full charge.

resolution_check.py executes this protocol with35-bit nodes, original IDs and distinct row-specific b values. It tests21small capped cases, six normalized actual P256 recoveries, partial replies, maximum home tags and nine output-tree configurations. Scalar word checks cover fifth-word carry in halving. The asynchronous blocking/acquire/release and completion machine stays an explicit analytical convention, not a hardware-timing claim. Every context's instructions and all storage are charged.


# Appendix: CARRY_MEAN.md

# Carry mean and collision conditioning: proposed completion

2026-11-02, operative analytical carry-tail argument requiring independent exact-version review.

Source theorem: Todd Cochrane, Exponential Sums Modulo Prime Powers, equation(1.2), author PDF https://www.math.k-state.edu/~cochrane/research/expprpowa.pdf . For a nonzero additive linear phase and a quadratic character of a squarefree cubic, it gives the bound3sqrt(p); the degree/order hypotheses are1,3,2. This is the only external theorem used below.

For a uniformly random nonidentity P256 point, its canonical Montgomery x-coordinate K has mass (1+chi(g(k)))/(r-1), where g is the nonsingular curve cubic after the invertible Montgomery rescaling. Its nonzero additive Fourier coefficients are bounded by3sqrt(p)/(r-1). For an interval, the normalized absolute Fourier coefficient sum away from zero is at most H_floor(p/2)<=1+ln(p): pair frequencies t and p-t, and use the geometric-sum bound p/(2t). Fourier inversion therefore bounds every interval discrepancy from uniform Fp by epsilon=3sqrt(p)(1+ln(p))/(r-1). With sqrt(p)<2^128 and ln(p)<180, epsilon<543*2^128/(r-1)<2^-117.

For independent K,L of this distribution, replace L by a uniform integer U in[0,p) inside the event K+L>=R0, R0=2^256. The conditional event is an interval, costing at most epsilon. The resulting function of K is monotone in[0,1], so summation by parts and the CDF discrepancy replace K by an independent uniform V at cost at most epsilon. The uniform probability is (2p-1-R0)(2p-R0)/(2p²)<1/2. Hence carry probability<=1/2+2^-116. Sign normalization does not change x. Identity draws can be assigned zero carry and handled by their existing defective-batch loss; conditioning all generators nonidentity instead would alter independence and is not assumed here.

For collision conditioning, first fix the target and each node's shared b scalar. Uniform independent a_i makes X_i=a_iP+bQ independent uniform group points even after this fixing; c_j gives independent Y_j. A cross-node witness depends on four generator values. Fix those generator values as well; every unmentioned variable stays independent with its original distribution. This is the setting of the anchored-tail lemma in the separate carry-conditioning audit. Do not condition on the entire execution avoiding every other defect. Instead charge such defect events separately. At most four anchored row and four anchored column indices is a conservative bound per node; all incident edges contribute at most B=4(q+u).

Thus the candidate's unconditional variance denominator slack² must conservatively become(slack-B)² for witness-conditioned tail accounting, unless a tighter conditional analysis is supplied. The exact ledger now incorporates this anchored correction; its conditional total is included in the complete declared cost bound, not a measured production score. This proof must be checked together with the exact leaf/collision event definition, zero-target path, scalar rejection caps and existing Poisson loss ledger. No claim that this note alone establishes the full attack.

Carry is indexed by original unsorted generator identities. The total is invariant under Y permutation and sign normalization; never condition on a random sorted rank. Capped scalar sampling is coupled to ideal independent draws, and exhaustion of eight attempts is charged separately by the global failure term.


# Appendix: CARRY_TAIL.md

# Witness-conditioned carry tail

Fix the target, shared node scalars and the witness's original generator values as in CARRY_MEAN.md. Remaining generators are independent. Delete every incidence touching an anchored row/column; at most B=4(q1+q2) such incidences exist per node and their carry contribution is at most B deterministically. The remaining carry indicators with disjoint endpoints are independent; any covariance is at most1/4 by Cauchy-Schwarz. Counting ordered pairs that share an endpoint (including diagonal pairs) bounds variance by q1*q2*(q1+q2-1)/4. Their total mean is at most q1*q2*(1/2+2^-116).

If cap=ceil(mean)+2^50 and cap-mean>B, Chebyshev gives conditional failure at most var/(cap-mean-B)^2. This bound is uniform in anchored values, so averaging over their collision event preserves it. Only original indices are anchored; sorting and sign normalization leave total carry invariant. Do not condition on global cleanliness. Identity X contributes zero in the ideal carry array and its independent defective-batch loss is separately charged. A meter failure stops further leaf generation only; all receiver/service/control phases continue as specified in proof.md. Two witness endpoints yield the mu*(1+1/r)*2*d_m loss term. The exact ledger checks the positive gap and charges all admitted corrections fully.


# Appendix: GENERATED_PROGRAM.md

# Fixed generated arithmetic program

The four specialized backward variants and forward entries are defined by range_witness.generated_assembly from the retained reference templates and fixed literal transformations. The rolled base is pinned to SHA256 07f98f57a417d080a5ccdf49fc8d1381b501b174251f89606e5ea69a48b95a29 (zero test deleted; X-hoist and post-indexed pointers); the generated unrolled listing is pinned to SHA256 94854ef8c93104b2f1f4da04bb099be5d99a325312da3a3a356b331eb30ef9b6 (103032 encoded bytes) and includes the short-Montgomery/cached-MC transformation, X-hoisting, pointer folding, emission pair-loads, ping-pong mains, the b16 first-operand reuse and the b17 TBZ/BFI emission. Every arithmetic experiment loads and executes that generated listing; the public fixed program depends on no challenge or searched advice. Expanded code is delivered/initialized under the per-node2^20 program allowance and code/stack reservation, not assumed to occupy compressed bytes at execution. A hash pins the exact source; it does not prove correctness.


# Appendix: BLOCK_PROTOCOL.md

## Fixed windows, complete delivery and progress

Each sender has one header port and a fixed4096-byte payload window per destination d, at payload_base+8192*d+32: an8192-byte stride from an8192-byte-aligned base, whose upper4096 bytes are unused padding so that bit12 of a slot address is set exactly on the128th slot. A header is (destination, destination slot address, fixed payload-window address,128). The emitted address always equals that destination's allocated window; it is not an arbitrary remote-pointer service. An active delivery context maps only this immutable header and that fixed payload window into source space. All payload stores use destination-local hop1 memory. Fixed source windows lie below2^51; destination slots start at2^52, so their numerical ranges do not alias in this mapping.

The completing header store waits until its own transaction has been completely copied and acknowledged, not merely enqueued. The destination's dedicated scalar context is independently and fairly scheduled; its parameters stay resident. The preinitialized receiver context stores its fixed ACK-endpoint address at header+32; this is context metadata, not a fifth transmitted header word. It reads a complete header, copies128records using explicit LDP/STP instructions, and issues a native STLR1 to its fixed device ACK endpoint. That ACK is an event consuming the single current transaction, never a predicate satisfied by a stale stored1. No new transaction may replace the active mapping before that ACK. The sender can reuse the payload after ACK because the handler performs no further payload read. The trailing callable RET costs the same as the resident branch back to its blocking header load.

There is one bounded current window per receiver and one outstanding header per sender. Arbitration is fair and backpressured; no unbounded software queue, polling, epoch arithmetic or software dispatcher is hidden. The receiver never waits on main-context work or another payload before acknowledging. Its fixed finite work therefore drains finite workloads even when producers block. Publication/acquire, payload writes/release ACK and sender completion are ordered device operations. Two charged in-flight updates bracket the transaction, with decrement after ACK. Existing phase drains/barriers precede any window/cursor reuse or hop1 consumer scan. These are explicit abstract-machine assumptions, not measured physical-network performance.

## Emission, partial blocks and prefix preservation

block.s initializes local cursors to fixed window starts and remote cursors to each sender's slot base. Each record is written with two STPs; only the last of128 aligned32-byte slots has bit12 set (equivalently bits5..11 zero), tested by one TBZ. On that slot, reset the local cursor, reserve4096 remote bytes, compare against the slot end and publish only if within capacity. A full but dropped block still pays the sender flush. Local buffers persist across all row descriptors in the phase.

At group end block_finish computes the remaining valid byte count. It publishes one padded full block if the remainder is nonzero and the remote cursor is below the slot end. Padding is initialized before reads. Then it sends the clamped valid cursor min(remote+remainder,slot_end) by the ordinary final-cursor port. The completing block store already waits for ACK before this cursor is sent. CAP_PG=5*2^17 is divisible by128. Therefore every admitted full block lies wholly within its slot, a final padded block stays in that slot, and the valid cursor excludes padding. In sender order these are exactly the first min(number_emitted,CAP1) records. After the drain, h2s traverses the same source-slot order as before. Record words and original IDs are unchanged. All prior capacity, intruder, pair-load and success bounds consequently apply without assuming a new distribution.

## Complete instruction and initialization bound

The native per-incidence emitter base is25 (one TBZ per record) rather than29. Its fixed13 rather than10 costs at most3*87=261 extra per row; the previous explicit84768 driver bound becomes85029<2^17. A full flush costs at most8 sender instructions. Native block_dlv costs551 including ACK; two in-flight updates give561 total. Charge600 per full flush, reserving39 additional instructions. Across all senders and all groups full flushes are at most S/128 because each consumes128 emitted leaves, so the hop1 row is600*S/128 plus the fused subroutine's per-firing excess stated in proof.md section3; the per-incidence arithmetic/emission coefficient is derived there, not here.

Per group block_init costs775*SIDE+2, at most SIDE partial deliveries cost553*SIDE including ACK and in-flight updates, block_finish costs at most40*SIDE+4, and ordinary final-cursor delivery costs8*SIDE. Allow64*SIDE for endpoint metadata construction and32*SIDE for group/completion setup. Their sum is1472*SIDE+6 per group; together with the per-group hop1 cursor initialization3*SIDE+1, finishing sweep17*SIDE+3, h2s slot traversal11*SIDE+2, h2done3+6*SIDE, DONE handling8*SIDE+2, DONE wait/reset5, hchg per-bucket work13*2^26+1 and stop test5 this is at most1517*SIDE+13*2^26+28, charged as1536*SIDE+14*2^26+64 per group per node for all64 groups (group_fixed row). The metadata loop writes eight fixed fields per entry (payload bounds, source/channel and ACK bindings plus reserved words), using preloaded constants, incremented addresses, stores and a counted loop;64 per entry exceeds four paired zero stores, eight field stores and bounded address/control work. Constants and resident-context startup are also covered by the retained2^14 per-node setup allowance. No repeated polling is included or permitted by this construction.

All other arithmetic, cold generation, sampling, carry cap, grouping, full phase4 resolution, verification and output rows are unchanged and remain charged. The expected total is the integer ceiling84460781846257086092550172748372930554314 instructions (proof.md section15). Rational logarithm enclosure gives outward125989621112938436/10^15 at C1000. The computed success lower bound is at least0.3900000022331415 at q1=377464499.

## Simultaneous memory and finite evidence

fixed_layout.py adds8192*SIDE+32 payload bytes on an8192-byte-aligned base, the8*SIDE local cursor table on a2^38-byte-aligned base,8*SIDE remote-cursor bytes,64*SIDE metadata bytes,8*2^26 read-cursor bytes and4096 context bytes below the2^51 former counter region. The pre-counter endpoint stays below2^51. Per-group hop1 slots (32*CAP_PG*SIDE), buckets, next-array and hash heads follow2^52; the final per-node endpoint is16340889727475384320<2^64 (CAP_3 follows q1); total memory is below2^134 bytes and the declared2^136 bound. The absolute attempted remote cursor is at most2^52+32*CAP_PG*(SIDE-1)+32*2*ceil(q1/64)*2^28<2^64, checked in fixed_layout.py. Cursors reset only after the stated drain. Every endpoint and context uses the explicit bases, not an overlapping implicit allocation.

routing_check.py executes the native full/partial block path, ACK events, initialization, second-hop receiver, hash-chain grouping and the toy stop over16nodes and3 executed groups with49344 synthetic leaves, caps and independent retained-prefix references. Its event count requires exactly one native ACK per publication without a Python flag reset. These finite synchronous tests support instruction/data-path claims; they do not prove physical asynchrony or authorize a weaker completion convention. The mathematical progress argument and fixed-device contract above are material review obligations. block.s is an independent P256 implementation inspired by the public secp block/header design e59be87d-4235-4f77-ae3f-54a6ea896975; that submission's acceptance is not a premise.


# Lossless reference-assembly factoring

Only repeated literal instruction blocks in p256_inv.s and pair_inline_reference.s are factored into ordinary assembly macros; expansion preserves every nonblank noncomment line, indentation, labels, directives, operand spelling and order, adding no loop, instruction, branch or register change. The canonical expanded SHA256 values ared3b9a146237ea972367016b6b66079aa4a45c2708c28836e923f811ac53b093d and4f44679819e9d1e50306d4f31ce86604c999f0bb599871599d8f8cbdb8f80787; probe.py checks both. All expanded delivery, initialization and execution costs stay charged.
