# Appendix: NORMALIZATION.md

# Signed normalization with original-ID range traversal

Normalization changes signs and traversal, not the underlying class multiset.

## 1. Exact point and scalar representation

Generate the original independent scalars exactly as in the baseline: X_k=[a_k]P+[b]Q, Y_l=[c_l]P. Preserve b during all X construction. Before hot work, for each nonidentity point A choose epsilon(A)=+1 if canonical Montgomery y_A lies in [1,(p-1)/2], otherwise -1. Odd prime order excludes nonidentity y=0. Replace A by epsilon(A)A. For X store a'_k=epsilon(X_k)a_k mod r and b'_k=epsilon(X_k)b mod r; for Y store c'_l=epsilon(Y_l)c_l mod r. Identity X retains the original abort/identity handling and is never sent to the arithmetic kernel. Normalization takes canonical p-y, not a redundant hot representation.

Keep source Y coordinates and scalar labels indexed by original l. Range-tree records carry that original l explicitly; sorting is a permutation used only for arithmetic traversal, never a redefinition of probabilistic anchors. X keeps original k and its phase. Tree and descriptor construction are specified in RANGE_DRIVER.md.

For every original incidence and s in {+1,-1},

X'_k+sY'_l = epsilon(X_k) [ X_k + s epsilon(X_k)epsilon(Y_l) Y_l ].

Thus each incidence's two negation classes are exactly unchanged as a multiset. The full ideal leaf matrix, indexed for analysis by original generators, has the same class multiset. Analyze independence and collision counts on those original generators, not on falsely uniform normalized points. Sorting permutes leaves within each full row; normalization may exchange its two signs. The actual weak raw labels may change with prefix order, but the deterministic U(K) bound in FILTER_FREE.md is independent of representation history and still dominates each label. All good-witness, load and pair-envelope statements therefore need only this unchanged class matrix and U(K), not execution-order independence. Capacity predicates are fixed point-only predicates; sorted-prefix order is not treated as independent randomness.

The event x_X=x_Y and identity X is invariant. Zero-denominator tests remain before inversion. Intruder/false-pair bounds cover all pairs, so reordering equal-key runs cannot increase their envelope. For good canonical witness labels, the adjacent useful-pair lemma applies to any ordering: if every adjacent coefficient pair has zero resolver denominator, transitivity makes every pair useless.

## 2. Hot kernels

RANGE_DRIVER.md partitions each row into at most87 segments with known x-difference direction, y-difference direction and x-sum carry. range_witness.generated_assembly defines the fixed loops. Across segments use a single prefix chain and inverse, then reverse the descriptors. Forward costs125 per incidence; backward628 plus4 for a high-carry incidence; explicit-ID block emitter27 plus separately charged full-flush work. Segment fixed terms and dispatch are covered by the additional85029 per-row driver charge. Normalized y sum remains below p. All arithmetic/source/ID invariants are retained. Planted algebraic fixtures exercise rare raw labels; actual curve fixtures are separate.

## 3. Resolution and uselessness

An emitted normalized plus/minus leaf has coefficients (a'_k+s c'_l,b'_k). Fetch row-specific b'_k rather than the node's global b. The reply remains224 bytes: X, Y, a', b', c'. A collision Z2=epsilon Z1 yields
x=(epsilon A1-A2)/(B2-epsilon B1) mod r,
followed by [x]P=Q verification. All operations and fallback/infinity cases remain charged.

Condition on original generator points and Q. Since a is uniform modulo r, original b is uniform and independent of X. The normalization signs are fixed by these conditioned points, so signed b'_k is uniform too. Across distinct nodes the coefficient pair remains independent. More strongly, if original collision sign is u and the corresponding normalized leaves have signs epsilon_X1 and epsilon_X2, then u'=epsilon_X2*u*epsilon_X1. Both resolver numerator and denominator become epsilon_X2 times the original ones. Vanishing and the recovered scalar are therefore preserved pointwise. This preserves the cross-node1/r uselessness bound. Same-node loads use the original matrix coupling, not independence of row-specific b' within one node.

## 4. Cold work and memory

Retain the previous1280q1+3968q2+4352 normalization/sort/search allowance in full, even where the new driver makes parts redundant. Additionally charge8192q2+85029 q1 per node for the range tree and driver. No saving is taken from the old reserve. RANGE_DRIVER.md itemizes the added construction. fixed_layout.py allocates source tables, labels, all tree records/nodes/sort scratch and batch storage simultaneously below the fixed-memory boundary; nothing relies on unproved reuse. The exact current total is in range_ledger.py and proof.md.


# Appendix: FILTER_FREE.md

# Filter-free labels and complete load domination

Write R0=2^256, delta=R0-p and k for the canonical Montgomery x-key of an ideal nonidentity leaf. Set U(O)=empty. The hot arithmetic and its denominator/identity exception charge are unchanged.

## 1. Exact output set

The raw Montgomery square satisfies 0<=T<R0+p. The common-sum representation c lies in[0,R0), hence z=T-c lies in(-R0,R0+p), and h=floor(z/R0) is -1,0 or1. The assembly stores u=(z-h*p) modulo R0. Its canonical field answer is k=z modulo p.

- For h=-1 and z>=-p, u=z+p=k. If z<-p, then u=z+p+R0=k+delta, with k=z+2p>p-delta.
- For h=0, u=z. If z<p this is k; otherwise u=k+p and k<delta.
- For h=1, u=z-p lies in[delta,R0). It is k if below p; otherwise it is k+p with k<delta.

Thus every output belongs to U(k)={k}, with k+p added when k<delta and k+delta added when k>p-delta. Since delta<p/2, the exceptional bands are disjoint and |U(k)|<=2. Every extra label is in[p,R0), not arbitrary machine garbage. This covers the noncongruent negative-wrap case as well as redundant congruent residues. The emitter now outputs u without a representation test, subject only to its fully charged slot cap.

## 2. Witness exception loss

Let Bbad=[0,delta) union(p-delta,p). A true negation collision shares its canonical k. This set has at most2delta keys and4delta curvepoints, each with at most two negation-collision partners. For two independent cross-node uniform points, a collision with k in Bbad therefore has probability at most8delta/r^2<2^227/r^2. Summing over Gamma gives the retained Gamma*2^227/r^2 allowance. It is applied unconditionally, before any capacity/survival conditioning. No global S*2^-32 term is omitted: exceptional competitors are handled by the set-valued load bounds below, not excluded globally.

## 3. Single-cell bounds

For a fixed f-bit pattern let h be the number of its fixed positions in[224,256). Because p>R0-2^224, every extra label has all its top32bits set. The containing interval[R0-2^224,R0) has at most2^(224-f+h) integers satisfying the pattern. Each extra raw label has at most two inverse translated canonical keys, and each canonical key corresponds to at most two curvepoints. Extra point preimages are therefore at most4*2^(224-f+h).

For the35-bit slot masks and96-bit home/bucket mask, h=0. FOURIER_CAPACITY.md bounds canonical point counts. In the96-bit case its error is below4721808*2^128, while the new extra fibers add at most4*2^128. Their sum4721812*2^128 is below2^151, the spare1/512 fraction of2^160. For35bits, canonical error plus4*2^189 is also below2^221/512. Consequently both refined coefficients remain513/512.

For home70, the crude canonical count2*2^186 plus4*2^154 gives coefficient2+2^-30. The group144 and matched182 masks include sixteen fixed positions in[224,256), so their coefficient is2+2^-14. Divide these counts by r-1 where a free generator excludes one group point. These are mask-specific bounds, not claims about arbitrary f-bit subsets.

Use the fixed indicator A_C(P)=1[U(k(P)) intersects C]. It is a deterministic function of the ideal point, even when the actual representation depends on a whole batch. Every actual record in C is dominated by this indicator. Define all ideal variables under the original independent generator distribution; do not condition on clean execution or chosen raw representatives.

## 4. Pair loads and dependence

For a fixed opposite ideal point P', its set U(k(P')) has at most two labels. A free point P whose possible labels can match P' on182bits lies in a union of at most two fixed cells. Hence the matched-degree bound is M=2*KB182. Home and relay-cell factors use the single-cell KB bounds. With rho=(r/(r-1))^4, retain the source's conservative shared-generator cases with

h=rho*(S^2*KB70*M+(2S+20)*M),
g=2*rho*n0*S*M,
v=18*rho*((S/SIDE)*S*KB35*M+(2S+20)*M).

The witness is charged separately if bad; otherwise its actual key is its unique canonical label and its home/group are fixed by the witness generators. Thus its intruder rate is rho*(S+12)*KB144, without another factor2. Slot and bucket indicators remain fixed point predicates, so the source covariance argument applies: independent incidences/free Y cases, the one-excluded-point bound for a common X, invertible two-sign forms for a shared incidence, and the conservative60-partner allowance for witness-generator incidences. The phase slot uses qphi, while the persistent bucket uses all q1.

Actual partial/full matches can be false, including between two exceptional labels. They are nevertheless members of the ideal all-unordered-pair envelopes above. Caps and verification account for them. Removing records can create adjacent pairs, so no subset claim about an earlier adjacency list is used. Source labels, signs, generator dependencies and exact final verification are unchanged.

## 5. Work and evidence

For the retained filter-free reference emitter, deleting the two CMP/B.HS pairs and one threshold initialization changes its upper bound from33n+12 to29n+11. The current block emitter instead follows BLOCK_PROTOCOL.md. The arithmetic is unchanged. Its full-record envelope already charges up to two records for every incidence; formerly filtered outputs are not free new work. Exact numerical success and cost rows, including the enlarged h/g/v bounds, are recomputed in range_ledger.py. Boundary and differential tests explicitly assert membership in U(k), not equality of every raw output with k. These finite tests supplement this argument; they do not independently prove production success.


# Appendix: FOURIER_CAPACITY.md

# Fixed-bit capacity lemma

Let p be the protected P-256 field prime, r its prime curve order, R0=2^256, and rho=(R0 mod p)^(-1) mod p. For canonical Montgomery keyk, the number of affine curvepoints is1+chi((rho*k)^3-3rho*k+b). Define g(k) to be this quadratic character for0<=k<p and0 forp<=k<R0.

Primary reference for the mixed character bound: Todd Cochrane, Exponential sums modulo prime powers, Acta Arithmetica101.2(2002), introduction equation(1.2), author PDF https://www.math.k-state.edu/~cochrane/research/expprpowa.pdf . A nonzero additive frequency has polynomial degree1 and the multiplicative cubic degree3, characterorder2, giving3sqrt(p). At zero frequency Hasse gives2sqrt(p). The cubic is nonsingular and not a square.

An arbitrary realfrequency exponential over0..p-1 can be expanded in the p additive characters. The geometric coefficients satisfy min(1,1/(2p||theta-t/p||)). Counting at mosttwo frequencygridpoints per distanceannulus bounds their normalizedL1 sum by2+H_floor(p/2)<3+lnp; the deliberately looser4+2lnp is also safe. Hence every R0-DFT coefficient of zero-extendedg has magnitude at most3sqrt(p)(4+2lnp).

For an interval of lengthL in Z/q, its normalized FourierL1 norm is at most2+lnL: use|D_L(t/q)|<=min(L,q/(2|t|)) and the decreasing-sum integral. Fixing a bitblock[a,b) is a periodic interval of length2^a and period2^b; periodic lifting preserves normalizedL1. The lowest-bit fixedblock has norm exactly1. The complete bitmask is a POINTWISE PRODUCT of these indicators, so the Fourier convolution inequality bounds its normalizedL1 by the product of these norms. Do not use the false normalizedL1 product claim for pointwise Fourier products of AP convolutions.

For SIDE35 and d3=26, bucket bits are[0,35),[64,99),[128,154). There are2^160 total integerkeys in the cell. Its indicator has norm at most(2+64ln2)(2+128ln2). Pairing againstg gives pointcount <=2^160+3sqrt(p)(4+2lnp)(2+64ln2)(2+128ln2). Zero extension avoids wrap or proper-GAP assumptions. The count of canonical integers in the cell is at most2^160.

Safe rational constants: ln2<7/10, lnp<180, sqrt(p)<2^128,2+64ln2<47,2+128ln2<92. Error<3*364*47*92*2^128=4721808*2^128<2^151. Thus use KB96=(513/512)*2^160/(r-1). For either35-bit slot, norm<=47, error<3*364*47*2^128<2^144, far below2^221/512; the same513/512 coefficient is safe. Retain crude two-point fiber bounds for all other patterns, particularly144and182bits. This explicit lemma is used by proof.md; its argument and cited character-sum theorem are material evidence.

Conditional leaves with one free independent generator are uniform or uniform excluding one group point; denominatorr-1 is retained. Recheck every dependency/variance step and every positive capgap when applying the new KBs. Do not merelysquare a negative gap.


# Appendix: RECEIVER_PROTOCOL.md

# Two-phase direct-bucket receiver machine

This specifies the changed transport component. It is an abstract scalar AArch64 instruction-work construction, not a claim about measured network throughput. All processors' instructions and the stated delivery/in-flight work are charged. No direct-home shortcut or uncharged key oracle is used.

## 1. Contexts, channels and bounded state

Each algorithmic node has its main scalar context and one dedicated scalar receiver context sharing that node's coherent private memory. The receiver is resident, not invoked by a per-packet software dispatcher. Its code and register parameters are fixed before hop2 begins and persist across both phases. It executes independently while the main context stalls on a send or waits for completion. Both contexts' work is included.

Final-cursor, barrier and phase4 messages retain the original addressed delivery handler. Hop1 is replaced by BLOCK_PROTOCOL.md, including its additional fixed payload windows and completion contract. Hop2 payloads have a separate fixed row-payload endpoint/port bound to the receiver window; control messages cannot enter it. This is a fixed channel distinction, not a content-dependent classification operation. The extra context, stack,48-byte packet window, valid/ack state and payload-port state require less than2KiB per node, included in the fixed-region reservation. The existing sender ports and their in-flight state remain charged separately.

Use the source blocking-port convention: a sender's completing STP waits if its bounded port/window cannot accept another message; the receiver's first LDP waits until a complete packet is atomically published. A waiting instruction counts once, not repeatedly as a poll. Arbitration is fair, backpressured and the same transport mechanism as the source delivery construction. There is one bounded receiver window, no hidden unbounded software queue, and no software context-switch loop. The published packet remains immutable until acknowledgment. Packet publication and the blocking device load have acquire ordering; the window is an explicitly ordered device endpoint, not an ordinary unsynchronized shared-memory buffer.

## 2. Receiver body and acknowledgment

Resident inputs are the packet address, d3 mask, cursor/end-table base, d4 mask, counter base, acknowledgment endpoint, and constants1,32,-32. The body reads the four payload words in two LDPs. It computes d3, addresses the interleaved cursor/end pair, reserves32bytes with LDADD, loads the end and compares. On a full bucket it undoes the increment. Otherwise it stores the two payload pairs, computes d4 and increments its counter.

Both exits execute STLR xzr,[ack] after those writes, then branch to the blocking first load. STLR is a base scalar release store with an unindexed address; its encoding is restricted in the interpreter. The finite callable test uses RET in place of that final branch, with the same one-instruction price. The accepted path costs15 instructions; a full-bucket drop costs11. Charge15 for every potential record, including ones that in fact never reach hop2 because of earlier drops.

The explicit release acknowledgment is charged in addition to the source's two per-message in-flight updates. In-flight decrement occurs only after the acknowledgment, after insertion/drop and the d4 update are visible. Delivery completion means processing completion, not inbox arrival. The main context's completion wait observes that acknowledgment with acquire semantics; its DSB/barrier cannot pass merely because the packet reached the window. This release/completion contract is a material part of the machine definition.

There is no per-message call, argument setup, register save/restore or polling: parameters stay in the dedicated context. Startup and public parameter/address setup are bounded by the separate2^14 per-node phase/context-control allowance. A final issued blocking receive after the last packet costs at most one additional instruction and is included there; with no further messages it executes no further work. Program, context and port initialization also remain within the charged2^20 node initialization row. These bounds do not excuse an implementation that adds a software scheduler; such an implementation would require its own extra ledger.

## 3. Phases and progress

Partition absolute X indices into two contiguous phases, each containing at mostceil(q1/2) indices. Keep absolute k in w2. For each phase initialize the block windows and local/remote cursor tables as in BLOCK_PROTOCOL.md; emit full and partial blocks and final valid cursors; drain and barrier. Then h2s streams those local slots through the dedicated row-payload port. Extract destination d2 before BFI replaces w3bits0..34 with original source column j. Sign/l and w2 are untouched.

The receiver never sends another payload or waits on main-context work while processing a packet: it only performs bounded local operations and acknowledges. Fair execution therefore drains finite sender workloads even if main contexts block on their ports. After h2s, each sender waits for every acknowledgment and all nodes barrier. Only then may hop1 slot storage be reused. Buckets and d4 counts persist through both phases. Prefix4, full-record pass4 and scanning start only after the second phase fully drains.

The two first-hop cursor initializations and finishing sweeps are charged. h2s costs7n+2floor(n/8)+11*SIDE+2 per phase. Bucket cursor/end initialization costs5*2^26+1. After the phases, full-record pass4 costs6n+2floor(n/8)+11*2^26+1. Prefix and scanner costs are unchanged. Sixty-four global barriers at1024instructions/node and2^14 additional phase/context-control instructions per node conservatively cover the finite public phase transitions, arguments and completion setup. Receiver record work is not hidden inside those fixed allowances.

## 4. Addresses and resource obligations

Reserve[0,2^52) for fixed memory; add the d3 end table and2KiB extra receiver state explicitly. Place phase receive slots, tight32*CAP_3 full-record buckets, then the8-byte pointer array in disjoint increasing regions. range_ledger.py and fixed_layout.py check the exact endpoint below2^64, fixed reservation, positive cap gaps and attempted hop1 cursor bound. The bucket increment/rollback cannot wrap because its end is strictly below the subsequent pointer-array region and the final address limit. No second-hop R2 receive array or final-cursor table is needed. The bound does not rely on reusing dead storage.

Finite tests run actual emit, first-hop delivery, final cursors, sender, receiver acknowledgment, prefix, scatter and scanner over sixteen nodes and two phases, with filtered-free high raw labels, slot/bucket caps and independent reference contents/counts. They establish the stated finite behavior and instruction counts, not a physical asynchronous implementation. The protocol and all ideal-machine assumptions remain explicit review obligations.


# Appendix: COLD_ACCOUNTING.md

# Conservative complete cold allowances

This analytic construction refines unused reserves, not measured production throughput. It uses the exact P-256 field primitives already supplied and checked: multiply116, square92, add23, subtract19, inverse39479, including each primitive's return. Scalars and field values are four64-bit words. No arbitrary-width operation is a unit-cost primitive.

## Affine add/double, including exceptional points

Cold temporary frames represent a point by two canonical Montgomery field coordinates and one infinity flag. Persistent source_X stores infinity in-band as(0,0); nonidentity points cannot have y=0 in the protected odd-order group. Encoding/decoding this sentinel is charged in the existing construction/row control allowances. Check the flags, then coordinate equality. Copy the other point for an identity input; return infinity for opposite points or doubling with y=0. Otherwise use the ordinary slope formulas. Addition uses two differences, one inverse, two multiplies, one square and the remaining coordinate differences. Doubling uses slope=(3*x^2-3)/(2*y): two squares, two multiplies, one inverse and at most eight additions/subtractions. The constants1 and3 are stored in Montgomery representation. Thus either nontrivial path uses at most16 field calls, of which at most one is an inverse. The other fifteen calls may all be charged at116, deliberately exceeding the actual mix.

A direct wrapper places at most16 temporary fields in a fixed frame, keeps live pointers/counters in frame slots across calls, and uses no dynamic allocation. Charge16 non-field instructions per call (pointer loads/address formation, BL and live-state management),128 for flag/coordinate comparisons, branches and copies, and128 for frame/constants/result handling. Four-word equality/zero tests use bounded loads, XOR/OR and branches; coordinates and frame offsets are below4096. These512 wrapper instructions are a conservative explicit construction allowance. Required field loads/stores and RET are already in each field primitive count and are not omitted.

One complete affine operation therefore costs at most39479+15*116+512=41731. The41731 allowance also applies to exceptional branches, which skip arithmetic. Storage is bounded by the fixed frame and does not alter the published node-memory bound.

## Scalar multiplication and generation

Left-to-right binary multiplication over at most256 scalar bits performs at most255 doublings and256 additions, starting from infinity. Allow8192 further instructions: 16 per bit for extraction/branch/loop state (4096 over 256 bits, the same per-bit rate used for modular products below) plus 4096 for initialization, conversion/copies and output, separately from the point wrappers. Its work is at most511*41731+8192=21332733.

Compute R_g=[b]Q once per node for this target, retaining original b and R_g. Each X then requires one [a_k]P and one affine addition with R_g; each Y requires one [c_l]P. Charge21332733+41731+16384=21390848 for each of the q1 X items, and21332733+16384=21349117 for each of the q2 Y items and the one R_g, which need scalar multiplication and bookkeeping but no affine addition. It does not amortize across targets, and row-specific signed b'_k never replaces original b during construction. Sampling is separate. Header/canonical conversions and fixed bookkeeping stay within16384 per item: at most 64 loads/stores for two coordinates, two scalar labels and headers, at most two field products for conversions, and at most 512 address/control instructions; the true total is below 2000.

## Scalar-order arithmetic

For modulus r, canonical four-word modular addition can add into five words, subtract r extended by a zero word, and select the canonical result. Even charging all loads/stores and control,64 scalar instructions suffice: twelve input/modulus word loads, five additions, five subtractions, four conditional selects, four stores and at most16 pointer/control instructions total46; the64 allowance is larger. The same construction doubles a canonical value. The fifth word handles overflow; inputs below r give a result below2r and one subtraction suffices.

A256-step shift-and-add modular product doubles the accumulator each step and conditionally adds the multiplicand. Charge16 additional instructions per bit for extraction/branch/loop state and64 fixed instructions. Worst work is64+256*(2*64+16)=36928, strictly below2^16. Each addition is canonical, so induction gives the exact product modulo r. Binary exponentiation to r-2 uses at most512 such products. Including32 controller instructions per product and512 fixed instructions costs at most512*(36928+32)+512=18924032, below2^25. Fermat applies because the protected scalar group has prime order r. A zero solve denominator is rejected before inversion.

## Candidate resolution

The source-service and four explicit bounded-routing rows remain separately charged in full. The resolver receives the stored generator points and scalar labels, reconstructs its two leaf points with at most two affine operations, tests the negation relation and nonzero scalar denominator, solves for x modulo r and verifies [x]P=Q. Allocate:
- scalar inverse:18924032;
- verification scalar multiplication:21332733;
- two affine reconstructions:83462;
- all remaining scalar operations, comparisons, canonical output, local message preparation and bookkeeping:2^18=262144.

The last row can cover four full modular products at36928 each, sixteen modular additions at64 each, and more than100000 instructions for the remaining fixed-width control/copy/output operations. The solve itself needs only one product after the inverse and a few modular differences/sign changes. Transport/service counts are not reduced or charged implicitly to this row.

The sum is40602371, below2^26. Charge40602371 for every allowed candidate, including false partial matches and failed verifications. Early exits cannot increase this bound. No claim depends on a prior submission's acceptance, or on the sharper54413312 resolution allowance in public submission1d6d2b74.

The exact integer identities are reproduced in range_ledger.py; scalar_check.py independently checks word-level products, carry-aware halving and inversion. Canonical cold field routines are p256_fields_reference.s, distinct from the weak hot routines; a separately namespaced copy is used for cold operations and inversion. Expanded program delivery remains within2^20 per node. Four coordinate decodings, sign/equality checks and exceptional dispatch are included in the262144 miscellaneous resolver allowance; no second charge is silently removed. Those finite tests supplement, not replace, the word-level and loop-invariant arguments above. A full protected package replay and external exact-version review are still required before publication.


# Appendix: RANGE_DRIVER.md

# Range-tree driver decomposition

2026-11-02. Analytical construction with tested native merge/search helpers. Python finite orchestration is not itself a counted production implementation.

Build an array of72-byte records (x,y,original ID), stable-sort by(y,ID), and use its fixed balanced binary tree. At each tree node store the same records sorted by(x,ID). Build node lists bottom-up by merging child lists, rather than independently sorting every node. This avoids the extra logarithm of the naive construction. Keep IDs unchanged and original scalar tables separately. Each level stores exactly n records; there are29 levels for n=2^28. Total tree records72*n*29; a64-byte descriptor per node and two72*n sorting buffers bring the explicit allocation below2^40 bytes. fixed_layout.py inserts these regions without overlap into the operative address map.

A scalar-instruction budget is128 per merged record plus64 per merge. A fixed five-word lexicographic comparison, selected72-byte record copy and pointer/count branches can fit this allowance, The supplied merge.s has at most64instructions per record plus8fixed, leaving additional space for caller setup in the128/64 allowance. A y sort and an x merge tree each cost at most128*28*n+64*n. Add64 per tree node and64 per leaf copy: total7488*n-64, within8192*n. No free sorting oracle, per-target amortization or construction of advice is assumed.

At query rank r, a single root-to-leaf split path leaves at most29 disjoint fully contained tree nodes covering all leaves, with a fixed y-side bit on each. In each node list, binary searches at X.x and2^256-X.x split it into at most3 segments. If X.x=0, the latter threshold is a sentinel beyond all field coordinates, not a wrapped zero. Thus at most87 descriptors suffice. Each record occurs exactly once because node intervals are disjoint and each local split is a partition. Empty segments are omitted.

There are at most29*(1+2*29)=1711 coordinate comparisons: one global rank search and two searches per cover node. With a 32-instruction comparison/search-step budget (search.s uses at most16 per step plus10 fixed),256 per segment for descriptor writes/argument setup/call overhead,64 per visited path/cover node, and4096 fixed initialization, the per-row driver bound is84768, below2^17. Each64-byte descriptor stores node-list base, start/count, batch offset, x/y direction, carry and original-ID cursor. At most48 memory/address instructions construct it; caller spills/restores and arguments cost at most128, and all three routine fixed terms total37. These fit256 per descriptor. Merge callers use at most64 instructions for public loop/address calculations, six argument moves, BL and bookkeeping. Tree nodes use at most64 instructions for eight field stores, address arithmetic and traversal. These reserves are deliberately separate from per-item arithmetic; call overhead already charged elsewhere must not be subtracted twice. The stated constants give the bound directly; BLOCK_PROTOCOL.md conservatively allows the additional261 per row.

Forward calls consume descriptors in order while carrying the prefix across segment boundaries. Reverse calls traverse the same descriptors backwards with the one row inverse. A zero denominator aborts the entire row before any key emission. Carry count is the sum of lengths of high-carry segments; the global cap must be checked before emission, and cannot reset per phase. Slot admission and source-ID recovery are separate obligations.

The tree can be generated by explicit bounded loops over levels and contiguous runs, with two reusable y-sort buffers and a separately allocated list for each tree node. No library sort or uncharged allocator is used. Child/node offsets are computed by multiplication/addition from level/range indices; store at most eight64-bit descriptor words. Stack depth is at most29 and fits the context reservation. Argument setup for x-key search uses node base; for y-rank search it uses the y-sorted array base+32 with the same72-byte stride. Call the second x threshold search only for X.x!=0; otherwise its result is list length. All addresses and multiplications fit the explicit layout. The instruction budgets above include these operations, including empty searches/segments, rather than relying on finite timing. The full attack proof separately treats carry caps, phase4 and source recovery.

Control suballocation, relative to the4096-byte-aligned65536-byte contexts_ports_barriers base: query descriptors[0,5568); traversal/call stack[6144,14336); resident contexts and endpoints[14336,22528); source scratch[22528,23296); eight isolated2048-byte barrier-word granules[24576,40960). The gaps are reserved padding. Barrier words are2048-byte aligned and no unrelated stores share their granules. Epoch reuse follows the drain/barrier protocol. All regions are disjoint and end below65536. The caller retains original IDs and treats sentinel(0,0) X as a defective row before search.


# Appendix: RESOLUTION_PROTOCOL.md

# Explicit bounded resolution routing

This protocol replaces the earlier count/offset exchange completely. It does not assume that payload conservation establishes grouping or allocation work. There is no count exchange, sender-side sorting, offset reply, per-destination linked list, or undisclosed allocator.

## 1. Record formats and endpoints

A phase-4 packet occupies an 80-byte endpoint window: destination index and unused zero address word (16 bytes), followed by an eight-word, 64-byte record. Publication occurs only when all five STPs have completed. The first receiving LDP blocks until the complete record is published; the window remains unchanged until STLR releases it. A busy sender store blocks. In-flight state is incremented once on acceptance and decremented only after the receiver's release. Both updates are charged.

All node coordinates are 35-bit integers. A request record contains:

- word 0: source row | (X index k << 35);
- word 1: source column | (leaf sign c << 35) | (Y index l << 36);
- word 2: home row | (candidate index j << 35);
- word 3: home column | (candidate side << 35);
- words 4 through 7: zero.

The fields fit: k < q1 < 2^29, l < 2^28, c is 0 or 1, j < CAND_CAP = 2^24, and side is 0 or 1. These are produced from already bounded leaf and scanner records, not untrusted network input. Each request is uniquely identified by (home row, home column, j, side).

A reply record contains:

- word 0: home row | (j << 35);
- word 1: home column | (side << 35) | (chunk << 36);
- word 2: valid flag 1, separate from every routing field;
- word 3: zero;
- words 4 through 7: one 256-bit data chunk.

There are exactly eight chunks numbered 0 through 7. Chunks 0 through 6 carry X.x, X.y, Y.x, Y.y, a'_k, b'_k and c'_l; chunk 7 is zero padding. The leaf type remains in the home's candidate record, so it need not be retransmitted.

Each phase has a distinct input endpoint and resident receiving context. All four contexts are initialized once. They are not reconfigured while blocked. Phases 1 and 3 have distinct cursor words although they reuse the same relay storage after a drain. The phase-2 source cursor and phase-4 indexed-home handler are separate. The earlier leaf receiver has its own endpoint. Ordinary 48-byte fixed-address delivery remains a different message class; it does not silently acquire a new ACK requirement.

Every context performs the blocking payload load BEFORE reserving a record. Thus an idle receiver cannot leave a phantom cursor advance after its last packet. After completion, its next blocked load performs no reservation or further work; that one attempted load and all context setup are included in the fixed control allowance. The algorithm's output/failure decision is made after all required data drains, without waiting for another packet at an idle context.

## 2. Four phases and buffer lifetimes

At a home, `make_requests` reads each candidate's two leaf-record pointers and their source metadata, then writes two 64-byte requests. It costs exactly 17 instructions per candidate plus 3. The output fits the existing 64*REQ_CAP region because 2*CAND_CAP = 2^25 < REQ_CAP = 2^29. This region is not the relay region.

1. Request hop 1: home (h_r,h_c) sends each request to relay (g_r,h_c), using word 0's low 35 bits. The relay appends to its 64*RELAY_CAP region.
2. Request hop 2: after all hop-1 senders drain and a barrier completes, each relay scans its accepted records once and sends to source (g_r,g_c), using word 1's low 35 bits. The source appends to the 64*REQ_CAP region. The old home request output is dead by this phase and may be overwritten. Relay and source regions are disjoint.
3. Reply hop 1: after phase 2 drains, each source retains its accepted requests while processing them sequentially. It gathers the seven requested 256-bit values, builds eight tagged replies in a 512-byte temporary and sends them to relay (h_r,g_c). The relay region from phase 1 is now dead and reused, with its separate phase-3 cursor initialized at the base.
4. Reply hop 2: after phase 3 drains, each relay scans each accepted reply exactly once and sends it to home (h_r,h_c). Home placement is direct indexed storage, not append allocation.

Each phase ends with all its sender ports drained and a global barrier before any receiving buffer is read or reused. No source discards its request before all eight replies have been copied into outgoing packet windows. The 512-byte output temporary can then be reused because the endpoint holds its own packet bytes, not a pointer into that temporary. At most eight replies per request are ever produced. Forwarding does not duplicate records.

The append receiver accepts the first CAP records in arrival order. It reserves 64 bytes by LDADD, compares the returned address with the end, and either stores all eight words or rolls the cursor back. Both paths release the packet and continue. Excess records are dropped individually; no whole node is aborted. The cursor remains at or below the end between packets, and its temporary end+64 address is within the reserved fixed-space padding.

## 3. Reply association and completeness

At home the receiving kernel decodes j, side and chunk from the first two words and computes

slot = ((2*j + side)*8 + chunk),
address = reply_base + 64*slot.

This is injective on the valid tag tuple. Since each request is produced once and each accepted request produces each chunk once, no two legitimate writes collide. Two candidate leaves from the same source still differ by side. Different candidates requesting the same source leaf still differ by j or home.

Before any request is emitted, all 16*CAND_CAP flag words are zeroed. Only reply records write their flag to 1. After the final drain/barrier, the `complete` kernel checks all 16 flags for a candidate before any payload is used. A missing request or any missing reply chunk therefore drops that candidate. A final chunk cannot falsely certify the earlier chunks. All flags are from this run; there is no epoch reuse without initialization. Tags have bounded widths, so the largest slot is 2^32-1 and the home reply region is exactly 16*CAND_CAP*64 bytes.

The home reconstructs the two leaves from X,Y and the retained leaf types, rejects either reconstructed identity leaf, tests equality up to negation, selects a matching sign with nonzero denominator, performs at most one scalar inversion and one target verification for that candidate, and otherwise rejects. For valid nonidentity witness leaves the matching sign is unique; exceptional or invalid candidate data may be rejected without trying additional inversions. It verifies [x]P=Q before output. The transport returns precisely the seven values stored at the requested source indices. Source X stride is64 bytes, with x and y at offsets 0 and 32. Source Y stride is64 bytes, with x at0 and y at32. Scalar arrays have 32-byte entries; b'_k is fetched from its row-specific32-byte scalar slot; original b is never substituted.

## 4. Executed instruction counts

The supplied phase4.s and resolution_check.py execute the following scalar routines:

- `send64`: 14 per record plus 2. It loads the chosen routing word, masks the destination, writes the routing header, copies all eight payload words, advances the input pointer and loops.
- `append64`: 13 on acceptance, 7 on a full-buffer drop. This includes acquire-before-reserve, cursor/cap handling, payload copies or rollback, release and the loop return/branch.
- `home64`: 16 per record, including tag decoding, indexed address calculation, all payload/flag stores and release.
- `clear_flags`: 3 per flag plus 2, using a postindexed store, SUBS and B.NE.
- `complete`: at most 67 per candidate, including all 16 loads and checks and its result.
- `make_requests`: 17 per candidate plus 3, including both leaf metadata loads and every request store.
- `gather224`: 31 per served request, including all seven 256-bit reads, corresponding writes and padding initialization.
- `build_replies`: 84 per request for all eight complete tagged records.

The finite code returns to a sentinel for counting. A production resident receiver branches back to its blocking first load instead of RET; this is the same one instruction. The model's blocking/acquire/release behavior is explicit, not something the finite emulator proves about real hardware.

For source service, decoding and checking k,l,type and locating X,Y,a,b,c require at most 64 instructions (shifts/masks, comparisons, fixed-stride MUL/ADD and constant loads). Gathering takes 31 and constructing replies 84. Allow another 256 for argument setup, spilled loop state, request iteration and scratch management. Total 435 is below the retained 1024/request allowance. The output transmission itself is separately charged below. Source processing uses only 256 bytes of gather scratch and 512 bytes of reply scratch, not an unallocated 8*REQ_CAP reply array.

## 5. Complete phase-4 work bound

For each phase, every input record is visited once by `send64`. Intermediate phases read contiguous append arrays; no sort or directory lookup is needed. Request creation and final completeness checking are explicit routines. Per record, the executed sender (14), receiving path (at most 16) and in-flight updates (2) total at most 32. Charge 256 per record, which also dominates per-request send-call fixed costs, record iteration and any fixed-width address/tag moves not already charged. Empty calls and context initialization are in fixed work.

The fixed allocation per node per phase is 2^37 instructions. This covers, conservatively, request creation (17*CAND_CAP+3), flag initialization (48*CAND_CAP+2), all candidate completeness checks (67*CAND_CAP), cursor resets, scalar argument/context initialization, buffer-length calculations, state spills and phase exits. Even 512*CAND_CAP is 2^33, within the 2^37 allocation; the explicitly listed per-candidate terms sum to 132*CAND_CAP, leaving substantial fixed-control space. The separate source-service allowance is not used twice to justify missing transport work. A final constant-degree tree reduction selects one verified scalar or reports failure. The parent of (i,j) with i>0 is ((i-1)//2,j); for i=0,j>0 it is (0,(j-1)//2). Each node has at most four children with public child-slot addresses in a separate16384-byte output-mailbox region. Each of four slots has a2048-byte data granule and a distinct2048-byte flag granule, so data writes cannot wake a flag wait. All four flags are initialized to zero at node setup BEFORE the first global barrier or any possible publication, never cleared lazily when waiting begins. They are disjoint from all earlier barrier words. Each child sends one 64-byte completion record containing its selected scalar and a final flag (1=found,2=empty). send_completion writes the explicit parent axis index and public child-slot address in its routing header, then all eight payload words, in10instructions; the ordinary send64 routine with zero address is not reused for this job. For a row child i>0, slot=(i-1)&1; for a column child i=0,j>0, slot=2+((j-1)&1). Choose OUTPUT_BASE=2^20, exactly after the2^20-byte code/stack reservation; it is4096-byte aligned. Reserve the following16384bytes before all later fixed regions, so no additional alignment padding is needed. The public destination address is OUTPUT_BASE+4096*slot. Parent and slot use shifts/parity from the public row/column, at most64setup instructions. publish64 stores the seven data words and release-publishes the last packet word at child_slot+2048 in12instructions. The parent waits on that isolated flag address. Each flag has exactly one post-initialization publisher. Parent waits acquire the flag before reading the scalar; every child writes once, so the existing bounded single-wakeup wait construction applies. Each node chooses at most one valid result and forwards once. Even4096instructions/node covers four child reads/selections/waits, one send/delivery and fixed addressing. This is included in the unused fixed allowance, not the global OUTPUT constant:132*CAND_CAP+2^20<2^37, with2^20covering this reduction and all remaining context/control constants. Root serialization alone uses OUTPUT. Thus many simultaneous successful candidates do not cause an uncharged all-to-root fan-in.

Let M be the total number of attempted records in a phase across all nodes, including those later dropped. Requests have M <= N*2*CAND_CAP; replies have M <= N*8*REQ_CAP. Drops can only reduce subsequent phase inputs. Therefore each phase costs at most

N*2^37 + 256*M <= N*(2^37 + 256*2^32) = N*9*2^37 < N*2^41.

The existing ledger's four EXCH rows, each N*2^41, thus bound this complete replacement. This is not an assertion that conservation alone establishes work: the routines and fixed-cost calculation above establish the work per record and per node; conservation then bounds how many times they execute. No per-candidate arithmetic slack is needed to cover this protocol.

## 6. Memory and success-loss compatibility

Existing fixed-space regions already include 64*RELAY_CAP relay bytes, 64*REQ_CAP source/request-output bytes, and 16*CAND_CAP*64 home reply bytes. The candidate-pointer list is separate. The relay is reused only after phase 2 drains. Source/request-output reuse occurs only after phase 1 drains. Home flags occupy their own metadata words inside reply records.

Reserve 8192 bytes, replacing the earlier 2048, for the leaf and four phase-4 resident contexts, their separate 80-byte windows and cursor/ACK words, phase-specific sender windows, 768 bytes of source scratch, and bounded stack/control state. Even six 31-register contexts (including fixed completion delivery), twelve 80-byte endpoint windows,768bytes of scratch and1024bytes of stack/control fit. Together with the separate16384-byte output-mailbox region, this increases the explicit F sum by22528; fixed_layout.py explicitly allocates these and a larger64KiB context/control region below2^52. The reserved bulk endpoint E and total memory bound do not change.

Candidate cap2^24 and source cap2^29 are retained. Allocate2^34 relay records per phase. For each witness endpoint put B=rho*((S/SIDE)*S*KB35*M+(2S+20)*M), where M=2*KB182 and SIDE=2^35. A request traverses relay(source row,home column); its two endpoint requests have envelope2B. A reply traverses relay(home row,source column); its sixteen candidate chunks have envelope16B. Each direction uses the corresponding35-bit marginal bound, identical by the symmetric bit-block bounds in FOURIER_CAPACITY.md/FILTER_FREE.md. The shared-generator cases in proof.md apply separately to these fixed marginals. Both dominate ALL unordered matching pairs, including pairs newly adjacent after deletions; no adjacency-subset assertion is used. Their union has conditional nonfixed-load expectation at most18B=v.

Four anchored original generators fix at most four sign-leaves (two incidences), hence at most six fully fixed unordered pairs. Reserve108=6*18 possible request/reply records, not merely the historical16. Markov therefore bounds either required-relay failure by v/(2^34-108), which is at most the retained numerical v/(2^30-16). The stricter diagnostic denominator remains safe; it is not the actual buffer size. Drops can only remove traffic. Missing home flags arise only from one of these capped request/reply paths. No additional completion probability or whole-node abort is introduced. Source-request and candidate-load bounds remain as in proof.md with their existing explicit reservations. Full work is bounded by native per-record routines plus fixed initialization and only THEN by packet-count conservation.

## 7. Target-specific resolver and evidence

COLD_ACCOUNTING.md supplies the complete40602371 P256 bound, including18924032 scalar inverse,21332733 verification multiplication, two41731 affine reconstructions and262144 miscellaneous control/scalar/coordinate-decode work. Payload coordinates are canonical Montgomery values; decode them with four field products or execute the same canonical Montgomery affine routines. Signed coefficients use only negation, additions and one final product after inversion. Zero denominators and reconstructed identities are rejected before inversion; false candidates still receive the full charge.

resolution_check.py executes this protocol with35-bit nodes, original IDs and distinct row-specific b values. It tests21small capped cases, six normalized actual P256 recoveries, partial replies, maximum home tags and nine output-tree configurations. Scalar word checks cover fifth-word carry in halving. The asynchronous blocking/acquire/release and completion machine remains an explicit analytical convention, not a hardware-timing claim. Every context's instructions and all storage are charged.


# Appendix: CARRY_MEAN.md

# Carry mean and collision conditioning: proposed completion

2026-11-02, research argument requiring review, not a submitted proof.

Source theorem: Todd Cochrane, Exponential Sums Modulo Prime Powers, equation(1.2), author PDF https://www.math.k-state.edu/~cochrane/research/expprpowa.pdf . For a nonzero additive linear phase and a quadratic character of a squarefree cubic, it gives the bound3sqrt(p); the degree/order hypotheses are1,3,2. This is the only external theorem used below.

For a uniformly random nonidentity P256 point, its canonical Montgomery x-coordinate K has mass (1+chi(g(k)))/(r-1), where g is the nonsingular curve cubic after the invertible Montgomery rescaling. Its nonzero additive Fourier coefficients are bounded by3sqrt(p)/(r-1). For an interval, the normalized absolute Fourier coefficient sum away from zero is at most H_floor(p/2)<=1+ln(p): pair frequencies t and p-t, and use the geometric-sum bound p/(2t). Fourier inversion therefore bounds every interval discrepancy from uniform Fp by epsilon=3sqrt(p)(1+ln(p))/(r-1). With sqrt(p)<2^128 and ln(p)<180, epsilon<543*2^128/(r-1)<2^-117.

For independent K,L of this distribution, replace L by a uniform integer U in[0,p) inside the event K+L>=R0, R0=2^256. The conditional event is an interval, costing at most epsilon. The resulting function of K is monotone in[0,1], so summation by parts and the CDF discrepancy replace K by an independent uniform V at cost at most epsilon. The uniform probability is (2p-1-R0)(2p-R0)/(2p²)<1/2. Hence carry probability<=1/2+2^-116. Sign normalization does not change x. Identity draws can be assigned zero carry and handled by their existing abort loss; conditioning all generators nonidentity instead would alter independence and is not assumed here.

For collision conditioning, first fix the target and each node's shared b scalar. Uniform independent a_i makes X_i=a_iP+bQ independent uniform group points even after this fixing; c_j gives independent Y_j. A cross-node witness depends on four generator values. Fix those generator values as well; every unmentioned variable remains independent with its original distribution. This is the setting of the anchored-tail lemma in the separate carry-conditioning audit. Do not condition on the entire execution avoiding every other abort. Instead charge such abort events separately. At most four anchored row and four anchored column indices is a conservative bound per node; all incident edges contribute at most B=4(q+u).

Thus the candidate's unconditional variance denominator slack² must conservatively become(slack-B)² for witness-conditioned tail accounting, unless a tighter conditional analysis is supplied. The exact ledger now incorporates this anchored correction; its conditional total is not a submitted score. This proof must be checked together with the exact leaf/collision event definition, zero-target path, scalar rejection caps and existing Poisson loss ledger. No claim that this note alone establishes the full attack.

Carry is indexed by original unsorted generator identities. The total is invariant under Y permutation and sign normalization; never condition on a random sorted rank. Capped scalar sampling is coupled to ideal independent draws, and exhaustion of eight attempts is charged separately by the global failure term.


# Appendix: CARRY_TAIL.md

# Witness-conditioned carry tail

Fix the target, shared node scalars and the witness's original generator values as in CARRY_MEAN.md. Remaining generators are independent. Delete every incidence touching an anchored row/column; at most B=4(q1+q2) such incidences exist per node and their carry contribution is at most B deterministically. The remaining carry indicators with disjoint endpoints are independent; any covariance is at most1/4 by Cauchy-Schwarz. Counting ordered pairs that share an endpoint (including diagonal pairs) bounds variance by q1*q2*(q1+q2-1)/4. Their total mean is at most q1*q2*(1/2+2^-116).

If cap=ceil(mean)+2^50 and cap-mean>B, Chebyshev gives conditional failure at most var/(cap-mean-B)^2. This bound is uniform in anchored values, so averaging over their collision event preserves it. Only original indices are anchored; sorting and sign normalization leave total carry invariant. Do not condition on global cleanliness. Identity X contributes zero in the ideal carry array and its independent abort loss is separately charged. A meter failure stops further leaf generation only; all receiver/service/control phases continue as specified in proof.md. Two witness endpoints yield the mu*(1+1/r)*2*d_m loss term. The exact ledger checks the positive gap and charges all admitted corrections in full.


# Appendix: GENERATED_PROGRAM.md

# Fixed generated arithmetic program

The four specialized backward variants and forward entries are defined by range_witness.generated_assembly from the retained reference templates and fixed literal transformations. The generated UTF-8 assembly is pinned to SHA256 7bf0bedffd2eb53652489011a8f2be7e840dccf911388a56d9c0a4364114f20a and is byte-identical to v15 range_batch.s. This removes a redundant expanded source listing from this add-only version, not any instruction from the attack. Every arithmetic experiment loads and executes that generated listing. The public fixed program is specified by this generator; it depends on no challenge or searched advice. Expanded code is delivered/initialized under the existing per-node2^20 program allowance and code/stack memory reservation, not assumed to occupy compressed bytes at execution. Independent finite behavior, instruction counts, source transformation and exact hash checks remain. This is not a claim that a hash proves correctness.


# Appendix: BLOCK_PROTOCOL.md

## Fixed windows, complete delivery and progress

Each sender has one header port and a fixed32768-byte payload window for each destination d, at payload_base+32768*d+32. A header is (destination, destination slot address, fixed payload-window address,1024). The emitted address always equals that destination's allocated window; it is not an arbitrary remote-pointer service. An active delivery context maps only this immutable header and that fixed payload window into source space. All payload stores use destination-local hop1 memory. Fixed source windows lie below2^51; destination slots start at2^52, so their numerical ranges do not alias in this mapping.

The completing header store waits until its own transaction has been completely copied and acknowledged, not merely enqueued. The destination's dedicated scalar context is independently and fairly scheduled; its parameters stay resident. The preinitialized receiver context stores its fixed ACK-endpoint address at header+32; this is context metadata, not a fifth transmitted header word. It reads a complete header, copies1024records using explicit LDP/STP instructions, and issues a native STLR1 to its fixed device ACK endpoint. That ACK is an event consuming the single current transaction, never a predicate satisfied by a stale stored1. No new transaction may replace the active mapping before that ACK. The sender can reuse the payload after ACK because the handler performs no further payload read. The trailing callable RET has the same cost as the resident branch back to its blocking header load.

There is one bounded current window per receiver and one outstanding header per sender. Arbitration is fair and backpressured; no unbounded software queue, polling, epoch arithmetic or software dispatcher is hidden. The receiver never waits on main-context work or another payload before acknowledging. Its fixed finite work therefore drains finite workloads even when producers block. Publication/acquire, payload writes/release ACK and sender completion are ordered device operations. Two charged in-flight updates bracket the transaction, with decrement after ACK. Existing phase drains/barriers precede any window/cursor reuse or hop1 consumer scan. These are explicit abstract-machine assumptions, not measured physical-network performance.

## Emission, partial blocks and prefix preservation

block.s initializes local cursors to fixed window starts and remote cursors to each sender's slot base. Each record is written with two STPs; only the last of1024 aligned32-byte slots has bits5..14 zero. On that slot, reset the local cursor, reserve32768 remote bytes, compare against the slot end and publish only if within capacity. A full but dropped block still pays the sender flush. Local buffers persist across all row descriptors in the phase.

At phase end block_finish computes the remaining valid byte count. It publishes one padded full block if the remainder is nonzero and the remote cursor is below the slot end. Padding is initialized before reads. Then it sends the clamped valid cursor min(remote+remainder,slot_end) by the ordinary final-cursor port. The completing block store already waits for ACK before this cursor is sent. CAP1=25*2^18 is divisible by1024. Therefore every admitted full block lies wholly within its slot, a final padded block remains in that slot, and the valid cursor excludes padding. In sender order these are exactly the first min(number_emitted,CAP1) records. After the drain, h2s traverses the same source-slot order as before. Record words and original IDs are unchanged. All prior capacity, intruder, pair-load and success bounds consequently apply without assuming a new distribution.

## Complete instruction and initialization bound

The native per-incidence emitter base is27 rather than29. Its fixed13 rather than10 costs at most3*87=261 extra per row; the previous explicit84768 driver bound becomes85029, charged at that explicit bound. A full flush costs at most9 sender instructions: mov x12,#32736 plus sub rewind, then the seven unchanged descriptor instructions; x12 is dead there because the record tag was already stored and is recomputed in the next half, and both new immediates are single AArch64 encodings. Native block_dlv costs4359 including ACK (3 prologue, 128 loop iterations of 34, 4 epilogue); two in-flight updates and the nine sender flush-path instructions give4370 total. Charge4370 per full flush with no additional reserve: the 4359 delivery/ACK body is exact static code asserted on every executed delivery, the two in-flight updates are protocol-defined, the nine sender flush-path instructions are exact static code, and both completion waits count once inside the completing store and the blocking header load by the stated transaction contract, leaving no unmodeled per-flush path. Across all senders and both phases full flushes are at most S/1024 because each consumes1024 emitted leaves. Thus the hop1 row is4370*S/1024 with the780-incidence arithmetic/emission coefficient retained.

Charge an additional21312*SIDE+12 per node for both phases' fixed block work, retaining the old cursor-initialization/finishing rows conservatively. In each phase block_init costs6151*SIDE+2, at most SIDE partial deliveries cost4361*SIDE including ACK and in-flight updates, block_finish costs at most40*SIDE+4, and ordinary final-cursor delivery costs8*SIDE. Allow64*SIDE for endpoint metadata construction and32*SIDE for phase/completion setup. Their sum is10656*SIDE+6 per phase, hence21312*SIDE+12 for SIDE=2^35. The metadata loop writes eight fixed fields per entry (payload bounds, source/channel and ACK bindings plus reserved words), using preloaded constants, incremented addresses, stores and a counted loop;64 per entry exceeds four paired zero stores, eight field stores and bounded address/control work. Constants and resident-context startup are also covered by the retained2^14 per-node setup allowance. No repeated polling is included or permitted by this construction.

All other arithmetic, cold generation, sampling, carry cap, grouping, full phase4 resolution, verification and output rows are unchanged and remain charged. The exact total is102627271243399308900963981540191989923840 instructions. Rational logarithm enclosure gives outward1262706817563496959/10^15 at C1000. The computed success lower bound remains at least0.3900000024051574611.

## Simultaneous memory and finite evidence

fixed_layout.py adds32768*SIDE+32 payload bytes,8*SIDE remote-cursor bytes,64*SIDE metadata bytes and4096 context bytes below the2^51 counter region. The pre-counter endpoint is1131627513255936<2^51. Existing bulk slots/caps and final per-node endpoint18214529204345634816<2^64 remain unchanged; total memory is below2^134 bytes and the declared2^136 bound. The absolute attempted remote cursor is at most13695129515217387520<2^64 even under the larger conservative phase leaf bound. Cursors reset only after the stated drain. Every endpoint and context uses the explicit bases, not an overlapping implicit allocation.

routing_check.py executes the native full/partial block path, ACK events, initialization, second-hop receiver, grouping and scanner over4nodes and2phases with17600 synthetic leaves, caps and independent retained-prefix references. It observes 10 full-flush publications and asserts a positive count, locking 1024-block flush-path coverage. Its event count requires exactly one native ACK per publication without a Python flag reset. These finite synchronous tests support instruction/data-path claims; they do not prove physical asynchrony or authorize a weaker completion convention. The mathematical progress argument and fixed-device contract above are material review obligations. block.s is an independent P256 implementation inspired by the public secp block/header design e59be87d-4235-4f77-ae3f-54a6ea896975; that submission's acceptance is not a premise.


# Lossless reference-assembly factoring

Only repeated literal instruction blocks in p256_inv.s and pair_inline_reference.s are factored into ordinary assembly macros. Expansion preserves every nonblank noncomment line, retaining leading indentation (instruction versus label), including labels, directives, operand spelling and instruction order. No new loop, instruction, branch or register change is introduced. The canonical expanded SHA256 values ared3b9a146237ea972367016b6b66079aa4a45c2708c28836e923f811ac53b093d and4f44679819e9d1e50306d4f31ce86604c999f0bb599871599d8f8cbdb8f80787 respectively; probe.py checks both before the retained complete tests. All expanded program delivery, initialization and execution costs remain charged as before. This reduces source text supplied to review, not computation or memory of the proposed attack. The exact changed source still requires fresh protected replay and external review.
