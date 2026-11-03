# Exact kernel body and retuned coverage budget for the powers-of-tau target

## Complete claim and source scope

For the protected BLS12-381 powers-of-tau target at n=2^20, the supplied paired class-walk Cheon construction has total work at most 2^116.285334254520161 G1-reference units (C=2000), conditional on explicitly disclosed H1-H6 of source-proof.md and s6 below. Exact score 116285334254520161/1000000000000000. Success at least 0.39 under those hypotheses. Memory below 2^122.0004 bytes, preprocessing below 2^115.2943 units and included in total work, advice zero. Only A_1 and A_(2^20) accessed, all accesses charged; same bound W(n) for 2^20<=n<r. No claim for smaller n; no auxiliary powers transferred to an ordinary track.

Algorithm, reduction, sieves, walk, planner, transport, collisions, recovery, verification and H1-H5 are in source-proof.md, the immutable source argument from 9c1b2d5d-98fb-42a4-99ab-593f6406a59f, commit d0a6e4cfea64902889b71847dbc4409709b073fc. 49e06a3d-e5c5-4fb4-8302-2106dc93104a supplies the s1 lazy-fused body + artifacts (lazy_witness.py, native_check.py, measure/side_lanes.S); its amendment prose is omitted for budget. 8d7d6d7c is cited for its H4 ladder only (remeasured by h6_ladder.py); s6 replaces its H6. Ancestor Callgrind recipes, measure/ generators, v26 log, parity-key amendment prose (a38cb931; rule restated in s1, audited by key_check.py/native_check.py) omitted (retained through k14; nothing reads them); side_k2.S comment-free except `// ----` markers (instructions identical); all audited artifacts included. Direct file base: k14 package 13efa912-0c66-4f25-a6b6-f103dcbf542d (validated), supplying sections 9-14 (GEN12/RELAY14, INS14, LOOK7+capped EXT3, CMP10/cap-2, four setup rows, chain schedule, 2^58 heads, SM77 restart, contiguous k10) plus kernel, H6, rows, bounds, repairs, budget, layout, grid (earlier k13 a17a54af, k12 96c560bb, k11 d15a39fb, k10 5cf4be66, k9 eba14ef2, k8 0779d86b, k7 e46615a7, k6 df0e2fad, k5 8dd778a3, k4 bc91c5ff, b4 780784b1; conceptual 2ed74a2c, 8d65035e, c5eeef54). This document reprices restarts exactly (section 15: SMX 647185 for SM77; restart row only, outputs identical); kernel (996), records, bookkeeping, F, N, else unchanged. Old constants/totals below are historical, not claims.

Source uses r as protected, d=2^20 dividing r-1, m=(r-1)/d, six-element GLV/negation classes, minimum pair fraction 49/50. Coverage target K=6N^2/m>=24757/50000 (0.49514). N=64239578809815195260251725009983698 is the least integer satisfying the retuned budget. Scheme and failure checks unchanged. No previous hosted decision proves H1-H6.

## 1. Exact k6 kernel body (996) with explicit fixed-call row

The side-output body is charged at 996 instructions worst-path (k6 = k5+T19; k2 file names retained), proven by a composed exact source transform (kern2_witness.py) from the measured kernel measure/side_lanes.S to side_k2.S (asserted byte equality), differentially audited vs whole-integer field arithmetic (native_check.py, 458 lanes, both branches of both fused blocks). Fixed calls funded through a separate, nonoverlapping ledger row.

The side call issues at most once per round per base, so (BASES+2)*LANE_OPS bounds invocations (retained 49e06a3d call bound). Each invocation costs at most 512 instructions incl. call preparation and return (deliberately above the 38 fixed function instructions: 36 audited + T19 prologue pair). The ledger adds an exact row of 512*(BASES+2)*LANE_OPS = 1066509928891664157045105880989696 instructions, split evenly between preprocessing and online work (asserted exact). This row is in the instruction total; fixed-call funding draws on no per-output slack, allowance capacity, or padding conservatism.

First stage (49e06a3d): p base-field prime, R=2^384>8p. For a,b<2p, word-serial Montgomery reduction returns S=(ab+qp)/R, q<R, hence S<4p^2/R+p<2p, accumulator below a+p<3p<R; so the first product (a=y1+y2<2p, b=inv<p) skips its final conditional subtraction (mu<2p); the square stays canonical. The second product takes both operands below 2p, raw S<2p. For final y, D=S-y1 in (-p,2p); one carry branch adds p once (negative) or subtracts p once (nonnegative), replacing 13+23 by at most 23. Second stage (k2,-11/lane): (T1,-2) two `mov x17,xzr` drop; `adc x4,x17,xzr` becomes `adc x4,xzr,xzr` (x17 dead; adc sets no flags); (T2,-2) n0 pinned in x18 per product (untouched, asserted); 10 reloads for spill/restore/loads via [x29,#112]; (T3,-6) branchy x' norm, `b.hs` over +2p adds (worst 28 vs 34); (T4,-1) out spill deleted, out loads [x18,#-32]. (T5,-1) deletes the second tail n0 reload (x4 unread after; `ldr` sets no flags). (T6,-1) folds `add x2,sp,#96`, mu loads direct [sp,#96..136] (x2 dead). (T7/T8/T12-T17,-8) delete both products' dead round-6 carry chains (4 links each; values feed only deleted links; flags overwritten first). (T9,-1) head `ldp x2,x3` becomes `ldp x3,x2`, y2 via x3, deleting `mov x2,x3` (x3 dies into `umulh x3,...`). (T10,-1) rides out in x0, P in x2 across the second product, deleting the y' P reload. (T11,-2) n0 stays in x18 across both products with mid-lane cursor via a frame tmp, deleting first restore, redc reload, second spill+pin. (T18,-1) that tmp rides in x4 (dead q-gap, asserted), P/T via [x4,#-24], deleting second cursor reload. (T19,-3) reassigns the two whole-lane live ranges: x30, absent from both Montgomery bodies and used once in the counted region (loop-back compare), pins n0 for the whole lane; loop-end to frame slot. Prologue (uncounted) stashes end to [x29,#112], pins n0 in x30; 18 n0 reads rename x18 to x30; spill+pin, cursor-to-x4 reload, cursor restore deleted while out/P-T revert to [x18] direct; loop-back reloads end into x4 (dead from second product's end to compare, written-first next lane, both asserted) and compares against x4, so the pin survives all lanes. Audit asserts byte-identical x', y' and parity keys on all 458 lanes (edges 0,1,p-1,p-2,R; 300 randoms; 32 genuine curve points), composed byte derivation, per-group saving in [55n+5neg,61n+5neg] after normalizing once-per-call +2 prologue pair (raw 3n-2 per group, call-row funded; T3 extras in multiples of 6), and cn<=996n+512.

A pair costs 2000+996+53=3049 instructions for two points; a root costs 2050. Primary 2000 reference charge retained untouched, incl. batch-inversion amortization. Pair bookkeeping (53) and root bookkeeping (50), incl. seven-instruction worst-path parity keys, unchanged; re-derived average record transport of s9-14 is 73/2. With pooled pair fraction F=197/200 of s6, weighted per-point charge is 627553/400 = 1568.8825 instructions.

## 2. Retuned coverage budget

Success floor is 1-exp(-K)-0.0005-2^-20-7*2^-40-2^-111, where K=6N^2/m is realized coverage and subtracted terms are H1 allowance, modelling budget, H3/H5/EXT cap-tail union, diagonal event (source-proof.md s7). The H1 allowance (0.0005) and every other term are retained exactly; only the target K changes, from 619/1250 (0.4952) to 24757/50000 (0.49514).

At K=0.49514 the exact floor is 1-exp(-0.49514)-0.0005-2^-20-7*2^-40-2^-111 = 0.39001347..., asserted above 0.39001. Declared success stays 0.39, margin 0.000013 of pure exact arithmetic on 1-exp(-K); no heuristic allowance narrowed. N is the least integer with 6N^2>=KT*m (asserted both directions), so realized K is minimal feasible, no coverage wasted.

## 3. Explicit pool-cap charge and allowance 1

The source charges an 8-instructions-per-point allowance row for unlisted bookkeeping, naming its largest known residual: the owner pool-cap compare, two instructions per baby record (source-proof.md s5.2). The base charges it explicitly: every baby record inserts once into its owner's pool; the insert path performs one pool-counter bound compare plus its not-taken branch (subs/cmp plus b.hi), for an exact row of N*2 = 128479157619630390520503450019967396 instructions. Giant records only read their head slot and walk the chain (REC_LOOK); never insert; no pool-cap check. The pool-cap row sits entirely in the baby phase, included in preprocessing in full.

With pool-cap residual and head-index combine explicit and pass control below 2^112.1 ((v) below), the allowance row drops from 2 to 1 instruction per point (2N*1 = 128479157619630390520503450019967396), removing 2N while keeping dominant residuals auditable.

Residual bounds, retained from the base. (i) Pool-cap compares: charged exactly (N x 2), residual zero. (ii) Key reloads: zero (key stored to slot, reloaded by listed ldp+ldr, travels in-record). (iii) Head-index combine: charged exactly (2N; s11 4+54 split; lsr+bfi; consumed next instruction). (iv) Cap enforcement: charged exactly (2 x cmp_cap x OWNERS, online only; counter in register, wave-end traffic in link row). (v) Pass control below 2^112.1: owner-init unrolled 2^20-fold (control below 2^99); all other bulk rows below 2^109 in sum with at most subs+b.ne control per level over four levels (below 8x sum, below 2^112); vs retained 2N (2^116.63). The allowance covers only purely-unknown misses.

## 4. Exact complete ledger

Unchanged source dimensions, quotas, capped schedules evaluated by evidence.py at (KT, F, KERNEL). Every row below is an exact integer instruction count, rounding at original locations. The complete instruction total is 202490038095811793825301726699459094871, divided by 2000 for W. Preprocessing 101862524537629152940278063598999940138 instructions, included (falls with restarts row). Baby-phase transport charged at 40/record vs combined 36.5 average; baby-phase walk part at conservative floor F_PRE=39/40 of s6, not pooled F.

- base interval sieve, characters, assignment: 14778874133111650197864576.
- candidate delivery and replies: 5488762045857792.
- candidate formation: 281474976710656.
- candidate re-walks: 1257895200023828017403475787776.
- candidate tests: 13858116485374300651520.
- candidate-link wave counts: 48961495694392481575600128.
- chain-id formation: 30635557439458550347906155000.
- character sieve, class lists, their transfer: 30225869176339564179694288896.
- comparison cap enforcement: 25085181676799472938181262097514496.
- comparisons at the per-owner cap: 125425908383997364690906310487572480.
- gap tables: 127161744382312499152959732592009770.
- global barriers: 6630756600764919022560032736.
- head-index word combine: 128479157619630390520503450019967396.
- listing allowance (1/point): 128479157619630390520503450019967396.
- nonempty-lookup extracts and cap: 34564008753675207535941252726915072.
- owner initialisation: 166153499473114557899952177373249536.
- owner pool-cap compares: 128479157619630390520503450019967396.
- padding ops at the pair price: 50807649054705996927880279846027264.
- per-wave link words, memory barriers: 48668920232217315802772306657280.
- plan distribution (copy tree): 121734115224496661965030475956224.
- planner (scans, matching, emission): 487728773955200596712816258318336.
- planner descendant traversals (adjacency lists): 10161011842454401296372011630592.
- records: 2N points: 201568702004179776155783753925951498155.
- root recovery and input distribution: 73872550223301050368.
- segment restarts: 4956840852593753310928802855538750.
- side-output fixed calls: 1066509928891664157045105880989696.
- table strand starts and base points: 7657983763737776997925928960000.
- window-table initialisation: 10880332376531662572355584.

Memory and all per-node layout assertions are retained under the s8 slot size (generating-node term updated, owner term doubled by s11): total below 2^122.0004 bytes and every node below 2^63 bytes. No free precomputation, future-target amortization, or wall-clock parallel discount. Exponentially large node/memory assumptions are analytical accounting conventions, not feasible hardware.

Score grid is 10^-15: evidence.py reports the least fifteen-decimal ceiling 116285334254520161/1000000000000000 of exact log2 W = 116.2853342545201603... (90-digit Decimal enclosure; on-grid assertion holds). log_bound.py encloses ln(x) by 64 rational terms of 2*sum(t^(2j+1)/(2j+1)), t=(x-1)/(x+1), plus the positive geometric tail; similarly ln2; divides in the safe direction. Applied to the exact total/2000, it proves the declared fifteen-decimal ceiling.

## 5. Success and honest limitations

Actual points, canonical stored keys, routing, candidate decisions, and recovered scalars match the preceding algorithm. Neither allowance-1 accounting and kernel micro-steps, the k7 layout amendment (s8), the s9 re-derivation (same caps/routing/counts/events; sequences change), the s10 memory schedule (same events; counts change), the s11 head doubling (comparison cap halved, H3 re-verified), the s12 conditional extracts (EXT cap failure in tail budget), the s13/s15 restart repricings (same outputs; price only), nor the s14 layout (same events; EXT at true offset) change any source success event; x', y' and parity keys are byte-identical on all audited paths. N is unchanged. H1-H5 keep definitions, evidence, failure budgets from source-proof.md s7 + planner sections; H6 is in s6. The allowance change is exact accounting with a stated residual enumeration. Those hypotheses remain extrapolations requiring review.

The recomputed floor is 1-exp(-0.49514)-0.0005-2^-20-7*2^-40-2^-111 > 0.39001347; declared success is only 0.39. Conditional on the disclosed model, not an empirical production rate or reviewer-confidence statement. All failures and capped work remain in the complete ledger.

The full combined test reruns the production-parameter integer ledger, planner controls, reduced-order curve recovery, synchronous transport with overflow controls, toy census, 12142 parity-key cases, plus the 458-lane native audit. Tests exercise finite cases only; algebra and accounting above carry the production extrapolation. No full-size attack run.

## 6. Pooled pair fraction 197/200: H4 class mean plus jointly measured H6 weighting

The scored records row prices all 2N emitted points at pooled pair fraction F=197/200. Two premises: production class-mean pair fraction at least 0.9862 (H4 mean); production selected-base-weighted pooled fraction at least 0.985 (new disclosed H6). The per-class floor 49/50 of source H4 is retained alongside both, bounding the worst case.

H4 mean. The shipped planner ran at production degree GMAX=16d^2 on real chi-class lists (h6_ladder.py mean-only specs, rerunnable; hosted run recomputes two small joint specs (h6_hosted); table quoted from the same code offline, as with source planner tables):

| d | span X2/2GMAX | classes | mean f | min f |
| --- | --- | --- | --- | --- |
| 16 | 16.0 | 4 | 0.98651 | 0.98601 |
| 32 | 8.0 | 4 | 0.98689 | 0.98589 |
| 64 | 4.0 | 2 | 0.98646 | 0.98502 |
| 128 | 4.0 | 2 | 0.98776 | 0.98751 |

Every row sits above 0.9862; source tables agree (0.9861 at d=16 span 16, 0.9874 at d=64 span 16), means rising in span and d (uncapped controls 0.9794, 0.9834, 0.9858, 0.9859, 0.9879 at d=8..128). Production (span 2^31, d=2^20) lies far beyond every measured row on both trends, and every deficit source (boundary GMAX/span, cycle 1/d, demotions) vanishes at production ratios. H4-mean is the disclosed extrapolation that level and trends continue; the production run is infeasible and unclaimed.

H6 (joint weighting). Production emitted-point weights depend jointly on class lengths, opposite-character base-prime counts, and final partial walks. h6_ladder.py measures this joint process at toy scale: real chi-class lists, shipped planner per class, bases as consecutive primes above X2, baby taking the first N_t elements and giant the next (last base partial), each base walking its opposite-character class plan. Twenty-two joint specs (d=8..128, 7..25 walked classes, 20..54 bases, uncapped + ratio-capped max_tree=d controls) split pooled-minus-class-mean into selection, length, base-count and partial-walk steps. `python3 h6_ladder.py < offline-specs.json` reproduces the 22-row table (about 100 seconds, stdlib); summary:

- pooled-minus-walked-mean, all 22 rows: range [-0.00185, +0.00282], mean +0.00015;
- pooled-minus-class-mean, all 22 rows: range [-0.00214, +0.00233], mean +0.00040;
- in the 15 rows with at least 12 walked classes, |pooled-minus-walked-mean| <= 0.00050;
- 3 rows above 0.0005 all walk at most 8 classes with up to 31% of emitted weight on one class;
- length-weighting step, all 22 rows: range [-0.00001, +0.00057], nonnegative in 20 of 22 (other two -0.00001, rounding-scale);
- base-count step: [-0.00133, +0.00191]; partial-walk step: [-0.00089, +0.00091];
- capped controls (cap binding, fractions 0.82-0.94): deltas -0.00043, +0.00023, +0.00282, same scale.

Production weights are 10^4x+ more diffuse than the most concentrated measured row: one class holds at most 2^42 bases x 2^57 elements = 2^99 of N ~= 2^115.63 emitted points per side (max single-class weight 2^-16.6 vs 1/7 worst measured); all 2^20 classes are walked (~2^39 bases per character per side), so selection vanishes. Measured deltas center on zero; spread shrinks in diffuseness; no systematic negative component across 22 diverse fixed weightings. H6 margin 0.9862-0.985=0.0012 covers 2.4x the worst drag measured at production-like diffuseness (0.00050). Partial-walk drag is analytic: at most one remainder per side, under 2^57 elements each, deficit under 0.3 vs N ~= 2^115.63 points, i.e. under 2^-60 pooled. Adversarial weight alignment with low-fraction classes would break H6, but is measured-absent and mechanism-opposed (H2 equidistribution over 2^20 diffuse classes); residual risk disclosed, not hidden.

Phase-safe accounting. The combined records row is exactly linear in F (2N x (F x PAIR/2 + (1-F) x ROOT + 73/2)), so the pooled H6 bound covers any baby/giant split with no per-phase fraction asserted for the scored total. Preprocessing funded separately and conservatively: baby-phase walk part uses F_PRE=39/40=0.975, below retained H4 per-class floor, no weighting claim (101862524537629152940278063598999940138 instructions, below 2^115.2943).

## 7. Segment-index encoding repair (33 bits)

The source specifies a 32-bit segment index in the record chain id (source-proof.md s4.2) while allowing classes with up to 2^33 segments (s3 under H2), so the index overflows on long classes. This package amends the encoding explicitly (verdict: procedural context only): segment-index field is 33 bits in the pool record (62-bit base + 33-bit segment + 56-bit link = 226 of 256 bits, within 31 spare bits) and the 128-bit transport chain id (62-bit base + 33-bit segment + spare = 95 of 128 bits). The segment directory already holds 2^33+1 eight-byte entries, and the 62-bit base index is unchanged (at most 2^61.6 bases). Zero instruction cost: same-word stores/loads; the re-walk index extract widens an immediate only (ubfx width 32 to 33, one instruction); record sizes, memory, occupancy and cap analysis unchanged.

## 8. Production lane-stack layout (128-byte slots)

This section supersedes the 96-byte lane-stack stride of source-proof.md s5.2 (k6 verdict's layout gap, verdict: procedural context only); other source listings stand. Every lane-stack slot is 128 bytes: [0,96) point coords (x [0,48), y [48,96)), [96,120) 24-byte parity key (limbs 1..3), [120,128) reserved pad no charged instruction touches. Slots are contiguous at stride 128 from a 16-aligned base, so every slot is 16-aligned and every listed ldp/stp/ldr/str keeps form and offsets. Kernel writes exactly [out,out+120): x' at 0/16/32, y' at 48/64/80, key at 96/112; out always at a slot base, next slot at +128, so no lane touches another slot. Gap tables stay 96 bytes per entry (read-only (x,y); kernel never writes T; strands store points only, never keys); 96-byte batch lane buffers stay (denominator/inverse halves, no keys).

Zero instruction delta: midpoint-load and P+T-store `madd` keep single-instruction hoisted form with constant 128 for 96; T[g] load keeps its 12; batch buffers, argument moves, P+T key 3-word store (in slot at +96), key-reload ldp+ldr (+96/+112), root load/store instruction-identical. Memory: lane-stack term LANES*(2*DEPTH+1)*128; total below 2^122.0004, every node below 2^63 (asserted). native_check.py adds a production-stride audit: 8 lanes, P and out in one contiguous 128-stride buffer, outputs vs the whole-integer reference, pad/trailing canaries intact, P inputs unchanged, write footprint exactly out ranges plus frame; a 96-stride negative control corrupts its trailing canary, proving the audit catches the k6 bug.

## 9. Record-path re-derivation: pointer buffers, honest matching, explicit setup

This section re-derives every per-record instruction from scratch, superseding the source s5.2 transport counts. The five charged sequences (GEN12/RELAY14/INS16/LOOK9/CMP12, GEN12 per s10) are asserted by evidence.py audit_record_path (lengths, bits, immediates, pool fit, division proofs); the proto toy executes the pointer mechanism on every record buffer (pre-write check, +40 advance, exact derivation, reset). Score delta vs k7: -0.001670367816893.

Pointer buffers (GEN 13 to 11, RELAY 16 to 14). Sender per-link state becomes (ptr, end): each record's three stores post-index ptr by 16+16+8 = 40 bytes, deleting slot-address madd and counter increment. Entries stay 16 bytes (one ldp; source "24-byte" prose counted a hoisted cap word), addressed by one add-with-shift (t/u below 2^30, lsl #4). The cap check precedes the write (cmp ptr, end; b.hs fail): identical failure event, never an out-of-bounds write. GEN11 is ubfx t, entry address, ldp, cmp, b.hs, three stores, ptr store-back (interleaved links, pointers spill per record), plus key-reload/loop pair. RELAY14 is the same output side; input loads and loop unchanged. At wave end the sender derives count = (ptr-base)/40 (sub, lsr #3, umull/lsr /5 with magic 0xCCCCCCCD, exact by m*5 = 2^34+1 plus asserted bound, count store, ptr reset: 8), writes the per-parity count word, and resets; receiver side (3) and double-buffer discipline unchanged. LINK_WAVE is 16 (8 sender + 3 receiver + 5 margin); the per-node term rises to NODES*24 for dmb, counter traffic and constants (16 used, 18 with s12 EXT traffic). Chain-id words are segment-invariant, formed once per segment (ubfx+bfi+lsr: 3) by the formation row and register-held across the segment's records, so GEN's five stored words are all in registers with no per-record formation cost.

Honest matching (INS 16 verified, LOOK 7 to 9, comparison 7 to 12, cap 1 to 2). Record key: owner 60 bits (implicit), head 57 bits 60..116 (straddling: listed lsr + combine row), compare key 75 bits (k11 bits 117..127 + k64 = w2). The pool record stores k64, link56 with k11low8, w3 verbatim, w4 verbatim with k11high3 in chain spare (226 of 256 bits): pack is two ubfx/bfi pairs (4, as listed), the counter add is the listed 1 (counter in register; wave-end traffic in link row), all other INS items stand: INS 16 verified. LOOK9 loads only the three key words (ldp [x],#16; ldr [x],#24 advances exactly 40, skipping chain words), head lsr, head load, a cbz skipping empty chains (pool slot 0 unused, counter from 1, so zero means empty/end; no sentinel compare; cap counts chain records only), and two k11 ubfx. Each comparison is 12: link scale (lsl #5; links below 2^56 by pool-cap row), three loads, two k11 ubfx, orr, cmp, ccmp, equal branch, link ubfx (no high-bits invariant), cbnz merging null check with loop branch. The cap row is subs plus exhaustion branch (2). Four new explicit rows: window-table init (8/link, 4 used), chain-id formation (3/segment, above), candidate-link wave counts (9 on each of 4 traversals: record /40 magic, reply /64 shifts, candidate /56 with proven /7 magic), candidate formation (16/match, 11 used). Caps, routing, counts, failure events unchanged, so H1-H6, N, success, memory untouched.

## 10. k9 repair: chain-word memory schedule

This section repairs the k8 verdict's chain-register gap (verdict: procedural context only) and supersedes two section-9 sentences (GEN11 listing; register-held chain words). Two chain-id words formed once per segment into a per-lane 16-byte slot (LANES slots per generating node, in memory term); every GEN reloads them. Formation 4 per segment (ubfx+bfi+lsr plus stp); GEN11 becomes GEN12 (ubfx t, entry address, ldp entry, ldp chain words, cmp, b.hs, three stores, ptr store-back, key-reload/loop pair). No chain word crosses records in a register, so lane interleaves and side-output clobbering are compatible; allowance untouched. Delta vs k8: +0.000913879237243; vs k7: -0.000756488579650. The manifest states GEN12, LOOK9 and four setup rows exactly as charged.

## 11. k10 heads-x2: 2^58 slots per owner

This section doubles the per-owner head array, superseding the s9 57-bit head width and the head-slot/H3-comparison sentences of source-proof.md sections 4.3/7. Head is key bits 60..117 (58 bits, straddling w0/w1 as 4+54; lsr + one combine-row bfi: width 53 to 54, same two instructions); compare key 74 bits (k10 bits 118..127 + k64). LOOK9/CMP12 extracts become k10 ubfx (widths only); pool record holds link56+k10low8 and k10high2 in chain spare (225 of 256 bits). All listing counts stand. cmp_cap halves to 10878963388471823 per owner; the H3 comparison tail above 1.01x mean (Bernstein) stays below 2^-40 over all nodes at the halved mean (exponent linear in mu, margin enormous); toy support model-level, unchanged. Owner init doubles (2^57-stp stream, control below 2^99); (v) stays below 2^112.1. Delta vs k9: -0.000658038009622.

## 12. k11 conditional extracts: LOOK7 base plus occupancy-capped EXT4

This section moves LOOK9's two k10 extracts behind the nonempty test, charging them at an occupancy cap. Every lookup executes LOOK7 (ldp, ldr, lsr, ldr, cbz, subs, b.ne); only fall-through lookups execute EXT4 (subs + b.lo on a per-owner down-counter, then two k10 ubfx). The counter lives in a register across the owner's lookups with wave-end traffic in the link row (s9 pattern; link-row use 16 to 18 of 24); exhaustion aborts into the tail budget: cap-plus-tail, never expectation. Per-owner cap is ceil(1.02 x mu_L x p) = 9993166234811374 with mu_L = N/2^60 exact, p = y-y^2/2+y^3/6 at y = N/(OWNERSxHEADS-1): under H3 each key is uniform over the OC cells, so P(nonempty) = 1-(1-1/OC)^N <= 1-e^-y <= p (ln(1-x) >= -x/(1-x); decreasing alternating Taylor, asserted), and E[nonempty] <= mu_L x p by the tower argument over independent baby/giant placements. Three 1.005 stages cover the cap (1.005^3 = 1.0150751 < 1.02, asserted): lookups/owner (Binomial, Bernstein), occupied/owner (negatively-associated, Chernoff, Dubhashi-Panconesi), probes-given-both (Binomial, stochastic dominance). Per-owner means exceed 1e15 (asserted); at delta >= 0.005 every Chernoff exponent exceeds 1e10 nats vs the 70.5-nat budget (2^-40 over 2^60 owners), so the success union grows to 7x2^-40. The toy counts nonempty probes per lookup phase. Base drops 2N, EXT row adds 4 x OWNERS x cap; capped EXT beats 2/lookup (asserted). Delta vs k10: -0.000586272451324.

## 13. k12 77-bit restarts

This section repriced restarts at 77-bit width (s15 supersedes the price; the [b]Q_a identification stands). Every restart computes a root from the 77-bit b in its plan entry (source-proof.md sections 2-4: b (77 bits) per entry; segments begin at restarts; each lane handoff executes one), so all ceil(2N/SEG)+BASESxLANES counted restarts are [b]Q_a with b < 2^77. Was SM77 = ceil(77xSM/255) = 2473663: same double-and-add program with loop bound 255 to 77 (77 doublings, <=77 additions), same per-bit rate incl. control, setup covered 1.2x by scaled 4x slack (4x77/255 = 1.208 > 1, asserted). +16 handoff unchanged. Strand starts, base points, re-walks, candidate tests stay at full SM (conservative). Preprocessing falls with the row's halved share. Delta vs k11: -0.000311727843010.

## 14. k14 contiguous k10: INS14/EXT3/CMP10 (giant EXT at true offset)

This section supersedes the split k10 layout of sections 9/11/12 (link56+k10low8, k10high2 in chain spare; two pack pairs; two extracts joined by orr). The k10 field (key bits 118..127, second-word bits 54..63) is stored contiguously: one ubfx+bfi pair deposits all 10 bits into chain-spare w4[31..40] (formation lsr zeroes w4-high; (31,.) audit immediate), and the link word carries link56 low with zero high, wasting 8 link bits. Pool use stands at 225 of 256 bits (64+56+10+62+33, asserted); occupancy, caps, H3, counts, memory, success untouched: layout only, no fresh crypto. Single-writer per word (pack before stores), so no tearing argument.

INS is 14 (ubfx (54,10) from the baby second word, bfi (31,10) into the chain word), EXT 3 (subs, b.lo, one ubfx (54,10) on the giant second-word register), CMP 10 (k10 ubfx (31,10) from the baby chain word, orr deleted; link ubfx (0,56) unchanged). Register schedule: LOOK7's ldp/ldr leave the three giant key words in registers; lsr writes a scratch, head ldr writes the slot register, cbz/subs/b.ne write no key-word register, so the second word is intact at EXT; extracted k10 + k64 word stay in registers across the chain walk (loads/compares only), feeding cmp/ccmp. Immediates ((0,30),(30,30),(54,10),(31,10),(0,56),(4,54)) with per-side mapping (INS (54,10)/(31,10), EXT (54,10), CMP (31,10)); baby-phase 40, average 73/2, per-point 627553/400. Records -2N, EXT x3/4, comparisons x10/12 (cap unchanged). Delta vs validated k12: -0.001175637155688.

Semantics audited, not just counted: 2^54 yields k10=1 on both sides; 8 adversarial + 1500 random second words agree through pack/extract vs the layout oracle; 1500 identical full keys match, k10-bit flips mismatch, other bits inert; the toy routes every comparison through packed ops, asserting packed iff full-key decision. Identical keys match, so H1 coverage implies recovery as in the validated base.

## 15. k15 exact restart price SMX

This section supersedes s13 SM77; only the restarts row changes. Each restart computes [b]Q_a, b in [1,X2] (bit-length <=77), by branchy Jacobian double-and-add, all-ones worst path charged uniformly. Qa enters canonical affine; outputs canonical affine proved equal to [b]Q_a, so success is unchanged; strands/tests stay full SM.

M=325 (lazy-first 317 + glue 8); Mc=337 (+12 canon); S=257 (square 92 + reload 3 + mul_by_1 150 + high-add 6 + glue 9); Sc=269; C1=158. M/S values below 3p/2 (R>8p asserted). Smalls (p in x5-x10, n0 in x30): 29/41/29/26/38; x3=55, x4=52, x8=78. DBL=2404, ADD=3814 (14 sub-feeding uses take _c form); ladder 76x6218=472568; loop 380, scan 233. R=[k]Q_a, 1<=k<2^77, never O (k<r; ADDs at m>=2 with 2m=+-1 needing m~2^254). Fermat (p-2 fixed, 380S+228M unrolled) directly yields Montgomery Z^-1: 171760. Conversion 1572, entry 656, setup 16. SMX=647185. Preprocessing falls; memory unchanged. Delta vs k14: -0.000099663495972. Machine-checked: emulated use-sites/smalls, ladder vs affine over Fp, recounts.
