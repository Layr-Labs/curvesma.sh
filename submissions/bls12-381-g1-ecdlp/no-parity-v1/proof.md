# BLS powers: balanced giant/baby ratio (115.887900713539257)

## Complete claim and source scope

For the protected BLS12-381 powers-of-tau target at n=2^20, the supplied paired class-walk Cheon construction has total work at most 2^115.887900713539257 G1-reference units (C=2000), conditional on explicitly disclosed H1-H5 of source-proof.md and H7 (s37), with H4 restated as the per-class floors of s30 (pair fraction F = 998/1000 and leaf share 193/200 since s36, both runtime-checked; s6 records the withdrawn pooled argument). Exact score 115887900713539257/1000000000000000. Success at least 0.39 under those hypotheses. Memory below 2^122.0010 bytes, preprocessing below 2^114.9084 units and included in total work, advice zero. Only A_1 and A_(2^20) accessed, all accesses charged; same bound W(n) for 2^20<=n<r. No claim for smaller n; no auxiliary powers transferred to an ordinary track.

The direct complete-file source is our published submission
3c7fc300-a7a6-466a-b221-ccc30ec9685a, commit
d0302d68c68a44d048ef3d2411673393681dbc76. It supplied sections1-39 and
the accompanying code, with section 40 removing unused parity keys.
Section 41 specifies the only production change here: analytically balance
the giant/baby ratio Q = 115917/100000.

Source uses r as protected, d=2^20 dividing r-1, m=(r-1)/d, six-element GLV/negation classes, Coverage target K=6 N_b N_g/m>=12377947373/(25x10^9) (0.49511789492, s2) with N_b = ceil(115917 N_g/100000) (s41): N_g=59664969621627809652901534016047992 is the least integer satisfying it, N_b=69161842836302308115353871185382351. No previous hosted decision proves H1-H6.

## 1. Inherited arithmetic derivation; parity output superseded by s40

The historical raw side_lanes(args,n) in measure/side_lanes.S kept the k16 interface: args[i] = {out, P, T, inv}; it writes fully reduced x' of P-T to out+0, y' to out+48, the parity key (limbs 1..3 of whichever of y' and p-y' is even, a38cb931) to out+96; args[n] = {NULL, NULL} ends the lane loop. Inputs: x1, y1, x2, y2 canonical Montgomery residues (P a kernel output, restart root or strand start; T a gap-table entry), y2 != 0 (G1 has odd prime order), inv = 1/(x2-x1) mod p from the batch round, below 1.17p but not necessarily canonical (s27). Each gap-table entry stores v2 = p - x2 at T+96, n2 = p - y2 at T+144 (192-byte entries; +21/entry in the gap-table step, TABLE_BK 48 <= 50; re-walk rebuild 2000+50; memory GMAX x 192; hoisted stride holds 192). Base-point/strand-start entries get v/n inside their scalar-mult charge (SM = 8192000/point, ~4x a double-and-add).

Arithmetic (unchanged; bounds at the s27 inverse range). R = 2^384, p/R < 0.1016; the word-serial Montgomery product returns S = (ab + qp)/R < ab/R + p before its final conditional subtraction, which the kernel omits. blst f7a88067's __mul_mont_384 (twice), __sqr_384, __mul_by_1_mont_384 inlined with n0 in x18, modulus in x5..x10, the 768-bit square in registers. Steps: t = y1 + y2 < 2p; mu = S(t, inv) < 1.24p; V = S(mu^2) < 1.16p; X = V + v2 - x1 in (-p, 2.16p) reduced to x' by a sign branch or one/two borrow-checked subtractions; w = x' + v2 in (0, 2p); z = S(w, mu) < 1.26p; y' = mu(x' - x2) + y2 = mu(x' - x1) - y1 (line through P and -T; mu(x2 - x1) = y1 + y2), i.e. D = z - n2 in (-p, 1.26p), reduced by one branch and stored with the key block (tbz, <= 4 borrows, stp, str: 7 worst path; unread since s30, which keys records by x). Code after the x' reduction is emitted twice, so no x' path jumps through shared code.

Derivation (measure/derive_k929.py of k31 78c29cd4 and earlier rungs; not shipped here, nor the 971 file it reads). The shipped file is produced byte for byte from the measured 971 listing of 99e03e91/fd385c98 (SHA-256 ba448d6b7d1623febeae2c4b9b9e7f343232e1651639d1c3cc5cde9e709f355a, the file of its Callgrind record; asserted by the script) by two text transforms, every edit asserted by count:

K1, value-bound carry deletion (-37). Lemma: in a word-serial Montgomery product with a < 2p and any 384-bit b, acc_j < a + p before round j, and every pre-shift sum acc_j + a b_j + q_j p is below 2^64 (a + p). Induction: acc_0 = 0; if acc_j <= a + p - 1 then acc_j + a b_j + q_j p <= (a + p - 1) + (2^64 - 1)(a + p) = 2^64 (a + p) - 1, so acc_{j+1} < a + p. With a + p < 3p < 2^383, every pre-shift sum is below 2^447, every shifted accumulator below 2^383; every in-round partial sum is a sub-sum of nonnegative pre-shift terms. So word 7 of every partial sum and word 6 of every shifted accumulator are zero. The 971 listing used this for the final round only; K1 applies it to rounds 0..4 of all three product copies, deleting the carry words adc x17,xzr,xzr; adc x17,x17,xzr; adc x4,{x17|xzr},xzr; adc x25,x4,xzr (18 per product), and turning the next round's first read of zero word 6, adcs x25,x25,xzr, into adc x25,xzr,xzr (same count; its carry-out fed only a deleted word; the next instruction resets flags). Both products satisfy a < 2p: t = y1 + y2 <= 2p - 2, w = x' + v2 <= 2p - 1. In the square, mu < 1.24p < 2^382, so its top word is below 2^62, the cross sum sum_{i<j} mu_i mu_j 2^(64(i+j)) is below 2^702 + 2^639 < 2^703 and doubling words 1..10 never carries into word 11: adc x1,xzr,xzr deleted, adc x16,x16,x1 reads xzr (-1).

K2, paired b-word loads (-5). After K1, x17 holds only the current b word in each product, so b words load in ldp pairs with a partner register dead inside that product (x29 in the first: square scratch, restored by the epilogue; x2 = out in the second, reloaded once after it); the transform asserts the partner appears nowhere else in its product.

Instruction counts per lane (static, from the shipped file; 971 listing in brackets): first product 295 [316], square 91 [92], square reduction 150 [150], second product 296 [316], driver worst path 97 [97]; total 929 [971]. Static enumeration gives 24 lane paths, worst 929, best 918; every K1/K2 edit lies on all 24. The prologue (22), epilogue (10) and terminator total 32 per call, constant, plus the last lane's exit branch; the fixed-call row charges the itemized 40 per invocation of section 20 against the measured 33 (kernel_check.py asserts the fixed part equals 32 on every call; 32 + 1 + 5 + 2 = 40). Since s30 each round issues eight calls (four side-output kernels, four beta_keys passes), at most once per round per base, so 8 x (BASES+2) x LANE_OPS bounds invocations and the row is (BASES+2) x (8 x 48 x LANE_OPS + 64) instructions (48 since s37) (value: s4), both phases charged, inside the instruction total and drawing on no allowance or padding capacity.

Inherited measurement: measure/measurement.log identifies the shipped
assembly (SHA-256 b4c7a7c2d17e6decfc8c0ee93f1def949602ad9d7e2ee0edff4a7ec3c8e4be72),
GCC 14.2.0, Valgrind 3.24.0 and blst f7a88067. It reports
(Ir(8192)-Ir(4096))/4096=922.925 instructions/lane, below 929,
and field multiplication 383. Its harness self-checks passed. This package
inherits that record without claiming a new native AArch64 measurement.

Finite audits are executable in kernel_check.py and nocarry_witness.py.
They check the shipped assembly against exact field and curve arithmetic,
all 24 enumerated lane paths, input/callee-saved preservation, the current
output stride and write footprint, and a failing old-stride negative control.
Cases include extreme operands, genuine curve pairs and unreduced inverses.
The carry witness checks the K1 recurrence and square bounds on extreme and
seeded random inputs. These tests support the analytic lemmas above; their
actual current results appear in the declared experiment output.

Pair and root prices are those of s26-s28 per child with the leaf kernel for leaf children, s30-s31 side prices and s32 kernels: baby point 3615580199499/3276800000 and giant point 4192706599499/3276800000 instructions (pinned; s37, s38), no reference addition in the walk. Root 243368215/131072 = 1856.7521...

## 2. Coverage budget at the exact floor limit

Success floor is 1-exp(-K)-0.0005-2^-20-7*2^-40-2^-111, where K=6 N_b N_g/m is realized coverage and the subtracted terms are the H1 allowance, modelling budget, H3/H5/EXT cap-tail union and diagonal event (source-proof.md s7, s11-s12 below). H1 (source-proof.md s7) is the hypothesis that |R|/m at the production construction is at least 1-e^-K-0.0005 at the operating K, with its census support and ln(X2)/d extrapolation disclosed there; it is applied directly at the realized K of this package. No monotonicity of the coverage deficit in N or K is asserted or needed.

The floor is evaluated exactly: evidence.py encloses exp(-K) by the alternating Taylor series (40 terms, both partial-sum brackets, width below 10^-50) in rational arithmetic and takes the lower end of the floor; the target KT = 12377947373/(25 x 10^9) = 0.49511789492 is the least multiple of 10^-12 whose floor lower bound is at least 0.39 + 10^-10 (asserted in both directions), the same tightening pepedesigner's 54ffdd65 made on a 10^-9 grid with a 10^-7 margin. At the realized K the exact lower bound is 0.39000000010024666..., asserted at or above 0.39 + 10^-10; the declared success is 0.39. N_g is the least integer with 6 ceil(29 N_g/25) N_g >= KT*m (asserted both directions; s39), 59643620144729401137567435202637066, N_b = 69186599367886105319578224835058997, so realized K >= KT and no coverage is wasted. Every N-dependent quantity is regenerated from N_b, N_g by the same code: records, pool-cap wave-check row (s18), allowance, restarts, bases, gap tables, strands, the s11 comparison cap (32312633143699954 per owner, 3 N_g lookups) and the s12 occupancy cap (29182932951402857 per owner) with their asserted exponents and ratios (below). The tail arguments of s11 and s12 are evaluated at these values, not inherited.

## 3. Explicit pool-cap charge and allowance 1

The one-instruction-per-point allowance covers residual bookkeeping; named pool-cap checks and key/index access costs are charged explicitly.

Residual bounds, retained from the base. (i) Pool-cap checks: charged exactly (wave row of s18), residual zero. (ii) Key reloads: zero (key stored to slot, reloaded by listed ldp+ldr, travels in-record). (iii) Head index: one lsr (s17), combine row deleted, residual zero. (iv) Cap enforcement: charged exactly (2 x cmp_cap x OWNERS, online only; counter in register, wave-end traffic in link row). (v) Pass control below 2^114.1: owner-init unrolled 2^20-fold (control below 2^99); all other bulk rows below 2^111 in sum (since s30 the largest are gap-table entries 2^108.65, planner scan visits 2^108.58 and augmentation ancestor steps 2^108.4, added on top of the former 2^109 budget) with at most subs+b.ne control per level over four levels (below 8x sum, below 2^114); vs retained N_b+N_g (2^116.65). The allowance covers only purely-unknown misses.

## 4. Exact complete ledger

Unchanged source dimensions, quotas, capped schedules evaluated by evidence.py at (KT, F, ELL, KERNEL=916 (s40), KX, BETA = 128, BATCH, the s28 share, N_b, N_g; s38 tables). Every row below is an exact integer instruction count, rounding at original locations. The complete instruction total is 153731998392942372608847136489421039335, divided by 2000 for W. Preprocessing 77986473067904468812759746943788259437 instructions, included. Baby transport 2180/64 per record, giant 3 x 1796/64 per point (s30); baby-phase walk part at conservative floor F_PRE=39/40 of s6, not F. Preprocessing takes the baby share QB = max(N_b/(N_b+N_g), BASES_B/BASES, (BASES_B+1)/(BASES+2)) of the per-point and per-base rows (padding, side calls, allowance, restarts, gap tables, strands, chain ids) and half of the wave rows (baby waves are under half of all waves, asserted).

- records: N_b + N_g points: 152654387983403776134121040908485828775.
- gap tables: 378674037237058483341550182416817480.
- comparisons at the per-owner cap: 260777507352905359601837545753673728.
- listing allowance (1/point): 128830219512615506457145660037696063.
- nonempty-lookup extracts and cap: 100936892901516346244119520184631296.
- owner initialisation: 83076749736557315843464236105728000.
- comparison cap enforcement: 74507859243687245600525013072478208.
- the other 28 rows, together 50807143554256217637453423364185785 (0.033% of W, each below 0.05%): names and values in the experiment output (ledger rows_instructions).

Memory and all per-node layout assertions are retained under the s8 slot size (generating-node term with 192-byte gap entries, owner term doubled by s11): total below 2^122.0009 bytes and each non-generating node below 2^63 bytes, generating node below 2^64 (s38). No free precomputation, future-target amortization, or wall-clock parallel discount. Exponentially large node/memory assumptions are analytical accounting conventions, not feasible hardware.

Score grid is 10^-15: evidence.py reports the least fifteen-decimal ceiling 115887900805657670/1000000000000000 of exact log2 W = 115.8879008056576693173898... (90-digit Decimal enclosure; on-grid assertion holds). log_bound.py encloses ln(x) by 64 rational terms of 2*sum(t^(2j+1)/(2j+1)), t=(x-1)/(x+1), plus the positive geometric tail; similarly ln2; divides in the safe direction. Applied to the exact total/2000, it proves the declared fifteen-decimal ceiling.

## 5. Success and honest limitations

Points, stored keys, routing, candidate decisions, recovered scalars match the algorithm. No section s1-s37 changes a source success event: each changes prices, layouts, bit-identical listings or superset failure events in the same union, as each section states; current full kernels preserve exact x',y'; beta passes produce the live routed keys. Parity output is removed by s40. s30's x-keys with three giant lookups test exactly the U-orbit match the parity key tested (orbit lemma in s30); leaf kernels and the s30 planner change prices, not points. N changes only through s2. H1-H5 keep source-proof.md s7 definitions, evidence, budgets; H7 (s37) is a new disclosed hypothesis; H6 is superseded by the s30 per-class floors (H4 as f_c >= 998/1000 and H4L as leaf share >= 193/200 (s36) per production class, runtime-checked before any group op; a failed check stops the run, challenge-independent). The allowance change is exact accounting with stated residuals. Those hypotheses remain extrapolations requiring review.

The exact lower bound of the floor at the realized K (s2) is 0.39000000010024666..., at or above 0.39 + 10^-10 (asserted in rational arithmetic); declared success is 0.39. Conditional on the disclosed model, not an empirical production rate or reviewer-confidence statement. All failures and capped work remain in the complete ledger.

Mode all of combined_lazy.py reruns the production integer ledger, planner controls, reduced-order curve recovery, transport with overflow controls, the toy census at 29/25, the executed record-path audits, the s40 full kernel, unchanged leaf and batch-round audits (s1, s30, s32), the s31 beta keys and the K1 witness. The toy curve search skips a prime only when none of its six j = 0 orders can pass the acceptance test (cm_ok; speed only). Finite cases only; algebra and accounting carry the production extrapolation. No full-size attack run.

## 6. Superseded pooled floor

H6 is withdrawn. Use the runtime-checked per-class floors of s30;
preprocessing uses F_PRE=39/40.

## 7. Segment-index encoding repair (33 bits)

The source specifies a 32-bit segment index in the record chain id (source-proof.md s4.2) while allowing classes with up to 2^33 segments (s3 under H2), so the index overflows on long classes. This package amends the encoding explicitly: segment-index field is 33 bits in the pool record (62-bit base + 33-bit segment + 56-bit link = 226 of 256 bits, within 31 spare bits) and the 128-bit transport chain id (62-bit base + 33-bit segment + spare = 95 of 128 bits). The segment directory already holds 2^33+1 eight-byte entries, and the 62-bit base index is unchanged (at most 2^61.6 bases). Zero instruction cost: same-word stores/loads; the re-walk index extract widens an immediate only (ubfx width 32 to 33, one instruction); record sizes, memory, occupancy and cap analysis unchanged.

## 8. Production lane-stack layout (128-byte slots; 192 since s30)

Supersedes the source's 96-byte lane-stack stride. Slots are 192 bytes (s30; 128 at k7), contiguous from a 16-aligned base: x [0,48), y [48,96), keys [96,160), pad; every listed ldp/stp keeps form and offsets, slot addresses use one hoisted madd. The s40 full kernel writes [out,out+96), leaving the key region to beta_keys; gap-table entries are 192 bytes (x, y, v, n; s1), read-only to kernels; batch buffers are separate (144 bytes since s32). Memory term LANES x (2 DEPTH + 3) x 192. kernel_check.py's stride audit (s1): 8 out slots in one contiguous buffer under random canaries, outputs exact, write footprint exactly the out ranges plus the frame; the current 80-byte-stride negative control trips it.

## 9. Record-path re-derivation: pointer buffers, honest matching, explicit setup

This section re-derives every per-record instruction from scratch, superseding the source s5.2 transport counts. The charged record sequences (GEN12 s10, RELAY s9, INS s22, LOOK/EXT3 s23, CMP7 s25, loops s29) are asserted by evidence.py audit_record_path (lengths, bits, immediates, pool fit, division proofs); the proto toy executes the pointer mechanism on every record buffer (pre-write check, +40 advance, exact derivation, reset).

Pointer buffers. Sender per-link state becomes (ptr, end): each record's three stores post-index ptr by 16+16+8 = 40 bytes, deleting slot-address madd and counter increment. Entries stay 16 bytes (one ldp; source "24-byte" prose counted a hoisted cap word), addressed by one add-with-shift (t/u below 2^30, lsl #4). The cap check precedes the write (cmp ptr, end; b.hs fail): identical failure event, never an out-of-bounds write. GEN listing: s10. The relay output side matches GEN's; relay input loads unchanged; s24 pairs the relay loop. At wave end the sender derives count = (ptr-base)/40 (sub, lsr #3, umull/lsr /5 with magic 0xCCCCCCCD, exact by m*5 = 2^34+1 plus asserted bound, count store, ptr reset: 8), writes the per-parity count word, and resets; receiver side (3) and double-buffer discipline unchanged. LINK_WAVE is 16 (8 sender + 3 receiver + 5 margin); the per-node term rises to NODES*24 for dmb, counter traffic and constants (16 used, 18 with s12 EXT traffic). Chain-id words are segment-invariant, formed once per segment by the formation row (register-holding superseded by s10; GEN reloads them).

Honest matching: this section's record key, pool layout and INS/LOOK/compare/cap counts are superseded by s11, s12, s17, s19 and s22-s25, s29 (current bodies in evidence.py, executed by audit_ins12/audit_look6). Four new explicit rows: window-table init (8/link, 4 used), chain-id formation (3/segment), candidate-link wave counts (9 on 4 traversals: /40 magic, /64 shifts, /56 with proven /7 magic), candidate formation (16/match, 11 used). Caps, routing, counts, events unchanged, so H1-H6, N, success, memory untouched.

## 10. k9 repair: chain-word memory schedule

This section repairs the k8 verdict's chain-register gap (s9's GEN11 and register-held chain words). Two chain-id words formed once per segment into a per-lane 16-byte slot (LANES slots per generating node, in memory term); every GEN reloads them. Formation 4 per segment (ubfx+bfi+lsr plus stp); GEN11 becomes GEN12 (ubfx t, entry address, ldp entry, ldp chain words, cmp, b.hs, three stores, ptr store-back, key-reload/loop pair). No chain word crosses records in a register (except within one giant point's three consecutive GEN copies, which no call separates; s34), so lane interleaves and side-output clobbering are compatible; allowance untouched.

## 11. k10 heads-x2: 2^58 slots per owner

This section doubles the per-owner head array, superseding the s9 57-bit head width and the head-slot/H3 sentences of source-proof.md s4.3/7. Head is key bits 60..117 (58 bits, 4+54 across w0/w1; lsr + combine bfi: width 53 to 54, same two instructions); compare key 74 bits (k10 bits 118..127 + k64). Extracts and pool layout: s17, s19. cmp_cap is 32312633143699954 per owner at the s2 sizes, i.e. 1.00001 x 3 N_g N_b/(OWNERS^2 x HEADS) since s30 (3 N_g giant lookups; was 1.01x). Support (exponents asserted): an owner's comparisons sum the chain lengths of its looked-up slots. Under H3 distinct giant points' keys are independent uniform and independent of the baby layout; a point's three keys are fixed bijections of x', so its k_g <= 3 lookups at the owner form one block in [0, 64 k_g], independent across points given the layout. Chains exceed 64 with probability below 2^-200 (Poisson(0.23) tail since s30, below 2^-324 over 2^118 slots; in the failure union). k_g >= 2 has probability 3 x 2^-120 per point (mean 0.14 per owner), so sum k_g^2 <= G_o + 1800 except with probability (e 0.14/300)^300 < 2^-2000. Hoeffding over blocks at relative deviation 9.7x10^-6: exp(-2 (9.7x10^-6 E_o)^2 / ((G_o + 1800) 64^2)), exponent above 94 nats/owner (2^-75 over 2^60 owners). B_o is Binomial(N_b, 2^-60); G_o = sum k_g is Binomial(N_g, ~3 x 2^-60) plus at most 600; both means exceed 5.5x10^16 and are exceeded by 10^-7 relatively with Chernoff exponent above 180 nats; 2 x 10^-7 + 9.7 x 10^-6 < 10^-5. The same blocks cover s12's lookups. All three sit far inside the retained 2^-40 budget. Toy support model-level, unchanged. Owner init doubles (2^57-stp stream, control below 2^99); (v) stays below 2^114.1 (s3).1.

## 12. k11 conditional extracts at an occupancy cap

This section charges the key extracts only on nonempty fall-throughs (EXT, listings since s23), at a per-owner occupancy cap. The counter lives in a register across the owner's lookups with wave-end traffic in the link row (s9 pattern; link-row use 16 to 18 of 24); exhaustion aborts into the tail budget: cap-plus-tail, never expectation. Per-owner cap is ceil(1.000031 x mu_L x p) = 29182932951402857 (at the s2 sizes; was 1.02) with mu_L = 3 N_g/2^60 exact (s30), p = y-y^2/2+y^3/6 at y = N_b/(OWNERSxHEADS-1): under H3 each key is uniform over the OC cells, so P(nonempty) = 1-(1-1/OC)^N <= 1-e^-y <= p (ln(1-x) >= -x/(1-x); decreasing alternating Taylor, asserted), and E[nonempty] <= mu_L x p by the tower argument over independent placements. Three 1.00001 stages cover the cap (1.00001^3 = 1.0000300003 < 1.000031, asserted): lookups/owner (Binomial, Chernoff), occupied/owner (negatively-associated, Chernoff, Dubhashi-Panconesi), probes-given-both (Binomial, stochastic dominance). Per-owner means exceed 1e16 (asserted: mu_L, occupied floor HEADS x 0.17); at delta = 10^-5 every Chernoff exponent mu delta^2/3 exceeds 3x10^5 nats, far above the asserted 70.5-nat budget (2^-40 over 2^60 owners), so the success union stays 7x2^-40. The toy counts nonempty probes per lookup phase. EXT row 3 x OWNERS x cap (EXT3, s23); capped EXT3 beats 1/lookup (asserted).

## 13. k12 77-bit restarts (price superseded by s15)

Pointer: restated and priced in s15 (SMX 630722); strand starts, base points, re-walks and candidate tests stay at full SM.

## 14. k14 contiguous k10 (layout superseded by s17)

Pointer: s17, s19 and s21-s23 supersede it.

## 16. Augmented planner (pointer; floors and planner since s30)

The s16 augmentation (1/3/5 alternating paths, partner cap 32, length-5 cap 8, exact rollback, OPS_BOUND = 10752 per element, asserted) runs unchanged inside the s30 planner, over centre pairs and in the fallback; its floor and tables are superseded by s30.

## 15. Exact restart price SMX from the shipped kernel bodies

Supersedes the s13 price; only the restarts row changes. Each restart computes [b]Q_a, b in [1, X2] (<= 77 bits, from the plan entry), by branchy Jacobian double-and-add (dbl-2007-bl, 2M+5S; madd-2007-bl with affine Q_a, 7M+4S) on an all-ones worst path in Montgomery residues. Strands, base points and candidate tests stay at SM = 8192000. SMX = 48 + 380 + 5 + 76 x (2340 + 11 + 3656) + 6 + 172520 + 1231 = 630722, an upper bound on every b; row count x (SMX + 16). Checked by smx_audit.py.

## 17. Permuted record key: the head index is one instruction

This section supersedes the s14 key layout (now a pointer) and restates INS/EXT/CMP. The 192-bit record key is permuted, keeping every field width: owner bits 0..59, k10 bits 60..69, head bits 70..127, k64 bits 128..191 (60+10+58+64 = 192, asserted). The head index is now w1 >> 6, a single listed lsr, so the 2N-instruction head-index word combine row is deleted. Head is still 58 bits (asserted), so HEADS, the s11 tail, caps, N, success and memory stand; owner routing (bits 0..59) is unchanged and H3 statistics are identical.

INS keeps its k10 count (insert body since s22): extr k,w1,w0,#60 (low 10 = k10) replaces the ubfx; bfi c4,k,#54,#10 deposits k10 into w4[54..63] (chain data [0..30] intact, spare [31..53] untouched; formation zeroes high as in s14). EXT3 keeps its count: subs, b.lo, extr (w0/w1 intact at EXT: the lookup lsr writes a scratch, others write no key-word register). CMP keeps its k10 count: and b,c4,#0xFFC0000000000000 (10 contiguous ones, asserted) replaces the ubfx; cmp b,k,lsl #54 compares full words (shift drops upper bits; both sides zero outside 54..63); ccmp on k64 unchanged. Register schedule as in s14 (extracted k10 + k64 stay in registers across the chain walk). Counts: s23.

Semantics audited, not just counted: audit_packed_match emulates the full 64-bit extr/and/shifted-cmp against the layout oracle (2^60 boundary; 10 straddles at 60/63/64/69/70; 1500 random keys with extract/head-lsr = oracle and spare-only writes; 1500 keys with match, k10-flip mismatch, other-flip inertness, head-flip chain change). The toy routes every comparison through the packed ops, asserting packed iff full-key decision. The permutation idea follows the public note of gbotrel's 12528fa8-ce59-4797-9e36-2849dab46c11 (s17 there); this section is an independent re-derivation (own verification, own audit code, no files taken).

## 18. Wave-granular pool-cap check

This section supersedes the per-insert pool-counter compare (N x 2, deleted) with one check per owner per wave. An owner inserts only through its 2^20 inbound relay links (source-proof.md s4.3), each carrying at most CAP2 = 2^21 records per wave, enforced by the relay's already-charged pre-write cap check (GEN/RELAY cmp + b.hs; overflow is a run failure in the same budget). So an owner inserts at most 2^41 records per wave (asserted); checking counter <= 2^56 - 2^41 at each wave start guarantees no wave overflows the pool (start bound + one wave = cap, asserted). The failure event (count above threshold at some wave start) supersedes the old overflow event, taking the same H3-tail union member: union and success formula stand. Tail: mean N_b/OWNERS exact, excess d in (0,1), Chernoff mu d^2/3 above the 70.5-nat budget (actual above 2.7x10^14 nats at s30's N_b; asserted).

Schedule, all counted (counter = pool pointer since s22): counter and threshold live in pool slot 0 (allocated, unused by records; zero memory delta). Per owner per wave: ldr, ldr, subs, b.hi at wave start plus str at wave end (5; no cross-wave register lifetime). Once per owner: movz counter, threshold in movz + 3xmovk (asserted), stp (6). Row = OWNERS x (5 x WAVES + 6) = 32160694642461793694253056 (asserted exact; WAVES with s30's 2^32-step giant waves). Both phases charged though giant waves never insert (conservative); preprocessing takes the row in full; wave-loop control (~2^82) in the (v) bound. After the public note of 12528fa8 (s18 there); this section itemizes 5+6 with nothing absorbed.

## 19. Full-address chain links

This section restates the link word (superseding the s7 56-bit-link sentence): head entries and chain links store full 64-bit record addresses, 0 = end or empty; records load as [addr,#imm] (comparison CMP7 since s25, lookup/EXT3 since s23, head store in the s22 loop). Layout: pool after the head array at nonzero base HEADS x 8, records at base+32 and above, addresses below 2^63 (asserted). Pool use 233 of 256 bits (64+64+10+62+33, asserted); word sizes, memory, owner init unchanged. Links never leave the owner (chain ids are base+segment).

## 20. Itemized side-call fixed price: 40 per invocation

This section replaces the s1 64/invocation bound with a complete per-call schedule: 32 fixed (prologue 22 + epilogue 10 with ret, asserted every call) + 1 exit branch (measured) + 5 staging (mov x0,x_args; mov x1,x_n; lsl x2,x_n,#5; stp xzr,xzr,[x_args,x2]; bl) + 2 round control (subs, b.ne) = 40 (asserted). The per-lane argument records sit in PAIR_BK/ROOT_BK (s26), so staging is exactly these 5. Zero spills: 8 loop variables in x19-x28 across the call; the kernel preserves callee-saved regs (asserted); nothing live in caller-saved regs across bl. Lanes are unrolled x n (n fixed; PADOPS pads to full calls), so no per-lane loop control exists; the 40 is n-independent anyway. Driver and kernel share one small text (bl in range, no veneer). Per-base-per-phase setup bounded by 32 (entry/exit 10 + materialize 8 + base-loop 2 + 12 margin; true cost <= 20), charged as 64 x (BASES+2). Row: s1 (eight calls per round since s30; 48 per call since s37), both phases charged. gbotrel's 03e66f0b-6d4a-4c97-bbdf-fe9963d0579e independently itemized 40/invocation (convergent corroboration, not a premise).

## 21. Buffer-granular EXT-cap check

The per-lookup b.lo moves to a buffer check (EXT3 since s23). A fall-through lookup runs subs + extr only; before each link buffer the owner checks remaining >= 2^21 (cmp-immediate + b.lo: 2). A buffer holds at most CAP2 = 2^21 records by the relay's charged pre-write check and EXT runs at most once per lookup, so a buffer adds at most 2^21 EXT (hard) and the total never exceeds cap_ne (asserted). The failure event (above cap_ne - 2^21 at some buffer start) supersedes the old one and takes the same s12 union member: cap_ne - 2^21 still covers the three 1.00001 stages (asserted, margin 2.8x10^10 at s30's sizes), so tails and success stand. Schedule, all counted: EXT counter in pool slot 0 after the s18 pair (24/32 bytes, asserted; zero memory delta); per owner per wave ldr + str (2); per buffer cmp + b.lo (2, in the buffer loop, control in link row); one-time movz + 3xmovk + str (5, asserted). Row = OWNERS x (WAVES x (2 + 2 x 2^20) + 5) = 13489182981024494619417473187840 (asserted; s30 WAVES); all waves charged, giant only needed (conservative); preprocessing takes it in full. Toy enforces per-link pre-write caps; thresholds are ledger asserts.

## 22. Pool pointer: insert pair loop

This section supersedes the s17/s19 insert. The owner's pool counter becomes a pointer p to the next free 32-byte record, deleting the slot-address add and the counter add (the s9 sender device at the owner). The insert body in order: ldp x0,x1,[x20],#16; ldp x2,x3,[x20],#16; ldr x4,[x20],#8; lsr x9,x1,#6 (head); ldr x11,[x21,x9,lsl #3] (old head); extr x10,x1,x0,#60; bfi x4,x10,#54,#10; stp x3,x4,[x22,#16]; str x22,[x21,x9,lsl #3] (head = record address); stp x2,x11,[x22],#32 (k64, link = old head; the only write of p, after both stores through it). s24 runs 64 bodies per iteration sharing one subs + b.ne (s29: 10 + 2/64 per record); setup and tail ladder add at most 10 per invocation. Record words, chains, links, comparisons and events are unchanged. s18 on the pointer: pool base PB = HEADS x 8 = 2^61, first record PB + 32, threshold TP = PB + 32 x THRESH = 2^62 - 2^46; p <= TP at a wave start plus one wave (at most 2^41 records) ends at PB + 32 x POOL = 2^62, the pool end (asserted); the wave check compares p with TP (same five instructions, same failure event); p and TP each initialise by movz + movk inside the charged 6 (asserted). audit_ins12 executes the s29 loop with a64.py on 600/601/1/0/639/38 records over 40 heads: 10n + 2(n div 64) + 10 steps, final p, and memory exactly the oracle's (records, links, heads; nothing else written).

## 23. Deferred k64: lookup pair loop + EXT3

This section supersedes the s12/s21 lookup and EXT. The lookup body is ldp x0,x1,[x20],#40 (w0,w1; ptr past the record), lsr, ldr head, cbz (empties skip free, ptr already advanced); fall-through runs EXT3 (subs, extr k10, ldur x2,[x20,#-24] = w2 of the passed record). s24 runs 64 bodies per iteration sharing one subs + b.ne (s29: 4 + 2/64 per lookup plus EXT); setup and tail ladder add at most 10 per invocation. LDP #40 / LDUR #-24 are one instruction each (ranges asserted). w2 feeds only CMP (s17 schedule); buffer, heads, events unchanged; s21's buffer cap still bounds EXT. audit_look6 runs the s29 loop in a64.py on 600/601/1/0/639/38 records/40 heads: 4N+2(N div 64)+3K+10 steps, final ptr, x24, single-iteration regs, untouched mem; wrong advance, offset or lsr and dropped ldur each fail it.

## 24. Transport loop unrolling (64-fold since s29)

This section supersedes the loop control of the s9/s10 RELAY14 listing and the s22/s23 single-record loops (bodies unchanged). The loop pair (subs + b.ne) is each transport loop's only cross-record instruction, so every loop runs UNT = 64 body copies back to back sharing one pair (s29). Each copy reuses the previous copy's registers after its stores complete, in order; pair-control scratch x5 is dead in the exhibited bodies (asserted absent) and exists for RELAY by the s9 dataflow (5 input words, routing, entry, ptr/end: at most 13 live of 31). Setup is and x5,x23,#63 + lsr x23,x23,#6 + cbz; the tail is a binary ladder: for j = 5..0, tbz x5,#j skips 2^j body copies. Per record RELAY 12 + 2/64 (body asserted verbatim vs the s9 mnemonics), INS 10 + 2/64, LOOK 4 + 2/64 plus unchanged EXT3; a tail record costs its body only. Per invocation at most setup 3 + six tests + ret = 10 (UNROLL_FIX, asserted), incl. empty loops. GEN12 has no per-record loop pair (driver control in s3(v), unchanged); CMP7's cbnz walks chains of mean length below 1; EXT3's subs is a per-item counter; lanes and owner-init are already unrolled (s20, s3(v)). Charged: 770/64, 642/64, 258/64 per record (average 497/16, asserted) plus the unroll row 10 x ((LINKS1+LINKS2) + 2 x OWNERS) x WAVES (s4; both phases charged, one used). No new heuristic or premise change: H1-H5, N, success, memory untouched.

## 25. Paired k64/link comparison load: CMP7

This section supersedes the s19 comparison listing. The executed insert stores k64 then link per record: `stp x2, x11, [x22], #32` (k64 +0, link +8; audit_ins12), so the first two comparison loads fuse: CMP7 is ldp,ldr,and,cmp,ccmp,b.eq,cbnz (ldp pair, no writeback; ldr [#24]; rest as s17/s19). Pool base 2^61 is 16-aligned (asserted), records 32 B, so every ldp address is aligned; the s24 non-unrolling reason stands. Same words/decisions/events/caps; H1-H5, N, success, memory, preprocessing untouched. Row = OWNERS x cmp_cap x 7 = 260777507352905359601837545753673728 (asserted; s30 cmp_cap).

## 26. Batch round: both children by the derived kernels, no reference addition in the walk

Supersedes the s1 pair/root prices and the pair, root and inversion-share items of source s5.2 (and their s8 restatement); only the records, padding, partial-plan and side-call rows change, plus the memory row. Points, keys, records and events are unchanged (canonical residues are unique): H1-H5, N, success stand.

Round (one step of the 2^18 lanes): the driver writes argsB[i] = {out, P, T, buf+144i} per lane op (out = the P+T child, or the root's point) and argsA = {P-T out, P, T, buf+144i} per pair (since s30 each split into a full and a leaf array by the leaf flags), each ended by a 32-byte NULL record; padded lanes use the fixed pair (Q_a, T[2]), x distinct, so every denominator is nonzero (n = 2^18, a multiple of the s28 unroll factor 32). batch_round(argsB, buf) (batch_audit.py, generated from the shipped MULNR body and s15's C2/Fermat): forward d_i = x1 + v2 in (0,2p) (v2 = p - x2 at T+96; s27) at buf_i+48 and c_i = S(c_(i-1), d_i) at buf_i (c_(-1) = R mod p at buf-144); the canonical c_(n-1) is inverted by s15's Fermat and negated once (t = p - t); backward inv_i = S(t, c_(i-1)) at buf_i, unreduced (s27), then t = S(t, d_i). side_lanes(argsA) then writes P-T and side_lanes_plus(argsB) writes P+T: side_lanes with its y2 loads (T+48/64/80) and n2 loads (T+144/160/176) exchanged (9 ldp offsets), i.e. side_lanes on -T = (x2, p-y2, v2, y2), so every kernel bound holds verbatim. With b < 2p and p/R < 0.1016 every MULNR operand a stays below 1.26p (K1 holds) and inv_i < 1.17p, congruent to 1/(x2 - x1) (s27).

Counts (current: s28 loops, s32 body): per element forward FB + backward 595, one cbnz and one cmp + b.ne per 32 elements, BATCH = (32 FB + 1 + 32 x 595 + 2)/32 per lane op, the only path. Round fixed 172585 (with the negation) plus 5 per call (two moves, bl, terminator address, stp) = 172590 per round of exactly 2^18 lane ops (s20), charged as the share s = 172590/2^18 per lane op (s28). PAIR = 2 KERNEL + BATCH + 28 + s for two points (plan entry 8, P address 1, T[g] address 6, buffer 1, two out addresses 2, two argument records 4, s30 routing6); ROOT = KERNEL + BATCH + 21 + s (gap entry 5, T[gap] address 6, P and out addresses 2, buffer 1, argument record 2, loop 5). P/T and inverse loads and coordinate stores are inside the executed listings. Side-call row: s1 (48 per call since s37). Re-walks keep 3300 + 232864.

Audits: batch_audit.run (every ledger evaluation) executes one round of 64 lanes of genuine points (s28): every inverse is congruent to the residue of 1/(X2-X1) and below 1.17p, one unreduced (s27), sp and x18-x30 restored, inputs untouched, writes only buffers and frame, count ROUND_FIXED + 29411n/32; both side calls give exact P-T/P+T coordinates with no parity writes; the composed beta pass verifies routed keys. kernel_check.py plus mode runs the 509-lane suite on side_lanes_plus with negated entries (the same operands at every instruction): 12 paths, worst916, all7 branches both ways, 24 genuine pairs, stride audit, fixed 32.

The reference addition (2000: 3M formula plus 3M batch share, A = 0) is no longer used anywhere in the walk. P+T, the batch products (three per lane op, 3n >= 3(n-1)) and the inversion are other work charged at their listed and executed counts, the rule s1 already applies to P-T; nothing is charged twice or omitted.

## 27. Lazy batch round: reflected difference and unreduced inverses

Changes only the batch listing of s26 (measure/side_lanes.S unchanged, SHA-256 b4c7a7c2... as in measurement.log). Prices only: points, keys, records, events, H1-H5, N, success and memory stand.

Reflected difference. The forward pass forms d_i = x1 + v2 from v2 = p - x2 at T+96, present in every entry the round reads (gap entries, base points, the T[2] pads; s1), by six adds; d_i is in (0, 2p) and d_i = -(x2 - x1) mod p, nonzero as x1 != x2. Montgomery's trick is linear in its running value t, so the round negates the canonical Fermat output once (t = p - t in (0, p), six instructions in the round-fixed part); every backward t and inv_i = S(t, c_(i-1)) is then the negative of its value for d_i, so inv_i = -1/d_i = 1/(x2 - x1) mod p, the inverse both kernel calls need.

Unreduced inverses. The C2 after inv_i is deleted. With t, c < 1.26p (s26), inv_i < 1.26^2 x 0.1016p + p < 1.17p. The kernel reads inv only as the b operand of its first product, and K1 (s1) holds for every 384-bit b. Its bounds become: t < 2p; mu < 2.34 x 0.1016p + p < 1.24p < 2^382; V < 1.24^2 x 0.1016p + p < 1.16p; z < 2 x 1.24 x 0.1016p + p < 1.26p; X in (-p, 2.16p) and D in (-p, 1.26p) stay inside the existing reduction arms (X < 3p, D < 2p), so x',y' stay canonical; s40 removes unused parity output. The square needs 2c < 2^704 for its cross sum c: mu < 2^382 puts the top word below 2^62, so c < 2^702 + 2^639. The comments in measure/side_lanes.S state the canonical-inverse bounds (inv fully reduced, mu < 1.21p < 2^381, V < 1.15p, z < 1.25p); this section supersedes them, and no instruction changes.

Evidence: the batch-round audit (unreduced inverses exercised and asserted) is in s26/s28. kernel_check.py (plain and plus) runs lanes with inv in [p, 1.17p), including inv = 1.17p - 1 at t = 2p - 2; nocarry_witness.py covers b up to 1.17p and 1.24p and square inputs up to 2^382 - 1.

## 28. Exact round share and unrolled batch loops

Changes the s26 round share and the s26/s27 batch loop control only. Exact share: a walk round is one lockstep step of a node's 2^18 lanes, each doing one pair or root operation, idle lanes padded and charged (source s2 step 4, s3 padding lemma; s20): exactly 2^18 lane ops and one fixed cost 172585 + 5 = 172590 (s26), so each lane op carries s = 172590/2^18, not 1. Every executed lane op sits in a row charged at PAIR or ROOT (records at the floor F, which charges at least the lane ops of any actual f >= F; padding; remainders), so the shares cover every round. Table strands keep their rounded share.

Unrolled loops: the forward loop's only per-element control is cbnz x0, the backward loop's cmp x2,x30 + b.ne. Each loop holds 32 copies of its unchanged body and one control; bodies address only through x30 (pre-indexed) and x2 (post-indexed); n = 2^18 is a multiple of 32, no tail (asserted). Counts and listing size: s32 (inside gen_mem's 2^20-byte code term). Evidence: batch_audit.run executes one 64-lane round, each loop branch both ways.

## 29. Transport loops 64-fold with a binary tail ladder

Changes only the s24 unroll factor (2 -> 64) and tail (cbz + one body -> six tbz over 32, 16, ..., 1 bodies); bodies, records, chains, events, H1-H5, N, success and memory stand. evidence.py generates the loops (unroll), so no static copies ship. Average record transport 497/16, baby phase 2180/64 per record; rows: records, candidate delivery (eight relay hops at 770/64), unroll setup and tails (10 per invocation). Evidence: audit_ins12 and audit_look6 execute the generated INS and LOOK loops with a64.py on 600/601/1/0/639/38 records (every tail test both ways): exactly 10n + 2(n div 64) + 10 and 4N + 2(N div 64) + 3K + 10 steps, memory equal to the oracle (INS) or untouched (LOOK), single-record registers; the four LOOK negatives still fail.

## 30. x-only leaves, beta-triple giant lookups, leaf planner and asymmetric sides

Supersedes the record key (source s1 parity key), source s3 matching, the s16 floor, GMAX 2^44 -> 2^47 and DEPTH 2^14 -> 2^17 (source s2-s3), N of s2, the 128-byte slot of s8, the per-point price of s26-s29, the two side calls per round of s20/s26 and the 2^34-step giant wave (source s4.3). Everything else stands.

Keys. Every record key is limbs 1..3 of the canonical Montgomery residue of x' (out+8..31). [-1](x, y) = (x, -y) and [lambda](x, y) = (beta x, y), beta fixed by [lambda]P0 = (beta x0, y0) (asserted in beta_run), so the U-orbit of Q has x-set {x, beta x, beta^2 x}; conversely x_P in that set gives y_P^2 = x_P^3 + 4 = y_Q^2, so P is in U Q. A baby x-key equal to one of a giant's three keys is exactly the orbit match the parity key tested, and at most one of the three can match (x = beta x forces x = 0, a 3-torsion point; 3 does not divide r). Babies insert one record per point; giants send three lookup records per point (keys of x', beta x', beta^2 x', one chain id; GEN unrolled three times per point at key offsets +8/+104/+136; a giant re-walk keeps positions where any of its three keys equals the candidate key). Coverage events, K, success and the 96 recovery tests are unchanged. Giant waves are 2^32 steps, so a generating-node link still carries at most 2^23 records per wave in the mean (six per step).

beta_keys (batch_audit.beta_run; tables since s38): per giant slot, limbs 1..3 of canonical beta X' at +104 and of beta^2 X' = -(X' + beta X') at +136, walking the side-call argument arrays (the driver stores nothing extra).

Leaves. A pair child that is no pair's midpoint is only keyed: it runs the leaf kernel kernel_check.leaf_text(), the kernel cut after the x' store (y' tail and mu spill deleted; prices since s32, s37), bit-identical in x', executed by run_leaf (365 lanes each copy). Midpoints and roots run the full kernel. Slots are 192 bytes: x' 0, y' 48, keys 96..159.

Driver and calls. Plan entries carry two leaf flags (kind 2 + gap 48 + three 19-bit slots + 2 = 109 of 128 bits); per pair lane op two ubfx and two pointer ldr/str route the children to the full or leaf call (ROUTE = 6, listed). Per round four side calls (full/leaf x minus/plus) and four beta_keys passes, charged 8 x 40 on both phases (s1, s20). Pads run full pairs plus two beta passes.

Planner (planner.py). One centre per run of 50 class-list elements (the one with most candidate pairs); greedy by increasing degree over centre pairs, hosting centres first; the s16 augmentation; a fallback pass over all midpoints; source s3 demotion and emission; every pair passes the shipped cycle, depth and tree rules with exact rollback. GMAX = 2^47 = 128 d^2 gives each midpoint about 128 candidate pairs; WIN = 2^29 keeps rel32 offsets at most 2^30 < 2^31. Charged: 6 window scans per element (was 3); matching at 64 kappa = 512 candidate pairs per element (twice the mean 2 GMAX/d^2), at most five tries each (centre pass twice, fallback three times); two augmentation passes, each with its partner cache; memory adds a 3-byte centre flag and count per element (planning node 2^62.79 < 2^63). Each class is planned before any group operation and the run stops unless f_c >= F = 998/1000 and the leaf share ell_c (leaf children / pair children) >= ELL = 193/200 since s36 (H4 restated; H4L, a cost floor). Measured at kappa = GMAX/(16 d^2) = 8, span X2/(2 GMAX) = 4: the declared run at d = 4; d >= 16 by a line-for-line C port of plan() that reproduces planner.py's forests exactly (pairs, roots, leaves) on six declared-type lists, capped and uncapped:

| d | classes | min F | min ell | mean ell | max height |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 4 | 2 | 0.99550 | 0.97349 | 0.97375 | (declared) |
| 16 | 4 | 0.99822 | 0.97278 | 0.97298 | 235 |
| 64 | 4 | 0.99863 | 0.97213 | 0.97296 | 931 |
| 128 | 4 | 0.99875 | 0.97305 | 0.97329 | 1062 |
| 256 | 8 | 0.99877 | 0.97306 | 0.97324 | 2172 |
| 512 | 4 | 0.99886 | 0.97311 | 0.97338 | 6931 |

ell is flat in d (0.972-0.974) and F rises with d and span; F = 998/1000 sits 0.00077 below every d >= 256 row, ELL (96/100 at k36; 193/200 since s36) below every row. Tree cap: at TREEMAX = d (the production cap/window ratio 2^-9) the leaf deficit is about 5/d (4.5-5.1 at d = 32..128) and the F excess 1.1/d to 2.2/d (d = 16..128); at fixed d = 64 the F excess falls 0.0193, 0.0100, 0.0027, 0.0011 for caps 64..4096. Both routes put the production cap effect below 10^-4 at TREEMAX = 2^20, inside the margins; that extrapolation is part of H4/H4L. Heights grow with tree size (6931 at 2^19 elements): DEPTH = 2^17, stack 2^18 + 1 slots.

Sides. With side prices c_b, c_g the optimum is N_b/N_g = c_g/c_b (1.169 at s31's 23/20). H1 keeps its form; census minus 1 - e^-K (m = 2^24, d_t = 4096, seeds 0-3): -0.00046, +0.00004, +0.00008, -0.00023 at 23/20. The ec toy (x-keys, triple lookups, solved == exact coverage every trial, all scalars correct): declared 27/60 at 0.3956 predicted (23/20). Its rate is a binomial sample of toy coverage at K = 0.50 (SE 0.063 at 60 trials, 0.028 at 300): it checks that success is exactly the coverage event; production success rests on H1 and the census (s2, source s7), not on toy rates. Offline, seeds 7-12 x 300 trials: 713/1800 = 0.396 at 23/20, 715 at 27/20, 717 at 1:1 (model 0.396-0.397). Caps re-derived with 3 N_g lookups and N_b inserts (s11, s12, s18, s21 asserts pass): pool mean N_b/2^60 < 2^55.74; false matches 3 N_b N_g/2^192 < 2^40.85 < 2^44 (Chebyshev < 2^-40). H3 now covers N_b + 3 N_g keys: a giant's three keys are images of x' under a fixed bijection, treated as three uniform keys; the declared toy's comparisons (28.9 vs 27.8 model) and candidates (0.85 per trial vs 0.217 false-match model plus true hits) follow the three-lookup model; the re-walk step price 3300 covers 2 x 929 + 2 x 343 and key compares.

Prices (pinned in audit_record_path; current values s1): pair child (BATCH + PAIR_BK + 6)/2 + ELL x KX + (1 - ELL) x KERNEL; baby point F x child + (1 - F) x ROOT + 2180/64; giant point F x child + (1 - F) x ROOT + 3 x 1796/64 + BETA. Records = N_b x baby + N_g x giant.

## 31. Historical 23/20 tuning, superseded by s39

Re-tunes N_b/N_g to 23/20 (optimum 1.169); keys, records, events, H1-H5 and H4L stand. Tables, lookup and evidence: s38 (k37k's in e5e8ae33).

## 32. Shared x-sum e = v2 - x1, once per lane op

Prices only. Both children of a lane op formed X = V + v2 - x1 from the same P, T (reloading x1, v2; 12 adds/subs). The forward body, which holds x1, v2 for d, now forms e = v2 - x1 (subs, 4 sbcs, sbc into the dead x0, x1) and stores it at buf+96 (three interleaved stp): FB 324, BATCH 29411/32. Slot 96 -> 144 bytes (immediates only: strides 144, c_(i-1) at x2-224; memory (LANES+1) x 272; listing 201996 instructions, inside the 2^20-byte code term).

Kernel: kernel_check.esum_text replaces one asserted block of the measured listing by `ldp x29,x0,[x30,#16]` (T, buf), three ldp of e, the v2 loads, X = V + e (adds, 5 adcs), `tbz x16,#63,.Lx_nonneg`. As integers e is in (-p, p), so X is the old integer, in (-p, 2.16p), inside (-2^382, 2^383): bit 383 of its two's complement is set exactly when X < 0, the old b.hs condition; x', y' are bit-identical (x19..x24 now hold e, rewritten before any read). KERNEL 923 (24 paths); the s30 leaf cut also drops the v2 loads: KX 587 (4 paths). At that historical rung, the ratio stayed 23/20; s39 now supersedes it.

Historical s32 evidence: kernel_check.run plain and plus (509 lanes, e at inv+96: exact, worst923, all branches both ways incl. the tbz, fixed 32, stride audit); run_leaf plain and plus (365 lanes, worst 587); batch_audit.run (64 genuine lanes: e = v2 - x1 mod 2^384 at every buf+96, inverses below 1.17p, both side calls exact; 10369/19042 per 32). Callgrind (s1) measured the 929 listing; 923 and 587 are derived listings priced by static path enumeration and executed by a64.py, like s30's leaf.

## 33. Gap tables by pair levels

Changes only how tables T[k] = [k]Q_a, k <= GMAX, are built (entries, layout, memory and use unchanged). With T[1..M] known (M = 3 from the strand-start scalar multiplications, row unchanged), one level computes the centre c = M + k + 1 as a root T[M] + T[k+1] and then the pairs T[c -+ j] = T[c] -+ T[j], j = 1..k, k <= M - 2: the s26 lane op on (T[c], T[j]) with both children full (916 kernel, minus and plus calls, then the 21-instruction v/n step), so T[1..M + 2k + 1] is known. Full levels (k = M - 2) give M' = 3M - 3; the last stops at GMAX: 30 levels. Incomplete additions are safe: c - j, c + j and c are distinct, nonzero and below GMAX < r, and x(T[c]) != x(T[j]) as c != +-j.

Price per base: every level is charged ceil(k/LANES) + 1 full rounds of LANES lane ops at PAIR (+1 for the centre; padding included) plus 128 per round for calls, and 21 per entry: BASES x (rounds x (LANES x PAIR + 128) + 21 GMAX) (strands: BASES x (GMAX - 1) x 2050; re-walks keep 2050). Row in s4.

Evidence: batch_audit.run(table) executes a 64-lane level (M = 66, centre [131]Q, gaps [j]Q): the s32 round and both kernel calls give exactly [131 -+ j]Q; evidence.py asserts the 30-level recursion.

## 34. Giant GEN chain-word reuse; head zeroing by dc zva

Prices only. (i) A giant point's three lookup records share one chain id (s30; GEN unrolled three times per point, key offsets +8/+104/+136). Copies 2 and 3 drop their chain-word ldp (copy 1's registers still hold the words; GEN only reads them). Giant transport 3 x 1796/64 - 2; the proto's generic per-record GEN 12 stays an upper bound. (ii) Owner initialisation zeroes each 64-byte-aligned head array by dc zva (64-byte blocks, DCZID_EL0.BS = 4, common AArch64) plus one add: 2 instructions per 64 bytes, not one stp per 16. Row OWNERS x (HEADS x 8/64 x 2 + 64).

## 35. EMIN = 2^56 - 2^52 (H2 restated)

Bases are charged at N/EMIN per side; EMIN was 2^55, half the mean X2/d = 2^56. H2 now has class lists between EMIN = 2^56 - 2^52 = 0.9375 x 2^56 and EMAX = 2^57, still sieve-checked before any chain starts (a short class stops the run: an H2 failure event). Bases per side at the bound: less than 2^59.84 per side (I's 2^64.3 primes: margin above 20); segments per class at least 2^31.9. Only per-base rows change. Support: class-length spread shrinks like (X/d)^(-1/2). Executed (evidence.py h2_sizes, mode all) on toy primes, [d, log2(X/d), min, max of |L_c| d/X]: [16, 10, 0.9648, 1.0391], [16, 14, 0.9908, 1.012], [64, 12, 0.9714, 1.0295], [256, 10, 0.9287, 1.1035], [256, 12, 0.9646, 1.0317] (offline, d = 16 at 2^18: min above 0.994). At 2^56 per class (d = 2^20) this extrapolates below 10^-6, far inside 6.25% (heuristic).

## 36. Leaf floor ELL = 193/200

H4L's runtime-checked floor rises from 96/100 to 193/200 (a class below it stops the run before any group operation, as before). Per-class minima: s30's table (0.97213-0.97349); new offline planner.py rows (kappa 8, four classes): d = 8 span 4: 0.97274; d = 8 span 8: 0.9728; d = 16 span 4: 0.97289; d = 4 span 8: 0.97524. One smaller two-class d = 8 run (lists of 7625 and 8734) had a class at 0.9691; 0.965 is below every measured class and ell rises with span (production 2^28). Prices only.

## 37. Leaf reduction at a capped share; H7

The leaf's X in (-p, 2.16p) needs the second subtraction (+7; paths 580/580/587/587) only if X >= p; the leaf counts those lanes (ldr/add/str: 590) in its free slot [sp,#96], zeroed in the prologue, returned in x0 (fixed 34). The driver adds each return to a node count C and its lanes to L (inside the 48 now charged per call, 40 + 8). At every wave barrier the tree sums (C, L) (16 per node per wave) and the run stops unless C <= 7L/32. Between barriers C grows by at most one wave of leaf lanes, so a leaf child costs at most 580 + 10 x 7/32 = 582.1875 plus the row GENS x WSTEP x LANES x 2 x 10 + WAVES x GENS x 16.

**H7 (new, in the claim's conditioning).** At every wave barrier, the cumulative share of leaf lanes with X >= p is at most 7/32. A violating run stops (excluded by H7, as H2's sieve checks). Support: exact arithmetic on 2 x 10^5 random lanes: 0.1878. batch_audit.xge_run (declared; one base Q, gap table [j]Q, walked centres, exact s32 batch-round inverses, 16 rounds of 256 lanes, both children) gives 0.1874, per-round shares 0.152-0.232 at 512 flags. Its model equals the executed listing's inverses (asserted in batch_audit.run); a barrier share averages >= 2^99 lanes. run_leaf (plain and plus): the returned count equals the exact X >= p count (61 of 365), body exactly 580 n + 10 count.

## 38. Seven-chunk beta tables (k37l)

Replaces s31's tables, lookup and evidence. Canonical Montgomery X < p < 2^381 splits at bits 0, 55, 110, 165, 219, 273, 327 into three 55-bit and four 54-bit chunks c_i; T_i[c] = beta c 2^shift_i mod p, canonical. Tables 0-2 interleave at A + 256c + 64i, tables 3-6 at B + 256c + 64(i - 3), B = A + 2^63 (one add-lsl per address); q p (q < 8) at A - 512 + 64q; production A = 512. BUILD is s31's 24-instruction recurrence at stride 256: row GENS x (5 x 2^55 x 24 + 7 x 2000 + 64) in preprocessing; storage 3 x 2^62 + 512 bytes.

Lookup, 128 per giant point: 3 loads, 12 chunk, 7 address, 21 table loads, 36 sum, 23 reduction, 2 stores, 22 for beta^2 X, 2 loop; fixed call part 26 of the 48 charged (s37). S < 7p and q = floor(s5/(p5 + 1)) <= 6 give q p <= S < (q + 1)(p + 2^320) <= q p + 2p (7 x 2^320 <= p): one conditional subtraction; keys and offsets unchanged.

Addresses: each node is one flat 64-bit space (source s4.1). Tables fill [0, 3 x 2^62 + 512), the other generating-node regions (windows included) follow, total 17609105416545370912 < 2^64 (asserted); other nodes stay < 2^63. Addresses are unsigned and never wrap: the ledger asserts that no shipped listing or mnemonic list has a signed condition, sign-extending load or asr, and that bit 63 is tested only on data (the kernel's X sign, the restart scan's scalar words); caps compare with hs; rel32 offsets add modulo 2^64 within one region. Analytical, not hardware.

Evidence: beta_run (66 slots incl. every chunk boundary, 663 reduction cases, 448 BUILD entries, layout asserts); batch_audit.run composes batch inversion, both side calls and beta on the production 192-byte slots (outputs, keys, registers, memory, counts). Finite checks. H1-H5 and H7 stand. Historical s38 gain:0.005605324156009 bits; current prices are in s4/s40.

## 39. Inherited ratio and evidence repairs

The ratio is29/25, unchanged from3c7fc300. Choose the least N_g with
6*ceil(29*N_g/25)*N_g>=KT*m, then set N_b=ceil(29*N_g/25).
Section2 and the ledger give current integers. Mode scan checks six ratios;
it is design work, not an attack stage, and makes no global-optimum claim.
H1 is asserted at these exact sets and realized K; historical23/20 toy
observations are supporting evidence, not a proof for these sets.
Each side uses fewer than2^60 bases, giving6*2^120/m<2^-111 for diagonals.
WAVES=5578990<2^23; cap_ne=29182932951402857 and cmp_cap=32312633143699954.

Planner depth and height use32-bit fields for DEPTH=2^17:28 bytes base,
8 adjacency,4 previous link,3 centre, totaling43 bytes/element. Wider
loads/stores retain the instruction prices. The evidence parser caches
16 listings and returns fresh lists. Planner toy degrees are powers of two.
Mode all adds30 recovery trials; arithmetic audits retain200 random samples,
boundary cases and negative controls. The character check certifies order d
for7^m and exhibits k=256 where log_g(x^m) is the wrong character base;
source section2 uses the correct log_(g^m). Earlier public ratio-scan and
cache ideas from unpublished cf957fb1 are credited in the direct parent's
note; no new import from that entry occurs here.

## 40. Delete the unused full-point parity key

Section30 already replaced every record key by x-coordinate limbs at+8,
with giant beta keys at+104 and+136. Full internal points still need y at+48
for future additions, but no current consumer reads the old parity key
at+96. The full kernel first stores canonical x and y, then each of four
branch exits conditionally negates y solely to emit that discarded key.
drop_parity_text removes exactly those four blocks, checking each mnemonic
sequence: tbz, subs, sbcs, sbcs, sbc, stp, str. Every lane formerly executed
one block: three instructions on even y, seven on odd y. All coordinate
arithmetic and stores precede it and are retained.

The next-lane ldp/cbnz use no flags. The next lane reloads its scratch
registers and executes adds before any carry consumer; exit restores the
callee-saved registers. Thus deleting the blocks preserves the live state.
The full output footprint is now[0,96); no bytes at96..191 are touched.
GEN reads x at8..31. beta_keys reads x at0..47 and independently writes its
keys at104..127 and136..159 before giant routing. Gap tables consume x,y
and explicitly construct their v,n fields. No removed key is an input to
planning, matching, transport, re-walks, scalar recovery or verification.

The derived full kernel has 12 paths and worst916 instructions/lane,
versus24 paths and worst923. Fixed-call overhead remains32. The leaf
generator is byte-identical to its parent, asserted for both signs; its
590/580 paths and H7 counter are unchanged. KERNEL=916 prices full roots,
nonleaf pair children, padded pairs and gap tables. Rounded re-walk budgets
are retained. The per-point saving on each side is exactly
7*(1-(998/1000)*(193/200))=0.25851 instructions. No success event, emitted
point, size, resource cap or heuristic changes.

Evidence runs both full kernels against exact x,y arithmetic on509 lanes
each, with all remaining conditional-branch successors reached. The192-byte
canary audit checks coordinate-only writes; an80-byte slot must fail.
The composed batch checks both signs and absence of writes to the deleted
key region, then executes beta_keys and verifies every routed key. Both
leaf kernels match their parent text exactly. The raw measured assembly
and its log are unchanged historical inputs;916 is an analytically counted
derived listing, not a new native measurement.

The supplied-total rational certificate now brackets the least10^-15 grid
ceiling, using ceiling-10^-15<low<=log2(W)<=high<ceiling and an asserted
interval width<10^-50. Shifting the supplied ceiling one grid step down or
up must fail. This strengthens the certificate, not the attack assumptions.

## 41. Balanced giant/baby ratio Q

The direct parent 3c7fc300 / d966957b inherited the ratio Q = 29/25 = 1.16 from earlier drafts without optimization. Under the exact per-point costs BABY_PT = 1103.3875120541381 and GIANT_PT = 1279.5125120541381, the ratio GIANT_PT / BABY_PT is approximately 1.159622. Since total point work is proportional to N_b * BABY_PT + N_g * GIANT_PT subject to 6 * N_b * N_g >= KT * m, balancing baby and giant effort motivates candidate ratios near GIANT_PT / BABY_PT. Evaluating the discrete ledger across candidate rational ratios shows that Q = 115917/100000 yields total work of 153731988576905836562981836590626336476 instructions (saving 98159936534509748680470287 instructions over 29/25). No claim of unrestricted global discrete minimization is declared or required; the declared claim is strictly the certified upper bound W(2^20) <= 2^115.887900713539257 G1-reference units (C=2000).

The rational certificate log_bound.certify verifies ceiling 115887900713539257/1000000000000000 with width < 10^-50, and rejects both +-1 grid perturbations. All algorithm mechanisms, kernels (916/590), H1-H5 and H7 hypotheses remain identical to the direct parent.

