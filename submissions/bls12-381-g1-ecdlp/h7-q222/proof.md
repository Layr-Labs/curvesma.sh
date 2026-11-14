# BLS powers: run gate, round-stopped giants (expected work 115.055045677162760)

## Complete claim and source scope

The complete claim is claims[0] of submission.json: one cold run that draws a gate coin (K = 7236988362665215240; a closed gate is a charged failure), walks its giants in rounds and stops at the first verified s; success at least 0.39; W(n), the run's expected total work, at most 2^115.055045677162760 units (exact score 115055045677162760/10^15); every run at most 2^118.032068814159745 units; conditional on H1' (s44), H2-H5 of source-proof.md and H7 (s37). s2, s4, s43 and s44 derive every figure in it. The flat unsigned 64-bit node address and communication model is analytical, not hardware.

Source scope. Sections 1-42 and the code are the construction of our submission d1586a3f-8498-41f1-9c8d-7a64c236e552 (declared dependency), carried complete and re-executed, its record comparison repaired (s17, s25); no acceptance premise rests on it. Sections 43 and 44 are the change.

Parameters: r as protected, d=2^20 dividing r-1, m=(r-1)/d, six-element GLV/negation classes; Q, N_b, N_g: s2, s41. No previous hosted decision proves H1', H1-H5 or H7.

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

Pair and root prices are those of s26-s28 per child with the leaf kernel for leaf children, s30-s31 side prices and s32 kernels: current baby and giant point prices are s42's (pinned), no reference addition in the walk; root 486650135/262144 at s42's round share.

## 2. Opened-run success floor (s43, s44)

An opened run's floor is lo_run = 1-exp(-K)-0.0005-2^-20-7*2^-40-2^-108, where K = 6 N_b N_g/m >= K_T = 261/50 counts the giant target N_g and the subtracted terms are the H1' allowance at the full set (H1), the modelling budget, the H3/H5/EXT cap-tail union and the diagonal event (source-s7, s11-s12; s43). The run emits N_g giant points plus the rest of its last class walk (at most EMAX more), which only adds ratios: H1' at the full set gives coverage at least 1 - e^(-K_full) - 0.0005 with K_full >= K.

The floor is evaluated exactly: evidence.py encloses exp(-K/8) by the alternating Taylor series (40 terms, both partial-sum brackets) in rational arithmetic, takes eighth powers (bracket width below 10^-48) and rounds lo_run down to the 10^-60 grid: lo_run = 0.99409171719287617267... K_T = 261/50 and Q = 222/1000 are design values chosen by a ledger scan, not minima. N_g = 442687687325533329476781238004094501 is the least integer with 6 ceil(Q N_g) N_g >= K_T m (asserted both directions), N_b = ceil(Q N_g) = 98276666586268399143845434836908980, so K = 5.220000000000000000000000000000... >= K_T. Every N-dependent quantity is regenerated from N_b and N_g + EMAX by the same code with its asserted exponents (list: s43); tails are evaluated at these values, not inherited.

## 3. Explicit pool-cap charge and allowance 1

The one-instruction-per-point allowance covers residual bookkeeping; named pool-cap checks and key/index access costs are charged explicitly.

Residual bounds, retained from the base. (i) Pool-cap checks: charged exactly (wave row of s18), residual zero. (ii) Key reloads: zero (key stored to slot, reloaded by listed ldp+ldr, travels in-record). (iii) Head index: one lsr (s17), combine row deleted, residual zero. (iv) Cap enforcement: charged exactly (2 x cmp_cap x OWNERS, online only; counter in register, wave-end traffic in link row). (v) Pass control below 2^115.1: owner-init unrolled 2^20-fold (control below 2^99); all other bulk rows below 2^112 in sum (at s43's N the largest are gap-table entries BASES x GMAX = 2^109.62, planner scan visits 2^108.58 and augmentation ancestor steps 2^108.4, on top of the former 2^109 budget) with at most subs+b.ne control per level over four levels (below 8x sum, below 2^115); vs retained N_b+N_g (2^118.52). The allowance covers only purely-unknown misses.

## 4. Exact complete ledger

Unchanged source dimensions, quotas, capped schedules evaluated by evidence.py at (KT, F, ELL, KERNEL=916 (s40), KX, BETA = 128, BATCH, the s28 share, N_b, N_g; s38 tables). Every row below is an exact integer instruction count, rounding at original locations. The every-run total of an opened run (giants never stopped) is W_run = 679552726052473388381488596980646223864 instructions; the bound on its mean is E_run = 219994460429955022321640963315313057745 (s44); the 64-instruction gate is added in s43. An opened run's preprocessing, the challenge-independent baby phase, is 110802912844150860813505957232403002701 instructions, included in both. Baby transport 2180/64 per record, giant 3 x 1796/64 per point (s30); baby-phase walk part at conservative floor F_PRE=39/40 of s6, not F. Preprocessing takes the baby share QB = max(N_b/(N_b+N_g), BASES_B/BASES, (BASES_B+1)/(BASES+2)) of the per-point and per-base rows (padding, side calls, allowance, restarts, gap tables, strands, chain ids) and half of the wave rows (baby waves are under half of all waves, asserted).

- records: N_b + N_g points: 674772463584609614743531502632711104250.
- gap tables: 1589885355070353484986889774144032668.
- comparisons at the per-owner cap: 1374681265158185150885979014735855616.
- listing allowance (1/point): 540964353911801728764741860916859353.
- nonempty-lookup extracts and cap: 547748808775090825143073991961870336.
- comparison cap enforcement: 392766075759481471681708289924530176.
- owner initialisation: 166153499473114557899952177373249536.
- the other 29 rows, together 150842572767680052542680929107387887 (0.025% of W_run, each below 0.05%): names and values in the experiment output (ledger rows_instructions).

Memory and all per-node layout assertions are retained under the s8 slot size (generating-node term with 192-byte gap entries, owner term doubled by s11): total below 2^123.0005 bytes; the generating node and the owner (s43) below 2^64 bytes, every other node below 2^63 (s38). No free precomputation, future-target amortization, or wall-clock parallel discount.

Score grid is 10^-15: evidence.py reports the least fifteen-decimal ceiling 115055045677162760/10^15 of the exact log2 of E_num/(2000 x 2^64), E_num = 64 x 2^64 + K E_run (s43, s44; 90-digit Decimal enclosure). log_bound.py encloses ln(x) by 64 rational terms of 2*sum(t^(2j+1)/(2j+1)), t=(x-1)/(x+1), plus the positive geometric tail; similarly ln2; divides in the safe direction. Applied to E_num/(2000 x 2^64) it proves the declared ceiling and rejects both one-step shifts; applied to (W_run + 64)/2000 it proves the every-run ceiling 118.032068814159745.

## 5. Success and honest limitations

Points, stored keys, routing, candidate decisions, recovered scalars match the algorithm. No section s1-s37 changes a source success event: each changes prices, layouts, bit-identical listings or superset failure events in the same union, as each section states; current full kernels preserve exact x',y'; beta passes produce the live routed keys. Parity output is removed by s40. s30's x-keys with three giant lookups test exactly the U-orbit match the parity key tested (orbit lemma in s30); leaf kernels and the s30 planner change prices, not points. N changes only through s2. H1-H5 keep source-s7 definitions, evidence, budgets; H7 (s37) is a new disclosed hypothesis; H6 is superseded by the s30 per-class floors (H4 as f_c >= 998/1000 and H4L as leaf share >= 193/200 (s36) per production class, runtime-checked before any group op; a failed check stops the run, challenge-independent). The allowance change is exact accounting with stated residuals. Those hypotheses remain extrapolations requiring review.

An opened run's floor at the realized K (s2) is lo_run = 0.99409171719287617267...; with the gate (s43) success is at least (K/2^64) lo_run = 0.39000000000000000005... >= 0.39 (exact rationals; K - 1 fails, asserted); declared success is 0.39. All failures and capped work remain in the complete ledger.

Mode all of combined_lazy.py adds the planner controls, reduced-order curve recovery, transport with overflow controls and the round-order prefix census to the ledger-mode audits. Finite cases only; no full-size attack run.

## 6. Superseded pooled floor

H6 is withdrawn. Use the runtime-checked per-class floors of s30;
preprocessing uses F_PRE=39/40.

## 8. Production lane-stack layout (128-byte slots; 192 since s30)

Supersedes the source's 96-byte lane-stack stride. Slots are 192 bytes (s30; 128 at k7), contiguous from a 16-aligned base: x [0,48), y [48,96), keys [96,160), pad; every listed ldp/stp keeps form and offsets, slot addresses use one hoisted madd. The s40 full kernel writes [out,out+96), leaving the key region to beta_keys; gap-table entries are 192 bytes (x, y, v, n; s1), read-only to kernels; batch buffers are separate (144 bytes since s32). Memory term LANES x (2 DEPTH + 3) x 192. kernel_check.py's stride audit (s1): 8 out slots in one contiguous buffer under random canaries, outputs exact, write footprint exactly the out ranges plus the frame; the current 80-byte-stride negative control trips it.

## 9. Record-path re-derivation: pointer buffers, honest matching, explicit setup

This section re-derives every per-record instruction from scratch, superseding the source s5.2 transport counts. The charged record sequences (GEN12 s10, RELAY s9, INS s22, LOOK/EXT3 s23, CMP7 s25, loops s29) are asserted by evidence.py audit_record_path (lengths, bits, immediates, pool fit, division proofs); the proto toy executes the pointer mechanism on every record buffer (pre-write check, +40 advance, exact derivation, reset).

Pointer buffers. Sender per-link state becomes (ptr, end): each record's three stores post-index ptr by 16+16+8 = 40 bytes, deleting slot-address madd and counter increment. Entries stay 16 bytes (one ldp; source "24-byte" prose counted a hoisted cap word), addressed by one add-with-shift (t/u below 2^30, lsl #4). The cap check precedes the write (cmp ptr, end; b.hs fail): identical failure event, never an out-of-bounds write. GEN listing: s10. The relay output side matches GEN's; relay input loads unchanged; s24 pairs the relay loop. At wave end the sender derives count = (ptr-base)/40 (sub, lsr #3, umull/lsr /5 with magic 0xCCCCCCCD, exact by m*5 = 2^34+1 plus asserted bound, count store, ptr reset: 8), writes the per-parity count word, and resets; receiver side (3) and double-buffer discipline unchanged. LINK_WAVE is 16 (8 sender + 3 receiver + 5 margin); the per-node term rises to NODES*24 for dmb, counter traffic and constants (16 used, 18 with s12 EXT traffic). Chain-id words are segment-invariant, formed once per segment by the formation row (register-holding superseded by s10; GEN reloads them).

Honest matching: this section's record key, pool layout and INS/LOOK/compare/cap counts are superseded by s11, s12, s17, s19 and s22-s25, s29 (current bodies in evidence.py, executed by audit_ins12/audit_look6). Four new explicit rows: window-table init (8/link, 4 used), chain-id formation (3/segment), candidate-link wave counts (9 on 4 traversals: /40 magic, /64 shifts, /56 with proven /7 magic), candidate formation (16/match, 11 used). Caps, routing, counts, events unchanged, so H1-H6, N, success, memory untouched.

## 10. k9 repair: chain-word memory schedule

Two chain-id words formed once per segment into a per-lane 16-byte slot (LANES slots per generating node, in memory term); every GEN reloads them. Formation 4 per segment (ubfx+bfi+lsr plus stp); GEN11 becomes GEN12 (ubfx t, entry address, ldp entry, ldp chain words, cmp, b.hs, three stores, ptr store-back, key-reload/loop pair). No chain word crosses records in a register (except within one giant point's three consecutive GEN copies, which no call separates; s34), so lane interleaves and side-output clobbering are compatible; allowance untouched.

## 11. Comparison cap per owner

This section sets the per-owner comparison cap. The head array has 2^59 slots per owner (key layout, extracts and compare: s17, s25), superseding the head-slot/H3 sentences of source-s4.3/7. cmp_cap is 170335133046814413 per owner at the s2 sizes, i.e. 1.00001 x 3 N_g N_b/(OWNERS^2 x HEADS) since s30 (3 N_g giant lookups; was 1.01x). Support (exponents asserted): an owner's comparisons sum the chain lengths of its looked-up slots. Under H3 distinct giant points' keys are independent uniform and independent of the baby layout; a point's three keys are fixed bijections of x', so its k_g <= 3 lookups at the owner form one block in [0, 64 k_g], independent across points given the layout. Chains exceed 64 with probability below 2^-200 (Poisson tail at slot load N_b/2^119 < 0.15, below 2^-324 over the 2^119 slots; in the failure union). k_g >= 2 has probability 3 x 2^-120 per point (mean below 0.86 per owner, s43), so sum k_g^2 <= G_o + 1800 except with probability (e 0.86/300)^300 < 2^-2000. Hoeffding over blocks at relative deviation 9.7x10^-6: exp(-2 (9.7x10^-6 E_o)^2 / ((G_o + 1800) 64^2)), exponent above 94 nats/owner (2^-75 over 2^60 owners). B_o is Binomial(N_b, 2^-60); G_o = sum k_g is Binomial(N_g, ~3 x 2^-60) plus at most 600; both means exceed 5.5x10^16 and are exceeded by 10^-7 relatively with Chernoff exponent above 180 nats; 2 x 10^-7 + 9.7 x 10^-6 < 10^-5. The same blocks cover s12's lookups. All three sit far inside the retained 2^-40 budget. Owner init: s34.

## 12. k11 conditional extracts at an occupancy cap

This section charges the key extracts only on nonempty fall-throughs (EXT, listings since s23), at a per-owner occupancy cap. The counter lives in a register across the owner's lookups with wave-end traffic in the link row (s9 pattern; link-row use 16 to 18 of 24); exhaustion aborts into the tail budget: cap-plus-tail, never expectation. Per-owner cap is ceil(1.000031 x mu_L x p) = 158365452920080162 (at the s2 sizes; was 1.02) with mu_L = 3 N_g/2^60 exact (s30), p = y-y^2/2+y^3/6 at y = N_b/(OWNERSxHEADS-1): under H3 each key is uniform over the OC cells, so P(nonempty) = 1-(1-1/OC)^N <= 1-e^-y <= p (ln(1-x) >= -x/(1-x); decreasing alternating Taylor, asserted), and E[nonempty] <= mu_L x p by the tower argument over independent placements. Three 1.00001 stages cover the cap (1.00001^3 = 1.0000300003 < 1.000031, asserted): lookups/owner (Binomial, Chernoff), occupied/owner (negatively-associated, Chernoff, Dubhashi-Panconesi), probes-given-both (Binomial, stochastic dominance). Per-owner means exceed 1e16 (asserted: mu_L, occupied floor HEADS x 0.17); at delta = 10^-5 every Chernoff exponent mu delta^2/3 exceeds 3x10^5 nats, far above the asserted 70.5-nat budget (2^-40 over 2^60 owners), so the success union stays 7x2^-40. The toy counts nonempty probes per lookup phase. EXT row 3 x OWNERS x cap (EXT3, s23); capped EXT3 beats 1/lookup (asserted).

## 16. Augmented planner (pointer; floors and planner since s30)

The s16 augmentation (1/3/5 alternating paths, partner cap 32, length-5 cap 8, exact rollback, OPS_BOUND = 10752 per element, asserted) runs unchanged inside the s30 planner, over centre pairs and in the fallback; its floor and tables are superseded by s30.

## 15. Exact restart price SMX from the shipped kernel bodies

Supersedes the s13 price; only the restarts row changes. Each restart computes [b]Q_a, b in [1, X2] (<= 77 bits, from the plan entry), by branchy Jacobian double-and-add (dbl-2007-bl, 2M+5S; madd-2007-bl with affine Q_a, 7M+4S) on an all-ones worst path in Montgomery residues. Strands, base points and candidate tests stay at SM = 8192000. SMX = 48 + 380 + 5 + 76 x (2340 + 11 + 3656) + 6 + 172520 + 1231 = 630722, an upper bound on every b; row count x (SMX + 16). Checked by smx_audit.py.

## 17. Record key layout: 59-bit head, 9-bit field

The 192-bit record key (words w0, w1, w2) is split as owner bits 0..59, field f9 = bits 60..68, head bits 69..127, k64 = bits 128..191 (60+9+59+64 = 192, asserted). The head index is w1 >> 5, one listed lsr (HEADS = 2^59 = 2^(128-69), asserted). Owner routing (bits 0..59) and the H3 statistics are those of s11 and s12 at 2^59 heads.

INS (body: s22): extr x10,x1,x0,#60 puts key bits 60..123 in x10; bfi x4,x10,#55,#9 deposits f9 into w4[55..63]. The chain word w4 holds the 35-bit segment in [0..34] (s42) and the pool link is its own word, so bits 55..63 are free (35 <= 55, asserted); bits 35..54 are never read.

Compare (routine: s25): and x14,x14,#0xFF80000000000000 keeps w4[55..63] of the stored record; cmp x14,x10,lsl #55 compares it with the giant's x10 shifted left by 55. The shift discards x10's bits 9..63, so exactly key bits 60..68 land on bits 55..63 and no head bit reaches the comparison: the two operands are equal exactly when the two f9 fields are equal, for every value of the head bits. ccmp then compares k64. With the owner (routing) and the head (chain) already equal, a comparison succeeds exactly when the two 192-bit keys are equal.

Audit (every ledger run): audit_packed_match evaluates this displayed compare word for word, with no mask on the giant operand, (w4 & 0xFF80000000000000) == (x10 << 55) mod 2^64: the key 2^69 (head index 1, f9 = 0, x10 = 512) and 3013 further keys with random junk in w4[0..54] compare equal to themselves; any single flip in bits 60..68 compares unequal; a flip elsewhere leaves the field compare equal; the former constants (field and mask at 54..62, shift 54) reject the key 2^69. The executed routine: s25.

## 18. Wave-granular pool-cap check

This section supersedes the per-insert pool-counter compare (N x 2, deleted) with one check per owner per wave. An owner inserts only through its 2^20 inbound relay links (source-s4.3), each carrying at most CAP2 = 2^21 records per wave, enforced by the relay's already-charged pre-write cap check (GEN/RELAY cmp + b.hs; overflow is a run failure in the same budget). So an owner inserts at most 2^41 records per wave (asserted); checking counter <= 2^57 - 2^41 (pool 2^57 since s43) at each wave start guarantees no wave overflows the pool (start bound + one wave = cap, asserted). The failure event (count above threshold at some wave start) supersedes the old overflow event, taking the same H3-tail union member: union and success formula stand. Tail: mean N_b/OWNERS exact, excess d in (0,1), Chernoff mu d^2/3 above the 70.5-nat budget (above 13.3x10^15 nats at s43's N_b; asserted).

Schedule, all counted (counter = pool pointer since s22): counter and threshold live in pool slot 0 (allocated, unused by records; zero memory delta). Per owner per wave: ldr, ldr, subs, b.hi at wave start plus str at wave end (5; no cross-wave register lifetime). Once per owner: movz counter, threshold in movz + 3xmovk (asserted), stp (6). Row = OWNERS x (5 x WAVES + 6) = 79761387861594776929304576 (WAVES = 13836394, s42/s43). Both phases charged though giant waves never insert (conservative); preprocessing takes the row in full; wave-loop control (~2^82) in the (v) bound.

## 19. Full-address chain links

This section restates the link word (superseding the s7 56-bit-link sentence): head entries and chain links store full 64-bit record addresses, 0 = end or empty; records load as [addr,#imm] (comparison CMP7 since s25, lookup/EXT3 since s23, head store in the s22 loop). Layout: pool after the head array at nonzero base HEADS x 8, records at base+32 and above, addresses below 2^63 (asserted). Pool use: s42; word sizes, memory, owner init unchanged. Links never leave the owner (chain ids are base+segment).

## 20. Itemized side-call fixed price: 40 per invocation

This section replaces the s1 64/invocation bound with a complete per-call schedule: 32 fixed (prologue 22 + epilogue 10 with ret, asserted every call) + 1 exit branch (measured) + 5 staging (mov x0,x_args; mov x1,x_n; lsl x2,x_n,#5; stp xzr,xzr,[x_args,x2]; bl) + 2 round control (subs, b.ne) = 40 (asserted). The per-lane argument records sit in PAIR_BK/ROOT_BK (s26), so staging is exactly these 5. Zero spills: 8 loop variables in x19-x28 across the call; the kernel preserves callee-saved regs (asserted); nothing live in caller-saved regs across bl. Lanes are unrolled x n (n fixed; PADOPS pads to full calls), so no per-lane loop control exists; the 40 is n-independent anyway. Driver and kernel share one small text (bl in range, no veneer). Per-base-per-phase setup bounded by 32 (entry/exit 10 + materialize 8 + base-loop 2 + 12 margin; true cost <= 20), charged as 64 x (BASES+2). Row: s1 (eight calls per round since s30; 48 per call since s37), both phases charged.

## 21. Buffer-granular EXT-cap check

The per-lookup b.lo moves to a buffer check (EXT3 since s23). A fall-through lookup runs subs + extr only; before each link buffer the owner checks remaining >= 2^21 (cmp-immediate + b.lo: 2). A buffer holds at most CAP2 = 2^21 records by the relay's charged pre-write check and EXT runs at most once per lookup, so a buffer adds at most 2^21 EXT (hard) and the total never exceeds cap_ne (asserted). The failure event (above cap_ne - 2^21 at some buffer start) supersedes the old one and takes the same s12 union member: cap_ne - 2^21 still covers the three 1.00001 stages (asserted at s43's sizes), so tails and success stand. Schedule, all counted: EXT counter in pool slot 0 after the s18 pair (24/32 bytes, asserted; zero memory delta); per owner per wave ldr + str (2); per buffer cmp + b.lo (2, in the buffer loop, control in link row); one-time movz + 3xmovk + str (5, asserted). Row = OWNERS x (WAVES x (2 + 2 x 2^20) + 5) = 33454379818480017084081720262656 (WAVES = 13836394); all waves charged, giant only needed (conservative); preprocessing takes it in full. Toy enforces per-link pre-write caps; thresholds are ledger asserts.

## 22. Pool pointer: insert pair loop

This section supersedes the s17/s19 insert. The owner's pool counter becomes a pointer p to the next free 32-byte record, deleting the slot-address add and the counter add (the s9 sender device at the owner). The insert body in order: ldp x0,x1,[x20],#16; ldp x2,x3,[x20],#16; ldr x4,[x20],#8; lsr x9,x1,#5 (head, s17); ldr x11,[x21,x9,lsl #3] (old head); extr x10,x1,x0,#60; bfi x4,x10,#55,#9; stp x3,x4,[x22,#16]; str x22,[x21,x9,lsl #3] (head = record address); stp x2,x11,[x22],#32 (k64, link = old head; the only write of p, after both stores through it). s24 runs 64 bodies per iteration sharing one subs + b.ne (s29: 10 + 2/64 per record); setup and tail ladder add at most 10 per invocation. Record words, chains, links, comparisons and events are unchanged. s18 on the pointer: pool base PB = HEADS x 8 = 2^62 (s42 heads), first record PB + 32, threshold TP = PB + 32 x THRESH = 2^63 - 2^46; p <= TP at a wave start plus one wave (at most 2^41 records) ends at PB + 32 x POOL = 2^63, the pool end (asserted; s43); the wave check compares p with TP (same five instructions, same failure event); p and TP each initialise by movz + movk inside the charged 6 (asserted). audit_ins12 executes the s29 loop with a64.py on 600/601/1/0/639/38 records over 40 heads: 10n + 2(n div 64) + 10 steps, final p, and memory exactly the oracle's (records, links, heads; nothing else written).

## 23. Deferred k64: lookup pair loop + EXT3

This section supersedes the s12/s21 lookup and EXT. The lookup body is ldp x0,x1,[x20],#40 (w0,w1; ptr past the record), lsr, ldr head, cbz (empties skip free, ptr already advanced); fall-through runs EXT3 (subs, extr x10, ldur x2,[x20,#-24] = w2 of the passed record). s24 runs 64 bodies per iteration sharing one subs + b.ne (s29: 4 + 2/64 per lookup plus EXT); setup and tail ladder add at most 10 per invocation. LDP #40 / LDUR #-24 are one instruction each (ranges asserted). w2 feeds only CMP (s17 schedule); buffer, heads, events unchanged; s21's buffer cap still bounds EXT. audit_look6 runs the s29 loop in a64.py on 600/601/1/0/639/38 records/40 heads: 4N+2(N div 64)+3K+10 steps, final ptr, x24, single-iteration regs, untouched mem; wrong advance, offset or lsr and dropped ldur each fail it.

## 24. Transport loop unrolling (64-fold since s29)

This section supersedes the loop control of the s9/s10 RELAY14 listing and the s22/s23 single-record loops (bodies unchanged). The loop pair (subs + b.ne) is each transport loop's only cross-record instruction, so every loop runs UNT = 64 body copies back to back sharing one pair (s29). Each copy reuses the previous copy's registers after its stores complete, in order; pair-control scratch x5 is dead in the exhibited bodies (asserted absent) and exists for RELAY by the s9 dataflow (5 input words, routing, entry, ptr/end: at most 13 live of 31). Setup is and x5,x23,#63 + lsr x23,x23,#6 + cbz; the tail is a binary ladder: for j = 5..0, tbz x5,#j skips 2^j body copies. Per record RELAY 12 + 2/64 (body asserted verbatim vs the s9 mnemonics), INS 10 + 2/64, LOOK 4 + 2/64 plus unchanged EXT3; a tail record costs its body only. Per invocation at most setup 3 + six tests + ret = 10 (UNROLL_FIX, asserted), incl. empty loops. GEN12 has no per-record loop pair (driver control in s3(v), unchanged); the s25 chain loop alternates two copies and is not unrolled further; EXT3's subs is a per-item counter; lanes and owner-init are already unrolled (s20, s3(v)). Charged: 770/64, 642/64, 258/64 per record (average 497/16, asserted) plus the unroll row 10 x ((LINKS1+LINKS2) + 2 x OWNERS) x WAVES (s4; both phases charged, one used). No new heuristic or premise change: H1-H5, N, success, memory untouched.

## 25. The chain comparison, as an executed routine

A nonempty lookup (s23) holds x10 (extr, key bits 60..123), x2 = k64 and x11 = the head's record address, and walks the chain inline before the next lookup record. A pool record is k64 at +0, link at +8, chain words at +16 and +24 (s22). The routine is two copies that alternate the record register (x11, x13), so no move is needed between records; copy A (evidence.py CMP_A; copy B swaps x11 and x13 and ends cbnz x11, A):

    subs x25, x25, #1          ; s11 comparison budget (CAP 2)
    b.lo fail                  ; budget exhausted: run failure
    ldp x12, x13, [x11]        ; k64, link
    ldr x14, [x11, #24]        ; w4
    and x14, x14, #0xFF80000000000000
    cmp x14, x10, lsl #55      ; f9 (s17)
    ccmp x12, x2, #0, eq       ; k64
    b.ne next
    str x11, [x26], #8         ; match: record address to the candidate list (candidate formation row, s9)
    next: cbz x13, done        ; end of chain: on to the next lookup record

Per comparison: CMP 7 (ldp, ldr, and, cmp, ccmp, b.ne, cbz or cbnz) plus the 2-instruction cap pair, the rows OWNERS x cmp_cap x 7 = 1374681265158185150885979014735855616 and OWNERS x cmp_cap x 2 (s11 cap 170335133046814413; asserted). Pool base 2^62 is 16-aligned and records are 32 bytes, so every ldp is aligned. No signed condition; addresses stay below 2^63.

Executed (audit_cmp, every ledger run, a64.py): the INS loop of s22 inserts 127 baby records, among them a twin of the key 2^69, equal heads with different f9 and equal f9 with different k64; the lookup loop of s23 with this routine inline then processes 315 giant records: the key 2^69, exact copies of baby keys, and for each copy one flip in bits 60..68, in k64, in the head and in the owner bits. Against the full-key oracle: the candidate list equals, in order, the records whose key bits 60..191 equal the giant's on its chain (113 matches, 54 with key bit 69 set); exactly 4967 comparisons; step count 4N + 2(N div 64) + 3 x nonempty + 9 x comparisons + matches + 10 = 47021; nothing else written; one unit short of budget takes the fail exit. The same audit with the field at bit 54, and with mask 0x7FC0000000000000 and shift 54, must fail, and does. The transport simulation (proto) uses the same field compare and replays a quarter of its giant keys from baby keys.

## 26. Batch round: both children by the derived kernels, no reference addition in the walk

Points, keys, records and events are unchanged (canonical residues are unique): H1-H5, N, success stand.

Round (one step of the 2^18 lanes): the driver writes argsB[i] = {out, P, T, buf+144i} per lane op (out = the P+T child, or the root's point) and argsA = {P-T out, P, T, buf+144i} per pair (since s30 each split into a full and a leaf array by the leaf flags), each ended by a 32-byte NULL record; padded lanes use the fixed pair (Q_a, T[2]), x distinct, so every denominator is nonzero (n = 2^18, a multiple of the s28 unroll factor 32). batch_round(argsB, buf) (batch_audit.py, generated from the shipped MULNR body and s15's C2/Fermat): forward d_i = x1 + v2 in (0,2p) (v2 = p - x2 at T+96; s27) at buf_i+48 and c_i = S(c_(i-1), d_i) at buf_i (c_(-1) = R mod p at buf-144); the canonical c_(n-1) is inverted by s15's Fermat and negated once (t = p - t); backward inv_i = S(t, c_(i-1)) at buf_i, unreduced (s27), then t = S(t, d_i). side_lanes(argsA) then writes P-T and side_lanes_plus(argsB) writes P+T: side_lanes with its y2 loads (T+48/64/80) and n2 loads (T+144/160/176) exchanged (9 ldp offsets), i.e. side_lanes on -T = (x2, p-y2, v2, y2), so every kernel bound holds verbatim. With b < 2p and p/R < 0.1016 every MULNR operand a stays below 1.26p (K1 holds) and inv_i < 1.17p, congruent to 1/(x2 - x1) (s27).

Counts (current: s28 loops, s32 body): per element forward FB + backward 595, one cbnz and one cmp + b.ne per 32 elements, BATCH = (32 FB + 1 + 32 x 595 + 2)/32 per lane op, the only path. Round fixed 172585 (with the negation) plus 5 per call (two moves, bl, terminator address, stp) = 172590 per round of exactly 2^18 lane ops (s20), charged as the share s = 172590/2^18 per lane op (s28). PAIR = 2 KERNEL + BATCH + 28 + s for two points (plan entry 8, P address 1, T[g] address 6, buffer 1, two out addresses 2, two argument records 4, s30 routing6); ROOT = KERNEL + BATCH + 21 + s (gap entry 5, T[gap] address 6, P and out addresses 2, buffer 1, argument record 2, loop 5). P/T and inverse loads and coordinate stores are inside the executed listings. Side-call row: s1 (48 per call since s37). Re-walks keep 3300 + 232864.

Audits: every ledger evaluation executes one 64-lane round of genuine points with both side calls (batch_audit.run) and kernel_check's 509-lane suites; results are in the experiment output.

The reference addition (2000: 3M formula plus 3M batch share, A = 0) is no longer used anywhere in the walk. P+T, the batch products (three per lane op, 3n >= 3(n-1)) and the inversion are other work charged at their listed and executed counts, the rule s1 already applies to P-T; nothing is charged twice or omitted.

## 27. Lazy batch round: reflected difference and unreduced inverses

Changes only the batch listing of s26 (measure/side_lanes.S unchanged, SHA-256 b4c7a7c2... as in measurement.log). Prices only: points, keys, records, events, H1-H5, N, success and memory stand.

Reflected difference. The forward pass forms d_i = x1 + v2 from v2 = p - x2 at T+96, present in every entry the round reads (gap entries, base points, the T[2] pads; s1), by six adds; d_i is in (0, 2p) and d_i = -(x2 - x1) mod p, nonzero as x1 != x2. Montgomery's trick is linear in its running value t, so the round negates the canonical Fermat output once (t = p - t in (0, p), six instructions in the round-fixed part); every backward t and inv_i = S(t, c_(i-1)) is then the negative of its value for d_i, so inv_i = -1/d_i = 1/(x2 - x1) mod p, the inverse both kernel calls need.

Unreduced inverses. The C2 after inv_i is deleted. With t, c < 1.26p (s26), inv_i < 1.26^2 x 0.1016p + p < 1.17p. The kernel reads inv only as the b operand of its first product, and K1 (s1) holds for every 384-bit b. Its bounds become: t < 2p; mu < 2.34 x 0.1016p + p < 1.24p < 2^382; V < 1.24^2 x 0.1016p + p < 1.16p; z < 2 x 1.24 x 0.1016p + p < 1.26p; X in (-p, 2.16p) and D in (-p, 1.26p) stay inside the existing reduction arms (X < 3p, D < 2p), so x',y' stay canonical; s40 removes unused parity output. The square needs 2c < 2^704 for its cross sum c: mu < 2^382 puts the top word below 2^62, so c < 2^702 + 2^639. The comments in measure/side_lanes.S state the canonical-inverse bounds (inv fully reduced, mu < 1.21p < 2^381, V < 1.15p, z < 1.25p); this section supersedes them, and no instruction changes.

## 28. Exact round share and unrolled batch loops

Changes the s26 round share and the s26/s27 batch loop control only. Exact share: a walk round is one lockstep step of a node's 2^18 lanes, each doing one pair or root operation, idle lanes padded and charged (source s2 step 4, s3 padding lemma; s20): exactly 2^18 lane ops and one fixed cost 172585 + 5 = 172590 (s26), so each lane op carries s = 172590/2^18, not 1. Every executed lane op sits in a row charged at PAIR or ROOT (records at the floor F, which charges at least the lane ops of any actual f >= F; padding; remainders), so the shares cover every round. Table strands keep their rounded share.

Unrolled loops: the forward loop's only per-element control is cbnz x0, the backward loop's cmp x2,x30 + b.ne. Each loop holds 32 copies of its unchanged body and one control; bodies address only through x30 (pre-indexed) and x2 (post-indexed); n = 2^18 is a multiple of 32, no tail (asserted). Counts and listing size: s32 (inside gen_mem's 2^20-byte code term). Evidence: batch_audit.run executes one 64-lane round, each loop branch both ways.

## 29. Transport loops 64-fold with a binary tail ladder

The s24 loops run 64 bodies per pair with a six-step binary tail ladder; evidence.py generates them (unroll). Average record transport 497/16, baby phase 2180/64 per record; candidate delivery is eight relay hops at 770/64; unroll setup and tails are 10 per invocation. Evidence: audit_ins12 and audit_look6 (s22, s23).

## 30. x-only leaves, beta-triple giant lookups, leaf planner and asymmetric sides

Supersedes the record key (source s1 parity key), source s3 matching, the s16 floor, GMAX 2^44 -> 2^47 and DEPTH 2^14 -> 2^17 (source s2-s3), N of s2, the 128-byte slot of s8, the per-point price of s26-s29, the two side calls per round of s20/s26 and the 2^34-step giant wave (source s4.3). Everything else stands.

Keys. Every record key is limbs 1..3 of the canonical Montgomery residue of x' (out+8..31). [-1](x, y) = (x, -y) and [lambda](x, y) = (beta x, y), beta fixed by [lambda]P0 = (beta x0, y0) (asserted in beta_run), so the U-orbit of Q has x-set {x, beta x, beta^2 x}; conversely x_P in that set gives y_P^2 = x_P^3 + 4 = y_Q^2, so P is in U Q. A baby x-key equal to one of a giant's three keys is exactly the orbit match the parity key tested, and at most one of the three can match (x = beta x forces x = 0, a 3-torsion point; 3 does not divide r). Babies insert one record per point; giants send three lookup records per point (keys of x', beta x', beta^2 x', one chain id; GEN unrolled three times per point at key offsets +8/+104/+136; a giant re-walk keeps positions where any of its three keys equals the candidate key). Coverage events, K, success and the 96 recovery tests are unchanged. Giant waves are a quarter of a pair wave (2^31 steps since s42), so a generating-node link still carries at most 2^23 records per wave in the mean (six per step).

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

Sides. With side prices c_b, c_g the optimum is N_b/N_g = c_g/c_b (1.169 at s31's 23/20). The ec toy (x-keys, triple lookups, solved == exact coverage every trial, all scalars correct): the declared mode-all input (K = 4.5, Q = 13/50) gives 8/8 at 0.9889 predicted. Its rate is a binomial sample of toy coverage: it checks that success is exactly the coverage event; production success rests on H1 and the census (s2, source s7), not on toy rates. Caps with 3 N_g lookups and N_b inserts: s43. H3 now covers N_b + 3 N_g keys: a giant's three keys are images of x' under a fixed bijection, treated as three uniform keys; the declared toy's comparisons (259.6 vs 248.2 model) and candidates (5.75 per trial: 1.939 false-match model plus true hits) follow the three-lookup model; the re-walk step price 3300 covers 2 x 929 + 2 x 343 and key compares.

Prices (pinned in audit_record_path; current values s1): pair child (BATCH + PAIR_BK + 6)/2 + ELL x KX + (1 - ELL) x KERNEL; baby point F x child + (1 - F) x ROOT + 2180/64; giant point F x child + (1 - F) x ROOT + 3 x 1796/64 + BETA. Records = N_b x baby + N_g x giant.

## 32. Shared x-sum e = v2 - x1, once per lane op

Prices only. Both children of a lane op formed X = V + v2 - x1 from the same P, T (reloading x1, v2; 12 adds/subs). The forward body, which holds x1, v2 for d, now forms e = v2 - x1 (subs, 4 sbcs, sbc into the dead x0, x1) and stores it at buf+96 (three interleaved stp): FB 324, BATCH 29411/32. Slot 96 -> 144 bytes (immediates only: strides 144, c_(i-1) at x2-224; memory (LANES+1) x 272; listing 201996 instructions, inside the 2^20-byte code term).

Kernel: kernel_check.esum_text replaces one asserted block of the measured listing by `ldp x29,x0,[x30,#16]` (T, buf), three ldp of e, the v2 loads, X = V + e (adds, 5 adcs), `tbz x16,#63,.Lx_nonneg`. As integers e is in (-p, p), so X is the old integer, in (-p, 2.16p), inside (-2^382, 2^383): bit 383 of its two's complement is set exactly when X < 0, the old b.hs condition; x', y' are bit-identical (x19..x24 now hold e, rewritten before any read). KERNEL 923 (24 paths); the s30 leaf cut also drops the v2 loads: KX 587 (4 paths).

## 33. Gap tables by pair levels

Changes only how tables T[k] = [k]Q_a, k <= GMAX, are built (entries, layout, memory and use unchanged). With T[1..M] known (M = 3 from the strand-start scalar multiplications, row unchanged), one level computes the centre c = M + k + 1 as a root T[M] + T[k+1] and then the pairs T[c -+ j] = T[c] -+ T[j], j = 1..k, k <= M - 2: the s26 lane op on (T[c], T[j]) with both children full (916 kernel, minus and plus calls, then the 21-instruction v/n step), so T[1..M + 2k + 1] is known. Full levels (k = M - 2) give M' = 3M - 3; the last stops at GMAX: 30 levels. Incomplete additions are safe: c - j, c + j and c are distinct, nonzero and below GMAX < r, and x(T[c]) != x(T[j]) as c != +-j.

Price per base: every level is charged ceil(k/LANES) + 1 full rounds of LANES lane ops at PAIR (+1 for the centre; padding included) plus 128 per round for calls, and 21 per entry: BASES x (rounds x (LANES x PAIR + 128) + 21 GMAX) (strands: BASES x (GMAX - 1) x 2050; re-walks keep 2050). Row in s4.

Evidence: batch_audit.run(table) executes a 64-lane level (M = 66, centre [131]Q, gaps [j]Q): the s32 round and both kernel calls give exactly [131 -+ j]Q; evidence.py asserts the 30-level recursion.

## 34. Giant GEN chain-word reuse; head zeroing by dc zva

Prices only. (i) A giant point's three lookup records share one chain id (s30; GEN unrolled three times per point, key offsets +8/+104/+136). Copies 2 and 3 drop their chain-word ldp (copy 1's registers still hold the words; GEN only reads them). Giant transport 3 x 1796/64 - 2; the proto's generic per-record GEN 12 stays an upper bound. (ii) Owner initialisation zeroes each 64-byte-aligned head array by dc zva (64-byte blocks, DCZID_EL0.BS = 4, common AArch64) plus one add: 2 instructions per 64 bytes, not one stp per 16. Row OWNERS x (HEADS x 8/64 x 2 + 64).

## 35. EMIN = 2^56 - 2^52 (H2 restated)

Bases are charged at N/EMIN per side; EMIN was 2^55, half the mean X2/d = 2^56. H2 now has class lists between EMIN = 2^56 - 2^52 = 0.9375 x 2^56 and EMAX = 2^57, still sieve-checked before any chain starts (a short class stops the run: an H2 failure event). Bases at the bound: less than 2^62.51 per side, 2^62.80 in all (I's 2^64.28 primes: margin above 2.78; s43); segments per class at least 2^31.9. Only per-base rows change. Support: class-length spread shrinks like (X/d)^(-1/2). Executed (evidence.py h2_sizes, mode all) on toy primes, [d, log2(X/d), min, max of |L_c| d/X]: [16, 10, 0.9648, 1.0391], [16, 14, 0.9908, 1.012], [64, 12, 0.9714, 1.0295], [256, 10, 0.9287, 1.1035], [256, 12, 0.9646, 1.0317] (offline, d = 16 at 2^18: min above 0.994). At 2^56 per class (d = 2^20) this extrapolates below 10^-6, far inside 6.25% (heuristic).

## 36. Leaf floor ELL = 193/200

H4L's runtime-checked floor rises from 96/100 to 193/200 (a class below it stops the run before any group operation, as before). Per-class minima: s30's table (0.97213-0.97349); new offline planner.py rows (kappa 8, four classes): d = 8 span 4: 0.97274; d = 8 span 8: 0.9728; d = 16 span 4: 0.97289; d = 4 span 8: 0.97524. One smaller two-class d = 8 run (lists of 7625 and 8734) had a class at 0.9691; 0.965 is below every measured class and ell rises with span (production 2^28). Prices only.

## 37. Leaf reduction at a capped share; H7

The leaf's X in (-p, 2.16p) needs the second subtraction (+7; paths 580/580/587/587) only if X >= p; the leaf counts those lanes (ldr/add/str: 590) in its free slot [sp,#96], zeroed in the prologue, returned in x0 (fixed 34). The driver adds each return to a node count C and its lanes to L (inside the 48 now charged per call, 40 + 8). At every wave barrier the tree sums (C, L) (16 per node per wave) and the run stops unless C <= 7L/32. Between barriers C grows by at most one wave of leaf lanes, so a leaf child costs at most 580 + 10 x 7/32 = 582.1875 plus the row GENS x WSTEP x LANES x 2 x 10 + WAVES x GENS x 16.

**H7 (new, in the claim's conditioning).** At every wave barrier, the cumulative share of leaf lanes with X >= p is at most 7/32. A violating run stops (excluded by H7, as H2's sieve checks). Support: exact arithmetic on 2 x 10^5 random lanes: 0.1878. batch_audit.xge_run (declared; one base Q, gap table [j]Q, walked centres, exact s32 batch-round inverses, 16 rounds of 256 lanes, both children) gives 0.1874, per-round shares 0.152-0.232 at 512 flags. Its model equals the executed listing's inverses (asserted in batch_audit.run); a barrier share averages >= 2^99 lanes. run_leaf (plain and plus): the returned count equals the exact X >= p count (61 of 365), body exactly 580 n + 10 count.

## 38. Seven-chunk beta tables (k37l)

Replaces s31's tables, lookup and evidence. Canonical Montgomery X < p < 2^381 splits at bits 0, 55, 110, 165, 219, 273, 327 into three 55-bit and four 54-bit chunks c_i; T_i[c] = beta c 2^shift_i mod p, canonical. Tables 0-2 interleave at A + 256c + 64i, tables 3-6 at B + 256c + 64(i - 3), B = A + 2^63 (one add-lsl per address); q p (q < 8) at A - 512 + 64q; production A = 512. BUILD is s31's 24-instruction recurrence at stride 256: row GENS x (5 x 2^55 x 24 + 7 x 2000 + 64) in preprocessing; storage 3 x 2^62 + 512 bytes.

Lookup, 128 per giant point: 3 loads, 12 chunk, 7 address, 21 table loads, 36 sum, 23 reduction, 2 stores, 22 for beta^2 X, 2 loop; fixed call part 26 of the 48 charged (s37). S < 7p and q = floor(s5/(p5 + 1)) <= 6 give q p <= S < (q + 1)(p + 2^320) <= q p + 2p (7 x 2^320 <= p): one conditional subtraction; keys and offsets unchanged.

Addresses: each node is one flat 64-bit space (source s4.1). Tables fill [0, 3 x 2^62 + 512), the other generating-node regions (windows included) follow, total 17609118817070351136 < 2^64 (asserted); owners stay < 2^64 (s43), other nodes < 2^63. Addresses are unsigned and never wrap: the ledger asserts that no shipped listing or mnemonic list has a signed condition, sign-extending load or asr, and that bit 63 is tested only on data (the kernel's X sign, the restart scan's scalar words); caps compare with hs; rel32 offsets add modulo 2^64 within one region. Analytical, not hardware.

Evidence: beta_run and batch_audit.run (finite checks; experiment output). Current prices: s4/s40.

## 39. Inherited ratio and evidence repairs

The ratio is Q of s41. H1 is asserted at these exact sets and realized K. Diagonals and WAVES: s43.

Planner depth and height use32-bit fields for DEPTH=2^17:28 bytes base,
8 adjacency,4 previous link,3 centre, totaling43 bytes/element. Wider
loads/stores retain the instruction prices. The character check certifies order d
for7^m and exhibits k=256 where log_g(x^m) is the wrong character base;
source section2 uses the correct log_(g^m).

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

The raw measured assembly
and its log are unchanged historical inputs;916 is an analytically counted
derived listing, not a new native measurement.

## 41. Feasible ratio Q, not a minimum

Q = 222/1000: N_g is the least integer with 6*ceil(Q*N_g)*N_g >= K_T*m, and N_b = ceil(Q*N_g). No search over all positive integers N_b, N_g is claimed. Kernels stay 916 and 590.

## 42. Lanes 2^19, segments 2^22, wave 2^33

Heads are 2^59 per owner (key layout: s17). Round share 172590/2^19. Baby 7230079554123/6553600000, giant 8384332354123/6553600000. SEG = 2^22 gives 2^35 segments on an EMAX class. Segment field 35 bits (29+35=64). Pool 64+64+9+63+35=235 of 256 (base index 63 bits since s43). Chain id 63+35=98 of 128. Directory charged at EMAX/SEG+2.

Wave 2^33 steps. Two records per step on 2^19 lanes is 2^53 per generating node, mean 2^23 over 2^30 links, half of CAP1=2^24. CAP2 and the 2^41 owner wave stay as in s18. Asserted: 2*WSTEP*LANES = 2^23*2^30.

Totals, score, preprocessing and memory: s4 and s43. H1-H5 and H7 remain open.

## 43. Run gate; sizes at K_T = 261/50

**One run.** The attack is one cold run on one target. It first draws one fresh uniform 64-bit algorithm coin u (one RNDR read at one coordinator node, before any other step of source-s2; all cold setup follows). If u >= K = 7236988362665215240 the run stops there and outputs failure, at charged cost GATE = 64 instructions. Otherwise it releases the first global barrier (one of WAVES' 64 spare barriers, charged in the barrier rows), which starts every other node, and runs to the first verified s (s44) or to the end of the giant phase. A closed gate is a failure, never a success. Apart from u the construction is deterministic.

**Success.** The opened construction's success event S is a function of the challenge t alone and u is independent of t, so P(success) = P(u < K) P(S) >= (K/2^64) lo_run, lo_run = 0.99409171719287617267... (s2). evidence.gate asserts that K satisfies (K/2^64) lo_run >= 0.39 and K - 1 does not; (K/2^64) lo_run = 0.39000000000000000005...

**W(n).** W(n) is the expected total work of this one run over the probability space the rules fix (uniform challenge, fresh algorithm coins), every trial, failure and closed gate included. Pointwise work(u, t) = 64 + 1[u < K] work_run(t), so W(n) = 64 + (K/2^64) E[work_run] <= 64 + (K/2^64) E_run = E_num/2^64 instructions, with E_run = 219994460429955022321640963315313057745 the bound of s44 on an opened run's mean work and E_num = 1592097349982397680518748905618285287073417755790385337224; in units W(n) <= 2^115.055045677162760 for 2^20 <= n < r. What is bounded is this mean. It is not a bound on one run: an opened run does at most W_run + 64 instructions, 2^118.032068814159745 units (s4: giants never stopped, every loop at its hard cap, every overflow a charged failure), certified separately (opened_run_log_bound), and that maximum is above W(n).

**Sizes re-derived at K_T = 261/50, Q = 222/1000** (asserted in evidence.py ledger; giant counts at N_g + EMAX). Pool 2^57 records per owner, mean N_b/2^60 < 2^56.25, wave threshold 2^57 - 2^41, s18 exponent above 13.3x10^15 nats; owner map heads [0, 2^62), pool [2^62, 2^63), so records and s19 links stay < 2^63 (s22), owner < 2^63.0001 bytes in s38's unsigned flat 64-bit model. Bases 2^62.80 in all, below 2^63 (the base index is 63 bits; s42), 2^42.80 per character < 2^43 (H2's cap, sieve-checked), at most 2^16 per generating node; I's 2^64.28 primes: margin above 3. False matches < 2^44.24 < 2^45 - 2^20 (candidate cap 2^45), Chebyshev keeps P(> 2^45) < 2^-40; owner candidate mean 2^-15.9, candidate-link mean 2^-35.9; s11 two-key blocks < 0.86 per owner. Diagonal < 2^-108. H5 conditional mean K/coverage < 5.25. WAVES = 13836394 < 2^24. Memory < 2^123.0005 bytes; opened-run preprocessing < 2^115.4155 units.

## 44. Round-stopped giant phase

Source-s2 steps 5-6 run per round, one candidate round (s4.5) per giant round; kernels, records, keys, routing, caps and planner are unchanged.

**Rounds.** A generating node walks its giant bases one after another (source s2); round j is the j-th giant base of every generating node, a union of complete class walks, fixed before the run and independent of the challenge. At most RMAX = 2^16 rounds; a node without a j-th base waits. Round j, in this order: (i) at most 161 emission waves of the unchanged wave protocol (source s4.3); (ii) 2 drain waves, after which every round-j giant record has been looked up at its owner; (iii) the status pass and the candidate and reply rounds of source s4.5 on the candidates appended in round j; (iv) a close barrier at which the found flag (set only by an owner whose scalar test passed and whose recovered s verifies against A_d and T, source s1) and the fail flag (any cap overflow or failed check) are OR-ed up the barrier tree and released down; every node halts if either is set. Round j+1's emission starts only after round j's close release. So a run executes rounds 1..T, each in full: nothing is in flight at a stop, nothing is emitted after it, and only a verified s or a failure stops the run, never a false candidate. Candidate state is cumulative: owner lists (2^16 slots per run) are append-only with a delivered cursor, the status pass checks the running total against 2^45, and a link carries per round at most its whole-run load, so every cap event is contained in the source's whole-run event and the failure union of s2 is unchanged.

**Charges, at RMAX rounds, never scaled.** WAVES_G = RMAX x 171 (161 + 2 + 8 per round), so WAVES = 13836394 < 2^24 (asserted); the per-wave rows (link words, barriers, s18, s21, s24) grow with it. Candidate-link wave counts are charged every round (RMAX x 4 traversals). Round control is the listing ROUND_CTL (batch_audit.py), executed on toy trees by audit_round_ctl in every ledger run: no flag, a found flag, found and fail flags at two leaves, a 1024-leaf tree, and the round cursor at RMAX - 2 and RMAX - 1; the root holds the OR, every leaf sees it, and a node costs at most 8 + 10 instructions per round. The row charges 24 for every leaf and every barrier node in every round. Per-candidate rows stay at the whole-run 2^45 cap. Every row except the giant records, the allowance, the table rows and the walk rows keeps its every-run value in E_run.

**Stopping bound, through found.** Let G_<j be the giant points of rounds before j, R_(j-1) = U.{b/g : b in B, g in G_<j} and T the last round run. If t is in R_(j-1) and no failure event of s2 occurs, some round before j holds a true pair; a true pair always yields a candidate (exact orbit key; s17, s25), its own round resolves it, s verifies and the run stops before round j. Hence T >= j implies that t is not in R_(j-1) or that a failure event occurred: P(T >= j) <= 1 - |R_(j-1)|/m + FAIL, with FAIL = 2^-20 + 7 x 2^-40 + 2^-108 the failure terms of s2. t is uniform on H and R_(j-1) is challenge independent; no coin, ordering or independence assumption enters.

**H1' (prefix coverage; in the claim's conditioning).** For every round prefix, |R_j|/m >= 1 - e^(-K_j) - 0.0005, K_j = 6 N_b |G_<=j|/m. H1 is its last case. Support and limits: source-s7.

**Expected giant work.** Rounds are fixed sets, so the work w_j of round j does not depend on T, and E[giant work] = sum_j w_j P(T >= j). A giant point costs a fixed part plus 10 per expensive leaf (s37). Every giant base walks a complete class plan (the run emits up to EMAX points beyond N_g for this, charged), each class meets the runtime-checked floors F and ELL (s30, s36), a root costs more than a pair child and a full child more than a leaf, so the fixed part is at most the floor price per point; transport (three lookups) and the beta keys are per point. H7 is cumulative at every barrier, so any conforming run's total leaf cost is at most KXC = 580 + 10 x 7/32 = 582.1875 per leaf; leaves enter linearly, so the mean leaf cost is at most KXC x E[leaves] <= KXC x EG, and E[giant records] <= EG x g with g = 8384332354123/6553600000 = 1279.3476, the every-run price (BABY_PT and W_run already use KXC; the predecessor's F x ELL x 7.8125 = 7.524 per-point premium is not part of the mean). The allowance (s3, 1 per walked point) is priced likewise: N_b + EG + 2^116 (s3(v) pass control < 2^115.1), 2^117.68 not 2^118.52, asserted. With x_j = |G_<j|/N_g and K >= K_T, P(T >= j) <= e^(-K_T x_j) + 0.0005 + FAIL. e^(-K_T x) decreases, a round holds at most MAXR = GENS x EMAX = 2^104 points and the run at most N_g + EMAX, so the left Riemann sum gives sum_j |G_j| P(T >= j) <= EG = ceil(N_g ((1 - lo)/K_T + 0.0005 + FAIL)) + MAXR + 2 EMAX = 84589544470433856447439874622616398 (0.2203 of N_g), lo the rational lower enclosure of e^(-K_T).

**Tables.** The two table rows (gap tables; strand starts and base points: s33, s4) are per base. A node computes a giant base's Q_a' = [a']A_d, strand starts and gap table only after the previous round's close release, inside that base's round (step 5 as step 4, per round; a node holds one table, gen_mem's single GMAX x 192 term), so a run stopped after round T builds none for later rounds. The build, TROUNDS = 134217779 lane rounds plus one strand start per lane (SM = 8192000 instructions), fits in the round's emission slack 161 x 2^31 - LANE_OPS = 2142240768 steps (giant waves 2^31 steps): wave counts are unchanged. Each giant base walks a whole class of at least EMIN points (H2), so round j has |B_j| <= |G_j|/EMIN giant bases and the mean number walked is sum_j |B_j| P(T >= j) <= EG/EMIN; baby bases are all walked. The table rows are charged at BE = BASES_B + ceil(EG/EMIN) = 2706963705881522956 bases (expected_bases), not BASES = 8007881082858055318: 2^118.69, not 2^120.08, instructions (their sum divides by BASES, asserted). **Walk rows.** Padding, side-output fixed calls (s20), restarts and chain ids (s9) are per base or segment, run while the base walks, inside its round; the same bound gives SG = ceil((N_b + EG)/SEG) + BE x LANES segments (expected_segments) and rows (BE LANES SEGLEN + 2 LANES LANE_OPS)(PAIR + 2 BETA) + (BE + 2)(384 LANE_OPS + 64) + SG (SMX + 20) (walk_rows_expected), asserted below their every-run sum. E_run = W_run - (records, allowance, table, walk rows) + N_b x baby + EG x g + (N_b + EG + 2^116) + (table rows) x BE/BASES + (walk rows at BE, SG) = 219994460429955022321640963315313057745 (asserted below W_run).
