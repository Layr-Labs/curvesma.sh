# BLS powers: run gate, interleaved round-stopped babies and giants, Montgomery rows, table reduction, s52 arithmetic tightening (expected work 114.555827210730)

## Complete claim and source scope

The complete claim is claims[0] of submission.json: one cold run that draws a gate coin (K = 7202280856715389586; a closed gate is a charged failure), walks its babies and giants in interleaved rounds and stops at the first verified s; success at least 0.39; W(n), the run's expected total work, at most 2^114.555827210730 units (exact score 114555827210730/10^12); every run at most 2^117.530378341483 units; conditional on H1' (s45) and H2-H5 of source-proof.md (H7 is retired by s47). s2, s4 and s43-s52 derive every figure in it. The flat unsigned 64-bit node address and communication model is analytical, not hardware.

Source scope. This candidate extends accepted submission e983f8e3 (kt-opt-v4-20261118, 114.579559759048) and its lineage cb4e8dc0 / 40ea9361 / 38c6ea6e / 65e66aa0 / b5b30ac6 (complete credits in source-proof.md). No prior review is used as a premise, and no prior success hypothesis is closed. The s52 revision of b574037 modified only five pure scalar AArch64 arithmetic listings (s52, s52-proof.md). This revision changes only the design point (K_T, Q) and the N-dependent quantities that follow from it; the arithmetic listings, schedules, memory, DC ZVA (s50) and hypotheses are unchanged.

Parameters: r as protected, d=2^20 dividing r-1, m=(r-1)/d, six-element GLV/negation classes; Q, N_b, N_g: s2, s45. No previous hosted decision proves H1' or H1-H5.

## 1. Inherited arithmetic derivation; parity output superseded by s40

The kernels derive from the measured listing measure/side_lanes.S (SHA-256 b4c7a7c2d17e6decfc8c0ee93f1def949602ad9d7e2ee0edff4a7ec3c8e4be72; measurement.log: GCC 14.2.0, Valgrind 3.24.0, blst f7a88067, (Ir(8192)-Ir(4096))/4096 = 922.925 instructions/lane for the then 929-instruction listing, field multiplication 383, harness self-checks passed; no new native measurement is claimed). side_lanes(args, n): args[i] = {out, P, T, buf}, args[n] = {NULL, NULL} ends the lane loop; it writes fully reduced x' of P-T to out+0 and y' to out+48 (rows: s46). Inputs are canonical Montgomery residues (P a kernel output, root or strand start; T a gap-table entry with v2 = p - x2 at T+96, n2 = p - y2 at T+144; base-point and strand-start entries get v/n inside their scalar-mult charge, SM = 8192000/point).

Arithmetic. R = 2^384, p/R < 0.1016; the word-serial Montgomery product returns S = (ab + qp)/R < ab/R + p without its final conditional subtraction (blst's __mul_mont_384, __sqr_384, __mul_by_1_mont_384 inlined, n0 in x18, p in x5..x10). From mu (s46): V = S(mu^2); X = V + v2 - x1 (as V + e, s32) reduced to x' by a sign branch or one or two borrow-checked subtractions; w = x' + v2 < 2p; z = S(w, mu); y' = mu(x' - x2) + y2 = z - n2, reduced by one branch. Code after the x' reduction is emitted twice, so no x' path jumps through shared code.

Derivation from the measured 971 listing of 99e03e91/fd385c98 by asserted text transforms. K1, value-bound carry deletion. Lemma: in a word-serial Montgomery product with a < 2p and any 384-bit b, acc_j < a + p before round j and every pre-shift sum acc_j + a b_j + q_j p is below 2^64 (a + p) (induction from acc_0 = 0: (a + p - 1) + (2^64 - 1)(a + p) = 2^64 (a + p) - 1). With a + p < 3p < 2^383, word 7 of every partial sum and word 6 of every shifted accumulator are zero, so their carry words are deleted (18 per product) and the next round's first read of word 6 reads xzr. The kernel's remaining product has a = w < 2p; in the square mu < 2^382 puts the top word below 2^62 and the cross sum below 2^703, so doubling never carries into word 11. K2: b words load in ldp pairs with a partner register dead inside that product. nocarry_witness.py checks the deleted words at extreme and random operands.

Fixed per call: s20 (48 charged per call); side-call row (BASES+2) x (8 x 48 x LANE_OPS + 64), both phases. Current kernel counts and point prices: s46, s47 (pinned).

## 2. Opened-run success floor (s43, s45)

An opened run's floor is lo_run = 1-exp(-K)-0.0005-2^-20-7*2^-40-2^-108, where K = 6 N_b N_g/m >= K_T = 73909/10000 counts the two side targets and the subtracted terms are the H1' allowance at the full sets (H1), the modelling budget, the H3/H5/cap-tail union (the former s12 member now covers the per-round comparison budgets of s45) and the diagonal event (source-s7, s11; s43). Each side emits its target plus the rest of its last class walk (at most EMAX more), which only adds ratios: H1' at the full sets gives coverage at least 1 - e^(-K_full) - 0.0005 with K_full >= K.

The floor is evaluated exactly: evidence.py encloses exp(-K/8) by the alternating Taylor series (60 terms, both partial-sum brackets) in rational arithmetic, takes eighth powers (bracket width below 10^-48) and rounds lo_run down to the 10^-60 grid: lo_run = 0.99888220577219534119... K_T = 73909/10000 and Q = 156921/125000 are design values chosen by the fresh joint ledger scan of s51b at this candidate's prices, not proved minima (s41). The giant pool binds: its per-owner mean 3 (N_g + EMAX)/2^62 is 2^56.99985 against the wave-check threshold 2^57 - 2^41 = 2^56.99997, a slack of 0.0000916 of the mean (many binomial standard deviations; asserted). A giant point's three keys are not independent, so the Chernoff bound takes each point as one block of at most 3: exponent 2^27.1 nats against 82 (2^28.6 for independent keys; asserted). The slack would still be 5.62e-08 after a further full round of MAXR = 2^104 points (asserted), and K_T + 0.0001 at this Q fails the pool assertion. N_g = 221514324027853722891620852279381995 is the least integer with 6 ceil(Q N_g) N_g >= K_T m (asserted both directions), N_b = ceil(Q N_g) = 278081993926198672399008286084263217, so K = 7.390900000000000000000000000000... >= K_T. Every N-dependent quantity is regenerated from N_b + EMAX and N_g + EMAX by the same code with its asserted exponents (list: s43, s45); tails are evaluated at these values, not inherited.

## 3. Explicit pool-cap charge and allowance 1

The one-instruction-per-point allowance covers residual bookkeeping; named pool-cap checks and key/index access costs are charged explicitly.

Residual bounds, retained from the base. (i) Pool-cap checks: charged exactly (wave row of s18, both pools), residual zero. (ii) Key reloads: zero (key stored to slot, reloaded by listed ldp+ldr, travels in-record). (iii) Head index: one lsr and one add (s45), residual zero. (iv) Cap enforcement: charged exactly (the pair of s25 per comparison, per-round budgets of s45). (v) Pass control below 2^115.1: owner-init unrolled 2^20-fold (control below 2^99); all other bulk rows below 2^112 in sum (at s45's N the largest are gap-table entries BASES x GMAX = 2^109.67, planner scan visits 2^108.58 and augmentation ancestor steps 2^108.4, on top of the former 2^109 budget) with at most subs+b.ne control per level over four levels (below 8x sum, below 2^115); vs retained N_b+N_g (2^118.58). The allowance covers only purely-unknown misses.

## 4. Exact complete ledger

Unchanged source dimensions, quotas, capped schedules evaluated by evidence.py at (KT, F, ELL, KERNEL=778 (s52), KX=211 (s52), BETA=138 (s52), BATCH=9209/8 (s52), the s28 share, N_b, N_g; s38 tables; s45 owners). Every row below is an exact integer instruction count, rounding at original locations. The every-run total of an opened run (both sides walked to their targets, every loop at its hard cap) is W_run = 479953627230568756514720629955568406370 instructions; the bound on its mean is E_run = 156393908636956402415832412106483204887 (s45, s52); the 64-instruction gate is added in s43. Challenge-independent preprocessing (sieves, planner, plan copies, beta tables, owner initialisation, windows) is below 2^105.3247 units, included in both; babies are no longer preprocessing (s45). Transport per record GEN 12 + RELAY 770/64 + owner body 770/64 (s45).

- records: N_b + N_g points: 475141853002071973802586736158879480632.
- gap tables: 1507665726991305433247600194045512135.
- beta tables: build (s38): 7605903601371637215993159221248.
- reduction tables (s47): 882284817761208846683140253548544.
- Every remaining row and its exact value is in the declared experiment output (ledger.rows_instructions); no row is dropped or amortized.

In E_run the records, allowance, table, walk and comparison rows are replaced by their expected values of s45 (expected_baby_points, expected_giant_points, expected_bases, walk_rows_expected, comparison_rows_expected in the output); every other row keeps its every-run value.

Memory and all per-node layout assertions are retained under the s8 slot size, owners per s45: total below 2^125.0002 bytes; the generating node and the owner below 2^64 bytes, every other node below 2^63 (s38). No free precomputation, future-target amortization, or wall-clock parallel discount.

Score grid is 10^-12: evidence.py reports the least twelve-decimal ceiling 114555827210730/10^12 of the exact log2 of E_num/(2000 x 2^64), E_num = 64 x 2^64 + K E_run (s43, s45, s52; 90-digit Decimal enclosure). log_bound.py encloses ln(x) by 64 rational terms of 2*sum(t^(2j+1)/(2j+1)), t=(x-1)/(x+1), plus the positive geometric tail; similarly ln2; divides in the safe direction. Applied to E_num/(2000 x 2^64) it proves the declared ceiling and rejects both one-step shifts; applied to (W_run + 64)/2000 it proves the every-run ceiling 117.530378341483.

## 5. Success and honest limitations

Points, stored keys, routing, candidate decisions and recovered scalars match the algorithm. No section s1-s50 changes a source success event: each preserves canonical point/key results and source success events, with internal representatives justified in s49; kernels output exact x', y'; s30's x-keys with three giant lookups test exactly the U-orbit match (orbit lemma in s30). N changes only through s2. H1-H5 keep the source-s7 definitions and budgets; H6 is superseded by the runtime-checked per-class floors (H4: f_c >= 998/1000, H4L: leaf share >= 193/200, s36), checked before any group operation. These hypotheses remain extrapolations requiring review.

An opened run's floor at the realized K (s2) is lo_run = 0.99888220577219534119...; with the gate (s43) success is at least (K/2^64) lo_run = 0.39000000000000000003870... >= 0.39 (exact rationals; K - 1 fails, asserted); declared success is 0.39. All failures and capped work remain in the complete ledger.

The declared run adds curve recovery with interleaved stopping, transport and the two-sided census to the ledger audits (mode all adds the planner and h2). Finite cases only; no full-size attack run.

## 6. Superseded pooled floor

H6 is withdrawn. Use the runtime-checked per-class floors of s30.

## 8. Production lane-stack layout (128-byte slots; 192 since s30)

Supersedes the source's 96-byte lane-stack stride. Slots are 192 bytes (s30; 128 at k7), contiguous from a 16-aligned base: x [0,48), y [48,96), keys [96,160), pad; every listed ldp/stp keeps form and offsets, slot addresses use one hoisted madd. The s40 full kernel writes [out,out+96), leaving the key region to beta_keys; gap-table entries are 192 bytes (x, y, v, n; s1), read-only to kernels; batch buffers are separate. Since s46 full-child slots are 672 bytes (rows around the 192-byte slot), entries 912 and batch slots 496; memory term LANES x (2 DEPTH + 3) x 672. kernel_check.py's stride audit (s1): 8 out slots in one contiguous buffer under random canaries, outputs exact, write footprint exactly the out ranges plus the frame; since s46 the audit uses 672-byte slots with rows checked, and an overlapping 640-byte stride trips it.

## 9. Record-path re-derivation: pointer buffers, honest matching, explicit setup

This section re-derives every per-record instruction from scratch, superseding the source s5.2 transport counts. The charged record sequences (GEN12 s10, RELAY s9, owner body s45, CMP7 s25, loops s29) are asserted by evidence.py audit_record_path (lengths, bits, immediates, pool fit, division proofs); the proto toy executes the pointer mechanism on every record buffer (pre-write check, +40 advance, exact derivation, reset).

Pointer buffers. Sender per-link state is (ptr, end): each record's three stores post-index ptr by 40 bytes; entries are 16 bytes (one ldp), addressed by one add-with-shift. The cap check precedes the write (cmp ptr, end; b.hs fail). At wave end the sender derives count = (ptr-base)/40 (sub, lsr #3, umull/lsr /5 with magic 0xCCCCCCCD, exact by m*5 = 2^34+1 and the asserted bound; 8 with the count store and reset), writes the per-parity count word and resets; receiver side 3. LINK_WAVE is 16 (8 + 3 + 5 margin); the per-node term is NODES*24 for dmb, counters and constants.

Honest matching: this section's record key, pool layout and owner/compare/cap counts are superseded by s11, s17, s19, s25, s29 and s45 (current bodies in evidence.py, executed by audit_il). Four new explicit rows: window-table init (8/link, 4 used), chain-id formation (3/segment), candidate-link wave counts (9 on 4 traversals: /40 magic, /64 shifts, /56 with proven /7 magic), candidate formation (16/match, 11 used). Caps, routing, counts, events unchanged, so H1-H6, N, success, memory untouched.

## 10. k9 repair: chain-word memory schedule

Two chain-id words formed once per segment into a per-lane 16-byte slot (LANES slots per generating node, in memory term); every GEN reloads them. Formation 4 per segment (ubfx+bfi+lsr plus stp); GEN11 becomes GEN12 (ubfx t, entry address, ldp entry, ldp chain words, cmp, b.hs, three stores, ptr store-back, key-reload/loop pair). No chain word crosses records in a register (except within one giant point's three consecutive GEN copies, which no call separates; s34), so lane interleaves and side-output clobbering are compatible; allowance untouched.

## 11. Comparison cap per owner

Superseded by s45's per-round budgets. What stands: a run prefix's comparisons at an owner are sum over its head slots of (baby records) x (giant keys) in that slot, each pair compared once (s45); chains above 64 have probability below 1/65! < 2^-302 per slot at the giant-table load 3 N_g/2^119 < 0.9999 (binomial tail at most load^65/65!), below 2^-183 over all 2^119 slots (inside the union); a giant point's keys at one owner form one block (k_g >= 2 has probability 3 x 2^-124 per point, mean below 0.04 per owner); B_o and G_o are binomial with Chernoff tails, and Hoeffding over blocks bounds the sum (s45 asserts the exponents at every prefix scale).

## 12. Extract step deleted

The s45 owner body holds the key field and k64 in registers before its lookup, so no extract step, extract cap or extract counter remains (s21 likewise deleted).

## 16. Augmented planner (pointer; floors and planner since s30)

The s16 augmentation (1/3/5 alternating paths, partner cap 32, length-5 cap 8, exact rollback, OPS_BOUND = 10752 per element, asserted) runs unchanged inside the s30 planner, over centre pairs and in the fallback; its floor and tables are superseded by s30.

## 15. Exact restart price SMX from the shipped kernel bodies

Supersedes the s13 price; only the restarts row changes. Each restart computes [b]Q_a, b in [1, X2] (<= 77 bits, from the plan entry), by branchy Jacobian double-and-add (dbl-2007-bl, 2M+5S; madd-2007-bl with affine Q_a, 7M+4S) on an all-ones worst path in Montgomery residues. Strands, base points and candidate tests stay at SM = 8192000. SMX = 48 + 380 + 5 + 76 x (2340 + 11 + 3656) + 6 + 172520 + 1231 = 630722, an upper bound on every b; row count x (SMX + 16). Checked by smx_audit.py.

## 17. Record key layout: 57-bit head, 9-bit field (s45 owners)

The 192-bit record key (words w0, w1, w2) is split as owner bits 0..61, field f9 = bits 62..70, head bits 71..127, k64 = bits 128..191 (62+9+57+64 = 192, asserted). The head index is w1 >> 7 (HEADS = 2^57, asserted); owner routing uses t = bits 0..29 and u = bits 30..61.

The owner body (s45): extr x10,x1,x0,#62 puts key bits 62..125 in x10; bfi x4,x10,#55,#9 deposits f9 into w4[55..63]. The chain word w4 holds the 35-bit segment in [0..34] (s42) and the pool link is its own word, so bits 55..63 are free; bits 35..54 are never read.

Compare (routine: s25): and x14,x14,#0xFF80000000000000 keeps w4[55..63] of the stored record; cmp x14,x10,lsl #55 discards x10's bits 9..63, so exactly key bits 62..70 meet, equal exactly when the two f9 fields are equal, whatever the head bits. ccmp then compares k64. With owner and head equal by routing and chain, a comparison succeeds exactly when the two 192-bit keys are equal.

Audit (every ledger run): audit_packed_match checks this compare word for word on 3001 keys with junk in w4[0..54], with single-bit flips inside and outside the field and in the head.

## 18. Wave-granular pool-cap check (two pools since s45)

One check per owner per wave and per pool. An owner inserts only through its 2^20 inbound relay links (source-s4.3), each carrying at most CAP2 = 2^21 records per wave, enforced by the relay's charged pre-write cap check (overflow is a run failure in the same budget). So an owner inserts at most 2^41 records per wave into either pool (asserted); checking each pool pointer against its threshold at each wave start guarantees no wave overflows it. Baby pool 2^56 records at [2^61, 2^62), threshold pointer 2^62 - 2^46; giant pool 2^57 records at [2^62, 2^63), threshold 2^63 - 2^46 (each one movz + movk, asserted). Per-owner means below 2^55.7430 and 2^56.9999 against 2^56 and 2^57, counted at N + EMAX with three keys per giant point; the wave check uses 2^56 - 2^41 and 2^57 - 2^41, and the Chernoff exponents are asserted above the 82-nat budget with a giant point's three keys counted as one block of at most 3 (2^49.5 and 2^27.1 nats). Per owner per wave and pool: ldr, ldr, subs, b.hi plus str at wave end (5); once per owner and pool 6. Row = 2 x OWNERS x (5 x WAVES + 6), WAVES = 11214912 (s45); in every-run value.

## 19. Full-address chain links

This section restates the link word (superseding the s7 56-bit-link sentence): head entries and chain links store full 64-bit record addresses, 0 = end or empty; records load as [addr,#imm] (comparison CMP7 since s25, owner body and head store in s45). Layout: pools after the head array of 16-byte head pairs (s45), records at pool base + 32 and above, addresses below 2^63 (asserted). Pool use: s42; word sizes, memory, owner init unchanged. Links never leave the owner (chain ids are base+segment).

## 20. Itemized side-call fixed price: 40 per invocation

Per call: 32 fixed (prologue 22 + epilogue 10, asserted every call) + 1 exit branch + 5 staging (two moves, lsl, terminator stp, bl) + 2 round control = 40 (asserted), no spills (loop variables in x19-x28, callee-saved preserved), lanes unrolled so no per-lane control; 48 per call charged since s37. Per-base-per-phase setup bounded by 32, charged as 64 x (BASES+2). Row: s1, both phases charged.

## 21. Deleted (s12, s45)

## 22. Pool pointer (owner body: s45)

The owner's pool counter is a pointer p to the next free 32-byte record (the s9 sender device at the owner): the s45 body writes the record with stp x3,x4,[x22,#16] and stp x2,x13,[x22],#32 (k64, link = old head; the only write of p, after both stores through it) and the head word with str x22. One pointer per pool; the per-round kind selects it. s18 compares p with its threshold at each wave start.

## 23. Superseded by s45 (one owner body inserts and looks up)

## 24. Transport loop unrolling (64-fold since s29)

Every transport loop runs UNT = 64 body copies back to back sharing one subs + b.ne, its only cross-record instructions; each copy reuses registers after its stores, and scratch x5 is absent from the bodies (asserted). Setup is and x5,x23,#63 + lsr x23,x23,#6 + cbz; the tail is a binary ladder (tbz x5,#j skips 2^j copies, j = 5..0). Per record RELAY 12 + 2/64 (asserted against the s9 mnemonics) and owner body 12 + 2/64 (s45); per invocation at most 10 (UNROLL_FIX, asserted), charged in the row 10 x ((LINKS1+LINKS2) + 2 x OWNERS) x WAVES. GEN12 has no per-record loop pair (driver control in s3(v)); the s25 chain loop alternates two copies.

## 25. The chain comparison, as an executed routine

After the s45 owner body, x10 (extr, key bits 62..125), x2 = k64 and x11 = the other table's head record address are in registers; a nonempty head walks the chain inline before the next record. A pool record is k64 at +0, link at +8, chain words at +16 and +24 (s22). The routine is two copies that alternate the record register (x11, x13), so no move is needed between records; copy A (evidence.py CMP_A; copy B swaps x11 and x13 and ends cbnz x11, A):

    subs x25, x25, #1          ; s45 per-round comparison budget (CAP 2)
    b.lo fail                  ; budget exhausted: run failure
    ldp x12, x13, [x11]        ; k64, link
    ldr x14, [x11, #24]        ; w4
    and x14, x14, #0xFF80000000000000
    cmp x14, x10, lsl #55      ; f9 (s17)
    ccmp x12, x2, #0, eq       ; k64
    b.ne next
    str x11, [x26], #8         ; match: record address to the candidate list (candidate formation row, s9)
    next: cbz x13, done        ; end of chain: on to the next lookup record

Per comparison: CMP 7 (ldp, ldr, and, cmp, ccmp, b.ne, cbz or cbnz) plus the 2-instruction cap pair, every-run rows OWNERS x cap x 7 = 1937888384957516273534622887040778240 and OWNERS x cap x 2 (cap 60030365628724705 per owner, s45), charged in E_run at their s45 expectation. Pool bases are 16-aligned and records are 32 bytes, so every ldp is aligned. No signed condition; addresses stay below 2^63.

Executed: audit_il (s45), every ledger run.

## 26. Batch round: both children by the derived kernels, no reference addition in the walk

Points, keys, records and events are unchanged. Round (one lockstep step of a node's 2^19 lanes, s42): the driver writes argsB[i] = {out, P, T, buf+496i} per lane op (out = the P+T child, or the root's point) and argsA = {P-T out, P, T, buf+496i} per pair (since s30 split into full and leaf arrays by the leaf flags), each ended by a 32-byte NULL record; padded lanes use the fixed pair (Q_a, T[2]), x distinct, so every denominator is nonzero. batch_round(argsB, buf) runs Montgomery's trick over d_i = x1 + v2 (s27) and leaves both slopes per lane op (s46); side_lanes(argsA) writes P-T and side_lanes_plus(argsB) writes P+T (side_lanes on -T = (x2, p-y2, v2, y2): its n2 loads read T+48/64/80). Counts: s46. PAIR = 2 KERNEL + BATCH + 28 + s for two points (plan entry 8, P address 1, T[g] address 6, buffer 1, two out addresses 2, two argument records 4, s30 routing 6); ROOT = KERNEL + BATCH + 21 + s (gap entry 5, T[gap] address 6, P and out addresses 2, buffer 1, argument record 2, loop 5). Side-call row: s1. Re-walks keep 3300 + 232864. The reference addition (2000) is not used anywhere in the walk; every product, the inversion and both kernels are charged at their listed and executed counts, nothing twice.

## 27. Lazy batch round: reflected difference and unreduced inverses

The forward pass forms d_i = x1 + v2 (v2 = p - x2 at T+96) in (0, 2p), d_i = -(x2 - x1) mod p, nonzero as x1 != x2. Montgomery's trick is linear in its running value, so the round negates the canonical Fermat output once (t = p - t) and every backward inverse is inv_i = 1/(x2 - x1); inv_i stays unreduced (bounds: s46).

## 28. Exact round share and unrolled batch loops

A walk round is one lockstep step of a node's 2^19 lanes, each doing one pair or root operation, idle lanes padded and charged (source s2 step 4, s3 padding lemma; s20): exactly 2^19 lane ops and one fixed cost 172597 + 5 = 172602 (s26, s46), so each lane op carries s = 172602/2^19. Every executed lane op sits in a row charged at PAIR or ROOT (records at the floor F; padding; remainders), so the shares cover every round. The forward loop's only per-element control is cbnz x0, the backward loop's ldr x30 (end slot, from the frame since s47) + cmp x2,x30 + b.ne, each over 32 copies of its body; n = 2^19 is a multiple of 32 (asserted). The listing (about 2.2 x 10^5 instructions with the Fermat chain) sits inside gen_mem's 2^20-byte code term.

## 29. Transport loops 64-fold (s24)

Candidate delivery is eight relay hops at 770/64; evidence: audit_il (s45).

## 30. x-only leaves, beta-triple giant lookups, leaf planner and asymmetric sides

Supersedes the parity record key, source s3 matching, the s16 floor, GMAX (2^47) and DEPTH (2^17), N, the slot of s8, the per-point price and the side calls per round of s20/s26.

Keys. Every record key is limbs 1..3 of the canonical Montgomery residue of x' (out+8..31). [-1](x, y) = (x, -y) and [lambda](x, y) = (beta x, y), beta fixed by [lambda]P0 = (beta x0, y0) (asserted in beta_run), so the U-orbit of Q has x-set {x, beta x, beta^2 x}; conversely x_P in that set gives y_P^2 = x_P^3 + 4 = y_Q^2, so P is in U Q. A baby x-key equal to one of a giant's three keys is exactly the orbit match the parity key tested, and at most one of the three can match (x = beta x forces x = 0, a 3-torsion point; 3 does not divide r). Babies insert one record per point; giants send three lookup records per point (keys of x', beta x', beta^2 x', one chain id; GEN unrolled three times per point at key offsets +8/+104/+136; a giant re-walk keeps positions where any of its three keys equals the candidate key). Coverage events, K, success and the 96 recovery tests are unchanged. Giant waves are a quarter of a pair wave (2^31 steps since s42), so a generating-node link still carries at most 2^23 records per wave in the mean (six per step).

beta_keys (batch_audit.beta_run; tables since s38): per giant slot, limbs 1..3 of canonical beta X' at +104 and of beta^2 X' = -(X' + beta X') at +136, walking the side-call argument arrays (the driver stores nothing extra).

Leaves. A pair child that is no pair's midpoint is only keyed and runs the leaf kernel (kernel_check.leaf_text: the kernel cut after the x' store; s47), bit-identical in x'. Midpoints and roots run the full kernel.

Driver and calls. Plan entries carry two leaf flags (kind 2 + gap 48 + three 19-bit slots + 2 = 109 of 128 bits); per pair lane op two ubfx and two pointer ldr/str route the children to the full or leaf call (ROUTE = 6, listed). Per round four side calls (full/leaf x minus/plus) and four beta_keys passes, charged 8 x 40 on both phases (s1, s20). Pads run full pairs plus two beta passes.

Planner (planner.py): one centre per run of 50 class-list elements; greedy by increasing degree over centre pairs, hosting centres first; the s16 augmentation; a fallback pass over all midpoints; source s3 demotion and emission; every pair passes the cycle, depth and tree rules with exact rollback. GMAX = 2^47 = 128 d^2; WIN = 2^29 keeps rel32 offsets below 2^31. Charged: 6 window scans per element, matching at 64 kappa = 512 candidate pairs per element with at most five tries, two augmentation passes with partner caches; memory adds 3 bytes per element. Each class is planned before any group operation and the run stops unless f_c >= F = 998/1000 and the leaf share ell_c >= ELL = 193/200 (H4, H4L). Measured at kappa = GMAX/(16 d^2) = 8, span X2/(2 GMAX) = 4: the declared payload's spec at d = 4 (mode planner, offline here); d >= 16 by a line-for-line C port of plan() that reproduces planner.py's forests exactly (pairs, roots, leaves) on six declared-type lists, capped and uncapped:

| d | classes | min F | min ell | mean ell | max height |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 4 | 2 | 0.99550 | 0.97349 | 0.97375 | (mode planner) |
| 16 | 4 | 0.99822 | 0.97278 | 0.97298 | 235 |
| 64 | 4 | 0.99863 | 0.97213 | 0.97296 | 931 |
| 128 | 4 | 0.99875 | 0.97305 | 0.97329 | 1062 |
| 256 | 8 | 0.99877 | 0.97306 | 0.97324 | 2172 |
| 512 | 4 | 0.99886 | 0.97311 | 0.97338 | 6931 |

ell is flat in d (0.972-0.974) and F rises with d; F sits 0.00077 below every d >= 256 row and ELL below every row. With the tree cap at TREEMAX = d the leaf deficit is about 5/d and the F excess 1.1/d to 2.2/d (d = 16..128); both put the production effect at TREEMAX = 2^20 below 10^-4, inside the margins (part of H4/H4L). Heights grow with tree size (6931 at 2^19 elements): DEPTH = 2^17.

Sides. With side prices c_b, c_g the unconstrained optimum is N_b/N_g = c_g/c_b, approximately1.23925 at the current prices. Q=156921/125000 is selected by the joint s51b ledger scan (s41), not claimed optimal. The ec toy (x-keys, triple keys, interleaved bases, solved == exact coverage every trial, stop round == first covering round, all scalars correct): the declared input (mode declared, K = 7.3909, Q = 59/50) gives 4/4 at 0.9994 predicted, 4 trials stopping before their last round. Production success rests on H1' and the census, not on toy rates. H3 covers N_b + 3 N_g keys (a giant's three keys are fixed bijections of x', treated as uniform); caps: s45.

Prices (pinned in audit_record_path; current values s45): pair child (BATCH + PAIR_BK + 6)/2 + ELL x KX + (1 - ELL) x KERNEL; baby point F x child + (1 - F) x ROOT + 12 + 2 x 770/64; giant point F x child + (1 - F) x ROOT + 3 x (12 + 2 x 770/64) - 2 + BETA.

## 32. Shared x-sum e = v2 - x1, once per lane op

The forward element stores e = v2 - x1 (at buf+336 since s46); kernel_check.esum_text replaces one asserted block of the measured listing so the kernel forms X = V + e. As integers e is in (-p, p), X in (-p, 2.11p) (s46), and bit 383 of its two's complement is the sign (tbz x16,#63); x', y' are bit-identical. Counts: s46.

## 33. Gap tables by pair levels

Changes only how tables T[k] = [k]Q_a, k <= GMAX, are built (entries, layout, memory and use unchanged). With T[1..M] known (M = 3 from the strand-start scalar multiplications, row unchanged), one level computes the centre c = M + k + 1 as a root T[M] + T[k+1] and then the pairs T[c -+ j] = T[c] -+ T[j], j = 1..k, k <= M - 2: the s26 lane op on (T[c], T[j]) with both children full (the s46 full kernel with rows, minus and plus calls, then the 21-instruction v/n step and the 60-instruction V rows), so T[1..M + 2k + 1] is known. Full levels (k = M - 2) give M' = 3M - 3; the last stops at GMAX: 30 levels. Incomplete additions are safe (c -+ j, c distinct, nonzero, below GMAX < r).

Price per base: every level is charged ceil(k/LANES) + 1 full rounds of LANES lane ops at PAIR (+1 for the centre; padding included) plus 128 per round for calls, and 21 per entry: BASES x (rounds x (LANES x PAIR + 128) + (21 + 60) GMAX) (strands: BASES x (GMAX - 1) x 2050; re-walks keep 2050). Row in s4.

Evidence: batch_audit.run(table) executes a 64-lane level (M = 66, centre [131]Q, gaps [j]Q): the s32 round and both kernel calls give exactly [131 -+ j]Q; evidence.py asserts the 30-level recursion.

## 34. Giant GEN chain-word reuse; head zeroing by dc zva

Prices only. (i) A giant point's three lookup records share one chain id (s30; GEN unrolled three times per point, key offsets +8/+104/+136). Copies 2 and 3 drop their chain-word ldp (copy 1's registers still hold the words; GEN only reads them). Giant transport 3 x (12 + 2 x 770/64) - 2 (s45); the proto's generic per-record GEN 12 stays an upper bound. (ii) Owner initialisation zeroes each 64-byte-aligned head array by dc zva (64-byte blocks, DCZID_EL0.BS = 4, common AArch64) plus one add: 2 instructions per 64 bytes, not one stp per 16. Historical row OWNERS x (HEADS x 16/64 x 2 + 64) with s45's16-byte head pairs is superseded by s50, which also explicitly charges the zeroing loop control.

## 35. EMIN = 2^56 - 2^52 (H2 restated)

H2 has class lists between EMIN = 2^56 - 2^52 and EMAX = 2^57, sieve-checked before any chain starts (a short class is an H2 failure event); bases are charged at N/EMIN per side. Bases at the bound: 2^61.82 baby and 2^61.51 giant, 2^62.69 in all (I's 2^64.28 primes: margin above 3.01; s45); segments per class at least 2^33.9. Support: class-length spread shrinks like (X/d)^(-1/2); h2_sizes (mode all) on toy primes gives min/max of |L_c| d/X from 0.9287/1.1035 (d = 16-256, X/d = 2^10) to 0.9908/1.012 (X/d = 2^14), which extrapolates far inside 6.25% at X/d = 2^56 (heuristic).

## 36. Leaf floor ELL = 193/200

H4L's runtime-checked floor is 193/200 (a class below it stops the run before any group operation). Per-class minima: s30's table (0.97213-0.97349); offline planner.py rows (kappa 8, four classes, d = 4-16, span 4-8) 0.97274-0.97524; one small two-class d = 8 run had 0.9691. Prices only.

## 37. Leaf reduction (H7 retired by s47)

Since s47 every leaf lane ends with x' canonical after one conditional subtraction (paths 213 and 214), so no X >= p count, cap, barrier sum or hypothesis H7 remains; the former leaf-cap row is deleted.

## 38. Eight-chunk beta tables (current table trade)

Canonical Montgomery X < p < 2^381 splits into seven 48-bit chunks and a final 45-bit chunk at shifts 0,48,96,144,192,240,288,336. T_i[c] = beta*c*2^shift_i mod p is canonical. Eight tables interleave at A+512*c+64*i; q*p for 0<=q<8 is at A-512+64*q. Production A=512. Allocate and initialize all eight 2^48-entry columns, including unused entries of the last column, so storage is 2^57+512 bytes and the charged construction is GENS*(8*2^48*24+8*2000+64). BUILD is the supplied 24-instruction recurrence at stride 512; the eight setup constants are conservatively charged 2000 each.

Lookup costs 139 instructions per giant point (128 in the immediate source): loads, exact chunks, addresses,24 table loads, accumulation, reduction, stores, beta^2*X and loop. The sum S<8p<2^384. Writing p=p5*2^320+p_low and S=s5*2^320+s_low, q=floor(s5/(p5+1)) gives q*p<=S and S<(q+1)(p+2^320). Since q<=7 and 8*2^320<=p, S-q*p<2p, so one conditional subtraction yields canonical beta*X. The unchanged beta^2 identity computes -(X+beta*X) mod p. Keys are bit-identical to the source; H3 is not strengthened or replaced.

The large TMR tables of s47 follow the small beta tables, then other generator regions. The exact generator total is 17854630968541906320 bytes, below 2^64 and2^57 bytes smaller than the immediate source. Owners, relays and their lookup/transport layout are unchanged. All addresses are unsigned and do not wrap.

Evidence: beta_run and the composed batch round execute boundary and curve-point cases against exact arithmetic. The new number of chunks changes costs and local table layout only.

## 39. Inherited evidence repairs

Planner fields: 28 bytes base, 8 adjacency, 4 previous link, 3 centre (43 bytes/element, DEPTH = 2^17). The character check certifies order d for 7^m and exhibits k = 256 where log_g(x^m) is the wrong character base; source s2 uses log_(g^m).

## 40. Delete the unused full-point parity key

Section 30 replaced every record key by x-coordinate limbs. kernel_check.drop_parity_text removes the four exit blocks that only emitted the old parity key (each asserted as tbz, subs, sbcs, sbcs, sbc, stp, str); coordinates are stored before them, no consumer reads +96..191 of a full output, and flags and scratch registers are dead at the next lane. Counts: s46.

## 41. Re-scanned K_T and Q; manifest checks

N_g is the least integer with6*ceil(Q*N_g)*N_g>=K_T*m and N_b=ceil(Q*N_g), asserted both ways. K_T=73909/10000 (7.3909) and Q=156921/125000 come from the fresh joint ledger scan of s51b at this candidate's prices, replacing the b574037 design point K_T=14717/2000, Q=5/4 (itself inherited from 40ea9361's K_T=727/100, Q=247/200). The algorithm, schedules, hypotheses H1'-H5 and census scope are unchanged; only N_b, N_g and the gate coin K follow from (K_T, Q), and every N-dependent quantity is regenerated by the same code (s2, s43, s45). No proved minimum is claimed: this is a finite ledger grid, not optimality.

**Manifest checks.** On the declared input combined_lazy.py reads submission.json and asserts that claims[0].score is the ledger's score numerator over 10^12, that the title and statement carry the ledger's W(n) exponent, every-run exponent, K, lo_run (60 digits), K_T, E_run, Q and the preprocessing and memory exponents (each exactly once), and that every hypothesis bound quoted from the ledger (pool means, false matches, owner-candidate and link means, the H5 means, bases, prime margin, slot load) holds against its cap and stands verbatim in proof.md or source-proof.md (outputs manifest_check, hypothesis_bounds). log_bound certifies both exponents from the same ledger integers. An input that overrides KT or ratio skips the manifest check and says so.

## 42. Lanes 2^19, segments 2^22, wave 2^33

Heads: s17, s45. Round share 172602/2^19 (s46). Point prices: s46. SEG = 2^22 gives 2^35 segments on an EMAX class. Segment field 35 bits (29+35=64). Pool 64+64+9+63+35=235 of 256 (base index 63 bits since s43). Chain id 63+35=98 of 128. Directory charged at EMAX/SEG+2.

Wave 2^33 steps. Two records per step on 2^19 lanes is 2^53 per generating node, mean 2^23 over 2^30 links, half of CAP1=2^24. CAP2 and the 2^41 owner wave stay as in s18. Asserted: 2*WSTEP*LANES = 2^23*2^30.



## 43. Run gate; sizes at K_T = 73909/10000

**One run.** The attack is one cold run on one target. It first draws one fresh uniform 64-bit algorithm coin u (one RNDR read at one coordinator node, before any other step of source-s2; all cold setup follows). If u >= K = 7202280856715389586 the run stops there and outputs failure, at charged cost GATE = 64 instructions. Otherwise it releases the first global barrier (one of WAVES' 64 spare barriers, charged in the barrier rows), which starts every other node, and runs to the first verified s (s44, s45) or to the end of its last round. A closed gate is a failure, never a success. Apart from u the construction is deterministic.

**Success.** The opened construction's success event S is a function of the challenge t alone and u is independent of t, so P(success) = P(u < K) P(S) >= (K/2^64) lo_run, lo_run = 0.99888220577219534119... (s2). evidence.gate asserts that K satisfies (K/2^64) lo_run >= 0.39 and K - 1 does not; (K/2^64) lo_run = 0.39000000000000000003870...

**W(n).** W(n) is the expected total work of this one run over the probability space the rules fix (uniform challenge, fresh algorithm coins), every trial, failure and closed gate included. Pointwise work(u, t) = 64 + 1[u < K] work_run(t), so W(n) = 64 + (K/2^64) E[work_run] <= 64 + (K/2^64) E_run = E_num/2^64 instructions, with E_run = 156393908636956402415832412106483204887 the bound of s45/s52 on an opened run's mean work and E_num = 1126392854282846724778895998043138370365995771849475410206; in units W(n) <= 2^114.555827210730 for 2^20 <= n < r. What is bounded is this mean. It is not a bound on one run: an opened run does at most W_run + 64 instructions, 2^117.530378341483 units (s4: never stopped, every loop at its hard cap, every overflow a charged failure), certified separately (opened_run_log_bound), and that maximum is above W(n).

**Sizes re-derived at K_T = 73909/10000, Q = 156921/125000** (asserted in evidence.py ledger; side counts at N + EMAX): owners, pools, heads and comparison budgets: s45. Bases 2^62.69 in all, below 2^63 (the base index is 63 bits; s42), at most 2^42 per character and side (H2's cap, sieve-checked), at most 2^15 per side per generating node; I's 2^64.28 primes: margin above 3.01. False matches 2^44.75 < 2^45 - 2^20 (candidate cap 2^45), Chebyshev keeps P(> 2^45) < 2^-40; owner candidate mean 2^-17.25, candidate-link mean 2^-37.2. Diagonal < 2^-108. H5 conditional mean K/coverage < 7.40. WAVES = 11214912 < 2^24. Memory < 2^125.0002 bytes.

## 44. Round-stopped run

Source-s2 steps 4-6 run per round, one candidate round (s4.5) per round; kernels, records, keys, routing and planner are unchanged. Round kinds and owner tables: s45.

**Rounds.** A generating node walks its bases one after another (source s2); round j is, for every generating node, its next base of the round's kind (s45), a union of complete class walks, fixed before the run and independent of the challenge. At most RMAX = 2^16 rounds (2^15 per side); a node without such a base waits. Round j, in this order: (i) at most 161 emission waves of the unchanged wave protocol (source s4.3); (ii) 2 drain waves, after which every round-j record has been inserted and looked up at its owner; (iii) the status pass and the candidate and reply rounds of source s4.5 on the candidates appended in round j; (iv) a close barrier at which the found flag (set only by an owner whose scalar test passed and whose recovered s verifies against A_d and T, source s1) and the fail flag (any cap overflow or failed check) are OR-ed up the barrier tree and released down, with the next round's comparison budget (s45); every node halts if either flag is set. Round j+1's emission starts only after round j's close release. So nothing is in flight at a stop, nothing is emitted after it, and only a verified s or a failure stops the run, never a false candidate. Candidate state is cumulative: owner lists (2^16 slots per run) are append-only with a delivered cursor, the status pass checks the running total against 2^45, and a link carries per round at most its whole-run load, so every cap event is contained in the whole-run event and the failure union of s2 is unchanged.

**Charges, at RMAX rounds, never scaled.** WAVES = RMAX x 171 (161 + 2 + 8 per round) + plan rounds + 64 = 11214912 < 2^24 (asserted); the per-wave rows (link words, barriers, s18, s24) grow with it. Candidate-link wave counts are charged every round (RMAX x 4 traversals). Round control is the listing ROUND_CTL (batch_audit.py), executed on toy trees by audit_round_ctl in every ledger run; a node costs at most 8 + 10 instructions per round and the row charges 24 for every node in every round. Per-candidate rows stay at the whole-run 2^45 cap.

**Stopping bound, through found.** Let B_<j, G_<j be the baby and giant points of rounds before j, R_(j-1) = U.{b/g : b in B_<j, g in G_<j} and T the last round run. If t is in R_(j-1) and no failure event of s2 occurs, some round before j holds a true pair; a true pair always yields a candidate in the later of its two rounds (exact orbit key; s17, s25, s45), that round resolves it, s verifies and the run stops before round j. Hence T >= j implies that t is not in R_(j-1) or that a failure event occurred: P(T >= j) <= 1 - |R_(j-1)|/m + FAIL, with FAIL = 2^-20 + 7 x 2^-40 + 2^-108 the failure terms of s2. t is uniform on H and R_(j-1) is challenge independent; no coin, ordering or independence assumption enters.

**Tables and walk rows.** A node builds a base's Q_a or Q_a' = [a']A_d, strand starts and gap table only inside that base's round (one table at a time), so a run stopped after round T builds none for later rounds; the build (134217779 lane rounds plus one strand start per lane) fits in the round's emission slack. Each base walks a class of at least EMIN points (H2), so the mean number of bases walked is at most (E_b + E_g)/EMIN + 2 = BE = 2414234128596394668 (expected_bases), not BASES = 7377509947778909855; padding, side-output calls, restarts and chain ids take BE and SG = ceil((E_b + E_g)/SEG) + BE x LANES segments (walk_rows_expected), asserted below their every-run sum.

## 45. Interleaved rounds: babies and giants grow together

Supersedes the baby-first order (source s2, s4.3), the giant-only rounds of s44, the owner layout, bodies and caps (s11, s12, s17, s18, s21-s23, source s4.2) and Q, N (s2, s41, s43). Kernels (prices: s46), planner, records, transport hops, beta tables and the gate are unchanged (K_T, lo_run: s2).

**Round kinds.** Bases are dealt as in source s2 step 3, at most 2^42 baby and 2^42 giant bases per character (H2, sieve-checked), so at most 2^15 of each per generating node. With y_j = |B_<=j|/N_b and x_j = |G_<=j|/N_g the cumulative shares after round j (exact class lengths from the sieves, before the run), round j+1 is a baby round if y_j <= x_j and babies remain, else a giant round. The round table (kind, n_b(j), n_g(j), cap_j) is accumulated during base assignment (one add per base, inside its 64-instruction charge), held at the barrier root and sent down in the close-release word (the owner's load is inside its 8 per round). In a baby round every generating node walks its next baby base; in a giant round its next giant base.

**Owners.** 2^62 owners (t, u) = key bits 0..29 and 30..61; per owner 2^57 head pairs of 16 bytes (baby head at +0, giant head at +8; [0, 2^61)), indexed by key bits 71..127 (s17); a baby pool of 2^56 and a giant pool of 2^57 32-byte records ([2^61, 2^62) and [2^62, 2^63)), records as s22. A baby record inserts its key into the baby table and looks up the giant table; each of a giant point's three records (keys of x', beta x', beta^2 x', s30) inserts into the giant table and looks up the baby table. One 12-instruction body per kind (evidence.il_body; unrolled as s24): three loads of the record, lsr (head), add and ldp of the head pair, extr and bfi (field), stp and str (record and head), stp (k64 and link = own old head), cbz on the other table's head, which enters copy A of s25 directly; no extract step remains.

**Each equal-key pair once.** Rounds are kind-homogeneous and a record is stored before the next record is processed, so for an equal-key (baby key, giant key) pair the record processed later finds the other on its chain and the earlier one did not. Every true pair yields exactly one candidate and every same-slot pair is compared exactly once, so a prefix's comparisons at an owner are sum over slots of B_s G_s. audit_il executes both loops (a64.py) over 12 rounds of both kinds against the full-key oracle (candidate lists, heads, links, pools, step count 12n + 2(n div 64) + 9C + matches + 10, each pair once, budget fail exit, four negative controls). proto runs the wave protocol over six rounds of both kinds and asserts the same exactly-once property.

**Prices and sizes.** Per baby point s30's walk price + 12 + 2 x 770/64 = 1145678673093/1310720000 = 874.0835; per giant point the walk price + 3 x (12 + 2 x 770/64) - 2 + BETA = 1419782993093/1310720000 = 1083.2085 (pinned). Q = 156921/125000 (s51b ledger scan, s41; c_g/c_b is approximately1.23925). Pool means: s18.

**H1' (prefix coverage; in the claim's conditioning).** For every round prefix j of the production schedule, |R_j|/m >= 1 - e^(-K_j) - 0.0005 with K_j = 6 |B_<=j| |G_<=j|/m. H1 is its last case. Support: source-s7 (two-sided census).

**Expected work.** Rounds are fixed sets, so E[work] = sum_j w_j P(T >= j) with w_j fixed, and by s44 and H1' P(T >= j) <= e^(-K x_(j-1) y_(j-1)) + 0.0005 + FAIL with K >= K_T. At a baby round y_(j-1) <= x_(j-1), so x y >= y^2; at a giant round x y >= x^2. A round adds at most MAXR = GENS x EMAX = 2^104 points and a side ends at most EMAX past its target, so the left Riemann sum gives expected baby points E_b = ceil(N_b (J + 0.0005 + FAIL)) + MAXR + 2 EMAX = 90791727675373716712186999345312124 (0.3266 of N_b) and likewise E_g = 72338331172536214682056472866667513 (0.3266 of N_g), with J = int_0^1 e^(-K_T y^2) dy < 0.325987 (alternating Taylor series in rationals, stopped on a negative term below 10^-50; evidence.gauss_int_up). Point prices are every-run prices (per-class floors) and enter linearly. The allowance is N_b + N_g at the every-run value and E_b + E_g + 2^116 in E_run (the 2^116 covers s3(v)'s pass control, below 2^115.1).

**Comparison budgets.** At each close release every owner receives cap_(j+1) = ceil((1 + 10^-5) 3 n_b n_g/(OWNERS^2 HEADS)) + 2^40, n_b, n_g the cumulative side counts after round j+1, and sets its s25 budget register to cap_(j+1) minus the comparisons it has done (8 instructions per owner per round, row); exhaustion takes the s25 fail exit. Tails: at every prefix the s11 argument (B_o, G_o binomial with deviation 10^-7 mean + 2^30, Hoeffding over giant blocks of at most 64 k_g comparisons with deviation 2^39) has each of its three exponents at least 143 nats (the minimum over the mean G of (G/10^7 + 2^30)^2/(3G) is 4 x 2^30/(3 x 10^7) > 143; Hoeffding above 700), asserted above 120 nats at prefix scales 2^-4k of each side; each owner-round event is below 2^-200, so 2^62 owners x 2^16 rounds stay far inside the 2^-40 union member. So the comparisons of a run that stops after round T are at most OWNERS x cap_T, deterministic in T, and with u = x y, P(T >= j) <= e^(-K_T u_(j-1)) + 0.0005 + FAIL gives E[cap_T] <= OWNERS x [(1 + 10^-5) 3 N_b N_g/(OWNERS^2 HEADS) (1/K_T + du + (0.0005 + FAIL) u_end) + 2^40 + RMAX + 1], du <= 2 MAXR/min(N_b, N_g): u-weight 0.13602, i.e. 0.136020 of the every-run comparison rows with the 2^40 + RMAX + 1 slack, comparison_rows_expected = 340348566647448052506376002380759040 instructions with the cap pair.

**Result.** E_run = W_run - (records, allowance, table, walk and comparison rows) + their expectations above = 156393908636956402415832412106483204887 (asserted below W_run); W(n) = (64 x 2^64 + K E_run)/(2000 x 2^64) <= 2^114.555827210730 (s43, s52); every opened run at most 2^117.530378341483 units.

**Evidence.** audit_il and proto above; ec (declared input) walks interleaved bases on a real toy curve and checks every trial that the run stops exactly at the first round whose prefix sets cover t; census runs the two-sided production order in the declared input at the adopted K = 7.3909 (m = 2102730, d_t = 64, N = [1750, 1484], 23 rounds: lowest prefix deficit -0.008854, -8.6e-05 after the literal-repeat term; shuffled control -0.008819, -7.9e-05; source-s7). All finite; production rests on H1'-H5.

## 46. Montgomery rows of reused operands

Prices and layouts only: points, keys, records, rounds and events stand (leaf and inverse: s47). Supersedes the s26-s28 and s32 counts, the s1/s27 first product and lane start, and the slot sizes of s8.

**Rows.** For v <= p, Row_0 = v and Row_(j+1) = (Row_j + q p)/2^64 (q = Row_j[0] n0 mod 2^64): Row_(j+1) <= (p + (2^64 - 1) p)/2^64 = p and Row_j = v 2^(-64 j) mod p. For any 384-bit a, S = sum_k a_k Row_(5-k) = a v 2^-320 mod p and one Montgomery word r = (S + q p)/2^64 is congruent to a v R^-1. Rows <= p give S < 6 x 2^64 p < 2^448 and r < 7p (RP6, 186 with the row loads); rows <= 2p take an 8th word, r < 13p < 2^385 (RP7). The correction subtracts m(idx) p, idx = r >> 376, m(idx) = floor(idx 2^376/p), from a 512-entry table: m(idx) p <= idx 2^376 <= r < (m(idx) + 1) p + 2^376, so the result lies in [0, p + 2^376), p + 2^376 < 1.04p < 2^382.

**Where operands repeat.** Every point that is ever a lane op's P (midpoints and roots, written by the full kernel; base points, strand starts and restarts) carries rows 1..5 of x at +192..431 and of y at -240..-1; gap-table entries also carry V rows of v2 = p - x2 at -480..-241 (V_j = p - X_j, rows.vrows_l, 60 per entry, executed on every table-level output). A midpoint serves all its pairs (s30).

**Batch round** (rows.py emitters, batch_audit.py listing). Slot of 496 bytes per lane op: D rows 0..287, c 288, e 336, mu- 384, mu+ 432, P and T addresses 480. Forward element, 292: store P, T; e = v2 - x1 (stored for the kernels); D_j = X_j(x1) + V_j(v2) <= 2p for j = 0..5 (D_0 = d_i of s27), formed pairwise just before row j of the product and stored; c_i = RP7(c_(i-1); D) corrected and stored, and kept in registers for the next element (two register-renamed copies alternate in the unrolled loop; the argument pointer lives on the round frame, ldr/str per element). End: C2 makes c_(n-1) canonical, s15's Fermat, one negation (s27). Backward element (875 in current s49): inv_i is the raw seven-table representative of t c_(i-1) R^-1, with p<=inv_i<8.110p<R; A = RP6(inv_i; Y rows of P), B = RP6(inv_i; Y rows of T); mu- = A + B and mu+ = A - B, each corrected by one table step and stored (the sum, below 14p < 2^385, indexes the unsigned table; the difference, in (-7p, 7p), indexes a signed table at 256 + the top byte of its 384-bit two's complement on a borrow, ranges [0, 182] and [329, 511] disjoint, entries the low 384 bits of floor((idx - 512) 2^376/p) p); t = RP7(t; D) corrected. mu- = (y1 + y2) inv_i R^-1 and mu+ = (y1 - y2) inv_i R^-1 mod p are the slopes of P - T and P + T, both below p + 2^376. t and A spill to the round frame (counted). BATCH = (32 x 292 + 1 + 32 x 875 + 3)/32 = 9337/8; the fixed part is s26's plus the two table addresses and the end slot (172597 since s47, with the frame copies of p, n0 and the table base); the two tables are 64 KiB per generating node, built once per node (row 'row-product tables (s46)', GENS x 2^15 = 2^62 instructions).

**Kernels** (kernel_check.mu_text, rows_text). The lane start (t = y1 + y2 and the first product, 309 instructions) becomes ldr buf and three ldp of mu (minus kernel at buf+384, plus kernel at buf+432); e is read at buf+336. Bounds: mu < p + 2^376 < 1.04p < 2^382 replaces s27's 1.24p, so V < 1.11p, X in (-p, 2.11p), z < 1.22p, inside the existing reduction arms (built for mu < 1.24p, z < 1.26p) (x', y' canonical); K1 holds for the remaining product (w < 2p). The full kernel stores the rows of y' (from registers) and of x' (three loads) before the next lane: 140 + 3 + 140; worst 894 over 12 paths. The leaf builds no rows (paths: s47, s49).

**Charges.** KERNEL = 894 (full children, roots, padded pairs, tables), leaf price s47, BATCH = 9337/8; gap entries + 60; restarts, base points, strand starts and re-walk rebuilds + 283 each; PAIR = 2 KERNEL + BATCH + 28 + s < 3300 keeps the re-walk step price (asserted). Point prices (pinned): baby 1145678673093/1310720000 = 874.0835, giant 1419782993093/1310720000 = 1083.2085. Memory: full-child slots 672 bytes, entries 912, batch slots 496, two 32 KiB tables per generating node (s4 total).

**Leaf.** Superseded by s47.

**Evidence.** batch_audit.run executes a 64-lane round of genuine points and a table level: every stored c, e, mu- and mu+ equals the exact big-integer model (rows.py), mu (x2 - x1) = (y1 +- y2) R mod p, and both side calls give exact P -+ T with exact rows; kernel_check runs 509 lanes per kernel with rows checked, a 672-byte footprint audit under canaries and an overlapping-stride control.

## 47. Seven-chunk table Montgomery reduction (current table trade)

For X=H*2^384+L, X*R^-1=H+L*R^-1 mod p, R=2^384. Split L at shifts0,55,110,165,220,275,330 into six 55-bit chunks and a 54-bit final chunk. T_i[c]=c*2^shift_i*R^-1 mod p, canonical for i>0; T_0 has an additional+p. The exact model is rows.tmr_val; emitters rows.tmr_l and kernel_check.tmr_leaf_block use the same chunks.

**Packed table layout.** Allocate 384*2^55=3*2^62 bytes. If allocation starts at A and base register is B=A+256, T_i[c] occupies B+384*c+48*i-256 through the next47 bytes. Each row has seven 48-byte entries followed by a 48-byte correction slot. For0<=idx<512, correction m(idx)*p mod R, m(idx)=floor(idx*2^376/p), occupies B+384*idx+80 through+127. Thus all eight slots exactly fill the 384-byte row. Every LDP offset is in[-256,112],16-byte aligned and encodable. For c=2^55-1 the final row ends exactly at A+3*2^62. The last column only indexes c<2^54, but its whole2^55-row allocation and initialization are paid.

**Addressing and registers.** MADD computes B+384*c in one scalar instruction. At the first row a temporary register receives 384. Chunk1 reads the last live bits of low product word L0; only after that extraction does L0 receive384 and remain the stride register through correction. The first two MOVs cost two instructions per reduction. Rows2–6 never read the original L0. The leaf uses x17; the backward inverse uses x0. No live product, accumulator, pointer or frame value is overwritten. The new MADD interpreter case is independently checked on input/output aliases and modular overflow; listings still pay one instruction each.

**Bounds and correction.** The slope satisfies0<=mu<p+2^376. With p/R<0.1016 and2^376/p<0.039, H=floor(mu^2/R)<0.110p. All seven table contributions total<8p, including the+p bias; with e=p-x2-x1 in(-p,p), the leaf sum is strictly positive and<9.110p<R. The inverse product has the same H bound and e=0, so its sum<8.110p<R. Section49 proves the effective high word is zero even after the signed addition, so its explicit bookkeeping is omitted. The correction index therefore fits nine bits (indeed below237); all 512 slots are supplied. By the definition of m(idx), the corrected value lies in[0,p+2^376). One conditional subtraction makes the leaf x canonical. The inverse skips that correction entirely: s49 proves the raw value belowR is sufficient for both unchanged RP6 consumers.

**Counts (current s49).** The seven table rows require84 instructions including two stride MOVs; seven zero carry-word ADCs are deleted. The leaf still uses the11-instruction correction, with LSR replacing an EXTR whose high word was zero. Three signed-e high-word bookkeeping instructions also disappear. Leaf paths are213 or214, worst214. The backward inverse uses its raw table sum and omits the11-instruction correction, so its multiplication costs230 instructions; the backward element is875. Forward292 and fixed172597 remain unchanged. BATCH=(32*292+1+32*875+3)/32=9337/8. The full-child kernel remains894. Final point prices are1145678673093/1310720000 per baby and1419782993093/1310720000 per giant.

**Construction, including the bias.** Each column has a canonical accumulator a_c=c*T_i[1] mod p. The recurrence stores its value (three STPs), advances the pointer by 384, adds its fixed step (six instructions), conditionally subtracts p (six subtracts and six CSELs), and executes SUBS/B.EQ/B control: 24 instructions per entry. Column 0 first computes a_c+p into a temporary six-word register vector, costing six additional instructions; it stores that biased value while leaving the canonical accumulator for the recurrence. Its cost is 30 per entry. table_trade_audit.build_body emits and checks both exact listings. All seven 2^55-entry columns are initialized, so the row is GENS*((30+6*24)*2^55+2^24). The per-node2^24 fixed setup bound includes every TMR step constant, correction value, pointer/counter initialization and the shared beta root. Explicitly, obtain beta from the x-coordinate ratio of [lambda]P0 and P0 (the inherited scalar-multiplication bound SM=8192000 plus one field inverse and multiplication); check its cube and the preserved y-coordinate. Compute R^-1 by384 modular halvings and successive reduction-table steps by330 modular doublings; even30 instructions per six-word half/double is below22000 instructions. The beta table step constants are separately covered by its eight2000-instruction setup entries. For all512 correction entries, increment idx*2^376, maintain q with q*p<=idx*2^376<(q+1)*p using seven-word comparisons/additions (at most one q increment because2^376<p), and store low384 bits of q*p;128 per entry bounds this loop including initialization. The scalar multiplication, field arithmetic, fewer than22000 half/double instructions,65536 correction instructions, setup loads/stores and checks together are well below2^24. Thus derived constants are constructed and charged, not free advice. No column, unused allocated table value, cold setup or future-target amortization is free.

## 48. Retained table-memory trade and provenance

The previous package introduced the beta/reduction table trade and paid all cold setup and table construction costs. This revision retains that layout and accounting. Point sets, class plans, emitted keys, round schedule, root recovery, gate and H1′/H2–H5 assumptions are inherited from40ea9361. Canonical output coordinates and keys are unchanged; internal inverse and slope representatives need not be bit-identical to the previous implementation. Their bounds and congruence are the subject of s49. The previous table construction and MADD/address audit remain in the declared run.

## 49. Zero high words and a lazy backward inverse

Let R=2^384, B=2^64, Delta=2^376 and h=(p+Delta)^2/R. The exact protected integers satisfy9p+h<R,7p<R and14p<2R; range_audit asserts these rational/integer inequalities directly.

**Table accumulation.** All seven table entries are nonnegative, and T0>=p. The high half H of any product in question satisfies0<=H<h. Therefore every partial sum after the first table lies in[p,8p+h), strictly belowR. Its seventh64-bit carry word is exactly zero. Removing its ADC after each of the seven rows neither changes the lower six words nor changes flags (ADC without S preserves flags). Each following addition starts with ADDS and each correction starts with SUBS, so no deleted flag or carry dependency survives.

**Signed leaf term.** Before e, S0 lies in[p,8p+h). Since -p<e<p, the mathematical result S0+e lies in(0,9p+h), still belowR. If e>=0, the represented addition has no outgoing carry. If e<0, its six-word two's-complement representation is R+e, so the represented sum is R+(S0+e), with carry exactly one; subtracting the sign extension would give high word1−1=0. Thus the previous LSR-sign, ADC-top and SUB-top instructions compute an identically zero value and are deleted. The six-word result is already the exact positive S0+e moduloR. Correction indexing uses LSR on its high limb instead of EXTR from a zero high word. The subsequent correction and canonical subtraction are unchanged.

**Lazy backward inverse.** The raw sum v=H+sum_i T_i belongs to[p,8p+h), belowR, and v is congruent to t*c*R^-1 modulo p. Its only arithmetic consumers are A=RP6(v;Yrows(P)) and B0=RP6(v;Yrows(T)); the running inverse t is separately saved and reloaded for its original update. For any scalar a<R, its digits satisfy0<=a_k<B, and every Y row lies in[0,p]. Thus S=sum_k a_k*Y_(5−k)<=6(B−1)p. One Montgomery word adds q*p with0<=q<B, yielding S+q*p<7Bp<BR=2^448. All seven-word intermediate accumulators therefore fit, and the output r=(S+q*p)/B is below7p and congruent to a*y*R^-1. The RP6 lemma never required a<p+Delta. The raw inverse v consequently needs no table correction.

The two RP6 outputs remain in[0,7p), their sum below14p<2R and difference in(-7p,7p). The unchanged signed/unsigned mu corrections produce representatives in[0,p+Delta). The resulting mu may differ from the previous implementation by a multiple of p, but its residue and accepted range are the same. The full and leaf kernels therefore produce the same canonical curve points and orbit keys. No success/coverage assumption or failure event changes.

**Evidence and scope.** range_audit executes116 RP6 emitter cases, including a=R−1,0,p,2p,7p,8p and boundary rows; checks216 raw inverse products, including208 above2p in the declared seed; and checks200 prefix/signed-e cases. The composed batch audit has64 raw inverses above2p and verifies its exact model, slope congruences, accepted slope ranges and canonical P±T outputs. Both leaf signs and full kernels remain in the complete declared experiment. These finite cases supplement the range proof, not replace it. The ledger, exact-log certificates, negative grid controls and manifest/hypothesis checks run on the final bytes.

## 50. Declared256-byte DC ZVA owner initialization

The complete source, instruction-model discussion, exact count, address proof, architecture excerpts and evidence for this section are in **init-proof.md**, supplied as part of the claim.

## 51. Fresh K_T/Q ledger scan at this candidate's prices

Retuned from b574037 (K_T = 14717/2000, Q = 5/4) by the s51b joint scan; full scan in **scan-proof.md**.

## 52. TMR full-child kernel, zero-spill batch round, pinned-stride leaf and aligned beta chunk

Complete scalar AArch64 derivations, exact range inequalities, regenerated ledger integers (E_run = 156393908636956402415832412106483204887, W(n) <= 2^114.555827210730, every-run <= 2^117.530378341483) and audit summary are in **s52-proof.md**, supplied as evidence.
