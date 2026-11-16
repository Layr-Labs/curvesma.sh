# BLS powers: run gate, interleaved round-stopped babies and giants, Montgomery rows (expected work 114.704949928220679)

## Complete claim and source scope

The complete claim is claims[0] of submission.json: one cold run that draws a gate coin (K = 7236538282107112258; a closed gate is a charged failure), walks its babies and giants in interleaved rounds and stops at the first verified s; success at least 0.39; W(n), the run's expected total work, at most 2^114.704949928220679 units (exact score 114704949928220679/10^15); every run at most 2^117.422738478825422 units; conditional on H1' (s45), H2-H5 of source-proof.md and H7 (s37). s2, s4 and s43-s46 derive every figure in it. The flat unsigned 64-bit node address and communication model is analytical, not hardware.

Source scope. This revises our b5b30ac6 (declared dependency): the s46 forward element forms the D rows inside its row product and keeps c in registers, and Q is re-selected at the new prices. Sections 1-44 and the code are the construction of e04c792b (aburan28, declared dependency; its lineage is credited there), carried complete and re-executed. Sections 45 and 46 are the change: babies and giants grow together in interleaved rounds, every owner record both inserts and looks up, and the run stops at the first verified s (s45); reused operands (midpoints, table entries) carry Montgomery rows, so the batch round forms both slopes by row products and the kernels start from them (s46). Where s45 or s46 supersedes an earlier section, the later text says so; no acceptance premise rests on the dependency's review.

Parameters: r as protected, d=2^20 dividing r-1, m=(r-1)/d, six-element GLV/negation classes; Q, N_b, N_g: s2, s45. No previous hosted decision proves H1', H1-H5 or H7.

## 1. Inherited arithmetic derivation; parity output superseded by s40

The kernels derive from the measured listing measure/side_lanes.S (SHA-256 b4c7a7c2d17e6decfc8c0ee93f1def949602ad9d7e2ee0edff4a7ec3c8e4be72; measurement.log: GCC 14.2.0, Valgrind 3.24.0, blst f7a88067, (Ir(8192)-Ir(4096))/4096 = 922.925 instructions/lane for the then 929-instruction listing, field multiplication 383, harness self-checks passed; no new native measurement is claimed). side_lanes(args, n): args[i] = {out, P, T, buf}, args[n] = {NULL, NULL} ends the lane loop; it writes fully reduced x' of P-T to out+0 and y' to out+48 (rows: s46). Inputs are canonical Montgomery residues (P a kernel output, root or strand start; T a gap-table entry with v2 = p - x2 at T+96, n2 = p - y2 at T+144; base-point and strand-start entries get v/n inside their scalar-mult charge, SM = 8192000/point).

Arithmetic. R = 2^384, p/R < 0.1016; the word-serial Montgomery product returns S = (ab + qp)/R < ab/R + p without its final conditional subtraction (blst's __mul_mont_384, __sqr_384, __mul_by_1_mont_384 inlined, n0 in x18, p in x5..x10). From mu (s46): V = S(mu^2); X = V + v2 - x1 (as V + e, s32) reduced to x' by a sign branch or one or two borrow-checked subtractions; w = x' + v2 < 2p; z = S(w, mu); y' = mu(x' - x2) + y2 = z - n2, reduced by one branch. Code after the x' reduction is emitted twice, so no x' path jumps through shared code.

Derivation from the measured 971 listing of 99e03e91/fd385c98 by asserted text transforms. K1, value-bound carry deletion. Lemma: in a word-serial Montgomery product with a < 2p and any 384-bit b, acc_j < a + p before round j and every pre-shift sum acc_j + a b_j + q_j p is below 2^64 (a + p) (induction from acc_0 = 0: (a + p - 1) + (2^64 - 1)(a + p) = 2^64 (a + p) - 1). With a + p < 3p < 2^383, word 7 of every partial sum and word 6 of every shifted accumulator are zero, so their carry words are deleted (18 per product) and the next round's first read of word 6 reads xzr. The kernel's remaining product has a = w < 2p; in the square mu < 2^382 puts the top word below 2^62 and the cross sum below 2^703, so doubling never carries into word 11. K2: b words load in ldp pairs with a partner register dead inside that product. nocarry_witness.py checks the deleted words at extreme and random operands.

Fixed per call: prologue 22 + epilogue 10 (asserted 32 on every call) + exit branch + 5 staging + 2 round control = 40, charged as 48 per call since s37; each round issues eight calls (four side-output kernels, four beta_keys passes), row (BASES+2) x (8 x 48 x LANE_OPS + 64), both phases, inside the instruction total. Current kernel counts and point prices: s46 (pinned).

## 2. Opened-run success floor (s43, s45)

An opened run's floor is lo_run = 1-exp(-K)-0.0005-2^-20-7*2^-40-2^-108, where K = 6 N_b N_g/m >= K_T = 10463/2000 counts the two side targets and the subtracted terms are the H1' allowance at the full sets (H1), the modelling budget, the H3/H5/cap-tail union (the former s12 member now covers the per-round comparison budgets of s45) and the diagonal event (source-s7, s11; s43). Each side emits its target plus the rest of its last class walk (at most EMAX more), which only adds ratios: H1' at the full sets gives coverage at least 1 - e^(-K_full) - 0.0005 with K_full >= K.

The floor is evaluated exactly: evidence.py encloses exp(-K/8) by the alternating Taylor series (60 terms, both partial-sum brackets) in rational arithmetic, takes eighth powers (bracket width below 10^-48) and rounds lo_run down to the 10^-60 grid: lo_run = 0.99415354528490547951... K_T = 10463/2000 and Q = 603/500 are design values chosen by a ledger scan, not minima. N_g = 190142035840427625271021522089034850 is the least integer with 6 ceil(Q N_g) N_g >= K_T m (asserted both directions), N_b = ceil(Q N_g) = 229311295223555716076851955639376030, so K = 5.231500000000000000000000000000... >= K_T. Every N-dependent quantity is regenerated from N_b + EMAX and N_g + EMAX by the same code with its asserted exponents (list: s43, s45); tails are evaluated at these values, not inherited.

## 3. Explicit pool-cap charge and allowance 1

The one-instruction-per-point allowance covers residual bookkeeping; named pool-cap checks and key/index access costs are charged explicitly.

Residual bounds, retained from the base. (i) Pool-cap checks: charged exactly (wave row of s18, both pools), residual zero. (ii) Key reloads: zero (key stored to slot, reloaded by listed ldp+ldr, travels in-record). (iii) Head index: one lsr and one add (s45), residual zero. (iv) Cap enforcement: charged exactly (the pair of s25 per comparison, per-round budgets of s45). (v) Pass control below 2^115.1: owner-init unrolled 2^20-fold (control below 2^99); all other bulk rows below 2^112 in sum (at s45's N the largest are gap-table entries BASES x GMAX = 2^109.43, planner scan visits 2^108.58 and augmentation ancestor steps 2^108.4, on top of the former 2^109 budget) with at most subs+b.ne control per level over four levels (below 8x sum, below 2^115); vs retained N_b+N_g (2^118.34). The allowance covers only purely-unknown misses.

## 4. Exact complete ledger

Unchanged source dimensions, quotas, capped schedules evaluated by evidence.py at (KT, F, ELL, KERNEL=894 and KX (s46), BETA = 128, BATCH (s46), the s28 share, N_b, N_g; s38 tables; s45 owners). Every row below is an exact integer instruction count, rounding at original locations. The every-run total of an opened run (both sides walked to their targets, every loop at its hard cap) is W_run = 445447422996986209108288716143471760003 instructions; the bound on its mean is E_run = 172603434747228806072920362889798020392 (s45); the 64-instruction gate is added in s43. Challenge-independent preprocessing (sieves, planner, plan copies, beta tables, owner initialisation, windows) is below 2^107.1114 units, included in both; babies are no longer preprocessing (s45). Transport per record GEN 12 + RELAY 770/64 + owner body 770/64 (s45).

- records: N_b + N_g points (every-run): 441385983949647917112351201368427942282.
- comparisons at the per-owner cap (every-run): 1377745271740761511439821508215570432.
- gap tables (every-run): 1398363181441903306596747108892572792.
- listing allowance (1/point, every-run): 419453331063983341636103853880122624.
- comparison cap enforcement (every-run): 393641506211646146125663288061591552.
- owner initialisation: 332306998946229263373856944422912000.
- the other 30 rows, together 139928757933768426765322071571048321 (0.081% of E_run, each below 0.05%): names and values in the experiment output (ledger rows_instructions).

In E_run the records, allowance, table, walk and comparison rows are replaced by their expected values of s45 (expected_baby_points, expected_giant_points, expected_bases, walk_rows_expected, comparison_rows_expected in the output); every other row keeps its every-run value.

Memory and all per-node layout assertions are retained under the s8 slot size, owners per s45: total below 2^125.0002 bytes; the generating node and the owner below 2^64 bytes, every other node below 2^63 (s38). No free precomputation, future-target amortization, or wall-clock parallel discount.

Score grid is 10^-15: evidence.py reports the least fifteen-decimal ceiling 114704949928220679/10^15 of the exact log2 of E_num/(2000 x 2^64), E_num = 64 x 2^64 + K E_run (s43, s45; 90-digit Decimal enclosure). log_bound.py encloses ln(x) by 64 rational terms of 2*sum(t^(2j+1)/(2j+1)), t=(x-1)/(x+1), plus the positive geometric tail; similarly ln2; divides in the safe direction. Applied to E_num/(2000 x 2^64) it proves the declared ceiling and rejects both one-step shifts; applied to (W_run + 64)/2000 it proves the every-run ceiling 117.422738478825422.

## 5. Success and honest limitations

Points, stored keys, routing, candidate decisions, recovered scalars match the algorithm. No section s1-s37 changes a source success event: each changes prices, layouts, bit-identical listings or superset failure events in the same union, as each section states; current full kernels preserve exact x',y'; beta passes produce the live routed keys. Parity output is removed by s40. s30's x-keys with three giant lookups test exactly the U-orbit match the parity key tested (orbit lemma in s30); leaf kernels and the s30 planner change prices, not points. N changes only through s2. H1-H5 keep source-s7 definitions, evidence, budgets; H7 (s37) is a new disclosed hypothesis; H6 is superseded by the s30 per-class floors (H4 as f_c >= 998/1000 and H4L as leaf share >= 193/200 (s36) per production class, runtime-checked before any group op; a failed check stops the run, challenge-independent). The allowance change is exact accounting with stated residuals. Those hypotheses remain extrapolations requiring review.

An opened run's floor at the realized K (s2) is lo_run = 0.99415354528490547951...; with the gate (s43) success is at least (K/2^64) lo_run = 0.39000000000000000000361... >= 0.39 (exact rationals; K - 1 fails, asserted); declared success is 0.39. All failures and capped work remain in the complete ledger.

Mode all of combined_lazy.py adds the planner controls, reduced-order curve recovery with interleaved stopping, transport with overflow controls and the two-sided round-order census to the ledger-mode audits. Finite cases only; no full-size attack run.

## 6. Superseded pooled floor

H6 is withdrawn. Use the runtime-checked per-class floors of s30.

## 8. Production lane-stack layout (128-byte slots; 192 since s30)

Supersedes the source's 96-byte lane-stack stride. Slots are 192 bytes (s30; 128 at k7), contiguous from a 16-aligned base: x [0,48), y [48,96), keys [96,160), pad; every listed ldp/stp keeps form and offsets, slot addresses use one hoisted madd. The s40 full kernel writes [out,out+96), leaving the key region to beta_keys; gap-table entries are 192 bytes (x, y, v, n; s1), read-only to kernels; batch buffers are separate. Since s46 full-child slots are 672 bytes (rows around the 192-byte slot), entries 912 and batch slots 496; memory term LANES x (2 DEPTH + 3) x 672. kernel_check.py's stride audit (s1): 8 out slots in one contiguous buffer under random canaries, outputs exact, write footprint exactly the out ranges plus the frame; since s46 the audit uses 672-byte slots with rows checked, and an overlapping 640-byte stride trips it.

## 9. Record-path re-derivation: pointer buffers, honest matching, explicit setup

This section re-derives every per-record instruction from scratch, superseding the source s5.2 transport counts. The charged record sequences (GEN12 s10, RELAY s9, owner body s45, CMP7 s25, loops s29) are asserted by evidence.py audit_record_path (lengths, bits, immediates, pool fit, division proofs); the proto toy executes the pointer mechanism on every record buffer (pre-write check, +40 advance, exact derivation, reset).

Pointer buffers. Sender per-link state becomes (ptr, end): each record's three stores post-index ptr by 16+16+8 = 40 bytes, deleting slot-address madd and counter increment. Entries stay 16 bytes (one ldp; source "24-byte" prose counted a hoisted cap word), addressed by one add-with-shift (t below 2^30, u below 2^32, lsl #4). The cap check precedes the write (cmp ptr, end; b.hs fail): identical failure event, never an out-of-bounds write. GEN listing: s10. The relay output side matches GEN's; relay input loads unchanged; s24 pairs the relay loop. At wave end the sender derives count = (ptr-base)/40 (sub, lsr #3, umull/lsr /5 with magic 0xCCCCCCCD, exact by m*5 = 2^34+1 plus asserted bound, count store, ptr reset: 8), writes the per-parity count word, and resets; receiver side (3) and double-buffer discipline unchanged. LINK_WAVE is 16 (8 sender + 3 receiver + 5 margin); the per-node term rises to NODES*24 for dmb, counter traffic and constants (16 used). Chain-id words are segment-invariant, formed once per segment by the formation row (register-holding superseded by s10; GEN reloads them).

Honest matching: this section's record key, pool layout and owner/compare/cap counts are superseded by s11, s17, s19, s25, s29 and s45 (current bodies in evidence.py, executed by audit_il). Four new explicit rows: window-table init (8/link, 4 used), chain-id formation (3/segment), candidate-link wave counts (9 on 4 traversals: /40 magic, /64 shifts, /56 with proven /7 magic), candidate formation (16/match, 11 used). Caps, routing, counts, events unchanged, so H1-H6, N, success, memory untouched.

## 10. k9 repair: chain-word memory schedule

Two chain-id words formed once per segment into a per-lane 16-byte slot (LANES slots per generating node, in memory term); every GEN reloads them. Formation 4 per segment (ubfx+bfi+lsr plus stp); GEN11 becomes GEN12 (ubfx t, entry address, ldp entry, ldp chain words, cmp, b.hs, three stores, ptr store-back, key-reload/loop pair). No chain word crosses records in a register (except within one giant point's three consecutive GEN copies, which no call separates; s34), so lane interleaves and side-output clobbering are compatible; allowance untouched.

## 11. Comparison cap per owner

Superseded by s45's per-round budgets. What stands: a run prefix's comparisons at an owner are sum over its head slots of (baby records) x (giant keys) in that slot, each pair compared once (s45); chains above 64 have probability below 2^-316 per slot at the giant-table load 3 N_g/2^119 < 0.87, below 2^-197 over all 2^119 slots (inside the union); a giant point's keys at one owner form one block (k_g >= 2 has probability 3 x 2^-124 per point, mean below 0.03 per owner); B_o and G_o are binomial with Chernoff tails, and Hoeffding over blocks bounds the sum (s45 asserts the exponents at every prefix scale).

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

Audit (every ledger run): audit_packed_match evaluates this compare word for word on the key 2^71 (head 1, f9 = 0, x10 = 512) and 3000 further keys with junk in w4[0..54]: equal to themselves; any flip in bits 62..70 unequal; a flip elsewhere leaves the field compare equal; a flip in the head bits changes the slot.

## 18. Wave-granular pool-cap check (two pools since s45)

One check per owner per wave and per pool. An owner inserts only through its 2^20 inbound relay links (source-s4.3), each carrying at most CAP2 = 2^21 records per wave, enforced by the relay's charged pre-write cap check (overflow is a run failure in the same budget). So an owner inserts at most 2^41 records per wave into either pool (asserted); checking each pool pointer against its threshold at each wave start guarantees no wave overflows it. Baby pool 2^56 records at [2^61, 2^62), threshold pointer 2^62 - 2^46; giant pool 2^57 records at [2^62, 2^63), threshold 2^63 - 2^46 (each one movz + movk, asserted). Means N_b/2^62 < 2^55.47 and 3 N_g/2^62 < 2^56.78 (giant keys, three per point); Chernoff exponents above the 82-nat budget (asserted). Per owner per wave and pool: ldr, ldr, subs, b.hi plus str at wave end (5); once per owner and pool 6. Row = 2 x OWNERS x (5 x WAVES + 6), WAVES = 11214912 (s45); in every-run value.

## 19. Full-address chain links

This section restates the link word (superseding the s7 56-bit-link sentence): head entries and chain links store full 64-bit record addresses, 0 = end or empty; records load as [addr,#imm] (comparison CMP7 since s25, owner body and head store in s45). Layout: pools after the head array of 16-byte head pairs (s45), records at pool base + 32 and above, addresses below 2^63 (asserted). Pool use: s42; word sizes, memory, owner init unchanged. Links never leave the owner (chain ids are base+segment).

## 20. Itemized side-call fixed price: 40 per invocation

This section replaces the s1 64/invocation bound with a complete per-call schedule: 32 fixed (prologue 22 + epilogue 10 with ret, asserted every call) + 1 exit branch (measured) + 5 staging (mov x0,x_args; mov x1,x_n; lsl x2,x_n,#5; stp xzr,xzr,[x_args,x2]; bl) + 2 round control (subs, b.ne) = 40 (asserted). The per-lane argument records sit in PAIR_BK/ROOT_BK (s26), so staging is exactly these 5. Zero spills: 8 loop variables in x19-x28 across the call; the kernel preserves callee-saved regs (asserted); nothing live in caller-saved regs across bl. Lanes are unrolled x n (n fixed; PADOPS pads to full calls), so no per-lane loop control exists; the 40 is n-independent anyway. Driver and kernel share one small text (bl in range, no veneer). Per-base-per-phase setup bounded by 32 (entry/exit 10 + materialize 8 + base-loop 2 + 12 margin; true cost <= 20), charged as 64 x (BASES+2). Row: s1 (eight calls per round since s30; 48 per call since s37), both phases charged.

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

Per comparison: CMP 7 (ldp, ldr, and, cmp, ccmp, b.ne, cbz or cbnz) plus the 2-instruction cap pair, every-run rows OWNERS x cap x 7 = 1377745271740761511439821508215570432 and OWNERS x cap x 2 (cap 42678697621513294 per owner, s45), charged in E_run at their s45 expectation. Pool bases are 16-aligned and records are 32 bytes, so every ldp is aligned. No signed condition; addresses stay below 2^63.

Executed: audit_il (s45), every ledger run.

## 26. Batch round: both children by the derived kernels, no reference addition in the walk

Points, keys, records and events are unchanged. Round (one lockstep step of a node's 2^19 lanes, s42): the driver writes argsB[i] = {out, P, T, buf+496i} per lane op (out = the P+T child, or the root's point) and argsA = {P-T out, P, T, buf+496i} per pair (since s30 split into full and leaf arrays by the leaf flags), each ended by a 32-byte NULL record; padded lanes use the fixed pair (Q_a, T[2]), x distinct, so every denominator is nonzero. batch_round(argsB, buf) runs Montgomery's trick over d_i = x1 + v2 (s27) and leaves both slopes per lane op (s46); side_lanes(argsA) writes P-T and side_lanes_plus(argsB) writes P+T (side_lanes on -T = (x2, p-y2, v2, y2): its n2 loads read T+48/64/80). Counts: s46. PAIR = 2 KERNEL + BATCH + 28 + s for two points (plan entry 8, P address 1, T[g] address 6, buffer 1, two out addresses 2, two argument records 4, s30 routing 6); ROOT = KERNEL + BATCH + 21 + s (gap entry 5, T[gap] address 6, P and out addresses 2, buffer 1, argument record 2, loop 5). Side-call row: s1. Re-walks keep 3300 + 232864. The reference addition (2000) is not used anywhere in the walk; every product, the inversion and both kernels are charged at their listed and executed counts, nothing twice.

## 27. Lazy batch round: reflected difference and unreduced inverses

The forward pass forms d_i = x1 + v2 (v2 = p - x2 at T+96) in (0, 2p), d_i = -(x2 - x1) mod p, nonzero as x1 != x2. Montgomery's trick is linear in its running value, so the round negates the canonical Fermat output once (t = p - t) and every backward inverse is inv_i = 1/(x2 - x1); inv_i stays unreduced (bounds: s46).

## 28. Exact round share and unrolled batch loops

A walk round is one lockstep step of a node's 2^19 lanes, each doing one pair or root operation, idle lanes padded and charged (source s2 step 4, s3 padding lemma; s20): exactly 2^19 lane ops and one fixed cost 172591 + 5 = 172596 (s26, s46), so each lane op carries s = 172596/2^19. Every executed lane op sits in a row charged at PAIR or ROOT (records at the floor F; padding; remainders), so the shares cover every round. The forward loop's only per-element control is cbnz x0, the backward loop's cmp x2,x30 + b.ne, each over 32 copies of its body; n = 2^19 is a multiple of 32 (asserted). The listing (about 2.2 x 10^5 instructions with the Fermat chain) sits inside gen_mem's 2^20-byte code term.

## 29. Transport loops 64-fold (s24)

Candidate delivery is eight relay hops at 770/64; evidence: audit_il (s45).

## 30. x-only leaves, beta-triple giant lookups, leaf planner and asymmetric sides

Supersedes the record key (source s1 parity key), source s3 matching, the s16 floor, GMAX 2^44 -> 2^47 and DEPTH 2^14 -> 2^17 (source s2-s3), N of s2, the 128-byte slot of s8, the per-point price of s26-s29, the two side calls per round of s20/s26 and the 2^34-step giant wave (source s4.3). Everything else stands.

Keys. Every record key is limbs 1..3 of the canonical Montgomery residue of x' (out+8..31). [-1](x, y) = (x, -y) and [lambda](x, y) = (beta x, y), beta fixed by [lambda]P0 = (beta x0, y0) (asserted in beta_run), so the U-orbit of Q has x-set {x, beta x, beta^2 x}; conversely x_P in that set gives y_P^2 = x_P^3 + 4 = y_Q^2, so P is in U Q. A baby x-key equal to one of a giant's three keys is exactly the orbit match the parity key tested, and at most one of the three can match (x = beta x forces x = 0, a 3-torsion point; 3 does not divide r). Babies insert one record per point; giants send three lookup records per point (keys of x', beta x', beta^2 x', one chain id; GEN unrolled three times per point at key offsets +8/+104/+136; a giant re-walk keeps positions where any of its three keys equals the candidate key). Coverage events, K, success and the 96 recovery tests are unchanged. Giant waves are a quarter of a pair wave (2^31 steps since s42), so a generating-node link still carries at most 2^23 records per wave in the mean (six per step).

beta_keys (batch_audit.beta_run; tables since s38): per giant slot, limbs 1..3 of canonical beta X' at +104 and of beta^2 X' = -(X' + beta X') at +136, walking the side-call argument arrays (the driver stores nothing extra).

Leaves. A pair child that is no pair's midpoint is only keyed: it runs the leaf kernel kernel_check.leaf_text(), the kernel cut after the x' store (y' tail and mu spill deleted; prices since s32, s37), bit-identical in x', executed by run_leaf (365 lanes each copy). Midpoints and roots run the full kernel. Slots are 192 bytes: x' 0, y' 48, keys 96..159.

Driver and calls. Plan entries carry two leaf flags (kind 2 + gap 48 + three 19-bit slots + 2 = 109 of 128 bits); per pair lane op two ubfx and two pointer ldr/str route the children to the full or leaf call (ROUTE = 6, listed). Per round four side calls (full/leaf x minus/plus) and four beta_keys passes, charged 8 x 40 on both phases (s1, s20). Pads run full pairs plus two beta passes.

Planner (planner.py). One centre per run of 50 class-list elements (the one with most candidate pairs); greedy by increasing degree over centre pairs, hosting centres first; the s16 augmentation; a fallback pass over all midpoints; source s3 demotion and emission; every pair passes the shipped cycle, depth and tree rules with exact rollback. GMAX = 2^47 = 128 d^2 gives each midpoint about 128 candidate pairs; WIN = 2^29 keeps rel32 offsets at most 2^30 < 2^31. Charged: 6 window scans per element (was 3); matching at 64 kappa = 512 candidate pairs per element (twice the mean 2 GMAX/d^2), at most five tries each (centre pass twice, fallback three times); two augmentation passes, each with its partner cache; memory adds a 3-byte centre flag and count per element (planning node 2^62.79 < 2^63). Each class is planned before any group operation and the run stops unless f_c >= F = 998/1000 and the leaf share ell_c (leaf children / pair children) >= ELL = 193/200 since s36 (H4 restated; H4L, a cost floor). Measured at kappa = GMAX/(16 d^2) = 8, span X2/(2 GMAX) = 4: the declared payload's spec at d = 4 (mode planner, offline here); d >= 16 by a line-for-line C port of plan() that reproduces planner.py's forests exactly (pairs, roots, leaves) on six declared-type lists, capped and uncapped:

| d | classes | min F | min ell | mean ell | max height |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 4 | 2 | 0.99550 | 0.97349 | 0.97375 | (mode planner) |
| 16 | 4 | 0.99822 | 0.97278 | 0.97298 | 235 |
| 64 | 4 | 0.99863 | 0.97213 | 0.97296 | 931 |
| 128 | 4 | 0.99875 | 0.97305 | 0.97329 | 1062 |
| 256 | 8 | 0.99877 | 0.97306 | 0.97324 | 2172 |
| 512 | 4 | 0.99886 | 0.97311 | 0.97338 | 6931 |

ell is flat in d (0.972-0.974) and F rises with d and span; F = 998/1000 sits 0.00077 below every d >= 256 row, ELL (96/100 at k36; 193/200 since s36) below every row. Tree cap: at TREEMAX = d (the production cap/window ratio 2^-9) the leaf deficit is about 5/d (4.5-5.1 at d = 32..128) and the F excess 1.1/d to 2.2/d (d = 16..128); at fixed d = 64 the F excess falls 0.0193, 0.0100, 0.0027, 0.0011 for caps 64..4096. Both routes put the production cap effect below 10^-4 at TREEMAX = 2^20, inside the margins; that extrapolation is part of H4/H4L. Heights grow with tree size (6931 at 2^19 elements): DEPTH = 2^17, stack 2^18 + 1 slots.

Sides. With side prices c_b, c_g the optimum is N_b/N_g = c_g/c_b (s45: Q = 603/500). The ec toy (x-keys, triple keys, interleaved bases, solved == exact coverage every trial, stop round == first covering round, all scalars correct): the declared input (mode declared, K = 5.2315, Q = 59/50) gives 4/4 at 0.9948 predicted, 4 trials stopping before their last round. Its rate is a binomial sample of toy coverage: it checks that success is exactly the coverage event; production success rests on H1 and the census (s2, source s7), not on toy rates. Caps with 3 N_g giant keys and N_b baby keys: s45. H3 now covers N_b + 3 N_g keys: a giant's three keys are images of x' under a fixed bijection, treated as three uniform keys; the toy's candidates follow the false-match model plus true hits; the re-walk step price 3300 covers 2 x 929 + 2 x 343 and key compares.

Prices (pinned in audit_record_path; current values s45): pair child (BATCH + PAIR_BK + 6)/2 + ELL x KXC + (1 - ELL) x KERNEL; baby point F x child + (1 - F) x ROOT + 12 + 2 x 770/64; giant point F x child + (1 - F) x ROOT + 3 x (12 + 2 x 770/64) - 2 + BETA.

## 32. Shared x-sum e = v2 - x1, once per lane op

The forward element stores e = v2 - x1 (at buf+336 since s46); kernel_check.esum_text replaces one asserted block of the measured listing so the kernel forms X = V + e. As integers e is in (-p, p), X in (-p, 2.11p) (s46), and bit 383 of its two's complement is the sign (tbz x16,#63); x', y' are bit-identical. Counts: s46.

## 33. Gap tables by pair levels

Changes only how tables T[k] = [k]Q_a, k <= GMAX, are built (entries, layout, memory and use unchanged). With T[1..M] known (M = 3 from the strand-start scalar multiplications, row unchanged), one level computes the centre c = M + k + 1 as a root T[M] + T[k+1] and then the pairs T[c -+ j] = T[c] -+ T[j], j = 1..k, k <= M - 2: the s26 lane op on (T[c], T[j]) with both children full (the s46 full kernel with rows, minus and plus calls, then the 21-instruction v/n step and the 60-instruction V rows), so T[1..M + 2k + 1] is known. Full levels (k = M - 2) give M' = 3M - 3; the last stops at GMAX: 30 levels. Incomplete additions are safe: c - j, c + j and c are distinct, nonzero and below GMAX < r, and x(T[c]) != x(T[j]) as c != +-j.

Price per base: every level is charged ceil(k/LANES) + 1 full rounds of LANES lane ops at PAIR (+1 for the centre; padding included) plus 128 per round for calls, and 21 per entry: BASES x (rounds x (LANES x PAIR + 128) + (21 + 60) GMAX) (strands: BASES x (GMAX - 1) x 2050; re-walks keep 2050). Row in s4.

Evidence: batch_audit.run(table) executes a 64-lane level (M = 66, centre [131]Q, gaps [j]Q): the s32 round and both kernel calls give exactly [131 -+ j]Q; evidence.py asserts the 30-level recursion.

## 34. Giant GEN chain-word reuse; head zeroing by dc zva

Prices only. (i) A giant point's three lookup records share one chain id (s30; GEN unrolled three times per point, key offsets +8/+104/+136). Copies 2 and 3 drop their chain-word ldp (copy 1's registers still hold the words; GEN only reads them). Giant transport 3 x (12 + 2 x 770/64) - 2 (s45); the proto's generic per-record GEN 12 stays an upper bound. (ii) Owner initialisation zeroes each 64-byte-aligned head array by dc zva (64-byte blocks, DCZID_EL0.BS = 4, common AArch64) plus one add: 2 instructions per 64 bytes, not one stp per 16. Row OWNERS x (HEADS x 16/64 x 2 + 64) with s45's 16-byte head pairs.

## 35. EMIN = 2^56 - 2^52 (H2 restated)

H2 has class lists between EMIN = 2^56 - 2^52 and EMAX = 2^57, sieve-checked before any chain starts (a short class is an H2 failure event); bases are charged at N/EMIN per side. Bases at the bound: 2^61.56 baby and 2^61.29 giant, 2^62.43 in all (I's 2^64.28 primes: margin above 3.58; s45); segments per class at least 2^33.9. Support: class-length spread shrinks like (X/d)^(-1/2); h2_sizes (mode all) on toy primes gives min/max of |L_c| d/X from 0.9287/1.1035 (d = 16-256, X/d = 2^10) to 0.9908/1.012 (X/d = 2^14), which extrapolates far inside 6.25% at X/d = 2^56 (heuristic).

## 36. Leaf floor ELL = 193/200

H4L's runtime-checked floor is 193/200 (a class below it stops the run before any group operation). Per-class minima: s30's table (0.97213-0.97349); offline planner.py rows (kappa 8, four classes, d = 4-16, span 4-8) 0.97274-0.97524; one small two-class d = 8 run had 0.9691. Prices only.

## 37. Leaf reduction at a capped share; H7

The leaf's X in (-p, 2.16p) (2.11p since s46) needs the second subtraction (+7) only if X >= p; the leaf counts those lanes (ldr/add/str; paths since s46 275 to 285) in its free slot [sp,#96], zeroed in the prologue, returned in x0 (fixed 34). The driver adds each return to a node count C and its lanes to L (inside the 48 now charged per call, 40 + 8). At every wave barrier the tree sums (C, L) (16 per node per wave) and the run stops unless C <= 7L/32. Between barriers C grows by at most one wave of leaf lanes, so a leaf child costs at most 275 + 10 x 7/32 = 277.1875 (s46) plus the row GENS x WSTEP x LANES x 2 x 10 + WAVES x GENS x 16.

**H7 (new, in the claim's conditioning).** At every wave barrier, the cumulative share of leaf lanes with X >= p is at most 7/32. A violating run stops (excluded by H7, as H2's sieve checks). Support (s46 round): batch_audit.xge_run (declared; one base Q, gap table [j]Q, walked centres, the exact s46 round model of rows, corrections and the mu combination, 16 rounds of 256 lanes, both children) gives 0.1848, per-round shares 0.1445-0.2285 at 512 flags; a barrier share averages >= 2^99 lanes. run_leaf (plain and plus): the returned count equals the exact X >= p count (66 of 365), body exactly 275 n + 10 count.

## 38. Seven-chunk beta tables (k37l)

Replaces s31's tables, lookup and evidence. Canonical Montgomery X < p < 2^381 splits at bits 0, 55, 110, 165, 219, 273, 327 into three 55-bit and four 54-bit chunks c_i; T_i[c] = beta c 2^shift_i mod p, canonical. Tables 0-2 interleave at A + 256c + 64i, tables 3-6 at B + 256c + 64(i - 3), B = A + 2^63 (one add-lsl per address); q p (q < 8) at A - 512 + 64q; production A = 512. BUILD is s31's 24-instruction recurrence at stride 256: row GENS x (5 x 2^55 x 24 + 7 x 2000 + 64) in preprocessing; storage 3 x 2^62 + 512 bytes.

Lookup, 128 per giant point: 3 loads, 12 chunk, 7 address, 21 table loads, 36 sum, 23 reduction, 2 stores, 22 for beta^2 X, 2 loop; fixed call part 26 of the 48 charged (s37). S < 7p and q = floor(s5/(p5 + 1)) <= 6 give q p <= S < (q + 1)(p + 2^320) <= q p + 2p (7 x 2^320 <= p): one conditional subtraction; keys and offsets unchanged.

Addresses: each node is one flat 64-bit space (source s4.1). Tables fill [0, 3 x 2^62 + 512), the other generating-node regions (windows included) follow, total 17609118817070351136 < 2^64 (asserted); owners stay < 2^64 (s43), other nodes < 2^63. Addresses are unsigned and never wrap: the ledger asserts that no shipped listing or mnemonic list has a signed condition, sign-extending load or asr, and that bit 63 is tested only on data (the kernel's X sign, the restart scan's scalar words); caps compare with hs; rel32 offsets add modulo 2^64 within one region. Analytical, not hardware.

Evidence: beta_run and batch_audit.run (finite checks; experiment output). Current prices: s4/s40.

## 39. Inherited evidence repairs

Planner fields: 28 bytes base, 8 adjacency, 4 previous link, 3 centre (43 bytes/element, DEPTH = 2^17). The character check certifies order d for 7^m and exhibits k = 256 where log_g(x^m) is the wrong character base; source s2 uses log_(g^m).

## 40. Delete the unused full-point parity key

Section 30 replaced every record key by x-coordinate limbs. kernel_check.drop_parity_text removes the four exit blocks that only emitted the old parity key (each asserted as tbz, subs, sbcs, sbcs, sbc, stp, str); coordinates are stored before them, no consumer reads +96..191 of a full output, and flags and scratch registers are dead at the next lane. Counts: s46.

## 41. Feasible ratio Q, not a minimum

Q = 603/500 (s45): N_g is the least integer with 6*ceil(Q*N_g)*N_g >= K_T*m, and N_b = ceil(Q*N_g). No search over all positive integers N_b, N_g is claimed. Kernels: s46.

## 42. Lanes 2^19, segments 2^22, wave 2^33

Heads: s17, s45. Round share 172596/2^19 (s46). Point prices: s46. SEG = 2^22 gives 2^35 segments on an EMAX class. Segment field 35 bits (29+35=64). Pool 64+64+9+63+35=235 of 256 (base index 63 bits since s43). Chain id 63+35=98 of 128. Directory charged at EMAX/SEG+2.

Wave 2^33 steps. Two records per step on 2^19 lanes is 2^53 per generating node, mean 2^23 over 2^30 links, half of CAP1=2^24. CAP2 and the 2^41 owner wave stay as in s18. Asserted: 2*WSTEP*LANES = 2^23*2^30.



## 43. Run gate; sizes at K_T = 10463/2000

**One run.** The attack is one cold run on one target. It first draws one fresh uniform 64-bit algorithm coin u (one RNDR read at one coordinator node, before any other step of source-s2; all cold setup follows). If u >= K = 7236538282107112258 the run stops there and outputs failure, at charged cost GATE = 64 instructions. Otherwise it releases the first global barrier (one of WAVES' 64 spare barriers, charged in the barrier rows), which starts every other node, and runs to the first verified s (s44, s45) or to the end of its last round. A closed gate is a failure, never a success. Apart from u the construction is deterministic.

**Success.** The opened construction's success event S is a function of the challenge t alone and u is independent of t, so P(success) = P(u < K) P(S) >= (K/2^64) lo_run, lo_run = 0.99415354528490547951... (s2). evidence.gate asserts that K satisfies (K/2^64) lo_run >= 0.39 and K - 1 does not; (K/2^64) lo_run = 0.39000000000000000000361...

**W(n).** W(n) is the expected total work of this one run over the probability space the rules fix (uniform challenge, fresh algorithm coins), every trial, failure and closed gate included. Pointwise work(u, t) = 64 + 1[u < K] work_run(t), so W(n) = 64 + (K/2^64) E[work_run] <= 64 + (K/2^64) E_run = E_num/2^64 instructions, with E_run = 172603434747228806072920362889798020392 the bound of s45 on an opened run's mean work and E_num = 1249051363171498192194173626250470135718421938389528468560; in units W(n) <= 2^114.704949928220679 for 2^20 <= n < r. What is bounded is this mean. It is not a bound on one run: an opened run does at most W_run + 64 instructions, 2^117.422738478825422 units (s4: never stopped, every loop at its hard cap, every overflow a charged failure), certified separately (opened_run_log_bound), and that maximum is above W(n).

**Sizes re-derived at K_T = 10463/2000, Q = 603/500** (asserted in evidence.py ledger; side counts at N + EMAX): owners, pools, heads and comparison budgets: s45. Bases 2^62.43 in all, below 2^63 (the base index is 63 bits; s42), at most 2^42 per character and side (H2's cap, sieve-checked), at most 2^15 per side per generating node; I's 2^64.28 primes: margin above 3. False matches 2^44.25 < 2^45 - 2^20 (candidate cap 2^45), Chebyshev keeps P(> 2^45) < 2^-40; owner candidate mean 2^-17.75, candidate-link mean 2^-37.7. Diagonal < 2^-108. H5 conditional mean K/coverage < 5.27. WAVES = 11214912 < 2^24. Memory < 2^125.0002 bytes.

## 44. Round-stopped run

Source-s2 steps 4-6 run per round, one candidate round (s4.5) per round; kernels, records, keys, routing and planner are unchanged. Round kinds and owner tables: s45.

**Rounds.** A generating node walks its bases one after another (source s2); round j is, for every generating node, its next base of the round's kind (s45), a union of complete class walks, fixed before the run and independent of the challenge. At most RMAX = 2^16 rounds (2^15 per side); a node without such a base waits. Round j, in this order: (i) at most 161 emission waves of the unchanged wave protocol (source s4.3); (ii) 2 drain waves, after which every round-j record has been inserted and looked up at its owner; (iii) the status pass and the candidate and reply rounds of source s4.5 on the candidates appended in round j; (iv) a close barrier at which the found flag (set only by an owner whose scalar test passed and whose recovered s verifies against A_d and T, source s1) and the fail flag (any cap overflow or failed check) are OR-ed up the barrier tree and released down, with the next round's comparison budget (s45); every node halts if either flag is set. Round j+1's emission starts only after round j's close release. So nothing is in flight at a stop, nothing is emitted after it, and only a verified s or a failure stops the run, never a false candidate. Candidate state is cumulative: owner lists (2^16 slots per run) are append-only with a delivered cursor, the status pass checks the running total against 2^45, and a link carries per round at most its whole-run load, so every cap event is contained in the whole-run event and the failure union of s2 is unchanged.

**Charges, at RMAX rounds, never scaled.** WAVES = RMAX x 171 (161 + 2 + 8 per round) + plan rounds + 64 = 11214912 < 2^24 (asserted); the per-wave rows (link words, barriers, s18, s24) grow with it. Candidate-link wave counts are charged every round (RMAX x 4 traversals). Round control is the listing ROUND_CTL (batch_audit.py), executed on toy trees by audit_round_ctl in every ledger run (no flag, a found flag, found and fail flags at two leaves, a 1024-leaf tree, the round cursor at RMAX - 2 and RMAX - 1); a node costs at most 8 + 10 instructions per round, the row charges 24 for every leaf and barrier node in every round. Per-candidate rows stay at the whole-run 2^45 cap.

**Stopping bound, through found.** Let B_<j, G_<j be the baby and giant points of rounds before j, R_(j-1) = U.{b/g : b in B_<j, g in G_<j} and T the last round run. If t is in R_(j-1) and no failure event of s2 occurs, some round before j holds a true pair; a true pair always yields a candidate in the later of its two rounds (exact orbit key; s17, s25, s45), that round resolves it, s verifies and the run stops before round j. Hence T >= j implies that t is not in R_(j-1) or that a failure event occurred: P(T >= j) <= 1 - |R_(j-1)|/m + FAIL, with FAIL = 2^-20 + 7 x 2^-40 + 2^-108 the failure terms of s2. t is uniform on H and R_(j-1) is challenge independent; no coin, ordering or independence assumption enters.

**Tables and walk rows.** A node computes a base's Q_a (baby) or Q_a' = [a']A_d (giant), strand starts and gap table only after the previous round's close release, inside that base's round (step 4 or 5, per round; a node holds one table, gen_mem's single GMAX x 912 term), so a run stopped after round T builds none for later rounds. The build, TROUNDS = 134217779 lane rounds plus one strand start per lane, fits in the round's emission slack 161 x 2^31 - LANE_OPS = 2142240768 steps. Each base walks a whole class of at least EMIN points (H2), so the mean number of bases walked is at most (E_b + E_g)/EMIN + 2 = BE = 2406608190315584774 (expected_bases), not BASES = 6209156612039369554; padding, side-output fixed calls, restarts and chain ids are per base or segment and take BE and SG = ceil((E_b + E_g)/SEG) + BE x LANES segments (walk_rows_expected), asserted below their every-run sum.

## 45. Interleaved rounds: babies and giants grow together

Supersedes the baby-first order (source s2, s4.3), the giant-only rounds of s44, the owner layout, bodies and caps (s11, s12, s17, s18, s21-s23, source s4.2) and Q, N (s2, s41, s43). Kernels (prices: s46), planner, records, transport hops, beta tables, gate, K_T and lo_run are unchanged.

**Round kinds.** Bases are dealt as in source s2 step 3, at most 2^42 baby and 2^42 giant bases per character (H2, sieve-checked), so at most 2^15 of each per generating node. With y_j = |B_<=j|/N_b and x_j = |G_<=j|/N_g the cumulative shares after round j (exact class lengths from the sieves, before the run), round j+1 is a baby round if y_j <= x_j and babies remain, else a giant round. The round table (kind, n_b(j), n_g(j), cap_j) is accumulated during base assignment (one add per base, inside its 64-instruction charge), held at the barrier root and sent down in the close-release word (the owner's load is inside its 8 per round). In a baby round every generating node walks its next baby base; in a giant round its next giant base.

**Owners.** 2^62 owners (t, u) = key bits 0..29 and 30..61; per owner 2^57 head pairs of 16 bytes (baby head at +0, giant head at +8; [0, 2^61)), indexed by key bits 71..127 (s17); a baby pool of 2^56 and a giant pool of 2^57 32-byte records ([2^61, 2^62) and [2^62, 2^63)), records as s22. A baby record inserts its key into the baby table and looks up the giant table; each of a giant point's three records (keys of x', beta x', beta^2 x', s30) inserts into the giant table and looks up the baby table. One body per kind (evidence.il_body; loops unrolled as s24), 12 instructions:

    ldp x0, x1, [x20], #16 ; ldp x2, x3, [x20], #16 ; ldr x4, [x20], #8 ; lsr x9, x1, #7
    add x12, x21, x9, lsl #4 ; ldp x13, x11, [x12]   (giant: ldp x11, x13, [x12])
    extr x10, x1, x0, #62 ; bfi x4, x10, #55, #9 ; stp x3, x4, [x22, #16]
    str x22, [x12, #0]   (giant: [x12, #8]) ; stp x2, x13, [x22], #32 ; cbz x11, next

x13 is the own table's old head (it becomes the link), x11 the other table's head, which enters copy A of s25's routine directly. The key field and k64 are in registers before the lookup, so no extract step remains (s12, s21 deleted).

**Each equal-key pair once.** Rounds are kind-homogeneous and a record is stored before the next record is processed, so for an equal-key (baby key, giant key) pair the record processed later finds the other on its chain and the earlier one did not. Every true pair yields exactly one candidate and every same-slot pair is compared exactly once, so a prefix's comparisons at an owner are sum over slots of B_s G_s. audit_il executes both loops (a64.py) over 12 rounds of both kinds (870 records over 10 heads, 1068 equal-key pairs, 485 of them found by baby records, 21791 comparisons): candidate lists, heads, links and pools equal the full-key oracle, the step count is 12n + 2(n div 64) + 9C + matches + 10, every pair is reported once, one unit short of a round's budget takes the fail exit, and wrong lsr, head offset, field position or compare shift fail it. proto runs the wave protocol over six rounds of both kinds and asserts the same exactly-once property.

**Prices and sizes.** Per baby point s30's walk price + 12 + 2 x 770/64 = 3153843925641/3276800000 = 962.4768; per giant point the walk price + 3 x (12 + 2 x 770/64) - 2 + BETA = 3803059925641/3276800000 = 1160.6018 (pinned). Q = 603/500 (ledger scan near c_g/c_b). Means per owner N_b/2^62 < 2^55.47 and 3 N_g/2^62 < 2^56.78 (s18).

**H1' (prefix coverage; in the claim's conditioning).** For every round prefix j of the production schedule, |R_j|/m >= 1 - e^(-K_j) - 0.0005 with K_j = 6 |B_<=j| |G_<=j|/m. H1 is its last case. Support: source-s7 (two-sided census).

**Expected work.** Rounds are fixed sets, so E[work] = sum_j w_j P(T >= j) with w_j fixed, and by s44 and H1' P(T >= j) <= e^(-K x_(j-1) y_(j-1)) + 0.0005 + FAIL with K >= K_T. At a baby round y_(j-1) <= x_(j-1), so x y >= y^2; at a giant round x y >= x^2. A round adds at most MAXR = GENS x EMAX = 2^104 points and a side ends at most EMAX past its target, so the left Riemann sum gives expected baby points E_b = ceil(N_b (J + 0.0005 + FAIL)) + MAXR + 2 EMAX = 88876914453953062401067662964710513 (0.3876 of N_b) and likewise E_g = 73699081783027707051866047749359910 (0.3876 of N_g), with J = int_0^1 e^(-K_T y^2) dy < 0.386993 (alternating Taylor series in rationals, stopped on a negative term below 10^-50; evidence.gauss_int_up). Point prices are every-run prices (s44's per-class floors and H7 cap) and enter linearly. The allowance is N_b + N_g at the every-run value and E_b + E_g + 2^116 in E_run (the 2^116 covers s3(v)'s pass control, below 2^115.1).

**Comparison budgets.** At each close release every owner receives cap_(j+1) = ceil((1 + 10^-5) 3 n_b n_g/(OWNERS^2 HEADS)) + 2^40, n_b, n_g the cumulative side counts after round j+1, and sets its s25 budget register to cap_(j+1) minus the comparisons it has done (8 instructions per owner per round, row); exhaustion takes the s25 fail exit. Tails: at every prefix the s11 argument (B_o, G_o binomial with deviation 10^-7 mean + 2^30, Hoeffding over giant blocks of at most 64 k_g comparisons with deviation 2^39) has each of its three exponents at least 143 nats (the minimum over the mean G of (G/10^7 + 2^30)^2/(3G) is 4 x 2^30/(3 x 10^7) > 143; Hoeffding above 700), asserted above 120 nats at prefix scales 2^-4k of each side; each owner-round event is below 2^-200, so 2^62 owners x 2^16 rounds stay far inside the 2^-40 union member. So the comparisons of a run that stops after round T are at most OWNERS x cap_T, deterministic in T, and with u = x y, P(T >= j) <= e^(-K_T u_(j-1)) + 0.0005 + FAIL gives E[cap_T] <= OWNERS x [(1 + 10^-5) 3 N_b N_g/(OWNERS^2 HEADS) (1/K_T + du + (0.0005 + FAIL) u_end) + 2^40 + RMAX + 1], du <= 2 MAXR/min(N_b, N_g): u-weight 0.19186, i.e. 0.191885 of the every-run comparison rows with the 2^40 + RMAX + 1 slack, comparison_rows_expected = 339902337058414704253741265923866624 instructions with the cap pair.

**Result.** E_run = W_run - (records, allowance, table, walk and comparison rows) + their expectations above = 172603434747228806072920362889798020392 (asserted below W_run); W(n) = (64 x 2^64 + K E_run)/(2000 x 2^64) <= 2^114.704949928220679 (s43); every opened run at most 2^117.422738478825422 units.

**Evidence.** audit_il and proto above; ec (declared input) walks interleaved bases on a real toy curve and checks every trial that the run stops exactly at the first round whose prefix sets cover t; census runs the two-sided production order in the declared input (m = 2102730, d_t = 64, 23 rounds: lowest prefix deficit -0.007697, -4.5e-05 after the literal-repeat term; shuffled control -0.007697, -4.5e-05; source-s7). All finite; production rests on H1'-H5, H7.

## 46. Montgomery rows of reused operands

Prices and layouts only: points, keys, records, rounds and events stand; H7 is re-measured. Supersedes the s26-s28 and s32 counts, the s1/s27 first product and lane start, and the slot sizes of s8.

**Rows.** For a residue v <= p, Row_0 = v and Row_(j+1) = (Row_j + q p)/2^64 with q = Row_j[0] n0 mod 2^64, j < 5: Row_(j+1) <= (p + (2^64 - 1) p)/2^64 = p and Row_j = v 2^(-64 j) mod p. For any 384-bit a = sum a_k 2^(64 k), S = sum_k a_k Row_(5-k) = a v 2^(-320) mod p and one Montgomery word r = (S + q p)/2^64 (q = S[0] n0) is congruent to a v R^-1, the residue MULNR returns. Rows <= p give S <= 6 (2^64 - 1) p < 2^448, so the 7-word accumulator never wraps and r < 7p (RP6, 186 instructions, row loads included); rows <= 2p (sums of two row sets) take an 8th word and r < 13p < 2^385 (RP7). The correction subtracts m(idx) p, idx = r >> 376 (extr of the top words), m(idx) = floor(idx 2^376/p), read from a 512-entry table of the low 384 bits of m(idx) p: since m(idx) p <= idx 2^376 <= r < (idx + 1) 2^376 < (m(idx) + 1) p + 2^376, the result lies in [0, p + 2^376), p + 2^376 < 1.04p < 2^382.

**Where operands repeat.** Every point that is ever a lane op's P (midpoints and roots, written by the full kernel; base points, strand starts and restarts) carries rows 1..5 of x at +192..431 and of y at -240..-1; gap-table entries also carry V rows of v2 = p - x2 at -480..-241 (V_j = p - X_j, rows.vrows_l: three loads, six subtractions and three stores per row, 60 per entry, executed on every table-level output of batch_audit.run). A midpoint serves all its pairs (s30), an entry all its uses.

**Batch round** (rows.py emitters, batch_audit.py listing). Slot of 496 bytes per lane op: D rows 0..287, c 288, e 336, mu- 384, mu+ 432, P and T addresses 480. Forward element, 292: store P, T; e = v2 - x1 (stored for the kernels); D_j = X_j(x1) + V_j(v2) <= 2p for j = 0..5 (D_0 = d_i of s27), formed pairwise just before row j of the product and stored; c_i = RP7(c_(i-1); D) corrected and stored, and kept in registers for the next element (two register-renamed copies alternate in the unrolled loop; the argument pointer lives on the round frame, ldr/str per element). End: C2 makes c_(n-1) canonical, s15's Fermat, one negation (s27). Backward element, 930: inv_i = MULNR(t, c_(i-1)) (unchanged; a = t < 2p, so K1 holds; inv_i < 1.04^2 x 0.1016p + p < 1.11p); A = RP6(inv_i; Y rows of P), B = RP6(inv_i; Y rows of T); mu- = A + B and mu+ = A - B, each corrected by one table step and stored (the sum, below 14p < 2^385, indexes the unsigned table; the difference, in (-7p, 7p), indexes a signed table at 256 + the top byte of its 384-bit two's complement on a borrow, ranges [0, 182] and [329, 511] disjoint, entries the low 384 bits of floor((idx - 512) 2^376/p) p); t = RP7(t; D) corrected. mu- = (y1 + y2) inv_i R^-1 and mu+ = (y1 - y2) inv_i R^-1 mod p are the slopes of P - T and P + T, both below p + 2^376. t and A spill to the round frame (counted). BATCH = (32 x 292 + 1 + 32 x 930 + 2)/32 = 39107/32; the fixed part is s26's plus the two table addresses and the end slot (172591); the two tables are 64 KiB per generating node, built once per node (row 'row-product tables (s46)', GENS x 2^15 = 2^62 instructions).

**Kernels** (kernel_check.mu_text, rows_text). The lane start (t = y1 + y2 and the first product, 309 instructions) becomes ldr buf and three ldp of mu (minus kernel at buf+384, plus kernel at buf+432); e is read at buf+336. Bounds: mu < p + 2^376 < 1.04p < 2^382 replaces s27's 1.24p, so V < 1.11p, X in (-p, 2.11p), z < 1.22p, inside the existing reduction arms (built for mu < 1.24p, z < 1.26p) (x', y' canonical); K1 holds for the remaining product (w < 2p). The full kernel stores the rows of y' (from registers) and of x' (three loads) before the next lane: 140 + 3 + 140; worst 894 over 12 paths. The leaf builds no rows: 275 to 285 over 4 paths.

**Charges.** KERNEL = 894 (full children, roots, padded pairs, tables), leaf KXC = 275 + 10 x 7/32 under H7, BATCH = 39107/32; gap entries + 60; restarts, base points, strand starts and re-walk rebuilds + 283 each; PAIR = 2 KERNEL + BATCH + 28 + s < 3300 keeps the re-walk step price (asserted). Point prices (pinned): baby 3153843925641/3276800000 = 962.4768, giant 3803059925641/3276800000 = 1160.6018. Memory: full-child slots 672 bytes, entries 912, batch slots 496, two 32 KiB tables per generating node (s4 total).

**H7, re-measured.** batch_audit.xge_run with the exact s46 round model gives the leaf share 0.1848, per-round 0.1445-0.2285 at 512 flags, below the 7/32 cap. Per-round spread is sampling at 512 correlated flags (both children share e); for uniform x1, x2 and Montgomery q the share is analytically about 0.1845, and the cap binds the cumulative share at barriers of at least 2^99 lanes. run_leaf returns the exact X >= p count (66 of 365).

**Evidence.** batch_audit.run executes a 64-lane round of genuine points and a table level: every stored c, e, mu- and mu+ equals the exact big-integer model (rows.py), mu (x2 - x1) = (y1 +- y2) R mod p, and both side calls give exact P -+ T with exact rows; kernel_check runs 509 lanes per kernel with rows checked, a 672-byte footprint audit under canaries and an overlapping-stride control.
