This revision keeps the deferred-tests kernel and its top-limb distinguished test (key < 2^196, theta = 2^196/p,
home slot (key >> 32) mod C, DP_MAX re-derived) and adds six changes: a 4-bit windowed batch inversion (11,814 fewer
instructions per round), a one-instruction 2^52 anchor cap test, the rare-path repricing implied by the shared `finv`
(F_MAX 41,190 -> 29,376; wL 2,789,700 -> 2,021,790), an A=2^32 adder table that cuts the two-cycle reserve, a
three-arm anchor schedule that omits the anchor load and comparison in 16334 of every 16384 rounds (period 64 -> 16384,
with the exact coverage condition proved below), 512-copy batched blocks (backward jump-table), and the 2P+64
persistent-state re-derivation of the collision-exponent margin. PLAIN_FIXED is 154,896,754,179 per 2^28-lane round (the measured 2^20-lane body price 605,065,449 held per lane), C_COMMON 575.041797,
the rare reserve 0.003030230 and the total per step at most 589.211918; CORR is 14.167091. The declared score is 126.728598631. 

# P-256 sigma rho with a 2^32-adder table, a three-arm anchor schedule and a concentrated two-cycle reserve

## Current claim and provenance

For a uniform scalar x in 0..r-1, given Q=[x]P on the protected NIST P-256 profile, the following classical
single-target algorithm has success probability at least 0.39 under the disclosed H1, M1/L-ori, fresh-class exposure
and H-trail premises. Its total charged computation is at most **2^126.728598631** target-group-operation equivalents,
normalized at C=1000 scalar AArch64 instructions per unit, not one unit per P-256 addition. Peak memory is below 2^77
bytes, preprocessing (the rejection-sampled adder table and record-table zeroing included) is at most 2^89.961190
units, nonuniform advice is zero bytes, and instruction accounting retains H3.

The exact declared score is 126728598631/1000000000 (126.728598631), the ledger's total rounded up at the ninth decimal: a
conditional analytic bound, not an unconditional theorem or an executed full-size attack.

The immediate sources are 8337718a (126.735049), supplying the windowed inversion, one-instruction
anchor cap, rare-path repricing and a smaller A table, and 32895fa2 (126.732733), supplying the three-arm anchor
schedule; both build on the top-limb-test revision c1e31e69 and its ancestors. The material
predecessor is b89d1d72 / commit 0e1d63e9, directory submissions/p256-ecdlp/rho-sigma2.

**Precedence:** the claim is this document and `submission.json`. The predecessor's Sections 1.2--1.3 define the same
algorithm except for A, L, the batch size and the rare-work threshold; 2.0--2.3 and Lemma P of 2.8 give the
probability-space, escape and replay arguments; `probability.md` replaces the abandonment/selected-merge bounds and
supplies report overflow; 3--4 the forced-path pricing, resource caps and storage. The batched step, its
measurement, the record-table layout, the two-filter rule and the adder rule are replaced below; the predecessor's
score and parameters are historical.

## Algorithm, fixed costs and changed parameters

Use the protected p, a=-3, b, r and generator P unaltered, return zero for Q=identity and otherwise validate the curve
point. A class is {W,-W}, key(W) the integer Montgomery representation of its affine x-coordinate, canon(W) the
representative with even Montgomery y, and sigma(W) the sign taking canon(W) to W. The walk is

    W -> W + sigma(W) R[key(W) mod A].

Here A=2^32 adders, each R_i=[c_i]P+[d_i]Q with fresh uniform scalar coins, redrawn until 2^192 <= key(R_i) < 2^231 and its
Montgomery y lies within p/2 +- 2^-8 p (at most 2^40 trials; "Small-key adders", "The y3 fold"); the table contains
both signs. K=2^32 walk slots are arranged into batches of 2^28 for simultaneous inversion. Distinguished points have
key < 2^196, so theta=2^196/p under the stated class model. Each trail starts from fresh independent a,b
with W=[a]P+[b]Q, ends at its first distinguished point, records it once with its coins, and a fresh trail starts; coefficients are reconstructed by replay only at a match.

Test ordering and replay are predecessor-proof.md Sections 1.2--1.3 with checkpoint-only cycle detection:
the distinguished-point test precedes the cycle test at every generated point, look-ahead and escape points included;
a detected cycle of length 2..64 escapes by doubling the canonical representative of its least-key member; the 64-bit
anchor filter only triggers exact checks. Every trail is limited to TRAIL_MAX=16384*(2^52+1)=2^66+16384 processed steps; the first
key match is resolved by two replays, each limited to 65 times its original live length. Degenerate matches and
exceeded caps fail the run, and each derived scalar is verified by [x]P=Q before output.

Set N=(r-1)/2+1 and

    tau_target = 4943073/10000000,
    T_min = ceil(sqrt(2*N*tau_target)),
    L = ceil(T_min/K), T=K*L,
    S = K*(L+TRAIL_MAX).

Every slot takes L main-phase steps and finishes its current trail, so S
bounds all processed steps including persistent steps and finishing tails. Sampling caps, maximum distinguished records, abandonment cap, probe cap,
restart costs, bounded replay, verification and output are unchanged formulas evaluated at this S; `accounting.py` recomputes the
point-order certificate proving ord_r(2)>2*TRAIL_MAX.

The common path costs at most C_COMMON=575.041797 instructions per processed step under H3 (anchor schedule below)
and corrections at most CORR=14.167091 more, metered (section "Blinded corrections and their meter"). A fresh-source two-cycle entry and a longer entry
each cost at most wL=2,021,790 more (an anchor check with its batch exit, look-ahead and escape); false filter episodes with their
exposure debt are charged 32834*(w2+wL), w2=71,704 the current two-cycle check, merged steps wL, plus the
inherited forced-path and control allowances. The runtime uses **one total rare-work meter** and need not
recognize whether a point is globally fresh or classify observed events as fresh versus merged; that partition is
solely for analysis.

### The adder table at A=2^32 and the measured step price

Every entry rate and reserve term is a function of A. At A=2^32:

    fresh-source two-cycle rate  p2 = 1/(2A) = 2^-33,
    longer-entry rate            <= 3/(4A^2) + 2/A^3 = 4.065758e-20,
    conditional entry bound      P' = 64/(2A) = 2^-27,
    collision exponent penalty   epsilon = 32836/A + 49248/A^2 + 131328/A^3 = 7.645227e-6
                                 (2P+64 = 32832 wasted states per entry, E2 below).

The two-cycle rate is exact in probability.md's coupled ideal-index experiment under L-ori; the real finite-key residue defect is paid in its 2^-80 coupling charge. The longer-entry rate is at most 3/(4A^2)+2/A^3 (Lemma E1(c)):
E[w'] <= 3/(2A)+4/A^2 over windows of length 3..64, so E[w']/(2A) <= 3/(4A^2)+2/A^3, which the ledger charges and `accounting.py` asserts equals (3/(2A)+4/A^2)/(2A). Epsilon's closure margin charges every short-cycle entry, two-cycles included, 2P+64 = 32832 wasted states (at most 2P-1 <= 32767 persistent steps from E2 plus 64 look-ahead states), with two-cycle term 2*32832/(2A) and longer-entry term
2*32832*(3/(4A^2)+2/A^3). The split-concentration threshold mu=S/(2A), the two-cycle reserve ceil((1+eta)*mu)*wL, E0 and
the reserves built on them are recomputed, not rescaled: the two-cycle charge is 0.000235368 and E0 is 7.452956e-9
instructions per step.

The A-dependence of the compiled step is checked on organizer-image objects (flags as measured, B=2^16,
ADDER_BITS=20, 24 and 27): the three whole-translation-unit mnemonic streams have the same 17,868 instructions and
differ only in the immediates of the DP test and the adder index. `batch_step` has 10,303 instructions at A=2^20 and
2^24 (9,935 before the deferred tests), followed by five (A=2^20) or one (A=2^24) alignment nops that never execute;
`step_rare` (98 instructions; three trailing padding nops at A=2^24), `new_trail`, `add_single`, `finv` and
`dp_restart` are identical up to immediates, and `rare_anchor` has one more instruction at A=2^24, inside the
per-event OVERHEAD. A=2^32 changes only the adder-index mask (2^32-1, a legal logical immediate) and the slot shift;
instruction forms, counts and branch geometry are unchanged (analytic, not measured). The A=2^32 table is 2^39 bytes
for both signs; it is not linked into the finite harness, so no A=2^32 execution is claimed, and its cost is charged
below. The 64-bit anchor filter
compares the whole low limb, so its false-match charge (2^-63 per step, twice a 64-bit filter's 2^-64) is independent
of A; the table is data.

The larger table is charged, not assumed free: construction (worst case of the small-key rejection sampling)
A*2^40*(START+2^10) = 2^89.961189 units; storage A*2*64 = 2^39 bytes for both signs, and 2*A*2*64 = 2^40 bytes for the two tables
(the block reads a second copy), inside the peak-memory formula; each 128-byte record's read-back (16 words) inside START's 2^14 constant. With the
record-table zeroing below, preprocessing is at most 2^89.961190 units, included in W.

### The batched step in this revision

The walk, its distinguished-point test and every rare-path primitive are the predecessor's. The implementation of the common batched step changes
in ways 1-13 below, of which 1-6 and 8-10 are value-preserving:

1. **Branchy corrections.** Every conditional correction of the field
   arithmetic, computed branch-free and always executed in the predecessor,
   is now a forward branch over its fix-up, and a taken fix-up adds its own
   extra instruction count to a meter (below).
2. **A register-resident forward block.** The prefix products PRE[w+1] =
   D[w]*PRE[w] (and Q[w] = dy[w]*PRE[w], item 12), w=1..B-2, run as one
   assembly loop (`FLOOP`) with PRE in x19..x22, never reloaded.
3. **Operand forms.** inv is canonical whenever read (item 9); PRE, Q,
   lambda and (item 13) the square, x3 and y3' stay below 2^256 (their consumers accept that; the tail tests x3, y3') and the products'
   last steps write straight into their destination registers.
4. **2^28-walk batches** and a copy of the two-sign table T2 at virtual address
   2^40, read by the blocks, so the adder address is one `orr` and x27 is
   free for the meter; the C code keeps the static T2.
5. **A blinded prefix chain.** Each global round draws one fresh uniform
   unit rho shared by all batches (rejection sampling, at most 8 draws) and
   every batch starts its chain at PRE[1] = rho*D[0]; the C path for w=0
   takes inv*rho = 1/D[0]. Inverses, slopes and points are unchanged,
   and every prefix PRE[w] is a uniform unit given the walk, which the
   meter's abort bound uses.

6. **Five hundred twelve copies per loop.** Both blocks hold 512
   labelled copies; the end test (`tst x23`/`b.ne`, or `subs x26,#512`/`b.ne`)
   runs once per 512 steps, copy j processing lanes w = -j mod 512. The
   backward block enters through a 512-entry jump table keyed by bits 6-14 of
   x23 (entry v is copy (-v) mod 512); the forward block enters at copy 2 with
   counter B, running 510 + 512(B/512-1) = B-2 steps. A local trampoline keeps
   the 2^52 cap's `tbz` target in the 32 KiB range. Every walk is processed as
   before (B a multiple of 512, asserted at compile time; evidence at the end of `predecessor-proof.md`).
   The 128 -> 512 step is analytic: the measured 128-copy term
   `4B/8 - 4B/128 - 3 = 491,517` at B=2^20 is a four-instruction end test once
   per 128 lanes plus a fixed three-instruction dispatch; 512 lanes halve that
   end-test count twice, adding `4*(2^20/128 - 2^20/512) = 24,576` per round, with the
   six-instruction dispatch charged unchanged.
7. **Checkpoint-only cycle detection** (ae463634). The two-steps-back
   filter is gone; the anchor filter alone catches every short cycle
   (section at the end of `predecessor-proof.md`). The lane record {sel, anchor, tper} keeps the 32-byte
   stride of the implicit-previous-key revision.
8. **Small-key adders, reversed D and t** (5a258d4a; section "Small-key adders" at the end of
   `predecessor-proof.md`).
9. **Deferred canonical test** (e2973980): inv's test r3 >= p3 runs
   in the step tail (section at the end of `predecessor-proof.md`).
10. **Paired reduction words** (21edf4bd; section at the end of `predecessor-proof.md`).
11. **The y3 fold and y-window adders** (y3-fold revision; section below).
12. **Shared rows** (shared-row revision; section below).
13. **Deferred x3 and y3' tests** (`measurement/gen51.py`, sent as d9337f20; section below).
14. **The distinguished test on the top limb** (this revision; section below).

`measurement/step.c` is the predecessor's harness with these changes; the backward block `BLOOP_FIELD_B`
(`measurement/bloop_branchy.h`) and the forward block `FLOOP` (`measurement/floop.h`) are generated by `gen10.py` from
the predecessor's measured `BLOOP_FIELD` text (the generator stays in 70c0711e, `measurement/gen10.py`,
omitted for the size limit; this revision's blocks are its output after 5a258d4a's `measurement/gen11.py`,
`measurement/gen12.py` (e2973980), `measurement/gen48.py` (21edf4bd), `measurement/gen13.py` (205a498f), `measurement/gen49.py`
(9116d890) and `measurement/gen51.py` (the deferred-tests revision, sent as d9337f20), all omitted likewise). Code changed: item 7 the step tails, rare exits, `step_rare`,
`new_trail`, `main` (`rare_2cycle` deleted); item 9 the inv finish, tails' last branch, rare exits; item 10 the product rows;
item 11 the backward block's end; item 12 the forward block, backward first products, `step.c`; item 13 the tails, slow blocks, asm operands.

### The corrections

With x28 = p1 = 2^32-1 and x29 = p3 = 2^64-2^32+1 set at block entry (p = (2^64-1) + p1*2^64 + p3*2^192), adding p
to four limbs is `subs r0,r0,#1; adcs r1,r1,x28; adcs r2,r2,xzr; adc r3,r3,x29` and subtracting it is
`adds r0,r0,#1; sbcs r1,r1,x28; sbcs r2,r2,xzr; sbc r3,r3,x29`.

* **Subtraction** a-b of canonical b and a < 2^256 (six per walk step: D, dy, two in x3, t, y3): `subs/sbcs` x4, then `b.cs`
  over [add p; meter += 5]. The common path is 5 instructions (was 11; y3's 7 with the fold), a borrow costs 5 more and gives
  a-b+p in [0, p), else a-b: below 2^256, congruent, canonical if a is. D, dy and the kept y3 fold take
  canonical a; since item 13 the x3 pair and t take a = S, S - x(W) and x3, each v or v + p below 2^256 (t <= 2^256 - 1 - 2^192).
* **Carry-only finish** (PRE, Q, lambda, inv; the square and y3' since item 13): T = r + c*2^256, c in {0,1}, T < 2^256+p, c the carry flag left by the
  product's top `adcs`. `b.cc` over [subtract p; meter += 5]: 1 instruction (was 7 or 9), 5 more when c=1, giving
  T - c*p < 2^256.
* **Canonical finish** (y3' = t*lambda before item 13, inv before item 9; T < 2p): `b.cc L` (c=0), else subtract p, meter += 4, `b done`;
  `L: cmn r3,x28; b.hs slow`. Common path 3 (was 9), c=1 costs 4 more. If c=0, T < p unless r3 >= p3 (`cmn` sets C
  exactly then); the out-of-line slow block computes T-p, keeps it with four `csel` on no borrow and adds 10.
* **lambda^2** before item 13, from lambda < 2^256: T < 2^256+p, so after a c=1 subtraction T-p may still be at least p; the fix-up
  (meter += 5) falls through into the same `cmn/b.hs` re-test, and one more conditional subtraction in the slow block
  suffices (T - 2p < p).

### Flag-carry finishes (from 1f942c41)

The product finishes branch on the carry flag C: c2b498ec ended each product with `adc c,xzr,xzr` (c = x13, x11 at
lambda^2) and `cbz c,L`, and 1f942c41 made each pair one `b.cc`, which branches exactly when C = 0; that keeps every
outcome if the branch follows the `adcs` whose carry is c and nothing later reads c (x13 and x11 are dead at the
square's and y3''s sites; gen49.py's products never formed a carry register). On this revision's four builds (B=2^20
at A=2^20 and 2^24, B=2^16 at A=2^16 and 2^20) all 80 `b.cc` sites (two per forward copy, four per plain and anchor
backward copy) directly follow the top `adcs` (16 `b.lo`, the same encoding, follow `ccmn`). Rows' domains: "Shared rows". Every value a step ends or
exits with (x3, y3, key, sign bit) is canonical (raw-x3 tests: item 13).

**Exactness checks** (5a258d4a's `evidence/measurement/branchy.json`, keys `previous_revision` and earlier,
omitted here for the size limit, quoted values in `ancestor_quotes`; on the builds that introduced these changes,
all with the two-steps-back filter; reruns: sections "Deferred canonical test", "Paired reduction words", "The y3 fold" and "Shared rows"): (i) n=8 dumps of 2^16 walks on
two input sets and stress builds with 2^4, 2^6 and 2^10 adders equal the predecessor kernel's byte for byte; (ii) the
backward block equals an exact big-integer model (Montgomery integers, branch outcomes, meter) on 240,000 random
and crafted inputs; (iii) the forward block on 6,000 runs; all three rerun for 1f942c41, 0 mismatches.

### The y3 fold and the y-window adder table (y3-fold revision)

**Kernel change.** The backward block ends with y3 = y3' - y(W), y3' = t*lambda and y(W) canonical, and a `b.cs`
borrow fix-up (extra 5). 205a498f's `measurement/gen13.py` (omitted here) derived that revision's `bloop_branchy.h`
from the paired-reduction revision's by one edit and its header line: `cmp x30,x11` compares the top limbs y3'[3] and y(W)[3]; when y3'[3]
is smaller a four-limb chain stores y(W) - y3' = -y3, which cannot borrow, and branches past the old code; otherwise
(`b.hs`) the old y3' - y(W) and its fix-up run, the fix-up firing only at equal top limbs (slow site `y3_tie`). The
stored y is y3 or p - y3: |y3' - y(W)|, or y3 at a tie. Each of the 16 backward copies gained these 7 instructions (`batch_step` +112
at all four shapes); either path runs 2 more instructions per block step.

**The walk is unchanged up to sign.** The stored y feeds only the sign bit (the parity of its Montgomery value, which
selects R_i or -R_i and flips with the negation, p being odd) and the next step's y(W). A point stored as -W' selects
the negated adder and yields the negated old successor, so every class, key, adder index, distinguished-point test
and rare-path outcome is the paired-reduction kernel's; x3, lambda, t and PRE are unchanged, and M1, L-ori and the entry
rates, which concern the class map, are untouched. Dumps (each walk at the checkpoint-only kernel's point or its
negation), mutants and block checks: 205a498f's `branchy.json`, `y_fold`.

**The adder rule.** Step 2 of the procedure accepts R when R != O, 2^192 <= key(R) < 2^231 (the top limb of the Montgomery x
in [1, 2^39); the floor bounds t, "Deferred x3 and y3' tests") and Y_LO <= Y <= Y_HI for its Montgomery y Y, with Y_LO = (p - (p>>7))/2 and Y_HI = p - Y_LO, i.e.
|2Y - p| <= p>>7: three compares inside each trial's 2^10. This replaces the test key < 2^224 of "Small-key adders",
whose arguments hold with S the set of logarithms s for which [s]P passes. As 2Y_LO - p = -(p>>7) and 2Y_HI - p = p>>7,
the window is closed under Y -> p - Y, so S = -S, and u = Y/p satisfies |u - 1/2| <= U_DEV <= 2^-8 (`accounting.py` asserts both). D and the second x3 subtraction borrow below 2^231/p < 2^-24.9999 for every adder, charged without pad; the floor only
removes keys, so every use of 2^231 as an upper bound stands.

**Acceptance count** (asserted in `accounting.py`). With I = [2^192, 2^231) and J = [Y_LO, Y_HI], |J| < 2^249, the 2-D
finite Fourier count of affine points with Montgomery coordinates in I x J has main term |I||J|(r-1)/p^2 >
2^224 - 2^185.0001. For (a,b) != 0, f = aRx + bRy has its only pole at O (order at most 3), so Bombieri's bound gives
|sum_P e(f(P)/p)| <= 6 sqrt(p); with sum_{a != 0} |Ihat(a)| <= pL for any interval, L = 1 + ln p, and the same for J, the error is at
most BOX_ERR = 6(isqrt(p)+1)(2L + L^2) < 2^145.5598 < 2^160/2^14. So at least M_ACC = 2^224 - 2^185 - 2^160 - 3 points pass
(it replaces 2^224 - 2^160 - 3 in predecessor-proof.md's "Small-key adders"); Q_ACCEPT_LOWER = M_ACC/r and DELTA_LOOP
move by about 2^-39 relative, LOG2_ADDER_EXHAUST by under 10^-9, and every bound quoted on them holds; the trial cap is unchanged and the table charge is 2^89.961189 units; a constant 2^14 times Bombieri's would still suffice (asserted).

**dy under the fold model.** dy = y(R) - y(W) (`ldp x3,x4,[x7,#32]`, the selected adder's y, minus `[x24]`, the
stored y) borrows iff y(W) > y(R). Let v = y(R)/p for the adder and sign the step selects, so v = u or 1 - u. Model
(support for H-trail, as uniform operands are for the other sites): each step's y3' is uniform on [0,p) given the
past, its parity independent of its magnitude. A stored y3 or p - y3 then has a fair-coin parity independent of every
magnitude, so v = u or 1 - u with probability 1/2 each; and given everything through the stored y two steps back
(c*p, m = min(c, 1-c)), the stored y one step back has density 2 on [0,m) and 1 on [m, 1-m] in units of p, so
y(W), its fold with a fresh y3', exceeds v with probability (1-v)^2 for every c. The rate is ((1-u)^2 + u^2)/2 = 1/4 + (u-1/2)^2 <= 1/4 + 2^-16 for every
such history (not given the stored y one step back) and every adder, not an average over the table; `accounting.py` asserts 1/4 + U_DEV^2 <= q_dy = 1/4 + 2^-16 and charges q_dy + 2^-17. With c
the stored y one step back, the model gives E[dy/p] = 1/2 - m^2 when m <= min(v, 1-v) and at most 1/4 + 2^-8 otherwise, so
in that revision lambda = dy*iw (iw a uniform unit given the walk) kept its fifth-word bound 1/4 given any
history, measured 0.208; here Q = PRE*dy has that form ("Shared rows"). The law needs two fold outputs with fresh y3': a trail's first two steps
(the start point's y is not one) and the steps around a short-cycle entry (up to 2P-1 <= 32767 persistent steps repeat values;
the escape point's y is canonical) are charged the full 5, DY_EXCESS = 5(2K_TRAILS/S + 32770(p2 + 3/(4A^2) + 2/A^3)) <
1.938e-5 per step, below the dy pad 5*2^-17 (both asserted). Walk 0 takes the compiled C step every round and lane B-1's dy is a C
subtraction ("Shared rows"), both branch-free and unmetered.

**The y3 charge.** The fix-up left at y3 fires only if y3'[3] = y(W)[3], for a uniform canonical y3' probability at
most 2^192/p < 2^-63.99 whatever y(W); it is charged 2^-24 with weight 5 in CORR_SLOW beside the three slow blocks.
The y3 site (5 at 1/2 + 2^-10) leaves CORR_SITES and dy moves from 1/2 without pad to 1/4 + 2^-16 + 2^-10. C_MAX stayed
86 (y3_tie's 5 replaced y3's); CORR falls from 17.709638 to 13.959715 and ceil(CORR) from 18 to 14.

**Measured** (205a498f's `branchy.json`, `y_fold`): FIXED +2(B-1) against the paired-reduction kernel in every
counted round, both input sets, plain and anchor (2,097,150 at 2^20); over 127 further rounds of the fold built on the
checkpoint-only kernel, -48 and +144 from net escape negations (24 instructions each). Site rates, rounds 2..64 on
four tables: dy 0.249782 to 0.250142 (checkpoint-only 0.50009 at rounds 2..8), lam 0.208314 to 0.208426; binned by
octile of the stored y two steps back, dy stays at 0.2456 to 0.2543, flat at the model's 0.250005 (about 2.2 binomial
standard errors at most); binned by the stored y one step back it follows |1/2 - beta| (0.4385 down to 0.0607). The y3-fold kernel's meter: 13.702 to 13.715 per walk step (rounds 2..8 and 2..128, two
tables), below CORR; later kernels': "Shared rows", "Deferred x3 and y3' tests".

### Shared rows (shared-row revision)

**Rows and products.** As -1/p = 1 mod 2^64, X*2^-64 mod p is (X + x0*p)/2^64 = X>>64 + x0*(2^32 + 2^128*p3), seven
instructions (`lsr`, `adds` with `lsl`, `adcs`, `mul` and `umulh` against x29, `adcs`, `adc`). The two products of a
pair with factor R use B1 = (R + r0*p)/2^64 and B2 = (B1 + b1_0*p)/2^64, formed once (14 instructions). A product with
scalar words a0..a3 sums S = a0*B2 + a1*B1 + a2*R + 2^64*a3*R (= a*R*2^-128 mod p; three rows at word 0, one at word
1) and divides by 2^128 with two Montgomery rounds (m0, m1 the low words, their p3 multiples by `mul` and `umulh`,
m0>>32 and m1<<32 joined by one `extr`): 78 instructions against the paired-reduction product's 90. The row idea is
secp256k1-ecdlp submission a00774d8's (rows of the running inverse shared by its two products), carried to Montgomery
form. `measurement/gen49.py` writes both blocks from this schedule with a liveness allocator (pool x0..x18, x30;
lambda in x0, x2, x18, x30 and INV', PRE' in x19..x22 as before), each edit asserted to match once.

**Dataflow.** iw = INV*PRE is gone. The forward block computes dy (the subtraction and fix-up the backward block had,
site dy), Q = PRE*dy (finish site q) into the pre slot, then D and PRE' = PRE*D, both products on PRE's
rows. The backward block computes lambda = INV*Q (site lambda) and INV' = INV*D (site inv) on INV's rows, then
`ldr x7,[x26]` and the y3-fold kernel's text from the square on; lambda = INV*PRE*dy = dy/D as before. Per walk step
the backward block's common path loses 110 instructions (287 -> 177) and the forward block's gains 90 (105 -> 195); a
backward copy has 420 instructions (540), a forward copy 215 (115). `step.c`: FLOOP clobbers x0, x2, x18, x30;
`batch_step` stores PRE(B-1) from the forward block's registers into DP[B-1].pre (the inverse's input, as before) and
forms Q(B-1) = PRE(B-1)*dy(B-1) with `k_sub` and `k_mul` (branch-free and unmetered, so lane B-1's dy fix-up leaves
the meter); the scalar tail's w > 0 branch, unreachable (RUN_BLOCK always leaves w0 = 0) and reading PRE from the pre
slot, is dropped.

**Values.** Expanding the rows, T = (a*R + N'*p)/2^256 with N' = L + 2^128*(m0 + 2^64*m1), L = a0*r0 +
2^64*(a0*b1_0 + a1*r0) < 2^193 + 2^128, while the exact Montgomery integer uses N = -a*R*p^-1 mod 2^256. As N' = N
mod 2^256 and N' < 2^257, N' = N + k*2^256 with k = [N < L], so T = T_exact + k*p. Domains: backward R = INV < p at
every block step (`finv` ends in the C kernel's canonical finish and the tail test canonicalises each INV'), scalars
below 2^256; forward R = PRE < 2^256, scalars dy, D < p (top limb at most p3). Then a row stays below p + 2^192 (below
p if R < p), S + m0*p < 2^384 (six words never overflow) and T < 2^256 + p, so the carry-only finish leaves T - c*p <
2^256 (all asserted in `accounting.py`). The walk is unchanged: lambda's integer, which may differ from the y3-fold
kernel's by p, enters only the square and y3', whose finishes were canonical before item 13; INV' is canonicalised before it is read;
Q and PRE' enter only as multipliers or rows (any value below 2^256) and PRE(B-1) the inverse. x3, y3, keys, signs and
each step's INV equal the y3-fold kernel's; no exit was added and C_MAX stayed 86.

**Corrections.** Every site keeps its weight and idealized rate (site q replaces iw), so CORR is unchanged: Q = PRE*dy
has fifth-word probability at most E[PRE]*dy/2^512 < dy/(2p) given dy (up to 1/2 as dy nears p), so at most
E[dy/p]/2 <= 1/4 given any history up to two steps back, as E[dy/p] <= 1/2 under the fold model ("The y3 fold",
the bound lambda = dy*iw had); the steps DY_EXCESS names add below 5e-6 per step to q (5*(1/2-1/4) times their
share DY_EXCESS/5), inside its 2^-17 pad; lambda = INV*Q has two blinded operands,
as iw = INV*PRE had; k = 1 needs N < L, below 2^-62.99 under the uniform-multiplier model, inside each site's 2^-17
pad. The held PRE integer is pi or pi + p as before, so the C_MAX*2^-31 term stands. Q, lambda and INV' are other
integers than the y3-fold kernel's iw, lambda and INV', so the meter is not that kernel's (per-site
rates: 9116d890's `branchy.json`, `shared_rows`); the totals stay within noise of it (meters there).

**Measured** (9116d890's `branchy.json`, `shared_rows`, with its exactness checks, mutants, dumps and meters;
organizer image). 80 `b.cc` sites, each right after its top `adcs`; memory operands encode. At B=A=2^16 a round falls by
1,310,630 = 110(2^16-1) - 90(2^16-2) - 160 (Q(B-1)'s `k_sub` 18, `k_mul` 115 and 27 compiled stores and moves), and 8
rounds by 8 times that, on a y-window table and both uniform input sets; 2^20 rounds: next section.

### Deferred x3 and y3' tests (this revision)

**Kernel.** `measurement/gen51.py` (sent as d9337f20) writes `bloop_branchy.h` and `step.c` from the shared-row revision's (9116d890),
each edit asserted to match once. The square's and y3''s finishes become carry-only (their `cmn r3,x28; b.hs`, slow
blocks Ls5, Ls9 and the `b` past y3''s test go); x3[3] stays in x17 and the fold writes y3 to x13..x16, and the tail
tests `tst`/`ccmp` (event: Z = 1), then `ccmn x17,x28,#10,ne`, `ccmn x30,x28,#10,lo`, `ccmn x22,x28,#10,lo`, `b.hs
Lz`: x3[3], y3'[3], INV'[3] >= p3, each only if nothing fired (else NZCV = 1010). INV' is last, so Lz's `b.mi` stays
exact (N = 0 iff only INV' fired; its sum less 2^64 is below 2^32) and its INV' block (11) is unchanged. N = 1 goes to
Lx, which re-tests x3[3] and y3'[3]: neither is the rare exit as before; else the fix block subtracts p from x3 and
y3' on no borrow (`csel`), stores x3, refolds and stores y3 (both folds from one generator function), zeroes x17 and
x30 and re-runs the tail. The anchor tail's cap test moves after the chain, so every exit sees
canonical x3, y3, INV'. `tst`/`ccmp` read the raw x3, but an event (NZCV = 1010) or x3 = X + p (top limb at least p3) reaches Lx,
and key and sign follow the chain: exits, keys and signs rest on x3 = X, y3' = Y, and a test of
X + p is dropped and re-run on X.

**Values.** The carry-only finish leaves the square S below 2^256: S_c or S_c + p. A borrow at x3 = S - x(W) - x(R)
or t = x3 - x(W) adds p to a value in [-p, 0), giving [0, p), and later subtractions of canonical operands keep it there;
so x3 = X or X + p (below 2^256 < 2p), and t >= p only with no borrow at all three: t = S - 2x(W) - x(R) <= 2^256 - 1 -
2^192 by the adder floor. The y3' product's five-word accumulator, whose top word takes non-flag `adc`s, then stays below
t*2^64 + p + 2^128 < 2^320 for every lambda < 2^256 (asserted in `accounting.py`; the p*2^64 bound of "Paired reduction
words" assumed t < p). Without the floor, t > 2^256 - B, B = (p + 2^128) >> 64 < 2^192, and a lambda word 2^64 - 1 can
lose a carry (y3' off by 2^192). So y3''s T < 2^256 + p and y3' = Y or Y + p. x3[3] < p3 gives x3 = X, y3'[3] < p3 gives
y3' = Y; a fired test restores both. So x3, y3, keys, signs, every acted-on test and exit, and INV' are the shared-row kernel's
on every path; only the meter differs, on rare steps. Checked (`branchy.json`, `tail_tests`, `domain`):
a solver proof finds no overflow at the product's 10 non-flag `adc`s for t[3] <= 2^64 - 2 (three if t[3] = 2^64 - 1 is
allowed) nor at the square's 8; 1,150,000 emulated and organizer-build cases up to t's maximum equal the shared-row
walk and the exact model; every measured table meets the floor.

**Charges.** The common path loses two `cmn`/`b.hs` pairs and gains two `ccmn`: 2(B-1) fewer per round. y3''s c=1
path runs 5 over it (was 4: no test left to skip). CORR_SLOW drops sq_slow and y3m_slow (10 each) and adds, at 2^-24,
xy_fix: Lz's `b.mi` to the re-run's `b.hs`, 37 + 7 + BL_FM_N + BL_DP_N (asm operand `wfix`; 46 measured, charged 47,
any build's most), and y3_tie_fix (5). x3[3] >= p3 needs a square at least p (the old square slow event, below
2^-31) or X in [p3*2^192, p) (length 2^96 - 1, below 2^-159 at density below 1/p); y3'[3] >= p3 needs the old y3'
slow event (below 2^-31), a c=1 result in that interval or t >= p (the x3 event): below 2^-29 per step. y3''s c=1
fix-up and the fix (after x3 = X + p), and both ties, can share a step, so C_MAX = 11*5 + 11 + 47 + 5 + 5 = 123 (was
86, the sum less 4). CORR < 14.167091, ceil 15; the meter bounds keep their structure ("Meter and cap").

**Exits.** A distinguished-point or anchor exit runs Lx's 3 instructions first, 12 above the clean step (was 9); a cap
exit at most 12 (was 7). Charged 13 per batch exit (was 10;
emulated, `branchy.json`); wL 2,021,790.

**Measured** (`branchy.json`, `tail_tests`; organizer image). `batch_step` 10,303 instructions at both A; finish `b.cc`
follow their `adcs`; other functions unchanged; operands encode. 156,000 emulated random and crafted cases (three
builds, both tails) and 24,000 injected states give the shared-row walk and the exact models, executed less metered
instructions one constant per kernel, tail and build off exits; 2,400,000 organizer-build cases equal; twelve
one-instruction mutants fail 16,016 to 332,076 of 100,000 per seed, the control none. A round falls by 2(B-1) at
B=A=2^16 (1 and 8 rounds; a y-window table, both uniform sets) and 2^20; ten dump configurations equal the shared-row
kernel's byte for byte but fix_meter. Meter per walk step, y-window adders, rounds 2..8 and 2..128: 13.963 and 13.954
(one table), 13.956 and 13.962 (the other; rounded up), below CORR.

### The distinguished test on the top limb (this revision)

**Kernel.** The base submissions' `gen55.py` writes the unrolled `bloop_branchy.h` and `step.c` from the deferred-tests revision's (d9337f20),
each edit asserted to match once. The DP test stops computing a straddling key: with bits 24..83 zero it read
`extr x8,x4,x3,#24; tst x8,#dpf`. Distinguished is redefined as the key below 2^196, i.e. x3's top
limb below 2^4 = 16: one `tst x17,#0xfffffffffffffff0` on the limb the deferred chain already tests. The adder index
stays key bits 0..26, so the walk's class sequence is unchanged. Two repairs follow from the fix path. The y3'
canonicalisation's temporary moves off x17 to x5 (dead after `stp x5,x17,[x25,#32]`), so x17 keeps the canonical top
limb; and Lx's `mov x17,#0` becomes `and x17,x17,#%c[dpf]`, which clears bit 0 of a value below p3, so it stays
below p3 (never the rare exit, no loop), preserves the DP bit pattern, and makes the re-run from Lt test canonical x3.
The stores precede the mask, so x3[3] is stored canonical.

**Set identity.** For canonical x3 = X < p both tests are functions of X alone, so a recorded class carries one
distinguished flag under either rule; what changes is the index at which trails end, and so theta. As 2^196 < p the set
is a full interval [0, 2^196) with theta = 2^196/p >= 2^-60 (asserted in `accounting.py`). In the assembly DP exit, the overlapping anchor bits cannot influence the decision: at a distinguished point the
`tst` leaves Z = 1 and C = 0, so the following `ccmp x3,x6,#4,ne` and the `ccmn ... ,ne` links are skipped and the
comparison the anchor filter would make never executes; the later CCMN links test x3, y3' and INV' and are skipped via fallback flags on a DP. A
distinguished key's bits 24..63 reach the walk only through the home slot (4,000 pairs, every frame word but the two
anchor passthroughs equal, a control flipping the exit flag 4,000 of 4,000 times).

**Home slot and the record bound.** The slot becomes (key >> 32) mod C. A distinguished key is uniform on [0, 2^196)
under M1, so key >> 32 is uniform on exactly M = 2^164 values, with no p-granularity correction, and
DELTA_HASH = 2*DP_MAX*(C+1)/2^164 < 2.9464e-8. The record bound is re-derived at the new theta:
DP_MAX = 2*ceil(S*theta) + 2*K, asserted at least 2*S*theta + 2*K, so the mean S*theta sits at DP_MAX/2;
DP_MAX is 2^68.491741, below the 2^69.492 slot count. The cap
CAP = 16384*2^52 = 2^66 stays fixed. `probability.md` replaces both old abandonment/merge losses by one
all-abandonment envelope and proves report overflow via distinct first exposures; the residual and correction terms retain their formulas.

**Charges.** The common path loses one instruction per step in every build: PLAIN_FIXED 606,641,931 falls to
605,581,542. CORR does not move: no fix-up fires at top limb = p3 (rate below 2^-32), so xy_fix keeps its charged 47
against the measured 45, and C_MAX stays 123. An exit still runs Lx first, charged 13.

**Measured** (`branchy.json`, `top_limb_dp`; organizer image). `batch_step` gains one instruction at every A; the DP
`tst` replaces the `extr`/`tst` pair in all sixteen copies and the slow blocks; operands encode.
~300,000 emulated random and crafted cases over six build configurations and both tails, with x3 = X + p and the
top-limb boundary set (p3, p3-1, 0xf, 0x10, 0), agree with their exact models (walks compared only where the DP predicates agree); 2,400,000 organizer-build
cases equal, including the 45-vs-46 fix weight; a round falls by exactly B-1 at B = A = 2^16 (three inputs, 1 and 8
rounds) and by (B-1) per round at 2^20 on both seeds; anchor rounds fall by the same 65,535, so ANCHOR_EXTRA stays 4(B-1)+9 for that change (and becomes 3(B-1)+9 with the one-instruction cap below).

### Windowed inversion, one-instruction anchor cap, and rare-path repricing (this revision)

`batch_step` performs one inversion per global round. The inherited binary square-and-multiply over p-2 (255
squarings, 127 multiplications) is replaced by a 4-bit windowed walk: the top nibble is 15 (accumulator a^15), the
odd table a^1..a^15 costs one squaring and seven multiplications, and each of the 63 lower nibbles costs four
squarings plus, when nonzero (32 are), one multiplication. That is 253 squarings and 39 multiplications, saving
11,814 emulated instructions per round in the organizer image; the saving is per round and B-independent.
`add_single` and all its rare-path callers share the same `finv`, so every unbatched addition is exactly 11,814
cheaper; with F_MAX = 41,190 - 11,814 = 29,376, the rare-entry price wL (65 F_MAX terms) falls from 2,789,700 to
2,021,790. ESC65_MAX and ESC2_MAX are left at their measured values, so wL is conservative. The 2^52 anchor cap
test takes one `tbz x12,#52` on the normal path, with a local `b 2f` trampoline on the cap path (same predicate: the
counter never exceeds 2^52+1), so ANCHOR_EXTRA falls from 4*(2^20-1)+9 to 3*(2^20-1)+9.

Paired organizer-image counts at 2^16 walks: plain 39,222,158 -> 39,210,344 (11,814); anchor 39,485,165 ->
39,407,816 (77,349 = 11,814 + 2^16-1); meters and event counters equal; walk dumps byte-identical. The ledger then
gives PLAIN_FIXED = 605,065,449 per 2^20-lane round at the measured batch. The production batch is 2^28 lanes, so the per-lane body price is unchanged and PLAIN_FIXED is 154,896,754,179 per production round; the per-round control terms below are re-amortized at that batch. Together with the anchor schedule, the 512-copy jump-table blocks and A=2^32, the ledger gives
C_COMMON = 575.041797, two-cycle charge 0.000235371, E0 7.452956e-9, reserve 0.003030230 and total per step at most 589.211918.

### Three-arm anchor schedule (this revision)

Anchor rounds are `t mod 16384 == 0`; rounds with
`1 <= t mod 16384 <= 16334` use a tail without the anchor load and comparison (saving `2*(B-1)-2` per skip round,
with dispatch charged +2); rounds 16335 through 16383 retain the plain comparison and pay the +3 dispatch, and the
anchor round compares against the previous anchor and then refreshes it.

**Coverage theorem.** Let the walk enter a cycle of length l at t_c and let t0 be the next anchor round,
0 <= t0-t_c <= P-1, so W(t0) is on the cycle. For offsets d = t-t0 in the tested window U+1..P the walk compares
W(t0+d) with W(t0). With P = q*l + r, 0 <= r < l, the largest multiple of l not exceeding P is P-r, and it lies in
the window iff r <= P-U-1 = 49. `accounting.py` checks P mod l <= 49 < 50 for every 1 <= l <= 64 (maximum 49 at
l=55), so every cycle of length l <= 64 is detected within the period starting at t0: at most (P-1)+P = 2P-1 <= 32767
persistent steps per entry, charged as 2P+64 = 32832 states in epsilon. P is a power of two, so the predicate is one `and`/`cmp`:
arm prices are unchanged, only the per-period counts (16334 skip, 49 plain, one anchor) move.

### Fixed path and metered corrections

Counts are organizer-image QEMU 8.2.2 guest instructions (GCC 14.2.0-19, -O3 -mgeneral-regs-only), as in the
base submissions; the flag-carry kernel was cross-checked against Callgrind.
PLAIN_FIXED=605,065,449, ANCHOR_EXTRA=3*(B-1)+9 and the 1,024 per-round allowance are charged.

### Blinded corrections and their meter

Each site j has extra cost w_j when taken and an idealized probability q_j for a uniform point and prefix: borrow 1/2 for the first x3 subtraction and t; below 2^231/p for D and the second x3 subtraction (small-key adders);
1/4 + 2^-16 for dy (history-conditional bound derived and measured in "The y3 fold and the y-window adder table");
fifth word set 1/4 for the products (PRE, Q, inv, lambda, y3') and 1/3 for the square; below 2^-31 for the inv slow block
(r3 >= p3 needs T in an interval of length below 2^224 + 2^192, and T has density below 1/p), below 2^-29 for the fix
block and below 2^-24 for the y3 ties. The walk step is charged

    CORR = sum_{j not in Z} w_j (q_j + 2^-17) + sum_{j in Z} w_j q_j + 68*2^-24 + 123*2^-31 = 14.167090692 instructions,


Z = {D, second x3 subtraction}, w_inv = 5 and the inv slow test 11 (section "Deferred canonical test"), dy at
1/4 + 2^-16 and y3 replaced by y3_tie (section "The y3 fold"), w_y3' = 5, xy_fix 47 and y3_tie_fix 5 (section "Deferred
x3 and y3' tests").

**What a lane's count depends on.** Let pi = PRE[w] mod p, the lane's prefix residue. With inv canonical when read, everything
a lane computes in a round is a function of its walk point W (with its adder) and of the integer held for PRE[w],
which is pi unless it is pi + p: dy, Q, D, PRE[w+1], inv = 1/(pi*D) up to a fixed factor, lambda, the square, x3, t,
y3' and y3 (lane B-1 forms dy and Q in C, unmetered). Write c(pi, W) in [0,123] for the lane's correction count when PRE[w] is held as pi. The value pi + p
(below 2^256) needs pi < 2^256 - p < 2^224, so for a uniform pi it occurs with probability below 2^-31 and changes the
count by at most 123, covered by the charge 123*2^-31 per step. The walk (x3, y3, keys, signs) does not depend on pi
or rho.

**Premise H-trail (average only).** For a trail started from fresh uniform coins, with the table fixed and any slot
and starting round, the expected sum over its batch steps of g(W) = E_{pi uniform}[c(pi, W)] is at most CORR
times its expected number of processed steps: average behavior only. Support: for a uniform point and prefix, uniform
canonical operands borrow with probability (p-1)/(2p); D = x(W)-x(R) and u - x(R) borrow only below x(R) < 2^231,
probability below 2^231/p for every adder, so no table-mean pad is needed (the predecessor's table-mean Hoeffding
term DELTA_TABLE = 4 exp(-2A 2^-20) is kept in DELTA_CORR although no site uses it); dy's rate is bounded given
any history up to two steps back under the fold model ("The y3 fold"); modelling the reduction multiple m = -ab p^-1 mod 2^256 as
uniform gives P(T >= 2^256) = E[ab]/2^512 (1+O(2^-32)) <= 1/4 (1/3 for a square), for Q = PRE*dy at most
E[PRE]*dy/2^512 < dy/(2p), at most 1/4 given any history up to two steps back as E[dy/p] <= 1/2 under the
fold model; lambda = INV*Q has two blinded operands, as iw = inv*PRE had; a shared-row
product exceeds the exact integer by p only when N < L, below 2^-62.99 under the same model, inside each site's 2^-17
pad; the inv slow block is below 2^-31, the fix block below 2^-29.
Along actual walks the measured per-step mean is 13.953 to 13.963 with y-window adders (13.703 to 13.712 on the
shared-row kernel, 13.702 to 13.715 on the y3-fold kernel, 17.662485 to 17.667348 on the deferred-inv kernel). H-trail bounds only an abort probability,
never the work.

**Meter and cap.** The blocks add each taken fix-up's extra count to x27;
each round adds the block totals to one 192-bit run-wide meter M, inside
the per-round control allowance. After round k, with n_k processed steps
so far, the run aborts if

    M > n_k*CORR + LAMBDA_CAP,
    LAMBDA_CAP = 2^112 + 2*(123+15)*K*TRAIL_MAX + 2^96.

So total correction work is at most S*CORR + LAMBDA_CAP plus one round's
overshoot K*123 (123 = most correction instructions in one step), all charged: CORR inside the step charge, the rest (below 2^102.059 units) as
`correction_meter_slack`; the blinding draws and products are charged in
the per-round control (8*2^16 per global round spread over the 2^4
batches, plus 2^11 per batch round). Work is bounded
even if H-trail fails.

**Abort probability without a walk model.** Write M = sum_s g(W_s) + sum_r N_r, s over processed batch steps and
N_r = sum over round r's lanes of (c(pi, W) - g(W)), round r's blinding noise.

* *Blinding noise.* The walk, hence every W and every PRE up to the factor rho, is independent of the rho draws.
  Given the whole walk, each global round's rho is a fresh uniform unit, so each lane's pi = rho*P[w] is uniform and
  E[N_r | walk] <= 2^28*123*2^-31 (the non-canonical prefix case); distinct rounds have independent rho, and one
  round's noise over all batches is at most 2^32*123. Hoeffding's maximal inequality over at most 2^95.492 rounds gives
  P(max partial sum >= 2^96) < 2^-1000.
* *Trails.* A trail's batch steps, and so its sum T_i of g, are a deterministic function of its start coins and the
  fixed table (the walk, its escapes, its distinguished-point and cap rules). Trails start from fresh, independent
  coins, so given the table (and each trail's slot and starting round, which fix its lane and anchor phase) T_i and
  its length len_i are independent of all earlier trails, whatever the walk's internal dependence; merges and cycles
  do not change this. Every slot finishes its current trail, a trail has at most TRAIL_MAX steps, and at most
  K_TRAILS = DP_MAX + 2K trails start, so D_i = T_i - CORR*len_i lies in an interval of width (123+15)*TRAIL_MAX (15 = ceil(CORR)) and,
  by H-trail, E[D_i] <= 0 given earlier trails. Ordering trails by start, Azuma-Hoeffding gives
  P(max_n sum_{i<=n} D_i >= 2^112 - 2^109) < 2^-1384 (exponent
  2(2^112-2^109)^2/(K_TRAILS*(138*TRAIL_MAX)^2) bits, with the drift debits below already deducted). At most K trails
  are in progress at a check; they move the sum by at most 2*(123+15)*K*TRAIL_MAX.
* *Drifts.* Two terms are not centered: the non-canonical-prefix case adds at most 123*2^-31 per step in expectation
  (below 2^104 over the run) and the per-round cap increment floor(2^28*CORR) loses below 1 per batch round (below
  2^108). Both are debited from the 2^112 trail slack.

So the cap is exceeded only on these two events or when a round's eight draws all fail (below 2^-160 over all rounds):

    DELTA_CORR = 2^-1384 + 2^-1000 + 2^-160 + DELTA_TABLE <= 5.07e-14,

subtracted from the success bound. The independence used is that of the algorithm's own coins (trail starts and rho),
not of walk states; only the average H-trail is a premise.

### Common path charged by the anchor schedule

An anchor is set only in rounds t of the global batch round counter with `t % 16384 == 0` (ANCHOR_PERIOD=16384; `main`
runs `batch_step(t)` for t=1,2,..., and all walks of a batch share t). Anchor rounds run the same backward field block
with the anchor tail, whose common path has 15 instructions against the plain tail's 12 in the compiled kernel (the
period counter, loaded with the anchor by one `ldp`, tested (1) and incremented (1), and one `stp` stores the new
anchor and counter), and the same forward block. The skip tail omits the anchor load and comparison, two instructions
per copy against the plain tail. Measured on this revision's kernel at 2^16 walks and 2^20 adders, a plain round falls
by 11,814 instructions and an anchor round by 11,814 + (2^16-1), so with the one-instruction cap the anchor-minus-plain
gap is 3*(2^20-1)+9 at that batch; production uses 3*(2^28-1)+9. The anchor-schedule section gives the paired skip deltas.

Of a slot's rounds t=1..R, R = L+TRAIL_MAX, write a = R // 16384 and rest = R - 16384a. The schedule has a anchor rounds,
16334a + min(16334, rest) skip rounds, and the remaining rounds plain; these are exact counts, not expectations. Hence every
run's common path costs at most

    CONTROL = K*(35072*L + 35840*TRAIL_MAX)/2^28  (includes 34816 per round for blinding)
    COMMON = K*(SKIP_ROUNDS*SKIP_ARM + PLAIN_ROUNDS*PLAIN_ARM + ANCHOR_ROUNDS*ANCHOR_ARM)/BATCH + CONTROL

with SKIP_ARM=154,359,884,295, PLAIN_ARM=154,896,755,206 and ANCHOR_ARM=155,702,061,576 instructions per 2^28-walk round,
including the 1,024 allowance. C_COMMON = COMMON/S < 575.041797 (exact rational in `accounting.py`), plus the metered
corrections, CORR per processed step. Replays are unbatched and priced separately, so their per-trail anchors do not
enter this count. The orientation enters only through the table index and the sign bit, so the fixed path costs the
same for all coin values; only the corrections depend on the data, and they are metered.

### Batch size and per-round control amortization

The production batch is **B=2^28 lanes** (up from 2^20), the unit of simultaneous inversion. The per-lane work is unchanged: the block bodies are the measured 605,065,449/2^20 plain instructions per lane, and field, correction, tail and rare-path costs do not depend on B. Only the per-round charges are spread over 256 times as many lanes: the 256-instruction wrapper, the 1,024 allowance and the fixed 2^11 blinding-product allowance give (1,280+2,048)*(2^-20-2^-28) = 0.003169 fewer instructions per step (7.76e-6 bits). The 8*2^16-per-global-round draws are *not* amortized: a global round spans 2^4 batches, so RHO_COST = 8*2^16/2^4+2^11 = 34,816 per batch round and the per-step draw total 8*2^16/K is unchanged; CONTROL = K*(35,072*L + 35,840*TRAIL_MAX)/2^28. The batch buffer grows to about 328*2^28 bytes (2^36.4), far below the 2^76.49-byte record store, and peak memory stays 2^76.491741. Fewer rounds only tighten the correction-meter and Hoeffding bounds; the detailed arithmetic is in `unroll.md`.

### Production counter widths and control allowance

The C file is a finite measurement harness (signed `int` loops, a 64-bit round argument, a diagnostic `ring_n`); the
specified algorithm uses unsigned 128-bit batch-round and replay-length counters and 192-bit instruction meters and
work counters. The per-trail period counter `tper` is one 64-bit word: `measurement/step.c` increments it once per 16384
processed steps and compares it with CAP_PERIODS=2^52, so the cap TRAIL_MAX=16384*(2^52+1)=2^66+16384 steps is enforced through that
word, which holds at most CAP_PERIODS+1=2^52+1 (53 bits); the step count is never materialised. The anchor total
ceil(R/64)~2^89.5 is analytic, not a runtime counter, and the diagnostic escape ring is not part of the attack.

Charge 256 extra scalar instructions per full main-phase batch round (all B=2^28 lanes active) and retain 1,024 for
every finishing round; finishing, rare-event and replay bookkeeping keep their allowances. The production wrapper
stores completed t, L and R_MAX as two-limb values at offsets 0,16,32 of state x19, which an AAPCS call preserves.
This listing bounds recurring main control: 15 instructions (the two-limb round counter at [x19], its increment, the
`bl batch_step`, the loop branch and the two `b.hs` guards against the finishing and abort limbs). The correction
meter adds 19 more per round: the 192-bit meter M at [x19,#48] and the cap at [x19,#72] each take a three-limb add of
the block total x20 and of the per-round increment floor(B*CORR) at [x19,#96], then the three-limb compare against the
cap and `b.lo abort`. The normal path has 15+19 instructions. The full counter is kept in state; passing its low limb preserves t modulo 16384
and hence the anchor schedule, and the t bound below 2^97 excludes overflow. The wrapper does not touch the kernel or
its register allocation. With 32 register save/restore pairs, 6 stack/flags operations, 20 address/setup
instructions and 32 phase-handoff operations it is at most 124 with the meter lines, so 256 per main round covers it,
the final main comparison and transition included (L is positive). No per-lane finishing check is moved into this bound; the 1,024 finishing allowance and all event-only wide-counter prices are retained. 

The record counter reaches DP_MAX=2*ceil(S*theta)+2*K, below 2^68.492 records; the trail counter DP_MAX+2*K (K
starts, at most DP_MAX+K restarts); the table index spans C=2*DP_MAX, below 2^69.492 slots. Each is a 128-bit
integer whose control is charged inside the table operation. `accounting.py` charges DP_MAX table operations at
DP_COST=PROBE_CAP*PROBE_COST+DP_INDEX_COST+DP_COUNTER_COST+DP_OVERFLOW_COST
=2^10*2^11+2^11+8+8=2,099,216=2^21+2^11+16 each. A probe's two-limb increment and comparison cost at most 8 of its 2^11.
Once per table operation, the index reduction (a long division) costs at
most 2^11, the two-limb record-counter update 8 and the overflow test 8.

The home slot is (key(W)>>32) mod C. Under M1 the distinguished key is uniform on [0,2^196),
so its slot input is uniform on M=2^164 values. Its slot is within total variation
rem(C-rem)/(C*M) < 8.8116e-30 of uniform, rem=M mod C. As proved in the top-limb section,
the DP path skips the anchor comparison: the slot bits 32..195 are not otherwise read there,
except in exact class equality. Coupling each recorded class's slot to an independent uniform one
changes only the probe-cap event, charged DELTA_HASH=2*DP_MAX*(C+1)/2^164 < 2.9464e-8.
The K initial starts and
the DP_MAX+K restarts are charged K*START and (DP_MAX+K)*START,
START=255,184,896, so the
trail counter adds at most 8*(DP_MAX+2*K) against START*(DP_MAX+2*K), a
factor of 31,898,112 for each eight-instruction counter update. No counter is
incremented per processed step: `tper` advances once per 64 steps and the
batch-round counter once per round, inside the phase-specific allowances
stated above.

### The record table: layout, empty test and charged initialization

Slot i (0 <= i < C) lies in bank i>>40 at byte offset 128*(i mod 2^40) < 2^47, one 64-bit register; a probe forms
bank and offset from the two-limb index inside its 2^11. A slot is 16 64-bit words: w0..w3 the full key(W), w4..w7
coin a, w8..w11 coin b, w12..w13 the trail's live step count n (< 2^67, for the 65n replay cut-off: one walk step per
round, so n is the batch-round counter less the start round START writes into the slot state, a two-limb
subtraction), w14 the occupancy tag, w15 zero. **Empty test: w14 = 0**, a load and `cbz` inside the probe allowance.
An insert stores the 16 words with w14 = 1 inside its operation's charge, and nothing clears a slot, so a stored
record never tests empty; a probe finding w14 = 1 compares w0..w3 with the full key.

Before the walk each bank is zeroed, x0 its base, x1 its slot count (1 <=
x1 <= 2^40):

a loop of eight `stp xzr,xzr,[x0],#16` then `subs x1,x1,#1; b.ne`, exactly 10 instructions per slot,
plus at most 64 per bank to select it, set x0, x1 and advance the bank counter. With N_B=ceil(C/2^40)=754,917,021 banks,
`accounting.py` charges record_table_init = 10*C + 64*N_B =
8,300,400,423,641,169,372,144 instructions, below 2^62.847887 units, as
preprocessing, so it enters both totals. Loads and stores are
instructions in the cost model; no other transfer term applies.
Preprocessing is at most 2^89.961190 units; its rounding slack, over 2^59
units, also covers setup and the 2^40 table copy (16 stores per adder).

Every other array is first written by a charged step: the adder table by its accepted trial inside
A*2^40*(START+2^10), whose 2^14 constant covers negation and stores; a walk slot's state (512 bytes) by START
at each trail start; the batch buffers each round before they are read, inside PLAIN and ANCHOR; rare-path scratch
inside the forced maxima; and the counters and meters in setup, inside the 2^24 input_validation allowance.

## Concentration for the dominant event

This uses the **conditional** fresh-source premise, not a measured average rate. In the deferred-exposure filtration
of predecessor-proof.md Section 2.0 and Lemma E1(b), before the key and orientation of a globally fresh source are
exposed, its successor closes the two-cycle with probability p2=1/(2A). Stop the exposure process at the first
violation of the good-event conditions, padding with zeros after it and after termination: at most S trials, each
indicator I_t of conditional mean at most p2; already exposed sources belong to the inherited residual categories. The
law is **not** conditioned on G: the stopped process lives on the original space, counts every fresh-source entry on
G, and outside G the separately charged DELTA_CUT applies.

For lambda>0, conditional expectation gives

    E[exp(lambda*I_t) | F_(t-1)]
       <= 1+p2*(exp(lambda)-1)
       <= exp(p2*(exp(lambda)-1)).

Iterating this inequality and applying Markov to the exponential yields,
with mu=S*p2 and eta=2^-16,

    Pr[sum I_t > (1+eta)*mu]
       <= exp(-mu*((1+eta)*ln(1+eta)-eta))
       <= exp(-mu*eta^2/(2+eta)).

Only the adapted conditional bound is used; the last step is (1+u)ln(1+u)-u >= u^2/(2+u), u>=0. The limit is
n2=ceil((1+eta)*mu), the reserve B2=n2*wL; the accounting checks the exponent exactly as a rational above 100, so the
tail below 2^-100 is charged as delta2=2^-100. This uses M1/L-ori more strongly: reduced-order measurements do not
establish the conditional law at production lengths, and a marginal rate alone would not give the concentration, a
material limitation.

## Residual reserve and a single enforceable meter

Write c2=B2/S. The expected residual event cost per processed step before within-trail exceptional work has three
parts, only the first depending on A:

    E0 = (3/(4A^2)+2/A^3)*wL         longer entry, depends on A
       + 2^-63*32834*(w2+wL)         false filter, A-independent
       + (tau_exact*TRAIL_MAX/S)*wL merged step, A-independent
    tau_exact = T^2/(2N).

The longer-entry term is charged the derived bound above; the merged term retains its predecessor price. The
false-filter term carries a generous debit for speculative look-ahead: a false anchor check exposes at most 65 future
classes (the debit was sized for two filters), which the concentration proof cannot then treat as unexposed, so the
originating episode is charged its direct work, at most 66 later exposed-source steps and at most 2P-1 <= 32767 persistent
steps through an ensuing short cycle; the debit 32834*(w2+wL) exceeds that sum, and repeated charges of overlapping
exposures only overcount. An unclosed speculative path that returns to visited classes merges with another trail
(charged), closes a short cycle (this debit and the 8192-step bound) or violates G (failure); otherwise continuing past
its end exposes a new source, so no known look-ahead key is treated as fresh. The runtime meters the same calls.

For the residual success argument, stop the analytical residual process at the first G violation, as for the
two-cycle process. With Eres=E0 and kM=375000,

    B_res = kM*Eres*S,
    B_rare = B2+B_res.

The factor kM is a free analysis parameter: raising it multiplies the reserve kM*Eres per step and shrinks the
abort allowance 1/kM. At A=2^32, E0 is at most 7.452956e-9 (longer part 8.2201e-14, false-filter part 7.4525652572e-9,
merged part 3.0823e-13). Tau and kM are tuned jointly; this revision keeps kM=375000 and tau=4943073/10000000 (below),
for a reserve of 0.003030230 instructions per step and an abort allowance of 1/375000; the threshold is public
algorithm data.

All within-trail coincidences and excessive-length formal closures are excluded by G and paid through DELTA_CUT, not
by the per-step allowance or the predecessor's within-trail fixed-point argument. On G the stopped residual process is
the actual residual cost, and by the retained rate bounds plus the speculative-exposure debit its expectation is at
most Eres*S, so Markov gives a 1/kM = 1/375000 tail.

    Pr[rare work > B_rare and G] <= 2^-100 + 1/375000.

No independence between residual and two-cycle work is asserted. Rare operations are metered at their worst-case
prices, as in the predecessor, with an abort before any priced operation that would exceed B_rare, so work is bounded
even when the premises fail; the decomposition adds no runtime counter, freshness oracle or uncharged lookup, and its
parameters are public algorithm constants.

## Complete cost and success ledger

The ledger charges (COMMON+S*CORR+B_rare)/1000 for all walk steps, metered corrections and metered rare work, then
**adds**, rather than amortizes away:

* the correction meter's slack LAMBDA_AZ+LAMBDA_NF and one round's overshoot K*123;
* 13 instructions for each trail-ending batch exit, K_TRAILS*13 (anchor exits are inside wL);
* A adder constructions at 2^40 small-key trials each, and K initial trail starts, each trial or start with the
  complete capped scalar-sampling and two-scalar-multiplication price (trials also 2^10 for the key test);
* up to DP_MAX+K restart constructions, DP_MAX=2*ceil(S*theta)+2*K;
* DP_MAX table operations at 2^21+2^11+16 instructions each (probing, key comparisons, addressing, index
  reduction, counter update, overflow test), and the record table's zeroing, 10*C + 64*ceil(C/2^40) instructions;
* two bounded coefficient replays at the generous unbatched operation count; solving, verification, serialization,
  input validation and setup allowances;
* the Lemma F closure series: the first 64 terms are exact Fractions. For s>=64 and s<=A-2,
  the consecutive-term ratio is at most (e^2/8)*(1+3/A)^2<93/100: bound
  (A-s)^2/(A-s-1)<=A-s+2 and (s+1)(A-s+2)<=(A+3)^2/4.
  The remaining partial-support terms are bounded by term_64*(93/100)/(1-93/100).
  The full-support term TRAIL_MAX*A*2^(1-A)<2^-600 is separately included.
  The exact rational sum-plus-tail is below 2^-528; S<2^128 gives DELTA_LONG<=2^-400.
  No lgamma differences or floating guard supply an upper probability bound.

All processors' work is summed; no future-target amortization or free advice is used. The storage formula is
unchanged, 2*DP_MAX*128 + K*512 + 2*A*2*64 + 2^24 bytes, covering records, walk states and batch buffers, both adder
signs (and the block's copy) and a fixed reserve; at A=2^32 it evaluates to 2^76.491741 bytes, dominated by the
2*DP_MAX*128 record table. `accounting.py` evaluates these formulas and asserts that the total, rounded up at the
ninth decimal, equals the declared score 126728598631/1000000000.

The collision exponent penalty, re-derived at the bound with the predecessor's 2.6 margins, is
epsilon = 32836/A + 49248/A^2 + 131328/A^3 = 7.645227e-6 (section "The adder table"). The conditional success lower
bound is

    1-exp(-tau_exact*(1-epsilon))
      - 1/kM - 2^-100 - 2^-59 - DELTA_CORR - DELTA_ADDER_CLASS
      - DELTA_ANY_ABANDON - DELTA_CUT - DELTA_CAP - DELTA_PROBE - DELTA_HASH.

`probability.md` derives DELTA_ANY_ABANDON=13/10^7, covering any abandonment, including a nonresolving
selected merge. It allocates the 2^-59 bucket to report overflow (2^-60), degenerate first match (1/r)
exceptional adder hits (2^-95) and finite-index coupling (2^-80); their sum fits the bucket. No known DP-free prefix is resampled.
Tau is 4943073/10000000. The model term is approximately 0.3900043914313467; the recomputed lower
bound exceeds 0.3900001 (approximately 0.3900001438642747). 
Other failure charges are DELTA_CUT<=2.5144e-7, DELTA_CAP<=3.0195e-17, DELTA_ADDER_CLASS<=1.7748e-29,
DELTA_PROBE<=2.4924e-44, DELTA_HASH<=2.9464e-8 and DELTA_CORR<=4.8107e-49.

Per processed step the ledger charges C_COMMON=575.041797 (fixed path,
anchor schedule, control with the blinding draws, re-amortized at 2^28 lanes per batch), CORR=14.167091
(metered corrections) and the rare reserve 0.003030230 (two-cycle
0.000235371, residual 0.002794859; rounded up), in total at most 589.211918
instructions. 

## Evidence scope, execution and limitations

`evidence/anchor-schedule.json` records this revision's paired schedule counts, identical correction meters and
byte-identical one-round walk dumps; the 128-copy rows are in `unroll.md`. The A=2^32 source changes only the adder-index
mask and slot-shift immediates against the measured A=2^24 object (section "The adder table"), but the static finite
harness is not linked at A=2^32, so no A=2^32 execution is claimed. Finite counts/models do not prove production success. Pooled harness spares do not
price cold starts, which are charged separately. The reduced-order experiment and replays retain every-step anchor comparisons and test a superset of the
production skip schedule, not the skipped rounds themselves; the inherited toy runs the period-64 schedule because
its per-walk budget is far below 16384 steps, so the period change rests on the accounting.py coverage check and the
coverage theorem, not on the toy run.

combined.py reads only accounting.py and experiment.py from the credential-free verifier snapshot.
The inherited toy solves/replays with the old DP predicate, full-key dictionary, key-only adder rule,
no y window, a zero toy floor, and no production meter/probe cap/zeroing. Its successes/failures have
only that scope. Assembly checks are separate offline evidence. Production mixing, orientations, fresh-window domination, long-trail behavior, kernel prices, conditioned
adder H1 and H-trail/dy fold remain heuristic premises. 
