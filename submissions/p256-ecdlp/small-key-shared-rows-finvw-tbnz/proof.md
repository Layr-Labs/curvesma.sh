This revision adds one change, shared Montgomery rows, to the y3-fold revision (submission 205a498f, folder
`small-key-y-fold`, declared 126.789581). Two products with the same factor R use R's rows B1 = R*2^-64 and
B2 = R*2^-128 mod p, formed once (7 instructions each), and take 78 instructions instead of 90. So that both blocks
hold such a pair, the forward block also forms Q = dy*PRE into the pre slot, and the backward block forms
lambda = INV*Q and INV' = INV*D (sharing INV) in place of iw = INV*PRE, INV' = D*INV, dy and lambda = dy*iw; the
forward products Q and PRE' = PRE*D share PRE. Per walk step the backward block runs 110 fewer instructions and the
forward block 90 more, with 160 more per round in C: C_COMMON falls from 600.604538 to 580.577731 and the declared
score from 126.789581 to 126.741585. The walk, CORR, the anchor schedule, the rare paths and every success term
are unchanged; the meter counts carries of other integers at two sites (section "Shared rows").

Earlier revisions (items 7-11; mechanisms and evidence in the sections below, scores in "Current claim and provenance"): the y3-fold revision stores y3 by magnitude and draws adders with Montgomery y within
p/2 +- 2^-8 p, so dy borrows with probability at most 1/4 + 2^-16 given any history up to two steps back (CORR
17.709638 -> 13.959715; 126.793683 -> 126.789581); CORR was 22.468426 -> 17.458661 at the small-key revision (predecessor-proof.md "Small-key adders").

Lineage before 205a498f (content in `predecessor-proof.md`): ebc47628 (checkpoint-only on an earlier kernel), 961165d1 (eight copies), 1f942c41 (flag-carry), c2b498ec (dy pad 0), b75c69bb, 778628fb and earlier kernels, 7169b17e..8a577856 (trail cap), material predecessor b89d1d72.

# P-256 sigma rho with a 2^30-adder table and a concentrated two-cycle reserve

## Current claim and provenance

For a uniform scalar x in 0..r-1, given Q=[x]P on the protected NIST P-256 profile, the following classical
single-target algorithm has success probability at least 0.39 under the disclosed H1, M1/L-ori, fresh-class exposure
and H-trail premises. Its total charged computation is at most **2^126.741585** target-group-operation equivalents,
normalized at C=1000 scalar AArch64 instructions per unit, not one unit per P-256 addition. Peak memory is below 2^77
bytes, preprocessing (the rejection-sampled adder table and record-table zeroing included) is at most 2^87.961190
units and is included in the total, nonuniform advice is zero bytes, and instruction accounting retains H3.

The exact declared score is 25348317/200000 (126.741585), the ledger's total rounded up at the sixth decimal: a
conditional analytic bound, not an unconditional theorem or an executed full-size attack.

The immediate source is submission 205a498f-2dc0-4aa2-a78f-a6027f7c4e16 (`submissions/p256-ecdlp/small-key-y-fold`,
declared 126.789581, commit `6bf75c4419327f2d31f6b8b81aa6466f7ee46567`), then 21edf4bd
(`submissions/p256-ecdlp/small-key-paired-reduction`, declared 126.793683), then e2973980 (`submissions/p256-ecdlp/small-key-checkpoint-deferred-inv`, declared
126.851040), then ae463634 (`submissions/p256-ecdlp/small-key-checkpoint-only`, declared
126.852724), then 5a258d4a-928d-4f44-9338-ff9a64df7257 (declared 126.857022, `submissions/p256-ecdlp/small-key-adders`,
commit `01f9429a6c53200980f7d7561bc13893eae96b9f`), whose source is submission 70c0711e-ed9a-4be7-a5f6-03165ce1e6ae
(declared 126.868201, `submissions/p256-ecdlp/blinded-implicit-prev`, commit
`306bc0f8b8492eaacaee3a9f3059d94bb77d78f9`), with the lineage above. The
material predecessor is submission `b89d1d72-d9a1-4ef9-9ed7-c6ef0eec1ed3`,
immutable commit `0e1d63e9afd0fdb04b311b1ae72a2618271b371b`, directory
`submissions/p256-ecdlp/rho-sigma2`, whose algorithm, arguments and
measurements `predecessor-proof.md` reproduces with a historical banner. No
prior platform decision is an accepted premise.

**Precedence:** the claim is this document and `submission.json`. The predecessor's Sections 1.2--1.3 define the same
algorithm except for A, L and the rare-work threshold; 2.0--2.4 and 2.8 give the probability-space, escape, replay,
abandonment and merge arguments, re-evaluated here; 3--4 the forced-path pricing, resource caps and storage. Replaced
below: the batched step, its measurement, the record-table layout and zeroing, the two-filter rule of 1.3 step 5 and
all that depends on it ("Checkpoint-only cycle detection"), and the adder rule of 1.3 step 2 ("Small-key adders"; both at the end of `predecessor-proof.md`;
"The y3 fold and the y-window adder table"; the inverse chain's products, "Shared rows"). Its score and parameters are historical; its `evidence/accounting.json`
stays in the predecessor submission (omitted for the size limit).

## Algorithm, fixed costs and changed parameters

Use the protected p, a=-3, b, r and generator P unaltered, return zero for Q=identity and otherwise validate the curve
point. A class is {W,-W}, key(W) the integer Montgomery representation of its affine x-coordinate, canon(W) the
representative with even Montgomery y, and sigma(W) the sign taking canon(W) to W. The walk is

    W -> W + sigma(W) R[key(W) mod A].

Here A=2^30 adders, each R_i=[c_i]P+[d_i]Q with fresh uniform scalar coins, redrawn until key(R_i) < 2^231 and its
Montgomery y lies within p/2 +- 2^-8 p (at most 2^40 trials; "Small-key adders", "The y3 fold"); the table contains
both signs. K=2^32 walk slots are arranged into batches of 2^20 for simultaneous inversion. Distinguished points have
bits 30..89 of the key zero, so theta=2^-60 under the stated class model. Each trail starts from fresh independent a,b
with W=[a]P+[b]Q, ends at its first distinguished point, records it once with its coins, and a fresh trail starts; coefficients are reconstructed by replay only at a match.

Test ordering and replay are predecessor-proof.md Sections 1.2--1.3 with checkpoint-only cycle detection:
the distinguished-point test precedes the cycle test at every generated point, look-ahead and escape points included;
a detected cycle of length 2..64 escapes by doubling the canonical representative of its least-key member; the 64-bit
anchor filter only triggers exact checks. Every trail is limited to TRAIL_MAX=64*2^60+64 processed steps; the first
key match is resolved by two replays, each limited to 65 times its original live length. Degenerate matches and
exceeded caps fail the run, and each derived scalar is verified by [x]P=Q before output.

Set N=(r-1)/2+1 and

    tau_target = 4943092/10000000,
    T_min = ceil(sqrt(2*N*tau_target)),
    L = ceil(T_min/K), T=K*L,
    S = K*(L+TRAIL_MAX).

Every slot takes L main-phase steps and finishes its current trail, so S
bounds all processed steps including persistent steps and finishing tails. Sampling caps, maximum distinguished records, abandonment cap, probe cap,
restart costs, bounded replay, verification and output are unchanged formulas evaluated at this S; `accounting.py` recomputes the
point-order certificate proving ord_r(2)>2*TRAIL_MAX.

The common path costs at most C_COMMON=580.577731 instructions per processed step under H3 (anchor schedule below)
and corrections at most CORR=13.959715 more, metered (section "Blinded corrections and their meter"). A fresh-source two-cycle entry and a longer entry
each cost at most wL=2,021,787 more (an anchor check with its batch exit, look-ahead and escape); false filter episodes with their
exposure debt are charged 256*(w2+wL), w2=83,518 the former two-cycle check, merged steps wL, plus the
inherited forced-path and control allowances. The runtime uses **one total rare-work meter** and need not
recognize whether a point is globally fresh or classify observed events as fresh versus merged; that partition is
solely for analysis.

### The adder table at A=2^30 and the measured step price

Every entry rate and reserve term is a function of A. At A=2^30:

    fresh-source two-cycle rate  p2 = 1/(2A) = 2^-31,
    longer-entry rate            <= 3/(4A^2) + 2/A^3 = 6.505213e-19,
    conditional entry bound      P' = 64/(2A) = 2^-25,
    collision exponent penalty   epsilon = 196/A + 288/A^2 + 768/A^3 = 1.825392e-7.

The two-cycle rate is exact under L-ori. The longer-entry rate is at most 3/(4A^2)+2/A^3 (Lemma E1(c)):
E[w'] <= 3/(2A)+4/A^2 over windows of length 3..64, so E[w']/(2A) <= 3/(4A^2)+2/A^3, which the ledger charges and `accounting.py` asserts equals (3/(2A)+4/A^2)/(2A). Epsilon's closure margin charges every short-cycle entry, two-cycles included, 192 wasted states, with two-cycle term 2*192/(2A) and longer-entry term
2*192*(3/(4A^2)+2/A^3). The split-concentration threshold mu=S/(2A), the two-cycle reserve ceil((1+eta)*mu)*wL, E0 and
the reserves built on them are recomputed, not rescaled: the two-cycle charge is 0.0009415 and E0 is 5.972949e-11
instructions per step.

The A-dependence of the compiled step is checked on the binaries (flags as measured, B=2^20, ADDER_BITS=20, 24, 26 and 30,
disassembled): `batch_step` has 9,935 instructions at every A (11,272 before the shared rows), followed by five (A=2^20) or one (A=2^24, 2^26 and 2^30) alignment nops
that never execute; the two differ only in the `extr`, `ubfiz` and `and` immediates of the DP test and the adder
index (34 instructions: sixteen `extr` and sixteen `ubfiz` in the sixteen block step tails, one `and` and one `extr` in
the compiled C step tail); `step_rare` (98 instructions; three trailing padding nops at A=2^24, 2^26 and 2^30), `new_trail`,
`add_single`, `finv` and `dp_restart` are identical up to immediates, and `rare_anchor` has one more instruction at
A=2^24, 2^26 and 2^30, inside the per-event OVERHEAD. The A=2^30 whole-binary mnemonic sequence equals
A=2^24's instruction for instruction (122,403 instructions, immediates aside), so the step price
is A-invariant. The 64-bit anchor filter compares the whole low limb, so its false-match
charge (2^-63 per step, twice a 64-bit filter's 2^-64) is independent of A; the table is data.

The larger table is charged, not assumed free: construction (worst case of the small-key rejection sampling)
A*2^40*(START+2^10) = 2^87.961189 units; storage A*2*64 = 2^37 bytes, twice with the block's copy at 2^40, inside the
peak-memory formula; each 128-byte record's read-back (16 words) inside START's 2^14 constant. With the
record-table zeroing below, preprocessing is at most 2^87.961190 units, included in W.

### The batched step in this revision

The walk, its distinguished-point test and every rare-path primitive are the predecessor's. The implementation of the common batched step changes
in the value-preserving ways 1-6 and 8-10, cycle detection is checkpoint-only (item 7), y3 is stored by magnitude (item 11) and the
inverse chain shares rows (item 12):

1. **Branchy corrections.** Every conditional correction of the field
   arithmetic, computed branch-free and always executed in the predecessor,
   is now a forward branch over its fix-up, and a taken fix-up adds its own
   extra instruction count to a meter (below).
2. **A register-resident forward block.** The prefix products PRE[w+1] =
   D[w]*PRE[w] (and Q[w] = dy[w]*PRE[w], item 12), w=1..B-2, run as one
   assembly loop (`FLOOP`) with PRE in x19..x22, never reloaded.
3. **Operand forms.** inv is canonical whenever read (item 9); PRE, Q and
   lambda stay below 2^256 (their consumers accept that) and the products'
   last steps write straight into their destination registers.
4. **2^20-walk batches** and a copy of the two-sign table T2 at virtual address
   2^40, read by the blocks, so the adder address is one `orr` and x27 is
   free for the meter; the C code keeps the static T2.
5. **A blinded prefix chain.** Each global round draws one fresh uniform
   unit rho shared by all batches (rejection sampling, at most 8 draws) and
   every batch starts its chain at PRE[1] = rho*D[0]; the C path for w=0
   takes inv*rho = 1/D[0]. Inverses, slopes and points are unchanged,
   and every prefix PRE[w] is a uniform unit given the walk, which the
   meter's abort bound uses.

6. **Eight copies per loop.** Both blocks hold eight labelled copies; the
   end test (`tst x23`/`b.ne`, or `subs x26,#8`/`b.ne`) runs once per eight
   steps, copy j processing lanes w = -j mod 8. The backward block enters
   through a `tbz`/`tbnz` tree on bits 6-8 of x23, the pointer to 64-byte DP
   record w (array aligned to its size): 3 instructions at w = B-1, at most
   4 on a re-entry after a rare exit (two-copy entry: 1); its end test
   follows the w = 1 mod 8 copy. The forward block enters at copy 2 with
   counter B, running 6 + 8(B/8-1) = B-2 steps. Every walk is processed as
   before (B a multiple of 8, asserted at compile time).
7. **Checkpoint-only cycle detection** (ae463634). The two-steps-back
   filter is gone; the anchor filter alone catches every short cycle
   (section at the end of `predecessor-proof.md`). The lane record {sel, anchor, tper} keeps the 32-byte
   stride of the implicit-previous-key revision.
8. **Small-key adders, reversed D and t** (small-key revision; section "Small-key adders", moved to the end of
   `predecessor-proof.md`).
9. **Deferred canonical test** (e2973980): inv's test r3 >= p3 runs
   in the step tail (section at the end of `predecessor-proof.md`).
10. **Paired reduction words** (paired-reduction revision; section below).
11. **The y3 fold and y-window adders** (y3-fold revision; section below).
12. **Shared rows** (previous revision; section below).
13. **Windowed batch inversion and a one-instruction anchor cap** (this revision; section below).

`measurement/step.c` is the predecessor's harness with these changes; the backward block `BLOOP_FIELD_B`
(`measurement/bloop_branchy.h`) and the forward block `FLOOP` (`measurement/floop.h`) are generated by `gen10.py` from
the predecessor's measured `BLOOP_FIELD` text (the generator stays in 70c0711e, `measurement/gen10.py`,
omitted for the size limit; this revision's blocks are its output after 5a258d4a's `measurement/gen11.py`,
`measurement/gen12.py` (e2973980), `measurement/gen48.py` (21edf4bd) and `measurement/gen13.py` (205a498f), all omitted
likewise, and `measurement/gen49.py`): gen48.py re-formed the reduction words (item 10), gen49.py replaces the forward
block and the backward block's first products (item 12), and the corrections are as before. Item 7 changed the step tails, the rare exits, `step_rare`, `new_trail` and `main` and deleted `rare_2cycle`
(section "Checkpoint-only cycle detection"); item 9 changed the inv finish, the tails' last branch and the rare exits; item 10 changed only the product rows;
item 11 changed the end of the backward block; item 12 edits `step.c` (section "Shared rows"); item 13 edits only `finv` and the 2^60 anchor cap test in `step.c` (section "Windowed inversion and the one-instruction anchor cap").

### Windowed inversion and the one-instruction anchor cap (this revision)

**Windowed batch inversion.** `batch_step` performs one inversion per global round,
`finv(&S[0], &acc)`, on the product of the lane prefixes. The inherited `finv` is binary
square-and-multiply over the 256-bit exponent p-2, whose top bit is absorbed by starting at a (255 squarings, 127 multiplications). This
revision walks the same exponent in 4-bit windows from the top: the top nibble of p-2 is 15, so
the accumulator starts at a^15; the table a^1, a^3, ..., a^15 costs one squaring and seven
multiplications; each of the 63 lower nibbles costs four squarings and, when nonzero (32 are),
one multiplication. That is 253 squarings and 39 multiplications against the binary method's 255 and 127: the
same field operations, re-associated. `finv` is a single shared function, so every caller gets
the windowed form. The once-per-round `batch_step` call is the charged fixed-path change; the
rare-path callers (`add_single`, and through it `f_single`, `escape`, `step_rare`,
`rare_anchor`) become cheaper, never more expensive, so the forced-path maxima measured with
the binary `finv` (F_MAX=41190, ADD_MAX=41185) remain upper bounds and the rare-path charge is
conservative. The other rare primitives (`new_trail`, `dp_restart`) are unchanged; only `finv` and
the cap test change. Value-preserving: the nibbles are recomputed from
p-2 and base and patched walk dumps are byte-identical.

**One-instruction anchor cap test.** The anchor tail tests `tper >= 2^60`; the production build
uses one `tbnz x12,#60` instead of `lsr`+`cbnz`. Same predicate on every reachable value: the
stored counter never exceeds 2^60+1, and every value in [2^60, 2^61) has bit 60 set. The generic
non-2^60 build keeps `cmp`/`b.hs`; rare-tail prices are unchanged.

**Measurement.** Organizer image (gcc 14.2.0-19, valgrind 3.24.0, arm64, `-O3
-fno-tree-vectorize -mgeneral-regs-only -static`), 2^16 adders/walks/spares, one round each.
Callgrind Ir: plain 39,402,364 -> 39,390,550; anchor 39,665,371 -> 39,588,022; meters, counters
and walk dumps equal. The inversion is per round (B-independent); the cap test is per lane per
anchor round (B-1). Hence `PLAIN_FIXED` = 608,739,081 - 11,814 = 608,727,267, `ANCHOR_EXTRA` =
3*(2^20-1)+9 = 3,145,734, C_COMMON = 580.577731 and the step charge 594.538411 (from
594.650484).

**Rare-path repricing.** Every unbatched addition (`add_single` and its callers) calls
the shared `finv`, whose count is data-independent, so each is exactly 11,814 cheaper
(`single` mode 40,884 -> 29,070 Ir). F_MAX 41,190 -> 29,376; wL (65 F_MAX terms; ESC
unchanged) 2,789,697 -> 2,021,787. The reserve scales with wL, so the two-cycle charge is 0.0009415 at A=2^30 (E0
5.972949335e-11).

### The corrections

With x28 = p1 = 2^32-1 and x29 = p3 = 2^64-2^32+1 set at block entry (p = (2^64-1) + p1*2^64 + p3*2^192), adding p
to four limbs is `subs r0,r0,#1; adcs r1,r1,x28; adcs r2,r2,xzr; adc r3,r3,x29` and subtracting it is
`adds r0,r0,#1; sbcs r1,r1,x28; sbcs r2,r2,xzr; sbc r3,r3,x29`.

* **Subtraction** a-b of canonical a, b (six per walk step: D, dy, two in x3, t, y3): `subs/sbcs` x4, then `b.cs`
  over [add p; meter += 5]. The common path is 5 instructions (was 11; y3's 7 with the fold), a borrow costs 5 more, and the result is a-b
  or a-b+p, canonical.
* **Carry-only finish** (PRE, Q, lambda, inv): T = r + c*2^256, c in {0,1}, T < 2^256+p, c the carry flag left by the
  product's top `adcs`. `b.cc` over [subtract p; meter += 5]: 1 instruction (was 7 or 9), 5 more when c=1, giving
  T - c*p < 2^256.
* **Canonical finish** (y3' = t*lambda, and inv before item 9; T < 2p): `b.cc L` (c=0), else subtract p, meter += 4, `b done`;
  `L: cmn r3,x28; b.hs slow`. Common path 3 (was 9), c=1 costs 4 more. If c=0, T < p unless r3 >= p3 (`cmn` sets C
  exactly then); the out-of-line slow block computes T-p, keeps it with four `csel` on no borrow and adds 10.
* **lambda^2**, from lambda < 2^256: T < 2^256+p, so after a c=1 subtraction T-p may still be at least p; the fix-up
  (meter += 5) falls through into the same `cmn/b.hs` re-test, and one more conditional subtraction in the slow block
  suffices (T - 2p < p).

### Flag-carry finishes (from 1f942c41)

The product finishes branch on the carry flag C: c2b498ec ended each product with `adc c,xzr,xzr` (c = x13, x11 at
lambda^2) and `cbz c,L`, and 1f942c41 made each pair one `b.cc`, which branches exactly when C = 0; that keeps every
outcome if the branch follows the `adcs` whose carry is c and nothing later reads c (x13 and x11 are dead at the
square's and y3''s sites; gen49.py's products never formed a carry register). On this revision's four builds (B=2^20
at A=2^20 and 2^24, B=2^16 at A=2^16 and 2^20) all 80 `b.cc` sites (two per forward copy, four per plain and anchor
backward copy) directly follow the top `adcs`. Rows' domains: "Shared rows". Every exposed value (x3, y3, key,
sign bit) is canonical.

**Exactness checks** (5a258d4a's `evidence/measurement/branchy.json`, keys `previous_revision` and earlier,
omitted here for the size limit, quoted values in `ancestor_quotes`; on the builds that introduced these changes,
all with the two-steps-back filter; reruns: sections "Deferred canonical test", "Paired reduction words", "The y3 fold" and "Shared rows"): (i) n=8 dumps of 2^16 walks on
two input sets and stress builds with 2^4, 2^6 and 2^10 adders equal the predecessor kernel's byte for byte; (ii) the
backward block equals an exact big-integer model (Montgomery integers, branch outcomes, meter) on 240,000 random
and crafted inputs; (iii) the forward block on 6,000 runs; all three rerun for 1f942c41, 0 mismatches.

### Eight copies per loop (eight-copy revision)

Copies equaled 1f942c41's instruction for instruction (forward 119, plain backward 570, anchor backward 575, slow
sets 30); `batch_step` differed only in the entry tree, the counter immediate and nops. Loop control per round
fell by exactly 3B/2-5 (1,572,859 at B=2^20); a re-entry after a rare exit costs at most 3, inside OVERHEAD.

### Paired reduction words (paired-reduction revision)

A product computes T = (a*b + M*p)/2^256 in four rows: row j adds a*b_j, takes m_j, the accumulator's low word
(-1/p = 1 mod 2^64), and adds m_j*(p+1)/2^64 = m_j*2^32 + m_j*p3*2^128 from the next word up. The paired-reduction
kernel forms m_j*p3 by `mul`, `umulh` of m_j and x29 = p3 and joins m_j>>32 and m_(j+1)<<32 (disjoint bits of one
word) by one `extr`: a product fell from 94 to 90 instructions and the square from 75 to 70, 21 fewer per backward
block step and 4 per forward step (21edf4bd's `measurement/gen48.py`, omitted here). Each row adds the same integers
in another grouping, so T, the carry and every finish and meter outcome were unchanged; the five accumulator words
stay below p*2^64 + p + 2^128 < 2^320 (asserted in `accounting.py`). Evidence (21edf4bd's and 205a498f's
`branchy.json`, `paired_rows`): `batch_step` -368 on five build pairs, 12,000,000 + 2,400,000 exact-model block cases,
four one-instruction mutants each failing 96,001 to 97,622 of 100,000, rounds falling by 21(B-1) + 4(B-2) on both
input sets, dumps equal. Here only the square and y3' keep this form, their text unchanged (from `mul x1,x0,x2` on,
the backward block is the y3-fold revision's).

### The y3 fold and the y-window adder table (y3-fold revision)

**Kernel change.** The backward block ends with y3 = y3' - y(W), y3' = t*lambda and y(W) canonical, and a `b.cs`
borrow fix-up (extra 5). 205a498f's `measurement/gen13.py` (omitted here) derived that revision's `bloop_branchy.h`
from the paired-reduction revision's by one edit and its header line, each asserted to match once: `cmp x30,x11` compares the top limbs y3'[3] and y(W)[3]; when y3'[3]
is smaller a four-limb chain stores y(W) - y3' = -y3, which cannot borrow, and branches past the old code; otherwise
(`b.hs`) the old y3' - y(W) and its fix-up run, the fix-up firing only at equal top limbs (slow site `y3_tie`). The
stored y is y3 or p - y3: |y3' - y(W)|, or y3 at a tie. Each of the 16 backward copies gained these 7 instructions (`batch_step` +112
at all four shapes); either path runs 2 more instructions per block step.

**The walk is unchanged up to sign.** The stored y feeds only the sign bit (the parity of its Montgomery value, which
selects R_i or -R_i and flips with the negation, p being odd) and the next step's y(W). A point stored as -W' selects
the negated adder and yields the negated old successor, so every class, key, adder index, distinguished-point test
and rare-path outcome is the paired-reduction kernel's; x3, lambda, t and PRE are unchanged, and M1, L-ori and the entry
rates, which concern the class map, are untouched. Dumps (205a498f's `branchy.json`, `y_fold`): the fold
built on the checkpoint-only kernel ends every walk at that kernel's point or its negation (flipped in 0.4858 to
0.5146 of the walks), event counters and record digests equal, in 16 pairs (8 to 600 rounds, 2^4 to 2^16 adders,
rare paths forced); the y3-fold kernel's equal the fold on the deferred-inv kernel's in 6 more; mutants of the edit
move 65,145 to 65,338 of 65,536 final points out of class; its block equals the fold block on 12,000,000 random and
crafted inputs.

**The adder rule.** Step 2 of the procedure accepts R when R != O, key(R) < 2^231 (the top limb of the Montgomery x
below 2^39) and Y_LO <= Y <= Y_HI for its Montgomery y Y, with Y_LO = (p - (p>>7))/2 and Y_HI = p - Y_LO, i.e.
|2Y - p| <= p>>7: two compares inside each trial's 2^10. This replaces the test key < 2^224 of "Small-key adders",
whose arguments hold with S the set of logarithms s for which [s]P passes. As 2Y_LO - p = -(p>>7) and
2Y_HI - p = p>>7, the window is closed under Y -> p - Y, so S = -S as the accidental-relation argument needs, and
u = Y/p satisfies |u - 1/2| <= U_DEV <= 2^-8 (`accounting.py` asserts both). D and the second x3 subtraction now
borrow with probability below 2^231/p < 2^-24.9999 for every adder, charged without pad.

**Acceptance count** (asserted in `accounting.py`). With I = [0, 2^231) and J = [Y_LO, Y_HI], |J| < 2^249, the 2-D
finite Fourier count of affine points with Montgomery coordinates in I x J has main term |I||J|(r-1)/p^2 >
2^224 - 2^94.07. For (a,b) != 0, f = aRx + bRy has its only pole at O (order at most 3), so Bombieri's bound gives
|sum_P e(f(P)/p)| <= 6 sqrt(p); with sum_{a != 0} |Ihat(a)| <= pL, L = 1 + ln p, and the same for J, the error is at
most BOX_ERR = 6(isqrt(p)+1)(2L + L^2) < 2^145.5598 < 2^160/2^14. So at least M_ACC = 2^224 - 2^160 - 3 points pass,
as before, and Q_ACCEPT_LOWER, the trial cap, LOG2_ADDER_EXHAUST, DELTA_LOOP and the 2^87.961189-unit table charge
are unchanged; a constant 2^14 times Bombieri's would still suffice.

**dy under the fold model.** dy = y(R) - y(W) (`ldp x3,x4,[x7,#32]`, the selected adder's y, minus `[x24]`, the
stored y) borrows iff y(W) > y(R). Let v = y(R)/p for the adder and sign the step selects, so v = u or 1 - u. Model
(support for H-trail, as uniform operands are for the other sites): each step's y3' is uniform on [0,p) given the
past, its parity independent of its magnitude. A stored y3 or p - y3 then has a fair-coin parity independent of every
magnitude, so v = u or 1 - u with probability 1/2 each; and given everything through the stored y two steps back
(c*p, m = min(c, 1-c)), the stored y one step back has density 2 on [0,m) and 1 on [m, 1-m] in units of p, so
y(W), its fold with a fresh y3', exceeds v with probability (1-v)^2 for every c. The rate is ((1-u)^2 + u^2)/2 = 1/4 + (u-1/2)^2 <= 1/4 + 2^-16 for every
such history (not given the stored y one step back) and every adder, not an average over the table; `accounting.py` asserts 1/4 + U_DEV^2 <= q_dy = 1/4 + 2^-16 and charges q_dy + 2^-10. With c
the stored y one step back, the model gives E[dy/p] = 1/2 - m^2 when m <= min(v, 1-v) and at most 1/4 + 2^-8 otherwise, so
in that revision lambda = dy*iw (iw a uniform unit given the walk) kept its fifth-word bound 1/4 given any
history, measured 0.208; here Q = PRE*dy has that form ("Shared rows"). The law needs two fold outputs with fresh y3': a trail's first two steps
(the start point's y is not one) and the steps around a short-cycle entry (up to 128 persistent steps repeat values;
the escape point's y is canonical) are charged the full 5, DY_EXCESS = 5(2K_TRAILS/S + 130(p2 + 3/(4A^2) + 2/A^3)) <
1.938e-5 per step, below the dy pad 5*2^-10 (both asserted). Walk 0 takes the compiled C step every round and lane B-1's dy is a C
subtraction ("Shared rows"), both branch-free and unmetered.

**The y3 charge.** The fix-up left at y3 fires only if y3'[3] = y(W)[3], for a uniform canonical y3' probability at
most 2^192/p < 2^-63.99 whatever y(W); it is charged 2^-24 with weight 5 in CORR_SLOW beside the three slow blocks.
The y3 site (5 at 1/2 + 2^-10) leaves CORR_SITES and dy moves from 1/2 without pad to 1/4 + 2^-16 + 2^-10. C_MAX stays
86 (asserted; y3_tie's 5 replaces y3's); CORR falls from 17.709638 to 13.959715 and ceil(CORR) from 18 to 14.

**Measured** (205a498f's `branchy.json`, `y_fold`): FIXED +2(B-1) against the paired-reduction kernel in every
counted round, both input sets, plain and anchor (2,097,150 at 2^20); over 127 further rounds of the fold built on the
checkpoint-only kernel, -48 and +144 from net escape negations (24 instructions each; rare-path work at its worst
case). Site rates of that build (rounds 2..8, two tables): dy 0.25008 and 0.25024 (checkpoint-only 0.50009), lam
0.20791 and 0.20861, y3 0; rounds 2..64 on four tables: dy 0.249782 to 0.250142, lam 0.208314 to 0.208426. Binned by
octile of the stored y two steps back, dy stays at 0.2467 to 0.2526 and 0.2456 to 0.2543, flat at the model's 0.250005
(largest deviation about 2.2 binomial standard errors); binned by the stored y one step back it follows |1/2 - beta|
(0.4385 down to 0.0607). The y3-fold kernel's meter: 13.702 to 13.715 per walk step (rounds 2..8 and 2..128, two
tables), below CORR; this revision's: "Shared rows".

### Shared rows (this revision)

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
kernel's by p, enters only the square and y3', whose finishes are canonical; INV' is canonicalised before it is read;
Q and PRE' enter only as multipliers or rows (any value below 2^256) and PRE(B-1) the inverse. x3, y3, keys, signs and
each step's INV equal the y3-fold kernel's; no exit is added and C_MAX stays 86 (asserted).

**Corrections.** Every site keeps its weight and idealized rate (site q replaces iw), so CORR is unchanged: Q = PRE*dy
has fifth-word probability at most E[PRE]*dy/2^512 < dy/(2p) given dy (up to 1/2 as dy nears p), so at most
E[dy/p]/2 <= 1/4 given any history up to two steps back, as E[dy/p] <= 1/2 under the fold model ("The y3 fold",
the bound lambda = dy*iw had); the steps DY_EXCESS names add below 5e-6 per step to q (5*(1/2-1/4) times their
share DY_EXCESS/5), inside its 2^-10 pad; lambda = INV*Q has two blinded operands,
as iw = INV*PRE had; k = 1 needs N < L, below 2^-62.99 under the uniform-multiplier model, inside each site's 2^-10
pad. The held PRE integer is pi or pi + p as before, so the 86*2^-31 term stands. Q, lambda and INV' are other
integers than the y3-fold kernel's iw, lambda and INV', so the meter is not that kernel's: rates per walk step over
rounds 2..8 on two tables are q 0.20789/0.20847 (the old lambda law), lambda 0.24930/0.24954, dy 0.25008/0.25023, inv
0.24977/0.24966, PRE 0.25103/0.24977 (y3-fold: iw 0.24972/0.24986, lambda 0.20791/0.20861), and the totals stay
within noise of it (below).

**Measured** (`branchy.json`, `shared_rows`; organizer image). Static: 80 `b.cc` sites, each right after its top
`adcs`; 505 or 506 memory operands per build encode; the new code reads only x19..x29, memory and registers it wrote
first; every other own function compiles as before (one alignment nop). The generated text equals the exact model on
6,000 backward-pair and 6,000 forward-copy cases (PRE >= p and solved k = 1 events included). One forward then one
backward copy (tail's inv test, Lz and exit) of each kernel equal their exact models, and each other's walk, on
2,400,000 random and crafted cases (k = 1 at all four shared-row sites; lambda, PRE' and raw INV' representative
classes, 23,637 to 31,584 each per 600,000): 0 mismatches; six one-instruction mutants (row shift, `extr` amounts, a
scalar offset, the Q store, a finish branch) move the walk in 87,899 to 100,000 of 100,000 cases per seed. Rounds: at
B=A=2^16 a round falls by 1,310,630 = 110(2^16-1) - 90(2^16-2) - 160, and 8 rounds by 8 times that, on a y-window
table and both uniform input sets; 160 is Q(B-1)'s `k_sub` (18) and `k_mul` (115) and 27 compiled stores and moves.
Dumps of ten configurations (8 to 600 rounds, 2^4 to 2^16 adders, small DP masks and one-period caps forcing escapes,
restarts and abandonment; two tables) equal the y3-fold kernel's byte for byte, event counters included; only
fix_meter differs. Meter per walk step with y-window adders, rounds 2..8 and 2..128: 13.711969 and 13.703247 (one table),
13.704984 and 13.711226 (the other), against 13.714181, 13.702654, 13.707284 and 13.710572 for
the y3-fold kernel; 2^20 rounds: next
section.

### Measured fixed path and metered corrections

Counts are exact guest instructions under QEMU 8.2.2 user mode over the harness's functions main..read_fe and
add_single..the end of batch_step (identical to `-one-insn-per-tb` counting; a 1,000-iteration three-instruction loop
adds exactly 3,000; Valgrind 3.22.0 Callgrind reproduced round differences on three earlier kernels, built with the
Ubuntu GCC 14.2.0 cross build and the organizer flags `-O3 -fno-tree-vectorize -mgeneral-regs-only`). A
METHOD_ALLOWANCE of 1,024 instructions per round is charged on every measured round for method and toolchain
differences.

**This revision's kernel** was built in the organizer's reference-count
image (GCC 14.2.0-19, same flags, `-static`) and counted by the same
QEMU method, filter ranges from each binary's symbol table. Every taken
correction adds exactly its extra instructions to the meter, so a round
without a rare event costs exactly FIXED + (meter increase); an escape
adds its rare path, priced by the rare reserve. Measured (round =
cumulative count at n less at n-1 rounds; FIXED = round minus meter
increase; clean = no escape in the round):

    B=2^20, A=2^20, input set 2, round 5 (clean)        608,739,053
    same build, input set 1, round 7 (clean)            608,739,053
    y3-fold kernel, same rounds                         629,710,483
    paired-reduction kernel, same rounds                627,613,333
    deferred-inv kernel, same rounds                    653,827,704
    checkpoint-only kernel, same rounds                 654,876,279
    1f942c41's kernel, same image, set 2, round 5       659,594,869
    B=2^16, A=2^20, plain round 2, both sets            38,084,331
    same build, anchor round 64 (clean), both sets      38,346,478
    y3-fold kernel, both sets                 39,394,961, 39,657,108
    paired-reduction kernel, both sets        39,263,891, 39,526,038

No counted round has an escape or unbatched operation. Against the y3-fold kernel a round falls by exactly
110(B-1) - 90(B-2) - 160, 20,971,430 at B=2^20 and 1,310,630 at B=2^16 in plain and anchor rounds, on both input sets
(section "Shared rows"); that kernel rose by 2(B-1) against the paired-reduction kernel, which fell by 21(B-1) +
4(B-2) against the deferred-inv kernel, which fell by B-1 against the checkpoint-only kernel (each batch exit gives
back 9, charged separately). 1f942c41 charged 659,594,897, 28 above this image's count of its round, and each later
revision charged its predecessor's charge less its exactly counted change: eight-copy 3B/2-5 =
1,572,859, implicit-previous-key its measured 1,048,577 (656,973,461, which the small-key revision
kept), checkpoint-only 2,097,154, deferred-inv 1,048,575, paired-reduction
26,214,371 (627,613,361), y3-fold -2,097,150 (629,710,511), keeping the offset of 28:

    PLAIN_FIXED = 629,710,511 - 20,971,430 - 11,814 = 608,727,267 instructions,
    per walk step 580.527561 (rounded up),

so with the allowance PLAIN = (608,727,267+1,024)/2^20 = 580.528537 per walk step (rounded up), against 600.539718 in
the y3-fold revision, 598.539720 in the paired-reduction revision, 629.039690 in 1f942c41 and 714.613617 in the
predecessor. On these counted rounds (uniform adders) the meter lies between 19.341 and 19.350 per walk step (y3
never fires; dy borrows at the uniform-table rate 1/3); y-window meters: section "Shared rows".
Anchor rounds: section "Common path charged by the anchor schedule".

**Callgrind.** The organizer image's Valgrind 3.24.0 Callgrind (`--collect-atstart=no`) on 1f942c41's B=2^16,
A=2^16 build, input set 1, gives Ir = 42,730,422 for one round and 85,459,279 for two: round 2 costs 42,728,857,
the QEMU count of that binary (90,122,446 - 47,393,589).

### Blinded corrections and their meter

Each site j has extra cost w_j when taken and an idealized probability q_j for a uniform point and prefix: borrow 1/2 for the first x3 subtraction and t; below 2^231/p for D and the second x3 subtraction (small-key adders);
1/4 + 2^-16 for dy (history-conditional bound derived and measured in "The y3 fold and the y-window adder table");
fifth word set 1/4 for the products (PRE, Q, inv, lambda, y3') and 1/3 for the square; below 2^-31 for each slow block
(r3 >= p3 needs T in an interval of length below 2^224 + 2^192, and T has density below 1/p) and below 2^-24 for the
y3 tie. The walk step is charged

    CORR = sum_{j not in Z} w_j (q_j + 2^-10) + sum_{j in Z} w_j q_j + 36*2^-24 + 86*2^-31 = 13.959715 instructions,

Z = {D, second x3 subtraction}, w_inv = 5 and the inv slow test 11 (section "Deferred canonical test"), dy at
1/4 + 2^-16 and y3 replaced by y3_tie (section "The y3 fold"), against 13.703 to 13.712 per walk step measured on this kernel
with y-window adders (`branchy.json`, `shared_rows`), below CORR.

**What a lane's count depends on.** Let pi = PRE[w] mod p, the lane's prefix residue. With inv canonical when read, everything
a lane computes in a round is a function of its walk point W (with its adder) and of the integer held for PRE[w],
which is pi unless it is pi + p: dy, Q, D, PRE[w+1], inv = 1/(pi*D) up to a fixed factor, lambda, the square, x3, t,
y3' and y3 (lane B-1 forms dy and Q in C, unmetered). Write c(pi, W) in [0,86] for the lane's correction count when PRE[w] is held as pi. The value pi + p
(below 2^256) needs pi < 2^256 - p < 2^224, so for a uniform pi it occurs with probability below 2^-31 and changes the
count by at most 86, covered by the charge 86*2^-31 per step. The walk (x3, y3, keys, signs) does not depend on pi
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
product exceeds the exact integer by p only when N < L, below 2^-62.99 under the same model, inside each site's 2^-10
pad; the slow blocks are below 2^-31.
Along actual walks the measured per-step mean is 13.703 to 13.712 with y-window adders (13.702 to 13.715 on the
y3-fold kernel), 17.662485
to 17.667348 on the deferred-inv kernel with small-key adders,
17.40 to 17.43 on the small-key kernel (5a258d4a's smallkey.json) and 22.409 to 22.422 on 70c0711e's with uniform adders. H-trail bounds only an abort probability,
never the work.

**Meter and cap.** The blocks add each taken fix-up's extra count to x27;
each round adds the block totals to one 192-bit run-wide meter M, inside
the per-round control allowance. After round k, with n_k processed steps
so far, the run aborts if

    M > n_k*CORR + LAMBDA_CAP,
    LAMBDA_CAP = 2^112 + 2*(86+14)*K*TRAIL_MAX + 2^96.

So total correction work is at most S*CORR + LAMBDA_CAP plus one round's
overshoot K*86 (86 = most correction instructions in one step), all charged: CORR inside the step charge, the rest (below 2^102.052 units) as
`correction_meter_slack`; the blinding draws and products are charged in
the per-round control (8*2^16 per global round spread over the 2^12
batches, plus 2^11 per batch round). Work is bounded
even if H-trail fails.

**Abort probability without a walk model.** Write M = sum_s g(W_s) + sum_r N_r, s over processed batch steps and
N_r = sum over round r's lanes of (c(pi, W) - g(W)), round r's blinding noise.

* *Blinding noise.* The walk, hence every W and every PRE up to the factor rho, is independent of the rho draws.
  Given the whole walk, each global round's rho is a fresh uniform unit, so each lane's pi = rho*P[w] is uniform and
  E[N_r | walk] <= 2^20*86*2^-31 (the non-canonical prefix case); distinct rounds have independent rho, and one
  round's noise over all batches is at most 2^32*86. Hoeffding's maximal inequality over at most 2^95.492 rounds gives
  P(max partial sum >= 2^96) < 2^-1000.
* *Trails.* A trail's batch steps, and so its sum T_i of g, are a deterministic function of its start coins and the
  fixed table (the walk, its escapes, its distinguished-point and cap rules). Trails start from fresh, independent
  coins, so given the table (and each trail's slot and starting round, which fix its lane and anchor phase) T_i and
  its length len_i are independent of all earlier trails, whatever the walk's internal dependence; merges and cycles
  do not change this. Every slot finishes its current trail, a trail has at most TRAIL_MAX steps, and at most
  K_TRAILS = DP_MAX + 2K trails start, so D_i = T_i - CORR*len_i lies in an interval of width (86+14)*TRAIL_MAX (14 = ceil(CORR)) and,
  by H-trail, E[D_i] <= 0 given earlier trails. Ordering trails by start, Azuma-Hoeffding gives
  P(max_n sum_{i<=n} D_i >= 2^112 - 2^109) < 2^-2635 (exponent
  2(2^112-2^109)^2/(K_TRAILS*(100*TRAIL_MAX)^2) bits, with the drift debits below already deducted). At most K trails
  are in progress at a check; they move the sum by at most 2*(86+14)*K*TRAIL_MAX.
* *Drifts.* Two terms are not centered: the non-canonical-prefix case adds at most 86*2^-31 per step in expectation
  (below 2^103 over the run) and the per-round cap increment floor(2^20*CORR) loses below 1 per batch round (below
  2^108). Both are debited from the 2^112 trail slack.

So the cap is exceeded only on these two events or when a round's eight draws all fail (below 2^-160 over all rounds):

    DELTA_CORR = 2^-2635 + 2^-1000 + 2^-160 + DELTA_TABLE <= 5.07e-14,

subtracted from the success bound. The independence used is that of the algorithm's own coins (trail starts and rho),
not of walk states; only the average H-trail is a premise.

### Common path charged by the anchor schedule

An anchor is set only in rounds t of the global batch round counter with
`t % 64 == 0` (ANCHOR_PERIOD=64, `main` runs `batch_step(t)` for
t=1,2,..., all walks of a batch share t). Anchor rounds run the same
backward field block with the anchor tail, whose common path has 14
instructions against the plain tail's 10 in the compiled kernel (the
period counter, loaded with the anchor by one `ldp`, tested (2) and
incremented (1), and one `stp` stores the new anchor and counter), and
the same forward block. Measured on this revision's kernel at 2^16 walks and 2^20 adders, the fixed path of a clean anchor round (t=64) is
38,346,478 instructions against 38,084,331 for a plain round on both input sets: 4*(2^16-1) + 7, as in the y3-fold
kernel (39,657,108 against 39,394,961), the paired-reduction kernel (39,526,038 against 39,263,891) and the earlier
kernels (the implicit-previous-key kernel's 41,361,018 against 41,098,871). The ledger charges an anchor round 4*(2^20-1) + 9 more
than a plain round, plus the same allowance.

Of a slot's rounds t=1..R, R = L+TRAIL_MAX, at most ceil(R/64) are anchor rounds, a schedule count, not an
expectation, so every run's common path costs at most

    CONTROL = K*(256*L + 1024*TRAIL_MAX)/2^20
    COMMON = S*PLAIN + CONTROL + K*ceil(R/64)*(ANCHOR - PLAIN)

instructions, C_COMMON = COMMON/S = 580.577731 (exact rational in `accounting.py`), plus the metered corrections,
CORR per processed step. Replays are unbatched and priced separately, so their per-trail anchors do not enter this
count. The orientation enters only through the table index and the sign bit, so the fixed path costs the same for all
coin values; only the corrections depend on the data, and they are metered.

### Production counter widths and control allowance

The C file is a finite measurement harness (signed `int` loops, a 64-bit round argument, a diagnostic `ring_n`); the
specified algorithm uses unsigned 128-bit batch-round and replay-length counters and 192-bit instruction meters and
work counters. The per-trail period counter `tper` is one 64-bit word: `measurement/step.c` increments it once per 64
processed steps and compares it with CAP_PERIODS=2^60, so the cap TRAIL_MAX=64*(2^60+1) steps is enforced through that
word, which holds at most CAP_PERIODS+1=2^60+1 (61 bits); the step count is never materialised. The anchor total
ceil(R/64)~2^89.5 is analytic, not a runtime counter, and the diagnostic escape ring is not part of the attack.

Charge 256 extra scalar instructions per full main-phase batch round (all B=2^20 lanes active) and retain 1,024 for
every finishing round; finishing, rare-event and replay bookkeeping keep their allowances. The production wrapper
stores completed t, L and R_MAX as two-limb values at offsets 0,16,32 of state x19, which an AAPCS call preserves.
This listing bounds recurring main control (scratch x9..x12 reloaded each round):

    ldp x9,x10,[x19]
    ldp x11,x12,[x19,#16]
    subs xzr,x9,x11
    sbcs xzr,x10,x12
    b.hs finish
    adds x9,x9,#1
    adc x10,x10,xzr
    stp x9,x10,[x19]
    ldp x11,x12,[x19,#32]
    subs xzr,x9,x11
    sbcs xzr,x10,x12
    b.hs abort
    mov x0,x9
    bl batch_step
    b main_loop

The correction meter adds, per round (x20 = the round's block meter totals, the 192-bit meter M at [x19,#48], the cap
at [x19,#72], its per-round increment floor(B*CORR) at [x19,#96]):

    ldp x9,x10,[x19,#48]; ldr x11,[x19,#64]
    adds x9,x9,x20; adcs x10,x10,xzr; adc x11,x11,xzr
    stp x9,x10,[x19,#48]; str x11,[x19,#64]
    ldp x12,x13,[x19,#72]; ldr x14,[x19,#88]; ldr x15,[x19,#96]
    adds x12,x12,x15; adcs x13,x13,xzr; adc x14,x14,xzr
    stp x12,x13,[x19,#72]; str x14,[x19,#88]
    subs xzr,x12,x9; sbcs xzr,x13,x10; sbcs xzr,x14,x11; b.lo abort

The normal path has 15+19 instructions. The full counter is kept in state; passing its low limb preserves t modulo 64
and hence the anchor schedule, and the t bound below 2^97 excludes overflow. The wrapper does not touch the kernel or
its register allocation. With 32 register save/restore pairs, 6 stack/flags operations, 20 address/setup
instructions and 32 phase-handoff operations it is at most 124 with the meter lines, so 256 per main round covers it,
the final main comparison and transition included (L is positive). No per-lane finishing check is moved into this bound; the 1,024 finishing allowance and all event-only wide-counter prices are retained. This is an explicit
production-control specification, not a new compiled measurement.

The record counter reaches DP_MAX=2*floor(S/2^60)+2*K, below 2^68.492 records; the trail counter DP_MAX+2*K (K
starts, at most DP_MAX+K restarts); the table index spans C=2*DP_MAX, below 2^69.492 slots. Each is a 128-bit
integer whose control is charged inside the table operation. `accounting.py` charges DP_MAX table operations at
DP_COST=PROBE_CAP*PROBE_COST+DP_INDEX_COST+DP_COUNTER_COST+DP_OVERFLOW_COST
=2^10*2^11+2^11+8+8=2,099,216=2^21+2^11+16 each. A probe's two-limb increment and comparison cost at most 8 of its 2^11.
Once per table operation, the index reduction (a long division) costs at
most 2^11, the two-limb record-counter update 8 and the overflow test 8.

The home slot is (key(W)>>90) mod C, above the 30 adder and 60
distinguished bits. As p = -1 mod
2^90, a distinguished key's bits from 90 up are uniform on M=(p+1)/2^90
values under M1, so a record's slot is within total variation
rem(C-rem)/(C*M) < 1.9801e-30 of uniform, rem = M mod C. The walk reads a
distinguished key only below bit 90 and in exact equality tests (class
identity), and lookups match the full key, so coupling each recorded
class's slot with an independent uniform one changes only the probe-cap
event, with probability at most DELTA_HASH = 2*DP_MAX*(C+1)/2^166 <
7.3658e-9, which covers M=2^166.
The K initial starts and
the DP_MAX+K restarts are charged K*START and (DP_MAX+K)*START,
START=255,184,896, so the
trail counter adds at most 8*(DP_MAX+2*K) against START*(DP_MAX+2*K), a
factor of 31,898,112 for each eight-instruction counter update. No counter is
incremented per processed step: `tper` advances once per 64 steps and the
batch-round counter once per round, inside the phase-specific allowances
stated above. With the itemised charge the exact ledger
total rounds up to the declared 2^126.741585 claim.

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

    1:  stp xzr, xzr, [x0], #16      // 8 times: w0..w15
        subs x1, x1, #1
        b.ne 1b

That is exactly 10 instructions per slot (8 general-register store pairs,
a subtract, a branch), plus at most 64 per bank to select it, set x0, x1
and advance the bank counter. With N_B=ceil(C/2^40)=754,917,021 banks,
`accounting.py` charges record_table_init = 10*C + 64*N_B =
8,300,400,421,708,581,799,384 instructions, below 2^62.847887 units, as
preprocessing, so it enters both totals. Loads and stores are
instructions in the cost model; no other transfer term applies.
Preprocessing is at most 2^87.961190 units; its rounding slack, over 2^59
units, also covers setup and the 2^40 table copy (16 stores per adder). The total stays 2^126.741585.

Every other array is first written by a charged step: the adder table (both signs) by its accepted trial inside
A*2^40*(START+2^10), whose 2^14 constant covers negation and stores; a walk slot's state (at most 512 bytes) by START
at each trail start; the batch buffers (D, PRE, selections) each round before they are read, inside PLAIN and ANCHOR;
rare-path scratch inside the forced maxima; and the record, trail and batch-round counters with the rare-work and
correction meters in setup, inside the 2^24 input_validation allowance.

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
       + 2^-63*256*(w2+wL)          false filter, A-independent
       + (tau_exact*TRAIL_MAX/S)*wL merged step, A-independent
    tau_exact = T^2/(2N).

The longer-entry term is charged the derived bound above; the merged term retains its predecessor price. The
false-filter term carries a generous debit for speculative look-ahead: a false anchor check exposes at most 65 future
classes (the debit was sized for two filters), which the concentration proof cannot then treat as unexposed, so the
originating episode is charged its direct work, at most 66 later exposed-source steps and at most 128 persistent
steps through an ensuing short cycle; the debit 256*(w2+wL) exceeds that sum, and repeated charges of overlapping
exposures only overcount. An unclosed speculative path that returns to visited classes merges with another trail
(charged), closes a short cycle (this debit and the 128-step bound) or violates G (failure); otherwise continuing past
its end exposes a new source, so no known look-ahead key is treated as fresh. The runtime meters the same calls.

For the residual success argument, stop the analytical residual process at the first G violation, as for the
two-cycle process. With Eres=E0 and kM=400000,

    B_res = kM*Eres*S,
    B_rare = B2+B_res.

The factor kM is a free analysis parameter: raising it multiplies the reserve kM*Eres per step and shrinks the
abort allowance 1/kM. At A=2^30, E0 is at most 5.972949335e-11 (longer 1.315215518e-12, false 5.810604775e-11,
merged 3.082300770e-13). Tau and kM are tuned jointly; this revision keeps kM=400000 and tau=4943092/10000000 (below),
for a reserve of 0.0009654 instructions per step and an abort allowance of 1/400000; the threshold is public
algorithm data.

All within-trail coincidences and excessive-length formal closures are excluded by G and paid through DELTA_CUT, not
by the per-step allowance or the predecessor's within-trail fixed-point argument. On G the stopped residual process is
the actual residual cost, and by the retained rate bounds plus the speculative-exposure debit its expectation is at
most Eres*S, so Markov gives a 1/kM = 1/400000 tail.

    Pr[rare work > B_rare and G] <= 2^-100 + 1/400000.

No independence between residual and two-cycle work is asserted. Rare operations are metered at their worst-case
prices, as in the predecessor, with an abort before any priced operation that would exceed B_rare, so work is bounded
even when the premises fail; the decomposition adds no runtime counter, freshness oracle or uncharged lookup, and its
parameters are public algorithm constants.

## Complete cost and success ledger

The ledger charges (COMMON+S*CORR+B_rare)/1000 for all walk steps, metered corrections and metered rare work, then
**adds**, rather than amortizes away:

* the correction meter's slack LAMBDA_AZ+LAMBDA_NF and one round's overshoot K*86;
* 10 instructions for each trail-ending batch exit, K_TRAILS*10 (anchor exits are inside wL);
* A adder constructions at 2^40 small-key trials each, and K initial trail starts, each trial or start with the
  complete capped scalar-sampling and two-scalar-multiplication price (trials also 2^10 for the key test);
* up to DP_MAX+K restart constructions, DP_MAX=2*floor(S/2^60)+2*K;
* DP_MAX table operations at 2^21+2^11+16 instructions each (probing, key comparisons, addressing, index
  reduction, counter update, overflow test), and the record table's zeroing, 10*C + 64*ceil(C/2^40) instructions;
* two bounded coefficient replays at the generous unbatched operation count; solving, verification, serialization,
  input validation and setup allowances;
* the Lemma F closure series, summed for s < 4096 with the tail bounded by the consecutive-term ratio
  (e^2/8)(1+3/A)^2 < 0.924, so the truncation is an upper bound, equalling the enumerated sum at A=2^16 and A=2^20.

All processors' work is summed; no future-target amortization or free advice is used. The storage formula is
unchanged, 2*DP_MAX*128 + K*512 + 2*A*2*64 + 2^24 bytes, covering records, walk states and batch buffers, both adder
signs (and the block's copy) and a fixed reserve; at A=2^30 it evaluates to 2^76.491743 bytes, dominated by the
2*DP_MAX*128 record table. `accounting.py` evaluates these formulas and asserts that the total, rounded up at the
sixth decimal, equals the declared score 25348317/200000.

The collision exponent penalty, re-derived at the bound with the predecessor's 2.6 margins, is
epsilon = 196/A + 288/A^2 + 768/A^3 = 1.825392e-7 (section "The adder table"). The conditional success lower
bound is

    1-exp(-tau_exact*(1-epsilon))
      - 1/kM - 2^-100 - 2^-20 - DELTA_CORR - DELTA_ADDER_CLASS
      - DELTA_ABANDON - DELTA_CUT - DELTA_MERGE
      - DELTA_CAP - DELTA_PROBE - DELTA_HASH.

Each term is recomputed from the predecessor's formula at the new S and T; DELTA_CORR is the correction meter's
abort term. kM=400000 was cheapest in an earlier 1e5-to-2e6 screen (1.0e-6 margin then), and since the
checkpoint-only revision tau is the smallest on the 1e-7 grid keeping the asserted 5e-7 floor; this revision keeps both. At tau=4943092/10000000 the model term is
0.390007800602308 and the exact bound is 0.39000401600773443, a margin of 4.016e-6 above 0.39; the named failure terms at this S are DELTA_ABANDON=7.2117e-8, DELTA_CUT=2.5144e-7 (accidental
relations bounded with M_ACC in place of r, section "Small-key adders"), DELTA_MERGE=1.2521e-14,
DELTA_CAP=9.3782e-18 (adder exhaustion included), DELTA_ADDER_CLASS=4.4370e-30, DELTA_PROBE=2.4924e-44,
DELTA_HASH=7.3658e-9 and DELTA_CORR=4.8106e-49. The bundled ledger asserts a margin of at least 5e-7 and the claimed
cost. The margin is a conditional-model margin, not a confidence interval for empirical trials: five times the
predecessor's 1.0e-7 floor, which only had to cover numerical error in the recomputed bound (below 1e-14).

Per processed step the ledger charges C_COMMON=580.577731 (fixed path,
anchor schedule, control with the blinding draws), CORR=13.959715
(metered corrections) and the rare reserve 0.0009654 (two-cycle
0.0009415, residual 0.0000239; rounded up), in total 594.538411
instructions (614.650398 in the y3-fold, 616.400323 in the paired-reduction, 641.400295 in the deferred-inv, 642.149317 in the checkpoint-only and 644.068668 in the small-key revision, rounded up), against
714.759990 in the predecessor (714.754495 + 0.005494). The correction meter's slack
and overshoot add below 2^102.052 units. The total is below 2^126.741585.

## Evidence scope, execution and limitations

Inherited forced rare paths, replica results and reduced-order records keep their attribution and scopes
(the sign-law and Lemma R records and 5a258d4a's `smallkey.json`, not shipped here, come from kernels or procedures
with the two-steps-back filter; the offline records `trail-paths.json`, `replay16-production.json`,
`trails/summary.json` and `trails-stress.json` are omitted for the size limit, byte-identical in 5a258d4a's,
ae463634's, e2973980's, 21edf4bd's and 205a498f's folders, as is 205a498f's `measurement/gen13.py`). New here:
`measurement/gen49.py` and its output, the `step.c` edits of "Shared rows" and `evidence/measurement/branchy.json` key
`shared_rows` (builds, static counts, emulator and block-harness exactness, mutants, paired, 2^20-walk and anchor
rounds, meters, site rates, dumps), with `ancestor_quotes` (earlier kernels' quoted values; their records are in
5a258d4a's file and ae463634's, e2973980's and 205a498f's, the last with `y_fold` and `paired_rows`). None
establishes a production success probability;
pooled restart points in the harness do not measure cold restarts, which
the attack pays in full.

The fixed path is the same for every input, and correction work is metered
and capped, so H-trail moves only DELTA_CORR, not the charged work.

One declared combined experiment uses the verifier's bounded credential-free runner and returns two separately
labelled results: `accounting.py` recomputes the ledger, the correction charge, its slack and abort terms, the
concentration threshold, failure terms and the point-order certificate; `experiment.py`, the inherited complete
reduced-order procedure with checkpoint-only cycle detection at the declared collision parameter on its 0.0001 grid,
keys its table by the full key (no home slot, probe cap or zeroing), draws its adders by the key rule alone (no y
window), implements no production meter and runs no full-size attack; its finite-run successes, failures and replay
checks are diagnostic, reported as observed. The assembly blocks cannot run in the Python runner; their checks are
recorded above. Parent files no executed path reads are omitted, the evidence list names exactly the files present,
and `combined.py` reads only `accounting.py` and `experiment.py`. The forced-path maxima link to the predecessor
kernel's recorded sources (`measurement/forced/summary.json`, stored without whitespace, content unchanged); the rare
primitives they price compile to the same instructions here.

The remaining material uncertainties are production-scale collision mixing, conditional fresh-source orientation,
the fresh-window domination model, analytic handling of long trails and rare returns, validity of the measured kernel
bounds, H1 for the small-key y-window adder table, and the average correction premise H-trail with the fold model
for dy (abort term only). The
exponential-tail arguments do not resolve them but are conditional on them. Memory remains exponential and
impractical, and no claim of practical compromise, unconditional improvement in P-256 security, or human
certification is made.
