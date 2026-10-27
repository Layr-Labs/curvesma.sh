This revision changes two coupled details of submission
70c0711e-ed9a-4be7-a5f6-03165ce1e6ae (declared 126.868201, folder
`blinded-implicit-prev`, commit 306bc0f8; its kernel is the "implicit
previous key" kernel below): **every adder is rejection-sampled until its
Montgomery x-coordinate is below 2^224**, and the batch forms
**D = x(W)-x(R) and t = x3-x(W)** instead of x(R)-x(W) and x(W)-x3 (item 8;
section "Small-key adders and reversed differences"). Walk states, keys,
signs, distinguished points, every rare path and the fixed path are
unchanged. Two metered corrections that together fired on average 5
instructions per step *for every table* (D borrowed with probability
1-mean x(R)/p, the second x3 subtraction with mean x(R)/p) now each fire
with probability below 2^224/p < 2^-31.99999, so CORR falls from 22.468426
to 17.458661 and the declared score from 126.868201 to 126.857022. The
table's worst-case construction, 2^81.961189 units, is charged in full.

Lineage (content kept below): 70c0711e (implicit previous key), 961165d1
(eight copies), 1f942c41 (flag-carry finishes), c2b498ec (dy pad 0),
b75c69bb, 778628fb and earlier kernels, 7169b17e..8a577856 (trail-level
correction cap), material predecessor b89d1d72.

# P-256 sigma rho with a 2^24-adder table and a concentrated two-cycle reserve

## Current claim and provenance

For a uniform scalar x in 0..r-1, given Q=[x]P on the protected NIST P-256
profile, the following classical single-target algorithm has success probability
at least 0.39 under the disclosed H1, M1/L-ori, fresh-class exposure and H-trail premises.
Its total charged computation is at most **2^126.857022** target-group-operation
equivalents. The normalization is C=1000 scalar AArch64 instructions per unit, not one
unit per P-256 addition. Peak memory is below 2^77 bytes, preprocessing
(the rejection-sampled adder table and record-table zeroing included) is
at most 2^81.961192 units and is included in the total, and nonuniform
advice is zero bytes; instruction accounting retains H3.

The exact declared score is 126857022/1000000 (126.857022), the ledger's total
rounded up at the sixth decimal. This is a conditional analytic bound, not an
unconditional complexity theorem or an executed full-size attack.

The immediate source of this revision is submission
`70c0711e-ed9a-4be7-a5f6-03165ce1e6ae` (commit
`306bc0f8b8492eaacaee3a9f3059d94bb77d78f9`, directory
`submissions/p256-ecdlp/blinded-implicit-prev`), with the lineage above.
The material predecessor is submission `b89d1d72-d9a1-4ef9-9ed7-c6ef0eec1ed3`, immutable
commit `0e1d63e9afd0fdb04b311b1ae72a2618271b371b`, directory
`submissions/p256-ecdlp/rho-sigma2`, whose algorithm, arguments and
measurements `predecessor-proof.md` reproduces with a historical banner. No
prior platform decision is an accepted premise.

**Precedence:** the claim is this document and `submission.json`. The
predecessor's Sections 1.2--1.3 define the same algorithm except for A, L
and the rare-work threshold replaced below; Sections 2.0--2.4 and 2.8 give
the unchanged probability-space, escape, replay, abandonment and merge
arguments, re-evaluated here; Sections 3--4 the forced-path pricing,
resource caps and storage; the batched step and its measurement are
replaced, and the record-table layout and zeroing added, below. Its score
and parameters are historical; its `evidence/accounting.json` remains in
the predecessor submission (omitted for the size limit).

## Algorithm, fixed costs and changed parameters

Use the protected p, a=-3, b, r and generator P without alteration. Handle
Q=identity by returning zero; otherwise validate the curve point. Let a
class be {W,-W}, key(W) the integer Montgomery representation of its affine
x-coordinate, and canon(W) the representative with even Montgomery y.
Let sigma(W) be the sign taking canon(W) to W. The walk is

    W -> W + sigma(W) R[key(W) mod A].

Here A=2^24 adders, each R_i=[c_i]P+[d_i]Q with fresh uniform scalar coins,
redrawn until key(R_i) < 2^224 (at most 2^40 trials; section "Small-key adders").
The table contains both signs. There are K=2^32 walk slots, arranged into
batches of 2^20 for simultaneous inversion. Distinguished points have
bits 24..83 of the key zero, so theta=2^-60 under the stated class model.
Each trail starts from fresh independent a,b with W=[a]P+[b]Q. It ends at
its first distinguished point, records it once with its coins, then starts
a fresh trail. Coefficients are reconstructed by replay only at a match.

The short-cycle rules, test ordering and replay algorithm are those
specified in predecessor-proof.md Sections 1.2--1.3. In particular,
the distinguished-point test precedes cycle tests at every generated
point, including look-ahead and escape points. A detected two-cycle or
cycle of length 3..64 escapes by doubling the canonical representative
of its least-key member. The two 64-bit filters only trigger exact checks.
Every trail is limited to TRAIL_MAX=64*2^60+64 processed steps; the first
key match is resolved by two replays, each limited to 65 times its original
live length. Degenerate matches and exceeded caps fail the run. Each
derived scalar is verified by [x]P=Q before output.

Set N=(r-1)/2+1 and

    tau_target = 4943041/10000000,
    T_min = ceil(sqrt(2*N*tau_target)),
    L = ceil(T_min/K), T=K*L,
    S = K*(L+TRAIL_MAX).

Every slot takes L main-phase steps and finishes its current trail. Thus S
bounds all processed steps, including persistent steps and finishing tails.
Sampling caps, maximum distinguished records, abandonment cap, probe cap,
restart costs, bounded replay, verification and output are unchanged
formulas, evaluated at this S. The point-order certificate proving
ord_r(2)>2*TRAIL_MAX is recomputed in `accounting.py`.

The common path costs at most C_COMMON=626.604513 instructions per
processed step under H3 (anchor schedule below), and corrections at most
CORR=17.458661 more, metered (below; dy carries pad 0 by sign symmetry, D and the second x3
subtraction pad 0 by the small-key bound, the other sites 2^-10). A fresh-source two-cycle
entry costs at most w2=83,518 more, a longer entry wL=2,789,687; false
filter episodes with their exposure debt are charged 256*(w2+wL), merged
steps wL, the inherited forced-path/control allowances.
The runtime still uses **one total rare-work meter**: it need not
recognize whether a point is globally fresh or classify observed events
as fresh versus merged. That partition is solely for analysis.

### The adder table at A=2^24 and the measured step price

Every entry rate and reserve term is a function of A. At A=2^24:

    fresh-source two-cycle rate  p2 = 1/(2A) = 2^-25,
    longer-entry rate            <= 3/(4A^2) + 2/A^3 = 2.664536e-15,
    conditional entry bound      P' = 64/(2A) = 2^-19,
    collision exponent penalty   epsilon = 10/A + 288/A^2 + 768/A^3 = 5.960475e-7.

The two-cycle rate is exact under L-ori. The longer-entry rate is at most
3/(4A^2)+2/A^3 (Lemma E1(c)): E[w'] <= 3/(2A)+4/A^2 over windows of length
3..64, so E[w']/(2A) <= 3/(4A^2)+2/A^3, which the ledger charges and
`accounting.py` asserts equals (3/(2A)+4/A^2)/(2A).
Epsilon's closure margin uses the same bound, each 3..64 window wasting 192
states: its third term is 2*192*(3/(4A^2)+2/A^3). The split-concentration threshold mu=S/(2A), the
two-cycle reserve ceil((1+eta)*mu)*w2, the residual expectation E0 and the
reserves built on them are recomputed, not rescaled: the two-cycle charge is
0.002489 and E0 is 7.51339e-9 instructions per step, against 0.637201 and
4.87164e-4 at A=2^16.

The A-dependence of the compiled step is checked on the binaries (flags as
measured, B=2^20, ADDER_BITS=20 and 24, disassembled): `batch_step` has 11,568 instructions at both A, followed by seven
(A=2^20) or five (A=2^24) alignment nops that never execute; the two
differ only in the `extr`, `ubfiz` and `and` immediates of the DP test and
the adder index (34 instructions: sixteen `extr` and sixteen `ubfiz` in
the sixteen block step tails, one `and` and one `extr` in the compiled C
step tail);
`step_rare`, `new_trail`, `add_single`,
`finv` and `dp_restart` are identical up to immediates, and `rare_anchor`
and `rare_2cycle` have one more instruction at A=2^24, inside the per-event OVERHEAD. The two 64-bit filters compare the whole low limb, so their
2^-63 false-match charge is independent of A, and the table is data.

The larger table is charged, not assumed free: construction (worst case
of the small-key rejection sampling) A*2^40*(START+2^10) = 2^81.961189
units, storage A*2*64 = 2^31
bytes, twice with the block's copy at 2^40, inside the peak-memory formula, and each 128-byte record's read-back
(16 words) lies inside START's 2^14 constant term. Adders and initial
starts are 2^81.961189 units of preprocessing; with the record-table zeroing
below, preprocessing is at most 2^81.961192 units, included in W.

### The batched step in this revision

The walk, its tests and every rare path are the predecessor's. Only the
implementation of the common batched step changes, in seven ways, each
value-preserving:

1. **Branchy corrections.** In the predecessor kernel every conditional
   correction of the field arithmetic was computed branch-free and always
   executed. They are now forward branches over their fix-ups, and a
   taken fix-up adds its own extra instruction count to a meter (below).
2. **A register-resident forward block.** The prefix products PRE[w+1] =
   D[w]*PRE[w], w=1..B-2, run as one assembly loop (`FLOOP`) with PRE kept
   in x19..x22: the product's last reduction step writes the new PRE over
   the multiplier registers it has consumed, so PRE is never reloaded.
3. **Operand forms.** inv is kept canonical (branchy finish), so as
   register factor it lets iw = inv*PRE use the five-word rows instead of
   six-word rows. PRE, iw and lambda are left below 2^256 (their consumers
   accept that); the products' last reduction steps write straight into
   their destination registers.
4. **2^20-walk batches** (per-round inversion and set-up amortized over
   16 times more walks) and a copy of the two-sign table T2 at virtual
   address 2^40, read by the blocks, so the adder address is one `orr`
   and x27 is free for the meter; the C code keeps the static T2.
5. **A blinded prefix chain.** Each global round draws one fresh uniform
   unit rho, shared by all batches (rejection sampling, at most 8 draws),
   and every batch starts its prefix chain at
   PRE[1] = rho*D[0]; the C path for w=0 takes iw = inv*rho = 1/D[0]. All
   inverses, slopes and points are unchanged; every prefix PRE[w] is then a
   uniform unit given the walk, which the meter's abort bound uses.

6. **Eight copies per loop.** Both blocks hold eight copies of the step
   (distinct labels); the end test (`tst x23`/`b.ne`, or `subs x26,#8`/
   `b.ne`) runs once per eight steps, copy j processing lanes w = -j
   mod 8. The backward block enters through a `tbz`/`tbnz` tree on bits
   6-8 of x23, the pointer to 64-byte DP record w (array aligned to its
   size): 3 instructions at w = B-1, at most 4 on a re-entry after a rare
   exit (two-copy entry: 1); its end test follows the w = 1 mod 8 copy.
   The forward block starts with one `b` to copy 2 and counter B, running
   6 + 8(B/8-1) = B-2 steps. Every walk is processed as before (B a
   multiple of 8, asserted at compile time).
7. **Implicit previous key.** The filter's previous key is read from
   word 0 of Wn[w].x, not a lane-record copy (section below).
8. **Small-key adders, reversed D and t** (this revision; section below).

`measurement/step.c` is the predecessor's harness with these changes; the
backward block `BLOOP_FIELD_B` (`measurement/bloop_branchy.h`) and the
forward block `FLOOP` (`measurement/floop.h`) are generated by
`gen10.py` from the predecessor's measured `BLOOP_FIELD` text (the
generator stays in submission 70c0711e, `measurement/gen10.py`, omitted
here for the size limit; this revision's blocks are its output after
`measurement/gen11.py`): every multiplication, square and reduction row
is that text, and only the corrections and the destinations of the last
reduction steps differ.
Item 7 changes the step tails, the rare exit, `step_rare` and
`new_trail` (section "Implicit previous key"); `rare_anchor`,
`rare_2cycle`, `add_single`, `finv` and `dp_restart` disassemble identically to the
eight-copy build up to address immediates, so the forced-path maxima of
`measurement/forced/summary.json` still price the rare primitives.

### The corrections

With x28 = p1 = 2^32-1 and x29 = p3 = 2^64-2^32+1 set at block entry
(p = (2^64-1) + p1*2^64 + p3*2^192), adding p to four limbs is
`subs r0,r0,#1; adcs r1,r1,x28; adcs r2,r2,xzr; adc r3,r3,x29` and
subtracting it is `adds r0,r0,#1; sbcs r1,r1,x28; sbcs r2,r2,xzr;
sbc r3,r3,x29` (carry chains as in the predecessor's finishes).

* **Subtraction** a-b of canonical a, b (six per walk step: D, dy, two in
  x3, t, y3): `subs/sbcs` x4, then `b.cs` over [add p; meter += 5]. The
  common path is 5 instructions (was 11); a borrow costs 5 more. The
  result is a-b or a-b+p, canonical.
* **Carry-only finish** (PRE, iw, lambda): T = r + c*2^256, c in {0,1},
  T < 2^256+p, where c is the carry flag left by the product's top `adcs`.
  `b.cc` over [subtract p; meter += 5]: 1 instruction (was 7 or 9), 5 more
  when c=1; the result is T - c*p < 2^256.
* **Canonical finish** (inv, y3' = t*lambda; T < 2p): `b.cc L` (c=0), else
  subtract p, meter += 4, `b done`; `L: cmn r3,x28; b.hs slow`. Common
  path 3 (was 9), c=1 costs 4 more. If c=0, T < p unless r3 >= p3
  (`cmn` sets C exactly then); the out-of-line slow block computes T-p,
  keeps it with four `csel` on no borrow, meter += 10, and returns.
* **lambda^2**, from lambda < 2^256: T < 2^256+p, so after a c=1
  subtraction T-p may still be at least p. The fix-up (meter += 5) falls
  through into the same `cmn/b.hs` re-test; one more conditional
  subtraction in the slow block always suffices (T - 2p < p).

### Flag-carry finishes (from 1f942c41)

The six product finishes (iw, inv, lambda, lambda^2, y3', PRE) branch on
the carry flag C. In c2b498ec each product's last reduction row ended
with `adc c,xzr,xzr` (c = x13, x11 at lambda^2) and the finish began with
`cbz c,L`; that and `b.cc L` both branch exactly when C = 0, so dropping
the `adc` and replacing `cbz c` keeps every branch outcome, value and
meter increment if (a) the branch follows the `adcs` whose carry is c
and (b) nothing later reads c. `gen10.py`'s `FCARRY=True` selects this
form, `FCARRY=False` c2b498ec's (six sites per copy apart). On each of
this revision's four builds (B=2^20 with A=2^20 and 2^24, B=2^16 with
A=2^16 and 2^20) all 88 `b.cc` sites (eight forward, eight plain and
eight anchor backward copies) follow the top `adcs`, and on every path
from each, 19 to 117 instructions deep, x13 (x11) is written before it
is read; both are declared clobbers, unread by C code. Against
c2b498ec's build, 1f942c41's two-copy `batch_step` lost exactly 22 `adc
..,xzr,xzr` and 22 `cbz` became `b.cc`, all else identical.

Five-word rows need a register factor below p (D, dy, t, inv) and accept
any 4-limb multiplier; iw and lambda < 2^256 enter only as multipliers and
the square. Every exposed value (x3, y3, key, sign bit) is
canonical, so walk states, dumps and rare-path inputs equal the
predecessor's.

**Exactness checks** (`evidence/measurement/branchy.json`): (i) n=8 dumps
of 2^16 walks on two input sets and stress builds with 2^4, 2^6 and 2^10
adders (every rare path) equal the predecessor kernel's byte for byte;
(ii) the backward block, wrapped as a function, equals an exact
big-integer model (Montgomery integers, branch outcomes, meter) on
240,000 random and crafted inputs covering every path (22,118
non-canonical lambdas, 10,063 to 41,984 cases of each slow path and its
subtracting branch, operands 0, 1, p-1, p, 2^256-1); (iii) the forward
block equals the model on 6,000 runs of 62 steps (PRE inputs up to
2^256-1, 2,403 non-canonical PRE outputs).

For 1f942c41 the three checks were rerun on the flag-carry blocks
(byte-identical dumps and stress builds, 0 model mismatches).

### Eight copies per loop (eight-copy revision)

Each compiled copy equaled 1f942c41's copy instruction for instruction
(forward 119, plain backward 570, anchor backward 575, slow sets 30);
`batch_step` differed only in the entry tree, the forward counter's
immediate, the entry `b` and nops, so checks (ii) and (iii) covered every
copy, and check (i) matched 1f942c41's build byte for byte. Loop control
per round fell by exactly 3B/2-5 (1,572,859 at B=2^20); a re-entry after
a rare exit costs at most 3 more, inside its OVERHEAD of 1,024.

### Implicit previous key (previous revision)

Earlier kernels kept a 40-byte lane record {prev, anchor, cur, sel,
tper}, and each step loaded cur, stored it as prev and stored the new
key as cur. Round t reads P_t from buffer Wc and writes P_{t+1} into Wn,
which are swapped after the round, so Wn[w] holds P_{t-1} as round t
begins; word 0 of Wn[w].x, the key the filter needs, is unchanged until
step w's tail loads it, on every path:

* Before the tail, step w stores only x3[2..3] and y3 into Wn[w]; the
  forward pass writes only the separate D/PRE array; every other step,
  rare exit and C call writes only its own walk's entries.
* First round: `main` writes ~0 = 2^64-1 into word 0 of every Wn entry,
  the old prev's initial value.
* Trail start (`dp_restart` after a distinguished point, abandonment at
  the period cap): `new_trail` writes ~0 into word 0 of Wc[w], next
  round's Wn[w], where the old code set prev = ~0. Wc[w] is not read
  again that round (`rare_2cycle` reads it first; the pass is past w).
* Otherwise word 0 of Wc[w] keeps P_t's key, as old prev := cur did; an
  escape replaces only P_{t+1}, in Wn[w].
* Anchor rounds use the same load (anchor and counter stay in the record).
* A rare exit (distinguished point, filter hit, period cap) leaves after
  the tail loaded the key into x5 and stores it as the record's prev,
  which `step_rare` reads; re-entry continues at walk w-1. The slow
  correction blocks rejoin their copy and never reach the exit.
* The C step (walk 0) sets prev from word 0 of Wn[w] before `k_lxy`
  writes the new point there.

The counted builds compare the whole low limb (FILTER_MASK = 2^64-1), so
every test sees its old operands; the 8-bit-filter stress build masks
the loaded key in three more tail instructions.

Static changes, all four builds: each block copy drops `stp
x3,x4,[x25,#-48]`; the plain tail keeps 12 instructions (two loads, the
moved pair, one `str` of sel for `stur`+`stp`) and the anchor tail falls
from 17 to 16, so a plain step loses 1 and an anchor step 2 (copies 569
and 573); forward copies read sel at stride 32, not 40. `batch_step` has
34 fewer instructions at B=2^20, 35 at B=2^16, nops aside. Checks (ii)
and (iii) cover the unchanged arithmetic; check (i), rerun, matches the
eight-copy build byte for byte (n=8 dumps, four stress builds, both
input sets). All 88 `b.cc` sites pass the flag check on all four builds,
and every memory operand (558 per compiled `batch_step`) has an
encodable AArch64 offset.

Rare paths: an exit adds the `str` of prev, and re-entry drops 3
instructions of a division by 40 (records are now 32 bytes);
`step_rare` compiles to 105 instructions instead of 121 (`rare_2cycle`,
still a separate call, is unchanged), `new_trail` to 49 instead of 47 (the ~0 store). Calls to
the priced primitives are unchanged; each body is far inside the
per-event OVERHEAD of 1,024, and a trail start's slot writes are inside
START.

### Small-key adders and reversed differences (this revision)

**Algorithm change.** Step 2 of the procedure (predecessor Section 1.3)
becomes: for each adder i < A, repeat at most ADDER_TRIALS = 2^40 times:
draw fresh c, d under the caps of Section 8.2, set R = [c]P+[d]Q, and
accept when R != O and key(R) < 2^224, i.e. the top 64-bit limb of the
Montgomery x-coordinate is below 2^32 (a load, a shift and a branch); then
tabulate R_i = R and -R_i with coins (c_i, d_i). If any adder exhausts its
trials the run fails. Nothing else in the walk, its tests, its
distinguished points or its resolution changes; the walk still reads
sigma(W) R[key(W) mod A] from the two-sign table.

**Kernel change.** The batch now forms D = x(W)-x(R) (forward block FLOOP,
and the C code for walks B-1 and 0) and t = x3-x(W) (backward block
BLOOP_FIELD_S and the C step `k_lxy` of walk 0). `measurement/gen11.py`
derives this revision's `floop.h`, `bloop_branchy.h`, `p256_fused.h` and
`step.c` from the previous revision's files by swapping the two operands
of exactly those five four-limb subtraction chains (each edit asserted to
match once; same instruction count and registers). The inverse chain uses
D' = -D, so iw' = -1/D, lambda' = dy*iw' = -lambda, lambda'^2 = lambda^2
and x3 is unchanged; t' = -t gives y3' = t'*lambda' = t*lambda, so y3 is
unchanged. Every exposed value (x3, y3, key, sign bit) is the previous
kernel's. PRE[w] becomes (-1)^(w-1) PRE[w], still a uniform unit given the
walk under the per-round blinding, so the product-carry sites keep their
law.

**Why two corrections disappear.** For a canonical u in [0, p),
u - x(R) borrows iff u < x(R). D' = key(W) - x(R) borrows iff key(W) <
x(R), and the second x3 subtraction (lambda^2 - x(W)) - x(R) iff its
canonical first difference is below x(R). With every x(R) < 2^224, each
has probability below 2^224/p < 2^-31.99999 for a uniform operand, for
*every* adder, so these two sites need no table-mean pad and no Hoeffding
term over the table. Before, D = x(R)-x(W) borrowed with probability
1 - x(R)/p and the second x3 subtraction with x(R)/p: 5 instructions per
step on average whatever the table. CORR changes by
-2*5*(1/2+2^-10) + 2*5*2^224/p = -5.009766 to 17.458661.

**Measured on this kernel** (`evidence/measurement/smallkey.json`; GCC
14.2.0 Ubuntu cross build 14.2.0-4ubuntu2, the organizer flags, `-static`,
QEMU 8.2.2 user-mode exact counting as described below, harness functions
only; inputs are random curve points from `random.Random(seed)`, adders
with Montgomery x below 2^224 for the small-key sets):

* Walk dumps (8 rounds of 2^16 walks, 2^16 adders) are byte-identical to
  the previous kernel's on a uniform-adder set and a small-key set; event
  counters equal, meter aside.
* Fixed path (round minus meter increase), B=2^16, plain round 2: 41,098,898
  for both kernels on both sets. This toolchain counts the previous
  kernel 27 above the organizer image's 41,098,871, the offset the method
  documents; the new kernel's fixed path equals the previous kernel's
  exactly. B=2^20, A=2^20, small-key set: clean round 1 656,973,460 for
  both kernels (organizer image 656,973,433 for the previous kernel: the
  same offset 27); round 2, one two-cycle escape, 657,055,820 for both.
* Static, B=2^20, A=2^24: `batch_step` has the same instruction count in
  both builds; the address-normalised diff is exactly 100 swapped operand
  pairs in 25 subtraction chains (8 forward copies, 16 backward copies,
  `k_lxy`) and 12 swapped base registers of loads in the three C
  subtraction helpers. `step_rare`, `new_trail`, `add_single`, `finv`,
  `rare_anchor`, `dp_restart` and `main` are identical, so the forced-path
  maxima still apply.
* Per-site counters (one build per site, 8 rounds of 2^16 walks, small-key
  set): D and the second x3 subtraction fire 0 times in 524,272 and
  524,280 opportunities (previous kernel: D every time); other sites keep
  their frequencies (subtractions 0.498-0.500, products 0.249-0.251, square
  0.333). Meter 17.40 per step (17.41-17.42 at 2^20 walks) against 22.41
  for the previous kernel on the same sets, below the charged 17.458661.

**Table cost and failure.** Each trial costs at most START + 2^10
instructions (START as before: capped draws, two scalar multiplications,
retries on O, stores; 2^10 for the Montgomery-x test and loop control).
The worst case A*2^40*(START+2^10) = 4,707,349,357,454,119,741,916,446,720
instructions = 2^81.961189 units is charged as preprocessing. Acceptance:
a trial point is uniform on the group (c is uniform and independent of d
and x). The affine points with Montgomery x = u < X number X + S_X up to
three roots of f, where S_X = sum_{u<X} chi(g(u)), chi the quadratic
character and g(u) = f(u 2^-256) = (u 2^-256)^3 - 3(u 2^-256) + b, a
squarefree cubic (the curve is nonsingular). The Weil bound with the
standard completion over an interval gives |S_X| <= 3 sqrt(p) (1+ln p) <
2^137.07; the ledger uses 2^160. So at least 2^224 - 2^160 - 3 group
elements are accepted, q_acc >= 2^-32.0000000004, and all A adders finish
except with probability A (1-q_acc)^(2^40) < 2^-345.3, added to DELTA_CAP
(2.5477e-18 with the 2*A*2^40 extra scalar draws). Expected trials per
adder are about 2^32.

**Probability space and accidental relations.** Every walk point is
[l]P for a logarithm l that is a linear form, with integer coefficients
fixed by the transcript, in the trail's start logarithm t = a + b x and the
adder logarithms s_i = c_i + d_i x; the walk, its keys and its tests are
functions of these logarithms only, and x enters only at resolution. For
fixed x, (c, d) -> (c + d x, d) is a bijection of Z_r^2, so: each start
logarithm is uniform on Z_r; each accepted s_i is uniform on the accepted
set S of logarithms (key([s]P) < 2^224), of size m >= M_ACC = 2^224 -
2^160 - 3 (count above), closed under negation; all of these are
independent of each other and of x; and given them, every d_i and every b
is still uniform (so the degenerate-match bound, which uses the b and d
coefficients, holds as stated). **Acceptance does constrain the s_i,** so
an accidental relation among adder logarithms is likelier than over Z_r:
a nonzero linear form in the s_i (and starts) vanishes with probability at
most 1/m. Fix the form and an index j with nonzero coefficient; given all
other variables it pins s_j to one value, which a uniform element of S
takes with probability at most 1/m (1/r if the coefficient is on a start).
For example two adders satisfy R_i = +-R_j with probability 2/m, not 2/r.

Every accidental-coincidence bound of predecessor-proof.md Section 2.2
(the in-trail pair bound 2|Γ|/r, its conditioning on earlier
non-coincidence, Lemma G's |Γ|t/(r - |Γ|t(t-1)/2) and the table term
A^2|Γ|/r) is a union of such events, so each holds with r replaced by
M_ACC; this is applied to all of them, including the terms that also
contain a fresh start logarithm, which only over-counts. Their side
condition becomes |Γ|TRAIL_MAX^2 <= M_ACC/2 (asserted). In the ledger
DELTA_LOOP, hence DELTA_CUT, rises by the factor r/M_ACC < 2^32 from
5.8542e-17 to 2.5144e-7; it is subtracted, and the success bound becomes
0.39000078761, a margin of 7.876e-7 above 0.39 (the asserted floor is now
5e-7; tau, S and every work term are unchanged). Formal closures (D1, D2,
Lemma F) are combinatorial facts about transcripts with i.i.d. class data
and do not involve the adder logarithms' distribution. Merges between
different trails are the intended collisions, governed by H1.

M1 and L-ori concern walk classes. The selection conditions only the 2^24
adder classes; the walk ever processes a point in one of them with
probability at most S*A/N < 6.94e-32 (`DELTA_ADDER_CLASS`, subtracted).
H1 (random mapping for the class walk with the A-adder penalty epsilon)
depends on the table only through its A elements, whose logarithms are
independent and uniform on S; S is defined by a coordinate bound with no
known relation to discrete logarithms. Extending H1 to this table is a
disclosed premise at the level of H1 itself, not a theorem. Toy support
(`evidence/measurement/smallkey.json`, `toy`): the complete procedure on
28-bit a=-3 curves at tau=0.4943, with adders redrawn until key < p>>4,
gives 513/1200 successes (three seeds) against 436/1000 with uniform
adders; two-cycle rates 2.37e-4 to 2.46e-4 against 2.41e-4 (model
1/(2*2^11) = 2.44e-4) and no replay mismatch. The declared experiment runs
the rule with `adder_key_shift` 2 (keys below p>>2), which fits the
runner's time limit (shift 4 did not).

### Measured fixed path and metered corrections

c2b498ec's counts were taken with GCC 14.2.0 (Ubuntu cross build
14.2.0-4ubuntu2), the organizer flags `-O3 -fno-tree-vectorize
-mgeneral-regs-only`, by exact guest-instruction counting under QEMU 8.2.2
user mode (block sizes from `-d in_asm`, executions from `-d
exec,nochain`, restricted to the harness's functions main..read_fe and
add_single..the end of batch_step; identical to `-one-insn-per-tb`
counting; a 1,000-iteration three-instruction loop adds exactly 3,000).
Valgrind 3.22.0 Callgrind reproduced its round differences on the
predecessor kernel (46,832,945), an earlier branchy kernel (43,253,957)
and 00db6aff's kernel (43,253,483) (2^16 walks and adders, input set 1,
plain round 2). The organizer image's 46,832,918 for the predecessor
round differs by 27 (toolchain build). A METHOD_ALLOWANCE of 1,024
instructions per round is charged on every measured round for method and
toolchain differences.

**The implicit-previous-key kernel** (whose fixed path this revision's
kernel equals exactly, section above) was built in the organizer's reference-count
image (GCC 14.2.0-19, the same flags, `-static`) and counted by the same
QEMU method, filter ranges from each binary's symbol table. Every taken
correction adds exactly its extra instructions to the meter, so a round
without a rare event costs exactly FIXED + (meter increase); an escape
adds its rare path, priced by the rare reserve. Measured (round =
cumulative count at n less at n-1 rounds; FIXED = round minus meter
increase; clean = no escape in the round):

    B=2^20, A=2^20, input set 2, round 5 (clean)        656,973,433
    same build, input set 1, round 7 (clean)            656,973,433
    eight-copy kernel, same image, both rounds          658,022,010
    1f942c41's kernel, same image, set 2, round 5       659,594,869
    B=2^16, A=2^16 and 2^20, plain round 2, both sets   41,098,871
    eight-copy kernel, same settings and sets           41,164,409
    B=2^16, A=2^20, anchor round 64 (clean), both sets  41,361,018
    eight-copy kernel, same round and sets              41,492,092

Meters and escapes are equal in every pair. Against the eight-copy
kernel a plain round falls by 65,538 = (2^16-1)+3 at B=2^16 and
1,048,577 = (2^20-1)+2 at B=2^20: one per backward block step and 3 or 2
in the once-per-round compiled C code (register allocation differs by
build). 1f942c41 charged 659,594,897, 28 above this image's count of its
round (c2b498ec's declared 665,886,346 less its exact saving); the
eight-copy revision charged that less its exact 3B/2-5 = 1,572,859, i.e.
658,022,038. The implicit-previous-key revision charged that less the
saving measured at B=2^20, keeping the offset of 28, and this revision
keeps that charge:

    PLAIN_FIXED = 658,022,038 - 1,048,577 = 656,973,461 instructions,
    per walk step 626.538717 (rounded up),

so with the allowance PLAIN = (656,973,461+1,024)/2^20 = 626.539693 per
walk step (rounded up), against 627.539694 in the eight-copy revision,
629.039690 in 1f942c41, 635.039683 in c2b498ec and 714.613617 in the
predecessor. On those counted rounds (previous kernel, uniform adders)
the meter lay between 22.409 and 22.422 per walk step; this revision's
meter is measured in smallkey.json. Anchor rounds: section "Common path charged by the anchor schedule".

**Callgrind (1f942c41's kernel).** The organizer image's Valgrind 3.24.0
Callgrind (`--collect-atstart=no`, the harness's toggles) on 1f942c41's
B=2^16, A=2^16 build, input set 1: Ir = 42,730,422 for one round and
85,459,279 for two, so round 2 costs 42,728,857, equal to the QEMU count
of that binary (90,122,446 - 47,393,589); this kernel's rounds are QEMU
counts.

### Blinded corrections and their meter

Each site j has extra cost w_j when taken and an idealized probability
q_j for a uniform point and a uniform prefix: borrow 1/2 for dy, the first x3
subtraction, t and y3, below 2^224/p for D and the second x3 subtraction
(small-key adders), fifth word set 1/4 for the products (PRE, iw, inv, lambda,
y3') and 1/3 for the square, and below 2^-31 for each slow block (r3 >= p3
needs T in an interval of length below 2^224 + 2^192, and T has density
below 1/p). The walk step is charged

    CORR = sum_{j not in Z} w_j (q_j + 2^-10) + sum_{j in Z} w_j q_j + 30*2^-24 + 80*2^-31 = 17.458661 instructions,

Z = {dy, D, second x3 subtraction}, against 17.40 to 17.43 per walk step measured
on this kernel (smallkey.json), below CORR.

**What a lane's count depends on.** Let pi = PRE[w] mod p, the lane's
prefix residue. With inv canonical, everything a lane computes in a round
is a function of its walk point W (with its adder) and of the integer
held for PRE[w], which is pi unless it is pi + p: D, PRE[w+1], inv =
1/(pi*D) up to a fixed factor, iw, dy, lambda, the square, x3, t, y3' and
y3. Write c(pi, W) in [0,80] for the lane's correction count when PRE[w]
is held as pi. The value pi + p (below 2^256) needs pi < 2^256 - p <
2^224, so for a uniform pi it occurs with probability below 2^-31 and
changes the count by at most 80: the charge 80*2^-31 per step covers it.
The walk itself (x3, y3, keys, signs) does not depend on pi or rho.

**Premise H-trail (average only).** For a trail started from fresh
uniform coins, with the table fixed and any slot and starting round, the
expected value of sum over its batch steps of g(W) = E_{pi uniform}[c(pi, W)]
is at most CORR times its expected number of processed steps. This is a
statement about average behavior only. Support: for a uniform point and
prefix, uniform canonical operands borrow with probability (p-1)/(2p); for
D = x(W)-x(R) and u - x(R) a borrow needs the uniform operand below
x(R) < 2^224, probability below 2^224/p for every adder, so no table-mean
pad is needed (the predecessor's table-mean Hoeffding term DELTA_TABLE =
4 exp(-2A 2^-20) is kept in DELTA_CORR although no site now uses it); dy borrows with probability exactly 1/2 by sign
symmetry (L-ori) with no table-mean term, so dy carries pad 0; modelling the reduction multiple m = -ab p^-1 mod
2^256 as uniform gives P(T >= 2^256) = E[ab]/2^512 (1+O(2^-32)) <= 1/4 (1/3
for a square); the slow blocks are below 2^-31. Along actual walks the
measured per-step mean is 17.40 to 17.43 on this kernel (smallkey.json, with
per-site frequencies; `branchy.json` keeps c2b498ec's, 22.409 to 22.422
per step for the previous kernel). H-trail bounds only an abort probability, never the work.

**Meter and cap.** The blocks add each taken fix-up's extra count to x27;
each round adds the block totals to one 192-bit run-wide meter M, inside
the per-round control allowance. After round k, with n_k processed steps
so far, the run aborts if

    M > n_k*CORR + LAMBDA_CAP,
    LAMBDA_CAP = 2^112 + 2*(80+18)*K*TRAIL_MAX + 2^96.

So total correction work is at most S*CORR + LAMBDA_CAP plus one round's
overshoot K*80 (80 = most correction instructions in one step), all
charged: CORR inside the step charge and the rest (2^102.05 units) as
`correction_meter_slack`; the blinding draws and products are charged in
the per-round control (8*2^16 per global round spread over the 2^12
batches, plus 2^11 per batch round). Work is bounded
even if H-trail fails.

**Abort probability without a walk model.** Write M = sum_s g(W_s) +
sum_r N_r, where s runs over processed batch steps and N_r = sum over
round r's lanes of (c(pi, W) - g(W)) is round r's blinding noise.

* *Blinding noise.* The walk, hence every W and every PRE up to the
  factor rho, is independent of the rho draws. Given the whole walk, each
  global round's rho is a fresh uniform unit, so each lane's pi = rho*P[w]
  is uniform and E[N_r | walk] <= 2^20*80*2^-31 (the non-canonical
  prefix case); distinct rounds have independent rho, and the noise of
  a round over all batches is at most 2^32*80. Hoeffding's maximal
  inequality over at most 2^95.49 rounds gives P(max partial sum >=
  2^96) < 2^-1000.
* *Trails.* A trail's batch steps are a deterministic function of its
  start coins and the fixed table (the walk, its escapes, its
  distinguished-point and cap rules), and so is its sum T_i of g. Trails
  start from fresh, independent coins, so given the table (and each
  trail's slot and starting round, which fix its lane and anchor phase)
  T_i and its length len_i are independent of all earlier trails, whatever
  the walk's internal dependence; merges and cycles do not change this.
  Every slot finishes its current trail, a trail has at most TRAIL_MAX
  steps, and at most K_TRAILS = DP_MAX + 2K trails start, so D_i = T_i -
  CORR*len_i lies in an interval of width (80+18)*TRAIL_MAX and, by
  H-trail, E[D_i] <= 0 given earlier trails. Ordering trails by start, the
  Azuma-Hoeffding maximal inequality gives P(max_n sum_{i<=n} D_i >=
  2^112 - 2^109) < 2^-2744 (exponent 2(2^112-2^109)^2/(K_TRAILS*(98*TRAIL_MAX)^2)
  bits, with the drift debits below already deducted). At most K trails are in progress at a check;
  they move the sum by at most 2*(80+18)*K*TRAIL_MAX.
* *Drifts.* Two terms are not centered: the non-canonical-prefix case adds
  at most 80*2^-31 per step in expectation (below 2^103 over the run) and
  the per-round cap increment floor(2^20*CORR) loses below 1 per batch
  round (below 2^108); both are debited from the 2^112 trail slack.

So the cap is exceeded only on these two events or when a round's eight
draws all fail (probability below 2^-160 over all rounds):

    DELTA_CORR = 2^-2744 + 2^-1000 + 2^-160 + DELTA_TABLE <= 5.07e-14,

subtracted from the success bound. The independence used is that of the
algorithm's own coins (trail starts and rho), not of walk states; only the
average H-trail is a premise.

### Common path charged by the anchor schedule

An anchor is set only in rounds t of the global batch round counter with
`t % 64 == 0` (ANCHOR_PERIOD=64, `main` runs `batch_step(t)` for
t=1,2,..., all walks of a batch share t). Anchor rounds run the same
backward field block with the anchor tail, whose common path has 16
instructions against the plain tail's 12 in the compiled kernel (the
period counter, loaded with the anchor by one `ldp`, is tested (2) and
incremented (1), and one `stp` stores the new anchor and counter), and
the same forward block. Measured on the implicit-previous-key kernel at 2^16 walks
and 2^20 adders, the fixed path of a clean anchor round (t=64) is
41,361,018 instructions against 41,098,871 for a plain round on both
input sets: 4*(2^16-1) + 7 (the eight-copy kernel: 41,492,092 against
41,164,409, 5*(2^16-1) + 8). The ledger charges an anchor round
4*(2^20-1) + 9 more than a plain round, plus the same allowance.

Of the rounds t=1..R of a slot, R = L+TRAIL_MAX, at most ceil(R/64) are
anchor rounds. This is a schedule count, not an expectation, so every
run's common path costs at most

    CONTROL = K*(256*L + 1024*TRAIL_MAX)/2^20
    COMMON = S*PLAIN + CONTROL + K*ceil(R/64)*(ANCHOR - PLAIN)

instructions, C_COMMON = COMMON/S = 626.604513 (exact rational in
`accounting.py`), plus the metered corrections, CORR per processed step.
Replays are unbatched and priced separately; their per-trail anchors do
not enter this count. The orientation enters only through the table index
and the sign bit, so the fixed path costs the same for all coin values;
only the corrections depend on the data, and they are metered.

### Production counter widths and control allowance

The C file is a finite measurement harness (signed `int` loops, a 64-bit
round argument, a diagnostic `ring_n`). The specified algorithm uses
unsigned 128-bit batch-round and replay-length counters and 192-bit
instruction meters and work counters. The per-trail
period counter `tper` is a single 64-bit word: `measurement/step.c`
increments it once per 64 processed steps and compares it with
CAP_PERIODS=2^60, so the cap TRAIL_MAX=64*(2^60+1) steps is enforced
through that 64-bit word, which holds at most CAP_PERIODS+1=2^60+1 (61
bits), and the step count is never materialised in a counter. The anchor
total ceil(R/64)~2^89.5 is analytic, not a runtime counter, and the
diagnostic escape ring is omitted from the attack (no test reads it).

Charge 256 extra scalar instructions per full main-phase batch round and
retain 1,024 for every finishing round. Only main rounds are repriced;
all B=2^20 lanes are active there. Finishing semantics, rare-event and
replay bookkeeping retain their existing allowances.

The production wrapper stores completed t, L and R_MAX as two-limb values
at offsets 0,16,32 of state x19. An AAPCS call preserves x19. This listing
bounds recurring main control (scratch x9..x12 is reloaded each round):

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

The correction meter adds, per round (x20 = the round's block meter
totals, the 192-bit meter M at [x19,#48], the cap at [x19,#72], its
per-round increment floor(B*CORR) at [x19,#96]):

    ldp x9,x10,[x19,#48]; ldr x11,[x19,#64]
    adds x9,x9,x20; adcs x10,x10,xzr; adc x11,x11,xzr
    stp x9,x10,[x19,#48]; str x11,[x19,#64]
    ldp x12,x13,[x19,#72]; ldr x14,[x19,#88]; ldr x15,[x19,#96]
    adds x12,x12,x15; adcs x13,x13,xzr; adc x14,x14,xzr
    stp x12,x13,[x19,#72]; str x14,[x19,#88]
    subs xzr,x12,x9; sbcs xzr,x13,x10; sbcs xzr,x14,x11; b.lo abort

The normal path has 15+19 instructions. The full counter is kept in state;
passing its low limb preserves t modulo 64 and hence the anchor schedule.
The t bound below 2^97 excludes overflow. The wrapper does not touch the
kernel or its register allocation. Even allowing 32 register save/restore pairs,
6 stack/flags operations, 20 address/setup instructions and 32
phase-handoff operations gives at most 124 with the meter lines; charge 256 per main round. The final
main comparison and transition are included in that envelope; L is positive.
No per-lane finishing check is moved into this bound. The 1,024 finishing
allowance and all event-only wide-counter prices are retained. This is
an explicit production-control specification, not a new compiled measurement.

The record counter reaches DP_MAX=2*floor(S/2^60)+2*K, below 2^68.492
records; the trail counter DP_MAX+2*K (K starts, at most DP_MAX+K
restarts); the table index spans C=2*DP_MAX, below 2^69.492 slots. Each is
a 128-bit integer whose control is charged inside the table operation. `accounting.py` charges DP_MAX table operations at
DP_COST=PROBE_CAP*PROBE_COST+DP_INDEX_COST+DP_COUNTER_COST+DP_OVERFLOW_COST
=2^10*2^11+2^11+8+8=2,099,216=2^21+2^11+16 each. A probe's two-limb increment and comparison cost at most 8 of its 2^11.
Once per table operation, the index reduction (a long division) costs at
most 2^11, the two-limb record-counter update 8 and the overflow test 8.

The home slot is (key(W)>>84) mod C, above the 24 adder and 60
distinguished bits. As p = -1 mod
2^84, a distinguished key's bits from 84 up are uniform on M=(p+1)/2^84
values under M1, so a record's slot is within total variation
rem(C-rem)/(C*M) < 3.0228e-32 of uniform, rem = M mod C. The walk reads a
distinguished key only below bit 84 and in exact equality tests (class
identity), and lookups match the full key, so coupling each recorded
class's slot with an independent uniform one changes only the probe-cap
event, with probability at most DELTA_HASH = 2*DP_MAX*(C+1)/2^172 <
1.1509e-10, which also covers M=2^172.
The K initial starts and
the DP_MAX+K restarts are charged K*START and (DP_MAX+K)*START,
START=255,184,896, so the
trail counter adds at most 8*(DP_MAX+2*K) against START*(DP_MAX+2*K), a
factor of 31,898,112 for each eight-instruction counter update. No counter is
incremented per processed step: `tper` advances once per 64 steps and the
batch-round counter once per round, inside the phase-specific allowances
stated above. With the itemised charge the exact ledger
total rounds up to the declared 2^126.857022 claim.

### The record table: layout, empty test and charged initialization

Slot i (0 <= i < C) lies in bank i>>40 at byte offset 128*(i mod 2^40) <
2^47, one 64-bit register; a probe forms bank and offset from the two-limb
index inside its 2^11. A slot is 16 64-bit words: w0..w3 the full key(W),
w4..w7 coin a, w8..w11 coin b, w12..w13 the trail's live step count n (<
2^67, for the 65n replay cut-off: one walk step per round, so n is the
batch-round counter less the start round START writes into the slot state,
a two-limb subtraction), w14 the occupancy tag, w15 zero. **Empty
test: w14 = 0** (a load and cbz, inside the probe allowance). An insert
stores the 16 words with w14 = 1, inside its operation's charge, and
nothing clears a slot, so a stored record never tests empty; a probe
finding w14 = 1 compares w0..w3 with the full key.

Before the walk each bank is zeroed, x0 its base, x1 its slot count (1 <=
x1 <= 2^40):

    1:  stp xzr, xzr, [x0], #16      // 8 times: w0..w15
        subs x1, x1, #1
        b.ne 1b

That is exactly 10 instructions per slot (8 general-register store pairs,
a subtract, a branch), plus at most 64 per bank to select it, set x0, x1
and advance the bank counter. With N_B=ceil(C/2^40)=754,913,127 banks,
`accounting.py` charges record_table_init = 10*C + 64*N_B =
8,300,357,602,202,811,136,168 instructions, below 2^62.847880 units, as
preprocessing, so it enters both totals. Loads and stores are
instructions in the cost model; no other transfer term applies.
Preprocessing is at most 2^81.961192 units; its rounding slack, over 2^59
units, also covers setup and the 2^40 table copy (16 stores per adder). The total stays 2^126.857022.

Every other array is first written by a charged step: the adder table
(both signs) by its accepted trial inside A*2^40*(START+2^10), whose 2^14
constant covers negation and stores;
a walk slot's state (at most 512 bytes) by START at each trail start; the
batch buffers (D, PRE, selections) each round before they are read,
inside PLAIN and ANCHOR; rare-path scratch inside the forced maxima; the
record, trail and batch-round counters, the rare-work and correction meters in setup,
inside the 2^24 input_validation allowance. The escape ring is not part
of the attack.

## Concentration for the dominant event

This uses the **conditional** fresh-source premise, not a measured
average rate. In the deferred-exposure filtration of predecessor-proof.md
Section 2.0 and Lemma E1(b), before the key and orientation of a globally
fresh source are exposed, its successor closes the two-cycle with
probability p2=1/(2A). Stop the exposure process at the first violation of
the good-event conditions, padding with zeros after it and after
termination: at most S trials, each indicator I_t of conditional mean at
most p2. Already exposed sources belong to the inherited residual
categories. The law is **not** conditioned on G: the stopped process lives
on the original space, counts every fresh-source entry on G, and outside G
the separately charged DELTA_CUT applies.

For lambda>0, conditional expectation gives

    E[exp(lambda*I_t) | F_(t-1)]
       <= 1+p2*(exp(lambda)-1)
       <= exp(p2*(exp(lambda)-1)).

Iterating this inequality and applying Markov to the exponential yields,
with mu=S*p2 and eta=2^-16,

    Pr[sum I_t > (1+eta)*mu]
       <= exp(-mu*((1+eta)*ln(1+eta)-eta))
       <= exp(-mu*eta^2/(2+eta)).

Only the adapted conditional bound is used, and the last step is
(1+u)ln(1+u)-u >= u^2/(2+u), u>=0. The limit is n2=ceil((1+eta)*mu), the
reserve B2=n2*w2; the accounting checks the exponent exactly as a rational
above 100, so the tail, below 2^-100, is charged as delta2=2^-100.

This uses M1/L-ori more strongly: reduced-order measurements do not
establish the conditional law at production lengths, and a marginal rate
alone would not give the concentration, a material limitation.

## Residual reserve and a single enforceable meter

Write c2=B2/S. The expected residual event cost per processed step before
within-trail exceptional work is the sum of three parts, of which only the
first depends on A:

    E0 = (3/(4A^2)+2/A^3)*wL         longer entry, depends on A
       + 2^-63*256*(w2+wL)          false filter, A-independent
       + (tau_exact*TRAIL_MAX/S)*wL merged step, A-independent
    tau_exact = T^2/(2N).

The longer-entry term is charged the derived bound above; the merged term
retains its predecessor price. The false-filter term has a new generous
debit for speculative look-ahead. A false first-filter check exposes
at most one future class; a false anchor check exposes at most 65. Those
classes cannot subsequently be treated as conditionally unexposed for the
concentration proof. Charge the originating false-filter episode for its
direct work, at most 66 later exposed-source steps, and at most 128
persistent steps through an ensuing short cycle before escape. The debit
256*(w2+wL) exceeds that sum. Overlapping exposures may be charged
repeatedly, which only overcounts. An unclosed speculative path that
returns to visited classes merges with another trail (charged), closes a
short cycle (this debit and the 128-step bound) or violates G (failure);
otherwise continuing past its end exposes a new source. So no known
look-ahead key is treated as fresh. The runtime meters the same calls.

For the residual success argument, stop the analytical residual process
at the first G violation, just as for the two-cycle process. Let Eres=E0
and kM=400000. Set

    B_res = kM*Eres*S,
    B_rare = B2+B_res.

The factor kM is a free analysis parameter: raising it multiplies the
reserve kM*Eres per step and shrinks the abort allowance 1/kM. At A=2^24, E0
is at most 7.513393311e-9 (false-filter part 7.974745864e-11, merged part
4.252975230e-13). Tau and kM are tuned jointly; tau=4943041/10000000 with
kM=400000 kept a 1.0e-6 margin before the small-key revision (7.876e-7 now). Its reserve is 0.005494 instructions per step and its abort
allowance 1/400000. Its threshold
is public algorithm data.

All within-trail coincidences and excessive-length formal closures are
excluded by G and paid through DELTA_CUT; we do not bound their realized
costs by the average per-step allowance, nor use the predecessor's
within-trail fixed-point argument. On G the stopped residual process is
the actual residual cost, and its expectation is at most Eres*S by the
retained rate bounds plus the new speculative-exposure debit. Markov
gives a 1/kM = 1/400000 tail.

    Pr[rare work > B_rare and G] <= 2^-100 + 1/400000.

No independence between residual and two-cycle work is asserted. The
actual rare operations are metered at their worst-case prices, as in the
predecessor, with an abort before any priced operation that would exceed
B_rare, so work is bounded even when the premises fail. The decomposition
adds no runtime counter, freshness oracle or uncharged lookup, and the
profile-dependent parameters are public algorithm constants.

## Complete cost and success ledger

The ledger charges (COMMON+S*CORR+B_rare)/1000 for all walk steps, metered
corrections and metered rare work, then **adds**, rather than amortizes away:

* the correction meter's slack LAMBDA_AZ+LAMBDA_NF and one round's overshoot K*80;

* A adder constructions at 2^40 small-key trials each, and K initial trail
  starts, each trial or start with the complete capped scalar-sampling and
  two-scalar-multiplication price (trials also 2^10 for the key test);
* up to DP_MAX+K restart constructions, DP_MAX=2*floor(S/2^60)+2*K;
* DP_MAX table operations at 2^21+2^11+16 instructions each (probing, key
  comparisons, addressing, index reduction, counter update, overflow test);
* the record table's zeroing, 10*C + 64*ceil(C/2^40) instructions;
* two bounded coefficient replays at the original generous unbatched
  operation count, followed by solving, verification and serialization;
* input validation and setup allowances;
* the Lemma F closure series, summed for s < 4096 with the tail bounded by
  the consecutive-term ratio (e^2/8)(1+3/A)^2 < 0.924, so the truncation is an
  upper bound, equalling the fully enumerated sum at A=2^16 and A=2^20.

All processors' work is summed. No future-target amortization or free
advice is used. The storage formula is unchanged:

    2*DP_MAX*128 + K*512 + 2*A*2*64 + 2^24 bytes.

It covers records, walk states and batch buffers, both adder signs (and
the block's copy) and a fixed reserve. At A=2^24 it evaluates to
2^76.491736 bytes, dominated by the 2*DP_MAX*128 record table. `accounting.py`
evaluates these formulas and asserts that the total, rounded up at the sixth
decimal, equals the declared score 126868201/1000000.

The collision exponent penalty re-derived at the bound, with the predecessor's 2.6 margins:

    epsilon = 4/A + 2*(3/A) + 2*192*(3/(4A^2)+2/A^3)
            = 10/A + 288/A^2 + 768/A^3 = 5.960475e-7.

The conditional success lower bound is

    1-exp(-tau_exact*(1-epsilon))
      - 1/kM - 2^-100 - 2^-20 - DELTA_CORR - DELTA_ADDER_CLASS
      - DELTA_ABANDON - DELTA_CUT - DELTA_MERGE
      - DELTA_CAP - DELTA_PROBE - DELTA_HASH.

Each term is recomputed from the predecessor's formula at the new S and
T, and DELTA_CORR is the correction meter's abort term. The parameter
screen (kM from 1e5 to 2e6) of the earlier revisions took the smallest tau
on its 1e-7 grid keeping a 1.0e-6 margin above 0.39, with kM=400000 the
cheapest; this revision keeps that tau and kM. At tau=4943041/10000000 the
model term is 0.3900045649523908 and the exact bound is
0.39000078761012763, a margin of 7.876e-7 above 0.39; the named failure
terms at this S are DELTA_ABANDON=7.2116e-8, DELTA_CUT=2.5144e-7 (accidental
relations bounded with M_ACC in place of r, section "Small-key adders"),
DELTA_MERGE=1.2520e-14, DELTA_CAP=2.5477e-18 (adder exhaustion included),
DELTA_ADDER_CLASS=6.94e-32, DELTA_PROBE=2.4923e-44, DELTA_HASH=1.1509e-10
and DELTA_CORR=5.0657e-14. The bundled ledger asserts a margin of at least
5e-7 and the claimed cost. The margin is a conditional-model margin, not a
confidence interval for empirical trials; it is five times the
predecessor's 1.0e-7 floor, which only had to cover numerical error in
the recomputed bound (below 1e-14).

Per processed step the ledger charges C_COMMON=626.604513 (fixed path,
anchor schedule, control with the blinding draws), CORR=17.458661
(metered corrections) and the rare reserve 0.005494 (two-cycle 0.002489,
residual 0.003005), in total 644.068667 instructions (649.078434 in the
implicit-previous-key revision), against 714.759990 in the predecessor
(714.754495 + 0.005494). The correction meter's slack and overshoot add
2^102.05 units. The total is below 2^126.857022.

## Evidence scope, execution and limitations

Inherited forced rare paths, replica results and reduced-order records
keep their attribution and scopes. New: the generated blocks, their
generator, and in `evidence/measurement/branchy.json` the measured fixed
paths, anchor difference, site frequencies, per-round meter statistics
and exactness checks. None establishes a production success probability;
pooled restart points in the harness do not measure cold restarts, which
the attack pays in full.

The counts are exact executed-instruction counts, by emulation, of builds
made in the organizer's image; that image's Callgrind equals the method
on 1f942c41's 2^16-walk round, the method reproduces the predecessor's
Callgrind plain round within 27 instructions, c2b498ec's declared round
exceeds this image's count of it by 28, and 1,024 per round are charged
for such differences. The fixed path is the same for every input;
correction work is metered and capped, so H-trail moves only DELTA_CORR,
not the charged work.

One declared combined experiment uses the verifier's bounded
credential-free runner and returns two separately labelled results:
`accounting.py` recomputes the ledger, the correction charge, its slack
and abort terms, the concentration threshold, failure terms and
point-order certificate; `experiment.py`, the inherited complete
reduced-order procedure at the declared collision parameter on its
0.0001 grid, keys its table by the full key (no home slot, probe cap or
zeroing; the production table is specified above), and its finite-run
success is diagnostic only. It does not implement the production meters
or run a full-size attack; its successes, failures and replay checks
must be reported as observed. The assembly blocks cannot run in the
Python runner; their checks are recorded above. Parent files no executed
path reads are omitted (the predecessor's
`evidence/measurement/seed1.json`, `measurement/walk_ref.py`,
`summarize.py`, `reproduce.sh`, `run.sh`, `evidence/sign_law.py` and the
block test harness, whose results are recorded); `combined.py` reads
only `accounting.py` and `experiment.py`, and the evidence list names
exactly the files present. The forced-path maxima link to the
predecessor kernel's recorded sources (`measurement/forced/summary.json`,
stored without whitespace, content unchanged); the rare primitives they
price compile to the same instructions here.

The remaining material uncertainties are production-scale collision
mixing, conditional fresh-source orientation, the validity of the
fresh-window domination model, analytic handling of long trails and
rare returns, universal validity of the measured kernel bounds, H1 for
the small-key adder table, and the average correction premise H-trail
(abort term only).
The exponential-tail arguments do not resolve them; they are conditional
on them. Memory remains exponential and impractical. No claim of
practical compromise, unconditional improvement in P-256 security, or
human certification is made.
