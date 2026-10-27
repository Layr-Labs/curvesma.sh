This revision derives from submission b75c69bb-d881-4d46-b801-e2f670ff7f34
(floor-1e7-joint-retune, candidate commit 77ccb326), whose kernel and
measurements come from 778628fb (anchor lane block, aburan28) and its
lineage 6c22dbfd, 427d9f07, edec9eab, 74599488, 214ec421, f9bd18fb,
dc8d1ebc and the material predecessor b89d1d72. It keeps that algorithm,
every rare-path price, every premise and the success-model form, and
changes the implementation of the common batched step: field corrections
become metered branches, each round's prefix chain is blinded, the forward pass becomes a
register-resident block, and batches hold 2^20 walks. It extends
submission 7169b17e (AI-validated at 126.927283), whose correction-cap
argument used the independence of trails' own start coins with an
average-only premise: here every correction is metered, and each round's
prefix chain is blinded by a fresh random unit so that the batch-dependent
corrections reduce to the same trail-level argument plus independent
per-round blinding noise.

# P-256 sigma rho with a 2^24-adder table and a concentrated two-cycle reserve

## Current claim and provenance

For a uniform scalar x in 0..r-1, given Q=[x]P on the protected NIST P-256
profile, the following classical single-target algorithm has success probability
at least 0.39 under the disclosed H1, M1/L-ori, fresh-class exposure and H-trail premises.
Its total charged computation is at most **2^126.896671** target-group-operation
equivalents. The normalization is C=1000 scalar AArch64 instructions per unit, not one
unit per P-256 addition. Peak memory is below 2^77 bytes, preprocessing
(record-table zeroing included) is at most 2^62.848071 units and is included
in the total, and nonuniform advice is zero bytes; instruction accounting
retains H3.

The exact declared score is 126896671/1000000 (126.896671), the ledger's total
rounded up at the sixth decimal. This is a conditional analytic bound, not an
unconditional complexity theorem or an executed full-size attack.

The immediate source is submission
`b75c69bb-d881-4d46-b801-e2f670ff7f34`, candidate commit
`77ccb326ae7f10bd72d1af3e769a379b8d4d3f2b`, directory
`submissions/p256-ecdlp/floor-1e7-joint-retune`; its kernel source is
`778628fb-264e-4865-9413-f934795125ab` (commit `1fcbf53a`). The material
predecessor is submission `b89d1d72-d9a1-4ef9-9ed7-c6ef0eec1ed3`, immutable
commit `0e1d63e9afd0fdb04b311b1ae72a2618271b371b`, directory
`submissions/p256-ecdlp/rho-sigma2`, whose algorithm, arguments and
measurements `predecessor-proof.md` reproduces with a historical banner. No
prior platform decision is an accepted premise, so the reviewer assesses
them anew.

**Precedence:** the claim is this document and `submission.json`. The
predecessor's Sections 1.2--1.3 define the same algorithm except for A, L
and the rare-work threshold replaced below; Sections 2.0--2.4 and 2.8 give the
unchanged probability-space, escape, replay, abandonment and merge arguments,
re-evaluated here; Sections 3--4 give the forced-path pricing, resource
caps and storage; the batched step and its measurement are replaced, and the
record-table layout and zeroing added, below. The predecessor's score and
parameters are historical, not current claims; its `evidence/accounting.json`
remains in the predecessor submission (omitted here for the size limit).

## Algorithm, fixed costs and changed parameters

Use the protected p, a=-3, b, r and generator P without alteration. Handle
Q=identity by returning zero; otherwise validate the curve point. Let a
class be {W,-W}, key(W) the integer Montgomery representation of its affine
x-coordinate, and canon(W) the representative with even Montgomery y.
Let sigma(W) be the sign taking canon(W) to W. The walk is

    W -> W + sigma(W) R[key(W) mod A].

Here A=2^24 adders, each R_i=[c_i]P+[d_i]Q with fresh uniform scalar coins.
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

The common path costs at most C_COMMON=639.619998 instructions per
processed step under H3 (anchor schedule below), and corrections at most
CORR=22.389326 more, metered (below). A fresh-source two-cycle
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
measured, B=2^20, ADDER_BITS=20 and 24, disassembled): `batch_step` has
2353 instructions in both and differs only in the `extr`, `ubfiz` and
`and` immediates of the DP test and the adder index (six instructions);
`step_rare`, `add_single`, `finv` and `dp_restart` are identical up to
immediates, and `rare_anchor` and `rare_2cycle` have one more instruction
at A=2^24, as in the predecessor's comparison, inside the per-event
OVERHEAD. The two 64-bit filters compare the whole low limb, so their
2^-63 false-match charge is independent of A, and the table is data.

The larger table is charged, not assumed free: construction A*START =
4,281,292,120,129,536 instructions = 2^41.961183 units, storage A*2*64 = 2^31
bytes, twice with the block's copy at 2^40, inside the peak-memory formula, and each 128-byte record's read-back
(16 words) lies inside START's 2^14 constant term. Adders and initial
starts are 2^49.966808 units of preprocessing; with the record-table zeroing
below, preprocessing is at most 2^62.848071 units, included in W.

### The batched step in this revision

The walk, its tests and every rare path are the predecessor's. Only the
implementation of the common batched step changes, in five ways, each
value-preserving:

1. **Branchy corrections.** In the predecessor kernel every conditional
   correction of the field arithmetic was computed branch-free and always
   executed. They are now forward branches over their fix-ups, and a
   taken fix-up adds its own extra instruction count to a meter (below).
2. **A register-resident forward block.** The prefix products PRE[w+1] =
   D[w]*PRE[w], w=1..B-2, run as one assembly loop (`FLOOP`) with PRE kept
   in x19..x22: the product's last reduction step writes the new PRE over
   the multiplier registers it has consumed, so PRE is never reloaded.
3. **Operand forms.** PRE and inv are kept canonical (branchy finish):
   inv as register factor lets iw = inv*PRE use the five-word rows instead
   of six-word rows, and a canonical PRE makes each lane's arithmetic a
   function of the prefix value alone. iw and lambda are left below 2^256
   (their consumers accept that); the products' last reduction steps write
   straight into their destination registers.
4. **2^20-walk batches** (per-round inversion and set-up amortized over
   16 times more walks) and a copy of the two-sign table T2 at virtual
   address 2^40, read by the blocks, so the adder address is one `orr`
   and x27 is free for the meter; the C code keeps the static T2.
5. **A blinded prefix chain.** Each batch round draws a fresh uniform unit
   rho (rejection sampling, at most 8 draws) and starts the prefix chain at
   PRE[1] = rho*D[0]; the C path for w=0 takes iw = inv*rho = 1/D[0]. All
   inverses, slopes and points are unchanged; every prefix PRE[w] is then a
   uniform unit given the walk, which the meter's abort bound uses.

`measurement/step.c` is the predecessor's harness with these changes; the
backward block `BLOOP_FIELD_B` (`measurement/bloop_branchy.h`) and the
forward block `FLOOP` (`measurement/floop.h`) are generated by
`measurement/gen7.py` from the predecessor's measured `BLOOP_FIELD` text:
every multiplication, square and reduction row is that text, and only
the corrections and the destinations of the last reduction steps differ.
The step tail, anchor tail, rare exits and rare primitives are unchanged:
`step_rare`, `rare_2cycle`, `rare_anchor`, `add_single`, `finv`,
`new_trail` and `dp_restart` disassemble instruction-identically to the
predecessor build up to address immediates, so the forced-path maxima of
`measurement/forced/summary.json` still price them.

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
* **Carry-only finish** (iw, lambda): T = r + c*2^256, c in {0,1},
  T < 2^256+p. `cbz c` over [subtract p; meter += 5]: 1 instruction (was 7
  or 9), 5 more when c=1; the result is T - c*p < 2^256.
* **Canonical finish** (PRE, inv, y3' = t*lambda; T < 2p): `cbz c,L`, else
  subtract p, meter += 4, `b done`; `L: cmn r3,x28; b.hs slow`. Common
  path 3 (was 9), c=1 costs 4 more. If c=0, T < p unless r3 >= p3
  (`cmn` sets C exactly then); the out-of-line slow block computes T-p,
  keeps it with four `csel` on no borrow, meter += 10, and returns.
* **lambda^2**, from lambda < 2^256: T < 2^256+p, so after a c=1
  subtraction T-p may still be at least p. The fix-up (meter += 5) falls
  through into the same `cmn/b.hs` re-test; one more conditional
  subtraction in the slow block always suffices (T - 2p < p).

Five-word rows need a register factor below p (D, dy, t, inv) and accept
any 4-limb multiplier; iw and lambda < 2^256 enter only as multipliers and
the square. Every exposed value (x3, y3, key, sign bit) is
canonical, so walk states, dumps and rare-path inputs equal the
predecessor's.

**Exactness checks** (`evidence/measurement/branchy.json`): (i) n=8 dumps
of 2^16 walks on two input sets, and stress builds with 2^4, 2^6 and 2^10
adders exercising every rare path, equal the predecessor kernel's byte
for byte; (ii) the whole backward block, wrapped as a function, equals an
exact big-integer model (Montgomery integers, every branch outcome, the
meter value) on 240,000 random and crafted inputs covering every path
(22,118 non-canonical lambdas, 10,063 to 41,984 cases of each slow path
and its subtracting branch, boundary operands 0, 1, p-1, p, 2^256-1);
(iii) the forward block equals the model on 6,000 runs of 62 steps
(4,607 slow-path cases and 2,085 of its subtracting branch).

### Measured fixed path and metered corrections

Instruction counts were taken with GCC 14.2.0 (Ubuntu cross build
14.2.0-4ubuntu2) and the organizer flags `-O3 -fno-tree-vectorize
-mgeneral-regs-only`, by exact guest-instruction counting under QEMU 8.2.2
user mode (translation-block sizes from `-d in_asm`, executions from
`-d exec,nochain`, restricted to the harness's own functions; identical to
`-one-insn-per-tb` counting; a 1,000-iteration three-instruction loop adds
exactly 3,000). Valgrind 3.22.0 Callgrind (arm64, `--collect-atstart=no`,
the harness's toggles) run under the same emulator gives the same round
differences: 46,832,945 for the predecessor kernel and 43,253,957 for an
earlier branchy kernel, equal to this method's counts. The organizer
image's 46,832,918 differs by 27 per round (toolchain build). A
METHOD_ALLOWANCE of 1,024 instructions per round is charged on top of
every measured round for the method and toolchain difference.

Because every taken correction adds exactly its extra instructions to the
meter, a round costs exactly FIXED + (meter increase): verified on every
clean measured round below; rounds containing a two-cycle escape exceed
it by that escape's rare path, priced by the rare reserve. The
predecessor kernel at the same settings costs 748,723,509 per round
(714.04 per walk step). At B=2^20 walks and 2^20 adders (A-dependence as
above):

    plain round fixed path  670,080,642 instructions, both input sets
    per walk step           639.038165
    meter per walk step     22.16 to 22.17 measured

so with the allowance PLAIN = (670,080,642+1,024)/2^20 = 639.039142 per walk step, against
714.613617 in the predecessor. The anchor block's tail is the
predecessor's (5 more instructions per block step than the plain tail);
measured at 2^16 walks on both input sets, the anchor round's fixed path is the
plain round's plus 5*(2^16-1)+8 (42,245,763 against 41,918,080); the ledger
charges 5*(2^20-1)+9 per anchor round.

### Blinded corrections and their meter

Each site j has extra cost w_j when taken and an idealized probability
q_j for a uniform point and a uniform prefix: borrow 1/2 for the six
subtractions, fifth word set 1/4 for the products (PRE, iw, inv, lambda,
y3') and 1/3 for the square, and below 2^-31 for each slow block (r3 >= p3
needs T in an interval of length below 2^224 + 2^192, and T has density
below 1/p). The walk step is charged

    CORR = sum_j w_j (q_j + 2^-8) + 40*2^-24 = 22.389326 instructions,

against 22.16-22.17 per step measured.

**What a lane's count depends on.** With PRE and inv canonical, everything
a lane computes in a round is a function of its walk point W (with its
adder) and of its prefix value pi = PRE[w]: D, PRE[w+1] = pi*D, inv =
1/(pi*D) up to a fixed factor, iw, dy, lambda, the square, x3, t, y3' and
y3. So the lane's correction count is c(pi, W) in [0,85]. The walk itself
(x3, y3, keys, signs) does not depend on pi or rho.

**Premise H-trail (average only).** For a trail started from fresh
uniform coins, with the table fixed and any slot and starting round, the
expected value of sum over its batch steps of g(W) = E_{pi uniform}[c(pi, W)]
is at most CORR times its expected number of processed steps. This is a
statement about average behavior only. Support: for a uniform point and
prefix, uniform canonical operands borrow with probability (p-1)/(2p); for
D = x(R)-x(W) and u - x(R) the adder average is 1/2 -/+ (mean x(R)/p -
1/2), within 2^-10 by Hoeffding over the A independent adders except with
probability DELTA_TABLE = 4 exp(-2A 2^-20); dy borrows with probability 1/2
by sign symmetry (L-ori); modelling the reduction multiple m = -ab p^-1 mod
2^256 as uniform gives P(T >= 2^256) = E[ab]/2^512 (1+O(2^-32)) <= 1/4 (1/3
for a square); the slow blocks are below 2^-31. Along actual walks the
measured per-step mean is 22.16-22.17 per step (`branchy.json`, with per-site
frequencies). H-trail bounds only an abort probability, never the work.

**Meter and cap.** The blocks add each taken fix-up's extra count to x27;
each round adds the block totals to one 192-bit run-wide meter M, inside
the per-round control allowance. After round k, with n_k processed steps
so far, the run aborts if

    M > n_k*CORR + LAMBDA_CAP,
    LAMBDA_CAP = 2^112 + 2*(85+23)*K*TRAIL_MAX + 2^96.

So total correction work is at most S*CORR + LAMBDA_CAP plus one round's
overshoot K*85 (85 = most correction instructions in one step), all
charged: CORR inside the step charge and the rest (2^102.05 units) as
`correction_meter_slack`; the blinding draws and products are charged in
the per-round control (8*2^16 + 2^11 per batch round). Work is bounded
even if H-trail fails.

**Abort probability without a walk model.** Write M = sum_s g(W_s) +
sum_r N_r, where s runs over processed batch steps and N_r = sum over
round r's lanes of (c(pi, W) - g(W)) is round r's blinding noise.

* *Blinding noise.* The walk, hence every W and every PRE up to the
  factor rho, is independent of the rho draws. Given the whole walk, each
  batch round's rho is a fresh uniform unit, so each lane's pi = rho*P[w]
  is uniform and E[N_r | walk] = 0; distinct batch rounds have independent
  rho, and |N_r| <= 2^20*85. Hoeffding's maximal inequality over at most
  2^107.49 batch rounds gives P(max partial sum >= 2^96) < 2^-1000.
* *Trails.* A trail's batch steps are a deterministic function of its
  start coins and the fixed table (the walk, its escapes, its
  distinguished-point and cap rules), and so is its sum T_i of g. Trails
  start from fresh, independent coins, so given the table (and each
  trail's slot and starting round, which fix its lane and anchor phase)
  T_i and its length len_i are independent of all earlier trails, whatever
  the walk's internal dependence; merges and cycles do not change this.
  Every slot finishes its current trail, a trail has at most TRAIL_MAX
  steps, and at most K_TRAILS = DP_MAX + 2K trails start, so D_i = T_i -
  CORR*len_i lies in an interval of width (85+23)*TRAIL_MAX and, by
  H-trail, E[D_i] <= 0 given earlier trails. Ordering trails by start, the
  Azuma-Hoeffding maximal inequality gives P(max_n sum_{i<=n} D_i >=
  2^112) < 2^-2951. At most K trails are in progress at a check; they move
  the sum by at most 2*(85+23)*K*TRAIL_MAX.

So the cap is exceeded only on these two events or when a round's eight
draws all fail (probability below 2^-148 over all rounds):

    DELTA_CORR = 2^-2951 + 2^-1000 + 2^-148 + DELTA_TABLE <= 5.07e-14,

subtracted from the success bound. The independence used is that of the
algorithm's own coins (trail starts and rho), not of walk states; only the
average H-trail is a premise.

### Common path charged by the anchor schedule

An anchor is set only in rounds t of the global batch round counter with
`t % 64 == 0` (ANCHOR_PERIOD=64, `main` runs `batch_step(t)` for
t=1,2,..., all walks of a batch share t). Anchor rounds run the same
backward field block with the predecessor's anchor tail, whose common path
has 19 instructions against the plain tail's 14 (load and test of the
period counter, its increment and store, and one more store pair), and the
same forward block. Measured at 2^16 walks and 2^20 adders, the fixed path of a clean anchor
round (t=64) is 43,425,155 instructions against 43,097,471 for a plain
round: 5*(2^16-1) + 9. The ledger charges an anchor round 5*(2^20-1) + 9
more than a plain round, plus the same allowance.

Of the rounds t=1..R of a slot, R = L+TRAIL_MAX, at most ceil(R/64) are
anchor rounds. This is a schedule count, not an expectation, so every
run's common path costs at most

    CONTROL = K*(256*L + 1024*TRAIL_MAX)/2^20
    COMMON = S*PLAIN + CONTROL + K*ceil(R/64)*(ANCHOR - PLAIN)

instructions, C_COMMON = COMMON/S = 657.117799 (exact rational in
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
total rounds up to the declared 2^126.896671 claim.

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
Preprocessing is at most 2^62.848071 units; its rounding slack, over 2^40
units, also covers setup and the 2^40 table copy (16 stores per adder). The total stays 2^126.896671.

Every other array is first written by a charged step: the adder table
(both signs) by A*START, whose 2^14 constant covers negation and stores;
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
kM=400000 keeps a 1.0e-6 margin. Its reserve is 0.005494 instructions per step and its abort
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

* A adder constructions and K initial trail starts, each with the complete
  capped scalar-sampling and two-scalar-multiplication price;
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
decimal, equals the declared score 126896671/1000000.

The collision exponent penalty re-derived at the bound, with the predecessor's 2.6 margins:

    epsilon = 4/A + 2*(3/A) + 2*192*(3/(4A^2)+2/A^3)
            = 10/A + 288/A^2 + 768/A^3 = 5.960475e-7.

The conditional success lower bound is

    1-exp(-tau_exact*(1-epsilon))
      - 1/kM - 2^-100 - 2^-20 - DELTA_CORR
      - DELTA_ABANDON - DELTA_CUT - DELTA_MERGE
      - DELTA_CAP - DELTA_PROBE - DELTA_HASH.

Each term is recomputed from the predecessor's formula at the new S and
T, and DELTA_CORR is the correction meter's abort term. The parameter
screen (kM from 1e5 to 2e6) takes the smallest tau on its 1e-7 grid at
which the conditional bound keeps a 1.0e-6 margin above 0.39; kM=400000
is the cheapest. At tau=4943041/10000000 the model term is
0.3900045649523908 and the exact bound is 0.39000103904670785, a margin of
1.0390e-6 above 0.39; the named failure terms at this S are
DELTA_ABANDON=7.2116e-8, DELTA_CUT=5.8542e-17, DELTA_MERGE=1.2520e-14,
DELTA_CAP=2.4393e-18, DELTA_PROBE=2.4923e-44, DELTA_HASH=1.1509e-10 and
DELTA_CORR=5.0657e-14. The bundled ledger asserts the 1.0e-6 margin and
the claimed cost. The margin is a conditional-model margin, not a
confidence interval for empirical trials; it is ten times the
predecessor's 1.0e-7 floor, which only had to cover numerical error in
the recomputed bound (below 1e-14).

Per processed step the ledger charges C_COMMON=639.619998 (fixed path,
anchor schedule, control with the blinding draws), CORR=22.389326
(metered corrections) and the rare reserve 0.005494 (two-cycle 0.002489,
residual 0.003005), in total 662.014817 instructions, against 714.759990 in the predecessor
(714.754495 + 0.005494). The correction meter's slack and overshoot add
2^102.05 units. The total is below 2^126.896671.

## Evidence scope, execution and limitations

The inherited forced rare paths, replica results and reduced-order
records are included with attribution at their stated scopes. New: the
generated blocks, their generator, and in `evidence/measurement/branchy.json`
the measured fixed paths, anchor difference, site frequencies, per-round
meter statistics and exactness checks. None establishes a production
success probability; pooled restart points in the harness are not taken
as a measurement of cold restarts, which the attack pays in full.

The counts are exact executed-instruction counts taken by emulation, not
by the organizer's Callgrind image; the method reproduces the
predecessor's Callgrind plain round to within 27 instructions, and 1,024
per round are charged for that. The fixed path is the same for every
input; correction work is metered and capped, so the charged work does
not depend on H-trail, which moves only DELTA_CORR.

One declared combined experiment uses the verifier's bounded credential-free
runner and returns two separately labelled results. `accounting.py`
recomputes the numerical ledger, the correction charge, its slack and abort
terms, the concentration threshold, failure terms and point-order
certificate. `experiment.py` is the inherited complete reduced-order
procedure, with the declared collision parameter on its 0.0001 input
grid; its table is a dictionary keyed by the full key, so it has no home
slot, probe cap or zeroing (the production table is specified above); its
finite-run success is diagnostic only. It does not implement the
production instruction meters or execute a full-size attack. Its actual
successes, failures and replay checks must be reported as observed. The
assembly blocks cannot run in the Python runner; their checks are the
recorded ones above. The package omits parent files no executed path reads
(the predecessor's `measurement/walk_ref.py`, `summarize.py`,
`reproduce.sh`, `run.sh`, `evidence/sign_law.py`, and the block test
harness, whose results are recorded); `combined.py` reads only
`accounting.py` and `experiment.py`. The evidence list names exactly the
files present. The forced-path maxima link to the predecessor kernel's
recorded sources (see `measurement/forced/summary.json`); the rare-path
functions they price compile to the same instructions here.

The remaining material uncertainties are production-scale collision
mixing, conditional fresh-source orientation, the validity of the
fresh-window domination model, analytic handling of long trails and
rare returns, universal validity of the measured kernel bounds, and the
average correction premise H-trail (abort term only).
The exponential-tail arguments do not resolve them; they are conditional
on them. Memory remains exponential and impractical. No claim of
practical compromise, unconditional improvement in P-256 security, or
human certification is made.
