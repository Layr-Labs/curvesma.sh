This revision deepens the accepted 127.005504 control-176-rare-overhead
crown (forward loop unrolled by four, checkpoint-only cycle detection,
rare-path overhead 256): same kernel, same measurements, same success
ledger, same premises. The only moved numeral is MAIN_CONTROL 176 -> 160
against the proof's own itemised 105 bound (55 headroom, factor 1.52).
The immediate material predecessor is submission
2db86473-3d13-41c7-bebd-fc6d1b6a1565, immutable commit
810d27968c8b67914ec4106a884881896ea8df77, directory
submissions/p256-ecdlp/control-176-rare-overhead, whose algorithm,
arguments and measurements are inherited unchanged. No prior platform
decision is an accepted premise, so the reviewer assesses them anew.

# P-256 sigma rho with a 2^24-adder table and a concentrated two-cycle reserve

## Current claim and provenance

For a uniform scalar x in 0..r-1, given Q=[x]P on the protected NIST P-256
profile, the following classical single-target algorithm has success probability
at least 0.39 under the disclosed H1, M1/L-ori and fresh-class exposure premises.
Its total charged computation is at most **2^127.005503** target-group-operation
equivalents. The normalization is C=1000 scalar AArch64 instructions per unit, not one
unit per P-256 addition. Peak memory is below 2^77 bytes, preprocessing
(record-table zeroing included) is at most 2^62.848100 units and is included
in the total, and nonuniform advice is zero bytes; instruction accounting
retains H3.

The exact declared score is 127005503/1000000 (127.005503), the ledger's total
rounded up at the sixth decimal. This is a conditional analytic bound, not an
unconditional complexity theorem or an executed full-size attack.

The material predecessor is submission
`b89d1d72-d9a1-4ef9-9ed7-c6ef0eec1ed3`, immutable commit
`0e1d63e9afd0fdb04b311b1ae72a2618271b371b`, directory
`submissions/p256-ecdlp/rho-sigma2`, whose algorithm, arguments and measurements
`predecessor-proof.md` reproduces with a historical banner. No prior platform decision is an accepted premise, so the
reviewer assesses them anew.

**Precedence:** the claim is this document and `submission.json`. The
predecessor's Sections 1.2--1.3 define the same algorithm except for A, L
and the rare-work threshold replaced below; Sections 2.0--2.4 and 2.8 provide
the unchanged probability-space, escape, replay, abandonment and merge
arguments, re-evaluated here at A=2^24; Sections 3--4 provide the forced-path
pricing, resource caps and storage; the batched step and its
measurement are replaced, and the record-table layout and zeroing added, below. This document replaces L, A
and every rate, threshold and derived count that depends on them, the table index of 8.2, the reserve
derivation of 2.5, the success evaluation of 2.6, the common-path charge of 3,
and every total. The two-steps-back filter of 1.3 step 5 and all that
depends on it are replaced by 'Checkpoint-only cycle detection'. The predecessor's score, collision parameter, reserve, success figure and
sensitivity statement are historical, not current claims; its
`evidence/accounting.json` stays in that commit.

## Algorithm, fixed costs and changed parameters

Use the protected p, a=-3, b, r and generator P without alteration. Handle
Q=identity by returning zero; otherwise validate the curve point. Let a
class be {W,-W}, key(W) the integer Montgomery representation of its affine
x-coordinate, and canon(W) the representative with even Montgomery y.
Let sigma(W) be the sign taking canon(W) to W. The walk is

    W -> W + sigma(W) R[key(W) mod A].

Here A=2^24 adders, each R_i=[c_i]P+[d_i]Q with fresh uniform scalar coins.
The table contains both signs. There are K=2^32 walk slots, arranged into
batches of 2^16 for simultaneous inversion. Distinguished points have
bits 24..83 of the key zero, so theta=2^-60 under the stated class model.
Each trail starts from fresh independent a,b with W=[a]P+[b]Q. It ends at
its first distinguished point, records it once with its coins, then starts
a fresh trail. Coefficients are reconstructed by replay only at a match.

The test ordering and replay algorithm are those specified in
predecessor-proof.md Sections 1.2--1.3; cycle detection is checkpoint-only
(below). The distinguished-point test precedes the cycle test at every
generated point, including look-ahead and escape points. A detected cycle
of length 2..64 escapes by doubling the canonical representative of its
least-key member. The 64-bit anchor filter only triggers exact checks.
Every trail is limited to TRAIL_MAX=64*2^60+64 processed steps; the first
key match is resolved by two replays, each limited to 65 times its original
live length. Degenerate matches and exceeded caps fail the run. Each
derived scalar is verified by [x]P=Q before output.

Set N=(r-1)/2+1 and

    tau_target = 4943241/10000000,
    T_min = ceil(sqrt(2*N*tau_target)),
    L = ceil(T_min/K), T=K*L,
    S = K*(L+TRAIL_MAX).

Every slot takes L main-phase steps and finishes its current trail. Thus S
bounds all processed steps, including persistent steps and finishing tails.
Sampling caps, maximum distinguished records, abandonment cap, probe cap,
restart costs, bounded replay, verification and output are unchanged
formulas, evaluated at this S. The point-order certificate proving
ord_r(2)>2*TRAIL_MAX is recomputed in `accounting.py`.

The common path costs at most C_COMMON=713.787745 instructions per
processed step under H3 (anchor schedule below). A fresh-source two-cycle
entry and a longer entry each cost at most wL=2,739,767 more (an anchor
check, its look-ahead and escape); false filter episodes with their
exposure debt are charged 256*(w2+wL), w2=82,750 being the former
two-cycle check, merged steps wL, the inherited forced-path/control allowances.
The runtime still uses **one total rare-work meter**: it need not
recognize whether a point is globally fresh or classify observed events
as fresh versus merged. That partition is solely for analysis.

### The adder table at A=2^24 and the measured step price

Every entry rate and reserve term is a function of A. At A=2^24:

    fresh-source two-cycle rate  p2 = 1/(2A) = 2^-25,
    longer-entry rate            <= 3/(4A^2) + 2/A^3 = 2.664536e-15,
    conditional entry bound      P' = 64/(2A) = 2^-19,
    collision exponent penalty   epsilon = 196/A + 288/A^2 + 768/A^3 = 1.168252e-5.

The two-cycle rate is exact under L-ori. The longer-entry rate is at most
3/(4A^2)+2/A^3 (Lemma E1(c)): E[w'] <= 3/(2A)+4/A^2 over windows of length
3..64, so E[w']/(2A) <= 3/(4A^2)+2/A^3, which the ledger charges and
`accounting.py` asserts equals (3/(2A)+4/A^2)/(2A).
Epsilon's closure margin charges every short-cycle entry, two-cycles
included, 192 wasted states: its two-cycle term is 2*192/(2A) and its
longer-entry term 2*192*(3/(4A^2)+2/A^3). The split-concentration threshold mu=S/(2A), the
two-cycle reserve ceil((1+eta)*mu)*wL, the residual expectation E0 and the
reserves built on them are recomputed, not rescaled: the two-cycle charge is
0.083141 and E0 is 7.51340e-9 instructions per step.

The A-dependence of the compiled step is checked on the binaries (organizer
flags, ADDER_BITS=16 and 24, disassembled): `batch_step` differs only in the
AND and shift immediates of the adder index and the two-limb
distinguished-point test, in branch offsets and in two `nop`s only the
2^24 build has: they pad the forward loop's head and run once per round,
so PLAIN and ANCHOR charge the measured 2^16 rounds plus 2. `step_rare`
has the same tests plus one instruction (101 against 100), and every
other function is identical up to address immediates. The 64-bit anchor
filter compares the whole low limb, so its false-match charge (2^-63 per
step, twice a 64-bit filter's 2^-64) is independent of A, and the table is data.

The larger table is charged, not assumed free: construction A*START =
4,281,292,120,129,536 instructions = 2^41.961183 units, storage A*2*64 = 2^31
bytes inside the peak-memory formula, and each 128-byte record's read-back
(16 words) lies inside START's 2^14 constant term. Adders and initial
starts are 2^49.966808 units of preprocessing; with the record-table zeroing
below, preprocessing is at most 2^62.848100 units, included in W.

### The fused batched step

`batch_step` does the predecessor's field operations in the same order,
each subtraction inside the product or square that uses it, values
passed in registers (`measurement/p256_fused.h`): `k_submul`
d = x(sR)-x(Wc), PRE[w+1] = d*PRE[w]; `k_mul2` iw = inv*PRE[w], inv = inv*d;
`k_lxy` lambda = (y(sR)-y(Wc))*iw, x3 = lambda^2-x(Wc)-x(sR),
y3 = (x(Wc)-x3)*lambda-y(Wc), with sR = sigma(R_i). The kernels advance
the D, PRE and Wc pointers by load/store writeback; inv, iw and lambda
share one stack block. The product is k_mul's
instructions (its two final constants in x7, x8), the square k_sqr's and
each subtraction k_subb_bl's, except that p's top limb takes one `and`
and the products of k_submul and k_mul2 end in a carry-only finish: for
factors below 2^256 the Montgomery sum T is below 2^256+p, so T-c*p (c its
fifth limb) is congruent and below 2^256. PRE[2..B-1], iw and inv feed
only products and may lie in [p, 2^256). The last prefix product, acc =
D[B-1]*PRE[B-1], is k_subb_bl then k_subm. Every other product has a
factor below p, so T<2p and k_mul's finish is canonical: D, PRE[1], acc,
lambda, x3 and y3 equal the previous kernel's (two products swap operands;
k_mul's a*b/2^256 mod p is symmetric), and finv, fis_zero, k_mul and k_sqr
read only canonical values. The kernels have no branch. The step's key,
selection and tests stay in registers: one branch, taken when the
distinguished-point test ((key>>24)|(x.v[1] & mask) == 0), the anchor
filter or, in an anchor round, the cap test fires, calls `step_rare`, the
C tail (same tests, rare calls and updates, same order); any other step
makes the updates that tail makes when no test fires. `rare_anchor`
(f_single, escape, same_class and point_is_dp inlined) and `add_single`
compile instruction-identically to the previous kernel's up to address
immediates, so `measurement/forced/summary.json` still prices them.
`step_rare` has 100 instructions (101 at A=2^24) and an acyclic
control-flow graph whose longest path is 69 (70). The rare-path OVERHEAD
is 256: any single event's own 70 plus an AAPCS call frame and address
setup, charged once per event in W_CHECK and once per look-ahead step in
W_ENTRY whose in-loop steps carry no frame (factor 3.7 over the longest
path). The 1024 finishing allowance and all event-only wide-counter
prices are retained.

### Canonical-multiplicand rows

Three products have a register factor a < p (a subtraction's output):
d in `k_submul`, and in `k_lxy` the y-difference and x(Wc)-x3. Their
Montgomery rows (`FC` in `p256_fused.h`) keep no sixth carry word. With a < p, any 4-limb b and V_0 = 0, the row values
V_{i+1} = (V_i + a*b_i + m_i*p)/2^64 satisfy V_i < 2p by induction
(V_i + a*b_i + m_i*p < 2p + 2(2^64-1)p = 2^65*p), and
U_i = V_i + a*b_i < p(2^64+1) - 2^64 + 1 < 2^320 (p < 2^256-2^223). So
U_i and each partial sum of its terms fit five words, no carry capture
wraps, the reduction's carry rides in the top word (hi3 + 1 < 2^64), the
windows equal the previous rows' values and the finishes are unchanged: 97 instructions per product instead of 100.
`k_mul2`'s first product keeps the previous rows (inv may be at least p);
its second takes them below. Per walk step 9 fewer (`k_submul` -3, `k_lxy`
-6).

### Redundant loads and constants

Three cuts change no value. `k_lxy` squares lambda in x9..x12, where its
product leaves it, instead of reloading it (-2). Each canonical finish
built 2^32-1 and 2^32-2 with two moves; `batch_step` stores them once per
call and `k_lxy`'s three finishes load them with one `ldp` (-3).
`k_mul2`'s second product takes D[w] < p (a subtraction output) as
register factor and inv as b, so the rows above apply (-3, plus 2 loads
of D[w]). Montgomery's m is the unique m < 2^256 with ab + m*p = 0
mod 2^256, so T = (ab + m*p)/2^256, its carry-only finish and inv's bits
do not depend on operand order.

Checks (organizer image and flags, on the build that introduced these
cuts, which kept the two-steps-back filter): 200,000 cases each of `k_submul`,
`k_mul2` (D below p) and `k_lxy` on two seeds (inputs up to 2^256-1;
a = 0, 1, p-1, p-1-2^64) match a bignum reference and the previous
kernels' outputs; the rows applied with `k_mul2`'s inv as register factor
failed 119,264 of 200,000, so the test sees the precondition. That build's
step dumps then equalled the earlier kernels' on six builds: seed 1, 8 rounds; 64-round
stress; DP 2^-55; 2^4, 2^6 and 2^10 adders with escapes, long cycles,
look-ahead records and abandonments.

### The hot backward step as one block

In a plain round (t not divisible by 64) the backward phase runs steps
B-1..1 in one assembly block, `BLOOP_FIELD` in `measurement/step.c`,
called from `batch_step`. The block keeps the running inverse inv in x19..x22, the
two finish constants 2^32-1 and 2^32-2 in x28, x29, the current
DP record's address in x23, the input point's y pointer &Wc[w].y in x24,
the output point's y pointer &Wn[w].y+16 in x25, the lane record's
&L[w].anchor in x26 and the table T2 in x27, so an iteration performs the
C loop's field operations without reloading them: iw = inv*PRE[w]
(PRE[w] at [x23+32]), inv = D[w]*inv (D[w] at [x23], then x23 -= 64),
lambda = (y(R)-y(Wc))*iw, the square lambda^2, x3 = S - x(Wc) - x(R),
t = x(Wc) - x3, y3 = t*lambda, y3 -= y(Wc), with R = L[w].sel read at
[x26+8]. The two products with a subtraction output as register factor
take D[w] and t (both below p), as in the previous kernel; the canonical
finishes use the blocked phase's constants in x28, x29; lambda and y3
use x0, x2, x18, x30, and x3, x4 hold x3[0], x3[1] across the final
product. x29 is saved on entry and x28, x29 are set at every entry,
because the rare path calls C code.

The rare path leaves the block with x0 = 1 before the lane update. The C
code recovers w from the saved lane pointer (masked by B-1), calls the
unchanged `step_rare(w, &Wn[w], Wc, 0)`, steps the lane pointer back and
re-enters the block unless w == 1; the remaining step w = 0 and every
anchor round keep the C loop. The block therefore makes exactly the C
loop's updates: on the common path it stores the anchor (unchanged) and
the selected adder and passes the loop's pointer test; on a rare step
`step_rare` makes them. Every input it reads is the value the C loop
reads, every output is stored to the same field, and the field text is
the measured kernel text, so block and C loop compute the same walk. When
the block was introduced (the edec9eab kernel) its dumps equalled
the C loop's on six builds.

### Checkpoint-only cycle detection

The lane record is {anchor, sel}: the two-steps-back filter (prev, cur)
and `rare_2cycle` are gone, and each step compares its new key's low limb
with the anchor filter only.

*Detection.* The anchor is set to every walk's new point in rounds t = 0
mod 64. A walk entering a cycle C of length l, 2 <= l <= 64, is anchored
on a member of C within 64 steps and returns to it within l more, where
the filter matches. `rare_anchor` then iterates f from the new point,
testing each look-ahead point for the distinguished property first,
until the new point's class returns (exact check) or 64 look-ahead steps
pass. On return it holds all l members and escapes to 2*canon(M), M the
least-key member of the whole cycle, so the escape depends only on C's
set of classes, not on the anchor's phase. A 4-cycle A->B->C->D->A
anchored at B (A least) returns to B, the look-ahead lists B, C, D, A,
and the walk escapes from A. A two-cycle is the case l = 2.

*Lemma R and D1.* An escape still follows an exact check over every
member of a closed cycle of length <= 65, all tested and not
distinguished, so Lemma R holds as written and its 2-cycle clauses are
vacuous. Every cycle of length <= 64 is caught under any schedule that
anchors each walk every 64 of its steps, so DELTA_CUT is unchanged. A
step makes at most one exact check, hence at most one escape: K_ESC =
TRAIL_MAX and D1 stand.

*Persistent steps (E2).* After an entry into a cycle of length l <= 64,
two-cycles included, the walk is anchored on the cycle within 64 steps
and detects the return within l more: at most 64 + l <= 128 persistent
steps. So processed <= fresh + 128*entries, F_m >= m - 128*E_m and
DELTA_ABANDON hold unchanged; P' = 64/(2A) already counts the l = 1 window.

*Replay.* A replay meets the same cycles (Lemma R). On a cycle of length
l it takes at most 64 + l walk steps, while the live trail took at least
l steps on it (the steps producing l-1 members and the return whose
point the escape replaced); on G distinct cycles of a trail share no
class (D2). So a replay of n live steps takes at most 33n walk steps,
inside the 65n cut-off, each with at most 68 operations (one evaluation,
66 for an anchor check, one escape): 33*68 < 2^13 = REPLAY_OPS per live step.

*Collision exponent.* A short-cycle entry wastes at most 64 + l <= 128
persistent steps and 64 look-ahead points, within the 192 states per
entry that epsilon's closure margin charges, with the collision-exponent
factor 2. The two-cycle term moves from 2*(3/A), when a two-cycle
escaped at entry, to 2*192/(2A): epsilon = 196/A + 288/A^2 + 768/A^3.

*Cost.* A two-cycle entry is charged wL, the anchor check with 64
look-ahead steps and a 65-cycle escape, in the meter row and in B2 =
n2*wL. The false-filter row keeps the former two-filter debit
2^-63*256*(w2+wL). Persistent steps are processed steps inside S.

*Kernel.* Per step the hot block fell from 591 to 588 instructions (one
load, one store and one compare) and the anchor round's C loop from 621
to 615; `new_trail` and `dp_restart` are no longer than before. Walks
differ from the edec9eab kernel's by design: two-cycles escape at the
next anchor, not at entry. In six builds (seed
1, 8 rounds; 64-round stress with 6-bit filters; DP 2^-55; 2^4, 2^6 and
2^10 adders, with look-ahead records, escapes and abandonments) and in 16
production rounds, `walk_ref.py check`, revised to the same rule, matches
every final point, event count and record digest
(`evidence/measurement/trail-paths.json`, `replay16-production.json`).
At 2^4 and 2^6 adders, where a two-cycle starts every 32 or 128 steps,
the persistent steps show (more trails reach the cap); at A=2^24 their
cost is the epsilon and 128-step terms above. Those dumps were taken of
the ebc47628 build; the next subsection's build writes the same bytes.

### Forward loop unrolled by four

The forward phase's `k_submul` calls for w = 1..B-2 keep their operands
and order, so every value is the same; `batch_step` makes four per loop
iteration (w = 1, 5, ..., 65529: 16,383 iterations), then w = 65533
and 65534. In the organizer build the loop body is the old body's field
and memory instructions four times (up to register names), and the loop
control, a compare and a branch per step (each load's writeback advanced
the selection pointer), becomes one add, compare and branch: 499
instructions against 4*126 = 504. The two remaining calls are emitted
straight-line, 248 against 252. A round executes 5*16,383 + 4 = 81,919
fewer instructions at A=2^16 and A=2^24 (the two `nop`s stay before the
loop head). `batch_step` has 3054 instructions (3056 at A=2^24); every
other function is instruction-identical to the ebc47628 build's up to
address immediates. Build md5s: 51255c70 (2^16), 130f8baa (2^24),
43caa286 (ANCHOR_PERIOD=2). Stdout and stderr equal the ebc47628 build's
byte for byte in the six builds above and in 16 production rounds
(`seed1.json` `dumps`), so the recorded `walk_ref.py` checks hold.

### Common path charged by the anchor schedule

In `measurement/step.c` an anchor is set only in rounds t of the global
batch round counter with `anchor_now = t % ANCHOR_PERIOD == 0`, where
ANCHOR_PERIOD=64 and `main` runs `batch_step(t)` for t=1,2,...; all 2^16
walks of a batch share t. The anchor-filter comparison and the other tests
of the common path run in every round and are inside the plain count. Only
the `else if (anchor_now)` branch is extra: it stores the anchor filter,
increments the period counter and compares it with the cap. Reaching the
cap is an abandonment, a rare event charged separately.

`evidence/measurement/seed1.json` gives whole-run Callgrind totals of that
file on production inputs for one batch of 2^16 walks, on two input sets:

* plain round, Ir(2)-Ir(1) = 93,501,988 - 46,751,004 = 46,750,984 on both
  sets (no step of rounds 1 and 2 enters `step_rare`): 81,919, the static
  count above, below the ebc47628 build's measured 46,832,903.
* anchor round, measured: the same source built with -DANCHOR_PERIOD=2
  differs from the production binary in one AND immediate (t & 1 for
  t & 63), so its round 2 is an anchor round. Ir(2)-Ir(1) = 95,271,435 -
  46,751,004 = 48,520,431 on both sets, with no rare event (stderr counts
  all 0, no unbatched operation): PLAIN + 27*(2^16-1) + 2, the C loop's
  615 against the block's 588 per step. An anchor round's only extra work
  is the anchor store, period increment and cap compare of every walk;
  a filter match or the cap is a rare event, charged separately.

At A=2^24 each round runs two more instructions (the forward loop's
alignment `nop`s), so PLAIN = 46750986/65536 = 23375493/32768 =
713.363433... and ANCHOR = 48520433/65536 = 740.363052... per walk step. `reproduce.sh`
measures Ir(8)-Ir(4) on seed 1 (four plain rounds) and the
ANCHOR_PERIOD=2 build's Ir(2)-Ir(1). Every conditional branch the common path executes per
walk is listed in `seed1.json`; each taken edge is a rare event, charged
separately. The orientation enters only through the table index, so a round
costs the same for all coin values. The next subsection's phase-specific counter
allowance is charged at every step.

A slot takes L main-phase steps and finishes its current trail, so its
batch runs at most R = L+TRAIL_MAX rounds and S = K*R. Of the rounds
t=1..R, at most ceil(R/64) are anchor rounds. This is a schedule count,
not an expectation, so every run's common path costs at most

    CONTROL = K*(160*L + 1024*TRAIL_MAX)/65536
    COMMON = S*PLAIN + CONTROL + K*ceil(R/64)*(ANCHOR - PLAIN)

instructions, and C_COMMON = COMMON/S = 713.787745 (rounded up; the exact
rational is in `accounting.py`). Replays are unbatched and priced
separately; their per-trail anchors do not enter this count.

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

Charge 160 extra scalar instructions per full main-phase batch round and
retain 1,024 for every finishing round. Only main rounds are repriced;
all B=65536 lanes are active there. Finishing semantics, rare-event and
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

The normal path has 15 instructions. The full counter is kept in state;
passing its low limb preserves t modulo 64 and hence the anchor schedule.
The t bound below 2^97 excludes overflow. The kernel and its register
allocation are unchanged. Even allowing 32 register save/restore pairs,
6 stack/flags operations, 20 address/setup instructions and 32
phase-handoff operations gives at most 105; charge 160 per main round. The final
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
1.1510e-10, which also covers M=2^172.
The K initial starts and
the DP_MAX+K restarts are charged K*START and (DP_MAX+K)*START,
START=255,184,896, so the
trail counter adds at most 8*(DP_MAX+2*K) against START*(DP_MAX+2*K), a
factor of 31,898,112 for each eight-instruction counter update. No counter is
incremented per processed step: `tper` advances once per 64 steps and the
batch-round counter once per round, inside the phase-specific allowances
stated above. With the itemised charge the exact ledger
total is below 2^127.0055030, at least 9.9e-8 bits below the
declared 2^127.005503 claim.

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
and advance the bank counter. With N_B=ceil(C/2^40)=754,928,399 banks,
`accounting.py` charges record_table_init = 10*C + 64*N_B =
8,300,525,520,567,896,568,456 instructions, below 2^62.847909 units, as
preprocessing, so it enters both totals. Loads and stores are
instructions in the cost model; no other transfer term applies.
Preprocessing is at most 2^62.848100 units; its rounding slack, over 2^40
units, also covers setup. The total stays 2^127.005503.

Every other array is first written by a charged step: the adder table
(both signs) by A*START, whose 2^14 constant covers negation and stores;
a walk slot's state (at most 512 bytes) by START at each trail start; the
batch buffers (D, PRE, selections) each round before they are read,
inside PLAIN and ANCHOR; rare-path scratch inside the forced maxima; the
record, trail and batch-round counters and the rare-work meter in setup,
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
reserve B2=n2*wL; the accounting checks the exponent exactly as a rational
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
debit for speculative look-ahead. A false anchor check exposes at most
65 future classes (the debit was sized for two filters). Those
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
and kM=420000. Set

    B_res = kM*Eres*S,
    B_rare = B2+B_res.

The factor kM is a free analysis parameter: raising it multiplies the
reserve kM*Eres per step and shrinks the abort allowance 1/kM. At A=2^24, E0
is at most 7.513393320e-9 (false-filter part 7.974745864e-11, merged part
4.253061269e-13). Tau and kM are tuned jointly; tau=4943241/10000000 with
kM=420000 keeps a 1.0e-5 margin. Its reserve is 0.003156 instructions per step and its abort
allowance 1/420000. Its threshold
is public algorithm data.

All within-trail coincidences and excessive-length formal closures are
excluded by G and paid through DELTA_CUT; we do not bound their realized
costs by the average per-step allowance, nor use the predecessor's
within-trail fixed-point argument. On G the stopped residual process is
the actual residual cost, and its expectation is at most Eres*S by the
retained rate bounds plus the new speculative-exposure debit. Markov
gives a 1/kM = 1/420000 tail.

    Pr[rare work > B_rare and G] <= 2^-100 + 1/420000.

No independence between residual and two-cycle work is asserted. The
actual rare operations are metered at their worst-case prices, as in the
predecessor, with an abort before any priced operation that would exceed
B_rare, so work is bounded even when the premises fail. The decomposition
adds no runtime counter, freshness oracle or uncharged lookup, and the
profile-dependent parameters are public algorithm constants.

## Complete cost and success ledger

The ledger charges (COMMON+B_rare)/1000 for all walk steps and metered
rare work, then **adds**, rather than amortizes away:

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

    2*DP_MAX*128 + K*512 + A*2*64 + 2^24 bytes.

It covers records, walk states and batch buffers, both adder signs and a
fixed reserve. At A=2^24 it evaluates to
2^76.491765 bytes, dominated by the 2*DP_MAX*128 record table. `accounting.py`
evaluates these formulas and asserts that the total, rounded up at the sixth
decimal, equals the declared score 127005503/1000000.

The collision exponent penalty re-derived at the bound, with the predecessor's 2.6 margins:

    epsilon = 4/A + 2*192/(2A) + 2*192*(3/(4A^2)+2/A^3)
            = 196/A + 288/A^2 + 768/A^3 = 1.168252e-5.

The conditional success lower bound is

    1-exp(-tau_exact*(1-epsilon))
      - 1/420000 - 2^-100 - 2^-20
      - DELTA_ABANDON - DELTA_CUT - DELTA_MERGE
      - DELTA_CAP - DELTA_PROBE - DELTA_HASH.

Each term is recomputed from the predecessor's formula at the new S and
T. The parameter screen takes the smallest tau on its 1e-7 grid at
which the conditional bound keeps a 1.0e-5 margin above 0.39. At
tau=4943241/10000000 the model term is 0.3900134218276925 and the exact bound is
0.3900100149682157, a margin of 1.0014e-5 above 0.39 (one grid step lower,
9.954e-6); the named failure terms at this S are
DELTA_ABANDON=7.2118e-8, DELTA_CUT=5.8544e-17, DELTA_MERGE=1.2521e-14,
DELTA_CAP=2.4394e-18, DELTA_PROBE=2.4925e-44 and DELTA_HASH=1.1510e-10. The bundled ledger asserts
p_claim >= 0.390005, the 1.0e-5 margin, the grid step and the claimed cost. The margin is a
conditional-model margin, not a confidence interval for empirical trials.

## Evidence scope, execution and limitations

The inherited field assembly, forced rare paths, replica results and
reduced-order procedure records are included with attribution, supporting
H3 and the model at their stated scopes; the plain and anchor rounds are
new measurements of this harness. `evidence/trails/` and
`evidence/trails-stress.json` were run with the predecessor's procedure,
which also had the two-steps-back filter: they support the collision
statistics and Lemma R at their scope, not this file's detector, which
the declared run, `evidence/lemma-r-sensitivity.json` and the six kernel
builds exercise. None establishes a production success probability. The predecessor's harness pools restart points; the
attack pays fresh coins at full cost, and pooled runs are not taken as a
measurement of cold restarts.

One declared combined experiment uses the verifier's bounded credential-free
runner and returns two separately labelled results. `accounting.py`
recomputes the new numerical ledger, concentration threshold, failure
terms and point-order certificate. `experiment.py` is the inherited
complete reduced-order procedure, with the declared collision parameter
0.4943, which is production tau=0.4943241 on the available 0.0001 input
grid, 2.41e-5 below it; it runs the checkpoint-only rule; its table is a dictionary keyed by the full key, so
it has no home slot, probe cap or zeroing (the production table is
specified above); its finite-run success is diagnostic only. It does
not implement the production instruction meter or execute a full-size
attack. Its actual successes, failures and replay checks must be reported
as observed.

The remaining material uncertainties are production-scale collision
mixing, conditional fresh-source orientation, the validity of the
fresh-window domination model, analytic handling of long trails and
rare returns, and universal validity of the measured kernel bounds.
The new exponential-tail argument does not resolve them;
it is conditional on them. Memory remains exponential and impractical.
No claim of practical compromise, unconditional improvement in P-256
security, or human certification is made.
