# P-256 sigma rho with a 2^24-adder table and a concentrated two-cycle reserve

## Current claim and provenance

For a uniform scalar x in 0..r-1, given Q=[x]P on the protected NIST P-256
profile, the following classical single-target algorithm has success probability
at least 0.39 under the disclosed H1, M1/L-ori and fresh-class exposure premises.
Its total charged computation is at most **2^127.142273** target-group-operation
equivalents. The normalization is C=1000 scalar AArch64 instructions per unit;
an actual P-256 addition is not assumed to cost one unit. Peak memory is below
2^77 bytes, all preprocessing is below 2^50 units and is included in the total,
and nonuniform advice is zero bytes. The instruction accounting retains H3,
the predecessor's implementation-measurement assumption.

The exact declared score is 127142273/1000000 (127.142273), the ledger's total
rounded up at the sixth decimal. This is a conditional analytic bound under
an exploratory acceptance policy. It is neither an unconditional complexity
theorem nor an executed full-size attack.

The material predecessor is submission
`b89d1d72-d9a1-4ef9-9ed7-c6ef0eec1ed3`, immutable commit
`0e1d63e9afd0fdb04b311b1ae72a2618271b371b`, directory
`submissions/p256-ecdlp/rho-sigma2`, whose algorithm, arguments, arithmetic
sources and measurements `predecessor-proof.md` reproduces with an explicit
historical banner. No prior platform decision is an accepted premise, so the
reviewer assesses them anew.

**Precedence:** the claim is this document and `submission.json`. The
predecessor's Sections 1.1 and 1.3 define the same algorithm except for A, L
and the rare-work threshold replaced below; Sections 2.0--2.4 and 2.8 provide
the unchanged probability-space, escape, replay, abandonment and merge
arguments, re-evaluated here at A=2^24; Sections 3--4 provide the forced-path
pricing, resource caps and storage layout; the batched step and its
measurement are replaced below. This document replaces L, A
and every rate, threshold and derived count that depends on them, the reserve
derivation of 2.5, the success evaluation of 2.6, the common-path charge of 3,
and every total. The predecessor's 127.370918 score, 0.5303 collision
parameter, 64-fold reserve, 0.395898 success figure and 124-fold sensitivity
statement are historical, not current claims; its `evidence/accounting.json`
is reproduced as `evidence/predecessor-accounting.json`.

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

    tau_target = 494321/1000000,
    T_min = ceil(sqrt(2*N*tau_target)),
    L = ceil(T_min/K), T=K*L,
    S = K*(L+TRAIL_MAX).

Every slot takes L main-phase steps and finishes its current trail. Thus S
bounds all processed steps, including persistent steps and finishing tails.
Sampling caps, maximum distinguished records, abandonment cap, probe cap,
restart costs, bounded replay, verification and output are unchanged
formulas, evaluated at this S. The point-order certificate proving
ord_r(2)>2*TRAIL_MAX is recomputed in `accounting.py`.

The common path costs at most C_COMMON=784.858576 instructions per
processed step under H3, derived from the anchor schedule below. A fresh-source two-cycle entry costs at most w2=83,518
additional instructions. Longer entries cost at most wL=2,789,687;
false filter episodes, including their speculative exposure debt, are charged 256*(w2+wL); merged steps are charged wL.
These are the same inclusive forced-path/control allowances as before.
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
`accounting.py` asserts equals (3/(2A)+4/A^2)/(2A), below the inherited
allowance 3/A^2, 4/(1+8/(3A)) (nearly 4) times larger and not charged.
Epsilon's closure margin uses the same bound, each 3..64 window wasting 192
states: its third term is 2*192*(3/(4A^2)+2/A^3), not the predecessor's
1152/A^2. The split-concentration threshold mu=S/(2A), the
two-cycle reserve ceil((1+eta)*mu)*w2, the residual expectation E0 and the
reserves built on them are recomputed, not rescaled: the two-cycle charge is
0.002490 and E0 is 7.51340e-9 instructions per step, against 0.637201 and
4.87164e-4 at A=2^16.

The A-dependence of the compiled step is checked on the binaries (organizer
flags, ADDER_BITS=16 and 24, disassembled): every other function is
identical, and `batch_step` differs only in the AND and shift immediates of
the adder index and the two-limb distinguished-point test, branch offsets
and alignment no-ops. The 2^24 build runs two more no-ops per batch round,
before its first loop; `accounting.py` adds 2 to both rounds.
`select_adder1`'s empty asm keeps the index mask its own instruction
(at 16 bits it could fold into the address add). The two 64-bit
filters compare the whole low limb, so their 2^-63 false-match charge is
independent of A, and the table is data.

The larger table is charged, not assumed free: construction A*START =
4,281,292,120,129,536 instructions = 2^41.961183 units, storage A*2*64 = 2^31
bytes inside the peak-memory formula, and each 128-byte record's read-back
(16 words) lies inside START's 2^14 constant term. Preprocessing is
2^49.966808 units, below the declared 2^50 bound and included in W.

### The fused batched step

`batch_step` does the predecessor's field operations in the same order,
each subtraction inside the product or square that uses it, values
passed in registers (`measurement/p256_fused.h`): `k_submul`
d = x(sR)-x(Wc), PRE[w+1] = d*PRE[w]; `k_mul2` iw = inv*PRE[w], inv = inv*d;
`k_subm` lambda = (y(sR)-y(Wc))*iw; `k_xy` x3 = lambda^2-x(Wc)-x(sR),
y3 = (x(Wc)-x3)*lambda-y(Wc), with sR = sigma(R_i). The product is k_mul's
instructions (its two final constants in x7, x8), the square k_sqr's and
each subtraction k_subb_bl's, except that p's top limb takes one `and`
and k_mul2's two products end in a carry-only finish: for factors below
2^256 the Montgomery sum T is below 2^256+p, so T-c*p (c its fifth limb)
is congruent and below 2^256. iw and inv feed only products. Every other
product has a factor below p, so T<2p and k_mul's finish is canonical:
D, PRE, lambda, x3 and y3 equal the previous kernel's (two products swap
operands; k_mul's a*b/2^256 mod p is symmetric). The kernels have no branch. The distinguished-point test is one compare of
(key>>24)|(x.v[1] & mask) with zero, inside the measured round.
`p256_field.h`, `p256_sub.h` and the rare paths are unchanged and compile
instruction-identically, so `measurement/forced/summary.json` (measured on
the previous `step.c`, sha256 listed) still prices them. On seed 1, 8
rounds of 2^16 walks, `walk_ref.py check` finds 0 mismatches and the dump
equals the previous kernel's byte for byte.

### Common path charged by the anchor schedule

The predecessor charged every processed step the larger of its measured
plain and anchor batch rounds, plus 3, rounded up to 847, with the same
1/64 wide-counter allowance. That charges an anchor round's extra work at
every step.
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

* plain round, Ir(2)-Ir(1) = 51,420,348 on both sets, 917,495 below the
  previous kernel's: 14 fewer per walk (k_submul 1, k_mul2 8, k_subm 1,
  k_xy 4), less 9 as k_submul and k_mul2 run 2^16-1 times;
* anchor round: the previous kernel's Ir(64)-Ir(63) is 53,140,932 and
  53,305,368. As the walk state matches it after every round, branches
  agree and a round runs 917,495 fewer, so the larger gives 52,387,873
  (derived), used without removing the one and three two-cycle escapes it
  contains, which are also metered as rare work and so are charged twice.

With the 2 no-ops of the 2^24 build, PLAIN = 51420350/65536 = 784.612274...
and ANCHOR = 52387875/65536 = 799.375534... per walk step. On seed 1,
Ir(8)-Ir(4) = 4*51,420,348 + 82,194 (derived), one two-cycle escape over the window
`reproduce.sh` uses. Every conditional branch the common path executes per
walk is listed in `seed1.json`; each taken edge is a rare event, charged
separately. The orientation enters only through the table index, so a round
costs the same for all coin values. The next subsection's 1/64 counter
allowance is charged at every step.

A slot takes L main-phase steps and finishes its current trail, so its
batch runs at most R = L+TRAIL_MAX rounds and S = K*R. Of the rounds
t=1..R, at most ceil(R/64) are anchor rounds. This is a schedule count,
not an expectation, so every run's common path costs at most

    COMMON = S*(PLAIN + 1/64) + K*ceil(R/64)*(ANCHOR - PLAIN)

instructions, and C_COMMON = COMMON/S = 784.858576 (rounded up; the exact
rational is in `accounting.py`). Replays are unbatched and priced
separately; their per-trail anchors do not enter this count. The rare
paths and their forced prices are unchanged; the reserve built on them is
re-derived below.

### Production counter widths and control allowance

The historical C file is a finite measurement harness. Its signed `int`
loop variables, 64-bit batch round argument and unsaturated diagnostic
`ring_n` are not a production implementation. In the specified algorithm,
use unsigned 128-bit batch-round and replay-length counters; use unsigned
192-bit instruction meters and aggregate work counters. The per-trail
period counter `tper` is a single 64-bit word: `measurement/step.c`
increments it once per 64 processed steps and compares it with
CAP_PERIODS=2^60, so the cap TRAIL_MAX=64*(2^60+1) steps is enforced
through that 64-bit word, which holds at most CAP_PERIODS+1=2^60+1 (61
bits), and the step count is never materialised in a counter. The per-slot
anchor total is not a runtime counter: it is the analytic 128-bit quantity
ANCHOR_ROUNDS=ceil(R/64)~2^89.5 used by the W bound above, computed from
the round schedule. The rare-path diagnostic escape ring is omitted from
the attack (it is not used by the cycle filters or escape decisions).

Charge 1,024 extra scalar instructions per batch of 2^16 point updates,
or 1/64 instruction per processed step, for widened batch-round increment,
comparison to L and finishing control, including loads, stores and branch
operations; a two-limb increment and comparison need far fewer. Wide
counters touched only at rare events fit the existing 1,024-instruction
control allowance per rare call; a three-limb increment, load/store and
threshold comparison fit well below it. This change alters no field
kernel, adder choice, coordinate test or per-point ordinary loop, and
replay bookkeeping stays inside its conservative unbatched price. The
common price above thus includes the 1/64 at every step, with all wide
arithmetic charged.

The distinguished-record, trail/restart and table-index widths follow
from their maxima. The record counter reaches DP_MAX=2*floor(S/2^60)+2*K,
below 2^68.492 records; the trail counter reaches DP_MAX+2*K trails started
(the K initial starts plus at most DP_MAX+K restarts), the same magnitude;
the table index spans C=2*DP_MAX, below 2^69.492 slots. Each is a two-word
(128-bit) integer; its control is charged inside the table
operation the ledger already names. The ledger charges "DP_MAX table
operations at 2^21+2^11+16 instructions each, including bounded probing,
loads/stores, key comparisons, addressing and control, the index reduction
from the canonical key, the record-counter update and its overflow test",
`accounting.py` implements it as
`led.add_instructions("distinguished_points", DP_MAX * DP_COST)` with
DP_COST=PROBE_CAP*PROBE_COST+DP_INDEX_COST+DP_COUNTER_COST+DP_OVERFLOW_COST
=2^10*2^11+2^11+8+8=2,099,216=2^21+2^11+16; the ledger also charges
"(DP_MAX+K) restart constructions", priced by `accounting.py` at
START, 255,184,896 instructions. Each probe's own two-limb increment and
comparison costs at most 8 instructions against the probe allowance of
2^11, a factor of 256. The once-per-operation work is added to the charge
rather than assumed to fit in the probe allowance: the index is reduced
once per table operation, by a long division of at most 2^11 instructions;
the record counter is incremented and written back as two limbs in at most
8 instructions; and the overflow test against DP_MAX costs at most 8 more.
The worst case L*2^11+2^11+16 therefore equals DP_COST exactly, with no
once-per-operation work outside the charge. The K initial starts and
the DP_MAX+K restarts are charged K*START and (DP_MAX+K)*START, so the
trail counter adds at most 8*(DP_MAX+2*K) against START*(DP_MAX+2*K), a
factor of 31,898,112 for each eight-instruction counter update. No counter is
incremented per processed step: `tper` advances once per 64 steps and the
batch-round counter once per round, both inside the 1,024-instruction
allowance per round stated above. With the itemised charge the exact ledger
total is below 2^127.1422728697, 1.30e-7 bits below the declared
2^127.142273 claim.

## Concentration for the dominant event

This step uses the existing **conditional** fresh-source premise, not
merely a measured average rate. In the deferred-exposure filtration of
predecessor-proof.md Section 2.0 and Lemma E1(b), before exposing the key
and orientation of a globally fresh source, the probability that its
successor closes the two-cycle is p2=1/(2A). Stop the analytical exposure
process at the first violation of the predecessor's good-event
conditions, padding with zeros after that stopping time and after normal
termination. There are at most S trials, and its indicator I_t has
conditional mean at most p2 at each trial. Already exposed sources are
not charged here; their additional costs belong to the inherited
residual categories, including merged steps and within-trail coincidences.

Crucially, we do **not** condition the probability law on the final event
G. The stopped adapted process is defined on the original probability
space. On G it counts every fresh-source two-cycle entry; outside G the
separately charged DELTA_CUT applies, which avoids assuming that future
conditioning preserves the conditional event bound.

For lambda>0, conditional expectation gives

    E[exp(lambda*I_t) | F_(t-1)]
       <= 1+p2*(exp(lambda)-1)
       <= exp(p2*(exp(lambda)-1)).

Iterating this inequality and applying Markov to the exponential yields,
with mu=S*p2 and eta=2^-16,

    Pr[sum I_t > (1+eta)*mu]
       <= exp(-mu*((1+eta)*ln(1+eta)-eta))
       <= exp(-mu*eta^2/(2+eta)).

No independence of the indicators is required beyond the stated adapted
conditional bound; the final inequality is the standard scalar bound
(1+u)ln(1+u)-u >= u^2/(2+u), u>=0. Use the integer limit
n2=ceil((1+eta)*mu), with reserve B2=n2*w2. The accounting checks the
exponent exactly as a rational and verifies it exceeds 100, so the tail
is less than exp(-100)<2^-100, charged as delta2=2^-100; the computed
exponent is vastly larger, and the claim needs no underflowing
floating-point estimate.

This uses the disclosed M1/L-ori premise more strongly. Bounded
reduced-order measurements alone do not establish it at production
lengths: if it is only a marginal-rate observation rather than the
conditional law assumed, the concentration conclusion does not follow.
This distinction is a material limitation.

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
returns to visited main-trail classes merges with a different trail
(already charged), closes a short cycle (covered by this debit and the
128-step bound), or violates G (charged as failure). Continuing past the
finite speculative path's end exposes a new source unless one of those
cases occurs. Thus the split never silently treats a known look-ahead key
as fresh randomness. This is an analytical assignment; the runtime still
meters the same actual calls.

For the residual success argument, stop the analytical residual process
at the first G violation, just as for the two-cycle process. Let Eres=E0
and kM=2^18. Set

    B_res = kM*Eres*S,
    B_rare = B2+B_res.

The factor kM is a free analysis parameter, not a measured constant: the
reserve charges kM*Eres at every processed step, so raising kM multiplies the
reserve and shrinks the abort allowance 1/kM. E0's three parts scale
differently: only the longer-entry part depends on A, as
3/(4A^2)+2/A^3 = (3/(4A^2))(1+8/(3A)), A^-2 up to the factor (1+8/(3A)), so
it falls 412333637632/6291457-fold (more than 65,538.65-fold) from A=2^16 to 2^24,
while the false-filter part is 7.974745864e-11 and the merged part 4.253047933e-13
at every A. With tau and S at their production values, E0 is 4.871633669e-4
at A=2^16, 1.902989172e-6 at A=2^20 and 7.513393318e-9 at A=2^24;
`accounting.py` prints the first two and every component. So E0 falls
more than 64,839-fold over A=2^16..2^24, less than 2^16-fold, since only the
longer-entry part falls. The declared pair keeps a margin of 1.0e-5: at
tau=494321/1000000 the factor 2^18 keeps 1.0033e-5 above 0.39 and is the
lowest score among the powers of two that keep at least 1.0e-5; 2^17 scores
1.8e-6 lower but keeps only 6.2e-6 above 0.39, and 2^19 keeps more margin at
3.6e-6 more work. Its reserve is 0.001970 instructions per step and its abort
allowance is 2^-18, 309 times smaller than the inherited 1/849. Its threshold
is public algorithm data.

All within-trail coincidences and excessive-length formal closures are
excluded by G and paid through DELTA_CUT; we do not bound their realized
costs by the average per-step allowance, nor use the predecessor's
within-trail fixed-point argument. On G the stopped residual process is
the actual residual cost, and its expectation is at most Eres*S by the
retained rate bounds plus the new speculative-exposure debit. Markov
gives a 1/kM = 2^-18 tail.

    Pr[rare work > B_rare and G] <= 2^-100 + 2^-18.

No independence between residual and two-cycle work is asserted. Meter
the actual rare operations at their charged worst-case prices, just as
in the predecessor. Abort before executing a priced operation that would
exceed the total threshold B_rare, so computation is bounded even when
the probabilistic premises fail. The event decomposition adds no runtime
counter, global-freshness oracle or uncharged lookup. The inherited
control allowance pays the same meter operations at a different constant
threshold. A constant-sized setup calculation can be absorbed by the
input-validation allowance; all profile-dependent parameters may also be
compiled as public algorithm constants.

## Complete cost and success ledger

The ledger charges (COMMON+B_rare)/1000 for all walk steps and metered
rare work, then **adds**, rather than amortizes away:

* A adder constructions and K initial trail starts, each with the complete
  capped scalar-sampling and two-scalar-multiplication price;
* up to DP_MAX+K restart constructions, DP_MAX=2*floor(S/2^60)+2*K;
* DP_MAX table operations at 2^21+2^11+16 instructions each, including
  bounded probing, loads/stores, key comparisons, addressing and control,
  and the once-per-operation index reduction, record-counter update and
  overflow test;
* two bounded coefficient replays at the original generous unbatched
  operation count, followed by solving, verification and serialization;
* input validation and setup allowances;
* the Lemma F closure series, summed for s < 4096 with the tail bounded by
  the consecutive-term ratio (e^2/8)(1+3/A)^2 < 0.924, so the truncation is an
  upper bound, equalling the fully enumerated sum at A=2^16 and A=2^20.

All processors' work is summed. No future-target amortization or free
advice is used. The storage formula is unchanged:

    2*DP_MAX*128 + K*512 + A*2*64 + 2^24 bytes.

It includes distinguished records, live walk states and batch buffers,
both signs of adders, and a fixed reserve. At A=2^24 it evaluates to
2^76.491761 bytes, dominated by the 2*DP_MAX*128 record table. `accounting.py`
evaluates these formulas and asserts that the total, rounded up at the sixth
decimal, equals the declared score 127142273/1000000. That is a rational
submission bound, not a timing result.

The collision exponent penalty re-derived at the bound, with the predecessor's 2.6 margins:

    epsilon = 4/A + 2*(3/A) + 2*192*(3/(4A^2)+2/A^3)
            = 10/A + 288/A^2 + 768/A^3 = 5.960475e-7.

The conditional success lower bound is

    1-exp(-tau_exact*(1-epsilon))
      - 2^-18 - 2^-100 - 2^-20
      - DELTA_ABANDON - DELTA_CUT - DELTA_MERGE
      - DELTA_CAP - DELTA_PROBE.

Each term is recomputed from the unchanged predecessor formula with the
new S and T; no historical numerical failure allowance is silently
reused. The parameter screen takes the smallest tau on its 1e-6 grid at
which the conditional bound keeps a 1.0e-5 margin above 0.39. At
tau=494321/1000000 the model term is 0.3900148737819886 and the exact bound is
0.39001003329294714, a margin of 1.0033e-5 above 0.39; the named failure terms at this S are
DELTA_ABANDON=7.2117e-8, DELTA_CUT=5.8543e-17, DELTA_MERGE=1.2520e-14,
DELTA_CAP=2.4393e-18 and DELTA_PROBE=2.4924e-44. The bundled ledger asserts
p_claim >= 0.390005, the 1.0e-5 margin and the claimed cost. The margin is a
conditional-model margin, not a confidence interval for empirical trials.
The legacy 124-fold residual-expectation sensitivity no longer applies:
this candidate trades that large reserve for a stronger conditional-tail
argument and for the residual meter the A=2^24 entry rates make affordable.

## Evidence scope, execution and limitations

The inherited field assembly, enumerated forced rare paths, replica
results and reduced-order full procedure records are included with
attribution, supporting H3 and the model at their stated scopes; the
plain round is a new measurement of the fused step by the same method,
the anchor round derived from the previous kernel's. None establishes a production success probability. The
predecessor's measurement harness uses pooled restart points; the
specified attack requires fresh coins and pays their full generation
cost, a distinction retained without treating pooled runs as a
measurement of cold production restarts.

One declared combined experiment uses the verifier's bounded credential-free
runner and returns two separately labelled results. `accounting.py`
recomputes the new numerical ledger, concentration threshold, failure
terms and point-order certificate. `experiment.py` is the inherited
complete reduced-order procedure, with the declared collision parameter
0.4943, which is production tau=0.494321 on the available 0.0001 input
grid, 2.1e-5 below it; its finite-run success is diagnostic only. It does
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
