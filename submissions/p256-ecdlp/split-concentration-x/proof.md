# P-256 sigma rho with a concentrated two-cycle reserve

## Current claim and provenance

For a uniform scalar x in 0..r-1, given Q=[x]P on the protected NIST P-256
profile, the following classical single-target algorithm has success probability
at least 0.39 under the disclosed H1, M1/L-ori and fresh-class exposure premises.
Its total charged computation is at most **2^127.241865** target-group-operation
equivalents. The normalization is C=1000 scalar AArch64 instructions per unit;
an actual P-256 addition is not assumed to cost one unit. Peak memory is below
2^77 bytes, all preprocessing is below 2^50 units and is included in the total,
and nonuniform advice is zero bytes. The instruction accounting retains H3,
the predecessor's implementation-measurement assumption.

The exact declared score is 25448373/200000 (127.241865), the ledger's total
rounded up at the sixth decimal. This is a conditional analytic bound under
an exploratory acceptance policy. It is neither an unconditional complexity
theorem nor an executed full-size attack.

The material predecessor is submission
`b89d1d72-d9a1-4ef9-9ed7-c6ef0eec1ed3`, immutable commit
`0e1d63e9afd0fdb04b311b1ae72a2618271b371b`, directory
`submissions/p256-ecdlp/rho-sigma2`. Its full algorithm, formal arguments,
arithmetic sources and historical measurements accompany this candidate.
`predecessor-proof.md` reproduces its argument with an explicit historical
banner. No prior platform decision is an accepted mathematical premise.
The reviewer must assess the reproduced premises and arguments anew.

**Precedence:** the current claim is this document and `submission.json`.
The predecessor's Sections 1.1 and 1.3 define the same algorithm except for
L and the rare-work threshold replaced below. Its Sections 2.0--2.4 and 2.8
provide the unchanged probability-space, escape, replay, abandonment and
merge arguments. Its Sections 3--4 provide the unchanged kernel pricing,
resource caps and storage layout. This document replaces the numerical
parameter L, the reserve derivation of Section 2.5, the numerical success
evaluation of Section 2.6, the common-path charge of Section 3, and every
resulting total. The predecessor's
127.370918 score, 0.5303 collision parameter, 64-fold reserve, 0.395898
success figure and 124-fold sensitivity statement are historical and are
not current claims. Original `evidence/accounting.json` is reproduced here
as `evidence/predecessor-accounting.json`.

## Algorithm, fixed costs and changed parameters

Use the protected p, a=-3, b, r and generator P without alteration. Handle
Q=identity by returning zero; otherwise validate the curve point. Let a
class be {W,-W}, key(W) the integer Montgomery representation of its affine
x-coordinate, and canon(W) the representative with even Montgomery y.
Let sigma(W) be the sign taking canon(W) to W. The walk is

    W -> W + sigma(W) R[key(W) mod A].

Here A=2^16 adders, each R_i=[c_i]P+[d_i]Q with fresh uniform scalar coins.
The table contains both signs. There are K=2^32 walk slots, arranged into
batches of 2^16 for simultaneous inversion. Distinguished points have
bits 16..75 of the key zero, so theta=2^-60 under the stated class model.
Each trail starts from fresh independent a,b with W=[a]P+[b]Q. It ends at
its first distinguished point, records it once with its coins, then starts
a fresh trail. Coefficients are reconstructed by replay only at a match.

The exact short-cycle rules, test ordering and replay algorithm are those
fully specified in predecessor-proof.md Sections 1.2--1.3. In particular,
the distinguished-point test precedes cycle tests at every generated
point, including look-ahead and escape points. A detected two-cycle or
cycle of length 3..64 escapes by doubling the canonical representative
of its least-key member. The two 64-bit filters only trigger exact checks.
Every trail is limited to TRAIL_MAX=64*2^60+64 processed steps; the first
key match is resolved by two replays, each limited to 65 times its original
live length. Degenerate matches and exceeded caps fail the run. Each
derived scalar is verified by [x]P=Q before output.

Set N=(r-1)/2+1 and

    tau_target = 496471/1000000,
    T_min = ceil(sqrt(2*N*tau_target)),
    L = ceil(T_min/K), T=K*L,
    S = K*(L+TRAIL_MAX).

Every slot takes L main-phase steps and finishes its current trail. Thus S
bounds all processed steps, including persistent steps and finishing tails.
Sampling caps, maximum distinguished records, abandonment cap, probe cap,
restart costs, bounded replay, verification and output remain unchanged
as formulas, evaluated at this S. The point-order certificate proving
ord_r(2)>2*TRAIL_MAX is recomputed in `accounting.py`.

The common path costs at most C_COMMON=836.842722 instructions per
processed step under H3, derived from the anchor schedule below. A fresh-source two-cycle entry costs at most w2=83,518
additional instructions. Longer entries cost at most wL=2,789,687;
false filter episodes, including their speculative exposure debt, are charged 256*(w2+wL); merged steps are charged wL.
These are the same inclusive forced-path/control allowances as before.
The runtime still uses **one total rare-work meter**. It does not need to
recognize whether a point is globally fresh, or to classify observed
events as fresh versus merged. That partition is solely for analysis.

### Common path charged by the anchor schedule

The predecessor charged every processed step the larger of its measured
plain and anchor batch rounds, plus 3, rounded up to 847; the 1/64 below
was then added. That charges an anchor round's extra work at every step.
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

* plain round, Ir(2)-Ir(1) = 54,631,599 on both sets, so
  PLAIN = 54631599/65536 = 833.612045... per walk step;
* anchor round, Ir(64)-Ir(63) = 55,369,153 and 55,533,591. The larger is
  used as measured, without removing the one and three two-cycle escapes it
  contains, which are also metered as rare work and so are charged twice:
  ANCHOR = 55533591/65536 = 847.375350... per walk step.

The predecessor's anchor figure with those escapes removed at their average
cost, 55286976/65536 = 843.612304..., is smaller and is not used. The other
edge of the first distinguished-point test (key>>16 = 0) evaluates the
second test and adds at most 3 instructions. It is charged at every step,
whether or not a measured step took it, as is the 1/64 counter allowance
of the next subsection. Section 3 of predecessor-proof.md (H3)
gives the uniform cost of a round over all coin values.

A slot takes L main-phase steps and finishes its current trail, so its
batch runs at most R = L+TRAIL_MAX rounds and S = K*R. Of the rounds
t=1..R, at most ceil(R/64) are anchor rounds. This is a count from the
schedule, not an expectation, so on every run the common path costs at most

    COMMON = S*(PLAIN + 3 + 1/64) + K*ceil(R/64)*(ANCHOR - PLAIN)

instructions, and C_COMMON = COMMON/S = 836.842722 (rounded up; the exact
rational is in `accounting.py`). Replays are unbatched and priced
separately, so their per-trail anchors do not enter this count. The rare
paths, the reserve and every other ledger term are unchanged.

### Production counter widths and control allowance

The historical C file is a finite measurement harness. Its signed `int`
loop variables, 64-bit batch round argument and unsaturated diagnostic
`ring_n` are not a production implementation. In the specified algorithm,
use unsigned 128-bit batch-round and replay-length counters; use unsigned
192-bit instruction meters and aggregate work counters. The per-slot
anchor-period count already fits in 64 bits because the cap is 2^60.
The rare-path diagnostic escape ring is omitted from the attack (it is
not used by the cycle filters or escape decisions). Historical source
and measurement records remain unchanged for reproducibility.

Charge 1,024 extra scalar instructions per batch of 2^16 point updates,
or 1/64 instruction per processed step, for widened batch-round increment,
comparison to L and finishing control, including loads, stores and branch
operations. A two-limb increment and comparison require far fewer than
this allowance. Wide counters touched only at rare events fit the existing
1,024-instruction control allowance per rare call; a three-limb increment,
load/store and threshold comparison fit well below it. This counter-width
change does not alter the field kernel, adder choice, coordinate tests,
or the per-point ordinary loop. Replay bookkeeping is likewise inside
its very conservative unbatched price. Therefore the common price above
includes the 1/64 at every step, with all wide arithmetic charged.

## Concentration for the dominant event

This step uses the existing **conditional** fresh-source premise, not
merely a measured average rate. In the deferred-exposure filtration of
predecessor-proof.md Section 2.0 and Lemma E1(b), before exposing the key
and orientation of a globally fresh source, the probability that its
successor closes the two-cycle is p2=1/(2A). Stop the analytical exposure
process at the first violation of the predecessor's good-event conditions.
Pad it with zeros after that stopping time and after normal termination.
There are at most S trials, and its indicator I_t has conditional mean
at most p2 at each trial. Sources already exposed are not charged to
this count; their additional costs belong to the inherited residual
categories, including merged steps and within-trail coincidences.

Crucially, we do **not** condition the probability law on the final event
G occurring. The stopped adapted process is defined on the original
probability space. On G it counts every fresh-source two-cycle entry;
outside G the separately charged DELTA_CUT applies. This avoids assuming
that future conditioning preserves the conditional event bound.

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
conditional bound. The final inequality follows by the standard scalar
calculus bound (1+u)ln(1+u)-u >= u^2/(2+u), u>=0.
Use the integer limit n2=ceil((1+eta)*mu), with reserve B2=n2*w2.
The accounting checks the exponent exactly as a rational and verifies it
exceeds 100. Thus the tail is less than exp(-100)<2^-100, and we charge
delta2=2^-100. In fact the computed exponent is vastly larger. The claim
does not require estimating an underflowing floating-point probability.

This is stronger use of the already disclosed M1/L-ori premise. Bounded
reduced-order measurements alone do not establish that premise at
production lengths. If it is only a marginal-rate observation rather
than the conditional law explicitly assumed, the concentration conclusion
does not follow. This distinction is a material limitation.

## Residual reserve and a single enforceable meter

Write c2=B2/S. The expected residual event cost per processed step before
within-trail exceptional work is

    E0 = (3/A^2)*wL
       + 2^-63*256*(w2+wL)
       + (tau_exact*TRAIL_MAX/S)*wL,
    tau_exact = T^2/(2N).

The longer-entry and merged terms retain their predecessor prices. The false-filter
term has a new generous debit for speculative look-ahead. A false first-filter
check exposes at most one future class; a false anchor check exposes at most 65.
Those classes cannot subsequently be treated as conditionally unexposed for the
concentration proof. Charge the originating false-filter episode for its direct
work, at most 66 later exposed-source steps, and at most 128 persistent steps
through an ensuing short cycle before escape. The debit 256*(w2+wL) exceeds
that sum. Overlapping exposures may be charged repeatedly, which only overcounts.
An unclosed speculative path that returns to previously visited main-trail
classes either merges with a different trail (already charged), closes a
short cycle (covered by this debit and the 128-step bound), or violates G
(charged as failure). Continuing beyond the end of the finite speculative
path exposes a new source unless one of those cases occurs. Thus the split
never silently treats a known look-ahead key as fresh randomness. This is
an analytical assignment; the runtime still meters the same actual calls.

For the residual success argument, stop the analytical residual process
at the first G violation, just as for the two-cycle process. Let Eres=E0
and kM=849. Set

    B_res = kM*Eres*S,
    B_rare = B2+B_res.

All within-trail coincidences and excessive-length formal closures are
already excluded by G and paid through DELTA_CUT. We do not bound their
individual realized costs by the average per-step meter allowance, and
we do not use the predecessor's within-trail fixed-point cost argument.
On G the stopped residual process is the actual residual cost, and its
expected cost is at most Eres*S by the retained rate bounds plus the new
speculative-exposure debit. Markov gives a 1/849 tail. Therefore

    Pr[rare work > B_rare and G] <= 2^-100 + 1/849.

No independence between residual and two-cycle work is asserted. Meter
the actual rare operations by their charged worst-case prices, just as
in the predecessor. Abort before executing a priced operation that would
exceed the total threshold B_rare. Therefore computation is bounded even
when the probabilistic premises fail. The event decomposition does not
add any runtime counter, global-freshness oracle or uncharged lookup.
The inherited control allowance pays the same meter operations with a
different constant threshold. A constant-sized setup calculation can
be absorbed by the input-validation allowance; all profile-dependent
parameters may also be compiled as public algorithm constants.

## Complete cost and success ledger

The ledger charges (COMMON+B_rare)/1000 for all walk steps and metered
rare work, then **adds**, rather than amortizes away:

* A adder constructions and K initial trail starts, each with the complete
  capped scalar-sampling and two-scalar-multiplication price;
* up to DP_MAX+K restart constructions, DP_MAX=2*floor(S/2^60)+2*K;
* DP_MAX table operations at 2^21 instructions each, including bounded
  probing, loads/stores, key comparisons, addressing and control;
* two bounded coefficient replays at the original generous unbatched
  operation count, followed by solving, verification and serialization;
* input validation and setup allowances.

All processors' work is summed. No future-target amortization or free
advice is used. The storage formula is unchanged:

    2*DP_MAX*128 + K*512 + A*2*64 + 2^24 bytes.

It includes distinguished records, live walk states and batch buffers,
both signs of adders, and a fixed reserve. `accounting.py` evaluates these
formulas and asserts that the total, rounded up at the sixth decimal, equals
the declared score 25448373/200000. That is a rational submission bound, not a
timing result.

Keep the original collision exponent penalty

    epsilon = 10/A + 1152/A^2.

The conditional success lower bound is

    1-exp(-tau_exact*(1-epsilon))
      - 1/849 - 2^-100 - 2^-20
      - DELTA_ABANDON - DELTA_CUT - DELTA_MERGE
      - DELTA_CAP - DELTA_PROBE.

Each named term is recomputed from the unchanged predecessor formula with
the new S and T. No historical numerical failure allowance is silently
reused. The analytic parameter screen selects tau with roughly 0.0001
probability margin above 0.39; the bundled ledger asserts both this
threshold and the claimed cost. The margin is conditional-model margin, not
a confidence interval for empirical trials. The legacy 124-fold residual
expectation sensitivity no longer applies. This candidate deliberately
trades that large legacy reserve for a stronger conditional-tail argument.

## Evidence scope, execution and limitations

The inherited field assembly, measured production-curve common paths,
enumerated forced rare paths, replica results and reduced-order full
procedure records are included with attribution. They support H3 and
the model at their stated scopes. They are not new measurements by this
candidate and do not establish a production success probability.
The predecessor's measurement harness uses pooled restart points; the
specified attack requires fresh coins and pays their full generation
cost. That distinction is retained without treating pooled runs as a
measurement of cold production restarts.

One declared combined experiment uses the verifier's bounded credential-free
runner and returns two separately labelled results. `accounting.py` recomputes the new numerical ledger, concentration
threshold, failure terms and point-order certificate. `experiment.py`
is the inherited complete reduced-order procedure, with the declared
collision parameter 0.4965 (the available input grid has spacing 0.0001).
The latter is close to but not identical to production tau=0.496471;
its finite-run success is diagnostic only. It does not implement the
production instruction meter or execute a full-size attack. Its actual
successes, failures and replay checks must be reported as observed.

The remaining material uncertainties are production-scale collision
mixing, conditional fresh-source orientation, the validity of the
fresh-window domination model, analytic handling of long trails and
rare returns, and universal validity of the measured kernel bounds.
The new exponential-tail argument does not resolve those uncertainties;
it is conditional on them. Memory remains exponential and impractical.
No claim of practical compromise, unconditional improvement in P-256
security, or human certification is made.
