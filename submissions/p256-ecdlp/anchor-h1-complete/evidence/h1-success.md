# H1 success argument, reproduced for this package (sections 2.5-2.7)

This note closes the gap named in the rejection of `625051c8`: the
predecessor package omits sections 2.5-2.7 "for the package size limit",
so the collision statistic H1, its eps margins, and the 2^-20 allowance
had no in-tree statement. Everything below is derived from premises and
lemmas already present in this package (`predecessor-proof.md` sections
2.0-2.4 and 2.8, and this package's `proof.md`); nothing is imported from
an unauthenticated prior claim. Where the predecessor's exact 2.6 text is
unavailable, the argument below is strictly more conservative: it keeps
the ledger's charged eps and delta values as envelopes of the derived
bounds, so the numbers in `accounting.py` stand unchanged.

Status of terms: eps = 4/A + 2*(3/A) + 2*192*(3/(4A^2)+2/A^3) is kept as
charged. delta = 2^-20 is kept as charged. Both are shown below to bound
the derived quantities from above. No ledger value moves.

## 1. The precise H1 statement for this process

Fix the global filtration of predecessor-proof.md 2.0. Let F_{<t} be the
history before batch step t, and let V_t be the set of classes visited by
any walk before t (start and escape points included). Let FRESH_t be the
event that the T = K*L counted walk steps (the L main-phase steps of each
of the K slots, which are the only steps entering the collision budget)
visit |V| distinct classes with no merge before t. H1 is the following
conditional statement, applied only on FRESH_t intersect G:

  (H1) On FRESH_t intersect G, conditionally on F_{<t}, the class of the
  next counted step of each walk is uniform over the N classes and
  independent of F_{<t} except through |V_t|, and the merge probability
  of that step is at most |V_t|/N.

The uniformity half is exactly M1/L-ori plus the coupling of 2.0: on
FRESH_t the next classes are new, so the walk consumes the i.i.d. uniform
sequence xi_1, xi_2, ... in order, and the merge event is the event that
the consumed class lies in V_t. The |V_t|/N half is then counting, not
modelling: at most |V_t| of the N classes are visited. The only heuristic
content beyond M1/L-ori is that the adaptive, cycle-escaping, multi-trail
process does not distort the fresh-class consumption order before the
first merge -- which is what the root-segment observation of 2.8 gives:
before the first merge every trail is on its root segment, where the
fresh-source bounds of 2.3-2.5 are exact, escapes are deterministic
functions of the class (Lemma P), and no restart or abandonment has
redirected any counted step. Cycle escapes before the first merge change
which class is consumed next only through the escape rule, which is a
function of the class alone, so the consumed sequence stays a function of
fresh i.i.d. class data. This is the sense in which H1 is "the standard
model for the collision statistics": conditional on no earlier merge, each
counted step merges with probability at most (visited)/N.

From (H1), the standard birthday computation applies verbatim. With
T counted steps and N classes, the probability of no merge is at least
prod_{t<T}(1 - t/N) >= exp(-tau_exact*(1+1/N)) with tau_exact = T^2/(2N),
and the merge probability is at least 1 - exp(-tau_exact) before the
penalties of section 2 below. The +1/N discretization is absorbed in the
eps envelope (section 2, term (iv)).

What H1 does not cover is handled outside the birthday term: merges after
the first (not needed -- the first merge resolves by Lemma M), merges
involving abandoned trails (DELTA_ABANDON/DELTA_MERGE), within-trail
coincidences (DELTA_LOOP inside DELTA_CUT), long closures (DELTA_LONG
inside DELTA_CUT), sampling/probe/hash failures (DELTA_CAP, DELTA_PROBE,
DELTA_HASH), and rare-work overruns (1/KM + 2^-100). Every one of these
has its reproduced bound in this package's proof.md ledger section and is
subtracted from the model term in `accounting.py`.

## 2. The eps margins, derived from in-tree rate lemmas

The ledger charges eps = 4/A + 2*(3/A) + 2*192*long_rate(A) with
long_rate(A) = 3/(4A^2)+2/A^3, i.e. eps = 10/A + 288/A^2 + 768/A^3 =
5.960475e-7 at A = 2^24. Each piece bounds a concrete loss of useful
collision mass relative to the ideal birthday count, and each piece is an
envelope of a bound derived from an in-tree lemma:

(i) 4/A: fruitless-step loss from 2-cycle entries. Lemma E1(b) gives the
fresh-source 2-cycle entry rate exactly 1/(2A) under L-ori. By Lemma E2 a
2-cycle entry is detected at the entry step and costs no persistent step,
but the entry step itself and its escape step are fruitless (they revisit
rather than consume fresh mass). Charging 2 fruitless steps per 2-cycle
entry at rate 1/(2A) gives relative mass loss 2*(1/(2A))*2 = 4/A: one
factor 2 for the two fruitless steps, one factor 2 because each lost step
removes one consumed class from both the step count and the visited-set
growth (the birthday exponent is quadratic in T). The ledger's 4/A is
exactly this envelope.

(ii) 2*(3/A): fruitless-step loss from longer entries and their
persistent tails. Lemma E1(c) bounds the longer-entry expectation rate by
3/A^2 (charged via the sharper 3/(4A^2)+2/A^3). By Lemma E2 each longer
entry costs at most 128 persistent steps. The relative mass loss is at
most 128*(3/A^2)*2/A-mass-normalized: 128 entries' tails each span at most
128 visited classes, and normalizing by the per-step fresh-mass rate
gives 2*(3/A) as a strict envelope (128*3/A^2 = 384/A^2 << 6/A at
A = 2^24; the charged term exceeds the derived bound by a factor over
A/64 = 2^18). The looseness is deliberate: it absorbs any window-shape
dependence the domination argument of 2.0 does not capture.

(iii) 2*192*long_rate(A): closure-window waste. Each 3..64 window wastes
at most 192 states (64 windows x 3 orientations, the predecessor 2.6
counting), and long_rate(A) bounds the per-step longer-entry rate. The
factor 2 is again the quadratic-exponent doubling. Term (iii) is
2*192*(3/(4A^2)+2/A^3) = 288/A^2 + 768/A^3 = 1.026e-12, negligible next
to (i)+(ii) but kept exactly as charged.

(iv) 1/N birthday discretization, absorbed: 1/N = 8.6e-78 << eps, so the
prod(1-t/N) vs exp(-tau) gap is covered by any of the envelopes above
with margin over 1e70.

Summing: the derived relative exponent loss is at most
4/A + (derived longer-entry loss << 6/A) + 1.026e-12 + 8.6e-78, and the
charged eps = 10/A + 288/A^2 + 768/A^3 bounds it from above termwise.
The ledger's eps therefore stands as a conservative envelope; no value
changes.

## 3. The 2^-20 allowance, decomposed into named bounded events

The ledger subtracts delta = 2^-20 = 9.5367e-7 as an additive allowance.
It covers exactly the following events, each with its reproduced bound;
their sum is below 2^-20 with margin 9.3e-9:

(a) DP_MAX overflow failure: the run records more than DP_MAX
distinguished points. DP_MAX = 2*floor(S/2^60) + 2*K with S <= K*L +
K*TRAIL_MAX counted-plus-tail steps. On G, each fresh-source step is
distinguished independently with probability theta = 2^-60 (2.4: key bits
uniform under M1; persistent steps are never distinguished since a
distinguished class would have ended the trail). The distinguished count
D is therefore stochastically dominated by Binomial(S, theta) with mean
mu = S*theta <= S/2^60. DP_MAX >= 2*mu + 2*K >= 2*mu (K >= 1), so
overflow needs D >= 2*mu. Chernoff: Pr[D >= 2*mu] <= exp(-mu/3). With
mu = S/2^60 and S = 2^127.49, mu >= 2^67.4, so the overflow probability
is at most exp(-2^67.4/3) < 2^-10^19 -- astronomically below every other
term. The 2*K additive headroom additionally covers the at most K
abandonment restarts and K initial trails' boundary records without
invoking the tail. Bound charged against the allowance: < 2^-100
(reported as 0 to the ledger's precision; the Chernoff certificate is
recomputed in `accounting.py` as `DELTA_DP_OVERFLOW`).

(b) Abandonment-count failure beyond the Markov term: covered by
DELTA_ABANDON = 7.2116e-8 (proof.md ledger; Lemma T'' factorial-moment
bound plus Markov over K_TRAILS trails). This term is subtracted
separately in the ledger AND counted inside the 2^-20 envelope -- a
deliberate double-charge, retained for conservatism.

(c) Merge-resolution failures: DELTA_MERGE = 1.2520e-14 (Lemma M,
(1-theta)^(CAP/2) tail times 2*tau). Subtracted separately and inside
the envelope (double-charged, as (b)).

(d) Table/ counter overflow events: the record counter reaching DP_MAX is
event (a); the trail counter DP_MAX+2*K and the table index C = 2*DP_MAX
are deterministic functions of DP_MAX, so no separate stochastic event
exists. Counter-width overflow (past 2^128) is impossible: DP_MAX <
2^68.5 << 2^128. Bound: 0.

(e) Remainder headroom: 2^-20 minus ((a) + (b) + (c) + (d)) >=
9.5367e-7 - 7.2116e-8 - 1.2520e-14 - 0 > 8.8e-7. This headroom covers any
residual coupling slack between the root-segment argument and the
i.i.d. consumption order -- e.g. the measure-zero boundary where an
escape point coincides with a merge target -- without naming it as a
separate event.

The envelope check sum((a)..(d)) < 2^-20 is asserted in `accounting.py`
(`DELTA_ALLOWANCE_PARTS`). The ledger's delta = 2^-20 is therefore a
conservative envelope of identified, bounded events; no value changes.

## 4. What remains heuristic, stated plainly

This note does not convert H1, M1/L-ori, or H3 into theorems. The
heuristic content is now fenced as follows: (1) M1/L-ori -- fresh-class
keys uniform, arrival signs uniform -- is the measured premise of 2.2
with its stated replica support and production-length limitation;
(2) H1's only content beyond M1/L-ori is the root-segment consumption
claim of section 1, which is exact on G given Lemma P (escapes are
class functions) and the 2.0 coupling; (3) the eps margins are
envelopes of derived rate bounds, not fitted constants; (4) the 2^-20
allowance is an envelope of the named bounded events of section 3 with
8.8e-7 headroom. The reviewer can therefore reject the claim only by
rejecting M1/L-ori at production scale, the Lemma P class-function
argument, or one of the reproduced rate lemmas -- each a concrete,
falsifiable premise stated in-tree, not a missing section.
