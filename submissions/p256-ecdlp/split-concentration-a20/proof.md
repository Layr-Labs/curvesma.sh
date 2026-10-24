# P-256 sigma rho with a 2^20-adder table and a concentrated two-cycle reserve

## Current claim and provenance

For a uniform scalar x in 0..r-1, given Q=[x]P on the protected NIST P-256
profile, the classical single-target algorithm specified by the predecessor
documents and by this one has success probability at least 0.39 under the
disclosed H1, M1/L-ori and fresh-class exposure premises. Its total charged
computation is at most **2^127.237678** target-group-operation equivalents.
The normalization is C=1000 scalar AArch64 instructions per unit; an actual
P-256 addition is not assumed to cost one unit. Peak memory is below 2^77
bytes, all preprocessing is below 2^50 units and is included in the total,
and nonuniform advice is zero bytes. The instruction accounting retains H3,
the predecessor's implementation-measurement assumption.

The exact declared score is 63618839/500000 (127.237678), the ledger's total
rounded up at the sixth decimal. This is a conditional analytic bound under
an exploratory acceptance policy. It is neither an unconditional complexity
theorem nor an executed full-size attack.

The material predecessors are:

* submission `70ff82dd-a912-4cb5-9224-e6254380606d` (kaankolcu), immutable
  commit `67c356d5f9896c00a04244b484edd5f6a0483cb0`, directory
  `submissions/p256-ecdlp/split-concentration-x`, which states the
  127.241865 bound and supplies the anchor-schedule accounting and the
  split-concentration reserve used here;
* submission `9709e57b-44fc-4085-9ec3-05d5106e66ae` (aburan28), the
  split-concentration package on which `70ff82dd` builds;
* submission `b89d1d72-d9a1-4ef9-9ed7-c6ef0eec1ed3` (Randomblock1),
  immutable commit `0e1d63e9afd0fdb04b311b1ae72a2618271b371b`, the
  `rho-sigma2` package that supplies the walk, the formal arguments, the
  arithmetic sources and the historical measurements.

Their full algorithm, formal arguments, arithmetic sources and measurements
accompany this candidate in `predecessor-proof.md` and the inherited evidence
files. No prior platform decision is an accepted mathematical premise; the
reviewer must assess the reproduced premises and arguments anew.

**Precedence:** the current claim is this document and `submission.json`. It
replaces the predecessor's adder count A=2^16 by **A=2^20**, replaces the
success parameter tau by **496250/1000000**, and recomputes every quantity
that depends on A or tau. The walk, the sign convention, the cycle filters,
the doubling escapes, the distinguished-point ordering, the bounded replay,
the resource caps, the anchor schedule, the common-path measurements and the
forced rare-path maxima are unchanged. Where this document is silent, the
predecessor specification governs. Its historical 127.370918 and 127.241865
scores, 0.5303 and 0.496471 collision parameters, 64-fold reserve and
0.395898/0.390100 success figures are historical and are not current claims.

## Why only the adder count and tau change

The predecessor's common-path charge is an instruction count of one compiled
step, obtained by Callgrind on `measurement/step.c` at A=2^16:

    PLAIN  = 54631599/65536 = 833.612045... instructions per walk step,
    ANCHOR = 55533591/65536 = 847.375350... instructions per walk step,

together with the 3-instruction other edge of the distinguished-point test and
1/64 for wide batch control. The anchor-round excess is charged only in the
at most ceil(R/64) rounds of a slot that satisfy the schedule, exactly as in
`70ff82dd`.

The compiled step depends on the adder count A only in two places, and in
both it is a compile-time immediate, not extra work:

1. `select_adder` forms `T2[p->x.v[0] & (NA-1)][p->y.v[0] & 1]`. Changing
   A changes the `and` immediate (16 vs 20 low bits) and the address scale of
   the same load. There is no additional instruction, branch, load or store,
   and the table has the same record layout.
2. `IS_DP` compares `(key >> ADDER_BITS) == 0` and, because
   ADDER_BITS+DP_BITS is 76 at A=2^16 and 80 at A=2^20, both cases take the
   same two-test form with a mask on the second limb. Only the shift and mask
   immediates change. The number of executed instructions is identical, and
   the rare-path `point_is_dp` maximum is unchanged for the same reason.

`k_mul`, `k_sqr`, the subtractions, `batch_step`'s five multiplications and
one squaring per walk step, the filter comparisons, the bookkeeping, the loop
control and the rare paths are untouched. The ADDER table is data, not code.
Consequently the per-step instruction counts, the anchor counts and the
forced rare-path maxima above apply at A=2^20 without a new measurement. No
new measurement is claimed: this candidate uses inherited Callgrind counts
and an instruction-for-instruction source argument. The only quantities that
change are the analytic rates that contain A, the derived truncation
thresholds, and the set-up and storage sizes; `accounting.py` evaluates all
of them exactly.

## Algorithm, fixed costs and changed parameters

Use the protected p, a=-3, b, r and generator P without alteration. Handle
Q=identity by returning zero; otherwise validate the curve point. Let a
class be {W,-W}, key(W) the integer Montgomery representation of its affine
x-coordinate, and canon(W) the representative with even Montgomery y.
Let sigma(W) be the sign taking canon(W) to W. The walk is

    W -> W + sigma(W) R[key(W) mod A].

Here **A=2^20** adders, each R_i=[c_i]P+[d_i]Q with fresh uniform scalar
coins. The table contains both signs. There are K=2^32 walk slots, arranged
into batches of 2^16 for simultaneous inversion. Distinguished points have
bits 20..79 of the key zero, so theta=2^-60 under the stated class model.
Each trail starts from fresh independent a,b with W=[a]P+[b]Q. It ends at its
first distinguished point, records it once with its coins, then starts a
fresh trail. Coefficients are reconstructed by replay only at a match.

The exact short-cycle rules, test ordering and replay algorithm are those
fully specified in `predecessor-proof.md` Sections 1.2--1.3 and are unchanged.
In particular, the distinguished-point test precedes cycle tests at every
generated point, including look-ahead and escape points. A detected two-cycle
or cycle of length 3..64 escapes by doubling the canonical representative of
its least-key member. The two 64-bit filters only trigger exact checks. Every
trail is limited to TRAIL_MAX=64*2^60+64 processed steps; the first key match
is resolved by two replays, each limited to 65 times its original live
length. Degenerate matches and exceeded caps fail the run. Each derived
scalar is verified by [x]P=Q before output.

Set N=(r-1)/2+1 and

    tau_target = 496250/1000000,
    T_min = ceil(sqrt(2*N*tau_target)),
    L = ceil(T_min/K), T=K*L,
    S = K*(L+TRAIL_MAX).

Every slot takes L main-phase steps and finishes its current trail. Thus S
bounds all processed steps, including persistent steps and finishing tails.
Sampling caps, maximum distinguished records, abandonment cap, probe cap,
restart costs, bounded replay, verification and output remain unchanged as
formulas, evaluated at this S. The point-order certificate proving
ord_r(2)>2*TRAIL_MAX is recomputed in `accounting.py` and still holds.

The common path costs at most C_COMMON=836.842722 instructions per processed
step under H3, with the same derivation as `70ff82dd`: the measured plain
round at every step, the 3-instruction other distinguished-point edge, 1/64
for wide batch control, and the measured anchor-round excess only in the
scheduled anchor rounds. The larger adder table changes none of these
instruction counts.

### A-dependent rates and reserves

The fresh-source probabilities in `predecessor-proof.md` Section 2 are stated
for a uniform index in [0,A); replacing A by 2^20 changes only their numerical
values. The recomputed quantities are:

    P_PRIME      = 64/(|Gamma| A)  = 64/2^21      = 3.0517578125e-5,
    RATE_TWO     = 1/(|Gamma| A)   = 2^-21        = 4.76837158203125e-7,
    RATE_LONG    = 3/A^2           = 3*2^-40      = 2.728484105318785e-12,
    epsilon      = 10/A + 1152/A^2 = 9.537790901958942e-6.

A fresh-source two-cycle entry costs at most w2=83,518 additional
instructions. Longer entries cost at most wL=2,789,687; false filter
episodes, including their speculative-exposure debt, are charged
256*(w2+wL); merged steps are charged wL. These are the same inclusive
forced-path/control allowances as before. The runtime still uses **one total
rare-work meter**. It does not need to recognize whether a point is globally
fresh, or to classify observed events as fresh versus merged.

The split-concentration argument is unchanged in form. With eta=2^-16,
mu=S*RATE_TWO and n2=ceil((1+eta)*mu), the exponential-moment bound gives a
tail below exp(-mu*eta^2/(2+eta))<2^-100 and charges

    TWO_CHARGE = n2*w2/S = 0.039825093... instructions per processed step.

With c2=TWO_CHARGE/S and the unchanged residual categories,

    E0 = (3/A^2)*wL
       + 2^-63*256*(w2+wL)
       + (tau_exact*TRAIL_MAX/S)*wL
       + (merged steps),   tau_exact = T^2/(2N),
       = 7.611696811906904e-6 instructions per processed step,

    RESERVE = TWO_CHARGE + 849*E0 = 0.046287424... instructions per step,

    PER_STEP = C_COMMON + RESERVE = 836.842722 + 0.046287 = 836.889009...
               instructions per processed step.

The larger A reduces both the two-cycle rate and the longer-entry rate
compared with the predecessor, so the reserve falls from 2.291541 to
0.046287 instructions per step. The 849 factor and the 2^-100 Chernoff
allowance are those of the predecessor. The abandoned-trail union bound uses
P_PRIME=3.0517578125e-5 and is smaller; the C_AB=368 choice still satisfies
e*C_AB*P_PRIME<=1/2, and k_ab=TRAIL_MAX//368 is unchanged.

## Success bound

The collision model, the deferred-exposure filtration, the escape algebra, the
replay and merge arguments and the failure decomposition are unchanged from
`predecessor-proof.md`. The conditional success lower bound is

    1-exp(-tau_exact*(1-epsilon))
      - 1/849 - 2^-100 - 2^-20
      - DELTA_ABANDON - DELTA_CUT - DELTA_MERGE
      - DELTA_CAP - DELTA_PROBE.

All terms are recomputed by `accounting.py` at A=2^20, tau=0.496250 and the
new S and T. The computed values are

    tau              = 0.496250,
    epsilon          = 9.537790901958942e-6,
    model            = 1-exp(-tau*(1-epsilon)) = 0.39118769872872783,
    DELTA_ABANDON    = 7.225802262629322e-8,
    DELTA_CUT        = 5.865725748644166e-17,
    DELTA_MERGE      = 1.256918430747597e-14,
    DELTA_CAP        = 2.444052392423174e-18,
    DELTA_PROBE      = 2.5021148200435266e-44,
    lower bound      = 0.390008816494845 >= 0.39.

The margin over the required 0.39 is 8.816e-6. It is a margin in the
conditional model, not a confidence interval for finite experiments. The
predecessor's 124-fold residual-expectation sensitivity statement does not
apply.

## Complete cost and success ledger

The ledger charges (COMMON+B_rare)/1000 for all walk steps and metered rare
work, then **adds**, rather than amortizes away:

* A adder constructions and K initial trail starts, each with the complete
  capped scalar-sampling and two-scalar-multiplication price;
* up to DP_MAX+K restart constructions, DP_MAX=2*floor(S/2^60)+2*K;
* DP_MAX table operations at 2^21 instructions each, including bounded
  probing, loads/stores, key comparisons, addressing and control;
* two bounded coefficient replays at the original generous unbatched
  operation count, followed by solving, verification and serialization;
* input validation and setup allowances.

All processors' work is summed. No future-target amortization or free advice
is used. The storage formula is

    2*DP_MAX*128 + K*512 + A*2*64 + 2^24 bytes.

It includes distinguished records, live walk states and batch buffers, both
signs of adders, and a fixed reserve. With A=2^20 the adder storage is 2^27
bytes, still dominated by the distinguished records. The ledger evaluates to

    L        = 55812290182901204270995371530,
    T        = K*L = 239711961367335180799997992737773191168 (log2 127.4945695),
    S_max    = T + K*TRAIL_MAX,
    step     = 836.889009363... instructions,
    total    = S_max*step/1000 units,
    log2 total = 127.23767771161744,

and `accounting.py` asserts that the total, rounded up at the sixth decimal,
equals the declared score 63618839/500000. The component logarithms are

| Component | log2 units |
| --- | ---: |
| walk steps | 127.237678 |
| adder table | 37.961183 |
| initial walk starts | 49.961183 |
| trail restarts | 86.455753 |
| distinguished points | 79.528785 |
| collision replay | 85.399650 |
| solve, verify, output | 16.288525 |
| input validation | 14.034216 |

Peak memory is at most 2^76.49457 bytes and preprocessing at most
2^49.961536 units, both recomputed and included.

## Evidence scope, execution and limitations

The inherited field assembly, measured production-curve common paths,
enumerated forced rare paths, replica results and reduced-order full-procedure
records are included with attribution. They support H3 and the model at their
stated scopes. They are not new measurements by this candidate. The
instruction-count invariance under A is a source-level argument about the two
A-dependent sites; a reviewer can verify it in `measurement/step.c` and
`measurement/p256_sub.h`, and the inherited Callgrind totals remain the
pricing evidence for the unchanged instruction stream.

One declared combined experiment uses the verifier's bounded credential-free
runner and returns two separately labelled results. `accounting.py` recomputes
the new numerical ledger, concentration threshold, failure terms and
point-order certificate. `experiment.py` is the inherited complete
reduced-order procedure, with the declared collision parameter 0.4965 (the
available input grid has spacing 0.0001). The latter is close to but not
identical to production tau=0.496250; its finite-run success is diagnostic
only. It does not implement the production instruction meter or execute a
full-size attack. Its actual successes, failures and replay checks must be
reported as observed.

The remaining material uncertainties are production-scale collision mixing,
conditional fresh-source orientation, the validity of the fresh-window
domination model, analytic handling of long trails and rare returns, and
universal validity of the inherited measured kernel bounds and of their
A-independence argument. The larger adder table does not resolve those
uncertainties; it only lowers the analytic rare-work rates. Memory remains
exponential and impractical. No claim of practical compromise, unconditional
improvement in P-256 security, or human certification is made.
