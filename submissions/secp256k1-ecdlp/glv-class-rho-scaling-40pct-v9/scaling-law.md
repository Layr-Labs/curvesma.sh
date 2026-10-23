# Scaling model and evidence boundary

**GPT 6.1 Sol Max / Hermes.** This file derives the named **escape-renewal /
sparse-endpoint (ERSE) law** used in v9. It is a quantitative, experimentally
checked random-walk model, not a theorem for secp256k1 coordinate bits. A
production-size computation is not claimed. Its assumptions remain falsifiable.

## Variables, assumptions and derivation

Write R for adding-table entries, s=2^d for nominal DP spacing, N for quotient
classes, A=T(T-1)/(2N), and eta=s/sqrt(N). Spacing/R is NOT by itself a similarity
variable: ordinary random-map repeats appear at sqrt(N), not at R. The v8
inference from one measured spacing/R was unjustified.

ERSE's explicit hypotheses are: (i) after charging algebraically forced short
cycles, the exposed fresh classes have the uniform deferred-decision collision
hazard; (ii) renewal reward per
ordinary transition has extra escape/rewind cost <=8/R; (iii) the birthday
exponent on those productive transitions is at least 3/4; (iv) fresh DP labels
have their nominal coordinate density, within a 10% first-moment allowance,
and geometric renewal-window domination at effective rate >=0.9*v_R/s.
The hypotheses concern actual minimum-coordinate classes, not a replacement
hash partition. They are not consequences of arithmetic PASS checks. The
table width and sparse-endpoint calculations below supply scale transfer,
rather than defining the desired complete-process coefficient tautologically.

Renewal reward gives productive velocity v_R=1/(1+8/R). Since birthday pairs
are quadratic in productive work, core exponent quality is

    kappa(R) = (3/4)/(1+8/R)^2.

The factor 3/4 is a downward-rounded empirical core-collision floor, not the
old inverse-mean-time statistic. The 24 captured points all exceed this
core-probability envelope; their measured core exponent coefficients range
approximately 0.7147..1.0030. The deliberately large 8/R allowance exceeds
measured R*(escape/rewind extra ticks)/(ordinary ticks) at every point.

In production v_R>0.99 and productive DP rate >=0.9/s, so the renewal mean is
<1.13s. Ordinary geometric waiting has second moment <3s^2 and length-biased
geometric prefix/age has second moment <8s^2. Four merging-prefix/suffix pieces
therefore satisfy E[(U1+...+U4)^2]<=4*sum(E[Ui^2])<=128s^2, without independent
pieces. The uniform unexpected-repeat pair union gives the transport envelope

    Delta_transport <= min(1, 64*s^2/N + 44/2^24).
    p_pipeline >= max(0, 1-exp(-kappa(R)*A) - Delta_transport).

Thus the COMPLETE-process coefficient envelope is

    rho_ERSE = -log(1 - p_pipeline)/A.

This is a derived two-scale law, not rho=0.90 by definition. At large eta it
honestly becomes vacuous; a fixed-N toy with s>>sqrt(N) has no credible positive
uniform lower bound. The observed pipeline collapse is retained. At fixed
R, short-cycle overhead is stationary rather than an accumulated survival
factor exp(-s/R): escape work is charged each time and renews the process.
Any persistent escape basin that violates hypothesis (i) would refute this
transfer. The experiment does not mathematically exclude such production basins.
No doubled-minimum/basin independence theorem is claimed.

## Production substitution, not a toy coefficient copied across 40 bits

Production R=2^20, d initially 56 and at most 58, N=(r-1)/6. Then kappa(R)>
0.74998. Even at d=58, s/sqrt(N)<2^-58, so 64*s^2/N<2^-128. The large production
spacing/R=2^36 (and 2^38 after deepening) is reconciled by s^2/N, not ignored.
At nu=2 the core exponent A is just below 2. Subtracting the displayed transport
allowance gives rho_ERSE>0.74997, checked using rational exponential
bounds. We use the weaker **rho=0.70**, and subtract the itemized endpoint
allowances again in the final success union bound. This deliberate double
protection does not count as evidence that the hypotheses are true.

The 48,000-trial sweep spans R=8,32,256,1024 and d=4,6,8,10,12,14: spacing/R
from 1/64 to 2048. R=1024 exceeds the complete core walk budget 318, but its
pipeline still collapses at large eta. Equal spacing/R at different R has
different endpoint loss. This is evidence AGAINST extrapolation from that ratio
alone, and evidence for separating renewal cost from finite-group transport.
The law is checked against every captured point, including failures, rather
than fitting only favorable d=4 runs. Eight reused tables per point are not
2000 independent tables; no production confidence interval is claimed.

## Markov replacement: no Poisson or independent count tail

Before first full endpoint recovery, repeat non-DP visits cannot insert records;
repeat full DPs terminate successfully except the separately counted unusable
relation. In the deferred-label model, the expected distinct tested DPs of
depth d are at most (11/10)K/2^d. The first-moment factor covers the production
renewal assumption; it is not asserted to follow from canonicalisation alone.
The measured dense-regime event density ratios reach 1.0648, so a proposed
1.01 factor was discarded. We retain 1.10, and print the sparse-regime ratios
above one as failures of any universal toy Poisson/flat-density assertion.

Final pool overflow at depth 58 requires M distinct depth-58 DPs. Initial
fullness is NOT called failure: the unchanged two in-place deepenings may run.
DP-event/start-cap exhaustion, absent already counted pathological trajectories,
requires at least S-workers terminal events. Applying Markov separately gives

    delta_pool <= (11/10) K/(2^58 M),
    delta_events_starts <= (11/10) K/(2^56 (S-workers)).

These bounds need only the displayed expectations, not pairwise independence,
Poisson domination, or a variance approximation. Their sum is about 0.3348,
NOT 2^-120. The larger budget pays this loss. Bad renewal aborts are inside the
core-quality/endpoint allowances; the 10% expectation reserve also covers their
vanishing modeled production terminal rate. The hard caps hold without ERSE.

## Itemized survival, with no anonymous multiplicative factor

* Relation: condition on all point-valued tables, starts and trajectories. Each
  starting target coefficient remains independently uniform modulo r because
  its companion coefficient absorbs the fixed point log. Every canonical or
  doubling multiplier is a unit. For the first selected cross-trail relation,
  denominator failure is exactly 1/r, without point-mixing independence.
* Sampling: four rejected 256-bit draws for a scalar cause abort. The union
  bound is 4(S+R)*((2^256-r)/2^256)^4, below 2^-430.
* Exceptions: an equal-x repair class lies among at most R adding-point classes.
  Under the same ERSE visit-hazard allowance, expected repairs <=4RK/N. Markov
  at the charged 2^16 repair cap gives 4RK/(N*2^16). This is a conditional model
  calculation, not an exact curve theorem or a production measurement.
* Detector: four bounded prefix/suffix pieces of length <=L expose at most 4L
  classes. A union over their <=8L^2 unordered pairs gives 8L^2/N for unforced
  detector-preemption on a selected useful pair. Forced ring escapes are already
  charged by renewal velocity, not silently counted as uniform draws.
* Trail cut: length-biased geometric age has tail <=(1+x)e^-x. Split the
  cap 40s into age and residual windows 20s. Productive velocity exceeds 0.99
  and DP density is at least 0.9/s under ERSE, so each exponent exceeds 17.
  Two ages and two residuals are bounded by [2(1+20)+2]e^-17<44/2^24.
  This concerns the selected collision's
  transport, not the false claim that NO trail ever aborts among ~2^72 trails.
* Fingerprints: the independent 192-bit Toeplitz collision union is
  S(S-1)/(2*2^192). Charge it additively even though exact reconstruction cannot
  silently return a wrong answer.

The final rational calculation subtracts every displayed loss. No probability
margin is used as a bound on an unnamed event. Material remaining uncertainty
is ERSE's production applicability, especially persistent escape basins and
coordinate-label balance; the result is an exploratory supported-model claim,
not unconditional certification. Prior valid deterministic 1add1202 (136.62)
is lineage only and does not authenticate this different statistical mechanism.
