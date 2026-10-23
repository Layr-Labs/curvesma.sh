# Fixed-spacing GLV class-rho v10b: conditional 128.21

## Claim, model and exact parameters

Under the explicitly conditional H-Rom/H-rate/H-window laws below, success is
at least 0.40 within 2^(12821/100) classical target-operation equivalents.
This is conditional analysis, not production execution.
C=1000 is fixed by tracks/docs/COST_MODEL.md; unrelated instructions cost R/C.

Let r be the protected secp256k1 order and N=(r-1)/6. Set
(d,M,S,nu)=(64,2^72,2^85,153/125), DP=2^64, R=2^20 and workers=1024.
T=ceil(sqrt((98859/100000)*nu^2*N)) is
169064858030058607474129302997935005868.
The actual-spacing cutoff L=40*DP=737869762948382064640 and
K=T+1024L=169064858030059363052766562141169197228.
No adaptive depth is used.

## Mechanism and faithful labels

The six nonidentity GLV units are ±1, ±lambda, ±lambda^2. Choose the least
x,beta*x,beta^2*x and even y. Each start has independent full scalars (a0,b0)
with point=[a0*x_secret+b0]G. The coordinate branch is bits 64..83 of canonical
x, disjoint from DP bits 0..63. Uniform full-scalar jump labels (aj,bj) form
R table points. For canonical alignment u after adding the selected jump,
labels update to (u*(a+aj),u*(b+bj)) modulo r. Equal canonical points give
x_secret=(b1-b2)/(a2-a1) when the denominator is nonzero. Full target verification
is mandatory. Identity and equal-x additions use complete exceptional branches.

Keep an ordered ring of eight canonical classes. On its closed-cycle repeat,
rewind to the minimum of that cycle, double there (two ticks), and clear the
ring. Test DP before escape and after doubling. Expanded-state Brent compares
(point,ordered ring), aborting other repeats; it does not choose escapes.
Stop starting after T ticks, then drain the in-flight trails to their L limits.
Escapes, rewinds, aborted trails and drains all consume K, not only useful steps.
Production defers labels to charged replay; toy labels update live. Records store
both full starting scalars, a replay length and an independent fingerprint.

## D1: no deletion, finite capacity and the actual event process

Use one fixed coordinate DP predicate throughout. A terminated endpoint never
becomes ineligible at a later filter switch: there is no switch. Insert unique
records into the capped pool; stop with failure BEFORE an insertion would exceed
M. Never evict, replace or deepen. Event/start/replay/exception cap excess also
fails rather than allocating hidden work or memory. This specifies the full-pool
policy that the prior draft lacked; no transport term is used to excuse eviction.

Let D count emitted endpoints under the K cap. H-rate states E[D]<=b*K/DP,
b=11/10, for THIS fixed-density process. The analogous terminal-event expectation
is also bounded by b*K/DP; this includes the terminal mechanisms and is an explicit
model premise, not inferred from DP density alone. At most D records can be live.
Thus pool abort is contained in {D>=M}, and Markov gives

  fill = K/(2^64*M) = 0.0019407693900088557...
  pool loss = b*fill = 0.0021348463290097415... <0.01.

After the initial workers, starts consume terminal events, so the corresponding
Markov debit is b*K/(2^64*(S-1024))=0.00000026060135852169696... <0.005.
Both denominators are DP, not old depths58/56; Markov needs no independence.
H-rate itself still needs production extrapolation.

H-window supplies the separate selected-pair endpoint debit 64*DP^2/N,
approximately 1.1284745768e-36. It concerns renewal/age and four transport pieces,
not pool replacement. The four <=L pieces yield the detector debit 8*L^2/N,
approximately 2.2569491536e-34. These terms are recomputed after L and d move.
Remaining debits are S(S-1)/(2*2^192), relation failure 1/r, bounded scalar
rejections, exceptional repairs, and trail cap 44/2^24. Fingerprints hash full
coordinates with an independent universal Toeplitz map; DP bits do not supply
fingerprint entropy. Reconstruction prevents silent false answers.

H-Rom is the useful-core collision law with coefficient 7/10. Directed rational
alternating series enclose exp(-0.7*T*(T-1)/(2N)). Subtracting ALL listed losses
gives 0.4023768860995998..., hence >=0.40. The margin over 0.40 is
0.0023768860995998..., in addition to the named 0.01 reserve above the protected
0.39 floor. Neither margin bounds error in transferring the hypotheses.

## Computation, counters and honest peak memory

The original extra-step rows total 925 instructions. Add a 24-instruction ceiling
for two-limb length increments/carries, comparisons, checkpoint copies and replay
bookkeeping: total 949<=E0=950. Reprice E=6*950/5=1140. DP64 is a one-limb zero
test; the branch still starts at bit 64. L no longer fits a 64-bit length: it uses
TWO limbs, explicitly. K fits two limbs. No arbitrary-width unit-cost counter is
assumed. Work main=K*(1000+1140)/1000. All other rows receive 6/5 too:

* lookup/transfers = S*2^18/1000; starts = S*4096; both use NEW S=2^85;
* retention initialization/maintenance =2^96; zero deepening rows;
* complete replay, scalar checks and key matches =2^110;
* exceptional arithmetic =2^40; table preprocessing =2^32;
* constant checks, sampling/output =2^20.

At most A=2^24 match requests plus 2^16 exceptional requests replay <=2 trails
of length L. Even a loose 2^10 equivalents per replay tick plus 2^12 per final
verification fits 2^110; the certificate checks this with the NEW L. Fresh-start
full scalars, complete candidate verification and all failure paths are charged.
The pool red-black tree height is <=146; 96-bit byte pointers cover its addresses.
The 2^18 instructions per event cover routed multiword accesses and balancing.

Raw repriced log2 work is 128.0884614260. Reserve W=(27/25)*raw has
log2 128.1994927383; claim 128.21 leaves 0.0105072617 bit and >0.10 bit above raw.
Exact power comparisons, not floating logarithms, check both gates and the
strict claim<log2(W)+0.06 window. Changing nu reprices T,K and every dependent row.

Each 192-byte record is fingerprint24 + full scalars64 + length16 + pointers36
+ flags2 + metadata50. The extra length limb consumes padding, not scalar entropy.
Memory rows are pool=M*192=3*2^78 bytes and auxiliary cap=2^30 bytes; peak is
906694364710972954771456 bytes, strictly below 2^80. Table, code, workers and scratch
share the auxiliary cap. Memory is unscored under COST_MODEL, but reported as an
actual hard allocation, never expected occupancy. Nonuniform advice is zero.
The prior 3*2^79-byte pool fits its published <2^81 cap, not this tighter cap;
2^79 RECORDS cannot be relabeled as 2^79 bytes.

## D2: rerun a binding three-scale sensitivity test

Toy orders1297,5449,523969 give N=216,908,87328. Each Group checks its own order,
all GLV classes and2000 independent EC additions. Every branch is live. Actual
x coordinates drive partition and DP; the oracle only accelerates addition.

The spacing matrix co-moves one sensitivity bit when production d widens, clamped
at the exact delta<0.30 feasibility limits (small0/medium1/large4). This is NOT a
literal simulation of a ten-bit production increase at tiny fixed N. Small stays
at spacing1; medium moves to spacing2; large tests spacings8/16 at widths8/64/256.
Variants at the same parameters are separate table batches, not extra densities.

All 24 configurations ran eight tables x250 trials: 48000 trials. T is the
smallest multiple of16 with E=1-exp(-0.7*T*(T-1)/(2N))-64*s^2/N >0.02.
Directed lower E is rounded DOWN to denominator10^12 and recomputed exactly.
Minimum E=0.025894125697; all 24 measured probabilities exceed their own E.
Toy/production transport ratios exceed10^30 in every row. There are18057 core recoveries,
zero wrong recoveries, 38945 escapes, 21688 rewind ticks, 489 detector aborts,
zero length aborts and17 degenerate endpoint pairs.

Mean bare-map repeat/sqrt(pi*N/2) spans0.119498..1.022378, NOT an exponent
coefficient. The separate full-pipeline fitted coefficient spans0.845673..1.129714.
Finite compatibility does not prove production mixing, H-rate/H-window, rare
basins, costs or large-pool behavior. scripts/verify_certificate.py certifies
conditional arithmetic and capture reconciliation, not organizer acceptance.
