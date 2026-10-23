# GLV fixed-retention rho v12: conditional 130.68

## Scope and unchanged retention

The claim is success >=0.40 within 2^(3267/25) target-operation equivalents,
conditional on the explicit production laws below, not a production execution.
Use protected secp256k1 p,r; N=(r-1)/6, C=1000 (COST_MODEL.md), R=2^20,
(d,M,S,workers)=(64,2^72,2^85,1024), s=2^64, L=40s, nu=137/100.
T=ceil(sqrt((98859/100000)*nu^2*N)); K=T+1024L.
All ordinary/rewind/escape/drain ticks consume K.
Retain all eligible records until failure BEFORE exceeding M: no eviction,
replacement, adaptive density or filter change. Event/start, exception and replay
caps also fail closed. This retention architecture has not been redesigned.

## Actual step and labels

Canonicalize nonidentity points using min(x,beta*x,beta^2*x), even y; alignment
u is one of +/-1,+/-lambda,+/-lambda^2. Branch=bits64..83 of canonical x;
DP is its zero low64 word. Start with independent full 256-bit a0,b0 and
point=[a0*x_secret+b0]G; jumps have independent full labels aj,bj.

```
STEP(state):
  check two-limb meters; unpack live x,y; load jump[branch(x)]
  batched_affine_add: six field subtractions, six reference products
  normalize three x images (two EXTRA beta products), normalize y
  two 256-bit argmin compares/selects; parity-select y and alignment
  if DP: emit full-start record, compare/replay/verify candidates, terminate
  if new class in ordered ring[8]:
    select cycle minimum; rewind using charged STEP ticks
    double minimum (two ticks), clear ring; test DP again
  append/shift ring; Brent_compare(point, ordered ring, ring length)
  if expanded repeat: terminate detector; else update checkpoint/meters
```

Replayed labels update (a,b)<-u*(a+aj,b+bj) mod r, or 2u*(a,b) at doubling.
A match yields (b1-b2)/(a2-a1); reject zero denominator and verify the target.
Complete exceptional repair is capped. Toy labels update live; production
replays. Independent 192-bit full-point fingerprints are not DP bits.

## Collision-selected joint renewal law: measured, then conditional transfer

Each toy trial selects its FIRST shared-class visit between distinct trails before
T, NOT a uniformly sampled pair. Keep both local tick positions, times since last
escape, all subsequent escape times, terminal paths and useful endpoint match.
Dyadic histogram bins are [0,s),[s,2s),...,[32s,infinity); paired age bins stay
joint. Completed inter-escape gaps are histogrammed; terminally censored exposure
is retained in gap_exposure. Gaps alone are NOT fitted as uncensored samples.

For B branches, the diagnostic null hazard is e=1/(12B), DP probability 1/s,
q=(1-1/s)*(1-e). DP precedes escape. At h=m*s for m=1,2,4,40, null joint
survival before min(h,DP) is
  J=[1-(1-1/s)*e*(1-q^h)/(1-q)]^2.
Actual paired survival counts no escape and no abort before that horizon; an
earlier DP ends exposure successfully, not as artificial right censoring. The
square is a NULL diagnostic, not a claim of independence for selected halves.
All ratios use a directed-UP denominator-10^12 enclosure of J.

For n selected trials use eps=ceil(10000*sqrt(5/n))/10000. Hoeffding, conditional
on the eight fixed tables AND selection masks, gives simultaneous coverage >99%
for 96 tail cells and24 completion cells:120*exp(-10)<0.01, checked by a positive
rational Taylor sum. It covers those finite weighted mixtures only, NOT fresh
production maps. Take the minimum of (k/n-eps)/J and useful-completion lower
fractions, clamp to [0,1], then round DOWN to denominator10000:
eta=8521/10000. Raw worst tail/null=0.9798861865; selected trials=18079.

H-Rom posits the useful-core first-collision lower law P0=1-exp(-0.7*T*(T-1)/(2N)).
H-window now explicitly requires conditional useful selected-pair completion
>=eta, before separately union-debited adverse transport paths. This is a
production-transfer hypothesis supported at finite scope by JOINT selected
counts, not an inference that coordinate DP density controls history. The ledger
uses eta*P0, NOT an unjustified discount inside the exponent.
H-window further requires conditional selected-window spatial domination for
four pieces of length<=L: at most8L^2 relevant pairs, total collision mass<=8L^2/N;
and endpoint window debit64s^2/N. Counting pairs alone does not prove either
probability for this walk. The production cap loss44/2^24 is also conditional.

Every terminal route is separately counted: DP, identity, detector, exception
cap, trail limit, pool full and event cap. All24 observed path/debit ratios<=1;
maximum is0.051. Terminal-event rate/(1.1*ticks/s)<=0.917897. There are489 detector
aborts,319 identity terminals,719 exceptional repairs, zero exception-cap,
trail-limit, pool-full or event-cap aborts. Detector ceilings at these tiny sizes
clip to1, so their domination test is VACUOUS as a rare-event calibration. No
zero count proves the production cap loss. Scaled assumed detector ceilings for
common width8 decrease59.2593,56.3877,37.5238 across scales; measured detector
probabilities do NOT decrease monotonically. We do not claim empirical monotone
convergence or transfer to spacing2^64. H-rate's endpoint AND all-terminal first
moments <=(11/10)*K/s remain conditional, separately tested at toy scales.

Markov then debits (11/10)*K/(s*M) for pool abort and
(11/10)*K/(s*(S-1024)) for starts, without independence. Additional debits are
S(S-1)/(2*2^192),1/r,bounded scalar rejection,4RK/(N*2^16),8L^2/N,64s^2/N,
and44/2^24. Exact union lower eta*P0-sum(losses)=0.404614009730 >=0.40.
Neither this margin nor finite confidence is a production transfer-error bound.

## Code-level scalar-64 ceiling (complete 34-row table in ledger.py)

This is a 256-bit secp implementation, NOT a 384-bit Montgomery kernel. COST_MODEL
pins EXTRA 5x52 multiplications at179 instructions; six reference products are
already in C and are not charged again. Scalar arithmetic, compare, branch,
select, load/store and addressing primitives cost one instruction each.
normalize52_pack ceiling80=2*5*3 carry instructions+2*5 top folds+2*5 subtraction/
select+12 pack+9 traffic+9 control. modsub52 ceiling96=5 sub+1 mask+5 AND+5 add+80
normalize. unpack52=5*3=15; cmp256=4 SUBS/SBC+1 CSET=5. These are constructive
analytical ceilings, not new hardware measurements.

Routine totals (34 primitive-count rows in ledger.py): add/repair594;
canonicalize801; add62; step406; escape279; Brent258; meters32. Sum2432<=8000,
margin5568. The complete table ties EVERY row to the pseudocode routine and its
primitive operation count; certificate derives the sum before testing8000.
Heavy escape/control rows are charged EVERY tick; doubling uses two reference
ticks. No production timing measurement is claimed.

Reprice ALL8000 (not only2432) at6/5: E=9600. Main=K*(1000+9600)/1000.
Other base rows, all repriced6/5: lookup S*2^18/1000, starts S*4096, maintenance
2^96, replay2^110, exceptional2^40, table2^32, constants/sampling/output2^20.
Replay: (2^24+2^16) requests*2L ticks*1024 plus4096 final verification fits2^110.
Tree height<=146; routed96-bit byte-pointer traffic is charged.
raw log2=130.5594102322; W=(27/25)*raw has log2=130.6704415446.
Claim130.68 is >log2(raw)+0.10 and <log2(W)+0.06; exact powers check the inequalities.

Records are192B=fingerprint24+full scalars64+length16+pointers36+flags2+metadata50.
Peak=M*192+2^30=3*2^78+2^30<2^80 BYTES; auxiliary allocation includes table,
code, workers and scratch. Memory is reported, unscored; nonuniform advice zero.

## Reproduction and scope

Run python3 -B scripts/toy_sweep.py, then scripts/verify_certificate.py --write
and without --write.24 configs*8 tables*250 trials use N=216,908,87328,
spacings1,2,8/16 and live widths8/64/256. Full GLV enumeration and2000 independent
EC additions/group check the coordinate oracle.18057 recoveries, no wrong
answers;24/24 binding envelopes pass. Arrays preserve all counters. Finite shape
tests do not prove production magnitudes, rare basins or instruction timings.
