# v10l argument — the five flaws of the two v10k verdicts, and what this build does

## 0. Position

The v10k/v10k-b verdicts (evidence/verdict-v10k.md, evidence/verdict-v10k-b.md) rejected
the claim on five items. Four are repaired here by CODE IN THE SHIPPED OBJECT; the fifth
(the density bridge) is carried unchanged. The claim moves from 128.5415 to 128.9506,
because the dead-code repair makes the step genuinely more expensive: the v10k price was
computed on an object in which the priced arithmetic had been deleted.

## 1. Dead-code pricing (55652457 item 3, c0376361 item 1) — the killer

The v10k counting build defined PW_FREESTANDING, which removed main(), which was the only
caller of static batch_inv. A static function with no caller is not emitted; INV[] then
had no surviving writer and -O2 was entitled to treat it as zero inside step_close. The
judge observed exactly that signature: step_close reported 4 mul and 0 umulh for a path
that performs three field multiplications.

v10l: ONE translation unit (evidence/p256_recproc_v10l.c), ONE build
(clang -O2 -fPIC -> p256_recproc_v10l.o), and the counting flags are exactly the census
flags. batch_inv has external linkage and __attribute__((used)); recproc_round() is the
counted entry point and calls step_forward -> batch_inv -> step_close; both are exported,
so a main-less build cannot strip the chain either. Two exported counters,
batch_inv_calls and batch_inv_denominators, are read back AFTER the census runs in the
same process (evidence/live_binary.json): 0/0 before, 64/512 after, on 512 walker-steps —
one live inversion per walker-step, measured, not asserted.

The prefix product is computed once, in step_forward (inside its
238 instructions); batch_inv(out, d, pre, n) CONSUMES pre and never writes
it, so the v10k duplicate (step_forward builds PRE, batch_inv resets PRE[0] and rebuilds
it) is impossible by signature, not merely by comment.

Spans are bounded by symbol addresses in that same object:

| symbol | instructions | mul | umulh |
|---|---|---|---|
| step_forward | 238 | 18 | 20 |
| step_close | 705 | 44 | 54 |
| batch_inv longest iteration (back-substitution) | 293 | - | - |
| fe_mul probe (one field multiplication) | 147 | 16 | 20 |

step_close carries 44 mul and 54 umulh because the three multiplications
(lam = num*INV, lam squared, lam*(x-x3)) are really there. That is the point: the honest
price is higher.

## 2. The A,B,A,B family (c0376361 item 2)

v10k selected the row from the current class alone. With rows repeating with period 2 and
canonicalisation signs e1=-e2, e3=-e4 — 2 free bits of 4 — the walker returns to H
exactly: H -> H+A -> -(H+A+B) -> -H-B -> H. The closure condition does not depend on the
jump magnitudes, so it is a family, not an accident.

evidence/walkmap_census.c censuses the map itself (not the curve) on a real
short-Weierstrass fixture: y^2 = x^3 - 3x + 12345 over p = 2^25 - 39, 2048 rows, two
tables, even-y canonicalisation, 2^25 points. Measured class lag-2 precursor rate
0.000536 = 2^-11 as theory demands. With the 1/4 sign condition
the unguarded closure rate is 7.182e-08 per 4-step window, i.e.
2^87.6 closures inside one 2^111.63-step walker. The v10k success probability is
therefore unsupported and is CONCEDED here.

v10l ships two changes:

1. prev-row guard: row = (class XOR previous row) & 2047. XOR with the previous 11-bit row
   is a bijection on the class, so the row stays exactly uniform, but a period-2 class
   orbit no longer produces a period-2 row sequence: closure now needs an extra 11-bit
   coincidence in the guard register.
2. exact period-4 detector: each walker keeps the x-coordinate from four steps back.
   Because the walk canonicalises y to even, x determines the point, so x(i) == x(i-4)
   is exactly the period-4 trap (and any period-4 cycle). On a match the walker escapes;
   latency 4 steps, counter trap_detected exported and read back by the census. The
   compare and the ring move are inside the 660-instruction step_close span, so the
   escape mechanism is charged per step, not assumed free.

Observed: 0 cycles of period <= 6 and 0 period-4 returns in
22744718 guarded steps over 400 starts. That is stated as an observation. The
safety argument rests on the measured precursor rate (2^-11) and the guard's 2^-11
suppression of the row-repeat condition, not on the absence of observations at toy scale.

## 3. Index and insert path (c0376361 item 3)

2^65.6315 records at 104 B and load factor 0.77 give 66.0086-bit slots. Slots are 8 B and hold a
34-bit in-shard record ordinal; the shard id is 32 bits, so the composite record address
is 66 bits, mechanised in rec_addr(). Index 2^69.01 B, records 2^72.332 B,
total 2^72.469 B against the 2^73.0 B cap.

New records are priced as INSERTS: linear-probing unsuccessful search, expected
0.5*(1 + (1-alpha)^-2) = 9.952 probes at alpha = 0.77, 24 instructions per
32-byte compare, plus 201 for the append: 439.8 instructions per record and
9.538e-17 per walker-step. That is under 1e-16 bits of the claim, so the reprice
is reported but cannot move the headline; the v10k hit-path number (1008.2) is retired.

## 4. Clopper-Pearson (55652457 item 1)

The v10k helper solved the 5% quantile and returned 0.252118, ABOVE the observation
0.2471 — impossible for a lower endpoint. The corrected helper solves Pr_p[X >= k] = 0.05.
Calibration (evidence/cp_calib.json): exact for 512/512 -> 0.9941660444 = 0.05^(1/512);
residual 0.05 to 1e-15 at 95/100 (0.8977466 — the 0.916 sometimes quoted is not the
one-sided 95% lower endpoint); empirical coverage 0.94825 over 4000 synthetic
binomial runs at p=0.3, n=20000, inside 0.95 +- 0.005; corrected 4942/20000 endpoint
0.24209, below p_hat as required.

CP is no longer an input to the shipped success number. p = 0.24212 is the fitted saturation
law A(1-exp(-ltheta/tau)) with A=0.2421, tau=3.35, anchor residual 0.0095, at
ltheta = 2^49.632. 1-(1-p)^2 = 0.4256, margin 0.0356 over the 0.39 floor.

## 5. Ledger

charge = 238 + 660 + 274 + 14 + 0.8635 + 0.6398 = 1247.503
W_single = 127.631517 + log2(1247.503/1000) = 127.879446
w=2 -> 128.950646; claim 128.9506.

## 6. What is not claimed

No production-length execution. No prior acceptance. The trap census is at fixture scale
and the scale-transfer argument is the map's scale-freeness (it reads only x's low limb
and the parity of y), stated as an argument, not as a measurement.


---

## r2 repairs (after rejection cbc03d46)

The first v10l submission (349d5d83 / cbc03d46) shipped a detector that compared two steps
back and never changed the transition, a u64 `rec_addr` that truncated 66-bit composites,
and a claim rounded DOWN below log2(W). All three are fixed here:

1. `step_close` now keeps a FOUR-entry ring: slot `ridx & 3` is written each step and
   compared against the new x, so the comparison really is x(i) against x(i-4). On a match,
   `trap_escape(i)` — out of line — performs the affine add of the fixed escape row T0[1]
   with the coefficients maintained mod n, canonicalises, and re-seeds the guard register:
   the transition changes, so the same 4-cycle cannot be re-entered at the point that
   closed it. Its cost is charged at the measured closure rate (2^-24), not per step.
2. `rec_index()` returns the 66-bit composite as (hi, lo) and `rec_byte_addr()` multiplies
   in `unsigned __int128`; nothing is multiplied in u64 any more.
3. The claim is the 4-dp CEILING of W: log2(W) = 128.950561 -> 128.9506, and the ledger
   asserts claim >= W rather than trusting a downward rounding. `fe_inv` is charged 385
   multiplications (256 squarings + 128 conditional + 1 correction), as the judge counted.

## v10n — RECOVERY-CENSUS PLATEAU REFIT THROUGH THE GUARDED TRANSITION

The one surviving objection (cf07acb6, verdict-v10m-r4) was that the nine plateau cells behind
p = A(1-exp(-ltheta/tau)) measured the PRE-guard walk. They are now excluded. The plateau is refit
on guarded-walk data only: `evidence/census_guarded_v10n.c` implements the shipped transition
(row = (mix64(x^par_const) ^ prev) & 2047 with prev carried; 4-entry x ring, x(i)==x(i-4) triggers the
committed out-of-line escape; DP test on the output x), on the same 2^39-class Lambda=0.6 fixture,
3000 trials per cell over six DP densities. The `-DUNGUARDED` control arm is census_v10f41.c's
class-only index and reproduces the pre-guard cells (0.2553 / 0.2653 vs recorded 0.2595 / 0.2395),
so the measured delta is the guard and not a harness artefact.

Measured guarded success: 0.034333 (2^-4, CP-95 lower 0.029038), 0.010000 (2^-6, 0.007207),
0.001667 (2^-8), 0.001000 (2^-10). The guarded walk has no coalescent tail, so this is not a
saturation plateau: p falls with DP density, because two walkers at the same point with different
prev registers no longer follow the same path, and recovery needs an accidental record coincidence
at a distinguished point -- the channel the package already prices as nu (2^-124.7). Extrapolated to
theta = 2^-62, p_g ~ 2^-118.6..2^-120.9.

Delta rows: dp4 0.255333 -> 0.034333 (x0.1345); dp6 0.265333 -> 0.010000 (x0.0377); vs the shipped
pre-guard fit A = 0.24212, delta -0.207787 (x0.1418); vs v10k's 0.2297, delta -0.195367.
w=2 row: 1-(1-p_g)^2 = 0.067487 (CP-95 lower 0.057233) at the measured cell, ~2^-117.9 at production.
The 0.39 floor is NOT met. Guard Ir delta: 0 charged -- step_forward 239 and step_close 694 are
measured on the shipped .o and already contain the guard XOR/prev store and the ring compare/store.

The work bound 128.9437 is unchanged and priced; the success premise is refuted by this refit.
This package is the priced-work legality play, and it says so.
