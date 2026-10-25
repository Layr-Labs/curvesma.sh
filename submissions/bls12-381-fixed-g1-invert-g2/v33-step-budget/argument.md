# Argument: the integrated protocol with a per-EXECUTED-STEP budget (v33)

v32 (ticket bca1e72a, claim 126.9558) was rejected on four flaws. v33 repairs
exactly those four and keeps v32's machinery: `run_protocol` is still the only
path, DP-first termination, the twelve-shift recovery, the sampler-order control,
the in-flight gate, the instruction-counted pairing/output rows, regen + sweep.

Claim **126.582** = 63291/500, the least 4-decimal log2 c with
2^c >= W/C_REF where W is the exact INSTRUCTION total of
`evidence/ledger33.json` (W = 2^140.869663 instructions, W/C_REF =
2^126.581950 target operations). Emitted by `scripts/regen_v33.py`; no claim
literal is typed.

## 1. The budget is EXECUTED STEPS, not launches (objection 1)

v32 charged theta per LAUNCH and clamped the launch count, so with a per-trail
cap of 128*theta the work bound was ~128*S and `t >= S and not inflight` was dead
code behind the empty-in-flight break. The same seed now runs twice:

    seed 4242, S = 20000, theta = 2^5
      v32 rule (charge at launch, no step cutoff)  completed 654518   consumed 1341440
      v33 rule (charge per executed step)          completed  20000   consumed   20000

654518/20000 = **32.70x** the declared budget on that seed (the v32 trace the
judge read is the same regime), and v33's stop is a trace row:

    {"round":360,"t":20000,"event":"halt_step_budget","inflight":103,
     "launches":360,"completed":20000,"consumed":20000,"remaining":0}

The loop tests `remaining = S - t` BEFORE every step and executes at most one
step per trail per round, so the executed total over ALL trails -- abandoned,
in-flight and bootstrap ones included -- cannot exceed S; `consumed == executed`
is asserted, not assumed, and `launches <= ceil(S/theta)` follows from RULE 2
(a launch needs `remaining >= theta`) and is also checked (360 <= 625). The
in-flight cap is a hard gate in the code; on the shipped trace the maximum
in-flight is 103, so it never binds, and a control at cap 32 binds it (541
deferrals, completed still exactly 20000).

SUCCESS IS RECOMPUTED UNDER THE REAL BUDGET. The v32 number (0.75) was measured
at ~33x the declared work. Sweep through `run_protocol`, 24 challenges each,
stopping at the first verified recovery:

    S        success   cp95-lo   mean completed (<= S)   mean launches/bound
    200      18/24     0.5653    120.0                 7.0 / 7
    400      22/24     0.7602    150.0                 12.2 / 13
    800      24/24     0.8827    138.6                 16.5 / 25
    1600     24/24     0.8827    138.7                 16.7 / 50
    4000     24/24     0.8827    138.7                 16.7 / 125
    20000    24/24     0.8827    138.7                 16.7 / 625

The success numbers moved UP, not down, because the fixed protocol stops at the
first merge inside a smaller executed budget and the fixture's walk is
degenerate: it collapses onto few classes, so merges are frequent. That
degeneracy is the fixture's, and it is disclosed in section 6 -- the fixture
measures the MECHANISM (protocol order, budget, detection), not the production
rate.

## 2. P(D > R) on the ACTUAL event (objection 2)

v32 used Markov on E[D] to bound P(D >= 128 theta), which does not bound P(D > R).
v33 bounds the event that actually loses a trail. `routing.dp_hit` reads a FRESH
5-bit field of the class key at every step, so the per-step hit indicator is
i.i.d. Bernoulli(1/theta) and

    P(D > R) = q^R exactly,  q = 1 - 1/theta,

with R the trail's remaining lifetime at its launch, `R = min(cap, S - t_launch)`,
MEASURED per launch and carried in `launch_log`. Integrated over the launch
schedule (one launch per theta of executed work at most, so R runs over
theta, 2*theta, ..., S) the mean is bounded by the geometric sum

    delta = P(D>R) <= (1/(e-1)) * theta/S = 0.582 * theta/S,

solved to a fixed point with S itself (production: 0.582*2^80/S = 9.41e-15).
`evidence/tail33.json` measures the same row on the shipped run and checks it
under the same ceiling:

    R (measured, 360 launches)      mean q^R
    R as shipped, 64..4096          3.790e-4      <= ceiling 9.312e-4  ok
    R/2  (audit)                    1.254e-3      row RISES when R falls
    2R   (audit)                    4.779e-5      row FALLS when R rises

The row is a function of R, not a constant, and no inequality crosses events.

## 3. P(detect | merge), conditioned on merges (objection 3)

The v32 charge, 6/10 = 0.600 "retention", was `distinct stored keys later shared
/ reports` with stopping disabled: a COVERAGE fraction of a finite run, already
inside the birthday term, not a conversion probability. v33 conditions on MERGE
events and counts detection, through the same `run_protocol`, with stopping off
so the batch can collect merges inside the budget (each merged trail is cheap):

    merge events          1054
    detected              1054        (twelve-shift witness verified against z)
    P(detect | merge)     1.000000
    CP-95 lower bound     0.997162    <-- CHARGED, not the point estimate
    old shared fraction   0.600000    (same protocol, labelled: a coverage number)

`evidence/detect33.json` carries both numbers, labelled, and control (c) asserts
they differ. The j-histogram of the twelve shifts is recorded. The correction
moves the claim DOWN by 0.366 bits (S scales as 1/sqrt(p_det)): 126.9484 with the
old 0.600 charge vs the shipped 126.582. Charging the wrong object was hiding a
credit; v33 takes the debit away because the measurement says the conversion
factor is ~1 and the coverage deficit is the (1-exp(-.)) term. mu = 11/10 stays
charged: it is a 10 % work allowance, i.e. a debit, not a credit.

## 4. One division by C_REF (objection 4)

v32's six rows were already target-operation equivalents and the statement
compared their sum against 20000*2^c, normalising twice. v33 rebuilds the chain
so there is exactly one division line, `ledger33.total33`:

    rows are INSTRUCTIONS (integers; every judgement row rounded UP)
    W = sum(six rows)                      = 2^140.869663 instructions
    W_ref = W / K.C_REF                    <-- THE ONE DIVISION LINE
    claim = least 4-decimal c, 2^c >= W_ref = 126.582 (2^126.581950)

Control (d) re-emits the claim after multiplying the RAW instruction ledger by
C_REF: 126.582 -> 140.8697, a shift of 14.287700 bits against log2(20000) =
14.287712 (the residual is one 1e-4 claim quantum). A double-normalised chain
would move by 0 or by 2*log2(C_REF); the audit asserts neither.

## 5. Rows (instructions, and target ops after the one division)

    row                 instructions        target ops
    w_effective         2^140.869663        2^126.581950   S * step price
    abandon_row         2^49.597536         2^35.31        cycle-abandon charge
    aux_sum             2^70.789040         2^56.50        stores, starts, solve, verify
    pairing_row         2^26.213926         2^11.926       3892 Fp12 mul (Miller 1533
                                                           + fexp 300 + G2 226 + verify 1833)
    overflow_row        0                   reports capped at ceil(S/theta) literally
    stop_loss_row       2^109.056976        2^95.77        in-flight population + abandoned work

S = 74712357095124621315373948935154181646 (2^125.81) solves
(1-delta)(1-exp(-S^2(1-f^2)p_det/(2 mu N))) >= 0.392 with delta = 0.582*theta/S
and p_det = 0.997162, mu = 11/10, f = 3/8, N = (r-1)/12. Reports = 2*ceil(S/theta)
= 1.236e14; the pairing row is instruction-counted (Miller 64 DBL x 9 + 87 ADD x
11 Fp12 = 1533; final exp 300; G2 reconstruction 254 DBL x 24 + 127 ADD x 48 Fp mul
/ 54 = 226 Fp12; verify 1833); labels are 255-bit / four u64 limbs and the step
counter is 87-bit.

## 6. Scope and limitations (disclosed, not papered over)

* The fixture is 785 classes at theta = 2^5, so theta ~ sqrt(N) there; production
  is theta = 2^80 against sqrt(N) = 2^125.6. The fixture therefore measures the
  protocol's mechanism and the conditional detection rate, not the production
  birthday regime; the production success is the analytic law in section 5, and
  no production run is claimed.
* The fixture walk is degenerate (it collapses onto few classes, so merge rates
  are far above a uniform birthday rate). Every conditional number above is
  fixture-scope; the charged factors are the CP-95 lower bound and the analytic
  tail ceiling, not the fixture point estimates.
* No production mixing bound is claimed. mu = 11/10 is charged as a work
  allowance. The claim is behind the live board minimum and this account's other
  tracks; it is submitted because the repaired accounting is the point.
* `evidence/` was regenerated by `scripts/measure_v33.py` (all numbers) and
  `scripts/regen_v33.py --write` (ledger + certificate); `--check` re-emits and
  agrees, `--control-ledger` runs control (d).
