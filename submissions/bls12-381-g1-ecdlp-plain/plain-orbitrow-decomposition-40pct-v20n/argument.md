# Argument - plain G1 v20n

## Mechanism

The walk is the mu_6 quotient walk: class(c) = min{+-c, +-lambda c, +-lambda^2 c}, one
seed u[m] per column, and the step m = col(canon), t = canon + u[m]*x, canon' = canon(t)
with (a, c) <- e*(a, c + u[m]) so a + c*x == canon (mod q) holds at every step. Two DP
records of one class from two different walkers solve x = (a1 - a2)/(c2 - c1) exactly.
The two-step repair substitutes the alternate orbit row when the successor class equals
the class one step back; the row is orbit-generated, so the repair is equivariant under
mu_6 and the escape decision is a class function.

## Evidence shipped in this revision

1. `evidence/absorb-direct.json` - the repaired predecessor-dependent process. prev is the
   step's RETURNED predecessor every step, kcap = 8 with a coin generator, and the
   repetition key is the full state (class, prev, escape-run). 16,384 starts at each of two
   sizes; period_2 = 0 at both; the three-cycle share of the residual falls 0.0263 ->
   0.0007 over a 4.0x class increase. `evidence/absorb-direct.json:mutation` is the re-pin
   control (the v20l defect) and moves the counts by 0.7810 / 0.4626.
2. `evidence/restart-count.json` - written by `scripts/restart_count.py` from itemised
   383-instruction field products: 2,681 insn per doubling (2M+5S), 4,213 insn per window
   addition (11M), 953,287 insn per 256-bit scalar action (476.644 ref), 2,970,320 insn per
   fresh start (1,485.1600 ref). `ledger.py` loads this file; the self-audit binds
   `ledger.FRESH_REF` to `evidence/restart-count.json:fresh_ref`.
3. `evidence/selfaudit.json` - the citation audit: every cited file must exist, every
   `file:field` citation must resolve, every evidence file must be cited, and a planted
   wrong path must be caught. This is the check whose absence let a write-path typo hide
   an entire absorption experiment for a round.
4. `evidence/walk-census-exact.json` and `evidence/walk-census.json` - the at-budget
   census (visits/sqrt(N) asserted within 1e-3) and the drain=0 production-shape arm whose
   CP-95 (0.376122) minus debits is the measured success 0.373911.

## What is measured, what is modelled, what is missing

Measured: absorption and its cycle-length split, the mutation control, the restart price,
the at-budget success rows. Modelled: the 0.3995 production decomposition - its (w-1)/w
walker factor and DP detection window are analytic at 2^20 walkers and 2^-60 spacing.
Missing, and stated as missing: the trajectory-matched d=0 control. `scripts/census_v20n.py`
shares one fixed start block across every arm and compares the escape and no-escape
trajectories byte-for-byte up to the branch step, but it aborted on its own assertion that
at least one trial fires the escape: zero fires in 86.4M lane-steps at the production-shaped
fixture against an analytic prediction of order 4e5. The absorption harness on the same step
fires at the model rate, so the divergence is in the invariant-consistent start block, not in
the repair. No proxy number is substituted for the missing control.

## The cap

The deterministic work meter bounds EXPENDITURE (METER_CAP 8.0021e31 ref). It is not a
lossless cap and no success-loss bound is derived from it; the success row carries no
cap-loss debit and the note says so.

## Position against the bar

log2 W = 126.26594404, claim 126.266, kappa product 1.099950711 <= 1.099993376. The claim
is 0.0248 above the 126.2412 frontier and the measured success is 0.0161 below the 0.39
floor; `ledger.py` prints FLOOR_MET=False and DECISION HOLD. The revision is an evidence
integrity revision, not a claim increase.
