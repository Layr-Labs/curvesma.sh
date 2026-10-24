# argument (v32) -- repairs verdict b30ddc21 (rejection of gt-instantiable-40pct-v31)

Claim **127.7213**. Every number below is emitted by `scripts/gen32.py` from the same evidence the
certificate audits; nothing is asserted that the generator does not recompute.

## 1. Capacity is derived at the production steady state, not transferred from the fixture

The occupancy law is `occupancy = min(capacity, arrival_rate x mean_residence)` with arrivals at one
insert per step, so `mean_residence = capacity` for an evicting store and `occupancy = min(capacity,
arrivals)`. Two routes are computed.

**Route B (v31's 2^21 cap) -- computed and REJECTED.** buckets 524288, capacity 2097152, occupancy 2097152 =
capacity: the store is SATURATED. Residence is capacity/rate = 2097152 steps. In a birthday-limited rho
search the number of candidate pairs at gap g is (S-g)/N, so the gap to a matching counterpart is
uniform on [1,S]; the probability a record survives that gap is (M/S)(1-e^-S/M) = 2.189e-32 = 2^-105.17.
Success would be 1.262e-32 -- four orders below the 0.39 floor. The v31 census could not see this because
its fixtures run at 1164/16384 occupancy, i.e. far UNDER capacity, so their 0.08% eviction loss is
not the production ratio.

**Route A (shipped here).** Capacity is sized from the arrivals the run actually makes: buckets =
next_pow2(2S) = 340282366920938463463374607431768211456, capacity = 4 x buckets = 1361129467683753853853498429727072845824 records, occupancy = S = 95803668470929860550224427498997735090. No record is ever
evicted, so the retention factor is 1 **by construction**, and the success resampling drops the
retention factor instead of importing a one-step, depth-5 measurement (the exact transfer the
verdict refused).

`trap_waste` and `store_init_route_a` are charged once, at production scale; the storage the route
needs is 194705035654402988429786479712611689348352 B (log2 137.16). That is the machine the claimed 142.01-bit computation requires; it is
stated plainly, charged for its one-time initialisation, and reported as an unweighted memory row
(the cost model scores computation, not memory).

## 2. Every auxiliary array is bounded

`bounds_table` in `evidence/ledger32.json` lists nine arrays with a bound formula from design
parameters and a byte value: records/alive/matched bounded by min(arrivals, capacity) = S, guard ring
by depth 6, probe chain by PROBES 4, export scratch by 10 words, sort work by 12x6 words, node pool
by 2^130, bucket array by 2^128, counters by 20. Nothing is left "grows with arrivals, unbounded".

## 3. The full-width state-key path is priced, at its EXECUTED maximum

v31 priced `state_key48w` (3398 Ir) and nothing before it, while `Group.state` first computes the
twelve canonical 381-bit residues, their min/argmin and a sort. `scripts/sort32.py` executes that
path 20000 times at production width with counting instrumentation on every comparison and move:
comparisons p50/p95/max = 53/64/75, executed Ir p50/p95/**max = 1095/1263/1458**. The ledger charges the
MAX (1458 Ir/step), not an aggregate divided by calls -- the verdict's objection to `exec31`'s
mean-per-call. `exec31`'s own charge rule (static span + per-loop structural maximum x body span)
upper-bounds every individual call, so it bounds the maximum a fortiori; the divided figure is
retained only as a labelled mean.

The coefficient-update row is now an implementation bound: `A += linv[kap]*a` and `D += linv[kap]*b`
are two 255-bit modular products at 64 Ir each, plus the 12-lane load -- 176 Ir/step, not 48.

## 4. Trap charge is a Markov union, with no independence premise

A trapping window requires a guard admission (the abort is armed by one and needs K_STALL = 32
squarings after it). The census counts that necessary precondition exactly: 4497 admissions in
5040000 steps. So the number of trapping steps at S is at most ceil(4497*S/5040000) = 85481963713049917240944295726784289 and
P(any trap) <= that count -- a union bound over steps, no independence, no Poisson limit, no
asserted overlap multiplier, no rule of three on a zero count. The 0-admission window census
(lengths 14..22, 360000 windows) is supporting evidence only. Charge = 85481963713049917240944295726784289 x 32 x STEP = 131037341647165465940094003989067234091060 Ir.

## 5. Exact integers

`trap_at_S` = 85481963713049917240944295726784289 and `trap_waste` = 131037341647165465940094003989067234091060 are Python ints (asserted by type in the certificate); v31
computed trap_at_S through float arithmetic and used floor division in trap_waste. Every row is an
exact integer ceil; the claim is the minimal n with 20000*2^(n/10000) >= W, verified in both
directions on the integers.

## 6. Certificate regenerated, tree swept

`certificate.json` is emitted by `scripts/cert32.py` from THIS ledger (v31 shipped a certificate for
the older 127.6293 claim). The same script sweeps the tree for literals that contradict the current
values (127.6293, 1276422, "524288 buckets", "2097152 records", "overlap factor of six") and reports
the sweep verdict inside the certificate.

## Success

Nonparametric bootstrap over the recorded fixture indicators (225/320, 200/320) with retention
removed: p5 = 0.578125, median = 0.625, floor 0.39 met. The v31 retention-resampled p5 (0.5767) and the Beta
model quantile (0.5762) are labelled comparisons.

## Rider from the same verdict, handled

`dependencies.json` stays empty and the "already accepted" collision result is no longer invoked
anywhere: the argument is self-contained, so nothing is claimed whose hypotheses are not reproduced
here.
