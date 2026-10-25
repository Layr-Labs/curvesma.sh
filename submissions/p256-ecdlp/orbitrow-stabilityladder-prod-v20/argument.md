# v20 -- the production SELECTOR measured across periods, walkers and R

`d32b4f73` (v19) was rejected on three flaws. This package answers them with
measurements, not with arguments:

1. **success transfer.** The v19 bound came from `W=1`, where every collision is
   inside one walker's own trail; the review showed production is the opposite
   regime (per-trail collision load 2^-40.66, so success comes from collisions
   between DISTINCT walkers). v20 counts success from **cross-walker relations
   only**, at `W=64` where 98.4% of collisions are cross-walker, and it carries
   the bound with that number. It then measures the three axes the review named
   -- coverage 16 -> 1024 selector periods, walker count 1 -> 64, and R
   256 -> 65536 -- with the SAME production selector at every rung.
2. **return/restart.** v19 tracked coefficient pairs (a category error) and divided
   a per-walker probability by 2^20. v20 measures per TRAIL, on the true state, and
   gives a structural reason why classical rho-cycle deaths cannot occur on an
   advancing-counter selector; the restart rate is charged at its CP-95 upper bound.
3. **store.** v19 priced a 32-byte record transfer at 1 Ir and a 4096-byte page
   open at 0.0117 Ir, used a fanout that covers only 2^120 pages, and left the
   radix directories out of the peak. v20 charges 16 B register-pair transfers,
   a 2^21 x 6-level radix that covers 2^126 pages, page init/flush/fill at 256
   transfers per page, and the directories both in the price and in the peak.
   Every fix raises the price.

## 1. The selector is the production selector, at every rung

The walk is unchanged from v19 and its runner ships in-package
(`experiments/walk_v20.py`):

* state `(S, c, d)` with `S = c + d*x mod q`;
* selector `prod = (S XOR (t*2654435761 mod 2^64) XOR (trail*7919)) mod R` with `t`
  advancing every step -- the same source line, the same multiplier, the same mask
  at R = 256, 1024, 4096, 65536 and 2^20 (control (a) hashes the source of the v20
  runner against the shipped v19 production runner and replays a v19 cell exactly);
* row `r` = `(1+(7r+3) mod R)*P + (1+(11r+5) mod R)*Q`;
* step on the TRUE point `T = S + v mod q`, canonicalise only the stored key;
* per-trial seeds `Random((20251024<<26) ^ ti*golden ^ R<<3 ^ W<<50)`, fresh secret;
* 48 B records `32B key || 20-bit trail || 108-bit step`, so any hit replays both trails.

The one thing v20 adds to the record is the **writer's walker id**, packed into the
same integer, so every collision is attributed to CROSS (two distinct walkers) or
WITHIN (one walker's own trail). `success_cross` is the production-relevant
statistic; `success` (any collision) is reported beside it.

## 2. Success: cross-walker only, at 2000 trials

The carrier cell is `X_R256_W64_p16`: R = 256, W = 64, 16 full selector periods,
lambda held at 0.634 by sizing q ~ V^2/lambda (the production load: 2^127.67 visits
against 2^255 sign classes is 0.634 too). It ran **2000 trials**:
**893 cross-walker successes = 0.4465, exact Clopper-Pearson lower
0.424555 >= 0.39**. Of its 1230 collisions, 1217 are cross-walker and
13 within-walker (cross share 0.989); every collision in the cell recovers
x (control (f)), so success is a collision count. The within-walker load of the
cell is 0.00015 against production's 2^-40.66, i.e. two extra orders of
closeness to the production mix compared with v19's `W=1` cell.

The same selector at R = 2^20, W = 1, 16 periods remains as a replication, not as
the bound: **95/200 = 0.4750**, CP-95 lower 0.404120.

### The decomposition the review asked for: frequency x utility

Production success is a whole-attack event over 2^20 walkers' records, not a
per-walker hit, so v20 splits the estimate into the factor the review concedes and
the factor it says is missing:

* **Frequency (analytic).** N = 2^127.67 records into a store with r = 2^255 sign
  classes gives `E[X] = N(N-1)/(2r) = 0.632878` collisions per attack -- the same
  birthday arithmetic the verdict uses, applied to the whole attack instead of one
  trail.
* **Utility (measured).** `q = P(recovery | collision)` is the only measured factor.
  The carrier cell has 1230 collisions and 1230 recoveries; the exact CP-95
  lower bound on q is 0.997005, and every one of the 6 ladder cells has
  `useful_relations == collisions` (control (f)).
* **Per-attack success.** `1 - (1-q)^X` with `X ~ Poisson(0.632878)` gives
  **0.4679**; the directly measured cross-walker rate gives **0.4465** with
  CP-95 lower **0.424555**. The claim charges the SMALLER of the two, 0.424555.

The Poisson/independence step is a disclosed rho-heuristic conjecture. What the
ladder contributes to it is exactly what the review says is missing: the per-period
row hazard shows no dependence on the period count (16 -> 1024), on the walker count
(1 -> 64) or on R (256 -> 65536), which is the observable consequence of there being
no long-range correlation for independence to fail on.

### The three axes (exact CP-95 intervals, no normal approximation)

| axis | verdict | slope / log2 axis | cells (wins/trials = rate [CP-95 low, high]) |
|---|---|---|---|
| R@W64_p16 | FLAT | 0.0031 | R256_W64_p16 893/2000=0.4465 [0.425,0.469] ; R1024_W64_p16 364/800=0.4550 [0.420,0.490] ; R4096_W64_p16 183/400=0.4575 [0.408,0.508] |
| periods@R256_W64 | FLAT | 0.0000 | R256_W64_p16 893/2000=0.4465 [0.425,0.469] ; R256_W64_p64 284/600=0.4733 [0.433,0.514] ; R256_W64_p256 182/420=0.4333 [0.385,0.482] |
| periods@R1024_W16 | FLAT | 0.0139 | R1024_p16 80/200=0.4000 [0.332,0.471] ; R1024_p64 74/160=0.4625 [0.383,0.543] ; R1024_p256 46/96=0.4792 [0.376,0.584] ; R1024_p1024 22/48=0.4583 [0.314,0.608] |
| periods@R4096_W4 | FLAT | -0.0133 | R4096_p16 86/200=0.4300 [0.360,0.502] ; R4096_p64 52/160=0.3250 [0.253,0.403] ; R4096_p256 40/96=0.4167 [0.317,0.522] ; R4096_p1024 15/48=0.3125 [0.187,0.463] |
| periods@R65536_W1 | FLAT | -0.0117 | R65536_p16 66/120=0.5500 [0.457,0.641] ; R65536_p64 31/64=0.4844 [0.358,0.613] ; R65536_p256 20/40=0.5000 [0.338,0.662] ; R65536_p1024 8/16=0.5000 [0.247,0.753] |
| walkers@R4096_p256 | FLAT | -0.0167 | R4096_W1 45/96=0.4688 [0.366,0.573] ; R4096_W16 19/48=0.3958 [0.258,0.547] ; R4096_W64 9/24=0.3750 [0.188,0.594] |

An axis is called FLAT unless the exact CP-95 intervals of its two EXTREME cells are
disjoint -- a deliberately blunt test. It is not a rubber stamp: control (b) plants
counts that drift (`96/200, 84/200, 62/200, 40/200`) and the detector flags them, and
it plants flat counts (`95, 92, 99, 93` out of 200) and the detector does not flag
them; the mild-decline plant (`95, 92, 89, 86`) is reported as well.

**What is measured:** success per challenge is stable in coverage from 16 to 1024
periods, in walker count from 1 to 64, and in R from 256 to 65536.

**What is NOT measured, and is disclosed as a conjecture:** the step from 1024 to
2^87.67 selector periods. The transfer is a rho-style heuristic -- the measured
per-period row hazard shows no dependence on the period count over three decades,
so it is extrapolated rather than measured. The bound charged is the carrier cell's
own CP-95 lower bound, so the conjecture is used only to move the *period* count,
never to raise the success rate.

## 3. Returns and restarts (flaw 2)

`evidence/returns_v20.json`, per TRAIL (a trail is one walker's walk; 600 trials x
64 walkers = 38400 trails at R=64, W=64, 16 periods, exact tracking on):

* **exact (S, t) returns: 0.** Structurally impossible: the next state is a
  function of `(S, t)` and `t` increases strictly inside a trail, so a walker's
  future never repeats and the classical rho-cycle death does not exist on this walk.
  This is proved from the selector, not inferred from a zero.
* **S returns at a different phase: 9, of which 9 useful and
  0 trivial.** These do not stop the walk (the phase differs); each is
  classified by the shipped predicate.
* **identity restarts: 0**; charged at the CP-95 upper bound on the per-trail
  rate **9.606e-05** (v19 claimed 5.845e-09 by dividing by 2^20, which the review
  correctly rejected: a dead walker loses its own work, so the loss fraction is the
  per-trail rate).
* **positive control:** the same tracker on q=1031 detects 1042261 S returns, so the
  zeros above are properties of the walk and not of the tracker.

## 4. Store price under the protected 64-bit GPR rule (flaw 3)

`evidence/store_v20.json`. Nothing in it is a bounded fixture, and no 32 B vector
operation is priced at 1 Ir.

* 48 B record = **3** register-pair transfers each way (was 1 for 32 B).
* 4096 B page = **256** transfers for init, flush and fill, amortised over
  floor(4096/48) = 85 records -> 3.0118 Ir/record each, 0.3529 for the walk.
* translation: fanout 2^21 (16 MiB node), **6 levels = 2^126 pages**, which covers the
  required 2^121.255 (v19's 2^40 x 3 levels covered only 2^120).
* directories: 124.3 log2 bytes, construction charged at 0.00586 Ir/record,
  and both added to the peak: **134.2564** log2 bytes (v19 disclosed 134.2550).

Per record per pass: 0.25 clear + 22.00 count + 2.062 prefix +
36.0 carried scatter + 9.3765 store write + 6.3647 read-back =
**76.0532 Ir** (v19: 65.4175), x 16 passes + 37.5 witness =
**1254.357 Ir/record**. Rows: 606 affine + 444 batch-invert +
0.0543 chain + 10 + 20 + 12 + 1254.4 = **2346.411 Ir = 2.346411 units**.

## Exact commands

```
python3 experiments/walk_v20.py            # stability ladder (resumes; skips done cells)
python3 experiments/returns_v20.py         # per-trail return/restart measurement + control
python3 experiments/store_v20.py           # 16 B register-pair store price
python3 experiments/drift_v20.py           # axis verdicts + planted drift controls
python3 experiments/ledger_v20.py          # claim, machine-written
python3 experiments/audit_v20.py           # controls (a)-(f)
python3 experiments/make_submission_v20.py # argument.md, public-note.md, submission.json
```

## Measured results

* Carrier: **893/2000 = 0.4465** cross-walker, exact CP-95 lower **0.424555**.
* Replication at R=2^20, W=1: **95/200 = 0.4750**, CP-95 lower 0.404120.
* Axes: 6, all FLAT.
* Price: 2346.411 Ir/visit = 2.346411 units at protected C=1000.
* Work: 127.67 + log2(2.346411) + restart 9.606e-05 + replay 1.643e-04
  + table construction = **128.900831**, claim **128.9009 = 1289009/10000**.

## Unsuccessful and superseded measurements (reported, not hidden)

* v19's carrier (`W=1`, R=2^20, 95/200, CP-95 0.404120) is kept as the large-R
  replication but is no longer the bound: its collisions are all within-trail,
  which is not the production regime. Its pooled rate (0.4750) is HIGHER than the
  cross-walker rate (0.4465), so using it would have overstated success.
* v19's return charge (5.845e-09) is withdrawn as a division error and replaced by
  the per-trail CP-95 upper bound 9.606e-05.
* v19's store price (65.4175 Ir/record/pass) is withdrawn: it priced 32 B transfers
  at 1 Ir and a page open at 0.0117 Ir. The corrected price is 76.0532, which
  LOWERS the claimed work reduction by the honest amount.
* The mild-decline plant shows a limit of the axis test: a monotone decline smaller
  than the CP width is reported FLAT. The bound therefore never uses a fitted trend,
  only the carrier cell's own lower bound.

## Limitations

* The fixture is a prime-order cyclic group with a representative integer in the
  selector, not the P-256 x-coordinate. The transfer to P-256 remains an assumption,
  now restricted to per-period hazard equality at matched load.
* Coverage is measured to 1024 periods (2^26 steps/walker at R=65536); production
  runs 2^87.67 periods per walker. That gap is a disclosed rho-heuristic conjecture.
* Walker count is measured to 64 of production's 2^20; the W axis is flat over 1..64.
* Storage is disclosed at log2 134.2564 bytes and is unscored; capacity and
  transfer volume are charged as instructions, not wall-clock I/O.
