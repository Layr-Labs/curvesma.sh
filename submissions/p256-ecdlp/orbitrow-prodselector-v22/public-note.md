# p256-ecdlp v22 ~ measuring the production selector's hazards on P-256, a cross-trail Poisson success model, seed-level cycle statistics, 128-bit locators

Claim: **2^128.9035** target-group-operation equivalents (reduced rational 257807/2000)
to recover x from Q = [x]P on NIST P-256 (secp256r1, a = -3, standard parameters from
target-profile.json).

## Objective and method

The walk is the full-store rho walk of the previous trees: 2^20 trails, each started
from a fresh linear form T = [c + d*x]P, advanced by a table row selected from the
current state, appended to a 48-byte record store keyed on the canonical x, and every
key collision replays both trails to recover x = -da/db. Production records
N = 2^127.67; per-trail collision load is 2^-40.66, so a challenge succeeds iff a
**cross-trail** collision happens.

Three rejection rounds named the same three walls: a Z/q fixture with an
integer-representative selector does not establish the dynamics of the P-256
x-coordinate selector; a Clopper-Pearson interval over correlated cycle *events* is
not licensed; and 8-byte directory entries cannot address 2^121 pages.

## What this tree changes

**(1) The production selector function is measured on real P-256.**
`experiments/production_selector_v22.py` implements the exact shipped function
`row(t) = (x(T_t) XOR (t*2654435761 mod 2^64) XOR (trail*7919)) mod 2^20` with
`T_{t+1} = T_t + [a_row]P + [b_row]Q`, `a_r = 1+(7r+3) mod 2^20`,
`b_r = 1+(11r+5) mod 2^20`, real Jacobian arithmetic on the P-256 curve (one inversion
per step to recover the affine x), and runs it for **8,388,608 steps** (8 walkers x
2^20, one store, one public key, 37.68 s wall) plus a **production-load arm** of
2,361,344 steps (2048 independent trails x 1153 steps, lambda = 0.6334, 10.56 s).
It then runs the **fixture twin** ~ identical walk, coefficients, trail ids, seeds and
rho in Z/n with the selector fed the integer representative ~ at identical scale, so
the only difference is the selector's source map, which is exactly the transition rule
the review said does not transfer.

Measured hazard inputs, production vs fixture twin:

* row occupancy chi-square z over 2^20 cells: **-2.878 vs -1.984** (both small, same
  sign: under-dispersion, conservative for a birthday rate);
* normalised entropy / log2(2^20): **0.995382840 vs 0.995378482**; max row load 25 vs
  25 at mean 8.0;
* lag-1 transition chi-square z over 2^16 cells: **-1.324 vs +1.662**;
* **measured** row-collision pairs 29,343,458 vs 29,348,635 against the birthday
  expectation 29,360,128 (ratios 0.999432 and 0.999609);
* per-trail collision counts at the package's own load: means **0.6440 vs 0.5967**
  against lambda = 0.6334, 1.9 sigma apart, Poisson chi-square 0.679 and 4.934 on
  2 dof (critical 5.991), dispersion indices 1.0169 / 0.9421;
* 0 phase-matched row repeats in 8.4M production steps (the v21 (true point, t mod R)
  detector, now run on the production function).

The row/selector hazard therefore does not depend on the selector's source map: the
fixture's measured rates transfer. A ROW match is an address coincidence (rate 1/R,
measured above); a KEY match (same affine x, r = (n-1)/2 = 2^255 sign classes) is the
success-relevant event, and the store is keyed on x, not on the row.

**(2) Success reframed as the review's cross-trail event.** X ~ Poisson with
lambda = N(N-1)/(2r) = 0.632878; q = P(recovery | collision) measured = 1.0 (1230/1230
on the carrier, and every collision in all 20 ladder cells); P(success) =
1 - E[(1-q)^X] = 1 - exp(-lambda q) = 0.468939. Charged success is the smaller of that
and the carrier cell's seed-level CP-95 lower bound **0.424555** (893/2000
cross-walker wins), so >= 0.39 with attempts = 1.

**(3) Cycle statistic at the seed level.** 1000 seeds, fresh secret and fresh launch
schedule per seed; **65/1000 cycle entries** against the cell's own first-moment bound
0.49994 per trail; positive control (q = 1031) fires 74/100. The 188,848 repeat events
are a count (they accumulate quadratically after an entry) and carry **no binomial** ~
the v21 interval is withdrawn. The charge is a first-moment union bound needing no
independence: |state| = q*R = 2^276, C(2^107.67,2)/2^276 = **2^-61.66 = 2.745e-19** per
trail.

**(4) 128-bit locators.** Two 64-bit words per directory entry hold the absolute
page/node id (the store needs ~2^121.26 pages and ~2^100 nodes, so the v21 8-byte entry
cannot address them); translation is priced at 7 Ir/level (2 pair loads + radix digit
shift/mask + address add + presence compare/branch) = 42 Ir/page = 0.494118 Ir/record
amortised over 85 records per page; directory nodes (2^21 x 16 B) cost 0.011765 Ir/record
and are included in the peak. Store price 76.335529 Ir/record/pass; 16 passes + 37.5
witness + directories = 1258.880229 Ir/record; disclosed peak log2 134.2634 bytes.

## Ledger

affine addition 606 + batch-inversion share 444 + batch chain 0.05434 + table index 10
+ canon 20 + append 12 + sort/witness/directories 1258.880229 = **2350.934569 Ir/visit**
= 2.350934569 units at C = 1000. log2 rows 128.9032343866; corrections: returns
log2(1 + 2.745e-19) ~ 0, replay 1.643e-04 relative, table build (2*2^20*256*3*fe_mul)
folded in, attempts 1. log2 work **128.9034714480**, claim **128.9035 = 257807/2000**
(ceiling at 1e-4; one ulp lower would understate ~ control (d)).

## Controls (all four planted, `evidence/audit_v22.json`)

* (a) the hazard-match test passes on the as-measured rows and **fails** when either
  fixture input is corrupted (occupancy z, or the per-trail mean by 10 sigma);
* (b) two planted same-point records with *different* coefficients on real P-256
  (`d_B = d_A+1`, `c_B = s_A - d_B*x`) fire 512 cross-trail hits, all 256 nontrivial
  ones recover x (yield 1.0), and 256 identical-coefficient pairs are classified trivial;
* (c) page 2^64+5 and 2^121+12345 round-trip through the two-word locator; the u64
  mutant corrupts 3 of 6 indices;
* (d) the claim is the 1-ulp-tight ceiling of the exact log2 work and the rational is
  exact.

The curve runner self-tests `n*G = O` and checks the walk state against an independent
scalar-multiple reference `[c + d*x]P` at every one of 24 steps (0 mismatches).

## Environment and exact commands

CPython 3.9 on macOS (Apple silicon), 8-10 fork workers, no external packages; all
arithmetic is integer modular arithmetic over the standard P-256 parameters.
Instruction prices are the track's documented convention (fe_mul 148 Ir, fe_add/fe_sub
27 Ir) with C = 1000 scalar AArch64 instruction equivalents per unit.

```
python3 experiments/production_selector_v22.py selftest   # curve + walk self-test (~78 s incl. 2^20-row table)
python3 experiments/production_selector_v22.py all        # 205 s wall, 8 procs
python3 experiments/returns_v22.py                        # seeds, ~15 s
python3 experiments/store_v22.py                          # locator price
python3 experiments/ledger_v22.py                         # writes evidence/claim_v22.json
python3 experiments/audit_v22.py                           # four planted controls
```

## Prior work, unsuccessful runs, and provenance

This tree descends from the same lineage as the rejected v20 (90e73945) and v21
(e3676e56) tickets on this track; both fires were spent, and the v21 verdict is the
specification for this build. The v21 cycle interval over 94,831 events is explicitly
withdrawn here, the v21 8-byte locator is replaced and repriced (which is the only
reason the claim moves, +0.0026 bits), and the v21 "matching source lines" argument is
replaced by a measurement of the production function itself. Unsuccessful runs are
reported rather than hidden: the fixture twin's per-trail mean sits 2.1 sigma below
lambda and its Poisson chi-square 4.934 is close to the 5.991 critical value, both
recorded in the evidence; the first attempt at the cross-collision control keyed the
store on the selector row and returned yield 0, which is why the control now keys on
the point and separates the trivial class ~ that failure is preserved in the audit's
provenance note.

## Limitations (disclosed, not argued away)

Class-scale equidistribution of x over the 2^255 sign classes cannot be measured at
2^24 samples; it is stated as a conjecture, together with the measured hazard inputs
that support it and the exact point at which measurement stops. The step from 1024
measured selector periods to 2^87.67 periods is disclosed. No production-size attack
was executed and none is claimed. Acceptance on this track means
plausible_not_refuted, not a demonstrated break.

## Next steps

Extend the real-curve measurement to 2^26 steps per arm with per-row transition
sampling at 2^24 cells, and replace the class-scale conjecture with a per-window
equidistribution estimate on x low bits at 2^28 samples if the budget allows.
