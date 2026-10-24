# Argument — plain-orbitrow-40pct-v20d

## Claim

On `bls12-381-g1-ecdlp-plain` (input `T = [x]P0`, `x` uniform in `[0, r)`, no
auxiliary powers), the algorithm specified here recovers and verifies `x` with
probability at least 0.39 using at most **2^126.262** reference additions at `C =
2000`, i.e. `W = (1 + 161/2000) * (1.015^2 S0^2 + 2^60) + EXTRAS` with `log2 W =
126.26195600`. The success claim is carried by a 192,000-trial census whose minimum
exact CP-95 lower bound over the declared band d = 5,7,9,11 is **0.398043**, and by
the ledger's model success 0.399004.

## 1. The step and the price (v20c, unchanged)

The step is the 448-bit-window twist index of v19d: `h = argmin` over `x, beta x,
beta^2 x` from the low 448-bit window of `x * floor(2^448 beta/p)` (12 partial
products with `i+j in {5,6}`), `K = min(y, p-y)`, `m = (K*R) >> (bitlen(p)-1)`,
`S <- S + Tab[h][m]`, `(a,c) <- (a + r[h][m], c + s[h][m])`. The ledger prices 161
itemised instructions on top of one reference addition: index 46, v0 10, v2 6, argmin
14, sign/key 17, slot/address 7, loads/stores 28, coefficient update 12, checkpoint
12, inversion share 1, counters 7, meter 1. `drain = workers*2*DP = 2^81` steps and
the one-spacing `+2^60` lag term are inside `W`; the near-tie fallback is an exact
union bound at `3*2^-64` per step and is 2^-64.9 of phase 1.

Nothing in this package changed the price. `ledger.py` is v20c's file byte for byte
in its arithmetic: `log2 S0 = 126.18...`, `log2 W = 126.26195600`, claim 126.262,
ship line 126.266, gate 126.274 -> margins +0.004 and +0.012. The kappa product is
1.096914315 <= 2^(claim - log2 S0) = 1.096947771.

## 2. The mechanism the census replays

`qwalk.Universe` is the quotient walk of the mu_6 action on `Z_q` (`q = 1 mod 6`,
`lambda^3 = 1`): `class(c) = min{+-lambda^k c}` and each step advances the class by
`u[m]*x` with the coefficient transform `(a,c) <- e(a, c+u[m])`, `e = canon/(canon+J)
in mu_6`. The invariant

    a + c*x == canon  (mod q)

holds at every step by construction, so two records in one class from two different
walkers solve `x = (a1-a2)/(c2-c1) mod q` exactly and are verified against the
trial's `x`. This is the surrogate for `S = [a]G + [c]T` with `T = [x]G`.

Protocol (van Oorschot-Wiener, what the ledger prices): `theta = 2^-d` gates
**reporting only**; a walker that hits a DP appends `(class, a, c)` to the pool and
restarts at a fresh `[a0]G + [c0]T`; total steps per trial = `kappa*S_ref`; pool =
records only, `records = budget/2^d`; detection = two pool records of one trial in
one class from two different walkers; recovery verified. The independent merge audit
is over in-budget visited classes (any levels, two different walkers), so
`success <= merge` is a genuine invariant, not a restatement.

## 3. The drain arm (the v20d change)

v20c stopped each lane at the budget. A merge in the last DP spacing of the budget
then produced no report, and the census read 0.360 against merge 0.418 at d=5. In
production the lost window is `2*DP/S_ops = 2^-65` (per walker, `2*DP` out of
`2^106.3` steps); the toy's budget is 0.18–11.7 spacings per walker, so the same lag
costs it 1–12% of its merges. The drain arm closes exactly that gap and nothing more:

1. After the budget, a walker whose trial has an audited in-budget merge but no
   in-budget duplicate report is stepped to its next DP; the report enters the pool.
2. The flush cap is 8 spacings per lane (charged drain: 2 spacings). A lane that has
   not reported is dropped, so no lane ever walks further than the flush bound.
3. A duplicate pair counts only on a trial with an audited in-budget merge.
   Drain-only collisions are counted (`drain_only_excluded`) and excluded: at d=11
   the toy's drain is ~11x its budget, and counting them would give success ~1.0 by
   manufacturing search that production does not pay for.

`success_budget` (0.3641 / 0.2422 / 0.0968 / 0.0296) and `success_drain` (0.0475 /
0.1629 / 0.3054 / 0.3748) are reported separately, so the reader can see that the
in-budget statistic is the v20c one and the drain adds only the lag.

## 4. Census

48,000 trials per density at d = 5, 7, 9, 11 (192,000 total), seed 20270 +
d*1000 + shard, 64 walkers per trial, `kappa = 1.015`, `N = 561,690,001` classes
(`2^29.07`).

| d | success | CP-95 exact | merge | records | expected | drain pairs excluded |
|---:|---:|---:|---:|---:|---:|---:|
| 5 | 0.402313 | 0.398625 | 0.4130 | 745.7 | 746.0 | 451 |
| 7 | 0.402646 | 0.398958 | 0.4108 | 186.4 | 186.5 | 118 |
| 9 | 0.401729 | 0.398043 | 0.4032 | 46.6 | 46.6 | 21 |
| 11 | 0.404250 | 0.400560 | 0.4042 | 11.6 | 11.7 | 5 |

CP-95 is the exact one-sided Clopper-Pearson lower bound in Decimal (the
`gen_cert_v15k.py` function, with the bracket opened to [0,1] because that census's
[0.30,0.45] is not a general bracket; the arithmetic is identical). The certificate
recomputes it from `(k,n)` and compares to the stored value at 1e-12, then checks it
against 0.39 at every density. Zero successes without a merge, zero solve mismatches,
zero zero-denominator solutions; 145,798 verified recoveries across 77,325 successful
trials (the excess is multiple verified duplicate pairs per trial).

The merge rate exceeds the 0.399 model by 0.3–1.3%, the same excess v20c measured
(0.418 vs 0.399 at 400 trials) and in the direction the audit predicts: it sees the
reported class as well as the visited ones, so it is slightly above a pure coverage
birthday.

**Low-n disclosure.** At 8,000 trials/density the same arm read CP-95
0.3921/0.3938/0.3872/0.3908 — d=9 below the floor by 0.0028. That file
(`evidence/prod-sweep-8000.json`) is shipped and is not the gate. At p = 0.399 the
decision rule needs n >= 8,012 for CP-95 >= 0.39, so 8,000 is one trial short of the
power it needs, and the declared census is 48,000 at all four densities.

## 5. Reproduce

    python3 ledger.py
    python3 scripts/gen_cert_v20d.py
    nice -n 10 python3 scripts/prod_sweep.py --trials 48000 --chunk 100 \
        --dens 5,7,9,11 --workers 8

The census is 2m35s wall on 8 workers at 732% CPU. Raw shard files are not shipped
(1,920 files would breach the 30-file archive); the reduced evidence carries every
per-density count the certificate needs, and the command regenerates the shards.
