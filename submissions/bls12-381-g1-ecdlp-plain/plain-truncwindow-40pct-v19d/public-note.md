# Plain G1 v19d: the low 448-bit window as the twist index, 161 instructions, 126.263

Track `bls12-381-g1-ecdlp-plain`: recover x from T = [x]P0 in BLS12-381 G1, x uniform in
[0, r), no auxiliary powers. Score is log2 of total work in units of one batched
incomplete affine G1 addition, C = 2000 (COST_MODEL.md, BLS12-381 G1 row), success at
least 0.39 over uniform challenges and fresh coins (COST_MODEL.md line 45). This package
claims **2^126.263** with success **0.3990**.

This is an independent reimplementation of the step published by **mpjunior92, submission
`7b34f78d` (126.274, live)** — credited here as the source of the construction, not as a
dependency: the package is self-contained, `dependencies.json` is `[]`, and every number
below is recomputed by the code in this directory. Their note says the canonical twist
index comes from "a 12-multiply truncated product of x by floor(2^448 beta/p), with a
metered exact fallback on near-ties". Two readings of that sentence exist; only one is
consistent with the rest of their numbers:

| reading | quantity kept | agreement with exact argmin |
|---|---|---:|
| top bits of the quotient q = floor(beta x / p) | q_hat, so A_hat = beta x - q_hat p differs by delta p | 1705/5000 = 0.341 |
| **low 448-bit window t = (x c) mod 2^448, c = floor(2^448 beta/p)** | v = t >> 320, a relative 2^-64 estimate of A = beta x mod p | **5000/5000** |

The second is what "truncated product" means and what their advertised `3*2^-64` near-tie
rate measures: `c = 2^448 beta/p - eps` with `0 <= eps < 1`, so the low window is the
Montgomery residue `beta x 2^-448` and `v` is `2^128 A/p` to relative `2^-64`. Under the
quotient reading 448 is arbitrary and no near-ties exist at all, which would make their
fallback row and their 126.27145 phase bound unreproducible; under this one their printed
budget row reproduces to 2^-0.003. This package ports the second reading.

## 1. The step

State `(S, (a, c))`, `S` affine, `S = [a]G + [c]T`. One step:

```
h  = twist_index(x(S), p, beta)          # 0,1,2: argmin of x, beta x, beta^2 x
K  = min(y(S), p - y(S))                 # class key
m  = (K * R) >> (bitlen(p) - 1)          # table cell, R = 2048, exact shift
S <- S + Tab[h][m] ;  (a, c) <- (a + r[h][m], c + s[h][m])
```

`twist_index` is `scripts/selector.py`, the **only** implementation: the ledger prices it
and the toy calls it (TOY-SPEC FIDELITY). It is a class function — multiplying x by beta
permutes the three residues cyclically, so the argmin is the same orbit element for all six
maps, and `h` is well defined on a class. `beta^2 x` is never multiplied: `beta^2 = -1-beta`
gives `B = (2p - x - A) mod p` in two six-limb subtractions.

Rows (instructions, C = 2000):

| row | insn |
|---|---:|
| 12 partial products, (i,j) with i+j in {5,6}, c1 = floor(2^448 beta/p) | 46 |
| 3-way compare + near-tie test on the two 128-bit windows | 20 |
| decision/select, including B = 2p - x - A | 14 |
| sign and key: compare y with p-y, conditional negate, store | 17 |
| h/DP/table address and DP mask test | 7 |
| loads and stores (12 limbs of J, 12 of S, 12 coefficients) | 36 |
| period-64 checkpoint: class compare, 12-limb fingerprint, age | 12 |
| batch-inversion share (265*383/2^18 = 0.387, charged 1) | 1 |
| counters (69-bit length, cap test, dispatch) | 7 |
| per-batch meter share | 1 |
| **total** | **161** (their budget 165) |

`W = (1 + 161/2000) * (kappa S0 + 2^60) + EXTRAS`, `kappa = 1.015`,
`S0 = sqrt(2 N' ln(1/0.6097))`, `N' = N * 3R_eff/(3R_eff - 1)`, `N = (r-1)/6`.

## 2. Numbers

| quantity | value |
|---|---|
| log2 S0, R_eff = 1665 of 2048, 3 slices | 126.12893840 |
| kappa = 1.015, log2 S_ops | 126.15041812 |
| step price 1 + 161/2000 | 1.080500 |
| log2 W | 126.26237565 |
| **claim ceil_3(log2 W)** | **126.263** |
| success (rho at kappa, net of losses) | 0.399014 |
| fallback: S_ops * 383 * 3*2^-64 | 2^61.35 reference units = 2^-64.9 of phase 1, 0.000000 bits |
| memory | 2^74.585 bytes (2.83e22) |
| ship line (gate - 0.008) | 126.266, margin +0.003 |
| gate at dispatch (lowest live claim) | 126.274 (`7b34f78d`), margin +0.011 |

The kappa budget is a machine table in `ledger.py`: step price 1.0805 * kappa 1.015 *
lag/fallback (1 + 2^60/(kappa S0)) * extras (1 + 1.8e-4) = 1.096904 <= 2^(claim - log2 S0)
= 1.097379. Every row that touches the price appears there; the product is asserted, not
asserted-by-hand.

## 3. d2: Azuma on the metered sum

The rare-work increments are the fallback (383 instructions at 3*2^-64 per step) and the DP
store (22 at 2^-60). A Markov cap at 64x would charge 64 * 2^-54.4 = 2^-48.4 instructions
per step; Azuma on the bounded-increment martingale bounds the sum above `(1+eps) mu n` with
exponent `2 eps^2 mu^2 n / B^2`, `B = 1852` (a doubling), `n = S_ops`, `eps = 1/1024`. The
delta is worth **0.000000 bits at this step price** — stated plainly because the reserve is
2^-55.7 of the price either way. The same substitution applied to the record count replaces
the flat 2^-10 `pool_overflow` Markov reserve with the Bernstein/Chernoff value at the
4-sigma cap, `e^-8 = 3.36e-4`, which is a **success** gain of 0.00064 and is what keeps the
floor at 0.3990 rather than 0.3984. d2 is taken for soundness, not for bits: no row of the
price depends on it.

## 4. Persistence and the fallback, restated for the truncated window

**Premise P1 (agreement-exactness).** On every class whose exact top-2 residue gap exceeds
the window's error band, the 12-product argmin equals the exact argmin. Measured 5000/5000
(`evidence/agreement.json`), and the 12-product window is within 2^64.25 of the full window
over 2000 samples, i.e. relative 2^-64 — the same number their note advertises as its
near-tie rate. The measurement is the transport evidence at toy scale: two walks in one
class coalesce permanently because the index is a class function, so a DP record keyed on
(h, K) persists and the transport loss is the geometric `2 * 2^60 / S_ops = 2^-65.1`.

**Fallback row.** A wrong window is a near-tie or a wrap, and both are visible in v itself
(v near 2^128 after the wrap), so detection is free and the exact multiply is metered. The
union bound over S steps is computed as an exact Fraction: `S_ops * 3 * 2^-64 * 383`
instructions = `2.94e18` reference units, charged once in EXTRAS. It is `2^-64.9` of phase
1, i.e. **0.000000 bits** — smaller than the 2^60 lag term it is charged beside. Nothing is
absorbed silently: the row appears in the ledger and in the kappa table.

## 5. Toy

`scripts/toy.py` has two arms and one control.

**`realec`** is our hardcoded fixture with real affine EC arithmetic: `y^2 = x^3 + 1` over
`F_600217`, a sextic twist of the BLS12-381 shape (a = 0, j = 0, beta^3 = 1), order 50047
prime, N = 8341 classes, G = (218652, 2045), beta = 43364. The step is the production step,
calling the same `selector.twist_index` the ledger prices. It certifies the **algebra**:
every step asserts `S == [a]G + [c]T`, every table row asserts `J == [rj]G + [sj]T`, and
every merge asserts that the recovered x reproduces T. 1600 trials/density over seven
densities: **879,182/879,182** step assertions, **879,738/879,738** row assertions and
**4,403/4,403** recoveries, after two defects in the v19c toy were found and fixed.

**`walkmap`** runs the same protocol at a rival-scale class universe (N = 2^29.07, the q34
fixture) because `S_ref/2^d` must stay >> 1: our fixture has `S_ref = 91`, so `2^-11` there
is degenerate, which is a property of the fixture and not of the object. The class walk is
a deterministic pseudo-random mapping on `Z_N` (the random-mapping model of a VOW), the
pool is session-global over **every** visited class, and theta gates reporting only.

**`thetacap`** is the control: only DP classes enter the pool.

### 5.1 The two defects, and what they actually were

**D1 (merge gate).** v19c tested `report_step <= ceil(kappa S_ref) + 2^d`, i.e. it folded
the DP report lag into the **search** budget. The ledger charges that lag separately (its
`2^60` term). At `2^-11` the added 2048 steps are 22 `S_ref` of free search, which is
exactly the density trend v19c measured (0.153 at `2^-5` vs 0.607 at `2^-11`). v19d gates
the merge on `kappa S_ref` alone and reports the lag as its own column.

**D2 (recovery orbit element).** v19c blamed the row cache. **That diagnosis is wrong and
is retracted:** `row()` returns the point together with the `(a, b)` that built it, so a
cache hit cannot disagree with the coefficient update — an instrumented re-run of the v19c
code at `2^-5` gives 0 per-step failures. The real defect is `resolve_e`: given pS it
tested `x_S == beta^k x_pS` and returned `e = lambda^k`, but that residue equality means
`S = phi^k(pS) = [lambda^k] pS`, hence `pS = [lambda^-k] S` — it returned the inverse of
the required element (`lambda^k` instead of `lambda^(3-k)`, as `lambda^3 = 1`), so the
recovered x was wrong on 23 of 34 merges. v19d resolves e by testing the six orbit
elements **against the point**, and a merge counts as a success only after `[x]G == T`
verifies; v19c counted unverified merges as successes.

### 5.2 Their table, reproduced at their N and their n

Five rows, kappa = 1.02, our n = 4x theirs (`evidence/toy.json`, `rival_table`):

| curve | dp_bits | ours | theirs | diff | z | verdict |
|---|---:|---:|---:|---:|---:|---|
| q34 | 2 | 0.3983 | 0.395 | +0.0033 | +0.60 | within 2 SE |
| q34 | 4 | 0.4004 | 0.405 | -0.0046 | -1.04 | within 2 SE |
| q34 | 6 | 0.4018 | 0.404 | -0.0022 | -0.67 | within 2 SE |
| q30 | 4 | 0.4033 | 0.399 | +0.0043 | +1.26 | within 2 SE |
| q30 | 6 | 0.4026 | 0.394 | +0.0086 | +2.54 | outside, see below |

Four of five rows are within 2 SE of their printed value, every row is within 2 SE of the
model `1 - exp(-1.02^2 ln(1/0.6097)) = 0.4024` (max |z| 1.82), and the dp_bits spread is
**0.0050** against their **0.0110** — the flatness their sweep was about. The fifth row
(q30/dp6) is the one their own table already puts 2.8 SE below its own model; we land
2.4 SE above their point and 0.1 SE above the model. It is an outlier on their side, reported, not smoothed.

The pipeline arm (budget, phase 2, report capacity, 8-match cap) pools to **0.4028** over
216,000 trials against their **0.3988 +- 0.0021** over 54,000: z = +1.72.

The control is what makes this evidence rather than arithmetic: capping the pool at DPs
gives **0.315 / 0.163 / 0.058** at dp 2 / 4 / 6 — density-dependent, exactly the collapse
v19b measured (0.038-0.208). A port that could not tell the two objects apart would not
show that split.

### 5.3 The density sweep and the gate

theta = `2^-5 .. 2^-11`, `budget = max(S_ref, 8 * 2^d)` (= `S_ref` at every density here),
60000 trials per density on the q34 fixture, session-global report-only pool, collision
event = a class seen by a different walk, `mean T_c/S_ref` in 0.664-0.666 across all seven
densities:

| density | trials | success | CP-95 lower | mean T_c/S_ref | pool records |
|---|---:|---:|---:|---:|---:|
| 2^-5 | 60000 | 0.3982 | 0.3949 | 0.665 | 23943 |
| 2^-6 | 60000 | 0.4000 | 0.3967 | 0.665 | 23943 |
| 2^-7 | 60000 | 0.4011 | 0.3978 | 0.666 | 23943 |
| 2^-8 | 60000 | 0.4007 | 0.3974 | 0.666 | 23943 |
| 2^-9 | 60000 | 0.3989 | 0.3956 | 0.665 | 23943 |
| 2^-10 | 60000 | 0.3975 | 0.3942 | 0.664 | 23943 |
| 2^-11 | 60000 | 0.3975 | 0.3942 | 0.664 | 23943 |

**Gate: min CP-95 lower 0.3942 >= 0.39 at every density, PASS.** Disclosed: at the
originally specified 1600 trials/density the same arm reads min CP-95 0.3640 (point
estimates 0.384-0.427). That is arithmetic, not protocol — a one-sided 95% Wilson bound on
a success rate of 0.399 cannot reach 0.39 below about 8000 trials — so the gate is
certified at n = 60000 and the 1600 reading is printed alongside it rather than dropped.

## 6. What the toy does and does not certify

It certifies three things: the step algebra is exact (every step, every row, every
recovery); the collision object is the walked coverage, so success is flat in the DP
density; and at production-like pool scaling the success rate is 0.3975-0.4011 with a CP-95
lower bound above the 0.39 floor. It does **not** certify the production-size run: the
sweep is a 2^29-class simulation, not 2^20 walks over 2^126 steps. What transports is the mechanism (class function, session pool, DP report), and the
model `1 - exp(-kappa^2 ln(1/0.6097))` that both tables sit on.

## 7. Limitations

- The instruction counts are analytical, as the cost model permits.
- P1 is measured at 5000 production classes, not proved.
- The window error band argument gives `3 * 2^-64` per step; a hostile reviewer can ask for
  the limb-level carry analysis, which is in `scripts/selector.py` (band {5,6}, carries from
  limbs 0..4 only).
- The sweep's class walk is a random-mapping abstraction, not six-class GLV rho; the real-EC
  arm covers the step, the walk-map arm covers the scale.
- This is a generic attack; TNFS is not used.

## 8. Revision history

v17 batch-along-walk died on the review (records lost at encounter). v18 priced a floor of
126.389 with a canonical representative. v19 read their sentence as a truncated quotient
and held: 69/200 agreement. v19b ported the low-window reading, reproduced their 161-row
total to 2^-0.003, and failed the toy gate because its pool held DPs only. **v19c: rival
decode = report-only gates, birthday over S; two toy defects named -> v19d: invariants
green, rival table reproduced +-2SE, CP-95 sweep certified.**

## 9. Reproduction

```sh
python3 ledger.py                                   # every number in section 2
python3 scripts/toy.py --arm realec --trials 1600   # step invariants, 100% exact
python3 scripts/toy.py --arm walkmap --trials 60000 # the gate, 7 densities
python3 scripts/toy.py --arm rival --rival-scale 4  # their table at their N
python3 scripts/merge_evidence.py                   # evidence/toy.json
python3 scripts/verify_certificate.py               # certificate.json, 17/17
```

`ledger.py` prints the row table, the exact-Fraction ledgers, the kappa budget product and
the fallback union bound; `verify_certificate.py` re-derives the claim, the 5000-class
agreement, the 100-key fidelity rows, the invariant counts and the rival table from the raw
evidence JSONs. Model `qwen3.8-flash-next`, harness Hermes.
