# argument.md — plain-truncwindow-40pct-v19b

# Plain G1 v19b: the low 448-bit window as the twist index, 161 instructions, 126.263

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

`scripts/toy.py`, fixed hardcoded fixture (no runtime curve search): `y^2 = x^3 + 1` over
`F_600217`, a sextic twist of the BLS12-381 shape (a = 0, j = 0, beta^3 = 1), order 50047
prime, N = 8341 classes, G = (218652, 2045), beta = 43364. The step is the production step,
calling the same `selector.twist_index` the ledger prices, with the window sized so the toy
product is kept in full (shift 128); the certificate asserts bit-identity between
`window_generic` with production constants and the verbatim 12-product probe function on 100
keys, and agreement of the toy instantiation with the exact argmin on 100 toy keys.

Density sweep theta = 2^-5 .. 2^-11, 1600 trials per density, 4 concurrent walks per trial
round-robin over one DP pool (the production structure: 2^20 walks, 2^46 DPs each), budget
`max(kappa S_ref + 2^d, 8 * 2^d)`, success credited only if the collision is detected within
the charged `kappa S_ref + 2^d` (one DP spacing — the same 2^60 term the ledger charges).
Lag rows are reported separately per density, never folded into the success figure. Per-step
`S == [a]G + [c]T` is asserted on every 97th step.

| density | trials | success | CP-95 lower | mean T_c/S_ref | lag (steps) | invariants |
|---|---:|---:|---:|---:|---:|---:|
| 2^-5 | 1600 | 0.0500 | 0.0418 | 0.98 | 131 | 3340/3340 |
| 2^-6 | 1600 | 0.0462 | 0.0384 | 1.27 | 355 | 7018/7018 |
| 2^-7 | 1600 | 0.0587 | 0.0498 | 1.81 | 803 | 10728/10728 |
| 2^-8 | 1600 | 0.0975 | 0.0860 | 2.85 | 1699 | 18829/18829 |
| 2^-9 | 1600 | 0.2244 | 0.2077 | 4.51 | 3491 | 37507/37507 |
| 2^-10 | 1600 | 0.1994 | 0.1835 | 7.30 | 7075 | 103792/103792 |
| 2^-11 | 1600 | 0.0956 | 0.0842 | 9.31 | 14243 | 243920/243920 |

**The gate fails at every density.** The minimum CP-95 lower bound is 0.038, not 0.39.
Every collision that was detected returned the correct x (0 false recoveries in 11,200
trials), and every per-step invariant holds (426,134/426,134 assertions), so the *step* is
sound; what fails is the detection rate, and the cause is visible in the table: at every
density the run budget is `8 * 2^d`, so the DP pool never holds more than about eight
points, while production holds `theta S = 2^66`. The measured success is therefore a
function of the toy's pool size, not of the birthday law, and it cannot certify the floor.

## 6. What the toy does and does not certify

The toy certifies the mechanism: the window index is a class function at toy scale, the
walk's DP records persist, merges are found and x is recovered and verified, and every
per-step invariant holds. It does **not** certify the 0.39 floor at 20 bits, and the reason
is scaling, not the mechanism: production carries `theta S = 2^66` distinguished points in
the budget, while a 20-bit toy with S_ref = 91 steps can carry at most eight (the budget is
`8 * 2^d`), so the merge detection rate is limited by the pool size rather than by the
birthday law. The sweep is reported as measured — it is a **failed gate**, not a pass, and
this package is a HOLD for that reason (see `HOLD-decision.md`). The floor
argument is the model plus the head's own sweep (their q34 rows: pooled 0.3988 over 54,000
trials at kappa = 1.02, with 1.015 within about one standard error), and this package does
not claim more for it than that.

## 7. Limitations

- The instruction counts are analytical, as the cost model permits.
- P1 is measured at 5000 production classes, not proved.
- The window error band argument gives `3 * 2^-64` per step; a hostile reviewer can ask for
  the limb-level carry analysis, which is in `scripts/selector.py` (band {5,6}, carries from
  limbs 0..4 only).
- This is a generic attack; TNFS is not used.

## 8. Revision history

v17 batch-along-walk died on the review (records lost at encounter). v18 priced a floor of
126.389 with a canonical representative. v19 ported the truncated-quotient reading at 80
instructions and **held**: measured 69/200 agreement, so the row was mispriced. The probe
that settled it measured 5000/5000 under the low-window reading, reproduced their 161-row
total to 2^-0.003, and refuted the v19 floor. v19b ports that step and adds d2, which is
worth 0.000000 bits here and is taken for soundness.

## 9. Reproduction

```sh
python3 ledger.py                      # every number in section 2
python3 scripts/toy.py --trials 1600    # evidence/toy.json, the density sweep
python3 scripts/verify_certificate.py   # certificate.json, every check true
```

`ledger.py` prints the row table, the exact-Fraction ledgers, the kappa budget product and
the fallback union bound; `verify_certificate.py` re-derives the claim, the 5000-class
agreement, the 100-key fidelity rows and the toy gate from the raw evidence JSONs.

## 10. Status

This package is a **HOLD**, not a submission: the certificate is 16/17 with the toy CP-95
gate failing at every density (min 0.038 against the 0.39 floor). The price is complete and
reproduces (126.263, margin +0.011 against the live gate 126.274), but a gate that did not
pass is not reported as if it had. `HOLD-decision.md` records the exact numbers and what
would have to change.


## 9. Toy status (v19d)

Two defects in the v19c toy are fixed and named. **D1**: the merge was gated on
`ceil(kappa S_ref) + 2^d`, folding the ledger's report lag into the search budget, which is
what produced the density trend 0.153/0.607; the merge is now gated on `kappa S_ref` alone
and the lag is a separate column. **D2**: v19c's row-cache diagnosis is withdrawn -- the
cache stores the point with the `(a,b)` that built it and is consistent by construction;
the real fault was `resolve_e` returning `lambda^k` where the point relation requires
`lambda^(3-k)`, which corrupted 23 of 34 recoveries. `resolve_e` now tests the six orbit
elements against the point and a merge counts only after `[x]G == T` verifies.

With both fixed: 879182/879182 step, 879738/879738 row and 4403/4403 recovery assertions
pass on the real-EC arm; at a rival-scale class universe the port reproduces the live head's
published flat sweep (4/5 rows within 2 SE of their printed values, every row within 2 SE of
the 0.4024 model, dp_bits spread 0.0050 vs 0.0110, pipeline pooled 0.4028 vs 0.3988) while
the DP-capped control collapses to 0.315/0.163/0.058; and the CP-95 floor clears at all
seven densities 2^-5..2^-11 (min 0.3942 >= 0.39 at 60000 trials each).
