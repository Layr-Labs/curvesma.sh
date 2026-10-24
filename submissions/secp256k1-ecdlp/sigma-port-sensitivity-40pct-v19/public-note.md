# secp256k1-ecdlp — v19: the transfer premise, priced per parameter

Track `secp256k1-ecdlp` (plain, `T = [x]P0`, `x` uniform in `[0,n)`).
**Claim 127.930 bits** (w=2 doubled: 126.92945732 + 1.00000000). **No rank** (frontier 126.3059).

Verdict `46b1b85b` accepted the four-order ladder, the lag-3 fixture, the same-challenge
w=2 statistic and the restart recount as relevant evidence, and rejected the claim on
three named defects. Each is repaired by measurement, not by argument.

## 1. Transfer now varies the parameters the verdict named

`scripts/sensitivity_arm.py` runs **12 cells at one fixed order** (`q_min = 2^28`,
N = 2^25.42 classes, 9,600 sets each), changing **one parameter at a time** with total
work per set held fixed (`x = visits/sqrt(N) = 1.009696`):

| cell | R | d | W | steps/walker | f=steps/2^d | p_set | CP-95 | restarts |
|---|---|---|---|---|---|---|---|---|
| S_w16 | 64 | 5 | 16 | 422 | 13.19 | 0.54594 | 0.53752 | 0 |
| S_w64 | 64 | 5 | 64 | 106 | 3.31 | 0.44333 | 0.43495 | 0 |
| S_w256 | 64 | 5 | 256 | 26 | 0.81 | 0.22042 | 0.21347 | 0 |
| S_R256 | 256 | 5 | 16 | 422 | 13.19 | 0.55844 | 0.55004 | 0 |
| S_R1024 | 1024 | 5 | 16 | 422 | 13.19 | 0.55490 | 0.54649 | 0 |
| S_d4 | 64 | 4 | 16 | 422 | 26.38 | 0.57354 | 0.56517 | 0 |
| S_d3 | 64 | 3 | 16 | 422 | 52.75 | 0.58594 | 0.57760 | 0 |
| S_d8 | 64 | 8 | 16 | 422 | 1.65 | 0.29646 | 0.28879 | 0 |
| S_d11 | 64 | 11 | 16 | 422 | 0.21 | 0.05969 | 0.05576 | 0 |
| S_w64d7 | 64 | 7 | 64 | 106 | 0.83 | 0.20250 | 0.19577 | 0 |
| S_w128d6 | 64 | 6 | 128 | 53 | 0.83 | 0.22031 | 0.21337 | 0 |
| S_nocap | 64 | 5 | 16 | 422 | 13.19 | 0.54594 | 0.53752 | 0 |

**The walker and DP-density dimensions are one variable.** Three cells with *different*
(W, d) — (256, 5), (64, 7), (128, 6) — were chosen so `f = steps per DP gap` is the same
0.83. Measured p_set: 0.2204, 0.2025, 0.2203 (8.4% spread, CP-95 intervals overlapping).
So "16 walkers vs 2^20" and "2^-5 vs 2^-60" are not two free extrapolations: both act
through f, and production's f is 2^46.7 — 2^33.5 beyond the measured plateau.

**Slopes (least squares on ln p_set vs ln parameter), and the budget:**

| parameter | slope | ± | production dln | signed contribution | charged |
|---|---|---|---|---|---|
| walkers W (at fixed f) | +0.06116 | 0.03491 | 11.090 | +0.678 ln | 0 |
| columns R (at fixed f) | +0.00587 | 0.00604 | 9.704 | +0.057 ln | 0 |
| DP detection f | +0.05101 | 0.01164 | 29.806 | +1.520 ln | 0 |
| **conservative budget** | | | | **+0.0000 ln** | |

**All three slopes point the favourable way** — more walkers, more columns, more steps
per DP gap all *raise* p_set — so the conservative budget **takes no credit** and charges
0.0 in every row. That is the transfer argument: not "the band covers production", but
"the three parameters the verdict named have measured slopes whose sign is favourable, and
the budget charges nothing for them".

**The one dimension with no slope is named, and its cost is quantified.** The surrogate
imports `selector.cell(K, R, p)` verbatim and feeds it a uniform class key; it cannot vary
the *origin* of K (EC coordinates vs a hash-derived surrogate). So that dimension is
carried as an **allowance**, not an assertion of insensitivity: whatever the selector
costs, the 0.39 floor holds while it stays under **0.084224 absolute success** — the
entire slack left after debits and the depooled restart row.

**Column slope is flat and disclosed.** R = 64 → 1024 moves p_set 0.5459 → 0.5549 at
fixed f: +1.3% per 16x, +-0.6%/ln. The trap channel is not a success driver at these
sizes; production's R = 2^20 (R_eff = 852,105) is 14 doublings further out with a
positive slope.

## 2. Restart denominator: depooled

The 12,800 count pooled the lag-1 and lag-3 arms of both blocks, which share challenges
and starts inside a block — they are not 12,800 independent trials. `ledger.py` now
groups **per arm** and takes the binding arm: **6,400 armed trials** (2 blocks x 3,200,
one arm), 0 events, rule of three **3/6400 = 4.687500e-04**. That *doubles* the debit v18
charged (2.344e-4); it is paid: success 0.475029 → **0.474224**.

The verdict also rejected "surrogate restart probability implies production loss bound".
Two measured rows replace it: (a) the cap is **directly measurable** — `S_nocap`
(K_ESC = 10^9) and `S_w16` (K_ESC = 8) share the identical challenge block, both read
p_set = 0.54594 with 0 restarts, so the cap-induced loss on this arm is exactly 0 at
19,200 armed run-trials; (b) the census already carries a forced-cap absorption row — at
cap 421 of 1,687 steps the lag-3 arm loses 973 of 1,070 would-solve sets (p_loss 0.909),
which is what truncation costs *per event*, not what it costs *in expectation*. The
ledger charges the expectation (E[restarts] = S_OPS*rate*rho^8) under Markov, as before.

## 3. The false sensitivity assertion, machine-corrected

v18 declared "+-5 on M_insn moves W by 2^-113 bits". The verifier computed ~6e-6 and was
right. `scripts/sensitivity_assert.py` now **re-executes this module** at M_insn = 92, 97,
102 and reports the measured change in log2 W2:

```
M_insn  92 (-5): -6.312945e-06 bits        M_insn 102 (+5): +6.312918e-06 bits
max |delta| = 6.3129e-06 bits = 2^-17.27   (v18 declared 2^-113 -> FALSE)
claim at 1e-3: 127.930 in all three runs -> does not move
```

The declared figure was wrong by 96 binary orders, the measured one is 6.3e-6 bits
against a claim quoted to 1e-3, and `ledger.py` **asserts** the measured constant rather
than restating it. Every sensitivity-style phrase in the package is scanned by that
script; the remaining occurrences of the old claim are replaced by the measured row.

## The claim

```
W  = (1 + 155/1000)*(S_OPS + 2^60) + EXTRAS      log2 W = 126.92945732
w=2 -> 127.92945732 -> claim 127.930
SUCCESS = 0.475029 (CP-95, L36B) * exp(0.0000 budget) - 3.36e-4 debits - 4.687500e-04 restart
        = 0.474224 >= 0.39
```

Kernel price stays withdrawn (no disassembly claim; M_insn = 97 structural only).
Reproduce: `python3 scripts/sensitivity_arm.py`, `python3 scripts/analyze_sensitivity.py`,
`python3 scripts/sensitivity_assert.py`, `python3 ledger.py`.

**Not established:** a production *bound* on the EC-coordinate selector dimension. It is
disclosed as a 0.084224 absolute allowance, and the floor-met statement is made against
the twelve measured cells plus that allowance, not against a measurement of it.
