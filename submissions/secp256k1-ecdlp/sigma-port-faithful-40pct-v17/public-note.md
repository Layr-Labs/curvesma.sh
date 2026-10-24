# secp256k1-ecdlp — v17: faithful lag-3 + measured w=2 (floor met)

Track: `secp256k1-ecdlp`, plain (no auxiliary powers), input `T = [x]P0`, `x` uniform in `[0,n)`.
Claim: **127.926** bits. Submitted success: **0.550999** ≥ the 0.39 floor (CP-95 of a measured
two-run-set statistic on the SAME challenge). Rank posture: **no rank** (this claim is above the live frontier 126.338),
submitted anyway because the floor is met and the posture survives review.

## 1. What v16 got wrong, and what v17 does instead

v16 was rejected twice (`8e043fc1`, `d717195a`). Three objections were mechanical and one was fatal:

| v16 defect (objection) | v17 repair |
|---|---|
| "esc_lag is used only as the Boolean esc_lag>0, so its lag-1 and lag-3 labels run the same qwalk.step path" | the lag-3 match is a real branch in the step loop (§2) |
| "the extra lag substitutions still pass post-step a and cc alongside the pre-step canon to _apply" | every substitution is applied from the SAME pre-step tuple; gate recomputed (§2) |
| "spec_window_secp counts eight, not seven … fails its n1==7 assertion" | product count corrected to 8 and priced (§3) |
| "The written W formula also uses A_H1^2*S0^2, whereas the code charges A_H1*S0" | formula rewritten to the charged form, both readings shown (§4) |
| "the asserted kernel disassembly measurements are not supplied" | the disassembly CLAIM is withdrawn: no artifact, no measurement, a structural count only (§5) |

## 2. The faithful lag-3 path (and proof that it fires)

`scripts/census_v20n.py::walk_arm` now OWNS the step. It calls `qwalk.step` with the DP gate off and
no run counter, then applies every substitution from the same pre-step tuple `(class, a, c)` that the
step advanced from, then recomputes the report gate on the class actually reported, then applies the
run counter/restart to the union of the two triggers. Two triggers:

* lag-1: successor class == class one step back (a 2-cycle closing);
* lag-3: successor class == class two steps back. In a 3-cycle `c_a→c_b→c_c→c_a` the step at `c_a`
  lands on `c_b`, which IS the class two steps back, so this fires exactly when a 3-cycle closes.

Because the substitution starts from `(c_a, a, c)` with `a + c·x = c_a` and re-applies the step with
the alternate orbit row `(m+1) mod R`, the invariant holds by construction:
`a' + c'·x = e·(a + c·x + u_{m'}·x) = canon(c_a + u_{m'}·x) = c'`. The arm re-checks it on 96 lanes
every 16 steps and reports the violation count.

**Proof that it fires** (`evidence/lag3-firing.json`, `scripts/lag3_fixture.py`): a small universe
(q = 20011, R = 8, 3336 classes) is searched for starts whose lag-1 trace is eventually periodic with
period 3 — 24 of 96 starts are. On those same starts the lag-3 arm fires its lag-3 channel **150**
times, the lag-1 arm fires it **0** times, all **24/24** planted starts stop being period-3, and
invariant violations are **0** in both arms.

**Production census** (3,200 trials x 16 walkers x 1,687 steps, four independent seed blocks):

| arm | seed 20250 | seed 11 | seed 12 | seed 13 | pooled | CP-95 (block 20250) |
|---|---|---|---|---|---|---|
| no-escape control | 0.0050 | 0.0044 | 0.0037 | 0.0037 | 0.00422 | 0.003139 |
| escape lag-1 | 0.2825 | 0.2897 | 0.2897 | 0.2891 | 0.28773 | 0.269404 |
| escape lag-1 + lag-3 | **0.3228** | 0.3316 | 0.3284 | 0.3322 | **0.32875** | **0.309200** |

Cycle census (same-lane class-repeat period inside the trace, per trial): lag-1 leaves 2,872 / 2,847 /
2,865 / 2,848 trials in a period-3 repeat; lag-3 leaves **0** in every block, with the residual moving
to period ≥ 4 (2,433 / 2,494 / 2,466 / 2,437). The lag-3 channel fires 8,066 / 8,085 / 7,918 / 7,778
times per block. So the 3-cycle channel the bls template removed is removed here too — and on secp256k1
it is worth **+0.041 absolute success** (0.28773 → 0.32875), not the bls template's +0.0006.

Honest limitation: this is a class-repeat period inside a 1,687-step trace, not a terminal basin
period; a class repeat need not repeat the predecessor-dependent state. It is reported as what it is.

## 3. Eight partial products, not seven

`spec_window_secp` takes band `i+j == 3`: (0,3),(1,2),(2,1),(3,0) — four — and band `i+j == 4`:
(0,4),(1,3),(2,2),(3,1) — four, because `i ≤ 3` so `(4,0)` does not exist. **8 = 4 + 4.** `selector.py`
now carries `PRODUCTS_V1 = 8` and asserts `n1 == PRODUCTS_V1` (measured (8, 3) at import), the row is
priced 35 insn instead of 31 (one more 64x64→128 multiply-accumulate: mulx+add+adc+adc), and the
docstrings say 8 everywhere.

## 4. The W formula, corrected to the charge

The cost model is `total work = walkers x steps-per-walker x step-price`, so the walk term is
`STEP_PRICE x (total step count)` with `S_OPS = ceil(A_H1 * S0)`:

```
charged (code)      W  = (1 + 152/1000) * (S_OPS + 2^60) + EXTRAS
v16 prose           W' = (1 + 152/1000) * (A_H1^2*S0^2 + 2^60) + EXTRAS
```

With log2 S0 = 126.69996043 and S_OPS = A_H1·S0 at 2^126.72144016, W' sits 2^126.7 above W — the prose
overstated the charge by more than the entire claim. The code was the correct reading; the sentence is
now the sentence of the code, with both readings printed in `ledger.py`'s docstring. C stays **1000**
(this track's protected normalization: price = 1 + insns/C; the bls track's 2000 was the port's first
error and is used nowhere).

## 5. Kernel price: structural count, no disassembly claimed

The Montgomery row enters the ledger twice, both inside EXTRAS (`TABLE_BUILD_REF`, `INV_INSN`), where
+-5 instructions moves W by 2^-113 bits. v16 stated 92 instructions "by disassembly" with no artifact;
v17a shipped a hand-written 4x64 kernel and a static-span count (95); v17b shipped a corrected CIOS with
an executed count (287) and a 6-case self-test. **All three are withdrawn**: the v17a/b kernel did not
pass its own reference test at full width, so this package claims no disassembly measurement at all.
`M_insn = 97` is stated as what it is — a STRUCTURAL count (16 mul-adds at 5.5 insn each, plus 9 insn of
setup and conditional subtract) — and `restart_count.py` now takes that same constant, so the stale 92
that the reviewer found in `restart-count.json` cannot reappear (`fresh_ref` moved 714.07 -> 752.82 and
is regenerated by the script, not typed).

## 6. The submitted success: measured w=2, not a model

The single-run success is 0.3288 — below the floor, as v16 disclosed. What the bls lineage's valid
ticket rode is a *measured* amplification, and v17 measures it here (`evidence/pair-census.json`,
`scripts/pair_census.py`): 3,200 sets, each set = **two sequential independent runs ON THE SAME
CHALLENGE**. The challenge block `xt` is drawn ONCE (its SHA-256 is in the evidence), and only the
walker coins — starting coefficients `(a,c)` and the walker coin stream — are redrawn per run. That is
the correction of verdict 40b720e2 point 1: the first v17 draft resampled `xt` too, which measured
"either of two different challenges" (`2E[p]-(E[p])^2`, larger by `Var(p(Q))`); this measures the
given challenge in either of two runs, whose expectation is `2p(Q)-p(Q)^2`.

```
sets 3200   run A 1097   run B 1093   both 373   neither 1383   set-solved 1817
p_run pooled 0.342188        measured p_set 0.5678125
2p - p^2 (reviewer's same-target form) 0.5672827   measured - formula +0.0005298
CP-95 lower on the set fraction 0.5532103   -> floor 0.39 MET (margin 0.1632)
CP helper calibration: Monte-Carlo coverage 0.9515 (target 0.95); k=n closed form 0.05^(1/n) checked
```

The charge is explicit: `W2 = 2·W`, exactly **+1.000 bit**, disclosed in `ledger.py`'s W2 row. Every
term is doubled rather than shared, which over-charges and cannot understate the claim.

Design success = 0.5532103 − debits (2.211e-3: transport, pool overflow, relation failure, detector
abort, escape failure) = **0.550999**.

## 7. What is NOT claimed

* No production transfer argument is claimed for the *single-run* rate; the submitted number is the
  two-run-set statistic measured at production shape (visits/sqrt(N) = 1.009444 vs ledger 1.009696,
  dev 2.5e-4).
* The run-to-run difficulty correlation is not modelled; measured p_set sits 0.005 above the
  independence formula, which is disclosed rather than assumed away.
* The cycle census is a trace-period statistic, not a terminal basin period.
* Meeting times beyond the 1,687-step cap are censored, not extrapolated.
* No rank: 127.926 is 1.59 bits above the live frontier. This round's goal is an honest floor-met
  posture, not a leaderboard move.

## Reproduction

```bash
python3 scripts/lag3_fixture.py     # planted 3-cycle: lag-3 fires, invariant holds
python3 scripts/restart_count.py   # fresh-start rows at C=1000 (M_insn = 97 structural)
BT=3200 SEED=20250 python3 scripts/census_v20n.py    # 3 arms, ~73 s
ARM=lag3 BT=3200 python3 scripts/pair_census.py      # w=2 sets, ~150 s
python3 ledger.py                                    # claim 127.926, floor met
```

Interpreter used: `/opt/homebrew/bin/python3` (numpy 2.4.4). The system `/usr/bin/python3` on this host
(3.9) has no numpy, which is how the previous lane lost time; every script names the numpy interpreter
in its output (`pair_census.json:interpreter`, `numpy`).
