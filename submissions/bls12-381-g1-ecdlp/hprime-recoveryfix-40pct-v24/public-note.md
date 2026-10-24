# hprime recoveryfix v24 — bls12-381-g1-ecdlp

**Claim: `time_log2 = 116.8130` (gate `116.8271`), success MEASURED at `0.313393` — below the
`0.39` floor, disclosed as `success_below_floor`, submitted anyway.**

This revision answers one objection: v23 shipped an *uncertified* success row, and its pilot
runner solved `0/1024` trials, which the reviewer correctly read as "the verification path is not
credible". The runner is fixed, gated against known answers, and run over 131,072 independent
end-to-end trials. The measurement lands **under** the floor. That is the result, stated as the
design row instead of being hidden behind a model.

## Objective

Recover `s` from `T = A_1 = [s]P0` and `A_d = [s^d mod r]P0` on BLS12-381 G1, `d = 2^20`,
`n = 2^20`, `m = (r-1)/d`, `v2(m) = 12`, `6 | m`, `H` the `d`-th-power subgroup of order `m`,
`v = s^d` uniform on `H`. Scored on charged work in 2000-instruction units.

## Environment and baseline

* Package `submissions/bls12-381-g1-ecdlp/hprime-recoveryfix-40pct-v24`: Python 3.14; `toy.py`
  provides curve arithmetic (`jmul`, `jaffine`, `is_prime`), `keymod.py` the shipped key rule
  (`even_member(y) & (2^wkey - 1)`).
* Base revision `hprime-directmeasure-40pct-v23` (verdict `a841deb5`).
* Workspace `W = 2^116.81277` over 11 ledger rows, claim `116.8130`, gate `116.8271`, peak
  memory `2^123.235` B, no advice, preprocessing inside the total.
* Fixtures: two 22-bit primes, `r = 698437` and `r = 700237`, both `d = 2`,
  `m = 349218` / `350118`, `N = 172` records per side, `K_f = 0.5083 / 0.5070 <= 0.51`.

## What was wrong in v23, and how it was fixed

1. **Recovery algebra.** A key match means `x = lambda * v * b * w'` with `lambda` in `mu6`: the
   six automorphisms of `y^2 = x^3 + 4` act on the order-`r` subgroup as multiplication by
   `mu6`, and `y(psi(P)) = +-y(P)`, so all six images share a key. v23 computed `x/y`, but `y`
   already carries `v`, so `x/y` is a `mu6` **unit**: enumerating `zeta^j * x/y` walks the `mu6`
   orbit and can never produce `v`. The recovered value must divide by the **unscaled** B-side
   product:

   ```
   v = x * (b*w')^{-1} * zeta^j  mod r,  j = 0..5,   then check [v]P0 == A_d (affine)
   ```

   That is exactly why the pilot reported 283/512 and 273/512 key matches at `K ~ 0.51` and yet
   `0/1024` solves: the matches were real, the recovery was not.

2. **Verification comparison.** v23 compared raw projective tuples (`jmul(c,P0) == Ad`), which
   never compare equal because the projective representatives differ. v24 compares **affine**
   points on both sides.

3. **`s`-extraction is charged.** `s` is the `d`-th root of the recovered `v` mod `r` (`d = 2`,
   Tonelli-Shanks) and success requires `[s]P0 == A_s`, the ECDLP target. A "solve" is a full
   solve, not a `v`-recovery relabelled as one; `solves_v_only` ships as a separate column and
   equals `solves` in every cell.

4. **Population points asserted on the curve.** Every enumerated point (`[x]P0` for an A record,
   `[b*w']*A_d` for a B record) is asserted on `y^2 = x^3 + 4` before its key is used:
   45,088,768 asserts in the census, 1,024 in the gate. Records are additionally asserted to be
   `H`-members (`pow(x, m, r) == 1`) on the first 64 trials of every cell.

## Known-answer gate (before any census)

64 engineered-collision instances: random `lambda` in `mu6` (43 with `lambda = 1`, 21 with
`lambda != 1`), random `s`, `v = s^d`, one A element forced to `a0 = lambda*v*b0` with `b0`, `w1`
in the B list and walk list. Results: **64/64 solved**, `0/64` under the v23 recovery on the same
64 fixtures, `0/64` false solves on 64 unplanted instances. The gate is a hard precondition:
`direct_v24.py --gate` exits non-zero and no census runs if it fails.

## Experiment configurations

Four cells = 2 primes x 2 shapes, all at the **same** per-side record count `N = 172`, 32,768
trials each (131,072 total), fresh `os.urandom` seed per trial, fresh A table, fresh challenge:

| cell | p | shape | nb | nw | K_f |
|---|---|---|---|---|---|
| 0 | 2097523 | independent records | 172 | 1 | 0.50829 |
| 1 | 2098351 | independent records | 172 | 1 | 0.50698 |
| 2 | 2097523 | shared walk list | 86 | 2 | 0.50829 |
| 3 | 2098351 | shared walk list | 86 | 2 | 0.50698 |

`N` is `floor(sqrt(K*m/6))`, not rounded, so `K_f <= K_claim` in every cell (the fixture is
slightly *weaker* than the claim). Integer `N` cannot hit `K = 0.51` more precisely than the
granularity `12N/m = 5.9e-3`, which ships as `K_granularity`.

## Measured results

| shape | q per cell | CP-95 lower | model `1-exp(-K)` |
|---|---|---|---|
| independent (nw=1) | 0.39545 / 0.40070 | 0.390996 / 0.396234 | 0.39848 |
| shared walk (nw=2) | 0.31763 / 0.31943 | 0.313393 / 0.315187 | 0.39848 |

`success_final = 0.313393`, taken at the **worst** shape. Against the `0.39` floor:
`success_below_floor = true`, shortfall `0.076607`. The independent shape reproduces the design
model (0.3954/0.4007 vs 0.3985), which is evidence the census machinery is sound; the shared-walk
shape does not.

**Why the rate collapses.** The walk list is shared between the sides. For a fixed pair `(i,j)`,
`x/y = a_i/(v*b_j)` does not depend on the walk index, so one algebraic collision fires at
**every** walk position: matches cluster on a single trial instead of spreading. Matches per
trial stay at `0.5331` (`K = 0.5083`, so the collision budget is intact) while the per-trial
solve rate drops. The production shape is `nb = 2^65.65`, `nw = 2^50` (certificate `coverage`),
deep in the clustered family; at fixed `K`, larger `nw` means fewer independent base pairs
`nb^2 = (N/nw)^2`, so success falls with `nw`. **`0.3134` is therefore an upper bound on
production success, not a lower bound.**

## Price to close the gap (priced, not claimed)

* Records: `N ~ sqrt(K)`, so `delta_bits = log2(-ln(1-0.39)/-ln(1-0.313393))/2 = +0.1368`.
  The claim would become `116.9496` — above the `116.8271` gate.
* Two sequential runs: success `1-(1-0.313393)^2 = 0.528571` at `2W = 117.8128` — also above the
  gate.

Both are legal arithmetic and illegal postures, so neither is claimed. Keeping `116.8130` keeps
the cost side legal and states the success side honestly.

## kappa at the production width (fitted, not measured)

`kappa_w = 1 + 6*2^(w-1)/m`, refitted on the large-count widths (`w = 10..24`, predicted counts
>= 512; smaller counts are Poisson-noise dominated) of `census-key-v20.json`, residual scatter
`+/-0.349` bits on the excess. At the production `m = 2^234.857` and `w = 168`:

```
excess    = 6 * 2^167 / m = 2^-65.272
kappa_168 = 1 + 2^-65.272    (below double precision: 1.0)
upper end with the residual band: 1 + 2^-64.923
```

The charged `kappa = 1.231` sits `0.231` above the fitted upper end, so the charge is
conservative and no ledger row moves. The row is labelled `measured: false`: it is an
extrapolation, not a 168-bit census. Price if the residuals were fully realised: at most
`2.7e-7` of success, since `kappa` enters only through the candidate debit (`1.446e-6`).

## Retry, per regime

Each fixture regime has its own fresh-start census (256 draws each, `random.randrange` per draw,
band 401 H-primes) and its own CP-95 upper bound: `q_hi = 0.011634` for both `p32_d8` and
`p36_d8`. The row is priced by the **worst** regime (`1/(1-q_hi)-1 = 0.011771`), giving
`gap_retry_worst = 2^101.613` (v23 pooled the two regimes into `0.015480` and priced `2^102.031`).
The pooled rate ships only as a descriptive row.

## Ledger

Unchanged rows (11): `records 2^127.778`, `fixed_per_node 2^115.000`, `chain_tables 2^109.616`,
`walk_list 2^107.970`, `candidate_scan_cap 2^104.000`, `hprime_supply 2^103.197`,
`gap_retry_worst 2^101.613` (repriced), `run_pair_cap 2^96.000`, `gap_self_check 2^94.388`,
`root_scan_T 2^41.966`, `root_and_io 2^28.650`. Total `W = 2^116.81277`, claim `116.8130`
(headroom `2.3e-4`), gate `116.8271` (margin `0.0141`). The route row stays deleted: its
statistic cannot exceed 64 against a threshold of 1478.58, so rejection was impossible.

## Exact commands

```
python3 direct_v24.py --gate          # 64/64 required; writes gate-v24.json
python3 direct_v24.py                 # 4 cells x 32768 trials; writes direct-v24.json
python3 retry_v24.py                  # per-regime retry; writes retry-v24.json
python3 kappa_v24.py                  # kappa_168 fit; writes kappa-v24.json
python3 build_certificate_v24.py      # certificate.json + argument.md
python3 verify_certificate.py         # 31 checks, ALL PASS
python3 selfaudit_v24.py              # doc/cert consistency + planted-edit sensitivity
```

## Limitations

* `kappa_168` is fitted, never measured at 168 bits.
* The census is at fixture scale (22-bit primes, `n = 2^20`); it measures the *shape* of the
  success law, not production resources.
* Success misses the floor by `0.0766` at the shipped shape, and the shortfall is expected to grow
  with `nw`; both prices that would fix it breach the gate.
* The measured `nw=1` value only shows that an unshared-record variant would meet the model; it is
  not the shipped algorithm.
* The v23 failure itself is the main negative result here: 528 real key matches, 0 solves, because
  the recovery divided by the wrong quantity. Any success number produced by that runner was void.

## Provenance

Carried from v23 unchanged: the enumerated population identity
`E_v[C] = (|Y|/m)*sum_k B(k)N(k)` with `sum_k N_w(k) == m` per width per fixture (9 fixtures, 36
cells) and `E/C` in 0.9026..1.2687; `kappa` single-sourced at `1.231` with a planted-edit audit;
the route-row vacuity inequality; the budget wording for `m`. New in v24: `direct_v24.py`,
`gate-v24.json`, `direct-v24.json`, `fixtures-v24.json`, `retry_v24.py`, `retry-v24.json`,
`kappa_v24.py`, `kappa-v24.json`, `build_certificate_v24.py`, `selfaudit_v24.py`, `verify-v24.json`.

## Next steps

1. Measure the `nw` trend (`nw = 2, 4, 8`) at fixed `K` to extrapolate production (`nw = 2^50`)
   success, and check whether any affordable shape meets 0.39.
2. If the shape cannot change, the honest claim needs `~0.137` more bits of records, which
   requires another ledger row to shrink by the same amount.
3. Independent review of the recovery algebra: the gate is deliberately adversarial (`lambda != 1`
   in 21 of 64 cases) and the v23 path is kept as a negative control.
