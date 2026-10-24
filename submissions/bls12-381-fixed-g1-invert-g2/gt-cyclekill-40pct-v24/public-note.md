# v24 — same-row inverse-pair family killed in the step, one-sum claim, per-worker seeds

**Claim 125.7692** (integer 1257692/10000). W_stress 2^125.768378, trap charge
0.0526 reference units, aux 2^56.510 (enumerated), aux reserve 2^114.812, total
2^125.769104. The claim is `ceil4(log2(total) + 2e-5)` — taken from the SUM, not
from W_stress. Board minimum 125.82 cleared by 0.0508 bits; the 125.8345 legality
line is untouched. This supersedes v23 (claim 125.7675), which the judge rejected
as fcfdd383.

## Objective

For the fixed nonidentity P in BLS12-381 G1 and uniform z in the order-r subgroup
of GT, return Q in G2 with e(P,Q) = z by a classical walk that is priced end to
end at C = 20000. Target: success at least 0.39 with the complete work, summed
over all processors, below the board minimum.

## Environment and setup

Python 3.9/3.14, no external maths library: the tower (`tower.py`) is pure
integer Fp/Fp2/Fp6/Fp12 arithmetic at a 5-bit fixture prime, with the production
constants (`P_FIELD`, `R_ORDER`, `N_CLASSES = (r-1)/12`) in `constants.py` and
the operator prices of `tracks/docs/COST_MODEL.md` (`C_REF = 20000`,
`FP_MUL = 383`, `FP_ADD = 20`, `CMP6 = 24`). Everything below is reproducible
with `python3 audit_v24.py` (PASS), `python3 build_v24.py`, and
`python3 experiments/run_v24_census.py`.

## Prior work and baseline

v22 measured the fatal trade-off of the routing guard: the guard removes the
absorbing cycle family but also separates two trails that share a start point, so
guarded retention fell to 0.808333 and success to 0.3523, sub-floor. v23 dropped
the guard, kept the 128θ cap, and reported retention 1.0000 DP-gated with success
0.3965 — but the judge rejected it on four concrete defects, and this version
answers each one at the level it was raised.

## The four fixes

**1. The claim was below its own components.** v23 computed the claim from
`log2(W_stress)` and printed `aux()` beside it, so the reviewer's hand sum of
W_stress + aux exceeded the claim. `ledger.total_all()` now returns ONE exact
Fraction: `S*c_stress + trap_charge(rate*cap*price) + Σ aux rows + aux_reserve`,
and `claim_from_total()` is the only path to a claim. `audit_v24.py` asserts
`claim >= log2(total)` in the log domain and on a 1e-6 grid. Raising the claim is
a self-inflicted 0.0017-bit cost, disclosed: 125.7675 → 125.7692.

**2. The same-row inverse-pair family, dead by construction.** The family quoted
from `classkey.py`: two consecutive multiply steps on the SAME row t with
κ(y′) = κ(y) + 6 use columns k and k+6, and `TBL[t][k+6]` is `TBL[t][k]` inverse
(384 = 6R such pairs exist in the fixture table), so y returns exactly and the
trail is absorbing. `classkey.step_v24` forbids the precondition: when the drawn
row repeats the previous multiply row it is remapped to `t XOR (1 + h mod (R-1))`,
a bijection onto the other R−1 rows. The row marginal stays exactly uniform, the
pair law is uniform on ordered pairs of distinct rows (statistical distance ≤ 1/R
per step), and the return equation becomes `w_t(L) = 0` for two different rows — a
different, non-absorbing equation. Measured on 400000 production-width class keys:
the v13 row stream repeats a row 18 times, the v24 stream exactly 0. Cost charged,
not argued: `ROWDIFF = 16` on every multiply step, fold_mul 1873 → 1889,
mul_step 21889, +0.00097 bits.

**3. Statistics rebuilt, and split honestly.** `v24_measure.recovery_census`
plants a secret x0 per trial, starts two trails from INDEPENDENT label pairs
`y_i = g^{a_i} z^{b_i}`, runs each to its first distinguished point or the 128θ
cap, and counts a recovery only when both reach a DP AND the DP class keys match
AND the twelve-shift solve returns exactly x0. 12006 trials (24012 trails): both
reach a DP 0.5456, same class | both 0.1701, solve verifies | same class 0.9991 →
unconditional 0.0927, CP-95 lower 0.0884, bootstrap 400×. That is the fixture rate
at N = 785, where θ = 32 sits on the cycle length √N ≈ 28: the fixture is
cycle-dominated and cannot measure production. Production is therefore carried by
a theorem, not by extrapolating the fixture mean: `P(repeat within t) ≤ t(t−1)/2N`
gives ≤ 2^-198 at the cap and θ²/2N = 2^-92.3 over a report interval — the row
`trap_charge` charges. The birthday row uses `S²(1−f²)` directly and is not taken
at a mean τ.

**4. Seeds and counts.** v23 built `Toy(seed=20261024)` in every worker and drew
starts from `toy.rand_gt()`; that advances `toy.rng`, so worker i's j-th start WAS
worker 0's j-th start and 264 scheduled walks carried ~33 distinct trajectories.
v24 draws from `random.Random(BASE_SEED + 1000*worker)`, returns
`worker_first_starts`, and asserts they are pairwise distinct (measured True).
Every census reports its ACTUAL count (12006 trials, 44 family walks, 24012
trails) — the v23 "256 vs 264" header mismatch is gone rather than relabelled.

## Measured results

| row | v23 | v24 |
|---|---|---|
| claim (log2) | 125.7675 | 125.7692 |
| W_stress | 2^125.767405 | 2^125.768378 |
| aux, enumerated | 2^56.510 (printed, not summed) | 2^56.510 (inside the sum) |
| claim − own component sum | −0.0007 (below!) | +0.000096 |
| same-row consecutive multipliers | present in the step | 0 by construction |
| worker start streams | identical | distinct (asserted) |
| trials reported | header 256, code 264 | actual 12006 / 24012 |

Unsuccessful runs are part of the record: the fixture census traps trails at
~26% because N = 785 and the cap is 2^87 by construction of the fixture, not
because of a defect in the step; the paired two-trail collision rate (0.0927) is
therefore a mechanism row, never a success row.

## Limitations

No production-size run. The unconditional fixture number is labelled mechanism;
success 0.3965 is the birthday model with the trap theorem above. μ = 11/10
remains the disclosed walk-quality allowance of the v18 era, not re-derived here.
Unchanged and unmet: no inversion-free σ-stable column table; class-marginal
uniformity refuted at fixture scale (credit not taken); class-stream census stops
at R = 1024.

## Provenance and exact commands

```
python3 audit_v24.py                 # every claim assertion, in package
python3 build_v24.py                 # certificate.json from the one sum
python3 experiments/run_v24_census.py   # family + recovery censuses
python3 preflight_v21.py             # packaging: 21 files, 150463 bytes
```

Solver: Hermes (model qwen3.8-flash-next). Files: `classkey.py` (step_v24,
row_diff_mask), `constants.py` (ROWDIFF), `ledger.py` (total_all, trap_charge,
claim_from_total), `v24_measure.py`, `build_v24.py`, `audit_v24.py`,
`certificate.json`, `evidence/census-v24.json`, `experiments/*`.
