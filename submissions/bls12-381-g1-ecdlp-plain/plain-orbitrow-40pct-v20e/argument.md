# Argument — plain-orbitrow-40pct-v20e

## Claim

On `bls12-381-g1-ecdlp-plain` (input `T = [x]P0`, `x` uniform in `[0, r)`, no
auxiliary powers), the algorithm specified here recovers and verifies `x` with
probability at least 0.39 using at most **2^126.266** reference additions at `C = 2000`,
i.e. `W = (1 + 167/2000) * (1.015^2 S0^2 + 2^60) + EXTRAS` with `log2 W =
126.26594270`. The success claim is carried by a 192,000-trial **unconditional** census
whose gated CP-95 lower bound is **0.3950** (coverage rate, d=5,7,9) and by the ledger's
model success 0.399014.

## 1. What this revision answers (the four defects of `e2bd4bef`)

| defect | answer | where |
|---|---|---|
| D1 audit double division: `rtr` = repeated trial array indexed at the lane, then `/w` again, so trial `t`'s report was filed under `floor(t/64)` | `trial = lane // w` once, in the single `audit()` function the census and the construction check share; `plant_check()` plants reports at known `(trial, level)` and asserts the attribution exact, that one planted walker marks nothing, and that the old key disagrees | `scripts/prod_sweep.py`, certificate `audit_plant_check_passed` |
| D2 `detect(only_merged=True)` enforced `success <= merge` | filter deleted; two unconditional statistics computed independently and compared post hoc, both directions reported; drain-only pairs counted and printed with their fraction | §3, tables |
| D3 trap aborts destroy coalescence at 2^37.4 opportunities per DP spacing; escape unspecified and unmeasured; 1.015 unjustified | coalescence-preserving equivariant escape designed, proved, priced (167 rows), measured at a fixture in the real regime, kappa split | §4 |
| D4 note/evidence conflict (145,798 vs 275,798; 19,760 vs 19,311) | note tables emitted by the certificate and byte-compared; decompositions are exact integer identities; the pair/trial distinction is printed as two fields | certificate `note_tables_byte_equal_to_certificate`, `sums_reconcile_exactly` |

## 2. The step and the price

The step is the 448-bit-window twist index: `h = argmin` over `x, beta x, beta^2 x` from
the low 448-bit window of `x * floor(2^448 beta/p)` (12 partial products, `i+j in
{5,6}`), `K = min(y, p-y)`, `m = (K*R) >> (bitlen(p)-1)`, `S <- S + Tab[h][m]`,
`(a,c) <- (a + r[h][m], c + s[h][m])`, then the escape row of §4. The ledger prices
167 itemised instructions on top of one reference addition: index 46, v0 10, v2 6,
argmin 14, sign/key 17, slot/address 7, loads/stores 28, coefficient update 12,
checkpoint 12, **trap escape 6**, inversion share 1, counters 7, meter 1. That is +6
over v20d's 161; the 64-step abort row is removed, not added to.

`drain = workers*2*DP = 2^81` steps and the one-spacing `+2^60` lag term are inside
`W`; the near-tie fallback is an exact union bound at `3*2^-64` per step.

## 3. Unconditional statistics (D2)

Protocol: `theta = 2^-d` gates reporting only; a walker that reports restarts at a
fresh `[a0]G + [c0]T`; budget `kappa*S_ref`; pool = records; detection = two records of
one class from two different walkers solving `x = (a1-a2)/(c2-c1) mod q`, verified
against the trial's `x`; merge audit = two different walkers in one class at any
visited level inside the budget. Both are per-trial indicator statistics computed from
the same replay with **no conditioning**: `detect()` takes no merge argument, and the
audit never looks at reports' recovery.

The production success statistic is the coverage rate: two walkers that ever occupy one
class have two records with probability `1 - 2*DP/S_ops = 1 - 2^-46`, and the ledger
charges that lag explicitly. The pool statistic needs *records*, and the toy's pool is
`budget/2^d` per walker (12 records at d=11, `2^46` in production), so at the coarse end
it under-reads; at the fine end the toy's drain window is a large fraction of its
budget and it over-reads. Both are printed. The gated band is the densities whose
walker restarts at a DP inside its budget (11.7, 2.9, 0.73 spacings at d=5,7,9); d=11
(0.18 spacings: no restart at all) is reported in full and excluded from the gate, with
its point estimate 0.3915 above the floor and its CP-95 0.3878 below it.

Drain-only pairs — verified detections on trials the audit did not mark — are counted
in `success_raw` and printed (`drain_only_pairs`, `drain_only_trials`,
`drain_only_fraction_of_merges`). At toy scale the drain/budget ratio is 0.17–0.69, so
these pairs are a real fraction; in production the drain is `2*DP/S_ops = 2^-45` of the
search and the channel is empty. The scaling is stated; the channel is not deleted.

## 4. The trap escape (D3)

**Trap.** `s -> s1 -> s` with `s1 = canon(s + u[m(s)]x)`; structural rate `1/(6R)` per
step, `R = 2^20`, so `2^37.4` opportunities per `2^60` DP spacing. A 64-step abort plus
an independent restart destroys coalescence at that rate, so it is not a repair.

**Escape.** At class `s1`, if the successor would be the previous class `s`, apply the
alternate orbit row at `s1` for one step — column `(m+1) mod R`, same orbit
construction — then continue. No abort, no restart, no re-randomisation.

**Equivariance.** The row is orbit-generated, so `J(sigma S) = sigma J(S)` for all six
`sigma in mu_6`, and the alternate row is the same construction one column over. The
decision compares class fingerprints, which are sigma-invariant, so all sigma-partners
decide identically. Measured with the escape armed and forced to fire on every step:
0 mismatches over 12,000 comparisons (2000 classes x 6 sigma; normal column, alternate
column, successor, alternate successor) and 0 violations of `a + c*x == canon` over
599,296 sampled lane-steps.

**Coalescence.** The walk state for coalescence is (class, previous class). The escape is
deterministic in that pair, so two walkers sharing it share the future exactly; two
walkers that merely *meet* at a class can differ only if exactly one is closing a
two-cycle there, an event at the trap rate itself, which is what the model debits.

**Measurement in the regime.** `scripts/medium_fixture.py`: `q ~ 2^28`, `R = 2^10`,
spacing `2^20`, so the trap period is 6,144 steps against a 1,048,576-step spacing —
**183.5 opportunities per spacing** (structural 170.7; production 2^37.4). Over 256
lanes x 16,384 steps: trap frequency `1.750e-4` per step vs structural `1.628e-4`
(ratio **1.075**), **734 escape fires, all of which left the two-cycle, 0 residual
traps, K_esc = 4** consecutive substitutions at most. The unarmed control spends
**56.1%** of all its steps inside a two-cycle — absorbing, hence the withdrawal of the
abort.

**Price and kappa.** Per-step: 6 instructions (class compare 3, alternate-slot select 2,
join 1). Per event: `S_ops * TRAP_RATE * K_esc = 6.5e31` reference units, i.e.
`2^-0.0000009` bits — the abort row is replaced. Escape slowdown
`1 + K_esc*rate = 1.000000636`, four orders inside the `1.015` allowance, so the trap
regime consumes none of kappa; the 1.015 remains the walk-length allowance
`T_c/S_ref`. Residual failure would require the substituted step to land back on the
predecessor class, probability `<= 1/#classes` per event, after which the escape
repeats; the model debits `S_ops*TRAP_RATE*K_esc/#classes = 6.9e-45`.

## 5. Census

48,000 trials per density at d = 5,7,9,11 (192,000 total), 64 walkers per trial,
`kappa = 1.015`, `N = 561,690,001` classes, seed 20290, escape armed throughout. Every
number in the note's tables is emitted by `scripts/gen_cert_v20e.py`; the certificate
re-reads the note and byte-compares the block, recomputes both CP-95 bounds from
`(k,n)`, and checks the exact decomposition identities.

## 6. Reproduce

    python3 ledger.py
    python3 scripts/medium_fixture.py
    python3 scripts/gen_cert_v20e.py
    nice -n 10 python3 scripts/prod_sweep.py --trials 48000 --chunk 100 \
        --dens 5,7,9,11 --workers 8

The census is ~2 min wall on 8 workers. Raw shard files are not shipped (1,920 files
would breach the archive cap); the reduced evidence carries every count the certificate
uses and the command regenerates the shards.
