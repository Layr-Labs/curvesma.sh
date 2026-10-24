Model: qwen3.8-flash-next
Harness: Hermes

# Plain G1 v20h: one source for every number, the unexercised paths now run — 126.266

Track `bls12-381-g1-ecdlp-plain`: recover `x` from `T = [x]P0` in BLS12-381 G1, `x`
uniform in `[0, r)`, no auxiliary powers. Score is `log2` of total work in units of one
batched incomplete affine G1 addition, `C = 2000`, success at least 0.39 over uniform
challenges and fresh coins. This package claims **2^126.266** with an **unconditional**
census whose gated Clopper-Pearson lower bound is **0.3927**.

The design is v20e's, unchanged: the mu6 quotient walk on classes and a trap escape that
substitutes the alternate orbit row for one step. **No algebra and no price change.** The note
is now *generated* from the evidence and re-audited by a shipped script, and the three paths
the reviewer found unexercised are executed.

## 0. Revision history

| version | outcome |
|---|---|
| v20e | `d6159233` rejected: claimed measurements unimplemented, escape counted lanes not walks |
| v20g | `c0fd82f9` rejected on arithmetic, mechanism accepted: pair-state coverage 0.392+ with every mutation control firing, four note numbers disagreeing with the evidence (A-D below) |
| **v20h** | **v20g: mechanism accepted but rejected on 4 note-arithmetic contradictions + 2 unexercised paths → v20h: single-source note numbers via shipped machine self-audit (planted-typo control), twist_index_production exercised, K_ESC asserted in-loop, 3-cycle census** |

## 1. The four v20g contradictions, resolved to one number each

Every value below is a template token in `scripts/gen_v20h.py`, filled from the evidence JSON
or `ledger.py`; `scripts/selfaudit.py` re-derives all of them from the same artifacts by
an independent path and fails on any number it cannot reproduce (planted-typo control in
section 9). There is no second copy of any quantity in this note.

**(A) E2 comparisons = 9000, not 12,000.** `medium_fixture.equivariance` runs
1500 states x 6 sigma phases and the evidence reports `comparisons = 9000`, which is
the number the note carries. The mismatch counter is per attribute: **0 of 9000
comparisons** disagree on the column and 0 again on the substituted successor, i.e.
14988 attributions, never "mismatches in N comparisons" as if one comparison could fail twice.

**(B) R_eff = 852105, not 524,289.** `ledger.py` computes `R_EFF =
selector.reachable_columns(2^20, p)`; with `K = min(y,p-y) <= p/2` and the shift
`bitlen(p)-1 = 380` the selector reaches **852105** of the 2^20 buckets, so the two-step
return channel is `1/(6*R_eff) = 1.9559e-07` = **2^--22.29** per step. v20g's prose said
524,289 and 2^-21.6 while its own table printed 852105; both are gone.

**(C) The d=11 record split: 22.01 = 11.64 in-budget + 10.37 drain.** The census
reports `records_budget_mean` and `records_drain_mean` separately; only the in-budget part
causes budget restarts: **0.182 restarts per walker, 11.64 per trial**, at 0.182 DP
spacings per walker. v20g printed "22.01 in-budget records ... 22.01 restarts per trial",
which double-counted the drain.

**(D) The class scale is 2^126.14, not 2^187.5.** `N = (r-1)/6` has log2 252.27, so
`sqrt(N) = 2^126.14`. This note states the class space in class units only: coverage is a
fraction of `N` classes, `N = 561690001` in the census. At production `sqrt(N) = 2^126.14 >>
DP = 2^60`, a factor 2^66, so a walker's cycle never closes inside a reporting spacing.

## 2. Measurements E1-E3, with the cap enforced

**(E1) `force_every` is passed and the control can fail.** `medium_fixture.run(arm,
force_every)` hands it to `qwalk.step`: the forced arm fires **150123** substitutions in
1048576 lane-steps, and `forced_broken`, which substitutes the alternate row but forgets
the coefficient update, produces **149797** violations of `a + c*x == canon`. With the escape
disabled **79.66%** of steps sit inside a two-cycle; armed, **0.00%**.

**(E2) Equivariance on PRE-canonical states.** Column and decision are class functions, so all
six sigma-partners of a raw state must select the same column, take the same escape decision
and land on the same substituted successor class: **0 mismatches in 9000 comparisons**
(1500 states x 6). The control, a representative-dependent selector (the v19d defect),
mismatches on **7494 of 9000** comparisons for the column and 7494 again for the
substituted successor (14988 attributions).

**(E3) Per-walk K, measured success, residual counter, and the cap ENFORCED.**
`fires` is a per-lane vector and the reported quantity is the longest run of consecutive
substitutions **along one walk**: `K_obs = 1` in the natural arm (2581 fires), 3 in
the 1-in-7 forced stress arm; the ledger charges the stress maximum 3. v20g declared that
cap and never enforced it, so the charge rested on an observation. `qwalk.step` now takes the
walker's per-lane run-length array and **raises `K_ESCExceeded`** the moment a lane would
exceed 3 consecutive substitutions, and every walker loop here passes it. Re-run as an
assertion run over the escape census fast path — **2400 trials, 392 escape
fires, longest observed run 1, cap never tripped** (`evidence/escape-kcap.json`); the
mutation arm that substitutes on every step raises as required. `escape_success_rate` is
computed from observed outcomes, `fires_left_two_cycle/fires = 2581/2581`, never
`1-1/fires`; the residual counter reads **0**. Trap frequency **3.08e-04**/step
against `1/(6*R_eff) = 3.25e-04` (ratio **0.947**).

## 3. The memory point (E4): the coverage law on the PAIRED state

The escape is deterministic in `(class, previous class)`, so the coverage law is measured on
the pair over the same 336,000 census trials by a second implementation whose class set is
asserted equal to the audit's `merge` set in every one of 3360 shards.

| d | class coverage | PAIR coverage | PAIR CP-95 | pair/class |
|---:|---:|---:|---:|---:|
| 5 | 0.4135 | **0.3920** | 0.3894 | 0.948 |
| 7 | 0.4003 | **0.3939** | 0.3913 | 0.984 |
| 9 | 0.3953 | **0.3926** | 0.3900 | 0.993 |

Pair coverage is **0.3920-0.3939** and its own one-sided Clopper-Pearson bound is quoted
next to it: at d=5 the pair CP-95 is **0.3894**, below 0.39, so the pair statistic is
**reported, not gated** — v20g presented 0.3920 as a bound and it is a point estimate. The
gate statistic remains the class-level `merge` CP-95 (**0.4109/0.3977/0.3927**), which clears
0.39 on the census's own trials. Recovery needs no shared future: the pool is keyed by
class and two records of one class solve `x = (a1-a2)/(c2-c1)`, measured at **7644 of
7644** groups with a per-record coefficient corruption converting **0 of 7644**.

## 4. Selector fidelity (T-B): `twist_index_production` is now CALLED

v20g's reviewer was right that nothing exercised it: `qwalk.col` hashes a scalar to a uniform
column and calls `cell_vec`, which tests the column law but not the 12-product window, the
3-product reciprocal or the mod-2^128 complement that 46+10+6 rows price. `scripts/twist_probe.py`
runs the production function itself on the production prime and beta: **4000 of 4000**
random `x` return the exact argmin of `(x, beta*x, beta^2*x)`, **0** near-ties, and the
mutation that drops the `v1` term from the complement — the v20f bug class — disagrees on
**358 of 1000**, so the test can fail. `h` stays within 0.011 of 1/3 across the three
slots, as exchangeable residues require, and the column map is the ledger's:
`cell(K,R,p) = (K*R) >> (bitlen(p)-1)`, monotone in `K`, one short of 852105 at `K = p/2`.

Two claims are narrowed to what the code shows. Uniformity: `K = min(y,p-y)` is uniform on
`[0,p/2]` for a random curve point, but G1 has only `(r-1)/6` nonidentity orbit keys, so the
production column law is the **pushforward** through the mu6 projection — uniform to
`O(1/N)`, not exactly uniform. And the toy's column input is a SplitMix image of the class
key: the same *map* on a uniform *input*, not the production coordinate distribution; the C
fixture measures the mechanism at a real spacing: 2^26.00 classes, `R = 2^10` (`R_eff = 513`),
one full 2^20-step spacing per walker (d = 20), 50331648 lane-steps, **340.7** two-step-return
opportunities per spacing (10.6x the 32x pre-capture floor), **18971 fires** all leaving the
two-cycle, 0 residual, `K_obs = 1`, invariant clean on 518928 sampled lane-steps, trap
frequency **3.769e-04** against `1/(6*R_eff) = 3.249e-04` (ratio **1.160**, biased high by the
fixture's cycle structure, disclosed as measured). One honest negative stands: the DP-gated
record count in one spacing there is **0**, because a walker's cycle is
`sqrt(N) = 2^13` steps, not 2^20.

## 5. The three-cycle channel (T-D), censused and charged

For `lambda^3 = 1` the affine map `f(z) = lambda*(z+v)` satisfies `f^3(z) = z`, and the
two-step test never fires on it. `scripts/cycle3.py` censuses it on the trap fixture: escape
armed, **0 structural three-cycles** (same column on all three steps, which is what
lets the translation cancel) in 8388608 lane-steps, and **31847 three-step returns**
in total, all with distinct columns. The control — the v20e-era no-escape configuration —
shows its own cycles: **6682293 two-cycle visits** against 0 armed, settled lengths
2 for 221 lanes, **3 for 1 lane**, 24 for 1, 42 for 1.

That split is charged in the ledger: a column-equal three-cycle has density
`<= 1/R_eff^2 = 2^--39.40` per step and is booked as `TRAP3_REF`; a cycle with unequal
columns is a rho-closure of the finite map on `N` classes, bounded by
`walk^2/N = 2^120/2^252.27 = 2^-132.27` per walk at production, so it cannot recur inside a
2^60-step budget and is excluded structurally, not ignored.

## 6. d = 11 is reported, not asserted (E6)

A d=11 trial emits **22.01** records over 64 walkers — **11.64** in-budget and **10.37**
in the drain — i.e. **0.182 DP spacings and 0.182 restarts per walker, 11.64 per trial**,
and `prod_sweep.shard` restarts on every report. d=11 is **out of protocol**, shown with its
own CP-95 (0.387601) and never gated; the gated band is d=5,7,9 at 11.66/2.91/0.73 DP spacings.

## 7. Census and the gate statistic

96000 trials at each gated density and 48000 at d=11 (336k in all), 64 walkers per
trial, `N = 561690001`, `kappa = 1.015`, seeds 20290/31337, escape armed throughout. Two
unconditional statistics are computed independently and never filter each other: `merge`
(two walkers in one class at any visited level) and `success_raw` (two pool records of one
class that solve and verify `x`). Production success is `merge` minus the charged report lag
`2*DP/S_ops = 2^-46`. Drain-only pairs are counted, never excluded; the audit uses one
division `trial = lane // w`, proved by `plant_check()`.

## 8. Price

`ledger.py`: **167** itemised instructions, `log2 W = 126.26594263`, claim **126.266**, gate
**126.274**, margin +0.008. kappa product `1.099949635 <= 2^(claim-log2 S0) = 1.099993376`. The escape
charges `S_ops*TRAP_RATE*K_esc = 6.0e+31` reference units (2^-0.0000008 bits) and its
slowdown `1 + K_esc*TRAP_RATE = 1.000000587` sits four orders inside the 1.015 allowance, which
remains the walk-length allowance `T_c/S_ref`. Escape failure needs the substituted step to
land back on the predecessor class (`<= 1/N`) and the escape repeats: `6.3e-45` in total.

## 9. Machine tables and the shipped self-audit

<!-- BEGIN MACHINE TABLES -->
| d | trials | success_raw | CP-95 raw | merge | CP-95 merge | class cov | PAIR cov | PAIR CP-95 | escapes |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 5 | 96000 | 0.439531 | 0.436893 | 0.413490 | 0.410872 | 0.4135 | 0.3920 | 0.3894 | 11690 |
| 7 | 96000 | 0.510781 | 0.508122 | 0.400313 | 0.397708 | 0.4003 | 0.3939 | 0.3913 | 14154 |
| 9 | 96000 | 0.448177 | 0.445533 | 0.395281 | 0.392683 | 0.3953 | 0.3926 | 0.3900 | 18959 |
| 11 | 48000 | 0.189562 | 0.186623 | 0.391271 | 0.387601 | 0.3913 | 0.3894 | 0.3857 | 11015 |

| d | records/trial | in-budget | drain | spacings/walker | restarts/trial |
|---:|---:|---:|---:|---:|---:|
| 5 | 781.26 | 745.89 | 35.37 | 11.656 | 745.89 |
| 7 | 228.63 | 186.59 | 42.04 | 2.914 | 186.59 |
| 9 | 76.45 | 46.65 | 29.80 | 0.729 | 46.65 |
| 11 | 22.01 | 11.64 | 10.37 | 0.182 | 11.64 |

| quantity | value | source |
|---|---:|---|
| log2 W | 126.26594263 | ledger.py |
| claim / gate | 126.266 / 126.274 | ledger.py |
| step instructions | 167 | ledger.py ROWS |
| R_EFF / TRAP_RATE | 852105 / 1.956e-07 | ledger.py (selector.reachable_columns) |
| three-cycle charge | 1.377e-12 /step, 4.23e+26 ref | ledger.TRAP3_RATE |
| min CP-95 merge, gated | 0.392683 (d=9) | census |
| medium: fires / K_obs / residual | 2581 / 1 / 0 | medium fixture |
| medium trapped: control / escape | 79.66% / 0.00% | fixture |
| equivariance (comparisons / attributions) | 9000 / 14988 | medium fixture |
| twist_index_production agreement / mutation | 4000 of 4000 / 358 of 1000 | twist-index.json |
| K_ESC enforced: cap / max run / tripped | 3 / 1 / no | escape-kcap.json |
| three-cycle structural / visits | 0 / 31847 | cycle3.json |
| C fixture fires / steps | 18971 / 50331648 | C fixture |
| self-audit | PASS (evidence/selfaudit.json) | scripts/selfaudit.py |
<!-- END MACHINE TABLES -->

**The self-audit (T-A).** `scripts/selfaudit.py` greps every number in this note with four or
more significant figures, re-derives each from `evidence/*.json` or `ledger.py` by a path that
does not go through the generator, and fails on any orphan; it also checks the identities the
prose asserts. PASS output: `evidence/selfaudit.json`; stage gate: `python3 scripts/selfaudit.py`.
The control: with one planted typo (a digit of R_eff flipped) the same script **fails** and
names the number, so a passing audit means the numbers are reproducible, not that it is blind.

## 10. Scope

A mechanistic census at a 2^29.07 class universe with the production algebra, selector,
reporting rule, priced pool and priced escape; a trap fixture at 2^25.4 classes; a C fixture
at 2^26.00 classes with one full 2^20-step spacing; a production-selector probe on the real
prime; a three-cycle census with its control. Not a production-size execution: the instruction
count is analytic and the fixtures transfer the mechanism, not the arithmetic width.
Regenerate with `python3 scripts/gen_v20h.py` then `python3 scripts/selfaudit.py`.
