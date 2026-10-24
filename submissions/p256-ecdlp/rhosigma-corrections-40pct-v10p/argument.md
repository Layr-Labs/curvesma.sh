# argument.md — v10o

## 0. What is claimed, precisely

For Q = [x]P on NIST P-256 with the SEC 2 §2.4.2 parameters (including a = −3), the two-table
affine sigma walk recovers x with

* **work**: at most 2^128.9437 target-group-operation equivalents under the C=1000 scalar AArch64
  instruction normalization (`evidence/ledger.json`), and
* **success**: p per instance ≤ 2^-61.55 at production parameters, w=2 amplification 2^-60.55.

The work bound is unchanged from v10n and re-audited here. The success bound is **new**: the θ²
extrapolation is withdrawn, the production figure is derived from the measured coincidence law, and
the collapse that v10n asserted on a flawed control is now attributed by a factorial experiment.
The claim **does not meet** the 0.39 floor; see §6.

## 1. The walk

Two 2048-row tables T0/T1 of points [t]G, t uniform mod n from four independent splitmix64 streams
(256-bit seed, `recproc_seed(seed[4])`). A walker holds a point P_i, coefficients (a_i, c_i) with
P_i = [a_i]G + [c_i]Q, a carried previous row `prev`, and a 4-deep ring of recent x-coordinates.

Step (shipped, one TU):

```
par  = step parity
c    = mix64(x ^ (par ? SA : SB))        SA = 0xa5a5..., SB = 0x5a5a...
j    = (c ^ prev) & 2047
prev = j
T    = (par ? T1 : T0)[j]
lambda = (y_T - y) / (x_T - x)  or  (3x^2)/(2y) on doubling, batched over 16 walkers
x' = lambda^2 - x - x_T ;  y' = lambda(x - x') - y
canonicalise y' even (negate a, c with it)
ring[(ridx++)&3] compared against x' BEFORE the store: x' == x(i-4) => do_escape
commit locals, then (out of line) trap_escape: add fixed row 1, a += TLOG[1], resalt prev
DP test: mix64(x' ^ DPSALT) & (2^62 - 1) == 0
```

Records carry (x, a, c). A cross-walker repeat of x at a distinguished point yields
x = (a1 − a2)/(c2 − c1) mod n, checked by [x]P == Q.

## 2. Why the guard exists, and what it costs

Pre-guard the row was a function of the class alone, so rows repeating with period 2 plus sign
conditions (2 free bits of 4) return H exactly — measured lag-2 precursor 0.000536 = 2^-11, closure
rate 2^-24 per 4-step window (`evidence/trap_census.json`), i.e. 2^25.6 closures in a 2^49.6-step
production walk. Without a guard every walk dies inside a short cycle long before it is useful. The
guard (`prev` XOR) breaks the row repetition. It is not optional.

Its cost is the subject of this package. Under the guard, walkers that reach the same point with
different `prev` registers take different next steps, so the coalescent tail that the birthday
argument relies on is gone: a merge no longer delivers a shared distinguished point. Recovery must
come from an accidental record–record coincidence instead, and that channel is much weaker.

v10n measured the weakness but tested it against a control built on a *fixed* salt instead of the
parity SA/SB rule the pre-guard cells used. The verdict caught that. §3 replaces it.

## 3. Factorial attribution (the correction the verdict demanded)

`evidence/census_v10o.c`: three mechanisms as independent compile-time switches — index salt rule,
carried prev, ring+escape — six arms, four DP densities, 3000 trials each, **one seed (424242)** for
all 24 cells. Identical fixture arithmetic to `census_v10f41.c`/`census_guarded_v10n.c`.

Fidelity: the shipped arm (P2) reproduces the v10n binary **exactly** on 200 trials at dp4
(6 successes, 10148289 records, 6 meets, 1 self, 8 ring hits).

| arm | dp4 | dp6 | dp8 | dp10 |
|---|---|---|---|---|
| P0 pre-guard, parity salt | 0.24833 | 0.23467 | 0.24100 | 0.23967 |
| F0 control, fixed salt | 0.26300 | 0.26633 | 0.25667 | 0.23933 |
| P1 guard only | 0.03333 | 0.00733 | 0.00300 | 0.00133 |
| F1 guard, fixed salt | 0.03533 | 0.00967 | 0.00200 | 0.00100 |
| P2 shipped | 0.03333 | 0.00733 | 0.00300 | 0.00133 |
| P2s shipped, sentinel ring init | 0.03333 | 0.00733 | 0.00300 | 0.00133 |

* salt rule: x1.059 (0.26300 → 0.24833) — real, small, and it made v10n's control *generous*;
* carried prev alone: **x0.134** (0.24833 → 0.03333), −2.90 bits;
* ring + escape: **x1.000 exactly** — identical success counts in all four cells (100/22/9/4), and
  identical for the sentinel-init variant, so neither the trap hardware nor its seeding quirk moves p.

P0 is a saturating plateau: mean 0.2410, spread 1.5 pct, flat in ℓθ over 49.6…3172.8, matching the
shipped pre-guard fit A = 0.24212. P1/P2 are not a plateau at any A — the best saturating fit sits on
the grid floor with rms 0.032, larger than the entire data range.

## 4. Production success without a θ² law

At fixed Λ, with w walkers, trail length ℓ and budget S = wℓ:

```
nu = R/sqrt(2N') = (w·ell·theta)/(w·ell/sqrt(Lambda)) = theta·sqrt(Lambda).
```

ν is an identity at fixed budget, not a fitted quantity. Fixture ν(dp4) = 2^-4.37; production
ν = 2^-62.37. Measured κ = p/ν on guarded arms: 0.6885, 0.6059, 0.9915, 1.7622 (dp 4/6/8/10).
Carrying the maximum κ = 1.7622 (the reading most favourable to the algorithm):

```
p_prod <= 1.7622 * 2^-62.3685 = 2.963e-19 = 2^-61.551
w=2      = 1-(1-p)^2          = 2^-60.551
floor 0.39 -> short by 59.2 bits
```

The same data as a density exponent: α_θ = 0.87 ± 0.09 (weighted LS of log p on log θ, 4 points,
max |residual| 0.63 octaves). The sign is the point — p **decreases** as density falls toward
production's 2^-62; it does not improve.

### 4.1 θ² withdrawal and the band

v10n wrote p ∝ θ². Slope 2.00 against a measured 0.87 ± 0.09 over 58 octaves is a 65.5-bit
over-penalty; the law is withdrawn, and with it the 2^-119…2^-121 bracket.

Replacement band over the defensible family {measured power law (2^-55.4), ν-linear (2^-61.55),
constant (2^-4.91)}: **log2 p ∈ [−61.55, −4.91]**. Every member misses the floor. The most generous
member is the constant p = 0.03333 — 6.57× short of the 0.218975 that w=2 requires (exact: 1−sqrt(0.61)); the linearised 0.195 understates the threshold. The claim therefore
rests on no choice of extrapolation: even refusing to extrapolate at all fails the floor.

## 5. Ring engagement, coalescence, witness validity

**Engagement.** Measured closure rate 5.09e-8…6.65e-8 per slot-step = 2^-24.03 across all four
densities, against the ledger's modelled 2^-24 (97.9 pct). At production ℓθ = 2^49.63: 5.09e7
closures per walk, each 57277 Ir, a 2.69e-6 relative surcharge = 3.9e-6 bits. The escape row is now
measured, not modelled, and needs no reprice. Limit disclosed: the fixture integrates 50764 steps, so
it measures the per-step rate, never the cumulative regime — the cumulative charge is computed.

**Counter defect.** v10n's `escapes` field counted table-x mismatches only and never ring escapes, so
it printed 0 while the ring fired 124–162 times. The instrumented build reports both (83 hits ⇒ 83
escapes in one run).

**Coalescence.** The walk state is augmented (class × prev × ring); the solve is not. Measured at dp4
over 2000 trials: 62 class-meets, 0 with equal prev, 0 with equal ring, 0 augmented-meets, 62/62
verified. Augmented coalescence is ≤4.7 pct of class coincides at 95 pct confidence, and irrelevant,
because recovery is x = (a1−a2)/(c2−c1) with a full-point check that consumes neither prev nor the
ring. All 24 matrix cells: verified == meets exactly (115…900375 meets/cell) — witness-usable 1.000.

## 6. Scorable subset — the honest statement

The track's score is log2 of total classical work, admitted only with success ≥ 0.39 over uniform
production challenges. There is **no separately scorable work-only component**; v10n's attempt to bank
the work side while conceding the success side was refused, correctly. This package states the same
thing without routing around it: work exact and unchanged, success 59 bits short, **claim not
scoreable, not promotable**. It is filed as a testable correction of the record.

## 7. Work derivation (unchanged, re-audited)

Measured on the shipped counted TU `p256_recproc_v10l.o` (sha256 96b736ef…, 32864 B, spans bounded by
TEXT symbols only; live counters 0→64 `batch_inv` calls and 0→512 denominators on 512 walker-steps):

```
 239  step_forward
 694  step_close (44 mul, 54 umulh)          <- guard XOR + prev store + ring compare/store/index live here
 293  batch_inv longest back-substitution
  14  denominator zero-scan
0.8636  385 fe_mul x 147/65536
0.6398  rare-event reserve
0.00341 escape, 57277 Ir x 2^-24 (now measured at 2^-24.03)
-------
1241.5067493023873 Ir per walker-step
```

Σ rows equals the stated charge exactly. W_single = 127.631517 + log2(1.2415067) = 127.9436093073;
w=2 → 128.9436093073; claim 128.9437 is the 4-dp ceiling, so the at-most direction holds.

## 8. Not claimed

Production-length execution; any prior acceptance; a guarded plateau; any extrapolation beyond the
band of §4.1; that the guard can be repaired without changing the transition. The required next
construction is a coalescence-preserving guard — `prev` derived from a canonical function of the
point rather than carried history — and it is not in this package.
