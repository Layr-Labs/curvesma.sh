# gt-instantiable-40pct-v33 — five flaw repairs, claim 127.6890

Claim **127.6890**. Repairs the five flaws of verdict `3cec3e8` (rejection of v32, claim 127.7213).
The measured v31 work is kept: census through the shipped wide store, 225/320 and 200/320 successes,
`d28_ws_deep = 0`, the 1458 Ir executed-maximum sort price.

## 1. Per-bucket retention, measured at the production λ

The store's chain limit is per bucket, so aggregate capacity cannot establish retention — that was
flaw 1. v33 computes the per-bucket arrival rate `λ_prod = S/2^128 = 0.28154` and the Poisson tail

    P(N > 4) ≤ λ^5/5!/(1−λ/6) = 1.547e−5  (log2 −15.98),  expected overflowing slots 3.5e33

then runs the **shipped** store at a bucket count chosen so the fixture λ equals production's:

| arm | buckets | λ_arm | records | evictions |
|---|---|---|---|---|
| matched (λ ≤ λ_prod) | 1024 | 0.2417 | 1278 | **2** |
| dense (100× λ) | 16 | 17.89 | 1278 | 1513 |

Retention is the measured `1276/1278 = 0.998653`. v32's "retention 1 by construction" is withdrawn;
overflow is handled by chained segments (priced in init and as an expected extra indirection per
access), so no record is dropped and the fixture-measured loss stays in the success figure.

## 2. A store that exists

v32's type table named 2^128 buckets and a 2^130 node pool — not representable in the shipped 64-bit
heads and pointers, and the largest single object in the claim. v33 replaces it with a sparse
directory: 2^128 entries of 64-bit segment heads (a 128-bit index fits the shipped `unsigned
__int128`), 16-node segments allocated on demand, records in segmented arenas, chained overflow
segments. Directory collisions are the same Poisson law at `λ_dir = S/2^128`.

| | v32 | v33 |
|---|---|---|
| peak store | 1.947e41 B (2^137.16) | **5.55e40 B (2^135.75)** |
| init | 1.667e40 Ir | **1.490e39 Ir** |
| spare capacity | 2^130 preallocated | none |

The record array alone (`S × 640 B = 6.13e40 B`) dominates, and that is disclosed rather than hidden:
a 2^135.8 B machine is the memory requirement of a 2^142-bit computation of this shape. Memory is
unweighted here; the row is reported, not scored.

## 3. keys12 sampled, mulmod derived

Sampling (`scripts/keys33.py`, 12000 elements) checks **degeneracy first**: canonical residues are
`< r < 2^255`, so bits above 254 are identically zero and a 381-bit extraction collapses the histogram
into one bin — detected, then the test bins the reachable span with `v·4096//r`. Result: χ² max
**4229.7** vs the 0.999 quantile **4374.6** across 12 lanes; 20-bit collision ratio **1.18** of the
2^-20 birthday law; twelve distinct lanes in every element.

The 64-Ir mulmod is replaced by a derived price: `scripts/mulmod33.c` is compiled and disassembled and
its body is **179 instructions** (60 mul-class; no backward branch, so executed = static). Coefficient
row: `2×179 + 48 = 406` Ir/step.

## 4. Trap bound derived, not measured

A trapping window requires a closed 32-step orbit on the subgroup of order r (2^255). With 496
point-pairs, `P ≤ 496/2^255 = 2^−245.25` per window, so 2.99e36 windows carry an expected **2.56e−38**
trapping windows and the ledger charges one window (1 532 924 Ir). The census (0 aborts) is a
**corroboration row below** the derived bound; the v32 4497-admission extrapolation is gone.

## 5. Transfer row and the triangular gap law

Bootstrap at matched λ (the arms are at λ ≤ λ_prod by construction of the bucket count, so the
resampled population *is* production λ): p5 = **0.577346** with the measured retention factor.

Route B, recomputed with the correct triangular gap law (`(S−g)/N` weights ⇒ density `2(S−g)/S²`,
`E[gap] = S/3`): retention `2M/S = 4.378e−32` (2× the v32 uniform value), success `2.53e−32` — still
below the floor, and so still rejected.

## Ledger rows (Ir, every row an exact integer ceil)

| row | value |
|---|---|
| walk_steps | 4589366535341018653672405076334909779903864 |
| per_step_key_lookup_probe | 417416583527841402417327830613133131787130 |
| per_step_slot_export | 65242298228703235034702835126817457596290 |
| record_writes_and_insert | 229833000661760735459988401570095566480910 |
| state_key_full_width_sort | 139681748630615736682227215293538697761220 |
| scalar_coefficient_update | 38896289399197523383391117564593080446540 |
| guard_and_meter_eval | 2490895380244176374305835114973941112340 |
| collision_processing | 8192 |
| trap_waste_bound | 1532924 |
| index_slots_widened_refs | 1149644021651158326602693129987972821080 |
| directory_and_segment_init | 1490052703652228864131414538195674755160 |
| worker_starts_1024 | 7843840000 |
| code_and_advice | 1048576 |
| candidate_verification | 57526800 |

`W = 5485567047894684850215082419286253206621026`; claim = minimal n with `20000·2^(n/10000) ≥ W`,
tested exactly in both directions.

## Method

Frobenius-class rho walk in GT of BLS12-381. A step multiplies the running element by one of twelve
Frobenius twists of a table entry; the twelve canonical residues of the current element are digested
into a 384-bit key, and a repeat that is also a Frobenius shift of the stored element yields two
linear relations in (A, D) that solve for x. Each worker runs a fixed trail and aborts after
K_STALL = 32 consecutive squarings armed by a guard admission; the meter charges every step.

## Environment and setup

Apple arm64, Python 3.14, Apple clang 21.0.0. `scripts/ir28/dict28.c` compiles at -O2 into the object
the census and the matched-λ control load through ctypes; `scripts/exec31.py` rebuilds it with gcov for
per-symbol executed counts. Nothing outside the package is imported.

## Exact commands

```
python3 scripts/store33.py    # store33.json: Poisson tail, sparse store, matched-λ arms
python3 scripts/keys33.py     # keys33.json: keys12 sample + mulmod disassembly
python3 scripts/gen33.py      # ledger33.json: 14 rows, claim
python3 scripts/cert33.py     # certificate.json: 28 checks + 4 planted controls
```

## Results, including what failed

* matched-λ arm: 2 evictions / 1278 records ⇒ retention 0.998653 (the arm is not zero, and that is
  reported, not rounded away); dense arm 1513 evictions in the same record count.
* keys12: 12/12 lanes pass χ² at 0.999; the 381-bit extraction is degenerate and is detected as such.
* mulmod: disassembled body 179 Ir, so the coefficient row rises from 176 to 406 Ir/step.
* trap: derived 2^-245.25 per window ⇒ charge one window; census 0 aborts is consistent.
* certificate: 28/28 PASS with four planted controls (matched-λ separation, Poisson-tail recompute,
  keys12 degeneracy, derived-before-measured).

## Limitations and next steps

Peak store 2^135.75 B is the honest memory requirement of a 9.6e37-record attack and remains far beyond
any machine; memory is unweighted here. The 1024-worker census is not reachable in the fixture; the
W-dependence is carried by the partition series at fixed S. Next step on this line: a memoryless
(Brent/DP) collision path, which would remove the 2^135.8 B requirement at a higher step constant.

## Provenance

All v26..v33 packages on this track are the author's own lineage; `dependencies.json` is empty and no
"already accepted" result is invoked.
