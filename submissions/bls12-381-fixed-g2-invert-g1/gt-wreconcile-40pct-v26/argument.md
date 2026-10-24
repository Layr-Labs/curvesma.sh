# gt-wreconcile-40pct-v26 — argument

This package answers verdict `d2be26ee` (v25) on the four material flaws. Every
number below is emitted by `scripts/gen26.py` from `evidence/*.json`; the package
fails its own certificate if any of them disagrees.

## 1. The claim now covers W, by minimal ceiling

v25 published `claim = floor(log2(W/20000))` and certified `2^claim <= W/20000`.
That is a lower bound on the work presented as an upper bound on the score, which
is exactly the direction error the verdict names. v26 publishes the opposite:

    claim = min{ n/10000 : 20000 * 2^(n/10000) >= W }

decided in exact integer arithmetic, with no floating point in the decision:

    assert 20000^10000 * 2^n      >= W^10000     # covers W
    assert 20000^10000 * 2^(n-1)  <  W^10000     # minimal

Published claim: **126.5429** (numerator 1265431 / 10000) against

    W = 2 478 733 097 941 441 351 663 450 346 708 083 015 908 845   (log2 140.8306)

The nine ledger rows sum to that W exactly, and the audit perturbs **every** row
by +1 in turn and requires the sum to change each time. The claim is larger than
v25's 126.4981 because two rows grew when the dictionary and the executed trap
were priced honestly (sections 2 and 3); no number was tuned to a target.

## 2. The dictionary is implemented, tested and priced

v25 cited a "full-store dictionary" and priced `lookup_probe`: a 64-bit hash, two
probes, an `int` return, no insertion, no overflow rule, no load factor, and a key
builder `Group.state` did not use.

`scripts/ir26/dict26.c` ships the structure:
- `state_key48(const u64 v[12], u64 out[6])` — insertion-sorts the twelve field
  values `Group.state` sorts, then folds them pairwise into the 48-byte key, so the
  key is invariant to the rotation `kappa`.
- `dict_lookup(key, tab, mask)` — open addressing, **four** probes at stride
  `0x9E3779B97F4A7C15`, six-word compare on an occupied slot, stop at the first
  empty slot.
- `dict_insert(key, tab, mask, &evicted)` — first empty of the four probes, else
  **evict** the occupant with the smallest leading key word (ties to the lowest
  probe index) and count the eviction.
- `dict_lf_ok(n, cap)` — the load-factor guarantee the caller enforces: `8n <= 5cap`.
- Layout: eight 64-bit words per slot (six key words, a record identifier, a tag), 64 B/slot.

Differential test (`scripts/verify26.py`, `evidence/dict26.json`): 20 000 ops
(11 289 inserts, 8 711 lookups, 8 458 hits, 253 misses, 306 evictions) against a
Python model — **20/20 table digests match, 0 structure mismatches, 0 counter
mismatches**, `state_key48` matches on 32 tuples, and a mutant whose stride
constant differs by one is **detected** (20 digest diffs, 3 counter diffs).

Prices are `otool -tvV` spans of the compiled object, not notation:
`state_key48` 493 Ir, `dict_lookup` 153 Ir, `dict_insert` 246 Ir, `dict_record` 4 Ir
(`evidence/ir26.json`). The ledger therefore charges

    per step  : LOOKUP    = 493 + 153 = 646 Ir  (build the key, then probe)
    per record: REC_WRITE = 688 + 246 = 934 Ir  (serialize 688 B, then insert)

`evidence/ledger26.json → symbol_audit` lists every cited symbol with its span and
the evidence file that prices it; the certificate fails unless each one exists.
Nothing named in the ledger lacks a symbol in the tree.

## 3. Trap evidence: protocol traffic for coverage, executed steps for the mechanism

The verdict was right that v25 counted only slot words of lengths 4, 6, 8 and that
its forced-admission branch incremented a counter without running `Group.step`.

**Coverage (protocol-generated traffic only).** `evidence/scan26.json`: 1 920 000
steps of the shipped step on protocol-generated walks (fresh coins per segment,
`q=523`, `rp=521041`, `Rg=65536`, depth 5), counting balanced alternating slot
words of lengths **4, 6, 8, 10, 12**: 0 admissions, 0 aborts, max squaring run 9.
`evidence/scan25.json` adds 5 760 000 protocol steps with 0 admissions for 4/6/8.
The rule-of-three bound used in the waste row is `3/5 760 000` per step.

**Mechanism (existence proof, labelled as such).** `evidence/scan_forced26.json`
forces an admission every 1000 steps in a mutant build, and the forced branch now
**executes** the shipped transition: `y <- cyclo_sqr(y)`, `A,D <- 2A,2D`, `hist`
untouched — the same path a guard-blocked multiply takes. The data shows it:
`multiply_per_step` falls from 0.6394 (protocol) to 0.6022 (mutant) because the
trap really consumes steps as squarings. Delay to the 32-squaring abort: **mean
31.269, max 32 over 312 events**. The control build with the abort disabled runs
2 436 stall steps past 32 and never abandons (`evidence/scan_abortctl26.json`).
Waste charged:

    S x (3/5760000 x 32 + (1-f_mul)^32 x 32) x STEP
      = 29 082 097 875 473 205 031 598 916 382 482 736 638 Ir

## 4. Survival: measured at production Rg, extrapolation charged

The production arm is not extrapolated. `evidence/survival25.json` contains the
`Rg = 65536` arm at **400 pairs**, gap-1 survival **0.995**, for both `true_merge`
and `indep_hist`; the `samerep` control sits at 0.075, so the measurement
separates the arms. What remains unmeasured is the shape beyond the fitted range,
so v26 **charges the fit's residual envelope** rather than assuming it away:

    success_final = CP-95 x measured gap-1 survival / band_factor
                  = 0.533116 x 0.995 / 1.358246
                  = 0.39054                       (floor 0.39, met)

with `max |log residual| = 0.3062` over the seven fitted gaps reported alongside.
The uncharged point estimate (0.5305) is published only as a point estimate and is
not what enters success. Route shipped: **band**, applied to the measured
production-Rg arm.

## Cost model

    STEP = f_mul x MUL_CHARGE + (1-f_mul) x SQR_CHARGE
         = 0.639418 x 22409 + 0.360582 x 27611 = 24284.746371875 Ir
    MUL_CHARGE 22409, SQR_CHARGE 27611 : disassembled corrected cyclotomic_sqr
    f_mul 0.639418, guard rate 0.001480729 : measured over 1 920 000 protocol steps
    S = 829/572 x sqrt(N), N = (q-1)/12, 1024 workers, one record per visited class

Ten kernels (f2/f6/f12 mul, sqr, inv, cyclotomic_sqr) pass 10 000 bigint vectors
each with per-kernel mutation controls (`evidence/kernel_verify25.json`).

## Ledger (exact integers)

| row | Ir | derivation |
|---|---|---|
| walk_steps | 2 326 567 790 311 729 259 904 662 356 923 284 227 425 719 | S x STEP |
| per_step_lookup_key_and_probe | 61 889 169 832 220 689 915 444 980 164 352 536 868 140 | S x 646 (ir26) |
| record_writes_and_insert | 89 480 626 351 848 489 753 909 615 284 063 884 574 060 | S x 934 (ir26) |
| collision_processing | 8 192 | ceil(S(S-1)/2N) x 4096 |
| trap_waste_bound | 29 082 097 875 473 205 031 598 916 382 482 736 638 | rule-of-three x executed delay |
| index_slots_one_per_record | 766 429 347 767 438 884 401 795 419 991 981 880 720 | S slots x 8 words |
| worker_starts_1024 | 7 843 840 000 | 1024 x 383 muls |
| code_and_advice | 1 048 576 | 2^20 |
| candidate_verification | 57 526 800 | 1200 x f12mul span + 18000 compare |
| **W** | **2 478 733 097 941 441 351 663 450 346 708 083 015 908 845** | exact sum of the nine rows |

Peak storage is charged explicitly: pool `S x 688 B`, table `S x 64 B`, worker
state `1024 x 256 B`, history `1024 x 5 x 8 B`, code 2^20.

## Reproduce

```
python3 scripts/verify26.py
python3 scripts/protocol25.py '{"mode":"scan","Rg":65536,"depth":5,"q":523,"rp":521041,"segs":16,"steps":120000,"force_admit":false,"abort_enabled":true,"procs":16,"out":"scan26.json"}'
python3 scripts/protocol25.py '{"mode":"scan","Rg":65536,"depth":5,"q":523,"rp":521041,"segs":8,"steps":40000,"force_admit":true,"force_period":1000,"abort_enabled":true,"procs":8,"out":"scan_forced26.json"}'
python3 scripts/gen26.py
```

## Withdrawn from v25

The `floor(log2)` claim convention and its `2^claim <= W` certificate are
withdrawn: they certify the wrong direction. `lookup_probe` is withdrawn as the
priced per-step lookup; `dict_lookup` + `state_key48` replace it. The
"coverage" reading of the forced-admission run is withdrawn; it is now labelled a
mechanism existence proof, and coverage is claimed only from protocol traffic.
