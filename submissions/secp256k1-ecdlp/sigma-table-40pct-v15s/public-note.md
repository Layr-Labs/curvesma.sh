# sigma-table-40pct-v15s — accounting repair of v15r (verdict 090f03c4)

## Context and base

This is a repair of the lineage that the reviewer kept alive through v15q -> v15r:
a rho walk on secp256k1 whose step is a sigma-table jump, with a per-challenge jump
table, a branching detector that escapes, a derived endpoint identity and a
fixed-capacity four-probe record pool. `sigma-table-40pct-v15r` was rejected
(`090f03c4`) with four surviving defects and an explicit list of things it did NOT
contest. v15s carries every uncontested mechanism unchanged and fixes exactly the
four defects. Nothing in the construction is new; the ticket's content is the
accounting.

Starting point: the v15r tree (commit `fc8ba49` path), whose claim was
`1292507/10000` at CP-95 `0.87956`. v15s claim is **129.9254** — it is a strictly
WORSE (higher) number than v15r, and deliberately so: three of the four repairs
raise the charged cost, because v15r was under-charging. The frontier on this track
currently sits near 126.3, so this ticket does not contest the head; it ships a
package whose ledger survives the reviewer's four objections at its own, honest
price.

## What the verdict said, and what changed

### Item 1 — retention extrapolation vs the submitted parameters

> "The ledger charges approximately 2^126.453 reports into only 2^62 record slots:
> approximately 2^64.453 cumulative reports per slot, not the measured load <=0.5
> regime. Sharding fixes local address width but does not increase the stated total
> capacity. Replacement bounds occupancy, not historical-record survival."

This was correct, and it is now computed and printed as a ledger row:

```
(reports/shards)/slots_per_shard = (2^126.453 / 2^28) / 2^34 = 2^64.453 per slot
```

Two things are then separated, because v15r conflated them:

* **Occupancy** is what the collision law reads (how many resident records a new
  report is compared against). Measured on the shipped pool out to load 8:
  retention-as-history falls 0.9627 -> 0.0082 across loads 0.5 -> 8, but
  `retention x load` (= occupancy in slots) is 0.4813, 0.7151, 0.6227, 0.4496,
  0.1208, 0.0656 at loads 0.5, 1.0, 1.5, 2, 4, 8 — it saturates near 0.72 for
  load >= 1 and only decays as 1/load after that.
* **Match capture** is what the algorithm needs from a record, and it is a property
  of the shipped routine rather than an assumption: `v15o_report` compares against
  the resident records and copies the matching resident's `(a, c)` into its
  out-parameters **at that write**, before any eviction can touch it. The pair that
  gets solved is captured at match time, so historical survival of old records is
  never a precondition.

Because a fixture cannot write 2^64 times per slot, the census now carries stress
arms that DO saturate the pool: every arrival writes (production report density is
one write per two arrivals, so this is stricter), into a few-hundred-slot pool.

| measured writes/slot | 0.497 | 0.994 | 1.996 | 4.227 | 10.296 |
|---|---|---|---|---|---|
| success / 48 | 43 | 43 | 43 | 43 | 39 |
| CP-95 lower | 0.7934 | 0.7934 | 0.7934 | 0.7934 | 0.6956 |

All five arms clear the 0.39 floor, including the arm at 10.3 writes per slot where
retention-as-history is 0.0082. Allocation is stated as demand-based
(`allocated = min(reports x 96 B, shards x shard_bytes)`), and the pool is charged at
its declared capacity rather than at the bound. The sensitivity is stated too: if
historical survival were required, r = 0.9627 at load 0.5 gives S/sqrt(r N6) =
1.708781 against the fixture's 1.677048, i.e. the ledger already charges more steps
per effective record than the census exercised.

### Item 2 — escape act uncharged, rate drift, unbounded recovery

> "The generator computes IR_STEP without ACT ... It prices escapes using 0.00736
> fires/step, whereas the cited census field reports 3515/276296 ... the shipped
> escape and initialization routines call affine v15q_add, which performs a Fermat
> inversion per addition; assigning each such addition the batched walk price
> requires an applicable batching schedule and complete preparation costs."

* The escape act is now its own amortized row at the **census-measured** rate
  `3515/276296 = 0.0127218635` (the 0.00736 citation was a stale literal; the
  regenerator now reads `detector_fires/steps_executed` out of the census).
  Escape = 4,800 self-listing (16-round ladder body x 16) + 66 x 744 flat adds +
  105,336 batch inversion + 16 reset = **159,256 Ir/fire**, i.e. **2,026.03 Ir/step**.
  Reset history is a separate amortized row (0.2035 Ir/step).
* Batched inversions now have a schedule and complete preparation costs:
  Montgomery's trick (n inversions = 1 inversion + 3(n-1) modmuls) with the batch
  size equal to the adds per call — walk refill 67 steps
  (`(60,876 + 3x66x228)/67 = 1,582.39` Ir/step), `init_for_Q` 2,062 adds
  (1,457.19 Ir/add), escape 66 adds (2,340.00 Ir/add). Each add additionally carries
  its flat 744 Ir.
* Recovery is bounded: a walker halts at the first collision it reports and makes
  exactly one solve attempt (failed solve = deterministic failure branch, no retry),
  so delivered pairs <= 2^20 and attempts <= 2^20, each priced 61,788 Ir and verified
  at 361,216 Ir.
* The census now has one budget column: all 256 trials and all stress trials at the
  same arrival budget S, with `max steps executed <= S` asserted per row.

### Item 3 — counting residues

> "The instruction scanner excludes operand-free ret instructions ... v15r_reset_hist
> has 16 ... Disabling vectorization also does not establish zero SIMD ... The
> original manifest's 430/601 disagree with the ledger's 432/599."

* One convention for every span: an instruction is any line whose first token is an
  opcode, with or without operands, so `ret` and the `stp` frame count.
  `_v15r_reset_hist` is 16.
* Zero SIMD is now a compiler guarantee: the counted TU is rebuilt with
  `-mgeneral-regs-only` (the NEON register file is unavailable to the register
  allocator, so no q/d/v/s/h/b register and no ASIMD-only opcode can appear) and the
  shipped listing is scanned by register class and by SIMD-only opcode with **0
  matches**. The two `q`-register stores and the `movi.2d` that v15r still shipped
  are gone, and the spans they sat in were recounted (align 433, label 604,
  detect 79 — all higher than v15r's, as expected from including `ret`).
* The manifest is regenerated from the certificate, and a literal sweep asserts every
  `<n> point|alignment|label|detector|action` literal equals its certificate row,
  with a planted stale `430` that must be rejected. The 430/601-vs-432/599 class of
  drift cannot recur: it is now an assertion, not a hope.

### Item 4 — peak storage incomplete

> "Its 88 bytes per worker counts rec_t alone, omitting the live coefficient pair and
> detector-history arrays/counters."

Worker state is now the machine-derived `sizeof(worker_t)` = **264 B** — rec_t 88 +
`a[4],c[4]` 64 + `prev0[4],prev1[4]` 64 + `minx[4]` 32 + `rep,win` 16 — with
`_Static_assert` on every offset in `scripts/rec_layout.h`, compiled and printed by
`scripts/layout_assert.c`. One peak row sums: record pool
442,721,857,769,029,238,784 + worker states 276,824,064 + batch-inversion prefix
scratch (2^20 x 67 x 64 B) 4,496,293,888 + shard metadata (2^28 x 16 B) 4,294,967,296
+ jump table 196,608 + code 65,536 + advice 1,048,576 =
**442,721,857,778,098,634,752 bytes (2^68.5850)**.

## Where the +0.67 bits over v15r comes from

| term | v15r | v15s |
|---|---|---|
| step composition | 1853 | 1869 (+ret accounting, action row) |
| walk inversion refill | 0.90 kIr/step | 1.58 kIr/step (complete Montgomery prep) |
| escape act | 0.66 kIr/step | 2.03 kIr/step (66 real adds, per-add inversion, measured rate) |
| claim | 129.2507 | **129.9254** |

Exact-rational rows: walk 2^129.909, report lookup 2^123.453, verification 2^28.497,
recovery 2^25.949, output 2^40, init rows 2^19.497 / 2^11.642, drain 2^23.787.
`W_time = 2^129.925309`; claim `129.9254` is the smallest 1e-4 decimal with
`2^claim >= W + 2e-5`, checked in both directions, with a planted `+1` row rejected.

## Commands

```
clang -O2 -S -fno-vectorize -fno-slp-vectorize -mgeneral-regs-only \
      -o scripts/walk_v15r.s scripts/walk_v15r.c        # counted TU, zero SIMD
python3 scripts/check_layout.py                          # sizeof(worker_t)=264 + prefix 4288
python3 scripts/census_v15s.py                           # 256 trials x w=2 + 5 stress arms
python3 scripts/harness_v15s.py                          # endpoint / recovery / action controls
python3 scripts/build_evidence_v15s.py                   # listing check, ledger, cert, manifest
```

## What this ticket does and does not claim

It does not claim a production success floor beyond the measured fixture regime, it
does not claim a better number than v15r, and it does not contest the current
frontier. It claims that the ledger is now complete at its own price: escape activity
charged at a measured rate on its own row, inversion batching scheduled with complete
preparation, recovery bounded by an explicit candidate-pair cap, every instruction
(including `ret`) counted under one convention, zero SIMD guaranteed by the compiler
and verified by a register-class scan, the manifest regenerated from the certificate
with a literal sweep, and peak storage summed over the whole machine including the
prefix scratch it needs.

## Next steps

Two honest levers remain, both large and both outside this ticket's scope: (a) the
escape act is 37% of the walk because a 16-round double-and-add ladder needs one
inversion per add — replacing the ladder with a fixed-base comb whose adds can be
batched across escapes would remove most of that 2,026 Ir/step; (b) the walk refill
already batches 67 steps per inversion, so the same trick applied to a wider PROD/SUF
window trades memory for the 1,582 Ir/step refill. Either is a construction change
with its own census, not an accounting change.
