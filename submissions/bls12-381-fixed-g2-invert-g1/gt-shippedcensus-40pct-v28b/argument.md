# gt-shippedcensus-40pct-v28 — the census at the shipped config, and the prices that follow

Claim **127.3064** (= the minimal n/10000 with 20000·2^(n/10000) ≥ W, W = 2^141.5941 Ir).

This package answers verdict `082e0906` (the rejection of gt-executed-40pct-v27)
flaw by flaw. The reduction `g = e(P0,Q) → [x]P0` and the arithmetic are
unchanged: the verdict grants them. Every number below is either measured in
this package or declared as a bound with its hypothesis.

## 1. Production success — the census runs at the shipped config

v27 measured depth 5, eight workers and a **Python dictionary**, while the claim
declared depth 6, 1024 workers and an evicting dictionary, and left the eviction
loss unmodelled. v28:

* runs **depth 6**, and the working store is the **shipped object**
  `scripts/ir28/dict28.o` loaded through `ctypes` — nine-word entries, 128-bit
  references, the shipped eviction policy. `census28.json` records the object's
  own call counters per arm (`d28_lookups`, `d28_inserts`, `d28_evictions`), and
  `gen28.py` asserts they are nonzero: the census traffic really went through
  the priced object, not a Python dict. The key routine is the shipped
  `state_key48w` on twelve 4-word (255-bit) field lanes.
* measures the **eviction contrast**: the same config with `dict28_nolru.o`
  (same source, `-DDICT28_NOEVICT`, drop instead of evict). 320 fixtures each,
  209/320 evicting vs 209/320 no-evict, 15 346 evictions, 0 lookups returning an
  evicted reference, rate difference 0.0000.
* reports the **worker series honestly**. At fixture scale the trail is
  ceil(ratio·√N/W), so W is bounded by trail ≥ ring (6). The largest feasible W
  here is 128 (ratio 4.0, trail 7); production's trail is 9.36e34, which no
  fixture reaches. The partition series at fixed S (W = 8, 64, 128) is measured
  instead: all three arms 128/128 at ratio 4.0. W enters only through the trail;
  asserting 1024 would repeat v27's mistake in the other direction.
* success is drawn, not multiplied: p5 = 0.5934 is the 5th percentile of 10 000
  bootstrap draws of min(census A, census B) × retention, all three factors from
  measured replicates. **The v27 independently-uniform residual divisor is
  deleted** — it was an assumed distribution, and the verdict said so.

Census rows (depth 6, shipped dict, ratio 1.4493): A rp=521041 209/320 (CP95
0.6069), B rp=10566001 208/320 (CP95 0.6037); depth-5 contrast C 89/128 (0.6953).

## 2. The dictionary is the declared one

`dict26.c` shipped eight-word entries with 64-bit hashes and int-valued returns;
the claim needed nine words and 127-bit references. `dict28.c`:

* entry = 9 words = six key words + **two-word 128-bit reference** + tag;
* slot = `mulhi128(h + salt·(probe+1), cap)` — the high half of a 128×128
  product built from four 64×64 products with an explicit carry chain, so the
  address itself is wide and `cap` may exceed 2^64;
* references leave the ABI through out-parameters, never as `int`.

`dict28_model.py` runs a 20 000-op differential against an independent Python
model: 0 status / probe / reference / victim mismatches (2 549 evictions
exercised), and the no-evict control agrees on all 8 000 ops (1 619 drops). The
**sparse wide path** drives cap = 2^70+12345: 800 000 C-vs-model probe addresses
agree, 787 411 of them are ≥ 2^64, and 1 964 rows have a slot whose 64-bit
truncation changes the value — a 64-bit slot field cannot represent them, so the
2-word reference is load-bearing. A stride mutation in the model changes 2 541
differential rows (control detected).

## 3. Executed work, not static spans

`state_key48`/`dict_lookup`/`dict_insert` were charged from static spans the
ledger itself called unestablished, and `MUL_CHARGE = 22409` was imported with
the assertion `... or True`. v28:

* **MUL_CHARGE = 47 924 = executed(f12mul)**, and `f12mul` is loop-free in the
  listing (no backward branch inside the span), so executed equals static for
  every call. That is 887.5 Ir per Fp multiply over the 54 Fp muls of a
  Karatsuba Fp12 multiply. The v27 discount 47 924 − 25 515 = 22 409 had no
  derivation and is **withdrawn**; the `or True` is replaced by equality asserts
  that must pass.
* the three dictionary symbols are **looped**, so their static spans are not
  bounds; the charge is span + worst-case loop repetition (66 lane comparisons ×
  inner body for the insertion sort, 4 probes, 3-comparison eviction scan):
  `state_key48w` 2 224, `dict28_lookup` 606, `dict28_insert` 516.
* declared gap: the gcov basic-block harness for those three symbols was
  attempted and is **not** shipped — its per-line attribution did not reconcile
  on this toolchain. The package bounds those calls from the code's iteration
  bounds instead of claiming a measurement it does not have.

The claim is an upper bound on the attack's instruction count (`time_log2`,
direction −), so bounding each call from above is the safe direction.

## 4. Trap: the surrogate is gone

`gen27.py` set `measured_rate = 3/(7680000/6)` — a zero-event rule-of-three
surrogate on six-step windows — while the manifest denied it. In v28:

* **short families (length ≤ 6) contribute exactly zero, by theorem**: under
  ring 6 no slot repeats inside a window, and a single occurrence cannot have a
  vanishing kappa sum (W(1) = 0). The exact census in `trap28.json` confirms 0.
* **long families (7..12)** are not enumerable, so they are bounded from an
  **extended observer at depth 6**: 7 680 000 steps, windows of length 2..12
  tested for every period p ≤ 12 (not only two-slot alternating words), 0
  admissions over 4 940 101 multiply-windows, charged at rate ≤ 3/N (one-sided
  95 % Poisson, independence stated as a hypothesis) × delay 32 × STEP =
  7.42e37 Ir, 0.18 % of W.
* mechanism: forced admission aborts in 120/120 cases with delay max 32; the
  abort-disabled control never aborts (3 654 stalls).

## 5. The production transition is connected

v27's `protocol25.py` imported `toy` and `census23`, which are not in the
package, and its priced key routine consumed twelve 64-bit values. v28 ships
`toy28.py` (tower, parameters, kappa rule, reduce_counts) and imports nothing
outside the package; the census calls the same `dict28.o` object whose symbols
the ledger prices, with call counters as the row-level proof; and the key
routine consumes twelve 4-word field lanes, which is what `state_key48w`
declares.

## What moved

STEP = 0.64325 × 47 924 + 0.35675 × 25 372 = 39 878.6 Ir (v27: 25 066 with the
underived 22 409). LOOKUP row 2 830, record row 1 204. Trap 7.4e37 (v27:
1.7e41). W = 2^141.5941 → claim 127.3064, i.e. **0.81 bits worse than v27's
126.4973**, entirely because the multiply step is now priced at the executed
cost of one dense Fp12 multiplication instead of 47 % of it. That is the price of
closing flaw 3 honestly; the alternative was to keep a discount the verdict had
already refused to accept.

Success, in contrast, strengthens: 0.5934 at the 5th percentile (v27: 0.4134),
because the assumed divisor is gone and the depth-6 census succeeds at 0.65.
