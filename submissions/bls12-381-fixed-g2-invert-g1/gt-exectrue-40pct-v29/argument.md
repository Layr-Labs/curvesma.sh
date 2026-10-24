

## v29: the six defects of verdict 3e668bfb, fixed

1. **Wide storage access, not just wide addresses.** `dict28_wide_lookup`/`dict28_wide_insert`
   key a sparse store by the FULL 128-bit slot; the fixture round-trips 4096/4096 references at
   `cap = 2^70+12345` with 8062 slots above 2^64 and 0 false hits (`evidence/exec28.json`,
   `wide_access`). Every `(u64)slot` truncation is gone from the routines and from the Python
   call path, which now passes `tab_slots` as two words.
2. **Key width.** v28's `state_key48w` ended in `out[j] = mix64(s ^ c_j)` — `mix64` is a
   bijection, so all six words carried the same 64 bits: 2^64 distinct keys against ~2^126.2
   records, i.e. 2^188 expected false-key pairs, and the store could retain 2^64 references at
   most. v29 emits two INDEPENDENT 128-bit accumulators (`out[0..3]`) plus two cross terms, so
   the key carries 256 bits: 2^-3.6 expected false keys at production scale. The dictionary
   compares all six words and lookup precedes insert, so a reference is overwritten only after
   an equal key has already been reported.
3. **The executed transition is the charged one.** `run_one` charged a squaring whenever
   `br OR guard` held (`f_mul = 0.6433`), while `Group.step` squares only on a guard hit and
   performs one dense Fp12 multiplication otherwise. The charge is now
   `guard → cyclotomic square, every other step → dense multiply`:
   STEP = 0.99680 × 47 924 + 0.00320 × 25 372 = **47 851.7 Ir**, the number the verdict computed.
4. **The executed-count harness was RUN.** `exec28.py` now invokes its counting builds:
   `state_key48w` is called 20 000 times in a driver linked against the shipped object and its
   per-loop iteration counters are read out (`lane_cmp`, `lane_move`, `outer`) and checked
   against the analytic maxima (66, 66, 11); the wide-storage fixture runs and its JSON is
   embedded. The loop allowance now sums EVERY loop (outer 11, comparisons 66, lane moves 66,
   4 probes, 3-scan) instead of `MAXIT[s][0] × bodies[0]`. `f12mul` is priced as executed work
   with FOUR Fp6 multiplications = **72** base-field multiplications (per-Fp row 665.6, not
   887.48). Charges: `state_key48w` 2 926, `dict28_lookup` 759, `dict28_insert` 1 224.
5. **The eviction counter measures losses.** The old predicate (`evicted index >= new record
   index`) could not fire for a genuine old record. v29 tracks per-record match status and
   counts a record evicted BEFORE it produced a match: **1 792 of 15 261 evictions** are lost
   references, and the measured success contrast is 225/320 (evicting) vs 224/320 (drop), rate
   difference 0.003125 — the loss is now disclosed, not asserted away. Library counters are
   read per fixture as deltas, not accumulated.
6. **No rule of three.** The `rate = 3/windows` literal is gone. One inequality block: zero
   admissions over 4 940 101 windows gives a one-sided 95 % POISSON limit 2.996/N per window
   under a Poisson-arrival hypothesis (overlapping windows are dependent, so this is a stated
   hypothesis, not a theorem); correlation can overstate precision by at most the window
   overlap multiplicity, bounded by the ring (6), so the charged per-window rate is 6 × 2.996/N
   = 3.49e-6 and the union over windows is carried by charging every one of the S steps. Trap
   charge 5.34e38 Ir (10.4 % of W), not 0.18 %.

Census at depth 6 with the shipped `dict28.o` and the new key routine: A 225/320 (CP95 0.6583),
B 200/320 (0.5782), depth-5 contrast 88/128, partition series D/E/F 128/128. Bootstrap p5 =
0.5762 (median 0.6194) over min(A,B) × retention, floor 0.39.

## What moved

STEP 39 878.6 → 47 851.7 (+20 %), lookup row 2 830 → 3 685, record row 1 204 → 1 912, trap
7.4e37 → 5.34e38. W = 2^141.8777, claim 127.5900 — 0.28 bits worse than v28's 127.3064 and 1.87
bits worse than our own 125.72 ticket. The number is what the audit says; we ship the honest
one rather than keep a transition the reviewer has already refuted.
