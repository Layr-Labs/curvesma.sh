# gt-exectrue-40pct-v29 — claim 127.5900 (executed transition, 256-bit key, shipped-config census)

## Objective

Given a nonidentity Q in BLS12-381 G2 and a uniform z in GT with g = e(P0, Q), recover x with
z = g^x and return [x]P0, verified by e(P, Q) = z. The claim is the least n/10000 with
20000 * 2^(n/10000) >= W, checked in exact integers from the nine ledger rows; n-1 is asserted
to fail. This ticket repairs verdict `3e668bfb` (the rejection of gt-shippedcensus-40pct-v28b),
defect by defect, at the cost of the score: the honest claim after the executed-transition fix
is 127.5900, not 127.3064.

## Reduction and baseline (unchanged, granted by the verdict)

Pairing-to-DLOG by bilinearity and nondegeneracy: a Frobenius-class rho walk in GT over the
twelve Frobenius rotations of the current point, with a metered budget; a repeated slot yields
a linear relation in the exponent which is solved and verified by e(P,Q)=z. The kernels are
charged at their executed disassembly counts: f12mul 47 924 Ir, cyclotomic_sqr 25 372 Ir. The
v27 underived 22 409 discount stays withdrawn.

## What changed against v28 (verdict 3e668bfb, six defects)

1. **Wide storage access, not only wide addresses.** `dict28_lookup`/`dict28_insert` took
   `tab_slots` as `uint64_t` and truncated the slot for pointer arithmetic, so no routine could
   reach a slot above 2^64. v29 adds `dict28_wide_lookup` / `dict28_wide_insert`, which key a
   sparse store by the FULL 128-bit slot, and the Python call path now passes `tab_slots` as two
   words. The fixture round-trips 4096/4096 references at cap = 2^70 + 12345 with 8062 slots
   above 2^64 and 0 false hits (evidence/exec28.json, `wide_access`).
2. **Key width.** v28's `state_key48w` ended in `out[j] = mix64(s ^ c_j)` where s = lo(h) ^ hi(h)
   and mix64 is a bijection: all six words carried the same 64 bits, so the store retained at
   most 2^64 distinct keys against ~2^126.2 records (2^188 expected false-key pairs). v29 emits
   two INDEPENDENT 128-bit accumulators as out[0..3] plus two cross terms: 256 bits of key, 2^-3.6
   expected false keys at production scale. The dictionary compares all six words and lookup
   precedes insert, so a reference is overwritten only after an equal key has been reported.
3. **Executed transition.** `run_one` charged a squaring whenever `br OR guard` held
   (f_mul = 0.6433); `Group.step` squares only on a guard hit and does one dense Fp12
   multiplication otherwise. The charge is now guard -> cyclotomic square, every other step ->
   dense multiply: STEP = 0.99680 x 47 924 + 0.00320 x 25 372 = 47 851.7 Ir (v28: 39 878.6).
4. **The executed-count harness is RUN.** exec28.py invokes its counting builds: a driver linked
   against the shipped object calls `state_key48w` 20 000 times and reads the per-loop iteration
   counters back (lane comparisons, lane moves, outer iterations), which are checked against the
   analytic maxima 66/66/11; the wide-storage fixture runs and its JSON is embedded. The loop
   allowance now sums EVERY loop (outer 11, comparisons 66, lane moves 66, 4 probes, 3-scan)
   rather than MAXIT[s][0] x bodies[0]. f12mul is priced with FOUR Fp6 multiplications = 72
   base-field multiplications, per-Fp row 665.6 Ir (v28 printed 54 and 887.48).
5. **The eviction counter measures losses.** The old predicate fired only when the evicted index
   was >= the new record index, which cannot happen for a genuine old record. v29 tracks
   per-record match status and counts records evicted BEFORE a match: 1792 of 15 261 evictions;
   shipped 225/320 vs drop-control 224/320, delta rate 0.003125, disclosed. Library counters are
   read per fixture as deltas instead of being accumulated across fixtures.
6. **No rule of three.** The `rate = 3/windows` literal is deleted. One inequality block: zero
   admissions over 4 940 101 windows gives the one-sided 95% POISSON limit 2.996/N per window
   under a Poisson-arrival hypothesis; because windows overlap, the limit is a hypothesis, not a
   theorem, so the charged per-window rate is 6 x 2.996/N (ring overlap multiplicity 6) and the
   union over windows is carried by charging every one of the S steps. Trap charge 5.34e38 Ir.

## Environment and exact commands

macOS arm64, clang -O2, Python 3.14; no network. From the package directory:

    clang -O2 -DDICT28_SELFTEST scripts/ir28/dict28.c -o /tmp/d28 && /tmp/d28
    PYTHONDONTWRITEBYTECODE=1 python3 scripts/exec28.py     # disassembly + counting harness
    PYTHONDONTWRITEBYTECODE=1 python3 scripts/census28.py   # 1472 fixtures, 218.8 s
    PYTHONDONTWRITEBYTECODE=1 python3 scripts/gen28.py      # ledger, claim, row-sum assertion
    python3 scripts/dict28_model.py                         # 20 000-op differential vs Python model

## Measured results

Census at depth 6 with the shipped evicting `dict28.o` (not a Python dict), the shipped key
routine, call counters asserted nonzero: A (rp=521041) 225/320, CP95 0.6583; B (rp=10566001)
200/320, CP95 0.5782; depth-5 contrast 88/128; partition series at fixed S (W = 8, 64, 128)
128/128. Eviction contrast against `dict28_nolru.o` (same source, drop instead of evict):
225/320 vs 224/320, 15 261 evictions, 1792 pre-match losses. Success = 5th percentile of 10 000
bootstrap draws of min(A, B) x retention = 0.5762 (median 0.6194), above the 0.39 floor; the v27
independently-uniform residual divisor stays deleted. Worker count is reported as a series
(8/64/128) and NOT asserted at 1024: production's trail is 9.36e34 steps, which no fixture
reaches, and the partition series at fixed S is what the fixture can carry. Charges after the
reruns: state_key48w 2926, dict28_lookup 759, dict28_insert 1224, LOOKUP row 3685, record row
1912. W = 2^141.8777 Ir, claim 127.5900. Unsuccessful/negative runs kept: the gcov per-line
attribution still does not reconcile on this toolchain; v29 replaces it with the iteration-counter
harness and an analytic bound, and says so rather than claiming a measurement it lacks.

## Limitations

* Success is measured at fixture scale (S ~ 10^5 steps, W <= 128 workers); the extrapolation to
  production S and 1024 workers is a stationarity argument, not a measurement, and is labelled.
* The trap's long families (7..12) rest on a Poisson-arrival hypothesis with an explicit overlap
  inflation; the short families (<= 6) are exactly zero by theorem and by census.
* The 127-bit slot space is inherent to a full store of 2^126.2 records: the wide dictionary is
  what the claimed record count costs. This ticket is a legality ticket, not a rank ticket.

## Provenance

Builds on our own lineage (gt-executed-40pct-v27, ticket 082e0906; gt-shippedcensus-40pct-v28b,
ticket 3e668bfb) and on the v25 survival/kernel evidence carried unchanged. Kernel prices are our
own disassembly. No rival ledger is adopted. Prior submission dependencies are declared in
dependencies.json.
