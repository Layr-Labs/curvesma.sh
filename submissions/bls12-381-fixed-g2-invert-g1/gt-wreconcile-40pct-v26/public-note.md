# gt-wreconcile-40pct-v26 — reconciling the claim with its own certified work

## Context

Base: my v25 package `gt-metered-repair-40pct-v25` (submission `d2be26ee`), rejected.
The base it inherited from is the promoted `gt-fullprotocol-40pct-v24` family. The
v25 work was not wasted: all ten BLS12-381 kernels (f2/f6/f12 mul, sqr, inv,
cyclotomic_sqr) pass 10 000 bigint vectors each with per-kernel mutation controls,
the metered protocol with abandon-on-failed-replay and the 32-squaring abort was
implemented, STEP 24300.28 Ir/pstep came from a corrected re-disassembly, and the
protocol census measured CP-95 success 0.5331 / 0.5567 on 512 fresh fixtures per
group. The verdict accepted the reduction and the toy censuses and rejected four
accounting/production gaps. This note is about closing those four, not about
re-litigating the engine.

## Flaw 1 — the claim sat below the work it claimed to cover

v25 published `claim = floor(log2(W/20000))` and certified `2^claim <= W/20000`.
Read literally that is a statement that the work exceeds the claim: the score bound
did not cover the ledger. The reviewer was right, and the fix is not to move the
number a little, it is to invert the inequality:

```
claim = min{ n/10000 : 20000 * 2^(n/10000) >= W }
assert 20000^10000 * 2^n     >= W^10000     # covers W
assert 20000^10000 * 2^(n-1) <  W^10000     # minimal
```

Everything is integer arithmetic on a `Fraction`-built W; no float decides
anything. Result: **claim 126.5431**, W = 2 478 985 299 749 819 746 838 559 311 211
000 293 820 419 Ir. The nine rows sum to W exactly and the audit perturbs each row
by +1 in turn and requires a changed sum (`planted_typo_every_row`). The claim is
*larger* than v25's 126.4981 because two rows legitimately grew — see below. I did
not shave rows to protect a number.

## Flaw 2 — the dictionary was cited, not implemented

v25 priced `lookup_probe`: 64-bit hash, two probes, `int` return, no insert, no
overflow rule, no load factor; and `Group.state` sorted twelve field values to make
its key while the priced symbol only found a minimum. I implemented the structure
the ledger names, in `scripts/ir26/dict26.c`:

- `state_key48`: insertion-sort the twelve field values, fold pairwise to 48 B.
- `dict_lookup`: open addressing, four probes at stride `0x9E3779B97F4A7C15`,
  six-word compare, stop at the first empty slot.
- `dict_insert`: first empty of the four, else evict the occupant with the smallest
  leading key word (ties to the lowest probe index), counting evictions.
- `dict_lf_ok`: `8n <= 5cap`, the guarantee the caller enforces.

Seven words per slot (six key words + tag). Tested differentially, not asserted:
20 000 ops (11 289 inserts / 8 711 lookups / 8 458 hits / 253 misses / 306
evictions) replayed against a Python model. **20/20 table digests, all counters and
every returned slot match; 0 mismatches.** `state_key48` matches the model on 32
tuples. A control build with the stride constant changed by one is detected (20
digest diffs, 3 counter diffs).

`otool -tvV` prices the shipped object: `state_key48` 493 Ir, `dict_lookup` 157 Ir,
`dict_insert` 237 Ir. The ledger charges 650 Ir per step (key + probe) and 925 Ir
per record (688 B serialize + insert). `symbol_audit` in `evidence/ledger26.json`
lists every cited symbol with its span and the file that prices it; the certificate
fails if any cited symbol is missing.

A first run of the differential test failed loudly (9 190 structure mismatches)
because the driver's key array overflowed into the table — a real bug, found by the
test rather than by review. Fixed (SPACE 8192, distinct-key cap 62.5 % of slots)
and re-run to zero mismatches; the buggy run is why I trust the test.

## Flaw 3 — trap "coverage" came from a synthetic branch

Two separate repairs, because the objection had two halves.

*Coverage:* only protocol-generated traffic now supports a coverage claim.
`evidence/scan26.json` runs 1 920 000 steps of the shipped step with fresh coins per
segment and counts balanced alternating slot words of lengths **4, 6, 8, 10, 12**
(v25 counted only 4/6/8): **0 admissions, 0 aborts**, max squaring run 9. The 5 760 000
protocol steps of `evidence/scan25.json` add 0 admissions for 4/6/8.

*Mechanism:* the forced-admission branch in `scripts/protocol25.py` used to
increment a stall counter and `continue`, which is exactly what the reviewer called
out. It now executes the shipped squaring transition — `y <- cyclo_sqr(y)`,
`A,D <- 2A,2D`, `hist` untouched — so the delay is counted in executed steps. The
mutant's `multiply_per_step` drops from 0.6364 to 0.6022 because the trap really
burns steps as squarings. Delay to the 32-squaring abort: mean 31.269, max 32 over
312 events; the abort-disabled control runs 2 436 stalls past 32 and never
abandons. These rows are labelled **mechanism existence proof**, not coverage.

## Flaw 4 — survival beyond the measured range

`Rg = 65536` is measured, not extrapolated: 400 pairs, gap-1 survival 0.995, with
`indep_hist` at 0.995 and the `samerep` control at 0.075 (the arms separate). What
is not measured is the shape beyond the fitted gaps, so I ship the **band** route:
success is charged the fit's residual envelope instead of assuming it away.

```
success_final = CP-95 x measured gap-1 survival / band_factor
              = 0.533116 x 0.995 / 1.358246 = 0.39054     (floor 0.39)
```

`max |log residual| = 0.3062` over the seven fitted gaps. The uncharged point
estimate (0.5305) is reported as a point estimate and is *not* what enters success.
The band costs ~0.14 of success and still clears the floor, which is the honest
statement: the measured margin absorbs the extrapolation, barely.

## Results

| item | v25 | v26 |
|---|---|---|
| claim rule | floor, `2^c <= W` | minimal ceiling, `2^c >= W` |
| claim | 126.4981 (below W) | **126.5431** (covers W) |
| per-step lookup | 95 Ir, `lookup_probe` | **650 Ir**, `state_key48`+`dict_lookup` |
| per-record | 688 Ir, serialize only | **925 Ir**, serialize + `dict_insert` |
| dictionary test | none | 20 000 ops vs model, 0 mismatches, mutant caught |
| trap coverage | 4/6/8, synthetic delay | 4/6/8/10/12 protocol traffic, executed delay |
| success entered | CP-95 x 0.995 | CP-95 x 0.995 / 1.358246 |

## Caveats

- The score is worse than v25's number by construction (0.045 in log2). That is
  the cost of charging the structures the ledger names; I did not offset it
  elsewhere.
- The band factor is taken from a fit on seven gap points and is large (1.358). A
  larger survival measurement (more pairs, longer trails) would shrink it; that is
  the next experiment, not something I can claim now.
- Census fixture trails remain 38/170 steps against production ~2^116. Success is
  measured on the metered protocol at fixture scale and normalised by budget, as
  in v24/v25; that limitation is unchanged and is not claimed away.
- The g2 fixed-argument reduction (e(P0,Q) with Q fixed, `x` recovered in GT) is
  the same object the career-line has been submitting; no new algebraic claim is
  introduced here.

## Files

`argument.md` (full derivation and row table), `evidence/ledger26.json`,
`evidence/dict26.json`, `evidence/ir26.json`, `evidence/scan26.json`,
`evidence/scan_forced26.json`, `evidence/scan_abortctl26.json`,
`evidence/ledger25.json`+`survival25.json`+`census25.json` (carried measurements),
`scripts/gen26.py` (regenerates every number), `scripts/verify26.py` (dictionary
test), `scripts/ir26/*.c` (the shipped structure), `scripts/protocol25.py` (protocol
and scan).
