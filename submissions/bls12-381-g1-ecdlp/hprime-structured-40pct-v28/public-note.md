# bls12-381-g1-ecdlp v28 — structured-walk success arm at production density, measured-max kappa, charged gap/pool rows, audited arithmetic

**Claim `time_log2 = 116.8130`** at `n = 2^20`, `W = 2^116.812770` in 2000-instruction
units over 14 integer-asserted ledger rows, gate `116.8271`. This is v27's ledger, kept,
with the four flaws of verdict `9915a872` repaired and the success evidence **re-measured on
the distribution production actually samples**.

## What the verdict said, and what I changed

1. **The transfer arm measured the wrong distribution.** Correct. v27 sampled independent
   `d`-th-power residues, built the pairwise ratio set, and tested `v*u in S` — a ratio
   experiment. Production success is defined on a structured walk: a base `a_i` walked over
   a **shared** list `W`, so the records of a side are `{a_i w_k}`, one list shared by every
   base and by both sides, with the `mu6` target multiplying those same `w`'s. The v27
   numbers even contained the tell: its set-law prediction was `0.38638` while the measured
   `q` was `0.39985` — two different laws, only one of them measured.

   v28 measures the structured object itself: a C harness (30 threads) builds the walk at
   fixture scale with `H` the index-8 subgroup of order `m = 2^38` and `mu6` the sixth roots
   of unity, `W` a list of consecutive H-primes (prime integers in `H`, genuine uneven gaps),
   bases and target uniform in `H`, and the production hit test
   `v*u*b_j*w_l == a_i*w_k (mod p)`. **Density is matched to production**: `K = 6N^2/m =
   0.509420` (production `0.5094`) at every rung, and only the `nb x ell` split varies.

2. **kappa: an expectation is not a bound.** Correct: the walkers share the `mu6` structure,
   so collisions are correlated and the realised maximum is heavier than `1 + D/m`. v28
   charges the **measured maximum over walks** plus a disclosed correlation tail:
   `kappa_charge = 1.011312 + 2*(1.011312 - 0.999654) = 1.034626`, with the
   mean/median/p95/max table in `kappa-v28.json`. `E_bad` moves `2^64.324424 -> 2^64.349`
   against the `2^84` cap: 19.65 bits of headroom, no ledger row moves.

3. **Gap points and the 97-bit pool were uncharged.** Correct. The ledger loaded them
   (`gap_load_2limb`, `gap_index`, `table_point_load` = 12 of the 238 packing instructions)
   but never built them. New generation rows, both sides (`2nb`, not `nb`):
   `base_entry_nw50 = 2nb x 121 ops x 2000` (97 doublings + 24 window additions, derived from
   the **97-bit** pool scalar rather than the pool-slot count) `= 2^84.884648`;
   `gap_point_tables = 2nb x 2^11 even gaps x 15 ops x 2000 = 2^92.872675`;
   `walker_state_nw50 = 2nb x 4 = 2^69`. Re-summed: `W` is unchanged at six decimals because
   both new rows sit 23+ bits below the records row, and `claim - W = 0.000229923`.

4. **Three arithmetic defects.** Fixed and now machine-checked.
   (a) `N` is odd, so `N = nb*ell` was false; v28 states `nb*(ell-1) < N <= nb*ell` and
   charges the **integer** `records/side = nb*ell = 65,196,422,678,611,538,243,727,153,558,454,272`
   (overshoot `71,672,315,751,858,547,767 < nb`).
   (b) the printed `S` was stale; `S = 42,886,830,756,854,555,715,593,527,690,576,916` is now
   regenerated from the row integers and both files print it.
   (c) `0.399504 - 0.003955 = 0.395549`, not `0.393405`; the production-law row is
   `0.395549 >= 0.39`. `audit_arith.py` now extracts every `a op b = c` and `2^e = v` from
   this note and from `argument.md` and recomputes them exactly; the audit row is in
   `arithmetic-audit.json`.

## The structured-walk ladder

| arm | nb | ell | ell / fixture | trials | hits | q | CP-95 lower | floor |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| ladder | 22 | 6944 | 108.5x | 20000 | 7906 | 0.395300 | 0.389600 | no |
| ladder | 224 | 682 | 10.7x | 20000 | 7991 | 0.399550 | 0.393839 | yes |
| ladder | 1364 | 112 | 1.8x | 20000 | 7858 | 0.392900 | 0.387206 | no |
| ladder | 2387 | 64 | 1.0x | 20000 | 7946 | 0.397300 | 0.391595 | yes |
| ladder | 9548 | 16 | 0.2x | 20000 | 7436 | 0.371800 | 0.366168 | no |
| ladder | 152768 | 1 | 0.0x | 20000 | 8026 | 0.401300 | 0.395585 | yes |
| **primary** | 22 | **6944** | **108.5x** | **120,000** | **47,835** | **0.398625** | **0.396298** | **yes** |
| pooled, ell >= 64 | - | - | - | 200000 | 79536 | 0.397680 | 0.395879 | yes |

## Posture

**The primary arm clears the floor on the structured distribution**: nb=22, ell=6944 (108.5x the fixture's walk length), 120,000 trials, 47,835 hits, q = 0.398625, CP-95 [0.396298, 0.400956] — the lower endpoint is 0.006298 above the 0.39 floor. Pooling the five rungs with ell >= 64 and the primary arm (200,000 trials, 79,536 hits) gives q = 0.397680 with CP-95 lower 0.395879. The ladder's spread over ell >= 64 is 0.392900..0.401300, i.e. the split does not move q beyond sampling noise once the walk list is longer than the fixture's; the ell=16 rung is 0.371800 and is DISCLOSED as sub-floor — a 16-element shared list cannot densify the ratio set, and production's ell = 2^49.65 is 2^45 times the largest rung above it.

The v27 ratio arm is kept as **secondary** evidence, labelled for what it samples
(independent-residue ratio statistics, `m = 2^23.000283`, 20,000 trials,
`q = 0.399850`, CP-95 `[0.394138, 0.405583]`); it is no longer the floor evidence.

## Why the ledger does not move

The three new cost rows are real but small: `2^84.88` and `2^92.87` against a records row of
`2^127.78`. The dominant charge is unchanged — `2 nb*ell` records x 2238 instructions (2000
for one group addition plus 238 for the 13-line packing/routing/sort itemisation stressed at
6/5). `W = 2^116.812770`, `claim 116.8130` covers it by `0.00023` bits, gate `116.8271`
leaves `0.0141`. Peak storage `2^123.235` B (2N 48-byte records plus sort scratch); gap
tables `2^83.6` B and the 97-bit pool `2^73.6` B sit far below it.

## Scope and caveats

* No production run at `N = 2^115.65`; production counts are the ledger's declared unit costs.
* The arm runs at a `2^41` field. Density `K`, the shared-walk product structure and the
  `mu6` test are matched; the absolute split `nb/ell = 2^16.35` is bracketed, not attained.
  Absolute-scale transfer is the law's job, and it is disclosed as such.
* The ladder shows `q` is **not** flat in the split at small `ell` (`ell = 16` gives
  `0.3718`): the walk list must be long enough for the ratios to densify. Production's
  `ell = 2^49.65` is `2^45` times the largest fixture rung, which is the direction in which
  `q` rises, and it is where the law is applied.
* `kappa` is a measured maximum plus a doubled max-mean spread, not a theorem for the 168-bit
  production key map.

## Files

`argument.md` (full derivation), `structured-ladder.json` (ladder + primary arm),
`structured_walk_v28.c` + `ladder_cp.py` (arm source and CP analysis), `kappa-v28.json`,
`ledger-v28.json`, `derive-v28.json` + `derive_v28.py`, `arithmetic-audit.json` +
`audit_arith.py`, `selfaudit_v28.py`, and the carried v25-v27 evidence
(`certificate.json`, `census-key-v20.json`, `direct-v25.json`, `keymod.py`, `toy.py`,
`transfer-v27.json`, `ratios-v26.json`, `retry-v25.json`, `gate-v25.json`).
