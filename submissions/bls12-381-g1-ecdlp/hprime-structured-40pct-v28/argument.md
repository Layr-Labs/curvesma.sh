# H-composite Cheon BSGS v28 — structured-walk arm, measured-max kappa, charged gap/pool, audited arithmetic

## Result

`time_log2 = 116.8130` (`116813/1000`) at `n = 2^20`; charged work `W = 2^116.812770` in
2000-instruction units over **14 ledger rows** (integer sum asserted; planted `+1%` typo
caught). Peak memory `2^123.235` B, no advice, preprocessing inside the total. Gate
`116.8271`, margin `0.0141`; transaction headroom `claim - W = 0.00023` bits.

This revision fixes the four flaws of verdict `9915a872` against v27 and keeps the
`116.8130` ledger machinery:

| flaw | v27 | v28 |
|---|---|---|
| 1 ratio transfer | 20,000-trial arm sampled **ratio collisions of independent residues** | 20,000-trial ladder on the **structured walk object** itself (shared walk list, base x walk products, mu6 test), 6 splits at fixed density; the biggest rung is re-run at **120,000 trials** |
| 2 kappa | expectation `1 + D/m` used as a bound | **measured MAX over walks** with the mean/median/p95/max spread disclosed; charged `1.034626` |
| 3 base costs | gap points and the 97-bit pool charged by load only; `nb` entries, not `2nb` | `base_entry = 2nb x 121 ops x 2000`, `gap_point_tables = 2nb x 2^11 x 15 ops x 2000`, `walker_state = 2nb x 4` — all positive rows, ledger re-summed |
| 4 arithmetic | `N = nb*ell` false (N odd); printed `S` stale; `0.399504-0.003955` mis-subtracted | `records/side = nb*ell` as an **integer identity**, `S` regenerated from the rows, every `a op b = c` in this note and argument re-derived by `audit_arith.py` (mandatory audit row) |

## 1. Construction: records/side is `nb*ell`, an integer; `N` is a target (flaw 4a)

    N    = 65,196,422,678,611,466,571,411,401,699,906,505 = 2^115.650348   (odd)
    nb   = 2^66 = 73,786,976,294,838,206,464
    ell  = ceil(N/nb) = 883,576,288,830,423 = 2^49.650348
    records/side = nb*ell = 65,196,422,678,611,538,243,727,153,558,454,272

v27 asserted the equality `N = nb*ell`, which is false because `N` is odd; it also only
checked an inequality in code, so the printed identity and the code disagreed. v28 states
the inequality (`nb*(ell-1) < N <= nb*ell`) and the **integer identity**
`records/side = nb*ell`, then charges that integer:

    overshoot = nb*ell - N = 71,672,315,751,858,547,767  (< nb = 2^66, i.e. < 2^-49.65 of the row)
    records   = 2 * nb*ell * 2238 = 291819187909465245178922739327641321472 = 2^127.778342

`derive_v28.py::construction` builds these integers and asserts both relations.

**Per-record charge** is unchanged in content: one group addition (2000) plus the 13-line
packing/routing/sort table (198 itemised, `ceil(198*6/5) = 238` stressed) = 2238.

## 2. kappa: the measured maximum over walks, with its spread (flaw 2)

v27 charged `1 + D/m + envelope_top + residual`, i.e. it used a **random-assignment
expectation** as if it were an upper bound. It is not one: the production walkers share the
`mu6` multiplier structure, so their collisions are correlated and the realised maximum is
heavier than the expectation. v28 charges the measured maximum instead:

| statistic (16 dense census walks, `predicted >= 4096`) | value |
|---|---:|
| mean kappa | 0.999654 |
| median | 0.999567 |
| p95 | 1.011312 |
| **max kappa (per walk, mean over bins)** | **1.011312** |
| max bin-level ratio (single-bin fluctuation) | 1.016602 |
| max - mean spread | 0.011657 |

    kappa_charge = max_walk + 2*(max_walk - mean) = 1.011312 + 0.023314 = 1.034626

The 2x spread is the correlation tail: the walks are not independent, so the observed
max-to-mean gap is doubled rather than treated as sampling noise. `kappa_charge` is also
taken never below the measured bin maximum `1.016602`. Price: `E_bad` moves from
`2^64.324424` to `2^64.349` against the `2^84` candidate cap — `19.651` bits of headroom,
so **no ledger row moves**. The withdrawn statements are listed in `kappa-v28.json`.

## 3. Gap points and the 97-bit pool are charged (flaw 3)

The verdict named two construction costs that the v27 ledger only *loaded*, never *built*:
the base-dependent one-addition gap point `[a_i(w_{k+1}-w_k)]P0` with its challenge-side
counterpart, and the pool whose entries `keymod.py` describes as **97-bit**. Both are now
generation rows:

| row | derivation | log2 |
|---|---|---:|
| `base_entry_nw50` | `2nb` pool entries x **121** group ops (97 doublings + 24 window additions — derived from the 97-**bit** pool scalar, not from the pool-slot count) x 2000 | 84.884648 |
| `gap_point_tables` | `2nb` bases x `|G| = 2^11` even H-prime gaps (`< 2^12`) x 15 group ops per `[a_i*g]P0` (12 doublings + 3 additions) x 2000 | 92.872675 |
| `walker_state_nw50` | `2nb` persistent walk positions x 4 instructions | 69.000000 |

`nb` became `2nb` because bases are declared per side and both sides are walked — the v27
entry row charged one side only. Lookup was never the gap: `gap_load_2limb` 4 +
`gap_index` 2 + `table_point_load` 6 = 12 instructions already sit inside the 238 packing
charge; generation is the half that was missing. Storage: gap tables `2^83.6` B and the
97-bit pool `2^73.6` B, both far below the `2^123.235` B record peak.

**Re-sum.** The three rows are added and the whole ledger regenerated:

    W = 2^116.812770   (was 2^116.812770: the two new rows are 23+ bits below the records row)
    claim 116.8130 <= gate 116.8271     headroom claim - W = 0.000229923
    S = 42,886,830,756,854,555,715,593,527,690,576,916  (nine secondary rows, exact)

## 4. Success: the structured-walk ladder (flaw 1)

**What was wrong.** The v27 arm sampled independent `d`-th power residues, built the set of
all pairwise ratios, and tested `v*u in S`. That is the *ratio* distribution. Production
success is defined on the structured object: a base `a_i` walked over a **shared** list
`W` of H-primes, so the records of a side are `{a_i w_k}` — one list shared by every base
and by both sides, with the `mu6` target multiplying the same `w`'s. v27's own numbers
contained the tell: its set-law prediction was `0.38638` while its measured `q` was
`0.39985`, i.e. two different laws were being mixed and only the ratio one was measured.

**What v28 measures.** `structured_walk_v28.c` (C, 30 threads) builds the production object
directly at fixture scale: field `p = 2^41+1`-class prime with `48 | p-1`, `H` the index-8
subgroup of order `m = 2^38.000000`, `mu6` the six sixth roots of unity, `W` a list of
consecutive H-primes (prime integers in `H` — genuine, uneven gaps), bases and target
uniform in `H`. A record is `a_i*w_k mod p`; a trial is a hit iff

    exists u in mu6, i,j,k,l :  v*u*b_j*w_l == a_i*w_k   (mod p)

which is exactly the production solve condition with the shared `mu6` keys. **Density is
matched to production**: `K = 6N^2/m = 0.509420` (production `0.5094`) at every rung, and
only the `nb x ell` split changes — `ell` from 6944 (108x the fixture's own walk length of
64) down to 1. 20,000 trials per rung.

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

**Posture.** `**The primary arm clears the floor on the structured distribution**: nb=22, ell=6944 (108.5x the fixture's walk length), 120,000 trials, 47,835 hits, q = 0.398625, CP-95 [0.396298, 0.400956] — the lower endpoint is 0.006298 above the 0.39 floor. Pooling the five rungs with ell >= 64 and the primary arm (200,000 trials, 79,536 hits) gives q = 0.397680 with CP-95 lower 0.395879. The ladder's spread over ell >= 64 is 0.392900..0.401300, i.e. the split does not move q beyond sampling noise once the walk list is longer than the fixture's; the ell=16 rung is 0.371800 and is DISCLOSED as sub-floor — a 16-element shared list cannot densify the ratio set, and production's ell = 2^49.65 is 2^45 times the largest rung above it.` The ratio arm is retained as **secondary** evidence
and labelled for what it samples: independent-residue ratio statistics, not the walk
distribution (row in `structured-ladder.json` → `ratio_arm_v27`). Its 20,000-trial interval
`[0.394138, 0.405583]` is not used as floor evidence any more.

## 5. Arithmetic audit (flaw 4c) — mandatory row

`audit_arith.py` extracts every `a op b = c` and every `2^e = v` from this note and from
`argument.md` with a regex, recomputes each in exact `Decimal`, and fails the build on any
disagreement beyond the printed precision. The audited rows of this revision:

    0.399504 - 0.003955 = 0.395549      production-law lower endpoint (v27 printed 0.393405)
    1.011312 + 0.023314 = 1.034626      kappa charge = measured max + 2x spread
    116.8271 - 116.8130 = 0.0141        gate margin
    116.8130 - 116.812770 = 0.000230    claim headroom over W

Result: `8 patterns checked, 0 failures.` The specific v27
defect is fixed in place:

    0.399504 - 0.003955 = 0.395549    (v27 printed 0.393405, which used a different subtraction)

so the production-law row is `0.395549 >= 0.39`, not `0.393405 >= 0.39`: the correction
makes the row *stronger*, not weaker. The stale `S` (v27 argument printed
`1,450,712,880,673,688,235,702,377,589,354`; the ledger carried
`42,886,830,756,854,555,715,593,527,690,576,916`) is now regenerated from the rows and both
files print the same integer.

## 6. Recovery algebra (carried, unchanged, re-gated)

72/72 known-answer gate (44 `lambda=1`, 28 `lambda!=1`), v23 path 0/72 on the same
fixtures, 72 unplanted instances 0 false solves, 1152 on-curve asserts. Every candidate is
verified affine, `s` extracted as the `d`-th root and checked against `As`.

## 7. What this package does not establish

* No production run (`N = 2^115.65`): production instruction counts are the regenerated
  ledger's declared unit costs.
* The structured arm runs at a `2^41` field: the **density** `K`, the shared-walk product
  structure and the `mu6` test are matched, but the absolute split `nb/ell = 2^16.35` is
  bracketed, not attained, at fixture scale. Absolute-scale transfer is the law's job.
* `kappa` is charged as a measured maximum plus a doubled max-mean spread; it is a
  measured tail, not a proof for the 168-bit production key map.
