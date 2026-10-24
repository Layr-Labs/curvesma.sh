# H-composite Cheon BSGS v29 — derived gap scale, exact-product 97-bit arm, per-walk kappa, extended label audit

`time_log2 = 116.8130` (`116813/1000`) at `n = 2^20`; charged work `W = 2^116.812794` in the
2000-instruction unit, every row a regenerated Python integer. Gate `116.8271`; the claim sits
0.0143 below the gate and 0.000206 above `W`. Peak storage `2^123.235` B, no advice,
preprocessing inside the total.

## 1. Geometry (unchanged, carried)

    r-1 = 2^32 u,  d*m = r-1,  d = 2^20,  m = 2^234.857   (certificate order_structure)
    N    = 65,196,422,678,611,466,571,411,401,699,906,505 = 2^115.650348   (odd)
    nb   = 2^66 bases/side,  ell = ceil(N/nb) = 883,576,288,830,423 = 2^49.650
    records/side = nb*ell (exact integer), 2*rps*2238 = 291819187909465245178922739327641321472 = 2^127.778342

Per-record charge 2238 = 2000 (one group addition) + 238 (packing, itemised 198 x 6/5 stress rule).

## 2. FLAW 1 — the gap scale is DERIVED now, and priced at 2^25.72 / support 2^30

v28 charged gap points at `|G| = 2^11` even values `< 2^12`, and `2^12` was read off
`keymod.py`'s `HOME_BITS` — the number of **recovery steps** between checkpoints, a different
object from the gap **magnitude**. The production walk is a list of consecutive H-primes
`w_0 < w_1 < ...` and a step is `g_k = w_{k+1} - w_k`.

Derivation (row, not assertion). An integer near `q` lies in the index-`d` subgroup `H` with
density `1/d` and is prime with density `1/ln q`, so consecutive H-primes near `q` are spaced

    mu(q) = d * ln q
          = 2^20 * ln(2^76)                      walk band q_w ~ 2^76 (97-bit pool scalar
          = 2^20 * 52.68 = 5.5238e7 = 2^25.719   x 77-bit walk element = the 174-bit key scalar)

The shared list holds `ell = 2^49.65` gaps; the maximum of that many geometric gaps is
`mu*(ln ell + gamma) = 2^25.719 * 35.0`, i.e. `2^30.848`. The even-gap support therefore runs to
`< 2^31`, i.e. `|G| = 2^30` entries per base — **19 bits above** the v28 charge of `2^11`, and
the mean gap itself is `2^25.72`, 13.7 bits above the `2^12` label the verdict named.

Construction (supplied, not assumed): `T_i[g] = [g*a_i]P0` for every even `g < 2^31` is built by
the **even-gap addition chain** `T_i[2] = 2*(a_i P0)`, `T_i[2j] = T_i[2(j-1)] + T_i[2]` — one
group addition per entry, charged at 15 ops/entry (12 doublings + 3 additions) as a stress rule.
Per base that is `2^30` additions against `ell = 2^49.65` walk steps, i.e. `2^-19.65` of the walk,
which is why a 19-bit larger table still moves `W` by `2.4e-5` bits.

Fallback (priced, so the row is a bound rather than an expectation): a gap above the support is
handled per base as a 32-doubling small-scalar multiple of `a_i P0` (62 group ops). Under the
same density law `P(gap > 2^31) = exp(-2^31/mu) = 2^-56.09`, so the list's expected number of
such gaps is `0.0115`; **one** worst-case occurrence over all `2nb` bases is charged
(`gap_tail_fallback = 2nb*62*2000 = 2^83.943`). Any further occurrence costs `2^83.9`, absorbed
by the `2^43` bits of headroom.

Ledger effect: `gap_point_tables` 2^92.873 -> 2^111.873, `W` 2^116.812770 -> 2^116.812794
(+2.35e-5 bits). The gap row is still `2^15.9` below the records row: **the re-charge does not
threaten the claim**, and the row now exists at the right scale with a construction and a tail.

## 3. FLAW 2 — what the 97-bit quantity IS, and a same-structure arm that measures it

`keymod.py` stores a 97-bit **integer exponent** per pool slot (`pool_load_2words`: 2 LDP from the
stored 97-bit pool entry) and does **not** store the scalar. The record scalar is
`a_i * w_k`: 97-bit base x 77-bit walk element = a 174-bit product, and `174 < log2(m) = 234.857`,
so the product **never wraps mod m**. Production ratios are ratios of exact integers, not residues
of a uniform draw. That is the mismatch the verdict named: the v28 arm drew `a_i` uniform in `H`,
whereas production draws `a_i` uniform in `[1, 2^97)`, a `2^-137.9` sliver of `H`.

What the hit test actually sees is the ratio set `{a_i w_k / (u b_j w_l)}`: `N^2` base-ratios
times walk-ratios times 6 (`mu6`) = `2^233.89` candidate ratios against `|H| = 2^234.857`, so the
pair law `1 - exp(-K)` with `K = 6N^2/m` is a statement about the ratio map, not about the base
support — and `K` is held at the production value in every arm below.

Arm (same walk object, same density, production pool law): bases drawn from `[1, 2^18)` instead of
uniform in `H`, with `span + walk_band = 18 + 21 < 41`, so `a_i * w_k` is an **exact product (no
mod-p wrap)** — the production regime at reduced absolute scale, and the only place a 64-bit
harness can hold both the exactness and the pool/group ratio.

    nb=22, ell=6944, K=0.509420246, trials=24000, span 18 (exact products, no wrap)
    hits=9504, q=0.396000, CP-95 [0.390796, 0.401223], pair law 1-exp(-K)=0.39916, floor 0.39

Matched control, identical geometry and seeds, only the base law reverted to full support:
24,000 trials, 9,524 hits, q = 0.396833, CP-95 lower 0.391627. The pool law moves q by
-0.00083 (0.2 sigma) — the base-exponent restriction is invisible at this density.

The pool-scalar arm clears the floor on its own CP interval (lower 0.390796 >= 0.39) and agrees
with the full-support control to 0.2 sigma, which is what the ratio-map argument predicts: the
hit test sees `{a_i w_k/(u b_j w_l)}`, and restricting a 97-bit exponent support does not shrink
that set at fixed `N`. Residual band and what is NOT claimed are stated in §7.

## 4. FLAW 3 — kappa is measured per WALK, not on census aggregates

v28's "max over 16 walks" was 12 big-fixture + 4 small-fixture **key-width aggregates** from
`kappa-v26.json` — mean-of-challenges rows, not walks. v29 measures collisions **inside each
walk**: one trial is one full record set at the primary rung; for each walk

    kappa_walk = (# colliding record pairs in that walk) / (uniform expectation C(N,2)/m)

| walks | expectation/walk | mean | p50 | p95 | p99 | MAX |
|---|---|---|---|---|---|---|
| 24,000 | 0.042451 | 0.939309 | 0.000000 | 0.000000 | 23.556344 | 70.669032 |

The per-walk count is integer-valued with expectation 0.042451, so one extra duplicate moves
`kappa_walk` by 23.56; the MAX is therefore a conservative charge by construction (p99 = 23.56 is
exactly one duplicate, the max is three). K is held at the production value and the duplicate
expectation depends only on K, so this per-walk tail is scale-free between fixture and production.
Charging `kappa_charge = 70.669` prices the debit at `E_bad = 2^70.443675` against
`cap_bad = 2^84.0` — headroom 13.56 bits; the working row does not move (kappa prices the
failure debit only). This retires the v28 aggregate-based 1.011312 charge and the `+2*spread`
heuristic with a measured per-walk maximum.

## 5. FLAW 4 — the label audit is extended, and it ships

`audit_labels.py` (in this tree, intake may re-run it) audits, over `argument.md`,
`public-note.md` and `submission.json`:

* every `2^y` token against its paired numeric magnitude **within one sentence**, both directions,
  recomputing `log2` of the number and asserting it against the printed exponent;
* every percentage claim against the ratio it names (hits/trials, secondary/dominant, ...);
* every `x`-multiple claim against the recomputed multiple;
* **counts** ("four rungs with ell >= 64"): `len([r for r in ladder if r["ell"] >= 64])`;
* direction words (`above`/`below`) against the operands' actual order;
* a **planted-flip control** per class (flipped exponent, `+1` count, `x1.05` percentage,
  `+0.5` exponent) which must FAIL — the audit exits non-zero if any control passes silently.

Result: `audit-labels-v29.json`, regenerated by `audit_labels.py` on this tree: every `2^y`/
magnitude pair, percentage, multiple, count and direction word in the three documents passes, and
all four planted flips are caught (the script exits non-zero if any control survives). The two v28 numerical errors are corrected at
the source and re-derived: `ell/6944 = 883576288830423/6944 = 2^36.89` (not `2^45`), and `ell >= 64`
selects **four** ladder rungs plus the primary arm (not five).

## 6. Success evidence (carried crown, plus the new arm)

PRIMARY (unchanged, the row the objection was resolved on): the structured ladder at the
production object — shared walk list, `records = a_i w_k`, shared `mu6` test, density
`K = 6N^2/m` held at every rung — **nb=22, ell=6944 (108.5x the fixture walk), 120,000 trials,
47,835 hits, q = 0.398625, CP-95 lower 0.396298 >= 0.39 floor**. Pooled over the four rungs with
`ell >= 64` plus the primary arm: 200,000 trials, 79,536 hits, q = 0.397680, CP-95 lower 0.395879.
The `ell=16` rung is disclosed as sub-floor at 0.371800. Fixture arm (v25, carried): 131,000
trials, q = 0.396844. Pool-scalar arm: §3.

## 7. Scope and disclosures

* No production run: instruction counts are regenerated from exact integer rows; the success
  evidence is measured at fixture scale and at the production *object* and *density*, not at
  `N = 2^115.65`.
* The 97-bit arm matches the pool law, the exact-product regime, the walk object and the density;
  it does **not** attain the absolute split `nb/ell = 2^16.35` (bracketed by the ladder), and its
  CP interval bounds sampling error, not the modelling residual. The residual is disclosed: if
  the ratio map de-clumped at production scale, the floor would move with it; nothing here
  extrapolates a fixture CP interval over `2^138`.
* kappa prices the failure debit at `2^E_BAD` (cap `2^84.0`); it does not enter the work row.
* Storage: records `2^122.235` B, sort scratch `2^123.235` B peak; gap tables `2^102.6` B,
  97-bit pool `2^70.6` B — all inside the peak.
* The pool arm's first run drew the bases as bare integers, which puts records outside `H` and
  collapsed q to 0.0627 for that reason alone; that run is recorded in `derive-v29.json` as a
  disclosed confound, not used, and the arm was re-run with the production structure (small
  discrete log, membership automatic). The fixture cannot attain the absolute
  `2^17`-of-`2^234.857` exponent sliver: 18 of 38 bits is the widest exact-product span a 41-bit
  fixture admits, and that gap is disclosed rather than papered over.

## 8. Ledger (every row a regenerated integer, `total == sum(rows)`)

| row | value | log2 |
|---|---|---|
| records | 291819187909465245178922739327641321472 | 127.778342 |
| chain_tables | 994818461282523598820751704064000 | 109.616132 |
| fixed_per_node | 41538374868278621028243970633760768 | 115.000000 |
| hprime_supply | 11626940766239494561217535541248 | 103.197242 |
| run_pair_cap | 79228162514264337593543950336 | 96.000000 |
| candidate_scan_cap | 20282409603651670423947251286016 | 104.000000 |
| root_and_io | 421315604 | 28.650326 |
| root_scan_T | 4294967296000 | 41.965784 |
| walk_list | 317745634036292242930980832149504 | 107.969571 |
| gap_self_check | 25924182448222977157188550656 | 94.388285 |
| gap_retry_worst | 3877290542265193293613612662784 | 101.612892 |
| base_entry_nw50 | 35712896526701691928576000 | 84.884648 |
| gap_point_tables | 4753689750855860255612637020160000 | 111.872675 |
| gap_tail_fallback | 18594318026299228028928000 | 83.943147 |
| walker_state_nw50 | 590295810358705651712 | 69.000000 |
| **W** | **2918620811728985960903916068037924725/20** | **116.812794** |
| claim | 116813/1000 | 116.813000 |
| gate | 1168271/10000 | 116.827100 |

`S = 42886830756854555715593527690576916` (exact integer sum of the nine secondary rows).
Control: `+1%` on the records row moves `W` to `2^116.827123`, above the gate — the ledger's
planted typo is detected, so the rows are load-bearing and the sum is not decoration.
