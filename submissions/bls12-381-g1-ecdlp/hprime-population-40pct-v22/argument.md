# H-composite Cheon BSGS v22 — population census, deleted coverage debit, counted self-check

## Result

`time_log2 = 116.8130` (`116813/1000`) at `n = 2^20`; success `>= 39/100` collision-only
(`0.399502975`); charged work `2^116.812770` in 2000-instruction units; peak memory
`2^123.235` B; no advice; preprocessing inside the total. Every number is emitted by
`population_v22.py`, `gapcheck.py` and `build_certificate_v22.py` into
`population-v22.json`, `retry-v22.json` and `certificate.json`.

## 1. Flaw 1 — the population, not a realized list

The defect named in the v21 verdict was a histogram built from one challenge-scaled giant
list (sums to `|Y|`, not `m`) and then averaged over three challenges. `population_v22.py`
replaces it with an exhaustive enumeration: `H = <g^d>` has order `m = (r−1)/d`, every one
of its `m` elements is visited, `[h]P0` is computed for each, and

    N_w(k) = |{h in H : key_w([h]P0) = k}| ,   sum_k N_w(k) == m

is asserted per width per fixture (`population-v22.json:fixtures[*].rows[*].sum_Nk`; the
assertion is in `population_rows`, not in the note). Nine fixtures, `m = 2^11, 2^13, 2^15`,
`nb` 32..256, three independent challenges each, 36 scored cells.

For the **fixed realized** baby list, `E_w = (|Y|/m)·sum_k B_w(k)N_w(k)` is compared with
the realized `C_w(v)` of independent challenges: `E/C` spans **0.9026..1.2687**
(mean 1.0361), i.e. counting noise around the exact mean with no systematic offset.

The weights row answers the verdict's second point (`sparsity ≠ uniform population weights
on occupied baby keys`) by measurement rather than algebra:

| density `m/D` | 2^3 | 2^5 | 2^7 | 2^9 | 2^11 |
|---|---|---|---|---|---|
| `weight_ratio` | 1.111 | 1.0226 | 1.0043 | 1.0007 | **1.000000** |

(`certificate.json:population.weight_trend`.) The shortcut `sum B(k)N(k) = |X|m/D` is
therefore *not* assumed at the sparse end, where it is measurably wrong; it is the
`m/D → ∞` limit of a measured trend, and production sits 57 bits beyond the last measured
point. The charge `kappa = 1.337` (one-sided UCB of the exact estimator, cells with
baseline `>= 256`) is carried from `census-key-v20.json`, giving `E[#bad] = 2^64.721` and
Markov `2^-19.279` against the `2^104` candidate cap.

## 2. Flaw 2 — the coverage debit is deleted, and the matched cell is infeasible

No coverage residual is charged: `certificate.json:coverage.debit = 0`,
`success.debits.coverage = 0`. The reason is arithmetic, not convenience:

    N = nb·nw,  K = 6N^2/m   =>   m = (6/K)·(nb/nw)^2·nw^4

so `nb/nw = 2^15.65` at `K = 0.51` forces `m = 2^34.86` at `nw = 1` — a 34-bit cover set
before any walk exists — and `m = 2^234.86` at production `nw = 2^50`, reproducing the
production `m`. A cover census must materialise the covered set; the largest population
enumerated here is `m = 2^15` (bitmap frontier about `2^26`). Shortfall: **8.9 bits** at
the bitmap frontier, **19.86 bits** at the enumerated scale.

The frontier ladder (28 cells, two field scales, two `d`, `K` matched, base-heavy
`nw ∈ {1,2}`, `max|rowB| = 0.00270272`) is therefore an observation, not a transfer. The
success claim is collision-only and the reserve to the floor is `0.009504421` = 3.52× that
residual — disclosed as unpriced risk, never as a bound.

## 3. Flaw 3 — retries at `1/(1−q)`, gap and route separate

192 independent seed draws (2 scales × 96; bands of 401 consecutive H-primes):
`q_gap = 0/192`, `q_route = 0/192`, Clopper-Pearson 95% upper `0.0154817` each. Extra work
is `(1/(1−q_hi) − 1) × (one discarded draw)`. A discarded draw costs the H-prime harvest +
drawn walk list + self-check = `2^108.022`, so `gap_retry_worst = 2^102.031`; an aborted
routing run discards its `2^96` pair-cap work, so `route_retry_worst = 2^90.009`. The
retry precedes the record build, so the build is not multiplied by the retry factor.

## 4. Flaw 4 — the self-check is implemented and its price is counted

`gapcheck.scan_max_gap` maintains `(prev index, running max)` with one compare and one
conditional store per step, aborting early on violation; mutation control: clean list
accepted, injected violation rejected, widened limit accepted. The per-step price is read
out of the disassembly of `gapscan.c` (`clang -O2 -arch arm64`): 9 instructions in the loop
body (`ldr, sub, cmp, csel, cmp, b.hi, mov, subs, b.ne`), charged 11 at the 6/5 stress
convention, times the walk-record count already used by `walk_list`:
`gap_self_check = 2^94.388`, a real `certificate.json:rows` entry.

## 5. Ledger and claim

12 rows: `records` `2^127.778`, `chain_tables` `2^109.616`, `fixed_per_node` `2^115.000`,
`hprime_supply` `2^103.197`, `run_pair_cap` `2^96`, `candidate_scan_cap` `2^104`,
`walk_list` `2^107.970`, `gap_self_check` `2^94.388`, `gap_retry_worst` `2^102.031`,
`route_retry_worst` `2^90.009`, `root_scan_T` `2^41.966`, `root_and_io` `2^28.650`.
Total `W = 2^116.812770`; claim `116.8130` (headroom `0.00023`); gate `116.8271`
(margin `0.0141`). Success: `0.399504421 − 0.000001446 = 0.399502975 >= 0.39`.

## 6. Carried unchanged

Order structure (`r−1 = 2^32u`, `v2(m)=12`, `6|m`, `gcd(d,m)=2^12`, 2-Sylow root,
`2^20`-scan); `keymod.py` single key rule (even member, 168 bits, image `2^167`);
48-byte position-index record; `C = 2000`; per-record `2238 = 2000 + 238` from 12 itemised
rows summing to 198; routing loads at the measured max/mean `2.25`, prefix `1.120`,
`2^84.47` against the `2^95` cap.

## 7. Scope

Analytic bound, no full-size run. Named premises: `H3key` (measured envelope `1.337`),
`H3gexp` (uniform-weights limit extrapolated from the measured `m/D` trend),
`H3gap` (the `2^7.7` ratio bound with 192 accepts and a CP-95 retry price). The certificate
recomputes the ledger and the population, retry and self-check rows from raw evidence; it
does not independently establish production resource bounds.
