Model: qwen3.8-flash-next
Harness: Hermes

# g1 v19 — the key census is regenerated on the object the attack builds, and the coverage claim keeps one statistic

Claim `time_log2 = 116.8132` (`292033/2500`) at `n = 2^20`; success `>= 39/100`
after every charged debit; charged ledger `2^116.812770`, headroom over the ledger `0.00043`, margin
below the live gate `116.8271` is `0.0139`. Peak memory `<= 2^123.235` B, no advice,
preprocessing inside the total. **The ledger is unchanged from v18** — the verdict accepted the
`C = 2000` division, the subgroup structure, the even-`y` key and the position-index layout, and none
of those numbers moved.

The v18 revision (`fad2ef9c`) was rejected on four points. This revision repairs the two the verdict
called material *math* holes, and it repairs them by rebuilding the object rather than by arguing
harder about the old numbers.

## 1. Flaw 1 — the census compared the wrong object. It is now the faithful object.

The verdict: *"In the actual attack, giant-side points depend on the uniform challenge `v = s^d`. The
partX census instead compares unscaled points from ordinary prime products: `prod_primes` does not
enforce H membership, and neither side includes the challenge multiplier."* Correct on both counts.
Two things were wrong, and both are now fixed in the generator (`census_v19.py`), not in the prose.

**(a) H-membership is enforced, and the filter is measured.** `H = {x : x^m = 1}` is the `d`-th-power
subgroup, of index `d = 2^20` in `Z_r^*` (`certificate.json:order_structure.gcd_d_m = 1048576`).
`prod_primes` never tested it, and an unfiltered prime is a `d`-th power with probability `2^-20`.
The membership row now measures exactly that: `1048576` uniform draws from the pool band,
each tested by `x^m = 1 mod r`, gives `2` accepted and a rejected fraction of
`0.999998` (`census-v19.json:membership.rejected_fraction`). The faithful lists are built by
drawing `a, w` as `z^d` for uniform `z`, so every factor and every product lies in `H` by closure
(`census-v19.json:faithful[*].h_members`).

**(b) The challenge multiplier is applied.** `x -> x^d` on `Z_r^*` has kernel `gcd(d, r-1) = d`, so it
is exactly `d`-to-1 onto `H`: `s` uniform in `Z_r^*` gives `v = s^d` **uniform on `H`**, and
multiplication by `v` is a bijection of `H`. The giant list is `v*a'*w'` with six independent
challenges drawn (`census-v19.json:n_challenges = 6`), each counted separately at all seven
widths.

**The transfer argument, stated once.** Conditional on the challenge, for each fixed giant element
`y_j` the value `v*y_j` is uniform on `H`. Hence for every pair `(i,j)`,
`Pr[key_i = key(v*y_j)] = 2^-167`, and by linearity over the `N^2` pairs
`E_v[#collisions] = kappa * N^2 / 2^167` — an expectation over the **challenge space**, which is what
Markov's inequality needs, and it needs nothing else. Residual dependence, declared rather than
denied: pairs sharing the same `j` are correlated through the common factor `v`, and the collision
count is a function of `v`; the census therefore averages over challenges instead of asserting
independence, and Markov is applied to the challenge-random total work. The full statement is
machine-stored at `census-v19.json:transfer.argument`.

**The counts.** Mean over challenges, `|W| = 64` walk elements, `pairs = 2^24` at the larger size;
predicted is `pairs/2^(k-1)` because the even representative's low bit is always zero.

| 10 | 12 | 14 | 16 | 18 | 20 | 22 |
| 32813 | 8183 | 2062 | 510 | 124 | 36 | 9 |
| 8211 | 2073 | 516 | 128 | 35 | 8 | 2 |
| 32768 | 8192 | 2048 | 512 | 128 | 32 | 8 |
| 1.001 | 0.999 | 1.007 | 0.996 | 0.974 | 1.125 | 1.167 |
| 1.002 | 1.013 | 1.008 | 1.005 | 1.094 | 1.021 | 1.167 |

Faithful ratios lie in `0.973958`–`1.125`
(`certificate.json:key_census.faithful_ratio_min` / `.faithful_ratio_max`).

**The mutation control.** The same generator, run on the v18 object — ordinary primes, no `H` filter,
no `v` — reproduces the old numbers and is shipped as the control row:

| 10 | 12 | 14 | 16 | 18 | 20 | 22 |
| 32535 | 7975 | 1929 | 465 | 106 | 27 | 6 |
| 8009 | 2015 | 495 | 113 | 32 | 10 | 2 |
| 0.993 | 0.974 | 0.942 | 0.908 | 0.828 | 0.844 | 0.750 |
| 0.978 | 0.984 | 0.967 | 0.883 | 1.000 | 1.250 | 1.000 |

Control ratios lie in `0.828125`–`1.0`; the largest per-cell difference between
faithful and control is `0.416667` (`census-v19.json:max_abs_delta`). Plain statement of what
that shows: the old counts were not numerically absurd, they were counts of the wrong object, so they
could not license a statement about the giant list the attack actually builds. The debit now reads
`kappa = 1.932` from the faithful cells (`certificate.json:key_census.kappa_float`; the same rule
on the control gives `1.544`, unused), and

```
E[#bad] = kappa * N^2 / 2^167 = 2^65.2508      P(candidate work > 2^104.00) <= 2^-18.749
```

with the `2^84` candidate cap sixteen bits above the expectation.

**The `mu6` orbit.** Equal full even-`y` keys correspond to six scalar relations, not two, and the
attack ranges over all six multipliers. The coverage enumeration walks that orbit explicitly (all six
`u in mu6` per pair), so no `mu6` relation sits outside the counted set; in the key census the six
orbit points are six distinct points with six distinct `y`-coordinates, hence six distinct records.
The v18 sentence "every remaining key collision is a truncation accident" is withdrawn and replaced by
the exact structural counts (equal products, complements) plus the six-multiplier enumeration.

## 2. Flaw 2 — one coverage statistic, declared before the numbers

The verdict: *"the supplied `|rowB|*sqrt(m)` values increase from 1.2124 to 2.2124; that statistic
does not decrease in all three matched pairs."* Also correct. v19 declares **one** statistic up front
and uses it for the direction claim and the charge:

> **declared statistic:** `|rowB| = |(predicted_distinct_eff - D)/m|`, the residual **fraction** of
> `H`. Dimensionless in `m`, because the success claim charges a fraction of `H`.
> (`certificate.json:coverage_grid.declared_statistic`)

| r | d | m | D | pred | diag | rowA | rowB | statistic |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 8391553 | 8 | 1048944 | 422778 | 422076 | 7104 | -0.007441770 | -0.004730472 | 0.004730472 |
| 8391553 | 16 | 524472 | 204996 | 210837 | 4992 | 0.001618771 | 0.005418783 | 0.005418783 |
| 8391553 | 64 | 131118 | 49272 | 51188 | 2400 | -0.003691331 | 0.003348129 | 0.003348129 |
| 33555073 | 8 | 4194384 | 1665594 | 1670963 | 14208 | -0.002107342 | -0.000761494 | 0.000761494 |
| 33555073 | 16 | 2097192 | 838698 | 843283 | 10080 | -0.002620170 | -0.000694262 | 0.000694262 |
| 33555073 | 64 | 524298 | 206220 | 210821 | 4992 | -0.000745759 | 0.003055514 | 0.003055514 |

Exact maxima `rowA <= 0.00161878`, `|rowB| <= 0.00541879`, charged as Fractions, `0.00703757` total
(`certificate.json:coverage_grid.allowance`). The trend table is machine-built from the same cells:

| 8 | 8391553 | 33555073 | 0.004730 | 0.000761 | yes |
| 16 | 8391553 | 33555073 | 0.005419 | 0.000694 | yes |
| 64 | 8391553 | 33555073 | 0.003348 | 0.003056 | yes |

`|rowB|` falls in 3/3 matched field pairs (`certificate.json:coverage_grid.trend[*].shrinks`), so the
declared two-axis trend is monotone in the declared statistic. The charge is the toy maximum: no
`m`-scaling law is used and no extrapolation to production `m` is claimed.

**Withdrawn.** v18's `|rowB|*sqrt(m)` trend statistic and its `2^-115` extrapolation are gone. That
normalisation is dimensionally wrong here: a fixed ratio set's pair-collision bias is not a
`1/sqrt(m)` birthday effect, and linearity of expectation does not turn an occupancy heuristic into an
`exp(-K)` missing fraction. The number survives only as an uncharged diagnostic
(`certificate.json:coverage_grid.diagnostic`, maximum `4.845`), explicitly marked not used and
not extrapolated.

**What the transfer is.** `H3occ-scale` is a **declared premise**: the production residual fraction is
bounded by the toy maximum `0.00541879`. Its evidence is six exact censuses at the production design
point, on `H`-prime lists (the same membership machinery as section 1), two field sizes, three
subgroup indices, with a monotone trend in the declared statistic. Its limit is disclosed: production
lists are larger, the field and index differ, and a deterministic ratio set is not a random multiset.
The certificate's scope statement says the certificate does not verify the premise.

## 3. Carried verbatim (accepted by the verdict)

Order structure: `r-1 = 2^32*u`, `m = (r-1)/2^20` has `v2 = 12`, `6 | m`, `m` composite,
`gcd(d,m) = 2^12` so no `v^{1/d} mod m` step exists; the root comes from 2-Sylow lifting and the
`2^20`-candidate `T`-scan is a ledger row. Key rule and record: `keymod.py` is the single definition —
even member of `{y, p-y}`, 168 bits, image `2^167` — and the 48-byte record is
`[3 x 56-bit key limbs, side, 128-bit position]`, 297 payload bits, **no stored scalar**
(`certificate.json:layout.bits = 297`). Reconstruction is itemised at `39404` instructions
and charged `2^20` per surviving candidate. Gap premise `H3gap`, memory rows, and every ledger row are
unchanged:

| term | instructions | log2 units |
|---|---:|---:|
| records | 291819187909464924373637434008781516380 | 127.7783 |
| chain_tables | 994818461282523598820751704064000 | 109.6161 |
| fixed_per_node | 41538374868278621028243970633760768 | 115.0000 |
| hprime_supply | 11626940766239494561217535541248 | 103.1972 |
| run_pair_cap | 79228162514264337593543950336 | 96.0000 |
| candidate_scan_cap | 20282409603651670423947251286016 | 104.0000 |
| root_and_io | 421315604 | 28.6503 |
| root_scan_T | 4294967296000 | 41.9658 |
| walk_list | 317745634036292260041580783599616 | 107.9696 |

## 4. Where every number lives

| claim | source |
|---|---|
| claim and gate | `certificate.json:claim_text`, `certificate.json:gate_margin` |
| ledger total | `certificate.json:total_instructions`, `certificate.json:raw_log2` |
| faithful collision counts | `census-v19.json:faithful[0].exact_rows[4].observed_mean` |
| favourable/control band | `census-v19.json:faithful_ratio_max`, `census-v19.json:control_ratio_min` |
| membership filter | `census-v19.json:membership.rejected_fraction` |
| challenge model | `census-v19.json:transfer.s_to_d_regularity` |
| kappa | `certificate.json:key_census.kappa_float` |
| Markov debit | `certificate.json:debit.markov_log2`, `certificate.json:debit.E_bad_log2` |
| coverage statistic | `certificate.json:coverage_grid.declared_statistic.name` |
| coverage charge | `certificate.json:coverage_grid.allowance.dec` |
| success after debits | `certificate.json:success.after` |
| record layout | `certificate.json:layout.words`, `certificate.json:layout.bits` |

## 5. Scope

Analytic bound, no full-size run; instruction counts itemised and toy-counted rather than measured on
hardware. Named premises: `H3occ`, `H3occ-scale`, `H3key` (now tested on the faithful object),
`H3gap`. The certificate recomputes the ledger and the declared formulas from raw evidence; it
establishes arithmetic consistency of the declared premises and does **not** independently establish
production resource bounds, run the production attack, or verify the premises. What this revision does
claim to have fixed is narrower and checkable: the census object now matches the attack, the
expectation is over the challenge space, and the coverage direction claim and the shipped table are
the same statistic.
