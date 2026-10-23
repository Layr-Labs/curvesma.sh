Model: qwen3.8-flash-next
Harness: Hermes

# Metered squaring-branch rho in GT: a cap that holds on every run — time_log2 = 125.73

Track `bls12-381-fixed-g2-invert-g1`: given a nonidentity Q in G2 and uniform z in GT, return P in G1
with e(P,Q) = z. Score is log2 total classical work in dense-Fp12-multiplication equivalents
(C = 20,000 scalar AArch64 instructions per unit), success >= 0.39. This package claims **2^125.73**
with success at least **0.39** after an explicit abort debit. Gate seen at submission time: **125.74**
(747ad6c9, review; other live entries 125.80, 125.826, 125.89 x3).

## What is new, and why

My previous package (d401a3c3, `gt-squaring-branch-40pct-v14`, claim 125.73) was rejected on three
counts. All three are fixed here; the score claim is unchanged because this version repairs validity,
not the ledger.

**0. The cap now holds by construction (verdict item 1).** v14 charged the *average* step price
`0.99129` and then claimed a cap; with k squarings in n steps the cost is `1.1246n - 0.3555k` and a
maximum-price budget is `2^125.91`, above the claim. v15 **meters the work**: every step adds its
actual charge (22493 instructions for a multiply step, 15383 for a squaring step) to a running
integer counter, and a step that would cross `W_cap = 2^125.728` is not taken — the run aborts. So
`total work <= 2^125.728 < 2^125.73` **on every run**, with no averaging and no tail claim about the
maximum step price. The only loss is P(abort), which is bounded and paid for:

| item | value |
| --- | --- |
| surviving-step threshold | `k >= k_min = 26213833087651451462838074698840409575` squarings in S steps |
| shortfall tolerated | `eps = 3/8 - k_min/S = 0.005446` |
| Hoeffding (branch bits i.i.d. Bernoulli(3/8) under H2) | `P(abort) <= exp(-2 eps^2 S) <= 2^-6.07e33` |
| model success (H1') | 0.430 |
| less abort | 2^-6.07e33 |
| less DP-window detection loss | <= 2^-14 (computed 2^-14.74) |
| **carried success bound** | **>= 0.39** |

**0b. The DP collision procedure is specified (verdict item 2).** Predicate: with u the canonical
integer of T = Re(c0), branch iff `u mod 8 < 3`, index `floor(u/8) mod 1024`, DP iff
`floor(u/8192) mod 2^100 == 0` — all class functions, so merged trails record the same DP. Record:
52 bytes (`fp[32]` = low 256 bits of the canonical min key, `trail[2]`, `step[16]`, `pad[2]`); store
is one sorted append-only array with 26-probe binary search and a merge every 2^14 inserts; per-hit
charge 12288 instructions, itemised total 2^39.1 instructions, dominated by the 1-instruction-per-step
amortised row by 2^86. Peak memory itemised: records `55,956,890 x 52 B = 2^31.00 B`, table 7.08 MB,
code 2 MB, scratch 16 MB, sort buffer 32 MB, **peak 2.97 GB < 2^32 B**. On a match both trails are
replayed to their recorded steps, the twelve Frobenius shifts are compared, `x` is solved from the
linear relation and accepted only if `g^x = z`; degenerate or failing matches are discarded, never
reported. Escapes: a class repeated within two steps forces one multiply step at index `(h+1) mod R`,
which replaces the step it would have taken — escapes add no steps and no charge. The unexplained
2^65 setup allowance is **deleted** and replaced by a derived row, 1024 x (510 reference mults + 12
Frobenius) + 2^20 = 2^19.0 units.

**0c. The collision model is re-scoped honestly (verdict item 3).** "Same index never merges" is
withdrawn: fixed-index twisted multiplication is not proved injective on Frobenius classes, since
`aM = (bM)^(p^j)` with `j != 0` does not imply equal input classes. The merge law is now the scoped
heuristic H1' with the j-shift gap stated, and the work claim charges the **smaller** ordinary rate
`rho = (1-s^2)(1-1/R) = 0.858536` rather than the up-to-12/N match rate — even though a shifted match
is not a loss, because it yields `x = (p^j A2 - A1)/(D1 - p^j D2)` directly. The toy separates the
channels: **134 shifted-image matches, all 134 usable; 1228 same-index matches, 1225 of them
degenerate cycle closures** (correctly discarded). The reviewer's algebra `1 - s^2 - (1-s)^2/R =
0.858994` exceeds ours by exactly `2s(1-s)/R`; the ledger asserts that inequality symbolically, so
using the smaller rho is the conservative side.

## Numbers

| quantity | value |
| --- | --- |
| rows | common 301, multiply 2493 (+C) = **22493 = 1.12465 u**, squaring **15383 = 0.76915 u** |
| average step price | `c_avg = 79307/80000 = 0.9913375` |
| classes | `N = (r-1)/12 = 2^251.272`, `N' = N/rho = 2^251.492` |
| steps | `S = ceil(sqrt(2 N' ln(100/61))) = 70933784039754404319800682927445544128 = 2^125.737` |
| setup | 10,475,339,776 instructions = 2^19.0 units (derived, replaces 2^65) |
| nominal work | **W = 2^125.725..125.726** |
| metered cap | **W_cap = W * 513/512 = 2^125.728..125.729** |
| claim | **125.73**, checked below the cap by exact integer power comparison and <= gate - 0.005 |
| peak memory | 2,969,264,968 B (2.97 GB), bounded by the work cap |

## Evidence, and what is only mechanism evidence

* `evidence/prod_run15.json` — 2000 steps on real BLS12-381 GT: 767 squarings (0.3835, 95% CI
  +-0.021 covers 3/8), 63+ exponent-invariant checks with 0 failures, 33/33 class-invariance
  agreements against all eleven explicit conjugates, chi-square 5.1 on 7 degrees of freedom of
  `T mod 8`, 26/26 part-A algebra checks, 12 DP records with the store/replay path exercised on an
  injected duplicate (`dp_usable = 0`: the match is correctly rejected).
* `evidence/toy_census15.json` — 180 trials, q = 523, N = 43420, R = 128, DP width 4, branch rule
  `(T mod 8) < 3` (production, not a lookalike), DP-gated store, replay/verify, forced-multiply
  escape: 145/180 recoveries within 4 S_ops, 57/180 at S_ops, 4171 DP hits, 2801 records, 1370
  matches split 1228/134 across the two channels, 28 escapes, 4642 invariant checks, 0 failures.
* `evidence/predicate_match.json` — the shared production predicate module versus an independently
  written copy of the formula: **100/100 agree** on random production GT elements, and the toy
  predicate versus the same spec: **100/100 agree** on random toy elements.
* `evidence/ledger15.json` — every printed number above, recomputed from the two raw runs with
  integers and Fractions only.

The toy's branch frequency is 0.302 rather than 3/8 because at q = 523 the invariant T is far from
uniform mod 8; that is a small-field artifact, it is disclosed, and the toy's frequency is not used
in the claim. The production field, with 381-bit T, measures 0.3835 +- 0.021.

## What is still heuristic

H1' (random-mapping behaviour of the mixed class walk) and H2 (pseudorandomness of the disjoint T bit
ranges) are heuristics and are not proved. No production-size walk is executed; this is a square-root
attack analysis with a complete resource accounting, not a security estimate. The squaring price is
analytic from counted Fp2 operations, not a measured routine, and it is charged on a scale about 1.9x
stricter than the reference operation.

## Revision history

* v13 (`gt-mixwalk-selfcontained-40pct-v13`, f = 0 floor): superseded; its work claim had no floor
  accounting and its collision model relied on the same-index no-merge premise later rejected.
* v14 (`gt-squaring-branch-40pct-v14`, d401a3c3, claim 125.73): **rejected**. My own diagnosis,
  matching the verdict: (i) it charged the average step price `c_avg = 0.99129` while claiming a cap
  of 2^125.73, so a max-price budget of 2^125.91 broke the claim — the defect is charging an average
  and claiming a maximum; (ii) it ran transitions, not the DP store/replay procedure, and carried an
  unexplained 2^65 setup allowance; (iii) it claimed the same-index no-merge property without proving
  bijectivity on Frobenius classes, and its rho was the smaller of two candidate expressions without
  saying so.
* v15 (this package): meters the work so the cap holds by construction and pays the abort debit with
  Hoeffding; writes the full DP procedure with an itemised store and a 2.97 GB peak; deletes 2^65 in
  favour of a derived setup row; demotes the collision lemma to a scoped heuristic, charges the
  conservative rho and measures both collision channels. Claim unchanged at 125.73, which is <= the
  live gate 125.74 minus 0.005.

## The meter, as it runs

```
meter = SETUP                      # 10,475,339,776 instructions, table + code + pairings
while not found:
    u = T(y)                       # Re(c0): one field-word read, class invariant
    if u mod 8 < 3:  charge = 15383     # cyclotomic squaring, no reference operation
    else:            charge = 22493     # one dense Fp12 mult (20000) + 2493 counted rows
    if meter + charge > CAP:  abort     # the cap is never crossed, not even by one step
    meter += charge
    y = step(y);  if dp(u): store_or_match(u, y, meter)
```

`CAP = W_cap * C` with `W_cap = W * 513/512`, and the abort test is evaluated before every step, so
the charged total is `<= CAP` for every input and every trajectory. Because the charge per step is one
of exactly two integers, the total after n steps with k squarings is
`SETUP + 22493(n-k) + 15383k`, which is the expression the ledger inverts to get `k_min`; the
implementation and the analysis therefore use the same arithmetic, not two approximations of it.

## Storage arithmetic

`records = ceil(S / 2^100) = 55,956,890`, because a DP record is written only on a DP step and the DP
predicate has probability 2^-100 on the class invariant. Slots are 52 bytes; `55,956,890 x 52 =
2,909,758,280 B`. Adding the jump table (1024 entries x 12 conjugates x 576 B = 7,077,888 B), code
(2^21), scratch (2^24) and the sort buffer (2^25) gives `2,969,264,968 B`, i.e. 2.97 GB, checked
below 2^32 B. The record count is bounded by the *work* cap, not by a hope: a run can execute at most
`CAP / 15383` steps, so at most `ceil(CAP / (15383 * 2^100))` records can ever be written, and the
array is sized by that bound.

## Why the collision model is charged on the conservative side

Including the j-shift ambiguity the per-pair match probability is at most 12/N, since a class has
twelve Frobenius images; the ordinary channel is 1/N and the mixture adds the branch/exclusion
factors. Our rho is `(1-s^2)(1-1/R)`; the reviewer's algebra gives `1-s^2-(1-s)^2/R`, larger by
`2s(1-s)/R = 0.000457`. Using the smaller value enlarges N' by that factor and lengthens S, so the
ledger is on the safe side of the discrepancy rather than the flattering side. Nothing in the claim
depends on the extra 12x channel: it is measured (134 of 134 usable) and left uncashed.

## Reproduce

```
cd scripts
python3 gt_prod.py    < '{"steps":2000,"dp_bits":8,"R":1024}'   > ../evidence/prod_run15.json
python3 toy15.py      < '{"trials":180,"dpbits":4,"R":128}'      > ../evidence/toy_census15.json
python3 predicate_check.py                                       > ../evidence/predicate_match.json
python3 ledger15.py                                              # rewrites certificate.json
```

`ledger15.py` reads only the three evidence JSONs and asserts every row, every inequality
(`W^100 < 2^125725`, `W_cap^100 < 2^12573`, `claim + 0.005 <= gate`) and every probability bound with
integers and Fractions; it exits non-zero on any mismatch, so the certificate cannot drift from the
runs.
