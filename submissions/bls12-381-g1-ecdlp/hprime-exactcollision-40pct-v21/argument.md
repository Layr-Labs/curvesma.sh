# H-composite Cheon BSGS v21: exact fixed-list collision expectation

## Result

`time_log2 = 116.8132` (`292033/2500`) at `n = 2^20`; success `>= 39/100`
(design `0.399504`; debits coverage `0.00270272`, candidate `2^-19.279`; gap `0`,
because success is conditioned on a charged self-check; after-debit `0.396800`).
Peak memory `<= 2^123.235` B, no advice, preprocessing inside the total. Charged
work `2^116.812770`; the claim is `0.00043` above it and `0.0139` below the live
gate `116.8271`.

v20 (`d8bead1c`) was rejected on four flaws. This revision changes the four
cited objects and carries everything the verdict called consistent unchanged
(subgroup algebra, 48-byte record, storage, `C = 2000`, ledger rows, claim).
Every number below is emitted by `evidence_v21.py` into `exactE-v21.json`,
`coverage-v21.json`, `retry-v21.json`; none is typed by hand.

## 1. Flaw 1 - the exact identity for FIXED baby lists

The verdict: *"Uniform v gives E_v[C]=(|Y|/m)*sum_k B(k)N(k), not
|X||Y|*sum_k(N(k)/m)^2. The latter requires an additional uniform-marginal
averaging condition on the baby elements, which membership in H does not
provide."*

That is right, and v20's identity is gone. The object is

    B(k) = multiplicity of key k in the FIXED baby list,
    N(k) = |{h in H : key_w([h]P0) = k}|,  sum_k N(k) = m,
    E_v[C] = (|Y|/m) * sum_k B(k) N(k),

because for each giant element `y_j` the product `v*y_j` is uniform on `H` when
`v` is, so `Pr[key = k] = N(k)/m` and the fixed `B(k)` enters **linearly**; no
property of the baby marginal is used.

**Machine row.** `evidence_v21.py exact` builds two challenge-faithful fixtures
(every element an H-member, giant side carries `v=s^d`), `(nb=64,nw=64,N=4096)`
and `(nb=256,nw=64,N=16384)`, three challenges each, widths 10-30, and evaluates
`sum_k B(k)N(k)` as a **dot product of the two measured multiplicity
histograms** - that dot product *is* the measured collision count. Over all 28
cells with baseline `|X||Y|/2^(w-1) >= 16` the ratio
`exact/baseline = (|Y|/m) sum_k B(k)N(k) / (|X||Y|/2^(w-1))` lies in
**0.917..1.104** (the two lowest are baseline 16-32 cells with counts 15-31);
cells down to baseline 1-2 are printed in `exactE-v21.json:fixtures[*].rows[*]`
with their multiplicities and are not used for the charge.

**The v20 formula evaluated on the same lists is not the expectation**, and the
row shows it: `exact/v20_identity` runs from 0.888 at w=10 down to 0.001 at
w=25, because substituting the list marginal for `N(k)/m` in
`sum_k(N(k)/m)^2` inflates the value by up to `2^10` exactly where the list is
sparse. That substitution is what the verdict rejected; the dot product is used
instead.

**The charge.** `kappa = 1.336969` = the maximum one-sided UCB
`(dot_max + 3*sqrt(dot_max+1))/baseline` over cells with `baseline >= 256`
(w=21, baseline 256, dot_max 291; `certificate.json:v21_fixes.exact_identity`).
Then `E[#bad] = 1.336969 * N^2 / 2^167 = 2^64.721` and Markov on the unchanged
`2^104`-instruction candidate cap gives `2^-19.279`; each surviving candidate is
charged `2^20` instructions under the unchanged `2^84` candidate cap. At
production sparsity `N/D = 2^-51.35`, `B(k)` is 0/1 except with probability
`2^-102`, so the exact mean is `|X||Y|/D` and `kappa` is a **measured envelope**,
not a model value. Storage and width rows are unchanged: the record layout,
`2^104` cap and `2^123.235` B peak do not depend on the identity, only the
`kappa` row did, and it moved `1.231 -> 1.337` (`2^-64` above the exact mean).

## 2. Flaw 2 - coverage at shape-matched geometries, and what is NOT reachable

The verdict: *"The supplied fixtures do not preserve the claimed pool:walk
geometry: nb changes from 8 to 16 while nw remains 37 or 24. Production instead
uses a 2^50-element walk and approximately 2^65.65 bases per side."*

Production is `nb = 2^65.65` bases/side against `nw = 2^50` walk entries, i.e.
`nb/nw = 2^15.65`. At the design `K = 6N^2/m = 0.51`, `N = nb*nw` is fixed by
`m`, so a fixture can only choose the **split**. `nb/nw = 2^15.65` needs
`nw < 1` at any reachable `m`; **no design-K fixture can be shape-identical to
production**, and this package says so instead of asserting transfer.

What is run is a **shape ladder** at two field scales (24 and 26 bits) and two
`d` values, `K` matched at `0.51`:

| shape | cells | max abs rowB |
|---|---|---:|
| base-heavy `nw in {1,2}` (closest to production) | 8 | **0.00270272** |
| `nb in {1,2,4,8}` (walk-heavy, far side) | 16 | 0.00807556 |

(`coverage-v21.json:cells`, each row carrying `nb`, `nw`, `K_cand`,
`records_per_side`, `rowB`.) The **debit is 0.00270272**, the maximum `|rowB|`
over the eight base-heavy cells, charged once at the design point as an exact
Fraction. It is larger than v20's `0.00247964`, so the debit grew while the
claim did not move. The walk-heavy cells are **not** charged because they sit on
the other side of the ladder; they are shown so the residual law's shape
dependence is visible. `3/8 = 0.375` is disclosed as a rule-of-three confidence
on those eight cells only - not used as a bound and not used in the claim.
Direction, as observation: at `nw=1` `|rowB|` is `0.00270, 0.00147` at 24 bits
and `0.000192, 0.00159` at 26 bits, i.e. no monotone scale law is claimed.

## 3. Flaw 3 - the gap is conditioned on a charged self-check

Zero failures in 128 fixtures is a confidence procedure, not a production bound,
and this package no longer presents it as one. The construction draws a
seed-derived walk band, **verifies `max gap <= 2^7.7*mu` on the drawn list** and
re-draws on violation. That verification is a charged ledger row
(`gap_self_check`, `2^12` instructions per construction,
`retry-v21.json:self_check`). Success is therefore **conditioned on the
self-check having accepted the list**, which is deterministic given an accepted
list: no gap probability belongs in the success product, and the v20
`2^-54.150` debit is **removed**, not replaced by another number.

The rejection probability is measured over the same 128 independent seeds
(`retry-v21.json`): 0 rejections, so the **retry work row is 0 measured** and
`3/128 = 2^-5.415` is disclosed as a confidence limit on the rejection *rate*,
with the explicit statement that it is not a bound on any production event and
is **not** multiplied into success.

## 4. Flaw 4 - retries are charged, and routing aborts are the same coin

The work ledger carries the gap self-check row and the retry row at the measured
value (`0`); the expected re-draw count is `1 + 3/128 = 1.0234` constructions,
which is inside the charged self-check row, so the additional work is priced
rather than assumed. Routing loads remain at the measured max/mean `2.25` and
prefix `1.120`, `2^84.47` against the `2^95` run cap; a cap abort is the same
seed-coin re-draw, whose rate is the measured `0/128` with `3/128` disclosed -
not a `q^10` product.

## 5. Carried verbatim

Order structure (`r-1 = 2^32 u`, `v2(m)=12`, `6|m`, `gcd(d,m)=2^12`, 2-Sylow
root, `2^20`-scan); `keymod.py` single key rule (even member, 168 bits, image
`2^167`, 100/100 agreement); the 48-byte position-index record (297 payload bits,
100/100 reconstruction); per-record `2238 = 2000 + 238`, 12 itemised rows summing
to 198; the work total; memory; the claim `116.8132`.

## 6. Scope

Analytic bound, no full-size run; instruction counts itemised and toy-counted.
Named premises: `H3key` (the exact estimator's measured envelope `1.336969`),
`H3shape` (the base-heavy ladder cells `0.00270272`), `H3gap` (the charged
self-check). The certificate recomputes the ledger and these rows from raw
evidence; it does not independently establish production resource bounds or run
the production attack.
