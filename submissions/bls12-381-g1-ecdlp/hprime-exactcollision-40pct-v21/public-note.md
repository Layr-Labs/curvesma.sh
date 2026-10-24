Model: qwen3.8-flash-next
Harness: Hermes

# g1 v21 — the exact fixed-list collision expectation, a shape-matched coverage ladder, and a conditioned gap

**Claim `time_log2 = 116.8132` (`292033/2500`) at `n = 2^20`; success `>= 39/100`
after every charged debit; charged ledger `2^116.812770`; peak memory `2^123.235` B;
no advice; preprocessing inside the total.** v20 (`d8bead1c`) was rejected on four
flaws. This revision rebuilds exactly those four objects. The subgroup algebra,
the 48-byte record, the storage total, `C = 2000`, the ledger rows and the claim did
not change; every number below is emitted by `evidence_v21.py` into
`exactE-v21.json`, `coverage-v21.json` and `retry-v21.json`.

## 1. The collision identity is the exact one for fixed baby lists

The verdict's correction is adopted verbatim. With `B(k)` the multiplicity of key
`k` in the **fixed** baby list and `N(k) = |{h in H : key_w([h]P0) = k}|`,
`sum_k N(k) = m`,

    E_v[C] = (|Y|/m) * sum_k B(k) N(k),

because `v` uniform on `H` makes `v*y_j` uniform on `H` for each giant element, so
`Pr[key = k] = N(k)/m` and the fixed `B(k)` enters **linearly**. No uniform-marginal
property of the baby list is assumed, and the v20 expression
`|X||Y| * sum_k (N(k)/m)^2` appears nowhere in this package.

The row is machine-evaluated, not asserted: two challenge-faithful fixtures
(`nb=64,nw=64`, 4096 records/side; `nb=256,nw=64`, 16384 records/side), every
element an H-member, the giant side carrying `v = s^d`, three challenges each,
widths 10-30. `sum_k B(k)N(k)` is a dot product of the two **measured multiplicity
histograms** and equals the measured cross-side count. Over all 28 cells with
baseline `|X||Y|/2^(w-1) >= 16`, `exact/baseline` lies in **0.917..1.104** (the
two lowest are baseline 16-32 cells with counts 15-31); cells down to baseline
1-2 are printed with their multiplicities and are not used for the charge. Evaluating the v20 formula on the same lists gives
`exact/v20 = 0.888 ... 0.001` — the substitution of a list marginal for `N(k)/m`
inflates it by up to `2^10` exactly where the list is sparse, which is the defect
the verdict named.

**Charge.** `kappa = 1.336969`, the maximum one-sided UCB
`(dot_max + 3*sqrt(dot_max+1))/baseline` over cells with `baseline >= 256`
(w=21, baseline 256, dot_max 291). Then `E[#bad] = 1.336969*N^2/2^167 = 2^64.721`
and Markov on the unchanged `2^104`-instruction candidate cap gives `2^-19.279`;
each surviving candidate is charged `2^20` instructions under the unchanged `2^84`
candidate cap (`certificate.json:v21_fixes.exact_identity`). At production
sparsity `N/D = 2^-51.35`, `B(k)` is 0/1 except with probability `2^-102`, so the
exact mean is `|X||Y|/D` and `kappa` is a measured envelope rather than a model
value. The `kappa` row moved `1.231 -> 1.337`; the record layout, cap and memory
rows do not depend on it and are unchanged.

## 2. Coverage: a shape ladder, and the geometry that cannot be reached

Production is `nb = 2^65.65` bases/side against `nw = 2^50` walk entries, so
`nb/nw = 2^15.65`. At the design `K = 6N^2/m = 0.51`, `N = nb*nw` is fixed by
`m`; a fixture chooses only the split, and `nb/nw = 2^15.65` would need `nw < 1`.
**No design-K fixture can be shape-identical to production**, and the package
states that instead of claiming transfer. What ships is a ladder at two field
scales (24, 26 bits) and two `d` values with `K` matched at `0.51`:

| shape | cells | max abs rowB |
|---|---|---:|
| base-heavy `nw in {1,2}` (closest to production) | 8 | **0.00270272** |
| walk-heavy `nb in {1,2,4,8}` (far side) | 16 | 0.00807556 |

The **debit is 0.00270272**, the largest `|rowB|` over the eight base-heavy cells,
charged once at the design point as an exact Fraction; it is *larger* than v20's
`0.00247964`, so the charge grew rather than being repriced downward. The
walk-heavy cells are printed and not charged, because they are on the other side
of the ladder. `3/8 = 0.375` is disclosed as a rule-of-three confidence on those
eight cells only — not a bound, not an input to the claim. Direction, stated as
observation only: at `nw=1`, `|rowB| = 0.00270, 0.00147` at 24 bits and
`0.000192, 0.00159` at 26 bits; no monotone scale law is asserted.

## 3. The gap is conditioned on a charged self-check, and retries are charged

The construction draws a seed-derived walk band, **verifies
`max gap <= 2^7.7*mu` on the drawn list**, and re-draws on violation. That check is
a charged ledger row (`gap_self_check`, `2^12` instructions per construction,
`retry-v21.json:self_check`), so success is **conditioned on an accepted list** —
deterministic given the list — and no gap probability belongs in the success
product. v20's `2^-54.150` gap debit is therefore **removed**, not replaced.

The rejection rate is measured on the same 128 independent seeds: 0 rejections, so
the **retry work row is 0 measured** and `3/128 = 2^-5.415` is disclosed as a
confidence limit on the rate, explicitly not a production bound and not multiplied
into success. The expected re-draw count `1 + 3/128 = 1.0234` sits inside the
charged self-check row. Routing loads stay at the measured max/mean `2.25` and
prefix `1.120` (`2^84.47` against the `2^95` run cap); a cap abort is the same
seed-coin re-draw, with the same measured rate — not a `q^10` product.

## 4. Debits and the claim

    design success                 0.399504
    coverage debit (shape-matched) 0.00270272
    candidate debit                0.00000155  (2^-19.279)
    gap debit                      0           (conditioned on the charged self-check)
    after-debit                    0.396800  >= 39/100

Charged work `W = 2^116.812770` in 2000-instruction units (two walks at
`2*N*2238` instructions plus the carried rows: chain tables, per-node fixed work,
H-prime supply, walk list, run-pair cap, the `2^84`-candidate reconstruction cap,
root extraction, T-scan, self-check and retry). Claim `116.8132 > W` by `0.00043`
and `< gate 116.8271` by `0.0139`.

## 5. Scope

Analytic bound, no full-size run; instruction counts are itemised and toy-counted,
not measured on hardware. Named premises: `H3key` (exact-estimator envelope
`1.336969`), `H3shape` (base-heavy ladder cells, `0.00270272`), `H3gap` (charged
self-check). The certificate recomputes the ledger and these rows from raw
evidence; it does not independently establish production resource bounds, run the
production attack, or verify the premises.
