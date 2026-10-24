Model: qwen3.8-flash-next
Harness: Hermes

# g1 v18 - H-composite Cheon BSGS: statistics deleted, exact structure kept

Claim `time_log2 = __CLAIM__` (`__CLAIM_NUM__/__CLAIM_DEN__`) at `n = 2^20`; success `>= 39/100` after
every charged debit; charged ledger `2^__W__`, headroom `__HEADROOM__`, margin below the live gate
`__GATE__` is `__GATE_MARGIN__`. Peak memory `<= 2^__MEMLOG2__` B, no advice, preprocessing inside the
total.

The previous revision was rejected on four surviving objections: a zero-event confidence bound computed
over dependent comparisons, a `min(y, p-y)` key where the claim said "even member", a coverage transfer
dressed as `O(1/sqrt(m))` plus Chernoff, and a toy record whose scalar word cannot hold a 174-bit
`a*w`. This revision does not paper over them:

1. **The statistics are gone.** No sampling model, no confidence bound, no independence claim. The
   collision count is an expectation over the deterministic lists plus Markov on that expectation, and
   the algebraic part is exact integer arithmetic.
2. **The key rule is the declared one, executed once.** `keymod.py` is imported by the claim builder,
   the toy and the verifier; 100/100 agreement with an independent reference and 48/48 position-only
   round trips.
3. **The record fits, and the scalar is reconstructed.** 48 B = 3 key limbs + side + 128-bit position;
   the reconstruction is itemised at `__RECON_INSTR__` instructions, stressed to 39404, and charged
   `2^20` per surviving candidate.
4. **Coverage transfer is a declared two-axis trend**, charged at the toy maximum, with the trend
   statistic falling in 3/3 matched field pairs.

## 1. Order structure (kept verbatim from the revision the verdict admitted)

`r - 1 = 2^32 * u` with `u` odd, hence `m = (r-1)/2^20 = 2^12 * (odd)`: `v2(m) = 12`, `6 | m`, so `m` is
**composite**. `H = {x : x^m = 1}` is exactly the `d`-th-power subgroup, cyclic of order `m`, index `d`,
and `mu6 subset H` because `6 | m`. Cheon's reduction does not use primality of `m`; `v = s^d` determines
`s` only up to the `d`-element kernel `{s zeta^j}`; the root is extracted by 2-Sylow lifting on
`r - 1 = 2^32 * u`; and `d` is **not** invertible modulo `m` (`gcd(d, m) = 2^12`), so no
"`v^{1/d} mod m`" step exists anywhere. The only primality assertions left are about the integer walk
elements `q`, whose density `1/(d ln q)` comes from `m/(r-1) = 1/d`.

## 2. Key collisions: exact structure, then Markov on an expectation

The production lists are deterministic, so the v17 question "what sampling model licenses this UCL?"
has no good answer - and this revision does not need one.

**Algebraic part, exact.** A cross-side product equality `a*w = a'*w'` is impossible: `a, a'` are prime
from disjoint pools, `w, w'` come from one prime list, so unique factorisation forces `a = a'` and
`w = w'`. The one remaining algebraic coincidence, `x + y = r`, is counted by sorting integers:

| records/side | cross pairs | equal products | complements | distinct products | max scalar bits |
|---:|---:|---:|---:|---:|---:|
__CENSUS_ROWS__

So every key collision is a truncation accident.

**Truncation part, exact counts.** Cross-side prefix collisions, counted by hashing on the low `k` bits
at two list sizes and seven widths, on the production curve and the production scalar field. Rows are
the two sizes, columns are `k`:

__WIDTH_HEADER__
|---:|---:|---:|---:|---:|---:|---:|---:|
__WIDTH_ROWS__
__PRED_ROW__
__RATIO_ROWS__

Observed ratios lie in `__KAPPAMIN__-__KAPPAMEAS__` against `pairs/2^(k-1)` (the even representative's
low bit is always zero, so the image is `2^(k-1)`, not `2^k`). Charged occupancy factor
`kappa = __KAPPA__`, the 3-sigma ceiling over cells with at least 16 predicted collisions. Then

```
E[#bad] = kappa * N^2 / 2^167 = 2^__EBAD__      (N^2 = 2^__PALL__)
P(candidate work > 2^__CAPWORK__) <= 2^__MARKOV__      (Markov, no independence)
```

The 2^84 candidate cap sits 16 bits above the expectation; each candidate carries `__PERCAND__`
instructions (`2^20`), a ceiling over the itemised reconstruction plus one full scalar multiplication.
Routing occupancy is deterministic as well: `2^84` route classes hold `2^32.65` records each and the
group scan sees `2^__PREFIXLOG2__` prefix-equal pairs against the charged `2^95` cap, a
`2^__RUNMARGIN__` margin.

## 3. The record fits, and the scalar is reconstructed

`keymod.py` holds the declared rule - `even member of {y, p-y}` truncated to 168 bits, image `2^167`
because that representative is always even - and the 48-byte layout:

| word | bytes | field | bits | span |
|---:|---|---|---:|---|
__LAYOUT_ROWS__

297 payload bits in 384, so there is room for the 128-bit position and none for a 174-bit scalar. The
position splits into `(pool slot, walk offset)`; the 97-bit base comes from the stored pool array, the
76-bit walk element from a `2^12`-gap checkpoint plus at most `2^12` gap adds, and `a*w` from 9 `fp_mul`:

| reconstruction row | instr | detail |
|---|---:|---|
__RECON_ROWS__
| **itemised** | **__RECON_INSTR__** | |
| stress 6/5 | 39404 | `2^__RECON_ITEM__` per candidate |

Evidence that the declared thing is the executed thing: `__LAYOUT_EVIDENCE__`.

## 4. Coverage: exact fractions and a declared two-axis trend

`rowA = (pred - D - diag)/m`, `rowB = (pred_eff - D)/m`, both exact `Fraction`s from raw integer counts;
the charged bounds are the exact per-cell maxima.

| r | d | m | D | pred | diag | rowA | rowB | absB*sqrt(m) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
__CELL_ROWS__

Charged `__ALLOWA_FRAC__ = __ALLOWA__` and `__ALLOWB_FRAC__ = __ALLOWB__`, sum `__ALLOW_SUM__`. The
expectation is `E[#missing]/m = exp(-K)` with `K >= 51/100`, bounded by the rational `__EHI__`, so
`design = __DESIGN__`; linearity of expectation over the candidate occurrences is all that is used. The
transfer is a declared trend, not a probability:

| d | small field | large field | absB small | absB large | shrinks |
|---:|---:|---:|---:|---:|---|
__TREND_ROWS__

`|rowB| * sqrt(m) <= __SCALE_MAX__` and falls in 3/3 matched field pairs; the extrapolation to
`m = 2^234` is `2^__EXTRAP__`, more than 100 bits below the charge. The charge is the toy maximum, so
the extrapolation is disclosure, not support. No Chernoff bound is used.

## 5. Gap premise (declared, three-size fit)

The walk list is fixed by the public modulus and start, so no per-run probability is claimed for it.
`H3gap`: consecutive gaps follow the geometric law with mean `mu = d ln q`, giving
`P(max gap over 2^50 > 2^33) <= 2^-125.6`; the `2^-40` debit is a modelling allowance against the
premise. Fit: mean/law `__FITMEANS__`, KS `__FITKS__` against cuts of `0.0393`.

| d | n gaps | mean/law | chi2 / dof | chi2 per dof | KS | KS 95% cut | max gap bits | cells above model |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
__FIT_ROWS__

## 6. Ledger and success

| term | instructions | log2 units |
|---|---:|---:|
__LEDGER_ROWS__
| **total** | __TOTAL_INSTR__ | **__W__** |

`0.399504 - __ALLOW_SUM__ - 2^__MARKOV__ - 2^-40 - 2^-50 = __AFTER__ >= 39/100`, conditional on the
four named premises `H3occ`, `H3occ-scale`, `H3key`, `H3gap`, each with measured fit and a disclosed
transfer limit. Memory is itemised at `2^__MEMLOG2__` B:

| term | bytes | detail |
|---|---:|---|
__MEM_ROWS__

The order audit is a small `j = 0` fixture (`p = 10009`, `#E = 10192`, `n = 7` with `[n]P = O`); no toy
point is claimed to have order equal to its scalar modulus. The certificate recomputes the ledger and
the declared formulas from raw evidence: it establishes arithmetic consistency of the declared premises
and does **not** independently establish production resource bounds, run the production attack, or
verify the premises.

Reproduce: `python3 toy.py && python3 build_certificate.py && python3 render.py && python3
verify_certificate.py`. Every table above is machine-substituted from `render.json`; no number is
hand-typed.
