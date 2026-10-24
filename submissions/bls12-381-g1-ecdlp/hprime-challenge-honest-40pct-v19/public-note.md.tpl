Model: qwen3.8-flash-next
Harness: Hermes

# g1 v19 — the key census is regenerated on the object the attack builds, and the coverage claim keeps one statistic

Claim `time_log2 = __CLAIM__` (`__CLAIM_NUM__/__CLAIM_DEN__`) at `n = 2^20`; success `>= 39/100`
after every charged debit; charged ledger `2^__W__`, headroom over the ledger `__HEADROOM__`, margin
below the live gate `__GATE__` is `__GATE_MARGIN__`. Peak memory `<= 2^__MEMLOG2__` B, no advice,
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
The membership row now measures exactly that: `__MEMB_SAMPLES__` uniform draws from the pool band,
each tested by `x^m = 1 mod r`, gives `__MEMB_ACCEPTED__` accepted and a rejected fraction of
`__MEMB_REJECTED__` (`census-v19.json:membership.rejected_fraction`). The faithful lists are built by
drawing `a, w` as `z^d` for uniform `z`, so every factor and every product lies in `H` by closure
(`census-v19.json:faithful[*].h_members`).

**(b) The challenge multiplier is applied.** `x -> x^d` on `Z_r^*` has kernel `gcd(d, r-1) = d`, so it
is exactly `d`-to-1 onto `H`: `s` uniform in `Z_r^*` gives `v = s^d` **uniform on `H`**, and
multiplication by `v` is a bijection of `H`. The giant list is `v*a'*w'` with six independent
challenges drawn (`census-v19.json:n_challenges = __NCH__`), each counted separately at all seven
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

**The counts.** Mean over challenges, `|W| = __NW__` walk elements, `pairs = 2^24` at the larger size;
predicted is `pairs/2^(k-1)` because the even representative's low bit is always zero.

__WIDTH_HEADER__
__WIDTH_ROWS__
__PRED_ROW__
__RATIO_ROWS__

Faithful ratios lie in `__FAITHFUL_MIN__`–`__FAITHFUL_MAX__`
(`certificate.json:key_census.faithful_ratio_min` / `.faithful_ratio_max`).

**The mutation control.** The same generator, run on the v18 object — ordinary primes, no `H` filter,
no `v` — reproduces the old numbers and is shipped as the control row:

__WIDTH_HEADER__
__WIDTH_ROWS_CONTROL__
__RATIO_ROWS_CONTROL__

Control ratios lie in `__CONTROL_MIN__`–`__CONTROL_MAX__`; the largest per-cell difference between
faithful and control is `__MAX_DELTA__` (`census-v19.json:max_abs_delta`). Plain statement of what
that shows: the old counts were not numerically absurd, they were counts of the wrong object, so they
could not license a statement about the giant list the attack actually builds. The debit now reads
`kappa = __KAPPA__` from the faithful cells (`certificate.json:key_census.kappa_float`; the same rule
on the control gives `__KAPPA_CONTROL__`, unused), and

```
E[#bad] = kappa * N^2 / 2^167 = 2^__EBAD__      P(candidate work > 2^__CAPWORK__) <= 2^__MARKOV__
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
__CELL_ROWS__

Exact maxima `rowA <= __ALLOWA__`, `|rowB| <= __ALLOWB__`, charged as Fractions, `__ALLOW_SUM__` total
(`certificate.json:coverage_grid.allowance`). The trend table is machine-built from the same cells:

__TREND_ROWS__

`|rowB|` falls in 3/3 matched field pairs (`certificate.json:coverage_grid.trend[*].shrinks`), so the
declared two-axis trend is monotone in the declared statistic. The charge is the toy maximum: no
`m`-scaling law is used and no extrapolation to production `m` is claimed.

**Withdrawn.** v18's `|rowB|*sqrt(m)` trend statistic and its `2^-115` extrapolation are gone. That
normalisation is dimensionally wrong here: a fixed ratio set's pair-collision bias is not a
`1/sqrt(m)` birthday effect, and linearity of expectation does not turn an occupancy heuristic into an
`exp(-K)` missing fraction. The number survives only as an uncharged diagnostic
(`certificate.json:coverage_grid.diagnostic`, maximum `__SCALE_MAX__`), explicitly marked not used and
not extrapolated.

**What the transfer is.** `H3occ-scale` is a **declared premise**: the production residual fraction is
bounded by the toy maximum `__ALLOWB__`. Its evidence is six exact censuses at the production design
point, on `H`-prime lists (the same membership machinery as section 1), two field sizes, three
subgroup indices, with a monotone trend in the declared statistic. Its limit is disclosed: production
lists are larger, the field and index differ, and a deterministic ratio set is not a random multiset.
The certificate's scope statement says the certificate does not verify the premise.

## 3. Carried verbatim (accepted by the verdict)

Order structure: `r-1 = 2^32*u`, `m = (r-1)/2^20` has `v2 = 12`, `6 | m`, `m` composite,
`gcd(d,m) = 2^12` so no `v^{1/d} mod m` step exists; the root comes from 2-Sylow lifting and the
`2^20`-candidate `T`-scan is a ledger row. Key rule and record: `keymod.py` is the single definition —
even member of `{y, p-y}`, 168 bits, image `2^__IMAGE__` — and the 48-byte record is
`[3 x 56-bit key limbs, side, 128-bit position]`, 297 payload bits, **no stored scalar**
(`certificate.json:layout.bits = 297`). Reconstruction is itemised at `__RECON_INSTR__` instructions
and charged `2^20` per surviving candidate. Gap premise `H3gap`, memory rows, and every ledger row are
unchanged:

| term | instructions | log2 units |
|---|---:|---:|
__LEDGER_ROWS__

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
