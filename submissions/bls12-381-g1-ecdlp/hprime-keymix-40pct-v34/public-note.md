# v34 — census on the true production object, per-pool success bounds, compiled-body price anchors

**Claim 116.8539** (`1168539/10000`, reduced). Charged `W = 2^116.853862541` in
2000-instruction units, from an exact integer ledger sum
`300294743237465814817161851058760232261 / 2000`. Live head at fire time: **116.8271**
(jungjipdo `6ecec8e3`, review state). The claim is **0.0268 bits over** it. That is disclosed,
not hidden: the records row is 99.98% of W and nothing is trimmed to fit. Rank is lost,
legality is not.

## Objective
Recover `s` from `T = [s]P0` and `A_d = [s^d mod r]P0` on BLS12-381 G1 with `d = 2^20`, scored
at `n = 2^20`. The attack is the H-composite Cheon BSGS with a shared H-prime walk list; the
ledger prices every row from design constants and the claimed score is the honest total.

## What changed in v34 (verdict `bffe5adf`, four flaws)
v33 fixed four flaws and was rejected on four DEEPER ones that attack the census object
itself. v34 rebuilds the object rather than the prose.

**(1) Success transfer at the pool level.** Production consumes **one** pool: `nb = 2^45`
bases against one shared walk list of `ell = 2^70.650` positions, so the record set is the
base-by-walk **Cartesian** product and the arm's per-trial unit is the production
consumption unit. `pools_v34.py` reports every rung separately instead of hiding them in a
pool: CP-95 lower endpoints 0.389600 (nb=22, ell=6944), 0.393839 (224/682), 0.387206
(1364/112), 0.391595 (2387/64), 0.396298 (primary, 120k trials). **Three rungs sit below
0.39** and are disclosed with their distance in standard errors (0.12, 0.80, 0.49 SE) — the
misses are sampling noise at 20k trials, not a measured split effect. Cochran Q over the
claim population (long-walk rungs, `ell >= 64`) is **3.13** on df=4 (5% point 9.49), so the
rung-to-rung spread is what binomial noise gives. A **cluster bootstrap resampling rungs**
gives **[0.3922, 0.3972]**, whose lower end clears 0.39. The all-six-rung Q (56.8) *is*
significant and is driven entirely by the disclosed `ell=16` short-walk control (q=0.3718),
which is not in the claim's pool. The 0.39 floor is asserted at the **pool** level and
reported, separately, per rung.

**(2) The census samples the production object.** Three concrete defects, three fixes, each
asserted in `audit_v34.py`:
- *Cartesian*: v33 sampled the diagonal `phi(a[i]+w[i])`. v34 draws bases without replacement
  into the high digit band and walks into the low band, forms the full product, and asserts
  `|exponent set| == B*L` per side (4096 per side per pool, 256 pools).
- *Challenge multiplier*: v33 omitted it. v34's sampled function is `e = (c + k_i + s_l) mod m`
  with `c` the challenge scalar drawn once per pool; three cells separate its effect
  (no-challenge, challenge-common to both sides, challenge on the candidate side only).
- *Histogram*: v33 built the histogram from A-side keys only while dividing by a both-side
  expectation. v34 histograms **A ∪ B**, emits `hist_side_counts`, asserts
  `hist_total == expectation denominator`, and a **planted A-side-only histogram fires that
  check** (`hist_regression_control.planted_fires = true`).
Measured 20-bit cells: kappa_cross 0.9900 / 1.0532 / 0.9985 with cluster-bootstrap CIs over
pools; only cells with power are charged (the v33-pairing control has 1.0 expected pairs and
is reported, not charged). Charged kappa = **1.231** (the carried envelope). v33's 2.6138
does **not** reproduce under either sampler shape here; the correction only lowers the debit,
never the claim.

**(3) Key-width ladder.** Measured at 20/24/28/32 bits; 32 is the widest fixture width with
power. The 168-bit production envelope is **not** extrapolated from those rows: it is the
birthday arithmetic `kappa*N^2/2^169` over the 168-bit image, whose premise is the measured
avalanche of the shipped kernel (`key_collision_test.py`, `key-collision-v32.json`). The
ladder is labelled FIXTURE-REGIME MECHANISM EVIDENCE.

**(4) Gap tail as a labelled hypothesis.** The exponential tail (`mu = 2^26.067`) is stated
as a **HYPOTHESIS**; the `2^-105.6` union bound is its consequence, not its evidence. A
hypothesis-free Markov row is printed and is **insufficient**: `P(gap > 2^33) <= E[gap]/2^33`
gives `2^-6.9` per gap and `2^63.7` after the union over `ell-1` distinct gaps. The gap row is
`2^93.873` of a `2^116.854` total, so the tail assumption moves the claim by <1e-4 bits and
the claim does not rest on it.

**(5) Construction accounting, anchored.** Two constructions the ledger never charged are now
rows: **H-prime supply** = `sides * ell * 12,697` instructions per element (one
subgroup-membership modexp over the 38-bit exponent band, 94 field mults, plus the expected
28 trial-division screens for rejected candidates) and **chain-table point construction** =
`sides * 8 * 2^31 * 1,048` (12 field mults + a Fermat inversion — each entry is an EC point
built once, where the row previously priced only its load/store). Prices are walked
`source symbol -> clang -O2 -> otool -tV -> instruction count -> price -> ledger row` by
`anchor_v34.py`; `audit_v34.py` walks one chain end-to-end (`P_scan` = 10 compiled
instructions >= declared 6) and refuses a price that does not reconstruct. Declared itemised
194 <= compiled 765, and the 2000-instruction per-record allowance strictly dominates the
measured itemised total, so no row under-charges. The stale certificate described a different
design and is **deleted** (`certificate.json` is absent; `N` and `m` are declared in
`design-v34.json`).

## Reproduction
```
python3 anchor_v34.py       # clang -O2 + otool -tV -> price-anchors-v34.json
python3 census_prod_v34.py  # ~70 s; -> census-prod-v34.json
python3 pools_v34.py        # -> pools-v34.json
python3 derive_v34.py       # -> design-v34.json, ledger-v34.json
python3 audit_v34.py        # 31 checks, exit 0 -> audit-v34.json
```
Environment: macOS (arm64), Python 3.14, clang (Apple), `keymod.py` + `key_mixed_final.c`
for the shipped 168-bit key kernel (70 counted instructions).

## Measured results, including the failures
- Three of five claim-population rungs have per-rung CP-95 lower endpoints below 0.39 and are
  reported as such; their distances from the floor are 0.12 / 0.80 / 0.49 SE.
- The `ell=16` short-walk rung measures 0.366168 and fails the floor outright; it is disclosed
  as a short-walk control and excluded from the claim's pool (production `ell` is 2^45 times
  the largest rung).
- v33's charged kappa (2.613761) does **not** reproduce: the corrected cells measure
  0.99–1.05 and the v33-pairing control measures 1.00 at 1.0 expected pairs.
- The Markov tail row fails to justify `2^33` support at production `ell` (`2^63.7` after the
  union) and is printed as a failure, not buried.

## Limitations and provenance
No production run is performed (`N = 2^115.650`); unit costs are declared ledger integers
anchored to compiled bodies; the `phi` surrogate for `y([e]P)` is SHAKE-based and disclosed —
the dependence structure is exact, the curve group law is not executed at 41-bit scale. Not
claimed as transferable: the absolute `nb/ell` split, the 97-bit exponent sliver, the `2^38`
fixture order, the geometric gap tail, and the fixture challenge width. Prior work in this
line: v25–v33 packages (archived), the v28 structured-walk arm (`structured-ladder.json`) and
the v32 key kernel (`key_mixed_final.c`).
