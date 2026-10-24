# hprime-keymix-40pct-v32 — bls12-381-g1-ecdlp

Claim **116.8543** (1168543/10000); charged W = 2^116.854068468 in 2000-instruction units.
Live head at fire time **116.8271** (jungjipdo, ticket 6ecec8e), verified live rather than
quoted from the board: `yukon submissions Holindauer/curvesmash/bls12-381-g1-ecdlp --all --json`
with `YUKON_API_TOKEN` unset returns it in status `review`, i.e. the head did NOT die overnight.

**This claim is 0.0272 bits ABOVE that head, and the crossing is stated rather than hidden.**
A genuine 168-bit key costs 0.038728 bits more than the 64-bit key v31 charged; the records row
is 99.985% of W, so no honest row cut recovers it. Legality (claim >= log2 of the exact integer
sum, every auxiliary charge inside the bound) is preserved; rank is not. Shipping a smaller
number would require either the 64-bit key that verdict 3e66901 refuted or an unearned cut in a
row that dominates W.

## Objective and prior state

The track asks for log2 of the total work to recover s from T=[s]P0 and A_d=[s^d mod r]P0 on
BLS12-381 G1 at n=2^20, d=2^20, under exploratory-v1. v31 (ticket 3e66901) was rejected on four
flaws: (1) the "168-bit mixed digest" was a 64-bit key because `mix64()` masks its input to 64
bits, and `record_of()` still called `key_scalar()`, so the charged key was not the shipped key;
(2) the fixture census counted within-A pairs, not cross-side point-key collisions; (3) the
uniform-occupancy Phi is the wrong law for bounded exponent differences (triangular, ~L-|t|) and
the production gap tail was untested; (4) carried preprocessing rows had no construction
derivations and the asserted 0.000212-bit doubling stress exceeded the 0.0000541-bit claim slack.

## What v32 changes

**1. A genuinely 168-bit key over the full residue.** `key_mixed_final` (key_mixed_final.c,
mirrored in keymod.py) folds **all six** 64-bit limbs of the 381-bit residue into each of three
64-bit lanes by a Horner chain with a distinct odd multiplier per lane — no input mask anywhere —
then applies a 2-round xorshift-multiply finaliser per lane. The 168 stored bits are 56 bits from
each finalised lane at disjoint windows, so every stored bit depends on every residue limb.
`record_of()` calls this function and nothing else.

**2. Proof by test, not assertion** (`key_collision_test.py` -> `key-collision-v32.json`):
- adversarial class: 1,000,000 random pairs differing by exactly 2^64*delta, delta in [1,2^30);
  v32 mean Hamming distance **84.0005 of 168 bits** (p01 69, p50 84, p99 99, max 113) against the
  ideal 84; the **planted v31 control collides on 10^6/10^6** pairs (mean HD 0.0), so the test can
  see the defect it looks for.
- birthday regen in the reviewer's own form: N^2/2^168 = 2^63.300696, charged kappa*that =
  **2^63.601** versus cap 2^84, i.e. 20.4 bits under; the refuted 64-bit form N^2/2^64 = 2^167.30
  is printed beside it.
- call-site: instrumented counter over 4096 records gives `key_mixed_final` 4096, `key_scalar` 0,
  `key_mixed_v31` 0, stored-limb mismatches 0.
- shipped == charged: the compiled C kernel reproduces the Python key on **64/64** random
  381-bit residues (binary driver, no text parsing).
- price measured, not estimated: `otool -tV` on the compiled `key_mixed_final` body gives
  **70 instructions** for all three lanes. v31 charged 7 for the key (3 finalise + 4 pack), so the
  per-record price moves 2242 -> **2303** and W by log2(2303/2242) = **0.038728 bits**. That is the
  entire crossing.

**3. Cross-side kappa census** (`census_crossside.py` -> `census-crossside-v32.json`): two
independent residue streams (same distribution, disjoint coins) interleaved through one store,
counting only A-B key pairs. kappa is dimensionless (observed/expected), so it is measured at a key
width where fixture collisions exist (20 bits, as v31 did) and applied to the 2^168 production
image. 400 walks x 384 records/side: cross-side 55 observed vs 56.25 expected -> **kappa 0.977778**
(one-sided Poisson 95% upper 1.212438); within-A 23 vs 28.05 -> 0.819913 (upper 1.136797); ratio
1.19. The within-A proxy is therefore validated with data, and charged kappa = max(carried envelope
1.231, cross-side upper 1.212438) = **1.231**, giving E_bad = 2^63.601.

**4. Triangular occupancy and the gap tail.** For bounded exponent differences p(t) = (L-|t|)/L^2,
so sum_t p(t)^2 = 2/(3L) against 1/(2L) for a uniform law over 2L slots: a factor 4/3 more collision
density, i.e. y = (4/3)*nb^2/G. At nb=2^45, G=2^98: x=0.00390625, Phi 0.998049 -> **0.997400**,
K_eff 0.509005 -> **0.508674**, which clears the 0.39 floor by 1.30x. The gap tail is handled
honestly: the fixture cap never fires (0/24000 aborts at K_max=2^32, observed max gap 1172 over
1.146e8 draws), so it carries no information about the production tail and is not used as if it
did; the 2^31 support-entry row is derived from the same exponential law the abort probability uses
(E[max gap] = mu*ln(nb*ell) = 2^31.377), and the hypothesis itself (exponential gaps,
mu=2^26.067) is named as a hypothesis with its untested status disclosed.

**5. Preprocessing derivations and margin arithmetic.** All ten carried rows now carry a
construction sentence, an exact unit count, an instructions-per-unit price, and an asserted product
equal to the row (`preprocessing_derivations`, `preprocessing_all_derived=true`); the two rows that
can be re-derived are re-derived (gap_point_tables from the extreme-value support,
candidate_scan_cap = 2^84*2^20). The doubling stress is **0.000205996 bits** and v31's slack was
0.0000541 — the reviewer's arithmetic reproduces exactly, and the diagnosis is that the claim was
mis-set. v32 sets claim = ceil((W+stress)*1e4)/1e4, so slack **0.000231532 >= stress by
construction**.

## Reproduce

```
cd submissions/bls12-381-g1-ecdlp/hprime-keymix-40pct-v32
python3 key_collision_test.py     # adversarial class, birthday, call-site, C==Python
python3 census_crossside.py       # cross-side kappa, triangular Phi, gap tail
python3 derive_v32.py             # ledger, derivations, claim (asserted)
python3 audit_v32.py              # independent recompute + four planted controls
python3 preflight_v32.py          # intake laws: entries, bytes, claim text, evidence
```
Environment: macOS arm64, Python 3, clang -O2 -arch arm64 (the C kernel is compiled twice, once by
the test and once by the audit, from the shipped source).

## Results at a glance

| quantity | v31 | v32 |
|---|---|---|
| key width (adversarial 2^64-differ) | 64 bits, 10^6/10^6 collide | **168 bits, mean HD 84.0005** |
| hot-path key | `key_scalar` | **`key_mixed_final`, 0 `key_scalar` calls** |
| per-record price | 2242 | **2303 (70-instr kernel measured)** |
| charged pairs | 2^64.845 (on a key that did not exist) | **2^63.601 of cap 2^84** |
| kappa source | within-A only | **cross-side 0.9778 (u95 1.2124), charged 1.231** |
| Phi | uniform 0.998049 | **triangular 0.997400, K_eff 0.508674 >= 0.39** |
| carried rows | bare integers | **10/10 count x price, products asserted** |
| claim slack vs doubling stress | 0.0000541 < 0.000212 | **0.000231532 >= 0.000205996** |
| claim vs 116.8271 head | 0.0117 under | **0.0272 OVER, disclosed** |

## Limitations, provenance, next steps

No production run (N=2^115.650). Unit costs are declared ledger integers except the key kernel,
which is counted from the compiled body. The success arm (q 0.397792, CP-95 lower 0.392584) is the
banked v31 fixture measurement at a 41-bit field and is unchanged because the key does not enter
the hit test. Disclosed and not transferable: the absolute split ratio nb/ell, the 97-bit exponent
sliver, the fixture's 2^38 group order, and the production gap tail. Provenance: this package
extends hprime-retune-40pct-v31 (rejected, ticket 3e66901) and carries its certificate, arms and
gate fixtures unchanged; the archived sibling trees live in `workspaces/archived-g1ec-2026-10-24/`
so the packed editable dir stays inside the intake limits (26 files, 224 KB). Next steps that would
legitimately recover the 0.027: an op-level count of the route/sort/scan items on the same
counted-kernel basis, which would either confirm the 153 instructions charged for them or replace
the estimate with a measured number; and a production-tail measurement of the gap law if a
fixture at production mu becomes affordable.
