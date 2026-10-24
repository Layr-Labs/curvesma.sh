# g2-invert v35 -- detection in the code and measured, bounds from source, trap predicate quoted

Claim **127.7536** (reduced 79846/625). This is a repair of v34 (ticket `ea688603`, claim 127.723),
which was rejected on four concrete flaws. Each is answered below in the order the verdict raised
them, and two answers *raise* the work rather than lower it.

## Objective

For nonidentity Q in BLS12-381 G2 and uniform z in GT, compute g = e(P0,Q) and recover x with
z = g^x by a Frobenius-class rho walk in GT on 1024 workers, returning [x]P0 verified by e(P,Q) = z.
The claim is the smallest n/10000 with 20000*2^(n/10000) >= W, where W is the exact integer sum of
sixteen rows in `evidence/ledger35.json`.

## Environment and setup

macOS arm64, clang -O2. The dictionary is the shipped `scripts/ir28/dict28.c`, compiled to a dylib at
import time by `protocol28.py` (no compiled object ships: objects are not text and intake rejects
them). The tower and field constants are `scripts/ir25/kernels.c`. Run:

    cd submissions/bls12-381-fixed-g2-invert-g1/gt-instantiable-40pct-v35/scripts
    python3 bounds35.py      # re-parse dict28.c, regenerate bounds_table, sweep + planted control
    python3 arm35.py         # retention and kappa_detect at the production lambda (900 trials)
    python3 smoke35.py       # run the patched protocol; asserts the trail gate both ways
    python3 regen35.py       # rebuild ledger35.json, submission.json, argument.md
    python3 cert35.py        # independent re-derivation of all sixteen rows + four planted controls

`cert35.py` exits 0 with `ok: true` and all four controls true.

## Prior work / baseline

v34 (claim 127.723) priced a store table it did not implement, claimed a detection mechanism the code
did not contain, derived an abort charge from a counter that measures something else, and divided the
Fp12 cost by 54. Earlier revisions were rejected on their own grounds; nothing from them is re-used
here except the censused per-symbol prices in `evidence/exec31.json`, `evidence/sort32.json` and the
step budget S = 9.580366847e37.

## 1. Collision detection: the comparison is now a step, and its power is measured

v34's manifest said collisions were found by per-worker trail comparison while the code went from the
store hit straight to `recover()`. That was false and it is withdrawn. `protocol28.py` now performs
the comparison: a lookup hit counts as a collision **only** when the stored trail id differs from the
arriving worker's. Same-trail hits are counted and discarded, because two records of one trail share
`(a,b)`, so `recover()` divides by `D2 - q^k D1 = 0` and returns `None`.

This does not lose cross-worker collisions: one `d28_ws` is created before the worker loop and shared
by all 1024 workers, so a cross-worker collision is a cross-trail hit. Detection is therefore
store-load-bearing -- and it is measured rather than assumed. `evidence/arm35.json` runs the shipped
four-probe FIFO store at the production load factor and **plants the collision**: it replays a
uniformly chosen earlier key through the shipped `dict28_wide_lookup` and counts a detection only when
the planted record's reference comes back (a wrong-reference hit is not counted). 215647 of 215858
probes detect: `kappa_detect = 0.999023` (se 6.726e-05). Success carries the debit:
`0.578125 x 0.998821 = 0.577443`, above the 0.39 floor.

## 2. Store bounds regenerated from `dict28.c`

`scripts/bounds35.py` parses the shipped source: `PROBES=4`, `ENTRY=9`, `sizeof(ws_node)=96 B`,
`u64 nb`, `if (nb & (nb-1)) return -1;`, and `ws_bucket` masking a 64-bit mix. One level of that type
addresses at most 2^63 buckets, so v34's declared head count is not representable (the disclosure file
`evidence/historical-disclosure.json` records what was withdrawn and why, and it is the only place
those values may appear).

Production is priced as **two levels of the same type**: 2^63 shards x 2^63 heads = 2^126 heads,
capacity 4x2^126, every field `u64` exactly as the file declares. Hence `lambda = S/2^126 = 1.1262`
(not 0.2815), with measured retention **0.993527** over 900 trials against the Poisson terminal limit
0.994015 -- the lower figure is the one the success debit uses. Peak storage is the sum of the
regenerated rows (7.22e40 B) and initialisation is **charged** (8.51e37 Ir), not deleted. An
unsuccessful control is reported too: the arm at the withdrawn v34 lambda cannot be run at all,
because no bucket count of that size can be initialised by the shipped field.

## 3. The Fp12 path is counted from the tower

`kernels.c` `f12mul` calls `f6mul` four times; each `f6mul` is six `f2mul`; each `f2mul` is three
`fpmul`. One Fp12 multiply is **72** Fp multiplies, so v34's "54" is withdrawn. The measured `f12mul`
span (47924 Ir) is charged as its own row (`fp12_dense_mulmod`) against the executed non-guard step
count, the guard branch is charged at the measured cyclotomic-square span (25372 Ir), and the
coefficient row becomes `2 x ceil(47924/72) + 12 = 1344` Ir/step. v34 charged 1788: this is a price
*cut*, taken because the divisor was wrong.

## 4. One trap event, three identical definitions, one quote each

| where | name | quote |
|---|---|---|
| ledger row | `trap_abort_waste` | `evidence/ledger35.json` rows_ir |
| census counter | `aborts` | `scripts/protocol28.py:351`: `aborts += 1; aborted = abandoned = True` under `if armed and stall >= K_STALL:` |
| Markov event | `#aborts <= #guard_conversions` | `evidence/trap31.json`: `guard_conversions = 4497` |

`br OR guard` is the *stall* predicate, not the abort predicate; 4497 guard conversions and 1813539
branch events over 5040000 steps are the step mix and are not charged as aborts. The bound needs no
independence and no p_s^32: an abort requires an armed counter, `armed` is set only on a guard step
and cleared on restart, so distinct aborts consume distinct arming guards. At S that bounds 8.548e34
aborts, each charged its 32-step wasted window at the ledger's own per-step price (58233 Ir) plus one
fresh trail start (7660 Ir): 1.599e41 Ir. The raw census abort count is zero, which is below the
bound as it must be.

## Measured results

| quantity | value | source |
|---|---|---|
| W | 5.737038977e42 (log2 142.041286) | ledger35.json, sixteen rows |
| claim | 127.7536 = 1277536/10000, minimal n verified at n and n-1 | cert35.py |
| kappa_detect | 0.999023 (se 6.7e-05, 3-sigma lower 0.998821) | arm35.json, 215858 planted probes |
| retention at lambda 1.1262 | 0.993527 (900 trials, se 3.9e-05) | arm35.json |
| p_used | 0.577443 >= 0.39 | ledger35 success block |
| gate | same-trail hit rejected, cross-trail hit accepted, shipped lookup returns the reference | protocol35.json |

## Limitations

The census runs 8/64/128 workers, not 1024; only the partition series at fixed S carries the
W-dependence, and that is stated rather than papered over. `kappa_detect` is measured at fixture
scale at the production load factor, not at 1024 workers. The trap bound is a counting bound: it is
valid and it is loose. The store's two-level instantiation is a specification whose per-level fields
are exactly the shipped types; no 2^126-head array is allocated anywhere here.

## Provenance

Everything is inside this directory: `scripts/ir28/dict28.c` and `scripts/ir25/kernels.c` are the
shipped arithmetic, `evidence/exec31.json` and `evidence/sort32.json` are the executed price anchors,
`evidence/trap31.json` is the census that supplies the 4497 count, and
`evidence/historical-disclosure.json` is the only place withdrawn v34 values are named. No file
outside this tree is referenced by path.
