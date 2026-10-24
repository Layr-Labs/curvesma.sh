# argument (v35) -- repair of verdict ea688603 (rejection of gt-instantiable-40pct-v34, claim 127.723)

Claim **127.7536** (reduced 79846/625). W = 5737038977428639211251683104658917243163054, summed in one function from the sixteen rows of
`evidence/ledger35.json`. The verdict named four defects; each is answered in code, and two of the
answers *cost* work (the trap row and the store initialisation) rather than saving it.

## 1. Detection is the code's step, and its power is measured

v34's manifest said "detection is per-worker trail comparison, so an evicted record cannot lose the
walk". The code did not implement that: `run_one` went from the `dict28_wide_lookup` hit straight to
`recover()`. The claim was false for the code, and v35 withdraws it.

The comparison is now a step of the shipped protocol (`scripts/protocol28.py`, `run_one`):

    if trail_ids[ref[0]] == wid:   same_trail_hits += 1        # discarded
    else:                          cross_trail_hits += 1; matches += 1; ... recover()

* It is not decoration. Two records of one trail start from the same (a, b), so
  `G.recover()` has `den = (D2 - q^k D1) = 0` and returns `None`; a same-trail hit is *unusable*.
  The gate discards exactly the hits that cannot produce x.
* It does not lose cross-worker opportunities. One `d28_ws` is created before the worker loop and
  shared by all `nw` workers, so a cross-worker collision *is* a cross-trail hit.
* Detection is therefore store-load-bearing, and the honest question is the verdict's: what is
  P(detect | a true collision)? That is measured, not assumed. `evidence/arm35.json` builds the
  shipped store at the production lambda and **plants the collision**: at arrival t it replays a
  uniformly chosen earlier key through the shipped `dict28_wide_lookup` and counts a detection only
  when the planted record's reference comes back. 215647 of 215858 probes detect,
  kappa_detect = 0.999023 (se 6.73e-05), 3-sigma lower bound 0.998821. Success carries the
  debit explicitly: p_used = 0.578125 x 0.998821 = 0.577443 >= 0.39.

## 2. bounds_table is regenerated from dict28.c, not transcribed

`scripts/bounds35.py` parses the source instead of restating a wish: `PROBES=4`, `ENTRY=9`,
`sizeof(ws_node)=96 B`, `typedef struct { u64 nb; ... }`, `if (nb & (nb - 1)) return -1;`, and
`ws_bucket` masks a **64-bit** `mix64`. A `u64` power of two reaches 2^63, and a 64-bit index cannot
address more, so one level of the shipped type has at most 2^63 buckets. the v34 bucket-count value and its lambda
0.2815 are not representable and are withdrawn (disclosed in `evidence/historical-disclosure.json`).

Production is priced as **two levels of that same type**: 2^63 shards x 2^63 bucket heads = 2^126
heads, capacity 4x2^126, every field `u64` exactly as the file declares. Hence

* lambda = S / 2^126 = 1.126167 (not 0.2815), load factor `8S <= 5*capacity` holds;
* measured retention at that lambda on the shipped four-probe FIFO chains: 0.993527 over 900
  trials, against the Poisson terminal limit 0.994015 (delta -4.88e-04, the fill-phase evictions
  the terminal limit cannot see; the lower of the two is what the success debit uses);
* peak storage 72150101413155552897541746081705523214244 B, and initialisation is **charged**, not deleted: 85070591730234615875067023894796828672 Ir.

## 3. The Fp12 path is counted from the tower and charged as its own row

`scripts/ir25/kernels.c` `f12mul` calls `f6mul` four times, each `f6mul` is six `f2mul`, each
`f2mul` is three `fpmul`: one Fp12 multiply is **72** Fp multiplies. The v34 "54" is withdrawn.
`MUL_CHARGE = 47924 Ir` is the executed `f12mul` span (the 6-limb Montgomery CIOS chain is inside
that span); it is charged as `fp12_dense_mulmod` = executed non-guard steps x 47924, with the guard
branch charged separately at the measured cyclotomic-square span 25372. The coefficient row uses
the same anchor: 2 x ceil(47924/72) + 12 = 1344 Ir/step (v34's 1788 becomes 1344,
a price *cut*).

## 4. The trap charges one event, and that event is quoted

`evidence/trap35.json` holds the three definitions side by side, each with one quoting line:

| where | name | quote |
|---|---|---|
| ledger row | `trap_abort_waste` | `evidence/ledger35.json` rows_ir |
| census counter | `aborts` | `scripts/protocol28.py:353` `aborts += 1; aborted = abandoned = True` under `if armed and stall >= K_STALL:` |
| Markov event | `#aborts <= #guard_conversions` | `evidence/trap31.json` `measured.guard_conversions = 4497` |

`br OR guard` is the **stall** predicate, not the abort predicate: 4497 guard conversions and
1813539 branch events in 5040000 observed steps are the step mix, and v34 charging 4497 as the abort
count was a mis-read of its own record. The bound needs no independence and no p_s^32: an abort
requires an armed counter, `armed` is set only on a guard step and cleared on restart, so distinct
aborts consume distinct arming guards, hence #aborts <= #guard_conversions = 85481963713049917240944295726784289 at S. Each
abort is charged its 32-step wasted window at the ledger's own per-step price (58233 Ir/step)
plus one fresh trail start: 159946670014907108948206726843117711696524 Ir. The raw census aborts observed in the fixture arms are
reported (zero) and are *below* the bound, as they must be.

## 5. What is reported rather than charged

Retention, kappa_detect and the success rate are measured values; they enter the success premise and
are labelled, and no Ir row is derived from them. Every charged term is a Python int (`cert35.py`
re-checks all sixteen rows, the exact claim inequality at n and n-1, and the reduced rational).
