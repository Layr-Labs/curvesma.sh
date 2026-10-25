# v19 -- production cell alone carries the success bound; analytic store; runner shipped

`eab1bb1f` (v18) was rejected on four points. v19 keeps the v18 machinery that
survived (carried 256-bit cursors, planted addressing control, in-package price
anchors, deleted unconditional 2-cycle claim) and answers each point.

## 1. The production walk is now specified and its runner ships (flaw 1)

`experiments/period_ladder_v17.py` is the census runner, **inside the package**.
The walk it runs is fully determined:

* **State.** A visit carries `(S, c, d)` with `S = c + d*x mod q`, i.e. the group
  element and its representation in the basis `(P, Q)`. `x` is the secret.
* **Selector.** `row(t) = (S XOR (t * 2654435761 mod 2^64) XOR (trail * 7919)) mod R`,
  where `t` is the step counter that advances every step and `trail` is the walker's
  index. Row `r` holds `(v, a, b) = ((1+(7r+3) mod R) + (1+(11r+5) mod R)*x mod q, a, b)`,
  so a row is an affine point whose coefficients are known in closed form.
* **Step.** `T = S + v mod q`; `(c, d) += (a, b)`; if `T = 0` the trail restarts
  (identity hit, charged); the **walk steps on the true point `T`** and
  canonicalises only the stored key: `key = min(T, q-T)` with `(c, d)` negated when
  `T` is negated. No walk-state canonicalisation.
* **Initialisation and replay.** Per trial, `x` and `(c, d)` come from
  `Random((20251024<<26) ^ (ti*golden) ^ (R<<3) ^ (W<<50))`, fresh secret per trial.
  A record stores `32B key || 20-bit trail || 108-bit step`, so any hit replays both
  trails from their stored identities: that replay is what decides triviality and,
  when nontrivial, recovers `x = -da/db`. This is why the census is reproducible
  from the record alone.

## 2. The production store is priced analytically, not as a fixture (flaw 2)

v18 shipped a C store whose page table was a bounded fixture (64 special + 256
pages, then NULL, with a linear search) and priced only the loop body. v19 does not
offer that program as the store. `evidence/store_v19.json` prices an analytic
external store and charges every access:

* capacity: `2^127.67` records x 48 B = **2^133.25 B**, page 4096 B -> `2^121.25`
  pages, addressed by a 3-level radix over the page number with fanout 2^40
  (2^41 pointers/node), so translation is one node walk per page, not per record;
  no fixed page count, no linear search, no NULL return;
* charged per record per pass: page-table walk 0.0469 + page open 0.0117 + record store
  1.0 + page flush (4096 B by 32 B vector stores) 1.5000 = **2.5586 Ir**;
* read-back materialisation (the objection that `main`'s per-pass `ext_read` was
  not counted): 0.0469 + 1.5000 + 1.0 = **2.5469 Ir** per record per pass, charged on all
  16 passes although only the count and witness passes read the stream back;
* the carried-scatter body is quoted at **one** number everywhere, 36.0 Ir/record
  (72 insns / 2-record unroll). v18 quoted 24.5 in the manifest while the review was
  granted 36.0; that inconsistency is gone.

Per pass 0.25 + 22.0 + 2.062 + 36.0 + 2.5586 + 2.5469 = **65.4175 Ir/record**; 16 passes +
37.5 witness = **1084.2 Ir/record** (v18: 818.5). The price went **up**: the claim rises
from 128.7369 to 128.7921 as the cost of the accounting, and that is stated, not
hidden. `experiments/radix_sort_v18.c` is retained as a **unit test of the 256-bit
addressing arithmetic only**: planted cursors 2^32+17 / 2^64+5 / 2^128+3 place all
4096 records exactly (`bad_new=0`), each read back from its own address
(`distinct_locations=1 readable_at_5=1 no_alias_at_2p64=1`), while the v17 low-word
path writes the same record at position 5 and faults at a 2^40 low word.

## 3. The success bound rides on the production cell alone (flaw 3)

The claimed object is the R = 2^20 cell **by itself**, extended from 100 to **200
trials** (same q, same seed rule; the first 100 are the v18 trials and are reproduced
exactly, which `evidence/ladder_R2p20_W1_16p.json` records as a bit vector and a
prefix check):

```
95/200 = 0.4750 verified successes, exact Clopper-Pearson lower 0.404120 >= 0.39
145 collisions, useful_per_collision 1.00, 16 full selector periods, fresh secret per trial
```

The bound as trials land: 100 -> 0.388644, 120 -> 0.399259, 140 -> 0.414390,
160 -> 0.420040, 180 -> 0.424712, 200 -> 0.404120. No cross-R, cross-order or
cross-walker-count pooling is used to carry the bound. The smaller rungs are
reported as **supporting** measurements and the pooled number is labelled
descriptive:

| cell (supporting) | R | periods | trials | wins | measured | CP-95 lower |
|---|---|---|---|---|---|---|
| R256_W64_16p | 256 | 16 | 400 | 163/400 | 0.4075 | 0.3589 |
| R1024_W64_16p | 1024 | 16 | 300 | 146/300 | 0.4867 | 0.4288 |
| R4096_W64_16p | 4096 | 16 | 240 | 107/240 | 0.4458 | 0.3819 |
| R64_W64_16p | 64 | 16 | 600 | 287/600 | 0.4783 | 0.4377 |

What is **assumed**, stated as an assumption because it is the one link not
measured: the per-period row hazard is the same in both regimes. Its justification
is structural (the row index is the state representative XOR an advancing counter,
and the representative is a uniform element of a prime-order cyclic group in both
regimes) and its **direction is conservative**: production affords 2^87.67 selector
periods per walker, not 16, and the within-trail collision probability is monotone
non-decreasing in the number of periods, so the measured 16-period rate is a lower
bound on the production per-walker probability. v19 does **not** claim a 49-bit
fixture is isomorphic to the P-256 subgroup.

## 4. The gap-2 congruence had the wrong signs (flaw 4)

The review is right. From `eps2*(eps1*(S+v1)+v2) = S` with `eps1, eps2` in
`{+1,-1}` and `1/eps2 = eps2`:

```
eps2 = eps1:   v2 = -eps1*v1
eps2 = -eps1:  v2 = -2*eps1*S - eps1*v1        (all mod q)
```

v18 printed the opposite signs. v19 states the congruence above, derives it in one
line, and control (c) in `audit_v19.json` checks all four `(eps1, eps2)` cases:
every corrected value returns to `S`, and **no** v18-signed value does. The row is
still not load-bearing: the cycle charge is now a measurement on the declared walk
(section 5).

## 5. The cycle row is measured on the declared walk

`evidence/cycle_declared_v19.json`: R = 64, W = 64, 16 periods, 600 trails, exact
true-state tracking on (the declared walk, not the v16 canonicalised-state walk):
**0 deaths in 600 trails**, 95% Clopper-Pearson upper bound on the per-trail death
rate 0.006129, and a death kills at most one of 2^20 walkers, so the charge is
**5.845e-09** relative (v18 charged 2.265e-07 from a different, canonicalised walk).
Positive control: the v16 walk cycles on every trail in the same harness, so the
zero is a property of the declared walk, not of the tracker.

## 6. Claim, machine-copied from the ledger

`experiments/ledger_v19.py` computes the bound once and writes
`evidence/claim_v19.json`; this script copies the reduced rational and asserts the
score is >= the computed log2 work.

Rows: 606 affine + 444 batch-invert + 0.05434 chain + 10 + 20 + 12 + 1084.2 sort and
witness = **2176.2 Ir** = 2.176234 units at protected C=1000. Work =
127.67 + log2(units) = 128.791834, plus cycle 5.845e-09 and replay 1.755e-04 relative, plus table
construction. **Claim = 128.7921 = 1287921/10000.** Storage log2 134.2550 bytes (96 V, external,
disclosed, unscored).

## Audits (`evidence/audit_v19.json`, all executed, all pass)

* (a) addressing unit test: 2^32 / 2^64 / 2^128 planted cursors exact; v17 path
  misplaced at 5; retrievable from each address.
* (b) success from the shipped outcome vector: 200 bits -> 95 wins -> 0.404120, and the
  first 100 bits -> 49 wins -> 0.388644, i.e. the prefix is intact and the new trials
  move the bound.
* (c) gap-2 signs: corrected values return, v18-signed values do not.
* (d) witness: collision-free fill emits zero witnesses.
* (e) store price: ledger row equals the analytic store total (1084.1795) and exceeds
  the body-only price (60.3120), so the store is charged, not skipped.
* (f) cycle: 0 deaths on the declared walk; v16 positive control cycles.

## Limitations

* The census fixture is a prime-order cyclic group with a representative-based
  selector, not the P-256 x-coordinate selector. Section 3 states the single
  assumption that bridges them and its conservative direction; it is not measured.
* The production cell runs W = 1 walker per trial at 16 periods. The 2^20-walker,
  2^87.67-period configuration is not executed -- it is the regime the measured
  hazard is charged against.
* Store capacity and transfer volume are disclosed at log2 134.2550 bytes and charged
  as instructions, not as wall-clock I/O.
