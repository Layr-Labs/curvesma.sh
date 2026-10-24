# gt-executed-40pct-v27 — argument

This package answers the five flaws of verdict `55e7de49` (v26). Every number is
emitted by `scripts/gen27.py` from `evidence/*.json`; the certificate fails if any
of them disagrees with the claim. v26 is left in the tree as published provenance.

## 1. The manifest is generated, not asserted

v26's claim text said the `Rg = 65536` survival arm was measured "at both fixture
scales"; `survival25.json` holds **one** scale, `rp = 521041`. v27 does not re-word
the sentence. `evidence/manifest27.json` is produced by walking the package: every
file, its byte count and its `sha256`. The survival row names the single measured
scale, `not_measured_at_rp10566001 = true` stays, and the certificate asserts that
every manifest path exists.

## 2. Success is a bootstrap of the product, not a product of quantiles

v26 charged `CP-95 x point retention / band` — a product of quantiles presented as
a quantile of the product. v27 draws the three factors:

    census A    ~ Beta(292+1, 512-292+1)   (292/512 at rp=521041)
    census B    ~ Beta(304+1, 512-304+1)   (304/512 at rp=10566001)
    retention   ~ Beta(398+1, 400-398+1)   (398/400 pairs at production Rg)
    divisor     ~ U[1, 1.358246]           (fit residual envelope, both ends)

10 000 draws, seed 20261024, product = min(census A, census B) x retention /
divisor. Result: **5th percentile 0.41336**, 1st percentile 0.39805, median
0.47582. The bound is the 5th percentile and it is above the 0.39 floor, so S needs
no resizing. The eviction effect on success stays **unmodelled** and is stated as
such: the shipped dictionary evicts (306 evictions in 20 000 ops), and a productive
hit lost to eviction is a success loss this bound does not charge.

## 3. One address layout, one memory formula

`records = S = 9.5804e37` visited classes, one record each. A store that keeps
occupancy `<= 5/8` cannot hold S records in S slots:

    slots       = ceil(S / 0.625) = 153285869553487776880359083998396376144
    slot ref    = ceil(log2 slots) = 127 bits -> 2 words
    slot layout = 6 key words + 2 ref words + 1 tag = 9 words = 72 B
    bytes       = S x 688 (pool) + slots x 72 (index) + 1024 x 256 (worker state)
                  + 1024 x 6 x 8 (ring) + 2^20 (code)

The v26 pair — 8 Ir per record on the index row while `gen26.py` counted
`65536*7*8` table bytes and the prose said `S*64` — is withdrawn. There is one
formula, machine-derived, and the index row now charges 9 Ir per visited record
(one touched slot, 1 Ir per word). The byte count is reported for consistency, not
scored: the bound is a meter bound on instructions, and a 2^127-slot store is not
claimed buildable. The fix removes the contradiction and makes the row dearer, not
cheaper.

## 4. Trap: ring 6, and a census instead of the rule of three

The guard ring keeps the last **6** multiply slots (v26 scanned 2-slot words; its
guard held 5). A trap needs some slot whose accumulated kappa sum vanishes in
Z[zeta], and vanishing needs that slot at least twice, so under ring 6 **no** slot
family of length <= 6 can trap: the exact census over canonical words of length
2..6 with no repeat within distance 6 returns zero surviving families
(`ledger27.json -> trap.family_census`; W(1) = 0 is the whole argument). The
period-6 family `abcdefabcdef` is blocked for the same reason: its second `a` sits
at distance exactly 6, inside the ring.

The first family the ring permits repeats a slot at distance 7, i.e. a window of
length >= 14. So the charge comes neither from the census (empty) nor from the rule
of three on correlated steps. It comes from the **measured protocol traffic**:
7 680 000 steps, 0 admissions, 0 aborts. With correlated steps the independent-trial
count is the number of ring windows, not the number of steps, so the rate is
`3/(7680000/6)` per window, and the charge is

    (3/(7680000/6) x 32 + (1-f_mul)^32 x 32) x S x STEP = 168 691 606 152 ... Ir

Period-6 is not claimed excluded by the census — it is excluded by the ring, and any
family beyond ring 6 is covered only by this measured charge. The forced-admission
build stays labelled a mechanism existence proof.

## 5. Executed counts: the priced spans are loop-free

v26's objection was that `otool` spans count static instructions, "not executed
instructions for looping routines". `evidence/loopfree27.json` tests the direction
the claim needs: inside each priced span, **no branch targets an address at or below
itself** — 0 of 25 372 instructions in `cyclotomic_sqr`, 0 of 47 924 in `f12mul`.
Control only moves forward, so every instruction executes at most once per call and
the static span is an **upper bound on the executed count**, which is what covering
W requires. It is not claimed tight: a forward branch can skip an arm and the true
executed count is smaller.

Also removed: v26 charged `SQR_CHARGE = 27611` against a measured span of 25 372,
an unexplained +2 239 Ir on every squaring. v27 charges the span. `STEP` falls from
24 284.746 to **23 477.404** Ir at the same measured `f_mul = 0.639418`.

The 12 x 64-bit key question: `state_key48` reads the twelve 64-bit words that
`Group.state` sorts and folds them into a 48-byte key. Production fields are
381-bit, so the key is a **hash of the state, not the state**: a key collision is
not a field collision, and the six-word compare in `dict_lookup` is what decides a
hit. The accounting is the differential's counters (8 458 hits, 253 misses, 306
evictions, 0 structure mismatches) plus the collision row
`ceil(S(S-1)/(2N))` replays, which charges every match at full verify cost.
Widening to 7 limbs (448 bits) would add 7 Ir per lookup (+0.00035 units/step); it
is not done, because the shipped key is a 48-byte fold and every row is priced on
the shipped symbol.

## Ledger (exact integers)

| row | Ir | derivation |
|---|---|---|
| walk_steps | 2 249 221 408 965 905 297 733 409 839 103 413 355 600 159 | S x STEP |
| per_step_lookup_key_and_probe | 61 889 169 832 220 689 915 444 980 164 352 536 868 140 | S x 646 |
| record_writes_and_insert | 89 480 626 351 848 489 753 909 615 284 063 884 574 060 | S x 934 |
| collision_processing | 8 192 | 2 matches x 4096 |
| trap_waste_bound | 168 691 606 152 510 875 773 251 867 992 491 670 823 | measured ring windows x executed delay |
| index_slots_widened_refs | 862 233 016 238 368 744 952 019 847 490 979 615 810 | S x 9 (layout row) |
| worker_starts_1024 | 7 843 840 000 | 1024 x 383 |
| code_and_advice | 1 048 576 | 2^20 |
| candidate_verification | 57 526 800 | 1200 x f12mul + 18000 |

W = 2 401 622 129 772 365 357 023 489 706 267 321 150 752 560 (2^140.7850).
Claim = smallest n/10000 with 20000 x 2^(n/10000) >= W, exact integer check with
n-1 asserted to fail: **126.4973**. Against the track frontier 125.72 this is
0.78 bits above; the package does not claim to beat its own promoted ticket. It
claims to be a valid bound with the five v26 flaws repaired.

## Cost model

    STEP = 0.639418 x 22409 + 0.360582 x 25372 = 23477.404 Ir
    S = 829/572 x sqrt(N), N = (r-1)/12, 1024 workers, one record per visited class

## Reproduce

```
python3 scripts/loopfree27.py
python3 scripts/verify27.py
python3 scripts/gen27.py
```

## Withdrawn from v26

The "measured at both fixture scales" survival sentence; the CP x point / band
success composition; the 8-Ir index row and the `65536*7*8` table count; the
2-slot trap words and the rule-of-three waste rate; the +2239 Ir squaring
surcharge.
