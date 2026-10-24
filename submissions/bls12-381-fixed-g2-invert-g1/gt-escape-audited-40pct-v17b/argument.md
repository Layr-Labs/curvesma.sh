# Argument — gt-escape-audited-40pct-v17b

Claim: `time_log2 = 125.732`, success at least 0.39, every-run work at most
`70648576363634561378687894185206506161` units (`2^125.731..125.732`).

## 1. Reduction

`Q` nonidentity in G2 makes `g = e(P0,Q)` a generator of GT's order-`r` subgroup, so
`z = g^x` has a unique `x` and `[x]P0` satisfies `e(P,Q) = z`. The walk carries the
label `(A,D)` with the invariant `y = g^A z^D`, so a class collision with distinct
labels yields a linear equation in `x` modulo `r`.

## 2. Disjoint-half slots (v17 repair, propagated here)

`T = Re(c0)` is a free class invariant, identical on all twelve Frobenius conjugates.
The production map is `prod_bits(t, 4096, 100)`:

    br = (t & 7) < 3                 bits 0..2
    h  = (t >> 3) & 2047             bits 3..13
    dp = ((t >> 14) & (2^100-1)) == 0   bits 14..113
    slot = h + 2048*[kappa >= 6]     kappa = canonicalisation index

The three ranges are disjoint and the half select is a class function, not a T bit.
`p^6 = -1 mod r` is why the halves exist: in a `p^6` pair, `kappa` differs by 6, so
the two steps read different halves and the inverse-multiplier cancellation that
killed v16's walk cannot occur. No escape rule is used anywhere.

The v16 defect was that `gt_prod.py` kept `idx_bits=10` and started the DP window at
bit 13, so `h` and the DP predicate shared a bit and every DP hit forced `h` into the
lower half. `gt_prod17b.py` passes R=4096 explicitly, asserts the effective map and
publishes it as `config`; evidence is regenerated at exactly that configuration
(120160 invariants vs two independent derivations, 0 mismatches, witnesses for
every index in both halves with a simultaneous DP hit).

## 3. Degeneracy

R=4096 slots are two disjoint halves. Residual risk would be two distinct slots that
are exact Frobenius multiples; the 49152 scaled label pairs
`(p^k a_u, p^k b_u) mod r` are pairwise distinct (0 collisions), so none exist.
Toy measurement of the repair: degenerate matches fell from 5816/400 trials under the
v16 rule to 12/1600 under this one, and every match converted to a verified `g^x = z`.

## 4. Accounting

`W = H + R/C` with C = 20000. Rows: common 301; multiply branch 2505 (including the
new 12-instruction half-select row); squaring branch 15042+40; setup
`4096 x (510 C + 12 x 2400) + 2^20 = 41898213376`; auxiliary rows solve
5340504 instructions for solve+verify+output. Nominal work
`2^125.731..2^125.732`, inside the cap `2^125.731..125.732`.

The meter charges every routine and checks the cap BEFORE the step runs. The executed
run's event log, its meter and the row formula agree as exact integers
(41921914768 each), and the cap-boundary run stops at exactly the admitted count.

## 5. Success floor

The tail `1 - exp(-S(S-1)/(2N/rho))` with `rho = (1-9/64)(1-1/4096)` at
S = 71192122588783531220144472073284087098 gives 0.392418. The design-budget census
(80000 trials, 30854 usable collisions, 0.385675) calibrates the model at
the toy's own design budget (0.384972), giving a one-sided 95% lower ratio of
994464441/1000000000. Discounting the tail and subtracting the abort, DP-detection and
storage debits leaves 0.390100 >= 0.39. The discount is why S is
0.00542 bits above v17's budget. Claiming 0.39 from the undiscounted tail was the v16
error; claiming less than 0.39 is not permitted by COST_MODEL.md line 45.

## 6. Storage and limits

Records 57209259 x 52 B = 2974881468 B, R=4096 table 28311552 B (12 keys x 576 B),
code/scratch/sort rows itemised, peak 3055621820 B < 2^32, DP waits charged as a
256-multiplication row. Overflow refuses new records and pays a Chernoff debit.

H1'', H2 and the random-mapping idealisation remain disclosed heuristics; the census is
a calibration at toy scale, not a production-size execution.
