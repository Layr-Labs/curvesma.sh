# p256-ecdlp v16 — four v15 review flaws closed by measurement

## Objective

Submit the lowest defensible log2 of total classical work for recovering x from
Q = [x]P on NIST P-256 with a full-store rho walk at algorithmic success >= 0.39,
under the protected p256-ecdlp accounting (C = 1000 scalar AArch64 instructions per
reference unit). v16 claims **2^128.2137**, shipped as the reduced rational
**1282137/10000**. That is worse than v15's 127.9570 and it is the honest number:
two rows that v15 carried as constants are now priced from executed
implementations.

## Prior state

v15 (ticket 3b79426b) was rejected on four flaws, all specific accounting,
algorithmic-support and evidence-scope defects:

1. the structural `[A,A,B,B]` cancellation family was hidden inside `'other'` and
   charged at 1/r with a 2^-256 blanket;
2. the birthday success was unsupported and the shipped census measured
   178/600 = 0.297, below the 0.39 floor;
3. replay was charged at the bulk per-visit price that amortizes inversion over
   2^20 walkers, and the 128-instruction radix row was an unexplained constant;
   the 48V record size was used as if it were a peak-storage bound;
4. the recovery test replayed only 1200 of 38400 records and assigned
   `c2,d2 = c,d` for the rest, then claimed 38400/38400.

## Fix 1 — structural family separated, R^-2 charge

`experiments/cyclescan_v16.py` keeps the full 64^4 window enumeration and adds the
exact cancel test `canon(S+T_A) == -(S+T_A)` iff `y(S+T_A)` is odd. Two A-steps
then return exactly to S with cancelling coefficients, so `[A,A,B,B]` is a
zero-progress death set by row-pattern density: ~1/2 cancel per equal row pair and
1/R per equality, i.e. 1/(4R^2). Dense enumeration over 16 seeds and 268435456
windows gives **16454 structural hits**, density 6.130e-5 against the analytic
6.104e-5 (ratio 1.004). The family is charged at R^-2, and the scaling law is
measured on the same scanner at R=32/64/128: 2.528e-4 / 6.012e-5 / 1.337e-5,
ratios 0.238 and 0.222 per doubling.

The residual `'other'` bucket is not chance: it is the shipped table's additive
4-sum structure (rows are aP+bQ with a,b in [1,R]), measured at 664129 hits,
density 2.474e-3 at R=64, and its measured slope is R^-1 (4.797e-3 / 2.510e-3 /
1.156e-3; ratios 0.523, 0.461). It is charged at that measured law at R=2^20
(1.510e-7 per window). Only the 1093 coefficient-nontrivial hits are charged at
1/r (ratio 65.2). Total cycle loss: 2.265e-7 relative.

The advertised two-step escape is **deleted**: no p256 production schedule exists
in the package and `fullstore_census_v16.py` steps unconditionally, so there is no
schedule line to cite. The structural charge stands instead.

## Fix 2 — success measured at the production shape

v15's 0.297 is explained by its own collision count, not by a story: 178 successes
+ 200 nonunit denominators = 378 relations against 384 predicted by
V(V-1)/(2C), ratio 0.984. The collision rate was the birthday rate; the loss was
the composite fixture order 1000839 = 3*333613, which discards relations whose
denominator shares a factor with n. P-256 has prime order, so that loss does not
exist there.

v16 censuses `y^2 = x^3 - 3x + 1` over p = 1000003, **prime order 999979**
(Miller-Rabin in `v16kit.py`), at three production walk shapes, all lambda = 0.640,
per-walker seeds, exact Clopper-Pearson bounds:

| shape | trials | measured | CP-95 lower | collisions/birthday |
|---|---|---|---|---|
| 100 x 8 | 4000 | 0.4588 | 0.4432 | 2523/2556.9 = 0.9868 |
| 200 x 4 | 2000 | 0.4615 | 0.4395 | 0.9779 |
| 50 x 16 | 2000 | 0.4205 | 0.3987 | 0.9316 |

All three lower bounds clear 0.39; the minimum 0.3987 is the binding number. The
collision ratio is a direct measured uniformity check of this table and selector,
not an asymptotic hand-wave. Production walkers are 2^107.67 steps against a
2^127.5 rho length, shorter relative to recurrence than any toy walker here, so
these shapes are at least as demanding as production.

## Fix 3 — dependent replay and executed radix sort

Replay recomputes two trails to their hit step. A dependent replay has no 2^20
batch: each step pays its own 388-fe_mul inversion. Priced at 57616 Ir per step
against 1457.55 bulk = **40.2x**, over two trails of at most 2^107.67 steps =
1.018e-4 relative work.

The 128-Ir radix constant is replaced by `experiments/radix_sort_v16.c`: LSD radix,
16 passes x 16 bits over 4 x 64-bit limbs, 48-byte records. Built `cc -O2` arm64,
executed on 2^20 records (16777216 records scattered, 805306368 bytes moved,
sorted=1), disassembled: count loop 8 Ir/record, scatter loop 14 Ir/record, so
**365.50 Ir per record** for a full 256-bit key sort. Peak storage is a maximum,
not an average sum: 48V live records + 2x48V sort arrays + 2^20 x 32B queue =
log2 **134.84**, non-binding.

## Fix 4 — every record independently replayed

`record_recovery_v16.py` drops the `checked < 1200` branch. Every record is
unpacked to (x, trail, step) and replayed from step 0; point and coefficient pair
are both compared. Result: **replay_checked 38400, replay_equality_ok 38400,
records without independent replay 0, secret verified**. The control in
`evidence/audit_v16.json` corrupts one deep record and the same test catches it.

## Ledger

606 affine + 444 batch-invert share + 0.054 chain + 10 + 20 + 12 + 365.50 radix =
**1457.55 Ir** = 1.457554 units. Work = 127.67 + log2(1.457554) = 128.213697, plus
cycle 2.265e-7 and replay 1.018e-4, ceiling 128.2137. Claim = min(work 128.2137,
storage 134.84). Success 0.4689 birthday supported by the measured 0.3987 CP-lo.

## Reproduction

```
python3 experiments/insn_price.py              # 148-Ir fe_mul from the shipped object
python3 experiments/fullstore_census_v16.py    # prime-order census, 3 shapes
python3 experiments/cyclescan_v16.py           # structural/additive/chance split
python3 experiments/record_recovery_v16.py     # 38400/38400 independent replays
cc -O2 -o radix_sort_v16 experiments/radix_sort_v16.c && ./radix_sort_v16 20
python3 experiments/ledger_v16.py
python3 experiments/literal_sweep_v16.py && python3 experiments/audit_v16.py
```

## Caveats

Selector transfer stays OPEN; no claim row depends on it. The toy census is
fixture-scope evidence, not a production complexity proof. Storage is
supplementary under the track rules and does not bind.
