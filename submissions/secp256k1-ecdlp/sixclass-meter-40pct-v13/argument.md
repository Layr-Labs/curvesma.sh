# Six-class metered rho v13: conditional 127.03

## Scope

The claim is success >= 0.40 within 2^(12703/100) target-operation equivalents,
conditional on the explicit laws below, not a production execution. Use the
protected secp256k1 p, r; N=(r-1)/6 = 2^253.41504 classes; C=1000
(tracks/docs/COST_MODEL.md); R=2^20 twisted table entries; DP spacing 2^69;
trail cutoff L=40*2^69; 2^20 workers; 2^16 storage nodes with 2^43-record caps;
192-byte records; S=2^85 events. Every arrival, re-walk, escape, emission and
verification is charged once; the report pool aborts before exceeding a cap, so
no record is evicted or filtered.

## Step, key and escapes

```
STEP(state):
  load 192-bit key K; branch = K bits 0..19, twist = 3 bits at 64..66
  Q = Q + T[branch][twist]            # one reference affine addition, 6M
  f0 = top 64 bits of x(Q)
  f1 = 64-bit fraction of beta*x(Q)    # 8 limb products, column 5/4
  f2 = 64-bit fraction of -x(Q)-beta*x(Q)
  if two of f0,f1,f2 are within 40 ulp: exact 6-orbit argmin (metered)
  else: argmin on the three fractions
  force even y by flipping the table twist (no field operation)
  store K; compare K with K_{-1} (2-cycle) and with the checkpoint key
  if either matches: double the minimum-key member, clear history, retest DP
  if K is a distinguished point: emit full-start record to node K bits 128..143
  tick the global meter; halt at S_ops
```

Only doubling escapes are used. A fixed additive escape is rejected because it
creates structural composite cycles through negation cancellation; the toy
control in section "Toy evidence" still traps at that design.

## Ledger (log2 units)

- S_ops = ceil(sqrt(2*N'*ln(1/0.58))) = 2^126.76932 arrivals, N'=N*R/(R-1).
  This is the hard phase-1 operation counter, not a sampling estimate.
- Phase 1 = 1.11*(S_ops + 2^60) = 2^126.91988. The 2^60 term pays the start
  and meter-block remainders; every one of the S_ops arrivals is counted.
- Hard extras 2^89.00: table 2^26, 2^16 node setup at 2^20 each, meter
  remainders 2^61, 2^68 report traffic, 2^89 index fallback allowance
  (rate <= 2^-50 at 2^12 units), 2^40 verification, 2^24 starts.
- Random extras, expected 2^107.30: the doubling surcharge
  (1/(6R) + 2^-20)*(2.35-1.11)*S_ops, Markov-debited below.
- raw = 2^126.91988; 27/25 reserve gives W = 2^127.03092; claim 127.03.

## Itemised control ceiling: 110 instructions, step price 1.11

Primitives: insn 1, cmp256 5. Nine rows pay two 4x64 loads or a 192-bit
compare; the index rows pay one 8-product truncated product (17) plus the
v0+v1+v2 fold that yields the third fraction without a second multiplication
(6).

| group | rows | instructions |
| --- | --- | ---: |
| key, branch/twist, table address, jump load, twist select | 5 | 24 |
| index product, fold, near-tie test, argmin, parity twist | 5 | 36 |
| key store, 2-cycle and checkpoint compares, DP test, rotate, flags | 6 | 19 |
| meter, loop, replay tag, spill/reload, un-itemised allowance | 5 | 31 |
| **total** | **21** | **110** |

Price = 1 + 110/1000 = 1.11 for one reference addition (6M = 1074 measured
instructions, already C=1000 and not charged twice) plus that ceiling.
Doubling is 7*200 + 6*20 + 3*200 + 200 = 2320 instructions = 2.35 units: the
Jacobian a=0 formula (7M), six additions, the 3M share of the batched inversion
that restores an affine x, and 200 control. M<=200 is a limb count: 4x4
schoolbook plus Montgomery reduction is 32 multiply, 32 multiply-high, about 56
add/carry and 8 load/store classes; the 5x52 kernel's 179 is not reused.
No timing measurement is claimed; these are analytic ceilings.

## Probability

H1: T_c counts every class arrival, and the useful collision is the first
arrival at a class already occupied by another trail. On that statistic
P[T_c <= S] = 1 - exp(-S^2/(2N')), so S_ops gives model success 0.420. The
Markov early-halt allowance is random_extras/headroom = 2.18e-5 with headroom
2^127.03 - phase1. Debits: pool Chernoff 2^-10 at twice the per-node mean,
fingerprint union S(S-1)/2^193, relation 1/r, bounded scalar rejection,
detector abort 8L^2/N, endpoint transport 64*2^138/N. Success >= 0.419 > 0.40.

Doubling cannot cancel structurally: a return after m doublings needs
2^m in U = {+-1, +-lambda, +-lambda^2} mod r. Part A verifies 2^m not in U for
every m <= 2^22 exhaustively, and a chain of 2^22 consecutive doublings has
probability (2^-19.78)^(2^22). The argument therefore does not depend on the
factorisation of r-1, which this package does not reproduce.

## Storage partition

Reports route on bits 128..143 of K, which are part of H1 together with the
partition and DP bits: the routing word is uniform, independent of the branch
and DP words, and identical for two trails that reach the same class. Mean load
is 2^42 records per node; the cap is 2^43 (2.0x). Per-node bytes are
2^43*192 = 2^50.58 < 2^56, total 2^59 records, peak 2^66.58 bytes including the
2^40 auxiliary cap and the 6R table. Node setup and index construction are
charged. Overflow is a failure, never a silent drop.

## Reproduction

```
echo '{"curve":"q25","trials":1200,"R":2048,"dp_bits":6,"seed":1,"production_checks":2,"seconds":14,"pipeline":true}' | python3 -I scripts/run_experiments.py
python3 -B scripts/run_experiments.py --offline      # captured suite
python3 -B scripts/verify_certificate.py --write
python3 -B scripts/verify_certificate.py             # regenerates, compares
```

## Toy evidence (actual scope only)

Four census configurations (14,400 trials total) give mean T_c/sqrt(N) within
+3.2% of 1.2533*sqrt(R/(R-1)) and pooled P[T_c <= S_ops] >= 0.40. The 2-cycle
rate tracks 1/(6R). A 16-node pipeline routes every report with zero overflows.
The long-trail test at q29 separates doubling from the additive control. The
truncated index rule gives zero mismatches against the exact six-orbit argmin
on two million random inputs and on every orbit-twisted partner checked, with
the fallback exercised. These are finite shape tests; they do not establish
production magnitudes, rare basins or instruction timings.

## Limitations

H1 and the twist-uniformity clause are heuristics supported only at toy scale.
Instruction counts are analytic, not measured. The 2^22 direct verification of
the doubling-return condition covers all chains that can occur, but it is not a
statement about the order of 2 mod r. Memory is reported and unscored.
