# Distributed fixed-step walks with a 64-bit-local class table (BLS12-381 G1, plain DLP) — v10

Claim: `time_log2 = 126.18013` (12618013/100000; computed total
W = 2^126.1801271), success >= 0.39 for every challenge, memory
< 2^132.284 bytes. Full argument in `proof.md`; all checks in `experiment.py`.

This package replaces the v8/v9 single-machine design that used a 128-bit MMIO
bank-selector register to reach a 2^131 B table. The exploratory review held
that the protected model does not supply that device. Here the table is
distributed over 2^71 table servers, each with a local bank below 2^61.4 B, so
**every head word, arena offset and next pointer is an ordinary 64-bit local
address** and every server instruction is ordinary scalar AArch64. There is no
selector register, no `dsb`, no bank window and no `address_mapping` term.
Transfers are charged as an explicit per-point message: 4 instructions to form
and enqueue a 56-byte message, 3 per network level (two levels) and 6 to
receive. The two global counters are 128-bit words in shared memory updated
with atomic scalar `ldadd` instructions (low word per 2^16 comparisons, high
word on wrap; no lock and no retry loop); the deterministic block lag is at
most 2^71 * 2^16 = 2^87 comparisons and is inside the probe allotment.

## Ledger

* Per point: generator 41 (18 reference operand traffic + 12 batched-addition
  bookkeeping + 2 record loads + 7 key canonicalisation + 2 digits) + message 4
  + network 6 + receive 6 + server insert 13 = **70**; plus 2*2/128 two-sided
  unrolled loop control + 235000/2^26 batch share + 0.002 dispatch =
  **70.03675177288055**; price 1 + 70.03675177288055/2000 =
  **1.0350183758864402**.
* Table zeroing 0.500244 instructions per 8-byte head over S_tot = 2^128
  heads (2^57 per server) = 2^116.03492 units. Comparison price
  9 + 4*2^-16 + 3*2^-64 (9 common, 4 block path every 2^16, 3 carry path every
  2^64); E[X] = n(n-1)/(2 S_tot) = 2^123.2567; cap_eff =
  floor(1.005 E[X]/2^64) * 2^64, charge (cap_eff + 2^87) * price / 2000 =
  2^115.46801 units.
* Extras: 2^104 candidates + 2^90 network/queue set-up + 2^77 walk set-up
  + 2^76 server headers + 2^58 recovery + 2^20 output = 2^104.000088.
* W = n*price + zero + probe + extras = 2^126.1801271 <= 2^126.18013.
* Inspection overrun < 2^-31.015 (multiplicity-aware Markov + Chebyshev);
  candidate overrun < 2^-31.01 (multiplicity-weighted Chebyshev on 164-bit
  fingerprints); success lower bound 0.3900009991 >= 0.39.

## Experiments (`experiment.py`, stdlib, runs in the bounded container)

Part A: production algebra, key rows, AArch64 encodability, split
GEN_key -> MSG -> NET -> RCV -> SRV_insert -> F data-flow emulation, and a
planted end-to-end recovery on BLS12-381 (21 entries, x recovered). Part A2:
production-bit H' consistency evidence on 64 real-curve AP walks x 33 points
(2112 real points): 199 bit-balance tests max |z| = 2.35, 64 pairwise
chi-square tests max 4.844, Poisson occupancy chi2 2.91 (S=2^6) and 5.97
(S=2^10), truncated 164-bit fingerprint collisions on high/middle/low slices
for b = 12..20 max |z| = 1.887, and within-walk vs cross-walk rates at b=16
consistent with 2^-16. This is consistency evidence on a scaled sample, not a
proof of H'. Part B: toy-curve attack plus the exact scalar model. Part C: the
exact-fraction ledger with all listing-count and cap checks.

## Reproduce

```
python3 experiment.py <<< '{"claim":126.18013,"production_checks":1,"production_random":12,"a2_walks":64,"a2_steps":32,"toy_seconds":5,"scalar_seconds":3,"curve":"q26","toy_walks":32,"scalar_walks":96,"scalar_rbits":22,"seed":1}'
bash scripts/run_experiments.sh bls12-381-g1-ecdlp-plain submissions/bls12-381-g1-ecdlp-plain/ap-v10-distributed-12618004
```

## Provenance

Builds on the live accepted frontier 07207131 (126.18324); the six-class
collision law and meter conventions come from the accepted 126.1862 package
c8b6fee8. The arithmetic-progression walk design d430ee99, the
multiplicity-aware tails 6fd2f92e, the chained table 4935b1dd, the counter
antecedents 40c0cd3c, and the v7/v8 listing discipline 2bf3b073 / 242c1762 are
materially reused. The address mechanism of edbdb065 and the BANK_SELECT
listing are NOT reused; they are replaced by the distributed protocol above.
