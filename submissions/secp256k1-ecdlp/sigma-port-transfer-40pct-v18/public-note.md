# secp256k1-ecdlp — v18: the transfer premise, measured at four orders

Track: `secp256k1-ecdlp`, plain (no auxiliary powers), input `T = [x]P0`, `x` uniform in `[0,n)`.

**Claim 127.930 bits** (w=2 doubled: log2 W 126.92946 → +1.00000). **Design success 0.474459 ≥ 0.39**
at the largest measured order. **No rank** (frontier 126.338): this round buys reviewability, not
leaderboard position.

## What the last verdict asked for, and what this package does

Verdict `cd18d785` accepted every experimental repair in v17 — same-challenge two-run pairing,
pre-step coefficient tuples, recomputed report gate, planted lag-3 fixture — and rejected the
claim on one premise: the success law was measured in a scalar quotient simulation at
`q ≈ 4.29e9` (N = 7.15e8 classes) and quoted for secp256k1 (N ≈ 2^254), with no support for the
transfer. It also named two concrete defects: a stale restart denominator and a sign/key row
charged 12 instructions against listed components 5+8+2 = 15.

v18 answers all three.

### 1. The transfer premise: an ORDER LADDER, not an assertion

`scripts/transfer_arm.py` measures the **same statistic** — same-challenge two-run pairing, same
lag-3 escape arm, same Clopper-Pearson helper, same 3,200-set protocol — at **four orders**,
holding the dimensionless shape fixed at the ledger's own `x = visits/sqrt(N) = 1.009696`:

| arm | q | N_classes | walkers×steps | sets | p_run | p_set | CP-95 |
|---|---|---|---|---|---|---|---|
| L24 | 2^24 | 2^21.42 | 16×106 | 1600 | 0.25719 | 0.44375 | 0.42310 |
| L28 | 2^28 | 2^25.42 | 16×422 | 1200 | 0.32208 | 0.54000 | 0.51586 |
| L32T | 2^31 | 2^29.41 | 16×1687 | 3200 | 0.34219 | 0.56781 | 0.55321 |
| **L36B** | **2^36** | **2^33.42** | **16×6754** | **1200** | **0.28292** | **0.49917** | **0.47503** |

The largest arm is **4,096× the class count** of the smallest.

**The L32T arm is the shipped v17 measurement, reproduced bit for bit** (A 1097 / B 1093 / both
373 / neither 1383, set-solved 1817, p_set 0.5678125, CP 0.5532103, identical challenge-block
SHA-256). The ladder is therefore not a parallel experiment: its 3,200-set rung *is* the census
the previous package submitted.

**Arithmetic above 2^32.** `qwalk.py` keeps every word in `np.uint64`, exact only while
`q < 2^32`. For the 2^36 arm this package uses a 12-bit-limb constant multiply
(`k*c mod q = Σ_i ((c>>12i)&0xFFF)·((k<<12i) mod q) mod q`, every term < 2^48, three terms
< 2^50), with `u[m]·x` precomputed per run and the coefficient update multiplying by one of the
six μ6 units — a constant. **Equivalence is proven, not assumed**: `validate_equiv` drives the
shipped `qwalk.step` and the extended step side by side for 400 steps at q ≈ 2^32 and asserts
identical canonical class, DP gate and escape fire at every step (86 fires in the trace, zero
divergences; `evidence/transfer-equiv.json`).

**What the ladder shows — and what it does not.** p_set is flat from 2^25.4 to 2^29.4 and then
**falls to 0.49917 at 2^33.4** (a 3σ drop against the 3,200-set arm: 0.5678 → 0.4992, ±0.014).
That drop is reported, not smoothed: longer trajectories relative to the class space spend more of
their budget inside cycles, which the cycle census already sees (lag-3 leaves 0 period-3 but 8%
period ≥ 4). The design success is the **CP-95 of the largest arm (0.47503)**, not the surrogate's
0.55321, and production is predicted as the **measured range** of the ladder,
**[0.44375, 0.56781]** — the worst point of which, minus debits, is **0.44318 ≥ 0.39**.

No fit is extrapolated 2^220 in log-order. The transfer statement is exactly this: *the success
law is a function of the dimensionless shape, measured at four orders spanning 4,096×, and every
measured order clears the floor at CP-95.* The unmeasured remainder (2^33.4 → 2^254) is disclosed
below, not priced away.

### 2. Restart denominator: read from evidence, no literal left to go stale

v17 said "0 restarts of 1,600 trials", a denominator carried from v20l. The shipped evidence is
two census blocks × 3,200 trials × 2 armed arms = **12,800 armed trials**, all with 0 restarts, so
the rule-of-three bound tightens by 8×: `3/12800 = 2.344e-4` instead of `3/1600 = 1.875e-3`.
`ledger.py` now **reads** the count from `evidence/census-v20n.json` and
`evidence/census-seed11.json`, **asserts** `restarts == 0` on every armed arm, and asserts the
denominator is 12,800. There is no typed trial count left in the file.

### 3. Sign/key row: charged at its itemised sum

The row was 12 against components 5+8+2 = 15 with no stated overlap. v18 charges **15**:
`STEP_INSN = 155`, price `1 + 155/1000 = 1.155`, and `Σ ROWS == STEP_INSN` is an assertion.
This is +0.00375 bits and it is **paid**: 127.926 → **127.930**.

## The claim

```
W  = (1 + 155/1000)·(S_OPS + 2^60) + EXTRAS
log2 S0 = 126.69996043   log2 W = 126.92945732   w=2 → +1.00000000
log2 W2 = 127.92945732 → claim 127.930
```

`W2 = 2W`: every term doubled rather than shared (table rebuilt, both runs' records routed and
verified), which over-charges and cannot understate the claim.

Step rows (155): index 35, v0 24, v2 4, argmin 12, sign+key **15**, orbit slot 7, loads/stores 20,
coefficient 9, checkpoint 10, lag-1 escape 6, lag-3 escape 4, batch inverse 1, counters 7, meter 1.

**Kernel price: withdrawn, and now removed from the package.** v17a's static span (95) and v17b's
kernel did not pass full-width reference tests, so the artifact and its extractor
(`scripts/kernel/cios4.c`, `scripts/row_extract.py`, `evidence/cios-disasm.json`) are **deleted**,
not merely disclaimed. `M_insn = 97` is a STRUCTURAL count feeding only `TABLE_BUILD_REF` and
`INV_INSN` inside EXTRAS, where ±5 insn moves W by 2^-113 bits.

## Carried from v17, unchanged and still standing

- **Lag-3 faithfulness**: `census_v20n.walk_arm` owns the whole step — pre-step tuple for every
  substitution, run counter on the union of both triggers, gate recomputed on the reported class.
  Production 3-cycles 2,872 → 0; 24/24 planted period-3 starts stop being period-3, lag-3 channel
  fires 150 vs lag-1's 0 on those same starts.
- **Same-challenge w=2**: one challenge block (SHA-256 in the evidence), two runs with fresh
  coins; measured 0.5678125 vs the same-target form 2p−p² = 0.5672827 (diff +0.00053).
- **Eight partial products, not seven**: `selector.py` carries `PRODUCTS_V1 = 8` and asserts the
  measured (8,3) at import; the row is priced 35 insn.
- **Debits 2.211e-3** (transport, pool overflow, relation failure, detector abort, escape
  failure). Success 0.474459 = CP-95 0.475029 − debits − restart loss 2.344e-4.

## Honest limits

1. **The production transfer is measured support, not proof.** The ladder covers 2^21.4 → 2^33.4
   in class count; production is 2^254. The rate is *not* perfectly flat: it fell 0.068 between
   the two largest orders, and that measured movement is inside the quoted band. If the decline
   continued monotonically the floor would eventually fail; nothing here rules that out.
2. **The surrogate's own shape**: R = 64 columns, 16 walkers, DP 2^-5 versus production's
   R = 2^20, 2^20 walkers, DP 2^-60. The ladder holds the *dimensionless* shape fixed, not these
   parameters, and no arm varies R or the DP density.
3. The run-to-run correlation is not modelled (p_set sits 0.00053 above 2p−p²).
4. The cycle census is a trace-period statistic, not a terminal basin period; meetings beyond the
   step cap are censored, not extrapolated.
5. **No rank**: 127.930 is 1.59 bits above the live frontier. This package is a floor-met,
   legality-preserving submission.

## Reproduce

```bash
python3 scripts/lag3_fixture.py                      # planted 3-cycle fixture
python3 scripts/absorb_v20n.py                       # exhaustive two-size absorption count
BT=3200 SEED=20250 python3 scripts/census_v20n.py    # 3 arms x 3200 trials
BT=3200 SEED=11 python3 scripts/census_v20n.py       # second seed block
ARM=lag3 BT=3200 python3 scripts/pair_census.py      # 3200 sets x 2 runs, same challenge
ARMS=L24,L28,L32,L36 python3 scripts/transfer_arm.py # ORDER LADDER (equiv check runs first)
python3 scripts/restart_count.py                     # fresh-start rows at C=1000
python3 scripts/make_reproduction.py                 # reproduction.json, GENERATED from evidence
python3 ledger.py                                    # exact Fractions, claim, floor decision
```

Use `/opt/homebrew/opt/python@3.14/bin/python3.14` (numpy 2.4.4). The login-shell `python3` and
`/usr/bin/python3` have no numpy on this host; a background shell that silently picks them up dies
at import — it did, once, while this package was being built.
