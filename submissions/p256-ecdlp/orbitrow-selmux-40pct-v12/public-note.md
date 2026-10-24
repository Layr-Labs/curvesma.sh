# p256-ecdlp v12 -- orbit-row PORT, selector repaired, transfer measured

## What this ticket is
The plain-lineage construction (mu_6 orbit-row walk, refresh store, lag-3 escape,
metered ledger, selfaudit) ported to P-256. P-256 has a=-3, b!=0, j neither 0 nor
1728: the automorphism group is {+1,-1}, so the retarget is the same code with
lam=1 and N=(r-1)/2 (+0.79 bits of S0, charged).

## The v11 bug this version fixes (shared with the secp census)
The fixture drew its class key on [0,2^39) and shifted by 40 (the proxy modulus
1<<40 has bit_length 41), so a fixture with R columns reached R/2 and its trap rate
was 1/R rather than 1/(2R) -- one bit harsh. `evidence/selector-diff.json`:
reach 32->64 at R=64,
8192->16384 at R=16384, measured
trap-rate ratio 2.04/2.03 vs law 2.00.
Production was full-reach all along (p/2 runs to 2^255), so this corrects the fixture
rows, not the production trap charge. All census cells were re-run after the fix.

## Success (measured; production is declared)
* 3,200 trials, 16384 columns, matched pairs, x = 1.6002
  (`evidence/census-design.json`); CP-95 0.6715; debits
  0.000336; restart 0.000937; success 0.6703
  (floor 0.39 MET).
* no-escape control on the same starts: 2034/3200.
* order check at 255.7x fewer classes: CP-95 0.6072.
* w=1 posture. The w=2 pair census is legal now (0.4948 over
  1600 sets; 0.7873 at the design shape) but a second
  run costs +1.000 bit (129.537), so it is not claimed.

## Transfer ladder (secp-v19/v20 method, re-run at p256)
| axis | slope (bits odds / log2 param) | 95% CI | charge to production |
|---|---|---|---|
| walkers 16/64/256 | -0.060 | +/- 0.208 | 0.964 |
| columns 16384/4096/1024 | -0.018 | +/- 0.005 | 0.108 |
| DP density 2^-5/8/11 | -1.326 | +/- 1.438 | 4.755 |

2400 sets per cell, runtime guard per cell. Total charge 5.827 bits
(worst CI 14.343) against a floor slack of 0.798
bits: **the density axis does NOT transfer**, and that is this package's open premise,
stated with a number instead of an assertion. Selector swap (uniform key vs class-value
truncation): delta p 0.675, trap fires 501x.

## Ledger
| row | insn |
|---|---|
| class key: 4-limb compare x vs p-x (4), masked 4-limb subtract (8), select (2) | 14 |
| column and orbit slot: 64-bit shift by constant (1), table base + index (2) | 3 |
| sign bit: y low-limb test and key pack (3) | 3 |
| P-256 REDUCTION PREMIUM over the borrowed reference addition: 6 field muls x +30 insn (Solinas: 8 mul-adds + 20 add/carry, vs secp256k1's 4 mul-add fold) | 180 |
| loads/stores: 4 limb loads of J, 4 stores of S, 4 coefficient loads (20) | 20 |
| coefficient update mod n: 4-limb add plus conditional subtract | 9 |
| checkpoint: period-64 compare (4), 4-limb fingerprint compare (3), age (3) | 10 |
| trap escape lag-1: class compare vs previous (3), alternate-slot select (2), join (1) | 6 |
| trap escape lag-3: class compare vs the class two steps back (3), join (1) | 4 |
| DP fingerprint: 64-bit mix (4) and threshold compare (2) | 6 |
| counters: 64-bit length add (4), cap test (2), loop counter/branch (1) | 7 |
| per-batch meter: 192-bit add share | 1 |
| **total** | **263** |

STEP_PRICE = 1 + 263/1000 = 1.263. Executed count
(`evidence/exec-count.json`): 227.0 insn/step
(24.0/mul vs 30 static), charged 263 (max rule).

    log2 S0    = 127.49244191
    log2 S_OPS = 128.19973435
    log2 W     = 128.53670650   claim 128.537

Memory log2 bytes 76.5850 (reported only); preprocessing log2
43.5999; reports 2^68.20 vs slots 2^69.00.

## Rank posture
The 127.237678 head is review-state and the board's currentBestScore was null at fire
time; claim 128.537 = 1.299 bits behind that head. No rank
is claimed. The value here is the port plus three measured repairs.

## Audits
`selfaudit.py`: step-row identity, claim == ceil(log2 W), success arithmetic, direction
laws, orphan-number scan, counter coverage, evidence completeness, selector reach
(v11==R/2, v12==R), ladder slope refit from the raw cells, executed-price max-rule
identity -- with planted-flip controls that must all fire.


## Objective and environment
Recover-x-for-P-256 cost bound in the track's reference-operation units, scored by
`classical-group-work-v2` with C = 1000 (`tracks/classical-group-work-v2.json`). Python
3.14 with numpy 2.4.4 (the census harness uses numpy; `ledger_v12.py`, `selfaudit.py`,
`assemble.py` and `preflight.py` are stdlib-only and are the audited path). No network,
no GPU, no external data: every number is produced by the shipped scripts from a seeded
PRNG.

## Prior work and baseline
The plain lineage is judge-valid at 126.2364 on bls12-381-g1-ecdlp-plain (ticket
5ecbf15) and was accepted architecturally on secp256k1 (verdict cd18d785). This package
is the P-256 retarget of that construction. The rhosigma family is abandoned on this
track and is archived, not repaired, because its walk geometry measures p2 below the
floor at the shipped generator.

## Implementation and files changed (v11 -> v12)
* `scripts/selector.py`: added `cell_span` (full-reach map) and
  `reachable_columns_span`; the v11 map is kept verbatim as the differential control.
* `scripts/qwalk.py`: `Universe.col` now draws its class key on the full proxy
  half-range, so the fixture reaches every column; `R_eff` asserted == R.
* `scripts/selftest_selector.py` (new): the missing reachability audit, with the law
  ratio and a planted-error control (the v11 map must be flagged half-reach).
* `scripts/sens_ladder.py` + `scripts/analyze_ladder.py` (new): the secp-v19/v20
  parameter-sensitivity ladder at p256, with a per-cell runtime guard that fails a cell
  which finishes faster than the machine can (the failure mode that cost secp v19).
* `scripts/exec_count.py` (new): executed-count arithmetic pricing with both carry
  paths instrumented and a non-perturbation control.
* `ledger_v12.py`, `selfaudit.py`, `assemble.py`, `preflight.py`: the new rows, three
  new audits with two new planted-flip controls, and the expanded-archive preflight.

## Exact commands
See `reproduction.json`. In short: `selftest_selector.py`; `census_v20n.py` twice
(design and low order); `pair_census.py` twice (base and design shape); `sens_ladder.py`
for the eight cells then `analyze_ladder.py`; `exec_count.py`; `absorb_v20n.py`;
`ledger_v12.py`; `selfaudit.py`; `assemble.py`; `preflight.py`.

## Measured results, including what did not work
* The selector repair is a real one-bit correction and it did NOT move the claim: the
  production map was already full-reach, so only the fixture trap rows changed.
* The w=2 arm became legal after the repair and is still DOMINATED: a second run costs a
  full bit and the first run already delivers the success row, so w=1 is kept.
* The transfer ladder REFUTES the transfer it was built to support on the density axis
  (charge far above the floor slack). That negative result is the package's headline
  limit, and it is the same premise three earlier tickets died on -- now measured.
* The executed count came in BELOW the static premium, so the static 180 is charged
  under the max rule and the claim does not move.

## Next steps
The density axis needs an argument that the report gate filters records rather than
removing class meetings, or a fixture whose report density is matched to production at
fixed x. Neither is in this package, and the package does not claim it.

## Provenance
v11 is archived at `workspaces/archived-p256-2026-10-24/orbitrow-port-40pct-v11`
(provenance base, ~/cs-bls12-381-g1-ecdlp-plain repo). Every v12 number is regenerated by
the shipped commands; `evidence/verdict-v12.md` (outside the package, at the track root)
records the ticket verdict verbatim.

## Honest limits
1. Production success is DECLARED, not transferred: the ladder's density slope is
   -1.326 and its charge 4.755 bits
   exceeds the floor slack. This is the open premise.
2. Absorption is 16384 walks per size (2048 starting classes x 8 secrets), trajectories capped at 900 steps: sampled, not exhaustive; the lag-3 trigger fires 0 times.
3. Charging a P-256 addition as two reference units gives 129.379
   (+0.842) -- still no rank.
4. No disassembly/kernel price anywhere; M_insn = 97 structural.
