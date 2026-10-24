# p256-ecdlp v11 -- orbit-row plain-lineage PORT (rhosigma family abandoned)

## What this ticket is
The rhosigma family is dead on this track at the physics level (shipped-generator
census p2 0.0705 vs 0.4558 control), so this lane stops repairing it and PORTS the
campaign's proven plain construction: mu_6 orbit-row walk, refresh store, lag-3
escape, metered ledger, selfaudit -- judge-valid at 126.2364 on bls12-381-g1-ecdlp-plain
and accepted architecturally on secp256k1 (verdict cd18d785).

## The P-256 retarget, and the one structural delta
P-256 has a = -3, b != 0 (j = 1728), so its automorphism group is {+1,-1}: there is
no mu_6 orbit. Verified in `ledger_v11.py` (every nontrivial cube root of unity in F_p
fails the curve equation). The retarget is the same code with lam = 1 --
canon(c) = min(c, q-c), e in {+1,-1} -- so the class space is N = (r-1)/2, which
costs +0.79 bits of S0 (log2 S0 = 127.4924). The walk, the escape, the
recovery and the ledger form are otherwise unchanged.

## Measured success (the claim's statistic)
* 3,200 trials, matched pairs, escape armed (lag-3), budget at the exact ledger shape
  x = visits/sqrt(N) = 1.6002 (`evidence/census-design.json`).
* CP-95 lower bound 0.6826; debits 0.000336; restart loss
  0.000469; **design success 0.6818** (floor 0.39:
  MET).
* Order transfer, design shape: same statistic at 255.7x fewer classes
  gives CP-95 0.6091;
* w=2 pair census at the base shape: CP-95 0.2808 over 3200 sets
  (comparison row only -- the claim is the single-run measured row, so no +1 bit).

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

STEP_PRICE = 1 + 263/1000 = 1.263 (leading 1 = the reference
group operation; the rows are the overhead, plus a 180-insn P-256 Solinas reduction
premium charged because this track prices P-256 work in instructions, not at the
borrowed secp addition price). C = 1000. S_OPS = ceil(A_H1 * x_meas * sqrt(N')).

    log2 S0   = 127.49244191
    log2 S_OPS= 128.19973435
    log2 PHASE1 = 128.53658899
    log2 EXTRAS = 114.95291103
    log2 W    = 128.53670650   claim 128.537

Preprocessing (table build) log2 = 43.5999; memory log2 bytes =
74.5850 (reported only). Fresh-start price: 492300 insn
(492.3 ref), restart allowance log2 = 109.54.

## Rank posture (stated, not spun)
Frontier 127.237678; claim 128.537; **1.299 bits behind,
NO RANK.** That is the arithmetic of the port: the 6-fold orbit is worth 0.79 bits of S0
and the P-256 reduction premium 0.34 bits of step price, so a legal P-256 instance of
this lineage cannot contest 127.2377. The ticket is submitted as an honest, floor-met,
auditable P-256 baseline of the proven architecture rather than as a frontier move.

## Counter coverage (reviewer point: counters must span ALL call sites)
`evidence/census-design.json` fires counters inside `qwalk.step`, which is called from
exactly one place in the shipped code (`census_v20n.walk_arm`); `pair_census` and
`transfer_arm` call `walk_arm`, never `step` directly. `selfaudit.py` A7 greps every
shipped script for `qwalk.step(`/`walk_arm(` call sites and fails if an unaudited path
appears. The escape counter, the restart counter and the DP counter are incremented in
the same function that executes the step, so table construction / init / verify cannot
run the counted step uncounted.

## Audits
`selfaudit.py`: row-sum identity, claim == ceil(log2 W) at 1e-3, success arithmetic,
floor/rank direction laws, orphan-number scan over this note, counter coverage, evidence
completeness -- plus planted-flip controls (step row +1, claim -0.002, success +0.01,
note floor 0.39 -> 0.3901), all of which must fire.

## Honest limits
1. Measured at fixture DP density 2^-5 (production 2^-60): the density axis is
   NOT claimed to transfer; the fixture is a stress test of the gate.
2. Fixture column count R=64 gives a trap rate 1/(2*64) per step against production
   1/(2*2^20): the escape channel is over-fired by ~13 bits in the fixture.
3. Residual absorption is counted exhaustively at two sizes (share falls with size); no
   production basin bound is claimed.
4. Run-to-run correlation of the paired runs is not modelled.
5. No disassembly/kernel price is claimed anywhere; M_insn = 97 is structural.
