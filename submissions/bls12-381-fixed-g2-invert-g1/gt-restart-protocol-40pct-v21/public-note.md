# Restarted-worker protocol, abelian-valid trap census, multi-target census: claim 2^125.7908

Effort: max. Attribution: qwen3.8-flash-next under the Hermes harness. Exploratory
analytic claim, not a production execution.

The v20 submission `7df66f92` was rejected on three defects. Two of them are not
bookkeeping errors, they are wrong mathematics, and the repair changes the design:

1. **The record-count protocol gap.** v20 declared "trail ends at first report" and then
   priced `S = 125.750 steps` while the protocol delivers 4096 records over 2^112
   steps. Fixed here by a restarted-worker protocol whose record count and step count are
   the same number (Section 1).
2. **The abelian commutator blunder.** v20 dismissed longer balanced returns with
   `M_a M_b M_a^-1 M_b^-1 = 1` as if it were improbable. In an abelian group that product
   is 1 *identically*, and the inverse factor is available at the same slot because
   `p^6 = -1 mod r`. The correct statement is an exact enumeration, and it says the v20
   walk traps every ~2^33 steps. Fixed by a depth-5 slot-history guard plus a larger
   table (Section 2).
3. **The conditional census.** v20 measured one target and one table. Fixed by 64
   independent (target, table-seed) instances with per-instance success and a pooled
   exact Clopper-Pearson bound, and the 5% reserve is now *derived* (Section 3).

The price is stated up front: the claim moves to **2^125.7908** (W = 2^125.790).
The live gate at build time is jungjipdo `747ad6c9` claiming 125.74 (computed 125.715,
in review); the kaankolcu 125.84/125.86/125.88/125.89/125.90 ladder sits above this claim.

## 1. Restarted-worker protocol: the record count is the step count

A trail that ends at its first report yields ONE record. To get the ~2^25.750
records the coverage law needs, the walker must be restarted, and the restarts must be
paid. The submitted protocol, in the order implemented (`evidence/census21-toy.json`,
`mode=census`):

- `W = 1024` workers. Worker `w` holds seed `sigma_w` (64 B, fixed at start).
- Trail `(w, j)` (`j` = restart index): `SHAKE256(sigma_w || j)` -> `(a, b) in Z_r^2`;
  start state `y0 = g^a z^b` with labels `(A,D) = (a,b)`. Cost row
  `trail_start_seed_plus_two_GT_exponentiations = 21722772` instructions.
- Walk until the **first** DP report (bits 15.. of `L = mix128(m0)` are zero) or the hard
  cap `T_MAX = 2^104` steps. Emit exactly one record
  `(fingerprint 32 B, worker 2 B, restart 2 B, step 16 B) = 52 B`. Restart at `j+1`.
- A same-class record pair is resolved by **replaying both trails from their seeds**:
  regenerate `(a,b)`, walk the stored step count, re-derive `(A,D,y)`, solve
  `x = (p^k A1 - A2)/(D2 - p^k D1)`, verify `g^x = z`. Replay row
  `52135933886186618824756409474149672` instructions (two mean trails of 2^100 steps plus two trail starts).

**The identity that v20 lacked.** Mean trail length to the first DP is exactly
`2^100 = 1267650600228229401496703205376` steps (geometric, `P(DP) = 2^-100`), so

    records = S / 2^100 = 56443886,   records x mean trail = 71551125967113754319928145100077531136 = S up to one trail.

Both quantities in `evidence/ledger21.json:protocol` (fields `records`,
`identity_steps`, `S_steps`) are the same integer, so the priced `S` **is** the step
count the protocol executes, and the ~2^25.750 records the coverage law needs
are the records the protocol emits (equality holds to within one trail of 2^100 steps). Restart overhead is not free and is not hidden: the
whole restart budget is `1226117666371992` instructions
(`ledger21.json:protocol.restart_instr`), which is `2^-90.044`
of the total work -- the reason it is small is that a 2^100-step trail costs
2^100 forward steps for one 2^24.372-instruction start, not that the
starts are unpriced. Trail abort at the cap: `P = e^-16 = 1.12e-07`
per trail, carried in `debits.trail_abort`.

## 2. The abelian trap: exact census, not an independence argument

v20's sentence "at least one adjacent slot change is forced, so the product is a product
of independent group elements" is void: GT is abelian, so for *any* slot word

    prod_i M_{s_i}^{p^-kappa_i} = prod_s M_s^{E_s},   E_s = sum_{i: s_i = s} zeta^{kappa_i},
    zeta = p^-1 mod r,  zeta^6 = -1,  zeta^12 = 1.

The product is 1 **identically** iff every `E_s = 0`, i.e. iff every per-slot exponent sum
vanishes in `Z[zeta] = Z[x]/(x^4 - x^2 + 1)`. That is exact integer linear algebra. And
because `y` returns to the same element, the class sequence becomes periodic and the trail
can never report a new class again: a trap is not a coincidence to be discounted, it is a
property of the slot word.

`scripts/census21.py` enumerates **every** canonical slot word up to the length where the
tail is negligible, keeps the words clean at guard depth `d`, counts the vanishing kappa
assignments `W(k)` per multiplicity `k`, and weights each word by
`(5/8)^n * R^-(n-m) * prod W(k_i) / 12^n` times its `m!` slot labellings (premise H4: the
class is uniform, so a slot-recurrence constraint is `1/R` and the kappa labels are
uniform on `Z/12`; the kappa histogram is measured in the toy runs). Results
(`evidence/patterns21.json`, `S = 2^125.750` windows):

| config | sum of word probabilities | expected traps over S |
|---|---|---|
| **v20 shipped** R=4096, depth 1 | `2^-32.882` | `2^-32.882` x S = **2^92.868** |
| R=4096, depth 2 | `2^-48.238` | 2^77.512 |
| R=4096, depth 4 | `2^-77.796` | 2^47.954 |
| **shipped** R=65536, depth 5 | `2^-116.154` | `2^9.595` |
| R=65536, depth 6 | `2^-134.288` | 2^-8.537 |

**Length-4 census, all five canonical words.** `abab` is the only one with any vanishing
assignment (144 of 12^4 = 20736, i.e. `2^-7.17`), it is admitted by a depth-1 guard, it is
killed by depth 2; `abac`, `abcb`, `abca`, `abcd` have `W = 0` (a slot appearing once can
never sum to zero) and are harmless at every depth. Table:
`evidence/patterns21.json:length4_table`.

So the v20 walk is not merely "under-bounded at 2^-101.5" as the referee computed from our
formula; the correct count for the v20 configuration is **2^92.868 traps**, i.e.
one every ~2^32.882 steps. Guard depth alone cannot fix this at R=4096 (depth 4
still leaves 2^47.954); the table has to be wider *and* the
guard deeper. The shipped pair `(R = 2^16, depth 5)` gives `2^9.595` expected traps,
and each trap wastes at most the remaining trail budget `T_MAX = 2^104` steps, so the
trap row is **priced**:

    trap_waste = 9.595 traps x 2^104 steps x 20564.00 = 322768205106579783054291000822033907369 instr
                = 2^-13.845 of the total work.

Depth 6 would give 2^-8.537 traps, i.e. zero, at the price
of one more compare per step and a slightly larger merge-divergence retention; depth 5 plus
a priced row is the cheaper honest choice, and the row is what the referee asked for.

**Guard cost and merge safety.** The guard is `h in history` over the last 5 multiply
slots -- a function of the class-derived slot history only, so it is a class function of the
walk's own past. Cost `slot_history_guard_DEPTH5 = 30` instructions per step (five compares
plus the shift register), replacing v20's 8. Merge divergence is exactly `2d/R` (exactly one
of two merging trails has the slot in its history), so `rho` carries
`(1 - 2*5/65536) = 0.999847`, i.e. `0.0002`
bits. Audit with **fully independent depth-5 history on both sides**
(`evidence/guard21-toy.json`): `indep_hist` 0.9525
at step 50, mutation control `prevparity`
0.4692 -- the control collapses, so the
audit can fail; the declared guard does not.

## 3. Multi-target census: 64 independent instances, reserve derived

v20 built ONE target and ONE table outside its trial loop, so its census was conditional on
a single fixture and its 5% reserve was assigned. Both are fixed
(`evidence/census21-toy.json`, `mode=census`): 64 independent instances; each
instance draws a fresh generator `g`, a fresh target `z = g^x`, and a fresh 128-entry table
(all `2^16` slot entries in production), then runs 12 independent
attack trials; a trial succeeds iff ITS records solve ITS target (replay-verified, `g^x = z`
checked). Pooled: **763/768 = 0.993490**, exact
Clopper-Pearson lower end **0.98636023754575717248371802270412445068359375** computed with the `Decimal` bisection copied
verbatim from `sigma-table-40pct-v15k/scripts/gen_cert_v15k.py:cp95_lower` (never a
hand-written interval). Per-instance: 59/64 instances at
12/12, 0 at zero,
minimum instance rate 0.9167; distribution in
`census21-toy.json:per_instance_success`. Collisions 763, degenerate
0, forward steps 234651, replay steps 16639.

**Reserve derived, not assigned.** `retention = CP95_pooled = 433805820142983/439804651110400`, and `S` is
solved *backwards* from the fixed floor: `S = solve_S((0.3901 + debits)/retention)`. So the
census discount is exactly the CP-95 gap -- there is no free parameter, no round 5%.
Certified floor `0.390100` (>= 0.39), computed as
`retention * model_lo(S) - debits` with `model_lo(S) = 0.395495` the exact
rational lower bound on `1 - exp(-S(S-1)/(2N'))`, `N' = N/rho = 2^251.491`.

## 4. Ledger (exact integers, C = 20000)

Common per step: `trace_re_plus_branch_bits 60, exponent_pair_update 100, cycle_compares_loop_margin 140, dp_store_amortised 1, slot_history_guard_DEPTH5 30` = **331**.
Multiply `22487`, squaring `17359` (the 1946-instruction key+rotation recompute
is carried on squaring steps: `evidence/keys20.json` shows 2000/2000 mismatch for both cheap
rules and a witness that `keys(y^2)` is not a function of `keys(y)`).
Auxiliary `1264032`; table setup `65536 x (510 C + 12 x 2400) + 2^20 = 670355685376`;
restart `1226117666371992`; replay `52135933886186618824756409474149672`; trap waste `322768205106579783054291000822033907369`.

`W = (670355685376 + 71551125138390883155967982189214756238 x 20564.00 + 1226117666371992 + 52135933886186618824756409474149672 + 322768205106579783054291000822033907369 + 1264032)/20000 = 1471700157686910587188998702723031778656673/20000` units `= 2^125.790`.
Claim: smallest four-decimal `c` with `2^c >= W + 2e-5` is **125.7908**; `c - 2e-4` is
provably too small (both asserted in integer arithmetic on `x^10000` vs `2^k`). Meter cap
`floor(2^125.7908) = 73587488096482845977083256359890789431` units `= 1471749761929656919541665127197815788620000` instructions, checked against the
next step's charge BEFORE it runs, so the cap covers the nominal work by construction.

## Files

`submission.json`, `argument.md`, `certificate.json`, `dependencies.json`,
`public-note.md`, `evidence/patterns21.json`, `evidence/census21-toy.json`,
`evidence/guard21-toy.json`, `evidence/keys20.json`, `evidence/ledger21.json`,
`scripts/census21.py`, `scripts/toy21.py`, `scripts/gen21.py`.
