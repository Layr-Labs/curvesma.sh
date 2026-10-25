# g2-invert v39 -- element-keyed rho at the collision floor

Recover `x` from `z = g^x` in the BLS12-381 GT group. Each walk step stores the exponent
coordinates `(A, D)` of `y = g^A z^D`; two records with the same element satisfy
`A1 + x D1 = A2 + x D2 (mod r)`, so `x = (A2 - A1)/(D1 - D2)` and the only failure is
`D1 = D2`, probability `1/r`.

Claim **128.9625** against `W = 13261162772257673514973661745255587682588424` Ir, the
exact sum of 17 rows in `evidence/ledger39.json`.

This revision answers ticket `6000bc0f`, which rejected v38 because (1) the production
success was 0.084, below the 0.39 floor, (2) the page table was not a wide-address
implementation, and (3) the schedule, the collision row and the certificate were
assertions rather than executed work.

## 1. Success at the floor -- the scale is arithmetic

v38 priced `S = 95803668470929860550224427498997735090 = 0.418 * sqrt(r)` steps. Under the
standard rho model (S draws without replacement from the order-r group) the collision
probability is

    p(S) = 1 - prod_{i<S} (1 - i/r)  >=  1 - exp(-S(S-1)/(2r))

which at v38's S is **0.0838**. Conditional recovery given a collision is `1 - 1/r`, so the
production success was 0.084 -- the verdict's first flaw.

v39 prices `S = 227678985548529776010180189721849206666 = 0.99428 * sqrt(r)`, the SMALLEST
integer whose lower bound reaches the floor:

    1 - exp(-S(S-1)/(2r)) = 0.390000  >=  0.39      (exact integer test, evidence/scale39.json)

`evidence/arm39.json` runs the SHIPPED protocol at the production-normalised budget
(718 = 1.00 * sqrt(521041) steps per run, 200 independent seeds, verified recoveries only):
**83/200 seeds solved = 0.415**, above the model's 0.3898; the 95 % Clopper-Pearson
interval is [0.356, 0.477], which contains the model. The arm is a consistency check on the
model; the floor is the exact bound above and no toy rate is imported into W. v38's arms ran
at `55 * sqrt(r)`, which is the objection answered.

Every S-dependent row moved with S: starts `ceil(S/65536) = 3474105614449001709139712367581928`,
cold restart 73,084,100 Ir each, occupancy `ceil(M(1-e^(-S/M))) = M` for M = 2^110.

## 2. Why the Frobenius class merge is NOT claimed

v38 withdrew the class merge and v39 does not smuggle it back. `evidence/invariance39.json`
measures both directions at production width:

| object | trials | Frobenius-invariant |
|---|---|---|
| the shipped element key (two Fp limbs of one coefficient) | 1200 | **0** |
| positive control: `min` over the twelve conjugates | 1200 | **1200** |
| the toy residue key (v37's claim) | 1200 | **0** |

The test bites: it passes a real class key and fails the shipped one, and v37's toy "one key
over 200 randoms" was a probe bug, not an invariant. A genuine class key costs 11 Fp12
Frobenius maps per step -- 635,406 Ir against a 57,082 Ir step, **12.13x** -- and buys 0.687x
the steps, i.e. **8.33x** the work against path (a)'s 2.379x at the same success. The merge
is rejected on price; the success above is earned by S.

## 3. The page table is real, and its overflow boolean is checked

`scripts/protocol37.py::PageTable`:

* a locator is FOUR opaque `u64` chunks in a `ctypes` array, assembled explicitly into a
  256-bit word (`c0 | c1<<64 | c2<<128 | c3<<192`) and disassembled on every access, so no
  wide pointer is ever dereferenced by host code;
* three 38-bit radix levels plus a 22-bit offset name **136 address bits**;
* `_idx` tests `addr >> 136` and the sign bit, sets `overflow`, and RAISES; `alloc()` opens
  with `if self.overflow: raise`, so a tripped table refuses further inserts;
* `walk()` returns `(present, locator)` and the shipped loop BRANCHES on `present` -- the
  Boolean v38 ignored is load-bearing.

Controls (`evidence/pagetable39.json`): a planted address `2^137` (138 bits) trips the
overflow and every later insert is refused; a 121-bit address round-trips its chunks; and
the same fixture against a table whose membership is always False recovers **0** matches
against **6** live.

## 4. Budget exits, collision row, certificate

* **Wasted work is charged.** `run_one` adds the steps executed inside an abandoned trail to
  `wasted_steps`; the ledger carries `abandoned_trail_waste` = 32 steps x 57,082 Ir from the
  forced-abort control (cap 200,000 Ir), so the row is live and W is an upper bound even if an
  abort were reachable. On the production schedule the trail is truncated at the remaining
  budget, so the counter is 0 there.
* **The collision row is executed.** `collision_processing` = 5 calls x 55,128,548 Ir (12
  Frobenius maps + one Fermat inversion + one verification), where 5 =
  `ceil(S(S-1)/(2r)) + 4` is the same birthday count the success probability uses. The v38
  constant 8,192 is gone; the row is 33,648x it, and doubling the stub's work doubles the row
  (`evidence/collision39.json`), which a constant cannot do.
* **The certificate checks an external invariant.** `scripts/cert39.py` does not re-sum rows:
  it (A) recomputes `x = (p^k A1 - A2)/(p^k D1 - D2) mod r` from RAW exponent inputs on 300
  constructed non-degenerate crossed collisions and checks `g^xx == z` by exponentiation
  (300/300), (B) checks the collision algebra at the 255-bit production r on 1000 random
  triples, including that a wrong shift index fails (919/919), (C) checks the floor and the
  minimality of the claim in exact integers, and (D) corrupts one input and requires the check
  to FAIL. The row sum is reported as bookkeeping and labelled as such.

## What is charged

Seventeen exact-integer rows: 9 per-step rows (dense Fp12 multiply, lookup probe, slot export,
record write and insert, six-limb key formation, six-limb coefficient update, ctypes pack,
guard and meter, widened references), address translation, the cold trail starts, the store
initialisation from the occupancy, worker starts, the executed collision row, the wasted-work
row, code and advice, and the final verification. No float enters a charged term.

## Limitations

* The uniform random element model is an assumption, stated where it is used; the arm at the
  normalised budget is consistent with it (0.415 vs 0.3898) but cannot, at n = 200, establish
  0.39 by interval arithmetic alone.
* The six-limb spans are read from the shipped object on this machine (`otool -tvV`).
* Toy fixtures use q = 523, r = 521041; the store, the detector, the recovery and the page
  table are the shipped code, only the field is small.

## Reproduce

    python3 scripts/v39.py quick     # invariance, scale, page-table controls
    V39_SEEDS=200 python3 scripts/v39.py arm
    python3 scripts/regen39.py       # ledger39.json, submission.json
    python3 scripts/cert39.py        # certificate.json (asserts ok)
