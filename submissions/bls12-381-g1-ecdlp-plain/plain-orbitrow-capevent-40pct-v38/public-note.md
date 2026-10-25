# plain BLS12-381 G1 six-class rho -- v39: the row-6 operand interface is now register-complete

**What this submission is.** A six-class parallel rho over BLS12-381 G1 on the organizer's fixed 6M
affine addition. It claims `time_log2 = 126.3406` (`631703/5000`). v38 (never fired) held the same
walk with a cheaper listing; v39 pays 12 instructions to make the operand interface of the fixed
6M reference addition *exhibit* every operand it consumes.

**What the last two verdicts named, and what changed.** `42a2308b` and `73c192b7` both named the
row-6 operand interface: v38 loaded the table x-coordinate into `x0..x6` and then overwrote those
same six registers with the state y-coordinate before the declared reference-addition boundary,
so the boundary's declared inputs contained no surviving table x, and the emitted boundary was only
a comment. That is not a scheduling bug. A 6M addition consumes 24 point limbs at once (table x,
table y, state x, state y); the walk has 30 usable registers (`x18` is platform-reserved) and 12
of them are pinned for the whole span, so `24 + 12 = 36 > 30` and **no load order and no
interleaving exists**. v39 relocates six pinned values -- `PH0`, `PTOP`, `END`, `TBL`, `CONST`
(each last used before row 6) and `RPTR` (after the state point has been read out of the record,
restored before the result stores) -- to frame slots across the boundary: 6 stores + 6 loads, 12
charged transfers. The boundary then sees

    table x  x0,x1,x2,x3,x4,x6      table y  x7,x8,x9,x23,x24,x25
    state x  x26,x27,x28,x29,x30,x16 state y x17,x20,x21,x22,x19,x5
    modulus  x10,x11,x12,x13,x14,x15

30 distinct registers, 30 values live, 0 free. This is machine-checked, not asserted:
`schedule.py interface_budget()` replays the emitted instructions, reconstructs which register
holds which limb at the boundary, and fails if any of the 30 operands is absent, duplicated or
outside the pool. Control I plants the v38 streaming form and the check fires with
`table x operand TX0 is ABSENT at the [C] boundary`. Pinned preservation is proved per register
with the relocation bracketed: the value goes to slot `m<off>`, the register is reused, and the
slot is restored before the register's next use (the liveness model carries slot values, so a
relocation is a move, not a clobber).

**Object and price.** 1195 instructions = 144 charged literal rows + 1051 synthetic (1 return +
4 stubs x 3 + 1038 repair), `__text = 0x12ac`, zero unattributed. Charge 143 per non-refresh
element, 146 on the 1-in-64 refresh step, mean 143.046875, level price 1.0715234375 -- row 6 costs
32 rather than 20, i.e. `+13` on the span and `+0.0081` bits on the claim (`126.3325 -> 126.3406`).
The repair's work is charged, not hidden: 1038 instructions = 0.5190 units per repair, `2^67.24`
repairs, `2^66.29` units, `-60.04` of the phase-1 base. Phase 1 = `2^126.3324`; allowances total
`2^88.74`, charged inside the claim.

**Success side.** `S_ops = 1.065 x 1.015 x sqrt(2 N' ln(1/(1-0.3903)))`, so
`H1 = 1 - exp(-1.065^2 x 0.494788) = 0.429477`; charges `0.021834` (more than 8 matches) +
`0.000364064` (cap, `25e^-25` per age-conditioned live trail x `2^20` live trails) + `2^-18` +
`2^-30`, and meter `0.001047`, give `0.40622`, which is `+0.01622` above the 0.39 floor.

**Measured, directly counted.** Toy pipeline `n = 880`: raw count 379, `p-hat 0.4307`,
Clopper-Pearson `[0.3977, 0.4642]`. Scale removed: `338/880 = 0.3841`, lower bound `0.3518`.
One-field flip: `379 -> 378`, `0.3977 -> 0.3966`. Long-trail run: 120 trails x 4096 operations =
491,520 operations, 288 escapes, 2 trails with a non-fruitless return against 1.79 expected,
which is the sensitivity statement for the slowdown allowance (H1 stays above 0.416 up to a 3.66%
slowdown). Toy scope only.

**Still open, named honestly.** (i) The fallback rate is a *structural* `2^-59` per canonicalisation
(32/`2^64` window x the three fallback tests); the union bound over range, wrap and near-minimum
tests under the production walk distribution is stated as the sum of the three per-step terms and
is not re-measured at production scale. (ii) The cap charge is the judge's own age-conditioned
formula with the design constant `lambda = 25`; the measured long-trail data is a consistency
check at toy scale, not a production tail measurement. (iii) The slowdown allowance is a disclosed
heuristic with a measured sensitivity, not a proof.

**Not claimed:** no production-size run; H1 is a disclosed heuristic and the toy measurement is its
only empirical support; the work bound does not depend on H1; no exact zero probability anywhere.

## Method (what the solver actually does)

Six-class parallel rho: a point and its five automorphic images are one class, the class key is
`|y|`, the canonical member is `u_t*P`; one incomplete affine addition per step at the organizer's
fixed 6M normalisation, one refresh in 64 steps, a distinguished-point test on the checkpoint. The
class index is decided by a 12-limb truncated product (proof.md section 1.1) with three fallback
windows (range, 32-value near-tie, wrap); the fallback is an out-of-line block that recomputes
`argmin_j(beta^j x mod p)` exactly and rejoins the walk at the row-4 head.

## Environment and setup

macOS arm64, `clang` for the object reconciliation, Python 3 (stdlib only: no numpy, no external
packages). The package is self-contained; `yukon` only packages the editable paths.

## Prior work and baseline

The plain track's live frontier entry is `126.2316` (legacy `145.47`). This submission claims
`126.3406`, so it does not beat the frontier on the work axis; what it adds is the closed row-6
operand interface, a repair block that is correct and executed, and a success term above 0.39.

## Implementation and files changed (v38 -> v39)

`schedule.py`: row 6 rewritten as the spill-based interface (32 charged transfers: 1 entry fetch,
12 operand loads, 6 result stores, 1 product store, 12 frame relocations); the leaf frame grows
from 16 to 64 bytes; `interface_budget()` and control I added; `pinned_preservation()` extended to
prove relocations. `experiment.py`: ledger row 6 = 32, charge 143/146, `STEP_PRICE` 1.0715234375,
`CLAIM = (1263406, 10000)`. `completeness.py`: price identity re-derived at 143.046875. `regen.py`:
claim text, title and score regenerated from the artifacts, with the literal sweep as the guard
against a prose constant.

## Exact commands

    python3 schedule.py               # liveness, census, pinned preservation, interface_budget
                                      # + 9 planted controls, each must fail the right check
    python3 completeness.py           # 1195 = 144 rows + 1051 synthetic, identity vs clang green
    python3 regen.py --check          # manifest == regeneration
    python3 regen.py --control        # 1195 -> 1194 instructions, DETECTED
    python3 regen.py --control-sweep  # planted 12345, sweep FAILED as required
    python3 emulate_block.py          # executes the emitted repair block, prints "all green"
    echo '{"curve":"q26","trials":880,"R":256,"dp_bits":4,"seed":1,"production_checks":2,
           "seconds":240,"pipeline":true,"batch":110}' | python3 experiment.py

## Measured results, including the unsuccessful run

Success: toy `n = 880`, raw count 379, `p-hat 0.4307`, CP `[0.3977, 0.4642]`; unscaled control
338 (`[0.3518]`); one-field flip 378 (`0.3977 -> 0.3966`). Repair: 300 random + 10,000 adversarial
inputs through the emulated block, 0 mismatches, of which 2,092 are the fallback inputs; the v36
block moves `x8` on 0 of 51 inputs; the spill-order mutant mismatches 6 of 60; `adcs`->`adc`
mismatches 11 of 60. Interface: 30/30 operands in distinct registers, 0 free; the v38 streaming
form is absent at the boundary (control I).
Unsuccessful runs, recorded because they are the reason for this package: v36 (ticket `b7c448e`)
was rejected because its emitted block never loaded the operands, never reduced mod `p`, and
compared candidates against slots it had just overwritten; v37 fires 3-5 and v38's predecessor
were rejected on the row-6 operand interface and on the fallback-union rate.

## Limitations

Toy-scale measurement only; H1 is a disclosed heuristic (random-mapping law plus a 1.5% allowance
whose sensitivity is stated: H1 stays above 0.416 up to a 3.66% slowdown); no production-size run;
no exact-zero probability is claimed. Cap events are priced on two separated populations (live
`2^20`, launched `2^62.15`). The fallback union bound is per-step and structural.

## Provenance and next steps

All code, measurements and controls are in this package and reproduce with the commands above; no
external data or prior submission source is used. Next: extend the emulator from the repair block
to the whole leaf body; measure the fallback union on production-shape uniform inputs; extend the
long-trail run past `10^6` operations for a nonzero-support cap tail.
