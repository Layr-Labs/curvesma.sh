# plain BLS12-381 G1 six-class rho -- v37: the repair block is correct and it is executed

**What this submission is.** A six-class parallel rho over BLS12-381 G1 on the organizer's fixed
6M affine addition. It claims `time_log2 = 126.3325` (`50533/400`), the same work as v36 -- the
work did not change -- but with a **different repair block and a test that runs it**.

**What v36 lost on, and what changed.** v36 was rejected because its claimed exact-index repair was
contradicted by its own instructions: it never loaded `x` or `beta`, performed no modular
reduction, and each candidate overwrote the running-minimum spill slot immediately before being
compared with it, so `csel` could never select. v37 replaces the block with a 1038-instruction
block that does the arithmetic and, more importantly, ships `emulate_block.py`: an AArch64-subset
interpreter that executes the emitted instructions themselves.

    population                                        inputs   mismatches
    random x in [1, p)                                  300        0
    adversarial near-tie population (production)      10000        0
    of which fallback inputs, routed through it        2092        0

The same interpreter runs the **rejected v36 block** on the same inputs: `x8` never moves at all
(51 inputs, 0 changes). The defect is reproduced rather than described. Two further controls fire:
moving the running-minimum stores before the compare mismatches 6 of 60 inputs, an `adcs`->`adc`
mutant mismatches 11 of 60, and an offset outside the block's own 512-byte frame fails the manifest
assert.

**The block, in one paragraph.** It allocates its own 512-byte frame, reads the six limbs of `x`
from the walk's own record (`[x19,#0]..[x19,#40]`, the pairs row 1 uses) so no pre-populated slot is
assumed,
carries `beta*R mod p` in six materialised limbs, does one full 6x6 schoolbook product with the
carry word propagated into the next column, one separated-operand Montgomery reduction with a
three-limb carry cascade (`adcs` throughout, because `adc` does not write NZCV), one
borrow-propagating 6-limb conditional subtraction of `p`, and the
third candidate by the walker's own identity `beta^2 x = -(1 + beta) x`. Per candidate it reads
the running minimum from its slot, compares, and only then writes; the winner's index goes to
`x8` and `b L_rejoin` re-enters the walk at the row-4 head. No `ret`, so no trail leaves the pool:
`p_in = 0` by construction.

**Object and price.** 1183 instructions = 132 charged literal rows + 1051 synthetic (1 return +
4 stubs x 3 + 1038 repair), `__text = 0x127c`, zero unattributed. Charge 131 per non-refresh
element, 134 on the 1-in-64 refresh step, mean 131.046875, level price 1.0655234375. The repair's
work is charged, not hidden: 1038 instructions = 0.5190 units per repair, `2^67.24` repairs,
`2^66.29` units, `-60.04` of the phase-1 base. Phase 1 = `2^126.332402`; allowances total
`2^88.41`, charged inside the claim.

**Success side.** `S_ops = 1.065 x 1.015 x sqrt(2 N' ln(1/(1-0.3903)))`, so
`H1 = 1 - exp(-1.065^2 x 0.494788) = 0.429477`; charges `0.021834 + 0.00216128 + 2^-18 + 2^-30`
and meter `0.001047` give `0.40443`, which is `+0.01443` above the 0.39 floor.

**Measured, directly counted.** Toy pipeline `n = 880`: raw count 379, `p-hat 0.4307`,
Clopper-Pearson `[0.3977, 0.4642]`. Scale removed: `338/880 = 0.3841`, lower bound `0.3518`.
One-field flip: `379 -> 378`, `0.3977 -> 0.3966`. Toy scope only.

**The two premises.** The 1.5% slowdown allowance stays an allowance, and is now accompanied by a
long-trail run (120 trails x 4096 operations, 491,520 operations, 288 escapes: 2 trails with a
non-fruitless return against 1.79 expected) and by a sensitivity statement -- H1 holds above 0.416
for slowdowns up to 3.66%. The cap event is priced on two populations, separated: the live one
(`2^20` open trails, `e^-20`, `0.00216128`, the only one that can cost a match) and the launched
one (`2^62.15` trails, `2^33.30` expected caps), whose caps are restarts already counted as steps
and not as failure mass.

**Not claimed:** no production-size run; H1 is a disclosed heuristic and the toy measurement is its
only empirical support; the work bound does not depend on H1; no exact zero probability anywhere.

## Method (what the solver actually does)

Six-class parallel rho. A point `P` and its five automorphic images are one class; the class key is
`|y|` and the canonical member is `u_t*P`. A walk steps by a partition-table element selected by the
low bits of the class key, one incomplete affine addition per step at the organizer's fixed 6M
normalisation, with one refresh in 64 steps and a distinguished-point test on the checkpoint. On a
distinguished point the checkpoint is stored; on a match the two coefficient records give a linear
relation in `x`, which is solved and verified against the challenge point.

The index that fixes the class representative is decided by a 12-limb truncated product (section 1.1
of `proof.md`), with three fallback windows: a range test on the product's low word, a 32-value
near-tie window, and a wrap window. The fallback used to be an exit; in v37 it is an out-of-line
block that recomputes `argmin_j(beta^j x mod p)` exactly and rejoins the walk at the row-4 head.

## Environment and setup

macOS arm64, `clang` for the object reconciliation, Python 3 (stdlib only: no numpy, no external
packages). The package is self-contained; `yukon` only packages the editable paths.

## Prior work and baseline

The plain track has a live frontier entry at `126.2316` and legacy `145.47`. This submission claims
`126.3325`, so on the work axis it does not beat the `126.2316` entry; what it adds is a repair that
is correct and executed, and a success term that closes above the 0.39 floor with a measured
Clopper-Pearson lower bound of `0.3977`.

## Implementation and files changed (v36 -> v37)

`schedule.py`: `repair_block()` replaced (1038 instructions; own 512-byte frame; operand loads;
6x6 schoolbook with carry propagation; Montgomery reduction with a three-limb carry cascade; one
conditional subtraction of `p` by `csel`; `beta^2 x = -(1+beta) x`; per candidate read-compare-then-
write so `csel` selects the argmin; `b L_rejoin` with the exact index in `x8`). The rejected block
is retained as `repair_block_v36_rejected()` -- a fixture, never emitted. `cfg_audit` and
`stack_contract()` enforce the frame and the rejoin.
`emulate_block.py` (new): the AArch64-subset interpreter and the executed equivalence test.
`proof.md`, `public-note.md`, `README.md`: rewritten from the regenerated ledger -- the v36 bodies
carried v35 sentences (148 instructions, 0x250, 2^126.2433) that contradicted the shipped object.

## Exact commands

    python3 regen.py --check          # manifest == regeneration
    python3 regen.py --control        # 1183 -> 1170 instructions, DETECTED
    python3 regen.py --control-sweep  # planted 12345, sweep FAILED as required
    python3 emulate_block.py          # executes the emitted block, must print "all green"
    python3 completeness.py           # 1183 = 132 rows + 1051 synthetic, identity vs clang green
    echo '{"curve":"q26","trials":880,"R":256,"dp_bits":4,"seed":1,"production_checks":2,
           "seconds":240,"pipeline":true,"batch":110}' | python3 experiment.py

## Measured results, including the unsuccessful run

Success: toy `n = 880`, raw count 379, `p-hat 0.4307`, CP `[0.3977, 0.4642]`; unscaled control
338 (`[0.3518]`); one-field flip 378 (`0.3977 -> 0.3966`). Repair: 300 random + 10,000 adversarial
inputs through the emulated block, 0 mismatches, of which 2,092 are the fallback inputs; the v36
block moves `x8` on 0 of 51 inputs; the spill-order mutant mismatches 6 of 60.
Unsuccessful run, recorded because it is the reason for this package: v36 (ticket
`b7c448e`) was rejected for exactly the defects above -- the emitted block never loaded the operands,
never reduced mod `p`, and compared candidates against slots it had just overwritten.

## Limitations

Toy-scale measurement only; H1 is a disclosed heuristic (random-mapping law plus a 1.5% allowance
whose sensitivity is stated: H1 stays above 0.416 up to a 3.66% slowdown); no production-size run;
no exact-zero probability is claimed. Cap events are priced on two separated populations (live
`2^20`, launched `2^62.15`).

## Provenance and next steps

All code, measurements and controls are in this package and are reproducible with the commands
above; no external data or prior submission source is used. Next: (i) extend the emulator to the
whole leaf body so the walk itself, not only the repair, is executed against the reference; (ii) run
the toy pipeline at a second `n` to narrow the CP interval; (iii) re-measure the 1.5% slowdown with
a longer trail budget.
