# proof.md — plain BLS12-381 G1 six-class parallel rho, v37

v37 is v36 with the exact-index repair REWRITTEN and EXECUTED. v36 was rejected because the
repair block it claimed did not compute what it said: no operand load, no carry propagation, no
modular reduction, and a compare that read the very spill slots the same candidate had just
overwritten, so `csel` never selected anything. The arithmetic of the attack did not change and
does not need to: the floor is met by the measurement and by the scaled schedule. What changed is
that the object now contains a block that can be run, and a test that runs it.

Nothing below is inherited text. Every number in this file is produced by `regen.py --write` from
the shipped tree (`regen-report.txt` is its extraction table), by `experiment.py` when it is run,
or by `emulate_block.py` when it executes the emitted instructions.

## 1. Classes, canonicalisation, and the index the repair must reproduce

A point `P = (x, y)` has six images under the automorphism orbit. The class key is `|y|` and the
canonical representative is `u_t · P` where `u_t` is the unit whose action moves the x-coordinate
to the smallest of `x`, `beta x`, `beta^2 x` (mod p). The walker computes that minimum directly:

    bx  = beta * x % q
    b2x = (-x - bx) % q          # beta^2 x = -(1 + beta) x
    j   = argmin(x, bx, b2x)     # canon(), experiment.py

Section 1.1 decides `j` from a 12-limb truncated product and falls back to that exact rule on the
three windows (range test, 32-value near-tie window, wrap window). The fallback is entered by
`b.hs L_fallback` / `b.ls L_fallback`; what happens there is the subject of section 5.

## 2. Walk, batch, link register

Unchanged in structure: the walk steps by a table element whose index is the low bits of the class
key, one 6M affine addition per step with the organizer's fixed normalisation, one refresh in 64
steps, and a per-batch prologue that allocates 16 bytes and saves the link at `[sp,#8]`. The
measurement (section 7) supports the step count; the object enumeration supports the price.

## 3. Ledger and prices

The charge is 131 instructions per non-refresh element (129 literal rows plus 2 amortised
per-batch rows, meter and Fermat share) and 134 on the 1-in-64 refresh step (131 plus 3 refresh
stores). The 64-cycle mean is

    131.046875 = (63 x 131 + 134) / 64

and the level price is `1.0655234375`. Phase 1 is `1.0655234375 x (S_ops + 2^60) = 2^126.332402`,
with headroom `2^112.48` and a meter early-halt bound of `0.001047`.

## 4. The shipped object

The span is one single-entry function of **1183** object instructions: **132** charged literal rows
in emit order plus **1051** synthetic lines — one standalone return, 4 three-instruction exit stubs
(12), and the **1038**-instruction repair block (`13 + 1038 = 1051`). `__text = 0x127c = 4 x 1183`
bytes, ending at the last instruction, with zero unattributed instructions (`completeness.py`;
three planted enumeration defects each detected). Every record field is read before its own first
write: offsets 0, 16, 32, 48, 64, 80, 96, 104, 120, 136 are loaded, in emitted order, before the
first store to each (10/10).

The repair block's own work is charged, not hidden: 1038 instructions = **0.5190 units** per
repair at the structural `2^-59` per canonicalisation, so `2^67.24` repairs over
`S_ops = 2^126.24084` steps cost `2^66.29` units, i.e. `-60.04` of the phase-1 base. That is an
EXPECTED-COUNT bound on the probability of at least one repair; it is not converted into a
probability of none, and no exact zero is claimed anywhere.

## 5. The repair block, as emitted (schedule.repair_block)

The block allocates its own 512-byte frame (`sub sp,sp,#512` … `add sp,sp,#512`), so the leaf
prologue's 16-byte link frame is untouched and every offset it uses is inside the frame it owns
(largest offset 496). It saves every register it touches except `x8`, and `x8` is the result.

Per repair:

1. **Load the operands.** The block's first action after its frame allocation and register saves is
   to read the six limbs of `x` from **the walk's own record** at `[x19,#0]..[x19,#40]` (the same
   three `ldp` pairs row 1 uses) and copy them into its frame: it does not assume any slot is
   pre-populated, and the interpreter test feeds `x` through that record, not into the frame. `beta`
   is carried as `beta·R mod p` in six materialised limbs (`movz`/`movk`), and `p` likewise.
2. **Full product with carry propagation.** One 6x6 schoolbook: for each multiplier limb, six
   column steps `(C,S) = T[i+j] + a[j]·x_i + C` with `mul`/`umulh`/`adds`/`adc`/`adc`, the carry
   word written into the next column (`cset`/`add`), so no limb is dropped or folded.
3. **Montgomery reduction.** Separated operand scan: `m = T[i]·n0' mod 2^64`, then `m·p` added at
   limb offset `i` with a three-limb carry cascade. Result `= beta·R·x·R^-1 = beta·x mod p`.
4. **Reduction to canonical form.** One borrow-propagating 6-limb `subs`/`sbcs` cascade against
   `p`, then ONE conditional subtraction selected by `csel` on the carry (the `op_point_construct`
   lesson: never a bare subtract). The stored candidate is in `[0, p)`. Every multi-limb carry
   chain uses `adcs` except its final limb: `adc` does not write NZCV on AArch64, so a chain of
   plain `adc` would not propagate the carry.
5. **Third candidate.** `beta^2 x = -(1 + beta) x` — the walker's own identity: one
   borrow-propagating 6-limb subtract `p - x - bx` and a masked add of `p` when it borrows.
6. **argmin, in the right order.** Per candidate: six loads of the running minimum from its frame
   slot, six loads of the candidate, the 6-limb `subs`/`sbcs` cascade, `cset` to a flag, six
   `csel` into the winner, six stores, and one `csel` for the index. The stores are AFTER the
   compare, so the spill slots cannot feed the comparison they belong to.
7. **Rejoin.** `ldr x8,[sp,#496]` puts the exact index in `x8` — the register row 4 reads — and
   `b L_rejoin` re-enters the walk at the row-4 head. There is no `ret` in the block, so no trail
   leaves the pool: `p_in = 0` by construction and no abandonment debit exists as a term.

`cfg_audit` enforces: no `ret` in the block, `L_rejoin` inside the leaf body, the block ending in
`b L_rejoin`, the block allocating its own frame, releasing it exactly once, and
`stack_contract()` finding every offset in-frame.

## 6. The equivalence is EXECUTED

`emulate_block.py` is a 64-bit AArch64 subset interpreter (`ldr`/`str`/`stp`/`ldp`, `mov`,
`movz`/`movk`, `mul`/`umulh`/`madd`, `add`/`adds`/`adc`, `sub`/`subs`/`sbcs`, `cmp`, `cset`,
`csel`, `neg`, `and`, `b`/`b.cond`) with ARM carry semantics (`sbcs` subtracts `NOT(C)`). It runs
the emitted instructions themselves on the walker's inputs and compares with `exact_j`:

| population | inputs | mismatches |
|---|---|---|
| random `x` in `[1, p)` | 300 | 0 |
| adversarial near-tie population (production generator) | 10000 | 0 |
| of which FALLBACK inputs routed through the block | 2092 | 0 |

Controls, all run in the same pass:

* **argmin control.** An input whose third candidate is strictly the smallest must move `x8`. The
  v37 block selects index 2 and `x8` changes. The **rejected v36 block**, run through the same
  interpreter on the same inputs, never moves `x8` at all (51 inputs, 0 changes): its compare read
  the slots its own candidate had just overwritten. The defect is reproduced, not described.
* **spill-order control.** The block with the running-minimum stores moved before the compare
  (v36's order) mismatches on 6 of 60 inputs, so the check is sensitive to the ordering and not
  only to the values.
* **ADC/ADCS control.** `adc` does not update NZCV; only `adcs` does. The block therefore uses
  `adcs` in every multi-limb carry chain except the last limb. A mutant block with `adcs` replaced
  by `adc` mismatches 11 of 60 inputs, so the test is sensitive to exactly that instruction
  semantics; an interpreter that (wrongly) let `adc` write flags was the fire-2 defect.
* **operand-transfer control.** The block reads `x` from `[x19,#0]..[x19,#40]`; if the record is
  not populated the block's output is zero, so the test exercises the transfer rather than
  pre-loading its frame.
* **stack-contract control.** A block emitting an offset at 1008 in a 512-byte frame makes the
  manifest assert fail.

This is the evidence that was missing in v35 and v36: a test that executes the artifact.

## 7. Measurement, directly counted

The toy pipeline at `n = 880` counts the event `recovered == secret` once per trial and reports the
raw count **379**: `p-hat 0.4307`, Clopper-Pearson at the two-sided five per cent level
`[0.3977, 0.4642]`, computed from `(k, n)` alone. Two controls run in the same pass: with the
scale removed (`x1.0`) the same `n` gives `338/880 = 0.3841` with lower bound `0.3518`, so the
scaling is not flat; flipping one recovered-s field drops the count by exactly one (`379 -> 378`)
and moves the lower bound `0.3977 -> 0.3966`, so the count is the only channel.

The toy walker canonicalises with exact Python arithmetic and the shipped object reproduces the
same argmin (section 6), so toy and object agree by execution rather than by assertion. Toy scope
only: there is no production-size run, and none is claimed.

## 8. Success arithmetic

`S_ops = 1.065 x 1.015 x sqrt(2 N' ln(1/(1-0.3903)))` buys 6.5% more steps than the 0.3903
requirement, at `log2(1.065) = 0.0909` bits of work. The random-mapping law, evaluated on the
scaled step count, gives

    H1 = 1 - exp(-(S_ops/1.015)^2/(2N')) = 1 - exp(-1.065^2 x 0.494788) = 0.429477

Charges, each line shown: more than 8 processed matches `0.62^8 = 0.021834`; length cap `e^-20`
applied per live trail and `x 2^20` open trails `= 0.00216128`; report table `2^-18`; recovery
`2^-30`; meter `0.001047`. Success term

    0.429477 - 0.0239991 - 0.001047 = 0.40443

which is `+0.01443` above the 0.39 floor.

## 9. The two premises the v36 verdict did not accept

**(a) The 1.5% slowdown allowance.** It is an allowance, not a measurement, and it is disclosed as
one. Two things are now on record. First, a long-trail run: 120 trails of 4096 operations each
(491,520 operations, 288 escapes) on the q34 toy curve gives **2** trails with a non-fruitless
return against **1.79** expected from `1 - exp(-L^2/(2N))` (`L^2/2N = 0.0149` per trail). The
random-mapping law is consistent with long-trail behaviour; there is no structural trap. Second,
the claim is not on a knife edge at 1.5%: holding every other term fixed, H1 stays above
`0.39 + 0.026 = 0.416` for any slowdown factor up to **3.66%**. The allowance is conservative, and
the margin is quantified rather than asserted.

**(b) The cap-event population.** Two different populations matter and are now separated. The
population that can turn a cap into a *failure* is the live one: `2^20` open trails, `e^-20` per
trail, union bound `0.00216128` — no independence assumed. The population over which caps actually
*happen* is every launched trail, `S_ops/theta = 2^62.15`, giving `2^33.30` expected cap events
over the attack. Those events are a restart, not a lost trail: a capped trail's work is already
inside `S_ops` (it is counted as steps), the trail is re-launched as a fresh trail counted in the
`2^62.15`, and the repair keeps `p_in = 0`, so no match-failure probability is manufactured by the
restarts. The `e^-20` factor is therefore applied to the live population only, where a cap can
still hold the matching trail; the launched population is priced as work, not as failure mass.

## 10. What is not claimed

No production-size run. H1 is a disclosed heuristic (random-mapping law plus the 1.5% allowance)
and the toy measurement is its only empirical support; the work bound does not depend on H1. No
exact zero probability is claimed anywhere. The v33 `esc2_events` statistic stays deleted: its
counter was initialised and never incremented, so its zero was not evidence of anything.
