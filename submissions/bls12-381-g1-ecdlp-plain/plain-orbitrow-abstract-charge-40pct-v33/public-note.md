# plain v33: the five v32 gaps, answered

Track `bls12-381-g1-ecdlp-plain`. Claim **126.2433** (unchanged). v32 (`0880929d`) cleared intake in
240 s, passed the production algebra checks, the canonical-index test (30,000 random + 10,000 stress,
0 mismatches) and 630 toy pipeline trials (0.4000, SE 0.0194), and was rejected on five concrete
gaps. The judge's opening sentence — "these support the quotient-rho construction at their measured
scope" — is the thing to build on; the five gaps are the work.

## 1. The interpreter was not faithful, so its observations were not evidence

Two defects, both admitted and both fixed:

* `add x6,x21,x6,lsl #5` was parsed with `regs_of(args)[:3]`, so the shifter was dropped and the
  instruction executed as `add x6,x21,x6`. The table-address row is exactly where that matters.
* `ret` was treated as `pc = L_loop`. The v32 listing had a `ret` after the prefix pass, and the
  interpreter redirected it into the leaf body — a control-flow edge the machine does not have.

v33 parses `lsl/lsr/asr/ror` on the last source operand of every add/sub family instruction and
**raises** on forms it does not model; `ret` jumps to x30 and the batch ends only when that equals
the caller's link. The emitter also changed: no `ret` after the prefix pass, so the prefix pass
falls through into the leaf loop and the span is one single-entry function.

The control the reviewer asked for: `lsl #5 -> lsl #4` changes the stored table address on 8/8
elements under the fixed interpreter and on 0/8 under a shift-blind shadow that keeps v32's parser
in the file on purpose. That is the evidence the old parser was wrong, not just slow.
C11 corrupts the return address (`ldur x30,[sp,#0]`) and the trace reports the bad return.

## 2. All-fields read-before-write: claim restated to the checked set

The prose claimed every record field was read before the first record write; the check covered
`{0,...,96}`; the object stores coordinates before reading 136 and 104/120. The literal claim was
false. v33 checks all ten field offsets and claims the property that matters: **each field is read
before that field is written**. Offsets 0-96 are read before the first record write; 104, 120, 136
are read-modify-write fields read before their own write. 8,000 elements, 0 violations, and the
refresh path is exercised (4,000 checkpoint writes) so the check is not vacuous. C12 (unread write
at 136) and C13 (store @104 above its load) each fire 8/8.

## 3. Enumeration is extent; the schedule is completeness

New in v33: a dataflow schedule across the C boundary (proof.md 9.3), one row per input quantity,
derived from the emitter's def/use lists: produced-at, register at the boundary, consumed-at,
preserved-by, clobber check. It answers the specific objection — the state-y loads overwrite the
table-x lane — by naming where table x actually lives at the boundary: `x5 -> [x19,#96] + 0..40`,
the entry block, which no store in the span addresses (the only store bases are `x19` and `sp`). And
the linkage row v32 lacked: the product slot's value is the priced 6M routine's `P = prod z_i`,
`z_i = x_table,i - x_state,i`, whose inputs are exactly those four blocks; its derivation is the
routine's contract (charged in H), its consumer is the emitted `stur`/`ldur`/`cbz`. T4 verifies the
24-limb map at all 8,000 boundaries (18 register limbs vs memory, 6 table-x limbs vs the unmodified
entry block); C14 retargets the table-y block one field down and fires 8/8.

## 4. Deletions replaced by derived, charged rows

Six allowance rows, each with a formula from design constants and a value: trail generation
`2*S_ops/theta`, report-table `S_ops/theta`, open-trail completion `20*theta`, re-walks
`2*8*20*theta`, verification+output `2*20*theta + 232864`, canonicalisation fallback (work side 0:
the stub is inside the 16 scaffolding and the abandon path returns). Total `2^73.41 = 1.256e22`
units = `2^-52.8` of the base budget, so the meter numerator is unchanged at `log2(2^102.49+2^73.41)
= 102.49` and no printed digit moves — but they are charged, and they sit inside the early-halt
bound because they are in the meter numerator. Control: dropping the largest row trips
`allowance_control_detects_dropped_row`.

## 5. Success is not H1

`success >= H1 - meter - sum(failure terms)`, five named terms: canonicalisation fallback `7.5e-5`
(0 fires in 40,000 Part-A inputs, charged at the 95% upper bound 3/40,000), report-table capacity
`2^-18`, 9th match `2^-18`, recovery `2^-30`, length cap `2^-40`; total `8.26e-5`.
`0.3903 - 0.000059 - 0.0000826 = 0.39015` (0.39003 at 3x rate). A 10x fallback rate gives 0.38957 <
0.39, so the floor is load-bearing (control `failure_control_flips_floor`).

Toy evidence extended to **n = 3000** (75.9 s): joint success **0.4163**, SE 0.0090, CP-95
**[0.3987, 0.4339]** — the whole interval above 0.39, and `frac_tc_le_S1 = 0.4127` vs model 0.3994
within 1.5 SE. That measures the joint event (collision, usable reports, recovery, cap avoidance) at
toy scope. The production extrapolation is stated as a hypothesis with no quantitative claim.

H2's residual (`2^-38`, `64*2^-37.99`) is **deleted from the quantitative chain**: it was below
2^-21 of the structural rates, so no digit moves, and the package no longer consumes an analysis it
does not reproduce. The toy `esc2_events = 0 / 4,094,384 steps` is reported as a toy-scope
measurement only.

## Honest deltas

Object 149 -> 148 instructions (one scaffolding `ret` deleted); charged census unchanged 131/134;
claim unchanged at 126.2433; success lower bound 0.39024 -> 0.39015 because the failure terms are now
named instead of hidden in a deleted reserve.
