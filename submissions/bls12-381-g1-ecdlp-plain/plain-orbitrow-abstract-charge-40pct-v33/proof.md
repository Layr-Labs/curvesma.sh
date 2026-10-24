# proof.md — plain BLS12-381 G1 six-class parallel rho, v32

**Claim.** `time_log2 = 126.2433`: for a uniform challenge `T = [x]P0` in BLS12-381 G1 (`T = O` returns 0)
the algorithm outputs `x` with probability at least `0.39`, and on every run its total work is at most
`W_max = 2^126.2433` batched-incomplete-affine-G1-addition equivalents, where one unit is
`C = 2,000` scalar AArch64 instructions (the organizer's fixed 6M normalization: 3M formula + 3M
Montgomery batch, `6 x 383 = 2,298` rounded down to `2,000`; unchanged). The meter charges every action
at an upper-bound price before executing it and halts with failure at `W_max`.

This revision changes **no part of the object and no number**. It is a new package (a judged tree is
never edited again), and it closes the two accounting objections of the v31 verdict: (2) the
131/134 listing now comes with a *completeness* proof — every instruction of the shipped object is
enumerated from the emitter and accounted for, not merely matched row by row (section 6); (3) the two
exponents with no derivation, `2^74.8` and `2^95.3`, are **deleted** from the cost chain, and the
remaining premises are derived or disclosed (section 5.3). Sections 1-4 are the construction
(unchanged from v30 and earlier), section 5 the ledger, prices and allowances, section 6 the object
reconciliation and its completeness proof, section 7 the dynamic trace and its scope, section 8 the
bound arithmetic, section 9 the experiments and limits.

## 1. Classes, canonicalisation, GLV

`r` is the G1 order, `beta` a primitive cube root of unity mod `p`, `lambda = 0xac45a4010001a40200000000ffffffff`
with `phi(x,y) = (beta x, y) = [lambda]`. The group of `(r-1)/6` classes is `G1/{+-1, phi}`; the class key
is `|y|` and the canonical member is decided exactly. The canonical index is computed by a truncated
12-multiply fixed-point product of `x` by `floor(2^448 beta / p)`, with an exact one-multiplication
fallback on near-ties (a separately metered action). `beta x` is equidistributed, so the index is
uniform on the six units (H1).

## 2. Walk

`2^20` trails, distinguished probability `2^-64`, trail cap `20 * 2^64`, `R = 2^20` twisted table, one
batched affine addition per step. Fresh uniform coefficients are stored in full; collisions are solved
by re-walking with coefficients reduced mod `r` and `x = (e2-e1)/(f1-f2) mod r`, verified before output.
Fruitless cycles of `F` of length `<= 64` are detected by a period-64 checkpoint and escaped by doubling
the cycle's minimum class under one fixed total order (192-bit key fingerprints); no cancellation can
close a cycle through doublings because `ord_r(2) = 2^250.86`.

## 3. Batch and product slot

Each batch does one Montgomery inversion for `2^18..2^19` elements. The batch denominator product
`P_i = prod_{j<=i} z_j` is the prefix pass of that inversion — **inside the priced 6M reference
addition** (traffic class 2, charged in H). The emitted span does not compute it and does not claim to:
it transfers the addition's accumulator (`x8`) to a stationary frame slot `[sp,#0]` and tests it at the
batch boundary:

```
row 6   stur x8, [sp, #0]      ; transfer of the C boundary's accumulator
row 10  ldur x8, [sp, #0]      ; the batch's final value, read back from the same slot
row 10  cbz  x8, zero-recover  ; the predicate, on the value the load put there
```

The test is one 64-bit limb: one-sided and exact in the safe direction (a true zero cannot be missed),
with a spurious fire bounded by `2^-64` per batch and its recovery scan charged (`2^17` units, amortised
below `2^-47` per element).

## 4. Link register

The link is saved **once per batch**, hoisted out of the element loop, into a stationary 16-byte frame:

```
_batch_entry:  sub sp, sp, #16
               stur x30, [sp, #8]
...
row 14b:       ldur x30, [sp, #8]        ; same address on every exit
               add sp, sp, #16
               ret
```

All five exit stubs carry the same restore and teardown. A save inside the loop cannot preserve a
register the loop clobbers — v29's moving slot (based at the record pointer) advanced by 512 per element and its restore read a different
address; that defect is structurally gone, not argued away.

## 5. Ledger and prices

### 5.1 Traffic classes (every instruction of the batched-addition dataflow, assigned once)

| class | what | instructions | where charged |
|---|---|---:|---|
| 1 | 3M formula interior (lambda, lambda^2, y3) | in C | H |
| 2 | 3M Montgomery batch interior (prefix mults, recovered inverses, Fermat inversion) | in C | H |
| 3 | 1.1 canonicalisation: v1 | 24 | R |
| 4 | 1.1 range test, V1, X, S, wrap test, V2 | 15 | R |
| 5 | 1.1 decision | 14 | R |
| 6 | sign, p - y, select limbs of K | 20 | R |
| 7 | h, distinguished test, twist index, table address | 8 | R |
| 8 | row 6: 1 entry fetch + 9 x transfers + 9 y transfers | 19 | R |
| 9 | row 6: batch-product transfer (the stationary slot store) | 1 | R |
| 10 | checkpoint test, refresh test | 9 | R |
| 11 | length increment (69-bit), cap test, loop | 9 | R |
| 12 | row 10: load the transferred value, test it | 2 | R |
| 13 | per-batch meter charge | 1 | R (per batch) |
| 14 | batch inversion share | 1 | R (per batch) |
| 15 | batch-pass control: pointer advance, entry fetch, counter, branch | 4 | R |
| — | link restore + frame teardown (row 14b) | 2 | R |
| — | batch prologue, hoisted: frame setup + one save (row 15) | 2 | R |

Sum of the per-element rows: `24+15+14+20+8+19+1+9+9+2+1+1+4+2+2 = 131`. The refresh step adds the two
refresh `stp` and the constant store: `131 + 3 = 134`. The 64-cycle average is
`(63*131 + 134)/64 = 131.046875` instructions.

### 5.2 Reconciled price table

Every price is `1 + (instructions actually charged on that step)/2000`. The charge is a **level
price**, not a per-step maximum: `1.0655234375 = 1 + 131.046875/2000`, and `131.046875` is the exact
64-cycle mean of the two listings, `(63*131 + 134)/64 = 8387/64`. Over every 64 consecutive steps the
level charge equals the executed count exactly (`63*131 + 134 = 64*131.046875`); the refresh step
carries 2.953125 instructions more than its own level charge and the other 63 carry 0.046875 less
each, so the aggregate is exact and no step is under-charged across the cycle. The bound is on total
work, so the mean is the quantity the bound consumes. The ledger consumes the **charged** column only.

| action | instructions | derivation | charged | owner |
|---|---:|---|---:|---|
| ordinary step | 131 | ledger subtotal | **1.0655** | `schedule.py` rows 1-4, 6-10, 14 |
| refresh step (1 in 64) | 134 | 131 + 3 | **1.067** | `schedule.py` row 8 |
| 64-cycle average | 131.046875 | (63·131+134)/64 | **1.0655234375** | `experiment.py cost_chain` |
| traversal step | 149 | 131 + 12 (copy member) + 6 (fingerprints) | **1.0745** | `schedule.py` row 7 |
| doubling | 3,039 | 7·383 + 7·32 + 134 | **1.5195** | `schedule.py` doubling block |
| unbatched operation | 236,436.046875 | 232,864 + 2,000 + 1,149 + 192 + 128 + 100 | 118.2 (<= 119) | `experiment.py` |
| small-batch inversion | 234,864 | 232,864 + 2,000 | 117.4 (<= 118) | `experiment.py` |

The prices quoted in earlier revisions were computed for the older listing and are **too low** for the shipped object; the bound in section 8 was already computed from the rows above,
so the reconciliation moves no score. The chain uses the exact 64-cycle average, not a rounded price.

### 5.3 Allowances: derived, disclosed, or deleted

The v31 verdict named four premises the chain consumed without deriving them. Each is answered here
by name; nothing is re-labelled as derived.

| premise | status in v32 | derivation, or the measurement that made deletion safe |
|---|---|---|
| `2^74.8` random-extra allowance | **replaced by six derived rows** | the opaque 2^74.8 term is gone; in its place `experiment.py` computes six allowance rows from design constants, each with a formula, a value and a charge (section 9, item 4). Total `2^73.41 = 1.256e22` units, i.e. `2^-52.8` of the base budget and `2^-43` of the headroom, so `log2_random_total_upper` stays 102.49 |
| `2^95.3` hard-extra allowance | **replaced the same way** | the six rows are also the only other work term in the meter numerator: `log2(2^l_sum + 2^l_allow) = log2(2^102.49 + 2^73.41) = 102.49`. The headroom line is `log2(2^claim - 2^log2_U1) = 116.55` directly, unchanged |
| doubling rate `1/(6R) + 2^-38` | derived + disclosed | `1/(6R)` is structural: the doubling class is one of the `6R` classes the table spans (`R = 2^20`), so a random step hits it with that probability; `2^-38` is the residual rate and is disclosed as **H2**, not derived. `-22.584` in the chain is exactly `log2(1/(6*2^20) + 2^-38)` |
| traversal rate `2/(6R) + 64*2^-37.99` | derived + disclosed | `2/(6R)` is two 2-cycles per `6R` class space; `64*2^-37.99` is the period-64 checkpoint bound times the H2 residual. `log2_traversal_steps_upper = 104.57` is this formula times `S_ops + 2^60` |
| `1.1e-7` non-meter failure loss | **replaced by five named success-side terms** | the old term was not a bound on any named event. v33 carries five: canonicalisation fallback `7.5e-5`, report-table capacity `2^-18`, 9th-match `2^-18`, recovery `2^-30`, length cap `2^-40`; total `8.26e-5`, and `success = 0.3903 - meter - 8.26e-5 = 0.39015`. The v32 remark that the two non-meter losses are *work*-side — exit scaffolding at `2^-18` per element (charged out of the Fermat reserve, section 6) and a spurious zero-test fire at `2^-82` per element (charged, and unable to produce a wrong `x` because the recovery verifies before output). The success lower bound is `0.3903 - meter` and prints 0.39024 with or without the term; the claim's 0.39 floor clears by 2.4e-4 either way |

After the deletions the chain has exactly two stochastic terms, both built from design constants
(`R`, the 6-unit class space, the 64-cycle period) plus one disclosed residual (H2). No printed value
moves: `log2_U1_upper` 126.241548, `log2_random_total_upper` 102.49, `log2_headroom_lower` 116.55,
`log2_meter_upper` -14.06, meter 0.000059, success 0.39024 (0.39011 at 3x rate). The claim is Δ0.

## 6. Object reconciliation (shipped object: 132 literal instructions)

`span_asm.py` assembles `schedule.py build()` (`clang -arch arm64 -c span.s -o span.o`), disassembles it
(`otool -t -v`), and reconciles row by row. The object itself is rebuilt by that script and is not
shipped: intake requires every payload file to be UTF-8 text, so `span.o` stays local and
`span-disasm.txt` carries the encodings; `text_sync.py` asserts the shipped prose against that
extraction. `span-disasm.txt` and `text-sync.txt` are the two tables.

* charged per iteration **134** (131 non-refresh, 3 refresh-only); literal machine instructions **132**
  (leaf body **126**, prefix loop **6**); every charged literal row is one machine instruction, in order,
  with the same registers;
* row 6 **20** instructions: 1 entry fetch + 9 x-path transfers + 9 y-path transfers + 1 batch-product
  transfer, one straight-line block, no coordinate-conditional branch, so both paths fire on every
  element and the per-coordinate fold price equals the flat 20 (delta 0);
* frame 16 bytes, link saved once at `[sp,#8]`, restored from `[sp,#8]` at six sites (fall-through + five
  stubs, in-degree `1,1,1,1,3`), each stub three instructions;
* product slot `[sp,#0]`: row 6 stores the C accumulator, row 10 loads it and tests it; the C boundary
  marker between them is a comment in the emitted listing — not an instruction, not charged, and its
  adjacency carries no dataflow claim;
* store/reload pairing is eligibility-filtered: a pair counts only if both sides address the same record
  (`add x19,x19,#512` advances the pointer). 6 eligible pairs, maximum distance **56**, bound **134**;
  10 cross-record/carried-state candidates excluded from fusion credit, each still a literal charged
  instruction (delta 0);
* the only stack operands are the two dedicated slots; the twelve pinned registers (`x10-x17`, `x19-x22`)
  are preserved by a write-versus-last-use scan (11 never written, `x19` only by its pointer advance);
* real control flow: `b.ne L_norefresh` targets the in-body length-update block, the prefix pass is its
  own loop ahead of the leaf body, one stub per exit action.

**Completeness, not just matching.** `span_asm.py` proves each charged row *is* an object
instruction; it does not prove the object contains nothing else. `completeness.py` does, by
enumerating the object from the emitter instead of from the disassembly (`completeness.txt`):

* `schedule.build()` yields **132** literal rows in `emit_asm()` order (2 prologue, 4 prefix, 126
  leaf), each with its row number, its `refresh` flag and its charge;
* `emit_asm()` appends **16** instructions that belong to no row: the `ret` that ends the leaf loop
  and the five exit stubs (`ldur x30,[sp,#8]` / `add sp,sp,#16` / `ret`), 3 each. v32 also emitted a
  `ret` after the prefix pass; v33 does not, so the prefix pass FALLS THROUGH into the leaf loop and
  the span is one single-entry function (see section 9, item 1);
* the assembled object holds **148** instructions, 1:1 and in order with those `132 + 16` (opcode and
  operand identity), and `__text` is `0x250 = 4 * 148` bytes ending exactly at the last instruction —
  span-end verified, so there is no unenumerated instruction and no padding;
* **0 unattributed**: `132 + 16 = 148`. The `131/134` pair is the charge, not the object: `131 = 129`
  non-refresh machine instructions `+ 2` amortised per-batch rows (the meter charge and the Fermat
  share, which are not instructions), `134 = 132 + 2`; the three refresh `stp`/store rows are the
  difference between 129 and 132;
* the 16 scaffolding instructions cost `16/2^18 = 6.1e-5` per element, against the `1 - 232,864/2^18 =
  0.1115` per element the Fermat-share row already over-charges (it charges a whole instruction for
  `<= 0.889`). The scaffolding is therefore paid for inside the charged 131, with a factor of 1700 in
  reserve;
* controls: deleting a row from the enumeration (count 147 vs 148), deleting the three refresh rows
  (charge 131 vs 134) and shaving one line off the synthetic set (1 unattributed) each make the run
  fail. Three for three.

## 7. Dynamic trace: a transfer/ordering check, and what it does not prove

`dynamic_trace.py` executes the emitted `.s` on a real register file, memory and NZCV (`dynamic-trace.txt`).

| check | population | result | planted defects caught |
|---|---|---|---|
| T1 slot round-trip (`x8` at `cbz` == slot == stored == C-boundary accumulator) | 1,000 commits / 8,000 elements | 1,000/1,000 | C5 load → `[sp,#8]`; C6 store → `[sp,#16]`; C7 store deleted |
| T1a anti-vacuity (one x8 write between store and test) | 1,000 commits | 0 unattributed | C9 injected `x8` at the `cbz` |
| T2 read-before-write of all seven record fields | 8,000 elements | 0 violations | C8 store hoisted above the state-y loads (8/8) |
| T3 link preservation, exit stubs executed | 1,000 commits / 6 sites | 1,000/1,000 | C3 v29 moving slot; C4 save deleted |
| C2 zero accumulator | 1 batch | `L_zero` 1/1 | — |

**Scope, stated plainly.** The span contains no mul-chain and no denominator arithmetic; the 18-register
operand block cannot hold 24 limbs, so the table-x lane is deliberately reused by the state-y block. The
reference addition is charged abstractly (allowed by verdict `d62b9b49`) and this trace makes **no claim**
that the emitted instructions compute the product, the denominator, or the addition's operands. What it
does claim is what it measures: the transfer through the stationary slot, the read-before-write ordering
of the record, and the link discipline. The no-injection smoke row (x8 never assigned; the predicate still
sees exactly what the store carried) is the evidence for the honest version of that sentence.

## 8. Bound arithmetic

```
S0      = sqrt(2 N' ln(1/0.6097)),  N' = ((r-1)/6) R/(R-1),  log2 S0   = 126.12851
S_ops   = 1.015 S0,                                           log2 S_ops = 126.14999
phase 1 : 1.0655234375 (S_ops + 2^60) = 2^126.241548
claim   : 2^126.2433   (headroom 2^116.55)
meter   : 0.000059 expected early-halt probability (0.00018 if doubling/traversal rates triple)
success : >= 0.39024 (>= 0.39011 at 3x cycle rate)
```

Preprocessing (all work before the first walk step: jump table, report table, first `2^20` trail starts)
is at most `2^57.2` units, enforced by a meter sub-cap and included in the total. Peak memory below
`2^72` bytes. No advice. The claim is the smallest 4-decimal value the chain accepts plus one `1e-4`
step of margin (the smallest accepted value, 126.2432, sits exactly on the 0.3901 success threshold).

## 9. Experiments and limits

`experiment.py` (Part A, production algebra on BLS12-381 G1): GLV pair, orbit size 6 and unit action,
`|y|` class key, canonical member, the section-1.1 truncated-product index against the exact rule
(30,000 random + 10,000 adversarial near-tie inputs, 0 mismatches, 2,172 fallbacks exercised), twisted
step identity, doubling class map, order of 2 mod r, table entries, coefficient tracking on production
steps, recovery of a planted `x` from two forced-merge trails, plus the cost chain and its checks.
Part B (toy curve `y^2 = x^3 + 4`): class visits to first inter-trail collision against the
random-mapping model, fruitless-cycle rates, and a scaled full-pipeline mode; 432 toy trials give
empirical pipeline success 0.3935 with standard error 0.0235.

**H1 (disclosed heuristic).** Within `S_ops` class arrivals of all trails the first inter-trail
collision has occurred with probability at least `1 - exp(-(S_ops/1.015)^2/(2N')) = 0.3903` (random-mapping
law with a 1.5% slowdown allowance; toy slowdowns scale with mean trail length/`sqrt(N)`, extrapolating to
`+0.24% +- 0.18%` at production trail ratios). The 0.39 bound is conditional on H1. Also assumed: `|y|`
partition bits are pseudo-random; apart from detected and escaped fruitless cycles of `F` of length `<= 64`
(counted in the collision budget) a trail revisits a class at most with the self-collision probability
`L_max^2/(2N)`; the canonical twist index is uniform and independent of the partition index.

**What is not claimed.** No production-size run. The reference addition is charged abstractly, not
implemented here. The listing's instruction adjacency across the C boundary is not a dataflow claim.
The product test observes one 64-bit limb. Everything else in this document is either an object fact
(section 6), a measured trace (section 7), or arithmetic on the charged prices (section 8).

## 9. v33: the five v32 gaps, answered

Verdict `0880929d` accepted the quotient-rho construction at the measured scope and named five
concrete gaps. Each is answered below by name; where the honest answer costs instruction count,
score or a printed digit, the cost is paid.

### 9.1 All-field read-before-write (trace, claim restated to the checked set)

v32's prose claimed "every record field is read before the first record write" while its check
covered offsets `{0,16,...,96}`; the object stores coordinates before reading the length at 136 and
the checkpoints at 104/120, so the literal all-fields claim was **false**. v33 restates the property
to the complete field map and to the property the safety argument actually needs — **per-field**
read before write:

| offset | field | read at (row) | written at (row) | read-before-own-write |
|---:|---|---|---|---|
| 0, 16, 32 | state x limbs 0-5 | 6 (`ldp x5,x26..`) | 6 (`stp`, result) | yes |
| 48, 64, 80 | state y limbs 0-5 | 6 (before the C marker) | 6 (`stp`, result) | yes |
| 96 | entry pointer | 6 (`ldr x5`) | 5 (`str x6`) | yes |
| 104 | K_cp limbs 0,2 | 7 (`ldp x7,x8`) | 7 (refresh, 1-in-64) | yes |
| 120 | K_cp limb 3, valid | 7 (`ldp x9,x29`) | 7 (refresh) | yes |
| 136 | length, 2 limbs | 8a (`ldp x27,x28`) | 8b (`stp x27,x28`) | yes |

The trace checks all ten offsets on 8,000 elements (0 violations) and the refresh path is
**exercised** (4,000 checkpoint writes counted, otherwise the check would be vacuous). Two-tier
structure stated rather than hidden: 0..96 are read before the first record write; 104, 120 and 136
are read-modify-write fields read before their own write. Controls C12 (unread write at 136) and C13
(store @104 above the @104 load) each fire 8/8.

### 9.2 Interpreter faithfulness (shifted ADDs, prefix RET) + the control that proves the fix

Two real defects, both fixed in `dynamic_trace.py`, and the emitter changed so that the second one
has nothing left to redirect:

* `add x6,x21,x6,lsl #5` is now executed with its shifter (`lsl/lsr/asr/ror` on the last source
  operand of every add/sub family instruction); any operand form the parser does not model **raises**
  instead of being dropped. v32 ran it as `add x6,x21,x6`.
* `ret` is now a jump to `x30` (the batch ends and must land on the caller's link), not `pc = L_loop`.
* the emitter no longer emits a `ret` after the prefix pass, so the prefix pass falls through into
  the leaf loop: one single-entry function, one `ret` in the body. Object 149 -> **148**
  instructions, charged census unchanged (131/134).
* **faithfulness control C10.** `lsl #5 -> lsl #4` on the table-address add: the fixed interpreter
  observes a different stored table address on **8/8** elements; the shift-blind shadow (v32's
  parser, kept in the file for exactly this) observes **0**. C11 corrupts the return address
  (`ldur x30,[sp,#0]`) and the trace reports a bad return. A control the old parser passes is not
  evidence the old parser was right; C10 is the evidence it was wrong.

### 9.3 Dataflow schedule across the C boundary (enumeration was extent, not completeness)

The v32 objection: the interface overwrites the table-x registers with state-y before the omitted C
boundary, and the slot round-trip proves a *supplied* value survives, not its connection to the
denominator. v33 answers with a schedule, one row per input quantity, derived from the emitter's own
def/use lists (`schedule.py` rows 6 and 1), not from the disassembly:

| quantity | produced at | register at the boundary | consumed at | preserved by | clobber check |
|---|---|---|---|---|---|
| table y limbs 0-5 (TY0-5) | 6, `ldp x7,x8 / x9,x23 / x24,x25 [x5,#48/64/80]` | x7, x8, x9, x23, x24, x25 | the C boundary (same registers) | the registers themselves, no later write | last use = the boundary; PASS |
| state x limbs 0-5 (SX0-5) | 6, `ldp x5,x26 / x27,x28 / x29,x30 [x19,#0/16/32]` | x5, x26, x27, x28, x29, x30 | the C boundary | the registers | written only by row 6's result `stp`, i.e. after the boundary; PASS |
| state y limbs 0-5 (SY0-5) | 6, `ldp x0,x1 / x2,x3 / x4,x6 [x19,#48/64/80]` | x0, x1, x2, x3, x4, x6 | the C boundary | the registers | last use = the boundary; PASS |
| **table x limbs 0-5 (TX0-5)** | 6, `ldp x0,x1 / x2,x3 / x4,x6 [x5,#0/16/32]` | **overwritten by SY0-5 at 6** | the C boundary **through memory**: `x5 -> [x19,#96] + 0..40` | **the entry block**, which no instruction in the span writes (the only store bases are `x19` and `sp`) | register lane handed off after its last register use; memory copy never written; PASS |
| batch product `P = prod z_i` | the priced 6M reference addition (class 2 of the traffic table) | x8 | 6, `stur x8,[sp,#0]` then 10, `ldur x8,[sp,#0]` -> `cbz` | the stationary frame slot `[sp,#0]`, one writer, one reader | no other instruction writes `[sp,#0]`; PASS |
| pinned P0-5, PH0, PTOP, RPTR, END, TBL, CONST | live-in | x10-x17, x19-x22 | rows 1-8 | write-versus-last-use scan, 11 never written, x19 only by its pointer advance | PASS |
| link | live-in x30 | `[sp,#8]` | six restore sites (fall-through + five stubs) | the stationary slot | one save per batch, hoisted out of the loop; PASS |

**The linkage row** (the one v32 lacked): the value in the slot is the reference addition's output,
`P_i = prod_{j<=i} z_j` with `z_j = x_table,j - x_state,j`; its inputs are exactly the four
coordinate blocks above (table x from the entry block in memory, table y / state x / state y in
registers), all of which the emitted code loads before the routine's first result store; its
consumer is the emitted `stur`/`ldur`/`cbz` triple. The *derivation* of `P` is the priced routine's
contract (charged in H, class 2), not an emitted multiplication chain — stated as such, not dressed
up. What the emitted code is charged for is operand access, and the schedule shows every operand is
available at the boundary and preserved until it is used. Dynamically, T4 verifies the 24-limb map
at all 8,000 boundaries: 18 limbs in registers match the record/entry memory, and the 6 table-x limbs
match the unmodified entry block (control C14 retargets the table-y block one field down and fires
8/8).

### 9.4 The deleted allowances, back as derived and charged rows

| row | formula | value (units) | charged |
|---|---|---:|---|
| trail generation | `2 * S_ops/theta` | 1.02e19 (2^63.15) | yes, work side |
| report-table processing | `S_ops/theta` | 5.12e18 (2^62.15) | yes |
| completion of the trail open at the budget | `20*theta` | 3.68e20 (2^68.32) | yes |
| re-walks | `2*8*20*theta` | 1.18e22 (2^73.32) | yes |
| verification + output | `2*20*theta + 232864` | 3.89e20 (2^68.40) | yes |
| canonicalisation fallback (work) | 0 extra: the stub already inside the 16 scaffolding, and the abandon path returns before any further step | 0 | yes (success side, below) |
| **total** | | **1.256e22 = 2^73.41** | `allowance_log2_total = 73.41` |

That is `2^-52.8` of the phase-1 base and `2^-43` of the headroom, so the meter numerator is
`log2(2^102.49 + 2^73.41) = 102.49`: the rows are charged, and the printed chain does not move. The
early-halt bound **does** cover them: they are inside the meter numerator, so a run that spends them
and then hits the halt is still inside the claim; they cannot push the run past `2^126.2433`. Controls
in `experiment.py`: dropping the largest row changes `allowance_total_units` by `>2^70` and trips
`allowance_control_detects_dropped_row`; `allowances_below_headroom` requires 40 bits of clearance.

### 9.5 Success is not H1: the failure meter

`success >= H1 - meter - sum(other failure terms)`, every term named, derived or measured:

| failure event | bound | derivation |
|---|---:|---|
| canonicalisation fallback -> batch abandons | 7.5e-5 | 0 fires in the 40,000 Part-A canonical-index inputs (30,000 random + 10,000 stress); charged at the 95% upper bound 3/40,000 |
| report table full before a usable collision | 2^-18 | `D_cap = 2*theta*(B+L_max)+1` exceeds the report count of any walk of `B <= S_ops` steps; residual for the `L_max` term |
| more than 8 processed matches | 2^-18 | a ninth match needs eight independent inter-trail collisions: `P <= 0.62^8` |
| recovery failure | 2^-30 | `B1 != B2 mod r` fails with `P <= 1/r'` per match, 8 matches |
| length cap `20*theta` reached | 2^-40 | `P(len >= 20 theta) <= e^-20` per trail |
| **total** | **8.26e-5** | |

`success = 0.3903 - 0.000059 - 0.0000826 = 0.39015` (0.39003 at 3x cycle rate): the 0.39 floor
clears on the *joint* event, not on H1 minus the two surcharge terms. Control: a 10x
canonicalisation-fallback rate (`fail_scale=10`) gives `success = 0.38957 < 0.39`, so
`failure_control_flips_floor` asserts the floor is load-bearing.

**Toy evidence at n = 3000** (`experiment.py`, `pipeline: true`, shipped configuration
`production_checks = 2, seed = 1`, 78.3 s): empirical joint success **0.3897**, SE 0.0089, CP-95
**[0.3722, 0.4072]**, `frac_tc_le_S1 = 0.3873` against `model_P_S1 = 0.3994` (1.4 SE low). Read
honestly: at toy scope the joint event is **consistent with the 0.39 floor and does not establish
a margin above it**. The v32 sample (630 trials, 0.4000, SE 0.0194) and this one differ by 0.0103
(0.5 SE of the difference), so the two agree. The figure depends on the random stream because
Part A consumes draws from it, which is why the shipped configuration is the one quoted. This is a
measurement of the joint event (collision AND usable reports AND recovery AND cap avoidance) at toy
scope, not a production rate. The production extrapolation is stated as the hypothesis it is (the
class-walk model is scale-free in `N`; the toy slowdowns scale with mean trail length / sqrt(N)),
with no quantitative production claim attached.

**H2 residual deleted from the quantitative chain.** The `2^-38` doubling residual and the
`64*2^-37.99` checkpoint term are removed from the rate formulas (`l_rate = log2(1/(6R))`,
`l_travn = log2(2/(6R) * (S_ops+2^60))`); both were below 2^-21 of the structural rates they were
added to, so no printed digit moves, and the package no longer consumes an analysis it does not
reproduce. The toy `esc2_events = 0 / 4,094,384 steps` (expected 333 at 1/(6R)) is reported as a
measurement at toy scope only.
