# proof.md — plain BLS12-381 G1 six-class parallel rho, v31

**Claim.** `time_log2 = 126.2433`: for a uniform challenge `T = [x]P0` in BLS12-381 G1 (`T = O` returns 0)
the algorithm outputs `x` with probability at least `0.39`, and on every run its total work is at most
`W_max = 2^126.2433` batched-incomplete-affine-G1-addition equivalents, where one unit is
`C = 2,000` scalar AArch64 instructions (the organizer's fixed 6M normalization: 3M formula + 3M
Montgomery batch, `6 x 383 = 2,298` rounded down to `2,000`; unchanged). The meter charges every action
at an upper-bound price before executing it and halts with failure at `W_max`.

This revision changes **no part of the object**: it changes the prose to match the object, replaces a
circular dynamic check with a transfer/ordering check that can fail, and reconciles the price rows.
Sections 1-4 are the construction (unchanged from v30 and earlier), section 5 is the ledger and the
reconciled price table, section 6 the object reconciliation, section 7 the dynamic trace and its scope,
section 8 the bound arithmetic, section 9 the experiments and limits.

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

Every price is `1 + (instructions actually charged on that step)/2000`; the charge is
`max(static listing, dynamic count)`, which coincide by construction (section 6). The ledger consumes
the **charged** column only.

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

## 6. Object reconciliation (shipped object: 132 literal instructions)

`span_asm.py` assembles `schedule.py build()` (`clang -arch arm64 -c span.s -o span.o`), disassembles it
(`otool -t -v`), and reconciles row by row; `text_sync.py` asserts the shipped prose against that
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
