# plain v29 — folds deleted (no dataflow claim), batch-product row 10, x30 link save/restore

Track: `bls12-381-g1-ecdlp-plain`. Claim **126.2405** (v28: 126.2392). The rise is a correction:
the two fixes that cost instructions are those the v28 verdict demanded, and the third fix removes
a claim rather than adding code.

## The three v28 verdict flaws (`28abe582`) and what changed

| verdict flaw | v28 state | v29 fix | evidence |
|---|---|---|---|
| the folds drop values and the emitter deletes them, so the audited sequence is not the executed one: fold 1 omitted state-y and the Montgomery recovery state, fold 2 omitted old x and the saved lambda | two cost-0 `[C]` rows; the object was loads and stores only, so 19 matching transfers witnessed nothing | **the fold rows and the claim are deleted.** Row 6 is 19 charged transfers = operand access for the fixed-price reference addition (6M, priced in H), with **no zero-extra-cost interface and no interior dataflow claim**. What the object proves is ordering, now correct: all eight input loads precede the first result store (v28 stored new x before reading state y). Audits run on the **emitted object** (`clang -c` -> `otool -t -v`), and the evidence carries the listing hash | `schedule.py`, `span_asm.py`, `span-disasm.txt` |
| row 10 tests `x1` after the leaf loop, where it holds the last record's key limb K1 — not the batch product, and the row-5 report branch already exited on it | `cbz x1,zero-recover` | **`cbz x8,zero-recover`**: the batch denominator accumulator Pi = prod z_i, the quantity the shared inversion consumes and the batch commit gates on. Planted control: K1 != 0 with Pi = 0 -> the old test passes, the new one fires | `schedule.py`, `audit.py` control |
| exit stubs `ret` on x30, but rows 6 and 4 clobber x30 (state-x limb 5, `sbcs` temp) | `ret` with no save | **`stur x30,[x19,#-8]` as the first leaf instruction; `ldur x30,[x19,#-8]` before every `ret`** (fall-through and all five stubs). Pinned-preservation table extended with an x30 row; planted-clobber control refires | `schedule.py` pinned table, `span-disasm.txt` |

## Price

`stur`+`ldur` cost +2 charged instructions per element (the five stub restores are exit-only, at
most one fire per pass, amortised < 2^-18). Ledger **127** non-refresh / **130** refresh step,
step price `1 + 127.046875/2000 = 1.0635234375`, traversal `1.0725` (127 + 18 = 145), doubling
`1.5175` (3,035 instructions), phase-1 bound `1.063523 (S_ops + 2^60) < 2^126.238842`, headroom
`2^126.2405 - 2^126.238842 - 2^95.3 > 2^116.47`, early-halt `< 2^-13.98 = 0.000062`, success
`0.3903 - 0.000062 - 1.1e-7 = 0.39024 > 0.39` (`0.39012` at 3x the cycle rate). Row 10's predicate
is a 1-for-1 swap and costs nothing.

## Board posture

Our own `5ecbf15` = **126.2364** (valid) stands and is not re-used; `67154fda` = **126.2337** was
in review state at fire time. This package claims 126.2405 — worse than both, deliberately: v28's
ledger omitted the link save/restore and mis-priced row 10, so 126.2392 was not a bound.

## Object reconciliation (re-run on the v29 object, not inherited)

| check | value |
|---|---|
| charged per iteration | 130 (127 non-refresh, 3 refresh-only) |
| literal instructions in `.o` | 128 (leaf 124 + prefix 4) |
| row-6 region | 19 instructions: 1 entry + 9 table/state loads + 6 state-y loads, then 6 result stores |
| row-by-row identity mismatches | 0 (EXACT) |
| back-edges | 1 leaf, 1 prefix loop |
| exit stubs | 5 distinct targets, in-degree 1,1,1,1,3, each `ldur x30; ret` |
| eligible fusion pairs | 6, max distance 56, bound 130 |
| controls | 8/8 schedule, 4/4 planted flips in audit |

## Controls

`schedule.py` runs eight planted violations, each of which must make the verifier fail: A
clobber-while-live, B missing preservation, C wrong output register, D census flip, E pinned
clobber, F fusion-bound flip, G a record operand made `sp`-relative, H an exit edge aliased into
the body. All eight fire. `audit.py` adds four planted flips (inequality reversal, equation
perturbation, a dropped listing row, a dropped fold transfer), all detected.

## Objective

Recover x from T = [x]P0 in BLS12-381 G1 with no auxiliary powers, priced as H + R/C where H
counts batched incomplete affine G1 additions (C = 2,000 scalar AArch64 instructions, the
organizer's fixed 6M convention) and R counts every other scalar instruction. The claim is an
unconditional work bound enforced by a global meter, plus a success probability >= 0.39 driven by
heuristic H1 (toy-supported random-mapping law with a 1.5% slowdown allowance).

## Environment and setup

macOS (arm64), Python 3.14 stdlib only, Apple clang for the object check:
`clang -arch arm64 -c span.s -o span.o` then `otool -t -v span.o`. No third-party packages, no
network access, no GPU. The algebra and toy experiments run on the production curve for part A and
on toy curves y^2 = x^3 + 4 over small prime fields for part B.

## Baseline and provenance

Lane history: v22 126.2364 (our own `5ecbf15`, the best *valid* published position), v26 126.2385,
v27 126.2385 (three listing repairs at zero cost), v28 126.2392 (prefix loop charged honestly,
rejected on the fold-witness, row-10 and x30 defects), this v29 126.2405 (+0.0013 for the link
save/restore, with the fold claim withdrawn). The v28 tree is superseded here; nothing in this
package re-uses the v28 folded-row claim.

## Exact commands

```
python3 schedule.py            # contract: liveness + occupancy, census 127/130, pinned table
                               # incl. x30, per-coordinate row 6, filtered pairs, 8 controls
python3 span_asm.py            # emit -> clang -> otool -> reconcile 1:1, enumerate branch edges
python3 audit.py               # arithmetic/direction/literal checks + planted flips
python3 experiment.py
```

## Measured results

* `schedule.py`: peak 30 occupied registers (0 free of 30 usable), 127 charged non-refresh / 130
  refresh, rows 24/15/14/20/8/19/12/9/1/1/1/4/2, 12/12 pinned preserved plus x30 save/restore,
  6 eligible pairs max distance 56, 10 excluded candidates, 8/8 controls fire.
* `span_asm.py`: 128 literal instructions (124 leaf + 4 prefix), 0 identity mismatches, 1 leaf
  back-edge, 1 prefix back-edge, 5 exit stubs (in-degree 1,1,1,1,3), `b.ne norefresh` -> the
  length-update block.
* `audit.py`: 0 findings, all planted flips detected.
* `experiment.py`: part A all true (30,000 random + 10,000 adversarial canonical-index inputs, 0
  mismatches, 0 fallbacks on random); toy pipeline at the claimed 0.39 floor.
* Unsuccessful runs: the v27 tree (rejected, `invalid`); v28 (rejected, 28abe582: fold rows that
  the emitter deletes, row 10 testing K1, `ret` on a clobbered x30).

## Limitations

H1 is a heuristic supported only at toy scale; instruction counts are analytic plus one
object-verified listing; no production-size run is claimed, and rejection does not require one.
Row 6 prices the operand access of the reference addition and claims nothing about its interior.

## Next steps

Keep the listing-level lane closed: the remaining objections are about the priced implementation,
and v29 answers them against the emitted object rather than the source listing.
