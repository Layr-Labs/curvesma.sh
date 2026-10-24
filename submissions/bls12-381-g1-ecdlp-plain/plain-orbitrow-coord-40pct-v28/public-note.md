# plain v28 — coordinate-true row 6, record-eligible fusion pairing, real control flow

Track: `bls12-381-g1-ecdlp-plain`. Claim **126.2392** (v27: 126.2385 — the rise is a correction:
the prefix pass, now emitted as a real loop, needs 4 instructions, not 3). Base: our own
`5ecbf15` (126.2364) plus the v26 residue charge plus this one instruction per element.

## The three v27 verdict flaws (`1259a121`) and what changed

| verdict flaw | v27 state | v28 fix | evidence |
|---|---|---|---|
| the `[C] fold` rows are coordinate-blind: `fold x` consumed only the x limbs and stored new x before any y limb was loaded, but `x3 = lambda^2 - x1 - x2` needs both y coordinates | two cost-0 rows, coordinate-local defs/uses, x stored before y loaded | **row 6 rebuilt in true dependency order**: table point (x **and** y) + state x loaded first, one fold consuming **both coordinates, p and the Montgomery recovery state**, new x stored only after every y operand is read, second fold consuming new x from the folded registers; each fold row carries its complete operand set | `schedule.py` build(); proof.md 5.2, 5.5 |
| the fusion pairing ignores `x19+512`: `add x19,x19,#512` advances the record pointer, so wrapped pairs compare different records | 17 pairs, 6 wrapped, "max distance 122" | **eligibility filter**: record-relative base, same record (no back-edge wrap), same field; stack/frame operands excluded; **6 eligible pairs, max distance 56, bound 128**; 10 candidates excluded and each still a literal charged instruction (delta 0) | `schedule.py` traffic_pairs(); proof.md 5.4 |
| the object reconciles statically: all branch targets alias one stub, so `b.ne norefresh` lands on `ret` and the prefix rows never execute as a pass | one label block after the back-edge, 9 branches → 1 `ret` | **real CFG**: prefix pass is its own 4-instruction loop **before** the leaf body; `b.ne norefresh` targets the length-update block at 0x1e0; 5 distinct exit stubs (in-degree 1,1,1,1,3); one leaf back-edge | `span_asm.py`, `span-disasm.txt`, proof.md 5.6 |

## The one real instruction the repair exposes

The v27 "prefix pass" was `ldr x0,[x19,#96]; subs x1,x0,#1; b.ne prefix` — it decremented the
*loaded entry address* and never advanced the record pointer, so it was not a loop over the batch
at all. A real prefix loop is four instructions (`ldr x2,[x0,#96]; add x0,x0,#512; subs x1,x1,#1;
b.ne L_prefix`). Row 13 goes 3 → 4, the step price 1 + 125.046875/2000 = 1.0625234375 (125
non-refresh, 128 on the 1-in-64 refresh step), and the claim moves **126.2385 → 126.2392**
(+0.0015 units per step). Dependent chain, regenerated in `experiment.py`:
`1.062523 (S_ops + 2^60) < 2^126.237481`, headroom 2^116.52, early-halt < 2^-13.94 = 0.000060,
success ≥ 0.3903 - 0.000060 - 1.1e-7 = 0.39023 > 0.39 (0.39011 at 3x the cycle rate).

## Object reconciliation (re-run, not inherited)

| check | value |
|---|---|
| charged per iteration | 128 (125 non-refresh, 3 refresh-only) |
| literal machine instructions in the object | 126 (122 leaf body + 4 prefix loop) |
| row-6 region | 19 instructions, split 1 entry + 9 x path + 9 y path |
| row-by-row identity mismatches | 0 |
| leaf-body back-edges / prefix-loop back-edges | 1 / 1 |
| exit stubs | 5 distinct `ret`, in-degree 1,1,1,1,3 |
| eligible store→reload pairs | 6, max distance 56, bound 128 |

## Controls

`schedule.py` runs eight planted violations, each of which must make the verifier fail: A
clobber-while-live, B missing preservation, C wrong output register, D census flip, E pinned
clobber, F fusion-bound flip, **G a record operand made `sp`-relative (stack pair)**, **H an exit
edge aliased into the body (the v27 defect)**. All eight fire.

`span_asm.py` re-derives the counts, the per-coordinate segmentation, the branch→target table and
the eligibility-filtered pair table from the object and fails on any mismatch.

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

Lane history: v22 126.2364 (our own `5ecbf15`, the best *valid* published position), v26 126.2385
(after the v26 verdict's cost-model reconciliation), v27 126.2385 (three listing repairs at zero
cost), this v28 126.2392 (+0.0007 for the one instruction the v27 prefix "pass" was missing).
The v27 tree is archived under `workspaces/archived-plain-2026-10-24/` and is superseded; nothing
in this package re-uses the v27 folded-row claim.

## Exact commands

```
python3 schedule.py            # contract: liveness, census 125/128, pinned table, per-coordinate
                               # row 6, eligibility-filtered pairs, branch/stub table, 8 controls
python3 span_asm.py            # emit -> clang -> otool -> reconcile 1:1, enumerate branch edges
python3 audit.py               # 110 arithmetic/direction/literal checks + 3 planted flips
python3 experiment.py <<< '{"seconds": 14}'
python3 ../../scripts/preflight-v28.py
```

## Measured results

* `schedule.py`: peak 30 live values (0 free of 30 usable), 125 charged non-refresh / 128 refresh,
  rows 24/15/14/20/8/19/12/9/1/1/1/4, 12/12 pinned preserved, 6 eligible pairs max distance 56,
  10 excluded candidates, 8/8 controls fire.
* `span_asm.py`: 126 literal instructions (122 leaf + 4 prefix), 0 identity mismatches, 1 leaf
  back-edge, 1 prefix back-edge, 5 exit stubs (in-degree 1,1,1,1,3), `b.ne norefresh` -> 0x1e0
  (the length-update block).
* `audit.py`: 110 checks, 0 findings, 3/3 planted flips detected.
* `experiment.py` (14 s): part A all true (30,000 random + 10,000 adversarial canonical-index
  inputs, 0 mismatches, 0 fallbacks on random); toy pipeline pooled 0.396 +- 0.002 over 58,000
  trials against the claimed 0.39.
* Unsuccessful runs: the v27 tree (rejected, `invalid`); a v28 draft with a 3-instruction prefix
  "loop" that did not advance the record pointer and therefore was not a loop (fixed by paying
  the fourth instruction); a draft that kept 19 charges for row 6 while loading all 24 operand
  limbs at once (rejected by the register checker: 36 > 30).

## Limitations

H1 is a heuristic supported only at toy scale; instruction counts are analytic plus one
object-verified listing; no production-size run is claimed, and rejection does not require one
(the v27 verdict says so explicitly). The fold's interior is inside the fixed reference price and
is not itemised per coordinate beyond the symmetric 1.5M split.

## Next steps

Keep the listing-level lane closed: the remaining objections in this campaign are about the
priced implementation, and v28 answers them against the object rather than the source listing.
