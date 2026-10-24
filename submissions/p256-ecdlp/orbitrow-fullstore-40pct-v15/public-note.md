# p256-ecdlp v15 — full-store rho walk with the v14 review flaws repaired

## Objective

Claim the lowest defensible log2 of total classical work for recovering x from
Q=[x]P on NIST P-256 by a coefficient-free full-store rho walk, with algorithmic
success at least 0.39, under the protected p256-ecdlp accounting (C = 1000 scalar
AArch64 instructions per reference unit). The claim is 2^127.9570, shipped as the
reduced rational 127957/1000.

## Prior work and what was rejected

v13 was rejected on four flaws (four-cycle channel unhandled, selector depth, step
price, ledger). v14 repaired several and cleared intake as ticket 7db962fb, then was
rejected at review with six flaws: (1) `v14kit.cls` used `min(x, p-x)`, conflating
points with opposite x-coordinates, which are not negatives of each other; (2) the
fixture order 1000839 is composite and the inverse test used Fermat
`pow(db, n-2, n)`, invalid for a composite modulus; (3) the four-cycle scan imposed
only `[A,B,A,B]` and then applied `min(1,ratio)`, discarding the measured 1.779x
inflation; (4) the selector experiment compared the successor with a variable that
the update order had already set to the current state, so it measured a one-step
return; (5) the 147-instruction `fe_mul` price had no object behind it; (6) the
ledger row said 133.25 while the 48-byte record requires 133.2549625. Every one of
those is repaired below against a measured artifact, and because fixes (1) and (6)
change the numbers, the claim moves from 127.8076 to 127.9570. That is a worse
claim and it is the honest one.

## 1. Class semantics fixed in code, not in prose

The shipped sign rule canonicalises a point by negating it iff y is odd. Negation
preserves x, so the class key in production is x itself: P and -P collide, and x
and p-x are different keys belonging to unrelated points. v15's `cls(Pt)` returns
`Pt[0]`, exactly the 32 bytes `record_pack` writes and exactly what the collision
algebra assumes. v14's `min(x, p-x)` merged x with p-x; that partition is gone.
Because the partition changed, the census and the cycle scan were both re-run from
scratch on the corrected code, and their measured counts differ from v14's.

## 2. Fixture primality disclosed, inverse switched to extended Euclid

The fixture group order is 1000839 and it is composite: 3 * 333613. The package
states this in `evidence/census_v15.json` along with a Miller-Rabin verdict and the
factorization; there is no primality assertion anywhere. Inversion is extended
Euclid (`pow(a, -1, n)`) on units, which needs no primality; the Fermat sentence is
deleted. Re-verified: 0 misses in 649969 inverse trials, with 278870 non-units
skipped and disclosed. The v14 count of 0/649969 was computed under the same
modulus with FLT; the corrected method reproduces the zero and removes the invalid
justification.

## 3. Four-cycle channel: full 4-tuple enumeration, ratio applied

`experiments/cyclescan.py` no longer imposes a shape. For each of 16 fixture seeds
it builds all 4096 two-step prefix states, both preimages of S under every one-step
map, and probes all 64 third rows against the preimage table: every 4-step window is
tested, 268435456 windows in total. Every hit is replayed with coefficient tracking,
so a return with zero coefficient delta (a death) is separated from a return with
nonzero delta (a usable relation). Measured: the alternating family reproduces at
ratio 4.024, the reversed family at 3.969, and a chance family of 680583 hits whose
production shape probability is 2^-256 rather than R^-2, so its charge is 2^-243.
The production charge is `per4 = (4.024 + 3.969) * 2^-4 * R^-2 + 40599.9 * 2^-4 *
2^-256`, i.e. a relative work loss of 6.8e-13. There is no `min(1, ratio)` anywhere
in the source; the ratio multiplies.

## 4. Two-step history test

The production schedule's trap is a revisit after two steps. v14 compared the
successor with `hist0` after the update had set `hist0 = S`, measuring a one-step
return. v15 compares `S_{t+1}` with `S_{t-1}` explicitly, with escape enabled and a
cluster bootstrap over whole walks: delta 1.70e-4, CI [-1.14e-3, 1.48e-3], traps 303
against 297. The interval contains zero, the arm is inconclusive, selector transfer
is declared OPEN and no claim row depends on it.

## 5. The multiplication price now has an object behind it

`evidence/p256_fe_v15.c` is the v10m field-kernel source verbatim plus the counted
wrappers; `evidence/insn_disasm.txt` is the objdump of the build
(`clang -O2 -fPIC -arch arm64`), with `evidence/insn_nm.txt` carrying the symbol
addresses. `fe_mul_obj` spans 148 instructions (0x0..0x250), `fe_add_obj` and
`fe_sub_obj` 27 each, so the affine-addition row is 3*148 + 6*27 = 606 Ir.
`evidence/trace_counts.json` records the dynamic call counts: one affine addition
executes 388 fe_mul with its own inversion (3 plus a 385-multiplication chain: 256
squarings, 128 multiplies for popcount(p-2), one fixup) and 3 fe_mul when the
inverse is supplied by the batch. v14's 147 is superseded by the measured 148, and
the per-visit total moves from 1100 to 1220.054 Ir.

## 6. Ledger regenerated from the shipped 48-byte record

Rows in `evidence/ledger_v15.json`: affine addition 606; batch-inversion share 444
(3 muls) plus 0.054 (the 385-mul chain over a 2^20 batch); table read 10; canonical
sign rule 20; 48-byte append 12; radix sort and adjacent detect 128. Total
1220.054 Ir per visit = 1.220054 units, so work = 127.67 + log2(1.220054) = 127.9569,
plus restart 6.8e-13 and replay 2^-19, ceilinged at four decimals to 127.9570.
Storage is recomputed honestly: 48 bytes times 2^127.67 visits is 133.2549625 log2
bytes, i.e. 133.2550 as an upper bound; the v14 row's 133.25 was rounded down and is
superseded. Storage does not bind, so claim = min(bounds) = 127.9570, and the
in-package assertion `claim >= log2(total)` holds. The required literal sweep is
`experiments/literal_sweep.py`: every integer of five or more digits in this package
must reconcile to an owning evidence value, and it does.

## Success, stated at its true scope

`tracks/p256-ecdlp/rules.md`: "Accept the complete claim only if it is concrete,
evaluable, supported and has no surviving fatal flaw. Supported plausible heuristics
may qualify." and "Algorithmic success must be at least 0.39 over uniform
full-production challenges and fresh algorithm coins." and "An analytic argument
need not run an exponentially large attack. Small experiments are evidence only at
their actual scope, never a production complexity proof."

The claim's success law is therefore the birthday heuristic: lambda = V^2/r with
V = 2^127.67 visits and r ~ 2^256, giving lambda = 0.6329 and success 0.4689, above
the 0.39 floor. The toy census is reported at its own scope and it does not reach
that regime: the toy walk recurs in roughly sqrt(1000839)/7 ~ 140 steps, so a
lambda = 0.63 toy store is cycle-dominated (measured 0.297, CP-95 0.260). v15
reports that number instead of hiding it, and does not use it as the production
value; the production walk takes 2^107.67 steps per walker in a 2^256-point group,
far below its 2^128 recurrence scale. Record capacity is demonstrated by replay
rather than by identity serialization: 38400 of 38400 reconstructions from the
unpacked 48-byte record equal truth, and a store collision between two replayed
coefficient pairs recovers the secret, verified by cand*P == Q.

## Files and reproduction

* `argument.md` — the six repairs with their measured numbers.
* `experiments/v15kit.py` — curve, x-class semantics, extended-Euclid inverse, 48-byte record.
* `experiments/insn_price.py` — parses the shipped objdump and nm, emits `evidence/insn_price.json`.
* `experiments/fullstore_census.py` — shared-store census, per-walker seeds, CP-95, inverse re-verification.
* `experiments/cyclescan.py` — full 4-tuple scan and ratio pricing, emits `evidence/cycle_channel_v15.json`.
* `experiments/record_recovery.py` — replay recovery through the 48-byte layout.
* `experiments/selectors.py` — two-step history arm, exploratory only.
* `experiments/ledger_v15.py` — the regenerated bound vector and the in-package assertions.
* `experiments/literal_sweep.py`, `experiments/audit.py` — literal reconciliation and independent controls.

Commands: `python3 experiments/insn_price.py && python3 experiments/fullstore_census.py
&& python3 experiments/cyclescan.py && python3 experiments/record_recovery.py &&
python3 experiments/selectors.py && python3 experiments/ledger_v15.py && python3
experiments/literal_sweep.py && python3 experiments/audit.py`.

## Negative results and open problems

* The generalized scan's chance family is large at fixture scale (680583 hits) but
  cannot bind at production because its shape probability is 2^-256.
* A toy group cannot test the birthday law at lambda = 0.63: the recurrence scale and
  the birthday scale are both sqrt(group size). This is a scope limit, stated.
* Selector transfer remains OPEN; v14 and v15 both fail to separate the arms.
* The success claim is a heuristic with fixture-scope support for the detector and
  the record, not a production demonstration.
