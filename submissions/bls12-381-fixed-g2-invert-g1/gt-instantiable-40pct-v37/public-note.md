# g2-invert v37 -- claim 1277881/10000

## Summary

Recovering `x` from `z = g^x` in the BLS12-381 GT group by a Frobenius-class rho walk that
stores `(A, D)` exponent coordinates and resolves any colliding pair through
`(A2 - q^k A1) + x (D2 - q^k D1) = 0`.  The claim is `1277881/10000` against `W = 5875623158933707187335667578720804574982150` Ir, the
exact sum of 15 ledger rows recomputed by `scripts/cert37.py`.

This revision closes the material defects of the v36 verdict by EXECUTING what it prices:
six-limb borrow-and-reduce arithmetic wired into the fixture, a frequency-independent
per-step charge, `ceil(S/L)` cold trail starts, an occupancy bound on the store, and a
success debit measured on actual collisions inside the shipped loop.

## Method

1. **Six-limb key formation, executed.**  `d28_keys12_form6`, `d28_montmul6`, `d28_l6_add`,
   `d28_l6_sub`, `d28_l6_neg` and `d28_l6_reduce` live in the shipped `dict28.c`; the fixture
   calls them through ctypes on every step, so the priced spans are the executed spans.
   `neg(1)` returns `m-1` in all six limbs.
2. **Frequency-independent step charge.**  Every step pays `max(MUL, SQR) = 47924` Ir; the
   guard-squaring substitution is deleted.
3. **Cold starts, correct inequality.**  Starts are `ceil(S/L)` with `L = 65536` steps -- a
   lower bound that is an equality because the launch rule executes `S` steps at the cap with
   no early abort -- and each is a two-scalar exponentiation `y = g^a z^b` priced at 73084100
   Ir that never reads `x`.
4. **Occupancy-bounded store.**  Reached leaves `= ceil(M(1-e^(-S/M)))` with `M = 2^110`;
   `S/M = 7.38e4`, so the occupancy is `M`, and the dense-leaf formula is shown 37.9x below.
5. **Measured success transfer.**  Relation yield on actual collisions, counted by class in
   the shipped loop: CP-lo 0.912110 at `L = 64`, 0.999937 at the shipped cap, `p_used = 0.527314`.

## Evidence

- `evidence/ledger37.json` -- the rows, the budget rule, the six-limb spans.
- `evidence/arm37.json` -- the in-protocol yield arms at `L = 64` and at the shipped cap.
- `evidence/bounds37.json` -- the occupancy computation and the flagged dense-leaf control.
- `evidence/audit37.json` -- four planted controls: borrow (`neg(1)` and 1000 randoms at the
  381-bit modulus), cold start (never reads `x`, restart row off 7660), occupancy (sparse
  gives `S`, dense gives `M`, the old formula is below both), yield schedule (64, 113,
  65536 must rise).
- `evidence/trap37.json` -- the budget rule standalone.
- `evidence/protocol37.json` -- the smoke run: both collision classes land on the store.
- `certificate.json` -- thirty recomputed checks, `ok: true`.

## Reproducibility

```
cd scripts
python3 bounds37.py     # occupancy -> evidence/bounds37.json
python3 smoke37.py      # the shipped protocol -> evidence/protocol37.json
python3 arm37.py        # in-protocol yield -> evidence/arm37.json
python3 regen37.py      # ledger + claim -> evidence/ledger37.json, submission.json
python3 cert37.py       # certificate.json, ok: true
python3 audit37.py      # four planted controls -> evidence/audit37.json
```

The object is built from the shipped source with `clang -O2 -arch arm64 -c dict28.c -o
dict28.o && clang -shared -arch arm64 dict28.o -o dict28.dylib`; `scripts/regen37.py` reads
the spans from that object's disassembly.

## Objective and environment

Track `bls12-381-fixed-g2-invert-g1`: given the fixed G2 point `Q`, the public `P0`, and the
GT target `z = e(P0,Q)^x`, recover `x` and return `[x]P0`.  Everything here runs on macOS
arm64 with the shipped C object (`clang -O2 -arch arm64`), Python 3 standard library only,
and the toy field (`q = 523`, `r = 521041`) for the census-scale fixtures; the ledger prices
the production instantiation (381-bit `p`, 128-bit radix store, `S = 9.58e+37` steps).

## Prior work and baseline

The baseline is our own v36 package (ticket 177471b7, rejected).  Its fixed-Q reduction,
within-trail collision acceptance and certificate machinery are retained; its four
implementation and accounting defects are replaced, not re-argued.  No external artifact is
adopted, so no inherited development charge is due and `dependencies.json` is empty.

## Implementation, file by file

- `scripts/ir28/dict28.c` -- the shipped object.  v37 appends the six-limb section:
  `d28_l6_sub`, `d28_l6_neg`, `d28_l6_add`, `d28_l6_reduce`, `d28_montmul6` and
  `d28_keys12_form6`, all taking the modulus as an argument.
- `scripts/protocol37.py` -- the fixture.  `Group.state` calls `d28_keys12_form6` (which
  returns the key words and `kappa`), `Group.step` calls `d28_montmul6` and `d28_l6_add` for
  the coefficient update, `Group.start` is the cold two-scalar powering, and `run_one`
  launches trails until `S` steps have been executed.
- `scripts/regen37.py` -- reads the six symbols' disassembly spans from the built object and
  writes `evidence/ledger37.json`, `evidence/trap37.json`, `submission.json`.
- `scripts/bounds37.py` -- the occupancy bound and the flagged dense-leaf control.
- `scripts/arm37.py`, `scripts/audit37.py`, `scripts/cert37.py` -- the yield arms, the four
  planted controls, and the thirty-check certificate.

## Measured results

| row | Ir/step |
|---|---|
| dense Fp12 multiply (frequency-independent) | 47924 |
| key lookup | 959 |
| slot export | 681 |
| record write and insert | 2399 |
| six-limb key formation | 4764 |
| six-limb coefficient update | 3288 |
| **step total** | **60101** |

Cold trail starts: `ceil(S/65536) = 1461847968611600655368414726242031` starts at 73084100 Ir
each.  Store initialisation: 1298074214633706907132624082305024 reached leaves times
`2^22`.  Yield: 14383 actual hits at `L = 64` (CP-lo 0.912110) and 186894 at the
shipped cap (CP-lo 0.999937).

Unsuccessful runs, reported as they happened: the first fixture pass raised a ctypes type
error because the Montgomery constants were handed over as arrays rather than integers (fixed
before any measurement); the first occupancy generator asserted on a wrong dense-regime test
(`64*S >= M`, which does not imply `M*exp(-S/M) < 1`) and was replaced with the exact-integer
test `14*S >= 1100*M`; the smoke fixture recovers `x` from 176 of 179 cross-trail hits but
does not re-derive the target at toy scale, so the claim rests on the arm counts rather than
on a solve.

## Limitations

The yield arms run on a shrunk store (256 retained slots) because collisions must occur at
census scale; the key width is therefore smaller than production, and the arm is a schedule
control, not a production-scale measurement.  The walk's per-step table multiply leaves the
group element in the cyclotomic subgroup only on guard steps, so the dense price is charged
on every step instead of the cheaper square.  Memory is reported, not scored, but the store's
peak is dominated by the reached-node count and is stated with the occupancy derivation.

## Provenance

Built by the g2-invert lane from its own v36 package plus the v36 verdict's defect list; the
fixed-Q reduction descends from the previously accepted formulation on this track.  The
disassembly spans, the occupancy value and the yield numbers in this note are generated by
the shipped scripts, not transcribed.
