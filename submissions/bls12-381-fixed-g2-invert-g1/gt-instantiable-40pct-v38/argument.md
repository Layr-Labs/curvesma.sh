# g2-invert v37 -- argument

Claim `1277881/10000` on `bls12-381-fixed-g2-invert-g1`: recover `x` from `z = g^x` in the GT
group of BLS12-381 by a Frobenius-class rho walk on 1024 workers, returning `[x]P0` verified
by `e(P,Q) = z`.  `W = 5875623158933707187335667578720804574982150` Ir is the exact sum of 15 rows in
`evidence/ledger37.json`; `scripts/cert37.py` recomputes every one and both directions of the
claim inequality.

## The fixed-Q reduction is unchanged

`g = e(P0, Q)` on the fixed G2 point `Q` turns the target into a discrete log in the
order-`r` subgroup of GT.  A walk step multiplies by a table element and tracks the exponent
as a linear form `(A, D)` in the public generators `g` and `z`, so two records that collide
in the store give `(A2 - q^k A1) + x (D2 - q^k D1) = 0`, and `x` follows from one division
when the denominator is nonzero.  This part was accepted before and is untouched.

## Flaw 1 closed: priced arithmetic is executed arithmetic

`Group.state` and `Group.step` CALL the shipped six-limb routines through ctypes on every
step: `d28_keys12_form6` forms the twelve lanes and returns both the six key words and
`kappa`; `d28_montmul6` performs `A <- A + linv[kap]*a (mod r)` once per coefficient with the
constant held in Montgomery form; `d28_l6_add` accumulates.  The Python path
(`Tower.keys12` then `min/index`) is gone from production and `cert37.py` asserts its absence.

The routines propagate borrow and reduce: `d28_l6_sub` borrows across all six limbs and adds
`m` back on underflow; `d28_l6_neg` returns `m - a` and maps `0` to `0`, so the planted case
`neg(1)` returns `m - 1` in all six limbs (v36 gave `2^64 - 1` with zero upper limbs);
`d28_l6_add` carries and conditionally subtracts `m`; `d28_l6_reduce` handles a carry limb;
`d28_montmul6` is a six-limb Montgomery CIOS with a six-limb final subtraction.  The modulus
is an ARGUMENT: the fixture runs the toy modulus and `scripts/audit37.py` runs the same
symbols at the 381-bit production modulus `p` -- 1000 random additions, subtractions and
negations match the Python reference exactly, the carry-limb reductions match, and the old
limb-local rule is evaluated beside the new one to show it fails the planted case.

## Flaw 3 closed: budget, inequality direction, cold start

Every executed step pays one price, `max(dense-multiply span, cyclotomic-square span)` =
47924 Ir: the guard-squaring substitution `ceil(4497*S/5040000)` is deleted.

A trail of length at most `L` forces AT LEAST `ceil(S/L)` starts; v36 charged at most that
many and started only 1024 trails, which cannot execute `S` steps.  The shipped loop runs
until `S` steps have been executed, every trail runs to the cap and no early-abort path
exists on the launched schedule, so starts = `ceil(S/L)` with equality by construction
(`L = 65536`).  Each start is cold: `y = g^a z^b`, a two-scalar exponentiation in the public
`g` and the target `z`, priced at 1525 full multiplies = 73084100 Ir.  v36 priced 7660 Ir for
a start that read the secret `x`.  `audit37.py` checks the start against `g^(a+xb)` computed
from `x` only afterwards and asserts the source never reads it.

## Flaw 4 closed: reached-node bound from occupancy

The 128-bit slot is consumed 22 bits per level over `h = 6` levels, node = 32 MiB.  The
reached-leaf count is the uniform-address occupancy

    reached(M, S) = ceil( M * (1 - exp(-S/M)) ),   M = 2^110

with `S/M = 7.3804e+04`, so `M*exp(-S/M) < 1` and the occupancy is `M` itself
(1298074214633706907132624082305024 leaves); the init row is reached nodes times `2^22`.  The
dense-leaf formula `6*ceil(S/(4*2^22))` gives 34262061764334390360197220146298, a factor 37.9
BELOW the occupancy, and is reported as a flagged control, not a bound.  The 96-byte record
node and the 32 MiB radix node are separate constants, and every calloc argument is below
`SIZE_MAX`.

## Flaw 5 closed: yield on actual collisions, in the shipped loop

`evidence/arm37.json` comes from `scripts/protocol37.run_one` -- the same walk, store code and
recovery the claim prices -- with the store shrunk so ACTUAL collisions occur at census
scale.  Each hit is a real pair of the walk; the arm records its class, whether
`D2 - q^k D1` is nonzero, and whether the recovered `x` re-derives the target.  No pair is
planted with a chosen denominator.  Yield at trail length 64: 14383 hits, CP-lo
0.912110; at 65536 (the shipped cap): 186894 hits, CP-lo 0.999937.  The yield RISES with
the trail length, which is the schedule sensitivity the debit needs, and the debit takes the
worse of the shipped value and the `L = 64` control: `p_used = 0.578125 * 0.912110 = 0.527314`,
above the 0.39 floor.

## What is charged

Fifteen exact-integer rows: the per-step dense multiply, lookup, slot export, record write
and insert, six-limb key formation, six-limb coefficient update, the ctypes pack and unpack,
guard and meter evaluation, the widened 128-bit references, the cold trail starts, the store
initialisation from the occupancy, worker starts, collision processing, code and advice
bytes, and the final verification.
