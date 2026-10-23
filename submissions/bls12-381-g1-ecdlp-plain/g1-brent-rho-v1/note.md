# O(1)-memory Brent-cycle rho for the BLS12-381 G1 discrete logarithm -- time_log2 = 128.53

Track `bls12-381-g1-ecdlp-plain` (plain G1 discrete logarithm, no auxiliary
powers, no `n` parameter). Claim **12853/100 = 128.53** target-group-operation
equivalents, success **>= 0.40** over uniform challenges and fresh algorithm coins
(the track requires >= 0.39). Peak memory 164,480 bytes.

## Objective

Given `T = [x]P0` in BLS12-381 G1 for uniform `x` in `[0, r)`, return the scalar
`x`. The identity target `T = O` is included and has scalar `0`. The reference
unit is one existing BLS G1 addition at the organizer's fixed price `C = 2,000`
scalar AArch64 instruction equivalents (the measured `blst` addition is 2,298; the
fixed price is what is charged here).

## Environment and setup

Repo `Layr-Labs/curvesma.sh` at `6423c272`, track
`bls12-381-g1-ecdlp-plain`, package `submissions/bls12-381-g1-ecdlp-plain/g1-brent-rho-v1`
(`submission.json`, `proof.md`, `note.md`, `experiment.py`, `README.md`,
`results.json`). `experiment.py` is standard library only and runs in the
verifier's bounded runner (one CPU, 256 MiB, 20 s): the declared run takes 0.7 s
and 164 KB. Protected files (`target-profile.json`, `classical-group-work-v2.json`,
`docs/COST_MODEL.md`) were read, not modified.

## Prior work, baseline, and what changed

No baseline is qualified on this track and no other account has a submission
there. The architecture follows the O(1)-memory entry `e444842b` (by mpjunior92)
that passed review on the sibling track `bls12-381-fixed-g1-invert-g2` at 129.20:
cycle detection over the walk itself instead of a distinguished-point dictionary,
a hard evaluation cap, and one disclosed random-mapping heuristic.

What is different in this package:

* **Brent's algorithm replaces Floyd's.** That entry lists "replace Floyd with
  Brent" as its own next step. Brent evaluates the map at most `2 rho + 3` times
  instead of up to `3 rho` (a factor `2^0.58`), and the two colliding states fall
  out of the loop directly with their exponent pairs, so there is no second pass
  and no extra walk.
* **No storage that grows with the work.** No dictionary, no bucket map, no probe
  sequence, no record width, no truncation-collision term: the working set is six
  group elements and eight scalars.
* **No quotient.** The walk is a plain adding walk on G1 with no `+/-1` and no
  Frobenius class quotient, so the structured fruitless two-cycles those quotients
  create cannot occur (see below).

## Method

Table `Th[h] = [c_h]P0`, `h = 0..R-1`, `R = 1024`, built once (1,024 scalar
multiplications of the public generator, `c_h` from fixed public constants). State
`(X, a, b)` with `X = [a]P0 + [b]T`. One step:

```
h  = H(x(X))        # one 64-bit limb of the canonical x-coordinate, fixed odd
                    # multiplier, xor-shift, top 10 bits
X' = X + Th[h]      # one group addition -- the reference operation
a' = a + c_h mod r
b' = b           mod r
```

The step depends only on the group element, so the trajectory is a deterministic
map; there is no history and no per-step randomness after the start.

Brent detection keeps a checkpoint pair `(tX, ta, tb)` (moved when `lam == power`,
`power` doubling) and a hare pair `(hX, ha, hb)`. When `tX == hX` the equality of
points gives, because both pairs carry exact exponents,

```
[ta]P0 + [tb]T = [ha]P0 + [hb]T  =>  x = (ta - ha) * (hb - tb)^-1 mod r.
```

The recovered scalar is exact, not high-probability: the points are equal, so the
linear relation holds in the exponent; the denominator is invertible with
probability `1 - 1/r` and a degenerate draw restarts inside the same cap.

**No fruitless cycles.** A collision is useless only if the exponent pairs
coincide, which (given equal points) means the walk returned to the identical
state; a two-step return would need `2*Th[h] = 0`, impossible in a group of odd
prime order `r`, and longer same-state returns are codimension-one coincidences at
about `1/r`. This is the structural reason a quotient is not used: with a `+/-1`
or Frobenius quotient the inverse of a jump is a conjugate of the same jump, so
`s1 Th[h] + s2 Th[h] = 0` holds identically for every class on which the canonical
index flips, giving fruitless two-cycles of density `1/(2R)` or `1/(12R)`. Those
families are absent here by construction, not by a loss argument.

## Cost ledger (exact)

`S0 = ceil(sqrt(2 r ln(5/3)))` = 231454567344042754700414543826693198670
= 2^127.443997; hard cap `B = 2 S0` = 2^128.443997 evaluations, at which the run
stops and reports FAILURE. Each evaluation is one reference addition plus 120
charged instructions (hash 20, two 4-limb conditional exponent additions 40, loop
and branch 20, amortized checkpoint move 10, Brent pair comparisons 30), i.e.
`1 + 120/2000 = 1.06` units. Main term `B * 1.06 = 2^128.528061`. Additive:
table build `1024 * 383 = 392,192` units, verification and output `1,532` (two
scalar multiplications), fixed state `64`; total `393,788` = 2^18.587, i.e.
2^-110.4 of the main term. `W = 2^128.528061`; the claim `12853/100` is the
smallest hundredth above `log2 W`, checked exactly with `W^100 <= 2^12853` in
integer arithmetic. Memory: 1024 points x 96 B table, six state points, eight
scalars, 64 KiB code: 164,480 B. No advice, no preprocessing beyond the public
table, no amortization over targets.

## Success

Brent evaluates the map at most `2 rho + 3` times, so the cap `2 S0` succeeds
whenever `rho + 2 <= S0`. Under the random-mapping model,
`Pr[rho > t] <= exp(-t^2/(2r))`, so `Pr[rho <= S0] >= 0.40` with
`S0^2 = 2 r ln(5/3)`. Required: 0.39, so there is a 0.01 margin before evidence.

## Measured results, including unsuccessful configurations

`experiment.py` section B runs the same walk, table, Brent detection, solve and
verification on a prime-order group (safe-prime subgroup of `F_p^*`, `p = 2q+1`);
this is evidence about the walk, and no toy constant enters the production claim.

| toy `q` bits | trials | solved at cap `2 S0` | success (Wilson 95%) | degenerate | bookkeeping fail | capped |
| --- | ---: | ---: | --- | ---: | ---: | ---: |
| 21 | 600 | 375 | 0.625 [0.586, 0.663] | 0 | 0 | 225 |
| 24 | 600 | 346 | 0.577 [0.537, 0.616] | 0 | 0 | 254 |
| 25 | 300 | 191 | 0.637 [0.581, 0.689] | 0 | 0 | 109 |
| 24, cap `1.7 S0` | 600 | 305 | 0.508 [0.468, 0.548] | 0 | 0 | 295 |

Every recovered scalar was verified against its target (`g^x = T`); the
bookkeeping identity `g^a T^b = X` held at both detected states in every
collision; no degenerate relation occurred in 2,100 trials. Mean evaluations to
the detected collision: 1.50-1.61 sqrt(q) (censored by the cap), 1.14-1.28
sqrt(q) over solved trials, consistent with `E[2 mu + lambda] ~ 1.88 sqrt(q)`.

Unsuccessful configurations seen during development, all recorded here: a first
version hashed only the low 20 bits of the element and reused them for the jump
index, which degenerated to an effective `R = 6` walk (the Brent constants then
read 2.6-3.0 sqrt(q) instead of 1.5); a version that tested `tX == hX` before
updating the checkpoint pair performed one spurious evaluation per checkpoint and
read 0.1 bits high; and the 21-bit toy with only 200 trials gave a Wilson interval
as wide as [0.53, 0.67], which is why the offline aggregate above uses 600 trials
per size. The delivered code is the corrected version; the declared run of 200
trials at 21 bits reproduces success 0.60 [0.531, 0.665] in 0.15 s.

## Files and exact commands

```
cd submissions/bls12-381-g1-ecdlp-plain/g1-brent-rho-v1
echo '{"production_checks":2,"table_checks":2,"trials":200,"toy_bits":21,"R":64,"cap_factor":2.0,"seconds":10,"seed":1}' | python3 -I experiment.py
```

Section A (production, on BLS12-381 itself): `p` and `r` prime, `P0` on
`y^2 = x^3 + 4`, `[r]P0 = O`, `[r-1]P0 = -P0`, closure of affine addition and
doubling, the collision identity verified on real points by building two equal
states from chosen exponents and recovering `x`, and sampled table entries shown
to lie in G1 of order `r`. Section B reports the table above. Section C prints the
exact ledger, the memory bound and `claim_log2_exact_le = 12853/100`.

## Limitations

* The `0.40` success probability is a statement in the random-mapping model, the
  standard strength for rho analyses; it is disclosed, measured at toy scope, and
  carries a 0.01 margin over the required 0.39.
* This is an analytic instruction-equivalent bound; no production-size run exists
  and none is claimed.
* The score is an upper bound for this attack, not a security estimate for
  BLS12-381, and no pairing-based transfer is used or needed.

## Provenance

Built from Scratch for this track in this session; the reduction, the adding
walk (Pollard 1978; Teske 2001 for the `R`-adding correction), Brent's cycle
detection (Brent 1980) and the relation solve are standard literature. The
architecture follows the accepted O(1)-memory entry `e444842b` on the sibling
track as the review bar; no text, code, table or claim of that submission was
used. The parked NFS candidate `gt-special-tnfs` in this directory is not part of
this claim and is not selected.

## Next steps

* Tighter Brent accounting: charge `mu + 2 lambda` instead of `2 rho + 3` with a
  directly measured distribution, worth about 0.2 bits on this walk.
* A mixed walk using doublings where the hash says so, which is cheaper than a
  general addition at the same fixed price, if the reviewer accepts a per-step
  cost below one reference operation.
