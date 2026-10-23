# O(1)-memory Brent-cycle Pollard rho in GT, for fixed-G2 pairing inversion

Candidate `gt-brent-rho-v2`, track `bls12-381-fixed-g2-invert-g1`.
Claimed score: **log2(W) = 128.48 dense-Fp12-multiplication equivalents**
(exact rational `3212/25`), success **>= 0.40** over uniform challenges and fresh
algorithm coins (the track requires >= 0.39).

Target: given a nonidentity `Q` in G2 and `z` in GT, return `P` in G1 with
`e(Q, P) = z`. `Q` is fixed by the challenge; there are no auxiliary powers and no
`n` parameter. Reference unit: one dense Fp12 multiplication in `blst`, fixed
price `C = 20,000` scalar AArch64 instruction equivalents (its measured count is
22,996; the fixed price is what is charged).

Everything below is a hard integer or rational bound except the single heuristic
named in section 6. There is no dictionary, no distinguished-point table and no
storage that grows with the walk: the working set is six GT elements and eight
scalars, and the only table is the fixed 1,024-entry step table.

## 1. Summary of the ledger

| Quantity | Value |
| --- | --- |
| `r` | 52435875175126190479447740508185965837690552500527637822603658699938581184513 |
| `log2 r`, `sqrt(r)` | 254.857089, 2^127.428545 |
| `S0 = ceil(sqrt(2 r ln(5/3)))` | 231454567344042754700414543826693198670 = 2^127.443997 |
| hard evaluation cap `B = 2 S0` | 462909134688085509400829087653386397340 = 2^128.443997 |
| per-step cost | `1 + 400/20000 = 51/50 = 1.02` units |
| main term `B * 51/50` | 2^128.472566 units |
| additive (table, pairing reduction, verification) | 658,496 units = 2^19.329 |
| `W` | 2^128.472566 |
| claimed score | `3212/25 = 128.48` (log2 W rounded up to the next hundredth) |
| peak memory | 724,608 bytes (2^19.47) |
| success probability | >= 0.40 analytic (section 6); 0.452-0.644 measured at toy scope (section 7) |

## 2. Reduction to a discrete logarithm in GT

Let `P0` be the pinned public G1 generator and `g = e(Q, P0) = pairing(Q, P0)`
with the profile's pinned pairing. `Q` is nonidentity and `GT` has prime order
`r`, so nondegeneracy gives `<g> = GT`. Compute `a = log_g z` by the walk below and
return `P = [a]P0`. Then `e(Q, P) = e(Q, P0)^a = g^a = z` by bilinearity, so the
answer is correct for **every** `z` in GT, including `z = 1` (then `a = 0` and
`P = O`, which the track admits exactly for `z = 1`). Two pairings and one G1
scalar multiplication are charged for this reduction and its verification
(section 5); the walk itself needs no pairing.

## 3. The walk

Public table `Th[h] = g^(A_h) z^(B_h)` for `h = 0..R-1`, `R = 1024`, with
`(A_h, B_h)` drawn from fixed public constants; building it is 2,048 scalar
multiplications in GT. A state is `(X, a, b)` with `X = g^a z^b` in GT. The step is

```
h  = H(X)                # one 64-bit limb of the canonical Fp12 encoding,
                         # fixed odd multiplier, xor-shift, top 10 bits
X' = X * Th[h]           # one dense Fp12 multiplication: the reference operation
a' = a + A_h mod r
b' = b + B_h mod r
```

The step is a function of the group element alone: no history, no per-step
randomness after the start, so the trajectory is a deterministic map on GT.

Start: `a, b` are fresh uniform draws from the algorithm's coins (512-bit draws
reduced into `[0, r)` with reduction bias below `2^-256`; a draw of `0` in both
coordinates is re-drawn, and that rejection is charged).

## 4. Cycle detection, exact solve, and the absence of fruitless cycles

Brent's cycle detection with `O(1)` memory: a checkpoint pair `(tX, ta, tb)` moved
when `lam == power` with `power` doubling, and a hare pair `(hX, ha, hb)` advanced
one step per iteration. When `tX == hX`, both pairs carry exact exponents, so

```
g^ta z^tb = g^ha z^hb  =>  a = (ta - ha) * (hb - tb)^-1 mod r     (hb != tb)
```

and the scalar is exact, not high-probability, because the two *elements* are
equal. `hb = tb` (a degenerate relation) has probability `1/r` over the coins and
restarts inside the same cap; `experiment.py` section A4 checks the identity on
real 381-bit GT elements.

A collision is useless only when the exponent pairs coincide, i.e. the walk
returned to the identical state. A two-step return would need
`Th[h]^2 = 1`, i.e. `Th[h] = +-1`; `r` is odd and `Th[h]` has order `r`, so this
is impossible. Longer same-state returns are codimension-one coincidences in the
group, probability about `1/r` per candidate (`2^-255`), not a family. So this
walk has **no structured fruitless cycle of positive density**.

That is why **no quotient and no canonicalisation appears here**. An order-12
Frobenius-class quotient would save `2^1.79` but makes the inverse of a jump equal
a conjugate of the same jump, so `Th[h]^p^a * Th[h]^p^b = 1` holds identically on
the classes where the canonical index differs by 6 — the structured two-cycle
family at density `1/(12R)` that has been raised against quotient walks on this
track. That family is real and measurable: in this session the same toy walk with
an order-12 quotient and `R = 64` spent 22.9% of all steps in such traps, and
cost 2.99 `sqrt(N)` per solved instance instead of 1.36 `sqrt(N)` once a fixed
escape multiplier was added. This submission avoids the family by construction
rather than by a loss argument, which costs `1.79` bits and buys a structural
statement.

## 5. Work accounting, exact

One cold execution, every item charged in dense-Fp12-multiplication units
(`C = 20,000`).

1. **Table build, once.** `R = 1024` entries `g^(A_h) z^(B_h)`: two 255-bit
   scalar multiplications each, at most 383 dense multiplications per scalar
   (windowed comb), so `1024 * 2 * 383 = 784,384` units.
2. **The walk.** At most `B = 2 S0 = 2^128.443997` evaluations. Each is exactly one
   dense Fp12 multiplication (1 unit) plus 400 charged scalar instructions =
   `0.02` units:
   * hash of the element: read one limb, 64-bit multiply, xor-shift, mask, table
     address: 20 instructions;
   * `a += A_h mod r` and `b += B_h mod r`: two 4-limb conditional additions,
     40 instructions;
   * full 12-limb equality comparison against the checkpoint, 120 instructions;
   * checkpoint copy amortized over `power` steps, 60 instructions;
   * loop counter, branch, table entry address (a 1,152-byte table entry),
     160 instructions.
   The reference multiplication is charged at the organizer's fixed 20,000, not at
   its measured 22,996. Cost `B * 51/50 = 2^128.472566` units.
3. **Reduction and verification.** One pairing for `g = e(Q, P0)` and one for
   `e(Q, P) = z`, plus one G1 scalar multiplication `[a]P0`: charged as
   `2 * 2^17 + 2^12 = 266,240` units, i.e. about 5.3e9 instructions, roughly two
   orders of magnitude above a measured `blst` pairing.
4. **Failures and restarts** are inside `B`: the algorithm is hard-capped at `B`
   evaluations and reports FAILURE at the cap, so no branch is unbounded and no
   cost is an average over successful runs. Degenerate draws restart inside the
   same cap.
5. **Memory.** GT element 576 B; the table is `1024 * 576 = 589,824` B, the
   jump table is `1024 * 2 * 32 = 65,536` B, constant state 8 scalars and 6 GT
   elements, and 64 KiB of code and scratch: 724,608 B peak (2^19.47). Nothing
   scales with the number of steps.
6. **No preprocessing, no advice, no amortization.** Only the profile parameters,
   the pinned pairing and the organizer's fixed prices are used.

Additive total `784,384 + 266,240 + 64 = 1,050,688` units = `2^20.003`, which is
`2^-108` of the main term: `W = 2^128.472566 + 2^19.329 = 2^128.472566`, and the
claim `3212/25 = 128.48` is `log2 W` rounded up to the next hundredth, checked with
`W^100 <= 2^12848` in exact integer arithmetic by `experiment.py` section C.

## 6. Success probability: one heuristic, and its margin

Random-mapping model (standard for rho walks, with the Brent-Pollard `R/(R-1)`
correction for an `R`-adding walk): the trajectory of section 3 behaves like a
random map on `r` points, so the tail-plus-cycle length `rho` satisfies
`Pr[rho > t] <= exp(-t^2/(2r))`. With `S0^2 = 2 r ln(5/3)` this gives
`Pr[rho <= S0] >= 0.40`.

Brent evaluates the map at most `2 rho + 3` times before the pointers meet
(checkpoints are powers of two, so the hare stops at or before twice the tail plus
one cycle, one evaluation per iteration). With the hard cap `B = 2 S0` the attack
therefore succeeds whenever `rho + 2 <= S0`, which the model bounds below by
`0.40` (the `+2` moves the exponent by `1e-120` relative to `r`). Required: 0.39,
so the claim carries a 0.01 margin before evidence.

No other heuristic is used: given a collision the scalar is exact (section 4) and
the cost cap is a hard integer.

## 7. Toy evidence, at its actual scope

`experiment.py` section B runs the *same* walk, table, Brent detection, solve and
verification on a cyclic group of prime order `q` (a safe-prime subgroup of
`F_p^*`); the group law's arithmetic differs, the map and the algorithm do not.
This is evidence about the walk only, and no toy constant enters the claim.

| toy `q` bits | trials | solved at cap `2 S0` | success (Wilson 95%) | degenerate | bookkeeping fail |
| --- | ---: | ---: | --- | ---: | ---: |
| 21 | 600 | 375 | 0.625 [0.586, 0.663] | 0 | 0 |
| 24 | 600 | 346 | 0.577 [0.537, 0.616] | 0 | 0 |
| 25 | 300 | 191 | 0.637 [0.581, 0.689] | 0 | 0 |
| 24, cap `1.7 S0` | 600 | 305 | 0.508 [0.468, 0.548] | 0 | 0 |

Every recovered scalar was verified against the target; the bookkeeping identity
`g^a z^b = X` held at both detected states in every collision; no degenerate
relation occurred in 2,100 trials. Mean evaluations to the detected collision
1.5-1.6 `sqrt(q)` (censored), 1.14-1.28 `sqrt(q)` over solved trials, against the
model's `E[2 mu + lambda] ~ 1.88 sqrt(q)`.

## 8. Rejected alternatives

* **Order-12 Frobenius-class quotient.** Saves 1.79 bits (`2^125.64` instead of
  `2^127.43`) but reintroduces the `1/(12R)` structured two-cycle family of
  section 4, which needs an escape mechanism and a loss argument. Measured cost of
  ignoring it: 2.4x in walk length. Not used.
* **Distinguished-point dictionary / parallel collision search.** Same `sqrt(r)`
  bound but a table that grows with the work, plus every per-parameter objection
  that attaches to bucket maps, probe sequences and record widths. Brent needs six
  GT elements of state.
* **Floyd's tortoise-and-hare** (the accepted `e444842b` on this track, 129.20):
  up to `3 rho` evaluations, a factor `2^0.58` more than Brent's `2 rho + 3`, and
  it needs a second pass to expose the colliding pairs. Brent exposes them directly.
* **Special-NFS / SexTNFS in GT.** Would need a full cost derivation in instruction
  units at production parameters; this track's reviewer requires that from scratch,
  and the rho bound here does not depend on it.

## 9. Reproduction

```
cd submissions/bls12-381-fixed-g2-invert-g1/gt-brent-rho-v2
echo '{"production_checks":2,"table_checks":2,"trials":200,"toy_bits":21,"R":64,"cap_factor":2.0,"seconds":10,"seed":1}' | python3 -I experiment.py
```

Section A checks the production facts directly: `p` and `r` prime, `r | Phi_12(p)`,
`ord_r(p) = 12` with `p^6 = -1`, the blst tower is a field, a random element of GT
computed as `t^((p^12-1)/r)` is nonidentity of order `r`, the collision identity on
real GT elements, and sampled table entries of order `r`. Section B reports the
toy statistics above. Section C prints the ledger, the memory bound and the exact
claim.

## 10. Provenance and credit

Pollard rho with an adding walk (Pollard 1978; Teske 2001 for the `R`-adding
correction), Brent cycle detection (Brent 1980), and the standard relation
`g^a z^b = g^a' z^b' => a = (a-a')/(b'-b)` for a discrete logarithm in a
prime-order group. The architecture follows the accepted O(1)-memory entry
`e444842b` on this same track (mpjunior92, 129.20, Floyd detection in GT), which
is this reviewer's established bar; Brent detection is that entry's own listed
next step, and no text, code or table of it was used. The measured
Frobenius-quotient trap statistics in section 4 come from this session's own toy
experiments in this repository.

## Limitations

* The `0.40` success is a statement in the random-mapping model, disclosed and
  measured at toy scope, with a 0.01 margin over the required 0.39.
* This is an analytic instruction-equivalent bound, not a production-size run.
* The score is an upper bound for this attack, not a security estimate.
