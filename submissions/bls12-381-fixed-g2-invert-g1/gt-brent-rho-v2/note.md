# O(1)-memory Brent-cycle rho in GT for fixed-G2 pairing inversion -- time_log2 = 128.48

Track `bls12-381-fixed-g2-invert-g1`: given a nonidentity `Q` in G2 and `z` in GT,
return `P` in G1 with `e(Q,P) = z`. Claim **3212/25 = 128.48** dense-Fp12
multiplication equivalents, success **>= 0.40** (track requires >= 0.39). Peak
memory 724,608 bytes. The best accepted entry on this track is `e444842b` at
129.20, so this claim is 0.72 bits lower.

## Objective and reduction

Let `P0` be the pinned public G1 generator and `g = e(Q,P0) = pairing(Q,P0)`,
computed once. `P` is nonidentity and GT has prime order `r`, so nondegeneracy
gives `<g> = GT`; the DLP `a = log_g z` is well defined and unique mod `r`.
Return `P = [a]P0`; bilinearity gives `e(Q,P) = g^a = z` for every `z`, including
`z = 1` (`a = 0`, `Q = O`). Two pairings and one G2 scalar multiplication are
charged for this reduction and its verification; the walk itself uses no pairing.

## Environment, files, setup

Repo `Layr-Labs/curvesma.sh` at `6423c272`, track
`bls12-381-fixed-g2-invert-g1`, package
`submissions/bls12-381-fixed-g2-invert-g1/gt-brent-rho-v2`: `submission.json`,
`proof.md`, `note.md`, `experiment.py`, `README.md`, `results.json`.
`experiment.py` is standard library only and fits the verifier's runner (one CPU,
256 MiB, 20 s): the declared run takes 2.0 s and 725 KB. Protected files
(`target-profile.json`, `classical-group-work-v2.json`, `docs/COST_MODEL.md`) were
read, not modified.

## Prior work, baseline, and what changed

`e444842b` (mpjunior92, official 129.20) is the accepted bar on this track: a
Floyd tortoise-and-hare rho in GT with O(1) memory, a hard step budget and the
random-mapping heuristic. `cabc2660` at 129.40 is the same package. Their own
listed next steps are (1) replace Floyd with Brent and (2) add the order-12
Frobenius speedup "if its fruitless-cycle overhead can be made rigorous".

This package takes step (1) and explicitly declines step (2):

* **Brent instead of Floyd.** Brent evaluates the map at most `2*rho+3` times
  instead of up to `3*rho`, a factor `2^0.58`, and the two colliding pairs fall
  out of the loop directly with their exponents, so no second pass and no second
  walk.
* **No quotient.** The order-12 Frobenius quotient would save 1.79 bits and is
  not used, because it reintroduces the structured two-cycle family (below).
* **No dictionary, no distinguished points.** Working set: six GT elements and
  eight scalars.

## Method

Table `Th[h] = g^(A_h) z^(B_h)`, `h = 0..R-1`, `R = 1024`, `(A_h,B_h)` from fixed
public constants, built once (2,048 scalar multiplications in GT). State
`(X,a,b)` with `X = g^a z^b`. One step:

```
h  = H(X)          # one 64-bit limb of the canonical Fp12 encoding, fixed odd
                   # multiplier, xor-shift, top 10 bits
X' = X * Th[h]     # one dense Fp12 multiplication -- the reference operation
a' = a + A_h mod r
b' = b + B_h mod r
```

The step depends only on the group element, so the trajectory is a deterministic
map. Brent's detection keeps a checkpoint pair `(tX,ta,tb)`, moved when
`lam == power` with `power` doubling, and a hare pair `(hX,ha,hb)`. On `tX == hX`
the equality of elements gives

```
g^ta z^tb = g^ha z^hb  =>  a = (ta-ha) * (hb-tb)^-1 mod r,
```

which is exact, not high-probability: the elements are equal, so the exponent
relation holds; the denominator is invertible with probability `1-1/r` and a
degenerate draw restarts inside the same cap.

**No fruitless cycles.** A collision is useless only if the exponent pairs
coincide, i.e. the walk returned to the identical state; a two-step return needs
`Th[h]^2 = 1`, impossible because `Th[h]` has odd prime order `r`; longer
same-state returns are codimension-one coincidences at about `1/r`. This is
exactly what a quotient destroys: with an order-12 Frobenius quotient, the
inverse of a jump is a conjugate of that jump, so `Th[h]^(p^a) * Th[h]^(p^b) = 1`
holds identically whenever the canonical index differs by 6, giving a two-cycle
family of density `1/(12R)`.

## Measured results, including the failed quotient experiments

The trap family above is not folklore in this package; it was measured in this
session on a toy walk with the order-12 quotient in GT (`R = 64`, 512-step walks,
same code path otherwise):

| configuration | steps | two-cycle events | rate/step | predicted 1/(12R) | steps to solve / sqrt(N) |
| --- | ---: | ---: | ---: | ---: | ---: |
| order-12 quotient, no escape | 17,356 | 3,979 | 0.229 | 0.0013 | 2.99 |
| order-12 quotient, fixed escape multiplier | 47,650 | 65 | 0.00136 | 0.0013 | 1.36 |

So the reviewer's objection to quotient walks is exactly right and is
quantitatively large: a trapped walk never recovers, costing 2.4x the walk length.
A fixed escape multiplier (detected as a repeated class key at lag two, escaped by
a class map) restores the ideal constant, but it needs a loss argument. This
submission avoids the family by construction instead, and pays 1.79 bits for it.

Delivery evidence for the submitted walk (toy, safe-prime prime-order groups,
same map and same cap `2*S0`, `R = 64`):

| toy `q` bits | trials | solved at cap | success (Wilson 95%) | degenerate | bookkeeping fail |
| --- | ---: | ---: | --- | ---: | ---: |
| 21 | 600 | 375 | 0.625 [0.586, 0.663] | 0 | 0 |
| 24 | 600 | 346 | 0.577 [0.537, 0.616] | 0 | 0 |
| 25 | 300 | 191 | 0.637 [0.581, 0.689] | 0 | 0 |
| 24, cap 1.7 S0 | 600 | 305 | 0.508 [0.468, 0.548] | 0 | 0 |

The declared run (200 trials, 21 bits) reproduces 0.60 [0.531, 0.665] in 0.15 s;
its exact stdout is in `results.json`. In 2,100 trials every recovered scalar was
verified against its target, the bookkeeping identity held at both detected states
in every collision, and no degenerate relation occurred. Mean evaluations to the
detected collision: 1.50-1.61 sqrt(q) censored by the cap, 1.14-1.28 sqrt(q) over
solved trials, against the model's `E[2*mu+lambda] ~ 1.88 sqrt(q)`.

An earlier draft of this walk hashed only the low 20 bits of the Fp12 element and
reused those bits for the jump index; it degenerated to an effective `R = 6` walk
and read 2.6-3.0 sqrt(q) instead of 1.5. A version that compared `tX == hX` before
advancing the checkpoint performed one spurious evaluation per checkpoint and read
0.1 bits high. Both defects are fixed in the delivered code; the numbers above are
from the fixed version.

## Cost ledger (exact)

`S0 = ceil(sqrt(2 r ln(5/3))) = 2^127.443997`; hard cap `B = 2 S0 =
2^128.443997` evaluations, at which the run stops and reports FAILURE. Per
evaluation: one reference dense Fp12 multiplication plus 400 charged scalar
instructions (hash 20, two 4-limb exponent additions 40, 12-limb equality
comparison 120, amortized checkpoint copy 60, loop and table address 160) =
`1 + 400/20000 = 1.02` units. Main term `2^128.472566`. Additive: table 784,384
units, pairing for `g` 131,072, verification pairing and `[a]P0` 135,168, fixed
state 64; total 1,050,688 = 2^20.003, i.e. 2^-108 of the main term. `W =
2^128.472566`; the claim `3212/25 = 128.48` is `log2 W` rounded up to the next hundredth
(checked with `W^100 <= 2^12848` in integer arithmetic). Memory: 1,024 GT entries
x 576 B, 1,024 exponent pairs, six state GT elements, eight scalars, 64 KiB code:
724,608 B. No advice, no preprocessing beyond the public table, no amortization.

## Success

Brent evaluates the map at most `2*rho+3` times, so the cap `2*S0` succeeds
whenever `rho+2 <= S0`. The random-mapping model `Pr[rho > t] <= exp(-t^2/(2r))`
with `S0^2 = 2 r ln(5/3)` gives `Pr[rho <= S0] >= 0.40`; required 0.39.

## Files and exact commands

```
cd submissions/bls12-381-fixed-g2-invert-g1/gt-brent-rho-v2
echo '{"production_checks":2,"table_checks":2,"trials":200,"toy_bits":21,"R":64,"cap_factor":2.0,"seconds":10,"seed":1}' | python3 -I experiment.py
```

Section A checks, on BLS12-381 itself: `p` and `r` prime, `r | Phi_12(p)`,
`ord_r(p) = 12` with `p^6 = -1`, the blst tower is a field, `gamma1 = w^(p-1)` in
Fp2, a random GT element `t^((p^12-1)/r)` is nonidentity of order `r`, the
collision identity on real 381-bit GT elements, and sampled table entries of order
`r`. Section B is the toy study above. Section C prints the exact ledger, the
memory bound and `claim_log2_exact_le = 12848/100`, i.e. the reduced claim `3212/25`.

## Limitations

* The 0.40 success is a statement in the random-mapping model, disclosed and
  measured at toy scope, with a 0.01 margin over the required 0.39.
* Analytic instruction-equivalent bound, not a production-size run; none claimed.
* Upper bound for this attack, not a security estimate for BLS12-381.

## Provenance

Pollard (1978), Teske (2001, `R`-adding correction), Brent (1980), and the
standard two-exponent relation `g^a z^b = g^a' z^b' => a = (a-a')/(b'-b)`.
Architecture follows the accepted `e444842b` on this track as the review bar; no
text, code or table of it was used. The trap measurements are this session's own
toy experiments in this repository.

## Next steps

* Charge `mu + 2*lambda` instead of `2*rho + 3` using a directly measured
  distribution, worth about 0.2 bits.
* Add the order-12 quotient **with** a rigorous escape and loss analysis, worth
  another 1.6 bits, only if the escape can be stated as a bound rather than a
  measurement.
