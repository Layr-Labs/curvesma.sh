# Accounted STNFS for fixed-G1 pairing inversion

Track `bls12-381-fixed-g1-invert-g2`. Given nonidentity `P` in G1 and `z` in GT, return `Q` in G2 with `e(P,Q)=z` for that exact `P`. `Q` is infinity exactly when `z=1`. The unit is one dense Fp12 multiplication, `C=20000` scalar AArch64 instructions. No auxiliary powers. This package keeps the same tower and the same `(A,B)=(650, 2^66)` as rejected submission `9cb00ec2`, and replaces the four pieces that rejection said were not a bound. The resulting total is higher. It is not a tightening of `116.337899`.

## 1. Reduction

Let `Q0` be the public G2 generator and `e` the pinned reduced pairing (`pairing(Q0,P)` in py_ecc order).

1. Reject `P` unless it is a nonidentity point of the order-`r` subgroup, and reject `z` unless `z^r=1`.
2. If `z=1`, return infinity.
3. Set `g=e(P,Q0)`. Nondegeneracy and prime order make `g` a generator of GT, so `z=g^a` for a unique `a` mod `r`.
4. Compute `a` by the STNFS below in the degree-12 field.
5. Return `Q=[a]Q0` only if a fresh evaluation gives `e(P,Q)=z`.

The field isomorphism and the final pairing are charged in the wrapper. A wrong logarithm fails that check.

## 2. What changed relative to `9cb00ec2`

1. **Filtering.** The old simulator replaced every row through a light column by one symmetric difference and reset coefficients to 1. That deletes constraints and is not elimination over `F_r`. The new `filter_sim` works over the prime `2^61-1`: a weight-1 column drops its row, and a weight-`k` column keeps `k-1` combinations `row - c*pivot`. The production ledger takes **no** reduction credit from that simulator. The matrix that is costed is the unfiltered square system on the factor-base orbits.
2. **Descent.** The schedule that fixed both cofactors at `300+12 log2(9)` bits is withdrawn. Each of the two target logarithms is charged one extra copy of relation collection and one extra copy of linear algebra, the same formulas as the main phase. That charge is deterministic. There is no Markov tail on an underived tree.
3. **Yield arithmetic.** Primitive density multiplies by `1/zeta=0.903208`. The old code divided by that constant. Alpha is charged at the measured value **plus** `1.5` standard errors, because a larger alpha lowers Dickman `rho`. The second-order term `(1-gamma) rho(u-1)/ln N` is evaluated in the Monte Carlo, per side. The sieved set is the exact integer domain below, not the non-integer real number `(2A+1)^12/4`.
4. **Kernels and prime ideals.** `kern.c` subtracts on one sign and adds on the other, propagates carry and borrow, and subtracts `r` until the accumulator is reduced. The ledger does not use a callgrind count of the old kernel. Hit and entry prices are analytical upper bounds above the old measured inner loops. The mean number of roots of `h` modulo primes below 3000 is computed in the experiment (`0.958` on 430 primes). The charge uses `1`, which is above that mean. The old `2/3` mixed a splitting probability with a mean and is not used.

## 3. Field, domain, norms

`h(Y)=Y^6-Y^2+1=g3(Y^2)` with `g3(Z)=Z^3-Z+1`, `f` and `g` as in the claim statement, `u=-0xd201000000010000`. The experiment checks `p` and `r` from `u`, and that `f` vanishes at `X^2=uy` modulo `(p,h)`.

The sieved set is

```
D = { (a,b) : coordinates in [-A,A], the first nonzero coordinate of a is positive, b_5 in 1..A }.
```

Its cardinality is the integer `((2A+1)^6-1)/2 * A * (2A+1)^5`. At `A=650` that is `5873995684991862158037873511093785000`. Pairs with `b_5=0` are not sieved. Norms are unchanged by `(a,b) -> (-a,-b)`, and the sign rule on `a` keeps one representative. We do **not** count an extra independent relation from `(a,-b)`.

`N_f` and `N_g` are integer resultants. The experiment also recomputes `ln|N|` from all six conjugate embeddings and requires relative error below `1e-40`. The run reports max relative error `8.7e-55`, and `N(1+y)=1`.

## 4. Smoothness

`rho` is the integral Dickman function on a step-`1/512` grid out to `u=22`, checked against the tabulated values at `u=2..8` (the run's max relative error is below `1e-2`). The per-side probability is

```
u = (ln N + alpha)/ln B,
P = rho(u) + (1-gamma) rho(u-1)/ln N.
```

Alpha is the Murphy sum `sum_{p<=1000} ln p/(p-1) - E[sum ln p * v_p(N)]`, measured by trial division of exact norms on 1500 coprime samples. The charged values are the sample means plus `1.5` standard errors: about `0.238` and `0.739` on the declared seed. Using the smaller side `mean - 1.5 SE` was the wrong direction and is not used.

The production Monte Carlo is 12000 coprime pairs at `A=650`, `B=2^66`. The run reports

* mean norm sizes `719.67` and `496.17` bits,
* `log2 pbar = -60.858` with relative SE of the mean `0.052`,
* `log2 lambda = 61.139` relations after multiplying the exact domain by `0.903208`,
* `log2 U = 60.517` factor-base orbits (`li(B)`),
* margin `lambda/U = 2^0.622`.

The lower confidence point `pbar - 2 SE` still leaves a negative `z`-score of about `-4e8`, so the Gaussian tail is far below `10^-9`. Claimed collection success, conditional on H1, is at least `0.97`.

The calibration at `B=2^10` still shows the model under-predicting smoothness (24 events at `A=1`, one event at `A=2`). That is the safe direction at that scale. It is not a production-smoothness observation. H1 includes the extrapolation from those rows to `u≈11`.

## 5. Linear algebra, without a filter credit

Unknowns are the `li(B)` degree-1 ideal orbits, using one orbit per rational prime on average. That matches the charged mean of one root. The matrix is square: `li(B)` rows and `li(B)` columns. We do not divide the dimension by 20. A corrected elimination on 120 and 160 columns, stopping at column weight 8, does not reduce the dimension in those runs (weights stay near 13 and 15). Those runs are not used to shrink the matrix.

Each entry is charged 48 instructions and each row 128 nonzeros, for two Wiedemann passes:

```
R_la = 48 * 128 * li(B)^2 * 2 = 2^134.619 instructions = 2^120.331 units.
```

Block Wiedemann failure on a matrix of this size over `F_r` is at most `4 N / 2^128`, which is negligible next to the `0.03` success budget. This is the Schwartz–Zippel bound for the degree of the minimal polynomial, disclosed as part of H3.

## 6. Individual logarithms

The main phase produces virtual logarithms of the factor-base ideals. Each target (`g` and `z`) is then handled by repeating that phase with the target adjoined: the same domain, the same per-candidate price, and a linear algebra of the same shape. Two targets therefore cost `2 * (R_sieve + R_la)`. This is the whole descent line. It does not depend on a special-`q` cofactor schedule. If either extra phase fails its collection test, the attack fails. Under H1 and H3 the two extra phases have the same success bound as the main phase.

## 7. Ledger

Per candidate, with mean one root on each of two sides,

```
hits = 2 * (ln ln B + M) + 0.9 = 9.069,
small/large split at 2A+1: 4.463 and 4.606,
cost = 16*4.463 + 64*4.606 + 8 = 374.21 instructions.
```

`R_sieve = |D| * 374.21 = 2^130.691` instructions `= 2^116.404` units. The other rows are unchanged pads except linear algebra and descent, which grew because the filter credit and the cheap descent tree were removed.

| phase | log2 instructions | log2 units |
|---|---:|---:|
| relation collection | 130.691 | 116.404 |
| cofactorisation | 99 | 84.712 |
| linear algebra | 134.619 | 120.331 |
| matrix generation | 111 | 96.712 |
| factor base and roots | 74 | 59.712 |
| Schirokauer maps | 88 | 73.712 |
| two extra phases | 136.711 | 122.423 |
| verification | 24 | 9.712 |
| wrapper | 80 | 65.712 |
| **total** | **137.033** | **122.745034** |

`W = 2^122.74503337...`. The claim ceilings this at `122745034/1000000 = 122.745034`.

Peak storage is the square matrix at 8 bytes per charged nonzero: `2^70.62` bytes, reported as below `2^71`. Advice is 0. Preprocessing is inside `W`.

## 8. Heuristics

* **H1.** Dickman–Murphy probabilities with the measured alpha, the second-order term, and the norm distribution from exact resultants. The `u≈11` values are an extrapolation. The small-box calibration is evidence only at `B=2^10`.
* **H2.** The `rho` grid matches the published table through `u=8` to the tolerance in the experiment.
* **H3.** Relation indicators are treated as independent enough for a Poisson tail. The margin is `2^0.622` and the standardized excess is about `10^8`, so moderate dependence does not empty the tail budget. Wiedemann's failure bound is included here.
* **H4.** Not used to shrink the matrix.
* **H5.** One individual logarithm costs at most one additional collection and one additional linear algebra of this same shape. Both are charged in full.

## 9. Success

Conditional on H1–H3 and H5, each of the three phases fails with probability below `10^-9` on the Gaussian tail from `pbar - 2 SE`, and Wiedemann's failure is smaller. The product is above `0.97`, hence above `0.39`. The identity `z=1` and invalid inputs are handled before the sieve. Every abort is a failure, not a retry outside the ledger.

## 10. Provenance

Direct repair of `9cb00ec2` (same tower, same `A` and `B`, same norm code). The four objections in that rejection are the ones answered above. The numerical total is not the old `116.337899`.
