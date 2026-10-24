# Argument (v20)

1. *Reduction* (unchanged). `g = e(P0,Q)`, `z = e(P,Q) = g^x`; a walk in GT with
   `y = g^A z^D`; a class collision with distinct labels gives one linear equation in `x`:
   `x = (p^k A1 - A2)/(D2 - p^k D1)`, `k` from the Frobenius shift between the two states.
2. *Class function* (unchanged). `kappa` from the minimal cyclic rotation of the twelve
   keys, `s = idx(H(m0))`, multiplier `M_s^(p^-kappa)`.
3. *Repair 1*. The twelve keys are linear forms of `d0 = y[1][0]`; `d0(y^2)` is not a
   function of `d0(y)`, so the keys are recomputed on squaring steps and charged
   (1802 + 144). Measured refutation of both candidate cheap rules and a witness that
   `keys(y^2)` is not a function of `keys(y)`: `evidence/keys20.json`.
4. *Repair 2*. All three label fields come from `L = mix128(m0)`: branch bits 0-2, index
   bits 3-14, DP bits 15+. Trails end at their first report; the pool stores
   (fingerprint, trail, step); recovery replays both trails. Census:
   `evidence/sparse20.json`.
5. *Repair 3*. Coalescence under independent random previous slots, plus a mutation
   control that fails when the guard reads `prev` in a non-class-consistent way; the
   residual merge-divergence rate 2/R is priced into rho: `evidence/coalesce20.json`.
6. *Trap census*. 0 length-2/3/4 returns per 1600 traces x
   200 steps under the declared step, against
   114 per 1600 under the guard-free control:
   `evidence/trap20.json`. Longer balanced products are bounded analytically (2^-105).
7. *Work*. Claim 125.8259 with cap `floor(2^125.8259)` covering the nominal work;
   rows and exact rational W in `evidence/ledger20.json`.
8. *Success*. Certified floor 0.390100 >= 0.39 after abort, DP-detection and storage
   debits, with the census retention `0.950000` applied as a discount.
