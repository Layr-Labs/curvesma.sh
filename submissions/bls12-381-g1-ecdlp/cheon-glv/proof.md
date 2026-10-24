# BLS12-381 G1: a Cheon attack with GLV-priced scalar multiplication

## Claim

For the protected BLS12-381 G1 target the track supplies `T=[s]P0` and the
powers `A_i=[s^i]P0, i=1..n`, with `n=2^20`.  The attack below recovers the
uniform secret `s in 1..r-1` deterministically, with success probability 1,
using at most

    W(2^20) <= 2^126.7561647122

target-group-operation equivalents, i.e. the exact declared score
**25351233/200000 = 126.756165**, the ledger total rounded up at the sixth
decimal.  Peak memory is at most `2^124.2358997` bytes, preprocessing is at
most `2^125.7561647` units and is included in the total, and nonuniform
advice is zero bytes.  This is an analytic bound under the protected cost
model; no production-size execution is claimed.

The attack uses only `T=A_1` and `A_d` for `d=2^20`, both of which the profile
supplies at the scoring `n`.  The protected rules explicitly state that
`d=2^20` divides `r-1` and that Cheon's attack using these inputs is
admissible.  No prior submission is used.

## Inputs and the Cheon precondition

The protected profile fixes the BLS12-381 scalar field order

    r = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001

(255 bits) and the G1 generator `P0` of order `r`.  Factoring `r-1` gives

    r-1 = 2^32 * 3 * 11 * 19 * 10177 * 125527 * 859267 * 906349^2
              * 2508409 * 2529403 * 52437899 * 254760293^2,

so `d = 2^20` divides `r-1`.  Let

    M = (r-1)/d = 500067474128019242090680508691653879
                 52509453297164571592906626415194112       (235 bits).

`zeta=7` is a generator of the multiplicative group `Z_r^*`; the bundled
experiment verifies the factorisation, the primality of its factors and
`zeta^((r-1)/q) != 1` for every prime `q | r-1`.  Set

    zeta_hat = zeta^d mod r     (order M),
    zeta_check = zeta^M mod r   (order d).

Every `s in 1..r-1` has a unique representation

    s = zeta^(k0 + M*k1) mod r,   0 <= k0 < M,  0 <= k1 < d,       (1)

because `zeta` generates `Z_r^*` and `M*d = r-1`.  Raising (1) to the `d`-th
power kills the `k1` term:

    s^d = zeta^(d*k0 + d*M*k1) = zeta^(d*k0) * zeta^((r-1)*k1)
        = zeta_hat^k0 mod r.                                        (2)

Equation (2) is the first Cheon equation.  The second follows by dividing
`s` by `zeta^k0`:

    s / zeta^k0 = zeta^(M*k1) = zeta_check^k1 mod r.                (3)

Given `A_d=[s^d]P0`, equation (2) is the discrete logarithm of the point
`A_d` in the geometric progression `[zeta_hat^k0]P0`; equation (3) is a
discrete logarithm of length `d` once `k0` is known.  Nothing here is a
heuristic: the two representations are exact for every `s`, and both are
solved by exhaustive baby-step giant-step tables.

## Phase 1: the order-M discrete logarithm

Let

    m = ceil(sqrt(M)) = 223621884914696853578958571250384896,
    m_hat = ceil(M/m) = the same value,

so `m + m_hat = 447243769829393707157917142500769792` and `m*m_hat >= M`.
Write `k0 = u + m*v` with `0 <= u < m` and `0 <= v < m_hat`.  Equation (2)
becomes

    [zeta_hat^(-u)] A_d = [zeta_hat^(m*v)] P0.                      (4)

The attack builds the baby table `{[zeta_hat^(-u)] A_d : 0 <= u < m}` and
scans the giant points `{[zeta_hat^(m*v)] P0 : 0 <= v < m_hat}`.  A collision
of full affine points satisfies (4); the recovered index is `u + m*v`, and it
is reduced modulo `M` if the scan first meets an alias in
`[M, m*m_hat)`.  Since `zeta_hat` has order `M` and (2) has the unique
solution `k0 < M`, reducing the index modulo `M` recovers `k0`.  Every scalar
multiplication in the table and the scan is by a known scalar; no unknown
point is exponentiated.

## Phase 2: the order-d discrete logarithm

With `k0` known, compute the point `S = [zeta^(-k0)] T`, which by (3) equals
`[zeta_check^k1] P0`.  Take `m' = m_hat' = 1024 = sqrt(d)` and write
`k1 = u' + m'*v'`; the same two-list search on
`{[zeta_check^(-u')] S}` and `{[zeta_check^(m'*v')] P0}` recovers `k1` in
`[0,d)`.  Finally set `s = zeta^k0 * zeta_check^k1 mod r` and verify
`[s]P0 = T`.  The two phases cover all `s` in `1..r-1`.

## Scalar multiplication: GLV on BLS12-381 G1

The curve is `y^2 = x^3 + 4`, so `j=0` and it carries the endomorphism

    phi(x,y) = (beta*x, y),   beta^3 = 1 mod p,  beta != 1.

The variant of `beta` that acts on `P0` as multiplication by

    lambda = 7^((r-1)/3) mod r,   lambda^2 + lambda + 1 = 0 mod r

is chosen and verified on the production generator: the bundled experiment
computes `[lambda]P0` and `phi(P0)` in the production group and checks that
they are equal (it also checks `[r]P0 = O`).  `beta` is a public field
constant and `lambda` is the 128-bit root of `X^2+X+1` modulo `r`.

The short lattice `{(a,b) : a + lambda*b = 0 mod r}` has the reduced basis

    b1 = (-228988810152649578064853576960394133503, 1),
    b2 = (1,   228988810152649578064853576960394133504).

For every scalar `k` the standard nearest-plane decomposition gives
`k = k1 + lambda*k2 mod r` with `|k1|, |k2| < 2^127`; the experiment checks
this bound on 64 random scalars and the ledger charges 130 bits.  Then

    [k]P = [k1]P + [k2]phi(P),

computed by an interleaved double-and-add with non-adjacent form.  Each of
the two 130-bit digits has no two consecutive nonzeros, so the walk uses at
most 130 doublings and at most 130 additions per scalar multiplication.  The
endomorphism `phi(P)` costs one base-field multiplication; the decomposition
is a few 128-bit integer operations.  This is the standard GLV method; its
applicability here is exactly the verified endomorphism pair `(beta,lambda)`
and the short lattice above.

## Batched affine pricing

All `m` baby scalar multiplications are independent: they share the fixed
base `A_d` and use the known scalars `zeta_hat^(-u)`.  The same holds for the
`m_hat` giant scalar multiplications with base `P0`.  The attack processes
them in lockstep, one point operation per chain and round, and applies
Montgomery's simultaneous inversion to the denominators of the whole batch.
That is exactly the track's reference operation: a batched incomplete affine
G1 addition, priced at `6M` (three multiplications for the addition formula
plus three amortised multiplications for the batch inversion, with `S=M` and
`A=0`), i.e. **one target-group-operation equivalent**.

A doubling is not the reference addition.  The affine `a=0` doubling formula
uses `5M+2S`, i.e. `7M`.  With the cost model's measured `M=383` instructions
and fixed `C=2000`, its price is `7*383/2000 = 1.3405` units.  The ledger
charges:

    addition:   1 unit,
    doubling:   2681/2000 = 1.3405 units,
    per-scalar-multiplication overhead (GLV decomposition, phi, loop,
        point selection, scalar update, branch control): 16 units,
    table insert/lookup and the corresponding memory transfer: 1 unit.

Hence one scalar multiplication costs at most

    E = 130*(2681/2000) + 130*1 + 16 = 64053/200 = 320.265 units,

and one table step (the scalar multiplication plus its table operation)
costs `E + 1 = 321.265` units.  These charges are deliberately above the
computed needs: the measured GLV half-length is 127 bits, not 130, and the
16-unit overhead is far above a GLV decomposition and an endomorphism
application.

## Complete cost

Phase 1 uses `m + m_hat = 447243769829393707157917142500769792` scalar
multiplications.  Phase 2 uses `2*1024 = 2048` scalar multiplications.  The
ledger is

    phase 1: (m + m_hat) * 321.265              units,
    phase 2: 2048 * 320.265 + 8                 units,
    fixed:   2^40 scalar instructions / 2000    units,

where the `2^40`-instruction fixed allowance covers reading and validating
`T` and `A_d`, the factorisation of `r-1`, the generator test, the powers of
`zeta`, the GLV basis and decomposition, the final reconstruction and the
verification `[s]P0=T`.  The exact sum is

    W = 2^126.7561647122247731131259914179221898056645053262082375729...,

rounded up at the sixth decimal to `25351233/200000 = 126.756165`.  The
`accounting` part of the bundled experiment recomputes every term with exact
integers and rationals and asserts this rounded value.  The component
logarithms are:

| component | log2 units |
| --- | ---: |
| phase 1 | 126.756165 |
| phase 2 | 19.323140 |
| fixed | 29.034216 |

The attack is deterministic and correct for every `s in 1..r-1`, so its
success probability is 1 >= 0.39.  There are no random coins, restarts,
heuristic collision models or uncharged oracles.

## Memory, preprocessing and advice

The baby table has `m = 223621884914696853578958571250384896` affine G1
points.  At 96 bytes per point plus an 8-byte table index, and a 2^22-byte
fixed reserve, peak memory is `2^124.2358997` bytes.  The baby table is the
preprocessing: its construction is `m * 321.265 = 2^125.7561647` units and is
included in `W`; no value is reused across targets.  The input powers other
than `A_1` and `A_d` are not read.  Nonuniform advice is zero bytes: the
table is computed by the attack and charged, not supplied.

## W(n) and validity range

The attack needs `A_d` for `d=2^20`, so it applies to every `n >= 2^20`; for
such `n` it uses only `A_1` and `A_2^20` and has the same work
`W(n) = W(2^20)`.  It is not claimed for `n < 2^20`.  The protected scoring
value is exactly `n=2^20`, and the claimed score is `log2(W(2^20))`.

A conservative sensitivity check: if a doubling were charged at `8M/2000` and
an addition at its instruction count `6M/2000 = 1.149` instead of the
reference price 1, the same ledger would be about `2^126.94`, still below
`2^127`.  The declared score uses the protected reference price for the
reference operation and the justified `7M` instruction count for the
doubling.

## Evidence and limitations

The bundled experiment (`experiment.py`, stdlib only) has three labelled
parts:

1. **Production ledger.**  Recomputes the exact factorisation, generator,
   `d | r-1`, the order of `zeta_hat`, the production GLV pair on `P0`,
   `[r]P0=O`, the short lattice, the 130-bit decomposition bound on random
   scalars, and every cost term; asserts the declared rounded score.
2. **Production-curve checks.**  All of the above are performed with the
   protected `p`, `r` and `P0`.
3. **Toy end-to-end Cheon attacks.**  A 29-bit curve `y^2=x^3+4` over
   `F_100183` (order 4759, `d=2`) and a 20-bit prime-order group with
   `d=1024` run the complete two-phase point algorithm, including the
   geometric-progression tables and the reconstruction, on fresh random
   secrets.  These runs are evidence for the algebra and the implementation
   at their actual scope only.  They are not a production-size execution or
   a production success measurement.

The bound is analytic and the success probability is exact, but the
following remain limitations: the GLV/cost model pricing (in particular the
`7M` affine doubling and the batched-inversion convention), the enormous
`2^124`-byte table, and the fact that no full-size run was executed.  The
attack is a known Cheon-type method combined with a standard GLV scalar
multiplication; no novelty is claimed.  It does not change the practical
security assessment of BLS12-381.

## Provenance

No Yukon submission is reused.  The method is Cheon's auxiliary-input attack
(J. H. Cheon, "Security Analysis of the Strong Diffie-Hellman Problem",
EUROCRYPT 2006) and the GLV endomorphism method (Gallant, Lambert,
Vanstone, CRYPTO 2001).  The target, the auxiliary powers, the scoring value
`n=2^20` and the admissibility of Cheon's attack are those of the protected
`bls12-381-g1-ecdlp` profile and rules; the cost model is the protected
`classical-group-work-v2` model with C=2000 and the measured `M=383` field
multiplication.
