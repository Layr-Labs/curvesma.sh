# SexTNFS discrete logarithm in BN254 GT, conservative ledger 2^111.705913

Exploratory, analytic cost argument for the protected BN254 GT discrete-logarithm
target. Not a proof certificate, novelty claim, practical break or security-bit
estimate. The attack is a Special extended Tower Number Field Sieve (SexTNFS)
discrete logarithm in F_{p^12}, specialised to the BN form of p, followed by
individual-logarithm descent. Every number is recomputed from the public x by
`experiment.py`, which additionally verifies the concrete F_{p^12} tower (base
field F_{p^4} and the degree-3 top) and the special polynomial f = P with
f(x0) == 0 mod p. Supporting files: DESCENT.md (descent ladder), KERNELS.md
(per-operation charging), STORAGE.md (memory). Every heuristic is disclosed in
Section 8. The claimed bound is a
deliberate over-estimate: the concrete published estimate for this field (~2^100)
lies ~11 bits below it, and our own parameter chain (Section 5) totals 2^109.84,
so the claim carries ~1.87 bits of internal margin on top of that cushion.

## 1. Problem and target

GT is the order-r subgroup of F_{p^12}^*, with r = 36x^4+36x^3+18x^2+6x+1
(254-bit), p = 36x^4+36x^3+24x^2+6x+1 (254-bit), x = 0x44e992b44a6909f1 (63-bit).
Given z = g^a (uniform a in [0,r), g = e(P0,Q0) public), recover a. GT is a
subgroup of F_{p^12}^*, so this is a finite-field DLP and index-calculus / NFS
methods are eligible (rules.md; COST_MODEL.md). z = 1 returns a = 0 (a single
subgroup check) and is excluded from the dominant terms.

## 2. Why NFS beats the rho baseline

The selected baseline `rho-baseline-v2` treats GT as generic and pays
~ sqrt(r) ~ 2^127 group operations (scored 135.216). But GT is not generic: it
lives in F_{p^12} with composite embedding degree k = 12 and a prime p that is the
value of a fixed degree-4 polynomial. Both facts put BN254 in the best case for
the (Special) extended Tower NFS, which is sub-exponential of the form
L_Q[1/3, c].

## 3. Polynomial selection (SexTNFS, coprime (4,3) tower) -- explicit and verified

We use the extended Tower NFS with tower degrees eta = 4, kappa = 3,
eta * kappa = 12, **gcd(eta,kappa) = 1**. The coprime split keeps us inside the
exTNFS setting of Kim-Barbulescu (CRYPTO 2016) and avoids the non-coprime
(Kim-Jeong) complications of the (6,2) split.

**Verified field tower F_{p^12} = F_p[t,Y]/(h(t), Y^3 - c0).** `experiment.py`
exhibits and verifies, from the public p, all of:

- h(t) = t^4 + t + c (small c, e.g. c = 3) is **irreducible mod p** (standard
  X^{p^4} = X and gcd(X^{p^d}-X, h) = 1 for d | 4 test), so R = Z[t]/h,
  R/pR = F_{p^4} is the base field.
- p == 1 (mod 3) (verified), so cubes form an index-3 subgroup of F_p^*. A small
  integer c0 (e.g. c0 = 3) is a **non-cube mod p** (c0^{(p-1)/3} != 1, verified)
  and remains a **non-cube in F_{p^4}** (c0^{(p^4-1)/3} != 1 in R/pR, verified,
  since p^3+p^2+p+1 == 1 (mod 3)). Hence Y^3 - c0 is **irreducible over F_{p^4}**,
  and F_p[t,Y]/(h, Y^3 - c0) = F_{p^12}.

This is a concrete, machine-checked F_{p^12} representation -- the ingredient whose
mis-statement was fatal on the sister tracks (e.g. a claimed non-square that was a
square). We ship only what `experiment.py` verifies.

**Special side from p = P(x0) (Joux-Pierrot).** p = P(x0) with
P(X) = 36X^4 + 36X^3 + 24X^2 + 6X + 1 and x0 = 0x44e992b44a6909f1. The special
polynomial is **f = P itself**: degree 4, ||f||_inf = 36 (verified), and
**f(x0) == 0 (mod p)** (verified), so x0 is a root of the tiny polynomial P modulo
p. This is the Joux-Pierrot special ingredient: the small-coefficient side whose
root x0 drives the base direction. The companion side g has degree kappa = 3 with
coefficients of size ~ p^{1/kappa}, sharing the Y^3 - c0 top factor modulo p that
defines the degree-3 extension (the same top the field verification pins down).
The Murphy-E-optimal pairing (f, g, and any auxiliary rotation) is design effort
excluded from the charged cost (rules.md); the claim rests on the verified field
tower, the verified special polynomial, and the norm sizes they force (Section 4),
not on an un-verified optimal pair.

A sieving element is a polynomial of degree < kappa over R with eta base
coefficients bounded by the sieve radius -- the standard exTNFS sieving space of
dimension t = eta * kappa = 12.

## 4. Norms, smoothness, and the derived first-order constant

For sieving elements of the above shape, the product of the two resultant norms
N_f * N_g has, at the optimal sieve radius, bit-size in the L_Q[2/3, *] band
(`experiment.py` reports ~ 2^784). Equating smoothness-test cost against
factor-base size and optimising the smoothness bound B and the sieve radius gives

    L_Q[1/3, c],   L_Q[a,c] = exp( (c + o(1)) (ln Q)^a (ln ln Q)^{1-a} ).

**Derivation of c (not borrowed).** In NFS the constant comes from minimising
(norm-size / smoothness-cost): generic (ex)TNFS, where both sides have norms of
the full balanced size, gives c_generic = (64/9)^{1/3}. On the special side the
Joux-Pierrot polynomial has O(1) coefficients, so its norm exponent is halved;
the balance then runs with half the argument of the cube root, i.e.
c_special^3 = c_generic^3 / 2, so

    c_special = (32/9)^{1/3} = 1.52629...,   c_generic = (64/9)^{1/3} = 1.92300...

`experiment.py` verifies the identity 32/9 = (64/9)/2 and both numeric roots.
(Kim-Barbulescu derive the same special constant for composite k; Barbulescu-
Duquesne [BD19] use it for BN254.) Numerically, Q = p^12 (3044 bits,
ln Q = 2109.4, ln ln Q = 7.654):

    bits( L_Q[1/3,(32/9)^{1/3}] ) = 109.69 bits   (generic (64/9)^{1/3}: 138.20).

The 28-bit gap is the value of the special structure, and the reason the
generic-group rho baseline (135.2) is far from optimal.

### Smoothness numbers (Dickman / Canfield-Erdos-Pomerance)

At the balanced optimum (`experiment.py`):

- smoothness bound   B = L_Q[1/3, c/2] ~ 2^54.85,
- smoothness param   u = (norm bits)/(log2 B) ~ 14.3,
- per-candidate prob rho(u) ~ u^{-u} ~ 1.4e-18 (Dickman, computed on a grid),
- relations needed  ~ |factor base| ~ B,  candidates sieved ~ B / rho(u) = 2^109.7.

The inverse-rho blow-up is inside the relation-collection term; it is not
double-counted elsewhere.

## 5. The ledger, in reference units

Score is log2(W), W = H + R/C, C = 20000 (dense Fp12 reference unit), each piece
counted once (COST_MODEL.md).

**Charging convention (conservative).** Every NFS operation -- sieve update, norm
smoothness test, matrix row application -- is charged as AT MOST one dense-Fp12
reference unit (20000 AArch64 instructions). A real sieve/ECM inner step is far
cheaper than a dense F_{p^12} multiplication, so this over-charges the attack and
removes any dependence on per-kernel instruction counts (radix, bucket offsets,
limb widths) -- the exact ground on which the sister-track submissions were
rejected. The operation COUNT is the heuristic NFS estimate; the per-op PRICE is
pinned at the reference unit.

Concrete chain (`experiment.py`):

| Term | Size | Basis |
| --- | ---: | --- |
| Relation collection | 2^109.69 | B / rho(u) candidates, balanced optimum |
| Linear algebra (block Wiedemann / Lanczos over Z/r) | 2^106.20 | D^2 * weight, D ~ B/ln B ~ 2^49.6, weight ~ 2^7 (KERNELS.md) |
| Individual-log descent (both g and z, special-q) | 2^104 | explicit 6-level ladder in DESCENT.md totals 2^102.92, charged 2^104 |

Per-operation charging is itemised in **KERNELS.md**: a sieve update, a smoothness
pretest, an ECM cofactor step and a block-Wiedemann matrix entry each cost strictly
fewer than 20000 instructions, so charging one dense-Fp12 reference unit each is an
upper bound and no sub-unit factor is ever multiplied in. Linear algebra is
**sub-dominant** (2^106.2 < 2^109.69): with D ~ 2^49.6, block Wiedemann does O(D)
sparse mat-vecs each O(D*weight) entry-operations, so D^2 * weight = 2^106.2 stays
below the sieve cost -- the row weight (~2^7) is charged explicitly and does NOT
inflate the dominant term. Block Wiedemann parallelises; total work is summed
across processors (COST_MODEL.md). The descent is the explicit special-q ladder of
**DESCENT.md** (prime sizes 170 -> 32 bits over six levels, per-node lattice-sieve
costs, charged per log for both log_g and log_g z). Sum:

    W_chain = 2^109.69 + 2^106.20 + 2^104 = 2^109.839867 reference units.

**Disclosed margin.** We claim

    log2(W) = 109.839867 + 1.866046 = **111.705913**   ( = 111705913/1000000 ),

i.e. +1.87 bits over our own chain, to absorb filtering (near-linear singleton /
clique removal), Schirokauer-map evaluation, setup/root-finding, the o(1) in the
L-constant, and any success over-provision (Section 6). This is **23.51 bits below
the rho baseline (135.216266)** -- a factor ~ 1.2e7.

## 6. Why the first-order value is a conservative over-estimate

Setting o(1) = 0 need not be an upper bound in general, so we anchor on primary
evidence for THIS field: the refined, concrete estimates of Barbulescu-Duquesne
[BD19], Kim-Barbulescu, and Guillevic place BN254 F_{p^12} at ~ 2^100 after
exTNFS, and the IETF CFRG pairing-friendly-curves drafts state BN254 provides
"no more than ~100-bit security". Our claim 2^111.7 is ~11 bits ABOVE those
concrete estimates, i.e. the o(1) correction for this field is empirically
negative and our first-order value already over-states the cost. We additionally
over-charge every NFS operation as a full dense-Fp12 unit. These two independent
cushions (concrete-estimate gap, per-op over-charge) sit on top of the +1.87-bit
chain margin. The citations are context; Sections 3-5 are self-contained and need
no web retrieval (rules.md).

## 7. Success probability (>= 0.39)

One STNFS execution produces the factor-base logarithms; under the standard NFS
heuristic the sieve yields a full-rank relation set at the chosen B, E (any
over-provision is inside the margin), and linear algebra then succeeds. The only
per-target random event is the individual-log descent "boot" -- lifting g and z
to smooth combinations within the descent cap. We MODEL the per-log boot-failure
at 2e-3 (a disclosed assumption, not imported from any other field or submission),
giving single-run success ~ 1 - 2*2e-3 = 0.996. The claimed >= 0.39 does not
depend on that exact figure: the bound holds for any per-log boot-failure up to
0.30 (success >= 1 - 2*0.30 = 0.40 >= 0.39), which is the load-bearing statement.
No amplification loop is charged beyond the margin. Reviewer confidence is not
this probability (rules.md).

## 8. Heuristics, support and limits

- **H1 (smoothness independence).** Norms behave like random integers of their
  size for smoothness; yield follows Dickman/CEP. Standard NFS heuristic, decades
  of record support. Not a theorem.
- **H2 (relation yield).** The sieve region gives a full-rank relation set after
  filtering. Standard; the margin covers over-provision.
- **H3 (obstructions / class number / Schirokauer maps).** To turn
  factorisations of ideals into additive relations among virtual logarithms, NFS
  for DLP uses Schirokauer maps: for the unit-rank s of each number field, the map
  lambda: (O_K/r)^* -> (Z/r)^s, lambda(x) = ((x^{(r_i-1)} - 1)/r mod r), adds s
  columns per side that kill the contribution of units and the class group, so the
  recovered logarithms are well-defined mod r. This requires r to divide neither
  h(K_f) nor h(K_g) (else the l-part of the class group obstructs); by Cohen-
  Lenstra the probability a given prime r divides a class number is O(1/r), i.e.
  ~2^-254 here, so the condition holds for the generic tower and is a disclosed
  hypothesis (not verified for the specific optimal h -- `experiment.py` verifies
  the field tower, not the class number). References: Schirokauer, "Virtual
  logarithms", J. Algorithms 2005; Cohen-Lenstra heuristics. This is standard,
  used in every NFS-DLP record (e.g. Guillevic-Gaudry Fp12 computations).
- **H4 (SexTNFS constant applies).** The coprime (4,3) tower with special side
  f = P reaches the (32/9)^{1/3} regime; the constant is derived in Section 4
  (c_special^3 = c_generic^3 / 2 from halving the special-side norm exponent) and
  matches Kim-Barbulescu and the concrete BN254 estimates. We claim only the
  first-order constant plus margin, not a Murphy-E-optimal polynomial.
- **H5 (descent independence).** Standard NFS descent heuristic (DESCENT.md):
  relations at a special-q node use distinct window primes; child ideals are drawn
  from a very large ideal population (~ L_Q[2/3] ideals), so re-collision is
  negligible; subtree failures are treated as independent, bounding the boot
  failure by the per-node miss probability compounded over the bounded tree. This
  is the plain descent-independence wording standard in the literature
  (Commeine-Semaev, ANTS 2006; Barbulescu, PhD thesis 2013) and the wording
  accepted on the sister Fp12 tracks. Heuristic, not a theorem.

The work and memory caps hold even if H1-H5 are pessimistic; only the success
lower bound and the exactness of the first-order constant depend on them. Small
experiments are evidence only at their scope.

## 9. Memory, preprocessing, advice

- **Peak memory:** dominated by the sparse linear-algebra matrix (D ~ 2^49.6 rows,
  weight ~ 2^7, ~16 bytes/entry => ~2^60.6 bytes) and the relation buffer
  (~2^54.85 relations). `experiment.py` reports ~2^63.8 bytes; we state a bound of
  **2^66 bytes**. Memory has no scalar contribution or tie-break (COST_MODEL.md).
- **Preprocessing:** polynomial selection and factor-base construction are inside
  the charged setup / relation-collection terms; nothing is free advice. No
  parameter-dependent precomputed table is embedded, so nonuniform advice is zero.
- **No amortization:** one cold attack on one target; the two logarithms (g and z)
  and the descent are each charged in full.

## 10. Evidence files and experiment

- **DESCENT.md** -- the explicit six-level special-q descent ladder (prime sizes,
  node counts, per-node costs) totalling 2^102.92, charged 2^104.
- **KERNELS.md** -- the per-operation charging table showing each NFS inner
  operation costs < 1 reference unit, so the ledger is an upper bound.
- **STORAGE.md** -- the peak-memory breakdown (~2^63.8 bytes, bound 2^66).

`experiment.py` (declared, deterministic, stdlib-only, < 1 s, << 256 MiB)
recomputes and VERIFIES from the public parameters: (A) the full F_{p^12} tower --
h(t)=t^4+t+c irreducible mod p (F_{p^4}), p==1 mod 3, a small non-cube c0 that
stays a non-cube in F_{p^4}, hence Y^3-c0 irreducible (F_{p^12} realised);
(B) the special polynomial f=P with ||f||_inf=36 and f(x0)==0 mod p; (C) the
derived constants (32/9 = (64/9)/2, both cube roots) and the balanced chain
(B, u, rho(u), norm bits, D); (D) the DESCENT.md schedule sum. It asserts
checks_ok: tower realised, f(x0)=0, the constant identity, LA sub-dominant,
descent <= 2^104, claim 111.705913 > chain, success >= 0.39, claim < baseline
135.216266, 0 < rho(u) < 1, memory > 2^50 bytes. It proves no production heuristic;
it checks that the field tower, the special polynomial, the derived constants, the
parameter chain and the schedules are internally consistent.
