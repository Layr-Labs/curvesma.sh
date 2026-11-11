# SexTNFS discrete logarithm in BN254 GT -- concrete finite-regime ledger 2^123.304316

Exploratory, analytic cost argument for the protected BN254 GT discrete-logarithm
target. Not a proof certificate, novelty claim, practical break or security-bit
estimate. The attack is a Special extended Tower Number Field Sieve (SexTNFS) in
F_{p^12}, specialised to the BN form of p, followed by individual-logarithm
descent. Every derived number is recomputed by `experiment.py`, which also exhibits
and verifies the concrete polynomial triple and computes the two sieve norms by
exact resultants on real sieve elements. Supporting files: DESCENT.md, KERNELS.md,
STORAGE.md.

This is a deliberately conservative FINITE-REGIME ledger. It does not borrow an
asymptotic L-value as a cost; it computes the candidate count directly from
explicitly computed norms and a self-checked Dickman rho, charges a full reference
unit per sieve candidate (a disclosed ~7.8-bit over-charge against the honest
instruction count), and adds a stated 4.5-bit margin. The honest instruction-level
cost of the same attack is ~2^117.4; concrete published estimates for BN254 are
near 2^100. The claim 2^123.304316 sits above both, 11.9 bits below the selected
Pollard-rho baseline 135.216266.

## 1. Problem and target

GT is the order-r subgroup of F_{p^12}^*, with r = 36x^4+36x^3+18x^2+6x+1
(254-bit), p = 36x^4+36x^3+24x^2+6x+1 (254-bit), x = x0 = 0x44e992b44a6909f1
(63-bit). Given z = g^a (uniform a in [0,r), g = e(P0,Q0) public), recover a. GT is
a subgroup of F_{p^12}^*, so this is a finite-field DLP and index-calculus / NFS
methods are eligible (rules.md; COST_MODEL.md). z = 1 returns a = 0 by a single
subgroup check and is excluded from the dominant terms. Q := p^12 has 3044 bits.

## 2. Why NFS beats the rho baseline

The selected baseline treats GT as generic and pays ~ sqrt(r) ~ 2^127 group
operations (scored 135.216266). GT is not generic: it lives in F_{p^12} with
composite embedding degree k = 12 and a prime p that is the value of the fixed
degree-4 BN polynomial P at x0. Both facts put BN254 in the special case of the
extended Tower NFS. Note that *generic* exTNFS (constant (64/9)^{1/3}) is already
near 2^138 at first order and does NOT beat the baseline; the special structure
(p = P(x0)) is load-bearing, and the ledger below is built from an explicit
special polynomial pair, not from an asymptotic constant.

## 3. Explicit, verified polynomial selection (coprime (4,3) tower)

We use the extended Tower NFS with tower degrees eta = 4, kappa = 3,
eta*kappa = 12, gcd(eta,kappa) = 1 (Kim-Barbulescu, CRYPTO 2016). `experiment.py`
exhibits and verifies, from the public p, the full triple (h, g, f):

- **Base tower h.** h(t) = t^4 + t + 3 is **irreducible mod p** (checked by
  t^{p^4} = t and gcd(t^{p^2}-t, h) = 1), so R = Z[t]/h and R/pR = F_{p^4}.
- **Companion g (Joux-Pierrot special form, S = 0).** g(X) = X^3 - x0. Since
  p = P(x0) = 0 mod p, and p == 1 (mod 3) (verified), we check that x0 is a
  **non-cube in F_p** and remains a **non-cube in F_{p^4}**
  (x0^{(p^4-1)/3} != 1 in R/pR, verified). Hence g = X^3 - x0 is **irreducible
  over F_{p^4}**, and F_p[t,X]/(h, X^3 - x0) = **F_{p^12}**. This is the concrete,
  machine-checked degree-12 field; g is the explicit companion polynomial that was
  missing from v1.
- **Special side f.** f(X) = P(X^3) = 36 X^12 + 36 X^9 + 24 X^6 + 6 X^3 + 1, of
  degree 4*kappa = 12 and infinity-norm 36. Modulo (g, p) we have X^3 == x0, so
  f == P(x0) == 0: **g divides f mod p** (verified). f and g therefore share the
  degree-3 factor g over the base F_{p^4}; combined with the degree-4 tower h this
  is exactly the common degree-12 structure that defines the target field. (f(x0)
  itself is non-zero for kappa > 1; the load-bearing identity is f == 0 mod (g,p),
  which is what `experiment.py` verifies.)

A sieving element a(t,X) is a polynomial of degree < kappa = 3 (sieving degree
tau; here tau = 3) over R with eta = 4 base coefficients each bounded by the sieve
radius A, i.e. n = eta*tau = 12 integer coefficients.

## 4. Norms from the actual polynomials, and smoothness

For a sieving element a, the two norms are the exact integer resultants

    N_f(a) = | Res_t( Res_X(a, f), h ) |,   N_g(a) = | Res_t( Res_X(a, g), h ) |.

`experiment.py` computes these as big integers by fraction-free (Bareiss)
elimination over Z[t] followed by an integer resultant with h, on 48 real random
sieve elements at the production radius A = 2^9. The sampled bit-sizes are

    Nf ~ 522 bits (max 552),   Ng ~ 612 bits (max 622).

These match the Kim-Barbulescu bound log2 N_f ~ eta*deg_f*log2 A + ... : the
measured slope in log2 A is 48 = eta*deg_f on the f-side and 12 = eta*deg_g on the
g-side, confirming the degrees. The special benefit is visible: the g-side carries
the size through ||g|| = x0 ~ 2^63 (= p^{1/4}), not p^{1/2} (generic) or p^{1/3}.

At the smoothness bound B = 2^62.4 (chosen offline to balance relation collection
against linear algebra; `experiment.py` recomputes the balance):

    u_f = log2 Nf / log2 B = 8.37,   u_g = log2 Ng / log2 B = 9.81,
    P_smooth = rho(u_f) * rho(u_g) = 2^-60.16   (two sides independent, heuristic H1).

rho is the Dickman function, built by the stable integral recurrence
u*rho(u) = int_{u-1}^{u} rho; `experiment.py` self-checks it (two step sizes agree
to 1.2%, the ratio rho(u+1)/rho(u) is monotone) instead of trusting a table.
P_smooth is the MEAN of rho(u_f)*rho(u_g) over the real samples; the per-sample max
norms and the resulting spread are inside the margin.

## 5. The finite-regime ledger (reference units, C = 20000)

Score is log2(W), W = H + R/C, each piece counted once (COST_MODEL.md).

**Relation collection.** Factor base size |FB| ~ 2*B/ln B = 2^57.97 primes. The
number of sieve candidates that must be tested to obtain a full-rank relation set
is

    candidates = |FB| / P_smooth = 2^57.97 / 2^-60.16 = 2^118.12.

This is the direct (relation-count)/(smoothness-probability) quantity the v1
review asked for; it is NOT an asymptotic L-value. Each candidate is charged a
FULL dense-Fp12 reference unit (20000 instructions). The honest per-candidate
instruction count is ~88 (a few Mertens sieve updates on each side, array init,
cofactor bookkeeping; KERNELS.md), so the honest relation-collection cost is
2^110.30 and the full-unit charge over-charges by +7.83 bits. We claim the
conservative value:

    relation collection = 2^118.12.

**Linear algebra (block Wiedemann over Z/r), instruction-priced.** The sparse
matrix has dimension D = |FB| = 2^57.97. Block Wiedemann does ~3D sparse
matrix-vector products; each nonzero is a small-integer x 254-bit multiply-add
(~20 instructions), and each side contributes a Schirokauer block of dense 254-bit
columns. We derive the average row weight ~18 from (log Nf + log Ng)/log(sqrt B)
and charge 60 Schirokauer columns (unit ranks up to eta*deg_f and eta*deg_g, an
over-count), each a 254-bit modmul (~300 instructions). Per matrix-vector-entry
cost is then weight*20 + schir*300 = 18364 instructions, and

    LA = 3 * D^2 * 18364 / C = 2^117.39.

LA is sub-dominant to relation collection; the Schirokauer columns, not the sparse
weight, dominate it. (This is honest instruction-level pricing; no "1 unit per
operation" inflation is applied to LA.)

**Individual-logarithm descent (both log_g and log_g z).** DESCENT.md gives a
capped Las-Vegas special-q ladder. Each target is first randomised to an algebraic
integer of norm ~ L_Q[2/3, 0.95] ~ 2^444, then a bounded special-q tree drives the
large primes down to the factor base; each level's success probability is rho of
the level norm over log2 B (same Dickman), and each node is run with a cap
c/P_level, c = 40, charged a full reference unit per candidate. The ladder totals
2^26 per log and 2^27 for both; it is sub-dominant by >90 bits, so the chain is
insensitive to it (even charging descent at 2^115 moves the chain < 0.5 bit, inside
the margin).

**Chain.**

    W_chain = 2^118.12 + 2^117.39 + 2^27 = 2^118.804.

## 6. Claim and margin

    log2(W) = 118.804 + 4.500 (disclosed margin) = 123.304316 = 123304316/1000000.

The 4.5-bit margin absorbs: the asymptotic o(1) and the variance of the sampled
norms around their mean (a few tenths in u shifts P_smooth by ~2 bits); filtering
(singleton/clique removal, a few near-linear passes over ~2^58 relations);
Schirokauer-map evaluation per relation; polynomial-selection / root-finding setup;
and any success amplification. These are individually far below the dominant term;
the margin is a single conservative envelope. The much larger ~7.8-bit per-candidate
over-charge is a SEPARATE cushion already inside the chain. The claim is
**11.91 bits below the Pollard-rho baseline 135.216266** and above both the honest
cost of this attack (~2^117.4) and the concrete published estimates for this field
(~2^100).

## 7. Success probability (>= 0.39)

One SexTNFS execution produces the factor-base logarithms; under the standard NFS
heuristic the sieve yields a full-rank relation set at the chosen B (over-provision
is inside the margin), and linear algebra then succeeds. The only per-target random
event is the descent. Each descent node is a capped Las-Vegas loop with cap
c/P_level, c = 40, so its failure probability is at most e^-40; with ~2^13 nodes
the expected failed-node count per log is ~2^13 * e^-40 ~ 10^-14, and success per
target is ~ 1 - 2*10^-14. The >= 0.39 bound does not depend on that exact figure:
it holds for any per-log boot-failure up to 0.30 (success >= 1 - 2*0.30 = 0.40).
Reviewer confidence is not this probability (rules.md).

## 8. Heuristics, support and limits

- **H1 (smoothness).** The two resultant norms behave like independent random
  integers of their (computed) sizes for smoothness; yield follows Dickman/CEP.
  Standard NFS heuristic; here the norm sizes are computed exactly, not assumed.
- **H2 (relation yield).** The candidate region gives a full-rank relation set
  after filtering; over-provision is inside the margin.
- **H3 (Schirokauer / class number).** Virtual-logarithm relations use Schirokauer
  maps; r divides neither class number with probability O(1/r) ~ 2^-254
  (Cohen-Lenstra), a disclosed hypothesis on the specific tower (the field tower is
  verified; the class number is not). Refs: Schirokauer 2005; Cohen-Lenstra.
- **H4 (special construction).** The explicit coprime (4,3) tower with f = P(X^3)
  of norm-free small coefficients and g = X^3 - x0 carrying size x0 = p^{1/4} is
  verified to realise F_{p^12}; the cost is computed from its actual norms, not
  from an asserted constant.
- **H5 (descent independence).** Standard NFS descent heuristic (Commeine-Semaev
  2006; Barbulescu 2013): child ideals are drawn from a very large population so
  re-collision is negligible and subtree failures are treated as independent.

The work and memory caps hold even if H1-H5 are pessimistic; only the success lower
bound and the exactness of the asymptotic constant depend on them. Small
experiments are evidence only at their scope: the norm computation is exact at the
production radius, but the full attack is not executed.

## 9. Memory, preprocessing, advice

Peak storage is dominated by the relation buffer (~|FB| = 2^57.97 relations, few
hundred bytes each, ~2^66 bytes) and the sparse matrix (D = 2^57.97, weight ~18,
~16 B/entry, ~2^66.3 bytes); `experiment.py` reports ~2^67.1 bytes and we state a
bound of **2^70 bytes**. Memory has no scalar contribution (COST_MODEL.md).
Polynomial selection and factor-base construction are inside the charged setup /
relation terms; no parameter-dependent precomputed table is embedded, so nonuniform
advice is zero. One cold attack on one target; both logarithms and the descent are
charged in full, with no amortization.

## 10. Experiment

`experiment.py` (declared, deterministic for the claim, stdlib-only, ~2 s,
<< 256 MiB) constructs and verifies (h, g, f), computes the sieve norms by exact
resultants on real samples, the self-checked Dickman rho, the factor base,
candidate count, both relation-collection costs, the linear-algebra cost, the
descent ladder and success bound, the chain, the disclosed margin, the claimed
rational 123304316/1000000, the storage bound, and an independent cross-check with
the caller's seed (asserted within 1 bit of the claim chain; the claim itself uses
a fixed internal seed and does not depend on the caller's seed). Its JSON ends with
`"checks_ok": true`. It proves no production heuristic; it checks that the field
tower, the explicit polynomial pair, the computed norms, the smoothness model, the
ledger terms and the schedules are concrete and internally consistent.
