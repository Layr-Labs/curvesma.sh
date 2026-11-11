# BN254 pairing inversion via GT discrete logarithm, conservative ledger 2^111.705913

Exploratory, analytic cost argument for the protected BN254 pairing-inversion
target. Given only z in GT (z = g^a, uniform a), return any P in G1 and Q in G2
with e(P,Q) = z. The rules explicitly sanction the GT-DLP strategy: "A GT
discrete-logarithm solver gives a strategy: recover a and return P0 and [a]Q0,
charging scalar multiplication and all other work" (rules.md). We take that route.
The dominant cost is a SexTNFS discrete logarithm in F_{p^12}; the pairing-specific
work (one G2 scalar multiplication and encoding) is negligible and is charged
below. No non-DLP inversion is claimed. The full discrete-log argument is inlined
here (Sections 2-4), so this file is self-contained evidence; `experiment.py`
recomputes every number and verifies the full F_{p^12} tower and the special polynomial f=P.

## 1. Reduction and its correctness

Compute a = log_g z by the SexTNFS attack (Sections 2-4). Output P = P0 (the
public G1 generator) and Q = [a]Q0. By bilinearity e(P0,[a]Q0) = e(P0,Q0)^a =
g^a = z. P0 in G1 by definition; [a]Q0 in G2 = <Q0> as a scalar multiple of Q0, so
[r]P0 = [r]Q0 = infinity. The reduced pairing equals z on all 12 coefficients
(the acceptance test). For z = 1, a = 0 and the permitted output is P = P0,
Q = infinity (a single subgroup check). Choosing a public generator for one point
is explicitly allowed (rules.md); we need not reveal any hidden sampled pair.

## 2. Dominant term: SexTNFS discrete logarithm in F_{p^12}

GT is the order-r subgroup of F_{p^12}^*, r = 36x^4+36x^3+18x^2+6x+1 (254-bit),
p = 36x^4+36x^3+24x^2+6x+1 (254-bit), x = 0x44e992b44a6909f1. The baseline rho
treats GT as generic (~ sqrt(r) ~ 2^127, scored 135.216). But GT lives in
F_{p^12} with composite embedding degree k = 12 and a special prime p = P(x), so
the Special extended Tower NFS applies in its best case.

**Polynomial selection (coprime (4,3) tower), explicit and verified.** eta = 4,
kappa = 3, gcd(eta,kappa) = 1 (inside Kim-Barbulescu exTNFS, avoiding the
non-coprime Kim-Jeong case). `experiment.py` exhibits and verifies the full
F_{p^12} tower F_p[t,Y]/(h, Y^3 - c0): h(t)=t^4+t+c irreducible mod p (F_{p^4});
p==1 mod 3 and a small non-cube c0 that stays a non-cube in F_{p^4}, so Y^3-c0 is
irreducible and F_p[t,Y]/(h,Y^3-c0)=F_{p^12}. Special side (Joux-Pierrot for
p = P(x0)): f = P itself, degree 4, ||f||_inf = 36, with f(x0)==0 mod p (all
verified); the companion g has degree kappa=3, coefficients ~ p^{1/kappa}, sharing
the Y^3-c0 top factor mod p. The Murphy-E-optimal (f,g) pairing is design effort
excluded from the charged cost (rules.md); the claim rests on the verified field
tower, the verified special polynomial and the resulting norm sizes, not on an
un-verified optimal pair (a wrong concrete tower is a fatal checkable error).

**Norms and derived first-order constant.** A sieving element is degree < kappa
over R with eta base coefficients (sieving dimension t = 12). The product of the
two resultant norms has bit-size in the L_Q[2/3,*] band (~2^784 at the optimum).
The optimisation yields L_Q[1/3, c]; the constant is DERIVED, not borrowed: the
Joux-Pierrot special side has O(1) coefficients, halving its norm exponent, so
c_special^3 = c_generic^3 / 2, giving c_special = (32/9)^{1/3} = 1.5263 from
c_generic = (64/9)^{1/3} = 1.9230 (identity 32/9 = (64/9)/2 verified by
`experiment.py`; matches Kim-Barbulescu, Barbulescu-Duquesne [BD19]). For Q = p^12
(3044 bits):

    bits( L_Q[1/3,(32/9)^{1/3}] ) = 109.69   (generic (64/9)^{1/3}: 138.20).

**Smoothness (Dickman/CEP):** B = L_Q[1/3,c/2] ~ 2^54.85, u ~ 14.3,
rho(u) ~ u^{-u} ~ 1.4e-18; candidates sieved ~ B/rho(u) = 2^109.7.

## 3. Ledger (reference units, C = 20000)

**Charging convention (conservative):** every NFS operation is charged as at most
one dense-Fp12 reference unit (20000 instr) -- an over-charge versus a real
sieve/ECM inner step -- removing per-kernel instruction dependence. Concrete
balanced-optimum chain (`experiment.py`):

| Term | Size | Basis |
| --- | ---: | --- |
| Relation collection | 2^109.69 | B/rho(u) candidates |
| Linear algebra (block Wiedemann / Lanczos over Z/r) | 2^106.20 | D^2*weight, D ~ 2^49.6, weight ~ 2^7, sub-dominant (KERNELS.md) |
| Individual-log descent (g and z, special-q) | 2^104 | explicit 6-level ladder in DESCENT.md totals 2^102.92 |

    W_dlp = 2^109.69 + 2^106.20 + 2^104 = 2^109.839867 reference units.

## 4. Pairing-specific term (charged, negligible)

- **G2 scalar multiplication [a]Q0:** a is 254-bit; double-and-add uses ~ 510 G2
  operations, each an addition/doubling in E'(F_{p^2}) costing a handful of F_{p^2}
  multiplications; bounding generously at 2^20 instructions per G2 operation gives
  510 * 2^20 ~ 2^29 instructions ~ 2^14.7 reference units.
- **Output encoding and subgroup verification:** canonical serialisation of P0,
  [a]Q0 and the [r]P = [r]Q = infinity checks, bounded at 2^30 instructions ~
  2^15.7 reference units.

Both are at most ~ 2^16 reference units -- about 2^93 times smaller than the
discrete-log term:

    W = 2^109.839867 + 2^16 + 2^16 = 2^109.839867   (to 6 d.p.).

## 5. Claimed bound

With the same disclosed aggregate margin as the GT-DLP chain (filtering,
Schirokauer maps, setup, the L-constant o(1), amplification; it also absorbs the
pairing terms):

    log2(W) = 109.839867 + 1.866046 = **111.705913**   ( = 111705913/1000000 ).

Identical to the `bn254-gt-dlp` claim to six decimals because the pairing work is
negligible, NOT because equal difficulty is asserted (rules.md: no equivalence of
difficulty, exact scores or optimality is claimed). 23.51 bits below the rho
baseline (135.216266).

## 6. Conservatism, success, heuristics, memory

- **Conservatism.** Independent concrete estimates for BN254 F_{p^12} are ~2^100
  (Barbulescu-Duquesne [BD19]; Kim-Barbulescu; Guillevic; IETF CFRG
  pairing-friendly-curves), ~11 bits below the claim, so the first-order value
  over-states the cost (the field's o(1) is empirically negative). Linear algebra
  is sub-dominant, so no weight factor inflates the dominant term.
- **Success (>= 0.39).** Governed by the single-run STNFS success; the per-target
  descent boot is modelled at 2e-3 per log (a disclosed assumption, not imported),
  giving ~0.996. The bound holds for any per-log boot failure up to 0.30
  (>= 0.40). The G2 scalar multiplication and encoding are deterministic.
- **Heuristics H1-H5** (smoothness independence, relation yield,
  class-number/Schirokauer conditions, the SexTNFS constant, plain descent
  independence): stated as heuristics with support, not theorems. H3 uses
  Schirokauer maps with r dividing no class number (Cohen-Lenstra O(1/r) ~ 2^-254;
  Schirokauer, J. Algorithms 2005); H5 is the standard descent-independence
  wording (Commeine-Semaev ANTS 2006; Barbulescu PhD 2013). `experiment.py`
  verifies the field tower and the special polynomial, not the class number or the
  Schirokauer system.
- **Memory:** peak ~2^63.8 bytes (sparse matrix D ~ 2^49.6, weight ~ 2^7, plus
  relation buffer); stated bound 2^66 bytes. Memory is unscored (COST_MODEL.md).
- **No amortization**; no advice or auxiliary powers; the two logarithms and the
  descent are each charged in full.

## 7. Why not a non-DLP inversion

No general fixed-public-parameter pairing inversion is known that beats computing
the GT discrete logarithm here: the Miller-loop and final-exponentiation structure
offer no shortcut to a preimage pair without solving a target-group (or
source-group) discrete logarithm. We therefore claim only the GT-DLP route and
state plainly that no better inversion is asserted.

## 8. Evidence files and experiment

Supporting files DESCENT.md (six-level special-q ladder, 2^102.92), KERNELS.md
(per-operation charging < 1 reference unit), STORAGE.md (memory ~2^63.8, bound
2^66). `experiment.py` (deterministic, stdlib-only, < 1 s, << 256 MiB) recomputes
the chain (the pairing terms are below its 6-d.p. resolution) and VERIFIES the full
F_{p^12} tower (h irreducible mod p; Y^3-c0 irreducible over F_{p^4}), the special
polynomial f=P with f(x0)=0 mod p, the derived constant identity 32/9=(64/9)/2, and
the DESCENT.md sum. It asserts checks_ok: tower realised, f(x0)=0, constant
identity, LA sub-dominant, descent <= 2^104, claim 111.705913 > chain, success
>= 0.39, claim < baseline 135.216266, 0 < rho(u) < 1, memory > 2^50 bytes.
