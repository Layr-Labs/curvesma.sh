# BN254 fixed-argument pairing inversion via GT-DLP -- finite-regime ledger 2^123.304316

Exploratory, analytic cost argument for the protected BN254 pairing-inversion
target. Not a proof certificate, novelty claim, practical break or security-bit
estimate. The task is: given z = g^a (g = e(P0,Q0) public), return P in G1 and Q in
G2 with e(P,Q) = z. We use the reduction the rules sanction: recover the scalar a
by a Special extended Tower Number Field Sieve (SexTNFS) GT-discrete-logarithm, then
output P = P0 and Q = [a]Q0, since e(P0,[a]Q0) = e(P0,Q0)^a = g^a = z. Every derived
number is recomputed by `experiment.py`. Supporting files: DESCENT.md, KERNELS.md,
STORAGE.md.

The GT-DLP cost chain is identical to the sister `bn254-gt-dlp` analysis; the only
addition is the output work (one G2 scalar multiplication and one verification
pairing), which is ~2^9.3 reference units and negligible against the sieve. The
claim 2^123.304316 is 11.91 bits below the Pollard-rho baseline 135.216266.

## 1. Problem and strategy

Given only z = g^a (uniform a in [0,r)), return any P in G1, Q in G2 with
e(P,Q) = z; both points must lie in the stated prime-order subgroups, and either
may be the identity exactly when z = 1 (target-profile.json). Finding an arbitrary
pair does not require revealing source-group logarithms, but a GT-DLP solver gives
a concrete strategy (rules.md): recover a, output P0 and [a]Q0. We charge the
scalar multiplication and all other work; no equivalence of difficulty with GT-DLP
is asserted, only that this particular strategy costs what the ledger states.

p = 36x^4+36x^3+24x^2+6x+1 and r = 36x^4+36x^3+18x^2+6x+1 are 254-bit,
x = x0 = 0x44e992b44a6909f1 (63-bit), Q := p^12 has 3044 bits. z = 1 returns
(P0, O) or (O, O) directly and is excluded from the dominant terms.

## 2. Why NFS beats the generic baseline

The selected `rho-baseline-v2` returns P0 and [x]Q0 after a generic Pollard-rho
GT-discrete-logarithm (~2^127 group operations, scored 135.216266). GT is not
generic: it lives in F_{p^12} with composite embedding degree k = 12 and a prime
p = P(x0) that is the value of the fixed degree-4 BN polynomial. Generic exTNFS is
near 2^138 and does not beat the baseline; the special structure p = P(x0) is
load-bearing, and the ledger is built from an EXPLICIT special polynomial pair.

## 3. Explicit, verified polynomial selection (coprime (4,3) tower)

Identical to the gt-dlp construction; `experiment.py` exhibits and verifies the
triple (h, g, f):

- h(t) = t^4 + t + 3, **irreducible mod p** (base field F_{p^4}).
- g(X) = X^3 - x0: x0 is a verified **non-cube in F_{p^4}** (p == 1 mod 3), so g is
  **irreducible over F_{p^4}** and F_p[t,X]/(h, X^3 - x0) = **F_{p^12}**.
- f(X) = P(X^3) = 36X^12 + 36X^9 + 24X^6 + 6X^3 + 1, degree 12, infinity-norm 36.
  Since P(x0) == 0 mod p, X^3 == x0 mod (g,p) gives f == 0, so **g | f mod p**: the
  explicit common degree-12 structure. (f(x0) itself is non-zero for kappa > 1; the
  verified identity is f == 0 mod (g,p).)

A sieving element a(t,X) has degree < kappa = 3 over R = Z[t]/h with eta = 4 base
coefficients each bounded by the radius A (n = eta*tau = 12 integer coefficients).

## 4. Norms from the actual polynomials, and smoothness

The two norms are the exact integer resultants
N_f(a) = |Res_t(Res_X(a,f),h)|, N_g(a) = |Res_t(Res_X(a,g),h)|. `experiment.py`
computes them as big integers (Bareiss over Z[t], then an integer resultant with h)
on 48 real random sieve elements at radius A = 2^9: Nf ~ 522 bits, Ng ~ 612 bits.
The measured slopes in log2 A (48 = eta*deg_f, 12 = eta*deg_g) confirm the degrees;
the g-side carries size through ||g|| = x0 ~ p^{1/4}. At B = 2^62.4:
u_f = 8.37, u_g = 9.81, and P_smooth = rho(u_f)*rho(u_g) = 2^-60.16 from a
self-checked Dickman rho (two step sizes agree to 1.2%, ratios monotone).

## 5. The finite-regime ledger (reference units, C = 20000)

- **Relation collection.** |FB| = 2*B/ln B = 2^57.97; candidates = |FB|/P_smooth =
  2^118.12. Each candidate is charged a full dense-Fp12 reference unit (honest
  per-candidate ~88 instructions -> 2^110.30; the full-unit charge over-charges by
  +7.83 bits). Claimed: 2^118.12.
- **Linear algebra.** Block Wiedemann, D = 2^57.97, derived row weight ~18, 60
  Schirokauer columns, instruction-priced (weight*20 + schir*300 per entry,
  3*D^2 entries): 2^117.39. (KERNELS.md.)
- **Descent.** Capped Las-Vegas special-q ladder for both log_g and log_g z, per-
  level Dickman probabilities: 2^27, sub-dominant by >90 bits. (DESCENT.md.)
- **Output work (pairing-inversion specific).** Having recovered a, output P0 and
  [a]Q0. The G2 scalar multiplication is ~381 point operations over Fp2 (~10 Fp2
  muls each, ~4 Fp muls per Fp2 mul, ~300 instructions per 254-bit Fp mul) ~ 2^22.1
  instructions; one verification pairing (Miller loop + final exponentiation) ~ 2^23
  instructions. Together ~2^9.34 reference units -- negligible against the sieve.

Chain = 2^118.12 + 2^117.39 + 2^27 + 2^9.34 = **2^118.804**.

## 6. Claim and margin

log2(W) = 118.804 + 4.500 (disclosed margin) = **123.304316 = 123304316/1000000**,
11.91 bits below the baseline 135.216266. The 4.5-bit margin absorbs the asymptotic
o(1) and norm-model variance, filtering, Schirokauer-map evaluation, setup and
amplification; the ~7.8-bit per-candidate over-charge is a separate cushion already
in the chain. The honest instruction-level cost of the attack is ~2^117.4 and
concrete published BN254 estimates are near 2^100; the claim is above both.

## 7. Success probability (>= 0.39)

One SexTNFS run yields the factor-base logarithms; the only per-target random event
is the descent, whose capped loops (cap c/P_level, c = 40) give per-target success
~ 1 - 2e-14. The scalar multiplication and pairing check are deterministic. The
>= 0.39 bound holds for any per-log descent failure up to 0.30. Both output points
are in the prime-order subgroups by construction (P0 in G1, [a]Q0 in G2).

## 8. Heuristics, support and limits

Same disclosed heuristics as gt-dlp: H1 (smoothness on the computed norms),
H2 (relation yield), H3 (Schirokauer / class-number condition, O(1/r) ~ 2^-254),
H4 (the verified special construction), H5 (descent independence; Commeine-Semaev
2006, Barbulescu 2013). The work and memory caps hold even if H1-H5 are pessimistic;
only the success lower bound and the exact asymptotic constant depend on them. The
norm computation is exact at the production radius; the full attack is not executed.

## 9. Memory, preprocessing, advice

Peak storage ~2^67.1 bytes (relation buffer and sparse matrix at D = 2^57.97),
bound 2^70 bytes; memory is not scored. No embedded searched advice; no
amortization; both logarithms and the descent are charged in full.

## 10. Experiment

`experiment.py` (deterministic for the claim, stdlib-only, ~2 s, << 256 MiB)
constructs and verifies (h, g, f), computes the sieve norms by exact resultants on
real samples, the self-checked Dickman rho, the ledger terms, the output work, the
chain, the disclosed margin, the claimed rational 123304316/1000000, the storage
bound, and an independent cross-check with the caller's seed (within 1 bit of the
claim chain; the claim uses a fixed internal seed). Its JSON ends with
`"checks_ok": true`. It checks concrete internal consistency, not a production run.
