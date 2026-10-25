# Additive Cheon v30: prime class lists — the exact energy decomposition, the residual's certified strata, one named joint-defect bound

**Lineage.** Provenance: the construction follows the published `additive-cheon-v1`
package (submission cbadde77, in review on this track, declared as a dependency).
It is restated in full below, and every identity, instruction itemisation, bound,
price and toy number in THIS proof is this package's own, re-derived here and
re-verified by the runs of THIS package's `evidence.py`; nothing is inherited by
citation. The accounting differences against the provenance package, carried
since v2:

1. canonical y: itemised below at ~25 instructions, charged **30**.
2. walk bookkeeping: itemised below at ~75 instructions, charged **100**.
3. hash insert/lookup: a complete per-transfer accounting, itemised below:
   endpoint work 61 plus six radix-1024 relay stages of 29 instructions over the
   ≤ 2^60 nodes, 235 itemised, charged **240**; evidenced by a 16,384-key probe
   and shard check on structured production-style walk keys.

**The v23 change (carried): the b's are the PRIMES of the χ-class lists.** With
every walked element a product a·b of two primes from disjoint ranges,
integer-level coincidences off the diagonal are impossible by unique
factorisation. Parameters: X1 = 2^82 (N_A ≈ 2^73.32 primes), X2 = 2^70, GMAX = 2^34.

**The v28 change (carried): the calibrated null.** The birthday law was tested
against its own null — multiple independent random-set realizations counted by
the same routine at the same parameters — at two scales spanning a 16-fold
change in d and both measured class ratios, the construction lying INSIDE both
null distributions while the structured baselines lie far outside.

**The v29 change (carried): the success argument restructured into proved
families.** The distinct-index collision energy decomposes exactly (§3.1);
the mixed-root integer identities are proved EMPTY by a cutoff-free
certificate over the six fixed constants (§3.2); the repeated-coordinate
family is bounded at 2.671×10⁻¹⁹ of the birthday baseline (§3.3); the
shared-prime wrapped collisions at ρ = ±1 are excluded by the retained
2·X1²·X2 < r, whose strength is demonstrably necessary (the r = 481,153
counterexample, §3.4); and the same-root swap is reduced to the emitted-mass
class second moments under P1 (§3.5).

**The v30 change: the residual's strata are certified and the premise
shrinks to one joint-defect bound.** The algorithm, walks, ledger,
parameters and score are unchanged. What changes is the second half of §3:
the residual is decomposed by common-factor structure, and every structured
stratum is now PROVED at production scale under P1 with no H2 clause — the
exactly-one-shared-large stratum at 0.002332911916573675 of the birthday
baseline B₀ (an affine-coset lattice theorem with the mod-30 wheel's exact
16/225 area factor, weighted by the actual P1 class masses, at
Dusart-certified block endpoints a_B ≤ 0.57·X1, a_G ≤ 0.64·X1), the
shared-small Selberg main at 0.038976659226098 and its list-prime exception
at 0.002641842 (finite optimal Selberg weights at R = 2^35 with an
integer-certified hyperbola denominator G(2^35) > 322.7629594, no analytic
input), the two-large-only correction at 1.63×10⁻¹⁷, the repeated-coordinate
charge at 2.671×10⁻¹⁹, and the coprime main term certified to 3.5×10⁻¹³.
A coverage-preserving pruning discharges the swap-energy charge entirely
(proved density loss ≤ 4.79×10⁻¹⁴, no H2). The five certified pieces sum to
α = 0.043951414 of B₀, and the single remaining unproved estimate is the
JOINT DEFECT t = (Δ₀ + R_SS)/B₀ — the coprime modular discrepancy plus the
signed Selberg remainder — on which the claim is conditional.

**Claim.** For every n with 2^20 ≤ n < r, the deterministic algorithm below recovers s with heuristic
success probability at least 0.390 (≥ 0.39 required; route A of §3.7 — the DETERMINISTIC
pruned-coverage implication: κ − (1+α+t)·K²/2 ≥ 0.39000000388245215 at the safe joint-defect bound
t ≤ 0.0138998, with α = 0.043951414 certified at production under P1; route B gives ≥ 0.4120276 at
the sharper target t = 0.0121794 additionally under P3/P4 — t itself is the single named unproved
estimate). Its total work is
W(n) = W(2^20) ≤ **2^116.958777** target-group-operation equivalents (exact score 116958777/1000000,
rounded up; unchanged from v28 — the cost ledger is identical).
- Preprocessing, meaning the challenge-independent work, is ≤ 2^115.968 (exactly the ledger's own
  preprocessing output) and is included in W.
- Peak memory is ≤ 2^122.9 bytes. There is no advice.
- The algorithm reads only A_1 = T and A_d, with d = 2^20.

The idea is Cheon's (p−1)-attack, run as a baby-step/giant-step (BSGS) inside the order-m subgroup H of
d-th powers, with one change. The baby and giant sets are **small products of two primes that lie in H**,
not geometric progressions ζ^i. Products a·b listed in increasing order can be walked with **one batched
affine addition per element**. A geometric step needs a scalar multiplication, which is 2 additions with
the Kozaki–Kutsuma–Matsuo (KKM) c = 3 tables.

## 1. Facts about the fixed target

These are checked by `evidence.py` and by direct computation.
- r ≈ 2^254.857 is prime and v₂(r−1) = 32, so d = 2^20 divides r−1.
- H = {x^d : x ∈ F_r^*} = ker(x ↦ x^m) has order m = (r−1)/d ≈ 2^234.857, and 6 | m.
- **The hidden value is uniform on H.** t := s^d is uniform on H, because s is uniform on F_r^* and
  x ↦ x^d is d-to-1 onto H. A_d = [t]P0 is public.
- **The symmetry group sits inside H, with production-grade structure.** Let
  L = x² − 1 = 228988810152649578064853576960394133503 ≈ 2^127.429, where x is
  the BLS12 parameter 0xd201000000010000. Then **L² + L + 1 = r exactly**, as an
  integer identity (r = x⁴ − x² + 1), so L is a primitive cube root of unity
  mod r, [L](x, y) = (βx, y) is the GLV endomorphism, and U = {±1, ±L, ±L²} ⊂ H.
  The six fixed representatives of U in (0, r) are exactly the polynomials

  | representative | polynomial | order in F_r^* |
  |---|---|---:|
  | 1 | 1 | 1 |
  | L | x² − 1 | 3 |
  | L + 1 | x² | 6 |
  | L² | (x² − 1)² | 3 |
  | L² + 1 | (x² − 1)² + 1 | 6 |
  | L(L + 1) = r − 1 | x²(x² − 1) | 2 |

  (three of 2^127.43 scale, three of r scale; the representative of −L̄ is L + 1
  and of −L is L² + 1, since −L² ≡ L + 1 and −L ≡ L² + 1). Also
  r − 1 = L·(L + 1) and L + 1 = x², facts the certificate of §3.2 consumes.
- **One key identifies a whole orbit.** All six points [u]Q, for u ∈ U, have y-coordinate ±y(Q),
  because λ fixes y. So canonical(Q) := min(y, p − y) identifies the U-orbit of Q.
- **Membership in H is a character test.** χ(x) := log_ζ(x^m) ∈ Z/d is a character with kernel H, for
  a fixed generator ζ of μ_d. So a·b ∈ H ⇔ χ(a) + χ(b) ≡ 0 (mod d).
- **2 ∉ H.** Indeed 2^m ≠ 1 (mod r); 2^m has order 2^16.

## 2. Algorithm

Parameters: K = 0.55, X1 = 2^82, X2 = 2^70 < X1/2, GMAX = 2^34 and N_B = N_G = N = ⌈√(K·m/6)⌉ ≈ 2^115.705.

**The deterministic set rule (fixes B and G given the input).** The primes a
are taken from (X1/2, X1] in increasing order, each used exactly once: the
first ⌈N_A/2⌉ form the baby side and the next ⌈N_A/2⌉ the giant side, where
N_A = (N_B+N_G)·2d·49/X2 + 2 ≈ 2^73.32 (the 49 = ln(2^70) + slack buys the
H2 minimum supply per class, below). The primes of (X1/2, X1] number at
least 2^74.716 by the Rosser–Schoenfeld interval bound — a proven
**2.63×** slack over N_A (`log2_primes_interval_lower_bound`,
`log2_prime_availability_ratio_lower` in the output). For each prime a, the b's
are **the primes of L_{−χ(a)} in increasing order** — exactly the order the
χ-sieve streams them — and the walk stops for that prime when its side
reaches N elements. Nothing is sampled or shuffled; B and G are fixed by the
input. The evidence's toy instantiation of this rule is the same formula
with the toy's ln X2 in place of 49: **ka = 2n·d·LN/X2 + 1 with
LN = ⌊ln X2⌋ + 1**, implemented with no fallback, and the attack uses the
same form at its own parameters (ka = 41 there). One attack configuration
runs in fixed order (16 trials) alongside the shuffled runs, which
robustness-check the rule against unfavourable prime blocks.

1. **Precompute** the characters, sorted into prime class lists.
   - Compute χ(q) for the primes q ≤ X2 and for the primes used in (X1/2, X1]. Each costs one 235-bit
     exponentiation in F_r and a Pohlig–Hellman logarithm in μ_{2^20}.
   - Sieve χ(b) for all b ≤ X2 by multiplicativity (the sieve necessarily processes every integer;
     **only the primes are stored to the class lists**).
   - Stream each prime b ≤ X2 into its class list L_c = {b ≤ X2 prime : χ(b) = c}, sorted in
     increasing order. Each list holds ≈ π(X2)/d ≈ 2^44.7 primes (H2).
2. **Baby walk.** For fresh primes a ∈ (X1/2, X1]:
   - set Q_a = [a]P0 and build T[k] = [k]Q_a for k ≤ GMAX (T[2] by doubling, then T[k] = T[k−1] + Q_a);
   - walk b over the primes of L_{−χ(a)} in increasing order, with [ab]P0 = [ab_prev]P0 + T[b − b_prev];
   - for a gap greater than GMAX, first add T[GMAX] repeatedly;
   - insert key canonical([ab]P0) with locator (a, b) into a hash table.

   Stop after N elements. Every element a·b lies in H. They are distinct integers below r: a is the
   unique prime factor of a·b above X2, and with b prime the integer a·b determines both factors.
3. **Giant walk.** Do the same with disjoint primes a′ and base A_d, which gives the points
   [a′b′·t]P0. Look up the canonical key of each point. Stop after N elements, or earlier on success.
4. **Collision.** A key match means [g·t]P0 = [u·ab]P0 for some u ∈ U, where g = a′b′. Test the six
   candidates t = u·ab/g by checking [t]P0 = A_d.
5. **Recover s.** Take one d-th root s0 of t: twenty square roots, keeping at each step the root that is
   still a 2^j-th power. Then s = s0·ω^i for some i < 2^20, where ω generates μ_d. Test [s0·ω^i]P0 = T
   for each i and output s.

**Batching and representation.** Each prime's walk, and its table build, is sequential. The 2^16
primes in flight are therefore advanced in lockstep. Every round performs 2^16 independent incomplete
affine additions, one per prime, sharing a single Montgomery batch inversion; this is the reference
operation's batched setting. The inversion costs ≈ 100 units per round, about 0.0015 units (3
instructions) per element, which fits inside the walk-bookkeeping budget. All stored points, meaning
tables, walk states and Q_a, are affine: 96 bytes each.
- Q_a = [a]P0 (or [a]A_d) is computed in Jacobian coordinates and normalised by one inversion, charged
  in the per-prime scalar budget.
- T[2] is one affine doubling, also charged there.

**Incomplete additions.** A walk step adds [b_prev]Q and [gap]Q, with Q ≠ O and b = b_prev + gap < r. The
x-coordinates can match only if b_prev ≡ ±gap (mod r).
- b_prev ≡ −gap would force b ≡ 0, which is impossible.
- b_prev = gap means b = 2·b_prev in the same class list: 2·b_prev is an even integer above 2, hence
  composite, and never in a prime list. (The χ(2) = 0 route is also closed: 2 ∉ H.)
- *Overflow chunks.* [b_prev + i·GMAX]Q + [GMAX]Q collides only when b_prev = GMAX and i = 0. That is at
  most once per prime; it is done by a doubling, charged to the per-prime scalar budget. The remainder
  addition is safe because b_prev + j·GMAX > GMAX ≥ rem for j ≥ 1.
- *Table build.* T[k] = [k−1]Q + Q is safe for k ≥ 3.
- Giant points are multiples of Q = [a′t]P0 ≠ O, so the same argument covers them.

## 3. Success probability

Steps 4–5 are deterministic. The algorithm succeeds exactly when t ∈ S := U·B·G^{-1}. Because t is
uniform on H, Pr[success] = |S|/m. Write N = 6·N_B·N_G for the total number of walked triples
(u, x, y), u ∈ U, x ∈ B, y ∈ G.

### 3.1 The exact distinct-index energy decomposition

Let E be the number of ordered pairs of **distinct** triples with equal coverage
values u·x·y^{-1} ≡ u′·x′·y′^{-1} (mod r). The diagonal (identical triples)
belongs to E + N, not to E. Split each congruence by k, the integer quotient in
u·x·y′ − u′·x′·y = k·r:

- **k = 0, same root (u = u′).** Unique factorisation over the two block primes
  (each > X1/2; the baby and giant blocks are disjoint) forces the block primes
  to match on both sides, and then the four class primes force either the
  identical pair or the swap b₁ = b′₂, b₂ = b′₁: a same-class row pair
  (a, a′) whose consumed prefixes share s = min(t_a, t_a′) elements. Every
  shared element maps both triples to the single value u·a/a′, so the pair
  contributes exactly s(s−1) ordered distinct-index pairs per root — the
  **shared-class-prefix family**, counted exactly: 6·Σ s(s−1). This includes
  the final partial prefixes exactly (s is the consumed length, not the list
  length).
- **k = 0, mixed root (u ≠ u′).** Empty — §3.2.
- **k ≠ 0.** The residual E_{k≠0}, bounded in §3.3–§3.6.

$$E \;=\; 6\sum_{\text{same-class row pairs}} s(s-1) \;+\; E_{k\ne0}.$$

Both sides of this identity are computed by `evidence.py` at every measured
point (`energy`: the swap term, the mixed term, and E_{k≠0}, with the family
split of §3.6's measurements).

### 3.2 Mixed-root integer identities are empty — cutoff-free

A mixed integer identity u·P = u′·Q with P = x·y′, Q = x′·y (four-prime
products) gives, after reducing u/u′ = v/w to lowest terms, v | Q and w | P,
hence the necessary condition

$$\Omega(v) \;=\; \Omega(w) \;\le\; 4$$

(Ω with multiplicity; P and Q carry exactly four prime factors). All **15
unordered mixed pairs** of the six constants of §1 fail it (`mixed_root_k0_certificate`):

- **eight** have a reduced coefficient equal to 1 against a coefficient > 1
  (Ω 0 against ≥ 1): the pairs {1, ·} together with {L, L²}, {L, L+1},
  {L+1, L(L+1)}, {L², L(L+1)}, whose shared parts strip an entire constant;
- **five** have a reduced coefficient divisible by 2^31: the shared parts strip
  only {1, L, L+1, 2}, and L + 1 = x² carries v₂ = 32 (Ω ≥ 31);
- the **remaining two** contain L² + 1, and 2·5·37·1637 | L² + 1 with
  quotient > 1, so Ω(L² + 1) ≥ 5 > 4.

No size hypothesis is used anywhere in this certificate — it is a statement
about the six fixed constants of the target alone. The corresponding count is
measured **0** at every toy point of the evidence.

### 3.3 The repeated-coordinate family is negligible at production

At small toy d, elevated family counts can come from **internal torsion
pairs** — elements x ≡ ρx′ (mod r) inside one side, which turn
repeated-coordinate collisions (x = x′ or y = y′) into whole families. At
production this mechanism is provably absent:

- For **real ρ** (±1): an internal pair x ≡ −x′ would need x + x′ = r with
  x, x′ ≤ Z := X1·X2 = 2^152 < r/2 — impossible; x = x′ is the same element.
- For **nonreal ρ**: an internal torsion pair must be **coprime** — a shared
  prime p | gcd(x, x′) would force r | p² ± pq + q² with
  0 < p² ± pq + q² ≤ 3·X1² = 2^165.58 < r (89.3 bits of room). Coprime
  torsion pairs are primitive lattice points: at most ⌊2Z²/r⌋ + 1 per side, and
  their entire repeated-coordinate contribution satisfies

$$\frac{D_{\rm repeat}}{N^2/m} \;\le\; \frac{8\,(2Z^2/r + 1)}{K\,N} \;=\; \mathbf{2.671\times10^{-19}}$$

  (`d_repeat`) — twenty orders of magnitude inside every budget below.

### 3.4 Shared-prime wrapped collisions at ρ = ±1: excluded, and the exclusion is tight

A k ≠ 0 collision with a prime shared between the two elements of one side
reduces, after cancelling the shared prime, to an integer relation among
products of at most three allowed primes:

- shared **block** prime (one row emitting two elements): the three-prime
  products are ≤ X1·X2² = 2^222 < r;
- shared **class** prime (a same-class row pair sharing consumed prefix): the
  products are ≤ X1²·X2 = 2^234.00 < r, and for ρ = −1 their sum is
  ≤ 2·X1²·X2 = 2^235.00 < r.

In every case the reduced relation is an integer of absolute value < r that
must be a nonzero multiple of r — impossible — or zero, which is an integer
identity already classified at k = 0 (§3.1–§3.2). So **every shared-prime
collision at ρ = ±1 is either impossible or counted at k = 0**, under the
retained inequality

$$2\,X_1^2 X_2 \;=\; 2^{235.00} \;<\; r \;=\; 2^{254.857} \qquad \textbf{(19.857 bits of slack).}$$

The stronger of the two monomials is **necessary**: at r = 481,153 with
d = 16, X1 = 400, X2 = 197 — where 3·X1² = 480,000 < r but
X1²·X2 = 4.30×10⁹ > r — the class-matched semiprimes

- P = (227·167)(347·197), Q = (241·73)(347·151),

all four in H (verified as (ab)^m ≡ 1) with disjoint baby/giant tiers, satisfy
P − Q = 3470·r with gcd(P, Q) = 347: a genuine shared-prime wrapped collision
(`counterexample_r481153`). We retain 2·X1²·X2 < r precisely because the weaker
condition does not exclude this family.

### 3.5 The same-root swap: emitted-mass class moments (premise P1)

Let t_a be the number of elements row a actually emits (the final partial
prefix included) and W_s(c) = Σ_{a on side s, χ(a)=c} t_a the emitted mass per
block class (Σ_c W_s(c) = N/6 = n per side). Since
min(t, t′)[min(t, t′) − 1] ≤ t·t′,

$$E_{\rm swap} \;\le\; 6\sum_c W_B(c)\,W_G(c) \;\le\; \frac{6n^2\bigl(1+\sqrt{\delta_B^{w}\delta_G^{w}}\bigr)}{d},$$

the second step by orthogonality and Cauchy–Schwarz on the class characters,
with the **emitted-mass class second moments**

$$\delta_s^{w} \;=\; \frac1{n^2}\sum_{j=1}^{d-1}\Bigl|\sum_a t_a\,\chi(a)^j\Bigr|^2 .$$

**Premise P1:** δ_B^w, δ_G^w ≤ 500. Then

$$\frac{E_{\rm swap}}{N^2/m} \;\le\; \frac{501}{K\,d} \;=\; \mathbf{8.69\times10^{-4}}$$

— 0.0869% of the birthday baseline (`swap_bounds.wform_fraction_of_baseline`).
P1 is a precise, checkable second-moment statement — not pointwise
equidistribution — and it is measured with 59–737× headroom at all three
evidence scales:

| point | d | δ_B^w | δ_G^w |
|---|---:|---:|---:|
| coverage, X2/d = 64 | 2^4 | 0.678 | 0.679 |
| structured-λ (§3.6) | 2^4 | 1.399 | 1.888 |
| coverage, X2/d = 16 | 2^8 | 7.225 | 8.490 |

Two corroborations, both verified in the same run: the **full-interval form**
— since the two row sets are disjoint subsets of I₀ = (X1/2, X1],
Σ_c n_B(c)n_G(c) ≤ (T² + Δ₀)/(4d) for the ACTUAL stopped subsets, so the
premise Δ₀ ≤ T² over I₀ gives E_swap/(N²/m) ≤ 1.452×10⁻³
(measured Δ₀/T² ∈ {0.066, 0.094, 0.505}); and the v28 **L2-form**
6·(⌈N_A/2⌉²/d)·(L2 − 1) = 2^172.6 lost values, unchanged in the ledger. The
list-size concentration bound L2 is no longer load-bearing for the swap; the
block rule's supply side of H2 is unchanged.

### 3.6 The residual's certified strata and the single remaining premise

Both applications consume moments of ONE object, the **labelled prime-product
Fourier transform** 𝒯_s(j,η) = Σ_a Σ_{b∈P_a} χ(a)^j η(ab) (slices
𝒯_s(j,1) = the emitted-mass class moments of §3.5, 𝒯_s(0,η) = Σ_{x∈B_s} η(x)),
and Fourier analysis on H gives the exact **U-filter identity**
(Û(η) = 6·1[η|_U = 1] for characters η of H):

$$E + N \;=\; \frac{N^2}{m} \;+\; \frac{36}{m}\sum_{\substack{\eta \ne 1 \text{ on } H\\ \eta|_U = 1}}
\bigl|\widehat B(\eta)\bigr|^2\,\bigl|\widehat G(\eta)\bigr|^2 ,$$

with ≈ m/6 − 1 ≈ 2^232.3 contributing characters; the identity is verified
EXACTLY at both energy points of the run (the U-coset tuple count and
(m(E+N) − N²)/36 agree as integers).

**The coverage-preserving pruning (no H2, no swap-energy premise).** Prune
each deterministic same-root shared-prefix cluster to one tuple, consistently
across all six roots: every deleted index has a retained equal value, so the
support is unchanged, and the clusters are disjoint (the two large factors
identify the row pair). The deletion count satisfies
D = 6·Σ_{(a,c)}(s(a,c)−1)⁺ ≤ 6·Σ_{(a,c)}s(a,c) = 6·Σ_p e_B(p)e_G(p) ≤ 6Mn
(using e_s(p) ≤ M and Σ_p e_B(p) = n), so **D/m ≤ 6Mn/m = 4.79×10⁻¹⁴** — a
proved density loss, no H2 clause — and the pruned energy E₀ ≤ E − E_swap
carries no swap charge at all (the 8.69×10⁻⁴ of §3.5 is discharged; P1
remains the input to the strata below).

**The covering of all shared-prime wrapped collisions (proved).** After
pruning, every collision falls into exactly one of:
(i) mixed-root integer k = 0 — absent (§3.2's certificate);
(ii) shared-prime wrapped at ρ = ±1 — absent (2·X1²·X2 < r, §3.4);
(iii) repeated coordinates — ≤ 48n·(⌊2Z²/r⌋+1)/B₀ < 2.671×10⁻¹⁹ (§3.3);
(iv) a nonrepeated collision with at least one common SMALL prime — covered
by the broadened Selberg bound E ≤ E_exc + V/G(R) + R_SS: the main term with
the finite optimal Selberg weights (λ_d at R = 2^35, λ₁ = 1, quadratic main
exactly 1/G(R)) is CERTIFIED ≤ **0.038976659226098·B₀** under P1, the
list-prime exception (S₃₅ = 1 + Σ_{j=2}^{35}⌊(2^j−1)/(j−1)⌋ = 2,086,667,226
by the dyadic binomial argument) ≤ **0.002641842·B₀**, and R_SS is the
SIGNED weighted remainder — retained exactly, not assumed away; the
denominator G(2^35) > 322.7629594 is certified by an integer hyperbola
computation with no analytic input (selberg-derivation.txt);
(v) a nonrepeated collision whose common factors are only LARGE: exactly
one common large — the affine-coset lattice theorem, CERTIFIED ≤
**0.002332911916573675·B₀** under P1 (the mod-30 wheel gives the exact
16/225 area factor; the weighted short-vector step
F_s(T) ≤ Σ_c W_s(c)√k_c(T) ≤ ‖W_s‖₂·√(Σ_c k_c(T)) preserves the actual P1
class masses; explicit root-lattice cells carry ≤ 9 coefficient antecedents
per primitive point; the block endpoints a_B ≤ 0.57·X1, a_G ≤ 0.64·X1 are
certified from the stated Dusart bounds — shared-large-derivation.txt); two
common larges — ≤ 96·M²·(⌊2X₂⁴/r⌋+1)/B₀ = 1.633×10⁻¹⁷;
(vi) coprime — E_coprime = 36·M₀/m + Δ₀ with the coprime tuple density
certified to **1 − M₀/n⁴ < 3.492×10⁻¹³** (the main term is B₀ to thirteen
decimals) and Δ₀ the modular discrepancy.

The self-moment route to these bounds is blocked, exactly as measured: the
diagonal contribution alone forces S_B, S_G ≥ h(2n²−n) − n⁴ on the sector
fourth moments (verified per run), so separate self-moment bounds plus
Cauchy–Schwarz cannot reach below 1/K − 1 ≈ 0.818 (measured Cauchy ratios
1.985 and 2.064) — direct control of the mixed correlation, stratum by
stratum, is what the certificates above do.

**The conditional reduction and the single remaining premise.** The five
certified pieces sum to the ledger bound α = **0.043951414** (a directed
9-place ceiling over the exact rationals, recomputed fresh by the run from
the three in-folder source modules; the only external analytic input is
Dusart arXiv:1002.0442v1 Theorem 6.9 eq. (6.6), stated with its validity
domains in prime-count-dependency.txt and transcription-checked against
independent sieves). With t bounded:

> **Premise (the joint defect).** t := (Δ₀ + R_SS)/B₀ ≤ 0.0138998 —
> the coprime modular discrepancy plus the signed Selberg remainder.

E₀ ≤ B₀·(1 + α + t). **t is named, not proved.** Naming it does not
discharge a demand for its proof; the exact maximal t in the certified
route-A formula is 0.013899825669104974.

**The structured-λ experiment** (`structured_lambda_point`): the cyclotomic
toy r = 48⁴ − 48² + 1 = 5,306,113, whose small cube root L_c = 48² − 1 = 2303
satisfies L_c² + L_c + 1 = r exactly — the production constant structure of
§1 in miniature, three small and three large representatives, asserted in
code. At d = 16, X2 = 1024, n = 175, ka = 39 (rows 15/19 consumed):
|S|/m = 0.426 against the law's 0.4254; E = 108,984, swap = 9,708, mixed = 0,
E_{k≠0} = 99,276 = **−2.4907%** of baseline; per-order family rates
(order 1 / 2 / 3 each / 6 each): 0.926415 / 0.996426 / 0.990062 / 0.973796 —
**every family below null**. The production-like root structure creates no
family-specific excess. At the d = 2^8 coverage point the rates are
1.026129 / 0.960423 / 0.974403 / 1.1156 — the only above-null families sit at
order 6 there, an order that is below null at the structured point: no
consistent order pattern exists at any measured scale.

### 3.7 The claimed success

Both routes are DETERMINISTIC coverage implications — no birthday law, no
coupling assumption, no Poisson approximation:

- **Route A — the claim; second Bonferroni on the pruned multiset.** For
  every finite occupancy multiset, coverage ≥ F₁ − F₂/2. With
  κ = K − 6Mn/m ≥ K − 4.79×10⁻¹⁴ and F₂ = E₀/m ≤ (1+α+t)·K²/2:
  |S|/m ≥ κ − (1+α+t)·K²/2. At the safe gap t ≤ 0.0138998 this is
  ≥ **0.39000000388245215**, floor-stated **0.390** ≥ 0.39 — the claim,
  conditional on P1 and the named t. At the sharper target t = 0.0121794:
  ≥ 0.390260.
- **Route B — fourth Bonferroni; additionally P3/P4.** Coverage ≥
  F₁ − F₂/2 + F₃/6 − F₄/24 on the pruned multiplicities; with P3
  (F₃ ≥ 0.95·κ³) and P4 (F₄ ≤ 1.20·κ⁴) — both unproved, measured 0.959/1.024
  and 0.900/0.976 at the two energy scales, the inequality verified per run —
  and u + g(u) increasing on the certified range (g(u) = 0.95u³/6 −
  1.20u⁴/24): |S|/m ≥ κ + g(κ) − (1+α+t)·K²/2. t ≤ 0.157816 ⟹ ≥
  0.39000007446578 (t = 0.158 fails: 0.38997); t ≤ 0.10492393035226562 ⟹
  ≥ 0.398; at t = 0.0121794 ⟹ ≥ 0.4120276.
- **The original-total caveat (wording discipline).** The pruned
  implications must NOT be described as proving the original total-energy
  P2(0.057): unpruned, ε_total ≤ α + t + β = 0.058720214 at t = 0.0138998,
  where β = 8.69×10⁻⁴ is P1's swap charge; only at the sharper t = 0.0121794
  does the unpruned total (0.056999814) stay under 0.057.

The claim is route A's **0.390** — the smallest assumption set that clears
the 0.39 requirement: P1 plus one named bound on the joint defect. The
law-based figures of the earlier lineage (0.42305, 0.415) are retained only
as the §3.8 empirical context.

### 3.8 The law against its calibrated null (carried from v28, unchanged)

| scale | law | prime construction | null range (draws) | null mean ± sd | z |
|---|---:|---:|---:|---:|---:|
| d = 2^4, X2/d = 64 (r ≈ 2^24.24, n = 338, ka = 74) | 0.4245 | **0.4219** | 0.4216–0.4252 (8) | 0.4235 ± 0.0012 | −1.4 |
| d = 2^8, X2/d = 16 (r ≈ 2^25.33, n = 124, ka = 140) | 0.4278 | **0.4242** | 0.4238–0.4303 (6) | 0.4282 ± 0.0023 | −1.7 |

At both scales the construction lies INSIDE the null distribution. The
structured baselines at d = 2^4 lie far outside every null range: the matched
integer-list control measures 0.3296 (+22.4%) and the smallest integers of H
measure 0.3952 (−6.9%). The counted shared-class-prefix lost values are 3,858
at d = 2^4 and 66 at d = 2^8, tracking the analytic L2-form.

**The end-to-end attack confirms the law.** 40 shuffled trials on the toy
curve: success 0.55 against the law's 0.4267 (+1.48 σ); 16 fixed-order trials
0.4375 (+0.11 σ); every recovered scalar correct in both runs (the toy attack
runs wrapped, unlike production — figures stated as measured).

## 4. Cost (unit: one batched incomplete affine G1 addition = 2,000 instructions)

Per-element budgets, in instructions on top of the reference addition.

| Item | Budget | Itemised count (this package's own ledger) |
|---|---:|---|
| canonical y | 30 | p − y (≈ 7), 6-limb compare (≈ 12), select (6): ≈ 25 |
| hash insert (baby) / lookup (giant) | 240 | endpoint 61 + 6 relay stages × 29 (18 word moves + 4 destination selection + 4 queue updates + 3 control); evidenced by the hash_probe toy check |
| walk bookkeeping | 100 | read b and the gap (4), overflow test (2), address T[gap] (3), load T[gap] (12), load and store points (24), batch-inversion buffers beyond the 6M (24), loop (5): ≈ 75 |

**Hash breakdown.** The 2^116.7-slot table, at load ≤ 1/2 with 72-byte records, cannot sit in one
address space. It is sharded over nodes by the high bits of the key. Each insert or lookup:
- computes the shard and slot from key bits: 4;
- appends the 9-word record to that shard's outbound buffer: 14;
- is forwarded through h hops, each loading and storing every word: 18 per hop;
- is read at the receiver: 9;
- is probed by linear probing at load ≤ 1/2, ≤ 2.5 expected probes of 6 instructions: 15;
- stores the record, for inserts: 9;
- costs loop and control: 10.

Apart from the hops, the work is 61 instructions. The budget of 300 therefore covers h ≤ 13 hops. A
radix-1024 switch fabric over the same ≤ 2^60 nodes needs ceil(60/10) = 6 hops, at the identical
per-stage relay cost (18 instructions load and store each of the 9 record words). Lookups follow the same
route, and a match reply is needed only on success. These are expected costs summed over 2^117
operations, where linear probing and balanced sharding concentrate sharply. The χ-sieve is segmented
over ranges of b, so each segment is local, and appends to the prime class lists are streamed. Walks
read one list sequentially.

This gives c_B = c_G = 1 + 370/2000 = 1.185. The full ledger, recomputed exactly by `evidence.py`
(every charge above its itemised count):

| Term | log₂ (units) |
|---|---:|
| baby walk N·c_B | 115.95 |
| giant walk N·c_G | 115.95 |
| gap-overflow additions (deterministic: Σ gaps ≤ X2 per prime, ≤ N_A·X2/GMAX) | 109.32 |
| per-prime tables T[k] (GMAX additions × N_A primes) | 107.32 |
| per-prime scalar multiplications (3 × 1,000 units per prime) | 84.87 |
| χ for primes ≤ 2^82 (100 units = 2·10⁵ instructions each: exponentiation plus Pohlig–Hellman) | 83.148 |
| prime sieve to X1 (64 instructions per integer); χ-sieve and prime class lists to X2 | 77.035 |
| collision checks (6 candidates) | 14.136 |
| d-th root (≤ 2^20 F_r mults of 200 instructions) | 16.678 |
| μ_d search (2^20 scalar mults) + final check | 31.551 |
| input decoding and subgroup checks | 13.551 |
| **total W(2^20)** | **116.958777** |

- A scalar multiplication is charged at 1,000 units for a ≤ 82-bit scalar and 3,000 units for a 255-bit
  scalar. These are generous Jacobian ladders at ≤ 5 units per step, plus one inversion at 100 units.
- The score is rounded up: 116.958777 = 116958777/1000000. K = 0.55 (rather than 0.52) costs 0.040 of
  score because every term scales as √K, and keeps both floors of §3.7 above 0.39 with margin.

**Memory.**
- Hash table: 2N slots × 72 bytes ≈ 2^122.9 bytes.
- Prime class lists: π(2^70) ≈ 2^65.03 primes × 9 bytes ≈ 2^68.7 bytes.
- Walk state and affine tables for 2^16 primes in flight: 2^16 · 2^34 · 96 bytes ≈ 2^56.6 bytes.
- Memory is not scored. Any BSGS form of Cheon stores one side, about 2^115.7 entries.
- Rho-type Cheon variants need little memory, but each of their steps is a multiplication in the
  exponent (≥ 7 additions with KKM tables).

## 5. W(n) and applicability

The algorithm uses only A_1 and A_{2^20}, so W(n) = 2^116.958777 for every 2^20 ≤ n < r. Nothing is
claimed for n < 2^20. More generally, for any d | r−1 with 6 | (r−1)/d and d ≤ n:
- W ≈ 2·c·√(K(r−1)/(6d)), plus the lower-order terms above;
- X2 = d·2^50.

## 6. Context (unscored estimates)

| Method | log₂ W |
|---|---:|
| Pollard rho on G1, order-6 classes: 0.5·log₂(πr/12) additions, plus overhead | ≥ 126.5 |
| Cheon BSGS, double-and-add scalar multiplication per step (≈ 383 additions) | ≈ 125.2 |
| Cheon BSGS with KKM c = 3 tables (2 additions plus one F_r multiplication per step); same K and per-element budgets as above, computed by `evidence.py` | 117.897 |
| **this work (v29)** | **116.959** |

- The gain over KKM-accelerated Cheon (recomputed under the same budgets) is a constant factor, about
  0.94 bit.
- **Relation to the generic bound.** Boneh–Boyen's generic lower bound Ω(√(ε·r/d)) has unstated
  constants. Here the walks use 2N ≈ 2^116.7 additions, and the generic coverage is
  N²/m = K/6 ≈ 0.092. The factor 6 that lifts success to 1 − e^{−K} comes from the cheap action of U,
  through negation and the GLV endomorphism. The model does not capture that action, so there is no
  conflict with the bound.

**Disclosed premises and heuristics:**
- **P1** (§3.5): the emitted-mass class second moments δ_s^w ≤ 500 —
  measured 0.678–8.49 at all three scales (59–737× headroom); the input to
  every stratum certificate of §3.6 (the class masses W_s(c) enter the
  lattice and Selberg weights).
- **The joint defect** t = (Δ₀ + R_SS)/B₀ ≤ 0.0138998 (§3.6) — THE single
  named unproved estimate behind the claim (route A): the coprime modular
  discrepancy plus the signed Selberg remainder, after every structured
  stratum is certified at production under P1 with no H2.
- **P3, P4** (route B only): the pruned factorial-moment bounds
  F₃ ≥ 0.95·κ³ and F₄ ≤ 1.20·κ⁴ — measured 0.959/1.024 and 0.900/0.976 at
  the two energy scales; not load-bearing for the claim, not proved.
- **H2** (supply only): the primes ≤ 2^70 and the primes of (2^81, 2^82]
  equidistribute over the χ-classes (each class at least X2/(2d·49) with
  ln-slack; proven 2.63× interval availability) — this affects how many
  primes are needed, not success; NO H2 clause enters any §3.6 stratum
  bound or the pruning.
- Expected-cost linear-probing hashing over a sharded table (evidenced by
  the probe/shard check).

No other assumptions.
