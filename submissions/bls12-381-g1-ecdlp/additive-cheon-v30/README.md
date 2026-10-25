# additive-cheon-v30

Cheon's (p-1) attack at d = 2^20 as a baby-step/giant-step in the order-m
subgroup H of d-th powers, with small-integer product sets (a prime in
(2^81, 2^82] times a PRIME b <= 2^70 with chi(a)+chi(b) = 0 — the b's are
the primes of the chi-class lists) walked at one batched incomplete affine
G1 addition per element, and the 6-element symmetry
U = {+-1, +-lambda, +-lambda^2} used through the hash key min(y, p-y).

The algorithm, parameters, walks and cost ledger are IDENTICAL to v28's;
the score is unchanged at 116.958777. What v30 adds is the residual's
CERTIFIED STRATA: every structured stratum of the wrapped-collision
residual is now PROVED at production scale under P1 with no H2 clause,
leaving one named joint-defect bound.

## What is proved

- The exact distinct-index energy decomposition
  E = 6*sum_{same-class row pairs} s(s-1) + E_{k!=0} (diagonal in E+N).
- Mixed-root integer identities EMPTY, cutoff-free: all 15 unordered mixed
  pairs of the six polynomial constants {1, L, L+1, L^2, L^2+1, L(L+1)}
  (L = x^2 - 1, L^2 + L + 1 = r exactly) fail the necessary condition
  Omega(v) = Omega(w) <= 4 — eight with a reduced coefficient 1 against
  >1, five with a reduced coefficient divisible by 2^31, two via
  Omega(L^2+1) >= 5 (2*5*37*1637 | L^2+1, quotient > 1).
- Repeated-coordinate family <= 8*(2Z^2/r + 1)/(K*n) = 2.671e-19 of the
  birthday baseline B0 (Z = X1*X2).
- Shared-prime wrapped collisions at rho = +/-1 excluded under the
  retained 2*X1^2*X2 = 2^235 < r (19.857 bits), whose strength is
  necessary: the class-matched counterexample r = 481,153 (P - Q = 3470*r,
  gcd(P,Q) = 347) exists with only 3*X1^2 < r.
- The coverage-preserving pruning: same-root shared-prefix clusters
  pruned to one tuple (support unchanged), with proved density loss
  D/m <= 6Mn/m = 4.79e-14 — no H2, no swap-energy premise (the 8.69e-4
  swap charge is discharged).
- The residual's strata, each certified at production under P1 with no
  H2: exactly-one-shared-large <= 0.002332911916573675*B0 (the
  affine-coset lattice theorem: mod-30 wheel factor 16/225, the weighted
  short-vector step F_s(T) <= ||W_s||_2*sqrt(sum_c k_c(T)) preserving the
  actual P1 class masses, explicit root-lattice cells with <= 9
  coefficient antecedents, Dusart-certified endpoints a_B <= 0.57*X1,
  a_G <= 0.64*X1); shared-small Selberg main <= 0.038976659226098*B0
  (finite optimal Selberg weights at R = 2^35; the denominator
  G(2^35) > 322.7629594 certified by an integer hyperbola computation,
  no analytic input); small-prime exception <= 0.002641842*B0
  (S35 = 2,086,667,226); two-large-only <= 1.633e-17; coprime main
  certified to 3.492e-13.
- The five pieces sum to alpha = 0.043951414 of B0 (a directed 9-place
  ceiling over exact rationals, recomputed fresh by the run from the
  three in-folder source modules; the only external analytic input is
  Dusart arXiv:1002.0442v1 Thm 6.9 eq (6.6), stated with validity
  domains and transcription-checked against independent sieves).

## The single named premise

> t := (Delta_0 + R_SS)/B0 <= 0.0138998 — the coprime modular discrepancy
> plus the signed Selberg remainder.

The U-filter identity (E + N = N^2/m + (36/m)*sum over the ~m/6 - 1
U-invariant nontrivial characters) is verified EXACTLY at both energy
points of the run; the self-moment route is blocked exactly as measured
(the diagonal alone forces S_B, S_G >= h(2n^2 - n) - n^4, so
Cauchy-Schwarz cannot reach below 1/K - 1 = 0.818; measured ratios
1.985/2.064), which is why the certificates control the mixed correlation
stratum by stratum. **t is named, not proved.**

## The claimed success

Both routes are DETERMINISTIC coverage implications — no birthday law,
no coupling, no Poisson approximation:

- **Route A (the claim)**: coverage >= F1 - F2/2 >= kappa - (1+alpha+t)*K^2/2
  >= **0.39000000388245215** at t <= 0.0138998, floor-stated **0.390**
  >= 0.39 (kappa = K - 4.79e-14; exact maximal t = 0.013899825669104974;
  at the sharper t = 0.0121794: >= 0.390260).
- **Route B**: additionally under P3 (F3 >= 0.95*kappa^3) and P4
  (F4 <= 1.20*kappa^4) — measured 0.959/1.024 and 0.900/0.976 —
  t <= 0.157816 gives >= 0.39000007446578, t <= 0.10492393035226562
  gives >= 0.398, and t = 0.0121794 gives >= 0.4120276.
- **Wording discipline**: the pruned implications are NOT the original
  total-energy P2(0.057) — unpruned, epsilon_total <= alpha + t + beta =
  0.058720214 at the safe gap (beta = 8.69e-4, P1's swap charge), under
  0.057 only at the sharper t.

## The law against its calibrated null (carried)

| scale | law | prime construction | null range (draws) | null mean ± sd | z |
|---|---:|---:|---:|---:|---:|
| d = 2^4, X2/d = 64 (r ~ 2^24.24, n = 338, ka = 74) | 0.4245 | 0.4219 | 0.4216–0.4252 (8) | 0.4235 ± 0.0012 | −1.4 |
| d = 2^8, X2/d = 16 (r ~ 2^25.33, n = 124, ka = 140) | 0.4278 | 0.4242 | 0.4238–0.4303 (6) | 0.4282 ± 0.0023 | −1.7 |

The structured-lambda cyclotomic point (r = 48^4 - 48^2 + 1, whose six
constants are the production polynomials in miniature) reproduces the
constant structure with every family below null and |S|/m = 0.426 vs the
law's 0.4254. The end-to-end attack: 0.55 (40 shuffled) and 0.4375 (16
fixed) against 0.4267 predicted, every recovered scalar correct.

## The proven supply bound

Rosser-Schoenfeld give at least 2^74.716 primes in (2^81, 2^82] — a
proven 2.63x over N_A ~ 2^73.32 — output as ledger fields.

Provenance: the declared dependency cbadde77; every identity, itemisation
and bound is restated and re-verified in this package's own proof.md and
evidence.py. The stratum certificates ship as in-folder source modules
executed fresh by the run, with their full derivations in-folder
(shared-large-derivation.txt, selberg-derivation.txt,
prime-count-dependency.txt). The 16,384-key probe check measures
1.486/2.512 mean probes vs the 1.5/2.5 bounds, shard loads within 2.55
sigma. The evidence runs in ~4.9s local inside the hosted limits.

Preprocessing <= 2^115.968 (the ledger's own output). Claimed
log2 W(2^20) = 116.958777 (frontier: 116.982765). Peak memory <= 2^122.9
bytes. Prepared with GLM 5.3 in the omp harness.
