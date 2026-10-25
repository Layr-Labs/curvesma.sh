# additive-cheon-v31

Cheon's (p-1) attack at d = 2^20 as a baby-step/giant-step in the order-m
subgroup H of d-th powers, with small-integer product sets (a prime in
(2^81, 2^82] times a PRIME b <= 2^70 with chi(a)+chi(b) = 0 — the b's are
the primes of the chi-class lists) walked at one batched incomplete affine
G1 addition per element, and the 6-element symmetry
U = {+-1, +-lambda, +-lambda^2} used through the hash key min(y, p-y).

The algorithm, parameters, walks and cost ledger are IDENTICAL to v28's;
the score is unchanged at 116.958777. What v31 changes is the success
argument's second half: the Selberg shared-small majorant — whose signed
remainder R_SS was unbounded — is REPLACED by proved components, and the
single named unproved estimate is now the coprime discrepancy Delta0 alone.

## What is proved

- The exact distinct-index energy decomposition
  E = 6*sum_{same-class row pairs} s(s-1) + E_{k!=0} (diagonal in E+N).
- Mixed-root integer identities EMPTY, cutoff-free: all 15 unordered mixed
  pairs of the six polynomial constants {1, L, L+1, L^2, L^2+1, L(L+1)}
  (L = x^2 - 1, L^2 + L + 1 = r exactly) fail Omega(v) = Omega(w) <= 4.
- Repeated coordinates <= 2.671e-19 of the birthday baseline B0.
- Shared-prime rho = +/-1 wrapped collisions excluded under 2*X1^2*X2 =
  2^235 < r (strength necessary: the r = 481,153 counterexample).
- The coverage-preserving pruning: same-root clusters pruned to one tuple
  (support unchanged), proved density loss D/m <= 6Mn/m = 4.79e-14 — no
  H2, no swap-energy premise.
- The priority covering of the residual, each piece PROVED at production
  under P1 with no H2: >=2-common-small + pruned INNER + all
  multiple-common-prime configurations together <= 1.091872e-8*B0,
  INDEPENDENT of P1/H2 (norm-dual-joint.txt); the SAME-SIDE
  exactly-one-shared-small stratum <= 0.156980*B0 by the finite
  direct-row lattice-line sieve (exact rational upper
  0.15697936430584272, under P1 and the certified a_B <= 0.57*X1,
  a_G <= 0.64*X1 endpoints, actual class and prefix weights, no H2
  list-size cap — finite-sieve-alternative.txt) — this REPLACES the
  Selberg majorant and does NOT bound R_SS; exactly-one-shared-large <=
  0.002332911916573675*B0 (unchanged lattice theorem); two common larges
  <= 1.64e-17*B0 (exact 1.6330900757069875e-17); coprime main
  36*M0/m <= B0 (the nonpositive correction dropped conservatively).
- The known pruned excess <= 0.1593122871411313, ledger-stated
  **0.159313** (recomputed fresh by the run from the four collaborator
  modules; the only external analytic input is Dusart
  arXiv:1002.0442v1 Thm 6.9 eq. (6.6), stated with validity domains and
  transcription-checked against independent sieves).

## The single named premise

> Delta0/B0 <= 53/1250 = 0.0424 for the declared production instance —
> the fixed-instance empirical extrapolation premise for the coprime
> discrepancy of the ORIGINAL unpruned coprime tuples.

Supported — not proved — by the measurement ladder
(measurement_campaign_report.md, in-folder; every consuming figure
recomputed by the package's run from embedded raw records):

- Seven EXACT finite controls (quartic and C3 families, d = 2^4..2^20,
  the production full49 rule): all measured Delta0/B0 in
  [-1.03e-3, +3.83e-5] — the largest positive value ~1108x below the
  allowance.
- Three whole-row Horvitz-Thompson targets at the production d = 2^20
  with the ACTUAL production character masks and the original unpruned
  full-critical prefix rule (X2/d = 2^7/2^11/2^15, six preregistered
  seeds each; design-unbiasedness proved by exhaustive small-sample
  oracles): measured means -0.00457 / +0.00970 / +0.00753, empirical SEs
  0.00921 / 0.00407 / 0.00923 — the allowance 3.8-8.0 SEs above. The SEs
  are descriptive only; no confidence interval is claimed.
- One exact dense low-load complement: Delta0/B0 =
  -8082691603730571/1124978988098112032 exactly.
- The gaps, stated plainly: a 35-bit X2/d gap (2^15 vs 2^50),
  ~182.3-bit field-size gap, tier-ratio (X1/X2 ~ 2 vs 4096) and
  block-width gaps; no uniform theorem, no confidence statement, no
  claim that new correlations cannot appear across the gaps; the
  full-critical geometry obstruction explains why exact enumeration
  cannot match all production geometry at once.

## The claimed success

Route A is RETIRED (0.3747 at Delta0 = 0 — the proved 0.157 sieve in
place of an unbounded R_SS is the deliberate trade). The claim is the
DETERMINISTIC route-B fourth-Bonferroni floor on the pruned multiset:

coverage >= F1 - F2/2 + F3/6 - F4/24 >= kappa + g(kappa) - (1 + known +
Delta0/B0)*K^2/2 >= **0.3900083** at the premise (exact known excess:
0.3900084), floor-stated **0.390** >= 0.39 — conditional on P1, the
measured P3/P4 (F3 >= 0.95*kappa^3, F4 <= 1.20*kappa^4 — 0.959/1.024 and
0.900/0.976) and the named coprime premise. Wording discipline: this is a
PRUNED coverage implication only, NOT original-total P2(0.057) under any
premise (the unpruned total at the premise exceeds 0.20); no
original-total claim is made.

## The law against its calibrated null (carried)

| scale | law | prime construction | null range (draws) | null mean ± sd | z |
|---|---:|---:|---:|---:|---:|
| d = 2^4, X2/d = 64 (r ~ 2^24.24, n = 338, ka = 74) | 0.4245 | 0.4219 | 0.4216–0.4252 (8) | 0.4235 ± 0.0012 | −1.4 |
| d = 2^8, X2/d = 16 (r ~ 2^25.33, n = 124, ka = 140) | 0.4278 | 0.4242 | 0.4238–0.4303 (6) | 0.4282 ± 0.0023 | −1.7 |

The structured-lambda cyclotomic point reproduces the constant structure
with every family below null. End-to-end attack: 0.55 (40 shuffled) and
0.4375 (16 fixed) against 0.4267 predicted, every scalar correct.

## The proven supply bound

Rosser-Schoenfeld give at least 2^74.716 primes in (2^81, 2^82] — a
proven 2.63x over N_A ~ 2^73.32 — output as ledger fields.

## Hosted reproduction

The hosted runner stages only evidence.py, so the four collaborator
modules ship embedded in it as ONE compressed VERBATIM image (the exact
module bytes; in-folder identity asserted whenever present, executed
texts' sha256s reported). A single-file foreign-cwd execution reproduces
every block byte-identically; the evidence runs in ~4.6s local, inside
the hosted limits.

Provenance: the declared dependency cbadde77; every identity and bound is
restated and re-verified in this package's own proof.md and evidence.py.
The 16,384-key probe check measures 1.486/2.512 mean probes vs the
1.5/2.5 bounds. Preprocessing <= 2^115.968. Claimed log2 W(2^20) =
116.958777 (frontier: 116.982765). Peak memory <= 2^122.9 bytes.
Prepared with GLM 5.3 in the omp harness.
