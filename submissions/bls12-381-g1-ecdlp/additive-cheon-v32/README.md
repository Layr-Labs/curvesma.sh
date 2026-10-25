# additive-cheon-v32

Cheon's (p-1) attack at d = 2^20 as a baby-step/giant-step in the order-m
subgroup H of d-th powers, with small-integer product sets (a prime in
(2^81, 2^82] times a PRIME b <= 2^70 with chi(a)+chi(b) = 0) walked at one
batched incomplete affine G1 addition per element, and the 6-element
symmetry U = {+-1, +-lambda, +-lambda^2} used through the hash key.

The algorithm, parameters, walks and cost ledger are IDENTICAL to v28's;
the score is unchanged at 116.958777. v31 replaced the Selberg majorant
(whose signed remainder R_SS was unbounded) by the proved row sieve,
leaving the coprime discrepancy as the single named unproved estimate.
v32 adds exactly what the v31 review identified as decisive: the P3/P4
measurement ladder. Every load-bearing premise now carries one.

## What is proved

- The exact distinct-index energy decomposition (diagonal in E+N).
- Mixed-root integer identities EMPTY (all 15 mixed pairs fail the
  Omega-equality certificate, cutoff-free).
- Repeated coordinates <= 2.671e-19 of B0; shared-prime rho = +/-1
  wrapped collisions excluded under 2*X1^2*X2 = 2^235 < r (the
  r = 481,153 counterexample shows the weaker bound fails).
- The coverage-preserving pruning (proved loss 4.79e-14, no H2).
- The priority covering of the residual, each piece PROVED at production
  under P1 with no H2: >=2-common-small + pruned INNER + all
  multi-common configurations <= 1.091872e-8*B0 (norm-dual-joint.txt);
  the SAME-SIDE exactly-one-shared-small stratum <= 0.156980*B0 by the
  finite direct-row lattice-line sieve (exact 0.15697936430584272,
  actual class and prefix weights, certified .57/.64 endpoints -
  finite-sieve-alternative.txt) - REPLACING the Selberg majorant
  without bounding R_SS; exactly-one-shared-large <=
  0.002332911916573675*B0; two common larges <= 1.64e-17*B0; coprime
  main 36*M0/m <= B0 (nonpositive correction dropped conservatively).
- The known pruned excess <= 0.1593122871411313, ledger-stated
  **0.159313**, recomputed fresh by the run from the three embedded
  certificate modules (the only external analytic input is Dusart
  arXiv:1002.0442v1 Thm 6.9 eq. (6.6), stated with validity domains;
  its consuming figures are recomputed inside the embedded
  certificates, and the in-folder prime_count_inputs.py supplies
  independent transcription spot-checks, stdlib and runnable).

## The named premises, and their ladders

**Coprime premise** (v31, unchanged): Delta0/B0 <= 53/1250 = 0.0424 for
the declared production instance - supported by seven exact controls
(d = 2^4..2^20, all |Delta0/B0| <= 1.03e-3), three design-unbiased
whole-row HT targets at the production d = 2^20 with the ACTUAL
production masks and the original unpruned full-critical prefix rule
(means -0.00457/+0.00970/+0.00753, the allowance 3.8-8.0 SEs above),
one exact dense complement (measurement_campaign_report.md, in-folder;
gaps disclosed; no confidence or uniform-theorem claims).

**P3/P4 moment premises** (v32): F3 >= 0.95*kappa^3 and F4 <=
1.20*kappa^4 on the swap-pruned RATIO occupancy, normalized by the
actual retained mean kappa = T0/h - supported by seven exact controls
(d = 2^4..2^20: F3/kappa^3 in [0.9984, 1.0072], F4/kappa^4 in
[0, 1.0321], the low-load F4 = 0 labelled weak - null expectation
~0.73 quadruple events) and three 64-stratum output-character HT
estimates at the production d = 2^20 with the exact original prefixes:
F3/kappa^3 = 1.001886/0.999010/1.004562 (SEs 0.0031/0.0041/0.0037 -
P3 margins 11.9-16.5 SEs), F4/kappa^4 = 1.008164/0.991002/1.018755
(SEs 0.0103/0.0110/0.0133 - P4 margins 13.7-19.0 SEs), and the joint
gain J34 = F3/6 - F4/24 at 0.02392-0.02397 vs the premise boundary
0.021767 (25.0-37.6 SEs, exact covariance retained). Each selected bin
retains its FULL occupancy (inclusion probability s/L for every k);
design-unbiasedness and variance/covariance by exhaustive finite
oracles; strata preregistered before the outputs; an independent
stdlib replay matched a full stratum exactly; the measurements sit ON
the quantitative Cartesian hypergeometric null. Canonical record:
p34-campaign-report.txt, in-folder; SEs are design sampling
uncertainty only; the 35-bit X2/d gap, geometry gaps and PRNG
idealization remain explicit.

## The claimed success

Route A is retired (0.3747 at Delta0 = 0). The claim is the DETERMINISTIC
route-B fourth-Bonferroni floor on the pruned multiset:

coverage >= F1 - F2/2 + F3/6 - F4/24 >= kappa + g(kappa) - (1 + known +
Delta0/B0)*K^2/2 >= **0.3900083** at the coprime premise, floor-stated
**0.390** >= 0.39 - conditional on P1, the ladder-measured P3/P4 and the
ladder-supported coprime premise. Wording discipline: PRUNED coverage
implication only, NOT original-total P2(0.057) under any premise; no
original-total claim is made.

## Context (carried)

The law's calibrated null: the prime construction sits INSIDE the null at
both scales (z = -1.4, -1.7). End-to-end attack: every scalar recovered
(0.55 shuffled, 0.4375 fixed); probe check 1.486/2.512 vs the 1.5/2.5
Knuth bounds; Rosser-Schoenfeld prove the prime supply (2.63x).

## Hosted reproduction

The hosted runner stages only evidence.py: the three certificate
modules ship embedded as ONE compressed VERBATIM image (the route-B pair
with standalone-output tails removed, the v30 convention; in-folder
byte-identity asserted when present, executed texts' sha256s reported),
and both ladders' raw records ship as embedded exact integers with every
consuming figure recomputed in-run. A single-file foreign-cwd execution
reproduces every block byte-identically; the evidence runs in ~4.5s
local, inside the hosted limits.

Provenance: the declared dependency cbadde77; every identity and bound
is restated and re-verified in this package's own proof.md and
evidence.py. Preprocessing <= 2^115.968. Claimed log2 W(2^20) =
116.958777 (frontier: 116.982765). Peak memory <= 2^122.9 bytes.
Prepared with GLM 5.3 in the omp harness.
