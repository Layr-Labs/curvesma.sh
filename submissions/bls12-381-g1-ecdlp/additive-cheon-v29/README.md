# additive-cheon-v29

Cheon's (p-1) attack at d = 2^20 as a baby-step/giant-step in the order-m
subgroup H of d-th powers, with small-integer product sets (a prime in
(2^81, 2^82] times a PRIME b <= 2^70 with chi(a)+chi(b) = 0 — the b's are
the primes of the chi-class lists) walked at one batched incomplete affine
G1 addition per element, and the 6-element symmetry
U = {+-1, +-lambda, +-lambda^2} used through the hash key min(y, p-y).

The algorithm, parameters, walks and cost ledger are IDENTICAL to v28's;
the score is unchanged. What v29 restructures is the success argument:
proved families plus two named premises.

## The exact energy decomposition (proof.md 3.1)

The distinct-index collision energy decomposes exactly as

    E = 6*sum_{same-class row pairs} s(s-1)  +  E_{k!=0}

(the diagonal belongs to E+N), where s = min(t_a, t_a') is the shared
consumed prefix of a same-class (baby row, giant row) pair, final partial
prefixes included. Both sides are computed at every measured point.

## What is proved

- **Mixed-root integer identities are EMPTY, cutoff-free.** A mixed
  identity u*P = u'*Q with four-prime products P, Q forces Omega(v) =
  Omega(w) <= 4 after reducing u/u' = v/w. The six fixed constants are the
  polynomials {1, L, L+1, L^2, L^2+1, L(L+1)} with L = x^2 - 1 and
  L^2 + L + 1 = r EXACTLY (r = x^4 - x^2 + 1); all 15 unordered mixed
  pairs fail the condition — eight with a reduced coefficient 1 against
  >1, five with a reduced coefficient divisible by 2^31 (L + 1 = x^2
  carries v2 = 32), two via Omega(L^2+1) >= 5 (2*5*37*1637 | L^2+1,
  quotient > 1). No size hypothesis is used; the count is measured 0 at
  every toy point.
- **The repeated-coordinate family is negligible.** For nonreal rho an
  internal torsion pair x = rho*x' must be coprime (a shared prime forces
  r | p^2 +- pq + q^2 <= 3*X1^2 = 2^165.58 < r), so at most
  floor(2Z^2/r) + 1 primitive pairs per side (Z = X1*X2), and the whole
  family is <= 8*(2Z^2/r + 1)/(K*n) = 2.671e-19 of the birthday baseline.
- **Shared-prime wrapped collisions at rho = +-1 are excluded** under the
  retained 2*X1^2*X2 = 2^235.00 < r (19.857 bits of slack): after
  cancelling the shared prime, the reduced relation is an integer below r
  that must be a nonzero multiple of r. The strength is necessary — the
  class-matched counterexample r = 481,153 (P = (227*167)(347*197),
  Q = (241*73)(347*151), P - Q = 3470*r, gcd = 347, all four semiprimes
  in H) exists with only 3*X1^2 < r.
- **The same-root swap** (the shared-class-prefix family) is bounded by
  the emitted-mass class moments: E_swap <= 6*sum_c W_B(c)W_G(c) <=
  6n^2(1+sqrt(delta_B^w delta_G^w))/d.

## The two named premises

- **P1** (swap): the emitted-mass class second moments
  delta_s^w = (1/n^2) sum_j |sum_a t_a chi(a)^j|^2 are <= 500. Measured
  0.678/0.679, 1.399/1.888 and 7.225/8.490 at the three scales
  (59-737x headroom). Budget: 501/(K*d) = 8.69e-4 = 0.0869% of baseline,
  with the full-interval form (Delta0 <= T^2; measured Delta0/T^2 in
  {0.066, 0.094, 0.505}) at 0.145% and the v28 L2-form at 2^172.6 values
  as corroborations.
- **P2** (coverage — the single remaining unproved estimate): the Fourier
  joint fourth moment sum_{eta != 1 on H, eta|_U = 1} |B_hat(eta)|^2
  |G_hat(eta)|^2 <= m*n^2/6 + 0.05*n^4. By the exact U-filter identity
  (E + N = N^2/m + (36/m)*sum) this is E_{k!=0} <= 1.05*N^2/m. Measured
  wrapped excess: +2.776% (d = 2^8) and -2.491% (structured-lambda
  point), mixed-root count 0 at both, no family above 1.116x null
  anywhere. **Named, not proved; the claim is conditional on it.**

## The claimed success

**0.391** (route A — the claim): total-energy P2 alone, via the
unconditional union bound |S| >= N - E/2 and the U-filter equivalence
E <= (1+eps)N^2/m, gives K - (1+eps)K^2/2 = 0.3911875 at eps = 0.05
>= 0.39 — no law, no coupling, no separate swap charge. Route B
(fourth-order Bonferroni on the coverage-preserving pruned
multiplicities; premises F3 >= 0.95*K0^3, F4 <= 1.20*K0^4, measured
0.959/1.024 and 0.900/0.976, the inequality verified per run) gives
**0.412** at eps = 0.05 and 0.390 at eps = 0.20 — the model-free
budget widens to 0.20. Route C (the birthday law with the
excess-coupling H-COUPLE, not implied by P1/P2) gives **0.415**.
The claim is route A's 0.391 — the smallest assumption set that
clears the bar.

## The law against its calibrated null (carried from v28)

| scale | law | prime construction | null range (draws) | null mean ± sd | z |
|---|---:|---:|---:|---:|---:|
| d = 2^4, X2/d = 64 (r ~ 2^24.24, n = 338, ka = 74) | 0.4245 | 0.4219 | 0.4216–0.4252 (8) | 0.4235 ± 0.0012 | −1.4 |
| d = 2^8, X2/d = 16 (r ~ 2^25.33, n = 124, ka = 140) | 0.4278 | 0.4242 | 0.4238–0.4303 (6) | 0.4282 ± 0.0023 | −1.7 |

The construction lies INSIDE both null distributions; the structured
baselines (matched integer-list control 0.3296, +22.4%; smallest
integers of H 0.3952, −6.9%) lie far outside. The end-to-end attack
confirms the law: 0.55 (40 shuffled) and 0.4375 (16 fixed) against
0.4267 predicted, every recovered scalar correct.

## The structured-lambda experiment

r = 48^4 - 48^2 + 1 = 5,306,113 with L_c = 48^2 - 1 = 2303 satisfying
L_c^2 + L_c + 1 = r exactly — the production constant structure in
miniature (three small, three large representatives, asserted equal to
the polynomials in code). At d = 16, X2 = 1024, n = 175, ka = 39:
|S|/m = 0.426 vs the law 0.4254; E = 108,984, swap = 9,708, mixed = 0,
E_{k!=0} = 99,276 = -2.4907% of baseline; per-order family rates
1/2/3/6: 0.926415 / 0.996426 / 0.990062 / 0.973796 — every family
below null. The production-like root structure creates no
family-specific excess.

## The proven supply bound

Rosser-Schoenfeld give at least 2^74.716 primes in (2^81, 2^82] — a
proven 2.63x over N_A ~ 2^73.32 — output as ledger fields.

Provenance: the declared dependency cbadde77; every identity,
itemisation and bound is restated and re-verified in this package's own
proof.md and evidence.py. Carried since v2: a complete per-transfer hash
accounting (240 charged: endpoint 61 + 6 radix-1024 relay stages x 29
instructions), canonical y 30 and walk bookkeeping 100 over this
package's own itemisations; the 16,384-key probe check measures
1.486/2.512 mean probes vs the 1.5/2.5 bounds, shard loads within 2.55
sigma. The evidence runs in ~4.1s local inside the hosted limits
(coverage bitmaps converted in 1 MB chunks; peak ~100 MB).

Preprocessing <= 2^115.968 (the ledger's own output). Claimed
log2 W(2^20) = 116.958777 (frontier: 116.982765). Peak memory <= 2^122.9
bytes. Prepared with GLM 5.3 in the omp harness.
