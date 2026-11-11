# Individual-logarithm descent (capped special-q ladder)

After the factor-base logarithms are known (relation collection + linear algebra),
each target logarithm -- both log_g and log_g z -- is obtained by special-q descent.
`experiment.py` recomputes this ladder; every probability is the same Dickman rho
used for relation collection, and every node is a capped Las-Vegas loop.

## Phases

1. **Randomisation (boot).** Replace the target w (= g, then = z) by w*g^s for
   random s and lattice-reduce in the number field until the representative is an
   algebraic integer of norm ~ L_Q[2/3, 0.95] ~ 2^444 bits (`experiment.py`
   recomputes this from Q). Smooth it over primes up to a medium bound.
2. **Tree.** Each remaining large prime (special-q) is lattice-sieved into a few
   smaller primes; the prime size falls geometrically until every factor is in the
   factor base (<= B = 2^62.4).

## Ladder (recomputed by experiment.py, field F_{p^12})

| Level | Norm bits | u = bits/log2 B | P_level = rho(u) | nodes (log2) | level cost (log2) |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | 444 | 7.12 | 2^-20.7 | 0.0 | 26.0 |
| 1 | 276 | 4.41 | 2^-9.2 | 2.2 | 16.7 |
| 2 | 171 | 2.74 | 2^-3.6 | 4.4 | 13.3 |
| 3 | 106 | 1.70 | 2^-1.1 | 6.6 | 13.0 |
| 4 | 66 | 1.05 | 2^-0.08 | 8.8 | 14.2 |
| 5 | 64 | 1.03 | 2^-0.05 | 11.0 | 16.4 |

Per-log cost = log2( sum 2^level ) = **2^26.0**; both logs = **2^27.0**. The boot
level dominates (2^26); subsequent levels are cheaper because the prime size and
hence u fall faster than the node count grows. Each node is charged a full
reference unit per candidate (KERNELS.md convention), and the per-node trial count
is the cap c/P_level with c = 40 (over-provision).

## Why descent does not affect the claim

Descent is **sub-dominant by more than 90 bits** (2^27 vs the relation term
2^118.12). The chain is therefore insensitive to the descent model: even if the
true descent cost were as large as 2^115 -- far above any standard NFS descent
estimate for this field -- the chain
log2(2^118.12 + 2^117.39 + 2^115) would rise by < 0.5 bit, which the 4.5-bit margin
already covers. We give the explicit ladder for completeness, not because the claim
depends on its exact value.

## Success probability

Each node runs its Las-Vegas loop up to cap c/P_level with c = 40, so a node fails
with probability at most e^-40. With ~2^13 nodes total the expected number of
failed nodes per log is ~2^13 * e^-40 ~ 10^-14, and the per-target success is
~ 1 - 2*10^-14 for the two logs. The >= 0.39 requirement holds for any per-log
descent failure up to 0.30 (success >= 0.40); see proof.md section 7. Heuristic H5
(descent independence; Commeine-Semaev 2006, Barbulescu 2013) is disclosed in
proof.md section 8.
