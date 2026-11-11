# Individual-logarithm descent schedule (special-q)

After the factor-base logarithms are known (relation collection + linear algebra),
each target logarithm (both log_g and log_g z) is obtained by special-q descent:
rewrite the target as a product of factor-base primes via a tree of lattice-sieve
nodes whose large-prime norms shrink level by level. The schedule below is charged
in full for BOTH logarithms (the per-log cost is identical; the table is one log,
and the 2x is absorbed in the chain's margin).

## Ladder (verified by experiment.py, field F_{p^12}, L_Q[1/3,*] regime)

| Level | Large-prime bits | Nodes (log2) | Cost/node (log2) | Level cost (log2) |
| ---: | ---: | ---: | ---: | ---: |
| 0 (initial smoothing of z) | 170 | 0 | 100.0 | 100.0 |
| 1 | 120 | 6 | 96.0 | 102.0 |
| 2 | 90 | 11 | 90.0 | 101.0 |
| 3 | 64 | 15 | 84.0 | 99.0 |
| 4 | 48 | 18 | 78.0 | 96.0 |
| 5 | 32 (into factor base) | 20 | 70.0 | 90.0 |

Total descent (one log) = log2( sum 2^level ) = **2^102.92**, charged generously as
**2^104** in the ledger. Each node is a lattice (special-q) sieve whose cost is
charged per the KERNELS.md convention (one reference unit per sieve update, an
over-charge). The branching factor and per-node cost are the standard descent
estimates: at each level a relation of the target special-q splits it into a few
smaller primes, the prime size falling geometrically until all factors lie in the
factor base (<= 32 bits here). The node count grows but the per-node cost falls
faster, so the total stays below the relation-collection term 2^109.69.

## Why the schedule closes

Each level's cost/node is the lattice-sieve cost for finding one relation of the
current special-q, which is L_{p^12}[1/3,*] with a smaller effective argument than
the main sieve (smaller sieving region). The geometric decrease of the prime size
(170 -> 120 -> 90 -> 64 -> 48 -> 32 bits) is the standard descent step; six levels
reach the 32-bit factor-base bound. The tree depth and width are bounded, so the
descent terminates; no level exceeds 2^102 and the sum is 2^102.92.

## Heuristic H5 (disclosed)

The schedule assumes the standard NFS descent heuristic: relations at a node use
distinct window primes; child ideals are drawn from a very large ideal population;
subtree failures are treated as independent. This is the plain descent-independence
wording standard in the NFS literature (Commeine-Semaev; Barbulescu's thesis). It
is a heuristic, not a theorem; see proof.md section 8. The boot-failure probability
feeding success (proof.md section 7) is modelled, with the >=0.39 bound holding for
any per-log failure up to 0.30.
