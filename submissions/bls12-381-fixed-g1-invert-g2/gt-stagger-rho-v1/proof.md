# Fixed-argument pairing inversion on BLS12-381 G1: staggered-checkpoint rho in GT

Track `bls12-381-fixed-g1-invert-g2`. Candidate `gt-stagger-rho-v1`.
Claim **127.67** (`12767/100`), success at least **0.39**, peak memory 656,704 B.

Everything numeric below is recomputed by `experiment.py` (standard library only)
from the declared inputs; the offline tools named in §9 reproduce the model
numbers. Nothing here is copied from an earlier candidate.

## 1. Reduction

Instance: nonidentity `P` in G1 and `z` in GT, with `r = |G1| = |G2| = |GT|`
prime. Choose a public fixed nonidentity `Q0` in G2 (for example the pinned
generator), and compute

    g = e(P, Q0) = pairing(Q0, P)          (the pinned library takes (Q, P) order)

`e` is nondegenerate and bilinear, `r` is prime and `P != O`, so `g != 1` and
`<g> = GT` (a group of prime order has no proper nontrivial subgroup). Hence
there is a unique `a` in `Z/r` with `z = g^a`, and

    e(P, [a]Q0) = e(P,Q0)^a = g^a = z,

which solves the challenge for that `P` and that `z`, including `z = 1`
(`a = 0`, `Q = O`). The problem is therefore exactly a discrete logarithm in the
target group, with a single instance and no auxiliary powers: what is left to
price is the DLP.

Verification after the solve is one pairing `e([a]Q0, P) =? z` plus one G2 scalar
multiplication; both are charged (§6).

## 2. The walk

Public jump table `Th[h] = g^(A_h) z^(B_h)`, `h = 0..R-1`, `R = 1024`, built once
from the instance (2 x 1024 scalar multiplications in GT). State `(X, a, b)` with
`X = g^a z^b`, starting from a fresh uniform `(a0, b0)`:

    h  = H(X)                    one 64-bit limb of the canonical Fp12 encoding,
                                 fixed odd multiplier, xor-shift, top 10 bits
    X <- X * Th[h]               one dense Fp12 multiplication = 1 reference unit
    a <- a + A_h mod r,  b <- b + B_h mod r

The next state is a function of the group element `X` alone. The trajectory is
therefore a deterministic map `X -> X * Th[H(X)]` on a set of size `r`, and it is
modelled as a random map on a cyclic group of prime order `r`. There is no
history, no value lookup, no distinguished-point set and no dictionary: the only
sequence data kept is the small checkpoint window of §3.

## 3. Detection: staggered checkpoints, and the working set is *counted*

The predecessor candidates charged `2*S0` evaluations on the strength of the
sentence "Brent evaluates the map at most `2*rho+3` times". That sentence is
false (a tail-0 cycle of length 9 is detected at evaluation 24) and it is the
single defect named in this account's three rho rejections. It is deleted, not
patched: this candidate charges no cycle-length bound at all.

The detector is a staggered checkpoint schedule. Let

    P_k = { ceil(2^j (k+i)/k) : j >= 0, 0 <= i < k }        k = 4

so that `P_4` has exactly 4 points in every dyadic level `[2^j, 2^(j+1))`
(`{2^j, 2^j+2^(j-2), 2^j+2^(j-1), 2^j+3*2^(j-2)}` for `j >= 2`). A checkpoint
created at position `p` stores `(word(X), a, b)`; it is compared with the hare
only at hare positions in the **time window** `[p+1, W*p]`, `W = 8`, and is then
dropped. Expiry is one test on the oldest entry, because the live list is in
increasing `p`.

**Working-set law.** A window of multiplicative span `W = 2^b` straddles exactly
`b` dyadic levels, and each level carries `k` points, so

    n_max(k, W) = k * log2(W)          (W a power of two)

and for `(k,W) = (4,8)` that is **12** checkpoints, each 112 B (8-byte word,
48-byte `a`, 48-byte `b`, 8-byte position). This is the correction that matters:
the rejected `gt-multicp-brent-v2` priced `k*(1-1/W)` — 3.5 states here, 56 at
`k=64,W=8` where the true count is 192. `experiment.py` §B1 therefore does not
assert the number, it **enumerates** it: it walks every hare position in
`[1, 2^20]`, takes the maximum (12), and reports the count at the three probe
positions 1024, 4096 and 2^20 (12, 12, 12). The old `k*(1-1/W)` value is printed
beside the enumerated one so the two cannot be confused again.

**Detection index.** For a trajectory with tail `mu` and cycle `lam`, the states at
positions `c < m` are equal iff `m = c` or both `c >= mu` and `m = c (mod lam)`.
A checkpoint stored at `p` and compared with the hare at `p+1 .. W*p` is matched at

    m = p + lam,   first for p = p* = min{ p in P_k : p >= max(mu, ceil(lam/(W-1))) }

because the checkpoint must lie on the cycle (`p >= mu`) and its window must still
be open when the hare reaches `p + lam` (`lam <= (W-1) p`). `m` is the number of
map evaluations before the match is seen, so it is the evaluation count.

Verification, at two independent scopes:

* inside `experiment.py` §B2, against a **direct simulation** of the comparison
  decisions (build the grid, walk the hare, test the live checkpoints) for every
  `(mu, lam)` with `mu + lam <= 60`: 1830 pairs, **0 mismatches**;
* offline, exhaustively over every pair with `mu <= 600, mu+lam <= 600` for
  `k=1,W=2` (180,300 pairs) and over 9 `(k,W)` schedules x 4,095 pairs:
  **0 mismatches**.

The theorem `m <= 3*rho` holds on the enumerated range and is attained
(`max m/rho = 2.0` for this schedule over the range of §B2, `2.9942` for `k=1` at
`rho = 513`), which is why charging `2*rho` was wrong.

## 4. Success probability: a computed quantile, not a work bound

Model (stated as a heuristic): `rho` satisfies `Pr[rho >= t] = exp(-t^2/(2r))` and
`mu` is uniform on `[0, rho)`, independent across instances and fresh coins.
Substituting the law of §3 into that model and integrating deterministically
(quadrature in `(rho, mu)`, no sampling error; cross-checked against 2,000,000
Monte Carlo draws of the identical law, agreement to 4 decimals) gives the
quantiles of `m/sqrt(r)` for `(k,W) = (4,8)`:

| quantile | 0.39 | 0.40 | 0.45 | 0.50 |
| --- | --- | --- | --- | --- |
| `c = Q_q(m)/sqrt(r)` | 1.0503 | 1.0682 | 1.1573 | 1.2430 |

**Phase robustness.** The grid is dyadic, so the rounding factor
`2^ceil(log2 x)/x` depends on the fractional part of `log2` of the position
scale. An earlier table in this session computed the model at `sqrt(N) = 2^35`
and so misplaced it by 0.43 in phase; that is why its Brent figure (1.537) does not
match this recomputation (1.437). Sweeping 20 phases, `c` at the 0.45 quantile is

| `k` | 1 | 2 | 4 | 8 | 16 | 64 |
| --- | --- | --- | --- | --- | --- | --- |
| phase spread of `c` | 9.0% | 2.4% | **0.8%** | 0.4% | 0.1% | 0.0% |

For `(4,8)` the 0.45 quantile lies in `[1.151, 1.159]` at every phase, so the cap
below does **not** depend on the lucky alignment of the production scale. The
candidate caps at `c = 1.16`, the top of that band, and declares only the 0.39
floor the track requires — a 0.06 absolute margin that also absorbs the model's
own error.

**Hard cap.** `B = ceil(1.16*sqrt(r)) = 265627019777073510555230149274057194865`
(exact integer square root, `B-1` verifiably fails). The run stops after `B`
evaluations and reports FAILURE. So

    P[success] = P[m <= B] = 0.45 >= 0.39.

No amplification or restart is charged, because none is needed: the requirement is
a single trial above 0.39 and the trial above is at 0.45.

**Toy validation of the model, at production `R`.** `experiment.py` §B3 runs the
complete algorithm — same grid, same window, same expiry, same solve, production
jump-table size `R = 1024` — on a 26-bit prime-order subgroup, 300 trials, cap
`ceil(1.16*sqrt(q)) = 6720`:

| quantity | value |
| --- | --- |
| solved inside the cap | 137/300 = **0.4567** (Wilson 95% [0.401, 0.513]) |
| solved + capped = trials | 300 = 300 |
| recovered log verified correct | 137/137 |
| bookkeeping failures | 0 |
| live-bound violations, max live observed | 0, **12** |
| degenerate relations | 0 |

The toy success fraction 0.4567 matches the 0.45 the model predicts for this
schedule, which is the check the earlier "aggregate evidence does not match the
declared cap" objection asked for. It is evidence at toy scope only.

## 5. Solve, degenerate draws, and the withdrawn argument

On a detected match, the two states are equal group elements with exponent pairs
`(a,b)` and `(ta,tb)`, so `(ta-a) + (tb-b) x = 0 mod r`. Hence

    den = (b - tb) mod r;   if den = 0 restart inside the same cap
    a   = (ta - a) * den^-1 mod r

and the value is exact because the elements are equal. `den = 0` means the two
exponent pairs coincide, i.e. the walk returned to the identical state. That
requires the product of the jumps over the gap to be the identity of GT:

    sum A_h = 0 AND sum B_h = 0  (mod r)

a codimension-one coincidence in a prime-order group, probability about `1/r`
per instance — far below the 0.06 margin of §4. The predecessor's argument here
(`T_i^2 = 1`, "two-step returns impossible") is **withdrawn**: a return needs
`T_i * T_j = 1`, not `T_i^2 = 1`, and the correct statement is the codimension-one
one above.

No quotient by the order-12 Frobenius class and no canonicalisation is used. That
quotient is what creates the structured two-cycle family of density `1/(12R)`
raised against quotient walks on this track; in this session the same walk with
that quotient spent 22.9% of all steps in such traps and cost 2.99 sqrt(N) per
solved instance against 1.36 sqrt(N) without it.

## 6. Cost ledger (exact)

Units: one dense `Fp12` multiplication at the organizer's fixed `C = 20,000`
instructions; `r = 2^254.857089`, `sqrt(r) = 2^127.428545`.

Per evaluation: `1 + 328/20000 = 2541/2500 = 2^0.023468`.

| item | instructions | note |
| --- | --- | --- |
| dense Fp12 multiplication | 20,000 | the reference operation |
| jump index from one 64-bit limb, xor-shift | 20 | |
| two 4-limb exponent additions mod r | 40 | |
| 12 live-word comparisons (load, cmp, jcc, loop) | 48 | 4 per live checkpoint |
| eviction test, loop, pointer and address overhead | 160 | |
| checkpoint creation, amortised over the grid | 60 | `k` insertions per `2^j` steps |
| **total per evaluation** | **20,328** | |

Prefix comparison discipline: the stored word is `X & (2^16-1)` of the canonical
representative, an *exact necessary condition* for equality (a bijective encoding
means equal elements have equal words), so the filter cannot miss a true match; it
can only produce false positives, which are resolved by a full 576-byte comparison
(120 instructions). Expected false positives over the whole run:
`2^127.7 evaluations x 12 x 2^-64 = 2^67.3`, costing `2^74.2` instructions
`= 2^59.9` units — charged in the additive line.

Additive work:

| item | units |
| --- | --- |
| jump table, 2048 GT scalar multiplications | 784,384 |
| one pairing for `g = e(P,Q0)` | 131,072 |
| verification `e([a]Q0,P) = z` and `[a]Q0` | 135,168 |
| expected prefix false positives | `2^60.0` |
| **total additive** | **1,152,921,504,607,897,600** `= 2^60.0` |

Total work

    W = B * 2541/2500 + additive = 2^127.666138,   claim = 12767/100 = 127.67

The additive term is `2^60.0` against a main term of `2^127.67`; it moves the
claim by less than `2^-67` bits and is carried anyway.

Peak memory: 656,704 B = 1024-entry jump table (589,824 B) + 12 checkpoints of
112 B (1,344 B) + 64 KiB of walk state and stack. Nothing grows with the walk
length, with the number of recorded points, or with the number of trials.

## 7. Heuristics disclosed

1. The random-map model for `X -> X * Th[H(X)]` on `GT`, with the Rayleigh
   collision time. Standard, supported by the toy run of §4 at production `R`
   (0.4567 against a predicted 0.45) and by the exactness of the detection index
   given `(mu, lam)`.
2. The uniform tail `mu` in the model. The toy run's median detection index
   (0.70 sqrt(q) over the solved trials) sits where the model puts it.
3. The per-item instruction counts of §6. They are stated as the derivation of
   each item, and the claim is insensitive to them: removing the 12-word
   comparison line entirely moves the claim from 127.666 to 127.663 (0.003 bits),
   and charging a full 576-byte comparison against all 12 states instead of one
   word moves it to 127.762 (0.096 bits).
4. The success probability is a model quantile, not a proof. That is why the cap
   is set at the 0.45 quantile while only the 0.39 floor is claimed.

## 8. Differences from this account's rejected rho candidates

| rejection | this candidate |
| --- | --- |
| `2*rho+3` evaluations is false (`d0f054cf`, `21c408d5`, `e5e6249f`) | no cycle-length bound is used anywhere; the index is the verified law of §3 |
| 56 claimed states against 192 actual (`624fef4`) | the count is enumerated in `experiment.py` §B1 at 1024, 4096 and 2^20 and the bound is `k*log2(W) = 12`; the old formula is printed beside it |
| "the aggregate evidence does not match the declared cap"; `R=64` not `R=1024` | `solved + capped = trials_run` is asserted, and the declared run uses `R = 1024` |
| "returning requires `T_i*T_j = 1`" | the wrong two-step argument is withdrawn; §5 gives the codimension-one statement |
| quotient walks trap on structured two-cycles | no quotient, no canonicalisation |

## 9. Reproduce

```
cd submissions/bls12-381-fixed-g1-invert-g2/gt-stagger-rho-v1
echo '{"production_checks":2,"table_checks":2,"trials":300,"toy_bits":26,
       "R":1024,"cap_factor":1.16,"seconds":12,"seed":1}' | python3 -I experiment.py
```

`results.json` is that exact command's stdout. Offline model tools (not shipped,
they only produce the tables of §4): an exact-law simulator, a deterministic
quadrature over the random-map model, a 2,000,000-draw Monte Carlo cross-check
and a 20-phase sweep.
