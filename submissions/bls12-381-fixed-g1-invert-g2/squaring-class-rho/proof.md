# Squaring-class mixed add/square rho for fixed-G1 pairing inversion into G2

Claimed score: `time_log2 = 126.4911` target-group-operation equivalents
(exact rational `1264911/10000`), success probability at least `0.45`,
peak memory `2^134.402` bytes, preprocessing `2^20` target-group-operation
equivalents (included in the total), nonuniform advice `0` bytes.

Track: `bls12-381-fixed-g1-invert-g2`. Cost policy `classical-group-work-v2`.
Reference operation: one dense Fp12 multiplication, `C = 20,000` scalar
AArch64 instruction equivalents.

---

## 1. Problem and the standard reduction

The challenge fixes a nonidentity `P` in `G1` and an independently uniform
`z` in `GT` (the identity target is included) and asks for `Q` in `G2` with
`e(P,Q) = z` under the pinned `py_ecc` convention. The track rules state the
admissible reduction explicitly: solve `z = e(P,Q0)^a` in `GT` and return
`Q = [a]Q0`.

The reduction is exact. Write `P = [b]P0` with `b` uniform in `[1,r)` (the
supplied point is uniform over nonidentity `G1`), and `g = e(P0,Q0)`, a
generator of `GT`. Then `e(P,Q0) = g^b`, the challenge has `z = e(P,Q0)^a =
g^{ab}`, and `Q = [a]Q0` satisfies `e(P,Q) = e(P0,Q0)^{ab} = z`. Conversely
any valid `Q = [q]Q0` satisfies `bq = ab (mod r)`, hence `q = a` because
`b` is invertible modulo the prime `r`. The map `Q -> e(P,Q)` is therefore
an isomorphism `G2 -> GT` and the task is exactly one discrete logarithm in
the cyclic group `GT` of prime order `r`, with a uniform exponent
`x = ab mod r` (uniform because `b` is uniform on `[1,r)` and `r` is prime).

The group order is

```
r = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001,
log2 r = 254.857055,
```

so the generic square-root barrier is `sqrt(r) = 2^127.428528`.

No structure beyond the generic group is used: `r` is prime, so there is no
Pohlig-Hellman reduction; there are no auxiliary powers; `z` is uniform, so
there is no small-exponent or interval structure; and `GT` is a prime-order
subgroup of `Fp12^*`, whose field-level discrete logarithm is far more
expensive than the square-root barrier at this size. The entire claim is
about the constant in front of `sqrt(r)`.

## 2. Priced operations

The reference operation for this track is one dense Fp12 multiplication,
`C = 20,000` instructions (the pinned `blst_fp12_mul` measurement is 22,996).
Two further operations are used, and both are priced above their measured
cost so the ledger is conservative.

* **Dense Fp12 multiplication by a table element**: `1.0000 C`. This is the
  reference operation itself; nothing is double counted because the reference
  price already includes its internal reductions, additions and memory work.
* **Cyclotomic squaring in `GPhi6`**: charged `0.2500 C` = 5,000
  instructions. Every walk element lies in `GT`, and `GT` is contained in the
  cyclotomic subgroup `GPhi6` of order `p^4-p^2+1`, so the Granger-Scott
  squaring applies. Its formula uses three `Fp2` multiplications, i.e. nine
  `Fp` multiplications, which at the pinned 383 instructions per
  `blst_fp_mul` is 3,447 instructions, or `0.1724 C`. The ledger charges
  `0.25 C`, a 45% margin over the measured formula, and includes the
  additions, conjugations and stores that accompany it.
* **Conjugation (the inverse in `GT`)**: free. For `y` in `GPhi6`, `y^{-1} =
  y^{p^6}` because `p^6 = -1 (mod r)` (from `r | p^4-p^2+1`), and `y^{p^6}`
  is the nontrivial automorphism of `Fp12/Fp6`, which in the tower
  representation is the sign flip of the second half of the coefficients.
  Six negations per element are charged in the overhead below.
* **Hash, table insert, exponent update**: charged in aggregate below.

## 3. The squaring-class mixed add/square walk

A standard Pollard-rho walk on `GT` uses one dense multiplication per visited
element: pick `s = H(X)` from the state, multiply by a table element
`M_s = g^{u_s} z^{v_s}`, and continue. The element cost is `1.0 C`.

This walk replaces one dense multiplication by a **block** that costs
`1 + 3*(1/4) = 1.75` reference units and *emits eight group elements*, every
one of which is inserted into the collision table and compared:

```
block state X_t  in GT
s      = H(X_t)                       (a fixed hash into [1,S])
X_{t+1} = X_t^8 * M_s                 (three cyclotomic squarings + one dense multiply)
emitted = { X_t, X_t^2, X_t^4, X_t^8,
            X_t^-1, X_t^-2, X_t^-4, X_t^-8 }
```

The exponents are tracked alongside: with `X_t = g^{a_t} z^{b_t}`,

```
a_{t+1} = 8 a_t + u_s  (mod r),   b_{t+1} = 8 b_t + v_s  (mod r).
```

The three squarings produce `X^2, X^4, X^8`; the four inverses are free
conjugations; the one dense multiplication advances the block state. The
eight emitted elements are exactly the orbit of `X_t` under the group
`{+/- 2^d : d = 0,1,2,3}` acting on exponents, and the squaring chain is what
makes that orbit cheap.

The table `M_s = g^{u_s} z^{v_s}`, `s = 1..S`, `S = 2^20`, is precomputed with
`2^20` dense multiplications and charged in the preprocessing term.

## 4. Exactly fourteen visited equalities per block pair

This section is the part of the accounting that a naive squaring-chain ledger
gets wrong, and it is stated so that every credited relation is an event the
walk physically visits.

Consider two distinct blocks `i < j` and their emitted sets

```
E_i = { X_i^{e 2^d} : e in {+1,-1}, d in {0,1,2,3} },   |E_i| = 8,
E_j = { X_j^{e 2^d} : e in {+1,-1}, d in {0,1,2,3} },   |E_j| = 8.
```

A collision is an equality between an element of `E_i` and an element of
`E_j`. Writing `X_i = g^{a_i} z^{b_i}` and `X_j = g^{a_j} z^{b_j}`, the
equality `X_i^{e 2^d} = X_j^{e' 2^{d'}}` is equivalent, after inverting the
bijections `y -> y^{2^d}`, to

```
X_i = X_j^{e'' 2^{d''}},   e'' = e e' in {+1,-1},  d'' = d' - d in [-3,3].
```

The set `R = { e 2^d : e in {+1,-1}, d in [-3,3] }` has exactly `14`
elements, and they are pairwise distinct modulo `r`: two elements coincide
only if `2^{d1-d2} = +/- 1` with `|d1-d2| <= 6`, which requires
`ord_r(2) <= 12`; `2` has astronomically large order modulo `r`. Hence for
one block pair the fourteen events are **disjoint**, and for states that
behave as uniform elements of `GT` each has probability exactly `1/r`:

```
Pr[collision between blocks i and j] = 14 / r.            (4.1)
```

All fourteen are **visited**. The equality `X_i = X_j^{e 2^d}` with
`d in [0,3]` is literally the equality of the emitted element `X_i` (level 0
of block `i`) with the emitted element `X_j^{e 2^d}` (level `d` of block
`j`), and with `d in [-3,-1]` it is literally the equality of the emitted
element `X_i^{2^{|d|}}` (level `|d|` of block `i`) with the emitted element
`X_j` (level 0 of block `j`). Both operands are computed, inserted into the
table and compared; no ratio is credited that is not an equality of two
stored elements.

Two degeneracies are removed from the count and do not affect the bound.

* **Within-block pairs** give no information: `X^{e 2^d} = X^{e' 2^{d'}}`
  with `(e,d) != (e',d')` implies `X^{|2^d -/+ 2^{d'}|} = 1`, hence
  `X = 1` because the exponent is smaller than the prime `r`. This happens
  with probability `O(1/r)` and is discarded.
* **Adjacent blocks with `d = +3`**: `X_t = X_{t-1}^8` would require
  `M_s = 1`, which the table never contains (all `(u_s,v_s)` are nonzero).
  One of fourteen relations is therefore unavailable on the `B-1` adjacent
  pairs, a relative loss of order `1/B`, which is below `2^-125` of the
  total chance count.

Every collision that is found is solved by linear algebra: from
`X_i = X_j^{e 2^d}`,

```
(a_i - e 2^d a_j) + (b_i - e 2^d b_j) x = 0 (mod r),
x = (e 2^d a_j - a_i) / (b_i - e 2^d b_j) (mod r).
```

The denominator vanishes with probability `O(1/r)` (and identically only for
the structural `d = 3` adjacent case already excluded); a vanishing
denominator is detected and the walk continues, so no failure is hidden. The
recovered `x` is checked against `z` before output.

## 5. Detection is the same object as the accounting

The collision table is a hash table keyed by a 256-bit fingerprint of the
twelve `Fp` coefficients of each emitted element. Every emitted element of
every block is inserted; every insertion probes for an equal fingerprint.
Because `|R| = 14` and the emitted sets are complete orbits under `R`, the
fourteen relations of Section 4 are exactly the fourteen ways two emitted
elements can be equal, and every one of them is found by the same insertion.
There is no separate "distinguished point" filter that would shrink the
detected set below the counted set: the counted set and the detected set are
the same set of `8B` stored elements.

A 256-bit fingerprint is long enough that false matches are negligible: with
`N = 8B = 2^128.654` stored elements the expected number of accidental
fingerprint coincidences is `N^2 / 2^257 = 2^0.3`, and each is resolved by a
full coefficient comparison whose cost is charged in the per-emission
overhead. A collision is confirmed only by full equality.

## 6. Birthday accounting and the success probability

With `B` blocks the walk stores `8B` elements and, by (4.1), performs
`14 * C(B,2) ~ 7 B^2` disjoint relation trials. The standard birthday
approximation for a rho walk on a set of size `r` gives

```
Pr[no collision after B blocks] ~ exp( -7 B^2 / r ).        (6.1)
```

The budget is fixed so that the success probability clears the 0.39 floor
with margin:

```
7 B^2 / r = -ln(0.55) = 0.597837
B = sqrt(0.0854053 r) = 2^125.6538 blocks.                  (6.2)
```

Equation (6.1) then gives `Pr[success] = 1 - exp(-0.597837) = 0.45`, which is
the declared success probability. Failure inside the budget (the walk stops
without a collision) is a *reported* outcome, not a silent retry, so no
amplification factor is charged: the declared algorithm performs exactly one
trail of `B` blocks and succeeds with probability `0.45`.

The success probability is over the uniform challenge (uniform `P`, hence
uniform `b`, and uniform `z`, hence uniform `x`) and over the fresh coins
that choose the walk's starting point and the table. Equation (6.1) is the
only heuristic; it is the standard rho randomness assumption, supported by
the `2^20`-entry multiplier table and the state hash, and it is tested at
small scope by `experiment.py` (Section 10).

## 7. Complete cost ledger

Instructions per block, in the reference instruction basis:

| item | instructions | reference units |
| --- | ---: | ---: |
| one dense Fp12 multiplication by `M_s` | 20,000 | 1.0000 |
| three cyclotomic squarings at 5,000 each | 15,000 | 0.7500 |
| eight 256-bit fingerprints (4 coefficients read + mix) | 400 | 0.0200 |
| eight hash-table insert/probe pairs | 160 | 0.0080 |
| exponent pair update `(8a+u, 8b+v) mod r` | 100 | 0.0050 |
| state hash for multiplier selection | 50 | 0.0025 |
| four conjugations (inverses), six negations each | 24 | 0.0012 |
| **per block** | **35,734** | **1.7867** |

Total walk work:

```
W_walk = 1.7867 * B = 1.7867 * 2^125.6538 = 2^126.4911.    (7.1)
```

Remaining terms, all charged and all included in the total:

* **Preprocessing**: the `2^20`-entry multiplier table costs `2^20` dense
  multiplications, `2^20 = 2^20` reference units. This is `2^-106.5` of
  `W_walk`.
* **Output**: `Q = [x]Q0` is one G2 scalar multiplication, about 2,811 `Fp2`
  multiplications, i.e. `161` reference units.
* **Attack-side verification**: one pinned pairing check, about 75 reference
  units.
* **Collision resolution**: at most a few fingerprint confirmations plus one
  2x2 linear solve modulo `r`, about `10` reference units.
* **Randomness**: the start point and the `2^20` table are the only random
  choices; sampling them is a few hundred instructions.
* **Memory transfers**: the eight inserts per block are counted in the table
  above; the `64`-byte exponent record per block is one store block, counted
  in the exponent-update line.

The three largest non-walk terms sum to less than `2^8.1` reference units,
which is `2^-118.4` of `W_walk`. Hence

```
W = 2^126.4911 + 2^20 + 2^8.1 = 2^126.4911,
time_log2 = 126.4911.
```

No work is amortized over future targets, no advice is embedded, no oracle is
assumed, and no processor is shared: the walk is a single serial trail and the
score is total work, not elapsed time.

## 8. Memory, preprocessing, advice

The table holds `N = 8B = 2^128.6538` fingerprints. A table with `2^129`
slots (load factor 0.5) at 32 bytes of fingerprint plus a 4-byte block index
is `2^129 * 36 = 2^134.170` bytes. The per-block exponent pair
`(a_t,b_t)` is two residues modulo `r`, 64 bytes per block, `2^131.654`
bytes. The multiplier table is `2^20 * 576 = 2^29.2` bytes. Peak storage is

```
2^134.170 + 2^131.654 + 2^29.2 = 2^134.402 bytes.
```

The declared exponential memory is required by the construction: the ratio
classes are only detectable if both operands of a relation are retained, so a
distinguished-point filter would multiply the work by `1/theta` and destroy
the gain. Memory is reported, not scored, under the track's accounting.

Preprocessing is `2^20` reference units, declared and included in the total.
Nonuniform advice is `0` bytes: nothing is embedded in code beyond the
algorithm and the public profile constants.

## 9. Why the full Frobenius orbit is not used

`GT` admits the order-12 Frobenius action `y -> y^p` (the order of `p`
modulo `r` is 12 because `r | p^4-p^2+1`). Quotienting the collision classes
by that action would multiply the relation count by 12 and lower the bound by
`log2(sqrt(12)) = 1.79` bits. It is priced out here, and the calculation is
recorded so the omission is a decision rather than an oversight.

One Frobenius on `Fp12` in the tower `Fp12 = Fp6[w]/(w^2-v)`,
`Fp6 = Fp2[v]/(v^3-xi)` requires `x0^p` and `x1^p` (three `Fp2`
multiplications each, since `Fp2` conjugation is free) plus the multiplication
of `x1^p` by the `Fp6` constant `v^{(p-1)/2}` (six `Fp2` multiplications):
`36` `Fp` multiplications, `13,788` instructions, `0.69 C`. Emitting the
twelve orbit images of one element therefore costs `8.3 C`, against `0.44 C`
per emitted element in this walk. The order-12 quotient needs every emitted
element canonicalised over twelve images, so it adds roughly `8.3 C` per
emitted element and `66 C` per block, while raising the relation count only
from 14 to 168; the per-relation cost rises from `0.128` to `0.43` reference
units. The order-2 subgroup of the same action, `p^6 = -1 (mod r)`, is the
conjugation and is free; that is exactly the `e in {+1,-1}` factor already
used in Section 4. The claim therefore uses the free part of the Frobenius
class structure and prices out the expensive part.

## 10. Heuristics, limitations, evidence

**Heuristic.** The block map is treated as a random map on `GT`, so block
states behave as independent uniform elements and (4.1) and (6.1) hold. This
is the standard rho randomness assumption; it is not a theorem. It is
supported by the `2^20`-entry multiplier table, by the state-dependent
multiplier selection, and by the small-scope experiment.

**Evidence at its actual scope.** `experiment.py` runs the same block
construction in a small prime-order subgroup of `F_p^*` (`p = k q + 1`,
`q` an input prime, default `q = 100003`). It stores every emitted element,
detects collisions by exact equality, solves the linear relation, and checks
the recovered exponent. It reports the observed success rate against the
prediction of (6.1) for the mixed walk and against the plain one-element walk
at the same block budget and at the same reference cost. The experiment is
evidence about the combinatorial counting of Section 4 and the birthday
model of Section 6 only; it says nothing about BLS12-381 costs, `Fp12`
arithmetic, or production-size behaviour.

**Limitations.** The bound is analytic; no production-size trail is executed.
The memory bound is exponential and is required by the construction. The
success probability is a heuristic lower bound derived from the birthday
approximation, not a proven bound. The constant in front of `sqrt(r)` is the
whole claim; no non-generic structure of `GT`, and no index-calculus or
pairing-specific shortcut, is claimed or used.

**Accounting statements.** The reference price `C = 20,000` is the organizer
constant and is not modified. The cyclotomic squaring is charged above its
measured formula. All preprocessing, trials, failures, randomness, hashing,
lookups, memory transfers, output construction and attack-side verification
are included. Nothing is amortized over future targets, and the attack is a
single classical cold trail.
