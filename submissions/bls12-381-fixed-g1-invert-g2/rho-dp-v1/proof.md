# bls12-381-fixed-g1-invert-g2: Pollard-rho discrete logarithm in GT

Status: ready for independent review. The first declared claim in
`submission.json` supplies the score `127798437/1000000 = 127.798437`.

This candidate is a self-contained generic attack. It makes no claim of
novelty, no claim of a practical break, and no claim that index-calculus
methods are unavailable. It is an admissible attack under the track's stated
approach: solve `z = e(P,Q0)^a` in GT and return `Q = [a]Q0`. All production
parameters are used; nothing is scaled down.

## Target and convention

The protected profile fixes BLS12-381, the subgroup order

```text
r = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
```

(255 bits, prime, `log2 r = 254.8570894130472`), the G1 generator P0, the G2
generator Q0, and the pairing

```text
e(P,Q) = py_ecc.optimized_bls12_381.pairing(Q, P, final_exponentiate=True)
```

at the pinned commit, with the library taking G2 first. The profile hash is
`af94a59f...abef888` and the cost-policy hash is `e0edc25a...d5e17294`;
`claim.json` binds both.

GT is the group of nonzero z in Fp12 with `z^r = 1`; it is cyclic of prime
order r. The identity is `[1,0,...,0]`. The challenge independently samples
P = [b]P0 for uniform b in `1..r-1` and z = e(P,Q0)^a for uniform a in
`0..r-1`, and releases only P and z.

## Reduction used by the attack

Given the supplied nonidentity P, compute

```text
g = e(P,Q0) = pairing(Q0, P, final_exponentiate=True).
```

Bilinearity gives `g = e(P0,Q0)^b`. Since P0 and Q0 have prime order r and the
pairing is nondegenerate, `e(P0,Q0)` generates GT, so `g != 1` and g generates
GT. The released target is `z = e(P,Q0)^a = g^a`. Therefore the unique
`a = log_g z` in `[0,r)` is the challenge exponent, and

```text
Q = [a]Q0  satisfies  e(P,Q) = e(P,Q0)^a = g^a = z.
```

For a uniform b in `1..r-1`, `g = e(P0,Q0)^b` is uniform over GT minus the
identity. Since z is independently uniform over GT (including the identity),
`a` is uniform over `[0,r)` and is independent of g. The case z = 1 is exactly
a = 0; the attack returns the G2 identity, which the profile accepts exactly
for z = 1. No auxiliary powers are used and P is not replaced.

Thus the production problem reduces, for every valid instance, to one discrete
logarithm in the prime-order group GT, followed by one G2 scalar
multiplication and one verification pairing.

## Algorithm

All group elements below are in GT. Write `X = g^u z^v` for known exponents
`u,v` in `Z_r`; initially the walk starts at a randomized point whose
exponents are known relative to g and z.

### Tables

Choose 128 small multipliers `c_j = j+1` for `j = 0..127`, and precompute
`h_j = g^{c_j}` by repeated multiplication (127 dense Fp12 multiplications;
`h_0 = g`). Also precompute `g^{2^i}` for `i = 0..254` (255 multiplications)
for the randomized walk starts. All tables total 383 dense Fp12
multiplications.

### Partition

Represent X canonically as its 12 base-field coefficients in `[0,p)`. Let
`w_i` be the low 64-bit word of coefficient i, and define

```text
chi(X) = (w_0 xor w_1 xor w_2) mod 256.
```

This is a deterministic function of the group element (equal elements have
equal canonical encodings). Split on `t = chi(X) mod 3`:

```text
t = 0:  X <- X * h_{chi(X)>>1}       (chi>>1 in 0..127),  u <- u + c_{chi>>1} mod r
t = 1:  X <- X * z,                                      v <- v + 1 mod r
t = 2:  X <- X^2,                                        u <- 2u mod r, v <- 2v mod r
```

Every iteration is one dense Fp12 multiplication or one Fp12 squaring. Both
are charged as one reference operation; the squaring is not discounted.

### Distinguished points and collision detection

Fix the distinguished-point density at `2^-80`. X is distinguished when the
low 64 bits of `w_0` and the low 16 bits of `w_1` are both zero. This predicate
is a deterministic function of X.

For each walk:

1. Sample s uniformly in `[0,r)` from a cryptographic PRNG (rejection
   sampling); set `(u,v) = (s,1)` and `X = z * g^s`, using the power table.
2. Iterate the partition step until X is distinguished, with a hard abort at
   `128 * 2^80 = 2^87` iterations (probability at most `e^-128` per walk).
3. Insert `(X,u,v)` into a hash table keyed by the canonical encoding of X.
   A key collision means two walks reached the same distinguished element, so
   the two trajectories intersect. With entries `(X,u1,v1)` and `(X,u2,v2)`,
   the collision gives

   ```text
   u1 + v1*a = u2 + v2*a (mod r),   so   a = (u1-u2)/(v2-v1) mod r
   ```

   provided `v1 != v2`. If `v1 = v2`, discard and continue (negligibly rare;
   see limitations).
4. Verify `g^a = z` (255 dense multiplications), set `Q = [a]Q0` (G2 scalar
   multiplication), verify `pairing(Q,P) = z`, and return Q. For z = 1 return
   the identity before starting a walk.

Every group state is stored exactly once per encountered distinguished point;
all table accesses, hashing, probing, record stores, PRNG calls, scalar
bookkeeping, the final scalar multiplication, the verification pairing and the
output are charged in the ledger below.

## Success argument

The walk is deterministic on the group element and mixes through the 3-way
adding/squaring step. Under the standard random-mapping heuristic, after total
work x the probability of at least one collision among visited states is
`1 - exp(-x^2/(2r)) + o(1)`; the expected first-intersection work is

```text
E[T] = sqrt(pi*r/2) = 1.2533141373 * sqrt(r) = 2^127.7542927713.
```

The attack budget is

```text
t = ceil(63*sqrt(r)/50) = 1.26*sqrt(r) = 2^127.7619684402,
```

which exceeds `E[T]` by a factor 1.0053345 (about `2^120.2` spare
iterations). Distinguished-point detection happens when the later colliding
walk reaches its endpoint; the extra uncharged-beyond-t work is at most the
walk abort bound `2^87` iterations, which is smaller than the spare budget
`2^120.2`. Hence the expected total number of walk iterations, including all
restarts and detection delay, is at most t. The walk succeeds with
probability 1 in the Las Vegas sense: on the negligible abort or degenerate
events it restarts with fresh coins, and the expected restart work is below
one iteration times a factor `1 + 2^-130`. Therefore the claimed success
probability over uniform P, z and fresh coins is 1 (at least 0.39), and all
trials, failures and recovery are already inside the expected total.

A collision yields the correct a unless `v1 = v2`. For a collision of an
effectively random mapping the two v-coefficients are effectively independent
uniform elements of Z_r, so this happens with probability about `1/r`; the
algorithm restarts in that case. The abort probability per walk is at most
`e^-128`, and there are at most `2^48` walks, so the total failure probability
before a restart is below `2^-130`. This is a heuristic statement about the
walk, disclosed below, not a theorem for this specific finite-field map.

## Cost ledger

The cost model fixes C = 20,000 scalar AArch64 instruction equivalents per
dense Fp12 multiplication, and scores `log2(H + R/C)` with all work charged
exactly once. The reference multiplication price includes its internal
arithmetic, loads, stores and reduction; the 512 instructions below are only
the surrounding walk bookkeeping, so nothing is double charged.

### Reference-operation work

| Component | Amount (target-group-operation equivalents) |
| --- | ---: |
| Walk iterations `t` (one dense Fp12 mul or squaring each) | 288,525,900,792,338,468,361,715,506,970,096,608,215 |
| Ancillary group work reserve (detection delay, starts, preprocessing, checks) | 309,485,009,821,345,068,724,781,056 = 2^88 |

The second row covers everything that is not a main walk iteration:

- detection-delay iterations at most `2^87`;
- walk starts: at most `2^48` walks times 128 dense multiplications for
  `g^s`, below `2^55`;
- preprocessing table `h_j` (383 dense multiplications) and the two pairings
  `e(P,Q0)` and `e(P,Q)`; a full pairing is charged at most `2^20` dense
  multiplication equivalents each, below `2^21` in total;
- the `g^a = z` check (255 dense multiplications) and the G2 scalar
  multiplication `[a]Q0` (below `2^22` scalar instruction equivalents,
  dominated by 2^88).

All of these are below `2^88` in total, so the reserve is a valid upper bound.
Preprocessing and advice are included in total work; there is no reuse credit.

### Scalar-instruction work

Each walk iteration performs at most 512 scalar AArch64 instructions outside
the reference kernel: the `chi` extraction (loads, xor, mask), the 3-way
branch, one or two 256-bit modular additions/doublings of u and v (with
conditional subtraction), the 80-bit distinguished-point test, pointer swaps
and loop control. A conservative breakdown is below 230 instructions; 512 is a
factor-2.2 margin. Thus

```text
R_main = 512 * t = 147,725,261,205,677,295,801,198,339,568,689,463,406,080
```

A reserve `R_other = 2^97 = 158,456,325,028,528,675,187,087,900,672`
instructions covers all non-reference work not proportional to t: 512
instructions on each of at most `2^88` ancillary iterations, hash-table
insertion and probing for at most `2^48` records (each at most `2^16`
instructions, below `2^64` in total), record transfers of `2^57.5` bytes,
PRNG output, collision resolution, the verification pairing, scalar
multiplication and serialization.

### Total and score

```text
W = H + R/C
  = 288,525,900,792,338,468,361,715,506,970,096,608,215
    + 309,485,009,821,345,068,724,781,056
    + (147,725,261,205,677,295,801,198,339,568,689,463,406,080
       + 158,456,325,028,528,675,187,087,900,672) / 20,000
  = 295,912,163,852,939,740,977,848,195,451,015,216,836.

log2(W) = 127.79843660731963, rounded upward to 127.798437.
```

The two reserves are `2^-39.76` of the main term and `< 2^-13` of it in
combination with the per-step overhead; they cannot hide a material cost.

### Memory, preprocessing, advice

The distinguished-point table holds at most `ceil(t/2^80)+2 = 238,663,031,354,822`
records. With 768 bytes per record (576 bytes for the canonical Fp12 element,
32 bytes each for u and v, and 128 bytes of metadata), an 11/10 open-addressing
load factor, the 383-entry multiplier tables and a 1 MiB code/scratch reserve,
peak memory is

```text
262,529,334,490,306 slots * 768 bytes + 1 MiB + 383*576 bytes
= 201,622,528,889,824,768 bytes = 2^57.4844344647.
```

Reported memory_log2_bytes is 57.484435. Preprocessing is at most `2^22`
target-group-operation equivalents (383 table multiplications plus two
pairings charged at `2^20` each) and is included in W; preprocessing_log2 is
22.000000. Nonuniform advice is zero; the only precomputed values are the
publicly computable powers of the instance base g.

### Arithmetic checks

The exact integers above were recomputed with Python's arbitrary-precision
integers from the profile's r. The rounding is upward at six decimal places,
as in the organizer control. The score uses the first claim only; no other
claim is declared.

## Experiments and local verification

`experiment.py` is a standard-library toy that runs the exact exponent
bookkeeping and distinguished-point collision logic in an additive group of
prime order (not a production curve and not a cost proof). It recovers random
discrete logarithms and checks the recovered exponent against the target. The
repository's `tests/test_bls_challenges.py` additionally checks the pinned
pairing convention and the fixed-G1 challenge semantics with py_ecc; the
candidate was reviewed against those checks. No full-size execution is
claimed: 2^127.76 walk steps are analytic, and the experiment is evidence only
at its stated toy scope.

## Heuristics, limitations and objections

- Random-mapping walk: the collision-time bound `sqrt(pi*r/2)` is a standard
  Pollard-rho / van Oorschot-Wiener heuristic, not a proof for this specific
  partition. The 3-way adding/squaring walk with a balanced 256-way partition
  is the standard construction; the 0.53% budget margin covers ordinary
  constant-factor variation.
- Distinguished-point uniformity: the predicate uses 80 low bits of canonical
  Fp12 encodings, assumed close to uniform for random GT elements. The
  expected number of walks and records depends on it; the memory and
  insertion reserves are generous.
- Eventual termination: a walk that enters a distinguished-point-free cycle is
  aborted at 128 times the mean walk length; the per-walk abort probability is
  at most `e^-128` and the total is below `2^-130`.
- Degenerate collision: if the two colliding trajectories have equal v
  coefficients the equation gives no a; the probability is about `1/r` per
  collision and the algorithm restarts. This is part of the disclosed
  random-mapping heuristic.
- Expected rather than worst-case cost: the attack is Las Vegas; it succeeds
  with probability 1 and the claimed total work is its expected work,
  including every walk, restart and detection delay. This matches the cost
  model's success, trial and recovery requirements; no uncharged amplification
  remains. A fixed-budget one-shot variant at the same t would succeed with
  probability 0.548, and charging recovery at that probability would only
  increase the cost, so the expectation bound is the tightest honest claim.
- Memory: 2^57.5 bytes is astronomically large. The track scores computation
  only, memory has no scalar contribution, and the organizer control itself
  declares 2^137.6 bytes. The value is analytically justified, not practical.
- Index calculus: finite-field NFS/TNFS in Fp12 is a different and potentially
  faster attack on the same GT discrete logarithm. This submission does not
  claim to be the best-known attack, does not use those methods, and charges
  no work for them. Its bound is valid for the algorithm described here.
- The walk uses canonical GT encodings only for the partition and
  distinguished-point predicate; an implementation is free to use any correct
  Fp12 arithmetic. The reference multiplication price is the track's fixed
  normalization, not a measured property of this candidate.

## References

- J. M. Pollard, "Monte Carlo methods for index computation (mod p)",
  Mathematics of Computation 32 (1978).
- P. C. van Oorschot and M. J. Wiener, "Parallel collision search with
  cryptanalytic applications", Journal of Cryptology 12 (1999), for the
  `sqrt(pi*N/2)` collision-search expectation and distinguished points.
- The protected `tracks/docs/COST_MODEL.md` and `target-profile.json` for the
  reference price, pairing convention and acceptance rules.
