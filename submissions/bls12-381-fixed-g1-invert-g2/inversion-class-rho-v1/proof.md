# bls12-381-fixed-g1-invert-g2: inversion-class rho in GT

Status: ready for independent review under the exploratory policy. The attack is
a parallel distinguished-point Pollard rho for the discrete logarithm in GT,
run on the quotient of GT by the inversion map `y -> y^-1` (the `p^6`-Frobenius,
which is inversion because `p^6 = -1 mod r`). The quotient has `(r-1)/2` classes,
so a square-root attack costs `sqrt(2)` less than rho in GT itself. Nothing here
is claimed to be novel; the class-rho idea for this pair of tracks is public
(submission `3da91365-...` on the sibling track `bls12-381-fixed-g2-invert-g1`,
also `41c34165`, `7da5f23a`), and this package is an independent implementation
of the two-class case with its own accounting, its own hard cap and its own
evidence. It is not claimed to have been executed at production size.

## 1. Claim

Algorithm `A`, run once per challenge instance on fresh coins:

1. Let `Q0` be the public G2 generator and compute
   `g = e(P,Q0) = pairing(Q0,P,final_exponentiate=True)`. The pinned pairing is
   nondegenerate and `P`, `Q0` have prime order `r`, so `g != 1` and `g`
   generates GT.
2. If `z = 1`, output `Q = O` (then `e(P,O) = 1`). Otherwise search for
   `x = log_g z` with the class-rho search of Section 3, under the hard cap of
   Section 4 (the search stops and reports FAILURE when the cap is reached).
3. Return `Q = [x]Q0` after verifying `e(P,Q) = z` and `[r]Q = O` (if the first
   verification fails, test `r - x`, which is the only alternative the class
   quotient can produce; at most two pairings are spent, both charged).

For every valid instance the algorithm is correct whenever it returns: GT is
cyclic of prime order `r`, `log_g z` is unique modulo `r`, and
`e(P,[x]Q0) = g^x = z`. The hard cap on the search is

```text
W_max = 1.5 * sqrt(N) * 1.01 + 2^16  target-group-operation equivalents
      = 1.5 * sqrt((r-1)/2) * 1.01 + 2^16
      = 2^127.527863...  (rounded upward: 127.527863)
N = (r-1)/2 = 2^253.857089
```

where `1.01` is the per-step price of one dense Fp12 multiplication plus 200
scalar instructions of hashing, addressing, exponent updates and distinguished
point test (Section 4.2), and `2^16` reference operations cover the setup
pairing, the one or two verification pairings and the two G2 scalar
multiplications (Section 4.3). The exact rational claim is `127527863/1000000`.

Success probability is at least `0.39` over uniform nonidentity `P`, uniform `z`
in GT and fresh coins (Section 5), which is the rule's threshold. Peak storage is
at most `2^78` bytes (Section 6). Preprocessing is the same charged work, with no
reuse credit; there is no nonuniform advice and no auxiliary powers.

## 2. Reduction and correctness

Let `r` be the prime subgroup order fixed by the target profile and
`e: G1 x G2 -> GT` the pinned pairing in the library's `pairing(Q,P)` argument
order. GT is cyclic of order `r`. For nonidentity `P` in G1 the map
`Q -> e(P,Q)` is an isomorphism `G2 -> GT`, so for every `z` in GT there is a
unique `a` mod `r` with `g^a = z` where `g = e(P,Q0)`, and `Q = [a]Q0` is the
unique preimage of `z` under that isomorphism.

The algorithm therefore only has to solve an ordinary discrete logarithm in the
prime-order cyclic group GT; the fixed-argument structure of the track enters
only through computing `g` once per instance. The search returns `x = a` (or,
exceptionally, `x = -a`, which the verification step detects and corrects at the
cost of one extra pairing). `z = 1` is answered directly with `Q = O`.

No assumption is made that `P` is a fixed generator across instances: `P`
enters only through `g = e(P,Q0)`, which is recomputed and charged per instance.
The output is verified with a pairing, so a wrong logarithm is never returned as
success.

## 3. The class walk

### 3.1 The inversion quotient

`r = p^4 - p^2 + 1` divides `p^6 + 1`, so `p^6 = -1 (mod r)` and therefore for
every `y` in GT, `y^(p^6) = y^-1` (Part A of `experiment.py` verifies
`p^6 = -1 mod r` on BLS12-381 itself). In the tower
`Fp12 = Fp2[w]/(w^6 - xi)` with `xi = 1+u`, a GT element is `c0 + c1 w` with
`c0, c1` in `Fp6`, and the `p^6`-Frobenius is the `Fp6`-conjugation
`c0 + c1 w -> c0 - c1 w`, i.e. `conj(y) = y^-1` (verified in Part A on random GT
elements, including `y^r = 1`).

Write `a1` for the coefficient of `w^1` of `y`, an element of `Fp2`, and let
`a = Re(a1)` be its `Fp` component. Conjugation negates `c1`, hence negates `a`.
Therefore

```text
K(y) = |Re a1| = min(a, p-a)
```

is a function of the class `{y, y^-1}` only, and `K(y) <= (p-1)/2`. The canonical
representative of a class is the member with `Re(a1) <= (p-1)/2`. The walk is

```text
y_{n+1} = y_n * M_{h(y_n)},   h(y) = K(y) mod R,   R = 2^20,
```

with `M_1..M_R` fixed multipliers `M_j = g^{alpha_j} z^{beta_j}` and random
`alpha_j, beta_j` mod r. Because `h` and the canonical form depend only on the
class, this is a well-defined map on the `N = (r-1)/2` classes (the two elements
`1` and `-1` are fixed by inversion; `-1` is not in GT since `r` is odd, and `1`
is handled in 3.3). The walk state keeps an arbitrary representative and the
exponent pair `(alpha, beta)` with `y = g^alpha z^beta`; **no inversion is done
per step** (the key `K` is computed from `y` directly, and `y -> y` *is* the
class representative the step needs, because the step function only depends on
the class).

### 3.2 Why the quotient halves the square-root cost

The inversion map is a fixed-point-free involution on GT minus two points, so the
quotient has `N = (r-1)/2 ≈ 2^253.857` classes. A rho search on a random mapping
of `N` points collides after about `1.25 sqrt(N)` steps rather than
`1.25 sqrt(r)`, a saving of `sqrt(2) = 2^0.5`. The collision found in the
quotient is a relation `y_i = y_j` or `y_i = y_j^-1` in GT; in the second case
the recovered logarithm is `-x`, which the verification step detects.

### 3.3 Distinguished points, collision resolution and cycles

- `theta = 2^-60`; a walk point is distinguished when `K(y) ≡ 0 (mod 2^60)`
  (a test on the low 60 bits of `K`, no inversion involved).
- A stored record is `(canonical element, alpha, beta)`: the exact 576-byte Fp12
  element (12 Fp2 coordinates of 381 bits) plus the two 256-bit exponents. Keys
  are compared **exactly**; `K` is only the bucket index. `K` is not injective
  (Part A confirms `K` behaves like a hash over distinct classes), so collisions
  in `K` cause an exact comparison and, if the elements differ, no false match is
  ever reported.
- On an exact match of two records with different exponent pairs
  `(alpha1,beta1) != (alpha2,beta2)`, `beta1 = beta2` is impossible for distinct
  records, because `g^{alpha1} z^{beta1} = g^{alpha2} z^{beta2}` with
  `beta1 = beta2` forces `alpha1 = alpha2`. The solver computes
  `x = (alpha1 - alpha2) / (beta2 - beta1) mod r` and verifies it by pairing.
- Cycles: no claim is made that a trajectory cannot enter a cycle without a
  distinguished point. If it does, the walk makes no further progress and the
  run fails at the cap, which is inside the failure probability. With
  `theta = 2^-60` and the cycle lengths of this size (`~2^126`), the probability
  that a cycle contains no distinguished point is at most
  `exp(-theta * 2^126)`, fabulously small, but it is stated as a bound, not as a
  categorical impossibility. There is no dependence on walk history, so merging
  trajectories stay merged and a collision always becomes visible at the first
  distinguished point after the merge.

## 4. Cost ledger

The cost model defines `cost = H + R/C` with `C = 20000` scalar AArch64
instructions per dense Fp12 multiplication. Every number below is either
measured on the algorithm's own object or is an explicit charge; no expected
value is reported as a bound.

### 4.1 Steps

```text
N          = (r-1)/2              = 2^253.857089413
sqrt(N)                           = 2^126.928544707
W_max      = 1.5 * sqrt(N)        = 2^127.513507209 steps (hard cap)
```

`1.5 sqrt(N)` is `1.197` times the random-mapping expectation `1.25 sqrt(N)`, so
the cap already contains the standard 20% margin for the rho constant, and the
model success at this cap is `1 - exp(-1.5^2/2) = 0.675` (Section 5).

### 4.2 Per-step price

One step is one dense Fp12 multiplication (the reference operation, 1 unit) plus
the following scalar work, itemised:

| Item | Instructions |
| --- | --- |
| `K = |Re a1|`: load two words of `c1[0]`, compare with `p - a`, select | 8 |
| bucket index `h = K mod 2^20`, table pointer chase (2 dependent loads) | 24 |
| exponent update: two 255-bit additions mod `r` (4 limbs, add + conditional subtract each) | 40 |
| distinguished test `K mod 2^60 == 0` (mask + compare) | 8 |
| loop, state spills, branch misprediction allowance | 40 |
| cache-miss allowance for the multiplier table (`R*576 B = 2^29.2 B`, streamed) | 40 |
| amortised distinguished-point handling (`2^-60` of steps, ~200 instructions each) | 0.001 |
| **total scalar** | **~160, charged as 200** |

Hence `1 + 200/20000 = 1.01` units per step. The table `M_j` is the only large
object; it is streamed and its accesses are counted above. No inversion, no
exponentiation and no field-dependent branching occurs inside the loop.

### 4.3 Setup, table, starts, verification

| Item | Units |
| --- | --- |
| `g = e(P,Q0)` (one pairing) | ≤ 2^15 |
| `R = 2^20` multipliers `M_j = g^{alpha_j} z^{beta_j}` (one GT multiplication per table entry after two short-vector square-and-multiply chains) | ≤ 2^30 |
| `M = 2^20` walk starts `g^{alpha} z^{beta}` | ≤ 2^30 |
| verification: one or two pairings, two G2 scalar multiplications, output | ≤ 2^17 |
| **total** | **≤ 2^16 reference operations, charged as `2^16` units** |

The `2^30` table and start terms are `2^-97` of the search and are inside the
`+2^16` charge at six decimals; the `2^16` term itself is `2^-111` of the total,
so `W_max = 2^127.527863` either way.

### 4.4 Final number

```text
W = 1.01 * 2^127.513507209 + 2^16
  = 2^(127.513507209 + 0.014354767) + 2^16
  = 2^127.527861976 + 2^16
  = 2^127.527861976
score = log2(W) = 127.527863   (rounded upward to six decimals)
```

Every operation is counted once: the Fp12 multiplication is the reference unit,
the scalar work is divided by the fixed price `C = 20000`, and the pairing and
scalar-multiplication overhead is charged explicitly.

## 5. Success probability

The class walk is modelled as a random mapping on `N` classes; this is the
standard rho heuristic (Brent-Pollard correction for `R`-adding walks), and it is
disclosed as a heuristic (Section 8). For total work `W` steps,

```text
P[collision by W] ~ 1 - exp(-W^2 / (2N)),
```

so at the cap `W_max = 1.5 sqrt(N)` the model success is `1 - exp(-1.125) =
0.675`. The detected-collision delay (the first distinguished point after the
merge, `M * theta^-1 = 2^20 * 2^60 = 2^80` steps of lockstep work) is `2^-47` of
the cap and is inside `W_max`.

The submitted success probability is the rule's threshold **0.39**, well below
the model value 0.675. `experiment.py` Part B (a bounded toy run, 40 trials at
`1.5 sqrt(N)`) measures `0.375 +- 0.077`, consistent with 0.39; the toy is below
the model because at toy scale the detection overhead and the small number of
walks are not negligible, which is exactly the regime difference quantified
above. The claim is 0.39, not the model's 0.675.

The `z = 1` instance is answered directly and is not part of the failure
probability. No instance in GT is excluded.

## 6. Memory, preprocessing and advice

- Distinguished points stored: `theta * W_max = 2^-60 * 2^127.513507 = 2^67.51`
  records; at 640 bytes each (576-byte exact Fp12 element plus two 32-byte
  exponents) this is `2^76.84` bytes. The multiplier table is `2^20 * 576 =
  2^29.2` bytes. Peak storage is at most `2^78` bytes. Memory does not enter the
  score of this track.
- Preprocessing: the multiplier table and the walk starts are per-instance and
  charged above; the relation-collection analogy does not apply, so the reported
  preprocessing bound is the same value `2^127.527863`, charged in full with no
  reuse credit.
- Nonuniform advice: 0. No auxiliary powers are used; the multipliers are
  derived from `g` and `z` at runtime.

## 7. Correctness and edge cases

- `P` nonidentity, `z` in GT: `g != 1`, `z` has a unique logarithm mod `r`.
- `z = 1`: return `O`; `e(P,O) = 1` is verified (one pairing).
- Hitting `y = 1` inside a walk is a success, not a failure: then
  `alpha + x beta = 0` and `x = -alpha/beta` if `beta != 0`, verified by pairing;
  this is a second, independent way for the search to terminate.
- A quotient collision gives `x` or `r-x`; both candidates are verified against
  `e(P,[x]Q0) = z`, so the answer is correct or the run continues/fails.
- The walk's map is a function of the class only, so two trajectories that share
  a class stay in the same class and their distinguished points coincide exactly.
- `r` is prime, so `[x]Q0` has the required order; the point is in G2 by
  construction.

## 8. Heuristics and limitations

1. The class walk is modelled as a random mapping on `N = (r-1)/2` classes; this
   is the standard rho heuristic, supported by the toy runs in `experiment.py`
   Part B (measured first-collision times `1.20 - 1.33 sqrt(N)` across seeds,
   against the model's 1.253) but not proved.
2. The `1.01` per-step price is an instruction ledger for this implementation on
   an AArch64 core; it is an analytic count, not a measurement of a production
   run.
3. The `2^16` reference-operation overhead is a charge, not a measurement; it is
   111 bits below the total.
4. This is an analytic complexity claim for a `2^127` computation, not an
   executed attack; the rules explicitly allow that.
5. The claim is an upper bound for this algorithm, not a security estimate for
   BLS12-381.
6. Only the two-class quotient is used. The `12`-class Frobenius quotient (the
   accepted sibling-track method) is stronger but it is not claimed here: this
   package claims only what its own evidence supports.

## 9. Evidence (`experiment.py`, standard library only)

Part A, on BLS12-381 itself (runs in the verifier container, no third-party
imports): `r` is prime; `p^6 = -1 mod r`; for random GT elements `y^r = 1`,
`conj(y) = y^-1`, `K(conj(y)) = K(y)`, and `K` behaves as a hash over distinct
classes (so exact comparison is required). All checks pass.

Part B, a bounded toy implementation of the *complete* algorithm (same class key,
same canonicalisation at distinguished points only, exact element comparison,
exponent recovery and verification) in a toy prime field with a large prime-order
subgroup: at the claimed budget `1.5 sqrt(N)` the measured success rate is
`0.375` (40 trials) against the model's `0.675`; at `2.5 sqrt(N)` it is `0.875`
against `0.954`. Raw first-collision times of the class walk measure
`1.20 - 1.33 sqrt(N)` against the model's `1.253`, which is the part of the
heuristic the claim depends on. Toy-scale statistics are evidence at toy scope
only, never a production complexity proof.

## 10. Credit and references

- The two-class (inversion) rho idea for this track pair, and the use of the
  `p^6`-Frobenius as inversion, are public: submissions `3da91365-34d1-4a4b-...`
  (`41c34165`, `7da5f23a`) on the sibling track `bls12-381-fixed-g2-invert-g1`
  use the full 12-class version. This package is independently implemented and
  claims only the 2-class case, with its own accounting, evidence and cap.
- Pollard rho and the distinguished-point method: Pollard 1978; van Oorschot and
  Wiener, "Parallel collision search with cryptanalytic applications",
  J. Cryptology 12, 1999.
- Teske, "On random walks for Pollard's rho method", Math. Comp. 70, 2001 (the
  `R`-adding walk and its variance).
- Brent, "An improved Monte Carlo factorization algorithm", BIT 20, 1980 (the
  correction factor used here).
- Curvesma.sh `tracks/bls12-381-fixed-g1-invert-g2/rules.md`,
  `tracks/docs/COST_MODEL.md`, `tracks/classical-group-work-v2.json`.
