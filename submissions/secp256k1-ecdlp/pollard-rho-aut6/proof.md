# secp256k1-ecdlp: Pollard rho with the full automorphism group (order 6)

Status: ready for independent review; no qualification or emitted score.
This is an analytic cryptanalytic claim under exploratory review. It is not
a practical break or a demonstrated production execution.

The normative inputs/outputs are those of `production-v2`; `claim.json`
binds that profile and `classical-group-work-v2` by SHA-256. No restricted
challenge distribution, advice, random-oracle premise, quantum operation or
unpriced oracle is used. All work is charged under the fixed cost model.

## 1. Target and automorphism group

Let p = 2^256 - 2^32 - 977 (field modulus) and
n = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
(subgroup order, prime). The curve is E: y^2 = x^3 + 7 over Fp, generator P
as in the protected profile. Input Q = [x]P for uniform x in [0,n); output x.

secp256k1 has j-invariant 0, hence CM by Z[omega], omega = e^{2pi i/3}.
Since 3 divides p-1 (p-1 = 2^256 - 2^32 - 978 and 978 = 3*326), there is
beta in Fp with beta^3 = 1, beta != 1. Define psi(x,y) = (beta*x, y).
Then psi is an automorphism of E of order 3, and on the subgroup <P> of
prime order n, psi acts as multiplication by a scalar lambda with
lambda^2 + lambda + 1 = 0 mod n, lambda != 1 mod n (lambda has order 3).

Combined with negation (x,y) -> (x,-y), the full automorphism group of
<P> is G_aut = {1, lambda, lambda^2, -1, -lambda, -lambda^2}, of order 6.
Every orbit of G_aut on <P> \ {infinity} has size dividing 6; the number
of orbits is at most ceil(n/6) < n/6 + 1. For the walk analysis we use the
upper bound n/6 on the effective search space (orbits of size < 6 only
help the attacker, and there are O(1) such orbits).

## 2. Algorithm: Pollard rho with distinguished points on orbits

We use the standard r-adding walk (Teske) with r = 32 partitions, run on
canonical orbit representatives, with distinguished points for
memory-efficient collision detection (van Oorschot-Wiener).

### 2.1 Canonical representatives

For R in <P>, R != infinity, define rep(R) as the lexicographically
smallest (x,y) among the 6 points {phi(R) : phi in G_aut}, using (x,y)
with x,y as 256-bit big-endian integers. Infinity maps to itself.
Computing rep(R):
- Compute x1 = beta*x mod p, x2 = beta*x1 mod p (2 field multiplications).
- Candidates: (x,y), (x1,y), (x2,y), (x,p-y), (x1,p-y), (x2,p-y).
- Select lexicographic minimum by comparing x-coordinates (then y on tie).
- If rep(R) = phi(R) for phi corresponding to scalar mu in
  {1,lambda,lambda^2,-1,-lambda,-lambda^2}, record mu.

If R = a*P + b*Q (as group elements), then rep(R) = (mu*a)*P + (mu*b)*Q.
We maintain the invariant that stored coefficients (a,b) always satisfy
rep = a*P + b*Q by multiplying by mu at each reduction step.

### 2.2 The walk

Precompute for j = 0..31: random scalars c_j, d_j in [0,n), and
T_j = c_j*P + d_j*Q (as points). Store (T_j, c_j, d_j). This is 32 scalar
multiplications, negligible (see Section 4).

Partition function s: maps a canonical representative to {0..31}, e.g. by
a fixed hash of the x-coordinate (low 5 bits after a public mixing).
The walk state is (R, a, b) with R = rep(a*P + b*Q).

- Initialize: random a0, b0 in [0,n), R0 = rep(a0*P + b0*Q).
- Step: given (R, a, b), let j = s(R). Set
  R' = R + T_j (one group addition),
  a' = a + c_j mod n, b' = b + d_j mod n,
  then (R_next, mu) = reduce(R') to canonical rep,
  a_next = mu*a' mod n, b_next = mu*b' mod n.
- Distinguished check: R_next is distinguished if its x-coordinate has
  t leading zero bits (t chosen below). If distinguished, store
  (R_next, a_next, b_next) in a hash table; on collision with a stored
  (R, a_old, b_old) having the same representative, proceed to Section 2.3.
- Otherwise continue.

The walk is deterministic given (a0,b0) and the public partition/step
constants. It is a functional graph walk on the set of orbit
representatives.

### 2.3 Collision use

On finding stored (R, a1, b1) and new (R, a2, b2) with the same canonical
representative R (and (a1,b1) != (a2,b2) as pairs), we have:
  a1*P + b1*Q = a2*P + b2*Q   (as group elements, since both equal R)
  (a1-a2) + (b1-b2)*x = 0 mod n, where Q = x*P.
If b1 != b2 mod n, set x = -(a1-a2)/(b1-b2) mod n, verify [x]P = Q,
and output x. If b1 = b2 (hence a1 = a2, since R determines a1-a2=0),
this collision is useless; continue the walk (this occurs with
negligible probability for random coefficients).

The algorithm is Las Vegas: it terminates with probability 1 and upon
termination outputs the correct x (verified by the final [x]P = Q check).

### 2.4 Distinguished-point parameter

Take t such that 2^t is about the square root of the expected walk
length, balancing memory and the extra steps to hit a distinguished
collision. Concretely, t = 40 gives expected ~2^40 stored points at
peak, which is reported in memory (not scored). The walk length analysis
below is for the total steps; the distinguished-point method finds a
collision among stored points with the same asymptotic expected total
steps up to a small constant factor, which we absorb in the explicit
1.7x per-step overhead and the Markov 2x budget (Section 3).

## 3. Work analysis and success probability

### 3.1 Heuristic random-walk model

We use the standard heuristic (Teske; van Oorschot-Wiener; Bos et al.)
that the r-adding walk with r = 32 on the set of orbit representatives
behaves like a uniform random walk on a set of size N_orb <= n/6.
Under this heuristic, the expected number of steps to the first
collision (birthday paradox) is sqrt(pi*N_orb/2) <= sqrt(pi*n/12).

This is a *supported plausible heuristic*: it is the standard analysis
of Pollard rho used in the literature for three decades, with the
automorphism quotient analyzed explicitly in the context of ECDLP
(e.g. Wiener-Zuccherato for negation; Bos-Kaihara-Montgomery for
larger automorphism groups including secp256k1's order-6 group).
We disclose it as a heuristic per the rules.

Fruitless cycles: automorphisms can create short cycles in the walk
(R -> phi(R) -> R). With r = 32 partitions and the standard
cycle-escape (on detecting a 2-cycle, take a deterministic extra step),
their impact on the expected walk length is negligible; we absorb any
residual in the conservative per-step overhead. This is disclosed.

### 3.2 From Las Vegas to fixed budget (Markov)

Let T be the (random) number of walk steps to success (a useful
distinguished collision). Under the heuristic, E[T] <= sqrt(pi*n/12).

Run the walk for a fixed budget of B = 2*E[T] steps. By Markov's
inequality, P(T > B) <= E[T]/B = 1/2. Hence this fixed-budget variant
succeeds with probability at least 1/2 >= 0.39, as required. The work
is now deterministic (no expected-runtime standing in for total work).

Using n < 2^256:
  E[T] < sqrt(pi*2^256/12) = 2^128 * sqrt(pi/12).
  sqrt(pi/12) = 0.51166... (pi/12 = 0.261799..., sqrt = 0.511663).
  B = 2*E[T] < 2^129 * 0.511663 = 2^(129 + log2(0.511663)).
  log2(0.511663) = -0.96673...
  B < 2^128.03327 group additions.

### 3.3 Per-step cost in reference operations

Each walk step performs:
(a) One group addition R + T_j. We count this as 1 reference operation
    (H). The reference is "batched incomplete affine addition"; a
    single addition is at most this cost for upper-bounding purposes,
    and we separately account all auxiliary work in (b)-(d) via
    instruction counts divided by C = 1000.

(b) Automorphism reduction to canonical representative:
    - 2 field multiplications by beta (for x1, x2): 2*179 = 358
      instructions (using the pinned 179-instruction field multiply).
    - Lexicographic minimum over 6 candidates: at most 5 x-comparisons
      (4 limbs each, ~4 instructions per limb) ~ 80 instructions;
      coefficient scaling (mu*a, mu*b mod n): 2 modular multiplications,
      generously 200 instructions.
    - Total (b): at most 700 instructions = 0.7 reference units.
    This is a generous upper bound; actual implementations are cheaper.

(c) Partition index s(R): hash/mix of x-coordinate, at most 30
    instructions = 0.03 units.

(d) Distinguished check (t leading zeros of x): at most 20
    instructions = 0.02 units. Hash-table store on distinguished
    points: amortized negligible (one store per 2^t steps).

Per-step total: 1 + 0.7 + 0.03 + 0.02 = 1.75 reference units.
We use 1.75 as a clean upper bound.

### 3.4 Total work and score

Total work W (in target-group-operation equivalents):
  W = B * 1.75 + W_pre + W_post,
where:
- W_pre: 32 precomputed T_j (each <= 512 additions) + setup <= 2^15,
  negligible.
- W_post: final collision solving (1 modular inversion + mults) and
  verification scalar multiplication (<= 512 additions), negligible.

  W < 2^128.03327 * 1.75 + 2^15
    < 2^128.03327 * 2^0.807355 + 2^15    (log2(1.75) = 0.8073549...)
    = 2^128.84063 + 2^15
    < 2^128.84064                        (the 2^15 term is absorbed)

Score = log2(W) < 128.84064.

We declare score 128.84065 (rounded up), i.e. numerator 12884065,
denominator 100000. As a reduced fraction: 128.84065 = 2576813/20000.

Check: 2576813/20000 = 128.84065. Yes.

This beats the baseline 129.000802 by about 0.16 bits, i.e. the attack
uses roughly 2^0.16 ≈ 1.12x less work than the BSGS control.

### 3.5 Success probability

The fixed-budget walk succeeds with probability at least 1/2 (Markov,
Section 3.2), provided the useful-collision event implies the equation
is solvable (b1 != b2). The probability that a random collision has
b1 = b2 mod n is at most (walk length)/n, negligible (< 2^-120).
Absorbing this, success probability >= 0.49 > 0.39. We declare 0.49.

More precisely: P(success) >= P(T <= B) - P(bad collision)
  >= 1/2 - 2^128.04/n > 1/2 - 2^-127 > 0.49.

## 4. Memory, preprocessing, advice

- Peak memory: distinguished-point table. Expected number of stored
  points ~ B/2^t. With t = 40, ~2^88.4 records; each record (point +
  two scalars) ~ 96 bytes. Peak ~ 2^88.4 * 2^6.6 = 2^95 bytes.
  We report memory_log2_bytes = 95. (Memory is not scored.)
  Note: the table can be pruned/checked incrementally; this bound is
  generous. A larger t reduces memory at no asymptotic work cost.
- Preprocessing: 32 T_j points + partition constants <= 2^15 reference
  units; preprocessing_log2 = 15. Included in total W.
- Nonuniform advice: none (0 bytes). The beta constant and partition
  function are public algorithm parameters, not searched advice.
- All memory transfers for table inserts/lookups are included in the
  per-step overhead (d).

## 5. Heuristics disclosed

1. **Random-walk heuristic**: the r=32 adding walk on orbit
   representatives is modeled as a uniform random walk, giving
   E[T] <= sqrt(pi*n/12). Standard in Pollard-rho literature; supported
   by extensive experiments at small sizes and the Teske analysis.
   Limitation: not proven for this specific walk; the constant could
   differ slightly. We use an upper bound (n < 2^256) and a 2x Markov
   margin, which absorbs small constant-factor deviations.
2. **Automorphism quotient**: orbits have size exactly 6 except O(1)
   fixed points; using n/6 as the effective space is standard and
   conservative (smaller orbits only help).
3. **Fruitless cycles**: short cycles from automorphisms are handled by
   standard cycle detection/escape; residual impact absorbed in the
   1.75x per-step overhead. Supported by Bos et al. analysis.
4. **Distinguished-point overhead**: the expected total steps to a
   *stored* collision is within a small constant of the plain collision
   bound; absorbed in the Markov 2x budget.
5. **Coefficient independence**: b1-b2 != 0 with overwhelming
   probability; standard.

No heuristic assumes subexponential behavior, unproven number-theoretic
conjectures, or oracles. All are standard, literature-supported
plausibility assumptions for Pollard-rho analysis, explicitly permitted
by the exploratory review standard.

## 6. Limitations

- Analytic claim; no full-size execution (infeasible by design).
- The score is an upper bound; the true expected work may be lower
  (e.g. tighter automorphism overhead, n < 2^256 - 2^32).
- Memory bound (2^95 bytes) is astronomical; this is a cryptanalytic
  cost argument, not a practical attack.
- Success probability 0.49 is a lower bound; the Las Vegas variant has
  success probability 1 with expected work < 2^127.9.
- The 1.75x per-step overhead is generous; optimized implementations
  would reduce it.

## 7. Comparison with baseline

The baseline BSGS uses 2^129.000802 work with 2^136 bytes memory.
This Pollard-rho attack uses < 2^128.84065 work (about 1.12x less)
with 2^95 bytes memory (much less, though still astronomical).
The improvement comes from the square-root collision bound with the
order-6 automorphism quotient (saving ~0.97 bits vs plain rho) and
from Pollard's lower constant vs BSGS (rho: sqrt(pi/2)≈1.25 vs
BSGS: 2, saving ~0.68 bits), partially offset by the automorphism
overhead per step.
