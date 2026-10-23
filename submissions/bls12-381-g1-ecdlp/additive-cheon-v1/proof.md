# Additive Cheon: one addition per baby/giant step

**Claim.** For every n with 2^20 ≤ n < r, the deterministic algorithm below recovers s with heuristic
success probability 0.4055 (≥ 0.39 required). Its total work is W(n) = W(2^20) ≤ **2^116.982765**
target-group-operation equivalents.
- Preprocessing, meaning the challenge-independent work, is ≤ 2^115.986 and is included in W. It
  consists of the baby walk with its hash table, the per-prime tables, the characters and the sieves.
- Peak memory is ≤ 2^122.9 bytes. There is no advice.
- The algorithm reads only A_1 = T and A_d, with d = 2^20.

The idea is Cheon's (p−1)-attack, run as a baby-step/giant-step (BSGS) inside the order-m subgroup H of
d-th powers, with one change. The baby and giant sets are **small integers that lie in H**, not geometric
progressions ζ^i. Small integers listed in increasing order can be walked with **one batched affine
addition per element**. A geometric step needs a scalar multiplication, which is 2 additions with the
Kozaki–Kutsuma–Matsuo (KKM) c = 3 tables.

## 1. Facts about the fixed target

These are checked by `evidence.py` and by direct computation.
- r ≈ 2^254.857 is prime and v₂(r−1) = 32, so d = 2^20 divides r−1.
- H = {x^d : x ∈ F_r^*} = ker(x ↦ x^m) has order m = (r−1)/d ≈ 2^234.857, and 6 | m.
- **The hidden value is uniform on H.** t := s^d is uniform on H, because s is uniform on F_r^* and
  x ↦ x^d is d-to-1 onto H. A_d = [t]P0 is public.
- **The symmetry group sits inside H.** Let U = {±1, ±λ, ±λ²}, where λ = x² − 1 is the GLV eigenvalue
  (λ² + λ + 1 = r) with [λ](x, y) = (βx, y). U ⊂ H because its order 6 divides m.
- **One key identifies a whole orbit.** All six points [u]Q, for u ∈ U, have y-coordinate ±y(Q),
  because λ fixes y. Conversely, the points with y-coordinate ±y₀ are exactly the six (βʲx₀, ±y₀). So
  canonical(Q) := min(y, p − y) identifies the U-orbit of Q.
- **Membership in H is a character test.** χ(x) := log_ζ(x^m) ∈ Z/d is a character with kernel H, for
  a fixed generator ζ of μ_d. So a·b ∈ H ⇔ χ(a) + χ(b) ≡ 0 (mod d).
- **2 ∉ H.** Indeed 2^m ≠ 1 (mod r); 2^m has order 2^16.

## 2. Algorithm

Parameters: K = 0.52, X1 = 2^76, X2 = 2^70 < X1/2, GMAX = 2^30 and N_B = N_G = N = ⌈√(K·m/6)⌉ ≈ 2^115.664.

1. **Precompute** the characters, sorted into class lists.
   - Compute χ(q) for the primes q ≤ X2 and for the primes used in (X1/2, X1]. Each costs one 235-bit
     exponentiation in F_r and a Pohlig–Hellman logarithm in μ_{2^20}.
   - Sieve χ(b) for all b ≤ X2 by multiplicativity.
   - Stream each b into its class list L_c = {b ≤ X2 : χ(b) = c}, sorted in increasing order.
2. **Baby walk.** For fresh primes a ∈ (X1/2, X1]:
   - set Q_a = [a]P0 and build T[k] = [k]Q_a for k ≤ GMAX (T[2] by doubling, then T[k] = T[k−1] + Q_a);
   - walk b over L_{−χ(a)} in increasing order, with [ab]P0 = [ab_prev]P0 + T[b − b_prev];
   - for a gap greater than GMAX, first add T[GMAX] repeatedly;
   - insert key canonical([ab]P0) with locator (a, b) into a hash table.

   Stop after N elements. Every element a·b lies in H. They are distinct integers below r, because a is
   the unique prime factor of a·b above X2.
3. **Giant walk.** Do the same with disjoint primes a′ and base A_d, which gives the points [a′b′·t]P0.
   Look up the canonical key of each point. Stop after N elements, or earlier on success.
4. **Collision.** A key match means [g·t]P0 = [u·ab]P0 for some u ∈ U, where g = a′b′. Test the six
   candidates t = u·ab/g by checking [t]P0 = A_d.
5. **Recover s.** Take one d-th root s0 of t: twenty square roots, keeping at each step the root that is
   still a 2^j-th power. Then s = s0·ω^i for some i < 2^20, where ω generates μ_d. Test [s0·ω^i]P0 = T
   for each i and output s.

**Batching and representation.** Each prime's walk, and its table build, is sequential. The 2^16
primes in flight are therefore advanced in lockstep. Every round performs 2^16 independent incomplete
affine additions, one per prime, sharing a single Montgomery batch inversion; this is the reference
operation's batched setting. The inversion costs ≈ 100 units per round, about 0.0015 units (3
instructions) per element, which fits inside the walk-bookkeeping budget. All stored points, meaning
tables, walk states and Q_a, are affine: 96 bytes each.
- Q_a = [a]P0 (or [a]A_d) is computed in Jacobian coordinates and normalised by one inversion, charged
  in the per-prime scalar budget.
- T[2] is one affine doubling, also charged there.

**Incomplete additions.** A walk step adds [b_prev]Q and [gap]Q, with Q ≠ O and b = b_prev + gap < r. The
x-coordinates can match only if b_prev ≡ ±gap (mod r).
- b_prev ≡ −gap would force b ≡ 0, which is impossible.
- b_prev = gap means b = 2·b_prev in the same class list, which would need χ(2) = 0. Since 2 ∉ H, this
  never happens.
- *Overflow chunks.* [b_prev + i·GMAX]Q + [GMAX]Q collides only when b_prev = GMAX and i = 0. That is at
  most once per prime; it is done by a doubling, charged to the per-prime scalar budget. The remainder
  addition is safe because b_prev + j·GMAX > GMAX ≥ rem for j ≥ 1.
- *Table build.* T[k] = [k−1]Q + Q is safe for k ≥ 3.
- Giant points are multiples of Q = [a′t]P0 ≠ O, so the same argument covers them.

## 3. Success probability

Steps 4–5 are deterministic. The algorithm succeeds exactly when t ∈ S := U·B·G^{-1}. Because t is
uniform on H, Pr[success] = |S|/m.

- **H1 (coverage).** Apart from the structural loss below, the 6N² products u·b·g^{-1} coincide no
  more often than random products would. Then |S|/m ≈ 1 − e^{−K} = 0.4055. Two kinds of coincidence
  are ruled out exactly:
  - A pair with the same baby prime a and the same giant prime a′ that coincides with another such pair
    would need an integer identity below r, which is impossible even with u = −1.
  - Products built from different primes cannot be equal as rationals, because a is determined.

  *Structural loss.* Rationals repeat when b/b′ = c/c′, which happens through a shared factor
  g = gcd(b, b′) and a cofactor k with χ(k) = χ(g). The expected number of repeats per pair is at most
  2·E[gcd]/d ≈ 2·(6/π²)(ln X2 + 0.15)/d ≈ 2^−14.1. That changes the success probability by less than
  10⁻⁴. The margin is K = 0.52 against the 0.494 needed for 0.39.
- **H2 (equidistribution).** The values χ(a) of the primes used are spread evenly over Z/d. So the
  average length of the class lists actually walked is ≈ X2/d = 2^50; the lists' total length is
  exactly X2. This affects only the number of primes needed (about 2^67.7, charged at twice the
  expected number). There are 2^69.3 primes available in (X1/2, X1], which is ≥ 3× slack. H2 does not
  affect success.

Evidence: `evidence.py` with the declared inputs. These results hold at toy scale only.

- **Coverage.** Exact |S|/m for a toy prime r ≈ 2^28 with d = 2^8, K = 0.52 and prediction 0.4056:
  - random sets: 0.4061;
  - smallest integers in H: 0.4069;
  - product construction at X2/d = 4, 16, 64 and 256: 0.4038, 0.4020, 0.4028 and **0.4056**.

  The measured shortfalls from the prediction are 0.0018, 0.0036, 0.0028 and 0.0000. They are small
  (≤ 0.0036, well inside the 0.0155 margin) and shrink as X2/d grows. The gcd formula is only a
  conservative upper bound on the loss: it is not a model of these deficits, and it grows with ln X2.
  The class-length spread also narrows as X2/d grows (0.76–1.22 × mean at X2/d = 256). Production has
  X2/d = 2^50.
- **Toy attack.** The whole algorithm on a 24-bit j = 0 curve, with d = 2^6 and GMAX = 16, so that the
  overflow path runs (1.3·10⁶ chunk additions; the doubling case b_prev = GMAX did not occur).
  - 1,500 trials gave success 0.3953 against 0.4095 predicted, about 1.1σ low, at X2/d = 16 where a
    small toy-scale loss is expected. The `conservative_gcd_bound_toy_scale` field in the output (0.366)
    applies at toy scale only; at production scale the same bound is 2^−14.1.
  - Every recovered scalar was correct.
  - Earlier runs of 500 trials on seeds 1–3 gave 0.404, 0.420 and 0.450.

## 4. Cost (unit: one batched incomplete affine G1 addition = 2,000 instructions)

Per-element budgets, in instructions on top of the reference addition.

| Item | Budget | Itemised count |
|---|---:|---|
| canonical y | 40 | p − y (≈ 7), 6-limb compare (≈ 12), select (6): ≈ 25 |
| hash insert (baby) / lookup (giant) | 300 | see the breakdown below: 61 + 18 per routing hop, ≤ 13 hops |
| walk bookkeeping | 150 | read b and the gap (4), overflow test (2), address T[gap] (3), load T[gap] (12), load and store points (24), batch-inversion buffers beyond the 6M (24), loop (5): ≈ 75 |

**Hash breakdown.** The 2^116.7-slot table, at load ≤ 1/2 with 72-byte records, cannot sit in one
address space. It is sharded over nodes by the high bits of the key. Each insert or lookup:
- computes the shard and slot from key bits: 4;
- appends the 9-word record to that shard's outbound buffer: 14;
- is forwarded through h hops, each loading and storing every word: 18 per hop;
- is read at the receiver: 9;
- is probed by linear probing at load ≤ 1/2, ≤ 2.5 expected probes of 6 instructions: 15;
- stores the record, for inserts: 9;
- costs loop and control: 10.

Apart from the hops, the work is 61 instructions. The budget of 300 therefore covers h ≤ 13 hops. A
radix-32 routing network over ≤ 2^60 nodes needs 12 hops. Lookups follow the same route, and a match reply is
needed only on success. These are expected costs summed over 2^117 operations, where linear probing and
balanced sharding concentrate sharply. The χ-sieve is segmented over ranges of b, so each segment is
local, and appends to the class lists are streamed. Walks read one list sequentially.

This gives c_B = c_G = 1 + 490/2000 = 1.245. The full ledger, recomputed exactly by `evidence.py`:

| Term | log₂ (units) |
|---|---:|
| baby walk N·c_B | 115.981 |
| giant walk N·c_G | 115.981 |
| gap-overflow additions (deterministic: Σ gaps ≤ X2 per prime) | ≤ 107.7 |
| per-prime tables T[k] (2^30 additions × 2^67.7 primes) | 97.7 |
| per-prime scalar multiplications (3 × 1,000 units per prime) | 79.2 |
| χ for primes ≤ 2^76 (100 units = 2·10⁵ instructions each: exponentiation plus Pohlig–Hellman) | 77.3 |
| prime sieve to X1 (64 instructions per integer); χ-sieve and class lists to X2 (200 per b) | 71.1 |
| collision checks, d-th root, μ_d search, input checks | ≤ 31.6 |
| **total W(2^20)** | **116.982765** |

- A scalar multiplication is charged at 1,000 units for a ≤ 76-bit scalar and 3,000 units for a 255-bit
  scalar. These are generous Jacobian ladders at ≤ 5 units per step, plus one inversion at 100 units.
- The score is rounded up: 116.982765 = 23396553/200000.

**Memory.**
- Hash table: 2N slots × 72 bytes ≈ 2^122.8 bytes.
- Class lists: ≈ 2^73.2 bytes.
- Walk state and affine tables for 2^16 primes in flight: 2^16 · 2^30 · 96 bytes ≈ 2^52.6 bytes.
- Memory is not scored. Any BSGS form of Cheon stores one side, about 2^115.7 entries.
- Rho-type Cheon variants need little memory, but each of their steps is a multiplication in the
  exponent (≥ 7 additions with KKM tables).

## 5. W(n) and applicability

The algorithm uses only A_1 and A_{2^20}, so W(n) = 2^116.982765 for every 2^20 ≤ n < r. Nothing is
claimed for n < 2^20. More generally, for any d | r−1 with 6 | (r−1)/d and d ≤ n:
- W ≈ 2·c·√(K(r−1)/(6d)), plus the lower-order terms above;
- X2 = d·2^50.

## 6. Context (unscored estimates)

| Method | log₂ W |
|---|---:|
| Pollard rho on G1, order-6 classes: 0.5·log₂(πr/12) additions, plus overhead | ≥ 126.5 |
| Cheon BSGS, double-and-add scalar multiplication per step (≈ 383 additions) | ≈ 125.2 |
| Cheon BSGS with KKM c = 3 tables (2 additions plus one F_r multiplication per step); same K and per-element budgets as above, computed by `evidence.py` | 117.894 |
| **this work** | **116.983** |

- The gain over KKM-accelerated Cheon is a constant factor, about 0.91 bit.
- **Relation to the generic bound.** Boneh–Boyen's generic lower bound Ω(√(ε·r/d)) has unstated
  constants. Here the walks use 2N ≈ 2^116.66 additions, and the generic coverage is
  N²/m = K/6 ≈ 0.087. The factor 6 that lifts success to 1 − e^{−K} comes from the cheap action of U,
  through negation and the GLV endomorphism. The model does not capture that action, so there is no
  conflict with the bound.
- **Context, not a premise.** KKM (Pairing 2007, paraphrased) set aside a table of small-integer
  multiples as hard to use, because such integers are scattered among the exponents. The product
  construction a·b, with χ(a) + χ(b) ≡ 0, is one way to use them. We found no prior description of
  it, but nothing in the claim depends on novelty.

**Disclosed heuristics:**
- H1, with its structural loss bounded;
- H2;
- expected-cost hashing and sharding.

No other assumptions.
