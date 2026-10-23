# Additive Cheon v25: prime class lists — integer-level coincidences eliminated by construction

**Lineage.** Provenance: the construction follows the published `additive-cheon-v1`
package (submission cbadde77, in review on this track, declared as a dependency).
It is restated in full below, and every identity, instruction itemisation, bound,
price and toy number in THIS proof is this package's own, re-derived here and
re-verified by the runs of THIS package's `evidence.py`; nothing is inherited by
citation. The accounting differences against the provenance package, carried
since v2:

1. canonical y: itemised below at ~25 instructions, charged **30**.
2. walk bookkeeping: itemised below at ~75 instructions, charged **100**.
3. hash insert/lookup: a complete per-transfer accounting, itemised below:
   endpoint work 61 plus six radix-1024 relay stages of 29 instructions over the
   ≤ 2^60 nodes, 235 itemised, charged **240**; evidenced by a 16,384-key probe
   and shard check on structured production-style walk keys.
4. the success statement is floor-rounded with a heuristic allowance (section 3);
   the deterministic set rule is stated (section 2) and exercised by a
   fixed-order attack run.

**The v23 change against the v22 lineage: the b's are the PRIMES of the
χ-class lists.** Twenty-two revisions measured the integer-list construction's
coverage deficit across every controlled dimension — d-ladders, population and
truncation controls, matched ratio pairs — and the allowance remained an
extrapolation to production's d = 2^20. v23 removes the deficit's dominant
family by construction instead of extrapolating it: with every walked element a
product a·b of two primes from disjoint ranges, integer-level coincidences are
impossible off the diagonal by unique factorisation (section 3), and the one
family that survives is counted exactly at toy scale and bounded — not
extrapolated — at production. The parameters rebalance for the thinner prime
lists: X1 = 2^82 (N_A ≈ 2^73.3 primes at 16× availability slack) and
GMAX = 2^34 (keeping the gap-overflow and table terms at 2^109.3 / 2^107.3).
The walks, batching, symmetry key, recovery and the one-addition-per-element
structure are unchanged.

**The v24 changes against the v23c judgment:** the coverage sweep implements
the declared block rule exactly (ka = 2n·d·LN/X2 + 1, no fallback), the modular
premise is measured as an exact-rule d-ladder instead of asserted, and the
hosted memory blow-up is fixed (the 78 MB coverage bitmap is converted in 1 MB
chunks, not one 78 MB integer; peak ≈ 100 MB against the 256 MB cap; the
evidence runs in 4.1s local, inside the 20-second limit that the v22-lineage
package passed at 5.0s local).

**The v25 changes against the v24 judgment (three items).**
1. **The shared-prefix accounting is corrected to lost distinct values.**
   A same-class (baby prime a, giant prime a′) pair whose consumed prefixes
   share t elements has its t diagonal pairs (x = a·b, y = a′·b) all map to the
   single value u·a/a′ — for each u — losing 6·(t−1) distinct values, not
   C(t,2) pairs of pairs. The per-run counter now reports exactly that
   (`prefix_lost_values`, as a fraction of the 6n² pair slots), and the
   production bound becomes 6·(⌈N_A/2⌉²/d)·(L2−1) = **2^172.6 values =
   2^−61 of the 6N_B·N_G pair slots** — 42 orders of magnitude below the
   v24 bound, because the diagonal collapse is linear, not quadratic, in the
   prefix length.
2. **Every measured deficit is at or below the exactly counted integer-level
   loss.** At each of the four prime-construction coverage points the measured
   deficit and the counted lost values are both reported, and the deficit never
   exceeds the count (at the two largest d the measurement is above prediction)
   — an arithmetic inequality between two reported numbers, with no statistical
   model invoked: there is no unexplained modular excess at any measured point.
   The ratio point below additionally shows X2/d does not drive the deficit.
3. **K moves 0.52 → 0.55 and the allowance widens 0.014 → 0.032.** The margin
   the old K left unused is converted into claim robustness: the model value
   1 − e^−0.55 = 0.42305 floor-rounds at 0.423, and 0.423 − 0.032 = **0.391**
   ≥ 0.39 with the allowance now **12.3×** above the worst measured deficit
   (+0.0026) rather than 5.4×. The score moves 116.918317 → **116.958777** —
   still 0.024 under the frontier's 116.982765, because every cost term scales
   as √K.

Every toy number below is reproduced by this package's own code and seeds.

**Claim.** For every n with 2^20 ≤ n < r, the deterministic algorithm below recovers s with heuristic
success probability at least 0.391 (≥ 0.39 required; the model value 1 − e^−0.55 = 0.42305 floor-rounded at
0.423, minus a 0.032 heuristic allowance for the production coverage deficit — section 3). Its total work
is W(n) = W(2^20) ≤ **2^116.958777** target-group-operation equivalents (exact score 116958777/1000000,
rounded up).
- Preprocessing, meaning the challenge-independent work, is ≤ 2^115.968 (exactly the ledger's own
  preprocessing output) and is included in W. It consists of the baby walk with its hash table, the
  per-prime tables, the characters and the sieves.
- Peak memory is ≤ 2^122.9 bytes. There is no advice.
- The algorithm reads only A_1 = T and A_d, with d = 2^20.

The idea is Cheon's (p−1)-attack, run as a baby-step/giant-step (BSGS) inside the order-m subgroup H of
d-th powers, with one change. The baby and giant sets are **small products of two primes that lie in H**,
not geometric progressions ζ^i. Products a·b listed in increasing order can be walked with **one batched
affine addition per element**. A geometric step needs a scalar multiplication, which is 2 additions with
the Kozaki–Kutsuma–Matsuo (KKM) c = 3 tables.

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

Parameters: K = 0.55, X1 = 2^82, X2 = 2^70 < X1/2, GMAX = 2^34 and N_B = N_G = N = ⌈√(K·m/6)⌉ ≈ 2^115.705.

**The deterministic set rule (fixes B and G given the input).** The primes a
are taken from (X1/2, X1] in increasing order, each used exactly once: the
first ⌈N_A/2⌉ form the baby side and the next ⌈N_A/2⌉ the giant side, where
N_A = (N_B+N_G)·2d·49/X2 + 2 ≈ 2^73.32 (the 49 = ln(2^70) + slack buys the
H2 minimum supply per class, below). For each prime a, the b's are **the
primes of L_{−χ(a)} in increasing order** — exactly the order the χ-sieve
streams them — and the walk stops for that prime when its side reaches N
elements. Nothing is sampled or shuffled; B and G are fixed by the input.
The evidence's toy instantiation of this rule is the same formula with the
toy's ln X2 in place of 49: **ka = 2n·d·LN/X2 + 1 with LN = ⌊ln X2⌋ + 1**,
implemented with no fallback, and the attack uses the same form at its own
parameters (ka = 41 there). One attack configuration runs in fixed order
(16 trials) alongside the shuffled runs, which robustness-check the rule
against unfavourable prime blocks.

1. **Precompute** the characters, sorted into prime class lists.
   - Compute χ(q) for the primes q ≤ X2 and for the primes used in (X1/2, X1]. Each costs one 235-bit
     exponentiation in F_r and a Pohlig–Hellman logarithm in μ_{2^20}.
   - Sieve χ(b) for all b ≤ X2 by multiplicativity (the sieve necessarily processes every integer;
     **only the primes are stored to the class lists**).
   - Stream each prime b ≤ X2 into its class list L_c = {b ≤ X2 prime : χ(b) = c}, sorted in
     increasing order. Each list holds ≈ π(X2)/d ≈ 2^44.7 primes (H2).
2. **Baby walk.** For fresh primes a ∈ (X1/2, X1]:
   - set Q_a = [a]P0 and build T[k] = [k]Q_a for k ≤ GMAX (T[2] by doubling, then T[k] = T[k−1] + Q_a);
   - walk b over the primes of L_{−χ(a)} in increasing order, with [ab]P0 = [ab_prev]P0 + T[b − b_prev];
   - for a gap greater than GMAX, first add T[GMAX] repeatedly;
   - insert key canonical([ab]P0) with locator (a, b) into a hash table.

   Stop after N elements. Every element a·b lies in H. They are distinct integers below r: a is the
   unique prime factor of a·b above X2, and with b prime the integer a·b determines both factors.
3. **Giant walk.** Do the same with disjoint primes a′ and base A_d, which gives the points
   [a′b′·t]P0. Look up the canonical key of each point. Stop after N elements, or earlier on success.
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
- b_prev = gap means b = 2·b_prev in the same class list: 2·b_prev is an even integer above 2, hence
  composite, and never in a prime list. (The χ(2) = 0 route is also closed: 2 ∉ H.)
- *Overflow chunks.* [b_prev + i·GMAX]Q + [GMAX]Q collides only when b_prev = GMAX and i = 0. That is at
  most once per prime; it is done by a doubling, charged to the per-prime scalar budget. The remainder
  addition is safe because b_prev + j·GMAX > GMAX ≥ rem for j ≥ 1.
- *Table build.* T[k] = [k−1]Q + Q is safe for k ≥ 3.
- Giant points are multiples of Q = [a′t]P0 ≠ O, so the same argument covers them.

## 3. Success probability

Steps 4–5 are deterministic. The algorithm succeeds exactly when t ∈ S := U·B·G^{-1}. Because t is
uniform on H, Pr[success] = |S|/m.

- **H1 (coverage; heuristic).** The 6N_B·N_G pair products u·b·g^{-1} coincide no more often than
  random products would, so |S|/m ≥ 1 − e^{−K} = 0.42305, stated as the floor 0.423. The production
  coverage deficit is charged at a **0.032 heuristic allowance**; the claimed success is at least
  0.423 − 0.032 = **0.391** ≥ 0.39. The deficit's possible sources split into an integer level (k = 0)
  and a modular level (k ≠ 0); the integer level is provable and the modular level is measured:

  **Theorem (integer-level coincidences, k = 0).** Let u·x₁y₂ = u′·x₂y₁ be an
  exact integer identity with x_i = a_i·b_i ∈ B, y_i = a′_i·b′_i ∈ G, where
  u, u′ are the six fixed representatives of U in (0, r). Each side carries
  exactly two block primes above X1/2 (a_i and a′_j, all > 2^81), and each of
  the six fixed constants u < r ≈ 2^254.9 has at most ⌊254.9/81⌋ = 3 prime
  factors above X1/2.
  - **If u = u′**, unique factorisation gives {a₁, a′₂} = {a₂, a′₁} as
    multisets of the primes above X1/2. The baby and giant blocks are disjoint,
    so a₁ = a′₁ and a₂ = a′₂ are impossible; hence a₁ = a₂ and a′₁ = a′₂.
    Cancelling, b₁b′₂ = b₂b′₁ with all four factors prime, so
    {b₁, b′₂} = {b₂, b′₁}: either the diagonal (the same pair twice), or the
    swap b₁ = b′₁, b₂ = b′₂. In the swap case **both elements of each triple
    use the same small prime**: x₁ = a·b, y₁ = a′·b and x₂ = a·c, y₂ = a′·c
    with b, c in the shared consumed prefix of the same class list (this
    requires χ(a) = χ(a′)) — and the value u·x·y^{-1} = u·a/a′ is then
    **independent of the shared element**. A same-class prime pair with a
    shared prefix of t elements therefore loses exactly **6·(t−1)** distinct
    values (the t diagonal values collapse to one, for each of the six u) —
    the *shared-class-prefix family*, linear in t.
  - **If u ≠ u′**, a block prime must equal one of the at most 6×3 = 18 fixed
    prime factors of the six constants; each such configuration admits at
    most (2π(X2)/d)³ completions in the b's, bounding this exceptional family
    at ≤ 2^151 values — 21 orders of magnitude below the prefix bound, and
    covered by it.
  - **The composite-b multiplicative family is empty.** With integer b's,
    b₁b′₂ = b₂b′₁ has many off-diagonal solutions (2·15 = 3·10 = 5·6); with
    prime b's, unique factorisation leaves only the two cases above. This is
    the family the v2–v22 campaign measured and extrapolated; v23 deletes it.

  **The prefix family, counted exactly at every toy point and bounded at
  production.** The evidence computes, for every coverage point, the exact
  number of lost distinct values `prefix_lost_values` = Σ 6·(t−1) over
  same-class prime pairs with shared consumed prefix t (reported also as a
  fraction of the 6n² pair slots). At production the ledger bounds it at
  6·(⌈N_A/2⌉²/d)·(L2−1) with L2 = 2·π(X2)/d (the H2 concentration bound,
  factor-2 slack): **2^172.6 values = 2^−61 of the 6N_B·N_G pair slots**
  (`log2_prefix_duplicate_values_bound`, `prefix_duplicate_fraction_log2` in
  the ledger output; the latter divides by N_B·N_G and reads 2^−58.8). This is
  a bound under H2, not an extrapolation in d, and it is 42 orders of
  magnitude below the random model's own duplicate rate.

  **Measured coverage under the exact rule** (`coverage_prime_sweep`; every
  point: ascending prime blocks by the declared rule ka = 2n·d·LN/X2 + 1, no
  fallback, X1·X2 < r asserted in code, exact bitmap counts of |U·B·G⁻¹| in
  1 MB chunks; the fraction column is lost values / 6n² pair slots):

  | point | r | n | ka | counted lost values (fraction) | prediction | measured | deficit |
  |---|---|---:|---:|---:|---:|---:|---:|
  | d = 2^4, X2/d = 64 | 2^24.24 | 338 | 74 | 3,858 (0.0056) | 0.4245 | 0.4219 | **+0.0026** |
  | d = 2^4, X2/d = 256 | 2^25.04 | 445 | 32 | 1,896 (0.0016) | 0.4242 | 0.4234 | **+0.0008** |
  | d = 2^6, X2/d = 64 | 2^25.27 | 241 | 68 | 420 (0.0012) | 0.4236 | 0.4244 | **−0.0008** |
  | d = 2^8, X2/d = 64 | 2^29.22 | 473 | 148 | 360 (0.0003) | 0.4238 | 0.4243 | **−0.0005** |

  Three facts read directly off this table, each an arithmetic comparison of
  two reported numbers with no statistical model invoked:
  - **No unexplained modular excess anywhere.** At every point the measured
    deficit is at or below the exactly counted integer-level loss (at the two
    largest d the measurement is *above* prediction). The deficits at the two
    smaller points are the counted structure itself; nothing additional is
    needed to explain them.
  - **X2/d does not drive the deficit.** At fixed d = 2^4, moving X2/d from
    64 to 256 (a 4× thicker class, a different modulus and n) *halves* the
    deficit: +0.0026 → +0.0008, with the counted loss falling in proportion.
  - **d does not drive it either.** |deficit| ≤ 0.0008 at every point with
    d ≥ 2^5 — the worst point is the smallest toy.
  At d = 2^4 the random-set baseline measures 0.4219 (against the 0.4245
  prediction) and the smallest-integers-in-H baseline 0.3952: structured sets
  are deficient, the prime construction sits with the random baseline.

  **The matched control.** At d = 2^4, the SAME r, m, n, X1, X2, ka and prime
  blocks with the class lists holding all integers instead of the primes
  measures **0.3296** (deficit +0.0949) against the prime construction's
  0.4219 (+0.0026): the composite multiplicative family the construction
  eliminates, isolated by direct control (its counted prefix loss is only 624
  values — the deficit is the uncounted multiplicative family).

  **Why the modular residual is heuristically birthday-random: the values are
  squarefree.** Every coverage value is u·x·y^{-1} with x = a·b, y = a′·b′ —
  up to the six fixed constants, a ratio of squarefree integers whose prime
  factors all occur once, from disjoint ranges. The birthday behaviour of
  such ratios in F_r^* is the standard premise of every low-memory DLP
  heuristic on this benchmark family; the claim does not rest on proving it,
  but on what the construction and the measurements contribute: (i) the
  family that measurably breaks the premise is the repeated-factor
  (composite-b) family, and the construction eliminates it (control:
  +0.0949 → +0.0026); (ii) at every measured point the deficit is at or below
  the counted integer-level loss — no unexplained modular excess at any
  point, and the worst total deficit anywhere is +0.0026; (iii) the
  integer-level production loss is bounded, not extrapolated, at 2^−61 of
  pair slots; (iv) the 0.032 allowance is 12.3× the worst measured deficit.
  (As clearly-labelled loose context, the output's integer-list
  repeated-ratio estimate — 0.3572 at the toy attack's scale, 2^−14.1 at
  production — upper-bounds the deleted family and is not used.)
- **H2 (equidistribution of the class primes).** The primes ≤ X2 spread evenly
  over the χ-classes: each class holds ≈ π(X2)/d ≈ 2^44.7 primes, at least
  X2/(2d·49) with ln-slack (the quantity the block rule budgets), so each
  prime a supplies that many walked elements and N_A ≈ 2^73.32 primes
  suffice. The primes of (X1/2, X1] number ≈ 2^77.3 — a 16× slack over N_A.
  The same concentration (factor-2 slack) enters the prefix bound above. H2
  affects only how many primes are needed, not success.

Evidence: `evidence.py` with the declared inputs; one run in 4.1s local. These results hold at toy scale
only.

- **Exact-rule coverage ladder** (`coverage_prime_sweep`): the four prime-construction points of the
  table above, plus at d = 2^4 the random-set baseline (0.4219), the smallest-integers-in-H baseline
  (0.3952) and the matched integer-list control (0.3296, deficit +0.0949).
- **Hash probe** (`hash_probe`): 16,384 canonical-y keys from the structured production-style walk
  (primes in increasing order, prime class lists walked in increasing order at one addition per
  element, on a 27.8-bit curve with 2 ∉ H — the production walk precondition), inserted into 16
  shards at load 1/2: **1.486** mean probes per insert and **2.512** per absent lookup against the
  1.5/2.5 Knuth bounds, shard loads within **2.55 σ** of Binomial.
- **End-to-end toy attack** on a 23.85-bit j = 0 curve with d = 2^6 and GMAX = 16, X2/d = 64, under
  the production block rule (ka = 41), so that the overflow path runs (270,020 chunk additions over
  40 trials; 0 overflow doublings — the b_prev = GMAX case did not occur). The toy attack runs wrapped
  (X1·X2 > r_toy), unlike production (2^82·2^70 = 2^152 < r): its success figures are stated as
  measured.
  - 40 shuffled trials: 22 successes (**0.55**) against 0.4267 predicted, 1.48 σ high, at 7,492.7
    walk additions per trial.
  - 16 fixed-order trials under the deterministic rule exactly as declared: 7 successes (**0.4375**),
    0.11 σ high.
  - Every recovered scalar was correct in both runs.

## 4. Cost (unit: one batched incomplete affine G1 addition = 2,000 instructions)

Per-element budgets, in instructions on top of the reference addition.

| Item | Budget | Itemised count (this package's own ledger) |
|---|---:|---|
| canonical y | 30 | p − y (≈ 7), 6-limb compare (≈ 12), select (6): ≈ 25 |
| hash insert (baby) / lookup (giant) | 240 | endpoint 61 + 6 relay stages × 29 (18 word moves + 4 destination selection + 4 queue updates + 3 control); evidenced by the hash_probe toy check |
| walk bookkeeping | 100 | read b and the gap (4), overflow test (2), address T[gap] (3), load T[gap] (12), load and store points (24), batch-inversion buffers beyond the 6M (24), loop (5): ≈ 75 |

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
radix-1024 switch fabric over the same ≤ 2^60 nodes needs ceil(60/10) = 6 hops, at the identical
per-stage relay cost (18 instructions load and store each of the 9 record words). Lookups follow the same
route, and a match reply is needed only on success. These are expected costs summed over 2^117
operations, where linear probing and balanced sharding concentrate sharply. The χ-sieve is segmented
over ranges of b, so each segment is local, and appends to the prime class lists are streamed. Walks
read one list sequentially.

This gives c_B = c_G = 1 + 370/2000 = 1.185. The full ledger, recomputed exactly by `evidence.py`
(every charge above its itemised count):

| Term | log₂ (units) |
|---|---:|
| baby walk N·c_B | 115.95 |
| giant walk N·c_G | 115.95 |
| gap-overflow additions (deterministic: Σ gaps ≤ X2 per prime, ≤ N_A·X2/GMAX) | 109.32 |
| per-prime tables T[k] (GMAX additions × N_A primes) | 107.32 |
| per-prime scalar multiplications (3 × 1,000 units per prime) | 84.87 |
| χ for primes ≤ 2^82 (100 units = 2·10⁵ instructions each: exponentiation plus Pohlig–Hellman) | 83.148 |
| prime sieve to X1 (64 instructions per integer); χ-sieve and prime class lists to X2 | 77.035 |
| collision checks (6 candidates) | 14.136 |
| d-th root (≤ 2^20 F_r mults of 200 instructions) | 16.678 |
| μ_d search (2^20 scalar mults) + final check | 31.551 |
| input decoding and subgroup checks | 13.551 |
| **total W(2^20)** | **116.958777** |

- A scalar multiplication is charged at 1,000 units for a ≤ 82-bit scalar and 3,000 units for a 255-bit
  scalar. These are generous Jacobian ladders at ≤ 5 units per step, plus one inversion at 100 units.
- The score is rounded up: 116.958777 = 116958777/1000000. Every cost term scales as √K; raising K
  from 0.52 to 0.55 costs 0.040 and buys the 12.3× allowance cushion of section 3.

**Memory.**
- Hash table: 2N slots × 72 bytes ≈ 2^122.9 bytes.
- Prime class lists: π(2^70) ≈ 2^65.03 primes × 9 bytes ≈ 2^68.7 bytes.
- Walk state and affine tables for 2^16 primes in flight: 2^16 · 2^34 · 96 bytes ≈ 2^56.6 bytes.
- Memory is not scored. Any BSGS form of Cheon stores one side, about 2^115.7 entries.
- Rho-type Cheon variants need little memory, but each of their steps is a multiplication in the
  exponent (≥ 7 additions with KKM tables).

## 5. W(n) and applicability

The algorithm uses only A_1 and A_{2^20}, so W(n) = 2^116.958777 for every 2^20 ≤ n < r. Nothing is
claimed for n < 2^20. More generally, for any d | r−1 with 6 | (r−1)/d and d ≤ n:
- W ≈ 2·c·√(K(r−1)/(6d)), plus the lower-order terms above;
- X2 = d·2^50.

## 6. Context (unscored estimates)

| Method | log₂ W |
|---|---:|
| Pollard rho on G1, order-6 classes: 0.5·log₂(πr/12) additions, plus overhead | ≥ 126.5 |
| Cheon BSGS, double-and-add scalar multiplication per step (≈ 383 additions) | ≈ 125.2 |
| Cheon BSGS with KKM c = 3 tables (2 additions plus one F_r multiplication per step); same K and per-element budgets as above, computed by `evidence.py` | 117.897 |
| **this work (v25)** | **116.959** |

- The gain over KKM-accelerated Cheon (recomputed under the same budgets) is a constant factor, about
  0.94 bit.
- **Relation to the generic bound.** Boneh–Boyen's generic lower bound Ω(√(ε·r/d)) has unstated
  constants. Here the walks use 2N ≈ 2^116.7 additions, and the generic coverage is
  N²/m = K/6 ≈ 0.092. The factor 6 that lifts success to 1 − e^{−K} comes from the cheap action of U,
  through negation and the GLV endomorphism. The model does not capture that action, so there is no
  conflict with the bound.
- **Context, not a premise.** KKM (Pairing 2007, paraphrased) set aside a table of small-integer
  multiples as hard to use, because such integers are scattered among the exponents. The product
  construction a·b, with χ(a) + χ(b) ≡ 0, is one way to use them. We found no prior description of
  it, but nothing in the claim depends on novelty.

**Disclosed heuristics:**
- H1, with the integer-level families eliminated or bounded, the modular residual measured with no
  unexplained excess at any point, and the 12.3× allowance cushion;
- H2;
- expected-cost hashing and sharding.

No other assumptions.
