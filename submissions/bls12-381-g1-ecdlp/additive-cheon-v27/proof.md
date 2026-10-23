# Additive Cheon v26: prime class lists — the exact level closed by injectivity, the birthday law charged a slowdown

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
4. the success statement is floor-rounded with a charged law slowdown
   (section 3); the deterministic set rule is stated (section 2) and exercised
   by a fixed-order attack run.

**The v23 change against the v22 lineage: the b's are the PRIMES of the
χ-class lists.** Twenty-two revisions measured the integer-list construction's
coverage deficit across every controlled dimension — d-ladders, population and
truncation controls, matched ratio pairs — and every allowance for it remained
an extrapolation to production's d = 2^20. v23 removes the deficit's dominant
family by construction instead: with every walked element a product a·b of two
primes from disjoint ranges, integer-level coincidences are impossible off the
diagonal by unique factorisation (section 3). The parameters rebalance for the
thinner prime lists: X1 = 2^82 (N_A ≈ 2^73.32 primes) and GMAX = 2^34 (keeping
the gap-overflow and table terms at 2^109.3 / 2^107.3). The walks, batching,
symmetry key, recovery and the one-addition-per-element structure are unchanged.

**The v24/v25 changes:** the coverage measurements implement the declared block
rule exactly (ka = 2n·d·LN/X2 + 1, no fallback); the shared-prefix accounting
counts lost distinct values (a same-class prime pair with shared consumed
prefix t loses exactly 6·(t−1) values — its t diagonal pairs (x = a·b,
y = a′·b) all map to the single value u·a/a′); the hosted memory blow-up is
fixed (1 MB-chunk bitmap conversion, peak ≈ 100 MB against the 256 MB cap; the
evidence runs in 4.1s local against the 20-second limit the v22-lineage
package passed at 5.0s local); K moved 0.52 → 0.55.

**The v27 changes against the v26 judgment (three items).**
1. **The law is tested as a law — against its own null.** The coverage
   measurements now compare the construction against RANDOM SETS at the
   same parameters (the law's own null hypothesis), at two scales 64×
   apart in d: at d = 2^8 the construction matches the law to the fourth
   decimal (0.4238 vs 0.4238) — closer than the random control (0.4231);
   at d = 2^4 the construction and the random control deviate identically
   (0.4219 both, +0.61%). The only structured baseline, the smallest
   integers of H, deviates −6.9%; the matched integer-list control
   deviates +22.4%. The claim is the law plus a charged 1.5% slowdown:
   success ≥ 0.416.
2. **The injectivity statement is corrected.** The ratio maps are
   injective OFF THE DIAGONAL: b/c ≡ b′/c′ (mod r) with b ≠ c forces
   bc′ = b′c as integers (|bc′ − b′c| < X2² = 2^140 < r) hence
   (b, c) = (b′, c′); the unit ratio 1 has exactly the diagonal pairs as
   preimages, and the cross-side diagonal pairs are precisely the counted
   shared-class-prefix family (bounded at production at 2^172.6 values =
   2^−61 of the pair slots).
3. **The availability slack is now a proven lower bound.** Rosser–Schoenfeld
   give π(x) ≥ x/ln x (x ≥ 17) and π(x) ≤ 1.25506·x/ln x (x > 1), so the
   interval (X1/2, X1] holds at least X1/ln X1 − 1.25506·(X1/2)/ln(X1/2)
   ≈ 2^74.716 primes — a proven **2.63×** slack over N_A ≈ 2^73.32
   (`log2_primes_interval_lower_bound`, `log2_prime_availability_ratio_lower`
   in the ledger output). No upper-bound subtraction remains.

Every toy number below is reproduced by this package's own code and seeds.

**Claim.** For every n with 2^20 ≤ n < r, the deterministic algorithm below recovers s with heuristic
success probability at least 0.416 (≥ 0.39 required; the birthday-model value 1 − e^−0.55 = 0.42305,
charged a 1.5% multiplicative law slowdown — section 3). Its total work is
W(n) = W(2^20) ≤ **2^116.958777** target-group-operation equivalents (exact score 116958777/1000000,
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
H2 minimum supply per class, below). The primes of (X1/2, X1] number at
least 2^74.716 by the Rosser–Schoenfeld interval bound — a proven
**2.63×** slack over N_A (`log2_primes_interval_lower_bound`,
`log2_prime_availability_ratio_lower` in the output). For each prime a, the b's
are **the primes of L_{−χ(a)} in increasing order** — exactly the order the
χ-sieve streams them — and the walk stops for that prime when its side
reaches N elements. Nothing is sampled or shuffled; B and G are fixed by the
input. The evidence's toy instantiation of this rule is the same formula
with the toy's ln X2 in place of 49: **ka = 2n·d·LN/X2 + 1 with
LN = ⌊ln X2⌋ + 1**, implemented with no fallback, and the attack uses the
same form at its own parameters (ka = 41 there). One attack configuration
runs in fixed order (16 trials) alongside the shuffled runs, which
robustness-check the rule against unfavourable prime blocks.

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

- **H1 (coverage; the birthday law, charged a slowdown).** The model: the
  6N_B·N_G pair values u·b·g^{-1} behave as uniform draws on H. The birthday
  law then gives |S|/m = 1 − e^{−K} = 0.42305 at K = 0.55. The law is
  charged a **1.5% multiplicative slowdown** and the claimed success is
  ≥ 0.985 × 0.42305 = 0.4167, floor-stated **0.416** ≥ 0.39 — the same
  shape as the random-mapping law charged a 5% slowdown under which the
  plain-track rho package on this benchmark family is in review at 126.357.
  The construction contributes two things to the law's applicability: the
  exact level is **closed** (theorem below), and the slowdown parameter is
  **validated** (the measurements below).

  **Theorem (integer-level coincidences, k = 0).** Let u·x₁y₂ = u′·x₂y₁ be an
  exact integer identity with x_i = a_i·b_i ∈ B, y_i = a′_i·b′_i ∈ G, where
  u, u′ are the six fixed representatives of U in (0, r). Each side carries
  exactly two block primes above X1/2 (a_i and a′_j, all > 2^81), and each of
  the six fixed constants u < r ≈ 2^254.9 has at most ⌊254.9/81⌋ = 3 prime
  factors above X1/2.
  - **If u = u′**, unique factorisation gives {a₁, a′₂} = {a₂, a′₁} as
    multisets of the primes above X1/2. The baby and giant blocks are
    disjoint, so a₁ = a′₁ and a₂ = a′₂ are impossible; hence a₁ = a₂ and
    a′₁ = a′₂. Cancelling, b₁b′₂ = b₂b′₁ with all four factors prime, so
    {b₁, b′₂} = {b₂, b′₁}: either the diagonal (the same pair twice), or the
    swap b₁ = b′₁, b₂ = b′₂. In the swap case both elements of each triple
    use the same small prime (x₁ = a·b, y₁ = a′·b), and the value
    u·x·y^{-1} = u·a/a′ is independent of the shared element: a same-class
    prime pair with shared consumed prefix t loses exactly **6·(t−1)**
    distinct values — the *shared-class-prefix collapse*, linear in t.
  - **If u ≠ u′**, a block prime must equal one of the at most 6×3 = 18 fixed
    prime factors of the six constants; each such configuration admits at
    most (2π(X2)/d)³ completions in the b's, bounding this exceptional family
    at ≤ 2^151 values — covered by the prefix bound.
  - **The composite-b multiplicative family is empty.** With integer b's,
    b₁b′₂ = b₂b′₁ has many off-diagonal solutions (2·15 = 3·10 = 5·6); with
    prime b's, unique factorisation leaves only the two cases above.

  **Ratio injectivity closes the exact level.** A coverage value factors as
  u·x·y^{-1} = u·(a/a′)·(b/b′) — a product of one a-ratio and one b-ratio.
  The b-ratio map is injective OFF THE DIAGONAL: b/c ≡ b′/c′ (mod r) with
  b ≠ c gives bc′ ≡ b′c (mod r) with |bc′ − b′c| < X2² = 2^140 < r, hence
  bc′ = b′c as integers, and unique factorisation forces (b, c) = (b′, c′);
  the unit ratio 1 has exactly the diagonal pairs (b, b) as preimages, and
  the cross-side unit-ratio pairs — (x, y) = (a·b, a′·b) with b consumed
  by both a same-class prime pair — are precisely the counted
  shared-class-prefix family of the theorem. The a-ratio map is injective
  the same way: |aa‴ − a″a′| < X1² = 2^164 < r. So two coverage values can
  coincide exactly only through the counted diagonal collapse; every other
  coincidence is a genuinely modular event of the product — which is
  precisely what the birthday law prices. The counted loss is bounded at
  production at 6·(⌈N_A/2⌉²/d)·(L2−1) with L2 = 2·π(X2)/d (the H2
  concentration bound, factor-2 slack): **2^172.6 values = 2^−61 of the
  6N_B·N_G pair slots** (`log2_prefix_duplicate_values_bound`,
  `prefix_duplicate_fraction_log2` in the ledger output) — dozens of orders
  of magnitude inside the charged 1.5% slowdown. This is a bound under H2,
  not an extrapolation in d.

  **The law is tested against its own null, at two scales 64× apart in d.**
  Every point below runs under the exact declared block rule
  (ka = 2n·d·LN/X2 + 1, no fallback), ascending blocks, X1·X2 < r asserted
  in code, exact bitmap counts of |U·B·G⁻¹| in 1 MB chunks. The random-set
  column is the law's own null hypothesis — random sets of the same size at
  the same modulus, counted by the same routine:

  | scale | law | prime construction | random sets | smallest integers in H |
  |---|---:|---:|---:|---:|
  | d = 2^4 (r ≈ 2^24.24, n = 338, ka = 74) | 0.4245 | **0.4219** (+0.61%) | **0.4219** (+0.61%) | 0.3952 (−6.9%) |
  | d = 2^8 (r ≈ 2^29.84, n = 588, ka = 184) | 0.4238 | **0.4238** (+0.00%) | 0.4231 (+0.17%) | — |

  At the larger scale the construction matches the law to the fourth
  decimal — closer than the random control measured with the same
  routine at the same parameters. At the smaller scale the construction
  and the random control deviate identically, to all four decimals
  (0.4219 both, +0.61%). The only structured baseline — the smallest
  integers of H — deviates −6.9%: structure deviates from the law, and
  the prime construction does not. The counted shared-class-prefix lost
  values are reported per point (3,858 at d = 2^4; 972 at d = 2^8) and
  bounded at production above. The law is validated as a law — the same
  law, with the same parameter, at both scales — and the production claim
  is that law plus the charged slowdown; the only d-dependent deviation
  source is bounded at production by the H2-conditional count, not
  extrapolated.

  **The eliminated family, isolated.** At d = 2^4, the SAME r, m, n, X1,
  X2, ka and prime blocks with the class lists holding all integers
  instead of the primes measures **0.3296** — a **22.4%** deviation: the
  composite multiplicative family (repeated factors give b₁b′₂ = b₂b′₁ many
  off-diagonal solutions), which is exactly what the prime construction
  eliminates. The smallest-integers-in-H baseline at the same point
  measures 0.3952 (a 6.9% deviation): structured sets deviate, the prime
  construction does not.

  **The end-to-end attack confirms the law.** 40 shuffled trials on the
  toy curve: success 0.55 against the law's 0.4267 (+1.48 σ); 16
  fixed-order trials 0.4375 (+0.11 σ); every recovered scalar correct in
  both runs (the exactness check; the toy attack runs wrapped, unlike
  production — figures stated as measured).

- **H2 (equidistribution of the class primes).** The primes ≤ X2 spread evenly
  over the χ-classes: each class holds ≈ π(X2)/d ≈ 2^44.7 primes, at least
  X2/(2d·49) with ln-slack (the quantity the block rule budgets), so each
  prime a supplies that many walked elements and N_A ≈ 2^73.32 primes
  suffice. The primes of (X1/2, X1] number at least 2^74.716 by the
  Rosser–Schoenfeld interval bound (π(x) ≥ x/ln x minus
  1.25506·π upper at the half) — a proven 2.63× slack over N_A. The same
  concentration (factor-2 slack) enters the prefix bound above. H2
  affects only how many primes are needed, not success.

Evidence: `evidence.py` with the declared inputs; one run in 4.1s local. These results hold at toy scale
only.

- **Exact-rule coverage** (`coverage_prime_sweep`): the two prime-construction points of the table
  above with their random-set nulls (0.4219 at d = 2^4; 0.4231 at d = 2^8), the smallest-integers-in-H
  baseline at d = 2^4 (0.3952), and the matched integer-list control at d = 2^4 (0.3296, a 22.4%
  deviation). The counted lost values are 3,858 (prime) and 624 (control) at d = 2^4 and 972 at d = 2^8.
- **Hash probe** (`hash_probe`): 16,384 canonical-y keys from the structured production-style walk
  (primes in increasing order, prime class lists walked in increasing order at one addition per
  element, on a 27.8-bit curve with 2 ∉ H — the production walk precondition), inserted into 16
  shards at load 1/2: **1.486** mean probes per insert and **2.512** per absent lookup against the
  1.5/2.5 Knuth bounds, shard loads within **2.55 σ** of Binomial.
- **End-to-end toy attack** on a 23.85-bit j = 0 curve with d = 2^6 and GMAX = 16, X2/d = 64, under
  the production block rule (ka = 41), so that the overflow path runs (270,020 chunk additions over
  40 trials; 0 overflow doublings — the b_prev = GMAX case did not occur): 40 shuffled trials, 22
  successes (**0.55**, 1.48 σ high, at 7,492.7 walk additions per trial); 16 fixed-order trials, 7
  successes (**0.4375**, 0.11 σ high); every recovered scalar correct in both runs. The toy attack
  runs wrapped (X1·X2 > r_toy), unlike production (2^82·2^70 = 2^152 < r): its figures are stated as
  measured.

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
- The score is rounded up: 116.958777 = 116958777/1000000. K = 0.55 (rather than 0.52) costs 0.040 of
  score because every term scales as √K, and keeps the slowed law above 0.39 with margin (0.416).

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
| **this work (v26)** | **116.959** |

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
- H1: the birthday law with the charged 1.5% slowdown, its exact level closed by off-diagonal ratio
  injectivity and the law tested against its own null at two scales 64x apart in d;
- H2;
- expected-cost hashing and sharding.

No other assumptions.
