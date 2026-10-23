# GT Frobenius-class walk: shared specification (revision v11)

Used by `curvesmash-core/cscore/gt.py` (toy model, Python), the blst step kernel
(`curvesmash-core/examples/blst-gt-rho/step.c`) and its Python replica. Written 2026-10-22; revisions v7 to v10
are summarised in the package history (proof.md, section 0). Revision v11 (2026-10-23, round 11) changes four things,
each answering a round-10 g2 objection:

1. Restart and escape exponents are outputs of a keyed PRF (section 6), not entries of a fixed table. Restarts use a
   fixed-base comb for g^e. Premise H8 is now the PRF premise (section 9.3).
2. The square/add decision and the adder index come from a keyed mix of the key (section 5), not from the key's
   low bits. The square probability is exactly s = 458/1024. The non-degenerate collision rate V_low is an assumption, not a
   derived law (section 9.3). The square share enters through its expectation under the random oracle (H9, section 9.3).
3. The trail age is not reset by an escape (section 6). The trail cap therefore bounds every trail, including escape
   loops. The trapped-trail bound follows from the cap and uniform landings (H1), not from an in-degree argument.
4. The additive cells are all 12 r_a (not the 12 (r_a - S) of v10), which changes the formal-closure union bound
   (section 9.5).

Any change here must be applied to all three implementations (kernel, replica, toy).

## 1. Group and Frobenius

- G = cyclic group of prime order m (production: m = r, the BLS12-381 subgroup order;
  toy: m = ℓ, a prime dividing Φ12(q) = q^4 − q^2 + 1 with ℓ ∤ 12).
- σ(x) = x^q, with q = p in production and a small prime in the toy model. λ = q mod m.
  Then σ(g) = g^λ, λ has order 12 and λ^6 = −1, so σ^6(x) = x^{-1}.
- Frobenius conjugates of x: σ^n(x) for n = 0..11. The orbit has 12 elements when x ≠ 1.
- Classes: the orbits. N = (m − 1)/12 classes (the identity excluded).

## 2. Scalar value and key

- For x in GT, V_n(x) is a fixed scalar of σ^n(x) (production: the real part of the
  w^1 coefficient in the blst basis Fp2[w]/(w^6 − ξ), ξ = 1+u; toy: the base-q integer
  encoding of the coefficient vector of σ^n(x) in F_q[x]/f).
- Key: K(x) = min_n V_n(x) and k(x) = argmin_n V_n(x) (smallest n on ties).
- K and the orbit are class invariants: K(σ^j x) = K(x). Ties are exceptional (probability
  about 2^-250 per step in production; they are resolved by smallest n and verified by the
  replica).
- K_lo(x) = the low 64-bit word of the integer K(x) (production: the first limb of the Montgomery
  representation of K that the kernel stores; toy: K mod 2^64). K_lo is a class invariant.
- Production implementation of the 12 values: with a = Re c1 and B = Re(conj(c1)·d1),
  where d1 = γ_{1,1} = ξ^{(p−1)/6}, the 12 values are ±a·ζ^j and ±B·ζ^j for ζ a primitive
  sixth root of unity in Fp, j = 0, 1, 2. Four Fp multiplications (B: two; aζ and Bζ: one
  each; aζ² = aζ − a and Bζ² = Bζ − B need none). The replica checks this formula against
  direct computation of V_n(σ^n(x)) on real GT elements.

## 3. Exponent invariant

- A walk state is (x, a, b) with x = g^{a + b·y} (mod m), where z = g^y is the target and
  g generates G (production: g = e(P, Q0), z = e(P, Q0)^y).
- The state x need NOT be the canonical element. The canonical element is σ^k(x) =
  g^{A + B·y} with (A, B) = (λ^k a, λ^k b), where k = k(x).
- Two states whose canonical elements are equal give y = (A2 − A1)/(B1 − B2) mod m, when
  B1 ≠ B2. Every candidate y is checked by g^y = z before it is output.

## 4. Twisted adder table (no Frobenius map per step)

- r_a = 2^adder_bits = 1024 adders. Adder h has exponents (c_h, d_h), chosen at random, with
  R_h = g^{c_h} z^{d_h}. Adder elements are built by the simulation shortcut in the toy
  model; the walk never reads y.
- For k in 0..11: T[h][k] = σ^{-k}(R_h) = g^{λ^{-k} c_h} z^{λ^{-k} d_h}. Store the pairs
  E_c[h][k] = λ^{-k} c_h mod m and E_d[h][k] = λ^{-k} d_h mod m, with the element.
  Table size 12 · r_a elements (production: 12 · 1024 · 576 bytes ≈ 7 MiB).
- Every adder is used by add steps (section 5); the additive cells (h, k) are all 12 r_a of them.

## 5. One step (mixed add/square walk)

Let S = 458 and r_a = 1024, so s = S / r_a = 458/1024 exactly. Let κ_b be a 64-bit branch key drawn at setup (uniform in production; the
measurements draw it per trial). Let mix be the splitmix64 finaliser (section 6).

1. Compute K = K(x), k = k(x) and K_lo (section 2).
2. u = mix(K_lo XOR κ_b). The step is a square slot iff (u >> 54) < S (the top ten bits of u). For a uniform
   u this has probability exactly s.
3. Square slot: x ← x²; (a, b) ← (2a, 2b). Squaring commutes with σ, so no twist is needed. Production uses
   blst_fp12_cyclotomic_sqr (charged at its instructions).
4. Add slot: h = u mod r_a (the low ten bits of u); x ← x · T[h][k]; a ← a + E_c[h][k]; b ← b + E_d[h][k]
   (mod m). Production uses one blst_fp12_mul (one reference unit).

Correctness of step 4: σ^k(x') = σ^k(x)·R_h, so the class of x' is the class of the canonical element times the
adder, as required. Step 3 is exact. The branch and the index depend only on K_lo and k, which are class invariants,
so the step is a function of the class (this is what Lemma R uses).

## 6. Trails, distinguished points, restarts, escapes and solving

- Trails. A trail begins at a walk's start element, at a DP restart, or at a trail-cap restart. It ends at its first
  distinguished point (DP) or at the trail cap. The trail age is the number of main steps since the trail began; it is
  reset only when a trail begins and is NOT reset by an escape (below). The cap is CAP = 2^66 steps of trail age. The age
  counts the walk's main steps (adder and square steps, escapes' targets included as the steps that follow them); the
  walk-through of a detected cycle (at most 64 points, each DP-tested) is not counted, and it is charged as an event
  (66 adder operations, section 9.3). An escape target continues its trail; it does not start one.
- Restart counter. Each walk has a counter t, incremented before each DP restart and each cap restart (not at escapes).
  The walk's index is w.
- PRF. Let κ_P be a 256-bit key fixed at setup (random in production; fixed in the measurements). For a message m,
  PRF(m) = (H(κ_P ‖ m ‖ 0x00) ‖ H(κ_P ‖ m ‖ 0x01)) mod r, read as a 512-bit big-endian integer, where H is SHA-256.
  Bias of the reduction is below 2^-260. Messages:
  - restart: m = 0x01 ‖ w (8 bytes, big-endian) ‖ t (8 bytes, big-endian);
  - escape: m = 0x02 ‖ K_s (the 48-byte big-endian integer of the least-key state's key; section below).
  PRF outputs for distinct messages are treated as independent and uniform on Z_r (premise H8).
- Comb. Preprocessing builds G[i][v] = g^{v · 2^{8i}} for i = 0..31 and v = 1..255 (8160 Fp12 elements, 8160
  multiplications, about 4.7 MiB). For e < r < 2^256, g^e = Π_i G[i][byte_i(e)] over the nonzero bytes, so g^e costs at
  most 32 multiplications.
- Restart (DP or cap). Increment t; e = PRF(restart message); x ← x · g^e (comb); (a, b) ← (a + e, b) mod r. Clear the
  detector state (previous key, anchor) and set the trail age to 0. The restarted state is itself a DP candidate:
  if it is a DP, its record is stored and the walk restarts again (chained restarts, probability 2^-60 each in
  production). A restart does not advance the anchor schedule, which refreshes every 64 walk steps (rounds, one step per
  walk per round). A walk stores no start element after loading.
- Escape (detected cycle, replay-invariant). On detection of a cycle C (2-cycle or anchor), the walk (i) walks C once
  more (at most 64 steps; the DP test applies to every produced point) and finds its least-key state s (minimal K over the
  members of C; K is class-invariant); (ii) takes the canonical representative c_s of s's class and its canonical
  exponents (A_s, B_s) = (λ^kn a_s, λ^kn b_s), kn the argmin index of s; (iii) sets e = PRF(0x02 ‖ K_s), where K_s is the
  48-byte integer of the key of s, and moves to x = c_s · g^e with exponents (A_s + e, B_s) mod r. The escape clears the
  detector state, does NOT reset the trail age, and does NOT advance t. Its target is a DP candidate (chained, as a
  restart, with the counter incremented for the DP restart only). The escape depends only on the cycle's least-key state
  (through K_s), not on the walk's detector state, previous key, anchor phase, or trail age. It is charged as one escape
  (66 adder operations) plus one restart (section 8).
- Replay lemma (Lemma R, v11). Two walks at the same class (any detector states, exponents or trail ages) produce the
  same element sequence until a DP, a detection or a cap in either walk, because each step is a function of the current
  element (section 5). Let a detection occur: the cycle C it closes contains no DP (the DP test precedes the detector),
  and every point of C has been produced by the detecting walk; the other walk, on the same element sequence, enters C
  at the same point and detects it within 127 steps of that entry, with the same least-key state s and the same escape
  target c_s g^e (e a function of K_s only). So after the detections both walks sit at the same element, with the exponent
  offset of the shared class preserved (the escape adds the same e; squares double both exponents). Consequently the
  first DP on the shared element sequence is reached by both walks and recorded by both as the same element with the same
  offset, provided no cap intervenes in either walk before it (a cap in one walk separates the pair; it is counted in R,
  section 9.4). A pair of such records with B1 ≠ B2 solves for y.
- Escape outputs and closure rate (v11). The escape target c_s g^e is a state of the same walk: it is counted in M (or in
  R when it is a restart, section 9.4). Its class is uniform (H8) for distinct cycles, so it is a fresh class for H7. The
  closure rate is argued per fresh entry: a fresh state's successor closes a cycle (a first revisit of the walk's window)
  with conditional probability at most p (section 9.3), by the adder structure. Revisit states are charged through the
  reserve w. No independence across revisits is claimed; the expected total is covered by the Markov reserve K p M.
- Revisits and loops (v11; replaces the v10 premise H10 and its in-degree argument). Every escape target is a PRF output
  (escape, above), so distinct cycles give independent uniform targets, and no table element or fixed power is reused. The
  only non-uniform landings are the closures of section 9.3 (short formal cycles, each detected and counted in E). Any other
  landing of a trail on a class already visited by that trail is a revisit with probability at most (visited)/N under H1
  (uniform landing on unvisited classes), so at most CAP/N for a trail of age at most CAP. A trail that loops (a self-trap
  through its own escape cycle, or a composite loop through several escape cycles) must contain such a revisit, because each
  escape is followed by a fresh uniform path until the next detected closure, and a loop closes only on a visited class.
  Over the whole run the expected number of own-trail revisits is at most M · CAP / N (section 9.3, a term of delta). A trail
  that never revisits is not trapped; it reaches its DP or its cap, and the cap bounds the steps of any trail. The table-reuse
  loop of round-10 g1 (A, B = (AR)^(-1), C = AY, D = (CR)^(-1) with a shared Y) needs one fixed element used by two escapes; v11
  has no such element, so its composite-loop analogue needs an accidental equality of PRF exponents or a uniform revisit, and
  both are bounded above. The toy checks the configuration directly (section 5 of the proof).
- Joins and duplicates. A join is a first landing of a walk on a state of another trail (its own earlier trail included).
  At a join the walk follows the joined trail to its DP (the joined trail is deterministic), stores one duplicate of
  that record, and restarts with its next counter value. A join is a collision: two records with equal fingerprints and,
  generically, different B. The analysis counts joins by the landing bound of section 9.3 (H1 for steps, H8 for restarts).
- Distinguished point: (K >> adder_bits) & (2^dp_bits − 1) = 0 (production dp_bits = 60). Store (fingerprint, A, B),
  80 bytes: the fingerprint is key bits 70..197 of the Montgomery integer of K (128 bits; bits 10..69 are zero for
  every distinguished point), and A, B are the canonical exponents A = λ^kn a, B = λ^kn b mod r. The record is stored
  before the restart.
- Two-cycle (checked every step, after the DP test of the new point): the new key equals the key two steps back and the
  following step returns to the older key. These 2-cycles are formal (the exponents return exactly; section 9.5).
- Longer cycles up to length 64 (per walk): an anchor key is recorded every 64 of the walk's own steps. When the new key
  equals the anchor key, the walk has entered a cycle of length L ≤ 64 (detected within 127 steps of its entry).
- Solving: records with equal fingerprints are compared pairwise. A pair from the same class with B1 ≠ B2 gives
  y = (A2 − A1)/(B1 − B2) mod r; a degenerate pair (B1 = B2) is discarded and costs one comparison. Every candidate is
  verified by g^y = z before output. No trail is replayed.
- Walks: WALKS = 2^20 walks run in parallel. Each walk is a separate state with its own start element, detector and
  trail age. M is the number of main steps (of all walks together) that the charge bounds cover (section 9).

## 7. Parameters

- Production: r_a = 1024, adder_bits = 10, S = 458 square slots, s = 458/1024 = 0.447265625 exactly (section 5).
  S is chosen by minimising the per-step charge over sqrt(V) (adder and square charges of section 8); the optimum is within
  a few parts in 10^6 of S = 458, which is the measured count.
- Toy: adder_bits, S and the trail cap are inputs of the experiment. The budget of the toy uses the assumed value V_low
  (section 9.3). The toy success statistic is the verified-output rate, not the first-useful rate.

## 8. Measured quantities required by the accounting

- Per-type step costs: the adder step (reference operation blst_fp12_mul, charged 1, plus the other instructions divided
  by C = 20,000) and the square step (blst_fp12_cyclotomic_sqr at kappa, plus the other instructions divided by C). Every
  main step is charged the adder cost (the larger of the two). The step code (v11) is straight-line (the mix, the branch
  comparison, the table index, comparisons over all limbs, branch-free modular addition and minimum), so each type has one
  instruction count, measured with every step of one type (forced inputs). The counts are in measurement/RESULT.json of the
  package.
- Restart (DP or cap): 32 multiplications (the comb, charged at the adder cost each, conservative) plus the PRF
  (two SHA-256 evaluations for the restart message, charged PRF_INST = 20,000 instructions, conservative). Escape: the
  closure charge (66 adder operations) plus one restart. Preprocessing: 8160 multiplications for the comb table.
- Distinguished-point costs: the store of one record (the dp path of the kernel, measured as the difference between the
  every-step-DP and the no-DP runs minus the restart; RESULT.json) and the bitonic sorting network of section 9.2: c_cmp
  instructions per compare-exchange and c_scan per slot (measured in the standalone kernel sortnet.c;
  measurement/sortnet/RESULT-sortnet.json).
- kappa: instructions of one blst_fp12_cyclotomic_sqr, divided by C.
- The operation counts per main step and the cycle-event counts come from uninstrumented runs of the same kernel
  (measurement/native_run.sh) over 2^20 rounds per seed (evidence/native-counts.json).
- Measurement: measurement/reproduce.sh (measure/measure.sh with LIB=blst, N1 = 256, N2 = 512, two seeds), and
  replica.py for the dump, distinguished-point, restart and count checks.

## 9. Procedure, work cutoff, capacity and the distinguished-point table

### 9.1 Work counter and cutoff

The procedure keeps a work counter W in units (charges as in section 8). Each operation adds its charge: main steps
(adder or square, each at the adder charge), restarts (section 8, at DP and at the trail cap), escapes (66 adder operations
plus one restart each), distinguished points (section 9.2 and 9.4), walk starts (2 exponentiations each), verified
fingerprint comparisons (one exponentiation each), the adder table, the comb table, and the final solve (2^29
instructions plus one pairing).

The walks run until W reaches W_MAX − reserve, where reserve is the solve and output charge. Then the procedure stops and
reports failure if no verified y was found. W_MAX is the claimed computation.

### 9.2 Distinguished-point table (buckets, v9)

- Buckets. The table is B = 2^37 independent buckets, each of n = 2^31 slots (64-bit indices), so each bucket fits a 64-bit
  address space. A record goes to bucket b = fp[1] >> 27 (the top 37 bits of the 128-bit fingerprint). Equal fingerprints
  therefore share a bucket, and matching needs no global sort. Total slots 2^68 (padding included).
- Store. A distinguished point appends its 80-byte record to its bucket (a cursor per bucket). If a bucket is full the
  procedure stops and reports failure (overflow rule); if the table holds D records before a store the same rule applies.
- Pad. Before sorting, the unused slots of each bucket are filled with sentinels (fingerprint all ones).
- Sort. Each bucket is sorted independently by a bitonic network on its n = 2^31 slots, ascending by the 128-bit
  fingerprint (fp[1] the high word, unsigned): for k = 2, 4, ..., n; for j = k/2, ..., 1; for i = 0..n−1 with l = i XOR j
  and l > i: compare-exchange(i, l), ascending if (i AND k) = 0. The compare-exchange is branch-free. Each bucket has
  n k (k+1)/4 compare-exchanges, k = 31.
- Match. Equal neighbouring fingerprints in a bucket are enumerated over the equal run, each pair charged as a fingerprint
  match (section 9.3). A pair with B1 = B2 is discarded at one comparison.
- Scratch and memory. A bucket sorts in place; its scratch is one 80-byte temporary. Total storage 2^68 slots x 80 bytes
  (2^74.3 bytes), spread over the buckets.
- Cost. Per record: store (measured), one restart (section 8), bucket scatter (measured), transfer and index widening (an
  explicit assumption, 100,000 instructions per record), and a 10-instruction overflow check. Per slot: the bucket sort
  (c_cmp per compare-exchange times the comparators per slot, k = 31) and the scan (c_scan per slot). The table is charged
  for its full slots (padding included) and its records at capacity D.

### 9.3 Premises, charges and reserves for the first M main steps (accounting.py, build())

Premise H1 (random map). Each landing of a walk (a main step or a restart) is uniform on the classes not yet visited by
that walk, given the history, except for the revisit structure described in section 6 (closures, detected cycles).

Premise H7 (fresh states). Every fresh state's key is uniform given the past. Escape targets are fresh (H8).

Premise H8 (PRF restarts and escapes, v11). The restart and escape exponents are outputs of a PRF with a random key
κ_P (section 6). Distinct messages give independent uniform exponents on Z_r, independent of the history. So a restart or
escape target x · g^e (or c_s · g^e) is uniform on G given the history, and it lands on a visited class with probability at
most (visited)/N. This is the standard PRF premise. It is used only for joins, duplicates and the DP overflow bounds.

Premise H9 (square share, expectation form). The square decision of a class x is the branch bit
b(x) = [(mix(K_lo(x) XOR κ_b) >> (64 − adder_bits)) < S], and the adder index is the low bits of the same word. Random-oracle
assumption: κ_b is a uniform key, and the words mix(K_lo(x) XOR κ_b) of distinct classes are independent uniform 64-bit words.
Then (a) a fresh step is a step TAKEN FROM a class the walk has never visited (its branch bit is decided at that class, and the
freshness of the class is fixed before the decision, so the bit is unrevealed); its square decision is Bernoulli(s) with s = S/r_a exactly,
conditional on the history, because the history reveals only the branch bits of visited classes; hence the expected number of
square decisions among the fresh steps is s times their expected number, with no conditioning. (b) The number of useful
collisions among fresh pairs is taken at the assumed rate V_low (the next paragraph), not derived from the share. (c) The
success model is the H1 Poisson approximation with this expected count. The key-level fraction f_N of square classes is an
average of N independent bits over the random key, non-adaptive (no walk quantity enters it): by Hoeffding,
|f_N − s| ≤ t_N = sqrt(ln(2 / δ_N) / (2N)) with probability at least 1 − δ_N, δ_N = 2^-60, and δ_N is a term of delta.
No conditional concentration (Azuma), no visited-set correction and no margin is used. The measurement of the expectation
draws a fresh κ_b per trial (section 5 of the proof).

Non-degenerate collision rate (assumption, v12). The budget uses V_low = 0.7996551 as the non-degenerate probability of a
fresh collision pair. This value is an assumption. It is the output of the earlier expression 1 − s² − (1 − s)²/r_a at s = 458/1024, kept as an assumed conservative lower bound; the expression is
withdrawn, because two square steps from distinct classes cannot collide (squaring is injective on the odd-order group),
and degenerate pairs are formal two-cycles. At N = 76880 the toy measures a degenerate fraction of 0.002 to 0.012 at
r_a = 1024 (mixed) and 0.02 (add-only), so the measured non-degenerate rate near 0.98 lies above V_low. That margin is
evidence at toy size, not a bound at production size (proof section 2).

Unconditional means. Let p = event_rate_bound (section 9.5 and accounting.py) be the unconditional probability per main step
that a closure event starts (a first revisit of a window state of the walk's own trail, lengths 2 to 64). Its expectation over
the first M steps is at most p M (the union over the windows, each closing with the probability of its formal or
non-formal pattern given fresh keys, by the tower property).

- Closure events (Markov). A closure event is a detected cycle. Its entry is a first revisit of a window state of the walk,
  made from a fresh state. Each fresh state's successor is such a revisit with conditional probability at most p (H7, the
  adder structure), and fresh states number at most M, so E[E] ≤ p M and P(E > K p M) ≤ 1/K. Independence across revisits is
  not used. The reserve is K = 50 (failure 1/50 = 0.02, in delta). Each event is charged 66 adder-type operations and one
  restart (section 8).
- Revisit steps. Each closure event is followed by at most 127 revisit steps before its detection (the anchor fires within
  127 steps of entry). These steps are main steps. So the revisit steps are at most 127 E ≤ 127 K p M, and the reserve
  fraction is w = 127 K p (+ the join suffix term and the trapped-trail term below).
- Own-trail revisits (H1, section 6). The expected number of revisits of a trail to its own earlier classes over the whole run
  is at most M · CAP / N (each landing at trail age a ≤ CAP revisits with probability at most a/N). Such a revisit is treated as
  a trap, its steps are at most CAP each, and the probability that any occurs is a term of delta (section 9.7). The reserve
  includes (M CAP / N) · CAP steps, a negligible fraction of M.
- Joins (Markov). A join is a first landing of a walk on a state of another trail (its own earlier trail included). Every new
  state, a main step or a restart (a DP restart, an escape or a cap reset), lands on a visited class with probability at most
  (visited)/N: H1 for main steps, H8 for restarts. With R the number of restarts, R ≤ E_reserve + M 2^-60 + M/CAP, the
  expectation is at most (M + R)² / (2N); since R ≤ M, this is at most M²/N (the charged bound). With the allowance
  J_MAX = 1000, P(J > J_MAX) ≤ E[J]/J_MAX. A join is followed by at most CAP revisit steps, so the join suffixes are at most
  J_MAX CAP steps, a fraction of M of order 2^-50.
- Detection delay. A collision is detected when both walks reach a DP. Each walk's residual to its DP after the collision is
  geometric with mean 2^60 (under H7, DP flags of fresh keys have probability 2^-60). The allowance is T_DELAY = 40 mean DP
  waits per walk, so lag = WALKS · 40 · 2^60 steps, and the residual exceeds the allowance with probability at most e^-40.
  The collision rate is multiplied by (1 − e^-40) in the success model.
- Fresh steps. The fresh steps are M minus the revisit steps; M is chosen so that M (1 − w) − J_MAX CAP ≥ T + lag. The success
  model needs T fresh steps, and F ≥ T holds on the reserve event.
- Restarts: one restart each (section 8). Cap resets are at most M/CAP. The comb table costs 8160 multiplications.
- Fingerprint matches: metered (one exponentiation each); the metered total is bounded by 2^60 times the expected count
  (Markov, failure 2^-60).
- Non-formal closures of length at most 64: expectation 12 · 64 · M / r (Markov).
- Trap events (undetected cycles longer than 64): the trap term of section 9.5.

### 9.4 Distinguished-point capacity

New states are the main steps (M) and the restarts (R: DP restarts, escapes and cap resets). Each new state is a fresh DP
candidate with flag probability 2^-60 (under H7 and H8), so the fresh flags X are dominated by a binomial with n = M + R
trials and mean mu = (M + R) 2^-60 (Chernoff: P(X > 2 mu) ≤ exp(−mu/3)). Records are fresh DP flags X plus duplicates, one
per join, so records ≤ X + J pathwise. The capacity is D = ceil(2 mu) + J_MAX, with R ≤ E_reserve + M 2^-60 + M/CAP ≤ M, so
mu ≤ (M + M) 2^-60 at most. Overflow needs X > 2 mu (failure exp(−mu/3)), or J > J_MAX (failure E[J]/J_MAX, section 9.3),
or a bucket above 2^31 records. Bucket occupancy is binomial with mean D/2^37 against capacity 2^31; the Chernoff bound per
bucket and the union over 2^37 buckets give a negligible failure (accounting.py, bucket_failure_log2). The procedure stops
when a bucket or the total capacity is exceeded (9.2), so only these events bind.

### 9.5 Formal cancellations and window lengths

A closed window is a window of steps whose composition returns the walk to an earlier class or state. Its exponent coefficient
is the sum of the step contributions. If the window contains a square, the coefficient carries a factor 2^q with q ≥ 1 that is
not a 12th root of unity (section 9.6), so it is nonzero: such closures are non-formal and bounded by the 1/r term (H5). If the
window is all additive, the coefficient of adder h is sum_k n_{h,k} λ^(−k) in F_r, where n_{h,k} counts the steps of adder h with
rotation k. Lemma (norm). For windows of length L < 2^63, this is zero in F_r only if sum_k n_{h,k} ζ^(−k) = 0 in Z[ζ12].
Proof: the image of α = sum n_{h,k} ζ^(−k) in F_r is zero, so α lies in a prime 𝔭 of Z[ζ12] above r (r = 1 mod 12, so r splits),
and a nonzero α in 𝔭 has norm divisible by r with |N(α)| ≤ L^4 < 2^252 < r, so α = 0.

Proposition G (vanishing sums of 12th roots of unity). Let n in Z_{≥0}^12 with sum_k n_k ζ^k = 0 for a primitive 12th root of
unity ζ. Then n is a nonnegative integer combination of the indicator vectors of the 2-gons {k, k+6} and the 3-gons {k, k+4, k+8}
(indices mod 12). The proof is in the package's proof.md, section 2, and is stated here without repetition.

Per adder, the closed window's rotation multiset is a nonnegative solution, so Proposition G gives its decomposition. The union
bound of accounting.py (union_log2, lengths 3 to 64 and 65 to 400) counts the gon decompositions of each length, with ADDITIVE_CELLS
= 12 r_a (v11: every adder is an additive cell; v10 used 12 (r_a − S)). For windows longer than 400 the probability that every step
is additive, (1 − s)^401/s, is used. The 2-cycle term (1 − s)/(12 r_a) is the exact 2-cycle rate: the first step additive, the
return step with the same adder and rotation +6, which returns the exact element and exponents (T[h][k+6] = T[h][k]^{-1}).

### 9.6 Order of 2 and the smallest q

r − 1 = 2^32 · 3 · 11 · 19 · 10177 · 125527 · 859267 · 906349^2 · 2508409 · 2529403 · 52437899 · 254760293^2 (cscore/number_theory.py:
the product, primality of each factor by Miller–Rabin with the first twelve primes as bases, deterministic below 3.18·10^23, and r
as a probable prime with the same bases). The multiplicative order of 2 modulo r is o = 3277242198445386904965483781761622864855659531282977363912728668746161324032,
251 bits, certified by 2^o = 1 and 2^(o/p) ≠ 1 for each prime divisor p of o. Since 2^q is a 12th root of unity iff o | 12q iff q0 | q,
with q0 = o/gcd(o,12) (gcd(q0, 12/gcd(o,12)) = 1), the smallest q ≥ 1 is q0, a 248-bit number larger than 2^247. Hence 2^q is not a
12th root of unity for 1 ≤ q < 2^247.

### 9.7 Failure terms

Each term is computed in accounting.py (build()) and summed into the failure total subtracted from the success probability:
the trap term (9.5), the own-trail revisit term (H1, section 6: M · CAP / N), the fingerprint Markov term 2^-60, the
key-level Hoeffding term for the square share (H9, 9.3: delta_N = 2^-60), the closure-event Markov term 1/K (9.3), the non-formal closure expectation,
the join Markov term (9.3), the fresh-flag overflow exp(−mu/3) and bucket overflow (9.4), and the detection delay e^-40 (9.3).
Success is then the model with the collision factor (1 − e^-40), with V evaluated at s + ε_s (section 9.3), minus the total.
