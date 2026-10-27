# bls12-381-g1-ecdlp: Cheon with twin-output steps on a three-term-progression walk list (twin-ap-v1)

**Claimed time_log2 = 116.538356** (exact rational `29134589/250000` = `116538356/1000000`; an upper bound on log2 W(2^20)). Heuristic success 0.39649 (K = 0.505), claimed 0.396 after every failure term (required 0.39). Peak memory 2^123.3242 bytes, every node counted at its full itemised allocation (exact sum 13313099851307310120424180719332360192 bytes = 2^123.324186..., asserted by the script as an integer inequality). No advice. Preprocessing 2^115.5427 units, included in W. `python3 -I evidence.py` recomputes everything: W, the record counts and memory in exact integers, the probability tails and K in floating point, and a 60-digit decimal check that log2 W does not exceed the claimed score (mode `ledger`). Mode `ec` is the toy run; mode `h2` gives the statistics of the list rule and of the character pairs on the real r.

Dependencies (earlier Yukon submissions on this benchmark, declared with the submit command):
- 0a45ddb2, the promoted entry at 116.793080. Its routing, bucket matching, caps, candidate handling, re-walk recovery, status tree and cluster model are all reused here with the same design.
- Through it, cbadde77: one batched affine addition per step over small elements of H, the mu6 y-key, six candidates per match.
- Through it, 6ecec8e3 / cf6d45a8: multiplier-times-list products, data-independent routing, capped data-dependent items.

**What is new: every step produces two useful points for one inversion.** In the dependencies each record costs one batched affine addition. Here each step computes both P + G and P - G from the same lane inverse 1/(x_G - x_P).
- P + G is exactly the reference operation: a batched incomplete affine addition, priced 1 unit = 2000 instructions.
- P - G adds 2 multiplications, 1 squaring and 5 additions/subtractions, listed at 1427 instructions using the measured 383 per multiplication.

Both outputs are useful because the walk list is built so that u - g, u and u + g are all in H at every step, a three-term arithmetic progression inside H. All chains share the list, so the 2^-40 density of such progressions costs only a one-off scan of at most 2^100 integers (2^107 units). The H-prime factorisation argument is replaced by a short-interval lemma: every list element lies in [2^210, 2^210 + 2^100 + 2^45), which keeps the ratios pairwise distinct as fractions.

Each twin step is charged 2000 + 1700 instructions for two records, i.e. 1850 per record against 2180 in 0a45ddb2. Everything after a record exists (routing, matching, caps, recovery) is the dependency's design, with the numbers re-derived for this N.

Cost accounting. Work is the sum of instructions executed over all processors (COST_MODEL.md: "Charge total work across processors, including preprocessing, advice construction, random sampling, all trials and failures, sorting, lookup, memory transfers, verification, output and success amplification"; rules.md: "Sum work across processors; do not amortize"). The cluster model is the one in 0a45ddb2, stated here so that every charge has a data path:
- Nodes have private memory and no remote-memory instruction. Owner and storage nodes hold 2^40 bytes, generating nodes 2^51 bytes.
- The only communication is posting a 512-byte block into another node's receive FIFO. Filling the block and the post are charged to the sender, the poll and the read to the receiver.
- A node waiting for a block executes no instructions and is charged nothing.

This is my model, not a rule of the challenge. Wall time is not claimed.

## 1. Problem and reduction

Given T = [s]P0 and A_i = [s^i]P0 for i = 1..n, n >= 2^20, recover s. Let d = 2^20. Since v_2(r-1) = 32, d | r-1 and 6d | r-1. H = {x^d : x in F_r^*} has order m = (r-1)/d, log2 m = 234.857089. z = s^d is uniform in H and A_d = [z]P0 is an input. Only T, A_1 and A_d are read, so the same W applies to every n >= 2^20. An integer c is "in H" when c mod r is in H, i.e. c^m = 1 mod r (the H-test).

Lemma 1 (recovery of s from z). Twenty successive square roots in F_r (Tonelli-Shanks) give one t with t^d = z. Then s = t * zeta for one of the 2^20 d-th roots of unity zeta = g^k, where g is a fixed primitive 2^20-th root. Write k = i + 2^10 j. Compute the 2^10 points [t g^i]P0 and the 2^10 points [g^(-2^10 j)]A_1 (2^11 scalar multiplications, 2^23 units), sort the first list by x and look up each point of the second in it. Charged at 2^24 units.

Lemma 2 (mu6 in H). The sixth roots of unity lie in H because 6d | r-1 (asserted by the script). On this j = 0 curve [omega](x, y) = (beta x, y) and [-1](x, y) = (x, -y), so the six points [u]X, u in mu6, share the unordered pair {y, p-y}.

## 2. The twin step

Let P = (x1, y1) and G = (x2, y2) be affine points with x1 != x2. Let i = 1/(x2 - x1) be the lane's inverse from a Montgomery batch over 2^16 independent lanes. The reference operation computes P + G = (x3, y3):
- lambda = (y2 - y1) i
- x3 = lambda^2 - x1 - x2
- y3 = lambda (x1 - x3) - y1

This is the batched incomplete affine addition of the cost model (3M batch share + 3M). It is charged 1 unit = 2000 instructions and nothing else.

P - G = P + (x2, -y2) has slope (-y2 - y1)/(x2 - x1) = -mu with mu = (y1 + y2) i, so

  x' = mu^2 - x1 - x2,  y' = mu (x' - x1) - y1.

Cost on top of the reference, step by step:

| operation | cost |
| --- | --- |
| t = y1 + y2 | 1 add |
| mu = t i | 1 mul |
| mu^2 | 1 sqr |
| x' = mu^2 - x1 - x2 | 2 sub |
| x' - x1 | 1 sub |
| mu (x' - x1) | 1 mul |
| y' | 1 sub |

Total: 2M + 1S + 5 add/sub. The listing prices are:
- 383 instructions per multiplication or squaring. This is the pinned blst_fp_mul measurement in COST_MODEL.md; squaring is charged as a multiplication.
- 48 per 384-bit modular addition or subtraction: 6 loads, 6 add/adc or sub/sbc, the carry, 6 conditional selects, 3 stores and call overhead. This is an estimate, not a measurement.
- 32 for the argument setup of 8 calls, and 6 for keeping the lane inverse (3 stp + 3 ldp).

Listed total: 1149 + 240 + 32 + 6 = 1427 instructions. Sensitivity: if an add/sub cost 64 instead of 48, the side arithmetic would be 1507 and the per-step baby listing in Section 6 1747, 47 above the 1700 charged; the score would rise by log2(3747/3700) = 0.018.

The only quantity shared with the reference operation is the lane inverse i, which the batch computes anyway. No internal instruction of the reference is charged twice, and none is omitted from the extra.

Sanity check against the convention: in field multiplications a twin step is 9M + 5A, against 12M for two separate batched additions (ratio 0.75 to 0.8). The ledger charges 3427 listed (3700 charged) against 4000, a ratio of 0.86 to 0.93, so the 2000-unit convention is not being stretched.

The toy asserts that the twin outputs equal [a(u + g)]Q and [a(u - g)]Q computed by independent scalar multiplications.

Lemma 3 (no degenerate addition). In a chain, P = [a u]Q and G = [a g]Q with 2^210 <= u < 2^210 + 2^101 and 1 <= g <= 2^44. P = +-G would need u = +-g mod r. Both u - g and u + g lie in (0, r), so this is impossible. Hence x1 != x2 at every twin step and neither output is the point at infinity. Table strands add [a]Q to [k a]Q with 2 <= k <= 2^44 (k != +-1 mod r), so they never degenerate either.

## 3. The twin walk list and the multipliers

U is the set of integers in H. X0 = 2^210, T = 2^44.

**List rule.** The list is built once, does not depend on the challenge, and is charged in full.
- **Start.** u_0 is the least element of U at or above X0 + T. It is found by the same scan; under H2 it is below X0 + 2T, since the probability that a 2^44 window holds no element of U is e^(-2^24).
- **Twin step.** Given u_j, scan c = u_j + 1, u_j + 2, ..., u_j + T and take the first c with all three of: c in U; v = 2 u_j - c in U; v not equal to a list element lying in the window [u_j - T, u_j). Then the gap is g_j = c - u_j, the main output is u_{j+1} = c, and the side output is v_j = u_j - g_j = 2 u_j - c.
- **Fallback step.** If no c in (u_j, u_j + T] qualifies, u_{j+1} is the least element of U in (u_j, u_j + T], with no side output. The scan's U-bits for (u_j, u_j + T] are kept in a T-bit bitmap, so nothing is retested. The ledger nevertheless charges a full rescan of T integers for each allowed fallback.
- **Storage and termination.** Per step the list stores the gap (44 bits) and a twin/fallback flag, 6 bytes. It also stores the absolute value u_j at the 2^16 segment boundaries (prefix sums, negligible). Construction runs on the storage nodes and stops, reporting fail, if the scan passes 2^100 integers or the list has more than 2^40 fallback steps.

S denotes the set of all outputs: u_1, u_2, ... and every side output v_j.

Properties enforced by the rule:
- The main outputs strictly increase, and every new main is above every earlier element.
- Every side output lies in [u_j - T, u_j). An earlier element below u_j - T cannot equal it, so the window check makes it differ from all earlier elements. Hence the elements of S are pairwise distinct.
- Every side satisfies v_j >= u_0 - T >= X0.
- Every element is below u_0 + 2^100 < X0 + 2^100 + 2T. So S lies in [X0, X0 + Delta) with Delta = 2^100 + 2^45, and Delta^2 < 2^210 = X0 (asserted).

**Construction cost.** Under H2 (Section 8) a candidate c qualifies with probability 2^-40 (more precisely, at most a factor 1 + 2^-20 less, see H2), so twin gaps are geometric with mean 2^40.
- A fallback has probability (1 - 2^-40)^(2^44) <= e^-16 = 2^-23.08 per step. That gives 2^35.56 expected fallbacks in the whole list against the reserve of 2^40.
- The L = 2^58.6432 twin steps scan about L 2^40 = 2^98.64 integers. The cap of 2^100 is 2.56 times that mean.

The list costs one H-test per scanned integer, plus one more per c in U together with a window duplicate check (at most 2^81). Each H-test is an exponentiation c^m mod r charged 2^18 instructions (235-bit exponent, about 291 4-limb Montgomery multiplications at no more than 383 instructions each, 2^16.8). Charged: (2^100 + 2^44 + 2^81 + 2^84) x 2^18 instructions = 2^107.03 units.

**This is not a probabilistic event.** The list does not depend on the challenge and uses no coins. Whether it fits under its caps is therefore a fixed property of r, not an event over challenges or coins. If H2 were false in a way that pushed the span past 2.56 times its predicted mean, or the fallbacks past 2^4.4 times theirs, the attack would fail on every challenge, not with a small probability. The tail numbers in Section 5 are heuristic estimates under H2, not probabilities over the challenge distribution. The same holds for the H-prime lists of the dependencies.

**Multipliers.** There are 2^56 baby multipliers and 2^56 disjoint giant multipliers. They are H-primes in (2^211, 2^211 + 2^86), all larger than every element of S. The range holds 2^58.81 H-primes in expectation (one prime in 2^20 is an H-prime, density 1/(2^20 ln 2^211)). They are found by an H-test on every integer of the range and a 16-base Miller-Rabin test on the 2^66 survivors (bounded): 2^93.03 units.

**Lemma 4 (distinct ratios, short interval).** Baby records are the points [a s]P0 and giant records are [a' s']A_d = [a' s' z]P0, with a in A, a' in A', and s, s' in S.

Suppose a s / (a' s') = b t / (b' t') as rationals, i.e. a s b' t' = b t a' s'.

1. The prime a divides the right side. It does not divide t or s', which are positive and smaller than a. Also a != a', since A and A' are disjoint. So a = b, and likewise a' = b'.
2. Then s t' = t s'. Write s = X0 + alpha, t' = X0 + delta, t = X0 + beta and s' = X0 + gamma, with alpha, beta, gamma, delta in [0, Delta).
3. Expanding gives X0 (alpha + delta - beta - gamma) = beta gamma - alpha delta. The right side has absolute value below Delta^2 < X0, so alpha + delta = beta + gamma and alpha delta = beta gamma, i.e. {alpha, delta} = {beta, gamma}.
4. So either (s, s') = (t, t'), the same pair, or s = s' and t = t', the diagonal.

So the N^2 ratios are pairwise distinct as fractions, except that for each multiplier pair (a, a') the 2L diagonal pairs share the ratio a/a'. That loses a fraction 1/(2L) = 2^-59.6 of the pairs.

**Lemma 4 is about rationals, not about residues mod r.** The products a s are about 2^421 > r. That the distinct fractions behave like independent uniform elements of H modulo r is heuristic H1, as in the dependencies. One quick check within one multiplier pair: a nontrivial mod-r coincidence needs r | (s t' - t s'). There are about 2^238 quadruples and the differences are about 2^310, so under any reasonable model the expected count of such coincidences per pair is negligible next to the 2^119 distinct ratios. This is a heuristic estimate, not a proof.

## 4. Algorithm B

**Chains.** Chain a holds the points [a u_j]Q, with Q = P0 on the baby side and Q = A_d on the giant side.
- At each twin step it emits records for [a v_j]Q and [a u_{j+1}]Q. The step is X_{j+1} = X_j + G_a[g_j], and the side is X_j - G_a[g_j].
- G_a[g] = [a g]Q for 1 <= g <= 2^44 is the chain's gap table: 96-byte affine entries, 2^50.58 bytes in all.
- The table is built by repeated addition of [a]Q in 2^16 strands, each started by one scalar multiplication. Strand 0 starts at [2a]Q, so every strand addition is [k a]Q + [a]Q with k >= 2. One extra scalar multiplication per chain computes the base [a]Q.
- The chain runs as 2^16 lockstep segments, each started by one scalar multiplication [a u_j mod r]Q at its stored first position. So each round is 2^16 independent twin steps sharing one Montgomery inversion.
- A fallback step is a plain reference addition with one record.

Each chain has L = 450167994536394810 = 2^58.6432 twin steps and emits 2L records. Per side there are 2^56 chains and N = 2 L 2^56 = 2^115.6432 records.

**Generating nodes.** There are 2^20 generating nodes. Each holds one gap table at a time (2^51 bytes) and runs 2^36 baby chains, then 2^36 giant chains. The gap list is streamed to them from 2^30 storage nodes in 512-byte blocks. At 6 bytes per step this is charged 4 instructions per step for transport, covering the sender copy, the post and the receiver poll amortised over 85 steps per block.

**Key.** The record for X = (x, y) is (key, id): key = bits 1..192 of the even one of y and p - y, and id = the 64-bit chain index. Records are 32 bytes. By Lemma 2 the six points [u]X share the key. Main and side outputs are keyed alike.

**Order, routing, matching, caps, candidates and recovery** follow the design of 0a45ddb2 exactly, with this N.
- **Routing.** The whole baby side runs first. There are 2^84 owner nodes indexed by key bits 0..83, reached by three data-independent hops of 28 bits. Records are scattered into 512-byte blocks of 16 records and relayed at (t,0,0) and (t,u,1).
- **Baby insert.** At the owner, baby records go into a direct-mapped head array of 2^33 slots (key bits 84..116), chained through a pool of 2^33 32-byte records (capacity 2.56 x the mean of 2^31.643). A slot holds at most 64 baby records; a 65th makes the node fail.
- **Giant lookup.** Giant records read their head slot and compare their 75 remaining key bits with each chained baby record. Each comparison costs 8 instructions and counts against a per-node counter preset to 1.01 x the per-node mean.
- **Candidates.** An equal key appends a 40-byte candidate (full key, both chain ids) to a 2^16-slot list. After the giant side the status tree totals the candidates; above 2^44 the run fails.
- **Forwarding.** Each candidate is copied to both generating nodes named by the top 20 bits of its chain ids, each with a 2^24-slot receive area.
- **Re-walk.** Each generating node re-walks its own chain once, from start to end: the gap table is rebuilt, the strands, segments and base are restarted, and every twin step is recomputed with both outputs keyed and compared with the retained key. Each re-walk step is charged as a full twin step. The node keeps the outputs whose key equals the retained key, at most 4 per chain; a fifth makes the run fail.
- **Testing.** The giant chain's node sends its kept values a' s' in a 168-byte block to the baby chain's node. That node tests the six z = a s / (u a' s') for at most 16 pairs (at most 96 scalar multiplications) against A_d, and on a hit runs Lemma 1.

Nothing is compared against a whole other chain, and a wrong candidate never produces a wrong answer. Status collection is a fan-in-256 tree over all nodes, run in three passes.

## 5. Sizes, caps and success

L is the least integer with 6 N^2 >= 0.505 m for N = 2 L 2^56. This gives N = 2^115.6432 and K = 6 N^2 / m = 0.505000. With the diagonal loss, K' = K (1 - 1/(2L)). Under H1 the probability that z is among the covered classes is 1 - exp(-K') = 0.39649.

The caps below have the same structure and derivations as in 0a45ddb2, recomputed for this N. Every tail is rounded towards the larger probability, and the charged work is the maximum under the caps.

- **Comparisons.**
  - 2^117 head slots at load lambda = N / 2^117 = 0.390. Mean comparisons are N^2 / 2^117 = 2^114.286, i.e. 2^30.286 per node, with a node cap of 1.01 x that (2^30.301).
  - A slot over 64 records: Poisson tail 2^-273.7 over all slots.
  - A node's baby or giant count more than 2^-10 off its mean: 2^-1453 over all nodes.
  - Given these, Bernstein (M = 64, variance margin 4) for a 0.8 % excess: 2^-10462 over all nodes.
- **Matches.**
  - False matches (equal 192-bit key, different mu6 orbit) have mean N^2 / 2^192 = 2^39.29 against a cap of 2^44; Chebyshev tail 2^-47.6.
  - True matches are Poisson with mean K under H1, with a budget of 2^20; exceeding it has probability below 2^-20492289.
  - Candidates per owner node have mean 2^-44.7 against 2^16 slots; overflow below 2^-2930260.
  - Candidates per generating node have mean 2^20.29 against 2^24 slots; overflow below 2^-39943108.
- **Re-walk key repeats.** One key on more than 4 of a chain's 2L outputs: below 2^-419 over all chains.
- **Node pool overflow.** Chernoff 2^-(3.3 x 10^9) per node; a budget of 2^-20 covers modelling error.
- **List construction.** These are heuristic estimates under H2, not probabilities over challenges (see Section 3). Span over 2^100: Chernoff bound for a sum of L exponentials exceeding 2.56 times its mean, below 2^-4e17. Fallback steps over 2^40: below 2^-3.37e12.

Success claimed: 0.39649 - 2^-20 - 2^-40 > 0.396 > 0.39 (asserted). The 2^-40 covers every cap and re-walk tail above, each of which is below 2^-47. The 2^-20 covers node-overflow modelling.

## 6. Instruction listing per twin step (AArch64 scalar, on top of the one 2000-instruction reference addition)

**Walk overhead, once per step (48 instructions):**

| item | instructions |
| --- | ---: |
| gap stream read and pointer advance | 5 |
| gap stream transport (6 bytes per step, blocked) | 4 |
| table address | 2 |
| operand loads and stores not counted inside the reference price | 18 |
| batch-inversion prefix and suffix loads and stores | 6 |
| lane control and loop | 8 |
| flow control and counters | 5 |

**Side output arithmetic, once per step (1427 instructions):** see Section 2.

**Record handling, once per record and twice per step.** Baby records cost 96 each:

| item | instructions |
| --- | ---: |
| parity select and key extraction | 23 |
| scatter at 3 hops (12 per hop) | 36 |
| post and poll of full blocks, amortised over 16 records, at 3 hops | 9 |
| 32-byte transport copy at 3 hops | 18 |
| insert | 10 |

Giant records cost 92 each: the same items without the insert, plus a lookup fixed cost of 6. Comparisons are charged separately, at 8 each, up to the node cap.

**Totals.** Listed per twin step: baby 48 + 1427 + 2 x 96 = 1667, giant 48 + 1427 + 2 x 92 = 1659. Charged: 1700 per twin step on both sides, which is 1.02 x the larger listing (33 instructions of margin for unlisted items). That makes 3700 instructions per twin step and 1850 per record. Gap-table additions are charged 2000 + 48, and re-walk steps are charged as full twin steps (3700).

## 7. Ledger (units of 2000 instructions; W = ceil of the instruction sum / 2000)

| term | log2 units |
| --- | ---: |
| records: 2^57 chains x L twin steps x (2000 + 1700), two records per step | 116.531 |
| per-node fixed work: (2^84 + 2^30) x 2^35 + 2^20 x 2^48 instructions | 108.034 |
| twin list S: (2^100 + 2^44 + 2^81 + 2^84) H-tests at 2^18 instructions, window duplicate checks included | 107.034 |
| comparisons: 2^84 nodes x cap x 8 | 106.335 |
| candidates at cap 2^44: re-walk both chains (L + 2^40 steps at 3700), rebuild both tables, 4 x 2^16 + 4 restarts | 104.531 |
| gap tables: 2^57 chains x 2^44 additions at 2000 + 48 | 101.034 |
| fallback reserve: 2^40 steps per chain x 2^57 chains at 3700 | 97.888 |
| multipliers: 2^86 H-tests + 2^66 Miller-Rabin tests | 93.034 |
| status tree: 3 passes x every node x 2^12 instructions | 86.619 |
| table strand starts: 2^57 x (2^16 + 1) scalar multiplications at 2^12 units | 85.000 |
| segment starts: 2^57 x 2^16 scalar multiplications | 85.000 |
| chain bases [a]Q: 2^57 scalar multiplications | 69.000 |
| candidates at cap: up to 6 x 4 x 4 = 96 scalar multiplications each | 62.585 |
| Lemma 1 up to 2^20 + 1 times, input access, output | 44.000 |
| candidates at cap: kept values to the verification node, 256 instructions | 41.034 |
| candidates at cap: copy 40 B to two generating nodes, 128 instructions | 40.034 |
| **W** | **116.5383553** |

W = 120653688096167569697692259864438375 units, log2 W = 116.538355337040093... (the 60-digit decimal check is in the script). Score = ceil(10^6 log2 W) / 10^6 = 116.538356.

For comparison with the dependency: its record row is 2N x (2000 + 180) at N = 2^115.6644, with K = 0.52. The twin step moves the per-record price from 2180 to 1850 instructions, which is log2(1850/2180) = -0.237. K = 0.505 accounts for another -0.021.

Preprocessing is the work that does not depend on the challenge: the baby side, the lists, the baby tables and node setup. It comes to 2^115.5427 units and is included in W.

Memory, with every node at its full itemised allocation, summed as integers:

| nodes | bytes each | itemised allocation |
| --- | --- | --- |
| 2^84 owner nodes | 2^39.3242 | heads 2^36, pool 2^38, send blocks 2^37, three receive FIFOs 3 x 2^36, candidate list 2^16 x 40 B, code and counters 2^30 |
| 2^20 generating nodes | 2^50.5851 | gap table 2^44 x 96 B, send blocks 2^37, stream buffers 2^36, candidate receive area 2^24 x 40 B; allocation 2^51 |
| 2^30 storage nodes | 2^40 | list 2^58.64 x 6 B plus multipliers 2^57 x 32 B, 2^62.67 B used |

The total is 2^123.324186 B. The claim is 2^123.3242 and the script asserts mem^10000 <= 2^1233242. The records alone are 2^120.64 B.

## 8. Heuristics

**H1 (coverage).** The 6 N^2 (1 - 1/(2L)) distinct fractions a s / (u a' s') behave like independent uniform elements of H for the purpose of hitting z. This is the same heuristic the dependencies use, applied to a different but equally structured set: products of a large prime and an element of a short interval. Lemma 4 removes the exact repeats. The toy runs the identical construction at d = 2^7 and agrees with 1 - exp(-K') within 0.7 standard errors at 400 trials (0.3625 against 0.3762, z = -0.56) and at 2000 trials (0.383, z = +0.63). The toy cannot separate 0.396 from 0.39; the production figure rests on H1.

**H2 (twin density).** For an element u of the list and nearby integers c, the events "c in H" and "2u - c in H" behave like independent events of probability 2^-20 each. This is a statement about the 2^20-th power residue character on pairs of nearby integers. It is of the same kind as H2 of the dependencies (H-prime density and exponential gaps).

The list rule itself adds a small bias. A side 2u_j - c lies below u_j, in integers that an earlier step already scanned and rejected as candidates. Rejection lowers their U-probability from q = 2^-20 to at most q/(1 + q). The mean gap therefore exceeds 1/q^2 by at most a factor 1 + q: 1 + 2^-20 in production, but 1.125 at level 2^3 and 1.0625 at 2^4.

Mode `h2` gives two kinds of evidence on the real r:

(a) **The exact list rule** applied to the integers 2^210 + i, i < 24000, with H replaced by the 2^k-th powers modulo r:

| level | mean gap / 4^k | prediction (1 + 2^-k) | twin steps | coefficient of variation |
| --- | ---: | ---: | ---: | ---: |
| 2^3 | 1.228 (standard error 0.061) | 1.125 | 281 | 1.03 |
| 2^4 | 1.054 (standard error 0.122) | 1.0625 | 59 | 0.94 |

Geometric gaps have a coefficient of variation of about 1. There were no fallbacks.

(b) **Direct sampling** of u uniform in [2^210, 2^210 + 2^100) and g uniform in [1, 2^44], with the low 8 bits of the 20-bit character logarithms of u, u + g and u - g computed exactly (Pohlig-Hellman in mu_{2^8}). With 4000 samples (hosted input):
- chi-square of (e+ mod 16, e- mod 16) over 256 cells: 263.3 (expected 255, sd 22.6);
- chi-square of (e0, e+, e-) mod 4 over 64 cells: 50.9 (expected 63, sd 11.2);
- P(e+ = e- = 0 mod 2^k) / 4^-k = 0.989, 0.960, 0.848, 0.960 for k = 1..4 (counts 989, 240, 53, 15).

A local run with 30000 samples (seed 7) gave chi-squares 258.4 and 65.3 and ratios 0.982, 0.991, 1.082, 1.118 (counts 7367, 1859, 507, 131).

The top 12 bits of the character logs, which matter for the full 2^-40 density at level 2^20, are not sampled. Their independence is the heuristic itself. The toy curve list (d = 2^7) has a mean gap / d^2 of 0.855 over 13 steps.

**H3 (key bits).** The 192-bit keys behave like independent uniform strings across mu6 orbits: occupancies are binomial and cross-orbit equality has probability 2^-192. This is the same as in the dependencies. The toy chi-square of pseudo-node occupancy is 30.79 ± 0.42 at 400 trials (hosted input) and 30.65 ± 0.18 at 2000 trials, against 31.

## 9. Evidence (`evidence.py`, standard library, about 7 s for mode all on one core)

**`ledger`** computes every number above. W, the record counts and memory are exact integers; the tails, K and success are floating point. It asserts d | r-1, 6d | r-1, mu6 in H, Delta^2 < X0, that the multipliers lie above the list, the listed versus charged per-step instructions, every tail, and the node memory sums with the exact memory inequality. It also runs a 60-digit check that log2 W does not exceed the claimed score.

**`h2`** produces the two experiments quoted in H2.

**`ec`** runs Algorithm B end to end on the j = 0 curve y^2 = x^3 + 3 over F_207670579, of prime order r = 207671809 (d = 128, m = 1622436):
- X0 = 2^40 and T = 2^18. The twin list is built by `build_twin_list`, the rule of Section 3 including the window duplicate check: 13 twin steps, no fallback. Delta^2 < X0 and distinctness of S modulo r are asserted.
- The multipliers are H-primes above the list. There are 14 chains per side with 26 records per chain (n = 364, K = 0.490).
- Each side runs in lockstep with batch inversion, every twin step computing P + G and P - G from one shared lane inverse (`batch_twin`). On the first trial, 52 twin outputs are asserted equal to independently computed scalar multiples [a s]P0.
- Then: the 16-bit parity key, 32 pseudo-nodes x 16 heads, the bucket table, candidates (key, id, id'), one linear re-walk of each candidate chain recomputing both outputs of every twin step, at most 4 kept outputs per chain, six z per pair tested against A_d, 7 square roots and the root-of-unity scan.

Differences from production:
- The gap table holds only the 13 gaps of the toy list, built from doublings, instead of all g <= T.
- The caps other than the re-walk cap are not enforced; the counts are compared with the model instead.
- It runs in one process, so nothing is routed.

Results with seed 2026:

| run | success | predicted (diagonal-corrected) | z | comparisons (exp. 258.8) | false matches (exp. 2.022) | true matches |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 400 trials (hosted input) | 0.3625 | 0.3762 | -0.56 | 258.2 | 1.998 | 0.515 |
| 2000 trials (local) | 0.383 | 0.3762 | +0.63 | 258.6 | 1.976 | 0.52 |

There were no re-walk cap failures, and every recovered s was correct.

## 10. Limitations

- **No full-size run.** This is an analytic bound.
- **Side-output price.** The instruction count uses the measured blst multiplication price and an itemised estimate of 48 instructions per modular addition. The charge of 1700 per step is 1.02 x the listing; at 64 per addition the listing would exceed it by 47 instructions per step (+0.018 in score).
- **Pricing rule.** The pricing relies on the cost model's rule that one batched incomplete affine addition costs one unit. The side output shares only that addition's lane inverse and is charged in instructions on top.
- **Heuristics.** H1, H2 and H3 are heuristics validated at toy scale. H2 is also checked on the real r at power levels up to 2^4 and on the low 8 bits of the character logarithm. The list and multiplier constructions are deterministic: a failure of H2 would make the attack fail on every challenge, not with a small probability.
- **Memory.** Memory is exponential (2^123.3242 B). The generating nodes need 2^51 bytes each for the 2^44-entry gap table. A smaller table (T = 2^42) raises the fallback rate to e^-4 per step and the score by about 0.013.
- **Matching stage.** The matching stage and its cluster model are the dependency's, with the same caveats. A reviewer who prices communication differently gets a different number for the routing items: scatter, post and poll, and the transport copy are 63 listed instructions per record, which is 126 of the 1700 charged per twin step. Section 6 has the full list.
