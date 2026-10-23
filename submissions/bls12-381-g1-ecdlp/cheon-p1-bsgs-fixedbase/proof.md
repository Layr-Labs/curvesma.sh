# bls12-381-g1-ecdlp: Cheon p-1 attack with fixed-base comb tables and a mu_6 class reduction

Status: ready for independent exploratory review; no qualification or emitted
score is claimed by this file. Costs are in `classical-group-work-v2` units:
one batched incomplete affine BLS12-381 G1 addition (6M, S=M, A=0) costs 1, and
any other work costs (scalar AArch64 instructions)/C with C=2000.

## Notation and number-theoretic facts (checked against the protected profile)

- r = 0x73eda753...00000001 is the G1 order, log2 r = 254.857.
- r-1 = 2^32 * 3 * 11 * 19 * 10177 * 125527 * 859267 * 906349^2 * 2508409 *
  2529403 * 52437899 * 254760293^2 (the product equals r-1 exactly).
- g = 7 generates F_r^*: 7^((r-1)/q) != 1 mod r for every prime q above.
- d = 2^20 divides r-1, and A_d = [s^d]P0 is supplied at the scoring n = 2^20.
- m = (r-1)/d, log2 m = 234.857. 6 | m. zeta = g^d generates the unique
  subgroup H < F_r^* of order m.
- mu_6 (the sixth roots of unity in F_r) lies in H, because 6 | m and F_r^* is
  cyclic. -1 = zeta^(m/2), and the GLV endomorphism phi(x,y) = (beta*x, y),
  with beta a primitive cube root of unity in F_p, acts on G1 as multiplication
  by a primitive cube root of unity lambda in F_r. So for every Q in G1 the set
  {[w]Q : w in mu_6} = {+-Q, +-phi(Q), +-phi^2(Q)} is exactly the set of points
  whose x-coordinate lies in {x, beta*x, beta^2*x}.
- Canonical class key: kappa(Q) = min(x, beta*x, beta^2*x), compared as
  canonical integers in [0,p). Since 1+beta+beta^2 = 0, beta^2*x = -x - beta*x,
  so kappa costs one F_p multiplication and two modular subtractions.

## Success argument

s is uniform on F_r^*, the d-th power map F_r^* -> H is exactly d-to-1, and
H is cyclic of order m. So s^d = zeta^k with k uniform on [0,m).

Phase 1 (recover s^d). Set M6 = m/6, choose integers u, G with u*G <= M6 and
6*u*G >= 0.39*m. Baby points are B_j = [zeta^j]P0 for 0 <= j < u. Giant points
are C_i = [zeta^(-i*u)]A_d = [s^d * zeta^(-i*u)]P0 for 0 <= i < G. Every
k in S = {f*M6 + i*u + j : 0<=f<6, 0<=i<G, 0<=j<u} yields
s^d zeta^(-iu) = zeta^(f*M6) * zeta^j with zeta^(f*M6) in mu_6, hence
kappa(C_i) = kappa(B_j). The 6uG elements of S are distinct because uG <= M6,
so Pr[k in S] = 6uG/m >= 0.39 exactly. (No heuristic enters this.)

Conversely, kappa(C_i) = kappa(B_j) means C_i = [w]B_j for some w in mu_6,
i.e. s^d = w * zeta^(iu+j). There are no false positives: the table stores the
full 381-bit canonical key. The six candidates w * zeta^(iu+j) are each checked
by computing [candidate]P0 and comparing it with A_d, and exactly one passes.
Baby keys are pairwise distinct: kappa(B_j) = kappa(B_j') would force
zeta^(j-j') in mu_6, i.e. j = j' mod M6, which is impossible for j, j' < u <= M6.

Phase 2 (recover s from y = s^d). Compute one d-th root s0 of y in F_r with a
Tonelli-Shanks/Adleman-Manders-Miller style 2-power root extraction. The
2-adicity is 32 >= 20, and y in H is guaranteed to be a d-th power. Then
s = s0 * rho^t for a primitive 2^20-th root rho = g^((r-1)/2^20) and a unique
t in [0, 2^20). Test [s0 rho^t]P0 = T for all t; exactly one passes and gives s.
Return it. Phase 2 is deterministic given phase-1 success, so the overall
success probability is >= 0.39 over the uniform challenge, with no algorithm
coins and no amplification needed.

Chosen sizes: u = G = 57012617742321962895928165416718183 (log2 = 115.4568).
These satisfy u*G >= ceil(0.39*m/6) and u*G <= m/6, and 6uG/m = 0.39000...
(>= 0.39).

## Fixed-base evaluation: two reference additions per point

Every baby or giant point is [e]Base for a known nonzero scalar e in [1,r) and
one of two fixed bases, P0 or A_d. Write e = a0 + 2^85 a1 + 2^170 a2 with
0 <= a_w < 2^85 (since r < 2^255). For each base, precompute three comb tables
T_w[v] = [v * 2^(85w)]Base, v = 1..2^85-1, w = 0,1,2, stored affinely.
Then [e]Base = T_0[a0] + T_1[a1] + T_2[a2], skipping zero digits:
at most two group additions per point.

The scalars are maintained incrementally in F_r: e_{j+1} = zeta * e_j for the
baby points and e'_{i+1} = zeta^(-u) * e'_i for the giant points. This is one
F_r Montgomery multiplication per point. The constants zeta and zeta^(-u) are
computed once by square-and-multiply.

Points are processed in batches of 2^20 independent scalars. Layer 1 forms
T_0[a0]+T_1[a1] for the whole batch with one Montgomery-trick inversion, and
layer 2 adds T_2[a2]. This is exactly the batched incomplete affine addition
of the reference operation.

Exceptional cases are handled exactly:
- Layer 1 can never hit a doubling or an inverse pair. Both digits are below
  2^85, so a0 +- 2^85 a1 lies in (-2^170, 2^170) and is nonzero mod r unless
  both digits are zero. Zero digits are skipped (infinity is never fed to the
  formula).
- Layer 2 needs a0 + 2^85 a1 = +-2^170 a2 mod r. The minus sign means e = 0,
  which is impossible. The plus sign is a doubling. Before batching, each pair
  compares x-coordinates (6 words). Equal-x pairs are removed from the batch
  and finished with a separate affine doubling, charged at <= 1 reference
  addition plus one F_p inversion (<= 2^18 instructions = 2^7 units). The
  equal-x condition pins e to at most 2^85 of the r-1 scalars per base, so its
  frequency is negligible. Even if all phase-1 points hit it, the analysis
  only needs the comparison cost, which is charged below; the fallback's own
  cost is covered by the per-batch slack (one fallback per 2^20 points would
  add only 2^-13 per point).

## Per-point ledger (baby insertion or giant lookup)

| Item | Reference adds | Instructions (upper bound) |
| --- | ---: | ---: |
| Two batched incomplete affine G1 additions (6M convention) | 2 | - |
| Batch-state loads/stores and Montgomery-trick bookkeeping beyond the 6M arithmetic, 2 layers | - | 200 |
| Per-batch F_p inversion amortised over 2^20 points, 2 layers | - | 1 |
| Exceptional-case x comparison, 2 layers | - | 40 |
| F_r scalar update (4-limb Montgomery multiplication) | - | 300 |
| Digit extraction (3 x 85 bits), multiword table addressing | - | 60 |
| Table loads (3 affine points, 36 words) and zero-digit tests | - | 60 |
| kappa: one F_p mul (383 measured), two modular subtractions, two 6-word compares and select | - | 500 |
| Hash insert or lookup (linear probing, load factor <= 1/2, 8-word records, two-word addresses) | - | 300 |
| Loop counters (two-word i, j) and control | - | 100 |
| **Total** | **2** | **1561 <= 2000** |

So each point costs at most 2 + 2000/C = 3 units.

Hash table. There are 2u slots of 8 words (6-word canonical key, 2-word j).
The slot is taken from key bits not used for comparison. Linear probing at
load <= 1/2 needs at most 2.5 expected probes for an unsuccessful search and
1.5 for an insertion (Knuth). The 300-instruction budget covers 3 probes of
<= 60 instructions each, plus hashing, addressing and the 8-word store.
Exceeding 3 probes on average over ~2^115 operations would need a gross
deviation from the mean, given the heuristic below.

## Total work W(n)

Applicability: n >= 2^20 with d = 2^20 (the attack reads only A_1 = T and
A_d). For the scoring n = 2^20:

- Phase 1 points: 3*(u+G) = 2^118.041799 (worst case, where all G giant
  steps run).
- Comb tables: 2 bases x 3 tables x (2^85 - 1) entries. Each entry costs one
  batched addition, <= 200 instructions of store/addressing, and doubling
  chains of 2 x 170 doublings for the window bases: < 2^88 units. This is the
  stated preprocessing, charged once, not amortised.
- Phase-1 match resolution (six candidates, each one F_r exponentiation
  < 2^19 instructions and one 2-add fixed-base check), phase 2 (root
  extraction < 2^24 instructions, plus 2^20 candidates at <= 3 units each),
  setup (decoding and validating only T and A_d, constants zeta, zeta^(-u),
  rho, beta) and output/verification ([s]P0 = T): <= 3d + 2^26 units.

W(2^20) <= 3(u+G) + 2^88 + 3*2^20 + 2^26, so
log2 W(2^20) <= 118.041799 (rounded up at the sixth decimal). For general
d = 2^e with 20 <= e <= 32 and n >= d, the same construction gives
W(n) <= 6*ceil(sqrt(0.39*(r-1)/(6d))) + 2^88 + 3d + 2^26.

## Memory and advice

- Comb tables: 6 * 2^85 affine points * 96 bytes = 2^94.17 bytes.
- Hash table: 2u * 64 bytes = 2^122.46 bytes. This dominates.
- Peak memory <= 2^122.457 bytes. No nonuniform advice: every table is built
  at run time from P0 and A_d and charged above.

## Heuristics and limitations

- The only heuristic is that hash-slot bits of the canonical keys of baby and
  giant points behave like uniform random bits for linear-probing cost. This
  is the standard assumption for BSGS on elliptic curves. It affects cost,
  not correctness: matches are exact on full keys.
- Success probability 0.39 follows exactly from the uniform distribution of
  k. It needs no heuristic.
- Exponential memory (2^122.5 bytes) and multiword addressing are part of
  the analytic model; memory carries no score. This is not a practical attack.
- This is the classical Cheon (Eurocrypt 2006) p-1 attack. The improvements
  over the textbook version are accounting choices: fixed-base comb tables
  bring scalar multiplication down to two batched additions, and the mu_6
  symmetry (negation plus the GLV endomorphism, both inside H) shrinks the
  search by 6. No novelty is claimed beyond that.
