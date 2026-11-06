# Success, section 6 of proof.md (no heuristic)

### 6.1 Model

Coins are independent. Q is uniform. Node g draws `b_g` and `a_j` (j < MX = 256) uniform in [0, r), giving its permanent points `X_j = [b_g]Q + [a_j]P`. Each of its q1 batches k draws `c_k, t_k, s_k` uniform in [1, r) and `rho_k` uniform in [1, p). The batch's L = 2^17 points are `Y_(m,n) = [c_k + m t_k + n s_k]P` (m < 16384, n < 8); incidence l pairs `X_(l mod 256)` with the point of number l >> 8, so each X meets each point once per batch. The running products start at `rho_k`: uniform on F_p^* while no factor is 0. Facts:
- (F1) Given b, the `X_j` are i.i.d. uniform; so all X's are i.i.d. uniform, independent of the batches. Given all points, `b_g` stays independent uniform: for each value of b_g exactly one tuple `(a_j)` fits the points.
- (F2) Given `(t_k, s_k)`, each Y point is uniform on G \ {one point} (`c_k` uniform in [1, r)), independent of the X's and of other batches.
- (F3) For index pairs v != w of one batch, `Y_w - Y_v = [dm t + dn s]P` with `(dm, dn) != (0, 0)`, |dm| < 16384, |dn| < 8, which takes each value with probability <= 1/(r - 1) over `(t, s)`. So `(Y_v, Y_w)` takes each pair of points with probability <= `1/((r-1)(r-1))` (the r/(r-1) factors are inside ρ = (r/(r-1))^4), and given one point of a batch, any other takes each value with probability <= 1/(r - 1) (up to ρ).
- (F4) For three points with affinely independent index pairs, `(c, t, s) -> (Y_u, Y_v, Y_w)` has a nonzero integer determinant of absolute value < r: the triple is uniform up to ρ.

Nodes, hence leaves of different nodes, are independent; a leaf `φ^i(X) + eY` (i = 0, 1, 2; e = ±1) is uniform on G (X uniform, independent of Y).

Γ is the set of unordered cross-node leaf pairs, `|Γ| = (S² - 2^72*n0²)/2` with `S = 2^72*6*q1*q2`, `n0 = 6*q1*q2`. For α = {L, L'} ∈ Γ let `I_α = 1[L' ∈ class(L)]`; then `P(I_α = 1) = P1 = (6(r-1) + 1)/r²` and `μ = |Γ|*P1 = 0.494318018` (`mu`).
Same-node pairs are omitted; their collisions only help. Collisions are single pairs, not clusters (cross-node: four distinct points): a match of `φ^i(X) + eY` with `u(φ^i'(X') + e'Y')` forces a second match among the same incidences' other leaves only through a second independent relation among X, Y, X', Y' (fresh coins), probability <= 1/(r - 1), already inside P1 and 6.2's sums.

### 6.2 Chen–Stein on cross-node pairs

Let `W = Σ_(α∈Γ) I_α` and `B_α` be the pairs sharing an X or a batch with α. I_α depends only on α's four generators (two X points and two batches, from two nodes), which are independent of all other generators (6.1), so I_α is independent of `{I_β : β ∉ B_α}` (b3 = 0). An X lies in `6*q1*L` leaves, a batch in `6*q2`, so `|B_α| <= 12(q1*L + q2)*S`.

For a unit u, the event `L' = u*L` with `L = φ^i(X) + eY`, `L' = φ^(i')(X') + e'Y'` is a linear relation in the eight coins `(x, c, t, s, x', c', t', s')` (x, x' the logarithms of X, X'; Y's is `c + m t + n s`): coefficients `-uλ^i, -ue, -uem, -uen, λ^(i'), e', e'm', e'n'`. I_α <= the sum of these 6 events. Take α != β in Γ sharing a generator, units u, v; E_β is β's relation.
- **Case 1: β has a generator outside α's four.** Its x (or c) enters E_β with a unit coefficient, not E_α; with an x or c of α entering E_α with a unit coefficient it forms a triangular, invertible 2 × 2 minor.
- **Case 2: β has only α's generators.** Restrict both relations to `(x, c, x', c')`, rows `(-uλ^i, -ue, λ^(i'), e')` and `(-vλ^(i2), -ve2, λ^(i2'), e2')`. If they are not proportional, some 2 × 2 minor is invertible. **If proportional,** the factor κ has `e' = κ e2'`, so κ = ±1, and `λ^(i') = κ λ^(i2')`, so κ is a power of λ; -1 is not (powers of λ cube to 1, (-1)³ = -1), so κ = 1, `i2' = i'`, `e2' = e'`; then `uλ^i = vλ^(i2)` and `ue = ve2` give `v/u = e/e2 = ±1 = λ^(i-i2)`, so v = u, `e2 = e`, `i2 = i`. E_β - E_α is `e'(dm' t' + dn' s') - ue(dm t + dn s) = 0`, `(dm, dn, dm', dn')` the index differences of β's points from α's, not all zero (else β = α): probability <= 1/(r - 1) over `(t, s, t', s')`; given these and x, x', c', E_α fixes c: <= 1/(r - 1).

So the relations fix two coins (uniform on >= r - 1 points given the others), and `E[I_α I_β] <= 36/(r-1)²`. Chen-Stein (Arratia, Goldstein and Gordon, Ann. Probab. 17 (1989), Theorem 1) gives `|P(W = 0) - e^-μ| <= b1 + b2 <= |Γ|*12(q1*L + q2)S*(P1² + 36/(r-1)²)`, certificate-checked < 2^-79.

### 6.3 Usefulness

A cross-node colliding pair is **detected** if both records carry their class's label (no batch defect, wrong class
index or label error), both survive the caps, no intruder enters their run, and no candidate, request or relay cap
is reached first. Its run is then resolved, a useful adjacent pair (section 4) gives x, and `[x]P = Q` is checked
before output. The pair is useful iff `n' != u*n`, where `n = λ^i*b_g` and `n' = λ^(i')*b_(g')` for its two nodes; given
all points, this fails with probability exactly 1/r (F1). The resolver tries the records' nine type pairs,
solving with the first relation having `n' - u*n != 0` (2.3): any relation it computes between the points, right or
wrong types, gives the logarithm, so a useful true pair's relation is never lost. So `P(fail) <= P(W = 0) + E #{α
: I_α = 1, α useless or not detected}`.

### 6.4 Defect rows

A **node-local defect** of L (3.4 batch defect, wrong class index, label error, fixup meter) is a function
of node g's coins, given which L' is uniform, so `P(I_α ∧ defect) <= (6/r)*P(defect)`; over Γ a row costs <=
`μ(1 + 1/r)` times its per-leaf rate.

**Caps and intruders.** Condition on α's four leaf points (its two X's and its two Y points), not on whole batches. A **free** generator is an X or a batch other than α's. Every leaf outside the 12 of α's two incidences (`FIXED_LEAVES`) is then uniform on >= r - 1 points: its X is free, or its Y point is of a free batch (F2) or of α's batch (F3). **Label lemma.** In a clean batch (3.4; defective batches anywhere are excluded by `g_def`) a record's fields (d1, d2, K2, d4) are a function of its incidence's fractions, which the six leaves share, but 3.3's error ranges hold for every value of them: for a leaf with exact point Z the emitted fields are the exact label's when Z is outside W and D', and one of <= 12 labels (three paths, two blocks per fraction) inside, all within Z's cell widened by 3.3's ranges. So every event on emitted fields fixed in advance is contained in an event on exact points: a leaf's fields hit a given value of f bits only if Z lies in a widened cell, of mass <= `KB(f)` over any point uniform on >= r - 1 values (3.3's lattice count, whatever index KEY returns). All probabilities below are of exact points; no fraction, representation or index enters them. **Same-incidence pairs.** Two leaves A, B of one incidence are an invertible linear image of (X, Y) (6.2), so given Z_A, Z_B is uniform on >= r - 1 points. If Z_A avoids W and D', A's fields are its exact label's and B's agree with them on f bits only if Z_B lies in that label's widened cell: <= `KB(f)`; Z_A lies in W or D' with probability <= `d_wj + d_lab`, and then A's fields are one of <= 12 values: <= `12 KB(f)`. So the pair agrees on f bits with probability <= `(1 + 11(d_wj + d_lab))KB(f)`, and jointly with I_α: in α's node take it over the node's coins (A, B are an invertible image of their incidence's two generators, and P(I_α | the node's coins) <= 6/r); elsewhere it is independent of α. These 5 pairs per leaf join the candidate, request and relay loads without the home condition (`sip`). **Joint with I_α.** In every case below I_α is a function of α's two leaf points, the pair event (by the lemma) one of A's and B's points, and the named generator is uniform given all the others; each bound is that generator's conditional probability, so it holds jointly with I_α in every overlap configuration listed (i), (ii). The six leaves of an incidence through exactly one of α's points are functions of the other point.
- **Pair slots** count single leaves. A slot has mean <= `m_p = n0*KB(36) + 1` from the free leaves (at hop 1
  one node's leaves of one group with a given d1, n0 the group's; at hop 2 one column's with a given (d1, d2), n0 all
  groups'); each cap leaves room for the <= 12 fixed leaves. Two leaves' indicators have covariance 0 unless they share an X or a batch. A common free X (different batches) gives <= `1/(4(r-1)²)` (<= `KB(f)/(r-1)` with α's Y point). Outside (a)-(c), two points of one batch with the same X or a common free X, or one incidence with both generators free, are uniform up to ρ (F3, F4; an invertible linear image of (X, Y) with Y != O): covariance <= `3/(r²(r-1))` times the two cells' sizes. The ordered pairs of (a)-(c) have covariance <= one leaf's cell mean and number <= `COVX(q) = 60(q*L + MX) + 72 L² + 12 MX*L`, for the two nodes of α:
  (a) the 30 ordered leaf pairs of an incidence through exactly one of α's points: `q*L` through an X per group (q its batches), MX through a Y point;
  (b) two leaves on α's X and batch, functions of `(t, s)` given α's point: `(6L)²` per node;
  (c) a free X with two points of α's batch symmetric about α's, `Y_v + Y_w = 2Y_u`, of one type and opposite signs: `6*MX*L` per node. No other pair through one free X and α's batch is dependent: its two leaves (types i, i2) are an invertible image of (X, c, t, s) unless `P_w - P_u = ±λ^(i2-i)(P_v - P_u)` in the index plane; for i2 != i that makes λ = ±a/b with 0 < |a|, |b| < YM+YN = 16392, which `LAMBDA_FREE` rules out (0 < a² ± ab + b² < r); for i2 = i the sign + gives P_w = P_v.

  So a count's variance is <= its mean plus `1 + COVX(q)*KB(f)`, and Chebyshev bounds a record's drop by `(m + 1 + COVX(q)KB(f))/(c - 12 - m)²`
  for cap c and mean bound m (`d_p1` per group with cap_pg at hop 1, q the group's batches; `d_p2` with cap_p at hop 2, q = q1).
- **Intruders and the candidate, request and relay loads** count leaf pairs, jointly with I_α. For a pair
  π = {A, B} != α with event E_π (keys equal on f bits; at the candidate list and the relay, also A's home or home column
  is α's), `E[I_α*1(E_π)] <= ρ*(6/r)*P*(π)`, where `P*(π)` is E_π's probability for independent uniform leaves and
  `ρ = (r/(r-1))^4` covers the exclusion Y != O and F3's, F4's r/(r-1) factors. (i) If A, say, has a free generator B lacks,
  condition on all else: A is uniform on >= r - 1 points while I_α, B and α's home are fixed, so the f-bit
  match costs <= `(r/(r-1))KB(f)` and the home condition falls on B. If B has a free generator, B is uniform given
  α's points, and its home is α's with probability <= `(r/(r-1))KB(72)` (KB(36) for a home column); if B is a
  partner of L (or L') (a leaf of its node through α's X and a point of α's batch), (L, B) (or (L', B)) is uniform up to ρ over that node's coins, given which `P(I_α) <= 6/r`; if
  B is L or L', π is a pair through L (below). (ii) Otherwise A and B have the same free generators. If their incidence avoids α's points, or they share a free X or batch and avoid α's node pair, condition on α's points: (A, B) is uniform up to ρ (an
  invertible linear image of the incidence's two generators; F3, F4). If both lie in α's node g and use α's X or α's batch there (α's own incidence and (c)'s symmetric pairs included), E_π (with α's home, which on I_α is L's) is a function
  of g's coins, given which L' is uniform, so `P(I_α | g's coins) <= 6/r`; over g's coins (A, B) is uniform up to ρ (F3, or an invertible pair in one incidence), and
  L is uniform given (A, B): through α's X when A and B avoid it (a free X and α's batch: (A, B) is a rank-two image of the free X and two points, and on (c)'s pairs A - B fixes only Y_u), and otherwise by F4, unless A is L (pairs through L, or L',
  are counted without the home condition, both at α's home) or α's point and the two leaves' points are collinear in the index plane: those pairs, at most `72(8 + 2)L` ordered, and the 20 pairs of two partners of L in an incidence of α (10 per incidence), are counted without it too: `PARTNER_PAIRS`. If A is in L's incidence or a partner of L and B likewise for L', `(L, A)` and `(L', B)` are independent and uniform
  up to ρ (F3, or an invertible pair), so each of the six events `L' = uL` meets E_π with probability <= `ρ*P*(π)/r`.

  An intruder stops the later member's walk, so shares its 182 compared bits (d1, d2, K2, d4; section 4): the intruder
  rate is <= `ρ(S + 12)KB(182)` per colliding pair; the candidate list at α's home has mean
  `other_h <= ρ(S²*KB(72)*KB(182) + (2S + PARTNER_PAIRS)*KB(182))` (pairs at the home equal on d1, d2, K2 and d4, the 182 label
  bits a listed pair matches; L and L' are there); a source receives <= `other_g <= 2ρ*n0*S*KB(182)` requests; a
  relay holds <= `other_v <= 18ρ((S/2^36)*S*KB(36)*KB(182) + (2S + PARTNER_PAIRS)*KB(182))` records (a request is one record
  and a reply <= 8), each plus the `sip` pairs. Markov with the caps less α's own entries (4 candidates, 8 requests, 16 records)
  gives `d_cl`, `d_rq`, `d_rl`. Leaves of defective batches elsewhere are not uniform; the whole-run row "a
  computational defect in any batch" excludes them.

| row (per leaf of the pair unless stated) | term in `certificate.py` | checked |
| --- | --- | --- |
| useless collision (per pair) | `1/r` | - |
| batch defect: 3.4's sites (`SITES`), D = 0, ygen's zero denominators (`YGEN_DEF`), equal or opposite points in a batch | `d_batch` = `D_BATCH` | < 2^-127.1 |
| wrong class index (3.3): KEY's window set W | `d_wj` | < 2^-58.47 |
| fixup meter, Markov on the mean `q1*E_FIX` with the cap `2^62 - FIX_BATCH` (28 q2 fixups, 3.4), plus `q1*d_batch` | `d_f` | < 2^-33.5 |
| the label-error set D' (3.3) | `d_lab` | < 2^-16.65 |
| pair slot at hop 1, per group (rounded up to a multiple of 2^-256), summed over groups (Chebyshev, `cap_pg`); at hop 2 (`cap_p`) | `d_p1`; `d_p2` | < 2^-18.58; < 2^-19.9 |
| intruder (per pair), `ρ(S + 12)KB(182)` | `lam_int` | < 2^-52.7 |
| candidate, request and relay caps (Markov; a reply counts 8 records) | `d_cl`; `d_rq`; `d_rl` | < 2^-23.4; < 2^-25.9; < 2^-25.2 |
| a computational defect in any batch (whole run) | `g_def = 2^72*q1*d_batch` | < 2^-26.5 |
| an identity `X_j` (`2^72*MX/r`); a coin whose two range-rejection draws fail (`2^-254`; 4 per batch, MX + 1 per node) (whole run) | in `glob` | all of `glob` < 2^-26.5 |
| RNDR failure (whole run) | `g_rndr` = 0 under the rules' ideal coins (section 7) | - |

`d_batch` is unconditional per batch. A site's true operand is a function of one incidence's points (X, Y) or a running product; for any fixed X, Y is uniform on >= r - 1 points (F2) and at most `mult` (the site table's last column) give one operand value, so a bad set of size b has probability <= `b*mult/(r-1)`; products are uniform on F_p^* by `rho_k` while no factor is 0 (D = 0: six values of Y per incidence). ygen's `YM - 1` start and `YM(YN - 1)` step additions each divide by a difference vanishing for two coin values (1/r each: `YGEN_DEF`); equal or opposite points need `dm t + dn s = 0` or one value of c (`L^2/(r-1)`). A Y point equal to O comes only from an addition with equal x-coordinates (C != O as c_k >= 1), which sets ygen's flag: it is inside `YGEN_DEF`.

### 6.5 The bound

`certificate.py` picks, by bisection from 2^27 (where it fails), the least q1 with

```text
1 - e^-mu - (b1 + b2) - mu(1 + 1/r)(1/r + intruder + 2 * sum of per-leaf rows) - whole-run rows >= 0.39,
```

the per-leaf rows being 6.4's (the hop-1 pair-slot row d_p1 for the 160 groups of 6.6, the hop-2 row d_p2 once each),
with an upper enclosure of e^-μ and exact Fractions: q1 = 145285144, with left side >= 0.3900000000 (ten places, rounded down) and the bound at q1 - 1 below 0.39;
the target is the rules' 0.39 (searched there, as in aburan28's 7516608f).

## 6.6 Grouped early abort: stopping time, hash chains, E-work

Before the run, the q1 batches are split into G = 160 contiguous ranges (groups) of Y-batches, group g ending at
floor(q1*(g/G)^(21/32)) (exact: `group_sizes`; big early groups, small late ones). All nodes run each group in lockstep. Group g generates its batches' Y points (the X's are permanent, set up once) and runs the leaf, route, hash-join and resolve
phases on them exactly as in 6.1-6.5 (same kernels, send path and coins-per-batch laws; fused
hop-1 arrival, proof.md 4; cumulative hop-2 slots of cap_p with per-slot
final cursors F; the cumulative hash table of 2^58 zeroed head words with a packed next-array
and uncapped walks); then resolves group g's new
candidates, stopping the run at the first group resolving to an x with [x]P = Q
(checked, latched). If none resolves, the run ends after group G without output. T
in 1..G+1 is the first resolving group (G+1 if none); group g's leaf/route work and hash fixed cost
are deterministic given g and run iff T > g - 1, while its main-pair walks have means
bounded uniformly in H_g (below), so with P(T > g - 1) <= PHI[g] the expected
work is E_pre + sum_g PHI[g] (W_g + W_K2,g) + W_sameX + W_K2,glob + W_list + W_def + W_rest (worst-case rows, with
candidates and fixups at their expected counts, proof.md 7; W_g includes the key-path mean below, W_K2,g the main pairs' K2 verifications, W_K2,glob the same-batch and dependent ones,
W_def row DEFJ) by
linearity and the domination below. No independence is claimed: chains may depend on the stopping history.

Non-stopping bound. P(T > g - 1) is P(no resolved useful collision among the first m leaves). The
bound of 6.5 at prefix m (same Gamma, S, n0, q1 at prefix values, global
caps unchanged, hop-1 drop rows summed over the prefix's groups)
bounds it: P(W_m = 0) <= e^{-mu(m)} + (b1+b2)(m) is the Chen-Stein upper tail,
and the usefulness/defect/cap rows count the colliding pairs each mechanism hides, valid at any
m (each row conditions only on coins of prefix batches). Uselessness stays 1/r. The fixup meter cap stays global (2^62) with prefix means: looser,
sound. PHI[g] = min(1, 1 - lo(m)), rounded up to a multiple of 2^-256 (still an upper bound), is nonincreasing, PHI[0] = 1. Fused arrival
keeps the sender-side cap test (proof.md 4), so a dropped record implies >= cap_pg
records addressed to its per-group slot (cursors reset per group, cap_pg a multiple of k, keep
lemma): P(drop) <= P(its slot's addressed count >= cap_pg), which 6.4's Chebyshev row bounds for any sender order. `d_p1` sums that row over the prefix's groups; `d_p2` is unchanged.

Hash chains (detection, no order). Each home keeps, per K2 mod 2^58, a prepend chain of R2 record addresses linked by the packed next-array (every receive pattern up to the slot caps fits); dlvh joins each arriving block's records in [max(F, c), c + count), F the slot's
final cursor after the previous group, so every record is inserted exactly once, in any arrival order (proof.md 4). It links each new record
to its latest same-key predecessor (same K2 and d4), lists that pair unless the list is full, always inserting: same-key records form one connected
chain per key, the listed pairs are among the 182-bit-matching pairs at their home (the d_cl/d_rq/d_rl rows apply unchanged), and a
colliding pair of one class gets its chain edge unless an intruder (same 182 bits, other label) arrives between its records (lam_int).

Uncapped walks (no overflow event). A walk runs to the first full (K2 and d4) match or the chain
end (chains hold <= HOME_MAX = 2^36 slot caps), so it terminates, with visits <= the insert's fold-matching priors, which pair counting (below)
bounds; the worst-case rows charge HOME_MAX visits with verification and listing on every record.

Walk means (E-rows, history-conditional domination with a validity split). H_g = (all X's, batches of groups < g, b, Q)
is the history before group g. {T > g - 1} (no earlier group resolved) is a function of H_g alone: group j's records
are functions of (the X's, batch j's coins); its chains, candidates, requests and replies are deterministic given the
records, and so are the resolver's tests, the usefulness of a listed pair (n' - u*n != 0, n = lambda^i*b_g) and
its [x]P = Q test given also Q. Group g's batch coins are independent of (X's, b, Q) and earlier batches (6.1): fresh given H_g. Group g's walk cost is walks_g = sum_u V_u over its inserts u, with V_u <= sum over priors v
of the fold-match indicator M_uv. Pairs are main (different batches) or same-batch (inside group g). For a
main pair, condition on H_g and on group g's batches except u's, and split by u's validity. If u's batch is valid, v's
fold is fixed, u's X is fixed by H_g and, given (t_k, s_k), the coin c_k is fresh, so u's leaf is
uniform on >= r - 1 points and u's valid-fold function takes the 3.3 mass bound: P(M_uv and u valid | .) <= KB(130) + 1/r
(certificate PAIRM; the 1/r is slack: X = O is excluded by `glob`). If u's batch is
defective (unconditional d_batch per batch, 3.4's marginal derivation), u walks <= HOME_MAX visits: E[V_u x 1_{u defective}]
<= HOME_MAX x d_batch per insert, charged globally (certificate DEFV). Hence E[walks_g^main x 1_{valid} | H_g] <= (pair
count) x (KB130 + 1/r) uniformly in H_g (the count s*gp + s(s-1)/2 - sameX_g is deterministic: arrivals never exceed
sends), and with 1_{T>g-1} H_g-measurable, E[1_{T>g-1} walks_g^main] <= PHI[g] x bound + DEFV. Same-batch pairs
(certificate SAMEX: a batch's 6 q2 leaves) are charged globally without PHI (E[1 x M_uv] <= E[M_uv]; the indicator
dropped, the X's not conditioned): different X's are independent uniform, so one leaf is uniform given the other; the
same X and different points take F3's pair bound; same-incidence leaves take the single-leaf marginal. Each costs <=
KB130 + 1/r (F3's r/(r-1) is inside the 1/r; certificate SXV), plus d_batch for a defective shared batch (inside DEFV).

K2-hit means (E-row, same domination). A K2 hit is a visited pair agreeing on d1, d2 (the home, 72 bits) and K2 (64
bits): 136 of the 200 label bits. For a main pair, condition as for walks; if u's batch is valid, v's 136 bits are fixed
while u's leaf is uniform on >= r - 1 points, so P(K2 hit on (u, v) and u valid | .) <= KB(136) + 1/r (3.3, f = 136; equal K2 implies equal K2 mod 2^58: no bucket factor), and E[1_{T>g-1} k2hits_g^main] <= PHI[g] x ρ x (pair
count) x (KB(136) + 1/r) (certificate `EK2G`). Same-batch pairs are global without PHI, at most their visits, <= KB130 + 1/r
each (`EK2SX`). The dependent terms stay global at KB(64) with the home assumed: the (2S + PARTNER_PAIRS) pairs and the sip
term (`EK2DEP`). Listings drop the stopping indicator and stay global at KB(182) (ELIST, whose sip term covers
same-incidence pairs); a 182-agreeing pair sits at exactly one home (the 182 bits include d1, d2), so ELIST's leading
term takes KB(182) alone, without the fixed-home KB(72) of 6.4's other_h. Both add
<= 2 x d_batch per pair with either side defective (DEFK2, DEFL) and a defective record's own-walk verifies/lists (<=
HOME_MAX each, inside DEFK2/DEFL). A walk of v visits costs at most 4.5v (4 per visit, one B at a chain end after an even
count): the E-rows charge 4.5 per expected visit; a K2 hit's verification (5) and a listing (at most 5) are in EK2 and ELIST. A request
is served from the source's permanent X's and stored batch coins: it regenerates the Y point [c_k + m t_k + n s_k]P (`REQ_SERVE = 2^10 + GEN_MUL` on every request-cap entry, proof.md 7), identically for every group.

Key-path means (E-row, same domination). A leaf takes JH or JT, 1 more instruction (`JPATH`; the sign's one B, `SJOIN`,
is charged always), only if the value its batch computes has exact index 1 or 2 or lies in W (proof.md 3.2, 3.3). If its
batch is valid, that value is the exact leaf x(φ^i(X) ± Y); given H_g (fixing X), group g's other batches and (t_k, s_k),
Y is uniform on G minus <= one point (c_k fresh), hence so is the leaf point (φ and the shift are bijections of G).
Exactly two thirds of G \ {O} have exact index 1 or 2 (φ acts freely on G \ {O}, its orbits carry the three
indices once each, proof.md 3.2), and <= d_wj (r - 1) points lie in W, so P(JH or JT and valid | .) <= 2/3 + d_wj + 1/(r - 1)
<= PJ = 2/3 + 2d_wj + 1/(r - 1) (d_wj > 1/(r - 1)), uniformly in the conditioning. Hence E[keypaths_g x 1_{valid} | H_g] <=
s_g (PJ + 1/3) for the group's s_g leaves, and E[1_{T>g-1} keypaths_g^valid] <= PHI[g] x s_g (PJ + 1/3) (certificate `WGJ`,
inside W_g). A defective batch takes at most 1 on each of its 6q2 leaves past the B charged always: E[keypaths x
1_{defective}] <= d_batch per leaf of the whole run, global (certificate `DEFJ`). There is no key-path meter (JH and JT only
form N'; R_worst charges every leaf 1 + 1/3, `JROW_W`), so 6.4 has no key-path row.

**The 6.1-6.5 argument with these substitutions (E1-E7).** E1: records as in 6.1 (groups are contiguous ranges of batches). E2: cumulative hop-2 slots (6.4's drop rows). E3: hash chains link same-key records (each inserted exactly once). E4: same-key classes are connected (uncapped walk). E5: coverage (a colliding pair with no intruder between its records gets its chain edge; a chain with a useful pair has a useful chain edge, proof.md 4; a listing past the list cap is covered by d_cl). E6: resolution as in 6.3. E7: the success event (stopping only skips groups after a confirmed x, which latches; prefix rows bound prefix defects and caps). So P(success) is at least the full grouped bound of 6.5 at the 0.39 target, and q1 is the least under it (certificate, 145285144).

## Appendix: kernel macro definitions and counts

## 10. Kernel evidence excerpts and count derivation

The macros below are copied verbatim from the `de987687` kernel listings. They are the instruction-level accounting for the dominant certificate constants. The full files and the bit-exact test code are published in `de987687`; this section supplies the relevant macro bodies so the constants can be checked without them.

### PFWD = 153 per pair

```asm
.macro PFWD
	EFWD
	st64b x20,[x3]
	add x3,x3,#1024
	P1X
	.endm
```

```asm
.macro EFWD
	ldp x7,x8,[x3,#-280]
	ldp x20,x21,[x1,#48]
	ldp x22,x23,[x1,#64]
	ldp x8,x9,[x1,#336]
	add x10,x7,x8
	ld64b x12,[x10]
	AC4F x20,x21,x22,x23,x12,x13,x14,x15,x11
	add x10,x7,x9
	ld64b x12,[x10]
	AC4 x20,x21,x22,x23,x16,x17,x18,x19,x11
	ldp x8,x9,[x1,#352]
	add x10,x7,x8
	ld64b x12,[x10]
	AC4 x20,x21,x22,x23,x12,x13,x14,x15,x11
	add x10,x7,x9
	ld64b x12,[x10]
	AC4 x20,x21,x22,x23,x16,x17,x18,x19,x11
	ldp x8,x9,[x1,#368]
	add x10,x7,x8
	ld64b x12,[x10]
	AC4 x20,x21,x22,x23,x12,x13,x14,x15,x11
	add x10,x7,x9
	ld64b x12,[x10]
	AC4 x20,x21,x22,x23,x16,x17,x18,x19,x11
	ldp x8,x9,[x1,#384]
	add x10,x7,x8
	ld64b x12,[x10]
	AC4 x20,x21,x22,x23,x12,x13,x14,x15,x11
	add x10,x7,x9
	ld64b x12,[x10]
	AC4 x20,x21,x22,x23,x16,x17,x18,x19,x11
	LF x20,x21,x22,x23,x11,x10
	.endm
```

```asm
.macro P1X
	RF x20,x24,x25,x26,x27,x13,x14,x15,x16,x17,x28,x19
	SHX x5,x24,x27,x25,x26
	RA x21,x5,x24,x25,x26,x6,x7,x8,x28,x19,x13,x14,x15,x16,x17,x29,xzr,x15,x16
	SHX x27,x5,x26,x24,x25
	RA x22,x27,x5,x24,x25,x6,x7,x8,x28,x19,x13,x14,x15,x16,x17,x29,x29,x15,x16
	SHX x26,x27,x25,x5,x24
	RA x23,x26,x27,x5,x24,x6,x7,x8,x20,x21,x13,x14,x15,x16,x17,x29,x29,x26,x27
	FDX x13,x14,x17,x29,x24,x25,x5,x6,x26,x27
	.endm
```

`PFWD` is `EFWD` + `ST64B` + `ADD` + `P1X`; expanding the three macros gives 153 instructions per pair. The certificate charges 16 PFWD per 32 incidences (4 per pass) plus 5 fixed and one per Y switch: `fwd_count(q2)=5+2452*(q2//32)+q2//256`, 76.6289 per incidence.

### BWDP = 2281 per pair

```asm
.macro BWDP
	ld64b x16,[x3]
	W2S
	stp x0,x1,[x3,#-496]
	stp x2,x24,[x3,#-512]
	C1BS
	INC
	sub x3,x3,#512
	ldp x0,x1,[x3,#32]
	ldp x2,x5,[x3,#48]
	INC
	ldp x25,x26,[x3,#16]
	ldp x27,x28,[x3],#-512
	.endm
```

```asm
.macro W2S
	RF x20,x25,x26,x27,x28,x9,x10,x11,x12,x13,x5,x6
	RF x16,x25,x26,x27,x28,x0,x1,x2,x24,x29,x5,x6
	SHX x5,x25,x28,x26,x27
	RA4 x21,x5,x25,x26,x27,x9,x10,x11,x12,x13,x14,xzr,x6,x7,x8,x28
	RA4 x17,x5,x25,x26,x27,x0,x1,x2,x24,x29,x30,xzr,x6,x7,x8,x28
	SHX x6,x5,x27,x25,x26
	RA4 x22,x6,x5,x25,x26,x9,x10,x11,x12,x13,x14,x14,x7,x8,x28,x27
	RA4 x18,x6,x5,x25,x26,x0,x1,x2,x24,x29,x30,x30,x7,x8,x28,x27
	SHX x7,x6,x26,x5,x25
	RA4 x23,x7,x6,x5,x25,x9,x10,x11,x12,x13,x14,x14,x8,x28,x27,x26
	RA4 x19,x7,x6,x5,x25,x0,x1,x2,x24,x29,x30,x30,x8,x28,x27,x26
	FDX x9,x10,x13,x14,x9,x10,x5,x13,x11,x12
	FDX x0,x1,x29,x30,x0,x1,x6,x29,x2,x24
	.endm
```

```asm
.macro C1BS
	ldp x13,x14,[sp,#16]
	ldp x15,x16,[sp,#32]
	ldp x17,x18,[x3,#192]
	ldp x19,x20,[x3,#208]
	FS x13,x14,x15,x16,x13,x14,x15,x16,x17,x18,x19,x20,x21
	RF x9,x13,x14,x15,x16,x0,x1,x2,x5,x6,x7,x8
	SHX x17,x13,x16,x14,x15
	RA4 x10,x17,x13,x14,x15,x0,x1,x2,x5,x6,x18,xzr,x7,x8,x19,x20
	SHX x16,x17,x15,x13,x14
	RA4 x11,x16,x17,x13,x14,x0,x1,x2,x5,x6,x18,x18,x7,x8,x19,x20
	SHX x15,x16,x14,x17,x13
	RA4 x12,x15,x16,x17,x13,x0,x1,x2,x5,x6,x18,x18,x7,x8,x19,x20
	FDX x0,x1,x6,x18,x0,x1,x7,x8,x2,x5
	stp x0,x1,[x3,#-480]
	stp x2,x5,[x3,#-464]
	ldp x13,x14,[x3,#-128]
	ldp x15,x16,[x3,#-112]
	bfi x13,x9,#8,#32
	ld64b x20,[x13]
	AC4F x0,x1,x2,x5,x20,x21,x22,x23,x6
	lsr x7,x9,#32
	bfi x13,x7,#8,#32
	ld64b x20,[x13]
	AC4 x0,x1,x2,x5,x24,x25,x26,x27,x6
	bfi x14,x10,#8,#32
	ld64b x20,[x14]
	AC4 x0,x1,x2,x5,x20,x21,x22,x23,x6
	lsr x7,x10,#32
	bfi x14,x7,#8,#32
	ld64b x20,[x14]
	AC4 x0,x1,x2,x5,x24,x25,x26,x27,x6
	bfi x15,x11,#8,#32
	ld64b x20,[x15]
	AC4 x0,x1,x2,x5,x20,x21,x22,x23,x6
	lsr x7,x11,#32
	bfi x15,x7,#8,#32
	ld64b x20,[x15]
	AC4 x0,x1,x2,x5,x24,x25,x26,x27,x6
	bfi x16,x12,#8,#32
	ld64b x20,[x16]
	AC4 x0,x1,x2,x5,x20,x21,x22,x23,x6
	lsr x7,x12,#32
	bfi x16,x7,#8,#32
	ld64b x20,[x16]
	AC4 x0,x1,x2,x5,x24,x25,x26,x27,x6
	LF x0,x1,x2,x5,x6,x7
	.endm
```

```asm
.macro INC
	C2UV
	SQ13 x17,x18,x19,x20,x21,x22,x23,x24,x5,x6,x7,x8,x9
	F5 x21,x22,x23,x24,x8,x21,x22,x9,x8
	SQ13 x27,x28,x29,x30,x25,x26,x0,x2,x5,x6,x7,x8,x9
	F5 x25,x26,x0,x2,x8,x25,x26,x9,x8
	PAIRB
	stp x0,x2,[sp,#464]
	str x1,[sp,#480]
	LKE
	TFP4X
	TOUTCX sp,240,256,320,336,128,144
	LKF
	TFP4X
	TOUTCX x3,160,176,288,304,96,112
	LKG
	TFP4X
	TFINS x11,x12
	TFINS x13,x14
	TFINS x15,x16
	TFINS x17,x18
	stp x12,x11,[x3,#256]
	stp x16,x15,[x3,#272]
	stp x14,x13,[x3,#64]
	stp x18,x17,[x3,#80]
	.endm
```

```asm
.macro TOUTCX r,c,d,p,q,m,n
	ldp x9,x10,[\r,#\c]
	ldp x5,x6,[\r,#\d]
	TFINC x11,x12,x9,x10
	TFINC x13,x14,x9,x10
	TFINC x15,x16,x5,x6
	TFINC x17,x18,x5,x6
	stp x12,x11,[x3,#\p]
	stp x16,x15,[x3,#\q]
	stp x14,x13,[x3,#\m]
	stp x18,x17,[x3,#\n]
	.endm
```

```asm
.macro TFP4X
	TROWX TR0,x21,x25,x29,x30,x10,x9,x6,x19
	SHX x19,x5,x10,x6,x9
	TROWX TRK,x22,x26,x27,x28,x9,x6,x5,x10
	SHX x10,x19,x9,x5,x6
	TROWX TRK,x23,x0,x7,x8,x6,x5,x19,x9
	SHXT x10,x6,x19,x5
	TROWX TRK,x24,x2,x1,x20,x5,x19,x10,x6
	.endm
```

```asm
.macro C2UV
	ldp x13,x14,[sp,#112]
	ldp x15,x16,[sp,#128]
	RF x0,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23
	ldp x13,x14,[sp,#144]
	ldp x15,x16,[sp,#160]
	RA4 x1,x13,x14,x15,x16,x17,x18,x19,x20,x21,x24,xzr,x22,x23,x25,x26
	ldp x13,x14,[sp,#176]
	ldp x15,x16,[sp,#192]
	RA4 x2,x13,x14,x15,x16,x17,x18,x19,x20,x21,x24,x24,x22,x23,x25,x26
	ldp x13,x14,[sp,#208]
	ldp x15,x16,[sp,#224]
	RA4 x5,x13,x14,x15,x16,x17,x18,x19,x20,x21,x24,x24,x22,x23,x25,x26
	FDX x17,x18,x21,x24,x17,x18,x22,x23,x19,x20
	ldp x9,x10,[x3,#352]
	ldp x11,x12,[x3,#368]
	bfi x9,x0,#8,#32
	ldp x13,x14,[x9]
	ldp x15,x16,[x9,#16]
	lsr x7,x0,#32
	bfi x9,x7,#8,#32
	ld64b x22,[x9]
	AC4F x13,x14,x15,x16,x26,x27,x28,x29,x6
	bfi x10,x1,#8,#32
	ld64b x22,[x10]
	AC4 x13,x14,x15,x16,x22,x23,x24,x25,x6
	lsr x7,x1,#32
	bfi x10,x7,#8,#32
	ld64b x22,[x10]
	AC4 x13,x14,x15,x16,x26,x27,x28,x29,x6
	bfi x11,x2,#8,#32
	ld64b x22,[x11]
	AC4 x13,x14,x15,x16,x22,x23,x24,x25,x6
	lsr x7,x2,#32
	bfi x11,x7,#8,#32
	ld64b x22,[x11]
	AC4 x13,x14,x15,x16,x26,x27,x28,x29,x6
	bfi x12,x5,#8,#32
	ld64b x22,[x12]
	AC4 x13,x14,x15,x16,x22,x23,x24,x25,x6
	lsr x7,x5,#32
	bfi x12,x7,#8,#32
	ld64b x22,[x12]
	AC4 x13,x14,x15,x16,x26,x27,x28,x29,x6
	LF x13,x14,x15,x16,x6,x7
	FADDW x27,x28,x29,x30,x17,x18,x19,x20,x13,x14,x15,x16,x5
	FS x17,x18,x19,x20,x17,x18,x19,x20,x13,x14,x15,x16,x5
	.endm
```

```asm
.macro LKE
	ldr x1,[x3,#224]
	ldp x5,x6,[sp,#80]
	ldp x9,x10,[sp,#96]
	ldp x0,x2,[sp,#272]
	add x0,x1,x0
	ld64b x12,[x0]
	AC4F x5,x6,x9,x10,x12,x13,x14,x15,x11
	add x2,x1,x2
	ld64b x12,[x2]
	AC4 x5,x6,x9,x10,x16,x17,x18,x19,x11
	ldp x0,x2,[sp,#288]
	add x0,x1,x0
	ld64b x12,[x0]
	AC4 x5,x6,x9,x10,x12,x13,x14,x15,x11
	add x2,x1,x2
	ld64b x12,[x2]
	AC4 x5,x6,x9,x10,x16,x17,x18,x19,x11
	ldp x0,x2,[sp,#304]
	add x0,x1,x0
	ld64b x12,[x0]
	AC4 x5,x6,x9,x10,x12,x13,x14,x15,x11
	add x2,x1,x2
	ld64b x12,[x2]
	AC4 x5,x6,x9,x10,x16,x17,x18,x19,x11
	ldp x0,x2,[sp,#320]
	add x0,x1,x0
	ld64b x12,[x0]
	AC4 x5,x6,x9,x10,x12,x13,x14,x15,x11
	add x2,x1,x2
	ld64b x12,[x2]
	AC4 x5,x6,x9,x10,x16,x17,x18,x19,x11
	LF x5,x6,x9,x10,x11,x12
	ldp x0,x2,[sp,#464]
	ldr x1,[sp,#480]
	.endm
```

```asm
.macro LKF
	ldr x1,[x3,#224]
	ldp x0,x2,[sp,#336]
	add x0,x1,x0
	ldp x5,x6,[x0]
	ldp x9,x10,[x0,#16]
	add x2,x1,x2
	ld64b x12,[x2]
	AC4F x5,x6,x9,x10,x16,x17,x18,x19,x11
	ldp x0,x2,[sp,#352]
	add x0,x1,x0
	ld64b x12,[x0]
	AC4 x5,x6,x9,x10,x12,x13,x14,x15,x11
	add x2,x1,x2
	ld64b x12,[x2]
	AC4 x5,x6,x9,x10,x16,x17,x18,x19,x11
	ldp x0,x2,[sp,#368]
	add x0,x1,x0
	ld64b x12,[x0]
	AC4 x5,x6,x9,x10,x12,x13,x14,x15,x11
	add x2,x1,x2
	ld64b x12,[x2]
	AC4 x5,x6,x9,x10,x16,x17,x18,x19,x11
	ldp x0,x2,[sp,#384]
	add x0,x1,x0
	ld64b x12,[x0]
	AC4 x5,x6,x9,x10,x12,x13,x14,x15,x11
	add x2,x1,x2
	ld64b x12,[x2]
	AC4 x5,x6,x9,x10,x16,x17,x18,x19,x11
	LF x5,x6,x9,x10,x11,x12
	ldp x0,x2,[sp,#464]
	ldr x1,[sp,#480]
	.endm
```

```asm
.macro LKG
	ldr x1,[x3,#224]
	ldp x0,x2,[sp,#400]
	add x0,x1,x0
	ldp x5,x6,[x0]
	ldp x9,x10,[x0,#16]
	add x2,x1,x2
	ld64b x12,[x2]
	AC4F x5,x6,x9,x10,x16,x17,x18,x19,x11
	ldp x0,x2,[sp,#416]
	add x0,x1,x0
	ld64b x12,[x0]
	AC4 x5,x6,x9,x10,x12,x13,x14,x15,x11
	add x2,x1,x2
	ld64b x12,[x2]
	AC4 x5,x6,x9,x10,x16,x17,x18,x19,x11
	ldp x0,x2,[sp,#432]
	add x0,x1,x0
	ld64b x12,[x0]
	AC4 x5,x6,x9,x10,x12,x13,x14,x15,x11
	add x2,x1,x2
	ld64b x12,[x2]
	AC4 x5,x6,x9,x10,x16,x17,x18,x19,x11
	ldp x0,x2,[sp,#448]
	add x0,x1,x0
	ld64b x12,[x0]
	AC4 x5,x6,x9,x10,x12,x13,x14,x15,x11
	add x2,x1,x2
	ld64b x12,[x2]
	AC4 x5,x6,x9,x10,x16,x17,x18,x19,x11
	LF x5,x6,x9,x10,x11,x12
	ldp x0,x2,[sp,#464]
	ldr x1,[sp,#480]
	.endm
```

`BWDP` is `LD64B` + `W2S` + `C1BS` + `INC` + `INC` + two `LDP` + `SUB`; expanding `W2S`, `C1BS`, `INC`, `TOUTCX`, `TFP4X`, `C2UV` and the three lookup blocks gives 2281 per pair. The certificate charges 8 BWDP per 16 incidences (3 per pass) plus 6 fixed and 5 per inner Y switch: `bwd_count(q2)=6+18251*(q2//16)+5*((q2+255)//256-1)`, 1140.7070 per incidence.

### Field inversion and ygen

`fe_inv` is the `secp_fe.s` routine above with its `fe_mul_f`/`fe_sqr_f` helpers; the certificate charges the field calls at its stated multiply/square/add/subtract prices and uses `INV=22549` as the routine allowance (the proof does not use the p256 39479 field-inverse price). `ygen` is the routine above; its count is `A+B*YM+C*YN+D*YM*YN+E*ceil(log2 YM)` with (-23285,5,22598,1606,45680): 211,380,571 at (256,16384,8) (1612.71 per point), and `YGEN=+3*GEN_MUL=213,316,801` per batch.

### Transport and LS64 delivery

```asm
.macro HCHGI u,v,a,l
	bfi x10,\u,#3,#58
	swp x13,x11,[x10]
	stp x11,\u,[x13],#16
	cbz x11,8f
9:
	ldp x5,x6,[x11]
	cmp x6,\u
	b.eq 10f
12:
	cbz x5,8f
	ldp x11,x6,[x5]
	cmp x6,\u
	b.eq 11f
13:
	cbnz x11,9b
	b 8f
10:
	add x6,x15,x11,lsl #1
	ldr x12,[x6,#8]
	eor x12,x12,\v
	tst x12,#0x3fffffffffff
	b.ne 12b
	cmp x3,x4
	b.hs 8f
	\l
	stp x6,\a,[x3],#16
	b 8f
11:
	add x6,x15,x5,lsl #1
	ldr x12,[x6,#8]
	eor x12,x12,\v
	tst x12,#0x3fffffffffff
	b.ne 13b
	cmp x3,x4
	b.hs 8f
	\l
	stp x6,\a,[x3],#16
8:
	.endm
```

```asm
.macro H2B a,b,c,d
	bfi x10,\d,#3,#36
	ldadd x17,x11,[x10]
	stp \a,\b,[x11]
	stp \c,\d,[x11,#16]
	tbz x11,#13,9f
	sub x9,x11,#4080
	sub x9,x9,#4080
	str x9,[x10]
	orr x4,x10,#0x80000000000
	ldadd x18,x11,[x4]
	cmp x11,x12
	b.lo 8f
	str x12,[x4]
	b 9f
8:
	stp x18,x19,[x9,#-16]
	stp x10,x11,[x9,#-32]
9:
	.endm
```

```asm
.macro H2P
	ld64b x20,[x1]
	add x1,x1,#64
	H2B x20,x21,x22,x23
	H2B x24,x25,x26,x27
	.endm
```

```asm
.macro DLVP
	ld64b x16,[x1]
	add x1,x1,#64
	st64b x16,[x2]
	add x2,x2,#64
	.endm
```

```asm
.macro DLVS
	ldp x16,x17,[x1],#32
	ldp x18,x19,[x1,#-16]
	stp x18,x19,[x2,#16]
	stp x16,x17,[x2],#32
	.endm
```

```asm
.macro DLVHS
	DLVS
	HCHGI x16,x17,x12,LSUB
	.endm
```

Counts: dlv 6 per one-record message; h2s 6 per record plus 8 fixed and the flush/copy terms; dlvh 4; join 4 per record; walk 4.5 per visit. The itemized charges are in proof.md section 4 and certificate.py.

### Counts and parameter consistency

fwd 76.6289 and bwd 1140.7070 per incidence; fe_inv 22549; ygen 211,380,571 per batch; q1=145285144, G=160, success 0.3900000000833567; q2=2^25, MX=256, L_Y=2^17. The certificate inserts these constants and checks the exact rational ledger; the macros above are the instruction-level source for them.

### Arithmetic primitives used by PFWD/BWDP

```asm
.macro RF o,a,b,c,d,e,f,g,h,i,t,u
	mul \e,\o,\a
	umulh \t,\o,\a
	mul \u,\o,\b
	adds \f,\t,\u
	umulh \t,\o,\b
	mul \u,\o,\c
	adcs \g,\t,\u
	umulh \t,\o,\c
	mul \u,\o,\d
	adcs \h,\t,\u
	umulh \t,\o,\d
	adc \i,\t,xzr
	.endm
```

```asm
.macro RA4 o,b0,b1,b2,b3,a0,a1,a2,a3,a4,a5,y,r,s,t,v
	mul \r,\o,\b0
	adds \a0,\a0,\r
	umulh \s,\o,\b0
	mul \r,\o,\b1
	adcs \s,\s,\r
	umulh \t,\o,\b1
	mul \r,\o,\b2
	adcs \t,\t,\r
	umulh \v,\o,\b2
	mul \r,\o,\b3
	adcs \v,\v,\r
	umulh \r,\o,\b3
	adc \r,\r,xzr
	adds \a1,\a1,\s
	adcs \a2,\a2,\t
	adcs \a3,\a3,\v
	adcs \a4,\a4,\r
	adc \a5,\y,xzr
	.endm
```

```asm
.macro RA o,a,b,c,d,r,s,t,v,u,e,f,g,h,i,j,y,m,n
	mul \r,\o,\a
	umulh \s,\o,\a
	mul \u,\o,\b
	adds \s,\s,\u
	umulh \t,\o,\b
	mul \u,\o,\c
	adcs \t,\t,\u
	umulh \v,\o,\c
	mul \u,\o,\d
	adcs \v,\v,\u
	umulh \u,\o,\d
	adc \u,\u,xzr
	adds \e,\e,\r
	adcs \f,\f,\s
	adcs \m,\g,\t
	adcs \n,\h,\v
	adcs \i,\i,\u
	adc \j,\y,xzr
	.endm
```

```asm
.macro SHX l,a,d,p,q
	mul \l,\d,x4
	umulh \d,\d,x4
	adds \a,\a,\d
	b.lo 9f
	adcs \p,\p,xzr
	adcs \q,\q,xzr
	stadd x4,[sp]
9:
	.endm
```

```asm
.macro SHXT a,d,p,q
	umulh \d,\d,x4
	adds \a,\a,\d
	b.lo 9f
	adcs \p,\p,xzr
	adcs \q,\q,xzr
	stadd x4,[sp]
9:
	.endm
```

```asm
.macro FDX e,f,i,j,o,p,l,h,y,z
	mul \l,\i,x4
	umulh \h,\i,x4
	madd \h,\j,x4,\h
	adds \o,\e,\l
	adcs \p,\f,\h
	b.lo 9f
	adcs \y,\y,xzr
	adc \z,\z,xzr
	stadd x4,[sp]
9:
	.endm
```

```asm
.macro SQ13 a0,a1,a2,a3,w0,c1,c2,c3,c4,c5,c6,c7,t
	mul \c3,\a0,\a3
	mul \t,\a1,\a2
	adds \c3,\c3,\t
	umulh \c4,\a0,\a3
	umulh \t,\a1,\a2
	adcs \c4,\c4,\t
	umulh \c5,\a1,\a3
	adc \c5,\c5,xzr
	umulh \c2,\a0,\a1
	mul \t,\a0,\a2
	adds \c2,\c2,\t
	umulh \t,\a0,\a2
	adcs \c3,\c3,\t
	mul \t,\a1,\a3
	adcs \c4,\c4,\t
	mul \t,\a2,\a3
	adcs \c5,\c5,\t
	umulh \c6,\a2,\a3
	adc \c6,\c6,xzr
	mul \c1,\a0,\a1
	lsr \c7,\c6,#63
	extr \c6,\c6,\c5,#63
	extr \c5,\c5,\c4,#63
	extr \c4,\c4,\c3,#63
	extr \c3,\c3,\c2,#63
	extr \c2,\c2,\c1,#63
	mul \w0,\a0,\a0
	umulh \t,\a0,\a0
	adds \c1,\t,\c1,lsl #1
	mul \t,\a1,\a1
	adcs \c2,\c2,\t
	umulh \t,\a1,\a1
	adcs \c3,\c3,\t
	mul \t,\a2,\a2
	adcs \c4,\c4,\t
	umulh \t,\a2,\a2
	adcs \c5,\c5,\t
	mul \t,\a3,\a3
	adcs \c6,\c6,\t
	umulh \t,\a3,\a3
	adc \c7,\c7,\t
	mul \t,\c4,x4
	adds \w0,\w0,\t
	umulh \c4,\c4,x4
	mul \t,\c5,x4
	adcs \c4,\c4,\t
	umulh \c5,\c5,x4
	mul \t,\c6,x4
	adcs \c5,\c5,\t
	umulh \c6,\c6,x4
	mul \t,\c7,x4
	adcs \c6,\c6,\t
	umulh \c7,\c7,x4
	adc \c7,\c7,xzr
	adds \c1,\c1,\c4
	adcs \c2,\c2,\c5
	adcs \c3,\c3,\c6
	adc \c7,\c7,xzr
	.endm
```

```asm
.macro F5 w0,w1,w2,w3,w4,o0,o1,ta,tb
	mul \ta,\w4,x4
	umulh \tb,\w4,x4
	adds \o0,\w0,\ta
	adcs \o1,\w1,\tb
	b.lo 9f
	adcs \w2,\w2,xzr
	adc \w3,\w3,xzr
	stadd x4,[sp]
9:
	.endm
```

```asm
.macro PAIRB
	mov x9,#0xfff8000000000000
	orr x5,x9,#64
	orr x6,x9,#128
	BLK0 x21,x22,x23,x24,x29,x27,x7,x1
	BLK0 x25,x26,x0,x2,x30,x28,x8,x20
	.endm
```

```asm
.macro TROWX op,a,b,c,d,r3,r2,r1,t
	\op \a,x11,x12,\r3,\r2,\r1,\t
	\op \b,x13,x14,\r3,\r2,\r1,\t
	\op \c,x15,x16,\r3,\r2,\r1,\t
	\op \d,x17,x18,\r3,\r2,\r1,\t
	.endm
```

```asm
.macro TR0 m,f,s,r3,r2,r1,t
	umulh \f,\m,\r2
	madd \f,\m,\r3,\f
	mul \t,\m,\r2
	lsr \s,\t,#3
	umulh \t,\m,\r1
	add \s,\s,\t,lsr #3
	.endm
```

```asm
.macro TRK m,f,s,r3,r2,r1,t
	umulh \t,\m,\r2
	add \f,\f,\t
	madd \f,\m,\r3,\f
	mul \t,\m,\r2
	add \s,\s,\t,lsr #3
	umulh \t,\m,\r1
	add \s,\s,\t,lsr #3
	.endm
```

```asm
.macro AC4 a0,a1,a2,a3,b0,b1,b2,b3,c
	adds \a0,\a0,\b0
	adcs \a1,\a1,\b1
	adcs \a2,\a2,\b2
	adcs \a3,\a3,\b3
	adc \c,\c,xzr
	.endm
```

```asm
.macro AC4F a0,a1,a2,a3,b0,b1,b2,b3,c
	adds \a0,\a0,\b0
	adcs \a1,\a1,\b1
	adcs \a2,\a2,\b2
	adcs \a3,\a3,\b3
	cset \c,hs
	.endm
```

```asm
.macro LF a0,a1,a2,a3,c,t
	mul \t,\c,x4
	adds \a0,\a0,\t
	adcs \a1,\a1,xzr
	b.lo 9f
	adcs \a2,\a2,xzr
	adcs \a3,\a3,xzr
	csel \t,x4,xzr,hs
	add \a0,\a0,\t
	stadd x4,[sp]
9:
	.endm
```

```asm
.macro FS d0,d1,d2,d3,a0,a1,a2,a3,b0,b1,b2,b3,t
	subs \d0,\a0,\b0
	sbcs \d1,\a1,\b1
	sbcs \d2,\a2,\b2
	sbcs \d3,\a3,\b3
	csel \t,x4,xzr,lo
	subs \d0,\d0,\t
	b.hs 9f
	sbcs \d1,\d1,xzr
	sbcs \d2,\d2,xzr
	sbc \d3,\d3,xzr
	stadd x4,[sp]
9:
	.endm
```

```asm
.macro FADDW d0,d1,d2,d3,a0,a1,a2,a3,b0,b1,b2,b3,t
	adds \d0,\a0,\b0
	adcs \d1,\a1,\b1
	adcs \d2,\a2,\b2
	adcs \d3,\a3,\b3
	csel \t,x4,xzr,hs
	adds \d0,\d0,\t
	b.lo 9f
	adcs \d1,\d1,xzr
	adcs \d2,\d2,xzr
	adc \d3,\d3,xzr
	stadd x4,[sp]
9:
	.endm
```

```asm
.macro W2S
	RF x20,x25,x26,x27,x28,x9,x10,x11,x12,x13,x5,x6
	RF x16,x25,x26,x27,x28,x0,x1,x2,x24,x29,x5,x6
	SHX x5,x25,x28,x26,x27
	RA4 x21,x5,x25,x26,x27,x9,x10,x11,x12,x13,x14,xzr,x6,x7,x8,x28
	RA4 x17,x5,x25,x26,x27,x0,x1,x2,x24,x29,x30,xzr,x6,x7,x8,x28
	SHX x6,x5,x27,x25,x26
	RA4 x22,x6,x5,x25,x26,x9,x10,x11,x12,x13,x14,x14,x7,x8,x28,x27
	RA4 x18,x6,x5,x25,x26,x0,x1,x2,x24,x29,x30,x30,x7,x8,x28,x27
	SHX x7,x6,x26,x5,x25
	RA4 x23,x7,x6,x5,x25,x9,x10,x11,x12,x13,x14,x14,x8,x28,x27,x26
	RA4 x19,x7,x6,x5,x25,x0,x1,x2,x24,x29,x30,x30,x8,x28,x27,x26
	FDX x9,x10,x13,x14,x9,x10,x5,x13,x11,x12
	FDX x0,x1,x29,x30,x0,x1,x6,x29,x2,x24
	.endm
```

```asm
.macro TFINC f,s,cl,ch
	add \f,\f,\s,lsr #61
	adds \s,\cl,\s,lsl #3
	adc \f,\f,\ch
	.endm
```

```asm
.macro TFINS f,s
	add \f,\f,\s,lsr #61
	.endm
```

These are the remaining macros invoked by PFWD/BWDP and the lookup/fold paths; with the macro bodies above they close the instruction-level expansions.

