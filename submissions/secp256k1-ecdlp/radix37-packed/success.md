# Success, section 6 of proof.md (no heuristic)
### 6.1 Model
Coins are independent. Q is uniform. Node g draws `b_g` and `a_j` (j < MX = 1024) uniform in [0, r), giving its permanent points `X_j = [b_g]Q + [a_j]P`. Each of its q1 batches k draws `c_k, t_k, s_k` uniform in [1, r) and `rho_k` uniform in [1, p). The batch's L = 2^15 points are `Y_(m,n) = [c_k + m t_k + n s_k]P` (m < 8192, n < 4); incidence l pairs `X_(l mod 1024)` with the point of number l >> 10, so each X meets each point once per batch. The running products start at `rho_k`: uniform on F_p^* while no factor is 0. Facts:
- (F1) Given b, the `X_j` are i.i.d. uniform; so all X's are i.i.d. uniform, independent of the batches. Given all points, `b_g` stays independent uniform: for each value of b_g exactly one tuple `(a_j)` fits the points.
- (F2) Given `(t_k, s_k)`, each Y point is uniform on G \ {one point} (`c_k` uniform in [1, r)), independent of the X's and of other batches.
- (F3) For index pairs v != w of one batch, `Y_w - Y_v = [dm t + dn s]P` with `(dm, dn) != (0, 0)`, |dm| < 8192, |dn| < 4, which takes each value with probability <= 1/(r - 1) over `(t, s)`. So `(Y_v, Y_w)` takes each pair of points with probability <= `1/((r-1)(r-1))` (the r/(r-1) factors are inside ρ = (r/(r-1))^4), and given one point of a batch, any other takes each value with probability <= 1/(r - 1) (up to ρ).
- (F4) For three points with affinely independent index pairs, `(c, t, s) -> (Y_u, Y_v, Y_w)` has a nonzero integer determinant of absolute value < r: the triple is uniform up to ρ.
Nodes, hence leaves of different nodes, are independent; a leaf `φ^i(X) + eY` (i = 0, 1, 2; e = ±1) is uniform on G (X uniform, independent of Y).
Γ is the set of unordered cross-node leaf pairs, `|Γ| = (S² - 2^72*n0²)/2` with `S = 2^72*6*q1*q2`, `n0 = 6*q1*q2`. For α = {L, L'} ∈ Γ let `I_α = 1[L' ∈ class(L)]`; then `P(I_α = 1) = P1 = (6(r-1) + 1)/r²` and `μ = |Γ|*P1 ≈ 6.749995172` (`mu`, an opened run's).
Same-node pairs are omitted; their collisions only help. Collisions are single pairs, not clusters (cross-node: four distinct points): a match of `φ^i(X) + eY` with `u(φ^i'(X') + e'Y')` forces a second match among the same incidences' other leaves only through a second independent relation among X, Y, X', Y' (fresh coins), probability <= 1/(r - 1), already inside P1 and 6.2's sums.
### 6.2 Chen–Stein on cross-node pairs
Let `W = Σ_(α∈Γ) I_α` and `B_α` be the pairs sharing an X or a batch with α. I_α depends only on α's four generators (two X points and two batches, from two nodes), which are independent of all other generators (6.1), so I_α is independent of `{I_β : β ∉ B_α}` (b3 = 0). An X lies in `6*q1*L` leaves, a batch in `6*q2`, so `|B_α| <= 12(q1*L + q2)*S`.
For a unit u, the event `L' = u*L` with `L = φ^i(X) + eY`, `L' = φ^(i')(X') + e'Y'` is a linear relation in the eight coins `(x, c, t, s, x', c', t', s')` (x, x' the logarithms of X, X'; Y's is `c + m t + n s`): coefficients `-uλ^i, -ue, -uem, -uen, λ^(i'), e', e'm', e'n'`. I_α <= the sum of these 6 events. Take α != β in Γ sharing a generator, units u, v; E_β is β's relation.
- **Case 1: β has a generator outside α's four.** Its x (or c) enters E_β with a unit coefficient, not E_α; with an x or c of α entering E_α with a unit coefficient it forms a triangular, invertible 2 × 2 minor.
- **Case 2: β has only α's generators.** Restrict both relations to `(x, c, x', c')`, rows `(-uλ^i, -ue, λ^(i'), e')` and `(-vλ^(i2), -ve2, λ^(i2'), e2')`. If they are not proportional, some 2 × 2 minor is invertible. **If proportional,** the factor κ has `e' = κ e2'`, so κ = ±1, and `λ^(i') = κ λ^(i2')`, so κ is a power of λ; -1 is not (powers of λ cube to 1, (-1)³ = -1), so κ = 1, `i2' = i'`, `e2' = e'`; then `uλ^i = vλ^(i2)` and `ue = ve2` give `v/u = e/e2 = ±1 = λ^(i-i2)`, so v = u, `e2 = e`, `i2 = i`. E_β - E_α is `e'(dm' t' + dn' s') - ue(dm t + dn s) = 0`, `(dm, dn, dm', dn')` the index differences of β's points from α's, not all zero (else β = α): probability <= 1/(r - 1) over `(t, s, t', s')`; given these and x, x', c', E_α fixes c: <= 1/(r - 1).
So the relations fix two coins (uniform on >= r - 1 points given the others), and `E[I_α I_β] <= 36/(r-1)²`. Chen-Stein (Arratia, Goldstein and Gordon, Ann. Probab. 17 (1989), Theorem 1) gives `|P(W = 0) - e^-μ| <= b1 + b2 <= |Γ|*12(q1*L + q2)S*(P1² + 36/(r-1)²)`, certificate-checked < 2^-73.
### 6.3 Usefulness
A cross-node colliding pair is **detected** if both records carry their class's label (no batch defect, wrong class
index or label error), both survive the caps, no intruder enters their run, and no candidate, request or relay cap
is reached first. Its run is then resolved (each record's batch k decoded exactly, proof.md 4), a useful adjacent pair (section 4) gives x, and `[x]P = Q` is checked
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
- **Canonical routing cells (slot loads only).** Let `v_m = β^m x mod p`, `X_m = 2^64 v_m/p`, ε = 2^-64. The wrap set is where some T(v_m) is below 1024 or at least 2^128 - 1024: at most 6*1024*2^128 values of x. Outside it, by 3.3's ranges (T(v_m) - W_m in [-132, 194], m = 0, 1) W_m does not wrap, `F_m` lies in `(X_m - 1 - 195ε, X_m + 132ε]`, and (T_0 + T_1 + T_2 = s 2^128 - {0, 1, 2}, s = 1 or 2) `σ = 2^64 - X_2 - δ`, δ in `(-265ε, 1 + 390ε)` (one unit wider than 3.3's interval, from the floor in T_2). KEY's branches give, for its index j and every m: j = 0 or 1 without the σ carry, `X_j < X_2 + 1 + 585ε` (`F_j <= X_2 + δ - 1`) and `X_j < X_(1-j) + 1 + 327ε` (the F comparison); j = 2 with the carry, `X_2 <= min(F_0, F_1) - δ < X_m + 397ε`. So **`X_j < X_m + 1 + 585ε` for every m** (`GAPE` = 585, derived in `certificate.py` from the ranges). The emitted pair estimates `(T(u), T(w))` for the point `(u, w) = (v_j, βv_j mod p)` of Λ; with `z = -u - w mod p` and macrocell `(i, l) = (⌊2^64 u/p⌋, ⌊2^64 w/p⌋)`, `X_u < X_w + 1 + 585ε` gives `l >= i - 2`; if u + w < p, `X_z = 2^64 - X_u - X_w` gives `2i + l <= 2^64 + 1`; else `X_z = 2^65 - X_u - X_w` gives `2i + l <= 2^65 + 1` and `i + l >= 2^64 - 1`. These lie in the pad-8 set whose size NC < (1 + 2^-58) 2^128/3 `certificate.py` counts exactly (`ncells`, brute-checked for sides <= 60): the canonical domain is a third of the square. A fixed 36-bit digit meets a macrocolumn in <= 2 T-intervals of 2^28 + 654 values in all (3.3's widening), so in at most `CELLU = (2^28 + 654) 2^128` values of u; a macrorow has <= 2^192 values of w. By 3.3's lattice lemma (each interval widened by D = max(A, 2A - B)), a macrocell holds <= `(CELLU + 2D)(2^192 + D)/p` points of Λ with a given d1 and <= `(CELLU + 2D)²/p` with a given (d1, d2). Summing over the three bijections `x -> (v_m, v_(m+1))` onto Λ and two points per x, a leaf uniform on >= r - 1 points emits a given d1 with probability <= `KB_R36 = 6 NC (CELLU + 2D)(2^192 + D)/(p(r-1)) + WRAPB` and a given (d1, d2) with probability <= `KB_H72 = 6 NC (CELLU + 2D)²/(p(r-1)) + WRAPB`, `WRAPB = 12*1024*2^128/(r-1)` (the wrap set): at most 2.0000049 2^220/(r-1) and 2.0000098 2^184/(r-1), about a third of KB(36), KB(72) (checked: below them / 2.99). This uses only the error ranges, so it holds in a clean batch whatever index KEY returns (W, D' included). It replaces KB(36) only in the slot means and covariances below; every other event (KB(64), KB(72), KB(136), KB(182), the candidate, request and relay loads) keeps the factor-6 KB, since a high-bit cell can lie wholly inside the canonical domain. `experiment.py` `canon_key` checks the gap and the macrocell set for every x of p = 32779 over KEY's full error box.
- **Pair slots** count single leaves. A slot has mean <= `m_p` from the free leaves: at hop 1 `n0*KB_R36 + 1` (one node's leaves of one group with a given d1, n0 the group's), at hop 2 `2^36 n0*KB_H72 + 1` (one column's with a given (d1, d2), n0 all
  groups'); each cap leaves room for the <= 12 fixed leaves. Two leaves' indicators have covariance 0 unless they share an X or a batch. A common free X (different batches) gives <= `1/(4(r-1)²)` (<= `KB(f)/(r-1)` with α's Y point). Outside (a)-(c), two points of one batch with the same X or a common free X, or one incidence with both generators free, are uniform up to ρ (F3, F4; an invertible linear image of (X, Y) with Y != O): covariance <= `3/(r²(r-1))` times the two cells' sizes. The ordered pairs of (a)-(c) have covariance <= one leaf's cell mean and number <= `COVX(q) = 60(q*L + MX) + 72 L² + 12 MX*L`, for the two nodes of α:
  (a) the 30 ordered leaf pairs of an incidence through exactly one of α's points: `q*L` through an X per group (q its batches), MX through a Y point;
  (b) two leaves on α's X and batch, functions of `(t, s)` given α's point: `(6L)²` per node;
  (c) a free X with two points of α's batch symmetric about α's, `Y_v + Y_w = 2Y_u`, of one type and opposite signs: `6*MX*L` per node. No other pair through one free X and α's batch is dependent: its two leaves (types i, i2) are an invertible image of (X, c, t, s) unless `P_w - P_u = ±λ^(i2-i)(P_v - P_u)` in the index plane; for i2 != i that makes λ = ±a/b with 0 < |a|, |b| < 8196, which `LAMBDA_FREE` rules out (0 < a² ± ab + b² < r); for i2 = i the sign + gives P_w = P_v.
  So a count's variance is <= its mean plus `1 + COVX(q)*KB_R36`, and Chebyshev bounds a record's drop by `(m + 1 + COVX(q)KB_R36)/(c - 12 - m)²`
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
  are counted without the home condition, both at α's home) or α's point and the two leaves' points are collinear in the index plane: those pairs, conservatively at most `72 L^2` ordered, and the 20 pairs of two partners of L in an incidence of α (10 per incidence), are counted without it too: `PARTNER_PAIRS`. If A is in L's incidence or a partner of L and B likewise for L', `(L, A)` and `(L', B)` are independent and uniform
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
| fixup meter, Markov on the mean `q1*E_FIX` with the cap `2^62 - FIX_BATCH` (28 q2 fixups, 3.4), plus `q1*d_batch` | `d_f` | < 2^-32.5 |
| the label-error set D' (3.3) | `d_lab` | < 2^-16.65 |
| pair slot at hop 1, per group (rounded up to a multiple of 2^-256), summed over groups (Chebyshev, `cap_pg`); at hop 2 (`cap_p`) | `d_p1`; `d_p2` | < 2^-18.3; < 2^-18 |
| intruder (per pair), `ρ(S + 12)KB(182)` | `lam_int` | < 2^-50.8 |
| candidate, request and relay caps (Markov; a reply counts 8 records) | `d_cl`; `d_rq`; `d_rl` | < 2^-19.6; < 2^-22.2; < 2^-21.4 |
| a computational defect in any batch (whole run) | `g_def = 2^72*q1*d_batch` | < 2^-26.1 |
| an identity `X_j` (`2^72*MX/r`); a coin whose two range-rejection draws fail (`2^-254`; 4 per batch, MX + 1 per node) (whole run) | in `glob` | all of `glob` < 2^-26.1 |
| RNDR failure (whole run) | `g_rndr` = 0 under the rules' ideal coins (section 7) | - |
`d_batch` is unconditional per batch (proof.md 3.4: a bad set of size b costs <= `b*mult/(r-1)`, Y uniform on >= r - 1 points by F2; products uniform on F_p^* by `rho_k` while no factor is 0; ygen's start and step divisions `YGEN_DEF`; equal or opposite points `L^2/(r-1)`). A Y point equal to O comes only from an addition with equal x-coordinates (C != O as c_k >= 1), which sets ygen's flag: it is inside `YGEN_DEF`.
### 6.5 The bound (an opened run)
`certificate.py` fixes the opened run's length, q1 = 536870720 (320*1677721 equal batches; k decoded by proof.md 4), and computes
```text
lo_run = 1 - e^-mu - (b1 + b2) - mu(1 + 1/r)(1/r + intruder + 2 * sum of per-leaf rows) - whole-run rows,
```
the per-leaf rows being 6.4's (the hop-1 pair-slot row d_p1 for the 320 groups of 6.6, the hop-2 row d_p2 once each), with an upper enclosure of e^-μ (μ approximately 6.749995172 > 1: the |x| <= 1 series bounds squared k times, exact Fractions rounded outward to 2^-512) and exact Fractions: lo_run >= 0.9985833838. With the canonical cells the hop-2 slots run at mean 2^21.585 of cap_p = 2^22 - 2^17 (d_p2 < 2^-18) and the hop-1 slots at 2^13.264 of cap_pg; the label-error row d_lab < 2^-16.65 is the largest; every row is charged.
## 6.6 Grouped early abort: stopping time, hash chains, E-work
Before the run, the q1 batches are split into G = 320 contiguous ranges (groups) of Y-batches, group g ending at
1677721*g (exact: `group_sizes`; all groups have 1677721 batches). All nodes run each group in lockstep. Group g generates its batches' Y points (the X's are permanent, set up once) and runs the leaf, route, hash-join and resolve
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
**The 6.1-6.5 argument with these substitutions (E1-E7).** E1: records as in 6.1 (groups are contiguous ranges of batches). E2: cumulative hop-2 slots (6.4's drop rows). E3: hash chains link same-key records (each inserted exactly once). E4: same-key classes are connected (uncapped walk). E5: coverage (a colliding pair with no intruder between its records gets its chain edge; a chain with a useful pair has a useful chain edge, proof.md 4; a listing past the list cap is covered by d_cl). E6: resolution as in 6.3. E7: the success event (stopping only skips groups after a confirmed x, which latches; prefix rows bound prefix defects and caps). So an opened run succeeds with probability at least the full grouped bound lo_run of 6.5 (certificate, q1 = 536870720), and 6.7 adds the gate.
### 6.7 Run gate
One coin u, uniform on [0, 2^64) and independent of the challenge and of every other coin, opens the run iff u < K (proof.md 2.4). The run's success event is a function of the challenge and the run's coins only, and conditioning on u < K leaves their joint law unchanged, so 6.1-6.6 hold for an opened run verbatim and P(success) >= P(u < K) lo_run = (K/2^64) lo_run. `certificate.py` takes K = 7204436109379059621, the least integer with (K/2^64) lo_run >= 0.39 (it checks K and K - 1): success >= 0.39 (exact rational lower bound; downward decimal). Expected work: a closed gate costs GATE = 2^7 instructions, an opened one GATE plus the run's work, whose expectation (over the run's coins and the challenge, the same convention as 6.6) is E_run; so E = 2^7 + (K/2^64) E_run (proof.md 7). The gate is independent of the stopping history H_g of 6.6, so PHI and every E-row are an opened run's unchanged.

### The 8192 x 4 rectangle
F1-F4 hold for the 8192 x 4 index set (L = 32768, MX*L = 2^25): index differences are nonzero integer vectors with |dm| <= 8191, |dn| <= 3, and affinely independent triples have determinants of absolute value at most 2*8191*3 < r. The lambda test uses the bound 8196 (3*8194^2 < r; lambda^2 + lambda + 1 = 0). So the rank exclusions and COVX hold: through a free X only (c)'s charged symmetric partner can lose rank; through α's own X every ordered point pair is charged (72 L^2: L^2 pairs, 36 leaf-type pairs, two endpoints; with the 20 incidence partners, PARTNER_PAIRS), overcounting every collinear line. No independence hypothesis is added. YGEN_DEF is exact ((M - 1) + M(N - 1) = MN - 1). PHI, every row, K and the gate are exact-Fraction recomputations: the gated bound is at least 39/100 and K - 1 fails; finite tests are not success evidence.

### Radix-37 representation change
The lookup decomposition in proof.md 3.1 changes only a weak representation of the same field operand. Seven canonical terms plus an optional base replace eight plus a base, all below the inherited nine-term LF bound. All eleven defect-site and seventeen fixup rows remain conservative; all uses of a computed operand already allow every congruent representative below 2^256. No bit-identical comparison with epoch29 is assumed. The independent batch/point laws, truncated-product error box and canonical routing-cell counts are unchanged. The smaller cap_p = 2^22 - 2^17 changes the hop-2 Chebyshev row at every prefix: certificate.py recomputes it, all PHI values and the least gate K from those prefixes; d_p2 < 2^-18 and lo_run >= 0.9985833838.
