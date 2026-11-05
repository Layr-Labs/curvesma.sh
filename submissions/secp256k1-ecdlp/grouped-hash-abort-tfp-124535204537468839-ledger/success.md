# Success, section 6 of proof.md (no heuristic)

### 6.1 Model

Coins are independent uniform; Q is uniform. Each `Y_l` is uniform on G \ {O}; given R_g, the
`X_k = R_g + [a_k]P` are independent uniform. So the X's are independent uniform points, independent of the Y's;
nodes, hence leaves of different nodes, are independent; a leaf is uniform on G. Given all generator
points, the `b_g` stay independent uniform: for each value of b_g exactly one tuple `(a_k)` fits the points.

Γ is the set of unordered cross-node leaf pairs, `|Γ| = (S² - 2^72*n0²)/2`. For α = {L, L'} ∈ Γ let
`I_α = 1[L' ∈ class(L)]`; then `P(I_α = 1) = P1 = (6(r-1) + 1)/r²` and `μ = |Γ|*P1 = 0.494468123` (the run's `mu`).
Same-node pairs are omitted; their collisions only help. Collisions are single pairs, not clusters (cross-node: four distinct points): a match of `φ^t(X) + sY` with `u(φ^t'(X') + s'Y')` forces a second match among the same two incidences' other leaves only through a second independent relation among X, Y, X', Y' (fresh coins), probability ≤ 1/(r − 1), already inside P1 and the neighbourhood sums of 6.2.

### 6.2 Chen–Stein on cross-node pairs

Let `W = Σ_(α∈Γ) I_α` and `B_α` be the pairs sharing a generator point with α. I_α depends only on α's four
generators (distinct, from two nodes), so is independent of `{I_β : β ∉ B_α}` (b3 = 0). A generator
lies in 6q2 (an X) or 6q1 (a Y) leaves, so `|B_α| <= 12(q1 + q2)*S`.

For a unit u, the event `L' = u*L` with `L = φ^t(X) + s*Y` and `L' = φ^(t')(X') + s'*Y'` is the linear relation
`λ^(t')x' + s'y' - uλ^t x - us*y = 0` in the generators' logarithms; I_α <= the sum of these 6 events. Take
α ≠ β in Γ sharing a generator, and units u, v. **The two relations are not proportional:** proportional relations
have the same support, so each node gives exactly one X and one Y to both, and β's leaf in each node uses both:
β is `{φ^(t2)(X) + s2*Y, φ^(t2')(X') + s2'*Y'}` with event `L2' = v*L2`. A factor c with `λ^(t2') = c*λ^(t')` and
`s2' = c*s'` is both ±1 and a power of λ, so c = 1 (-1 is not: powers of λ cube to 1, (-1)³ = -1),
`t2' = t'` and `s2' = s'`; then `vλ^(t2) = uλ^t` and `v*s2 = u*s` force `λ^(t2 - t) = s2/s`, so `t2 = t`, `s2 = s` and
β = α. So some 2 × 2 minor is invertible: given all other generators, the relations fix two generators, each
uniform on >= r - 1 points, so `E[I_α I_β] <= 36/(r-1)²`. Chen-Stein (Arratia, Goldstein and Gordon, Ann. Probab. 17 (1989), Theorem 1)
gives `|P(W = 0) - e^-μ| <= b1 + b2 <= |Γ|*12(q1 + q2)S*(P1² + 36/(r-1)²)`, certificate-checked < 2^-95.

### 6.3 Usefulness

A cross-node colliding pair is **detected** if both records carry their class's label (no batch defect, wrong class
index or label error), both survive the caps, no intruder enters their run, and no candidate, request or relay cap
is reached first. Its run is then resolved, a useful adjacent pair (section 4) gives x, and `[x]P = Q` is checked
before output. The pair is useful iff `n' ≠ u*n`, where `n = λ^t*b_g` and `n' = λ^(t')*b_(g')` for its two nodes; given
all points, this fails with probability exactly 1/r. The resolver tries the records' nine type pairs,
solving with the first relation having `n' - u*n ≠ 0` (2.3): any relation it computes between the points, right or
wrong types, gives the logarithm, so a useful true pair's relation is never lost. So `P(fail) <= P(W = 0) + E #{α
: I_α = 1, α useless or not detected}`.

### 6.4 Defect rows

A **node-local defect** of L (3.4 batch defect, wrong class index, label error, fixup meter) is a function
of node g's coins, given which L' is uniform, so `P(I_α ∧ defect) <= (6/r)*P(defect)`; over Γ a row costs <=
`μ(1 + 1/r)` times its per-leaf rate.

**Caps and intruders.** Condition on α's four generators. A leaf with a generator outside them (a **free** generator)
is then uniform on >= r - 1 points, so <= `KB(f)` of its mass falls on given values of f label bits,
whatever index KEY returns (3.3's lattice count). A leaf's emitted fields depend on its incidence's components but lie
in its value's cell widened by 3.3's error ranges: each count below is <= the number of leaves with value in
the widened cell, and a pair event with one leaf's fields fixed lies in an event of the other leaf's value, of
mass <= KB(f) whatever the components; we bound these value events. For two leaves of one
incidence, the second emits its class label outside 3.3's sets W and D', else (probability < `d_wj +
d_lab`) one of <= 12 labels near it (three paths, two blocks per fraction), so their pair event has mass <=
`(1 + 11(d_wj + d_lab))KB(f)`; these 5 pairs per leaf join the candidate, request and relay loads without the
home condition (`sip`). The 12 leaves of α's two incidences are fixed; the six leaves of an incidence through one of α's
generators and a free one are functions of the free one.
- **Pair slots** count single leaves. A pair slot has mean <= `m_p = n0*KB(36) + 1` from the free leaves (at hop 1
  one node's leaves of one group with a given d1, n0 the group's; at hop 2 one column's with a given (d1, d2), n0 all
  groups'); each cap leaves room for the <= 12 fixed leaves. Two leaves'
  indicators have covariance 0 without a common free generator or with a common free Y; <= `1/(4(r-1)²)` with a
  common free X (<= `KB(f)/(r-1)` if one also uses α's Y); <= `1/(r²(r-1))` times the two cells'
  sizes for two leaves of one incidence with both generators free, whose pair is uniform on r(r - 1) pairs (an
  invertible linear image of (X, Y), with Y ≠ O); and <= one leaf's mean for two leaves of one incidence through
  α's generators, <= 60(q1 + q2) ordered pairs. So a count's variance is <= its mean plus
  `1 + 60(q1 + q2)*KB(f)`, and Chebyshev bounds a record's drop by `(m + 1 + 60(q1 + q2)KB(f))/(c - 12 - m)²`
  for cap c and mean bound m (`d_p1` per group with cap_pg at hop 1, `d_p2` with cap_p at hop 2).
- **Intruders and the candidate, request and relay loads** count leaf pairs, jointly with I_α. For a pair
  π = {A, B} ≠ α with event E_π (keys equal on f bits; at the candidate list and the relay, also A's home or home column
  is α's), `E[I_α*1(E_π)] <= ρ*(6/r)*P*(π)`, where `P*(π)` is E_π's probability for independent uniform leaves and
  `ρ = (r/(r-1))^4` covers the exclusion Y ≠ O. (i) If A, say, has a free generator B lacks,
  condition on all else: A is uniform on >= r - 1 points while I_α, B and α's home are fixed, so the f-bit
  match costs <= `(r/(r-1))KB(f)` and the home condition falls on B. If B has a free generator, B is uniform given
  α's generators, and its home is α's with probability <= `(r/(r-1))KB(72)` (KB(36) for a home column); if B is a
  partner of L (or L'), (L, B) (or (L', B)) is uniform up to ρ over that node's coins, given which `P(I_α) <= 6/r`; if
  B is L or L', π is a pair through L (below). (ii) Otherwise A and B lie in
  one incidence or in α's incidences. If their incidence avoids α's generators, condition on α's: (A, B) is an
  invertible linear image of the incidence's two generators, so uniform up to ρ. If both lie in α's node g and their
  incidence uses α's X or Y there (α's own incidence included), E_π (with α's home, which on I_α is L's) is a function
  of g's coins, given which L' is uniform, so `P(I_α | g's coins) <= 6/r`; over g's coins (A, B) is uniform up to ρ, and
  L is uniform given (A, B) through a generator they lack, unless A is L (pairs through L, or L',
  are counted without the home condition, both at α's home) or A and B are both partners of L (20 pairs, 10 per
  incidence of α, also counted without it). If A is in L's incidence and B in L''s, `(L, A)` and `(L', B)` are uniform
  up to ρ, so each of the six events `L' = uL` meets E_π with probability <= `ρ*P*(π)/r`.

  So the intruder rate is <= `ρ(S + 12)KB(142)` per colliding pair; the candidate list at α's home has mean
  `other_h <= ρ(S²*KB(72)*KB(182) + (2S + 20)*KB(182))` (pairs at the home equal on d1, d2, K2 and d4, the 182 label
  bits a listed pair matches; L and L' are there); a source receives <= `other_g <= 2ρ*n0*S*KB(182)` requests; a
  relay holds <= `other_v <= 18ρ((S/2^36)*S*KB(36)*KB(182) + (2S + 20)*KB(182))` records (a request is one record
  and a reply <= 8). Markov with the caps less α's own entries (4 in the candidate list, 8 requests, 16 records)
  gives `d_cl`, `d_rq` and `d_rl`. Leaves of defective batches elsewhere are not uniform; the whole-run row "a
  computational defect in any batch" excludes them.

| row (per leaf of the pair unless stated) | term in `certificate.py` | checked |
| --- | --- | --- |
| useless collision (per pair) | `1/r` | - |
| batch defect: the sites of 3.4, D = 0, X = O | `d_batch` | < 2^-126.05 |
| wrong class index (3.3): KEY's window set W | `d_wj` | < 2^-58.47 |
| fixup meter, Markov on the mean `q1*E_FIX` with the cap `2^45 - 33q2` (3.4) | `d_f` | < 2^-16.32 |
| the label-error set D' (3.3) | `d_lab` | < 2^-16.65 |
| pair slot at hop 1, per group, summed over the groups (Chebyshev, `cap_pg = cap_p/8`) | `d_p1` | < 2^-16.5 |
| pair slot at hop 2 (Chebyshev, `cap_p`) | `d_p2` | < 2^-19.9 |
| intruder (per pair), `ρ(S + 12)KB(142)` | `lam_int` | < 2^-12.7 |
| candidate, request and relay caps (Markov; a reply counts 8 records) | `d_cl`; `d_rq`; `d_rl` | < 2^-23.4; < 2^-25.9; < 2^-25.2 |
| a computational defect in any batch (whole run) | `g_def = 2^72*q1*d_batch` | < 2^-26.5 |
| two coin range rejections (whole run) | in `glob` | all of `glob` < 2^-26.5 |
| RNDR failure (whole run) | `g_rndr` = 0 under the rules' ideal coins (section 7) | - |

### 6.5 The bound

`certificate.py` picks, by bisection from 2^27 (where it fails), the least q1 for which

```text
1 - e^-mu - (b1 + b2) - mu(1 + 1/r)(1/r + intruder + 2 * sum of per-leaf rows) - whole-run rows >= 0.39,
```

the per-leaf rows being 6.4's (the hop-1 pair-slot row d_p1 for the 114 groups of 6.6, the hop-2 row d_p2 once each),
with an upper enclosure of e^-μ and exact Fractions: q1 = 145307201, with left side >= 0.3900000005;
the target is the rules' 0.39 (searched there, as in aburan28's 7516608f).

## 6.6 Grouped early abort: stopping time, hash chains, E-work

Before the run, the q1 X-batches are split into G = 114 contiguous groups (sizes differ
by <= one; certificate `group_sizes`). All nodes run each group in lockstep
(barriers between groups). Group g generates its X's and runs the leaf, route, hash-join and resolve
phases on them exactly as in 6.1-6.5 (same leaf kernels and keys/fin send path, same coins-per-batch laws, fused
hop-1 arrival (records processed on block arrival, never written to the allocated,
unused R1 slots, proof.md 4), cumulative hop-2 slots of cap_p with per-slot
final cursors F (E2 by group parity; group 1 the slot bases), and the cumulative hash table of 2^58 zeroed head words with a packed next-array
(one 16-byte entry, next address and K2, per R2 slot position) and uncapped walks); then resolves group g's new
candidates, stopping the run at the first group resolving to an x with [x]P = Q
(checked; the output latches). If none resolves, the run ends after group G without output. T
in 1..G+1 is the first resolving group (G+1 if none); group g's leaf/route work and hash fixed cost
are deterministic given g and run iff T > g - 1, while its main-pair walks have history-conditional
means bounded uniformly in H_g (below), so with P_g = P(T > g - 1) <= PHI[g] the expected
work is E_pre + sum_g PHI[g] (W_g + W_K2,g) + W_sameX + W_K2,glob + W_list + W_def + W_rest (worst-case rows, with
candidates and fixups at their expected counts, proof.md 7; W_g includes the key-path mean below, W_K2,g the main pairs' K2 verifications, W_K2,glob the same-X and dependent ones,
and W_def its defective-batch row DEFJ) by
linearity and the history-conditional domination below. No independence
is claimed: chains may depend on stopping's history.

Non-stopping bound. P(T > g - 1) is P(no resolved useful collision among the first m leaves). The
bound of 6.5 at prefix m (same Gamma, with S, n0, q1 at prefix values, global
caps unchanged, per-group hop-1 drop rows summed over the prefix's groups; no overflow row: walks
are uncapped) bounds it: P(W_m = 0) <= e^{-mu(m)} + (b1+b2)(m) is the Chen-Stein upper tail,
and the usefulness/defect/cap rows count the colliding pairs each mechanism hides, valid at any
m as each row conditions only on coins of prefix batches. Uselessness stays 1/r
given the points. The fixup meter cap stays global (2^45) with prefix means: looser,
sound. The certificate's PHI[g] is min(1, 1 - lo(m)), nonincreasing, PHI[0] = 1. Fused arrival
keeps the sender-side cap test: `leaf.s` FLUSH and `route.s` fin still drop a block whose RC cursor
is at/past the slot end (proof.md 4), so on every run a dropped record implies >= cap_pg
records addressed to its per-group slot (cursors reset per group, cap_pg a multiple of k, keep
lemma). Hence P(a record is dropped) <= P(its slot's addressed count >= cap_pg), the
count event 6.4's Chebyshev row bounds: emit, cursor rule and per-group reset are
unchanged, so the mean/variance argument holds verbatim for any sender order. The
certificate's `d_p1` sums that per-group row over the prefix's groups; the hop-2 row
`d_p2` is unchanged.

Hash chains (detection, no order). Bucket b at home h holds all stored records with K2 mod 2^58 = b
as a prepend chain of R2 record addresses linked by the packed next-array (one 16-byte entry, next
address and K2, per R2 slot position: every receive pattern up to the slot caps fits by construction,
no bump, no capacity check); dlvh joins each arriving block's records in [max(F, c), c + count), F the slot's
final cursor after the previous group, so every record is inserted exactly once, at arrival, in any arrival order
in the group (proof.md 4). It links each new record
to its latest same-key predecessor (same K2 and d4; the walk is uncapped, so on every run),
lists that pair unless the list is full, always inserting: same-key records form one connected
arrival-order chain per key, no sorting needed, the listed pairs are among the
182-bit-matching pairs at their home (the d_cl/d_rq/d_rl pair-counting rows apply unchanged), and a
colliding pair of one class gets its chain edge unless an intruder (same 142 bits, other label)
arrives between its records (the lam_int row).

Uncapped walks (no overflow event). The walk runs to the first full (K2 and d4) match or the chain
end; chains are finite (one home's records of one fold, <= HOME_MAX = 2^36 slot caps), so every
walk terminates on every run, with no overflow row. Each visit examines a same-fold prior,
so an insert's visits are <= its fold-matching priors, and pair counting (below)
bounds uncapped walks; the worst-case rows charge HOME_MAX visits with verification and
listing on every record of every group.

Walk means (E-rows, history-conditional domination with a validity split). H_g = (X_<g, Y, b, Q)
is the history before group g: earlier groups' X-batches, all Y generators, every node's
coefficient b_g' (so every R_g' = [b_g']Q) and Q. {T > g - 1} (no earlier group
resolved) is a function of H_g alone: group j's records are functions of (X_j, Y); its chains,
candidates, requests and replies are deterministic given the records, and the resolver's unit tests,
the usefulness of a listed pair (n' - u*n != 0 with n = lambda^t*b_g, a function of the coefficients)
and its [x]P = Q test given also the scalars and Q; so H_g decides
stopping by g. Conditioning on b and Q is free: given (b_g, Q), node g's group-g points
X_k = R_g + [a_k]P (fresh independent uniform a_k) stay independent uniform and
independent of the rest of H_g (6.1), all the fresh-X step uses. Group g's walk cost is walks_g
= sum_u V_u over its inserts u, with V_u <= sum over priors v of the fold-match indicator M_uv on
every run (as above). Pairs are main (different
X-batches) or same-X (same batch, so inside group g). For a main pair, condition on H_g and
on all of group g's X-batches except u's, and split by u's validity. If u's batch is valid, v's fold is fixed (whether v is valid or corrupt) while u's X is fresh uniform (a_k is independent of H_g and
the other batches), so u's valid-fold function takes
the 3.3 mass bound: P(M_uv and u valid | .) <= KB(130) + 1/r (KB(130) for X != O, 1/r for X = O
where the leaf is fixed). If u's batch is defective (unconditional d_batch per batch, 3.4's fresh-Y
derivation, never conditioned on Y), u walks <= HOME_MAX visits: E[V_u x 1_{u defective}] <=
HOME_MAX x d_batch per insert, charged globally (certificate DEFV). Hence E[walks_g^main x
1_{valid} | H_g] <= (pair count) x (KB130 + 1/r) uniformly in H_g — the pair count s*gp +
s(s-1)/2 - sameX_g is deterministic, since arrivals never exceed sends — and 1_{T>g-1} is
H_g-measurable, so E[1_{T>g-1} walks_g^main] <= PHI[g] x bound + DEFV: the domination is established
before multiplying by PHI[g]. Same-X pairs are charged globally without PHI (E[1 x M_uv] <= E[M_uv]; given
shared X, different-Y leaves are independent and same-incidence leaves take the single-leaf
marginal): <= KB130 + 1/r each (certificate SXV) plus d_batch for a defective shared batch
(inside DEFV). K2-hit means (E-row, the same domination, joint mass bound). A K2 hit is a visited pair agreeing
on d1, d2 (the home, 72 bits) and K2 (64 bits): 136 of the 200 label bits. For a main pair, condition as for walks;
if u's batch is valid, v's 136 bits are fixed while u's X is fresh uniform, so P(K2-hit on (u, v) and u valid | .)
<= KB(136) + 1/r (3.3's mass bound with f = 136; 1/r for X = O, as for walks; equal K2 implies equal K2 mod 2^58,
so the bucket adds no factor). Main pairs take walks' deterministic pair counts (s*gp + s(s-1)/2 - sameX_g) and
1_{T>g-1} is H_g-measurable, so E[1_{T>g-1} k2hits_g^main] <= PHI[g] x ρ x (pair count) x (KB(136) + 1/r)
(certificate `EK2G`). Same-X pairs are charged globally without PHI: a K2 hit is a visit, so their hits are at
most their visits, <= KB130 + 1/r per pair as above (certificate `EK2SX`). The dependent terms stay global as
before, at KB(64) with the home assumed: the (2S + 20) pairs and the same-incidence sip term (certificate
`EK2DEP`). Listings drop the stopping indicator (E[1*X] <= E[X]) and stay global at KB(182) (ELIST, whose sip
term covers same-incidence pairs). Both add <= 2 x d_batch per pair with either side defective (DEFK2, DEFL) and
a defective record's own-walk verifies/lists (<= HOME_MAX each, inside DEFK2/DEFL). ELIST's leading term is global: a 182-agreeing pair sits at exactly one home
(the 182 bits include d1, d2), so takes KB(182) alone with no fixed-home KB(72) (6.4's other_h
validly keeps it for a specified collision's home). A walk of v visits costs at most 4.5v (4 per visit and one B
at a chain end after an even count v >= 2), so the E-rows charge
4.5 per expected visit; a K2 hit's verification (5: ADD, LDR, EOR, TST, B.NE) and a listing (at most 5: 4, + 1 forming an address off x2) are in EK2 (EK2G, EK2SX, EK2DEP) and ELIST. Old-X
requests serve from each source's retained all-groups X affine store (proof.md 2.3, 5), identically
for every group, so cross-group candidates need no regeneration row.

Key-path means (E-row, the same domination). A leaf takes JH or JT, 1 more instruction (`JPATH`; the sign's one B, `SJOIN`, is charged always), only if the
value its batch computes has exact index 1 or 2 or lies in W (proof.md 3.2, 3.3). If its batch
is valid, the computed value is the exact leaf x(φ^t(X_k) ± Y_l); given H_g (fixing R_g) and
group g's other X-batches, X_k = R_g + [a_k]P is fresh uniform on G minus <= one point, hence so is the leaf
point (φ and the shift are bijections of G). Exactly two thirds of
G \ {O} have exact index 1 or 2 (φ acts freely on G \ {O}; its orbits {P, φP, φ²P} carry the three indices once
each, the uniformity of proof.md 3.2), and <= d_wj (r - 1) points lie in W, so P(JH or JT and valid | ·) <=
2/3 + d_wj + 1/(r - 1) <= PJ = 2/3 + 2d_wj + 1/(r - 1), uniformly in the conditioning (d_wj > 1/(r - 1)); the domination is established before multiplying by PHI[g]. Hence
E[keypaths_g x 1_{valid} | H_g] <= s_g (PJ + 1/3) for the group's s_g leaves, and as 1_{T>g-1} is H_g-measurable,
E[1_{T>g-1} keypaths_g^valid] <= PHI[g] x s_g (PJ + 1/3) (certificate `WGJ`, inside W_g). A defective batch (unconditional d_batch per batch, 3.4's fresh-Y derivation, never
conditioned on Y) takes at most 1 on each of its 6q2 leaves past the B charged always, so E[keypaths x 1_{defective}] <= d_batch per leaf of
the whole run, charged globally without PHI (certificate `DEFJ`). There is no key-path meter: JH and JT only form N'
and a sign takes at most one B, so R_worst charges every leaf 1 + 1/3 (certificate `JROW_W`) and 6.4 has no key-path row.

**The 6.1-6.5 argument with these substitutions (E1-E7).** E1: records as in 6.1 (the same kernels and
X/Y order; groups are contiguous X-ranges). E2: cumulative hop-2 slots (hop-2 cursors persist across groups; 6.4's
drop rows). E3: hash chains link same-key records (dlvh joins each block's
records at or past F, so the joined ranges partition each slot and each record is inserted exactly once; no drops, buckets or capacity check). E4: same-key classes are connected
(each new record links to its latest same-key predecessor, reached by the uncapped walk). E5:
coverage (every colliding pair with no intruder between its records gets its chain edge; a
chain with a useful pair has a useful chain edge, proof.md 4; a listing past the list cap is
covered by the d_cl rows). E6: resolution as in 6.3. E7: the success
event (stopping only skips groups after a confirmed x, which latches; prefix rows bound prefix
defects and caps). So P(success) is at least the full grouped bound of 6.5 at the 0.39
target, and q1 is the least under it (certificate, 145307201).
