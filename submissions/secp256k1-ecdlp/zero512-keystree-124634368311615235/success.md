# Success, section 6 of proof.md (no heuristic)

### 6.1 Model

Coins are uniform and independent, and Q is uniform. Each `Y_l` is uniform on G \ {O}; given R_g, the
`X_k = R_g + [a_k]P` are independent and uniform. So the X's are independent uniform points, independent of the Y's, and
nodes are independent; a leaf is uniform on G, and leaves of different nodes are independent. Given all generator
points, the `b_g` are still independent and uniform: for each value of b_g exactly one tuple `(a_k)` fits the points.

Γ is the set of unordered pairs of leaves from different nodes, `|Γ| = (S² - 2^72*n0²)/2`. For α = {L, L'} ∈ Γ let
`I_α = 1[L' ∈ class(L)]`; then `P(I_α = 1) = P1 = (6(r-1) + 1)/r²` and `μ = |Γ|*P1 = 0.494455702` (the run's `mu`).
Same-node pairs are left out; their collisions only add to success.

### 6.2 Chen–Stein on cross-node pairs

Let `W = Σ_(α∈Γ) I_α`, and let `B_α` be the pairs that share a generator point with α. I_α depends only on α's four
generators, distinct because they come from two nodes, so I_α is independent of `{I_β : β ∉ B_α}` (b3 = 0). A generator
lies in 6q2 (an X) or 6q1 (a Y) leaves, so `|B_α| <= 12(q1 + q2)*S`.

For a unit u, the event `L' = u*L` with `L = φ^t(X) + s*Y` and `L' = φ^(t')(X') + s'*Y'` is the linear relation
`λ^(t')x' + s'y' - uλ^t x - us*y = 0` in the generators' logarithms; I_α is at most the sum of these 6 events. Take
α ≠ β in Γ sharing a generator, and units u and v. **The two relations are not proportional:** proportional relations
have the same support, so each node contributes exactly one X and one Y to both, and β's leaf in each node uses both:
β is `{φ^(t2)(X) + s2*Y, φ^(t2')(X') + s2'*Y'}` and its event is `L2' = v*L2`. A factor c with `λ^(t2') = c*λ^(t')` and
`s2' = c*s'` is both ±1 and a power of λ, so c = 1 (-1 is not: every power of λ cubes to 1, and (-1)³ = -1),
`t2' = t'` and `s2' = s'`; then `vλ^(t2) = uλ^t` and `v*s2 = u*s` force `λ^(t2 - t) = s2/s`, so `t2 = t`, `s2 = s` and
β = α. So some 2 × 2 minor is invertible: given all
other generators, the relations fix two generators, each uniform on at least r - 1 points, and
`E[I_α I_β] <= 36/(r-1)²`. The Chen–Stein method (Arratia, Goldstein and Gordon, Ann. Probab. 17 (1989), Theorem 1)
gives `|P(W = 0) - e^-μ| <= b1 + b2 <= |Γ|*12(q1 + q2)S*(P1² + 36/(r-1)²)`, which the certificate checks is below
2^-95.

### 6.3 Usefulness

A cross-node colliding pair is **detected** if both records carry their class's label (no batch defect, wrong class
index or label error), both survive the caps, no intruder enters their run, and the candidate, request and relay caps
are not reached first. Its run is then resolved, a useful adjacent pair (section 4) gives x, and `[x]P = Q` is checked
before output. The pair is useful iff `n' ≠ u*n`, where `n = λ^t*b_g` and `n' = λ^(t')*b_(g')` for its two nodes; given
all points, this fails with probability exactly 1/r. The resolver tries the nine type pairs of the two records and
solves with the first relation that has `n' - u*n ≠ 0` (2.3): a relation between the points it computes, of right or
wrong types, gives the logarithm, so the true pair's relation, if useful, is never lost. So `P(fail) <= P(W = 0) + E #{α
: I_α = 1, α useless or not detected}`.

### 6.4 Defect rows

A **node-local defect** of L (a batch defect of 3.4, a wrong class index, a label error, the fixup meter) is a function
of node g's coins, and L' is uniform given them, so `P(I_α ∧ defect) <= (6/r)*P(defect)`; over Γ a row costs at most
`μ(1 + 1/r)` times its per-leaf rate.

**Caps and intruders.** Condition on α's four generators. A leaf with a generator outside them, a **free** generator,
is then uniform on at least r - 1 points, so at most `KB(f)` of its mass falls on given values of f label bits,
whatever index KEY returns (3.3's lattice count). A leaf's emitted fields depend on its incidence's components, but lie
in its value's cell widened by 3.3's error ranges: each count below is at most the number of leaves whose value lies in
the widened cell, and a pair event in which one leaf's fields are fixed lies in an event of the other leaf's value, of
mass at most KB(f) whatever the components; the bounds below are for these value events. For two leaves of one
incidence, the second emits its class label outside 3.3's sets W and D', and otherwise (probability below `d_wj +
d_lab`) one of at most 12 labels near it (three paths, two blocks per fraction), so their pair event has mass at most
`(1 + 11(d_wj + d_lab))KB(f)`; these 5 pairs per leaf are added to the candidate, request and relay loads without the
home condition (`sip`). The 12 leaves of α's two incidences are fixed; the six leaves of an incidence through one of α's
generators and a free one are functions of that free generator.
- **Pair slots and d3 buckets** count single leaves. A pair slot has mean at most `m_p = n0*KB(36) + 1` from the free
  leaves (at hop 1 it holds one node's leaves with a given d1, at hop 2 one column's leaves with a given (d1, d2), the
  same bound), and a d3 bucket `m_3 = S*KB(96) + 1`; each cap leaves room for the at most 12 fixed leaves. Two leaves'
  indicators have covariance 0 without a common free generator or with a common free Y; at most `1/(4(r-1)²)` with a
  common free X (at most `KB(f)/(r-1)` when one of them also uses α's Y); at most `1/(r²(r-1))` times the two cells'
  sizes for two leaves of one incidence with both generators free, whose pair is uniform on r(r - 1) pairs (an
  invertible linear image of (X, Y), with Y ≠ O); and at most one leaf's mean for two leaves of one incidence through
  α's generators, at most 60(q1 + q2) ordered pairs. So a count's variance is at most its mean plus
  `1 + 60(q1 + q2)*KB(f)`, and Chebyshev bounds the drop of a given record by `(m + 1 + 60(q1 + q2)KB(f))/(c - 12 - m)²`
  for a cap c and mean bound m (`d_p`, `d_3`).
- **Intruders and the candidate, request and relay loads** count pairs of leaves, jointly with I_α. For a pair
  π = {A, B} ≠ α with event E_π (keys equal on f bits; at the candidate list and the relay, also A's home or home column
  is α's), `E[I_α*1(E_π)] <= ρ*(6/r)*P*(π)`, where `P*(π)` is E_π's probability for independent uniform leaves and
  `ρ = (r/(r-1))^4` covers the exclusion Y ≠ O. (i) If one of them, say A, has a free generator that the other lacks,
  condition on everything else: A is uniform on at least r - 1 points while I_α, B and α's home are fixed, so the f-bit
  match costs at most `(r/(r-1))KB(f)` and the home condition falls on B. If B has a free generator, B is uniform given
  α's generators, and its home is α's with probability at most `(r/(r-1))KB(72)` (KB(36) for a home column); if B is a
  partner of L (or L'), (L, B) (or (L', B)) is uniform up to ρ over that node's coins, given which `P(I_α) <= 6/r`; if
  B is L or L', π is a pair through L (below). (ii) Otherwise A and B lie in
  one incidence or in α's incidences. If their incidence is free of α's generators, condition on α's: (A, B) is an
  invertible linear image of the incidence's two generators, so uniform up to ρ. If both lie in node g of α and their
  incidence uses α's X or Y there (α's own incidence included), E_π (with α's home, which on I_α is L's) is a function
  of g's coins, given which L' is uniform, so `P(I_α | g's coins) <= 6/r`; over g's coins (A, B) is uniform up to ρ, and
  L is uniform given (A, B) through a generator they lack, unless A is L (the pairs through L, and likewise L',
  are counted without the home condition, as both are at α's home) or A and B are both partners of L (20 pairs, 10 per
  incidence of α, also counted without it). If A is in L's incidence and B in L''s, `(L, A)` and `(L', B)` are uniform
  up to ρ, so each of the six events `L' = uL` meets E_π with probability at most `ρ*P*(π)/r`.

  Hence the intruder rate is at most `ρ(S + 12)KB(142)` per colliding pair; the candidate list at α's home has mean
  `other_h <= ρ(S²*KB(72)*KB(182) + (2S + 20)*KB(182))` (pairs at the home equal on d1, d2, K2 and d4, the 182 label
  bits `scan_ranges` matches; L and L' are there); a source receives at most `other_g <= 2ρ*n0*S*KB(182)` requests; a
  relay holds at most `other_v <= 18ρ((S/2^36)*S*KB(36)*KB(182) + (2S + 20)*KB(182))` records (a request is one record
  and a reply at most 8). Markov with the caps less α's own entries (4 in the candidate list, 8 requests, 16 records)
  gives `d_cl`, `d_rq` and `d_rl`. Leaves of a defective batch elsewhere are not uniform; the whole-run row "a
  computational defect in any batch" excludes them.

| row (per leaf of the pair unless stated) | term in `certificate.py` | checked |
| --- | --- | --- |
| useless collision (per pair) | `1/r` | - |
| batch defect: the sites of 3.4, D = 0, X = O | `d_batch` | < 2^-126.05 |
| wrong class index (3.3): KEY's window set W | `d_wj` | < 2^-58.47 |
| fixup meter, Markov on the mean `q1*E_FIX` with the cap `2^45 - 33q2` (3.4) | `d_f` | < 2^-16.32 |
| the label-error set D' (3.3) | `d_lab` | < 2^-16.65 |
| pair slot at hop 1, and again at hop 2 (Chebyshev, `cap_p`) | `d_p`, counted twice | < 2^-19.9 each |
| d3 bucket (Chebyshev, `cap_3 = ceil(m_3) + 2^32`) | `d_3` | < 2^-30 |
| intruder (per pair), `ρ(S + 12)KB(142)` | `lam_int` | < 2^-12.7 |
| candidate, request and relay caps (Markov; a reply counts 8 records) | `d_cl`; `d_rq`; `d_rl` | < 2^-23.4; < 2^-25.9; < 2^-25.2 |
| a computational defect in any batch (whole run) | `g_def = 2^72*q1*d_batch` | < 2^-26.5 |
| two coin range rejections (whole run) | in `glob` | all of `glob` < 2^-26.5 |
| RNDR failure (whole run) | `g_rndr` = 0 under the rules' ideal coins (section 7) | - |

### 6.5 The bound

`certificate.py` picks the least q1 for which

```text
1 - e^-mu - (b1 + b2) - mu(1 + 1/r)(1/r + intruder + 2 * sum of per-leaf rows) - whole-run rows >= 0.390000001,
```

with an upper enclosure of e^-μ and exact Fractions: q1 = 145305376, where the left side is at least 0.3900000023, and
the target is 10^-9 above the rules' 0.39.

## 6.6 Grouped early abort: stopping time, hash chains, E-work

The q1 X-batches are partitioned once, before the run, into G = 114 contiguous groups (sizes differ
by at most one; certificate `group_sizes`). Every node runs the same group at the same time
(barriers between groups). Group g: generate its X's; run the leaf, route, hash-join and resolve
phases for them exactly as in 6.1-6.5 (same leaf kernels and keys/fin send path, same coins-per-batch laws, fused
hop-1 arrival (records processed on block arrival, never written to R1 slots; the allocated
slots stand unused, proof.md 4), cumulative hop-2 slots of cap_p with per-slot
hash-read cursors, and the cumulative hash table of 2^57 zeroed head words with a packed next-array
(one 16-byte entry, next address and K2, per R2 slot position) and uncapped walks); then resolve group g's new
candidates and stop the whole run at the first group whose candidates resolve to an x with [x]P = Q
(checked; the output latches). If no group resolves, the run ends after group G with no output. T
in 1..G+1 is the first resolving group (G+1 if none); group g's leaf/route work and hash fixed cost
are deterministic given g and run iff T > g - 1, while its main-pair walks have history-conditional
means bounded uniformly in the history (below), so with P_g = P(T > g - 1) <= PHI[g] the expected
work is E_pre + sum_g PHI[g] W_g + W_sameX + W_K2 + W_list + W_def + W_rest (worst-case rows; W_g includes the
key-path mean below and W_def its defective-batch row DEFJ) by
linearity over indicators and history-conditional domination with a validity split. No independence
is claimed: the chains may depend on the same history as stopping.

Non-stopping bound. P(T > g - 1) is P(no resolved useful collision among the first m leaves). The
bound of 6.5 re-evaluated at the prefix m (same Gamma with S, n0, q1 at their prefix values, global
caps unchanged, per-group hop-1 drop rows summed over the groups inside; no overflow row: walks
are uncapped) is an upper bound: P(W_m = 0) <= e^{-mu(m)} + (b1+b2)(m) is the Chen-Stein upper tail,
and the usefulness/defect/cap rows count colliding pairs hidden by each mechanism, all valid at any
m since every row conditions only on coins of batches inside the prefix. Uselessness stays 1/r
given the points. The fixup meter cap is global (2^45) with prefix means, hence looser and
sound. The certificate's PHI[g] is min(1, 1 - lo(m)), nonincreasing, PHI[0] = 1. Fused arrival
keeps the sender-side cap test: `leaf.s` FLUSH and `route.s` fin still drop a block whose RC cursor
is at or past the slot end (proof.md 4), so on every run a dropped record implies at least cap_pg
records addressed to its per-group slot (cursors reset per group, cap_pg a multiple of k, keep
lemma). Hence P(a given record is dropped) <= P(its slot's addressed count >= cap_pg), the same
count event 6.4's Chebyshev row bounds: the emit, the cursor rule and the per-group reset are
unchanged, so the mean/variance argument applies verbatim, whatever the sender order. The
certificate's `d_p1` sums that per-group row over the groups inside the prefix; the hop-2 row
`d_p2` applies unchanged.

Hash chains (detection, no order). Bucket b at home h holds all stored records with K2 mod 2^57 = b
as a prepend chain of R2 record addresses linked by the packed next-array (one 16-byte entry, next
address and K2, per R2 slot position: every receive pattern up to the slot caps fits by construction,
no bump, no capacity check); each group call reads its slots from their read cursors to E2 and advances the
cursors, so every record is read and inserted exactly once (proof.md 4). hchg links each new record
to its most recent same-key predecessor (same K2 and d4; the walk is uncapped, so on every run),
lists that pair unless the list is full, and always inserts: same-key records form one connected
arrival-order chain per key, no sorted order is needed, the listed pairs are a subset of the
182-bit-matching pairs at their home (the d_cl/d_rq/d_rl pair-counting rows apply unchanged), and a
colliding pair of one class gets its chain edge unless an intruder (same 142 bits, other label)
arrives between its records, the event the lam_int row covers.

Uncapped walks (no overflow event). The walk runs to the first full (K2 and d4) match or the chain
end; chains are finite (one home's records of one fold, at most HOME_MAX = 2^36 slot caps), so every
walk terminates on every run with no overflow event or row. Every visit examines a same-fold prior,
so an insert's visit count is at most its number of fold-matching priors, and pair counting (below)
bounds uncapped walks directly; the worst-case rows charge HOME_MAX visits with verification and
listing on every record of every group.

Walk means (E-rows, history-conditional domination with a validity split). Write H_g = (X_<g, Y, b, Q)
for the history before group g: the X-batches of all earlier groups, all Y generators, every node's
coefficient b_g' (so every R_g' = [b_g']Q) and the input Q. The event {T > g - 1} (no earlier group
resolved) is a function of H_g alone: group j's records are functions of (X_j, Y); its chains,
candidates, requests and replies are deterministic given the records; the resolver's unit tests,
the usefulness of a listed pair (n' - u*n != 0 with n = lambda^t*b_g, a function of the coefficients)
and its [x]P = Q test are deterministic given the records, the scalars and Q; so stopping by g is
decided by H_g. Conditioning on b and Q costs nothing below: given (b_g, Q), node g's group-g points
X_k = R_g + [a_k]P with their fresh uniform independent a_k are still independent uniform points,
independent of everything else in H_g (6.1), which is all the fresh-X step uses. Group g's walk cost is walks_g
= sum_u V_u over its inserts u, with V_u <= sum over priors v of the fold-match indicator M_uv on
every run (every visit examines a same-fold prior). Split the pairs: main pairs (different
X-batches) and same-X pairs (same batch, so inside group g). For a main pair, condition on H_g and
on all of group g's X-batches except u's, and split by u's validity. If u's batch is valid, v's fold is fixed (a fixed
value whether v is valid or corrupt) while u's X is fresh uniform (a_k is independent of H_g and of
the other batches), so u's valid-fold function takes
the 3.3 mass bound: P(M_uv and u valid | .) <= KB(129) + 1/r (KB(129) for X != O, 1/r for X = O
where the leaf is fixed). If u's batch is defective (unconditional d_batch per batch, 3.4's fresh-Y
derivation, never conditioned on Y), u walks at most HOME_MAX visits: E[V_u x 1_{u defective}] <=
HOME_MAX x d_batch per insert, charged globally (certificate DEFV). Hence E[walks_g^main x
1_{valid} | H_g] <= (pair count) x (KB129 + 1/r) uniformly in H_g — the pair count s*gp +
s(s-1)/2 - sameX_g is deterministic, since arrivals never exceed sends — and 1_{T>g-1} is
H_g-measurable, so E[1_{T>g-1} walks_g^main] <= PHI[g] x bound + DEFV: the domination is established
before multiplying by PHI[g]. Same-X pairs are carved out and charged globally without PHI (E[1 x M_uv] <= E[M_uv]; given
the shared X, different-Y leaves are independent and same-incidence leaves take the single-leaf
marginal): at most KB129 + 1/r each (certificate SXV) plus d_batch for a defective shared batch
(inside DEFV). K2 hits and listings
drop the stopping indicator (E[1*X] <= E[X]) and are charged globally by KB(64)/KB(182) pair counting (EK2, ELIST,
whose sip terms cover same-incidence pairs), plus at most 2 x d_batch per pair with either side
defective (DEFK2, DEFL) and a defective record's own-walk verifies/lists (at most HOME_MAX each,
inside DEFK2/DEFL). ELIST's leading term is global: a 182-agreeing pair sits at exactly one home
(the 182 bits include d1, d2), so it takes KB(182) alone with no fixed-home KB(72) (6.4's other_h
keeps it, legitimately, for a specified collision's home). A walk of v visits costs 5v plus one B to
the end after an even number of visits without a match (at most floor(v/2)), so the E-rows charge
5.5 per expected visit; a K2 hit's verification (4: LDR, EOR, TST, B.NE) and a listing (4) are in EK2 and ELIST. Old-X
requests serve from each source's retained all-groups X affine store (proof.md 2.3, 5), identically
for every group, so cross-group candidates resolve with no regeneration row.

Key-path means (E-row, the same domination). A leaf takes JH or JT, 3 more instructions (`JPATH`; the sign's
one B, `SJOIN`, is charged always), only if the value its batch computes has exact index 1 or 2 or lies in W
(proof.md 3.2, 3.3). Split by its batch's validity. If the batch is valid, the computed value is the exact leaf
x(φ^t(X_k) ± Y_l); given H_g (which fixes R_g) and the other X-batches of group g, X_k = R_g + [a_k]P is fresh
uniform on G minus at most one point, hence so is the leaf point (φ and the shift are bijections of G). Exactly two
thirds of the points of G \ {O} have exact index 1 or 2 (φ acts freely on G \ {O}; its orbits {P, φP, φ²P} carry the
three indices once each, the uniformity of proof.md 3.2), and at most d_wj (r - 1) points lie in W, so P(JH or JT and
valid | ·) <= 2/3 + d_wj + 1/(r - 1) <= PJ = 2/3 + 2d_wj + 1/(r - 1), uniformly (d_wj > 1/(r - 1)). Hence
E[keypaths_g x 1_{valid} | H_g] <= s_g (3 PJ + 1/3) for the group's s_g leaves, and as 1_{T>g-1} is H_g-measurable,
E[1_{T>g-1} keypaths_g^valid] <= PHI[g] x s_g (3 PJ + 1/3) (certificate `WGJ`, inside W_g): the domination is
established before multiplying by PHI[g]. A defective batch (unconditional d_batch per batch, 3.4's fresh-Y
derivation, never conditioned on Y) takes at most 3 on each of its 6q2 leaves past the B charged always, so
E[keypaths x 1_{defective}] <= 3 d_batch per leaf of the whole run, charged globally without PHI (certificate `DEFJ`).
There is no key-path meter: JH and JT form N' and run into the tree copies ending in one B, so R_worst charges every
leaf 3 + 1/3 (certificate `JROW_W`) and 6.4 has no key-path row.

Equivalence E1-E7 (grouped with no early stop = base success). E1: same records (same kernels and
X/Y order; groups are contiguous X-ranges; cursors persist by the same rule). E2: same hop-2 slots
(cumulative cursors, same drops). E3: hash chains link same-key records (each group call reads
[read-cursor, E2) and advances the cursor, so the ranges partition each slot and every record is
inserted exactly once; no drops, buckets or capacity check). E4: same-key classes are connected
(each new record links to its most recent same-key predecessor; the uncapped walk reaches it). E5:
same coverage (every colliding pair with no intruder between its records gets its chain edge; a
chain with a useful pair has a useful chain edge as for the base adjacency; candidate sets agree up
to the list cap, else the per-prefix d_cl rows cover it). E6: same resolution. E7: same success
event (stopping only skips groups after a confirmed x, which latches; prefix defects and caps are
bounded by the prefix rows). Therefore P(success) is the grouped full-prefix bound at the same 0.39
target; q1 is re-searched least under it (certificate, 145307571).
