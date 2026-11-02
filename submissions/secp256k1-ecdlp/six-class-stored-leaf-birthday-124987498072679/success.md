# Success, section 6 of proof.md (no heuristic)

### 6.1 Model

Coins are uniform and independent, and Q is uniform. Each `Y_l` is uniform on G \ {O}; given R_g, the
`X_k = R_g + [a_k]P` are independent and uniform. So the X's are independent uniform points, independent of the Y's, and
nodes are independent; a leaf is uniform on G, and leaves of different nodes are independent. Given all generator
points, the `b_g` are still independent and uniform: for each value of b_g exactly one tuple `(a_k)` fits the points.

Γ is the set of unordered pairs of leaves from different nodes, `|Γ| = (S² - 2^72*n0²)/2`. For α = {L, L'} ∈ Γ let
`I_α = 1[L' ∈ class(L)]`; then `P(I_α = 1) = P1 = (6(r-1) + 1)/r²` and `μ = |Γ|*P1 = 0.494456934` (the run's `mu`).
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

A **node-local defect** of L (a batch defect of 3.4, a wrong class index, a label error, either meter) is a function
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
| batch defect: the sites of 3.4, D = 0, X = O | `d_batch` | < 2^-125.6 |
| wrong class index (3.3): KEY's window set W | `d_wj` | < 2^-58.47 |
| j = 1, 2 meter, Chebyshev with slack `2^42 - 6q2 - 1 - 2n0*d_wj` | `d_j` | < 2^-30 |
| fixup meter, Markov on the mean `q1*E_FIX` with the cap `2^45 - 37q2` (3.4) | `d_f` | < 2^-15.85 |
| the label-error set D' (3.3) | `d_lab` | < 2^-17.38 |
| pair slot at hop 1, and again at hop 2 (Chebyshev, `cap_p`) | `d_p`, counted twice | < 2^-19.9 each |
| d3 bucket (Chebyshev, `cap_3 = ceil(m_3) + 2^32`) | `d_3` | < 2^-30 |
| intruder (per pair), `ρ(S + 12)KB(142)` | `lam_int` | < 2^-12.6 |
| candidate, request and relay caps (Markov; a reply counts 8 records) | `d_cl`; `d_rq`; `d_rl` | < 2^-23.4; < 2^-25.9; < 2^-25.2 |
| a computational defect in any batch (whole run) | `g_def = 2^72*q1*d_batch` | < 2^-26.5 |
| two coin range rejections (whole run) | in `glob` | all of `glob` < 2^-26.5 |
| RNDR failure (whole run) | `g_rndr` = 0 under the rules' ideal coins (section 7) | - |

### 6.5 The bound

`certificate.py` picks the least q1 for which

```text
1 - e^-mu - (b1 + b2) - mu(1 + 1/r)(1/r + intruder + 2 * sum of per-leaf rows) - whole-run rows >= 0.390000001,
```

with an upper enclosure of e^-μ and exact Fractions: q1 = 145305557, where the left side is at least 0.3900000023, and
the target is 10^-9 above the rules' 0.39.
