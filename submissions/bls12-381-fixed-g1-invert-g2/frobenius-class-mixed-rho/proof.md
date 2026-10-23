# Frobenius-class mixed add/square rho on GT

## 1. Reduction

The challenge fixes a nonidentity `P` in G1 and an independently uniform target
`z` in GT. Put

```
g = e(P, Q0)        (one pairing, P is the supplied argument, Q0 the fixed G2 generator)
h = z
```

`e(P, Q0)` is nondegenerate because `P != O` and `Q0` generates G2, so `g`
generates the prime-order group GT. Hence there is a unique `x` in `[0, r)` with
`h = g^x`, and `Q = [x]Q0` satisfies

```
e(P, Q) = e(P, Q0)^x = g^x = z.
```

`Q` is returned canonically; `Q = O` arises exactly when `x = 0`, i.e. when
`z = 1`. Nothing else is needed from the challenge: `P` is used only through the
single pairing `e(P, Q0)`, and no auxiliary powers, advice or oracle are assumed.

Everything below is a discrete logarithm in GT with known base `g` and target
`h`, so the reference unit is the dense Fp12 multiplication with
`C = 20,000` scalar AArch64 instruction equivalents.

## 2. Cheap maps on GT

GT = mu_r is the group of `r`-th roots of unity inside Fp12, and
`r | p^4 - p^2 + 1`, so

```
p^6 + 1 = (p^2 + 1)(p^4 - p^2 + 1) = 0  (mod r).
```

Therefore for every `Y` in GT

```
Y^{p^6} = Y^{-1}.
```

In the tower `Fp12 = Fp6[w]/(w^2 - v)` the map `Y -> Y^{p^6}` is the nontrivial
automorphism of Fp12 over Fp6, i.e. negation of the `w`-half: 6 Fp2 negations,
12 Fp negations, about 24 instructions. **Conjugation is free.**

Three maps generate the emission set:

| map | action on exponents | cost |
| --- | --- | --- |
| walk step `X -> X * S_j` | `a + b x -> a + c_j + (b + d_j) x` | 1.0 C (dense Fp12 mul) |
| cyclotomic squaring `X -> X^2` | `e -> 2 e` | 0.53 C |
| Frobenius `X -> X^p` | `e -> p e` | 0.36 C |
| conjugation `X -> X^{p^6}` | `e -> -e` | free |

Costs, from the pinned measurement `383` instructions per `blst` Fp
multiplication (`COST_MODEL.md`):

* cyclotomic (Granger-Scott) squaring for `k = 12`: 3 Fp2 squarings and 3 Fp2
  multiplications, priced here as 9 Fp2 multiplications = 27 Fp muls
  = `27 * 383 = 10341` instructions, plus additions and memory traffic, call it
  `10541` instructions = `0.527 C`; we use **0.53 C**.
* Frobenius `pi^k` (`k` not a multiple of 6): 6 Fp2 coefficients, the `i = 0`
  coefficient of the Fp6 half has multiplier 1, so 5 Fp2 constant
  multiplications = 15 Fp muls = `5745` instructions, plus 12 Fp2 loads/stores
  and 12 Fp conjugations, call it `7100` instructions = `0.355 C`; we use
  **0.36 C**.
* table insert or lookup of one element (64-bit hash over the 12 Fp
  coefficients, one probe, 48-byte payload, 64-byte exponent pair): about
  `130` instructions = `0.0065 C`; we use **0.007 C**.

All constants used by these maps (Frobenius coefficients, `p^k mod r`) are
public profile constants and cost nothing to construct.

## 3. Block emission

The walk keeps a table of `(X_i, a_i, b_i)`. A block is

1. one dense multiplication `X_i = X_{i-1} * S_{j_i}`, where the 32 step
   elements `S_j = g^{c_j} h^{d_j}` are precomputed once (32 C total), and
   `(a_i, b_i) = (a_{i-1} + c_{j_i}, b_{i-1} + d_{j_i})`;
2. `D` cyclotomic squarings giving `Y_d = X_i^{2^d}` for `d = 0..D`;
3. for each `Y_d`, five Frobenius maps `Y_d^{p^k}`, `k = 1..5`.

The emitted set of one block is

```
{ Y_d^{eps p^k} : d = 0..D, k = 0..5, eps = +-1 } ,   |set| = K = 12 (D + 1).
```

The `eps = -1` elements are conjugates of the `eps = +1` elements and cost
nothing. Note `Y_d^{p^{k+6}} = Y_d^{-p^k}`, so the twelve elements per level are
exactly `{Y_d^{p^j} : j = 0..11}`: the Frobenius orbit of `Y_d`, whose size is 12
because the Frobenius has order 12 on mu_r.

Cost of one block, `D = 31`:

```
1.0                      walk multiplication
0.53 * 31 = 16.43        cyclotomic squarings
0.36 * 5 * 32 = 57.60    Frobenius maps (5 per level, 32 levels)
0.007 * 384 = 2.69       table inserts
------------------------------
total 77.72 C            for K = 384 elements, i.e. 0.2024 C per element
```

The `D`-dependence is `(2.884 + 2.414 D)/(12 (D + 1))`, decreasing to `0.2012`
as `D -> infinity`; `D = 31` is within 0.6% of the limit and keeps the per-block
scratch small.

## 4. Collision algebra

Write `m_u = eps_u p^{k_u} 2^{d_u} mod r` for the public multiplier of emitted
element `u` of block `i`. Element `u` is the group element

```
V_u = g^{m_u a_i} h^{m_u b_i} = g^{m_u (a_i + b_i x)}.
```

If `V_u = V_v` for two elements of **different** blocks `i != j`, then

```
m_u (a_i + b_i x) = m_v (a_j + b_j x)  (mod r)
```

and, provided `m_u b_i - m_v b_j != 0 (mod r)`,

```
x = (m_v a_j - m_u a_i) * (m_u b_i - m_v b_j)^{-1}  (mod r).
```

One modular inverse and three multiplications mod r, negligible against the
search. The candidate is checked by the pairing `e(P, [x]Q0) = z` (about `2^7`
C); a failed check means the collision was degenerate, and the next collision is
used. Degeneracy needs `m_u b_i = m_v b_j (mod r)`, probability about `1/r` per
collision, and there is on average `0.6` collision per attempt, so this costs
nothing.

Within one block all pairs are *deterministic* relations (every emitted element
is a fixed power of the same `X_i`); two of them are equal only if `X_i = 1`,
which would immediately give `x = -a_i / b_i`. They are excluded from the count
and from the solver. This is the point where earlier candidates in this family
went wrong: they counted *ratio classes* (potential equalities) rather than
*visited equalities*. Here every emitted element is inserted into the table and
every table hit is a real value equality, so the counted set and the detected
set are identical.

## 5. Birthday accounting

Let `E = K N` be the total number of emitted elements over `N` blocks. Cross
block unordered pairs number `C(E,2) - N C(K,2) = E^2/2 (1 - 1/N)`.

Model (standard rho heuristic, supported below): for `i != j` the walk values
`X_i` and `X_j` are independent uniform elements of GT, the multipliers
`m_u, m_v` are fixed public constants, and `Y -> Y^m` is a bijection of GT for
every `m != 0 (mod r)`. Hence `(V_u, V_v)` is uniform on `GT x GT` and

```
Pr[V_u = V_v] = 1/r        for every cross-block pair.
```

By linearity, the expected number of collisions is `lambda = E^2 / (2 r)`, and
the Poisson approximation gives

```
Pr[at least one collision] = 1 - exp(-E^2 / (2 r)).
```

Choose success probability `0.45`:

```
lambda = -ln(0.55) = 0.59784,   E = sqrt(2 * 0.59784 * r) = 2^127.5574
```

using `log2 r = 254.8568`. Then `N = E / 384 = 2^118.972` blocks and

```
W = N * 77.72 C = 2^118.972 * 2^6.280 = 2^125.2526 C.
```

Adding the one-off work (initial pairing `e(P,Q0)`, 32 step elements, final
`[x]Q0`, verifying pairing, modular inverse) worth about `2^8` C changes nothing
at the seventh binary digit. The submitted bound is

```
time_log2 = 125.30.
```

The attempt is charged once. Amplification to higher success would cost
`1/0.45` times as much and is not claimed or needed; `0.45 >= 0.39` satisfies the
track minimum with margin for model error.

Peak storage: `2^127.5574` table entries of 176 bytes (12 Fp coefficients,
`a_i`, `b_i`, block and multiplier index) gives `log2 = 135.0` bytes; with code
and scratch we report `2^135.1` bytes. Memory has no scalar contribution.

Preprocessing: one pairing and 32 step elements, `2^8.0` target-group-operation
equivalents, already included in `W`. Advice: 0 bytes.

## 6. Heuristics, support and limitations

1. **Uniformity of the walk.** Steps are chosen uniformly from a fixed set of 32
   elements `g^{c_j} h^{d_j}` with `(c_j, d_j)` uniform in `[1, 2^16)^2`; the
   exponent pair after `i` steps is a sum of `i` independent uniform increments,
   which is uniform on `Z_r x Z_r` after about `r^{1/2}` steps by the standard
   mixing heuristic. The walk is not required to be a function of the point
   because the whole table is stored, so no functional-graph bias enters.
2. **Independence across blocks.** `X_j = X_i * (product of steps i+1..j)` and
   the later steps are independent of `X_i`, so `X_i` and `X_j` are independent
   uniform. This is the same heuristic used by every parallel-collision-search
   claim for this track; it is the weakest link and it is what experiment.py
   tests at toy scale.
3. **Cheap maps.** The Frobenius and cyclotomic-squaring prices are analytic
   instruction counts derived from the pinned `383`-instruction Fp
   multiplication, not from a new hardware measurement; they are the second
   weakest link. If a Frobenius were priced as a full dense multiplication the
   bound would rise by about 0.9 bits, which is a cost-model dispute, not a
   correctness failure.
4. **Free conjugation.** `p^6 = -1 (mod r)` is exact and is checked numerically
   by experiment.py on the real profile constants.
5. **Toy-scale evidence.** experiment.py runs the same multiplier set, the same
   block structure and the same solver in `Z_q` for a 30-bit prime `q`, and
   compares the recovered-exponent rate and the collision count with
   `1 - exp(-E^2/(2q))` and `E^2/(2q)`. It supports the combinatorial model only;
   it says nothing about production-size behaviour or field-arithmetic cost.
6. **Not claimed.** Novelty, practical executability at `2^135` bytes, or any
   improvement of the generic-group exponent. The gain over plain
   distinguished-point rho is a constant factor: the emission gives 12 elements
   per level and `2(D+1)` elements per walk multiplication, worth about 2.3 bits
   in the reference unit.

## 7. Accounting summary

```
log2 r                 254.8568
success probability    0.45      (minimum required 0.39)
E (emitted elements)   2^127.5574
K (elements per block) 384
N (blocks)             2^118.972
cost per block         77.72 C   (1 walk mul + 31 squarings + 160 Frobenius + 384 table ops)
log2 W                 125.2526
submitted time_log2    125.30
memory_log2_bytes      135.1
preprocessing_log2     8.0       (included in W)
advice bytes           0
```
