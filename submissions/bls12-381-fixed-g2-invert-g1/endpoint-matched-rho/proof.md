# Endpoint-matched multi-walk rho in GT: complete attack and cost bound

Status: ready for independent exploratory review. No qualification, execution or
emitted score is claimed. This document states the algorithm, its success bound
and its complete classical work bound for `bls12-381-fixed-g2-invert-g1` under
`classical-group-work-v2` and the protected `target-profile.json`
(SHA-256 `c18402d08303b68d536a4c1ea2241420f04455e7e90b970cf9310d12374788a9`);
`claim.json` binds both protected files.

Submitted score: **time_log2 = 127.602539**, a rational 127602539/1000000,
rounded upward from log2(2^127.60253841) so that the declared bound is never
below the computed one.

## 1. Target and reduction

The instance is a nonidentity `Q` in G2 and a target `z` in GT. Both groups are
cyclic of prime order

```
r = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
log2 r = 254.8570894130
```

with `P0` the pinned G1 generator of `target-profile.json`. For fixed `Q` the map
`P -> e(P,Q)` is a group homomorphism `G1 -> GT`. It is injective because
`e(P,Q) = 1` with `Q` nonidentity forces `P = infinity`, and both groups have
order `r`, so it is an isomorphism. Hence

```
g := e(P0,Q)  generates GT,   g != 1,
```

and the challenge guarantees `z = e(P0,Q)^a = g^a` for a unique `a` in `[0,r)`.
Recovering `a` from `(g,z)` and returning `P = [a]P0` is admissible under
`rules.md`; bilinearity gives `e(P,Q) = e(P0,Q)^a = z` for that exact `Q`. The
output is then verified by the attack itself with one pairing evaluation. If
`z = 1` then `a = 0` and `P = infinity`, which the acceptance rule permits
exactly in that case.

The remaining problem is therefore a discrete logarithm in the prime-order group
GT. No representation shortcut is used: the attack is a generic birthday search,
priced in GT group operations, each of which is one dense Fp12 multiplication,
the reference operation of this track.

## 2. Algorithm

Parameters, all fixed before the instance is seen:

```
k = 2^20                       number of walks
L = 240218821685709510680521902710375   steps per walk   (log2 L = 107.5660482303)
N = k*L = 251887691167914535871338934656434176000   (log2 N = 127.5660482303)
N / sqrt(r) = 1.1000000
```

The canonical encoding of a GT element is the twelve Fp coefficients of
`target-profile.json` (increasing degree, each a 96-hex-character integer below
`p`). Write `c11` for the integer coefficient of `w^11` and define the class

```
e(x) = 1 + (c11 mod 32)        in {1,...,32}
```

which depends only on the group element, not on the walk. The iteration is

```
step(x) = x * g^{e(x)}
```

a single dense Fp12 multiplication by a table entry.

```text
1.  Compute g = e(P0,Q) with the pinned pairing.        (one pairing)
2.  Build the multiplier table T[j] = g^j for j = 1..32. (31 multiplications)
3.  For i = 1..k, sample u_i, v_i uniform in [0,r) and set x_i = g^{u_i} z^{v_i}.
4.  For each i, iterate x <- step(x) exactly L times, keeping only the current
    value; store the final endpoint y_i together with (i, u_i, v_i).
5.  Sort the k records by the canonical endpoint encoding and scan for equal
    adjacent endpoints. Every equal group yields a merged pair (i,j).
6.  For the first pair with v_i != v_j, re-run walks i and j from their starts,
    accumulating the class sums c_i and c_j, and locate the first step at which
    the two runs coincide.
7.  Set a = (u_i + c_i - u_j - c_j) * (v_j - v_i)^{-1} mod r, return P = [a]P0,
    and accept only if e(P,Q) = z.
```

The walk is a deterministic function of the current group element, so two walks
that ever occupy the same element are identical from that step onward. Therefore

* equal endpoints are necessary and sufficient for the two paths to have
  intersected: no false positives, because equal final elements mean the paths
  met at or before step `L`; no false negatives, because a meeting at step
  `m <= L` forces equality at every later step including step `L`.

Detection is exact; only the event "some pair meets" is probabilistic.

## 3. Correctness of the recovered scalar

Walk `i` after `c_i` steps holds `g^{u_i + c_i} z^{v_i}`. At the meeting point
the two values are equal:

```
g^{u_i + c_i} z^{v_i} = g^{u_j + c_j} z^{v_j}
g^{u_i + c_i - u_j - c_j} = z^{v_j - v_i} = g^{a (v_j - v_i)}
u_i + c_i - u_j - c_j = a (v_j - v_i)  (mod r)
```

`v_j - v_i` is invertible modulo the prime `r` unless `v_i = v_j`; that single
excluded case has probability `1/r` and step 6 skips such pairs. The equation
then determines `a` uniquely, and step 7 checks the result with the pairing, so
a returned point is always a correct answer. Degenerate pairs can only cause a
failed run, never a wrong output. The case `z = 1` is covered: the equation
gives `a (v_j - v_i) = 0`, hence `a = 0` and `P = infinity`, as required.

## 4. Success probability

The only probabilistic event is the existence of a meeting between two of the
`k` walks. Under the standard random-mapping model for r-adding walks, a pair of
walks of length `L` meets with probability `1 - exp(-L^2/r)` for `L^2/r << 1`,
and the `C(k,2)` pairs are treated as independent trials, giving

```
E[merges] = C(k,2) * L^2 / r  ~=  (k L)^2 / (2 r) = 1.1^2 / 2 = 0.605
P(at least one merge) ~= 1 - exp(-0.605) = 0.45393
```

The declared bound is the conservative **success_probability = 0.39**, which the
model supports with a 16 per cent relative margin. Three model errors are
absorbed by that margin and are disclosed in section 8: the walk is not literally
a random mapping; the Poisson approximation to the number of merges is not a
proof; and pairs are not independent. Each of these affects the constant in
`1 - exp(-(kL)^2/(2r))`, not the exponent, and the design point `N = 1.10 sqrt(r)`
sits 10.6 per cent above the `0.9943 sqrt(r)` that the same model needs for 0.39.

Within-walk self-intersections are irrelevant here: `L^2/r = 2^-47.7`, and such
an event would only shorten a path, never produce a wrong answer.

The work bound does not depend on the coins at all. `k` walks of exactly `L`
steps are always executed, so the total number of GT multiplications is the
constant `N`, not an expectation; the randomness affects only whether the run
succeeds. This is why the score is a deterministic function of the fixed
parameters and cannot be improved by a lucky coin.

## 5. Complete work ledger

The track price is one dense Fp12 multiplication = 1 target-group-operation
equivalent, `C = 20000` scalar AArch64 instruction equivalents. Measured
`blst_fp12_mul` is 22,996 instructions, so a multiplication is charged as one
reference operation and all surrounding work is charged as instructions divided
by 20,000. Nothing is charged twice: the multiplication's own internal loads,
stores and reduction are inside the reference operation, and only the caller-side
work is added.

| Item | Work |
| --- | ---: |
| Walk steps: `N` dense Fp12 multiplications | 2^127.5660482 |
| Caller-side per-step work: 512 instructions/step | `N * 512/20000` = 2^122.2783 |
| Pairing `g = e(P0,Q)` (charged 4x a blst pairing) | 2^14 |
| Multiplier table `g^1..g^32` | 31 |
| `k` random starts `g^{u_i} z^{v_i}`: 2 scalar muls x 383 ops each | 803,209,216 = 2^29.5812 |
| Endpoint stores: `k` records of 644 bytes | 2^12 |
| Sort of `k` endpoints: `k log2 k` comparisons x ~200 instructions | 2^17.6781 |
| Collision resolution: 16 walk re-runs of `L` steps | 2^111.6025 |
| Randomness: `k` pairs of 255-bit words | 2^8.7123 |
| Modular inversion and output construction | 1 |
| Attack-side verification `e(P,Q) = z` (one pairing) | 2^14 |

Per-step caller work is itemised as: one 64-bit load and mask of the low limb of
`c11` plus the increment for `e(x)` (4 instructions), multiplier address
computation (4), loop counter and branch (4), and the 576-byte reload/store of
the state around the multiplication (72 loads, 72 stores and 144 address
operations, 288). That totals about 300 instructions; 512 is charged, a 1.7x
margin, and it is 2.2 per cent of the 22,996-instruction multiplication it
surrounds.

Summing the table:

```
W = N * (1 + 512/20000) + 2^111.6025 + 2^29.5812 + 2^18 + 2^14 + 2^12 + 2^9 + 31
  = 2^127.60253841
```

The walk term is larger than every other term by at least 2^16, so the total is
dominated by the multiplications and the smaller lines are listed only for
completeness. Dividing by the fixed price is already done: the ledger is in
target-group-operation equivalents, so

```
time_log2 = log2(W) = 127.60253841  ->  submitted 127.602539 (rounded upward)
```

Rounding upward is deliberate: the submitted number is never smaller than the
computed bound.

## 6. Memory, preprocessing and advice

Peak storage is the endpoint table, `k = 2^20` records of 1024 bytes each
(576 bytes of canonical Fp12 endpoint, 64 bytes for `u_i,v_i`, 4 bytes of walk
index, and slack for padding and allocator headers), plus one sorted copy of the
same size:

```
peak = 2 * 2^20 * 1024 bytes = 2^31 bytes
memory_log2_bytes = 31
```

The walk itself keeps one 576-byte state, the 32-entry multiplier table (18,432
bytes) and the 64-byte random pair per walk, all far below the table. Every
record is retained until the sort completes; nothing is streamed away.

Preprocessing is the start construction plus the pairing and the table:
`803,209,216 + 2^14 + 31` equivalents, `preprocessing_log2 = 29.581231`. It is a
subset of the total above and is not charged a second time. There is no stored
advice: the multiplier table and the starts are derived from the public input at
run time, so `nonuniform_advice_log2_bytes = 0`. No result of any previous
instance is reused, and nothing is amortised over future targets.

## 7. Why the endpoint decomposition is used

The usual distinguished-point variant terminates each walk at a sparse set and
then charges a random walk length with mean `1/theta` plus a table sized
`N*theta`. That makes the work bound an expectation over walk lengths and forces
the collision-detection analysis to carry the distinguished-point density. Here
walks have a fixed length and collisions are found by comparing final elements,
so the number of group operations is the constant `N`, the table is exactly `k`
records, and detection is exact rather than probabilistic. The cost of that
choice is one sort of `2^20` records, which is `2^17.7` equivalents, and the
resolution re-runs, `2^111.6` equivalents; both are more than 2^16 times smaller
than the walk term. Choosing `k = 2^20` keeps the table at 2 GiB, so the memory
claim is a commodity machine rather than an unphysical one.

## 8. Heuristics, support and limitations

Disclosed heuristics:

1. *Random-mapping model.* The r-adding walk with a 32-class alphabet is assumed
   to behave like a random mapping on GT for the purpose of the birthday
   estimate. This is the standard model for Pollard rho and r-adding walks
   (van Oorschot-Wiener 1999; Teske's analyses of r-adding walks), and the class
   function depends on a single canonical coefficient of the state, which is the
   usual construction. It is a heuristic, not a theorem; the 16 per cent margin
   over 0.39 and the 10.6 per cent margin in `N` are the allowance for it.
2. *Poisson approximation.* `P(no merge) ~= exp(-E[merges])` assumes the pair
   events are close to independent. The true count is a sum of weakly dependent
   indicators; the deviation is a constant-factor effect in the exponent of the
   success expression.
3. *Uniformity of `u_i,v_i`.* Sampling from `[0,r)` by rejection sampling makes
   the starts uniform over GT; the analysis uses only that, and a repeated
   `(u,v)` pair has probability `~2^40/r`, negligible.

Limitations stated plainly:

* No production-size execution exists; the bound is analytic, which the cost
  model explicitly permits. Small experiments are not offered as evidence for
  the production complexity.
* The bound is a generic birthday bound for GT and is not a security estimate
  for BLS12-381; it makes no claim about index-calculus or number-field methods
  over Fp12, which are outside this submission.
* The bound is not claimed to be tight. It is a stated upper bound with a
  deliberate success margin and upward rounding.
* The 2^31-byte table and the `2^20` records are required by this decomposition;
  they are reported as supplementary information and carry no score weight.

## 9. Accounting summary

| Field | Value |
| --- | ---: |
| time_log2 | 127.602539 |
| time_unit | target-group-operations (dense Fp12 multiplication, C = 20,000) |
| memory_log2_bytes | 31 |
| preprocessing_log2 | 29.581231 (included in the total) |
| success_probability | 0.39 (model gives 0.45393) |
| nonuniform_advice_log2_bytes | 0 |
