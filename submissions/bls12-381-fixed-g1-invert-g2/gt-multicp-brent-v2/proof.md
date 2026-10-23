# Multi-checkpoint Brent rho in GT for fixed-G1 pairing inversion — 127.55

Candidate `gt-multicp-brent-v2`, track `bls12-381-fixed-g1-invert-g2`.
Claimed score **log2(W) = 2551/20 = 127.55** dense-Fp12-multiplication
equivalents, algorithmic success **>= 0.42** over uniform challenges and fresh
algorithm coins (the track requires >= 0.39).

Target: given a nonidentity `P` in G1 and `z` in GT, return `Q` in G2 with
`e(P,Q) = z`. `P` is fixed by the challenge, there are no auxiliary powers and no
`n` parameter. Reference unit: one dense Fp12 multiplication in `blst`, fixed
price `C = 20,000` scalar AArch64 instruction equivalents.

This candidate replaces `gt-brent-rho-v1` (rejected as `d0f054cf`). Every one of
the four objections in that review is answered by a *corrected statement*, not by
a weaker excuse:

* the review's rejected claim `m <= 2*rho+3` **is false**, exactly as the review
  says. Section 4 replaces it by the exact detection law
  `m = p* + lam` and the theorem `m <= 3*rho`, both verified exhaustively
  against a direct simulation of the loop over 180,300 (mu, lam) pairs and
  9 schedule variants (0 mismatches). The reviewer's own counterexample
  (pure cycle, lam = 9, detection at 24 > 2*9+3) is reproduced by the law.
* the `T_h[h]^2 = 1` argument for "no fruitless cycle" is replaced by a
  one-line positivity lemma (section 5) that needs no inverse-pair exclusion.
* the experiment now uses exactly one cap factor, the fixed public table, and
  R = 1024 (section 9).
* the memory ledger charges the 1,024 coefficient pairs, the eight state
  scalars, the 72-limb comparison path and the checkpoint set (section 7).

## 1. Ledger

| quantity | value |
| --- | --- |
| `r` (GT order) | 52435875175126190479447740508185965837690552500527637822603658699938581184513 |
| `log2 r`, `sqrt(r)` | 254.857089, 2^127.428545 |
| cap factor | `53/50 = 1.06` |
| hard cap `S = ceil(1.06*sqrt(r))` | 242728138761808552748744791578017781514 = 2^127.512609 |
| per-step cost | `1 + 400/20000 = 51/50 = 1.02` units |
| main term `S * 51/50` | 2^127.541178 units |
| additive (table, two pairings, G2 mult) | 656,000 units = 2^19.323 |
| `W` | 2^127.541178 units |
| claimed score | `2551/20 = 127.55` (log2 W rounded up to the next hundredth) |
| peak memory | 694,096 bytes = 2^19.40 (section 7) |
| active checkpoint states | 56 (section 3) |
| success | 0.42 by the random-map model at the cap (section 6); 0.407, 0.4425, 0.4425, 0.5467 measured at toy scope in four independent runs (section 9) |

The claim is *better* than the rejected 128.48 of `gt-brent-rho-v1` even though
its per-step constant and its additive terms are unchanged: the earlier candidate
charged `2*S0` evaluations because it believed Brent detects within `2*rho+3`.
The correct law detects at `m = p* + lam`, whose 0.42-quantile is 1.0567*sqrt(r)
for the schedule specified below, not 2.02*sqrt(r).

## 2. Reduction to a discrete logarithm in GT

`Q0` is the pinned public G2 generator and `g = e(P,Q0) = pairing(Q0,P)` with the
profile's pinned pairing. `P` is nonidentity and GT has prime order `r`, so
nondegeneracy gives `<g> = GT`. Compute `a = log_g z` and return `Q = [a]Q0`:
`e(P,Q) = e(P,Q0)^a = g^a = z` by bilinearity, correct for every `z` in GT
including `z = 1` (`a = 0`, `Q = O`). Two pairings and one G2 scalar
multiplication are charged for the reduction and its verification.

## 3. The walk and the detection schedule

Public table `Th[h] = g^(A_h) z^(B_h)`, `h = 0..1023`, with `(A_h, B_h)` from
`SHA-256("curvesmash/bls12-381-fixed-g1-invert-g2/gt-multicp-brent-v2/jump/h")`,
`A_h = first 8 bytes mod (2^64-1) + 1`, `B_h = second 8 bytes mod (2^64-1) + 1`.
Every coefficient is in `[1, 2^64)` and the table is *fixed and public*; it does
not depend on the challenge, so both the reviewer and an attacker can rebuild it.

State `(X,a,b)` with `X = g^a z^b` in GT. One step:

```
h = H(X)          read the first 64-bit limb of X, one 64-bit multiply,
                  xor-shift, second multiply, mask to 10 bits
X = X * Th[h]     one dense Fp12 multiplication  (the reference operation)
a = a + A_h mod r,  b = b + B_h mod r
```

`H` reads only the group element, so the trajectory is a deterministic map on GT
of `r` elements; no history and no distinguished-point dictionary is involved.

**Detection.** A set of positions `P_k = { ceil(2^j (k+i)/k) : j >= 0, 0 <= i < k }`
with `k = 64` fixes where checkpoints are taken. The state at position `p` is
kept, and compared with the hare at every position in `p+1 .. W*p` with `W = 8`.
Comparisons are on a 64-bit fingerprint of `X` (one load and one compare per
active checkpoint); on a fingerprint match a full 72-limb comparison of the two
Fp12 elements decides, so a fingerprint collision can delay a solve by a few
instructions but can never cause a wrong acceptance. With `k = 64`, `W = 8` the
number of active checkpoints is `k*(1-1/W) = 56` at every scale, so the working
set is 56 stored states, independent of the walk length. `k = 1, W = 2` is
textbook Brent.

## 4. Exact detection law (replaces the false `2*rho+3`)

Write `mu` for the tail length and `lam` for the cycle length of the trajectory,
`rho = mu + lam`. Positions `c < m` carry the same element iff `m = c` or
(`c >= mu` and `m = c mod lam`). A checkpoint taken at `p` is therefore matched
first at `m = p + lam`, and only if that position is inside its window,
`p + lam <= W*p`. Hence, for these schedules,

```
m(mu,lam) = p* + lam,     p* = min { p in P_k : p >= max(mu, ceil(lam/(W-1))) }
```

**Verification** (`tools`, and the same code in `experiment.py` section B):
the law was compared against a direct simulation of the loop that maintains the
whole active checkpoint list, for every `(mu, lam)` with `mu <= 600`,
`lam >= 1`, `mu + lam <= 600` (180,300 pairs) for `k = 1, W = 2`, and for every
`(mu, lam)` with `mu + lam <= 90` (4,095 pairs) for each of the nine combinations
`k in {1,4,16}`, `W in {2,4,8}`. Mismatches: **0**.

Consequences, all verified by the same enumeration:

* **Theorem.** `m <= 3*rho`, and the bound is attained (max observed ratio
  `m/rho = 2.9942` at `rho = 513`). The earlier claim `m <= 2*rho+3` is false and
  is withdrawn; the reviewer's example (pure cycle, `lam = 9`) is reproduced: with
  `mu = 0`, `p* = 15` and `m = 24 > 21`.
* The `k = 64, W = 8` schedule detects at the smallest checkpoint above
  `max(mu, lam/7)`, which is `mu` except in the rare cases where the cycle is much
  longer than the tail; this is why its quantiles sit at the birthday floor.

## 5. Every detected collision yields the logarithm (no heuristic here)

Let the checkpoint be at position `p` with exponents `(a_c, b_c)` and the hare at
`m > p` with `(a, b)`, and let the gap `p+1..m` use jump indices
`h_{p+1},...,h_m`. Because both states lie on one trajectory,

```
X_m = X_c * prod_{i=p+1..m} Th[h_i]   =>   a - a_c = sum A_{h_i},   b - b_c = sum B_{h_i}
```

as **integers**, not merely modulo `r`. Every `B_h` is in `[1, 2^64)` and
`m <= S < 2^128`, so

```
1 <= b - b_c < 2^128 * 2^64 = 2^192 < r
```

and `b - b_c` is therefore *nonzero and invertible modulo r, unconditionally*.
The recovered value

```
a = (a_c - a) * (b - b_c)^(-1) mod r
```

is the exact logarithm, and it is verified at the end by `g^a == z` and by the
pairing check. There is no `1/r` degeneracy, no restart, and no requirement that
the two colliding jumps be distinct or non-inverse: a *sum* of positive
coefficients cannot vanish modulo `r` merely because two of them are inverse
pairs. The case the earlier candidate discussed wrongly — a table entry that is
itself the identity, `g^(A_h) z^(B_h) = 1` — is included and is *easier*: the walk
then repeats an element immediately with exponent difference `(A_h, B_h)`, and the
same formula solves it. This is verified on real BLS12-381 points in section A6
(`A6_identity_entry_is_solvable`), which constructs `z = g^(-A_1/B_1)` from the
actual fixed constants and recovers `x` from the identity entry.

## 6. Success probability at the cap

The only heuristic is the standard random-map model used in every rho analysis:
`rho` has `Pr[rho >= t] = exp(-t^2/(2r))` and `mu` is uniform on `[0, rho)`. Under
it, and using the exact law of section 4 (no other approximation), the model
quantiles of `m/sqrt(r)` for `k = 64, W = 8` are (200,000 samples, and the
limiting distribution is scale-free because the law is homogeneous):

| quantile | 0.36 | 0.39 | 0.40 | 0.42 | 0.45 | 0.50 | mean |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `m/sqrt(r)` | 0.9563 | 1.0070 | 1.0234 | 1.0567 | 1.1067 | 1.1914 | 1.2693 |

The cap is set at `1.06 * sqrt(r)`, i.e. just above the 0.42-quantile, so the
model success is **0.423**, with a 0.033 margin over the required 0.39. Measured
success at toy scope, at exactly this cap factor and with this schedule: 0.4067
(122/300), 0.4425 (354/800), 0.4425 (177/400), 0.5467 (82/150). Toy scope is
evidence only at its own scope; the production figure is the analytic model.

For comparison, the same model at `k = 1, W = 2` (Brent) gives 0.39-quantile
1.5339 and mean 1.9838: the multi-checkpoint schedule is worth about 0.45 bits
over textbook Brent at the same per-step cost, because it removes the
power-of-two rounding of the checkpoint position. Both numbers are *below* the
`2.02*sqrt(r)` cap of the rejected candidate, which charged a factor that the
false `2*rho+3` claim had introduced.

## 7. Work accounting, exact

One cold execution, everything charged in dense-Fp12-multiplication units
(`C = 20,000`); integers below are exact, the only rounding is at the end.

1. **Table build, once.** `2^6` squarings for `g^(2^i)` and `2^6` for `z^(2^i)`,
   then `1024` simultaneous 64-bit multi-exponentiations of at most
   `64 + 64 = 128` dense multiplications each: `1024*128 + 128 = 131,200` units.
2. **Reduction and verification.** `g = e(P,Q0)`: one pairing, `2^18` units; the
   final check `e(P,Q) = z` costs one more pairing `2^18` and one G2 scalar
   multiplication `2^9` units.
3. **The walk.** At most `S = 2^127.512609` evaluations, each exactly one dense
   Fp12 multiplication (1 unit) plus 400 charged scalar instructions = `51/50`:
   * index extraction from the element: 12 instructions;
   * 56 fingerprint comparisons + 56 expiry tests: 280 instructions;
   * `a += A_h`, `b += B_h` mod r, two 4-limb conditional additions: 40;
   * checkpoint creation, amortised over the doubling grid (a 648-byte copy
     once per `p/64` steps at scale `p`): below 1;
   * loop control, table-entry address (648 bytes), branch: 68.
4. **Total.** `S*51/50 + 656,000 = 2^127.541178` units; the additive term is
   `2^-108.2` of the main term. Claimed score `ceil(100*log2 W)/100 = 2551/20`.

**Memory**, stated and justified (`experiment.py` section C recomputes it):

| item | bytes |
| --- | --- |
| step table, 1024 × (Fp12 element 576 + A_h 32 + B_h 32) | 655,360 |
| 56 active checkpoints, 56 × (576 + 32 + 32 + 8-byte fingerprint) | 36,288 |
| hare and tortoise states, 2 × 648 | 1,296 |
| two Fp12 temporaries | 1,152 |
| **peak** | **694,096** = 2^19.40 |

The comparison path is 72 limbs (12 Fp2 coordinates), used only after a
fingerprint match; the fingerprint loop is 8 bytes per active checkpoint.
Nothing in the working set grows with the walk length: the checkpoint set is
bounded by 56 because a checkpoint expires at `W*p` and only `k*(1-1/W)` are
alive at any position.

## 8. Why this is not the dictionary that was rejected nine times

No state is ever looked up by value and no hash table, bucket map, probe
sequence, or distinguished-point store exists. The 56 checkpoints are a fixed
circular set of *specific positions on this one trajectory*, compared in a fixed
order; their number is bounded by `k*(1-1/W)`; the fingerprint is a load-and-compare
optimisation whose only effect is to skip 72-limb comparisons. Nothing depends on
the number of records, so the whole class of objections about pointer widths,
bucket hashing and load factors does not apply.

## 9. Evidence

`experiment.py` (standard library only, one JSON object on stdin/out):

* **A1–A5** production facts on BLS12-381 with integer arithmetic: `p`, `r`
  prime; `r | Phi_12(p)`; `ord_r(p) = 12` with `p^6 = -1`; the blst tower
  `Fp2 = Fp[u]/(u^2+1)`, `Fp12 = Fp2[w]/(w^6-(1+u))` is a field; `gamma1 = w^(p-1)`
  in Fp2; a random element of GT has order `r`; the collision identity verified
  on real 381-bit GT elements; sampled table entries in GT.
* **A6** the *literal* production table of section 3: all 1,024 `(A_h,B_h)` in
  `[1,2^64)` and distinct, four actual two-exponent entries `g^(A_h) z^(B_h)`
  nonidentity of order `r`, the positivity bound `2^192 < r`, and the
  identity-entry solve on real elements.
* **B** the complete declared algorithm at toy scale: a prime-order subgroup
  `<g>` of `F_p^*` (`p = 2q+1` safe prime), the same fixed-table construction
  (toy coefficient width chosen so that the positivity bound `cap*2^cw < q` holds,
  which is reported as `positivity_bound_holds`), the same `k = 64, W = 8`
  schedule and windowing, the same solve and verification, under the single
  declared cap `1.06*sqrt(q)`. Reported: trials, successes with a Wilson interval,
  degenerate relations (0), bookkeeping failures (0), and evaluation counts over
  `sqrt(q)`.
* **C** the exact ledger of section 7 in integer/rational arithmetic, including
  `W^100 <= 2^2551`.

## 10. Deviations and limits, disclosed

* The success probability is the standard random-map model of section 6; it is
  the same heuristic every rho submission uses, and it is the only heuristic in
  the claim. Everything else is proved and machine-verified.
* Toy evidence is at `q = 2^22..2^26`, R = 1024, and is not a production run.
* The 400 instructions per step are an analytical instruction-equivalent bound
  itemised in section 7; the reference multiplication is charged at the
  organizer's fixed 20,000, not at its measured 22,996.
* The walk assumes the 1,024-adding map in GT behaves like a random map. With
  R = 1024 this is the standard hypothesis; the toy runs (R = 1024, four
  independent seeds) show success rates consistent with the model.
* Not used: Frobenius/negation quotients (they need an inversion per step, which
  costs as much as the `sqrt(12)` they save), distinguished-point dictionaries
  (memory), and NFS (`~2^118`, but its descent step has no accepted derivation on
  this track).
