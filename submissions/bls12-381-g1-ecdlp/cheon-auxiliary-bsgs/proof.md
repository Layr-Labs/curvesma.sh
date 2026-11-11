# bls12-381-g1-ecdlp: Cheon auxiliary-input attack with BSGS and b-adic scalar tables

Status: ready for independent review; no qualification or emitted score.
This is an analytic cryptanalysis argument, not a claim of practical execution.
The normative inputs/outputs are those of the track's `target-profile.json`
and `classical-group-work-v2.json`. No restricted challenge distribution,
advice, random-oracle premise or quantum operation is used.

## 1. The claim

Let `r` be the BLS12-381 scalar-field prime (order of G1), `P0` the G1
generator from the target profile, and `s` uniform in `[1, r)`. The input
is `T = [s]P0` and powers `A_i = [s^i mod r]P0` for `i = 1..n`. The output
is `s`. The attack below uses only `A_1 = T` and `A_d` where `d = n`; it
may use a subset of the powers, as the rules allow.

For every `n` with `n | (r-1)` and `1 <= n < r`, the attack recovers `s`
deterministically (success probability 1, hence `>= 0.39`) with total work
`W(n)` target-group-operation equivalents given in Section 6. At the
scoring value `n = 2^20` (which divides `r-1`, see Section 2),

- `log2(W(2^20)) <= 122.759104` (exact rational score `61379552/500000`),
- preprocessing `<= 2^122.402707` group-operation equivalents (included),
- peak memory `<= 2^125.128985` bytes,
- nonuniform advice `0` bytes.

All logarithms below are rounded **up** to six decimal places; the exact
integers in Section 6 are the audit authority. The reference price is
`C = 2000` (batched incomplete affine G1 addition, 6M convention).

**A correction to a naive reading, stated up front.** The textbook
statement of Cheon's attack is `O(sqrt(r/d) + sqrt(d))` *group operations*,
which at `d = 2^20` is about `2^117.5`. That figure counts only the
baby-step/giant-step group operations. Under this challenge's cost model
*all* work must be charged once, including the scalar multiplications that
produce the baby/giant steps, the field arithmetic that produces their
scalars, the sorting and lookup, and the table construction. When all of
that is charged honestly (Sections 5-6), the bound is `2^122.759104`,
not `2^117.5`. This submission claims the honest `2^122.759104`. It still
beats the generic Pollard rho bound of about `2^127.754` (Section 8) by
about five bits, and it sets the first score on this unscored track.

## 2. Number-theoretic facts (all verified by direct computation)

```
r = 52435875175126190479447740508185965837690552500527637822603658699938581184513
  = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
r - 1 = 2^32 * 3 * 11 * 19 * 10177 * 125527 * 859267 * 906349^2
        * 2508409 * 2529403 * 52437899 * 254760293^2
```

The product above was verified to equal `r-1` exactly. In particular
`2^20 | (r-1)`, as the rules state. Set

```
d = 2^20 = 1048576
m = (r-1)/d = 50006747412801924209068050869165387952509453297164571592906626415194112
```

`m` has bit-length 235. The element `zeta = 7` is a generator of the
cyclic group `F_r^*`: for each of the 12 distinct prime factors `l` of
`r-1`, `7^((r-1)/l) != 1 (mod r)` (checked directly; `7` is the smallest
positive integer with this property). Define

```
eta   = zeta^d mod r = 0x628fc184374a1219e85d9e943d57d0c20bad4d8f813859f6121b51c927eaa4dc
theta = zeta^m mod r = 0x3e1c54bcb947035a57a6e07cb98de4a2f69e02d265e09d9fece7e0e39898d4b
```

Then `eta` has exact order `m` in `F_r^*` (it is `zeta^d` with `zeta` a
generator and `d*m = r-1`), and `theta` has exact order `d`: the order
of `zeta^m` is `(r-1)/gcd(r-1, m) = (r-1)/m = d`.

Useful sizes (exact):

```
B1 = ceil(sqrt(m)) = 223621884914696853578958571250384896   (log2 ~ 117.428545)
B2 = ceil(sqrt(d)) = 1024                                  (exact, since d is a square)
b  = ceil(r^(1/4)) = 15132376222941642752                   (log2 ~ 63.714272, b^4 >= r)
```

## 3. The algorithm

Notation: `[a]P` is scalar multiplication of a G1 point `P` by a field
element `a` (taken mod `r`). Group operations are G1 additions/doublings.
All field arithmetic is in `F_r`.

### 3.1 The two phases

Write `s = zeta^k` for the (unknown) discrete logarithm `k` of `s` base
`zeta` in `F_r^*`, `0 <= k < r-1`. Write `k = k0 + k2*m` with
`0 <= k0 < m` and `0 <= k2 < d` (division of `k` by `m`). Then

- `s^d = zeta^(k*d) = zeta^(k0*d) = eta^k0`, because `zeta^(k2*m*d) =
  zeta^(k2*(r-1)) = 1`. So `k0` is determined by `s^d`, and `s^d` is
  exactly the exponent of the supplied point `A_d = [s^d]P0`.
- `zeta^k0` is one `d`-th root of `s^d` in `F_r`, and the `d` roots are
  `zeta^k0 * theta^j` for `j = 0..d-1`, since `theta` has order `d`.
  Hence `s = zeta^k0 * theta^k2` for a unique `k2` in `[0, d)`.

**Phase 1 — find `k0`.** Baby-step/giant-step over `k0` in `[0, m)`:

- Baby table: `B_j = [eta^j]P0` for `j = 0..B1-1`, where
  `B1 = ceil(sqrt(m))`.
- Giant steps: for `i = 0..B1-1`, compute `C_i = [eta^(-i*B1)]A_d`.
  Note `ceil(m/B1) <= B1` because `B1^2 >= m`, so `i` ranges far enough.
- Sort the baby table by `(x, y)` and binary-search each `C_i` in it.
  A match `B_j = C_i` gives `eta^j = s^d * eta^(-i*B1)`, i.e.
  `s^d = eta^(i*B1 + j)`, so `k0 = i*B1 + j`. Exactly one match exists:
  baby points are distinct (distinct scalars `eta^j` mod `r`, `P0` of
  order `r`), and `k0 < m <= B1^2` is covered by the `i, j` ranges.

Then `t := s^d = eta^k0 mod r` is computed by fast powering, and
`S := [zeta^k0]P0` is computed by one scalar multiplication.

**Phase 2 — find `k2`.** Baby-step/giant-step over `k2` in `[0, d)`:

- Baby table: `D_j = [theta^j]S` for `j = 0..B2-1`, `B2 = 1024`.
- Giant steps: for `i = 0..B2-1`, compute `E_i = [theta^(-i*B2)]A_1`.
- A match `D_j = E_i` gives `theta^j * zeta^k0 = s * theta^(-i*B2)`,
  i.e. `s = zeta^k0 * theta^(i*B2 + j)`, so `k2 = i*B2 + j`.

Return `alpha = zeta^k0 * theta^k2 mod r`.

**Verification.** Compute `[alpha]P0` by scalar multiplication and
accept only if it equals `A_1`; also check `alpha^d mod r = t`. If
either check fails, the attack reports failure (it never happens for a
correct execution; see Section 4).

### 3.2 Evaluating the scalar multiplications: b-adic tables

Every baby/giant step above is a scalar multiplication `[delta]G_t`
where `delta` is a full 255-bit field element and `G_t` is a fixed base
point for that sequence. Charging each such step as a fresh
double-and-add scalar multiplication (about 512 group operations) would
make Phase 1 cost about `2^126.4` and would lose to generic rho. The
attack instead uses the standard table method of Kutsuma-Matsuo (SCIS
2007, Section 4, improving Cheon 2006), which evaluates many scalar
multiplications with the *same* base point at `c-1` group additions each
after a one-time table build. This is the same idea as a fixed-base
comb / windowed scalar multiplication, with the window size chosen so
the table fits the accounting below.

Fix `c = 4` and `b = ceil(r^(1/4))` (Section 2), so `b^4 >= r`. Every
`delta` in `[0, r)` has a unique base-`b` expansion
`delta = d0 + d1*b + d2*b^2 + d3*b^3` with `0 <= di < b`. For a base point
`G_t`, build the table

```
T(i, j) = [j * b^i]G_t     for i = 0..3, j = 0..b-1.
```

Construction: `T(0, 1) = G_t`; `T(i, 1) = [b]T(i-1, 1)` for `i = 1..3`
(a scalar multiplication by the 64-bit integer `b`, at most 128 group
operations each); then `T(i, j) = T(i, j-1) + T(i, 1)` for `j = 2..b-1`
(one group addition each). Then for any `delta`,

```
[delta]G_t = T(0, d0) + T(1, d1) + T(2, d2) + T(3, d3),
```

which is 3 group additions (terms with `di = 0` are the point at
infinity and contribute nothing). The identity is exact integer
arithmetic on the base-`b` expansion, not a heuristic.

Four sequences need tables, with four distinct base points, so four
tables are built in total: base `P0` (Phase 1 baby steps), base `A_d`
(Phase 1 giant steps), base `S` (Phase 2 baby steps), base `A_1`
(Phase 2 giant steps). The same table serves every step of its sequence.
The scalar `delta` for each step is produced by one `F_r`
multiplication from the previous step's scalar (e.g.
`eta^(j+1) = eta^j * eta mod r`), plus base-`b` digit extraction;
both are charged in Section 5.

Setup also includes: finding `zeta` (testing `a = 2, 3, ...` against the
12 prime factors of `r-1`, each test a fast powering), computing `eta`,
`theta`, and the four tables. All of it is charged in Section 6 under a
single generous lump sum plus the exact table-build cost.

## 4. Correctness and success probability

The algorithm is deterministic. Every step is exact arithmetic:

- `zeta` is a generator of `F_r^*` (verified in Section 2), so the
  decompositions `s = zeta^k`, `k = k0 + k2*m` exist and are unique in
  their ranges.
- Phase 1: `k0 < m <= B1^2`, so `k0 = i*B1 + j` for some
  `i, j < B1`; the corresponding baby and giant points are equal, so a
  match is found; it is unique because baby points are distinct and the
  giant points are distinct (distinct scalars, base points of order `r`).
  The recovered `k0` satisfies `eta^k0 = s^d` by construction.
- Phase 2: `k2 < d = B2^2`, so `k2 = i*B2 + j` for some `i, j < B2`;
  the same distinctness argument gives a unique match, and the recovered
  `k2` satisfies `s = zeta^k0 * theta^k2`.
- The final verification `[alpha]P0 = A_1` then holds by construction.

Hence the attack recovers `s` for **every** uniform `s` in `[1, r)`,
including the edge cases `s = 1` (`k0 = k2 = 0`) and `s` with
`s^d = 1` (`k0 = 0`, Phase 2 alone determines `s`). Success probability
is **1**, hence `>= 0.39`. There is no randomness, no rejection sampling,
no restart and no amplification to charge. The only heuristic-like
ingredient is none: the b-adic table identity and the BSGS coverage are
theorems, and the toy experiment (`experiment.py`) checks the two-phase
logic and the table identity on small groups (see Section 10).

## 5. Per-unit cost charges (instructions, C = 2000)

Work is counted in scalar AArch64 instruction equivalents, then divided
by the fixed `C = 2000` to give group-operation equivalents, exactly as
the cost model prescribes. One G1 addition/doubling is one reference
operation, i.e. 2000 instructions. All other work is charged at
justified instruction counts; nothing is charged twice (a reference
operation is never also charged through its internal instructions).

| Work item | Charge | Justification |
| --- | --- | ---: |
| G1 addition / doubling | 2000 instructions (1 unit) | Reference operation, fixed `C`. |
| Scalar multiplication `[a]P`, `a` 255-bit | at most 512 group ops = 1,024,000 instructions | Double-and-add: at most 2 additions per bit (the same `2b` bound used in the secp256k1 baseline proof for a scalar multiplication). Used only for `S = [zeta^k0]P0`, the final verification `[alpha]P0`, and the three `[b]T(i-1,1)` table steps (there `b` is 64-bit, so at most 128 group ops each). |
| `F_r` multiplication | at most 1000 instructions (0.5 units) | The organizer's pinned measurement is 383 instructions per `blst_fp_mul` (381-bit field). `F_r` is 255-bit; a Montgomery multiplication plus reduction costs no more, and 1000 is a safe upper bound covering the multiply, the reduction, loads/stores and loop control. See the sensitivity analysis in Section 9 for the much larger secp-style bound. |
| `F_r` addition / subtraction | at most 256 instructions | The secp256k1 baseline's `A` bound for a modular addition on 256-bit limbs. |
| One b-adic table evaluation `[delta]G_t` | 3 group additions + 4 `F_r` multiplications = 3*2000 + 4*1000 = 10,000 instructions | 3 additions from Section 3.2; the 4 multiplications cover updating the running scalar (`delta <- delta * eta` or the analogous factor) and the base-`b` digit extraction (3 divisions by the 64-bit `b`, each bounded by one multiplication's worth of work, plus the update itself). This is charged for *every* baby and giant step. |
| Point comparison (sorting / binary search) | 64 instructions | Comparing `(x, y)`: 8 limb compares plus loads and a branch; `y` is compared only if `x` is equal. |
| Record move (sorting) | 128 instructions | Copying a record (two field elements + index, 104 bytes) plus addressing. |
| Merge-sort loop overhead | 128 instructions per record per pass | The same allowance as the secp256k1 baseline's sort formula. |
| Table-entry store | 128 instructions | Storing one 96-byte point plus addressing. |
| Setup lump sum | 2^30 instructions | Covers finding `zeta` (at most 6 candidates x 12 powerings x ~255 multiplications), computing `eta`, `theta`, fast powerings for `t`, and all small field work. The true cost is orders of magnitude smaller; the lump is a generous upper bound. |
| Verification lump sum | 2^22 instructions | Covers the final scalar multiplication `[alpha]P0` (at most 1,024,000 instructions), the `alpha^d` field powering, encoding checks and output. The lump is headroom on top of those items, all inside it. |

Sorting is by deterministic bottom-up merge sort with two buffers,
following the secp256k1 baseline's formula: sorting `N` records costs
`N * L * (CMP + 2*MOVE + LOOP)` instructions with `L = ceil(log2 N)`.
A hash table is deliberately **not** used: its total work is only an
expected value (collision-dependent), while the claim must be a
deterministic upper bound, and merge sort gives one. (A hash table or a
radix sort would tighten the bound; they are not claimed here, and the
reviewer is asked to review the submitted bound, not to tighten it.)

Memory transfers for table lookups are included in the per-step charge:
each table evaluation reads 4 table entries (96 bytes each); the loads
are inside the 10,000-instruction step budget, which has large headroom
over the 3 group additions (6,000 instructions) plus 4 multiplications
(4,000 instructions) it nominally covers.

## 6. Total work as a function W(n), and the score at n = 2^20

For general `n` with `n | (r-1)`, `1 <= n < r`, put `d = n`,
`m(n) = (r-1)/n`, `B1(n) = ceil(sqrt(m(n)))`, `L1(n) = ceil(log2 B1(n))`,
`B2(n) = ceil(sqrt(n))`, `L2(n) = ceil(log2 B2(n))`. With
`STEP = 10000` instructions (Section 5), one table build
`TAB = 3*128*2000 + 4*(b-2)*2000 + 4*b*128` instructions (with
`b = ceil(r^(1/4))`, independent of `n`), and `C = 2000`,

```
W(n) = (1/C) * [ 2^30                                    (setup)
               + 4*TAB                                  (four b-adic tables)
               + B1(n)*STEP                             (Phase 1 baby steps)
               + B1(n)*L1(n)*448                        (sort Phase 1 baby table)
               + B1(n)*STEP + B1(n)*(L1(n)+1)*64         (Phase 1 giant steps + binary searches)
               + 512*2000                               (S = [zeta^k0]P0)
               + B2(n)*STEP                             (Phase 2 baby steps)
               + B2(n)*L2(n)*448                        (sort Phase 2 baby table)
               + B2(n)*STEP + B2(n)*(L2(n)+1)*64         (Phase 2 giant steps + binary searches)
               + 2^22 ]                                 (verification)
```

Notes on the formula: the Phase 1 giant-step count is charged as
`B1(n)`, an upper bound on `ceil(m(n)/B1(n))`; the Phase 2 giant-step
count is charged as `B2(n)`, an upper bound on `ceil(n/B2(n))`. The
four tables are charged in full even though the two Phase 2 tables are
built after the Phase 1 tables could be freed. No work is amortized
over future targets: every term above is for this one challenge.
`W(n)` is claimed only for `n | (r-1)`; for other `n` the attack as
stated does not apply (the published score uses `n = 2^20`, which
satisfies the condition, and `A_d` is supplied at that `n`).

Exact integers at `n = 2^20` (`B1 = 223621884914696853578958571250384896`,
`L1 = 118`, `B2 = 1024`, `L2 = 10`, `b = 15132376222941642752`,
`TAB = 128806786409679263857024` instructions):

| Component | Instructions (exact integer) |
| --- | ---: |
| setup | 1073741824 |
| four b-adic tables | 515227145638717055428096 |
| Phase 1 baby steps | 2236218849146968535789585712503848960000 |
| sort Phase 1 baby table | 11821547324130534467598065910580347142144 |
| Phase 1 giant steps + binary searches | 3939323124657299772646934191146780327936 |
| `S = [zeta^k0]P0` | 1024000 |
| Phase 2 baby steps | 10240000 |
| sort Phase 2 baby table | 4587520 |
| Phase 2 giant steps + binary searches | 10960896 |
| verification | 4194304 |
| **Total instructions** | **17997089297934803291261731452949136606720** |
| Preprocessing instructions (setup + tables + both baby builds + both sorts) | 14057766173277503518614797261802340099584 |

Dividing by `C = 2000`:

- `log2(total instructions) = 133.724888` (rounded up),
- **`log2(W(2^20)) = 122.759104`** (rounded up; exact value
  `122.75910310655439699...`), i.e. the submitted rational score
  `61379552/500000`. Direct check: `2^(61379552/500000) >= W(2^20)`.
- `log2(preprocessing) = 122.402707` (rounded up).

The dominant term is sorting the Phase 1 baby table
(`2^122.15` group-operation equivalents of the `2^122.76` total);
the baby and giant steps together contribute about `2^121.2`. The
four b-adic tables together cost about `2^67.8` group-operation
equivalents — negligible at this scale, which is why the table method
is worth using: it converts each 512-group-operation scalar
multiplication into a 5-group-operation-equivalent table evaluation.

## 7. Memory, preprocessing, advice

- **Peak memory** (upper bound): two merge-sort buffers of `B1` records
  (`2 * B1 * 104` bytes, a record being two 48-byte field elements plus
  an 8-byte index) + all four b-adic tables simultaneously
  (`4 * 4 * b * 96` bytes) + a `2^20`-byte reserve for code, constants
  and scratch, totalling `46513352062256968787753261258444374016` bytes,
  i.e. `log2 = 125.128985` (rounded up). The tables could be freed
  between phases, so this is an upper bound, not an exact peak.
- **Preprocessing**: setup + four table builds + both baby-table builds
  + both sorts = `14057766173277503518614797261802340099584`
  instructions, i.e. `2^122.402707` group-operation equivalents,
  included in the total above. No preprocessing is reused across
  targets and no credit is taken for it.
- **Nonuniform advice**: `0` bytes. The b-adic tables are built by the
  attack itself from the public inputs (`P0`, `A_d`, `A_1`, and the
  computed point `S`) and their construction is charged in full above;
  they are not embedded precomputed answers. No other advice is used.

## 8. Comparison with the generic attack

For a group of prime order `r`, Pollard's rho algorithm needs about
`sqrt(pi*r/2)` group operations, i.e. `log2 = 127.754293` here; with the
negation map (which applies to this curve, points and their negatives
being distinct except at infinity) about `sqrt(pi*r/4)`, i.e.
`log2 = 127.254293`. This submission's `122.759104` beats generic rho
by about 5.0 bits (4.5 bits against negation-map rho), using only the
supplied auxiliary powers — exactly the advantage Cheon's attack is
known for, reduced from the naive `2^117.5` figure by the honest
charging of sorting, field arithmetic and table construction
(Section 1). The BSGS baseline control on the secp256k1 track
(`129.000802`) is a different track and a deliberately conservative
control; it is not the comparison point here. There is no prior score
on this track: the selected template is unscored.

## 9. Sensitivity: what if field arithmetic costs far more?

The most contestable charge is the `F_r` multiplication at 1000
instructions. The organizer's own pinned measurement is 383
instructions per 381-bit `blst_fp_mul`; 1000 for a 255-bit multiply is
already generous. But suppose instead the much larger secp-style
bit-serial bound were used: 147,712 instructions per `F_r`
multiplication (73.856 group-operation equivalents). Then each
b-adic table evaluation costs `3*2000 + 4*147712 = 596,848`
instructions, and recomputing the Section 6 total with that single
change gives `log2(W(2^20)) = 126.721069` — still about 0.5 bits better
than generic rho with negation (`127.254293`) and 1.0 bit better than
plain rho (`127.754293`). So the claim that Cheon beats generic rho on this track
does not depend on the 1000-instruction charge; it survives even a
~148x more expensive field multiplication. The submitted bound
`122.759104` does depend on the 1000-instruction charge, which is
justified in Section 5.

Similarly, if merge sort were charged at twice the per-record cost used
here, the total would rise to about `2^123.5`, still beating generic
rho. The margin is not balanced on any single charge.

## 10. Heuristics, limitations, and what the experiment shows

- **Heuristics: none.** The attack is deterministic; every step is
  exact arithmetic, and success probability is 1 by the argument in
  Section 4. No random-oracle, uniformity, smoothness or independence
  assumption is used. (The uniformity of `s` matters only for the
  success-probability statement, which is 1 for every `s`.)
- **Limitations.** The bound is an analytic upper bound, not a measured
  runtime; no full-size execution is claimed or possible. Memory is
  astronomical (`2^125.13` bytes), as with the secp256k1 BSGS control
  (`2^136` bytes); this is not a practical attack. The attack needs
  `A_d` with `d = n = 2^20` supplied, which the track provides at the
  scoring `n`. The `W(n)` formula is claimed only for `n | (r-1)`.
- **Experiment.** `experiment.py` implements the exact two-phase logic
  of Section 3 (with direct scalar multiplications in place of the
  b-adic tables, which are a production-scale optimization) on additive
  groups `Z/qZ` for `(q, d) = (61, 12), (101, 10), (271, 18)`, and
  separately checks the b-adic table identity
  `[delta]G_t = sum_i T(i, d_i)` on the same toy groups. It recovers
  the secret in every case, including the edge cases `s = 1` and
  `s = q-1` (checked locally during preparation). It does not establish
  production cost bounds, memory bounds, or G1 arithmetic; those rest
  on the analytic accounting in Sections 5-6.
- **What could still be wrong.** The main risks are (a) an arithmetic
  slip in the exact integers of Section 6 (they were computed by
  script from the formulas in Section 5; the formulas, not just the
  numbers, are stated so they can be rechecked), (b) a charge the
  reviewer judges too low (the sensitivity analysis in Section 9 bounds
  the consequence of the largest such charge), and (c) a subtle
  mismatch between the track's input encoding and the algorithm's
  inputs (the algorithm uses only `A_1` and `A_d`, both supplied, and
  returns a 32-byte scalar `s` with `[s]P0 = T`, as required).

## References

- J. H. Cheon, "Security analysis of the strong Diffie-Hellman
  problem", EUROCRYPT 2006, LNCS 4004, pp. 1-11. The auxiliary-input
  DLP algorithm and its `O(sqrt(p/d) + sqrt(d))` group-operation
  analysis; see also J. Cryptology 23(3), 2010, pp. 457-476.
- T. Kutsuma and K. Matsuo, "Remarks on Cheon's algorithms for
  pairing-related problems", SCIS 2007. The b-adic table method of
  Section 3.2 and the improved `O(sqrt(p/d) + sqrt(d))` analysis
  without the `log p` factor; Theorem 1' in that paper is the version
  used here.
- A. Menezes, P. van Oorschot and S. Vanstone, Handbook of Applied
  Cryptography, Section 3.6 (generic DLP algorithms: BSGS, Pollard
  rho, Pohlig-Hellman).
- The track's own `target-profile.json`,
  `classical-group-work-v2.json` and `tracks/docs/COST_MODEL.md` for
  the target, the reference price `C = 2000`, and the accounting rules
  applied in Sections 5-6. The `r-1` factorization, the generator
  `zeta = 7`, and the values of `eta` and `theta` in Section 2 were
  verified by direct computation from the profile's `r`.
