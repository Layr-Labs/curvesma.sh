# Stored six-leaf birthday search on secp256k1 (v2, 127.622149460486554)

Status: conservative generic attack, complete accounting, no qualification implied.
This is v2 of a rejected v1 package. The four defects named in that rejection are
fixed here and are listed first, so a reviewer can check them directly.

## 0. What changed from v1

| v1 defect | v2 fix |
| --- | --- |
| Claim `15887222037891/125000000000000` equals 0.127, not 127.1 | The claim is now emitted by `ledger.py` as the exact reduced rational `63811074730243279/500000000000000 = 127.622149460486554`, bracketed on both sides |
| Collision pairs occur in dependent pairs, so `mu` overcounts | The success bound now uses `mu/2`: `P(no collision) <= exp(-mu/2)`, and the parameters are chosen so `mu/2 = 0.5`, giving `0.3935 >= 0.39` |
| The 2^128-slot table overwrites slots, so not every leaf is retained | The table is a 2^129-head bucket array with **chained records and no eviction**; every inserted leaf stays reachable |
| 32-byte slots cannot hold the record, and axis/buffer storage was uncharged | Records are 64 bytes (`x`, `k`, `l`, `tag`, `next`); axis scalars, axis points, batch buffers and head zeroing are itemised and charged |

## 1. Target

`tracks/secp256k1-ecdlp/target-profile.json`: given `Q = [x]P` for the fixed
secp256k1 generator `P` (order `r`, cofactor 1) and uniform `x` in `[0, r)`,
return `x`. Input is `Q` only; no auxiliary powers, no advice, no preprocessing
shared with other targets. Success must be at least 0.39 over the uniform
challenge and fresh algorithm coins.

```
p = 2^256 - 2^32 - 977
r = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
```

## 2. The six x-coordinates of one batch

`E: y^2 = x^3 + 7` over `F_p` has an automorphism group of order 6. Let `beta`
be the element of order 3 in `F_p`:

```
beta = 0x7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee
beta^2 + beta + 1 = 0 (mod p), beta^3 = 1, phi(x, y) = (beta x, y)
```

`phi` is multiplication by an element `lambda` of order 3 in `F_r`.

For affine `X = (x_X, y_X)` and `Y = (x_Y, y_Y)` with `x_X != beta^t x_Y`, the
six points `phi^t(X) +- Y` have x-coordinates

```
x = (u_+- / d_t)^2 - x_Y - beta^t x_X,   d_t = x_Y - beta^t x_X,
u_+- = y_Y -+ y_X
```

because `phi^t(X)` shares `y_X` and has x-coordinate `beta^t x_X`. All six share
the three denominators `d_0, d_1, d_2`, so one inversion of
`D = d_0 d_1 d_2` gives all six:

```
w_t = D^-1 d_j d_k      (j, k != t)
x_{t,+-} = (u_+- w_t)^2 - x_Y - beta^t x_X
```

`beta^t x_X` is computed once per `X` and reused. This is the standard
Galbraith–Wang–Zhang sharing of one inversion between `X + Y` and `X - Y`,
extended to the three twists of `X`.

## 3. Construction

All parameters are exact integers printed by `ledger.py`:

```
n         = 7530851732716320753                (2^62.7075)
records   = n^2 = 56713727820156410591308303813178487009
leaves S  = 6 n^2 = 340282366920938463547849822879070922054  (2^128.0000)
```

1. **Axes.** Draw `a_k, b_k` for `k = 1..n` and `c_l` for `l = 1..n`, each
   uniform in `[1, r)` by rejection on 256-bit `RNDR` words (64 instructions per
   draw, charged). Set `X_k = [a_k]P + [b_k]Q` and `Y_l = [c_l]P`. The axis
   points and all `3n` scalars are stored (section 6); they are not free advice.
2. **Leaves.** For every pair `(k, l)` and every `t` in `{0,1,2}`, `s` in
   `{+1,-1}`, the leaf

   ```
   L = phi^t(X_k) + s Y_l = [m]P + [q]Q,
   m = lambda^t a_k + s c_l,   q = lambda^t b_k
   ```

   is inserted. Six leaves per record, `n^2` records.
3. **Table.** `H` is an array of `2^129` head pointers and `R` an array of `S`
   64-byte records. A record holds `x` (32 B), `k` (8 B), `l` (8 B), a tag word
   with `t` and `s` (8 B), and `next` (8 B). `H` is zeroed; `R` needs no
   initialisation because every record is written before it is linked.
4. **Insert.** For a leaf with coordinates `(x, k, l, t, s)`: `h = x mod 2^129`;
   walk the chain `H[h]` comparing `x`; on a hit resolve (section 4); otherwise
   write the record, point it at the old head, and store it into `H[h]`.
   Nothing is ever overwritten or evicted, so every leaf remains reachable.

The expected chain length is `S / 2^129 = 2^-1`, so an insert touches one stored
record on average.

## 4. Resolution

If two leaves have equal x-coordinates then their points are equal or opposite,
so exactly one of

```
m_2 - m_1 + (q_2 - q_1) x = 0        (L_2 = +L_1)
m_2 + m_1 + (q_2 + q_1) x = 0        (L_2 = -L_1)
```

holds at the true `x`. Both coefficients are known from the stored tags and the
axis scalars; solve modulo `r`, then check `[x]P = Q` and output `x`. A zero
coefficient is a charged failure row (probability `O(1/r)`); the corresponding
line is then skipped and the other solved. Correctness is exact: the check
accepts precisely when the returned scalar is the discrete logarithm.

## 5. Success

Two distinct leaves are equal or opposite with probability `2/r`, and there are
`C(S,2)` pairs, so

```
mu = C(S,2) * 2/r = S(S-1)/r = 1.000000000   (exact in ledger.py)
```

`mu` alone would give `1 - exp(-mu) = 0.632`. It cannot be used, because the
construction creates **dependent collision pairs**: for `k != j`, `l != m`, if
`phi^t(X_k) + sigma Y_l = eps (phi^u(X_j) + tau Y_m)` with `eps = +-1`, then
applying `eps` to both sides and rearranging gives

```
phi^t(X_k) - eps tau Y_m = eps (phi^u(X_j) - eps sigma Y_l),
```

and both sides are other declared leaves. Every cross-index collision is
therefore accompanied by a partner, so the number of independent collision
orbits is `mu/2`, not `mu`; shared-index collisions (`k = j` or `l = m`) are
`O(n^3/r) = 2^-68` and need no partner. With the orbit-Poisson bound

```
P(no collision) <= exp(-mu/2) = exp(-0.5) = 0.606531,
P(success) >= 1 - exp(-0.5) = 0.393469 >= 0.39.
```

The orbit-Poisson bound is the only heuristic and is disclosed as such. The
construction has no other heuristic: the six-leaf formulas are exact algebra, the
resolution is exact linear algebra, and the final test `[x]P = Q` is exact.

## 6. Accounting

Units: one reference target operation is `C = 1000` scalar AArch64 instruction
equivalents. The fixed calibration measures a 5x52 base-field multiplication as
`M = 179` instructions; squaring is charged as `M`; a modular add/sub is charged
at the conservative `A = 30`.

Per record (one `(k, l)` pair, six leaves), all items executed:

| item | instructions |
| --- | ---: |
| `d_t = x_Y - beta^t x_X`, t = 0,1,2 | 3 A |
| `D = d_0 d_1 d_2` | 2 M |
| batched inverse of `D`, amortised 3M per element | 3 M |
| `w_t = D^-1 d_j d_k`, six products | 6 M |
| `u_+- = y_Y -+ y_X` | 2 A |
| `u_+-^2` (2) and `w_t^2` (3) | 5 M |
| `x = u^2 w^2 - x_Y - x_t`, six products | 6 M |
| the twelve subtractions | 12 A |
| pointer and loop work | 24 |
| six chained inserts (index, head load, one compare, record store, link) | 6 x 24 = 144 |
| **per record** | **4616** |

Non-record rows:

| row | instructions | log2 |
| --- | ---: | ---: |
| records `4616 n^2` | 261790567617841991289479130401631896033544 | 137.587 |
| generators `2n (645410 + 642930)` | 19404595042655489357840040 | 84.005 |
| coin draws `3n * 64` | 1445923532681533584576 | 70.292 |
| axis stores (points and scalars) | 52715962129014245271 | 65.515 |
| head zeroing `2^129 * 8 / 64` | 85070591730234615865843651857942052864 | 126.000 |
| resolution and `[x]P = Q` | 2097152 | 21.000 |
| **total T** | **261875638209572245311438656203789745853447** | **137.588** |

Nothing is amortised over other targets; preprocessing (heads, axes, generators)
is charged in full. `T` counts every executed instruction of the construction
plus the itemised allowances.

## 7. Storage

| object | bytes |
| --- | ---: |
| 2^129 head pointers (8 B) | 2^132 |
| S records (64 B) | 2^134.6 |
| 3n scalars (32 B) + 2n points (64 B) | 2^68.3 |
| batch-inversion buffer (one limb array per lane, 2^10 lanes) | 2^15 |
| **peak** | **2^134.32** |

Memory is reported, not scored, and is construction storage, not advice. Every
stored object is charged when written; no object is assumed free.

## 8. Claim

`W = T / 1000`. `ledger.py` brackets `2^c` with 90-digit decimal arithmetic:

```
log2 W          = 127.622149460486554...
claim           = 63811074730243279 / 500000000000000 = 127.622149460486554
2^claim         >= W      (verified)
2^(claim-1e-15) <  W      (verified)
```

The submitted score is this reduced rational; `ledger.py` re-derives it from the
integer ledger, so the spelling and the value cannot drift apart.

## 9. Evidence

* `ledger.py` — prints every integer above, the exact `mu`, the success bound,
  and both bracket checks; runs in under a second.
* `experiment.py` — on `y^2 = x^3 + 7` over `F_127` (prime order 127) it checks
  `beta` and `lambda`, draws the Y axis once and reuses it for every X, inserts
  leaves into a chained bucket table with no overwrite, resolves a collision by
  the two linear equations and verifies `[x]P = Q`. Toy scale only.

## 10. Limitations and provenance

* Generic attack: about 3.0 bits behind the promoted frontier and no threat to
  secp256k1. The parameters are astronomically infeasible.
* The orbit-Poisson success bound is a supported heuristic, disclosed above.
* `A = 30` and the 3M amortised inversion are conservative upper bounds, not
  measurements.
* Derived from the published rules, the cost model, and standard literature
  (HAC 3.6; Galbraith–Wang–Zhang, ePrint 2015/605). No prior submission's code,
  ledger or text is reused.
