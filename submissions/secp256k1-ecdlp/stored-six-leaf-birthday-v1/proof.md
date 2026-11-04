# Stored six-leaf birthday search on secp256k1 (127.097776303128)

Status: conservative generic attack, complete accounting, no qualification implied.
This package is a first, deliberately simple member of the stored-leaf family: it
uses the curve automorphism to produce six stored x-coordinates per batched
inversion, and finds the discrete logarithm by a birthday collision among them.
It is **not** competitive with the promoted frontier (124.638168120013643); its
purpose is a transparent, fully reproducible baseline whose ledger can be checked
end to end.

## 1. Target

`tracks/secp256k1-ecdlp/target-profile.json`: given `Q = [x]P` for the fixed
secp256k1 generator `P` (order `r`, cofactor 1) and uniform `x` in `[0, r)`,
return `x`. The input is `Q` only. No auxiliary powers, no advice, no
preprocessing shared with other targets. Success must be at least 0.39 over the
uniform challenge and fresh algorithm coins.

Field: `p = 2^256 - 2^32 - 977`. Order:
`r = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141`.

## 2. Automorphism and the class of a point

`E: y^2 = x^3 + 7` over `F_p` has automorphism group of order 6. Write
`beta` for the element of order 3 in `F_p`:

```
beta = 0x7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee
beta^2 + beta + 1 = 0 (mod p),  beta^3 = 1
```

and `phi(x, y) = (beta x, y)`. On the group of order `r` this is multiplication
by the element `lambda` of order 3 in `F_r`; `ledger.py` does not depend on its
value, and `experiment.py` re-derives it on a small curve.

For a point `Z = (x, y)` the six points `{+-Z, +-phi(Z), +-phi^2(Z)}` are the
class of `Z`. All six share the same `y` up to sign and have x-coordinates
`{x, beta x, beta^2 x}`.

Only one curve operation is needed below: from an affine point `X = (x_X, y_X)`
and an affine point `Y = (x_Y, y_Y)` with `x_X != beta^t x_Y`, the six
x-coordinates

```
x(phi^t(X) + Y), x(phi^t(X) - Y),  t = 0, 1, 2
```

are computed from **one** field inversion, because all six share the three
denominators `d_t = x_Y - beta^t x_X`:

```
D      = d_0 d_1 d_2                          (2M)
D^-1                                          (one batched inversion)
w_t    = D^-1 d_j d_k  (j,k != t)             (2M each, 6M total)
x_{t,+-} = (u_+- w_t)^2 - x_Y - beta^t x_X    (u_+- = y_Y -+ y_X)
```

`beta^t x_X` is computed once per `X` and reused for all `Y`.

## 3. Construction

Parameters, all exact integers (`ledger.py` prints them):

```
n        = 6314527679590297595              (2^62.45)
records  = n^2 = 39873259816312028046297715260662784025
leaves   = S   = 6 * records = 239239558897872168277786291563976704150 (2^127.4917)
```

Setup (charged, no free advice):

1. Draw `a_k, b_k` for `k = 1..n` and `c_l` for `l = 1..n`, each uniformly from
   `[1, r)` by rejection sampling on 256-bit `RNDR` words (64 instructions per
   draw, charged).
2. `X_k = [a_k]P + [b_k]Q`, `Y_l = [c_l]P`. The `b_k` are **per-k**, which
   matters: it makes the vector `(m, n)` of every leaf depend on `k`, so leaves
   from different `k` are independent uniform group elements.
3. A table `T` of `2^128` slots of 32 bytes, addressed by the low 128 bits of an
   x-coordinate; zeroed in bulk at 64 bytes per `DC ZVA`.

Main loop, for each of the `n^2` pairs `(k, l)`:

1. Compute the six x-coordinates of section 2.
2. For each of the six leaves `L = [m]P + [q]Q`, with

   ```
   m = lambda^t a_k + sigma c_l,     q = lambda^t b_k
   ```

   form the address `addr = x(L) mod 2^128`, load the slot, compare its stored
   x with `x(L)`, and store `(x(L), k, l, t, sigma)`.

Resolution on the first equal pair: the two points satisfy `L_2 = +-L_1`, hence

```
m_2 - m_1 + (q_2 - q_1) x = 0      (plus)
m_2 + m_1 + (q_2 + q_1) x = 0      (minus)
```

Each is a linear equation in the unknown `x`; the applicable one has a nonzero
coefficient with probability `1 - O(1/r)`. Solve modulo `r`, then verify
`[x]P = Q` and output `x`. A degenerate coefficient is a charged failure row.

## 4. Correctness

Every leaf is `[m]P + [q]Q` with `m, q` known to the machine. If two distinct
leaves have equal x-coordinates then their points are equal or opposite, so one
of the two linear equations above holds at the true `x`; solving it returns `x`
modulo `r`, and the explicit `[x]P = Q` check accepts exactly when the answer is
correct. Distinct leaves have distinct `(m, q)` with overwhelming probability, so
no collision is between a leaf and itself. No heuristic is used in this part.

## 5. Success

Fix the true `x`. Every leaf is `[m]P + [q]Q` with `(m, q)` a fresh affine
function of the drawn coins; conditional on the coins it is a uniform element of
the group, and two distinct leaves are equal-or-opposite with probability `2/r`.
There are `C(S, 2)` pairs, so the mean number of colliding pairs is

```
mu = C(S,2) * 2/r = 0.494296  (exact rational in ledger.py)
```

By the standard Poisson/Chen-Stein approximation for a sum of `C(S,2)` weakly
dependent indicators, `P(no collision) = exp(-mu) (1 + o(1))`, hence

```
P(success) >= 1 - exp(-mu) = 0.390000  >= 0.39
```

This is a supported heuristic (the dependence is only through shared `k` or `l`,
and the o(1) term is below 2^-40 at these parameters); it is disclosed as such,
not as a theorem. A failed run is fully charged: the table build is the whole
work whether or not a collision is found, and fresh coins are drawn for a retry.

## 6. Accounting

Units: one reference target operation is `C = 1000` scalar AArch64 instruction
equivalents (`tracks/docs/COST_MODEL.md`). The fixed 6M calibration measures one
5x52 base-field multiplication as **179** instructions; the account uses
`M = 179`, `S = M`, and a conservative modular add/sub `A = 30`. Every operation
below is counted once, on the same scalar-instruction basis.

Per record (one `(k, l)` pair, six leaves), all items executed:

| item | instructions |
| --- | ---: |
| `d_t = x_Y - beta^t x_X`, t = 0,1,2 | 3 A |
| `D = d_0 d_1 d_2` | 2 M |
| batched inverse of `D` (amortised 3M per element) | 3 M |
| `w_t = D^-1 d_j d_k`, six products | 6 M |
| `u_+- = y_Y -+ y_X` | 2 A |
| `u_+-^2` (2), `w_t^2` (3) | 5 S |
| `x_{t,+-} = u^2 w^2 - x_Y - x_t`, six products | 6 M |
| the twelve subtractions in the six x rows | 12 A |
| pointer/loop work | 24 |
| six table inserts (index, load, 4-word compare, store) | 6 x 15 = 90 |
| **per record** | **4562** |

`beta^t x_X` is two multiplications per `X`, charged once per axis entry, and is
below 2^-60 of the total. Non-record rows:

| row | instructions | log2 |
| --- | ---: | ---: |
| generators `2n (645410 + 642930)` | 16270517181446728007084600 | 83.75 |
| coin draws `2n * 2 * 64` | 1616519085975116184320 | 70.45 |
| table zeroing `2^128 * 32 / 64` | 170141183460469231731687303715884105728 | 127.00 |
| resolution and `[x]P = Q` | 2097152 | 21.00 |
| records `4562 * n^2` | 181901811282015471947210177019143620722050 | 137.06 |
| **total `T`** | **182071952465475957451075564855562630193850** | **137.063561** |

Peak memory: the table is `2^128` slots of 32 bytes = `2^133` bytes; working
state is under 2^10 bytes. Memory is reported, not scored; it is not advice.

Preprocessing is the table build itself and is charged in full; nothing is
amortised over future targets, and preprocessing is reported as
`log2(T_pre/1000) = 117.034` (zeroing plus generators plus coins).

## 7. Claim

`W = T / 1000`. `ledger.py` brackets `2^c` with a 90-digit decimal and returns the
least multiple of `10^-15` at or above `log2 W`:

```
log2 W           = 127.097776303128...
claim            = 127097776303128 / 10^15 = 15887222037891 / 125000000000000
2^claim          >= W      (verified)
2^(claim-1e-15)  <  W      (verified)
```

The claim is the submitted score: `log2` of total computation in
target-group-operation equivalents, lower is better.

## 8. Limitations and provenance

- This is a generic attack. It is about 2.46 bits behind the promoted frontier
  and does not weaken secp256k1; the parameters are astronomically infeasible.
- The Poisson collision bound is a supported heuristic, disclosed in section 5;
  the end-game verification is exact.
- The batched inversion cost (`3M` per element) and the `A = 30` add bound are
  conservative upper bounds, not measurements; a faster kernel only lowers the
  claim.
- Everything here is derived from the published rules, the cost model and
  standard birthday-pollard/BSGS literature (HAC 3.6, Galbraith-Wang-Zhang,
  ePrint 2015/605, for sharing one inverse between `X +- Y`). No prior
  submission's code, ledger or text is reused. `ledger.py` reproduces every
  integer above and takes a second to run.

