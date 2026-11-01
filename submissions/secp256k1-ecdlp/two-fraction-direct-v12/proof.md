# Secp256k1 two-fraction labels with two-phase direct bucketing

Proposed complete bound125.142304476874=125142304476874/10^12 target-group-operation equivalents at C1000, total46944787124746853008406832072376781373440 scalar instructions. Success at least0.390000001, total peak memory below2^136 bytes, preprocessing below2^111 units included in total, no advice or auxiliary powers. The exact numerical certificate and complete source argument are operative together. DIRECT_PROTOCOL.md replaces the original two-slot transport with a fully charged receiver construction. No prior score or acceptance is a mathematical premise. This complete analytical candidate requires external exact-version review and user authorization before submission; finite checks are not independent production certification.

The two-fraction label, rare arithmetic fixups, all label-error and node-local meter losses, six leaves per incidence, field kernels and paired scanner are inherited unchanged from public bbe1cd85-58b0-4100-88bc-a095f03ae1e6, commit287f9790b7285be84c0bf913c055d3b069d7839e. New work is compatible routing/layout/cost integration, not addition of our alternative orbit-min key to that key. q1=145289002, q2=2^25,N=2^72, CAP1=7*2^19,CAP3=14755795394. Both meters remain persistent over all q1 batches across both phases.

## 1. Overview

Each node computes the **leaves** `x(φ^t(X) ± Y)`, t = 0, 1, 2, for independent points X and Y, six per pair from one
batched inverse: Galbraith, Wang and Zhang (ePrint 2015/605) share one inverse between `X ± Y`, and φ(x, y) = (βx, y)
extends this to the three twists of X. Labelled by two fractions of its class's canonical x, each leaf goes in two
bulk-synchronous hops to a home node, where equal keys give the logarithm. This is an analytic bound, not a practical
break.

## 2. The algorithm

### 2.1 Machine

Use the explicit scalar machine, ordered blocking ports, separate resident receiver context, charged release acknowledgement, in-flight accounting and barrier construction in DIRECT_PROTOCOL.md. Each node has private memory below2^64 bytes; all instructions of every context/processor are summed. FEAT_RNG implements the stated ideal coin source and FEAT_LSE the scalar atomics. No polling, callback, interrupt or message dispatch is assumed free. No physical hardware throughput claim follows from this instruction model.

### 2.2 Leaves and records

Node g draws fresh uniform `b_g` and `a_1..a_q1` in `[0, r)` and `c_1..c_q2` in `[1, r)`, with `q1 = 145289002` and `q2
= 2^25`, and stores them with `R_g = [b_g]Q`, `X_k = R_g + [a_k]P` and `Y_l = [c_l]P`. Its leaves are the x-coordinates
of `φ^t(X_k) + s*Y_l` for all k, l, t in {0, 1, 2} and s = ±1: `n0 = 6*q1*q2` per node and `S = 2^72*n0 = 2^126.699317`
in all. As φ = [λ], leaf `φ^t(X_k) + s*Y_l` is `[m]P + [n]Q` with `m = λ^t*a_k + s*c_l` and `n = λ^t*b_g`. A point Z ≠ O
has the class `{±Z, ±φ(Z), ±φ²(Z)}`, its x-coordinates `{x, βx, β²x}` (mod p) have the canonical member `v* = min(x, βx
mod p, β²x mod p)`, and its **label** is `(comp(v*), comp(βv* mod p))`, `comp(v) = floor(2^100 v/p)`: two 100-bit
components shared by the class (3.3). A record is `[K2, K3, w2, w3]` with K2 and K3 the top 64 bits of the first and the
second component, `w2 = i | k<<36` and `w3 = d2 | c<<36 | l<<39`, where `c = 1 + 2t + [s = -1]` in 1..6 is the leaf
type; its digits are d1 and d2, the components' low 36 bits (hop-1 row, hop-2 column), `d3 = K2 mod 2^24` and `d4 = K3
mod 2^48`. At its home a record's source is `(w2 mod 2^36, j)`, where j=w3 mod2^36 after hop2 replaces the destination field, as specified in2.3.

### 2.3 Phases

Preprocess the complete Y table, headers, constants and scalar labels as before. Initialize d4 counters and tight full-record buckets once. Partition the original X-index sequence into two consecutive phases, each of at most ceil(q1/2)=72644501 batches, without renumbering k. Arithmetic and j-path counters and admission guards span both phases; never reset either meter at a phase boundary.

For each phase initialize hop-1 cursor tables to the sender-row slot bases. The unchanged key emit advances its destination cursor by32 even on a full slot. At phase end send clamped final cursors with fin; drain and barrier. h2s streams all retained first-hop records directly to their home row/column, extracts d2 before replacing low36 bits of w3 with original source column j, and sends the four-word payload. The receiver reserves a full32-byte bucket entry, checks its explicit end, stores/counts or rolls back, then releases the payload window by STLR. Drain and barrier before reusing hop-1 storage. There are no hop-2 pair slots or finishing cursors.

Buckets and d4 counts persist across both phases. Afterwards prefix4 computes d4 starts and pass4 scatters pointers to full records in d3 order. scan_ranges lists K2-adjacent pairs within d4 ranges, capped at2^28. Source row and k remain in w2; source j now comes from low36 bits of w3, whose type bits36..38 and l bits39..63 are unchanged. Candidate processing includes64 additional instructions for decoding this layout.

Resolution retains four count-exchange subphases, source cap2^29 and relay cap2^34. Requests and replies use the stored points/scalars, recompute the two leaves and six-unit relation, solve the nonzero denominator and verify[x]P=Q. See the aggregate conservation accounting in DIRECT_PROTOCOL.md: a relay's local cap is not confused with the average total exchange bound.64barriers and2^14phase-control instructions per node cover all transitions and receiver setup; fixed kernel returns are charged separately.

## 3. Leaf kernels (`kernels/leaf.s`)

### 3.1 One batch

A batch is one X and all q2 Y's. With `D_l = y(Y_l)^2 - y(X)^2`, the forward pass stores `D_l` and `acc_(l-1)` in entry
l and forms `acc_l = acc_(l-1)*D_l` (`acc_(-1) = 1`); fe_inv gives `R = 1/acc_(q2-1)`; the backward pass, for
`l = q2-1..0`, sets `inv = R*acc_(l-1)`, `R <- R*D_l`, `u = inv*(y_Y - y_X)` and `v = inv*(y_Y + y_X)`, and for t = 0,
1, 2 writes

    x(φ^t(X) + Y) = Q_t²u² - β^t x_X - x_Y,    x(φ^t(X) - Y) = Q_t²v² - β^t x_X - x_Y,
    Q_t = x_Y² + β^t x_X x_Y + β^(2t) x_X².

**Why.** As R = 1/acc_l before step l, inv = 1/D_l. The curve gives `y² = x³ + 7`, so `D = x_Y³ - x_X³`, and as β³ = 1,
`(x_Y - β^t x_X)*Q_t = x_Y³ - x_X³ = D`. The slope of `φ^t(X) ± Y` with `φ^t(X) = (β^t x_X, y_X)` is
`(±y_Y - y_X)/(x_Y - β^t x_X) = (±y_Y - y_X)*Q_t/D`, which is `Q_t*u` or `-Q_t*v`; its square, less the two
x-coordinates, is the leaf. D ≠ 0 unless `x_Y³ = x_X³`, that is `Y = ±φ^t(X)` (3.4).

**Arithmetic per pair**, after W2S (`inv`, the new R) and W2L (u and v, on the shared rows of inv):
- u² and v² (SQ5 and F5: a square and its fold);
- `Q_t = (x_Y + β^t x_X/2)² + (3/4)β^(2t) x_X²` for t = 0, 1 (FADDW, SQ5, AC4 adds the constant, F5), and
  `Q_2 = 3x_Y² - Q_0 - Q_1` (two FS), as `1 + β + β² = 0`;
- for each t: Q_t² (SQ5, F5), the two products Q_t²*u² and Q_t²*v² on the shared rows of Q_t² (the pair schedule of
  W2L), and from each product's six-word sum `Z` the leaf `Z - s'' + 2^256`, folded (SB6, FDX), where
  `s'' = β^t x_X + (x_Y + C0)` (FADDW), so that the fold's 2^256 ≡ C0 cancels.

The header holds `x_X, y_X, βx_X, β²x_X, x_X/2, βx_X/2, 3x_X²/4, 3β²x_X²/4`; a Y entry holds `y², y, x, 3x², x + C0`;
the six leaves overwrite the 256-byte entry in the order `c - 1 = 2t + [s = -1]`. The arithmetic is F's row schedule;
the rare carries and borrows of SHX, FDX, F5, FS and FADDW are taken by branches over fixups that give the inline
chains' results (3.4), and its only errors are the sites of 3.4.

### 3.2 Key pass, emit and counts

`keys` walks the entries downwards and takes each entry's six leaves. For each leaf it forms KEY's column F1 and the
next column's shifted sum s (3.3), then decides KEY's index j: with `σ = F0 + F1`, if F0 ≥ F1 (CMP, B.HS) it branches to
the block JH, where j = 2 if `F1 + σ` carries (ADDS, B.HS to JT) and j = 1 otherwise; else j = 2 if `F0 + σ` carries
(ADDS, B.HS to JT) and j = 0 otherwise, which goes on inline. The j = 0 path forms `c6 = F1 + (s >> 61)` and emits the
pair (x's, A's); JH forms c6 and `N' = ~((x >> 128) + A)` (ADDS, ADC, MVN) and emits (A's, N''s); JT forms N' (ADDS,
ADC, ADD of `s >> 61`, MVN) and emits (N''s, x's). The **emit** forms d1 (LSR, UBFX or EOR) and its cursor's address in
the table at port + 64, advances the cursor by 32 with LDADD, which returns the old one, and drops the leaf if that is
at or past the slot end in x12 (5 instructions); otherwise it forms w3 from d2 and writes the message with three STPs
(5). The leaf pointer, advanced 32 bytes per leaf by a post-indexed LDP, gives w3's high bits: entries are 256 bytes at
a 2^34-aligned base, so `pointer << 31` is `c<<36 | l<<39`. JH and JT add 5 to the key-path count in x18 and branch
back. **Meter.** The exact index is uniform on {0, 1, 2} for a leaf other than O, and KEY's index equals it outside a
set of probability at most `d_wj` per leaf (3.3), so a leaf takes JH or JT with probability at most `2/3 + 2d_wj`. A
node starts a batch only while its count is at most `5(cap_j - 6q2)`, `cap_j = ceil(2n0/3) + 2^42`. Two leaves of one
incidence are an invertible linear image of (X, Y), and leaves sharing one generator are independent given it up to one
excluded point, so a node's number of such leaves has mean at most `n0(2/3 + 2d_wj)` and variance at most `(2/9 +
2d_wj)n0 + 1` (the dependent pairs' covariances are O(1/r) each; 6.4's meter row, while no batch of the node is
defective). Executed counts per incidence (fixed part in parentheses): `fwd` 99 (5) and `bwd` 1403 (5), each plus 3 or 4
per fixup (3.4), and `keys` 266 for six leaves on the j = 0 path, including their six emits of 10 (8), plus the count it
adds to x18; fe_inv (`secp_fe.s`) is 23646. Per leaf `keys` executes two LDPs, F1 (11), the column (15), the index's
ADD, CMP, B.HS, ADDS and B.HS, the ADD of c6 and the emit: 44, plus the loop's SUB and TBZ per six leaves; a leaf on JH
or JT executes 5 more. Kernels may overwrite x0-x29 and the flags and preserve sp; `bwd` and fe_inv save and restore
x30, which `fwd`, `keys` and the transport kernels do not write; `fwd` and `bwd` return their fixup word in x0 (3.4);
the driver keeps no registers across calls and loads the key-path count into x18 for `keys`. `bwd` keeps R in x25..x28
from pass to pass; `bwd` and `keys` end with a TBZ on bit 33 of the entry pointer: the 2^25 entries of 256 bytes fill
the lower 2^33 bytes of their 2^34-aligned block, so bit 33 is clear at every entry and set one entry below the block.

### 3.3 The class label (the index, the column, JH and JT)

**Fractions and the label.** For v in [0, p) let `T(v) = floor(2^128 v/p)` and `comp(v) = floor(2^100 v/p) = T(v) >>
28`. For a leaf value x in [1, p) let `v_j = β^j x mod p`; the canonical member `v* = min v_j` has index j = argmin, and
the class's **label** is `(comp(v*), comp(βv* mod p))`: the member `β^m x` has index `j - m` and the same pair. As
`2^128 v/p = v/2^128 + C0*v/(2^128 p)` with the second term in [0, 2^-95), `T(v) = (v >> 128) + e0(v)`, e0 in {0, 1};
and `comp(v) = v >> 156` unless `v mod 2^156 > 2^156 - 2^34`.

**The index.** This is F's canonical-index row without its windows. `G = floor(β*2^384/p)` has limbs G0..G5 (x16,
x7..x11), F0 is the top limb of x, `F1 = Σ_i [lo(x_i*G_(5-i)) + hi(x_i*G_(4-i))] mod 2^64` (column 5 of `x*G` without
its carry-in), `σ = F0 + F1`, `d = F0 - F1` and `e = min(F0, F1) + σ`. Then j = 2 if e carries, else 1 if F0 ≥ F1, else
0; 3.2's branches decide exactly this without forming the minimum. The **window set** is where (W1) `F1 >= 2^64 - 16`,
(W2) `σ >= 2^64 - 17`, (W3) `e mod 2^64 >= 2^64 - 33` or (W4) `d < 16`; it only describes where j may be wrong.

**Exactness of the index.** Let `V_i = 2^64*v_i/p` for `v0 = x`, `v1 = βx mod p` and `v2 = β²x mod p`, which are
distinct and sum to p or 2p. Then `F0 = V0 - a` with `0 <= a < 2`, and `a < 1` unless `x mod 2^192 > 2^192 - C0`. Bits
320..383 of `x*G` equal `F1 + c` and `floor(V1 - ε)` mod 2^64, with `0 <= c <= 8` (the terms left out of F1 sum below
8.01*2^320 mod 2^384) and `0 <= ε < 2^-64` (as `0 <= β*2^384/p - G < 1`). So outside (W1) `F1 = V1 - b` with `0 <= b <
10`, and outside (W1) and (W2) `2^64 - σ = V2 + a + b`. Hence j = 2 gives `V0, V1 >= V2`; j = 1 gives `V0 - V1 > 16 -
10` and `V2 - V1 > 34 - 2 - 20` by (W4) and (W3); j = 0 gives `V2 - V0 > 34 - 4 - 10` and `V1 - V0 >= 1 - a > 0` when `a
< 1`. Outside the window set j is therefore `argmin v_j`, except when j = 0 and `x mod 2^192 > 2^192 - C0`. **A wrong
index.** Where j is wrong the leaf emits another member's pair: a missed collision charged as the row `d_wj` of 6.4, or,
if it equals another class's label, a false candidate, which the mass bound below covers. **Rarity.** (W1) puts v1
within 2^196 of 0 or p; (W2) and (W4) put `(1±β)x mod p` within 2^197 of 0; (W3) puts `(2+β)x` or `(1+2β)x mod p` within
2^198 of 0. As these maps are bijections, at most `13*2^197` values of x lie in the window set and at most `2^64*C0`
more give the wrong index; a uniform point's x takes each value with probability at most `2/(r-1)`, so `d_wj <=
(13*2^197 + 2^64*C0)*2/(r-1) < 2^-54.29` per leaf.

**The column (15).** Column 4 of `x*G` has the eight words `lo(x_i*G_(4-i))` and `hi(x_i*G_(3-i))`. The kernel adds
`hi(x_0*floor(G3/8))` (x15 = G3 >> 3) and the other seven words shifted right by 3 (ADD with LSR #3) into s: each term
is below 2^61, so `s < 2^64`. Then `A = 2^64*F1 + 8s mod 2^128`, with top word `c6 = F1 + (s >> 61)`. **Bound.** As
`x*G/2^384 = xβ/p - xδ/2^384` with δ in [0, 1), `floor(2^128 {x*G/2^384})` is `T(v1)` or `T(v1) - 1` mod 2^128, and it
is columns 4 and 5 of `x*G` with the carry from columns 0 to 3. A leaves out that carry (columns 0 to 3 have 1, 3, 5 and
7 words, so their carries out are at most 0, 2, 4 and 6) and the bits the shifts drop: `Σ (t mod 8)` over seven words,
at most 49, and `hi(x_0 G3) - 8 hi(x_0 floor(G3/8))`, at most 12 as `G3 mod 8 = 5`. Hence `T(v1) - A ≡ e1 (mod 2^128)`,
e1 in [0, 68], for every x. **The third fraction.** As `v0 + v1 + v2` is p or 2p, `T(v2) ≡ -T(x) - T(v1) - k (mod
2^128)` with `k = ceil(f0 + f1)` in {0, 1, 2}, f0 and f1 the fractional parts of `2^128 x/p` and `2^128 v1/p` (the
fractional part of `2^128 v2/p` is in [0, 1) and the three multiples sum to a multiple of 2^128). `N' = ~((x >> 128) +
A) = -(x >> 128) - A - 1 mod 2^128`, and `x >> 128 = T(x) - e0(x)`, so `N' - T(v2) ≡ e0(x) + e1 + k - 1`, in [−1, 70].

**The paths.** A leaf of index j emits `(x >> 156, A >> 28)` for j = 0, `(A >> 28, N' >> 28)` for j = 1 and
`(N' >> 28, x >> 156)` for j = 2: its fields are each component's top 64 bits (K2, K3) and low 36 (d1, d2), so for j = 0
K2 = x3, d1 = x2 >> 28, K3 = c6 and d2 = bits 25..60 of s; for j = 1 K2 = c6, d1 = bits 25..60 of s, K3 = N''s top word
and d2 its bits 28..63 (EOR of the sum's bits 28..63 with 2^36 - 1); for j = 2 K2 = N''s top word, d1 its bits 28..63,
K3 = x3, d2 = x2 >> 28. With the index exact, `v_j = v*`, so the pair is the class's label unless a computed value's
error carries it across a block of 2^28 or x's own bits differ from comp: the **label-error set** D is where `T(v1) mod
2^28 < 70`, `T(v2) mod 2^28 >= 2^28 - 71` or `= 0`, or `x mod 2^156 > 2^156 - 2^34`. Outside D a leaf whose index is
exact emits its class's label. As each value of T is taken by at most 2^128 values and `x -> v_j` is a bijection, D has
at most `142*2^228 + 2^134` values of x, each of probability at most `2/(r-1)`: the row `d_lab < 2^-19.85` of 6.4, a
missed collision like `d_wj`; no leaf is dropped.

**Mass bound.** secp256k1's Frobenius is `π = A + Bω` with `B = 0xe4437ed6010e88286f547fa90abfe4c3` and
`A = 0x114ca50f7a8e2f3f657c1108d9d44cfd9`: `A² - AB + B² = p`, `B*β ≡ A - B` and `(A - B)*β ≡ -A (mod p)` (checked by
the run and the certificate). So `Λ = {(v, u) ∈ Z²: u ≡ βv (mod p)}` has the basis `(B, A - B), (A - B, -A)` of
determinant -p, and as `A > B > 0` its half-open parallelogram P spans [0, A) in v and (-A, A - B) in u. For each point
λ of Λ in a box `[v0, v0 + a) × [u0, u0 + b)` of integer sides a and b, `λ + P` lies in `[v0, v0 + a - 1 + A) × (u0 - A,
u0 + b - 1 + A - B)`; these translates are disjoint, as P tiles the plane under Λ, and each has area p, so the box holds
at most `(a + A)(b + 2A - B)/p` points of Λ. Fixing f of a label's 200 bits leaves `2^(200-f)` pairs of component
values; a component value is an interval of at most 2^156 values of v (of u = βv mod p for the second), widened at each
end by 72 units of `p/2^128` for a computed component (one that wraps mod 2^128 may be two end pieces of [0, p), which
count as one interval as `pZ² ⊂ Λ`). So on each path at most `2^(200-f)(2^156 + w + A)(2^156 + w + 2A - B)/p < (1 +
2^-19)2^-f p` values of `v_j`, and of x, give those bits, with w the widening (the certificate checks the factor),
whatever index KEY returns. With three paths and two points per x, the emitted pair takes given values on f of its bits
with probability at most `KB(f) = 6*2^(256-f)(1 + 2^-19)/(r-1)`.

### 3.4 Arithmetic sites

**Steps that can lose a carry.** A product on the rows of a weak b (below 2^256) uses `B_0 = b` and row shifts
`B_(k+1) = SHX(B_k)`; RF and RA sum `Z = Σ a_k*B_k < 2^322` exactly (each high word is at most 2^64 - 2), and a fold
(FDX) computes `t = (Z mod 2^256) + floor(Z/2^256)*C0` and keeps t mod 2^256. SQ5 forms the square of a weak a as `V =
L + H*C0` exactly (`a² = L + 2^256*H`; V's top limb is at most C0 + 1), AC4 may add a constant below p (top limb at
most C0 + 2), and F5 folds V the same way. SB6 forms `Z' = Z - s'' + 2^256` exactly (top word at most 4). Only these
steps drop a carry or a borrow: SHX's last ADCS, the final ADC of FDX and of F5, and FS when its subtrahend is a
fold's output at or above p. Without a dropped carry each returns a value below 2^256 congruent to the exact one:
- **row shift** of `B = (b0, b1, b2, b3)`: the carry is lost iff `(B mod 2^192) + floor(b3*C0/2^64) >= 2^192`, which
  needs `B mod 2^192 > 2^192 - 2^33`: at most 2^97 rows. Otherwise `SHX(B) = B*2^64 - b3*p`.
- **fold** of `Z < 2^322`: lost iff `t >= 2^256`, which needs `Z mod p = t - p` in `[C0, 2^66*C0)`; the fused fold of
  `Z' < 2^322 + 2^256` likewise needs the leaf in `[C0, (2^66 + 2)*C0)`.
- **square fold**: lost iff its sum is at least 2^256, which needs the true square in `[C0, C0² + 2C0)`, or Q_t in
  `[C0, C0² + 3C0)` with the constant added.
- **FS** `a - b` with `a < p` and `b < 2^256`: wrong iff `b - a > p`, which needs `b >= p`; the subtrahends Q_0 and Q_1
  are fold outputs, at or above p only when the true value is below C0.

A fold returns a value `>= p` only when the true value is below C0. That matters only where a canonical value is needed:
the leaf entering KEY and the key paths, and FS's subtrahend (above). Every other operation is exact on its inputs: RF
and RA; SB6; the label's column sums and path additions (3.3); FS
and FADDW on canonical inputs; FADDW of a canonical value and `x_Y + C0` (the sum is
below `p + 2^256`, so one conditional subtraction of p leaves it below 2^256); KEY's integer arithmetic (3.3 covers its
index); the emit; and F's `fe_mul_f` and `fe_sqr_f`, which fold the 512-bit product, fold its fifth limb again and add
C0 exactly when that last sum carries past 2^256, so fe_inv (255 squarings, 15 multiplications) is exact for every input
below 2^256. Generator construction, the table entries, the header values, the resolver and `[x]P = Q` run on
`fe_mul_f`, `fe_sqr_f` and fe_inv, each coordinate brought below p (an ADDS chain of C0 detects `t >= p`, CSELs keep
`t - p`) before it is subtracted, compared or stored, with exact arithmetic mod r for the solve. So they have no sites;
the tables hold exact canonical values, and `x_Y + C0 < 2^256` as `x_Y < p`.

**Distribution.** Given X ≠ O and every other generator, `Y_l` is uniform on the r - 1 points of G \ {O}. A value of
x_Y comes from at most 2 points, and of y_Y from at most 3. So `D_l = x_Y³ - x_X³`, `acc_l = acc_(l-1)*D_l`,
`R_l = 1/acc_l` and `inv_l = 1/D_l` take each value with probability at most `6/(r-1)` (the rest are bijections while
no `D_m = 0` for m ≤ l, the degenerate case below); `u = 1/(y_Y + y_X)` and `v = 1/(y_Y - y_X)` (as
`D = (y_Y - y_X)(y_Y + y_X)`) at most `3/(r-1)`; u² and v² at most `6/(r-1)`; Q_t, quadratic in x_Y, at most `4/(r-1)`;
Q_t² at most `8/(r-1)`; and each leaf at most `2/(r-1)`, as `φ^t(X) ± Y_l` is uniform on r - 1 points. A row
`B_k ≡ v*2^(64k)` determines its operand v, so at most 2^97 values of v make a given row shift carry.

| site, at incidence l | runs per batch | fresh generator: operand | bad values | points per value |
| --- | --- | --- | --- | ---: |
| fwd: row shifts of `acc_(l-1)`, l ≥ 1 | 3(q2 - 1) | `Y_(l-1)`: `acc_(l-1)` | 2^97 | 6 |
| fwd: fold of `acc_l`, l ≥ 1 | q2 - 1 | `Y_l`: `acc_l` | 2^66*C0 | 6 |
| bwd: row shifts of `R_l` | 3q2 | `Y_l`: `1/acc_l` | 2^97 | 6 |
| bwd: fold of `inv_l = R_l*acc_(l-1)`, l ≥ 1 | q2 - 1 | `Y_l`: `1/D_l` | 2^66*C0 | 6 |
| bwd: fold of `R_(l-1) = R_l*D_l`, l ≥ 1 | q2 - 1 | `Y_(l-1)`: `1/acc_(l-1)` | 2^66*C0 | 6 |
| bwd: row shifts of `inv_l` | 3q2 | `Y_l`: `1/D_l` | 2^97 | 6 |
| bwd: folds of u and v | 2q2 | `Y_l`: u, v | 2^66*C0 | 3 |
| bwd: square folds of u² and v² | 2q2 | `Y_l`: u², v² | C0² + 2C0 | 6 |
| bwd: square folds of Q_0 and Q_1 | 2q2 | `Y_l`: Q_t | C0² + 3C0 | 4 |
| bwd: FS forming Q_2 | 2q2 | `Y_l`: Q_0, Q_1 | C0 | 4 |
| bwd: square folds of Q_t², t = 0, 1, 2 | 3q2 | `Y_l`: Q_t² | C0² + 2C0 | 8 |
| bwd: row shifts of Q_t² | 9q2 | `Y_l`: Q_t² | 2^97 | 8 |
| bwd: fused folds of the six leaves | 6q2 | `Y_l`: the leaves | (2^66 + 2)*C0 | 2 |


A run fails with probability at most its bad count times its points per value over r - 1. At l = 0 three row shifts
and three folds cannot fail, and are not in the table:
- fwd's rows of `acc_(-1) = 1` are 1, 2^64 and 2^128, and its fold of `D_0*1` has `Z < 2^256`;
- bwd's fold of `inv_0 = R_0*1` has `Z = R_0 < 2^256`;
- bwd's last update `R_0*D_0 ≡ 1` folds to 1 or p + 1, as nothing else congruent to 1 lies below `2^256 + 2^66*C0`. It
  is left in x25..x28 after the last pass and never read.

The degenerate cases are X = O (probability 1/r) and `D_l = 0`, that is `Y_l = ±φ^t(X)` for some t (at most `6/(r-1)`
per incidence), which zeroes acc and spoils the batch; when `D_l ≠ 0`, `Q_t = Π_(t'≠t)(x_Y - β^(t') x_X) ≠ 0`. The
first site of a batch whose predicate holds has exact operands, so a batch has a defect with probability at most
`d_batch`: the table's sum plus `6q2/(r-1) + 1/r` (`SITES` in `certificate.py`; 6.4).

**Rare carries by branch.** In SHX the carry out of its first ADDS, in FDX and F5 the carry into limb 2, in FS the
borrow of its C0 subtraction and in FADDW the carry of its C0 addition are rare. Each macro branches over its **fixup**
(B.LO; B.HS for FS) unless that carry or borrow occurred. The fixup is the rest of the inline chain, two instructions
(ADCS and ADCS or ADC) or three for FS and FADDW, and `STADD x4, [sp]`. Where the branch is taken the skipped
instructions would add or subtract 0, so on every input every register equals the inline chain's, and so does every
deviation of the table above; only the flags differ, and every path from the join sets them before reading them or
returns (`fixups.flags_dead_at_every_join`; kernels may leave the flags changed). A fixup adds 3 instructions (SHX, FDX,
F5) or 4 (FS, FADDW) and C0 to the word at [sp], which opens a 16-byte frame pushed as zero by `fwd` (STR) and by
`bwd` beside its saved x30 (STP); each returns the word in x0, at most `46q2*C0 < 2^63` per call.
The driver adds x0/C0 to the node's fixup count N and starts a batch only while `N <= 2^45 - 46q2` (46 sites per
incidence: 18 SHX, 11 FDX, 7 F5, 4 FS and 6 FADDW). So N is at most 2^45 on every run, and section 7 charges `4*2^45`
per node.

**How often.** A fixup needs a word in a small set:
- SHX's carry needs the row B it shifts to have its low word above `2^64 - C0`, as the word it adds is below C0: at most
  `C0*2^192` rows.
- FDX's carry into limb 2 leaves its output T with `T mod 2^128 < F*C0 <= 2^66*C0` and `T >= 2^128`, where `F =
  floor(Z/2^256) <= 2^66` (for the fused fold too, as `Z' < 2^322 + 2^256`): at most `2^194*C0` words.
- F5's carry into limb 2 leaves `T mod 2^128 < 2^65`, as its fifth word is at most C0 + 2, so that word times C0 has a
  high word of at most 1: at most 2^193 words.
- FS's borrow leaves an output whose low word is at least `2^64 - C0`, and FADDW's carry one whose low word is below C0:
  at most `C0*2^192` words each.

While no site of the batch up to and including this one has deviated, the row or output is a word below 2^256 congruent
to the site's operand value (times `2^(64k)` for the k-th row), which therefore determines that value. So, given X and
every other generator, a site takes its fixup with probability at most its set's size times the points of `Y_l` per
operand value over r - 1, with the fresh generator and multiplicities of the table above; `y_Y ± y_X` have 3 points per
value, and `x_Y + β^t x_X/2` and `s''` have 2. At l = 0 the six sites the table leaves out cannot take a fixup either:
rows 1, 2^64 and 2^128 have a low word of at most 1, folds of values below 2^256 add nothing, and `R_0*D_0` folds to 1
or p + 1, neither of which is at least 2^128 with its low 128 bits below `2^66*C0`. Summed over a batch (`FIXUPS` and
`E_FIX` in `certificate.py`), the expected number of fixups a batch takes before its first deviation, and so of a batch
without a defect, is at most `E_FIX`, and over a node's batches at most `q1*E_FIX < 2^28.35`.

| fixup sites, at incidence l | runs per batch | words | points per value |
| --- | --- | --- | ---: |
| fwd: FS forming `D_l` | q2 | `C0*2^192` | 6 |
| fwd: row shifts of `acc_(l-1)`, l ≥ 1 | 3(q2 - 1) | `C0*2^192` | 6 |
| fwd: fold of `acc_l`, l ≥ 1 | q2 - 1 | `2^194*C0` | 6 |
| bwd: row shifts of `R_l` | 3q2 | `C0*2^192` | 6 |
| bwd: folds of `inv_l` and `R_(l-1)`, l ≥ 1 | 2(q2 - 1) | `2^194*C0` | 6 |
| bwd: FS of `y_Y - y_X`, and FADDW of `y_Y + y_X` | q2 each | `C0*2^192` | 3 |
| bwd: row shifts of `inv_l` | 3q2 | `C0*2^192` | 6 |
| bwd: folds of u and v | 2q2 | `2^194*C0` | 3 |
| bwd: square folds of u² and v² | 2q2 | 2^193 | 6 |
| bwd: FADDW of `x_Y + β^t x_X/2`, t = 0, 1 | 2q2 | `C0*2^192` | 2 |
| bwd: square folds of Q_0 and Q_1 | 2q2 | 2^193 | 4 |
| bwd: two FS forming Q_2 | 2q2 | `C0*2^192` | 4 |
| bwd: square folds of Q_t² | 3q2 | 2^193 | 8 |
| bwd: row shifts of Q_t² | 9q2 | `C0*2^192` | 8 |
| bwd: FADDW of `s''` | 3q2 | `C0*2^192` | 2 |
| bwd: fused folds of the six leaves | 6q2 | `2^194*C0` | 2 |

**The meter row.** A node stops early only if N exceeds `2^45 - 46q2` at the start of some batch. On the event that none
of its batches has a defect, N is that clean count, so by Markov's inequality a node stops with probability at most
`q1*E_FIX/(2^45 - 46q2) + q1*d_batch < 2^-16.6` (`d_f`, a node-local row of 6.4).

## 4. Transport and grouping

DIRECT_PROTOCOL.md gives each actual kernel formula and all completion costs. The full per-record charges are hop1=18; hop2=7+15+2+2/8=97/4; grouping=(6+2/8)+195/64=595/64, plus explicit range, cap, table, return and phase rows. The receiver's full-drop path is11, below the15 charged for every delivered record. No second-hop slot loss is charged because no such slot exists. Legacy bkt code is not called.

**Scan** (`scan_ranges`: the range layout of submissions 80e9e8af
and 1efcef88, the paired form of 892a3df8 and b00400fa). `prefix4` sets each d4 counter to its range's start and `pass4`
advances it per pointer it writes, so after pass4 the counter is the range's exclusive end. `scan_ranges` visits the
2^48 ranges in order and never compares across one. In a range of n >= 1 records, `n - 1 = 64q + 2r + o` with r < 32 and
o in {0, 1}, it loads the first record's K2, makes one single comparison (two LDRs, CMP, B.EQ) when o = 1, enters a
block of 32 pairs through a five-level TBZ tree (and a B when r is odd) for the r remaining pairs, and then runs the
block q times. A pair loads two pointers with one LDP and their K2s, compares the older adjacent pair with CMP and the
newer with `CCMP ...,#4,NE` (Z when either is equal) and branches once: 6 per two records, and 195 per 64 with the
block's end test (MOV, CMP, B.LO). Three key registers and four pointer pairs rotate through the block; the entry copies
the range's first record into every register a remainder can start from, so each comparison is between a record and the
one before it in the range, and as 32 mod 3 = 2 and 32 mod 4 = 0 one MOV realigns the keys for the next block, whose CMP
with the range's end (x18) stops it there, reading no pointer past the range. With range control a range without
candidates costs `25 + 6o + (r mod 2) + 6r + 195q`, which exceeds `(195/64)n` by `1405/64 + 189o/64 - 3r/32 + (r mod
2)`, at most 25.8125 (o = 1, r = 1): 26 per range is charged; an empty range costs 5. A hit re-tests both comparisons on
its out-of-line path (CMP, B.NE each) and appends each equal pair in order after its list-end test (CMP, B.HS, STP),
then branches back: 7, 8 or 11 for the first, the second or both, 4 for the single comparison, so at most 8 per listed
candidate; the attempt that finds the list full ends the scan. So a node's scan is at most `(195/64)n + 26*2^48 +
8(2^28 + 1) + 2`, the row of section 7 (S is a multiple of 64). A hop's receiving node executes nothing on arrival; it
reads the records in the next phase, with the LDPs of `h2s` (hop 1) or `bkt` (hop 2). Per-node costs are rows of section
7, and every address a kernel loads or stores lies in the regions of section 5.

**One local pass suffices.** A home's records share d1 and d2; the d3 buckets and the stable d4 pass make records with
equal (d3, d4) adjacent, grouping 144 label bits. The records of one class form one run unless an **intruder**, a record
with the same 144 bits that does not carry the run's label, sits inside (6.4); a run lies in one d4 range with one K2,
so its adjacent pairs are all listed unless the candidate list is full. `scan_ranges` lists every adjacent pair of a
range with equal K2, which agrees on d1, d2, K2 and d4, 184 label bits: besides a run's pairs it can list pairs of
different classes, which the resolver rejects (no unit relates their points) and 6.4 charges to the candidate, request
and relay caps through KB(184). Writing the members `u_i*Z` with Q coefficients `n_i`, a pair is useless iff `n_i/u_i =
n_j/u_j`, an equivalence relation, so a run with a useful pair has a useful adjacent pair.

## 5. Node layout and memory

DIRECT_PROTOCOL.md lists the complete fixed allocation F<2^52, including alignment padding, receiver context/window, relay2^34records and interleaved cursor/end table. The disjoint bulk regions are R1=7*2^60bytes,32*2^24*CAP3 full-record bytes, and8*2^24*CAP3 pointer bytes. Their final address E=17977400794953023488<2^64; total N*E<2^136. No bucket storage overlays live first-hop memory. Attempted hop1 cursors stay below2^63, bucket reservation temporaries and d4 counts do not wrap. The old2^38 mask/stride restriction is absent;32*CAP3 exceeds2^38 and is handled by explicit ends.

## 6. Success (no heuristic)

### 6.1 Model

Coins are uniform and independent, and Q is uniform. Each `Y_l` is uniform on G \ {O}; given R_g, the
`X_k = R_g + [a_k]P` are independent and uniform. So the X's are independent uniform points, independent of the Y's, and
nodes are independent; a leaf is uniform on G, and leaves of different nodes are independent. Given all generator
points, the `b_g` are still independent and uniform: for each value of b_g exactly one tuple `(a_k)` fits the points.

Γ is the set of unordered pairs of leaves from different nodes, `|Γ| = (S² - 2^72*n0²)/2`. For α = {L, L'} ∈ Γ let
`I_α = 1[L' ∈ class(L)]`; then `P(I_α = 1) = P1 = (6(r-1) + 1)/r²` and `μ = |Γ|*P1 = 0.494347007` (the run's `mu`).
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
all points, this fails with probability exactly 1/r. So `P(fail) <= P(W = 0) + E #{α : I_α = 1, α useless or not
detected}`.

### 6.4 Defect rows

A **node-local defect** of L (a batch defect of 3.4, a wrong class index, a label error, either meter) is a function
of node g's coins, and L' is uniform given them, so `P(I_α ∧ defect) <= (6/r)*P(defect)`; over Γ a row costs at most
`μ(1 + 1/r)` times its per-leaf rate.

**Caps and intruders.** Condition on α's four generators. A leaf with a generator outside them, a **free** generator,
is then uniform on at least r - 1 points, so at most `KB(f)` of its mass falls on given values of f label bits,
whatever index KEY returns (3.3's lattice count). The
12 leaves of α's two incidences are fixed; the six leaves of an incidence through one of α's generators and a free one
are functions of that free generator.
- **Pair slots and d3 buckets** count single leaves. A first-hop phase slot has mean at most `m_p = 6*ceil(q1/2)*q2*KB(36) + 1` from the free leaves; and a d3 bucket `m_3 = S*KB(96) + 1`; each cap leaves room for the at most 12 fixed leaves. Two leaves'
  indicators have covariance 0 without a common free generator or with a common free Y; at most `1/(4(r-1)²)` with a
  common free X (at most `KB(f)/(r-1)` when one of them also uses α's Y); at most `1/(r²(r-1))` times the two cells'
  sizes for two leaves of one incidence with both generators free, whose pair is uniform on r(r - 1) pairs (an
  invertible linear image of (X, Y), with Y ≠ O); and at most one leaf's mean for two leaves of one incidence through
  α's generators, at most60(qphi+q2) ordered pairs for a phase slot and60(q1+q2) for a full bucket. A count's variance is at most its mean plus `1+60(qphase+q2)*KB(f)`, with qphase=qphi for slots and q1 for buckets, and Chebyshev bounds the drop of a given record by `(m + 1 + 60(qphase + q2)KB(f))/(c - 12 - m)²`
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

  Hence the intruder rate is at most `ρ(S + 12)KB(144)` per colliding pair; the candidate list at α's home has mean
  `other_h <= ρ(S²*KB(72)*KB(184) + (2S + 20)*KB(184))` (pairs at the home equal on d1, d2, K2 and d4, the 184 label
  bits `scan_ranges` matches; L and L' are there); a source receives at most `other_g <= 2ρ*n0*S*KB(184)` requests; a
  relay holds at most `other_v <= 18ρ((S/2^36)*S*KB(36)*KB(184) + (2S + 20)*KB(184))` records (a request is one record
  and a reply at most 8). Markov with the caps less α's own entries (4 in the candidate list, 8 requests, 16 records)
  gives `d_cl`, `d_rq` and `d_rl`. Leaves of a defective batch elsewhere are not uniform; the whole-run row "a
  computational defect in any batch" excludes them.

| row (per leaf of the pair unless stated) | term in `certificate.py` | checked |
| --- | --- | --- |
| useless collision (per pair) | `1/r` | - |
| batch defect: the sites of 3.4, D = 0, X = O | `d_batch` | < 2^-126 |
| wrong class index (3.3): KEY's window set, or x mod 2^192 > 2^192 - C0 | `d_wj` | < 2^-54.29 |
| j = 1, 2 meter, Chebyshev with slack `2^42 - 6q2 - 1 - 2n0*d_wj` | `d_j` | < 2^-30 |
| fixup meter, Markov on the mean `q1*E_FIX` with the cap `2^45 - 46q2` (3.4) | `d_f` | < 2^-16.6 |
| the label-error set D (3.3) | `d_lab` | < 2^-19.8 |
| predetermined first-hop phase slot (Chebyshev, CAP1) | d_p once per witness leaf | exact rational certificate |
| d3 bucket (Chebyshev, `cap_3 = ceil(m_3) + 2^32`) | `d_3` | < 2^-30 |
| intruder (per pair), `ρ(S + 12)KB(144)` | `lam_int` | < 2^-14.6 |
| candidate, request and relay caps (Markov; a reply counts 8 records) | `d_cl`; `d_rq`; `d_rl` | < 2^-25.4; < 2^-27.9; < 2^-27.2 |
| a computational defect in any batch (whole run) | `g_def = 2^72*q1*d_batch` | < 2^-26.9 |
| two coin range rejections (whole run) | in `glob` | all of `glob` < 2^-26.9 |
| RNDR failure (whole run) | `g_rndr` = 0 under the rules' ideal coins (section 7) | - |

### 6.5 The bound

Use exactly success_terms(q1) in certificate.py. All unchanged Chen-Stein, arithmetic-site, wrong-index, two-fraction label-error, j-path meter, arithmetic-fixup meter, intruder, candidate, request and relay formulas are derived above. For each of the two witness leaves add one phase-slot loss d_p and one full-bucket loss d_3; do not add a nonexistent hop2-slot event. The unconditional union bound is
lo=1-exp(-mu)_upper-cs-mu*(1+1/r)*(1/r+lambda_intruder+2*(d_batch+d_wj+d_j+d_f+d_lab+d_p+d_3+d_cl+d_rq+d_rl))-global_loss.
The exact rational evaluation gives lo>=0.390000001 atq1=145289002, whileq1-1 fails that stated margin. No conditional cleanliness or empirical yield is substituted. Both meter caps are global across phases, as their probability and work formulas require.

## 7. Complete cold work ledger

All rows count scalar instructions across every node. Generator construction is669889 per generator including256doublings/mixed additions, affine conversion, bit/loop control and table/header creation: GEN_ITEM+GEN_CTRL+max(Y_ENTRY,X_ENTRY)=669889 in certificate.py. This removes only the unused old2^20allowance. Candidate cost retains the inherited1538117 bound (more than twice769058 itemized work) plus64 for changed source-column decoding. The same full scalar solve, point reconstruction, six-unit test, target verification and output remain charged.

Node setup2^20, header27, inversion23646+16, batch setup128, ideal coin draws64/scalar, source service2^10/request and output2^20 retain the explicit source itemizations; no quantity is omitted because it is cold or a failed path. Phase control increases to2^14/node and64barriers at2^10/node. Message delivery and ACK work, all fixed returns and aggregate exchange conservation are detailed in DIRECT_PROTOCOL.md. Scalar and generator draws have the same capped failure rules as above. All other numerical constants and site counts are unchanged.

| row | total instructions | preprocessing |
|---|---:|---|
| program delivery, node and port state | 4951760157141521099596496896 | yes |
| coins | 54052111560017664463916907888640 | yes |
| Y generators [c]P with their table entries | 106148149117036045692403126693265408 | yes |
| R_g = [b]Q and X = R_g + [a]P with their header values | 459616147145911973527422396878815232 | no |
| header constants | 127503895037480420769792 | yes |
| initial hop1 cursor tables | 973555660980002546832337706942464 | yes |
| second phase cursor reset | 973555660980002546832337706942464 | no |
| Q broadcast | 377789318635069175234560 | no |
| d4 counters zeroed | 1329227995784920595270289929925558272 | yes |
| d3 cursor table | 396140817293688170837364965376 | yes |
| forward | 2279174174504812288757006519075933782016 | no |
| inversion | 16234685446264599653580360336277504 | no |
| backward | 32299811740978061069950676863958585966592 | no |
| arithmetic fixups | 664613997892457936451903530140172288 | no |
| key pass without emits | 4742524038229269172735764654767210496000 | no |
| j = 1 and j = 2 key paths | 460543072416821751478853897424784588800 | no |
| batch set-up | 87821812911920748696571977138176 | no |
| hop 1, the emits, in-flight counts and deliveries | 2486371822990115696601929752490281009152 | no |
| hop 2, h2s, in-flight counts and deliveries | 3349702512047641928866166077460439891968 | no |
| hop1 final cursors, two phases | 11033630824414842909524522570874880 | no |
| prefix sums of d4 | 3323069989462294404626000520346075136 | no |
| local pass pass4 | 863323550520855520366055360457754017792 | no |
| scan_ranges | 455430168485668704138554841716190347264 | no |
| candidates | 1949876067909658129023600433148461056 | no |
| requests served | 2596148429267413814265248164610048 | no |
| phase-4 count exchanges | 41538374868278621028243970633760768 | no |
| barriers | 309485009821345068724781056 | no |
| phase control | 77371252455336267181195264 | no |
| output | 1048576 | no |

The sum is46944787124746853008406832072376781373440. The exact rational pow2_interval enclosure verifies total/1000<=2^(125142304476874/10^12). Preprocessing is included and below2^111units. Memory is independent of scalar score but fully bounded in section5. This is total computation, not elapsed time, throughput or an empirical attack cost.

### Itemized cold construction retained

- **Generator construction,669889 per point:** a complete Jacobian
  double-and-add, exceptional-case tests, two inversions with the affine conversion, all on the exact arithmetic of 3.4:
  `fe_mul_f` and `fe_sqr_f` at 115 and 82 instructions as executed, plus 4 and 3 for the call as in fe_inv, and an
  addition or subtraction with its operands brought below p at 27; 665,910 for the 256 doublings, 256 mixed additions
  and two conversions. The scalar's bit tests, the pointer step and the return add at most 6 per iteration, 1,600 in all, comprising1536per-bit control and64fixed setup/return. A Y's table entry (y², 3x², x + C0 and ten STPs, 239) or an X's mixed addition `R_g + [a_k]P` with its eight
  header values, y² and eighteen STPs (2,379) brings the total to at most 669,889, exactly the charged allowance.
- **Node set-up, 2^20:** a program of at most 2^16 bytes (4,050 kernel instructions, generator construction, resolver,
  phase control), 2^11 messages over 2 hops at 40, plus node state.
- **A message outside the hops, 40 per hop** (program, Q broadcast, phase 4): a broadcast message executes 8 (two LDPs,
  three STPs, ADD, SUBS, B.NE), plus 2 for its in-flight count and 6 for its delivery: 16.
- **Candidate, 1,538,117:** two requests and two replies of at most 8 records over 2 hops at 40, two exact additions
  with inversion, the six-unit test, the solve mod r and `[x]P = Q`, on the same exact arithmetic; 769,058, charged at
  twice that plus 1 (as in 1fbd26e8).
- **Phase-4 exchange:** the aggregate conservation proof in DIRECT_PROTOCOL.md replaces the old local message bound.
- **Barrier, 2^10 per node:** at most 5 waits of LDAXR, CMP, B.EQ and at most one WFE, B and re-check (8 each), at most
  5 messages of 14 (6, plus 2 for the in-flight count and 6 for delivery), DSB and 20 of set-up: 131.
- **Phase control:**2^14/node, including both phases and resident receiver setup.
- **d4 prefix sums, 10 per four counters:** two LDPs, two STPs, four ADDs, SUBS and B.NE; the RET is the row's +1.
- **Header 27** (G0..G5: six words of MOVZ and three MOVK, and three STPs),
  **batch set-up 128** (X's eight header values copied into H by 16 LDPs and 16 STPs, y_X² into registers, acc = 1,
  pointers, the index word, the slot end for `keys`, both meters' updates and tests and the loop: at most 119, itemized
  in `certificate.py`), **inversion call 16, request served 2^10, output 2^20.**
- **Coins, 64 per scalar:** at most two draws of 4 words `MRS Xt, RNDR` (FEAT_RNG), each with a B.EQ retrying once on
  the failure flag, plus a range test. RNDR is this machine's coin source: it returns 64 fresh uniform bits, the
  algorithm coins of the rules' probability space, and never reports failure, so the charged B.EQ is never taken and
  coins add no failure row.


## 8. Evidence and limitations of tests

probe_direct.py runs the unchanged frontier arithmetic/site/rare-path tests and class-label identities, the changed two-phase direct route with exact emit/delivery/receiver/reference comparisons, both narrow and wider slot configurations, source metadata and released ACKs, independent scanner paths including full-list exits, the small six-leaf toy, and exact success/cost/memory/cursor/exchange inequalities. The declared protected budget stays20seconds/256MB/oneCPU. The lower-level numerical source and complete proof, not small toy success, support the analytic bound. No historical long-run census is claimed as newly executed. Failed v1 intake due totalbytecap and earlier intermediate evidence remain preserved outside this final candidate.

## 9. Limitations and provenance

- **Model.** The bound is for the machine of §2.1, whose explicit blocking endpoints and resident receiver context are described in DIRECT_PROTOCOL.md. The sender's
  instructions (including those that form the address), each ordinary delivery (as `dlv`) or direct receiver including release ACK and the ports' in-flight counts are
  charged. Each hop is one link from the sender's port to its destination; a network that forwarded a message through
  further nodes would add their work at each extra hop, and that machine is not claimed.
- **Not executed:** the allowances of section 7 (generators, resolver, broadcasts, barriers, phase control, batch set-up
  with the driver's meter tests, and phase-4 exchanges), the in-flight counts, and the production attack, machine and
  links.
- **Provenance**, not a premise: package F (the original source author's 8f733826) for the row products (RF, RA, and the pair schedule of W2S and
  W2L) and FS's CSEL form, with kaankolcu's c8a6819e, its base, for the kernel macro set and `kernels/a64.py` (whose
  extensions its docstring lists), aburan28's 02feb284 for fe_inv and the F0 and F1 rows with their error analysis (the
  column of 3.3 extends F1), kaankolcu's 4595bad7 for SB6's subtraction of a stored offset, and aburan28's fce6e7ef for
  charging link-side counter updates as instruction-equivalents (the in-flight count); submissions 3828873f, 6ef3e7b3,
  d637fb0a, 9798f5f4 and d9e3dbf4 for STADD in `h2s` and `bkt`, paired scan comparisons, `prefix4`, the generator and
  candidate allowances and a claim on a decimal grid; submissions 472226d8, 4fdb8a9c and e1ecd7cc for a success target
  nearer 0.39 and a 12-decimal claim; submissions 80e9e8af and 1efcef88 for `scan_ranges`' range layout (K2-only
  comparisons within d4 ranges, a TBZ entry tree, the 2^28 and 2^29 caps; its first form was developed on 3828873f), and
  kaankolcu's 892a3df8 and b00400fa for its paired CMP/CCMP block of 64 records with the bound of 26 per range, taken
  unchanged; aburan28's 1fbd26e8 for charging a candidate at twice its itemized 769,058 plus 1; submission f9d5a673
  (the original source author's) for the LDADD cursors and the mask-tested buckets, and b3716a20 (the original source author's) for the six-leaf `fwd` and `bwd`, the
  site table and the success argument, all unchanged here but for the key, the scan and the fixups; kaankolcu's 807e4bf5
  for rare exits whose stubs run the original carry instructions and are charged, with a9297fb0's per-walk meter, which
  our 8f733826 used in SH and FD and our 34a30cf5 in fe_sub; our a3ca3235, the base (its fraction key formed the third
  fraction from x and the second; this package keeps its scanner port and charges), and our 502d7155, eb0aea59 and
  5b6cd013 for the windowless KEY, secp256k1's Frobenius coordinates A and B (here the lattice basis of 3.3), a class
  key with the range scan, and the metered fixups; Galbraith, Wang and Zhang for X + Y and X - Y from one inverse;
  Arratia, Goldstein and Gordon for Chen–Stein. The six-leaf sharing of 3.1 and the fixups of 3.4 are ours (b3716a20,
  5b6cd013); new here is the two-fraction label of 3.3: the shifted-sum column, the one-word complemented N' with its
  bound, the label-error set D and the lattice count behind its mass bound.


New integration provenance: the immediate public two-fraction source is bbe1cd85-58b0-4100-88bc-a095f03ae1e6. Direct-bucket ACK, memory and two-phase construction were developed here and in our P256v18 submission3f2ecd82-ce72-4111-bb03-fffb05041f20. The alternative orbit-min v8 key is not stacked onto the incompatible two-fraction key. Actual tools are OpenAI dot direct preparation, bounded Python/AArch64 interpreter execution and internal read-only mathematical audits. Underlying drafting model variant was not exposed; external reviewer is not sole author. Changed source requires fresh exact-version review and user submission instruction.


The final generator-control refinement follows public7e8d13c7-d6e9-49e7-8a2d-e4c33e0aaba1. A fixed four-limb/64-bit loop uses at most6instructions per bit for bit extraction/test, calls and loop control, plus64for four limb loads, outer-loop setup and return;256*6+64=1600. Scalar arithmetic and table generation remain in their separate itemized rows. Public acceptance is not the proof of this bound.
