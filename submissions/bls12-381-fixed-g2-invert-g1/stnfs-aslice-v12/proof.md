# STNFS in F_{p^12} with a packed-row slice sieve, certified thresholds and an exact fundamental domain

Track: `bls12-381-fixed-g2-invert-g1` (fixed-argument pairing inversion, FAPI-2). Unit: one dense Fp12
multiplication = C = 20,000 scalar AArch64 instructions (`tracks/classical-group-work-v2.json`, `tracks/docs/COST_MODEL.md`).

## Claim

For every production instance (uniform nonidentity Q in G2, uniform z in GT including z = 1), the algorithm below
returns P in G1 with e(P, Q) = z (P = infinity exactly when z = 1) with success probability at least **0.39** over the
uniform challenge and the algorithm's coins (section 8: >= 0.95 under the stated heuristics), using total classical
computation

  **W <= 2^122.343600 target-group-operation equivalents.**

W charges relation collection (every sieved cell, record, threshold block, setup and cofactorisation), the canonical
list, the factor base, filtering, linear algebra, Schirokauer maps, both individual-logarithm descents at their
metered cap, verification and randomness. There is no preprocessing or advice beyond this charged work; nothing is
amortised over other targets. Peak memory is reported in section 12 (not scored). Score: **122.343600**
(= 1223436/10000), the full-ledger logarithm 122.343555 rounded up.

This is an analytic cost claim supported by exact-norm sampling of the production region (regenerated in the sandbox
from per-index seeds), production-width kernels checked against brute force and measured in executed AArch64
instructions (`KERNELS.md`, `measurements12.txt`), a finite Euler-product model of the threshold loss calibrated
against exact smooth-number samples, a filtering simulation with the production occurrence profile, and an explicit
descent dynamic programme. `evidence12.py` recomputes every number of the ledger and prints `total_log2W_full_ledger`;
its parts are referenced below as [k].

### What changed since stnfs-aslice-v10 (rejected)

The v9/v10 reviews made three concrete objections and the v11 workbench found four more defects. Each is replaced
in this package; the honest price of the replacements is ~1.3 bits:

| objection / defect | v12 |
|---|---|
| global offsets of ~2^71 bytes, 64 fixed segments, 32-bit record counters | `slice12.c`: separately allocated segments <= 2^62 bytes, segment-local cells, 64-bit counters, runtime segment count; checked on every cell for 1-10 segments (4.3) |
| setup6 truncated 2^126/p to 64 bits; the check used `%` | one Barrett routine for every modular step, validated at every bit length 3-63 and exercised by the brute-force checks |
| scan priced per point, bucket bytes 14.9x larger | **packed rows**: one 4-byte cell per point, dense per-line suffix metadata, no bounding-box bytes |
| unpriced "norm at the block centre" thresholds | `nbound12.c`: certified interval bounds per block of <= 64 cells, evaluation counts and price measured; retention uses the lower bound, survivors the width (4.4) |
| survivor prime recovery unspecified | records (segment-local cell, global 64-bit ideal id), sentinel marking, one recovery pass per segment (priced per record) |
| finite unit test for the fundamental domain | canonical a = least (HNF key of aO_h, Q_a, lex) of its unit orbit: exact, no unit certificate (4.2); 2^90 units |
| Poisson "k primes above X2" keep factor | finite Euler-product retention of the actual threshold, exact-sample calibrations (4.6) |
| ECM cost a formula | measured 5-limb Montgomery multiplication (`ecm12.c`) (4.5, 7) |

v10 inherited NBAR = 1.947 degree-1 ideals per prime; [8] measures 1.996 below 2^12 and the model charges 2.0.

## 1. Reduction (admissible per rules.md)

Let g = e(P0, Q) with the pinned convention (py_ecc `pairing(Q, P0)`). Q is a nonidentity point of the prime-order
group G2 and e is nondegenerate, so g generates GT = mu_r < F_{p^12}^* (r prime, 255 bits). The challenge satisfies
z = g^a for a unique a mod r; P = [a]P0 answers it since e([a]P0, Q) = g^a = z, and z = 1 gives a = 0 and
P = infinity. The problem is one prime-order DLP in F_{p^12}^*, solved as two individual logarithms (z and g) over
one relation-collection/linear-algebra phase: a = log z / log g mod r.

## 2. Baselines

GT is a prime-order group whose only cheap automorphisms are the 12 Frobenius powers, so orbit-12 Pollard rho needs
sqrt(pi r / 24) = 2^125.96 group operations. No baseline is qualified on this track.

## 3. The tower, the polynomials and the automorphism

Special Joux-Pierrot tower for the BLS12-381 seed u = -0xd201000000010000 (|u| = 2^63.58):

- K_h = Q[y]/h with **h = y^6 + y^3 + y^2 - y - 1**. [0c] certifies exactly (`ideal_canonical.field_certificate`):
  Disc(h) = 58673 = 23 x 2551 is squarefree and h is irreducible mod p (Rabin), hence **O_h = Z[y]** and h is
  irreducible over Q. Signature (2, 2).
- **f = y^6 3p(x^2/y) mod h = x^12 - 2y x^10 + 2y^3 x^6 + y^5 x^2 + (1 + y - y^2 - y^3)**, degree 12 in x over O_h
  (3p(X) = X^6 - 2X^5 + 2X^3 + X + 1 is the BLS12 characteristic polynomial); **g = x^2 + |u| y = x^2 - u y**.

**The tower is F_{p^12}**: [0b] verifies in F_p[y]/h that h is irreducible mod p and that u y is a non-square in
F_{p^6} ((u y)^((p^6-1)/2) = -1; N(u y) = u^6 h(0) = -u^6 is a non-square in F_p because p = 3 mod 4), so g is
irreducible over F_{p^6} and F_p[y, x]/(h, g) = F_{p^12}; [0] checks exactly that Res_y(h, f(x^2 -> u y)) = 0 mod p,
so f and g share a root generating F_{p^12}.

**Automorphism.** x -> -x (f and g are even in x) is the Frobenius^6 conjugation: it maps the relation of (a, b) to
that of (a, -b) and log sigma(P) = -log P mod r, so candidates and unknowns are both quotiented by sigma (aut 2).

alpha (log2): the package charges **alpha_f = 5, alpha_g = 1**, conservative roundings of the Guillevic-Singh
`alpha_tnfs_2d` values recomputed in v10 (4.868, 0.844; a larger alpha means fewer relations).

## 4. Relation collection

Candidates are pairs (a, b) in O_h^2 (degree < 6 integer polynomials). A candidate is a relation when the ideals
(a - bx) O_f and (a - bx) O_g are both B-smooth: every prime ideal dividing them has norm <= B (exact smoothness, no
large primes). For a primitive pair (aO_h + bO_h = O_h) every such prime ideal has **relative degree 1 over K_h**
(a - bx in P with x generating a proper residue extension forces a, b in P cap O_h), so the **factor base** is the set
of prime ideals of relative degree 1 over K_h with norm <= B on both sides. Its K_h-degree-1 part ((p, y - s) with a
root of f or g mod p) is the sieved part and holds (NBAR/2) Li(B) sigma-orbits; the part over K_h-primes of degree
f >= 2 (norm p^f <= B) has fewer than 2 Li(B^(1/2)) < 2^31 ideals, is never sieved (a point lies in such an ideal with
probability p^-f) and is recovered at cofactorisation. The primes above 2 are of this kind: h has no root mod 2
(h(0) = -1, h(1) = 1; h = (y^2 + y + 1)(y^4 + y^3 + 1) mod 2), so K_h has primes P2 (norm 4) and P2' (norm 16); on the
g side g = x^2 + |u| y = x^2 mod P2, P2' (|u| = 2^16 x odd) gives one ramified relative-degree-1 prime of norm 4 and
one of norm 16, and N_g is even exactly when a in P2 or P2' (29.7% of pairs); on the f side there is none, so N_f is odd
for primitive pairs. O_h[x] is not 2-maximal on the g side (x/2^8 is integral): the factor-base construction computes
the maximal order locally at the index primes (p = 2 explicitly; every sieved prime passes Dedekind's p-maximality test
once, inside the 2^72 factor-base budget) and valuations at their primes are taken in that order (CADO-NFS's "bad
ideal" treatment). The yield model (4.6) credits a B-smooth norm pair whatever its small-prime part, consistently with
this factor base.

### 4.1 Region and counts (exact)

`config.json` holds the two exact rational positive-definite forms Q_a and Q_b (the yield-weighted block-diagonal
ellipsoid of v10, invariant under a -> -a and b -> -b separately). With A^2 = 2^21 and theta = 1/2: eligible a are
**a5 > 0, Q_a(a) <= theta A^2**; the half-slice of a is **b5 >= 1, Q_b(b) <= A^2 - Q_a(a)**. Continuous volumes:
eligible a 2^64.850, largest half-slice 2^67.472 cells ([3]). `lattice_count_bounds.py` proves with exact rational
Machin pi and sqrt/determinant intervals that for every radius >= 1024 (all radii here) the integer count of a
half-ellipsoid is within [0.98, 1.02] of its half-volume for each metric; the model multiplies counts by **2601/2500**
when they are work and by **2401/2500** when they are yield ([3] recomputes the factors).

### 4.2 Unit orbits and the fundamental domain (exact, no unit certificate)

For a unit eps of O_h, (eps a, eps b) is again a candidate and eps(a - bx) generates the same principal ideal as
a - bx in O_f and O_g, so both give **the same relation** (the Schirokauer-map entries differ by Lambda(eps), and eps
in F_{p^6} has r-logarithm 0: a trivial unit relation). Unit multiples are duplicates and must not be credited
(v1-v4 credited them, an 11x over-count; v5-v10 removed them with a finite unit test whose completeness the v9 review
did not accept). v12 uses an exact criterion:

**Canonical a.** The principal ideal aO_h has a unique row Hermite normal form key (O_h = Z[y] by [0c]); two eligible
a, a' have the same key iff aO_h = a'O_h iff a' = eps a for a unit eps. Among the eligible a of one key, the
**canonical** one is the least in the order (Q_a(a), lexicographic a): exactly one eligible a per unit orbit meeting
the eligible set, with an exact tie-break; nothing about the unit group, its regulator or embedding sizes is
assumed. The sieve runs one job per canonical a.

**One representative per orbit of <units, sigma>.** -1 is a unit (a5 > 0 selects the sign); sigma acts as b -> -b
and b5 >= 1 selects one of (a, b), (a, -b) (b5 = 0 is skipped by the kernel: a measure-zero loss). Hence every
relation orbit has exactly one point in F = {(a, b) : a canonical, b in the half-slice of a}, and the relation count
of F counts distinct orbits. Pairs with aO_h + bO_h != O_h ("non-primitive") are excluded from the yield (their
relations are multiples of smaller ones); the flag is computed exactly by HNF in the pool (2,227 of 2,401 canonical
rows are primitive).

**Price of the canonical list.** Enumerate the a-box (|a_i| <= 2907 by diagonal dominance, < 2^75 vectors, 2^14
instructions each for the exact quadratic test); key the < 2^66 eligible ones (|N(a)| < 2^106; HNF with an implicit
D Z^6 generator set keeps intermediates within 213 bits; < 2^34 instructions per key with bit-serial 4-limb arithmetic);
sort the keyed 512-byte records by two-way external mergesort in <= 2^62-byte segments with two-word ranks (66 levels,
2^14 instructions per comparison and move); one pass emits the minimum of each key group. The sum is below
20,000 x 2^90: the ledger charges **2^90 units** and 2^78 bytes of staging. (The pool's flags come from a different,
exhaustive oracle: complete enumeration of the vectors of aO_h with Q_a(v) <= Q_a(a) after an exact metric LLL, equal
absolute norms identifying unit associates; [2] recomputes them.)

### 4.3 The slice kernel `slice12.c` (production-width, packed rows)

One job per canonical a sieves the half-slice S(a) = {b : Q_b(b) <= A^2 - Q_a(a), b5 >= 1}:

1. **Layout.** Lines k = b5 = 1..K (K <= 2047) in segments of 2^gs lines; a **row** is the set of b0 of a fixed suffix
   (b1..b4, k); rows are packed consecutively inside their segment, so a cell index is a segment-local 64-bit integer.
   Per line a dense metadata array over the 4-dimensional suffix box (L4 <= 2^48.33 int64 entries: the cell of the
   row's b0 = 0, or poison for an empty row) gives the cell of a point as meta[k][idx4] + b0, idx4 affine in the walk
   coordinates; the fold b -> -b maps idx4 to L4 - 1 - idx4 and b0 to -b0. A cell holds two uint16 sums (f, g) in
   1/16 bit (16 x 1302 < 65535: no overflow).
2. **Line path (2^8 <= p < 2^20).** For each sigma-pair over (p, y - s), each line k and each of the two cosets (the
   ideal at +k; the ideal at -k folded, i.e. the conjugate's coset at +k), the 5-dimensional coset of
   L5(p, s) = {b' : sum_{i<5} b_i s^i = 0 mod p} inside the line's ellipsoid is enumerated by a Schnorr-Euchner walk
   with direct updates (+ floor(16 log2 p), computed exactly once per prime at factor-base construction). The coset origin is recentred exactly (integer Babai rounding in the
   reduced basis, wrapping uint64 arithmetic whose exact result is small), so walk coordinates are O(1) for every p.
3. **Record path (2^20 <= p < 2^63 = X2).** One 6-dimensional walk of the recentred coset of L(p, s) over the full
   ellipsoid per sigma-pair with affine (idx4, b0, b5) states; a hit with b5 < 0 is folded, b5 = 0 skipped, and the
   record (segment-local cell uint64, **global ideal id uint64**) is appended to its segment's retained store (64-bit
   counters; grown, never wrapped). The **apply pass** of a prime batch adds the batch log to the recorded cells.
   Sieved ideals have p < 2^63, so ids (< 2^60) and cell indices fit 64 bits.
4. **Modular arithmetic.** Every modular step of both paths (a(s) mod p by Horner, c = a(s) r^-1, c_k = c -+ k s^5)
   uses one Barrett routine valid for every p < 2^63 (mu = floor(2^2k/p) < 2^(k+1); two conditional subtractions);
   `slice12 barrett` compares it with the C remainder on 1,220,000 reductions at every bit length 3-63 (0 mismatches).
   There is no `%` in any sieving path.
5. **Scan.** Per row and per threshold block (4.4): a cell is a *preliminary survivor* iff f-sum >= T_f and
   g-sum >= T_g; survivors are listed and their f-sum set to the sentinel 0xFFFF (unreachable by a true sum).
6. **Recovery.** One pass over each segment's retained records copies every record whose cell holds the sentinel to
   the survivor list: a survivor obtains all its sieved ideals >= 2^20 with their global ids.

**Correctness.** `slice12 check` runs both paths (several primes, both sides), apply, scan and recovery on a slice
and compares every cell of every segment (digests of the complete arrays) and the recovered record sets with a
brute-force reference that tests each point of the box against Q_b and the coset congruences directly:
`measurements12.txt` shows `MATCH` on slices of 2^18-2^24 cells with 1 to 10 segments. The check also caught a defect
of the v10-style setup (it assumed |a_i| < p), which the Barrett path fixes. Primes up to 2^62 are exercised by the
`recentre` and `barrett` modes (membership of every recentred origin in its coset: 0 failures). The kernel aborts on
p >= 2^63 by design; X2 = 2^63 keeps every sieved prime in range.

### 4.4 Thresholds: certified per-block norm bounds (`nbound12.c`)

For a row (a and b1..b5 fixed) and a block b0 in [lo, hi], `nbound12` bounds 16 log2 N_f and 16 log2 N_g from the four
places v of K_h (two real, two complex: rational Rouche disks of radius 2^-100 around the roots of h, each certified to
contain exactly one root), evaluating N = prod_v |F_v| (F_v the x-homogenised f at (a(v), b(v))) and G_v = a(v)^2 +
|u| v b(v)^2 in outward binary64 interval arithmetic (every operation widened by one ulp on each side; no libm log:
a 17-entry dyadic table of 2^(j/16) bounds). A block is **accepted** when both widths are <= W = 128 (8 bits);
otherwise it is halved down to single points, and an unresolved single point is **rejected** (dropped from the sieve,
never a survivor; its yield weight is charged: 0 in the pool). Production blocks are the aligned 64-cell blocks of
each row (<= 127 evaluations per block); an accepted leaf gets the thresholds T = lower bound - 16 x slack, a
rejected one T_f = 0xFFFF.

`nbound12` is an operation-for-operation port of `embedding_bounds.py` (the Python reference, in the package); its
`check` mode replays the reference's leaves and bounds on all 2,401 pool rows (16-blocks: 4,221 evaluations, 3,311
leaves, 0 mismatches), and [2] recomputes the 64/128 statistics of the regenerated rows with the reference and checks
that the leaf bounds contain the exact norms. Over the pool (slice-weighted) the 64/128 blocks cost **0.04232
evaluations and 0.02932 leaves per cell**, with no rejected point and no non-finite value. Per row, the embedded
statistics (evaluations, leaves, block length, row length) are charged individually.

**Two deterministic implications** (used in 4.6-4.7). Let S be the true log2 of the product of the distinct sieved
primes dividing N (valuations 1 at p >= 256 in the retained subevent), S_16 the kernel's sum of floored scaled logs,
L <= 16 log2 N <= L + 128 the accepted bounds, D = log2 N - S the omitted mass (p < 2^8 with their valuations, primes
>= X2, unsieved (K_h-degree >= 2) and repeated factors). (i) *Retention*: with c the number of distinct sieved primes
(c <= floor(ceil(log2 N)/8)), S_16 >= 16 S - c, so D <= slack - (c + 1)/16 implies S_16 >= L - 16 slack: the cell
passes. (ii) *Survivors*: a passing cell has 16 S >= S_16 >= L - 16 slack >= 16 log2 N - 128 - 16 slack, i.e.
D <= slack + 8 bits, where every sieved valuation may be counted in full (a super-event of passage).

### 4.5 Cofactorisation of survivors

A preliminary survivor costs ([5]): exact norms N_f, N_g (two resultants over O_h, <= 1302 and 520 bits: 2^22
instructions); trial division of both by the 82,025 primes below 2^20 with valuations (<= 21 limbs x ~5 instructions
per limb per prime, both sides: 2^24.5), the recorded ideals' primes divided out with their valuations and each
p-part resolved into factor-base ideals (roots of a(s) - b(s) x in O_h/P for the K_h-primes P above p; the index primes
through the maximal order); what remains per side is a cofactor
of at most slack + 8 bits (4.4 (ii)), factored completely by ECM with the full curve budget for primes up to B = 2^69
(GMP-ECM table: B1 = 50,000, 214 curves for 25 digits, x ln 100 for 99% per-prime success; stage 1 = 1.44 B1 bits of
a Montgomery ladder at 11 modular multiplications per bit, stage 2 charged as a second stage 1; one 5-limb Montgomery
multiplication = **367.76 instructions**, the measured `ecm12.c` routine, checked on 20,000 products against an exact
reference): 2^39.06 instructions per side, 2^40.1 per survivor. The survivor is a relation iff every prime found is
<= B (every prime ideal is then in the factor base) and the pair is primitive (HNF test, 2^34 instructions, on relations
only). The 1% per-prime ECM miss is charged on the yield as 0.99^(lambda_f + lambda_g), lambda the expected number of
prime factors in (X2, B] of a B-smooth integer of the row's size (Dickman-conditioned; E[0.99^K] >= 0.99^E[K]).

### 4.6 Yield: retained relation orbits

The pool (`pool12-1.b85` + `pool12-2.b85`) is 20,000 uniform eligible pairs: a uniform in the eligible a-set, then b
uniform in the half-slice of a (per-index seeded, `evidence12.regenerate`; 2,000 with seed 1, 18,000 with seed 2); a
pair represents va / 20,000 eligible a and the weight w = (1 - Q_a/A^2)^3 of its slice. The 2,401 canonical pairs
(exhaustive oracle) carry Q_a (1/2^16), the exact norm sizes (1/16 bit), primitivity and the norm-bound block
statistics (evaluations, leaves, block length, row length rounded down to 16); [2] regenerates seed-1 indices from
scratch and verifies all of these. The retained relation orbits are

  rels = (va / 20,000) x vb x (2401/2500) x sum over canonical primitive non-rejected rows of
         w x rho((log2 N_f + alpha_f)/log2 B) x rho((log2 N_g + alpha_g)/log2 B) x r_f x r_g x 0.99^(lambda_f + lambda_g),

where rho is Dickman's function ([1]) and r_f, r_g are the **retentions** of 4.4 (i): the probability, for a B-smooth
integer of the row's size, that D <= slack - (c + 1)/16 with all valuations at p >= 256 equal to 1 (a conservative
sub-event: non-squarefree smooth integers are not credited). They are computed by `threshold_model.Retention` as
conditional probabilities of a finite Euler-product model: exact primes below 2^12, prime-density quadrature above
(step 1/4 bit, exponential tilt 0.05 for range, powers to the 8th), the generating-function recurrence
n p_n = sum j lambda_j p_{n-j} with signed coefficients for the squarefree selection, the small primes and the primes
>= X2 forming the omitted mass. Support: (a) on exact uniform smooth integers (tilted exact factor lists, exact
threshold decisions; `threshold_model.exact_smooth_samples`) the model's sufficient event matches the counts:
B = 2^10, window [2^100, 2^101), X = 2^9, slack 40: model 0.1206 vs 241/2000 = 0.1205 (actual threshold and squarefree
272/2000); B = 2^12, [2^120, 2^121): 325/1000 sufficient, 356/1000 actual; B = 2^14, [2^154, 2^155): 96/300, 108/300;
(b) grid convergence: at norm 770 / slack 296 the retention is 0.79028, 0.79020, 0.79033 at steps 1/4, 1/8, 1/16;
(c) conditional ratios are tilt-invariant on a fixed grid. At the operating point r_f = 0.785 and r_g = 0.895 at the
mean norms ([4]).

Result ([4]): **rels = 2^64.219** against the unknowns U = (NBAR/2) Li(B) + 2^31 = **2^63.451** sigma-orbits of the
factor base, i.e. excess x U = 2^63.830 at excess 1.3: a margin of **0.389 bits** (3.3 pool standard errors of 0.118
bits; the package requires >= 0.30).

### 4.7 Cost per cell (measured constants, per pool row)

All constants are executed-instruction counts (callgrind Ir, GCC 12 -O3 -fno-tree-vectorize, linux/arm64) of the
shipped kernels by the track's differential method (`KERNELS.md`). For the row of canonical a with slice S cells and
K lines, the model integrates over the primes in 1/4-bit bins (sum of 1/p over a bin = ln of the ratio of its log2
endpoints), with NBAR = 2 ideals per prime:

| term | charged | per cell (mean) |
|---|---|---|
| layout (allocation, metadata, rows, first zero pass) | 12.325 per cell (whole-program figure at a geometry with far more metadata and rows per cell than production) | 12.3 |
| zero pass | 0.625 per cell | 0.6 |
| line path, p in [2^8, 2^20): NBAR/p updates per cell | per update incl. its per-(pair, line, coset) setup, by h = S/(K p) points per line coset: 24.10 (h = 2^13.8), 31.11 (2^9.8), 88.31 (2^5.8), 310.70 (2^1.8); production h >= 2^34 is charged the densest measured value | 44.2 |
| record path, p in [2^20, 2^63): NBAR/p hits per cell | per hit incl. walk, fold, store and per-coset recentring, by h = S/p points per slice coset: 44.90 (2^18.9), 78.21 (2^10.9), 266.97 (2^3.9), 938.08 (2^0.9), 2231.66 (2^-1.1), linear in log2 h between; 1600 per (slice, pair) below h = 2^-1.1 (not reached: the smallest slice, 2^64.5 cells, has h >= 2^1.5 at X2) | 153.1 |
| apply + recovery passes | 7.00 + 7.003 per record (2.29 records per cell) | 32.2 |
| per-(slice, pair) setup (a(s) mod p, c, 36 streamed words) | 500 x 2^57.59 pairs / S | 1.9 |
| scan | 8.0 per cell + 21 per threshold leaf (from 9.313 at 16-cell and 13.250 at 4-cell blocks) | 8.6 |
| norm-bound blocks | 5,900 per evaluation (measured 5,862.7) x the row's evaluations per cell + 4,200 per row (measured 4,120) + 15,300 per slice | 252.5 |
| cofactorisation | survivors per cell (4.4 (ii): generic model with every sieved valuation counted, slack + 8 bits, B = 850 bits so that every pool norm is covered) x 2^40.1 | 11.1 |
| **total** | | **516.5** |

The population is **2^127.389 cells** (va / 20,000 x the sum of the upper-count slice sizes of the canonical rows;
mean slice 2^65.60, largest 2^67.53, 185 segments). Relation collection: 2^127.389 x 516.5 / 20,000 = **2^122.114
units** ([5]).

### 4.8 Operating point

B = 2^69, X2 = 2^63, slacks (292, 220) bits, excess 1.3, acceptance width 8 bits, 64-cell blocks, alpha (5, 1),
A^2 = 2^21, theta = 1/2. The slack grid 284-300 x 212-228 was searched offline with the same module; the neighbours
printed by [5] move the total by < 0.02 bits. Lower X2 is worse (the survivor rate explodes); B = 2^70 would push the
per-side ideal count past 2^64 ids.

## 5. Filtering and linear algebra

**Unknowns.** The sigma-orbits of the factor base (4), U = (NBAR/2) Li(B) + 2^31 = 2^63.451 (indices fit 64 bits),
plus k <= 82 Schirokauer-map columns (unit ranks <= 71 for K_f, <= 11 for K_g; all k are kept as dense
right-hand sides and the k dropped rows give a k x k dense system, so dependent map columns cost nothing but
verification rows).

**Filtering (H3).** Purge (singleton removal, clique removal to the target excess) and Markowitz merge with unit pivots
and coefficients bounded by 255, as in CADO-NFS, **simulated** on a random relation hypergraph with the production
occurrence profile (`filtsim.py`: a column of norm q = B^t on side s appears in a relation orbit with probability
2 kappa_s(t)/q, kappa_s(t) = rho(u_s - t)/rho(u_s); the small-prime columns are explicit with valuation coefficients
and the dense mid band is propagated through the merge). The simulation makes the small-prime columns explicit
(valuation coefficients, |coefficient| <= 255 enforced, max 29 observed) and carries the mid-band columns (norms
2^12.4-2^55) as a propagated per-row weight; a coefficient of a mid-band column in a merged row is bounded by the row's
combination multiplicity, and production rejects from the merge any row whose coefficient would exceed 255 (CADO's
rule), which costs excess, not correctness. At this profile (u_f = 11.12, u_g = 7.29) and excess 1.3, two seeds
at scale 2^17 give N/U = 0.0538 and 0.0519 at full row weight 100 (max |coefficient| 26); [6] runs scale 2^13 in the sandbox and embeds the
table of all excesses. The model charges **N = 2 x 0.0538 x U = 2^60.23** (a factor-2 pad on the dimension, x4
on the linear algebra) and row weight 100.

**Square system and cost.** Block Wiedemann on the square sparse part with the k map columns as right-hand sides
(blocking factor >= k, no extra products): 3 N single-vector products of N x 100 entries at **36 instructions per
entry** (`la_mac.s`: 33 per entry with 64-bit column indices, [9] re-counts the branch-free loop; plus the per-row
Barrett finish and the combination of the four column-block accumulators amortised over the row weight; the vector
of 2^60.23 elements of 32 bytes is held in four blocks of <= 2^63 bytes): **2^119.579 units**. The block
Berlekamp-Massey step is O(k^2 N log^2 N) << N^2 x 100. Rank: a random sparse square matrix over F_r with every column
of weight >= 2 (purge guarantees it) and small non-degenerate coefficients is singular with probability O(1/r) beyond
the structural dependencies merge removes as 0 = 0 rows; the dropped excess rows verify the solution and a failure is
repaired by collecting the 10% more relations that the excess 1.3 already contains (the yield margin is not used).

**Schirokauer maps.** <= 2^37 instructions per relation (exponent <= 72 x 255 bits in O_f / r^2) for 2^64.3 relations:
<= 2^88 units.

## 6. No large primes; nothing pairing-specific

Large primes above B would put columns of expected weight < 1 into a matrix already governed by the top of the factor
base; the slack tier of 4.3-4.5 uses the primes between X2 and B more cheaply. Nothing pairing-specific is used:
final-exponentiation preimages (Kim-Cheon 2012/657) leave a fibre of ~2^4314 elements; Miller inversion (Satoh
2019/385) needs the unreduced Miller value; hidden-root attacks (Vercauteren 2008/261) have degree >= 2^359;
Cheon-type restriction needs auxiliary powers.

## 7. Individual logarithms: boot and descent ([7])

Two logarithms (z and g); their ratio is a. The target is represented in F_p[y, x]/(h, g).

**Boot (g side).** For a randomiser t, T' = T gamma^t; the lattice L_T' = {(a, b) : a - b theta in T' F_{p^6}^*} has
index p^6 (six F_p-linear conditions on 12 coefficients; v10 verified an embedded example exactly); every nonzero
point lifts T' times a subfield element (r-log 0). After LLL (2^30 per t) the points of a ball of 2^37.6 candidates
have g-side norms of ~2^2725 bits (measured in v10 over random targets, 2661-2723 by radius; the model uses 2725); a
candidate is tested by exact norms and a full ECM run to B1 = 2^260 (GMP-ECM table extrapolated beyond 80 digits as
stated in `descent12.py`; the ECM multiplication is the measured 5-limb figure scaled by (L/5)^2).

**Tree.** A node is a prime ideal of norm 2^l on side s (2-bit bins, representative at the bin's midpoint above B);
its special-q lattice (12-dim) is reduced once (2^30). The node targets 2^n (n at or below the bin's lower edge, so
every child is in a strictly lower bin whose subtree cost is already known and is added in full with a x1.5 pad on the
child count): among 2^S lattice points with coordinates 2^(l/12 + S/12) one with both norms 2^n-smooth is needed (S the
least with 2^S P >= 4, a x2 node pad for failure); candidates have norms 72 c + 54 (f; measured for box-uniform coefficients of size 2^c) and 12 c + 6 log2 |u| + 2 (g)
bits with the special-q divided out; they are tested by exact norms + full ECM to 2^n on both sides, or, when n <= 80,
by sieving with the ideals <= 2^65 in the special-q lattice (per-(node, ideal) reduction 2^13, 250 instructions per
point, slack k n per side; the survivor rate is a Dickman estimate here, covered by the x2 node pad) and ECM. Children are
descended on their own side; the tree ends at B exactly; a terminal ideal with no determined logarithm (expected
fraction < 2^-10) is re-descended from level B + 1 inside the node pad. A child bin that is not costed makes the
programme fail (no silent zero).

Result: per logarithm **2^98.82 units**, both **2^99.82**. The attack runs at most **2^106 units** of descent per
challenge, else aborts: by Markov, P(abort) <= 2^(99.82 - 106) = 0.014.

## 8. Success probability

Relation collection and linear algebra do not involve the challenge; given the relation count (H1-H3 with the margin)
they succeed up to the rank failure of section 5. The descent succeeds within its meter
with probability >= 0.986; z = 1 is free. Overall success >= 0.95 >= **0.39** over uniform (Q, z) and the algorithm's
coins (the challenge is exactly uniform over GT). The 0.39 requirement tolerates a descent expectation 5.8 bits above
the modelled one, and the yield margin of 0.389 bits (3.3 pool standard errors) is the package's own acceptance
condition, not a success-probability reserve: a yield shortfall of the whole margin would leave excess 0.99 and no
square system, and the claim would simply be false rather than succeed with lower probability.

## 9. Total work

| component | log2 units |
|---|---|
| relation collection (2^127.389 cells x 516.5 instructions) | 122.114 |
| linear algebra (block Wiedemann, N = 2^60.23, weight 100, 36 instr/entry) | 119.579 |
| descent, two logarithms, metered cap | <= 106 |
| canonical a-list (box scan, HNF keys, segmented external sort) | <= 90 |
| Schirokauer maps | <= 88 |
| factor base (2^64.5 ideals x (roots 2^20 + reduction 2^13 + tables)) | <= 72 |
| filtering passes, relation/matrix storage, challenge pairing, randomness | <= 64 |
| verification (one pairing, subgroup checks, [a]P0, serialization) | <= 11 |
| **total** | **122.343555** |

The submitted bound **122.343600** is the full-ledger logarithm rounded up; `evidence12.py` prints
`total_log2W_full_ledger = 122.343555`.

## 10. Sensitivity (printed by [5])

- all sieve constants x1.5: total 122.847; matrix pad 4 instead of 2: total 122.871;
- slacks (288, 216): relation collection 122.104, margin 0.325; (296, 224): 122.128, margin 0.422;
- alpha +1 bit on either side changes the yield by < 0.05 bits (the log-derivative of rho at u = 11 is -0.048 per bit);
- NBAR 1.947 (v10's value) instead of 2.0: total -0.04 bits;
- without the unit quotient (crediting every point of the ellipsoid) the same model gives ~3.5 bits less: the v1-v4 error.

## 11. Heuristics (all disclosed, with support and limitations)

- **H1 (smoothness).** Norms are exact resultants ([2] recomputes them from seeds). B-smoothness probability rho(u)
  with GS21 alpha rounded up; a model, not an execution.
- **H1' (threshold retention).** The finite Euler-product model of 4.6 transfers generic-integer threshold statistics
  to the NFS norms (the alpha correction shifts the size; small-prime valuations of norms are governed by alpha rather
  than by the generic model); its calibrations are on exact generic smooth integers at B = 2^10-2^14, not on NFS norms
  at B = 2^69. The squarefree restriction is a conservative sub-event, the 0.30-bit margin and the ECM-miss factor are
  extra conservatism (the ECM-miss factor is 0.99 per expected prime in (X2, B], 0.972 at the mean norms), and the generic super-event of 4.4 (ii) is a model of the survivor rate, not a deterministic
  bound on it (the cofactorisation term is 2% of the cost).
- **H2 (region and pool).** 20,000 exactly sampled pairs, canonical flags by an exhaustive oracle, lattice counts
  bounded rationally; the only extrapolation is from pool statistics (2^11.2 canonical rows) to the region, with
  per-row weights and the standard error reported.
- **H3 (filtering).** Random hypergraph with the production occurrence profile: the standard heuristic behind every
  NFS filtering prediction, simulated at scale 2^17 with the dense columns propagated, padded x2 on the dimension.
- **H4 (kernel constants).** Every sieve constant is an executed-instruction count of a shipped kernel checked against
  brute force on every cell; production slices (2^65.6 cells) exceed the benchmark ones (<= 2^26), which the model
  handles by indexing each cost on points per coset (the quantity fixing the walk's regime), charging the densest
  measured value where production is denser still, and charging the layout at a geometry with more metadata per cell
  than production. Loads and stores cost one instruction each, as the metric prescribes. No full-size run exists.
- **H5 (descent).** Explicit dynamic programme with the measured ECM multiplication, table curve counts (99%
  per-factor), x2 node and x1.5 child pads; the boot norm is v10's measurement.
- **H6 (pairing transfer).** e(P0, Q) generates GT for nonidentity Q; standard.
- **H7 (field).** O_h = Z[y] is proved ([0c]); class-number or regulator statements are not used anywhere.

Limitations: no NFS-DLP at this field size has been executed; the claim is an analytic bound grounded in exact-norm
sampling, measured kernels, a calibrated retention model and simulated filtering, all reproducible from the package.

## 12. Memory (reported, not scored)

One a-slice job (largest slice, 2^67.53 cells): cells 4 bytes each = 2^69.53 bytes in 185 segments of <= 2^62 bytes;
suffix metadata 8 bytes x 2^48.33 entries x 2047 lines = 2^62.33 bytes; retained records 16 bytes x 2.29 per cell =
2^72.73 bytes; total **2^72.88 bytes per job**. Jobs are independent: peak = (concurrent jobs) x 2^72.88 bytes plus the
read-only streamed pair data (2^57.59 x 288 bytes = 2^65.8). Canonical-list staging 2^78 bytes; relation store
2^64.3 x ~48 = 2^69.9 bytes; filtered matrix 2^60.23 x 100 x 5 = 2^68.9 bytes; block-Wiedemann vectors < 2^69 bytes. No
64-bit index spans more than one segment.

## Attribution

Research, kernels, validation and this package: Claude Fable 5.1 (harness Claude Code), v5-v12; the v11 workbench
modules (`ideal_canonical.py`, `embedding_bounds.py`, `threshold_model.py`, `lattice_count_bounds.py`, the sampler) were
written with Codex in the same project; the exact-norm machinery and the boot lattice descend from v1-v4 (Kimi K3 /
Claude Opus 5.5, Claude Code). Polynomial selection and alpha: Guillevic-Singh (Math. Cryptology 2021); filtering
methodology: CADO-NFS, De Micheli-Gaudry-Pierrot (ASIACRYPT 2021); FAPI-2 context: Galbraith-Hess-Vercauteren
(2007/256), Kim-Cheon (2012/657), Satoh (2019/385), Vercauteren (2008/261). No credentials or private data.
