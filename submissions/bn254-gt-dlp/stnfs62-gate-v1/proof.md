# BN254 GT discrete logarithm by a twisted (6,2) special tower number field sieve

Track `bn254-gt-dlp`. Unit: C = 20,000 scalar AArch64 instruction equivalents. This is an analytic,
exploratory claim (`plausible_not_refuted`): no production attack was run, and success and work rest on the
heuristics H1-H7 of section 10, each stated with its support and its limits. Nothing is amortized over other
targets, no advice is used, and no earlier verdict is a premise. Every number used by the score is recomputed
by `tnfs.py` from the production-size samples stored in the four `data_*.json` files; `experiment.py` reruns that computation
and the exact checks cited as [E1]-[E6].

## Claim

W is the total classical work of the one attack on the one target, in units. Step 0 is a run gate: draw five
fresh uniform bits; unless their value is below 13, stop as a failure (below one unit). Past the gate the
algorithm below returns x with g^x = z with probability at least 0.96 under H1-H7 (section 7). The gate is independent of the
challenge and of all other coins, so success is at least (13/32) x 0.96 = 0.39 over uniform challenges and
fresh coins, and the expected total work is

    E[W] <= 1 + (13/32) x 2^100.5759 < 2^99.28 units (C = 20,000).

The score 99.28 is an expectation over the challenge and the coins of the one attack, not a bound on every
run. A run past the gate has expected work at most 2^100.5759 units (section 8); the only rows that are
expectations are named there, and the individual logarithms are additionally metered at 2^116.24
instructions (section 6).

## 1. Algorithm overview

Input z in GT (order r). Output: x in [0, r) with g^x = z, g = e(P0, Q0) the protected generator.

0. Run gate (Claim): draw five bits; if their value is at least 13, stop (failure).
1. Decode z; check z != 0 and z^r = 1. If z = 1 return x = 0.
2. Field map (section 2): move g and z into the representation F_p[y,x]/(h, x^2 - u y).
3. Main phase, independent of z: factor base (section 3), relation collection (section 4), linear algebra
   (section 5). It yields virtual logarithms of the factor-base ideals that occurred.
4. Individual logarithms of g and of z (section 6), giving vlog(g) and vlog(z); x = vlog(z)/vlog(g) mod r.
5. Check g^x = z in the protected representation and output x as 64 hex characters. Any wrong value is caught by this check and counted as a failure.

## 2. Tower, polynomials, field map ([E1])

u = 0x44e992b44a6909f1, p = 36u^4 + 36u^3 + 24u^2 + 6u + 1, r = 36u^4 + 36u^3 + 18u^2 + 6u + 1 (the protected
constants). K_h = Q[y]/h with

    h = y^6 + y^4 - 2y - 1.

Facts about K_h (computed with PARI/GP 2.19; the ones marked [E1] are rechecked by the experiment):
disc(h) = 340416 = 2^6 x 3^3 x 197 equals the field discriminant, so O_h = Z[y]; signature (2, 2), unit rank
3; class number 1 and fundamental units y, e2 = y^5 + y^3 + y^2 - 1, e3 = y^5 + y^3 + y - 1 (certified by
`bnfcertify`; [E1] rechecks the inverses y^5 + y^3 - 2, y^3 - 1, y^4 - y^3 + y^2 - 1 and N(y) = -1); the only
roots of unity are +-1. h is irreducible modulo p ([E1], Rabin's test), so F_p[y]/h = F_(p^6).

Side polynomials over O_h, with P(X) = 36X^4 + 36X^3 + 24X^2 + 6X + 1 the BN polynomial (P(u) = p):

    f(x) = y^4 P(x^2 / y) = 36x^8 + 36y x^6 + 24y^2 x^4 + 6y^3 x^2 + y^4,      g(x) = x^2 - u y.

At a root of g, f = y^4 P(u) = y^4 p = 0 modulo p ([E1]), so f and g share the root x modulo (p, h). Twist:
x^2 - u is reducible over F_(p^6) (u is a square in F_(p^2)), which is why g carries the unit y. u y is a
non-square in F_(p^6): (u y)^((p^6-1)/2) = -1 ([E1]; it follows from N(y) = -1 and p = 3 mod 4). Hence
F_p[y,x]/(h, g) = F_(p^12). K_f = K_h[x]/f has absolute degree 48 and K_g = K_h[x]/g absolute degree 12; the
absolute polynomials Res_y(h, f), Res_y(h, g) are irreducible over Q (PARI `polisirreducible`; not rechecked
by the experiment). K_f is totally complex (P has no real root, so f has none in any real embedding of K_h):
unit rank 23. K_g has two real embeddings (u y > 0 only at the real root y = 1.13..): unit rank 6. The method
uses 29 = 23 + 6 Schirokauer maps.

r divides p^4 - p^2 + 1 and gcd(r, p^6 - 1) = 1 ([E1]): every element of F_(p^6)* has trivial component in the
order-r subgroup, so factors from F_(p^6) do not change logarithms modulo r.

Field map. The protected field is F_p[w]/(w^12 - 18w^6 + 82). The wrapper finds a root of h in it and then a
square root of u y (equal-degree splitting with random elements, each trial splitting with probability > 1/4),
which defines the isomorphism; the 12 x 12 coordinate matrix is inverted by Gaussian elimination. One trial
costs fewer than 2^17 F_(p^12) multiplications of fewer than 2^14 instructions each; trials are capped at 2^10
(failure < 2^-400). With decoding, the group check z^r = 1, the final exponentiation or pairing, the scalar
multiplication and output, the wrapper stays far below the 2^70 instructions charged for it.

## 3. Factor base

B = 2^54. A degree-one prime ideal of K_f (resp. K_g) is (q, y - s, x - t) with q prime, s a simple root of h
mod q and t a simple root of f(s, x) (resp. g(s, x)) mod q. The factor base holds all such ideals with q <= B,
plus the finitely many exceptional ideals met in kept relations (ideals above 2, 3, 197, non-simple roots,
ideals of residue degree >= 2 of norm < 2^47; fewer than 2^25 of them, since such an ideal lies above a prime
below 2^23.5). By the prime ideal theorem each side has li(B)(1 + o(1)) prime ideals of norm <= B, and
li(B) = 1.028 B/ln B at B = 2^54. The ledger charges

    N = 2.2 B / ln B = 2^49.91 columns      (a 7% allowance over 2 li(B); H3),

and the linear system has dimension D = N + 29 + 64. Building the base (root finding modulo each prime, at most
2^24 instructions per prime, and a sieve for the primes) costs below 2^72.8 instructions.

## 4. Relation collection

Cells. (a, b) in Z[y]^2 with all twelve coefficients in [-A, A], A = 232, one of each pair +-(a, b):
V = (2A + 1)^12 / 2 = 2^105.333 cells. A cell is a relation when both integers

    N_f(a, b) = |N_(K_h/Q)(sum_i f_i a^i b^(8-i))|,     N_g(a, b) = |N_(K_h/Q)(a^2 - u y b^2)|

are B-smooth. `tnfs.norm_f`, `tnfs.norm_g` compute them exactly; [E2] replays 40 stored cells and [E3] draws
fresh uniform cells. Sizes (2,500 uniform draws, 1,981 primitive): mean log2 N_f = 442.2, log2 N_g = 471.4.

Sieve. For each degree-one factor-base ideal (q, s, t) of either side, the cells it divides are those with
a(s) = t b(s) mod q. They are enumerated by matching: the residues a(s) mod q of all (2A+1)^6 vectors a are
generated in Gray-code order (one modular addition each) and stored in a table keyed by residue; for each
of the (2A+1)^6 / 2 vectors b the residue t b(s) is looked up and each matching a gives one hit. Every hit adds
a 16-bit fixed-point log2 q (unit 1/64 bit) to that side's counter of the cell. No lattice reduction and no
assumption on the shape of the hit set is used; the number of hits of an ideal is the number of matching
pairs, (2A+1)^12 / (2q) on average.

Storage and addressing. Counters (two 16-bit words per cell) are split into chunks of 2^32 cells, each with
its own owner and 64-bit local addresses; a hit is a record (two-word cell index, value) appended to its
chunk's queue and applied by the owner. The per-ideal table has 2^53.2 entries and fits one address space.

Price per cell (`tnfs.model`, row `instr_per_cell` = 182.9):

- hits: 32 instructions per hit (match loop 6, two-word index 4, queue append 8, owner load/add/store and loop
  8, slack 6). Hits per cell are sum 1/q over the ideals: the exact sums over all degree-one simple-root ideals
  with q < 4096 are 0.7533 (f) and 1.1870 (g) (PARI enumeration, stored in `data_main.json`), and above 4096 each
  side is charged 1.05 x ln(54/12) (Mertens with one ideal per prime on average, 5% allowance): 5.099 hits.
- tables: 32 instructions per stored or looked-up residue, (1 + 1/2)(2A+1)^6 residues per ideal:
  96 N / (2A+1)^6 = 10.1 per cell.
- initialisation and threshold scan: 8 per cell.
- stage 2 (below): 1.7 per cell.

Thresholds. Slack s = 47 bits. Stage 1 keeps a cell when its f counter is at least n_f,min - 48 and its g
counter at least n_g,min - 48, where n_f,min = 394.6 and n_g,min = 450.0 are the 0.5% quantiles of the sampled
norm sizes (cells with smaller norms are not counted in the yield). Stage 2 computes both exact norms of a
kept cell (2^16 instructions each) and keeps it when each counter is at least log2 N - 47 (one bit absorbs
fixed-point rounding over at most 64 hits). Lemma: each hit ideal contributes at least its log2 q to log2 N
(distinct ideals above q each divide the norm), so a counter never exceeds the log of the B-smooth part of N.
A stage-2 survivor therefore has a cofactor below 2^47 < B on each side: it is B-smooth on both sides. There
are no false survivors. Stage-1 pass rate (work only): an integer n <= X has a B-smooth divisor >= Y for at
most X ln(X/Y) rho(ln Y / ln B) values of n; `tnfs.model` evaluates this for each side and charges the smaller
one times 2^17.

Yield ([E6]; `tnfs.model`). For each sampled primitive cell with exact norm sizes n_f, n_g, exact unsieved
masses X_f, X_g and unit multiplicity m,

    y = rho((n_f ln 2 + alpha_f) / ln B) rho((n_g ln 2 + alpha_g) / ln B) / m    if X_f, X_g <= 47 and
        n_f >= n_f,min, n_g >= n_g,min, else 0;   non-primitive draws count 0.

- rho is Dickman's function (table by the identity u rho(u) = int_(u-1)^u rho, step 1/512; [E6] prints
  rho(5) = 3.547e-4, rho(8) = 3.232e-8, rho(10) = 2.770e-11).
- alpha (Murphy): alpha = sum_(q < 1000) ln q / (q - 1) - E[ln of the 1000-smooth part of N], measured on the
  sampled exact norms and charged at + 3 SE: alpha_f = -0.60 (charged 0: no credit), alpha_g = +1.09.
- X is the exact mass of N not credited by the sieve among primes below 4096 (all powers of 2, 3, 197, excess
  powers, ideals of higher degree and non-simple roots); in all 1,981 primitive draws X_f <= 44.7 and X_g <= 23.1 bits. Above
  4096 an uncredited factor needs q^2 | N; a factor 0.98 is charged for it. The masses are measured on
  uniform cells; relation cells carry somewhat more small-prime mass, and the largest value seen is 2.3 bits
  under the slack, so a small relation-conditioned loss is possible; it is covered by the 47% relation surplus
  below, not by a separate factor.
- m is the number of units eps = y^i e2^j e3^k (|i|, |j|, |k| <= 4, up to sign) with eps a and eps b in the
  box: these cells give the same principal ideal, and 1/m counts it once. Mean 1/m = 0.790. [E2] replays m.
- primitive: the ideal (a, b) is O_h (PARI `idealadd`); fraction 0.792.

Mean y = 2^-54.724 per cell, relative SE 2.6%. Charged distinct relations: V x (mean - 3 SE) x 0.98 =
2^50.463 = 1.466 D. The algorithm needs D of them (section 5); the 47% surplus is the margin of H1/H2.

Kept relations. Each survivor is divided by every factor-base ideal (one residue evaluation, 64 instructions
per ideal), its exceptional part (below 2^47 per side) is factored and valued (2^30 instructions: trial division by
the 2^20 primes below 2^23.5 at 64 instructions, then at most 47 valuations at each of at most 48 ideals
above a prime from two-element representations, below 2^16 each), and
duplicates are removed by sorting (2^12): at most V x (raw mean + 3 SE) x (64 N + 2^30 + 2^12) = 2^107.78
instructions. Schirokauer maps of the D relations used (29 maps; one exponentiation in (Z/r^2)[y,x]/(h, f) with an
exponent below r^48, 12,192 bits at 1.5 products per bit, a product being 48^2 coefficient products of 8 limbs
at 1,088 instructions plus reduction, below 2^22: 2^36.2; the g side is below 2^30; charged 2^38 per
relation): 2^87.9.

Collection: V x 182.94 = 2^112.848 instructions.

## 5. Linear algebra

Relations are taken in order until their number equals the number of ideal columns that occur in them plus
29 Schirokauer columns plus 64; this happens within the first D = N + 29 + 64 distinct relations because at
most N ideals exist. The matrix over Z/r has these rows, the occurring ideal columns, the 29 Schirokauer
columns and 64 zero pad columns: square, of dimension at most D. H3: its kernel restricted to the ideal and
Schirokauer columns is the line of virtual logarithms. A kernel vector is found by Wiedemann's algorithm: 2D
products for the scalar sequence, Berlekamp-Massey, D products for the solution; the pad coordinates of the
solution are discarded. Ideals that do not occur in these rows have no logarithm (H7a). The solution is checked on the unused relations; a
failed check is a failure (probability < 2^-200 over the random projections, given H3).

Price per row and product: 96 sparse entries at 24 instructions (index load, four-limb add or subtract with
carry into the accumulator; a norm below 2^516 has at most 43 prime ideals of norm >= 2^12, the largest sampled
norm is 2^488.7, and sieved ideals below 2^12 average 1.94 per cell: an average weight), 29 dense four-limb products at 128 plus one reduction 256, row overhead
256, the dot product 130 per row: below 7,200. Three passes plus Berlekamp-Massey (below 800 D^2):

    22,400 D^2 = 2^114.274 instructions.

Storage: matrix 2^60.6 bytes, vectors 2^55 bytes.

## 6. Individual logarithms (expected work)

Both g and z need one; g being public does not make it free. Every search below is a sequence of independent,
identically distributed trials made with the algorithm's own fresh coins, so the number of trials is geometric
with mean exactly 1/pi, pi the success probability of one trial. This removes any assumption on how relations
arrive inside one lattice or on independence between its cells. It does not remove the smoothness heuristic:
for a node, pi is the realized fraction of relation cells in one fixed ball, and the ledger charges the
model value of that fraction. That the realized fraction of every encountered ideal is not below the charged
one is premise H7(c), stated with its support in section 10. If it fails for some ideal (in the extreme, a
ball without any relation), the search does not end and the meter below stops the run as a failure: the
effect is on success probability and work, never on correctness.

Boot (target T). One trial: draw e uniform in [1, r), T' = T^e = t0 + t1 x; LLL-reduce the lattice
{(a, b) in O_h^2 : a + (t0/t1) b = 0 in F_(p^6)} (dimension 12, determinant p^6, b weighted by 2^31); take its
first vector. Then a - b x = c T' with c in F_(p^6), which does not change the logarithm (section 2). The trial
succeeds when N_g(a, b) is 2^165-smooth, the ECM schedule finds all its primes, and every factor-base ideal in
it has a logarithm. 300 boot lattices were reduced: log2 N_g has mean 1708.6, maximum 1715.8 ([E5] rechecks
three stored vectors). Charged per trial: 2^44 (LLL: at most 1.1e5 swaps by the potential argument, each at most 576
operations on integers below 7,000 bits at 16 x 110^2 + 64 instructions: 2^43.5) + 2^24 (exponentiation) + the ECM schedule on a 1,724-bit modulus.
Success probability 2^-37.70 (`log2_boot_p`: mean of rho over the 300 norms with their measured alpha, - 2 SE,
times the factors of H6/H7).

Nodes. A prime ideal Q of norm q in (B, 2^165] met in a relation is a node. Bins by size, top edges
56, 60, 66, 74, 86, 104, 130, 165 bits. The lattice L_Q = {(a, b) : a(s) = t b(s) mod q} has determinant q;
it is LLL-reduced once (2^44 charged: entries below 2^165 give at most 6.2e4 swaps, 576 operations each on
integers below 3,960 bits, 2^41). One trial: draw c uniform in [-R, R]^12, form (a, b) = basis x c, compute
both exact norms, and test whether both (with q removed from Q's side) are smooth for the bound of the next
lower bin edge (B for the lowest bin). The relation then involves only smaller nodes and factor-base ideals.
Trials are independent draws with replacement, so trials until success are geometric with mean 1/p(Q), p(Q)
the fraction of relation cells in the ball; a repeated cell or a unit multiple is harmless. Smoothness test per
trial: for bounds up to 2^66 the product/remainder-tree batch test (H5; about 2^26 instructions per number,
`tnfs.batch_test_instr`; trials are tested in batches, and one extra batch is charged per node), above that
the ECM schedule (H4). The one relation found is then fully factored by ECM (charged).

Sampled lattices ([E4] verifies four stored bases: membership and determinant q). For every bin and side, six
random degree-one ideals at the bin's top edge were reduced and 60 uniform ball cells per ideal evaluated
exactly: 5,760 cells in `data_desc_f.json` and `data_desc_g.json`. For each ideal p_i is the mean of rho rho over its cells (non-primitive
cells count 0), and the bin is charged

    1/p = 1.5 x mean_i (1/p_i)

directly estimating E[1/p(Q)] (an estimate of 1/p_i from a finite sample is biased upward, the safe side),
then divided by the factors of H6/H7. Within a bin the six values of log2 p_i differ by at most 2.7 bits
(`ideal_spread_bits`; this includes the sampling noise of 60 cells). The radius R is admissible only if the worst sampled ideal has at least
64 expected relation cells in the ball. Table (`tnfs.model`, `nodes`): radius, log2 p, log2 expected work of a
node including its whole subtree:

| top bits | child bound | R | log2 p (f / g) | log2 work (f / g) |
|---:|---:|---:|---:|---:|
| 56 | 54 | 32 | -61.06 / -60.62 | 88.29 / 87.89 |
| 60 | 56 | 32 | -59.01 / -58.98 | 89.69 / 89.67 |
| 66 | 60 | 32 | -55.32 / -55.35 | 91.91 / 91.91 |
| 74 | 66 | 16 | -45.89 / -46.13 | 94.29 / 94.31 |
| 86 | 74 | 16 | -41.50 / -41.92 | 96.76 / 96.78 |
| 104 | 86 | 8 | -33.43 / -33.84 | 99.30 / 99.27 |
| 130 | 104 | 8 | -28.80 / -29.74 | 101.85 / 101.86 |
| 165 | 130 | 4 | -22.63 / -23.79 | 104.30 / 104.34 |

Bins are evaluated at their top edge: for the side of Q the norm with q removed, and for the other side the
norm itself, grow with the determinant q of the lattice, so the top edge is the least favourable size of the
bin (and only 2^(n-1) is removed for an n-bit ideal).

Recurrence. W(bin) = 2^44 + extra batch + (trial cost)/p + sum over lower bins of (expected number of prime
factors of the relation in that bin) x W(that bin, worse side). The expected number of prime factors of an
L-smooth norm of n bits in a window is the integral of kappa(t) dt/t, kappa(t) = rho(n/L - t)/rho(n/L), t the
prime's size relative to L (H6), evaluated at the largest sampled norm of the bin: between 2.4 and 12 children
per relation, 14.3 for the boot. By linearity of expectation and Wald's identity (a node's cost depends on the
past only through its ideal, and its bin's W bounds its conditional expectation), the expected work of one
logarithm is at most boot trials x boot trial cost + sum over children of W. Model value for both logarithms:
2^108.23 instructions. The ledger charges eight times that, 2^111.233, as the allowance for H5-H7; this row
is 8% of the total.

Meter. The two logarithms run under a meter set at 32 times the charged row (2^116.23 instructions); reaching
it is a failure. By Markov's inequality, if the expected work is at most the charged row, P(meter) <= 1/32.
Metered work is at most the unmetered work, so its expectation is at most the charged row.

## 7. Success

Past step 0 the run fails only if (i) the main phase falls short (fewer than D distinct relations, H3 fails,
Wiedemann fails): never under H1-H3, charged 0.005; (ii) the meter is reached: <= 1/32; (iii) the wrapper
fails or a logarithm of g is 0 mod r: < 2^-100. So P(success) >= 1 - 0.005 - 0.03125 - 2^-100 > 0.96 for every
z != 1 (the bound does not depend on z), and z = 1 is answered directly. With the gate, success is at least (13/32) x 0.96 = 0.39 exactly as required.

## 8. Ledger (`tnfs.model`, [E6])

| row | log2 instructions | kind |
|---|---:|---|
| relation collection (section 4) | 112.848 | fixed by the region |
| factoring and deduplicating kept relations | 107.784 | expected (H1), + 3 SE |
| Schirokauer maps | 87.911 | fixed |
| factor base | 72.774 | fixed |
| linear algebra (section 5) | 114.274 | fixed bound |
| two individual logarithms (section 6) | 111.233 | expected, 8 x model |
| wrapper (section 2) | 70.000 | bound |
| one run, total | 114.8636 | expected |

One run: 2^114.8636 instructions = 2^100.5759 units. With the gate, E[W] <= 1 + (13/32) x 2^100.5759 = 2^99.2763 < 2^99.28 units.

## 9. Memory, preprocessing, advice

Peak storage: sieve counters 4 bytes per cell, 2^107.4 bytes, plus hit queues, tables (2^57.2 bytes), relations
and matrix (2^61): below 2^108 bytes. The rules give memory no scalar contribution; the counters can be
processed in slices of a at the price of repeating table look-ups, which this bound does not use. All
preprocessing (everything before the individual logarithms, 2^114.742 instructions = 2^100.454 units) is for
this one target and charged in full. Nonuniform advice: none; h, f, g and the parameters are design choices.

## 10. Heuristics, support, limits

- H1 (smoothness). A norm of a sampled cell is B-smooth with the probability of a random integer of its size
  shifted by Murphy's alpha, and the two sides are independent given the sizes. Support: the standard NFS
  yield model; alpha is measured on exact production-size norms and charged at + 3 SE with no credit for a
  negative value; the exact local masses below 4096 are measured. Limit: no production-size smooth pair was
  exhibited; a yield error up to the 47% surplus changes nothing, a larger one lowers success.
- H2 (sampling). The uniform sample of 2,500 cells represents the region: yield at mean - 3 SE. [E3] redraws
  fresh cells for comparison of sizes and local masses.
- H3 (matrix). At most 2.2 B/ln B factor-base ideals; r divides neither class number; the 29 Schirokauer maps
  are independent on units modulo r-th powers; D distinct relations have a one-dimensional kernel on occurring
  columns. Standard for NFS discrete logarithms; not tested at toy scale here.
- H4 (ECM, work and detection). A schedule of 8 x exp(sqrt(2 ln P ln ln P)) modular multiplications finds a
  given prime below P with probability >= 0.99; a multiplication on k limbs costs 16k^2 + 64 instructions.
  Support: the standard ECM analysis, whose expected cost is exp((sqrt 2 + o(1)) sqrt(ln P ln ln P))
  multiplications; the factor 8 covers detection 0.99 (4.6 expected schedules) and the o(1). Charged: 2^50.5
  multiplications at 165 bits, 2^26.7 at 54 bits. Each trial carries 0.99 per possible
  prime factor above B.
- H5 (FFT, work only). A product of two N-bit integers costs at most 8 N log2 N instructions for N >= 2^16;
  the batch test of K numbers against the product of the primes below 2^L uses at most L + 6 log2 K + 16 such
  products of the full size.
- H6 (local density). Given that an n-bit norm is L-smooth, its prime factors of relative size t have density
  kappa(t) dt/t. Used for the expected number of children and for missing ideals.
- H7 (occurrence, encountered ideals). (a) A factor-base ideal of norm q occurs in each of the D relations
  independently with probability kappa/q, kappa the H6 ratio at the main phase's mean sizes; an ideal that
  never occurs has no logarithm, and a trial whose relation contains one is rejected. The expected number of
  such ideals per relation is 0.16 at the lowest bin and below that elsewhere; 1 - that number is charged
  (union bound). (b) The ideals met in the descent have the lattice statistics of random degree-one ideals of
  their size: the model value of p(Q) depends on Q only through the reduced lattice of the root data (s, t)
  mod q. (c) For every encountered ideal the realized fraction of relation cells in its ball is at least the
  charged one, i.e. at least the H1 mean over the ball divided by the allowances (1.5 and the row factor 8,
  together 12). Support: the admissible balls hold at least 64, and in the lowest bins about 2^12, expected
  relation cells, so (c) asks that the realized count of a ball not fall below one twelfth of its H1 mean;
  the six sampled ideals per bin differ by at most 2.7 bits in log2 p_i, sampling noise included; E[1/p] is
  estimated directly. Limits: no descent was executed and no ball was sieved, so the realized counts are not
  measured at any size; encountered ideals were not sampled from actual parent relations; six ideals per bin
  do not bound the worst encountered ideal. A failure of (b) or (c) raises work or triggers the meter. Without
  the factor 8 the row is 1% of the run; a further factor 4 beyond the charge would move the total by 0.3 bit.

What the experiment does and does not show. [E1]-[E6] check exact arithmetic at production parameters, replay
stored samples and recompute the ledger. They do not find a relation, run a toy sieve, or validate H1-H7.

## 11. Context

The method is the special tower number field sieve of Kim-Barbulescu with a twisted pair (f, g) for the
composite case gcd(6, 2) != 1 in the manner of Guillevic-Singh; published estimates for this field are context
only and no published total is a premise. The run gate and the convention of scoring expected total work
follow accepted practice on this benchmark's BLS12-381 tracks; no file or number of those submissions is
reused.
