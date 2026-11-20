# recover-P (bls12-381) through two GT logarithms with a (4,3) tower over Q(zeta_15) behind a run gate: expected work 2^109.179063486123

Track bls12-381-fixed-g2-invert-g1. Input: a nonidentity Q in G2 and z in GT; output P in G1 with e(P, Q) = z for the
exact supplied Q (P = infinity for z = 1). The score unit is one dense F_p12 multiplication, C = 20,000 scalar AArch64
instruction equivalents (ARMv8.0-A, general registers only). No advice, smaller targets or amortization are used.
No prior verdict is a premise.

## Claim

W is the total classical work of the one attack on the one target. Step 0 is a run gate: draw a 64-bit word u; unless
u < K = 7208358571546957167, stop as a failure (< 1 unit). K is the least integer with K x 0.998040 >= 0.39 x 2^64. Under H1-H5,
H-EA, H-R, H-X and H-SM (section 9) the algorithm below returns such a P with probability >= 0.998040 once past the gate
(section 7); u is independent of the challenge and the other coins, so success is at least 0.39 over uniform challenges
and fresh coins, and

    E[W] <= 1 + (K/2^64)(2^110.534687 + 2^66) < 2^109.179063486123  dense-Fp12 units (C = 20,000).

2^110.534687 is the sum of section 8's rows: collection and linear algebra are the same in every run; the two logarithms
and the ECM split are at expected work. The score is an expectation, not a bound on every run: a hard work meter stops
every opened run at 2^123.0 units (section 6), so all work of a run is at most 2^123.0 + 2^66 units.
Analytic, exploratory claim. Bounds round in the safe direction; `audit_cm.py` checks each against its unrounded value.

## 1. Reduction (fixed Q: two GT logarithms)

With the pinned reduced pairing e and the public generator P0 of G1: (0) run gate; (1) check Q in G2 \ {infinity} and
z^r = 1, return infinity if z = 1; (2) g = e(P0, Q); map g, z into F_p[y, X]/(h_pi, X^3 - omega^2 u) (field-map.md);
(3) draw uniform e_1, e_2 in [0, r), solve the core for log(g^e_1) and log(z^e_2) in a common base, set
a = (e_1 log z^e_2)(e_2 log g^e_1)^-1 mod r; (4) output P = [a]P0 iff e(P, Q) = z. Wrapper < 2^66 units, failure < 2^-48.

## 2. Tower and polynomials

u = -0xd201000000010000, r = u^4 - u^2 + 1, p = (u - 1)^2 r/3 + u. Base field K_H = Q(zeta_15), O_H = Z[zeta_15]
(class number 1), omega = y^5; Phi_15 = h_pi h_pibar mod p and P = (p, h_pi) has residue field F_p4.
g(X) = X^3 - omega^2 u is irreducible over F_p4 and divides f(X) = X^9 + (omega + 1) X^6 + X^3 - omega modulo P
(both checked in `algebra_bls12`), so F_p12 = F_p4[X]/(g) (field-map.md). sigma(y) = y^7 fixes P, f and g, and tau(y, X) = (y^7, omega^2 X) acts as Frobenius
(chi = p mod r, order 12): ideals are taken in tau-orbits of 12, and the matrix is over Z[zeta_12]. Only full-orbit
relative-degree-1 ideals enter the factor base; short orbits have log 0 mod r (ramified primes 3, 5, 19, 67 and the
primes of u are checked). Murphy alpha: f +2.74, g +1.28 bits. Schirokauer maps: Phi_15 = h_(r,1) h_(r,2) mod r; f stays
irreducible of degree 9 over both, g splits [4, 4, 4, 12]; P_chi Lambda has 3 (f) + 1 (g) chi-unit coordinates, and
uniform W_f (3 x 6), W_g (1 x 2) drawn before collection compress them to 4 columns (loss <= 4/r under H-SM).
`toy_e2e_dlp` runs the whole pipeline on u_0 = -5, p_0 = 7207, r_0 = 601: all relations' dependencies annihilate the
BSGS logarithms and 5/5 unseen targets are recovered.

## 3. Relation collection

Factor bases B_f = 2^71.29, B_g = 2^71.69 (U' = 2^63.317373 orbits), excess 1.09: target 1.148 U' = 2^63.516495
relations. The sieve region is the product of a Voronoi fundamental domain for O_H^* (10 units) and a 1/12 domain
for <omega, sigma>, cut by the exact 2-adic class rule (v_2(N_g) from the valuations at the two primes above 2;
three classes of radii 1.00, 1.20, 1.36). Pool: 200,000 samples, kept volume 2^119.5919; retained yield 2^63.516495
with standard error 0.0102 bits and a margin >= 0.060 bits above 1.09 U'. `verify_pool_2adic_live` redraws 400 cells
and recomputes volume, class shares, norms and the de Bruijn yield.

Sieve (KERNELS.md): per-row lattice sieve over replayed tables of the a-orbits (build 2^106.239198 units, storage
2^119.65 bytes), band 1 below 2^55, band 2 up to the factor-base bound with a 64-entry super-bucket gate, 16-bit
accumulators, then cofactorisation (schoolbook norms, Montgomery ECM, 4-prime NTT products for stage B). 24.1105
instructions per cell in total; collection <= 2^109.907365 units. Band split <= 2^100.543. ECM split of kept
relations: curves are drawn from the run's fresh coins and a schedule is repeated until the cofactor splits, so the
row is an expectation over the coins: 2^95.09 expected, charged 2^96.2.

## 4. Filtering and linear algebra

Z[zeta_12] = Z[omega] + chi Z[omega]; filtering (purge, clique removal, merge to weight 200) runs on the tau-quotient.
`filtsim_z12.py` simulates it at U' = 2^16..2^23 (6 seeds each): N_l/U' falls from 0.003260 to 0.001889 and N_sw/U'
from 0.006462 to 0.003846 at every doubling; unit share 97.45% at 2^23. Caps: N_l <= 0.00225 U', N_sw <= 0.00450 U',
k5 <= 17.50, table ratio <= 0.85 (above them the run aborts; H3).

Block Wiedemann on one owner (x, y and 22 tables in 2^63.913 bytes), (m, n) = (8192, 110), 4 Schirokauer columns:
<= 2^107.789554 units (unit entries 8.80, non-unit entries <= 4 basis updates of 29 instructions, Barrett row
reduction 68). Kernel recovery and back-substitution: RANK-BRIDGE.md; `hx_z12.py` runs the filtered operator and the
reconstruction on the live system at 2^12..2^14.

Matrix generator. The generator of the m x n sequence of length L = (1 + 2^-7)(N + 4)(1/m + 1/n) + 6 is an order basis of the k x k
(k = m + n) polynomial matrix [A(z) | I] to order L, computed by the divide-and-conquer PM-basis: solve to order L/2,
multiply the residual by the half basis, solve the second half, multiply the two bases. Each of the log2 L levels
multiplies polynomial matrices whose scalar products have total length <= 4 k^3 L; each scalar product of length l is
one integer FFT product of Kronecker-packed 255-bit coefficients (18 words each: product coefficients < 2^(510 + log2 l))
priced by fft_F (KERNELS.md section 3, 4-prime NTT) plus 640 instructions per coefficient for Barrett reduction mod r:
<= 4 k^3 log2(L) (fft_F(18 L) + 640 L) instructions = 2^98.584 units. `bm_toy.py` (audit) runs the recursion and the
iterative M-basis on random toy inputs: identical bases, F P = 0 mod z^d, and length-weighted product mass 9664, 46275, 156456
against the charged 16384, 80000, 331776.

## 5. Main-phase success

Region, polynomials and the sieve's ECM seeds are fixed, so the relation set and the filtered matrix are fixed (the
split of kept relations into ideals uses fresh coins, section 3). Under H-R the
kernel of [M | S] is the line of virtual logarithms; under H-X Block Wiedemann fails with probability
<= 4 N/2^128 < 2^-66 (Schwartz-Zippel; entries of X, Y uniform on >= 2^128 values).

## 6. Individual logarithms and the meter

DESCENT.md: 78 nodes in 39 bins and the boot; expected two-log work <= 2^107.8098 units. One meter counts all
work of an opened run and stops it at t = 2^123.0 units (checked every <= 2^66 units). The meter changes no expected
cost: E[min(W_open, t)] <= E[W_open] <= 2^110.534687, and P(meter) <= 2^(110.534687 - 123.0) <= 0.00017684 (Markov).

## 7. Success

Opened run: failure <= 2^-66 (BW) + 2^-80 (Miller-Rabin) + 0.0017831 (node counter and
repeats, DESCENT.md) + 0.00017684 (meter) + 6/r (degenerate targets, compression) + 2^-48 (wrapper); success
>= 0.998040. Subtrees fail independently (H5). With the gate: P(success) >= (K/2^64) x 0.998040 >= 0.39, and
E[W] <= 1 + (K/2^64) E[W_run].

## 8. Ledger

| component | log2 units (C = 20,000) |
|---|---:|
| relation collection (section 3) | <= 109.907365 |
| replay tables build (section 3) | <= 106.239198 |
| linear algebra (section 4) | <= 107.789554 |
| matrix generator (section 4) | <= 98.584 |
| both individual logarithms (expected; DESCENT.md) | <= 107.8098 |
| band split of kept relations | <= 100.543 |
| ECM split of kept relations (expected) | <= 96.2 |
| canonical a list and sort | <= 86.0 |
| Schirokauer maps (4 columns) | <= 84.0 |
| filtering, back-substitution, misc | <= 80.0 |
| P_g, P_f product trees | <= 75.5 |
| factor base, roots, ideals | <= 75.0 |
| numerical certificate | <= 50.0 |
| verification | <= 11.0 |
| wrapper (section 1) | <= 66 |

Sum <= 2^110.534687.

## 9. Heuristics

H1 smoothness (de Bruijn rho, Murphy alpha, exact v_2 at the two primes above 2). H1' retention (two-stage Euler
product). H-EA (work only): the survivor law prices stage B. H2 region (pool estimate, 0.774% allowance). H3 matrix
(r unramified, r prime to the class numbers, H-SM; the section-4 caps hold at production size). H4 kernel cost
(KERNELS.md, STORAGE.md). H-R / H-X (RANK-BRIDGE.md): rank and Block Wiedemann minors for this fixed system.
H5 descent (DESCENT.md): the negative-binomial recurrence of section 2 with independent subtrees of distinct ideals and
outcome-independent order, measured in DESCENT.md section 3.

## 10. Memory

Peak < 2^120.0 bytes (replay tables 2^119.65; linear algebra 2^63.913 on one owner; descent < 2^105).

## Provenance

Descent recurrence, run gate, meter and ledger layout follow the accepted c36p (50b04508) and its lineage
(kaankolcu, mpjunior92, gbotrel). Polynomials after Guillevic-Singh; tower and filtering after Kim-Barbulescu and
CADO-NFS. (4,3) tower over Q(zeta_15), Z[zeta_12] linear algebra: gbotrel with Claude.
