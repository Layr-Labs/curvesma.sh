# SexTNFS (6,2) over Q(zeta_7) for BN254 pairing inversion

This is an exploratory cryptanalytic claim, not a formal proof certificate, a
practical break or a security-bit estimate. The track asks for any P in G1, Q in G2
with e(P,Q) = z for a uniform z in GT. We use the strategy the rules permit: solve the
GT discrete logarithm z = g^x by the special extended tower number field sieve, then
output **P = P0 and Q = [x]Q0**, so that e(P0, [x]Q0) = e(P0, Q0)^x = g^x = z by
bilinearity. The scalar multiplication [x]Q0 and all output work are charged below
(Sections 8b and the ledger). The GT-DLP analysis is identical to our bn254-gt-dlp
package; this document adds the pairing output. It gives a concrete, self-contained
parameterisation and a complete instruction-equivalent cost ledger. All numbers
are recomputed from the protected constants by `experiment.py`; every heuristic is
disclosed with its support and its limits in Section 10. No production-scale attack
was executed.

Submitted bound: **time_log2 = 108.000000** (W = R/20000). The ledger of `machine.md`,
recomputed by `experiment.py`, derives a total **2^102.437** units; the claimed 108 is
that derivation plus a disclosed **5.56-bit margin** (Section 9). Success probability is
at least **0.50** (the model derives 0.92), over uniform challenges and fresh coins.

## 1. Target and why index calculus applies

The profile fixes BN254 (Ethereum alt_bn128): p = 36u^4+36u^3+24u^2+6u+1 with
u = 0x44e992b44a6909f1 = 4965661367192848881 (p is 254 bits), r = 36u^4+36u^3+18u^2+6u+1,
and GT the order-r subgroup of F_{p^12}^*. `experiment.py` checks p, r, and r | Phi_12(p)
with r coprime to p^k-1 for k in {1,2,3,4,6}: GT is the degree-12 cyclotomic subgroup, so
a discrete log in GT is a discrete log in F_{p^12}^* restricted to the order-r subgroup.

F_{p^12} is a composite-degree (k=12) extension of a prime field whose prime p is the value
of a fixed low-degree polynomial (the BN polynomial P(x)=36x^4+36x^3+24x^2+6x+1 at x=u).
Composite k and a special p put BN254 in the best case for the extended/special tower NFS
(exTNFS/SexTNFS; Kim-Barbulescu 2016; Barbulescu-Duquesne, J. Cryptology 2019; Guillevic
2020). Pollard rho needs ~sqrt(r) ~ 2^127 group operations; index calculus is far below that.

## 2. Tower and polynomial selection (concrete)

We use the tower F_{p^12} = F_{p^6}[x]/(g), F_{p^6} = F_p[y]/(h), i.e. eta=6, kappa=2.

- **Base field.** h = Phi_7(y) = y^6+y^5+y^4+y^3+y^2+y+1, the 7th cyclotomic polynomial.
  Because ord_7(p) = 6 (checked), h is irreducible mod p, so F_p[y]/(h) = F_{p^6}. The number
  field K_h = Q(zeta_7) has class number 1, ring of integers O_h = Z[zeta_7] (a PID), and a
  fully known unit group (Washington, *Cyclotomic Fields*: h^+(7)=1, units are
  +/- zeta_7^a (1+zeta_7)^b (1+zeta_7+zeta_7^2)^c). This known unit group is used to measure,
  not assume, the per-cell unit multiplicity (Section 4).

- **Special twist.** Put eps = 2 + y in O_h (a non-unit, N(eps)=Phi_7(-2)=43). We require
  u*eps to be a non-square in F_{p^6} so that g = x^2 - u*eps is irreducible over F_{p^6};
  `experiment.py` checks (u*eps)^((p^6-1)/2) = -1 and also u^((p^6-1)/2) = +1, i.e. the
  untwisted x^2 - u would *not* be irreducible — the twist is necessary.

- **Special side (Joux-Pierrot).** f = eps^4 P(x^2/eps) = 36x^8 + 36 eps x^6 + 24 eps^2 x^4
  + 6 eps^3 x^2 + eps^4, a degree-8 polynomial over K_h with coefficients of controlled size
  (|N(eps^k)| = 43^k). f and g share the root x with x^2 = u*eps: substituting, f reduces to
  eps^4 P(u) = eps^4 * p = 0 in F_{p^6}[x]/(g); `experiment.py` verifies this common root mod p.

- **Irreducibility of f over K_h (degree 8).** P is irreducible over K_h because its quadratic
  subfield is Q(sqrt(-3)) while the only quadratic subfield of K_h = Q(zeta_7) is Q(sqrt(-7));
  since sqrt(-3) is not in Q(zeta_7), Q(theta) ∩ K_h = Q and P stays degree 4 over K_h. Then
  f is the degree-2 Capelli extension X^2 - eps*theta over K_h(theta); it is irreducible iff
  eps*theta is not a square in K_h(theta). We ship a Capelli certificate: a degree-one prime
  (q, y-s, theta-theta0) = (43, 16, 11) of K_h(theta) at which the Legendre symbol
  (eps(s)*theta0 / q) = -1, witnessing that eps*theta is a non-square; `experiment.py` checks
  it with one Legendre symbol (q=1 mod 7, s a primitive 7th root, theta0 a root of P mod q).
  A reviewer sanity check — "f factors into two quartics modulo every prime" — is expected and
  not a counterexample: the product of the four values x^2 = eps*theta_i equals eps^4/36 =
  (eps^2/6)^2, a square in K_h, which forces the Galois group into the even-sign-flip subgroup
  of C_2^4 >< Gal(P) that contains no 8-cycle, so no prime can give an irreducible degree-8
  reduction even when f is globally irreducible. The degree-one Capelli test is the correct
  decider and it certifies irreducibility.

The isomorphism from the profile basis F_p[w]/(w^12-18w^6+82) (with u -> w^6-9) to our tower
representation is a fixed linear/rooted map computed once (<2^40 instructions, in the fixed
cost F of `machine.md`); the attack reads only z and the public profile.

## 3. Factor base, region and smoothness bound

Smoothness bound **B = 2^54** on both sides. The factor base is *all* prime ideals of norm
<= B on each side, of any residue degree (in K_h only primes q = 1 mod 7 split completely,
2 has residue degree 3 etc.; the factor base includes these so no cell is silently lost).
By the prime-ideal theorem the count per side is ~ Li(2^54) = 2^48.81, so the matrix dimension
is D = 2 * 2^48.81 + 28 = 2^49.81 (the 28 Schirokauer columns of Section 7).

Relation search is over the box of (a,b), a,b in O_h with the 12 integer coordinates in
[-A,A], A = 255, up to the global sign: 2^106.97 cells. For each cell we form the two norms
N_f(a,b) = Res_y(Res_x(f, a - b x)) and N_g(a,b) = N_{K_h/Q}(a^2 - u eps b^2); both are exact
integers computed in `tnfs_common.norm_pair`. On 250 uniformly sampled production cells (shipped
in `data_cells.json`) the measured means are log2 N_f = 440.46, log2 N_g = 475.74 (maxima 477.16,
485.75).

## 4. Relation yield (measured norms + disclosed rho, calibrated)

Per-cell smoothness probability is rho(u_f) * rho(u_g) with u_f = log N_f / log B,
u_g = log N_g / log B (Dickman rho; heuristic H1, independence H2). Averaged over the 250
shipped cells, with the charged Murphy alpha below added to each log norm, this joint yield is
2^-54.05. We then charge, all measured on the production cells:

- **Murphy alpha (H3).** The small-prime valuation excess of each norm, averaged over cells;
  measured mean +1.24 nats (f) and +1.41 nats (g); we charge the mean plus two standard errors,
  +2.20 nats (f) and +2.18 nats (g). Positive alpha *increases* the effective norm,
  so we charge it (shrinks yield); we take no credit when alpha is negative.
- **Primitivity.** Fraction of cells with no prime ideal of O_h dividing both a and b: 0.924
  (measured). A non-primitive cell gives a duplicate, so we multiply yield by 0.924.
- **Unit multiplicity.** Each relation is counted by the whole unit orbit inside the box; the
  mean inverse multiplicity is 0.558 (measured with the known unit group; the largest unit
  exponent reached inside the box is 2, so the orbit is fully enumerated). We multiply by 0.558.
- **Calibration (H1/H2 support).** On a toy instance with the *same construction* at
  u = 6 (p = 55333, 16-bit), 65 regions x 3000 cells, the realised joint-smooth count over the
  predicted rho count is 2.157 (2-sigma lower bound 2.109) at B = 2^20, i.e. rho *under*-predicts
  the real yield; we therefore charge the full rho yield (calibration factor min(1, ratio) = 1.0)
  and never take the 2x bonus. The dispersion index of per-region counts is 1.12 and no region
  had zero relations (both used in Section 8).
- A further 0.5-bit sampling margin. Charged yield per cell: 2^-55.51.

Relations obtained: 2^106.97 * 2^-55.51 = 2^51.46, against the 1.25 * D = 2^50.14 needed
(25% excess for filtering/rank, H4). The region is oversized on purpose; the surplus factor
2^1.32 is margin.

Because the relation count is a sum of 2^106.97 near-independent indicators of mean 2^51.46,
its relative standard deviation is ~2^-25.7, so the probability of collecting fewer than D
relations is negligible (second moment; this needs only the measured norm distribution and H1/H2,
not any assumed count law).

## 5. Relation collection cost

Both sides are lattice-sieved with their complete factor bases (every prime ideal of norm <= B,
including the large ones). The per-cell
work is the number of factor-base hits per side times the per-hit update, plus the scan: the
expected number of factor-base ideals dividing a cell's norm is sum_{N(Q)<=B} 1/N(Q) =
ln ln B + M + o(1), charged as ln ln B + 0.2615 + 0.3 ~ 4.7 per side; large ideals
(N(Q) > 2^36) carry the enumeration-tree surcharge. This gives 136.9 scalar instructions per
cell (derived in `machine.md`), hence sieve work 2^106.97 * 136.9 = 2^114.21 instructions
(the coset/root setup per (ideal, block) is included and is 2^110.8, subdominant).

Because every factor-base ideal is sieved on both sides, a cell's two sieve sums equal the logs
of its B-smooth parts up to rounding, so a cell passing both thresholds is smooth up to that
rounding slack (every prime above B exceeds 2^54); survivors are therefore essentially the
relations, and charging 16x the relation count is slack. Each survivor
is cofactored by a **fixed** ECM schedule — trial division, then stage-1 ECM with bound B1 = 2^20,
78 curves (sized so a 54-bit factor is missed with probability <= 2^-10), discarding any survivor
that does not split (a charged failed trial, never "repeat until success"). The measured per-curve
detection rate on 300 planted 24-bit factors (B1 = 2^10) is 0.321 against the rho prediction 0.156
(2.05x; H7; 300/300 found), so charging rho curves is conservative. Cofactor work: 2^100.51
instructions. Schirokauer-map evaluation on the kept relations: 2^83.91. Total relation collection (dominated by
the sieve): **2^114.21 instructions = 2^99.91 units**.

## 6. Linear algebra

A sparse system of dimension D = 2^49.81 over Z/r (r is 254 bits, four 64-bit limbs) is solved by
block Wiedemann. We take **no filtering credit** (dimension stays the full factor-base size,
conservative). Average row weight is 130 small-prime entries plus 28 dense Schirokauer columns;
`machine.md` prices one row's matrix-vector contribution at 19560 instructions (130*20 for the
small entries, 28*600 for the dense Schirokauer columns at one full 254-bit modmul each, 160
overhead). Block Wiedemann does ~3D matrix-vector products: 3 * D^2 * 19560 = **2^115.47
instructions = 2^101.17 units**. (The
Berlekamp-Massey/Thome generator step is O(D^2 log D) bit-operations, subdominant and included.)

## 7. Schirokauer maps and virtual logarithms

r is unramified in K_f and K_g; the unit ranks give 23 Schirokauer columns on the degree-8 side
and 5 on the degree-2 side (28 total), appended as dense columns mod r. K_h has class number 1,
so there is no class-group obstruction on the base; the number fields K_f, K_g may have nontrivial
class groups, handled in the standard way by including a basis of their r-torsion characters in the
Schirokauer system (H5, disclosed). The kernel vector gives the virtual logarithms of all
factor-base ideals; a final check (Section 8) catches any map deficiency as a failed run, never a
wrong answer.

## 8. Individual logarithm: boot and special-q descent (concrete, bounded)

**Boot.** Multiply the target z (and separately g) by random powers of g and reduce the lift in
the g-side number field K_g via LLL on the dimension-12 lattice of determinant p^6 (the ideal of
elements congruent to the target class mod p). A rigorous bound on the resulting g-side norm is
6 * (2*(log2 6 + 2.75 + 127) + log2(1+3u)) = 1970.2 bits (LLL gives coordinates <= 2^2.75 * p^(1/2);
each of the 6 conjugates of a^2 - u eps b^2 is bounded accordingly); the measured maximum over 8
production boot lattices is 1892.1 bits, below the bound. We smooth this ~1970-bit number to the
bound 2^200 by a **fixed** budget of 2^36.76 reduced candidates, each tested by a fixed stage-1 ECM
schedule (near-optimal B1), sized so the whole-boot failure probability is <= 2^-8 under rho
(H1, H8). Each candidate is charged 37 rounds of that schedule. A successful boot expresses the target
via at most 37 prime ideals of norm in (2^54, 2^200] (each exceeds 2^54 and the norm has at most
1970.2 bits), and every one of them is descended; smaller ideals are in the factor base. Boot work
(both boots): 2^111.14 instructions = 2^96.85 units.

**Descent.** Each of the at most 2*37 = 74 top special-q ideals is rewritten over four levels
200 -> 150 -> 110 -> 80 -> 54. At each level and ideal we build the special-q lattice
L(q,s,t) = { (a,b) : a(s) + t b(s) = 0 mod q }, an explicit dimension-12 integer lattice of
determinant q (rows: q e_0; e_i - s^i e_0; e_{6+j} - t s^j e_0), LLL-reduce it, and sieve a box of
side E in the reduced basis. The special-q divides the chosen side's norm exactly (removed before
smoothness testing). We sampled, at production size, the exact norm distributions of these balls at
q-sizes 80, 110, 150, 200 bits on both sides and several E (shipped in `data_spq.json`): e.g. at a
200-bit g-side special-q (E = 1) the f-norm averages 892 bits and the g-cofactor 388 bits. For each level the
ledger picks the harder side and the smallest E giving a mean smooth count mu >= 2^10. It charges
sieving the whole box with the factor base (at the bottom level, whose bound is B itself, the sieve
identifies the smooth cells directly and 16*mu survivors are cofactored, as in Section 5) and then,
at every level above 2^54, a fixed ECM schedule
on **every** cell of the box: no survivor filter is assumed, and each cell gets
ceil((log2 N_f + log2 N_g)/54) ECM rounds, each stripping one factor above 2^54. A node makes at most
**4 attempts** (hard per-node work bound): attempt r = 0..3 sieves the box translated by r(2E+1) along
the first reduced basis vector, so the attempts cover disjoint cell sets. Per-attempt failure to find
a relation is bounded by the second moment, P(count = 0) <= Var/mu_a^2 <= (charged dispersion)/mu_a,
with the charged dispersion 2.24 = 2 x the measured toy dispersion index 1.12 and mu_a = mu/4 (a
translated box has larger coordinates than the centred one, so each attempt is charged a quarter of
the centred mean count). Treating the four disjoint attempts as independent given q (H6), per-node
failure is <= (4 x 2.24/mu)^4. The number of children per node is at most
ceil((log2 N_f + log2 N_g)/54), because every prime ideal above the factor-base bound 2^54 is
descended and each exceeds 2^54; every child is charged a full node at the next level, although
smaller children are cheaper. With these worst-case counts the descent has at most 2^19.95 nodes,
the total descent failure is < 2^-11.6, and the descent cost is 2^115.35 instructions = 2^101.07
units (per-level figures in `machine.md`).

Nothing in boot or descent repeats until success: every loop (candidates, curves, retries) is a
fixed integer cap, and every cap's failure probability is in the success budget of Section 11.

## 8b. Pairing output P0, [x]Q0 (charged)

Having recovered x = log_g z mod r, we return P = P0 (the public G1 generator) and Q = [x]Q0,
computed by binary double-and-add in affine Fp2: at most 254 doublings and 254 additions, each a
constant number of Fp2 multiplies and one Fp2 inversion (Fp2 = Fp[u]/(u^2+1) in the profile's notation, where u is the Fp2 generator and not the BN
parameter; inversion is (a-bu)/(a^2+b^2), one Fp inverse). Each point operation is < 2^24 instructions, so [x]Q0 and the
optional on-curve and [r]Q0 = O subgroup checks are < 2^34 instructions. Q0 has order r, so [x]Q0
is in the prescribed G2 subgroup; P0 is in G1. e(P0,[x]Q0) = g^x = z holds by bilinearity, so no
pairing evaluation is needed by the attack (the organizer's pairing check is excluded by the rules).
For z = 1 we have x = 0 and output P0 and Q = infinity, permitted by the profile. This < 2^34 output
cost is inside the fixed term F = 2^50 of the ledger, so the pairing-inversion total equals the
GT-DLP total to the ledger's precision.

## 9. Total, score and disclosed margin

Summing all instruction counts (linear algebra 2^115.47, descent 2^115.35, sieve 2^114.21,
boot 2^111.14, cofactoring 2^100.51, Schirokauer 2^83.91, and the fixed setup/verification/
output F = 2^50): total R = 2^116.725 instructions. Dividing by C = 20000:

    W = R / 20000 = 2^102.437  units,   log2 W = 102.437453.

We submit **108.000000**. The 5.56-bit margin above the derivation covers, conservatively and
without claiming any of it is needed: the first-order rho model versus the true smoothness law
(the toy ratio shows rho under-predicts, so this is slack in our favour), the independence
heuristic H2, the class-group/Schirokauer conditions H5, any under-pricing of large-ideal sieve
updates or the generator step, and model error in the descent ball sizing. Independent literature
estimates for BN254 are ~2^100 in their own operation counts (Barbulescu-Duquesne 2019; Guillevic
2020); those units differ from ours, so this is only a plausibility check against our derived
2^102.44 units, and the claimed 108 sits ~8 bits above it.

## 10. Heuristics, support and limits

- **H1 (Dickman rho).** Smoothness probability of each norm is rho(log N / log B). Support:
  exact norm measurement on production cells and special-q balls; toy calibration ratio 2.16
  (rho conservative). Limit: rho is asymptotic; we disclose it and charge the full (not bonus) yield.
- **H2 (two-sided independence).** Joint smoothness = product of per-side rho. Support: toy joint
  realised/predicted ratio ~2.1 (no anticorrelation seen). Limit: not a theorem; margin covers it.
- **H3 (Murphy alpha).** Measured per side; positive charged, negative not credited.
- **H4 (rank/excess).** 25% relation excess yields a full-rank system mod r. Standard NFS heuristic.
- **H5 (Schirokauer/class group).** r unramified; K_h has class number 1; K_f, K_g class/unit
  conditions handled by the 28-column Schirokauer system. Disclosed as a hypothesis; a run that
  violates it fails the final check rather than returning a wrong x.
- **H6 (descent yield & lower tail).** Each node's ball contains a relation with the measured yield;
  the zero-probability is bounded by the measured dispersion via the second moment, and a node
  makes four attempts on disjoint translates, each charged a quarter of the centred mean count and
  treated as independent given q. Support: toy dispersion index 1.12 and zero-count fraction 0 over
  65 regions. This is the measured-distribution form the exploratory policy asks for, not an
  assumed count law. Limits: the conditional independence of disjoint attempts is a heuristic, and
  the per-node yield is the mean over the sampled special-q of each size, side and E (30 samples
  each), not a bound for every q.
- **H7 (ECM detection).** Per-curve success >= rho prediction. Support: 300/300 planted 24-bit
  factors found; measured per-curve rate 2.05x rho. Charged at rho (conservative).
- **H8 (boot smoothing).** The target's g-side lift reaches a 2^200-smooth value within 2^36.76
  candidates. Support: measured boot-norm size; rho yield. Limit: rho heuristic, disclosed.

Reproducible seeds in `data.json` are experimental evidence at their stated scope, not a proof of
ideal independent coins.

## 11. Success probability

The work bound is a hard cap on every run (all loops are fixed integer caps). The run returns the
correct x unless a charged cap is exhausted: relation shortfall (< 2^-25, Section 4), a 0.05
allowance for linear-algebra rank deficiency (H4) and missed ECM factors, boot failure
(<= 2 * 2^-8), descent-node failure (< 2^-11.6, Section 8), plus a 0.02 allowance for
Schirokauer-map deficiency (caught by the final g^x = z check). The derived success lower bound
is 0.922; we claim **>= 0.50**, comfortably above
the required 0.39. The identity z = 1 returns x = 0 directly.

## 12. Memory

Peak storage is dominated by the relation set and the sparse matrix: ~2^51.5 relations x ~130
entries x 8 bytes ~ 2^61.5, the factor base ~2^49.8 x 16 bytes ~ 2^53.8, and block-Wiedemann
vectors ~D x (few) x 32 bytes ~ 2^55. Peak total <= **2^62 bytes**. Memory is reported, not an
objective; it carries no scalar contribution. Nonuniform advice is zero (all precomputation is the
charged attack itself; public profile constants and algorithm design are excluded).

## 13. Correctness and output

The virtual logarithms of the factor base, composed with the descent of z, give x = log_g z mod r.
The attack verifies g^x = z before forming the output (one F_{p^12} exponentiation, in F); a mismatch
is a failed run, never an accepted wrong answer. This track then outputs P0 and [x]Q0 (Section 8b),
a valid pair with e(P0,[x]Q0) = z in the prescribed subgroups. All arithmetic is over O_h = Z[zeta_7],
modulo r, and in affine Fp2 with explicit limb counts in `machine.md`.
