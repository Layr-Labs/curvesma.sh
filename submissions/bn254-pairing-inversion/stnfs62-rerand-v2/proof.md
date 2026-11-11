# BN254 pairing inversion through the GT discrete logarithm: (6,2) special tower NFS with a re-randomised individual logarithm and fully swept nodes

Track `bn254-pairing-inversion`. Unit: C = 20,000 scalar AArch64 instruction equivalents. This is an analytic, exploratory
claim (`plausible_not_refuted`): no production attack was run. Success and expected work rest on four named
heuristics (H-main, H-boot, H-desc, H-ECM; section 11), each stated with its support and limits. Nothing is
amortized over other targets and no advice is used. Every number below is recomputed by `model.py` from the
production-size samples in `data.json`; `experiment.py` reruns that computation, compares it with `declared.json`,
with this text and with the manifest, and reruns the checks cited as [E1]-[E9].

## Claim

W is the total classical work of the one attack on the one target, in units. Step 0 is a run gate: draw six fresh
uniform bits; unless their value is below 27, stop as a failure (below one unit). Past the gate the
algorithm returns P in G1 and Q in G2 with e(P, Q) = z with probability at least 0.9374 under the heuristics (section 8). The gate is
independent of the challenge and of all other coins, so success is at least (27/64) x 0.9374 = 0.3954 >= 0.39
over uniform challenges (the identity included, answered directly) and fresh coins, and

    E[W] <= 1 + (27/64) x 2^108.288 = 2^107.043 < 2^107.05 units.

The score 107.05 is an EXPECTATION over the challenge and the coins, not a bound on every run. No every-run bound
is claimed for the whole attack; the only hard cap is the meter of section 7 (a run never spends more than
2^126.455 instructions on the two individual logarithms).

What is new against the rejected tower-NFS claims on this track (the construction of the field, of f and g and of
the matching sieve follows submission f8a14ba8-9a8a-4c47-91f8-29b129f2dbd6, declared as a dependency; every
number here is recomputed from our own samples):

1. The individual logarithm assumes nothing about a particular ideal. Every trial re-randomises the TARGET with
   a fresh uniform exponent, so a trial (boot, then one descent) is an independent draw from one
   challenge-independent law. A node without a usable relation is not a non-termination and is not repaired: it
   ends the trial and a new uniform target is drawn. The expected work is an exact Wald identity in three numbers
   attached to uniform targets (section 7). The premise on the descent is a RATE: a visited node has an empty
   ball in at most a fraction eps = 1/1000 of the cases, an average over uniform targets; empty balls are
   allowed and paid for (the ledger charges that only one descent in two succeeds), instead of assuming that
   every reached ideal has a relation. The descent has two levels and at most 644 nodes by a
   deterministic count, so pi_d >= (1 - eps)^644 with no model of the tree.
2. A node does not sample until success. It sweeps one fixed ball completely at a fixed price, so its cost does
   not depend on the realized relation density of its ideal. Every plan has at least 2^10 model-expected usable
   relations in the ball of the worst of 16 sampled ideals. Realized relation counts of completely swept balls,
   with the same predicate, are measured on a toy tower ([E9]).
3. No FFT, no asymptotic ECM cost, and no ECM in the descent at all: every node decides smoothness by the exact
   sieve itself. ECM appears only in the boot, as a stage-1 schedule (B1, number of curves) with an exact
   instruction count whose failures only drop the candidate; its capture probability is measured ([E7]).
4. No exceptional ideals: cells are restricted so that every prime ideal in a relation is given by a simple root,
   and Schirokauer maps are used in full degree, so no unit rank or maximal order computation is a premise.

## 1. Algorithm overview

Input z in GT (order r). Output P in G1 and Q in G2 with e(P, Q) = z; the algorithm recovers x with g^x = z, g = e(P0, Q0) the protected generator, and returns (P0, [x]Q0).

0. Run gate.
1. Decode z, check z^r = 1. If z = 1 return (P0, infinity).
2. Field map (section 2) to F_p[y,x]/(h, x^2 - u y).
3. Main phase, independent of z: factor base, relation collection (section 5), linear algebra (section 6). Output:
   the table KNOWN of virtual logarithms.
4. Individual logarithms of g and of z (section 7): vlog(g), vlog(z); x = vlog(z)/vlog(g) mod r.
5. Check g^x = z by exponentiation in the protected representation (a wrong value is a failure). Compute
   Q = [x]Q0 by double-and-add in G2 (at most 508 group operations, each fewer than 2^12 instructions: below
   2^21 instructions) and output (P0, Q) in the canonical encoding; e(P0, [x]Q0) = g^x = z by bilinearity, and Q
   lies in the order-r subgroup because Q0 does. This output step is inside the 2^70 instructions of the wrapper row.

## 2. Tower, polynomials, field map ([E1])

u = 0x44e992b44a6909f1, p = 36u^4 + 36u^3 + 24u^2 + 6u + 1, r = 36u^4 + 36u^3 + 18u^2 + 6u + 1.
K_h = Q[y]/h, h = y^6 + y^4 - 2y - 1, O = Z[y] (we only use the order Z[y]; disc(h) = 2^6 3^3 197, so Z[y] is
maximal at every prime other than 2, 3, 197, which are excluded below). h is irreducible modulo p ([E1]).
With P(X) = 36X^4 + 36X^3 + 24X^2 + 6X + 1 (P(u) = p):

    f(x) = y^4 P(x^2/y) = 36x^8 + 36y x^6 + 24y^2 x^4 + 6y^3 x^2 + y^4,        g(x) = x^2 - u y.

At a root of g, f = y^4 P(u) = y^4 p = 0 mod p ([E1] checks the identity). u y is a non-square in F_(p^6) ([E1]),
so F_p[y,x]/(h, g) = F_(p^12). K_f = K_h[x]/f has degree 48, K_g = K_h[x]/g degree 12.
r divides p^4 - p^2 + 1 and gcd(r, p^6 - 1) = 1 ([E1]): elements of F_(p^6), the image of Z[y], have trivial
component in GT, so they do not change logarithms modulo r.

Field map: find a root of h and a square root of u y in the protected field F_p[w]/(w^12 - 18w^6 + 82) by
equal-degree splitting (each trial fewer than 2^17 field multiplications of at most 2^15 instructions; at most
2^10 trials, failure below 2^-400), then invert a 12 x 12 matrix. With decoding, z^r, the final check and output
this is far below the 2^70 instructions charged as "wrapper".

## 3. Admissible cells and their prime ideals (no exceptional ideals)

A cell is (a, b) in Z[y]^2. N_f(a,b) = |N_(K_h/Q)(sum_i f_i a^i b^(8-i))| and N_g(a,b) = |N_(K_h/Q)(a^2 - u y b^2)|
are the absolute norms of a - b x in K_f (up to the leading coefficient) and K_g (`lat.norm_f`, `lat.norm_g`).
Bad primes: S_f = {2, 3, 197}, S_g = {2, 3, 197, 13, 29, 983, 11003, 405928799} (u = 3 x 13 x 29 x 983 x 11003 x
405928799 [E1]; disc(P) = 2^10 3^7; the leading coefficient of f is 36). A cell is ADMISSIBLE when

    gcd(N_h(a), N_h(b)) = 1,   gcd(N_f, S_f) = 1,   gcd(N_g, S_g) = 1        (`lat.admissible`).

Only admissible cells are used, in the main phase, in the boot and in every node. For such a cell and a rational
prime q dividing N_f (resp. N_g): q is unramified in K_h and does not divide the index of Z[y]; for a prime Q of
Z[y] above q (given by an irreducible factor of h mod q), Q divides the relative norm F(a,b) iff b is a unit
mod Q and t = a/b is a root of f (resp. g) in the residue field Z[y]/Q; that root is simple, because the roots of
f are the square roots of y rho_i (rho_i the roots of P), which collide only modulo primes dividing
2 x 36 x disc(P) or, for g, 2u. Hence the only prime ideal of K_f (resp. K_g) above Q dividing a - b x is
(Q, x - t), with valuation v_Q(F(a,b)); the leading coefficient 36 is a unit at Q. The valuations v_Q are read
from v_q(N) (one prime above q divides) or by lifting the factor of h modulo q^k (several). So the ideal
factorisation of an admissible cell is elementary, and no maximal order of K_f, class group or unit computation
is used. Measured admissible fraction of uniform cells: 0.246 ([E4] re-measures it).

## 4. Factor base, Schirokauer maps, virtual logarithms

B = 2^54. The factor base holds, for each side, every prime ideal (Q, x - t) as in section 3 of norm N(Q) <= B
(residue degree 1 or more). By the prime ideal theorem each side has li(B)(1 + o(1)) of them and
li(B) = 1.028 B/ln B at 2^54; the ledger charges N = 2.2 B/ln B = 2^49.912 columns in total (H-main).

Schirokauer maps in full degree: for a side with defining field of degree n (48 and 12), lambda(gamma) is the
coordinate vector of (gamma^E - 1)/r mod r computed in (Z/r^2)[y,x]/(h, f) (resp. g), E a fixed exponent below r^n
killing the unit group of the quotient modulo r: 60 maps in total. No unit rank is used. A vector
(v, chi) = (one value per ideal, one per map) is a virtual logarithm when for every admissible cell

    sum_Q v_Q val_Q(a - b x_f) + chi_f . lambda_f(a - b x_f) = sum_Q v_Q val_Q(a - b x_g) + chi_g . lambda_g(a - b x_g)   (mod r),

and then (standard: Schirokauer; needs that r divides neither class number and that the maps are injective on
units modulo r-th powers, both part of H-main) each side equals the logarithm, to a fixed base, of the GT
component of the image of a - b x. With full-degree maps the solution space has dimension larger than one, but two
solutions differ by a functional vanishing on all values of actual elements, so any solution whose logarithm
component is nonzero serves (a random kernel vector fails this with probability 1/r).

## 5. Relation collection

Cells: all twelve coefficients in [-A, A], A = 330, one of each pair +-(a, b): V = (2A+1)^12/2 = 2^111.423.
A cell is a relation when it is admissible and both norms factor over the factor base.

Sieve (residue matching; no lattice enumeration, no assumption on the shape of the hit set). For each factor-base
ideal with root data (q, s, t), the cells it divides are those with a(s) = t b(s) in the residue field. The residues
a(s) of all (2A+1)^6 vectors a are generated in Gray-code order (one modular addition each) into a table keyed by
residue; for each of the (2A+1)^6/2 vectors b the residue t b(s) is looked up; every match is a hit and adds a
16-bit fixed-point log2 N(Q) to the cell's counter of that side. Counters live in chunks of 2^32 cells with
64-bit local addresses; a hit is a record (two-word cell index, value) appended to the chunk's queue.
Primes of S_f, S_g are not sieved.

Thresholds. Slack 47 bits. Stage 1 keeps a cell when its f counter is at least n_f,min - 48 and its g counter at
least n_g,min - 48 (n_f,min = 414.6, n_g,min = 456.0 bits are the 0.5% quantiles of the sampled norm sizes;
smaller norms are not counted in the yield). Stage 2 computes both exact norms of a kept cell, tests
admissibility, and keeps it when each counter is at least log2 N - 47. Lemma: distinct hit ideals each contribute
at least their log2 norm to log2 N, so a counter never exceeds the log of the B-smooth part of N; a stage-2
survivor has a cofactor below 2^47 < B on each side whose prime ideals, by section 3, have norm at most the
cofactor: it factors over the factor base. There are no false survivors. A relation is lost only if repeated
factors hide more than 47 bits; [E4] measures this mass on fresh admissible cells and the ledger charges a
retention of 0.95. Stage 2 is capped: at most 2^-20 V cells are examined (model expectation
2^82.2 against the cap 2^91.4); that the cap suffices is part of H-main.

Yield (`model.ledger`). Sample: 3,000 uniform cells of the box, 737 admissible, exact norms. For each
admissible sampled cell with norm sizes n_f >= n_f,min, n_g >= n_g,min,

    y = rho((n_f ln 2 + alpha_f)/ln B) x rho((n_g ln 2 + alpha_g)/ln B) x 0.50,     y = 0 for all other draws.

- rho: Dickman's function (table from u rho(u) = int_(u-1)^u rho, step 1/512; [E2] prints spot values).
- alpha (Murphy), measured on the admissible sampled norms and charged at + 3 SE: alpha = sum_(q<1000) ln q/(q-1)
  - E[ln of the 1000-smooth part of N]: alpha_f = 3.63, alpha_g = 2.83 (large because admissible norms
  carry no factor 2 or 3: the restriction of section 3 is paid here, and in the admissible fraction).
- 0.50 is the charged mean inverse unit multiplicity: cells eps (a, b), eps a unit of Z[y], give the same
  relation; measured mean 1/m = 0.791 (A = 232) and 0.787 (A = 290) over admissible cells with the units
  y^i e2^j e3^k, |i|,|j|,|k| <= 4 ([E4] re-measures it on fresh cells of the production box). The table of units
  is a subgroup count (see section 11), so the measured value is an upper estimate of the true mean; the
  charge 0.50 leaves a factor 1.5 for that and for the weighting by smoothness.

Mean y = 2^-59.828 per cell (relative SE 7.7%). Charged distinct relations:
V x (mean - 3 SE) x 0.95 x 0.98 = 2^51.111 = 2.29 D, with D the matrix dimension of section 6.
The surplus factor 2.29 is the margin of H-main.

Price per cell (`instr_per_cell` = 172.38; section 10): 32 per hit, with sum over ideals of 1/N(Q) =
5.099 hits per cell (exact sums over degree-one ideals above 5 <= q < 4096: 0.7533 and 1.1870 [E3]; above
4096 each side 1.05 x ln(ln B / ln 4096) by Mertens' theorem for prime ideals); tables 32 x 1.5 (2A+1)^6 per ideal
= 96 N/(2A+1)^6 per cell; 8 for initialisation and scan. Stage 2: 2^-20 V x 2 x 2^19 instructions.
Kept relations: at most 8 D of them are factored over the factor base by one residue test per ideal (128
instructions) and two exact norms: 8 D (128 N + 2^20) instructions. Schirokauer maps of the D relations used: 2^38
each (section 10).

## 6. Linear algebra

Relations are taken in order until their number equals the number of ideal columns occurring in them plus 60
Schirokauer columns plus 64 (at most D = N + 124 = 2^49.912 rows). The matrix over Z/r has these rows, the
occurring ideal columns, 60 Schirokauer columns and 64 zero pad columns. A kernel vector is found by Wiedemann's
algorithm (2D products for the scalar sequence, Berlekamp-Massey, D products for the solution), and its pad
coordinates are discarded. KNOWN = the ideals with a column; an ideal that never occurred has no logarithm.
The solution is checked on the unused surplus relations; a failed check is a failure (H-main).

Price: one row of one product costs at most 96 sparse entries x 24 + 60 dense four-limb products x 128 + one
reduction 512 + row overhead 64 + the dot product 128 = 10688 instructions; three passes and Berlekamp-Massey
(544 D^2): 32608 D^2 = 2^114.816 instructions (section 10).

## 7. Individual logarithms: independent trials on a re-randomised target

Both g and z need one (g being public does not make it free). For a target T != 1 in GT:

    repeat (trial):
        e <- uniform in [1, r);  T' = T^e                         # T' is uniform in GT \ {1}, whatever T is
        BOOT(T'): if it fails, next trial                          # section 7.1
        DESCENT of the prime ideals of the boot: if any node fails, next trial      # section 7.3
        return vlog(T) = vlog(T') / e

T' is uniform in GT \ {1} and independent of everything before it, and BOOT and DESCENT use only T', the table
KNOWN and fresh coins. So, given the main phase, the trials are independent and identically distributed, with the
SAME law for T = g, for T = z and for every challenge. Let q_b = P(boot succeeds), pi_d = P(descent succeeds |
boot succeeded) and C_d = E[work of the descent | boot succeeded]. These are three numbers attached to the
uniform distribution on GT, not properties of any ideal. The number of trials is geometric with parameter
q_b pi_d, and by Wald's identity the expected work of one logarithm is exactly

    (c_b + q_b C_d) / (q_b pi_d) = c_b / (q_b pi_d) + C_d / pi_d,                      (*)

c_b being the fixed price of one boot trial. There is no backtracking and no retry inside a trial: a node that
finds no usable relation ends the trial, and the next trial starts from a fresh uniform target. No statement of
the form "every reached ideal has a relation" is used anywhere.

### 7.1 Boot

T' = t0 + t1 x with t0, t1 in F_(p^6). LLL-reduce the lattice {(a, b) in Z[y]^2 : a + (t0/t1) b = 0 in F_(p^6)}
(dimension 12, determinant p^6, b weighted by 2^31; `lat.boot_reduced`, exact integer LLL) and take its first
vector. Then a - b x = c T' with c in F_(p^6), which does not change the logarithm (section 2). The trial succeeds
when (a, b) is admissible, N_g(a, b) < 2^1729 (the largest sampled norm plus 8 bits; larger norms are
rejected, so all arithmetic of a trial has a fixed size) and SPLIT (section 7.2) writes N_g(a, b) as a product of
primes of at most 130 bits. Sample: 300 boot lattices of UNIFORM ELEMENTS OF GT (`lat.gt_element`: w^((p^12-1)/r) for
uniform w in the tower representation, exactly the law of T'), each reduced ([E6] reruns one from its seed):
log2 N_g has mean 1709.3, maximum 1721.6; admissible fraction 0.303. Model:
q_model = (mean over the 300 draws of [admissible] x rho((n ln 2 + alpha_b)/(130 ln 2)), minus 2 SE) x
(1 - missing) x 0.8 = 2^-55.088, where 0.8 is the charged probability that SPLIT succeeds on a smooth norm
(H-ECM). CHARGED: q_b >= q_model/2 = 2^-56.088 (H-boot). Price of one trial: c_b = LLL 2^46 + exponentiation
2^25 + norms + trial division + SPLIT on a 1792-bit modulus = 2^57.758 instructions.

### 7.2 SPLIT: stage-1 ECM with checkpoints and an exact cost

SPLIT(N, bound 2^c): (i) remove the primes below 4096 by trial division; (ii) run a FIXED number of curves of
stage-1 ECM (Montgomery form, Suyama parametrisation, Montgomery ladder, `ecm.py`) with B1 = 2^(c/6) on the
full modulus N: for each prime q <= B1 in increasing order multiply the point by the largest power of q not
exceeding B1 and take d = gcd(Z, N). A point that is the identity modulo a prime l stays the identity, so d only
grows along the checkpoints; every increase d/d_prev is recorded. (iii) Refine the recorded factors and the
remaining cofactor to a pairwise coprime set of parts by gcds. (iv) A part m <= B is factored completely by trial
division up to 2^27 (at most 2^22.6 primes, below 2^34 instructions, charged once per descent); a part m > B must be
at most 2^c and pass 40 rounds of Miller-Rabin. SPLIT succeeds when every part passes; otherwise the boot
trial is dropped. There is no stage 2, no retry and no early abort: every call costs the same, and
incomplete factorisations only lower yields. A composite that passes Miller-Rabin (probability below 2^-80) or a
prime ideal of residue degree above one yields no root pair (s, t) in the next step and is dropped there.

For an admissible boot vector and a prime part l > B of N, section 3 gives one prime Q of Z[y] above l dividing the
relative norm and one root t: the node is the degree-one prime ideal (l, y - s, x - t), found by the roots of h
and of f or g modulo l (2^24 instructions).

Schedule for bound 2^c: B1 = 2^(c/6), curves = ceil(3/rho(6)) = 152,673. One curve costs at
most 20 B1 + 11 modular products (11 per ladder step; at most pi(B1) log2 B1 <= 1.8107 B1 ladder steps by
pi(x) <= 1.25506 x/ln x; measured 1.40 B1), pi(B1) checkpoint gcds, and 2^21 for the setup inversion; a product
modulo a k-limb N costs 24 k^2 + 64 k instructions and a gcd 2 bits (8k + 8) (section 10). For the boot:
B1 = 2^21.7, 2^25.99 products per curve.

Capture model (H-ECM): one curve captures a prime l (its point becomes the identity modulo l by the end of stage
1) with probability at least rho(ln l / ln B1). With ceil(3/rho(6)) curves and l <= 2^c a given prime is missed
with probability at most e^-3 = 0.05, and smaller primes much more rarely; SPLIT succeeds as soon as all primes
but one are captured at distinct checkpoints or by distinct curves, since the last cofactor is then a single
part. Measured ([E7], `ecm_table.json`): the capture frequency of single curves, ladder run modulo l itself, for
primes of 24 to 64 bits at u = ln l/ln B1 in {3, 4, 5, 6, 6.5}: in every cell with captures the frequency is at
least rho(u); the smallest ratios with good statistics are 1.53 +- 0.18 (48 bits, u = 3) and 1.97 +- 0.19
(48 bits, u = 4), and at u = 6 the ratio is between 2.7 and 5.7 (24 to 48 bits).
176 planted products of 1,649 to 1,680 bits (143 with prime factors of 12-24 bits, 30 with 16-32 bits, 3 with
22-48 bits; B1 from u = 3.5 at the largest prime, about 240 curves) were all split completely with every prime
isolated. Limit: primes above 64 bits were not measured; the production
bound uses primes up to 130 bits at the same u = 6, which is an extrapolation in the size of l at
fixed u (the ratio decreases slowly with l in the measured range and is expected to tend to 1 from above, the
group orders of Suyama curves being multiples of 12). The ledger takes no credit for ratios above 1 and charges
0.8 per call instead of 1 - (number of top-size primes) x 0.05. SPLIT is used only in the boot. Since it has a
fixed cost and its failures only drop the trial, H-ECM enters the claim through the single number q_b and
never through a cost or a termination argument.

### 7.3 Descent: two levels, every node sweeps one fixed ball completely, no retry inside a trial

A NODE is a degree-one prime ideal Q = (l, y - s, x - t) with l in (B, 2^130], on side f or g. Bins by
the size of l, top edges 56, 60, 66, 74, 86, 104, 130 bits. A node with l <= 2^74 is a LEAF: its relation
must lie entirely in KNOWN. A node with l > 2^74 is BIG: its relation may contain leaves (child bound
2^c <= 2^74). So the tree has two levels below the boot. Each (bin, side) has a fixed plan (c, radius R):

1. LLL-reduce L_Q = {(a, b) in Z^12 : a(s) = t b(s) mod l} (`lat.specialq_reduced`; determinant l).
2. SWEEP the whole ball {(a, b) = basis x v : v in [-R, R]^12, up to sign}: (2R+1)^12/2 cells, all of them, with
   the residue-matching sieve of section 5 over all degree-one ideals of norm up to 2^c of both sides (c = 54
   for a leaf); the divisibility condition is a linear form in v, split into two halves of six coordinates
   (table over one half, look-ups over the other): 1.5 (2R+1)^6 residues per ideal. Same thresholds and lemma:
   a survivor has both norms (Q removed) 2^c-smooth with every prime ideal identified. Survivors whose norm
   exceeds the plan's size limit (the largest sampled norm of the plan plus 8 bits; none of the sampled cells
   does) are skipped. The sweep has the same two caps as the main phase, both part of the fixed node price: at
   most 2^-20 of the cells go to the exact-norm stage and at most 2^16 survivors are examined, in sweep order.
   A cell is a USABLE RELATION when it is admissible, a survivor examined within these caps, and every prime
   ideal of norm <= B in it is in KNOWN. Its CHILDREN are its prime ideals of norm in (B, 2^c] (none for a leaf).
3. If the capped sweep outputs no usable relation the node FAILS and the trial ends (section 7); "empty ball"
   below always means this event, cap exhaustion included. Otherwise take the usable
   relation with the fewest children (first in sweep order among those), descend its children, and set vlog(Q)
   from the relation (both sides, Schirokauer maps included).

The descent uses no ECM, no cofactorisation and no sampling: every node is one complete sieve of one fixed ball
at a fixed price, which does not depend on how many relations the ball holds.

Plans (`model.ledger`, rows `nodes`; cheapest radius among 16, 24, 32, 48, ..., 192, 256 with lambda >= 2^10 and
lambda' >= 2^8):

| side | top bits | child bound c | R | log2 p worst / mean | log2 lambda | log2 lambda' | missing | children mean / max | log2 node | log2 subtree | nodes |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| f | 56 | 54 | 64 | -70.22 / -67.72 | 12.60 | 9.09 | 0.150 | 0.00 / 0 | 97.57 | 97.57 | 1.0 |
| g | 56 | 54 | 64 | -69.80 / -67.76 | 13.02 | 9.53 | 0.148 | 0.00 / 0 | 97.57 | 97.57 | 1.0 |
| f | 60 | 54 | 96 | -75.66 / -72.67 | 14.13 | 10.55 | 0.158 | 0.00 / 0 | 101.17 | 101.17 | 1.0 |
| g | 60 | 54 | 96 | -76.24 / -73.03 | 13.55 | 9.99 | 0.157 | 0.00 / 0 | 101.17 | 101.17 | 1.0 |
| f | 66 | 54 | 128 | -81.42 / -77.13 | 13.31 | 9.68 | 0.167 | 0.00 / 0 | 104.11 | 104.11 | 1.0 |
| g | 66 | 54 | 128 | -82.01 / -78.54 | 12.72 | 9.11 | 0.166 | 0.00 / 0 | 104.11 | 104.11 | 1.0 |
| f | 74 | 54 | 192 | -87.15 / -83.68 | 14.55 | 10.87 | 0.179 | 0.00 / 0 | 109.74 | 109.74 | 1.0 |
| g | 74 | 54 | 192 | -87.22 / -83.39 | 14.49 | 10.83 | 0.178 | 0.00 / 0 | 109.74 | 109.74 | 1.0 |
| f | 86 | 60 | 64 | -70.29 / -67.51 | 12.61 | 9.46 | 0.103 | 6.78 / 23 | 103.41 | 104.44 | 7.8 |
| g | 86 | 60 | 64 | -71.32 / -69.01 | 11.58 | 8.44 | 0.105 | 6.92 / 23 | 103.41 | 104.45 | 7.9 |
| f | 104 | 66 | 48 | -65.45 / -62.99 | 12.57 | 9.75 | 0.072 | 10.38 / 24 | 106.81 | 107.77 | 11.4 |
| g | 104 | 66 | 48 | -65.76 / -64.52 | 12.26 | 9.47 | 0.073 | 10.50 / 24 | 106.81 | 107.78 | 11.5 |
| f | 130 | 66 | 96 | -78.41 / -75.94 | 11.50 | 8.57 | 0.083 | 12.19 / 27 | 112.76 | 112.79 | 13.2 |
| g | 130 | 66 | 128 | -82.96 / -80.11 | 11.90 | 8.99 | 0.085 | 12.51 / 27 | 115.24 | 115.25 | 13.5 |

- Samples: for each bin and side 16 degree-one prime ideals drawn uniformly among all those above the primes
  with exactly the bin's top number of bits (`lat.draw_ideal`; the top edge is the least favourable size of the
  bin: the norms grow with l), LLL-reduced, 48 uniform ball cells each at R = 16 with exact norms and
  admissibility (`data.json`; [E5] regenerates ideals, lattices and cells from their seeds and compares). Other
  radii (all at least 16) use the measured scaling of the norm sizes: per doubling of R, +48.25 bits for N_f and
  +12.5 for N_g, the conservative side of the shifts measured at R = 64 and 256 ([E5]).
- log2 p: model probability that a ball cell is admissible with both norms 2^c-smooth (Q removed), rho rho with
  the alphas of section 5, for the WORST and the mean of the 16 sampled ideals.
- lambda = cells x p(worst) x (1 - missing) x 0.95: model-expected number of usable relations in the ball of the
  worst sampled ideal. Every plan has lambda >= 2^10 = 1024. missing: expected number of prime ideals of norm in
  (2^42, B] of a relation that are not in KNOWN (an ideal of norm l occurs in each of the D matrix rows with
  probability kappa/l, kappa the local density ratio of smooth numbers).
- lambda': the same with both norms of every sampled cell enlarged by 18 bits (ideal-to-ideal spread, see
  "Where eps comes from" below). Every plan has lambda' >= 2^8.
- children: model mean number of prime ideals of norm in (B, 2^c] of a relation (local density of the prime
  factors of smooth numbers, at the largest sampled norms), and after the slash the MAXIMUM possible number:
  each child has norm above B = 2^54 and the two norms are below the plan's size limits, so a relation has at
  most floor((limit_f + limit_g)/54) children. The algorithm takes the relation with the fewest children.
- node: the fixed price of the node (LLL 2^46 + sweep + factoring and Schirokauer maps of up to 2^16 survivors).
  subtree, nodes: model means of the work and of the number of nodes of the node with its children, each child
  charged at the worse side of its bin.

Number of nodes: a deterministic bound. A boot norm is below 2^1729 (larger ones are rejected, section
7.1) and every boot part exceeds 2^54; a big part exceeds 2^74 and has at most 27 children. With x big
parts and y leaves among the boot parts, 74 x + 54 y <= 1729, a descent visits at most
x (1 + 27) + y <= m = 644 nodes. No model enters this bound.

Work of the descent: a deterministic bound by the same count. A leaf costs at most 2^109.741 instructions and
a big node at most 2^115.241 (the largest node prices of the table), so a big part with all its children
costs at most 2^115.241 + 27 x 2^109.741, and maximising x (that) + y 2^109.741 over
74 x + 54 y <= 1729 gives, with 2^34 + 2^38 for the boot's own factorisation and maps,

    C_d <= 2^120.440 instructions for EVERY descent (a failed descent costs less).

For comparison the model mean (15.7 boot parts on average, local density of prime factors of smooth
numbers) is 2^118.280 instructions and 166.9 nodes; the ledger does not use it.

Success of the descent. A descent visits its nodes one after the other in a fixed (depth-first) order, at most
m = 644 of them, and fails exactly when a visited node has no usable relation.

H-desc: for every k <= m, among the trials that make a k-th visit, the fraction whose k-th visited node has
a ball without usable relation is at most eps = 1/1000. (An average over uniform targets, conditional only on
the EVENT that a k-th visit is made. Given a particular ideal and the table KNOWN a ball is deterministic, so
no statement conditional on the full history is made or needed, and ideals with an empty ball may exist.)

Let A_k be the event that none of the first min(k, V) visits met an empty ball (V the number of visits made) and
R_(k+1) the event that a (k+1)-th visit is made; R_(k+1) is contained in A_k. A_(k+1) is A_k minus the event
that the (k+1)-th visit is made and is empty, whose probability is at most eps P(R_(k+1)) <= eps P(A_k) by H-desc.
Hence P(A_(k+1)) >= (1 - eps) P(A_k), and since V <= m always, the descent succeeds on A_m:

    pi_d >= (1 - eps)^m = (1 - 1/1000)^644 = 0.5250 >= 1/2      (charged pi_d = 1/2).

No independence between nodes, no random exponent and no Markov step is used: m is a constant of the design.

Where eps comes from (section 11 gives the support). The premise allows one visited node in a thousand to have
an empty ball. H-indep: for an ideal Q visited by the descent, the events "cell is a usable relation" for cells
of its ball that are not unit multiples of one another are independent, each with probability at least
1/8 of its model value (Dickman at the cell's norm sizes with the measured alpha, times the factors of
the table). This is the independence heuristic on which every NFS yield estimate rests, applied inside one
ball, with a factor 8 of slack on the probabilities. Under H-indep, with lambda_Q the model mean of the
ball of Q and E[1/mult] >= 1/4 for unit multiples,

    P(ball of Q empty) = prod over non-associate cells (1 - p_i) <= exp(-sum p_i) <= exp(-lambda_Q/(8 x 4)).

How lambda_Q varies from ideal to ideal is measured at production size (`lambda_spread.json`, [E5]): for
4,200 further uniformly drawn degree-one ideals (300 per bin and side) the same model mean lambda_Q of the plan's ball
was estimated from 48 exact ball cells each. The smallest of the 4,200 estimates is 2^8.87 = 467, the
smallest 1% quantile of a bin is 2^10.07 and the smallest median 2^12.98 (0 of the estimates are below 2^8; the
estimates carry the sampling noise of 48 cells, which only widens the lower tail). So:

- for an ideal with lambda_Q >= 467, as all 4,200 sampled ones, H-indep gives
  P(empty) <= exp(-467/32) = 4.6e-07;
- the fraction of ideals with lambda_Q below every sampled value is not observed (0 of 4,200; a zero count
  in 4,200 draws excludes fractions above 7.1e-04 at the 95% level).

The CHARGED eps = 1/1000 covers the sum of the two, 7.1e-04, and is of the order that this direct evidence can
resolve. In addition every plan was required to keep lambda' >= 2^8 for the worst of its 16 plan ideals when
both norms of every cell are enlarged by 18 further bits (column lambda'): a sensitivity margin of 36 bits on
log2 N_f + log2 N_g, not used in the bound above.

By (*), one logarithm costs at most 2^57.758 / (2^-56.088 x 1/2) + 2^120.440 / 1/2 =
2^121.455 instructions in expectation; two logarithms 2^122.455.

Global meter. The two logarithms together run under a meter of 16 x 2^122.455 = 2^126.455 instructions;
reaching it is a failure. By Markov's inequality P(meter) <= 1/16 when the expectation is at most the charged
value (H-boot, H-desc, H-ECM), and metered work is at most unmetered work.


## 8. Success

Under H-main the main phase is deterministic up to Wiedemann's random vectors and the choice of a kernel vector
(a kernel vector with zero logarithm component, a failed Wiedemann run: below 2^-100 in total, re-drawn at most
twice). Past the gate the run therefore fails only if (i) those coins fail: allowance 0.0001; (ii) the global meter of section 7 is reached: at most 1/16
(H-boot, H-desc, H-ECM, Markov); (iii) the wrapper fails or vlog(g) = 0 mod r: below 2^-100, inside the same
allowance. A wrong x is caught
by the final check. So P(success) >= 1 - 1/16 - 0.0001 = 0.9374 for every z != 1, and z = 1 is
answered directly. With the gate: (27/64) x 0.9374 = 0.3954 >= 0.39.

## 9. Ledger (`model.ledger`, [E8])

| row | log2 instructions | kind |
|---|---:|---|
| relation collection (section 5) | 118.860 | fixed by the region and the stage-2 cap |
| factoring kept relations | 109.823 | cap (8 D relations) |
| Schirokauer maps | 87.912 | fixed |
| factor base | 72.774 | fixed |
| linear algebra (section 6) | 114.816 | fixed bound |
| two individual logarithms (section 7) | 122.455 | expected, formula (*) with charged q_b, pi_d and the bound C_d |
| wrapper (section 2) | 70.000 | bound |
| one run past the gate | 122.576 | expected |

One run: 2^122.576 instructions = 2^108.288 units. With the gate:
E[W] <= 1 + (27/64) x 2^108.288 = 2^107.043 < 2^107.05 units.

## 10. Instruction prices (scalar AArch64, counted, not measured)

All arithmetic is schoolbook on 64-bit limbs. Counts are of instructions of the classes ldr/str, add/adds/adcs,
sub/subs/sbcs, mul, umulh, cmp, csel, lsl/lsr, b.cond; each price is a count plus stated slack.

- Sieve hit, 32: bucket walk and key compare 6 (ldr entry, cmp, b.ne, ldr next, cmp, b); two-word cell index from
  the a-index and b-index 4 (mul, umulh, adds, adc); queue append 8 (ldr tail, str x3, add, str tail, cmp, b);
  owner pass 8 (ldr x2, ldr counter, add, str, add, cmp, b); slack 6.
- Residue, 32: Gray-code step of a(s) (index update 4, ldr step, add, cmp, csel/sub 3) 9; for look-ups the product
  by t and its reduction (mul, umulh, two-multiplication Barrett reduction 6, correction 3) 11; bucket address 2,
  insert or probe 6; slack 4. Ideals of residue degree d > 1 (norm q^d <= B, so q <= 2^27 and they are fewer
  than 2^24) use d-coordinate residues at d^2 times the price, inside the 7% allowance on N.
- Counter initialisation and threshold scan, 8 per cell (str; ldr x2, cmp x2, b x2, add).
- Exact norm, 2^19: about ten products in Z[y]/h (36 coefficient products of at most 4 x 4 limbs, 128 each) and
  a 6 x 6 fraction-free determinant (at most 110 products and 55 exact divisions of numbers of at most 8 limbs,
  below 1,500 each): below 2^18.3.
- (Relation, ideal) test, 128: a(s), b(s) by Horner (5 steps of mul, umulh, reduction 6: 40 each), t b(s) 8,
  compare 2, loads and loop 10; slack 28.
- Modular product on k limbs (Montgomery, coarsely integrated operand scanning), 24 k^2 + 64 k: two inner loops
  of k steps per outer step, each step ldr, mul, umulh, adds, adcs, adc, str, loop 2 = 9, so 18 k^2; per outer
  step ldr, mul, carries and loop 13; final conditional subtraction 4k + 6: 18 k^2 + 17 k + 6. The slack
  6 k^2 + 47 k per product also pays the at most ten k-limb modular additions (5k + 8 each) of one ladder step
  (xDBL = 2M + 2S + 1 product by a constant, xADD = 4M + 2S: 11 products per ladder step; `ecm.py` counts them).
- Binary gcd of n-bit numbers, 2n (8k + 8): each iteration removes at least one bit of the 2n operand bits and
  costs a k-limb shift (3k), compare (k) and subtract (3k) plus control.
- Trial division by the 564 primes below 4096: one k-limb by one-limb division (4k + 8) per prime, charged three
  times (valuations).
- 12-dimensional integer LLL, 2^46: by the potential argument at most log_(4/3) of the product of the eleven
  leading Gram determinants, below 66 x 574 bits for the boot basis (entries below 2^285), i.e. 1.1 x 10^5 swaps;
  at most 2 swaps + 12 passes of the main loop, each at most 576 operations on integers below 7,000 bits
  (110 limbs) at 24 x 110^2 + 64 x 110: 2^45.1. Special-q lattices have smaller entries.
- Schirokauer maps of one relation, 2^38: one exponentiation in (Z/r^2)[y,x]/(h, f) with an exponent below r^48
  (12,192 bits, at most 2 products per bit), a product being 48^2 coefficient products modulo r^2 (8 limbs,
  24 x 64 + 64 x 8 = 2,048) plus the reductions modulo h and f (below 20%): 2^37.1; the g side is below 2^31.
- Linear algebra row, 10688: sparse entry 24 (ldr index, address, 4 ldr, adds/adcs x4, adc, loop 2 = 13 for a
  coefficient +-1; larger valuations repeat the addition; slack to 24), at most 96 entries on average: a norm
  below 2^516 has at most 43 prime ideals of norm at least 2^12 per side and the sieved ideals below 2^12
  contribute under 2 per cell; dense product 128 (16 x (mul, umulh, adds, adcs, adc) = 80, loads 8, loop 8; slack
  32) for each of the 60 Schirokauer columns into a nine-limb accumulator; one Barrett reduction of the row
  accumulator 512 (two 5 x 4-limb products and corrections, below 300); overhead 64; dot product for the scalar
  sequence 128. Berlekamp-Massey on a sequence of length 2D: at most 2D^2 + 4D modular products of 256 + 16: 544 D^2.
- Root data of one prime, 2^24: roots of h and of f(s, .), g(s, .) modulo q <= 2^54 by gcd with x^q - x: fewer
  than 2^8 polynomial products of degree below 8 per root, 64 coefficient products of 32 instructions each.
- Exponentiation in F_(p^12), 2^25: at most 508 products of 2^15 (144 four-limb Montgomery products plus additions).

## 11. Heuristics, support, limits

- H-main. A statement about the one fixed field and the one fixed cell region, not a probability: (a) the
  region contains at least D distinct relations and at most 2^-20 V stage-1 survivors; (b) at most 2.2 B/ln B
  factor-base ideals exist; r divides neither class number; the Schirokauer maps are injective on units modulo
  r-th powers; the linear system has the virtual logarithms as its solutions up to functionals vanishing on
  actual elements. Support for (a): the standard NFS yield model (norms of admissible cells are B-smooth with
  the probability of random integers of their size shifted by Murphy's alpha, independently on the two sides
  given the sizes) evaluated on exact production-size norms, alpha measured on the admissible cells at + 3 SE,
  and a 2.29-fold surplus of charged relations over the dimension: the model would have to
  overestimate the yield by that factor for (a) to fail. (b) is standard for NFS discrete logarithms and not
  exercised here at any size. Because the main phase has no coins except Wiedemann's (failure below 2^-100,
  section 8), H-main is a premise of the whole claim rather than a failure probability: if it is false the
  attack fails for every target. Limit: no production-size relation is exhibited.
- H-boot. q_b >= q_model/2: the fraction of uniform elements of GT whose reduced boot vector is admissible with
  a 2^130-smooth norm that SPLIT handles is at least half the Dickman value computed on the sampled
  norms. The 300 sampled lattices are those of uniform elements of GT, the exact law of the trial targets. This
  is an average over uniform targets, the heuristic of Guillevic (Math. Comp. 88, 2019, section 6.3: the norms
  of the reduced preimages behave like random integers of the same size), not a statement about any target.
- H-desc. Emptiness rate: among trials making a k-th visit, at most a fraction eps = 1/1000 meet an empty ball
  there (section 7.3). This is the ONLY premise on the descent: the number of nodes and the work of a descent
  are bounded deterministically (section 7.3), and no model of the tree, of the children or of their sizes
  enters the ledger. Support, in three layers. (i) Derivation: under H-indep (section 7.3) an ideal with
  lambda_Q >= 467 has an empty ball with probability at most 4.6e-07.
  (ii) Production-size measurements behind the constants of that derivation: the plans' lambda and lambda' are
  computed from exact norms of 16 ideals x 48 cells per bin and side; the ideal-to-ideal variation of lambda_Q
  from 4,200 further ideals x 48 cells (`lambda_spread.json`: none below 467); unit multiples inside
  production balls from 336 cells of those lattices: the number of units y^i e2^j e3^k keeping a ball cell
  inside its ball has mean at most 1.46 per bin (maximum 8) with exponents up to 6, and the count
  with exponents up to 4 differs in 0 of the 336 cells (units with larger exponents have coefficients
  too large to keep a cell inside). The ledger charges lambda/4 for the expected number of
  non-associate relations, more than twice the largest measured bin mean. Limits of this measurement: y, e2, e3
  are the fundamental units reported for K_h by the dependency (computed there with PARI/GP; we verify only
  that they are units), so the table counts multiples by the group they generate; and the average is over
  uniform ball cells, not weighted by smoothness.
  (iii) Realized counts, measured where they can be measured: a toy tower with the same shape and the same
  node predicate, [E9]:
  Toy: the same tower shape with a quadratic base field (h = y^2 - y - 1, f = x^4 + y x^2 + y^2,
  g = x^2 - 62 y, p = 3907) and the SAME node predicate: degree-one prime ideal of norm l in (B, 8B], B = 2^16,
  LLL-reduced lattice, complete sweep of the ball [-R, R]^4, relation = admissible cell with both norms
  B-smooth after removing l (child bound below l), counted exactly (`ballcount.py`; unit multiples counted).
  "random" = uniformly random degree-one ideals; "reached" = the large prime of an otherwise smooth cell of a
  swept parent ball, i.e. the ideals a descent actually visits. 44,735 completely swept balls in 17
  rows (`ballcount.json`): cells per ball, balls, empty balls, mean count, D = variance/mean of the counts,
  observed empty fraction, exp(-mean/D), realized total / plain Dickman prediction.

  | ideals | side | R | cells | balls | empty | mean | D | empty fraction | exp(-mean/D) | realized/Dickman |
  |---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
  | random | f | 3 | 1200 | 3051 | 288 | 7.7 | 4.17 | 0.0944 | 0.1567 | 0.35 |
  | random | f | 4 | 3280 | 3051 | 88 | 12.4 | 4.48 | 0.0288 | 0.0632 | 0.32 |
  | random | f | 5 | 7320 | 3051 | 13 | 19.9 | 4.99 | 0.0043 | 0.0185 | 0.31 |
  | random | f | 6 | 14280 | 6893 | 13 | 26.4 | 5.74 | 0.0019 | 0.0100 | 0.29 |
  | random | f | 7 | 25312 | 3356 | 0 | 36.0 | 6.44 | 0.0000 | 0.0038 | 0.28 |
  | random | f | 8 | 41760 | 3356 | 0 | 44.8 | 6.87 | 0.0000 | 0.0015 | 0.27 |
  | random | g | 4 | 3280 | 3051 | 604 | 5.6 | 4.97 | 0.1980 | 0.3261 | 0.32 |
  | random | g | 5 | 7320 | 3050 | 271 | 8.8 | 5.18 | 0.0889 | 0.1820 | 0.31 |
  | random | g | 6 | 14280 | 6893 | 356 | 11.5 | 5.58 | 0.0516 | 0.1274 | 0.29 |
  | random | g | 7 | 25312 | 2288 | 63 | 15.2 | 5.78 | 0.0275 | 0.0721 | 0.27 |
  | random | g | 8 | 41760 | 1832 | 17 | 19.3 | 6.28 | 0.0093 | 0.0462 | 0.26 |
  | random | g | 9 | 65160 | 1221 | 6 | 23.2 | 6.73 | 0.0049 | 0.0318 | 0.25 |
  | random | g | 10 | 97240 | 916 | 3 | 27.4 | 7.27 | 0.0033 | 0.0231 | 0.24 |
  | reached | f | 6 | 14280 | 1524 | 2 | 28.7 | 5.45 | 0.0013 | 0.0051 | 0.29 |
  | reached | f | 8 | 41760 | 245 | 0 | 49.6 | 6.56 | 0.0000 | 0.0005 | 0.27 |
  | reached | g | 6 | 14280 | 712 | 19 | 12.8 | 5.31 | 0.0267 | 0.0889 | 0.28 |
  | reached | g | 9 | 65160 | 245 | 2 | 25.4 | 7.89 | 0.0082 | 0.0403 | 0.25 |

  Reading. Counts are NOT Poisson: D is between 4 and 8 and empty balls DO occur when the mean is small.
  Two known sources inflate D: unit multiples of a relation inside the same ball (the toy's unit group is
  large relative to its small balls) and the pooling of ideals with different means; we did not separate them. In every row the observed empty fraction is at most exp(-mean/D), the
  value for independent clusters of the measured size (within three standard errors; [E9] asserts it), and it
  falls as the mean grows, on both sides and for both populations: from 9.4% at mean 7.7 to 0.19% at mean 26
  on side f, and no empty ball at all among the 6,957 balls of the three rows with mean above 35; reached
  ideals behave like random ones. The largest toy mean is 50; production plans have model means above 1,000. The realized totals are
  0.24 to 0.35 of the plain Dickman prediction at these tiny bounds (no alpha is applied in the toy, and
  admissible norms carry no small primes), inside the allowance 8 of H-indep.
  Limits: no production ball was swept, so realized counts are measured only at toy size, where smoothness
  bounds are tiny; H-indep itself is a heuristic (it excludes, for instance, a hidden mechanism that would make
  a fixed fraction of ideals have no relation at all; the toy shows no such mechanism: its empty fractions fall
  with the mean count as fast as the independent-cluster law predicts); the sampled ideals are uniform
  degree-one ideals at the top edge of their bin, not ideals drawn from production parent relations (the
  conditioning of a reached ideal is on the smoothness of the rest of its parent norm, which concerns other
  primes; the toy's "reached" population tests this at toy size); the model for missing ideals (inside
  lambda) is the local-density heuristic for smooth numbers.
  If H-desc or H-boot fails the effect is on expected work and on the meter, never on correctness.
- H-ECM. Capture probability of one stage-1 curve for a prime l is at least rho(ln l/ln B1), and SPLIT succeeds
  on a smooth norm with probability at least 0.8. Support and limit: section 7.2 and [E7].

What the experiment does and does not show. [E1]-[E8] check exact arithmetic at production parameters,
regenerate stored samples from their seeds, re-measure ECM captures at small prime sizes, and recompute the
ledger, the manifest score and the numbers of this text. [E9] sweeps fresh toy balls and audits the stored toy
campaign. Nothing is evidence for production cost beyond its stated scope: no production-size relation, boot or
descent was run, and no end-to-end run of the algorithm at any size is included.

## 12. Memory, preprocessing, advice

Peak storage: sieve counters (two 16-bit words per cell, 2^113.4 bytes) plus hit queues, tables, relations,
matrix (below 2^62 bytes) and the node sweeps (smaller): below 2^114 bytes. Storage beyond one 64-bit address
space is addressed by (chunk, offset) pairs, which is the two-word index priced in every hit. The rules give
memory no scalar contribution. All preprocessing is for this one target and is charged in full. Advice: none;
h, f, g and the plans are design choices and the public constants of the profile.

## 13. Context and credits

Method: special tower number field sieve (Kim-Barbulescu, Crypto 2016) with the twisted pair (f, g) for
gcd(6, 2) != 1 in the manner of Guillevic-Singh; boot by Guillevic's subfield lattice (Math. Comp. 2019);
Schirokauer maps; Wiedemann; Murphy's alpha; Lenstra's ECM in Montgomery form. Published estimates for this
field (about 2^100 to 2^103 for the main phase) are context only; no published total is a premise. The tower,
the polynomial pair and the matching sieve with its exact threshold are taken from submission
f8a14ba8-9a8a-4c47-91f8-29b129f2dbd6 (rejected for its descent, ECM and FFT premises, which this package
replaces); the run gate and the expected-work convention follow accepted entries of this benchmark. All samples,
the model and the ECM measurements are new.
