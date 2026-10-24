# Argument — v12b canonicalising Frobenius-class walk (bls12-381-fixed-g1-invert-g2)

Every number below is either an exact integer, an exact `Fraction`, or a field of
one of the package's raw JSON files. Nothing is retyped from a previous revision:
`verify_certificate.py` regenerates `certificate.json` from `constants.py`,
`ledger.py`, `a-*.json`, `d-*.json`, `e-equivariance.json`,
`modmul-count.json`, `routing-out.json` and `part_a-out.json`, and fails
on any moved number.

## 1. Reduction

P is the fixed nonidentity G1 point of BLS12-381 supplied by the challenge; z is
uniform in the order-r subgroup of GT. With g = e(P, Q0), non-degeneracy of the
pinned pairing makes g a generator of the order-r subgroup, so z = g^a for a
unique a mod r and Q = [a]Q0 satisfies e(P,Q) = z by bilinearity. `part_a.py`
verifies the tower, ord_r(p) = 12, N = (r-1)/12, the selector identities
(g1 = -g0, g0^6 = -1/8, zeta = 2 g0^2), and a planted instance through a real
pairing evaluation including the negative control: 24/24 PASS.

## 2. The walk

State: an element y of the order-r subgroup of GT, held at the canonical
Frobenius representative of its class, plus labels (a,b) with y = g^a z^b.

Class key: `kappa(y)` = the twelve Frobenius keys rotated so the minimum is
first; k = index of that minimum, so sigma^k(y) is the canonical representative.
The key is Frobenius-invariant by construction (it is a min-rotation over the
orbit), which audit A confirms at 528/528 per fixture.

Step, decided by the class key alone (`routing.fields`): the step type word picks
a dense Fp12 multiply on five of eight classes and a Granger-Scott cyclotomic
squaring on three of eight; the table word picks one of 2^14 rows and the column
word picks one of 12 columns. Then

    mul:  raw = y * tbl[t][col]           delta = u[t][col],  (a,b) -> (a+delta, b)
    sqr:  raw = cyclo_sqr_fast(y)                             (a,b) -> (2a, 2b)
    both: k = kappa(raw);  y <- sigma^k(raw);  (a,b) <- (a p^k, b p^k) mod r

This is v12b's change. The previous revision applied sigma^(-k)(m) at index
(tw-k) mod 12 to a non-canonical hold. That is a class function, and it is what
created the rejected two-cycle: p^6 = -1 (mod r) makes sigma^6 inversion, hence
sigma^(-(k+6))(m) = sigma^(-k)(m^-1); as k ranges over Z/12 the applied
multiplier set is inversion-closed for every table, so m followed by m^-1 returns
the element and the label sum delta_i + delta_{i+6} = 0 vanishes. v12 measured
this at 384/384 ladder pairs per fixture and showed that restricting the routing
modulus does not remove it, because k still spans Z/12.

In v12b the applied multiplier is the table entry itself, and the table is built
from independent exponents 1 <= u <= (r-1)/2 (`ClassWalk._build_table`). The
inverse of any entry has exponent r-u in [(r+1)/2, r-1]; the code checks
`u_neg_u_overlap == 0` on the actual list, so the (m, m^-1) family is
unreachable by a range identity rather than by a probability claim. The remaining
exact two-cycles are the O(1) states with y^3 = 1 and the single y = m^-1 state
per configuration, in a 2^251-node space.

## 3. Equivariance

The class group is mu_12 = <sigma>; the generators and their orders are
enumerated in `equivariance.py` and the orders are verified at run time
(sigma:12, sigma^2:6, sigma^3:4, sigma^4:3, sigma^6:2), with p^6 = -1 checked on
every fixture. The reviewer's law -- class(S + R_c) == class(sigma S + R_c) for
every sigma in the class group -- is tested in the form this design has to
satisfy: E1 the successor class through the shipped path is identical from all
twelve members of a class (0 splits required), and E2 no sigma^j of an applied
multiplier is again a table entry (0 hits required; this is the (m,m^-1) test,
and sigma^6 is exactly the inversion that made the old design fail). E3 is the
control: the same two tests on v10f's sigma-ladder and held-offset rule, where a
nonzero count proves the test detects the defect. Both control counts are nonzero
on every fixture; the shipped object reads 0/0.

## 4. Ledger

Prices: C = 20,000 (organiser reference dense Fp12 multiply), FP_MUL = 383 (blst
base-field multiply), FP_ADD = 20, CMP6 = 24, HASH_W = 24.

Fold, nine rows summing to 1813 = 403 + 766 + 180 + 264 + 72 + 8 + 32 + 32 + 56:

  * selector_ab 403 = FP_MUL + FP_ADD: a-b, then g0*(a-b), because gamma = g0*(1-u)
    exactly, so d = gamma0 a + gamma1 b collapses to one base-field multiply;
  * selector_zeta 766 = 2 FP_MUL: zeta*a and zeta*d; the zeta^2 forms are one
    subtraction from the zeta forms;
  * fold_adds 180 = 9 FP_ADD; fold_cmp 264 = 11 CMP6 (six |.| tests and a
    five-way maximum with cmov); hashing 72 = 3 HASH_W; dp_test 8; key traffic
    32; table/column/loop 32; witness add 56.

Cyclotomic squaring: 24 FP_MUL + 122 FP_ADD = 11,632, the exact operation count of
`tower.cyclo_sqr_fast`.

Renormalisation: RENORM = 286 = 2 x 143. One label's multiply is
`mulmod_r.c`, a 4-limb CIOS Montgomery multiply that computes
REDC(a * (p^k 2^256 mod r)) = a p^k mod r, with the twelve p^k pre-scaled once at
table build (auxiliary, charged). Compiled with `clang -O2 -arch arm64` the body is
143 instructions; the differential harness runs 255 cases, all twelve k against
random labels plus a in {0,1,2,r-2,r-1}, 0 mismatches. Both step types pay it.

So mul_step = 20,000 + 1,813 + 286 = 22,099 and sqr_step = 11,632 + 1,813 + 286 =
13,731. With f = 3/8, c_raw = (5/8)(22099/20000) + (3/8)(13731/20000) =
94805/100,000. W_stress reprices only the fold at 13/10 and leaves the counted
rows alone: c_stress = (5/8)((20000 + 23569/10 + 286)/20000) +
(3/8)((11632 + 23569/10 + 286)/20000).

S is the least integer with S^2 (1-f^2) >= 2 mu N ln(1/(1-T)), mu = 11/10,
T = 392/1000, N = (r-1)/12; log2 S = 125.810636. Then

    W_raw = S c_raw = 2^125.733671;  W_stress = S c_stress = 2^125.774473.

Claim 12578/100 = 125.78 is 0.005527 over W_stress and 0.040802 over W_raw; the
exact test 2^12578 >= W_stress^100 passes as integers
(`verify_certificate.claim_ge_W_stress_exact_integer`).

Auxiliary rows, all exact and all below 2^-60 of W_stress: table build
2^14 x 12 = 2^17.585 reference operations; trail starts 4; report writes
2^38.693 (S/theta expected reports x 12 stores of 8 bytes); the exponent solve
2^9 (it is now one cancellation x = (a1-a2)/(b2-b1) mod r, charged at the old and
larger value); verification 2^9.

Storage: 96 bytes per record (32-byte class key + two 32-byte limbs). Reports are
capped at (11/10) S theta + 2^25 with refusal on overflow, the Chernoff tail is
debited at 1e-20, and peak memory is below 2^53 bytes. Memory is unscored and
stated.

## 5. Success

Four trails of w = max(S/4, 8 theta) = 2^123.810636 steps, one shared DP store,
no restart at a report. Orbit coalescence gives

    Lambda = 4 C(w,2)/N + C(4,2) w^2/N = 8 w^2/N = S^2/(2N) = 0.6369029082,

the DP density enters as detection on the shared suffix, E[L] = w/3 = 2^122.226
so 1 - e^-(E[L]/theta) = 1 - 2^-42.2 at theta = 2^-80, and

    1 - e^-Lambda >= Lambda - Lambda^2/2 = 0.4340802,
    P(success) = 0.4710720 at w/sqrt(N) = 0.2821575.

The floor 0.39 is a declared 17.2 per cent reserve under that model, fixed by
COST_MODEL.md line 45.

Loss rows: cutoff waste 5e-7, cap overflow 1e-20, fingerprint false confirm
1e-25, degenerate denominator 1e-70, column bias 1e-6 (12 x 2^-32), undisclosed
random-map reserve 297e-6; total 2.985e-4, absorbed by the reserve and not by the
claim.

## 6. Toy evidence and its scope

`toy12d.py` executes the shipped rule: canonical state, renormalisation on every
step, inverse-free lower-half table, decision-only cache (no cached field element
is ever re-Frobeniused), four trails and one store per group, verified recovery
(every merged x is checked by exponentiation against the planted secret), and an
orbit census that flags any exact state repeat shorter than 8 steps. The fast
sigma^k is a counted object: 288/288 agreement with the literal tower Frobenius.

Audits per fixture: class-id invariance 528/528; successor class from all twelve
members through the shipped path 528/528; inverse of no applied multiplier in the
table 528/528; the v10f control splits on every fixture. Table checks:
`table_inverse_free` and `table_all_lower_half` true, zero kappa ties.

Scope, stated rather than glossed: at r ~ 2^15 a lap using one row's complete index
set returns element and labels exactly, because p has order 12 and p != 1 gives
sum_{i<12} p^i = 0 (mod r). That family is negligible at production and dominant at
toy size, so toy cycle lengths are fixture-limited and the production (m,m^-1)
statement is the range identity, not a measured rate. The toy is also not at
production's dimensionless point (w/sqrt(N) 2.286 against 0.2822), so what
transports is the event chain and the scaling variable, not the absolute level.
The 2^-5 cell is disclosed as the regime boundary of the closed form.

## 7. Disclosed hypotheses

H1 the class walk behaves as a random map to within the 10 per cent inefficiency
allowance already carried as mu = 11/10; H2 the routing fields derived from the
class key are pseudorandom (120,000-key histogram: fields bit-disjoint, no
collapse, table chi-square consistent with uniform, no checkpoint field, no
window); H3 the algorithm's coins are ideal. No production run, no unconditional
theorem, no novelty claim and no practical break of BLS12-381 is claimed.
