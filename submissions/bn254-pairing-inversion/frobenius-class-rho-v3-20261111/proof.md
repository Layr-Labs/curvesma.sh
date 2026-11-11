# BN254 pairing inversion via Frobenius-class Pollard rho with a stack cycle detector

This is a generic-group-style rho for the protected 254-bit prime-order GT,
improving the organizer's reviewed small-memory baseline
`submissions/bn254-pairing-inversion/rho-baseline-v2-20261110` (score 135.216266). It is
not a best-known finite-field attack, not a practical break and not a
security-bit estimate. It seeks the exploratory `plausible_not_refuted`
standard. Its work bound is a deterministic every-run cap; its success
probability rests on the disclosed random-mapping heuristic, as in the
baseline.

**What is reused (credited).** From the baseline: the scalar machine model and
its primitive bounds (Montgomery word arithmetic, the 2^20 Fp12 product, the
2^16 SHA-256 step, the 4096 table step), `arithmetic.c`, the bounded scalar
sampler, the degeneracy-conditioning argument, and the pairing output
procedure. All are restated in `machine.md`/below so the submission is
self-contained; no prior acceptance is assumed. **What is new:** (1) the walk
runs on classes under the order-12 Frobenius group (Duursma-Gaudry-Morain
1999; Gallant-Lambert-Vanstone 2000; Wiener-Zuccherato 1998), (2) Floyd's
three evaluations per iteration are replaced by Nivasch's stack detector
(one evaluation per step), (3) fruitless cycles are handled by exact
detection and a squaring escape (cf. Bos-Kleinjung-Lenstra 2010,
Bernstein-Lange-Schwabe 2011), (4) the step cap is set for a predicted
success of about 0.64 instead of 0.86.

## 1. Group facts used (all checked in experiment.py on the real profile)

phi(y)=y^p is a field automorphism of Fp12 of order 12 and maps GT to GT.
On GT (order r), phi^k(y) = y^(lambda_k) with lambda_k = p^k mod r. p has
order exactly 12 modulo r (p^12=1, p^k!=1 for k=1,2,3,4,6), so the lambda_k are
the 12 distinct 12th roots of unity mod r and lambda_6 = r-1: inversion is
phi^6 and in the protected basis is M_6 = diag((-1)^i). If phi^d(Y)=Y with
0<d<12 then Y^(p^d-1)=1, and r does not divide p^d-1, so Y=1. Hence every
non-identity class {phi^k(Y)} has exactly 12 elements and the identity is a
singleton: n = (r-1)/12 + 1 classes.

The representative rep(Y) is the minimum of the 12 conjugates under the fixed
total order of machine.md (lexicographic on canonical Montgomery coefficients,
priority 1,7,0,2,...). It is a function of the class. The cheap key uses only
coefficients 1 and 7 (rows 1,7 of M_k read only z1,z7); ties on the key are
resolved by full comparison (budgeted path), so rep is exact.

## 2. Complete algorithm (one trial; `rho.py:class_rho` is the reference)

Read only the profile and z; validate z!=0, z^r=1; if z=1 output x=0 (pairing
track: the permitted infinity output). Let R_T=1024 and

    E = isqrt(floor((9/5)^2 * n * 1024/1023)) + 1   (about 2^125.854566).

Sample a 256-bit salt and 2*1025 scalars by bounded rejection. Build
M_j = g^(c_j) z^(d_j), j<1024. Start: Y = rep(g^a0 z^b0) = (g^a0 z^b0)^lambda_k,
labels (a,b) = lambda_k(a0,b0). Invariant: Y = g^a z^b, Y canonical. Push
(Y,a,b) on an empty stack. Then for steps 1..E:

* normal step: j = SHA256(salt||encode(Y)) index mod 1024;
  Z = Y*M_j, (a,b) <- (a+c_j, b+d_j);
  escape step (flag set): Z = Y^2, (a,b) <- (2a, 2b);
* Y = rep(Z) = Z^(lambda_k), (a,b) <- lambda_k(a,b) mod r;
* pop while top > Y; if top = (Y,A,B):
  if B != b output x = (a-A)/(B-b) mod r after verifying g^x=z;
  else (then also a=A) set the escape flag (the next step squares this Y);
  stack unchanged;
* otherwise push (Y,a,b); abort if depth > 4096.

Abort (failure) also if more than 1024 tie events occur or E steps pass.
Verification failure is reported as failure, never as an answer. There is no
restart, amplification or amortisation.

**Correctness.** Labels are updated exactly as the exponents: multiplication
adds (c_j,d_j), squaring doubles, conjugation by phi^k multiplies by lambda_k
(valid because Y has order dividing r). Equal canonical elements with labels
(a,b),(A,B) give a+xb = A+xB mod r; with B!=b the unique solution is
x=(a-A)/(B-b); B=b forces a=A since g has order r. Every output is verified.

## 3. Stack detector (Nivasch 2004)

Lemma 1. Let y_0,y_1,... be produced by a deterministic map on a finite set,
tail mu, cycle lambda, and let m=y_(mu+J) be the cycle minimum (0<=J<lambda).
The stack procedure finds its first equality exactly at index mu+J+lambda,
and the equal stack entry is the first occurrence of m.
Proof. The stack is strictly increasing bottom to top. An equality needs a
repeated value, i.e. index >= mu+lambda. A cycle element y_(mu+i), i<J, pushed
before m, is popped when m arrives (m is smaller); its repeat at mu+i+lambda
finds it absent. m itself is never popped after its push (nothing smaller
occurs on the cycle), and on its repeat everything above it is larger and is
popped, so the top equals m. Elements y_(mu+i), i>J, repeat only after
mu+J+lambda. QED. Hence detection time T = mu + lambda + J <= mu + 2*lambda,
and the two equal entries are exactly one loop apart.

## 4. Fruitless cycles and the escape

Along a path, labels transform affinely: a loop of the class walk maps
(a,b) -> (mu a + sum e_j c_j, mu b + sum e_j d_j), with mu a product of
lambda's (and a power of 2 if escapes occur) and integer coefficients e_j
fixed by the path. A loop is *fruitless* if mu=1 and all e_j=0 mod r: then it
closes for every x and gives identical labels. Because of the inversion
(lambda_6=-1), Y -> rep(Y M_j)=(Y M_j)^(-1) -> rep(Y^(-1)) = Y is a fruitless
2-cycle, expected about once per 12*1024 steps (measured below); longer
fruitless cycles exist too. Every index occurring once in a loop contributes
a unit to its e_j, so a fruitless loop of length t uses each of its indices
at least twice.

Define the ideal map F'(Y) = rep(Y^2) if Y is the minimum of a fruitless
cycle of the class walk f, and F'(Y) = f(Y) otherwise. Fruitlessness is a
formal property of the path (indices and lambda choices), so F' is a fixed
deterministic map on classes.

Lemma 2. Until the algorithm stops, its element sequence is the F'-sequence
with, after each arrival at a fruitless minimum m_C that is not on the stack
at that moment, one *excursion*
(the remaining |C|-1 elements of C and m_C again) inserted; after each
excursion the stack equals the stack of the stack procedure run on the
F'-sequence. Proof. At an arrival at m_C that is not already on the stack,
m_C is pushed. The excursion elements are larger than m_C, so they never
equal entries below it, and they are distinct from each other (one loop of a
cycle). By Lemma 1 the only equality is the return to m_C, whose labels are
identical (fruitless): the algorithm escapes, leaving the stack as after
pushing m_C, and produces rep(m_C^2) = F'(m_C). If m_C is already on the
stack on arrival, the path since its earlier occurrence contains the escape
F'(m_C): that is a genuine F'-detection (below). QED.

Lemma 3 (no degenerate genuine detection except with probability exactly
1/r). Fix x, salt, s_j = c_j + x d_j and s_0 = a0 + x b0. The F'-sequence,
all hash indices, canonical choices and escape points depend only on these,
while b0 and all d_j remain independent and uniform (for each d_j exactly one
c_j realises s_j; bounded rejection sampling is uniform conditional on
success). Every label b is then a linear form in (b0, d_j) with coefficients
fixed by the walk; b0 has a nonzero coefficient (a product of lambda's and
2's). At the first F'-detection, B-b = (mu-1)*b_prev + sum e_j d_j for the
one loop between the equal entries (Lemma 1). If the loop contains an escape,
mu = 2^s * zeta with 1<=s<=E and zeta a 12th root of unity; mu=1 would need
ord_r(2) | 12s, impossible because ord_r(2) ~ 2^247.3 > 12E (computed from the
factorisation of r-1; even without assuming its 94-bit factor prime, the
computed order divided by that factor, ~2^153.8, is a proven lower bound and
exceeds 12E ~ 2^129.4; experiment.py asserts this). Otherwise the loop is an f-cycle not escaped
under F', hence non-fruitless. Either way B-b is a nonzero linear form, so it
is uniform and B=b has probability exactly 1/r. QED.

The excursion steps are the only *waste* W. Each escape costs its cycle
length (2 for the dominant 2-cycles).

## 5. Success probability (heuristic, disclosed)

Success holds if (i) the F'-detection time T' satisfies T' + W <= E, (ii) no
stack or tie abort occurs, (iii) the detection is nondegenerate (1-1/r).
Hence P(success) >= P(T' <= (1-delta)E) - P(W > delta E) - P(abort) - 1/r with
delta = 1/64.

H1 (random mapping on classes, the baseline's heuristic transferred to
classes). F' behaves like a random map on n classes started at a uniform
class, with the r-adding correction n_eff = n*R_T/(R_T-1); and the order
position J of the cycle minimum is uniform on the cycle. Under a random map,
rho = mu+lambda has P(rho > k) ~ exp(-k^2/(2n)), lambda is uniform on
1..rho, and T = rho(1 + v w) with v,w uniform on (0,1). With
c = E(1-delta)/sqrt(n_eff) = (9/5)(63/64):

    P(T <= c sqrt(n_eff)) = int u e^(-u^2/2) Q(c/u) du,
    Q(y) = 0 (y<=1), (y-1)(1-ln(y-1)) (1<y<2), 1 (y>=2),

evaluated numerically by experiment.py: 0.6448 (0.6558 at c=9/5). The mean
detection time is 1.25*sqrt(pi/2)*sqrt(n_eff).

H2 (waste). The waste fraction stays below delta = 1/64. Fruitless 2-cycles
occur at rate 1/(12 R_T) = 8.1e-5 per step and cost 2 steps (expected waste
~1.6e-4 per step, about 1/100 of delta). A fruitless loop of length t needs
every index to repeat; for t << R_T this has probability at most
(e t/(2 R_T))^(t/2), and for long loops it requires 1024 independent root-of-
unity sums to vanish simultaneously. Toy runs below measure the waste.

H3 (aborts). The stack holds the right-to-left minima of the sequence; under
H1 its depth is a sum of independent Bernoulli(1/i) with mean ~ln E ~ 88, so
P(depth > 4096) is astronomically small; tie events are expected ~2^-121.7
times per run against a cap of 1024.

Prediction: P(success) ~ 0.645 - negligible. **Claimed success >= 0.50**,
the protected minimum is 0.39. This is an exploratory heuristic, not a
production theorem; the work cap holds even if H1-H3 fail.

**Toy evidence (experiment.py, seed 20261111, 512 trials per row).** Prime
additive groups Z/q with q = 1 mod 12, generator 1, and the automorphism
group y -> lambda^k y of order 12: the same abstract structure as GT with its
Frobenius group, the same class_rho code, SHA-256 on a 384-byte encoding,
1024-entry table, stack detector and squaring escape; canonical choice by
integer minimum. Results (development-evidence.json):

| q | classes n | success at E=9/5*sqrt(n_eff) (pred 0.656) | mean T / predicted | 2-cycles per step (pred 8.1e-5) | waste fraction |
|---:|---:|---:|---:|---:|---:|
| 10009 | 835 | 0.688 | 0.976 | 1.3e-4 | 2.7e-4 |
| 100057 | 8339 | 0.641 | 0.994 | 9.6e-5 | 1.9e-4 |
| 1000033 | 83337 | 0.684 | 0.973 | 1.4e-4 | 2.8e-4 |
| 10000141 | 833346 | 0.686 | 0.981 | 7.7e-5 | 1.5e-4 |

All 2048 class-walk trials returned the correct logarithm (no degenerate,
cap, stack or verification outcomes in the uncapped runs); maximum stack
depth 22. The same detector without automorphisms (n=q) needs 3.51-3.57 times
more steps (sqrt(12)=3.46). A stress run with a 32-entry table measures
2.61e-3 fruitless 2-cycles per step against the predicted 1/(12*32)=2.60e-3,
waste 5.3e-3, all trials solved. These are walk statistics at 2^10-2^23
scale only; extrapolation to n ~ 2^250 is part of H1-H2.

On the real field, experiment.py checks the Frobenius matrices against
exponentiation, class invariance of rep on random GT elements (all 12
conjugates map to the same representative), and a 24-step production-field
class walk whose final state satisfies Y = g^a z^b and is canonical, and a
48-step production-field run with a one-entry table that exercises the
squaring escape (seed 20261111: 4 fruitless 2-cycles detected and escaped,
final invariant holds).

## 6. Work and score

Any run of s <= E steps costs at most s*I with I = 1,392,640 scalar
instruction equivalents (machine.md: G+H+A+K+S+control; the stack share S is
amortised over the run because pops never exceed pushes), plus setup, start
element, tie budget, final verification and output at most F = 2^42
(itemised and asserted in experiment.py, tally about 2^40.08). The step cap is
E = 76912954246949046915525912469710505284:

    R_total = E*1,392,640 + 2^42,   W = R_total/20,000,
    log2 W <= 131.976245;  claimed = 131.976245 + 0.02 = 131.996245.

The 0.02-bit margin covers rounding and bookkeeping slack; it is a stated
upward adjustment, not an estimate. Against the baseline: the class count
gives sqrt(12) (1.79 bits), one evaluation per step instead of three, and the
lower cap give the rest, minus 0.29 bits for canonicalisation and stack work (I = 1.22 times the baseline per-evaluation bound 3,420,160/3).

Memory < 2^22 bytes; preprocessing (all inside F) < 2^41 instructions, about
2^26.7 reference units; advice zero. Public profile constants and algorithm
design are the only excluded items.

## 7. Pairing-inversion track

The same run recovers and verifies x, then returns P=P0 and Q=[x]Q0 by affine
double-and-add (at most 508 affine G2 operations plus an optional [r]Q check,
<2^35 instructions, inside F). By bilinearity e(P0,[x]Q0)=g^x=z; both points
lie in the required subgroups; z=1 uses the permitted infinity output. The
reduction adds no work beyond F, so the same bound applies.

## References

Pollard, Math. Comp. 32 (1978). Teske, Math. Comp. 70 (2001) 809-825.
Wiener, Zuccherato, SAC 1998, LNCS 1556. Duursma, Gaudry, Morain,
ASIACRYPT 1999, LNCS 1716. Gallant, Lambert, Vanstone, Math. Comp. 69 (2000).
Bos, Kleinjung, Lenstra, ANTS-IX 2010, LNCS 6197. Bernstein, Lange,
Schwabe, PKC 2011, LNCS 6571. Nivasch, Information Processing Letters 90(3)
(2004) 135-140. Flajolet, Odlyzko, EUROCRYPT 1989 (random mapping
statistics). The arguments used are reproduced above; no retrieval is needed.
