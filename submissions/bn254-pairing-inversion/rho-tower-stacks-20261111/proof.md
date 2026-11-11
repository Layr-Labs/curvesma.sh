# BN254 unrestricted pairing inversion by a tower-arithmetic, stack-detected Pollard rho

This is a generic-group attack on the protected 254-bit prime-order GT,
followed by the permitted output (P0, [x]Q0). It is not a best-known
finite-field attack, a novelty claim, a practical break, or a security-bit
estimate. Its work bound is a uniform **every-run cap** (no expectation, no
run gate). Its success argument uses the standard, explicitly disclosed
random-mapping heuristic for r-adding walks; this seeks the protected
exploratory `plausible_not_refuted` standard, not formal certification.
No earlier result is assumed accepted.

It keeps the reviewed baseline's structure (one bounded trial, exact
scalar-label recovery, verified output, degeneracy argument) and changes
three things, each justified below and in `machine.md`:

1. **Arithmetic.** GT products use the standard BN254 tower
   Fp2 = Fp[u]/(u^2+1), Fp6 = Fp2[v]/(v^3-xi), Fp12 = Fp6[W]/(W^2-v),
   xi = 9+u, with Karatsuba at each level: exactly 54 Fp Montgomery
   products and 280 Fp additive operations per Fp12 product (counted by
   instrumenting `rho_ref.py`; checked against the protected-basis
   schoolbook product in `experiment.py`).
2. **Partition.** The step index is the top 10 bits of h1*w mod 2^64, where
   w is one fixed 64-bit word of the canonical state and h1 is a fresh odd
   64-bit salt; the walk is 1024-adding. No SHA-256 per step.
3. **Cycle detection.** One walk, one multiplication per step, with
   Nivasch's stack algorithm partitioned into K=256 stacks, instead of
   Floyd's three multiplications per iteration.

## Basis

The protected basis is Fp[w]/(w^12-18w^6+82) with u = w^6-9. Since
(w^6-9)^2 = -1, w^6 = 9+u = xi. Put W = w, v = w^2. A protected element
sum c_i w^i equals sum_{k=0..5} A_k w^k with
A_k = (c_k + 9c_{k+6}) + c_{k+6} u in Fp2, and grouping even and odd k
gives the tower element ((A0,A2,A4),(A1,A3,A5)). This is an Fp-linear
bijection with explicit inverse c_{k+6} = Im(A_k), c_k = Re(A_k) - 9Im(A_k).
It is applied once to z, g and on output only. `experiment.py` checks the
round trip and 57 full-width products (random and edge-valued) against the
protected-basis product. Element equality in the tower is equality of the
protected element.

## Complete algorithm

Input z (12 canonical coefficients). Decode, check z != 0 and z^r = 1,
convert to Montgomery tower form. If z = 1 return (P0, infinity), which the
profile permits. Fix the public parameters

    n'  = r * 1024/1023,
    T   = ceil(1.15 * sqrt(n'))      (exact integer: least T with
                                      400*1023*T^2 >= 529*1024*r),
    T   = 170221906881474042593545034831872327888 < 2^128.

One trial, fresh coins:

* Sample two odd 64-bit salts h1, h2 and 2*1024+2 uniform scalars modulo r
  by bounded 254-bit rejection (at most 1024 draws each).
* For j = 0..1023 build M_j = g^{c_j} z^{d_j}; set Y_0 = g^{a_0} z^{b_0}.
* State (Y, a, b) with invariant Y = g^a z^b. Let w(Y) be the low word of the
  Montgomery residue of Re(A_0) of Y. Index j(Y) = (h1*w mod 2^64) >> 54,
  class k(Y) = (h2*w mod 2^64) >> 56. Total order on states: lexicographic
  order of the 48 Montgomery words, w first.
* K = 256 empty stacks S_0..S_255, each holding at most 512 entries
  (Y, a, b).
* For t = 0, 1, ..., T: with Y = Y_t:
  * In stack S_{k(Y)}: pop while the top is strictly greater than Y.
  * If the top equals Y as a group element: collision found (below).
  * If the stack holds 512 entries: report failure. Else push (Y, a, b).
  * If t = T: report failure (cap). Else (Y, a, b) <- (Y M_j, a+c_j, b+d_j)
    with j = j(Y), labels modulo r.
* On a collision between the stored (Y, a', b') and current (Y, a, b): if
  b - b' = 0 mod r report failure; else x = (a' - a)/(b - b') mod r. Verify
  g^x = z (failure otherwise). Output P = P0 and Q = [x]Q0.

`rho_ref.walk` is this control flow verbatim for an abstract group; the toy
runs execute it.

## Correctness

The invariant holds initially and after each step by exponent addition.
Equal group elements g^{a'} z^{b'} = g^a z^b with z = g^x give
a' - a = x(b - b') mod r; with b - b' nonzero the quotient is the unique
discrete logarithm. The output check g^x = z is exact; a returned x is never
wrong. Then e(P0, [x]Q0) = e(P0,Q0)^x = g^x = z by bilinearity of the
protected reduced pairing and the profile definition g = e(P0,Q0). [x]Q0
lies in the order-r subgroup <Q0>; P0 generates G1. Identity: (P0, infinity).

**Detection (Nivasch).** The map Y -> Y M_{j(Y)} is deterministic on group
elements, so Y_t is eventually periodic with tail mu and cycle length
lambda, rho = mu + lambda (first repeated index). Each stack is kept sorted
(strictly increasing towards the top) because an entry is pushed only after
all larger ones are popped. Take any class c met on the cycle and its
minimum cycle element m_c, first reached at time t_c in [mu, mu+lambda-1].
After t_c, every element of class c is either a cycle element of class c
(not smaller than m_c, so m_c is never popped; equal only at its next
occurrence) or nothing else, since only cycle elements occur after mu.
At time t_c + lambda the walk returns to m_c, every entry above it has been
popped, and equality is detected. Hence detection happens no later than
mu + lambda + D, where D = min over classes c of (t_c - mu) < lambda. An
earlier equality is also a genuine repeat of a group element and is equally
usable. Stacks never contain equal elements before detection.

**Degeneracy (exact, no heuristic).** Fix the true x, the salts, the
increment logarithms s_j = c_j + x d_j and s_0 = a_0 + x b_0. Conditional on
these, the d_j remain independent and uniform (each has exactly one c_j
consistent with s_j), and the entire group walk, all stack operations, the
classes, the detection time and the two detection times t' < t depend only on
(s_j, s_0, salts), not on the d_j. Then b - b' = sum_j n_j d_j where n_j counts
uses of index j in steps t'..t-1; sum n_j = t - t' with 0 < t - t' <= T < r, so
some n_j is nonzero mod r, and the combination is exactly uniform: the
degeneracy probability is exactly 1/r < 2^-253 per run. Rejection sampling
fails with probability < 2050 * 2^-1024.

## Success heuristic, support and limits

**Model.** Under the random-mapping heuristic on n' points (the standard
Brent-Pollard/Teske correction for an r_a-adding walk replaces n by
n/(1 - 1/r_a); r_a = 1024 here), the first-repeat time satisfies
Pr[rho > t] = exp(-t^2/(2n')) to leading order, the repeated earlier index
is uniform so lambda = rho*V with V uniform on (0,1], and the cycle-minimum
offsets of the K class minima are close to independent uniform positions on
the cycle, so D = lambda * W with Pr[W > w] = (1-w)^K. Success requires
rho + D <= T. With s = T/sqrt(n') = 1.15 the integral
`experiment.py:model` gives **0.4825** (the simple bound without the D
term is 1 - exp(-s^2/2) = 0.4838; D costs only 0.0013 because K = 256).
The declared success is **at least 0.39**.

**Margin against model error.** If the effective number of points were 20%
larger than n' (i.e. the walk were noticeably worse than the random-map
model with the r-adding correction), the same cap gives 0.4224 > 0.39.
Teske's experiments (Math. Comp. 70 (2001) 809-825) found 20-adding walks
already statistically indistinguishable from random mappings; 1024-adding is
far closer. The factor 1/sqrt(1-1/r_a) is the usual heuristic r-adding
penalty (Brent-Pollard; discussed by Bernstein and Lange, e.g. "Two grumpy
giants and a baby", ANTS 2012); at r_a = 1024 it is 1.0005, so even an error
of many times this size is inside the 20% margin. The references are
context; the argument above is self-contained.

**Toy evidence (actual control flow).** `experiment.py` runs the exact walk,
1024-adding table, salted multiplicative hash, 256 partitioned stacks with
the 512 cap, hard step cap ceil(1.15 sqrt(n*1024/1023)), recovery and
verification on fresh prime-order groups (additive Z_q at 16, 20, 22 bits and
a multiplicative subgroup of F_P^* at 16 bits). Example output (seed
20261111, toy_trials 256): success fractions 0.430, 0.473, 0.500, 0.500,
pooled 0.459 over 608 trials, versus the model's 0.4825; all recovered
logarithms were correct; no degenerate, overflow or verification failures;
maximum stack depth 8. These are finite-scope statistics, not a production
proof: the toy groups are far smaller than 2^127, the toy table has more
entries than the toy walk has steps in the 16-bit case, and the toy hash
word is the group element itself rather than one Montgomery limb.

**Hash word.** j and k depend on one 64-bit Montgomery limb of one Fp
coordinate. For a GT element produced by a long random-looking walk, this
word is an essentially uniform value; distinct elements may share it, which
is harmless (full 48-word comparison decides order and equality). The
material assumption is that the salted index on this word gives a walk whose
collision time follows the r-adding random-map model at production scale,
and that cycle-class minima are positioned as in the model. These are
standard Pollard-rho heuristics, disclosed, not proved.

**Stack overflow.** Before detection, the elements of a class are distinct,
and the entries of its stack are exactly the suffix minima of that class's
sequence so far. Under the same heuristic (random relative order), their
number for m <= 2^128 items is a sum of independent Bernoulli(1/i) with mean
H_m < 89.7; Chernoff gives Pr[>= 512] <= (e*89.7/512)^512 < 2^-548, and a
union bound over all times and 256 stacks leaves < 2^-412.

**Total.** Pr[success] >= 0.4825 - 2^-253 - 2^-412 - 2050*2^-1024 under the
model; >= 0.42 with a 20% model penalty. Declared: >= 0.39.

## Computation, memory, output

`machine.md` bounds each step by S = 49,152 scalar AArch64
instruction-equivalents (Fp12 product 42,752 including 4,096 overhead; stack,
hash, labels and control 2,048; slack 4,352), and all setup, randomness,
table construction, subgroup check, recovery, verification and output by
F = 2^40. The cap gives the uniform every-run bound

    R_total = T * 49,152 + 2^40,   W = R_total / 20,000,
    log2 W = 128.297934444969765...,  submitted 128.297935 (rounded up).

Peak memory <= 2^26 bytes (stacks 60.8 MB worst case, table 0.46 MB).
No precomputation is reused across targets, no advice, no amortization.
