# Tag-tracing 8-adding distinguished-point rho in BN254 GT, then (P0,[x]Q0)

Lineage. This package is version 2 of e-man07's tag-tracing rho
73463878-0c19-4309-a991-749ea525d51c (bn254-pairing-inversion,
121.285447), which itself builds on
1da4a817-8e05-49ad-87ae-feabe0fecf72 and
cdae944e-4e31-4a5b-ac52-be49ad2618a0. The algorithm, the tag certificate,
the step-loop listing and its count (273), the output (P0,[x]Q0), the toy
experiments and the success analysis are those of 73463878. New in v2:
the window re-base (one Fp12 product and two label additions per 256
steps) is now built from the fused lazy-reduction Fp2 kernels f2mul and
f2mulxi and the base-field leaves fp_add and fp_sub of e-man07's
7e58d165-ebf6-41aa-8003-2f0a42370b6a, copied unchanged, together with
that package's idea of storing the table operand's Karatsuba sums next
to each table entry. Each table entry therefore grows from 1,280 to
2,048 bytes. A new straight-line driver `rebase` calls these leaves.
`experiment.py` executes the whole re-base instruction by instruction on
a simulator of the AArch64 subset used (11,730 instructions, the same on
every input), so the window end costs 11,923 instead of 24,576
instructions and the claim falls from 121.285447 to 121.077977. Also
carried over from our gt-dlp port of 73463878
(76bccbcf-95f2-46da-bce2-481011ac5601): machine.md states the width of
every meter, a store probe is priced at 64 instructions instead of 32,
and the self-cycle bound of section 6 is written out exactly. The walk,
the step cap and every other count are unchanged.

A generic square-root attack on the protected 254-bit prime-order GT. The
walk is an ordinary 8-adding Pollard-rho walk, but almost no step performs
a field multiplication: following Cheon, Hong and Kim's *tag tracing*
("Accelerating Pollard's rho algorithm on finite fields", J. Cryptology 25
(2012)), a step only computes a 64-bit approximation of one Fp coordinate
of the current element, as an inner product with precomputed table data,
and follows a stored pointer; one real Fp12 product is made every 256
steps. The method and every fact needed are reproduced here; the citation
is context, not an imported premise. This is not index calculus/NFS, not a
practical break (the table alone needs about 2^59.90 bytes) and not a claim
of optimality. The work bound is an analytical hard cap; the success bound
uses disclosed random-walk heuristics and seeks the exploratory
`plausible_not_refuted` standard. No prior result is assumed accepted.

Notation: p, r, g are the protected modulus, group order and generator;
z=g^x is the input. An element Y of Fp12 has twelve Fp coordinates
y_0..y_11 in a fixed basis (the tower order of machine.md); t(Y) in [0,p)
denotes its coefficient of w^0 in the protected basis, a fixed Fp-linear
functional of the coordinates.

## 1. The walk as a function of the group element

Fix eight multipliers M_s=g^(c_s)*z^(d_s), s=0..7. Define

    s(Y) = floor(8*t(Y)/p) in {0..7},      f(Y) = Y*M_{s(Y)},
    Y is distinguished  iff  t(Y)*2^80 < p.

f, s and the distinguished predicate depend on the group element only. The
algorithm below evaluates exactly these functions; the table machinery is
only a way to evaluate s(Y) cheaply. Hence the walk is a standard 8-adding
walk with a deterministic iteration function, and two trails that reach
the same element coincide from there on.

## 2. Tag certificate (exact statement, proved)

Row. For a fixed T in Fp12, Y -> t(Y*T) is Fp-linear, so there is a row
v(T) in Fp^12 with t(Y*T) = sum_k y_k*v_k mod p (explicit formula in
`reference.py`, `BN254GT.row`; checked on 40 random pairs).

Stored data. Coordinates are kept as 4-limb Montgomery residues
Y~_k = y_k*2^256 mod p with limbs Y~_{k,i}, i=0..3. For T put
u_{k,i} = v_k(T)*2^(64i-256) mod p and store the 128-bit integers
F_{k,i} = floor(2^128*u_{k,i}/p) = Fhi*2^64+Flo.

Tag. A = sum over the 48 pairs (k,i) of
( UMULH(Y~_{k,i},Flo_{k,i}) + Y~_{k,i}*Fhi_{k,i} ) mod 2^64.

**Lemma.** With t=t(Y*T), the real number tau = 2^64*t/p satisfies
(tau - A) mod 2^64 in [0,96).

Proof. sum_i Y~_{k,i}*u_{k,i} = Y~_k*v_k*2^-256 = y_k*v_k (mod p), so
t/p = sum_{k,i} Y~_{k,i}*u_{k,i}/p (mod 1). For one pair write y=Y~_{k,i}
<2^64 and u/p = F/2^128+eps, 0<=eps<2^-128. Then
y*u/p = y*Fhi/2^64 + y*Flo/2^128 + y*eps. Modulo 1 the first term equals
(y*Fhi mod 2^64)/2^64 exactly; the second equals UMULH(y,Flo)/2^64 plus a
remainder in [0,2^-64); the third lies in [0,2^-64). So each pair
contributes its two computed words plus an error in [0,2*2^-64), and the
48 errors add up to less than 96*2^-64. Multiply by 2^64. QED
(`experiment.py` checks the lemma on 260 full-width cases including edge
values; largest observed offset 41 in the declared run.)

**Corollary (certified step).** Call A *certified* if
((A+96) mod 2^64) mod 2^61 > 96. Then the low 61 bits of A lie in
[1, 2^61-97], so [A,A+96) contains no multiple of 2^61 and does not wrap:
tau and A have the same top three bits, i.e. s(Y*T) = A>>61 exactly, and
tau>=A>=1, i.e. t/p>=2^-64>2^-80: the element is not distinguished.

## 3. Complete algorithm

Parameters: R=8 multipliers; window length 256; distinguished density
2^-80; step cap S=floor(sqrt(1409*r/1000))+1; at most 2^50 walks, each at
most 2^86 steps; at most 2^90 exact evaluations; at most 2^70 store probes.
Reaching any cap ends the run as a failure. One trial, no restart.

1. Decode z, check z!=0 and z^r=1. If z=1 output P0 and infinity.
2. Draw c_s,d_s uniformly mod r (bounded rejection sampling) and form the
   eight M_s. Build the table: one entry e for every multiset of at most
   256 multipliers (C(264,8)=525,783,425,977,953 entries), holding
   T_e = product of the multiset, the 96 words F(T_e), alpha_e=sum c_s and
   beta_e=sum d_s mod r, eight child pointers (entry of the multiset
   plus one copy of M_s; null at size 256) and twelve stored Fp2 sums of
   the coordinates of T_e used by the re-base product (machine.md). The
   root is the empty multiset.
3. For each walk: draw a,b uniformly mod r, set Y0=g^a*z^b (skip if 1),
   e=root. The pair (Y0,e) represents the element Y=Y0*T_e with label
   (a+alpha_e, b+beta_e). Repeat:
   a. A = tag(Y0,e).
   b. If A is certified: s=A>>61. If e has a child pointer for s, set
      e to it; otherwise (window end) *re-base*: Y0<-Y0*T_e,
      (a,b)<-(a+alpha_e,b+beta_e), e<-child of the root for s.
      This is one step Y<-Y*M_s.
   c. If A is not certified: re-base to e=root, read t=t(Y0) exactly,
      and evaluate s(Y) and the distinguished predicate from t by their
      definitions. If Y is distinguished leave the loop with (t,a,b);
      otherwise set e to the root's child for s (one step).
4. At a distinguished point look t up in the store. If an entry (t,a',b')
   exists, put x=(a-a')/(b'-b) mod r (failure if b=b'); if g^x=z go to 5,
   otherwise treat the hit as a false match. Insert (t,a,b); next walk.
5. Output P=P0 and Q=[x]Q0 (affine double-and-add over Fp2).

`reference.py` (`solve`) is the executable reference with this control
flow, generic over the group, except that it omits the two meters (exact
evaluations, store probes); these only add failure exits and are charged
in section 5. `ttstep.c` is the priced loop for 3a-3b; the re-base of 3b
and 3c is the listing `rebase` in `listings.txt`.

## 4. Correctness

By the corollary, step 3b applies exactly f of section 1 and never skips a
distinguished point; step 3c uses the definitions directly. So every walk
is a trail of f from its start, stopped at its first distinguished point,
whatever representation (Y0,e) the element happens to have.
`experiment.py` drives the real `solve` loop at full BN254 width for 120
steps and checks at every step, against independently multiplied
elements, that the selector equals floor(8t/p), that consecutive elements
differ by exactly the selected multiplier, and (every 40 steps) that
Y=g^a*z^b for the represented label; it repeats the same walk with a
deliberately narrow certificate (54 protected bits, 14 exact evaluations)
and obtains the identical element sequence.

The re-base is the only place where the stored window base changes. Its
listing reads and writes the coordinates as fully reduced 4-limb
Montgomery residues (y*2^256 mod p, least significant limb first) in the
coordinate order of section 2, which is exactly the input of `tt_run`.
`experiment.py` executes it on 32 inputs (real GT elements, arbitrary
Fp12 elements and edge values such as all coordinates p-1, labels r-1,
in place and into a second buffer) and checks that the 56 output words
equal the Montgomery limbs of Y0*T_e computed independently in the
protected basis, followed by the two label sums mod r, and that the tag
of the output certifies coordinate 0 of Y0*T_e.

Labels: Y=g^a*z^b holds at each start and is preserved by multiplying with
M_s and adding (c_s,d_s). Two distinguished points with equal t that are
the same element give g^a z^b = g^a' z^b', so x=(a-a')/(b'-b) whenever
b!=b' (r is prime). The check g^x=z makes a wrong output impossible; a
false match (equal t, different elements) fails it and costs one
exponentiation, metered in the fixed allowance. By bilinearity
e(P0,[x]Q0)=g^x=z, and both points are in the prescribed subgroups; for
z=1 the pair (P0, infinity) satisfies e(P0, infinity)=1.

Degeneracy is exactly 1/r, without heuristics. Fix x and condition on the
logarithms m_s=c_s+x*d_s of the multipliers and s_i=a_i+x*b_i of all walk
starts. Every group element, tag, exact evaluation, distinguished point
and store hit is a function of these. Conditionally the b_i and d_s remain
independent and uniform. For a hit between walk i and an earlier walk j,
b-b' = b_i-b_j + (an integer combination of the d_s), which is uniform on
Z/r because b_i is; it vanishes with probability exactly 1/r.

## 5. Work bound (hard cap, independent of all heuristics)

From `machine.md`. Both priced pieces are counts of executed scalar
AArch64 instructions (every instruction, load, store, branch, call and
return counts one; no SIMD), which is the unit of the cost model:

| item | instructions | source |
|---|---:|---|
| one certified step (one trip of the `tt_run` loop: 48 UMULH, 48 MADD, pointer, exits, counter) | 273 | static count of the clang listing `ttstep.disasm.txt`, recounted by `experiment.py` |
| one window end: `rebase` 11,730 (Fp12 product by the table entry and two label additions mod r), `tt_run` entry/exit 65, driver <= 128 | <= 11,923 | `rebase`: executed on the simulator in `experiment.py` (listing `listings.txt`); 65: clang listing; 128: itemised in machine.md |
| everything else, by meters (table build, <=2^50 walk starts, <=2^90 exact evaluations, store, recovery, output) | <= 2^110 | machine.md |

At most S steps are taken and at most floor(S/256)+1 window ends occur
(a window end is preceded by 256 steps taken since the previous re-base).
Hence, for every run,

    R_total <= 273*S + 11,923*(floor(S/256)+1) + 2^110,
    W = R_total/20,000,   log2 W <= 121.077977.

That is 273 + 11,923/256 = 319.57 instructions (0.015979 reference units)
per step, against 369 in v1. Nothing is priced by substituting the
normalization unit for a BN254 operation. Peak storage is below 2^60.4
bytes (table 2^59.901 at 2,048 bytes per entry, store 2^58.585, the rest
negligible); every address is an ordinary 64-bit pointer below 2^63.
Preprocessing (the table) is below 2^69 instructions and is included.
There is no advice and no amortization over other targets.

## 6. Success probability: heuristics, support and limits

H1 (collision statistics). The trails of f, started at independent uniform
elements, meet as in the random-map model corrected for an m-adding walk:
after s steps in total the probability that some trail has touched an
earlier trail is at least 1-exp(-s^2/(2*kappa*r)) with kappa=m/(m-1)=8/7
(the Brent-Pollard/Teske correction). With S^2>=1.409*r this gives

    1-exp(-1.409*7/16) = 0.4601.

The declared 0.39 needs less: 1-exp(-1.409/(2*kappa)) >= 0.39 holds for
every kappa <= 1.425, against the assumed 8/7 = 1.143. So the claim
tolerates a collision law 25% slower than the corrected model.

H2 (tag bits). The low 61 bits of A, and the leading bits of t/p at
distinguished scale, behave like uniform bits along the walk.

Losses, all explicit; each event ends the run as a failure:
* A touch in the last 2^86 steps may not reach its distinguished point
  before the cap. The touch rate per step is at most S/r, so this is at
  most 2^86*S/r < 2^-40.
* Under H2 a walk fails to distinguish within 2^86 steps with probability
  exp(-64); a walk closes a cycle on itself before distinguishing with
  probability about E[len^2]/(2r) = 2^161/(2r) = 2^160/r. At most 2^50 walks: < 2^-40.
* The 2^50 walk cap binds only if 2^50 walks average under 2^77.05 steps
  against a mean of 2^80 (expected number of walks S*2^-80=2^47.05): under
  H2 this is a large-deviation event, < 2^-40.
* Exact evaluations: a step is uncertified with probability 97/2^61 under
  H2, so their expected number is below 2^74; by Markov the 2^90 meter is
  exceeded with probability < 2^-16.
* Store probes: a table of 2^52 slots holding at most 2^50 keys with
  linear probing has expected O(1) probes per operation; by Markov the
  2^70 meter is exceeded with probability < 2^-16 under H2.
* Degeneracy 1/r (exact); scalar rejection sampling fails with probability
  below 2^-1000; a start equal to 1 has probability 1/r per walk.
These sum to less than 2^-14. No other event prevents success once a
touch has happened: the later trail follows the earlier one to its stored
distinguished point. The heuristic prediction is therefore above 0.4600,
against the declared **success >= 0.39** (the protected minimum). The
run stops at the step cap in the remaining cases (about 54% under the
model); those runs are failures and are inside the same work cap. The
identity target succeeds directly.

Support for H1. (a) The selector is the classical one for rho in finite
fields (leading bits of a coordinate), and the walk is the textbook
r-adding walk; Teske ("On random walks for Pollard's rho method", Math.
Comp. 70 (2001)) studies exactly this family, and 8 <= m <= 20 is the
range where the 1/(1-1/m) correction is in use. (b) `experiment.py` runs
the complete algorithm (same `solve`, same tag arithmetic on 64-bit words,
same table, selector, cap formula, recovery and failure rules) on
subgroups of prime order n of 61-bit prime fields, with a fresh field,
target and coins in every trial, and reports all outcomes. Development
runs with the shipped cap formula are recorded in
`development-evidence.json`: declared input (seed 20261111, 256/256/128
trials at n of 16/20/24 bits) success 0.520, 0.457, 0.469 (two earlier
layouts of the same experiment, with only 64 trials in the last row, gave
0.438 and 0.391 there; a 64-trial row has standard error 0.06); a larger
offline run of the same code (seed 20261111, 4000/3000/2000/400 trials at
16/20/24/28 bits) gave 0.511, 0.452, 0.457, 0.498. All are above 0.39;
they scatter around the model value 0.46 rather than clearly above it, so
the margin over 0.39 is about 0.07 in the model and no extra slack is
claimed.

Limits. Toy groups have 16 to 28 bits and a one-coordinate field, so the
selector there reads the whole element; toy distinguished densities are
2^-3 to 2^-6, where the certificate is taken over the distinguished bits
as well (`cbits=dbits` in `solve`), while production uses the three
selector bits and finds distinguished points only through exact
evaluations. In production the selector reads one of
twelve coordinates of an element of a 2^254 subgroup of a 3048-bit field,
and trails have 2^80 steps. Equidistribution of that coordinate on GT and
the extrapolation of the 8-adding collision law to 2^127 steps are the
disclosed heuristics, not theorems. No production discrete logarithm was
computed. The work cap does not depend on H1 or H2. Sieve methods for
the degree-12 field are expected to be cheaper; this is not a security
estimate.
