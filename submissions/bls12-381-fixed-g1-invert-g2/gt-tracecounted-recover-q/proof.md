# Frobenius-class parallel rho with a cyclotomic-squaring branch read from a trace invariant

**Claim.** For uniform nonidentity P in G1 and uniform z in GT the algorithm below
recovers Q in G2 with e(P,Q) = z with probability at least 0.39, using at most
**2^125.62** dense-Fp12-multiplication equivalents (C = 20,000), including setup, the two
pairings, table construction, walk starts, distinguished-point (DP) handling, cycle
escapes, the DP detection delay, output and verification. Peak memory is below 2^74
bytes, preprocessing (included) at most 2^31 units, nonuniform advice 0 bytes. The success
bound rests on the usual random-mapping heuristic for rho walks, adjusted for the branch
mix (section 4), and is supported by toy-field experiments and by running the production
step rule in the real GT. No production-size run is claimed.

The class construction (Frobenius classes of GT, the 4-multiplication key, the twisted
multiplier table that never applies a Frobenius to the walk element) is the one published
in submission 7da5f23a-6633-46b1-9c03-11f6cfedeeec (mpjunior92), which this package builds
on. New here: (i) a cyclotomic-squaring branch that is cheaper than the table multiply;
(ii) a fresh instruction ledger; (iii) an analysis of how the branch mix changes the
collision law and the fruitless cycles; (iv) execution of the production step rule in the
real BLS12-381 GT; (v) explicit memory and report caps; (vi) the branch, multiplier and DP bits are read from the free class
invariant T = Re(c_0) instead of the minimum key K, so a squaring step needs no key computation at all; (vii) a bounded
distinguished-point report store (flat buffer, fixed eight-pass radix sort, linear scan), priced in section 3.2.
The package `frobenius-rho-squaring-mix` (submission 6d4e57d1, accepted as valid at 125.80)
is the same design with the bits read from K; this one differs only by (vi) and its ledger.

## 0. What this revision changes

**Port to this track.** This package is the fixed-G2, recover-G1 package of submission `56fc259f`
(`gt-tracecounted-exact-v25`, 125.72) carried to the fixed-G1, recover-G2 target. Both tracks pose the
same discrete logarithm in GT; only the fixed argument and the output group differ. The port changes
section 1 (g = e(P, Q0) for the supplied P, output Q = [x]Q0 in G2, the identity target), the setup and
output lines of section 5 (the two pairings and the G2 output, now counted operation by operation in
`experiment.py` part F), and adds part F. The walk, branch, keys, escapes, caps, prices, H1, H1+ and every
other figure below are those of `56fc259f`; the rest of this section describes how that package revised
its own parent, and "this package" there refers to the shared price vector.

**Bounded report store, and a measured baseline in the package.** A further revision replaces the
inherited chained multimap with one bounded structure, so that the per-report charge does not rest on a
load factor or a hash-table expectation (the defect that made the 125.89 package fail review). Reports
are appended to a flat buffer of 80-byte records and, after the walk budget is spent, sorted by the
128-bit fingerprint with eight fixed 16-bit radix passes and scanned for equal neighbours. Section 3.2
itemises the per-report worst case at 3,831 instructions, charged at 4,000. With the report cap
D at most 2^66.12 the DP-handling term is at most 2^63.80 units, the non-walk bucket stays below
2^64.5, and every figure of the walk ledger, including W = 2^125.6197624 and the claim 125.62, is
unchanged. `experiment.py` part D prints the new term and part G exercises the same eight-pass sort and
scan.

The measured library counts this package quotes are recorded in `measured-kernel.md`, in the package:
blst at commit f7a88067, the Callgrind harness in the grader's arm64 image, `blst_fp_mul` = 383,
`blst_fp12_cyclotomic_sqr` = 9,803, `blst_fp12_sqr` = 16,358 and `blst_fp12_mul` = 22,996 instructions,
with the arm64 primitive lengths (`sqr_382x`, `__add_mod_384x384`, `__redc_tail_mont_384` and the rest),
the compiled kernel written over those primitives, and its measured per-part split. Sections 0, 5 and 6
and `experiment.py` cite that file for the 383, the compiled kernel charged at 9,709, and the 9,803 control.

Everything below — the class walk, the 4-multiplication key, the twisted table, the
cyclotomic-squaring branch read from T, the escapes, the caps, H1, H1+, the toy and
production-field evidence — is the package of submission `747ad6c9`
(`frobenius-rho-trace-branch`, jungjipdo, 125.74), which is itself a revision of `6d4e57d1`
(125.80) and of the Frobenius-class construction of `7da5f23a` (mpjunior92). Only the three
charged prices, the report store and the figures that follow from them change.

That package itemises its three per-step charges and then charges each one above its own
counted listing:

| Work | counted or measured listing | parent's charge | this package |
| --- | ---: | ---: | ---: |
| every step | 270 | 300 | **270** |
| multiply steps, keys only | 2,119 | 2,200 | **2,119** |
| cyclotomic squaring | 9,708.0022 compiled kernel, measured (11,369 by the narrow 9 + 34 + 200 listing) | 11,400 | **9,709** |

The parent's charges are roundings of its own listing (300 for 270, 2,200 for 2,119, 11,400
for 11,369): 30, 81 and 31 instructions of slack, 111 instructions on a multiply step and 30
on a squaring step. This package charges the audited listing itself for the first two; for
the squaring step it replaces the narrow listing (11,369 instructions) by the same
Granger-Scott squaring written as one straight-line wide-lazy C kernel over blst's exported
arm64 primitives, compiled and measured in the grader's arm64 image: 9 wide Fp2 squarings
(3,690), 24 wide Fp2 add/subs (1,473), 12 Montgomery reductions (2,448), the narrow sums,
chains and `xi*C1` (1,708), the C body (379) and the loop (10), i.e. **9,708.0022
instructions per call**, rounded upward to **9,709 instructions = 0.48545 unit** for the charge. That is 95.0 instructions below the control
`blst_fp12_cyclotomic_sqr` (9,803.0018) in the same run, and the kernel is bit-exact against
that control on 1,000 random cyclotomic elements; the primitive bodies are blst f7a88067's
own `sqr_382x`, `__add_mod_384x384`, `__sub_mod_384x384`, `redc_mont_384` and
`mul_by_1_plus_i_mod_384x`, and `measured-kernel.md` reproduces the harness, the image and
the per-part split. Every dependent figure is regenerated: multiply step 1.125 -> **1.11945**
units, squaring step 0.585 -> **0.49895** units, and the squaring probability re-tuned from
3/8 to **35/64** (the model optimum at the new price; the toy at 35/64 gives 0.419 +- 0.011
against the model 0.43), so the average step cost is c 0.9225 -> **0.7801140625**
(= 0.846 x 0.9225). N' rises because more squaring steps make the walk less collision-prone
per pair at the same step count: sum q^2 0.141006 -> **0.299273**, log2 N' 251.4914 ->
**251.7852**, S(0.43) 2^125.8302 -> **2^125.9771**, and W_nominal 2^125.7147 ->
**2^125.6198**, claim 125.74 -> **125.62**. The algebra, eps, E, the 2^64.5 non-walk bucket, the report and start
caps and the success chain are unchanged; `experiment.py` part D recomputes the parent's
ledger (2^125.7147, c = 0.9225, at its own s = 3/8) and this one through the same code
path. Section 3.4 restates the abort headroom under the new cushion.

## 0a. Representation: which limbs the rule reads

The parent states its rule on ordinary integers in [0,p): T = Re(c_0) "read as an integer", and the
twelve keys ordered as integers. Its arithmetic prices are those of the pinned blst kernels, and those
kernels store field elements in Montgomery form x~ = x*2^384 mod p, so reading a stored limb does not
yield that integer. Converting costs work the parent's 270 / 2,119 / 11,369 charges do not contain —
one reduction per T read and twelve per multiply step, and the twelve alone are +0.14 unit on a
multiply step, which would put the parent's own bound above its claim.

This package therefore states the rule on the stored representation, which is what a limb-level
implementation of the design runs at zero cost:

* T is the six limbs of Re(c_0) **as stored**. The branch bits are (stored T) mod 64, the table index is
  (stored T >> 6) mod R, and the distinguished-point predicate is (stored T >> 16) mod 2^60 = 0. The
  fingerprint is read above them, (stored T >> 76) mod 2^128, so the branch, multiplier, DP and
  fingerprint bits are disjoint.
* kappa is the minimum of the twelve **stored** key products a*zeta^j and d*zeta^j, and the twist index
  is the smallest index attaining it.

This is the parent's rule composed with a fixed bijection of Z/p, so every structural property carries
over verbatim: the twelve stored keys are distinct; their multiset is the same for all twelve conjugates
(Frobenius permutes it), so kappa is class-determined and the twist index is covariant; stored T is
exactly invariant across conjugates because Re(c_0) is. What changes is which multiplier index is
selected, and that index is a pseudorandom function of the class under either rule, which is all H1
uses. No conversion is performed anywhere in the charged path, so the price table of section 5 is a
bound on the stated algorithm rather than on a representation it does not use.

Evidence, part C, both rules, 6,000 production-field steps each from the same seed: the stored rule
gives 0 class-check failures in the first 150 steps, the exponent invariant y = g^alpha z^beta holds at
the end, squaring fraction 0.5470, chi-square of the low 3 bits 3.71 (7 dof) and of the multiplier bits
(the 4 bits above the branch) 17.91 (15 dof); the integer rule gives 0 failures, 0.5484, 4.23 and 15.70.
The selected index differs in 2,483 of 2,716 multiply steps, so the two rules are different walks and
the claim is about the one the charges can pay for. `experiment.py` reports both blocks
(`partC_production_walk`, `partC_production_walk_integer_rule`).

## 1. Reduction

The challenge supplies a nonidentity P in G1, fixed within the instance (P = [b]P0 for a uniform b in
1..r-1, so not the generator across instances), and z uniform in GT. Compute g = e(P, Q0) once per
instance with the pinned convention pairing(Q0, P). P has prime order r and the pairing is
nondegenerate, so g != 1 generates GT (prime order r), and z = g^x for a unique x in [0, r).
The identity target is included: z = 1 exactly when x = 0, and the answer is then Q = [0]Q0 = O, the
point at infinity, since e(P, O) = 1 = z. The algorithm compares z with 1 (twelve coefficient
comparisons) before anything else and returns O without starting a walk, so this case (probability 1/r)
only adds to the success probability. Otherwise it finds x with the walk below, outputs Q = [x]Q0 in G2
and checks e(P, Q) = z, which holds by bilinearity: e(P, [x]Q0) = g^x = z. The answer is for the supplied
P: the table, the starts and every walk are built from this instance's g and z and are not reused.

`experiment.py` part F computes g with a transcription of the pinned library's pairing (the Miller loop
over the pseudo-binary encoding of 0xd201000000010000, its projective doubling, addition and line
functions over the twist, f_num / f_den and the final power (p^12-1)/r), checks g^r = 1 and g != 1 for a
sampled P, forms Q = [x]Q0 for a sampled x, checks that Q is on the G2 curve with [r]Q = O and that
e(P, Q) = g^x, and runs the identity case (Q = O, e(P, O) = 1). Outside the package the same g was
compared with py_ecc 8.0.0's pairing(Q0, P), whose four source files match the profile's sha256 values:
equal in all 12 coefficients for the sampled P.

## 2. Arithmetic facts (checked in `experiment.py`, parts A and C)

Tower: Fp2 = Fp[u]/(u^2+1), Fp12 = Fp2[w]/(w^6 - xi), xi = 1+u (the blst tower). With
u = w^6 - 1 this is the profile field, w^12 - 2w^6 + 2 = 0, and conversion of coefficients
is free: c_m = a_m - b_m, c_{m+6} = b_m.

**F1 Frobenius classes.** (sum c_m w^m)^p = sum conj(c_m) d_m w^m with d_m = xi^{m(p-1)/6}.
p has order exactly 12 modulo r, p^6 = -1 and 1 + p^4 + p^8 = 0 (mod r), so the maps
y -> y^{p^k} form a group of order 12 acting on GT, containing inversion. Every y != 1
lies in a class [y] of exactly 12 elements; N = (r-1)/12 classes, log2 N = 251.2721.

**F2 Keys.** Put zeta = xi^{(p^2-1)/6}; it lies in Fp and is a primitive 6th root of unity
(zeta^2 = zeta - 1). Let c_1 = a + b u be the w^1 coefficient of y and B = Re(conj(c_1) d_1)
= a Re(d_1) + b Im(d_1). The w^1 coefficient of y^{p^(k+1)} is conj(c^{(k)}) d_1, and
conj(d_1) d_1 = zeta, hence
`key_k := Re(w^1 coeff of y^{p^k}) = a zeta^j (k = 2j), B zeta^j (k = 2j+1)`.
With zeta^j in {1, zeta, zeta-1, -1, -zeta, 1-zeta} all 12 keys cost **4 Fp multiplications**
(a zeta, a Re d_1, b Im d_1, B zeta) plus 8 subtractions/negations and one addition. The
12 keys are distinct unless two of them coincide by a linear accident (probability <= 66/p
per element; ties are charged as failure). K([y]) = min_k key_k is a class invariant,
kappa(y) = argmin_k key_k, and y^{p^kappa} is the same element for all 12 members of the
class (the canonical representative y^). K is the minimum of 12 values that behave like
independent uniform residues (they come in negation pairs, so K lies in [0, p/2]); its low 200
bits are uniform up to statistical distance about 2^200/p = 2^-181 (the density of the minimum is smooth).

**F2b Trace invariant.** The w^0 coefficient of y^{p^k} is conj^k(c_0) (no twist), so T([y]) := Re(c_0) is the same
for all 12 conjugates, and 12 Re(c_0) is the trace of y to Fp. T is read directly from the coefficients (no
multiplication) and, for distinct classes, is an essentially uniform residue (collisions of T between distinct classes
have probability about 1/p per pair). It plays the role of the class identity for the bits below; K is needed only to
find kappa. `experiment.py` part A checks that Re(c_0) is equal over all 12 conjugates of random GT elements.

**F3 Cyclotomic squaring (Granger-Scott).** GT lies in the cyclotomic subgroup, where
squaring is cheaper than a product. Write y = A + B w + C w^2 with A=(c0,c3), B=(c1,c4),
C=(c2,c5) in Fp4 = Fp2[v]/(v^2-xi), v = w^3, and (A0,A1) = A^2, (B0,B1) = B^2,
(C0,C1) = C^2 computed in Fp4 with 3 Fp2 squarings each. Then
`c0' = 3A0 - 2c0, c3' = 3A1 + 2c3, c1' = 3 xi C1 + 2c1, c4' = 3C0 - 2c4, c2' = 3B0 - 2c2, c5' = 3B1 + 2c5.`
`experiment.py` verifies y^2 = formula on real BLS12-381 GT elements and counts the
operations of the implementation: **9 Fp2 squarings and 34 Fp2 additions/subtractions**
(the multiplications by xi are counted as one Fp2 addition).

**F4 Orders.** r - 1 = 2^32 * 3 * 11 * 19 * 10177 * 125527 * 859267 * 906349^2 * 2508409 *
2529403 * 52437899 * 254760293^2 (verified by multiplication in `experiment.py`). The order
of 2 modulo r is (r-1)/16 (about 2^250.9), so 2^t lies in the order-12 group <p> only if t is a
multiple of about 2^247.3.

## 3. Algorithm

Parameters: R = 1024 table multipliers, squaring probability s = 35/64, theta = 2^-60 (DP
density), M = 2^20 parallel walks. Bits of the class invariant T = Re(c_0) (read as the stored limbs, section 0a): T mod 64 selects
the branch (squaring iff T mod 64 < 35), h = (T >> 6) mod 1024 selects the multiplier, the class
is distinguished iff (T >> 16) mod 2^60 = 0. The fingerprint of a class used for matching reports is
F = (T >> 76) mod 2^128, the 128 bits above the branch, multiplier and DP bits, so the four fields are disjoint and H2 treats them as pseudorandom and independent.

Precompute M_h = g^{a_h} z^{b_h} for h < R (fresh uniform a_h, b_h mod r) and
T[h][t] = M_h^{p^t}, t = 0..11, with the scalar tables a_h p^t, b_h p^t mod r
(12,288 entries of 576 bytes, 7.1 MB).

A walk keeps some member y of the current class together with (alpha, beta) mod r such that
y = g^alpha z^beta. Each walk starts at fresh uniform (alpha, beta). One step:

```
1  T <- Re(c_0) of y (no multiplication); branch, h and DP bits from T
2  if class distinguished: compute the 12 keys (4 Fp mult) -> kappa; append the report
   (F, p^kappa*alpha, p^kappa*beta) as an 80-byte record to one flat buffer
3  2-cycle / checkpoint tests on T (section 3.3)
4  if T mod 64 < 35 (squaring branch):  y <- y^2 (Granger-Scott);  alpha,beta <- 2*alpha, 2*beta
   [no keys, no twist]
   else (multiply branch): compute the 12 keys (4 Fp mult) -> kappa; h <- (T>>6) mod 1024;
        t <- (12 - kappa) mod 12;  y <- y * T[h][t]   (ONE dense Fp12 multiplication);
        alpha += a_h p^t;  beta += b_h p^t
```

**3.1 The step is a map on classes.** Write y^ = y^{p^kappa}. Multiply branch:
y * M_h^{p^{-kappa}} = (y^ * M_h)^{p^{-kappa}} lies in the class [y^ M_h]. Squaring branch:
y^2 = (y^^2)^{p^{-kappa}} lies in [y^^2], because squaring commutes with Frobenius. Both
depend only on the class (T is a class function; kappa and y^ enter only through the conjugation
absorbed by the twist), so walks that reach the same class stay together and reach the same DP.

**3.2 Collision and the report store.** Reports go into one flat buffer of 80-byte records (a 128-bit
fingerprint F and the two 256-bit canonical exponents), in arrival order; the buffer holds at most
D = 2^66.12 records, and a report that would exceed the cap trips the total-work abort with failure
(section 3.4, inside its 0.01 budget). No report is ever replaced or dropped, so a useful class cannot
be shadowed. After the walk budget is spent, the buffer is sorted by F with a fixed number of stable
16-bit radix passes, least-significant digit first: eight passes cover the 128-bit fingerprint, each
pass counts, prefix-sums and scatters a 2^16-entry table of 128-bit counters and permutes a 128-bit
index array between two buffers. A single linear scan then compares the fingerprints of adjacent
records, so every pair of equal-F reports is adjacent; each run's first record is its representative
and the scan verifies it against the run's next two records, stopping at the first verified pair. Two
reports with equal F from different walks are the same class, up to a false match of probability 2^-128
per pair of distinct classes: the false equal-F pairs number at most D^2/2^129 = 2^3.23 in expectation,
so at most 2^16 false equal-F runs are verified (Markov, probability at most 2^-12.8), at most two
verifications each, one modular inverse and g^x = z per verification, <= 600 units. The canonical exponents (alpha^, beta^) = p^kappa (alpha, beta) of a true
pair satisfy alpha^1 + x beta^1 = alpha^2 + x beta^2 mod r, so x = (alpha^1 - alpha^2)/(beta^2 - beta^1).
Two records of the same class solve unless beta^1 = beta^2, probability 1/r, so a clean true run fails
only with probability 1/r^2. The true class's run holds at most 2^24.3 reports (the 20 M repeated
reports cap), so the expected number of foreign records that share its F falsely is at most
2^66.12 * 2^24.3 / 2^128 = 2^-37.58; the scan therefore finds the true pair except with probability at
most 2^-37.5. A run created by a false pair is mixed and is skipped after its two verifications; the
number of such false runs is at most the number of false pairs, 2^3.23 in expectation, below the 2^16
cap. Output Q = [x]Q0 after checking g^x = z.

*Report-store worst case.* The per-report instructions are an itemised sum, not an expectation:

| Work | Instructions |
| --- | ---: |
| 12 keys for kappa (charged again even on a multiply step that already paid them) | 2,119 |
| two 255-bit modular multiplications by p^kappa mod r, at 512 each | 1,024 |
| append the 80-byte record: 10 word stores, the 128-bit address 80*d + base, the 128-bit cap check and counter update | 32 |
| radix sort, eight passes: counting sweep 29 and scatter sweep 42, budgeted 80 per pass | 640 |
| final equal-fingerprint scan of the adjacent records | 16 |
| **itemised 3,831, charged** | **4,000** |

The fixed pass count and the per-pass cost do not depend on the data: every record is appended once,
read and written exactly eight times, and compared with one neighbour, so the worst case per report is
this sum whatever the fingerprints are. The 2^16 counters of each pass are zeroed and prefix-summed
once, 8 x 2^16 entries at 16 instructions = 2^23 instructions in total, at most 2^8.8 units, inside the 2^31
setup bucket. All addresses, offsets and counters are 128-bit, because the peak record address
80 D = 2^72.44 exceeds one 64-bit word; the append and both sweeps price them as such.

**3.3 Fruitless cycles.** A closed walk of the class graph is *useless* if the canonical
exponent pair returns unchanged. After l steps the pair is Lambda (alpha^_0, beta^_0) + tau
where Lambda = 2^t p^m (t squaring steps, m a sum of twists) and tau = sum_h c_h (a_h, b_h)
with c_h in Z[p]. Because (alpha^_0, beta^_0) is uniform and independent of the table, the pair
returns unchanged with probability about 1/r unless Lambda = 1 and every c_h = 0.
* **No squaring step can occur in a useless cycle.** Lambda = 1 with t >= 1 needs
  2^t in <p> (F4), i.e. t >= 2^247.
* A useless cycle made of multiply steps needs every table index used to occur at least twice
  with coefficients p^{e_1} + p^{e_2} + ... = 0 in Z[p]/(Phi_12(p)): for two occurrences
  e_2 = e_1 + 6 (probability 1/12), for three {e, e+4, e+8}, and so on; for any number of occurrences
  the probability is at most 1/12 (the last twist must hit one value). For 2-cycles this gives the
  model rate (29/64)^2/(12 R) <= 1/(12 R) = 2^-13.6 per step (measured 2.58e-4 per step pooled over
  2,200 runs at R = 64 versus 2.67e-4 predicted, section 7). For l = 2k indices drawn independently and
  uniformly (heuristic H2) the probability is at most (2k)!/(2^k k!) (12R)^-k (pairings dominate), which
  is below 2^-288 for 65 <= l <= 147. For l >= 148 all l consecutive steps must be multiply steps, which
  has probability (29/64)^l < 2^-209 (the branch bits are bits of T, H2). Such cycles are never met in
  practice; if one is, the stuck-walk rule restarts the walk after 20/theta = 2^64.3 steps, which adds
  less than 2^-30 relative work at that rate.
* **Detection.** The step before last is remembered; if its T equals the current T the two
  classes form a 2-cycle and the walk replaces the representative of the class with the smaller
  T by its square (y -> y^2, exponents doubled), cost <= 3 units. This depends on the cycle
  only, so merged walks stay merged. A checkpoint T is stored every 64 steps and compared at
  every step; a repeat means a cycle of length <= 64, which is traversed once more to its minimum-T
  class that is then squared (waste <= 200 steps). Every other walk that sees no DP within
  20/theta = 2^64.3 steps is restarted from a fresh uniform pair (stuck-walk rule).
* Total overhead factor on the step count: eps = 5/(12R) + 250/R^2 = 6.5e-4 (charged).

**3.4 Caps and the total-work abort.** The algorithm meters its own charged work with the prices of
section 5 (1.11945 units per multiply step, 0.49895 per squaring step, the stated costs of escapes,
starts, report handling and verifications) and **aborts with failure before executing any next item
if the metered total plus that item's charge would exceed W_cap = 2^125.62**, the submitted bound. Setup, report processing, verification and output obey the same pre-charge check;
unused output work is not free. The claimed bound therefore holds for every run by construction;
what remains is the probability that the abort fires, which is part of the failure budget below.

*Heuristic H1+ (no drift of the walk statistics).* The cap exceeds the nominal ledger
W_nominal = 2^125.6198 by 0.0164% (0.0001285 units per step; `experiment.py` part D), where the parent
of this package, which charged 300/2,200/11,400 for the same listing, had 1.7% and 0.0163 units
per step. A multiply step costs 1.11945 units and a squaring step 0.49895, so the cap is not
exceeded while the fraction of squaring steps stays within 0.0207 percentage points of 35/64. A
larger drop does not break the bound either: it makes the meter halt a corresponding number of
steps early, and the success probability of the shorter walk is still given by the section-4
collision law. A drop of D percentage points raises the average step price by 0.006205 units per
percentage point, delta = 0.006205 D / c, and the meter then buys n/(1 + delta) steps of the nominal n
(part D prints the table from this same formula):

| drop of the squaring fraction | success value | probability spent |
| ---: | ---: | ---: |
| 0.5 point | 0.4280 | 0.0020 |
| 1 point | 0.4255 | 0.0045 |
| 1.5 points | 0.4230 | 0.0070 |
| 2 points (the 0.01 budget) | 0.4205 | 0.0095 |
| 3 points | 0.4156 | 0.0144 |

H1+ states that over the whole run the frequencies of the branch bits, of the 2-cycle events and
of the DPs stay inside a drift of this size — the 0.01 budget is not exhausted until about 2
percentage points, and the production-field measurement of part C is 0.5469 +- 0.004 over five
seeds of 6,000 steps, within 0.8 points of 35/64 — and that the number
of reports, false fingerprint candidates and restarts stays within the caps below, with total
probability of any cap or the abort firing at most 0.01. This is a much weaker statement than
independence of the steps; it excludes only a drift by percentage points over 2^125 steps. Support:
(i) the marginal probability of the squaring branch is 35/64 for every class up to a bias below
2^-374 (T mod 64 is uniform on [0,64) up to the O(64/p) edge effect, F2b), so the expected charge of a
step is the nominal one within that bias and no independence
(ii) in the toy runs the fluctuation of the mean charge of a
run scales as the independent-steps prediction (section 7; 0.394 pooled against 0.392 over 2,200
runs), so there is no excess correlation of
branch labels along the walk at that scale; (iii) the production-field walk of part C has squaring fraction
0.5470 at 6,000 steps here (0.5484 under the integer rule) and 0.5399 to 0.5518 across five
6,000-step seeds (seed 1 to 5 of `experiment.py`), i.e. 0.5469 +- 0.004 measured and 35/64 by
construction; the largest deviation, 0.70 points, is inside the 2-point drift H1+ allows. The toy runs cannot test
concentration at 2^125 steps; H1+ extrapolates it, exactly as H1
extrapolates the collision law. The success analysis (section 4) charges 0.01 for it.

*Caps.* reports D <= 1.1 S theta + 2^25 (S theta = 2^65.98; 2^25 covers the at most 20 M repeated reports of
the merged trails), at most 2^66.12; at most 2^16 false equal-F runs, two verifications each, plus two for the true run
(the expectation of false pairs is at most 2^3.23, so Markov gives probability at most 2^-12.8; the verification
work is at most (2^17 + 2) * 600 < 2^26.23 units); at most
2^50 restarts (expected about 2^37, Markov 2^-13; each restart is a fresh start of about 2^10 units, 2^60 in
total); the total number of steps <= S(1+eps) + M 20/theta. Repeat steps after the first collision are at
most E = 2^84.32 steps (<= 1.11945 units each), negligible against the cap slack 2^113.05. Fruitless cycles
longer than 64 steps need all their steps to be multiply steps (rate below 2^-100 per step, about 2^25.8
events, each wasting at most 20/theta = 2^64.3 steps before the stuck-walk restart: 2^90.1 steps in total).
These caps bound the work and memory stated below.

## 4. Success probability

**H1 (random mapping with branch mix).** The class sequence behaves like iterates of a
random mapping on N classes in which two distinct classes have the same image with
probability (1 - sum_b q_b^2)/N, where q_b are the branch probabilities. Justification: for
c1 != c2 the equation F(c1) = F(c2) has, for each branch pair (b1,b2) with b1 != b2, about one
solution among the N classes (12 candidate representatives, each canonical with probability
1/12); for b1 = b2 = squaring it has none (squaring is a bijection of GT of odd order that
commutes with Frobenius), and for b1 = b2 = a table index it has about 11/12 of one. We
use the conservative value sum q^2 = s^2 + (1-s)^2/R = 0.299273 for all pairs; the exact
value (s^2 + (1-s)^2/(12R)) is smaller by 1.8e-4. Part E of `experiment.py` tests this law directly (section 7). This is the same correction as the Brent-Pollard
factor R/(R-1) of the pure adding walk (s = 0), and for Pollard's original three-branch walk
it is the mechanism that makes Pollard's original square-and-add walk slower than an r-adding walk
(Teske 1998). Hence the effective class count is
N' = N / (1 - sum q^2), log2 N' = 251.7852.

Under H1 the probability that the first S steps contain a class collision is
1 - exp(-S^2/(2N')). Choose P_model = 0.43:
`S = sqrt(2 N' ln(1/0.57)) = 2^125.9771.`
Losses before the first useful collision is detected: detection failure e^-20 = 2e-9, a
collision inside a single walk (fraction about 1/M = 1e-6), beta^1 = beta^2 (2^-255), cap
cap and abort failures (at most 0.01 together, H1+, section 3.4), steps lost to restarts (20 e^-20 = 4.1e-8 of all steps; about 2^37 restarts,
below the cap 2^40, each costing a fresh start of about 2^10 units, 2^47 in total, inside the 2^64.5
bucket of section 5). Their sum is below 0.011. The model value 0.43 therefore leaves 0.03 absolute
slack over the required 0.39 for the inaccuracy of H1. Toy measurements (section 7) of the fraction of
runs with a useful collision within S(0.43) are 0.393-0.443 per setting, pooled 0.419 +- 0.011 over 2,200
runs, at 0.406 for S(0.42); individual settings of a few hundred runs scatter by up to 0.04.

## 5. Cost (units of one dense Fp12 multiplication = 20,000 instructions)

Primitive costs, scalar AArch64, using the measured blst_fp_mul = 383 instructions (commit f7a88067,
harness and counts in `measured-kernel.md`): a
modular add/sub/neg/double of 6 limbs called as a routine, **45** instructions (6 add/sub, 6
conditional-subtract, 6 select, 12 loads, 6 stores, 9 for call, return and address set-up); Fp2
squaring (complex method: (a0+a1)(a0-a1) and 2 a0 a1) = 2 Fp multiplications + 3 additions = 901;
Fp2 addition = 90.

**Table multiply:** the product y*T[h][t] is the dense blst_fp12_mul reference operation,
charged exactly 1.

**Squaring branch:** the straight-line wide-lazy kernel compiled over blst's exported arm64
primitives and measured in the grader's image: 79,527,954 Callgrind `Ir` over 8,192 in-place calls
= **9,708.0022 instructions per call**, rounded upward and charged at **9,709 instructions = 0.48545 unit**.
No fractional instruction remainder is discarded. The per-part split is 9 wide
Fp2 squarings (3,690), 24 wide Fp2 add/subs (1,473), 12 Montgomery reductions (2,448), the narrow
sums, chains and `xi*C1` (1,708), the C body (379) and the loop (10). The control
`blst_fp12_cyclotomic_sqr` measures 9,803.0018 in the same run; the kernel is bit-exact against it
on 1,000 random cyclotomic elements. The primitives are blst f7a88067's `sqr_382x`,
`__add_mod_384x384`, `__sub_mod_384x384`, `redc_mont_384` and `mul_by_1_plus_i_mod_384x`;
`measured-kernel.md` reproduces the harness, the image and the counts. The narrow listing that the
parent priced (9 Fp2 squarings + 34 Fp2 additions + 200 moves = 9*901 + 34*90 + 200 = 11,369) stays
the derivation of the Fp2 item counts used above and on multiply steps. A squaring step needs no key
computation.

**Per-step work.** Every step (both branches):

| Work | Instructions |
| --- | ---: |
| read Re(c_0) (6 limbs), branch bits, DP test, table address | 50 |
| (alpha, beta) update mod r (two 4-limb add + conditional subtract, or doubling) | 100 |
| 2-cycle and checkpoint comparisons on T, history | 60 |
| loop control | 30 |
| margin for loads and moves | 30 |
| **total counted 270, charged** | **270 = 0.0135 units** |

Multiply steps in addition (the keys give kappa):

| Work | Instructions |
| --- | ---: |
| 4 Fp multiplications (a zeta, a Re d_1, b Im d_1, B zeta) | 1,532 |
| 1 Fp addition (B) + 8 Fp subtractions/negations, 45 each | 405 |
| minimum of the 12 keys (12 top-limb loads, 11 x 3 compare/select, tie test, index) | 150 |
| load of the w^1 coefficient, kappa, twist index | 32 |
| **total counted 2,119, charged** | **2,119 = 0.10595 units** |

So a multiply step costs 1 + 0.0135 + 0.10595 = 1.11945 units and a squaring step 0.48545 + 0.0135 = 0.49895 units.
Average step cost c = (29/64) * 1.11945 + (35/64) * 0.49895 = **0.7801140625 units** (= 0.846 x the parent's 0.9225).

Setup and output (all <= 2^31 units; `experiment.py` part D sums them to at most 2^29.9973): the
per-instance pairing g = e(P, Q0) and the verification pairing e(P, Q) (<= 2^14 each; part F counts
10,331 units for each, charging every Fp2 and Fp12 operation of the pinned Miller loop, every product of
the final power by the 4,314-bit exponent (p^12-1)/r, and every Fp product of the Fermat inversion behind
f_num / f_den at one full unit; the two counts are equal because the branch pattern does not depend on the
input: a doubling always takes the tangent branch of the line function, and the running point [k]Q with
1 < k < 2^64 < r is never +-Q, so every addition takes the generic branch), check z^r = 1 (389 units), the
table (1024 multipliers x (2 exponentiations <= 1,021 units + 11 Frobenius maps) <= 2^20.1), M = 2^20
starts g^alpha z^beta (<= 2^30.1), and the output Q = [x]Q0 in G2 at its hard maximum. For the output,
x <= r - 1 < 2^255, so the left-to-right ladder takes at most 254 doublings and 254 additions with the
pinned projective formulas over Fp2: at most 20 Fp2 operations per doubling and at most 28 per addition
(the largest branch; the generic branch is 25), plus 617 for the affine encoding (one Fermat inversion in
Fp with every Fp product counted, and the Fp2 products around it), each operation charged one full unit:
at most 254 x (20 + 28) + 617 = 12,809 units. DP handling: at most 2^66.12 reports (the cap) x 4,000 instructions (section 3.2: the 12 keys for
kappa, 2,119; two 255-bit modular multiplications by p^kappa, 1,024; append, 32; eight radix passes,
640; final scan, 16; itemised 3,831, charged 4,000) / 20,000, at most 2^63.80 units. All non-walk
work is therefore below 2^63.80 + 2^31 < 2^63.81 < 2^64.5 units; the 2^31 covers setup, output
(2^29.9973) and the at most 2^26.23 units of false-candidate verification.

Detection delay: after the first useful collision all M walks continue for 20/theta steps:
E = M * 20 / theta = 2^84.32 steps.

`W = c * S * (1 + eps) + c * E + 2^64.5 = 0.7801140625 * 2^125.9771 * 1.0006453 + 0.7801140625 * 2^84.32 + 2^64.5 = 2^125.6198`.

Submitted bound: **time_log2 = 125.62** (the least hundredth at or above the computed value; the
computed value is 2^125.6197624, so 0.0002 bits (rounded down) above it, a 0.0164% allowance over the whole ledger).
`experiment.py` (part D) recomputes this ledger, and the parent's 2^125.7147 with the parent's own
charge vector at its own s = 3/8, with exact decimal arithmetic in one code path.

Memory: the flat store holds at most 2^66.12 reports x 80 bytes (128-bit fingerprint, two 256-bit
exponents), 2^72.44 bytes, plus two 128-bit index buffers of 2^66.12 entries, 2^71.12 bytes, and 2^16
128-bit counters, 2^20 bytes: at most 2^72.93 bytes in total. With the 7.1 MB table, the 0.8 MB exponent
tables and 2^20 walk states (at most 2 KB each, including the remembered previous state) the peak stays
below 2^72.94 bytes, hence below 2^74.
Memory is reported, not scored.

## 6. Heuristics and limitations

* H1 above (random-mapping behaviour of the mixed walk on classes), H2 (pseudorandom branch and DP bits
  of T), independence of table indices (used for cycle probabilities). These are standard in
  rho analyses (Teske 1998/2001; van Oorschot-Wiener 1999; Duursma-Gaudry-Morain 1999;
  Bos-Kleinjung-Lenstra 2010), supported here by toy and production-field runs, not proved.
* The per-step bookkeeping and key prices follow from the measured fp_mul and a stated primitive model
  (modular add 45 instructions); they are analytic instruction bounds, not measured routines. The
  squaring price is the compiled straight-line kernel measured at 9,708.0022 instructions
  (`measured-kernel.md`: the harness, the grader's image, the per-part split, and the bit-exactness
  check against the 9,803.0018 control in the same run). A pessimistic
  squaring price of 1.0 unit (a
  full dense multiplication) would remove the gain: at this branch mix the ledger would read 2^126.06,
  above the plain walk's 2^125.88; each +0.1 unit on the squaring price adds 0.10 bits, each +100
  instructions of per-step work 0.009 bits.
* Collision statistics are validated in toy fields only. At s = 35/64 the measured mean of the mixed
  walk is within 3.4% of the model in all settings (pooled +0.7% over 2,200 runs); the plain-walk
  controls are -0.5% and -1.2%. The pooled fraction of runs inside S(0.43) is 0.419 against the model's
  0.43, and 0.406 inside S(0.42) against 0.42.
* This is a generic square-root attack. It does not use index-calculus (TNFS) methods in Fp12.
  It is a valid upper bound, not a security estimate.

## 7. Experiments (`experiment.py`, standard library only)

* **Part A (BLS12-381):** GT membership; Frobenius formula = p-th power; 4-multiplication key
  formula = keys of the 12 explicitly computed conjugates, all distinct; kappa picks one identical
  canonical element from every conjugate; p^6 inverts; ord_r(p) = 12; 1+p^4+p^8 = 0; r-1
  factorisation; ord_r(2); the Granger-Scott formula equals y*y with 9 Fp2 squarings and 34 Fp2
  additions.
* **Part C (BLS12-381, production step rule):** a real walk of thousands of steps in GT with the
  production branch rule (branch and multiplier from the bits of Re(c_0), twisted table, squaring,
  2-cycle escape). Checks: every one of the first 150 steps lands in the class of y^ M_h or
  y^^2 as computed with explicit Frobenius; the exponent invariant y = g^alpha z^beta at the
  end; the squaring fraction (35/64 = 0.5469 expected) and chi-square of the low 3 bits and of the multiplier bits of Re(c_0).
* **Part F (BLS12-381, reduction and output):** P = [b]P0 for a sampled b; g = e(P, Q0) by the
  transcribed pinned pairing; g^r = 1 and g != 1; Q = [x]Q0 for a sampled x, which stands in for the
  walk's logarithm (the toy fields of part B have no curve, so the output stage is checked on BLS12-381
  itself); Q on the G2 curve, [r]Q = O and [r]Q0 = O; e(P, Q) = g^x; the identity target returns O with
  e(P, O) = 1; the operation counts of section 5 (pairings, ladder branches, affine encoding).
* **Part E (pair-collision law, toy fields):** sample m distinct classes, apply one step F to each and
  count coincidences of the images; N times the pair-collision rate is the quantity H1 predicts as
  1 - s^2 - (1-s)^2/(12R). Offline results (30,000 sampled classes per row; relative standard error 1.7% to 5.0%):

  | q, R, s | observed | exact prediction | conservative value used | random mapping |
  | --- | ---: | ---: | ---: | ---: |
  | 67, 64, 35/64 | 0.7186 | 0.7007 | 0.6977 | 1 |
  | 67, 8, 35/64 | 0.6854 | 0.6988 | 0.6753 | 1 |
  | 211, 16, 35/64 (4.6%) | 0.7453 | 0.6999 | 0.6881 | 1 |
  | 67, 64, 0 | 0.9827 | 0.9987 | 0.9844 | 1 |
  | 67, 8, 0 | 1.0011 | 0.9896 | 0.8750 | 1 |

  The squaring branch removes its s^2 of the collision probability, and the table branches
  keep 11/12 of theirs (rows with s = 0 and R = 8 separate the exact from the conservative value).
  The largest deviation is the q = 211, R = 16, s = 35/64 row, 1.31 absolute standard errors (1.31
  relative to the prediction, from its 4.64% relative standard error); every observed value is within
  1.4 standard errors of the exact prediction.
* **Part G (report store):** the same eight-pass stable 16-bit radix sort over 128-bit fingerprints,
  followed by the adjacent-record scan, on 4,096 seeded records with eight planted three-record
  duplicate groups; it checks the fixed pass count (8), that the sorted order is ascending, and that
  all 8 planted groups and every dictionary-reference duplicate group are found (`complete`).
* **Part B (toy fields of the same tower shape, q = 67, 139, 163, 211, 379):** the complete
  algorithm (keys with tie-breaking by full conjugates, kappa trick, twisted table, squaring,
  2-cycle and checkpoint escapes, DPs, collision solving, verification g^x = z). Toy keys take
  few values, so the branch/DP bits are taken from a hash of the canonical conjugate (a
  stand-in for the bits of the 381-bit K; Part C tests those bits in the real field). Part B also reports the
  fluctuation of the mean charge per run against the independent-steps prediction (section 3.4).

Offline aggregate of the same script (all discrete logs recovered and verified; `trials` runs per
setting, `dp_bits` 3 and `walks` 4 unless noted; S(x) is the step budget at model success x):

| q, N classes, R, seed | squaring | runs | mean T_c / sqrt N (model) | P[T_c <= S(0.40)] | P[T_c <= S(0.42)] | P[T_c <= S(0.43)] |
| --- | --- | ---: | --- | ---: | ---: | ---: |
| 67, 129,145, 64, seed 11 | 35/64 | 400 | 1.543 (1.500) | 0.367 | 0.383 | 0.393 |
| 211, 721,279, 256, seed 21 | 35/64 | 300 | 1.509 (1.498) | 0.383 | 0.413 | 0.423 |
| 163, 4,524,905, 64, seed 31 | 35/64 | 800 | 1.517 (1.500) | 0.384 | 0.400 | 0.414 |
| 163, 4,524,905, 64, seed 22 | 35/64 | 100 | 1.551 (1.500) | 0.370 | 0.370 | 0.400 |
| 163, 4,524,905, 64, seeds 61-63, walks 8, dp_bits 4 | 35/64 | 600 | 1.474 (1.500) | 0.418 | 0.433 | 0.443 |
| **pooled mixed walk** | 35/64 | 2,200 | 1.510 (1.500) | 0.389 | 0.406 | 0.419 |
| 67, 129,145, 64, seed 12 (control) | 0 | 400 | 1.256 (1.263) | 0.378 | 0.395 | 0.415 |
| 163, 4,524,905, 64, seed 32 (control) | 0 | 600 | 1.248 (1.263) | 0.388 | 0.412 | 0.423 |

The control rows are plain adding walks (their model values use the plain-walk effective class count).
Measured 2-cycle rate at R = 64 with squaring: 2.51e-4 to 2.82e-4 per step across these settings
(predicted (29/64)^2/(12R) = 2.67e-4, bound 1/(12R) = 1.3e-3);
without squaring 1.28e-3 and 1.30e-3 (1.30e-3 predicted). Reproduce a row with e.g.
`echo '{"q":163,"trials":800,"R":64,"sq64":35,"dp_bits":3,"walks":4,"seed":31,"production_checks":0,"seconds":150}' | python3 -I experiment.py`.

References: Pollard 1978; van Oorschot & Wiener, J. Cryptology 12 (1999); Teske, Math. Comp.
70 (2001); Duursma, Gaudry & Morain, ASIACRYPT 1999; Granger & Scott, PKC 2010; Aranha et
al., EUROCRYPT 2011; Bos, Kleinjung & Lenstra, ANTS 2010; mpjunior92, submission 7da5f23a.
