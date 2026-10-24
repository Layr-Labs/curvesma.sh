# Frobenius-class parallel Pollard rho in GT with a cyclotomic-squaring branch (success target 0.39)

Claim: time_log2 = **125.667** (total work W <= 2^125.667 dense-Fp12-multiplication
equivalents, C = 20,000), success probability >= 0.39, peak memory < 2^73.4 bytes
(a hard bound, section 5); preprocessing and table construction at most 2^57 units,
including the 2^56.30-unit log initialization, all included in W; no advice.
Heuristics are listed in section 6.

Design provenance. The reduced-torus class walk, the four-multiplication class-key
formula and the twisted multiplier table are the published Frobenius-class rho design
for this track. The report log, the offline matcher and the per-step instruction ledger
this package extends are the ones accepted on this track at 125.826. The squaring
branch, the trace-invariant decision bits and the branch-adjusted collision law are the
mechanism accepted on this track at 125.74. This package re-derives all of it here and
differs in that the branch is stacked on the accepted report machinery: (1) every
per-step decision (partition, branch, distinguished test, cycle comparisons) is read
from one free class invariant T0, so a squaring step computes no class key at all and
the multiply step's bit reads are recounted; (2) the collision law, the step budget and
the fruitless-cycle rates are re-derived with the squaring branch in place (sections
3.2, 4); (3) because the charged work per step is now a random variable, the walk phase
carries a hard per-slot work meter whose abort rule and success allowance are stated
(section 5); (4) the report log, matcher, pair budget and memory bound are those of the
125.826 package, with the report margin restated for the larger step budget and the
append charge raised to carry the class key computed at a report. No earlier result is
used as evidence; every bound below is derived or computed here.

## 1. Reduction

Input: nonidentity Q in G2 and z in GT. Compute g = e(P0, Q) with the pinned
convention. Q has prime order r and the reduced Tate/ate pairing is nondegenerate,
so g != 1 and g generates GT (prime order r). Every valid z equals g^x for a unique
x in [0, r). Recover x by the rho algorithm below, output P = [x]P0, check
e(P, Q) = z. Then e(P, Q) = g^x = z. If z = 1, output P = infinity directly.
The answer is for the supplied Q; no other Q is used.

## 2. Arithmetic facts used (checked numerically in `experiment.py`, part A)

Represent Fp12 as blst does: Fp2 = Fp[u]/(u^2+1), Fp12 = Fp2[w]/(w^6 - xi),
xi = 1+u, basis w^0..w^5 (blst's v^i w^j with m = 2i+j). This is the profile
field: w^6 = 1+u gives u = w^6 - 1 and w^12 - 2w^6 + 2 = 0, so conversion to the
profile's coefficient list is c_m = a_m - b_m, c_{m+6} = b_m (free).

* Frobenius: (sum c_m w^m)^p = sum conj(c_m) d_m w^m with d_m = xi^{m(p-1)/6}, and
  d_0 = 1.
* p has multiplicative order exactly 12 modulo r (embedding degree 12), so the
  automorphism group {y -> y^{p^k}} of GT has order 12; it contains inversion
  (p^6 = -1 mod r). Every y != 1 lies in a class [y] of exactly 12 elements.
  Number of classes N = (r-1)/12, log2 N = 251.2721.
* **The free class invariant T0.** Write y = sum_{m=0..5} c_m w^m with c_m in Fp2 and
  let T0(y) = Re(c_0) be the Fp part of the w^0 coefficient (in the profile's flat
  coefficient list this is c_0 + c_6, since w^6 = 1+u). Because d_0 = 1 and the
  w^0 coefficient of y^{p^k} is conj^k(c_0), T0 is **equal on all 12 members of a
  class**, and reading it costs one load of the low limbs of c_0 — no field
  multiplication. Part A checks T0(y^{p^k}) = T0(y) for k = 0..11 on random GT
  elements. Outside cycle handling every per-step decision below is a function of T0
  and therefore of the class alone; the two cycle escapes are functions of the cycle
  reached, not of the walker (section 3.1).
* Let zeta = xi^{(p^2-1)/6}. It lies in Fp and is a primitive 6th root of unity
  (zeta^2 = zeta - 1). Writing c_1 = a + b'u for the w^1 coefficient of y, the
  Re(w^1 coefficient) of the conjugate y^{p^k} is
  `key_k = a*zeta^j` for k = 2j and `key_k = B*zeta^j` for k = 2j+1, where
  `B = Re(conj(c_1) d_1) = a*Re(d_1) + b'*Im(d_1)`.
  Since zeta^j in {1, zeta, zeta-1, -1, -zeta, 1-zeta}, the twelve keys are the two
  sign-symmetric triples {+-a, +-a*zeta, +-(a*zeta - a)} and {+-B, +-B*zeta,
  +-(B*zeta - B)}: all 12 keys come from **4 Fp multiplications** (a*zeta, B as two
  products, B*zeta) together with 8 subtractions/negations if every one of them is
  formed. No Frobenius of the whole element is ever computed.
  Note that inversion (k = 6) maps key a to -a, so the 12 keys are pairwise
  distinct for a generic element (checked in part A); there is no built-in tie.
  The minimum over a triple's six residues is min(m, p - M), with m and M the
  minimum and maximum of the triple's three magnitudes (a zero magnitude, probability
  <= 2/p per step, is a degenerate key treated like a tie, i.e. a failure). So K and
  kappa cost the two differences a*zeta - a and B*zeta - B, the two reflections p - M
  and the min/max reductions of the triples, not a twelve-way scan (section 5); part A
  checks that it reproduces the 12-way argmin and key.

Define kappa(y) = argmin_k key_k (keys as canonical residues in [0,p), compared as
integers) and the class key K([y]) = min_k key_k. Then y^{p^{kappa(y)}} is the same
element for all 12 members of a class (canonical representative y^), and K is a
class invariant. Ties need a 1-in-p linear coincidence (probability <= 66/p per
step, total <= 2^-249); they are counted as failure.

* **Cyclotomic squaring.** On GT the Granger-Scott formula squares y with three Fp4
  squarings and the linear maps ret = 3t - 2a on the a[0] half and 3t + 2a on the
  a[1] half (with xi*t2[1] for ret[1][0]). Written in Fp2 primitives that is
  18 Fp multiplications, 64 Fp additions and 31 Fp subtractions, counted by running
  the routine: part A checks that it equals y*y on random GT elements (both the fast
  and the counted version) and prints those operation counts. Squaring commutes with
  Frobenius ((y^{p^k})^2 = (y^2)^{p^k}), so it is a map on classes, and it is
  injective on classes because |GT| = r is odd: y1^2 = y2^2 implies y1 = +-y2, and
  y2 != -y1 unless y1^2 = 1, i.e. y1 = +-1.
* **ord_r(2).** r - 1 factors as 2^32 * 3 * 11 * 19 * 10177 * 125527 * 859267 *
  906349^2 * 2508409 * 2529403 * 52437899 * 254760293^2 (verified in part A), and
  ord_r(2) = (r-1)/16 = 2^250.857 (computed in part A by removing factors). The
  smallest number of squaring steps that can occur in a fruitless cycle is therefore
  the smallest t with ord_r(2) | 12t, i.e. t = ord_r(2)/gcd(ord_r(2),12) = 2^247.27
  (part A prints it). Section 3.2 uses this.

## 3. Walk (R-adding walk on classes with a squaring branch)

Parameters: R = 1024 partitions, distinguished-point fraction theta = 2^-60,
M = 2^20 parallel walks, each started at g^alpha z^beta with fresh uniform alpha, beta mod r.
Squaring-branch probability s = 3/8, from a fixed 3-bit field of T0.

Precompute M_h = g^{a_h} z^{b_h} (h < R, random a_h, b_h) and the table
T[h][k] = M_h^{p^k}, k = 0..11 (12,288 entries of 576 bytes, 7.1 MB), with exponent
tables a_h p^k, b_h p^k mod r.

**The four decision fields of T0** (disjoint bit ranges of one Fp element, read once
per step; T0 has 381 bits, so they fit with 308 bits to spare):

| bits of T0 | width | decision |
| --- | ---: | --- |
| 0 .. 9 | 10 | partition h = T0 mod R (1024 values) |
| 10 .. 12 | 3 | branch: square iff the field's value is below 3 (probability 3/8) |
| 13 .. 72 | 60 | distinguished iff all 60 bits are zero (probability 2^-60) |
| 73 .. 380 | 308 | unused |

The three fields are disjoint, so h, the branch and the distinguished flag are
functions of independent bits of one class invariant, and the marginal probability of
each is exact: 3/8 for the branch, 2^-60 for the distinguished flag. Part A2 measures
the frequency of each of the 73 used bits on real BLS12-381 GT elements (section 7).

Class map: F([y]) = [ y^ * M_h ] on a multiply step, F([y]) = [ y^2 ] on a squaring
step, with h and the branch read from T0([y]). A walk keeps an arbitrary member y of
the current class with its exponents (alpha, beta), y = g^alpha z^beta. Because
y^{p^kappa} * M_h = (y * M_h^{p^{-kappa}})^{p^kappa}, one multiply step is

```
read T0 (one load) and its three fields       (12 instructions, section 5)
if the class two steps back is this class (2-cycle, T0 compare):
    escape: square the stored representative of the smaller-T0 class of the pair,
    drop both key histories, continue        (1 or 2 steps, charged at the step price)
if the checkpoint class (recorded every 64 steps) is this class (T0 compare):
    escape: traverse the cycle once to its smallest-T0 class, square it, continue
if the distinguished field of T0 is all zero:
    compute the class key K (4 Fp mul) and append the 80-byte record
    (bits 70..197 of K, p^kappa*alpha, p^kappa*beta) unless it equals this walk's last
    record; increment this walk slot's cumulative report counter r; if r = R_w: stop
if (squaring branch):                        [3/8 of steps]
    y = y^2                                  (Granger-Scott, section 5)
    alpha = 2 alpha ; beta = 2 beta (mod r)
else:                                        [5/8 of steps]
    compute class keys, kappa, K             (4 Fp mul + 407 other instructions)
    h = T0 mod R ; t = (12 - kappa) mod 12
    y = y * T[h][t]                          (ONE dense Fp12 multiplication)
    alpha += a_h p^t ; beta += b_h p^t (mod r)
```

So a multiply step costs one reference operation plus the key work, and a squaring
step costs the Granger-Scott squaring (0.585 reference operations) with no key work.
Outside cycle handling the step is a function of the current class: the partition h,
the branch and the distinguished flag are read from the class invariant T0, and the
multiplier table is fixed, so two walks at the same class take the same step. The two
cycle escapes are the only history-dependent actions, and each is a function of the
cycle reached rather than of the walker (section 3.1). So two walks that reach the same
cycle compute the same escape class and the same escape action; their class sequences
coincide from their first common class onward except possibly for the iteration on which
each detects the cycle, and they coalesce within one cycle length of the later
detection. This is eventual coalescence, not a memoryless merge: equal current classes
alone do not imply equal next classes while the stored history differs, and section 3.1
states exactly what a class collision then still gives.

**Cycle comparisons compare T0.** On a squaring step no class key exists, so the 2-cycle
comparison, the checkpoint comparison and the escape's ordering all use T0 (plus nothing
else: 381 bits are enough). A false match needs two distinct classes with equal T0,
probability about 1/p = 2^-380.70 per comparison (log2 p = 380.7007 for the 381-bit
modulus); over at most 2 x 2^126 = 2^127 comparisons (2 per step, for fewer than 2^126
steps in all: the walk budget S = 2^125.780073 plus the detection window 2^84.33 and the
restart waste 2^95.66) the total is below 2^127/p = 2^-253.70, i.e. below 2^-253
(a loss term of section 4). Its effect is only a spurious escape
(a wasted step) and never a false logarithm: the report key is derived from K, not from
T0, so two classes that share T0 but differ in K are never matched, and the escape's
choice of class is a function of the T0 multiset of the cycle (section 3.1), so a
spurious escape only wastes steps and cannot desynchronise two walkers that reach the
same cycle.

**Report log (deterministic capacity, direct addressing, per-slot cap).** Reports are
appended to a log: an array of exactly C_rep = M R_w = 2^20 x 2^47 = 2^67 record slots
of 80 bytes each, allocated before the walks start and filled with one fixed sentinel
record (all bytes 0xFF, whose exponent fields exceed r, so the matcher skips it; 12
instructions per record, charged in section 5). Walk slot i owns the R_w = 2^47
consecutive records that start at record 2^47 i, and the slot's cumulative report counter
r addresses its next record directly: the record goes to address base + 80 (2^47 i + r),
a 128-bit offset (two 64-bit words, since the offset reaches 80 x 2^67, below 2^73.321929,
which does not fit in one) formed by a shift, an add and a multiplication by 80. An append is one 80-byte
compare, one counter compare, the offset computation, one 80-byte store and one
increment, a fixed instruction count with no hashing, probing, chaining or eviction, and
no shared counter: no two slots write the same record, so the append needs no atomic
operation. It is charged inside the append's instruction count (section 5). The total
number of appends is at most M R_w by construction: no restart resets a slot's counter r
and the slot is stopped when r reaches R_w = 2^47, so a slot never writes outside its own
region, the log can never overflow and no report-count probability enters the memory or
processing bound. The slot's report counter is one 64-bit word; 2^47 < 2^63, so its cap
test can fire. (The toy of part B keeps each log as a list of the filled records, which
is the same set of records the matcher sees.)

**Report margin for the larger budget.** Under H1 the number of distinguished classes a
slot visits is geometric with parameter theta: its mean over the S_slot + E_slot =
2^105.780073 + 2^64.33 walk steps is 2^45.7801 reports, so the cap R_w = 2^47 is
a = 2^1.2199 = 2.3293 times the mean and the Bernoulli-KL (Chernoff) bound gives
log2 P(slot reaches the cap) <= -mu D(a||1)/ln 2 = -5.58 x 10^13. That is a heuristic
consequence of H1 (distinguished status is a class property, so the visits are not
independent), so it is not used as a hard tail: it is carried in the loss list of
section 4 as the same order of magnitude the accepted 125.826 package carried, and the
capability itself (log <= M R_w) is deterministic and does not rest on it. The append
rule -- a report equal to this walk's last stored record is not stored, and a report is
only stored after both cycle checks -- is what stops a fruitless return from filling
the log with duplicates of one record.

**Matching by an in-place heapsort of the log (worst-case bound).** The walk phase
ends after its work-metered step budget (section 5). Then the n <= C_rep = 2^67
records of the log are sorted in place by the WHOLE 80-byte record
(the 16-byte key T first, then A, then B) with heapsort: Floyd's heap construction,
then n - 1 extractions, using one 80-byte temporary and no other memory. Sorting the
whole record, not the key alone, puts identical records adjacent: it is what makes the
duplicate removal below a single adjacent comparison instead of a search inside an
equal-key run. Its work does not depend on the data: construction performs at most n
sift-down levels, and the extractions at most sum_{k=2..n} floor(log2 k) <= 65 n levels
(exactly 65 x 2^67 + 69 for n = 2^67), so at most 67 n <= 67 x 2^67 levels in all. A
level compares the two children (two 10-word compares: 20 loads, 10 compares, 10
conditional selects), compares the larger child with the parent (10 loads, 5 compares,
selects), moves one 80-byte record on a swap (10 loads, 10 stores), and does index
arithmetic and loop control (50 instructions): 130 in all, so <= 8,710 instructions per
record, budgeted at 10,000. Duplicate removal
is then exactly one adjacent 80-byte comparison per record (<= 100 instructions), with
no search, no bucket and no probe.

**Equal-key runs and the hard pair budget.** After that dedup the log is a sequence of
maximal runs of records sharing the 16-byte key T with distinct (A, B) inside a run.
Such a run of d >= 2 records is the only place pairwise work occurs: every one of its
d(d-1)/2 pairs is enumerated, and each pair is charged <= 200 instructions (the
576-byte element comparison plus the canonical recomputation check), with one
recomputation of the canonical element g^A z^B (<= 1,021 units) for each record a
charged pair touches, so at most 2 x (number of pairs) recomputations.

The pair count is bounded by a HARD budget P = 2^41: the matcher keeps a pair counter
over the whole log, and when it would charge pair P + 1 it aborts the match and the
attempt is a failure. This is a deterministic count, not an expectation. Its failure
probability is bounded by Markov on the number of false-key pairs: under H1 (bits
70..197 of K of distinct classes independent and uniform) two records of distinct
classes agree on T with probability <= 2^-128, the log holds at most 2^67 records, so
there are at most (2^67)^2/2 = 2^133 pairs of distinct records and the expected number
of false-key pairs is <= 2^133 x 2^-128 = 2^5. Hence P(pair count above 2^41) <=
2^5/2^41 = 2^-36, which is a loss term of section 4.

There is no separate per-record confirmation count: the pairs are the whole charge. A
GENUINE repeat -- the same canonical element reached with different exponent pairs, which
is the useful collision -- is found inside such a run and is charged as pairs like any
other: it is not exempt, and no record of a run is exempt. A run whose dedup leaves one
record (in particular an all-identical run) has no pair, so it charges nothing and is not
recomputed. What the rule buys is that the matcher's work is bounded by the pair budget
whatever the mix of false key matches and repeats: the pairs are the only quadratic
term, and they are capped at P.

Collision: the first charged pair whose canonical elements are equal with different
exponent pairs is solved as p^k1(alpha1 + beta1 x) = p^k2(alpha2 + beta2 x) mod r, that
is x = (A1 - A2)/(B2 - B1) mod r (one modular inversion, <= 2^14 instructions). Section
3.1(iii) fixes what "different exponent pairs" covers: a genuine class collision is
always usable, a byte-identical pair is removed by the dedup before pairing, and only a
false-key pair whose B values coincide (a 1/r term of the loss list) gives no x from
that pair, in which case the matcher continues to the next pair. Verify g^x = z
(<= 510 units) before output.

**Cycle rules (deterministic, cycle-determined escapes).** Each step compares T0 with
the T0 two steps earlier; on equality the two classes {C1, C2} form a 2-cycle and the
walk squares the stored representative of the one with the smaller T0 (y -> y^2 maps a
class to a class). A checkpoint T0 is stored every 64 steps and compared at every step;
a repeat means a cycle of length L <= 64, which is traversed once more (L steps) to find
its smallest-T0 class, which is squared. Both escapes depend only on the multiset of T0
values of the cycle's classes, not on where or when the walk entered it, so two walks
that reach the same cycle escape at the same class by the same action; only the
detection iteration can differ, and section 3.1 proves they coalesce within one cycle
length. Detection of an L-cycle (L <= 64) happens
within 64 + L steps of entering it, so the waste per event is <= 64 + 2L + 1 <= 193 <=
200 steps; a 2-cycle event wastes <= 3 steps. An escape is one squaring step (cheap) or
one squaring plus one re-step, and section 5 charges every one of those steps at the
multiply-step price, i.e. above the squaring price.

**Stuck-walk rule.** Each walk counts steps since its last distinguished class in a
128-bit counter (reset at each report and at restart; this is the trajectory counter, not
the slot's report counter r, which no restart resets; 40/theta = 2^65.33 exceeds 2^64). When
the count exceeds 40/theta = 2^65.33 the walk is abandoned and restarted at a fresh
g^alpha z^beta. Its steps are real steps and are inside the slot's work-meter cap of
section 5. It is the safety net for cycles longer than 64 that contain no distinguished
class (section 3.2 bounds the probability that cycles matter to success). Under H1 the
steps between distinguished classes are geometric with parameter theta, so a run exceeds
40/theta with probability (1 - theta)^(40/theta) < e^-40; over at most 2^66 runs the
expected number of restarts is < 2^66 e^-40 < 2^8.3, and by Markov's inequality more than
2^30 restarts occur with probability < 2^-21. 2^30 restarts waste at most
2^30 x 2^65.33 = 2^95.33 steps, a fraction < 2^-30 of S.

### 3.1 Eventual coalescence (what a class collision still gives)

The walk is **not** memoryless. The two-step T0 ring and the checkpoint slot hold stored
values, so a walk's next step can depend on how it reached its current class, and equal
current classes do not by themselves imply equal next classes. What survives is the
weaker invariant below, and it is what the matcher needs.

**(i) Outside cycle handling the step is a function of the current class.** The partition
h = T0 mod R, the branch and the distinguished flag are bit fields of the class
invariant T0, the multiplier table T[h][k] is fixed, and the twist t = (12 - kappa) mod 12
is a class invariant because kappa is the argmin over the class keys. Two walkers at the
same class therefore take the same ordinary step, however they arrived.

**(ii) The escape actions are functions of the cycle, not of the walker.** The 2-cycle
escape squares the stored representative of the class with the smaller T0 among the two
classes of the pair; the checkpoint escape traverses the cycle once to its smallest-T0
class and squares it. Which class is squared therefore depends only on the T0 values of
the cycle's classes, so the escape class C_min^2 is the same for every walker that
reaches that cycle, and squaring any representative of that class lands in C_min^2 (squaring
is a map on classes), so the escape action is the same at the level of classes. Two
walkers that reach the same cycle compute the same escape class and the same escape
action; their class sequences coincide from their first common class
onward except possibly for the iteration on which each detects the cycle, and they
coalesce within one cycle length of the later one's detection: a 2-cycle is detected
within two steps of entering it and an L <= 64 cycle within 64 + L steps, and after the
escape both walkers have an empty T0 ring and an empty checkpoint slot, so by (i) their
subsequent steps are functions of the class again until a later checkpoint detection,
which can differ only in timing by the same argument. Part P of the experiment measures
exactly this: two walkers at one class with different histories, one holding an
all-multiply 2-cycle in its ring and one freshly started at a different conjugate, escape
at the same class within one cycle length (2 steps) and share their class sequence from
there. The walks are not merged in the memoryless sense, and the proof no longer claims
that they are; what matters is that a collision still gives the usable pair (iii).

**(iii) A class collision is still a usable pair.** A report at class C stores
(T, A, B), where T is the truncated class key K(C) and (A, B) = (p^kappa alpha,
p^kappa beta) are the exponents of that walker's element. Two walkers that meet at C
produce the same T. If their exponent pairs differ, the pair solves
p^k1(alpha1 + beta1 x) = p^k2(alpha2 + beta2 x), so x = (A1 - A2)/(B2 - B1) mod r, as in
section 3. If the pairs coincide (A1 = A2 and B1 = B2), the two records are byte-identical:
then alpha1 = alpha2 and beta1 = beta2 mod r, so both walkers hold the same group element
g^alpha z^beta at class C, and by (i) and (ii) their class sequences agree from C onward
(the history rings may differ for at most one cycle length). A pair with coincident
records is not useful, but it cannot hurt: the full-record heapsort puts identical
records adjacent and the adjacent full-record dedup removes the duplicate before any
pair is enumerated, and the append rule separately stores no repeat of a walker's own
last record. So a division by zero from A1 = A2, B1 = B2 cannot arise. Note that for two
records in one class, B1 = B2 holds only in that coincident case: a Frobenius shift
multiplies B by p^k, which is 1 mod r only for k = 0. The remaining unsolvable case is a
false-key pair (distinct classes whose truncated keys agree) that also has B1 = B2; it
gives no x from that pair, is charged as the 1/r term of the loss list of section 4, and
the matcher simply continues with the next pair of the run.

### 3.2 Fruitless-cycle rates (derivation), with the squaring branch

Let a walk hold y_0 and take L steps, of which t are squarings with multiplier 2^t and
the rest have indices h_1..h_m and twists t_1..t_m. Then
y_L = y_0^{2^t} prod_i M_{h_i}^{p^{t_i}}. The walk is back in the class of y_0 with the
same exponent pair (a fruitless cycle, i.e. one that occurs for every x, not by the
1/r chance of a random mapping) iff y_L = y_0^{p^e} as formal expressions in the
independent uniform exponents (alpha_0, beta_0, a_h, b_h). Comparing coefficients
forces 2^t = p^e and, for every index value h that occurs,
sum_{i: h_i = h} p^{t_i} = 0 mod r.

* **No squaring step lies in a fruitless cycle.** 2^t = p^e mod r needs
  2^{12t} = 1 mod r, i.e. ord_r(2) = 2^250.857 divides 12t, i.e. t >= 2^247.27
  (part A computes both numbers). A cycle has at most 2^126 steps in all, so
  t = 0: every fruitless cycle is all-multiply, and the branch removes exactly its
  (1-s)^L share of them.
* Two consequences, as in the base derivation: a single occurrence of an index is
  impossible (p^t != 0), so every index that occurs appears >= 2 times and the L
  steps use d <= floor(L/2) distinct indices; and for a block of k >= 2 occurrences,
  sum p^{t_i} = 0 mod r is the event that a sum of k 12th roots of unity vanishes in
  Z[x]/Phi_12. With uniform twists this is the return probability of a symmetric
  random walk; it is 1/12 at k = 2 and, by P(S_{2m} = 0) = ||mu^{*m}||_2^2
  (non-increasing in m) and Cauchy-Schwarz for odd k, at most 1/12 for every k >= 2
  (part C computes k = 2..10 exactly: max 0.083333 = 1/12).

Under H1 (indices and twists of successive distinct classes behave as independent
uniform values, the same assumption that gives the 2-cycle rate), the probability
that an all-multiply walk enters a fruitless L-cycle at a given step is

  P_L <= sum_{d=1}^{min(R, floor(L/2))} C(R, d) (d/R)^L 12^-d,

and with the branch the same event also requires all L steps to be multiply steps,
probability (1-s)^L, so the mixed rate is P_L (1-s)^L. Part C computes
sum_{L=3..64} P_L (1-s)^L = 0.0289 R^-2 for R = 1024 (against 0.1391 R^-2 without the
branch), and sum_{L=65..2R} P_L (1-s)^L <= 2^-281.6. For L > 2R every term with
d < R shrinks by d/R per extra step, so the tail is negligible as in the base
derivation. The cost ledger nonetheless charges the **unbranched** R^-2 per step with
200 steps per event (section 5), which is more than 34 times the mixed rate and more than
7 times the unbranched one: the squaring branch can
only reduce fruitless-cycle work, and the charge does not depend on that reduction.

**Class cycles with a squaring step are useful, not fruitless.** If a cycle of length L
contains a squaring step, then after one lap the class returns but the exponent pair is
multiplied by 2^t p^m with 2^t != p^e, so the pair is different: two records with the
same key T and different (A, B), which is exactly a solved instance. This has three
consequences. (i) If such a cycle contains a distinguished class, the matcher solves the
instance from those records. (ii) If it contains none, the walk laps without reporting
and the stuck rule restarts it (waste <= 2^65.33 steps, inside the restart accounting
above). (iii) The two cycle checks escape every cycle of length L <= 64, and a 2-cycle
is escaped before it can report, so only cycles with L > 64 can lap. A lapping cycle can
hurt in exactly one place: an equal-key run in the report log, i.e. the hard pair budget
P = 2^41 of section 3, because a lap that passes through a distinguished class stores one
more record of that class, with a different exponent pair each time. A lapping cycle
presupposes a within-walk return, a walk reaching a class it has already visited, so
charge that event itself rather than one consequence of it. Without any within-walk
return a walk stores at most one record per class, so a class holds at most M = 2^20
records (one per walk) and even the whole run gives at most C(M, 2) = 2^19 (2^20 - 1)
< 2^39 pairs at one key, inside P with room for the about 2^5 false-key pairs that the
2^-36 Markov charge of section 3 already covers. By a union bound over the step pairs of
all walks, a return to a visited class has probability at most
M S_slot^2/(2N) = 2^(20 + 2 x 105.780073 - 1 - 251.2721) = 2^-20.71, and that is the
loss term charged (section 4): it covers every consequence of lapping at once, because
the report cap R_w and the slot work meter are hard bounds and no other probabilistic
quantity enters from a cycle. A distance-limited version is consistent -- a return at
distance <= 2^85.3 has probability at most S x 2^85.3/N = 2^(125.780073 + 85.3 -
251.2721) = 2^-40.19 -- but 2^-20.71 holds at every distance and is what the loss list
carries. The accepted 125.826 package carried the same exposure inside its pair-budget
term; here it is charged as a stated loss term, because the pair budget is global over
the log while a lap adds records per key, so no per-key lap count bounds it.

## 4. Success probability

Heuristic (H1): the class sequence of the mixed walk behaves like a random mapping on
N classes in which only a (1-s) fraction of the steps can merge two walks, with the
Brent-Pollard correction for R-adding walks. Two steps merge two distinct classes when
both are multiply steps with different multiplier indices (probability
1 - s^2 - (1-s)^2/(12R) per pair of steps, section 4.1), so the effective class count is

  N' = N / (1 - s^2 - (1-s)^2/R) = 2^251.491408

where the charge uses 1/R instead of 1/(12R) for the same-index case: same-partition
pairs are charged as never merging even when their twist index differs, which only
enlarges N' (conservatively by 4.07 x 10^-4 relative, 5.87 x 10^-4 in log2, since the two
factors differ by (1-s)^2 (1/R - 1/(12R)) = 3.49 x 10^-4). Then the probability that no useful collision occurs among the first S
distinct-step points of all walks is exp(-S^2/(2N')). Choose

  S0' = sqrt(2 N' ln(1/0.61)) = 2^125.737428,   S = 1.03 S0' = 2^125.780073,

where the 3% is explicit slack for walk non-randomness (the ceil of the per-slot budget
adds another 1e-9 relative, 1.4e-9 bits); at S0' the model gives exactly 1 - 0.61 =
0.39, and at S it gives 1 - 0.61^(1.03^2) = 1 - 0.61^1.0609 >= 0.40808 (rounded down). The algorithm
stops after S walk steps (plus the detection window below); S is a fixed budget, not an
expectation.

### 4.1 The branch-adjusted collision law

Two steps of two distinct classes produce the same image class when
C1^{b1} M_1 = (C2^{b2} M_2)^{p^k} for some k, where b_i = 2 on a squaring step (M_i = 1)
and b_i = 1 with a table multiplier otherwise.

* **squaring-squaring pairs never merge.** Squaring is injective on classes (|GT| is
  odd), so C1^2 = C2^2 forces C1 = C2.
* **squaring-multiply pairs merge with probability about 1/N**: C1^2 = C2 M is one
  equation between two independent uniform classes.
* **multiply-multiply pairs merge with probability about (1 - 1/(12R))/N**: the 12R
  multiplier indices (h, t) are a 12R-element set, and equal indices (probability
  1/(12R)) force C1 = C2.

So the per-pair merge probability is (1 - s^2 - (1-s)^2/(12R))/N: exactly the factor
above. Part E of the experiment measures this law directly on toy class sets (distinct
classes, one step each, coincidence counts by pair type), and part B measures the
resulting success fraction.

**Detection.** After the first useful collision, the two walkers coalesce within one cycle
length (section 3.1) and from their first common class onward share a class sequence, so
the wait to a distinguished class is that of one trajectory: within 20/theta further
steps except with probability e^-20 (H1). The window is applied to a trail selected
through a collision, not to a fresh trail; the argument does not need the fresh-trail
form, because under H1 the distinguished indicator of each step is independent of the
past with probability theta, so the residual wait of the selected trail is geometric with
the same parameter (the geometric tail of a fixed-parameter wait) and its tail is
e^-20 as well. No length-biased trail allowance is claimed. Charging all M walks in
lockstep over that window costs E = M*20/theta = 2^84.33 steps.
The losses, each a probability bound: detection tail e^-20 < 2^-28; trailing walk
restarted before reaching the distinguished class (its distinguished-free run must
exceed 40/theta: under H1 the run containing a given step exceeds 40/theta with
probability (1 + 40) e^-40 < 2^-52); self-collisions of a single walk (fraction
<= 2^-19 of collisions, with M = 2^20 walks); a false-key pair with B1 = B2 (1/r); key ties (< 2^-249);
long fruitless cycles (< 2^-100, now only all-multiply and at the mixed rate of section
3.2, so the bound only improves); more than 2^30 stuck-walk restarts (< 2^-21);
more false-key pairs than the hard pair budget (< 2^-36 by Markov, section 3); a
within-walk class return, the only way a lapping cycle can reach the pair budget
(< 2^-20.71, section 3.2); a cycle-handling step count above its budget (probability 0:
cycle handling
is inside the slot's hard work meter, whose abort is the 0.01 allowance); a walk
slot stopped by the per-slot report cap, i.e. a healthy slot reaching R_w = 2^47 reports
(the KL figure of section 3, < 2^-300); two distinct classes with equal T0 causing a
spurious escape (< 2^-253). Their sum is below 2.98 x 10^-6
< 2^-18 < 0.4080 - 0.39 - 0.01 = 0.0080.

**The work-meter allowance.** Total charged work is a random variable (a squaring step
costs 0.585 units, a multiply step 1.1125), so the walk phase is metered and a slot is
stopped at W_cap = 1.01 x its nominal charged work (section 5). The abort is charged as
a stated heuristic allowance of 0.01 of success, with no tail-probability claim: a slot
stopped at 1.01 of its nominal work is the same as running the phase with
1.03/1.01 = 1.0198 times the 0.39 budget, whose model success is
1 - 0.61^((1.03/1.01)^2) >= 0.40194, a model loss of 0.0062. The 0.01 charged is 1.6 times
that, and it also covers the steps a stopped slot never takes.

Success >= 0.40808 - 0.01 - 2^-18 > 0.40808 - 0.01 - 0.00001 = 0.39807, so
success >= 0.3980 >= 0.39 (every figure rounded down), with 0.0080 of margin.

## 5. Cost (units of one dense Fp12 multiplication = 20,000 instructions)

Per-step non-reference instructions, scalar AArch64, using blst's measured
blst_fp_mul = 383 instructions. A 384-bit (6-limb) operation is counted from its
actual instruction sequence: an add or subtraction mod p is 6 ADD/ADC (or SUBS/SBCS)
plus 6 SUBS/SBCS plus 6 CSEL inside 18 instructions of frame, PAC hints, modulus loads
and stores, charged 45 (add) and 44 (subtract); a 6-limb compare with min, max and
index selection is 6 SUBS/SBCS plus 12 CSEL plus 2 conditional moves, charged 24;
p - v is 6 loads plus 6 SUBS/SBCS, charged 20; a 4-limb addition (or doubling) mod r
is 4 ADD/ADC plus 4 SUBS/SBCS plus 4 CSEL plus its load and store, charged 26; a
128-bit counter increment is LDP, ADDS, ADC, STP (4) and a 128-bit compare with a
constant is CMP, CCMP, B.cond (3). The measured 383 for blst_fp_mul already includes
its call, loop and reduction work ("executed instructions, including call and loop
work"), so no separate call charge is added anywhere; the add and subtract prices above
include their own frame and loads.

**Multiply step** (the accepted ledger, recounted for the new bit reads, the meter and
the T0 comparisons):

| Work | Instructions |
| --- | ---: |
| 4 Fp multiplications (a*zeta, two products for B, B*zeta) | 1,532 |
| 1 Fp addition mod p (B = a*Re(d_1) + b'*Im(d_1)) | 45 |
| 2 Fp subtractions mod p (a*zeta - a, B*zeta - B) | 88 |
| min3 and max3 of the two triples (compare, csel, index) | 144 |
| p - max_a, p - max_B (the negated minima of the two triples) | 40 |
| final minimum of the four candidates, kappa = 2j + {0,1} + 6*negated | 74 |
| load the w^1 coefficient (12 limbs) | 16 |
| **common rows (both step types):** read T0's low limbs and its three bit fields | 12 |
| common: 2-cycle T0 comparison | 18 |
| common: checkpoint T0 comparison | 15 |
| common: T0 history ring and the 64-step checkpoint test | 8 |
| common: slot work-meter add and 128-bit cap test | 7 |
| common: steps-since-report counter and stuck-walk test | 7 |
| common: slot step counter (128-bit) | 4 |
| common: loop control, output buffer, slot pointer | 30 |
| h = T0 mod R, t = (12 - kappa) mod 12, table address | 10 |
| alpha, beta += table exponents mod r (2 x 4-limb add, conditional subtract) | 52 |
| **total counted (budgeted 2,250)** | **2,102** |

The budget of 2,250 instructions is 7.0% above the count and is the number charged per
step; the twelve keys are two sign-symmetric triples, so the canonical minimum needs the
two triples' minima and maxima and their negations rather than a scan over twelve keys.
The product y*T[h][t] is the dense blst_fp12_mul reference operation (operands read
through pointers, output written to a rotating buffer; no copies), charged 1.
Multiply-step cost c_mul = 1 + 2,250/20,000 = 1.1125 units.

**Squaring step** (Granger-Scott, counted by running the routine in part A):

| Work | Instructions |
| --- | ---: |
| 18 Fp multiplications (9 Fp2 squarings of 2 multiplications each) | 6,894 |
| 64 Fp additions mod p | 2,880 |
| 31 Fp subtractions mod p | 1,364 |
| alpha, beta doubling mod r (2 x 4-limb doubling, conditional subtract) | 52 |
| common rows of the table above (T0 read, bit fields, cycle comparisons, meter, counters, control) | 101 |
| **total counted (budgeted 11,700)** | **11,291** |

The charge of 11,700 instructions is the accepted figure for a cyclotomic squaring step
(11,400 for the squaring itself plus 300 of bookkeeping), and our own count of the same
formula with the same instruction prices is 409 instructions below it, so the budget is
charged. A squaring step computes no class key: the four Fp multiplications of the key
appear only on multiply steps and at reports. Squaring-step cost
c_sq = 11,700/20,000 = 0.585 units. Average step cost, with s = 3/8, is
c = (5/8) c_mul + (3/8) c_sq = 0.9146875 units (18,293.75 instructions).

**Work meter (a hard cap on a random per-step cost).** Because a step costs 1.1125 or
0.585 units, the walk phase's work is a random variable, so it is metered and capped
rather than left as an expectation.

* Layout. Each of the M = 2^20 slots owns a private meter in **instructions**, held as
  a 128-bit pair of 64-bit words (low word, high word) in the slot record. Its update is
  LDP, ADDS, ADC, STP (4 instructions) and its cap test is a 128-bit compare with the
  constant W_cap (CMP, CCMP, B.cond: 3), charged together as 7 instructions per step in
  the common rows above. The meter is private to the slot, which one processor owns for
  the slot's lifetime, so there is no synchronisation, no atomic and no cache-line
  contention on it. The report log is split into per-slot regions (section 3), so the
  walk phase has no shared counter and no contended structure at all.
* Values. W_cap = ceil(1.01 x nominal), nominal = c x C x (S_slot + E_slot) +
  I_MUL x ceil(OVER x S_slot) = 2^119.939899 instructions, and
  W_cap = 1287740786101123107786939607107502080 = 2^119.954254.
  The meter's largest value is one over-cap update past the test, W_cap - 1 + I_MUL =
  2^119.954254 (log2 119.954254), which fits two 64-bit words with the cap test at
  2^119.954254 < 2^127. A slot is stopped when its meter reaches W_cap, so
  W_walk <= M x (W_cap - 1 + I_MUL) = 2^139.954254 instructions = 2^125.666542 units is a
  **hard bound**, not an expectation; it is what the claim is computed from.
* Abort allowance. The abort is charged as the 0.01 success allowance of section 4, with
  no tail claim.

**Counters and indices.** Every counter and index, its largest value at the end of the
attack (including one over-cap update per processor), and its width. Part C recomputes
this table and asserts every row.

| counter or index | largest value (log2) | words | capacity test |
| --- | ---: | ---: | --- |
| slot work meter (instructions) | 119.954254 | 2 x 64 | >= W_cap with W_cap = 2^119.954254 < 2^127 |
| slot step counter | 106.440033 | 2 x 64 | none (informational) |
| steps since last report | 65.321929 | 2 x 64 | > 40/theta (2^65.321929 rounded up), a 128-bit compare |
| slot report counter r | 47 | 1 x 64 | >= 2^47 with 2^47 < 2^63 |
| log offset 80 x (2^47 slot + r) | 73.321929 | 2 x 64 | < 2^128 |
| slot restart counter | 41.118105 | 1 x 64 | < 2^63 |
| matcher pair counter | 41 | 1 x 64 | > 2^41 with 2^41 + 1 < 2^63 |
| heapsort index | 67 | 2 x 64 | < 2^128 |

No counter's capacity test sits at its word's own capacity: the two cap tests that can
fire are at 2^119.954254 of a 2^128 word and at 2^47 of a 2^64 word.

**Fruitless-cycle overhead factor (a hard cap on cycle-handling steps).** The step meter
also carries the cycle-handling steps: the nominal charged work adds
OVER = 3/(12R) + 200/R^2 = 4.3488 x 10^-4 extra steps per step at the multiply-step
price (2-cycles at rate 1/(12R), <= 3 steps each; long cycles charged at rate R^-2,
<= 200 steps each). The branch makes the 2-cycle rate (1-s)^2/(12R) and the long-cycle
rate 0.0289 R^-2, so the charged event rates are 2.56 times and more than 34 times the
mixed rates; cycle handling is inside the same hard meter and the nominal charge is at most
1.000529 times the plain step work,
not an average. Section 3.2 bounds the all-multiply long-cycle
rate the charge is built on, and the branch can only lower it.

**Setup, output and preprocessing.** Setup and output are all <= 2^31 units
together: pairing g (Miller loop + final exponentiation, < 2^12 units, bounded by 2^14),
input check z^r = 1 (<= 510), table: 1024 multipliers x (2 x 255-bit exponentiations
<= 1,021 units + 11 Frobenius maps <= 1 unit each) <= 2^20.02, M = 2^20 random walk starts
g^alpha z^beta <= 2^30.0, output [x]P0 (<= 510 G1 operations, each < 1 unit) and
verification pairing (<= 2^14). **Preprocessing and table construction are at most 2^57
group-operation equivalents**: the multiplier table (<= 2^20.02, inside the setup term),
the walk starts (<= 2^30.0, inside it) and the initialization of all C_rep = 2^67 log
records with the sentinel record at 12 instructions each, 2^67 x 12 = 2^56.30 units.
Every one of those terms is charged inside W: the sentinel fill is listed in the
distinguished-report group below as 2^67 x 12 instructions, and the
setup/output and table terms are that group's 2^31 entry. No preprocessing is
amortised over future targets and no nonuniform advice is used (0 bytes).

**Restart starts**, charged at the worst case the meter allows: a slot's step counter can
reach 2^106.440033, and every restart needs 40/theta = 2^65.33 steps, so a slot has at
most 2^41.119 restarts; each costs one exponentiation (1,021 units, 2^24.29 instructions)
plus 100, so M x 2^41.119 x (1,021 x 20,000 + 100) = 2^71.12 units.

**Distinguished-report handling**, charged at its worst case (the deterministic cap, not
the mean), all with fixed instruction counts per record (section 3):

| Work per record (2^67 records at the cap) | Instructions |
| --- | ---: |
| append: the class key computed at the report (4 Fp multiplications 1,532, the addition and the two subtractions 133, the min/kappa rows 258 and the w^1 load 16, together 1,939, rounded up to 2,000), two modular multiplications by p^kappa (512-bit products, reduction, <= 400 each), 80-byte compare against the walk's last record, 128-bit record offset 80 (2^47 slot + r) (shift, add, multiplication by 80), 80-byte store, per-slot counter increment | 3,000 |
| full-record heapsort: <= 67 levels x 130 (10-word compare-and-swap of the 80-byte record), charged up from the computed 8,710 | 10,000 |
| adjacent dedup: one 80-byte comparison per record | 100 |
| **total per record** | **13,100** |

2^67 x 13,100 = 2^80.68 instructions = 2^66.39 units. The pair term is bounded by the
hard pair budget P = 2^41: P x 200 = 2^48.65 instructions = 2^34.36 units, plus one
canonical-element recomputation (1,021 units) per record a charged pair touches, at most
2 x P records = 2^42 x 1,021 < 2^52 units. Filling the log with the sentinel record before
the walks is 2^67 x 12 instructions = 2^56.3 units, which is the preprocessing/table
bound's dominant term above. The restart-starts term is 2^71.12 units,
which dominates this group, and the setup/output (2^31) and the log initialization
(2^56.3) are far below it: together **2^71.168 units**. The whole group is charged in
full whether or not the reports occur,
and no report-count probability enters it (the log's capacity is deterministic).

W <= 2^125.666542 + 2^71.168 units <= 2^125.666542 (every figure rounded up: the walk
bound is 2^125.6665413 to seven decimals, and the second term raises log2 W by
2^-54.49/ln 2 < 2^-53.9). In closed form, each factor rounded up,
W <= 1.01 x 0.9146875 x 2^125.780073 x 1.000529 <= 2^125.666543 (meter cap x average step
cost x step budget x cycle overhead; the log2 of the four factors, each rounded up, is
0.014356 - 0.128649 + 125.780073 + 0.000763 = 125.666543, above the exact ledger value).

Submitted bound: **time_log2 = 125.667** = 125667/1000 >= 125.666542 >= log2 W (part C of
experiment.py checks it as an exact integer inequality, W^1000 <= C^1000 x 2^125667, and
also checks that 125666 fails, so 125.667 is exactly the ceiling of log2 W at three
decimals).

Memory (a hard bound, every live array counted): the report log of C_rep = 2^67
records of 80 bytes is 2^67 x 80 < 2^73.3220 bytes, allocated in full and never
exceeded; the heapsort works in place with one 80-byte temporary; the multiplier table
is 12,288 x 576 bytes = 7.1 MB (2^22.76) plus its exponent tables (2 x 12,288 x 32
bytes < 2^20); the walk states are 2^20 x <= 1 KiB = 2^30 bytes (an Fp12 element 576
bytes, two 32-byte exponents, a T0 history ring of 48 bytes, an 80-byte last-record
buffer, the 128-bit work meter and the step/report/restart counters, < 900 bytes).
Total < 2^73.4 bytes. There is no hash table, no chained bucket and no eviction, so no
load factor arises. Memory is reported, not scored.

## 6. Heuristics and limitations

* **H1 (random-mapping behaviour of the mixed walk).** The class sequence is treated as a
  random mapping on N classes in which only a (1-s) fraction of the steps can merge two
  walks, with the Brent-Pollard correction. This is the standard heuristic behind every
  rho estimate (Teske 1998/2001; Brent-Pollard; van Oorschot-Wiener 1999 for parallel
  distinguished points; Duursma-Gaudry-Morain 1999 and Bos-Kleinjung-Lenstra 2010 for
  automorphism classes and fruitless cycles), applied here with the branch-adjusted
  collision law of section 4.1. It is supported by the toy runs of section 7 (part E
  measures the collision law directly, part B the success fraction and the collision
  time at three field sizes); it is not proved.
* **The decision bits of T0.** The branch probability is exactly 3/8 if the three bits
  of T0 are uniform, and the partition and distinguished fields are uniform if the other
  bits are. T0 = Re(c_0) is not a uniform Fp element, so these are heuristic. Support:
  part A2 measures the frequency of each of the 73 used bits on real BLS12-381 GT
  elements, part P runs the production rule on real GT elements and measures the
  squaring fraction, and part B's toy walks use the same rule and measure the branch and
  distinguished fractions against their nominal values. Because the three fields are
  disjoint bits of one value, a defect in one field cannot alias into another. The same
  uniformity assumption prices the 128-bit truncated report key (false matches, section
  3) and the probability that a healthy slot is stopped by the per-slot report cap.
* **Eventual coalescence, not a memoryless merge.** The ordinary step is a function of
  the class, but the two cycle escapes read stored history (the two-step T0 ring and the
  checkpoint slot). Section 3.1 proves the weaker invariant the matcher needs: the escape
  class and action are functions of the cycle, so two walkers reaching one cycle coincide
  from their first common class onward and coalesce within one cycle length; a class
  collision therefore still yields a usable pair, and coincident records are removed by
  the full-record dedup. Part P measures coalescence from different histories on a
  constructed all-multiply 2-cycle. This replaces the memoryless-merge claim the walk
  does not support.
* **Preprocessing and advice.** Preprocessing and table construction are at most 2^57
  group-operation equivalents, including the 2^56.30-unit initialization of all 2^67 log
  records with the sentinel record, the multiplier table (<= 2^20.02) and the walk starts
  and setup/output (<= 2^31). All of it is charged inside W (section 5), so no
  preprocessing is amortised; nonuniform advice is 0 bytes.
* **The work-meter allowance.** The abort's effect on success is charged as 0.01 with a
  model argument (section 4), not as a tail bound: the meter's own bound is hard, but no
  statement is made about the distribution of the work around its mean.
* **No production run.** The toy experiments are evidence only at their own scale: the
  toy's distinguished fraction is 2^-4 or 2^-3 rather than 2^-60, its field is near 2^20
  rather than 2^381, its class count is 2^21.4 to 2^22.4 rather than 2^251.3, and its T0
  is 20 bits rather than 381, so in part E (q = 67, T0 only 7 bits) the decision bits are
  derived from a hash of the canonical conjugate instead of T0, as section 7 states. The
  production cycle rate is the derivation of section 3.2, not a toy row.
* **The squaring price is an analytic instruction count** built from the measured
  blst_fp_mul = 383 and our own counted modular add and subtract, not a Callgrind count
  of a cyclotomic squaring routine. The count is 409 instructions below the charged
  11,700.
* This is a generic attack (about the square root of the class count). It ignores
  index-calculus methods in Fp12 (TNFS family); it is a valid upper bound, not a
  security estimate.

## 7. Experiments (`experiment.py`, standard library only)

Part A (BLS12-381 itself): for random GT elements (random Fp12 element raised to
(p^12-1)/r), checks y^r = 1; Frobenius formula = y^p; the 4-multiplication key formula
equals the keys of all 12 explicitly computed conjugates and they are distinct; the
sign-set minimum of section 2 equals the twelve-way argmin kappa and the class key K;
kappa selects the identical canonical element from each of the 12 conjugates; p^6 is
inversion; ord_r(p) = 12; 1 + p^4 + p^8 = 0 mod r; the factorisation of r-1 and
ord_r(2) = (r-1)/16 with the smallest t that can appear in a fruitless cycle,
ord_r(2)/gcd(ord_r(2),12) = 2^247.27; **T0 is equal on all 12 conjugates** of a GT
element; and the **Granger-Scott squaring equals y*y** on random GT elements, run
through Fp primitives so its operation counts (18 multiplications, 64 additions,
31 subtractions) are measured rather than asserted. Part A runs on every invocation with
at least one production element (the code floors `production_checks` at 1), so no run
prints uncomputed, initialized-true flags.

Part A2 (BLS12-381 itself, bit uniformity): it raises one random GT element to 64
random exponents mod r (a pool of random GT elements), forms up to 2,000 samples as
products of three pool members (each a random GT element; time-limited), reads T0 and
the class key K, and reports, for each of the low 73 bits of T0 (the 10 partition bits,
the 3 branch bits and the first 60 distinguished bits) and for the low 128 bits of K,
the frequency of a one bit and the largest deviation from 1/2 in standard errors, with
the mean deviation (0.80 expected for uniform bits) and the measured fractions of T0
whose branch field is below 3 and whose low 4 distinguished bits are zero, against 3/8
and 2^-4. A run at the initial inputs gave 2,000 samples: largest deviation 3.26
standard errors over 73 T0 bits and 3.13 over 128 K bits (about what uniform bits
give), mean 0.99 and 1.01, branch-field fraction 0.395 against 0.375 (1.85 binomial
standard errors) and low-4-zero fraction 0.0645 against 0.0625. A rerun differs in the
last digits because the sample count depends on the time limit. It supports the claim
that those bits behave as uniform; it cannot measure a 2^-60 event, which is why the
distinguished fraction is also measured in the toy.

Part P (BLS12-381 itself, the production step rule): two walkers start at two different
conjugates y0 and y0^{p^5} of one GT element and run the step rule in lockstep with a
16-entry table: the branch, partition and distinguished fields come from T0, keys are
computed only on multiply steps and at reports, the cycle checks compare T0, and the
escapes use the same mixed step function. The two walkers must produce the same T0
sequence (the ordinary step is a map on classes), the invariant y = g^alpha z^beta is
checked for both at the end, and the squaring fraction, the multiply steps and the number
of key computations are reported (key computations must equal multiply steps plus
reports). A run of 4,000 lockstep steps gave a squaring fraction of 0.3705, 2,746 key
computations against 2,515 multiply steps and 231 reports, identical T0 sequences for
the two walkers, and the exponent invariant holding for both.

Part P also runs the coalescence test of section 3.1 with different histories. A
2-cycle needs F(C1) = C2 and F(C2) = C1, so the test builds an all-multiply pair by
setting the two table cells the cycle uses: the rows are rebuilt from the same base
multipliers, so the exponent invariant and the charge are unchanged. Walker D starts at
C1 with the pair already in its two-step T0 ring (as if it had just come round the cycle)
and walker F starts at a different conjugate of C1 with an empty ring, so F must walk the
cycle once before detecting the return. Both must escape at C1_min^2, the smaller-T0
class of the pair squared, and their recorded class sequences (T0 and class key K) must
agree from that escape onward. A run gave: table cells holding M and M^-1 exactly, D
escaping at iteration 0 and F at iteration 2 (one cycle length apart), the same escape
class for both, 10 compared class steps after coalescence with identical sequences,
different exponent pairs at the first common class (so the meeting is a usable collision
rather than a coincident record), and the exponent invariant holding for both.

Part E (new: the branch-adjusted collision law, measured directly). In the toy field
q = 67 (N = 129,145 classes) distinct classes are sampled uniformly, one step is applied
to each, and the image classes are counted by pair type: squaring-squaring,
squaring-multiply, multiply-multiply with the same partition index and
multiply-multiply with different partition indices. At q = 67 T0 has 7 bits, too few for
the partition and branch fields, so the toy derives the decision bits from a hash of the
canonical conjugate of the class (a function of the class alone, as the production
partition, branch and distinguished fields are); this is the one place where the toy's
decision bits are not T0's. The pooled
statistic is N x (coincidences / pairs), whose expectation is 1 - s^2 - (1-s)^2/(12R);
the charged value 1 - s^2 - (1-s)^2/R is reported next to it, and squaring-squaring
coincidences must be exactly zero (squaring is injective on classes).

Part B (toy, same tower shape): q = 7 mod 12 near 2^20 with xi = 1+u neither a square
nor a cube (asserted), r' a prime factor of q^4 - q^2 + 1 with ord_r'(q) = 12
(asserted), N = (r'-1)/12. It runs the production walk of section 3: T0's three bit
fields drive the branch, the partition and the distinguished test (T0 has 20 bits at
this q, so the fields fit exactly as in production), the four-multiplication keys and
kappa are computed only on multiply steps, the squaring steps use the Granger-Scott
formula, the 2-cycle escape, checkpoint escape and stuck-walk rule work as in section 3
(both cycle checks before the distinguished test), and the append rule is the production
one. Toy-only differences: the 60 distinguished bits become dp_bits (the toy's T0 has 20
bits), the table has R <= 1024 entries, and ties of the minimum key (about 66/q per step
in the toy) are broken by the full conjugate.

Every trial runs its full fixed step budget 1.03 S0', where S0' is the branch-adjusted
budget of section 4 (1.079 to 1.083 times the plain S0 at s = 3/8) and no trial stops at its first
recovery: the walk phase is metered with the section 5 rule at toy scale (a slot is
stopped at 1.01 of its nominal charged work, and the toy's nominal work is the same
c x steps with c = 5/8 x 1.1125 + 3/8 x 0.585), so a trial's steps are a random number
bounded by its meter. The collected report log is then matched by the production
matcher, not by an online index: the log holds 80-byte records (the truncated key T and
the two canonical exponents), the matcher sorts it by the FULL record (T, then A, then
B) so identical records are adjacent, dedups by one adjacent full-record comparison,
enumerates every pair of the remaining distinct records of each equal-key run and
charges it against the toy's hard pair budget, recomputes the canonical element g^A z^B
for each record a charged pair touches, solves x from the first pair of equal elements
with different exponent pairs, and accepts x only when exponentiation verifies g^x = z.
A measurement-only class index timestamps each trial's first useful class collision so
that the collision-time distribution can be compared with the model; no success flag
reads it. Stress passes exercise the report cap (a cap of 4 reports per slot), the
stuck-walk restart (a threshold of 2 steps) with the slot's report counter carried
across, the work meter (a meter cap of 40 nominal steps) and the pair budget (a real log
and a two-record equal-key log at budgets 1 and 0).

Part C recomputes, without randomness, the zero-sum probabilities, the fruitless-cycle
bound of section 3.2 with and without the branch, the branch-adjusted N' and S0', the
deterministic report-log capacity M R_w, the KL bound on a healthy slot reaching the
cap, the memory bound, the per-record report-handling charge, the hard pair budget and
the charge it bounds, the work meter and its exact integer cap, the counter table of
section 5, and the W ledger of section 5 including the exact integer claim check. It also
checks the preprocessing bound of sections 5 and 6 from integers: the sentinel fill
(2^67 x 12 instructions = 2^56.30 units) together with the multiplier table and the
setup/output term is within 2^57 units, the fill is a term of W, the claim holds with it
in W, and leaving it out would change log2 W by less than 2^-53 (so W's ledger already
contains it and the claim 125.667 dominates log2 W exactly either way).

**Part E result** (2,600 distinct classes, q = 67, N = 129,145, R = 64, s = 3/8; the
decision bits come from a hash of the class because the toy's T0 has 7 bits):

| pair type | pairs | coincidences | N x rate | predicted |
| --- | ---: | ---: | ---: | ---: |
| squaring-squaring | 495,510 | 0 | 0 | 0 |
| squaring-multiply | 1,597,584 | 10 | 0.8084 | 1 |
| multiply-multiply, same partition | 20,073 | 0 | 0 | 11/12 |
| multiply-multiply, different partitions | 1,265,533 | 14 | 1.4287 | 1 |

The pooled statistic is N x coincidences/pairs = 0.9174 +- 0.1806, against 0.8528
predicted at the sample's own branch mix (branch fraction 0.3831 +- 0.0095), i.e.
z = +0.36; the nominal-law value is 0.8589 and the charged value 0.8533. Squaring-squaring
coincidences are exactly zero, as injectivity requires. The measurement cannot resolve the
difference between the exact and the charged collision factor (4.07 x 10^-4 of the factor
at the production R = 1024), and it does not need to: the charge is the larger N', so any
excess of the true rate over it makes the charge more conservative. Part E alone also
cannot separate 1 - s^2 (0.859) from 1 (no branch): its pooled statistic is less than 0.5
standard errors from 1. That separation comes from part B's collision time below.

**Part B result** (600 trials per field, 1,800 in all, seed 1, every trial run to the
fixed step budget 1.03 S0' and then the 20/theta detection window, with the slot work
meter at 1.01 of the nominal charged work; S0' is the branch-adjusted budget of section
4, so the plain budget would be up to 1.083 times smaller here):

| q, r', R, dp_bits | trials | P[T_c <= S0'] | P[T_c <= 1.03 S0'] | mean T_c/sqrt(N') | recovered and verified | 2-cycles/step (predicted (1-s)^2/(12R)) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1048963, 48674953, 64, 4 | 600 | 0.395 | 0.413 | 1.2541 +- 0.0273 | 0.468 | 6.42e-4 (5.09e-4) |
| 1049611, 64668529, 256, 4 | 600 | 0.363 | 0.383 | 1.2851 +- 0.0261 | 0.460 | 1.68e-4 (1.27e-4) |
| 1054003, 35456053, 1024, 3 | 600 | 0.378 | 0.392 | 1.2597 +- 0.0276 | 0.428 | 2.09e-4 (3.18e-5) |

Pooled over the 1,800 trials: P[T_c <= 1.03 S0'] = 0.396 +- 0.0115 against the model
0.40808 (1.1 standard errors below), and the pooled mean collision time is
T_c/sqrt(N') = 1.2664 against the random-mapping model sqrt(pi/2) = 1.2533, i.e. 1.05%
above the model and 0.9 standard errors from it (pooled standard error 0.0156). That last comparison is the load-bearing
one for the branch: without the branch adjustment the model mean would be sqrt(pi/2) =
1.2533 on sqrt(N) rather than on sqrt(N'), so a model that ignored the branch would sit
8.3% low, more than 6.5 standard errors below the measurement; the measurement is within
1.05% of the adjusted model. The recovered-and-verified
fraction (0.468, 0.460, 0.428; pooled 0.452 +- 0.0117) is *above* the in-budget collision
fraction because the walk phase then runs the detection window, during which further
collisions are found and the in-budget ones are reported; the production window is
E = 2^84.33 steps, i.e. a 2^-41.45 fraction of S, so in production this effect is
negligible, while the toy's 320-step window is 14% of its budget. Every logarithm the
matcher returned verified by exponentiation.

Other measured quantities, per row: squaring fraction 0.3747, 0.3758 and 0.3750 (nominal
0.375); distinguished fraction 0.06231, 0.06349 and 0.12901 (nominal 2^-4, 2^-4 and 2^-3);
long-cycle events 26, 25 and 20 at 1.70e-5, 1.44e-5 and 1.62e-5 per step, against the
charged R^-2 = 2.44e-4, 1.53e-5 and 9.54e-7 (the R = 64 and R = 256 rows are inside the
charge; the R = 1024 row is outside the derivation's regime, since R = 2^10 against
S0 = 2^10.9 steps there, while in production R/S0 = 2^-115.6); the report log stayed
within its cap in every trial (largest 209, 229 and 324 records against 328, 376 and 556)
and no slot was stopped by the report cap; the work meter stopped 499, 384 and 532 of the
2,400 slots, always within its bound (the stress pass's maximum slot meter 751,700 sits
below cap - 1 + 22,250 = 753,999, and all four slots stop at a 40-step cap); the meter's
relative deviation over a slot is 0.2914, 0.2841 and 0.2839 times sqrt(steps) against the
0.2792 predicted for independent steps; the pairs charged stayed at most 3 in a trial
against a toy budget of 2^36 times the worst-case pair count, with no abort, and the
pair-budget stress aborted a two-record log at budget 0 and charged one pair at budget 1.
Toy key ties (about 66/q per step) occurred 9, 16 and 3 times and are broken by the full
conjugate; the production tie term is 2^-249. No walk was abandoned by the stuck-walk
rule in the rows; the restart stress forced 68 restarts and kept the log at exactly its
16-record cap with the slot counters carried across.

The 2-cycle rate needs a word, because the toy's is 1.27 and 1.33 times the derivation
(1-s)^2/(12R) at R = 64 and R = 256 and 6.6 times it at R = 1024, where it also exceeds the
charged rate 1/(12R) (that row has R = 2^10 against S0 = 2^10.9 steps, far outside the
regime): at R = 64 and R = 256 the toy's R/S0' is 2^-5.1 and 2^-3.3, so a walk of
about 560 or 640 steps uses only 64 or 256 partitions and revisits indices by necessity,
and the 2-cycle and long-cycle rates are the two places where the toy is outside the
independence regime that production (R/S0 = 2^-115.6) satisfies. The cost charge does not
depend on that regime: it charges 3 steps per 2-cycle at rate 1/(12R), 2.56 times the
mixed rate (1-s)^2/(12R), and 200 steps per long cycle at R^-2, more than 34 times the mixed
0.0289 R^-2, both inside the meter's hard cap.

The tables above come from a longer run of the same code at seed 1 with a larger time
limit. The packaged experiment at its initial inputs (200 trials requested, 10 s budget,
seed 1) runs parts A2, P, E and B against a wall-clock limit, so the number of trials it
completes in each field, the part E sample size and every figure computed from them vary
with the machine's speed; its rows are smaller samples of the same measurement, reported
for reproducibility, not as the evidence.
