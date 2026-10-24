# Frobenius-class parallel Pollard rho in GT with a cyclotomic-squaring branch (success target 0.39)

Claim: time_log2 = **125.585** (total work W <= 2^125.585 dense-Fp12-multiplication
equivalents, C = 20,000), success probability >= 0.39, peak memory < 2^73.4 bytes
(a hard bound, section 5); preprocessing and table construction at most 2^57 units,
including the 2^56.30-unit log initialization, all included in W; no advice.
Heuristics are listed in section 6.

Design provenance. The reduced-torus class walk, the four-multiplication class-key
formula and the twisted multiplier table are the published Frobenius-class rho design
for this track. The report log, the offline matcher and the per-step instruction ledger
this package extends are the ones accepted on this track at 125.826; the squaring branch,
the trace-invariant decision bits and the branch-adjusted collision law are the mechanism
accepted at 125.74. This package re-derives all of it and differs in that the squaring step
is priced by measurement and the branch fraction re-optimised for that price: (1) every
per-step decision (partition, branch, distinguished test, cycle comparisons) is read
from one free class invariant T0, so a squaring step computes no class key at all and
the multiply step's bit reads are recounted; (2) the collision law, the step budget and
the fruitless-cycle rates are re-derived with the squaring branch in place (sections
3.2, 4); (3) because the charged work per step is now a random variable, the walk phase
carries a hard per-slot work meter whose abort rule and success allowance are stated
(section 5); (4) the report log, matcher, pair budget and memory bound are those of the
125.826 package, with the report margin restated for the larger step budget and the
append charge raised to carry the class key computed at a report.

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
Squaring-branch probability s = 35/64 = 0.546875, from a fixed 6-bit field of T0 (the
step squares iff that field's value is below 35). This is the W-minimising fraction the
package can state exactly: the continuous minimiser of log2 W is s = 0.5461, and on the
grids 2^-k the best value is 35/64 (s = 273/500, which a bit field cannot express, gives
log2 W lower by 1 x 10^-6, the same three-decimal claim). The reason there is an interior
optimum is that raising s shrinks the average step cost linearly but raises N' through
the factor 1 - s^2 - (1-s)^2/R, which flattens towards 1 - s^2.

Precompute M_h = g^{a_h} z^{b_h} (h < R, random a_h, b_h) and the table
T[h][k] = M_h^{p^k}, k = 0..11 (12,288 entries of 576 bytes, 7.1 MB), with exponent
tables a_h p^k, b_h p^k mod r.

**The four decision fields of T0** (disjoint bit ranges of one Fp element, read once
per step; T0 has 381 bits, so they fit with 308 bits to spare):

| bits of T0 | width | decision |
| --- | ---: | --- |
| 0 .. 9 | 10 | partition h = T0 mod R (1024 values) |
| 10 .. 15 | 6 | branch: square iff the field's value is below 35 (probability 35/64) |
| 16 .. 75 | 60 | distinguished iff all 60 bits are zero (probability 2^-60) |
| 76 .. 380 | 305 | unused |

The three fields are disjoint, so h, the branch and the distinguished flag are
functions of independent bits of one class invariant, and the marginal probability of
each is exact: 35/64 for the branch, 2^-60 for the distinguished flag. Part A2 measures
the frequency of each of the 76 used bits on real BLS12-381 GT elements (section 7).

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
if (squaring branch):                        [35/64 of steps]
    y = y^2                                  (cyclotomic squaring, section 5)
    alpha = 2 alpha ; beta = 2 beta (mod r)
else:                                        [29/64 of steps]
    compute class keys, kappa, K             (4 Fp mul + 407 other instructions)
    h = T0 mod R ; t = (12 - kappa) mod 12
    y = y * T[h][t]                          (ONE dense Fp12 multiplication)
    alpha += a_h p^t ; beta += b_h p^t (mod r)
```

So a multiply step costs one reference operation plus the key work, and a squaring
step costs the cyclotomic squaring (0.5049 reference operations, measured: section 5)
with no key work.
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
modulus); over at most 2 x 2^20 x 2^106.570101 < 2^127.571 comparisons (2 per step, for
fewer than 2^20 x SLOT_STEPS_MAX = 2^126.570101 steps in all: the hard step bound the
slot work meter allows, which already contains the walk budget S = 2^125.926969, the
detection window 2^84.33 and the restart waste 2^95.66) the total is below
2^127.571/p = 2^-253.13, i.e. below 2^-253
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
2^105.926969 + 2^64.321928 walk steps is 2^45.926969 reports, so the cap R_w = 2^47 is
a = 2^1.0738 = 2.1038 times the mean and the Bernoulli-KL (Chernoff) bound gives
log2 P(slot reaches the cap) <= -mu D(a||1)/ln 2 = -4.45 x 10^13. That is a heuristic
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
  (part A computes both numbers). A cycle has at most 2^126.6 steps in all, so
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
sum_{L=3..64} P_L (1-s)^L = 0.0101 R^-2 for R = 1024 (against 0.1391 R^-2 without the
branch), and sum_{L=65..2R} P_L (1-s)^L <= 2^-311.8. For L > 2R every term with
d < R shrinks by d/R per extra step, so the tail is negligible as in the base
derivation. The cost ledger nonetheless charges the **unbranched** R^-2 per step with
200 steps per event (section 5), which is 99 times the mixed rate and more than
7 times the unbranched one: the squaring branch can
only reduce fruitless-cycle work, and the charge does not depend on that reduction.

**Class cycles, within-walk returns, and what a lap can cost.** If a cycle of length L
contains a squaring step, then after one lap the class returns but the exponent pair is
multiplied by 2^t p^m with 2^t != p^e, so the pair is different: two records with the
same key T and different (A, B), which is exactly a solved instance. A lap can therefore
enlarge an equal-key run of the report log, and the only place that hurts is the hard
pair budget P = 2^41 of section 3. The accounting below distinguishes every kind of
revisit explicitly, because within-walk returns are common and must not be mistaken for
a rare failure.

(i) **Algebraic fruitless 2-cycles are handled events, not failures.** They occur at rate
(1-s)^2/(12R) = (29/64)^2/12288 = 1.67 x 10^-5 per step (part B measures 1.92 x 10^-5 at
R = 1024, part C computes the rate), i.e. of order 2^110.1 occurrences over the
S = 2^125.93 steps of the run: they are not an event of probability 2^-20, they are the
normal behaviour of the walk. Each is DETECTED by the T0 comparison against the class two
steps back and ESCAPED by squaring the stored representative of the smaller-T0 class of
the pair. Its step cost (at most 3 steps) is already inside the overhead factor
ceil(OVER x S_slot) charged at c_mul in section 5, and the escape is a function of the
class, not of the walker: the class that is squared is the one with the smaller T0, so the
destination is a function of the two-cycle's T0 multiset and two walkers on one cycle
escape at the same class (section 3.1). A fruitless 2-cycle is therefore a charged,
handled event.

(ii) **Checkpoint cycles of length 3..64 are handled the same way.** The 64-step
checkpoint ring detects them, the walk traverses the cycle once to its smallest-T0 class
and squares it, and their step cost is inside the same overhead factor. Their per-step
entry rate is the section-3.2 sum 0.0101 R^-2 = 2^-189.1 for R = 1024, which the ledger
charges at the unbranched R^-2 with 200 steps per event.

(iii) **Every other revisit is a tail collision, and with different exponent pairs it is
useful.** A revisit that is not one of (i) or (ii) reaches a class already visited with a
different exponent pair, i.e. two records of one key T with different (A, B): that is the
collision the algorithm is looking for, x = (A1 - A2)/(B2 - B1) mod r, not a failure. The
one case that is neither escaped nor useful is a revisit with an equal exponent pair, and
that means the walk is on a class cycle longer than the checkpoint window: the two
exponent conditions 2^t = p^e mod r and sum_i p^{t_i} = 0 mod r of section 3.2 force such
a cycle to be all-multiply with length above 64, and the section-3.2 rate bounds the
chance that a walk enters one at all by S x 2^-311.8 < 2^-185, far inside the
long-fruitless-cycle loss term of section 4 (below 2^-100, which is the figure carried).
Even then the record cap and the pair budget are hard caps that stop the slot.

(iv) **Consequence for the failure ledger: there is no birthday term for within-walk
returns.** Cases (i) and (ii) are handled inside the charged cycle overhead, case (iii)
with different exponent pairs is the useful collision, and case (iii) with equal exponent
pairs is a long class cycle, already below 2^-100. Nothing is left to charge, so the
failure list of section 4 carries no term of the form M S_slot^2/(2N) for returns, and the
earlier claim that the record-count argument needs the absence of a within-walk return is
withdrawn: the report log's capacity is deterministic and rests on no such assumption.
Each slot owns R_w = 2^47 record positions and stops at its cap, so the log holds at most
C_rep = M R_w = 2^67 records whatever the walk does, with no overflow, no probability and
no return assumption (section 3); and a lap that adds records of one class to one run is
charged as pairs of that run like any other, against the hard pair budget P = 2^41 whose
own abort probability is the 2^-36 Markov charge of section 3. The accepted 125.826 package carried these exposures in its pair-budget term.

## 4. Success probability

Heuristic (H1): the class sequence of the mixed walk behaves like a random mapping on
N classes in which only a (1-s) fraction of the steps can merge two walks, with the
Brent-Pollard correction for R-adding walks. Two steps merge two distinct classes when
both are multiply steps with different multiplier indices (probability
1 - s^2 - (1-s)^2/(12R) per pair of steps, section 4.1), so the effective class count is

  N' = N / (1 - s^2 - (1-s)^2/R) = 2^251.785202

where the charge uses 1/R instead of 1/(12R) for the same-index case: same-partition
pairs are charged as never merging even when their twist index differs, which only
enlarges N' (conservatively by 2.62 x 10^-4 relative, 3.78 x 10^-4 in log2, since the two
factors differ by (1-s)^2 (1/R - 1/(12R)) = 1.84 x 10^-4 against the factor 0.7007). Then
the probability that no useful collision occurs among the first S
distinct-step points of all walks is exp(-S^2/(2N')). Choose

  S0' = sqrt(2 N' ln(1/0.61)) = 2^125.884326,   S = 1.03 S0' = 2^125.926969,

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
more false-key pairs than the hard pair budget (< 2^-36 by Markov, section 3); a slot
whose charged work exceeds 1.01 of its nominal work, i.e. the work meter aborting before
the slot takes its S_slot productive steps (< 2^-100 by Hoeffding, see the concentration
bound below); a cycle-handling step count above its budget (probability 0: cycle handling
is inside the slot's hard work meter); a walk slot stopped by the per-slot report cap,
i.e. a healthy slot reaching R_w = 2^47 reports (the KL figure of section 3, < 2^-300);
two distinct classes with equal T0 causing a spurious escape (< 2^-253). There is no
birthday term for within-walk returns: section 3.2 shows that a return is either escaped
inside the charged cycle overhead or is the useful collision. Their sum is below
2.39 x 10^-6 < 2^-18 < 0.4080 - 0.39 - 0.01 = 0.0080.

**The work-meter concentration bound.** Total charged work is a random variable (a
squaring step costs c_sq = 0.50490 units, a multiply step c_mul = 1.1125), so the walk
phase is metered and a slot is stopped at W_cap = 1.01 x its nominal charged work
(section 5). The nominal work already covers the slot's S_slot + E_slot steps at the mean
step cost c = (29/64) c_mul + (35/64) c_sq = 0.78021875 units plus its cycle overhead, so
the only way the cap costs the run a productive step is a slot whose step costs happen to
sum above 1.01 of that nominal. Costs are bounded (c_sq <= cost <= c_mul) and independent
of the past, so Hoeffding bounds that event exactly:

  P(total work of a slot > 1.01 c n) <= exp(-2 (0.01 c)^2 n / (c_mul - c_sq)^2),

with n = S_slot + E_slot = 2^105.93 steps and c_mul - c_sq = 0.6076. The exponent is
2 (0.01 x 0.78021875)^2 / 0.6076^2 = 3.298 x 10^-4 nats per step, i.e. 2^-94.89 for a
slot and 2^-114.89 for the run's S = 2^125.93 steps: the bound is far below 2^-100 (part
C computes both exponents and asserts the first exceeds 100). Over the M = 2^20 slots the
union bound stays below 2^-100, which is the figure in the loss list. The hard cap
therefore retains all S productive steps except with that probability; there is no
"We retain S/1.01 steps" equivalence and none is used. A slot that does abort has its
loss inside the same 2^-100 term.

Success >= 0.40808 - 0.01 - 2^-18 > 0.40808 - 0.01 - 0.00001 = 0.39807, so
success >= 0.3980 >= 0.39 (every figure rounded down), with 0.0080 of margin. The 0.01
allowance is kept as slack for the heuristic parts of H1 and for the steps a stopped slot
never takes; it is not used to cover the concentration bound, which is a probability
below 2^-100.

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

**Squaring step** (measured: the shipped cyclotomic squaring routine):

| Work | Instructions |
| --- | ---: |
| cyclotomic squaring of the GT element, measured (see below) | 9,798 |
| alpha, beta doubling mod r (2 x 4-limb doubling, conditional subtract) | 52 |
| common rows of the table above (T0 read, bit fields, cycle comparisons, meter, counters, control) | 101 |
| **total counted (budgeted 10,098)** | **9,951** |

The squaring is **measured, not modelled**. The routine is blst's
`blst_fp12_cyclotomic_sqr` at the pinned blst revision f7a88067: it costs **9,798
executed instructions per call**, counted with Callgrind on the track's own COST_MODEL
toolchain (Debian cc 14.2.0, valgrind-3.24.0, `-O3 -fno-tree-vectorize
-mgeneral-regs-only`), identical at five slice offsets over five different GT elements
(one-call inclusive count 9,798; the marginal over 64 contiguous calls is 9,803, the 5
being the call plumbing outside the counted arc). The count cannot depend on the input:
objdump shows the routine's call tree is 18 functions and 1,017 static instructions with
no conditional, backward or indirect branch, so the executed count is the same for every
element. The routine is the legitimate squaring here because every walk element is in GT
(and every report key comes from a GT element), which lies in the cyclotomic subgroup,
where `cyclotomic_sqr_fp12` is defined and equals `blst_fp12_sqr` (checked on GT elements
and after 997 further squarings). The same harness on the same toolchain reproduces
blst_fp_mul = 383.0 and blst_fp12_mul = 22,996.0 instructions exactly (0.00% on both),
which are the two prices the rest of this ledger is built on. Counting the Granger-Scott
formula line by line through Fp primitives instead (18 Fp multiplies, 64 adds, 31
subtractions, priced as standalone blst_fp_mul calls, as part A counts them) gives 11,190
and so over-charges the shipped routine by 1,340: the shipped code squares through a
double-width `sqr_382x`/`redc_mont_384` path rather than 18 separate full-precision
multiplies. The measurement is what is charged.

The charge is 10,098 instructions = 9,798 measured + 300 for the alpha/beta doubling, the
shared rows and slack (the same convention the earlier 11,700 = 11,400 + 300 used); the
counted rows above are 9,951, so the budget is 147 instructions above the count. A
squaring step computes no class key: the four Fp multiplications of the key appear only on
multiply steps and at reports. Squaring-step cost c_sq = 10,098/20,000 = 0.5049 units.
Average step cost, with s = 35/64, is
c = (29/64) c_mul + (35/64) c_sq = 0.78021875 units (15,604.375 instructions).

**Work meter (a hard cap on a random per-step cost).** Because a step costs 1.1125 or
0.5049 units, the walk phase's work is a random variable, so it is metered and capped
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
  I_MUL x ceil(OVER x S_slot) = 1204233582164114075842234780252897280 = 2^119.857527
  instructions, and W_cap = 1216275917985755216600657128055426253 = 2^119.871882.
  The meter's largest value is one over-cap update past the test, W_cap - 1 + I_MUL =
  1216275917985755216600657128055448502 = 2^119.871882 (log2 119.871882), which fits two
  64-bit words with the cap test at 2^119.871882 < 2^127. A slot is stopped when its meter
  reaches W_cap, so
  W_walk <= M x (W_cap - 1 + I_MUL) = 2^139.871882 instructions = 2^125.584170 units is a
  **hard bound**, not an expectation; it is what the claim is computed from.
* Abort bound. A slot's nominal work is c x C x (S_slot + E_slot) plus the cycle
  overhead, and its actual work is a sum of S_slot + E_slot independent step costs in
  [c_sq, c_mul], so Hoeffding gives P(work > 1.01 nominal) <= exp(-2 (0.01 c)^2 n /
  (c_mul - c_sq)^2) = 2^-94.89 for the slot and below 2^-100 over the M slots (section 4
  carries it in the loss list). The hard cap therefore stops a slot only in that event,
  and it is the 0.01 success allowance of section 4 that covers the heuristic parts of H1,
  not an assumed distribution of the work.

**Counters and indices.** Every counter and index, its largest value at the end of the
attack (including one over-cap update per processor), and its width. Part C recomputes
this table and asserts every row.

| counter or index | largest value (log2) | words | capacity test |
| --- | ---: | ---: | --- |
| slot work meter (instructions) | 119.871882 | 2 x 64 | >= W_cap with W_cap = 2^119.871882 < 2^127 |
| slot step counter | 106.570101 | 2 x 64 | none (informational) |
| steps since last report | 65.321929 | 2 x 64 | > 40/theta (2^65.321929 rounded up), a 128-bit compare |
| slot report counter r | 47 | 1 x 64 | >= 2^47 with 2^47 < 2^63 |
| log offset 80 x (2^47 slot + r) | 73.321929 | 2 x 64 | < 2^128 |
| slot restart counter | 41.248172 | 1 x 64 | < 2^63 |
| matcher pair counter | 41 | 1 x 64 | > 2^41 with 2^41 + 1 < 2^63 |
| heapsort index | 67 | 2 x 64 | < 2^128 |

No counter's capacity test sits at its word's own capacity: the two cap tests that can
fire are at 2^119.871882 of a 2^128 word and at 2^47 of a 2^64 word.

**Fruitless-cycle overhead factor (a hard cap on cycle-handling steps).** The step meter
also carries the cycle-handling steps: the nominal charged work adds
OVER = 3/(12R) + 200/R^2 = 4.3488 x 10^-4 extra steps per step at the multiply-step
price (2-cycles at rate 1/(12R), <= 3 steps each; long cycles charged at rate R^-2,
<= 200 steps each). The branch makes the 2-cycle rate (1-s)^2/(12R) and the long-cycle
rate 0.0101 R^-2, so the charged event rates are 4.87 times and 99 times the
mixed rates; cycle handling is inside the same hard meter and the nominal charge is at most
1.000621 times the plain step work,
not an average. Section 3.2 bounds the all-multiply long-cycle
rate the charge is built on, and the branch can only lower it.

**Setup, output and preprocessing.** Setup and output are all <= 2^31 units together:
pairing g (< 2^12, bounded by 2^14), input check z^r = 1 (<= 510), table: 1024
multipliers x (2 exponentiation <= 1,021 + 11 Frobenius maps) <= 2^20.02, M = 2^20 walk
starts <= 2^30.0, output [x]P0 (<= 510 G1 operations) and verification pairing (<= 2^14).
**Preprocessing and table construction are at most 2^57 group-operation equivalents**: the
table (<= 2^20.02, inside the setup term), the walk starts (inside it) and the
initialization of all C_rep = 2^67 log records with the sentinel record, 2^67 x 12 =
2^56.30 units. Every term is charged inside W (the fill in the report group below, the
rest in its 2^31 entry); nothing is amortised and no nonuniform advice is used.

**Restart starts**, charged at the worst case the meter allows: a slot's step counter can
reach 2^106.570101, and every restart needs 40/theta = 2^65.33 steps, so a slot has at
most 2^41.248172 restarts; each costs one exponentiation (1,021 units, 2^24.29
instructions) plus 100, so M x 2^41.248172 x (1,021 x 20,000 + 100) = 2^71.244 units.

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
canonical recomputation (1,021 units) per record a charged pair touches, at most 2 x P
records < 2^52 units. The sentinel fill is 2^67 x 12 = 2^56.3 units. The restart-starts
term (2^71.244) dominates the group, and setup/output (2^31) and the fill are far below
it: together **2^71.293 units**, charged in full whether or not the reports occur, with no
report-count probability (the capacity is deterministic).

W <= 2^125.584170 + 2^71.293 units = 2^125.584170 (every figure rounded up: the walk
bound is 2^125.5841696 to seven decimals, and the second term raises log2 W by
2^-54.29/ln 2 < 2^-53.8, below the last printed digit). Closed form, each factor rounded
up: W <= 1.01 x 0.78021875 x 2^125.926969 x 1.000621 <= 2^125.584170 (meter cap x average
step cost x step budget x cycle overhead; the log2 of the four is 0.014356 - 0.358049 +
125.926969 + 0.000895 = 125.584171, above the exact 125.5841696).

Submitted bound: **time_log2 = 125.585** = 125585/1000 >= 125.584170 >= log2 W (part C of
experiment.py checks it as an exact integer inequality, W^1000 <= C^1000 x 2^125585, and
also checks that 125584 fails, so 125.585 is exactly the ceiling of log2 W at three
decimals).

Memory (a hard bound, every live array counted): the report log of C_rep = 2^67 records
of 80 bytes is 2^67 x 80 < 2^73.3220 bytes, allocated in full and never exceeded; the
heapsort works in place with one 80-byte temporary; the multiplier table is 12,288 x 576
bytes = 2^22.76 plus its exponent tables (< 2^20); the walk states are 2^20 x < 1 KiB =
2^30 bytes. Total < 2^73.4 bytes. There is no hash table, chained bucket or eviction, so
no load factor arises. Memory is reported, not scored.

## 6. Heuristics and limitations

* **H1 (random-mapping behaviour of the mixed walk).** The class sequence is treated as a
  random mapping on N classes in which only a (1-s) fraction of the steps can merge two
  walks, with the Brent-Pollard correction: the standard heuristic behind rho estimates
  (Teske 1998/2001; Brent-Pollard; van Oorschot-Wiener 1999; Duursma-Gaudry-Morain 1999;
  Bos-Kleinjung-Lenstra 2010), applied with the branch-adjusted collision law of section
  4.1. It is supported by section 7 (part E measures the collision law, part B the success
  fraction and collision time); it is not proved.
* **The decision bits of T0.** The branch probability is exactly 35/64 if the six branch
  bits of T0 are uniform, and the partition and distinguished fields are uniform if the
  other bits are. T0 = Re(c_0) is not a uniform Fp element, so these are heuristic.
  Support: part A2 measures the frequency of each of the 76 used bits on real BLS12-381 GT
  elements, part P runs the production rule on real GT elements, and part B's toy walks
  measure the branch and distinguished fractions against their nominal values. Because the
  fields are disjoint bits of one value, a defect in one cannot alias into another. The
  same assumption prices the 128-bit report key and the report-cap margin.
* **Eventual coalescence, not a memoryless merge.** The ordinary step is a function of
  the class, but the two cycle escapes read stored history. Section 3.1 proves the weaker
  invariant the matcher needs: the escape class and action are functions of the cycle, so
  two walkers reaching one cycle coalesce within one cycle length, and a class collision
  still yields a usable pair. Part P measures coalescence from different histories on a
  constructed 2-cycle.
* **Preprocessing and advice.** Preprocessing and table construction are at most 2^57
  group-operation equivalents, including the 2^56.30-unit initialization of all 2^67 log
  records with the sentinel record, the multiplier table (<= 2^20.02) and the walk starts
  and setup/output (<= 2^31). All of it is charged inside W (section 5), so no
  preprocessing is amortised; nonuniform advice is 0 bytes.
* **The work meter and its concentration bound.** The meter's own bound on total work is
  hard (section 5); the event that a slot is stopped before its productive step budget is
  bounded by Hoeffding, because step costs are bounded and history-independent (section
  4): below 2^-100 over the run, so the cap retains all S productive steps except with
  that probability. Nothing rests on the work's distribution around its mean, and no
  "S/1.01 steps" equivalence is used. The separate 0.01 success allowance is slack for the
  heuristic parts of H1.
* **No production run.** The toys are evidence only at their own scale: distinguished
  fraction 2^-4 or 2^-3 rather than 2^-60, field near 2^20 rather than 2^381, class count
  2^21.4 to 2^22.4 rather than 2^251.3, T0 20 bits rather than 381 (and 7 in part E, where
  the bits come from a hash). The production cycle rate is the derivation of section 3.2,
  not a toy row.
* **The squaring price is a measured instruction count of one shipped routine.** It is
  the Callgrind count of blst f7a88067 `blst_fp12_cyclotomic_sqr` (9,798 per call, plus
  300), not a model of the formula. Support: CALIBRATION on the same harness, which
  reproduces blst_fp_mul = 383.0 and blst_fp12_mul = 22,996.0 to the instruction, and a
  straight-line body (no conditional, backward or indirect branch in its 18-function call
  tree), so the count is input-independent by construction. The residual risk is different
  code generation on the grader's compiler, which is why the charge keeps 147 instructions
  of slack and why COST_MODEL's own toolchain was used. The line-by-line formula would give
  11,190, which the shipped double-width `sqr_382x` path beats by 1,340.
* This is a generic attack (about the square root of the class count). It ignores
  index-calculus methods in Fp12 (TNFS family); it is a valid upper bound, not a
  security estimate.

## 7. Experiments (`experiment.py`, standard library only)

Part A (BLS12-381 itself): on random GT elements, checks y^r = 1; Frobenius formula =
y^p; the 4-multiplication key formula equals the keys of all 12 conjugates, distinct; the
sign-set minimum equals the twelve-way argmin kappa and the class key; kappa selects the
same canonical element from all 12 conjugates; p^6 is inversion; ord_r(p) = 12;
1 + p^4 + p^8 = 0 mod r; the factorisation of r-1 and ord_r(2) = (r-1)/16 with the smallest
t a fruitless cycle could contain, 2^247.27; T0 equal on all 12 conjugates; and the
cyclotomic squaring equal to y*y, run through Fp primitives so its operation counts (18
multiplications, 64 additions, 31 subtractions) are measured, not asserted.

Part A2 (bit uniformity): up to 2,000 samples, each a product of three random GT
elements (time-limited), read T0 and the class key K and report, for each of the low 76
bits of T0 (10 partition, 6 branch, 60 distinguished) and the low 128 bits of K, the
frequency of a one bit and the largest deviation from 1/2 in standard errors, with the
mean deviation and the measured fractions whose branch field is below 35 and whose low 4
distinguished bits are zero, against 35/64 and 2^-4. A run at the initial inputs gave
2,000 samples: largest deviation 3.26 standard errors over 76 T0 bits and 3.13 over 128 K
bits, mean 1.03 and 1.01, branch-field fraction 0.5565 against 0.546875 (0.86 binomial
standard errors) and low-4-zero fraction 0.0615 against 0.0625. A rerun differs in the
last digits because the sample count depends on the time limit. It supports uniformity of
those bits; it cannot measure a 2^-60 event, which the toy measures instead.

Part P (the production step rule on real GT elements): two walkers start at two conjugates
y0 and y0^{p^5} of one GT element and run the step rule in lockstep with a 16-entry table,
all fields from T0, keys only on multiply steps and at reports, cycle checks on T0. They
must produce the same T0 sequence, keep y = g^alpha z^beta, and report key computations
equal to multiply steps plus reports. A run of 4,000 steps gave squaring fraction 0.5524
(2,209 + 1,790 steps), 2,042 key computations against 1,790 multiply steps and 252
reports, identical T0 sequences and the invariant holding for both.

Part P also runs the coalescence test of section 3.1 on a constructed all-multiply
2-cycle: two walkers start at class C1 with different histories, one with the cycle already
in its two-step T0 ring and one freshly started. A run gave escapes at iterations 0 and 2
(one cycle length apart), the same escape class (the smaller-T0 class squared), 10
identical class steps after coalescence, and different exponent pairs at the first common
class (a usable collision, not a coincident record).

Part E (the branch-adjusted collision law, measured directly). In the toy field q = 67
(N = 129,145 classes) distinct classes are sampled uniformly, one step is applied to each,
and image coincidences are counted by pair type: squaring-squaring, squaring-multiply,
multiply-multiply with the same partition index and with different indices. At q = 67 T0
has 7 bits, too few for the fields, so the decision bits come from a hash of the canonical
conjugate (a function of the class alone, as the production fields are); this is the one
place where the toy's bits are not T0's. The pooled statistic is N x coincidences/pairs,
whose expectation is 1 - s^2 - (1-s)^2/(12R); the charged value 1 - s^2 - (1-s)^2/R is
reported next to it, and squaring-squaring coincidences must be exactly zero.

Part B (toy, same tower shape): q = 7 mod 12 near 2^20, r' a prime factor of q^4 - q^2 + 1
with ord_r'(q) = 12, N = (r'-1)/12. It runs the production walk of section 3 at toy scale
(T0 has 20 bits, so the fields fit as in production), every trial to its full fixed budget
1.03 S0' with the section-5 meter and no stop at the first recovery, and then matches the
collected log with the production matcher (full-record sort, adjacent dedup, charged pairs
against a hard budget, canonical recomputation, verification of g^x = z). A
measurement-only class index timestamps each trial's first useful collision; no success
flag reads it. Stress passes exercise the report cap, the stuck-walk restart with the
carried slot counter, the work meter and the pair budget.

Part C recomputes, without randomness, the zero-sum probabilities, the fruitless-cycle
bounds, N' and S0', the deterministic log capacity, the KL report-cap bound, the memory
bound, the per-record report charge, the pair budget, the work meter's exact integer cap
and its Hoeffding exponent, the counter table, and the W ledger with the exact integer
claim check and the preprocessing bound (all within 2^57 units, the fill inside W and worth
less than 2^-53 of it).

**Part E result** (2,500 distinct classes, q = 67, N = 129,145, R = 64, s = 35/64):

| pair type | pairs | coincidences | N x rate | predicted |
| --- | ---: | ---: | ---: | ---: |
| squaring-squaring | 905,185 | 0 | 0 | 0 |
| squaring-multiply | 1,553,284 | 10 | 0.8314 | 1 |
| multiply-multiply, same partition | 10,213 | 0 | 0 | 11/12 |
| multiply-multiply, different partitions | 655,068 | 6 | 1.1829 | 1 |

The pooled statistic is 0.6615 +- 0.1713 against 0.7100 predicted at the sample's own
branch mix (branch fraction 0.5384 +- 0.0100), i.e. z = -0.28; the nominal-law value is
0.7007 and the charged value 0.6977. Squaring-squaring coincidences are exactly zero, as
injectivity requires. The measurement cannot resolve the difference between the exact and
the charged collision factor (2.62 x 10^-4 of the factor at R = 1024), and it does not need
to: the charge is the larger N', so any excess of the true rate makes it more conservative.
Part E alone also cannot separate 1 - s^2 (0.7009, the s^2 term is 0.299) from 1 (no
branch): its statistic is 1.98 standard errors from 1, and that separation comes from part
B's collision time.

**Part B result** (600 trials per field, 1,800 in all, seed 1, a 600 s limit; every trial
run to 1.03 S0' and then the 20/theta detection window with the meter at 1.01 of the
nominal charged work; S0' is the branch-adjusted budget, so the plain budget is 1.195 times
smaller):

| q, r', R, dp_bits | trials | P[T_c <= S0'] | P[T_c <= 1.03 S0'] | mean T_c/sqrt(N') | recovered and verified | 2-cycles/step (predicted) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1048963, 48674953, 64, 4 | 600 | 0.415 +- 0.020 | 0.425 | 1.1840 +- 0.0262 | 0.475 | 3.16e-4 (2.67e-4) |
| 1049611, 64668529, 256, 4 | 600 | 0.395 +- 0.020 | 0.407 | 1.2257 +- 0.0264 | 0.458 | 9.54e-5 (6.68e-5) |
| 1054003, 35456053, 1024, 3 | 600 | 0.355 +- 0.020 | 0.372 | 1.2766 +- 0.0267 | 0.397 | 3.46e-5 (1.67e-5) |

Pooled over the 1,800 trials: P[T_c <= 1.03 S0'] = 722/1,800 = 0.401 +- 0.0116 against the
model 0.40808 (0.6 standard errors below), and the pooled mean collision time is
T_c/sqrt(N') = 1.2287 +- 0.0153 against sqrt(pi/2) = 1.2533, i.e. 2.0% below it and 1.6
standard errors below. A model that ignored the branch adjustment would place the same
sqrt(pi/2) on sqrt(N) instead of sqrt(N'), i.e. 1.2533/sqrt(1.428) = 1.0488 in these units,
and the measurement is 17% above that, 11.8 standard errors away: it sits on the
branch-adjusted scale and excludes the unadjusted one, while running slightly early (2%),
which is the safe direction for the success claim. The recovered-and-verified fraction
(0.475, 0.458, 0.397; pooled 0.443 +- 0.0117) is above the in-budget collision fraction
because the walk phase then runs the detection window, during which further collisions are
found; the production window is E = 2^84.33 steps, a 2^-41.65 fraction of S, so there this
effect is negligible, while the toy's 320- and 160-step windows are 13% and 8% of their
budgets. Every logarithm the matcher returned verified by exponentiation.

Other per-row quantities: squaring fraction 0.5474, 0.5471, 0.5487 (nominal 0.546875);
distinguished fraction 0.06276, 0.06334, 0.1291 (nominal 2^-4, 2^-4, 2^-3); long-cycle
events 45, 42, 34; the log stayed within its cap (largest 214, 264, 352 records against
328, 376, 556) with no slot stopped by the report cap; the meter stopped 681, 643, 619 of
2,400 slots, always within its bound (stress maximum 634,808 < 646,424); its relative
deviation over a slot is 0.3955, 0.4094, 0.3864 sqrt(steps) against 0.3877 predicted; at
most 2, 3, 3 pairs were charged with no abort (the stress aborted a two-record log at
budget 0 and charged one pair at budget 1); toy key ties occurred 12, 21, 14 times (the
production tie term is 2^-249); no walk was abandoned, and the restart stress forced 130
restarts at the 16-record cap with the slot counters carried across.

The toy's 2-cycle rate is 1.18, 1.43 and 2.07 times the derivation (1-s)^2/(12R) at
R = 64, 256 and 1024, and its long-cycle rate at R = 256 and R = 1024 (2.21e-5, 2.50e-5 per
step against the charged 1.53e-5 and 9.54e-7) exceeds the charge: at R = 64 and R = 256
the toy's R/S0' is 2^-5.1 and 2^-3.3, so a walk of about 560 or 640 steps uses only 64 or
256 partitions and revisits indices by necessity, which is outside the independence regime
production satisfies (R/S0 = 2^-115.6). The R = 64 long-cycle row (2.69e-5 against 2.44e-4)
is inside the charge. The cost charge does not depend on that regime: 3 steps per 2-cycle
at rate 1/(12R), 4.87 times the mixed rate, and 200 steps per long cycle at R^-2, 99 times
the mixed 0.0101 R^-2, both inside the meter's hard cap.

These tables come from a longer run of the same code at seed 1 with a larger time limit. The
packaged experiment at its initial inputs (200 trials requested, 10 s budget, seed 1) runs
parts A2, P, E and B against a wall-clock limit, so the trials it completes in each field,
the part E sample size and every figure computed from them vary with the machine's speed;
its rows are smaller samples of the same measurement, reported for reproducibility, not as
the evidence.
