# Frobenius-class parallel Pollard rho in GT with a cyclotomic-squaring branch (success target 0.39)

Claim: time_log2 = **125.571** (total work W <= 2^125.571 dense-Fp12-multiplication
equivalents, C = 20,000), success probability >= 0.39, peak memory < 2^73.4 bytes
(a hard bound, section 5); preprocessing and table construction at most 2^57 units,
including the 2^56.30-unit log initialization, all included in W; no advice.
Heuristics are listed in section 6.

Design provenance. The reduced-torus class walk, the four-multiplication class-key formula
and the twisted multiplier table are the published Frobenius-class rho design for this
track; the report log, the offline matcher and the per-step instruction ledger are the ones
accepted here at 125.826, and the squaring branch, the trace-invariant decision bits and the
branch-adjusted collision law the mechanism accepted at 125.74. This package re-derives all
of it; it differs in that the squaring step is priced by measurement and the branch fraction
re-optimised for that price (1), and it carries the hard meter (2) and the attempt-charged
report cap (3) that the accounting needs. (1) every per-step decision is read from one free
class invariant T0, so a squaring step computes no class key; (2) step costs are random, so
the walk carries a per-slot work meter that also charges each distinguished-point attempt,
with its abort rule and success allowance (section 5); (3) the report cap counts attempts,
and the log, matcher, pair budget and memory bound are those of 125.826 with the report
margin restated.

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
* **The decisions are read in the Montgomery domain.** blst stores an Fp element as six
  64-bit limbs in Montgomery form, a R mod p with R = 2^384, and its kernels return their
  outputs in that form. Every per-step decision in section 3 is therefore read from
  T0M = T0 R mod p, the word the measured kernels produce: the partition, branch and
  distinguished fields, both cycle comparisons, the "smaller-T0" and "smallest-T0"
  orderings that pick an escape class, the class-key ordering that defines kappa, and the
  report-key slice of K. No canonical-residue conversion is claimed or charged: the
  decisions are bits of the representation the arithmetic already holds. x -> x R mod p is
  a bijection of Fp, so a uniform residue has a uniform Montgomery image and section 6's
  bit-field probabilities transfer unchanged; and T0M is a function of T0, so it is still
  equal on all 12 members of a class (part A checks this). The only ordering used compares
  two Montgomery residues as six-limb integers, a total order on Fp that is a function of
  the class, so section 3.1(ii) holds unchanged.
* Let zeta = xi^{(p^2-1)/6}, in Fp and a primitive 6th root of unity (zeta^2 = zeta - 1).
  Writing c_1 = a + b'u for the w^1 coefficient of y, the Re(w^1) of the conjugate y^{p^k}
  is `key_k = a*zeta^j` for k = 2j and `key_k = B*zeta^j` for k = 2j+1, where
  `B = Re(conj(c_1) d_1) = a*Re(d_1) + b'*Im(d_1)`. Since zeta^j in
  {1, zeta, zeta-1, -1, -zeta, 1-zeta}, the twelve keys are the two sign-symmetric triples
  {+-a, +-a*zeta, +-(a*zeta - a)} and {+-B, +-B*zeta, +-(B*zeta - B)}: all 12 keys come from
  **4 Fp multiplications** (a*zeta, B as two products, B*zeta) plus 8 subtractions or
  negations. No Frobenius of the whole element is computed. Inversion (k = 6) maps key a to
  -a, so the 12 keys are pairwise distinct for a generic element (part A); the minimum over
  a triple's six residues is min(m, p - M) with m, M the triple's minimum and maximum
  magnitude (a zero magnitude, probability <= 2/p per step, is treated as a tie, i.e. a
  failure). K and kappa therefore cost the two differences a*zeta - a and B*zeta - B, the
  two reflections p - M and the min/max reductions, not a twelve-way scan (section 5), and
  part A checks that this reproduces the 12-way argmin and key. Define
  kappa(y) = argmin_k key_k (over the Montgomery representatives, section 2) and the class
  key K([y]) = min_k key_k; then y^{p^{kappa(y)}} is the same element for every member of a
  class (the canonical representative y^), so K and kappa are class invariants. Ties need a
  1-in-p linear coincidence (probability <= 66/p per step, total <= 2^-249) and count as
  failure.

* **Cyclotomic squaring.** On GT the Granger-Scott formula squares y with three Fp4
  squarings and the linear maps ret = 3t - 2a on the a[0] half and 3t + 2a on the
  a[1] half (with xi*t2[1] for ret[1][0]); the shipped routine is measured, not
  modelled, in section 5, and part A checks that it equals y*y on random GT elements.
  Squaring commutes with
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
per step from the Montgomery value T0M of section 2; T0 has 381 bits, so they fit with
308 bits to spare). All later mentions of T0 mean T0M:

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
if the distinguished field of T0 is all zero:              [one DP attempt]
    increment this walk slot's cumulative DP-attempt counter r (every distinguished
    hit, stored or not) and add I_REP = 3,000 to the slot's work meter; stop the
    slot if r = R_w or if the meter reached W_cap
    compute the class key K (4 Fp mul) and append the 80-byte record
    (bits 70..197 of K, p^kappa*alpha, p^kappa*beta) unless it equals this walk's last
    record (duplicate suppression is a storage rule only)
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
the branch and the distinguished flag are read from the class invariant T0M, and the
multiplier table is fixed, so two walks at the same class take the same step. The two
cycle escapes are the only history-dependent actions, and each is a function of the
cycle reached rather than of the walker (section 3.1). So two walks that reach the same
cycle compute the same escape class and the same escape action, and from a common class
their class sequences are the same sequence. What is **not** bounded by one cycle is the
iteration on which each walk DETECTS a cycle: checkpoints are installed on the retained
step counter mod 64, which no escape resets, so a walk that enters a cycle just after
installing one detects it up to 64 + L <= 128 steps later than a walk that entered with
a checkpoint already inside it. Two walks on one cycle can therefore be out of phase by
up to one checkpoint period plus one traversal; their class sequences agree from their
first common class up to that time offset, and section 3.1 states exactly what a class
collision then still gives.

**Cycle comparisons compare T0.** On a squaring step no class key exists, so the 2-cycle
comparison, the checkpoint comparison and the escape's ordering all use T0 (plus nothing
else: 381 bits are enough). A false match needs two distinct classes with equal T0,
probability about 1/p = 2^-380.70 per comparison (log2 p = 380.7007 for the 381-bit
modulus); over at most 2 x 2^20 x 2^106.556025 < 2^127.556 comparisons (2 per step, for
fewer than 2^20 x SLOT_STEPS_MAX = 2^126.556025 steps in all: the hard step bound the
slot work meter allows, which already contains the walk budget S = 2^125.912895, the
detection window 2^84.33 and the restart waste 2^95.66) the total is below
2^127.556/p = 2^-253.14, i.e. below 2^-253
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
consecutive records that start at record 2^47 i, and its cumulative counter r counts
distinguished-point ATTEMPTS, not stored records. Every visit to a distinguished class
increments r, whether its record is stored or suppressed as the walk's last record, and
the slot is stopped at r = R_w, so a slot makes at most R_w attempts and the run at most
M R_w = 2^67. A suppressed attempt stores nothing, so a slot's stored records are at most
its attempts and the log of 2^67 record slots still cannot overflow. Duplicate suppression
is therefore a storage rule only and cannot hide work: the attempt counter bounds the
attempts, and each attempt is charged I_REP = 3,000 instructions to the slot's work meter
(section 5), which is what charges the K and the two p^kappa exponent multiplications that
an attempt computes before it decides whether the record is a duplicate. A stored record
goes to address base + 80 (2^47 i + r),
a 128-bit offset (two words: 80 x 2^67 exceeds 2^64) formed by a shift, an add and a
multiplication by 80. An append is one 80-byte
compare, one counter compare, the offset computation, one 80-byte store and one
increment, a fixed instruction count with no hashing, probing, chaining or eviction, and
no shared counter, so the append needs no atomic operation. The total number of appends is
at most the number of attempts, at most M R_w by construction: no restart resets r and the
slot stops at r = R_w = 2^47, so a slot never writes outside its own region, the log can
never overflow and no report-count probability enters the memory or processing bound. The slot's DP-attempt counter is one
64-bit word; 2^47 < 2^63, so its cap test can fire. (The toy of part B keeps each log as a
list of the filled records, which is the same set of records the matcher sees.)

**Report margin at the step budget.** Under H1 the number of distinguished classes a
slot visits is geometric with parameter theta: its mean over the S_slot + E_slot =
2^105.912895 + 2^64.321928 walk steps is 2^45.912895 attempts, so the attempt cap
R_w = 2^47 is a = 2^1.0871 = 2.1244 times the mean and the Bernoulli-KL (Chernoff) bound
gives log2 P(slot reaches the cap) <= -mu D(a||1)/ln 2 = -4.55 x 10^13. That is a
heuristic consequence of H1 (distinguished status is a class property, so the visits are
not independent), so it is not used as a hard tail: it is carried in the loss list of
section 4 at the same order of magnitude the accepted 125.826 package carried, and the
capability (attempts, hence log <= M R_w) is deterministic. Counted as attempts, the cap also bounds a repeated report's computation: a
fruitless cycle that revisits one distinguished class spends an attempt per pass, so at most
R_w are spent per slot and each is metered, whatever the restart counter does.

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
detection iteration can differ, by up to 64 (the checkpoint period) + L steps, because
checkpoints are installed on the retained step counter mod 64, which no escape resets.
Section 3.1 states that bound and section 4.1 charges it; no one-cycle coalescence is
claimed. Detection of an L-cycle (L <= 64) happens
within 64 + L steps of entering it, so the waste per event is <= 64 + 2L + 1 <= 193 <=
200 steps; a 2-cycle event wastes <= 3 steps. An escape is one squaring step (cheap) or
one squaring plus one re-step, and section 5 charges every one of those steps at the
multiply-step price, i.e. above the squaring price.

**Stuck-walk rule.** Each walk counts steps since its last distinguished class in a
128-bit counter (reset at each distinguished hit and at restart; this is the trajectory
counter, not the slot's DP-attempt counter r, which no restart resets; 40/theta = 2^65.33
exceeds 2^64). When
the count exceeds 40/theta = 2^65.33 the walk is abandoned and restarted at a fresh
g^alpha z^beta. Its steps are real steps and are inside the slot's work-meter cap of
section 5. It is the safety net for cycles longer than 64 that contain no distinguished
class (section 3.2 bounds the probability that cycles matter to success). Under H1 the
steps between distinguished classes are geometric with parameter theta, so a run exceeds
40/theta with probability (1 - theta)^(40/theta) < e^-40; over at most 2^66 runs the
expected number of restarts is < 2^66 e^-40 < 2^8.3, and by Markov's inequality more than
2^30 restarts occur with probability < 2^-21. 2^30 restarts waste at most
2^30 x 2^65.33 = 2^95.33 steps, a fraction < 2^-30 of S.

### 3.1 Coalescence up to the checkpoint period (what a class collision gives)

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
action. Their class sequences are the same sequence from their first common class, but
they are TIME-SHIFTED: a walk detects a cycle only when its own checkpoint ring notices
the return, and the 64-step checkpoint schedule runs on the retained step counter, which
neither a ring nor a checkpoint reset resets. A walker that enters an L <= 64 cycle just
after installing a checkpoint detects the return at the next ring pass, up to 64 + L <=
128 steps later than a walker that entered with a checkpoint already inside the cycle;
one escapes and continues while the other is still lapping. Two walkers on one cycle are
therefore the same class sequence up to an offset of at most one checkpoint period plus
one traversal (<= 128 steps); after each escape both hold an empty ring and an empty
checkpoint slot, so by (i) their steps are class functions again and a later cycle can
again differ only by that offset. That is the only coalescence claimed: no one-cycle
detection bound and no simultaneity, and part P measures the offset on constructed cycles.
The recovery does not need the timing at all (iii): two records of one class with different
exponent pairs already give x, and the matcher is offline, so an offset costs no pair and no
work beyond the steps already inside each walker's own meter.

**(iii) A class collision is still a usable pair.** A report at class C stores (T, A, B),
where T is the truncated class key K(C) and (A, B) = (p^kappa alpha, p^kappa beta) are the
exponents of that walker's element. Two walkers that meet at C produce the same T. If their
exponent pairs differ, the pair solves p^k1(alpha1 + beta1 x) = p^k2(alpha2 + beta2 x), so
x = (A1 - A2)/(B2 - B1) mod r, as in section 3. If the pairs coincide the two records are
byte-identical: then both walkers hold the same element g^alpha z^beta at class C.
Coincident records cannot hurt: the
heapsort puts identical records adjacent and the adjacent dedup removes them before any
pair is enumerated, and the append rule stores no repeat of a walker's own last record, so
a division by zero from A1 = A2 cannot arise. For two records in one class, B1 = B2 holds
only in that case, since a Frobenius shift multiplies B by p^k = 1 mod r only for k = 0.
The remaining unsolvable case is a false-key pair (distinct classes whose truncated keys
agree) with B1 = B2; it gives no x, is charged as the 1/r term of section 4's loss list,
and the matcher continues with the next pair.

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
S = 2^125.91 steps of the run: they are not an event of probability 2^-20, they are the
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
useful.** A revisit that is not (i) or (ii) reaches a class already visited; if its
exponent pair differs, that is two records of one key T with different (A, B), i.e. the
collision the algorithm wants, x = (A1 - A2)/(B2 - B1) mod r. A revisit with an equal
exponent pair means the walk is on a class cycle longer than the checkpoint window; the
conditions 2^t = p^e mod r and sum_i p^{t_i} = 0 mod r of section 3.2 force such a cycle to
be all-multiply with length above 64, and bound the chance a walk enters one by
S x 2^-311.8 < 2^-185, inside the long-fruitless-cycle term of section 4.

(iv) **Consequence for the failure ledger: a within-walk return is escaped, useful, or a
wasted slot.** (i) and (ii) are inside the charged cycle overhead and case (iii) with
different exponent pairs is the useful collision; in the other cases the walk can spend the
rest of the slot's budget on repeats instead of fresh classes. Section 4(ii) charges that
whole self-collision event at n^2/(2N') per slot as an upper bound, so no birthday term is
needed anywhere else, and the earlier claim that the record-count argument needs the
absence of a within-walk return stays withdrawn: the log's capacity is deterministic. Each
slot owns R_w = 2^47 record positions and stops at its cap, so the log holds at most
C_rep = M R_w = 2^67 records whatever the walk does, with no overflow and no probability
(section 3); a lap that adds records of one class to one run is charged as pairs of that
run like any other, against the hard pair budget P = 2^41 whose own abort probability is
the 2^-36 Markov charge of section 3.

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

  S0' = sqrt(2 N' ln(1/0.61)) = 2^125.884326,   S = 1.02 S0' = 2^125.912895,

where the 2% is explicit slack for walk non-randomness (the ceil of the per-slot budget
adds another 1e-9 relative, 1.4e-9 bits); at S0' the model gives exactly 1 - 0.61 =
0.39, and at S it gives 1 - 0.61^(1.02^2) = 1 - 0.61^1.0404 >= 0.40206 (rounded down). The algorithm
stops after S walk steps (plus the detection window below); S is a fixed budget, not an
expectation.

**Why 2% is enough slack.** The 0.01 H1 allowance is the slack for the heuristic parts of
H1 and is kept; the branch-adjusted collision law that sets N' is measured in part E, part
B's toy fractions support the model at S0' and just above it, and the losses sum to
3.342 x 10^-6, three orders below what the 2% contributes. The 2% raises the model from
0.39 to 0.40206, covering the allowance and the 2^-18 of losses with 0.0020 left; a larger
slack buys back part of the allowance at 0.0141 bits per 2% of steps.

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

**Detection.** After the first useful collision the two walkers are the same class
sequence up to an offset of at most one checkpoint period plus one traversal (<= 128 steps,
section 3.1), and from their first common class the wait to a distinguished class is that of
one trajectory: within 20/theta steps except with probability e^-20 (H1). The window is applied to a trail selected through a collision, not to a fresh
trail; the argument does not need the fresh-trail form, because under H1 the distinguished
indicator of each step is independent of the past with probability theta, so the residual
wait of the selected trail is geometric with the same parameter (the geometric tail of a
fixed-parameter wait) and its tail is e^-20 as well. No length-biased trail allowance is
claimed. Charging all M walks in lockstep over that window costs E = M (20/theta + 128) =
2^84.33 + 2^27 steps: the 128 is the offset above, charged although the matcher is offline
and the offset needs no step of its own; it is 2^-98.9 of S, below 2^-40 of S and far below
the 1% meter band.
The losses, each a probability bound: detection tail e^-20 < 2^-28; trailing walk
restarted before reaching the distinguished class (its distinguished-free run must
exceed 40/theta: under H1 the run containing a given step exceeds 40/theta with
probability (1 + 40) e^-40 < 2^-52); self-collisions of a single walk (fraction
<= 2^-19 of collisions, with M = 2^20 walks); a false-key pair with B1 = B2 (1/r); key ties (< 2^-249);
long fruitless cycles (< 2^-100, now only all-multiply and at the mixed rate of section
3.2, so the bound only improves); more than 2^30 stuck-walk restarts (< 2^-21);
more false-key pairs than the hard pair budget (< 2^-36 by Markov, section 3); a slot
that collides with itself early and wastes its budget (<= 2^-20 over the run, section
4(ii)), whose success loss is <= 2^-21.66; a slot aborted by its meter, i.e. its
productive work above 1.01 of nominal or its cycle charge exceeded (< 2^-100 over the
run, section 4(v)); a walk slot stopped by the per-slot DP-attempt cap,
i.e. a healthy slot reaching R_w = 2^47 attempts (the KL figure of section 3, < 2^-300);
two distinct classes with equal T0 causing a spurious escape (< 2^-253). Within-walk
returns are charged twice over: an escaped return sits in the cycle overhead of section
3.2, an unescaped one is the wasted slot above, and a return with a different exponent
pair is the useful collision. Their sum is below
3.342 x 10^-6 < 2^-18 < 0.4020 - 0.39 - 0.01 = 0.0020.

**The work meter, accounted per step type.** A step costs c_sq = 0.50490 or c_mul = 1.1125
units, so the walk phase is metered and a slot stops at W_cap (section 5). The 1% band is
applied to the slot's productive steps alone, because the branch bit is a deterministic
function of the current class: a revisit repeats an earlier decision, step costs are not
i.i.d. across time, and Hoeffding cannot be applied to the whole sum. The slot's work is
productive-step costs + repeat-step costs + cycle-handling steps.

(i) **Productive steps are the i.i.d. part.** A productive step is one on a class this walk
has not visited. Under H1 the class sequence is a random mapping, so a newly exposed class
is not a function of the walk's history and its branch bit - a field of T0, measured on
real GT elements in part A2 - is a fresh draw: a Bernoulli(35/64) choosing between two
bounded costs in [c_sq, c_mul] with mean c = 0.78021875. Conditional on the step being
productive the costs are therefore i.i.d., and this is the only concentration inequality
used. The fresh-class law is part of H1, not a proof: part A2's bit frequencies, part P's
production rule on real elements and part B's toy walks support it, and the 0.01 allowance
below is slack for it.

(ii) **A slot that collides with itself can waste its budget; it is charged twice.** Under
H1 a slot's own class sequence can return to a class it already holds. If that happens
early enough the walk enters a cycle longer than the 64-step checkpoint window: the
checkpoint rule does not escape it, a lap need not leave it, matching is offline, and the
no-report restart may be held off by reports from the cycle, so the walk keeps revisiting
classes instead of exposing fresh ones. Under H1 the probability that a walk of n =
2^105.913 steps meets its own class again at all is at most n^2/(2N') =
2^211.826/2^252.785 = 2^-40.959 per slot (the union over its step pairs), and a
self-collision within its first O(1) steps has probability O(1/N') = 2^-251.785,
negligible. Such a slot wastes at most its remaining n steps, and then its own meter, its
DP-attempt cap or a restart ends it. It is charged explicitly on both counts:
  (a) **Work.** The wasted steps are the slot's own step budget, already inside S and
  inside its meter, and the meter stops the slot at W_cap, so the hard walk bound
  W_walk <= M x W_cap of section 5 holds unchanged and charging the steps again would
  count them twice. As a figure: at most M x n^2/(2N') = 2^-20.959 slots in expectation,
  so the wasted work is at most (M n^2/(2N')) n = 2^-40.959 S = 2^84.95 steps in
  expectation, against the run's M x 2^80 = 2^100-step revisit allowance, and at most
  2^-20.96 S with probability >= 1 - 2^-20 (Markov).
  (b) **Success.** A slot that spends its budget on repeats exposes no fresh classes, so
  the productive count is T >= S - L with E[L] <= 2^-40.959 S, and L < 2^-20.96 S except
  with probability <= 2^-20 (Markov). A shortfall of that fraction lowers the collision
  parameter by (1 - 2^-20.96)^2 and the success probability by at most 2^-21.66, the
  entry below; the 2^-19 self-collision loss already reserved is kept, since a wasted slot
  also loses the pairs it would have formed. Per slot the loss is at most its own n steps
  of fresh classes, which is what the two charges bound.

(iii) **Cycle handling is escaped work, not a productive step.** Fruitless 2-cycles and
checkpoint cycles are detected by the T0 comparisons, which run before the report append,
and escaped (section 3.2); they are charged in the overhead factor OVER = 3/(12R) + 200/R^2
= 4.35 x 10^-4 steps per step. The measured mixed rates are 3 (1-s)^2/(12R) = 5.01 x 10^-5
and 200 x 0.0101 R^-2 = 1.93 x 10^-6, i.e. 5.21 x 10^-5 actual against 4.35 x 10^-4
charged: an 8.35x safety factor. A slot's expected cycle work is 2^91.68 steps, and a
Chernoff bound for the count to beat the 8.35x charge gives 2^(-2^94.74) per slot (part C).
That count is history-dependent, so it is stated as H1-conditional and is not in the i.i.d.
sum.

(iv) **The productive sum.** With n = N_slot = 2^105.91 steps per slot and
c_mul - c_sq = 0.6076, Hoeffding gives P(productive work > 1.01 c n) <=
exp(-2 (0.01 c)^2 n / (c_mul - c_sq)^2), with exponent 3.298 x 10^-4 nats per step: the
probability is 2^(-2^94.88) for a slot and 2^(-2^114.88) over the run's S steps. Those
94.88 and 114.88 are the exponents of the exponent, not probabilities; part C prints both
and checks the form.

(v) **The failure ledger for aborts.** A slot aborts only if productive work exceeds 1.01
of nominal (2^(-2^94.88) per slot) or the cycle charge is beaten (2^(-2^94.74) per slot);
their union over the M = 2^20 slots is below 2^-100, the figure in the loss list. An
aborted slot loses at most its remaining productive steps, so the success loss is at most
the fraction of aborted slots, which is zero outside that event. The wasted-slot event of
(ii) is charged separately, at 2^-20 over the run. There is no "S/1.01 steps" equivalence
and none is used.

Together: (i) and (iv) bound the productive work on fresh classes, (a) and the hard meter
cap bound the work whatever the walk does, and (b) bounds the productive count the success
model runs on.

Success >= 0.40206 - 0.01 - 2^-18 > 0.40206 - 0.01 - 0.00001 = 0.39205, so
success >= 0.3920 >= 0.39 (every figure rounded down), with 0.0020 of margin. The 0.01
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

The squaring is **measured, not modelled**. `measure/` ships the harness, the raw Callgrind
records, the exact totals at n = 1 and n = 65 over five slice offsets and a disassembly
audit; `measure/README.md` documents the toolchain, the instrument and the reproduction, and
`experiment.py` re-derives every number below from those records and fails if one differs.
blst's `blst_fp12_cyclotomic_sqr` at blst f7a88067 costs **9,803.0 executed instructions per
call** ((637,215 - 9,823)/64), identical at all five offsets on five GT elements. The ledger
charges the body, 9,803 - 5 = **9,798**; the 5 are harness loop control already inside the
common rows. The count cannot depend on the input: the call tree is 18 functions and 1,017
static instructions with no conditional, backward or indirect branch (its 16 indirect
instructions are `ret`s), and on GT the routine is `blst_fp12_sqr`. The same harness
reproduces blst_fp_mul = 383.0 and blst_fp12_mul = 22,996.0 exactly, which calibrates the
instrument; a static sum over the disassembled bodies (9,619) is 1.88% low because it omits
call and return plumbing. The measurement is what is charged.

The charge is 10,098 = 9,798 measured + 300 for the doubling, the shared rows and slack (the
convention the earlier 11,700 = 11,400 + 300 used); the counted rows are 9,951, so the
budget is 147 above the count. A squaring step computes no class key. Squaring-step cost
c_sq = 10,098/20,000 = 0.5049 units, and with s = 35/64 the average step cost is
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
  contention on it, and the log's per-slot regions leave the walk phase no shared counter.
* Values. W_cap = ceil(1.01 x nominal) + 2^80 x I_MUL, with nominal = c x C x (S_slot +
  E_slot) + I_MUL x ceil(OVER x S_slot) + I_REP x R_w =
  1192541993987766958416474245181635120 = 2^119.843452 instructions: the third term is the
  worst-case DP-attempt charge, 2^47 attempts at 3,000 instructions, 2^-61.29 of the
  nominal; the 2^80-step revisit allowance adds 2^94.44 instructions, so
  W_cap = 1204467440826244114426138124845867472 = 2^119.857807. The meter's largest value
  is W_cap - 1 + I_MUL + I_REP = 1204467440826244114426138124845892721 < 2^127, so it fits
  two 64-bit words (a step that is also a DP attempt adds both charges). A slot is stopped
  when its meter reaches W_cap, so W_walk <= M x (W_cap - 1 + I_MUL + I_REP) =
  2^139.857807 instructions = 2^125.570095 units is a **hard bound**, not an expectation;
  it contains every step and every DP attempt, and it is what the claim is computed from
  (the attempt term raises the bound by 2^-108.31 relative, part C checks the integers).
* Report attempts are inside that meter. The charge I_REP = 3,000 instructions is added to
  the slot's meter at every distinguished hit, stored or suppressed, so a fruitless cycle
  that keeps revisiting one distinguished class pays 3,000 instructions per pass inside the
  same hard cap as its steps, whether or not the no-report restart counter is reset; at
  most R_w = 2^47 attempts are possible per slot, and even that worst case is 2^-54.65 of
  the 1% band the cap carries above the nominal, so report work cannot consume the band.
  The same work is charged again in the capacity group below (3,000 of the 13,100 per cap
  unit): the group is a worst case over the log whether or not the meter stops the slot
  first, and charging it twice is deliberate.
* Abort bound. Section 4(v) gives the productive Hoeffding bound 2^(-2^94.88) per slot
  and 2^(-2^114.88) over the run (94.88 is the exponent of the exponent, not the
  probability) and the cycle charge's Chernoff bound 2^(-2^94.74) per slot; their union
  over the M slots is below the 2^-100 of the loss list. A slot that collides with itself
  cannot pass the meter either: it stops at W_cap with fewer productive steps, so the
  bound above is a work bound whatever the walk does, and the 2^-40.959-per-slot event
  costs fresh classes, charged in section 4(ii)(b). The 0.01 success allowance covers the
  heuristic parts of H1, not an assumed distribution of the work.

**Counters and indices.** Every counter and index, its largest value at the end of the
attack (including one over-cap update per processor), and its width. Part C recomputes
this table and asserts every row.

| counter or index | largest value (log2) | words | capacity test |
| --- | ---: | ---: | --- |
| slot work meter (instructions) | 119.857807 | 2 x 64 | >= W_cap with W_cap = 2^119.857807 < 2^127 |
| slot step counter | 106.556025 | 2 x 64 | none (informational) |
| steps since last DP attempt | 65.321929 | 2 x 64 | > 40/theta (2^65.321929 rounded up), a 128-bit compare |
| slot DP-attempt counter r | 47 | 1 x 64 | >= 2^47 with 2^47 < 2^63 |
| log offset 80 x (2^47 slot + r) | 73.321929 | 2 x 64 | < 2^128 |
| slot restart counter | 41.234097 | 1 x 64 | < 2^63 |
| matcher pair counter | 41 | 1 x 64 | > 2^41 with 2^41 + 1 < 2^63 |
| heapsort index | 67 | 2 x 64 | < 2^128 |

No counter's test sits at its word's own capacity: the two that can fire are at
2^119.857807 of a 2^128 word and 2^47 of a 2^64 word.

**Fruitless-cycle overhead factor (a hard cap on cycle-handling steps).** The step meter
also carries the cycle-handling steps: the nominal charged work adds
OVER = 3/(12R) + 200/R^2 = 4.3488 x 10^-4 extra steps per step at the multiply-step price
(2-cycles at rate 1/(12R), <= 3 steps each; long cycles at rate R^-2, <= 200 steps each).
The branch makes the 2-cycle rate (1-s)^2/(12R) and the long-cycle rate 0.0101 R^-2, so the
charged rates are 4.87 and 99 times the mixed ones; cycle handling is inside the same hard
meter, and the nominal charge is at most 1.000621 times the plain step work, not an
average. Section 3.2 bounds the all-multiply long-cycle rate it is built on.

**Setup, output and preprocessing.** Setup and output are all <= 2^31 units together:
pairing g (bounded by 2^14), input check z^r = 1 (<= 510), table: 1024 multipliers x (2
exponentiations <= 1,021 + 11 Frobenius maps) <= 2^20.02, M = 2^20 walk starts <= 2^30.0,
output [x]P0 (<= 510 G1 operations) and verification pairing (<= 2^14). **Preprocessing and
table construction are at most 2^57 group-operation equivalents**: the table, the walk
starts and the initialization of all C_rep = 2^67 log records with the sentinel record,
2^67 x 12 = 2^56.30 units. Every term is charged inside W (the fill in the report group
below, the rest in its 2^31 entry); nothing is amortised and no advice is used.

**Restart starts**, charged at the worst case the meter allows: a slot's step counter can
reach 2^106.556025, and every restart needs 40/theta = 2^65.33 steps, so a slot has at
most 2^41.234097 restarts; each costs one exponentiation (1,021 units, 2^24.29
instructions) plus 100, so M x 2^41.234097 x (1,021 x 20,000 + 100) = 2^71.234 units.

**Distinguished-report handling**, charged at its worst case (the deterministic
ATTEMPT cap 2^67, not the mean), all with fixed instruction counts per cap unit (section 3):

| Work (2^67 cap units: attempts) | Instructions |
| --- | ---: |
| one DP attempt, stored or suppressed: the class key computed at the hit (4 Fp multiplications 1,532, the addition and the two subtractions 133, the min/kappa rows 258 and the w^1 load 16, together 1,939, rounded up to 2,000), two modular multiplications by p^kappa (512-bit products, reduction, <= 400 each), 80-byte compare against the walk's last record, 128-bit record offset 80 (2^47 slot + r) (shift, add, multiplication by 80), 80-byte store, per-slot counter increment | 3,000 |
| full-record heapsort: <= 67 levels x 130 (10-word compare-and-swap of the 80-byte record), charged up from the computed 8,710 | 10,000 |
| adjacent dedup: one 80-byte comparison per record | 100 |
| **total per record** | **13,100** |

2^67 x 13,100 = 2^80.68 instructions = 2^66.39 units: the 2^67 cap units are ATTEMPTS, the
3,000 covers each attempt including the ones whose record is suppressed, and the
10,000 + 100 cover the sort and dedup of at most 2^67 stored records. The pair term is bounded by the hard
pair budget P = 2^41: P x 200 = 2^48.65 instructions = 2^34.36 units, plus one canonical
recomputation (1,021 units) per record a charged pair touches, at most 2P records < 2^52
units; the sentinel fill is 2^67 x 12 = 2^56.3 units. The restart-starts term (2^71.234)
dominates the group and setup/output (2^31) is far below it: together **2^71.293 units**,
charged in full whether or not the reports occur (the capacity is deterministic).

W <= 2^125.570095 + 2^71.293 units = 2^125.570095 (every figure rounded up: the walk
bound is 2^125.5700944 to seven decimals, and the second term raises log2 W by
2^-54.29/ln 2 < 2^-53.8, below the last printed digit). Closed form, each factor rounded
up: W <= 1.01 x 0.78021875 x 2^125.912895 x 1.000621 + 2^94.44 <= 2^125.570095 (meter cap
x average step cost x step budget x cycle overhead, plus the revisit allowance; the log2 of
the four factors is 0.014356 - 0.358049 + 125.912895 + 0.000895 =
125.570097, and the allowance adds 2^-25.4 relative, so the exact value is 125.5700944).

Submitted bound: **time_log2 = 125.571** = 125571/1000 >= 125.570095 >= log2 W (part C
checks it as an exact integer inequality, W^1000 <= C^1000 x 2^125571, and that 125570
fails, so 125.571 is exactly the ceiling of log2 W at three decimals).

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
  other bits are; T0 = Re(c_0) is not a uniform Fp element, so these are heuristic. The
  fields are read from the Montgomery value T0M = T0 R mod p (section 2), the bijective
  image of T0, so the assumption is the same in either representation.
  Support: part A2 measures the frequency of each of the 76 used bits on real GT elements,
  part P runs the production rule on real GT elements, and part B's toy walks measure the
  fractions against their nominal values. The fields are disjoint bits of one value, so a
  defect in one cannot alias into another, and the same assumption prices the 128-bit
  report key and the report-cap margin.
* **Coalescence only up to the checkpoint period.** The ordinary step is a function of
  the class, but the two cycle escapes read stored history, and the checkpoint phase is a
  property of the walker's step counter, not of the cycle. Section 3.1 proves what survives:
  escape class and action are functions of the cycle, and two walkers on one cycle follow
  the same class sequence up to a time offset of at most 64 + L <= 128 steps, so section 4.1
  charges M (20/theta + 128) steps for the detection window (2^-98.9 of S) rather than
  claiming a one-cycle coalescence. The recovery needs none of it: two records of one class
  with different exponent pairs already solve for x (section 3.1(iii)), and the matcher is
  offline. Part P measures both cases on constructed cycles: a 2-cycle detected by the
  two-step ring (escape gap 2) and a 3-cycle detected only by the checkpoint (gap 63 of the
  67 the bound allows, both walkers detecting the same cycle).
* **Every distinguished hit is a charged attempt.** The attempt counter is capped at R_w per
  slot and each attempt costs 3,000 instructions inside the slot's meter and again inside
  the 2^67-cap capacity group, so a repeated report inside an unescaped fruitless cycle pays
  per pass and neither the no-report restart counter nor duplicate suppression can move it
  out of either bound: duplicate suppression is a storage rule only.
* **Preprocessing and advice.** At most 2^57 group-operation equivalents, including the
  2^56.30-unit initialization of all 2^67 log records, the multiplier table (<= 2^20.02) and
  the walk starts and setup/output (<= 2^31). All of it is inside W (section 5), so no
  preprocessing is amortised; nonuniform advice is 0 bytes.
* **The work meter's dependent-process accounting.** Step costs are not i.i.d. across time,
  because the branch bit is a function of the current class. Section 4 splits the slot's
  work into productive steps (i.i.d. costs on newly exposed classes, where H1's
  random-mapping law does supply a fresh draw, bounded by Hoeffding), cycle handling (escaped
  work, charged at 8.35 times its measured mixed rate) and the self-collision event, which
  spends a slot's remaining budget on repeats with probability at most n^2/(2N') = 2^-40.959
  per slot and is charged as a 2^-20 lost-slot and 2^-21.66 success-loss entry. The meter's
  bound on total work is hard (section 5) and holds under that event too, and section 4(v)'s
  abort events are below 2^-100 over the run. What remains heuristic - the fresh-class law
  and the cycle count - is supported by part A2's measured bit frequencies and the toy's
  measured per-slot work deviation, 0.3955/0.4094/0.3864 sqrt(steps) against 0.3877 for
  independent Bernoulli step costs. No "S/1.01 steps" equivalence is used.
* **No production run.** The toys are evidence only at their own scale: distinguished
  fraction 2^-4 or 2^-3 rather than 2^-60, field near 2^20 rather than 2^381, class count
  2^21.4 to 2^22.4 rather than 2^251.3, T0 20 bits rather than 381. The production cycle
  rate is the derivation of section 3.2, not a toy row.
* **The squaring price is a measured instruction count of one shipped routine**, not a
  model of the formula: the Callgrind count of blst f7a88067 `blst_fp12_cyclotomic_sqr`
  (9,798 per call, plus 300), on the toolchain the cost model names, calibrated against the
  cost model's own 383.0 and 22,996.0 and re-derived by part C from the shipped record. Its
  straight-line body makes the count input-independent by construction; the residual risk is
  different code generation on the grader's compiler, which is why the charge keeps 147
  instructions of slack.

* This is a generic attack (about the square root of the class count). It ignores
  index-calculus methods in Fp12 (TNFS family); it is a valid upper bound, not a
  security estimate.

## 7. Experiments (`experiment.py`, standard library only)

`experiments.md` gives the per-part detail; `measure/README.md` and
`measure/experiment-results.md` hold the measurement record. In summary:

* **Part A** (real BLS12-381 GT elements): GT membership, the Frobenius formula, the key
  formula against all twelve conjugates, the sign-set minimum against the twelve-way
  argmin, p^6 = inversion, ord_r(p) = 12, the factorisation of r-1 and
  ord_r(2) = (r-1)/16 with the 2^247.27 squarings a fruitless cycle needs, T0 equal on all
  twelve conjugates of the Montgomery value, and cyclotomic squaring = y*y on GT.
* **Part A2** (Montgomery domain): one-bit frequencies of the 76 used bits of T0M and of
  the 128 report-key bits of K on up to 2,000 random GT elements, the branch and
  distinguished fractions against 35/64 and 2^-4, and the fraction of elements whose
  branch decision differs between the canonical and the Montgomery reading.
* **Part P** (production step rule and both cycle escapes on real GT elements): identical
  T0M sequences from two conjugate starts, the exponent invariant, key computations equal
  to multiply steps plus attempts, a constructed two-cycle detected by the two-step ring,
  and a constructed three-cycle whose detection times differ by the 64 + 3 steps section
  3.1 allows.
* **Part E**: the branch-adjusted collision law measured directly in the toy field q = 67
  by pair type, with squaring-squaring coincidences exactly zero.
* **Part B**: the complete mixed walk, attempt cap, meter and offline matcher in three toy
  fields near 2^20, with the pooled success fraction against the model.
* **Part C**: the bound arithmetic of sections 3-5 recomputed from integers, including the
  integer claim check, the attempt charge, the detection offset and the squaring price
  re-derived from `measure/`.
