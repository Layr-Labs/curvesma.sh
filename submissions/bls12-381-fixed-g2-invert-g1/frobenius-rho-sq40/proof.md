# Frobenius-class parallel Pollard rho in GT with a cyclotomic-squaring branch (success target 0.39)

Claim: time_log2 = **125.612** (total work W <= 2^125.612 dense-Fp12-multiplication
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
  outputs in that form. Every per-step decision in section 3 is read from the word the
  measured kernels produce: T0M = T0 R mod p for the partition, branch and distinguished
  fields, both cycle comparisons, the "smaller-T0" and "smallest-T0" escape orderings and the
  report-key slice, and the twelve keys in that same form for the class-key ordering that
  defines kappa. The four kernel multiplications that form a*zeta, the two for B and B*zeta
  read Montgomery operands and write Montgomery outputs, so key_k R mod p is what the
  ordering sees, at the same four Fp multiplications; no canonical-residue conversion is
  claimed or charged. R-multiplication is a bijection of Fp, so the keys are still twelve
  distinct values in the same sign-symmetric triples and section 5's min/max reduction is
  unchanged, but the selected kappa (hence the canonical representative and the twist t of
  every multiply step) can differ from the canonical-residue choice. T0M is a function of T0,
  equal on all 12 members of a class (part A checks it), so comparing two Montgomery residues
  as six-limb integers is a class function, which is what 3.1(ii) needs. The keys'
  selected bits are not uniform (section 6).
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
package can state exactly: the continuous minimiser of log2 W is s = 0.5464, and on the
grids 2^-k the best value is 35/64 (s = 273/500, which a bit field cannot express, gives
log2 W lower by 2.04 x 10^-7, the same three-decimal claim). The reason there is an interior
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
each is exact: 35/64 for the branch, 2^-60 for the distinguished flag. Part A2 measures their frequencies (section 7).

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
modulus); over at most 2 x 2^20 x 2^106.597846 < 2^127.598 comparisons (2 per step, for
fewer than 2^20 x SLOT_STEPS_MAX = 2^126.597846 steps in all: the hard step bound the
slot work meter allows, which already contains the walk budget S = 2^125.954715, the
detection window 2^84.33 and the restart waste 2^95.66) the total is below
2^127.598/p = 2^-253.10, i.e. below 2^-253
(a loss term of section 4). Its effect is only a spurious escape
(a wasted step) and never a false logarithm: the report key is derived from K, not from
T0, so two classes that share T0 but differ in K are never matched, and the escape's
choice of class is a function of the T0 multiset of the cycle (section 3.1), so a
spurious escape only wastes steps and cannot desynchronise two walkers that reach the
same cycle.

**Report log (deterministic capacity, direct addressing, per-slot cap).** The log is exactly
C_rep = M R_w = 2^20 x 2^47 = 2^67 record slots of 80 bytes, allocated before the walks start
and filled with one fixed sentinel record (all bytes 0xFF, exponent fields above r, so the
matcher skips it; 12 instructions per record, charged in section 5). Walk slot i owns the
2^47 consecutive records from record 2^47 i, and its cumulative counter r counts
distinguished-point ATTEMPTS, not stored records: every visit to a distinguished class
increments r, stored or suppressed, and the slot stops at r = R_w, so a slot makes at most
R_w attempts and the run at most M R_w. Suppression is therefore a storage rule only and
cannot hide work: a slot's stored records are at most its attempts, the 2^67 slots cannot
overflow, no restart resets r, every attempt is charged I_REP = 3,000 instructions to the
slot's meter (section 5) for the K and the two p^kappa exponent multiplications it computes
before the duplicate test, and no report-count probability enters the memory or processing
bound. A stored record goes to address base + 80 (2^47 i + r), a 128-bit offset (80 x 2^67
exceeds 2^64) formed by a shift, an add and a multiplication by 80. An append is one 80-byte
compare, one counter compare, the offset computation, one 80-byte store and one increment, a
fixed count with no hashing, probing, chaining or eviction and no shared counter, so it needs
no atomic. The DP-attempt counter is one 64-bit word and 2^47 < 2^63, so its cap test fires.
(The toy of part B keeps each log as a list of the filled records, the same set the matcher
sees.)

**Report margin at the step budget.** Under H1 the number of distinguished classes a
slot visits is geometric with parameter theta: its mean over the S_slot + E_slot =
2^105.954715 + 2^64.321929 walk steps is 2^45.954715 attempts, so the attempt cap
R_w = 2^47 is a = 2^1.0452 = 2.0637 times the mean and the Bernoulli-KL (Chernoff) bound
gives log2 P(slot reaches the cap) <= -mu D(a||1)/ln 2 = -4.24 x 10^13. That is a
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
70..197 of K of distinct classes independent, uniform under H1) two records of distinct
classes agree on T with probability <= 2^-128, the log holds at most 2^67 records, so
there are at most (2^67)^2/2 = 2^133 pairs of distinct records and the expected number
of false-key pairs is <= 2^133 x 2^-128 = 2^5. Hence P(pair count above 2^41) <=
2^5/2^41 = 2^-36, which is a loss term of section 4.

A GENUINE repeat -- the same canonical element reached with different exponent pairs, which
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
detection bound and no simultaneity; part P measures the offset on constructed cycles.
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
only reduce fruitless-cycle work, and the charge does not depend on that.

**Class cycles, within-walk returns, and what a lap can cost.** If a cycle of length L
contains a squaring step, then after one lap the class returns but the exponent pair is
multiplied by 2^t p^m with 2^t != p^e, so the pair is different: two records with the
same key T and different (A, B), which is exactly a solved instance. A lap can therefore
enlarge an equal-key run of the report log, and the only place that hurts is the hard
pair budget P = 2^41 of section 3. The accounting below distinguishes every kind of
revisit explicitly, because within-walk returns are common and must not be mistaken for
a rare failure.

(i) **Algebraic fruitless 2-cycles are handled events, not failures.** They occur at rate
(1-s)^2/(12R) = (29/64)^2/12288 = 1.671 x 10^-5 per step (part B measures 1.92 x 10^-5 at
R = 1024, part C computes the rate), i.e. of order 2^110.1 occurrences over the
S = 2^125.96 steps of the run: they are not an event of probability 2^-20, they are the
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

  N' = N / (1 - s^2 - (1-s)^2/R) = 2^251.785203

where the charge uses 1/R instead of 1/(12R) for the same-index case: same-partition
pairs are charged as never merging even when their twist index differs, which only
enlarges N' (by 2.624 x 10^-4 relative, 3.784 x 10^-4 in log2: the factors differ by
(1-s)^2 (1/R - 1/(12R)) = 1.839 x 10^-4 against 0.7007). Then
the probability that no useful collision occurs among the first S
distinct-step points of all walks is exp(-S^2/(2N')). Choose

  S0' = sqrt(2 N' ln(1/0.61)) = 2^125.884326,   S = 1.05 S0' = 2^125.954715,

where the 5% is explicit slack for walk non-randomness (the ceil of the per-slot budget
adds another 1e-9 relative, 1.4e-9 bits); at S0' the model gives exactly 1 - 0.61 =
0.39, and at S it gives 1 - 0.61^(1.05^2) = 1 - 0.61^1.1025 >= 0.42013 (rounded down). The
algorithm stops after S walk steps (plus the detection window below); S is a fixed budget,
not an expectation, and every toy trial was given exactly this budget.

**Why 5% is enough slack.** Part E measures the collision law that sets N'. Part B's three toy
fields pool 1,224 of 3,000 trials colliding within 1.05 S0' = 0.4080 +- 0.0090: 0.971 of the
model, 1.35 sigma below it. Taking the 0.01 allowance off the measured value leaves 0.3980 >=
0.39, with 0.0080 to spare (0.9 sigma), and each field alone leaves 0.391, 0.404 and 0.399.
That is evidence at the toy's scope (fields near 2^20), not a bound. The losses consume
0.007453 of the allowance (4(v)), and the claimed success after it is
0.42013 - 0.01 - 2^-18 = 0.41012 >= 0.39, margin 0.02012. The R = 64 row alone gives 0.401 +-
0.015 (0.391 after the allowance), within its error of 0.39, so nothing is claimed there.

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

So the per-pair merge probability is (1 - s^2 - (1-s)^2/(12R))/N, exactly the factor
above. Part E measures this law on toy class sets (distinct classes, one step each,
coincidence counts by pair type) and part B the resulting success fraction.

**Detection.** After the first useful collision the two walkers are the same class
sequence up to an offset of at most one checkpoint period plus one traversal (<= 128 steps,
section 3.1), and from their first common class the wait to a distinguished class is that of
one trajectory: within 20/theta steps except with probability e^-20 (H1). The window is applied to a trail selected through a collision, not to a fresh trail:
under H1 the distinguished indicator of each step is independent of the past with
probability theta, so the residual wait of the selected trail is geometric with the same
parameter and its tail is e^-20 too. No length-biased trail allowance is claimed. Charging all M walks in lockstep over that window costs E = M (20/theta + 128) =
2^84.33 + 2^27 steps: the 128 is the offset above, charged although the matcher is offline;
it is 2^-98.9 of S.
The losses, each a probability bound: detection tail e^-20 < 2^-28; trailing walk
restarted before reaching the distinguished class (its distinguished-free run must
exceed 40/theta: under H1 the run containing a given step exceeds 40/theta with
probability (1 + 40) e^-40 < 2^-52); self-collisions of a single walk (fraction
<= 2^-19 of collisions, with M = 2^20 walks); a false-key pair with B1 = B2 (1/r); key ties (< 2^-249);
long fruitless cycles (< 2^-100, now only all-multiply and at the mixed rate of section
3.2, so the bound only improves); more than 2^30 stuck-walk restarts (< 2^-21);
more false-key pairs than the hard pair budget (< 2^-36 by Markov, section 3); a slot
that collides with itself early and wastes its budget (<= 2^-19.91 over the run, section
4(ii)), whose success loss is <= 2^-21.62; a slot aborted by its meter, i.e. productive
work above (1 + 10^-5) of nominal (< 2^-100) or cycle work above its allowance plus the
cycle band (<= 2^-7.07, 0.00745, section 4(v)); a walk slot stopped by the per-slot DP-attempt cap,
i.e. a healthy slot reaching R_w = 2^47 attempts (the KL figure of section 3, < 2^-300);
two distinct classes with equal T0 causing a spurious escape (< 2^-253). Within-walk
returns are charged twice over: an escaped return sits in the cycle overhead of section
3.2, an unescaped one is the wasted slot above, and a return with a different exponent
pair is the useful collision. Their sum is 3.399 x 10^-6 < 2^-18; with the 0.00745 charge
of (v) the ledger consumes 0.0074534 of the 0.01 allowance and leaves 0.00254.

**The work meter, accounted per step type.** A step costs c_sq = 0.50490 or c_mul = 1.1125
units, so the walk phase is metered and a slot stops at W_cap (section 5). The 1% band is
split once between productive and cycle work (4(v)), and Hoeffding is applied to the
productive steps alone: the branch bit is a function of the current class, a revisit repeats
an earlier decision, and step costs are not i.i.d. across time. The slot's work is
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
2^105.955 steps meets its own class again at all is at most n^2/(2N') =
2^211.910/2^252.785 = 2^-40.875 per slot (the union over its step pairs), and a
self-collision within its first O(1) steps has probability O(1/N') = 2^-251.785,
negligible. It wastes at most its remaining n steps, and its meter,
DP-attempt cap or a restart ends it. It is charged on both counts:
  (a) **Work.** The wasted steps are the slot's own step budget, already inside S and
  inside its meter, and the meter stops the slot at W_cap, so the hard walk bound
  W_walk <= M x W_cap of section 5 holds unchanged. As a figure: at most M n^2/(2N') = 2^-20.875 wasted slots in
  expectation, i.e. 2^-40.875 S = 2^85.08 steps, against the run's M x 2^80 = 2^100-step
  revisit allowance, and at most 2^-20.96 S except with probability 2^-19.91 (Markov).
  (b) **Success.** A slot that spends its budget on repeats exposes no fresh classes, so
  the productive count is T >= S - L with E[L] <= 2^-40.875 S, and L < 2^-20.96 S except
  with probability <= 2^-19.91 (Markov). A shortfall of that fraction lowers the collision
  parameter by (1 - 2^-20.96)^2 and the success probability by at most 2^-21.62, the
  entry below; the 2^-19 self-collision loss already reserved is kept, since a wasted slot
  also loses the pairs it would have formed. Per slot the loss is at most its own n steps
  of fresh classes, which is what the two charges bound.

(iii) **Cycle handling is escaped work, charged as a declared allowance.** Fruitless
2-cycles and checkpoint cycles are detected by the T0 comparisons, which run before the
report append, and escaped (section 3.2). Under H1 the expected cycle-handling work is
3 (1-s)^2/(12R) + 200 x 0.0101 R^-2 = 5.206 x 10^-5 steps per step, and the package
charges OVER = 3/(12R) + 200/R^2 = 4.3488 x 10^-4 steps per step instead, i.e. 1.000435 x S,
which is 8.35 times that mean. This is a **deterministic allowance in the ledger, not a tail
bound**: no Chernoff or weighted-MGF argument is claimed for the count of cycle events, whose
steps are not i.i.d. (history-dependent: a cycle's classes and the walker's checkpoint phase
decide it). What covers those steps in work is the hard meter of section 5, which charges
each of them at its step price and stops the slot at W_cap. The allowance's job is only to
keep the nominal from being an under-charge; the toy of part B measures the mixed rates it is
built on. Its 8.35x factor and the 2^91.73-step expected per-slot cycle work are figures of
H1, listed as such in section 6; (v) turns that mean into a Markov term at the meter's band.

(iv) **The productive sum.** With n = N_slot = 2^105.96 steps per slot and
c_mul - c_sq = 0.6076, Hoeffding at the excess a = 10^-5 gives P(productive work >
(1 + a) c n) <= exp(-2 (a c)^2 n / (c_mul - c_sq)^2), exponent 3.298 x 10^-10 nats per
step: the probability is 2^(-2^74.98) for a slot and 2^(-2^94.98) over the run's S steps,
under (i)'s fresh-class law, i.e. H1, not unconditionally. The 74.98 and 94.98 are
exponents of the exponent, not probabilities; part C prints both.

(v) **The abort ledger: the 1% band is split once.** The cap is W_cap = 1.01 nominal +
2^80 steps; one meter sums productive, cycle and report work (report work is inside
nominal, bounded by the DP-attempt cap). A slot aborts only if (i) its productive work
exceeds (1 + 10^-5) of its productive nominal, 2^(-2^74.98) per slot under H1, the union
over the M = 2^20 slots below 2^-100; or (ii) its cycle work exceeds its allowance plus the
cycle band, 0.01 nominal + 2^80 steps less 10^-5 of the productive nominal (0.00999 of
nominal) = 2^98.80 steps per slot (part C, `cycle_abort`). With neither true the meter
stays below W_cap. H1's expected cycle work is 2^91.73 steps per slot, so Markov gives
P(cycle work > band) <= 2^91.73/2^98.80 = 2^-7.07 (0.0074251), charged as 0.00745. A slot
stopped this way is charged as if it stored no DP (DPs stored before the stop stay in the
log), so its success loss is at most that fraction of the walk population, inside the 0.01
allowance. Event (i)'s slot loses at most its remaining productive steps; a self-colliding
slot's repeat work is the wasted-slot case of (ii), charged at 2^-19.91 over the run.

The 0.01 allowance, line by line: 0.01000 total, less 0.00745 for event (ii), less
0.0000034 for the other named losses (3.399 x 10^-6), leaves 0.00254 (rounded down) as H1
slack; success >= 0.42013 - 0.01 - 2^-18 = 0.41012 >= 0.39, margin 0.02012.

Success >= 0.42013 - 0.01 - 2^-18 > 0.42013 - 0.01 - 0.00001 = 0.41012 >= 0.39 (every
figure rounded down), with 0.02012 of margin. The 0.01
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
| common: steps-since-attempt counter and stuck-walk test | 7 |
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
  its lifetime, so there is no synchronisation, atomic or cache-line contention on it; the
  log's per-slot regions leave the walk phase no shared counter.
* Values. W_cap = ceil(1.01 x nominal) + 2^80 x I_MUL, with nominal = c x C x (S_slot +
  E_slot) + I_MUL x ceil(OVER x S_slot) + I_REP x R_w =
  1227616758516808452599261592668764720 = 2^119.885272 instructions: the third term is the
  worst-case DP-attempt charge, 2^47 attempts at 3,000 instructions, 2^-61.33 of the
  nominal; the 2^80-step revisit allowance adds 2^94.45 instructions, so
  W_cap = 1239892953000576023550753345807868368 = 2^119.899627. The meter's largest value
  is W_cap - 1 + I_MUL + I_REP = 1239892953000576023550753345807893617 < 2^127, so it fits
  two 64-bit words (a step that is also a DP attempt adds both charges). A slot is stopped
  when its meter reaches W_cap, so W_walk <= M x (W_cap - 1 + I_MUL + I_REP) =
  2^139.899627 instructions = 2^125.611915 units is a **hard bound**, not an expectation;
  it contains every step and every DP attempt, and it is what the claim is computed from
  (the attempt term raises the bound by 2^-108.34 relative, part C checks the integers).
* Report attempts are inside that meter. I_REP = 3,000 instructions are added at every
  distinguished hit, stored or suppressed, so a fruitless cycle that keeps revisiting one
  distinguished class pays per pass inside the same hard cap, whether or not the no-report
  restart counter is reset; at most R_w = 2^47 attempts are possible per slot, and even that
  worst case is 2^-54.66 of the 1% band the cap carries above the nominal, so report work
  cannot consume the band. The same work is charged again in the capacity group below (3,000
  of the 13,100 per cap unit), a worst case over the log whether or not the meter stops the
  slot first; charging it twice is deliberate.
* Abort bound. Two events, the 1% band split once (section 4(v)). (i) Productive work above
  (1 + 10^-5) of nominal: 2^(-2^74.98) per slot, 2^(-2^94.98) over the run (exponents of the
  exponent, under H1), below 2^-100 over the M slots. (ii) Cycle work above its allowance
  plus the cycle band, 0.01 nominal + 2^80 steps less 10^-5 productive nominal, is
  2^98.80 steps per slot against H1's 2^91.73-step expected cycle work, so Markov gives at
  most 2^91.73/2^98.80 = 2^-7.07, charged as 0.00745 inside the 0.01 allowance, which
  leaves 0.00254. A slot that collides with itself cannot pass the meter either: it stops
  at W_cap with fewer productive steps, so (i) is a work bound whatever the walk does; the
  2^-40.875-per-slot event costs fresh classes, charged in section 4(ii)(b). The 0.01
  allowance covers H1's heuristic parts and these two fractions.

**Counters and indices.** Every counter and index, its largest value at the end of the
attack (including one over-cap update per processor), and its width. Part C recomputes
this table and asserts every row.

| counter or index | largest value (log2) | words | capacity test |
| --- | ---: | ---: | --- |
| slot work meter (instructions) | 119.899627 | 2 x 64 | >= W_cap with W_cap = 2^119.899627 < 2^127 |
| slot step counter | 106.597846 | 2 x 64 | none (informational) |
| steps since last DP attempt | 65.321929 | 2 x 64 | > 40/theta (2^65.321929 rounded up), a 128-bit compare |
| slot DP-attempt counter r | 47 | 1 x 64 | >= 2^47 with 2^47 < 2^63 |
| log offset 80 x (2^47 slot + r) | 73.321929 | 2 x 64 | < 2^128 |
| slot restart counter | 41.275917 | 1 x 64 | < 2^63 |
| matcher pair counter | 41 | 1 x 64 | > 2^41 with 2^41 + 1 < 2^63 |
| heapsort index | 67 | 2 x 64 | < 2^128 |

No counter's test sits at its word's own capacity: the two that can fire are at
2^119.899627 of a 2^128 word and 2^47 of a 2^64 word.

**Fruitless-cycle overhead factor (a declared allowance, not a tail bound).** The step meter
also carries the cycle-handling steps: the nominal charged work adds
OVER = 3/(12R) + 200/R^2 = 4.3488 x 10^-4 extra steps per step at the multiply-step price
(2-cycles at rate 1/(12R), <= 3 steps each; long cycles at rate R^-2, <= 200 steps each), i.e.
it charges 1.000435 x S in step count, which appears in the closed form below as 1.0006201
times the plain step work. The H1 mean is 3 (1-s)^2/(12R) + 200 x 0.0101 R^-2 = 5.206 x 10^-5
per step, so the charge is 8.35 times it. This is charged work under the same hard meter;
section 4(v) turns its H1 mean into a Markov term, 2^-7.07 charged as 0.00745, and section
3.2 bounds the all-multiply long-cycle rate it is built on, which part B measures mixed.

**Setup, output and preprocessing.** Setup and output are all <= 2^31 units together:
pairing g (bounded by 2^14), input check z^r = 1 (<= 510), table: 1024 multipliers x (2
exponentiations <= 1,021 + 11 Frobenius maps) <= 2^20.02, M = 2^20 walk starts <= 2^30.0,
output [x]P0 (<= 510 G1 operations) and verification pairing (<= 2^14). **Preprocessing and
table construction are at most 2^57 group-operation equivalents**: the table, the walk
starts and the initialization of all C_rep = 2^67 log records with the sentinel record,
2^67 x 12 = 2^56.30 units. Every term is charged inside W (the fill in the report group
below, the rest in its 2^31 entry); nothing is amortised and no advice is used.

**Restart starts**, charged at the worst case the meter allows: a slot's step counter can
reach 2^106.597846, and every restart needs 40/theta = 2^65.32 steps, so a slot has at
most 2^41.275917 restarts; each costs one exponentiation (1,021 units, 2^24.29
instructions) plus 100, so M x 2^41.275917 x (1,021 x 20,000 + 100) = 2^71.271692 units.

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
units; the sentinel fill is 2^67 x 12 = 2^56.3 units. The restart-starts term (2^71.271692)
dominates the group and setup/output (2^31) is far below it: together **2^71.320 units**,
charged in full whether or not the reports occur (the capacity is deterministic).

W <= 2^125.611915 + 2^71.320 units = 2^125.611915 (every figure rounded up: the ancillary
term raises log2 W by 2^-54.29/ln 2 < 2^-53.8, below the last printed digit). Closed form,
each factor rounded up: W <= 1.01 x 0.78021875 x 2^125.954715 x 1.000621 + 2^94.45 (meter
cap x average step cost x step budget x cycle overhead, plus the revisit allowance); the log2
of the four factors is 0.014356 - 0.358049 + 125.954715 + 0.000895 = 125.611917 (the sum of
the rounded-up terms, so above the exact 125.611915 that part C prints), and the revisit
allowance raises log2 W by 3.13e-08.

Submitted bound: **time_log2 = 125.612** = 125612/1000 >= 125.611915 >= log2 W (part C
checks it as an exact integer inequality, W^1000 <= C^1000 x 2^125612, and that 125611
fails, so 125.612 is exactly the ceiling of log2 W at three decimals).

Memory (a hard bound, every live array counted): the report log of C_rep = 2^67 records
of 80 bytes is 2^67 x 80 < 2^73.3220 bytes, allocated in full and never exceeded; the
heapsort works in place with one 80-byte temporary; the multiplier table is 12,288 x 576
bytes = 2^22.76 plus its exponent tables (< 2^20); the walk states are 2^20 x < 1 KiB =
2^30 bytes. Total < 2^73.4 bytes; no hash table, chaining or eviction, so no load factor arises.
Memory is reported, not scored.

## 6. Heuristics and limitations

* **H1 (random-mapping behaviour of the mixed walk).** The class sequence is treated as a
  random mapping on N classes in which only a (1-s) fraction of the steps can merge two walks,
  with the Brent-Pollard correction: the standard heuristic behind rho estimates (Teske
  1998; Brent-Pollard; van Oorschot-Wiener 1999; Bos-Kleinjung-Lenstra 2010), with the
  branch-adjusted law of section 4.1. Part E measures the
  collision law and part B the success fraction; it is not proved.
* **The decision bits and the report key.** The branch probability is exactly 35/64 if the
  six branch bits of T0 are uniform, and the partition and distinguished fields are uniform
  if the other bits are; T0 = Re(c_0) is not a uniform Fp element, so these are heuristic. The
  fields are read from the Montgomery value T0M (section 2), the bijective image of T0, so the
  assumption is the same in either representation. Support: part A2 measures the frequency of
  each of the 76 used bits on real GT elements, part P runs the production rule on real GT
  elements, and part B's toy walks measure the fractions against their nominal values.
* **The report key's bits are taken as pseudo-random under H1, not proved uniform.** K is the
  minimum of the twelve Montgomery keys, so it is not a uniform Fp element, and x -> x R mod p
  does not make selected bits of it uniform: a bijection preserves a uniform distribution on
  the whole field, not the selected bits of a nonuniform one. The package claims no
  uniform-bits preservation. It takes bits 70..197 of the Montgomery-ordered minimum K as
  pseudo-random under H1, the same random-mapping heuristic the rest of the ledger rests on,
  with part A2 as evidence at its scope: 2,000 real GT elements give per-bit largest deviation
  4.83 and mean 1.12 standard errors, which is what 128 such bits give and is not a proof of
  uniformity. Those bits price the report-cap margin (section 3) and the 1/r false-key term
  (section 4), both carried as heuristics.
* **Coalescence only up to the checkpoint period.** The ordinary step is a function of
  the class, but the two cycle escapes read stored history. Section 3.1 proves what survives:
  escape class and action are functions of the cycle, and two walkers on one cycle follow
  the same class sequence up to a time offset of at most 64 + L <= 128 steps, so section 4.1
  charges M (20/theta + 128) steps for the detection window (2^-98.9 of S) rather than
  claiming a one-cycle coalescence. The recovery needs none of it: two records of one class
  with different exponent pairs already solve for x (section 3.1(iii)). Part P runs both
  escapes on constructed cycles (a 2-cycle by the ring, a 3-cycle by the checkpoint with the
  full traverse-and-square escape) and compares the walkers' class sequences afterwards.
* **Every distinguished hit is a charged attempt.** The counter is capped at R_w per slot and
  each attempt costs 3,000 instructions inside the slot's meter and again inside the 2^67-cap
  capacity group, so a repeated report inside an unescaped fruitless cycle pays per pass and
  neither the restart counter nor duplicate suppression can move it out of either bound.

* **Preprocessing and advice.** At most 2^57 group-operation equivalents, including the
  2^56.30-unit initialization of all 2^67 log records, the multiplier table (<= 2^20.02) and
  the walk starts and setup/output (<= 2^31), all inside W (section 5), so nothing is
  amortised; nonuniform advice is 0 bytes.

* **The work meter's dependent-process accounting.** Step costs are not i.i.d. across time (the
  branch bit is a function of the current class). Section 4 splits the slot's
  work into productive steps (i.i.d. costs on newly exposed classes, where H1's
  random-mapping law does supply a fresh draw, bounded by Hoeffding), cycle handling (escaped
  work, charged at 8.35 times its H1 mean rate as a declared allowance; no concentration
  inequality is claimed for the count, but the mean gives the Markov abort term of 4(v)) and
  the self-collision event, which
  spends a slot's remaining budget on repeats with probability at most n^2/(2N') = 2^-40.875
  per slot and is charged as a 2^-19.91 lost-slot and 2^-21.62 success-loss entry. The meter's
  bound on total work is hard (section 5) and holds under that event too, and of section 4(v)'s two
  abort events, (i) is below 2^-100 over the run and (ii) is charged at 2^-7.07 per slot. What remains heuristic - the fresh-class law,
  the cycle count and the report-key bits - is supported by part A2 and part B. No "S/1.01
  steps" equivalence is used.
* **No production run.** The toys are evidence only at their own scale: distinguished
  fraction 2^-4 or 2^-3 rather than 2^-60, field near 2^20 rather than 2^381, class count
  2^21.4 to 2^22.4 rather than 2^251.3, T0 20 bits rather than 381. The production cycle
  rate is the derivation of section 3.2, not a toy row.
* **The squaring price is a measured instruction count**, not a model: the Callgrind count of
  blst f7a88067 `blst_fp12_cyclotomic_sqr` (9,798 per call, plus 300) on the toolchain the
  cost model names, calibrated against the cost model's own 383.0 and 22,996.0 and re-derived
  by part C from the shipped record. Its straight-line body makes the count
  input-independent; the residual risk is different code generation on the grader's compiler,
  which is why the charge keeps 147 instructions of slack.

* This is a generic attack (about the square root of the class count). It ignores
  index-calculus methods in Fp12 (TNFS family); it is a valid upper bound, not a
  security estimate.

## 7. Experiments (`experiment.py`, standard library only)

`experiments.md` gives the per-part detail and `measure/` the measurement record. In summary:
Part A checks on real BLS12-381 GT elements the identities of section 2, the twelve-way argmin
and class key in the Montgomery domain, and cyclotomic squaring = y*y; Part A2 measures the
frequencies of the 76 used bits of T0M and of the 128 report-slice bits 70..197 of the
Montgomery-ordered K; Part P runs the production step rule and both cycle
escapes on real GT elements; Part E measures the branch-adjusted collision law in a toy field;
Part B runs the full mixed walk, attempt cap, meter and matcher in three toy fields near 2^20;
Part C recomputes the integer bound arithmetic of sections 3-5, including the exact claim
check and the squaring price re-derived from `measure/`, found from the source file and printed
absolutely, so the entrypoint runs anywhere.
