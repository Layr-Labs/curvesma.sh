# Frobenius-class parallel Pollard rho in GT (success target 0.39)

Claim: time_log2 = **125.826** (total work W <= 2^125.826 dense-Fp12-multiplication
equivalents, C = 20,000), success probability >= 0.39, peak memory < 2^73.4 bytes
(a hard bound, section 5), preprocessing <= 2^31 units (included), no advice.
Heuristics are listed in section 6.

Design provenance: the algorithm, the four-multiplication key formula and the twisted
table are those of public submission
3da91365-2d58-40e5-ac7d-415e3061f7d8 (Frobenius-class rho, 125.89). This package keeps
that algorithm unchanged and differs only in: (1) the success target is the track's
0.39 instead of 0.40 (same 3% step slack, same margin structure); (2) a written
derivation of the fruitless-cycle rates (section 3.2); (3) an explicit stuck-walk rule
and record format; (4) toy experiments that run the production partition and
distinguished-point predicates on the raw class key, the production step order and append
rule, and the production offline matcher (full-budget walk, full-80-byte-record sort,
adjacent dedup, charged equal-key pairs under a hard pair budget, exponentiation check); (5) a directly
addressed distinguished-report log whose capacity is the deterministic product
M R_w = 2^67 records, with a sorted matching phase whose per-record work and whose only
quadratic term (the equal-key pairs) are both hard-bounded, so its memory and processing are
worst-case bounds with no report-count tail claim (section 3, section 5); (6) the per-step non-reference work is
itemized instruction by instruction, with the class-key minimum taken from the sign
structure of the twelve keys instead of a 12-way scan (section 5). No earlier result is used as
evidence; every bound below is derived or computed here.

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

* Frobenius: (sum c_m w^m)^p = sum conj(c_m) d_m w^m with d_m = xi^{m(p-1)/6}.
* p has multiplicative order exactly 12 modulo r (embedding degree 12), so the
  automorphism group {y -> y^{p^k}} of GT has order 12; it contains inversion
  (p^6 = -1 mod r). Every y != 1 lies in a class [y] of exactly 12 elements.
  Number of classes N = (r-1)/12, log2 N = 251.2721.
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

## 3. Walk (R-adding walk on classes, Duursma-Gaudry-Morain / Wiener-Zuccherato style)

Parameters: R = 1024 partitions, distinguished-point fraction theta = 2^-60,
M = 2^20 parallel walks, each started at g^alpha z^beta with fresh uniform alpha, beta mod r.

Precompute M_h = g^{a_h} z^{b_h} (h < R, random a_h, b_h) and the table
T[h][k] = M_h^{p^k}, k = 0..11 (12,288 entries of 576 bytes, 7.1 MB), with exponent
tables a_h p^k, b_h p^k mod r.

Class map: F([y]) = [ y^ * M_h ], h = K([y]) mod R. A walk keeps an arbitrary
member y of the current class with its exponents (alpha, beta), y = g^alpha z^beta.
Because y^{p^kappa} * M_h = (y * M_h^{p^{-kappa}})^{p^kappa}, one step is

```
compute keys, kappa, K                      (4 Fp mul + 334 other instructions, section 5)
if the class two steps back is this class (2-cycle):
    escape: square the smaller-K class of the pair, drop both key histories, continue
if the checkpoint class (recorded every 64 steps) is this class:
    escape: traverse the cycle once to its smallest-K class, square it, continue
if bits 10..69 of K are all zero (distinguished):
    if (bits 70..197 of K, p^kappa*alpha, p^kappa*beta) differs from this walk's last record:
        append that 80-byte record to the log
        increment this walk slot's cumulative report counter r
        if r = R_w = 2^47: stop this walk and leave its slot idle (rule below)
h = K mod R (bits 0..9) ; t = (12 - kappa) mod 12
y = y * T[h][t]                             (ONE dense Fp12 multiplication)
alpha += a_h p^t ; beta += b_h p^t (mod r)  (2 additions mod r)
```

So each step costs exactly one reference operation plus the key work. The walk is
a deterministic function of the class, so once two walks reach the same class they
follow the same class sequence and reach the same distinguished class.

**Report log (deterministic capacity, direct addressing, per-walk cap).** Reports are
appended to a log: an array of exactly C_rep = M R_w = 2^20 x 2^47 = 2^67 record slots
of 80 bytes each, allocated before the walks start and addressed directly by a report
counter c (a 128-bit counter, two 64-bit words, because 2^67 does not fit in one): the
record goes to address base + 80 c. An append is one 80-byte compare, one counter compare,
one 80-byte store and one increment, a fixed instruction count with no hashing, probing,
bucket or load factor.

*Order of the step, and the append rule.* The 2-cycle check and the checkpoint-cycle check
run before the distinguished test, so a walk that comes back to a distinguished class
through a fruitless 2-cycle has already escaped at the return step and appends nothing
there. A report is stored only when it differs from the last record this walk appended
(the same 128-bit key and the same two canonical exponents); an identical record is
dropped. The identical records a fruitless 2-cycle would otherwise append -- at rate
1/(12R) per step, about S/(12R) x 2^-60 = 2^52.1 events over the run -- therefore never
reach the log: the check fires at the return step, before the test. The append rule costs
one 80-byte compare against the walk's own last record, inside the 1,000-instruction append
charge of section 5; the reordering is free, because the two checks are already two of the
step's counted compares. Neither rule changes a step: a walk whose record was dropped keeps
walking. They are narrow on purpose: a later report with the same key but different
exponents is the useful collision this algorithm exists to find, and both rules store it.

Two observations keep the account exact. First, a class determines the canonical
ELEMENT, not the exponent pair: the stored pair (A, B) is a (g, z) representation of that
element, y = g^alpha z^beta, and a walk that revisits a distinguished class normally
arrives with a DIFFERENT pair, because the accumulated exponents depend on the path taken.
Two such records are the useful collision this algorithm exists to find: the record holds
the element g^A z^B = g^{A + B x}, so equality of two distinct records gives
A1 - A2 + (B1 - B2) x = 0 mod r, i.e. x = (A1 - A2)/(B2 - B1) mod r, which is exactly the
recovery the matcher below performs. The matcher skips a pair with B1 = B2, which for two
such distinct records cannot happen unless their elements differ, so it never relies on
the pairs coinciding.

Only a FRUITLESS revisit -- the walk returns to the same element with the same exponents,
which is what a multiplier cycle is (section 3.2: the accumulated multipliers are the
identity as formal expressions) -- reproduces an identical record, and the two cycle checks
fire before the append on exactly those returns.

Second, the two rules remove the fruitless volume but not every exact duplicate: a cycle of
length 3..64 is escaped within 64 + L steps yet can report a distinguished class twice
before that, and a cycle longer than 64 (probability <= 2^-100 per step, section 3.2) is
never seen by the checkpoint, so its slot can append an exact duplicate of a record it
already stored before the per-walk cap stops it. Every such record equals an earlier record
of the same walk in key and in both exponents, so it costs one log slot, and the full-record
heapsort below puts it adjacent to its mate (identical records compare equal in all 80
bytes); the adjacent dedup drops it as the same record, one 80-byte comparison, before
anything is paired. It is not a pair and it is not recomputed. A short-cycle repeat is at
most theta x 193 per event, i.e. <= 2^50 over the run; a long-cycle slot stops at
R_w = 2^47 and its steps are charged to its walk. Neither changes the log's capacity,
which is the hard count above and does not depend on their number.

*Per-walk report cap (the capacity rule).* Each of the M walk slots keeps a cumulative
report counter r that no restart resets. When a slot's counter reaches R_w = 2^47 the
walk is stopped: its trail is discarded, the reports already stored stay, and the slot
takes no further steps (a slot that has appended 2^47 reports is cycling on a class
cycle and cannot reach a new class, so restarting it could only produce more identical
reports). Its steps count against the global step budget, so the waste is charged, not
free. Therefore the log holds at most C_rep = M R_w = 2^67 records. This is a
deterministic count; no probability statement about the report count is used anywhere
in the memory or processing bound.

The counter belongs to the slot state and is initialised once, when the slot starts. A
stuck-walk restart (rule below) replaces the walk's trajectory, its class history, its
checkpoint and its last-record buffer, and does not touch the slot's report counter, so a
slot's lifetime appends are bounded by R_w however many times its walk restarts.

*Why the cap is safe for healthy walks.* A healthy walk's report count over its
S/M = 2^105.671 steps is a sum of independent Bernoulli(theta = 2^-60) indicators under
H1, with mean 2^45.672 and relative standard deviation 2^-22.8. R_w = 2^47 is a factor
a = 2.511 above that mean, i.e. 2^23.4 standard deviations above it, and the Bernoulli
tail exponent D(a||1) = a ln a - a + 1 = 0.8008 gives P(a slot reaches the cap)
<= exp(-0.8008 x 2^45.672) = 2^-(2^45.88), far below 2^-300 (the union over 2^20 slots
is still far below 2^-300). A walk stopped by this rule wastes at most the steps it had
taken, and those steps are already inside the step cap of section 5, so the loss is a
step-budget displacement, not extra work: the run makes the same N_steps = 1.000435 S + E
= 2^125.672 steps either way.

Record format: a report is stored as 80 bytes = bits 70..197 of K (16 bytes, the key
T: the 128 bits directly above the distinguished bits, since bits 10..69 of every
reported K are zero by construction) plus the canonical exponents A = p^kappa alpha,
B = p^kappa beta mod r (32 bytes each, reduced).

**Matching by an in-place heapsort of the log (worst-case bound).** The walk phase
ends after its step budget (section 5: S steps plus the lockstep window E). Then the
n <= C_rep = 2^67 records of the log are sorted in place by the WHOLE 80-byte record
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

Collision: the first charged pair whose canonical elements are equal with B1 != B2
(probability 1 - 1/r per true cross-walk collision, for random starting exponents) is
solved as p^k1(alpha1 + beta1 x) = p^k2(alpha2 + beta2 x) mod r, that is
x = (A1 - A2)/(B2 - B1) mod r (one modular inversion, <= 2^14 instructions). Verify
g^x = z (<= 510 units) before output.

**Fruitless cycles (rule).** Each step compares K with the key two steps earlier;
on equality the two classes {C1, C2} form a 2-cycle and the walk squares the stored
representative of the one with smaller K (y -> y^2 maps a class to a class:
(y^{p^k})^2 = (y^2)^{p^k}), doubling its exponents mod r (alpha, beta <- 2 alpha, 2 beta
mod r). A checkpoint key is stored every
64 steps and compared at every step; a repeat means a cycle of length L <= 64, which
is traversed once more (L steps) to find its minimum-K class, which is squared. Both
escapes depend only on the set of classes in the cycle, not on where or when the walk
entered it, so two walks that have merged leave the cycle at the same class and stay
merged. Detection of an L-cycle (L <= 64) happens within 64 + L steps of entering it,
so the waste per event is <= 64 + 2L + 1 <= 193 <= 200 steps; a 2-cycle event wastes
<= 3 steps.

**Stuck-walk rule.** Each walk counts steps since its last distinguished class in a
128-bit counter (reset at each report and at restart; this is the trajectory counter, not
the slot's report counter r, which no restart resets; 40/theta = 2^65.33 exceeds 2^64). When the count exceeds 40/theta = 2^65.33 the walk is
abandoned and restarted at a fresh g^alpha z^beta. Its steps are real steps and are
inside the walk-step cap of section 5. It is only a safety net for cycles longer than
64 (section 3.2 bounds their total probability by 2^-100). Under H1 the steps between
distinguished classes are geometric with parameter theta, so a run exceeds 40/theta
with probability (1 - theta)^(40/theta) < e^-40; over at most 2^66 runs the expected
number of restarts is < 2^66 e^-40 < 2^8.3, and by Markov's inequality more than
2^30 restarts occur with probability < 2^-21. 2^30 restarts waste at most
2^30 x 2^65.33 = 2^95.33 steps, a fraction < 2^-30 of S.

**Word widths, no overflow.** The exponents alpha, beta of every walk, the table
exponents a_h p^t, b_h p^t and the stored A, B are always held reduced in [0, r),
r < 2^255, as four 64-bit limbs. Each update is followed by a reduction mod r: the sum
of two reduced values is < 2r < 2^256 and fits in four limbs; an escape doubling is such
a sum of a value with itself, so any number of doublings leaves the exponent reduced
and it never accumulates unreduced bits; the product p^kappa x alpha of a report is a
512-bit product reduced mod r (inside the 1,000-instruction append charge of
section 5). Because g and z have order r, the reduced pair is the exponent of y mod r
at every step. The counters are 128-bit (two limbs, add and add-with-carry): steps
since the last report (< 2^66), steps per walk (< 2^107), the global step count
(< 2^126), the report counter c (<= 2^67) and the record address 80 c (< 2^73.33).
Three of them change at every step (steps since the last report, steps of this walk, global
step count): each is one load pair, ADDS, ADC and one store pair (4 instructions),
12 per step; the stuck-walk test and the global step-cap test are two 128-bit compares
with a constant (CMP, CCMP, B.cond), 6 instructions. Both are lines of the section 5
ledger. The report counter c, the record address, the walk slot's cumulative report
counter r (<= 2^47, one 64-bit word in the walk state) and the reset of the
steps-since-report counter change only on a report, inside the 1,000-instruction append
charge of section 5.

### 3.2 Fruitless-cycle rates (derivation)

Let a walk hold y_0 and take L steps with indices h_1..h_L and twists t_1..t_L. Then
y_L = y_0 * prod_i M_{h_i}^{p^{t_i}}. The walk is back in the class of y_0 with the
same exponent pair (a fruitless cycle, i.e. one that occurs for every x, not by the
1/r chance of a random mapping) iff y_L = y_0^{p^e} as formal expressions in the
independent uniform exponents (alpha_0, beta_0, a_h, b_h). Comparing coefficients
forces p^e = 1 and, for every index value h that occurs,
sum_{i: h_i = h} p^{t_i} = 0 mod r. Two consequences:

* A single occurrence is impossible (p^t != 0), so every index that occurs appears
  >= 2 times; the L steps use d <= floor(L/2) distinct indices.
* For a block of k >= 2 occurrences, sum p^{t_i} = 0 mod r. p mod r is a primitive
  12th root of unity, whose minimal polynomial is Phi_12, so this is the event that a
  sum of k 12th roots of unity vanishes in Z[x]/Phi_12 (for k < r^{1/4} a nonzero sum
  has norm < r, so "= 0 mod r" equals "= 0"). With uniform twists this is the return
  probability of a symmetric random walk; it is 1/12 at k = 2 and, by the identity
  P(S_{2m} = 0) = ||mu^{*m}||_2^2 (non-increasing in m) and Cauchy-Schwarz for odd k,
  at most 1/12 for every k >= 2 (part C computes k = 2..10 exactly: max 1/12).

Under H1 (indices and twists of successive distinct classes behave as independent
uniform values, the same assumption that gives the 2-cycle rate 1/(12R)), the
probability that the walk enters a fruitless L-cycle at a given step is

  P_L <= sum_{d=1}^{min(R, floor(L/2))} C(R, d) (d/R)^L 12^-d,

Here D is the number of distinct indices among the L steps; the cycle needs D = d
with every one of the d blocks vanishing, so P_L <= sum_d P(D = d) 12^-d, and
P(D = d) <= P(all L indices lie in some d-set) <= C(R, d) (d/R)^L (union bound).
Computed exactly in part C of experiment.py for R = 1024 (part C also prints
R = 64 and 256: 0.1422 R^-2 and 0.1397 R^-2, used for the toy comparison in section 7):

* sum_{L=3..64} P_L = 0.139 R^-2 <= R^-2 (L = 3 alone is R^-2/12, L = 4 is 0.056 R^-2).
* sum_{L=65..2R} P_L <= 2^-237.3. For L > 2R every term with d < R shrinks by d/R
  per extra step, so sum_{L>2R} <= R * P_{2R} + 2^127 * 12^-1024 < 2^-226.

Hence the long-cycle checkpoint fires at rate <= R^-2 per step, and over all
<= 2^126 steps a fruitless cycle longer than 64 occurs with probability < 2^-100.
Squaring steps cannot be part of a fruitless cycle (2^t is not in <p> mod r for small
t). This derivation -- the union bound over index sets and the Phi_12 zero-sum
probability above -- is the production cycle-rate bound; the cost in section 5 charges
the round R^-2 with 200 steps per event. The toy measurements of section 7 are a sanity
check only where the toy's regime holds, R << S0. At R = 64 (R/S0 = 2^-4.7) and R = 256
(R/S0 = 2^-2.7) the measured rates are 0.063 and 0.55 R^-2 (51 and 32 events), inside the
charged R^-2 but 0.44x and 3.9x the sharp 0.139 R^-2 derived here, in opposite directions
and outside Poisson noise. At R = 1024
(R/S0 = 2^-0.7) the toy is outside the regime: R = 2^10 is comparable to its S0 = 2^10.7
steps, partition indices must repeat inside one short walk and the independence premise of
this derivation does not hold, so that row (13 R^-2) neither supports nor contradicts the
bound and is not used for it. In production R/S0 = 2^-115.6. Section 7 gives the counts
and states the scope.

## 4. Success probability

Heuristic (H1): the class sequence of the walk behaves like a random mapping on
N classes, with the Brent-Pollard correction for R-adding walks (effective
N' = N * R/(R-1)). Then the probability that no useful class collision occurs among
the first S distinct-step points of all walks is exp(-S^2/(2N')). Choose

S0 = sqrt(2 N' ln(1/0.61)) = 2^125.6285,   S = 1.03 * S0 = 2^125.6711,

where the 3% is explicit slack for walk non-randomness; at S0 the model gives
exactly 1 - 0.61 = 0.39, and at S it gives 1 - 0.61^(1.03^2) = 1 - 0.61^1.0609 =
0.4081. The algorithm stops after S walk steps (plus the detection window below), so
S is a fixed budget, not an expectation. Toy evidence (section 7) measures the
fraction of instances with a useful collision within S0 and within 1.03 S0 using the
production predicates.

Detection: after the first useful collision, the merged trajectory reaches a
distinguished class within 20/theta further steps except with probability e^-20.
Charging all M walks in lockstep over that window costs E = M*20/theta = 2^84.3
steps. The losses, each a probability bound: detection tail e^-20 < 2^-28; trailing
walk restarted before reaching the distinguished class (its distinguished-free run
must exceed 40/theta: under H1 the run containing a given step exceeds 40/theta with
probability (1 + 40) e^-40 < 2^-52); self-collisions of a single walk (fraction
<= 2^-19 of collisions, with M = 2^20 walks); B1 = B2 (1/r); key ties (< 2^-249);
long fruitless cycles (< 2^-100); more than 2^30 stuck-walk restarts (< 2^-21, section 3);
more than P = 2^41 charged equal-key pairs, i.e. more false-key pairs than the hard pair
budget (< 2^-36 by Markov, section 3); a cycle-handling
step count above its budget (< 2^-40, section 5); a walk slot stopped by the per-walk
report cap, i.e. a healthy slot reaching R_w = 2^47 reports (< 2^-300, section 3; the
cap itself is a deterministic bound, not a tail). Their sum is
< 2^-19 + 2^-21 + 2^-28 + 2^-36 + 2 x 2^-37 + 2^-300 < 2^-18 < 0.4081 - 0.39 = 0.0181.
Success >= 0.4081 - 2^-18 > 0.4081 - 0.0001 = 0.4080 >= 0.39.

## 5. Cost (units of one dense Fp12 multiplication = 20,000 instructions)

Per-step non-reference instructions, scalar AArch64, using blst's measured
blst_fp_mul = 383 instructions. A 384-bit (6-limb) operation is counted from its
actual instruction sequence: an add or subtraction mod p is 6 ADD/ADC (or SUBS/SBCS)
plus 6 SUBS/SBCS plus 6 CSEL, charged 20; a 6-limb compare with min, max and index
selection is 6 SUBS/SBCS plus 12 CSEL plus 2 conditional moves, charged 24; p - v is
6 loads plus 6 SUBS/SBCS, charged 20; a 4-limb addition mod r is 4 ADD/ADC plus
4 SUBS/SBCS plus 4 CSEL plus its load and store, charged 26; a 128-bit counter
increment is LDP, ADDS, ADC, STP (4) and a 128-bit compare with a constant is CMP,
CCMP, B.cond (3).

| Work | Instructions |
| --- | ---: |
| 4 Fp multiplications (a*zeta, two products for B, B*zeta) | 1,532 |
| 1 Fp addition mod p (B = a*Re(d_1) + b'*Im(d_1)) | 20 |
| 2 Fp subtractions mod p (a*zeta - a, B*zeta - B) | 40 |
| min3 and max3 of {a, a*zeta, a*zeta - a} and of {B, B*zeta, B*zeta - B} (compare, csel, index) | 144 |
| p - max_a, p - max_B (the negated minima of the two triples) | 40 |
| final minimum of the four candidates, kappa = 2j + {0,1} + 6*negated | 74 |
| load the w^1 coefficient (12 limbs) | 16 |
| h = K mod R, t = (12-kappa) mod 12, distinguished test on bits 10..69, table address | 20 |
| alpha, beta += table exponents mod r (2 x 4-limb add, conditional subtract) | 52 |
| 2-cycle and checkpoint key comparisons, key history stores | 40 |
| loop control, rotating output buffer | 30 |
| three 128-bit step counters (since last report, this walk, global), 3 x 4 | 12 |
| stuck-walk test and global step-cap test, 2 x 3 | 6 |
| **total counted (budgeted 2,250)** | **2,026** |

The budget of 2,250 instructions is 11% above the count and is the number charged
per step; by section 2 the twelve keys are two sign-symmetric triples, so the
canonical minimum needs the two triples' minima and maxima and their negations
rather than a scan over twelve keys.

The product y*T[h][t] is the dense blst_fp12_mul reference operation (operands
read through pointers, output written to a rotating buffer; no copies), charged 1.
Step cost c = 1 + 2,250/20,000 = 1.1125 units.

Fruitless-cycle overhead factor: 1 + 3/(12R) + 200/R^2 = 1.000435 (2-cycles at rate
1/(12R), <= 3 steps each; long cycles charged at rate R^-2, <= 200 steps each).
The walk-step cap of the algorithm is 1.000435 S steps in total, cycle handling
included, so this factor is a hard cap, not an average. By section 3.2 the actual
mean cycle-handling fraction is at most 3/(12R) + 200 x 0.139/R^2 = 0.000271 <
0.000435. The cycle events are sums over ~2^125.7 steps of rare events (means >=
2^100), so by a Chernoff bound the cycle-handling steps exceed 0.000435 S - 2^-30 S
(the second term reserves room for restart waste, section 3), leaving fewer than S
useful steps, only with probability < 2^-40 (counted in section 4).

Setup and output, all <= 2^31 units together: pairing g (Miller loop + final
exponentiation, < 2^12 units, bounded by 2^14), input check z^r = 1 (<= 510),
table: 1024 multipliers x (2 x 255-bit exponentiations <= 1,021 units + 11
Frobenius maps <= 1 unit each) <= 2^20.02, M = 2^20 random walk starts
g^alpha z^beta <= 2^30.0, output [x]P0 (<= 510 G1 operations, each < 1 unit) and
verification pairing (<= 2^14).

Distinguished-report handling, charged at its worst case (the deterministic cap, not the
mean), all with fixed instruction counts per record (section 3):

| Work per record (2^67 records at the cap) | Instructions |
| --- | ---: |
| append: 2 modular multiplications by p^kappa (512-bit products, reduction, <= 400 each), key extraction, 80-byte compare against the walk's last record, 80-byte store, counter compare, per-slot counter increment | 1,000 |
| full-record heapsort: <= 67 levels x 130 (10-word compare-and-swap of the 80-byte record), charged up from the computed 8,710 | 10,000 |
| adjacent dedup: one 80-byte comparison per record | 100 |
| **total per record** | **11,100** |

2^67 x 11,100 = 2^80.44 instructions = 2^66.15 units. The pair term is bounded by the
hard pair budget P = 2^41: P x 200 = 2^48.6 instructions = 2^34.36 units, plus one
canonical-element recomputation (1,021 units) per record a charged pair touches, at most
2 x P records = 2^42 x 1,021 < 2^52 units. The setup and output (2^31) and these three
terms total < 2^66.16 units, and are bounded by the round figure **2^66.2 units**. The
whole term is charged in full whether or not the reports occur. No separate charge
arises from the per-walk cap: the slot counter is incremented inside the append charge
above, and a stopped slot takes no further steps.

W = c * S * 1.000435 + c * E + 2^66.2
  = 1.1125 * 2^125.6711 * 1.000435 + 1.1125 * 2^84.3 + 2^66.2
  = 2^125.825569 (log2 S = 125.671137, log2 1.1125 = 0.153805, log2 1.000435 = 0.000627;
    125.671137 + 0.153805 + 0.000627 = 125.825569; the other two terms add < 2^-41.3
    relative: c*E/2^125.825569 = 2^-41.37 and 2^66.2/2^125.825569 = 2^-59.63, so log2 W
    rises by less than 2^-41.3/ln 2 < 2^-40.8 and stays 125.825569 to six decimals).

Submitted bound: **time_log2 = 125.826** = 125826/1000 >= 125.825569 (part C of
experiment.py recomputes log2 W = 125.8256 and checks claim >= log2 W; at 80-digit
precision 2^125.826 / W = 1.0002986 > 1 and 2^125.825 / W < 1, so 125.826 is exactly
the ceiling of log2 W at three decimals).

Memory (a hard bound, every live array counted): the report log of C_rep = 2^67
records of 80 bytes is 2^67 x 80 = 2^73.3219 bytes, allocated in full and never
exceeded; the heapsort works in place with one 80-byte temporary; the multiplier table
is 12,288 x 576 bytes = 7.1 MB (2^22.76) plus its exponent tables (2 x 12,288 x 32
bytes < 2^20); the walk states are 2^20 x <= 1 KiB = 2^30 bytes (an Fp12 element 576
bytes, two 32-byte exponents, three 48-byte key histories, an 80-byte last-record
buffer and counters, < 900 bytes).
Total < 2^73.4 bytes. There is no hash table, so no load factor and no table
beyond the log. Memory is reported, not scored.

## 6. Heuristics and limitations

* H1 (random-mapping behaviour of an R-adding walk on Frobenius classes) is the
  standard heuristic behind every rho estimate (Teske 1998/2001; Brent-Pollard;
  van Oorschot-Wiener 1999 for parallel distinguished points; Duursma-Gaudry-Morain
  1999 and Bos-Kleinjung-Lenstra 2010 for automorphism classes and fruitless
  cycles). It is supported here by the toy runs; it is not proved.
* Distinguished-point and partition bits are taken from K; K is the minimum of
  12 values, which is not uniform in [0,p), but its low 70 bits are pseudo-random
  functions of the class (K is about p/13 > 2^377, so the bits used are far below
  its size; heuristic, standard). The same assumption (bits 70..197 of K of distinct
  classes independent and uniform) prices the 128-bit truncated report key (false
  matches, section 3) and the probability that a healthy walk slot is stopped by the
  per-walk report cap (section 3). Support: the
  toy runs use exactly this rule on the raw toy key (K about 2^16.3 there, bits
  used <= 14) and measure the distinguished fraction equal to its nominal value, and
  part A2 of the experiment measures the low 70 bits of K on real BLS12-381 GT
  elements (section 7). The collision-rate model itself (random mapping with the
  Brent-Pollard correction) is supported by the toy success fractions at three
  field sizes, not proved.
* No production run. Toy experiments are evidence only at their own scale.
* This is a generic attack (about the square root of the class count). It
  ignores index-calculus methods in Fp12 (TNFS family); it is a valid upper
  bound, not a security estimate.

## 7. Experiments (`experiment.py`, standard library only)

Part A (BLS12-381 itself): for random GT elements (random Fp12 element raised to
(p^12-1)/r), checks y^r = 1; Frobenius formula = y^p; the 4-multiplication key
formula equals the keys of all 12 explicitly computed conjugates and they are
distinct; the sign-set minimum of section 2 (two differences, two reflections,
min/max of triples) equals the twelve-way argmin kappa and the class key K (the
structured rule used for the cost ledger); kappa selects the identical canonical
element from each of the 12 conjugates; p^6 is inversion; ord_r(p) = 12; 1 + p^4 + p^8 = 0 mod r. Part A runs on
every invocation with at least one production element (the code floors
`production_checks` at 1), so no run prints uncomputed, initialized-true flags.

Part A2 (BLS12-381 itself, bit uniformity): it raises one random GT element to 64
random exponents mod r (a pool of random GT elements), forms up to 2,000 samples as
products of three pool members (each a random GT element; time-limited to about
2.2 s), computes the class key K = min of the 12 keys with the four-multiplication
formula, and reports, for each of the low 70 bits of K (the 10 partition bits and the
first 60 distinguished bits), the frequency of a one bit and the largest deviation
from 1/2 in standard errors, with the mean deviation (0.80 expected for uniform bits) and
the fractions of K whose 4 and 6 bits above the partition bits are zero against 2^-4 and
2^-6. A run at the initial inputs gave 2,000 samples, largest deviation 2.77 standard
errors over 70 bits (about what 70 uniform bits give), mean deviation 0.90, and zero
fractions 0.0705 and 0.0175 against 0.0625 and 0.0156; a rerun differs in the
last digits because the sample count depends on the time limit. It supports the claim that
those bits behave as uniform; it cannot measure a 2^-60 event, which is why the
distinguished fraction is also measured in the toy.

Part B (toy, same tower shape): q = 7 mod 12 near 2^20 with xi = 1+u neither a square
nor a cube (asserted), r' a prime factor of q^4 - q^2 + 1 with ord_r'(q) = 12
(asserted), N = (r'-1)/12. It runs the production walk: four-multiplication keys,
kappa, twisted table, partition h = K mod R and distinguished test "bits
log2 R .. log2 R + dp_bits - 1 of K are zero" on the raw key, 2-cycle escape,
checkpoint escape and the stuck-walk rule at 40/theta, in the production step order
(both cycle checks before the distinguished test) and with the production append rule.
Toy-only differences: the 60 distinguished bits become dp_bits (the toy key has about
20 bits); ties of the minimum key (about 66/q per step in the toy) are broken by the
full conjugate.

Every trial runs its full fixed step budget 1.03 S0: no trial stops at its first
recovery. The collected report log is then matched by the production matcher, not by an
online index. The log holds 80-byte records (the truncated key T = K >> (log2 R +
dp_bits), which is the toy analogue of bits 70..197 of K, and the two canonical
exponents). The matcher sorts the log by the FULL record (T, then A, then B), so identical
records are adjacent; dedup is one adjacent full-record comparison per record; each maximal
equal-key run of d remaining distinct records then has all d(d-1)/2 of its pairs enumerated
and charged, with the canonical element g^A z^B recomputed for each record a charged pair
touches; x is solved from the first pair whose elements are equal with different exponent
pairs; and x is accepted only when exponentiation verifies g^x = z. Pairs are charged
against the toy's hard pair budget, which carries the production margin: P = 2^41 is 2^36
times the about 2^5 false-key pairs production expects over its 2^67-record log, and the
toy's budget is 2^36 times its own worst-case pair count (every pair in one run, about
(log capacity)^2/2). No healthy toy run can approach that, exactly as no healthy
production run can approach P; the abort rule is exercised by a separate tight-budget
stress, described below, not by the rows. The
recovery the matcher can return is a pair whose two elements are equal and whose exponent
pairs differ, which is the production condition; a repeat is charged as a pair like any
other. Two numbers are reported apart: the fraction of trials with a first useful class collision inside the budget (from a
measurement-only class index that no success flag reads) and the fraction whose
logarithm the matcher recovered and verified. The per-walk report cap of section 3 is
scaled to toy size (R_w / mean = 2.511 for each walk slot) and stops a slot that reaches
it; each row reports the log length against its cap, the slots stopped, the pairs charged
against the pair budget, and the element recomputations; a separate stress pass runs the
same walk code with a tiny cap (4 reports per slot) so the rule and the log bound are
exercised, not asserted, and a pair-budget stress matches a real collected log and a
two-record log with one equal-key pair at budgets 1 and 0, so the abort fires and is seen
to fire. The report
counter lives in the slot record and a stuck-walk restart carries it over, exactly as
section 3 states; a second stress pass sets the stuck-walk threshold to 2 steps so that
restarts certainly occur, and it reports the restart count and the log length against
M x cap, which is meaningful only because the counter persists (a counter reset by a
restart would let one slot append more than its cap).

Part C recomputes, without randomness, the zero-sum probabilities, the fruitless-cycle
bound of section 3.2, the deterministic report-log capacity M R_w, the probability that a
healthy walk slot reaches the cap, the memory bound, the per-record report-handling charge
(append, full-record heapsort, adjacent dedup), the hard pair budget and the charge it
bounds, and the W ledger of section 5.

Offline runs of this file (seeds 101-104, 400 trials per field per seed, "seconds": 230;
1,600 instances per field, 4,800 in all, each run to its full step budget and then matched
offline). A trial with no useful collision inside the budget counts as a failure for
P[T_c <= S0] and the mean collision time is censored; S0 is the 0.39 budget, model 0.39 at
S0 and 0.4081 at 1.03 S0. The last column of measured cycles is per step, with the derived
value in parentheses.

| q, r', R, dp_bits | trials | P[T_c <= S0] | P[T_c <= 1.03 S0] | recovered and verified | 2-cycles/step (1/(12R)) | long-cycle events, per step (x R^-2) | log max / cap, cap stops |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1048963, 48674953, 64, 4 | 1600 | 0.378 +- 0.012 | 0.396 +- 0.012 | 0.371 | 1.35e-3 (1.30e-3) | 50, 1.50e-5 (0.062) | 174 / 328, 0 |
| 1049611, 64668529, 256, 4 | 1600 | 0.385 +- 0.012 | 0.410 +- 0.012 | 0.385 | 3.23e-4 (3.26e-4) | 39, 1.02e-5 (0.67) | 191 / 376, 0 |
| 1054003, 35456053, 1024, 3 | 1600 | 0.385 +- 0.012 | 0.400 +- 0.012 | 0.393 | 9.36e-5 (8.14e-5) | 46, 1.63e-5 (17.1) | 277 / 556, 0 |

Measured distinguished fractions: 0.06271, 0.06397 and 0.12962 (nominal 2^-4 = 0.0625,
2^-4 and 2^-3 = 0.125). No walk was abandoned by the stuck-walk rule. The report log
stayed within its cap in every trial (largest 174, 191 and 277 records against caps 328,
376 and 556) and the per-walk report cap stopped no walk; the stress pass, run 12 times
with a cap of 4 reports per slot, stopped all four slots each time and left the log at
exactly 16 records, its cap. Both new rules are exercised: the append rule dropped 44,
0 and 0 repeats (the R = 64 field is the one whose small R produces longer cycles), and
the adjacent full-record dedup dropped 48, 18 and 19 exact duplicates before any pair was
formed (they are the same record, not a pair). The pairs charged stayed far below the
toy's hard pair budget in every trial: at most 106, 520 and 2,180 pairs in one trial
against budgets of 2^36 times the worst-case pair count (3.7e15, 4.8e15 and 1.1e16), with
no abort; 63,363, 336,272 and 1,446,428 pairs were charged in all, and the element
recomputations were at most 111, 180 and 276 in one trial. The tight-budget stress matched
a real collected log at a budget of one pair (one pair charged, the match aborted) and a
two-record equal-key log at budgets 1 and 0 (one pair charged without abort; zero pairs
charged and the abort fired), so the hard rule is exercised, not asserted.

The recovered-and-verified fraction is below the first-collision fraction by 0.025, 0.025
and 0.007 (40, 40 and 11 trials). Those are collisions so near the end of the step budget
that no distinguished class followed them, not matcher losses: every logarithm the matcher
recovered verified by exponentiation, and no trial recovered without a collision. In
production that gap is covered by the detection window E = M x 20/theta = 2^84.3 steps
that all walks run after the budget (section 4), whose failure probability is e^-20.

Long fruitless cycles: 50, 39 and 46 events were seen, at 0.062, 0.67 and 17.1 R^-2. The
cost ledger charges R^-2 per step with 200 steps per event, so the R = 64 and R = 256 rows
are inside the charge (0.062 and 0.67 of it). Against the sharper derivation of section
3.2 (0.139 R^-2, i.e. 115 and 8.1 events expected) they are 0.45x and 4.8x, in opposite
directions and both outside Poisson noise; that scatter at toy scale is why the cost uses
the round R^-2 with 200 steps per event rather than the sharp constant. The R = 1024 row
is 17x above the charged R^-2 and outside the toy's regime: there
R = 2^10 against S0 = 2^10.7 steps, so partition indices repeat inside one short walk by
necessity and the independence premise of the derivation cannot hold, while in production
R/S0 = 2^-115.6. This is a sanity check at the toy's own scale, not a validation of the
production rate: the production cycle-rate bound is the derivation of section 3.2 (the
union bound over index sets and the Phi_12 zero-sum probability), and no toy row is used
in its place.

The R = 1024 row is also the weakest for the success fraction: its class count is only
N_toy = 2^21.8, S0 = 2^10.7 steps and R = 2^10 is a large fraction of S0, so the R-adding
mixing model (which presumes R << S0) is strained there; in production S0 = 2^125.6 is
astronomically larger than R = 2^10. Its 0.400 at 1.03 S0 (1,600 trials) is close to the
0.4081 model. The quantitative toy support for the success fraction is the R = 64 and 256
rows and the larger offline table; the R = 1024 row is reported but not leaned on.

References: Pollard 1978; van Oorschot & Wiener, J. Cryptology 12 (1999);
Teske, Math. Comp. 70 (2001); Duursma, Gaudry & Morain, ASIACRYPT 1999; Wiener &
Zuccherato, SAC 1998; Bos, Kleinjung & Lenstra, ANTS 2010.
