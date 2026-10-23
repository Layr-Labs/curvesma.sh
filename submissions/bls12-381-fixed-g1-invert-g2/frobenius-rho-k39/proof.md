# Frobenius-class parallel Pollard rho in GT (fixed G1 argument, output in G2)

Claim: time_log2 = **125.87** (total work W <= 2^125.87 dense-Fp12-multiplication
equivalents, C = 20,000), success probability >= 0.39 (the ledger below gives
0.4013), peak memory < 2^64 bytes (section 5), preprocessing <= 2^31 units
(included), no advice. Every array, index and **byte address** of the algorithm is
reached by a single 64-bit word: the report cap is 2^54 records of 640 bytes and the
whole peak is 11,817,445,424,400,702,208 bytes = 2^63.36 < 2^64, so each address,
offset and index is one native 64-bit load or store. The **counters are not one
word each, and this package does not claim they are**: a counter is as wide as the
cap it must reach, so the global slot counter (cap L_max + 2^61 = 2^125.6851) is 126
bits and the stuck-walk counter (cap 40 x 2^80) is 86 bits, and each is a two-word
register pair whose carry path is charged per increment in section 5. Storage below
2^64 bytes does not by itself make a counter one word; it only makes every *address*
one word. All randomness comes from one
fresh uniformly random 256-bit seed K drawn from the platform randomness source
at start (one 512-byte draw of 64 64-bit words, 2^16 instructions, itemised and
charged in section 5). Every sampled scalar is a 384-bit output of Keccak-f[1600]
in counter mode under K, reduced modulo r by the exact Barrett reduction at 8,192
instructions;
the sampling sites are hard-capped and nothing is rejected. That tape is a
computational premise and not a statistical one: it is under 2^49.98 bits
against a 256-bit seed, so no closeness to independent uniform draws is claimed
here, and section 6 states the premise as H3, the assumption that the tape is
indistinguishable from independent uniform draws at the claimed cost. The
success bound below is conditional on H3. The distinguished store is keyed on the
full exact canonical class key -- the canonical representative itself, all 12 Fp
coefficients, 4,608 bits -- so a key match is exactly a class match and no
distributional assumption about key bits is used anywhere (section 3). The twelve
keys the step compares are the Montgomery representatives blst_fp_mul returns, so
no conversion is needed or charged (section 2). Fruitless cycles are handled with
**no long-cycle-rate assumption anywhere** (section 3): a look-back window of the last
W = 24 class keys with A = 8 deterministic probes detects and escapes every repeat
inside the window; one far checkpoint key refreshed every 2^20 slots detects and
escapes every return within that many slots; the residual longer cycles are covered
by an explicitly charged hard budget line of at most 2^40 escapes of at most 2^21
slots each (2^61 slots, charged in W) and by the stuck-walk restart; and a walk
state that is the identity is detected by an eight-instruction test (six ORs of the
minimum key's limbs into an accumulator, a compare with zero and a branch) and
restarted, with its probability charged as the 2^-214.49 term of section 4. The success bound rests
on H1, the generic random-mapping model of rho walks, and on H3, the randomness
assumption: both are stated and scoped in section 6.
Every number below is reproduced as an exact rational by part C of
`experiment.py`, including the final integer certificate that W <= 2^125.87 and
that 125.86 is not enough, and a self-check that the ledger's quoted figures are
rounded in the safe direction (upper bounds up, lower bounds down).

## 1. Reduction

Input: nonidentity P in G1 (uniform, fixed for the instance) and z in GT. Compute
g = e(P, Q0) with the pinned convention (library call pairing(Q0, P)). P has prime
order r and the reduced pairing is nondegenerate, so g != 1 and g generates GT
(prime order r). Every valid z equals g^x for a unique x in [0, r). If z = 1,
output Q = infinity (e(P, infinity) = 1). Otherwise recover x by the rho algorithm
below, output Q = [x]Q0 and check e(P, Q) = z. By bilinearity
e(P, [x]Q0) = e(P, Q0)^x = g^x = z. The answer is for the supplied P; no other
P is used, and Q0 is the public profile generator.

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
  Since zeta^j in {1, zeta, zeta-1, -1, -zeta, 1-zeta}, all 12 keys come from
  **4 Fp multiplications** (a*zeta, a*Re(d_1), b'*Im(d_1), B*zeta) plus one Fp
  addition and 8 subtractions/negations. No Frobenius of the whole element is
  ever computed.

Define kappa(y) = argmin_k key_k and the step key K([y]) = min_k key_k. Then
y^{p^{kappa(y)}} is the same element for all 12 members of a class (canonical
representative y^), and both K and y^ are class invariants. K drives the walk
(partition index and distinguished test, section 3); the store key is y^ itself,
which is exact (section 3). A tie needs an equality of two of the 12 keys, which
for a uniform a (or B) is a linear coincidence modulo p: at most 66/p per step,
and at most 66 * L_max / p <= 2^-248.9 for the whole run (log2 p = 380.7007). A
tie is counted as failure.

**Montgomery form of the keys (no conversion is charged, and none is needed).**
blst keeps every Fp element internally in Montgomery form, x~ = x * 2^384 mod p
fully reduced in [0, p), and the four products above are exactly the
blst_fp_mul calls the step makes: blst_fp_mul(x~, y~) returns (x y)~ = x y
2^384 mod p, fully reduced, at the published 383 instructions (COST_MODEL, which
counts the final conditional subtraction that makes the result reduced). The one
Fp addition and the eight subtractions/negations commute with the Montgomery map,
because that map is multiplication by the fixed unit 2^384. So the twelve values
the step computes are the Montgomery representatives of the twelve keys, with no
reduction to canonical form, no export and no conversion anywhere, and the keys
are compared as the integers they are, i.e. as k * 2^384 mod p in [0, p).
x -> x * 2^384 mod p is a bijection of Fp (2^384 is invertible modulo the prime
p), so the minimum is attained at the same index, with the same comparisons and
the same ties, as it would be on canonical residues: K is a class invariant, the
tie bound above is unchanged, and the step key's bits (partition index and
distinguished test) are a fixed function of the class in either form. Nothing in
sections 3 to 5 depends on which representative is compared; the four products
are charged at 383 instructions each, exactly as before.

## 3. Walk (R-adding walk on classes, Duursma-Gaudry-Morain / Wiener-Zuccherato style)

Parameters: R = 4096 partitions, distinguished-point fraction theta = 2^-80,
M = 2^20 parallel walks, each started at g^alpha z^beta with alpha, beta from the
bounded sampler of section 5. All walks share one table and one class map.

Precompute M_h = g^{a_h} z^{b_h} (h < R, random a_h, b_h from the same sampler) and
the table T[h][k] = M_h^{p^k}, k = 0..11 (49,152 entries of 576 bytes, 28,311,552
bytes), with exponent tables a_h p^k, b_h p^k mod r (49,152 entries of 64 bytes).

Class map: F([y]) = [ y^ * M_h ], h = K([y]) mod R. A walk keeps an arbitrary
member y of the current class with its exponents (alpha, beta), y = g^alpha z^beta.
Because y^{p^kappa} * M_h = (y * M_h^{p^{-kappa}})^{p^kappa}, a slot is

```
compute keys, kappa, K                      (4 Fp mul + <= 2,316 other instructions)
if the minimum key is 0 (identity rule below): restart the walk at a fresh start
if bits 12..91 of K are all zero (distinguished): report (K, p^kappa*alpha, p^kappa*beta)
h = K mod R (bits 0..11) ; t = (12 - kappa) mod 12
for a = 0, 1, ..., A-1:                     (A = 8, deterministic probes)
    cand = y * T[(h+a) mod R][t]            (ONE dense Fp12 multiplication)
    reject if cand is the identity, or its class key is the far checkpoint key,
      or its class key is one of the W = 24 window keys: charge one slot, next a
    else accept: y = cand ; alpha += a_h p^t ; beta += b_h p^t (mod r) (one slot)
if all A probes were rejected: restart the walk at a fresh start
```

A slot is one reference operation plus the key work, accepted or rejected. Accepted
slots are the walk's steps; a rejected slot is one wasted slot, and the slot counter
of the hard step cap counts both. The walk is a deterministic function of its class
and its window, so once two walks reach the same class with the same window they
follow the same class sequence and reach the same distinguished class (the window
exception is charged in section 4). A walk continues after reporting; it is never
restarted after a distinguished point.

**Distinguished store (batched, deterministic, exact key).** The walks never
consult a store. Each report is appended to a flat array. Its key is the full
canonical class key: the canonical representative y^ = y^{p^kappa} written out in
all 12 Fp coefficients (4,608 bits, stored in 72 64-bit limbs, 576 bytes), plus
the canonical exponents A = p^kappa alpha and B = p^kappa beta mod r as two
32-byte residues, 640 bytes per record. Nothing is truncated, hashed or
projected; a key match is therefore exactly a class match (two classes share a
key only if they are the same class). The canonical representative is computed
with one coefficient-wise Frobenius pass using the 12 stored constants
d_{m,kappa} = xi^{m(p^kappa-1)/6} for the selected kappa: 12 conjugations and 12
Fp multiplications, budgeted 5,120 instructions per report. A global counter
caps the array at 2^54 records; the algorithm stops with failure at the next
report. When the step budget is spent, the array is ordered by the full key with
a fixed LSD radix sort: 288 stable passes of 16 bits over the 4,608-bit key, each
pass a counting sweep and a scatter sweep over one 2^16-entry counter array,
permuting an array of 8-byte indices that is read from one buffer and written to
a second buffer of the same size (the two buffers swap roles every pass; the
first is filled with the identity order). **Every address, offset and index is one
64-bit word**: the largest address is the last byte of the record array,
640 x 2^54 < 2^63.33, the record index needs 54 bits, each of the 2^16 digit counters
holds a count of at most 2^54 records and so is 55 bits in one 8-byte word, and the
peak of the whole algorithm is 11,817,445,424,400,702,208 bytes =
2^63.36 < 2^64 (section 5), so each of them is one native AArch64 load/store and no
multiword pointer arithmetic is used anywhere. The *counters* are sized by their
caps and not by the address width: the report counter (cap 2^54, 55 bits), the digit
counters (cap 2^54 each, 55 bits), the partition index (12 values, 4 bits), the probe
counter (4 bits), the far-checkpoint counter (21 bits) and the escape counter
(cap 2^40) are one 64-bit word each, while the global slot counter (126 bits) and
the stuck-walk counter (86 bits) are two words each, and section 5 charges their
carry path. The cap is 2^54
and not 2^56 because a 640-byte record at 2^56 would need 640 x 2^56 = 2^65.32
bytes for the records alone, beyond a 64-bit address space; 2^54 is the largest
cap with the exact 576-byte key under 2^64 bytes, and section 4 pays its
report-cap tail. The pass count, the per-pass instruction count and the counter
array are constants, so the handling work is a fixed function of the record
count: there is no probing, no load factor, no rejected insert and no discarded
report (section 5). One linear scan then reads the sorted indices in order, keeps
the first record of each run of equal key, and compares every later record of the
run with it and with its predecessor. If the exponents equal the stored ones the
record is dropped (same walk history); otherwise x = (A1 - A2)/(B2 - B1) mod r is
computed and checked against g^x = z. On success the algorithm outputs Q = [x]Q0;
on failure it keeps the first record.

Detection is complete and exact: two reports of one class carry the same key,
because the canonical representative is a class invariant, so they lie in one
run, and every record of a run is compared with the run's first record, so the
pair is tested and no report is ever skipped. Different classes never share a
key, so no run mixes two classes and no comparison is wasted or misleading. The
relation is exact: with equal keys, p^k1(alpha1 + beta1 x) = p^k2(alpha2 + beta2
x) mod r, and the stored exponents are the canonical ones. For two distinct
exponent pairs of one class, B1 = B2 would force A1 = A2 as well, because
g^{A1} z^{B1} = g^{A2} z^{B2} means A1 + B1 x = A2 + B2 x mod r; so the division
never has a zero denominator and such a record pair is a duplicate that is
dropped. Detection after the fact cannot miss a useful collision either: no walk
ever reads the store, so no walk stops or diverges early and each walk takes
exactly the steps charged in section 4.

**Hard report cap:** a global counter counts every report, stored or matched;
the algorithm processes at most 2^54 reports and stops with failure at the next
one. Hence at most 2^54 records are ever stored, and no report is lost to the
store before that cap.

**Hard global step cap:** every slot of every kind (accepted steps, rejected
probes, steps of walks that are later abandoned) increments one global counter. At
L_max = 1.01 S + E = 2^125.6851 slots (S, E in section 4) the algorithm stops with
failure. L_max needs 126 bits, so this counter is a two-word register pair (section
5 charges the carry path of every increment).

**Fruitless cycles: no rate is assumed.** A twisted adding walk can return to a
class with formally cancelled exponents. This package handles every repeat by a
deterministic rule and never by a probability estimate about cycle lengths.

*Look-back window.* Each walk holds the class keys of its last W = 24 accepted
states in a sorted array. Every slot's candidate class is compared with those keys
before it is accepted, so a candidate that repeats any of the last 24 classes is
rejected and the next probe index is tried. Hence **any repeat of a class visited
within the last 24 accepted steps is seen at the next slot and escaped there**, for
every class and every history, with no distributional assumption and no rate. In
particular every cycle of length <= 24 is escaped at its first return (its next
class lies in the window) at the cost of one rejected slot. Resource cost, per walk:
24 keys of 48 bytes = 1,152 bytes, plus 24 one-byte slot ids and two one-byte ring
pointers (26 bytes, carried in a 32-byte block), so 1,152 + 32 = 1,184 bytes, which
is carried inside the 2,048-byte walk state of section
5 (2^20 x 1,184 = 1,241,513,984 bytes in total). Per slot: five binary-search probes
over the 24 sorted slot ids, each 26 instructions (address 2, id load 1, key address
3, six limb loads 6, six-limb compare 12, branch 2) = 130; inserting the accepted
candidate's class = shifting at most 23 one-byte ids (69), one id store (1), one key
write (11) and two ring-pointer updates (4) = 85; so 130 + 85 = 215 instructions,
itemised in section 5.

*Far checkpoint.* Each walk also holds one far checkpoint key, refreshed every
F = 2^20 slots to the walk's current class key, and a candidate whose class equals
it is rejected. A cycle of period p <= F is therefore escaped: within F slots of the
refresh the checkpoint key lies on the cycle, the walk returns to it at most p <= F
slots later, and that candidate is rejected, one slot. So every cycle of period up
to 2^20 slots is escaped within at most 2F = 2^21 slots of its entry, deterministically,
and the cost is one rejected slot per escape. The compare is six limbs per slot, 12
instructions, itemised in section 5.

*Hard budget line for the escapes.* The far-checkpoint escapes are capped and
charged as their own budget line, separate from the slot cap: **at most 2^40
escapes, each costing at most 2^21 slots of walk (at most 2F slots of detection
delay and one rejected slot), so at most 2^61 slots, which is 1.13 x 2^61 = 2^61.18
units charged in W** (section 5). A global escape counter enforces the cap and the
run stops with failure if it is reached. The cap is never approached: under H1 the
expected number of within-walk repeats in the whole run is at most (L_max^2)/(2 N M)
<= 2^-20.9, so the probability of reaching 2^40 escapes is at most 2^-60.9, which is
inside the 2^-20 self-collision charge of section 4 and is not counted twice.

*What covers anything longer.* A cycle whose period exceeds 2^20 slots is seen by
neither the window nor the far checkpoint, so it is bounded by the stuck-walk rule
instead: it can cost a walk at most 40/theta = 40 x 2^80 slots before the walk is
abandoned and restarted, each such restart consumes one of the n_r < 2^40.37
restarts, and the waste is charged by the abandonment tail 2^-35.05 of section 4.
**No long-cycle rate is claimed or used.**

*Identity rule.* A walk state is the identity when all twelve of its keys are zero;
that state has twelve equal keys and the minimum rule below is undefined on it, so it
is handled explicitly rather than left to a tie. The identity test is eight
instructions (six ORs of the minimum key's limbs, a compare with 0 and a branch); when it fires
the walk is abandoned and restarted at a fresh start g^alpha z^beta, and the event is
charged: a fresh start is the identity with probability exactly 1/r, because its two
exponents are uniform and alpha + x beta = 0 mod r is one linear condition with x !=
0, so over the M + n_r + 1 start and restart sites of the run the probability is at
most (M + n_r + 1)/r = 1,414,115,699,111/r = 2^-214.49; the 2^20 walk starts alone
give 2^20/r = 2^-234.85. The restart costs 1,021 units and two samples, which the
restart cap n_r and the sampler site cap of section 5 already cover, and the walk
restarts from a fresh start. A zero minimum at a *step* state is either the identity
or a tie among keys, which is the tie term 66 L_max/p <= 2^-248.9 of section 2; the
rule turns it into a restart, so that term remains a valid (conservative) failure
bound. The identity event enters the failure tail once, as 2^-214.49.

**Stuck-walk rule:** a walk that has taken 40/theta = 40 x 2^80 slots since its
last distinguished point (or its start) is abandoned and restarted at a fresh
g^alpha z^beta (two samples of section 5). Since every abandoned walk consumed
40 x 2^80 counted slots, restarts number at most n_r = floor(L_max/(40 x 2^80)) + 1
= 1,414,114,650,534 < 2^40.37, deterministically. The restart cap n_r is the cap on
restarts of every kind, stuck-rule and identity alike, so the sampler site count and
the restart term of section 5 cover both.

## 4. Success probability

Heuristic (H1): the class sequence of the walk behaves like a random mapping on
N classes, with the Brent-Pollard correction for R-adding walks (effective
N' = N * R/(R - 1), log2 N' = 251.2725). Then the probability that no useful class
collision occurs among the first S points of all walks is exp(-S^2/(2N')). Choose

S0 = sqrt(2 N' ln(1/0.61)) <= 2^125.6280,   S = 1.03 * S0 <= 2^125.6707,

where the 3% is explicit slack for walk non-randomness. At S0 the model gives
exactly 0.39; at S it gives 1 - 0.61^(1.03^2) = **0.408**. The track requires
0.39 at the submitted bound, and the losses below are charged against the 0.018
of margin between 0.408 and 0.39.

Detection: a class collision is only visible once both merged walks report the
same distinguished class. After the first useful collision the merged trajectory
reaches a distinguished class within 40/theta further slots except with
probability e^-40, and both walks then report that class, so the pair appears in
the array and the final scan compares it. Charging all M walks in lockstep over
that window costs E = M*40/theta = 40 x 2^100 < 2^105.33 slots. This is why the
slot cap is 1.01 S + E: after paying the detection window and at most 0.010 S of
waste the run has S useful slots.

**Both tails are Markov, not Chernoff.** The per-slot indicators below are not
independent: the class at each slot is a deterministic function of the previous one
and of the window, the partition index and the distinguished test are fixed
functions of that class, and a walk can repeat classes (section 3), so a Chernoff
bound over L_max slot indicators is not available. Linearity of
expectation is, and it absorbs all of that looping: whatever the dependence
between slots, the expected total waste is the sum of the per-slot expected
waste, and a bounded random variable with mean mu exceeds t with probability at
most mu/t (Markov).

Slot waste (this is the only waste rate the ledger uses): useful slots are at least
L_max - E - waste, and every wasted slot is a rejected probe. The only structural
repeat a probe can reject at a rate above 24/N is the exact two-cycle, which needs
h' = h and kappa' = 6 (inverse conjugate) and therefore occurs at rate 1/(12R) =
2^-15.58 per slot at R = 4096, measured 1.267e-3 against 1.302e-3 predicted at
R = 64 in the toy; it costs exactly one rejected slot, since the next probe index
leaves the cycle. A structured cycle of length 3 to 24 needs, beyond the index
repeat the two-cycle already needs, at least one further index coincidence or a
vanishing twist sum, so the allowance 8/R^2 charges at most A - 1 = 7 further
rejected slots at a rate a factor R below the two-cycle's; it is a conservative
charge and not an estimate of any cycle length. A candidate that hits the window or
the far key for no structural reason occurs at 25/N < 2^-246. So a slot rejects
with probability at most
1/(12R) + 8/R^2 + 25/N = 2.082e-5 whatever its neighbours do, the run's mean number
of wasted slots is at most 2.082e-5 L_max = 2.082e-5 x 1.01 S, and Markov gives

  P(rejected slots > 0.008 S) <= 2.082e-5 x 1.01 / 0.008 = 0.00263  (2^-8.57),

which is the allowance charged below. Nothing here is a claim about how long cycles
are: the window, the far checkpoint and the escape budget of section 3 fix the
waste, and this tail only bounds how many rejected probes that costs.

Report cap: under H1 each slot lands on a distinguished class with probability
theta = 2^-80, whatever the dependence between slots, so the run's expected number
of reports is at most theta L_max <= 2^45.69, only 2^8.31 below the cap of 2^54.
Markov then gives P(reports >= 2^54) <= 2^45.69 / 2^54 < 2^-8.31, a loss of
0.0032 charged below. (The cap is 2^54 rather than 2^56 because the store must
stay below 2^64 bytes with 640-byte records; section 3.)

Abandoned walks: a 40/theta stretch with no distinguished point has probability
at most e^-40 <= 2^-57.7. Every stretch begins at a report, a start or a restart,
so there are at most 2^54 + 1 + 2^20 + n_r < 2^54.02 stretches, and the expected
number of abandonments is at most 2^54.02 x 2^-57.7 = 2^-3.68 <= 2^-3.67 (0.077).
Each abandoned walk strands at most 40 x 2^80 counted slots, so the abandonment
waste exceeds 0.002 S only if at least 0.002 S/(40 x 2^80) >= 2^31.38 walks are
abandoned; by Markov's inequality that has probability at most
2^-3.67/2^31.38 = 2^-35.05.

Coalescence: a meeting is only useful if the two walks stay together. Their windows
differ in at most W entries, so at one of the next W slots the probe rule can send
them to different classes -- once for each of the two walks, at the two-cycle rate
used above -- and after W further slots the windows are identical again. The
probability that a meeting is lost that way is at most 2W/(12R) = 48/49,152 = 2^-10,
charged below.

No report is lost to the store: every report is written to the array before any
comparison and every record of a run is compared against the run's first record
and its predecessor. Because the key is the exact canonical representative, a run
never contains two classes, so the store contributes no success loss and needs no
distributional premise; the array and its sort also remove the load-factor and
discard analysis of the previous design entirely.

Losses charged against the 0.408, each a probability that the run fails:
probe waste 0.00263 (Markov, above), report cap 2^-8.31, abandonment 2^-35.05,
the e^-40 detection tail 2^-57.7, ties 2^-248.9, the identity term 2^-214.49
(section 3), the coalescence term 2^-10, and the 1/M = 2^-20 self-collision
fraction. No sampling term is charged, because
under H3 (section 6) the sampled scalars are indistinguishable from independent
uniform draws at this cost, so sections 2 to 5 compute the failure probability of
the uniform-coin algorithm; H3 is a computational assumption and the bound is
conditional on it. Every one of them is deducted, including the 2^-20 that
the 3% step slack was originally for:

  success >= 0.408089 - 0.002629 - 0.003140 - 2^-35.05 - 2^-57.7
             - 2^-248.9 - 2^-214.49 - 2^-10 - 2^-20
           = **0.4013** >= 0.39.

The only randomness in the algorithm is consumed at the sampling sites of section
5; every slot is deterministic given its class and its window, so no slot consumes
coins and the slot cap cannot be outrun by retries. The one fresh seed is drawn once, at
start, and its draw is charged as an additive term in section 5.

## 5. Cost (units of one dense Fp12 multiplication = 20,000 instructions)

Per-slot non-reference instructions, scalar AArch64, using the published 383
instructions per blst_fp_mul. A slot is one dense Fp12 multiplication plus this
work, whether the probe is accepted (a walk step) or rejected (one wasted slot):

| Work | Itemised hard cap |
| --- | ---: |
| 4 Fp multiplications (a*zeta, a*Re(d_1), b'*Im(d_1), B*zeta) at 383 | 1,532 |
| 1 Fp addition (B), 6 limbs with carry and conditional subtract | 24 |
| 8 Fp subtractions/negations, the six zeta^j multiples from a and from B, 24 each | 192 |
| min of 12 six-limb keys with index select: 12 comparisons with 6 loads, 6 subtracts and 2 selects each, plus 16 for the running minimum | 192 |
| identity/tie test: OR the six limbs of the minimum key, compare with 0, branch | 8 |
| distinguished test (bits 12..91 of K), partition index, t = (12-kappa) mod 12, table address | 40 |
| look-back window search: 5 binary-search probes over the 24 sorted slot ids at 26 (address 2, id load 1, key address 3, 6 limb loads, six-limb compare 12, branch 2) | 130 |
| look-back window insert: shift at most 23 one-byte ids (69), one id store (1), one 48-byte key write (11), two ring-pointer updates (4) | 85 |
| far checkpoint compare: six-limb compare of the current class key with the stored key | 12 |
| counters: 126-bit slot counter (two words, carry path) 8, 86-bit stuck-walk counter (two words) 7, 21-bit far-checkpoint counter 4, probe counter 2 | 21 |
| alpha, beta += table values mod r: 2 loads of 32 bytes, 2 four-limb adds with conditional subtract | 64 |
| loop control | 16 |
| **total itemised** | **2,316** |
| **budgeted per slot** | **2,600** |

The 284-instruction difference between the itemised 2,316 and the budgeted 2,600
is an explicit margin for load and store scheduling, the register traffic of the
table entry and the multiply operands, and any micro-operation the table does not
name. It is a budget, not a measurement: each row above is a worst-case count for
a specified operation, and every row is a hard cap. The window rows replace the
48-instruction 2-cycle/checkpoint/history row of the previous version of this
package (the look-back window replaces the checkpoint escape), so the itemised
total rises from 2,108 to 2,316 and the per-slot budget stays at 2,600.

The counters are sized by their caps, not by the address width. The slot counter
reaches L_max + 2^61 = 2^125.69, so it is 126 bits: two 64-bit words, incremented
with ADDS on the low word and ADC on the high word and compared with a two-word
compare (8 instructions). The stuck-walk counter reaches 40 x 2^80 = 2^85.33, so
it is 86 bits, also two words, with the same carry path (7 instructions). The
far-checkpoint counter (21 bits), the probe counter and the escape counter fit one
word each. The carry path of both wide increments is charged per slot in the
counter row above; a counter is not made one word by living in an address space
below 2^64, only an address is.

The product y*T[(h+a)][t] is the dense blst_fp12_mul reference operation (operands
read through pointers, output written to a rotating buffer; no copies), charged 1.
Every counted slot therefore costs at most c = 1 + 2,600/20,000 = **1.13** units,
and the walk costs at most c * L_max; the long-cycle budget line of section 3 adds
c x 2^61 = 2^61.18 units to W.

**Randomness: one fresh seed, then a bounded sampler.** No slot draws randomness.
At start the algorithm draws one fresh uniformly random 256-bit seed K from the
platform randomness source: one 512-byte draw = 64 64-bit loads and stores = 128
32-bit words, of which the first 8 (four 64-bit words) are the 256-bit seed. The
draw is charged whole, because the charge covers any platform API granularity:

| Fresh seed draw | Itemised hard cap |
| --- | ---: |
| one getrandom(2) syscall for 512 bytes, entry, validation and return | 4,096 |
| the kernel CSPRNG filling 512 bytes at <= 64 instructions per byte | 32,768 |
| 64 64-bit word loads and 64 64-bit word stores with pointer updates and loop control | 768 |
| absorb the four 64-bit seed words, the counter base and the padding into the Keccak state | 1,920 |
| sampler buffer, counters and control | 8,192 |
| **total itemised** | **47,744** |
| **budgeted for the one seed draw** | **65,536 = 2^16** |

Coins after that are consumed only at these sites, and the count at each is a hard cap:

| Site | Scalars sampled |
| --- | ---: |
| multiplier table entries, 4096 x 2 | 8,192 |
| walk starts, 2^20 x 2 | 2^21 |
| stuck-walk restarts, at most n_r x 2, n_r < 2^40.37 | < 2^41.37 |
| one spare scalar pair | 2 |
| identity or tie restarts (inside the n_r restart cap above) | 0 |
| cycle escapes (a rejected probe draws nothing) | 0 |
| **total, hard cap: 2 x (4097 + 2^20 + n_r) = 2,828,231,406,414** | **< 2^41.39** |

The parenthesised count is the number of scalar pairs: 4,096 multiplier-table pairs
(the 4,096 x 2 entries above), 2^20 walk-start pairs, n_r restart pairs and the one
spare pair, so the site list above reads 8,192 + 2^21 + 2 n_r + 2 = 2 x (4,096 + 2^20
+ n_r + 1) = 2 x (4097 + 2^20 + n_r) = 2,828,231,406,414 samples < 2^41.37 < 2^41.39.

One scalar is one 384-bit output of Keccak-f[1600] in counter mode under K (the
seed absorbed once, a 128-bit block counter, six squeezed words) reduced modulo r.
There is no rejection, no retry and no data-dependent branch, so the sampler's
work is a fixed product of the site counts:

| Work per sample | Itemised hard cap |
| --- | ---: |
| Keccak-f[1600], per round: theta 55 (20 column-parity XORs, 5 rotations and 5 XORs for D, 25 lane XORs); rho and pi 49 (24 rotations, 25 moves); chi 100 (NOT, AND and XOR for each of 25 lanes, plus 25 row-copy moves); iota 2; loads, stores and spills 90 (the 25-lane state reloaded and stored once, 20 temporaries spilled and refilled); round loop 4. That is 300 per round, x 24 rounds | 7,200 |
| absorb the seed K and the 128-bit counter, pad, squeeze six words | 64 |
| 6 x 64-bit loads of the 48-byte draw, buffer and loop control | 48 |
| Barrett: y * mu, 6 x 3 limb products at 6, 9-limb carry propagation (mu = floor(2^384/r), 130 bits) | 126 |
| qhat * r, 3 x 4 limb products at 6, carry | 80 |
| y - qhat r over 7 limbs: subtract, borrow, select | 21 |
| two conditional subtractions of r, 4 subtractions, 4 selects and a borrow flag each | 18 |
| shifts, final select, store, loop control | 24 |
| **total itemised** | **7,581** |
| **budgeted per sample** | **8,192** |

Every Keccak row counts each logical operation as its own instruction and charges
a full reload and store of the state every round, so it is an upper bound: scalar
AArch64 fuses a rotation into EOR and NOT into BIC, which only lowers the count.

The reduction is exact. With mu = floor(2^384/r) and y < 2^384, qhat =
floor(y mu / 2^384) satisfies q - 1 <= qhat <= q = floor(y/r), so y - qhat r < 2r
and two conditional subtractions return y mod r.

**H3, the randomness premise.** The tape produces at most 2^41.39 x 384 < 2^49.98
bits from a 256-bit seed, so it is far longer than any seed-to-output statistical
distance argument could cover, and no such argument is made: this package claims
no statistical closeness between the tape and independent uniform draws, and no
distance figure appears anywhere. What it states instead is H3: under a fresh
uniform seed K, the Keccak-f[1600] counter-mode output is computationally
indistinguishable from independent uniform bit strings at the cost of this attack
(2^125.87 work), and the platform seed is uniform. Under H3 every sampled scalar
is an independent uniform 384-bit value, the exact reduction above applies as
stated, and the success bound of section 4 is the bound computed for uniform coins.
H3 is a scoped computational assumption, not a proved property: the bound is
conditional on it, and if the platform source or the permutation fell short of it
the bound would not apply. The deterministic per-sample cost, the site caps and
the absence of rejection are charged facts, not assumptions, and nothing in
sections 2 to 5 depends on the sampling device beyond H3.

Additive terms (all deterministic bounds, each charged once):

* Setup and output, <= 2^31 together. Pairing g = e(P, Q0): the Miller loop is 64
  iterations of one Fp12 squaring, two line multiplications counted as dense and
  one G2 doubling at <= 1.44 units, plus one G2 addition for each of the 6 set bits
  of the ate loop, so <= 64 x 4.44 + 6 x 1.44 = 292.8 units; the final
  exponentiation even by plain square-and-multiply on the 4,314-bit exponent
  (p^12-1)/r costs <= 8,628 units; the pairing is <= 8,921 units. Table: 4096
  multipliers x (2 x 255-bit exponentiations <= 1,020 units + 1 product + 11
  Frobenius maps <= 1 unit each + the 22 multiplications mod r of the exponent
  tables, 22 x 512 = 11,264 instructions, <= 1 unit) = 4096 x 1,033 = 4,231,168 units
  < 2^22.02. M = 2^20 walk starts g^alpha z^beta at <= 1,021 units each (the two
  exponentiations and the product; the two draws are charged in the sampler term
  below), 1,070,596,096 <= 2^30. Output Q = [x]Q0 by double-and-add,
  <= 510 Jacobian G2 operations,
  each <= 20 Fp2 multiplications (<= 3 Fp multiplications + 5 Fp additions each)
  plus <= 20 Fp2 additions and control, <= 28,780 instructions = 1.44 units, so
  <= 735 units, plus one Fp2 inversion (Fermat, <= 762 Fp multiplications,
  <= 15 units). Verification pairing e(P, Q) <= 8,921 < 2^14. The itemised sum is
  8,921 + 4,231,168 + 1,070,596,096 + 735 + 15 + 8,921 = 1,074,845,856 < 2^31;
  the declared cap 2^31 is what enters W.
* Restarts: at most n_r x 1,021 units with n_r < 2^40.37, so <= 2^50.37.
* Sampler: at most 2,828,231,406,414 < 2^41.39 samples x 8,192
  instructions / 20,000 <= 2^40.11 units.
* Distinguished-point handling: at most 2^54 + 1 reports (report cap), stored as a
  flat array of 640-byte records, sorted by the full key and scanned (section 3).
  **Every address, offset and index below is one 64-bit word**, because the
  peak is 2^63.36 bytes and the cap is 2^54 records (section 3): a record address
  is 640 x index + base (a shift-add chain, 6), an index is 8 bytes, and a digit
  counter is one 8-byte word holding a count up to 2^54 (55 bits). The report
  counter reaches 2^54, so it is 55 bits and fits one word; the round
  counter of the scan and the pass counter are 8 bits. So each is one load or
  store, and the itemisation is:
  A report first forms its canonical exponents A = p^kappa alpha and B = p^kappa
  beta mod r, and the canonical representative y^{p^kappa} with one
  coefficient-wise Frobenius pass. One multiplication mod r is itemised as: the
  4 x 4-limb product by the stored constant p^kappa mod r, 16 limb products at 6
  plus carries (112); a Barrett reduction with mu = floor(2^512/r), the high four
  limbs times the five-limb mu, 20 limb products at 6 (120); qhat x r, 16 limb
  products at 6 (96); the subtraction and two conditional subtractions (40);
  loads, stores and control (40); 408 itemised, budgeted 512, so 1,024 per report
  for the two exponents (itemised 2 x 408 = 816). The canonical representative is
  12 Fp values, each a conjugation (12) and one Fp multiplication by the stored
  constant d_{m,kappa} (383) plus store and select (13): 12 x 408 = 4,896
  itemised, budgeted 5,120. The append writes the record: 80 word stores (72 key
  limbs and 8 exponent words), the 64-bit record address (6), the bound check
  against 2^54 (3), the counter increment and store (3) and loop control (6):
  98 itemised, budgeted 112 instructions per report. The first index buffer is
  filled with the identity order, 1 store plus the counter step and loop control,
  8 itemised, budgeted 12 per record. Per record per pass the counting sweep reads
  the 8-byte index (1 load), forms the record address as 640 x index plus the base
  (6), extracts the 16-bit digit from the 576-byte key (offset shift and add, one
  load, one zero-extend: 8), loads, adds and stores one 8-byte counter with its
  address (5) and runs the loop (3): 23. The scatter sweep reloads the index,
  address and digit (15), reloads the counter (1), forms the 8-byte output address
  8 x counter plus the base (6), stores one index word (1), increments and stores
  the counter (3) and runs the loop (3): 29. That is 52 itemised, budgeted 64, per
  pass. The scan reads 72 key limbs and compares them with the run's first
  record, keeps the run position, and compares the 8 exponent words on a
  duplicate: 480 itemised, budgeted 512. So one record costs
  112 + 1,024 + 5,120 + 12 + 288 x 64 + 512 = 25,212 budgeted instructions,
  itemised 98 + 816 + 4,896 + 8 + 288 x 52 + 480 = 21,274. Each 2^16-entry
  counter array (8-byte counters) is zeroed and prefix-summed once per pass, at
  most 16 instructions per entry, 288 x 16 x 2^16 <= 2^28.17 instructions in
  total. The whole handling term is (2^54 + 1) x 25,212 + 2^28.17
  <= 2^68.63 instructions <= 2^54.35 units. Neither the pass count nor the
  per-pass cost depends on the data.
* Candidate solves: at most two per record (the run's first record and the
  predecessor are each a comparison), <= 520 units each (inversion mod r plus
  exponentiation g^x, <= 510 dense multiplications, plus the comparison):
  <= (2^54 + 1) x 1,040 <= 2^64.03 units.
* Fresh seed: the one 512-byte (64 64-bit word) platform draw of section 5,
  <= 2^16 instructions / 20,000 <= 2^1.72 units (itemised 47,744). It is drawn
  once, before any table entry or walk start.

Additive total T <= 2^31 + 2^54.35 + 2^50.37 + 2^64.03 + 2^40.11 + 2^1.72
<= 2^64.04. Each term is quoted as a power of two rounded up, and the sum of the
quoted terms is itself at most 2^64.04; part C evaluates the same terms exactly
(T <= 2^64.0243, 18,759,128,372,271,595,423 units by the printed ceiling, dominated
by the candidate solves).

W <= c * (1.01 S + E) + c * 2^61 + T
  <= 1.13 * (1.01 * 2^125.6707 + 2^105.33) + 2^61.18 + 2^64.04
  <= 2^125.8614 < 2^125.87.

The second term is the long-cycle budget line of section 3 (2^40 escapes of at most
2^21 slots, charged as slots at c = 1.13); it is 2^61.18, more than 62 bits below
c x 1.01 S, so it does not move the grid point. The exact value is
log2 W = 125.861287 (part C).

Exact form: with u = 0.494296321814780115... an exact rational upper bound on
ln(1/0.61) from the atanh series (2 atanh(39/161), 40 terms plus a geometric
tail), S0 = sqrt(2N'u) rounded up, and E and the restart count exact integers,
part C of `experiment.py` evaluates A = c * 1.01 * 1.03 and B = c * (E + 2^61) + T as
exact rationals and checks (A S0 + B)^100 <= 2^12587, which is true, and
(A S0 + B)^100 > 2^12586, which is also true, so (A S0 + B)^100 <= 2^12586 fails
and 125.87 is exactly the smallest multiple of 0.01 at or above the bound, with
0.0087 bits of slack. Nothing is rounded down anywhere in that chain; the two
exact checks are printed by part C.

Submitted bound: **time_log2 = 125.87**.

Memory (hard bound by the record cap; supplementary on this track). Every array
live at the peak, which is the scatter phase of the radix sort (the walk states
are never needed after the step budget is spent but are counted anyway):

| Array | Bytes |
| --- | ---: |
| report records, 2^54 x 640 (exact 576-byte class key + two 32-byte exponents) | 11,529,215,046,068,469,760 |
| scatter input index buffer, 2^54 x 8 | 144,115,188,075,855,872 |
| scatter output index buffer, 2^54 x 8 | 144,115,188,075,855,872 |
| digit counters, 2^16 x 8 | 524,288 |
| multiplier table T, 49,152 x 576 | 28,311,552 |
| exponent tables a_h p^k and b_h p^k, 49,152 x 64 | 3,145,728 |
| Frobenius constants d_{m,kappa}, 144 x 48 | 6,912 |
| walk states, 2^20 x 2,048 (element 576, window block 1,184 = 24 class keys 1,152 + ids and byte ring pointers 32, far checkpoint key 48, exponents 64, counters 64, scratch 112) | 2,147,483,648 |
| scalars, caps, checkpoints, sampler tape and reduction, code and stack slack | 1,048,576 |
| **peak** | **11,817,445,424,400,702,208** |

The peak is 5 x 2^61 + 2^58 + 2,180,520,704 bytes =
11,817,445,424,400,702,208 bytes = 2^63.36, below 2^64, so **every byte of the
store is reachable by a single 64-bit address**: the last byte of the record array
is 640 x 2^54 < 2^63.33, the largest index needs 54 bits, and the whole run's
slot counter needs 126 bits only in its two-word register pair, never as an
address into memory. That is the point of the 2^54 cap: a 2^56 or 2^67 cap with
640-byte records would need 2^65.32 or 2^74.32 bytes of records alone, which no
AArch64 load or store can address. A peak below 2^64 makes the addresses one word;
it does not make the counters one word, and section 5 sizes every counter by its
own cap.
The 2^20 miscellaneous allowance (scalars, caps, code, stack and tape buffer) and
the 2^19-byte digit-counter array are both in that sum, together with the walk
states, the multiplier and exponent tables and the Frobenius constant set: those
non-record arrays are 2,180,520,704 bytes, below 2^31.03. The records and both index buffers are 656 x 2^54 bytes.
The exact class key is what makes detection exact, so the record is 640 bytes
instead of a 112-byte projected key; memory is supplementary on this track and is
not scored, so the exact key costs nothing in the score and the cap it forces
(2^54 rather than 2^56) costs the 0.0032 report-cap tail charged in section 4.
This is a hard bound by the report cap.

## 6. Heuristics and limitations

* H1 and H3 are the two premises that carry the success bound. H1: the class
  sequence of the walk behaves like a random mapping on N classes, with the
  Brent-Pollard correction for R-adding walks, and the partition function
  h = K mod R and the distinguished test on bits 12..91 of K behave as the fixed
  pseudo-random functions of the class that this model assumes. H1 is the standard
  heuristic behind every rho estimate (Teske 1998/2001; Brent-Pollard; van
  Oorschot-Wiener 1999 for parallel distinguished points; Duursma-Gaudry-Morain
  1999 and Bos-Kleinjung-Lenstra 2010 for automorphism classes and fruitless
  cycles). It is supported here by the toy runs; it is not proved. No separate
  claim is made about the distribution of any individual bit of K, and none is
  used: the store key is the exact canonical representative, so detection never
  depends on how K or any of its bits is distributed. H1 also covers the one place
  the probe rule changes the walk: a rejected probe sends the walk to the next
  table index, which happens at the two-cycle rate 2^-15.58 per slot, far inside
  the 3% step slack and inside the 0.008 S rejected-slot allowance.
* H3 is the randomness premise, stated as a scoped computational assumption with
  its own label rather than a statistical one. The algorithm draws one fresh
  uniformly random 256-bit seed K from the platform randomness source at start,
  in one charged 512-byte (64 64-bit word) draw, and the sampler's tape is
  Keccak-f[1600] in counter mode under K. H3: that tape is computationally indistinguishable from
  independent uniform draws at the cost of this attack, and the platform seed is
  uniform. The tape is under 2^49.98 bits against a 256-bit seed, so no
  statistical closeness to independent uniform draws is claimed anywhere and no
  distance figure is quoted; the success bound of section 4 is conditional on H3,
  the way published packages on this track state the same premise. The
  per-sample cost, the site caps and the absence of rejection are charged facts,
  and the tape is not assumed to interact with the class map.
* The rejected-slot, record-count and abandonment estimates in section 4 use H1,
  and their tails are Markov bounds on those means rather than Chernoff bounds over
  slots (section 4), because the per-slot indicators are not independent; the caps
  themselves, the restart count, the sampler charge, the per-record instruction
  cost and every cost term in section 5 are deterministic. Cycle handling rests on
  no rate and on no heuristic: the look-back window, the far checkpoint, the
  escape budget and the stuck-walk restart are deterministic rules with a hard
  slot charge, and the only rate the ledger uses (the exact two-cycle, 1/(12R)) is
  a structural condition that the toy measures. Detection is exact and rests on no
  heuristic at all: the store key is the canonical representative.
* The keys are the Montgomery representatives of section 2, compared as integers.
  That is not a heuristic: multiplication by 2^384 is a bijection of Fp, so the
  minimum, the tie condition and the class invariant are identical to the
  canonical-residue formulation, and nothing is converted or exported.
* No production run. Toy experiments are evidence only at their own scale.
* This is a generic attack (about the square root of the class count). It ignores
  index-calculus methods in Fp12 (TNFS family); it is a valid upper bound, not a
  security estimate.

## 7. Experiments (`experiment.py`, standard library only)

Part A (BLS12-381 itself): for random GT elements (random Fp12 element raised to
(p^12-1)/r), checks y^r = 1; Frobenius formula = y^p; the 4-multiplication key
formula equals the keys of all 12 explicitly computed conjugates and they are
distinct; kappa selects the identical canonical element from each of the 12
conjugates; p^6 is inversion; ord_r(p) = 12; 1 + p^4 + p^8 = 0 mod r. The pairing
reduction of section 1 is pure bilinearity and is not part of the toy.

Part B (toy, same tower shape): q = 163, r' = 54,298,861 dividing q^4 - q^2 + 1,
N = (r'-1)/12. Runs the full walk (key formula, kappa trick, twisted table,
**the production cycle rule** -- a look-back window of the last W = 24 class keys
with A = 8 deterministic probes, the identity reject, and one far checkpoint key
refreshed every FAR slots -- distinguished points, stuck-walk restart after
40/theta slots as in production, collision solving and verification) over many
random instances. The probe set is a Python set rather than a sorted array: the
rule (reject a candidate that repeats one of the last 24 classes) is identical and
the data structure is an implementation detail; FAR is 64 slots instead of 2^20 so
that a toy trail of about 640 slots can exercise a refresh. Reported: slots to
first useful class collision / sqrt(N) (random mapping predicts mean 1.2533), the
fractions of trials within the 0.39 and 0.40 step budgets, the rejected-slot
fraction, the measured two-cycle rate against 1/(12R), and the counts of
window rejections, far-checkpoint escapes, all-probes-rejected restarts and
identity starts. The toy does not apply the two production caps; its counter
counts slots, not setup work. The
store replay uses the production map exactly: the report key is the full canonical
representative of the toy field, every Fp coefficient serialised with no hash, so
equal toy keys mean equal toy classes just as in production. The toy's partition
index and distinguished test come from a fixed hash of the canonical element, while
production takes them from bits of the step key K (at q = 163 the step key is a
residue below 163, fewer than the 12 partition bits production uses). Both are
fixed functions of the class and both are covered by the same random-mapping
heuristic H1 of section 6; no distributional property of K's bits is claimed or
needed anywhere.

Part C (exact ledger): recomputes N, N', the exact rational upper bound on
ln(1/0.61), S0, S, L_max, the per-slot itemisation (2,316 against the 2,600 budget)
with its window, far-checkpoint and counter rows, the counter widths (126 and 86
bits, two words each), the long-cycle budget line (2^40 escapes of at most 2^21
slots, 2^61.18 units), the setup, record-array, restart, solve, sampler and
fresh-seed terms, W and the memory total (the 2^20 miscellaneous allowance, the
walk states with their 1,184-byte windows, the 2^16 x 8-byte digit counters, the
multiplier and exponent tables and the Frobenius constant set included) from exact
integers and fractions, prints log2 W,
the priced 512-byte (64 64-bit word) seed draw and the sampling-site count it
drives, and checks the
final integer inequality (W^100 <= 2^12587 and W^100 > 2^12586), and prints the
rejected-slot, report-cap and abandonment Markov tails, the identity and
coalescence terms and the final success sum. It computes no statistical-distance
figure, because the randomness premise is the scoped assumption H3 and not a
distance bound. It is the machine-checked version of sections 4 and 5. Its
self-check tests the 44 log2 bounds that part C itself lists, under the names it
lists them as: S0, S, L_max, E, the restart cap, the samples, the sampler term,
the table setup, the setup sum, the non-record array total, the handling
instructions and the handling term, the counter work, the restarts term, the
solves term, the additive total T and the same total at 4 decimals, W at
2^125.8614, at the claim 2^125.87 and at 2^125.86, the fresh-seed draw and its
2^16 instruction budget, the mean reports, the report-cap tail, the probe-waste
tail, the reject-rate bound, the mean wasted slots, the stretches, e^-40, the
expected abandonments, the abandonment tail and the abandonment threshold, the
restart samples, the ties, the two-cycle rate, the coalescence, the identity
term over all start sites and over the 2^20 walk starts, the long-cycle slot cap
and its charge in units, the window bytes per walk, and the three memory rows.
Each upper bound (42 of the 44) is tested so that its exact value is at most 2
raised to the printed exponent, and each lower bound (the 2 left: W at 2^125.86
and the abandonment threshold) so that it is at least, by an exact
integer-power comparison. Three of
the 44 are printed in this text in other units rather than as 2^x: the
reject-rate bound as 2.082e-5, the mean wasted slots as the product
2.082e-5 x 1.01 S, and the window as 1,184 bytes. The list is not the set of
every 2^x form in the text: the tape length 2^49.98 and the 2^-20 self-collision
fraction are not tested, and neither are the exact powers the text uses as inputs
(2^20 starts, the 2^54 report cap, theta = 2^-80 and 2^384 outputs and the like),
which the plain-number checks cover instead. With 3 descriptive constants at 4
decimal places that makes 47 figures in all; it re-evaluates 47 plain-number
claims (sums, itemisations, ratios, the two counter widths) and recomputes 17
chains from the printed, rounded operands as a reader would; it prints all_safe
true.

The toy store replay uses the production store shape and nothing else: every
distinguished report of a trial is appended to a flat array keyed by the full
canonical class key (all 12 coefficients, serialised exactly, no hash), the array
is ordered by a fixed LSD radix sort (8-bit digits, 24 passes) and one linear scan
compares the records of each run of equal keys.
Its first 20 trials recover the same logarithm from the sorted array that the
online walk finds, 20 of 20, over 3,251 records.

The initial input (q = 163, R = 64, 4 walks, seed 1, 40 trials, 14 s limit) is
deterministic given its seed, takes about 4 s of CPU and completes all 40 trials:
every logarithm was recovered and verified; the mean first-collision time is
1.1530 sqrt(N) (model 1.2533); the fraction within the 0.39 slot budget is 0.475
and the fraction within the 0.40 budget is 0.475 (models 0.39 and 0.40); 2-cycles
run at 1.2673e-3 per slot (predicted 1.3021e-3); of 100,213 slots 127 were
rejected, all of them exact two-cycles, with no far-checkpoint escape, no walk lost
to eight rejections, no identity start, no repeat report and no abandoned
walk; the store replay recovers 20 of 20 trials. The 0.39 budget is 1.0021 sqrt(N)
and the 0.40 budget 1.0188 sqrt(N). A runner slower than about three times this one
hits the 14 s limit, completes fewer trials and reports different sample values; the
ledger of part C is unaffected either way.

Larger offline runs of this same file (time limit lifted, so these inputs are
outside the runner's 14 s schema; production_checks 0), all logs recovered and
verified; standard errors from the random-mapping spread 0.655/sqrt(n) and
sqrt(0.24/n):

| q, r', R, dp_bits, walks, seed | trials | mean T_c/sqrt(N) (model) | P[T_c <= S0] (model 0.40) | 2-cycles/slot (1/(12R)) | far escapes |
| --- | ---: | ---: | ---: | ---: | ---: |
| 163, 54298861, 64, 4, 4, 1 | 2400 | 1.2862 +- 0.013 (1.2533) | 0.388 +- 0.010 | 1.3004e-3 (1.3021e-3) | 175 |
| 787, 64571173, 256, 4, 4, 2 | 600 | 1.1983 +- 0.027 (1.2533) | 0.440 +- 0.020 | 3.178e-4 (3.255e-4) | 14 |
| 331, 923346397, 1024, 6, 8, 3 | 96 | 1.1738 +- 0.067 (1.2533) | 0.479 +- 0.050 | 9.24e-5 (8.14e-5) | 1 |

The three larger runs are quoted at the 0.40 slot point, S0(0.40), which is the
column the same file has always printed; the 0.39 point is 0.9838 times smaller
and the initial-input run above reports both fractions (0.475 each, models 0.39
and 0.40). The claim of this package is the 0.39 point S0(0.39) = 2^125.6280 with
the 1.03 slack on top, so what these runs evidence is the random-mapping model
itself, at toy scale, not the value of one quantile.

The measured two-cycle rate tracks 1/(12R) at all three sizes (1.3004e-3 against
1.3021e-3, 3.178e-4 against 3.255e-4, 9.24e-5 against 8.14e-5), which is the one
rate the ledger uses. The far checkpoint fires in all three runs (175, 14 and 1
times): every one of those is a return the 24-key window could not see, escaped by
one rejected slot, so the rule of section 3 is exercised and not only argued. No
walk was lost to eight rejections and none was abandoned in any run.

References: Pollard 1978; van Oorschot & Wiener, J. Cryptology 12 (1999);
Teske, Math. Comp. 70 (2001); Duursma, Gaudry & Morain, ASIACRYPT 1999; Wiener &
Zuccherato, SAC 1998; Bos, Kleinjung & Lenstra, ANTS 2010.
