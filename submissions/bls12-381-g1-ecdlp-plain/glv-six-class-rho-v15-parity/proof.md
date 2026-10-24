# Six-class parallel Pollard rho on BLS12-381 G1: negation + GLV, parity class key

**Claim.** time_log2 = 126.241 = 126241/1000. The algorithm recovers x from
T = [x]P0 with probability >= 0.39, and on every run its total work is
<= W_max = 2^126.241 units, one unit = one batched affine G1 addition =
C = 2,000 scalar AArch64 instructions. The bound is unconditional: a global meter
charges every action at an upper-bound price before it runs (section 3).
Memory < 2^72 bytes; preprocessing <= 2^57.2 units by a meter sub-cap, inside W;
no advice. Success rests on H1 (section 4). The work bound is
W_up = 2^126.240632 with W_up^1000 <= 2^126241 and W_up^1000 > 2^126240 (exact
integers, section 5).

## 1. Algebra (`experiment.py` part A checks every line on BLS12-381)

p = field prime, r = order (log2 r = 254.85709), E: y^2 = x^3 + 4;
beta = 0x1a0111ea...00000000aaac (primitive cube root of 1);
phi(x, y) = (beta x, y) = [lambda] on G1 with
lambda = 0xac45a4010001a40200000000ffffffff (phi(P0) = [lambda]P0 checked).
u_t = (-1)^s lambda^j (t = 3s + j) acts by [u_t](x, y) = (beta^j x, (-1)^s y);
no nonidentity G1 point is fixed by a unit != 1, so G1 \ {O} has N = (r-1)/6
classes of size 6 (log2 N = 252.27213).

**Parity class key.** K = y if y is even, K = p - y if y is odd. p is odd, so the
pair {y, p - y} holds one even and one odd element: "take the even one" is
invariant under y -> p - y and unchanged by phi (which leaves y). So K is constant
on every class and unchanged by every unit. Equal keys give y' = +-y, hence
x'^3 = x^3, x' = beta^j x: the same class. K is a **complete** class key and
**bit 0 of K is identically 0**. It replaces the accepted package's
|y| = min(y, p - y), whose exact six-limb comparison (rows 1-3 and 5-8 of
section 5 are unchanged from that package, including its truncated-product
canonical index, re-checked by part A on 40,000 production inputs with 0
mismatches) the parity bit makes unnecessary; that rule stays in `experiment.py`
as the reference rule, checked to induce the same partition of the orbits.

**Canonical member and bits.** c = (beta^j x, K), j = argmin(beta^j x mod p),
s = [y odd], t = 3s + j: c = [u_t]P (checked for all six members). Bit 0 of K is
forced, so h = bits 1..20 of K and the distinguished test is limb 1 (bits
64..127) of K = 0, theta = 2^-64; both lie above bit 0, so under H1 they are
uniform. fp(K) = limbs 0, 2, 3 (192 bits, 191 free).

## 2. Walk

R = 2^20; theta = 2^-64; 2^20 walks; L_max = 20 * 2^64; checkpoint period 64.
**Table:** c_h, d_h uniform in [0, r), M_h = [c_h]P0 + [d_h]T (redraw if O),
Tab[h][t] = [u_t^{-1}]M_h = (-1)^s phi^{-j}(M_h). **Step:** F(C) = class(c(C) + M_h),
h = bits 1..20 of K; the walk keeps a member P and sets P <- P + Tab[h][t], and
c + M_h = [u_t](P + Tab[h][t]), so this is F: one batched affine addition (the
reference operation) plus at most 117 instructions.

**Trails.** Start S = [a]P0 + [b]T, a, b uniform from fresh coins, stored in full
(redraw if S = O). At each class arrival, in order: if K is distinguished, report
(K, a, b) and end the trail; if length >= L_max, abandon; if fp(K) == K_cp,
traverse the cycle once more (<= 64 steps, else abandon) and resume at the member
of its minimum-fp class (no extra operation), then ESC; if length mod 64 == 0,
K_cp <- fp(K); finally P <- P + Tab[h][t]. Cycles are detected only by the
checkpoint (the valid flag is 0 at every start and after an escape). ESC is the
**doubling escape** P <- [2]P, a class map: doubling commutes with every unit, so
it sends class C to class(2c(C)). One fixed total order on classes (fingerprint
limbs 3, 2, 0) drives the checkpoint test and the traversal minimum; `length`
counts every operation (<= L_max + 65); equal x-coordinates abandon the trail
(< 2^-114 in total).

**Merging lemma.** Let A and B arrive at the same class. Until one escapes both
apply F and visit the same classes. A trail escapes only after K = K_cp, a repeat
within 64 ordinary steps, so only on a cycle Z of F of length <= 64 (a spurious
fingerprint match is a failure event); Z holds no distinguished class. A trail on
Z reaches a length that is a multiple of 64 within 63 arrivals, setting K_cp to a
class of Z, and is back there |Z| arrivals later, so the test fires; the traversal
visits all of Z once and doubles the member of min Z, reaching
E(Z) = class(2 c(min Z)), which depends on Z only, after <= 192 operations. So A
and B are again at a common class and by induction reach the same distinguished
class, unless one first hits L_max or an exceptional addition. No cycle closes
through doublings: the order of 2 mod r is 2^250.86 (part A), so 1 - u 2^m != 0
mod r for m < 2^248 and the table terms cannot cancel (that pins the class; H1
bounds coincidences).

**Solving.** On a key match with another trail's report, re-walk both from their
stored (a, b) with P = [A]P0 + [B]T, reducing mod r after every operation;
(e, f) = u_t (A, B). If the canonical points agree and f_1 != f_2:
x = (e_2 - e_1)/(f_1 - f_2) mod r; verify [x]P0 = T and output x. T = O gives 0.

## 3. Phases and the global work meter

Every action (random words, sampling attempts, starts, steps, doublings,
traversals, reports, probes, table build, hash initialisation, re-walks, solving,
verification) is charged its price (section 5) before it runs; each price
upper-bounds its instructions / C. A charge above W_max = 2^126.241 halts with
failure, so W <= 2^126.241 on every run. Preprocessing has its own cap 2^57.2.
The meter is a 192-bit counter; each of at most 2^20 processors reserves blocks of
2^40 units atomically (remainders count as spent, < 2^61 units, in the hard
extras). **Phase 1** runs until a global operation counter (steps, doublings,
traversal steps; re-walks separate) reaches S_ops, so at least S_ops and fewer than
S_ops + 2^60 operations run. **Phase 2** finishes open trails, no new trails.
Report capacity D_cap = 2^63 exactly (section 4) and at most 8 processed matches
(overflow = failure). **Batches:** 2^18..2^19 ready walks per Montgomery inversion
(a smaller batch pays a 118-unit inversion); the meter is charged once per batch
and the zero-denominator test runs once per batch on the product of the
denominators (a 2^17-unit recovery action then scans and abandons the offenders;
< 2^-114).

## 4. Success probability

N' = N * R/(R-1). T_c = the number of class arrivals of all trails (starts,
steps, doublings, traversal steps, repeats inside fruitless cycles: nothing
excluded) up to the first arrival at a class already visited by another trail.

**H1 (heuristic, at the budget).** P(T_c <= S_ops) >= 1 - exp(-(S_ops/kappa)^2/
(2N')) with kappa = 1.02 (a 2% slowdown allowance on the random-mapping law),
i.e. >= 0.3903. Toy slowdowns are a harness effect scaling with (mean trail
length)/sqrt(N); section 6 measures it: the measured mean T_c/sqrt(N) is
+0.08% +- 0.16% above the model (110,000 trials over five trail lengths) and
P(T_c <= 1.02 S0) = 0.4004 +- 0.0015, 6.7 standard errors above 0.3903, while the
earlier fit a = +0.09% +- 0.23% puts kappa = 1.02 about 8 errors above the value
extrapolated to production's trail ratio 2^-62. Also part of H1: partition and
distinguished bits of K are pseudo-random functions of the class; **excluding the
fruitless cycles of F of length <= 64** (2-cycles at rate 1/(6R) per step and
rarer short cycles, whose repeats are counted inside T_c), a trail makes a
*non-fruitless return* with probability at most L_max^2/(2N) (F-cycles of length
65..200 need every partition index repeated: rate < 2^-206 per step); and the
canonical twist index of a new class is uniform on the 6 units and independent of
its partition index.

S0 = sqrt(2 N' ln(1/0.6097)) = 2^126.12851; S_ops = 1.02 S0 = 2^126.15707 (H1
gives >= 0.3903; pure model 0.4024); 1.0585 (S_ops + 2^60) < 2^126.239096.

**Phase 1 covers the budget.** Phase 1 performs at least S_ops operations (hard
counter) and each completed operation is one class arrival, so if T_c <= S_ops the
first inter-trail collision happens in phase 1.

**Doubling and traversal cost.** Only the doublings' surcharge
(1.511 - 1.0585 = 0.4525) and the traversal steps' 0.009 are not fixed by the
operation count. 2-cycles need the next class to share the partition index and
have canonical twist -1 (probability 1/(6R) per step); longer fruitless cycles
have rate <= 2^-37.999. So the expected doublings are
< (1/(6R) + 2^-38)(S_ops + 2^60) = 2^103.57 (linearity over the pseudo-random
partition and twist of each new class), surcharge < 2^102.43, and the expected
traversal steps are < (2/(6R) + 64 * 2^-37.99)(S_ops + 2^60) < 2^104.57,
surcharge < 2^97.78.

**Report capacity, exactly.** Under H1 the distinguished bits of distinct classes
are independent 2^-64 coins, so new-class reports are Binomial(n, 2^-64) with
n <= S_ops + 2^60 + 2^104.6 (traversal arrivals) + 2^62.2 (starts), mean
< 2^62.17. D_cap = 2^63 is 1.79 times the mean and Chernoff gives
exp(-0.2522 * 2^62.17) < 2^-60; repeats of a class are matches, at most 8.

**Failure events** (unconditional union bound; if T_c <= S_ops the merging lemma
gives one report and the attack succeeds unless one occurs): trail cap hit first
(< 1e-7); a non-fruitless return in a colliding trail (< 2^-53); report overflow
(< 2^-60); f_1 = f_2 (1/r); exceptional additions (< 2^-114); 8 matches failing
(< 2^-60); spurious fingerprint match over 191 free bits (< 2^-63); meter halting
early (surcharges < 2^102.49 against the allowance 2^116.36: < 0.000054 by
Markov, 0.000161 at 3x cycle rates); preprocessing cap (< 2^-31).

**Success >= 0.3903 - 0.000161 - 1.1e-7 > 0.39013 > 0.39** (0.39024 at the model
cycle rate).

## 5. Cost

**Per ordinary step** (instructions besides the reference addition; blst Fp
multiplication = 383):

| Work | Instr. |
| --- | ---: |
| 1.1 v1: 3+4 loads, 6 `umulh` (column 5), 5 `add`, 6 `madd` (column 6) | 24 |
| 1.1 range test, V1, X, S, wrap test, V2 (3+5+1+6) | 15 |
| 1.1 decision (4 `cmp`, 6 `csel`, `cset`, `mov`, `sub`, `b.ls`) | 14 |
| parity bit, p - y limbs 0-3, select limbs 0-3 of K | 11 |
| h (bits 1-20), distinguished test, twist index, table address | 8 |
| loads of table point and state, store of state | 19 |
| checkpoint test, refresh test, refresh | 12 |
| length increment (69-bit), cap test, loop | 9 |
| zero-denominator test, per-batch meter charge, batch inversion share | 3 |
| **total, a per-step maximum** | **115** |
| **budget with 2 instructions of stated slack** | **117** |

**Row 4 (11)** replaces the accepted package's 20-instruction key row (rows 1-3
and 5-8 are its rows, unchanged; row 5 changes only `and x4,x0,#0xfffff` ->
`ubfx x4,x0,#1,#20`, same count):

```
row 4 (11)  ldp x0,x1,[x19,#48]; ldp x2,x3,[x19,#64]                  (y limbs 0-3)
            subs x7,x10,x0; sbcs x8,x11,x1; sbcs x9,x12,x2; sbcs x25,x13,x3  (p - y, limbs 0-3)
            ands x28,x0,#1                                            (s = [y odd], Z = [s = 0])
            csel x0,x7,x0,ne; csel x1,x8,x1,ne; csel x2,x9,x2,ne; csel x3,x25,x3,ne
```

The earlier row formed p - y over six limbs and compared it with y over six limbs
for the exact sign s = [p-y < y]. With the parity key s is bit 0 of y (one
instruction) and a borrow never propagates downward, so limbs 0-3 of p - y are
exact after four subtractions. A step reads only limbs 0-3 of K (limb 0 in row 5,
limb 1 in the distinguished test, limbs 0, 2, 3 in the fingerprint), so limbs 4-5
are not formed; the report recomputes all six inside its 0.5-unit price.
`work/notes/par2-listing-check.py` executes the row on a 64-bit register model
against exact big-integer arithmetic (200,000 random and 100,000 adversarial
inputs, 0 mismatches in the key limbs, the parity bit, the index, the
distinguished test, the table address).

**Prices** (units): step 1 + 117/2000 = 1.0585; traversal 1.0675 (135
instructions); doubling 1.511 (7 x 383 + 7 x 32 + 117 = 3,022); fallback 0.25
(454 <= 500); word 0.005; sampling attempt 0.025; start 2^12; report 0.5; probe
0.015; hash init 2^57.1 (once; unchanged, covers the smaller table); table 2^32.2
(once); phase-2/re-walk operation 119 (236,422 <= 238,000); small-batch inversion
118 (234,864); solving 2^9; verification 2^13.

**Budget.** Phase-1 operations 1.0585 (S_ops + 2^60) <= 2^126.239096 plus the
expected doubling and traversal surcharges 2^102.43 + 2^97.78. Hard extras: meter
remainders < 2^61, block reservations < 2^82.3, phase 2
(2^20 (L_max + 65) 119 = 2^95.217), re-walks 2^79.217, reports 2^62, hash init
2^57.1, table and verification < 2^32.3: together < 2^95.22. Other random charges
(starts 2^62.2 * 2^12, fallbacks 2^72.2, words/attempts/probes 2^63) < 2^74.4.
With the meter allowance A = ceil(3 (2^102.43 + 2^97.78 + 2^74.4)/0.000199) =
2^116.36, W_up = 1.0585 (S_ops_up + 2^60) + 2^95.22 + A = 2^126.240632 units,
where S_ops_up = ceil(1.02 sqrt(2 N' ln(1/0.6097))) =
94855758132035518167548215207775477484 (integer square root of an upper bound for
ln(1/0.6097)). In exact integers, with
W_2000 = 2117 (S_ops_up + 2^60) + 2000 (2^95.22 + A) = 2000 W_up:
W_2000^1000 <= 2^126241 2000^1000 and W_2000^1000 > 2^126240 2000^1000, so
126.241 is the exact ceiling of log2 W_up in thousandths and W_max = 2^126.241
covers W_up. Headroom at the cap: 2^126.241 - 1.0585 (S_ops_up + 2^60) - 2^95.22
> 2^116.67 = 2.05 A, so the meter halts early with probability < E/A <= 0.000161
even if the doubling and traversal rates are three times the model (Markov).

**Preprocessing and memory.** Jump table 2^32.2, report-table initialisation
2^57.1, first 2^20 starts 2^32: expected < 2^57.11, cap 2^57.2 (the random part is
< 2^21, far below the 2^53 gap: halts a run with probability < 2^-31); all inside
W. Memory: D_cap records x 112 bytes at load 1/2 = 2^70.81, plus table and walk
states: < 2^72.

## 6. Experiments and limitations

`experiment.py` (stdlib). **Part A** (production BLS12-381; all true): curve,
order, beta/lambda, orbit and unit action, parity key invariance and bit 0 zero,
the reference |y| rule's partition, canonical member, twisted step, doubling map,
order of 2 mod r, table entries, coefficient tracking, recovery of a planted x
from two forced-merge trails. **Part B** (toy y^2 = x^3 + 4, q = 58493047): the
algorithm of section 2; `"pipeline": true` applies the scaled production rules,
`"sweep": true` reports the T_c distribution and per-kappa fractions. Trails run
sequentially; [a]G + [b]T is evaluated as [a + bx]G; toy keys are below 2^64.

**Walk-margin sweep** (T_c mode, R = 2048, every recovered x verified,
binomial standard error 0.0033 per row, columns P[T_c <= kappa S0],
S0 = sqrt(2 N' ln(1/0.6097))):

| dp | trail/sqrt(N) | trials | mean T_c/sqrt(N) | 1.05 | 1.04 | 1.03 | 1.02 | 1.015 |
| ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: |
| 0 | 0.056% | 22000 | 1.2525 +- 0.0044 | 0.4170 | 0.4107 | 0.4052 | 0.3995 | 0.3968 |
| 1 | 0.111% | 22000 | 1.2533 +- 0.0044 | 0.4160 | 0.4101 | 0.4047 | 0.3989 | 0.3961 |
| 2 | 0.222% | 44000 | 1.2565 +- 0.0031 | 0.4188 | 0.4130 | 0.4069 | 0.4012 | 0.3983 |
| 3 | 0.444% | 22000 | 1.2541 +- 0.0044 | 0.4209 | 0.4146 | 0.4079 | 0.4013 | 0.3982 |

Model mean T_c/sqrt(N) = sqrt(pi N'/2)/sqrt(N) = 1.2536. Pooled (110,000
trials): P(T_c <= 1.02 S0) = 0.4004 +- 0.0015, 6.7 standard errors above the H1
floor 0.3903 and at the 0.40 gate; at 1.015 it is 0.3975 +- 0.0015 (4.9 se above
the floor, below the gate), so kappa = 1.02 is the smallest grid value the
measurement supports. The pooled slowdown is +0.08% +- 0.16%, consistent with the
earlier fit a = +0.09% +- 0.23%. Seeds 7001-7010, 7101-7110, 7501-7510,
7601-7610, 7701-7710, 7801-7806; `work/notes/par2.sweep.json`.

Traps and limits: the `long_trails` mode (q30, 150 trails of 2,048 operations
without distinguished points) gives 1/150 and 1/150 trails with a non-fruitless
return for doubling escapes at R = 4 and 8 (ordinary self-collision predicts 0.7
and 1.9) against 46/150 and 14/150 for the rejected additive escape. H1 is a
heuristic (random mapping with Brent-Pollard correction; Wiener-Zuccherato, GLV,
Duursma-Gaudry-Morain; Bos-Kleinjung-Lenstra and Bernstein-Lange-Schwabe escapes),
supported only at toy scale; instruction counts are analytic; no TNFS.