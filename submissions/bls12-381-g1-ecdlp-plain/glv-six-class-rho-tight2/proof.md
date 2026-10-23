# Six-class parallel Pollard rho on BLS12-381 G1 (negation + GLV)

**Claim.** time_log2 = **126.293** = 126293/1000, the exact three-decimal ceiling
of log2 W (section 5). x is recovered from T = [x]P0 with probability >= 0.39 for
every uniform instance and fresh coins; total work <= W_max = 2^126.293
units on every run (1 unit = 1 batched affine G1 addition, C = 2,000 scalar
AArch64 instructions), because a global work meter charges every action an
upper-bound price before it runs and halts at W_max. Peak memory < 2^72 bytes;
preprocessing <= 2^57.2 units under a meter sub-cap, inside W; no advice. Success
rests on heuristic H1 (section 4); README.md has the twist-index windows and the cycle argument.

**Four numbers change** from public submission
1c600df4-b531-401a-b9d8-5124c6bceeb2 (six-class rho, 126.36), whose design,
algebra, walk, lemma, meter, caps, prices and H1 are restated here: (1) the
per-visit ledger falls from 196 to 128 instructions, a step costing
1 + 128/2000 = 1.064 units, not 1.10; (2) the success target is the floor 0.3905,
so S_ops = 1.05 sqrt(2N' ln(1/0.6095)) = 2^126.19938; (3) R is 2^16, not 2048, so
1/(6R) = 2^-18.585, the expected doubling surcharge is 2^106.417, the headroom
Markov needs is 2^117.46, and the table costs 2^28 units and 6R entries of memory;
(4) the trail cap stays 20 theta^-1.

## 1. Algebra

p = 0x1a0111ea...ffffaaab, r = 0x73eda753...00000001 (log2 r = 254.857089), E:
y^2 = x^3 + 4. beta is a primitive cube root of 1 in Fp; phi(x,y) = (beta x, y) =
[lambda] on G1, lambda^2 + lambda + 1 = 0 mod r. Units u_t = (-1)^s lambda^j,
t = 3s + j, act by [u_t](x,y) = (beta^j x, (-1)^s y); no non-identity G1 point is
fixed by a unit, so G1\{O} has exactly N = (r-1)/6 classes of 6
(log2 N = 252.2721269). K = |y| = min(y, p-y) is a complete class invariant
(equal |y| gives y' = +-y, hence x' = beta^j x); the canonical member
c = (beta^j x, |y|) with j minimal is unique. t is decided exactly from bits
384..447 of x*floor(2^448 beta/p): 12 limb products give V1 within 3 of
(beta x mod p)/2^320 and V2 within 5 of (beta^2 x mod p)/2^320, so an argmin
separated by more than 16 is exact and a near-tie (rate < 2^-52 per step if
beta x mod p is equidistributed) falls back to the exact rule, a separately
metered action of 0.25 units. Part A of `experiment.py`: 0 mismatches.

## 2. Walk

R = 2^16; theta = 2^-64 (limb 1 of K zero); M = 2^20 walks; L_max = 20 theta^-1;
checkpoint 64. Table: c_h, d_h uniform in [0,r) from fresh coins,
M_h = [c_h]P0 + [d_h]T (redraw if O), Tab[h][t] = [u_t^{-1}]M_h
= (-1)^s phi^{-j}(M_h), coefficients u_t^{-1}(c_h,d_h) mod r. Step:
F(C) = class(c(C) + M_h), h = K mod R; the walk keeps any member P and sets
P <- P + Tab[h][t], which is F since c + M_h = [u_t](P + Tab[h][t]). A step is one
batched affine addition (the reference operation) plus the ledger of section 5.
Trails start at S = [a]P0 + [b]T with fresh uniform a, b, stored in full (redraw
if S = O). At each class arrival:

```
if K distinguished: report (K,a,b); end trail
if length >= L_max: abandon trail
if pending_escape: ESC; clear; continue
if fp(K) == K_{-2}: if fp(K) < K_{-1}: ESC; clear; continue
                    else pending_escape <- true
elif fp(K) == K_cp: traverse the cycle once more (<= 64 steps, else abandon),
     resume at its minimum-fp member (no extra operation), ESC; clear; continue
if length mod 64 == 0: K_cp <- fp(K)
P <- P + Tab[h][t]; K_{-2} <- K_{-1}; K_{-1} <- fp(K)
```

"Clear" resets K_{-1}, K_{-2}, K_cp, pending_escape. ESC is the doubling escape
P <- [2]P, a class map: doubling commutes with every unit, so C -> class(2c(C)). fp(K) = limbs 0, 2, 3 of K as an integer; every test and minimum uses
this one order, so trails' histories are comparable, and full keys are never
compared. `length` counts every operation, so a trail runs at most L_max + 65.

**Merging lemma.** A and B at one class with any histories follow F until one
escapes. An escape needs a repeat (K = K_{-2}, or K = K_cp after <= 64 steps),
i.e. a cycle Z of F of length <= 64, which then has no distinguished class (a
spurious fingerprint match is a failure event); it leaves Z's minimum-fp class,
directly or after the one pending step, for E(Z) = class(2 c(min Z)), which
depends on Z only. If A escapes from Z, B enters Z, detects it within 3 arrivals
(2-cycle) or 128 + 64 (checkpoint) and escapes to the same E(Z); they reunite
after <= 192 extra operations of one, and by induction reach the same
distinguished class unless one hits L_max or an exceptional addition. A fixed
point gives the same escape class by both routes. Escapes close no structural
cycle: README.md shows the order of 2 mod r, 2^250.8571, forbids cancellation.

**Solving.** On a key match re-walk both trails from their stored (a,b), tracking
P = [A]P0 + [B]T mod r; if the canonical points agree and f_1 != f_2,
x = (e_2-e_1)/(f_1-f_2) mod r, verified by [x]P0 = T. T = O gives 0.

## 3. Phases and the global work meter

Every action in the price list of section 5 is charged before it runs, each price
upper-bounding its instructions / C; a charge over W_max = 2^126.293 halts the
attack with failure, so W <= 2^126.293 always. The meter is a 192-bit
counter; at most 2^20 processors (one per walk) reserve 2^40-unit blocks
atomically, unused remainders counting as spent: < 2^61 units. **Phase 1** runs
until a global operation counter (steps, doublings, traversal steps; re-walks
separate) reaches S_ops, so S_ops <= operations < S_ops + 2^60. **Phase 2**
finishes the open trails, no new trail. The report table holds <=
D_cap = 2^63.25 records; at most 8 matches are processed (overflow = failure).
Batches: 2^15..2^16 ready walks per Montgomery inversion (a smaller batch pays
132).

## 4. Success probability

N' = N R/(R-1) = 2^252.2721489. T_c counts class arrivals of all trails (starts,
steps, doublings, traversal steps, repeats inside fruitless cycles: nothing
excluded) up to the first arrival at a class already visited by another trail.

**H1 (heuristic, at the budget).** P(T_c <= S_ops) >=
1 - exp(-(S_ops/1.05)^2/(2N')), the 5% slowdown allowance against the
random-mapping law (toy class walks -0.2% to +2.1%, a control matching it). Also H1: |y|'s bits are pseudo-random in the class; apart from fruitless cycles of
F of length <= 64 (2-cycles at rate 1/(6R) per step, escaped and counted in T_c;
longer ones need every partition index repeated, rate < 2^-30 per step here), a
trail's *non-fruitless return* to an earlier class has probability
<= L_max^2/(2N), the random-mapping self-collision rate (long-trail tests of the
source package); and a new class's canonical twist index is uniform on the 6
units, independent of its partition index.

At the target: S0 = sqrt(2N' ln(1/0.6095)) = 2^126.12899, S_ops = 1.05 S0
= 2^126.19938 (= 97678699983940186850688299877457044734 truncated),
U1 = 1.064 S_ops = 2^126.28888; the pure model gives exactly 0.3905. Phase 1
performs >= S_ops class arrivals (hard counter), so T_c <= S_ops is reached inside
phase 1; fruitless-cycle arrivals lie inside T_c.

**Surcharges.** 0.436 per doubling at rate 1/(6R) + 4R^-2 = 2^-18.585 + 2^-30 per
step (2-cycles need the next class to share the partition index with canonical
twist -1; a longer fruitless cycle needs every partition index twice), and 0.012
per traversal step at rate < 2^-30 with <= 64 steps. Linearity
of expectation over each new class's pseudo-random partition and twist gives
expected surcharges 2^106.417 and 2^99.819 (2^4 a margin on the Stirling sum). Toy
2-cycle counts agree with 1/(6R).

**Failure events** (union bound, unconditional; if T_c <= S_ops the lemma gives
both trails one report, so the attack succeeds unless one occurs): cap first (Gamma(2) tail, 21e^-20 = 4.3e-8); a non-fruitless return in a colliding
trail (H1: 2 L_max^2/(2N) = 2^-115.6); report overflow (new-class reports are
Binomial(n, 2^-64), n <= S_ops + 2^60 + 2^88.32 + starts, mean < 2^62.1994;
D_cap = 2^63.25 is 2^1.0506 times it, Chernoff < exp(-2^60.78); repeat reports of
a class are matches, <= 8 processed); f_1 = f_2 (1/r); exceptional additions
(< 2^-114); 8 matches failing (< 2^-60); a spurious fingerprint match (< 2^-64 by
union over < 2^127.5 comparisons); preprocessing cap (< 2^-32); early halt
(section 5). All but the early halt total < 1.2e-7.

**Success.** The headroom above the deterministic charges at W_max = 2^126.293 is
2^117.838 and the expected random charges are 2^106.432, so by Markov the meter
halts early with probability <= 2^106.432/2^117.838 = 2^-11.4063 = 0.000368437.
Hence success >= 0.3905 - 0.000368437 - 0.00000012 = 0.3901314 > 0.39 (margin
1.3e-4, exact arithmetic).

## 5. Cost

**Per ordinary step**, besides the reference addition (`ldp`/`stp`: two 64-bit
words per instruction):

| Work | # |
| --- | ---: |
| point load (6 ldp) | 6 |
| p - y (6 subs), s = bit 380, K select (6 csel) | 13 |
| L': 12 limb products | 12 |
| accumulate L' | 12 |
| L' range test | 4 |
| V1 (umulh, mul, add, lsr) | 4 |
| X, S, wrap test | 6 |
| V2 (sub, sub, csel) | 3 |
| argmin, two gap tests, fallback | 13 |
| h, t, table address | 5 |
| distinguished test | 2 |
| addend load (6 ldp) | 6 |
| exceptional-addition test | 7 |
| history: 2 tests, k2<-k1, k1<-fp, checkpoint | 23 |
| point store (6 stp) | 6 |
| length++, cap test, meter add, loop | 6 |
| **total, price 1 + 128/2000 = 1.064** | **128** |

Nothing is omitted and nothing is charged twice: the reference addition's 6M
covers the field arithmetic, the operand accesses and the batch inversion
(COST_MODEL: 3M formula + 3M batch), so only loads, stores, addressing, key,
index, history and loop control appear. Coefficients are not updated per step; a
re-walk recomputes them at 133 per operation.

**Prices** (units): step 1.064; traversal 1.076 (128 + <= 24 for copying the
remembered member); doubling 1.5 (7M + 128 = 2,809); fallback 0.25; random word
0.005; sampling attempt 0.025; trail start 2^12; report 0.5; hash probe 0.015;
hash initialisation 2^57.1; table 2^28 (2^16 rows at 2^12, once, plus sampling);
phase-2 or re-walk operation 133; small-batch inversion 132; solving 2^9;
verification 2^13.

**Budget.** Phase-1 operations: 1.064 (S_ops + 2^60) < 2^126.28889, plus the
surcharges (random, expectation < 2^106.417 + 2^99.819). Hard extras: meter
remainders < 2^61, phase 2 (2^20 (L_max+65) 133 = 2^95.377), re-walks
(8*2*((L_max+65)133 + 2^12) = 2^79.377), reports (D_cap 0.5 = 2^62.25), hash init
2^57.1, table 2^28, verification 2^13: together < 2^95.401. Random extras
(expected): starts 2^74.4, words/attempts/probes < 2^63, small batches < 2^8,
fallbacks 2^72.3: together < 2^74.8. So
W = 1.064(S_ops+2^60) + 2^95.401 + 2^106.432/0.00048 = 2^126.292044 (the third
term is the headroom) and W <= W_up =
104158187191707310815349162299993028690 with log2 W_up = 126.2920439. The claimed
score is that ceiling: 2^126293 >= W_up^1000 and 2^126292 < W_up^1000, both
checked with integers in `experiment.py`.

**Memory.** D_cap records of 104 bytes in 2^64.25 slots (load 1/2) = 2^70.95
bytes, the table 6*2^16 entries of 160 bytes = 2^25.91, walk
states 2^20*512 = 2^29: total 2^70.9504 < 2^72. **Preprocessing** (inside W):
table 2^28, report-table initialisation 2^57.1, first starts 2^32, sampling 2^21,
so expected < 2^57.11 against cap 2^57.2; the random part < 2^21 gives a halt
probability < 2^-32 (Markov).

## 6. Experiments

`experiment.py` (stdlib; JSON config on stdin) checks in Part A, on production
BLS12-381, the identities of sections 1-2 and the twist index against the exact
rule, and in Part B the walk on a toy curve over F_q (58493047 or 861101287). The declared run (q26, R = 64, 4 distinguished bits,
150 trials) reports pipeline success 0.46, fraction with T_c <= S_ops 0.4533 and
2-cycle rate 0.00245/step against 1/(6R) = 0.002604; its `cost` block prints every
number above and the exact integer check.
Offline, the same file at R = 2048 with 5 distinguished bits over seeds 6101-6103
(q26, 1500 trials each) and 6201-6202 (q30, 700 each) gives 0.4247, 0.4093, 0.4127
and 0.4086, 0.3871: pooled 0.4114 +- 0.0064 over 5900 trials against the H1 floor
0.3905, the mean T_c 0.8% above the model. README.md lists the checks and caveats.
Generic attack; no pairing, no NFS.
