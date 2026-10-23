# bls12-381-g1-ecdlp-plain: parallel Pollard rho on size-6 classes

Claim: x is recovered from T=[x]P0 with probability >= 0.39 using at most 2^126.81 target-group-operation
equivalents (C=2000 instructions per reference operation). Peak memory 2^36.5 bytes, preprocessing 2^29 units
(included), advice 2^0. Classical, under the three heuristics H1-H3 of the section "What is assumed". Scores are not
comparable with other tracks. Every number below was recomputed with exact arithmetic; W = 2^126.80963 before rounding up.

## Group structure used

p = field modulus, r = prime subgroup order, log2 r = 254.857, G1=<P0>. Since p = 1 mod 3 there is beta in Fp of
order 3 (beta^2+beta+1=0) and phi(x,y)=(beta*x,y) maps G1 to itself with phi(Q)=[lambda]Q, lambda = z^2-1 with
z=0xd201000000010000 (lambda^2+lambda+1=0 mod r). The six maps u in {+-1,+-phi,+-phi^2} (the group mu_6, the units of
Z[omega]) act on G1 minus {infinity} (r prime and phi has no fixed point of order r), so every class has exactly 6
points and there are N'=(r-1)/6 classes, log2 N' = 252.2721. The pairing of beta with lambda is fixed by one scalar
multiplication on P0 (about 2^10 units). If T=infinity return 0 (checked before anything else).

## Algorithm

1. Table: nT=2^16 pairs (a_j,c_j) uniform mod r and points Tj=[a_j]P0+[c_j]T (two scalar multiplications each, at most
   4000 units per entry, 2^28 in total). For each j store the six images u^-1 Tj (x-coordinates beta^k*x, y or -y) and
   the scalar pairs u^-1(a_j,c_j) mod r. Because a_j is uniform, the Tj are independent uniform points of G1.
2. Walkers: m=2^14 independent walkers, each started at [a]P0+[c]T with fresh uniform (a,c) (at most 4000 units per start).
   A walker holds the actual affine point R, scalar accumulators (a,c) in 6 limbs (no reduction between distinguished
   points), a step counter t (mod alpha=64) and one saved canonical x-coordinate S (the checkpoint register).
3. Canonical class point. From (x,y) of R compute x1=beta*x (one Fp multiplication) and x2=-(x+x1) mod p (since
   beta^2=-1-beta). The canonical point C of the class has x_C = min(x,x1,x2) as 6-limb integers and the y that is < p/2;
   this fixes the unit u in mu_6 with u*R=C. x_C identifies the class (two points with the same x differ by sign).
4. Hash. Let H(x_C) be the 128-bit value (a,b) computed from the six little-endian 64-bit limbs l_0..l_5 of x_C by two
   lanes, a <- (a xor l_i)*M1, a <- a xor (a>>32), and b <- (b + l_i)*M2, b <- b xor (b>>29), for i=0..5 (start
   values the digits of pi, M1=0x9E3779B97F4A7C15, M2=0xBF58476D1CE4E5B9), then b <- b xor a and the splitmix64
   finaliser z <- (z xor z>>30)*M2, z <- (z xor z>>27)*M3, z <- z xor z>>31 on each lane (M3=0x94D049BB133111EB).
   The table index is h = a mod 2^16. The distinguished test is: bits 16..63 of a and bits 0..49 of b all zero
   (theta = 2^-98). Index and test read disjoint bit fields of one hash value. The hash is a function of the class only.
5. Step. R <- R + u^-1 T_h (the table image u^-1 T_h) and (a,c) <- (a,c)+u^-1(a_h,c_h). R is always the known combination
   aP0+cT, and the walk on classes, class(C) -> class(C+T_h), is a function of the class only.
6. Distinguished point (tested first). The walker reports (x_C, u*(a,c) mod r) to a hash store keyed by x_C and restarts
   at a fresh random point. If x_C is already stored with scalars (a',c') != (a,c) and c != c' mod r, output
   x = (a-a')/(c'-c) mod r and verify [x]P0=T by one scalar multiplication; otherwise continue.
7. Fruitless cycles (Bos-Kaihara-Kleinjung-Lenstra-Montgomery 2010, Bernstein-Lange-Schwabe 2011). The class walk is a
   function of the class, so a walker that reaches a class on a cycle of the function circles it forever. Detection at a
   fixed period alpha=64: after the distinguished-point test of step 6, (i) if S is set and x_C = S (one 64-bit
   compare, then a full 6-limb compare), a cycle was found, go to (iii); (ii) otherwise if t is a multiple of alpha,
   set S = x_C. Compare comes before save and S is compared at every step, so a cycle of length L <= alpha is detected at
   most (alpha-1)+L steps after the walker enters it, and L is exactly the number of steps since S was saved.
   (iii) Escape. The walker walks the cycle once, L further steps, remembering the member with the smallest x_C (unique,
   a function of the cycle alone) with the state (R,a,c) held there, restores that state and replaces R by 2R, (a,c)
   by (2a,2c), t=0, S cleared. Since R = u^-1 C* for the walker's own representative u, 2R = u^-1(2C*), so the class of
   the new R is the class of 2C* whichever representative or entry point the walker had: every walker on that cycle
   continues on the same class sequence, so merged trails stay merged. Cycles with L > alpha are not detected; the
   section "Fruitless cycles" bounds them for every length. The doubling runs in the shared batch like an addition. A
   cycle containing a distinguished class never traps a walker because step 6 comes first.
8. Budget: stop when the total number of class visits (steps) reaches S = k_B + 24*m/theta, with
   k_B = 1.15 * sqrt(1+1/nT) * k0 and k0 = sqrt(2*ln(1/0.61)*N'), log2 k0 = 126.127787, log2 k_B = 126.329432.
   Steps spent detecting and walking cycles are not counted in S; they are charged in the cost section.

## What is assumed

No theorem covers this exact walk. Three heuristics are used, each stated with its support and its charge.

H1 (hash behaves as a random oracle). On the classes the run visits (at most 2^127), the pair (table index, distinguished
flag) taken from H(x_C) is independent and uniform with the probabilities of step 4, independent across classes and of the
table points, the start points and the scalars. This replaces any claim about bits of the coordinate itself: the
canonical x is the minimum of three field elements and we do not claim anything about its bits. H is an ordinary
non-cryptographic multiply-xorshift mixer, and random-oracle behaviour of such mixers on structured inputs is the usual
assumption of Pollard rho implementations. Its instruction cost is charged below (52 itemised instructions for hashing, index
extraction and the distinguished test).

H2 (walk quality). Under H1 the class walk is an additive r-adding walk with a random partition into nT=2^16 table entries,
acting on the quotient set of G1 minus infinity by mu_6 (negation, and the endomorphism of order 3). Its collision rate
is modelled by a uniform random function on an effective class count N_eff = rho^2 (1+1/nT) N' with rho = 1.10. Here
(1+1/nT) is the standard additive-walk allowance, about sqrt(1+1/r) more steps than a random function for an r-entry
table (Teske 2001; Bernstein-Lange-Schwabe 2011), charged in the direction that costs us, and rho = 1.10 is an explicit
allowance for everything else, chosen from the evidence below. This is the standard heuristic, not a proof. Teske 2001
analyses and measures r-adding walks with r >= 16 in groups, not this quotient; walks with the negation map and with
automorphism groups of small order are treated with the same random-mapping model (Duursma-Gaudry-Morain 1999; Bos et al.
2010; Bernstein-Lange-Schwabe 2011, which also analyses the fruitless cycles the negation map creates). The toy
measurements below run this exact walk (same canonical form, same hash and bit fields, same cycle escape) as a consistency
check at toy scale; they are far too small to measure rho, which is taken from the standard model above.

H3 (cycle patterns). (a) The canonicalising units v_i met along a trail are independent uniform on mu_6 and independent of
the table indices; this is used only for lengths above 8192, where even "each vanishing sum has probability at most 1/2"
is enough. (b) If the formal combination of table points in a closing window does not vanish, the walk closes with
probability at most 6/r for each starting step and length: the closing equation is one linear equation in an independent
uniform table point whose coefficient is a nonzero Gaussian integer of norm below r.

## Success probability

A useful collision is two visits to the same class by different (a,c) lines. Then a + c x = u(a' + c' x) with u in mu_6,
solvable for x unless the scalars cancel, which has probability about 2^-253 per collision (starts are uniform). The
probability over fresh coins is the same for every instance T, so it holds over uniform challenges.

From P = 1 - prod_{i<k}(1 - i/N_eff) >= 1 - exp(-k(k-1)/(2 N_eff)). With k = k_B and N_eff = rho^2 (1+1/nT) N' the
exponent is (1.15/1.10)^2 * ln(1/0.61) = 1.09298 * 0.494296 = 0.54025, so the collision probability is at least
1 - exp(-0.54025) = 0.4174. (The factor 1+1/nT cancels between k_B and N_eff.) In words: the budget meets 0.39 for every
true effective class count up to 1.15^2 (1+1/nT) N', and the model charged 1.10^2 (1+1/nT) N'.

Losses subtracted from 0.4174:
- both trails of a collision must reach a distinguished point. After a collision the two trails coincide (the class
  sequence is a function of the class, escapes included, step 7) and the next distinguished point appears after a
  geometric number of steps of mean 1/theta; every walker is given 24/theta further steps, miss probability
  e^-24 = 3.8e-11 per collision. The 24*m/theta = 2^116.585 steps are inside S (0.12% of k_B);
- a walker trapped in an undetected cycle: at most 2^-27.9 over the whole run (section "Fruitless cycles");
- incomplete-addition exceptions (equal x coordinates): below 2^-120 over the run;
- cancelling scalars: below 2^-250.

So P >= 0.4174 - 2^-27.9 - 3.8e-11 - 2^-120 > 0.4173 >= 0.39. The margin is 0.027 in probability, or a factor
1.15 in steps against a charged 1.10. The run may stop early at the first useful collision, but the full budget is charged.

## Fruitless cycles

Write the walker's position R_i and R_{i+1} = R_i + v_i T_{h_i}, v_i in mu_6 (v_i = u^-1 of step 5, a product of the
canonicalising units, acting as a scalar). Take the first time t* a trail revisits one of its own classes:
C_{t*} = C_{t'}, L = t*-t'. The classes C_{t'},...,C_{t*-1} are distinct and so are their oracle values: the L table
indices h_{t'},...,h_{t*-1} are independent uniform on [nT] (H1, deferred decisions). Closure means R_{t'+L} = u R_{t'}
for some u in mu_6, that is Z = (u-1)R_{t'} with Z = sum_i v_i T_{h_i} = sum_j m_j T_j and m_j = sum over the window's
steps with index j of v_i, a Gaussian integer.

Closure is either formal (all m_j = 0, and u=1) or accidental (some m_j != 0). All lengths are covered:

1. Formal closures. m_j = 0 forces every used index to appear at least twice in the window, since a single appearance
   leaves m_j = v_i != 0, whatever the units are. So the probability is at most Q(L) = P(no index drawn exactly once among
   L independent uniform draws from nT) = sum_b A(L,b) (nT)_b / nT^L, A(L,b) the number of partitions of L steps into b
   blocks of size at least 2 (exact for any units). Values: Q(2)=2^-16, Q(3)=2^-32, Q(4)=2^-30.4, Q(5)=2^-44.7,
   Q(6)=2^-44.1, Q(8)=2^-57.3, Q(16)=2^-107, Q(64)=2^-365.6. If the units are uniform, a block of size s >= 2
   vanishes with probability at most 1/6, giving the derived 2-cycle rate 1/(6 nT) = 2^-18.58 (the canonical form of
   C+T_h is -(C+T_h) with probability 1/6 and the table index is again h with probability 1/nT) and 2^-18.58 for all
   lengths up to 64 together.
2. Detected lengths L <= 64. The rate of cycle entries per step is at most sum_{L=2}^{64} Q(L) + 64*6/r = 2^-15.99991.
   The charge below uses 4/nT = 2^-14, a factor 3.9998 above this bound (and 2^4.4 above the derived 2^-18.58).
3. Undetected lengths L >= 65 (a walker on such a cycle would circle until the end of the run).
   (i) 65 <= L <= 8192. Q(L) <= sum_{b<=L/2} C(L,b) (b/nT)^{L-b} <= (2L/nT)^{L/2}, because a partition into b blocks
   is fixed by b representatives and an assignment of the other L-b steps, has probability nT^{-(L-b)}, and b <= L/2,
   L-b >= L/2. For L <= nT/8 the bound decreases in L, from 2^-291.8 at L=65, and sums to below 2^-291.6.
   (ii) L > 8192 (closure term, formal part). A formal closure needs D >= 1 distinct indices, each used at least twice,
   with every block sum zero. P(D <= 4095) <= C(nT,4095) (4095/nT)^L <= 2^(22110-4L) <= 2^-10600 for L >= 8192, and
   with at least 4096 distinct indices H3(a) gives at most 2^-4096. A closure at lag L also needs the trail to survive L
   steps without a distinguished point, probability at most e^{-theta L}, and sum_L e^{-theta L} <= 1/theta = 2^98. This
   part is below 2^-3990 per step for all L > 8192 together.
   (iii) Accidental closures for every L >= 65 (closure term, accidental part). Probability at most 6/r per starting step and
   length (H3(b); the Gaussian coefficient has norm at most L^2 < r for every L < 2^127, and lags beyond have survival weight
   below e^(-2^29)), times the survival e^{-theta L}, summed over all L: at most (6/r)/theta = 2^-154.27 per step. This is
   the generic birthday self-collision of a trail of mean length 1/theta.
   The rate of undetected-cycle entries is at most 2^-291.6 + 2^-3990 + 2^-154.27 < 2^-154.27 * (1 + 2^-130) per step.
   Over the S = 2^126.33 steps of the run the expected number of trapped walkers is at most 2^-27.9. A trapped walker
   would waste the rest of its run, so this probability is subtracted from success (above) and no work is credited.

## Cost, per walk step

One batched incomplete affine addition = 1 reference unit (6M with S=M, A=0, including its share of batch
inversion). All other work is charged in AArch64 instructions / 2000:

| Item | Instructions |
| --- | ---: |
| Fp multiplication beta*x (measured blst_fp_mul) | 383 |
| x2 = -(x+x1) mod p: 6-limb add, conditional subtract p, subtract from p | 45 |
| load x and y of R | 12 |
| min of three 6-limb values (two compare chains, selects of index and low limbs) | 40 |
| y sign test against (p+1)/2 | 14 |
| class index u from the two decisions, inverse index | 6 |
| hash lanes: 6 limbs x 2 lanes x 3 (xor or add, multiply, xor with shifted operand) | 36 |
| cross mix (1) and two finalisers (2 x 5) | 11 |
| table index extraction (1) and distinguished test (4) | 5 |
| table address | 6 |
| loads of table x, y (12) and scalar increments (8) | 20 |
| two 6-limb scalar accumulations (load, 6 adds with carry, store) | 40 |
| Cycle check: compare with S, counter, conditional save of S | 10 |
| loop and batch bookkeeping | 12 |
| Itemised total | 640 |
| Reserve (batch inversion share about 12, hash constants and moves, spills, branch fixes, 17% margin) | 110 |
| Charged | 750 |

Per-step charge 1 + 750/2000 = 1.375 units. Cycle escapes: rate at most 2^-14 events per step (section "Fruitless
cycles", item 2). Cost per event: from entering a cycle until detection at most (alpha-1)+L <= 127 steps, then L <= 64
steps of walking the cycle once, so at most 191 steps at 1.375 units = 262.6 units, plus the stored-state restore and
6-limb compares (at most 10 units) and one batched doubling (at most 20 units): at most 292.6 units, charged 300.
Charge per step 300*2^-14 = 0.018310546875 units. Per step total 1.393310546875 units, log2 = 0.478517.

## Total

| Component | Units |
| --- | ---: |
| Walk steps: S = k_B + 24*2^14*2^98 | 2^126.33111 |
| Steps times 1.393310546875 (including hash and cycle escapes) | 2^126.80963 |
| Table construction (2^16 entries at <= 4000, six images) | < 2^28.3 |
| Walker starts (distinguished points 2^28.33 plus 2^14) at <= 4000 | < 2^40.3 |
| Distinguished-point processing: 5000 instructions = 2.5 units each | < 2^29.7 |
| Final solve and verification, constants | < 2^31 |
| All non-walk work | < 2^41 |

W = 2^126.80963 + 2^41 < 2^126.8097, so time_log2 = 126.81 (rounded up from 126.8097). Compared with a 1.03 k0 budget at 1.3433 units per step this adds 0.21 to the exponent: 0.16 for the 1.15 allowance and the penalty, 0.05 for the hash. Preprocessing is the table
(2^28.3) plus constant checks, stated as 2^29 units. Transfers of the 2^28.33 records of 128 bytes are counted in the
5000 instructions each. A tenfold error in the cycle rate would move the score by about 0.16.

## Memory

Distinguished points expected: S*theta = 2^28.33; each record is x_C (48 bytes), two scalars (64) and padding to 128 bytes,
so 2^35.33 bytes, 2^36.33 with a hash table at load factor 1/2. Table 2^16*960 bytes = 2^25.9, walkers 2^14*320 bytes = 2^22.3,
code and scratch under 2^20. Peak under 2^36.5 bytes (stated; memory carries no score). With theta=2^-98 the store is large
by design: it keeps the detection delay at 0.12% of the work.

## Toy-scale experiment (supporting evidence only)

experiment.py runs this exact walk (steps 3-7, with the table index taken from the low bits of the same 128-bit hash and
the distinguished test from bits 16 and up, disjoint from the index as in the production walk) on three toy groups
y^2 = x^3 + b over F_p with p = 1 mod 3 and prime order n = 1 mod 3: (p,b,n) = (2251,7,2341), (30091,2,29881),
(130003,3,130261), with 390, 4980 and 21710 classes. Four configurations are run, with 4 or 8 walkers, tables of 16 or 64
entries and the checkpoint period at the production value 64. A control with the same walkers and the same collision rule
moves to a uniformly random class at every step. For each configuration the script prints the mean number of class visits
to the first collision between two walkers, that mean as a ratio to sqrt(pi N'/2), the delta-method 95% interval of the
ratio, the ratio to the uniform control with its 95% interval, the number of cycle events, the longest cycle and the escape
steps per visit. It then recovers x through distinguished points on the middle group and checks the two hash fields over
every class of the largest group. There is no success-at-budget column in its output.

At the default input trials=20 per configuration (four trials on each of five fixed seeds, under a second):

| classes | table | walkers | mean visits | ratio to sqrt(pi N'/2) (95%) | ratio to control (95%) |
| ---: | ---: | ---: | ---: | --- | --- |
| 390 | 16 | 4 | 30.35 | 1.2262 (0.9319 to 1.5206) | 1.1036 (0.7407 to 1.4665) |
| 4980 | 16 | 8 | 94.0 | 1.0628 (0.691 to 1.4346) | 1.1975 (0.7055 to 1.6894) |
| 4980 | 64 | 8 | 98.95 | 1.1188 (0.9192 to 1.3183) | 1.24 (0.8582 to 1.6217) |
| 21710 | 64 | 8 | 222.05 | 1.2024 (0.9007 to 1.5041) | 1.149 (0.7575 to 1.5405) |

The mean of the four ratios to sqrt(pi N'/2) is 1.1525 and every interval covers 1.0. At 20 trials the intervals are 0.40
to 0.74 wide, so this run cannot resolve the 10% allowance in N_eff: it is a consistency check of the exact walk, not a
measurement of rho. A run at trials=1000 (1.7 s) narrows the intervals to 0.07 to 0.17 and moves the four ratios to
1.3794, 1.1907, 1.1094 and 1.0948, while the ratio to the control falls from 1.1866 and 1.1127 at 16 entries to 1.0388
and 1.0203 at 64. The excess over the control therefore falls as the table grows, the direction the additive-walk penalty
predicts, and the charged table penalties at those sizes are 1.031 and 1.008. The toy tables are 4096 and 1024 times
smaller than the production table of 2^16 entries and the largest toy group has 2^14.4 classes against 2^252.3, so we take
rho = 1.10 and charge 1.15; the measurements do not prove the production ratio.

Cycle escapes were exercised: 7 detections over the four configurations, longest cycle 16, at most 0.056 escape steps per
visit, so the escape rule is part of what runs. The recovery part solves 6 random targets through distinguished points
(64-entry table, 8 walkers) and recovers x in all 6, with mean 143.0 visits against sqrt(4980) = 70.6. Over all 21710
classes of the largest group the table index gives chi-square 16.65 on 15 degrees of freedom, the distinguished rate is
0.0614 against 1/16 = 0.0625, and the joint index-times-test chi-square is 13.33 on 15 degrees of freedom.

## Limits

No full-size execution. The 383-instruction Fp multiplication is the organiser's fixed measurement; the other per-step
items are counted as instruction sequences and padded with a reserve. Memory latency is not an instruction and carries no
charge in this model. The bound is the textbook automorphism rho; it does not use the pairing, index calculus or any
algorithm faster than the square-root class search, and does not claim optimality. If the true walk factor exceeded 1.15,
the budget would have to grow by that excess, about 0.014 in the exponent per 1% above. If the hash were not oracle-like on
the visited classes, the walk analysis would not apply; this is H1 and is not proved.

References: van Oorschot and Wiener, Parallel collision search with cryptanalytic applications, J. Cryptology 1999; Teske, On random
walks for Pollard's rho method, Math. Comp. 2001; Duursma, Gaudry, Morain, Speeding up the discrete log computation on curves
with automorphisms, Asiacrypt 1999; Bos, Kaihara, Kleinjung, Lenstra, Montgomery, On the security of 1024-bit RSA
and 160-bit elliptic curve cryptography (2010); Bernstein, Lange, Schwabe, On the correct use of the negation map in the Pollard rho
method, PKC 2011; Gallant, Lambert, Vanstone, Faster point multiplication on elliptic curves with efficient endomorphisms, 2001.
