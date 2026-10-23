# bls12-381-g1-ecdlp-plain: parallel Pollard rho on size-6 classes

Claim: x is recovered from T=[x]P0 with probability >= 0.39 using at most 2^126.60 target-group-operation
equivalents (C=2000 instructions per reference operation). Peak memory 2^36.5 bytes, preprocessing 2^29 units
(included), advice 2^0. Classical, with the heuristics H1-H3 below. Scores are not comparable with other tracks.

## Group structure used

p = field modulus, r = prime subgroup order, log2 r = 254.857, G1=<P0>. Since p = 1 mod 3 there is beta in Fp of
order 3 (beta^2+beta+1=0) and phi(x,y)=(beta*x,y) maps G1 to itself with phi(Q)=[lambda]Q, lambda = z^2-1 with
z=0xd201000000010000 (lambda^2+lambda+1=0 mod r). The six maps u in {+-1,+-phi,+-phi^2} act on G1 minus {infinity}
(r prime and phi has no fixed point of order r), so every class has exactly 6 points and there are N'=(r-1)/6 classes,
log2 N' = 252.2721. The pairing of beta with lambda is fixed by one scalar multiplication on P0 (about 2^10 units).
If T=infinity return 0 (checked before anything else).

## Algorithm

1. Table: nT=2^16 pairs (a_j,c_j) uniform mod r and points Tj=[a_j]P0+[c_j]T (two scalar multiplications each, at most
   4000 units per entry, 2^28 in total). For each j store the six images u^-1 Tj (x-coordinates beta^k*x, y or -y) and
   the scalar pairs u^-1(a_j,c_j) mod r.
2. Walkers: m=2^14 independent walkers, each started at [a]P0+[c]T with fresh uniform (a,c) (at most 4000 units per start).
   A walker holds the actual affine point R, scalar accumulators (a,c) in 6 limbs (no reduction between distinguished
   points), a step counter t (mod alpha=64) and one saved canonical x-coordinate S (the checkpoint register).
3. One step. From (x,y) of R compute x1=beta*x (one Fp multiplication) and x2=-(x+x1) mod p (since beta^2=-1-beta).
   The canonical class point C has x_C = min(x,x1,x2) as 6-limb integers and the y that is < p/2 (Montgomery form is a
   bijection, so any fixed total order works); this fixes u in the group of 6 with u*R=C. Because of the cube-root
   structure, x_C identifies the class. The step is R <- R + u^-1 T_h with h = low 16 bits of x_C, i.e. the class of
   C+T_h, using the table image u^-1 Tj. The scalars update as (a,c) <- (a,c)+u^-1(a_h,c_h). The point R is therefore always
   a known combination aP0+cT, and the walk on classes is a function of the class only.
4. Distinguished point: bits 16..113 of x_C all zero (theta=2^-98). The walker then reports (x_C, u*(a,c) mod r) to a hash
   store keyed by x_C and restarts at a fresh random point. If x_C is already stored with scalars (a',c') != (a,c) and
   c != c' mod r, output x = (a-a')/(c'-c) mod r and verify [x]P0=T by one scalar multiplication; otherwise continue.
5. Fruitless cycles (Bos-Kaihara-Kleinjung-Lenstra-Montgomery 2010, Bernstein-Lange-Schwabe 2011). The class walk is a
   function of the class alone, so a walker that reaches a class lying on a cycle of the function circles it forever.
   Detection at a fixed period alpha=64: the walker counts steps t since its last start or escape. After computing the
   canonical x_C of each visited class (and after the distinguished-point test of step 4, which comes first): (i) if
   S is set and x_C = S (one 64-bit compare, then a full 6-limb compare), a cycle was found, go to (iii); (ii) otherwise if
   t is a multiple of alpha, set S = x_C. Compare comes before save, S is compared at every step, so S is a window of
   alpha steps. (iii) Escape. The current class lies on a cycle, and S is the first class that returned, so the cycle
   length L is exactly the number of steps since S was saved and L <= alpha. The walker walks the cycle once, L further
   steps, remembering the member with the smallest x_C (x_C identifies the class, so the members have distinct x_C and
   the minimum is unique; it is a function of the cycle alone) together with the state (R,a,c) the walker had there. It
   then restores that state, so R lies in the class of C*, the minimum member, and replaces R by 2R, (a,c) by
   (2a,2c), and sets t=0 and clears S. Since R = u^-1 C* for the walker's own representative u, 2R = u^-1(2C*), so the
   class of the new R is the class of 2C* whichever representative the walker held. The class sequence after an escape is
   therefore the same for every walker that was on the cycle, whatever its checkpoint phase or its point of entry into
   the cycle: the rule keeps merged trails merged (step "Success probability"). Cycles with L > alpha are not detected;
   their rate is negligible (H3). The doubling is performed in the shared batch inversion like any addition. If any member of a cycle were
   distinguished, the walker would already have reported it and restarted on its first visit (step 4 is first).
6. Budget: stop when the total steps reach S = 1.03*k0 + 24*m/theta (worst case charged), with
   k0 = sqrt(2*ln(1/0.61)*N') = sqrt(0.98859*N'), log2 k0 = 126.12779. Steps spent detecting and walking cycles are
   not counted in S; they are charged separately in the cost section (3*alpha steps per event).

## Success probability

A useful collision is two distinct visits to the same class by different (a,c) lines: then a + c x = u(a' + c' x) with
u in the six maps, solvable for x unless the scalars cancel, which has probability about 2^-253 per collision (starts are
uniform). Under H1 the N' classes behave as random: with k class visits the probability of at least one collision is
1-exp(-k^2/(2N')). At k=k0 this is 1-0.61=0.39; the budget uses k=1.03*k0 steps, giving 1-exp(-0.5244)=0.408, so the
extra 3% of walk length covers walk non-randomness (H1) with room 0.018 in probability. Both colliding trails coincide after
the collision (the class sequence is a function of the class, so one walker is the other delayed by a fixed number of steps)
and the next distinguished point appears after a geometric number of steps of mean 1/theta; the allowance
24/theta per walker gives miss probability e^-24 = 4e-11 per collision. Fruitless cycles do not break this: if the trail
reaches a cycle, every walker on it detects the cycle within 2*alpha steps of entering it, walks it once and escapes to the
class of 2*C* with C* the cycle's smallest member (step 5), a class that depends on the cycle only; hence the two walkers
are again on one common trail after their escapes and the argument above applies from the escape class on. During the
detection the lagging walker only repeats classes the leading walker visited, none of which is distinguished. Between
distinguished points (about 2^98 steps) the trail passes through many cycle events; this is harmless because every escape
is determined by the class sequence and not by the walker (the rule does not look at where the walker's own checkpoint
fell, nor at when it started). Miss probability per collision is therefore unchanged, e^-24. Incomplete-addition exceptions (equal x
coordinates) have probability below 2^-120 over the whole run. The bound used is 0.408 - 4e-11 - 2^-120 >= 0.39.
The expected run is shorter (the search stops at the first useful collision), but the cost charged is the full budget.

## Cost, per walk step

One batched incomplete affine addition = 1 reference unit (6M with S=M, A=0, including its share of batch
inversion; walkers are processed in batches of 2^14, so the single Fp inversion per batch adds about 12
instructions per addition, which is covered by the reserve below). All other work is charged in AArch64
instructions / 2000:

| Item | Instructions |
| --- | ---: |
| Fp multiplication beta*x (measured blst_fp_mul) | 383 |
| x2 = -(x+x1) mod p: 6-limb add, conditional subtract p, subtract from p | 45 |
| load x and y of R | 12 |
| min of three 6-limb values (two compare chains, selects of index and low limbs) | 40 |
| y sign test against (p+1)/2 | 14 |
| class index u from the two decisions, inverse index | 6 |
| distinguished-point test and table index extraction | 10 |
| table address | 6 |
| loads of table x, y (12) and scalar increments (8) | 20 |
| two 6-limb scalar accumulations (load, 6 adds with carry, store) | 40 |
| Cycle check: compare with S, counter, conditional save of S | 10 |
| loop and batch bookkeeping | 12 |
| Itemised total | 598 |
| Reserve (batch inversion share, spills, branch fixes) | 52 |
| Charged | 650 |

So the per-step charge is 1 + 650/2000 = 1.325 units. Cycle escapes. Rate: a class C lies on a 2-cycle when the
canonical representative of C+T_h is -(C+T_h) (the automorphism that maps it to canonical form is the negation, probability
1/6 among the six) and the table index of that class is again h (probability 1/nT, since the class has the same x-class
and the index is the low 16 bits of the new x_C): then the walker adds T_h to the negated point and lands on -C, the
class of C. Every other closing pattern needs a vanishing relation between the independent random table points beyond
this one (probability about 2^-254 per class) or more table coincidences. Hence the derived 2-cycle rate is
1/(6 nT) = 2^-18.58 per step. A cycle of length L closes only if the L table points used cancel formally in Z[omega]T_j, which forces
at least ceil(L/2) coincidences among the L table indices, each of probability 1/nT, over at most Bell(L) coincidence patterns:
the rate of L-cycles is at most Bell(L)*nT^-ceil(L/2) (L=3: 2^-29.7, L=4: 2^-28.1, L=64: 2^-295). The charge uses
the rate bound 4/nT = 2^-14 events per step for all cycles together, a margin of 24 over the derived 2-cycle rate.
Cost per event: from the walker's entry into a cycle until detection at most (alpha-1)+L steps (wait for the next
checkpoint, then at most L steps to return to it), then L steps of walking the cycle once, so at most 3*alpha-1 = 191 steps
at 1.325 units = 253.1 units, plus the stored-state restore and 6-limb compares (at most 10 units) and one batched doubling (at
most 20 units): at most 284 units, charged 300 units per event. Charge per step 300*2^-14 = 0.01831 units. Per step total:
1.343310546875 units, log2 = 0.425793.

## Total

| Component | Units |
| --- | ---: |
| Walk steps: S = 1.03*k0 + 24*2^14*2^98 | 2^126.17231 |
| Steps times 1.343310546875 (including cycle escapes) | 2^126.59810 |
| Table construction (2^16 entries at <= 4000, six images) | < 2^28.3 |
| Walker starts (distinguished points 2^28.2 plus 2^14) at <= 4000 | < 2^40.3 |
| Distinguished-point processing: 5000 instructions each (normalise scalars mod r, hash insert) | < 2^29.5 |
| Final solve and verification, constants | < 2^31 |
| All non-walk work | < 2^41 |

W = 2^126.59810 + 2^41 < 2^126.5982, so time_log2 = 126.60 (rounded up from 126.5982). The delay allowance
24*m/theta is 2^116.585, which is 0.13% of 1.03*k0. Preprocessing is the table (2^28.3) plus constant checks, stated as
2^29 units. Sums over the 2^14 walkers are included; no coordination work beyond transferring 2^28.2 records of 128 bytes
(counted in the 5000 instructions each).

## Memory

Distinguished points expected: S*theta = 2^28.17; each record is x_C (48 bytes), two scalars (64) and padding to 128 bytes,
so 2^35.17 bytes, 2^36.2 with a hash table at load factor 1/2. Table 2^16*960 bytes = 2^25.9, walkers 2^14*320 bytes = 2^22.3,
code and scratch under 2^20. Peak under 2^36.5 bytes (stated; memory carries no score). With theta=2^-98 the store is large by
design: it keeps the detection delay at 0.13% of the work.

## Heuristics

- H1. The class walk is a random function on N' classes; Teske (2001) reports that r-adding walks with table size at least 16-20
  match random-function collision counts within a few percent, and nT=2^16 here. The 3% margin is charged. Collision counts
  for walks on classes under a group of size 6 are standard (Bos et al. 2010; Bernstein-Lange-Schwabe 2011); nonrandomness
  near-collisions only raise the collision rate. Not proved; supported by that literature and the toy run in experiment.py.
- H2. Hash bits for table index (bits 0..15) and distinguished-point test (bits 16..113) of x_C in Montgomery form are
  independent and uniform; choosing the minimum of three x values changes only the high bits.
- H3. Fruitless cycles start with probability at most 2^-14 per step (derived 2-cycle rate 2^-18.58; L-cycle rate at most
  Bell(L)*nT^-ceil(L/2), which assumes the 2^16 table points are independent random points), and cycles with L > 64 do not
  occur in the run: their summed rate over L=65..200 is below 2^-302 per step, so their probability over 2^126.2 steps is
  below 2^-175. If the true rate were 10 times the charged 2^-14 the score would move by about 0.17 (to 126.77); the 300-unit
  charge per event is a worst-case count. A cycle that is not detected would trap a walker until the end of the run;
  the bound above excludes this.

## Limits

No full-size execution. The 383-instruction Fp multiplication is the organiser's fixed measurement; the other per-step
items are counted as instruction sequences and padded with a reserve. Memory latency is not an instruction and carries no
charge in this model. The bound is the textbook automorphism rho; it does not use the pairing, index calculus or any
algorithm faster than the square-root class search, and does not claim optimality.

References: van Oorschot and Wiener, Parallel collision search with cryptanalytic applications, J. Cryptology 1999; Teske, On random
walks for Pollard's rho method, Math. Comp. 2001; Bos, Kaihara, Kleinjung, Lenstra, Montgomery, On the security of 1024-bit RSA
and 160-bit elliptic curve cryptography (2010); Bernstein, Lange, Schwabe, On the correct use of the negation map in the Pollard rho
method, PKC 2011; Gallant, Lambert, Vanstone, Faster point multiplication on elliptic curves with efficient endomorphisms, 2001.
