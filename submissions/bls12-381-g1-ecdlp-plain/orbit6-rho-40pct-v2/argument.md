# bls12-381-g1-ecdlp-plain: parallel Pollard rho on size-6 automorphism classes

Claim: x is recovered from T=[x]P0 with probability at least 0.39 using at most
2^(633/5)=2^126.60 target-group-operation equivalents at the fixed C=2000.
Peak memory 2^36.5 bytes, preprocessing 2^29 units (included in total),
nonuniform advice 2^0 bytes. Classical; heuristics H1-H3 disclosed below.
Every number below was re-derived independently in plain integer/float Python
against the protected target-profile.json and ../classical-group-work-v2.json;
see certificate.json for machine-checked identities.

## Group structure (verified at production constants)

p = 0x1a0111ea...aaab (381 bits), p ≡ 1 mod 3. r = 0x73eda753...00000001,
254.85709 bits, prime (Miller-Rabin, bases 2..37). lambda = z^2-1 mod r with
z = 0xd201000000010000 satisfies lambda^2+lambda+1 ≡ 0 mod r, lambda ≠ 1,
so lambda has order 3. Both roots beta of x^2+x+1 ≡ 0 mod p were computed;
direct scalar multiplication on the standard generator P0 (blst constants)
shows (beta·x, y) = [lambda]P0 for one of the two roots (the other root
gives lambda^2). Naming either root beta, the six maps u ∈ {±1, ±phi, ±phi^2}
with phi(x,y) = (beta·x, y) act freely on G1∖{infinity} (r prime, r ∤ 6,
phi fixed-point-free mod r), so every class has exactly 6 points and there
are N' = (r-1)/6 classes, log2 N' = 252.27213. The sqrt(6) time gain is
log2(6)/2 = 1.29248 bits versus a negation-free walk. T = infinity returns 0
after validation, at negligible cost.

## Algorithm

1. Table: nT = 2^16 pairs (a_j,c_j) uniform mod r, points T_j=[a_j]P0+[c_j]T
   (≤4000 units each: 2^28 total), storing for each the six class-preimage
   images u^-1 T_j with scalar pairs u^-1(a_j,c_j).
2. m = 2^14 independent walkers from fresh random [a]P0+[c]T. Each holds the
   affine point R, accumulators (a,c) in unreduced 6-limb form, a step
   counter t mod alpha = 64 and one saved canonical x-coordinate S.
3. Step: from x compute x1 = beta·x (one Fp mul), x2 = -(x+x1) mod p (since
   beta^2 = -1-beta). Canonical class key x_C = min of the three x-values of
   {R, phi R, phi^2 R} with y < p/2; x_C identifies the class. Set
   h = low 16 bits of x_C, R <- R + u^-1 T_h, (a,c) <- (a,c)+u^-1(a_h,c_h).
   The walk on classes is a function of the class alone.
4. Distinguished point: bits 16..113 of x_C zero (theta = 2^-98). Report
   (x_C, u(a,c) mod r) to a hash store keyed on x_C. A repeated key with
   different (a,c), c ≠ c' mod r, yields x = (a-a')/(c'-c) mod r, verified
   by one scalar multiplication [x]P0 = T.
5. Fruitless cycles (Bos-Kaihara-Kleinjung-Lenstra-Montgomery 2010;
   Bernstein-Lange-Schwabe 2011): fixed-period-64 checkpoint detection on
   class keys only; on detection, walk the cycle once, remember the smallest
   class key C* (unique, a function of the cycle alone), restore the state
   there and replace R by 2R, (a,c) by (2a,2c). The escape class is
   cycle-determined, not walker-determined, so merged trails stay merged and
   van Oorschot-Wiener collision-search logic applies unchanged.

## Success probability

Budget S = 1.03·k0 + 24·m/theta steps, k0 = sqrt(2·ln(1/0.61)·N'),
log2 k0 = 126.12779 (re-derived: 2·ln(1/0.61) = 0.988593;
1.03^2·ln(1/0.61) = 0.524400). Under H1 the collision probability at
1.03·k0 visits is 1-exp(-0.5244) = 0.408089. A collided trail pair merges
and the next DP appears geometrically with mean 1/theta; the 24/theta
per-walker allowance bounds the miss probability at e^-24 = 3.8e-11 per
collision; incomplete-addition exceptions < 2^-120. Bound:
0.408089 - e^-24 - 2^-120 >= 0.39. The 3% length margin absorbs walk
non-randomness with 0.018 probability headroom; the floor is carried at
exactly 0.39 as the historical branch argued, and the 0.408 point-estimate
is disclosed. The delay allowance 24·m/theta = 2^116.585 is 0.130% of
1.03·k0 (verified).

## Per-step cost (charged in instructions / C=2000)

One batched incomplete affine addition = 1 reference unit (6M convention,
S=M, A=0, batch inversion share inside the 6M). Other work, itemised:

| Item | Instructions |
| --- | ---: |
| Fp mul beta·x (measured blst_fp_mul = 383, COST_MODEL) | 383 |
| x2 = -(x+x1) mod p (6-limb add, cond. subtract) | 45 |
| load x,y of R | 12 |
| min of three 6-limb values | 40 |
| y sign test vs (p+1)/2 | 14 |
| class index u and inverse index | 6 |
| DP test + table index extraction | 10 |
| table address | 6 |
| table loads (12) + scalar increments (8) | 20 |
| two 6-limb accumulations (load/add-carry/store) | 40 |
| cycle check (compare S, counter, cond. save) | 10 |
| loop/batch bookkeeping | 12 |
| itemised | 598 |
| reserve (batch inversion share, spills) | 52 |
| charged | 650 |

Cycle events: derived 2-cycle rate 1/(6·nT) = 2^-18.585 per step (verified);
L-cycle rate ≤ Bell(L)·nT^-ceil(L/2) (L=3: 2^-29.68, L=4: 2^-28.06, verified);
charged rate bound 4/nT = 2^-14 per step, a factor-24 margin over the derived
2-cycle rate. Per event ≤ 3·alpha = 191 steps × 1.325 units + ≤30 units
restore/doubling ≤ 284, charged 300. Per-step total = 1 + 650/2000 +
300·2^-14 = 1.343310546875 units exactly, log2 = 0.425793 (verified).

## Total

| Component | Units |
| --- | ---: |
| Walk steps S = 1.03·k0 + 24·2^14·2^98 | 2^126.17231 |
| S × 1.343310546875 | 2^126.59810 |
| Table build (2^16 × ≤4000 + six images) | < 2^28.3 |
| Walker starts + DP records (2^28.2 + 2^14 events × ≤4000) | < 2^40.3 |
| DP processing (≤5000 instr each) | < 2^29.5 |
| Solve, verification, constants | < 2^31 |
| All non-walk work | < 2^41 |

W = 2^126.59810 + 2^41 < 2^126.5982 (re-derived log2 W = 126.59810),
time_log2 = 126.60 rounded up, claimed exactly as the rational 633/5.
Preprocessing 2^29 (table + constant checks: verify beta, lambda, the
pairing identity on P0, ≤2^10 units). All cross-walker sums charged; record
transfer included in the 5000-instruction DP charge.

## Memory

DPs stored: S·theta = 2^28.17 records × 128 B = 2^35.17 B; hash table at
load factor 1/2: 2^36.2 B. Table 2^16×960 = 2^25.9 B; walkers 2^14×320 =
2^22.3 B; scratch < 2^20 B. Peak < 2^36.5 bytes. Memory unscored.

## Heuristics (disclosed)

- H1: the class walk behaves as a random function on the N' classes up to
  the charged 3% length margin (Teske 2001 for adding walks with ≥2^16
  tables; BKL LM 2010, BLS 2011 for class walks). Nonrandomness raises the
  collision rate, not lowers it. Not proved.
- H2: DP bits (16..113) and index bits (0..15) of x_C are independent and
  uniform.
- H3: cycle rates as above; L>64 cycles summed over L=65..200 contribute
  < 2^-302 per step, so < 2^-175 over the run. A 10× true rate would move
  the score to about 126.77; the 300-unit charge is worst-case.

## Limits

Analytic accounting; no full-size execution. The 383-instruction Fp
multiply is the organizer's pinned measurement; other itemised counts are
instruction-sequence estimates with a reserve. No pairing, index calculus
or sub-square-root method claimed. This is the textbook automorphism rho,
not a novelty claim.
