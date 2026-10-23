# Six-class parallel Pollard rho on BLS12-381 G1 (negation + GLV)

**Claim.** time_log2 = 126.41: the algorithm recovers x from T = [x]P0 with
probability >= 0.40, and on every run its total work is <= W_max = 2^126.41
units (one unit = one batched affine G1 addition = C = 2,000 scalar AArch64
instructions). The work bound is unconditional (global work meter, section 3).
Memory < 2^72 bytes; preprocessing (everything before the first walk step:
jump table, report-table initialisation, the first 2^20 trail starts) <= 2^57.2
units, enforced by the meter and included in W; no advice. Success
rests on the random-mapping heuristic H1 (section 4), supported by toy runs.

## 1. Algebra (all checked on BLS12-381 in `experiment.py`, part A)

p = field prime, r = subgroup order (log2 r = 254.85709), E: y^2 = x^3 + 4.

* beta = 0x1a0111ea...00000000aaac (primitive cube root of 1 in Fp),
  phi(x, y) = (beta x, y) = [lambda] on G1 with
  lambda = 0xac45a4010001a40200000000ffffffff (phi(P0) = [lambda]P0 checked).
* Units u_t = (-1)^s lambda^j (t = 3s + j) act by [u_t](x, y) = (beta^j x, (-1)^s y).
  No nonidentity G1 point is fixed by a unit != 1, so G1 \ {O} has exactly
  N = (r-1)/6 classes of size 6 (log2 N = 252.27213).
* **Key** K = |y| = min(y, p - y) is a complete class invariant: equal |y|
  gives y' = +-y and x'^3 = x^3, so x' = beta^j x.
* **Canonical member**: j minimises beta^j x mod p, s = [p - y < y];
  c = (beta^j x, |y|) = [u_t]P. The exact rule costs one Fp multiplication
  (beta x; beta^2 x = p - ((x + beta x) mod p)); the walk needs only the index
  t, which section 1.1 decides exactly without that multiplication.

### 1.1 Exact canonical index from a truncated product

G = floor(2^448 beta/p), PH = floor(p/2^256), PTOP = floor(p/2^320),
x = sum x_i 2^(64i) (6 limbs), v1 = beta x mod p, v2 = beta^2 x mod p.

1. L6' = bits 384..447 of the partial product of x and G formed from the high
   halves of the 5 limb products x_i G_j with i+j = 4, both halves of the 6
   with i+j = 5 and the low halves of the 6 with i+j = 6 (23 mul/umulh). The
   omitted part is < 2^324, so the true bits are L6' or L6'+1, and
   frac(xG/2^448) = v1/p - e with 0 <= e < 2^-67. Fall back if L6' < 2 or
   L6' > 2^64 - 3; otherwise v1/p is in [L6'/2^64, (L6'+2)/2^64 + 2^-67).
2. V1 = floor(PH L6'/2^128) has 0 <= v1/2^320 - V1 < 1.3; X = x_5 is the floor
   of x/2^320; S = X + V1 is within 2.3 of (x + v1)/2^320. Fall back if
   |S - PTOP| <= 8; otherwise V2 = PTOP - S (or 2 PTOP - S when x + v1 > p) is
   within 3.3 of v2/2^320.
3. jm = argmin(X, V1, V2); fall back if another value is within 8 of the
   minimum. Otherwise errors < 3.3 cannot reverse an approximate gap > 8, so jm
   is the exact argmin.

The fallback is the exact rule, a separately metered action (0.25 units). Its
windows total about 94 units out of > 2^60.7, so assuming beta x mod p is
equidistributed along the walk its rate is < 2^-54 per step (this affects only
the meter's early-halt probability). The decision always equals the exact rule.
Part A emulates the limb algorithm against the exact rule on 30,000 random and
10,000 adversarial near-tie production inputs (0 mismatches).

## 2. Walk

Parameters: R = 2048; distinguished probability theta = 2^-64 (limb 1 of K is
zero); M = 2^20 walks; L_max = 20 * 2^64; checkpoint period 64.

* **Table.** c_h, d_h uniform in [0, r) (fresh coins, rejection sampling),
  M_h = [c_h]P0 + [d_h]T (redraw if O), Tab[h][t] = [u_t^{-1}]M_h
  = (-1)^s phi^{-j}(M_h) with coefficients u_t^{-1}(c_h, d_h) mod r.
* **Step.** F(C) = class(c(C) + M_h), h = K mod R. The walk keeps any member P
  and sets P <- P + Tab[h][t]; since c + M_h = [u_t](P + Tab[h][t]) this is
  F. Each ordinary step is exactly one batched affine addition (reference
  operation) plus at most 240 instructions (section 5).
* **Trails.** Start S = [a]P0 + [b]T, a, b uniform from fresh coins, stored
  in full (redraw if S = O). At each class arrival, in order:

```
if K distinguished: report (K, a, b, length); end trail
if length >= L_max: abandon trail
if pending_double: P <- [2]P; clear history; continue
if fp(K) == K_{-2}: 2-cycle {current, previous}:
     if K < K_{-1}: P <- [2]P; clear history; continue
     else pending_double <- true
elif fp(K) == K_cp: traverse the cycle once more (<= 64 steps, else abandon),
     double the member of its minimum-key class; clear history; continue
if length mod 64 == 0: K_cp <- fp(K)
P <- P + Tab[h][t]; K_{-2} <- K_{-1}; K_{-1} <- fp(K)
```

  fp = limbs 0, 2, 3 of K (192 bits); `length` counts every operation, so a
  trail has <= L_max + 65 operations. Doubling commutes with units, so every
  escape is a class map determined by the cycle alone: trails that meet a
  common class stay merged. If an addition has equal x-coordinates
  (P = +-Tab), the trail is abandoned (probability < 2^-114 in total).
* **Solving.** On a key match with another trail's report, re-walk both from
  their stored (a, b), tracking P = [A]P0 + [B]T with every value reduced
  mod r after every operation (doubling doubles A, B). Canonical coefficients
  (e, f) = u_t (A, B). If the canonical points agree and f_1 != f_2:
  x = (e_2 - e_1)/(f_1 - f_2) mod r (sign checked on production data); verify
  [x]P0 = T and output x. T = O outputs 0.

## 3. Phases and the global work meter

* **Meter.** Every action (random words, sampling attempts, starts and redraws,
  steps, doublings, traversals, reports, hash probes, table build and hash
  initialisation, re-walks, solving, verification) is charged its price
  (section 5) before it runs; each price upper-bounds its instructions / C.
  If a charge would exceed W_max = 2^126.41 the attack halts with failure.
  So W <= 2^126.41 on every run; no loop needs its own limit. Preprocessing
  (all charges before the first walk step) has its own meter cap 2^57.2: if a
  preprocessing charge would exceed it, the attack halts with failure. The meter is a
  192-bit fixed-point counter; processors reserve blocks of 2^40 units by an
  atomic update (unused remainder counts as spent).
* **Phase 1** runs until the step charges reach U1 = 1.12 S_ops (within
  M * 2^40 = 2^60 because of blocks). **Phase 2** finishes open trails, no new
  trails. Report table capacity D_cap = 2^63.25 and at most 8 processed
  matches (overflow = failure); these bound memory and re-walk work.
* **Batches.** Ready walks (next action a step or doubling) are split into
  batches of 2^15..2^16 per Montgomery inversion; a smaller batch (fewer than
  2^15 ready walks) pays its inversion as a separate 132-unit action.

## 4. Success probability

N' = N * R/(R-1). T_c = fresh class visits (all trails) before the first
visit of a class already visited by another trail.

**H1 (heuristic, at the budget).** P(T_c <= S) >= 1 - exp(-(S/kappa)^2/(2N'))
for S = S_fresh, with kappa = 1.05 (a 5% slowdown allowance: toy class walks
are 0.5-2% slower than the random-mapping model, while a no-automorphism
control matches it; section 6). Also: partition and distinguished bits of |y|
are pseudo-random functions of the class.

S0 = sqrt(2 N' ln(1/0.6)) = 2^126.15187; S_fresh = 1.06 S0 = 2^126.23593
(H1 gives >= 0.4058; pure model 0.4367); S_ops = S_fresh/0.998 = 2^126.23882;
U1 = 1.12 S_ops = 2^126.40232.

**Fresh visits.** 2-cycles occur at rate 1/(6R) per step (need h' = h and
twist -1) and waste <= 2 visits; other fruitless cycles of length L <= 64 need
every partition index to occur >= 2 times (this includes identities such as
one index thrice with twists 1, lambda, lambda^2), so their rate is at most
sum_L sum_b S2(L, b) R^-(L-b) = 2^-19.995 (S2 = partitions into blocks of
size >= 2, summed exactly for L <= 64), each wasting <= 256 visits. Expected
waste < 0.00045 of operations; Bernstein's inequality (per-step contributions
<= 256, independent under H1) gives waste < 0.0019 and doublings
< 2^-13 S_ops except with probability < exp(-2^100). Then phase 1 performs
>= (U1 - 2^60)(1 - 0.88 * 2^-13/1.12)/1.12 >= 0.9999 S_ops operations and
>= 0.9999 * 0.9981 S_ops >= S_fresh fresh visits.

**Failure events** (given T_c <= S_fresh the merged trails reach the same
report and are solved unless): trail cap hit first (Gamma(2) tail, < 1e-7);
first collision inside one trail (L_max^2/(2N) per trail, < 2^-53 total);
report overflow (mean 2^62.24 reports, D_cap = 2^63.25: < exp(-2^60));
f_1 = f_2 (b uniform given the start point: 1/r); exceptional additions
(< 2^-114); more than 8 matches (< 2^-60); spurious 192-bit fingerprint match
(union bound over < 2^127.5 fingerprint comparisons: < 2^-64); meter halting early (< 2^-44, section 5); preprocessing cap
(< 2^-32, section 5); concentration
(< 2 exp(-2^100)). **Success >= 0.4058 - 2e-7 > 0.40.**

## 5. Cost

**Per ordinary step** (instructions besides the reference addition; blst
Fp multiplication = 383):

| Work | Instr. |
| --- | ---: |
| 1.1 step 1: 23 mul/umulh, 35 add/adc, 22 loads | 80 |
| 1.1 step 2: range test, V1 (3 mul, 2 add), S, wrap test, V2 | 20 |
| 1.1 step 3: differences, tolerance tests, argmin, fallback branch | 14 |
| sign, p - y, select limbs of K | 17 |
| h, distinguished test, twist index, table address | 7 |
| loads of table point and state, store of state | 36 |
| fingerprint tests, history and checkpoint stores | 26 |
| zero-denominator test | 8 |
| length, cap test, meter additions, loop | 14 |
| batch inversion share (< 2^18 per >= 2^15 elements) | 8 |
| **total (budget 240 = 0.12 units)** | **230** |

The rare fallback of section 1.1 is not in this price (separate action).

**Prices** (units): ordinary or traversal step 1.12; doubling 2 (7M + 230
instr. < 1.5); section-1.1 fallback 0.25 (<= 500 instr.); random 64-bit word (RNDR) 0.005; sampling attempt 0.025;
trail start [a]P0 + [b]T (two Jacobian scalar multiplications + inversion,
about 4.8e6 instr.) 2^12; report 0.5; hash probe 0.015; hash initialisation
2^57.1 (once); table 2^24 (once, plus its sampling); phase-2 or re-walk
operation (unbatched, incl. inversion and coefficient updates) 133; small-batch
inversion 132; solving 2^9; verification 2^13.

**Budget.** Phase-1 steps <= U1 + 2^60 < 2^126.40233. Hard extras: phase 2
(2^20 (L_max + 65) * 133 = 2^95.4), re-walks (8 * 2 * ((L_max + 65) * 133 +
2^12) = 2^79.4), reports (D_cap * 0.5 = 2^62.3), hash init 2^57.1, table and
verification < 2^25: together < 2^95.5. Random extras (expected): starts
< 2^62.4 * 2^12 = 2^74.4, words/attempts/probes < 2^63, small batches < 2^8:
fallbacks < 2^70.3: < 2^74.5. Headroom 2^126.41 - 2^126.40233 - 2^95.5
> 2^118.8, so the meter halts early with probability < 2^74.5/2^118.8 < 2^-44
(Markov).

**Preprocessing:** jump table (2048 entries x 2^12.1 = 2^23.1, plus sampling),
report-table initialisation 2^57.1, first 2^20 starts (2^32, plus sampling):
expected < 2^57.11, cap 2^57.2 (the expected random sampling part is < 2^21,
far below the 2^53 gap, so the cap halts a run with probability < 2^-32,
Markov); all of it is inside W.

**W <= 2^126.41 on every run.** Memory: D_cap records x 104 bytes at load 1/2
= 2^70.95, plus table (2 MB) and walk states (2^20 x 512 bytes): < 2^72.

## 6. Experiments and limitations

`experiment.py` (stdlib). **Part A** (production BLS12-381; all true): curve,
order, beta/lambda pairing, orbit and unit action, |y| invariance, unique
canonical member, twisted-step class identity, doubling class map, table
entries, coefficient tracking on production steps, and recovery of a planted
x from two production trails forced to merge. **Part B** (toy y^2 = x^3 + 4
over F_q, q = 58493047 or 861101287): the per-trail algorithm of section 2
(exact canonical rule, which section 1.1 always reproduces, escapes, caps, restarts, re-walk solving). Default mode
reports T_c; `"pipeline": true` also applies the scaled production rules
(phase-1 budget 1.06 S0/0.998 class visits, phase 2, report capacity, 8-match
cap) and reports success. Trails run sequentially (same visit distribution
under H1); [a]G + [b]T is evaluated as [a + bx]G for speed (same point).

Offline results (every recovered x verified):

| mode / classes | toy, R | trials | mean T_c vs model | P[T_c <= 1.06 S0] or success |
| --- | --- | ---: | --- | ---: |
| T_c, 6 | q26, 256 | 8000 | +1.2% | 0.431 |
| T_c, 6 | q26, 2048 | 7000 | +1.3% | 0.422 |
| T_c, 6 | q30, 256 | 2400 | +1.2% | 0.421 |
| T_c, 6 | q30, 2048 | 2400 | +0.5% | 0.443 |
| T_c, 6 (repeat) | q26, 256 | 4000 | +2.1% | 0.413 |
| T_c, 1 (control, no automorphism) | q26, 256 | 4000 | +0.0% | 0.428 |
| T_c, 2 (control, negation only) | q26, 256 | 4000 | +1.4% | 0.427 |
| pipeline success, 6 | q26 256 / q26 2048 / q30 2048 | 6000 / 12000 / 3000 | - | 0.449 / 0.441 / 0.425 |

Pooled class-6 P[T_c <= 1.06 S0] = 0.4255 +- 0.0032 (H1 needs 0.4058); pooled
pipeline success 0.441 +- 0.003. Measured 2-cycle rates match 1/(6R) (1/(2R)
for the negation control). A 14-second run has standard error ~0.026.

Limitations: H1 is a heuristic (random mapping with Brent-Pollard correction;
van Oorschot-Wiener; automorphism classes as in Wiener-Zuccherato, GLV,
Duursma-Gaudry-Morain; fruitless-cycle escapes as in Bos-Kleinjung-Lenstra and
Bernstein-Lange-Schwabe), supported only at toy scale; instruction counts are
analytic; generic attack (no TNFS).
