# Six-class (negation + GLV) parallel Pollard rho in BLS12-381 G1

Parallel distinguished-point Pollard rho on the N = (r-1)/6 orbits of G1 under
negation and phi(x, y) = (beta x, y) = [lambda]. The class key is |y|, computed
by the exact six-limb minimum; the canonical twist index is decided from a
12-limb-product truncated fixed-point product with a metered exact fallback on
near-ties; an R = 2^16 twisted table makes every step one batched affine G1
addition plus a 142-instruction ledger (139 work instructions and 3 of stated
rounding slack); fruitless cycles of length <= 64 are escaped by doubling the
minimum-fingerprint class. A global work meter charges every action before it
runs and halts at W_max = 2^126.302 units, so the stated work bound holds on
every run.

Claim: time_log2 = 126.302 (126302/1000), success >= 0.39 (0.3901659), peak
memory < 2^72 bytes (2^71.8074), preprocessing <= 2^59 units under a meter
sub-cap, no advice.

- [Argument, ledger, failure budget](proof.md)
- [Experiment: production identity checks, rule checks, toy runs](experiment.py)

## What changes relative to the accepted design

The algorithm is that of public submission
1c600df4-b531-401a-b9d8-5124c6bceeb2 (six-class rho, 126.36). Everything below
is re-derived and re-charged.

1. The per-visit ledger is 139 instructions plus 3 of stated slack: the
   |y| key row is the exact 18-instruction minimum, the V1 row is the full
   five-instruction carry-exact sequence of section 1 below, the history row
   carries the 192-bit fingerprint compares, and the distinguished test on
   limb 1 of K costs its own row. An ordinary step is 1 + 142/2000 = 1.071
   units and a traversal step 1.083.
2. The success floor is 0.3906, so S_ops = 1.05 sqrt(2N' ln(1/0.6094)) =
   2^126.19962, and at that target the meter's early-halt bound is 0.000433998
   and success is 0.3901659, a margin of 1.66e-4 above 0.39.
3. The table is R = 2^16. The 2-cycle rate 1/(6R) becomes 2^-18.585, the
   expected doubling surcharge 2^106.394238 and the headroom the Markov bound
   needs 2^117.579312; the table build is charged 2^28 units and 6R entries of
   memory.
4. The report store is hard-bounded: 2^65 slots of 112 bytes with double
   hashing, at most 96 probes per insertion at 0.15 units each (a hard
   2^67.09 units), overflow below 2^-40 under the stated random-oracle model.
   Nothing about the store is an expectation.
5. Trail starts are hard-capped, not expected: at most 2^63.25 starts, charged
   2^75.25 units as a hard extra, with "start cap exceeded" below 2^-40 by
   Chernoff. The expected bucket now holds only words, sampling attempts and
   twist-index fallbacks (2^73.9).
6. The doubling's 2y denominator joins the same Montgomery inversion batch as
   the additions, so the 1.5-unit price is the batch share (7M + 128 = 2,461.3
   instructions, 2,464.3 with slack, below 3,000 = 1.5 units). Doublings never
   form a batch of their own; the only partial flushes are one per processor at
   the run's end, charged 2^28 units.
7. The trail length cap stays 20 theta^-1 = 20*2^64, with the two-trail Gamma(2)
   tail 2*21e^-20 = 8.66e-8.

## The |y| key rule (proof.md section 1)

The class key is K = |y| = min(y, p - y). On six 64-bit limbs with y in
[1, p-1] the rule is exactly

```
d  = p - y                        # 6 subs, borrow propagated
br = borrow(y - ceil(p/2))        # 6 sbcs; br = 1 iff y < ceil(p/2)
K  = br ? y : d                   # 6 csel
```

p is odd, so y >= ceil(p/2) iff y > p/2 iff p-y < y: the rule is exact for every
y and needs no tie fallback, at 6 + 6 + 6 = 18 instructions. A top-bit test is
not available, because p is not a power of two: the top bit of y differs from the
top bit of p-y on 23.1% of y, so two members of a class would get different keys
and different partition indices. `limb_key` in `experiment.py` is this rule,
parameterised by the limb count; Part A checks it at the six-limb production
width on 100,000 real production trail positions and on the orbit members, part B
on 666,672 key comparisons at real toy arrivals.

## The exact twist index and the V1 row (proof.md section 1)

Only t is needed, i.e. j = argmin(v0, v1, v2) with v0 = x, v1 = beta x mod p,
v2 = beta^2 x mod p, three distinct integers in [1, p). No canonical coordinates
are formed and no Fp multiplication is done.

* G = floor(2^448 beta/p) is a 7-limb constant; frac(xG/2^448) = v1/p - e with
  0 <= e < 2^-67.
* L' = (the six high halves of x_i G_{5-i} plus the six low halves of
  x_i G_{6-i}) mod 2^64 is 12 limb products (6 umulh, 6 mul) and pins bits
  384..447 of xG to [L', L' + 12]. If L' < 16 or L' > 2^64 - 16, fall back.
* With PH = floor(p/2^256) = PH1*2^64 + PH0, the row is exactly

  ```
  umulh t0, PH0, L'
  mul   t1, PH1, L'
  umulh t2, PH1, L'
  adds  t3, t0, t1       # carry c out
  adc   V1, t2, xzr      # V1 = t2 + c
  ```

  which computes ((PH0*L' >> 64) + PH1*L') >> 64 exactly: PH1*L' = t2*2^64 + t1,
  so the sum is (t2 + c)*2^64 + t3 and its top word is t2 + c. Both products are
  full 64x64->128 products and the carry of adding the high half of PH0*L' to
  the low half of PH1*L' is added by `adc`; the final shift is folded in, so no
  `lsr` is needed. Five instructions, not four. The truncation is exact, since
  the sum is below 2^128.
* V1 is within 3 of v1/2^320. X = the top limb of x errs by < 1. S = X + V1; if
  |S - PTOP| <= 16 (PTOP = floor(p/2^320)) the wrap of x + v1 mod p cannot be
  decided, so fall back; else V2 = PTOP - S or 2 PTOP - S is within 5 of
  v2/2^320.
* If the smallest of X, V1, V2 is below both others by more than 16 it is the
  exact argmin. Otherwise the exact rule runs, a separately metered action
  (0.25 units, <= 500 instructions).

The fallback windows cover < 2^7 of > 2^60 units, so with beta x mod p
equidistributed the fallback rate is < 2^-52 per step; that assumption only
moves the meter's early-halt probability, never correctness. Part A runs the
6-limb rule against the exact rule on random and adversarial near-tie production
inputs and on 100,000 real production trail positions.

## The report store (proof.md section 3)

Records (K, a, b) are 112 bytes (14 64-bit limbs). The table has 2^65 slots;
h = H1(K), g = H2(K)|1, and an insertion probes (h + i g) mod 2^65 for i = 0..95
in order, writing into the first free slot. If all 96 are occupied the attack
halts with failure, so an insertion is deterministically bounded at 96 probes,
charged 0.15 units each (300 instructions for a two-word address add, the slot
load and the compare), and the whole run's probe charge is at most
D_cap * 96 * 0.15 = 2^67.09 units. The table cannot fill (D_cap = 2^63.25 <
2^65). At load alpha = 2^-1.75 the probability that one insertion sees 96
occupied slots is at most alpha^96 = 2^-168 under the random-oracle model, so the
union over D_cap insertions is below 2^-104.7; the proof claims only 2^-40. The
2^71.81-byte table needs 128-bit byte addresses; page faults and zeroing over
2^59.8074 pages are charged 1 unit each as a hard extra, and the write pass is
the 2^57.84-unit preprocessing row under the 2^59 sub-cap.

## The distinguished rule (proof.md sections 1 and 4)

The production rule tests limb 1 of K (bits 64..127) for zero, theta = 2^-64, at
3 instructions (compare and branch). It is disjoint from the partition index
h = K mod R, which reads bits 0..15, so under H1 (K's bits are pseudo-random in
the class) the distinguished indicator and h are independent without any further
hash assumption. Part A tests the exact bit range, the branch on constructed
keys with h unchanged, and the rate (0 hits on 100,000 production positions, as
expected for a 2^-64 event).

## Why doubling escapes close no structural cycles

A trail that returns to class c after m >= 1 doublings needs
u 2^m c + (unit and power-of-2 multiples of table points met on the way) = c for
a unit u. The order of 2 mod r is 2^250.8571 (computed in Part A from the full
factorisation of r - 1), so 1 - u 2^m != 0 mod r for every m < 2^248. The
relation therefore cannot hold by cancellation of the table terms, as it would
for an additive escape multiplier whose negation pairs annihilate: it pins c to
one point. Such returns are coincidences no more frequent than ordinary
self-collisions of a trail, which is part of H1.

## How the experiment is run

```sh
echo '{"curve":"q26","trials":12000,"R":64,"dp_bits":4,"seed":1,
       "production_checks":1,"canon_random":2000,"canon_stress":1000,
       "samples":1000000,"prod_positions":100000,"hash_slots":16384,
       "pipeline":true,"seconds":280}' | python3 -I experiment.py
```

The declared input runs in about 46 s. Output blocks:

* `part_a`: production BLS12-381 checks, all true in the declared run: the curve
  and order, [r-1]P0 = -P0, beta and lambda as cube roots of unity,
  phi(P0) = [lambda]P0, orbit size 6 with the unit action, |y| as the class key,
  uniqueness of the canonical member, the twisted step and the doubling as class
  maps, table entries equal to [u_t^-1]M_h, the full factorisation of r-1 and
  the order of 2 mod r, coefficient tracking along production steps, a second
  trail started in the same class, the two merged trails reaching one canonical
  point, recovery of the planted x, and all three production rules on 100,000
  real production trail positions at the six-limb width (0 key mismatches, 0
  twist mismatches, 0 twist fallbacks, 0 distinguished-rule mismatches).
* `rules`: the same rule code at one 64-bit limb against the exact rules along
  the real toy trail positions, over 1,000,008 sampled arrivals (333,336
  twist-index, 666,672 key), with the fallback count. Declared run: 1 fallback,
  0 mismatches.
* `report_store`: two-hash insertion at load 1/2 with the 96-probe cap: 8,192
  insertions into 16,384 slots, at most 9 probes used, 0 overflows.
* `cost`: the whole cost model of proof.md sections 3-5 from exact integers and
  90-digit decimals, including the exact-integer check of the claimed score
  (2^126302 >= W_up^1000 and 2^126301 < W_up^1000) and asserts on the doubling
  cap (2,464.3 <= 3,000), the small-batch count, memory and the success margin.
* `toy`: the walk of proof.md section 2 on a toy curve y^2 = x^3 + 4 over F_q
  (q = 58493047 or 861101287), with the exact canonical rule, escapes, caps,
  restarts and re-walk solving; `pipeline` also applies the scaled production
  rules (phase-1 budget with sequential trails, report capacity, 8-match cap).

Declared run (q26, R = 64, 4 distinguished bits, 12,000 trials): pipeline
success 0.4198 with standard error 0.0045, fraction with T_c <= S_ops 0.4146,
both above the 0.3906 floor; 2-cycle rate 0.002434/step against
1/(6R) = 0.002604, 387 checkpoint escapes, 5 abandoned trails, 3 exceptional
additions. Seeds 2 and 3 of the same file (12,000 trials each) give 0.4125 and
0.4203, so the three runs pool to 0.4175 over 36,000 trials, 10.5 sigma above
the floor.

## Caveats

* H1 is a heuristic. Its support is toy scale only, and the toy's R (64) is far
  below the production 2^16; only N' and the cycle rates depend on R.
* H2 is a random-oracle model for the two store hashes; the toy store measures
  the probe mechanics at load 1/2, not the 2^-168 insertion-overflow event,
  which the proof bounds conservatively at 2^-40.
* Instruction counts are analytic, not measured. Every reported store cost and
  every probe is metered or hard-bounded; nothing in the ledger or the store is
  an expectation.
* The margin of the success bound above 0.39 is 1.66e-4, arithmetic on the
  stated bounds at 90 significant digits, not a sampling margin.
* The toy's cycle detection uses the exact key K, not the 192-bit production
  fingerprint: it tests the exact-key variant, which is the stronger test. The
  fingerprint is a production optimisation whose only extra event is a spurious
  match (< 2^-64, in the failure list), which costs one re-walk.
* The toy runs the production rule code at one 64-bit limb, the limb count a
  26-bit toy field needs; the six-limb instance is exercised in Part A on
  production inputs and on 100,000 real production trail positions.
* This is a generic attack: it uses neither the pairing nor a number field
  sieve.
