# Six-class (negation + GLV) parallel Pollard rho in BLS12-381 G1

Parallel distinguished-point Pollard rho on the N = (r-1)/6 orbits of G1 under
negation and phi(x, y) = (beta x, y) = [lambda]. The class key is |y|; the
canonical twist index is decided exactly from a truncated fixed-point product,
with a metered exact fallback on near-ties; an R-entry twisted table makes every
step one batched affine G1 addition plus a 128-instruction ledger; fruitless
cycles of length <= 64 are escaped by doubling the minimum-fingerprint class. A
global work meter charges every action before it runs and halts at
W_max = 2^126.293 units, so the stated work bound holds on every run.

Claim: time_log2 = 126.293 (126293/1000), success >= 0.39 (0.3901314), peak
memory < 2^72 bytes, preprocessing <= 2^57.2 units, no advice.

- [Argument, ledger, failure budget](proof.md)
- [Experiment: production identity checks, twist-index check, toy runs](experiment.py)

## What changes relative to the accepted design

The algorithm is that of public submission
1c600df4-b531-401a-b9d8-5124c6bceeb2 (six-class rho, 126.36). Four numbers are
re-derived, each charged:

1. The per-visit ledger falls from 196 to 128 instructions (proof.md section 5),
   so an ordinary step costs 1 + 128/2000 = 1.064 units instead of 1.10.
2. The success target is the floor 0.3905 (was 0.395), i.e.
   S_ops = 1.05 sqrt(2N' ln(1/0.6095)) = 2^126.19938, and the failure budget is
   re-derived at that target: the meter's early-halt bound is 0.000368437 and
   success is 0.3901314.
3. The table is R = 2^16 (was 2048). The 2-cycle rate 1/(6R) becomes 2^-18.585,
   the expected doubling surcharge 2^106.417, the headroom the Markov bound needs
   2^117.46; the table build is charged 2^28 units and 6R entries of memory.
4. The trail cap stays 20 theta^-1 = 20*2^64, with the Gamma(2) tail 21e^-20.

## The exact twist index (proof.md section 1)

Only t is needed, i.e. j = argmin(v0, v1, v2) with v0 = x, v1 = beta x mod p,
v2 = beta^2 x mod p, three distinct integers in [1, p). No canonical coordinates
are formed and no Fp multiplication is done.

* G = floor(2^448 beta/p) is a 7-limb constant; frac(xG/2^448) = v1/p - e with
  0 <= e < 2^-67.
* L' = (the six high halves of x_i G_{5-i} plus the six low halves of
  x_i G_{6-i}) mod 2^64 is 12 limb products (6 umulh, 6 mul) and pins bits
  384..447 of xG to [L', L' + 12]. If L' < 16 or L' > 2^64 - 16, fall back.
* V1 = floor(PH L' / 2^128) with PH = floor(p/2^256) satisfies
  0 <= v1/2^320 - V1 < 3.
* X = the top limb of x errs by < 1. S = X + V1; if |S - PTOP| <= 16
  (PTOP = floor(p/2^320)) the wrap of x + v1 mod p cannot be decided, so fall
  back; else V2 = PTOP - S or 2 PTOP - S is within 5 of v2/2^320.
* If the smallest of X, V1, V2 is below both others by more than 16 it is the
  exact argmin. Otherwise the exact rule runs, a separately metered action
  (0.25 units, <= 500 instructions).

The fallback windows cover < 2^7 of > 2^60 units, so with beta x mod p
equidistributed the fallback rate is < 2^-52 per step; that assumption only moves
the meter's early-halt probability, never correctness. Part A of `experiment.py`
emulates the limb rule and compares it with the exact rule.

## Why doubling escapes close no structural cycles (proof.md section 2)

A trail that returns to class c after m >= 1 doublings needs
u 2^m c + (unit and power-of-2 multiples of table points met on the way) = c for
a unit u. The order of 2 mod r is 2^250.8571 (computed in Part A from the full
factorisation of r - 1), so 1 - u 2^m != 0 mod r for every m < 2^248. The
relation therefore cannot hold by cancellation of the table terms, as it does for
an additive escape multiplier whose negation pairs annihilate: it pins c to one
point. Such returns are coincidences no more frequent than ordinary
self-collisions of a trail, which is part of H1.

## How the experiment is run

```sh
echo '{"curve":"q26","trials":150,"R":64,"dp_bits":4,"seed":1,
       "production_checks":1,"canon_random":2000,"canon_stress":1000,
       "pipeline":true}' | python3 -I experiment.py
```

The declared input runs in about 1.3 s. Output blocks:

* `part_a`: production BLS12-381 checks, all true in the declared run: the curve
  and order, [r-1]P0 = -P0, beta and lambda as cube roots of unity,
  phi(P0) = [lambda]P0, orbit size 6 with the unit action, |y| as the class key,
  uniqueness of the canonical member, the twisted step and the doubling as class
  maps, table entries equal to [u_t^{-1}]M_h, the full factorisation of r-1 and
  the order of 2 mod r, coefficient tracking along production steps, a second
  trail started in the same class, the two merged trails reaching one canonical
  point, recovery of the planted x, and the truncated-product rule against the
  exact rule on random and adversarial near-tie inputs.
* `cost`: the whole cost model of proof.md sections 3-5 from exact integers and
  80-digit decimals, including the exact-integer check of the claimed score
  (2^126293 >= W_up^1000 and 2^126292 < W_up^1000).
* `toy`: the walk of proof.md section 2 on a toy curve y^2 = x^3 + 4 over F_q
  (q = 58493047 or 861101287), with the exact canonical rule, escapes, caps,
  restarts and re-walk solving; `pipeline` also applies the scaled production
  rules (phase-1 budget, phase 2, report capacity, 8-match cap).

Declared run (q26, R = 64, 4 distinguished bits, 150 trials): pipeline success
0.46, fraction with T_c <= S_ops 0.4533, 2-cycle rate 0.00245/step against
1/(6R) = 0.002604, 4 checkpoint escapes, 0 abandoned trails, 0 exceptional
additions.

Offline runs of the same file (R = 2048, 5 distinguished bits, seeds 6101-6103 on
q26 with 1500 trials each and 6201-6202 on q30 with 700 each): q26 gives 0.4247,
0.4093, 0.4127; q30 gives 0.4086, 0.3871. Pooled 0.4114 +- 0.0064 over 5900
trials against the H1 floor 0.3905, and the pooled mean T_c is 0.8% above the
random-mapping model, inside the 5% allowance.

## Caveats

* H1 is a heuristic. Its support is toy scale only, and the toy's R (<= 2048) is
  far below the production 2^16; only N' and the cycle rates depend on R.
* Trails are kept short relative to sqrt(N) in the toy, as in production; a
  150-trial sample carries a standard error near 0.04 on a success rate near 0.4.
* Instruction counts are analytic, not measured. The report store's probe cost is
  an expectation under H1's hash spreading.
* The margin of the success bound above 0.39 is 1.3e-4, exact arithmetic on the
  stated bounds, not a sampling margin.
* This is a generic attack: it uses neither the pairing nor a number field sieve.
