# plain-rho-v6: secp256k1-ecdlp candidate

The negation walk of v5 is replaced by the plain additive walk on the
full group, with constants that mix properly at production scale. The
fruitless-cycle problem disappears as a THEOREM rather than a census,
and the cycle detection machinery of v5 (look-back, ratchet, escapes,
tail bounds) is deleted entirely. Every figure below is computed by
`experiment.py` at production parameters (Part A) and cross-checked at
toy scale (Part B); the declared run records both.

## Summary

| Field | Value |
| --- | --- |
| time_log2 | **16058573/125000 = 128.468584** |
| success_probability (lower bound) | 0.39 (modeled 0.4066) |
| memory_log2_bytes | 105.1 (over-report; see below) |
| preprocessing_log2 | 22.0 (included in W) |
| nonuniform_advice_log2_bytes | 0 |

## Algorithm

- Group: secp256k1: y^2 = x^3 + 7 over the 256-bit prime field, r prime, standard generator. Target T = [x]P for uniform secret
  x in [0, r); recover x.
- Constants: A_j = [c_j]P with
  **c_j = j*Q + 2^j mod r**, j = 0..31, Q = (r - rho)/32,
  rho = r mod 32 = 1. Part A verifies all 32 constants are
  distinct and that the constant differences span 0.96875 of
  the modulus.
- K = 2^16 walkers, each seeded as C = [a_0]P + [b_0]T with fresh
  uniform a_0, b_0, carrying the lazy representation (a, b) with
  C = [a]P + [b]T and b fixed for the walker's lifetime.
- Step: read the state's index j = x_C mod 32 (5 bits of the
  x-coordinate), update C <- C + A_j and a <- a + c_j mod r.
- Budget: t = 1.25 * sqrt(r-1) walker steps total (log2 t =
  128.32193), shared round-robin across the walkers.
- Distinguished points: a state is distinguished when the low 32 bits
  of its x-coordinate are zero (density 2^-32). On reaching one, the
  walker inserts (x_C, y-parity) -> (a, b) into a shared table of
  capacity 2^98 slots. A key hit with a different representation
  is a collision: x = (a2 - a1)/(b1 - b2) mod r is extracted and the
  trial succeeds. A key hit with the identical representation is a
  cycle revisit (see the theorem: impossible at production) and the
  walker reseeds; the branch exists for the toy and as a safety net.
- The table key is the pair (x-coordinate, y-parity). On the full
  group +X and -X are DISTINCT points that share an x-coordinate;
  the parity bit is what makes a key hit a true point collision, so
  the extraction is exact algebra with no false positives (the toy
  measured 0 wrong extractions in
  200 trials plus 100 census trials).

## The no-cycle theorem

**Claim.** No fruitless cycle of length L <= L* =
r/(32*2^31) = 2^220.00 exists. The per-walker budget is
t/2^16 = 2^112.32 steps, which is below L* by a factor of
2^107.7.

**Proof.** The walk is the deterministic map f(C) = C + A_{j(C)}.
A cycle of length L is a word w_1..w_L with
sum_u c_{w_u} = 0 mod r (the points return) where the word is the
sequence of indices visited. Write

- A = sum_u w_u  (each index in 0..31, so 0 <= A <= 31L),
- B = sum_u 2^(w_u)  (so L <= B <= 2^31 L).

Since c_j = jQ + 2^j, the cycle condition is Q*A + B = 0 mod r.
Multiply by 32 and use 32Q = r - rho:

  32B = rho*A  (mod r).

For L <= L* we have 32B <= 32*2^31 L < r and rho*A <= 31*31L < r,
so the congruence is the integer equation 32B = rho*A. But
2^w >= 2w for every w >= 0, hence B >= 2A term by term, and

  32B >= 64A > 31A >= rho*A

whenever A >= 1 - impossible. If A = 0 then every index is 0 and
B = L >= 1, so 32B = 0 is impossible. Part A verifies
32Q = r - rho, the bound arithmetic, and
per_walker_budget < L* numerically.

Consequences: every walker's trail is self-avoiding, so all
t(t-1)/2 state pairs across walkers are cross-walker pairs and
potentially reportable; no look-back, ratchet, epoch, escape step or
tail bound is needed, and the identical-representation branch never
fires at production.

## Mixing, and why the constants have the j*Q term

The lazy representation makes the walker's scalar k = a + b*x mod r
perform a walk with steps c_j. For two walkers A, B the collision
condition at trail positions (u, v) is

  Q*(A_A(u) - A_B(v)) + (B_A(u) - B_B(v)) = -Delta mod r

where Delta is the (uniform) difference of the seeds. If the
constants were only 2^j (as in v5), the left side would be a signed
sum of powers of two confined to a band of width ~2^31 * L, which at
production scale (L = 2^112.3 per walker, band ~2^143)
does NOT wrap the 2^255-size modulus: the collision probability
would collapse to ~2^-111 per walker pair instead of the birthday
rate. The j*Q term repairs exactly this: the constant differences
c_i - c_j = (i-j)Q + (2^i - 2^j) span 0.9688 of the
modulus, the difference walk wraps within ~5 steps, and
the standard parallel-rho birthday model applies:

  P(collision within the budget) = 1 - exp(-t^2 / (2(r-1))) in the
  ideal random-walk limit.

This is the one genuinely heuristic step of the analysis (the same
heuristic that underlies all practical parallel rho, van
Oorschot-Wiener), and Part B measures the realized fraction of the
ideal directly: 94/200 = 0.47
against the same-parameter ideal 0.5420
(ratio 0.867). The toy ratio
additionally absorbs two strictly toy-scale losses that production
does not have: about 1/8 of the toy's pairs are walker-self pairs
(which production excludes by the theorem, and which the toy
realizes at ~zero rate: 1 cycle
reseed in 200 trials), and toy-scale wrap cycles.

## Distinguished points, reporting, and memory

A collision is detected when the later visitor inserts a
distinguished point whose key is already in the table. If the
collision point itself is not distinguished, the later walker has
merged onto the earlier walker's trail; because the step map is
deterministic, the merged pair CANNOT separate, and it reaches the
next shared distinguished point together. The only way a merge goes
unreported is that the walker's personal budget ends first: the
expected gap to the next distinguished point is 2^32 steps against
a personal budget of t/2^16 = 2^112.32 steps, so the
reporting loss is

  2^32 / (t/2^16) = 2^-80.32

plus the table-overflow probability. The supplementary census in the
declared run exercises this machinery at a production-shaped
distinguished-point density (2^-6 at toy scale): it detected
27/100 trials with
0 wrong extractions, and its deficit
relative to the dbits=0 run is exactly the budget-end race scaled
down: at toy scale the gap (2^6 = 64 steps) is comparable to the
per-walker budget (159 steps), and the
observed undetected fraction matches the scaling law
gap/per-walker-budget; at production the same law gives
2^-80.32. (This answers the v5 review item that the
declared run never exercised the production-density machinery.)

Table sizing: expected inserts = t * 2^-32 = 2^96.32; the
capacity is 2^98 slots (3.2x
the mean), so by a Chernoff bound
P(overflow) <= exp(-2^95.9); on overflow the attack
aborts and fails, and this probability is subtracted from the
success probability (it is far below the printed precision).

Memory. Record layout (32-byte x-coordinate, 8 bytes parity/flags,
32-byte a, 32-byte b, 8 bytes metadata): 112
bytes <= 128 bytes per slot, padded to 128. (This corrects the v5
review item that a 32-byte field cannot hold the 381-bit BLS
x-coordinate: here the x-coordinate gets 32 bytes.) Peak memory is
the 2^98 x 128-byte table plus the per-walker state (2^16 walkers
x a few hundred bytes) plus the constants: strictly above
2^105 bytes. As with the partial-cheon submissions, the
reported memory_log2_bytes = 105.1 is the rounded-up over-report of
that peak (98+7 = 2^105 bytes for the table alone); memory
is a supplementary metric and the over-report direction is
conservative.

## Work ledger

Per walker step (all figures are AArch64 instruction counts, the
"other" category of the classical-group-work-v2 model; one reference
operation = batched incomplete affine addition is charged as 1.000
and every walk step performs exactly one):

| Item | Instructions |
| --- | --- |
| index j = x_C mod 32 (5-bit extract) | 2 |
| representation add a += c_j mod r (full-width constant add + reduce) | 32 |
| distinguished-point predicate (low-32-bit OR-reduce of x) | 4 |
| walker state, loop, branch, counters | 69 |
| **other per step** | **107** |
| reference add (the step itself) | 1.000 |

Per-step price: 1 + 107/1000 =
1.107 target-group-operation units. One-offs, all included in W:

- preprocessing: 32 constant scalar multiples + seeding 2^16 walkers
  = 2^22 units (preprocessing_log2 = 22.0),
- distinguished-point processing: 2^96.32 inserts at ~0.02
  units = <= 2^90.7 units,
- one extraction with a modular inverse (~700 instructions).

W = 1.107 * 1.25 * sqrt(r-1) + 2^22 + 2^90.7 + 700,
log2(W) = 128.468584, reported exactly as
time_log2 = 16058573/125000.

## Success probability

- Ideal collision probability within the budget:
  1 - exp(-(1.25 * (1 - reporting loss))^2 / 2) = 0.5421666382.
- Every collision is cross-walker (theorem), every cross-walker key
  hit is a true collision (parity key), the extraction is exact
  algebra, and the toy verified 0 wrong
  extractions.
- Efficiency of the r-adding walk relative to the ideal random walk:
  charged at a hard floor of **0.75**. Support: Teske's analysis of
  r-adding walks with r = 32 partitions; the toy measurement
  0.867 of the all-pair ideal (0.99 of
  the cross-pair ideal, see Mixing) in a regime strictly harder than
  production; and the floor sits 25% below the measurement.
- Claimed: 0.5421666382 * 0.75 = 0.4066249787 >= 0.39.
  The claim is made at the floor; it does not depend on the toy
  measurement beyond the floor's support.

## Response to the v5 review

1. *The v5 census missed the carry families 2^j + 2^j = 2^(j+1), so
   all rates were invalid.* v6 needs no census: cycles of length up
   to 2^220.00 are excluded by an exact integer argument
   (above), so there are no rates to enumerate.
2. *The tail bound covered only even lengths 18..40.* Moot: no
   cycles of any length or parity below L* exist.
3. *The ratchet is not a proven escape.* The ratchet is deleted, as
   are the look-back and the escape step.
4. *DP coalescence/separation was unjustified.* The plain map is
   deterministic: a merged pair provably never separates. The loss
   is the budget-end race 2^32/(t/2^16) = 2^-80.32,
   measured at toy scale by the census and covered by the scaling
   law above.
5. *The declared run does not test the production-density
   machinery.* The declared run now also executes and records the
   dbits=6 sparse-distinguished-point census (density 2^-6), which
   exercises the sparse table, lagged reporting and the cycle-escape
   reseed, and whose undetected-merge fraction matches the derived
   scaling law.
6. *Memory was reported as exactly 2^103 while stated strictly
   above.* Now: peak is strictly above 2^105 bytes and
   memory_log2_bytes is its rounded-up over-report (105.1).
7. *A 32-byte record field cannot hold the 381-bit BLS
   x-coordinate.* Fixed at the layout level: the x-coordinate gets
   32 bytes on this track (48 on BLS12-381 G1), plus 8 bytes
   parity/flags, 32-byte a, 32-byte b, 8 bytes metadata: 120 <= 128.

## Reproduce

```
echo '{"trials": 200, "dbits": 0, "walkers": 8}' | \
  python3 experiment.py
```

Stdlib only, deterministic (fixed seeds), single JSON object on
stdout. `trials` scales the toy Monte Carlo; `dbits` sets the
distinguished-point density of the main attack (0 reproduces the
declared run); `walkers` sets the walker count. The production
checks (Part A) run in every execution and assert the theorem
premises, the mixing span, the ledger arithmetic and the success
bound >= 0.39.
