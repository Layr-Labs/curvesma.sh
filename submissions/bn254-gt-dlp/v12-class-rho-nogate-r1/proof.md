# bn254-gt-dlp: distinguished-point rho on Frobenius classes

Submitted bound: **time_log2 = 124.763426** (exact rational 62381713/500000).
This is a hard cap on every run: every row of section 4 is enforced by a counter. There is
no run gate, and no row is an average.

Success probability is at least 0.39 over uniform x and fresh algorithm coins,
under the disclosed heuristics H1-H3 (section 6). Peak memory below 2^47.000001
bytes, in one ordinary 48-bit address space. Preprocessing 2^23.478998 units,
included. No advice, nothing amortised, no production-size run.

This is a generic square-root method with explicit scalar AArch64 accounting.
It is not a finite-field sieve attack and is not a security estimate.
Notation: p, r, g are the protected modulus, group order and generator;
C = 20,000; N = (r-1)/12; lambda = p mod r.

## 1. Facts used (all checked by `experiment.py`, mode `algebra`)

- p and r are prime, p = 7 (mod 12), and p has order exactly 12 modulo r. So
  lambda is a primitive 12th root of unity mod r and lambda^6 = -1.
- On GT the Frobenius sigma(y) = y^p acts as y -> y^lambda. It has order 12 and
  sigma^6 is inversion. Every y != 1 has a class {sigma^k(y)} of exactly 12
  elements (a fixed point would need lambda^k = 1 with 0 < k < 12). There are
  N = (r-1)/12 classes.
- Tower: Fp2 = Fp[u]/(u^2+1), xi = 9+u, Fp6 = Fp2[v]/(v^3-xi),
  Fp12 = Fp6[w]/(w^2-v), so y = sum_{k<6} e_k w^k with e_k in Fp2. The pinned
  basis maps to it bijectively: e_k = (c_k + 9 c_{k+6}) + c_{k+6} u.
- sigma(y) = sum conj(e_k) gamma^k w^k with gamma = xi^((p-1)/6) = d + e u in Fp2.
  zeta = xi^((p^2-1)/6) lies in Fp and satisfies zeta^2 - zeta + 1 = 0.

**Class keys.** Write e_1 = a + b u and B = a d + b e. Then the real part of the
w^1 coefficient of sigma^k(y) is

    V_{2j} = a zeta^j,   V_{2j+1} = B zeta^j   (j = 0..5).

Proof: sigma(y) has w^1 coefficient conj(e_1) gamma, whose real part is B;
sigma^2 multiplies that coefficient by gamma^(1+p) = zeta, and induction gives
the rest. Since zeta^3 = -1, zeta^2 = zeta - 1, zeta^4 = -zeta and
zeta^5 = 1 - zeta, all twelve values come from **4 multiplications**
(a d, b e, a zeta, B zeta) and 9 additions or negations. They also satisfy
V_k(sigma(y)) = V_{k+1}(y) (indices mod 12).

**Selector.** Values are compared as their stored Montgomery residues (a
bijection of Fp, so the 12-element multiset is still a class invariant). If the
minimum is attained once, K(y) is that value and k*(y) its index; otherwise y has
a **tie**. Then K is constant on a class and k*(sigma^j y) = k*(y) - j. So
c(y) = sigma^{k*(y)}(y) is the same element for every member of the class: its
canonical representative.

## 2. Algorithm

Parameters: R = 2^14 table indices, window w = 24, A = 8 probes, distinguished
bits b = 87, trail cap L_cap = 2^94 slots. Caps: S_max accepted steps (section 4),
slot cap M_max + L_cap, ST_max = 2^40 trail starts, D_max = 2^40 reports,
REC_max = 16 recoveries.

1. **Input.** Decode z and check z != 0 and z^r = 1. If z = 1, output x = 0.
2. **Table.** For h < R draw uniform scalars u_h, v_h from the algorithm's ordinary fresh independent coins (64 bits per load, one RNDR with its status branch, 16 instructions per 64-bit word everywhere) (254-bit draws of 4 words, rejection,
   at most 128 draws each) and set T_h = g^{u_h} z^{v_h}. Store
   TW[h][k] = sigma^{-k}(T_h) and EX[h][k] = (lambda^{-k} u_h, lambda^{-k} v_h)
   mod r for k < 12.
3. **Hash.** Read two 64-bit words of fresh coins as salts s1, s2. khash(K, s) is the splitmix64 finalizer chained
   over the four words of K starting from s. Word w1 = khash(K, s1) gives the index
   h(K) = w1 mod 2^14. K is **distinguished** when bits 14..63 of w1 are zero and the
   low 37 bits of w2 = khash(K, s2) are zero (87 bits).
4. **Trails** (run one after another). Draw (alpha, beta) uniform mod r from fresh coins, as in step 2, and set
   y = g^alpha z^beta. On a tie, drop the trail. Otherwise put the window = (K(y)).
   Then repeat:
   - If K(y) is distinguished, **report** it (step 5) and end the trail. If the trail
     has used L_cap slots, end it without a report.
   - For a = 0..7: compute y' = y * TW[(h(K)+a) mod R][k*(y)] (one slot).
     - On a tie at y', end the trail.
     - If the hash word w1 = khash(K(y'), s1) equals one of the 24 hash words stored in
       the window, reject this probe. Equal keys give equal hash words, so every true
       window repeat is caught; a 64-bit hash collision adds a spurious rejection with
       probability at most 24 x 2^-64 per slot (H2), which is counted in nu.
     - Otherwise accept: add EX[same entry] to (alpha, beta) mod r, set y = y',
       append K(y') to the window (keep the last 24) and count one accepted step.

     If all 8 probes are rejected, end the trail.
5. **Report and recovery.** Set (alpha_c, beta_c) = lambda^{k*}(alpha, beta) mod r.
   These are the exponents of the canonical representative c(y). Look K up in an
   ordered store keyed by the full 4-word K.
   - If K is absent, insert (K, alpha_c, beta_c).
   - If K is present with (alpha', beta') and beta' != beta_c, set
     x = (alpha_c - alpha')/(beta' - beta_c) mod r and check g^x = z on all 12
     coefficients. If the check holds, output x and stop.
   - Otherwise continue.
6. **Stopping.** New trails start only while the accepted-step count is below S_max.
   The trail in progress then runs to its end. Any cap exceeded means failure. No
   restart or amplification is performed.

**Invariant.** y = g^alpha z^beta throughout, because
sigma^{-k}(T_h) = T_h^{lambda^{-k}}. Also c(y) = g^{alpha_c} z^{beta_c}.

**The walk is on classes.** Let c = c(y) and k = k*(y), so y = sigma^{-k}(c). Then
y' = sigma^{-k}(c * T_h), and the class of y' is the class of c * T_h. That class
depends only on the class of y and the index used. So whenever the first probe is
accepted, the next class is F(C) = class(c(C) * T_{h(K(C))}), a function of the
current class C alone. A later probe is used only when F(C) is among the last 24
classes of the trail (section 7).

**Correctness.** Suppose two reports have equal K and come from the same class.
Then their canonical representatives are the same element, so
alpha_c + x beta_c = alpha' + x beta' (mod r), which gives the formula in step 5.
Every output is checked against z, so a false key match (equal K, different
classes) only costs a recovery and can never produce a wrong answer.

## 3. What happens after two trails meet

Suppose a later trail reaches a class C already visited by an earlier trail. That
earlier trail ended (trails are sequential), normally with a report. From C both
trails apply F, except when a window rejection happens in one of them and not the
other. The windows hold different older keys for at most 24 further steps, so the
trails stay together unless such a rejection occurs in that stretch (loss L2). The
later trail then reaches the earlier trail's distinguished class, the store finds
the key, and recovery runs.

**Degeneracy.** Given the true x, the walk depends on each start only through
alpha_0 + x beta_0, and on the table only through u_h + x v_h. Given those
values, each start's beta_0 is still uniform, because (alpha_0, beta_0) ->
(alpha_0 + x beta_0, beta_0) is a bijection. So beta_c of the later trail is
uniform and independent of beta', and beta' = beta_c has probability exactly 1/r.

## 4. Work ledger (every row is a meter-enforced cap)

The unit is 20,000 scalar AArch64 instructions. No operation is priced at one
unit. Instruction counts come from `machine.md`:
- Fp12 product: 13,843 executed instructions.
- Selector: 1,067 executed instructions.
- Per-slot glue (hashing, window test, table addressing, exponent update, meters, control): 357
  executed instructions, the same on every path through a slot.

So one slot costs 15,267 instructions.

| Row | Every-run bound (instructions) | log2 units, rounded up |
| --- | --- | ---: |
| walk slots: 15,267 x (M_max + L_cap) | 722043730438187757203339685725163556407498 | 124.7635 |
| input check z^r: 508 products + 2^16 | 7,097,780 | 8.4713 |
| derived constants (zeta, gamma, lambda^k, Montgomery constants) | 2^22 | 7.7123 |
| table: R x (1,028 products + 2^13 + 2 x 16,384) | 233,825,304,576 | 23.4790 |
| trail starts: 2^40 x (1,017 products + selector + 2 x 16,384 + 2^10) | 15,517,616,510,011,965,440 | 49.4629 |
| reports: 2^40 x 2^15 (store lookup and insert, canonical exponents) | 2^55 | 40.7123 |
| recoveries: 16 x (508 products + 512 Montgomery mod-r products + 2^10) | 113,834,816 | 12.4747 |
| output 2^16, fixed overhead 2^32 (parsing, allocation, code loading) | 4,295,032,832 | 17.7124 |
| coins: 2 salt words x 16 | 32 | -9.2877 |
| **total R** | **722043730438187757203355239370708832801246** | 124.7635 |

- **Main cap.** S_max = 47288635671054831211067821185688338613 (2^125.1528). It is the least
  integer with S^2 >= 2 (6/5) N ln(5/3). M_max = S_max + floor(S_max/2^13) + 1.
- **Score.** W = R/20,000. The claim **124.763426** is log2 W rounded up at the sixth decimal;
  the experiment prints the unrounded value and checks the rounding.
- **Exponentiations.** Each uses square-and-multiply with at most 508 products (squarings are
  charged as products). Each Frobenius map is charged as one product; it really costs 15
  multiplications and 6 negations.
- **Preprocessing.** Constants, input check and table total 2^23.478998 units, included.
- **Coins.** Each 64-bit word of fresh coins costs 16 instructions. A scalar mod r costs at most 128 draws of
  4 words plus 64 for masking and comparison: 128 x (4 x 16 + 64) = 16,384 instructions. Two scalars per table
  entry and per trail start, two salt words in the coins row; there is no other randomness.
- **Memory.**
  - store: 2^40 nodes of 128 bytes (K, alpha_c, beta_c, two child pointers, parent and colour);
  - TW: 75,497,472 bytes; EX: 12,582,912 bytes; stack, scratch, code and constants: under 2^20.
  - Total: 140,737,577,484,288 bytes = 2^47.000001. All addresses are ordinary 64-bit pointers
    below 2^48, and every access is an ordinary counted load or store.

## 5. Success probability

**Main event E.** Count new states (trail starts and accepted steps) in the order they are made.
E is the event that, by the time S_max steps have been accepted, some state's class equals a class
produced earlier by a different trail. Under H1 with mu = 6/5 (section 6),

    P(E) >= 1 - exp(-S_max^2 / (2 (6/5) N)) >= 1 - exp(-ln(5/3)) = 2/5.

For an ideal random map (no walk-quality factor) the same S_max gives 0.4583.

**Losses.** Each row bounds the probability that E happens and recovery still fails.

| | Cause | Bound |
| --- | --- | ---: |
| L1 | the earlier trail reaches L_cap without a distinguished class ((1-2^-87)^(2^94) <= e^-128, H2), or the later trail is cut: its length is its own part before the meeting plus the earlier trail's remaining path, so it exceeds L_cap = 2^94 only if one part exceeds 2^93 (2 e^-64) | 2^-91 |
| L2 | a window rejection in one trail within 24 steps after meeting: <= 2 x 24 x nu, nu <= 1.0002/(12R) + 24/N (section 7) | 2.5e-4 |
| L3 | one of the two trails is caught by an undetected trap (2^94 x 2^-163.5) or loses all 8 probes (2^94 x nu^8); and **any** trail of the run repeats a class within itself beyond the window, wasting budget: over all trails, sum L_i^2/(2N) <= L_cap (S_max + L_cap)/(2N) = 2^-31.86 (H1) | 2^-31 |
| L4 | a tie in either trail: at most 2^95 states at the H3 rate | 2^-100 |
| L5 | beta' = beta_c (exact, 1/r); a stored key from a different class shadows the true one (2^80 pairs x 144/p, H2) | 2^-160 |
| L6 | sampling failure: at most 2^42 scalars, each failing with probability (1 - r/2^254)^128 = 0.2439^128 < 2^-260 | 2^-218 |
| L7 | a cap exceeded: starts or reports (mean 2^38.15, cap 2^40); rejections (mean 2^-17.58 per slot, allowance 2^-13); recoveries (needs 16 false matches) | 2^-1000 |

Total loss < 3e-4. **Success >= 0.4 - 0.0003 = 0.3997 > 0.39.**

The identity input returns 0 directly. For z != 1 the walk depends on z only through
the uniform u_h + x v_h, so the bound holds for every z, not only for uniform z.

## 6. Heuristics, support and limits

**H1 (class walk behaves like a random map, with factor mu = 6/5).** The claim
uses H1 at exactly one point: the bound P(E) >= 2/5 at S_max.
- *Support:* the random-mapping product bound (no repeat among s uniform values
  has probability at most exp(-s(s-1)/(2N))), and the measurements below.
- *Margin:* the ideal map gives 0.458. The protected minimum 0.39 plus the losses
  needs P(E) >= 0.3903, which holds for every mu <= 1.2388. This is the allowance of
  this attack.
- *Context only, not a premise:* adding walks are studied in Teske, Math. Comp. 70
  (2001) 809-825.

**H2 (hash fields and twists).** For distinct classes, the index, the 87
distinguished bits and the twist increment k*(c T_h) mod 12 behave as independent
uniform values. It is used for trail lengths, cap tails, nu and the trap sums.
- The twist increment is exactly uniform for a uniform GT element, because
  sigma(y) is uniform and shifts k* by one. Its independence from the index is the
  heuristic.

**H3 (ties).** A tie needs e_1 of the state to satisfy one of 8 linear conditions:
a = 0, B = 0, or B = zeta^t a for t = 0..5.
- A uniform element of Fp12 meets them with probability at most 8/p, but GT is a
  subgroup, so a rate of at most 2^-195 per state is assumed, not proved.
- It is used only in L4, with a margin above 2^90. Ties end the trail. They never
  produce output.

**Toy evidence** (`toy_offline.json`, and mode `toy` of `experiment.py`). The same
`ClassRho` code runs the whole algorithm in two families of BN-shaped fields. Both
have q = 7 (mod 12), xi = a0 + u neither a square nor a cube in Fq2, and a prime
r' with q of order 12 mod r'.
- **BN family:** q = p(-41) = 99,286,339 and r' = r(-41) = 99,276,253 from the BN
  polynomials, so q is close to r' as in production. xi_toy = 1+u.
- **Generic family:** r' = 67,108,933 with q = 17,600,215,022,347, so q is much
  larger than r'. xi_toy = 2+u.

**What the toy shares with production.** The selector, the window of 24, the 8
probes, the canonical-exponent recovery with its check, the hash function, and the
trail cap ratio L_cap = 128/theta.

**What differs, stated so the scope is clear:**
- (a) R, the distinguished density (6 bits), and the field sizes.
- (b) Keys are ordered as plain residues, not Montgomery residues. Both orders are
  bijections of Fq, so both give class invariants.
- (c) The distinguished test reads the 6 bits of hash word 1 just above the index
  bits. Production reads 50 bits of word 1 and 37 of word 2, so the second word is
  not exercised in the toy.
- (d) The store is a hash map keyed by the full key, not a red-black tree.
- (e) The toy window compares full keys; production compares 64-bit hash words of the
  keys. They differ only on a 64-bit hash collision (at most 24 x 2^-64 per slot).
- (f) The toy enforces the trail cap only. It has no slot, start, report or recovery
  caps. It runs past the budget to record the first repeat, and counts success only
  for trails started before the budget. The sandbox toy measures
  the conditional success of a run.

Seeds were fixed before the runs (instance i of size R uses seed 1000 i + R), and
every instance is reported. Each run went to 4 x ceil(sqrt(2 (6/5) N' ln(5/3))) new
states, so it covers the budget of this claim, S_budget = ceil(sqrt(2 (6/5) N' ln(5/3))).
The table is at S_budget, which is 3,185 for the BN family and 2,619 for the generic
family.

| Field | R | instances | P[first cross-trail repeat <= S_budget] | end-to-end success within S_budget | wrong outputs | rejections per slot (nu bound) | instances with a capped trail |
| --- | ---: | ---: | ---: | ---: | ---: | --- | ---: |
| BN | 16 | 1000 | 0.438 | 0.459 | 0 | 5.21e-3 (5.81e-3) | 23 |
| BN | 64 | 1000 | 0.443 | 0.463 | 0 | 1.26e-3 (1.34e-3) | 15 |
| BN | 256 | 300 | 0.427 | 0.447 | 0 | 3.30e-4 (3.31e-4) | 6 |
| generic | 16 | 1000 | 0.449 | 0.482 | 0 | 5.22e-3 (5.81e-3) | 23 |
| generic | 64 | 1000 | 0.438 | 0.459 | 0 | 1.37e-3 (1.34e-3) | 19 |
| generic | 256 | 300 | 0.460 | 0.487 | 0 | 2.99e-4 (3.32e-4) | 4 |

How to read the toy runs:
- *Every row exceeds the requirement 0.3903 for P(E).* Rows at 1,000 instances sit at
  least 3 binomial standard errors (0.016) above it. The weakest row, 0.427 (BN, R = 256,
  300 instances), is 1.3 standard errors above. The ideal value is 0.458.
- *Capped trails.* The first column also counts the measured toy trails trapped in
  random cycles with no distinguished class. They end at the cap and use up budget.
  Their rate is theta^-2/(2N') per trail, about 2% of instances here, and about
  2^-76 in production.
- *Rejections* are compared with the full printed bound nu (section 7, at each toy's R and N'). Five rows sit
  below it; generic R = 64 is 1.4 binomial standard errors above (4,333 rejections in 3,161,690 slots), and the
  pooled R = 64 rows are consistent with it.
- *Limits:* the runs support H1 and H2 only at toy scale. Extrapolation to N =
  2^250 is the heuristic itself.

## 7. Deterministic cycles and the window

Write x_n = log_g(y_n) and m_h = log_g(T_h) = u_h + x v_h, which are uniform and
independent. A step gives x_{n+1} = x_n + lambda^{-k_n} m_{h_n}. A return to the
class of x_t after l accepted steps means

    lambda^s x_t = x_t + sum_h P_h(lambda) m_h,   P_h = sum over steps using h of zeta12^{-k_i}.

There are two cases.
- **A random coincidence.** If some P_h(lambda) != 0 mod r, the return is one linear
  condition on uniform independent m_h: probability 1/r per pattern. H1 covers it.
- **A deterministic trap.** If all P_h(lambda) = 0, the return happens whatever the
  m_h are. If moreover s != 0, then x_t = 0, which is the identity and is excluded.

**From mod r to Z[zeta12].** Map Z[zeta12] -> Z/r by zeta12 -> lambda. The kernel is
a prime of norm r. Each P_h is a sum of at most l roots of unity, so every complex
embedding has absolute value at most l, and |Norm(P_h)| <= l^4. A nonzero P_h in the
kernel therefore needs l^4 >= r, that is l >= 2^63.399. So for l < 2^63 a trap needs
P_h = 0 in Z[zeta12] for every h used. A single root of unity is not zero, so **every
index in a trap is used at least twice**.

**Trap probability.** Let a trap of length l use d distinct indices (d <= l/2). Under
H2 the indices fall in a fixed set of d values with probability (d/R)^l. For each
index, the twist at its last use is a fresh uniform shift, so the sum over that index
vanishes with probability at most 1/12. Hence

    E(l) = sum_{d=1}^{min(l/2,R)} C(R,d) (d/R)^l 12^-d

bounds the probability that a trap of length l starts at a given step.
- E(2) = 1/(12R). This is the two-cycle: same index twice with twists 6 apart.
- sum_{l=3..24} E(l) = 1.02e-4 x E(2).
- sum_{l=25..512} E(l) = 2^-163.54.
- For l > 512: C(R,d) <= (eR/d)^d gives each term at most (e/12)^d (d/R)^(l-d)
  <= 2^-(l/2) when d <= R/2. Terms with d > R/2 are at most (e/12)^(R/2) < 2^-17000.
  So E(l) <= (l/2) 2^-(l/2) and the tail is below 2^-246.
- Patterns of length 63 bits or more are within-trail repeats (L3).

**How the window handles traps.**
- *Up to 24 steps:* closed within 24 steps, the trap is caught by the window and
  escaped by the next probe.
- *Per-step rejection rate:* nu <= sum_{l=2..24} E(l) + 24/N + 24 x 2^-64 <= 1.0002/(12R) + 24/N + 24 x 2^-64.
- *Longer traps:* a trapped trail uses at most L_cap slots, all inside the slot cap. A trap
  of length 25 to 2^63 starts at a given step with probability at most 2^-163.5 (the sums
  above), so at most 2^94 x 2^-163.5 = 2^-69.5 per trail; longer patterns are within-trail
  repeats. Both are counted in L3.

## 8. Evidence and reproduction

`experiment.py` (standard library) is the one declared experiment.
- **kernels:** generates the priced programs, executes them, checks them against
  big-integer arithmetic and prints the counts used above.
- **algebra:** checks section 1 and the step's class invariance on the pinned
  constants.
- **ledger:** recomputes every integer in section 4, S_max, the success bound, the score
  and the trap sums.
- **stated_figures** (every mode): reads `submission.json` and this file, and checks that
  the claim score equals the computed score, that every figure above appears as computed,
  and that every log2 in the ledger table is rounded up. Any stated-figure or kernel
  mismatch makes the experiment exit with status 1.
- **toy:** runs the algorithm in the BN-family toy q = p(-15), r' = r(-15), with 8
  to 48 instances (48 initially) and R in {16, 64}. It prints the success count with
  its Wilson 95% interval, and an estimate of P(E) from the measured walk constant
  (mu_hat, which uses every instance's first-repeat time, not only a yes/no). At 48
  instances a single run's success fraction still varies by about 0.07, so one run
  below 0.39 is not evidence against the bound. The decisive statistics are the
  offline table of 4,600 instances.

Offline, the table above was produced with the same `ClassRho` code. A small driver
called experiment.py's `bn_toy` or `find_toy`, `Tower` and `ClassRho` with the seeds
stated in `toy_offline.json`, using a longer budget to record the first repeat. Per-instance results are in
`toy_offline.json`.

**Limitations.**
- The success bound is heuristic (H1-H3).
- The instruction counts are executed counts of generated straight-line programs on an
  interpreter, not timings.
- Nothing runs at production size.
- Finite-field sieve methods are expected to be cheaper and are not used here.
