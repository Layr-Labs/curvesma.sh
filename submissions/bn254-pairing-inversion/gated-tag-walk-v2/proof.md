# bn254-pairing-inversion: gated tag-walk rho, output ([x]P0, Q0), 143-instruction step, capped-run H1, averaged H2

Submitted score: **time_log2 = 119.030737** (exact rational 119030737/1000000).

**The score is an expectation, not an every-run bound.** It is the expected total work over uniform
challenges and the algorithm's fresh coins, the same probability space in which the rules measure
success. One fresh 64-bit coin decides before any other work whether the run opens.
- A closed run costs 64 instructions and counts as a failure.
- An opened run costs 2^120.387428 units in expectation.
- On every run the work is at most 2^132.677182 units, enforced by counters (section 5).

The algorithm finds the exponent x of z = g^x and outputs P = [x]P0 in G1 and Q = Q0 in G2, so
e(P, Q) = z. Success probability is at least 0.39 over uniform x and fresh algorithm coins, closed runs
included as failures, under the disclosed heuristics H1 and H2 (section 7). Peak memory is at most
2^61.007988 bytes in one address space with 64-bit pointers. Preprocessing 2^56.249959 units,
included. No advice, nothing amortised, no production-size run.

x is found by a distinguished-point rho walk with 6 multipliers in GT. A walk step makes no field
multiplication: it is an inner product of 72 stored digits with 72 table words on 64-bit registers
(the *tag*), which gives the next multiplier together with a proof that it is the right one, and a
pointer load. One Fp12 product is made every 1,024 steps. This is a generic square-root method with
explicit scalar AArch64 accounting. It is not a finite-field sieve attack and is not a security estimate.

Notation: p, r, P0, Q0 and g = e(P0, Q0) are the protected modulus, group order and generators;
C = 20,000; m = 6 multipliers; L = 1,024 steps per window; D = 90 distinguished bits.

## 1. The walk is a function of the group element

An element Y of Fp12 is stored as 12 coordinates y_0..y_11 in Fp (the tower order of `machine.md`), each as
its 4-word Montgomery residue Y~_k = y_k 2^256 mod p. Let t(Y) in [0, p) be the stored residue of
coordinate 1. For multipliers M_0..M_5 define

    s(Y) = floor(6 t(Y) / p)  in {0..5},      f(Y) = Y * M_{s(Y)},
    Y is distinguished  iff  6 t(Y) - s(Y) p >= p - floor(p / 2^90).

s, f and the distinguished predicate depend on the element Y only. The algorithm below evaluates exactly
these functions; the table is only a way to evaluate s cheaply. So the walk is an ordinary 6-adding walk
with a deterministic iteration function, and two walks that reach the same element coincide from there on.

## 2. Tags (proved; checked by `experiment.py`, modes `algebra` and `kernels`)

**Row.** For a fixed T in Fp12 the map Y -> (coordinate 1 of Y * T) is Fp-linear, so there is a row
v(T) in Fp^12 with (Y T)_1 = sum_k y_k v_k mod p. Each v_k is a fixed combination of the coordinates of T
with coefficients in {0, 1, -1, 9, -9} (coordinate 1 of the products of two basis elements). Multiplying
by 2^256 gives t(Y T) = sum_k Y~_k v_k mod p.

**Digits.** Write each stored residue in 6 digits of 43 bits: Y~_k = sum_{j<6} d_{k,j} 2^(43 j) with
0 <= d_{k,j} < 2^43 (6 x 43 = 258 bits cover p < 2^254). There are 72 digits.

**Table words.** For T put u_{k,j} = v_k(T) 2^(43 j) mod p and G_{k,j} = floor(2^64 u_{k,j} / p), a
64-bit word. There are 72 words.

**Tag.** A = sum over the 72 pairs (k, j) of d_{k,j} G_{k,j}, taken mod 2^64.

**Lemma.** Let tau = 2^64 t(Y T) / p (a real number in [0, 2^64)) and E = 72 x 2^43. Then
(tau - A) mod 2^64 lies in [0, E).

*Proof.* sum_{k,j} d_{k,j} u_{k,j} = sum_k Y~_k v_k = t(Y T) (mod p), so t(Y T)/p = sum d u / p (mod 1).
For one pair write u/p = G/2^64 + eps with 0 <= eps < 2^-64. Then d u/p = d G/2^64 + d eps. Modulo 1 the
first term equals (d G mod 2^64)/2^64 exactly, because d G is an integer; the second lies in
[0, 2^43 x 2^-64). Summing the 72 pairs: tau/2^64 = A/2^64 + delta (mod 1) with 0 <= delta < 72 x 2^43/2^64.
Multiply by 2^64. QED

**Corollary (certified step).** Let

THR = 18442944161523957759

which is 2^64 - 6 E - 1. Call A *certified* when (6 A mod 2^64) < THR. If A is certified then
floor(6 A / 2^64) = s(Y T), and Y T is not distinguished.

*Proof.* By the lemma tau = A + delta' - w 2^64 with 0 <= delta' < E and w in {0, 1}. If w = 1 then
A = 2^64 - a with 1 <= a < E, and 6 A mod 2^64 = 2^64 - 6 a > 2^64 - 6 E > THR: not certified. So w = 0
and tau = A + delta'. Write 6 A = s' 2^64 + c with c = 6 A mod 2^64 < THR. Then
s' 2^64 <= 6 tau < s' 2^64 + c + 6 E < s' 2^64 + 2^64 - 1. Since 6 tau / 2^64 = 6 t(Y T)/p, this gives
s(Y T) = s' = floor(6 A / 2^64). It also gives 6 t - s p < p (2^64 - 1)/2^64 = p - p/2^64, which is below
p - floor(p/2^90), so Y T is not distinguished. QED

So every distinguished element has an uncertified tag, and a certified tag gives the exact selector.
For a uniform 64-bit tag the uncertified share is at most (6 E + 2)/2^64, which is below 2^-12
(6 A mod 2^64 is even, and each even value comes from two tags).

## 3. Algorithm

Caps: X_cap exits from the step loop (section 5), 2^97 loop trips in one walk, 2^50 walk starts,
2^50 reports, 16 recoveries.

0. **Gate.** Before any other work, read one 64-bit word u of fresh coins. If u >= K, stop and
   return no answer; this costs 64 instructions and is a failure. Otherwise the run is **opened**
   and continues with step 1.

K = 7203029819100908579

   K is the least integer with (K/2^64) x 1047295/1048576 >= 0.39; the experiment checks K and K - 1
   in exact integers. K/2^64 = 0.390478 (rounded up).
1. **Input.** Decode z and check z != 0 and z^r = 1. If z = 1, output P = infinity (null) and Q = Q0 and stop.
2. **Multipliers and table.** For s < 6 draw uniform scalars c_s, d_s from the algorithm's ordinary fresh independent coins (64 bits per load, one RNDR with its status branch, 16 instructions per 64-bit word everywhere) (254-bit draws of 4 words, rejection,
   at most 128 draws each; a scalar that fails 128 times fails the run) and set M_s = g^{c_s} z^{d_s}.
   Build one table entry for every multiset of at most L multipliers:

The table has 1,634,386,955,681,025 entries (the binomial coefficient C(1030, 6)).

   The entry e of a multiset with product T_e holds the 72 words G(T_e), the element T_e, the labels
   alpha_e = sum of its c_s and beta_e = sum of its d_s mod r, and six child pointers: the entry of the
   multiset with one more copy of M_s, or null when the multiset already has L elements. The root is
   the empty multiset, T = 1.
3. **Walks** (run one after another). Draw (a, b) uniform mod r as in step 2 and set Y0 = g^a z^b; if
   Y0 = 1 skip this walk. Put e = root. The pair (Y0, e) represents the element Y = Y0 T_e with label
   (a + alpha_e, b + beta_e). Repeat:
   - (a) Compute the tag A of the digits of Y0 and the words of e.
   - (b) If A is certified, s = floor(6 A / 2^64). If e has a child for s, put e = that child.
     Otherwise (a *window end*) **re-base**: Y0 <- Y0 T_e, (a, b) <- (a + alpha_e, b + beta_e), new
     digits, and e = the root's child for s. Either way this is one step Y <- Y M_s.
   - (c) If A is not certified, re-base to e = root, read t = t(Y0) from the stored coordinate and
     evaluate s(Y) and the distinguished predicate by their definitions (section 1). If Y is
     distinguished, **report** it (step 4) and end the walk. Otherwise put e = the root's child for s
     (one step).

   (a) and (b) are the step loop `tt_run`; a window end or an uncertified tag is an **exit** from it.
   Every exit adds 1 to the run's exit counter; the walk's trip counter counts loop trips. A walk
   that reaches 2^97 trips ends without a report.
4. **Report and recovery.** Look the element Y0 (all 12 coordinates, 48 words) up in an ordered store.
   - If it is absent, insert (Y0, a, b).
   - If it is present with (a', b') and b' != b, set x = (a - a')/(b' - b) mod r, check g^x = z on
     all 12 coefficients and go to step 6.
   - If it is present with b' = b, continue with the next walk.
5. **Stopping.** Walks follow one another until step 4 finds x. There is no step budget. The run stops
   without an answer, a failure, as soon as a cap is reached: X_cap exits (tested at every exit),
   2^50 walk starts, 2^50 reports or 16 recoveries (the last cap is a safeguard: it is never reached,
   see section 4). No restart or amplification is performed.
6. **Output.** Compute P = [x]P0 by left-to-right double-and-add from the top bit of x in
   Jacobian coordinates (`outstep.py`: a doubling uses 7 products, a mixed addition of P0 uses 11),
   convert P to affine with one inversion by exponentiation to p - 2, and output P and Q = Q0. Before
   each step the accumulator is [k]P0 with 1 <= k < x < r, so a doubling never meets infinity and the
   added P0 is never equal or opposite to it. As attack-side verification, evaluate the reduced pairing
   e(P, Q0) (optimal ate Miller loop over the profile's 6x + 2, then the final exponentiation by binary
   powering of (p^12 - 1)/r) and compare it with z on all 12 coefficients; they agree because step 4
   checked g^x = z.

## 4. Correctness

- **The walk.** By the corollary, step 3(b) applies exactly f and never passes a distinguished element;
  step 3(c) uses the definitions directly. So every walk is the trail of f from its start to its first
  distinguished element, whatever representation (Y0, e) the element has at the moment.
- **Labels.** Y = g^a z^b holds for the represented label at each start and is preserved by a step
  (multiplying by M_s adds (c_s, d_s)) and by a re-base.
- **Recovery.** The store key is the element itself. A key found in the store is the same element of GT
  reached twice: g^a z^b = g^a' z^b', so x = (a - a')/(b' - b) whenever b' != b (r is prime). There are
  no false matches, so every recovery outputs the true x; the check g^x = z is kept as a safeguard.
- **Output.** By bilinearity e([x]P0, Q0) = e(P0, Q0)^x = g^x = z; P lies in G1 (cofactor 1) and Q0 is the
  profile generator of G2. The identity z = 1 is answered in step 1 with the infinity output the rules
  permit exactly for z = 1.
- **After two walks meet.** Suppose a later walk reaches an element that an earlier walk visited (a
  **meeting**). f is deterministic, so the later walk follows the earlier one to its distinguished
  element, finds it in the store, and step 4 gives x. A meeting fails to end the run only in two cases:

| | Cause | Bound |
| --- | --- | ---: |
| M1 | the earlier walk was cut by the walk cap (2^97 trips) before its distinguished element, or the later walk is cut: its own part before the meeting plus the earlier walk's remaining path exceeds 2^97 trips. Such a walk has more than 2^96 steps (a walk makes at most one trip more than steps per exit), against about 2^90 steps for a walk that reaches its distinguished element; the effect on the stopping time is inside H1, which is stated for the capped run | inside H1 |
| M2 | b' = b: exactly 1/r, see below | 2^-253 |

- **Degeneracy (exact).** Fix x and condition on the values m_s = c_s + x d_s of the multipliers and
  sigma_i = a_i + x b_i of all walk starts. Every element, tag, exit and store match is a function of
  these values. Given them, each b_i is still uniform and independent of the d_s and of the other b_j,
  because (a_i, b_i) -> (a_i + x b_i, b_i) is a bijection. For a match between walk i and an earlier
  walk j, b - b' = b_i - b_j + (an integer combination of the d_s), which is uniform on Z/r because b_i is.
  It vanishes with probability exactly 1/r.

## 5. Work ledger

The unit is 20,000 scalar AArch64 instructions. No operation is priced at one unit. From `machine.md`
(generated programs, executed by the experiment):
- one walk step (one trip of the step loop): **143** executed instructions (20 of the 72 digits stay in
  registers for a whole call, `machine.md`);
- one exit: **15,162** instructions: the re-base program (one Fp12 product of 13,843, two label
  additions, 72 digits, both counters: 14,115 executed), the loop's entry and exit (23) and an
  allowance of 2^10 for the dispatch, the exact evaluation of step 3(c) and the pointer load.

**Counts.**

E_steps = 215308847512201500215711236325313214167

X_cap = 1076544237561007501078556181626566074368

- E_steps (below 2^127.3397) is the least integer that is at least (106/100) sqrt(pi r m / (2 (m - 1))) with
  m = 6. By H1 it bounds the expected number of walk steps of an opened run.
- Expected exits: at most ceil(E_steps/1,024) + floor(E_steps/2^12) + 1, by two terms.
  - *Window ends.* After an exit the entry has size 1, and a window end needs an entry of size 1,024. So
    1,023 steps lie between an exit and the next window end, and the window end completes one more
    step: there are at most steps/1,024 window ends.
  - *Uncertified tags.* By H2 (a) the expected number of uncertified tags in an opened run is at most
    E_steps/2^12. For a tag that behaves like a uniform 64-bit value the uncertified share is at most
    (6 E + 2)/2^64 < 2^-12.24 (section 2); the margin between 2^-12.24 and 2^-12, above 2^112 tags, absorbs
    tags that repeat earlier values (section 7).
- H1 is stated for the run as executed, with all its caps (section 7), so E_steps bounds the expected
  number of steps of the run that is charged; its exits are bounded by the two terms above.
- X_cap (below 2^129.6617) is 2^12 times that bound, the every-run cap on exits. `rebase` counts every exit
  in a 3-word counter and compares it with X_cap.
- Every call of the step loop makes at most 1,025 trips (its trip budget, `machine.md`) and ends in an
  exit. So an opened run makes at most 1,025 X_cap trips.
- An opened run therefore costs on average 161.51 instructions per walk step in the two walk rows.

**Rows.** The second and third columns are the expectation of an opened run (instructions, log2 units
rounded up). The last two are the bound that holds on every run. Rows other than the two walk rows
are charged at their every-run caps in both.

| Row | Opened run, expectation | log2 units | Every run, at most | log2 units |
| --- | ---: | ---: | ---: | ---: |
| gate | 64 | -8.2877 | 64 | -8.2877 |
| steps | 30,789,165,194,244,814,530,846,706,794,519,789,625,881 | 120.2119 | 157,794,471,620,504,674,470,589,372,321,913,922,350,489,600 | 132.5352 |
| exits | 3,985,000,910,620,116,145,349,870,318,804,197,953,996 | 117.2621 | 16,322,563,729,899,995,731,353,068,825,821,994,819,567,616 | 129.2621 |
| input_check | 7,097,780 | 8.4713 | 7,097,780 | 8.4713 |
| constants | 4,194,304 | 7.7123 | 4,194,304 | 7.7123 |
| table | 1,713,778,936,440,271,136,994 | 56.2500 | 1,713,778,936,440,271,136,994 | 56.2500 |
| walk_starts | 15,888,837,971,051,651,530,752 | 59.4628 | 15,888,837,971,051,651,530,752 | 59.4628 |
| reports | 295,147,905,179,352,825,856 | 53.7123 | 295,147,905,179,352,825,856 | 53.7123 |
| recoveries | 113,834,816 | 12.4747 | 113,834,816 | 12.4747 |
| output_fixed | 4,295,032,832 | 17.7124 | 4,295,032,832 | 17.7124 |
| point_output | 4,194,304 | 7.7123 | 4,194,304 | 7.7123 |
| pairing_check | 268,435,456 | 13.7123 | 268,435,456 | 13.7123 |
| **opened run** | **34,774,166,104,864,930,694,094,341,925,999,955,863,035** | 120.3875 | **174,117,035,350,404,670,201,960,338,912,548,593,138,340,374** | 132.6772 |

Row contents: gate: the coin, compare and return. steps: 143 x steps (expectation) or x 1,025 X_cap
(every run). exits: 15,162 x exits. input_check: 508 products + 2^16. constants: Montgomery and row
constants. table: 2^20 per entry (itemised in `machine.md`) plus 6 multipliers at 1,017 products and two
scalars each. walk_starts: 2^50 x (1,017 products + 2 x 16,384 + 2^10). reports: 2^50 x 2^18 (store
lookup and insert on 48-word keys). recoveries: 16 x (508 products + 512 Montgomery mod-r products +
2^10). output_fixed: output 2^16, parsing, allocation and code loading 2^32.

- *point output:* 2^22 for step 6's scalar multiplication: at most 2 x 253 Jacobian steps for x < 2^254,
  each at most 11 products and 14 sums, every field operation charged at the executed `fp_mul` count, plus
  128 per step for the bit test, branch and call; the inversion (at most 2 x 254 products) and 4 products
  of the affine conversion; 4,096 for encoding and control:

Scalar multiplication: 2 x 253 x ((11 + 14) x 155 + 128) + (2 x 254 + 4) x 155 + 4,096 = 2,108,974 <= 2^22.

- *pairing check:* 2^28 for the attack-side check of step 6. Fp12 products, each charged at the executed
  Fp12 product: the Miller loop over the 65 bits of 6x + 2 (per bit a squaring and a line product, per
  set bit one more line product, two final line products for the Frobenius images of Q0) and the final
  exponentiation by binary powering of (p^12 - 1)/r, 2,790 bits, with 128 per exponent bit for the bit
  test, branch and call. G2 steps (a doubling per bit, an addition per set bit, two final additions), each
  in affine coordinates with its line coefficients: at most one Fp2 inversion (two squarings, one sum, an
  Fp inversion by exponentiation, two products, one negation), 12 Fp2 products (4 products and 2 sums
  each, enough for the Frobenius images of Q0 too) and 40 Fp2 sums (2 sums each), every field operation
  charged at `fp_mul`, plus 4,096 per step for embedding the line and control; 4,096 for the comparison
  with z:

One G2 step: 2 x 254 + 6 + 12 x 6 + 40 x 2 = 666 field operations.

Pairing check: (3 x 65 + 2 + 2 x 2,790) x 13,843 + 2,790 x 128 + (2 x 65 + 2) x (666 x 155 + 4,096) + 4,096 = 94,499,259 <= 2^28.

  The scalar multiplication is executed in every mode with counted field operations (machine.md); the
  pairing check is priced from the executed Fp12 product and `fp_mul` counts and is not executed.
- **Expected total work.** The work of a run is 64 + 1[u < K] x (work of the opened run after the gate).
  The coin u is independent of the challenge and of every other coin, so the expectation is
  64 + (K/2^64) x (opened-run expectation - 64). W = that/20,000, and the claim **119.030737** is
  log2 W rounded up at the sixth decimal; the experiment prints the unrounded value and checks the
  rounding.
- **Every run.** A closed run costs 64 instructions. An opened run costs at most the last column's
  total, 2^132.677182 units: every loop of the algorithm is bounded by one of the counters.
- **Exponentiations.** Each uses square-and-multiply with at most 508 products (squarings are
  charged as products).
- **Preprocessing.** Constants, input check and table total 2^56.249959 units, included.
- **Coins.** Each 64-bit word of fresh coins costs 16 instructions. A scalar mod r costs at most 128 draws of
  4 words plus 64 for masking and comparison: 128 x (4 x 16 + 64) = 16,384 instructions. Two scalars per
  multiplier and per walk start, one word in the gate; there is no other randomness.
- **Memory.** The table (1,088 bytes per entry) and the store (2^50 nodes of 480 bytes: element, labels,
  tree pointers) total 2,318,644,963,066,463,296 bytes = 2^61.007988 at the report cap. All addresses are
  ordinary 64-bit pointers, and every access is an ordinary counted load or store.

## 6. Success probability

**An opened run fails only by reaching a cap**, because it otherwise continues until step 4 finds x,
and a found x is the true one (section 4). The probability that it reaches a cap:

| | Cap | Bound |
| --- | --- | ---: |
| F1 | X_cap exits. By H1 and H2 the expected number of exits is at most X_cap/2^12 (section 5) (Markov's inequality) | 2^-12 |
| F2 | 2^50 walk starts or 2^50 reports. By H2 (b) the expected number of reports is at most E_steps/2^89 <= 2^38.34. Each walk ends with a report, at the walk cap (more than 2^96 steps, so at most E_steps/2^96 such walks in expectation) or at a start equal to 1 (probability 1/r per start), and the last walk ends at the output; so the expected number of walk starts is at most E_steps/2^89 + E_steps/2^96 + 2. Markov's inequality for both caps gives 2^-10.65 | 2^-10 |
| F3 | a scalar draw fails 128 times: at most 2^52 draws, each (1 - r/2^254)^128 = 0.2439^128 < 2^-260. The recovery cap is never reached: every recovery outputs x | 2^-20 |

So an opened run succeeds with probability at least 1 - 2^-12 - 2^-10 - 2^-20 = 1047295/1048576.
The gate coin is independent of the run, so the algorithm succeeds with probability at least
(K/2^64) x 1047295/1048576 >= 0.39. A closed run is a failure and is counted as one.

The cap is the only way an opened run fails, and a run that reaches no meeting stops at the cap: this is
event F1 and is inside the 2^-12. No part of the argument assumes that a cap is reached rarely for any
other reason. A meeting that fails (M1, M2) does not fail the run: the run continues to the next meeting,
and this is part of what H1 measures.

The identity input is answered directly in an opened run. For z != 1 the walk depends on z only through
the uniform values c_s + x d_s and a_i + x b_i, so the bounds hold for every z, not only for uniform z.

## 7. Heuristics, support and limits

**H1 (mean stopping time).** An opened run, executed as in section 3 with all its caps, makes at most
(106/100) sqrt(pi r m / (2 (m - 1))) walk steps in expectation, m = 6.
- *Where it is used:* the two walk rows of the expectation column, and F1 (Markov's inequality is applied
  to the exits of this capped run).
- *Finite by construction.* Every opened run makes fewer than 1,025 X_cap < 2^140 steps (section 5), so
  the expectation exists whatever the coins. Exceptional draws are inside the bound: some multiplier equals
  1 with probability at most 6/r < 2^-251, and such a run costs fewer than 2^140 steps, under 2^-111 steps in
  expectation. H1 is not stated for the process without caps: its stopping time has infinite expectation
  (for instance when every multiplier equals 1 and a walk never reaches a distinguished element).
- *Why this form.* A random map on r points produces its first repeated point after sqrt(pi r/2) new
  points in expectation. In an m-adding walk two different elements can have the same image only if
  they use different multipliers, which removes a fraction 1/m of the pairs: the walk behaves like a
  random map on r m/(m - 1) points. After the meeting the later walk needs the rest of the earlier
  walk's path, about 2^90 steps against 2^127.
- *Why 106/100.* The toy runs below measure the stopping time of the complete algorithm, with its toy exit
  and walk caps, directly: over
  9,600 runs the pooled mean is 1.0241 times sqrt(pi r' m/(2 (m - 1))) with standard error 0.0054: the
  walk is slightly slower than the simple model above, and the ledger uses the allowance, not the model.
  The allowance 1.06 is above that mean plus 5 standard errors (1.0511); the experiment checks this.
  The toy figure contains the toy's larger overhead: its distinguished rate is 1/64, so the path from a
  meeting to the distinguished element is about 0.5% of a run, against 2^-37 of a run in production.
- *Limit.* It is measured at r' near 2^26 and used at r near 2^254. If the true mean were larger by a
  factor f, the expected work would grow by f and F1 would become f 2^-12.
- *Context only, not a premise:* adding walks are studied in Teske, Math. Comp. 70 (2001) 809-825; the
  tag method follows Cheon, Hong and Kim, J. Cryptology 25 (2012) 195-242.

**H2 (averaged rates).** In an opened run, executed with its caps: (a) the expected number of steps whose tag
is not certified is at most E_steps/2^12; (b) the expected number of reports (distinguished elements reached)
is at most E_steps/2^89.
- *Where it is used:* the exits row and F1 (a), F2 (b).
- *Why these figures.* For a tag that behaves like a uniform 64-bit value, the uncertified share is at most
  (6 E + 2)/2^64 < 2^-12.24 (section 2), and an element whose coordinate t behaves like a uniform value on the
  tested bands is distinguished with probability floor(p/2^90)/p, between 2^-90.0001 and 2^-90 (the number of
  residues t with 6 t mod p in the top floor(p/2^90) values, over p). The figures of H2 allow 2^0.24 = 1.18
  times the first share and twice the second. The allowances are for steps that repeat values seen before,
  which are not fresh draws: after the meeting that ends the run, the later walk retraces the earlier walk's
  path to its distinguished element (about 2^90 steps, against 2^127 steps in a run); a walk can close a
  cycle of its own before a distinguished element (probability about (2^90)^2/r = 2^-74 per walk, and the
  walk cap then ends it); and an element can share t with a visited element of its Frobenius orbit (below),
  which happens about as often as a collision, a few times per run. None of these approaches the margins:
  the first margin is above 2^112 tags and the second above 2^37 reports.
- *H2 is a statement about averages over the run*, not about every element conditional on the whole history:
  a revisited element repeats its tag and status exactly, and the walk's algebraic structure ties some
  elements together. It is supported by the measured rates of the capped toy runs below, which include all of
  these effects at toy scale.
- *Frobenius structure of t.* t is coordinate 1, the imaginary part of the constant Fp2 coefficient e_0 of
  Y in Fp12 = Fp2[w]/(w^6 - xi). The p^2-power Frobenius fixes Fp2 and maps w to a sixth root of unity
  times w, so it fixes e_0: t(Y^(p^2)) = t(Y). On GT it acts as Y -> Y^(p^2 mod r), of order 6, so t is
  constant on six-element orbits and takes at most 1 + (r - 1)/6 values; t is therefore not a uniform
  residue mod p, and H2 is stated as averaged rates on the tested bands, not as uniformity of t. The selector and the distinguished predicate
  are constant on each orbit, but the walk does not commute with the Frobenius: Y and Y^(p^2) != Y
  have the same s and the different images Y M_s and Y^(p^2) M_s, and (Y M_s)^(p^2) = Y^(p^2) M_s^(p^2)
  differs from Y^(p^2) M_s because M_s^(p^2) != M_s for M_s != 1, so the walk does not map an orbit onto
  an orbit. Two different elements still have the same image only if they use different multipliers, as
  in section 1, so the structure creates no extra collision; it only restricts the values t takes. The experiment checks t(Y^(p^2)) = t(Y) on production elements of GT (mode `algebra`).
- *Support.* The toy fields below have the same structure: q = 7 mod 12 of order 12 mod r', the same
  tower and the same coordinate, so t is fixed by the q^2-power Frobenius there too (checked in mode
  `algebra`), and in the BN family r' is close to q, so t takes about q/6 values as in production. The
  measured uncertified and report rates per step match the band shares within 1.5% (H2 allows factors 1.18
  and 2), and every checked certified step agrees with the exact selector. This supports H2 at toy scale and
  toy band widths only; production equidistribution at 2^-90 is extrapolated, not shown.

**Toy evidence** (`toy_offline.json`, and mode `toy` of `experiment.py`, which runs `toy.py`). The same
`TagRho` code runs the whole algorithm, to its output, in two families of BN-shaped fields Fq12 with a
subgroup of prime order r' (q of order 12 mod r', q = 7 mod 12, the same tower shape):
- **BN family:** q = p(-41) = 99,286,339 and r' = r(-41) = 99,276,253 from the BN polynomials, so q is
  close to r' as in production. One 43-bit digit per coordinate.
- **Generic family:** r' = 67,108,933 with q = 17,600,215,022,347, so q is much larger than r'. Two digits
  per coordinate.

**What the toy shares with production.** The walk of section 1 (selector floor(6 t/q) from stored
coordinate 1, six multipliers g^c z^d), the tag arithmetic on 64-bit words with 43-bit digits, the
certificate rule with the same threshold formula, table entries with child pointers and windows, the
re-base, the exact evaluation of uncertified elements, the store keyed by the whole element, the
recovery and the stopping rule (run until x is found or the exit cap is reached).

**What differs, stated so the scope is clear:**
- (a) The field sizes, and windows of 4 or 8 steps instead of 1,024.
- (b) The distinguished density is 2^-6, not 2^-90. The threshold formula widens the uncertified region
  to contain the distinguished region, so in the toy most uncertified tags are distinguished elements;
  in production the tag error dominates and almost none are.
- (c) Stored residues are plain residues, not Montgomery residues (a fixed bijection of Fq).
- (d) Table entries are built at their first use, with the same contents as the full table.
- (e) The store is a hash map keyed by the whole element (as in step 4), not a red-black tree.
- (f) The toy has no gate and no start or report caps. Its exit cap is 64 (floor(sqrt(r')) + 1) and its
  walk cap is 128 x 2^6 steps; no run reached the exit cap.

Seeds were fixed before the runs (instance i of a row uses seed 4 x 10^7 + 1000 i + 10 x window + digits
per coordinate; `toy.offline_instance` reproduces any instance), and every instance is reported. The runs
use the shipped `toy.py`, whose store is keyed by the whole element. Instances 0..1,199 of each row use the
same seeds as an earlier set of runs whose toy store was keyed by the coordinate t alone: 4,794 of those
4,800 instances give identical step counts, and the 6 that differ (3 in each BN row, where q is close to r',
so unequal distinguished elements share t more often) stop earlier with the whole-element key. Each run went on until it found x. "Model" is
sqrt(pi r' m/(2 (m - 1))) with m = 6. In every instance the first 100 certified steps are compared with the
exact selector and distinguished test computed from the multiplied element (last column: failures of
checks).

| Field | m | window | instances | model | mean steps at the stop / model | standard error | largest | correct outputs | wrong outputs | runs at the toy exit cap | exits per step | uncertified per step | reports per step | selector check failures |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| bn | 6 | 4 | 2400 | 13,680 | 1.0120 | 0.0107 | 3.79 | 2400 | 0 | 0 | 0.2562 | 0.01584 | 0.01580 | 0 of 240,000 |
| bn | 6 | 8 | 2400 | 13,680 | 1.0069 | 0.0106 | 3.07 | 2400 | 0 | 0 | 0.1322 | 0.01582 | 0.01578 | 0 of 240,000 |
| generic | 6 | 4 | 2400 | 11,247 | 1.0460 | 0.0111 | 3.58 | 2400 | 0 | 0 | 0.2562 | 0.01592 | 0.01585 | 0 of 240,000 |
| generic | 6 | 8 | 2400 | 11,247 | 1.0314 | 0.0108 | 3.09 | 2400 | 0 | 0 | 0.1323 | 0.01599 | 0.01591 | 0 of 240,000 |

Pooled over 9,600 runs: mean 1.0241, standard error 0.0054, largest 3.79.

Negative control, the same runs over sqrt(pi r'/2): mean 1.1218.

How to read the toy runs:
- *Stopping time.* The pooled mean plus 5 standard errors is below the allowance 106/100 of H1.
- *Negative control.* Measured against sqrt(pi r'/2), the model without the adding-walk factor, the same
  runs average above the allowance by more than 5 standard errors (last line under the table; the
  experiment checks this). So the factor m/(m - 1) is needed, and with it the allowance holds. This is
  the same set of runs read against another model, not a second experiment. H2 is supported by the
  three rate columns and the selector checks; it has no separate control.
- *Every run ends with the correct x.* No run reached the toy exit cap and no output was wrong.
- *Selector.* No certified step disagreed with the exact selector, and none was a distinguished
  element, in all checked steps: the corollary of section 2 at work on real walks.
- *Exits per step* equal 1/window plus the uncertified rate, as section 5 counts them.
- *Uncertified and reports per step* sit near the toy's 2^-6 = 0.0156, as H2 predicts for the toy
  threshold (a walk's last element is counted without a step, which adds one part in 64).
- *Limits:* the runs support H1 and H2 only at toy scale. Extrapolation to r near 2^254, with the
  selector reading one of twelve coordinates of an element of a subgroup of a 3,048-bit field, is the
  heuristic itself.

## 8. Evidence and reproduction

`experiment.py` (standard library, with `kern.py`, `toy.py` and `outstep.py`) is the one declared experiment.
- **kernels:** generates the priced programs, executes them on the interpreter and checks them against
  big-integer arithmetic: the field kernels, the Fp12 product, the step loop on three real windows with
  the lemma, selector, exit and registers checked at each step (including the registers the loop must
  preserve), the re-base with its counters on both
  sides of both caps, both sides of the certificate threshold, and the trip budget.
- **algebra:** checks p, r and g, the row identity, the lemma and the corollary on 60 production-size
  pairs including extreme digits, the size of the structure constants, the digit cover, the
  uncertified fraction, and the Frobenius invariance t(Y^(p^2)) = t(Y) on GT elements (and its toy analogue).
- **ledger:** recomputes every integer in section 5, E_steps, X_cap, K (and that K - 1 fails), the
  failure budget and the score.
- **stated figures and output step** (every mode): the output step is executed for 1, r - 1 and a seeded
  x; every figure and table row of the claim, this file and `machine.md` is compared with the computed
  value; the pooled toy mean plus 5 standard errors must be within the allowance of H1 and the negative
  control outside it. Any mismatch makes the experiment exit with status 1.
- **toy:** runs the algorithm to its output in the BN-family toy q = p(-15), r' = r(-15), with 8 to 32
  instances (32 initially) and windows of 4 or 6. At 32 instances the standard error of the mean is about
  0.08, so one run does not test the allowance. The decisive statistics are the offline runs above.

Offline, the table above was produced with the same `TagRho` code (`toy.offline_instance`) and the seeds
stated in `toy_offline.json`, which holds the step count of every instance.

**Limitations.**
- The score is an expectation. An opened run costs more than the score: 2^120.387428 units in
  expectation and up to 2^132.677182 units.
- The success bound and the expectation are heuristic (H1, H2).
- The table needs 2^61.007988 bytes; nothing of that size can be built today.
- The instruction counts are executed counts of generated programs on an interpreter, not timings. The
  table construction, the store, the exact evaluation and the pairing check are priced by itemised
  bounds and are not executed.
- Nothing runs at production size.
- Finite-field sieve methods are expected to be cheaper and are not used here.
