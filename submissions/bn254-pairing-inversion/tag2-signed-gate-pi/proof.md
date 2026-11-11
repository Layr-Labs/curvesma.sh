# bn254-pairing-inversion: gated tag-walk rho with a two-level signed tag, output ([x]P0, Q0)

Submitted score: **time_log2 = 118.828325** (exact rational 4753133/40000).

**The score is an expectation, not an every-run bound.** It is the expected total work over uniform
challenges and the algorithm's fresh coins, the same probability space in which the rules measure
success. One fresh 64-bit coin decides before any other work whether the run opens.
- A closed run costs 64 instructions and counts as a failure.
- An opened run costs 2^120.183957 units in expectation.
- On every run the work is at most 2^131.353421 units, enforced by counters (section 5).

The algorithm finds the exponent x of z = g^x and outputs P = [x]P0 in G1 and Q = Q0 in G2, so
e(P, Q) = z. Success probability is at least 0.39 over uniform x and fresh
algorithm coins, closed runs included as failures, under the disclosed heuristics H1 and H2 (section 7).
Peak memory is at most 2^61.353586 bytes in one address space with 64-bit pointers.
Preprocessing 2^56.249959 units, included. No advice, nothing amortised, no production-size run.

**Credit and relation to earlier packages.** This package is a modification of kaankolcu's package
70974bf7 ("bn254-gt-dlp: gated tag-walk rho, 153 instructions per step"): the gated distinguished-point
walk, the table of multisets with child pointers, the run gate, the counters and caps, the ledger layout,
the success analysis, the interpreter and the generated field kernels (`kern.py`), the Fp12 product and the
toy harness are that package's, and the text below follows its structure. The changes are listed in
section 0. Several of them were published earlier by others and are credited here:
- kaankolcu's package 4cea9e8c ("gated tag-walk rho, 121 instructions per step", 118.863948, the accepted
  front of this track when this package was prepared; first published as db7fa0ef) has a second tag level
  for this walk and 16 digits kept in registers, at 121 and 178 instructions. kaankolcu's later submission
  a86327f8 (pending when this package was prepared) keeps 18 digits in registers, at 120 and 177
  instructions, as this package does. This package was written in parallel and differs in its signed
  digits and signed second-level words, its read coordinate, its Frobenius and revisit accounting, its
  multiplier check and its exit cap; the two-level tag and the register-resident digit counts are not
  claimed as new.
- e-man07's submission e65a41db ("20 register-resident digits") published keeping digits in registers
  across the step loop and the exit cap of 2^10 times the expected exits with the failure budget
  1 - 2^-10 - 2^-10 - 2^-20.
- The reviews of kaankolcu's 271f78c0 and 13640b78 and of kshitij-hash's 801ac4a8 raised the points
  answered in sections 1, 3, 5 and 7: the read coordinate is no longer fixed by a Frobenius power; H1 is
  stated for the algorithm as specified, caps included; H2 is stated as a comparison with an explicit
  idealized process on fresh steps; steps at elements related to earlier ones by a Frobenius power
  (inverses included) are counted and charged in full, and a check on the multipliers stops such a step from
  forcing the next one; the toy store is keyed by the whole element.

The walk, its ledger and its evidence are those of the bn254-gt-dlp package submitted alongside this one
(the same files except `outstep.py` and the output rows). The output step (`outstep.py`, the point
output and the priced pairing check, section 3 step 6 and section 5) is kaankolcu's, from 271f78c0,
unchanged.

x is found by a distinguished-point rho walk with 6 multipliers in GT. A walk step makes no field
multiplication: it is an inner product of 60 signed digits with 60 table words on 64-bit registers (the
*level-1 tag*), which gives the next multiplier together with a proof that it is the right one, and a
pointer load. When the level-1 tag cannot certify the selector (about 1 step in 170 measured, at most 1 in 47.8) a
second inner product with 60 more words (the *level-2 tag*) certifies it with a 122-unit error interval. One
Fp12 product is made every 1,024 steps. This is a generic square-root method with explicit scalar AArch64
accounting. It is not a finite-field sieve attack and is not a security estimate.

Notation: p, r, P0, Q0 and g = e(P0, Q0) are the protected modulus, group order and generators; C = 20,000; m = 6 multipliers;
L = 1,024 steps per window; D = 90 distinguished bits.

## 0. Changes from 70974bf7

| Part | 70974bf7 | this package |
| --- | --- | --- |
| digits | 72 unsigned 43-bit | 60 signed 51-bit (5 per coordinate, balanced recoding) |
| tag | one level, error 72 x 2^43 | level 1 with a per-window error bound E1 <= Smax + 1 < 2^55.84; level 2 with error 122 |
| walk step | 153 instructions | 120 instructions, plus 177 more on the level-2 path |
| digits in registers | none | 18 of the 60 for the whole call of the step loop |
| read coordinate | stored coordinate 1 (fixed by the p^2 Frobenius) | stored coordinate 6 (the p^2 Frobenius multiplies it by a primitive sixth root of unity) |
| multipliers | drawn once | redrawn when two are related by a power of the p^2 Frobenius (at most 16 draws) |
| H1 | stated for a run without caps | stated for the algorithm as specified, caps included |
| H2 | uniform given the whole history | a comparison with an idealized process at fresh steps; revisits and Frobenius encounters bounded and charged at the full price of a step, a level-2 step and an exit |
| exit cap and budget | 2^8, 1 - 2^-8 - 2^-10 - 2^-20 | 2^10, 1 - 2^-10 - 2^-10 - 2^-20 |
| table entry | 1,088 bytes | 1,472 bytes (60 G words, 60 V words) |

## 1. The walk is a function of the group element

An element Y of Fp12 is stored as 12 coordinates y_0..y_11 in Fp (the tower order of `machine.md`;
coordinate 6 i + 2 j + k is part k, 0 real and 1 imaginary, of the Fp2 coefficient of v^j w^i), each as its 4-word
Montgomery residue Y~_k = y_k 2^256 mod p. Let t(Y) in [0, p) be the stored residue of **coordinate 6**, the
real part of the constant Fp2 coefficient of w. For multipliers M_0..M_5 define

    s(Y) = floor(6 t(Y) / p)  in {0..5},      f(Y) = Y * M_{s(Y)},
    Y is distinguished  iff  6 t(Y) - s(Y) p >= p - floor(p / 2^90).

s, f and the distinguished predicate depend on the element Y only. The algorithm below evaluates exactly
these functions; the table is only a way to evaluate s cheaply. So the walk is an ordinary 6-adding walk
with a deterministic iteration function, and two walks that reach the same element coincide from there on.

**The read coordinate and the Frobenius powers** (checked by mode `algebra`). For every Y in Fp12 and
i = 0..5, coordinate 6 of Y^(p^(2 i)) is zeta^i t(Y) mod p, with zeta = 21888242871839275220042445260109153167277707414472061641714758635765020556617, a primitive sixth root of
unity in Fp; i = 3 is the inversion, t(Y^-1) = p - t(Y). No other Frobenius power maps t(Y) to +-t(Y). For
coordinate 1, which 70974bf7 read, the same factor is 1: t(Y^(p^2)) = t(Y), so t takes at most 1 + (r - 1)/6
values on GT. For coordinate 6 the values are not collapsed (the toy census of section 7 finds as many
distinct values as a random map), but they are related: the six elements of an orbit
O(Z) = {Z^(p^(2 i))} have the coordinates zeta^i t(Z). The walk is not Frobenius-equivariant, because the
multipliers are random elements: an element related to an earlier one gets a selector and a distinguished
status that are fixed by the earlier element, and then the walk moves on to an element that is in general
unrelated (section 4). Such steps are counted and charged in full (section 5); they are the *encounters*.

## 2. Tags (proved; checked by `experiment.py`, modes `algebra` and `kernels`)

**Row.** For a fixed T in Fp12 the map Y -> (coordinate 6 of Y * T) is Fp-linear, so there is a row v(T) in
Fp^12 with (Y T)_6 = sum_k y_k v_k mod p. Each v_k is a fixed combination of the coordinates of T with
coefficients in {0, 1, -1, 9, -9}. Multiplying by 2^256 gives t(Y T) = sum_k Y~_k v_k mod p.

**Signed digits.** Write each stored residue in 5 digits of 51 bits by balanced recoding: d_0..d_3 in
[-2^50, 2^50) and 0 <= d_4 <= 851317936231195, with Y~_k = sum_{j<5} d_{k,j} 2^(51 j) exactly (carry c' = (c + d + 2^50)
>> 51 from each digit to the next; 5 x 51 = 255 bits cover p < 2^254). There are 60 digits. For the digits of
one window let N = sum of |d| over the negative digits, P = sum of the positive digits, and E1 = N + P + 1.
Then E1 <= Smax + 1 with

Smax = 64259010763220292

which is 12 (4 x 2^50 + 851317936231195) (below 2^55.835).

**Table words.** For T put u_{k,j} = v_k(T) 2^(51 j) mod p, G_{k,j} = floor(2^64 u / p) (64 bits) and
V_{k,j} = floor(2^127 u / p) mod 2^63 (63 bits, so it is a non-negative signed 64-bit integer).

**Level-1 tag.** A1 = (sum of d G over the 60 pairs) - N, taken mod 2^64.

**Level-2 tag.** A2 = A1 + 2 (sum of SMULH(d, V) over the 60 pairs) + N - 1, taken mod 2^64, where
SMULH(d, V) = floor(d V / 2^64) for the signed 64-bit d and V.

**Lemma (two levels).** Let tau = 2^64 t(Y T) / p (a real number in [0, 2^64)). Then (tau - A1) mod 2^64 lies
in [0, E1), and (tau - A2) mod 2^64 lies in [0, 122).

*Proof.* sum d u = sum_k Y~_k v_k = t(Y T) (mod p), so tau = sum d (2^64 u / p) (mod 2^64). For one pair
write 2^64 u/p = G + eps with 0 <= eps < 1. Then d (2^64 u/p) = d G + d eps, and d G is an integer.
*Level 1:* d eps lies in (d, 0] for d < 0 and in [0, d) for d >= 0, so the sum of the d eps lies in (-N, P).
Hence tau - sum d G mod 2^64 lies in (-N, P), and tau - A1 in (0, N + P) within [0, E1).
*Level 2:* 2^63 eps = 2^127 u/p - 2^63 G, and 2^63 G is a multiple of 2^63 below 2^127 u/p, so
V = floor(2^63 eps) and eps = V/2^63 + eta with 0 <= eta < 2^-63. Then d eps = 2 (d V / 2^64) + d eta, and
d V / 2^64 lies in [SMULH(d, V), SMULH(d, V) + 1). |d eta| < 2^50 x 2^-63 = 2^-13. Summing the 60 pairs, the sum
of d eps lies in [S - 60 x 2^-13, S + 120 + 60 x 2^-13) with S = 2 sum SMULH(d, V), a subset of (S - 1, S + 121).
So tau - (sum d G + S - 1) mod 2^64 lies in (0, 122), and sum d G + S - 1 = A2 (mod 2^64). QED

**Corollary (certified step).** For a level with error width E (E1 or 122) put THR = 2^64 - 6 E - 1. Call A
*certified* when (6 A mod 2^64) < THR. If A is certified then floor(6 A / 2^64) = s(Y T), and Y T is not
distinguished. For level 2,

THR2 = 18446744073709550883

and for level 1 the threshold THR1 = 2^64 - 6 E1 - 1 is computed once per window by the re-base.

*Proof* (as in 70974bf7). By the lemma tau = A + delta' - w 2^64 with 0 <= delta' < E and w in {0, 1}. If
w = 1 then A = 2^64 - a with 1 <= a < E, and 6 A mod 2^64 = 2^64 - 6 a > 2^64 - 6 E > THR: not certified. So
w = 0 and tau = A + delta'. Write 6 A = s' 2^64 + c with c = 6 A mod 2^64 < THR. Then s' 2^64 <= 6 tau <
s' 2^64 + c + 6 E < s' 2^64 + 2^64 - 1. Since 6 tau / 2^64 = 6 t(Y T)/p, this gives s(Y T) = s'. It also gives
6 t - s p < p (2^64 - 1)/2^64 = p - p/2^64, which is below p - floor(p/2^90), so Y T is not distinguished. QED

So every distinguished element has a tag that neither level certifies, and a certified tag at either
level gives the exact selector. The step uses level 1 when it certifies, else level 2, else the exact
evaluation of section 3 (c).

**Uncertified shares for a uniform tag.** 6 A mod 2^64 is even and each even value comes from two tags, so
a uniform tag fails a threshold 2^64 - 6 E - 1 with probability at most (6 E + 2)/2^64.
- Level 1: E1 <= Smax + 1 deterministically, so at most q1 = (6 (Smax + 1) + 2)/2^64 = 0.020901 (rounded
  down; the ledger uses the exact rational). The production-size walk of mode `algebra` measures 0.006.
- Level 2: at most (6 x 122 + 2)/2^64 < 2^-64 x 734. With the distinguished elements (at most 2^-90 of the
  elements), an element needs the exact evaluation with probability below 2^-54.

## 3. Algorithm

Caps: X_cap exits from the step loop (section 5), 2^97 loop trips in one walk, 2^50 walk starts, 2^50
reports, 16 recoveries.

0. **Gate.** Before any other work, read one 64-bit word u of fresh coins. If u >= K, stop and
   return no answer; this costs 64 instructions and is a failure. Otherwise the run is **opened**
   and continues with step 1.

K = 7208315804938894124

   K is the least integer with (K/2^64) x 1046527/1048576 >= 0.39; the experiment checks K and K - 1 in exact
   integers. K/2^64 = 0.390764 (rounded up).
1. **Input.** Decode z and check z != 0 and z^r = 1. If z = 1, output P = infinity (null) and Q = Q0 and stop.
2. **Multipliers and table.** For s < 6 draw uniform scalars c_s, d_s from the algorithm's ordinary fresh
   independent coins (64 bits per load, one RNDR with its status branch, 16 instructions per 64-bit word
   everywhere) (254-bit draws of 4 words, rejection, at most 128 draws each; a scalar that fails 128 times
   fails the run) and set M_s = g^{c_s} z^{d_s}. **Multiplier check:** if M_a^(p^(2 i)) = M_b for some
   a, b < 6 and i < 6 with (a, i) != (b, 0) (in particular if some M_s = 1 or M_a = M_b^-1), draw all six
   again; the 16th degenerate draw fails the run. Build one table entry for every multiset of at most L
   multipliers:

The table has 1,634,386,955,681,025 entries (the binomial coefficient C(1030, 6)).

   The entry e of a multiset with product T_e holds the 60 words G(T_e) and the 60 words V(T_e), the element
   T_e, the labels alpha_e = sum of its c_s and beta_e = sum of its d_s mod r, and six child pointers: the
   entry of the multiset with one more copy of M_s, or null when the multiset already has L elements. The
   root is the empty multiset, T = 1.
3. **Walks** (run one after another). Draw (a, b) uniform mod r as in step 2 and set Y0 = g^a z^b; if
   Y0 = 1 skip this walk. Compute the digits of Y0 and the window constants -N, THR1, N - 1. Put e = root.
   The pair (Y0, e) represents the element Y = Y0 T_e with label (a + alpha_e, b + beta_e). Repeat:
   - (a) Compute A1 from the digits of Y0 and the G words of e. If A1 is certified, s = floor(6 A1 / 2^64).
     Otherwise compute A2 with the V words of e; if A2 is certified, s = floor(6 A2 / 2^64).
   - (b) If a level certified s: if e has a child for s, put e = that child. Otherwise (a *window end*)
     **re-base**: Y0 <- Y0 T_e, (a, b) <- (a + alpha_e, b + beta_e), new digits and window constants, and
     e = the root's child for s. Either way this is one step Y <- Y M_s.
   - (c) If neither level certified, re-base to e = root, read t = t(Y0) from the stored coordinate and
     evaluate s(Y) and the distinguished predicate by their definitions (section 1). If Y is
     distinguished, **report** it (step 4) and end the walk. Otherwise put e = the root's child for s (one
     step).

   (a) and (b) are the step loop `tt_run`; a window end or a tag that neither level certifies is an
   **exit** from it. Every exit adds 1 to the run's exit counter; the walk's trip counter counts loop trips.
   A walk that reaches 2^97 trips ends without a report.
4. **Report and recovery.** Look the element Y0 (all 12 coordinates, 48 words) up in an ordered store.
   - If it is absent, insert (Y0, a, b).
   - If it is present with (a', b') and b' != b, set x = (a - a')/(b' - b) mod r, check g^x = z on all 12
     coefficients and go to step 6.
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

- **The walk.** By the corollary, step 3(b) applies exactly f and never passes a distinguished element,
  whichever level certified; step 3(c) uses the definitions directly. So every walk is the trail of f from
  its start to its first distinguished element, whatever representation (Y0, e) the element has.
- **Labels.** Y = g^a z^b holds for the represented label at each start and is preserved by a step
  (multiplying by M_s adds (c_s, d_s)) and by a re-base.
- **Recovery.** The store key is the element itself. A key found in the store is the same element of GT
  reached twice: g^a z^b = g^a' z^b', so x = (a - a')/(b' - b) whenever b' != b (r is prime). There are no
  false matches, so every recovery outputs the true x; the check g^x = z is kept as a safeguard.
- **Output.** By bilinearity e([x]P0, Q0) = e(P0, Q0)^x = g^x = z; P lies in G1 (cofactor 1) and Q0 is the
  profile generator of G2. The identity z = 1 is answered in step 1 with the infinity output the rules
  permit exactly for z = 1.
- **Encounters.** Suppose the walk reaches Y' = Z^(p^(2 i)), i != 0, for an element Z evaluated earlier with
  selector s. The step at Y' is correct (it uses the definitions or a certified tag) and its outcome is fixed
  by t(Z). Its next element Y' M_s' equals the image (Z M_s)^(p^(2 i)) of the earlier next element only if
  M_s' = M_s^(p^(2 i)), which the multiplier check excludes; so the walk does not follow the image of the
  earlier trail, and a further relation to an earlier element is a new coincidence, counted again as an
  encounter or a meeting. An encounter does not report anything it should not: the store compares whole
  elements.
- **After two walks meet.** Suppose a later walk reaches an element that an earlier walk visited (a
  **meeting**). f is deterministic, so the later walk follows the earlier one to its distinguished
  element, finds it in the store, and step 4 gives x. A meeting fails to end the run only in two cases:

| | Cause | Bound |
| --- | --- | ---: |
| M1 | the earlier walk reached 2^97 trips without a distinguished element, or the later walk is cut: its own part before the meeting plus the earlier walk's remaining path exceeds 2^97 trips only if one part exceeds 2^96. A walk makes at most one trip more than steps per exit, so 2^96 trips contain more than 2^95 elements. If these are fresh steps, each is distinguished with probability above 2^-91 in the process of H2: at most 2 (1 - 2^-91)^(2^95) <= 2 e^-16. Otherwise the walk has closed a cycle (at most 2^-70, section 5) or met an orbit (at most 2^4 such steps, each one step long, section 5) | 2^-22 |
| M2 | b' = b: exactly 1/r, see below | 2^-253 |

- **Degeneracy (exact).** Fix x and condition on the values m_s = c_s + x d_s of the multipliers and
  sigma_i = a_i + x b_i of all walk starts. Every element, tag, exit and store match is a function of
  these values. Given them, each b_i is still uniform and independent of the d_s and of the other b_j,
  because (a_i, b_i) -> (a_i + x b_i, b_i) is a bijection. For a match between walk i and an earlier walk j,
  b - b' = b_i - b_j + (an integer combination of the d_s), which is uniform on Z/r because b_i is. It
  vanishes with probability exactly 1/r.

## 5. Work ledger

The unit is 20,000 scalar AArch64 instructions. No operation is priced at one unit. From `machine.md`
(generated programs, executed by the experiment):
- one walk step (one trip of the step loop, level 1): **120** executed instructions;
- a step that needs level 2: **177** more executed instructions (60 SMULH and their sums, the
  second certificate);
- one exit: **15,619** instructions: the re-base program (one Fp12 product of 13,843, two label
  additions, 60 signed digits and the window constants, both counters: 14,574 executed), the loop's
  entry and exit (20, with the save and restore of the six callee-saved registers that hold digits) and an
  allowance of 2^10 for the dispatch, the exact evaluation of step 3(c) and the pointer load.

**Fresh and non-fresh steps.** A step of the run evaluates the iteration function at an element Y. It is a
*revisit* when the run evaluated Y before, an *encounter* when it evaluated some Z with Y = Z^(p^(2 i)),
i != 0, before, and *fresh* otherwise. H2 (section 7) is a statement about fresh steps only. At revisits and
encounters the selector and the distinguished status are fixed by the earlier element, and the tag level
reached may differ (it depends on the window's digits). The ledger therefore charges every non-fresh step
the full price of a step, a level-2 step and an exit, in a row of its own, and uses H2 only for fresh steps.

**Counts.**

E_steps = 217340063054769438896991531007627489772

E_level2 = 4542608948867007142795422443021408713

X_cap = 217340073195986595075066199665591810048

- E_steps (below 2^127.3533) is the least integer that is at least (107/100) sqrt(pi r m / (2 (m - 1))) with
  m = 6. By H1 it bounds the expected number of walk steps of an opened run.
- E_level2 = ceil(q1 E_steps) bounds the expected level-2 steps among fresh steps: in the process of H2 a
  fresh step reaches level 2 with probability at most q1 = (6 (Smax + 1) + 2)/2^64 given everything before
  it, and the number of steps is a stopping time (Wald's identity).
- Expected exits at fresh steps: at most E_exits = ceil(E_steps/1,024) + floor(E_steps/2^54) + 1 =
  212246155326935344969952136417935710, by two terms.
  - *Window ends.* After an exit the entry has size 1, and a window end needs an entry of size 1,024. So
    1,023 steps lie between an exit and the next window end, and the window end completes one more
    step: there are at most steps/1,024 window ends.
  - *Exact evaluations.* By H2 a fresh step needs the exact evaluation with probability below 2^-54
    (section 2); Wald's identity as above.
- *Non-fresh steps.* The expected number of non-fresh steps of an opened run is at most 2^93 (H2, part
  (b)). A revisit happens only after a meeting (with another walk's trail or with the walk's own trail).
  (i) A walk that meets the trail of an earlier walk that ended at a distinguished element follows it to that
  element and the run then stops (unless M2). The trail met is chosen in proportion to its length, and trail
  lengths of fresh steps are geometric with parameter above 2^-91, so the expected length followed is at
  most 2 x 2^91 = 2^92; further meetings happen only after M1 or M2. (ii) A walk meets its own trail (a cycle)
  with probability at most E[l^2]/r <= 2^-70 for a trail of l fresh steps, over at most 2^37.36 expected
  walks, and then runs to the walk cap, 2^97 trips: at most 2^64.4 expected revisits. (iii) Meeting the trail
  of a walk that was cut at the walk cap (expected number of such walks below 2^-8.8, each trail at most
  2^97): at most 2^88.2. (iv) Encounters: an element is related to an earlier one through 5 nontrivial
  Frobenius powers for every way it can equal one, so encounters occur at 5 times the rate of meetings;
  the run stops at its first successful meeting, so the expected number of encounters is about 5 (at most 2^4
  is used), and each is a single step (section 4). The sum is below 2^93. The every-run bound does not use
  this: there every trip is charged at the level-2 price and every exit is counted by the cap.
- H1 is stated for the algorithm as specified, with all its caps, so its steps and exits are those of the
  run that the ledger prices.
- X_cap (below 2^127.3533) is 2^10 times (E_exits + 2^93), the every-run cap on exits; it bounds the
  expected exits, those at non-fresh steps included, by X_cap/2^10. `rebase` counts every exit in a 3-word counter
  and compares it with X_cap.
- Every call of the step loop makes at most 1,025 trips (its trip budget, `machine.md`) and ends in an
  exit. So an opened run makes at most 1,025 X_cap trips.
- An opened run therefore costs on average 138.96 instructions per walk step in the three walk rows
  (steps, level2, exits).

**Rows.** The second and third columns are the expectation of an opened run (instructions, log2 units
rounded up). The last two are the bound that holds on every run. Rows other than the walk rows are charged
at their every-run caps in both.

| Row | Opened run, expectation | log2 units | Every run, at most | log2 units |
| --- | ---: | ---: | ---: | ---: |
| gate | 64 | -8.2877 | 64 | -8.2877 |
| steps | 26,080,807,566,572,332,667,638,983,720,915,298,772,640 | 119.9724 | 26,732,829,003,106,351,194,233,142,558,867,792,635,904,000 | 129.9739 |
| level2 | 804,041,783,949,460,264,274,789,772,414,789,342,201 | 114.9529 | 39,430,922,779,581,868,011,493,885,274,329,994,137,958,400 | 130.5346 |
| exits | 3,315,072,700,051,403,153,085,682,418,711,737,854,490 | 116.9966 | 3,394,634,603,248,114,628,477,458,972,576,878,481,139,712 | 126.9966 |
| non_fresh | 157,624,429,322,128,899,642,355,689,193,472 | 92.6705 | 0 | - |
| input_check | 7,097,780 | 8.4713 | 7,097,780 | 8.4713 |
| constants | 4,194,304 | 7.7123 | 4,194,304 | 7.7123 |
| multiplier_draws | 1,371,442,720 | 16.0654 | 1,371,442,720 | 16.0654 |
| table | 1,713,778,936,440,271,136,994 | 56.2500 | 1,713,778,936,440,271,136,994 | 56.2500 |
| walk_starts | 15,888,837,971,051,651,530,752 | 59.4628 | 15,888,837,971,051,651,530,752 | 59.4628 |
| reports | 295,147,905,179,352,825,856 | 53.7123 | 295,147,905,179,352,825,856 | 53.7123 |
| recoveries | 113,834,816 | 12.4747 | 113,834,816 | 12.4747 |
| output_fixed | 4,295,032,832 | 17.7124 | 4,295,032,832 | 17.7124 |
| point_output | 4,194,304 | 7.7123 | 4,194,304 | 7.7123 |
| pairing_check | 268,435,456 | 13.7123 | 268,435,456 | 13.7123 |
| **opened run** | **30,199,922,208,197,625,425,026,120,367,074,854,888,681** | 120.1840 | **69,558,386,385,936,333,834,222,384,570,587,342,594,727,990** | 131.3535 |

Row contents: gate: the coin, compare and return. steps: 120 x steps (expectation) or x 1,025 X_cap
(every run). level2: 177 x E_level2 (expectation) or x 1,025 X_cap (every run: every trip charged at
level 2). exits: 15,619 x exits. non_fresh: 2^93 x (120 + 177 + 15,619) (expectation;
on every run these steps are already inside the steps, level2 and exits rows). multiplier_draws: 16 draws
of the six multipliers (1,017 products and two scalars each) with the check (2^20 each, `machine.md`). input_check: 508 products +
2^16. constants: Montgomery and row constants. table: 2^20 per entry (itemised in `machine.md`) plus 6
multipliers at 1,017 products and two scalars each. walk_starts: 2^50 x (1,017 products + 2 x 16,384 +
2^10). reports: 2^50 x 2^18 (store lookup and insert on 48-word keys). recoveries: 16 x (508 products + 512
Montgomery mod-r products + 2^10). output_fixed: output 2^16, parsing, allocation and code loading 2^32.

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
  64 + (K/2^64) x (opened-run expectation - 64). W = that/20,000, and the claim **118.828325** is
  log2 W rounded up at the sixth decimal; the experiment prints the unrounded value and checks the
  rounding.
- **Every run.** A closed run costs 64 instructions. An opened run costs at most the last column's
  total, 2^131.353421 units: every loop of the algorithm is bounded by one of the counters.
- **Exponentiations.** Each uses square-and-multiply with at most 508 products (squarings are
  charged as products).
- **Preprocessing.** Constants, input check and table total 2^56.249959 units, included.
- **Coins.** Each 64-bit word of fresh coins costs 16 instructions. A scalar mod r costs at most 128 draws of
  4 words plus 64 for masking and comparison: 128 x (4 x 16 + 64) = 16,384 instructions. Two scalars
  per multiplier and per walk start, one word in the gate; there is no other randomness.
- **Memory.** The table (1,472 bytes per entry) and the store (2^50 nodes of 480 bytes: element, labels,
  tree pointers) total 2,946,249,554,047,976,896 bytes = 2^61.353586 at the report cap. All addresses are ordinary
  64-bit pointers, and every access is an ordinary counted load or store.

## 6. Success probability

**An opened run fails only by reaching a cap**, because it otherwise continues until step 4 finds x,
and a found x is the true one (section 4). The probability that it reaches a cap:

| | Cap | Bound |
| --- | --- | ---: |
| F1 | X_cap exits. By H1 and H2 the expected number of exits, non-fresh steps included, is at most E_exits + 2^93 = X_cap/2^10 (section 5); Markov's inequality | 2^-10 |
| F2 | 2^50 walk starts or 2^50 reports. A fresh step is distinguished with probability at most 2^-90 in the process of H2, and a walk reports at most once, so the expected number of reports is at most E_steps/2^90 + 64 <= 2^37.36 (the 64 covers the walks that end at a non-fresh distinguished element: after a meeting, fewer than 2 in expectation by M1 and M2, or at an encounter, at most 2^4); each walk ends with a report, at the walk cap (below 2^-22, M1) or at a start equal to 1 (1/r), so the same figure bounds the expected walk starts. Markov's inequality for both caps gives 2^-11.64 | 2^-10 |
| F3 | a scalar draw fails 128 times: at most 2^52 draws, each (1 - r/2^254)^128 = 0.2439^128 < 2^-260. Sixteen degenerate draws of the multipliers: each draw is degenerate with probability at most 216/r < 2^-245. The recovery cap is never reached: every recovery outputs x | 2^-20 |

So an opened run succeeds with probability at least 1 - 2^-10 - 2^-10 - 2^-20 = 1046527/1048576.
The gate coin is independent of the run, so the algorithm succeeds with probability at least
(K/2^64) x 1046527/1048576 >= 0.39. A closed run is a failure and is counted as one.

The cap is the only way an opened run fails, and a run that reaches no meeting stops at the cap: this is
event F1 and is inside the 2^-10. No part of the argument assumes that a cap is reached rarely for any
other reason. A meeting that fails (M1, M2) does not fail the run: the run continues to the next meeting,
and this is part of what H1 measures.

The identity input is answered directly in an opened run. For z != 1 the walk depends on z only through
the uniform values c_s + x d_s and a_i + x b_i, so the bounds hold for every z, not only for uniform z.

## 7. Heuristics, support and limits

**H1 (mean stopping time).** An opened run of the algorithm of section 3, with all its caps, makes at most
(107/100) sqrt(pi r m / (2 (m - 1))) walk steps in expectation, m = 6 (all steps, fresh or not). Because the run
is capped, its number of steps is bounded on every run, so the expectation is finite by construction; no
uncapped comparison process is used.
- *Where it is used:* the walk rows of the expectation column, and F1.
- *Why this form.* A random map on r points produces its first repeated point after sqrt(pi r/2) new
  points in expectation. In an m-adding walk two different elements can have the same image only if
  they use different multipliers, which removes a fraction 1/m of the pairs: the walk behaves like a
  random map on r m/(m - 1) points. After the meeting the later walk needs the rest of the earlier
  walk's path, about 2^90 steps against 2^127.
- *Why 107/100.* The toy runs below measure the stopping time of the complete algorithm, with the read
  coordinate 6 and two tag levels, directly: the pooled mean is 1.0250 times sqrt(pi r' m/(2 (m - 1)))
  with standard error 0.0076. The allowance 1.07 is above that mean plus 5 standard errors; the
  experiment checks this. The toy figure contains the toy's larger overhead: its distinguished rate is
  1/64, so the path from a meeting to the distinguished element is about 0.5% of a run, against 2^-37 of a
  run in production. The allowance is the one of 70974bf7.
- *Limit.* It is measured at r' near 2^26 and used at r near 2^254. If the true mean were larger by a
  factor f, the expected work would grow by f and F1 would become f 2^-10.
- *Context only, not a premise:* adding walks are studied in Teske, Math. Comp. 70 (2001) 809-825; the
  tag method follows Cheon, Hong and Kim, J. Cryptology 25 (2012) 195-242.

**H2 (idealized process for fresh steps).** Let P* be the process that runs the algorithm of section 3,
caps included, except that at every fresh step (section 5) the coordinate t of the element is replaced by a
fresh uniform residue mod p, independent of everything before, and the tags are formed from it with the
errors tau - A1 and tau - A2 that the digits of Y0 and the words of e fix; at revisits and encounters P*
uses the values fixed by the earlier element, as the run does. In P* the following hold exactly, given the
whole history: a fresh step reaches level 2 with probability at most q1, needs the exact evaluation with
probability below 2^-54, and is distinguished with probability floor(p/2^90)/p, between 2^-91 and 2^-90.
- (a) The heuristic: the expected numbers of steps at each level, exits and reports of an opened run, and
  the walk-length tails used in M1 and F2, are at most what they are in P*.
- (b) Non-fresh steps: an opened run makes at most 2^93 of them in expectation, by the random-mapping count
  of section 5.
- *Where it is used:* E_level2, the expected exits, the non-fresh row, the expected reports and M1.
- *Why it is a model and not a conditional probability of the real run:* in the real run t is a
  deterministic function of the history, so its conditional law is a point mass; the statement compares
  the run with P*, in which the probabilities above are exact. What is assumed is that the real walk is not
  more expensive than P* on these counts.
- *What it does not need:* uniformity over all of Fp. Only the frequencies of the tested intervals matter,
  and only at fresh steps. The orbit relations of section 1 are excluded from fresh steps by definition and
  are counted separately.

**Toy evidence** (`toy_offline.json`, and mode `toy` of `experiment.py`, which runs `toy.py`). The same
`TagRho` code runs the whole algorithm, to its output, in two families of BN-shaped fields Fq12 with a
subgroup of prime order r' (q of order 12 mod r', q = 7 mod 12, the same tower shape):
- **BN family:** q = p(-41) = 99,286,339 and r' = r(-41) = 99,276,253 from the BN polynomials, so q is
  close to r' as in production. Two signed 14-bit digits per coordinate, 26-bit words.
- **Generic family:** r' = 67,108,933 with q = 17,600,215,022,347, so q is much larger than r'. Two signed
  22-bit digits per coordinate, 35-bit words.
- The toy runs the multiplier check of step 2 without using the random generator, so a run reproduces from
  its seed; none of the 4,800 offline instances drew a degenerate set (`multiplier_check` in `toy_offline.json`).

**What the toy shares with production.** The walk of section 1 (selector floor(6 t/q) from stored
coordinate 6, six multipliers g^c z^d), the two-level tag arithmetic with signed digits, the G and V words,
SMULH-style floor products and the same certificate and threshold formulas, table entries with child
pointers and windows, the re-base with new window constants, the exact evaluation of elements that
neither level certifies, the store keyed by the whole element (all 12 coordinates), the recovery and the
stopping rule (run until x is found or the exit cap is reached).

**What differs, stated so the scope is clear:**
- (a) The field sizes, and windows of 4 or 8 steps instead of 1,024.
- (b) The word size is chosen per field (26 and 35 bits instead of 64) so that the level-1 share is 1 to
  2%, as in production, and both levels occur often enough to be tested.
- (c) The distinguished density is 2^-6, not 2^-90. The threshold formula widens the uncertified region
  to contain the distinguished region, so in the toy most exact evaluations are distinguished elements;
  in production almost none are.
- (d) Stored residues are plain residues, not Montgomery residues (a fixed bijection of Fq).
- (e) Table entries are built at their first use, with the same contents as the full table.
- (f) The store is a hash map keyed by the 12 coordinates, not a red-black tree.
- (g) The toy has no gate and no start or report caps. Its exit cap is 64 sqrt(r') and its walk cap is
  128 x 2^6 steps; no run reached the exit cap.

Seeds were fixed before the runs (instance i of a row uses seed 5 x 10^7 + 1000 i + 10 x window + 2), and
every instance is reported. Each run went on until it found x. "Model" is sqrt(pi r' m/(2 (m - 1))) with
m = 6. In every instance the first 100 certified steps are compared with the exact selector and
distinguished test computed from the multiplied element (last column: failures of checks). Rates are per
walk step: exits, level-1 not certified, level-2 certified, exact evaluations, reports.

| Field | m | window | instances | model | mean steps at the stop / model | standard error | largest | correct outputs | wrong outputs | runs at the toy exit cap | exits | level 1 not certified | level 2 certified | exact evaluations | reports | selector check failures |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| bn | 6 | 4 | 1200 | 13,680 | 1.0191 | 0.0150 | 2.96 | 1200 | 0 | 0 | 0.2561 | 0.02363 | 0.00784 | 0.01579 | 0.01579 | 0 of 120,000 |
| bn | 6 | 8 | 1200 | 13,680 | 1.0365 | 0.0157 | 3.21 | 1200 | 0 | 0 | 0.1322 | 0.02360 | 0.00776 | 0.01584 | 0.01583 | 0 of 120,000 |
| generic | 6 | 4 | 1200 | 11,247 | 1.0250 | 0.0150 | 2.90 | 1200 | 0 | 0 | 0.2562 | 0.02248 | 0.00664 | 0.01584 | 0.01584 | 0 of 120,000 |
| generic | 6 | 8 | 1200 | 11,247 | 1.0194 | 0.0148 | 3.18 | 1200 | 0 | 0 | 0.1322 | 0.02238 | 0.00665 | 0.01572 | 0.01572 | 0 of 120,000 |

Pooled over 4,800 runs: mean 1.0250, standard error 0.0076, largest 3.21.

Negative control, the same runs over sqrt(pi r'/2): mean 1.1228.

**Census of the read coordinate over a whole toy GT** (`toy_offline.json`, key `census`; mode `toy` repeats
it for q = p(-7)). Every one of the r' = 1,704,961 elements of the order-r' subgroup for q = p(-15) =
1,706,311 was enumerated and read at coordinate 1 (70974bf7's) and at coordinate 6 (this package's):

| Coordinate | distinct values of t | selector counts (expected r'/6 = 284,160 each) | distinguished, D = 4, 6, 8, 10 (expected r'/2^D = 106,560, 26,640, 6,660, 1,665) |
| --- | ---: | --- | --- |
| 1 | 261,789 (at most 1 + (r' - 1)/6 = 284,161) | 284,635, 283,302, 284,544, 284,544, 283,302, 284,634 | 106,722, 26,580, 6,768, 1,638 |
| 6 | 1,078,147 (a random map of r' points to q values: 1,078,097) | 285,264, 283,993, 283,224, 283,224, 283,993, 285,263 | 106,369, 26,636, 6,530, 1,616 |

Coordinate 1 shows the six-element Frobenius orbits collapsed to one value; coordinate 6 (scaled by zeta^i
on an orbit) has as many distinct values as a random map, and its selector and distinguished frequencies are within 0.4% and two standard deviations of the
uniform ones (the selector counts come in equal pairs s, 5 - s because t(Y^-1) = q - t(Y)).

How to read the toy runs:
- *Stopping time.* The pooled mean plus 5 standard errors is below the allowance 107/100 of H1.
- *Negative control.* Measured against sqrt(pi r'/2), the model without the adding-walk factor, the same
  runs average above the allowance by more than 5 standard errors (last line under the table; the
  experiment checks this). So the factor m/(m - 1) is needed, and with it the allowance holds. This is
  the same set of runs read against another model, not a second experiment.
- *Every run ends with the correct x.* No run reached the toy exit cap and no output was wrong.
- *Selector.* No certified step, at either level, disagreed with the exact selector, and none was a
  distinguished element, in all checked steps: the corollary of section 2 at work on real walks.
- *Exits per step* are below 1/window plus the exact-evaluation rate, the count of section 5.
- *Level rates.* Level 1 fails on 2.2 to 2.4% of the steps. This includes the toy's distinguished elements
  (2^-6 = 0.0156 of the elements), which no level may certify and which go to the exact evaluation; level 2
  certifies the rest (0.66 to 0.78% of the steps). The exact-evaluation and report columns agree to within
  0.00001: level 2 almost never failed on an element that was not distinguished, as H2 predicts for the toy
  thresholds.
- *Production-size rates* (mode `algebra`): on 2,000 steps of a production-size walk, level 1 failed on
  0.6% of the steps (bound q1 = 2.09%), level 2 certified all of them, and the six selector counts give a
  chi-square of 3.45 on 5 degrees of freedom.
- *Limits:* the runs support H1 and H2 only at toy scale. Extrapolation to r near 2^254 is the heuristic
  itself.

## 8. Evidence and reproduction

`experiment.py` (standard library, with `kern.py`, `toy.py` and `outstep.py`) is the one declared experiment.
- **kernels:** generates the priced programs, executes them on the interpreter and checks them against
  big-integer arithmetic: the field kernels, the Fp12 product, the step loop on four real windows (one with
  a level-2 step) with both lemmas, selector, exit and registers checked at each step, the re-base with
  its digits, window constants and counters on both sides of both caps, both sides of both certificate
  thresholds, 240 entries with random words on the level-2 path, the trip budget, and that the six
  callee-saved registers holding digits and sp come back unchanged.
- **algebra:** checks p, r and g, the row identity, the signed recoding, both lemmas and the corollary on
  60 production-size pairs including extreme digits, the structure constants, the read coordinate under
  the Frobenius powers (the factor zeta^i on two random elements of GT), and the level rates and selector counts on a 2,000-step production-size walk.
- **ledger:** recomputes every integer in section 5, E_steps, E_level2, X_cap, K (and that K - 1 fails), the
  failure budget and the score.
- **stated figures and output step** (every mode): the output step is executed for 1, r - 1 and a seeded
  x; every figure and table row of the claim, this file and `machine.md` is
  compared with the computed value; the pooled toy mean plus 5 standard errors must be within the
  allowance of H1 and the negative control outside it. Any mismatch makes the experiment exit with
  status 1.
- **toy:** runs the algorithm to its output in the BN-family toy q = p(-15), r' = r(-15), with 8 to 32
  instances (32 initially) and windows of 4 or 6, compares the element sequence of the tag walk with a
  plain reference walk (full products and the exact selector, same coins) on the first 4 instances, and
  repeats the census for q = p(-7). At 32 instances the standard error of the mean is about 0.1, so one
  run does not test the allowance. The decisive statistics are the offline runs above.

Offline, the table above was produced with the same `TagRho` code (`toy.offline_instance`) and the seeds
stated in `toy_offline.json`, which holds the step count of every instance.

**Limitations.**
- The score is an expectation. An opened run costs more than the score: 2^120.183957 units in
  expectation and up to 2^131.353421 units.
- The success bound and the expectation are heuristic (H1, H2), and the revisit bound of section 5 uses H2.
- The table needs 2^61.353586 bytes; nothing of that size can be built today.
- The instruction counts are executed counts of generated programs on an interpreter, not timings. The
  table construction, the store, the exact evaluation and the pairing check are priced by itemised
  bounds and are not executed.
- Nothing runs at production size.
- Finite-field sieve methods are expected to be cheaper and are not used here.
