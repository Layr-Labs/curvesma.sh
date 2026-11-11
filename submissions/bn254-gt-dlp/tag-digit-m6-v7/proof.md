# Gated digit-tag 6-adding rho in BN254 GT: 144 instructions per certified step, window 1280

**Claim.** For uniform z = g^x in the protected BN254 GT, the algorithm below returns x with g^x = z
with probability at least 0.39 over uniform challenges and fresh coins, at **expected** total work

    E[W] <= 1 + (K/2^64) (144 E[T] + 9,216 (E[T]/1280 + E[walks] + 1) + 1,792 E[X] + 2^110) / 20,000 < 2^118.98

reference units (W = R/20,000, R the total scalar AArch64 instruction-equivalent work). Here T is the number
of walk steps and X the number of exact evaluations of a run past the gate. This is an expectation over
the challenge and the coins of the one attack on the one target, **not a bound on every run**. Every run is
capped at 2^132.96 units (section 5). Peak storage is below 2^63 bytes, preprocessing (the table,
included) is below 2^73 instructions, and there is no advice. This is a generic square-root attack, not
index calculus or NFS, not a practical break and not a security estimate. The success and expected-work
bounds use disclosed heuristics (H1, H2r, H-mean); the tag certificate and the 1/r degeneracy are proved.

**Lineage and what is new.**
- **Certified-tag method.** The certified-tag walk method and the table of all products of at most ell
  multipliers are those of submission 73463878-0c19-4309-a991-749ea525d51c (e-man07, after Cheon, Hong and
  Kim, J. Cryptology 25 (2012)). Its `reference.py`, `experiment.py`, `ttstep.c`, `ttstep.disasm.txt`,
  `tower.c`, `tower.disasm.txt` and `development-evidence.json` are included byte for byte except DP_BITS = 84 in `experiment.py` (its full-width walk and predicate check) and the dbits=84
  docstring in `reference.py` (the base used 80). They remain the reference implementation of the
  128-bit-word tag; their walk selects on the w^0 coefficient, and this package selects on another one.
- **Digit tag and 6-adding configuration.** Holding the window base as 72 digits of 43 bits, with one
  64-bit entry word per digit, and the 6-adding walk with long windows, are those of submission
  70974bf7-180b-47e1-a09b-c88363cfa805 (kaankolcu), adopted with credit.
- **Gate.** The run gate and the run-to-collision expected-work accounting follow submission
  8b901b3d-a617-4b1f-a442-ff370e6ec4f7 (terrapinelf).
- **Loop form.** The twice-unrolled loop without a per-trip MOV follows f86cff50-80c5-43cf-a773-67417a7c29e4
  (e-man07).
- **Window product.** The window re-base uses this account's generated, encoding-validated
  lazy-reduction Fp12 product of 8,679 instructions (accepted in e3771b98 and a698079e).

New here:
(i) **the selector coordinate t = Re(w^1 coefficient)**, which no nontrivial Frobenius power fixes, in
place of the w^0 coefficient, which p^2-Frobenius fixes. This answers objection (2) of the review of
70974bf7's pairing-inversion twin 271f78c0 (section 6);
(ii) **a generated digit loop of 144** executed instructions per certified trip (`walks.py` part S): 18
resident digits, the constant 6 in a register, and a 4-instruction certificate test plus a 1-instruction
selector for 6 multipliers;
(iii) **exact evaluation from a stored row.** Every entry stores the exact row v(T_e), and an uncertified
tag is resolved by a generated 1,294-instruction routine computing t = sum_k y_k v_k(T_e), with no
re-base. Exact evaluation is charged 1,792, at rate (14 2^48 + 2)/2^64 = 2^-12.19 per step;
(iv) **compact entries.** Non-leaf entries hold 72 words, 6 pointers, the row and the labels (1,072 bytes);
only leaves hold T_e's 54 product operands. The window length is therefore **1280** within 2^63 bytes;
(v) **a generated 171-instruction digit conversion** at each re-base;
(vi) **new toy evidence for the 6-adding walk**, with the mean allowance set to 1.08 above the measured
large-n ratios (section 6).

A step costs 144 + 9,216/1280 + 1,792 x (14 2^48 + 2)/2^64 <= **151.59** instructions in expectation. `lazy.py`
part G recomputes every figure below.

Notation: p, r, g are the protected modulus, group order and generator; z = g^x is the input. An element
Y of Fp12 has twelve Fp coordinates y_0..y_11 in the tower layout of `lazy.py`: coordinate 2m + b is the
b-th Fp part of the w^m coefficient over Fp2, with Fp12 = Fp2[w]/(w^6 - xi). t(Y) in [0, p) is coordinate
2, the real part of the w^1 coefficient: a fixed Fp-linear functional.

## 1. The walk as a function of the group element

Fix six multipliers M_s = g^(c_s) z^(d_s), s = 0..5. Define

    s(Y) = floor(6*t(Y)/p) in {0..5},      f(Y) = Y*M_{s(Y)},
    Y is distinguished  iff  t(Y)*2^84 < p.

f, s and the distinguished predicate depend on the group element only. The algorithm evaluates exactly
these functions; the table machinery is only a way to evaluate s(Y) cheaply. Hence the walk is a standard
6-adding walk with a deterministic iteration function, and two trails that reach the same element
coincide from there on.

## 2. Digit tag certificate (exact statement, proved)

**Row.** For a fixed T in Fp12, Y -> t(Y*T) is Fp-linear: t(Y*T) = sum_k y_k v_k(T) mod p (`lazy.py` `row`).

**Stored data.** A coordinate is kept as its Montgomery residue Y~_k = y_k 2^256 mod p (canonical). It is
split into six 43-bit digits d_(k,i) = bits [43i, 43i + 43) of Y~_k, so Y~_k = sum_i d_(k,i) 2^(43i). Here
d_(k,5) < 2^39 since Y~_k < p < 2^254. For a table product T put u_(k,i) = 2^(43i) v_k(T) 2^-256 mod p and
store the 64-bit words w_(k,i) = floor(2^64 u_(k,i)/p).

**Tag.** A = sum over the 72 pairs of d_(k,i) w_(k,i) mod 2^64.

**Lemma.** With t = t(Y*T), the real number tau = 2^64 t/p satisfies (tau - A) mod 2^64 in [0, E), E = 72 2^43.

Proof. sum_i d_(k,i) u_(k,i) = Y~_k v_k 2^-256 = y_k v_k (mod p), so t/p = sum_(k,i) d_(k,i) u_(k,i)/p
(mod 1). For one pair write d = d_(k,i) < 2^43 and u/p = w/2^64 + eps, 0 <= eps < 2^-64. Then
d u/p = d w/2^64 + d eps. Modulo 1 the first term equals (d w mod 2^64)/2^64 exactly; the second lies in
[0, 2^43 2^-64). Summing 72 pairs, 2^64 t/p - A mod 2^64 lies in [0, 72 2^43). QED. Part T checks this on
window bases produced by the generated product and digit conversion. Part I checks it at every step of a
full-width walk (largest observed offset 19.296 x 2^43 in the declared run).

**Corollary (certified step, 6 multipliers).** Let L = 6A mod 2^64 and THR6 = 2^64 - 14 2^48. Call A
*certified* if 1 <= L <= THR6. Since 6E = 13.5 2^48 < 14 2^48, the interval [L, L + 6E) does not reach
2^64. So 6 tau in [6A, 6A + 6E) stays within one multiple of 2^64, floor(6 tau/2^64) = floor(6A/2^64),
and s(Y*T) = UMULH(A, 6) exactly. This also excludes a wrap of [A, A + E) past 2^64, which would put L in
(THR6, 2^64). L >= 1 excludes A = 0, so tau >= A >= 1, i.e. t/p >= 2^-64 > 2^-84: the element is not
distinguished. In the loop: `mul x9, x3, x17` (L, x17 = 6), `sub x9, x9, #1`, `cmp x9, x4`, `b.hs EXACT`
(x4 = THR6, one MOVN), and `umulh x9, x3, x17` for the selector.

## 3. Complete algorithm

Parameters: R = 6 multipliers; window length 1280; distinguished points as defined in section 1
(t(Y) 2^84 < p, density 2^-84 under H2r); step cap S_cap = 2^140 steps in total; at most 2^52 walks, each at
most 2^90 steps; at most 2^127 exact evaluations; at most 2^70 store probes. Reaching any cap ends the run
as a failure.

0. **Gate.** Draw a fresh uniform 64-bit word v. Unless v < K = 7237656125499723472 (the least K with
   (K/2^64) x 0.994 >= 0.39), stop and fail (below one unit of work).
1. Decode z, check z != 0 and z^r = 1. If z = 1 output 0.
2. Draw c_s, d_s uniformly mod r (rejection sampling of 254-bit words with at most 1,024 draws per scalar;
   exhausting them ends the run as a failure, probability below 2^-200 per scalar, charged in section 6)
   and form M_s = g^(c_s) z^(d_s), s = 0..5. Build the table: one entry e for every multiset of at most 1280
   multipliers (C(1286,6) = 6,209,268,427,492,417 entries). Each entry holds the 72 words of section 2 for
   T_e, the six child pointers, the exact row v(T_e) (12 Montgomery residues) and alpha_e = sum c_s,
   beta_e = sum d_s mod r. The C(1285,5) = 28,970,148,184,257 leaves (depth 1280) also hold T_e as its 54
   canonical Karatsuba operands (the form the generated product reads). Entries are addressed by
   P = entry + 512 (`machine.md`, table entry).
3. For each walk: draw a, b uniformly mod r (same bounded sampling), set Y0 = g^a z^b (skip if 1), convert
   Y0 to its 72 digits, e = root. Repeat:
   a. A = tag(Y0, e) (one trip of the generated loop).
   b. If A is certified: s = UMULH(A, 6). If e has a child for s, set e to it; otherwise (window end, a
      leaf) re-base: Y0 <- Y0 * T_e by the generated product, digits by the generated conversion,
      (a, b) <- (a + alpha_e, b + beta_e), e <- the root's child for s. This is one step Y <- Y * M_s.
   c. If A is not certified (an *exact evaluation*): compute t = t(Y0 T_e) = sum_k y_k v_k(T_e) exactly by
      the generated routine from Y0's residues and the stored row, with no re-base, and evaluate s(Y)
      and the distinguished predicate by their definitions. Leave the loop at a distinguished point (the
      labels are (a + alpha_e, b + beta_e)). Otherwise take the step as in b with this s.
4. At a distinguished point look t up in the store (2^54 slots, open addressing, **keyed by t**, the
   selector coordinate, not by the whole element). If an entry (t, a', b') exists, put x = (a - a')/(b' - b)
   mod r (failure if b = b'). If g^x = z output x; otherwise treat it as a false match. Insert (t, a, b) and
   start the next walk.

There is no fixed step budget: past the gate the run continues until it outputs x or a cap fires. The
global step count is a three-word counter of the driver; the per-walk length and the exact-evaluation
meter are two-word counters (`machine.md`); no counter wraps.

## 4. Correctness

By the corollary, step 3b applies exactly f of section 1 and never skips a distinguished point; step 3c
uses the definitions directly. So every walk is a trail of f from its start, stopped at its first
distinguished point, whatever representation (Y0, e) the element happens to have.

The base package's `experiment.py` (changed only to the 84-bit predicate, section 7) drives its
reference `solve` loop, which selects on its own coordinate, at full BN254 width for 120 steps
(`steps_wanted` = 120; its output field `steps` counts 121 traced positions). Against independently
multiplied elements, it checks at every step that the selector is exact and that consecutive elements
differ by exactly the selected multiplier, and (every 40 steps; 4 label checks in the declared run) that
Y = g^a z^b for the represented label. It repeats the same walk with a deliberately narrow certificate (54 protected bits, 14 exact evaluations) and obtains the identical element sequence.

The priced implementation of this package is checked the same way (`lazy.py` part I). A full-width walk of
120 steps (window length 8 in the test; production 1280) uses 6 real multipliers g^c z^d and real table
entries: 72 words from the row of t(Y T_e) with t = Re of the w^1 coefficient, and child pointers.
- Every step is one trip of the generated 144-instruction digit loop. Every window end is the generated
  product followed by the generated digit conversion. At every 10th step the generated exact evaluation is
  also run on the stored row and compared with the true t.
- At every step the tag equals its definition, its error is inside [0, 72 2^43), the selector equals
  floor(6 t(Y)/p), and the loop's exit and next entry match the definition.
- Each re-based Y0 equals Y0 T_e, and its digits equal the definition.
- Every 100 steps the label invariant Y = g^a z^b and the step Y -> Y M_s are checked.
- Results (declared run: 120 steps, 14 window ends, 1 label check, 0 mismatches, 0 product failures).
- The generated exact evaluation agreed with the true t in all 12 of its checks.

**Labels and recovery.** Y = g^a z^b holds at each start and is preserved by multiplying with M_s and
adding (c_s, d_s). Two distinguished points with equal t that are the same element give
g^a z^b = g^a' z^b', so x = (a - a')/(b' - b) whenever b != b' (r is prime). The check g^x = z makes a
wrong output impossible. Because the store is keyed by t, a hit can also be a different element with the
same t. Such a false match fails the check, costs one exponentiation, is metered in the fixed allowance,
and the walk continues with the next start; it never produces an output. Under H2r its probability per
stored pair is about 1/p, so with at most 2^52 entries it is below 2^104/2^253 = 2^-149.

**Degeneracy** is exactly 1/r, without heuristics. Fix x and condition on the logarithms m_s = c_s + x d_s
of the multipliers and s_i = a_i + x b_i of all walk starts. Every group element, tag, exact evaluation,
distinguished point and store hit is a function of these. Conditionally the b_i and d_s remain
independent and uniform. For a hit between walk i and an earlier walk j, b - b' = b_i - b_j + (an integer
combination of the d_s), which is uniform on Z/r because b_i is; it vanishes with probability exactly 1/r.

## 5. Work

Per-item prices (`machine.md`). The loop, the product, the digit conversion and the exact evaluation are
generated, encoding-validated and executed by `lazy.py`/`walks.py`; `fr_add` is recounted from
`tower.disasm.txt`.

| item | instructions |
|---|---:|
| one certified step (one trip of the generated digit loop: 71 MADD + 1 MUL, 63 LDP, 4-instruction certificate test, UMULH selector, child load, exit test, counter) | 144 |
| one window end: re-base 8,938 (generated product 8,679 + digit conversion 171 + two `fr_add` 68 + glue 20), loop entry 20 + exit 13, driver with the three-word counter <= 160 | 9,131 <= **9,216** |
| one exact evaluation: generated routine 1,294, selector and distinguished test, aborted trip, loop exit and re-entry, driver, meter (itemised in machine.md) | <= **1,792** |
| everything else, by meters and caps: table build, <= 2^52 walk starts, store, recoveries | <= 2^110 |

A window end occurs after 1280 steps since the previous re-base, or at a walk start, so the number of
window ends is at most T/1280 + walks + 1.

**Every run** (independent of all heuristics): T <= 2^140, walks <= 2^52 and exact evaluations <= 2^127, so
R <= 144 x 2^140 + 9,216 (2^140/1280 + 2^52 + 1) + 1,792 x 2^127 + 2^110 and W <= 2^132.96 units.

**Expected work.** Let T_touch be the number of steps (all walks together) until some step reaches an
element visited by an earlier trail. After it, the touching walk reaches its distinguished point within its
remaining length: at most 2^90 steps, 2^84 in expectation under H2r. Hence E[T] <= E[T_touch] + 2^85. The
number of walks is at most one more than the number of distinguished points reached, which is T/2^84 in
expectation under H2r (charged twice: E[walks] <= 2 E[T]/2^84 + 1). A step's tag is uncertified, with L
outside [1, THR6], with probability (14 2^48 + 2)/2^64 = 2^-12.19 under H2r, so
E[X] <= E[T] (14 2^48 + 2)/2^64. Every certified step costs the fixed 144, so

    E[R_run] <= 144 E[T] + 9,216 (E[T]/1280 + E[walks] + 1) + 1,792 E[X] + 2^110.

**H-mean.** E[T_touch] <= 1.08 sqrt(pi kappa r / 2), kappa = 6/5 (the 6-adding correction of H1). Then
E[T] <= 2^127.3667, E[X] <= 2^115.18, E[R_run]/20,000 <= 2^120.3229 units. With the gate,
W <= 1 + 1[v < K] W_run with v independent of the run, so

    E[W] <= 1 + (K/2^64) x 2^120.3229 <= 2^118.9732.

Submitted bound: **time_log2 = 118.98** (rational 5949/50); the computed value is 118.97312057..., so the
claim is at least 0.0068 bits above it. Per step this is 144 + 7.2 + 0.39 <= 151.59 instructions in
expectation.

## 6. Success probability: heuristics, support and limits

Past the gate the run fails only if a cap fires or a recovery degenerates. Under H1, H2r and H-mean
(`lazy.py` part G, `losses`):
* the step cap: P(T > 2^140) <= E[T]/2^140 <= 0.000158 (Markov; exact ledger value 0.00015739...);
* the walk cap: P(walks > 2^52) <= E[walks]/2^52 <= 0.005037 (Markov; exact ledger value 0.0050365...);
* the exact-evaluation meter: P(X > 2^127) <= E[X]/2^127 <= 0.000276 (Markov; exact ledger value
  0.00027543...);
* the store-probe meter (2^54 slots, at most 2^52 keys, linear probing, O(1) expected probes): below 2^-20;
* a walk not distinguishing within 2^90 = 64 x 2^84 steps, closing a cycle on itself first, or meeting
  the inverse of an element of its own trail, in expectation over E[walks] <= 2^44.4 walks: below 2^-40
  (see the revisitation argument below);
* degeneracy b = b' (exactly 1/r per hit, section 4), a false t-collision, scalar sampling failure, a start
  equal to 1: below 2^-140.

The sum is below 0.0055 (ledger: 0.0054702...), so success past the gate is at least 0.9945; the claim uses
0.994. With the gate, success is at least (K/2^64) x 0.994 >= 0.39 and K is the least such integer (part G:
`success_ok`, `K_is_least`). Once a touch has happened, no other event prevents success: the later trail
follows the earlier one to its stored distinguished point.

**H1** (collision law of the 6-adding walk with kappa = 6/5, the r-adding correction of Brent-Pollard and
Teske, Math. Comp. 70 (2001)) is used only through H-mean.

**H2r (restricted equidistribution), stated so as to answer the review of 271f78c0.** The selector
coordinate t is the real part of the w^1 coefficient. Under the p^k-Frobenius the w^1 coefficient becomes
conj^k(e_1) xi^((p^k - 1)/6), with a twist different from 1 for every k = 1..11. So no nontrivial
Frobenius power fixes t, and the only exact symmetry is inversion (k = 6), t(Y^-1) = p - t(Y). Part T
checks this on a GT element: the eleven Frobenius images have t different from t(Y), and the 6th has
p - t(Y).

The premise: over the walk's *first visits to inversion pairs {Y, Y^-1}*, the events used (selector value,
certificate failure with probability (14 2^48 + 2)/2^64, distinguished with probability 2^-84) occur at
their measure on Fp, as fresh events conditional on the history of first visits. No uniformity
conditional on revisits is assumed. The other lane of this account measured the selector on BN254 itself
for this coordinate: over a 3,000-step walk the selector counts were [505, 496, 505, 489, 506, 499]
(chi-square 0.45 on 5 dof). Kaankolcu's coordinate, which p^2-Frobenius fixes, gave chi-square 8.53 on
the same test. That check is not shipped here; Part T ships the Frobenius check.

**Revisitation-aware tail argument.** A walk's trajectory until its first distinguished point consists of
first visits to inversion pairs, except possibly for its first revisit. That revisit is one of three kinds:
- (a) an element of an earlier trail, a *touch*, after which the walk follows that trail to its stored
  distinguished point, which is success;
- (b) an element of its own trail, a *self-cycle*;
- (c) the inverse of an element of its own trail.

A non-distinguished element revisited is certainly not distinguished, and none of these revisits is
treated as a fresh event; only (b) and (c) are charged.
- Under H2r, the distinguished indicator of a first visit is fresh with probability 2^-84, so 2^90 first
  visits without one have probability at most exp(-64) = 2^-92.3.
- (b) and (c) within a walk of geometric length (mean 2^84, capped at 2^90) have probability at most
  about 2 E[len^2]/(2r) = 2^169/r <= 2^-84.6.
- The charged term is E[walks] (2^-92.3 + 2^-84.6) <= 2^44.4 x 2^-84.59 = 2^-40.19 < 2^-40.

**H-mean** (mean touch time <= 1.08 of the model) is a statement about the mean only. Support: the walk
itself at toy scale (`walks.py` part H with m = 6: fresh 61-bit prime field, prime n-bit subgroup, target
and coins per run; 6-adding walk with selector floor(6t/q); distinguished points; sequential walks to their
distinguished point; store keyed by the whole toy element; run to the first verified recovery; per-run
ratios to sqrt(pi kappa n/2), kappa = 6/5). Offline (`development-evidence-mean6.json`, eighteen processes,
seeds 20261301..20261318):

| n bits | dp bits | runs | mean touch / model | s.e. | mean steps to recovery / model |
|---:|---:|---:|---:|---:|---:|
| 16 | 3 | 3,000 | 0.923 | 0.009 | 0.967 |
| 20 | 4 | 2,000 | 0.966 | 0.011 | 0.990 |
| 24 | 5 | 2,000 | 0.984 | 0.011 | 0.996 |
| 28 | 6 | 1,800 | 1.009 | 0.012 | 1.014 |
| 32 | 7 | 1,800 | 1.015 | 0.013 | 1.018 |
| 36 | 8 | 240 | 1.020 | 0.032 | 1.021 |

Pooled at n >= 28 bits: 3,840 runs, mean 1.0129, s.e. 0.0085; mean + 3 s.e. = 1.0384.

**Limit, disclosed.** The ratio rises with n: from 0.92 at 16 bits to about 1.01-1.02 at 28-36 bits.
The increments shrink: +0.043, +0.018, +0.025, +0.006 and +0.005 per 4 bits. That is consistent with
convergence to a constant near 1.02, the r-adding walk's known slowdown beyond the r/(r-1) correction, but
it is not a proof of it. The charged allowance 1.08 is 0.06 above the large-n data and about 5 standard
errors above the pooled large-n mean. The premise is that the ratio stays below 1.08 at n = 2^254, an
extrapolation from groups of at most 36 bits. The toy groups have one coordinate; in production the
selector reads t, one Fp-linear functional of an element of a 3048-bit field.

## 7. Experiments

The declared experiment is `lazy.py` (standard library). Each part's time is reported in its output
(`part_seconds`), and a guard skips the base experiment if 11 s have elapsed. Parts H, I and S are in
`walks.py`, which `lazy.py` executes in its own namespace to keep each file under the 64 KiB limit.
- **Part C** reads `submission.json`, `proof.md` and `machine.md` and asserts that they agree with the code
  and the ledger:
  - the 84-bit predicate is the only distinguished-point predicate stated (no 80-bit density statement
    anywhere);
  - the code constants are 84 (`DP_BITS` in `lazy.py` and `experiment.py`, the `reference.py` docstring);
  - every parameter (144, 151.59, 8,679, 9,216, 1,792, 2^140, 2^52, 2^90, 2^54, 2^127, 1280, K) appears
    where it is used, and every itemised total is at least the sum of its items;
  - every rounded upper bound is at least, and every lower bound at most, the computed value;
  - the manifest score equals the ledger claim.

  (49 checks; `all_consistent`)
- **Part S** generates the digit loop, validates the encoding of every instruction (0 violations) and
  executes it on 40 synthetic entries from 43
  starts (1,878 trips). The entries have random 64-bit words, 43-bit digits, a random child structure with
  window ends, and some entries forced to fail the certificate. Every trip's A, selector, child and exit
  agree with the definition, every certified in-window trip executes 144 instructions (both trips), and
  the BUDGET, REBASE and EXACT exits all occur.
- **Part L** generates the window product, validates the encoding of every instruction (0 violations) and
  executes it (20 cases: random, edge-valued and all-limbs-large inputs, 0 failures, canonical outputs,
  8,679 instructions).
- **Part T** checks the digit tag lemma on 6 window bases produced by the generated product and digit
  conversion: 0 offsets outside [0, 72 2^43), digits equal to the definition, and the exact selector on
  every certified tag. It also runs the Frobenius check of section 6.
- **Part H** is the toy mean of the 6-adding walk. Declared input: up to 300 / 150 / 37 runs at 16 / 20 /
  24 bits within a budget of at most 1 s, split equally across the three settings and shortened if needed
  so that the whole run stays inside the 20 s limit. The output reports the budget used and the runs
  completed. The coin seed is 20261113, different from the offline seeds, and the samples are not pooled.
- **Part I** is the integrated full-width walk of section 4.
- **Part G** is the ledger.
- **Part X** runs last, for the declared input only. It reads this run's own output and the evidence file,
  and asserts that every experimental figure stated in `submission.json`, `proof.md` and `machine.md`
  equals the measured value.

Then `lazy.py` executes the base package's `experiment.py` (84-bit predicate) in-process with
`{"seed": 20261111, "trials": 16}`. That experiment covers tower checks, the row identity, its 128-bit-word
certificate on 260 full-width cases (largest observed offset 41 in the declared run), the 120-step
full-width walk through its reference `solve` loop (23 re-bases), toy runs at the cap (16 / 16 / 16 trials)
and the recount of its reference listings (273 per trip; `fr_add` 34). Its printed score (121.285447) is the
base package's hard-capped ledger, superseded here by part G.

References: Pollard 1978; Teske, Math. Comp. 70 (2001); Cheon, Hong and Kim, J. Cryptology 25 (2012); van
Oorschot and Wiener, J. Cryptology 12 (1999); Aranha, Karabina, Longa, Gebotys and Lopez, EUROCRYPT 2011
(eprint 2010/526); submissions 73463878 and f86cff50 (e-man07), 70974bf7 and 271f78c0 (kaankolcu),
8b901b3d and e2dfd4cf (terrapinelf), 67df70ff, e3771b98, a698079e.
