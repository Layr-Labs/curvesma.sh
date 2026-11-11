# Gated tag-tracing 8-adding rho in BN254 GT with a generated step loop and lazy window product

**Claim.** For uniform z = g^x in the protected BN254 GT, the algorithm below returns x with g^x = z
with probability at least 0.39 over uniform challenges and fresh coins, at **expected** total work

    E[W] <= 1 + (K/2^64) (215 E[T] + 9,216 (E[T]/312 + E[walks] + 1) + 2^110) / 20,000 < 2^119.59

reference units (W = R/20,000, R the total scalar AArch64 instruction-equivalent work), where T is the number
of walk steps of a run past the gate. This is an expectation over the challenge and the coins of the one
attack on the one target, **not a bound on every run**; every run is capped at 2^132.65 units (section 5).
Peak storage is below 2^63 bytes; preprocessing (the table, included) below 2^72 instructions; no advice.
This is a generic square-root attack, not index calculus or NFS, not a practical break and not a security
estimate. The success and expected-work bounds use disclosed heuristics (H1, H2, H-mean); the tag
certificate and the 1/r degeneracy are proved.

**Lineage and what is new.** This package is e-man07's revision of kshitij-hash's accepted submission
a698079e (`tag-gated-lazy-v4`, 119.61); every file is carried from it, and only the step loop (215
instead of 217 instructions, from e-man07's accepted f86cff50-80c5-43cf-a773-67417a7c29e4), the
step-cap value (2^139 instead of 2^136) and the charged success past the gate (0.994 instead of 0.99)
change. In a698079e: the walk, tag certificate and table are those of
submission 73463878-0c19-4309-a991-749ea525d51c (e-man07, tag tracing after Cheon, Hong and Kim, J. Cryptology
25 (2012); accepted at 121.285447 on bn254-pairing-inversion). Its `reference.py`, `experiment.py`,
`ttstep.c`, `ttstep.disasm.txt`, `tower.c`, `tower.disasm.txt` and `development-evidence.json` are included
byte for byte except DP_BITS = 84 in `experiment.py` (its full-width walk and predicate check) and the
dbits=84 docstring in `reference.py` (the base used 80). The step-loop
listing is kept as the reference implementation and is still recounted by that experiment. The run gate and the run-to-collision expected-work accounting follow submission
8b901b3d-a617-4b1f-a442-ff370e6ec4f7 (terrapinelf; the same accounting was accepted in its
pairing-inversion twin e2dfd4cf). New here:
(i) **a generated step loop** (`lazy.py` part S; kshitij-hash's generator as revised in f86cff50) of **215**
executed instructions per certified trip (a698079e: 217; base listing: 273). It computes the same
A = sum_j (UMULH(y_j, Flo_j) + y_j Fhi_j) mod 2^64 with 3 arithmetic instructions per limb (UMULH, MADD, ADD),
one LDP per (Flo, Fhi) pair, 19 limbs resident in registers and the other 29 reloaded by LDP/LDR. It is
unrolled twice so that the entry pointer alternates between x0 and x10 (no MOV per trip), the UMULH
result overwrites the dead Flo register, freeing a register that holds 768, and the certificate test is
`add x9, x4, x3, lsl #3; cmp x9, #768; b.ls` (3 instructions instead of 4), the same test since
768 + 8A = 8 ((A + 96) mod 2^61) mod 2^64. The loop is encoding-validated and executed against the
definition at every trip, and all exits are exercised;
(ii) the window re-base uses a **generated lazy-reduction Fp12 product** of 8,679 executed instructions
(`lazy.py` part L; kshitij-hash's corrected version of the routine of their rejected submission 67df70ff:
encodable frame, LDP operand loads, encoding validator) instead of the 24,044-instruction `fp12_mul` of
`tower.c`. A window end costs at most 9,216 instead of 24,576. With the window length raised from 256 to
ell = 312 (multisets of at most 312 multipliers, C(320,8) table entries; the larger proven multiset bound
follows submission 2da98c02, which raised it to 312 for the pairing-inversion twin), a step costs
215 + 9,216/312 <= **244.54** instructions instead of 369 (a698079e: 246.54). This is the same window routine as kshitij-hash's submission f990c954
(119.97), which used the base loop;
(iii) the run is **gated and taken to its first recovery**, scored at its expected work under the premise
H-mean, which is measured here on the walk itself (`lazy.py` part H; `development-evidence-mean.json`);
(iv) the track asks for x, so there is no output map.
`lazy.py` part G recomputes every figure below.

Notation: p, r, g are the protected modulus, group order and generator;
z=g^x is the input. An element Y of Fp12 has twelve Fp coordinates
y_0..y_11 in a fixed basis (the tower order of `tower.c`); t(Y) in [0,p)
denotes its coefficient of w^0 in the protected basis, a fixed Fp-linear
functional of the coordinates.

## 1. The walk as a function of the group element

Fix eight multipliers M_s=g^(c_s)*z^(d_s), s=0..7. Define

    s(Y) = floor(8*t(Y)/p) in {0..7},      f(Y) = Y*M_{s(Y)},
    Y is distinguished  iff  t(Y)*2^84 < p.

f, s and the distinguished predicate depend on the group element only. The
algorithm below evaluates exactly these functions; the table machinery is
only a way to evaluate s(Y) cheaply. Hence the walk is a standard 8-adding
walk with a deterministic iteration function, and two trails that reach
the same element coincide from there on.

## 2. Tag certificate (exact statement, proved)

Row. For a fixed T in Fp12, Y -> t(Y*T) is Fp-linear, so there is a row
v(T) in Fp^12 with t(Y*T) = sum_k y_k*v_k mod p (explicit formula in
`reference.py`, `BN254GT.row`; checked on 40 random pairs).

Stored data. Coordinates are kept as 4-limb Montgomery residues
Y~_k = y_k*2^256 mod p with limbs Y~_{k,i}, i=0..3. For T put
u_{k,i} = v_k(T)*2^(64i-256) mod p and store the 128-bit integers
F_{k,i} = floor(2^128*u_{k,i}/p) = Fhi*2^64+Flo.

Tag. A = sum over the 48 pairs (k,i) of
( UMULH(Y~_{k,i},Flo_{k,i}) + Y~_{k,i}*Fhi_{k,i} ) mod 2^64.

**Lemma.** With t=t(Y*T), the real number tau = 2^64*t/p satisfies
(tau - A) mod 2^64 in [0,96).

Proof. sum_i Y~_{k,i}*u_{k,i} = Y~_k*v_k*2^-256 = y_k*v_k (mod p), so
t/p = sum_{k,i} Y~_{k,i}*u_{k,i}/p (mod 1). For one pair write y=Y~_{k,i}
<2^64 and u/p = F/2^128+eps, 0<=eps<2^-128. Then
y*u/p = y*Fhi/2^64 + y*Flo/2^128 + y*eps. Modulo 1 the first term equals
(y*Fhi mod 2^64)/2^64 exactly; the second equals UMULH(y,Flo)/2^64 plus a
remainder in [0,2^-64); the third lies in [0,2^-64). So each pair
contributes its two computed words plus an error in [0,2*2^-64), and the
48 errors add up to less than 96*2^-64. Multiply by 2^64. QED
(`experiment.py` checks the lemma on 260 full-width cases including edge
values; largest observed offset 41 in the declared run.)

**Corollary (certified step).** Call A *certified* if
((A+96) mod 2^64) mod 2^61 > 96. Then the low 61 bits of A lie in
[1, 2^61-97], so [A,A+96) contains no multiple of 2^61 and does not wrap:
tau and A have the same top three bits, i.e. s(Y*T) = A>>61 exactly, and
tau>=A>=1, i.e. t/p>=2^-64>2^-84: the element is not distinguished.

## 3. Complete algorithm

Parameters: R = 8 multipliers; window length 312; distinguished points as defined in section 1
(t(Y) 2^84 < p, density 2^-84 under H2; this package uses 84 bits everywhere, the base package used 80);
step cap S_cap = 2^139
steps in total; at most 2^52 walks, each at most 2^90 steps; at most 2^90 exact evaluations; at most
2^70 store probes. Reaching any cap ends the run as a failure.

0. **Gate.** Draw a fresh uniform 64-bit word v. Unless v < K = 7237656125499723472 (the least K with
   (K/2^64) x 0.994 >= 0.39), stop and fail (below one unit of work).
1. Decode z, check z != 0 and z^r = 1. If z = 1 output 0.
2. Draw c_s, d_s uniformly mod r (rejection sampling of 254-bit words with at most 1,024 draws per scalar;
   exhausting them ends the run as a failure, probability below 2^-200 per scalar, charged in section 6)
   and form M_s = g^(c_s) z^(d_s), s = 0..7. Build the table: one entry e for
   every multiset of at most 312 multipliers (C(320,8) = 2,496,767,727,736,920 entries), holding the 96 words
   F(T_e) of section 2, the eight child pointers, alpha_e = sum c_s and beta_e = sum d_s mod r, and T_e stored
   as its 54 canonical Karatsuba operands (the 12 coordinates and their sums, `machine.md`), the form the
   generated product reads. Entries are addressed by P = entry + 512 (`machine.md`, table entry).
3. For each walk: draw a, b uniformly mod r (same bounded sampling), set Y0 = g^a z^b (skip if 1),
   e = root. Repeat:
   a. A = tag(Y0, e) (one trip of the generated loop).
   b. If A is certified: s = A >> 61. If e has a child for s, set e to it; otherwise (window end) re-base:
      Y0 <- Y0 * T_e by the generated routine, (a, b) <- (a + alpha_e, b + beta_e), e <- the root's child
      for s. This is one step Y <- Y * M_s.
   c. If A is not certified: re-base to e = root, read t = t(Y0) exactly and evaluate s(Y) and the
      distinguished predicate by their definitions; leave the loop at a distinguished point, otherwise
      set e to the root's child for s (one step).
4. At a distinguished point look t up in the store (2^54 slots, open addressing). If (t, a', b') exists,
   put x = (a - a')/(b' - b) mod r (failure if b = b'); if g^x = z output x, otherwise treat it as a
   false match. Insert (t, a, b); start the next walk.

There is no fixed step budget S: past the gate the run continues until it outputs x or a cap fires. The
global step count is a three-word (192-bit) counter of the driver, the per-walk length a two-word counter
(`machine.md`); 2^139 < 2^192, so neither wraps.

## 4. Correctness

By the corollary, step 3b applies exactly f of section 1 and never skips a
distinguished point; step 3c uses the definitions directly. So every walk
is a trail of f from its start, stopped at its first distinguished point,
whatever representation (Y0,e) the element happens to have.
`experiment.py` (base package; changed only to the 84-bit predicate, section 7) drives the real `solve` loop at full BN254 width for 120
steps (`steps_wanted` = 120; its output field `steps` counts 121 traced positions) and checks at every step, against independently multiplied
elements, that the selector equals floor(8t/p), that consecutive elements
differ by exactly the selected multiplier, and (every 40 steps; 4 label checks in the declared run) that
Y=g^a*z^b for the represented label; it repeats the same walk with a
deliberately narrow certificate (54 protected bits, 14 exact evaluations)
and obtains the identical element sequence.

The priced implementation is checked the same way (`lazy.py` part I). A full-width walk of 120 steps
(window length 8 in the test; production 312) uses 8 real multipliers g^c z^d and real table entries
(F words from the row of t(Y T_e) in this layout, child pointers). Every step is one trip of the
generated 215-instruction loop and every window end is the generated product. At every step the
selector equals floor(8 t(Y)/p) and the loop's exit and next entry match the definition. Each re-based
Y0 equals Y0 T_e. Every 100 steps the label invariant Y = g^a z^b and the step Y -> Y M_s are checked
(declared run: 120 steps, 14 window ends, 1 label check, 0 mismatches, 0 product failures).

Labels: Y=g^a*z^b holds at each start and is preserved by multiplying with
M_s and adding (c_s,d_s). Two distinguished points with equal t that are
the same element give g^a z^b = g^a' z^b', so x=(a-a')/(b'-b) whenever
b!=b' (r is prime). The check g^x=z makes a wrong output impossible; a
false match (equal t, different elements) fails it and costs one
exponentiation, metered in the fixed allowance. The output is the scalar x.

Degeneracy is exactly 1/r, without heuristics. Fix x and condition on the
logarithms m_s=c_s+x*d_s of the multipliers and s_i=a_i+x*b_i of all walk
starts. Every group element, tag, exact evaluation, distinguished point
and store hit is a function of these. Conditionally the b_i and d_s remain
independent and uniform. For a hit between walk i and an earlier walk j,
b-b' = b_i-b_j + (an integer combination of the d_s), which is uniform on
Z/r because b_i is; it vanishes with probability exactly 1/r.

## 5. Work

Per-item prices (`machine.md`; the step loop and the window product are generated, encoding-validated
and executed by `lazy.py` parts S and L; `fr_add` is recounted from `tower.disasm.txt` by `experiment.py`):

| item | instructions |
|---|---:|
| one certified step (one trip of the generated loop: 48 UMULH, 47 MADD + 1 MUL, 48 ADD, 62 LDP + 1 LDR, certificate test (3), pointer, exits, counter; unrolled twice, no MOV) | 215 |
| one window end: generated product 8,679 + two `fr_add` 68 + glue 20 (re-base 8,767), loop entry 19 + exit <= 13, driver with the three-word counter <= 160 | 8,959 <= **9,216** |
| everything else, by meters and caps: table build, <= 2^52 walk starts, <= 2^90 exact evaluations, store, recoveries | <= 2^110 |

A window end occurs after 312 steps since the previous re-base, at a walk start, or within an exact
evaluation (priced in that allowance), so the number of window ends is at most T/312 + walks + 1.

**Every run** (independent of all heuristics): T <= 2^139 and walks <= 2^52, so
R <= 215 x 2^139 + 9,216 (2^139/312 + 2^52 + 1) + 2^110 and W <= 2^132.65 units.

**Expected work.** Let T_touch be the number of steps (all walks together) until some step reaches an
element visited by an earlier trail. After it, the touching walk reaches its distinguished point within its
remaining length: at most 2^90 steps, 2^84 in expectation under H2. Hence E[T] <= E[T_touch] + 2^85.
The number of walks is at most one more than the number of distinguished points reached, which is T/2^84
in expectation under H2 (charged twice: E[walks] <= 2 E[T]/2^84 + 1). The cost of every step is the fixed
215 (the loop does not depend on the branch), so no Wald argument on step prices is needed:

    E[R_run] <= 215 E[T] + 9,216 (E[T]/312 + E[walks] + 1) + 2^110.

**H-mean.** E[T_touch] <= 1.05 sqrt(pi kappa r / 2), kappa = 8/7 (the 8-adding correction of H1). Then
E[T] <= 2^127.2909, E[R_run]/20,000 <= 2^120.9371 units, and with the gate, W <= 1 + 1[v < K] W_run with v
independent of the run, so

    E[W] <= 1 + (K/2^64) x 2^120.9371 <= 2^119.5873.

Submitted bound: **time_log2 = 119.59** (rational 11959/100); the computed value is 119.58723912..., so the
claim is at least 0.0027 bits above it. Against a698079e (119.60480778...) the changes contribute 0.0118
bits (loop 217 -> 215) and 0.0058 bits (gate K/2^64 = 0.39394 -> 0.39235). For
comparison, the base package's hard-capped run at S = sqrt(1.409 r) steps and 369 instructions per step is
2^121.285447 on every run.

## 6. Success probability: heuristics, support and limits

Past the gate the run fails only if a cap fires or a recovery degenerates. Under H1, H2 and H-mean
(`lazy.py` part G, `losses`):
* the step cap: P(T > 2^139) <= E[T]/2^139 <= 0.000299 (Markov; exact ledger value 0.00029866...);
* the walk cap: P(walks > 2^52) <= E[walks]/2^52 <= 0.004779 (Markov; exact ledger value 0.0047786...);
* the exact-evaluation meter: a step is uncertified with probability 97/2^61 under H2, so the expected
  number is below 2^72.9, so P(> 2^90) < 2^-17.1;
* the store-probe meter (2^54 slots, at most 2^52 keys, linear probing, O(1) expected probes): below 2^-20;
* a walk not distinguishing within 2^90 = 64 x 2^84 steps (exp(-64) per walk), or closing a cycle on itself
  first (about E[len^2]/(2r) = 2^168/2^253.6 = 2^-85.6 per walk), in expectation over E[walks] <= 2^44.3 walks:
  below 2^-41;
* degeneracy b = b' (exactly 1/r per hit, section 4), scalar sampling failure, a start equal to 1: below
  2^-200.
The sum is below 0.0051 (ledger: 0.0050853...), so success past the gate is at least 0.9949; the claim uses
0.994. With the
gate, success is at least (K/2^64) x 0.994 >= 0.39 and K is the least such integer (part G: `success_ok`,
`K_is_least`). Once a touch has happened no other event prevents success: the later trail follows the
earlier one to its stored distinguished point.

**H1** (collision law of the 8-adding walk with kappa = 8/7) and **H2** (uniform tag bits, and the leading
bits of t/p at distinguished scale) are those of the base package, with its support: the selector is the
classical leading-bits selector of finite-field rho and the walk is Teske's r-adding walk (Math. Comp. 70
(2001)); the base experiment runs the complete algorithm on toy prime-field subgroups at the cap (success
fractions around the model value 0.46).

**H-mean** is new and is a statement about the mean only. Support, from the walk itself (`lazy.py` part H:
fresh 61-bit prime field, prime n-bit subgroup, target and coins per run; 8-adding walk with selector
floor(8t/q); distinguished points; sequential walks to their distinguished point; run to the first
verified recovery; per-run ratios to sqrt(pi kappa n/2)). Offline (`development-evidence-mean.json`):

| n bits | dp bits | runs | mean touch / model | s.e. | mean steps to recovery / model |
|---:|---:|---:|---:|---:|---:|
| 16 | 3 | 4,000 | 0.927 | 0.008 | 0.990 |
| 20 | 4 | 3,000 | 0.948 | 0.009 | 0.971 |
| 24 | 5 | 1,500 | 0.971 | 0.013 | 0.984 |
| 28 | 6 | 400 | 0.996 | 0.026 | 0.999 |
| 32 | 7 | 100 | 0.964 | 0.053 | 0.966 |

Every row is at or below 1, and the allowance 1.05 is about two standard errors above the 28-bit row.
**Limit, disclosed:** the touch ratio rises with n from 0.93 towards 1, as finite-size effects of the
8-adding walk fade. Small groups favour early collisions: short trails of a 16-bit group revisit a
noticeable fraction of the group. The model sqrt(pi kappa n/2) is the large-n asymptotic value of the
corrected random-mapping law, which the ratios approach from below. The 32-bit row (0.964, s.e. 0.053) is
the least precise. The premise is that it stays below 1.05 at n = 2^254. The toy groups have one
coordinate, so the selector reads the whole element there. In production it reads one of twelve
coordinates of an element of a 3048-bit field, as in the base package. The declared run repeats smaller
samples of the same settings (their scatter is larger).

## 7. Experiments

The declared experiment is `lazy.py` (standard library; each part's time is reported in its output (`part_seconds`), and a guard skips the base experiment if 11 s have elapsed; parts H and I are in `walks.py`, which `lazy.py`
executes in its own namespace to keep each file under the 64 KiB limit). Part C reads `submission.json`, `proof.md` and
`machine.md` and asserts that they agree with the code and the ledger: the 84-bit predicate is the only
distinguished-point predicate stated (no 80-bit density statement anywhere), the code constants are 84 (`DP_BITS` in `lazy.py` and
`experiment.py`, the `reference.py` docstring), every parameter (215, 244.54, 8,679, 9,216, 2^139, 2^52, 2^90,
2^54, 312, K) appears where it is used, every rounded upper bound is at least
and every lower bound at most the computed value, and the manifest score equals the ledger claim
(40 checks; `all_consistent`). Part S generates the step loop, validates the
encoding of every instruction (0 violations) and executes it on 40 synthetic entries (random 128-bit F
words, random child structure with window ends, some entries forced to fail the certificate) from 43
starts (3,196 trips): every trip's A, selector, child and exit agree with the definition, every certified in-window trip
executes 215 instructions, and the BUDGET, REBASE and EXACT exits all occur. Part L generates and executes the window
product, validates the encoding of every instruction (0 violations) and executes it (20 cases: random,
edge-valued and all-limbs-large inputs, 0 failures, canonical outputs, 8,679 instructions). Part T checks the tag certificate in this package's coordinate layout on 6 window bases
produced by that routine (0 offsets outside [0, 96), all 6 certified with the exact selector), with F words computed from the row of t(Y T). Part H is the toy mean (declared input: up to 300 / 150 / 37 runs at 16 / 20 / 24 bits within a budget of at most 1 s, split
equally across the three settings and shortened if needed so that the whole run stays inside the 20 s
limit (the output reports the budget used and the runs completed), coin seed 20261113, a different stream from the offline seed 20261112; the output reports the
runs completed; the two samples are separate runs, not pooled). Part I is the integrated full-width walk of section 4. Part G is
the ledger. Part X (run last, for the declared input only) reads this run's own output and
`development-evidence-mean.json` and asserts that every experimental figure stated in `submission.json`,
`proof.md` and `machine.md` equals the measured value: step, case, trial, start, trip and window-end counts,
mismatch and failure counts, the largest certificate offset and the offline table.

Then `lazy.py` executes the base package's `experiment.py` (84-bit predicate) in-process with
`{"seed": 20261111, "trials": 16}`. That experiment covers tower checks, the row identity, the
certificate on 260 cases, the 120-step full-width walk through the real `solve` loop (23 re-bases), toy runs at the cap
(16 / 16 / 16 trials) and the recount of its reference listings (273 per trip; `fr_add` 34). Its
printed score (121.285447) is the base package's hard-capped ledger, superseded here by part G.

**Parameter note.** The base package's experiment fixed dbits = 80 in its full-width walk and in its
predicate assertion. Here both read the constant `DP_BITS = 84`
(`ref.solve(G, z, 5, DP_BITS, ...)`, `dp == ((t << DP_BITS) < p)`), so the algorithm, the ledger, the caps
and every executed check use the one predicate of section 1. Part C asserts this in the code as well
as in the documents.

References: Pollard 1978; Teske, Math. Comp. 70 (2001); Cheon, Hong and Kim, J. Cryptology 25 (2012); van
Oorschot and Wiener, J. Cryptology 12 (1999); Aranha, Karabina, Longa, Gebotys and Lopez, EUROCRYPT 2011
(eprint 2010/526); submissions 73463878 (e-man07), 8b901b3d and e2dfd4cf (terrapinelf), 67df70ff.
