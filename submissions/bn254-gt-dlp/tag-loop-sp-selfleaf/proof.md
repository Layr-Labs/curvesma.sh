# Gated tag-tracing 8-adding rho in BN254 GT: 213-instruction step loop (SP-relative limbs, self-pointer leaves) and lazy window product

**Claim.** For uniform z = g^x in the protected BN254 GT, the algorithm below returns x with g^x = z
with probability at least 0.39 over uniform challenges and fresh coins, at **expected** total work

    E[W] <= 1 + (K/2^64) (213 E[T] + 9,216 (E[T]/256 + E[walks] + 1) + 2^110) / 20,000 < 2^119.65

reference units (W = R/20,000, R the total scalar AArch64 instruction-equivalent work), where T is the number
of walk steps of a run past the gate. This is an expectation over the challenge and the coins of the one
attack on the one target, **not a bound on every run**; every run is capped at 2^126.68 units (section 5).
Peak storage is below 2^61.7320 bytes; preprocessing (the table, included) below 2^69 instructions; no advice.
This is a generic square-root attack, not index calculus or NFS, not a practical break and not a security
estimate. The success and expected-work bounds use disclosed heuristics (H1, H2, H-mean); the tag
certificate and the 1/r degeneracy are proved.

**Parent and changes.** This package is the gated tag-tracing package e3771b98-7faa-4277-a559-c5feb10a85f4
(source b3aa659ab7393ebaaee785a4e583e2be8396aa03, review 119.68, not promoted). Parent code and evidence are retained; the proof and manifest describe the new loop, store and corrected mean. The scripts `lazy.py` and `walks.py` have unchanged syntax trees and comments, with compact indentation. Their recorded window-product (part L), certificate-in-layout (part T), toy-mean (part H) and integrated-walk (part I) checks concern only unchanged components, not the new store or corrected ledger. The parent's own lineage is the tag-tracing walk of submission 73463878 and the gate accounting of
submissions 8b901b3d and e2dfd4cf. Its description of that lineage is kept in section 7 and in
`development-evidence*.json`.

**Step-loop implementation change.** The parent's generated loop executed 217
instructions per certified trip. The new loop (`ttloop2.py` part S2) executes **213**. It computes the same
A = sum_j (UMULH(y_j, Flo_j) + y_j Fhi_j) mod 2^64, the same certificate and the same selector on the same
entry layout. There are three sources of the saving:
1. **SP-relative limbs.** The driver keeps the window base Y0 in its own frame, at the caller's SP. The
   kernel's frame is 112 bytes, so limb j sits at [SP, #112 + 8j], at most #488, which is encodable. The
   kernel reads non-resident limbs SP-relative, so x1 (the parent's y pointer) is free.
2. **No scratch registers.** The certificate and the selector use x4, which is free after the last
   accumulation, so x9 and x10 are free. Together with x1 this gives 22 resident limbs instead of 19, and
   26 reloaded limbs instead of 29: 13 LDP per trip instead of 14 LDP + 1 LDR. That saves 2 instructions.
3. **Counted window, no null test.** The child is loaded straight into x0 (`ldr x0, [x0, x4, lsl #3]`), and
   the parent's `cbz x10` and `mov x0, x10` are removed. That saves 2 instructions. To make this possible:
   - each window leaf (multiset size 256) stores its own pointer in all eight child slots instead of null;
   - the driver calls with budget = min(steps left, walk left, 257 - d0), where d0 is the depth at the call:
     d0 = 0 at a walk start (root) and d0 = 1 after a re-base or an exact evaluation;
   - each certified trip moves from an entry of size k to one of size k + 1, so a call that no counter
     limits makes exactly 257 - d0 trips, and the last trip is the leaf's own trip;
   - the driver recognises that case (budget equal to 257 - d0 and BUDGET status, a one-word compare and a
     branch) and re-bases with s = A >> 61 from the cursor, exactly as the parent's REBASE exit did;
   - a counter-limited call ends the walk or the run, as in the parent;
   - the leaf's own trip counts as one step, as the parent's REBASE trip did.

Per step: 213 + 9,216/256 = **249** instructions instead of 253. The window end is re-stated as follows:
product 8,679 + labels 88 + loop entry 19 (11 LDP for 22 resident limbs) + one exit 12 + driver 168 (the
parent's 160 plus the leaf test) = 8,966, still within the charged **9,216**. The table changes only in the
leaf child slots (self instead of 0), at the same build cost. Every other price, cap, constant, heuristic,
the mean allowance, the gate constant and the rounding rule (ceiling at two decimals) are the parent's.
Under the same ledger formula the computed value falls from 119.67157621 to **119.64858457** (0.02299 bits).
The claim is **119.65** against the parent's 119.68, which included 0.0084 bits of rounding slack.

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

Parameters: R = 8 multipliers; window length 256; distinguished points as defined in section 1
(t(Y) 2^84 < p, density 2^-84 under H2; this package uses 84 bits everywhere, the base package used 80);
step cap S_cap = 2^133
steps in total; at most 2^52 walks, each at most 2^90 steps; at most 2^90 exact evaluations. Reaching any
cap ends the run as a failure. The distinguished-point store is deterministic and has no probe cap (step 4).

0. **Gate.** Draw a fresh uniform 64-bit word v. Unless v < K = 7416732153347139310 (the least K with
   (K/2^64) x 0.97 >= 0.39), stop and fail (below one unit of work).
1. Decode z, check z != 0 and z^r = 1. If z = 1 output 0.
2. Draw c_s, d_s uniformly mod r (rejection sampling of 254-bit words with at most 1,024 draws per scalar;
   exhausting them ends the run as a failure, probability below 2^-200 per scalar, charged in section 6)
   and form M_s = g^(c_s) z^(d_s), s = 0..7. Build the table: one entry e for
   every multiset of at most 256 multipliers (C(264,8) = 525,783,425,977,953 entries), holding the 96 words
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
4. At a distinguished point, take the full element Y. It is the re-based Y0 produced by the exact
   evaluation, 12 canonical Montgomery residues, 48 limbs, 3,072 bits. Look it up in the store, a
   **deterministic binary Patricia (crit-bit) trie keyed by the full canonical element** (machine.md
   "Distinguished-point store").
   - If a record (Y, a', b') for the same element exists, put x = (a - a')/(b' - b) mod r (failure if
     b = b'), check g^x = z and output x. Key equality is element equality, because Montgomery form is a
     bijection on canonical residues. A stored hit is therefore always a true collision, and there are no
     false matches.
   - Otherwise, insert (Y, a, b) and start the next walk.

   Store properties:
   - Branch bit positions strictly decrease along every path, so any path has at most 3,072 internal
     nodes.
   - N distinct elements use exactly N leaves and N - 1 internal nodes.
   - Worst case per operation: call 64 + two traversals 2 x 3,072 x 16 + 48-limb compare 192 + first
     differing bit 200 + link 16 + leaf write and element materialisation 96 = 98,872 instructions,
     charged 2^17.
   - There are at most 2^52 operations (one per walk), so at most 2^69 instructions.
   - The store is deterministic. It adds no failure term and no cap.

There is no fixed step budget S: past the gate the run continues until it outputs x or a cap fires. The
global step count is a three-word (192-bit) counter of the driver, the per-walk length a two-word counter
(`machine.md`); 2^133 < 2^192, so neither wraps.

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
(window length 8 in the test; production 256) uses 8 real multipliers g^c z^d and real table entries
(F words from the row of t(Y T_e) in this layout, child pointers). Every step is one trip of the
generated 213-instruction loop (`ttloop2.py`; part I of the parent ran its 217-instruction loop, whose tag, certificate and selector are identical) and every window end is the generated product. At every step the
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
| one certified step (one trip of the generated loop: 48 UMULH, 47 MADD + 1 MUL, 48 ADD, 61 LDP, child LDR into x0, certificate test, selector, counter) | 213 |
| one window end: generated product 8,679 + two `fr_add` 68 + glue 20 (re-base 8,767), loop entry 19 + exit 12, driver with the three-word counter and the leaf test <= 168 | 8,966 <= **9,216** |
| everything else, by meters and caps: table build, <= 2^52 walk starts, <= 2^90 exact evaluations, store, recoveries | <= 2^110 |

A window end occurs after 256 steps since the previous re-base, at a walk start, or within an exact
evaluation (priced in that allowance), so the number of window ends is at most T/256 + walks + 1.

**Every run** (independent of all heuristics): T <= 2^133 and walks <= 2^52, so
R <= 213 x 2^133 + 9,216 (2^133/256 + 2^52 + 1) + 2^110 and W <= 2^126.68 units.

**Expected work.** Let T_touch be the number of steps (all walks together) until some step reaches an
element visited by an earlier trail. After it, the touching walk reaches its distinguished point within its
remaining length: at most 2^90 steps, 2^84 in expectation under H2. Exceptional paths are charged the
whole step cap 2^133:
- a walk closing a cycle on itself before its distinguished point, probability at most
  2^52 x 2^(2 x 84 + 1)/r;
- a walk with no distinguished point within 2^90 steps, probability at most 2^52 e^-64.

T_touch is a nonnegative integer, so E[T_touch] = sum over s >= 0 of P[T_touch > s] <= 1 + (the
integral), which is the discrete +1. Hence

    E[T] <= M = 1.05 sqrt(pi kappa r/2) + 1 + 2^85 + (2^52 x 2^169/r + 2^52 e^-64) x 2^133

(this corrected mean bound is shared with this author's pairing package).
The number of walks is at most one more than the number of distinguished points reached, which is T/2^84
in expectation under H2 (charged twice: E[walks] <= 2 E[T]/2^84 + 1). The cost of every step is the fixed
213 (the loop does not depend on the branch), so no Wald argument on step prices is needed:

    E[R_run] <= 213 E[T] + 9,216 (E[T]/256 + E[walks] + 1) + 2^110.

**H-mean.** E[T_touch] <= 1.05 sqrt(pi kappa r / 2), kappa = 8/7 (the 8-adding correction of H1). Then
E[T] <= M <= 2^127.2909, E[R_run]/20,000 <= 2^120.9631 units, and with the gate, W <= 1 + 1[v < K] W_run with v
independent of the run, so

    E[W] <= 1 + (K/2^64) x 2^120.9631 <= 2^119.6486.

Submitted bound: **time_log2 = 119.65** (rational 11965/100); the computed value is 119.64858457..., so the
claim is at least 0.0014 bits above it. All displayed upper bounds are rounded upward (`ttloop2.py` part G2,
`displayed_upper_bounds_rounded_up`). 119.65 is a claim until it is evaluated. With the parent's 217-instruction loop the same ledger gives
119.67157621... With the base loop (273) the same ledger gives 119.96004566...; for
comparison, the base package's hard-capped run at S = sqrt(1.409 r) steps and 369 instructions per step is
2^121.285447 on every run.

## 6. Success probability: heuristics, support and limits

Past the gate the run fails only if a cap fires or a recovery degenerates. Under H1, H2 and H-mean
(`ttloop2.py` part G2, `losses`):
* the step cap: P(T > 2^133) <= E[T]/2^133 <= 0.019115 (Markov; exact ledger value 0.0191144...);
* the walk cap: P(walks > 2^52) <= E[walks]/2^52 <= 0.004779 (Markov; exact ledger value 0.0047786...);
* the exact-evaluation meter: a step is uncertified with probability 97/2^61 under H2, so the expected
  number is below 2^72.9, so P(> 2^90) < 2^-17.1;
* the store: deterministic, keyed by the full element; no probe meter, no false match, no failure term;
* a walk not distinguishing within 2^90 = 64 x 2^84 steps (exp(-64) per walk), or closing a cycle on itself
  first (about E[len^2]/(2r) = 2^169/r per walk), bounded over at most 2^52 walks: below 1.6e-10
  (also charged in M);
* degeneracy b = b' (exactly 1/r per hit, section 4), scalar sampling failure, a start equal to 1: below
  2^-200.
The sum is below 0.0240 (ledger: 0.02390010...), so success past the gate is at least 0.9760; the claim uses
0.97. With the
gate, success is at least (K/2^64) x 0.97 >= 0.39 and K is the least such integer (part G: `success_ok`,
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

**Declared initial run: `ttloop2.py`** (new, standard library, under 1 s). It has two parts.

Part S2 generates the 213-instruction loop and validates the encoding of every instruction:
- LDP/STP offsets are multiples of 8 in [-512, 504], including the SP-relative limb loads at #112..#488;
- the 61-ones AND mask is a valid bitmask immediate;
- x18 is not used.

It then executes the loop on 64 synthetic entries: random 128-bit F words, random children, window leaves
whose children point to themselves, and some entries forced to fail the certificate. There are 64 calls
with budgets of 256, random values and 1 or 2. At every trip, A, the certificate, the selector, the next
entry, the exit status and the stored cursor (entry, A, remaining budget) are compared with the definition.
The declared run (seed 20261112) gives:
- 6,672 trips checked, 0 mismatches, 0 encoding violations;
- every certified trip executes exactly 213 instructions;
- 45 BUDGET and 19 EXACT exits, 15 calls ending on a leaf self-pointer;
- Y0 in the caller's frame and SP are unchanged after every call.

**Part R (full driver)** runs the generated kernel under the driver of machine.md on a complete
multiset table with window W = 4. The production window is 256; the loop body does not depend on W. The
table has 495 entries, self-pointer leaves, and some entries with all-zero F words, so that their tag is
uncertified for every Y0. The driver uses budget = min(steps left, W + 1 - d0), re-bases after a
counter-unlimited BUDGET return, and takes exact evaluations. Its sequence of mathematical steps (child
step, leaf re-base, exact step, with selectors) is compared one by one with the parent's semantics (null
leaves, REBASE exit, EXACT exit), where every trip is one step. Declared run: 6,000 steps, an identical
sequence, 1,335 leaf re-bases and 321 exact steps.

**Part D (store)** runs the crit-bit trie with 3,072-bit keys built from 12 coordinates below p, on 980
operations against a dictionary reference:
- 601 random elements;
- 139 single-bit neighbours of one element, at bits 0, 1, 62, 63, 64, 65, 127, 128, 191, 192, 252 and 253
  of every coordinate (every limb boundary);
- 36 nested shared-prefix elements;
- 4 edge elements;
- 200 exact duplicates.
It finds 0 mismatches and exact node counts (778 leaves, 777 internal). The maximum depth seen is 175; the
bound of 3,072 holds because branch positions strictly decrease. The maximum counted cost is 5,880, below
the analytical 98,872, below the charged 2^17. The store is 2^52 x 512 + (2^52 - 1) x 32 bytes =
2^61.0875 bytes. Peak storage, including the table (2^60.2590) and 2^30 bytes of other data, is below
2^61.7320 bytes.

**Part G2** computes the ledger with the corrected mean M: 119.64858459 (with the parent's 217 under the
same formula, 119.67157623), window end 8,966 <= 9,216, success past the gate 0.97610, K least, every-run
bound below 2^126.6723.

**Historical, not declared: `lazy.py` and `walks.py`.** These are the parent's experiment, kept byte for
byte as attributed evidence. They are **removed from the declared experiments**, because their parts C and X
assert the parent's document figures (217, 253, 119.68, the open-addressing store) and would fail against
this package. Their recorded outputs are in `development-evidence.json` and
`development-evidence-mean.json`. They support the unchanged window product (part L, 8,679 instructions),
the certificate in this layout (part T), the toy mean (part H) and the integrated full-width walk with the
parent's loop (part I). Nothing in this package's current validation depends on running them.

The parent's experiment description follows for reference.



The declared experiment is `lazy.py` (standard library; about 10 s in the local Docker runner, every part's
time is reported, and a guard skips the base experiment if 11 s have elapsed; parts H and I are in `walks.py`, which `lazy.py`
executes in its own namespace to keep each file under the 64 KiB limit). Part C reads `submission.json`, `proof.md` and
`machine.md` and asserts that they agree with the code and the ledger: the 84-bit predicate is the only
distinguished-point predicate stated (no 80-bit density statement anywhere), the code constants are 84 (`DP_BITS` in `lazy.py` and
`experiment.py`, the `reference.py` docstring), every parameter (217, 253, 8,679, 9,216, 2^133, 2^52, 2^90,
2^54, 256, K) appears where it is used, every rounded upper bound is at least
and every lower bound at most the computed value, and the manifest score equals the ledger claim
(39 checks; `all_consistent`). Part S generates the step loop, validates the
encoding of every instruction (0 violations) and executes it on 40 synthetic entries (random 128-bit F
words, random child structure with window ends, some entries forced to fail the certificate) from 43
starts (3,196 trips): every trip's A, selector, child and exit agree with the definition, every certified in-window trip
executes 217 instructions, and the BUDGET, REBASE and EXACT exits all occur. Part L generates and executes the window
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
(eprint 2010/526); submissions 73463878 , 8b901b3d and e2dfd4cf , 67df70ff.
