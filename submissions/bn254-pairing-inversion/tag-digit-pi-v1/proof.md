# BN254 pairing inversion: gated digit-tag 8-adding rho in GT, 143 instructions per certified step, output (P0, [x]Q0)

**Claim.** For uniform z = g^x in the protected BN254 GT, the algorithm below finds x with g^x = z and returns
P = P0 and Q = [x]Q0 with e(P, Q) = z with probability at least 0.39 over uniform challenges and fresh coins, at **expected** total work

    E[W] <= 1 + (K/2^64) (143 E[T] + 9,216 (E[T]/312 + E[walks] + 1) + 10,240 E[X] + 2^110) / 20,000 < 2^119.11

reference units (W = R/20,000, R the total scalar AArch64 instruction-equivalent work). Here T is the number
of walk steps and X the number of exact evaluations of a run past the gate. This is an expectation over
the challenge and the coins of the one attack on the one target, **not a bound on every run**; every run is
capped at 2^133.16 units (section 5). Peak storage is below 2^63 bytes; preprocessing (the table, included)
below 2^72 instructions; no advice. This is a generic square-root attack, not index calculus or NFS, not a
practical break and not a security estimate. The success and expected-work bounds use disclosed
heuristics (H1, H2r, H-mean); the tag certificate and the 1/r degeneracy are proved.

**Lineage and what is new.**
- **Walk and certificate.** The walk, the certified-tag method and the table of all products of at most ell
  multipliers are those of submission 73463878-0c19-4309-a991-749ea525d51c (e-man07, after Cheon, Hong
  and Kim, J. Cryptology 25 (2012)). Its `reference.py`, `experiment.py`, `ttstep.c`, `ttstep.disasm.txt`,
  `tower.c`, `tower.disasm.txt` and `development-evidence.json` are included byte for byte, except DP_BITS = 84 in `experiment.py` (its full-width walk and predicate check) and the dbits=84 docstring in
  `reference.py` (the base used 80). They remain the reference implementation of the 128-bit-word tag.
- **Digit tag.** The window base is held as 72 digits of 43 bits, and an entry holds one 64-bit word per
  digit, so a digit costs one multiply-add and no high product. This is submission
  70974bf7-180b-47e1-a09b-c88363cfa805 (kaankolcu), adopted with credit.
- **Gate.** The run gate and the run-to-collision expected-work accounting follow submission
  8b901b3d-a617-4b1f-a442-ff370e6ec4f7 (terrapinelf).
- **Unrolled loop.** The twice-unrolled loop without a per-trip MOV and the 3-instruction certificate form
  `add x9, x4, x3, lsl #3; cmp; b.ls` follow f86cff50-80c5-43cf-a773-67417a7c29e4 (e-man07).
- **Window length.** Window length 312 follows 2da98c02-9350-4cb1-83c1-b25d1a07a4b6 (0xAJPanda).
- **Generated window product.** The window re-base uses this account's generated, encoding-validated
  lazy-reduction Fp12 product of 8,679 instructions (accepted in e3771b98 and a698079e).

New here:
(i) **a generated digit loop** (`walks.py` part S) of **143** executed instructions per certified trip
(70974bf7: 153). Its 19 resident digits mean only 53 digits are reloaded (26 LDP + 1 LDR) besides the 36
LDP of entry words, and it uses 72 multiply-adds, the 3-instruction certificate test with the constant
8E in a register, the selector and child load, and the budget count. It is unrolled twice;
(ii) **a generated digit conversion** of 171 instructions (part I executes it on every re-base), so a
re-base is 8,679 + 171 + 88 = 8,938 instructions and a window end at most 9,216;
(iii) **8 multipliers**, keeping the 8-adding walk whose collision law this account measured (H-mean,
section 6), and the gate with success charged 0.994;
(iv) explicit answers to the three objections raised against 70974bf7's pairing-inversion twin 271f78c0:
a revisitation-aware distinguished-point argument (section 6), a restricted equidistribution premise
for t, which is constant on p^2-Frobenius orbits (H2r), and an exact description of the store, which is
keyed by t (section 3);
(v) **pairing-inversion output**: this package is the pairing-inversion port of this account's GT-DLP package
05b9ff85 (`tag-digit-v6`, same walk, programs, ledger and checks); after the verified x the output is P = P0 and
Q = [x]Q0, as in 73463878 and this account's accepted pairing-inversion packages 0bb2e6af, 2ad0108b and d768ffde.
The G2 scalar multiplication, the check g^x = z and serialization are charged inside the fixed 2^110-instruction
allowance (`machine.md`, below 2^30 instructions), so the expected-work bound is the GT package's.
A step costs 143 + 9,216/312 + 10,240 x (E + 1)/2^61 <= **175.36** instructions in expectation.
`lazy.py` part G recomputes every figure below.

Notation: p, r, g are the protected modulus, group order and generator; z = g^x is the input. An element
Y of Fp12 has twelve Fp coordinates y_0..y_11 in the tower layout of `lazy.py` (coordinate 2m + b is
the b-th Fp part of the w^m coefficient over Fp2). t(Y) in [0, p) denotes its coefficient of w^0 in the
protected basis, t = Re(y_(w^0)) - 9 Im(y_(w^0)) mod p, a fixed Fp-linear functional of the coordinates.

## 1. The walk as a function of the group element

Fix eight multipliers M_s = g^(c_s) z^(d_s), s = 0..7. Define

    s(Y) = floor(8*t(Y)/p) in {0..7},      f(Y) = Y*M_{s(Y)},
    Y is distinguished  iff  t(Y)*2^84 < p.

f, s and the distinguished predicate depend on the group element only. The algorithm evaluates exactly
these functions; the table machinery is only a way to evaluate s(Y) cheaply. Hence the walk is a standard
8-adding walk with a deterministic iteration function, and two trails that reach the same element
coincide from there on.

## 2. Digit tag certificate (exact statement, proved)

**Row.** For a fixed T in Fp12, Y -> t(Y*T) is Fp-linear: t(Y*T) = sum_k y_k v_k(T) mod p. The row is
computed in `lazy.py` (`row`) and checked through part T.

**Stored data.** A coordinate is kept as its Montgomery residue Y~_k = y_k 2^256 mod p (canonical), split
into six 43-bit digits d_(k,i) = bits [43i, 43i + 43) of Y~_k, so Y~_k = sum_i d_(k,i) 2^(43i). Here
d_(k,5) < 2^39, since Y~_k < p < 2^254. For a table product T put u_(k,i) = 2^(43i) v_k(T) 2^-256 mod p and
store the 64-bit words w_(k,i) = floor(2^64 u_(k,i)/p).

**Tag.** A = sum over the 72 pairs of d_(k,i) w_(k,i) mod 2^64.

**Lemma.** With t = t(Y*T), the real number tau = 2^64 t/p satisfies (tau - A) mod 2^64 in [0, E), E = 72 2^43.

Proof. sum_i d_(k,i) u_(k,i) = Y~_k v_k 2^-256 = y_k v_k (mod p), so t/p = sum_(k,i) d_(k,i) u_(k,i)/p
(mod 1). For one pair write d = d_(k,i) < 2^43 and u/p = w/2^64 + eps, 0 <= eps < 2^-64. Then
d u/p = d w/2^64 + d eps. Modulo 1 the first term equals (d w mod 2^64)/2^64 exactly; the second lies in
[0, 2^43 2^-64). Summing 72 pairs, 2^64 t/p - A mod 2^64 lies in [0, 72 2^43). QED. Part T checks this on
window bases produced by the generated product and digit conversion, and part I checks it at every step of
a full-width walk (largest observed offset 19.929 x 2^43 in the declared run).

**Corollary (certified step).** Call A *certified* if ((A + E) mod 2^64) mod 2^61 > E. Then the low 61 bits
of A lie in [1, 2^61 - E - 1], so [A, A + E) contains no multiple of 2^61 and does not wrap. tau and A
have the same top three bits, i.e. s(Y*T) = A >> 61 exactly, and tau >= A >= 1, i.e.
t/p >= 2^-64 > 2^-84: the element is not distinguished. In the loop, x9 = 8E + 8A mod 2^64 =
8((A + E) mod 2^61), and certified iff x9 > 8E (8E = 18 2^48, one MOVZ).

## 3. Complete algorithm

Parameters: R = 8 multipliers; window length 312; distinguished points as defined in section 1
(t(Y) 2^84 < p, density 2^-84 under H2r); step cap S_cap = 2^140 steps in total; at most 2^52 walks, each
at most 2^90 steps; at most 2^127 exact evaluations; at most 2^70 store probes. Reaching any cap ends the
run as a failure.

0. **Gate.** Draw a fresh uniform 64-bit word v. Unless v < K = 7237656125499723472 (the least K with
   (K/2^64) x 0.994 >= 0.39), stop and fail (below one unit of work).
1. Decode z, check z != 0 and z^r = 1. If z = 1 output P0 and infinity (allowed exactly when z = 1).
2. Draw c_s, d_s uniformly mod r (rejection sampling of 254-bit words with at most 1,024 draws per scalar;
   exhausting them ends the run as a failure, probability below 2^-200 per scalar, charged in section 6)
   and form M_s = g^(c_s) z^(d_s), s = 0..7. Build the table: one entry e for every multiset of at most 312
   multipliers (C(320,8) = 2,496,767,727,736,920 entries). Each entry holds the 72 words of section 2 for
   T_e, the eight child pointers, alpha_e = sum c_s and beta_e = sum d_s mod r, and T_e as its 54 canonical
   Karatsuba operands (the form the generated product reads). Entries are addressed by P = entry + 512
   (`machine.md`, table entry).
3. For each walk: draw a, b uniformly mod r (same bounded sampling), set Y0 = g^a z^b (skip if 1), convert
   Y0 to its 72 digits, e = root. Repeat:
   a. A = tag(Y0, e) (one trip of the generated loop).
   b. If A is certified: s = A >> 61. If e has a child for s, set e to it; otherwise (window end) re-base:
      Y0 <- Y0 * T_e by the generated product, digits by the generated conversion, (a, b) <- (a + alpha_e,
      b + beta_e), e <- the root's child for s. This is one step Y <- Y * M_s.
   c. If A is not certified (an *exact evaluation*): re-base as in b to e = root, read t = t(Y0) exactly
      and evaluate s(Y) and the distinguished predicate by their definitions. Leave the loop at a
      distinguished point; otherwise set e to the root's child for s (one step).
4. At a distinguished point look t up in the store (2^54 slots, open addressing, **keyed by t**, the
   value of the selector coordinate, not by the whole element). If an entry (t, a', b') exists, put
   x = (a - a')/(b' - b) mod r (failure if b = b'); if g^x = z output P = P0 and Q = [x]Q0
   (affine double-and-add over Fp2), otherwise treat it as a false match. Insert (t, a, b) and start the next walk.

There is no fixed step budget: past the gate the run continues until it outputs x or a cap fires. The
global step count is a three-word counter of the driver, the per-walk length and the exact-evaluation
meter are two-word counters (`machine.md`); no counter wraps.

## 4. Correctness

By the corollary, step 3b applies exactly f of section 1 and never skips a distinguished point; step 3c
uses the definitions directly. So every walk is a trail of f from its start, stopped at its first
distinguished point, whatever representation (Y0, e) the element happens to have.

The base package's `experiment.py` (changed only to the 84-bit predicate, section 7) drives its
reference `solve` loop at full BN254 width for 120 steps (`steps_wanted` = 120; its output field `steps`
counts 121 traced positions). Against independently multiplied elements, it checks at every step that the
selector equals floor(8t/p) and that consecutive elements differ by exactly the selected multiplier, and
(every 40 steps; 4 label checks in the declared run) that Y = g^a z^b for the represented label. It
repeats the same walk with a deliberately narrow certificate (54 protected bits, 14 exact evaluations) and
obtains the identical element sequence.

The priced implementation is checked the same way (`lazy.py` part I). A full-width walk of 120 steps
(window length 8 in the test; production 312) uses 8 real multipliers g^c z^d and real table entries:
72 words from the row of t(Y T_e) in this layout, and child pointers. Every step is one trip of the
generated 143-instruction digit loop, and every window end is the generated product followed by the
generated digit conversion.
- At every step the tag equals its definition, its error is inside [0, 72 2^43), the selector equals
  floor(8 t(Y)/p), and the loop's exit and next entry match the definition.
- Each re-based Y0 equals Y0 T_e, and its digits equal the definition.
- Every 100 steps the label invariant Y = g^a z^b and the step Y -> Y M_s are checked.
- Results (declared run: 120 steps, 14 window ends, 1 label check, 0 mismatches, 0 product failures).

**Labels and recovery.** Y = g^a z^b holds at each start and is preserved by multiplying with M_s and
adding (c_s, d_s). Two distinguished points with equal t that are the same element give
g^a z^b = g^a' z^b', so x = (a - a')/(b' - b) whenever b != b' (r is prime). The check g^x = z makes a
wrong output impossible; by bilinearity e(P0, [x]Q0) = e(P0, Q0)^x = g^x = z, and both points lie in the
prescribed order-r subgroups (P0 is the profile generator, Q a multiple of Q0). Because the store is keyed by t, a hit can also be a different element with the
same t, for instance a p^2-Frobenius conjugate (t is constant on those orbits, section 6). Such a false
match fails the check, costs one exponentiation, is metered in the fixed allowance, and the walk continues
with the next start. It never produces an output.

**Degeneracy** is exactly 1/r, without heuristics. Fix x and condition on the logarithms m_s = c_s + x d_s
of the multipliers and s_i = a_i + x b_i of all walk starts. Every group element, tag, exact evaluation,
distinguished point and store hit is a function of these. Conditionally the b_i and d_s remain
independent and uniform. For a hit between walk i and an earlier walk j, b - b' = b_i - b_j + (an integer
combination of the d_s), which is uniform on Z/r because b_i is; it vanishes with probability exactly 1/r.

## 5. Work

Per-item prices (`machine.md`; the loop, the product and the digit conversion are generated,
encoding-validated and executed by `lazy.py`/`walks.py`; `fr_add` is recounted from `tower.disasm.txt`):

| item | instructions |
|---|---:|
| one certified step (one trip of the generated digit loop: 71 MADD + 1 MUL, 62 LDP + 1 LDR, 3-instruction certificate test, selector, child load, exit test, counter) | 143 |
| one window end: re-base 8,938 (generated product 8,679 + digit conversion 171 + two `fr_add` 68 + glue 20), loop entry 20 + exit 13, driver with the three-word counter <= 160 | 9,131 <= **9,216** |
| one exact evaluation: re-base 8,938, exact t and selector, loop exit and re-entry, driver, meter (itemised in machine.md) | <= **10,240** |
| everything else, by meters and caps: table build, <= 2^52 walk starts, store, recoveries | <= 2^110 |

A window end occurs after 312 steps since the previous re-base or at a walk start (re-bases inside exact
evaluations are priced in those), so the number of window ends is at most T/312 + walks + 1.

**Every run** (independent of all heuristics): T <= 2^140, walks <= 2^52 and exact evaluations <= 2^127, so
R <= 143 x 2^140 + 9,216 (2^140/312 + 2^52 + 1) + 10,240 x 2^127 + 2^110 and W <= 2^133.16 units.

**Expected work.** Let T_touch be the number of steps (all walks together) until some step reaches an
element visited by an earlier trail. After it, the touching walk reaches its distinguished point within its
remaining length: at most 2^90 steps, 2^84 in expectation under H2r. Hence E[T] <= E[T_touch] + 2^85. The
number of walks is at most one more than the number of distinguished points reached, which is T/2^84 in
expectation under H2r (charged twice: E[walks] <= 2 E[T]/2^84 + 1). A step's tag is uncertified with
probability (E + 1)/2^61 = 2^-11.83 under H2r, so E[X] <= E[T] (E + 1)/2^61. Every certified step costs the
fixed 143, so

    E[R_run] <= 143 E[T] + 9,216 (E[T]/312 + E[walks] + 1) + 10,240 E[X] + 2^110.

**H-mean.** E[T_touch] <= 1.05 sqrt(pi kappa r / 2), kappa = 8/7 (the 8-adding correction of H1). Then
E[T] <= 2^127.2909, E[X] <= 2^115.47, E[R_run]/20,000 <= 2^120.4572 units, and with the gate,
W <= 1 + 1[v < K] W_run with v independent of the run, so

    E[W] <= 1 + (K/2^64) x 2^120.4572 <= 2^119.1075.

Submitted bound: **time_log2 = 119.11** (rational 11911/100); the computed value is 119.10742308..., so the
claim is at least 0.0025 bits above it. Per step this is 143 + 29.54 + 2.82 <= 175.36 instructions in
expectation. The same ledger gives 2^119.5858 with this account's 215-instruction 128-bit-word loop.

## 6. Success probability: heuristics, support and limits

Past the gate the run fails only if a cap fires or a recovery degenerates. Under H1, H2r and H-mean
(`lazy.py` part G, `losses`):
* the step cap: P(T > 2^140) <= E[T]/2^140 <= 0.000150 (Markov; exact ledger value 0.00014933...);
* the walk cap: P(walks > 2^52) <= E[walks]/2^52 <= 0.004779 (Markov; exact ledger value 0.0047786...);
* the exact-evaluation meter: P(X > 2^127) <= E[X]/2^127 <= 0.000336 (Markov; exact ledger value
  0.00033599...);
* the store-probe meter (2^54 slots, at most 2^52 keys, linear probing, O(1) expected probes): below 2^-20;
* a walk not distinguishing within 2^90 = 64 x 2^84 steps, or closing a cycle on itself first, in expectation over E[walks] <= 2^44.3 walks: below 2^-41 (see the revisitation argument below);
* degeneracy b = b' (exactly 1/r per hit, section 4), scalar sampling failure, a start equal to 1: below
  2^-200.

The sum is below 0.0053 (ledger: 0.0052648...), so success past the gate is at least 0.9947; the claim uses
0.994. With the gate, success is at least (K/2^64) x 0.994 >= 0.39 and K is the least such integer (part G:
`success_ok`, `K_is_least`). Once a touch has happened, no other event prevents success: the later trail
follows the earlier one to its stored distinguished point.

**H1** (collision law of the 8-adding walk, kappa = 8/7) is the base package's. The selector is the
classical leading-bits selector of finite-field rho, and the walk is Teske's r-adding walk (Math. Comp. 70
(2001)).

**H2r (restricted equidistribution), stated so as to answer the review of 271f78c0.** t(Y) is the w^0
coefficient, fixed by p^2-Frobenius: t(Y^(p^2)) = t(Y), and p^2 has order 6 on GT minus 1. So t takes at
most 1 + (r - 1)/6 values on GT and is *not* uniform on Fp. The premise is restricted to what the algorithm
uses: over the run's *first visits to p^2-Frobenius orbits* (a step is fresh only if no element of its
orbit {Y^(p^(2k))} was visited before in the run, since orbit-mates share t), the events "top three bits of tau equal s", "certificate
fails" (probability (E + 1)/2^61) and "t 2^84 < p" (probability 2^-84) occur at their measure on Fp. That
is, the image of t on GT (one value per Frobenius orbit) is balanced for these intervals of lengths 2^61
and 2^-84 p, which are far longer than the spacing of the image. No uniformity conditional on revisits or on orbit-mates
is assumed.

**Orbit-mate meetings (charged).** A step that lands on Y' = Y^(p^(2k)), k != 0 mod 6, for an earlier
element Y of the run is not fresh: t(Y') = t(Y), so its selector, certificate status and distinguished status
repeat Y's. Under H1 a step lands on an orbit-mate of one of the earlier visited elements
with probability about 5 x (number visited)/r, so up to step T the expected number of meetings is at most
5 E[T^2]/(2r). With the Rayleigh law of H1, E[T^2] = (4/pi) E[T]^2 <= (4/pi) 1.05^2 (pi kappa r/2) =
2 x 1.05^2 kappa r, so the expected number per run is at most 5 x 1.05^2 x kappa = 6.30 < 2^3. Each one makes a single step non-fresh: the next element Y' M_s is not an orbit-mate of
Y M_s, because the multipliers are not Frobenius-fixed. Its effects are charged as follows. (i) The walk-length
tail below counts first visits of orbits; at most m non-fresh steps lengthen a walk by m, changing the bound
exp(-L 2^-84) by a factor exp(m 2^-84), which is 1 to within 2^-23 for every m < 2^60 (and m exceeds
2^60 with probability below 2^-57 by Markov's inequality on its mean 6.30). (ii) A non-fresh step may exit
(uncertified) even if Y certified; this adds at most 2^3 exits to E[X] = 2^115.46 and at most
2^3 x 10,240 instructions to the run, inside the fixed allowance. (iii) If Y was distinguished, Y' is too; it
is then a store hit with equal t from a different element, a false match caught by g^x = z (section 3), which
costs one exponentiation and ends that walk. None of these changes a bound by more than the rounding already
applied.

**Revisitation-aware tail argument.** A walk's trajectory, until its first distinguished point, consists
of first visits, except possibly for its first revisit. That revisit is either (a) to an element of an
earlier trail, a *touch*, after which it follows that trail to its stored distinguished point, which is
success, or (b) to an element of its own trail, a *self-cycle*. A non-distinguished element revisited is
certainly not distinguished. The argument therefore never treats a revisit as a fresh event; it only
charges (b).
- Under H2r, the distinguished indicator of a first visit is a fresh event of probability 2^-84. A walk
  that has made L first visits without a distinguished point therefore has probability at most
  (1 - 2^-84)^L <= exp(-L 2^-84) of doing so; for L = 2^90 this is exp(-64) = 2^-92.3.
- Under H1, a self-cycle within a walk of length at most 2^90 has probability at most about
  E[len^2]/(2r) = 2^169/(2r) <= 2^-85.6, with E[len^2] = 2 2^168 for the geometric length.
- A walk is capped by length only if one of these happens, so the charged term is
  E[walks] (2^-92.3 + 2^-85.6) <= 2^44.3 x 2^-85.5 < 2^-41.

**H-mean** (mean touch time <= 1.05 of the model) is a statement about the mean only. Support: the walk
itself at toy scale (`walks.py` part H: fresh 61-bit prime field, prime n-bit subgroup, target and coins per
run; 8-adding walk with selector floor(8t/q); distinguished points; sequential walks to their
distinguished point; run to the first verified recovery; per-run ratios to sqrt(pi kappa n/2)). Offline
(`development-evidence-mean.json`):

| n bits | dp bits | runs | mean touch / model | s.e. | mean steps to recovery / model |
|---:|---:|---:|---:|---:|---:|
| 16 | 3 | 4,000 | 0.927 | 0.008 | 0.990 |
| 20 | 4 | 3,000 | 0.948 | 0.009 | 0.971 |
| 24 | 5 | 1,500 | 0.971 | 0.013 | 0.984 |
| 28 | 6 | 400 | 0.996 | 0.026 | 0.999 |
| 32 | 7 | 100 | 0.964 | 0.053 | 0.966 |

A second, independent offline run (`development-evidence-mean2.json`: ten processes, seeds
20261201..20261210, disjoint from the first file's 20261112 and the declared run's 20261113) adds larger
sizes:

| n bits | dp bits | runs | mean touch / model | s.e. | mean steps to recovery / model |
|---:|---:|---:|---:|---:|---:|
| 24 | 5 | 2,000 | 0.992 | 0.011 | 1.002 |
| 28 | 6 | 1,800 | 1.002 | 0.012 | 1.008 |
| 32 | 7 | 600 | 1.005 | 0.022 | 1.008 |
| 36 | 8 | 100 | 0.982 | 0.055 | 0.989 |

Pooled at n >= 28 bits over both files: 3,000 runs, mean 1.0002, s.e. 0.0096; mean + 3 s.e. = 1.0291.

**Limit, disclosed.** The touch ratio rises with n from 0.93 at 16 bits to about 1.00 at 28-36 bits. Small
groups favour early collisions; the model sqrt(pi kappa n/2) is the large-n value of the corrected
random-mapping law. The charged allowance 1.05 leaves headroom above the large-n data. Independent
larger runs of the same walk by another lane of this account (not shipped here) read slightly above 1 at
32 bits (about 1.01). The premise is that the ratio stays below 1.05 at n = 2^254, an extrapolation from
groups of at most 36 bits. The toy groups have one coordinate; in production the selector reads t, one
Fp-linear functional of an element of a 3048-bit field, as in the base package.

## 7. Experiments

The declared experiment is `lazy.py` (standard library). Each part's time is reported in its output
(`part_seconds`), and a guard skips the base experiment if 11 s have elapsed. Parts H, I and S are in
`walks.py`, which `lazy.py` executes in its own namespace to keep each file under the 64 KiB limit.
- **Part C** reads `submission.json`, `proof.md` and `machine.md` and asserts that they agree with the code
  and the ledger:
  - the 84-bit predicate is the only distinguished-point predicate stated (no 80-bit density statement
    anywhere);
  - the code constants are 84 (`DP_BITS` in `lazy.py` and `experiment.py`, the `reference.py` docstring);
  - every parameter (143, 175.36, 8,679, 9,216, 10,240, 2^140, 2^52, 2^90, 2^54, 2^127, 312, K) appears
    where it is used;
  - every rounded upper bound is at least, and every lower bound at most, the computed value;
  - the manifest score equals the ledger claim.

  (47 checks; `all_consistent`)
- **Part S** generates the digit loop, validates the encoding of every instruction (0 violations) and
  executes it on 40 synthetic entries from 43
  starts (3,204 trips). The entries have random 64-bit words, 43-bit digits, a random child structure with
  window ends, and some entries forced to fail the certificate. Every trip's A, selector, child and exit
  agree with the definition, every certified in-window trip executes 143 instructions (both trips), and
  the BUDGET, REBASE and EXACT exits all occur, from both trips.
- **Part L** generates the window product, validates the encoding of every instruction (0 violations) and
  executes it (20 cases: random, edge-valued and all-limbs-large inputs, 0 failures, canonical outputs,
  8,679 instructions).
- **Part T** checks the digit tag lemma on 6 window bases produced by the generated product and digit
  conversion: 0 offsets outside [0, 72 2^43), digits equal to the definition, and the exact selector on
  every certified tag.
- **Part H** is the toy mean. Declared input: up to 300 / 150 / 37 runs at 16 / 20 / 24 bits within a
  budget of at most 1 s, split equally across the three settings and shortened if needed so that the
  whole run stays inside the 20 s limit. The output reports the budget used and the runs completed. The
  coin seed is 20261113, a different stream from the offline seeds, and the samples are separate runs,
  not pooled.
- **Part I** is the integrated full-width walk of section 4.
- **Part G** is the ledger.
- **Part X** runs last, for the declared input only. It reads this run's own output and the evidence files,
  and asserts that every experimental figure stated in `submission.json`, `proof.md` and `machine.md`
  equals the measured value: step, case, trial, start, trip and window-end counts, mismatch and failure
  counts, the largest offsets and the offline tables.

Then `lazy.py` executes the base package's `experiment.py` (84-bit predicate) in-process with
`{"seed": 20261111, "trials": 16}`. That experiment covers tower checks, the row identity, its 128-bit-word certificate on 260 full-width cases (largest observed offset 41 in the declared run), the 120-step full-width walk
through its reference `solve` loop (23 re-bases), toy runs at the cap (16 / 16 / 16 trials) and the
recount of its reference listings (273 per trip; `fr_add` 34). Its printed score (121.285447) is the base
package's hard-capped ledger, superseded here by part G.

References: Pollard 1978; Teske, Math. Comp. 70 (2001); Cheon, Hong and Kim, J. Cryptology 25 (2012); van
Oorschot and Wiener, J. Cryptology 12 (1999); Aranha, Karabina, Longa, Gebotys and Lopez, EUROCRYPT 2011
(eprint 2010/526); submissions 73463878 and f86cff50 (e-man07), 70974bf7 and 271f78c0 (kaankolcu),
8b901b3d and e2dfd4cf (terrapinelf), 2da98c02 (0xAJPanda), 67df70ff, e3771b98, a698079e.
