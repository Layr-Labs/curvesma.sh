# Fixed-G1 pairing inversion via O(1)-memory Floyd rho in GT

Track: `bls12-381-fixed-g1-invert-g2`.

* Claimed score: **time_log2 = 129.40** (647/5) target-group-operation
  equivalents (dense Fp12 multiplication at C = 20,000 scalar AArch64 instruction
  equivalents).
* Success probability: at least 0.39 under the textbook random-mapping heuristic,
  supported by experiments.

This is a **small-memory, deterministic-budget** Pollard rho with **Floyd
tortoise–hare cycle detection**. Its entire working set — two group elements with
their exponent pairs, and a modest adding table — is a few gigabytes, well within
the 2^64-byte native address space, so there is **no large/segmented/distributed
memory** and its associated access cost (the obstacle for a baby-step
giant-step attack with ~2^138-byte storage). The only modelling assumption is the
random-mapping collision-time heuristic, which the review accepts as qualifying.

## 1. Reduction

Input: nonidentity P in G1 (order r), z in GT. Compute g = e(P, Q0) (pinned
`pairing(Q0, P)`); nondegeneracy and prime order make g a generator of GT. Find a
with z = g^a (section 2) and output Q = [a]Q0; e(P,Q) = g^a = z. If z = 1, output
infinity. P is never replaced.

## 2. Floyd rho for z = g^a

GT is cyclic of prime order r. Use an R-adding walk with R = 2^20 precomputed
entries M_j = g^{alpha_j} z^{beta_j} (random alpha_j, beta_j in [0,r)). The walk
is the **deterministic function** f(x) = x·M_{H(x)}; a walk point x = g^A z^B
carries its exponents (A,B), updated A += alpha_j, B += beta_j mod r at each step.

**Partition function H (concrete).** Represent x by its canonical little-endian
byte encoding: the 12 Fp coefficients c_0,...,c_11 of the flat w-basis (section 2
of the target profile), each the 48-byte canonical integer in [0,p) exported from
the internal representation (section 4). Define

    H(x) = (the low 20 bits of c_0's canonical encoding) ∈ [0, 2^20).

Reading them is one load of the low 3 bytes of c_0 and a 20-bit mask (≈ 4
instructions). H is a fixed deterministic function of the group element, so f is a
function; modelling it as a near-uniform R-way hash is the random-mapping
assumption (R = 2^20 ≥ 16, well in Teske's regime). The canonical encoding is
computed once per walk point as part of the multiplication's output (section 4),
so H adds only the mask.

**Floyd detection.** Maintain a tortoise (x_t, A_t, B_t) advanced one step per
iteration and a hare (x_h, A_h, B_h) advanced two steps per iteration, both
starting at a common random point. After i iterations the tortoise is at x_i and
the hare at x_{2i}. Because f is a function on the finite set GT, the trajectory
x_0, x_1, ... is eventually periodic with tail μ and cycle λ; the tortoise and
hare coincide at the first iteration i that is a multiple of λ with i ≥ μ, and
that i satisfies **i ≤ μ + λ = ρ** (the rho length). At the meeting
x_i = x_{2i}, so g^{A_t} z^{B_t} = g^{A_h} z^{B_h}, giving
a = (A_h − A_t)/(B_t − B_h) mod r provided B_t ≠ B_h. B_t = B_h occurs with
probability ~1/r (a negligible exceptional event counted as a failure below).

**Hard budget.** Run at most **T = ceil(1.25·k0)** tortoise iterations, k0 =
sqrt(2 r ln(1/0.61)); each iteration is exactly 3 group multiplications (1
tortoise + 2 hare) and one element comparison, then stop. The total is a hard
deterministic ≤ 3T + table-build multiplications.

## 3. Success probability (within the hard budget)

Model the R-adding walk as a random mapping on the r elements. The rho length ρ
(first value-repeat) satisfies P(ρ > t) ≈ exp(−t^2/(2r)). Floyd detects by
tortoise iteration i ≤ ρ, so the bounded success probability is at least

    P(ρ ≤ T) = 1 − exp(−(1.25·k0)^2/(2r)) = 1 − 0.61^{1.5625} = 0.538 ≥ 0.39,

with margin ~0.148. (If ρ > T the run fails; the ~1/r unsolvable meeting is
likewise counted as a failure. Both are exceptional events already inside the
1 − 0.538 shortfall.) The toy experiment measures the actual bounded success.

### Experimental support (toy scale, Python 3.12, one core)

`experiment.py` runs this exact O(1)-memory Floyd rho (tortoise one step, hare
two, meet → solve) with the hard budget T = ceil(1.25·k0), in the prime-order-ell
subgroup of F_q*. It reports the fraction of trials finding a usable collision
within T.

| bits | seeds | success within T |
| ---: | --- | ---: |
| 24 | 1–5 | 0.645–0.785 |
| 28 | 5 | 0.54 (at budget factor 1.0; ≥ 0.6 at 1.25) |

Success is comfortably ≥ 0.39 and exceeds the 0.538 lower bound (Floyd detects at
i ≤ ρ, often earlier). Toy-scale evidence; standard support for R-adding walks as
random mappings: Teske (2001).

## 4. Cost (hard deterministic bound)

* Each tortoise iteration = 3 dense Fp12 multiplications (1 unit each). Other
  per-iteration work: three exponent updates (each two modular additions of
  255-bit values ≈ 160 instr; ≈ 480), three H masks (≈ 12), three
  native-addressable table loads of the M_j element used by the multiplication
  (≈ 60), one full-element tortoise/hare comparison (72 words ≈ 220), and control
  (≈ 120): ≤ 892 instructions over 3 multiplications ≈ 297 instr/mult ≤ 0.015
  units/mult. We **declare a per-multiplication bound of 1.02 units**, which
  covers up to 400 extra instructions per multiplication — comfortably above the
  ~297 estimate.
* Multiplications ≤ 3·T + (one-time build). 3·T = 3·1.25·k0 = 3.75·k0.

    W_walk = 3.75·k0 · 1.02
    log2(W_walk) = 127.420269 + log2(3.825) = 127.420269 + 1.935525 = 129.355794.

**One-time work (all bounded, charged, negligible).**
* Adding-table build: 2^20 entries, each needs two uniform exponents in [0,r)
  and two exponentiations. Sampling a uniform value in [0,r) is one 256-bit draw
  from a counter-mode PRNG (one block ≈ 400 instr) with rejection probability
  < 2^-128 (at most one retry in expectation); 2^20·2 draws ≤ 2^21·400 ≤ 2^30
  instr ≤ 2^16 units. The 2^20·2 exponentiations (≤ 510 mults each) ≤ 2^30.3
  units. Total table build ≤ 2^30.3 units.
* Start sampling: two uniform draws as above ≤ 2^10 instr, negligible.
* Collision solving at the meeting: one modular inverse of (B_t − B_h) mod r
  (Fermat, ≤ 2·255 mults ≤ 2^9 units) and two multiplications; ≤ 2^9.1 units.
* g = e(P,Q0) ≤ 2^15 units; Q = [a]Q0 ≤ 2^10 units (≤ 510 G2 point operations,
  ≤ 64 Fp multiplications each at the pinned 383 instructions); output encoding
  ≤ 2^−2.

All one-time terms sum to ≤ 2^30.5 units, so

    W = 2^129.355794 + 2^30.5 ≤ 2^129.355794 · (1 + 2^−98.9) ≈ 2^129.355794.

The deterministic iteration cap bounds the walk work; the randomness/rejection,
sampling and collision-solving work are the explicit one-time terms above, not
left to the iteration cap. Submitted bound **time_log2 = 129.40**, rounded up
with ~0.044-bit margin.

## 5. Memory, preprocessing, advice

* Peak memory: the 2^20-entry adding table (each entry: a 576-byte element and
  two 32-byte exponents ≈ 640 bytes) ≈ **2^29.3 bytes**, plus O(1) walk state
  (two elements and two exponent pairs ≈ 1.3 KB). This is far below 2^64 bytes,
  so it is held in a single native 64-bit address space with ordinary
  loads/stores — no segmented, tree, banked, or distributed memory is needed, and
  no such access cost arises. Declare memory_log2_bytes = 30 (reported only; no
  scalar contribution).
* Preprocessing: the table build (≤ 2^30.3 units) is included in W.
* Nonuniform advice: none (0 bytes); the adding table is charged, not supplied.

## 6. Limitations and better directions

* Single heuristic: the R-adding walk behaves like a random mapping on the r
  elements (Teske 2001; toy success ≥ 0.6). Not proven at production size.
* Floyd's 3× factor can be reduced to ~1.5–2× with Brent's cycle detection
  (re-walking to extract the exponent relation), lowering the bound by ~0.5–1
  bit, at the cost of a slightly more involved analysis; left as the next step.
* The order-12 Frobenius automorphism could give a further √12 but reintroduces
  fruitless cycles; a costed special-TNFS/SexTNFS ledger for DLP in Fp12 (~2^120–
  2^126 in those papers' own units) could be lower still.
* Analytic claim; no production-size execution.

References: Pollard (1978); Floyd cycle detection (Knuth, TAOCP vol. 2); Teske,
Math. Comp. 2001.
