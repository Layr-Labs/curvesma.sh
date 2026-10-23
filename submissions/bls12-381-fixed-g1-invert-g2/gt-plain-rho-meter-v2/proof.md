# Fixed-G1 pairing inversion by plain parallel Pollard rho in GT, with a round meter

**Claim.** time_log2 = 12804/100 = 128.04. Given a nonidentity P in G1 and
z in GT, the algorithm below returns Q in G2 with e(P,Q) = z. It succeeds with
probability at least 0.39 over uniform P, uniform z and fresh algorithm coins.
On **every** run, its total work is at most W_max < 2^128.0319 units. One unit
is one dense Fp12 multiplication, C = 20,000 scalar AArch64 instructions.
Peak memory is < 2^76 bytes, spread over 2^16 storage nodes of < 2^62 bytes
each. Preprocessing is < 2^32.1 units and is included. Nonuniform advice: none.

The work bound is deterministic. A round counter, event caps and capped samplers
stop every loop (section 3). Success uses one disclosed random-mapping heuristic,
H1, with a factor-2 margin in the birthday constant (section 4). End-to-end toy
runs of the same walk in prime-order multiplicative groups support it (section 7).

**Design choice.** The walk is a fixed function on GT. It uses no Frobenius or
inversion classes, no escapes and no cycle detector, so two trails that reach the
same element continue identically. The rejected entries on this track failed on
four things, and each is addressed here:

- the partition function and Fp12 representation were undefined (section 2
  defines them exactly);
- restart and verification costs were left out (everything is metered);
- the expected memory was declared as a maximum (memory is capped);
- native addressing beyond 2^64 bytes (the store is split over nodes, each
  with < 2^62 bytes and its own 64-bit address space, and every message is
  charged).

## 1. Reduction

The inputs are validated: canonical encodings, P on the G1 curve and P ≠ O,
[r]P = O, z ≠ 0 and z^r = 1. If z = 1, the attack outputs Q = O, which is valid
exactly when z = 1. Otherwise it computes g = e(P, Q0) with the protected
pairing. Q0 is the fixed generator of G2 and P ≠ O, so non-degeneracy makes g a
generator of GT, a group of prime order r. Then z = g^a for a unique a in
[0, r). The attack finds a and returns Q = [a]Q0, since e(P, [a]Q0) = g^a = z.
It verifies e(P, Q) = z before output. The parameter checks are in
`experiment.py`: p and r probable primes, r | Φ12(p) = p^4 − p^2 + 1, and P0 and
Q0 on their curves.

## 2. Algorithm

**Representation and keys.** GT elements are held as the reference arithmetic
holds them: the blst tower Fp2/Fp6/Fp12, with Fp coefficients in Montgomery
form. Each Fp12 element has exactly one tower representation, and each step's
reference operation is one dense blst Fp12 multiplication. Write κ0(X) for the
canonical integer in [0, p) of the first Fp coefficient of X in this fixed tower
basis. It is obtained by Montgomery export (multiplication by 1) followed by a
conditional subtraction of p. κ0 is a fixed function of the field element and
hence of the group element. Its cost is charged on every step (section 5).

- *Index:* h(X) = κ0(X) mod 2048 (bits 0..10).
- *Distinguished:* bits 11..74 of κ0(X) are zero (θ = 2^-64).
- *Storage node:* bits 75..90 of κ0(X) select one of 2^16 storage nodes.
- *Stored key:* all 12 coefficients, exported canonically. Two keys are equal
  exactly when the elements are equal.

**Bounded samplers.** No sampling loop is unbounded.

- *Scalar draw.* Take four 64-bit generator words, keep the low 255 bits as v
  (r < 2^255), and accept if v < r. Make at most 64 attempts; if all fail,
  halt with failure.
- *Element draw.* Draw scalars u, w and form U = g^u·z^w. Accept if U ≠ 1.
  Make at most 4 attempts; if all fail, halt with failure.

One rejection has probability 1 − r/2^255 < 2^-3.4, so a scalar-draw cap is
reached with probability < 2^-217. U = 1 has probability 1/r per attempt. Over
≤ 2^92 draws, the total sampler-halt probability is < 2^-125. The toy experiment
uses the same rule: draw bitlength(r) bits, with ≤ 64 and ≤ 4 attempts.

**Table.** For h = 0..2047, take M_h = g^{c_h}·z^{d_h} from one element draw,
and store (M_h, c_h, d_h).

**Walk.** F(X) = X·M_{h(X)}. F is a fixed function on GT, determined by g, z and
the table.

**Trails.** M = 2^20 walkers advance in lockstep rounds; each round applies F
once to every walker. A trail starts at U = g^s·z^t from one element draw. Only
(s, t) is stored. After each step, test the new X in this order:

1. If X is distinguished, export its key and send (key, s, t) to the storage
   node chosen by bits 75..90. The trail ends and the walker restarts.
2. Otherwise, if the trail has made L_MAX = 64/θ = 2^70 steps, it is abandoned
   and the walker restarts.

The start element itself is not tested. No inversion is needed in GT, so there
is no batching and no exceptional input.

**Storage nodes.** Each of the 2^16 nodes is a separate machine with its own
64-bit address space. It keeps a red-black tree of its records and accepts at
most 2^52 records. A node that would exceed this halts the attack with failure.
One record is ≤ 1024 bytes: a 576-byte key, s, t, tree links and colour. So a
node uses < 2^62 bytes, natively addressable. On a report, the node looks the
key up:

- if the key is absent, insert the record;
- if the key is present with the same (s, t), ignore the report;
- otherwise, reply with the stored (s', t') and run recovery.

**Recovery.** Replay both trails from their starts with the same F, tracking
A ← A + c_h and B ← B + d_h mod r, so every replayed element equals g^A·z^B.
Both replays end at their first DP, which is the key. If B1 ≠ B2, set
a = (A2 − A1)/(B1 − B2) mod r, compute Q = [a]Q0, verify e(P, Q) = z, and
output Q. If B1 = B2, continue the walk.

## 3. Meter and the deterministic work bound

The attack runs exactly N = ⌈S1/M⌉ + L_TAIL rounds, with L_TAIL = 50/θ, and then
halts with failure. S1 is defined in section 4. Event caps apply; if an event
would exceed its cap, the attack halts with failure before performing it.

| event | cap |
|---|---|
| DP reports | D_CAP = 4·⌈S/2^64⌉ + 2^30 < 2^65.95 |
| trail starts (initial, after DP, after abandon) | M + D_CAP + 2^20 |
| recoveries | 2^20 |
| records per storage node | 2^52 |

S = S1 + M·(L_TAIL + 1) ≥ M·N is the number of steps. Every event is charged
at its worst-case price, which includes capped sampler loops (section 5):

    W_max = S·c_step + D_CAP·c_DP + STARTS·c_start + PRE + 2^20·c_rec + c_final.

The exact rational evaluation in `experiment.py` (mode production) gives:

| term | log2 units |
|---|---|
| walk steps S·c_step | 128.031872 |
| DP reports | 68.66 |
| starts | 77.95 |
| recoveries | 91.11 |
| preprocessing (validation, two pairings, table, init) | 23.0 |
| final (G2 scalar multiplication, pairing, inversion mod r) | 15.2 |
| **total W_max** | **128.0318724** |

The script checks the exact inequality W_max^100 ≤ 2^12804, and it holds. So
every possible execution, including every sequence of sampler rejections and
every early halt, satisfies W ≤ W_max.

## 4. Success probability

**H1 (walk quality).** In the parallel process above, let T be the total number
of steps until some element of GT is produced as a step output by two
**different** trails. Then for t ≤ 2^129:

    Pr[T > t] ≤ exp(−t²/(2κr)),  with κ = 2.

H1 also states that the index, DP and routing fields of κ0 (bits 0..10, 11..74
and 75..90) are pseudo-random functions of the element, so their values on
fresh F-outputs behave as independent uniform coins.

For an ideal random mapping, κ = 1 (van Oorschot–Wiener parallel collision
search). Multiplicative R-adding walks in prime-order groups are the setting of
Teske's experiments (Math. Comp. 2001), which found random-like behaviour for
R ≥ 16. The heuristic correction for R-adding walks is of order 1/R ≈ 0.0005.
The toy runs in section 7 give κ̂ near 1. H1 includes the usual model in which
distinguishedness of fresh F-outputs behaves like θ-coins.

**Main phase.** S1 = ⌈√(2κr·0.51083)⌉ + 2, so log2 S1 = 127.9440 (log2 r =
254.857). Since ln(1/0.6) < 0.51083, Pr[T ≤ S1] ≥ 0.4.

**From the first repeat to output.** Let X be produced first by trail A and
later (or later in the same round) by trail B ≠ A.

1. **Merging.** From X onward both trails apply F and visit the same elements up
   to A's first DP D after X. B reaches D the same number of steps after X, so
   A reports D first. Reports within a round are processed sequentially, so B
   finds A's record at D's node.
2. **Timing.** Pr[B needs more than 50/θ further steps] ≤ e^-50. The tail gives
   every walker 50/θ rounds after the main phase.
3. **No abandonment.** A trail is abandoned only after 64/θ steps without a DP,
   with probability ≤ e^-64 per trail. A structural cycle of F of length k
   requires Π M_h^{n_h} = 1 with n_h ≥ 0 and Σ n_h = k. That is,
   Σ n_h·(c_h + d_h·a) ≡ 0 (mod r). For each fixed relation, uniform c_h makes
   this hold with probability exactly 1/r. A union bound over all k ≤ 28 gives
   < 2^-44.5 (computed in `experiment.py`). Longer self-cycles are ordinary
   self-collisions under H1, and their expected number over all trails is
   negligible (≲ 2^(90+92−254)).
4. **Distinct coefficients.** Conditioned on every group element (g, z, the
   table elements, the starts, and hence every walk), each trail's t remains
   uniform and independent: given U, t determines s. Each d_h is independent of
   M_h because c_h is uniform. The replayed B equals t plus a sum of d_h along
   a path fixed by the group elements, so B1 = B2 has probability exactly 1/r.
   Identical starts have probability ≤ 2^180/2^254.
5. **Caps.** D_CAP is at least 4 times the expected number of DPs. Under H1
   the routing field (bits 75..90 of κ0, disjoint from the index bits 0..10 and
   the DP bits 11..74) is a pseudo-random function of the element. A node's
   share of the n ≤ D_CAP reports is then Binomial(n, 2^-16) with mean
   < 2^47.95. Its cap is 2^52 records, about 17 times the mean. Chernoff gives
   < exp(−2^50) per node and < exp(−2^49) over all nodes. The toy pipeline
   routes reports to 16 nodes with a cap of only twice the expected share and
   records 0 overflows (section 7). The
   recovery and start caps are far above their expected values. Exceeding any
   cap has negligible probability under the DP model.

Hence Pr[success] ≥ 0.4 − e^-50 − 2^-44 − 2^-60 − 2^-125 > 0.399 ≥ 0.39. When
z = 1, the answer is found directly.

## 5. Prices

**Per step.** One dense Fp12 multiplication X·M_h is the reference operation:
1 unit, including its own operand loads and stores (it was measured on GT
elements in memory). On top of it, in instructions:

- **Key export, FP = 1000.** One Montgomery multiplication by 1 plus a
  conditional subtraction. Itemised: the 6×6-limb product is 36 MUL/UMULH pairs
  at ≤ 7 instructions each (252); 6 reduction rounds at ≤ 30 each (180); the
  final subtraction and selection ≤ 30. That is ≤ 462. The protected
  measurement gives 383 for blst_fp_mul.
- **Control, CTRL = 256.** The itemised count is about 40: index AND and table
  address (3), DP test (3), node test (2), trail counter (3), round and loop
  control and meter compares (8), and the table-pointer load (2). The allowance
  also covers storing the result pointer.

c_step = 1 + (1000 + 256)/20000 = 1.0628 units.

**DP report, c_DP = 2^17 instructions.**

| part | instructions |
|---|---|
| export of all 12 coefficients | 12·1000 |
| record build | ≤ 128 |
| message send and receive, network-stack instructions on both sides | ≤ 2^14 |
| node-local red-black tree search and insert (height ≤ 2·52 + 2 levels, a 72-limb compare per level ≤ 160 instructions, plus rebalancing) | ≤ 2^15 |
| reply | ≤ 2^12 |
| **sum** | < 2^16 |

**Group operations used outside the walk** (each a worst case):

- **GT exponentiation by a scalar < 2^256: 512 units.** ≤ 256 squarings and
  ≤ 256 multiplications, each charged as a dense multiplication.
- **Element draw: ≤ 4 attempts.** Each attempt is 2 exponentiations, 1
  multiplication, an identity test and 2 scalar draws. One scalar draw is ≤ 64
  attempts × (4 words × 64 instructions + 16).
- **c_start = 4·(1025 + small) units.**
- **Pairing: 2^15 units allowance.**
  - Miller loop over |x| = 0xd201000000010000, which has 64 bits and Hamming
    weight 6. Each of ≤ 70 doubling or addition steps costs ≤ 1 Fp12 squaring,
    ≤ 1 line multiplication, and the G2 step with line coefficients (≤ 30 Fp2
    multiplications, i.e. ≤ 120 Fp multiplications at ≤ 1000 instructions, = 6
    units). That is ≤ 8 units per step and ≤ 560 in total.
  - Final exponentiation. The easy part (conjugation, one Fp12 inversion
    through the tower norms using ≤ 40 dense-equivalents plus one Fp inversion
    of < 2^20 instructions, a Frobenius and 2 multiplications) is ≤ 100 units.
    The hard part, an exponent (p^4 − p^2 + 1)/r of ≤ 1270 bits, is ≤ 2540 units
    by square-and-multiply.
  - The total is < 2^12, so the allowance is 8 times the itemised count. Two
    pairings are charged: g and the final verification.
- **Input validation.** [r]P in G1 costs ≤ 510 Jacobian operations at ≤ 16
  Fp multiplications each, < 2^11 units. z^r costs 512 units.
- **Output.** [a]Q0 in G2 costs ≤ 512 Jacobian operations at ≤ 30 Fp2
  multiplications each (≤ 3072 units), plus an Fp2 inversion: < 2^12 units.
  The inversion mod r costs 2^19 instructions.
- **Replay step:** 1 + (1000 + 256 + 4·64)/20000 units (two additions mod r).
  **c_rec** = 2·(c_start + (L_MAX + 2)·replay step).
- **PRE** = 2 pairings + [r]P + z^r + 2048·c_start + 2^16 nodes × 2^20
  instructions of node set-up (process start, allocation, routing table) +
  2^30 instructions for decoding and initialisation. With the first 2^20 starts,
  preprocessing is < 2^32.1 units.

## 6. Memory

There are ≤ D_CAP < 2^65.95 records of ≤ 1024 bytes, so < 2^75.95 bytes in
total. Each node stores < 2^62 bytes. Walker state is ≤ 2^20 × 640 bytes, and
the table is 2048 × 640 bytes. Peak total < 2^76 bytes, so
memory_log2_bytes = 76. Memory has no score contribution.

## 7. Evidence (`experiment.py`)

**Production mode** checks the parameters, the exact rational ledger, the score
inequality and the short-relation union bound.

**Toy mode** runs the same walk in the order-r subgroup of F_q^*, for primes
q = k·r + 1 with r ≈ 2^20..2^32. Elements are canonical integers, so κ0 is the
integer itself. It uses the same index and DP bit layout (index = low 11 bits,
DP = the next θ bits, routing = the next 4 bits to 16 nodes with a cap of
twice the expected per-node share, overflow = failure), the capped samplers,
restarts, L_MAX, the store, replay
recovery and verification. A full visited-element map gives the exact first
inter-trail repeat T, and κ̂ = mean(T²)/(2r). The budget S1 uses κ = 2, as in
production. The initial run (r ≈ 2^20, 40 trials) takes about 2 s.

The offline rows below come from the same script with the inputs listed in
README.md. "Merge within S1" means a verified recovery whose first inter-trail
repeat occurred within S1. The proof needs 0.4; for κ = 1 the model value is
about 0.64. The tail is long compared with S1 in toy runs, so almost every trial
eventually recovers a.

| log2 r | R | walkers | θ | trials | verified recovery | merge within S1 | κ̂ | κ̂ excl. trapped-trail steps | trapped trails | wrong | B1=B2 | node overflows | mean max/mean node load |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 20 | 2048 | 16 | 2^-4 | 400 | 400/400 | 0.695 | 0.897 | 0.895 | 9 | 0 | 0 | 0 | 2.10 |
| 24 | 2048 | 32 | 2^-5 | 300 | 300/300 | 0.653 | 1.006 | 1.006 | 2 | 0 | 0 | 0 | 1.73 |
| 24 | 16 | 32 | 2^-5 | 300 | 300/300 | 0.583 | 1.097 | 1.097 | 1 | 0 | 0 | 0 | 1.70 |
| 28 | 2048 | 128 | 2^-6 | 200 | 200/200 | 0.655 | 0.913 | 0.913 | 0 | 0 | 0 | 0 | 1.47 |
| 28 | 2048 | 128 | 2^-6 | 200 | 200/200 | 0.670 | 0.878 | 0.878 | 0 | 0 | 0 | 0 | 1.47 |
| 32 | 2048 | 256 | 2^-6 | 80 | 80/80 | 0.713 | 1.050 | 1.050 | 0 | 0 | 0 | 0 | 1.24 |
| 32 | 2048 | 256 | 2^-6 | 80 | 80/80 | 0.562 | 1.117 | 1.117 | 0 | 0 | 0 | 0 | 1.23 |

At r ≈ 2^20 with R = 2048, R² > r, so short table relations Π M_h^{n_h} = 1
exist. Some trails become trapped, and they are counted in the table. The
union bound in section 4 excludes such relations up to length 28 in production.

**Limitations.** The toy groups are F_q^* subgroups, not GT. The walk only uses
group multiplication, a fixed key function and uniform table elements, so the
abstract process is the same. Nothing at 2^128 scale is tested; the claim rests
on H1 with a factor-2 margin.

## 8. Score

time_log2 = 12804/100. W_max < 2^128.03188 ≤ 2^128.04, success ≥ 0.399,
memory < 2^100 bytes, preprocessing < 2^32.1 units, advice 0 bytes.
