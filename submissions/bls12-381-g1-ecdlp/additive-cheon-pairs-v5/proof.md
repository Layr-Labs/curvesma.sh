# bls12-381-g1-ecdlp: additive Cheon with planned pair walks and two-hop bucket matching (additive-cheon-pairs-v3)

**Claimed time_log2 = 116.485098** (exact rational 58242549/500000 in lowest terms; an upper
bound on log2 W(2^20), recomputed from exact integers by `evidence.py` mode `ledger`, which also
asserts the production algebra identities below). Heuristic success at least 0.398 (required
0.39): 1 - e^{-K} = 0.399504 at the realized K = 6N^2/m >= 0.51, minus the 0.0005 coverage
allowance of H1, a 2^-20 modelling budget, a 3*2^-40 cap budget (the union of the cap tails
below) and the 2^-112.9 diagonal event, leaving 0.399003, claimed 0.398. Peak memory
2^122.6439 bytes with every node counted at its full itemised allocation. No advice.
Preprocessing at most 2^115.5030 units, included in W. Every cost, size and failure figure is
rounded up, every success figure rounded down, by the script.

This is the construction of the rejected submission `additive-cheon-pairs-v1` (20a8ca63) with
every defect of its review repaired, and no premise inherited from it: the pair operation is
priced 2000 + 3*383 = 3149 instructions per two points plus charged second-formula arithmetic
(the reviewer's own accounting, replacing the invalid 9M = 3000 scaling); the candidate-defined
global-address machine is replaced by an explicit bounded cluster model with no
remote-memory instruction, adapted from the bucket-matching architecture of the competing
frontier family (whose submissions' published accounting I mirror; no acceptance of theirs is
used as a premise here); the coverage argument no longer excludes modular coincidences by
integer size — H1 below is an explicit heuristic supported by an exact census of the
construction on toy primes, in which every modular and rational coincidence is measured, not
excluded; the experiment selects primes deterministically in increasing order and emits exactly
the sets the probability statement uses; and no prior submission is used as a proof premise.
The useful core of the earlier work is kept: cheap point generation from planned three-term
progressions with one shared inversion per two points, and the six-element symmetry key
min(y, p-y).

Cost accounting. Work is the sum of instructions executed over all processors
(tracks/docs/COST_MODEL.md, classical-group-work-v2). Cluster model, stated so that every charge
below has a data path: nodes with private memory and no remote-memory instruction; the only
communication is posting a 512-byte block into another node's receive FIFO; filling the block and
the post are charged to the sender, the poll and the read to the receiver; a node with nothing to
do executes no instructions. Wall time is not claimed.

## 1. Problem and reduction

BLS12-381 G1, order r, powers A_i = [s^i]P0, i = 1..n. At the scoring n = 2^20 let d = 2^20.
Since v_2(r-1) = 32 (asserted), d | r-1 and 6d | r-1, so m = (r-1)/d = 2^234.857089... is a
multiple of 6 (asserted). H = {x^d} has order m, t := s^d is uniform on H and A_d = [t]P0 is an
input. Only T = A_1 and A_d are read, so the same W applies to every n >= 2^20.

Let lambda = x^2 - 1 mod r with x the curve parameter; then lambda^2 + lambda + 1 = 0 mod r and
lambda^m = 1 (both asserted), so U = {+-1, +-lambda, +-lambda^2} is the sixth-root subgroup and
lies in H. On this curve [lambda](x, y) = (beta x, y) and [-1](x, y) = (x, -y), so the six
points [u]Q, u in U, have y-coordinates +-y(Q), and min(y, p-y) is an exact key for the U-orbit
of any finite point. The production implementation stores Montgomery y = y*2^384 mod p; on
stored residues the negation map is Y -> p - Y, so min(Y, p - Y) and its low 192 bits are an
exact orbit key of the stored value with no normal-form conversion. The toy exercises the min
rule on affine y; the stored-form rule follows from the residue identity, which is the only
claim made for it.

Recovery of s from t is the standard two-step, both charged: (i) twenty Tonelli-Shanks square
roots in F_r give one s0 with s0^d = t; (ii) s = s0*omega^i for one of the d roots of unity
omega = g^((r-1)/d); write i = i0 + 2^10 i1, compute the 2^10 points [s0 omega^{i0}]P0 and the
2^10 points [omega^{-2^10 i1}]T, sort the first by x and look the second up in it, one scalar
multiplication each. Charged at 2^24 units per run in the ledger (the listing: 20 square roots
at 2369 F_r multiplications each, 2^11 scalar multiplications at 2^12 units, 3*2^10 products,
2^11 comparisons of 48-byte x-coordinates at 64 instructions).

## 2. Algorithm

Parameters: K = 0.51; N per side = the least integer with 6N^2 >= K m, N = 2^115.6503; X2 = 2^76
(the class integers b); GMAX = 2^44 = 16 d^2 (the pair-gap and table bound); base interval
I = (2^76, 2^76 + 2^70]; f = 0.98, the claimed minimum pair fraction (H4).

1. **Characters and class lists** (challenge independent). chi : F_r^* -> Z/d, chi(x) = x^m
   projected by a fixed generator of mu_d, has kernel H (and chi(2) != 0: 2^m has order 2^16 on
   this r; the script asserts 2^m != 1). Sieve the primes q <= X2, compute chi(q) for each (one
   235-bit exponentiation, 366 multiplications at 383 instructions, plus a Pohlig-Hellman
   logarithm in mu_{2^20} (bit-by-bit, about 400 squarings), charged together at 400000 instructions =
   200 units each), then chi(b) for all b <= X2 by multiplicativity into a 20-bit array. The
   class list L_c = {b in [1, X2] : chi(b) = c}, in increasing order, has about X2/d = 2^56
   elements (H2); the exact lengths are computed by the sieve.
2. **Planning** (challenge independent, section 3). For each class c a plan produces every
   element of L_c exactly once: each element is a root (one addition from the previous root) or
   one child of a pair (x-g, x, x+g) with midpoint x and gap g <= GMAX.
3. **Base primes**. Sieve I and compute chi(a) for its primes in increasing order. The baby side
   takes them in order until its N elements are emitted; the giant side continues with the next
   primes. A base a of character chi(a) walks the plan of class -chi(a), so every emitted scalar
   a*b lies in H. Baby and giant base sets are disjoint, and a > X2 >= b, so a is the unique
   prime factor of a*b above X2: (a, b) is recovered from the product by one factorisation.
4. **Baby walk.** For each base a in order: Q_a = [a]P0 by scalar multiplication; the table
   T[k] = [k]Q_a for 2 <= k <= GMAX is built in 2^16 strands by repeated addition (T[2] by one
   doubling, strand starts by scalar multiplication); the plan is executed, emitting the record
   (192-bit truncated key, 64-bit chain id) for every point. The whole baby side completes before
   any giant chain starts.
5. **Giant walk.** The same with Q_{a'} = [a']A_d, emitting nothing: each giant record is routed
   to the owner of its key and compared against the stored baby records there.
6. **Collision.** A truncated-key match is a candidate (key, id, id'), appended at the owner node.
   After the giant side, the status tree collects candidate counts; if the global count exceeds
   2^44 the run outputs fail. Otherwise each candidate's two chains are re-walked once, linearly,
   from their starts by their generating nodes (tables rebuilt, all charged), collecting the at
   most 4 positions per chain whose key equals the retained key; a fifth is a failure. The at
   most 16 pairs of kept positions give at most 96 values z = u*b/g, each tested by [z]P0 = A_d.
   On a hit, section 1 recovers s. A wrong candidate never yields a wrong answer: every scalar
   is verified against A_d and then T before output.

**The pair operation.** For P = (x1, y1) and T = (x2, y2) with x1 != x2, both P+T and P-T need
the single inverse delta = (x2 - x1)^{-1}, provided by Montgomery's trick shared across the
round. Beyond one reference addition (6M = 2000 instructions: 3M formula + 3M batch share, with
the reference's own field additions inside the A = 0 normalization), the second formula needs
three field multiplication/squaring operations, priced at the supplied 383-instruction count:
**pair arithmetic = 2000 + 3*383 = 3149 instructions for two points**. The second formula's
field additions and subtractions are not part of the reference's A = 0 convention and are
charged separately (section 5).

**Incomplete additions are safe.** Every point is a multiple of Q = Q_a != O. For a pair with
midpoint x and gap g: the children are list elements, so 0 < g < x and x + g <= X2 < r; neither
x = g nor x = -g mod r. For a root step u -> v = u + gap with gap <= GMAX: u = -gap would force
v = 0; u = gap means v = 2u, but chi(2) != 0, so chi(v) != chi(u) and v would not be in the
class. The table's k = 2 step is the one doubling. The giant side is the same with
Q = [a't]P0 != O. The toy asserts a non-zero denominator at every step.

## 3. The planner

Fix a class list L = (b_0 < ... < b_{l-1}), l ~ 2^56. (i, j), i < j, is a pair with midpoint k if
b_k - b_i = b_j - b_k = g, 1 <= g <= GMAX. **Edges:** for each u in L, the class members v in
(u, u + 2 GMAX] with u = v mod 2 whose midpoint (u+v)/2 is in L; under H2 an element lies in
about 2 GMAX/d^2 = 32 such pairs as a child. **Matching:** greedy by increasing degree, each
pair accepted only if it creates no cycle (parent chain walk, at most D = 2^14 steps) and no
depth overflow (depth(k) + 1 + max(height(i), height(j)) <= D with never-lowered heights, so the
check only over-rejects). No augmenting paths and no unmatching during matching; only the final
root-chain pass demotes (an element more than GMAX past the last root has its predecessor
demoted), so depths only increase during matching (at most l*D updates) and only decrease after.
Consecutive roots are at most GMAX apart. **Emission:** list order; each root, then its subtree
depth-first, so the midpoint of a pair is always on the current path; the ancestor stack holds at
most 2D + 1 = 2^15 + 1 points, asserted by the toy census. A plan entry is 32 bytes: kind, gap
(45 bits, covering g = GMAX), source and two destination slots (16 bits each, covering the
2^15 + 1 stack), and the list indices of the one or two elements produced (57 bits each); the
spare bits cover control. The plan of class c is challenge independent, shared by every base of
character -c on both sides, and is streamed entry by entry from a storage node through the block
protocol into a ring at the generating lane that uses it.

The planner cost (window scan 2^24 tests per element at 16 instructions, matching with bounded
ancestor checks, emission 64) is a ledger row at 2^93.037 units — ten orders below the records
row — and every access names its data path (section 5). The pair fraction f_c is known before
any group operation; the claim uses f = 0.98 for every class (H4).

## 4. Matching and transport

**Owner nodes.** 2^56 owner nodes, indexed by key bits 0..55, read as two 28-bit fields (t, u).
A generating node scatters its records into 2^28 send blocks of 16 records (512 B) by t and posts
full blocks to the hop-1 relay (t, 0); that relay re-scatters by u and posts to the owner (t, u).
Two hops; each node's send array is 2^28 blocks = 2^37 B. Within an owner, baby records are
inserted into a direct-mapped head array of 2^62 slots indexed by key bits 56..117, chained
through a pool of 2^61 32-byte records (2.58 x the mean load 2^59.63): a head slot is an 8-byte
link, and the slot counts (for the at most 64 records per slot check) live in a separate byte
array, one byte per slot; a pool record stores key bits 118..191 (74 bits), the 64-bit chain id
and a 61-bit link. Giant records are not stored: each reads its head slot and compares its 74
remaining key bits with each pooled record on the chain (7 instructions per comparison: 24-byte
load, compare, branch, strided advance, counter decrement). An equal key appends (full 192-bit
key, both chain ids, 40 B) to the node's candidate list of 2^16 slots; a 2^16+1st candidate fails
the node. Comparisons are charged per node up to 1.01 x the mean (2^57.31), all 2^56 nodes: the
row is 2^105.157 units.

**Cluster inventory.** 2^20 generating nodes, each with 2^16 lanes in lockstep (one batch
inversion per round of 2^16 lane operations, all lanes local to the node); a lane holds one
base's table (2^44 affine points, 96 B each, 2^50.6 B) and its ancestor stack, so a generating
node is 2^67 B, plus a candidate receive area of 2^24 40-byte slots (mean forwarded load 2^20.3 per node; a
2^24+1st candidate fails the node, a tail inside the 3*2^-40 budget). 2^56 owner nodes at 2^66.644 B (head links
2^65 + 2^62, pool 2^66, send arrays and FIFOs, candidate list, code). 2^30 storage nodes holding
the class lists (10 B per b <= X2), the plans (32 B per entry, 0.51 entries per element) and the
sieving primes, at their ceil shares of the 2^80.5 B total. The status tree is a fan-in-256 tree
over all nodes, three passes (after baby, after giant, after recovery), 2^12 instructions per
node per pass. Ordering: baby insertion completes before giant probing begins (enforced by the
first status pass), candidate forwarding and re-walks happen after the second pass, and the
global cap check precedes any forwarding.

Every queue is bounded (one send block in flight per destination per credit, credits returned with
the block's completion), every frame is a fixed 512 B, and no node polls: a node with an empty
FIFO executes nothing until the next post. Partial blocks are flushed at the end of each side,
charged. The routing graph is acyclic (key digits strictly advance per hop), so no deadlock is
possible; the 3*2^-40 budget covers the union of queue, termination and cap tails (the leading
terms: comparisons 2^-40, re-walk repeats 2^-45, false matches 2^-47).

## 5. Instruction listing per point (AArch64 scalar, 64-bit registers)

**Pair op, two points** (the walk items beyond the fixed 3149 arithmetic; listed 190):

| item | instructions |
| --- | ---: |
| plan entry decode (4 loads, kind, gap, slots, indices) | 8 |
| load midpoint P (6 ldp, 2 address) | 8 |
| load T[g] (6 ldp, 2 address, 3 index) | 12 |
| batch-inversion lane buffers (6 per record, 2 records) | 12 |
| second formula's modular adds/subs: 4 six-limb operations at ~12 each with borrow correction | 48 |
| store two points (12 stp, 2 address) | 14 |
| F_p inversion shared by the round, 232864/2^16 | 4 |
| lane control and loop | 10 |
| key extraction, 2 x 25 (p - y: 13, compare and select: 12) | 50 |
| **listed** | **166** |
| charged | **190** |

The 48 is listed for four six-limb operations (negation-sum, two subtractions with conditional
correction, final subtraction); charged 72 within the 190 as margin for operand spills.

**Root op, one point** (gap stream read 5, table address 2, T[gap] load 12, load/store point 14,
batch buffers 6, inversion 4, loop 5, key 25; listed 73, charged 130).

**Record transport, per record:** scatter at 2 hops (block index, slot, bounds) 24; post and poll
of full blocks amortised over 16 records, 2 hops, 6; transport copy of 32 bytes at 2 hops
(4 ldp/stp pairs plus addressing per hop) 12; insert (head load, count check, link write, record
store, head update) 10 at the baby side, lookup fixed cost 6 at the giant side. Listed 50,
charged 58.

**Plan streaming**, per entry: storage-node read and post 6, transport copy 6, ring write at the
lane 3 (the lane's read of the entry is the 8 of the pair/decode listings); charged 15 per entry,
0.51 entries per point.

**Gap tables** are charged by their own ledger rows (2^60.65 bases x 2^44 additions at 2000 + 100,
plus one inversion per round of 2^16 strand additions and the strand starts), not amortised into
the records row.

**Per point, both sides.** Arithmetic f*3149/2 + (1-f)*2000 = 1583.01 at f = 0.98. Everything
else (walk bookkeeping 95.7 — which already contains the key extraction — record transport 50,
plan stream 7.65) listed 153.35, charged x1.15 = 176.35. **Total 1759.36 instructions per point;
the records row is 2N x 1759.36.**

## 6. Ledger (units of 2000 instructions; every row an exact integer ceiling; the script recomputes each and asserts the totals)

| term | log2 units |
| --- | ---: |
| records: 2N x 1759.36 instructions | 116.465 |
| owner init (head and count zeroing) | 110.204 |
| comparisons | 105.157 |
| gap tables | 104.721 |
| candidate re-walks | 101.073 |
| gap-table inversions | 95.514 |
| planning | 93.037 |
| table strand starts | 88.650 |
| chi of primes <= X2 | 78.252 |
| per-base setup | 74.235 |
| chi sieve + class lists | 72.678 |
| chi of base primes | 72.252 |
| plan storage write | 71.034 |
| base window sieve | 66.034 |
| candidate z tests | 62.586 |
| status tree | 58.619 |
| Lemma 1 + io | 44.000 |
| sieving primes sieve | 43.122 |
| candidate forwarding | 40.034 |
| **W** | **116.485098** |

W = 2^116.4850979... units; the score is ceil(10^6 log2 W)/10^6 = 116.485098 = 58242549/500000
(the script computes log2 W to 80 decimal digits and asserts the 10^-6 grid point is more than
10^-20 away, so the ceiling is exact). For context, the same budgets with every point computed as
a root (no pair ops) give 2^116.813835; the pair trick is worth 0.329 bits. Preprocessing (the
baby half of the records row, the baby tables, the lists, characters, plans, sieves, owner init
and the first status pass) is at most 2^115.5030 units, included. Memory, every node at its full
itemised allocation: 2^56 owner nodes at 2^66.644 B (head links 2^65 and count bytes 2^62, pool 2^66,
send arrays and FIFOs, candidate list, code), 2^20 generating nodes at 2^67 B plus candidate
receive areas, and 2^30 storage nodes holding lists, plans and sieving primes at their shares:
2^122.6439 bytes in all. Charge accounting: the key extraction is inside the pair/root
bookkeeping of section 5 (nowhere else); the gap tables, their inversions and strand starts are
only the four table rows; the walk inversions are the per-point share (4 per pair op) inside the
records row; the re-walk inversions are inside the re-walk row; candidate forwarding and
re-walking are separate from the records row. Each instruction is charged exactly once across
these categories.

## 7. Success probability

Steps 1-3 of section 2 are challenge independent; the run fails and outputs no scalar if the base
searches come up short (H2) or the plan breaks the gap bound (never observed; the walk checks its
gaps and stops before any chain starts). The algorithm succeeds exactly when t lies in
R = U * {b/g} over the emitted sets B, G (|B| = |G| = N) and no cap breaks. Because t is uniform
on H and R is a fixed set (the lists do not depend on the challenge), Pr[t in R] = |R|/m exactly.

**H1 (coverage).** |R|/m at the production construction is at least 1 - e^{-K} - 0.0005, with
K = 6N^2/m >= 0.51, so at least 0.399004. Support:

- *Exact census of the construction.* `evidence.py` mode `census` builds the construction on toy
  primes r = 6 c 2^k + 1 (c odd, smooth, deterministic; skips r with 2 in H_t) at the predeclared
  scaling N = ceil(sqrt(K m/6)), class-size target e = round(N^0.4843), X2 = d e, bases = the
  primes of (X2, X2 + X2^0.921] in order, each walking its full class list, the last base partial
  so each side emits exactly N — the production shape. It computes chi for every b <= X2 by a
  Pohlig-Hellman dlog sieve, forms the emitted sets, and marks the U-coset of every ratio
  (a b)/(a' b') in a byte array of m/6: the count is exact, with every modular and rational
  coincidence included, no sampling. The marking was validated against direct field arithmetic
  (an independent brute-force reconstruction of R) on 18 primes at m up to 2^18.5, agreeing to
  10^-9. The census rows:

| toy m | d_t = 2^k | primes | mean diff | lowest | repeat frac (uniform value) |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 2^24.0 | 2^4 | 8 | -0.028615 | -0.059002 | 0.273 (0.217) |
| 2^24.0 | 2^8 | 8 | -0.002269 | -0.007363 | 0.221 (0.217) |
| 2^24.0 | 2^12 | 8 | -0.000201 | -0.000834 | 0.217 (0.217) |
| 2^24.0 | 2^16 | 8 | -0.000051 | -0.000381 | 0.217 (0.217) |
| 2^26.0 | 2^4 | 8 | -0.029485 | -0.043529 | 0.275 (0.217) |
| 2^26.0 | 2^8 | 8 | -0.003670 | -0.008165 | 0.224 (0.217) |
| 2^26.0 | 2^12 | 8 | -0.000148 | -0.000335 | 0.217 (0.217) |
| 2^26.0 | 2^16 | 8 | +0.000050 | -0.000028 | 0.217 (0.217) |
| 2^28.0 | 2^4 | 8 | -0.026865 | -0.050595 | 0.269 (0.217) |
| 2^28.0 | 2^8 | 8 | -0.003045 | -0.013193 | 0.223 (0.217) |
| 2^28.0 | 2^12 | 8 | -0.000240 | -0.000412 | 0.217 (0.217) |
| 2^28.0 | 2^16 | 8 | +0.000038 | -0.000003 | 0.217 (0.217) |
| 2^30.0 | 2^4 | 8 | -0.027476 | -0.045332 | 0.271 (0.217) |
| 2^30.0 | 2^8 | 8 | -0.008283 | -0.037778 | 0.233 (0.217) |
| 2^30.0 | 2^12 | 8 | -0.000194 | -0.000411 | 0.217 (0.217) |
| 2^30.0 | 2^16 | 8 | -0.000011 | -0.000068 | 0.217 (0.217) |
| 2^32.0 | 2^4 | 8 | -0.025991 | -0.031271 | 0.268 (0.217) |
| 2^32.0 | 2^8 | 8 | -0.007463 | -0.036635 | 0.231 (0.217) |
| 2^32.0 | 2^12 | 8 | -0.000236 | -0.000414 | 0.217 (0.217) |
| 2^32.0 | 2^16 | 8 | -0.000012 | -0.000034 | 0.217 (0.217) |
Top-window control (products a*b near r, maximal modular wraparound):

| toy m | d_t | primes | mean diff | lowest |
| ---: | ---: | ---: | ---: | ---: |
| 2^26.0 | 2^8 | 4 | -0.002824 | -0.004590 |
| 2^26.0 | 2^12 | 4 | -0.000289 | -0.000403 |
| 2^28.0 | 2^8 | 4 | -0.001502 | -0.001665 |
| 2^28.0 | 2^12 | 4 | -0.000248 | -0.000441 |
| 2^30.0 | 2^12 | 4 | -0.000196 | -0.000219 |

  The deficit shrinks with the toy degree. Empirically, across the measured rows the mean deficit
  is roughly proportional to ln(X2)/d_t — an empirical ansatz, not a derivation: the k = 4 rows
  (driver about 0.35-0.4) measure means -0.026 to -0.029 (worst -0.059); the k = 8 rows (driver
  about 0.04) means -0.002 to -0.008 with worst rows to -0.038 (the small-d rows drift with m;
  every k >= 12 row is flat in m); the k = 12 rows (driver about 0.003) have means -0.00015 to
  -0.00024 and worst -0.00083, at the smallest m; the k = 16 rows (driver about 0.0002, four
  times production's) have means within +-0.00005 of zero and worst -0.00038. Production has
  d = 2^20 with driver ln(2^76)/2^20 = 5.0e-5 per pair, 60 times below the k = 12 driver and
  4.4 times below the k = 16 driver. A *top-window control* reruns the census with the base
  window placed at (r/X2, r/X2 + (r/X2)^0.921] instead, so that the products a*b are near r and
  the ratio cross-products wrap modulo r as much as they can; its rows are reported alongside and
  sit at the same deficits as the low-window rows at the same (m, k), so the measured trend does
  not depend on the window height. The analytic bound on the production rational-repeat rate
  (E[gcd]-based, 2*(6/pi^2)(ln X2 + 0.15)/d < 1.1e-4 per pair) is consistent with the ansatz; it
  bounds the literal-equality repeats, which unique factorisation does constrain, and no modular
  coincidence is excluded anywhere. The repeat-fraction diagnostic shows the k >= 12 maps are
  uniform-like: the measured fraction of pairs hitting an already-marked coset is within 0.01 of
  the value 1 - (1 - e^{-K})/K that independent uniform draws would give (0.2167 at K = 0.51),
  while the k = 4 rows sit well above it. H1 remains a heuristic: the census reaches m = 2^32 at
  k = 12 (X2 = 2^18.9, against production m = 2^234.86, X2 = 2^76, and products a*b = 2^152 below
  r = 2^254.9 while the toy products sit near or below r depending on the window), so the
  production statement is the explicit 0.0005 allowance. Driver-scaled to production (the ansatz
  applied to the worst row of every k), no measured worst row exceeds -0.0001, five times inside
  the allowance; the claim's 0.398 tolerates a coverage deficit up to 0.0015 before the other
  terms, and the required 0.39 is robust to a shortfall of 0.0095, nineteen times the allowance.
- *End-to-end agreement.* The ec toy of section 8 checks, every trial, whether the sampled t lies
  in the exactly-enumerated R of its emitted sets, and every observed candidate/miss event agrees.

**H2 (equidistribution).** The characters of the primes of I and of the integers b <= X2 spread
evenly over Z/d, so class lists have about 2^56 elements and the walks need about 2^59.65 bases
per side. The interval I = (2^76, 2^76 + 2^70] holds about 2^64.3 primes by the PNT in intervals
(Rosser-Schoenfeld bounds; the sieving charge covers all 2^70 + 2^43 integers), a 12 x margin
over the 2^60.65 bases charged. Exact class lengths and base counts are computed by the sieves at
runtime; a shortage stops the run before any chain starts and is a failure event under H2, which
excludes it. Toy evidence: the class lengths of the ec toy's 32 classes and the census's class
counts at every scale. H2 affects only charged small terms, never the sets.

**H3 (keys).** The 192-bit keys of the 2N records behave as independent uniform strings across
U-orbits: sharding and head occupancy are binomial, and a cross-orbit key equality has
probability 2^-192. Support: the toy's occupancy chi-square matches the binomial expectation
(scaled statistic about 1), and its comparison and false-match counts match the model. Tails
under H3, each already multiplied by the number of nodes or slots: a head slot above 64 records
is Poisson(0.196) below 2^-340 over all 2^118 slots; the node baby-count Chernoff and the
comparison-count Bernstein bounds are below 2^-40 over all nodes (the 1.01 cap is charged as the
maximum); false key matches have mean N^2/2^192 = 2^39.26 against the 2^44 cap, Chebyshev below
2^-47 under H3; a re-walk chain carrying one key twice is L^2/2^193 < 2^-81, and one key at more
than 4 steps is below 2^-45 over all chains; the generating nodes' candidate receive areas
(mean forwarded load 2^20.3 per node against the 2^24 slots, cap-checked) carry tails inside the
same budget. The success line subtracts 3*2^-40 for the union of these.

**Diagonal event.** If t = u a/a' for base primes a, a' of the same character (chi(a) = chi(a')),
then at every common position b the baby and giant points share an orbit, so one base pair gives
more than 2^44 true matches and the run outputs fail at the candidate-count check (which precedes
any forwarding). These t form at most 6 (2^59.65)^2 = 2^121.9 elements of H, so the event has
probability below 2^-112.9. Off the diagonal, true matches are Poisson-like under H1 with mean at
most K; the 2^20 true-match budget inside the 2^44 cap and the per-node 2^16 candidate slots
carry tails inside the 3*2^-40 budget.

**Claimed success, under H1-H3:** 0.399504 - 0.0005 - 2^-20 - 3*2^-40 - 2^-112.9 = 0.399003,
claimed **0.398** (margin 0.008 over the required 0.39). The failure budget covers: base-search
shortage (H2), no collision (1 - e^{-K} and the H1 allowance), sharding and comparison caps
(H3), false truncated-key matches and candidate caps (H3), diagonal events, re-walk key repeats
and cap overflows, queue and termination tails. Every cap, abort and failure path is charged at
its maximum in the ledger.

## 8. Evidence (`evidence.py`, standard library)

`ledger`: every number of section 6 from exact integers, the exact score with the boundary
assert, the preprocessing and memory inequalities, the success inequality with every failure term,
and the production algebra identities (v_2(r-1) = 32, 6 | m, lambda = x^2 - 1 with
lambda^2 + lambda + 1 = 0 mod r, lambda^m = 1, chi(2) != 0). Arithmetic only; 0.03 s.

`planner`: the plan of section 3 on real chi-class lists of toy primes at the production degree
parameter GMAX = 16 d^2, sweeping the span; evidence for H4 (f) and the planner bounds:

| d | span/2GMAX | classes | list lengths | min f | mean f | max depth (cap 2^14) |
| ---: | ---: | ---: | --- | ---: | ---: | ---: |
| 16 | 4.0 | 8 | 1993-2119 | 0.9815 | 0.9845 | 71 |
| 16 | 16.0 | 8 | 8095-8328 | 0.9839 | 0.9861 | 110 |
| 16 | 64.0 | 8 | 32565-33065 | 0.9847 | 0.9862 | 125 |
| 64 | 1.0 | 8 | 1962-2091 | 0.9686 | 0.9721 | 90 |
| 64 | 4.0 | 8 | 8144-8220 | 0.9859 | 0.9883 | 128 |
| 64 | 16.0 | 8 | 32636-33002 | 0.9861 | 0.9874 | 269 |

The measured f rises with the span ratio X2/(2 GMAX) (production value 2^31) and with d; every
deficit source (boundary ~ GMAX/span, cycle rejections ~ 1/d, demotions ~ GMAX/span) vanishes at
the production ratios, and the claim (f = 0.98) sits 0.0015 below the smallest measured value
at span/2GMAX >= 4. The depth cap (2^14) never bound in any toy run; the dependency chain length is
geometric with mean 1/(1 - f).

`ec`: Algorithm A end to end on a j = 0 curve y^2 = x^3 + c over F_p with prime order
r = 6 c' 2^5 + 1 and 2 not in H_t: real class lists, the planner, the paired walk with the
shared-inversion pair formula, bases in increasing order with disjoint baby/giant sets, truncated
keys binned into 32 x 64 buckets, insert/probe with false matches counted, candidates retained,
the bounded linear re-walk of both chains, the at most 4x4 cross pairs, the six z tested against
A_d, the d-th root and the root-of-unity scan. Differences from the production procedure (the
toy is not the production at small size): one process, so nothing is routed, posted or polled;
toy tables are built by scalar multiplication; the occupancy, pool and comparison caps are
compared with the model instead of enforced as run failures; the re-walk cap is enforced and
counted. Hosted (the declared 60-trial run, seed 2026): success 0.45 against 0.4017 predicted;
every recovered scalar correct; every sampled collision event agrees with the exact coverage of
the emitted sets, recomputed independently from the plan's emission order; comparisons 9.8
against 9.4 model; occupancy chi-square scaled 0.9842; 4.62 field multiplications per element
against 6 for the reference walk. Offline rows at larger sizes and other seeds:

| toy | bases (baby+giant) | size/side | pair frac | K_eff | predicted | trials | success | comps (model) | mults/el (ref 6) | all scalars correct |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| r=2^25.7 | 2+2 | 383 | 0.8921 | 0.5125 | 0.401 | 400 | 0.395 | 72.7 (71.6) | 4.625 | yes |
| r=2^25.9 | 2+2 | 408 | 0.8879 | 0.512 | 0.4007 | 400 | 0.3675 | 81.4 (81.3) | 4.592 | yes |
| r=2^25.3 | 2+2 | 334 | 0.8975 | 0.5127 | 0.4011 | 400 | 0.3975 | 55.1 (54.5) | 4.653 | yes |
| r=2^22.8 | 1+1 | 139 | 0.8833 | 0.5136 | 0.4017 | 1000 | 0.408 | 10.1 (9.4) | 4.619 | yes |
| r=2^23.0 | 1+1 | 147 | 0.8926 | 0.5105 | 0.3998 | 1000 | 0.394 | 10.9 (10.6) | 4.571 | yes |
| hosted (r=2^22.8) | 1+1 | 139 | 0.8833 | 0.5136 | 0.4017 | 60 | 0.45 | 9.8 (9.4) | 4.619 | yes |

`census`: section 7's exact count at any toy size and K the input names, with the optional
top-window control. Hosted rows (m = 2^21, k = 4) show the method and the small-m deficit; the
offline table above carries the trend. Larger sizes (m above 2^30, or the production
m = 2^234.86) do not fit the hosted budget; H1 remains a heuristic supported by the census.

## 9. W(n), limitations, context

The algorithm reads only A_1 and A_d, so W(n) = W(2^20) = 2^116.485098 for every 2^20 <= n < r;
nothing is claimed for smaller n. Analytic bound, no full-size run; the instruction counts are
itemised, not measured on hardware. The cluster model is mine, not a rule of the challenge.
H1-H4 are disclosed heuristics with their support and extrapolations stated; the census reaches
m = 2^32 at k = 12 and at k = 16, against production 2^234.86; the ln(X2)/d trend is
an empirical ansatz; the planner fraction is measured at span ratios up to 64 against production
2^31; the ec toy is 2^22.8 in size. Memory is exponential, as in any BSGS-form Cheon attack.
Nothing here is a practical break: the attack is 2^116.5 work and 2^122.6 storage against a
255-bit group order.

Negative results from my own search, not claimed: two-for-one P+Q/P-Q inside the H-prime walk
list (three-term progressions of H-primes are too rare at density 2^-20); x-only or projective
steps (no saving against the batched affine reference); three-output triplet formulas
(common-y-denominator: the scalar children leave H, and 14M/3 > 9M/2 per output); implicit
x-bundles as intersection detection (coefficient hashing has counterexamples over F_101 and
roots cost more than the explicit pair); a second routing hop on the frontier's 2^84-owner
layout (send blocks do not fit a 2^40 B node — mine uses 2^56 owners so they do); fewer than 192
key bits (false matches break the 2^44 cap below 2^188).
