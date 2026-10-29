# bls12-381-g1-ecdlp: additive Cheon with planned pair walks, a measured side-output kernel and synchronous wave transport

**Claimed time_log2 = 116.375** (exact rational 931/8; an upper
bound on log2 W(2^20), recomputed from exact integers by `evidence.py` mode `ledger`, which also
asserts the production algebra identities below). Heuristic success at least 0.3905 (required
0.39): 1 - e^{-K} > 0.391038 at the realized K = 6N^2/m >= 0.496, minus the 0.0005 coverage
allowance of H1, a 2^-20 modelling budget, a 4*2^-40 cap budget (the union of the cap tails
below, carried by H3 and H5) and the 2^-111 diagonal event, leaving more than 0.390537, claimed 0.3905.
Peak memory 2^121.5854 bytes with every node counted at its full itemised allocation; every
per-node structure and every per-node total is below 2^63 bytes, so every byte offset is a native
64-bit quantity. No advice. Preprocessing at most 2^115.3798 units, included in W. Every cost,
size and failure figure is rounded up, every success figure rounded down, by the script.

This is the construction of the predecessor author’s rejected submission `additive-cheon-pairs-v19` (1d94fe4e) with
the two fatal defects of its review repaired and nothing inherited from it as a premise:

1. *Payload movement is charged at every hop.* There is no router fabric and no declared physical
   axiom. Records move through exactly two charged hops (generating node -> relay -> owner), and
   at each hop every word is written by a counted store of the sender and read by a counted load
   of the receiver (section 4). The transport is a synchronous wave protocol over preallocated
   per-link receive buffers with fixed capacities; an overflow is a run failure, bounded under H3.
2. *There are no credits.* Flow control is replaced by double-buffered waves separated by global
   barriers; every barrier wake, sweep and release is a charged ledger row and the barrier tree is
   explicit (section 4.4). A sequential simulation of the protocol (`evidence.py` mode `proto`)
   delivers every record exactly once, checks that every buffer is read exactly in the wave after
   it was written, has every receiver read the record count from the parity-specific count word
   (a single shared count word is shown to be clobbered, as a negative control), stops on cap
   overflow, and reproduces the ledger's closed-form instruction count
   from counted events.

Two further changes lower the bound. The second output of the pair operation (P - T from the
shared inverse) is no longer itemised as three multiplications plus additions: it is one AArch64
function, measured at **1127 instructions per lane** in the organizers' own Callgrind pipeline
(section 5.1, `measure/`), with blst's Montgomery routines copied verbatim and self-checked
against blst's exported API and against genuine curve points. And each generating node now runs
2^16 lanes over one base (one gap table per node), so the per-lane share of the batched field
inversion drops from 57 to 4 instructions. K is 0.496 (the predecessor used 0.5); the success floor exceeds 0.390537 with the same H1 allowance, conditional on a new coverage assumption at this budget.

Cost accounting. Work is the sum of instructions executed over all processors
(tracks/docs/COST_MODEL.md, classical-group-work-v2). Cluster model (section 4.1): nodes with
private memory; the only inter-node mechanism is a store into a link buffer that lives in the
receiver's memory and is mapped write-only into the sender's address space, plus a doorbell
store that raises the receiver's interrupt; every such store is a scalar store instruction of
the sender, every read of the buffer a scalar load of the receiver. The one physical premise is
the one every instruction-count model uses, including the organizers' own kernel measurements:
a store reaches the memory it addresses and the work of the memory system that carries it is
not an instruction; every instruction any processor executes is charged, and an explicit
allowance of 8 instructions per point (a ledger row) covers bookkeeping the listings may have
missed. A node with nothing to do waits for its interrupt and executes no instructions. Wall
time is not claimed.

## 1. Problem and reduction

BLS12-381 G1, order r, powers A_i = [s^i]P0, i = 1..n. At the scoring n = 2^20 let d = 2^20.
Since v_2(r-1) = 32 (asserted), d | r-1 and 6d | r-1, so m = (r-1)/d = 2^234.857089... is a
multiple of 6 (asserted). H = {x^d} has order m, t := s^d is uniform on H and A_d = [t]P0 is an
input. Only T = A_1 and A_d are read, so the same W applies to every n >= 2^20.

Let lambda = x^2 - 1 mod r with x the curve parameter; then lambda^2 + lambda + 1 = 0 mod r and
lambda^m = 1 (both asserted), so U = {+-1, +-lambda, +-lambda^2} is the sixth-root subgroup and
lies in H. On this curve [lambda](x, y) = (beta x, y) and [-1](x, y) = (x, -y), so the six
points [u]Q, u in U, have y-coordinates +-y(Q), and min(y, p-y) is an exact key for the U-orbit
of any finite point. The implementation stores Montgomery residues Y = y*2^384 mod p; on
residues the negation map is Y -> p - Y, so min(Y, p - Y) and its low 192 bits are an exact
orbit key of the stored value with no normal-form conversion. The measured kernel computes
exactly this key for the side output (section 5.1); the toy exercises the min rule on affine y;
the stored-form rule follows from the residue identity, which is the only claim made for it.

Recovery of s from t is the standard two-step, both charged: (i) twenty Tonelli-Shanks square
roots in F_r give one s0 with s0^d = t; (ii) s = s0*omega^i for one of the d roots of unity
omega = g^((r-1)/d); write i = i0 + 2^10 i1, compute the 2^10 points [s0 omega^{i0}]P0 and the
2^10 points [omega^{-2^10 i1}]T, sort the first by x and look the second up in it, one scalar
multiplication each. Charged at 2^24 units per run in the ledger (the listing: 20 square roots
at 2369 F_r multiplications each, 2^11 scalar multiplications at 2^12 units, 3*2^10 products,
2^11 comparisons of 48-byte x-coordinates at 64 instructions).

## 2. Algorithm

Parameters: K = 0.496; N per side = the least integer with 6N^2 >= K m, N is approximately 2^115.6303; X2 = 2^76
(the class integers b); GMAX = 2^44 = 16 d^2 (the pair-gap and table bound); base interval
I = (2^76, 2^76 + 2^70]; f = 0.98, the claimed minimum pair fraction (H4); class lengths between
EMIN = 2^55 and EMAX = 2^57 and at most 2^42 bases per character (H2, checked by the sieves).

1. **Characters and class lists** (challenge independent). chi : F_r^* -> Z/d, chi(x) = x^m
   projected by a fixed generator of mu_d, has kernel H (and chi(2) != 0: the script asserts
   2^m != 1). A segmented sieve over b <= X2 (2^20 sieve nodes, 2^56 integers each) factors b
   over the primes below 2^38 (their characters computed once per node, 400000 instructions
   each: one 235-bit exponentiation and a bit-by-bit Pohlig-Hellman logarithm in mu_{2^20}) and
   computes the character of the at most one remaining prime factor above 2^38 by the same
   400000-instruction routine; the ledger charges that routine for every b. The class list
   L_c = {b in [1, X2] : chi(b) = c}, in increasing order, has about X2/d = 2^56 elements (H2);
   the exact lengths are computed by the sieve. Each sieve node streams its members of class c,
   in order, into a per-link buffer of the planning node of class c (2^20 links per planning
   node, 2^40 bytes each), so the planning node holds its sorted class list in those buffers.
2. **Planning** (challenge independent, section 3). For each class c a plan produces every
   element of L_c exactly once: each element is a root (one addition from the previous root), one
   child of a pair (x-g, x, x+g) with midpoint x and gap g <= GMAX, or a *restart* (a root computed
   by scalar multiplication, carrying its 77-bit value b). The plan is emitted as a sequence of
   16-byte entries and distributed to the 2^27 generating nodes of class -c through a copy tree.
3. **Base primes**. Sieve I and compute chi(a) for its primes in increasing order (one
   400000-instruction character per prime). The baby side takes them in order until its N elements
   are emitted; the giant side continues with the next primes. A base a of character chi(a) walks
   the plan of class -chi(a), so every emitted scalar a*b lies in H. Baby and giant base sets are
   disjoint, and a > X2 >= b, so a is the unique prime factor of a*b above X2: (a, b) is recovered
   from the product by one factorisation. Bases of character c are assigned in order, round-robin,
   to the 2^27 generating nodes of class -c (at most 2^15 per node under the 2^42 bound), as a
   list of at most 2^15 (base, index) pairs per node, transferred once.
4. **Baby walk.** A generating node processes its bases one after another. For each base a:
   Q_a = [a]P0 by scalar multiplication; the table T[k] = [k]Q_a for 2 <= k <= GMAX is built in
   2^16 strands by repeated addition (strand starts by scalar multiplication) with one batched
   inversion per round of 2^16 strand additions; then the node's 2^16 lanes walk the plan in
   lockstep, segments being handed out in plan order to whichever lane finishes its current
   segment (greedy assignment, section 3: every lane finishes within one maximal segment of the
   others, the difference padded and charged), each step of every lane being one pair
   operation or one root operation, with one batched inversion per round of 2^16 lane
   operations. Every produced point emits a 40-byte record (192-bit truncated key, 128-bit chain
   id = base index, segment index, spare) into the transport (section 4). The
   whole baby side completes, and every baby record is inserted at its owner, before any giant
   chain starts (a global barrier).
5. **Giant walk.** The same with Q_{a'} = [a']A_d. Giant records are routed to the owner of their
   key and compared against the stored baby records there; they are not stored.
6. **Collision.** A truncated-key match is a candidate (key, baby id, giant id), appended at the
   owner node (2^16 slots; a 2^16+1st fails the node). After the giant side, the barrier tree
   sums the candidate counts; if the global count exceeds 2^44 the run outputs fail. Otherwise
   each candidate is delivered to the two generating nodes that own its records (section 4.5),
   each re-walks the one plan segment named by the record's chain id once (section 6, charged
   at the cap), keeps the at most 4 positions whose key equals the retained key (a fifth is a
   failure), and returns them with the base to the owner, which
   forms the at most 16 pairs, tests the at most 96 values z = u*b/g by [z]P0 = A_d, and on a hit
   recovers s by section 1. A wrong candidate never yields a wrong answer: every scalar is
   verified against A_d and then T before output.

**The pair operation.** For P = (x1, y1) (the midpoint) and T = (x2, y2) = T[g] with x1 != x2,
both P+T and P-T need the single inverse i = (x2 - x1)^{-1}, provided by Montgomery's trick
shared across the round. P+T is the reference batched incomplete affine addition (6M = 2000
instructions: 3M formula + 3M batch share, with the reference's own field additions inside the
A = 0 normalization). P-T = (x', y') with mu = (y1 + y2) i, x' = mu^2 - x1 - x2,
y' = mu (x' - x1) - y1 is computed, together with its orbit key, by the measured kernel of
section 5.1: **1127 instructions**. (Check: lambda' = (-y2 - y1) i = -mu, x' = lambda'^2 - x1 - x2,
y' = lambda'(x1 - x') - y1 = mu (x' - x1) - y1.) Pair arithmetic = 2000 + 1127 = 3127
instructions for two points; nothing of the second formula is inside the reference's A = 0
convention, and nothing of it is charged twice: the kernel is priced at its measured count, the
reference addition at its fixed price.

**Incomplete additions are safe.** Every point is a multiple of Q = Q_a != O. For a pair with
midpoint x and gap g: the children are list elements, so 0 < g < x and x + g <= X2 < r; neither
x = g nor x = -g mod r. For a root step u -> v = u + gap with gap <= GMAX: u = -gap would force
v = 0; u = gap means v = 2u, but chi(2) != 0, so chi(v) != chi(u) and v would not be in the
class. The table's k = 2 step is the one doubling. The giant side is the same with
Q = [a't]P0 != O. The toy asserts a non-zero denominator at every step.

## 3. The planner

Fix a class list L = (b_0 < ... < b_{l-1}), l ~ 2^56. (i, j), i < j, is a pair with midpoint k if
b_k - b_i = b_j - b_k = g, 1 <= g <= GMAX. **Edges:** for each u in L, the class members v in
(u, u + 2 GMAX] with u = v mod 2 whose midpoint (u+v)/2 is in L (a hash-set probe); under H2 an
element lies in about 2 GMAX/d^2 = 32 such pairs as a child, and the window holds about
2 GMAX/d = 2^25 members. **Matching:** greedy by increasing degree, each pair accepted only if
it creates no cycle (parent chain walk, at most D = 2^14 steps), no depth overflow
(depth(k) + 1 + max(height(i), height(j)) <= D with never-lowered heights, so the check only
over-rejects) and no tree above TREEMAX = 2^20 elements (the size is kept at the tree's root,
found by the same ancestor walk). No augmenting paths and no unmatching during matching; only
the final root-chain pass demotes (an element more than GMAX past the last root has its
predecessor demoted), so depths only increase during matching (at most l*D updates) and only
decrease after. Consecutive roots are at most GMAX apart. **Emission:** list order; each root,
then its subtree depth-first (children re-enumerated by a window scan, so no adjacency is
stored), so the midpoint of a pair is always on the current path; the ancestor stack holds at
most 2D + 1 = 2^15 + 1 points. **Segments and restarts:** the first root at or after each
multiple of SEG = 2^25 emitted entries is emitted as a restart (its value b in the entry); a
segment is the stretch between consecutive restarts, so it has between SEG and SEG + TREEMAX
entries (every tree lies within one segment); a class has at least 2^29 and at most 2^32
segments under the H2 length bounds. **Lane assignment.** While its plan copy arrives, a
generating node records the plan offset of every restart entry in a segment directory (at most
2^32 + 1 entries of 8 bytes, 2^35 bytes). For each base, a node-local cursor hands segments out in
plan order: initially segment j to lane j, and whenever a lane's segment ends, that lane takes the
segment at the cursor (directory load, cursor increment, restart; charged per segment). This is
greedy list scheduling, and its padding bound is deterministic. *Lemma (padding).* Let tau be the
step at which the last segment of a base is handed out. Before tau no lane is idle, so tau is at
most the base's total operations divided by 2^16; at tau every lane is starting or inside a
segment, so it finishes at a step in [tau, tau + SEG + TREEMAX]. Hence every lane finishes within
SEG + TREEMAX steps of every other, the idle (padded) lane steps of a base are at most
2^16 (SEG + TREEMAX), and the longest lane performs at most 2^57 / 2^16 + SEG + TREEMAX <=
(2^57 / 2^25 / 2^16 + 1)(SEG + TREEMAX) = 2^41.047 operations per base (the script's LANE_OPS).
The bound uses only the per-segment maximum SEG + TREEMAX, which holds by construction. The
segment index in a record's chain id names the segment whatever lane walked it. The last base of a side walks only the
prefix the side still needs (its idle lanes are padded and charged).

A plan entry is 16 bytes: kind (2 bits), gap (45 bits, covering g = GMAX), source and two
destination slots (16 bits each, covering the 2^15 + 1 stack); a restart entry carries b
(77 bits) instead of gap and slots. The plan of class c is challenge independent and shared by
every base of character -c on both sides. It is not stored at the planning node: as it is
emitted, in 2^13 rounds of 2^48 bytes, it is copied through a tree (planning node -> 2^6 ->
2^10 -> 2^10 copy nodes -> 2 generating nodes each, 2^27 generating nodes per class) into the
local plan array of every generating node of the class (2^61 bytes per node), each hop being a
charged copy (section 5.3). The planner cost (three window scans per element at 16
instructions per member, matching with bounded ancestor checks, emission) is a ledger row at
2^97.6 units, nineteen orders below the records row, and every access names its data path
(section 5). The pair fraction f_c is known before any group operation; the claim uses f = 0.98
for every class (H4).

## 4. Matching and transport

### 4.1 Cluster model

Node classes, all counted in the memory total: 2^47 generating nodes (2^27 per class), 2^50
relays R(t, s) indexed by the 30-bit owner digit t and a 20-bit sender digit s, 2^60 owners
(t, u) indexed by the first 60 key bits, 2^20 planning nodes, 2^20 sieve nodes, 2^46 copy nodes
and 2^50 barrier nodes. Every node has private memory and a 64-bit address space holding its
own structures and its *windows*: a window is a region of another node's memory (a link buffer)
mapped write-only into this node's address space at a fixed virtual address, established once
at initialisation (the window tables are written once, charged in the setup rows). A
(sender, receiver) pair with a window is a *link*; links are fixed by the design (listed below)
and never change. The sender moves data by ordinary store instructions into the window; the
memory system delivers them to the receiver's buffer; no processor executes anything for the
delivery. The receiver reads its buffer with ordinary loads. Each link also has, in the
receiver's memory, two *count words*, one per buffer parity (the count of parity q is written
by the sender only at the end of a wave in which it wrote the parity-q buffer, and is read by the
receiver only in the next wave, together with that buffer), and the
receiver has one *doorbell* word whose store raises the receiver's interrupt (a message-signalled
interrupt: a store to a doorbell address is translated by the receiver's interrupt controller
into a local interrupt — GICv3 ITS/LPI; the doorbell carries no data). There is no remote load,
no remote atomic, no shared read-modify-write, no credit, no polling: every word of every buffer
has exactly one writer (its link's sender) and exactly one reader (the receiver), and a node
without work waits for its interrupt (wfi) executing nothing. Memory ordering uses dmb barriers,
charged where they occur.

This is the model's whole inter-node vocabulary. In v19 the same stores were routed through a
seven-hop switch fabric whose buffers were not charged; here there are no intermediate buffers
of any kind between a store and the receiver's buffer, and the two relay hops of the record path
are charged processor work (loads, stores, addressing, cap checks) at the relay node.

### 4.2 Owners and matching

The 60-bit owner key is split (t, u) = (30, 30). Within an owner, baby records are inserted into a
direct-mapped head array of 2^57 slots (8-byte links, 2^60 bytes) indexed by key bits 60..116,
chained through a pool of 2^56 32-byte records (2^61 bytes; more than 1.29 times the mean load below 2^55.636,
H3 tail below): a pool record stores key bits 117..191 (75 bits), the chain id (62-bit base
index and 32-bit segment index, 94 bits) and a 56-bit link: 225 of 256 bits. Giant records are not stored: each reads its head
slot and compares its 75 remaining key bits with each pooled record on the chain (7 instructions
per comparison: two loads, compare, conditional compare, branch, link load, loop). An equal key
appends (full 192-bit key, both chain ids, 56 bytes) to the node's candidate list of 2^16 slots.
Comparisons are charged per node up to 1.01 times the mean, all 2^60 owners (row 2^106.13).

### 4.3 Record path: two charged hops

A record is 40 bytes (five words): the 192-bit truncated key and the 128-bit chain id. The
generating node g has 2^30 links, one to each relay R(t, s_g) with s_g = g mod 2^20; the relay
R(t, s) has 2^27 incoming links (the generating nodes with that s) and 2^30 outgoing links, one
to each owner (t, u); the owner (t, u) has 2^20 incoming links (the relays R(t, .)). Every
receiver's fan-in and every sender's fan-out is fixed and the per-link buffers are preallocated
at their caps for both parities of the double buffering:

| link | per-link buffer (records) | mean load per wave (H3) | buffers at the receiver | windows at the sender |
| --- | ---: | ---: | ---: | ---: |
| generating node -> relay | CAP1 = 2^24 | <= 2^23 | 2^27 x 2 x 2^24 x 40 B = 2^57.3 B | 2^30 x 2 x 2^24 x 40 B = 2^60.3 B |
| relay -> owner | CAP2 = 2^21 | <= 2^20 | 2^20 x 2 x 2^21 x 40 B = 2^47.3 B | 2^30 x 2 x 2^21 x 40 B = 2^57.3 B |

(The sender maps both parities of every link buffer; the window columns are address space at the
sender, counted in its per-node total.)

A **wave** is 2^36 lockstep steps of a generating node's 2^16 lanes, hence at most 2^53 records
per generating node per wave (two per pair step), at most 2^23 per link in the mean (t uniform
under H3), and at most 2^50 records per relay per wave, 2^20 per owner link in the mean. Wave w
of the whole machine, between two global barriers: every generating node produces its wave-w
records and stores each into the relay link buffer of parity w mod 2 (the record's t selects the
window, the link's fill counter the slot; a counter reaching CAP1 fails the run); every relay
reads its parity (w-1) mod 2 input buffers (the wave w-1 records) and stores each into the owner
link buffer of parity (w-1) mod 2 (u selects the window; CAP2 likewise); every owner reads its
parity w mod 2 input buffers (the wave w-2 records) and inserts or looks them up. At the end of
its wave work a sender executes dmb, stores every link's parity-(w mod 2) count word, dmb, and
arrives at the barrier; after its release wake a receiver executes dmb and, in wave w+1, reads the
parity-(w mod 2) count words, never the other parity, which the senders are writing in that wave.
Each count word therefore has the same write/read discipline as its buffer. Hence a buffer
of parity q is written in wave w (w = q mod 2), read in wave w+1, and next written in wave w+2,
after a barrier that follows the completion of the read; the simulation asserts this discipline
on every buffer. Nothing is retried, nothing waits for space, and no sender ever observes the
receiver's state: the only failure mode is a cap overflow, counted in the failure budget.

Phases: 2^13 plan-copy rounds (section 3), the baby waves (at most 2^20.09 + 3 per the H2
bounds: 2^15 bases per node, 2^41.047 steps per base per lane, 2^36 steps per wave, two drain
waves and one rounding wave), a phase barrier, the giant waves, the candidate rounds (4.5). The
ledger bounds the number of barriers by 2 (2^20.09 + 3) + 2^13 + 64 = 2170982 (WAVES in the
script).

### 4.4 Global barrier

A barrier tree of fan-in 2^10 spans all nodes (2^60.0017 leaves; 2^50.0 internal barrier nodes
over six levels, computed by the script's `barrier_tree`). An internal node holds one arrival
word per child (written by that child through its window) and each node holds a release word
(written by its parent). Arrival: dmb, store the barrier generation into the parent's arrival
word, dmb, store the parent's doorbell (4 instructions). The parent's interrupt handler (wake
charged 256 instructions) sweeps its 2^10 arrival words (2 instructions per slot: load, compare)
and, when all carry the current generation, arrives at its own parent; the worst case — one wake
and one full sweep per child arrival, no coalescing — is what the ledger charges: per barrier,
every non-root node arrives once (4), every arrival wakes its parent once (256) and triggers one
full sweep (2 per child slot), and every non-root node is released by one wake (256) plus one
load and dmb (2). Release: the root stores the generation into each child's release word and its
doorbell, recursively. Generation numbers make clearing unnecessary; a stale doorbell is harmless
(the handler re-checks the words). The barrier arrival word also carries the node's candidate
count during the status pass, so the status reduction is the same mechanism. The toy
simulation computes its barrier cost by this exact rule and asserts equality with the ledger's
formula at its own parameters.

### 4.5 Candidate rounds and replies

After the giant phase and the status pass (candidate total <= 2^44, else fail), one *candidate
round* delivers each candidate twice along the reverse direction of the record path — the owner
(t, u) stores the 56-byte candidate into its link buffer at relay R(t, s_b) where s_b is the baby
chain's generating node mod 2^20 (reverse links: 2^30 per relay, cap 2^8 candidates each), and
likewise at R(t, s_g) for the giant chain; after a barrier each relay forwards to the generating
node named by the chain id (reverse links: 2^30 per generating node, cap 2^8 each); after a
barrier the generating nodes re-walk the named segment (section 6). The *reply round* returns,
along the forward links (generating node -> R(t_key, s_g) -> owner (t_key, u_key)), the base a
and the at most four kept values b of each record's segment (one 64-byte reply per record, caps
2^8 per link); the owner then has both sides and forms the at most 16 products (a b)/(a' b').
Every node received P0 and A_d (192 bytes) with the first barrier release (charged in the
input-distribution row), so owners can run the scalar tests.
Per-link loads are tiny in the mean (2^39.3 false candidates over 2^80 owner-relay pairs), and
any link above its cap fails the run (H3 tail inside the 4*2^-40 budget). The simulation runs
both rounds and asserts that every candidate is resolved exactly once.

## 5. Instruction listings (AArch64 scalar, 64-bit registers)

### 5.1 The measured side-output kernel

`measure/side_fused.S` is one function `side_fused(out, P, T, inv)`: it loads y1, y2, forms
t = y1 + y2 without reduction (t < 2p), computes mu = t * inv with blst's `__mul_mont_384`
(copied verbatim from blst's generated `build/elf/mul_mont_384-armv8.S` at the pinned commit
f7a88067 and renamed `sk_mul_mont_384`), squares mu with blst's `__sqr_384` +
`__mul_by_1_mont_384` + the inlined accumulation and conditional subtraction of
`__redc_tail_mont_384` (blst's own `sqr_mont_384` path, result stores dropped), computes
x' = (mu^2 - x1) - x2 by two full modular subtractions and stores it, forms
w = x' + p - x1 without reduction (0 < w < 2p), computes mu * w with `sk_mul_mont_384`,
subtracts y1 modularly, stores y', and computes and stores the key: the low 192 bits of
min(y', p - y'). `measure/gen_side_fused.py` regenerates the file from a blst checkout
(`reproduce.sh` checks the shipped file against the generator with `cmp`).

**Lemma 2a (lazy operands).** blst's `__mul_mont_384` is word-serial Montgomery multiplication
with one final conditional subtraction: for inputs a, b with a < 4p and b < p it returns
exactly a b 2^{-384} mod p, fully reduced. Proof: after each of the six word steps the
accumulator is below a + p < 5p < 2^384 (the invariant acc_{j+1} = (acc_j + b_j a + q_j p)/2^64
< (acc_j + 2^64 a + 2^64 p)/2^64), so the six-limb accumulator with its one carry limb never
overflows; the loop computes S = (a b + q p)/2^384 with q < 2^384, so S < (a b)/2^384 + p
< 4p^2/2^384 + p < 2p because 4p < 2^383 < 2^384; S is then compared with p and p subtracted
once if S >= p, giving S mod p in [0, p). Both lazy operands here satisfy a < 2p. The harness checks the lemma directly
(600000 products with a = r + k p, k = 1..3, against the reduced product; including
r = b = p - 1): 0 mismatches.

**Measurement** (`measure/measurement.log`, `measure/run_measure.sh`): the organizers' image
and flags, blst built exactly as in `scripts/reference_costs/run.sh`, Callgrind Ir toggled around
a loop of one call per lane:

| counted loop | Ir(4096) | Ir(8192) | per iteration |
| --- | ---: | ---: | ---: |
| organizers' `bls_field.c`, rebuilt (binary SHA-256 80a325b8..., as published) | 1568788 | 3137556 | 383 |
| **`side_fused`, one call per lane (charged)** | 4616213 | 9232405 | **1127** |

Breakdown from the disassembly (straight-line code, each instruction executed once per call):
driver 181 (prologue/epilogue 17, constant loads 6, frame stores 2, t 12, two call set-ups 4,
mu spill 3, square set-up and modulus reloads 7, reduction tail 23, two modular subtractions
for x' 45, w 12, modular subtraction for y' 23, nine stores 6, key 17, four `bl` 4),
`sk_mul_mont_384` 341 x 2, `sk_sqr_384` 99, `sk_mul_by_1_mont_384` 154, C lane loop 11:
total 1127. Self-checks after the counted region, all passing at 4096 and 8192 lanes: (1) every
lane against the same formula computed through blst's exported, fully reduced `blst_fp_*`
API, the first 625 lanes covering every combination of {0, 1, p-1, p-2, random} for
x1, y1, x2, y2; (2) Lemma 2a as above; (3) 64 pairs of genuine G1 points [a]G, [b]G with
inv = 1/(x2 - x1) from `blst_fp_inverse`: (x', y') equals blst's own affine P + (-T). The
charged price is the measured 1127 per lane; it includes the key extraction of the side point,
the stores of both coordinates and the C loop.

### 5.2 Per-operation listings (listed; charged x 21/20)

**Pair operation, two points** — beyond the reference addition (2000) and the kernel (1127):

| item | instructions |
| --- | ---: |
| plan entry load and decode (ldp, kind, gap, three slots) | 8 |
| load midpoint P from the lane stack (6 ldp, 2 address) | 8 |
| load T[g] (6 ldp, 2 address, 3 index) | 12 |
| batch-inversion lane buffers: denominator store, inverse load and keep (for the kernel) | 12 |
| store P+T (6 stp, 1 address); the kernel stores P-T itself | 7 |
| F_p inversion shared by the round of 2^16 lane operations, ceil(232864 / 2^16) | 4 |
| lane control and loop | 10 |
| key of P+T: p - Y (6), compare (6), select (3), pack (10) | 25 |
| **listed** | **86** |

**Root operation, one point** — beyond the reference addition: gap entry 5, table address 2,
T[gap] load 12, load/store point 14, batch buffers 6, inversion share 4, loop 5, key 25:
listed 73. **Table strand step** — beyond the reference addition: load 8, store 7, batch
buffers 6, inversion share 4, loop 5: listed 30.

**Record transport, per record** (five words, section 4.3):

| hop | instructions | items |
| --- | ---: | --- |
| generating node | 14 | extract t (1), window base load (1), fill counter load (1), address (2), two stp + one str (3), counter increment and store (2), cap compare and branch (2), reload of the kernel-stored key for the side record or loop control for the main record (2) |
| relay | 17 | two ldp + one ldr (3), extract u (1), window base load (1), counter load (1), address (2), three stores (3), counter update (2), cap check (2), input loop control (2) |
| owner, baby insert | 17 | three loads (3), head index (1), head load (1), pool slot (1), pack the pool record (4), two stp (2), head store (1), pool counter (1), loop (3) |
| owner, giant lookup | 8 | three loads (3), head index (1), head load (1), chain entry and loop (3) — comparisons are their own row |

Listed per record: 48 baby, 39 giant, average 43.5; per link per wave 5 (sender: count store,
counter reset; receiver: count load, loop set-up, bounds). The kernel's key of the side point
replaces the 25-instruction key item for that point; the pair listing above charges one
25-instruction key (for P+T) per pair operation.

### 5.3 Setup, planning and distribution listings

Character of a prime: 400000 instructions (366 multiplications at 383 plus the bit-by-bit
Pohlig-Hellman logarithm), charged for every b <= X2 and every prime of I, plus once per prime
below 2^38 per sieve node. Sieve: 16 per integer; class-list transfer: 20 per member (store into
the planning node's window, the planning node's read). Planner: three window scans of
2 x 2 GMAX / d = 2^26 members per element at 16 (the factor 2 is H2 slack), 64 candidate pairs
per element at 2 x 2^14 x 3 + 20 (ancestor walks at the depth cap), emission 64, degree
counting-sort 16. Plan distribution: per entry per hop one ldp, one stp and two address
instructions (4), over 2^6 + 2^16 + 2^26 + 2^27 copies per class (the copy tree). Base
assignment: 64 per prime of I and 64 per base. Owner initialisation: head array zeroed with
stp (one instruction per 16 bytes) plus 64 per owner. Barriers: section 4.4. Scalar
multiplication (strand starts, restarts, Q_a, candidate tests): 8192000 instructions (2^12
units; 255 doublings and at most 255 additions at the 2000 reference price would be 2040000,
charged four times over).

### 5.4 Per point, both sides

Arithmetic f (2000 + 1127)/2 + (1-f) 2000 = 1572.23 at f = 0.98. Bookkeeping f 86/2 + (1-f) 73
= 43.60 listed, records 43.5 listed; listed 87.10, charged x 21/20 = 91.455. **Total 1663.685
instructions per point** (exactly 332737/200 in the script); the records row is 2N x 1663.685,
and a separate row charges 8 further instructions per point as an explicit allowance for
bookkeeping the listings may have missed (0.5 per cent of the point).

## 6. Ledger (units of 2000 instructions; row values are rounded, while the script recomputes exact integer totals)

| term | log2 units |
| --- | ---: |
| records: 2N points at 1663.685 instructions | < 116.362 |
| listing allowance: 8 instructions per point | < 108.665 |
| comparisons at the per-owner cap (1.01 x mean, all 2^60 owners, 7 each) | 106.128 |
| gap tables: BASES x (2^44 - 1) additions at 2000 + 32 | 105.659 |
| owner initialisation (2^57 head slots per owner, one stp per 16 bytes) | 105.034 |
| segment restarts: (2N / 2^25 + BASES x 2^16) scalar multiplications, each with its 16-instruction lane hand-off | 103.636 |
| padding operations at the pair price (one segment per lane per base, the partial last bases) | 103.366 |
| planner (three window scans, matching at the depth cap, emission) | 97.622 |
| plan distribution (copy tree, four instructions per entry per hop) | 95.620 |
| per-wave link count words and memory barriers (2^77 + 2^80 links x WAVES x 5, 2 dmb per node per wave) | 92.576 |
| table strand starts and base points: BASES x (2^16 + 1) scalar multiplications | 89.636 |
| candidate re-walks at the 2^44 cap (two segments: 3300 + one inversion per step, tables rebuilt, two scalar multiplications) | 89.023 |
| character sieve, class lists and their transfer | 83.644 |
| global barriers: WAVES x (arrivals, wakes, sweeps, releases over 2^60.002 nodes) | 81.413 |
| base interval sieve, characters of its primes, assignment | 72.646 |
| candidate tests (96 scalar multiplications each at the cap) | 62.587 |
| root recovery and input distribution (192 bytes of inputs to every node) | 55.036 |
| candidate delivery and replies (eight relay-priced hops and 200 per candidate) | 41.427 |
| **W (conservative bound, see below)** | **< 116.375** |

The displayed rows are rounded presentations of the predecessor’s exact integer charges at K = 0.5, except the record and listing rows, which use the new N. The monotonic comparison below uses the exact integer charges, not their rounded displays. The predecessor total is 216274957927844535493409341692259263733 instructions. Keeping every other predecessor row at its original bound and replacing only these two rows yields at most 215409911693892974638256554088505215771 instructions, whose log2 in units of 2000 instructions is approximately 116.37456796 < 116.375. All other original rows remain valid upper bounds because they are constant or monotone nondecreasing in N; comparisons, base count, restarts and padding can only decrease. The ledger program recomputes the tighter costs and prints their own six-decimal ceiling; 116.375 is the deliberately looser claimed bound. BASES = 2 (ceil(N / 2^55) + 1)
< 2^61.636 (the EMIN bound; the sieve computes the exact number, which is about half of this
under H2 — the rows that scale with BASES are charged at the bound). For context, the same
budgets with every point computed as a root give 2^116.73; the pair operation is worth 0.35
bits, the measured kernel 0.039 bits against the 3149 + 72 itemisation of v19, and the 2^16-lane
inversion share 0.023 bits. Preprocessing (the baby half of the records, allowance, padding, restarts,
tables, strand starts, link words and barriers; owner initialisation; plan distribution; planner;
sieves) is at most 2^115.3798 units, included. Memory and address space, every node at its full
itemised allocation: 2^60 owners at 2^61.59 B (heads, pool, two parities of 2^20 link buffers,
candidate list, reply buffers, code), 2^50 relays at 2^58.32 B (two parities of 2^27 input
buffers, both parities of 2^30 owner windows, window tables, candidate and reply buffers), 2^47
generating nodes at 2^61.70 B (plan array 2^61, segment directory 2^35, both parities of 2^30 relay windows and count words 2^60.3, gap
table 2^50.6, 2^16 lane stacks of 2^15 + 1 points, window entries and fill counters, candidate
buffers, base list),
2^20 planning nodes at 2^62.33 B (class list, hash set of 2^58 slots, parent/depth/height/size
arrays, copy windows), 2^20 sieve nodes at 2^60.32 B, 2^46 copy nodes at 2^58 B and the barrier
nodes: 2^121.5854 bytes in all, with the per-node totals asserted below 2^63 by the script.
Charge accounting: the key extraction of P+T is inside the pair bookkeeping and that of P-T
inside the kernel (nowhere else); the gap tables, their strand starts and inversions are only
the two table rows; the walk inversions are the per-lane share inside the records row; the
re-walk inversions are inside the re-walk row. Each instruction is charged exactly once across
these categories.

## 7. Success probability

Steps 1-3 of section 2 are challenge independent; the run fails and outputs no scalar if the base
searches come up short or a class or base count exceeds its bound (H2), or the plan breaks the
gap bound (never observed; the walk checks its gaps and stops before any chain starts). The
algorithm succeeds exactly when t lies in R = U * {b/g} over the emitted sets B, G (|B| = |G| = N)
and no cap breaks. Because t is uniform on H and R is a fixed set (the lists do not depend on the
challenge), Pr[t in R] = |R|/m exactly.

**H1 (coverage).** |R|/m at the production construction is at least 1 - e^{-K} - 0.0005, with
K = 6N^2/m >= 0.496, so more than 0.390538. This is a new heuristic assumption at the changed K and is not implied by acceptance of the predecessor. The following census was conducted at K = 0.5 and provides indirect, not conclusive, support at K = 0.496. Support:

- *Exact census of the construction.* `evidence.py` mode `census` builds the construction on toy
  primes r = 6 c 2^k + 1 (c odd, smooth, deterministic; skips r with 2 in H_t) at the predeclared
  scaling N = ceil(sqrt(K m/6)), class-size target e = round(N^0.4843), X2 = d e, bases = the
  primes of (X2, X2 + X2^0.921] in order, each walking its full class list, the last base partial
  so each side emits exactly N — the production shape. It computes chi for every b <= X2 by a
  Pohlig-Hellman dlog sieve, forms the emitted sets, and marks the U-coset of every ratio
  (a b)/(a' b') in a byte array of m/6: the count is exact, with every modular and rational
  coincidence included, no sampling. The marking was validated against direct field arithmetic
  (an independent brute-force reconstruction of R) on 18 primes at m up to 2^18.5, agreeing to
  10^-9 (v19's validation; the census code is unchanged). The census rows at K = 0.5:

| toy m | d_t = 2^k | primes | mean diff | lowest | repeat frac (uniform value) |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 2^24 | 2^4 | 8 | -0.028426 | -0.058765 | 0.270 (0.213) |
| 2^24 | 2^8 | 8 | -0.002240 | -0.007191 | 0.218 (0.213) |
| 2^24 | 2^12 | 8 | -0.000212 | -0.000861 | 0.214 (0.213) |
| 2^24 | 2^16 | 8 | -0.000025 | -0.000283 | 0.213 (0.213) |
| 2^26 | 2^4 | 8 | -0.029179 | -0.043471 | 0.271 (0.213) |
| 2^26 | 2^8 | 8 | -0.003658 | -0.008197 | 0.220 (0.213) |
| 2^26 | 2^12 | 8 | -0.000150 | -0.000338 | 0.213 (0.213) |
| 2^26 | 2^16 | 8 | +0.000046 | -0.000024 | 0.213 (0.213) |
| 2^28 | 2^4 | 8 | -0.026602 | -0.049898 | 0.266 (0.213) |
| 2^28 | 2^8 | 8 | -0.003035 | -0.013144 | 0.219 (0.213) |
| 2^28 | 2^12 | 8 | -0.000257 | -0.000607 | 0.214 (0.213) |
| 2^28 | 2^16 | 8 | -0.000027 | -0.000113 | 0.213 (0.213) |
| 2^30 | 2^4 | 8 | -0.027199 | -0.044774 | 0.267 (0.213) |
| 2^30 | 2^8 | 8 | -0.008180 | -0.037002 | 0.229 (0.213) |
| 2^30 | 2^12 | 8 | -0.000191 | -0.000329 | 0.213 (0.213) |
| 2^30 | 2^16 | 8 | -0.000006 | -0.000045 | 0.213 (0.213) |
| 2^32 | 2^4 | 8 | -0.025740 | -0.031054 | 0.265 (0.213) |
| 2^32 | 2^8 | 8 | -0.007391 | -0.036289 | 0.228 (0.213) |
| 2^32 | 2^12 | 8 | -0.000233 | -0.000409 | 0.214 (0.213) |
| 2^32 | 2^16 | 8 | -0.000013 | -0.000036 | 0.213 (0.213) |

  Top-window control (products a*b near r, maximal modular wraparound), K = 0.5:

| toy m | d_t | primes | mean diff | lowest |
| ---: | ---: | ---: | ---: | ---: |
| 2^26 | 2^8 | 4 | -0.002780 | -0.004537 |
| 2^26 | 2^12 | 4 | -0.000278 | -0.000398 |
| 2^28 | 2^8 | 4 | -0.001469 | -0.001728 |
| 2^28 | 2^12 | 4 | -0.000228 | -0.000370 |
| 2^30 | 2^12 | 4 | -0.000163 | -0.000222 |

  The deficit shrinks with the toy degree. Empirically, across the measured rows the mean deficit
  is roughly proportional to ln(X2)/d_t — an empirical ansatz, not a derivation: the k = 4 rows
  (driver about 0.35-0.4) measure means of a few hundredths; the k = 8 rows (driver about 0.04)
  means of a few thousandths with worst rows a few hundredths at small d (the small-d rows drift
  with m; every k >= 12 row is flat in m); the k = 12 rows (driver about 0.003) have means of a
  few 10^-4, and the k = 16 rows (driver about 0.0002, four times production's) have means within
  10^-4 of zero. Production has d = 2^20 with driver ln(2^76)/2^20 = 5.0e-5 per pair, 60 times
  below the k = 12 driver and 4.4 times below the k = 16 driver. The top-window control reruns the
  census with the base window at (r/X2, r/X2 + (r/X2)^0.921], so that the products a*b are near r
  and the ratio cross-products wrap modulo r as much as they can; its rows sit at the same deficits
  as the low-window rows at the same (m, k), so the measured trend does not depend on the window
  height. The analytic bound on the production rational-repeat rate (E[gcd]-based,
  2*(6/pi^2)(ln X2 + 0.15)/d < 1.1e-4 per pair) is consistent with the ansatz; it bounds the
  literal-equality repeats, which unique factorisation does constrain, and no modular coincidence
  is excluded anywhere. The repeat-fraction diagnostic shows the k >= 12 maps are uniform-like:
  the measured fraction of pairs hitting an already-marked coset is within 0.01 of the value
  1 - (1 - e^{-K})/K that independent uniform draws would give (0.2131 at K = 0.5), while the
  k = 4 rows sit well above it. H1 remains a heuristic: the census reaches m = 2^32 at k = 12 and
  k = 16 (X2 = 2^18.9, against production m = 2^234.86, X2 = 2^76, and products a*b = 2^152 below
  r = 2^254.9 while the toy products sit near or below r depending on the window), so the
  production statement is the explicit 0.0005 allowance. Driver-scaled to production (the ansatz
  applied to the worst row of every k), no measured worst row exceeds -0.0001, five times inside
  the allowance; the new 0.3905 claim tolerates a coverage deficit of about 0.000537 before the other
  terms, and the required 0.39 tolerates about 0.001037; the heuristic allowance is 0.0005.
- *End-to-end agreement.* The ec toy of section 8 checks, every trial, whether the sampled t lies
  in the exactly-enumerated R of its emitted sets, and every observed candidate/miss event agrees.

**H2 (equidistribution).** The characters of the primes of I and of the integers b <= X2 spread
evenly over Z/d: class lists have between 2^55 and 2^57 elements, no character has more than 2^42
bases, and the walks need about 2^59.64 bases per side (charged at the 2^61.64 EMIN bound). The
interval I = (2^76, 2^76 + 2^70] holds about 2^64.3 primes by the PNT in intervals
(Rosser-Schoenfeld bounds; the sieving charge covers all 2^70 + 2^43 integers), a 12 x margin
over the 2^60.64 bases per side at the bound. Exact class lengths and base counts are computed
by the sieves at runtime; a shortage or an exceeded bound stops the run before any chain starts
and is a failure event under H2, which excludes it. Toy evidence: the class lengths of the ec
toy's 32 classes and the census's class counts at every scale. H2 affects only charged small terms
and the wave/buffer bounds, never the sets.

**H3 (keys).** The 192-bit keys of the 2N records behave as independent uniform strings across
U-orbits: routing digits t, u and head slots are binomial, and a cross-orbit key equality has
probability 2^-192. Support: the toy's occupancy chi-square matches the binomial expectation
(scaled statistic about 1), and its comparison and false-match counts match the model. Tails
under H3, each already multiplied by the number of links, nodes or slots: a generating-node link
above CAP1 = 2^24 in a wave (mean at most 2^23) is below (e/4)^(2^23) over 2^77 links and 2^22
waves; a relay link above CAP2 = 2^21 (mean at most 2^20) below (e/4)^(2^20) over 2^80 links and
2^22 waves; an owner pool above 2^56 (mean below 2^55.636) below exp(-2^55.6 x 0.027) over 2^60
owners; the owner candidate list above 2^16 (mean 2^-20.7), a candidate or reply link above 2^8
(mean below 2^-40) and the node comparison count above 1.01 x mean (Bernstein) are each below
2^-40 over all nodes; false key matches have mean N^2/2^192 < 2^39.27 against the 2^44 cap,
Chebyshev below 2^-47 under H3; a re-walk chain carrying one key twice is below 2^-81, and one
key at more than 4 steps is below 2^-45 over all chains. The success line subtracts 4*2^-40 for
the union of these.

**Diagonal event.** If t = u a/a' for base primes a, a' of the same character (chi(a) = chi(a')),
then at every common position b the baby and giant points share an orbit, so one base pair gives
more than 2^44 true matches and the run outputs fail at the candidate-count check (which precedes
any delivery). These t form at most 6 (2^60.64)^2 = 2^123.9 elements of H, so the event has
probability below 2^-111, far inside the 2^-20 budget.

**H5 (candidate multiplicity).** For a uniform t in R, the number nu(t) of emitted pairs whose
ratio lies in t's coset — the number of candidates a true hit produces. Two means enter, in the
correct directions: the expected total candidate count is E[#false] + E[nu(t)] = N^2/2^192 +
N^2/M < 2^39.27 + 0.5, a conservative upper bound on the exact unconditional mean of the sum, and the cap tails for the total
are evaluated at it; the multiplicity of the true hit itself, given that t is covered, has
conditional mean N^2/(M x coverage) <= 1.28 at the claimed coverage, and the per-run budget of
2^20 true matches is Markov-bounded by it (Pr[nu > 2^20 | covered] <= 1.28/2^20 < 2^-19.6).
Support: the census's multiplicity diagnostic counts nu exactly (no sampling) and reports the
second-moment ratio sum(nu^2)/(T + T^2/M) (the value 1.000 for independent uniform pairs), K = 0.5:

| toy m | d_t | primes | energy ratio (1.000 = independent) |
| ---: | ---: | ---: | ---: |
| 2^22 | 2^4 | 4 | 1.2089 (worst 1.2832) |
| 2^24 | 2^8 | 4 | 1.0235 (worst 1.0425) |
| 2^24 | 2^12 | 4 | 1.0010 (worst 1.0026) |
| 2^26 | 2^8 | 4 | 1.0225 (worst 1.0385) |
| 2^26 | 2^12 | 4 | 1.0009 (worst 1.0015) |
| 2^28 | 2^12 | 4 | 1.0021 (worst 1.0034) |

tracking the same ln(X2)/d driver as H1 (production's driver is 60 times below k = 12's); the
ec toy's candidate counts match the false-plus-true model at every size run. Under H5's
Poisson-like multiplicity, the 2^20 true-match budget inside the 2^44 cap, the per-node 2^16
candidate slots and the 2^8 per-link candidate caps carry tails below 2^-40 each, and the union is
subtracted as 4*2^-40; the distribution-free fallback is the Markov bound 2^-19.6 above, and even
replacing every 2^-40 by it leaves the floor above 0.390532, still above the claimed 0.3905.

**Claimed success, under H1-H5:** 1 - exp(-0.496) - 0.0005 - 2^-20 - 4*2^-40 - 2^-111 > 0.390537,
claimed **0.3905** (margin 0.0005 over the required 0.39). The failure budget
covers: base-search shortage and exceeded H2 bounds, no collision (1 - e^{-K} and the H1
allowance), link-cap overflows, pool and candidate-list overflows, comparison caps (H3), false
truncated-key matches and the global candidate cap, diagonal events, re-walk key repeats and the
kept-value cap. Every cap, abort and failure path is charged at its maximum in the ledger.

## 8. Evidence (`evidence.py`, standard library; `measure/`, the kernel)

`ledger`: updated costs from exact integers, their six-decimal ceiling with the boundary
assert, the preprocessing and memory inequalities with every per-node total (memory plus
mapped windows) below 2^63, the
success inequality with every failure term, and the production algebra identities (v_2(r-1) = 32,
6 | m, lambda = x^2 - 1 with lambda^2 + lambda + 1 = 0 mod r, lambda^m = 1, chi(2) != 0).
Arithmetic only; 0.03 s.

`planner`: the plan of section 3 (with the depth and tree-size caps) on real chi-class lists of
toy primes at the production degree parameter GMAX = 16 d^2, sweeping the span; evidence for H4
(f) and the planner bounds:

| d | span/2GMAX | classes | list lengths | min f | mean f | max depth (cap 2^14) | max tree size (cap 2^20) | size rejections |
| ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: |
| 16 | 4.0 | 8 | 1993-2119 | 0.9815 | 0.9845 | 71 | 1751 | 0 |
| 16 | 16.0 | 8 | 8095-8328 | 0.9839 | 0.9861 | 110 | 5491 | 0 |
| 16 | 64.0 | 8 | 32565-33065 | 0.9847 | 0.9862 | 125 | 4893 | 0 |
| 64 | 1.0 | 8 | 1962-2091 | 0.9686 | 0.9721 | 90 | 1895 | 0 |
| 64 | 4.0 | 8 | 8144-8220 | 0.9859 | 0.9883 | 128 | 5757 | 0 |
| 64 | 16.0 | 8 | 32636-33002 | 0.9861 | 0.9874 | 269 | 17147 | 0 |

The measured f rises with the span ratio X2/(2 GMAX) (production value 2^31) and with d; every
deficit source (boundary ~ GMAX/span, cycle rejections ~ 1/d, demotions ~ GMAX/span) vanishes at
the production ratios, and the claim (f = 0.98) sits below the smallest measured value at
span/2GMAX >= 4. Neither cap ever bound in any toy run (maximum depth seen in the hundreds
against 2^14, maximum tree size in the thousands against 2^20, zero size rejections), and both
are enforced deterministically by the planner: no probabilistic claim about depth or tree size is
made anywhere, the stack bound (2D + 1 points) and the tree bound are asserted on every emitted
plan, and cap rejections are counted and charged in the planner row. As a rough sanity check (not
a premise): if chain lengths were geometric with mean 1/(1 - f) = 50, the maximum over l = 2^56
elements would be about 50 ln(2^56) = 1950, far below either cap; the claim f >= 0.98 nevertheless
remains a disclosed extrapolation from the measured sizes.

`ec`: Algorithm A end to end on a j = 0 curve y^2 = x^3 + c over F_p with prime order
r = 6 c' 2^5 + 1 and 2 not in H_t: real class lists, the planner, the paired walk with the
shared-inversion pair formula, bases in increasing order with disjoint baby/giant sets, truncated
keys binned into 32 x 64 buckets, insert/probe with false matches counted, candidates retained,
the bounded linear re-walk of both chains, the at most 4x4 cross pairs, the six z tested against
A_d, the d-th root and the root-of-unity scan. Differences from the production procedure (the
toy is not the production at small size): one process, so nothing is routed or buffered; toy
tables are built by scalar multiplication; the occupancy, pool and comparison caps are compared
with the model instead of enforced as run failures; the re-walk cap is enforced and counted.
Hosted (the declared 60-trial run, seed 2026, K = 0.5): success 0.45 against 0.3972 predicted;
every recovered scalar correct; every sampled collision event agrees with the exact coverage of
the emitted sets, recomputed independently from the plan's emission order; comparisons 9.7
against 9.3 model; occupancy chi-square scaled 0.9859; 4.59 field multiplications per element
against 6 for the reference walk. Offline rows at larger sizes and other seeds:

| toy | bases (baby+giant) | size/side | pair frac | K_eff | predicted | trials | success | comps (model) | mults/el (ref 6) | all scalars correct |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| r=2^25.7 | 2+2 | 379 | 0.8921 | 0.5018 | 0.3946 | 400 | 0.39 | 71.2 (70.1) | 4.619 | yes |
| r=2^25.9 | 2+2 | 404 | 0.8879 | 0.502 | 0.3947 | 400 | 0.355 | 79.8 (79.7) | 4.593 | yes |
| r=2^25.3 | 2+2 | 330 | 0.8975 | 0.5005 | 0.3938 | 400 | 0.3925 | 53.8 (53.2) | 4.655 | yes |
| r=2^22.8 | 1+1 | 138 | 0.8833 | 0.5062 | 0.3972 | 1000 | 0.401 | 10.0 (9.3) | 4.587 | yes |
| r=2^23.0 | 1+1 | 146 | 0.8926 | 0.5036 | 0.3957 | 1000 | 0.393 | 10.8 (10.4) | 4.572 | yes |

`proto`: the section-4 transport at toy scale (8 generating nodes of 4 lanes, 8 relays, 16 owners,
3 waves per phase plus two drain waves, fan-in-4 barrier tree), sequentially: records with random
192-bit keys are stored into per-link buffers with caps, relays forward the previous wave's
buffers, owners insert and look up the wave before that, a global barrier (arrivals, worst-case
wakes and full sweeps, releases) separates waves, candidates are delivered along the reverse links
and replies along the forward links. Asserted: every record arrives at the owner of its key
exactly once; every buffer is read exactly in the wave after it was written and never reused
before; every candidate is resolved once; the simulated instruction total (the section-5 listing
applied to counted events: records per hop, inserts, lookups, comparisons, link count words,
barrier events) equals the closed form the ledger evaluates at production size (`barrier_instr`,
`LINK_WAVE`, `REC_*` are shared between the two). Three patterns: uniform keys (2275 records, 10
barriers, closed form equal), keys skewed to two owners (wide caps, closed form equal), all keys
to one owner (the cap overflows at the first relay link and the run fails, as section 4.3
specifies). Every receiver obtains the number of records from the count word of the buffer's
parity and the simulation checks that this word describes exactly the buffer being read; a
negative control with one count word for both parities is detected as clobbered. A lane-schedule
check runs the greedy segment assignment of section 3 on 20 random and adversarial length
sequences (lengths at both extremes, sorted either way, alternating) and verifies the padding
lemma in each: finishing steps within one maximal segment, padding at most LANES (SEG + TREEMAX)
(the stride assignment, for comparison, exceeds that bound on the adversarial sequences). The simulation models the event ordering the protocol constrains, not hardware;
instruction counts are the listing applied to counted events; evidence at simulation scope.

`census`: section 7's exact count at any toy size and K the input names, with the optional
top-window control and the optional exact multiplicity diagnostic of H5. Hosted row (m = 2^21,
k = 4, K = 0.5) shows the method and the small-m deficit; the offline tables above carry the
trend. Larger sizes (m above 2^32, or the production m = 2^234.86) do not fit the hosted budget;
H1 remains a heuristic supported by the census.

`measure/`: `side_fused.S` (the kernel; blst routines verbatim, renamed), `side_driver.S.in` and
`gen_side_fused.py` (regenerate the file from a blst checkout), `side_kernel.c` (measurement
loop and the three self-checks), `run_measure.sh` (the in-container build and Callgrind runs,
mirroring `scripts/reference_costs/run.sh`), `reproduce.sh` (clones the pinned blst, builds the
organizers' image, runs, asserts the published bls_field totals, the self-checks and prints the
per-iteration counts), `measurement.log` (the run used for the claim: SHA-256 of every input,
the Collected totals, the self-check lines and the breakdown).

## 9. W(n), limitations, context

The algorithm reads only A_1 and A_d, so W(n) = W(2^20) <= 2^116.375 for every 2^20 <= n < r;
nothing is claimed for smaller n. Analytic bound, no full-size run; the kernel is measured in the
organizers' pipeline, the remaining instruction counts are itemised, not measured on hardware.
The cluster model is the predecessor author’s, not a rule of the challenge: its one inter-node mechanism (a store
into a write-mapped link buffer, a doorbell interrupt) is stated in section 4.1 and every use of
it is a charged store, load, wake or sweep, and the listing allowance is explicit. H1-H5 are disclosed heuristics with their support and
extrapolations stated; the census reaches m = 2^32 at k = 12 and at k = 16, against production
2^234.86; the ln(X2)/d trend is an empirical ansatz; the planner fraction is measured at span
ratios up to 64 against production 2^31; the ec toy is 2^22.8-2^25.9 in size. The success margin
above 0.39 exceeds 0.000537 at the floor and rests on H1's 0.0005 allowance (five times the driver-scaled
worst census row). Memory is exponential, as in any BSGS-form Cheon attack. Nothing here is a
practical break: the attack is 2^116.4 work and 2^121.6 storage against a 255-bit group order.

Negative results from the predecessor author’s search, not claimed: two-for-one P+Q/P-Q inside the H-prime walk
list (three-term progressions of H-primes are too rare at density 2^-20); x-only or projective
steps (no saving against the batched affine reference); three-output triplet formulas
(common-y-denominator: the scalar children leave H, and 14M/3 > 9M/2 per output); implicit
x-bundles as intersection detection (coefficient hashing has counterexamples over F_101 and
roots cost more than the explicit pair); fewer than 192 key bits (false matches break the 2^44
cap below 2^188); a single routing hop (the per-link count words over 2^107 generating-node/owner
pairs cost more than the relay hop they would save); a one-wave design (the generating node's
windows would exceed 2^63 bytes).

## Revision and provenance

This revision inherits the algorithm, measurements, evidence and ledger from submission 2c102f66-1237-441b-a9fd-a04cc60e53b7. It changes only the production N and declared K to 0.496, updates success and the conservative work cap, and retains the original toy experiments at their actual K = 0.5. The coverage heuristic at the new K is a new assumption; no previous verdict is a premise. No full-size execution or new algorithm is claimed.
