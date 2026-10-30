# bls12-381-g1-ecdlp: additive Cheon with planned pair walks, a measured side-output kernel and synchronous wave transport

**Claimed time_log2 = 116.335493** (exact rational 116335493/1000000; an upper
bound on log2 W(2^20), recomputed from exact integers by `evidence.py` mode `ledger`, which also
asserts the production algebra identities below). Heuristic success at least 0.3900 (required
0.39): 1 - e^{-K} > 0.390550 at the realized K = 6N^2/m >= 0.4952, minus the 0.0005 coverage
allowance of H1, a 2^-20 modelling budget, a 4*2^-40 cap budget (the union of the cap tails
below, carried by H3 and H5) and the 2^-111 diagonal event, leaving more than 0.390050, claimed 0.3900.
Peak memory 2^121.5854 bytes with every node counted at its full itemised allocation; every
per-node structure and every per-node total is below 2^63 bytes, so every byte offset is a native
64-bit quantity. No advice. Preprocessing at most 2^115.3389 units, included in W. Every cost,
size and failure figure is rounded up, every success figure rounded down, by the script.

The transport design repairs two defects of the predecessor construction and inherits nothing from it as a premise:

1. *Payload movement is charged at every hop.* There is no router fabric and no declared physical
   axiom. Records move through exactly two charged hops (generating node -> relay -> owner), and
   at each hop every word is written by a counted store of the sender and read by a counted load
   of the receiver (section 4). The transport is a synchronous wave protocol over preallocated
   per-link receive buffers with fixed capacities; an overflow is a run failure, bounded under H3.
2. *There are no credits.* Flow control is replaced by double-buffered waves separated by global
   barriers; every barrier wake, sweep and release is a charged ledger row and the barrier tree is
   explicit (section 4.4). A sequential simulation of the protocol (`evidence.py` mode `proto`,
   section 8) delivers every record exactly once, enforces the buffer discipline, stops on cap
   overflow and reproduces the ledger's closed-form instruction count from counted events.

Two further changes lower the bound. The second output of the pair operation (P - T from the
shared inverse) is no longer itemised as three multiplications plus additions: it is one AArch64
lane-loop function, **1076 instructions per lane including its own lane loop**,
measured in the organizers' Callgrind pipeline
(section 5.1, `measure/`), with blst's Montgomery routines copied verbatim and self-checked
against blst's exported API and against genuine curve points. And each generating node now runs
2^18 lanes over one base (one gap table per node), so the per-lane share of the batched field
inversion drops from 57 to 1 instruction. K = 0.4952; the success floor exceeds 0.390050 with the same H1 allowance, conditional on coverage at this budget.

Cost accounting. Work is the sum of instructions executed over all processors
(tracks/docs/COST_MODEL.md, classical-group-work-v2). Cluster model (section 4.1): nodes with
private memory; the only inter-node mechanism is a store into a link buffer that lives in the
receiver's memory and is mapped write-only into the sender's address space, plus a doorbell
store that raises the receiver's interrupt; every such store is a scalar store instruction of
the sender, every read of the buffer a scalar load of the receiver. The physical premise, as in every instruction-count model here: a store reaches the memory it
addresses and the memory system that carries it executes no counted instruction; every
instruction any processor executes is charged, and explicit
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

Parameters: K = 0.4952; N per side = the least integer with 6N^2 >= K m, N is approximately 2^115.6291; X2 = 2^76
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
   2^18 strands by repeated addition (strand starts by scalar multiplication) with one batched
   inversion per round of 2^18 strand additions; then the node's 2^18 lanes walk the plan in
   lockstep, segments being handed out in plan order to whichever lane finishes its current
   segment (greedy assignment, section 3: every lane finishes within one maximal segment of the
   others, the difference padded and charged), each step of every lane being one pair
   operation or one root operation, with one batched inversion per round of 2^18 lane
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
section 5.1: **1076 instructions per lane**, lane loop included. (Check: lambda' = (-y2 - y1) i = -mu, x' = lambda'^2 - x1 - x2,
y' = lambda'(x1 - x') - y1 = mu (x' - x1) - y1.) Pair arithmetic = 2000 + 1076 = 3076
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
then its subtree depth-first. The planner keeps its state in the planning node's private memory:
parent (32 bits), partner (32), depth (16), height (16), pair-slot (32) and tree size (64) in a
24-byte per-element record, and the accepted child lists as adjacency arrays head (32) and next
(32), 8 further bytes per element; the planning-node memory row (section 6) charges these 32
bytes per element. The four pointer fields are **signed 32-bit offsets from the element's own
absolute class-list index**, not absolute indices: for the element at index i, `parent` stores
i_parent - i, `partner` stores i_partner - i, `head` stores i_first - k for the parent k, and
`next` stores i_next - a for the smaller child a of the previous pair, with -2^31 as the null
sentinel; the smaller child a of a pair decodes its sibling as a + partner[a], and a parent's
child list is a = k + head[k], then a = a + next[a]. Every stored target lies in the same planner
window of WIN = 2(2 GMAX)/d = 2^26 members that section 5.3 scans (the H2 window cap, a
stop-and-fail rule), so |i_parent - i|, |i_partner - i| <= 2^27 < 2^31 and every field fits;
`rel32()` in the executable `plan()` encodes exactly this representation and asserts the bound on
every store. Decoding is one 32-bit load, a sign extension and a 64-bit add (3 instructions) and
encoding is a subtract and a store (2); section 5.3 charges the decodes inside the 16 instructions
per traversal visit and the ancestor walks' 3 per step, and the encodes inside the 64 instructions
per emitted element. `set_depth` and the height updates traverse those child lists (the executable
`plan()` in `evidence.py` does exactly this), and section 5.3 charges that traversal explicitly,
so it is not hidden in the 8/point allowance. The midpoint of a pair is always on the current
path; the ancestor stack holds at most 2D + 1 = 2^15 + 1 points. **Segments and restarts:** the first root at or after each
multiple of SEG = 2^24 emitted entries is emitted as a restart (its value b in the entry); a
segment is the stretch between consecutive restarts, so it has between SEG and SEG + TREEMAX
entries (every tree lies within one segment); a class has at least 2^31 and at most 2^33
segments under the H2 length bounds. **Lane assignment.** While its plan copy arrives, a
generating node records the plan offset of every restart entry in a segment directory (at most
2^33 + 1 entries of 8 bytes, 2^36 bytes). For each base, a node-local cursor hands segments out in
plan order: initially segment j to lane j, and whenever a lane's segment ends, that lane takes the
segment at the cursor (directory load, cursor increment, restart; charged per segment). This is
greedy list scheduling, and its padding bound is deterministic. *Lemma (padding).* Let tau be the
step at which the last segment of a base is handed out. Before tau no lane is idle, so tau is at
most the base's total operations divided by 2^18; at tau every lane is starting or inside a
segment, so it finishes at a step in [tau, tau + SEG + TREEMAX]. Hence every lane finishes within
SEG + TREEMAX steps of every other, the idle (padded) lane steps of a base are at most
2^18 (SEG + TREEMAX), and the longest lane performs at most 2^57 / 2^18 + SEG + TREEMAX <=
(2^57 / (2^24 * 2^18) + 1)(SEG + TREEMAX) = 584133378048 = 2^39.0875 operations per base
(the script's LANE_OPS).
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

There are no intermediate buffers between a store and the receiver's buffer; the relay hops are
charged processor work.

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

A **wave** is 2^34 lockstep steps of a generating node's 2^18 lanes, hence at most 2^53 records
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
No per-link doorbell is rung in the data path: a receiver's wake for a wave is the barrier release of
section 4.4, charged there (one wake of 256 plus a load and dmb per node per barrier, in the
global-barrier row); the doorbell word is used only by the barrier protocol, charged in the same row.
Each count word has the same write/read discipline as its buffer: parity q is written in wave
w = q mod 2, read in wave w+1, and next written in wave w+2 after a barrier; the simulation asserts
this on every buffer. Nothing is retried, nothing waits for space, and no sender observes the
receiver's state; the only failure mode is a cap overflow, counted in the failure budget.

Phases: 2^13 plan-copy rounds (section 3), the baby waves (at most 2^20.09 + 3 per the H2
bounds: 2^15 bases per node, 2^39.0875 steps per base per lane, 2^34 steps per wave, two drain
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

`measure/side_lanes.S` is one function `side_lanes(args, n)` run once per round over the lanes
after the batched inversion; lane i reads its four pointers (out, P, T, inv) from args[i] (the
pointers the pair listing's lane-buffer and argument items prepare). Per lane it loads y1, y2,
forms t = y1 + y2 without reduction (t < 2p), computes mu = t * inv with blst's `__mul_mont_384`,
squares mu with blst's `__sqr_384` + `__mul_by_1_mont_384` (the bodies copied verbatim from
blst's generated `build/elf/mul_mont_384-armv8.S` at the pinned commit f7a88067 and inlined,
only their final `ret` dropped), adds the high half without the final subtraction (V < 2p),
computes x' = V - (x1 + x2): D = V - S lies in (-2p, 2p), 2p is added under the borrow mask and
p is subtracted once conditionally, giving x' in [0, p), and stores it, forms
w = x' + p - x1 without reduction (0 < w < 2p), computes mu * w with `sk_mul_mont_384`,
subtracts y1 modularly, stores y', and computes and stores the key: the low 192 bits of
min(y', p - y'). Inlining frees x18 and x30, which hold the argument cursor and end across
lanes; p, 2p and n0 live in the frame. `measure/gen_side_lanes.py` regenerates the file from a blst checkout
(`reproduce_lanes.sh` checks the shipped file against the generator with `cmp`).

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

**Measurement** (`measure/measurement.log`, `measure/run_measure_lanes.sh`): the organizers' image
and flags, blst built exactly as in `scripts/reference_costs/run.sh`, Callgrind Ir toggled around
one `side_lanes` call over n lanes:

| counted loop | Ir(4096) | Ir(8192) | per iteration |
| --- | ---: | ---: | ---: |
| organizers' `bls_field.c`, rebuilt (binary SHA-256 80a325b8..., as published) | 1568788 | 3137556 | 383 |
| predecessor `side_fused` (v23 package), one C call per lane, re-run | 4616213 | 9232405 | 1127 |
| **`side_lanes`, one call over n lanes (measured)** | 4407354 | 8814650 | **1076** |

Breakdown of the loop body from the disassembly (straight-line, each instruction once per lane,
1076 from the `.Llane` label to its `b.ne`): inlined blst bodies 340 x 2 + 98 + 153 = 931
(each is the predecessor's count minus its `ret`); driver 145 = lane loop 7 (two `ldp` of the
four pointers, cursor `add`, `str` out and `stp` P, T for the reloads inside the
multiplication, `cmp`, `b.ne`), t 12, mu set-up 2, mu spill 3, square set-up 1, p reload 3 +
`mov` 1, V 9, x' 47 (pointer and 6 coordinate loads 7, S 6, D 7, 2p load 3, mask 6, add 6,
conditional subtraction 12), x' stores 3, w 12, y' set-up 2, y' subtraction 23, y' stores 3,
key 17. The prologue/epilogue and constant set-up run once per call and cancel in the
difference. The charged constant is the measured **1076** per lane, lane control included;
the pair listing no longer has a lane-control item (section 5.2). Self-checks after the counted region, all passing at 4096 and 8192 lanes: (1) every
lane against the same formula computed through blst's exported, fully reduced `blst_fp_*`
API, the first 625 lanes covering every combination of {0, 1, p-1, p-2, random} for
x1, y1, x2, y2; (2) Lemma 2a as above; (3) 64 pairs of genuine G1 points [a]G, [b]G with
inv = 1/(x2 - x1) from `blst_fp_inverse`: (x', y') equals blst's own affine P + (-T),
each through a one-lane call (n = 1).

### 5.2 Per-operation listings (listed; charged at the listed count)

**Pair operation, two points** — beyond the reference addition (2000) and the kernel (1076):

| item | instructions |
| --- | ---: |
| plan entry load and decode (1 ldp, one kind test/branch pair, one gap extract, three slot extracts) | 8 |
| load midpoint P from the lane stack (6 ldp, 1 madd with a hoisted 96-byte stride) | 7 |
| load T[g] (6 ldp, 2 address, 3 index, 1 g = 1 select) | 12 |
| batch-inversion lane buffers: 3 stp denominator, 1 madd lane-buffer base, 3 ldp inverse at +48/+64/+80 (adjacent 48-byte halves), 4 kernel argument moves | 11 |
| store P+T (6 stp, 1 address); the kernel stores P-T itself | 7 |
| F_p inversion shared by the round of 2^18 lane operations, ceil(232864 / 2^18) | 1 |
| lane control and loop: inside the measured 1076 (the kernel's own loop) | 0 |
| key of P+T: measured kernel sequence (p - Y, 6-limb compare, 3 csel, 3-word store) | 17 |
| **listed** | **63** |

**Root operation, one point** — beyond the reference addition: gap entry 5, T[gap] load 12
(6 ldp, 2 address, 3 index, 1 g = 1 select; the table base and the 96-byte stride are hoisted,
so the two address instructions are inside the 12 and no separate table-address item is charged),
load/store point 14, batch buffers 6, inversion share 1, loop 5, key 17: listed 60. The pair midpoint address is one `madd` (slot x 96 + hoisted lane-stack base); the T[g] load keeps its 12. The root's only table access is T[gap], whose two address instructions are inside the 12, so the old separate "table address 2" is not charged. The lane-buffer base is one `madd` held across the inversion round; denominator and inverse are adjacent 48-byte halves of one 96-byte per-lane buffer, so the three inverse `ldp` reach +48/+64/+80 from it and the second address is redundant. **Table strand step** — beyond the reference addition: load 8, store 7, batch
buffers 6, inversion share 1, loop 5: listed 27.

**Record transport, per record** (five words, section 4.3):

| hop | instructions | items |
| --- | ---: | --- |
| generating node | 13 | extract t (1), ldp of the window base and fill counter from one 24-byte per-link window-table entry (base, counter, cap) (1), window-table entry address and 40-byte slot address (2), two stp + one str (3), counter increment and store (2), cap compare and branch (2), reload of the kernel-stored key for the side record, ldp + ldr (2), or the main record's loop control (2) |
| relay | 16 | two ldp + one ldr (3), extract u (1), ldp of the owner-window base and fill counter from one 24-byte per-link entry (1), window-table entry address and 40-byte slot address (2), three stores (3), counter update (2), cap check (2), input loop control (2) |
| owner, baby insert | 16 | three loads (3), head index (1), head load (1), pool slot (1), pack the pool record (4), two stp (2), head store (1), pool counter (1), post-indexed record loop control, subs + b.ne (2) |
| owner, giant lookup | 7 | three loads (3), head index (1), head load (1), chain entry and post-indexed record loop control (2) — comparisons are their own row |

The transport and every data path are unchanged; two listings are re-counted at instruction
level. (i) The window base and fill counter of a link are adjacent 64-bit fields
of the same 24-byte window-table entry (base, counter, cap), so they are one `ldp` at the
generating node and at the relay, not two `ldr`; the entry address (t or u times the 24-byte
stride) and the 40-byte record-slot address (one `madd` with the constant stride 40 kept in a
register) remain the two address instructions of the listing, so no address work is dropped. (ii) The owner receive loop reads its
record count from the per-link count word already charged in the per-link-per-wave row, and the
five-word record load is post-indexed (`ldp a,b,[x],#16; ldp c,d,[x],#16; ldr e,[x],#8`), so the loop
control is `subs` + `b.ne` (2) instead of pointer advance + compare + branch (3); the relay's
input loop is unchanged. No word, cap check or counter update is removed: each hop still loads
five words and stores five words, executes the counter update and the cap check. The only other
listing that changes is the candidate/reply row, which reuses REC_RELAY. Residual per-record
bookkeeping not itemised above -- the owner pool-cap compare (two instructions per baby record)
and any head-index work beyond the single extraction already listed -- remains inside the explicit
8-instructions-per-point allowance row, which is charged for both sides; the recount itself only
merges or re-forms instructions already present in the listing.

Listed per record: 45 baby, 36 giant, average 40.5; per link per wave 5 (sender: count store,
counter reset; receiver: count load, loop set-up, bounds). The kernel's key of the side point is the 17-instruction sequence of section 5.1 (6-limb
p - Y, 6-limb compare, 3 csel, 3-word store); the pair listing above charges the same measured
17-instruction sequence for P+T. Both are the identical operation on a 6-limb Montgomery
residue already in registers, so the earlier hand count of 25 is replaced by the measured count;
the loads of Y, if the caller needs them, are part of the unchanged 8-instructions-per-point
allowance.

### 5.3 Setup, planning and distribution listings

Character of a prime: 400000 instructions (366 multiplications at 383 plus the bit-by-bit
Pohlig-Hellman logarithm), charged for every b <= X2 and every prime of I, plus once per prime
below 2^38 per sieve node. Sieve: 16 per integer; class-list transfer: 20 per member (store into
the planning node's window, the planning node's read). Planner: three window scans of
2 x 2 GMAX / d = 2^26 members per element at 16 (the factor 2 is H2 slack), 64 candidate pairs
per element at 2 x 2^14 x 3 + 20 (ancestor walks at the depth cap; the 3 includes the parent
relative-offset decode), emission 64 (including the four relative-offset encodes, 2 sub + 2 str
for the two parent and two partner fields, and the head/next stores), degree counting-sort 16,
and the repeated set_depth/height traversals over the stored child lists: at most 4 l max_tree
subtree visits per class (l/2 merges, each attaching two subtrees of at most max_tree, plus the
root-chain demotion pass with two subtrees per demoted pair) at 16 instructions per visit and
l DEPTH height increments at 8, i.e. 64 max_tree + 8 DEPTH = 2^26.03 instructions per element,
charged as its own ledger row. The 16 per visit include the ldr/sxtw/add decode of the parent
and of each head/next pointer, and the 8 per height step include the parent decode. Child-list
maintenance is inside the matching lump sum: an accepted merge walks at most WIN = 2^26 child
pairs to find its degree/tail and a demotion walks at most WIN to its slot, while the three
window scans charge 3 x 2^26 x 16 = 48 x 2^26 instructions per element, ten times a 5-instruction
per-step walk, so no maintenance walk is uncharged. Plan distribution: per entry per hop one ldp, one stp and two address
instructions (4), over 2^6 + 2^16 + 2^26 + 2^27 copies per class (the copy tree). Base
assignment: 64 per prime of I and 64 per base. Owner initialisation: head array zeroed with
stp (one instruction per 16 bytes) plus 64 per owner. Barriers: section 4.4. Scalar
multiplication (strand starts, restarts, Q_a, candidate tests): 8192000 instructions (2^12
units; 255 doublings and at most 255 additions at the 2000 reference price would be 2040000,
charged four times over).

### 5.4 Per point, both sides

Arithmetic f (2000 + 1076)/2 + (1-f) 2000 = 1547.24 at f = 0.98. Bookkeeping f 63/2 + (1-f) 60
= 32.07 listed, records 40.5 listed; listed 72.57, charged exactly. **Total 1619.81 instructions per point** (exactly 161981/100
in the script); the records row is 2N x 1619.81, and a separate row charges 8 further instructions
per point as an explicit allowance for bookkeeping the listings may have missed (0.49 per cent
of the point).

## 6. Ledger (units of 2000 instructions; row values are rounded, while the script recomputes exact integer totals)

| term | log2 units |
| --- | ---: |
| records: 2N points at 1619.81 instructions | 116.325 |
| listing allowance: 8 instructions per point | 108.664 |
| comparisons at the per-owner cap (1.01 x mean, all 2^60 owners, 7 each) | 106.115 |
| gap tables: BASES x (2^44 - 1) additions at 2000 + 27 | 105.649 |
| padding operations at the pair price (one segment per lane per base, the partial last bases) | 104.367 |
| owner initialisation (2^57 head slots per owner, one stp per 16 bytes) | 105.035 |
| segment restarts: (2N / 2^24 + BASES x 2^18) scalar multiplications, each with its 16-instruction lane hand-off | 104.630 |
| planner (three window scans, matching at the depth cap, emission) | 97.622 |
| planner descendant traversals (adjacency lists, 64 max_tree + 8 DEPTH per element) | 92.038 |
| plan distribution (copy tree, four instructions per entry per hop) | 95.620 |
| per-wave link count words and memory barriers (2^77 + 2^80 links x WAVES x 5, 2 dmb per node per wave) | 92.619 |
| table strand starts and base points: BASES x (2^18 + 1) scalar multiplications | 91.630 |
| candidate re-walks at the 2^44 cap (two segments: 3300 + one inversion per step, tables rebuilt, two scalar multiplications) | 89.024 |
| character sieve, class lists and their transfer | 83.644 |
| global barriers: WAVES x (arrivals, wakes, sweeps, releases over 2^60.002 nodes) | 81.456 |
| base interval sieve, characters of its primes, assignment | 72.646 |
| candidate tests (96 scalar multiplications each at the cap) | 62.588 |
| root recovery and input distribution (192 bytes of inputs to every node) | 55.036 |
| candidate delivery and replies (eight relay-priced hops and 200 per candidate) | 41.392 |
| **W (exact recomputation, see below)** | **< 116.335493** |

The displayed rows are the ceiling to three decimals of every exact integer charge recomputed by `evidence.py` mode `ledger` (so their displayed sum is not the exact total). The exact integer total is 209653892360123246030034579761312949459 instructions, whose log2 in units of 2000 instructions is 116.3354929360572..., and the claimed 116.335493 is its six-decimal ceiling, asserted by the script. Relative to the base 34705e18 (v23, 116.357427) one change applies: the pair kernel is the measured lane-loop kernel at 1076 per lane, lane loop included, replacing v23's 1116 + the 11-instruction lane-control item (1127). The exact total falls by 3211738576180679419676022142607002935 instructions. N, K = 0.4952, f = 0.98, GMAX = 2^44, ALLOW = 8, LANES = 2^18, WSTEP = 2^34, SEG = 2^24, the 17-instruction key and all other rows are inherited unchanged; the planning-node total is 2^62.59 B with the adjacency arrays (still below 2^63) while the total memory stays 2^121.5854 B. BASES = 2 (ceil(N / 2^55) + 1)
< 2^61.630 (the EMIN bound; the sieve computes the exact number, which is about half of this
under H2 — the rows that scale with BASES are charged at the bound). Preprocessing (the baby half of the records, allowance, padding, restarts,
tables, strand starts, link words and barriers; owner initialisation; plan distribution; planner;
sieves) is at most 2^115.3389 units, included. The baby records are charged at their own
transport listing 13+16+16 = 45 instructions per record, not at the two-phase average 40.5
used for the combined total (the giant records are 13+16+7 = 36); the script's baby-phase
term is ceil(N (f PAIR/2 + (1-f) ROOT + 45)), so the asymmetric preprocessing is
charged at its own listing while the combined records row remains exact. Memory and address space, every node at its full
itemised allocation: 2^60 owners at 2^61.59 B (heads, pool, two parities of 2^20 link buffers,
candidate list, reply buffers, code), 2^50 relays at 2^58.32 B (two parities of 2^27 input
buffers, both parities of 2^30 owner windows, window tables, candidate and reply buffers), 2^47
generating nodes at 2^61.70 B (plan array 2^61, segment directory 2^36, both parities of 2^30 relay windows and count words 2^60.3, gap
table 2^50.6, 2^18 lane stacks of 2^15 + 1 points, window entries and fill counters, candidate
buffers, base list),
2^20 planning nodes at 2^62.59 B (class list, hash set of 2^58 slots, 32-byte per-element
planner records with adjacency head/next arrays, copy windows), 2^20 sieve nodes at 2^60.32 B, 2^46 copy nodes at 2^58 B and the barrier
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
K = 6N^2/m >= 0.4952, so more than 0.390550. This is the same heuristic assumption as the accepted frontier's at the changed K and is not implied by acceptance of the predecessor. The following census was conducted at K = 0.5 and provides indirect, not conclusive, support at K = 0.4952. Support:

- *Exact census of the construction.* `evidence.py` mode `census` builds the construction on toy
  primes r = 6 c 2^k + 1 (c odd, smooth, deterministic; skips r with 2 in H_t) at the predeclared
  scaling N = ceil(sqrt(K m/6)), class-size target e = round(N^0.4843), X2 = d e, bases = the
  primes of (X2, X2 + X2^0.921] in order, each walking its full class list, the last base partial
  so each side emits exactly N — the production shape. It computes chi for every b <= X2 by a
  Pohlig-Hellman dlog sieve, forms the emitted sets, and marks the U-coset of every ratio
  (a b)/(a' b') in a byte array of m/6: the count is exact, with every modular and rational
  coincidence included, no sampling. The marking was validated against direct field arithmetic on 18 primes at m <= 2^18.5 (agreement 10^-9). The census rows at K = 0.5:

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
  is roughly proportional to ln(X2)/d_t — an empirical ansatz, not a derivation: the means fall from a few hundredths (k = 4) to a few 10^-4 (k = 12) and within 10^-4 of zero
  (k = 16, four times production's driver). Production's driver is ln(2^76)/2^20 = 5.0e-5 per pair, 60x below k = 12 and 4.4x below k = 16. The top-window control (base near r/X2, products a*b near r, maximal cross-product wrap) matches
the low-window rows, so the trend is window-independent. The E[gcd]-based bound 2*(6/pi^2)(ln X2 + 0.15)/d < 1.1e-4 per pair is consistent with the
ansatz; it bounds literal-equality repeats, and no modular coincidence is excluded anywhere. The repeat-fraction diagnostic is uniform-like for k >= 12: the measured fraction of pairs
  hitting an already-marked coset is within 0.01 of the independent-draw value
  1 - (1 - e^{-K})/K (0.2131 at K = 0.5), while the k = 4 rows sit above it. H1 remains a heuristic: the census reaches m = 2^32 at k = 12 and
  k = 16 (X2 = 2^18.9, against production m = 2^234.86, X2 = 2^76, and products a*b = 2^152 below
  r = 2^254.9 while the toy products sit near or below r depending on the window), so the
  production statement is the explicit 0.0005 allowance. Driver-scaled to production (the ansatz applied to the worst row of every k), no measured worst
  row exceeds -0.0001, five times inside the 0.0005 allowance (the claim tolerates a coverage
  deficit of about 0.00055 before the other terms).
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
any delivery). These t form at most 6 (2^60.64)^2 = 2^123.9 elements of H, probability below 2^-111.

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

tracking H1's ln(X2)/d driver (production is 60x below k = 12); the ec toy's candidate counts
match the model at every size. Under H5's
Poisson-like multiplicity, the 2^20 true-match budget inside the 2^44 cap, the per-node 2^16
candidate slots and the 2^8 per-link candidate caps carry tails below 2^-40 each, and the union is
subtracted as 4*2^-40; the distribution-free fallback is the Markov bound 2^-19.6 above, and even
replacing every 2^-40 by it leaves the floor above 0.390044, still above the claimed 0.3900.

**Claimed success, under H1-H5:** 1 - exp(-0.4952) - 0.0005 - 2^-20 - 4*2^-40 - 2^-111 > 0.390050,
claimed **0.3900** (margin 0.00005 over the required 0.39). The failure budget
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
toy primes at the production degree parameter GMAX = 16 d^2, sweeping the span and, for the
matched cap-active controls below, the tree cap; evidence for H4 (f) and the planner bounds. The
declared run includes the d = 16 uncapped row and the d = 8 and d = 16 capped/uncapped pairs; the
d = 32..128 rows were run offline with the same code and inputs:

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
the production ratios, and the claim (f = 0.98) sits below the smallest uncapped measured value
at d >= 16, span/2GMAX >= 4. These uncapped runs never approach the tree cap; the matched controls next do.

**The tree cap at production: matched capped/uncapped controls.** Production has d = 2^20 and
GMAX = 16 d^2, so a pairing window holds 2 GMAX/d = 32 d = 2^25 class members while TREEMAX =
2^20 = d: the window is exactly 32 caps wide. To exercise this ratio with the cap binding, the
planner is rerun on the *same* real chi-class list at GMAX = 16 d^2 for d = 8, 16, 32, 64, 128
with `max_tree = d` (window/cap = 32, the production ratio) and with the matching uncapped
control `max_tree = 2^20`; seed, list, greedy order, depth cap and emission are identical.
`evidence.py` mode `planner` reports (d = 8, 16 are in the declared run; d = 32..128 offline with
the same code and the same `max_tree` spec field):

| d | list length | tree cap | size rejections | largest tree | f capped | f uncapped | cap excess | B = excess*d |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 8 | 1021-1050 | 8 | 4635 | 7 | 0.8110 | 0.9794 | 0.1684 | 1.35 |
| 16 | 2031-2111 | 16 | 3660 | 15 | 0.8991 | 0.9834 | 0.0843 | 1.35 |
| 32 | 4094-4181 | 32 | 2675 | 31 | 0.9428 | 0.9858 | 0.0430 | 1.38 |
| 64 | 8184-8214 | 64 | 2229 | 63 | 0.9655 | 0.9859 | 0.0204 | 1.31 |
| 128 | 16210-16297 | 128 | 2135 | 127 | 0.9750 | 0.9879 | 0.0129 | 1.65 |

In every capped row the cap is exactly what limits the trees (the largest tree is cap - 1, sizes
being odd) and 2135-4635 candidate merges are rejected; every uncapped control has zero size
rejections. The cap excess approximately halves at each doubling, with B = excess*d in
[1.31, 1.65]. A linear continuation with the largest measured B inflated by 33% to 2.2 bounds
the production cap excess by 2.2/2^20 = 2.1e-6. The same-list uncapped fraction rises with d at
fixed span: 0.9794, 0.9834, 0.9858, 0.9859, 0.9879 at d = 8..128, so the production no-cap
baseline is at least the d = 128 value 0.9879 (production's d = 2^20 and span 2^31 are far beyond
it; the measured deficit terms there are boundary GMAX/span = 2^-32 and cycle 1/d = 2^-20, both
smaller than in the d = 128, span-4 row). Hence f_prod >= 0.9879 - 2.1e-6 > 0.9878, which exceeds
the claimed 0.98 by 7.8e-3; even if the 1/d constant B were 100 times the measured value
(B = 165), the cap excess 1.6e-4 would still leave a 7.7e-3 margin. H4 is therefore the disclosed
extrapolation that the measured 1/d decay of the cap excess continues to d = 2^20 and that the
no-cap baseline does not fall with d; the cap mechanism the earlier evidence did not exercise is
now measured where it binds, and its production contribution is bounded by 2.1e-6. The depth cap
2^14 still never bound (maximum depth 269); it is enforced deterministically, the stack bound
(2D + 1 points) is asserted on every emitted plan, and all cap rejections are counted and charged
in the planner row.

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
buffers, owners insert and look up the wave before that, a global barrier separates waves, and
candidates and replies travel the reverse and forward links. Asserted: every record arrives at the owner of its key
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

`measure/`: `side_lanes.S` (the kernel; blst bodies verbatim, inlined), `side_lanes.S.in` and
`gen_side_lanes.py` (regenerate the file from a blst checkout), `side_lanes_kernel.c` (the
counted call and the three self-checks; the predecessor's `side_fused` files are in the v23
package); `run_measure_lanes.sh` (the in-container build and Callgrind runs,
mirroring `scripts/reference_costs/run.sh`), `reproduce_lanes.sh` (clones the pinned blst, builds the
organizers' image, runs, asserts the published bls_field totals, the self-checks and prints the
per-iteration counts), `measurement.log` (the run used for the claim: SHA-256 of every input,
the Collected totals, the self-check lines and the breakdown).

## 9. W(n), limitations, context

The algorithm reads only A_1 and A_d, so W(n) = W(2^20) <= 2^116.335493 for every 2^20 <= n < r;
nothing is claimed for smaller n.
The cluster model is the predecessor author's (section 4.1); every use is a charged store, load, wake or sweep. H1-H5 are disclosed heuristics with their support and extrapolations stated: the census reaches
m = 2^32 (production 2^234.86), the ln(X2)/d trend is an empirical ansatz, the planner fraction
is measured at span ratios up to 64 (production 2^31), and the ec toy is 2^22.8-2^25.9 in size. The success margin above 0.39 is 0.000050 after H1's 0.0005 allowance and the other failure
terms. Memory is exponential, as in any BSGS-form Cheon attack; this is 2^116.4 work and 2^121.6 storage
against a 255-bit group.

## Revision and provenance

Base 34705e18 (v23, 116.357427), itself on the accepted d421a23d (v22); algorithm, measurements, ledger, H1-H5 and parameters (SEG 2^24, LANES 2^18, WSTEP 2^34, N, K = 0.4952, f = 0.98, GMAX = 2^44) inherited. One change: the pair kernel is the measured lane-loop kernel `side_lanes` (blst bodies inlined, loop inside the kernel, lazy x'), 1076 per lane with lane control included, replacing 1116 + 11. Preprocessing at most 2^115.3389. Per point 1619.81; exact total 209653892360123246030034579761312949459; ceiling 116.335493. No full-size run.
