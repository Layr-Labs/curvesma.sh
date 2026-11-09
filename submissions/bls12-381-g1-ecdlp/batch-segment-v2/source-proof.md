# bls12-381-g1-ecdlp: additive Cheon with planned pair walks, a measured side-output kernel and synchronous wave transport

Superseded source figures (claim 116.3227458, memory, preprocessing, success arithmetic) are omitted; the current claim, ledger, memory and success floor are proof.md's. Every cost, size and failure figure is rounded up, every success figure rounded down, by the script.

Transport: two charged hops (generating node -> relay -> owner), every word stored by the sender and loaded by the receiver, synchronous double-buffered waves between charged global barriers, no credits; an overflow is a run failure (H3). The `proto` simulation (section 8) checks delivery, buffer discipline and the closed-form count.

The pair operation's kernels and prices are proof.md's (s1, s26-s32).

Cost accounting. Work is the sum of instructions executed over all processors (tracks/docs/COST_MODEL.md); the cluster model is section 4.1: every inter-node word is a counted store of the sender and a counted load of the receiver, the memory system executes nothing, an idle node waits (wfi) and executes nothing, and an explicit allowance row (now 1 per point, proof.md s3) covers unlisted bookkeeping. Wall time is not claimed.

## 1. Problem and reduction

BLS12-381 G1, order r, powers A_i = [s^i]P0, i = 1..n. At the scoring n = 2^20 let d = 2^20.
Since v_2(r-1) = 32 (asserted), d | r-1 and 6d | r-1, so m = (r-1)/d = 2^234.857089... is a
multiple of 6 (asserted). H = {x^d} has order m, t := s^d is uniform on H and A_d = [t]P0 is an
input. Only T = A_1 and A_d are read, so the same W applies to every n >= 2^20.

Let lambda = x^2 - 1 mod r with x the curve parameter; then lambda^2 + lambda + 1 = 0 mod r and
lambda^m = 1 (both asserted), so U = {+-1, +-lambda, +-lambda^2} is the sixth-root subgroup and
lies in H. On this curve [lambda](x, y) = (beta x, y) and [-1](x, y) = (x, -y), so the six
points [u]Q, u in U, have y-coordinates +-y(Q). As p is odd, exactly one of y, p-y is even,
and that *even member* is an exact key for the U-orbit of any finite point. The implementation
stores Montgomery residues Y = y*2^384 mod p; on residues the negation map is Y -> p - Y, so the
even one of Y, p - Y is an exact orbit key of the stored value with no normal-form conversion;
records keep its limbs 1..3 (bits 64..255). The kernel computes this key (section 5.1); the toy
uses the parity rule on affine y (even member >> 1); the stored-form rule follows from the
residue identity.

Recovery of s from t is the standard two-step, both charged: (i) twenty Tonelli-Shanks square
roots in F_r give one s0 with s0^d = t; (ii) s = s0*omega^i for one of the d roots of unity
omega = g^((r-1)/d); write i = i0 + 2^10 i1, compute the 2^10 points [s0 omega^{i0}]P0 and the
2^10 points [omega^{-2^10 i1}]T, sort the first by x and look the second up in it, one scalar
multiplication each. Charged at 2^24 units per run in the ledger (the listing: 20 square roots
at 2369 F_r multiplications each, 2^11 scalar multiplications at 2^12 units, 3*2^10 products,
2^11 comparisons of 48-byte x-coordinates at 64 instructions).

## 2. Algorithm

Parameters: X2 = 2^76 (the class integers b); base interval I = (2^76, 2^76 + 2^70]; class lengths between
EMIN = 2^56 - 2^52 (proof.md s35) and EMAX = 2^57 and at most 2^42 bases per character (H2, checked by the sieves). K, N_b, N_g, GMAX, F: proof.md (s2, s30, s31).

1. **Characters and class lists** (challenge independent). chi : F_r^* -> Z/d, chi(x) = x^m
   with chi(x)=log_omega(x^m), omega=g^m of exact order d, has kernel H (and chi(2) != 0: the script asserts
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
   Q_a = [a]P0 by scalar multiplication; the table T[k] = [k]Q_a for 2 <= k <= GMAX is built (proof.md s33); then the node's 2^18 lanes walk the plan in
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

**The pair operation.** For P = (x1, y1) (the midpoint) and T = (x2, y2) = T[g] with x1 != x2, P+T and P-T share the inverse (x2 - x1)^{-1} from Montgomery's trick across the round; kernels and prices: proof.md s1, s26-s28, s30, s32.

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
2 GMAX/d = 2^25 members. **Matching** (order: proof.md s30): each pair accepted only if
it creates no cycle (parent chain walk, at most D = 2^17 steps (current s30)), no depth overflow
(depth(k) + 1 + max(height(i), height(j)) <= D with never-lowered heights, so the check only
over-rejects) and no tree above TREEMAX = 2^20 elements (the size is kept at the tree's root,
found by the same ancestor walk). Only
the final root-chain pass demotes (an element more than GMAX past the last root has its
predecessor demoted), so depths only increase during matching (at most l*D updates) and only
decrease after. Consecutive roots are at most GMAX apart. **Emission:** list order; each root,
then its subtree depth-first. The planner keeps its state in the planning node's private memory:
parent (32 bits), partner (32), depth (32), height (32), pair-slot (32) and tree size (64) in a
28-byte per-element record, and the accepted child lists as adjacency arrays head (32) and next
(32), 8 further bytes per element; the planning-node memory row (section 6) charges these 36
bytes per element before the s16 previous link and s30 centre fields. The four pointer fields are **signed 32-bit offsets from the element's own
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
path; the ancestor stack holds at most 2D + 1 = 2^18 + 1 points. **Segments and restarts:** the first root at or after each
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

Plan entries: proof.md s30. The plan of class c is challenge independent and shared by
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
direct-mapped head array of 2^59 slots (8-byte links, 2^62 bytes) indexed by key bits 69..127 (proof.md s17),
chained through a pool of 2^56 32-byte records (2^61 bytes; more than the mean load, below 2^55.74 at s39,
H3 tail below): a pool record stores the 9-bit key field (bits 60..68), key bits 128..191, the chain id and a full-address link (proof.md s17, s19, s42). Giant records are not stored: each reads its head
slot and compares its 73 remaining key bits with each pooled record on the chain (proof.md s25). An equal key
appends (full 192-bit key, both chain ids, 56 bytes) to the node's candidate list of 2^16 slots.
Comparisons are charged per node at the proof.md s11 cap, all 2^60 owners (row: proof.md s25).

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

Omitted from this copy: the source's 1054-instruction kernel, its Lemma 2a and measurement are superseded by proof.md s1 (929 instructions, measure/, kernel_check.py).

### 5.2 Per-operation listings (superseded; pointers)

The source's pair, root and transport listings are superseded and omitted here: pair and root items by proof.md s26-s28 and s30 (PAIR_BK, ROOT_BK, ROUTE), the key block by s30's x-keys, record transport by s9-s25 and s29 (GEN12, RELAY, INS, LOOK, CMP7, EXT3, generated and executed). The table strand step stands for re-walk rebuilds (tables: proof.md s33): beyond the reference addition, load 8, store 7, batch buffers 6, inversion share, loop 5 (TABLE_BK, with the s1 +21 for v and n). Residual unlisted bookkeeping stays inside the explicit allowance row (proof.md s3).

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

Superseded: per-point prices are proof.md's (s1, s30-s32); the 8-per-point allowance became 1 (proof.md s3).

## 6. Ledger

Omitted from this copy: the source's rounded table, totals, preprocessing and memory figures are superseded by proof.md s4 (exact rows, totals and the per-node memory checks in evidence.py). Retained:
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
evenly over Z/d: class lists have between EMIN = 2^56 - 2^52 and 2^57 elements (proof.md s35), no character has more than 2^42
bases, and the walks need fewer than 2^60 bases per side at the EMIN bound. The
interval I = (2^76, 2^76 + 2^70] holds about 2^64.3 primes by the PNT in intervals
(Rosser-Schoenfeld bounds; the sieving charge covers all 2^70 + 2^43 integers), a margin above 20. Exact class lengths and base counts are computed
by the sieves at runtime; a shortage or an exceeded bound stops the run before any chain starts
and is a failure event under H2, which excludes it. Toy evidence: the class lengths of the ec
toy's 32 classes and the census's class counts at every scale. H2 affects only charged small terms
and the wave/buffer bounds, never the sets.

**H3 (keys).** The 192-bit keys of the 2N records behave as independent uniform strings across
U-orbits: routing digits t, u and head slots are binomial, and a cross-orbit key equality has
probability 2^-192. Support: the toy's occupancy chi-square matches the binomial expectation
(scaled statistic about 1), and its comparison and false-match counts match the model. Tails
under H3, each already multiplied by the number of links, nodes or slots: a generating-node link
above CAP1 = 2^24 in a wave (mean at most 2^23) is below (e/4)^(2^23) over 2^77 links and 2^23
waves; a relay link above CAP2 = 2^21 (mean at most 2^20) below (e/4)^(2^20) over 2^80 links and
2^23 waves; an owner pool above 2^56 (mean below 2^55.74 at proof.md s39; checked per wave, s18) over 2^60
owners; the owner candidate list above 2^16 (mean 2^-20.7), a candidate or reply link above 2^8
(mean below 2^-40) and the node comparison count above 1.01 x mean (Bernstein) are each below
2^-40 over all nodes; false key matches: proof.md s30 (3 N_b N_g/2^192 < 2^40.85 against the 2^44 cap, Chebyshev
below 2^-40); a re-walk chain carrying one key twice is below 2^-81, and one
key at more than 4 steps is below 2^-45 over all chains. The success line subtracts 4*2^-40 for
the union of these.

**Diagonal event.** If t = u a/a' for base primes a, a' of the same character (chi(a) = chi(a')),
then at every common position b the baby and giant points share an orbit, so one base pair gives
more than 2^44 true matches and the run outputs fail at the candidate-count check (which precedes
any delivery). The current bound is 6*2^120/m < 2^-111 (proof.md s39).

**H5 (candidate multiplicity).** For a uniform t in R, the number nu(t) of emitted pairs whose
ratio lies in t's coset — the number of candidates a true hit produces. Two means enter, in the
correct directions: the expected total candidate count is E[#false] + E[nu(t)] = 3 N_b N_g/2^192 +
N_b N_g/M < 2^40.85 + 0.5 (proof.md s30), a conservative upper bound on the exact unconditional mean of the sum, and the cap tails for the total
are evaluated at it; the multiplicity of the true hit itself, given that t is covered, has
conditional mean N_b N_g/(M x coverage) <= 1.28 at the claimed coverage, and the per-run budget of
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
candidate slots and the 2^8 per-link candidate caps carry tails below 2^-40 each, and the union is charged in proof.md section 2. The distribution-free
Markov alternative does not establish the current tightly margined floor;
the present claim explicitly relies on H5.

The current conditional success floor and every failure deduction are in
proof.md section 2. Those exact values supersede this source's older K=0.4952
numbers. All capped work, including aborts, remains in the full ledger.

## 8. Evidence scope

The shipped evidence.py exposes ledger, planner, ec, proto and census modes.
Ledger checks exact accounting and identities, not the production hypotheses.
Planner uses real character lists on toy fields. The ec mode uses actual
small-curve points and verifies every recovered scalar and coverage event;
it has one process, scalar-built tables and model-only occupancy/pool/
comparison caps, while enforcing its re-walk cap. Proto simulates charged
parity-buffer transport, overflow and wrong-generation controls; it is not
a physical cluster. Census counts all ratios at its actual toy parameters.
Measure files are inherited raw assembly and its log; the current derived
listing is generated and checked by kernel_check.py. Section40 of proof.md
specifies the coordinate-only change. All finite evidence has only its stated
scope; the complete code and declared inputs are supplied.

## 9. W(n), limitations, context

The algorithm reads only A_1 and A_d, so W(n) = W(2^20) (bound: proof.md) for every 2^20 <= n < r; nothing is claimed for smaller n.
The cluster model is the predecessor author's (section 4.1); every use is a charged store, load, wake or sweep. H1-H5 are disclosed heuristics with their support and extrapolations stated: the census reaches
m = 2^32 (production 2^234.86), the ln(X2)/d trend is an empirical ansatz, the planner fractions are measured at toy spans (production far larger), and the ec toy is 2^22.8-2^25.9 in size. Memory is exponential, as in any BSGS-form Cheon attack.

