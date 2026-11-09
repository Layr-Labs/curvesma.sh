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
most the base's total operations divided by LANES = 2^19 (proof.md s42); at tau every lane is starting or inside a
segment, so it finishes at a step in [tau, tau + SEG + TREEMAX]. Hence every lane finishes within
SEG + TREEMAX steps of every other, the idle (padded) lane steps of a base are at most
LANES (SEG + TREEMAX), and the longest lane performs at most 2^57 / LANES + SEG + TREEMAX <=
(2^57 / (SEG LANES) + 1)(SEG + TREEMAX) = 343602626560 operations per base at SEG = 2^22
(the script's LANE_OPS).
The bound uses only the per-segment maximum SEG + TREEMAX, which holds by construction. The
segment index in a record's chain id names the segment whatever lane walked it. The last base of a side walks only the
prefix the side still needs (its idle lanes are padded and charged).

Plan entries: proof.md s30. The plan of class c is challenge independent and shared by
every base of character -c on both sides. It is not stored at the planning node: as it is
emitted, in 2^13 rounds of 2^48 bytes, it is copied through a tree (planning node -> 2^6 ->
2^10 -> 2^10 copy nodes -> 2 generating nodes each, 2^27 generating nodes per class) into the
local plan array of every generating node of the class (2^61 bytes per node), each hop being a
charged copy (section 5.3). The planner cost (window scans, matching with bounded ancestor checks,
emission; section 5.3) is in ledger rows below 0.01% of W (proof.md s4), and every access names its data path
(section 5). The pair fraction f_c is known before any group operation; the claim uses the proof.md s30
floors (H4).

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
chained through a pool of 2^57 32-byte records (2^62 bytes; more than the mean load, below 2^56.41 at proof.md s43,
H3 tail below): a pool record stores key bits 117..191 (75 bits), the chain id (62-bit base
index and 32-bit segment index, 94 bits) and a 56-bit link: 225 of 256 bits. Giant records are not stored: each reads its head
slot and compares its 75 remaining key bits with each pooled record on the chain (7 instructions
per comparison: two loads, compare, conditional compare, branch, link load, loop). An equal key
appends (full 192-bit key, both chain ids, 56 bytes) to the node's candidate list of 2^16 slots.
Comparisons are charged per node up to 1.01 times the mean, all 2^60 owners (row: proof.md s25).

### 4.3 Record path: two charged hops

A record is 40 bytes (five words): the 192-bit truncated key and the 128-bit chain id. The
generating node g has 2^30 links, one to each relay R(t, s_g) with s_g = g mod 2^20; the relay
R(t, s) has 2^27 incoming links (the generating nodes with that s) and 2^30 outgoing links, one
to each owner (t, u); the owner (t, u) has 2^20 incoming links (the relays R(t, .)). Every
receiver's fan-in and every sender's fan-out is fixed and the per-link buffers are preallocated
at their caps for both parities of the double buffering:

CAP1 = 2^24 records on each generating node -> relay link and CAP2 = 2^21 on each relay -> owner link (mean loads per wave at most 2^23 and 2^20 under H3); the sender maps both parities of every link buffer as windows, counted in its per-node total (evidence.py memory terms).

A **wave** is 2^33 lockstep steps of a generating node's 2^19 lanes (proof.md s42; giant waves 2^32 steps, s30), hence at most 2^53 records
per generating node per wave (two per pair step, six per giant step), at most 2^23 per link in the mean (t uniform
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

Phases: 2^13 plan-copy rounds (section 3), the baby waves (at most ceil(2^15 LANE_OPS / 2^33) + 2 under the H2
bounds: 2^15 bases per node, LANE_OPS steps per base per lane, 2^33 steps per wave, drain
and rounding waves; giant waves of 2^32 steps), a phase barrier, the giant waves, the candidate rounds (4.5). The
ledger bounds the number of barriers by WAVES = 6561960 (both phases, plan rounds, 64 spare).

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
Per-link loads are tiny in the mean (2^42.19 false candidates over 2^80 owner-relay pairs), and
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
the planning node's window, the planning node's read). Planner: the proof.md s30 rows in evidence.py (six window scans of WIN = 2^29 members per element at 16, 5 x 64 kappa candidate tries per element at 2 DEPTH x 3 + 20 (ancestor walks at the depth cap), emission 64 with the four relative-offset encodes, degree counting-sort 16; descendant traversals at most 4 l max_tree visits at 16 and l DEPTH height steps at 8 per class, the 16 and 8 including the pointer decodes; s16 augmentation rows). Child-list maintenance is inside the matching lump sum: an accepted merge or a demotion walks at most WIN child pairs, while the six scans charge 96 x 2^29 instructions per element, above a 5-instruction per-step walk, so no maintenance walk is uncharged. Plan distribution: per entry per hop one ldp, one stp and two address
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

**H1 (coverage).** |R|/m at the production construction is at least 1 - e^{-K} - 0.0005 at the realized K = 6 N_b N_g/m >= 251/200 (proof.md s2, s43), so above 0.7144. The census below was run at K = 1.255 and supports, not proves, this at production. Support:

- *Exact census of the construction.* `evidence.py` mode `census` builds the construction on toy
  primes r = 6 c 2^k + 1 (c odd, smooth, deterministic; skips r with 2 in H_t) at the predeclared
  scaling N_g least with 6 ceil(Q N_g) N_g >= K m (Q = 29/25), class-size target e = round(N_b^0.4811) (production log(X2/d)/log N_b), X2 = d e, bases = the
  primes of (X2, X2 + X2^0.921] in order, each walking its full class list, the last base partial
  so each side emits exactly N — the production shape. It computes chi for every b <= X2 by a
  Pohlig-Hellman dlog sieve, forms the emitted sets, and marks the U-coset of every ratio
  (a b)/(a' b') in a byte array of m/6: the count is exact, with every modular and rational
  coincidence included, no sampling. The marking was validated against direct field arithmetic on 18 primes at m <= 2^18.5 (agreement 10^-9). The census rows at K = 1.255 (112 runs; 'top': products a*b near r, maximal wraparound; last column: lowest diff times production's ln(X2)/d over the row's):

| toy m | d_t | primes | mean diff | lowest | scaled lowest |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 2^24 | 2^4 | 4 | -0.028762 | -0.045615 | -5.7e-06 |
| 2^24 | 2^8 | 8 | -0.002882 | -0.008123 | -1.1e-05 |
| 2^24 | 2^12 | 8 | -0.000209 | -0.000554 | -9.5e-06 |
| 2^24 | 2^16 | 8 | -0.000078 | -0.000511 | -1.1e-04 |
| 2^26 | 2^4 | 4 | -0.037817 | -0.053888 | -6.4e-06 |
| 2^26 | 2^8 | 8 | -0.004157 | -0.008438 | -1.1e-05 |
| 2^26 | 2^12 | 8 | -0.000168 | -0.000339 | -5.7e-06 |
| 2^26 | 2^16 | 8 | -0.000075 | -0.000265 | -5.8e-05 |
| 2^28 | 2^8 | 8 | -0.003737 | -0.016661 | -2.2e-05 |
| 2^28 | 2^12 | 8 | -0.000212 | -0.000473 | -7.7e-06 |
| 2^28 | 2^16 | 8 | -0.000009 | -0.000139 | -3.0e-05 |
| 2^30 | 2^8 | 4 | -0.015753 | -0.045937 | -5.8e-05 |
| 2^30 | 2^12 | 4 | -0.000186 | -0.000203 | -3.2e-06 |
| 2^30 | 2^16 | 4 | -0.000037 | -0.000061 | -1.3e-05 |
| 2^32 | 2^12 | 4 | -0.000292 | -0.000472 | -7.3e-06 |
| 2^32 | 2^16 | 4 | -0.000036 | -0.000058 | -1.2e-05 |
| 2^26 top | 2^8 | 4 | -0.003446 | -0.005615 | -7.6e-06 |
| 2^28 top | 2^12 | 4 | -0.000295 | -0.000530 | -8.6e-06 |
| 2^30 top | 2^12 | 4 | -0.000189 | -0.000268 | -4.3e-06 |

  The deficit shrinks with the toy degree. Empirically the mean deficit is roughly proportional to ln(X2)/d_t, an empirical ansatz, not a derivation: the means fall from a few hundredths (k = 4) to a few 10^-4 (k = 12) and within 10^-4 of zero (k = 16). Single k = 8 rows scatter widely (lowest -0.0459); each stays inside the allowance after driver scaling. Production's driver is ln(2^76)/2^20 = 5.0e-5 per pair, 58x below the k = 12 rows and 4.5x below the k = 16 rows. The top-window control (base near r/X2, products a*b near r, maximal cross-product wrap) matches the low-window rows. No modular coincidence is excluded anywhere. The repeat-fraction diagnostic is uniform-like for k >= 12: the measured fraction of pairs hitting an already-marked coset is within 0.01 of the independent-draw value 1 - (1 - e^{-K})/K (0.4304 at K = 1.255), while the k = 4 rows sit above it. H1 remains a heuristic: the census reaches m = 2^32 at k = 12 and k = 16 (X2 up to 2^23.2, against production m = 2^234.86, X2 = 2^76, and products a*b = 2^152 below r = 2^254.9 while the toy products sit near or below r depending on the window), so the production statement is the explicit 0.0005 allowance. Driver-scaled to production (the ansatz applied to every row with k >= 8), no measured row is below -1.1e-04, inside the 0.0005 allowance (at m = 2^24 the k = 16 scatter matches the binomial noise of a coverage fraction over m/6 cosets, sd 2.7e-04); with the gate, an allowance larger than needed costs bits, not validity.
- *End-to-end agreement.* The ec toy of section 8 checks, every trial, whether the sampled t lies
  in the exactly-enumerated R of its emitted sets, and every observed candidate/miss event agrees.

**H2 (equidistribution).** The characters of the primes of I and of the integers b <= X2 spread
evenly over Z/d: class lists have between EMIN = 2^56 - 2^52 and 2^57 elements (proof.md s35), no character has more than 2^42
bases, and the walks need fewer than 2^61 bases per side at the EMIN bound (proof.md s43). The
interval I = (2^76, 2^76 + 2^70] holds about 2^64.3 primes by the PNT in intervals
(Rosser-Schoenfeld bounds; the sieving charge covers all 2^70 + 2^43 integers), a margin above 7. Exact class lengths and base counts are computed
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
2^23 waves; an owner pool above 2^57 (mean below 2^56.41 at proof.md s43; checked per wave, s18) over 2^60
owners; the owner candidate list above 2^16 (mean 2^-17.8), a candidate or reply link above 2^8
(mean below 2^-37.8) and the node comparison count above 1.01 x mean (Bernstein) are each below
2^-40 over all nodes; false key matches: proof.md s43 (3 N_b N_g/2^192 < 2^42.19 against the 2^44 cap, Chebyshev
below 2^-40); a re-walk chain carrying one key twice is below 2^-81, and one
key at more than 4 steps is below 2^-45 over all chains. The success line subtracts 4*2^-40 for
the union of these.

**Diagonal event.** If t = u a/a' for base primes a, a' of the same character (chi(a) = chi(a')),
then at every common position b the baby and giant points share an orbit, so one base pair gives
more than 2^44 true matches and the run outputs fail at the candidate-count check (which precedes
any delivery). The current bound is 6*2^122/m < 2^-110 (proof.md s43).

**H5 (candidate multiplicity).** For a uniform t in R, the number nu(t) of emitted pairs whose
ratio lies in t's coset — the number of candidates a true hit produces. Two means enter, in the
correct directions: the expected total candidate count is E[#false] + E[nu(t)] = 3 N_b N_g/2^192 +
N_b N_g/M < 2^42.19 + 1.26 (proof.md s43), a conservative upper bound on the exact unconditional mean of the sum, and the cap tails for the total
are evaluated at it; the multiplicity of the true hit itself, given that t is covered, has
conditional mean N_b N_g/(M x coverage) <= 1.76 at the claimed coverage, and the per-run budget of
2^20 true matches is Markov-bounded by it (Pr[nu > 2^20 | covered] <= 1.76/2^20 < 2^-19.18).
Support: the census's multiplicity diagnostic counts nu exactly (no sampling) and reports the
second-moment ratio sum(nu^2)/(T + T^2/M) (the value 1.000 for independent uniform pairs), K = 1.255:

| toy m | d_t | primes | energy ratio (1.000 = independent) |
| ---: | ---: | ---: | ---: |
| 2^24 | 2^4 | 4 | 1.1645 (worst 1.2228) |
| 2^24 | 2^8 | 8 | 1.0129 (worst 1.0260) |
| 2^24 | 2^12 | 8 | 1.0009 (worst 1.0020) |
| 2^24 | 2^16 | 8 | 1.0001 (worst 1.0008) |
| 2^26 | 2^8 | 8 | 1.0174 (worst 1.0282) |
| 2^26 | 2^12 | 8 | 1.0008 (worst 1.0014) |
| 2^26 | 2^16 | 8 | 1.0003 (worst 1.0009) |

tracking H1's ln(X2)/d driver (production is 58x below k = 12); the ec toy's candidate counts
match the model at every size. Under H5's
Poisson-like multiplicity, the 2^20 true-match budget inside the 2^44 cap, the per-node 2^16
candidate slots and the 2^8 per-link candidate caps carry tails below 2^-40 each, and the union is charged in proof.md section 2. The distribution-free
Markov alternative does not establish the current tightly margined floor;
the present claim explicitly relies on H5.

The current conditional success floor and every failure deduction are in
proof.md section 2. Those exact values supersede this source's older
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

