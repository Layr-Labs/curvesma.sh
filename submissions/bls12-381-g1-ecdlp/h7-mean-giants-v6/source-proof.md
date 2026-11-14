# bls12-381-g1-ecdlp: additive Cheon with planned pair walks, a measured side-output kernel and synchronous wave transport

The current claim, ledger, memory and success floor are proof.md's. Every cost, size and failure figure is rounded up, every success figure rounded down, by the script.

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
EMIN = 2^56 - 2^52 (proof.md s35) and EMAX = 2^57 and at most 2^43 bases per character (H2, checked by the sieves). K, N_b, N_g, GMAX, F: proof.md (s2, s30, s31).

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
   to the 2^27 generating nodes of class -c (at most 2^16 per node under the 2^43 bound), as a
   list of at most 2^16 (base, index) pairs per node, transferred once.
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
   key and compared against the stored baby records there; they are not stored. The giant side runs
   in rounds and stops at the first verified s (proof.md s44).
6. **Collision.** A truncated-key match is a candidate (key, baby id, giant id), appended at the
   owner node (2^16 slots; a 2^16+1st fails the node). After each giant round (proof.md s44), the barrier tree
   sums the candidate counts; if the running global count exceeds 2^45 the run outputs fail. Otherwise
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
direct-mapped head array of 2^59 slots (8-byte links, 2^62 bytes) indexed by key bits 69..127 (proof.md s17),
chained through a pool of 2^57 32-byte records (2^62 bytes; more than the mean load, below 2^56.25 at proof.md s43,
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

CAP1 = 2^24 records on each generating node -> relay link and CAP2 = 2^21 on each relay -> owner link (mean loads per wave at most 2^23 and 2^20 under H3); the sender maps both parities of every link buffer as windows, counted in its per-node total (evidence.py memory terms).

A **wave** is 2^33 lockstep steps of a generating node's 2^19 lanes (proof.md s42; giant waves 2^31 steps), hence at most 2^53 records
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

Phases: 2^13 plan-copy rounds (section 3), the baby waves (at most ceil(2^16 LANE_OPS / 2^33) + 2 under the H2
bounds: 2^16 bases per node, LANE_OPS steps per base per lane, 2^33 steps per wave, drain
and rounding waves; giant waves of 2^31 steps), a phase barrier, then the giant rounds, each with its waves, candidate round and close barrier (proof.md s44). The
ledger bounds the number of barriers by WAVES = 13836394 (both phases, plan rounds, 64 spare).

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

After each giant round (proof.md s44) and its status pass (running candidate total <= 2^45, else fail), one *candidate
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
Per-link loads are tiny in the mean (2^44.03 false candidates over 2^80 owner-relay pairs), and
any link above its cap fails the run (H3 tail inside the 4*2^-40 budget). The simulation runs
both rounds and asserts that every candidate is resolved exactly once.

## 5. Instruction listings (AArch64 scalar, 64-bit registers)

### 5.2 Per-operation listings (pointers)

Kernel, pair, root and transport listings: proof.md s1, s9-s30 (generated and executed). The table strand step stands for re-walk rebuilds (tables: proof.md s33): beyond the reference addition, load 8, store 7, batch buffers 6, inversion share, loop 5 (TABLE_BK, with the s1 +21 for v and n). Residual unlisted bookkeeping stays inside the explicit allowance row (proof.md s3).

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

## 6. Ledger

Rows, totals, preprocessing and memory: proof.md s4.
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

**H1' (prefix coverage; H1 is its last case).** For every round prefix j of the production schedule (proof.md s44), |R_j|/m >= 1 - e^(-K_j) - 0.0005 with K_j = 6 N_b |G_<=j|/m; the full set has K >= 9/2 (proof.md s2). This is a hypothesis, supported as follows and not proved.

- *Exact census in production round order.* `evidence.py` mode `census` with `nbar` builds the construction on toy primes r = 6 c 2^k + 1 (c odd, smooth, deterministic; skips r with 2 in H_t): N_g least with 6 ceil(Q N_g) N_g >= K m (Q = 13/50), class lists over X2 = d e with e = N_g/(nbar d), so each character has about nbar giant bases; baby bases are the primes after X2 in increasing order until N_b points (last base partial); giant bases are the next primes, complete class walks, until N_g points. chi is computed for every b <= X2 by a Pohlig-Hellman dlog sieve. The giant bases are then split as production splits them: the bases of character c are dealt in increasing order to gpc generating nodes of class -c and round j is the j-th base of every node, so the round-j prefix holds each character's first min(n_c, j gpc) bases. After every round the U-cosets of all ratios (a b)/(a' b') marked so far are counted in a byte array of m/6: exact, every modular and rational coincidence included, no sampling. The same final set is also swept in plain increasing-prime order with the same number of bases per step (control). Offline rows at K = 4.5 = K_T (52 runs, 4 primes per row; prefixes run from the first round to the full set; diff = coverage minus 1 - e^-K, 'full' its mean at the full set, 'lowest' over all rounds of all primes; 'max gap' = largest difference between the two orders at one round):

| toy m | d_t | nbar/gpc | rounds | full diff | lowest diff, prod / incr | eps | lowest + eps K e^-K, prod / incr | max gap |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 2^28 | 2^4 | 8/1 | 12-17 | -0.00575 | -0.06068 / -0.05906 | 0.0932 | -0.00696 / -0.00671 | 0.00328 |
| 2^28 | 2^6 | 8/1 | 13-17 | -0.00115 | -0.01188 / -0.01167 | 0.0242 | -0.00010 / -0.00016 | 0.00021 |
| 2^28 | 2^8 | 8/1 | 15-18 | -0.00013 | -0.00169 / -0.00170 | 0.0047 | +0.00004 / +0.00004 | 0.00009 |
| 2^28 | 2^8 | 16/2 | 13-17 | -0.00017 | -0.00167 / -0.00174 | 0.0047 | -0.00000 / +0.00003 | 0.00021 |
| 2^28 | 2^10 | 8/1 | 20-26 | -0.00012 | -0.00131 / -0.00135 | 0.0028 | -0.00002 / -0.00009 | 0.00012 |
| 2^30 | 2^4 | 8/1 | 13-15 | -0.00597 | -0.05720 / -0.05611 | 0.0951 | -0.00615 / -0.00540 | 0.00214 |
| 2^30 | 2^6 | 8/1 | 13-16 | -0.00219 | -0.02382 / -0.02439 | 0.0413 | -0.00078 / -0.00107 | 0.00100 |
| 2^30 | 2^8 | 8/1 | 16-18 | -0.00246 | -0.04519 / -0.04440 | 0.0408 | -0.00378 / -0.00353 | 0.00156 |
| 2^30 | 2^8 | 16/2 | 13-15 | -0.00246 | -0.04483 / -0.04431 | 0.0410 | -0.00354 / -0.00341 | 0.00102 |
| 2^30 | 2^10 | 8/1 | 17-19 | -0.00014 | -0.00158 / -0.00160 | 0.0030 | -0.00004 / -0.00006 | 0.00007 |
| 2^32 | 2^8 | 8/1 | 16-18 | -0.00045 | -0.00537 / -0.00528 | 0.0095 | +0.00000 / +0.00000 | 0.00010 |
| 2^32 | 2^10 | 8/1 | 17-19 | -0.00010 | -0.00125 / -0.00120 | 0.0022 | -0.00003 / +0.00000 | 0.00008 |
| 2^32 | 2^12 | 4/1 | 12-14 | -0.00003 | -0.00047 / -0.00045 | 0.0008 | -0.00003 / -0.00001 | 0.00004 |

- *What the deficit is.* Two emitted pairs give the same integer ratio only inside one base pair (a and a' are primes above X2 >= b), and there exactly when b/b' repeats: a pair (b, b') with gcd g repeats an earlier one if and only if some t < g has chi(t) = chi(g). The column eps is the measured share of such pairs (exact per class pair, 3000 sampled class pairs per prime). If the distinct ratios fall like uniform draws, coverage is 1 - e^(-(1 - eps) K): a deficit of about eps K e^-K, at most eps/e, under either prefix rule, because eps is a property of class pairs and not of the order of bases. Adding eps K e^-K back leaves at most 0.0001 in absolute value in every row with eps < 0.03 and at most 0.00003 at m >= 2^32; rows with eps above 0.03 (all d_t = 2^4 rows, and toy primes where a small prime has a low-order character) keep a second-order remainder up to 0.007. Production order and increasing order differ by at most 0.00010 at any round at m >= 2^32.
- *Production eps.* chi(g) = chi(t) if and only if g^m = t^m mod r, so the set D of g with an earlier equal character is computable. `eps_lit` (ledger output eps_lit_bound) sums 1/g^2 over D up to 2^15 and adds the tail 2^-15: eps <= 3.57e-05; the same sum taken offline to 2^22 is below 8.9e-6 (the first such g is 803: chi(803) = chi(8)). The share of class members divisible by g is taken as 1/g (H2 at scale X2/g). The literal-repeat deficit at production is then below 1.32e-05 at every prefix, against the 0.0005 allowance.
- *What this does not show.* The toy runs have 12 to 26 rounds, production up to 2^16; toy m reaches 2^32 against 2^234.86 and toy X2 2^15 against 2^76. That the distinct ratios cover like uniform draws at production size, for every prefix, is the heuristic content of H1'; the census tests it at toy size only, and the 0.0005 allowance is the stated margin.
- *End-to-end agreement.* The ec toy of section 8 checks, every trial, whether the sampled t lies
  in the exactly-enumerated R of its emitted sets, and every observed candidate/miss event agrees.

**H2 (equidistribution).** The characters of the primes of I and of the integers b <= X2 spread
evenly over Z/d: class lists have between EMIN = 2^56 - 2^52 and 2^57 elements (proof.md s35), no character has more than 2^43
bases, and the walks need fewer than 2^63 bases in all at the EMIN bound (proof.md s43). The
interval I = (2^76, 2^76 + 2^70] holds about 2^64.3 primes by the PNT in intervals
(Rosser-Schoenfeld bounds; the sieving charge covers all 2^70 + 2^43 integers), a margin above 3. Exact class lengths and base counts are computed
by the sieves at runtime; a shortage or an exceeded bound stops the run before any chain starts
and is a failure event under H2, which excludes it. Toy evidence: the class lengths of the ec
toy's 32 classes and the census's class counts at every scale. H2 affects only charged small terms
and the wave/buffer bounds, never the sets.

**H3 (keys).** The 192-bit keys of the 2N records behave as independent uniform strings across
U-orbits: routing digits t, u and head slots are binomial, and a cross-orbit key equality has
probability 2^-192. Support: the toy's occupancy chi-square matches the binomial expectation
(scaled statistic about 1), and its comparison and false-match counts match the model. Tails
under H3, each already multiplied by the number of links, nodes or slots: a generating-node link
above CAP1 = 2^24 in a wave (mean at most 2^23) is below (e/4)^(2^23) over 2^77 links and 2^24
waves; a relay link above CAP2 = 2^21 (mean at most 2^20) below (e/4)^(2^20) over 2^80 links and
2^24 waves; an owner pool above 2^57 (mean below 2^56.25 at proof.md s43; checked per wave, s18) over 2^60
owners; the owner candidate list above 2^16 (mean 2^-15.9), a candidate or reply link above 2^8
(mean below 2^-35.9) and the node comparison count above 1.01 x mean (Bernstein) are each below
2^-40 over all nodes; false key matches: proof.md s43 (3 N_b N_g/2^192 < 2^44.03 against the 2^45 cap, Chebyshev
below 2^-40); a re-walk chain carrying one key twice is below 2^-81, and one
key at more than 4 steps is below 2^-45 over all chains. The success line subtracts 4*2^-40 for
the union of these.

**Diagonal event.** If t = u a/a' for base primes a, a' of the same character (chi(a) = chi(a')),
then at every common position b the baby and giant points share an orbit, so one base pair gives
more than 2^45 true matches and the run outputs fail at the candidate-count check (which precedes
any delivery). The current bound is 6 x (baby bases) x (giant bases)/m < 2^-108 (asserted; proof.md s43).

**H5 (candidate multiplicity).** For a uniform t in R, the number nu(t) of emitted pairs whose
ratio lies in t's coset — the number of candidates a true hit produces. Two means enter, in the
correct directions: the expected total candidate count is E[#false] + E[nu(t)] = 3 N_b N_g/2^192 +
N_b N_g/M < 2^44.03 + 4.50 (proof.md s43), a conservative upper bound on the exact unconditional mean of the sum, and the cap tails for the total
are evaluated at it; the multiplicity of the true hit itself, given that t is covered, has
conditional mean N_b N_g/(M x coverage) <= 4.56 at the claimed coverage, and the per-run budget of
2^20 true matches is Markov-bounded by it (Pr[nu > 2^20 | covered] <= 4.56/2^20 < 2^-17.81).
Support: the census's multiplicity diagnostic counts nu exactly (no sampling) and reports the
second-moment ratio sum(nu^2)/(T + T^2/M) (the value 1.000 for independent uniform pairs); at K = 4.5, Q = 13/50 (means over 4 primes): m = 2^24, d_t = 2^4: 1.0748 (worst 1.1083), m = 2^24, d_t = 2^8: 1.0072 (worst 1.0121), m = 2^24, d_t = 2^12: 1.0006 (worst 1.0009), m = 2^24, d_t = 2^16: 1.0001 (worst 1.0004), m = 2^26, d_t = 2^8: 1.0068 (worst 1.0095), m = 2^26, d_t = 2^12: 1.0002 (worst 1.0005); falling with the toy degree, as eps does; the ec toy's candidate counts
match the model at every size. Under H5's
Poisson-like multiplicity, the 2^20 true-match budget inside the 2^45 cap, the per-node 2^16
candidate slots and the 2^8 per-link candidate caps carry tails below 2^-40 each, and the union is charged in proof.md section 2. The distribution-free
Markov alternative does not establish the current tightly margined floor;
the present claim explicitly relies on H5.

The current conditional success floor and every failure deduction are in
proof.md section 2. Those exact values supersede this source's older
numbers. All capped work, including aborts, remains in the full ledger.

## 8. Evidence scope

evidence.py exposes ledger, planner, ec, proto and census modes. Ledger checks exact accounting and identities, not the production hypotheses. Planner uses real character lists on toy fields. The ec mode uses actual small-curve points and verifies every recovered scalar and coverage event, with model-only occupancy, pool and comparison caps. Proto simulates charged parity-buffer transport with overflow and wrong-generation controls; it is not a physical cluster. Census counts all ratios at its toy parameters. All finite evidence has only its stated scope.

## 9. W(n), limitations, context

The algorithm reads only A_1 and A_d, so W(n) = W(2^20) (bound: proof.md) for every 2^20 <= n < r; nothing is claimed for smaller n.
The cluster model is the predecessor author's (section 4.1); every use is a charged store, load, wake or sweep. H1-H5 are disclosed heuristics with their support and extrapolations stated: the census reaches
m = 2^32 (production 2^234.86), the ln(X2)/d trend is an empirical ansatz, the planner fractions are measured at toy spans (production far larger), and the ec toy is 2^22.8-2^25.9 in size.

