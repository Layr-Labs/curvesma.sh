# Conservative 64-bit computation and storage certificate

This file prices the search part of the auxiliary-powers G1 algorithm in proof.md. The extra input-dependent term is appended below. Let m=2^127 and h=2^126.
The only priced primitive is one scalar AArch64 instruction equivalent.
No SIMD, quantum operation, arbitrary-width arithmetic primitive,
unpriced lookup, or wall-clock parallel speedup is used. We bound total
work, not expected work, and include every unsuccessful target.

## 1. Local finite arithmetic

Each field value uses six 64-bit limbs; a seventh accommodates addition
carry and modular correction. Bound modular addition/subtraction by A=2048
instructions: at most six passes over seven limbs at sixteen instructions
per limb (672), and a further 1376 for all fixed copying, calls, addresses,
control and correction. Canonical input sums need at most one subtraction
of p; a difference needs at most one addition. No unbounded reduction occurs.

Use 384-bit shift-and-add multiplication, including leading zero bits.
Each bit allows two modular additions and 64 control/extraction instructions;
allow 2048 fixed setup and transfer instructions. Then

    M=384*(2*2048+64)+2048=1599488 < 2^21.

Binary exponentiation to p-2 uses at most 768 such multiplications and
384 iterations of exponent-bit handling:

    I=768*M+384*64+2048 < 2^31.

The complete affine point law in proof.md requires at most one inverse,
16 other field multiplications and 32 field additions, plus 8192 instructions
for six-limb comparisons, exceptional cases, local copies and control:

    G=I+16*M+32*A+8192 < 2^32.

This deliberately loose bound covers a=0 and b=4, identity and opposite
points, and zero doubling denominators. Fermat inversion applies only to
nonzero field elements. All actual instructions enter R; H=0, so no
reference kernel is charged twice. Local scratch and fixed code fit in
1 MiB. Field arithmetic has six/seven limbs, while scalars have four/five;
the subgroup order has 255 bits. Signed scalar recovery fits five limbs,
and its range requires at most one modular correction. Scalar comparison,
addition and correction use at most 4096 instructions. Scalar multiplication
uses at most 512 point operations for a padded 256-bit scalar. The giant
index fits 128 bits. No arbitrary-width arithmetic primitive is used.

## 2. Global memory without a 128-bit unit-cost address

There are two table buffers, each of m records of 128 bytes. A record
contains six limbs of x, six of y, two of j, and padding. An ordinary
flat single-processor AArch64 address space cannot directly name these
astronomical arrays. We do not assume it can. The hypothetical classical
machine can instead distribute the arrays over a binary tree of storage
processors and charge ALL participating processors' instructions.
This is the permitted total-work, rather than elapsed-time, parallel model.

Use one leaf per record in both buffers, hence 2m=2^128 leaves and fewer
than 2^129 nodes in a complete binary tree. A leaf has its 128-byte
record and bounded scratch. At an internal node there are local links
to its parent and two children. An address is a buffer bit followed by
a 127-bit index. Route it one bit per edge, at most 128 hops in each
direction. Address manipulation uses two or three 64-bit words; no
large address is put in one architectural register. Each node uses
only a small local address space. There is no precomputed record index
structure that contains a secret answer.

The fixed classical interconnect assumption is a bounded-degree tree
whose edges expose two directional, bounded memory-mapped mailboxes.
Each mailbox has a valid flag and a local level-triggered receive event.
The controller serializes record transactions, so no unbounded queue,
congestion resolution, global scheduling oracle or free index service is
needed. A storage node waits on its three fixed local ports with a
charged blocking/wait instruction and processes one message when its
local event is set. Send, receive, acknowledgment and return follow the
encoded address bits. No processor spins or polls while another works.
This is a hypothetical instruction-equivalent parallel machine with
ordinary local word transfers; it is not a claim that an existing cluster
has this size. Physical latency is not a scalar-instruction metric.

Messages use at most twenty-four 64-bit words: the three-word routing
label, operation/return flags, and at most sixteen data words. For each
hop charge at most sixteen instructions per word for sender/receiver
loads, stores, addressing and copying, or 384 instructions. Additionally
allow 640 instructions for valid-flag/acknowledgment changes, memory
barriers, bounded wait/event entry and exit, next-bit extraction, one of
three local-port selections, and the fixed forwarding/return control.
Thus a hop uses at most 1024 instructions TOTAL across its endpoints,
including coordination and notification. The synchronous protocol needs
at most one request/response event per directed message hop; its fixed
control has no data-dependent retry loop. Mailbox storage is counted at
both endpoints. Startup, tree-label assignment and local-port state are
included in the node-initialization budget below. No new arithmetic or
arbitrary-width global-memory primitive is assigned unit cost.

Routing a full record or a request-and-record response therefore costs
less than 2^20 instructions (256 hops * 1024 plus request/header
allowance). Reserve B=2^22 for a complete lexicographic x-key comparison
requiring two record reads, and D=2^22 for a complete record move with
a read and write. B,D include local six-limb field comparisons, four-limb
indices, global transfers, and branch/control overhead. Lookup never
receives a free hash-table oracle. Field arithmetic is local after the
point record is transferred, and its own limb traffic was charged in
Section 1; these table-transfer bounds do not replace that traffic.

Reserve 16 KiB at every storage-tree node for fixed forwarding code,
links, labels, record(s) and scratch. The small fixed routing program
needs fewer than 1024 instructions, occupying less than 4 KiB, and its
bounded buffers/labels fit in the remaining 12 KiB. With fewer than
2^129 nodes, this is less than 2^143 bytes. Separately add 1 MiB for
the controller. The total is below 2^143+2^20<2^144.
This includes code replicas and memory addressing structures, not just
the point coordinates. No physically available memory of this size
is claimed.

Initialize all nodes and both buffers explicitly with one serialized
depth-first traversal from the controller. The fixed less-than-4-KiB
routing program and each child's label are sent locally down every
tree edge once; acknowledgments return along that edge. Nodes clear
their bounded buffers and install local-port state before passing the
program/label to their children. Even 512 code words copied at sixteen
instructions per word, plus label/control transfer, bounded zeroing
and acknowledgments, fits in 2^20 instructions per node summed across
senders and receivers. No global deployment service is free. The total tree/array
initialization cost is less than 2^149. This budget is charged in full
as cold preprocessing; storage/code deployment work is not free advice.

### Finite boot and mailbox protocol within the existing reserve

The declared initial machine includes a uniform fixed receive-loader in
each node's reset state. Its at most 64 instructions occupy 256 bytes of
the already reserved code area. This is public fixed machine/program
state, not advice. To avoid crediting free installation, charge another
1024 instructions per node for its initial placement within the existing
2^20-per-node budget. The runtime routing program together with this
loader is at most 1024 instructions (4 KiB), already in the 16 KiB node cap.

A parent sends the at-most-512-word program as 32 sequential packets of
at most 16 data words. Each packet includes its fixed chunk number and
destination offset. The child checks that number, writes the data into
its bounded local program area, and acknowledges before the next packet.
There are no retries or queue growth in the stipulated error-free machine.
Allow four directed one-edge transfers per packet, each charged at 1024
instructions across endpoints: 32*4*1024=131072. Charge another 8192 for
label and final-ready messages, 32768 to clear at most 2048 local words,
and 65536 for loop control, loader completion, state setup, and all other
fixed steps. With the 1024 loader-placement charge, this is 238592 per
node, below 2^20. Count each parent-child installation against the child;
the root's corresponding fixed work is in the controller reserve.

At runtime a storage processor follows a finite receive state, checks a
local valid flag, reads one bounded message, decodes its next address bit,
forwards or accesses its record, acknowledges, and re-enters the blocking
receive state. The fixed wait/event hardware assumption supplies no
spurious wakeups and no polling; changed hardware would require new
cost accounting. The existing 640 control instructions per hop cover these
bounded branches, flag operations and wait entry/exit in addition to the
384 data-transfer instructions. Serialization permits one in-flight
transaction, so neither an unbounded queue nor a free scheduler is needed.
This elaborates the declared hypothetical interface; it does not certify
an actual available machine or remove that explicit assumption.

## 3. Sorting and lookup

Use bottom-up merge sort on two arrays of m=2^127 records. It takes
exactly L=127 merge passes. One pass uses at most m comparisons,
at most 2m full record moves, and 4096m instructions of remaining
multiword loop and interval bookkeeping. With B=D=2^22,

    S <= m*127*(B+2D+4096)+4096*127 < 2^158.

The comparison includes both inputs; the moves cover temporary/output
traffic. Empty tails and copy-only final runs use no more work. This
implementation trades practical efficiency for transparent accounting.

Binary search has at most 128 iterations. Each iteration can read
the chosen record, compare the full x-coordinate, and update its
two-word endpoints inside B+4096 instructions. Reserve one D for
reading the returned full record and 4096 for the result flag:

    V <= 128*(B+4096)+D+4096 < 2^30.

Here a central controller chooses the next multiword record address
and uses the routed read of Section 2. There is no uncharged random
access over a beyond-64-bit address space. The extra y comparison
and scalar sign/reconstruction fit the final-check reserve below.

## 4. Total budget

Cold baby generation and writes cost at most m*(G+D), with m=2^127,
G<2^32 and D=2^22; therefore the term is below 2^160. The giant
traversal costs at most h*(G+V), with h=2^126 and V<2^30, hence it is
also below 2^160. Sorting is below 2^158, conservatively below 2^160.
Node/bootstrap initialization is below 2^149. Repeated baby/giant loop
bookkeeping costs 4096*(m+h)<2^140, including identity/tag checks and
failed iterations. These repeated costs are not assigned to fixed setup.

Reserve another 2^52 instructions for canonical input parsing, curve and
subgroup checking, constructing and negating [2m]P, final scalar checking,
output, scalar setup and fixed scratch. Two 256-bit scalar multiplications
take at most 1024 point operations, below 2^42 instructions; even twice
that amount plus all listed fixed work fits the 2^52 reserve. Coordinates
have 96 hex digits but the scalar output still has 64. At most one final
candidate check occurs, and correctness ensures a genuine hit passes.

Thus every valid generated challenge, including a full failure, has

    R < 3*2^160 + 2^149 + 2^140 + 2^52 < 2^162.

All arithmetic is charged in instructions, H=0. With the protected C=2000
and 2000>2^10, W=R/2000<2^152. The submitted rational bound is 152/1;
this is conservative, not a measured runtime or a claim of equality.
Preprocessing is a subset of total W and is bounded by 2^152 as well.
Memory remains below 2^144 bytes: fourteen 64-bit payload words still fit
each 128-byte record, m and the two-buffer routing tree are unchanged,
and six-limb local scratch fits the controller's 1 MiB. No random coins,
restarts, amplification or searched nonuniform advice are used. All
profile constants are public. No table is free or reused across targets.

## 5. Auxiliary-powers input and the n-dependent total

All preceding costs price the search itself and remain conservative after
removing the zero secret. For each of n powers and T, allow 2^45 additional
instructions. A canonical point has fewer than 256 bytes, each visited
at most twice at 64 instructions per byte; fixed parser work is at most
2^16. Curve checking costs at most 32 field multiplications/additions.
A subgroup check has at most 512 G1 operations at G<2^32, below 2^41;
the remaining margin covers transfers, counters and per-record control.
The first retained point and n-counter fit the existing 1 MiB controller
allowance; even a four-limb counter is sufficient for n<r<2^255. T is
checked against the saved A_1. Actual powers are promised by the honest
input generation; they are not recomputed or treated as free searchable advice.

The streaming representation and bounded local port are explicit interface
assumptions. There is no free random read of a distant array: input bytes
are consumed in order, with their transfers charged in the per-point term.
The full compact input, if buffered by the interface, is separately counted
as 512(n+1) bytes. The M(n) bound in proof.md includes this plus all working
memory. No claim of actual storage at this size is made.

Thus total R(n)<3*2^160+2^149+2^140+2^52+(n+1)*2^45. For n=2^20 the
additional term is less than 2^66, keeping R(n)<2^162 and W(n)<2^152.
The earlier 2^144 memory estimate is for working/search storage; the
complete scoring-point input-inclusive bound is now 2^145 bytes.
