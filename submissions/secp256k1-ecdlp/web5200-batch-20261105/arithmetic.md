# Conservative 64-bit computation and storage certificate

This file prices the exact secp256k1 algorithm in proof.md. Let m=h=2^127.
The only priced primitive is one scalar AArch64 instruction equivalent.
No SIMD, quantum operation, arbitrary-width arithmetic primitive,
unpriced lookup, or wall-clock parallel speedup is used. We bound total
work, not expected work, and include every unsuccessful target.

## 1. Local finite arithmetic

Each secp256k1 field value occupies four 64-bit limbs. Canonical operands
are in [0,p); a fifth limb accommodates addition carry and correction.
Bound a modular addition or subtraction by A=1024 scalar instructions.
One implementation performs at most six linear passes over five limbs
for loads, stores, ADC/SBC carry or borrow, comparison, selection and
correction, allowing sixteen instructions per limb per pass. That is
480 instructions; an additional 544 covers copying operands, fixed
addresses, loop counts, local calls and returns. Adding two canonical
values requires at most one subtraction of p; subtraction requires at
most one addition of p. Reduction is not an unbounded loop.

For multiplication, initialize result=0 and running=u. Scan all 256
bits of v; if the bit is one, add running to result modulo p, and
double running modulo p. Two A calls per iteration cover both paths.
Allow 64 instructions per iteration for extraction, branching, and loop
control, plus 1024 for initialization, copies and final local transfer:

    M = 256*(2*1024+64)+1024 = 541696 < 2^20.

This algorithm maintains result and running in canonical form, and its
binary expansion proves multiplication modulo p. Leading zeros do not
alter the fixed iteration bound.

Inversion by exponent p-2 uses at most 512 multiplications with a known
256-bit exponent. Allow 64 extra instructions per bit and 1024 setup:

    I = 512*M+256*64+1024 < 2^29.

The point law in proof.md uses at most one such inverse, sixteen
other field multiplications (squaring included), and thirty-two modular
additions/subtractions. A further 4096 instructions covers coordinate
comparisons, the identity/opposite/equal tests, local copying, tag
handling and control. Therefore a COMPLETE affine addition or doubling
has cost at most

    G = I+16*M+32*A+4096 < 2^30.

These bounds include exceptional cases by taking the maximum branch.
The actual coefficient a=0 enters the doubling numerator. Computing
3*x^2 can use additions already within the allowance. No per-addition
inversion is omitted. No reference 6M convention is being used to charge
this implementation: every instruction above enters R, and H=0.

All arithmetic in these routines uses only fixed-size local scratch,
not a table-length allocation. Allow 1 MiB of local code, registers,
stack, point scratch, parsing and scalar-control buffers. The code is
fixed nested loops on four/five limbs, not an exponentially large
unrolled program. Four-limb scalar additions, comparisons and reduction
are bounded by 4096 instructions each. The giant-step index has 128
bits; the signed reconstructed scalar has at most 257 bits and fits
five limbs before the single modular correction.

## 2. Global memory without a 128-bit unit-cost address

There are two table buffers, each of m records of 128 bytes. A record
contains four limbs of x, four of y, two of j, and padding. An ordinary
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
a read and write. B,D include local 256-bit comparisons, four-limb
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

The cold baby table requires at most m point additions and m record
writes, charged by m*(G+D). The full giant traversal requires at most
h point additions and h searches, charged by h*(G+V). It is legitimate
to charge an extra unused last addition: we are proving an upper bound.
Since G,V<2^30 and D=2^22,

    m*(G+D) < 2^158,
    h*(G+V) < 2^158.

The sorting bound is S<2^158. The routing/array initialization is
less than 2^149. Explicitly charge an additional
T_loop=4096*(m+h)=2^140 instructions for the repeated baby/giant loop
increments, identity/tag checks, loop bounds and controller bookkeeping,
including iterations with no candidate. These recurring costs are not
assigned to a fixed setup allowance. Reserve an additional 2^50
instructions for everything
not already in those four terms: canonical parsing, checking the curve
equation, [r]Q, constructing and negating [2m]P, final [c]P, comparison
and serialization of output, scalar/index setup, initial identity handling,
and any fixed local scratch initialization. Two 256-bit scalar
multiplications need at most 1024 point operations and thus fewer than
2^40 instructions; even doubling that bound plus the listed fixed
work is far below 2^50. At most one final verification occurs: its
unexpected failure terminates the algorithm, and the correctness proof
shows it cannot fail for a genuine table hit on a valid input.

Therefore the complete worst-case instruction count R satisfies

    R < 3*2^158 + 2^149 + 2^140 + 2^50 < 2^160.

Divide by the PROTECTED normalization C=1000. With H=0 separately
priced group operations, W=R/1000 < 2^151. There is no double-counted
reference kernel or candidate-chosen price. The exact declared score
is the rational 151/1, a conservative log2 work bound. No floating-point
logarithm or experimental timing is required for this certificate.

Preprocessing is at most the total budget W<=2^151 and is not amortized.
Peak storage is at most 2^144 bytes as shown above, including code and
routing. Searched nonuniform advice is zero; public p,r,a,b,P and the
fixed algorithm are the excluded profile constants/design. The arrays
and their construction ARE charged, regardless of their independence
from Q. No randomness, retries or amplification cost is omitted because
none is used. The probability statement concerns exactly the deterministic
coverage in proof.md, independent of this accounting upper bound.
