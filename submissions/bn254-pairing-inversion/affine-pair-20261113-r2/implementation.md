# Halfword implementation with explicit distributed storage

All arithmetic below is scalar and bounded; all processor work is summed.
The memory requirement is theoretical and deliberately enormous. It is not
an instruction to construct the machine or execute the attack on the author's
computer. Communication latency and hardware construction are not scalar
instruction counts; all software transfer, initialization and coordination are.

### Residues and polynomial arithmetic

A field residue is stored in sixteen 16-bit digits, each in an ordinary 64-bit
word. Digits are always in [0,65535]. For addition, make a sixteen-digit pass
with carry, then compare with p from the high digit, and subtract p once if
needed. Since p<2^254, an addition fits in sixteen digits. A subtract uses
borrow and conditional p addition. An explicit pass requires at most eight
scalar operations per digit: two loads, arithmetic/mask/shift or comparison,
a store and bounded addressing/control. Three passes and at most 128 setup
operations fit in 512 instructions. Doubling is the same modular add.

For multiplication, maintain q=0,d=A, and for all 256 bits of B add d to q
if the bit is one and double d, reducing after each operation. Always charge
two 512-instruction adds and 64 instructions of bit/control overhead per bit,
plus 4096 setup/copy instructions; this is less than 2^19. The elementary
invariant q=A*(low bits of B), d=A*2^k proves modular multiplication exactly.
The same procedure works modulo the protected prime r, with sixteen digits.
No machine instruction performs a 254-bit multiply for unit price.

Multiply Fp12 elements by convolution: 144 residue multiplications and modular
accumulations, forming degrees zero through 22. For k=22 down to 12, replace
c*w^k by 18c*w^(k-6)-82c*w^(k-12). The small constants use fewer than 24
doublings/adds per high term. Canonical reductions occur throughout. Total
arithmetic is below 2^27; permitting all halfword-array loads/stores, scratch
initialization and loops, use 2^28 for a complete Fp12 product. Binary powers
with at most 256 exponent bits use at most 512 products. Comparisons use all
192 halfword words, not a truncated digest. Converting the 16-digit coefficient
to/from a protected 64-hex-digit string is finite and explicitly accounted.

### Network and sorting schedule

Allocate n=2^128 record processors in a Boolean cube of dimension 128. Colocate
the controller with coordinate zero and share its cube endpoints; its additional
root heap/work is counted separately. Root routes therefore need at most 128
cube hops and do not use an unpriced extra attachment link. A
processor coordinate uses two 64-bit words; local addresses are ordinary 64-bit
addresses. Per-node allocation, including register/control state, code and
neighbor queues, is at most 16 MiB. Neighbor dimension d is a local fixed
link endpoint, selected by a bounded index. This does not posit a giant array
accessible by a free 128-bit pointer. Physical links are part of the declared
machine; they contain no precomputed field answers.

A key has 192 digit words. With side, four-word counters/indices and fixed
headers, the record and its route information occupy at most 256 words.
Serialize and reconstruct each packet using at most 32 scalar operations per
bit across the sender and receiver, plus fixed bounded headers and endpoint
control. A full hop is below 2^20 instruction equivalents. These operations
include loads, digit/bit tests, stores, counters and software transport; link
hardware is not asked to compare keys or perform arithmetic. Endpoints sleep
or halt between scheduled events instead of executing unbounded polling loops.
Code/data loading, wake/control and the finite schedule are counted. These are
claims for this specified machine, not arbitrary network stacks or OS kernels.

A spanning-tree bootstrap assigns coordinates, initializes all local memory and
loads a bounded loop program (less than 256 KiB). Distribute that program over
the spanning tree: 128 packets of 2048 bytes per edge suffice, at fewer than
n*128*2^20=n*2^27 aggregate instructions. Include the fixed uniform bootstrap
receiver's loading and state, local zeroing, identifier construction, port setup,
all control and acknowledgments in the larger n*2^31 initialization allowance.
This explicitly charges distribution as well as local copying. No preexisting key list is
part of that initialized state. Generate the two progressions at the root and
send each record along at most 128 coordinate-changing hops. Serialize these
deliveries to avoid uncontrolled congestion. One GT product, one route, all
packet handling and index arithmetic per record fit in 2^30. This charges
the complete n-record generation and population of storage from scratch.

Sort using the standard explicit bitonic comparator schedule: for each
k=2,4,...,n and j=k/2,...,1, form t=i XOR j and handle the pair only at the
endpoint i<t. When (i AND k)=0, retain the smaller record at i and larger at t;
otherwise retain the larger at i and smaller at t. Both endpoints receive
their explicitly assigned records; they do not both select a minimum.
The last mask k=n uses three
words. There are 8256 layers. An ordered comparator uses full keys, then side,
then index; equal values cause no ambiguous identification. A layer exchanges
records over disjoint cube edges, then uses a bounded spanning-tree barrier
and next-layer broadcast. At most n full packet hops plus 2n much smaller
control hops are needed. Charge n*2^21 per layer: full hops use <n*2^20,
local comparisons/copies are <n*2^16, and the smaller control messages are
bounded by n*2^19 in aggregate. Tree depth changes time, not the number of
edge messages; all control work is counted. There is no cost-free barrier.

Why the network sorts: two oppositely sorted equal-size lists form a bitonic
list. Pairing positions half a list apart puts the minima in one half and
maxima in the other; in a thresholded zero-one sequence each half remains
bitonic and all first-half values are <= all second-half values. Recursing
sorts both halves. Thresholding commutes with each comparator, so a general
output inversion would imply a binary inversion. Induction over list size
proves the schedule sorts every ordered input, including duplicates. Thus
the group-key/side boundary argument of proof.md applies to the actual sort.

For the final scan the root schedules adjacent-position queries in index
order. For each adjacent pair, two requests and two returned records form four
routed packets; allow four more full-packet equivalents for completion/acks.
Each follows at most 128 cube edges, costing at most 8*128*2^20=2^30.
Charge 2^31 per position including comparisons, address/counter operations and
copies. Queues remain bounded by serialization. The root stores the first
cross-side equality while continuing its bounded scan. A further 2^28 allowance
per position covers selection bookkeeping and final record reporting.
No ideal global lookup, random access or free comparison is used.

### Finite sampling, root work, parser and storage

At most 4096 candidate random scalars of 254 bits are read for four capped
rejection samplers. Charge up to 64 scalar instructions per random bit,
plus comparison and state updates; this is below 2^28. Fresh random bits
are the ordinary algorithmic coin source allowed by the probability space,
not an oracle providing logs. Failed samplers return failure within this cap.
Four starting/step exponents, input subgroup validation, scalar arithmetic and
final verification are bounded binary powers and four-limb operations.
All fixed work, including the point-output extension in output.md when present,
fits in the deliberately generous 2^52 instruction allowance.

Parse the sole z object with a fixed byte-state machine. Verify exactly twelve
64-hex-digit coefficients, their ranges, complete syntax, no duplicate/extra
keys and EOF. Whitespace is allowed at legal JSON positions and consumed in
full; charge 8192*(L+1) for actual length L, including transfers, digit
conversion, state changes and output encoding. No buffer of length L is
needed: only the fixed twelve coefficients and parser state are retained.
For malformed data fail rather than return a fabricated answer. This parser
does not invent a maximum L. Its fixed numerical score uses only the
explicit canonical-input interpretation, while R(L) retains arbitrary padding.

Sixteen MiB per record node includes code, arithmetic scratch, record buffers,
incoming/outgoing link storage, registers/counters and all retained randomness;
root allocation is below 64 MiB. Thus n*2^24+2^26<2^153 bytes, below the
claimed 2^154. Every processor's program copy is counted, as are all packet
buffers. Preprocessing is constructed and charged anew for every target.
