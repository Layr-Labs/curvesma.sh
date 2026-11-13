# Four-limb scalar implementation and charged finite network

This file specifies one intentionally slow implementation. All bounds are upper
bounds on scalar 64-bit instruction equivalents, including local loads/stores,
addressing and branch/loop control. No arbitrary-width operation has unit cost.
The target's polynomial field representation and group properties are protected
premises. Normal residues are used, not an unspecified Montgomery oracle.

## Local arithmetic

Represent a base-field residue by four little-endian 64-bit limbs. Because p<2^254,
the sum of two residues fits in four limbs. A modular add consists of four
carry-propagating adds, a full four-limb comparison with p from the most
significant limb downward, and at most
one four-limb subtraction. Modular subtraction uses the corresponding borrow
chain and conditional addition of p. Implement each limb using scalar loads,
adds-with-carry or comparisons, and stores; loop bounds are four. Fewer than
32 instructions per limb suffice for each of these three passes, plus fewer
than 128 for setup, branching, copying and return. Thus 512 instructions covers
one complete add or subtract, including alias-safe scratch. Doubling uses add.

To multiply a,b modulo p, run exactly 256 iterations of binary multiplication:
start t=0,v=a; for each bit b_i from least significant to most significant, add
v to t modulo p when b_i=1, then double v modulo p. Charge two full modular
adds regardless of the bit. Bit extraction, loop control and copying are bounded
by 64 more scalar instructions per iteration. Hence
256*(2*512+64)+4096 < 2^19 instructions cover a modular multiplication.
Inputs and outputs are canonical residues, and the invariant after k bits is
t=a*(b mod 2^k), v=2^k*a modulo p, proving the result.

For Fp12 multiply, compute the 144 coefficient products into the 23 coefficient
polynomial using modular multiplication and modular addition. Reduce degrees
22 down through 12: if the current high coefficient is c at degree k, erase it,
add 18c at k-6, and subtract 82c at k-12. This is exactly
w^12=18*w^6-82. Descending order processes all newly created high terms. Compute
18c and 82c using binary doubling/addition, at most 24 modular adds/subtracts
per reduced degree. 144*2^19 plus 144*512 plus 11*24*512 is below 2^27;
2^28 covers these operations, all coefficient transfers/loops, zeroing scratch
and an alias-safe final copy. Equality compares all 48 words; zero and identity
tests use those same full words. Scalar arithmetic modulo r uses the same four
limbs and bounds, since r<2^254. Exponentiation at most 256 bits uses at most
512 such products. No modular inverse is needed in the GT search.

## Actual address and communication representation

There are N=2^128 record processors arranged as a 128-dimensional Boolean cube.
The root/controller is colocated with cube vertex zero, sharing its link ports;
its extra local heap and work are separately bounded. Thus every root route
is a cube route of at most 128 hops, with no uncharged attachment link. Each has only local
64-bit addresses and at most 4 MiB of allocated memory. Its logical coordinate
is two 64-bit words; counters that can equal N use three words. A neighbor in
dimension d is reached through local port d, never through a 128-bit hardware
pointer or an uncharged global-array access. Cube links are part of the stated
finite classical hardware, not advice containing answers. This is an instruction
work model, not a wall-clock or physical construction-time claim.

A record contains 48 residue words, a tag, a two-word index and fixed headers,
at most 64 words total. Each port has bounded incoming/outgoing packet buffers.
Describe a hop as scalar software transfer/serialization of at most 4096 bits:
sender extracts bits/words, emits to the fixed adjacent link endpoint; receiver
reconstructs/checks/stores them, with explicit counters and header operations.
Charge up to 32 scalar operations per transferred bit across both endpoints,
plus 4096 fixed operations, bounded by 2^18 instruction equivalents per hop.
The chosen communication hardware performs no arithmetic or remote lookup;
link travel latency does not reduce or add scalar work. Endpoints block without
executing unbounded polling instructions. Protocol handling, instructions for
blocking/waking and control packets are included in the fixed allowance. If a
different machine requires more work, this is not a cost claim for that machine.

Initialization creates local zeroed buffers, loads a bounded loop program,
constructs the node coordinates/ports by cube broadcasts and clears state. Code
is less than 64 KiB (the stated loops, arithmetic and packet routines need no
large unrolling). Broadcast the program from the root over a cube spanning tree:
128 full 512-byte packets per tree edge suffice for 64 KiB, and all N-1 edges
are charged, below N*2^25 packet instructions. A fixed small bootstrap receiver
is part of each processor's uniform program; its loading and state are counted,
not a searched advice string. Local zeroing, endpoint setup, bootstrap copying,
all identifiers and acknowledgments are included in an aggregate N*2^30 budget.
This allowance covers program distribution, not merely local writes. Physical
capacity is stated, not treated as a precomputed answer. There are no initialized
key tables until the root generates and distributes them during this attack.

Generate the two progressions sequentially at the root, at most N GT products
and N record deliveries. Route a packet by clearing differing coordinate bits
one at a time, at most 128 hops. Deliveries are serialized, avoiding unbounded
queues; blocking time is not counted as a fictitious speedup. Each delivery costs
at most 128*2^18=2^25 plus control; each record's complete generation/delivery
is below 2^29 instructions, since the GT product is bounded by 2^28. All N
deliveries are charged. Random access is not used to construct the table.

In one sorting layer each cube link involved in a compare/exchange is disjoint
from the other pairs. Exchange full records in both directions; a pair sends
two bounded packets. Compare canonical keys, tags and indices locally, then
retain the correct record. Local full comparisons, bounded copies and global
layer-control messages are charged. Across all processors one layer is below
N*2^20: at most N full-packet hops (N*2^18), local comparison/copy below N*2^15,
and a spanning-tree completion/next-stage protocol below 2N*2^18. A tree has
N-1 edges; its depth is not another multiplicative N factor in total work.
The few root coordinates/counters are bounded separately. No free barrier is
assumed. The total across 8256 layers is exactly charged by 8256*N*2^20.

After sorting, the controller schedules adjacent-key comparisons in index order.
For each adjacent pair the root requests slot i and receives its record, then
requests slot i+1 and receives its record: four routed packets. Allow four more
full-packet equivalents for acknowledgments/completion, each path at most 128
hops. Thus packet work is at most 8*128*2^18=2^28 per adjacent pair; comparison,
address construction, counters and copying are covered by a total 2^29 per slot.
Processing these requests sequentially bounds queues and retained messages.
The root retains the first equal cross-tag pair while continuing the bounded
scan; no second global random access is required. The separate 2^26-per-slot
selection allowance is conservative slack for local flag/index handling and
final reporting. The root verifies the resulting scalar.

These resource bounds intentionally count aggregate work even if processors
are mostly idle. A single processor with an impossible wide address, an ideal
hash table, prebuilt data table, cost-free broadcast, or a free communication
oracle is not part of the construction.

## Input, output, fixed work and memory

The parser is a byte state machine for one z key and exactly twelve fixed-length
hex strings, allowing JSON whitespace at legal boundaries, checking all syntax,
duplicate/extra keys, bounds and end of document. It retains only its current
state, counters and twelve coefficients. Charge 4096 per consumed byte plus
4096 overhead, including input transfers, hex conversion and limb stores. On
malformed input stop with failure. L is the entire received serialization,
including leading, trailing and inter-token whitespace: no whitespace is free.
The canonical scored convention in proof.md has L=811; unlimited padding is
handled by the separate R(L) statement, not secretly discarded from its cost.

Fixed root work includes h computation, subgroup testing, index recovery and
reduction modulo r, final exponentiation g^x=z, canonical output, loading code,
small-field operations and failure handling. A bound of 2^50 covers all of it.
For the point-output variant the additional group arithmetic described in
output.md is included in the same fixed allowance; that document bounds it by
2^46. No organizer verification is substituted for attack-side work actually run.

Each record processor's code, incoming/outgoing ports, register state, local
record/scratch and counters fit in 4 MiB. Root state is below 16 MiB and has no
full table. Cube link packet storage is assigned to its endpoint buffers and
counted once. Peak total is below N*2^22+2^24 < 2^151 bytes, so the declared
2^152-byte bound is conservative. It includes all code copies, arithmetic
scratch, processor state and retained data. No precomputation is amortized.
