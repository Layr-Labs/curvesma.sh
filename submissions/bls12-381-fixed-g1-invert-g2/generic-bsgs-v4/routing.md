# Explicit storage routing supplement

This supplement specifies the abstract multi-processor implementation of the
storage mechanism in proof.md. The score counts total work on all processors,
not wall-clock time. Each processor has a local ordinary 64-bit address space.
No processor needs to map the entire table into one address space. One leaf
processor owns one 8192-byte record, with fixed code and scratch below 65536
bytes. Each array uses a binary routing tree with 2^128 leaves, height 128,
and fewer than 2^129 nodes. Internal nodes route messages according to the next
bit of the 128-bit record index. Node identities are three 64-bit words (a 128-bit prefix and a depth). They
never become single-machine virtual addresses.

Use explicit point-to-point word transfers, charged as loads, stores and
control rather than free communication. For a fixed-length record stream,
charge 8 scalar instruction equivalents per 64-bit word per edge: load, send
store, receive load, destination store, two address updates, decrement and
branch. Network delivery imposes no arithmetic primitive beyond copying the
word; latency is not the metric. Conservatively use 142 edges including request
headers and endpoints. A copy reads the source stream and writes the destination
stream, hence at most 2*1024*142*8=2326528 word-transfer instructions. Fixed
index/header routing, 128 decisions and bookkeeping take fewer than 2^18
additional instructions. Comparing two records reads two streams with that
same bound; scanning 4572 coefficient bits adds fewer than 2^18 instructions.
Thus each record comparison or copy is below 2^22 instructions. Arithmetic
operands use the same streams; the 2^30 multiplication bound has ample space
for its constant number of input/output transfers. Routing is deterministic;
there are no random failures or distributed hash-table assumptions.

Allocate at most 65536 bytes of code, local scratch and metadata per node,
including leaves, for each of two trees. This costs less than 2^146 bytes.
Consequently this supplement replaces the payload-only peak-storage bound
with a conservative 2^147-byte bound. This is a reported memory bound only;
it does not alter the scalar score. Code replication is counted in storage.
Construct and initialize every node's local state and routing labels using
at most 2^20 instructions per node, allowing code copying and zeroing, for
less than 2^150 instructions for both trees. Teardown is also below 2^150.
This fits the extra 2^160 allowance already in proof.md. Cold initialization
is paid in full; no prebuilt target-dependent tree or table is supplied as advice.
The total remains less than 2^162 instructions and 2^148 reference equivalents.

This is an analytic storage/transfer model in the challenge's instruction
accounting convention, not an assertion about an available physical network
or a production allocation. Total work, all processors and transfers are counted.
No execution of the exponential attack is claimed. The actual local preparation
requires no such hardware because the candidate is a permitted analytic argument.
