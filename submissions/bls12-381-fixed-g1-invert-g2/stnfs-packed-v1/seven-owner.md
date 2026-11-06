# Addressed thirteen-owner revision: complete resource supplement

This file supersedes the older linear-algebra layout and prices inherited from
276a1e92-22da-47f1-8bd2-027287da78de (commit91af88d731afa87d30cccd37b398ba88e309b7a4)
and the seven-owner idea of6734bd5a-26e4-42f4-b40e-31722d0d7bca. The latter schedule
was integrated in this session's429e9aa0-b313-439b-a37c-a86ff2df86a0 and rejected
because a256-byte slot times2^56 entries fills the64-bit address space. That
rejection is not evidence against the inherited smoothness model, nor a premise
of this revision. This revision explicitly replaces that layout and accounting.
All original NFS authors retain their provenance. The mailbox accounting also
builds on this session's ef18b93e-6a35-4d61-aaa5-906db3c7de66 and its source
f9b7a390-0898-4369-a97c-3d752133ef4d. No prior verdict is a premise.

## Machine interface, separately from H1-H5

Assume scalar64-bit instruction-equivalent local memories and error-free direct
indexed one-hop mailboxes, fixed reset loaders and blocking receive events with
no spinning, retries or spurious wakeups. This is an explicit hypothetical
instruction-accounting interface, not authenticated physical cluster availability.
Every transfer and all endpoint/queue/reset construction are charged below.
There is no free multi-hop routing tree: each ordered sender/receiver pair has
its own bounded mailbox slot, looked up by64-bit endpoint ID.

All LA backing storage needs fewer than2^34 nodes: sparse entries<2^68 records;
projection inputs<2^64 values; generator arrays have dimension k<=32, length
<2^57, at most8*k^2*length32-byte coefficients including Karatsuba scratch,
below2^75 bytes; retained sequences<2^70 bytes. Sparse-entry shards hold2^43 payload bytes, coefficient shards up to2^45,
row-descriptor shards up to2^45 and dense-row shards up to2^50. Every full
shard holds at least2^43 bytes; using that smaller capacity to upper-bound
the node count is conservative. Dedicated array tails and directories add
fewer than2^20 shards. Depth-first Karatsuba reuses scratch; simultaneously live array tails are
below8*k^2*57<2^20. Family payload totals are below2^76 bytes, so full shards
are fewer than2^33; adding the tail bound remains below2^34 nodes. These terms
and their headers fit below2^78 bytes, below the source main-phase2^83 allowance.
Use at most2^34 endpoints including the13 vector owners. Allow1024 bytes of
mailbox/queue/label state per ordered node pair: at most2^78 bytes globally and
2^44 bytes locally. Code/reset state at most64KiB per node is below2^50 bytes
in all. Together main storage remains below2^83 bytes; descent peak2^105 and
overall2^106 are unchanged. Initialize endpoint tables, labels, queues, code and
partitions in fewer than2^84 scalar instructions, far below the existing misc
row2^80 reference equivalents (more than2^94 instructions). This is a cold
construction, not free target-dependent advice. No future-target amortization.

For bulk streams allow8 scalar instructions per64-bit payload word across the
sender/receiver: load, mailbox store, receive load, destination store, two pointer
updates, loop control and handoff. Endpoint choice, page headers and blocking
flags use the per-entry/per-row/batch allowances below, not an extra free action.
One transaction at a time per link bounds queues. Time spent blocked is latency,
not executed instructions; there are no uncharged spin loops in this interface.

## Vector buffers and local addresses

Since N<=2^58.6841<13*2^55, assign global element j to owner j>>55 (0..12), with
local t=j mod2^55. A slot holds five four-limb field elements (160 bytes) in
192 bytes. Each owner holds BOTH current and next Krylov arrays, each192*2^55
=3*2^61 bytes; together3*2^62 bytes. Start them at2^20 and2^20+3*2^61. Their
last byte is below2^20+3*2^62, leaving more than2^61 bytes for code, bounded
buffers, endpoint tables and scratch. Use a working window above that endpoint;
no live object wraps a64-bit address. A row's next values are written by its
root owner into the next local array, after reducing all partial sums. Swap
array bases after each full matvec. Thus no in-place corruption or hidden
second full vector array is assumed. Current/next Horner vectors use separate
32-byte slots in reused phase-local arenas; they are not simultaneously placed
over either live Krylov array. X and retained sequences live in backing shards.

For stride192, replace base+(t<<8) by3*t followed by base+(3*t<<6). t<2^55 means
3*t<2^57, so both instructions and final addresses fit64 bits. This adds one
instruction per unit/general Krylov entry. The original Horner stride32 is
unchanged. The updated listing contains sixteen copies of the39-instruction
unit body, then subtract/branch:626 instructions per16 entries, or39.125/entry.
A saved group length at[sp] and four tbz decisions select8/4/2/1 tail copies.
[sp,#8] and saved registers at[sp,#16] remain untouched. Additions, carries,
entry order, overflow bounds and general arithmetic are preserved. The general
Krylov loop is35 instructions. Branches fit their AArch64 ranges.

## Variable rows, backing and bounded batches

301 bounds TOTAL sparse entries divided by N, not maximum row weight. Store
packed8-byte sparse records in a separate stream sharded by global STREAM-SLOT ordinal
(two64-bit limbs, below2^68). A shard contains2^40 entries, so its local offset
is below2^43; its endpoint ID is below2^28. A32-byte row descriptor stores a
128-bit start ordinal,64-bit length and reserved word. Row descriptors are
sharded by row>>40, local offset(row mod2^40)*32<2^45. Dense21-column rows use
1024-byte records and the same row shard split, offsets<2^50. A row can exceed
any local buffer: it is processed as successive group-relative batches of at
most2^16 records. Each batch uses512KiB. Groups are(row,vector-owner,sign,type).
There are at most26 unit groups and26 general groups per row; a row's counts
and partial sums require bounded scratch only. Indices and counters are explicit
multiword arithmetic, not arbitrary-width unit-cost addresses.

Charge8 extra instructions per packed sparse word at the8/word mailbox envelope.
The destination link stays fixed for each group-relative batch. Endpoint/header
lookup and its fixed handoff/control are charged to first-group or extra-batch
allowances, not repeated as an uncharged per-entry selection. Payload copies,
source/destination pointer updates, counted loop and word handoff fit8/word.
The first26 groups of each type are charged2048 per vector-row for calls,
register saves/restores, headers and control (the listings need well below100
instructions per unit call and20 per general call; five calls per general group
still fit this allowance). Additional batches beyond the first in a group are
at most total entries/2^16. Bound each extra batch's full setup/teardown/control
by2^16 instructions per vector, and charge ONE additional instruction per
vector-entry. Batches are group-relative, not aligned to arbitrary backing pages;
a page crossing alone never creates an extra kernel call.

Projection X uses separate16-byte records sharded by global value ordinal;
vector/sequence/generator coefficients use32-byte records, each at most2^40
values per shard (local offsets below2^45). Buffers are bounded and independent
of row length. Projection terms receive a conservative128 instructions instead
of60: six words of newly charged transport at8 each plus20 for reductions,
headers and addressing. Final owner reductions are O(N), absorbed by the extra
20*(m+n)*N^2 allowance at this full-production N. Generator work is doubled:
each original charged arithmetic unit costs at least97+420=517 instructions,
while reading two four-limb inputs and writing one output plus sharding/control
costs at most192 extra instructions. The extra factor covers that movement,
recursive scratch copies and fixed endpoint control. All generator arrays use
the above backing rather than fitting inside a vector owner's reserved window.

## Complete new prices and ledger

Thirteen owners require12 sends/combinations. Base row work, including the
expanded first-group allowance and unroll group controls, is
12*(15+96)+12*(62+2)+92+180+183+2+13*10+40+2048+32=4807.
The32 is ceil(26*6/5): new unit controls str/lsr/cbz/ldr/four tbz replace the
old cmp/branch, adding at most6 per group. Add2048 per vector-row for dense
transport, row descriptors, synchronization, bounded buffer control and final
local output. Two84-word dense operands at8/word use1344; descriptors and the
remaining704 cover headers/output. Total6855 is rounded to C_ROW6858.

The single-vector entry price is(20.9+1)+8+1=30.9. For five-vector Krylov,
0.87*(4+3+39.125)/5+0.13*(4+3+35)=13.48575; add8/5 for sparse streaming and1
for extra batches, giving16.08575, rounded to16.086. H4u remains0.87. Use m=26,
n=5, C_PROJ128, and the doubled generator bound108.389731954. All other upper
ledger rows remain unchanged. With N<=2^58.6841 the LA formula gives
117.68442829489047950448170183378054996638899602811391748498066471603059591458653,
charged117.684428298. Upper rows plus wrapper sum below119.4345700829; core
meter119.434570083 plus wrapper66 fits the exact final score119.434570084.

This raises the bound relative to the rejected seven-owner package because
missing address/streaming work is now paid. It improves the valid generic148
control by using the inherited conditional STNFS method. No false claim to
beat every existing public heuristic bound is made. All H1-H5 success premises
remain unchanged; the new prices bound the same logical computations. The
source0.44542 success minus wrapper failure2^-48 and its narrow >=0.39 bound
at10% extra descent work remain conditional on those premises and this interface.

store16.check now tests the actual thirteen-owner192-byte DOUBLE-buffer layout,
entry/row/dense shard boundaries and non-wrapping endpoint windows. The hosted
audit checks those boundaries, new row tail,626-instruction body, prices and
Decimal ledger. Full-run records are inherited, not a run of this revision.
No candidate code is executed on the host. The complete optional reproducer is
losslessly LZMA/base85 encoded with a standard-library loader. This is analytic
exploratory evidence, not production execution or platform acceptance.

## Record and partial-sum protocol details

Cold partitioning packs word=(coefficient<<55)|LOCAL55-bit column index.
The coefficient is an absolute unsigned8-bit value (the source merged bound<=255);
bit63 is zero. The lower55bits and upper8bits occupy disjoint fields. Each row
has26 fixed(owner,sign) groups with a16-byte header of unit/general lengths.
Headers occupy TWO8-byte stream slots, entries ONE. Row descriptors use128-bit
slot ordinals. At most52N header slots plus301N entries remain below2^68 slots.
Global columns recover as(owner<<55)+local. Repacking and bitfield construction
stay below2^90 cold instructions within misc80. Unit/general type and sign are
known from group metadata; no input property or coefficient bound is weakened.

Both Horner and Krylov index loads postincrement8. An AND with the encodable
55-one-bit AArch64 logical immediate extracts the index. General kernels extract
the coefficient by LSR55 BEFORE masking/reusing the raw index register; this
replaces the old byte load, and its old pointer argument remains conservatively
paid but unused. The added index mask costs one instruction per entry in every
kernel. Therefore Horner loops cost12/34 for unit/general, and Krylov unit bodies
cost39, general35. No deinterleaving or separate coefficient stream is assumed.

Header/descriptor/dense and partial-sum transfers remain charged exactly as in
the addressed predecessor. The decoding cost is included before reducing the
sparse transport charge from16 to8. The new bitfields do not change arithmetic,
row order, H4u, batching, double buffers, success or machine interface assumptions.
This exact record compression is this session's new contribution; the addressed
predecessor score119.451090414 is not an acceptance premise.

Reset FIVE-limb kernel accumulators to zero for every batch. Promote and merge
each batch result into SIX-limb per-row positive/negative accumulators before
processing another batch. Each type batch has at most2^16 entries; the initial
Horner call may combine one unit and one general batch in separate512KiB buffers,
so its sum is below2^(255+8+17)=2^280, within320bits. Additional calls handle one
type batch. An unrestricted row sum may exceed320bits but is below2^323 and fits
the384-bit row accumulators. This avoids carrying a narrow accumulator through
an arbitrarily long row. Reset/promotion/merge need at most20 instructions per
vector: three zero stores, six input/row paired-or-single loads, six carry adds,
three result stores and two pointer controls. The2048 first-group allowance covers
these operations too: Horner has at most26 initial combined calls, each below64
instructions of entry/exit, arguments and reset/merge; five-vector paths have
26 unit calls at below64/5 plus26 general calls at below12+6+20=38 per vector,
with the shared unit reset/merge20 also charged per vector, in all below2048.
Additional batches reset/promote inside the existing <=2^16-overhead/vector
allowance; the +1 per-vector-entry amortized charge therefore remains sufficient.

Each row-tree partial is a12-word positive/negative sum. Charge96 transfer
instructions at8/word, PLUS the inherited15-instruction send kernel (conservative
double counting), not the old unexplained18-instruction transfer. Twelve edges
therefore add936 to the previous row budget. The resulting row tail is6855,
charged6858. Combine's extra2 remains conservatively charged. The explicit
mailbox interface applies to these transfers too; no specialized free path is used.

The packed resource variant's current bound is119.434570084, improving the
addressed predecessor by0.016520330 bits (about1.14% of modeled total work).
store16.la_check now checks coefficient/index packing at boundary values; sparse
slot shards use2^40 words,2^43 bytes. With the smaller minimum shard size, family
payload<2^76 still gives fewer than2^33 full shards plus tails, below2^34 nodes.
