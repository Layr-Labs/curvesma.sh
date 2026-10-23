# Partial BSGS with explicit distributed table accounting

This is a conservative analytic classical cold attack, not a practical recovery,
a new discrete-log algorithm, a best-known estimate, or a qualified result.
The submitted upper bound is 2^146 target-group-operation equivalents; the
score is exactly the rational 146/1. The actual derived upper bound below is
smaller. Preprocessing is also bounded by 2^146, included once in total work.
Peak aggregate storage is bounded by 2^143 bytes. There is zero external advice.
Success on the unchanged full-production uniform challenge is strictly greater
than 1/2. No heuristic premise, randomness, retries or amplification is used.

## Target and algorithm

Use exactly p,r,P and the canonical encodings in the protected secp256k1
profile. In particular r is the specified 256-bit prime order, not a toy
order, and Q=[x]P is the only challenge-dependent input. Let m=2^128 and
h=2^127. Decode Q, check the curve equation, and check [r]Q=O. Handle infinity
as a tagged point with zero coordinate fields. All arithmetic uses 64-bit limbs.

Starting at O, store (jP,j) for j=0,...,m-1 and perform m repeated complete
affine additions by P, retaining mP. Set D=-(mP). Sort the m records by the
full canonical key (infinity tag,x,y), using bottom-up merge sort and two
buffers. Start gamma=Q. For i=0,...,h-1, binary-search the table for gamma.
On a match j, form a=im+j, verify [a]P=Q, and output a as 64 lowercase hex
characters. On a miss update gamma=gamma+D. Return failure after h misses.
Charge even the final giant update and the full loops when early return saves
work. No scalar reduction is necessary because im+j<mh<r.

The records are distinct since m<r. A hit implies Q=[im+j]P and hence x=im+j,
since both scalars are canonical and below r. Every x in [0,mh) has exactly
one decomposition im+j with these bounds, and is found. No x outside that
interval can hit. Therefore success is exactly 2^255/r>1/2>0.39. This does
not condition the input distribution on that interval: the attack fails on
its complement. The deterministic algorithm also satisfies the fresh-coins
probability space. Attack-side verification prevents an incorrect output.

## Field and point arithmetic in scalar instructions

Field elements use four 64-bit limbs, with extra carry/scratch limbs. Bounds
below include operand loads/stores, addressing and local control, and use no
SIMD or quantum instruction. They are deliberately loose analytic instruction
bounds, not measurements of compiled code. Modular add/subtract A=256 covers
six passes over at most five limbs at eight instructions per limb plus sixteen
instructions of control. Extend, add/subtract with carry/borrow, compare, and
conditionally correct once. Negation and equality use no more than A.

Field multiplication M=8192: form the eight-limb product using sixteen
64x64 partial products. Each uses MUL and UMULH and carry propagation across
at most eight result limbs. A budget of 128 instructions per pair (2048 total)
includes loads/stores, address arithmetic, clearing and loop overhead. Reduce
using p=2^256-c, c=2^32+977. If t=L+2^256 H, replace it by L+cH. Initially
H has four limbs; multiplication by the single-limb c and addition needs at
most ten limb products/additions including carries. After this first fold,
t<2^256*(c+1), so H<=c. The second fold is below 2^256+c^2+c; the third
has H<=1; if the third overflows, the fourth has H<=1 and its result is below
2^67. Thus four fixed folds suffice, followed by at most one subtraction of
p from a value below 2^256. Allocate 1024 instructions per fold, plus 2048
for final correction, operand preparation and result storage. Total <=8192.
Each limb-product carry chain is bounded in length; no arbitrary precision
operation is charged as one instruction. Squaring uses this same bound.

For a nonzero z compute z^(p-2) by binary square-and-multiply. At most 512
field multiplications and 65536 other instructions suffice. Thus
I=512*M+65536=4259840. Fermat's theorem supplies the inverse. Standard affine
chord-and-tangent formulas, including doubling, use one inverse and at most
16 multiplications plus 32 modular additions/subtractions and 4096 local
instructions. Explicit branches handle infinity, equal points, opposites and
y=0; inversion of zero is never requested. Hence every point update is
bounded by G=2^23 instructions, since I+16*M+32*A+4096=4403200<G.
For distinct finite x: s=(v-y)/(u-x), x'=s^2-x-u, y'=s*(x-x')-y. For doubling
s=3x^2/(2y). These formulas justify both the algorithm and its exceptional
cases. No reference group operation is charged in addition to these instructions.

## Storage, addresses and transfers

An exponential flat array cannot be addressed by one 64-bit load. Use a finite
binary routing tree of depth 128 for each of the two buffers, with one leaf
per 128-byte record. A record stores an eight-byte infinity tag, four limbs
for each coordinate and four limbs for j (104 bytes), padded to 128 bytes.
Each tree has fewer than 2^129 nodes. Each node is a distinct scalar worker
with at most 4096 bytes for code, node state, record/payload, links and scratch.
Workers have ordinary local 64-bit addresses. A logical index uses two limbs;
scalar reconstruction uses four. Loop endpoints, counts and merge widths can
equal m=2^128 and therefore use three 64-bit limbs (at least 129 bits),
including the final width and half-open endpoint. Actual leaf indices remain
128-bit values. Bounds checks prevent a fetch at the endpoint m. This
three-limb control arithmetic is included in the per-access and per-pass
allowances. Links connect parent and child workers.
The attack pays creation/initialization of all workers and all retained state.
Unused workers do not perform background computation.

A fetch or store travels at most 128 links down and 128 back. Packets contain
at most 32 64-bit words, including payload, logical index and message metadata.
At each directed link, at most 256 instructions pay copying with loads/stores,
address-bit selection, local addressing, branches and protocol bookkeeping.
This is achievable by a fixed bounded packet loop (64 payload load/stores,
at most 128 loop/control instructions and 64 for headers/routing). Charge
these instructions at BOTH communication endpoints by splitting the same
per-link budget: the packet copy uses <=128 instructions at each endpoint.
An alternative implementation may use a shared mailbox on each link; its
reads/writes are included, never a free remote array access. R=2^17 instructions
per logical record access exceeds 256*256=65536, leaving another 65536 for
request creation, replies, local record work and synchronization. Transfer
latency is irrelevant because the metric sums work rather than wall time.
This model assumes arbitrarily many finite linked workers, as permitted by
the parallel-work and aggregate-memory convention; it assumes no free global
addressing, oracle, simultaneous broadcast or physical wall-time bound.

Build the two trees recursively and initialize nodes. At most 2^130 nodes,
with <=2^16 instructions per node to assign a path, establish three links,
initialize local storage and copy the finite worker program. The total
initialization/input/verification/output reserve F=2^148 instructions covers
this (<2^146), fixed decoding, curve/order checks, <=1024 point operations
for scalar multiplications, scalar reconstruction, serialization and shutdown.
No background per-worker idle polling is used: messages schedule bounded work.
Link buffers are included in node storage. The algorithm needs no random bytes.

Aggregate node storage is <2^142 bytes. An additional <2^142 bytes covers the
controller, input/output, arithmetic code/scratch, merge state, tree construction
state and bounded communication scheduling. Peak memory is therefore <2^143
bytes. All buffers and processor code/state are charged. The limit is
astronomical, but finite; this is an analytic bound, not available hardware.
There are no searched embedded constants or other nonuniform advice. Public
profile constants and algorithm design effort are the permitted exclusions.

## Complete work inequality

Let R=131072, G=8388608, F=2^148, m=2^128, h=2^127. For each merge pass,
there are <=m comparisons and <=2m record moves. Fetch both operands for
each comparison (<=2m accesses) and fetch/store for each move (<=4m accesses).
Allow another 2Rm for comparisons, loop work, multiword counters and messages.
The pass bound is 8Rm. There are exactly 128 bottom-up merge passes; partial
runs cannot increase this bound. Sorting cost S=8*R*128*m=2^155 includes
all record transfers and control. Binary search uses at most 129 comparisons
and fetches. Each comparison, branch and multiword midpoint calculation fits
1024 local instructions in addition to its R fetch charge. Bounds are:

- initialization, input, final verification/output and shutdown: F;
- baby generation and storage: m*(G+2R+1024);
- merge sort: S;
- giant traversal and search: h*(G+129*(R+1024)+R+1024).

The extra baby access pays any initialization/store handoff. The extra giant
access pays matched record retrieval; the fixed F pays its final scalar check.
There is at most one hit and one check, since keys are exact points and a hit
is always correct. Exhaustion costs no more than the same fixed bounds.
Add these blocks to form the exact instruction upper bound T. In particular,
baby coefficient <2^24 and giant coefficient <2^25, giving
T < 2^148 + 2^152 + 2^155 + 2^152 < 1000*2^146.
All preprocessing (initialization, babies and sorting) is included and is
bounded by the same 2^146 target-equivalent upper bound. Dividing once by
protected C=1000 yields W=T/1000<2^146. The sole submitted score is 146,
not a confidence, measured time, tightened logarithm or altered reference price.
Memory is supplementary, with memory_log2_bytes=143; preprocessing_log2=146;
nonuniform_advice_log2_bytes=0 denotes the max(1,0)-byte reporting convention.

## Evidence, uncertainty and provenance

This proof is self-contained. The accompanying research.md records primary
sources for curve parameters and standard BSGS. No unfamiliar external theorem
or another submission's acceptance is used as a premise. The construction is
based on textbook BSGS and the organizer's BSGS control. Public Yukon submission
1d657384-384f-4941-87e9-4900ace5265e materially informed the need for a separate
64-bit arithmetic ledger and partial coverage; it is declared as a dependency.
Its signed-window method and claimed bound are not assumed. We do not reuse
any claimed score, evidence result, code or validity as proof. Listing other
public scores supplied context only, not algorithmic premises.

No experiment or production execution is claimed. The declared experiment
plan is empty; no executable candidate file is provided. Numerical inequalities
are independently reproducible from literal protected constants and the
formulas above. The analytic arithmetic and transfer bounds, and the distributed
machine abstraction, require exploratory review. This submission does not claim
to improve the public frontier. No heuristic, random-walk mixing assumption,
GLV speedup or free preprocessing underlies it. Acceptance, if obtained, means
plausible_not_refuted under the selected policy, not proof certification.
