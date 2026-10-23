# Threshold-sized BSGS with padded batches and charged record routing

This is a classical analytic upper bound for the unchanged production target,
not a practical break or a production execution. The sole submitted score is
545/4 = 136.25: total work W is at most 2^136.25 target-group-operation
equivalents. Success is at least 39/100 over the full uniform challenge space.
Peak memory is below 2^136 bytes. Preprocessing is below 2^136.25 equivalents
and is already included in W. External nonuniform advice is zero bytes (the
max(1, advice_bytes) reporting convention gives log2=0). There are no random
walk, occupancy, independence, asymptotic extrapolation or success heuristics.
Instruction bounds and the finite distributed-memory realization are analytic
implementation obligations, not measured hardware costs.

## Scope, provenance and the actual change

Use exactly the protected p, prime r, generator P, input Q, canonical encoding
and acceptance predicate. C=1000 scalar AArch64 instruction equivalents per
target-group equivalent; every instruction below is counted as R, so H=0 and
W=R/1000. No invented group-operation price or 256-bit unit instruction is used.

Public submission 1add1202-af75-4f66-8add-3beee9ee738b, commit
90e5ac979ecba4aa783fbd8860133e6a826ee1df, supplies the interleaved affine
batching, full-key radix sorting and native-leaf routing construction adapted
here. Its recorded score is 136.62, but its acceptance is NOT a premise. The
present argument restates the needed method and bounds. The partial-coverage
framing was also informed by fb9b62cc-bf39-4bb9-a5b3-964f39f51a63. Its price
and implementation are not reused. The organizer baseline and standard BSGS
are additional provenance. This is parameter and accounting engineering, not
a newly discovered DLP algorithm or a claim to beat heuristic rho entries.

Relative to the full-coverage radix construction, choose the smallest square
grid meeting 39%, charge ALL padded batch slots, and increase routing, local
record handling, bucket initialization and router-state allowances. None of
the predecessor's captured tests are imported as validation of this candidate.

## Parameters and exact success argument

Let c=ceil(39*r/100), m=ceil(sqrt(c)), h=m, B=2^16, L=ceil(m/B), N=2m.
These are computable public parameters; their construction is charged below.
For the protected r:

    m = 212506270031153001307260662621714878127
    L = 3242588348864028950611277200648726
    100*m^2 - 39*r = 17906974590669998841605620727838004533757
    (m-1)^2 < c <= m^2 < r

Decode Q canonically, reject duplicate/extra keys and off-curve points, and
check [r]Q=O. Null denotes O. Build the baby records (key([j]P),j), 0<=j<m,
starting with (key(O),0). Separately compute D=-[m]P. Build giant records
(key(Q+[i]D),i), 0<=i<m. Stable-sort both streams by their FULL point keys;
merge them. For an equal key, form a=im+j, verify [a]P=Q, and emit a using the
profile's 64 lowercase hex scalar encoding. If the merge ends without a
match, return failure. The cost charges complete generation, sorts and scans
even if early termination could save work.

The full key is 0 for O and (2+(y mod 2))*2^256+x for a finite point. It has
258 significant bits and distinguishes all valid points. Each stream contains
distinct points: its increment is nonzero in a prime-order group and its
length is below r. An equal key implies x=im+j mod r. Since 0<=im+j<m^2<r,
the candidate is already canonical and unique. Conversely every x<m^2 has
that division into i,j and is found. Success is EXACTLY m^2/r, at least .39.
All remaining uniform targets fail and are included in the same work budget.
The input distribution is not narrowed to the covered prefix. Infinity is
recovered through i=j=0. No randomness, restarts or amplification are needed.

## Actual batch schedule and exceptional points

For each family, chain t, 0<=t<B, owns the labels tL+k<m, 0<=k<L.
Baby seeds are [tL]P; giant seeds are Q+[tL]D. Compute each independently by
complete binary scalar multiplication and, for giants, a complete addition.
All production starts tL are below m; for a general small test, empty chains
are masked. At layer k emit the current labeled state, then advance by P or
D. This produces exactly the logical streams above, although their original
order need not be increasing label order. Labels are retained, and sorting
does not require label order. The next layer depends on the completed prior
layer; this is NOT batching inversions along one serial dependent chain.

For ordinary addition U+S form d=x(S)-x(U), n=y(S)-y(U). Replace denominators
of inactive or exceptional lanes with 1. Prefix-product/reverse-product
batch inversion takes at most 3B multiplies and one inverse per layer. The
affine formulas slope=n/d, x'=slope^2-x(U)-x(S),
y'=slope*(x(U)-x')-y(U) use at most three more multiplies and eight modular
adds/subtracts. Charge every one of 2BL slots, including inactive slots, at
6M+8A+1024 instructions; charge all 2L product inversions separately. Ordinary
slot overhead includes local loads/stores beyond the kernels, masks, extraction
of denominator/numerator, state and multiword label updates, and loop control.

Exceptional current states are O,S,-S. Both progressions are injective, so
at most six exceptional updates exist across both families. Substituting 1
keeps inversion valid. Complete affine handling supplies the exceptional
output. Charge 16 complete operations in addition to the ordinary padded
charges. Complete seed operations are charged independently. Both families
perform their final unused updates; they are included.

## Four-limb arithmetic and normalization

All stored field elements are canonical four-64-bit-limb integers in [0,p).
Temporaries carry a fifth limb or an eight-limb product as appropriate. Modular
add/subtract A=256 covers six passes of at most five limbs, eight instructions
per limb, plus sixteen instructions: loads/stores, ADD/ADCS or SUB/SBCS,
carry/borrow, compare/correct, selection and local control. Results are reduced
canonically. Negation, equality and scalar limb operations fit such passes.

For M=4096, multiply using sixteen 64x64 schoolbook limb pairs. Each uses MUL
and UMULH; accumulate the two words and propagate carries across at most eight
limbs. Allow 128 instructions per pair: up to twelve per visited limb plus
32 for product, address/index/control work. This costs at most 2048.
Use p=2^256-c0, c0=2^32+977, to fold z -> low256(z)+c0*high256(z) three times.
Initially z<p^2. The first high part is at most four limbs; each later one is
smaller. A fold needs at most four low/high constant products, bounded carry
addition and transfers; 256 instructions per fold includes eight multiply
instructions, loads, carries, stores, fixed address and loop work.

After fold one z<2^290; after fold two z<2^256+2^67. After fold three, either
there was no carry and z<2^256, or the high part was 1 and z<2^67+c0.
Thus z<2^256+c0<2p. KEEP a fifth limb through compare/subtract/select; one
subtraction of p canonicalizes it. Reserve 256 for that operation and 256
for call setup, clearing, input/output pointers and transfers. The sum
2048+3*256+256+256=3328 fits M. Squaring uses the same bound. Neither modular
congruence nor limb repacking alone is treated as canonical normalization.

Inversion of nonzero z is binary powering z^(p-2). At most 512 multiply/square
calls, 64 overhead instructions per call and 256 setup give
I=512*(4096+64)+256=2130176. No inverse of zero is requested.

For seeds, exceptional cases and scalar verification deliberately retain the
looser bit-serial multiplication bound Mslow=256*(2A+64)+256=147712:
scan all 256 multiplier bits, conditionally add and double mod p, with fixed
bit/control work. Islow=512*(Mslow+64)+256=75661568. The complete affine
bound G=Islow+16*Mslow+32*A+4096=78037248 covers ordinary, doubling,
infinity, opposite-point and zero-y branches. Doubling uses 3x^2/(2y), with
y=0 returning O. Binary scalar multiplication uses at most 512G. We do not
charge both the reference group operation and any of these instructions.

## Full-key records, radix sorting and local counters

A 56-byte record has x in four words at offsets 0..31, a tag word at 32..39
(0,2,3), and a two-word index at 40..55. Infinity has x=0 and tag=0. Index
labels fit 128 bits. Lengths, exclusive endpoints, histogram counts, prefix
sums and scatter cursors use THREE words, as do logical record addresses.

Each stream is sorted by six stable LSD passes of 43 bits. K=2^43 buckets
use 64 bytes each: three-word count, three-word cursor, sixteen padding bytes.
This native-addressable 2^49-byte array is reused. Clear it; count digits;
scan buckets in order, recording exclusive prefix sums; scan source records
in increasing physical position, incrementing cursors and scattering each
record to the other array; swap arrays. All records may be in one bucket.
No occupancy or hashing assumption is made; a stable pass preserves earlier
digits, and six passes sort all 258 bits. Each destination is written before
it is read. Empty capacity padding is never read or cleared as payload.

Per record per pass: THREE routed accesses (count read, scatter read, scatter
write), plus 512 local instructions. Digit extraction fits six operations;
three-limb counter loads, add/carry, stores take nine; metadata addressing and
loop/control fit sixteen; remaining allowance covers split key digits, bounds,
cursor use and bookkeeping for both phases. Multiword route formation and
payload copying are inside routed-access cost F. Per bucket per pass allow
192 instructions for clearing eight words, prefix loading/copying, three-limb
cumulative addition, bounds checks and loop work. All-one-bucket counts stay
below 2^128 here and are safe in three limbs. No arbitrary-width counter costs
one instruction.

## Storage addressing and work at both transfer endpoints

Use four record arrays, two per stream, each with logical capacity below 2^128.
Split into leaves holding at most 2^58 records. A leaf's 56*2^58 bytes plus
local state fits below 2^64; offsets use native 64-bit addresses. Four arrays
need at most 2^72 leaves. A binary routing tree of depth 72 suffices, with
fewer than 2^73 nodes. Each node is a scalar worker with 8192 bytes of code,
mailboxes and state, separate from payload. This is a finite theoretical
distributed machine, not available physical hardware. Ordinary neighbor links
use bounded mailbox transfers. All endpoint/router work is summed. No processor
spins while idle; one outstanding request uses event-driven handlers. Latency
is irrelevant to this computation metric, not an uncharged arithmetic resource.

A read sends four words downward and eleven upward; a write reverses those
packet sizes. Four scalar operations per transferred word cover sender load,
mailbox write, mailbox read and receiver store; packets are fixed-size unrolled
copies, so no per-word loop is omitted. Sixteen additional instructions per
direction per edge cover path-bit selection, event/acknowledgement handling,
dispatch and control. A whole access costs at most

    ((4+11)*4 + 2*16)*72 + 512 = 7136 < F=8192.

The 512 endpoint allowance covers three-limb address/index composition,
56*index via shifts/subtraction with carries, request packing, leaf byte
addressing, seven-word payload copies, bounds and return handling. Neighbor
mailbox access is local; remote indices are routed, never native pointers.
Control words and in-flight payloads fit the node-state allowance. All traffic
from the engine to leaves is included in F. Work to distribute code and create
all nodes is bounded by 8192 instructions per node, total 2^86: up to 1024
words per node copied at four instructions per word plus equally much reserve
for path/link setup, clears and traversal. Payload is not initialized separately
because it is entirely overwritten before its first read.

The engine holds both B-chain families and batch scratch in native-addressable
memory. B*1024 bytes plus the histogram and 2^28 bytes of additional engine
code, control, scratch and I/O suffice. For a simple conservative total below,
the 2^28 engine reserve INCLUDES B*1024 (2^26) bytes of chain scratch. All
router and leaf programs/state are counted, as are both source/destination
buffers. This abstraction follows the protected allowance for aggregate
parallel work and memory; realizing its astronomical resources is not claimed.

## Merge and full cold-work ledger

Each generated record costs F+1024 for its write and canonical encoding,
tag/index construction and local bookkeeping. A two-head sorted merge reads
at most N records and charges F+1024 per record for transfer, five-word key
comparison, endpoint/cursor updates and branches. Equality tests all words;
it does not infer full equality from the final high-limb subtraction flag.
Streams have no repeated keys, so there is at most one match. At most one
candidate needs a scalar verification. The 4096G reserve pays scalar product
im+j, range checks, verification and output work, even though fewer than
512G are needed for the group check alone.

All rows below are scalar instruction equivalents. Constants are
M=4096, A=256, I=2130176, G=78037248, F=8192, K=2^43.

| Work | Bound |
| --- | ---: |
| All padded affine slots | 2BL*(6M+8A+1024) |
| Every batch inversion | 2L*I |
| Complete seeds and D | B*1025G+512G |
| Initial record writes | N*(F+1024) |
| Both six-pass sorts | N*6*(3F+512) |
| Bucket clearing and prefixes | 12K*192 |
| Complete merge including misses | N*(F+1024) |
| Complete exception reserve | 16G |
| Candidate reconstruction/verification | 4096G |
| Routing initialization, code/state | 2^73*8192 |
| Input, parameter calculation, setup/output | 2^44 |

Parameter calculation (integer square root, ceiling divisions), input validation,
[r]Q, negation of D, fixed setup, all short loops and final serialization fit
2^44; scalar validation costs below 512G and fixed 256-bit operations are
tiny relative to this reserve. Invalid finite encodings are rejected before
field/curve operations. There is no unbounded parser allocation requirement:
bounded input for the exact profile is checked before parsing, duplicate keys
and noncanonical forms are rejected. On valid uniform challenges this is a
fixed cost. No research/design effort or organizer challenge generation is
charged, as those are the explicit protected exclusions.

The sum is

    T = 83574680044327118202874855850367101169472512.

The claim follows from the exact integer inequality T^4 < 1000^4 * 2^545.
For numerical orientation only, log2(T/1000) is about 135.94019. The submitted
136.25 is a conservative upper bound, not this floating-point diagnostic.
All preprocessing is a subset of T, so preprocessing_log2=136.25 is a valid
included ceiling. No future-target amortization is used. Peak bytes are at most

    4m*56 + 2^73*8192 + K*64 + 2^28
    = 47601404486978349664078844326481535752480 < 2^136.

No advice bytes are supplied; runtime tables are paid construction and storage.
The deterministic failure paths cost no more than the complete loops above.

## Evidence boundaries and primary sources

experiment.py is a bounded specification check for the shared credential-free
runner. It checks the literal production ledger and exact threshold inequality;
exhausts a small elliptic curve with a partial grid; exercises masked batch
lanes, exceptions, stable sorting and full-key merge; checks small grids over
the real secp256k1 field on planted in/out-of-range targets; and checks field
fold/carry boundaries. These are actual mechanism tests when executed, not an
empirical estimate of production success or a native instruction measurement.
The production success proof is counting; the cost proof is the analytic ledger.
Python big integers and small radix buckets in tests are not given unit-cost
status in the attack. No production search, physical routing implementation
or compiled AArch64 trace is supplied.

Primary references, for provenance only: SEC 2, section 2.4.1,
https://www.secg.org/sec2-v2.pdf ; Handbook of Applied Cryptography, section
3.6, https://cacr.uwaterloo.ca/hac/about/chap3.pdf . The essential algebra,
scope and prices are stated above; the reviewer need not retrieve a URL.
Exploratory acceptance, if obtained, means plausible_not_refuted, not proof
certification, human promotion or a practical weakening of secp256k1.
