# Signed threshold BSGS with charged routing

This is an analytic classical upper bound for the unchanged production secp256k1
ECDLP. The submitted computation bound is **135.75 = 543/4** in log2 of
target-group-operation equivalents. Success is at least 39/100 over the full
uniform challenge distribution. Peak storage is below 2^135 bytes. Cold
preprocessing is bounded by 2^135.75 equivalents and is included in total work.
External nonuniform advice is zero bytes, reported as log2(max(1,0))=0.

There are no probabilistic heuristics. The success argument is exact counting.
The implementation costs are analytical instruction allowances, not measured
AArch64 traces. The storage is theoretical and astronomical. No production
attack, practical break, novel DLP algorithm, formal certification, or human
approval is claimed.

## Protected target and provenance

Use only the protected point Q=[x]P, with x uniform in 0..r-1, and the public
secp256k1 constants in target-profile.json. In particular,

```
p = 0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f
r = 0xfffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141
```

Use the protected generator, equation y^2=x^3+7, canonical JSON encoding and
acceptance predicate. No auxiliary powers, restricted-input promise, quantum
step, external oracle, or future-target amortization is supplied. The fixed
normalization is C=1000 scalar AArch64 instruction equivalents per reference
unit. All work below is R, with H=0 and W=R/1000. A complete affine addition is
not assigned the organizer's batched-reference price.

Material sources are public submissions
161020a8-da21-4676-a3c6-41c89fa4ee2a (threshold-sized BSGS, source commit
8369fa05652f13e581a4c6d4107033437c7a5b35) and
1add1202-af75-4f66-8add-3beee9ee738b (streamed, batched BSGS, source commit
90e5ac979ecba4aa783fbd8860133e6a826ee1df). Their public notes were read. The
former supplies the threshold framing, padded schedule and conservative
ledger; the latter supplies the underlying arithmetic, radix and routing
design, which was also read directly. Both are declared dependencies. Their
reviews and scores are provenance, never mathematical premises. The needed
construction is restated here; no external file is required to check it.

The substantive changes are signed windows, x-orbit keys with separately
retained parity, exact sign recovery, explicit treatment of repeated giant
orbit keys, and smaller leaves fitting 52-bit native address spaces. The
square threshold grid is replaced by the least m with m(2m-1)>=ceil(.39r).
All batch padding and unsuccessful scans remain charged. Negation-map BSGS
and simultaneous inversion are established techniques; the contribution is
their concrete integration and accounting under these challenge rules.

## Parameters, algorithm and success

Set c=ceil(39r/100). Choose the smallest positive integer m with
m(2m-1)>=c, and put n=2m-1, h=m, B=2^16, L=ceil(m/B), N=2m. Computation
of these public parameters is charged. Exact production values are:

```
c = 45158914802553316215192684153388284062606650068839212709216013625192082982792
m = 150264624583687887066834814587324391195
n = 300529249167375774133669629174648782389
L = 2292856210078245347089154275319281
m*n = 45158914802553316215192684153388284062614466497944190435557297090567462664855
100*m*n - 39*r = 781642910497772634128346537537968206357
```

The inequalities (m-1)(2m-3)<c<=m(2m-1)<r and m<r/2 hold. Decode Q
canonically and validate the curve and subgroup. For infinity use the explicit
point tag; do not assign it an ordinary x-coordinate.

1. Construct baby records for [j]P, 0<=j<m. Each retains x, infinity/finite
   status, y parity, and j.
2. Compute D=-[n]P and C0=Q-[(m-1)]P. Construct giant records for
   Gamma_i=C0+[i]D=Q-[in+(m-1)]P, 0<=i<m, retaining the same data and i.
3. Sort each stream by the orbit key: 0 for infinity, and 2^256+x for a
   finite point. Parity is payload, not part of the sort or merge key.
4. Merge sorted streams with two heads. Advance the smaller-key head. At an
   equal key, let d=in+(m-1). If the baby is infinity, set a=d. Otherwise
   set a=d+j when the two y parities agree, and a=d-j when they differ.
   Verify 0<=a<r and [a]P=Q, then return its 64 lowercase hexadecimal digits.
   Exhaustion without equality returns failure. Verification failure also
   returns failure; it does not trigger an uncharged search or retry.

For a short Weierstrass curve over an odd prime, two finite points with the
same x-coordinate have equal or opposite y-coordinates. The protected group
has odd prime order, so a finite group point never has y=0. Negation therefore
flips y parity, and the sign test is exact. Infinity corresponds to j=0.

For giant index i, the recoverable scalar values are precisely

```
in+(m-1) + {-(m-1), ..., -1, 0, 1, ..., m-1}
    = [in, (i+1)n-1].
```

These disjoint windows tile [0,mn-1]. Each candidate a is already in [0,r),
so a match cannot be an unintended wraparound. Conversely, every x in this
prefix gives exactly one index i and one signed baby offset. Thus success is
EXACTLY mn/r >= .39 on uniform full-production challenges. Targets outside
the prefix fail and are charged the complete budget. The deterministic
algorithm does not need fresh coins, amplification, or restarts.

Baby orbit keys are unique: jP=+/-kP in this range implies j=k; j+k<r
excludes a nontrivial negative equality. Actual giant points are distinct,
since n is nonzero modulo prime r and m<r. Giant ORBIT keys may repeat when
two points are negatives. This is allowed. The usual sorted merge reaches an
equality if one exists even with repeated giant keys; it returns immediately
at that equality. Disjoint scalar windows imply at most one cross-stream
match and therefore at most one candidate verification. No claim that the
giant orbit keys are all unique is used.

## Dependency-respecting affine batching

Generate each family with B independent chains. Chain t owns labels tL+k<m,
0<=k<L. Baby seeds are [tL]P. Giant seeds are C0+[tL]D. Each seed is computed
independently with complete binary scalar multiplication and, for giants,
one complete addition. At layer k, emit every active labeled state, then
advance it by P or D. Mask empty or exhausted chains. This gives exactly
the logical sets above, in a different enumeration order. Labels are stored.
Sorting is independent of the original enumeration order.

At an ordinary addition U+S, form denominator d=x(S)-x(U) and numerator
v=y(S)-y(U). Substitute denominator 1 for inactive or exceptional lanes.
For the whole layer form prefix products, invert the total product once,
then traverse backwards to recover each inverse. This costs at most 3B
multiplications and one inversion. For each ordinary lane compute
s=v/d, x'=s^2-x(U)-x(S), y'=s*(x(U)-x')-y(U), at most three further
multiplications and eight modular additions/subtractions. Outputs are
canonical because every field kernel is canonical.

Charge every one of the 2BL slots, including inactive padding and final
unused updates, at 6M+8A+1024 instructions. Separately charge all 2L
inversions. The 1024 includes lane masks, local scratch transfers beyond
kernel transfers, numerator/denominator preparation, loop control, state
and multiword label updates. At most six exceptional ACTIVE updates occur:
for each family the injective progression can visit O, S and -S at most
once each. Complete affine handling supplies those outputs after the dummy
denominator protects the batch. Charge 16 complete additions in addition to
the ordinary padded cost. Seeds are separately charged, not hidden inside
this exception bound. Total padded slots beyond the 2m emitted records are
16842. No inversion is batched across a serial dependent chain.

## Canonical four-limb arithmetic

The field representation is four 64-bit limbs for canonical integers in
[0,p). Carry temporaries retain a fifth limb; unreduced products retain eight.
Every operation includes loads, stores, addressing, carry/borrow and control.

Modular add/subtract allowance A=256 permits six passes over at most five
limbs, eight scalar instructions per limb, plus sixteen fixed instructions.
Add into five limbs, compare/subtract p and select the canonical result;
for subtraction correct a borrow by adding p. No high carry is truncated.
Equality, negation and the fixed-size scalar operations have similarly
bounded passes. All dimensions here are fixed small limb counts.

Use M=4096 for field multiplication. Sixteen 64x64 schoolbook limb products
use MUL and UMULH, accumulating their two words into the eight-word product.
For each pair allow 128 instructions: carry propagation through at most
eight limbs at up to twelve instructions each, plus 32 for product,
indexing, addresses and control. This contributes at most 2048.

Write c0=2^32+977, so p=2^256-c0. Fold the unreduced product three times
using z := low256(z)+c0*high256(z). Each fold uses at most four low/high
constant products plus bounded carries and transfers; 256 instructions per
fold suffice. Initially z<p^2. After the first fold z<2^290; after the second
z<2^256+2^67. In the third fold either the high part is zero or it is one,
in which case the result is below 2^67+c0. Consequently the result is below
2^256+c0<2p. Retain the fifth limb through comparison/subtraction/selection;
one subtraction of p makes the output canonical. Reserve 256 for that step
and another 256 for initialization, clearing, calling and input/output.
Total 2048+3*256+256+256=3328 fits M. Squaring uses the same bound. This is
an explicit canonical construction, not an assumption about the organizer's
179-instruction calibration or an uncharged five-to-four-limb conversion.

Nonzero inversion by binary powering to p-2 needs at most 512 multiply/square
calls, 64 extra instructions per call and 256 setup, giving
I=512*(4096+64)+256=2130176. Inversion of zero is never requested.

For seeds, exceptions and verification retain conservative slow kernels:
Mslow=256*(2A+64)+256=147712 via bit-serial conditional-add-and-double;
Islow=512*(Mslow+64)+256=75661568; and
G=Islow+16*Mslow+32A+4096=78037248 for complete affine addition/doubling.
Handle infinity, opposite points, doubling slope 3x^2/(2y), and y=0 before
any inverse. These bounds include canonical corrections and control.
Binary scalar multiplication uses at most 512G. These instruction costs
are divided by 1000, never by G. Reference operations are not double charged.

## Records, stable sorting and counters

Each 56-byte record uses four x words, one tag word (0 for infinity, otherwise
2|(y&1)), and two words of index. For infinity x=0. The sort key uses the
finite bit tag>>1, not the parity bit. Thus it has 257 significant bits,
although six stable LSD radix passes of 43 bits cover a 258-bit envelope.
The unused highest bit is zero. Retaining parity in the payload costs no
additional record word. Full-key equality tests every coordinate word and
the finite flag; a high-limb subtraction's zero flag is insufficient.

Let K=2^43. One native-addressable histogram of K 64-byte cells is reused.
Each cell holds a three-word count, three-word cursor and padding. In each
pass: clear every cell; read records and count the selected digits; scan
buckets to compute exclusive prefixes; read the source in physical order
and scatter records to the next bucket cursor, incrementing it; swap buffers.
This is stable even if every record falls into the same bucket. Six passes
produce the orbit-key ordering without a hash/occupancy assumption.

Every pass uses THREE routed record accesses per input record: count-read,
scatter-read and scatter-write. In addition allow 512 local instructions per
record for digit extraction including masking parity, three-word counter
increments, metadata access, endpoint checks and loop control. Per bucket
per pass allow 192 instructions for clearing eight words, prefix calculation,
cursor writes and loop work. All lengths, endpoints, counts, cursors and
logical addresses use three words. They do not overflow at 2^64 or 2^128.
Indexes retained in individual records need only two words. Every payload
destination is written before reading; untouched capacity is not payload
and needs no separate zeroing. Histogram clearing is explicitly charged.

## Addressing and transfers on a finite distributed machine

The storage cannot fit on one native-addressed machine. Use four record
arrays, two buffers per stream, each of capacity below 2^128 records. Split
each array into leaves with at most 2^46 records. A full leaf's
56*2^46 bytes plus local code/state is below 2^52. Four arrays need at most
2^84 leaves, served by a binary routing tree of depth 84 with fewer than
2^85 total nodes. Each node has 8192 bytes of code, state and mailboxes.
Use workers capable of 52-bit physical addressing with the MMU disabled;
no omitted OS or page tables are required. The histogram occupies 2^49
bytes on the engine and also fits. Instructions and registers remain scalar
64-bit AArch64 throughout. This smaller-leaf construction avoids treating
an entire 64-bit virtual space as necessarily implementable.

One request at a time is routed from the engine through neighbor mailboxes.
A read transmits a four-word request downward and an eleven-word reply
upward; a write uses the reverse sizes. Four instructions per transferred
word cover source load, mailbox store, mailbox load and destination store.
Fixed-size copies can be unrolled. Sixteen extra instructions per direction
per edge cover next-path-bit selection, dispatch, event/acknowledgement and
control. Add 512 instructions for engine/leaf endpoint packing, three-word
logical address formation, native offset arithmetic and record movement.
The resulting whole-access cost is

```
((4+11)*4 + 2*16)*84 + 512 = 8240 < F=8448.
```

All work at both ends and every intermediate node is summed. No worker
spins while idle: event-driven handlers are inactive until work arrives.
An event/wakeup and its acknowledgment are included in the control allowance.
Local physical mailboxes implement bounded neighbor transfers; no arbitrary
remote pointer or free global lookup is assumed. Transit state fits each
node's reserve. Communication latency and physical distance are not this
instruction-equivalent metric; every scalar data move is charged.

Charge 8192 instructions per node to distribute code, initialize state and
establish links, hence at most 2^85*8192. Recursive distribution sends each
child its program across one edge (not from the root for every child): up to
1024 words at four instructions each, with another 4096 for control, clearing
and construction. Code is replicated and stored in the memory bound. A
single control path of bounded depth suffices. No searched advice is placed
in code. The engine reserve of 2^28 bytes INCLUDES the two families of
batch state and scratch, with at most 1024 bytes per lane in aggregate,
plus program, constants, input/output and control. Histogram storage is
separate. The enormous hardware/storage assumption is explicit; physical
construction of this machine is not asserted by an instruction-work bound.

## Complete cold-work ledger

Record creation costs F+1024 per emitted record for its routed write,
canonical-key/tag formation, index and loop work. The final two-head merge
reads at most N records at F+1024 each, including comparisons, duplicate-key
handling, advancing cursors and failure detection. Finding a match ends it;
charging the full scan bounds every outcome. The 4096G extraction reserve
pays formation of in+(m-1)+/-j, range/equality checks, one [a]P verification
and output. Sign selection is not a second scalar multiplication.

The seeds cost B*1025G: B baby scalars at 512G and B giant scalars plus
addition at 513G. An extra 1536G covers D=-[n]P and C0=Q-[(m-1)]P,
their complete addition and negations with slack. All costs are scalar
instruction equivalents. Constants are M=4096, A=256, I=2130176,
G=78037248, F=8448, K=2^43.

| Work category | Instruction upper bound |
| --- | ---: |
| Every padded affine slot | 2BL(6M+8A+1024) |
| Every batch product inversion | 2LI |
| Complete seeds, D and C0 | B*1025G+1536G |
| Initial record writes | N(F+1024) |
| Both six-pass radix sorts | 6N(3F+512) |
| Histogram clears and prefixes | 12K*192 |
| Full merge, including failed targets | N(F+1024) |
| Complete exceptions | 16G |
| Recovery, verification and output | 4096G |
| Every routing node's code/state initialization | 2^85*8192 |
| Input validation, parameters and fixed work | 2^44 |

The fixed reserve covers bounded canonical input parsing, duplicate/extra-key
rejection, on-curve and [r]Q checks, parameter square root/ceiling calculations,
negations, and fixed control/output. A [r]Q check uses at most 512G, well
inside 2^44; other fixed computations use only small, bounded limb arrays.
Only valid profile inputs matter to the uniform-challenge probability, but
malformed encodings are rejected before unbounded allocation or field use.

Summing the rows gives

```
R = 60634932750577016349807926009615425840505344
R^4 < 1000^4 * 2^543.
```

Thus W=R/1000 < 2^(543/4). No floating-point rounding is needed. For
orientation log2(R/1000) is approximately 135.477273; 135.75 is the submitted
conservative bound, not a request for the judge to tighten it. Preprocessing
is a subset of R, so the included bound 2^135.75 covers it, including all
tables, clears, router construction and code distribution. There are no
randomness, restart, amplification or coordination terms left implicit.

Peak bytes are at most

```
4m*56 + 2^85*8192 + K*64 + 2^28
 = 33659275907062999353028056380885061285792 < 2^135.
```

This counts both buffers of both streams, every processor's program/state,
in-flight packets, retained table indexes/parity, the entire histogram,
chain scratch, constants, input and output. Runtime tables have charged
construction and are not free advice. External advice is zero. Only public
profile constants, algorithm design effort and organizer challenge generation
are excluded, as allowed by the rules.

## Evidence and limits

The declared standard-library Python experiment is a bounded specification
check for the shared credential-free Docker runner. Its exact production
integer arithmetic checks coverage, rounding, padding, work and memory
inequalities. Toy elliptic-curve enumeration checks both signs, infinity,
out-of-prefix failure, actual batched inversion, masked slots and repeated
giant orbit keys. Small planted secp256k1 cases exercise the same mechanism
over the real field. Field-fold/carry and stable-radix checks cover relevant
representation errors. Deliberately incorrect sign/key controls must fail.
Python big integers are used to check formulas, never priced as unit-cost
production instructions. Test timings are not attack costs.

Neither these finite checks nor source inspection measures the native
instruction schedule, constructs the distributed machine, proves every
resource allowance formally, or executes a uniform full-size attack. The
success proof and analytical cost argument, with their explicit implementation
assumptions, remain necessary. Exploratory review means plausible_not_refuted
only; it is distinct from a formal proof or practical cryptanalytic advance.

Primary references for context are Standards for Efficient Cryptography 2
(SEC 2), version 2, section 2.4.1, https://www.secg.org/sec2-v2.pdf; and
Galbraith, Wang and Zhang, Computing elliptic curve discrete logarithms with
improved baby-step giant-step algorithm (2017),
https://www.aimsciences.org/article/doi/10.3934/amc.2017038 . Negation-map
BSGS and simultaneous inversion predate this submission. The complete
argument above is supplied so the reviewer need not retrieve external URLs.
