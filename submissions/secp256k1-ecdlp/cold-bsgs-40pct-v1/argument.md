# Deterministic partial-coverage BSGS with full scalar accounting

This is an analytic accounting control for the protected secp256k1 target.
It is not a practical attack, a new discrete-log method, or a measured execution
on the production group. The algorithm uses ordinary baby-step/giant-step
(BSGS). It has no heuristic premise. Its large memory requirement is mandatory.

The submitted bounds are total work at most **2^146 target-group-operation
equivalents**, preprocessing at most **2^145** in the same units, peak memory
at most **2^146 bytes**, and success probability at least **2/5**. There is no
nonuniform advice and no random sampling by the attack. The fixed conversion
is **1,000 scalar AArch64 instruction equivalents per target-group-operation
equivalent**, as specified in the protected cost model. This document does not
assign a new reference price to point addition.

## Algorithm and exact success probability

Use the protected prime field modulus p, subgroup order r, and generator P.
Set

```text
m = 2^127
h = ceil(2*r / (5*m))
L = 127
```

In this instance, the exact values are:

```text
m = 170141183460469231731687303715884105728
h = 272225893536750770770699685945414569164
m*h = 46316835694926478169428394003475163141171880919487850230397683760192544571392
```

The input remains the full production challenge Q=[x]P, with x uniform in
0..r-1. No interval restriction is imposed on the challenge. The attack covers
an interval and reports failure on the other challenges.

```text
1. Decode and validate Q. Use an explicit representation for infinity.
2. Start U at infinity. For j=0..m-1:
     Store (U,j) in the baby table.
     Set U=U+P with complete handling of exceptional cases.
   Thus U=mP after the loop. Set D=-U.
3. Sort the baby table by its full point key, with bottom-up merge sort.
4. Set gamma=Q. For i=0..h-1:
     Binary-search the table for the full point key of gamma.
     If a record with index j is found:
       Set a=i*m+j.
       Check 0<=a<r and [a]P=Q.
       If the check passes, encode a as the required 32-byte scalar and return.
     Set gamma=gamma+D.
5. Return failure after the fixed loop ends.
```

The full point key is (infinity tag, x coordinate, y coordinate). Infinity has
tag 1 and zero coordinates. Finite points have tag 0 and canonical coordinates.
Coordinate comparisons use their unsigned integer values. Thus a key match is
point equality, without a hash-collision premise.

We have 0<m*h<r. For each x in 0..m*h-1, division gives unique i,j with
x=i*m+j, 0<=i<h, and 0<=j<m. At giant step i, gamma=Q-[i*m]P=[j]P, so the
table contains the needed record. Conversely, any match produces exactly such
an x, since i*m+j<r and P has order r. The final scalar multiplication verifies
the result. Infinity gives x=0 and is included.

The exact success probability is therefore m*h/r. The definition of h gives
m*h/r>=2/5>0.39. The upper bound m*h<2*r/5+m<r also holds for the protected r.
There are no attack coins, restarts, or amplification. The work bound charges
all h giant steps, including the full unsuccessful searches. Early success can
only reduce work.

## Local field arithmetic in 64-bit words

All arithmetic below uses scalar 64-bit instructions. A field value uses four
64-bit limbs. A temporary sum uses five limbs when needed for the carry.
An index uses two limbs. The scalar i*m+j uses four limbs. A logical byte
address uses three limbs. No 128-bit or 256-bit instruction has unit cost.

The following deliberately loose bounds include local operand loads, stores,
address calculation, calls, loop control, and scratch use. They apply to the
simple kernels specified here, rather than to a measured optimized library.

| Kernel | Instruction bound |
| --- | ---: |
| Canonical modular addition or subtraction | A = 256 |
| Field multiplication | M = 147712 |
| Inverse of a nonzero field element | I = 75661568 |
| Complete affine point addition or doubling | G = 78037248 |

For modular addition, load the two four-limb operands and form the five-limb
sum using carry instructions. Subtract the five-limb encoding of p. Select the
original sum if subtraction borrows, or the reduced sum otherwise. Store the
four result limbs. The sum is below 2p, so one subtraction suffices. Modular
subtraction uses a four-limb subtraction, then adds p if there is a borrow.
The correction can use a borrow mask and four AND instructions. Each kernel
needs fewer than 64 arithmetic, selection, and operand-transfer instructions.
The A=256 allowance also covers fixed constant loads, saved registers,
scratch transfers, pointer setup, entry, and return. The kernels have fixed
size and use local memory, within one 64-bit address space.

Multiplication uses binary shift-and-add. Start an accumulator at zero and a
running value at the first operand. Scan the 256 bits of the second operand.
For each bit, conditionally add the running value to the accumulator, then
double the running value modulo p. Both values stay canonical. Charge two A
calls and 64 other instructions per bit, plus 256 instructions for setup and
output. Thus

```text
M = 256*(2*A+64)+256 = 147712.
```

The 64-instruction allowance covers bit selection from a four-limb operand,
tests, branches, loop counters, and argument pointers. This kernel uses only
64-bit arithmetic; it does not use an uncharged large multiplication.

Inversion uses binary exponentiation to p-2. Fermat's identity applies because
p is the protected prime and the input is nonzero. At most 512 M calls and
512 blocks of 64 other instructions, plus 256 setup instructions, suffice:

```text
I = 512*(M+64)+256 = 75661568.
```

Point addition first tests the infinity, equal-coordinate, and opposite-point
cases. Return the other operand for an infinity input. Return infinity for
opposite points. For distinct x coordinates, use slope (v-y)/(u-x). For
doubling with nonzero y, use slope 3*x^2/(2*y). A zero y in doubling returns
infinity. After computing the slope t, set x'=t^2-x-u and
y'=t*(x-x')-y, all modulo p. These complete branches use at most one nonzero
inverse, 16 multiplications, and 32 modular additions or subtractions. Reserve
4096 further instructions for coordinate tests, record copies, selection,
and control:

```text
G = I+16*M+32*A+4096 = 78037248.
```

This overcounts the ordinary affine formulas. It also covers all exceptional
branches. Each scalar multiplication uses binary double-and-add and at most
512 G calls. There is no batch-inversion assumption and no incomplete-formula
exception left untreated.

## Storage, routing, sorting, and search

A baby record has twelve 64-bit words (96 bytes): one infinity tag, four x
limbs, four y limbs, two index limbs, and one padding word. The key has nine
words. The algorithm has two arrays of m records for merge sort.

The table is too large for a single 64-bit address space. To avoid a free
large-address RAM assumption, use a binary tree of storage nodes. Its leaves
hold the 2m records of the two buffers. There are fewer than 4m nodes. Each
node has its own local 64-bit address space, two child communication links,
and one parent link. A node's code, routing data, and fixed message scratch
use at most 64 KiB. Route by a 128-bit record number: one buffer bit followed
by the 127-bit record index. Keep the byte or word offset within a record
local to its leaf. Store routing addresses in three-word values; each routing
decision examines the next record-number bit. The actual depth is at most
128; the following bound allows 192 links each way.

Use a serial protocol with one outstanding request. A request follows the
address bits to its leaf. A read returns a twelve-word record. A write sends
a twelve-word record and returns an acknowledgement. At each link, copy the
bounded packet, update the routing bit position, and select the next link.
The packet needs at most sixteen words, including its address and control.
Charge at most 256 scalar instructions per link in each direction, including
all sender and receiver loads, stores, control, and coordination. For example,
four word transfers per packet word use 64 instructions; address tests,
counter changes, and handshake state fit in the remaining 192. No polling
loop runs while a request is absent. This is instruction-equivalent work;
network latency and physical distance do not reduce the charged work.

A complete read or write therefore costs at most

```text
F = 131072 = 2^17 instructions,
```

which covers 2*192*256=98304 link instructions and 32768 endpoint
instructions. This is a conservative serial transfer protocol, not a free
lookup oracle. All routing processors' work is included in F. No parallel
speedup is credited.

Build the tree by recursive initialization. A parent copies the fixed router
program to each child and initializes its state. The parent then delegates
initialization of the child's subtree. Thus each code word crosses each edge
at most once during construction; code is not sent from the root separately
to every leaf. Reserve 2^20*m instructions for all node creation,
code transfers, zero initialization, links, array initialization, and local
state. This is more than 2^18 instructions per node. A node has at most
8192 local 64-bit words, so initialization and program copies fit this bound
even with 16 instructions per word and a further equal control allowance.
The reserve is charged in full before the target search.

Use the following bounds for table operations:

```text
B = 2*F+256 = 262400               # compare two full keys
D = 2*F+256 = 262400               # read and move a full record
V = (L+1)*(B+256)+D+512 = 33882880 # one binary search and returned record
```

B permits two full record reads, even when a search key is already local.
Nine word comparisons, in order, plus three-word address and counter work
fit the extra 256 instructions. D permits one full read and one full write,
including local record construction. V permits L+1 key comparisons, another
256 instructions per comparison for index and interval updates, and a full
record move. This also bounds an unsuccessful binary search.

Bottom-up merge sort takes exactly L passes because m=2^L. Each pass has at
most m key comparisons and at most 2m record moves. Charge another 128m
instructions for loop and address work. Fixed pass setup gets 4096L:

```text
S = m*L*(B+2*D+128)+4096*L.
```

The separate B and D charges can read the same record more than once. This
is intentional overcounting. No expected constant-time hash lookup or free
word-width extension is used.

## Total work and memory

The baby loop takes at most m*(G+D) instructions. This includes its final
addition that obtains mP. The giant loop takes at most h*(G+V), including
an unused final giant update. Charge 2^40 fixed setup instructions and a
separate 2^40 final validation, verification, and output instructions.
These reserves include canonical decoding, curve and subgroup checks,
constant setup, D=-mP, scalar construction, and serialization. Two full
256-bit scalar multiplications need at most 1024G instructions, which is
less than 2^40; the remaining reserve is ample for these fixed-size tasks.
If a task is also charged in a loop, this only increases the bound.

| Component | Scalar instruction equivalents |
| --- | ---: |
| Tree and memory initialization | 178405961588244985132285746181186892047843328 |
| Baby table generation | 13321994775258162759421546327022817487297183744 |
| Merge sort | 17012528546828668548667429485422489351716597760 |
| Giant search | 30467556849547518968755367500570599593899732992 |
| Fixed setup | 1099511627776 |
| Final checks and output | 1099511627776 |

The exact total is

```text
T = 60980486133222595261976629059197095523984613376
T < 1000 * 2^146.
```

Preprocessing includes initialization, the complete baby table, its sort,
and the fixed setup allowance. It has no amortization credit:

```text
T_pre = 30512929283675076293221261558626494830573252608
T_pre < 1000 * 2^145.
```

The inequalities are integer comparisons. They do not depend on rounded
logarithms. Dividing the scalar counts by the protected C=1000 proves the
submitted total and preprocessing bounds. Point arithmetic is charged once
through G; no extra reference-operation charge is added to it.

For memory, charge 64 KiB for each of 4m nodes, both 96-byte record arrays,
and a separate 1 MiB for root code, input/output, field scratch, scalar
scratch, and control state:

```text
peak_bytes <= 4*m*65536 + 2*m*96 + 1048576
           = 44634157504285656375563920507610172761309184
           < 2^146.
```

This includes code copies, table and routing metadata, advice, and all
processor state. The leaf allowance and record arrays can overlap in an
implementation; charging both is safe. No randomness is retained. There is
no nonuniform advice, so advice_bytes=0 and the required
log2(max(1,advice_bytes)) is 0.

## Evidence, limits, and provenance

The construction and cost bounds above are analytic. No full-size run, wall
time measurement, reviewer acceptance, or practical key recovery is claimed.
The finite instruction and memory bounds describe an enormous classical
computation. The target's field, subgroup order, generator, input encoding,
and success test remain unchanged. A conservative bound does not show an
improvement over the best known attacks.

BSGS is a standard method; see Algorithm 3.56 in the primary textbook source,
[Handbook of Applied Cryptography, Chapter 3](https://cacr.uwaterloo.ca/hac/about/chap3.pdf).
The fixed curve parameters are specified in Section 2.4.1 of
[Standards for Efficient Cryptography 2, version 2.0](https://www.secg.org/sec2-v2.pdf).
The algorithm, success proof, and complete claimed cost are given here so
neither reference is needed as an uncharged algorithmic oracle.

The repository's original generic baseline was read. Its BSGS design and
the structure of its conservative arithmetic discussion informed this
accounting control. Its stated word width and final normalization are not
used: this argument uses four 64-bit field limbs, explicit large-storage
transfers, and the current protected C=1000. Existing submission files are
unchanged. There are no heuristic assumptions to extrapolate from a toy
curve or an ideal random walk.
