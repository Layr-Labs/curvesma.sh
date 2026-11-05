# Complete 64-bit cold-work and storage certificate

This prices proof.md with m=h=2^127 under classical-group-work-v2. All
arithmetic is charged through scalar AArch64 instruction equivalents,
H=0. The protected C=20000 is a normalization constant, not a candidate
operation price. No 256-bit or arbitrary-width arithmetic is unit cost.

## 1. Base, extension-field and point operations

A base-field element uses six 64-bit limbs; a seventh handles carry.
Canonical modular addition/subtraction costs at most A=2048: six passes
over seven limbs at sixteen instructions (672), with 1376 additional
instructions for copies, calls, fixed addresses and correction. A canonical
sum needs one subtraction of p; a difference needs at most one addition.

For a product u*v, use a result initially zero and running value u. Scan
all 384 padded bits of v, conditionally add the running value, then double
it modulo p. At most two A calls and 64 bit/control instructions per bit,
plus 2048 setup instructions, give

    M=384*(2*2048+64)+2048=1599488 <2^21.

Binary powering by p-2 for a nonzero value takes at most 768 products,
384*64 control instructions and 2048 setup, so I<2^31. Fermat inversion
is valid because nonzero multiplication permutes Fp*, yielding z^(p-1)=1.
Leading zero bits are retained in the finite bounds. There is no division
or inversion oracle, and zero denominators are handled before inversion.

Multiply Fp12 polynomials by 144 base-field products and 144 coefficient
additions; then eliminate degrees 22 down to 12 using
w^12=2w^6-2. Each of eleven eliminations uses at most four modular
additions/subtractions, for 188 total. This includes possible later
high-degree terms: descending elimination consumes each coefficient once.
Allow another 65536 instructions for coefficient addressing, copies,
zeroing, loops and transfer of at most 23 six-limb coefficient cells.
Thus 144M+188A+65536=230776832<2^28. We use the looser envelope

    E=2^30

for every Fp12 product or sum/subtraction, including data movement and
control. Fp2 multiplication (four M and two A) and addition (two A),
with 4096 extra instructions, also fit E. Fp2 inverse uses
(a+bu)^(-1)=(a-bu)/(a^2+b^2), charging one base inverse plus products,
copies and control; its cost is below 2^33. Canonical comparisons inspect
all 72 Fp12 limbs or 12 Fp2 limbs, never a 256-bit unit word.

Complete projective add/double formulas, branch cases, twist, rational
line evaluation and the exact Miller schedule are in pairing-kernels.md.
Each Fp/Fp2 projective point operation fits 64E. Scalar multiplication of
a 256-bit nonnegative integer takes at most 512 point operations, below
2^45 instructions. Normalizing a projective point uses one Fp/Fp2 inverse
and a fixed number of products; infinity is tested first. It fits 2^35.
All subgroup tests and returned-point formation therefore have finite
priced implementations for both G1 and G2.

The appendix proves an exact pairing evaluation costs

    J=(128*512+48*384)*E+2^28 <2^47.

No precomputed free pairing or field table is used. The final Fp12
quotient and exponent, decoding and exponent construction are included.
Static formulas are implemented directly; Python module import-time
tests, precomputed exponent tables or a Python interpreter are not part
of the defined algorithm. Their omitted work supplies no premise.

## 2. Concrete storage and communication

Each record is 1024 bytes: 72 coordinate limbs (576 bytes), two index
limbs (16 bytes), and padding. There are two m-record arrays, one for
merge scratch. Their payload is 2*m*1024=2^138 bytes. A scalar index,
one-past endpoint and recovered product i*m+j fit four 64-bit limbs.

As a single architectural address cannot name the arrays, use a declared
bounded-degree binary storage tree with 2m=2^128 leaves and fewer than
2^129 nodes. Addresses consist of a buffer bit and 127-bit record index;
route one bit per local edge for at most 128 hops each way. A controller
serializes logical transactions. All processors' work is summed, not
parallel elapsed time. Each node uses a small local address space.

An edge exposes two directional bounded memory-mapped mailboxes with
valid flags and blocking receive events. There are no spurious wakeups,
errors, retries or spinning in this explicit hypothetical interface.
Each elementary packet has at most 24 scalar words: sixteen payload words
and at most eight header/address/control words. Charge 16 instructions
per word across sender/receiver for loads, stores and transfer (384),
and 640 for flags, barriers, waits, acknowledgment, address-bit selection
and fixed branches, totaling 1024 per elementary hop. This is the stated
local-event machine assumption, not proof of a real cluster's timing.

A full record is eight 16-word payload packets. Request, fragmentation,
sequence checking and acknowledgment fit sixteen elementary packet-hop
allowances per logical edge, at most 2^14 instructions. A path and
return of at most 256 edges plus endpoint assembly and controller work
therefore fits F=2^24 instructions for a full read or write. There is
one in-flight transaction and no unbounded queue or free scheduler.
All address manipulation is multiword arithmetic. Bound a comparison of
two remote keys and a full record move by B=D=2^26 each: two F transfers
plus local 72-limb comparisons, endpoint updates and control fit this.

Reserve 64 KiB per node, including fixed reset loader, routing code,
labels, links, all mailbox endpoints, a record and scratch/processor
state. A routing/fragmentation program of at most 2048 fixed instructions
occupies 8 KiB; the remaining 56 KiB is ample for the bounded packets,
record and local state. The controller gets a separate 1 MiB for code,
field/point scratch, parsing and exponent buffers. Thus total memory is
below 2^129*2^16+2^20=2^145+2^20<2^146 bytes, including code replicas.

The initial machine has a uniform fixed receive-loader, at most 64
instructions already within that code allowance. Charge initial placement
as 1024 instructions per node rather than free deployment. A parent
installs up to 1024 code words in 64 sequential 16-word packets with
chunk numbers and acknowledgment. Allow four directed hops per packet
at 1024 instructions (262144), 131072 for zeroing at most 8192 local
words, and 262144 for labels, control, buffers and completion. With
loader placement the total is 656384, below the declared 2^23 per-node
bootstrap reserve. Count each parent-child installation against the child,
including both endpoints' work; root work is in fixed controller setup.
Total node/table initialization is below 2^129*2^23=2^152 instructions.
This explicit reset/stream/event interface is a material assumption.

## 3. Sorting, lookup, failure and total work

Bottom-up merge sort uses 127 passes, each with at most m comparisons,
2m record moves and 4096m remaining multiword loop/boundary instructions:

    S<=127*m*(B+2D+4096)+4096*127 <2^163.

Binary search takes at most 128 iterations, each with a routed read,
full-key comparison and endpoint updates inside B+4096. Allow one D
to recover the full returned record and 4096 result-control instructions:

    V<=128*(B+4096)+D+4096 <2^34.

The baby phase costs at most m*(E+D)<2^158. The giant phase costs at
most h*(E+V)<2^162, including a conservatively charged last multiply.
All three search terms are individually below 2^163. Recurring sample/
search loop overhead costs 4096*(m+h)=2^140, including failed iterations.
Every lookup compares all coefficients; no collision oracle is assumed.

Reserve 2^70 instructions for all remaining fixed work: parsing compact
input and installing controller constants/code; field/curve/prime-subgroup
validation; computing g; checking its order; computing D by binary powering;
forming the recovered point; final subgroup and pairing verification;
normalizing and serializing output; and scalar-index setup. At most two
pairings (each J<2^47), fewer than eight 256-bit scalar multiplications
(each <2^45), fewer than eight 256-bit Fp12 powers (each <2^40), and
fixed decoding/control well below 2^40 fit this reserve. Pairing exponent
construction and both field inversions are already charged inside J.
There is at most one candidate verification, justified by correctness.

Consequently the full worst-case instruction count, success or failure, is

    R<3*2^163+2^152+2^140+2^70 <2^165.
    W=R/20000 <2^151, since 20000>2^14.

The declared rational score is 151/1, a conservative log2 upper bound.
Preprocessing is a subset of total W and shares this bound, without
amortization. The two tables are regenerated for the supplied input g.
No randomness, amplification, external searched advice or quantum step
is used. The compact representation has <4096 input bytes and constant
output size; its two-pass parsing at 64 instructions/byte plus 2^16
fixed work fits setup. Expanded spellings have the separately declared
128L+2^16 instruction charge and do not inherit the fixed parsing bound.
