# Finite arithmetic, memory and communication certificate

This is an analytic scalar-64-bit implementation certificate, not a run of
an attack, a measurement of a device, or an assertion that existing hardware
can hold the table. All steps below are charged, including cold initialization.

## A. Integers and prime-field operations

All scalar and base-field values are arrays of at most six 64-bit limbs
(the selected base prime is at most 384 bits); arithmetic is not unit-cost
on an unlimited integer. A 256-bit loop index uses four limbs. Compare,
copy, add, subtract, shift by one, carry/borrow propagation, zero testing,
load/store and loop control on at most six limbs use fewer than 2^10 scalar
instructions, including scratch writes. Modular addition/subtraction uses
one unreduced limb operation and a comparison plus at most two corrections,
with an extra carry limb. Charge 2^12 instructions to the complete operation.

Use 384-step binary double-and-add modular multiplication: start c=0,
scan the multiplier bits least-significant first, conditionally add the current multiplicand to c,
and double the multiplicand modulo the exact prime. Each step has at most
two modular additions plus bit selection and fixed-length loop handling;
2^14 instructions per step is a conservative bound. Its work is <2^23
instructions, including input/output copies. For a nonzero denominator
use binary square-and-multiply to exponent p-2; at most 768 modular
multiplications and their control give <2^34 instructions. The field
identity a^(p-2)=a^-1 follows because p is the exact prime in the protected
profile. Zero denominators are tested before inversion.

For G1 use complete affine branching. O+A=A; A+O=A. For finite A=(x,y),
B=(u,v), if x=u and y=-v, return O. This includes doubling with y=0.
Otherwise, for x!=u set lambda=(v-y)/(u-x); for A=B with y!=0 set
lambda=(3*x^2+a)/(2*y), where a is the selected profile's actual curve
coefficient (a=-3 for P-256, a=0 for secp256k1 and BLS G1).
Return (lambda^2-x-u, lambda*(x-x_new)-y), reduced canonically.
The prime is odd; x=u on a nonsingular curve forces v=y or v=-y,
so these branches cover every valid pair. Fewer than 32 field
multiplications, one inversion, and 128 other limb operations suffice.
Charge <2^36 instructions per complete G1 addition, including dispatch.
This deliberately does not borrow an incomplete 6M reference price.

Fp2 is c0+c1*u with u^2=-1. Multiplication is the expanded four base
products with signed combination; inverse of nonzero a+b*u is
(a-b*u)/(a^2+b^2). The norm is nonzero because the protected Fp2 is a
field. Complete G2 addition uses the same affine branches with these
Fp2 operations. Charge <2^40 instructions per G2 addition. Scalar
multiplication by any scalar below 2^256 uses at most 512 complete
additions plus bit handling; charge <2^50 including canonical output.
Input validation checks encoding, the actual curve equation and [r]A=O;
all such tests are charged, never mistaken for organizer work.

## B. GT arithmetic and the fixed reduced pairing

GT elements are the 12 canonical base-field coefficients in the protected
polynomial basis Fp[w]/(w^12-2*w^6+2), not a free Fp12 primitive.
A dense multiply forms 144 coefficient products and adds them into
23 coefficients. Reduce degrees d=22 down to 12 using
w^d=2*w^(d-6)-2*w^(d-12). This needs at most 144+44 coefficient
additions and small-constant operations in addition to the products.
Charge <2^33 instructions per dense multiply, including coefficient
copies. Fixed exponents have fewer than 4608 bits; binary exponentiation
uses fewer than 9216 dense multiplies and therefore <2^48 instructions.
The inverse of nonzero a in Fp12 is a^(p^12-2). All 12 coefficients
are compared; a zero element is not a GT input.

For pairing tracks use exactly the reduced pairing and argument order in
the protected profile, not raw Miller equality or a substituted pairing.
Its pinned py_ecc definition uses the given 64-bit ate loop, a rational
line-function Miller evaluation in Fp12, twist u=w^6-1 into that same
basis, and final exponent (p^12-1)/r. No result of the pairing is treated
as advice: computing the generator for the supplied fixed point and the
final returned-point verification are included. A once-per-attack budget
of 2^70 scalar instructions is reserved for each complete pairing and
its fixed-exponent arithmetic. The method requires the profile's exact
pinned pairing evaluator and does not assert that a different bilinear
pairing has identical coefficient outputs. For non-pairing tracks this
routine is not used. The arithmetic certificate's pairing coverage must
be independently checked against the protected reference before claiming
those two candidates ready.

## C. Explicit local-node packet protocol and cold bootstrap

The M=2^128 baby entries do not fit a 64-bit address space. Model a finite
binary tree of M leaves and height 128, using separate processors with
only 64-bit local addresses and at most 2^20 bytes of local memory. Two
such trees provide source and destination merge buffers: 4M-2 nodes
<2^130. Every node has three fixed ports (parent, left child, right child)
whose local 64-bit registers expose one incoming/outgoing 64-bit word,
valid/ready status, and an acknowledgment. This is the explicitly assumed
local-link machine interface; it is not a claim that such hardware exists,
that physical latency is small, or that global memory is unit-time. A
register load/store or a link-word load/store costs one scalar AArch64
instruction equivalent. If a deployment does not supply this interface,
this analytic machine assumption does not establish its physical cost.
Only one root transaction is active, so no queues grow with the table.
Blocking ready/ack instructions suspend the node; they are not unbounded
busy-wait loops. No background node spins. All active-processor work,
not blocked elapsed time, is summed. This stipulation is a machine
assumption, not an experimentally verified hardware capability.

A request is at most 256 64-bit words: four index limbs, an opcode,
flags, 128 record words, and bounded control fields. A full record is
1,024 bytes; every exact key fits. A routed transaction has at most
256 request crossings and 256 reply/control crossings. At each crossing
copy at most 256 words. For each word a fixed relay routine performs
at most 64 instructions: <=8 local load/stores, <=8 port reads/writes,
<=8 valid/ack/control operations, <=8 index/mask/length operations,
<=8 branches/counter operations, plus <=24 spare operations for fixed
loop entry/exit and temporary state. Port blocking charges the issue,
completion and acknowledgment, not an unexecuted waiting loop. At each
node an at most 2^12-instruction fixed setup selects the next port using
the appropriate index bit (a shift/mask on one of the four limbs), writes
its 256-word local request/reply buffer, and verifies the opcode.
A leaf copies <=128 words into/out of its fixed record array; comparison
of full 12-field or affine keys is <2^12 additional instructions.
The complete transaction is bounded by

    512*256*64 + 512*2^12 + 2^13 < 2^25 instructions.

Charge the larger 2^26 bound per routed access/comparison, including
both directions and complete record traffic. A binary-search comparison
may use two accesses; the later bound of 16 accesses per step covers it.
All indices, interval endpoints and counters use at most four 64-bit
limbs; arithmetic is bounded-limb, not unlimited-word unit cost.

Cold start assumes this fixed hardware topology with uninitialized table
data and one root controller. Every non-root processor starts at a fixed
read-only <=2^10-byte generic child-loader image, with a <=2^12-byte
loader buffer and reset program counter. This image is algorithm code,
not searched advice; it only receives length-delimited code/control
words from the parent, stores them, forwards the same words to child
ports using the relay above, acknowledges completion and jumps to the
received controller entry. Its instructions are charged exactly like
all other relay instructions. No routing index, target-dependent record
or initialized table is provided by that loader. Loader code/buffer are
counted in the node code/control allowance. If the machine requires
loading rather than ROM reset, explicitly write this image to each node:
<=2^7 words times <=64 charged store/transfer instructions and <=2^12
reset/control instructions per node, totaling <2^144 across <2^130 nodes;
this reserved finite installation work is included in cold start.
The presence of such a program-load/reset interface is a declared
fixed-machine assumption, not an uncharged ready-made search service.
The root broadcasts a <=2^18-byte fixed code/control image (<=2^15
words) down the tree using those child loaders. Each link forwards each
image word with the same <=64-instruction relay. At most 2^130 child links exist,
so the whole image broadcast costs <2^151 instructions. Each node
then clears <=2^20 local bytes (<=2^17 words, <=16 instructions per
word for store/address/loop control), costing <2^21 per node. Include
four-limb logical-prefix/depth setup, port configuration, end-of-image
acknowledgment and invocation: allow a further 2^22 per node. All-node
clearing/setup costs <2^153, including the root; charge the entire
broadcast, initialization and acknowledgment process <2^156. Baby
records, program copies and buffers are therefore computed/initialized
in this cold-work accounting, not available free from an earlier attack.
The fixed root controller and generic child loaders are algorithm-design
code, not searched advice. Their storage, active instructions and the
explicit installation allowance above are all included.
Physical manufacture/topology acquisition is not instruction work in the
protected model; its fixed presence is a declared modeling assumption.

Per-node live-state budget: code/control image <=2^18 bytes; request and
reply buffers <=2^12 each; leaf record <=2^10; up to 256 frames of 64
bytes <=2^14 (the proposed iterative protocol actually uses fewer);
register/port/scratch/control <=2^16 bytes; reserved bounded arithmetic
workspace <=2^17 bytes. The sum is less than 2^20 bytes, including flags,
acknowledgments, stack, counters and processor state. Root scratch is
bounded by the same allowance. Both trees have <2^130 nodes, hence all
replicated code, buffers, records and processors occupy <2^150 bytes;
the submitted 2^151-byte bound has a factor-two margin. No RNG history
is retained, because the Shanks algorithm uses no randomness.

Bottom-up merge sort has exactly 128 passes. A sequential root
controller routes exact records and maintains four-limb source/destination
indices; <=16M routed accesses per pass cover reading both heads,
comparison, copying one winner, final tails and pointer control. Thus
128*16M*2^26=2^165 bounds routed work. Charge all sort work <2^166.
Each of M giant lookups has <=129 binary-search iterations and <=16
routed accesses per iteration. 129*16*M*2^26<2^166 covers its routed
work; charge the whole lookup phase <2^167, including scalar endpoints.
Full digests/coordinates and infinity tags are compared, never a hash
fingerprint or sign-only quotient. All processors' active instructions
are included in these bounds. End-of-sort buffers are reused, not added
as another unbounded memory allocation.

## D. Complete cold-work inequality

At most M baby iterations and M giant iterations occur. An additive
iteration costs <2^36 instructions; GT multiplication costs <2^33.
Thus additive table/giant arithmetic totals <2^166. Initial scalar
multiplications, input parsing/subgroup checks, precomputation, final
verification and output fit below a further 2^166 bound. Section C
separately charges sorting <2^166, lookup <2^167 and cold start <2^156.
The sum of those bounds and arithmetic <2^166 is <2^169, leaving
margin below the declared 2^170 total. For an auxiliary-power input charge (n+1)*2^52 instructions
as well to stream, decode and validate every supplied point; the auxiliary
powers are not used to reduce search work. In every selected track use

    R(n) < 2^170 + (n+1)*2^52 scalar instruction equivalents;

omit the second term when there is no n. Each operation is charged in R
rather than also charging reference group operations. Divide the entire
R by the protected C: 1000 for P-256/secp256k1, 2000 for BLS G1,
20000 for pairing tracks. At n=2^20 the second term is <2^73.
Therefore R/C<2^162 for every selected track, and the declared time_log2
162 is a conservative upper bound. This integer declaration is not an
estimate of a fast implementation. Memory has no scalar tie-break in
these tracks, but its immense bound is explicitly reported. Probability
is 1 on every valid challenge, so no trials, random coins, failures or
amplification are omitted. Arithmetic invalid-input rejection is included
but the challenge distribution supplies valid inputs.
