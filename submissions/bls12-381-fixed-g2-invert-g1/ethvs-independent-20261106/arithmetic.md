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

For pairing tracks use exactly the protected reduced py_ecc pairing at
commit 04151f01f59f902ab932a51e0ca0ebce3883fc51. This is a translation of
the fixed algebraic expressions into the charged limb/field operations,
not a claim that importing Python has the same instruction count. We do
not execute module-import fixtures or compute the unused exp_by_p table.
The caller uses pairing(Q,P,final_exponentiate=True), which selects the
direct exponentiation path in miller_loop, not the separate optimized
final_exponentiate function. Identity inputs produce the 12-coefficient
GT identity after the charged curve/subgroup check.

The source's homogeneous coordinates mean (x/z,y/z). Set cast_P to its
three coefficients embedded at degree zero. For Q in Fp2, the exact
twist maps x=(x0+x1*u), y and z as follows, with u=w^6-1:

  nx has coefficients (x0-x1) at degree 1 and x1 at degree 7;
  ny has (y0-y1) at degree 0 and y1 at degree 6;
  nz has (z0-z1) at degree 3 and z1 at degree 9;

all other coefficients are zero. Retain this exact scaling for every
Miller point, rather than substituting an arbitrary twist. Start R=Q,
twist_R=twist_Q=twist(Q), f_num=f_den=1. Use the fixed 63 entries
pseudo_binary_encoding[62::-1] from the source whose bit expansion is
the protected ate_loop_count. At each entry call the exact line
formula below on twist_R,twist_R,cast_P, multiply f_num by itself
and its line numerator, and similarly for f_den/denominator. Double R
with the source homogeneous double formula and retwist. If the bit is
1, evaluate the same line formula on twist_R,twist_Q,cast_P, multiply
both fractions, add Q to R with the source homogeneous addition and
retwist. There are <=126 line calls, <=63 doubles and <=63 additions.

For line inputs (x1,y1,z1),(x2,y2,z2),(xt,yt,zt), set
mn=y2*z1-y1*z2 and md=x2*z1-x1*z2. If md!=0 return
(mn*(xt*z1-x1*zt)-md*(yt*z1-y1*zt), md*zt*z1).
If md=0 and mn=0, replace mn=3*x1*x1, md=2*y1*z1 and
use that same return expression. Otherwise return
(xt*z1-x1*zt,z1*zt). These are all source branches; no omitted
zero-denominator shortcut, raw Miller equality, Frobenius correction
or conjugation is silently introduced.

The homogeneous double sets W=3*x^2, S=y*z, B=x*y*S,
H=W^2-8*B; return (2*H*S, W*(4*B-H)-8*y^2*S^2,8*S^3).
For add, identity z=0 returns the other input; define
U1=y2*z1,U2=y1*z2,V1=x2*z1,V2=x1*z2. If V1=V2 and U1=U2,
double; if V1=V2 only, return (1,1,0). Otherwise U=U1-U2,
V=V1-V2, W=z1*z2, A=U^2*W-V^3-2*V^2*V2; return
(V*A, U*(V^2*V2-A)-V^3*U2, V^3*W). These formulas and
line outputs are reproduced from the exact source, so the operation
count is bounded directly rather than treating an external call as free.

After the 63 steps, compute f=f_num*f_den^(p^12-2), then
f^((p^12-1)/r) with the bounded binary exponentiation above.
Each line uses <32 dense Fp12 multiplications and at most 256
coefficient operations/comparisons, including tangent scalar factors,
zero comparisons and complete copies/control. Each homogeneous point step uses <64 Fp2
multiplications and <256 other bounded operations; each twist uses
<256 coefficient copies/subtractions. These ceilings include small
scalar factors, exceptional branches and full copied coordinates. Even pricing ALL <2^17 Miller,
point,twist,comparison operations as a dense multiply bound 2^33
leaves <2^50 work. The two bounded Fp12 exponentiations cost <2^48
each; loop/state, input and final copies fit <2^50, so the whole
nonidentity evaluator is <2^52 instructions. Conservatively reserve
2^70 for each complete pairing including subgroup/encoding validation;
only a fixed number of them is used, and each is included in Section D.
For nonidentity inputs the denominator is nonzero, with the following
explicit lemma. Each loop prefix k and doubled prefix lies between 1 and
the positive ate_loop_count, which is <r-1; r is an odd prime. Thus [k]Q
is neither infinity nor 2-torsion, so its homogeneous z and its tangent
y factors are nonzero. In an addition step the doubled prefix 2k is
between 2 and ate_loop_count<r-1; hence [2k]Q is neither Q nor -Q, so
the homogeneous x-difference is nonzero. Twist multiplies these factors
only by nonzero powers of w in the field, and the finite G1 evaluation
point has nonzero z. Every line denominator and their product is
therefore nonzero. This proves exponentiation gives the exact source
division for every valid nonidentity pair. The identity case is
dispatched first exactly as the source.

Reference file hashes are included in the protected profile reproduced
in proof.md. The four downloaded source files were checked against those
hashes; a source match does not by itself prove mathematical acceptance.

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
Connect the two roots via their previously unused parent ports: one
root is the controller, the other root relays requests to the second
buffer tree. A one-bit buffer selector plus the four-limb leaf index
chooses the tree. The controller handles tree0 locally or sends tree1
requests across the root bridge; the other root replies over that same
bridge. Thus one traversal has at most 129 edges, and the later allowance
of 256 request and 256 reply/control crossings covers both trees. The
same bridge carries the child-loader code/image at cold start: controller
boots its own children and the second root, which forwards to its own
children. Degree remains at most3. Only one root transaction is active,
so no queues grow with the table.
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
so the whole image broadcast costs <2^151 instructions. Each node then clears only its designated data/scratch/record regions,
excluding the read-only child-loader image, received controller code,
current code pointers, and live framing/ack state. The code image is
installed in a fixed code region; the loader passes its complete-state
flag and entry address before jumping, so no executing bytes or live
handoff metadata are erased. At most <=2^20 local bytes are cleared (<=2^17 words, <=16 instructions per
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
