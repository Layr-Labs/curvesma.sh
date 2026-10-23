# A deterministic accounting control for unrestricted secp256k1

## Claim and scope

The protected profile, rules, price C=1000, and classical-group-work-v2 policy
apply unchanged. Let p and r be their field modulus and prime subgroup order.
The attack receives only Q=[x]P, with x uniform in 0..r-1. Set

    m = 2^127, h = 2^128, B = 4096, K = m*h = 2^255.

This is an analytic finite algorithm; no production execution is claimed. It
has no heuristic assumption, precomputed answer, free table, random oracle,
quantum component or amortization across targets. The declared bound is
2^143 total target-group-operation equivalents. All computation below is
charged as scalar 64-bit AArch64 instruction equivalents and divided by the
fixed C, so H=0 and W=R/1000. No reference addition is counted a second time.

## Algorithm and exact success

Represent a point table record by a 64-bit infinity flag, four 64-bit limbs
of canonical x, four of canonical y, two of the 127-bit index j, and five
padding words: 16 words = 128 bytes. For infinity, both coordinates are zero.
Order keys lexicographically by (flag,x,y), comparing coordinate limbs in
big-endian significance. Record indices and counters use two or four 64-bit
limbs as needed; there is no arbitrary-width unit-cost instruction.

1. Canonically decode Q. Return failure on malformed input. For valid input,
   check the curve equation and [r]Q=O. Infinity is a valid input.
2. Starting at O, generate jP for j=0..m-1 by Jacobian mixed additions with
   the fixed affine P. Retain the endpoint mP. Normalize points in blocks of
   exactly B, writing records to banked array A. This makes m additions,
   conservatively including the endpoint addition. Normalize mP separately.
3. Sort A by the complete point key with bottom-up stable merge sort, using
   a second m-record array. There are exactly 127 merge passes. A run cursor,
   its two input positions and output position have bounded multiword size.
   No hashing, fingerprint, collision assumption or free sorting is used.
4. Set affine D=-(mP). Generate gamma_i=Q+[i]D for i=0..h-1 in Jacobian
   blocks, normalize each block and binary-search every resulting key in
   the sorted baby table. Retain the endpoint of a block for the next block.
   Charge h additions, including the unneeded final update. At a match j,
   form a=i*m+j (a 255-bit integer), independently check [a]P=Q, encode it
   as the required 32-byte scalar and return it. Exhaustion returns failure.

For 0<=x<K, Euclidean division gives unique 0<=i<h and 0<=j<m with
x=im+j, and gamma_i=[j]P. The complete table therefore finds x. Conversely,
a matching pair produces a in [0,K), with [a]P=Q. Since K<r, uniqueness
modulo r implies a=x. Thus valid inputs with x>=K have no match; returning
failure there is deliberate and fully charged. In particular x=0 and x=K-1
are covered, and x=K is not. Success is exactly K/r>1/2, since the protected
order satisfies 2^255<r<2^256. This is success on the original uniform full
range, not a promise that the secret lies in a shortened interval. There are
no algorithm coins, restarts or amplification and no probabilistic stopping.

## Elementary arithmetic on actual 64-bit words

The following are deliberately loose implementation ceilings, not measured
compiler output. Field values use four 64-bit limbs; a carry uses a fifth.

* A modular add or subtract of canonical values costs A0=256 instructions.
  Load operands into a five-limb zero-extended temporary, add/subtract with
  carry/borrow, perform the one needed conditional correction by p, and
  store four canonical limbs. At most 20 limb loads/stores for the operands,
  temporary and result, 10 carry/borrow arithmetic instructions, 20 for
  correction/comparison/selection, and 128 for addressing, branches and
  fixed loop control fit below 256. Constants may be loaded afresh. An
  implementation may unroll these fixed five-limb passes; no 256-bit word
  is assumed. The possible 257th bit is retained until correction.
* Multiply modulo p by the elementary bit-serial method: initialize z=0,
  u=a; for each of the 256 bits of b, conditionally set z=z+u mod p and
  set u=2u mod p. Two modular additions per bit, plus 64 instructions per
  bit for bit access, limb/index manipulation, moves and control, and 512
  for initialization/finalization, give

      256*(2*256+64)+512 = 147968 < M = 2^18.

  Squaring uses this same bound. Division is never a primitive.
* Invert a nonzero field value by binary exponentiation to p-2. At most
  512 multiplications and 256*64+512 other instructions cost less than
  I=2^28. This includes the multiplication calls. No zero is inverted.

Here is a complete mixed-Jacobian addition for the short Weierstrass a=0
curve. A Jacobian (X,Y,Z) denotes (X/Z^2,Y/Z^3), with Z=0 for infinity.
For a nonidentity affine addend (u,v), compute

    z2=Z^2; U=u*z2; S=v*Z*z2; H=U-X; t=S-Y.
    If H=0: double if t=0, otherwise return infinity.
    HH=H^2; J=4*H*HH; V=4*X*HH; q=2*t;
    X'=q^2-J-2*V;
    Y'=q*(V-X')-2*Y*J;
    Z'=(Z+H)^2-z2-HH.

Return the affine addend (u,v,1) when the accumulator is infinity, and
return the accumulator when the addend is infinity. Doubling uses

    A=X^2; C1=Y^2; C2=C1^2;
    D=2*((X+C1)^2-A-C2); E=3*A; F=E^2;
    X'=F-2*D; Y'=E*(D-X')-8*C2; Z'=2*Y*Z.

Y=0 doubles to infinity. These formulas, including every exceptional case,
take no more than 16 field multiplications/squarings, 64 modular additions,
and 4096 scalar instructions for tests, copying, counters, local memory and
control. The resulting ceiling is below 4,215,000 instructions per point.
Mixed addition may compute some unused intermediates before a branch; this
still fits the allowance. Since p>3, the formulas apply to this profile.

Normalize a B-point block by multiplying the nonzero Z's into prefix
products, inverting the final product once, and traversing backwards to
recover every Z inverse. Replace each zero Z by 1 in the product and retain
its infinity flag. At most 3B multiplications suffice for the prefix and
backward pass. Compute zinv^2, X*zinv^2, zinv^3 and Y*zinv^3 using four more
multiplications per point. Empty nonzero sets need no inversion. Thus the
normalization charge is at most 7M+I/B per point. The endpoints of the baby
and giant loops are powers of two and B divides both exactly. The one
separate endpoint normalization is charged in fixed work.

Including mixed addition, block scratch transfers, counters, record-key
construction, flag handling and an additional 8192 local instructions per
point, use the convenient combined ceiling

    G = 2^23 instructions per generated, normalized point.

Indeed 16M+64*A0+4096+7M+I/B+8192 < G. The state needed to generate a block
fits in 2^24 bytes, including prefix products and output points. Record
access to the huge table is additional and charged separately below.
Verification can use ordinary double-and-add with these complete cases;
its at most a few thousand point operations and standalone inversions fit
far below the fixed allowance below.

## Addressing and data transfers beyond a 64-bit address space

Two arrays occupy 2*m*128=2^135 bytes. They are not placed in one processor's
64-bit virtual address space. Use 2^75 storage banks, each containing at
most 2^60 bytes of record slots. Two sets of 2^74 banks hold the two arrays.
The high 74 bits of an array index select a bank; its low 53 bits give a
record offset, multiplied by 128. A separate array-selector bit makes the
75-bit bank identifier. Local offsets and all local addresses fit 64 bits.

Connect these banks as leaves of a complete binary routing tree of depth
75. There are fewer than 2^76 total leaf/internal nodes. A controller at the
root issues one table access at a time. This deliberately serial schedule
needs no contention or bandwidth heuristic. Nodes have fixed local buffers,
two child links and one parent link; one request and one reply pass through
at most 150 edges in total. Use a fixed 20-word packet: the 16-word record
plus four words for command, bank identifier and offset/return information.
An access returns a record or an acknowledgement. A write sends the record
in the request. A read returns it in the reply. Even a small request is
charged at the full packet size.

For each traversed edge, charge 128 scalar instructions: at most 80 loads
and stores for copying the 20 words through the sender/receiver link buffers
(count both endpoints), and 48 for packet header handling, routing-bit
selection, addresses, ready/acknowledgement bookkeeping, branches and
control. Link storage is point-to-point shared message storage; this
charges the memory operations on both sides. There is no DMA or uncharged
copy engine. Latency and wall time are not the metric, and idle processors
need execute no polling loop: the controller schedules the finite forward
and reverse route explicitly, with the per-edge bookkeeping charged above.
Node selection is one bit of the two-limb bank ID at each tree level.

At most 8192 extra instructions per access cover the root's multiword
index calculation, node-local record read/write, protocol setup/teardown,
and local address arithmetic. Consequently the uniform per-access ceiling

    T = 2^15 > 150*128 + 8192 = 27392

includes all processors' work and transfers. This is an analytical storage
network in the fixed scalar-instruction model, not a claim that such a
physical computer exists. There is no unit-time large-address dereference.
Fixed program and topology construction are also charged below. The message
storage model performs explicitly bounded copies and does not assume a free
remote lookup or conceal the routing costs in elapsed parallel time.

## Sorting, lookup and complete total

Each merge pass has at most m key comparisons and 2m record moves (a safe
overcount of a standard out-of-place merge). A comparison reads two records,
a move reads and writes one, so at most 6m accesses per pass. Allow another
4096m instructions per pass for local key comparisons, two-limb run/index
arithmetic, moves, loop control and run-boundary handling. There are 127
passes. Initial table construction writes each record once. These are the
only large arrays; generation/normalization uses local block scratch.

Binary search among m=2^127 keys uses at most 128 record probes. Each probe
costs T+512, where 512 additional instructions cover the nine-word key
comparison, two-limb midpoint and interval update, branches and control.
Another 4096 per query covers a matched index, scalar reconstruction and
surrounding control. The matched record is already in the last probe.

Charge F0=2^104 scalar instructions once for ALL fixed startup work:
program/topology construction and distribution, bank initialization, root
scratch and local buffers, public constants, parsing, membership check,
endpoint normalization, step negation, final verification and output, even
on failed inputs. A concrete allowance of less than 2^28 instructions per
node for initialization, constant/program copying and a root phase of less
than 2^60 instructions fits F0 because the number of nodes is less than
2^76; more specifically a full depth-75 tree has 2^76-1 nodes, and reserving
2^27 instructions per node leaves over 2^103 instructions for root work.
Unwritten record slots are not read and need no eager zeroing. The record
slots actually touched are charged in the table access/merge costs.

The complete hard instruction bound is therefore

    R = F0 + (m+h)*G + m*T
        + 127*m*(6*T+4096)
        + h*(128*(T+512)+4096).

With h=2m this simplifies exactly to

    R = 2^104 + 59215872 * 2^127.

The test for the DECLARED bound uses integers, without floating-point
rounding: R < 1000 * 2^143. All work on failures is bounded by these same
loops. Thus total classical computation W=R/1000 is less than 2^143 target
units. This is a deliberately loose submitted upper bound, not an automatic
tightening request. Preprocessing is fully included and bounded by

    R_pre = F0 + m*G + m*T + 127*m*(6*T+4096)
          = 2^104 + 33910784 * 2^127 < R,

so its reported bound of 2^143 units is also valid. No table is reused for
another target at discounted cost. The deterministic algorithm uses no
randomness, so its randomness cost is zero. Trials, failures, verification,
output, transfers and the sum of processor work are all accounted for.

## Peak memory and advice

The two record arrays use exactly 2^135 bytes. Reserve 2^32 additional bytes
per leaf/internal node, plus 2^32 bytes at the root for code, program copies,
public constants, packet buffers, processor registers/stacks, local block
scratch, endpoints, inputs/outputs, and every cursor. Fewer than 2^76 nodes
give total peak below 2^135+2^108+2^32 < 2^136 bytes. This includes retained
programs and all routing storage. The code and the described finite kernels
are small compared with the reserve; unused reserve is not eagerly written.
There is no stored advice. Advice bytes are exactly zero, and the reporting
convention log2(max(1,advice_bytes)) gives 0. Tables are constructed during
charged preprocessing and included in memory, not treated as free advice.

## Evidence and research context

The included audit checks exact integer inequalities and identities only.
Exhaustive small additive-group tests check coverage and intentional failure
outside the searched interval. Small samples of production-curve arithmetic
check mixed addition/doubling and batch normalization against an independent
affine formula, including infinities and inverse pairs. They do not execute
the large table, measure AArch64 instruction counts or prove the production
bound. The bound is the analytic argument above, open to independent review.

Heuristic premises: none. Engineering limitations: astronomical processor
count, storage and computation; analytical, unmeasured kernel/routing bounds;
no production-size execution; no claim to beat the published rho estimates.
The explicitly charged network is part of the model implementation, not an
assumption that an actual AArch64 installation has this address space.

The standard BSGS construction is described in Handbook of Applied
Cryptography, Algorithm 3.56, section 3.6.2:
https://cacr.uwaterloo.ca/hac/about/chap3.pdf . SEC 2 section 2.4.1 gives
the public curve parameters: https://www.secg.org/sec2-v2.pdf . The necessary
construction and formulas are reproduced here so review needs no web access.
The repository's organizer BSGS baseline motivated this accounting-control
exercise, but its old price and word-width ledger are not premises here.
Public six-class-rho and plain-rho notes were read for context; no walk,
meter, constants, empirical success data or code from them is used. No
previous acceptance or reported score is invoked as mathematical evidence.
