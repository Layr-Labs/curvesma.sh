# Deterministic half-coverage BSGS: complete analytic argument

Status: candidate for hosted independent review; no qualification is asserted.
The first claim's exact submitted score is **146/1**. This is a deliberately
conservative bound, not a best-known attack estimate or practical curve break.

## 1. Protected instance, probability space and constants

Use exactly `tracks/secp256k1-ecdlp/target-profile.json`: the standard curve
Y^2=X^3+7 over its 256-bit prime p, standard generator P, prime order r,
and cofactor one. Only Q=[x]P is supplied. The unknown x is uniform on
0,...,r-1. No bits, interval promise, auxiliary powers, secrets, or sampling
coins are supplied. Canonical JSON, point and scalar rules are unchanged.

Set m=2^127, h=2^128 and K=m*h=2^255. These are algorithm parameters,
not restrictions on the challenge distribution. Direct comparison of the
protected r with powers of two gives K<r<2^256.

The submitted supplementary bounds are:

| Quantity | Declared logarithmic upper bound / value |
|---|---:|
| Total computation, target-group equivalents | 146 |
| Peak total bytes | 150 |
| Preprocessing equivalents, included in total | 145 |
| Success lower bound | 1/2 |
| Nonuniform advice bytes | 0 |
| log2(max(1, advice bytes)) | 0 |

## 2. Exact algorithm

Represent infinity separately from ordinary points. Use complete case handling
for affine addition as specified below; do not apply an incomplete formula to
exceptional inputs. All coordinates used as internal keys are canonical
Montgomery residues, an injective representation of field elements.

1. Parse the bounded-size input, check its canonical encoding, curve equation
   and subgroup membership. Infinity is allowed. Reject malformed input.
2. Starting with A=infinity, for j=0,...,m-1 write the record (key(A),j),
   then set A=A+P. At the end A=[m]P. Retain D=-A.
3. Sort the m records by the complete point key (infinity flag, X, Y), using
   bottom-up stable two-buffer merge sort, exactly 127 passes.
4. Set gamma=Q. For i=0,...,h-1, use deterministic binary search for key(gamma).
   On a match with j, form a=i*m+j, check a<r, independently compute [a]P,
   compare it with Q, and output the 32-byte canonical scalar if equal.
   Otherwise advance gamma=gamma+D. The last advance may be omitted.
5. If no match is found, report failure; no scalar is claimed for that target.

Only a successful scalar has the protected output encoding. A failure is an
unsuccessful attempt, not an alternative accepted output. The cost bound covers
all h iterations on every valid target, even when the attack fails. There is
one cold attack, with no retry, amplification or future-target amortization.

## 3. Exact success proof; no mixing assumption

The baby points [j]P, 0<=j<m, are distinct because m<r. At giant iteration i,
gamma=Q-[im]P. A table match is therefore equivalent to
x congruent to im+j modulo r for some 0<=j<m. Every such im+j lies in
[0,K), and K<r, so the congruence is equivalent to equality for canonical x.
Conversely every x in [0,K) has its unique quotient/remainder decomposition
x=im+j with i<h and j<m. Thus the algorithm succeeds exactly on K scalars.
It handles x=0 through the explicit infinity key. Outside [0,K) there is no
match, not even a modular wraparound match.

Consequently success is exactly K/r > 1/2 > 0.39 over the ORIGINAL uniform
production distribution. No independent-coin assumption is needed because the
algorithm is deterministic; averaging over fresh coins leaves this probability
unchanged. The finite exhaustive evidence tests decomposition only, not this
production-size execution. The algebraic proof supplies production coverage.

## 4. Scalar AArch64 instruction bounds for field and point arithmetic

All numbers are arrays of 64-bit limbs, never arbitrary-width unit-cost words.
Field elements occupy four limbs. Five-limb carry handling and nine-limb
products are explicitly allowed. No measured optimized multiplication price is
substituted for a complete affine addition.

For modular addition/subtraction use at most **A=256 instructions**, including
loads/stores, carry/borrow, a five-limb temporary, conditional correction and
control. A fixed small number of passes over at most five limbs suffices.

For multiplication use four-limb schoolbook multiplication followed by
Montgomery reduction with radix B=2^64, R=B^4. For an odd p, precompute
n0=-p^(-1) mod B and R^2 mod p during charged setup. Form the eight-limb
product in a nine-limb scratch area. For i=0,...,3, compute
u_i = t_i*n0 mod B, add u_i*p*B^i to the scratch, and propagate carries.
The four low limbs become zero; divide by R and conditionally subtract p.
For operands below p, T<p^2<pR; the standard reduction result is below 2p,
so a single final correction suffices.

There are 16 full limb products in the initial product, 16 in the correction
products, and four low-word products for u_i: at most 36 word-product terms.
A full product uses MUL and UMULH. Budget 128 scalar instructions per term,
including accumulation and propagation through at most nine limbs. This can
be implemented with at most 12 instructions per propagated limb and at most
20 other instructions per term: 9*12+20=128. This deliberately loose bound
includes operand loads, result stores, carry tests and loop/address work.
Reserve another 2,048 instructions for scratch initialization, shifts,
comparison/correction and return. Thus 36*128+2048=6656 < **M=8192**.
Fixed-size loops can be unrolled; no wide-integer library is assumed.
Montgomery multiplication and canonical add/sub preserve the representation.
This is an analytical instruction budget, not a measured compiler result.

A nonzero inverse by binary powering to p-2 requires at most 512 modular
multiplications, plus 2^14 instructions for fixed exponent bits and loop work:
**I=512*M+2^14**. Start with the Montgomery representation of one. Conversion
of public constants and input coordinates is charged to setup.

Complete affine addition handles infinity, equal X coordinates, opposite
points and doubling before using a denominator inverse. For distinct X use
lambda=(v-y)/(u-x); for doubling with y!=0 use lambda=3*x^2/(2*y).
Return infinity for an opposite pair or y=0 doubling. The remaining formulas
are X'=lambda^2-x-u, Y'=lambda*(x-X')-y. One inversion, at most 16 modular
multiplications, 32 modular adds/subtracts, and 2^14 other instructions suffice.
The full point-step bound is therefore

    G = I + 16*M + 32*A + 2^14
      = 4,366,336 < 2^23 instructions.

This includes arithmetic scratch loads/stores and all exceptional branches.
Scalar multiplication uses ordinary binary double-and-add, at most 512 such
steps for a 256-bit scalar. No division, inversion or scalar multiplication is
free. We use H=0 in H+R_other/C: all work is charged as scalar instructions,
then divided ONCE by the fixed **C=1000**.

## 5. Records, addresses and paid memory transfers

Each record occupies 128 bytes: two four-limb coordinates (64 bytes), an
8-byte infinity flag, a 16-byte j and 40 padding bytes. The full key is nine
64-bit words. A local key comparison is bounded by 512 scalar instructions.
Indices are two limbs; array byte offsets and global addresses use three limbs.
Binary search uses at most 128 probes. It does not use a hash oracle.

A 2^127-record array cannot be addressed by one 64-bit pointer. To make the
large-memory convention explicit, use a binary routing tree of depth 127,
with one leaf per table slot, two record buffers at each leaf, and one
controller. Each node has only small locally 64-bit-addressable memory.
There are less than 2^128 tree nodes. Only one memory transaction is outstanding
at a time. Links are a classical communication topology, not advice or a
computation oracle. Sum all active nodes' instructions; there is no elapsed-time
or parallelism discount. Idle nodes perform no polling.

A request/response carries at most 32 words: the three-limb destination,
control, record payload and reply information fit. There are at most 256
node-to-node hops for a transaction. At each hop allow 16 instructions per
word for explicit loads, stores, address arithmetic and control, or 512
instructions for a packet. Together with leaf work and the controller's
multiword address setup, a complete read OR write is below **2^18** scalar
instructions. No transfer is free. The construction needs neither bandwidth
nor latency assumptions to obtain a total-work bound: links just copy words
and their instruction-equivalent copying work has been charged. Tree traversal
follows successive bits of the public slot number. Small per-node forwarding
programs, queues and state fit within 2^14 bytes per node. This is a conventional
abstract distributed-memory work model, not a claim that the hardware exists.

Initialize the entire tree, its buffers and replicated small routing program.
Recursively assign node-prefix addresses and broadcast the same bounded program
once per edge, rather than sending a separate root-to-leaf image for every
node. Charge at most 2^20 instructions per node, so **2^148** instructions
cover tree construction, code distribution and all memory initialization.
No parameter-dependent table or machine state is supplied as free advice.
The one arithmetic controller's program is included in fixed setup.

## 6. Complete preprocessing, search and fixed-work ledger

During a merge pass, retain two current head records locally. Each record is
read once from its source buffer and written once to its destination buffer;
there are at most m key comparisons. Charge per record per pass

    2*(2^18) + 2^13 instructions,

including both routed transfers, comparisons, two-/three-limb counters,
run boundaries and local buffer work. Across 127 passes, with table generation
and its routed initial writes, the per-baby total is bounded by

    2^23 + 2^18 + 127*(2^19+2^13) < 2^27.

Loop initialization, step negation and any final result-buffer selection fit
in the fixed reserve below. Thus table generation AND sorting cost at most
**m*2^27** instructions, not 2^27 per pass. All preprocessing is charged cold.

Each giant iteration uses one point step (charged even if omitted on the last
iteration), at most 128 routed record reads, local comparisons and wide-index
binary-search work. Its bound is

    2^23 + 128*(2^18+2^13) + 2^16 < 2^26.

Thus the whole search costs at most **h*2^26** instructions. This includes
unsuccessful lookups and all iterations on failed targets.

An additional **2^48** instruction reserve covers parsing with a fixed input
size cap, public-constant derivation, curve/subgroup checking, D construction,
[scalar]P verification, the canonical scalar computation/output and controller
code initialization. Even 4096 point steps at the above bound cost below
2^35; actual scalar membership/verification needs far fewer. Hence this fixed
reserve is ample. No random bits are generated. No trials, restarts, stored
randomness or amplification occur.

Total instruction upper bound (fixed reserve counted once) is

    R_total = 2^148 + m*2^27 + h*2^26 + 2^48
            = 2^155 + 2^148 + 2^48.
    W_total = R_total/1000.

Since 1000*2^146 = 250*2^148, whereas
R_total = 129*2^148 + 2^48 < 250*2^148, we obtain
**W_total < 2^146**, proving the submitted score bound without a logarithm
rounding assumption. Total time is a worst-case bound, not expected time.

All fixed work can conservatively be assigned to preprocessing:

    R_pre = 2^148 + m*2^27 + 2^48 = 65*2^148+2^48.

Since 1000*2^145=125*2^148, **W_pre < 2^145**. This is already included
in W_total, not an extra unpaid cost. Some verification occurs after search;
assigning its reserve to preprocessing only makes that supplementary bound looser.

## 7. Peak memory, advice, limitations and provenance

The two record arrays require 2*m*128=2^135 bytes. Even assigning 2^14
bytes of routing code, buffers and state to each of fewer than 2^128 nodes
uses below 2^142 bytes. Add a 2^20-byte controller/code/input/output/scratch
reserve. The sum is below 2^143 and therefore below the declared **2^150 bytes**.
The node allowance includes edge buffers and prefix addresses. Sort scratch
is the second array; local heads fit the controller reserve. There is no
retained randomness. No nonuniform advice is used: computed records and
machine initialization are paid preprocessing and included in peak storage.

Heuristics: **none**. Production assumptions: the protected curve/order/profile
and ordinary correctness of the explicitly bounded classical arithmetic and
routing operations. No random-walk mixing, collision law or independence model
is invoked. No production attack, AArch64 binary, router implementation, or
astronomical-memory machine was executed or supplied. The instruction budget is
an analytical implementation bound permitted by the rules, and should be reviewed
as such. The routing topology makes accounting explicit; physical construction
costs beyond the protected instruction convention are not a claimed engineering
feasibility result. The score is deliberately loose; acceptance would mean only
the challenge's exploratory plausible_not_refuted status.

The method is the standard Shanks baby-step/giant-step decomposition, limited
to an algorithm-chosen coverage interval to exceed the success threshold, not
an interval-DLP promise. HAC Algorithm 3.56 supplies the standard decomposition;
HAC Algorithm 14.32 supplies Montgomery reduction. The proof above is
self-contained and adapts both to this fixed target and complete cold ledger.
The organizer's existing BSGS package motivated the accounting-control scope,
but its stale divisor, wide-word counts and score are NOT adopted. No participant
submission's negation-map, rho walk, code, measurements or claimed score is used.

References:
- Menezes, van Oorschot, Vanstone, Handbook of Applied Cryptography, Algorithm
  3.56: https://cacr.uwaterloo.ca/hac/about/chap3.pdf
- Same, Algorithm 14.32: https://cacr.uwaterloo.ca/hac/about/chap14.pdf
- Protected target: https://www.secg.org/sec2-v2.pdf#page=13
- Repository rules and cost model at source commit
  467130d867499e7bbca7479f0f72af0df115b282 (hashes in evidence.md).
