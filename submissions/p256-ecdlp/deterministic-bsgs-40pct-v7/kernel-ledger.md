# Instruction derivation and executable-mechanism scope

This is an analytical scalar-AArch64 schedule, not a measurement of Python,
assembly or production hardware. `mechanism.py` checks arithmetic and the
stream protocol at bounded sizes. `build_certificate.py` checks every row
sum, exact power comparison, profile pin and captured-result count.

## D3: the price source and one consistent four-limb kernel

`tracks/docs/COST_MODEL.md`, Computation accounting and Complete attack work,
requires arithmetic, logical operations, memory, addressing and control on
scalar 64-bit AArch64 registers, all counted as instructions. Its target table
sets P-256 C=1000. Its kernel measurement table reports 179, 22996 and 383 for
the three named *other* kernels. None is a P-256 bound and none is used to
interpolate one here. The per-instruction price used in EVERY row below is
exactly 1, from protected `tracks/classical-group-work-v2.json`,
`primitive_operations.scalar_AArch64_instruction`. COST_MODEL.md contains no
opcode-by-opcode price table; inventing one would be a false citation.
MUL and UMULH each count separately, just like ADDS, ADC, LDR, STR and branch.

Let B=2^64, R=B^4. Four canonical limbs store each Montgomery residue.
Use four CIOS rounds, initially t=0:

```
for i in 0..3:
    t += a * b[i]            # four limb-pair accumulations
    q = low_word(t)          # p'=-p^(-1) mod B=1 for this exact P-256 p
    t += q * p               # four MORE limb-pair accumulations
    assert low_word(t)==0
    t /= B                  # register renaming/copies, not an arbitrary shift
return t if t<p else t-p
```

The exact modulus is p=2^256-2^224+2^192+2^96-1. There is no free
pseudo-Mersenne reduction. In particular the reduction is not omitted:
all sixteen q*p limb products and their carry propagation are paid below.
No special sparse-modulus speedup is used except the exact p'=1 identity.
The invariant t<2p after each division follows from
`(2p + p*(B-1) + p*(B-1))/B = 2p`; strict input bounds make it strict.
Before division t<2pB, so a nine-word scratch accumulator is more than enough.
At exit one five-word subtraction and selection gives canonical output.

A limb-pair accumulator uses MUL+UMULH, then ADDS for the low word, ADCS
for the high word and at most seven ADCs to ripple across nine scratch words.
That is two multiply instructions plus at most nine carry instructions.
Keep the scratch in registers: nine t words, eight operand words, four
modulus words, two product words, three argument pointers, the return
address and at most three scratch/control words fit within thirty registers.
q can alias the low t word. Explicit save/restore allowance is also present.
An unrolled schedule needs no arbitrary-width arithmetic instruction.

| M kernel row | Count x price 1 |
|---|---:|
| Schoolbook products, 16 MUL and 16 UMULH | 32 |
| Schoolbook carries, 16 x 9 | 144 |
| Four REDC passes, 4 x 4 MUL/UMULH pairs | 32 |
| REDC carries, 16 x 9 | 144 |
| Quotient aliases/copies | 4 |
| Operand/modulus loads, output stores, argument traffic | 24 |
| Scratch zeroing | 9 |
| Final five-word subtract, select and moves | 16 |
| Callee save/restore and spill/load allowance | 128 |
| Four shift/copy schedules and round control | 64 |
| Entry, return, pointer construction and moves | 40 |
| Raw ceiling | 637 |
| INSIDE-bound safety, ceil(6*637/5) | **M=765** |

The wide carry allowance is deliberately larger than the six-word CIOS
accumulator requires. The spill allowance is not relied upon to hold the
entire accumulator in memory: the register schedule just given establishes
that it can remain resident. Squaring is the same kernel with equal inputs,
charged S=M. A conversion out of Montgomery form is M, not free.

| A kernel row | Count x price 1 |
|---|---:|
| Two operands and modulus: 12 limb loads | 12 |
| Sum with carry-out, at most 5 ADDS/ADCS | 5 |
| Conditional-reduction subtract or correction | 5 |
| Four-word selection | 4 |
| Output stores | 4 |
| Calls, pointers and control | 10 |
| Raw ceiling | 40 |
| INSIDE-bound safety, ceil(6*40/5) | **A=48** |

For subtraction, four SUBS/SBCS and four masked-modulus correction additions
fit the same arithmetic rows. For addition the fifth word preserves overflow;
selection uses the full five-word comparison, not a wrapped low-word test.
The same bound covers modular negation. All code uses canonical residues.

## D2: projective generation is not affine batch addition

The fixed addend is converted to Montgomery affine form once in setup.
Store a Jacobian chain state (X,Y,Z) and copy its current value to a local
batch; advance the *projective* state using complete case handling. The baby
chain starts at O, emitting the record for index zero before advancing.
The giant starts at Q. There is no dependent-denominator batching trick.

For a finite unequal mixed addition use the schedule:

```
z2=Z^2; u=x2*z2; s=y2*Z*z2
h=u-X; d=s-Y
hh=h^2; ii=4hh; jj=h*ii; rr=2d; vv=X*ii
X3=rr^2-jj-2vv
Y3=rr*(vv-X3)-2Y*jj
Z3=(Z+h)^2-z2-hh
```

It uses eleven M/S calls and at most twenty A calls. If h=0,d!=0, return O;
if h=d=0, double using the exact a=-3 schedule:

```
delta=Z^2; gamma=Y^2; beta=X*gamma
alpha=3*(X-delta)*(X+delta)
X3=alpha^2-8beta
Z3=(Y+Z)^2-gamma-delta
Y3=alpha*(4beta-X3)-8gamma^2
```

This is eight M/S calls and at most twenty A calls. Four M calls precede
branch selection, so even the double branch costs at most twelve M and
thirty-two A, including branch pretests. Identity inputs and Y=0 return the
appropriately tagged state. Charge EVERY step `J=16M+32A+512`; the final
allowance covers coordinate copies, full-limb tests, branches and pointers.
This is an upper bound for the complete mixed operation, not a charge for
both mutually exclusive formulas. It also bounds doubling used in scalar
multiplication. These schedules are the a=-3 Jacobian EFD mixed-addition and
doubling identities; `mechanism.py` tests them against an independent oracle.

For each BATCH=1024 finite-or-tagged states, replace zero Z by Montgomery one.
Prefix multiplication and the backward product walk cost at most 3M per
point. Invert the product by fixed exponent p-2, using at most 512 calls of M,
64 control instructions per call and 256 setup instructions:
`I=512*(M+64)+256`. Separately charge `ceil(I/BATCH)` per point.
Both production chain lengths m and h are divisible by BATCH (checked from
the exact profile); there is no uncharged residual inversion.

Then charge x scaling as `z^-2` and `X*z^-2` (2M), y scaling as a freshly
computed `z^-2`, `z^-3` and `Y*z^-3` (3M), and conversion of both resulting
coordinates out of Montgomery form (2M). Recomputing z^-2 is intentional:
no shared operation is charged and then silently reused as a substitute for
another charge. Infinity still goes through the bounded masked path.
Array traffic and mask/control receive 256 instructions per element; record
assembly, index/tag/padding stores and local send preparation receive 512.

Thus the displayed generation rows sum EXACTLY to
`K_rec = J + 3M + ceil(I/1024) + 2M + 3M + 2M + 256 + 512`.
The generated ledger shows each integer and its exact sum. Inversion is PLUS
the three multiplications, never hidden inside them. Full scalar verification
costs `V=512J+I+7M+512`: complete double-and-add, standalone normalization,
both conversions and point comparison. R_verify pays for two V calls.

## D1: bounded-local-memory bit-routed streaming

There are two full binary trees of height 128, each with m=2^128 leaves.
No processor addresses any global array. Each has at most 2^14 bytes, native
local pointers and three neighbor ports. A port is a local small integer,
not a global processor identifier. A record carries its two-word index;
bit routing selects one child per level. The central engine adds one edge,
so a record traverses at most 129 edges between engine and leaf. This
explicit finite network is a theoretical resource, not a practical machine.

Initialize fewer than 4m nodes. Per node, 16384 instructions cover 4096 for
8KiB code distribution at both edge endpoints, 4096 for bounded local zeroing,
2048 for neighbor/link and prefix setup, 1024 for metadata and 5120 for local
control. Recursively distribute code and initialize links: no free global
allocator or free code replication. Each local 16KiB allocation includes
code, record caches, neighbor queues, stack, tags, counters and links. The
central batch state/input/code reserve is separate.

One record placement per edge costs 128: 48 for unrolled twelve-word transfer
loads/stores at both endpoints, plus five rows of 16 for route-bit selection,
queue state, local pointers, tags, and call/control. Every index bit is read
through the two native words; no 128-bit pointer is used.

Each leaf then streams its one record or END. Each internal node caches its
two child heads, emits the smaller (ties left first), and requests replacement
only for that child. Thus induction gives a stable sorted stream keyed by
(tag,x); no record is dropped. One record across one edge/merge level costs
512: request+response word transfers 52, full-key comparison 64, cached-record
copies 48, state and end checks 96, scheduling/calls 128, local pointers 64,
and flags/counters 60. Includes both endpoints. There is no busy-waiting.
Every record passes each level once, and every edge emits END only once;
charge a further 4m*512 for all END requests/responses. All traversal and
coordination work is summed, not divided by the processor count.

The final two-stream intersection costs 512 per consumed record: record
loads, key comparisons, y-sign selection, counters and control. It may stop
on the first equal key; charging both complete streams is conservative.
There are no radix buckets, scatter arrays, unbounded file descriptors,
wide-address constant-time memory accesses or free network transfers.

## Scope of the captured checks

The generator checks CIOS on bounded full-width operands, all ordinary,
identity, opposite and double Jacobian branches, nontrivial projective
rescalings, normalization, progression labels, sorted-stream multiplicities
and duplicate keys. It exhausts every scalar of a counted prime-order toy
curve, and tests selected actual P-256 fixtures including Q=-P. It also
checks negative controls for the predecessor's sign rule. These are executed
mechanism checks, not estimates of production success, instruction prices,
network throughput or physical feasibility. The production coverage proof
and the analytical instruction schedules carry those claims independently.
