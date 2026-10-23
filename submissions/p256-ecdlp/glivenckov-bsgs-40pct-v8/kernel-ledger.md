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
chain starts at P, emitting [2j+1]P for index j=0 before advancing by [2]P.
The giant starts at Q-[L]P. There is no dependent-denominator batching trick.

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
Both production chain lengths m=2^127 and h=ceil(r/2^128) are divisible by BATCH (checked from
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

## D1: bounded-local-memory bit-routed streaming (v8 counts)

Use the SAME protocol and prices as v7, but a baby tree with m=2^127
leaves and height 127 and a giant tree with L=2m=2^128 leaves and height
128. Actual giant length h=ceil(r/L)<L; unused leaves emit END. No global
array exists. Each node has at most 2^14 local bytes and at most three
neighbor ports; the central engine adds one edge. The engine-to-leaf bounds
are thus 128 and 129 edges, respectively. Record indices fit two words,
whereas counts and exclusive endpoints use three words.

Initialize fewer than 2(m+L) nodes. Per node, 16384 instructions cover 4096 for
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
charge a further 2(m+L)*512 for all END requests/responses. All traversal and
coordination work is summed, not divided by the processor count.

The final two-stream intersection costs 512 per consumed record: record
loads, key comparisons, y-sign selection, counters and control. It may stop
on the first equal key; charging both complete streams is conservative.
There are no radix buckets, scatter arrays, unbounded file descriptors,
wide-address constant-time memory accesses or free network transfers.
In particular v7 had NO radix passes to halve: it used a comparison merge
network. Recompute E=128m+129L record-edge operations, not an imagined radix
pass count over 2^128 records. Actual records m+h < 3*2^127; they are NOT
2^128. Charged records m+L=3*2^127, versus v7's 2^129.

## New cold setup, labels and recovery

Charge [1]P, [3]P and P0=[m]P as THREE full V scalar actions even though [1]P
is public. A two-entry affine table stores P and 3P (128 bytes). Derive 2P
by a complete mixed update 3P-P; double P0 to the initial center LP=2P0,
and again to the giant stride [4]P0=2LP; negate and seed Q-LP. These are four
complete updates plus FOUR standalone normalizations and conversions.
The expanded count is 3V+4J+4*(I+7M+512)+2^20. The charged envelope is
4V+4J+2^20, and V >= 4*(I+7M+512) is checked exactly. No per-giant doubling
is needed: the fixed affine increment is computed once, as in ordinary BSGS.

Enumerate the positive odd radix-4 lattice IN NUMERIC ORDER, using one mixed
Jacobian addition of cached 2P per baby. This is not 64 fresh digit actions
per record. It emits exactly the same 64-digit signed set; the induction
and label identity are in argument.md. Do not claim dependent affine batching.
The complete K_rec=23121 schedule above applies unchanged to each chain.
Additionally charge 64 per baby: 8 increment loads, 16 index/odd-label work,
16 two-entry table-selection/pointers and 24 control. The table need not be
searched globally. All fixed addends fit the central 2^20-byte reserve.

Keep 96-byte records: x32+y32+tag8+index16+padding8. Full-y resolution gives
at most TWO candidates (2i+1)L +/- (2j+1), not four candidates from digit
magnitudes. Magnitude 3 is NOT a curve automorphism. Raw candidates may
need 257 bits, so construct them with FIVE 64-bit limbs; they are positive
and below 2hL < 3r. Two explicit five-limb conditional subtractions by r
suffice. The 2^20 final-recovery allowance covers shifts, additions,
comparisons, these reductions, odd-lift filtering and serialization; a fixed
unrolled path uses fewer than 4096 instructions, independent of the streams.
Both 256-bit reduced candidates still receive their own V verification.

Peak <=2(m+L)*2^14+2^20 < 2^144 bytes. Every byte of node code/state,
queues, links, heads, tags, counters and scratch is inside the allocations.
No physically realizable deployment or production execution is claimed.
