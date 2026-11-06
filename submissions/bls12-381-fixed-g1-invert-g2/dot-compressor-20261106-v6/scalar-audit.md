# Complete scalar work outside one reference multiplication

All counts are scalar AArch64 instructions, including loads, stores, addressing
and control. No SIMD, arbitrary-width unit primitive or free oracle is used.
C remains the protected20000. Reference multiplication is charged as H=1;
its own instructions, stack and ABI preservation are not charged again.

## Norm-one injective encoding and sign:128

For X=a+bw, conjugation X->a-bw is p^6 Frobenius. Public exponentiation gives
p^6=-1 mod r, so a²-vb²=1 for every GT element. Given b, at most two a roots
exist, a and-a; if a=0 there is one. The unsigned canonical Montgomery
coefficient tuple's lexicographic sign distinguishes them. Thus six48-byte b
coefficients plus one sign byte are injective. No square root or decompression
is performed; full states stay stored and compared.

At the first nonzero a coefficient x, sign is x>floor(p/2), equivalently
x>=H=(p+1)/2; earlier zero coefficients equal their negations. Fully unroll six
coefficient tests, each six LDR, five ORR, one CBNZ: at most72 instructions.
Keep the selected six limbs, load six H limbs, compare by SUBS/five SBCS and
CSET CS:13. All-zero case outputs zero. Entry pointers, optional BL, sign-byte
store, return and control fit32. Total<=117; charge128. Only caller-saved
registers and immediate offsets below288/48 are needed. b serialization is
charged in SHA preparation, not treated as free. Canonical Montgomery residues
are unique modulo p, so tuple and field equality are consistent.

## Five exact polynomial dot products:1472 including leaf control

q=2^61-1 is prime. Its complete factorization q-1 is
2*3²*5²*7*11*13*31*41*61*151*331*1321. Each listed factor is prime by trial
division through37; base37 satisfies37^(q-1)=1 mod q and
 gcd(37^((q-1)/f)-1,q)=1 for each distinct factor f.
The full-factor Lucas criterion certifies q. Certificate/constants construction
are charged in setup. Distinct39-digit encodings define different polynomials
over Fq because digits<2^60<q; degree38 gives the same five-key root bound.

During cold setup construct190 canonical powers T_(i,l)=k_i^l mod q for
five keys and l=1..38, by successive modular multiplication. Store them in
1520 bytes with8-byte alignment;190 small fixed products/control are within
2^90. There is no precomputed answer/advice or free table. The powers of zero
keys are valid zeros, with no division or special failure. T_(i,0)=1 is implicit.

Initialize five128-bit accumulators (hi,lo)=(0,c_0). For each l=1..38, extract
one shared digit and execute this five-instruction scalar sequence per key:

    LDR power,[powers,#constant_offset]
    MUL product_lo,digit,power
    UMULH product_hi,digit,power
    ADDS accumulator_lo,accumulator_lo,product_lo
    ADC accumulator_hi,accumulator_hi,product_hi

All are64-bit scalar registers, not arbitrary-width primitives. Fully unroll
38*5=190 updates with fixed table offsets below1520, costing950 instructions.
No runtime indexing, loop counter, spill or pointer increment occurs per term.
Each product is<2^121 because digit<2^60 and power<q<2^61. Each complete
accumulator is<38*2^121+2^60<2^127, so high-word ADC cannot overflow64 bits.
Intermediate partial sums are nonnegative and smaller than that bound.

For each final A=lo+2^64*hi, form its three61-bit chunks and fold using
2^61=1 mod q. With mask58=2^58-1, fully unroll these12 instructions:

    AND t0,lo,q
    LSR t1,lo,#61
    AND t2,hi,mask58
    ORR t1,t1,t2,LSL #3
    LSR t2,hi,#58
    ADD t0,t0,t1
    ADD t0,t0,t2
    LSR t1,t0,#61
    AND t0,t0,q
    ADD t0,t0,t1
    SUBS t1,t0,q
    CSEL output,t1,t0,HS

The first chunks are A mod2^61, floor(A/2^61) mod2^61, floor(A/2^122).
The last is<32 because A<2^127. Their sum is<2q+32<2^63 and fits one word.
One further61-bit fold is at most q+2<2q, so one conditional subtraction gives
the canonical exact A mod q. This handles all boundary and zero-key cases.
Five reductions cost60. No wide division, hidden reduction loop or carry is free.

Retain the original8-instruction allowance per shared60-bit digit,312 for39
digits: two aligned loads, up to two constant shifts, ORR, AND and edge/control
slack. Load c_0 once and initialize all ten accumulator registers within the
retained128 setup/output allowance. Ten accumulator registers, digit, power,
product halves, source/power/output pointers, q/mask58 and temporary use at most20
GPRs. Use x0..x17 plus x19/x20, preserve the latter pair with STP/LDP. This is a
leaf, preserving LR with no internal calls. The128 covers that pair, pointers,
q/mask initialization, ten initializations, five byte reversals/stores and RET.
No register stack spill occurs in the fully unrolled inner updates.

Complete cost950+60+312+128=1450; charge1472 with22 spare instructions. The old
Horner calculation and this integer dot product are equal modulo q by the
polynomial identity, so U and every actual hash output are unchanged.

Aligned input preparation is unchanged and paid in SHA2048. For the big-endian
encoding b_bytes||sign, read/reverse36 b words B_0..B_35 once and write37 aligned
little-indexed64-bit limbs:
 d_0=(B_35<<8)|sign;
 d_j=(B_(35-j)<<8)|(B_(36-j)>>56), j=1..35;
 d_36=B_0>>56.
This costs<512 including pointers/control. All37 limbs are initialized; digit
extraction uses fixed valid offsets below296. No unaligned access or uncharged
whole-object copy is assumed. The extra1520-byte powers table fits inherited
2^50 peak storage and all reads are charged in the950 instructions above.

## Standard SHA-256:4256 including2048 preparation

The remaining walk block has40 payload bytes,0x80 padding and8-byte original
length832. The fixed64-byte prefix's chaining value is precomputed in charged
setup; no hash output changes. Copying the32-byte chaining state, forming the
64-byte final block and extracting40 label bits costs<512 instructions. Together
with<512 aligned-limb preparation above this uses<1024, leaving>1024 in the
retained2048 preparation bound. No whole-object copy is displaced or omitted.

A fully unrolled expansion word for i=16..63 uses14 instructions, with all
symbolic offsets below replaced by immediate bytes0..252 at assembly time:

    LDR w0,[sp,#4*(i-15)]
    ROR w1,w0,#7
    EOR w1,w1,w0,ROR #18
    EOR w1,w1,w0,LSR #3
    LDR w0,[sp,#4*(i-2)]
    ROR w2,w0,#17
    EOR w2,w2,w0,ROR #19
    EOR w2,w2,w0,LSR #10
    LDR w0,[sp,#4*(i-16)]
    ADD w1,w1,w0
    LDR w0,[sp,#4*(i-7)]
    ADD w1,w1,w0
    ADD w1,w1,w2
    STR w1,[sp,#4*i]

It computes sigma0(W_(i-15))+sigma1(W_(i-2))+W_(i-16)+W_(i-7) modulo2^32.
No branch, counter, pointer increment or spill occurs. Three temporaries suffice;
round-state initialization follows expansion and does not compete for registers.

For logical32-bit state(a,b,c,d,e,f,g,h), fully unroll this22-instruction round.
All named state/temporary operands below are W registers, LDR/STR are32-bit,
wp/kp/sp remain64-bit address registers, and additions wrap modulo2^32:

    ROR t0,e,#6
    EOR t0,t0,e,ROR #11
    EOR t0,t0,e,ROR #25
    EOR t1,f,g
    AND t1,t1,e
    EOR t1,t1,g
    ADD t2,h,t0
    ADD t2,t2,t1
    LDR t3,[wp,#word_offset]
    ADD t2,t2,t3
    LDR t3,[kp,#constant_offset]
    ADD t2,t2,t3
    ROR t0,a,#2
    EOR t0,t0,a,ROR #13
    EOR t0,t0,a,ROR #22
    ORR t1,a,b
    AND t1,t1,c
    AND t3,a,b
    ORR t1,t1,t3
    ADD t0,t0,t1
    ADD d,d,t2
    ADD h,t2,t0

t2 is T1, t0 is T2; the final writes give new e in old d and new a in old h.
The next compile-time register tuple is(old h,old a,old b,old c,old d,old e,
old f,old g). All old values are read before overwrite. No runtime moves occur.
The majority identity is(a&b)|(c&(a|b)); choice is g^(e&(f^g)). Eight state,
four temporary and two pointer registers fit x0..x13, all caller-saved. All
word/constant offsets are<=252. After64 rounds the original physical mapping
returns because64 is divisible by8. No internal call or loop exists.

Fixed block overhead128 covers48 initial32-bit LDR/REV/STR word preparation,8 working
state loads and24 feed-forward LDR/ADD/STR, plus48 for entry/exit,256-byte frame
and pointer setup. There is no callee-saved spill with14 caller-saved GPRs using their W views.
Thus64*22+48*14+128+2048=4256. SHA semantics are FIPS180-4 sections4.1.2/6.2;
rotated EOR and LDP are ordinary scalar A64 instructions. This is an analytical
instruction schedule, not a CPU timing or Python-instruction measurement.

## Full-state search, coefficient updates and counters:4992

Keys remain all72 canonical64-bit words. Each of36 fully unrolled word pairs is:

    LDP x1,x2,[x0],#16
    LDP x3,x4,[x5],#16
    CMP x1,x3
    B.ne mismatch
    CMP x2,x4
    B.ne mismatch

x0/x5 are distinct temporary bases, not destinations. Post-index16 is encodable
and includes pointer advancement; no oversized immediate offset or separate
increment is needed.640-byte record strides give16-byte alignment. Both pairs
are loaded within initialized record bytes. The first unequal word's flags
preserve the same unsigned lexicographic result. On equality all72 words agree.
No compressed-key or hash comparison replaces full equality.

36*6=216; retain96 for caller BL/RET, base resets, lookup, lo/hi/m computation,
branch/return result and preserved registers. Charge320 including each binary-
search iteration, eight more than312. Six iterations cost1920, even if keys
agree through71 words. Copies of immutable pointers preserve coefficients.

Both four-limb scalar updates together cost<=512: each uses4 LDP for operands,
ADDS/three ADCS,2 LDP r, SUBS/three SBCS,4 CSEL,2 STP,20 instructions plus
128 caller/control. Since r<2^255, adding two residues loses no fifth carry.
Retain256 index/address,512 multiword counters/checkpoint control,512 occasional
640-byte state-record copy,256 equality/dispatch,1024 spills/control.
Sum1920+512+256+512+512+256+1024=4992. Sorted insertion is in fixed2^90, not
hidden in a per-step average. There is no visited-state dictionary.

## Reference-call interface:64, explicit complete schedule

Keep current/next640-byte state-record pointers in x19/x20 and selected jump
record pointer in x21. They are initialized in charged setup; selecting x21 and
j*640 addressing belong to the256 index/address allowance. At the reference call:

    MOV x0,x20
    MOV x1,x19
    MOV x2,x21
    BL blst_fp12_mul
    MOV x9,x19
    MOV x19,x20
    MOV x20,x9

These seven scalar instructions pass output/current/jump pointers and swap
buffers. The reference multiplication writes its576-byte result directly to
next storage. It is the pinned ordinary ABI function and preserves x19..x21;
its own return, spills and memory work are inside the once-counted reference
operation. Coefficient updates write next record's64 coefficient bytes directly
within their512 allowance. Other helpers preserve these callee-saved pointers
or pay their own saved-register traffic in their explicit allowances. Their
argument preparations remain in the unchanged named budgets above.

The complete transition loop's counters, flags and dispatch already have1024
outer control instructions. No extra per-step wrapper call, tower dispatcher,
state copy or scratch allocation is introduced by this reference call. Charge64
rather than the exact seven instructions, leaving57 for any additional interface
moves, pointer materialization or caller scheduling. Function-entry preservation
of long-lived registers is fixed setup, paid in2^90; it does not recur each step.
This tighter explicit interface, rather than an assumed free call, is the new
change from the previous384 allowance. No callee internals are charged twice.

Total nonreference cost1472+4256+128+4992+64=10912 per transition. The fixed
setup/finish budget, peak memory, exact total work and sole heuristic are fully
specified in proof.md and probability.md. The optional Python fixture checks
finite hash semantics only and has not been executed on the host.
