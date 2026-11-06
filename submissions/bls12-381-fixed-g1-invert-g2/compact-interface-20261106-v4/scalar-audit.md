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

## Five Horner compressions:3072, including leaf control

q=2^61-1 is prime: q-1=2*3²*5²*7*11*13*31*41*61*151*331*1321.
All factors are prime by trial division through37. Base37 satisfies
37^(q-1)=1 mod q and gcd(37^((q-1)/f)-1,q)=1 for each distinct prime factor f;
the full-factor Lucas criterion proves primality. These small constants and
certificate work are charged in setup. Distinct289-byte encodings have distinct
39-digit degree<=38 polynomials because each60-bit digit is below q.
Independent uniform five keys give the fixed-pair root probability(38/q)^5.

Five accumulators and keys, q, shared digit, lo/hi/tmp, pointers and extraction
temporary fit19 GPRs. On ABIs reserving x18, use x19 with a saved/restored pair
charged in128 setup/output. Each h<-h*k+c mod q, h,k<q,c<2^60, is:

    MUL lo,h,k
    UMULH hi,h,k
    ADDS lo,lo,c
    ADC hi,hi,xzr
    LSR tmp,lo,#61
    ORR tmp,tmp,hi,LSL #3
    AND lo,lo,q
    ADD lo,lo,tmp
    LSR tmp,lo,#61
    AND lo,lo,q
    ADD lo,lo,tmp
    SUBS tmp,lo,q
    CSEL h,tmp,lo,HS

The two-word product is<2^122 plus a small digit. First61-bit fold is<2q;
second is at most q+1, and the final subtraction canonicalizes. There is no
free high product half or division. Fully unroll39 digits descending and five
keys:195 updates at13 instructions=2535. Each60-bit digit is shared across keys
and uses at most two aligned LDR, two immediate shifts, ORR and AND; charge8,
including edge handling that omits a shift by64:312. Initialization, keys/q,
pointers, saved pair if needed, five output reversals/stores and RET fit128.
Total2975<3072; there is no runtime loop or per-update spill.

The aligned input is37 little-indexed64-bit limbs of the big-endian encoding
b_bytes||sign. With B_0..B_35 the big-endian64-bit b words, form
 d_0=(B_35<<8)|sign;
 d_j=(B_(35-j)<<8)|(B_(36-j)>>56), j=1..35;
 d_36=B_0>>56.
Load/reverse each B word once, retain a neighbor and store each initialized d;
this costs<512 including control. The scalar digit offsets are fixed below296.
This preparation is included in the unchanged SHA2048 allowance below. Padding
and all37 initialized limbs make every extraction valid, including leading zeros.

## Standard SHA-512:4832 including2048 preparation

The remaining walk block has40 payload bytes,0x80 padding and16-byte original
length1344. The fixed128-byte prefix's chaining value is precomputed in charged
setup; no hash output changes. Copying the64-byte chaining state, forming the
128-byte final block and extracting40 label bits costs<512 instructions. Together
with<512 aligned-limb preparation above this uses<1024, leaving>1024 in the
retained2048 preparation bound. No whole-object copy is displaced or omitted.

A fully unrolled expansion word for i=16..79 uses14 instructions, with all
symbolic offsets below replaced by immediate bytes0..632 at assembly time:

    LDR t0,[sp,#8*(i-15)]
    ROR t1,t0,#1
    EOR t1,t1,t0,ROR #8
    EOR t1,t1,t0,LSR #7
    LDR t0,[sp,#8*(i-2)]
    ROR t2,t0,#19
    EOR t2,t2,t0,ROR #61
    EOR t2,t2,t0,LSR #6
    LDR t0,[sp,#8*(i-16)]
    ADD t1,t1,t0
    LDR t0,[sp,#8*(i-7)]
    ADD t1,t1,t0
    ADD t1,t1,t2
    STR t1,[sp,#8*i]

It computes sigma0(W_(i-15))+sigma1(W_(i-2))+W_(i-16)+W_(i-7) modulo2^64.
No branch, counter, pointer increment or spill occurs. Three temporaries suffice;
round-state initialization follows expansion and does not compete for registers.

For logical state(a,b,c,d,e,f,g,h), fully unroll this22-instruction round:

    ROR t0,e,#14
    EOR t0,t0,e,ROR #18
    EOR t0,t0,e,ROR #41
    EOR t1,f,g
    AND t1,t1,e
    EOR t1,t1,g
    ADD t2,h,t0
    ADD t2,t2,t1
    LDR t3,[wp,#word_offset]
    ADD t2,t2,t3
    LDR t3,[kp,#constant_offset]
    ADD t2,t2,t3
    ROR t0,a,#28
    EOR t0,t0,a,ROR #34
    EOR t0,t0,a,ROR #39
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
word/constant offsets are<=632. After80 rounds the original physical mapping
returns because80 is divisible by8. No internal call or loop exists.

Fixed block overhead128 covers48 initial LDR/REV/STR word preparation,8 working
state loads and24 feed-forward LDR/ADD/STR, plus48 for entry/exit,640-byte frame
and pointer setup. There is no callee-saved spill with14 caller-saved registers.
Thus80*22+64*14+128+2048=4832. SHA semantics are FIPS180-4 sections4.1.3/6.4;
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

Total nonreference cost3072+4832+128+4992+64=13088 per transition. The fixed
setup/finish budget, peak memory, exact total work and sole heuristic are fully
specified in proof.md and probability.md. The optional Python fixture checks
finite hash semantics only and has not been executed on the host.
