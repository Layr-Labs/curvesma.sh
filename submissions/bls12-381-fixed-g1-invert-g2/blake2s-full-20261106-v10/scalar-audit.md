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
is performed; full states stay stored; exact injective keys are compared.

At the first nonzero a coefficient x, sign is x>floor(p/2), equivalently
x>=H=(p+1)/2; earlier zero coefficients equal their negations. Fully unroll six
coefficient tests, each six LDR, five ORR, one CBNZ: at most72 instructions.
Keep the selected six limbs, load six H limbs, compare by SUBS/five SBCS and
CSET CS:13. All-zero case outputs zero. Entry pointers, optional BL, 64-bit sign
store, return and control fit32. Total<=117; charge128. Only caller-saved
registers and immediate offsets below288/48 are needed. b serialization is
charged in hash preparation, not treated as free. Canonical Montgomery residues
are unique modulo p, so tuple and field equality are consistent.

## Five exact polynomial dot products:1312 including leaf control

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

Each shared digit c_l is extracted by four scalar instructions. Set constant
j=floor(60l/64), t=(60l) mod64 and fully unroll:

    LDR lower,[input,#8*j]
    LDR upper,[input,#8*(j+1)]
    EXTR digit,upper,lower,#t
    AND digit,digit,#0x0fffffffffffffff

EXTR produces the low64 bits of the concatenated128-bit(upper:lower) shifted
right by t. Its immediate t is0..63, including valid t=0. The logical-immediate
mask has60 contiguous one bits and is encodable in one AND. Both loads are
aligned with fixed offsets<=288; for the last digit j=35,j+1=36, inside the
initialized37-limb buffer. No dynamic division/multiplication, pointer increment,
branch, spill or shift-by64 case occurs. Upper/lower registers reuse the product
half temporaries before those products are computed, preserving the20-GPR bound.
39 digits cost156, with each extracted once and shared across five keys.
Load c_0 once and initialize all ten accumulator registers within the
retained128 setup/output allowance. Ten accumulator registers, digit, power,
product halves, source/power/output pointers, q/mask58 and temporary use at most20
GPRs. Use x0..x17 plus x19/x20, preserve the latter pair with STP/LDP. This is a
leaf, preserving LR with no internal calls. The128 covers that pair, pointers,
q/mask initialization, ten initializations, five byte reversals/stores and RET.
No register stack spill occurs in the fully unrolled inner updates.

Complete cost950+60+156+128=1294; charge1312 with18 spare instructions. The old
Horner calculation and this integer dot product are equal modulo q by the
polynomial identity, so U and every actual hash output are unchanged.

The exact aligned input preparation is charged separately in the1024 hash
preparation allowance below. The1520-byte powers table fits the inherited2^50
peak storage; its reads are charged in the950 term instructions above.

## Standard full-round BLAKE2s-256:2432 including1024 preparation

The walk's fixed64-byte prefix WALK||seed32||zero16 is compressed once in paid
setup. The state-dependent second block contains U's40 bytes and24 zero padding bytes
with cumulative byte counter104 and final flag true. All preparations below use
ordinary scalar instructions and preallocated aligned buffers; no allocation,
memcpy, memory oracle or implicit per-step work is omitted.

For clarity, define the physical b coefficients in Fp6 order(v^0,v^1,v^2),
Fp2 real then imaginary, each as six canonical little-indexed64-bit Montgomery
limbs. The numerical big-endian b words B_i have fixed physical offsets
48*floor(i/6)+8*(5-(i mod6)), i=0..35, from b's base. A native scalar LDR from
that offset already gives the numerical value B_i; it does not reverse the value.
The conceptual big-endian bytes specify the encoding, not an uncharged copy.

Form the37 little-indexed64-bit limbs of M=b_bytes||sign:
 d_0=(B_35<<8)|sign;
 d_j=(B_(35-j)<<8)|(B_(36-j)>>56), j=1..35;
 d_36=B_0>>56.
A fully unrolled schedule starts with LDR B_35, LDRB sign, LSL, ORR, STR (5).
For each next B, keep previous B in a register and use LDR current, LSL temp,#8,
ORR temp,temp,previous,LSR #56, STR temp, MOV previous,current (5 each).
The last limb is LSR previous,#56 and STR (2). Sum5+35*5+2=182 instructions.
For every loaded B word, add one STR to its fixed slot in the exact-key buffer,
36 additional instructions. The sign routine writes its0/1 as one64-bit STR at
key offset288 instead of a byte STR, the same instruction count; LDRB for M
reads its low byte correctly. Thus182+36+64=282 including pointers/control,
within the retained512 preparation allowance with230 spare. Key b offsets are
0..280, all fixed/aligned. Key padding8 bytes is initialized zero once in setup. All LDR offsets<=280 and
STR offsets<=288 are fixed, aligned, valid and initialized. This implements the
same injective289-byte encoding and supplies the same39 digits.

Before compression, zero the64-byte final-block buffer with four STP xzr,xzr at
0,16,32,48. The compressor's five already charged REV/STR outputs write U's40
bytes directly into0..39; bytes40..63 remain zero. There is no SHA padding marker
or bit-length field. Copy the32-byte precomputed prefix h into a mutable state
buffer with two LDP/STP pairs(4). The final kernel uses counter104, final true,
last-node false; its initialization charges those constants below.

BLAKE2s serializes its numeric state words little-endian. The low40 digest bits
under the specified big-endian digest interpretation are
 j=((state[6]>>24)<<32)|REV32(state[7]).
Two LDR,LSR,REV,ORR use5 instructions; all values are zero-extended before the
64-bit shift/or. Charge256 for final-block/state/label work including its4 zero
stores,4 state-copy instructions,5 label instructions,pointers and constants.

Retain256 more for wrapper argument formation, its own caller bookkeeping,
additional pointer setup, buffer selection and any wrapper spills/control not
inside the block function's128 overhead. This wrapper budget does not replace
sign/compressor/reference interfaces, which retain their separate named costs.
Total preparation512+256+256=1024. Buffers and stack frames are constructed in
cold setup and reused; their accesses are all charged here or in the actual
leaf routines. No production timing or free call assumption enters.

Standard BLAKE2s-256 uses10 rounds, each8 G operations. For each specified G on
distinct32-bit state registers a,b,c,d, this complete scalar A64 schedule costs16:

    ADD a,a,b
    LDR message_word,[block,#x_offset]
    ADD a,a,message_word
    EOR d,d,a
    ROR d,d,#16
    ADD c,c,d
    EOR b,b,c
    ROR b,b,#12
    ADD a,a,b
    LDR message_word,[block,#y_offset]
    ADD a,a,message_word
    EOR d,d,a
    ROR d,d,#8
    ADD c,c,d
    EOR b,b,c
    ROR b,b,#7

All data operands are W registers, modulo2^32; block is an X address register.
The two message loads are included. RFC7693's10 SIGMA permutations and8 G index
quadruples are unrolled at assembly time into fixed registers and byte offsets
0..60. There is no runtime schedule lookup,counter,loop,permutation or register
shuffle. Standard little-endian32-bit message words are read directly by LDR on
the selected little-endian AArch64 machine. Full cost80 G*16=1280, not a timing.

Use W0..W15 for16 working values, W16 for message data,X17 for the block,
X19 for h andX20 for IV constants. Twenty GPRs suffice without x18. Save/restore
x19/x20 as one STP/LDP pair; this is a leaf preserving LR. Before clobbering
argument registers, move h/block bases into x19/x17. Load8 h words and8 IV words,
16 instructions. The final walk block's t=104, high counter0 and final flag true
need MOV104,EOR into v12,MVN v14:3. v13 XOR0 and v15 XOR0 are omitted, not free
nonconstant operations. After all rounds, feed-forward each h_i XOR v_i XOR v_(i+8)
with LDR,EOR,EOR,STR,32 instructions. Prologue,pointer moves,IV ADRP/ADD,
epilogue and RET keep the complete fixed overhead below128; the enumerated
16+3+32 arithmetic/memory overhead leaves77 instructions for these fixed tasks.
There is no per-G call or callee-saved spill. Both halves of each G are fully
counted; no SIMD or cryptographic extension is used.

The unkeyed32-byte hash parameter state is initialized in cold setup as
h_0=IV_0 XOR0x01010020, h_i=IV_i fori1..7. Prefix compression at counter64,
final false,is also paid in setup. The repeated final-block kernel starts from
that saved prefix h, not an already finalized prefix digest. The code computes
standard full10-round BLAKE2s-256, not a reduced function. No timing claim or
formal security theorem is imported from the specification.

Total full kernel1280+128+1024=2432, including every charged preparation above.
The new hash-model bridge is explicit in proof.md and is not inferred from an
earlier SHA review. The mathematical first-repeat proof needs only the idealized
initial values/labels and therefore transfers inside this new model.

## Exact norm-one key search, coefficients and counters:4640

The key is37 native numerical64-bit words: B_0..B_35 followed by sign0/1.
Its injectivity on full GT follows from a²-vb²=1: given b there are at most two
roots a,-a in the odd-characteristic field; lexicographic sign distinguishes
these unless a=0, when they are the same root. The key is not a probabilistic
fingerprint. Sorting keys depends only on group elements, preserving conditional
uniformity of the denominator and the first-equality stopping rule.

For each of18 pairs, unroll:

    LDP x1,x2,[x0],#16
    LDP x3,x4,[x5],#16
    CMP x1,x3
    B.ne mismatch
    CMP x2,x4
    B.ne mismatch

Then compare the final sign word with LDR from each advanced base, CMP,B.ne.
The first18 pairs cost108 and the last word4,112 total even on full equality.
Post-index16 handles pointer advancement; no oversized offset or increment is
hidden. Bases are distinct from loaded registers. Flags preserve unsigned word
order; all37 words are initialized. The unused8 padding bytes are not compared.
Retain96 for argument bases/key offsets, BL/RET, lookup, lo/hi/m, interval update,
return result and saved registers. Total208; charge224 including each search
iteration,16 spare. At most127 saved keys need7 iterations:1568 instructions.
No square root, decompression, hash collision or free equality oracle enters.

Both four-limb scalar updates retain512: each4 LDP operands,ADDS/three ADCS,
2 LDP r,SUBS/three SBCS,4 CSEL,2 STP=20 plus128 caller/control. r<2^255 prevents
loss of a fifth carry. Retain256 index/address including cached-label load/store,
512 multiword counters/checkpoint control,512 state-record copying,256 equality/
dispatch,1024 spills/control. Sum1568+512+256+512+512+256+1024=4640.

Saved/current records are944 bytes: the original640 state/scalar bytes plus304
aligned key bytes. A complete save needs59 LDP/STP pairs=118 instructions using
post-index16, well within512 with pointers/control. All padding is initialized
zero in charged setup. Up to127 sorted insertions and shifts remain below2^30
fixed instructions, charged within2^90. Jump-table records stay640 bytes and
indexing j*640 is unchanged.127 checkpoint records total119888 bytes, still
negligible within2^50. State arithmetic uses full X as before, not reconstructed
values. The extra initial key/hash preparation is fixed work inside2^90; each
subsequent transition builds one key/hash label for its new state before lookup.

## Reference-call interface:64, explicit complete schedule

Keep current/next944-byte state-record pointers in x19/x20 and selected jump
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

Total nonreference cost1312+2432+128+4640+64=8576 per transition. The fixed
setup/finish budget, peak memory, exact total work and sole heuristic are fully
specified in proof.md and probability.md. The optional Python fixture checks
finite hash semantics only and has not been executed on the host.

Arm Instruction Set Reference Guide D2.83 documents scalar EXTR and its64-bit
immediate range0..63: https://documentation-service.arm.com/static/6245c734b059dc5ff9a8bdab
