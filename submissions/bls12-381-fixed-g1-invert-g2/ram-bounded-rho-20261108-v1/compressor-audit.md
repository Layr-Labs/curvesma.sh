# Complete unrolled compressor schedule and preparation

The mathematical U is unchanged from universal-proof.md. This revision tightens only its complete scalar cost by spelling out its leaf implementation. The source remains evidence; no assembly or Python experiment executes on the host.

Use registers for five Horner accumulators, five key residues, q, the shared digit, lo, hi, one reduction temporary, input/output pointers and one extraction temporary. These fit in19 caller-saved registers x0 through x18. The routine is a leaf and preserves LR; no calls,extra spills,stack frame or implicit big-integer operations occur. On platforms reserving x18, use x19 with one saved/restored pair; the128 setup/output allowance pays those instructions. Key residues are loaded once, accumulators initialized to zero once, and q is loaded in at most four instructions.

Every update h<-h*k+c mod q follows this exact13-instruction arithmetic pattern, with the obvious operand register names:

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

All operands are64-bit scalar registers. Correctness follows from2^61=1 mod q and h,k<q,c<2^60. The128-bit product is represented by two words; neither half is free. First folding sum is<2q<2^62; second is<=q; final conditional subtraction canonicalizes. This is13 instructions per update,195 updates for five Horner evaluations over39 digits:2535. Fully unroll the39 iterations with digits from38 down to0; each digit is extracted once and shared by all five updates. There is no loop-counter work per update.

The aligned input is37 little-indexed64-bit limbs d_0,...,d_36 of M, the integer represented by b_bytes||sign. Each60-bit digit uses at most two aligned LDRs,two immediate shifts,one ORR and one AND mask. Allow8 instructions per digit,312 total, including the edge case of a zero shift, which omits rather than executes a shift by64. Every offset is a fixed immediate below296 bytes. Initialization,keys,constant construction,pointers,at most one preserved-register pair,output's five byte reversals/stores and RET cost<128 more. Thus the leaf needs at most2535+312+128=2975; charge3072, leaving97 instructions of slack. All loads,stores,addressing and control in this routine are included.

The input limb preparation and final SHA message preparation use the separate
512-instruction schedule in preparation-audit.md. The 3072 compressor allowance
above includes its digit extraction, five Horner accumulators, key loads, output
byte reversals and stores, all ABI work and RET. None of those costs are moved
into the smaller preparation allowance. No experiment is declared.
