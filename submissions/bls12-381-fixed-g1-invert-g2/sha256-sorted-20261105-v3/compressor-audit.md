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

The aligned input is37 little-indexed64-bit limbs d_0,...,d_36 of M, the integer represented by b_bytes||sign. Each60-bit digit uses at most two aligned LDRs,two immediate shifts,one ORR and one AND mask. Allow8 instructions per digit,312 total, including the edge case of a zero shift, which omits rather than executes a shift by64. Every offset is a fixed immediate below296 bytes. Initialization,keys,constant construction,pointers,at most one preserved-register pair,output's five byte reversals/stores and RET cost<128 more. Thus the leaf needs at most2535+312+128=2975; charge4096, leaving1121 instructions of slack. All loads,stores,addressing and control in this routine are included.

The input limb buffer preparation is charged in the SHA's existing2048 preparation allowance, which now processes only40 hash-payload bytes. Let B_0,...,B_35 be the six canonical Fp coefficients of b read as36 big-endian64-bit words, in the protected encoding order, and let bit be its separately charged sign byte. Form d_0=(B_35<<8)|bit; d_j=(B_(35-j)<<8)|(B_(36-j)>>56), 1<=j<=35; d_36=B_0>>56. Load and reverse each B word once, retain the neighboring word in a register, and store each d as one aligned machine word. This uses fewer than512 instructions including pointers and local control. All37 words are initialized, so extraction never reads uninitialized padding. It reconstructs M exactly, including leading zero bits. No unaligned accesses are assumed.

Within the same2048 preparation allowance,copy the precomputed32-byte SHA-256 chaining state,form one64-byte block containing40-byte U output,0x80 padding,15 zero bytes and an eight-byte length field832,and extract the final24-bit label. These fixed loads/stores,reversals and control cost below512 more. The complete aligned compressor input plus hash preparation is below1024. sha256-audit.md separately charges64*40+48*20+128=3648 for all rounds,schedule and internal control. Adding2048 preparation gives5696;charge6144. No hidden whole-object copy occurs.

A full transition is45312+4096+6144+2048+8192=65792;charge66048,with sorted-checkpoint bookkeeping audited in checkpoint-proof.md. For the current shorter L and K,the exact inequality (K*66048+2^90)^50<20000^50*2^6469 supports score6469/50=129.38. All probability bounds,including compression loss below2^-13,and memory below2^35 remain as in the revised finite proof. No executable experiment is declared; these are analytic bounds for shared review.
