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

Within the same2048 preparation allowance, copy the precomputed64-byte SHA chaining state, form one128-byte block containing the40-byte U output,0x80 padding and16-byte length field for1344 bits, and extract the final24-bit label. These fixed loads/stores,reversals and control cost<512 more. The total preparatory work<1024 leaves over1024 further instructions within that2048 charge. The per-block80*64+64*24+128 round/schedule bound already includes internal loads,stores,addressing and block feed-forward. The tighter unrolled SHA audit in dense-sha-audit.md uses80*40+64*20+128+2048=6656;charge7168. There is no displaced work or uncharged whole-object copy.

Consequently a full transition is39168+4096+7168+2048+8192=60672;charge60928, including the original nine-checkpoint accounting. Exact public arithmetic yields (K*60928+2^90)^40<20000^40*2^5179, supporting rational score5179/40=129.475. All probability bounds,including the key-compression coupling loss<2^-13,and peak memory below2^35 bytes are unchanged. The bounded arithmetic experiment checks the mathematical reduction,digit reconstruction,Lucas certificate and public work inequality only. It neither measures this assembly nor proves the SHA heuristic.
